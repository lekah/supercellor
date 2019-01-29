
import numpy as np

import math
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice

try:
    from pymatgen.util.coord import lattice_points_in_supercell
except ImportError:
    from pymatgen.util.coord_utils import lattice_points_in_supercell

from pymatgen.core.sites import Site, PeriodicSite
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from spglib import standardize_cell

#from supercellor.lib.optimal_supercell import utils, optimal_supercell_hnf


EPSILON = 1e-6 # The precision when comparing floats!

def rotate(R, lattice_vecs, frac_coords):
    new_lattice_vecs = np.dot(lattice_vecs.T, R).T
    new_frac_coords = np.dot(frac_coords, R)
    return new_lattice_vecs, new_frac_coords

def get_angles(cell_matrix):
    a, b, c = [v/np.linalg.norm(v) for v in cell_matrix]
    k_a = np.arccos(np.dot(b, c))
    k_b = np.arccos(np.dot(a, c))
    k_c = np.arccos(np.dot(a, b))
    return k_a, k_b, k_c


def standardize_cell(cell, wrap):
    frac_coords = np.empty((len(cell), 3))
    cellvecs = cell._lattice.matrix
    for i, site in enumerate(cell):
        frac_coords[i,:] = site.frac_coords

    # The code below sorts the cell and positions by lattice vector length,
    # shortest first:
    veclens = sorted([(np.linalg.norm(v), i) for i, v in enumerate(cellvecs)])
    M1 = np.zeros((3,3))
    for row, (_, idx) in enumerate(veclens):
        M1[idx, row] = 1
    cellvecs, frac_coords = rotate(M1, cellvecs, frac_coords)

    k_a, k_b, k_c = get_angles(cellvecs)
    right_angle = 0.5*np.pi
    if ((k_a <= right_angle and k_b <= right_angle and k_c <= right_angle ) or 
        (k_a > right_angle and k_b > right_angle and k_c > right_angle )):
        M2 = np.eye(3)
    elif ((k_a <= right_angle and k_b > right_angle and k_c > right_angle ) or 
        (k_a > right_angle and k_b <= right_angle and k_c <= right_angle )):
        M2 = np.diag([1,-1,-1])
    elif ((k_a > right_angle and k_b <= right_angle and k_c > right_angle ) or 
        (k_a <= right_angle and k_b > right_angle and k_c <= right_angle )):
        M2 = np.diag([-1,1,-1])
    elif ((k_a > right_angle and k_b > right_angle and k_c <= right_angle ) or 
        (k_a <= right_angle and k_b <= right_angle and k_c > right_angle )):
        M2 = np.diag([-1,-1,1])
    else:
        raise RuntimeError("Unrecognized case for k_a={}, k_b={}, k_c={}".format(k_a,k_b, k_c))
    cellvecs, frac_coords = rotate(M2, cellvecs, frac_coords)
    # Now applying the rules layed out in  http://arxiv.org/abs/1506.01455
    # to get the standardized triclinic cell.
    # Since not all cells give to me are triclinic < -> <= with respect to
    # the paper. So it's not truly standardized!
    metric = np.dot(cellvecs, cellvecs.T)

    a = np.sqrt(metric[0, 0])
    b = np.sqrt(metric[1, 1])
    c = np.sqrt(metric[2, 2])
    alpha = np.arccos(metric[1, 2] / b / c)
    beta = np.arccos(metric[0][2] / a / c)
    gamma = np.arccos(metric[0][1] / a / b)

    cg = np.cos(gamma)
    cb = np.cos(beta)
    ca = np.cos(alpha)
    sg = np.sin(gamma)

    cellvecs = np.zeros((3,3))
    cellvecs[0, 0] = a
    cellvecs[1, 0] = b * cg
    cellvecs[2, 0] = c * cb
    cellvecs[1, 1] = b * sg
    cellvecs[2, 1] = c * (ca - cb * cg) / sg
    cellvecs[2, 2] = c * np.sqrt(1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg) / sg
    
    # And some checks:
    if (
        cellvecs[0,0] < -1e-12 or
        abs(cellvecs[0,1]) > 1e-12 or
        abs(cellvecs[0,2]) > 1e-12):
        raise ValueError("First lattice vector not aligned with x-axis")
    if  (
        cellvecs[1,1] < -1e-12 or
        abs(cellvecs[1,2]) > 1e-12):
        raise ValueError("Second lattice vector not in X-Y plane in first quadrant")
    if  (cellvecs[2,2] < 0):
        raise ValueError("Third lattice vector not in positive Z-direction")
    
    new_sites = []
    new_lattice = Lattice(cellvecs)

    for i, site in enumerate(cell):
        s = PeriodicSite(site.species_and_occu, frac_coords[i],
                        new_lattice, properties=site.properties,
                        coords_are_cartesian=False, to_unit_cell=wrap)
        new_sites.append(s)
    
    cell = Structure.from_sites(new_sites)
    return cell

def get_diagonal_solution(cell, r_inner):
    # getting the vectors:
    a,b,c = [np.array(v) for v in cell]
    R_diag = np.zeros((3,3), dtype=int)
    #repetitions =[None]*3
    for index, (n1, n2, n3) in enumerate(((a,b,c), (b,c,a), (c,a,b))):
        d23 = np.cross(n2, n3) 
        d1 = np.linalg.norm(np.dot(n1, d23)) / np.linalg.norm(d23)
        R_diag[index, index] = int(np.ceil(2*r_inner/d1)) # x2 to get diameter
    return R_diag, np.dot(R_diag, cell).T


def get_possible_solutions(cell, r_inner, verbosity=1):
    """
    Returns all possible vectors (as integer crystal coordinates)
    that lie outside a sphere of radius min_dist_soll, but within a sphere given by a
    maximum radius that is determined by the upper bound, the diagonal solution.
    The trial vectors could be iteratively reduced with better solutions found, while still being exhaustive
    """
    def apply_sym_ops(G, cell):
        """
        Apply possible symmetry operations.
        For now, i am reducing by removing all negative indices in the third direction
        """
        return G[G[:,2] >= 0, :]
        
    # Would be good to first LLL/Niggle reduce to make the upper bound as low as possible.
    
    #atoms = structure.get_ase()
    R_diag, C_diag = get_diagonal_solution(cell, r_inner)

    v_diag = np.abs(np.dot(np.cross(C_diag[0], C_diag[1]), C_diag[2]))
    r_outer = v_diag / r_inner**2 / 8.0

    S = np.matrix(np.diag([1,1,1,-r_outer**2]))

    cellT = cell.T
    cellI = np.matrix(cell).I
    cellTI = np.matrix(cellT).I
    #  I describe the move from atomic to crystal coordinates with an affine transformation M:
    M=  np.matrix(np.r_[np.c_[np.matrix(cellTI), np.zeros(3)], [[0,0,0,1]]])
    # Q is a check, but not used. Check is orthogonality
    # Q is the sphere transformed by transformation M
    Q =  M.I.T * S * M.I
    # Now, as defined in the source, I calculate R = Q^(-1)
    R = M * S.I *M.T
    # The boundaries are given by:
    boundaries = np.zeros((3,2), dtype=int) # results truncated to integer
    boundaries[0,0] = (R[0,3] + np.sqrt(R[0,3]**2 - R[0,0]*R[3,3])) / R[3,3] - EPSILON
    boundaries[0,1] = (R[0,3] - np.sqrt(R[0,3]**2 - R[0,0]*R[3,3])) / R[3,3] + EPSILON
    

    boundaries[1,0] = (R[1,3] + np.sqrt(R[1,3]**2 - R[1,1]*R[3,3])) / R[3,3] - EPSILON
    boundaries[1,1] = (R[1,3] - np.sqrt(R[1,3]**2 - R[1,1]*R[3,3])) / R[3,3] + EPSILON

    boundaries[2,0] = (R[2,3] + np.sqrt(R[2,3]**2 - R[2,2]*R[3,3])) / R[3,3] - EPSILON
    boundaries[2,1] = (R[2,3] - np.sqrt(R[2,3]**2 - R[2,2]*R[3,3])) / R[3,3] + EPSILON

    # I create the first reduced grid, that is all the grids within the bounding box:
    Gc_r0 =  np.array([a.flatten() for a 
            in np.meshgrid(*[np.arange(lbound, ubound+1, 1)
                    for lbound, ubound in boundaries])]).T

    if verbosity > 1:
        print ('Gc_r0')
        for item in Gc_r0:
            print '  {:>3} {:>3} {:>3}'.format(*item)
    # I reduced the grid further by apply symmetry operations to the grid.
    # I do this before calculating any distances, since distances are expensive
    Gc_r1 = apply_sym_ops(Gc_r0, cell)
    # getting the gridpoints in real space, no longer the grid
    Gr_r1  = np.dot(Gc_r1, cell)
    # calculating the norm:
    norms_of_Gr_r1 = np.linalg.norm(Gr_r1, axis=1)
    msk = (norms_of_Gr_r1 > ( r_inner -EPSILON) ) & (norms_of_Gr_r1 < (r_outer + EPSILON) ) 

    Gr_r2 = Gr_r1[msk, :]
    Gc_r2 = Gc_r1[msk, :]
    norms_of_Gr_r2 = norms_of_Gr_r1[msk]


    if verbosity > 1:
        print ('Gc_r01')
        for item, norm in zip(Gc_r1, norms_of_Gr_r1):
            print '  {:>3} {:>3} {:>3} {norm}'.format(*item, norm=norm)
        print ('Gc_r02')
        for item, norm in zip(Gc_r2, norms_of_Gr_r2):
            print '  {:>3} {:>3} {:>3} {norm}'.format(*item, norm=norm)
        print r_inner, r_outer
        print R_diag
    # Now I am sorting everything via the norm:
    sorted_argindex = norms_of_Gr_r2.argsort()
    sorted_Gr_r2 = Gr_r2[sorted_argindex, :]
    sorted_Gc_r2 = Gc_r2[sorted_argindex, :]
    norms_of_sorted_Gr_r2 = norms_of_Gr_r2[sorted_argindex]

    return norms_of_sorted_Gr_r2, sorted_Gc_r2, sorted_Gr_r2, r_outer, v_diag



def make_supercell(structure, r_inner, method='bec', wrap=True, standardize=True,
        do_niggli_first=True, verbosity=1):
    """
    Creates from a given structure a supercell based on the required minimal dimension
    :param structure: The pymatgen structure to create the supercell for
    :param float r_inner: The minimum image distance as a float,
        The cell created will not have any periodic image below this distance
    :param str method: The method to get the optimal supercell. For now, the only
        implemented option is *best enclosing cell*
    :param bool wrap: Wrap the atoms into the created cell after replication
    :param bool standardize: Standardize the created cell.
        This is done based on the rules in Hinuma etal, http://arxiv.org/abs/1506.01455
        However, only rules for the triclinic case are applied, so further
        standardization using spglib is recommended, if a truly standardized
        cell is required.
    :param bool do_niggli_first: Start with a niggli reduction of the cell,
        to enable a faster search. Disable if there are problems with the reduction
        of if the cell is already Niggli or LLL reduced.
    :param int verbosity: Sets the verbosity level.

    :returns: A new pymatgen core structure instance and the used scaling matrix
    :returns: The scaling matrix used.
    """
    if not isinstance(structure, Structure):
        raise TypeError("Structure passed has to be a pymatgen structure")
    try:
        r_inner = float(r_inner)
        assert r_inner>1e-12, "Non-positive number"
    except Exception as e:
        print "You have to pass positive float or integer as r_inner"
        raise e

    if not isinstance(wrap, bool):
        raise TypeError("wrap has to be a boolean")
    if not isinstance(standardize, bool):
        raise TypeError("standardize has to be a boolean")

    # I'm getting the niggli reduced structure as first:
    if verbosity > 1:
        print("given cell:\n", structure._lattice)
    if do_niggli_first:
        starting_structure = structure.get_reduced_structure(reduction_algo=u'niggli')
    else:
        starting_structure = structure
    if verbosity > 1:
        print( "starting cell:\n", starting_structure._lattice)
        for i, v in enumerate(starting_structure._lattice.matrix):
            print( i, np.linalg.norm(v))
    
    # the lattice of the niggle reduced structure:
    lattice_cellvecs = np.array(starting_structure._lattice.matrix, dtype=np.float64)
    # trial_vecs are all possible vectors sorted by the norm
    if method == 'bec':
        reduced_solutions = get_possible_solutions(lattice_cellvecs, r_inner,
                            verbosity=verbosity)
        if verbosity:
            print( "I received {} possible solutions".format(len(reduced_solutions[0])))
        # I pass these trial vectors into the function to find the minimum volume:
        scale_matrix, supercell_cellvecs = get_optimal_solution(
                *reduced_solutions, r_inner=r_inner, verbosity=verbosity)
    elif method == 'hnf':
        raise NotImplementedError("HNF has not been fully implemented")
        lattice_cellvecs = np.array(lattice_cellvecs)
        scale_matrix, supercell_cellvecs = optimal_supercell_hnf(lattice_cellvecs, r_inner, verbosity)
    else:
        raise ValueError("Unknown method {}".format(method))
    # Constructing the new lattice:
    new_lattice = Lattice(supercell_cellvecs)
    # I create f_lat, which are the fractional lattice points of the niggle_reduced:
    f_lat = lattice_points_in_supercell(scale_matrix)
    # and transforrm to cartesian coords here:
    c_lat = new_lattice.get_cartesian_coords(f_lat)
    #~ cellT = supercell_cellvecs.T

    if verbosity > 1:
        print("Given lattice:\n", new_lattice)
        for i, v in enumerate(new_lattice.matrix):
            print(i, np.linalg.norm(v))

    new_sites = []
    if verbosity:
        print("Done, constructing structure")
    for site in starting_structure:
        for v in c_lat:
            new_sites.append(PeriodicSite(site.species_and_occu, site.coords +v,
                             new_lattice, properties=site.properties,
                             coords_are_cartesian=True, to_unit_cell=wrap))

    supercell = Structure.from_sites(new_sites)

    if standardize:
        supercell = standardize_cell(supercell, wrap)
        if verbosity > 1:
            print("Cell after standardization:\n",  new_lattice)
            for i, v in enumerate(new_lattice.matrix):
                print (i, np.linalg.norm(v))
    return supercell, scale_matrix



def get_optimal_solution(norms_of_sorted_Gr_r2, sorted_Gc_r2, sorted_Gr_r2, r_outer, v_diag,
        r_inner, verbosity=1):

    min_volume = np.inf
    max_min_inter_face_dist = 0
    #max_radius = trial_vecs[-1][0]
    indices = range(len(norms_of_sorted_Gr_r2))
    R_best = np.zeros((3,3), dtype=int)
    C_best = np.zeros((3,3), dtype=float)
    
    for i1, norm1 in enumerate(norms_of_sorted_Gr_r2, start=0):
        if norm1 > r_outer + EPSILON:
            # At this point I have finished!
            break
        vector1 = sorted_Gr_r2[i1]
        for i2, norm2 in enumerate(norms_of_sorted_Gr_r2[i1+1:], start=i1+1):
            if norm2 > r_outer + EPSILON:
                # I can stop iterating over the possible second vectors
                break
            # Checking the dot product, so that I continue if the vectors have an angle < 60
            vector2 = sorted_Gr_r2[i2]
            if np.abs(np.dot(vector1, vector2)) / (norm1*norm2) >= 0.5:
                if verbosity > 1:
                    print '   -> Angle < 60, continue'
                continue
            for i3, norm3 in enumerate(norms_of_sorted_Gr_r2[i2+1:], start=i2+1):
                if norm3 > r_outer:
                    if verbosity > 1:
                        print '     -> Max radius surpassed, break'
                    break
                vector3 = sorted_Gr_r2[i3]
                if np.abs(np.dot(vector2, vector3)) / (norm2*norm3) >= 0.5:
                    if verbosity > 1:
                        print '     -> Angle < 60, continue'
                    continue
                elif np.abs(np.dot(vector1, vector3)) / (norm1*norm3) >= 0.5:
                    if verbosity > 1:
                        print '     -> Angle < 60, continue'
                    continue
                # checking intersections of each plane
                cross23 = np.cross(vector2, vector3)
                d1 = np.abs(np.dot(cross23/np.linalg.norm(cross23), vector1))
                if d1 < r_inner:
                    if verbosity > 1:
                        print '     -> d1 {} < r_inner, continue'.format(d1)
                    continue
                cross13 = np.cross(vector1, vector3)
                d2 = np.abs(np.dot(cross13/np.linalg.norm(cross13), vector2))
                if d2 < r_inner:
                    if verbosity > 1:
                        print '     -> d2 {} < r_inner, continue'.format(d2)
                    continue
                cross12 = np.cross(vector1, vector2)
                d3 = np.abs(np.dot(cross12/np.linalg.norm(cross12), vector3))
                if d3 < r_inner:
                    if verbosity > 1:
                        print '     -> d3 {} < r_inner, continue'.format(d3)
                    continue

                volume = np.abs(
                        sorted_Gc_r2[i1, 0]*sorted_Gc_r2[i2, 1]*sorted_Gc_r2[i3, 2]+
                        sorted_Gc_r2[i1, 1]*sorted_Gc_r2[i2, 2]*sorted_Gc_r2[i3, 0]+
                        sorted_Gc_r2[i1, 2]*sorted_Gc_r2[i2, 0]*sorted_Gc_r2[i3, 1]-
                        sorted_Gc_r2[i1, 1]*sorted_Gc_r2[i2, 0]*sorted_Gc_r2[i3, 2]-
                        sorted_Gc_r2[i1, 2]*sorted_Gc_r2[i2, 1]*sorted_Gc_r2[i3, 0]-
                        sorted_Gc_r2[i1, 0]*sorted_Gc_r2[i2, 2]*sorted_Gc_r2[i3, 1])

                min_inter_face_dist = min((d1,d2, d3))

                if volume < min_volume or (
                        volume == min_volume and 
                        min_inter_face_dist > max_min_inter_face_dist):
                    min_volume = volume
                    r_outer = min_volume/r_inner**2
                    max_min_inter_face_dist = min_inter_face_dist

                    if verbosity:
                        print "New optimal supercell {} & {}".format(volume, max_min_inter_face_dist)
                    R_best[0,:] = sorted_Gc_r2[i1]
                    R_best[1,:] = sorted_Gc_r2[i2]
                    R_best[2,:] = sorted_Gc_r2[i3]
                    C_best[0,:] = vector1
                    C_best[1,:] = vector2
                    C_best[2,:] = vector3

    if verbosity:
        print('FINAL\n')
        print(R_best)
    return R_best, C_best

