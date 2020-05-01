import math
import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.core.lattice import Lattice
from pymatgen.core.sites import  PeriodicSite
try:
    from pymatgen.util.coord import lattice_points_in_supercell
except ImportError:
    from pymatgen.util.coord_utils import lattice_points_in_supercell


from supercellor.lib.optimal_supercell import (
    utils, fort_optimal_supercell_hnf, fort_optimal_supercell_bec)


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


def standardize_cell(structure, wrap=True):
    """
    :param structure: A pymatgen structure instance
    :param bool wrap: Whether to wrap the positions into the cell after standardizing, defaults to True.
    """
    frac_coords = np.empty((len(structure.sites), 3))
    cellvecs = structure._lattice.matrix
    for i, site in enumerate(structure.sites):
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

    for i, site in enumerate(structure.sites):
        s = PeriodicSite(site.species_and_occu, frac_coords[i],
                        new_lattice, properties=site.properties,
                        coords_are_cartesian=False, to_unit_cell=wrap)
        new_sites.append(s)
    
    standardized_structure = Structure.from_sites(new_sites)
    return standardized_structure

def get_diagonal_solution_bec(cell, d_inner):
    # getting the vectors:
    a,b,c = [np.array(v) for v in cell]
    R_diag = np.zeros((3,3), dtype=int)
    #repetitions =[None]*3
    for index, (n1, n2, n3) in enumerate(((a,b,c), (b,c,a), (c,a,b))):
        d23 = np.cross(n2, n3) 
        d1 = np.linalg.norm(np.dot(n1, d23)) / np.linalg.norm(d23)
        R_diag[index, index] = int(np.ceil(d_inner/d1)) # x2 to get diameter
    return R_diag, np.dot(R_diag, cell)

def get_diagonal_solution_hnf(cell, dmpi):
    a,b,c = [np.array(v) for v in cell]
    norms = np.linalg.norm(cell, axis=1)
    R_diag = np.diag(list(map(int, np.ceil(dmpi / norms ))))
    return R_diag, np.dot(R_diag, cell)

def get_possible_solutions(cell, d_inner, verbosity=1):
    """
    Returns all possible vectors (as integer crystal coordinates)
    that lie outside a sphere of radius r_inner=0.5 d_inner, but within a sphere given by a
    maximum radius that is determined by the upper bound, the diagonal solution.
    The trial vectors could be iteratively reduced with better solutions found,
    while still being exhaustive.
    """
    def apply_sym_ops(G, cell):
        """
        Apply possible symmetry operations.
        For now, i am reducing by removing all negative indices in the third direction
        """
        return G[G[:,2] >= 0, :]


    # Important: I work with everything times 2 here, therefore transforming
    # a radius to a diameter, and the 'half integer grid G' becomes an integer grid
    # G. This is done for having the convenience of having largely integer math.

    R_diag, C_diag = get_diagonal_solution_bec(cell, d_inner)

    # calculating the volume of the diagonal solution:
    v_diag = np.abs(np.dot(np.cross(C_diag[0], C_diag[1]), C_diag[2]))

    # r_outer is taken from the diagonal solution.
    d_outer = v_diag / d_inner**2 #/ 8.0

    S = np.matrix(np.diag([1,1,1,-d_outer**2]))

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
            print ('  {:>3} {:>3} {:>3}'.format(*item))
    # I reduced the grid further by apply symmetry operations to the grid.
    # I do this before calculating any distances, since distances are expensive
    Gc_r1 = apply_sym_ops(Gc_r0, cell)
    # getting the gridpoints in real space.
    # The dot product of Gc_r0 and cell gives the matrix of npoints x 3,
    # where every point is a coordinate.
    Gr_r1  = np.dot(Gc_r1, cell)
    # calculating the norm:
    norms_of_Gr_r1 = np.linalg.norm(Gr_r1, axis=1)


    msk = (norms_of_Gr_r1 > ( d_inner -EPSILON) ) & (norms_of_Gr_r1 < (d_outer + EPSILON) )

    Gr_r2 = Gr_r1[msk, :]
    Gc_r2 = Gc_r1[msk, :]
    norms_of_Gr_r2 = norms_of_Gr_r1[msk]

    if verbosity > 1:
        print ('Gc_r01')
        for item, norm in zip(Gc_r1, norms_of_Gr_r1):
            print ('  {:>3} {:>3} {:>3} {norm}'.format(*item, norm=norm))
        print ('Gc_r02')
        for item, norm in zip(Gc_r2, norms_of_Gr_r2):
            print ('  {:>3} {:>3} {:>3} {norm}'.format(*item, norm=norm))

    # Now I am sorting everything via the norm:
    sorted_argindex = norms_of_Gr_r2.argsort()
    sorted_Gr_r2 = Gr_r2[sorted_argindex, :]
    sorted_Gc_r2 = Gc_r2[sorted_argindex, :]
    norms_of_sorted_Gr_r2 = norms_of_Gr_r2[sorted_argindex]
    return norms_of_sorted_Gr_r2, sorted_Gc_r2, sorted_Gr_r2, d_outer, v_diag


def get_optimal_solution_bec(norms_of_sorted_Gr_r2, sorted_Gc_r2, sorted_Gr_r2, r_outer, v_diag,
        r_inner, verbosity=1):
    """
    Gettung the optimal solution for the best enclosing cell for a sphere of a given
    radius.
    """
    min_volume = np.inf
    max_min_inter_face_dist = 0
    #max_radius = trial_vecs[-1][0]
    indices = range(len(norms_of_sorted_Gr_r2))
    R_best = np.zeros((3,3), dtype=int)
    C_best = np.zeros((3,3), dtype=float)

    for i1, norm1 in enumerate(norms_of_sorted_Gr_r2, start=0):
        if norm1 > r_outer - EPSILON:
            # At this point I have finished!
            break
        vector1 = sorted_Gr_r2[i1]
        if verbosity > 1:
            print ('  Setting vector1', vector1)
        for i2, norm2 in enumerate(norms_of_sorted_Gr_r2[i1+1:], start=i1+1):
            if norm2 > r_outer - EPSILON:
                # I can stop iterating over the possible second vectors
                break

            # Checking the dot product, so that I continue if the vectors have an angle < 60
            vector2 = sorted_Gr_r2[i2]
            if verbosity > 1:
                print ('    Setting vector2', vector2)
            if np.abs(np.dot(vector1, vector2)) / (norm1*norm2) - EPSILON > 0.5:
                if verbosity > 1:
                    print ('   -> Angle < 60, continue')
                continue
            for i3, norm3 in enumerate(norms_of_sorted_Gr_r2[i2+1:], start=i2+1):
                if norm3 > r_outer - EPSILON:
                    if verbosity > 1:
                        print ('     -> Max radius surpassed, break')
                    break
                vector3 = sorted_Gr_r2[i3]
                if verbosity > 1:
                    print ('      Setting vector3', vector3)

                if np.abs(np.dot(vector2, vector3)) / (norm2*norm3) - EPSILON  > 0.5:
                    if verbosity > 1:
                        print ('     -> Angle < 60, continue')
                    continue
                elif np.abs(np.dot(vector1, vector3)) / (norm1*norm3) - EPSILON > 0.5:
                    if verbosity > 1:
                        print ('     -> Angle < 60, continue')
                    continue
                # checking intersections of each plane
                cross23 = np.cross(vector2, vector3)
                d1 = np.abs(np.dot(cross23/np.linalg.norm(cross23), vector1))
                if d1 < r_inner:
                    if verbosity > 1:
                        print ('     -> d1 {} < r_inner, continue'.format(d1))
                    continue
                cross13 = np.cross(vector1, vector3)
                d2 = np.abs(np.dot(cross13/np.linalg.norm(cross13), vector2))
                if d2 < r_inner:
                    if verbosity > 1:
                        print('     -> d2 {} < r_inner, continue'.format(d2))
                    continue
                cross12 = np.cross(vector1, vector2)
                d3 = np.abs(np.dot(cross12/np.linalg.norm(cross12), vector3))
                if d3 < r_inner:
                    if verbosity > 1:
                        print ('     -> d3 {} < r_inner, continue'.format(d3))
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
                        print("New optimal supercell {} & {}".format(volume, max_min_inter_face_dist))
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


 
def get_optimal_solution_hnf(prim_cell, dmpi, verbosity=0):
    """
    :param prim_cell: The primitive cell as a numpy array
    :param float dmpi: The minimum image distance
    :param int verbosity: The verbosity level
    """
    R_diag, C_diag = get_diagonal_solution_hnf(prim_cell, dmpi)
    v_diag = np.abs(np.dot(np.cross(C_diag[0], C_diag[1]), C_diag[2]))
    # calculating the maximal volume as number of unit cell
    # not using linalg, as it returns float
    V_max = np.prod(np.diag(R_diag)) # np.linalg.det(R_diag)
    V_prim = v_diag / V_max
    hnf = np.zeros((3,3), dtype=int)
    best_dmpi = dmpi
    inv_prim_cell = np.linalg.inv(prim_cell)
    found = False
    # I assume that the minimal volume is a close packed fcc structure,
    # with angles of 60 degrees. This gives that the total volume as a function
    # of the cube volume
    # Is there a better way to estimate V_min?
    V_min = int(dmpi**3 / np.sqrt(2) /V_prim) or 1 # in case the first gives 0
    if verbosity > 0:
        print("Prim volume is {}".format(V_prim))
        print("Testing HNF from {} to {}".format(V_min, V_max))
    # Important. Rows and colums are exchanged, so I have to build the upper
    # Hermite normal form!
    count = 0
    for V_s in range(V_min, V_max+1):
        if verbosity >1:
            print("ENTERING VOLUME", V_s)
        for a in range(1, V_max+1):
            if V_s % a:
                # If there is a modulo, I need to continue
                continue
            # this remains an integer, division integer by integer:
            quotient = int(V_s / a)
            hnf[0,0] = a
            for c in range(1, quotient+1):
                if quotient % c:
                    continue
                f = int(quotient/c)
                hnf[1,1] = c
                hnf[2,2] = f
                for b in range(0, c):
                    hnf[0,1] = b
                    for d in range(0, f):
                        hnf[0,2] = d
                        for e in range(0, f):
                            hnf[1, 2] = e
                            count += 1
                            C = np.dot(hnf, prim_cell)
                            reduced = Lattice(C).get_lll_reduced_lattice(delta=0.75)
                            C_red = reduced.matrix
                            shortest_dist = np.linalg.norm(reduced.matrix, axis=1).min()
                            
                            R = np.rint(np.dot(C_red, inv_prim_cell)).astype(int)
                            assert np.sum((np.dot(R, prim_cell)-C_red)**2) < 1e-12, "cannot do inversion properly"
                            #~ if verbosity > 1 and count==114:
                            if verbosity > 1:
                                print( '--------------------------')
                                for key, val in (('HNF', hnf), ('C', C),
                                        ('R', R), ('Reduced', C_red),
                                        ('shortest', shortest_dist)):
                                    print ('{}:'.format(key))
                                    print (val)
                            if shortest_dist > best_dmpi:
                                best_abs_norm = np.sum(np.abs(hnf))
                                R_best = R.copy()
                                C_best = C_red.copy()
                                # Setting found for true, which means I will not
                                # go to the next volume
                                found = True
                            elif abs(shortest_dist-best_dmpi) < EPSILON:
                                abs_norm = np.sum(np.abs(hnf))
                                if abs_norm < best_abs_norm:
                                    R_best = R.copy()
                                    C_best = C_red.copy()
                                    best_abs_norm = abs_norm
                                    found = True
        if found:
            break
    return R_best, C_best



def make_supercell(structure, distance, method='bec', wrap=True, standardize=True,
        do_niggli_first=True, diagonal=False, implementation='fort', verbosity=1):
    """
    Creates from a given structure a supercell based on the required minimal dimension
    :param structure: The pymatgen structure to create the supercell for
    :param float distance: The minimum image distance as a float,
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
    :param bool diagonal: Whether to return the diagonal solution, instead
        of the optimal cell.
    :param str implementation: Either fortran ('fort') or python-implementation ('pyth'),
        defaults to 'fort'
    :param int verbosity: Sets the verbosity level.

    :returns: A new pymatgen core structure instance and the used scaling matrix
    :returns: The scaling matrix used.
    """
    if not isinstance(structure, Structure):
        raise TypeError("Structure passed has to be a pymatgen structure")
    try:
        distance = float(distance)
        assert distance>1e-12, "Non-positive number"
    except Exception as e:
        print ("You have to pass positive float or integer as distance")
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
        if diagonal:
            lattice_cellvecs = np.array(lattice_cellvecs)
            # I get the diagonal solutions
            scale_matrix, supercell_cellvecs = get_diagonal_solution_bec(lattice_cellvecs, distance)
        else:
            # I get all possible midpoint vectors, based on the distance,
            # which for BEC method is the diameter of the sphere
            (norms_of_sorted_Gr_r2, sorted_Gc_r2, sorted_Gr_r2, r_outer, 
                v_diag ) = get_possible_solutions(lattice_cellvecs, distance,
                                verbosity=verbosity)
            if verbosity:
                print( "I received {} possible solutions".format(len(norms_of_sorted_Gr_r2)))
            # I pass these trial vectors into the function to find the minimum volume:
            if implementation == 'pyth':
                scale_matrix, supercell_cellvecs = get_optimal_solution_bec(
                    norms_of_sorted_Gr_r2, sorted_Gc_r2, sorted_Gr_r2, r_outer,
                    v_diag, r_inner=distance, verbosity=verbosity)
            elif implementation == 'fort':
                scale_matrix, supercell_cellvecs = fort_optimal_supercell_bec(norms_of_sorted_Gr_r2, sorted_Gc_r2, sorted_Gr_r2, r_outer,
                    v_diag, distance, verbosity,len(norms_of_sorted_Gr_r2))
            else:
                raise RuntimeError("Implementation {}".formt(implementation))
    elif method == 'hnf':
        if diagonal:
            lattice_cellvecs = np.array(lattice_cellvecs)
            scale_matrix, supercell_cellvecs = get_diagonal_solution_hnf(lattice_cellvecs, distance)
        else:
            if implementation == 'pyth':
                scale_matrix, supercell_cellvecs = get_optimal_solution_hnf(lattice_cellvecs, distance, verbosity)
            elif implementation == 'fort':
                scale_matrix, supercell_cellvecs = fort_optimal_supercell_hnf(lattice_cellvecs, distance, verbosity)
            else:
                raise RuntimeError("Implementation {}".formt(implementation))
    
        #raise NotImplementedError("HNF has not been fully implemented")
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
        print("Given Scaling:\n")
        print(scale_matrix)
        print("Given lattice:\n")
        print (new_lattice)
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

