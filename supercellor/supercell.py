
import numpy as np
from time import sleep
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse
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


def make_supercell(structure, min_image_distance, wrap=True, standardize=True,
        do_niggli_first=True, verbosity=1):
    """
    Creates from a given structure a supercell based on the required minimal dimension
    :param structure: The pymatgen structure to create the supercell for
    :param float min_image_distance: The minimum image distance as a float,
        The cell created will not have any periodic image below this distance
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
        min_image_distance = float(min_image_distance)
        assert min_image_distance>1e-12, "Non-positive number"
    except Exception as e:
        print "You have to pass positive float or integer as min_image_distance"
        raise e

    if not isinstance(wrap, bool):
        raise TypeError("wrap has to be a boolean")
    if not isinstance(standardize, bool):
        raise TypeError("standardize has to be a boolean")

    # I'm getting the niggli reduced structure as first:
    if verbosity > 1:
        print "given cell:"
        print structure._lattice
    if do_niggli_first:
        starting_structure = structure.get_reduced_structure(reduction_algo=u'niggli')
    else:
        starting_structure = structure
    if verbosity > 1:
        print "starting cell:"
        print starting_structure._lattice
        for i, v in enumerate(starting_structure._lattice.matrix):
            print i, np.linalg.norm(v)
    
    # the lattice of the niggle reduced structure:
    lattice_cellvecs = starting_structure._lattice.matrix
    # trial_vecs are all possible vectors sorted by the norm
    trial_vecs = get_trial_vecs(lattice_cellvecs, min_image_distance, verbosity=verbosity)
    # I pass these trial vectors into the function to find the minimum volume:
    scale_matrix, supercell_cellvecs = find_min_vol(trial_vecs, lattice_cellvecs, min_image_distance, 
            verbosity=verbosity)

    # Constructing the new lattice:
    new_lattice = Lattice(supercell_cellvecs)
    # I create f_lat, which are the fractional lattice points of the niggle_reduced:
    f_lat = lattice_points_in_supercell(scale_matrix)
    # and transforrm to cartesian coords here:
    c_lat = new_lattice.get_cartesian_coords(f_lat)
    #~ cellT = supercell_cellvecs.T

    if verbosity > 1:
        print "Given lattice:"
        print new_lattice
        for i, v in enumerate(new_lattice.matrix):
            print i, np.linalg.norm(v)

    new_sites = []
    if verbosity:
        print "Done, constructing structure"
    for site in starting_structure:
        for v in c_lat:
            new_sites.append(PeriodicSite(site.species_and_occu, site.coords +v,
                             new_lattice, properties=site.properties,
                             coords_are_cartesian=True, to_unit_cell=wrap))

    supercell = Structure.from_sites(new_sites)

    if standardize:
        frac_coords = np.empty((len(supercell), 3))
        for i, site in enumerate(supercell):
            frac_coords[i,:] = site.frac_coords

        # The code below sorts the cell and positions by lattice vector length,
        # shortest first:
        veclens = sorted([(np.linalg.norm(v), i) for i, v in enumerate(supercell_cellvecs)])
        M1 = np.zeros((3,3))
        for row, (_, idx) in enumerate(veclens):
            M1[idx, row] = 1
        supercell_cellvecs, frac_coords = rotate(M1, supercell_cellvecs, frac_coords)

        k_a, k_b, k_c = get_angles(supercell_cellvecs)
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
        supercell_cellvecs, frac_coords = rotate(M2, supercell_cellvecs, frac_coords)
        # Now applying the rules layed out in  http://arxiv.org/abs/1506.01455
        # to get the standardized triclinic cell.
        # Since not all cells give to me are triclinic < -> <= with respect to
        # the paper. So it's not truly standardized!
        metric = np.dot(supercell_cellvecs, supercell_cellvecs.T)

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

        supercell_cellvecs = np.zeros((3,3))
        supercell_cellvecs[0, 0] = a
        supercell_cellvecs[1, 0] = b * cg
        supercell_cellvecs[2, 0] = c * cb
        supercell_cellvecs[1, 1] = b * sg
        supercell_cellvecs[2, 1] = c * (ca - cb * cg) / sg
        supercell_cellvecs[2, 2] = c * np.sqrt(1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg) / sg
        
        # And some checks:
        if (
            supercell_cellvecs[0,0] < -1e-12 or
            abs(supercell_cellvecs[0,1]) > 1e-12 or
            abs(supercell_cellvecs[0,2]) > 1e-12):
            raise ValueError("First lattice vector not aligned with x-axis")
        if  (
            supercell_cellvecs[1,1] < -1e-12 or
            abs(supercell_cellvecs[1,2]) > 1e-12):
            raise ValueError("Second lattice vector not in X-Y plane in first quadrant")
        if  (supercell_cellvecs[2,2] < 0):
            raise ValueError("Third lattice vector not in positive Z-direction")
        
        new_sites = []
        new_lattice = Lattice(supercell_cellvecs)

        for i, site in enumerate(supercell):
            s = PeriodicSite(site.species_and_occu, frac_coords[i],
                            new_lattice, properties=site.properties,
                            coords_are_cartesian=False, to_unit_cell=wrap)
            new_sites.append(s)
        
        supercell = Structure.from_sites(new_sites)

        #~ view(ad.get_atoms(supercell))
        if verbosity > 1:
            print "Cell after standardization:"
            print new_lattice
            for i, v in enumerate(new_lattice.matrix):
                print i, np.linalg.norm(v)

    return supercell, scale_matrix
    



def get_trial_vecs(cell, min_image_distance, verbosity=1):
    """
    Returns all possible vectors (as integer crystal coordinates)
    that lie outside a sphere of radius min_dist_soll, but within a sphere given by a
    maximum radius that is determined by the upper bound, the diagonal solution.
    The trial vectors could be iteratively reduced with better solutions found, while still being exhaustive
    """
    # Would be good to first LLL/Niggle reduce to make the upper bound as low as possible.
    
    #atoms = structure.get_ase()
    a,b,c = [np.array(v) for v in cell]
    repetitions =[None]*3
    for index, item in enumerate(((a,b,c), (b,c,a), (c,a,b))):
        n1, n2, n3 = item
        d23 = np.cross(n2, n3) 
        d1 = np.linalg.norm(np.dot(n1, d23)) / np.linalg.norm(d23)
        repetitions[index] = int(math.ceil(min_image_distance/d1))
    maxcell = np.array([rep*v for rep, v in zip(repetitions, cell)])
    diagvol = np.abs(np.dot(np.cross(maxcell[0], maxcell[1]), maxcell[2]))
    if verbosity:
        print "Volume of diagonal supercell:", diagvol
    maxradius =  diagvol /  min_image_distance**2
    # Now I find all vectors that lie within the maxradius but outside the minradius
    # Inversion symmetry allows me to look only in the upper half.
    trials = []
    for ia in range(0, repetitions[0]+1):
    #for ia in range(-repetitions[0], repetitions[0]+1):
        for ib in range(-repetitions[1], repetitions[1]+1):
        #for ib in range(0, repetitions[1]+1):
            for ic in range(-repetitions[2], repetitions[2]+1):
            #for ic in range(0, repetitions[2]+1):
                vector = ia*cell[0]+ib*cell[1]+ic*cell[2]
                #print vector, np.linalg.norm(vector)
                veclen = np.linalg.norm(vector)
                if min_image_distance <= veclen <= maxradius:
                    trials.append((veclen, ia, ib, ic, vector))
    trials = sorted(trials)
    return trials

def find_min_vol(trial_vecs, cell, min_image_distance, verbosity=1):
    min_volume = np.inf
    max_min_per_image_distance = 0
    max_radius = trial_vecs[-1][0]
    for i1, (veclen1, ia1, ib1, ic1, vector1) in enumerate(trial_vecs, start=0):
        
        if veclen1 > max_radius:
            break
        for i2, (veclen2, ia2, ib2, ic2, vector2) in enumerate(trial_vecs[i1+1:], start=i1+1):
            if veclen2 > max_radius:
                break
            if np.abs(np.dot(vector1, vector2)) / (veclen1*veclen2) >= 0.5:
                continue
            for i3, (veclen3, ia3, ib3, ic3, vector3) in enumerate(trial_vecs[i2+1:], start=i2+1):
                if veclen3 > max_radius:
                    break
                if np.abs(np.dot(vector2, vector3)) / (veclen2*veclen3) >= 0.5:
                    continue
                elif np.abs(np.dot(vector1, vector3)) / (veclen1*veclen3) >= 0.5:
                    continue
                # checking intersections of each plane
                d1 = np.abs(np.dot(np.cross(vector2/veclen2, vector3/veclen3), vector1))
                if d1 < min_image_distance:
                    continue
                d2 = np.abs(np.dot(np.cross(vector1/veclen1, vector3/veclen3), vector2)) 
                if d2 < min_image_distance:
                    continue
                d3 = np.abs(np.dot(np.cross(vector1/veclen1, vector2/veclen2), vector3))
                if d3 < min_image_distance:
                    continue
                volume = np.abs(np.dot(np.cross(vector1,vector2), vector3))
                min_per_image_distance = min((d1,d2, d3))
                if volume < min_volume:
                    min_volume = volume
                    max_radius = min_volume/min_image_distance**2
                    max_min_per_image_distance = min_per_image_distance
                    is_best = True
                elif volume == min_volume and min_per_image_distance > max_min_per_image_distance:
                    max_min_per_image_distance = min_per_image_distance
                    is_best = True
                else:
                    is_best = False
                if is_best:
                    if verbosity:
                        print "New optimal supercell {} & {}".format(volume, max_min_per_image_distance)
                    chosen_coords = np.array(((ia1,ib1,ic1),(ia2,ib2,ic2),(ia3,ib3,ic3)))
                    cellvecs = np.array([vector1, vector2, vector3])
    if verbosity:
        print 'FINAL'
        print chosen_coords
    return np.array(chosen_coords, dtype=int), np.array(cellvecs)

