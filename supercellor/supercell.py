
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
from scipy.linalg import expm


def make_supercell(structure, min_image_distance, wrap=True, refine=True):
    """
    Creates from a given structure a supercell based on the required minimal dimension
    :param structure: The pymatgen structure to create the supercell for
    :param float min_image_distance: The minimum image distance as a float,
        The cell created will not have any periodic image below this distance
    :param bool wrap: Wrap the atoms into the created cell after replication
    :param bool refine: Use spglib to refine the created cell, use with care for big structures

    :returns: A new pymatgen core structure instance and the used scaling matrix
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
    if not isinstance(refine, bool):
        raise TypeError("wrap has to be a boolean")

    # I'm getting the niggle reduced structure as a first:
    niggli_reduced = structure.get_reduced_structure(reduction_algo=u'niggli')
    # the lattice of the niggle reduced structure:
    lattice_cellvecs = niggli_reduced._lattice.matrix
    # trial_vecs are all possible vectors sorted by the norm
    trial_vecs = get_trial_vecs(lattice_cellvecs, min_image_distance)
    # I pass these trial vectors into the function to find the minimum volume:
    scale_matrix, supercell_cellvecs = find_min_vol(trial_vecs, lattice_cellvecs, min_image_distance)

    # Constructing the new lattice:
    new_lattice = Lattice(supercell_cellvecs)
    # I create f_lat, which are the fractional lattice points of the niggle_reduced:
    f_lat = lattice_points_in_supercell(scale_matrix)
    # and transforrm to cartesian coords here:
    c_lat = new_lattice.get_cartesian_coords(f_lat)
    #~ cellT = supercell_cellvecs.T
    #~ cellI = np.array(np.matrix(cellT).I)

    new_sites = []
    print "Done, constructing structure"
    for site in niggli_reduced:
        for v in c_lat:
            pos = site.coords + v
            if wrap:
                pos = np.dot(cellT, np.dot(cellI, pos) % 1.0 )
            s = PeriodicSite(site.species_and_occu, pos,
                             new_lattice, properties=site.properties,
                             coords_are_cartesian=True, to_unit_cell=False)
            new_sites.append(s)

    supercell = Structure.from_sites(new_sites)
    #~ view(ad.get_atoms(supercell))
    if refine:
        from spglib import refine_cell
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        print "Refining cell using spglib"
        sa = SpacegroupAnalyzer(supercell, symprec=1e-21, angle_tolerance=-1)
        oldvolume = supercell._lattice.volume
        supercell = sa.get_refined_structure()
        newvolume = supercell._lattice.volume
        assert abs(oldvolume - newvolume) / newvolume < 1e-10, "Lost volume in refinement {:.5f} {:.5f}".format(oldvolume, newvolume)
    return supercell, scale_matrix


    



def get_trial_vecs(cell, min_image_distance):
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

def find_min_vol(trial_vecs, cell, min_image_distance):
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
                    print "New optimal supercell {} & {}".format(volume, max_min_per_image_distance)
                    chosen_coords = np.array(((ia1,ib1,ic1),(ia2,ib2,ic2),(ia3,ib3,ic3)))
                    cellvecs = np.array([vector1, vector2, vector3])
    print 'FINAL'
    print chosen_coords
    return np.array(chosen_coords, dtype=int), np.array(cellvecs)
    

if __name__ == '__main__':
    from pymatgen.io.cif import CifParser
    #~ parser = CifParser("LiGePS.cif")
    parser = CifParser("MyBaseFileNameCollCode421725.cif")
    structure = parser.get_structures()[0]
    print 'OLD VOLUME', structure._lattice.volume
    #~ structure._lattice = Lattice(5*np.eye(3) + np.random.random((3,3)))
    #~ structure._lattice = Lattice()
    radius=20
    make_supercell(structure, radius, wrap=False, refine=True)
    #~ structure
    #~ radius = 10
    #~ cell = np.array(((1, 0.5, -0.1) , (0, 1,-0.2), (00.1,-0.3,1))) 
    #~ trial_vecs = get_trial_vecs(cell, radius)
    #~ print "Number of trial vectors:", len(trial_vecs)
    #for _, i, j, k, _ in trial_vecs:
    #
    #~ print i, j, k
