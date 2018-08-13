import unittest

class TestRefinement(unittest.TestCase):
    def test_refinement(self):
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        from pymatgen.core.structure import Structure
        from pymatgen.core.lattice import Lattice
        from pymatgen.core.sites import PeriodicSite
        from supercellor.supercell import make_supercell
        import numpy as np
        import itertools
        lattice = Lattice( 5*(np.random.random((3,3))-0.5))
        
        sites = []
        for pos in itertools.product([-0.5,0.5], repeat=3):
            sites.append(PeriodicSite("H", pos, lattice, coords_are_cartesian=True))

        structure = Structure.from_sites(sites)

        #~ with open('structure.json','w') as f:
            #~ json.dump(structure.as_dict(), f)
        #~ with open('structure.json', 'r') as f:
            #~ d = json.load(f)
            #~ structure = Structure.from_dict(d)

        supercell1, scale1 = make_supercell(structure, min_image_distance=3, verbosity=2, wrap=True, standardize=True, do_niggli_first=False)
        supercell2, scale2 = make_supercell(structure, min_image_distance=3, verbosity=2, wrap=True, standardize=True, do_niggli_first=True)
        #~ view(ad.get_atoms(supercell))
        sa = SpacegroupAnalyzer(supercell1, symprec=1e-21, angle_tolerance=-1)
        supercell_refine = sa.get_refined_structure()
        for i in range(3):
            for j in range(3):

                self.assertTrue(abs(supercell1._lattice.matrix[i,j] - supercell2._lattice.matrix[i,j]) < 1e-1, '{} != {} for i,j={},{}'.format(
                    supercell1._lattice.matrix[i,j], supercell2._lattice.matrix[i,j],i,j))
                self.assertTrue(abs(supercell1._lattice.matrix[i,j] - supercell_refine._lattice.matrix[i,j]) < 1e-1, '{} != {} for i,j={},{}'.format(
                    supercell1._lattice.matrix[i,j], supercell_refine._lattice.matrix[i,j],i,j))




if __name__ == '__main__':
    unittest.main()
