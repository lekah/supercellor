import unittest

class Test(unittest.TestCase):
    def test_niggli(self):
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        from pymatgen.core.structure import Structure
        from pymatgen.core.lattice import Lattice
        from pymatgen.core.sites import PeriodicSite
        from supercellor.supercell import make_supercell
        import numpy as np
        import itertools
        return
        lattice = Lattice( 5*(np.random.random((3,3))-0.5))
        
        sites = []
        for pos in itertools.product([-0.5,0.5], repeat=3):
            sites.append(PeriodicSite("H", pos, lattice, coords_are_cartesian=True))

        structure = Structure.from_sites(sites)

        #~ with open('structure.json','w') as f:
            #~ json.dump(structure.as_dict(), f)
        
        supercell1, scale1 = make_supercell(structure, min_image_distance=3,
                verbosity=2, wrap=True, standardize=True, do_niggli_first=False)
        supercell2, scale2 = make_supercell(structure, min_image_distance=3,
                verbosity=2, wrap=True, standardize=True, do_niggli_first=True)

        sa = SpacegroupAnalyzer(supercell1, symprec=1e-21, angle_tolerance=-1)
        supercell_refine = sa.get_refined_structure()
        for i in range(3):
            for j in range(3):
                self.assertTrue(abs(supercell1._lattice.matrix[i,j] - supercell2._lattice.matrix[i,j]) < 1e-1, 
                    '{} != {} for i,j={},{}'.format(
                        supercell1._lattice.matrix[i,j],
                        supercell2._lattice.matrix[i,j], i,j))
                self.assertTrue(abs(supercell1._lattice.matrix[i,j] - supercell_refine._lattice.matrix[i,j]) < 1e-1,
                    '{} != {} for i,j={},{}'.format(
                        supercell1._lattice.matrix[i,j],
                        supercell_refine._lattice.matrix[i,j],i,j))
    def test_methods_compatibility(self):
        from pymatgen.core.structure import Structure
        from pymatgen.core.lattice import Lattice
        from pymatgen.core.sites import PeriodicSite
        from supercellor.supercell import make_supercell
        from supercellor.lib.optimal_supercell import utils
        import json, numpy as np, itertools, os
        create_new = False
        if create_new:
            lattice = Lattice( 1.1*np.eye(3) + (np.random.random((3,3))-0.5))
            sites = []
            for pos in itertools.product([-0.5,0.5], repeat=3):
                sites.append(PeriodicSite("H", pos, lattice, coords_are_cartesian=True))


            structure = Structure.from_sites(sites)
            counter = 1
            while True:
                fname = 'data/structure-{}.json'.format(counter)
                if os.path.exists(fname):
                    counter += 1
                else:
                    break
                if counter > 100:
                    raise Exception("Reached 100 files")
            with open(fname,'w') as f:
                json.dump(structure.as_dict(), f)
        else:
            with open('data/structure-1.json', 'r') as f:
                d = json.load(f)
                structure = Structure.from_dict(d)
            print structure._lattice

        supercell1, scale1 = make_supercell(structure, min_image_distance=3, method='hnf', verbosity=2, wrap=True, standardize=False, do_niggli_first=False)
        det1 = utils.determinant33_int(scale1)
        supercell2, scale2 = make_supercell(structure, min_image_distance=3, method='bfv', verbosity=1, wrap=True, standardize=False, do_niggli_first=False)
        det2 = utils.determinant33_int(scale2)
        print det1, det2
        print 'LATTICE HNF'
        print supercell1._lattice
        print 'LATTICE BFV'
        print supercell2._lattice
        self.assertEqual(abs(det1), abs(det2))


if __name__ == '__main__':
    unittest.main()
