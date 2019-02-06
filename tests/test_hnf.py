import unittest

class Test(unittest.TestCase):

    @unittest.skipIf(True,"")
    def test_methods_compatibility(self):
        from pymatgen.core.structure import Structure
        from pymatgen.core.lattice import Lattice
        from pymatgen.core.sites import PeriodicSite
        from supercellor.supercell import make_supercell
        #from supercellor.lib.optimal_supercell import utils
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

        
        supercell1, scale1 = make_supercell(structure, r_inner=10, method='hnf', verbosity=1, wrap=True, standardize=True, do_niggli_first=False)
        supercell2, scale2 = make_supercell(structure, r_inner=10, method='bec', verbosity=1, wrap=True, standardize=True, do_niggli_first=False)
        self.assertTrue(supercell1._lattice.volume <= supercell2._lattice.volume)


    #~ @unittest.skipIf(True,"")
    def test_hnf_dmpi(self):
        from pymatgen.core.structure import Structure
        from pymatgen.core.lattice import Lattice
        from pymatgen.core.sites import PeriodicSite
        from supercellor.supercell import make_supercell
        #from supercellor.lib.optimal_supercell import utils
        import json, numpy as np, itertools, os
        np.random.seed(50)
        print '\n'
        for trial in range(10):
            lattice = Lattice( 1.0*np.eye(3) + 0.2*(np.random.random((3,3))-0.5))
            sites = []
            for pos in itertools.product([-0.5,0.5], repeat=3):
                sites.append(PeriodicSite("H", pos, lattice, coords_are_cartesian=True))
            structure = Structure.from_sites(sites)

            for rad in np.linspace(3, 10, 4):
                supercell1, scale1 = make_supercell(structure, r_inner=rad, method='hnf', verbosity=0, wrap=True, standardize=False, do_niggli_first=False)
                reduced_supercell1 = supercell1.get_reduced_structure(reduction_algo=u'LLL')
                print 'distances not red.:', np.linalg.norm(supercell1._lattice.matrix, axis=1)
                print 'angles not red.:', supercell1._lattice.angles
                print 'distances reduced :', np.linalg.norm(reduced_supercell1._lattice.matrix, axis=1)
                print 'angles reduced:', reduced_supercell1._lattice.angles
                for dim in range(3):
                    self.assertTrue(np.linalg.norm(reduced_supercell1._lattice.matrix[dim]) >= rad)



if __name__ == '__main__':
    unittest.main()
