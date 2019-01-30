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

        np.random.seed(4)

        lattice = Lattice( 1*np.eye(3) + 0.5*(np.random.random((3,3))-0.5))
        sites = []
        for pos in itertools.product([-0.5,0.5], repeat=3):
            sites.append(PeriodicSite("H", pos, lattice, coords_are_cartesian=True))

        structure = Structure.from_sites(sites)

        supercell1, scale1 = make_supercell(structure, r_inner=3,
                verbosity=1, wrap=True, standardize=True, do_niggli_first=False)
        supercell2, scale2 = make_supercell(structure, r_inner=3,
                verbosity=1, wrap=True, standardize=True, do_niggli_first=True)

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

    #~ @unittest.skipIf(True,"")
    def test_structures(self):
        import itertools
        import numpy as np

        from pymatgen.core.structure import Structure
        from pymatgen.core.lattice import Lattice
        from pymatgen.core.sites import PeriodicSite

        from supercellor.supercell import make_supercell

        NTESTS = 100 #100
        RADIUS = 100.0
        EPS = 1e-3
        DIAG = 2
        NOISE_R = 1
        tests_run = 0
        np.random.seed(10)
        while (tests_run < NTESTS):
            R = DIAG*np.eye(3, dtype=int) - np.random.randint(NOISE_R, size=(3,3)) + NOISE_R/2
            #np.random.randint(10, size=(3,3)) -5 
            S = RADIUS* np.eye(3)
            try:
                P = np.dot(np.linalg.inv(R), S)
            except np.linalg.LinAlgError:
                continue

            lattice = Lattice(P)
            if lattice.volume < 0.01*RADIUS**3/DIAG:
                print 'skipping', lattice.volume
                continue
            sites = []
            try:
                for pos in itertools.product([-0.5,0.5], repeat=3):
                    sites.append(PeriodicSite("H", pos, lattice, coords_are_cartesian=True))
            except np.linalg.LinAlgError:
                continue
            structure = Structure.from_sites(sites)
            supercell, scale = make_supercell(structure, r_inner=RADIUS-EPS,
                    verbosity=0, wrap=True, standardize=True, do_niggli_first=True)
            self.assertTrue(np.sum(np.abs(supercell._lattice.matrix - S))< EPS)
            tests_run += 1

    def test_equivalence_fort_py(self):
        from datetime import datetime
        import itertools
        import numpy as np

        from pymatgen.core.structure import Structure
        from pymatgen.core.lattice import Lattice
        from pymatgen.core.sites import PeriodicSite

        from supercellor.supercell import make_supercell



        # tests per radius
        NTESTS = 10 #100
        VERBOSITY = 0
        for radius in np.arange(1.0, 3.1, 1.0):
            #RADIUS = 5.0
            tests_run = 0
            np.random.seed(10) # reinitialize seed!
            timings_f = 0.0
            timings_p = 0.0
            while (tests_run < NTESTS):
                P = np.eye(3)+ 0.1*(np.random.random((3,3)) -0.5)
                lattice = Lattice(P)
                if lattice.volume < 0.1:
                    print 'skipping', lattice.volume
                    continue
                sites = []
                try:
                    for pos in itertools.product([-0.5,0.5], repeat=3):
                        sites.append(PeriodicSite("H", pos, lattice, coords_are_cartesian=True))
                except np.linalg.LinAlgError:
                    continue
                structure = Structure.from_sites(sites)
                n = datetime.now()
                supercell_p, scale_p = make_supercell(structure, r_inner=radius,
                        verbosity=VERBOSITY, wrap=True, standardize=True, do_niggli_first=True,
                        implementation='pyth')
                timings_p += (datetime.now()-n).microseconds
                n = datetime.now()
                supercell_f, scale_f = make_supercell(structure, r_inner=radius,
                        verbosity=VERBOSITY, wrap=True, standardize=True, do_niggli_first=True,
                        implementation='fort')
                timings_f += (datetime.now()-n).microseconds
                self.assertTrue(np.sum((scale_p-scale_f)**2) ==0)
                self.assertTrue(np.sum((supercell_f._lattice.matrix -supercell_p._lattice.matrix)**2) < 1e-6)

                #~ self.assertTrue(np.sum(np.abs(supercell._lattice.matrix - S))< EPS)
                tests_run += 1
            print('Avg timing fortran impl rad={} {:.2e}'.format(radius, 1e-6*timings_f/ tests_run))
            print('Avg timing python  impl rad={} {:.2e}'.format(radius, 1e-6*timings_p/ tests_run))

if __name__ == '__main__':
    unittest.main()
