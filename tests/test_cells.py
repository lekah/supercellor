import unittest

class Test(unittest.TestCase):

    @unittest.skipIf(True,"")
    def test_niggli(self):
        from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
        from pymatgen.core.structure import Structure
        from pymatgen.core.lattice import Lattice
        from pymatgen.core.sites import PeriodicSite
        from supercellor.supercell import make_supercell
        import numpy as np
        import itertools

        np.random.seed(4)
        DIST = 0.99

        lattice = Lattice( 1*np.eye(3) + 0.5*(np.random.random((3,3))-0.5))
        sites = []
        for pos in itertools.product([-0.5,0.5], repeat=3):
            sites.append(PeriodicSite("H", pos, lattice, coords_are_cartesian=True))

        structure = Structure.from_sites(sites)

        supercell1, scale1 = make_supercell(structure, distance=DIST,
                verbosity=1, wrap=True, standardize=True, do_niggli_first=False)
        supercell2, scale2 = make_supercell(structure, distance=DIST,
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

    @unittest.skipIf(True,"")
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
        DIAG = 4
        NOISE_R = 3
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
            supercell, scale = make_supercell(structure, distance=RADIUS-EPS,
                    verbosity=0, wrap=True, standardize=True, do_niggli_first=True)
            self.assertTrue(np.sum(np.abs(supercell._lattice.matrix - S))< EPS)
            tests_run += 1

    # ~ @unittest.skipIf(True,"")
    def test_supercells(self):
        from pymatgen.io.cif import CifParser
        import numpy as np
        from supercellor.supercell import make_supercell

        EPS = 1e-6

        def get_shortest_distance(system):
            nat = len(system.sites)
            shortest_distance = np.inf
            for i in range(nat):
                for j in range(i+1, nat):
                    d = system.get_distance(i,j)
                    shortest_distance = min((d, shortest_distance))
            return shortest_distance


        filenames = ('data/Al-fcc.cif', )
        for fname in filenames:
            parser = CifParser(fname)
            structure = parser.get_structures()[0]
            shortestd_unit = get_shortest_distance(2*structure)
            for method in ( 'bec', 'hnf',):
                for wrap in (False, True):
                    for standardize in (False, True):
                        for do_niggli in (False, True):
                            for D in np.arange(5.0, 15.1, 5.0):
                                for implementation in ('fort', 'pyth'):
                                    if D >= 10 and implementation == 'pyth':
                                        continue
                                    S1, sm = make_supercell(structure, distance=D,
                                        method=method, verbosity=0, implementation='fort', wrap=wrap,
                                        standardize=standardize, do_niggli_first=do_niggli)
                                    S2 = structure*sm
                                    shortestd_super = get_shortest_distance(S1)
                                    # print method, wrap, standardize, do_niggli, D, implementation, shortestd_super, shortestd_unit
                                    # Checking that atoms are not moved closer together, for any reasons
                                    self.assertTrue(abs(shortestd_super-shortestd_unit)<EPS)
                                    if not standardize:
                                        angl1 = np.array(S1._lattice.angles)
                                        angl2 = np.array(S2._lattice.angles)
                                        angl1.sort() # sort in place
                                        angl2.sort()
                                        # checking angles are the same
                                        self.assertTrue(np.sum((angl1-angl2)**2) < EPS)
                                    if not (standardize or do_niggli):
                                        # check consistency between scaling matrix and lattice vecs:
                                        self.assertTrue(np.sum((S1._lattice.matrix - S2._lattice.matrix)**2) < 1e-6)


if __name__ == '__main__':
    unittest.main()
