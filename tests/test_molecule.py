import doctest
import unittest

import numpy as np
import numpy.testing as npt

from mast import molecule
import mast.selection as mastsel

class TestMolecule(unittest.TestCase):
    def setUp(self):
        pass
    def tearDown(self):
        pass


if __name__ == "__main__":

    from mast import molecule

    # doctests
    print("\n\n\n Doc Tests\n-----------")
    nfail, ntests = doctest.testmod(molecule, verbose=True)

    # unit tests
    print("\n\n\n Unit Tests\n-----------")
    unittest.main()
