import doctest
import unittest

import mast.selection as mastsel
import mast.molecule as mastmol
import mast.system as mastsys
import mast.interactions as mastinx

import mast.tests.data as mastdata

import mast.config.molecule as mastmolconfig
import mast.config.system as mastsysconfig
import mast.config.interactions as mastinxconfig



if __name__ == "__main__":

s    from mast import interactions

    # doctests
    print("\n\n\n Doc Tests\n-----------")
    nfail, ntests = doctest.testmod(interactions, verbose=True)

    # unit tests
    print("\n\n\n Unit Tests\n-----------")
    unittest.main()

