import doctest
from mast import selection

nfail, ntests = doctest.testmod(selection, verbose=True)
