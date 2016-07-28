import doctest
import unittest

from mast import selection
import mast.selection as mastsel

class TestSelectionMember(unittest.TestCase):
    def setUp(self):
        self.member = 'a'
        self.selection_member = mastsel.SelectionMember(self.member)

    def tearDown(self):
        pass

    def test_constructor(self):
        pass

    def test_member(self):
        self.assertEqual(self.selection_member.member, self.member)

    def test_registry(self):
        self.assertEqual(self.selection_member.registry, [])
        sel = mastsel.Selection([self.selection_member], [0])
        self.assertIn((0, sel), self.selection_member.registry)

    def test_repr(self):
        pass

    def test_get_selections(self):
        sel = mastsel.Selection([self.selection_member], [0])
        self.assertIn(sel, self.selection_member.get_selections())

    def test_register_selection(self):
        pass



if __name__ == "__main__":
    # doctests
    print("\n\n\n Doc Tests\n-----------")
    nfail, ntests = doctest.testmod(selection, verbose=True)

    # unit tests
    print("\n\n\n Unit Tests\n-----------")
    unittest.main()
