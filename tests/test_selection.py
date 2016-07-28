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

class TestGenericSelection(unittest.TestCase):

    def setUp(self):
        self.member = 'a'
        self.selection_member = mastsel.SelectionMember(self.member)
        self.container = [self.selection_member]
        self.generic_selection = mastsel.GenericSelection(self.container)

    def tearDown(self):
        pass

    def test_constructor(self):
        with self.assertRaises(AssertionError):
            mastsel.GenericSelection(['a'])
            mastsel.GenericSelection([])


class TestSelection(unittest.TestCase):

    def setUp(self):
       self.members = ['a', 'b', 'c']
       self.selection_members = [mastsel.SelectionMember(sel) for sel in self.members]
       self.sel_idxs = [0,2]
       self.selection = mastsel.Selection(self.selection_members, self.sel_idxs)

    def tearDown(self):
        pass

    def test_slicing(self):
        # the second element is the third from members
        self.assertEquals(self.selection[1], self.selection_members[self.sel_idxs[1]])


    def test_register_selection(self):
        # simply make sure the selection was registered in each member
        for member in self.selection:
            self.assertIn(self.selection, member.get_selections())
            # and that we can retrieve it
            for key, selection in member.registry:
                self.assertEquals(member, selection[key])

class TestChainedSelection(unittest.TestCase):
     # set up Selection -0-> [Selection -0-> [SelectionMember]]
    def setUp(self):
       self.members = ['a', 'b', 'c']
       self.selection_members = [mastsel.SelectionMember(sel) for sel in self.members]
       self.selection = mastsel.Selection(self.selection_members, [0])
       self.selection_container = [self.selection]
       self.meta_selection = mastsel.Selection(self.selection_container, [0])

    def test_chained_registry_assignment(self):
        self.assertIn(self.meta_selection, self.selection_members[0].get_selections())
        self.assertNotIn(self.meta_selection, self.selection_members[1].get_selections())
        

if __name__ == "__main__":
    # doctests
    print("\n\n\n Doc Tests\n-----------")
    nfail, ntests = doctest.testmod(selection, verbose=True)

    # unit tests
    print("\n\n\n Unit Tests\n-----------")
    unittest.main()
