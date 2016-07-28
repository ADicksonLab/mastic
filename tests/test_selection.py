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


    def test_unselected_registry(self):
        self.assertEqual(self.selection_member.registry, [])

    def test_repr(self):
        pass

    def test_get_selections(self):
        sel = mastsel.Selection([self.selection_member], [0])
        # recursive
        self.assertIn(sel, self.selection_member.get_selections())
        # non-recursive
        self.assertIn(sel, self.selection_member.get_selections(level=0))

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

    def test_getitem(self):
        # the second element is the third from members
        self.assertEqual(self.selection[1], self.selection_members[self.sel_idxs[1]])

    def test_selection_member_self_retrieval(self):
        for sel_memb in self.selection_members:
            for key, selection in sel_memb.registry:
                self.assertEqual(selection[key], sel_memb)


class TestChainedSelection(unittest.TestCase):
     # set up Selection -0-> [Selection -0-> [SelectionMember]]
    def setUp(self):
       self.members = ['a', 'b', 'c']
       self.selection_members = [mastsel.SelectionMember(sel) for sel in self.members]
       self.selection = mastsel.Selection(self.selection_members, [0])
       self.selection_container = [self.selection]
       self.meta_selection = mastsel.Selection(self.selection_container, [0])

    def test_chained_registry_assignment(self):
        # the first level is in the recursive get
        self.assertIn(self.selection, self.selection_members[0].get_selections())
        # the first level is in the level=0 get
        self.assertIn(self.selection, self.selection_members[0].get_selections(level=0))
        # specifying too many levels is ignored silently
        self.assertIn(self.selection, self.selection_members[0].get_selections(level=3))

        # the second level selection is in the recursive get
        self.assertIn(self.meta_selection, self.selection_members[0].get_selections())
        # the second level selection is not in the level=0 get
        self.assertNotIn(self.meta_selection, self.selection_members[0].get_selections(level=0))
        # the second level selection is in the level=1 get
        self.assertIn(self.meta_selection, self.selection_members[0].get_selections(level=1))

        for unselected_member in self.selection_members[1:]:
            # neither selection is in any other SelectionMember
            self.assertEqual(unselected_member.get_selections(), [])

class TestIndexedSelection(unittest.TestCase):

    def setUp(self):
       self.members = ['a', 'b', 'c']
       self.selection_members = [mastsel.SelectionMember(sel) for sel in self.members]
       self.selection_idxs = [0, 2]
       self.idx_selection = mastsel.IndexedSelection(self.selection_members, self.selection_idxs)

    def tearDown(self):
        pass

    def test_getitem(self):
        for idx in self.selection_idxs:
            self.assertEqual(self.idx_selection[idx], self.selection_members[idx])

if __name__ == "__main__":
    # doctests
    print("\n\n\n Doc Tests\n-----------")
    nfail, ntests = doctest.testmod(selection, verbose=True)

    # unit tests
    print("\n\n\n Unit Tests\n-----------")
    unittest.main()
