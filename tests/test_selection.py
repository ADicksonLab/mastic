import doctest
import unittest

import numpy as np
import numpy.testing as npt

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
        import ipdb; ipdb.set_trace()
        sel0 = mastsel.Selection([self.selection_member], [0])
        sel1 = mastsel.Selection([self.selection_member], [0], flags='other_selection')
        sel2 = mastsel.IndexedSelection([self.selection_member], [0])
        # 
        meta_selection = mastsel.Selection()
        meta_meta_selection = mastsel.Selection([meta_selection])
        sel_list = mastsel.SelectionsList([sel0, sel1, sel2], flags=['list-selection'])
        meta_sel_list = mastsel.Selection(sel_list, [0], flags=['meta-list-selection'])

        # level 0
        self.assertIn(sel0, self.selection_member.get_selections(level=0))
        self.assertIn(sel1, self.selection_member.get_selections(level=0))
        self.assertIn(sel2, self.selection_member.get_selections(level=0))
        self.assertNotIn(sel_list, self.selection_member.get_selections(level=0))
        self.assertNotIn(meta_sel_list, self.selection_member.get_selections(level=0))
        # level 1
        self.assertIn(sel0, self.selection_member.get_selections(level=1))
        self.assertIn(sel1, self.selection_member.get_selections(level=1))
        self.assertIn(sel2, self.selection_member.get_selections(level=1))
        self.assertIn(sel_list, self.selection_member.get_selections(level=1))
        self.assertNotIn(meta_sel_list, self.selection_member.get_selections(level=1))
        # level 2
        self.assertIn(sel0, self.selection_member.get_selections(level=2))
        self.assertIn(sel1, self.selection_member.get_selections(level=2))
        self.assertIn(sel2, self.selection_member.get_selections(level=2))
        self.assertIn(sel_list, self.selection_member.get_selections(level=2))
        self.assertIn(meta_sel_list, self.selection_member.get_selections(level=2))
        # recursive
        # explicit syntax
        self.assertIn(sel0, self.selection_member.get_selections(level=None))
        self.assertIn(sel1, self.selection_member.get_selections(level=None))
        self.assertIn(sel2, self.selection_member.get_selections(level=None))
        self.assertIn(sel_list, self.selection_member.get_selections(level=None))
        self.assertIn(meta_sel_list, self.selection_member.get_selections(level=None))
        # implicit syntax
        self.assertIn(sel0, self.selection_member.get_selections())
        self.assertIn(sel1, self.selection_member.get_selections())
        self.assertIn(sel2, self.selection_member.get_selections())
        self.assertIn(sel_list, self.selection_member.get_selections())
        self.assertIn(meta_sel_list, self.selection_member.get_selections())

        # with selection criteria
        # level 0 only Selections
        self.assertIn(sel0,
                      self.selection_member.get_selections(
                          selection_types=[mastsel.Selection],
                          level=0))
        self.assertIn(sel1,
                      self.selection_member.get_selections(
                          selection_types=[mastsel.Selection],
                          level=0))
        self.assertNotIn(sel2,
                      self.selection_member.get_selections(
                          selection_types=[mastsel.Selection],
                          level=0))
        self.assertNotIn(sel_list,
                      self.selection_member.get_selections(
                          selection_types=[mastsel.Selection],
                          level=0))
        self.assertNotIn(meta_sel_list,
                      self.selection_member.get_selections(
                          selection_types=[mastsel.Selection],
                          level=0))

        # level 0 only IndexedSelection
        self.assertNotIn(sel0,
                      self.selection_member.get_selections(
                          selection_types=[mastsel.IndexedSelection],
                          level=0))
        self.assertNotIn(sel1,
                      self.selection_member.get_selections(
                          selection_types=[mastsel.IndexedSelection],
                          level=0))
        self.assertIn(sel2,
                      self.selection_member.get_selections(
                          selection_types=[mastsel.IndexedSelection],
                          level=0))
        self.assertNotIn(sel_list,
                      self.selection_member.get_selections(
                          selection_types=[mastsel.IndexedSelection],
                          level=0))
        self.assertNotIn(meta_sel_list,
                      self.selection_member.get_selections(
                          selection_types=[mastsel.IndexedSelection],
                          level=0))

        # level 0 only 'other_selection' flag
        self.assertNotIn(sel0,
                      self.selection_member.get_selections(
                          flags=['other_selection'],
                          level=0))
        self.assertIn(sel1,
                      self.selection_member.get_selections(
                          flags=['other_selection'],
                          level=0))
        self.assertNotIn(sel2,
                      self.selection_member.get_selections(
                          flags=['other_selection'],
                          level=0))
        self.assertNotIn(sel_list,
                      self.selection_member.get_selections(
                          flags=['other_selection'],
                          level=0))
        self.assertNotIn(meta_sel_list,
                      self.selection_member.get_selections(
                          flags=['other_selection'],
                          level=0))

        # recursive only Selections
        self.assertIn(sel0,
                      self.selection_member.get_selections(
                          selection_types=[mastsel.Selection],
                          level=None))
        self.assertIn(sel1,
                      self.selection_member.get_selections(
                          selection_types=[mastsel.Selection],
                          level=None))
        self.assertNotIn(sel2,
                      self.selection_member.get_selections(
                          selection_types=[mastsel.Selection],
                          level=None))
        self.assertNotIn(sel_list,
                      self.selection_member.get_selections(
                          selection_types=[mastsel.Selection],
                          level=None))
        self.assertIn(meta_sel_list,
                      self.selection_member.get_selections(
                          selection_types=[mastsel.Selection],
                          level=None))

        # recursive only SelectionsList
        self.assertNotIn(sel0,
                      self.selection_member.get_selections(
                          selection_types=[mastsel.SelectionsList],
                          level=None))
        self.assertNotIn(sel1,
                      self.selection_member.get_selections(
                          selection_types=[mastsel.SelectionsList],
                          level=None))
        self.assertNotIn(sel2,
                      self.selection_member.get_selections(
                          selection_types=[mastsel.SelectionsList],
                          level=None))
        self.assertIn(sel_list,
                      self.selection_member.get_selections(
                          selection_types=[mastsel.SelectionsList],
                          level=None))
        self.assertNotIn(meta_sel_list,
                      self.selection_member.get_selections(
                          selection_types=[mastsel.SelectionsList],
                          level=None))

        # recursive only 'list-selection' flag
        self.assertIn(sel0,
                      self.selection_member.get_selections(
                          flags=['list-selection'],
                          level=None))
        self.assertIn(sel1,
                      self.selection_member.get_selections(
                          flags=['list-selection'],
                          level=None))
        self.assertIn(sel2,
                      self.selection_member.get_selections(
                          flags=['list-selection'],
                          level=None))
        self.assertIn(sel_list,
                      self.selection_member.get_selections(
                          flags=['list-selection'],
                          level=None))
        self.assertNotIn(meta_sel_list,
                      self.selection_member.get_selections(
                          flags=['list-selection'],
                          level=None))

        # recursive only 'meta-list-selection' flag
        self.assertIn(sel0,
                      self.selection_member.get_selections(
                          flags=['meta-list-selection'],
                          level=None))
        self.assertNotIn(sel1,
                      self.selection_member.get_selections(
                          flags=['meta-list-selection'],
                          level=None))
        self.assertNotIn(sel2,
                      self.selection_member.get_selections(
                          flags=['meta-list-selection'],
                          level=None))
        self.asserttIn(sel_list,
                      self.selection_member.get_selections(
                          flags=['meta-list-selection'],
                          level=None))
        self.assertIn(meta_sel_list,
                      self.selection_member.get_selections(
                          flags=['meta-list-selection'],
                          level=None))


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

    def test_constructor(self):
        with self.assertRaises(AssertionError):
            mastsel.Selection(self.selection_members, 'a')
            mastsel.Selection(self.selection_members, ['a','b'])
            mastsel.Selection(self.selection_members, -1)
            mastsel.Selection(self.selection_members, [-1,-2])
            mastsel.Selection(self.selection_members, [-1,2])
            mastsel.Selection(self.selection_members, [])
            mastsel.Selection(self.members, [1])

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

    def test_selection_member_self_retrieval(self):
        for sel_memb in self.selection_members:
            for key, selection in sel_memb.registry:
                self.assertEqual(selection[key], sel_memb)

class TestCoordArray(unittest.TestCase):

    def setUp(self):
        self.array = np.array([[0,0,0], [1,1,1], [2,2,2]])
        self.coords = mastsel.CoordArray(self.array)
        self.new_coord = np.array([3,3,3])

    def test_add_coord(self):
        target_array = np.array([[0,0,0], [1,1,1], [2,2,2], [3,3,3]])
        self.assertEqual(self.coords.add_coord(self.new_coord), 3)
        npt.assert_equal(self.coords.coords, target_array)
        with self.assertRaises(AssertionError):
            self.coords.add_coord(np.array([4,4,4,4]))
            self.coords.add_coord(np.array([2,2]))
            self.coords.add_coord(np.array([]))
            self.coords.add_coord([])
            self.coords.add_coord({'a', 1})

    def test_coord_setter(self):
        with self.assertRaises(AssertionError):
            self.coords.add_coord([])
            self.coords.add_coord({'a', 1})


class TestCoordArraySelection(unittest.TestCase):

    def setUp(self):
        self.array = np.array([[0,0,0], [1,1,1], [2,2,2]])
        self.coords = mastsel.CoordArray(self.array)
        self.sel_idxs = [0,2]
        self.coord_selection = mastsel.CoordArraySelection(self.coords, self.sel_idxs)

    def tearDown(self):
        pass

    def test_constructor(self):
        with self.assertRaises(AssertionError):
            mastsel.CoordArraySelection(self.coords, 'a')
            mastsel.CoordArraySelection(self.coords, ['a','b'])
            mastsel.CoordArraySelection(self.coords, -1)
            mastsel.CoordArraySelection(self.coords, [-1,-2])
            mastsel.CoordArraySelection(self.coords, [-1,2])
            mastsel.CoordArraySelection(self.coords, [])
            mastsel.CoordArraySelection({}, [1])

    def test_getitem(self):
        for i, idx in enumerate(self.sel_idxs):
            npt.assert_equal(self.coord_selection.container[idx], self.array[idx])
            npt.assert_equal(self.coord_selection.data[i], self.array[idx])
            npt.assert_equal(self.coord_selection[i], self.array[idx])

    def test_coords(self):
        target_coords = np.array([[0,0,0], [2,2,2]])
        npt.assert_equal(target_coords, self.coord_selection.coords)

class TestPoint(unittest.TestCase):

    def setUp(self):
        self.point1_coord = np.array([0,1,0])
        self.point1 = mastsel.Point(self.point1_coord)

        self.array = np.array([[0,0,0], [1,1,1], [2,2,2]])
        self.coord_array = mastsel.CoordArray(self.array)

        self.point2_coord = self.coord_array[0]
        self.point2 = mastsel.Point(self.point2_coord)

        self.point3_coord = np.array([0,1,0])
        self.point3 = mastsel.Point(self.point3_coord)

        self.bad_point_2d = mastsel.CoordArray(np.array([0,1]))
        self.bad_point_4d = mastsel.CoordArray(np.array([0,1,2,3]))
    def tearDown(self):
        pass

    def test_constructor(self):
        with self.assertRaises(AssertionError):

            # wrong dimension points
            mastsel.Point(self.bad_point_2d)
            mastsel.Point(self.bad_point_4d)

            # from existing CoordArray
            mastsel.Point(coord_array=np.array([1,2,3]))
            mastsel.Point(coord_array=self.coord_array, array_idx=3)
            mastsel.Point(coord_array=self.coord_array, array_idx='b')

    def test_overlaps(self):
        self.assertFalse(self.point1.overlaps(self.point2))
        self.assertTrue(self.point1.overlaps(self.point3))
        with self.assertRaises(AssertionError):
            self.point1.overlaps(np.array([0,1,0]))
            self.point1.overlaps([0,1,0])

class TestSelectionType(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

class TestSelectionTypeLibrary(unittest.TestCase):

    def setUp(self):
        pass

    def tearDown(self):
        pass

if __name__ == "__main__":
    from mast import selection
    # doctests
    print("\n\n\n Doc Tests\n-----------")
    nfail, ntests = doctest.testmod(selection, verbose=True)

    # unit tests
    print("\n\n\n Unit Tests\n-----------")
    unittest.main()
