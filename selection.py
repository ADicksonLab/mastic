import numpy as np
import collections as col
import typing as ty
import collections.abc as colabc
from copy import copy
from functools import reduce, partial

__all__ = ['SelectionMember', 'GenericSelection', 'IndexedSelection',
           'SelectionDict', 'SelectionList', 'CoordArray',
           'CoordArraySelection', 'Point']

SELECTION_REGISTRY = {}
sel_reg_counter = 0

# TODO might want to make the keys the id(selection) however would
# then need to write stuff to tear down selections and remove them
# from the registry.
def register_selection(selection):
    """Register a selection in the global SELECTION_REGISTRY and return
it's key in the registry

    """
    global sel_reg_counter
    global SELECTION_REGISTRY
    sel_reg_counter += 1
    sel_reg_id = sel_reg_counter
    SELECTION_REGISTRY[sel_reg_id] = selection

    return sel_reg_id


class SelectionMember(object):
    def __init__(self, member):
        super().__init__()
        self.member = member
        # {selection_registry_id : id_in_selection}
        self.registry = {}

    def __repr__(self):
        return self.member.__repr__()

    def get_selections(self):
        global SELECTION_REGISTRY
        return [SELECTION_REGISTRY[sel_id] for sel_id in self.registry.keys()]

Selected = ty.TypeVar('Selected')
SelectionIDs = ty.TypeVar('SelectionIDs')
class GenericSelection(SelectionMember, col.UserDict, ty.Mapping[SelectionIDs, Selected]):
    def __init__(self, container: ty.Container):
        super().__init__(self)
        assert '__getitem__' in dir(container), \
            "container must implement `__getitem__`, {} does not".format(
                container)
        assert container, "container must have at least one SelectionMember element"

        self.container = container
        self.sel_ids = SelectionIDs

    def __repr__(self):
        return "{0}[{1}]".format(self.container, self.sel_ids)

class IndexedSelection(GenericSelection[int, Selected]):
    def __init__(self, container: ty.Sequence, sel: ty.Sequence[int]):
        super().__init__(container)
        assert issubclass(type(container[0]), SelectionMember), \
            "container members must be a subclass of SelectionMember, not {}".format(
                type(container[0]))

        # register this selection
        self.sel_reg_id = register_selection(self)

        # make the selections from container
        self.sel_ids = sel
        for sel_idx in sel:
            self[sel_idx] = container[sel_idx]
            # set this selection in the SelectionMember registry
            self[sel_idx].registry[self.sel_reg_id] = sel_idx

    def __repr__(self):
        return str(dict(self))

class SelectionDict(col.UserDict, SelectionMember):
    def __init__(self, selection_dict=None):
        super().__init__(self)
        if not selection_dict:
            self.data = {}
        else:
            assert issubclass(type(selection_dict), col.Mapping), \
                "selection_dict must be a subclass of collections.Mapping, not {}".format(
                    type(selection_dict))
            self.data = selection_dict

    def __repr__(self):
        return str(self.data)


class SelectionList(col.UserList, SelectionMember):
    def __init__(self, selection_list=None):
        super().__init__(self)
        if not selection_list:
            self.data = []
        else:
            assert issubclass(type(selection_list), col.Sequence), \
                "selection_dict must be a subclass of collections.Sequence, not {}".format(
                    type(selection_list))
            self.data = selection_list

    def __repr__(self):
        return str(self.data)

class CoordArray(SelectionMember):
    def __init__(self, array):
        assert isinstance(array, np.ndarray), \
            "array must be a numpy.ndarray, not {}".format(
                type(array))

        super().__init__(array)

    def __getitem__(self, idx):
        return self.member[idx]

    @property
    def coords(self):
        return self.member

    @coords.setter
    def coords(self, new_coords):
        assert isinstance(new_coords, np.ndarray), \
            "array must be a numpy.ndarray, not {}".format(
                type(new_coords))
        self.member = new_coords

    @property
    def shape(self):
        return self.coords.shape

    def add_coord(self, new_coord):
        """Given a 1-D coordinate array will add this coordinate to the array
and return the index of the new coordinate in the array.

        """
        assert isinstance(new_coord, np.ndarray), \
            "array must be a numpy.ndarray, not {}".format(
                type(new_coord))
        assert new_coord.shape[0] == self.coords.shape[-1], \
            "new_coord must be the same number of dimensions as the current coords ({0}), not {1}".format(
                self.coords.shape[-1], new_coord.shape)

        # add the coordinate to the coordinates
        self.coords = np.concatenate((self.coords, [new_coord]), axis=0)

        # return the index of the added coordinate
        return self.shape[0] - 1

class CoordArraySelection(GenericSelection[int, np.ndarray]):
    def __init__(self, array: CoordArray, sel: ty.Sequence[int]):
        super().__init__(array)

        # global register_selection
        self.sel_reg_id = register_selection(self)

        # TODO add support for slices
        assert type(sel) in [int, list], \
            "sel must be either a positive int, list of positive ints or a slice"

        # handle sel inputs
        if isinstance(sel, int):
            assert sel >=0, "an integer sel must be non-negative, {}".format(sel)
            sel = [sel]
        elif isinstance(sel, list):
            assert sel, "a list sel must be nonempty, {}".format(sel)
            assert all([(lambda x: False if x < 0 else True)(b) for b in sel]), \
                "all values in selection keys must be non-negative"

        # make selections as views from the array
        self.sel_ids = sel
        for sel_idx in sel:
            # slices like this are views into the array
            self[sel_idx] = self.container[sel_idx:(sel_idx+1)]
            # set this selection in the CoordArray registry
            self.container.registry[self.sel_reg_id] = sel_idx

    def __repr__(self):
        return str(dict(self))

    @property
    def coords(self):
        return reduce(lambda x,y: np.concatenate((x,y), axis=0), list(self.values()))

class Point(CoordArraySelection):
    def __init__(self, coords=None, coord_array=None, array_idx=None):

        # if not given a CoordArray just make our own using coords kwarg
        if not coord_array:
            assert isinstance(coords, np.ndarray), \
                "coords must be a numpy.ndarray, not type {}".format(type(coords))
            assert len(coords.shape) == 1, \
                "coords must be 1-dimensional, not {}".format(len(coords.shape))
            self._coord_array = CoordArray(np.array([coords]))

            # use the CoordArraySelection constructor
            super().__init__(self._coord_array, 0)

        # use coord_array
        else:
            assert issubclass(type(coord_array), CoordArray), \
                "coord_array must be type CoordArray, not {}".format(
                    type(coord_array))
            # set the private coord_array to None because it is not used
            self._coord_array = None

            # if an array index is given don't add a new entry
            if isinstance(array_idx, int):
                assert array_idx < coord_array.shape[0], \
                    "array_idx must be within the coord_array {0}, not {1}".format(
                        coord_array.shape[0], array_idx)
                point_idx = array_idx
            # add the coord record to the CoordArray
            else:
                point_idx = coord_array.add_coord(coords)

            # use the CoordArraySelection constructor
            super().__init__(coord_array, point_idx)

if __name__ == "__main__":
    # test GenericSelection
    gensel = GenericSelection([SelectionMember(None)])
    print(gensel)

    # test SelectionMember
    string_selmember = SelectionMember('a')
    print(string_selmember)

    # test IndexedSelection
    selmembers = [SelectionMember(i) for i in [0,1,2]]
    print("selmembers[0] is a part of these selections", selmembers[0].get_selections())
    idxsel = IndexedSelection(selmembers, [0,2])
    print("idxsel", idxsel)
    print("selmembers[0] is a part of these selections", selmembers[0].get_selections())
    idxsel2 = IndexedSelection(selmembers, [0,1])
    print("selmembers[0] is a part of these selections", selmembers[0].get_selections())

    # test a selection from a list of SelectionMembers as a new type
    strings = [SelectionMember('a'), SelectionMember('b'), SelectionMember('c')]
    class StrSelection(IndexedSelection[int, str]):
        def __init__(self, strings: ty.Sequence[str], sel: ty.Sequence[int]):
            super().__init__(strings, sel)

    stringsel = StrSelection(strings, [0])
    print(stringsel)

    array = np.array([[0,0,0], [1,1,1], [2,2,2]])
    print("making CoordArray")
    coords = CoordArray(array)
    print("making CoordArraySelection")
    coordsel = CoordArraySelection(coords, [0,1])
    print(coordsel)
    print(coordsel.coords)

    print("SELECTION_REGISTRY:", SELECTION_REGISTRY)

    print("making Point with it's own CoordArray")
    point_coord = np.array([0,1,0])
    point1 = Point(point_coord)
    print(point1)
    print(point1.coords)
    print("making Point in previous CoordArray")
    point_coord = np.array([0,1,0])
    point2 = Point(point_coord, coord_array=coords)
    print(point2)
    print(point2.coords)
    print("SELECTION_REGISTRY:", SELECTION_REGISTRY)
