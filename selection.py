"""Abstract basic classes for building classes that are selection
capable and selection aware.

"""

import numpy as np
import collections as col
import collections.abc as colabc
from copy import copy

__all__ = ['SelectionMember', 'GenericSelection', 'IndexedSelection',
           'SelectionDict', 'SelectionList', 'CoordArray',
           'CoordArraySelection', 'SelectionType', 'SelectionTypeLibrary',
           'Point',]

class SelectionMember(object):
    """The base class that allows an object to be part of a
    selection. Implements a registry which is used to keep track of
    which selections it is a part of. Selection type classes must
    implement a mechanism to add themselves to this registry.

    Examples
    --------

    SelectionMember simply wraps an object:

    >>> SelectionMember('a')
    <class 'mast.selection.SelectionMember'>

    Which you can access:

    >>> a = SelectionMember('a')
    >>> a.member
    'a'

    And can access the selections which have registered themselves as
    containing this SelectionMember:

    >>> a.registry
    []

    """

    def __init__(self, member):
        super().__init__()
        self.member = member

        # list of selections
        self._registry = []

        # TODO move this to a SystemMember class for mixing in later
        self._in_system = False
        # TODO turn this into a general method where the selections in
        # the selection members registry are recursively searched for
        # one that is true for _in_system. Currently implemented
        # induvidually for each class I need for it.
        self._system = None

    def __repr__(self):
        return str(self.__class__)

    def get_selections(self):
        return self._registry

    @property
    def registry(self):
        return self._registry

    def register_selection(self, key, selection):
        self._registry.append((key, selection))

class GenericSelection(SelectionMember, col.UserDict):
    """The most basic class for making selections of SelectionMember
    objects. Requires only the container to be selected from, but not
    which are selected, see subclasses for this implementation.

    Examples
    --------

    The container of SelectionMembers from which to select from:
    >>> container = [SelectionMember('a'), SelectionMember('b')]
    >>> gensel = GenericSelection(container)
    >>> gensel
    <class 'mast.selection.GenericSelection'>

    The container can be accessed:
    >>> gensel.container
    [<class 'mast.selection.SelectionMember'>, <class 'mast.selection.SelectionMember'>]

    See IndexedSelection, and CoordArraySelection for classes with the
    selection mechanism implemented.


    """
    def __init__(self, container):
        super().__init__(self)
        assert '__getitem__' in dir(container), \
            "container must implement `__getitem__`, {} does not".format(
                container)
        assert container, "container must have at least one SelectionMember element"

        self.container = container
        # self.sel_ids = SelectionIDs

    def __repr__(self):
        return str(self.__class__)
        # return "{0}[{1}]".format(self.container, self.sel_ids)

class IndexedSelection(GenericSelection):
    """ A selection of a container by indices.

    Example
    -------

    Make non-contiguous selections from a container:

    >>> container = [SelectionMember(str) for str in ['a','b','c']]
    >>> idxsel = IndexedSelection(container, sel=[0,2])
    >>> idxsel
    <class 'mast.selection.IndexedSelection'>

    Access via dictionary syntax:

    >>> idxsel[0]
    <class 'mast.selection.SelectionMember'>
    >>> idxsel[0].member
    'a'

    The registry for the SelectionMember knows it is selected:

    >>> idxsel[0].registry
    [(0, <class 'mast.selection.IndexedSelection'>)]

"""
    def __init__(self, container, sel):
        super().__init__(container)
        assert issubclass(type(container[0]), SelectionMember), \
            "container members must be a subclass of SelectionMember, not {}".format(
                type(container[0]))

        # make the selections from container
        self.sel_ids = sel
        for sel_idx in sel:
            self[sel_idx] = container[sel_idx]
            # set this selection in the SelectionMember registry
            self[sel_idx].register_selection(sel_idx, self)

    def __repr__(self):
        return str(self.__class__)

    def register_selection_members(self, key, selection):
        """ Register this object in the child selections of another selection."""

        # TODO currently doesn't store the key for this selection,
        # only the toplevel one
        for this_key, selmemb in self.data.items():
            selmemb.register_selection(key, selection)

    def register_selection(self, key, selection):
        """ Register this IndexedSelection in a selection's registry."""
        self._registry.append((key, selection))
        # let the children selection members know they are now a part
        # of a higher-order selection
        self.register_selection_members(key, selection)

class SelectionDict(SelectionMember, col.UserDict):
    """ A dictionary of collections of SelectionMembers.
e.g. {'strings' : [StringSelection, StringSelection] 'ints' :
[IntSelection, IntSelection]}

    """
    def __init__(self, selection_dict=None):
        if not selection_dict:
            self.data = {}

        super().__init__(selection_dict)

        # add the selection_dict to the data
        if selection_dict:
            assert issubclass(type(selection_dict), col.Mapping), \
                "selection_dict must be a subclass of collections.Mapping, not {}".format(
                    type(selection_dict))
            self.data = selection_dict

        # if values in the selection_dict are SelectionMembers update
        # their registries
        for key, value in self.data.items():
            # if there are multiple
            try:
                for member in value:

                    if issubclass(type(member), SelectionMember):
                        member.register_selection(key, self)
            # if there is only one
            except TypeError:
                if issubclass(type(value), SelectionMember):
                    value.register_selection(key, self)

    def __repr__(self):
        return str(self.__class__)


class SelectionList(SelectionMember, col.UserList):
    def __init__(self, selection_list=None):
        if not selection_list:
            self.data = []

        super().__init__(selection_list)

        if selection_list:
            assert issubclass(type(selection_list), col.Sequence), \
                "selection_dict must be a subclass of collections.Sequence, not {}".format(
                    type(selection_list))
            self.data = selection_list

        # if values in the selection_list are SelectionMembers update
        # their registries
        for idx, member in enumerate(self.data):
            if issubclass(type(member), SelectionMember):
                member.register_selection(idx, self)

    def __repr__(self):
        return str(self.__class__)

class CoordArray(SelectionMember):
    """A numpy array that is SelectionMember.

    Examples
    --------

    Just wraps some array:
    >>> from numpy import array

    e.g. the 3D coordinates of 3 points
    >>> arr = array([[0,0,0], [1,1,1], [2,2,2]])
    >>> coords = CoordArray(arr)
    >>> coords
    <class 'mast.selection.CoordArray'>
    >>> coords.coords
    array([[0, 0, 0],
           [1, 1, 1],
           [2, 2, 2]])

    Except we can add coords like records in a table, and the index of
    the new record is returned:

    >>> coords.add_coord(array([3,3,3]))
    3

    """
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
        """Adds 1-D coordinate array and returns the index of the new
        coordinate in the array.

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

class CoordArraySelection(GenericSelection):
    """ A selection of coordinates records from a CoordArray.

    Examples
    --------

    >>> import numpy as np
    >>> arr = np.array([[0,0,0], [1,1,1], [2,2,2]])

    >>> coords = CoordArray(arr)

    >>> coordsel = CoordArraySelection(coords, [0,2])
    >>> coordsel
    <class 'mast.selection.CoordArraySelection'>
    >>> coordsel[0]
    array([[0, 0, 0]])
    >>> coordsel[2]
    array([[2, 2, 2]])

    """
    def __init__(self, array, sel):
        super().__init__(array)

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
            self.container.register_selection(sel_idx, self)

    def __repr__(self):
        return str(self.__class__)

    @property
    def _coords(self):
        from functools import reduce
        return reduce(lambda x,y: np.concatenate((x,y)), [value for value in self.values()])

    @property
    def coords(self):
        return self._coords

class Point(CoordArraySelection):
    """An n-dimensional point possibly drawn from an existing CoordArray.

    Examples
    --------

    Without a pre-existing CoordArray the point will make it's own:
    >>> import numpy as np

    >>> point_coord = np.array([0,1,0])

    >>> point1 = Point(point_coord)
    >>> point1
    <class 'mast.selection.Point'>
    >>> point1.container.coords
    array([[0, 1, 0]])

    However, many points made this way will be fragmented and will not
    benefit from fast numpy operations on arrays. So if you have a lot of
    related points, make a CoordArray first and then make selections
    from that:

    >>> import numpy as np
    >>> arr = np.array([[0,0,0], [1,1,1], [2,2,2]])
    >>> coords = CoordArray(arr)

    To make a new point and add it to an existing CoordArray:

    >>> point2 = Point(point_coord, coord_array=coords)
    >>> point2.container.coords
    array([[0, 0, 0],
           [1, 1, 1],
           [2, 2, 2],
           [0, 1, 0]])

    >>> coords.registry
    [(3, <class 'mast.selection.Point'>)]

    """

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

    @property
    def coords(self):
        return self._coords[0]

    def overlaps(self, other):
        """ Test if this point overlaps with another Point.

        Examples
        --------

        >>> import numpy as np
        >>> point1 = Point(np.array([0.0, 0.0, 0.0]))
        >>> point2 = Point(np.array([0.0, 0.0, 1.0]))
        >>> point1.overlaps(point2)
        False

        """
        assert issubclass(type(other), Point), \
            "Other must be a subclass of Point, not {}".format(type(other))
        return np.all(np.isclose(self.coords, other.coords))


class SelectionType(object):
    """Base type for other Types."""

    def __init__(self, attr_dict=None):
        if attr_dict:
            assert isinstance(attr_dict, dict), \
                "The attr_dict must be a dictionary, not {}".format(
                    type(attr_dict))
            # if there is no 'name' in the attr_dict set it to the empty string
            if 'name' not in attr_dict.keys():
                attr_dict['name'] = ''
            self.__dict__.update(attr_dict)

    def __eq__(self, other):
        return self.__dict__ == other.__dict__
    def __ne__(self, other):
        return self.__dict__ != other.__dict__
    def __lt__(self, other):
        return set(self.__dict__.keys()) < set(other.__dict__.keys())
    def __gt__(self, other):
        return set(self.__dict__.keys()) > set(other.__dict__.keys())
    def __le__(self, other):
        if self == other or self < other:
            return True
        else:
            return False
    def __ge__(self, other):
        if self == other or self > other:
            return True
        else:
            return False

class SelectionTypeLibrary(col.UserDict):
    """Class that keeps track of a collection of types with methods for
querying matches to attributes of types in the library, to reduce
duplication, promote standardization of attributes, and keep clear
naming distinctions.

    """

    def __init__(self):
        super().__init__()
        self._names_counter = {}

    def add_type(self, a_type, type_name, rename=False):
        """adds a SelectionType to the SelectionTypeLibrary using the
type_name.  If you try to add a duplicate with the same name it is
silently ignored. If you try to add a type with the same name as one
in the library but with different attributes it will raise an
error.

        """

        assert issubclass(type(a_type), SelectionType), \
            "added types must be a subclass of SelectionType, not {}".format(
                type(a_type))
        assert isinstance(type_name, str) or isinstance(type_name, int), \
            "type_name must be type str or int, not {}".format(type(type_name))
        if rename:
            assert isinstance(rename, bool), "rename should be type bool, not {}".format(
                type(rename))
            if type_name in ['',None]:
                type_name = 'NoName'
        if type_name not in self.data.keys():
            self.data[type_name] = a_type
            self._names_counter[type_name] = 1
        elif self.data[type_name] == a_type:
            pass
        elif rename is False:
            raise ValueError(
                "{0} is already in the {2}, {1}, "
                "cannot redefine attributes under the same type_name".format(
                    type_name, id(self), type(self)))
        else:
            self.data[type_name+str(self._names_counter[type_name])] = a_type
            self._names_counter[type_name] += 1

    def attributes_match(self, a_type):
        """Check if the attributes of an AtomType are equivalent to any
AtomType already in the library.

        """

        from itertools import product
        for pair in product(self.data.values(), [a_type]):
                if pair[1] == pair[0]:
                    return True
        return False

if __name__ == "__main__":

    from mast.interactions import Association

    from mast.selection import *

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
    class StrSelection(IndexedSelection):
        def __init__(self, strings, sel):
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

    print("testing point overlaps")
    print(point1.overlaps(point2))

    print("Making SelectionDict")
    seldict = SelectionDict()
    print(seldict)
    print(seldict.registry)
    seldict2 = SelectionDict({'points' : [point1, point2],
                              'strings' : strings,
                              'coords' : [coordsel, CoordArraySelection(coords, [1,2])]})
    print(seldict2)
    print(seldict2.registry)
    print("retrieving a seldict from it's members registries")
    selected_member = seldict2['points'][0]
    registry_selector = selected_member.registry[0][1]
    print(registry_selector)
    print("is the original point in this list?")
    print(selected_member in registry_selector['points'])
    print("getting the other points in this selection")
    print("other_points in seldict2")
    other_points = [p for p in registry_selector if p is not selected_member]

    print("registries from seldict2 selections")
    print("point1: ", point1.registry)
    print("point2: ", point2.registry)
    print("strings[0]: ", strings[0].registry)
    print("coordsel: ", coordsel.registry)

    print("Making SelectionList")
    sellist = SelectionList()
    print(sellist)
    print(sellist.registry)
    sellist2 = SelectionList([point1, point2])
    print(sellist2)
    print(sellist2.registry)

    print("Making Association")
    assoc = Association(association_list=[point1, point2], association_type=None)
    print(assoc)

    print("testing Type base classes")
    print("making SelectionTypes")
    seltype = SelectionType({'name':'a'})
    seltype2 = SelectionType({'name':'b'})

    print("making SelectionTypeLibrarys")
    seltype_lib = SelectionTypeLibrary()
    seltype_lib.add_type(seltype, type_name=seltype.name)
    seltype_lib.add_type(seltype2, type_name=seltype2.name)
    print("the library with two types")
    print(seltype_lib)
    seltype_dup = SelectionType({'name': 'a'})
    print("A duplicate type, attributes_match result:")
    print(seltype_lib.attributes_match(seltype_dup))
    seltype_lib.add_type(seltype_dup, type_name=seltype_dup.name)
    print("Adding it with the same name makes no new entry")
    print(seltype_lib)
    seltype_lib.add_type(seltype_dup, type_name='a_dup')
    print("Using a different name will add it though")
    print(seltype_lib)
    print("If you add a type with the same name but different attributes,"
          " it default raises an error. Set rename flag to true if you want to"
          " allow it to rename the type. Default is just numbering")
    new_seltype = SelectionType({'name':'a', 'newattr':0})
    seltype_lib.add_type(new_seltype, type_name='a', rename=True)
    print(seltype_lib)
