"""Basic classes for building classes that are selection
capable and selection aware.

"""

import numpy as np
import collections as col
import collections.abc as colabc
from copy import copy

__all__ = ['SelectionMember', 'GenericSelection', 'IndexedSelection',
           'SelectionsDict', 'SelectionsList', 'CoordArray',
           'CoordArraySelection', 'Point', 'Selection']

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

    def __init__(self, member, flags=None):

        if flags is not None:
            assert (issubclass(type(flags), colabc.Set) or \
                issubclass(type(flags), colabc.Sequence)) and \
                not isinstance(flags, str), \
                "flags must be a container and not a string"
            assert all([isinstance(flag, str) for flag in list(flags)]), \
                "all flags must be strings, given{}".format(flags)

        super().__init__()
        self.member = member

        # list of selections
        self._registry = []
        # list of flags for specific kinds of selections
        self._flags = set()
        if flags:
            self._flags.update(flags)


    def __repr__(self):
        return str(self.__class__)

    # def get_selections(self, level=None):
    #     """Returns all the selection objects referred to in the registry,
    #     recursively to the level specified, default is all.

    #     Examples
    #     --------

    #     Get all the selections this SelectionMember is a part of:

    #     >>> container = [SelectionMember('a'), SelectionMember('b'), SelectionMember('c')]
    #     >>> selection = Selection(container, [0])
    #     >>> container[0].get_selections()[0] is selection
    #     True

    #     If these selections then become part of other selections we
    #     want to get all of them recursively:

    #     >>> len(container[0].get_selections())
    #     1
    #     >>> meta_selection = Selection([selection], [0])
    #     >>> len(container[0].get_selections())
    #     2

    #     For a specific depth of selections:

    #     >>> len(container[0].get_selections(level=0))
    #     1

    #     >>> len(container[0].get_selections(level=1))
    #     2

    #     """

    #     # get the selections on this level
    #     selections = [tup[-1] for tup in self._registry]
    #     # return them if the level is at 0
    #     if level == 0:
    #         return selections
    #     # otherwise get new selections from lower levels
    #     new_selections = []
    #     # from each selection
    #     for selection in selections:
    #         # if the level is not infinite reduce the level depth and continue
    #         if level is not None:
    #             new_selections.extend(selection.get_selections(level=(level - 1)))
    #         else:
    #             new_selections.extend(selection.get_selections(level=level))

    #     selections.extend(list(new_selections))

    #     return selections

    def get_selections(self, selection_types=None, flags=None, level=None):
        """Recursively searches for and finds all selections of this
        SelectionMember of the specified types and flags to the given level.

        If any option is None (default) then that criteria is ignored.

        Examples
        --------

        >>> container = [SelectionMember('a'), SelectionMember('b'), SelectionMember('c')]
        >>> selection = Selection(container, [0])
        >>> container[0].get_selections()[0] is selection
        True

        If these selections then become part of other selections we
        want to get all of them recursively:

        >>> len(container[0].get_selections())
        1
        >>> meta_selection = Selection([selection], [0])
        >>> len(container[0].get_selections())
        2

        For a specific depth of selections:

        >>> len(container[0].get_selections(level=0))
        1

        >>> len(container[0].get_selections(level=1))
        2

        When we add more selections on the SelectionMembers we can
        select specific ones:

        >>> indexed_selection = IndexedSelection(container, [1])
        >>> selections_list = SelectionsList([selection, indexed_selection], flags=['meta-selection'])
        >>> container[0].get_selections(selection_types=[IndexedSelection], level=0)[0] is indexed_selection
        True
        >>> container[0].get_selections(flags=['meta-selection'])[0] is selections_list
        True
        """

        if flags:
            flags = set(flags)
        # get the selections on this level if they match the criteria
        # if both selection_types and flags are given
        if selection_types is not None and flags is not None:
            # add them if
            selections = [tup[-1] for tup in self._registry if
                          # the selection's type is in the selection types
                          type(tup[-1]) in selection_types and
                          # and at least one of the flags is the same
                          not flags.isdisjoint(tup[-1].flags)]

        # only flags are given
        elif selection_types is None and flags is not None:
            selections = [tup[-1] for tup in self._registry if
                          not flags.isdisjoint(tup[-1].flags)]

        # only the selection_types are given
        elif selection_types is not None and flags is None:
            selections = [tup[-1] for tup in self._registry if
                          type(tup[-1]) in selection_types]

        # neither are given add everything
        elif selection_types is None and flags is None:
            selections = [tup[-1] for tup in self._registry]

        # return them if the level is at 0
        if level == 0:
            return selections
        # otherwise get new selections from lower levels
        new_selections = []
        # from each selection
        for selection in selections:
            # if the level is not infinite reduce the level depth and continue
            if level is not None:
                new_selections.extend(
                    selection.get_selections(selection_types=selection_types,
                                              flags=flags,
                                              level=(level - 1)))
            else:
                new_selections.extend(
                    selection.get_selections(selection_types=selection_types,
                                              flags=flags,
                                              level=level))

        selections.extend(list(new_selections))

        return selections

    @property
    def registry(self):
        """The selections this SelectionMember is a part of.
        Each entry is a tuple (key, selection), where the key is the
        way of accessing this SelectionMember from the selection with
        __getitem__.

        Examples
        --------

        If the SelectionMember is not part of a selection it will be empty:
        >>> a = SelectionMember('a')
        >>> a.registry
        []

        When a properly implemented selection is made:
        >>> idxsel = IndexedSelection([a], sel=[0])
        >>> a.registry
        [(0, <class 'mast.selection.IndexedSelection'>)]

        >>> a.registry[0][1][a.registry[0][0]] is a
        True

        Be careful with making temporary selections because the
        reference is still referenced in the registry and can cause a
        "memory leak":
        >>> len(a.registry)
        1
        >>> IndexedSelection([a], sel=[0])
        <class 'mast.selection.IndexedSelection'>
        >>> len(a.registry)
        2

        """
        return self._registry

    @property
    def flags(self):
        """A set of flags that certain selections may raise when they register
        themselves in this SelectionMember.
        """
        return self._flags

    def register_selection(self, key, selection, flags=None):
        """Adds a selection and this SelectionMembers key to __getitem__ from
        that selection to the registry.

        """
        self._registry.append((key, selection))
        if flags:
            self._flags.update(flags)

class GenericSelection(SelectionMember):
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
    def __init__(self, container, flags=None):

        super().__init__(self, flags=flags)
        assert '__getitem__' in dir(container), \
            "container must implement `__getitem__`, {} does not".format(
                type(container))
        assert container, "container must have at least one SelectionMember element"

        self.container = container
        # self.sel_ids = SelectionIDs

    def __repr__(self):
        return str(self.__class__)
        # return "{0}[{1}]".format(self.container, self.sel_ids)

class Selection(GenericSelection, col.UserList):
    """ A simple selection of a container.

    Examples
    --------

    >>> container = [SelectionMember(str) for str in ['a','b','c']]
    >>> selection = Selection(container, sel=[0,2])
    >>> selection[0].member
    'a'
    >>> selection[1].member
    'c'

    """

    def __init__(self, container, sel, flags=None):
        super().__init__(container, flags=flags)

        # handle sel inputs
        # TODO add support for slices
        assert type(sel) in [int, list], \
            "sel must be either a positive int, list of positive ints or a slice"

        if isinstance(sel, int):
            assert sel >=0, "an integer sel must be non-negative, {}".format(sel)
            sel = [sel]
        elif isinstance(sel, list):
            assert sel, "a list sel must be nonempty, {}".format(sel)
            assert all([(lambda x: False if x < 0 else True)(b) for b in sel]), \
                "all values in selection keys must be non-negative"

        assert all([issubclass(type(member), SelectionMember) for member in container]), \
            "container members must be a subclass of SelectionMember"

        # make the selections from container
        self.sel_ids = sel
        idx = 0
        for sel_idx in sel:
            member = container[sel_idx]
            self.append(member)
            # set this selection in the SelectionMember registry
            member.register_selection(idx, self, flags=flags)
            idx += 1

class KeySelection(GenericSelection, col.UserDict):
    pass

class IndexedSelection(GenericSelection, col.UserDict):
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
    def __init__(self, container, sel, flags=None):
        super().__init__(container, flags=flags)
        assert all([issubclass(type(member), SelectionMember) for member in container]), \
            "container members must be a subclass of SelectionMember"

        # make the selections from container
        self.sel_ids = sel
        for sel_idx in sel:
            self[sel_idx] = container[sel_idx]
            # set this selection in the SelectionMember registry
            self[sel_idx].register_selection(sel_idx, self, flags=flags)

    def __repr__(self):
        return str(self.__class__)


class CoordArray(SelectionMember):
    """A numpy array that is SelectionMember.

    Examples
    --------

    Just wraps some array:
    >>> import numpy as np

    e.g. the 3D coordinates of 3 points
    >>> arr = np.array([[0,0,0], [1,1,1], [2,2,2]])
    >>> coords = CoordArray(arr)
    >>> coords
    <class 'mast.selection.CoordArray'>
    >>> coords.coords
    array([[0, 0, 0],
           [1, 1, 1],
           [2, 2, 2]])

    Except we can add coords like records in a table, and the index of
    the new record is returned:

    >>> coords.add_coord(np.array([3,3,3]))
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
    array([0, 0, 0])
    >>> coordsel[1]
    array([2, 2, 2])

    """

    def __init__(self, array, sel, flags=None):

        assert issubclass(type(array), CoordArray), \
            "array must be a subclass of CoordArray, not {}".format(
                type(array))

        super().__init__(array, flags=flags)

        # handle sel inputs
        # TODO add support for slices
        assert type(sel) in [int, list], \
            "sel must be either a positive int, list of positive ints or a slice"

        if isinstance(sel, int):
            assert sel >=0, "an integer sel must be non-negative, {}".format(sel)
            sel = [sel]
        elif isinstance(sel, list):
            assert sel, "a list sel must be nonempty, {}".format(sel)
            assert all([(lambda x: False if x < 0 else True)(b) for b in sel]), \
                "all values in selection keys must be non-negative"


        # make selections as views from the array
        self.sel_ids = sel
        self.data = []
        for sel_idx in sel:
            self.data.append(self.container[sel_idx:(sel_idx+1)][0])
            # set this selection in the CoordArray registry
            self.container.register_selection(sel_idx, self, flags=flags)

    def __getitem__(self, idx):
        return self.data[idx]

    def __repr__(self):
        return str(self.__class__)

    @property
    def _coords(self):
        from functools import reduce
        return reduce(lambda x,y: np.concatenate((x,y)), [[value] for value in self.data])

    @property
    def coords(self):
        """ Array of only the selected coordinates. """
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

class SelectionsDict(SelectionMember, col.UserDict):
    """ A dictionary of collections of SelectionMembers.
e.g. {'strings' : [StringSelection, StringSelection] 'ints' :
[IntSelection, IntSelection]}

    """
    def __init__(self, selection_dict=None, flags=None):
        if not selection_dict:
            self.data = {}

        super().__init__(selection_dict, flags=flags)

        # add the selection_dict to the data
        if selection_dict:
            assert issubclass(type(selection_dict), col.Mapping), \
                "selection_dict must be a subclass of collections.Mapping, not {}".format(
                    type(selection_dict))
            self.data = selection_dict

        # if values in the selection_dict are SelectionMembers update
        # their registries and flags
        for key, value in self.data.items():
            # if there are multiple
            try:
                for member in value:

                    if issubclass(type(member), SelectionMember):
                        member.register_selection(key, self, flags=flags)
            # if there is only one
            except TypeError:
                if issubclass(type(value), SelectionMember):
                    value.register_selection(key, self, flags=flags)

    def __repr__(self):
        return str(self.__class__)


class SelectionsList(SelectionMember, col.UserList):
    def __init__(self, selection_list=None, flags=None):
        if not selection_list:
            self.data = []

        super().__init__(selection_list, flags=flags)

        if selection_list:
            assert issubclass(type(selection_list), col.Sequence), \
                "selection_dict must be a subclass of collections.Sequence, not {}".format(
                    type(selection_list))
            self.data = selection_list

        # if values in the selection_list are SelectionMembers update
        # their registries
        for idx, member in enumerate(self.data):
            if issubclass(type(member), SelectionMember):
                member.register_selection(idx, self, flags=flags)

    def __repr__(self):
        return str(self.__class__)

if __name__ == "__main__":
    pass
