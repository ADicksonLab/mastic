import copy

__all__ = ['TrackedMember', 'TrackedList', 'Selection']

# the magic of the TrackedList class, used to decorate functions which
# should propogate a signal of changes.
def _changes(func):
    """ Decorator to indicate the list has changed """
    def update_func(self, *args, **kwargs):
        self.changed = True
        self.indexed = False
        return func(self, *args, **kwargs)
    return update_func

def _type_check(func):
    """ Decorator to do a type check incoming members to the TrackedList. """

    def type_check_func(self, *args, **kwargs):
        return_value = func(self, *args, **kwargs)
    return type_check_func

# doesn't work because you have to rturn the values in the getitem etc
def _check(func):
    """Check members of list to see if they have the `changed` flag set to
`True`, if so set this `TrackedList` flag to `True` as well.

    """
    def check_func(func, *args, **kwargs):
        # allow resolution of function
        func(self, *args, **kwargs)

        # check to see if any members changed
        if self.check_members() == True:
            self.changed = True

def check_members(tlist):
    # iterate each member until a True one is found for efficiency
    List_iter = iter(self)
    try:
        member = next(List_iter)
        while member.changed == False:
            member = next(List_iter)
        # if StopIteration not reached set changed = True
        return True
    # if none were found return False
    except StopIteration:
        return False

class TrackedMember(object):
    """Template reference for a tracked and indexed object used in a
TrackedList. Indices (self.idx) are stored in each TrackedMember and
can only be set when using the TrackedMember constructor, when added
to a TrackedList when idx is None, and privately by TrackedList
methods (i.e. index_members).

The ids property is a dict {id_type : id} so that members can be
    members of many different types of Selection types.

    """

    def __init__(self, idx=None, ids=None):
        self._changed = False

        # idx
        if idx is None:
            self._idx = None
        elif isinstance(idx, int):
            if idx < 0:
                raise ValueError("idx can't be negative")
            else:
                self._idx = idx
        else:
            raise TypeError("idx must be a positive integer not type {}".format(type(idx)))

        # ids
        if ids is None:
            self._ids = None
        elif isinstance(ids, dict):
            self._ids = ids
        else:
            raise TypeError("ids in constructor must be type dict, not type {}".format(type(ids)))

    def __copy__(self):
        return TrackedMember(idx=self._idx, ids=self._ids)
    
    @property
    def changed(self):
        return self._changed

    @changed.setter
    def changed(self, value):
        if isinstance(value, bool):
            self._changed = value
        else:
            raise TypeError("Must be type bool not type {}".format(type(value)))

    @property
    def idx(self):
        return self._idx

    @property
    def ids(self):
        return self._ids

class TrackedList(list):
    """This creates a list with members of a single type that allows you
    to see if anything has changed in those members implementing the
    TrackedMember template.

    Attributes
    ----------
    changed : ``bool``
        Determines if something has been done to fundamentally change the
        underlying topology defined by this list such that the topology needs to
        be rebuilt
    indexed : ``bool``
        A flag to determine whether or not the items in a tracked list need to
        be indexed or not.
    type
        The type of the tracked members in the list.
    List : ``list[tracked_member_type]``
        The actual list of the TrackedMembers

    Examples
    --------
    >>> tl = TrackedList()
    >>> tl.append(Atom())
    >>> tl.append(Atom())
    >>> tl.append(Atom())
    >>> tl.indexed, tl.changed
    (True, True)
    >>> tl.index_members()
    >>> tl.indexed, tl.changed
    (False, True)
    >>> tl.changed = False # Must do when changes have been incorporated
    >>> tl.indexed, tl.changed
    (False, False)

    """
    def __init__(self, List=None):
        if List is None:
            self._type = None
            self._List = []
            self._changed = False
            self._indexed = True
        elif not List:
            self._type = None
            self._List = []
            self._changed = False
            self._indexed = True
        elif isinstance(List, list):
            # check to make sure the members are inherited from TrackedMembers
            if issubclass(type(List[0]), TrackedMember):
                self._type = type(List[0])
                self._List = List
                self._changed = True
                self._indexed = True
            else:
                raise TypeError(
                    "Elements in list constructor must be a subclass of TrackedMember, not type {}".format(type(List[0])))

        elif isinstance(List, TrackedList):
            self._type = type(List[0])
            self._List = List
            self._changed = True
        else:
            raise TypeError(
              "Must provide List, TrackedList, or None not {} type".format(type(List)))

    def __str__(self):
        return self._List.__str__()

    def __repr__(self):
        str = "type={3}\nindexed={0}\nchanged={1}\n{2}".format(\
                                                               self.indexed, self.changed, self._List.__repr__(), self.type)
        return str


    # Immutable container protocol methods

    def __len__(self):
        return len(self._List)

    # think about whether or not a shallow copy is the right thing
    # here
    def __getitem__(self, index):
        """ Returns a shallow copy of the List"""
        print("TrackedList __getitem__ no copy", self._List[index])
        print("TrackedList __getitem__ copy", copy.copy(self._List[index]))
        return copy.copy(self._List[index])

    # Mutable container protocol methods including slicing
    @_changes
    def __setitem__(self, index, value):
        if self._type is None:
            self._type = type(value)

        if isinstance(value, self._type):
            self._List[index] = value
        else:
            raise TypeError("Members must be of type {0} not type {1}".format(self.type, type(value)))

    @_changes
    def __delitem__(self, index):
        del self._List[index]

    # iterable protocol
    def __iter__(self):
        for i in self._List:
            yield i

    # Descriptor Protocol for the List attribute

    # useful for if someone manually sets the List attribute instead
    # of just modifying it e.g.:
    # l = TrackedList()
    # l.append('a')
    # p = [1,2,3]
    # l.List = p
    #
    # in which case you would want to update the flags
    @property
    def List(self):
        """The List attribute of the TrackedList"""
        return self._List


    @List.setter
    @_changes
    def List(self, value):
        self._List = value
        self.type = type(value[0])

    @List.deleter
    @_changes
    def List(self):
        del self._List

    # Descriptor for the changed property
    @property
    def changed(self):
        return self._changed

    @changed.setter
    def changed(self, value):
        if isinstance(value, bool):
            self._changed = value
        else:
            raise TypeError( \
              "changed attribute must be type bool not {} type".format(type(value)))

    # Descriptor for the indexed property
    @property
    def indexed(self):
        return self._indexed

    @indexed.setter
    def indexed(self, value):
        if value == True:
            raise ValueError("Cannot manually set indexed to True, use the index_members function")
        elif value == False:
            self._indexed = value
        else:
            raise TypeError("Only valid value is bool False, not type {}".format(type(value)))

    # Descriptor for the type property
    @property
    def type(self):
        return self._type

    @type.setter
    def type(self, new_type):
        if self._type is None:
            self._type = new_type
        else:
            raise AttributeError("The type is already set to {}".format(self._type))

    # List emulation functions
    @_changes
    def pop(self, idx=-1):
        return self._List.pop(idx)

    @_changes
    def remove(self, thing):
        self._List.remove(thing)

    @_changes
    def append(self, thing):
        if type(thing) == self.type:
            self._List.append(thing)
        else:
            raise TypeError("Must append members of the same type not {}".format(type(thing)))

    @_changes
    def extend(self, tlist):
        if tlist and isinstance(tlist, TrackedList):
            if type(tlist[0]) == self.type:
                if tlist[0].idx is None:
                    self._List.extend(tlist._List)
                    self.index_members()
                else:
                    self._List.extend(tlist._List)
            else:
                raise TypeError("Members of extension list must be type {0}, not {1}".format(self._type, type(tlist)))
        elif isinstance(tlist, list):
            if type(tlist[0]) == self.type:
                self._List.extend(tlist)
        else:
            raise TypeError(
                "Must extend with another TrackedList or list, not {} type".format(type(tlist)))


    @_changes
    def insert(self, index, thing):
        self._List.insert(index, thing)

    @_changes
    def clear(self):
        self._List = []

    @_changes
    def reverse(self):
        self._List.reverse()

    @_changes
    def sort(self, key=None, reverse=False):
        self._List.sort(key=key, reverse=reverse)

    def copy(self):
        return self._List[:]

    # will have to rethink this one later
    def index(self, value):
        return self._List.index(value)


    # Methods that return another instance
    def __add__(self, other):
        """ Return new TrackedList instance of the concatenation of two TrackedLists."""
        tlist = TrackedList()
        tlist.extend(self)
        tlist.extend(other)
        return tlist

    def index_members(self):
        """
        Assigns the idx variable for every member of this list to its place in
        the list, if the members of this list implement the idx attribute.
        """
        for i, item in enumerate(self):
            # try to set the private idx property if available.
            # public idx has no setter.
            try:
                item._idx = i
            # if the object doesn't implement the idx protocol no
            #action necessary and TrackedList will not fail
            except AttributeError:
                pass
        self._indexed = True

    def member_idx(self):
        """Return the indices of the members in List.

        """
        return [member.idx for member in self._List]

    def check_members(self):
        """ Attribute convenience for check_members function. """
        return check_members(self)


class Selection(object):
    """ Base class for containing a TrackedList and indices for that list. 

    container : the container implementing a TrackedList like interface
    sel : The selection indices for the container
    container_type : the type of the container
    member_type : the type of the things in the container
"""

    def __init__(self, container=None, sel=None):

        if container is None:
            self._container = None
            if sel:
                raise ValueError("Cannot set selection while container is None")

        elif issubclass(type(container), TrackedList):
            self._container = container
            self._container_type = type(container)
            self._member_type = container.type

            # set the selection
            if sel is None:
                self._sel = None
            elif isinstance(sel, slice):
                self._sel = sel
            elif isinstance(sel, list):
                it = iter(sel)
                index = next(it)
                try:
                    while True:
                        if isinstance(index, int):
                            if index < 0:
                                raise ValueError(
                                    "Selection indices cannot be negative, {}".format(index))
                            elif index >= len(self._container):
                                raise ValueError(
                                    "Selection out of bounds of the Selection's TrackedList")
                            else:
                                next(it)
                        else:
                            raise TypeError(
                                "Selection indices must be ints, not type {}".format(type(index)))
                except StopIteration:
                    self._sel = sel

            elif isinstance(sel, int):
                # TODO should negative indices work like a python
                # slice with negative numbers?
                if sel < 0:
                    raise ValueError(
                            "Selection indices cannot be negative, {}".format(index))
                elif sel >= len(self._container):
                    raise ValueError(
                        "Selection out of bounds of the Selection's TrackedList")
                else:
                    self._sel = sel
            else:
                raise TypeError(
                    "Selection indices must be type int or list of ints, not type {}".format(type(sel)))

        else:
            raise TypeError(
                "container type must be a subclass of TrackedList, not type {}".format(type(container)))

    def __len__(self):
        return len(self.sel)

    def __getitem__(self, index):
        selection = self.sel
        print("Selection __getitem__", selection[0])
        return self.sel[index]

    @property
    def container(self):
        return self._container

    @property
    def sel(self):
        """Returns a TrackedList shallow copy of the Selection indexed members
of the container TrackedList.

        """
        if isinstance(self._sel, int) or isinstance(self._sel, slice):
            print("in sel int/slice", self._container[0])
            return TrackedList([self._container[self._sel]])
        elif isinstance(self._sel, list):
            print("in sel list", self._container[0])
            return TrackedList([member for member in self._container if member.idx in self._sel])

    @property
    def sel_idx(self):
        return self._sel

    @property
    def container_type(self):
        return self._container_type

    @property
    def member_type(self):
        return self._member_type

class SelectionList(TrackedList):
    """ Collection class for multiple TrackedMembers contains """

    def __init__(self, List=None):
        if List is None:
            self._List = []
        super().__init__(List=None)
        self._type = Selection
