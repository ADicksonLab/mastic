from copy import copy

__all__ = ['TrackedMember', 'TrackedList', 'TrackedDict', 'Selection', 'SelectionList', 'Association', 'AssociationType']

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
    members_iter = iter(self)
    try:
        member = next(members_iter)
        while member.changed == False:
            member = next(members_iter)
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

class TrackedList(TrackedMember):
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
    members : ``list[tracked_member_type]``
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
    def __init__(self, members=None, idx=None, ids=None):
        super().__init__(idx=idx, ids=ids)
        if members is None:
            self._type = None
            self._members = []
            self._changed = False
            self._indexed = True
            self._member_type = None
        elif not members:
            self._type = None
            self._members = []
            self._changed = False
            self._indexed = True
            self._member_type = None
        elif isinstance(members, list):
            # check to make sure the members are inherited from TrackedMembers
            if issubclass(type(members[0]), TrackedMember):
                self._type = type(members[0])
                self._members = members
                self._changed = True
                self._indexed = True
                self._member_type = type(members[0])
            else:
                raise TypeError(
                    "Elements in list constructor must be a subclass of TrackedMember, not type {}".format(type(members[0])))

        elif issubclass(type(members), TrackedList):
            self._type = type(members[0])
            self._members = members
            self._changed = True
            self._indexed = False
            self._member_type = members.member_type
        else:
            raise TypeError(
              "Must provide List, TrackedList, or None not {} type".format(type(members)))

    def __str__(self):
        return self._members.__str__()

    def __repr__(self):
        str = "member_type={3}\nindexed={0}\nchanged={1}\n{2}".format(\
                                                               self.indexed, self.changed, self._members.__repr__(), self.member_type)
        return str

    def __copy__(self):
        return type(self)(self)

    # Immutable container protocol methods

    def __len__(self):
        return len(self._members)

    def __getitem__(self, index):
        return self._members[index]

    # Mutable container protocol methods including slicing
    @_changes
    def __setitem__(self, index, value):
        if self._type is None:
            self._type = type(value)

        if isinstance(value, self._type):
            self._members[index] = value
        else:
            raise TypeError("Members must be of type {0} not type {1}".format(self.type, type(value)))

    @_changes
    def __delitem__(self, index):
        del self._members[index]

    # iterable protocol
    def __iter__(self):
        for i in self._members:
            yield i

    # Descriptor Protocol for the members attribute

    # useful for if someone manually sets the members attribute instead
    # of just modifying it e.g.:
    # l = TrackedList()
    # l.append('a')
    # p = [1,2,3]
    # l.members = p
    #
    # in which case you would want to update the flags
    @property
    def members(self):
        """The members attribute of the TrackedList"""
        return self._members


    @members.setter
    @_changes
    def members(self, value):
        self._members = value
        self.type = type(value[0])

    @members.deleter
    @_changes
    def members(self):
        del self._members

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
    def member_type(self):
        return self._member_type

    @member_type.setter
    def member_type(self, new_type):
        if self._member_type is None:
            self._member_type = new_type
        else:
            raise AttributeError("The member_type is already set to {}".format(self._member_type))

    # List emulation functions
    @_changes
    def pop(self, idx=-1):
        return self._members.pop(idx)

    @_changes
    def remove(self, thing):
        self._members.remove(thing)

    @_changes
    def append(self, thing):
        if type(thing) == self.type:
            self._members.append(thing)
        else:
            raise TypeError("Must append members of the same type not {}".format(type(thing)))

    @_changes
    def extend(self, other):
        other = copy(other)
        if len(other) == 0:
            pass
        elif isinstance(other, type(self)):
            if other.member_type is self.member_type:
                    self._members.extend(other._members)
                    self.index_members()
            else:
                raise TypeError(
                    "Members of extension list must be type {0}, not {1}".format(
                                                   self.member_type, other.member_type))
        else:
            raise TypeError(
                "Must extend with another {0}, not {1} type".format(type(self), type(other)))


    @_changes
    def insert(self, index, thing):
        self._members.insert(index, thing)

    @_changes
    def clear(self):
        self._members = []

    @_changes
    def reverse(self):
        self._members.reverse()

    @_changes
    def sort(self, key=None, reverse=False):
        self._members.sort(key=key, reverse=reverse)

    def copy(self):
        return self._members[:]

    # will have to rethink this one later
    def index(self, value):
        return self._members.index(value)


    # Methods that return another instance
    def __add__(self, other):
        """ Return new TrackedList instance of the concatenation of two TrackedLists."""
        tlist = copy(self.members)
        other = copy(other.members)
        new = type(self)(tlist + other)
        new.index_members()

        return new

    __radd__ = __add__

    def index_members(self):
        """
        Assigns the idx variable for every member of this list to its place in
        the list, if the members of this list implement the idx attribute.
        """
        for i, item in enumerate(self):
            item._idx = i
        self._indexed = True

    def member_idx(self):
        """Return the indices of the members in members.

        """
        return [member.idx for member in self._members]

    def check_members(self):
        """ Attribute convenience for check_members function. """
        return check_members(self)

class TrackedDict(TrackedMember):
    """ Class to contain a dict of TrackedMembers. """
    def __init__(self, tdict):
        # check to see if values in the dict are TrackedMembers
        for key, value in tdict.items():
            if not issubclass(type(value), TrackedMember):
                raise TypeError(
                    "Values in dictionary must be of type TrackedMember, not {0} for key {1}".format(
                        type(value), key))

        self._dict = tdict

    def __len__(self):
        return len(self._dict)

class Selection(TrackedMember):
    """ Base class for containing a TrackedList and indices for that list.

    container : the container implementing a TrackedList like interface
    sel : The selection indices for the container
    container_type : the type of the container
    member_type : the type of the things in the container
"""

    def __init__(self, container=None, sel=None, idx=None, ids=None):

        super().__init__(idx=idx, ids=ids)

        if container is None:
            self._container = None
            if sel:
                raise ValueError("Cannot set selection while container is None")

        elif issubclass(type(container), TrackedList):
            self._container = container
            self._container_type = type(container)
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

    def __copy__(self):
        return Selection(container=self._container, sel=self._sel)

    def __len__(self):
        return len(self.sel)

    def __getitem__(self, index):
        selection = self.sel
        return self.sel[index]

    @property
    def container(self):
        return self._container

    @property
    def sel(self):
        """Returns a TrackedList shallow copy of the Selection indexed members
of the container TrackedList.

        """
        if isinstance(self._sel, int):
            return TrackedList([self._container[self._sel]])
        elif isinstance(self._sel, slice):
            return TrackedList(self._container[self._sel])
        elif isinstance(self._sel, list):
            return TrackedList([member for member in self._container if member.idx in self._sel])

    @property
    def sel_idx(self):
        return self._sel

    @property
    def container_type(self):
        return self._container_type

    @property
    def member_type(self):
        return self.container.member_type

class SelectionList(TrackedList):
    """ Collection class for multiple Selection objects contains """

    def __init__(self, members=None, idx=None, ids=None):
        if not members:
            self._members = []
        elif isinstance(members, list) or issubclass(type(members), SelectionList):
            if not issubclass(type(members[0]), Selection):
                raise TypeError(
                    "members elements must inherit from Selection, not {}".format(type(members[0])))
        else:
            raise TypeError("members must be type list, not {}".format(type(members)))

        super().__init__(members=members, idx=idx, ids=ids)
        self._member_type = Selection

class SelectionDict(TrackedDict):
    def __init__(self, tdict):
        # check to make sure values are either inherited from
        # Selection or SelectionList
        pass

class Association(TrackedList):
    """Base class for associating data to a SelectionList.  """

    def __init__(self, members=None, association=None, idx=None, ids=None):
        super().__init__(members=members)
        self._association = association

    @property
    def association(self):
        return self._association

class AssociationType(object):
    """ Base class for association attributes. """

    def __init__(self, description=None):
        self._description = description

    @property
    def description(self):
        return self._description
