import copy


__all__ = ['Angle', 'AngleType', 'Atom', 'AtomList', 'Bond', 'BondType',
           'ChiralFrame', 'Cmap', 'CmapType', 'Dihedral', 'DihedralType',
           'DihedralTypeList', 'Improper', 'ImproperType', 'MultipoleFrame',
           'OutOfPlaneBend', 'PiTorsion', 'Residue', 'ResidueList',
           'StretchBend', 'StretchBendType', 'TorsionTorsion',
           'TorsionTorsionType', 'TrigonalAngle', 'TrackedList', 'UreyBradley',
           'OutOfPlaneBendType', 'NonbondedException', 'NonbondedExceptionType',
           'AmoebaNonbondedExceptionType', 'AcceptorDonor', 'Group', 'AtomType',
           'NoUreyBradley', 'ExtraPoint', 'TwoParticleExtraPointFrame',
           'ThreeParticleExtraPointFrame', 'OutOfPlaneExtraPointFrame',
           'RBTorsionType', 'UnassignedAtomType']

# Create the AKMA unit system which is the unit system used by Amber and CHARMM

scale_factor = u.sqrt(1/u.kilocalories_per_mole * (u.daltons * u.angstroms**2))
scale_factor = scale_factor.value_in_unit(u.picoseconds)
akma_time_unit = u.BaseUnit(u.picosecond_base_unit.dimension, 'akma time',
                            symbol='aks')
akma_time_unit.define_conversion_factor_to(u.picosecond_base_unit, scale_factor)
akma_unit_system = u.UnitSystem([
        u.angstrom_base_unit, u.dalton_base_unit, akma_time_unit,
        u.elementary_charge_base_unit, u.kelvin_base_unit,
        u.mole_base_unit, u.radian_base_unit]
)

def _strip_units(value, unit=None):
    """
    Strips any units from the given value by casting them into the AKMA unit
    system (or the requested unit)
    """
    if u.is_quantity(value):
        # special-case angles, since pure angles are always in degrees
        if unit is None:
            if value.unit.is_compatible(u.degrees):
                return value.value_in_unit(u.degrees)
            return value.value_in_unit_system(akma_unit_system)
        else:
            return value.value_in_unit(unit)
    return value

def _getstate_with_exclusions(exclusions=None):
    """ Serializes based on all attributes except requested exclusions

    Parameters
    ----------
    exclusions : list of str, optional
        List of all attributes to exclude from serialization (should be
        descriptors and 'list'). Default is None

    Notes
    -----
    If exclusions is None, it defaults to excluding 'list'. If this is not
    desired, set exclusions to the empty list.
    """
    if exclusions is None:
        exclusions = ['list']
    def __getstate__(self):
        return {key : val for (key, val) in iteritems(self.__dict__)
                    if key not in exclusions}
    return __getstate__


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

    """

    def __init__(self, idx=None):
        self._changed = False
        if idx is None:
            self._idx = None
        elif isinstance(idx, int):
            if idx < 0:
                raise ValueError("idx can't be negative")
            else:
                self._idx = idx
        else:
            raise TypeError("idx must be a positive integer not type {}".format(type(idx)))

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


class Atom(TrackedMember):
    """
    An atom. Only use these as elements in AtomList instances, since AtomList
    will keep track of when indexes and other stuff needs to be updated. All
    parameters are optional.

    Parameters
    ----------
    atomic_number : ``int``
        The atomic number of this atom
    name : ``str``
        The name of this atom
    type : ``str``
        The type name of this atom
    charge : ``float``
        The partial atomic charge of this atom in fractions of an electron
    mass : ``float``
        The atomic mass of this atom in daltons
    nb_idx : ``int``
        The nonbonded index. This is a pointer that is relevant in the context
        of an Amber topology file and identifies its Lennard-Jones atom type
    solvent_radius : ``float``
        The intrinsic solvation radius of this atom.
    screen : ``float``
        The Generalized Born screening factor for this atom.
    occupancy : ``float``
        The occupancy of the atom (see PDB file)
    bfactor : ``float``
        The B-factor of the atom (see PDB file)
    altloc : ``str``
        Alternate location indicator (see PDB file)

    Other Parameters
    ----------------
    list : :class:`AtomList`
        The AtomList that this atom belongs to. If None, this atom does not
        belong to any list. This can be any iterable, but should typically be an
        AtomList or None
    tree : ``str``
        The tree chain identifier assigned to this atom. Relevant in the context
        of an Amber topology file, and not used for very much.
    join : ``int``
        The 'join` property of atoms stored in the Amber topology file. At the
        time of writing this class, `join` is unused, but still oddly required.
        Add support for future-proofing
    irotat : ``int``
        The `irotat` property of atoms stored in the Amber topology file.
        Unused, but included for future-proofing.
    number : ``int``
        The serial number given to the atom (see PDB file)
    rmin : ``float``
        The Rmin/2 Lennard-Jones parameter for this atom. Default evaluates to 0
    epsilon : ``float``
        The epsilon (well depth) Lennard-Jones parameter for this atom. Default
        evaluates to 0
    rmin14 : ``float``
        The Rmin/2 Lennard-Jones parameter for this atom in 1-4 interactions.
        Default evaluates to 0
    epsilon14 : ``float``
        The epsilon (well depth) Lennard-Jones parameter for this atom in 1-4
        interactions. Default evaluates to 0

    Other Attributes
    ----------------
    element : ``int``
        This is an alias for atomic_number
    atom_type : :class:`AtomType`
        In some cases, "type" is an ambiguous choice of an integer serial number
        or a string descriptor. In this case, atom_type is an AtomType instance
        that disambiguates this discrepancy.
    anisou : numpy.ndarray(float64) (or list of floats)
        Anisotropic temperature scaling factors. This is a 6-element numpy array
        They are the 3x3 symmetric matrix elements U(1,1), U(2,2), U(3,3),
        U(1,2), U(1,3), U(2,3). If no factors available, it is None.
    idx : ``int``
        The index of this atom in the list. Set to -1 if this atom is not part
        of a list or the index cannot otherwise be determined (i.e., if the
        containing list does not support indexing its members)
    residue : :class:`Residue`
        The Residue that this atom belongs to. This is assigned when this atom
        is passed to `Residue.add_atom` -- see below for more information. Until
        it is set there, it is None
    other_locations : ``dict`` of :class:`Atom`
        A dict of Atom instances that represent alternate conformers of this
        atom. The keys are the `altloc` characters for those Atoms.
    bonds : ``list`` of :class:`Bond`
        list of Bond objects in which this atom is a member. This attribute
        should not be modified.
    angles : ``list`` of :class:`Angle`
        list of Angle objects in which this atom is a member. This attribute
        should not be modified.
    dihedrals : ``list`` of :class:`Dihedral`
        list of Dihedral objects in which this atom is a member. This attribute
        should not be modified.
    urey_bradleys : ``list`` of :class:`UreyBradley`
        list of UreyBradley objects in which this atom is a member (CHARMM,
        AMOEBA). This attribute should not be modified.
    impropers : ``list`` of :class:`Improper`
        list of Improper objects in which the atom is a member (CHARMM). This
        attribute should not be modified.
    cmaps : ``list`` of :class:`Cmap`
        list of Cmap objects in which the atom is a member (CHARMM, AMOEBA).
        This attribute should not be modified.
    tortors : ``list`` of :class:`TorsionTorsion`
        list of TorsionTorsion objects in which the atom is a member (AMOEBA).
        This attribute should not be modified.
    bond_partners : ``list`` of :class:`Atom`
        list of Atoms to which this atom is bonded. Do not modify this
        attribute -- it will almost certainly not do what you think it will
    angle_partners : ``list`` of :class:`Atom`
        list of Atoms to which this atom forms an angle, but not a bond. Do not
        modify this attribute -- it will almost certainly not do what you think
        it will
    dihedral_partners : ``list`` of :class:`Atom`
        list of Atoms to which this atom forms an dihedral, but not a bond or
        angle. Do not modify this attribute -- it will almost certainly not do
        what you think it will
    tortor_partners : ``list`` of :class:`Atom`
        list of Atoms to which this atom forms a coupled Torsion-Torsion, but
        not a bond or angle (AMOEBA). Do not modify this attribute -- it will
        almost certainly not do what you think it will
    exclusion_partners : ``list`` of :class:`Atom`
        list of Atoms with which this atom is excluded, but not bonded, angled,
        or dihedraled to. Do not modify this attribute -- it will almost
        certainly not do what you think it will
    marked : ``int``
        Mainly for internal use, it is used to indicate when certain atoms have
        been "marked" when traversing the bond network identifying topological
        features (like molecules and rings)
    children : ``list`` of :class:`ExtraPoint`
        This is the list of "child" ExtraPoint objects bonded to this atom
    number : ``int``
        The serial number of the atom in the input structure (e.g., PDB file).
        If not present in the original structure, a default value of -1 is used
    rmin : ``float``
        The Rmin/2 Lennard-Jones parameter for this atom. Default value is 0.
        If not set, it is taken from the `atom_type` attribute (if available)
    epsilon : ``float``
        The epsilon (well depth) Lennard-Jones parameter for this atom. Default
        value is 0. If not set, it is taken from the `atom_type` attribute (if
        available)
    rmin_14 : ``float``
        The Rmin/2 L-J parameter for 1-4 pairs. Default value is `rmin` (see
        above). If not set, it is taken from the `atom_type` attribute (if
        available).
    epsilon_14 : ``float``
        The epsilon L-J parameter for 1-4 pairs. Default value is `epsilon` (see
        above). If not set, it is taken from the `atom_type` attribute (if
        available).

    Possible Attributes
    -------------------
    xx : ``float``
        The X-component of the position of this atom. Only present if
        coordinates have been loaded. Otherwise raises AttributeError
    xy : ``float``
        The Y-component of the position of this atom. Only present if
        coordinates have been loaded. Otherwise raises AttributeError
    xz : ``float``
        The Z-component of the position of this atom. Only present if
        coordinates have been loaded. Otherwise raises AttributeError
    vx : ``float``
        The X-component of the velocity of this atom. Only present if the
        velocities have been loaded. Otherwise raises AttributeError
    vy : ``float``
        The Y-component of the velocity of this atom. Only present if the
        velocities have been loaded. Otherwise raises AttributeError
    vz : ``float``
        The Z-component of the velocity of this atom. Only present if the
        velocities have been loaded. Otherwise raises AttributeError
    type_idx : ``int``
        The AMOEBA atom type index. Only present if initialized with the AMOEBA
        force field. Otherwise raises AttributeError.
    class_idx : ``int``
        The AMOEBA class type index Only present if initialized with the AMOEBA
        force field. Otherwise raises AttributeError.
    multipoles : ``list(float)``
        The list of the 10 multipole moments up through quadrupoles Only present
        if initialized with the AMOEBA force field. Otherwise raises
        AttributeError.
    polarizability : ``list(float)``
        The polarizability of the atom. Only present if initialized with the
        AMOEBA force field. Otherwise raises AttributeError.
    vdw_parent : :class:`Atom`
        In the AMOEBA force field, this is the parent atom for the van der Waals
        term
    vdw_weight : ``float``
        In the AMOEBA force field, this is the weight of the van der Waals
        interaction on the parent atom

    Notes
    -----
    The bond_partners, angle_partners, dihedral_partners, and exclusion_partners
    arrays are actually generated as properties by taking differences of sets
    and sorting them. As a result, accessing this attribute constantly can be
    less efficient than you would expect. Iterating over them in a loop requires
    minimal overhead. But if frequent access is needed and these sets are
    guaranteed not to change, you should save a reference to the object and use
    that instead.

    Binary comparisons are done by atom index and are used primarily for sorting
    sublists of atoms. The == operator is not defined, so Atom equality should
    be done using the `is` operator. This allows Atom instances to be hashable
    (and so used as `dict` keys and put in `set`s)

    Examples
    --------
    >>> a1 = Atom(name='CO', type='C', charge=0.5, mass=12.01)
    >>> a2 = Atom(name='OC', type='O', charge=-0.5, mass=12.01)
    >>> a1.bond_to(a2)
    >>> a1 in a2.bond_partners and a2 in a1.bond_partners
    True
    >>> a1.idx # Not part of a container
    -1

    This object also supports automatic indexing when it is part of a container

    >>> atom_list = []
    >>> atom_list.append(Atom(list=atom_list, name='CO', charge=0.5))
    >>> atom_list.append(Atom(list=atom_list, name='OC', charge=-0.5))
    >>> atom_list[0].idx
    0
    >>> atom_list[1].idx
    1
    """
    #===================================================

    def __init__(self, idx=None, atomic_number=0, name='', type='',
                 charge=None, mass=0.0, nb_idx=0, solvent_radius=0.0,
                 screen=0.0, tree='BLA', join=0.0, irotat=0.0, occupancy=0.0,
                 bfactor=0.0, altloc='', number=-1, rmin=None, epsilon=None,
                 rmin14=None, epsilon14=None):

        super().__init__(idx=idx)

        self.atomic_number = atomic_number
        self.name = name.strip()
        try:
            self.type = type.strip()
        except AttributeError:
            self.type = type
        self._charge = _strip_units(charge, u.elementary_charge)
        self.mass = _strip_units(mass, u.dalton)
        self.nb_idx = nb_idx
        self.solvent_radius = _strip_units(solvent_radius, u.angstrom)
        self.screen = screen
        self.tree = tree
        self.join = join
        self.irotat = irotat
        self.bfactor = bfactor
        self.altloc = altloc
        self.occupancy = occupancy
        self._bond_partners = []
        self._angle_partners = []
        self._dihedral_partners = []
        self._tortor_partners = []
        self._exclusion_partners = [] # For arbitrary/other exclusions
        self.residue = None
        self.marked = 0 # For setting molecules
        self.bonds, self.angles, self.dihedrals = [], [], []
        self.urey_bradleys, self.impropers, self.cmaps = [], [], []
        self.tortors = []
        self.other_locations = {} # A dict of Atom instances
        self.atom_type = UnassignedAtomType
        self.number = number
        self.anisou = None
        self._rmin = rmin
        self._epsilon = epsilon
        self._rmin14 = rmin14
        self._epsilon14 = epsilon14
        self.children = []

    #===================================================

    @classmethod
    def _copy(cls, item):
        new = cls(atomic_number=item.atomic_number, name=item.name,
                  type=item.type, charge=item.charge, mass=item.mass,
                  nb_idx=item.nb_idx, solvent_radius=item.solvent_radius,
                  screen=item.screen, tree=item.tree, join=item.join,
                  irotat=item.irotat, occupancy=item.occupancy,
                  bfactor=item.bfactor, altloc=item.altloc)
        new.atom_type = item.atom_type
        new.anisou = copy(item.anisou)
        for key in item.other_locations:
            new.other_locations[key] = copy(item.other_locations[key])
        _safe_assigns(new, item, ('xx', 'xy', 'xz', 'vx', 'vy', 'vz',
                      'type_idx', 'class_idx', 'multipoles', 'polarizability',
                      'vdw_parent', 'vdw_weight'))
        return new

    def __copy__(self):
        """ Returns a deep copy of this atom, but not attached to any list """
        return type(self)._copy(self)

    #===================================================

    @property
    def bond_partners(self):
        """ Go through all bonded partners """
        bp = set(self._bond_partners)
        for p in self._bond_partners:
            for c in p.children:
                bp.add(c)
        return sorted(list(bp))

    @property
    def angle_partners(self):
        """ List of all angle partners that are NOT bond partners """
        bp = set(self._bond_partners)
        ap = set(self._angle_partners) - bp
        toadd = set()
        for p in ap:
            for c in p.children:
                toadd.add(c)
        ap |= toadd
        return sorted(list(ap))

    @property
    def dihedral_partners(self):
        " List of all dihedral partners that are NOT angle or bond partners "
        bp = set(self._bond_partners)
        ap = set(self._angle_partners)
        dp = set(self._dihedral_partners) - ap - bp
        toadd = set()
        for p in dp:
            for c in p.children:
                toadd.add(c)
        dp |= toadd
        return sorted(list(dp))

    @property
    def tortor_partners(self):
        """
        List of all 1-5 partners that are NOT in angle or bond partners. This is
        *only* used in the Amoeba force field
        """
        bp = set(self._bond_partners)
        ap = set(self._angle_partners)
        dp = set(self._dihedral_partners)
        tp = set(self._tortor_partners) - dp - ap - bp
        toadd = set()
        for p in tp:
            for c in p.children:
                toadd.add(c)
        tp |= toadd
        return sorted(list(tp))

    @property
    def exclusion_partners(self):
        """
        List of all exclusions not otherwise excluded by bonds/angles/torsions
        """
        # A little expensive, but the only way to ensure this is completely
        # correct easily
        bp = set(self.bond_partners)
        ap = set(self.angle_partners)
        dp = set(self.dihedral_partners)
        tp = set(self.tortor_partners)
        ep = set(self._exclusion_partners) - tp - dp - ap - bp
        toadd = set()
        for p in ep:
            for c in p.children:
                toadd.add(c)
        ep |= toadd
        return sorted(list(ep))

    #===================================================

    # Various parameters that can be taken from the AtomType if not set on the
    # atom directly.

    @property
    def charge(self):
        if self._charge is None:
            if (self.atom_type is UnassignedAtomType or
                    self.atom_type.charge is None):
                return 0.0
            return self.atom_type.charge
        return self._charge

    @charge.setter
    def charge(self, value):
        self._charge = value

    @property
    def rmin(self):
        """ Lennard-Jones Rmin/2 parameter (the to Lennard-Jones radius) """
        if self._rmin is None:
            if (self.atom_type is UnassignedAtomType or
                    self.atom_type.rmin is None):
                return 0.0
            return self.atom_type.rmin
        return self._rmin

    @rmin.setter
    def rmin(self, value):
        """ Lennard-Jones Rmin/2 parameter (the Lennard-Jones radius) """
        self._rmin = value

    @property
    def sigma(self):
        """ Lennard-Jones sigma parameter -- directly related to Rmin """
        return self.rmin * 2**(-1/6) * 2

    @sigma.setter
    def sigma(self, value):
        self._rmin = value * 2**(1/6) / 2

    @property
    def epsilon(self):
        """ Lennard-Jones epsilon parameter (the Lennard-Jones well depth) """
        if self._epsilon is None:
            if (self.atom_type is UnassignedAtomType or
                    self.atom_type.epsilon is None):
                return 0.0
            return self.atom_type.epsilon
        return self._epsilon

    @epsilon.setter
    def epsilon(self, value):
        """ Lennard-Jones epsilon parameter (the Lennard-Jones well depth) """
        self._epsilon = value

    @property
    def rmin_14(self):
        """ The 1-4 Lennard-Jones Rmin/2 parameter """
        if self._rmin14 is None:
            if (self.atom_type is UnassignedAtomType or
                    self.atom_type.rmin_14 is None):
                return self.rmin
            return self.atom_type.rmin_14
        return self._rmin14

    @rmin_14.setter
    def rmin_14(self, value):
        """ The 1-4 Lennard-Jones Rmin/2 parameter """
        self._rmin14 = value

    @property
    def sigma_14(self):
        """ Lennard-Jones sigma parameter -- directly related to Rmin """
        if self._rmin14 is None:
            if (self.atom_type is UnassignedAtomType or
                    self.atom_type.rmin_14 is None):
                return self.sigma
            return self.atom_type.rmin_14 * 2**(-1/6) * 2
        return self._rmin14 * 2**(-1/6) * 2

    @sigma_14.setter
    def sigma_14(self, value):
        self._rmin14 = value * 2**(1/6) / 2

    @property
    def epsilon_14(self):
        """ The 1-4 Lennard-Jones epsilon parameter """
        if self._epsilon14 is None:
            if (self.atom_type is UnassignedAtomType or
                    self.atom_type.epsilon_14 is None):
                return self.epsilon
            return self.atom_type.epsilon_14
        return self._epsilon14

    @epsilon_14.setter
    def epsilon_14(self, value):
        """ The 1-4 Lennard-Jones epsilon parameter """
        self._epsilon14 = value

    #===================================================

    def nonbonded_exclusions(self, only_greater=True, index_from=0):
        """
        Returns the total number of nonbonded atom exclusions for this atom. The
        general rules for building the exclusion list is to include both
        exceptions AND exclusions (i.e., the Amber scaling of 1-4 interactions
        means that the 1-4 terms are excluded and a special pairlist is built to
        handle those exceptions).

        All atoms in the `_partners` arrays are nonbonded exclusions.

        Parameters
        ----------
        only_greater : ``bool``, optional
            If True (default), only atoms whose `idx` value is greater than this
            `Atom`s `idx` will be counted as an exclusion (to avoid double-
            counting exclusions). If False, all exclusions will be counted.
        index_from : ``int``, optional
            This is the index of the first atom, and is intended to be 0 (for C-
            and Python-style numbering, default) or 1 (for Fortran-style
            numbering, such as that used in the Amber and CHARMM topology files)

        Returns
        -------
        ``list of int``
            The returned list will be the atom indexes of the exclusion partners
            for this atom (indexing starts from ``index_from``)

        Notes
        -----
        If this instance's `idx` attribute evaluates to -1 -- meaning it is not
        in an AtomList -- a IndexError will be raised. If you have two extra
        points (i.e., those with atomic numbers of 0) bonded to each other, this
        routine may raise a ``RuntimeError`` if the recursion depth is exceeded.
        """
        if self.idx < 0:
            raise IndexError('Cannot find exclusions of an unindexed Atom')
        if only_greater:
            baseline = self.idx
        else:
            baseline = -1
        excl = []
        for atm in self.bond_partners:
            i = atm.idx + index_from
            if i > baseline:
                excl.append(i)
        for atm in self.angle_partners:
            i = atm.idx + index_from
            if i > baseline:
                excl.append(i)
        for atm in self.dihedral_partners:
            i = atm.idx + index_from
            if i > baseline:
                excl.append(i)
        for atm in self.tortor_partners:
            i = atm.idx + index_from
            if i > baseline:
                excl.append(i)
        for atm in self.exclusion_partners:
            i = atm.idx + index_from
            if i > baseline:
                excl.append(i)
        return sorted(excl)

    #===================================================

    # Make 'element' an alias for 'atomic_number'

    @property
    def element(self):
        return self.atomic_number
    @element.setter
    def element(self, value):
        self.atomic_number = value

    #===================================================

    @property
    def segid(self):
        return self.residue.segid
    @segid.setter
    def segid(self, value):
        self.residue.segid = value

    #===================================================

    def bond_to(self, other):
        """
        Log this atom as bonded to another atom.

        Parameters
        ----------
        other : :class:`Atom`
            An atom that will be added to `bond_partners`

        Notes
        -----
        This action adds `self` to `other.bond_partners`. Raises
        :class:`MoleculeError` if `other is self`
        """
        if isinstance(other, ExtraPoint):
            self.children.append(other)
        elif isinstance(self, ExtraPoint):
            other.children.append(self)
        if self is other:
            raise MoleculeError("Cannot bond atom to itself!")
        self._bond_partners.append(other)
        other._bond_partners.append(self)

    #===================================================

    def angle_to(self, other):
        """
        Log this atom as angled to another atom.

        Parameters
        ----------
        other : :class:`Atom`
            An atom that will be added to `angle_partners`

        Notes
        -----
        This action adds `self` to `other.angle_partners`. Raises
        :class:`MoleculeError` if `other is self`
        """
        if self is other:
            raise MoleculeError("Cannot angle an atom with itself!")
        self._angle_partners.append(other)
        other._angle_partners.append(self)

    #===================================================

    def dihedral_to(self, other):
        """
        Log this atom as dihedral-ed to another atom.

        Parameters
        ----------
        other : :class:`Atom`
            An atom that will be added to `dihedral_partners`

        Notes
        -----
        This action adds `self` to `other.dihedral_partners`. Raises
        :class:`MoleculeError` if `other is self`
        """
        if self is other:
            raise MoleculeError("Cannot dihedral an atom with itself!")
        self._dihedral_partners.append(other)
        other._dihedral_partners.append(self)

    #===================================================

    def tortor_to(self, other):
        """
        Log this atom as 1-5 partners to another atom

        Parameters
        ----------
        other : :class:`Atom`
            An atom that will be added to `tortor_partners`

        Notes
        -----
        This action adds `self` to `other.tortor_partners`. Raises
        :class:`MoleculeError` if `other is self`
        """
        if self is other:
            raise MoleculeError('Cannot coupled-dihedral atom to itself')
        self._tortor_partners.append(other)
        other._tortor_partners.append(self)

    #===================================================

    def exclude(self, other):
        """
        Add one atom to my arbitrary exclusion list

        Parameters
        ----------
        other : :class:`Atom`
            An atom that will be added to `exclusion_partners`.

        Notes
        -----
        This action adds `self` to `other.exclusion_partners`. Raises
        :class:`MoleculeError` if `other is self`
        """
        if self is other:
            raise MoleculeError("Cannot exclude an atom from itself")
        self._exclusion_partners.append(other)
        other._exclusion_partners.append(self)

    #===================================================

    def prune_exclusions(self):
        """
        For extra points, the exclusion partners may be filled before the bond,
        angle, dihedral, and tortor partners. Since we don't want memory of
        these exclusions if any of those topological features were to break, we
        want to *remove* those from the exclusion list. This function makes sure
        that nothing in the bond, angle, dihedral, and tortor lists appears in
        the exclusion list.
        """
        excludes = (set(self._exclusion_partners) - set(self._tortor_partners) -
                    set(self._dihedral_partners) - set(self._angle_partners) -
                    set(self._bond_partners))
        self._exclusion_partners = sorted(list(excludes))

    #===================================================

    # Comparisons are done by comparing indexes

    def __gt__(self, other):
        return self.idx > other.idx

    def __lt__(self, other):
        return self.idx < other.idx

    def __ge__(self, other):
        return not self < other

    def __le__(self, other):
        return not self > other

    def __repr__(self):
        start = '<Atom %s [%d]' % (self.name, self.idx)
        if self.residue is not None and hasattr(self.residue, 'idx'):
            return start + '; In %s %d>' % (self.residue.name, self.residue.idx)
        elif self.residue is not None:
            return start + '; In object %r>' % self.residue
        return start + '>'

    #===================================================

    # For pickleability

    def __getstate__(self):
        retval = dict(name=self.name, type=self.type, atom_type=self.atom_type,
                      _charge=self._charge, mass=self.mass, nb_idx=self.nb_idx,
                      solvent_radius=self.solvent_radius, screen=self.screen,
                      tree=self.tree, join=self.join, irotat=self.irotat,
                      bfactor=self.bfactor, altloc=self.altloc,
                      occupancy=self.occupancy, number=self.number,
                      anisou=self.anisou, _rmin=self._rmin,
                      _epsilon=self._epsilon, _rmin14=self._rmin14,
                      _epsilon14=self._epsilon14, children=self.children,
                      atomic_number=self.atomic_number,
        )
        for key in ('xx', 'xy', 'xz', 'vx', 'vy', 'vz', 'multipoles',
                    'type_idx', 'class_idx', 'polarizability', 'vdw_weight',
                    'weights', '_frame_type'):
            try:
                retval[key] = getattr(self, key)
            except AttributeError:
                continue

        return retval

    def __setstate__(self, d):
        self._bond_partners = []
        self._angle_partners = []
        self._dihedral_partners = []
        self._tortor_partners = []
        self._exclusion_partners = []
        self.residue = None
        self.marked = 0
        self.bonds, self.angles, self.dihedrals = [], [], []
        self.urey_bradleys, self.impropers, self.cmaps = [], [], []
        self.tortors = []
        self.other_locations = {}
        self.__dict__.update(d)

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
        elif isinstance(List, list) or isinstance(List, TrackedList):
            self._type = type(List[0])
            self._List = List
            self._changed = True
            # index the new list if it's members have idx is None
            if List:
                if List[0].idx is None:
                    self.index_members()
                    self._indexed = True
                else:
                    self._indexed = False
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

    # returns a deepcopy of the items, so they cannot be modified
    def __getitem__(self, index):
        """ Return a deepcopy of the member so it cannot be modified. """
        return copy.deepcopy(self._List[index])

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

class AtomList(TrackedList):
    """TrackedList of Atoms, that propogates signals to all Selections
    and Structures it is a part of if changed.
    """

    def __init__(self, List=None):
        if List is None:
            self._List = []
        super().__init__(self)
        print(type(self))
        self._type = type(self)



class Structure(TrackedList):
    """A tracked list of Atoms connected through covalent bonds, with
methods for making selections."""

    def __init__(self, atom_list=None):
        if atom_list is None:
            # setter sets the changed and indexed attributes as well
            self.List = []
        elif isinstance(atom_list, AtomList):
            self.List = atom_list.List.copy()
        else:
            raise TypeError("Constructor type {} is not None or AtomList".format(atom_list))

