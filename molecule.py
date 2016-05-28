from copy import copy
from itertools import product

import networkx as nx

from mast.datastructures import TrackedMember, TrackedList, Selection, SelectionList
import mast.unit as u

__all__ = ['Atom', 'Bond', 'Angle',
           'AtomList', 'BondList', 'AngleList',
           'MoleculeTopology', 'Molecule', 'MoleculeList']


### Directly from ParmEd
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

class _UnassignedAtomType(object):
    """
    This raises the appropriate exceptions (ParameterError) when you try to
    access its properties
    """
    _OBJ = None

    def __new__(cls):
        if cls._OBJ is None:
            cls._OBJ = super(_UnassignedAtomType, cls).__new__(cls)
        return cls._OBJ

    def __int__(self):
        raise ParameterError('Atom type is not defined')

    def __str__(self):
        raise ParameterError('Atom type is not defined')

    def __eq__(self, other):
        return isinstance(other, _UnassignedAtomType) # Behave like a singleton

    def __reduce__(self):
        return 'UnassignedAtomType'

UnassignedAtomType = _UnassignedAtomType()
# Make sure it's a singleton
assert UnassignedAtomType is _UnassignedAtomType(), "Not a singleton"

def _delete_from_list(list, item):
    """
    Deletes a requested item from a list. If the item does not exist in the
    list, a ValueError is raised

    Parameters
    ----------
    list : ``list``
        The list from which an item will be deleted
    item : ``object``
        The object to delete from the list
    """
    list.pop(list.index(item))

def _safe_assigns(dest, source, attrs):
    """
    Shallow-copies all requested attributes from `source` to `dest` if they
    are present. If not present, nothing is done
    """
    for attr in attrs:
        if not hasattr(source, attr): continue
        myattr = getattr(source, attr)
        setattr(dest, attr, myattr)

#### End of unkown function ParmEd stuff


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

    def __init__(self, atom=None, idx=None, ids=None,
                 coordinate=None, velocity=None,
                 atomic_number=0, name='', parmed_type='',
                 charge=None, mass=0.0, nb_idx=0, solvent_radius=0.0,
                 screen=0.0, tree='BLA', join=0.0, irotat=0.0, occupancy=0.0,
                 bfactor=0.0, altloc='', number=-1, rmin=None, epsilon=None,
                 rmin14=None, epsilon14=None):

        # copy an atom in the arguments if given
        if atom:
            self = copy(atom)
        # otherwise make another fresh one
        else:
            super().__init__(idx=idx, ids=ids)

            # Mast defined attributes

            self.coordinate = coordinate
            # check to see if velocity tuple is correct
            self.velocity = velocity

            # ParmEd defined attributes
            self.atomic_number = atomic_number
            self.name = name.strip()
            try:
                self.parmed_type = parmed_type.strip()
            except AttributeError:
                self.parmed_type = parmed_type
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

    def __copy__(self):
        new_atom = Atom()
        new_atom._idx = self.idx
        new_atom._ids = self.ids
        new_atom._coordinate = self._coordinate
        new_atom._velocity = self._velocity

        new_atom.atomic_number=self.atomic_number
        new_atom.name=self.name
        new_atom.parmed_type=self.parmed_type
        new_atom.charge=self.charge
        new_atom.mass=self.mass
        new_atom.nb_idx=self.nb_idx
        new_atom.solvent_radius=self.solvent_radius
        new_atom.screen=self.screen
        new_atom.tree=self.tree
        new_atom.join=self.join
        new_atom.irotat=self.irotat
        occupancy=self.occupancy
        new_atom.bfactor=self.bfactor
        new_atom.altloc=self.altloc
        new_atom.atom_type = self.atom_type
        # is this okay
        new_atom.anisou = copy(self.anisou)
        for key in self.other_locations:
            new_atom.other_locations[key] = copy(self.other_locations[key])

        _safe_assigns(new_atom, self, ('type_idx', 'class_idx',
                                       'multipoles', 'polarizability',
                                       'vdw_parent', 'vdw_weight'))
        return new_atom

    # I implement a shallow copy
    # def __copy__(self):
    #     """ Returns a deep copy of this atom, but not attached to any list """
    #     return type(self)._copy(self)

    #===================================================

    @property
    def coordinate(self):
        return self._coordinate

    @coordinate.setter
    def coordinate(self, coordinate):
        if not coordinate or not all(coordinate):
            self._coordinate = None
        elif isinstance(coordinate, tuple):
            if len(coordinate) != 3:
                raise TypeError(
                    "Coordinate tuples must be length 3, not {}".format(len(coordinate)))
            elif not all([isinstance(dim, float) for dim in coordinate]):
                raise TypeError(
                    "Coordinate tuple values must be type float, not {}".format(
                                                          type(coordinate[0])))
            else:
                self._coordinate = coordinate
        else:
            raise TypeError(
                "Coordinate tuples must be type tuple, not {}".format(type(coordinate)))

    # convenience names
    coord = x = coordinate

    @property
    def xx(self):
        if self._coordinate:
            return self._coordinate[0]
        else:
            return None
    @property
    def xy(self):
        if self._coordinate:
            return self._coordinate[1]
        else:
            return None
    @property
    def xz(self):
        if self._coordinate:
            return self._coordinate[2]
        else:
            return None

    @property
    def velocity(self):
        return self._velocity

    @velocity.setter
    def velocity(self, velocity):
        if not velocity:
            self._velocity = None
        elif isinstance(velocity, tuple):
            if len(velocity) != 3:
                raise TypeError(
                    "Velocity tuples must be length 3, not {}".format(len(velocity)))
            elif not all([isinstance(dim, float) for dim in velocity]):
                raise TypeError(
                    "Velocity tuple values must be type float, not {}".format(
                                                          type(velocity[0])))
            else:
                self._velocity = velocity
        else:
            raise TypeError(
                "Velocity tuples must be type tuple, not {}".format(type(velocity)))

    # convenience names
    vel = v = velocity

    @property
    def vx(self):
        if self._velocity:
            return self._velocity[0]
        else:
            return None
    @property
    def vy(self):
        if self._velocity:
            return self._velocity[1]
        else:
            return None
    @property
    def vz(self):
        if self._velocity:
            return self._velocity[2]
        else:
            return None

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
        start = '<Atom {0} [{1}]'.format(self.name, self.idx)
        if self.residue is not None and hasattr(self.residue, 'idx'):
            return start + '; In %s %d>' % (self.residue.name, self.residue.idx)
        elif self.residue is not None:
            return start + '; In object %r>' % self.residue
        return start + '>'

    #===================================================

    # For pickleability

    def __getstate__(self):
        retval = dict(name=self.name, parmed_type=self.parmed_type, atom_type=self.atom_type,
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



# Bond
class Bond(Selection):
    """
    A covalent bond connecting two atoms.

    Parameters
    ----------
    atom1 : :class:`Atom`
        The first atom involved in the bond
    atom2 : :class:`Atom`
        The other atom involved in the bond
    bond_type : :class:`BondType`
        The bond type that defines the parameters for this bond

    Notes
    -----
    You can test whether an :class:`Atom` is contained within the bond using the
    `in` operator. A `MoleculeError` is raised if `atom1` and `atom2` are identical.
    This bond instance is `append`ed to the `bonds` list for both `atom1` and
    `atom2` and is automatically removed from those lists upon garbage
    collection

    Examples
    --------
    >>> a1, a2 = Atom(), Atom()
    >>> bond = Bond(a1, a2)
    >>> a1 in bond and a2 in bond
    True
    """

    def __init__(self, atoms, atom1_idx, atom2_idx):
        """ Bond constructor """

        # check to see if container is correct type
        if not isinstance(atoms, AtomList):
            raise TypeError("atoms must be type AtomList, not type {}".format(type(atoms)))

        # check that the indices are the correct type and values
        if not (isinstance(atom1_idx, int) or isinstance(atom2_idx, int) ):
            raise TypeError(
                "atom indices must be type int, not {0} and {1}".format(type(atom1_idx), type(atom2_idx)))

        elif (atom1_idx < 0 or atom2_idx < 0):
                raise ValueError(
                    "atom indices can't be negative, given {0} and {1}".format(atom1_idx, atom2_idx))

        # inherited contructor
        super().__init__(container=atoms, sel=[atom1_idx, atom2_idx])


    def __copy__(self):
        return Bond(self._container, self.atom1.idx, self.atom2.idx)
    
    @property
    def atom1(self):
        return self[0]

    @property
    def atom2(self):
        return self[1]

    def __contains__(self, thing):
        """ Quick and easy way to see if an Atom is in this Bond """
        return thing is self.atom1 or thing is self.atom2

    def delete(self):
        """
        Deletes this bond from the atoms that make it up. This method removes
        this bond from the `bonds` list for both atom1 and atom2 as well as
        removing atom1 from atom2.bond_partners and vice versa.
        """
        _delete_from_list(self.atom1.bonds, self)
        _delete_from_list(self.atom2.bonds, self)
        _delete_from_list(self.atom1._bond_partners, self.atom2)
        _delete_from_list(self.atom2._bond_partners, self.atom1)

        self.atom1 = self.atom2 = self.bond_type = None

    def __repr__(self):
        return '<{0} {1}--{2}; bond_type={3}>'.format(type(self).__name__,
                self.atom1, self.atom2, type(self))

# Angle
class Angle(Selection):
    """
    A valence angle between 3 atoms separated by two covalent bonds.

    Parameters
    ----------
    atom1 : :class:`Atom`
        An atom one end of the valence angle
    atom2 : :class:`Atom`
        The atom in the middle of the valence angle bonded to both other atoms
    atom3 : :class:`Atom`
        An atom on the other end of the valence angle to atom1
    angle_type : :class:`AngleType`
        The AngleType object containing the parameters for this angle

    Notes
    -----
    An Angle can contain bonds or atoms. A bond is contained if it exists
    between atoms 1 and 2 or atoms 2 and 3.

    Examples
    --------
    >>> a1, a2, a3 = Atom(), Atom(), Atom()
    >>> angle = Angle(a1, a2, a3)
    >>> Bond(a1, a2) in angle and Bond(a3, a2) in angle
    True
    >>> a1 in angle and a2 in angle and a3 in angle
    True
    >>> Bond(a1, a3) in angle # this is not part of the angle definition
    False
    """

    def __init__(self, atom1, atom2, atom3, angle_type=None):
        # Make sure we're not angling me to myself
        if atom1 is atom2 or atom1 is atom3 or atom2 is atom3:
            raise MoleculeError('Cannot angle atom to itself!')
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        # Log these angles in each atom
        atom1.angles.append(self)
        atom2.angles.append(self)
        atom3.angles.append(self)
        # Load the force constant and equilibrium angle
        self.angle_type = angle_type
        atom1.angle_to(atom2)
        atom1.angle_to(atom3)
        atom2.angle_to(atom3)
        self.funct = 1

    def __contains__(self, thing):
        """ Quick and easy way to see if an Atom or a Bond is in this Angle """
        if isinstance(thing, Atom):
            return (thing is self.atom1 or thing is self.atom2 or
                    thing is self.atom3)
        return ((self.atom1 in thing and self.atom2 in thing) or
                (self.atom2 in thing and self.atom3 in thing))

    def delete(self):
        """
        Deletes this angle from the atoms that make it up. This method removes
        this angle from the `angles` list for atom1, atom2, and atom3 as well as
        removing each atom form each others' angle partner lists
        """
        _delete_from_list(self.atom1.angles, self)
        _delete_from_list(self.atom2.angles, self)
        _delete_from_list(self.atom3.angles, self)

        _delete_from_list(self.atom1._angle_partners, self.atom2)
        _delete_from_list(self.atom1._angle_partners, self.atom3)
        _delete_from_list(self.atom2._angle_partners, self.atom1)
        _delete_from_list(self.atom2._angle_partners, self.atom3)
        _delete_from_list(self.atom3._angle_partners, self.atom1)
        _delete_from_list(self.atom3._angle_partners, self.atom2)

        self.atom1 = self.atom2 = self.atom3 = self.angle_type = None

    def __repr__(self):
        return '<%s %r--%r--%r; angle_type=%r>' % (type(self).__name__,
                self.atom1, self.atom2, self.atom3, self.angle_type)



class Molecule(TrackedMember):
    """An object containing minimum information necessary to specify a 3D
molecule in space with internal coordinates. Also contains 3D
coordinates in each atom.

    atoms : AtomList
    bonds : BondList

    angles : AngleList :: stub, only useful for if atoms don't have coordinates or if you
care about parameters


    """

    def __init__(self, atoms=None, bonds=None, angles=None,
                 idx=None, ids=None):
        super().__init__(idx=idx, ids=ids)
        if atoms is None:
            self._atoms = AtomList()
        elif isinstance(atoms, AtomList):
            self._atoms = atoms
        else:
            raise TypeError(
                "Constructor argument for atoms {} is not None or AtomList".format(atoms))

        if bonds is None:
            self._bonds = BondList()
        elif isinstance(bonds, BondList):
            # check all the bonds are between atoms that are in this molecule
            self._bonds = bonds
        else:
            raise TypeError(
                "Constructor argument for bonds {} is not None or BondList".format(bonds))

        if not angles:
            self._angles = None
        elif isinstance(angles, AngleList):
            # check that all angles are between bonds of this molecule
            self._angles = angles
        else:
            raise TypeError(
                "Constructor argument for bonds {} is not None or AngleList".format(angles))

        # construct topology
        if self._bonds and self._atoms:
            self._topology = MoleculeTopology(bonds=self._bonds)

    def __copy__(self):
        return Molecule(atoms=self.atoms, bonds=self.bonds,
                        angles=self.angles, idx=self.idx, ids=self.ids)

    @property
    def atoms(self):
        return self._atoms

    @property
    def bonds(self):
        return self._bonds

    @property
    def angles(self):
        return self._angles

    @property
    def topology(self):
        return self._topology

    top = topology

    def overlaps(self, other):
        """Check whether this molecule overlaps with another.
        Checks whether any two atoms in each molecule have the same coordinates.

        bool : returns True if any overlaps detected

        """
        # TODO replace with assert?? Then you would lose the Error types
        if __debug__:
            if not isinstance(other, Molecule):
                raise TypeError("Other must be type Molecule, not {}".format(type(other)))

        pairs = product(self.atoms, other.atoms)
        try:
            pair = next(pairs)
        # if it is empty no overlaps
        except StopIteration:
            return False
        flag = True
        while flag:
            if pair[0] == pair[1]:
                return True
            else:
                try:
                    pair = next(pairs)
                except StopIteration:
                    flag = False
        return False




class BondList(SelectionList):
    """ TrackedList of Bonds"""

    def __init__(self, members=None):
        if not members:
            self._members = []
        else:
            if issubclass(type(members), SelectionList) or isinstance(members, list):
                if isinstance(members[0], Bond):
                    super().__init__(members=members)
                else:
                    raise TypeError(
                        "members elements must be type Bond, not {}".format(type[members[0]]))
            else:
                raise TypeError(
                    "members must be type list or SelectionList, not {}".format(type(members)))

        self._member_type = Bond

    @property
    def bonds(self):
        return self._members

    @property
    def atoms(self):
        atoms = []
        for bond in self.bonds:
            if not bond.atom1 in atoms:
                atoms.append(bond.atom1)
            if not bond.atom2 in atoms:
                atoms.append(bond.atom2)

        return atoms


class AngleList(SelectionList):
    """ TrackedList of Angles"""

    def __init__(self, members=None):
        if members is None:
            self._members = []
        super().__init__(members=members)

        self._member_type = Angle

class MoleculeList(TrackedList):
    """ TrackedList of Molecules."""

    def __init__(self, members=None):
        if not members:
            self._members = []
        else:
            if isinstance(members, list):
                if isinstance(members[0], Molecule):
                    super().__init__(members=members)
                else:
                    raise TypeError(
                        "members elements must be type {0}, not {1}".format(Molecule, type(members[0])))
            else:
                raise TypeError(
                    "members must be type {0}, not {1}".format(list, type(members)))

        self._member_type = Molecule

        @property
        def molecules(self):
            return self._members

class MoleculeTopology(BondList):
    """Class to store a molecular topology which is a set of connected
bonds.

    atoms : AtomList
    bonds : BondList
    bond_graph : networkx.Graph

    """

    def __init__(self, bonds=None):
        super().__init__(bonds)

        # put the bonds in as edges to a graph
        self._bond_graph = nx.Graph()
        for bond in self.bonds:
            self._bond_graph.add_edge(bond.atom1, bond.atom2)
        # check to make sure all the bonds are connected
        num_components = len(list(nx.connected_component_subgraphs(self._bond_graph, copy=False)))
        if num_components != 1:
            raise ValueError(
                "All bonds must be connected, there are {0} graph components in the BondsList {1}"\
                   .format(num_components, bonds))

    @property
    def bond_graph(self):
        return self._bond_graph

class AtomList(TrackedList):
    """TrackedList of Atoms """

    def __init__(self, members=None):
        if not members:
            self._members = []
        else:
            if issubclass(type(members), TrackedList) or isinstance(members, list):
                if isinstance(members[0], Atom):
                    super().__init__(members=members)
                else:
                    raise TypeError(
                        "members elements must be type Atom, not {}".format(type(members[0])))
            else:
                raise TypeError("members must be type list or TrackedList, not {}".format(type(members)))

        self._member_type = Atom

    @property
    def atoms(self):
        return self._members
