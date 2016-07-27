
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
