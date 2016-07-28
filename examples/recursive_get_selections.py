from mast.selection import SelectionMember, Selection

container = [SelectionMember('a'), SelectionMember('b'), SelectionMember('c')]
selection = Selection(container, [0])
meta_selection = Selection([selection], [0])

selections_0 = container[0].get_selections(level=0)
print(selections_0)

selections_1 = container[0].get_selections(level=1)
print(selections_1)

selections_None = container[0].get_selections(level=None)
print(selections_None)
