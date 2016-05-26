from mast.TrackedStructures import *

tms = [TrackedMember(idx=i) for i in range(5)]

tl = TrackedList(tms)

tlsel = Selection(container=tl, sel=0)
