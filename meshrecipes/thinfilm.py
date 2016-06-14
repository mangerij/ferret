reset
set developer commands on
undo on
create brick x 20 y 20 z 12
create brick x 20 y 20 z 10
volume 2 move 0 0 11.0
compress all
merge surface all with surface all
volume 1 2 size 0.9 
#mesh size depends h-level and material DW width!
volume 1 2 scheme sweep
mesh volume 1 2
block 1 volume 1
block 2 volume 2

sideset 1 surface 1

sideset 2 surface 2

sideset 3 surface 3
sideset 4 surface 4
sideset 5 surface 5
sideset 6 surface 6

sideset 7 surface 8
sideset 8 surface 9
sideset 9 surface 10
sideset 10 surface 11
sideset 11 surface 12

block all element type hex8
set large exodus file off
export Genesis "/home/john/projects/ferret/mesh/./exodus_thinfilm_09_20_12_10.e" dimension 3 block all overwrite
