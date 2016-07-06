proc draw_interactions { mol pairs_file color radius} {

    # read in contact pairs
    set file [open $pairs_file]

    # read in each line from the file
    while {[gets $file cline] >= 0} {
        # split on whitespace
        set inds [split $cline " "]
        # set variables from each row
        set ind1 [lindex $inds 0]
        set ind2 [lindex $inds 1]

        # make two selections for each atom
        set sel1 [atomselect top "serial $ind1"]
        set sel2 [atomselect top "serial $ind2"]
        # assign the coordinates to two variables
        lassign [$sel1 get {x y z}] start
        lassign [$sel2 get {x y z}] end

        # set the color for the graphics
        graphics $mol color $color
        # draw the cylinder
        graphics $mol cylinder $start $end radius $radius
    }
    close $file
}
