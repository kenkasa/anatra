#-------------------------------------------------------------------------------
proc print_title {} {
#-------------------------------------------------------------------------------

  puts "============================================================"
  puts ""
  puts "                    Free-Volume Analysis"
  puts ""
  puts "============================================================"
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc show_usage {arglist} {
#-------------------------------------------------------------------------------

  set help false
  set help [parse_arguments $arglist \
      "-h" "flag" $help]

  if {$help} {
    puts "Usage:"
    puts "anatra free_volume                                                                        \\"
    puts "  -stype        <structure file type>                                                     \\"
    puts "  -sfile        <structure file name.>                                                    \\"
    puts "  -tintype      <input trajectory file type>                                              \\"
    puts "  -tin          <input trajectory file name>                                              \\"
    puts "  -flist_traj   <trajectory file list (neccesary if tin is not specified)>                \\"
    puts "  -fhead        <header of output file name>                                              \\"
    puts "  -fprmtop      <Amber parm7 file. In current implementation, only parm7 is acceptable>   \\"
    puts "  -xsta         <value of initial window point>                                           \\"
    puts "  -dx           <window interval>                                                         \\"
    puts "  -ngrid        <# of windows>                                                            \\"
    puts "  -nins         <# of particle insertion per window>                                      \\"
    puts "  -sel0         <VMD selection>                                                           \\"
    puts "  -sel1         <VMD selection>                                                           \\"
    puts "  -target_selid <selection id for target species>                                         \\"
    puts "                (default: 0)                                                              \\"
    puts "  -center_selid <selection id for system center>                                          \\"
    puts "                (default: 1)                                                              \\"
    puts "  -coord_type   <coordinate type (z)>                                                     \\"
    puts "                (default: z)                                                              \\"
    puts "  -prep_only    <whether analysis is performed or not (true or false)>                    \\"
    puts "                (default: false)>"
    puts ""
    puts "Usage:"
    puts "anatra free_volume              \\"
    puts "  -stype        psf             \\"
    puts "  -sfile        complex.psf     \\"
    puts "  -tintype      dcd             \\"
    puts "  -tin          run.dcd         \\"
    puts "  -fhead        run             \\"
    puts "  -fprmtop      complex.prmtop  \\"
    puts "  -xsta         0.0             \\"
    puts "  -dx           0.5             \\"
    puts "  -ngrid        100             \\"
    puts "  -nins         1000            \\"
    puts "  -sel0         all             \\"
    puts "  -sel1         resname POPC    \\"
    puts "  -target_selid 0               \\"
    puts "  -center_selid 1               \\"
    puts "  -coord_type    z              \\"
    puts "  -prep_only    false"
    puts ""
    exit
  }

}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc define_optinfo {} {
#-------------------------------------------------------------------------------
  
  global opt

  set opt(fhead)         "run" 
  set opt(fprmtop)       "complex.prmtop" 
  set opt(xsta)          0.0
  set opt(dx)            0.5
  set opt(ngrid)         100
  set opt(nins)          1000
  set opt(target_selid)  0 
  set opt(center_selid)  0
  set opt(coord_type)   "z" 

}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc read_optinfo {arglist} {
#-------------------------------------------------------------------------------

  global opt

  set opt(fhead)        [parse_arguments $arglist \
      "-fhead"          "value" $opt(fhead)]
  set opt(fprmtop)      [parse_arguments $arglist \
      "-fprmtop"        "value" $opt(fprmtop)]
  set opt(xsta)         [parse_arguments $arglist \
      "-xsta"           "value" $opt(xsta)]
  set opt(dx)           [parse_arguments $arglist \
      "-dx"             "value" $opt(dx)]
  set opt(ngrid)        [parse_arguments $arglist \
      "-ngrid"          "value" $opt(ngrid)]
  set opt(nins)         [parse_arguments $arglist \
      "-nins"           "value" $opt(nins)]
  set opt(target_selid) [parse_arguments $arglist \
      "-target_selid"   "value" $opt(target_selid)]
  set opt(center_selid) [parse_arguments $arglist \
      "-center_selid"   "value" $opt(center_selid)]
  set opt(coord_type)   [parse_arguments $arglist \
      "-coord_type"     "value" $opt(coord_type)]
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc show_optinfo {} {
#-------------------------------------------------------------------------------

  global opt

  puts "<< option info >>"
  puts "fhead        = $opt(fhead)"
  puts "fprmtop      = $opt(fprmtop)"
  puts "xsta         = $opt(xsta)"
  puts "dx           = $opt(dx)"
  puts "ngrid        = $opt(ngrid)"
  puts "nins         = $opt(nins)"
  puts "target_selid = $opt(target_selid)"
  puts "center_selid = $opt(center_selid)"
  puts "coord_type   = $opt(coord_type)"
  puts ""

}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc analyze {} {
#-------------------------------------------------------------------------------
  # in
  global str
  global traj
  global seltxt
  global opt
  global common

  global sel 

  set anatra_path $::env(ANATRA_PATH);list
  set fort       "${anatra_path}/f90/bin/free_volume.x";list

  # read trajectory
  #
  puts ""
  puts "--------------------"
  puts " Read trajectory"
  puts "--------------------"
  puts ""

  #set mol 0;
  #read_traj $mol $str(stype) $str(sfile) $traj(tintype) $traj(tin) $traj(stride)
  #set nf   [molinfo $mol get numframes]
  set mol [mol load $str(stype) "$str(sfile)"]
  set nsel $seltxt(nsel)

  # setup selection
  #
  puts ""
  puts "--------------------"
  puts " Setup selection"
  puts "--------------------"
  puts ""
  for {set isel 0} {$isel < $seltxt(nsel)} {incr isel} {
    puts [format "selection %5d : %s" $isel $seltxt($isel)]
    set sel($isel) [atomselect $mol "$seltxt($isel)"]
  } 

  set rnam [$sel(0) get resname]
  set res  [$sel(0) get resid]
  set mass [$sel(0) get mass]
  set anam [$sel(0) get name]
  set chg  [$sel(0) get charge]
  set ind  [$sel(0) get index]
  set segn [$sel(0) get segname]
  set natm [llength $res]

  # Convert 
  #
  puts ""
  puts "--------------------"
  puts " Start analysis"
  puts "--------------------"
  puts ""

  set nf [molinfo $mol get numframes] 
  puts ""

  puts ">> Start Free-Volume calculation"
  puts ""

  set ffvinp   [format "%s.fv.inp" $opt(fhead) ]
  set ffvout   [format "%s.fv.out" $opt(fhead) ]

  for {set isel 0} {$isel < $nsel} {incr isel} {
    set fmolinfo($isel) [format "%s.fv.%i.molinfo" $opt(fhead) $isel]
  }

  for {set isel 0} {$isel < $nsel} {incr isel} {
    set rnam [$sel($isel) get resname]
    set res  [$sel($isel) get resid]
    set mass [$sel($isel) get mass]
    set anam [$sel($isel) get name]
    set chg  [$sel($isel) get charge]
    set ind  [$sel($isel) get index]
    set segn [$sel($isel) get segname]
    set natm [llength $res]
    set nf   [molinfo $mol get numframes]
    set nres [llength [lsort -unique [$sel($isel) get residue]]]

    set f [open $fmolinfo($isel) "w"]
    for {set iatm 0} {$iatm < $natm} {incr iatm} {
      puts $f [format "%10d  %6s  %6s  %15.7f  %15.7f  %d  %6s  %3s" \
        [lindex $res  $iatm]           \
        [lindex $rnam $iatm]           \
        [lindex $anam $iatm]           \
        [lindex $mass $iatm]           \
        [lindex $chg  $iatm]           \
        [expr [lindex $ind $iatm] + 1] \
        [lindex $segn $iatm]           \
        "END"]
    } 
    close $f
  }

  set ntraj [llength $traj(tin)]

  set st $opt(target_selid)
  set sc $opt(center_selid) 

  set f [open $ffvinp "w"]
  puts $f " &input_param"
  puts $f "   flist_traj = \"$traj(flist_traj)\""

  if {$ntraj > 0} {
    puts $f "   ftraj ="
    for {set i 0} {$i < $ntraj} {incr i} {
      set t [lindex $traj(tin) $i]
      puts -nonewline $f "    \"$t\" "
    }
  }
  puts $f "   fprmtop = \"$opt(fprmtop)\""
  puts $f ""
  puts $f " /"
  puts $f " &output_param"
  puts $f "   fhead   = \"$opt(fhead)\""
  puts $f " /"

  puts $f " &trajopt_param"
  puts $f "   molinfo = \"$fmolinfo($st)\" \"$fmolinfo($sc)\""
  puts $f " /"

  puts $f " &option_param"
  #puts $f "   mode       = \"whole\" \"whole\""
  puts $f "   xsta       = $opt(xsta)"
  puts $f "   ngrid      = $opt(ngrid)"
  puts $f "   dx         = $opt(dx)"
  puts $f "   coord_type = \"$opt(coord_type)\""
  puts $f "   nins       = $opt(nins)"
  puts $f " /"

  close $f

  #animate write dcd $fdcdtmp beg 0 end -1 waitfor all sel $sel(0) $mol

  if {!$common(prep_only)} {
    puts "Free-Volume is calculated with ANATRA fortran program:"
    puts "$fort ..."
    puts "=== INPUT ==="
    set content [exec cat $ffvinp]
    puts $content
    puts "============="
    exec $fort $ffvinp > $ffvout
    puts ""
    puts "=== OUTPUT ==="
    set content [exec cat $ffvout]
    puts $content
  }
  puts "=============="
  puts ">> Finished"

  exit
}
#-------------------------------------------------------------------------------
