#-------------------------------------------------------------------------------
proc print_title {} {
#-------------------------------------------------------------------------------

  puts "============================================================"
  puts ""
  puts "            Get Species Information (spec) file"
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
    puts "anatra species_info                      \\"
    puts "  -stype   <structure file type>         \\"
    puts "  -sfile   <structure file name>         \\"
    puts "  -fhead   <header of output file name>  \\"
    puts "  -sel0    <VMD selection> (X=0,1,2...)  \\"
    puts "  -sel1    <VMD selection> (X=0,1,2...)"
    puts ""  
    puts "Usage:"
    puts "anatra species_info          \\"
    puts "  -stype   parm7             \\"
    puts "  -sfile   str.prmtop        \\"
    puts "  -tintype rst7              \\"
    puts "  -tin     inp.inpcrd        \\"
    puts "  -fhead   out               \\"
    puts "  -sel0    name C32  H2X H2Y \\"
    puts "  -sel1    water"
    puts ""
    exit
  }

}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc define_optinfo {} {
#-------------------------------------------------------------------------------
  
  global opt

  set opt(fhead)         "out"
}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
proc read_optinfo {arglist} {
#-------------------------------------------------------------------------------

  global opt

  set opt(fhead)      [parse_arguments $arglist \
      "-fhead"       "value" $opt(fhead)]
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc show_optinfo {} {
#-------------------------------------------------------------------------------

  global opt

  puts "<< option info >>"
  puts "fhead = $opt(fhead)"
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

  global sel 


  # check control parameter check
  #

  # read trajectory
  #
  #puts ""
  #puts "--------------------"
  #puts " Read trajectory"
  #puts "--------------------"
  #puts ""

  set mol 0;
  set mol [mol load $str(stype) $str(sfile) ]
  #read_traj $mol $str(stype) $str(sfile) $traj(tintype) $traj(tin) $traj(stride)
  #set nf   [molinfo $mol get numframes]
  #set nsel $seltxt(nsel)

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


  # Convert 
  #
  puts ""
  puts "--------------------"
  puts " Start analysis"
  puts "--------------------"
  puts ""

  set nsel $seltxt(nsel)
  for {set isel 0} {$isel < $nsel} {incr isel} {
    set fmolinfo($isel) [format "%s_%i.spec" $opt(fhead) $isel]
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

  puts "=============="
  puts ">> Finished"

  exit
}
#-------------------------------------------------------------------------------
