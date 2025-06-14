#-------------------------------------------------------------------------------
proc print_title {} {
#-------------------------------------------------------------------------------

  puts "============================================================"
  puts ""
  puts "                  Lipid-Order Analysis"
  puts ""
  puts "============================================================"
  puts "Remark: Current version only supports CHARMM forcefield"
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
    puts "anatra lipid_order                                                        \\"
    puts "  -stype      <structure file type>                                       \\"
    puts "  -sfile      <structure file name>                                       \\"
    puts "  -tintype    <input trajectory file type>                                \\"
    puts "  -tin        <input trajectory file name>                                \\"
    puts "  -flist_traj <trajectory file list (neccesary if tin is not specified)>  \\"
    puts "  -fhead     <header of output file name>                                 \\"
    puts "  -cindex    <index for analyzed carbons>                                 \\"
    puts "  -selX      <VMD selection> (X=0,1,2...)                                 \\"
    puts "  -prep_only <whether analysis is performed or not (true or false)>       \\"
    puts "             (default: false)>"
    puts ""  
    puts "Usage:"
    puts "anatra lipid_order                            \\"
    puts "  -stype     psf                              \\"
    puts "  -sfile     complex.psf                      \\"
    puts "  -tintype   dcd                              \\"
    puts "  -tin       inp.dcd                          \\"
    puts "  -fhead     run                              \\"
    puts "  -cindex    1 2 3                            \\"
    puts "  -sel0      name C32  H2X H2Y and segid MEMB \\"
    puts "  -sel1      name C33  H3X H3Y and segid MEMB \\"
    puts "  -sel2      name C34  H4X H4Y and segid MEMB \\"
    puts "  -prep_only false"
    puts ""
    exit
  }

}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc define_optinfo {} {
#-------------------------------------------------------------------------------
  
  global opt

  set opt(fhead)    "out"
  set opt(cindex)    " 1  2  3  4  5  6  7  8  9 10 \
                        11 12 13 14 15 16 17 18 19 20"
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc read_optinfo {arglist} {
#-------------------------------------------------------------------------------

  global opt

  set opt(fhead) [parse_arguments $arglist \
      "-fhead"     "value" $opt(fhead)]
  set opt(cindex) [parse_arguments $arglist \
      "-cindex"     "value" $opt(cindex)]
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc show_optinfo {} {
#-------------------------------------------------------------------------------

  global opt

  puts "<< option info >>"
  puts "fhead      = $opt(fhead)"
  puts "cindex     = $opt(cindex)"
  puts ""

}

proc analyze {} {
  # in
  global str
  global traj
  global seltxt
  global opt
  global common

  global sel 

  set anatra_path $::env(ANATRA_PATH);list
  set fort       "${anatra_path}/f90/bin/lipid_order.x";list

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

  # Convert 
  #
  puts ""
  puts "--------------------"
  puts " Start analysis"
  puts "--------------------"
  puts ""

  set floinp   [format "%s.lo.inp" $opt(fhead) ]
  set floout   [format "%s.lo.out" $opt(fhead) ]

  for {set isel 0} {$isel < $nsel} {incr isel} {
    set fmolinfo($isel) [format "%s.lo.%i.molinfo" $opt(fhead) $isel]
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

  set ntraj   [llength $traj(tin)]

  set f [open $floinp "w"]
  puts $f " &input_param"
  puts $f "   flist_traj = \"$traj(flist_traj)\""

  if {$ntraj > 0} {
    puts $f "   ftraj ="
    for {set i 0} {$i < $ntraj} {incr i} {
      set t [lindex $traj(tin) $i]
      puts -nonewline $f "    \"$t\" "
    }
  }
  puts $f ""
  puts $f " /"
  puts $f " &output_param"
  puts $f "   fhead   = \"$opt(fhead)\""
  puts $f " /"
  
  puts $f " &trajopt_param"
  puts $f "   dt      = 1.0 ! Not used in the analysis"
  puts $f "   molinfo = "
  for {set i 0} {$i < $nsel} {incr i} {
    puts -nonewline $f "    \"$fmolinfo($i)\" "
  }
  puts $f ""
  puts $f " /"
  
  puts $f " &option_param"
  puts $f "   carbon_id = $opt(cindex)"
  puts $f " /"
  close $f

  #animate write dcd $fdcdtmp beg 0 end -1 waitfor all sel $sel($isel) $mol
  #

  if {!$common(prep_only)} {
    puts "Lipid order is calculated with ANATRA fortran program:"
    puts "$fort ..."
    puts "=== INPUT ==="
    set content [exec cat $floinp]
    puts $content
    puts "============="
    exec $fort $floinp >& $floout
    puts ""
    puts "=== OUTPUT ==="
    set content [exec cat $floout]
    puts $content
  }

  puts "=============="
  puts ">> Finished"

  exit
}
#-------------------------------------------------------------------------------
