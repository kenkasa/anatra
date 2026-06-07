#-------------------------------------------------------------------------------
proc print_title {} {
#-------------------------------------------------------------------------------

  puts "============================================================"
  puts ""
  puts "                 Cluster-size Analysis"
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
    puts "anatra cluster_size                                                        \\"
    puts "  -stype      <structure file type>                                        \\"
    puts "  -sfile      <structure file name>                                        \\"
    puts "  -tintype    <input trajectory file type>                                 \\"
    puts "  -tin        <input trajectory file name>                                 \\"
    puts "  -flist_traj <trajectory file list (neccesary if tin is not specified)>   \\"
    puts "  -fhead      <header of output file name>                                 \\"
    puts "  -sel0       <VMD selection> (X=0,1,2...)                                 \\"
    puts "  -mode       <analysis mode (residue or whole or atom)>                   \\"
    puts "              (default: residue)                                           \\"
    puts "  -rcut       <cutoff distance for judging (Angstrom)>                     \\"
    puts "  -prep_only  <where analysis is performed or not (true or false)          \\"
    puts "              (default: false)"
    puts ""
    puts "Usage:"
    puts "anatra cluster_size     \\"
    puts "  -stype     parm7      \\"
    puts "  -sfile     str.prmtop \\"
    puts "  -tintype   dcd        \\"
    puts "  -tin       inp.dcd    \\"
    puts "  -fhead     out        \\"
    puts "  -sel0      not water  \\"
    puts "  -mode      residue    \\"
    puts "  -rcut      3.5        \\"
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

  set opt(fhead)      "out"
  set opt(nparallel)  1
  set opt(rcut)       3.5 
  set opt(dt)         1.0 
  set opt(mode)       "residue"

}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc read_optinfo {arglist} {
#-------------------------------------------------------------------------------

  global opt

  set opt(fhead)      [parse_arguments $arglist \
      "-fhead"      "value" $opt(fhead)]
  set opt(mode)       [parse_arguments $arglist \
      "-mode"       "value" $opt(mode)]
  set opt(nparallel)  [parse_arguments $arglist \
      "-nparallel"  "value" $opt(nparallel)]
  set opt(rcut)       [parse_arguments $arglist \
      "-rcut"       "value" $opt(rcut)]
  set opt(dt)         [parse_arguments $arglist \
      "-dt"         "value" $opt(dt)]
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc show_optinfo {} {
#-------------------------------------------------------------------------------

  global opt

  puts "<< option info >>"
  puts "fhead      = $opt(fhead)"
  puts "mode       = $opt(mode)"
  puts "nparallel  = $opt(nparallel)"
  puts "rcut       = $opt(rcut)"
  #puts "dt         = $opt(dt)"
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
  set fort       "${anatra_path}/f90/bin/cluster_size.x";list

  # read trajectory
  #
  puts ""
  puts "--------------------"
  puts " Read trajectory"
  puts "--------------------"
  puts ""

  set mol [mol load $str(stype) "$str(sfile)"]

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
  puts ""

  puts ">> Start Cluster-size calculation"
  puts ""
  
  set fcsinp   [format "%s.cs.inp" $opt(fhead) ]
  set fcsout   [format "%s.cs.out" $opt(fhead) ]
  set fmolinfo [format "%s.cs.molinfo" $opt(fhead) ]

  set f [open $fmolinfo "w"]
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

  set ntraj [llength $traj(tin)] 

  set f [open $fcsinp "w"]
  puts $f " &input_param"
  puts $f  "  flist_traj = \"$traj(flist_traj)\""

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
  puts $f "   dt      = $opt(dt)"
  puts $f "   molinfo = \"$fmolinfo\""
  puts $f " /"

  puts $f " &option_param"
  puts $f "   mode      = \"$opt(mode)\""
  puts $f "   dt        = $opt(dt)"
  puts $f "   rcut      = $opt(rcut)"
  puts $f "   nparallel = $opt(nparallel)"
  puts $f " /"

  close $f

  if {!$common(prep_only)} {
    puts "Cluster-size is calculated with ANATRA fortran program:"
    puts "$fort ..."
    puts "=== INPUT ==="
    set content [exec cat $fcsinp]
    puts $content
    puts "============="
    exec $fort $fcsinp >& $fcsout
    puts ""
    puts "=== OUTPUT ==="
    set content [exec cat $fcsout]
    puts $content
  }

  puts "=============="
  puts ">> Finished"

  exit

}
#-------------------------------------------------------------------------------
