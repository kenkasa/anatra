#-------------------------------------------------------------------------------
proc print_title {} {
#-------------------------------------------------------------------------------

  puts "============================================================"
  puts ""
  puts "                  CoM Coordinate Analysis"
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
    puts "anatra center_of_mass                                                      \\"
    puts "  -stype      <structure file type>                                        \\"
    puts "  -sfile      <structure file name>                                        \\"
    puts "  -tintype    <input trajectory file type>                                 \\"
    puts "  -tin        <input trajectory file name>                                 \\"
    puts "  -flist_traj <trajectory file list (neccesary if tin is not specified)>   \\"
    puts "  -fhead      <header of output file name>                                 \\"
    puts "  -sel0       <VMD selection> (X=0,1,2...)                                 \\"
    puts "  -mode       <analysis mode (residue or whole or atom)>                   \\"
    puts "              (default: residue)                                           \\"
    puts "  -out_com    <whether com file is generated (true or false)>              \\"
    puts "              (default:false)>                                             \\"
    puts "  -out_msd   <whether msd file is generated (true or false)>               \\"
    puts "              (default:false)>                                             \\"
    puts "  -unwrap     <whether unwrapping trajecotry is performed (true or false)> \\"
    puts "              (default:false)                                              \\"
    puts "  -msddim     <dimension used for MSD analysis (2 or 3)>                   \\"
    puts "              (default:3)                                                  \\"
    puts "  -dt         <time interval>                                              \\"
    puts "  -t_sparse   <output time interval>                                       \\"
    puts "  -t_range    <analysis time range>                                        \\"
    puts "  -prep_only  <where analysis is performed or not (true or false)          \\"
    puts "              (default: false)"
    puts ""
    puts "Usage:"
    puts "anatra center_of_mass   \\"
    puts "  -stype     parm7      \\"
    puts "  -sfile     str.prmtop \\"
    puts "  -tintype   dcd        \\"
    puts "  -tin       inp.dcd    \\"
    puts "  -fhead     out        \\"
    puts "  -out_com   true       \\"
    puts "  -out_msd   true       \\"
    puts "  -sel0      not water  \\"
    puts "  -mode      residue    \\"
    puts "  -unwrap    false      \\"
    puts "  -msddim    3          \\"
    puts "  -dt        0.1        \\"
    puts "  -t_sparse  1.0        \\"
    puts "  -t_range   10         \\"
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
  set opt(out_com)    false
  set opt(out_msd)    false
  set opt(onlyz)      false
  set opt(unwrap)     false
  set opt(msddim)     3
  set opt(dt)         0.1 
  set opt(t_sparse)   -1.0
  set opt(t_range)    -1.0
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
  set opt(out_com)    [parse_arguments $arglist \
      "-out_com"    "value" $opt(out_com)]
  set opt(out_msd)    [parse_arguments $arglist \
      "-out_msd"    "value" $opt(out_msd)]
  set opt(unwrap)     [parse_arguments $arglist \
      "-unwrap"     "value" $opt(unwrap)]
  set opt(onlyz)      [parse_arguments $arglist \
      "-onlyz"      "value" $opt(onlyz)]
  set opt(msddim)     [parse_arguments $arglist \
      "-msddim"     "value" $opt(msddim)]
  set opt(dt)         [parse_arguments $arglist \
      "-dt"         "value" $opt(dt)]
  set opt(t_sparse)   [parse_arguments $arglist \
      "-t_sparse"   "value" $opt(t_sparse)]
  set opt(t_range)    [parse_arguments $arglist \
      "-t_range"    "value" $opt(t_range)]
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc show_optinfo {} {
#-------------------------------------------------------------------------------

  global opt

  puts "<< option info >>"
  puts "fhead      = $opt(fhead)"
  puts "fts        = $opt(fts)"
  puts "out_com    = $opt(out_com)"
  puts "out_msd    = $opt(out_msd)"
  puts "onlyz      = $opt(onlyz)"
  puts "msddim     = $opt(msddim)"
  puts "dt         = $opt(dt)"
  puts "t_sparse   = $opt(t_sparse)"
  puts "t_range    = $opt(t_range)"
  puts "mode       = $opt(mode)"
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
  set fort       "${anatra_path}/f90/bin/center_of_mass.x";list

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

  puts ">> Start CoM calculation"
  puts ""
  
  set fcminp   [format "%s.cm.inp" $opt(fhead) ]
  set fcmout   [format "%s.cm.out" $opt(fhead) ]
  set fmolinfo [format "%s.cm.molinfo" $opt(fhead) ]

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

  set f [open $fcminp "w"]
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
  puts $f "   mode     = \"$opt(mode)\""
  puts $f "   out_com  = .$opt(out_com)."
  puts $f "   out_msd  = .$opt(out_msd)."
  puts $f "   onlyz    = .$opt(onlyz)."
  puts $f "   unwrap   = .$opt(unwrap)."
  puts $f "   msddim   = $opt(msddim)"
  puts $f "   dt       = $opt(dt)"
  puts $f "   t_sparse = $opt(t_sparse)"
  puts $f "   t_range  = $opt(t_range)"
  puts $f " /"

  close $f

  if {!$common(prep_only)} {
    puts "CoM is calculated with ANATRA fortran program:"
    puts "$fort ..."
    puts "=== INPUT ==="
    set content [exec cat $fcminp]
    puts $content
    puts "============="
    exec $fort $fcminp >& $fcmout
    puts ""
    puts "=== OUTPUT ==="
    set content [exec cat $fcmout]
    puts $content
  }

  puts "=============="
  puts ">> Finished"

  exit

}
#-------------------------------------------------------------------------------
