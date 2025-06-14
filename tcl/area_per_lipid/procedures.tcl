#-------------------------------------------------------------------------------
proc print_title {} {
#-------------------------------------------------------------------------------

  puts "============================================================"
  puts ""
  puts "                  Area Per Lipid Analysis"
  puts ""
  puts "============================================================"

}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc show_usage {arglist} {
#-------------------------------------------------------------------------------

  set help false
  set help [parse_arguments $arglist "-h" "flag" $help]

  if {$help} {
    puts "Usage:"
    puts "anatra area_per_lipid                        \\"
    puts "  -stype      <structure file type>          \\"
    puts "  -sfile      <structure file name>          \\"
    puts "  -tintype    <input trajectory file type>   \\"
    puts "  -tin        <input trajectory file name>   \\"
    puts "  -fhead      <header of output file name>   \\"
    puts "  -dt         <time interval>                \\"
    puts "  -da         <delta area per lipid>         \\"
    puts "  -sel0       <VMD selection> (X=0,1,2...)"
    puts ""
    puts "Usage:"
    puts "anatra area_per_lipid           \\"
    puts "  -stype      psf               \\"
    puts "  -sfile      complex.psf       \\"
    puts "  -tintype    dcd               \\"
    puts "  -tin        inp.dcd           \\"
    puts "  -fhead      out               \\"
    puts "  -dt         0.01              \\"
    puts "  -da         0.5               \\"
    puts "  -sel0       segid MEMB"
    puts ""
    exit
  }

}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc define_optinfo {} {
#-------------------------------------------------------------------------------
  
  global opt

  set opt(fhead)  "run" 
  set opt(dt)     0.01
  set opt(da)     0.50 
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc read_optinfo {arglist} {
#-------------------------------------------------------------------------------

  global opt

  set opt(fhead)      [parse_arguments $arglist \
      "-fhead"       "value" $opt(fhead)]
  set opt(dt)         [parse_arguments $arglist \
      "-dt"          "value" $opt(dt)]
  set opt(da)         [parse_arguments $arglist \
      "-da"          "value" $opt(da)]
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc show_optinfo {} {
#-------------------------------------------------------------------------------

  global opt

  puts "<< option info >>"
  puts "fhead      = $opt(fhead)"
  puts "dt         = $opt(dt)"
  puts "da         = $opt(da)"
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


  set anatra_path $::env(ANATRA_PATH);list
  set fort       "${anatra_path}/f90/bin/histogram.x";list

  # read trajectory
  #
  puts ""
  puts "--------------------"
  puts " Read trajectory"
  puts "--------------------"
  puts ""

  set mol 0;
  read_traj $mol $str(stype) $str(sfile) $traj(tintype) $traj(tin) $traj(stride)
  set nf   [molinfo $mol get numframes]
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

  set fapinp   [format "%s.ap.inp" $opt(fhead) ]
  set fapout   [format "%s.ap.out" $opt(fhead) ]

  set nres [llength [lsort -unique [$sel(0) get resid]]]
  package require pbctools
  set box [pbc get -all]
  for {set istep 0} {$istep < $nf} {incr istep} {
    set bnow [lindex $box $istep]
    set apl($istep) [expr [lindex $bnow 0] * [lindex $bnow 1] * 2 / $nres ] 
  }

  set aplfile [format "%s.apl" $opt(fhead) ]
  set f [open $aplfile "w"]
  for {set istep 0} {$istep < $nf} {incr istep} {
    set t [expr $istep * $opt(dt)]
    puts $f [format "%20.10f  %20.10f" $t $apl($istep)]
  }
  close $f

  set f [open $fapinp "w"]
  puts $f " &input_param"
  puts $f "   fcv = \"$aplfile\""
  puts $f " /"
  puts $f " &output_param"
  puts $f "   fhead          = \"$opt(fhead)\""
  puts $f "   file_extension = \"apldistr\""
  puts $f " /"
  puts $f " &option_param"
  puts $f "   dx       = $opt(da)"
  puts $f "   xsta     = 30.0"
  puts $f "   nx       = 1000"
  puts $f "   ncol     = 1"
  puts $f " /"

  close $f

  puts "APL is calculated with ANATRA fortran program:"
  puts "$fort ..."
  puts "=== INPUT ==="
  set  content [exec cat $fapinp]
  puts $content
  puts "============="
  exec $fort $fapinp >& $fapout
  puts ""
  puts "=== OUTPUT ==="
  set  content [exec cat $fapout]
  puts $content
  puts "=============="
  puts ">> Finished"

  exit
}
#-------------------------------------------------------------------------------
