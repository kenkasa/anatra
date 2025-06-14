#-------------------------------------------------------------------------------
proc print_title {} {
#-------------------------------------------------------------------------------
  puts "============================================================"
  puts ""
  puts "         Root Mean Square Deviation (RMSD) Analysis"
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
    puts "anatra rmsd                                               \\"
    puts "  -stype        <structure file type>                     \\"
    puts "  -sfile        <structure file name>                     \\"
    puts "  -tintype      <input trajectory file type>              \\"
    puts "  -tin          <input trajectory file name>              \\"
    puts "  -fout         <output file name>                        \\"
    puts "  -selX         <X-th VMD selection> (X=0,1,2...)         \\"
    puts "  -fit          <fit is performed or not (true or false)> \\"
    puts "                (default: false)                          \\"
    puts "  -refpdb       <reference pdb file name>                 \\"
    puts "  -fitselid     <selection id for fitting>                \\"
    puts "  -refselid     <selection id for reference>              \\"
    puts "  -rmsdselid    <selection id for rmsd>                   \\"
    puts "  -rmsdrefselid <selection id for rmsd reference>"
    puts ""
    puts "Usage:"
    puts "anatra rmsd                                   \\"
    puts "  -stype        parm7                         \\"
    puts "  -sfile        str.prmtop                    \\"
    puts "  -tintype      dcd                           \\"
    puts "  -tin          inp.dcd                       \\"
    puts "  -fout         out.rmsd                      \\"
    puts "  -sel0         resid 1 to 275 and name CA    \\"
    puts "  -sel1         resid 1 to 275 and name CA    \\"
    puts "  -sel2         resid 1 to 275 and name CA    \\"
    puts "  -sel3         resid 1 to 275 and name CA    \\"
    puts "  -fit          true                          \\"
    puts "  -refpdb       ref.pdb                       \\"
    puts "  -fitselid     0                             \\"
    puts "  -refselid     1                             \\"
    puts "  -rmsdselid    2                             \\"
    puts "  -rmsdrefselid 3 "
    puts ""
    exit
  }

}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc define_optinfo {} {
#-------------------------------------------------------------------------------
  
  global opt

  set opt(fout)          "out.rmsd"
  set opt(fit)           false
  set opt(wrap)          false
  set opt(centering)     false

  set opt(wrapcenter)    origin
  set opt(wrapcomp)      fragment
  set opt(fitselid)      0 
  set opt(refselid)      0 
  set opt(rmsdselid)     0 
  set opt(rmsdrefselid)  0 
  set opt(refpdb)        ""
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc read_optinfo {arglist} {
#-------------------------------------------------------------------------------

  global opt

  set opt(fout)          [parse_arguments $arglist \
      "-fout"         "value" $opt(fout)]
  set opt(fit)           [parse_arguments $arglist \
      "-fit"          "value" $opt(fit)]
  set opt(fitselid)      [parse_arguments $arglist \
      "-fitselid"     "value" $opt(fitselid)]
  set opt(refselid)      [parse_arguments $arglist \
      "-refselid"     "value" $opt(refselid)]
  set opt(rmsdselid)     [parse_arguments $arglist \
      "-rmsdselid"    "value" $opt(rmsdselid)]
  set opt(rmsdrefselid)  [parse_arguments $arglist \
      "-rmsdrefselid" "value" $opt(rmsdrefselid)]
  set opt(refpdb)        [parse_arguments $arglist \
      "-refpdb"       "value" $opt(refpdb)]

  # Combination error check
  #

  if {$opt(refpdb) == ""} {
    puts "ERROR: refpdb should be specified."
    exit
  }

}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
proc fitting {molid seltxt_fit seltxt_out refid seltxt_ref} {
#-------------------------------------------------------------------------------

  set refsel [atomselect $refid $seltxt_ref]
  set fitsel [atomselect $molid $seltxt_fit]
  set outsel [atomselect $molid $seltxt_out]

  set nf [molinfo $molid get numframes]
   
  for {set i 0} {$i < $nf} {incr i} {
    if {[expr $i % 100] == 0 || [expr $i + 1] == $nf} {
      puts [format "%10d / %10d" $i [expr $nf - 1]] 
    } 
    $fitsel frame $i
    $outsel frame $i
    $outsel move [measure fit $fitsel $refsel weight mass]
  } 

}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc show_optinfo {} {
#-------------------------------------------------------------------------------

  global opt

  puts "<< option info >>"
  puts "fout         = $opt(fout)"
  puts "fit          = $opt(fit)"
  puts "fitselid     = $opt(fitselid)"
  puts "refselid     = $opt(refselid)"
  puts "rmsdselid    = $opt(rmsdselid)"
  puts "rmsdrefselid = $opt(rmsdrefselid)"
  puts "refpdb       = $opt(refpdb)"
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

  # setup reference 
  #
  puts ""
  puts "--------------------"
  puts " Setup reference"
  puts "--------------------"
  puts ""
  set ref        [mol load pdb "$opt(refpdb)"]
  set refsel     [atomselect $ref "$seltxt($opt(refselid))"]
  set rmsdrefsel [atomselect $ref "$seltxt($opt(rmsdrefselid))"]
  set comref     [measure center $refsel weight mass]

  set fitsel     [atomselect $mol "$seltxt($opt(fitselid))"]

  # Convert 
  #
  puts ""
  puts "--------------------"
  puts " Start analysis"
  puts "--------------------"
  puts ""

  set nf [molinfo $mol get numframes] 

  if {$opt(fit)} {
    puts "o Start fitting"
    fitting $mol \
            "$seltxt($opt(fitselid))" \
            "all"                         \
	          $ref                          \
	          "$seltxt($opt(refselid))"
    puts "> Finish fitting" 
     
  }

  for {set i 0} {$i < $nf} {incr i} {
    $sel($opt(rmsdselid)) frame $i
    set rmsval($i) [measure rmsd                 \
                    $sel($opt(rmsdselid))    \
                    $rmsdrefsel                  \
                    weight mass]
  }

  set fid [open "$opt(fout)" w]
  for {set i 0} {$i < $nf} {incr i} {
    puts $fid [format "%10d  %15.7f" $i $rmsval($i)]
  }
  close $fid

  puts ""
  puts ">> Finish all Analysis"
  puts ""

}
#-------------------------------------------------------------------------------
