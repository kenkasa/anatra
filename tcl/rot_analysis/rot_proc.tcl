proc print_title {} {
  puts "============================================================"
  puts ""
  puts "              Rotational Correlation Analysis"
  puts ""
  puts "============================================================"
}

proc show_rot_usage {arglist} {

  set help false
  set help [parse_arguments $arglist \
      "-h" "flag" $help]

  if {$help} {
    puts "Usage:"
    puts "anatra rot                                     \\"
    puts "  -stype    <structure file type>              \\"
    puts "  -sfile    <structure file name>              \\"
    puts "  -tintype  <input trajectory file type>       \\"
    puts "  -tin      <input trajectory file name>       \\"
    puts "  -fhead    <header of output file name>       \\"
    puts "  -dt       <time interval>                    \\"
    puts "  -tcfrange <number of steps analyzed for TCF> \\"
    puts "  -sel0     <VMD selection>                    \\"
    puts "  -sel1     <VMD selection>                    \\"
    puts "  -mode     <analysis mode (residue or whole)> \\"
    puts "            (default: residue)"
    puts ""
    puts "Usage:"
    puts "anatra rot              \\"
    puts "  -stype    parm7       \\"
    puts "  -sfile    str.prmtop  \\"
    puts "  -tintype  dcd         \\"
    puts "  -tin      inp.dcd     \\"
    puts "  -fhead    out.rot     \\"
    puts "  -dt       0.1         \\"
    puts "  -tcfrange 1000        \\"
    puts "  -sel0     name P      \\"
    puts "  -sel1     name N      \\"
    puts "  -mode     residue"
    puts ""
    exit
  }

}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      define_rotoptinfo
#! @brief         define ROT option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc define_rotoptinfo {} {
  
  global rotopt

  set rotopt(fhead)      "out" 
  set rotopt(dt)         0.1
  set rotopt(tcfrange)   1000
  set rotopt(mode)       "residue" 
}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      read_rotoptinfo
#! @brief         read Dipole option paramerters 
#! @authors       KK
#! @param[in]  arglist : argument list
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc read_rotoptinfo {arglist} {

  global rotopt

  set rotopt(fhead)       [parse_arguments $arglist \
      "-fhead"        "value" $rotopt(fhead)]
  set rotopt(mode)        [parse_arguments $arglist \
      "-mode"         "value" $rotopt(mode)]
  set rotopt(dt)          [parse_arguments $arglist \
      "-dt"           "value" $rotopt(dt)]
  set rotopt(tcfrange)    [parse_arguments $arglist \
      "-tcfrange"     "value" $rotopt(tcfrange)]
}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      show_rotoptinfo
#! @brief         show ROT option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc show_rotoptinfo {} {

  global rotopt

  puts "<< option info >>"
  puts "fhead      = $rotopt(fhead)"
  puts "dt         = $rotopt(dt)"
  puts "tcfrange   = $rotopt(tcfrange)"
  puts "mode       = $rotopt(mode)"
  puts ""

}

proc rot_analysis {} {
  # in
  global str
  global traj
  global seltxt
  global rotopt

  global sel 


  set anatra_path $::env(ANATRA_PATH);list
  set rotfort    "${anatra_path}/f90/bin/rot_analysis.x";list

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

  set nf [molinfo $mol get numframes] 
  puts ""

  puts ">> Start CoM calculation"
  puts ""

  set rand [expr int((100000*rand()))]
  set frotinp   [format "rot%06d.inp"     $rand]
  set frotout   [format "rot%06d.out"     $rand]

  for {set isel 0} {$isel < $nsel} {incr isel} {
    set fdcdtmp($isel)  [format "rot%06d_%i.dcd"     $rand $isel]
    set fmolinfo($isel) [format "rot%06d_%i.molinfo" $rand $isel]
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

    animate write dcd $fdcdtmp($isel) \
      beg 0 end -1 waitfor all sel $sel($isel) $mol 
  }


  set f [open $frotinp "w"]
  puts $f " &input_param"
  puts $f "   fdcd = \"$fdcdtmp(0)\" \"$fdcdtmp(1)\""
  puts $f " /"
  puts $f " &output_param"
  puts $f "   fhead = \"$rotopt(fhead)\""
  puts $f " /"

  puts $f " &trajopt_param"
  puts $f "   dt      = $rotopt(dt)"
  puts $f "   molinfo = \"$fmolinfo(0)\" \"$fmolinfo(1)\""
  puts $f " /"

  puts $f " &option_param"
  puts $f "   tcfrange = $rotopt(tcfrange)"
  puts $f "   mode     = \"$rotopt(mode)\""
  puts $f " /"

  close $f


  puts "Rotation is calculated with ANATRA fortran program:"
  puts "$rotfort ..."
  puts "=== INPUT ==="
  set content [exec cat $frotinp]
  puts $content
  puts "============="
  exec $rotfort $frotinp >& $frotout
  puts ""
  puts "=== OUTPUT ==="
  set content [exec cat $frotout]
  puts $content
  puts "=============="
  puts ">> Finished"

  exec rm -f $fmolinfo(0) $fmolinfo(1) \
             $frotinp $frotout $fdcdtmp(0) $fdcdtmp(1) 
  exit
}
