#-------------------------------------------------------------------------------
proc print_title {} {
#-------------------------------------------------------------------------------

  puts "============================================================"
  puts ""
  puts "                  Gaussian QM Prepper"
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
    puts "anatra z_profile                                                                          \\"
    puts "  -stype        <structure file type>                                                     \\"
    puts "  -sfile        <structure file name>                                                     \\"
    puts "  -tintype      <input trajectory file type>                                              \\"
    puts "  -tin          <input trajectory file name>                                              \\"
    puts "  -flist_traj   <trajectory file list (neccesary if tin is not specified)>                \\"
    puts "  -fhead        <header of output file name>                                              \\"
    puts "  -sel0         <VMD selection>                                                           \\"
    puts "  -sel1         <VMD selection>                                                           \\"
    puts "  -prep_only    <whether analysis is performed or not (true or false)>                    \\"
    puts "                (default: false)>"
    puts ""
    puts "Usage:"
    puts "anatra gauqm_prepper         \\"
    puts "  -stype        psf          \\"
    puts "  -sfile        complex.psf  \\"
    puts "  -tintype      dcd          \\"
    puts "  -tin          run.dcd      \\"
    puts "  -fhead        run          \\"
    puts "  -sel0         water        \\"
    puts "  -sel1         not water    \\"
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
  set opt(dz)            0.5
  set opt(mode0)         "whole"
  set opt(mode1)         "whole"

}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc read_optinfo {arglist} {
#-------------------------------------------------------------------------------

  global opt

  set opt(fhead)        [parse_arguments $arglist \
      "-fhead"          "value" $opt(fhead)]
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc show_optinfo {} {
#-------------------------------------------------------------------------------

  global opt

  puts "<< option info >>"
  puts "fhead        = $opt(fhead)"
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
  set fort       "${anatra_path}/f90/bin/gauqm_prepper.x";list

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

  puts ">> Start Gaussian Input preparation"
  puts ""

  set fgauinp   [format "%s.gau.inp" $opt(fhead) ]
  set fgauout   [format "%s.gau.out" $opt(fhead) ]

  for {set isel 0} {$isel < $nsel} {incr isel} {
    set fmolinfo($isel) [format "%s.gau.%i.molinfo" $opt(fhead) $isel]
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


  set f [open $fgauinp "w"]
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
  puts $f "   fhead = \"$opt(fhead)\""
  puts $f " /"

  puts $f " &trajopt_param"
  puts $f "   molinfo = \"$fmolinfo(0)\" \"$fmolinfo(1)\""
  puts $f " /"

  puts $f " &option_param"
  puts $f "   mode       = \"$opt(mode0)\" \"$opt(mode1)\""
  puts $f " /"

  close $f

  #animate write dcd $fdcdtmp beg 0 end -1 waitfor all sel $sel(0) $mol

  if {!$common(prep_only)} {
    puts "Gaussian QM prepper is done with ANATRA fortran program:"
    puts "$fort ..."
    puts "=== INPUT ==="
    set content [exec cat $fgauinp]
    puts $content
    puts "============="
    exec $fort $fgauinp > $fgauout
    puts ""
    puts "=== OUTPUT ==="
    set content [exec cat $fgauout]
    puts $content
  }
  puts "=============="
  puts ">> Finished"

  exit
}
#-------------------------------------------------------------------------------
