#-------------------------------------------------------------------------------
proc print_title {} {
#-------------------------------------------------------------------------------

  puts "============================================================"
  puts ""
  puts "                    Z-Profile Analysis"
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
    puts "anatra z_profile                                                           \\"
    puts "  -stype        <structure file type>                                      \\"
    puts "  -sfile        <structure file name>                                      \\"
    puts "  -tintype      <input trajectory file type>                               \\"
    puts "  -tin          <input trajectory file name>                               \\"
    puts "  -flist_traj   <trajectory file list (neccesary if tin is not specified)> \\"
    puts "  -fhead        <header of output file name>                               \\"
    puts "  -dz           <z-grid interval>                                          \\"
    puts "  -sel0         <VMD selection>                                            \\"
    puts "  -sel1         <VMD selection>                                            \\"
    puts "  -mode0        <analysis mode (residue or whole or atom)>                 \\" 
    puts "                (default: residue)                                         \\"
    puts "  -mode1        <analysis mode (residue or whole or atom)>                 \\" 
    puts "                (default: residue)                                         \\"
    puts "  -target_selid <selection id for target species>                          \\"
    puts "                (default: 0)                                               \\"
    puts "  -center_selid <selection id for system center>                           \\"
    puts "                (default: 1)                                               \\"
    puts "  -denstype     <density type (number or electron)>                        \\"
    puts "                (default: number)                                          \\"
    puts "  -symmetrize   <whether symmetrized w.r.t. z=0 or not (true or false)>    \\"
    puts "                (default: false)                                           \\"
    puts "  -prep_only    <whether analysis is performed or not (true or false)>     \\"
    puts "                (default: false)>"
    puts ""
    puts "Usage:"
    puts "anatra z_profile             \\"
    puts "  -stype        psf          \\"
    puts "  -sfile        complex.psf  \\"
    puts "  -tintype      dcd          \\"
    puts "  -tin          run.dcd      \\"
    puts "  -fhead        run          \\"
    puts "  -dz           0.5          \\"
    puts "  -sel0         water        \\"
    puts "  -sel1         segid MEMB   \\"
    puts "  -mode0        residue      \\"
    puts "  -mode1        whole        \\"
    puts "  -target_selid 0            \\"
    puts "  -center_selid 1            \\"
    puts "  -denstype     number       \\"
    puts "  -symmetrize   false        \\"
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
  set opt(mode0)         "residue"
  set opt(mode1)         "residue"
  set opt(target_selid)  0 
  set opt(center_selid)  0
  set opt(denstype)      "number" 
  set opt(symmetrize)    "true" 

}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc read_optinfo {arglist} {
#-------------------------------------------------------------------------------

  global opt

  set opt(fhead)        [parse_arguments $arglist \
      "-fhead"          "value" $opt(fhead)]
  set opt(dz)           [parse_arguments $arglist \
      "-dz"             "value" $opt(dz)]
  set opt(mode0)        [parse_arguments $arglist \
      "-mode0"          "value" $opt(mode0)]
  set opt(mode1)        [parse_arguments $arglist \
      "-mode1"          "value" $opt(mode1)]
  set opt(target_selid) [parse_arguments $arglist \
      "-target_selid"   "value" $opt(target_selid)]
  set opt(center_selid) [parse_arguments $arglist \
      "-center_selid"   "value" $opt(center_selid)]
  set opt(denstype)     [parse_arguments $arglist \
      "-denstype"       "value" $opt(denstype)]
  set opt(symmetrize)   [parse_arguments $arglist \
      "-symmetrize"     "value" $opt(symmetrize)]
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc show_optinfo {} {
#-------------------------------------------------------------------------------

  global opt

  puts "<< option info >>"
  puts "fhead        = $opt(fhead)"
  puts "dz           = $opt(dz)"
  puts "mode0        = $opt(mode0)"
  puts "mode1        = $opt(mode1)"
  puts "target_selid = $opt(target_selid)"
  puts "center_selid = $opt(center_selid)"
  puts "denstype     = $opt(denstype)"
  puts "symmetrize   = $opt(symmetrize)"
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
  set fort       "${anatra_path}/f90/bin/z_profile.x";list

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

  puts ">> Start Z-profile calculation"
  puts ""

  set fzpinp   [format "%s.zp.inp" $opt(fhead) ]
  set fzpout   [format "%s.zp.out" $opt(fhead) ]

  for {set isel 0} {$isel < $nsel} {incr isel} {
    set fmolinfo($isel) [format "%s.zp.%i.molinfo" $opt(fhead) $isel]
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

  set f [open $fzpinp "w"]
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
  puts $f "   molinfo = \"$fmolinfo($st)\" \"$fmolinfo($sc)\""
  puts $f " /"

  puts $f " &option_param"
  puts $f "   mode       = \"$opt(mode$st)\" \"$opt(mode$sc)\""
  puts $f "   dz         = $opt(dz)"
  puts $f "   denstype   = \"$opt(denstype)\""
  puts $f "   symmetrize = .$opt(symmetrize)."
  puts $f " /"

  close $f

  #animate write dcd $fdcdtmp beg 0 end -1 waitfor all sel $sel(0) $mol

  if {!$common(prep_only)} {
    puts "Z-profile is calculated with ANATRA fortran program:"
    puts "$fort ..."
    puts "=== INPUT ==="
    set content [exec cat $fzpinp]
    puts $content
    puts "============="
    exec $fort $fzpinp > $fzpout
    puts ""
    puts "=== OUTPUT ==="
    set content [exec cat $fzpout]
    puts $content
  }
  puts "=============="
  puts ">> Finished"

  exit
}
#-------------------------------------------------------------------------------
