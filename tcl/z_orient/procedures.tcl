#-------------------------------------------------------------------------------
proc print_title {} {
#-------------------------------------------------------------------------------

  puts "============================================================"
  puts ""
  puts "                  Z-orientation Analysis"
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
    puts "anatra z_orient                                                                             \\"
    puts "  -stype        <structure file type>                                                       \\"
    puts "  -sfile        <structure file name>                                                       \\"
    puts "  -tintype      <input trajectory file type>                                                \\"
    puts "  -tin          <input trajectory file name>                                                \\"
    puts "  -flist_traj   <trajectory file list (neccesary if tin is not specified)>                  \\"
    puts "  -fhead        <header of output file name>                                                \\"
    puts "  -judgeup      <how to judge up region (nmolup or coord or none)>                          \\"
    puts "                (default: none)                                                             \\"
    puts "  -nmolup       <number of molecules in upper leaflet>                                      \\"
    puts "  -zdef_type    <definition type of z-direction>                                            \\"
    puts "                (default: coord)                                                            \\"
#    puts "  -orient_type  <coordinate for distribution (costheta or theta)                            \\"
    puts "                (default: costheta)>                                                        \\"
    puts "  -dx           <costheta or theta grid interval>                                           \\"
    puts "  -dt           <time interval>                                                             \\"
    puts "  -sel0         <VMD selection>                                                             \\"
    puts "  -sel1         <VMD selection>                                                             \\"
    puts "  -sel2         <VMD selection>                                                             \\"
    puts "  -mode0        <analysis mode for sel0 (residue or whole or segment)>                      \\"
    puts "                (default: residue)>                                                         \\"
    puts "  -mode1        <analysis mode for sel1 (residue or whole or segment)>                      \\"
    puts "                (default: residue)>                                                         \\"
    puts "  -mode2        <analysis mode for sel2 (residue or whole or segment)>                      \\"
    puts "                (default: residue)>                                                         \\"
    puts "  -bottom_selid <selection id for bottom>                                                   \\"
    puts "  -head_selid   <selection id for arrowhead>                                                \\"
    puts "  -center_selid <selection id for system center>                                            \\"
    puts "  -vec0_bottom_selid <selection id for vec0 bottom (necessary if zdef_type = outerprod>     \\"
    puts "  -vec0_head_selid   <selection id for vec0 arrowhead (necessary if zdef_type = outerprod>  \\"
    puts "  -vec1_bottom_selid <selection id for vec1 bottom (necessary if zdef_type = outerprod>     \\"
    puts "  -vec1_head_selid   <selection id for vec1 arrowhead (necessary if zdef_type = outerprod>  \\"
    puts "  -prep_only    <whether analysis is performed or not (true or false)>                      \\" 
    puts "                (default: false)"
    puts ""
    puts "Usage:"
    puts "anatra z_orient                      \\"
    puts "  -stype              psf            \\"
    puts "  -sfile              complex.psf    \\"
    puts "  -tintype            dcd            \\"
    puts "  -tin                inp.dcd        \\"
    puts "  -fhead              run            \\"
    puts "  -judgeup            none           \\"
    puts "  -orient_type        costheta       \\"
    puts "  -zdef_type          coord          \\"
    puts "  -dx                 0.1            \\"
    puts "  -dt                 0.1            \\"
    puts "  -sel0               name P         \\"
    puts "  -sel1               name N         \\"
    puts "  -sel2               segid MEMB     \\"
    puts "  -mode0              residue        \\"
    puts "  -mode1              residue        \\"
    puts "  -mode2              whole          \\"
    puts "  -bottom_selid       0              \\"
    puts "  -head_selid         1              \\"
    puts "  -center_selid       2              \\"
    puts "  -vec0_bottom_selid  3              \\"
    puts "  -vec0_head_selid    4              \\"
    puts "  -vec1_bottom_selid  5              \\"
    puts "  -vec1_head_selid    6              \\"
    puts "  -prep_only       false"
    puts ""
    exit
  }

}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc define_optinfo {} {
#-------------------------------------------------------------------------------
  
  global opt

  set opt(fhead)             "run"
  set opt(judgeup)           "none"
  set opt(nmolup)            0 
  set opt(orient_type)       "costheta"
  set opt(zdef_type)         "coord"
  set opt(dx)                0.1 
  set opt(dt)                1 
  set opt(mode0)             "residue"
  set opt(mode1)             "residue"
  set opt(mode2)             "residue"
  set opt(bottom_selid)      0
  set opt(head_selid)        1 
  set opt(center_selid)      2
  set opt(vec0_bottom_selid) 3 
  set opt(vec0_head_selid)   4 
  set opt(vec1_bottom_selid) 5 
  set opt(vec1_head_selid)   6 

}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc read_optinfo {arglist} {
#-------------------------------------------------------------------------------

  global opt

  set opt(fhead)         [parse_arguments $arglist \
      "-fhead"           "value" $opt(fhead)]
  set opt(judgeup)       [parse_arguments $arglist \
      "-judgeup"         "value" $opt(judgeup)]
  set opt(nmolup)        [parse_arguments $arglist \
      "-nmolup"          "value" $opt(nmolup)]
  set opt(orient_type)   [parse_arguments $arglist \
      "-orient_type"     "value" $opt(orient_type)]
  set opt(zdef_type)     [parse_arguments $arglist \
      "-zdef_type"       "value" $opt(zdef_type)]
  set opt(dcost)         [parse_arguments $arglist \
      "-dx"              "value" $opt(dx)]
  set opt(dt)            [parse_arguments $arglist \
      "-dt"              "value" $opt(dt)]
  set opt(mode0)         [parse_arguments $arglist \
      "-mode0"           "value" $opt(mode0)]
  set opt(mode1)         [parse_arguments $arglist \
      "-mode1"           "value" $opt(mode1)]
  set opt(mode2)         [parse_arguments $arglist \
      "-mode2"           "value" $opt(mode2)]
  set opt(bottom_selid)  [parse_arguments $arglist \
      "-bottom_selid"    "value" $opt(bottom_selid)]
  set opt(head_selid)    [parse_arguments $arglist \
      "-head_selid"      "value" $opt(head_selid)]
  set opt(center_selid)  [parse_arguments $arglist \
      "-center_selid"    "value" $opt(center_selid)]
  set opt(vec0_bottom_selid)  [parse_arguments $arglist \
      "-vec0_bottom_selid" "value" $opt(vec0_bottom_selid)]
  set opt(vec0_head_selid)    [parse_arguments $arglist \
      "-vec0_head_selid" "value" $opt(vec0_head_selid)]
  set opt(vec1_bottom_selid)  [parse_arguments $arglist \
      "-vec1_bottom_selid" "value" $opt(vec1_bottom_selid)]
  set opt(vec1_head_selid)    [parse_arguments $arglist \
      "-vec1_head_selid" "value" $opt(vec1_head_selid)]
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc show_optinfo {} {
#-------------------------------------------------------------------------------

  global opt

  puts "<< option info >>"
  puts "fhead             = $opt(fhead)"
  puts "judgeup           = $opt(judgeup)"
  puts "nmolup            = $opt(nmolup)"
  puts "orient_type       = $opt(orient_type)"
  puts "zdef_type         = $opt(zdef_type)"
  puts "dx                = $opt(dx)"
  puts "dt                = $opt(dt)"
  puts "mode0             = $opt(mode0)"
  puts "mode1             = $opt(mode1)"
  puts "mode2             = $opt(mode2)"
  puts "bottom_selid      = $opt(bottom_selid)"
  puts "head_selid        = $opt(head_selid)"
  puts "center_selid      = $opt(center_selid)"
  puts "vec0_bottom_selid = $opt(vec0_bottom_selid)"
  puts "vec0_head_selid   = $opt(vec0_head_selid)"
  puts "vec1_bottom_selid = $opt(vec1_bottom_selid)"
  puts "vec1_head_selid   = $opt(vec1_head_selid)"
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
  set fort     "${anatra_path}/f90/bin/z_orient.x";list

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

  set nf [molinfo $mol get numframes] 
  puts ""

  puts ">> Start CoM calculation"
  puts ""

  set forinp   [format "%s.or.inp" $opt(fhead) ]
  set forout   [format "%s.or.out" $opt(fhead) ]

  for {set isel 0} {$isel < $nsel} {incr isel} {
    set fmolinfo($isel) [format "%s.or.%i.molinfo" $opt(fhead) $isel]
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

    #animate write dcd $fdcdtmp($isel) \
    #  beg 0 end -1 waitfor all sel $sel($isel) $mol 
  }

  set ntraj [llength $traj(tin)] 

  set sb   $opt(bottom_selid)
  set sh   $opt(head_selid)
  set sc   $opt(center_selid)
  set sv0b $opt(vec0_bottom_selid)
  set sv0h $opt(vec0_head_selid)
  set sv1b $opt(vec1_bottom_selid)
  set sv1h $opt(vec1_head_selid)

  set f [open $forinp "w"]
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
  puts $f "   dt      = $opt(dt)"

  set ztype [string tolower $opt(zdef_type)]
  if {$ztype == "coord"} {
    puts $f "   molinfo = \"$fmolinfo($sb)\" \"$fmolinfo($sh)\" \"$fmolinfo($sc)\""
  } else {
    puts $f "   molinfo = \"$fmolinfo($sb)\"   \"$fmolinfo($sh)\" \"$fmolinfo($sc)\""
    puts $f "             \"$fmolinfo($sv0b)\" \"$fmolinfo($sv0h)\""
    puts $f "             \"$fmolinfo($sv1b)\" \"$fmolinfo($sv1h)\""
  }

  puts $f " /"

  puts $f " &option_param"
  puts $f "   judgeup       = \"$opt(judgeup)\""
  puts $f "   nmolup        = $opt(nmolup)"
  puts $f "   orient_type   = \"$opt(orient_type)\""
  puts $f "   zdef_type     = \"$opt(zdef_type)\""
  puts $f "   dx            = $opt(dx)"
  puts $f "   mode          = \"$opt(mode$sb)\" \"$opt(mode$sh)\"  \"$opt(mode$sc)\""
  puts $f " /"

  close $f


  if {!$common(prep_only)} {
    puts "Z-orient is calculated with ANATRA fortran program:"
    puts "$fort ..."
    puts "=== INPUT ==="
    set content [exec cat $forinp]
    puts $content
    puts "============="
    exec $fort $forinp >& $forout
    puts ""
    puts "=== OUTPUT ==="
    set content [exec cat $forout]
    puts $content
  }

  puts "=============="
  puts ">> Finished"

  exit
}
#-------------------------------------------------------------------------------
