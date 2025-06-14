#-------------------------------------------------------------------------------
proc print_title {} {
#-------------------------------------------------------------------------------

  puts "============================================================"
  puts ""
  puts "               Spatial-Distribution Analysis"
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
    puts "anatra spatial_distr                                                                           \\"
    puts "  -stype             <structure file type>                                                     \\"
    puts "  -sfile             <structure file name>                                                     \\"
    puts "  -tintype           <input trajectory file type>                                              \\"
    puts "  -tin               <input trajectory file name>                                              \\"
    puts "  -flist_traj        <trajectory file list (neccesary if tin is not specified)>                \\"
    puts "  -flist_cv          <input list file that contains list of cv files>                          \\"
    puts "                     (necessary if use_conditional is true)                                    \\"
    puts "  -flist_weight      <input list file that contains list of weight files>                      \\"
    puts "                     (necessary if use_weight is true)                                         \\"
    puts "  -fhead             <header of output file name>                                              \\"
    puts "  -sel0              <VMD selection> (X=0,1,2...)                                              \\"
    puts "  -mode              <analysis mode (residue or whole or atom)                                 \\"
    puts "                     (default: residue)>                                                       \\"
    puts "  -ng3               <number of grids for x, y, z axes>                                        \\"
    puts "  -del               <grid spacing for x, y, z axes>                                           \\"
    puts "  -origin            <origin of 3d-grids>                                                      \\"
    puts "  -use_spline        <whether spline is performed or not (true or false)>                      \\"
    puts "                     (default: false)>                                                         \\"
    puts "  -spline_resolution <spline resolution (integer)>                                             \\"
    puts "                     (default: 4)                                                              \\"
    puts "  -use_weight        <whether weights of each configuration is treated or not>                 \\"
    puts "                     (true or false) (default: false)                                          \\"
    puts "  -use_conditional   <whether conditional sampling is used or not>                             \\"
    puts "                     (true or false) (default: false)                                          \\"
    puts "  -ndim              <dimensions of reaction coords (neccesary if use_conditional is true)>    \\"
    puts "  -react_range       <range of sampled reaction coords (neccesary if use_conditional is true)> \\"
    puts "  -fit               <fit is performed or not (true or false)>                                 \\"
    puts "  -refpdb            <reference pdb file name>                                                 \\"
    puts "                     (necessary if fit = true)                                                 \\"
    puts "  -fitselid          <selection id for fitting or centering>                                   \\"
    puts "  -refselid          <selection id for reference>                                              \\"
    puts "  -prep_only         <whether analysis is performed or not (true or false)>                    \\"
    puts "                     (default: false)"
    puts ""
    puts "Usage:"
    puts "anatra spatial_distr                         \\"
    puts "  -stype             parm7                   \\"
    puts "  -sfile             str.prmtop              \\"
    puts "  -tintype           dcd                     \\"
    puts "  -tin               inp.dcd                 \\"
    puts "  -fhead             out                     \\"
    puts "  -sel0              not water               \\"
    puts "  -mode              residue                 \\"
    puts "  -ng3               50 50 50                \\"
    puts "  -del               0.4 0.4 0.4             \\"
    puts "  -origin            0.0 0.0 0.0             \\"
    puts "  -use_spline        false                   \\"
    puts "  -spline_resolution 4                       \\"
    puts "  -prep_only         false"
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
  set opt(mode)              "residue" 
  set opt(ng3)               "50 50 50" 
  set opt(del)               "0.5 0.5 0.5" 
  set opt(origin)            "0.0 0.0 0.0"
  set opt(use_pbcwrap)       false
  set opt(centertype)        "ZERO"
  set opt(use_weight)        false 
  set opt(use_conditional)   false
  set opt(ndim)              1
  set opt(react_range)       "0.0 0.0"
  set opt(fcv)               ""
  set opt(flist_cv)          ""
  set opt(fweight)           ""
  set opt(flist_weight)      ""
  set opt(use_spline)        false 
  set opt(spline_resolution) 4 
  set opt(count_threshold)   1.0e-10
  set opt(fit)               false
  set opt(refpdb)            ""
  set opt(fitselid)          0 
  set opt(refselid)          0


}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc read_optinfo {arglist} {
#-------------------------------------------------------------------------------

  global opt

  set opt(fhead)                [parse_arguments $arglist \
      "-fhead"                    "value" $opt(fhead)]
  set opt(mode)                 [parse_arguments $arglist \
      "-mode"                     "value" $opt(mode)]
  set opt(ng3)                  [parse_arguments $arglist \
      "-ng3"                      "value" $opt(ng3)]
  set opt(del)                  [parse_arguments $arglist \
      "-del"                      "value" $opt(del)]
  set opt(origin)               [parse_arguments $arglist \
      "-origin"                   "value" $opt(origin)]

  # Hidden options
  set opt(use_pbcwrap)          [parse_arguments $arglist \
      "-use_pbcwrap"              "value" $opt(use_pbcwrap)]
  set opt(centertype)           [parse_arguments $arglist \
      "-centertype"               "value" $opt(centertype)]

  set opt(use_weight)           [parse_arguments $arglist \
      "-use_weight"               "value" $opt(use_weight)]

  set opt(use_conditional)      [parse_arguments $arglist \
      "-use_conditional"          "value" $opt(use_conditional)]
  set opt(ndim)                 [parse_arguments $arglist \
      "-ndim"                     "value" $opt(ndim)]
  set opt(react_range)          [parse_arguments $arglist \
      "-react_range"              "value" $opt(react_range)]
  set opt(fcv)                  [parse_arguments $arglist \
      "-fcv"                      "value" $opt(fcv)]
  set opt(flist_cv)             [parse_arguments $arglist \
      "-flist_cv"                 "value" $opt(flist_cv)]
  set opt(fweight)              [parse_arguments $arglist \
      "-fweight"                  "value" $opt(fweight)]
  set opt(flist_weight)         [parse_arguments $arglist \
      "-flist_weight"             "value" $opt(flist_weight)]
  set opt(use_spline)           [parse_arguments $arglist \
      "-use_spline"               "value" $opt(use_spline)]
  set opt(spline_resolution)    [parse_arguments $arglist \
      "-spline_resolution"        "value" $opt(spline_resolution)]
  set opt(count_threshold)      [parse_arguments $arglist \
      "-count_threshold"          "value" $opt(count_threshold)]
  set opt(fit)                  [parse_arguments $arglist \
      "-fit"                      "value" $opt(fit)]
  set opt(refpdb)               [parse_arguments $arglist \
      "-refpdb"                   "value" $opt(refpdb)]
  set opt(fitselid)             [parse_arguments $arglist \
      "-fitselid"                 "value" $opt(fitselid)]
  set opt(refselid)             [parse_arguments $arglist \
      "-refselid"                 "value" $opt(refselid)]

}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc show_optinfo {} {
#-------------------------------------------------------------------------------

  global opt

  puts "<< option info >>"
  puts "fhead             = $opt(fhead)"
  puts "mode              = $opt(mode)"
  puts "ng3               = $opt(ng3)"
  puts "del               = $opt(del)"
  puts "origin            = $opt(origin)"

  if {$opt(use_pbcwrap)} {
    puts "use_pbcwrap       = $opt(use_pbcwrap)"
    puts "centertype        = $opt(centertype)"
  }

  puts "use_weight        = $opt(use_weight)"
  puts "use_conditional   = $opt(use_conditional)"
  puts "ndim              = $opt(ndim)"
  puts "react_range       = $opt(react_range)"
  puts "fcv               = $opt(fcv)"
  puts "flist_cv          = $opt(flist_cv)"
  puts "fweight           = $opt(fweight)"
  puts "flist_weight      = $opt(flist_weight)"
  puts "use_spline        = $opt(use_spline)"
  puts "spline_resolution = $opt(spline_resolution)"
  puts "count_threshold   = $opt(count_threshold)"
  puts "fit               = $opt(fit)"
  puts "fitselid          = $opt(fitselid)"
  puts "refselid          = $opt(refselid)"
  puts "refpdb            = $opt(refpdb)"
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
  set fort       "${anatra_path}/f90/bin/spatial_distr.x";list

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

  #set nf [molinfo $mol get numframes] 

  set fsdinp   [format "%s.sd.inp"     $opt(fhead)]
  set fsdout   [format "%s.sd.out"     $opt(fhead)]

  for {set isel 0} {$isel < $seltxt(nsel)} {incr isel} {
    set fmolinfo($isel) [format "%s.sd.%i.molinfo" $opt(fhead) $isel]
  }

  for {set isel 0} {$isel < $seltxt(nsel)} {incr isel} {
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
  
  if {$opt(fit)} {
    puts ""
    puts "--------------------"
    puts " Setup reference"
    puts "--------------------"
    puts ""
    set ref    [mol load pdb "$opt(refpdb)"]
    set refsel [atomselect $ref "$seltxt($opt(refselid))"]

    set refxyz [format "%s.sd.xyz" $opt(fhead)]
    
    animate write xyz $refxyz beg 0 end -1 waitfor all sel $refsel $ref
  } 

  set ntraj [llength $traj(tin)] 

  set f [open $fsdinp "w"]
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
  puts $f "   fcv          = \"$opt(fcv)\""
  puts $f "   flist_cv     = \"$opt(flist_cv)\""
  puts $f "   fweight      = \"$opt(fweight)\""
  puts $f "   flist_weight = \"$opt(flist_weight)\""
  puts $f " /"
  puts $f " &output_param"
  puts $f "   fhead     = \"$opt(fhead)\""
  puts $f " /"

  puts $f " &trajopt_param"
  puts $f "   molinfo   = \"$fmolinfo(0)\" \"$fmolinfo($opt(fitselid))\""
  puts $f " /"

  puts $f " &option_param"
  puts $f "   mode              = \"$opt(mode)\""
  puts $f "   ng3               = $opt(ng3)"
  puts $f "   del               = $opt(del)"
  puts $f "   origin            = $opt(origin)"
  puts $f "   use_pbcwrap       = .$opt(use_pbcwrap)."
  puts $f "   centertype        = \"$opt(centertype)\""
  puts $f "   use_weight        = .$opt(use_weight)."
  puts $f "   use_conditional   = .$opt(use_conditional)."
  puts $f "   ndim              = $opt(ndim)"
  puts $f "   react_range       = $opt(react_range)"
  puts $f "   use_spline        = .$opt(use_spline)."
  puts $f "   spline_resolution = $opt(spline_resolution)"
  puts $f "   count_threshold   = $opt(count_threshold)"
  puts $f "   fit               = .$opt(fit)."
  puts $f "   "
  puts $f " /"

  close $f

  
  #animate write dcd $fdcdtmp beg 0 end -1 waitfor all sel $sel(0) $mol

  if {!$common(prep_only)} {
    puts "SDF is calculated with ANATRA fortran program:"
    puts "$fort ..."
    puts "=== INPUT ==="
    set content [exec cat $fsdinp]
    puts $content
    puts "============="
    exec $fort $fsdinp >& $fsdout
    puts ""
    puts "=== OUTPUT ==="
    set content [exec cat $fsdout]
    puts $content
  } 
  puts "=============="
  puts ">> Finished"

  exit
}
#-------------------------------------------------------------------------------
