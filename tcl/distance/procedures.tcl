#-------------------------------------------------------------------------------
proc print_title {} {
#-------------------------------------------------------------------------------

  puts "============================================================"
  puts ""
  puts "                   Distance Analysis"
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
    puts "anatra distance                                                               \\"
    puts "  -stype           <structure file type>                                      \\"
    puts "  -sfile           <structure file name>                                      \\"
    puts "  -tintype         <input trajectory file type>                               \\"
    puts "  -tin             <input trajectory file name>                               \\"
    puts "  -flist_traj      <trajectory file list (neccesary if tin is not specified)> \\"
    puts "  -fhead           <header of output file name>                               \\"
    puts "  -pbc             <treat pbc or not (true or false)>                         \\"
    puts "                   (default: false)                                           \\"
    puts "  -mode0           <analysis mode of sel0 (residue or whole or atom>          \\"
    puts "  -mode1           <analysis mode of sel1 (residue or whole or atom>          \\"
    puts "  -distance_type   <standard or minimum or intra>                             \\"
    puts "                   (default: standard)>                                       \\"
    puts "  -mindist_type0   <mindist type for species 0 (site or com)>                 \\" 
    puts "                   (default: site)>                                           \\"
    puts "  -mindist_type1   <mindist type for species 1 (site or com)                  \\"
    puts "                   (default: site)>                                           \\"
    puts "  -dt              <time step>                                                \\"
    puts "  -t_sta           <time at which analysis start (default: 0)>                \\"
    puts "  -t_end           <time at which analysis stop (default: 0 (till the end))>  \\"
    puts "  -sel0            <VMD selection> (X=0,1,2...)                               \\"
    puts "  -sel1            <VMD selection> (X=0,1,2...)                               \\"
    puts "  -prep_only       <where analysis is performed or not (true or false)        \\"
    puts "                   (default: false)"
    puts ""
    puts "Remark : "
    puts "  o if distance_type = standard => standard distance is calculated"
    puts "  o if distance_type = minimum  => minimum distance between pairs is calculated"
    puts ""
    puts "  o mindist_typeX specifies the minimum distance type"
    puts "    if mindist_typeX = site => atomic sites are analyzed"
    puts "    if mindist_typeX = com  => CoMs are analyzed"
    puts ""  
    puts "Usage:"
    puts "anatra distance                                    \\"
    puts "  -stype          parm7                            \\"
    puts "  -sfile          str.prmtop                       \\"
    puts "  -tintype        dcd                              \\"
    puts "  -tin            inp1.dcd inp2.dcd                \\"
    puts "  -fhead          out                              \\"
    puts "  -pbc            true                             \\"
    puts "  -distance_type  standard                         \\"
    puts "  -mode0          residue                          \\"
    puts "  -mode1          residue                          \\"
    puts "  -dt             0.1                              \\"
    puts "  -sel0           name C32  H2X H2Y and segid MEMB \\"
    puts "  -sel1           water                            \\"
    puts "  -prep_only      false"
    puts ""
    exit
  }

}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc define_optinfo {} {
#-------------------------------------------------------------------------------
  
  global opt

  set opt(fhead)          "out"
  set opt(dt)             0.1
  set opt(t_sta)          0.0
  set opt(t_end)          0.0
  set opt(pbc)            false
  set opt(distance_type)  "standard" 
  set opt(mindist_type0)  "site" 
  set opt(mindist_type1)  "site" 
  set opt(mode0)          "residue"
  set opt(mode1)          "residue"
}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
proc read_optinfo {arglist} {
#-------------------------------------------------------------------------------

  global opt

  set opt(fhead)         [parse_arguments $arglist \
      "-fhead"         "value" $opt(fhead)]
  set opt(dt)            [parse_arguments $arglist \
      "-dt"            "value" $opt(dt)]
  set opt(t_sta)         [parse_arguments $arglist \
      "-t_sta"         "value" $opt(t_sta)]
  set opt(t_end)         [parse_arguments $arglist \
      "-t_end"         "value" $opt(t_end)]
  set opt(dt)            [parse_arguments $arglist \
      "-dt"            "value" $opt(dt)]
  set opt(pbc)           [parse_arguments $arglist \
      "-pbc"           "value" $opt(pbc)]
  set opt(distance_type) [parse_arguments $arglist \
      "-distance_type" "value" $opt(distance_type)]
  set opt(mindist_type0) [parse_arguments $arglist \
      "-mindist_type0" "value" $opt(mindist_type0)]
  set opt(mindist_type1) [parse_arguments $arglist \
      "-mindist_type1" "value" $opt(mindist_type1)]
  set opt(mode0)         [parse_arguments $arglist \
      "-mode0"         "value" $opt(mode0)]
  set opt(mode1)         [parse_arguments $arglist \
      "-mode1"         "value" $opt(mode1)]
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc show_optinfo {} {
#-------------------------------------------------------------------------------

  global opt

  puts "<< option info >>"
  puts "fhead          = $opt(fhead)"
  puts "dt             = $opt(dt)"
  puts "t_sta          = $opt(t_sta)"
  puts "t_end          = $opt(t_end)"
  puts "pbc            = $opt(pbc)"
  puts "distance_type  = $opt(distance_type)"
  puts "mindist_type0  = $opt(mindist_type0)"
  puts "mindist_type1  = $opt(mindist_type1)"
  puts "mode0          = $opt(mode0)"
  puts "mode1          = $opt(mode1)"
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
  set fort       "${anatra_path}/f90/bin/distance.x";list

  # read trajectory
  #
  puts ""
  puts "--------------------"
  puts " Read trajectory"
  puts "--------------------"
  puts ""

  set mol 0;
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

  set nsel $seltxt(nsel)

  set finp   [format "%s.ds.inp" $opt(fhead) ]
  set fout   [format "%s.ds.out" $opt(fhead) ]

  for {set isel 0} {$isel < $nsel} {incr isel} {
    set fmolinfo($isel) [format "%s.ds.%i.molinfo" $opt(fhead) $isel]
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

  set f [open $finp "w"]
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
  puts $f "   dt        = $opt(dt)"
  puts $f "   molinfo   = \"$fmolinfo(0)\" \"$fmolinfo(1)\""
  puts $f " /"
  
  puts $f " &option_param"
  puts $f "   pbc            = .$opt(pbc)."
  puts $f "   mode           = \"$opt(mode0)\" \"$opt(mode1)\""
  puts $f "   distance_type  = \"$opt(distance_type)\""
  puts $f "   mindist_type   = \"$opt(mindist_type0)\" \"$opt(mindist_type1)\""
  puts $f "   t_sta          = $opt(t_sta)"
  puts $f "   t_end          = $opt(t_end)"
  puts $f " /"
  close $f

  if {!$common(prep_only)} {
    puts "PBC Distance is calculated with ANATRA fortran program:"
    puts "$fort ..."
    puts "=== INPUT ==="
    set content [exec cat $finp]
    puts $content
    puts "============="
    exec $fort $finp >& $fout
    puts ""
    puts "=== OUTPUT ==="
    set content [exec cat $fout]
    puts $content
  }

  puts "=============="
  puts ">> Finished"

  exit
}
#-------------------------------------------------------------------------------
