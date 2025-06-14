#-------------------------------------------------------------------------------
proc print_title {} {
#-------------------------------------------------------------------------------

  puts "============================================================"
  puts ""
  puts "             Radial-Distribution Analysis"
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
    puts "anatra radial_distr                                                           \\"
    puts "  -stype         <structure file type>                                        \\"
    puts "  -sfile         <structure file name>                                        \\"
    puts "  -tintype       <input trajectory file type>                                 \\"
    puts "  -tin           <input trajectory file name>                                 \\"
    puts "  -flist_traj    <trajectory file list (neccesary if tin is not specified)>   \\"
    puts "  -fhead         <header of output file name>                                 \\"
    puts "  -sel0          <VMD selection> (X=0,1,2...)                                 \\"
    puts "  -sel1          <VMD selection> (X=0,1,2...)                                 \\"
    puts "  -mode0         <analysis mode of sel0>                                      \\"
    puts "                 (residue or whole or atom)                                   \\"
    puts "                 (default: residue)>                                          \\"
    puts "  -mode1         <analysis mode of sel1>                                      \\"
    puts "                 (residue or whole or atom)                                   \\"
    puts "                 (default: residue)                                           \\"
    puts "  -dr            <delta r (angstrom)>                                         \\"
    puts "  -identical     <true or false>                                              \\"
    puts "                 (default: false)                                             \\"
    puts "  -normalize     <true or false>                                              \\"
    puts "                 (default: false)                                             \\"
    puts "  -separate_self <true or false>                                              \\"
    puts "                 (default: false)                                             \\"
    puts "  -prep_only     <whethere analysis is performed or not>                      \\"
    puts "                 (true or false)                                              \\"
    puts "                 (default: false)>"
    puts ""  
    puts "Usage:"
    puts "anatra radial_distr                               \\"
    puts "  -stype         parm7                            \\"
    puts "  -sfile         str.prmtop                       \\"
    puts "  -tintype       dcd                              \\"
    puts "  -tin           inp.dcd                          \\"
    puts "  -fhead         out                              \\"
    puts "  -sel0          name C32  H2X H2Y and segid MEMB \\"
    puts "  -sel1          water                            \\"
    puts "  -mode0         residue                          \\"
    puts "  -mode1         residue                          \\"
    puts "  -dr            0.4                              \\"
    puts "  -identical     false                            \\"
    puts "  -normalize     false                            \\"
    puts "  -separate_self false                            \\"
    puts "  -prep_only     false"
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
  set opt(mode0)          residue
  set opt(mode1)          residue
  set opt(dr)             0.1
  set opt(identical)      false
  set opt(normalize)      true 
  set opt(separate_self)  false
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc read_optinfo {arglist} {
#-------------------------------------------------------------------------------

  global opt

  set opt(fhead)     [parse_arguments $arglist \
      "-fhead"      "value" $opt(fhead)]
  set opt(mode0)     [parse_arguments $arglist \
      "-mode0"      "value" $opt(mode0)]
  set opt(mode1)     [parse_arguments $arglist \
      "-mode1"      "value" $opt(mode1)]
  set opt(dr)       [parse_arguments $arglist \
      "-dr"         "value" $opt(dr)]
  set opt(identical) [parse_arguments $arglist \
      "-identical"  "value" $opt(identical)]
  set opt(normalize) [parse_arguments $arglist \
      "-normalize"  "value" $opt(normalize)]
  set opt(separate_self) [parse_arguments $arglist \
      "-separate_self"  "value" $opt(separate_self)]
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc show_optinfo {} {
#-------------------------------------------------------------------------------

  global opt

  puts "<< option info >>"
  puts "fhead          = $opt(fhead)"
  puts "mode0          = $opt(mode0)"
  puts "mode1          = $opt(mode1)"
  puts "dr             = $opt(dr)"
  puts "identical      = $opt(identical)"
  puts "normalize      = $opt(normalize)"
  puts "separate_self  = $opt(separate_self)"
  puts ""

}

proc analyze {} {
  # in
  global str
  global traj
  global seltxt
  global opt
  global common

  global sel 

  set anatra_path $::env(ANATRA_PATH);list
  set fort       "${anatra_path}/f90/bin/radial_distr.x";list

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

  set nsel $seltxt(nsel)

  set frdinp   [format "%s.rd.inp" $opt(fhead) ]
  set frdout   [format "%s.rd.out" $opt(fhead) ]

  for {set isel 0} {$isel < $nsel} {incr isel} {
    set fmolinfo($isel) [format "%s.rd.%i.molinfo" $opt(fhead) $isel]
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

    #animate write dcd $fdcdtmp($isel) \
    #  beg 0 end -1 waitfor all sel $sel($isel) $mol 
  }

  set ntraj [llength $traj(tin)] 

  set f [open $frdinp "w"]
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
  puts $f "   fhead  = \"$opt(fhead)\""
  puts $f " /"
  
  puts $f " &trajopt_param"
  puts $f "   molinfo = \"$fmolinfo(0)\" \"$fmolinfo(1)\""
  puts $f " /"
  
  puts $f " &option_param"
  puts $f "   mode      = \"$opt(mode0)\" \"$opt(mode1)\""
  puts $f "   dr        = $opt(dr)"
  puts $f "   identical = .$opt(identical)."
  puts $f "   normalize = .$opt(normalize)."
  puts $f "   separate_self = .$opt(separate_self)."
  puts $f " /"
  close $f

  if {!$common(prep_only)} {
    puts "RDF is calculated with ANATRA fortran program:"
    puts "$fort ..."
    puts "=== INPUT ==="
    set content [exec cat $frdinp]
    puts $content
    puts "============="
    exec $fort $frdinp >& $frdout
    puts ""
    puts "=== OUTPUT ==="
    set content [exec cat $frdout]
    puts $content
  }

  puts "=============="
  puts ">> Finished"

  exit

}
#-------------------------------------------------------------------------------

