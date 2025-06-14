#-------------------------------------------------------------------------------
proc print_title {} {
#-------------------------------------------------------------------------------

  puts "============================================================"
  puts ""
  puts "                     Energy Analysis"
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
    puts "anatra interaction_energy                                                \\"
    puts "  -stype      <structure file type>                                      \\"
    puts "  -sfile      <structure file name>                                      \\"
    puts "  -tintype    <input trajectory file type>                               \\"
    puts "  -tin        <input trajectory file name>                               \\"
    puts "  -flist_traj <trajectory file list (neccesary if tin is not specified)> \\"
    puts "  -fhead      <header of output file name>                               \\"
    puts "  -parmformat <parameter file format (prmtop or anaparm)>                \\"
    puts "              (default: prmtop)                                          \\"
    puts "  -fanaparm   <anaparam file>                                            \\"
    puts "              (necessary if parmformat = anaparm)                        \\"
    puts "  -pbc        <treat pbc or not (true or false)>                         \\"
    puts "              (default: false)                                           \\"
    puts "  -calc_vdw   <whether vdw is calculated (true or false)>                \\"
    puts "              (default: true)                                            \\"
    puts "  -calc_elec  <whether elec is calculated (true or false)>               \\"
    puts "              (default: false)                                           \\"
    puts "  -vdw_type   <vdw interaction type (standard or attractive)>      \\" 
    puts "  -elec_type  <elec interaction type (bare or pme)>                \\"
    puts "              (default: bare)>                                     \\" 
    puts "  -dt         <time step>                                          \\"
    puts "  -rvdwcut    <LJ cutoff distance (A)>                             \\"
    puts "  -relcut     <ELEC cutoff distance (A)>                           \\"
    puts "  -pme_alpha  <PME screening parameter (A^-1)>                     \\"
    puts "  -pme_grids  <PME grids>                                          \\"
    puts "  -pme_rigid  <PME solute rigidity (true or false)>                \\"
    puts "              (default: false)                                     \\"
    puts "  -sel0       <VMD selection> (X=0,1,2...)                         \\"
    puts "  -sel1       <VMD selection> (X=0,1,2...)                         \\"
    puts "  -mode0      <analysis mode of sel0 (residue or whole or atom>    \\"
    puts "  -mode1      <analysis mode of sel1 (residue or whole or atom>    \\"
    puts "  -prep_only  <where analysis is performed or not (true or false)> \\"
    puts "              (default: false)"
    puts ""  
    puts "Usage:"
    puts "anatra interaction_energy                        \\"
    puts "  -stype        parm7                            \\"
    puts "  -sfile        str.prmtop                       \\"
    puts "  -tintype      dcd                              \\"
    puts "  -tin          inp.dcd                          \\"
    puts "  -fhead        out                              \\"
    puts "  -parmformat   prmtop                           \\"
    puts "  -pbc          true                             \\"
    puts "  -sel0         resname APR                      \\"
    puts "  -sel1         resname WAT                      \\"
    puts "  -mode0        whole                            \\"
    puts "  -mode1        residue                          \\"
    puts "  -calc_vdw     true                             \\"
    puts "  -calc_elec    false                            \\"
    puts "  -vdw_type     standard                         \\"
    puts "  -elec_type    bare                             \\"
    puts "  -dt           0.1                              \\"
    puts "  -rvdwcut      12.0                             \\"
    puts "  -relcut       1.0e10                           \\"
    puts "  -pme_alpha    0.35e0                           \\"
    puts "  -pme_grids    64 64 64                         \\"
    puts "  -pme_rigid    false                            \\"
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

  set opt(fhead)         "out"
  set opt(parmformat)    "prmtop"
  set opt(fanaparm)      "complex.anapram"
  set opt(dt)             0.1
  set opt(pbc)            false
  set opt(mode0)          "residue"
  set opt(mode1)          "residue"
  set opt(calc_vdw)       true 
  set opt(calc_elec)      false 
  set opt(vdw_type)       "standard"
  set opt(elec_type)      "bare"
  set opt(rvdwcut)        12.0 
  set opt(relcut)         1.0e10 
  set opt(pme_alpha)      0.35e0
  set opt(pme_grids)      "64 64 64"
  set opt(pme_rigid)      false
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc read_optinfo {arglist} {
#-------------------------------------------------------------------------------

  global opt

  set opt(fhead)      [parse_arguments $arglist \
      "-fhead"       "value" $opt(fhead)]
  set opt(parmformat) [parse_arguments $arglist \
      "-parmformat"  "value" $opt(parmformat)]
  set opt(fanaparm)   [parse_arguments $arglist \
      "-fanaparm"    "value" $opt(fanaparm)]
  set opt(dt)         [parse_arguments $arglist \
      "-dt"          "value" $opt(dt)]
  set opt(pbc)        [parse_arguments $arglist \
      "-pbc"         "value" $opt(pbc)]
  set opt(calc_vdw)   [parse_arguments $arglist \
      "-calc_vdw"    "value" $opt(calc_vdw)]
  set opt(calc_elec)  [parse_arguments $arglist \
      "-calc_elec"   "value" $opt(calc_elec)]
  set opt(vdw_type)   [parse_arguments $arglist \
      "-vdw_type"    "value" $opt(vdw_type)]
  set opt(elec)       [parse_arguments $arglist \
      "-elec_type"   "value" $opt(elec_type)]
  set opt(mode0)      [parse_arguments $arglist \
      "-mode0"       "value" $opt(mode0)]
  set opt(mode1)      [parse_arguments $arglist \
      "-mode1"       "value" $opt(mode1)]
  set opt(rvdwcut     [parse_arguments $arglist \
      "-rvdwcut"     "value" $opt(rvdwcut)]
  set opt(relcut)     [parse_arguments $arglist \
      "-relcut"      "value" $opt(relcut)]
  set opt(pme_alpha)  [parse_arguments $arglist \
      "-pme_alpha"   "value" $opt(pme_alpha)]
  set opt(pme_grids)  [parse_arguments $arglist \
      "-pme_grids"   "value" $opt(pme_grids)]
  set opt(pme_rigid)  [parse_arguments $arglist \
      "-pme_rigid"   "value" $opt(pme_rigid)]
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc show_optinfo {} {
#-------------------------------------------------------------------------------

  global opt

  puts "<< option info >>"
  puts "fhead          = $opt(fhead)"
  puts "parmformat     = $opt(parmformat)"
  puts "fanaparm       = $opt(fanaparm)"
  puts "dt             = $opt(dt)"
  puts "pbc            = $opt(pbc)"
  puts "mode0          = $opt(mode0)"
  puts "mode1          = $opt(mode1)"
  puts "calc_vdw       = $opt(calc_vdw)"
  puts "calc_elec      = $opt(calc_elec)"
  puts "vdw_type       = $opt(vdw_type)"
  puts "elec_type      = $opt(elec_type)"
  puts "rvdwcut        = $opt(rvdwcut)"
  puts "relcut         = $opt(relcut)"
  puts "pme_alpha      = $opt(pme_alpha)"
  puts "pme_grids      = $opt(pme_grids)"
  puts "pme_rigid      = $opt(pme_rigid)"
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
  set fort       "${anatra_path}/f90/bin/interaction_energy.x";list


  # check control parameter check
  #
  #if {$str(stype) != "parm7"} {
  #  puts "Error: only parm7 structure file is supported in this analysis."
  #  exit
  #}

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

  set feninp   [format "%s.en.inp" $opt(fhead) ]
  set fenout   [format "%s.en.out" $opt(fhead) ]

  for {set isel 0} {$isel < $nsel} {incr isel} {
    set fmolinfo($isel) [format "%s.en.%i.molinfo" $opt(fhead) $isel]
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

  set f [open $feninp "w"]
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
  puts $f "   fprmtop  = \"$str(sfile)\""
  puts $f ""
  puts $f "   fanaparm = \"$opt(fanaparm)\""
  puts $f " /"
  puts $f " &output_param"
  puts $f "   fhead = \"$opt(fhead)\""
  puts $f " /"
  
  puts $f " &trajopt_param"
  puts $f "   dt        = $opt(dt)"
  puts $f "   molinfo   = \"$fmolinfo(0)\" \"$fmolinfo(1)\""
  puts $f " /"
  
  puts $f " &option_param"
  puts $f "   parmformat   = \"$opt(parmformat)\""
  puts $f "   pbc          = .$opt(pbc)."
  puts $f "   mode         = \"$opt(mode0)\" \"$opt(mode1)\""
  puts $f "   calc_vdw     = .$opt(calc_vdw)."
  puts $f "   calc_elec    = .$opt(calc_elec)."
  puts $f "   vdw_type     = \"$opt(vdw_type)\""
  puts $f "   elec_type    = \"$opt(elec_type)\""
  puts $f "   rvdwcut      = $opt(rvdwcut)"
  puts $f "   relcut       = $opt(relcut)"
  puts $f "   pme_alpha    = $opt(pme_alpha)"
  puts $f "   pme_grids    = $opt(pme_grids)"
  puts $f "   pme_rigid    = .$opt(pme_rigid)."
  puts $f " /"
  close $f

  if {!$common(prep_only)} {
    puts "Interaction Energy is calculated with ANATRA fortran program:"
    puts "$fort ..."
    puts "=== INPUT ==="
    set content [exec cat $feninp]
    puts $content
    puts "============="
    exec $fort $feninp >& $fenout
    puts ""
    puts "=== OUTPUT ==="
    set content [exec cat $fenout]
    puts $content
  }

  puts "=============="
  puts ">> Finished"

  exit
}
#-------------------------------------------------------------------------------

