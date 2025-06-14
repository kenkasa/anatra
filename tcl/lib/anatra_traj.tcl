package provide anatra_traj 1.0

# --- Procedure list ---
# * read_traj
# * cut_traj
# * get_com
# * gen_pdbline
# * gen_xyzline

#-------------------------------------------------------------------------------
proc read_traj {molid stype sfile trajtype trajlist stride} {
#-------------------------------------------------------------------------------
  set molid [mol load $stype "$sfile"]
  foreach d $trajlist {
    animate read $trajtype $d beg 0 end -1 skip $stride waitfor all $molid 
  }

  if {$stype == "pdb" || $stype == "gro" || $stype == "xyz"} {
    animate delete beg 0 end 0 $molid
  }
}

proc read_traj_begend {molid stype sfile trajtype trajlist stride beg end} {
  set molid [mol load $stype "$sfile"]
  #foreach d $trajlist {
  #  animate read $trajtype $d beg $beg end $end skip $stride waitfor all $molid 
  #}
  foreach d $trajlist {
    animate read $trajtype $d beg 0 end -1 waitfor all $molid 
  }

  if {$stype == "pdb" || $stype == "gro" || $stype == "xyz"} {
    animate delete beg 0 end 0 $molid
  }

  if {$end != -1} {
    set delright [expr $end + 1]
    animate delete beg $delright end -1 $molid
  }

  if {$beg != 0} {
    set delleft [expr $beg - 1]
    animate delete beg 0 end $delleft $molid
  }

  if {$stride != 1} {
    set nf [molinfo $molid get numframes]
    for {set i $nf} {$i > 0} {set i [expr $i - 1] } {
      set chk [expr $i % $stride]

      if {$chk != 0} {
        set j [expr $i - 1]
        animate delete beg $j end $j $molid
      }
    }
  }

}

#-------------------------------------------------------------------------------
proc cut_traj {molid varname} {
#-------------------------------------------------------------------------------

  upvar $varname reactive
  # check the consistency
  set nf [molinfo $molid get numframes];list

  set ncut 0
  for {set istep 0} {$istep < $nf} {incr istep} {
    if {!$reactive($istep)} {
      set isnap [expr $istep - $ncut];list
      animate delete beg $isnap end $isnap $molid
      incr ncut;list
    }
  }
}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
proc get_com {sel} {
#-------------------------------------------------------------------------------

  global com_xyz;list

  set com [measure center $sel weight mass];list
  set com_xyz(0) [lindex $com 0];list
  set com_xyz(1) [lindex $com 1];list
  set com_xyz(2) [lindex $com 2];list
}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
proc gen_pdbline {aid anam rnam rid x y z} {
#-------------------------------------------------------------------------------
  if {$aid <= 99999} then {
    set aids [format "%5d" $aid]
  } else {
    set aids "*****"
  }
  #                          1   1   2    2
  #                          1   6   0    6
  set pline [format "ATOM  %5s %3s %4s  %4s    %8.3f%8.3f%8.3f" \
            $aids $anam $rnam $rid $x $y $z]
  return $pline
}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
proc gen_xyzline {aid anam x y z} {
#-------------------------------------------------------------------------------

  set xyzline [format "%8d  %5s  %15.7f  %15.7f  %15.7f" \
            $aid $anam $x $y $z]
  return $xyzline
}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
proc define_trajinfo {type} {
#-------------------------------------------------------------------------------

  global traj

  if {$type == "in"} {
    set traj(tintype)    "dcd";list
    set traj(tin)        "";list
    set traj(flist_traj) "";list
    set traj(stride)     1;list
    set traj(beg)        1;list
    set traj(end)        0;list

    set traj(beg_zerosta)  0;list
    set traj(end_zerosta) -1;list

  } elseif {$type == "out"} {
    set traj(totype)  "dcd";list
    set traj(to)      "out.dcd";list
  }

}
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
proc read_trajinfo {arglist type} {
#-------------------------------------------------------------------------------

  global traj 

  if {$type == "in"} { 
    set    traj(tintype) [parse_arguments $arglist \
                         "-tintype"    "value" $traj(tintype)];list
    set    traj(tin)     [parse_arguments $arglist \
                         "-tin"        "value" $traj(tin)];list
    set    traj(flist_traj)  [parse_arguments $arglist \
                         "-flist_traj" "value" $traj(flist_traj)];list
    set    traj(stride)  [parse_arguments $arglist \
                         "-stride"     "value" $traj(stride)];list
    set    traj(beg)     [parse_arguments $arglist \
                         "-beg"        "value" $traj(beg)];list
    set    traj(end)     [parse_arguments $arglist \
                         "-end"        "value" $traj(end)];list

    # determined from "beg" and "end" variables
    #
    set traj(beg_zerosta) [expr $traj(beg) - 1];list

    if {$traj(end) > 0} {
      set traj(end_zerosta) [expr $traj(end) - 1];list
    } else {
      set traj(end_zerosta) -1;list 
    }

  } elseif {$type == "out"} {
    set    traj(totype) [parse_arguments $arglist \
                         "-totype"  "value" $traj(totype)];list
    set    traj(to)     [parse_arguments $arglist \
                         "-to"      "value" $traj(to)];list
  }
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc show_trajinfo {type} {
#-------------------------------------------------------------------------------

  global traj 

  if {$type == "in"} {
    puts "<< input trajectory info >>"
    puts "tintype    = $traj(tintype)"
    puts "tin        = $traj(tin)"
    puts "flist_traj = $traj(flist_traj)"
    puts "stride     = $traj(stride)"
    puts "beg        = $traj(beg)"
    puts "end        = $traj(end)"
    puts ""
  } elseif {$type == "out"} {
    puts "<< output trajectory info >>"
    puts "totype     = $traj(totype)"
    puts "to         = $traj(to)"
    puts ""
  }
}
#-------------------------------------------------------------------------------
