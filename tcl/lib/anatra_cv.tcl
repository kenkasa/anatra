package provide anatra_cv 1.0

#-------------------------------------------------------------------------------
proc define_cvinfo {} {
#-------------------------------------------------------------------------------

  global cv

  set cv(cvfile)  ""
  set cv(ncv)      1
  set cv(cvstride) 1
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc read_cvinfo {arglist} {
#-------------------------------------------------------------------------------

  global cv 

  set cv(cvfile)   [parse_arguments $arglist "-cvfile"   "value" $cv(cvfile)]
  set cv(ncv)      [parse_arguments $arglist "-ncv"      "value" $cv(ncv)]
  set cv(cvstride) [parse_arguments $arglist "-cvstride" "value" $cv(stride)]
  } 

}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc show_cvinfo {} {
#-------------------------------------------------------------------------------

  global cv 
 
  puts "<< cv info >>"
  puts "cvfile   = $cv(cvfile)"
  puts "ncv      = $cv(ncv)"
  puts "cvstride = $cv(cvstride)"

}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc read_cvfile {molid} {
#-------------------------------------------------------------------------------

  global cv

  set ncv $cv(ncv)
  set fcv [open "$cvfile" "r"];list
  set istep 0;list
  set jstep 0;list
  while {![eof $fcv]} {
    gets $fcv line
    if {[expr $istep % $cv(cvstride)] == 0} {
      for {set icv 1} {$icv <= cv} {incr icv} {
        set icv2 [expr $icv - 1]
        set cv(data($jstep,$icv2)) [lindex $line $icv]
      }
      incr jstep;list
    }
    incr istep
  }
  set jstep [expr $jstep - 1];list
  if {$cv(data($jstep,0)) == ""} {
    for {set icv 0} {$icv < $ncv} {incr icv} {
      unset cv(data($jstep,$icv));list
    }
  }
  #set nstep  [expr $jstep - 1];list
  #set nstep  $jstep;list
  set nstep  [expr [array size cv] / $ncv];list
  set nread  [array size cv];list
  set nf    [molinfo $molid get numframes];list
  if {$nstep != $nf} then {
    puts "ERROR: # of lines in cvfile is different"
    puts "       from # of frames"
    puts "       in cvfile: $nstep, in trajectory: $nf"
    exit
  }

  puts "# of CV    : $ncv"
  puts "# of lines : $nstep"

  close $fcv
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc read_cvfiles {cvfile ncv cvstride} {
#-------------------------------------------------------------------------------

  # out
  global cvs;list
  global cvids;list
  global ncvfile;list
  global nstepc;list


  set ncvfile [llength $cvfile];list
  for {set ifile 0} {$ifile < $ncvfile} {incr ifile} {
    set fnam [lindex $cvfile $ifile];list
    set fcv [open "$fnam" "r"];list
    set istep 0;list
    set jstep 0;list
    while {![eof $fcv]} {
      gets $fcv line
      if {[expr $istep % $cvstride] == 0} {
        set cvids($ifile,$jstep) $jstep;list
        for {set icv 1} {$icv <= $ncv} {incr icv} {
          set icv2 [expr $icv - 1]
          set cvs($ifile,$jstep,$icv2) [lindex $line $icv]
        }
        incr jstep;list
      }
      incr istep
    }
    set jstep [expr $jstep - 1];list
    if {$cvs($ifile,$jstep,0) == ""} {
      unset cvids($ifile,$jstep);list
      for {set icv 0} {$icv < $ncv} {incr icv} {
        unset cvs($ifile,$jstep,$icv);list
      }
      set jstep [expr $jstep - 1];list
    }

    set nstepc($ifile) $jstep;list

    puts "cvfile     : $fnam"
    puts "# of CV    : $ncv"
    puts "# of lines : $nstepc($ifile)"
    puts ""

    close $fcv
  } 

}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc search_reactive_conf {molid varname ncv cvminmax} {
#-------------------------------------------------------------------------------

  global reactive   ;list
  upvar $varname cv;list

  set nf [molinfo $molid get numframes];list

  set imm 0;list
  for {set icv 0} {$icv < $ncv} {incr icv} {
    set cvmm($icv,0) [lindex $cvminmax $imm];list
    incr imm;list
    set cvmm($icv,1) [lindex $cvminmax $imm];list
    incr imm;list
  }

  for {set istep 0} {$istep < $nf} {incr istep} {
    set is_reactive true;list
    for {set icv 0} {$icv < $ncv} {incr icv} {
      if {$cv($istep,$icv) < $cvmm($icv,0) || \
          $cv($istep,$icv) > $cvmm($icv,1)} then {
        set is_reactive false;list
      }
    }
    set reactive($istep) $is_reactive
  }
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc cut_cv {var_cv var_cvid var_reactive ncv} {
#-------------------------------------------------------------------------------

  global cv_cut  ;list
  global cvid_cut;list 
  upvar  $var_cv       cv      ;list
  upvar  $var_cvid     cvid    ;list
  upvar  $var_reactive reactive;list

  set nf [array size reactive];list
  set isnap 0;list
  for {set istep 0} {$istep < $nf} {incr istep} {
    if {$reactive($istep)} {
      set cvid_cut($isnap) $cvid($istep);list
      for {set icv 0} {$icv < $ncv} {incr icv} {
        set cv_cut($isnap,$icv) $cv($istep,$icv);list
      } 
      incr isnap;list 
    }
  }
}
#-------------------------------------------------------------------------------
