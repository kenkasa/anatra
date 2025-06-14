package provide anatra_arg 1.0

#-------------------------------------------------------------------------------
proc parse_arguments {arglist optname opttype optval} {
#-------------------------------------------------------------------------------

  set is_value true 
  if {$opttype == "value"} {
    set is_value true
  } elseif {$opttype == "flag"} {
    set is_value false
  }

  foreach opt $arglist {
    incr i
    if {$opt == $optname} {
      if {$is_value} {
        set optval [get_chunk $arglist $i];list
      } else {
        set optval 1
      }
    }
  }
  return $optval
}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc get_chunk {arglist ista} {
#-------------------------------------------------------------------------------

  set na [llength $arglist]; list
  set chunk ""; list
  for {set i $ista} {$i < $na} {incr i} {
    set word   [lindex $arglist $i];list
    set wlen   [string length $word];list
    set first  [string range $word 0 0];list
    if {$first != "-" || $wlen == 1} {
      lappend chunk [lindex $arglist $i];list
    } elseif {$first == "-"} {
      set second [string range $word 1 1];list
      set is_num false;
      for {set inum 0} {$inum <= 9} {incr inum} {
        if {$second == $inum} {
	  set is_num true;
	}
      }
      if {$is_num} {
        lappend chunk [lindex $arglist $i];list
      } else {
        break 
      } 
    }
  }
  return $chunk
}
#-------------------------------------------------------------------------------
