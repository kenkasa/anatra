package provide anatra_selection 1.0

#-------------------------------------------------------------------------------
proc define_selinfo {} {
#-------------------------------------------------------------------------------

  global seltxt 
  global nselmax 
  
  set nselmax 10000

  set seltxt(nsel) 0
  for {set i 0} {$i < $nselmax} {incr i} {
    set seltxt($i) ""
  }

}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc read_selinfo {arglist} {
#-------------------------------------------------------------------------------

  global seltxt 
  global nselmax

  set seltxt(nsel) 0  
  for {set i 0} {$i < $nselmax} {incr i} {
    set selopt     [format "-sel%d" $i]
    set seltxt($i) [parse_arguments $arglist "$selopt" "value" $seltxt($i)];list
    if {$seltxt($i) != ""} {
      incr seltxt(nsel)
    }
  } 

}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc show_selinfo {} {
#-------------------------------------------------------------------------------

  global seltxt 
 
  puts "<< selection info >>"
  set i 0
  while {$seltxt($i) != ""} {
    set sels [format "sel%5d = " $i]
    puts "$sels $seltxt($i)"
    incr i
  }

}
#-------------------------------------------------------------------------------
