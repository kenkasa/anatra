package provide anatra_common 1.0

#-------------------------------------------------------------------------------
proc define_commoninfo {} {
#-------------------------------------------------------------------------------

  global common 

  set common(prep_only)  false;list

}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc read_commoninfo {arglist} {
#-------------------------------------------------------------------------------

  global common 

  set common(prep_only)   [parse_arguments $arglist "-prep_only"   "value" $common(prep_only)];list

}
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
proc show_commoninfo {} {
#-------------------------------------------------------------------------------

  global common 
 
  puts "<< common info >>"
  puts "prep_only  = $common(prep_only)"

}
#-------------------------------------------------------------------------------
