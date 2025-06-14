#-------1---------2---------3---------4---------5---------6---------7---------8
#
#> Script     
#! @brief    Orientation Analysis 
#!
#! @authors  Kento Kasahara (KK)
#
#  (c) Copyright 2020 Osaka Univ. All rights reserved.
#
#-------1---------2---------3---------4---------5---------6---------7---------8


# Load ANATRA library
set     anatra_path $::env(ANATRA_PATH);list
set     libpath     "${anatra_path}/tcl/lib";list
lappend auto_path   $libpath;list
package require anatra 1.0;list


# Loaad sub-programs
set     fpath "${anatra_path}/tcl/rot_analysis";list
source  ${fpath}/rot_proc.tcl;list

print_title
show_rot_usage $argv

puts ""
puts "\[Step 1\] Read Structure parameters"
puts ""
define_strinfo
read_strinfo $argv
show_strinfo

puts ""
puts "\[Step 2\] Read Trajectory parameters"
puts ""
# in
define_trajinfo in
read_trajinfo   $argv in
show_trajinfo   in
# out
#define_trajinfo out
#read_trajinfo   $argv out
#show_trajinfo   out 

puts ""
puts "\[Step 3\] Read Selection parameters"
puts ""
define_selinfo
read_selinfo $argv
show_selinfo

puts ""
puts "\[Step 4\] Read Option parameters"
puts ""
define_rotoptinfo
read_rotoptinfo $argv
show_rotopt

puts ""
puts "\[Step 5 \] Start Analysis"
puts ""
rot_analysis


exit
