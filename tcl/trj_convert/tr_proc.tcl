proc print_title {} {
  puts "============================================================"
  puts ""
  puts "                    Trajectory Convert"
  puts ""
  puts "============================================================"
}

proc show_tr_usage {arglist} {

  set help false
  set help [parse_arguments $arglist \
      "-h" "flag" $help]

  if {$help} {
    puts "Usage:"
    puts "anatra trjconv                                                            \\"
    puts "  -stype      <structure file type>                                       \\"
    puts "  -sfile      <structure file name>                                       \\"
    puts "  -tintype    <input trajectory file type>                                \\"
    puts "  -tin        <input trajectory file name>                                \\"
    puts "  -totype     <output trajectory file type>                               \\"
    puts "  -to         <output trajectory file name>                               \\"
    puts "  -beg        <first frame to be read> (default: 1)                       \\"
    puts "  -end        <last frame to be read>                                     \\"
    puts "              (default: 0, correspoding to last frame))                   \\"
    puts "  -selX       <X-th VMD selection> (X=0,1,2...)                           \\"
    puts "  -fit        <fit is performed or not (true or false)>                   \\"
    puts "  -centering  <centering is performed or not (true or false)>             \\"
    puts "              (default: false)                                            \\"
    puts "  -wrap       <wrap is performed or not (true or false)>                  \\"
    puts "              (default: false)                                            \\"
    puts "  -wrapcenter <wrapping center (origin or com) (default: fragment)>       \\"
    puts "  -wrapcomp   <how to wrap molecules (residue, segid,chain, or fragment)> \\"
    puts "              (residue or segid or chain or fragment or nomp)             \\"
    puts "              (default: fragment)                                         \\"
    puts "  -refpdb     <reference pdb file name>                                   \\"
    puts "              (necessary if fit = true)                                   \\"
    puts "  -outselid   <selection id for output molecules>                         \\"
    puts "  -fitselid   <selection id for fitting or centering>                     \\"
    puts "  -refselid   <selection id for reference>"
    puts ""
    puts "Remark:"
    puts "o If you specify fit = true & wrap = true, wrapcomp is automatically changed to 'com'"
    puts ""
    puts "Example (fitting):"
    puts "anatra trjconv                           \\"
    puts "  -stype      parm7                      \\"
    puts "  -sfile      str.prmtop                 \\"
    puts "  -tintype    dcd                        \\"
    puts "  -tin        inp.dcd                    \\"
    puts "  -totype     dcd                        \\"
    puts "  -to         out.dcd                    \\"
    puts "  -beg        1                          \\"
    puts "  -end        150                        \\"
    puts "  -sel0       not water                  \\"
    puts "  -sel1       resid 1 to 275 and name CA \\"
    puts "  -sel2       resid 1 to 275 and name CA \\"
    puts "  -fit        true                       \\"
    puts "  -centering  false                      \\"
    puts "  -wrap       true                       \\"
    puts "  -wrapcenter origin                     \\"
    puts "  -wrapcomp   fragment                   \\"
    puts "  -refpdb     ref.pdb                    \\"
    puts "  -outselid   0                          \\"
    puts "  -fitselid   1                          \\"
    puts "  -refselid   2"
    puts ""
    exit
  }

}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      define_baoptinfo
#! @brief         define bond angle option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc define_troptinfo {} {
  
  global tropt

  set tropt(fit)          false
  set tropt(wrap)         false
  set tropt(centering)    false

  set tropt(wrapcenter)   origin
  set tropt(wrapcomp)     fragment
  set tropt(outselid)     0
  set tropt(fitselid)     0 
  set tropt(refselid)     0

  # hidden options
  #
  set tropt(output_sfile) false

  # optional
  set tropt(refpdb)     ""
}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      read_troptinfo
#! @brief         read trajectory option paramerters 
#! @authors       KK
#! @param[in]  arglist : argument list
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc read_troptinfo {arglist} {

  global tropt

  set tropt(fit)        [parse_arguments $arglist \
      "-fit"        "value" $tropt(fit)]
  set tropt(wrap)       [parse_arguments $arglist \
      "-wrap"       "value" $tropt(wrap)]
  set tropt(centering)  [parse_arguments $arglist \
      "-centering"  "value" $tropt(centering)]
  set tropt(wrapcenter) [parse_arguments $arglist \
      "-wrapcenter" "value" $tropt(wrapcenter)]
  set tropt(wrapcomp)   [parse_arguments $arglist \
      "-wrapcomp"   "value" $tropt(wrapcomp)]
  set tropt(outselid)   [parse_arguments $arglist \
      "-outselid"  "value" $tropt(outselid)]
  set tropt(fitselid)   [parse_arguments $arglist \
      "-fitselid"  "value" $tropt(fitselid)]
  set tropt(refselid)   [parse_arguments $arglist \
      "-refselid"  "value" $tropt(refselid)]
  set tropt(refpdb)     [parse_arguments $arglist \
      "-refpdb"    "value" $tropt(refpdb)]
  set tropt(output_sfile) [parse_arguments $arglist \
      "-output_sfile" "value" $tropt(output_sfile)]

  # Combination error check
  #
  if {$tropt(fit) && $tropt(centering)} {
    puts "ERROR: fit and centering can not be used at the same time."
    exit
  }

  if {$tropt(fit)} {
    if {$tropt(refpdb) == ""} {
      puts "ERROR: refpdb should be specified if fit option is used."
      exit
    }
  }

}

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      fitting 
#! @brief         perform fitting the structure to reference structure 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc fitting {molid seltxt_fit seltxt_out refid seltxt_ref} {
  set refsel [atomselect $refid $seltxt_ref]
  set fitsel [atomselect $molid $seltxt_fit]
  set outsel [atomselect $molid $seltxt_out]

  set nf [molinfo $molid get numframes]
   
  for {set i 0} {$i < $nf} {incr i} {
    if {[expr $i % 100] == 0 || [expr $i + 1] == $nf} {
      puts [format "%10d / %10d" $i [expr $nf - 1]] 
    } 
    $fitsel frame $i
    $outsel frame $i
    $outsel move [measure fit $fitsel $refsel] 
  } 

}	

#=======1=========2=========3=========4=========5=========6=========7=========8
#
#> Procedure      show_troptinfo
#! @brief         show trajectory option paramerters 
#! @authors       KK
#
#=======1=========2=========3=========4=========5=========6=========7=========8

proc show_troptinfo {} {

  global tropt

  puts "<< option info >>"
  puts "fit        = $tropt(fit)"
  puts "wrap       = $tropt(wrap)"
  puts "centering  = $tropt(centering)"
  if {$tropt(wrap)} {
    puts "wrapcenter = $tropt(wrapcenter)"
    puts "wrapcomp   = $tropt(wrapcomp)"
  }
  puts "outselid   = $tropt(outselid)"
  puts "fitselid   = $tropt(fitselid)"
  if {$tropt(fit) || $tropt(centering)} {
    puts "refselid   = $tropt(refselid)"
    puts "refpdb     = $tropt(refpdb)"
  }

  # hidden option
  #
  # puts "output_sfile = $tropt(output_sfile)"
  puts ""

}

proc tr_convert {} {
  # in
  global str
  global traj
  global seltxt
  global tropt

  global sel 

  # read trajectory
  #
  puts ""
  puts "--------------------"
  puts " Read trajectory"
  puts "--------------------"
  puts ""

  set mol 0;
  read_traj_begend $mol               \
                   $str(stype)        \
                   $str(sfile)        \
                   $traj(tintype)     \
                   $traj(tin)         \
                   $traj(stride)      \
                   $traj(beg_zerosta) \
                   $traj(end_zerosta)

  set nf   [molinfo $mol get numframes]

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

  # setup reference 
  #
  if {$tropt(fit)} {
    puts ""
    puts "--------------------"
    puts " Setup reference"
    puts "--------------------"
    puts ""
    set ref    [mol load pdb "$tropt(refpdb)"]
    set refsel [atomselect $ref "$seltxt($tropt(refselid))"]
    set comref [measure center $refsel weight mass]

    set fitsel [atomselect $mol "$seltxt($tropt(fitselid))"]
  } elseif {$tropt(centering)} {
    puts ""
    puts "--------------------"
    puts " Setup centering"
    puts "--------------------"
    puts ""
   
    set fitsel [atomselect $mol "$seltxt($tropt(fitselid))"]
  }

  # Convert 
  #
  puts ""
  puts "--------------------"
  puts " Start convert"
  puts "--------------------"
  puts ""

  set nf [molinfo $mol get numframes] 
  set outsel $sel($tropt(outselid))

  if {$tropt(wrap)} {
    package require pbctools

    if {$tropt(fit)} {
      for {set istep 0} {$istep < $nf} {incr istep} {
        $fitsel frame $istep
        $outsel frame $istep
        set com   [measure center $fitsel weight mass]
        set vdiff [vecsub $com $comref]
        $outsel moveby [vecscale -1.0 $vdiff]	
      }

      puts "o Start wrapping"
      if {$tropt(wrapcomp) != "nocomp"} {
        pbc wrap -molid      $mol                        \
                 -compound   $tropt(wrapcomp)            \
	               -center     com                         \
	               -centersel  "$seltxt($tropt(fitselid))" \
	               -all
      } else {
        pbc wrap -molid      $mol                        \
	               -center     com                         \
	               -centersel  "$seltxt($tropt(fitselid))" \
	               -all
      }
      puts "> Finish wrapping"
      puts ""

      puts "o Start fitting"
      fitting $mol \
              "$seltxt($tropt(fitselid))" \
              "$seltxt($tropt(outselid))" \
	      $ref \
	      "$seltxt($tropt(refselid))"
      puts "> Finish fitting"
      puts ""

    } elseif {$tropt(centering)} {

      for {set istep 0} {$istep < $nf} {incr istep} {
        $fitsel frame $istep
        $outsel frame $istep
        set com   [measure center $fitsel weight mass]
        $outsel moveby [vecscale -1.0 $com]	
      }
      puts "o Start wrapping"
      if {$tropt(wrapcomp) != "nocomp"} {
        pbc wrap -molid      $mol                        \
                 -compound   $tropt(wrapcomp)            \
	               -center     com                         \
	               -centersel  "$seltxt($tropt(fitselid))" \
                 -all
      } else {
        pbc wrap -molid      $mol                        \
	               -center     com                         \
	               -centersel  "$seltxt($tropt(fitselid))" \
	               -all

      }
      puts "> Finish wrapping"
      puts ""

    } else {
      puts "o Start wrapping"
      if {$tropt(wrapcomp) != "nocomp"} {
        pbc wrap -molid    $mol               \
                 -compound $tropt(wrapcomp)   \
                 -center   $tropt(wrapcenter) \
                 -all

      } else {
        pbc wrap -molid    $mol               \
                 -center   $tropt(wrapcenter) \
                 -all
      }
      puts "> Finish wrapping"
      puts ""
    }

  } else {

    if {$tropt(fit)} {
      puts "o Start fitting"
      fitting $mol \
              "$seltxt($tropt(fitselid))" \
              "$seltxt($tropt(outselid))" \
	      $ref \
	      "$seltxt($tropt(refselid))"
      puts "> Finish fitting" 
       
    } elseif {$tropt(centering)} {
      for {set istep 0} {$istep < $nf} {incr istep} {
        $fitsel frame $istep
        $outsel frame $istep
        set com   [measure center $fitsel weight mass]
        $outsel moveby [vecscale -1.0 $com]	
      }

    }

  }

  puts ""

  puts ">> Finish all conversion"
  puts ""

  # generate converted trajectory 
  #
  if {$traj(totype) != "xtc"} {
    animate write $traj(totype) \
                  $traj(to) \
                  beg 0 end -1 \
                  waitfor all \
                  sel $outsel $mol
  } elseif {$traj(totype) == "xtc"} {
    set rand [expr int((100000*rand()))]
    set ftmp [format "tr%06d.dcd" $rand]
    set catcrdlog [format "catcrd%06.log" $rand]

    animate write dcd              \
                  $ftmp            \
                  beg 0 end -1     \
                  waitfor all      \
                  sel $outsel $mol

    puts ""
    puts "Convert to XTC format file using CATCRD ..."
    puts ""
    exec catcrd -i $ftmp -o $traj(to) >& $catcrdlog
    set  content [exec cat $catcrdlog]
    exec rm -f $ftmp $catcrdlog
  }

  if {$tropt(output_sfile)} { 
    animate write pdb         \
            structure.pdb     \
	          beg 0 end 0       \
	          waitfor all       \
	          sel $outsel $mol
  }

}
