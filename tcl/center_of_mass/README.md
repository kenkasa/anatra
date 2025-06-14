# Trajectory Convert (cm_analysis.tcl)

## はじめに

 分子重心の座標の計算や，MSDの計算を行うスクリプト． 
  
## 使い方

  具体例を下記に示す．
  他のスクリプトに共通するオプションについては，[このページ](../README.md)を参照すること．
  ```
  vmd -dispdev text -e tr_convert.tcl        \
      -stype   parm7                         \
      -sfile   str.prmtop                    \
      -tintype dcd                           \
      -tin     inp.dcd                       \
      -comfile out.com                       \
      -outmsd  true                          \
      -msddim  3                             \
      -dt      0.1                           \
      -msdrange  10                          \
      -msdfile out.msd                       \
      -msdavefile out.msdave                 \
      -sel0    not water                     \
      -mode    residue
  ```
  
  * -outcom <whether com file is generated (true or false) (default:true)>
 
  重心の時系列データを出力するかどうか． 
   
  * -comfile <output CoM file name>

  重心の時系列データを出力するファイル名． 
    
  * -outmsd <whether msd file is generated (true or false) (default:false)>
 
  MSDの解析を行うかどうか．
 
  * -msdfile <output msd file name>
 
  MSD（各分子ごと）を出力するファイルの名前．
  
    * -msdavefile <output averaged msd file name>
 
  MSD（全ての分子に対する平均）を出力するファイル名前．
  
  * -msddim <dimension for MSD analysis>
 
  MSDを計算する際の空間次元．3の場合は3次元空間に対するMSD，2の場合はxy平面に対するMSDを計算する．
  脂質膜系などで側方拡散(Lateral diffusion)を計算する場合は2にする．
  
  * -dt <time interval>
 
  時間刻み．msdfileの1列目(x軸，時間)を出力するのに用いられる．
 
  * -msdavefile <output msd ave file name>

  MSDの平均を出力するファイル名．  
  
  * -mode <analysis mode (residue or whole) (default: residue)>
 
  selectionで選んだ分子集団の残基ごとの重心を解析するか，全体の重心を解析するか． 
