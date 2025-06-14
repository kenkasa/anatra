# Scd Order Parameter Analysis (sc_analysis.tcl)

## はじめに

 脂質膜系におけるOrder-Parameter Scdを計算するスクリプト．
 現時点では，CHARMM力場のみをサポート．
 (AMBER力場でも利用できるかもしれないが未確認) 
  
## 使い方

  具体例を下記に示す．
  他のスクリプトに共通するオプションについては，[このページ](../README.md)を参照すること．
  ```
  vmd -dispdev text -e tr_convert.tcl        \
      -stype   parm7                         \
      -sfile   str.prmtop                    \
      -tintype dcd                           \
      -tin     inp.dcd                       \
      -header  scd                           \
      -dt      0.1                           \
      -cindex  1 2 3                         \
      -sel0    name C32  H2X H2Y and segid MEMB \
      -sel1    name C33  H3X H3Y and segid MEMB \
      -sel2    name C34  H4X H4Y and segid MEMB

  ```
 
  * -selX <selection>

  Order-parameterを計算したい炭素及び水素原子に対するセレクション． 
  * -dt <time interval>
 
  スナップショット間の時間差

  * -header <header of output file name>

  出力するファイル名のヘッダー

  * -cindex <index for analyzed carbons>

  Scdを計算する炭素原子のラベル．計算結果には影響を与えないが，${header}ave.datのx軸として用いられる．
