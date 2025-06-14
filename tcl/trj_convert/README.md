# Trajectory Convert (tr_convert.tcl)

## はじめに

  トラジェクトリのwrap，fittingや特定の分子種や分子集団の抽出を行うスクリプト．
  
## 使い方

  具体例を下記に示す．
  他のスクリプトに共通するオプションについては，[このページ](../README.md)を参照すること．
  ```
  vmds tr_convert.tcl \
    -stype parm7 -sfile str.prmtop \
    -tintype dcd -tin inp.dcd \
    -totype  dcd -to  out.dcd \
    -sel0 all \
    -sel1 resid 1 to 275 and name CA \
    -sel2 resid 1 to 275 and name CA \
    -fit  true \
    -wrap true \
    -refpdb ref.pdb \
    -outselid 0 \
    -fitselid 1 \
    -refselid 2 \
  ```
  
  * -fit <whether fitting is performed or not (true or false)>
  
    フィッティングを行うかどうか．trueの場合は，`-refpdb`，`-refsel`，`-fitselid`を指定する必要がある．
    
  * -wrap <wheter wrapping is performed or not (true or false)>
  
    ラッピングを行うかどうか．
    
  * -refpdb <reference pdb file name>
  
    参照構造のpdbファイル．フィッティングを行う場合に指定する必要あり．
    このファイルに含まれている構造（`-refsel`で指定）に重なるように，トラジェクトリ内の分子を並進回転させる．
    
  * -outselid <selection id>
  
    `-to`で指定したトラジェクトリファイルに出力させたい分子集団を指定したセレクションID．
    上の例では0を指定しているので，`-sel0 all`で指定した通り，系に含まれる全ての原子を出力させる．
    
  * -fitselid <selection id>
  
    トラジェクトリに含まれる分子集団のうち，フィットする部分を指定したセレクションID．
    上記の例では1を指定しているので，系の中の`resid 1 to 275 and name CA`がフィッティングに用いられる．
  
  * -refselid <selection for reference>
  
    参照構造の中で，フィッティングに使う部分を指定したセレクションID．フィッティングを行う場合に指定する必要あり．
