# ANATRA 1.0

[Manual (English)](./docs/anatra_manual_en.pdf)  
[Manual (Japanese)](./docs/anatra_manual_jpn.pdf)

**ANATRA** (*Ana*lyze *Tra*jectories) is a collection of Tcl/Fortran90 programs for analyzing trajectories obtained from Molecular Dynamics (MD) simulations.

We often need to analyze specific atoms, molecules, or amino acid residues of proteins in a system. However, writing a general-purpose program to extract such specific parts is not an easy task. This is because it is necessary to extract atom groups that satisfy multiple complex conditions simultaneously, such as:

> "residues numbered from 11 to 50, with residue name ALA, and only heavy atoms."

Creating a program that can interpret and extract such sets defined by intersections and unions of multiple conditions can be quite difficult.

On the other hand, **VMD** (*Visual Molecular Dynamics*), a tool widely used by MD users around the world, offers an excellent feature called **Atomselection** for selecting specific parts of a system. For instance, the above condition can be expressed in VMD$^{\prime}$s Atomselection syntax as:

```
resid 11 to 50 and resname ALA and noh
```

which is easy to understand intuitively.

**ANATRA** adopts a design that leverages VMD$^{\prime}$s Atomselection feature. As a result, users can use exactly the same syntax as in Atomselection for selecting specific parts of the system, eliminating the need to learn a new selection language—especially convenient for existing VMD users.

**ANATRA** consists of numerous Tcl/Fortran programs tailored to different types of analyses. To manage these programs in a unified way, a command-line interface named `"anatra"` is provided. Users can perform various analyses by executing:

```
$ anatra analysis_mode -option1 -option2 ...
```

In this workflow, the corresponding Tcl script is executed on VMD according to the specified analysis mode, utilizing Atomselection to identify the region of interest. The selected information is then passed as input to a Fortran program for analysis. All of this is handled internally, and users do not need to be aware of the intermediate steps—everything runs in a streamlined manner.
At the same time, since the Fortran programs themselves do not implement Atomselection directly, the source code remains simpler and easier to maintain. This also allows users to create new analysis programs with ease by building upon the libraries provided by **ANATRA**.
 
 ---
* **Project Leader/Main Developer**   
    * Kento Kasahara (Graduate School of Engineering Science, Univ of Osaka)

* **Contributor**   
    * Ren Masayama (Graduate School of Engineering Science, Univ. of Osaka)
    * Kazuya Okita (Graduate School of Engineering Science, Univ. of Osaka)  
    * Yuya Matsubara (Graduate School of Engineering Science, Univ. of Osaka)
    * Yusei Maruyama (Graduate School of Engineering Science, Univ. of Osaka)
    * Ryo Okabe (Graduate School of Engineering Science, Univ. of Osaka)
    * Yuki Yamashita (Graduate School of Engineering Science, Univ. of Osaka)
    * Nobuyuki Matubayasi (Graduate School of Engineering Science, Univ. of Osaka) 

---
ANATRA is distributed under the **GNU General Public License version 2.0 (GPL v2.0)**.  

**ANATRA**  
Copyright (C) 2025  Kento Kasahara (Univ. of Osaka) 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see https://www.gnu.org/licenses/.


This program has been developed with dependencies on several external libraries, which are included in the package.  
Each of these libraries is subject to its own license terms as described below.

* **NetCDF**  
    Copyright 2025 Unidata

    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

    3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission. 
 
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

* **xdrfile-1.1.4**  
    Copyright (c) 2009-2014, Erik Lindahl & David van der Spoel  
    All rights reserved.

    Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    https://ftp.gromacs.org/pub/contrib/


* **xdf.F90**  
    This routine was originally developed by Wes Barnett and distributed under the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.  
https://github.com/wesbarnett/libgmxfort  
Modified version of the routine was developed by Kai-Min Tu.  
https://github.com/kmtu/xdrfort


* **mt19937.f**  
    This routine was developed by Makoto Matsumoto and Takuji Nishimura and was distributed under the GNU Library General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.  
https://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/VERSIONS/FORTRAN/fortran.html
