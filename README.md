# NSFinder
A small tool for identifying nonsynonymous amino acid substitutions from vcf files (output of Freebayes)    
NSFinder : Nonsynonymous Substitutions Finder

# Download
### Enter the following command in cmd: 
```
C:\Users\XXX\NSFinder>git clone https://github.com/RoderickNi/NSFinder.git
```

# How to use
### FIRST: Enter the following command in cmd:    
```
C:\Users\XXX\NSFinder>python NASFinder.py
```
### SECOND: Enter the appropriate file path based on the input prompt：
```
 __          __         _                                       _
 \ \        / /        | |                                     | |
  \ \  /\  / /    ___  | |   ___    ___    _ __ ___     ___    | |_    ___
   \ \/  \/ /    / _ \ | |  / __|  / _ \  | '_ ` _ \   / _ \   | __|  / _ \
    \  /\  /    |  __/ | | | (__  | (_) | | | | | | | |  __/   | |_  | (_) |
     \/  \/      \___| |_|  \___|  \___/  |_| |_| |_|  \___|    \__|  \___/


          _   _    _____   ______   _               _
         | \ | |  / ____| |  ____| (_)             | |
         |  \| | | (___   | |__     _   _ __     __| |   ___   _ __
         | . ` |  \___ \  |  __|   | | | '_ \   / _` |  / _ \ | '__|
         | |\  |  ____) | | |      | | | | | | | (_| | |  __/ | |
         |_| \_| |_____/  |_|      |_| |_| |_|  \__,_|  \___| |_|


                                          -----Author: Ruoyao Ni


#  NASFinder 1.0  (Nonsynonymous Substitution Finder)
#
#  Copyright 2021 Ruoyao Ni,Zoology,CAS <niruoyao@ioz.ac.cn>
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

<Enter q to Quit!>
#  Please input the path of reference sequence file (.fasta): C:\Users\XXX\NSFinder\TestReferenceSeq.fasta
#  Please input the path of SNP file (.vcf): C:\Users\XXX\NSFinder\TestVcf.vcf
#  Frequency greater than or equal to : 0.005
#  Result Output: C:\Users\XXX\NSFinder\TestResult.tab

<Enter q to Quit!>
#  Please input the path of reference sequence file (.fasta): q
THANKS FOR USING NASFinder!
```
# Additional information
Frequency greater than or equal to : (Only nonsynonymous substitutions greater than or equal to this frequency are displayed.)    

