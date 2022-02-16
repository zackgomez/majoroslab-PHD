#!/usr/bin/env python
#=========================================================================
# Copyright (C)William H. Majoros (bmajoros@alumni.duke.edu)
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import ProgramName

def getFile(chrom,files,outDir):
    fh=files.get(chrom)
    if(fh is None):
        filename=outDir+"/"+chrom+".gff"
        files[chrom]=fh=open(filename,"wt")
    return fh

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <in.gff> <out-dir> <skip-chr,skip-chr,...>\n")
(inFile,outDir,skipChrs)=sys.argv[1:]

files={}
skipChrs=set(skipChrs.split(","))
with open(inFile,"rt") as IN:
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)<1): continue
        chrom=fields[0]
        if(chrom in skipChrs): continue
        fh=getFile(chrom,files,outDir)
        print(line,end="",file=fh)



