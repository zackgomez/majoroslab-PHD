#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Author: William H. Majoros (bmajoros@alumni.duke.edu)
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import math
import ProgramName
import numpy as np

#=========================================================================
# This script simulates the accuracy of pseudophasing
#=========================================================================

def pseudophase(mat,pat):
    A=[]; B=[]
    for i in range(len(mat)):
        x=mat[i]; y=pat[i]
        if(x<y):
            A.append(x)
            B.append(y)
        else:
            B.append(x)
            A.append(y)
    return (A,B)

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <#sites-per-gene> <#genes> <reads-per-site>\n")
(sitesPerGene,numGenes,readsPerSite)=sys.argv[1:]
sitesPerGene=int(sitesPerGene)
numGenes=int(numGenes)
readsPerSite=int(readsPerSite)

for i in range(numGenes):
    Mat=[]; Pat=[]
    for j in range(sitesPerGene):
        mat=np.random.binomial(readsPerSite,0.5)
        pat=readsPerSite-mat
        Mat.append(mat)
        Pat.append(pat)
        #print(mat,pat,sep="\t")
    trueTheta=float(sum(Mat))/float(sum(Pat))
    #if(trueTheta>1): trueTheta=1/trueTheta
    (A,B)=pseudophase(Mat,Pat)
    theta=float(sum(A))/float(sum(B))
    print(math.log2(trueTheta),math.log2(theta),sep="\t")
