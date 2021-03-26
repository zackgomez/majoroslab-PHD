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
# This script simulates the effect of different switching error rates on
# estimation of theta
#=========================================================================

def simulateError(mat,pat,pi):
    A=[]; B=[]
    N=len(mat)
    switching=False
    for i in range(N):
        if(np.random.uniform(0,1)<pi):
            switching=not switching
        if(not switching):
            A.append(mat[i]); B.append(pat[i])
        else:
            A.append(pat[i]); B.append(mat[i])
    return (A,B)

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=6):
    exit(ProgramName.get()+" <sites-per-gene> <num-genes> <reads-per-site> <theta> <pi>\n")
(sitesPerGene,numGenes,readsPerSite,theta,pi)=sys.argv[1:]
sitesPerGene=int(sitesPerGene)
numGenes=int(numGenes)
readsPerSite=int(readsPerSite)
theta=float(theta)
pi=float(pi)
p=theta/(1+theta)

for i in range(numGenes):
    Mat=[]; Pat=[]
    for j in range(sitesPerGene):
        mat=np.random.binomial(readsPerSite,p)
        pat=readsPerSite-mat
        Mat.append(mat)
        Pat.append(pat)
    #trueTheta=float(sum(Mat))/float(sum(Pat))
    (A,B)=simulateError(Mat,Pat,pi)
    theta=float(sum(A))/float(sum(B))
    print(math.log2(theta),sep="\t")
