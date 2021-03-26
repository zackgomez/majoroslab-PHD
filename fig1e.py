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
# This script simulates the effect of increasing effective sample size
# by doing joint inference across pedigrees of increasing size
#=========================================================================


#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=6):
    exit(ProgramName.get()+" <sites-per-gene> <num-genes> <reads-per-site> <theta> <#indivs>\n")
(sitesPerGene,numGenes,readsPerSite,theta,numIndiv)=sys.argv[1:]
sitesPerGene=int(sitesPerGene)
numGenes=int(numGenes)
readsPerSite=int(readsPerSite)
theta=float(theta)
p=theta/(1+theta)
numIndiv=int(numIndiv)

for i in range(numGenes):
    Mat=[]; Pat=[]
    for j in range(sitesPerGene):
        N=readsPerSite*numIndiv
        mat=np.random.binomial(N,p)
        pat=N-mat
        Mat.append(mat)
        Pat.append(pat)
    #trueTheta=float(sum(Mat))/float(sum(Pat))
    theta=float(sum(Mat))/float(sum(Pat))
    print(math.log2(theta),sep="\t")
