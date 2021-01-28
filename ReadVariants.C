/****************************************************************
 ReadVariants.C
 Copyright (C)2021 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "ReadVariants.H"
using namespace std;
using namespace BOOM;

ReadVariants::ReadVariants()
{
  // ctor
}



Vector<VariantInRead> &ReadVariants::getVariants()
{
  return variants;
}



bool ReadVariants::consistentWithPhase()
{
  const int L=variants.size();
  for(int i=0 ; i<L-1 ; ++i) {
    const VariantInRead &thisVar=variants[i], nextVar=variants[i+1];
    const Allele a1=thisVar.allele, a2=nextVar.allele;
    const VariantPhase readPhase=a1==a2 ? IN_PHASE : ANTI_PHASED;
    const VariantPhase phase=thisVar.v->getPhase();
    //if(phase==UNPHASED) ...?
    if(readPhase!=phase) return false;
  }
  return true;
}


