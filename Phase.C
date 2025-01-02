/****************************************************************
 Phase.C
 Copyright (C)2021 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "Phase.H"
using namespace std;

VariantPhase swapPhase(VariantPhase p) 
{
  switch(p) {
  case UNPHASED: throw "Error in VariantPhase swapPhase()";
  case IN_PHASE: return ANTI_PHASED;
  case ANTI_PHASED: return IN_PHASE;
  }
  throw "No case in VariantPhase swapPhase()";
}




