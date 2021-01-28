/****************************************************************
 Phase.C
 Copyright (C)2021 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "Phase.H"
using namespace std;
using namespace BOOM;

VariantPhase swap(VariantPhase p) 
{
  switch(p) {
  case UNPHASED: throw "Error in VariantPhase swap()";
  case IN_PHASE: return ANTI_PHASED;
  case ANTI_PHASED: return IN_PHASE;
  }
  throw "No case in VariantPhase swap()";
}




