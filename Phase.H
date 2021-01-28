/****************************************************************
 Phase.H
 Copyright (C)2021 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_Phase_H
#define INCL_Phase_H
#include <iostream>
using namespace std;
using namespace BOOM;

enum VariantPhase {
  UNPHASED,
  IN_PHASE,   // this var is in phase with next var (ref-ref, alt-alt)
  ANTI_PHASED // this var is anti-phased with next var (ref-alt, alt-ref)
};

VariantPhase swap(VariantPhase);

#endif

