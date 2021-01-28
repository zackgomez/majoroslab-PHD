/****************************************************************
 VariantInRead.H
 Copyright (C)2021 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_VariantInRead_H
#define INCL_VariantInRead_H
#include <iostream>
#include "Variant.H"
#include "Allele.H"
using namespace std;
using namespace BOOM;

struct VariantInRead {
  Variant *v;
  int pos; // position in read
  Allele allele;
  float probCorrect; // Probability that the base was called correctly
  VariantInRead(Variant &v,int pos,Allele a,float probCorrect);
  VariantInRead() : v(NULL) {}
};

#endif

