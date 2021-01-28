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
using namespace std;
using namespace BOOM;

enum SNP_ALLELE {
  REF=0,
  ALT=1
};

struct VariantInRead {
  Variant &v;
  int pos; // position in read
  SNP_ALLELE allele;
  VariantInRead(Variant &v,int pos,SNP_ALLELE a);
};

#endif

