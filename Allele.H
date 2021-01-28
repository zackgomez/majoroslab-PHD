/****************************************************************
 Allele.H
 Copyright (C)2021 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_Allele_H
#define INCL_Allele_H
#include <iostream>
using namespace std;

enum Allele {
  REF=0,
  ALT=1
};

Allele swap(Allele);

#endif

