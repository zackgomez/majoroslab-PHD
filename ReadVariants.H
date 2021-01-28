/****************************************************************
 ReadVariants.H
 Copyright (C)2021 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_ReadVariants_H
#define INCL_ReadVariants_H
#include <iostream>
#include "VariantInRead.H"
using namespace std;
using namespace BOOM;

// This class represents the alleles in a single read.
class ReadVariants {
public:
  ReadVariants();
  void push_back(const VariantInRead &v) { variants.push_back(v); }
  int size() { return variants.size(); }
  Vector<VariantInRead>::iterator begin() { return variants.begin(); }
  Vector<VariantInRead>::iterator end() { return variants.end(); }
  VariantInRead &operator[](int i) { return variants[i]; }
  Vector<VariantInRead> &getVariants();
protedted:
  Vector<VariantInRead> variants; // The alleles in this read
};

#endif

