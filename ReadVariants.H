/****************************************************************
 ReadVariants.H
 Copyright (C)2021 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_ReadVariants_H
#define INCL_ReadVariants_H
#include <iostream>
#include "BOOM/String.H"
#include "VariantInRead.H"
using namespace std;
using namespace BOOM;

// This class represents the alleles in a single read.
class ReadVariants {
public:
  ReadVariants();
  ReadVariants(const String &);
  const String &getID() const { return id; }
  bool isFirstOfPair() const { return firstOfPair; }
  void setFirstOfPair(bool);
  void setMate(ReadVariants *r) { mate=r; }
  ReadVariants *getMate() const { return mate; }
  void skip() { wantToSkip=true; }
  bool shouldSkip() { return wantToSkip; }
  void push_back(const VariantInRead &v) { variants.push_back(v); }
  bool hasVariants() { return size()>0; }
  int size() { return variants.size(); }
  Vector<VariantInRead>::iterator begin() { return variants.begin(); }
  Vector<VariantInRead>::iterator end() { return variants.end(); }
  VariantInRead &operator[](int i) { return variants[i]; }
  Vector<VariantInRead> &getVariants();
  bool consistentWithPhase();
protected:
  String id;
  bool firstOfPair;
  ReadVariants *mate;
  bool wantToSkip; // For paired reads, only one should be counted
  Vector<VariantInRead> variants; // The alleles in this read
};

#endif

