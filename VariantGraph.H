/****************************************************************
 VariantGraph.H
 Copyright (C)2021 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_VariantGraph_H
#define INCL_VariantGraph_H
#include <iostream>
#include "BOOM/Vector.H"
#include "Variant.H"
using namespace std;
using namespace BOOM;

class VariantGraph {
public:
  VariantGraph();
  Vector<Variant> &getVariants();
  void getComponents(Vector<VariantGraph> &into);
  int size() { return variants.size(); }
  void push_back(const Variant &v) { variants.push_back(v); }
  Variant &operator[](int i) { return variants[i]; }
  Vector<Variant>::iterator begin() { return variants.begin(); }
  Vector<Variant>::iterator end() { return variants.end(); }
protected:
  Vector<Variant> variants;
};


#endif

