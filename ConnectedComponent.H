/****************************************************************
 ConnectedComponent.H
 Copyright (C)2021 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_ConnectedComponent_H
#define INCL_ConnectedComponent_H
#include <iostream>
#include "BOOM/Vector.H"
#include "Variant.H"
using namespace std;
using namespace BOOM;


class ConnectedComponent {
public:
  ConnectedComponent();
  Vector<Variant*> getVariants() { return variants; }
  int size() { return variants.size(); }
  void push_back(Variant &v) { variants.push_back(&v); }
  Variant &operator[](int i) { return *variants[i]; }
  Vector<Variant*>::iterator begin() { return variants.begin(); }
  Vector<Variant*>::iterator end() { return variants.end(); }
  void clear() { variants.clear(); }
protected:
  Vector<Variant*> variants;
};

#endif

