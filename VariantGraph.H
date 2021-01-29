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
#include "BOOM/IlluminaQual.H"
#include "Variant.H"
#include "ReadVariants.H"
#include "ConnectedComponent.H"
using namespace std;
using namespace BOOM;

class VariantGraph {
public:
  VariantGraph();
  Vector<Variant> &getVariants();
  void phase(const IlluminaQual &,float minQual);
  void getComponents(Vector<ConnectedComponent> &into,
		     const IlluminaQual &Q,float minQual);
  int size() { return variants.size(); }
  void push_back(const Variant &v) { variants.push_back(v); }
  void clear() { variants.clear(); }
  Variant &operator[](int i) { return variants[i]; }
  Vector<Variant>::iterator begin() { return variants.begin(); }
  Vector<Variant>::iterator end() { return variants.end(); }
  Vector<ReadVariants> &getReads() { return reads; }
  void phaseComponents(Vector<ConnectedComponent> &);
  //void assignReads(Vector<ConnectedComponent> &);
  void assignReads();
protected:
  void phaseComponent(ConnectedComponent &);
  //void assignReads(ConnectedComponent &);
  Vector<Variant> variants;
  Vector<ReadVariants> reads;
};


#endif

