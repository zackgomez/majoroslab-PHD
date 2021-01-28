/****************************************************************
 Variant.H
 Copyright (C)2021 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_Variant_H
#define INCL_Variant_H
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/Array2D.H"
#include "BOOM/IlluminaQual.H"
#include "BOOM/Vector.H"
#include "Phase.H"
using namespace std;
using namespace BOOM;

class Variant {
public:
  Variant();
  Variant(String ID,int pos,char ref,char alt,int g[2]);
  bool concordant() const;
  bool nonzero() const;
  float probInPhase(const IlluminaQual &);
  void setPhase(const IlluminaQual &Q,float confidence);
  bool isPhased() { return phase!=UNPHASED; }
  VariantPhase getPhase() { return phase; }
  VariantPhase getComponentPhase() { return componentPhase; }
  void setComponentPhase(VariantPhase p) { componentPhase=p; }
  const String &getID() { return ID; }
  int getPos() { return pos; }
  char getRef() { return ref; }
  char getAlt() { return alt; }
  int getGenotype(int i) { return genotype[i]; }
  Array2D<Vector<pair<float,float> > > &getProbCorrect() 
    { return probCorrect; }
  Array2D<int> &getEdges() { return edges; }
  int &getCount(int allele) { return counts[allele]; }
protected:
  String ID;
  int pos;
  char ref, alt;
  int counts[2]; // indexed by Allele: REF=0, ALT=1
  int genotype[2];
  VariantPhase phase; // Relative to the *next* site!
  VariantPhase componentPhase; // Relative to first site in connected comp
  Array2D<Vector<pair<float,float> > > probCorrect; // P(base is correct)
  Array2D<int> edges; // Edges to next variant; 0=ref, 1=alt
  float logLikInPhase(const IlluminaQual &);
  float logLikAntiPhased(const IlluminaQual &);
  float logProd(int i,int j,const IlluminaQual &);
  float logProdSwapped(int i,int j,const IlluminaQual &);
};

#endif

