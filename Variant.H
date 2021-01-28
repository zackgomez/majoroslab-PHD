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
using namespace std;
using namespace BOOM;

enum VariantPhase {
  UNPHASED,
  IN_PHASE,   // this var is in phase with next var (ref-ref, alt-alt)
  ANTI_PHASED // this var is anti-phased with next var (ref-alt, alt-ref)
};

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
  const String &getID() { return ID; }
  int getPos() { return pos; }
  char getRef() { return ref; }
  char getAlt() { return alt; }
  int getGenotype(int i) { return genotype[i]; }
  Array2D<Vector<pair<float,float> > > &getProbCorrect() 
    { return probCorrect; }
  Array2D<int> &getEdges() { return edges; }
protected:
  String ID;
  int pos;
  char ref, alt;
  int genotype[2];
  VariantPhase phase;
  Array2D<Vector<pair<float,float> > > probCorrect; // P(base is correct)
  Array2D<int> edges; // Edges to next variant; 0=ref, 1=alt
  float logLikInPhase(const IlluminaQual &);
  float logLikAntiPhased(const IlluminaQual &);
  float logProd(int i,int j,const IlluminaQual &);
  float logProdSwapped(int i,int j,const IlluminaQual &);
};

#endif

