/****************************************************************
 VariantGraph.C
 Copyright (C)2021 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "VariantGraph.H"
using namespace std;
using namespace BOOM;

VariantGraph::VariantGraph()
{
  // ctor
}


Vector<Variant> &VariantGraph::getVariants()
{
  return variants;
}



void VariantGraph::getComponents(Vector<VariantGraph> &into,
				 const IlluminaQual &Q,float confidence)
{
  VariantGraph comp;
  const int numVariants=size();
  for(int i=0 ; i<numVariants ; ++i) {
    Variant &v=variants[i];
    comp.push_back(v);
    if(!v.isPhased()) {
      into.push_back(comp);
      comp.clear();
    }
  }
  if(comp.size()>0) into.push_back(comp);
}



void VariantGraph::phase(const IlluminaQual &Q,float minQual)
{
  for(Vector<Variant>::iterator cur=begin(), end=this->end() ; cur!=end ; 
      ++cur) {
    Variant &v=*cur;
    if(v.nonzero()) v.setPhase(Q,minQual);
  }
}
