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



void VariantGraph::getComponents(Vector<ConnectedComponent> &into,
				 const IlluminaQual &Q,float confidence)
{
  ConnectedComponent comp;
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



void VariantGraph::phaseComponents(Vector<ConnectedComponent> &components)
{
  for(Vector<ConnectedComponent>::iterator cur=components.begin(),
	end=components.end() ; cur!=end ; ++cur)
    phaseComponent(*cur);
}



void VariantGraph::phaseComponent(ConnectedComponent &component)
{
  const int N=component.size();
  VariantPhase phase=IN_PHASE;
  for(int i=0 ; i<N ; ++i) {
    Variant &v=component[i];
    v.setComponentPhase(phase);
    v.setComponent(&component);
    switch(v.getPhase()) {
    case UNPHASED: 
      if(i!=N-1) throw "Unphased site in connected component";
      break;
    case IN_PHASE: break;
    case ANTI_PHASED: phase=swap(phase); break;
    default: throw "No case in VariantGraph::phaseComponent";
    }
  }
}



void VariantGraph::assignReads()
{
  for(Vector<ReadVariants>::iterator cur=reads.begin(), end=reads.end() ; 
      cur!=end ; ++cur) {
    ReadVariants &read=*cur;
    if(read.size()<1) return;
    VariantInRead &firstVar=read[0];
    VariantPhase compPhase=firstVar.v->getComponentPhase();
    Allele allele=firstVar.allele; // REF or ALT

    // Assign the read counts to the first variant in this connected component
    ConnectedComponent &component=*firstVar.v->getComponent();
    Variant &v=component[0];
    if(compPhase==IN_PHASE)
      ++v.getCount(allele);
    else
      ++v.getCount(swap(allele));
  }
}







