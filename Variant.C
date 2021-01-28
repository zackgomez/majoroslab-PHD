/****************************************************************
 Variant.C
 Copyright (C)2021 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <math.h>
#include "BOOM/SumLogProbs.H"
#include "Variant.H"
using namespace std;
using namespace BOOM;

Variant::Variant()
  : phase(UNPHASED)
{
  // ctor

  counts[0]=counts[1]=0;
}



Variant::Variant(String ID,int pos,char ref,char alt,int g[2])
  : ID(ID), pos(pos), ref(ref), alt(alt), edges(2,2),
    probCorrect(2,2), phase(UNPHASED)
{ 
  genotype[0]=g[0]; 
  genotype[1]=g[1]; 
  edges.setAllTo(0); 
  counts[0]=counts[1]=0;
}



bool Variant::concordant() const
{
  const int diag1=edges[0][0]+edges[1][1], diag2=edges[0][1]+edges[1][0];
  return diag1==0 && diag2>0 || diag1>0 && diag2==0;
}



bool Variant::nonzero() const
{
  int sum=0;
  for(int i=0 ; i<2 ; ++i)
    for(int j=0 ; j<2 ; ++j)
      if(edges[i][j]>0) return true;
  return false;
}



// This function uses Bayes' Thm with a uniform prior to compute the
// posterior probability of being in-phase rather than anti-phased
float Variant::probInPhase(const IlluminaQual &illumina)
{
  const float lik1=logLikInPhase(illumina);
  const float lik2=logLikAntiPhased(illumina);
  const float prior1=log(1), prior2=log(1); // uniform prior (for now)
  const float likPrior1=lik1+prior1;
  const float likPrior2=lik2+prior2;
  const float posterior=exp(likPrior1-sumLogProbs(likPrior1,likPrior2));
  return posterior;
}



float Variant::logLikInPhase(const IlluminaQual &Q)
{
  const float logLik=logProd(0,0,Q)+logProd(1,1,Q)+
    logProdSwapped(0,1,Q)+logProdSwapped(1,0,Q);
  //cout<<"  LOG LIK IN PHASE = "<<logLik<<endl;
  return logLik;
}



float Variant::logLikAntiPhased(const IlluminaQual &Q)
{
  const float logLik=logProd(0,1,Q)+logProd(1,0,Q)+
    logProdSwapped(0,0,Q)+logProdSwapped(1,1,Q);
  //cout<<"  LOG LIK ANTI-PHASED = "<<logLik<<endl;
  return logLik;
}



float Variant::logProd(int i,int j,const IlluminaQual &illumina)
{
  float sum=0;
  const Vector<pair<float,float> > &pairs=probCorrect[i][j];
  for(Vector<pair<float,float> >::const_iterator cur=pairs.begin(), 
	end=pairs.end() ; cur!=end ; ++cur) {
    const pair<float,float> &p=*cur;
    sum+=log(p.first*p.second+(1-p.first)*(1-p.second));
    //cout<<"\t\t\t("<<i<<","<<j<<") first="<<p.first<<" second="<<p.second<<endl;
  }
  //cout<<"\t\tLOGPROD="<<sum<<endl;
  return sum;
}



float Variant::logProdSwapped(int i,int j,const IlluminaQual &illumina)
{
  float sum=0;
  const Vector<pair<float,float> > &pairs=probCorrect[i][j];
  for(Vector<pair<float,float> >::const_iterator cur=pairs.begin(), 
	end=pairs.end() ; cur!=end ; ++cur) {
    const pair<float,float> &p=*cur;
    sum+=log(p.first*(1-p.second)+(1-p.first)*p.second);
    //cout<<"\t\t\t("<<i<<","<<j<<") first="<<p.first<<" second="<<p.second<<endl;
  }
  //cout<<"\t\tLOGPROD SWAPPED="<<sum<<endl;
  return sum;
}



void Variant::setPhase(const IlluminaQual &Q,float confidence)
{
  const float P=probInPhase(Q);
  if(P>=confidence) phase=IN_PHASE;
  else if(1-P>=confidence) phase=ANTI_PHASED;
  else phase=UNPHASED;
}
