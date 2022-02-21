/****************************************************************
 VcfStream.C
 Copyright (C)2022 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "VcfStream.H"
using namespace std;
using namespace BOOM;


VcfStream::VcfStream(const String &filename)
  : reader(filename), buffered(false)
{
  // ctor
}



VcfStream::~VcfStream()
{
  close();
}



void VcfStream::close()
{
  reader.close();
}



void VcfStream::getVariants(const Interval &interval,Vector<::Variant> &into)
{
  // PRECONDITION: VCF file contains only variants on the required chromosome
  // PRECONDITION: VCF file is sorted by chrom position
  
  const int begin=interval.getBegin(), end=interval.getEnd();
  BOOM::VariantAndGenotypes vg;
  while(true) {
    if(buffered) { vg=buffer; buffered=false; }
    else if(!reader.nextVariant(vg)) break;
    if(vg.genotypes.size()!=1)
      throw RootException("Too many genotypes in VcfStream::getVariants()");
    const BOOM::Variant &bv=vg.variant;
    const int pos=bv.getPos()-1; // Because they are 1-based in the file
    if(pos<begin) continue;
    if(pos>=end) { buffer=vg; buffered=true; break; }
    if(bv.numAlleles()!=2)
      throw RootException("#alleles not 2 in VcfStream::getVariants()");
    // ### This code limited to SNPs!  i.e., alleles of length 1 only
    if(bv.getAllele(0).length()!=1 || bv.getAllele(1).length()!=1) continue;
    const BOOM::Genotype &bg=vg.genotypes[0];
    if(!bg.isHet()) continue;
    int genotype[2];
    genotype[0]=bg.getAllele(0); genotype[1]=bg.getAllele(1);
    ::Variant v(bv.getID(),pos,bv.getAllele(0)[0],bv.getAllele(1)[0],
		genotype);
    into.push_back(v);
  }
}

