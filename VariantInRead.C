/****************************************************************
 VariantInRead.C
 Copyright (C)2021 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "VariantInRead.H"
using namespace std;
using namespace BOOM;


VariantInRead::VariantInRead(Variant &v,int pos,Allele a,float p)
  : v(&v), pos(pos), allele(a), probCorrect(p)
{
  // ctor
}



