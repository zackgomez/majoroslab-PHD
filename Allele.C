/****************************************************************
 Allele.C
 Copyright (C)2021 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "Allele.H"
using namespace std;
using namespace BOOM;

Allele swap(Allele allele)
{
  return allele==REF ? ALT : REF;
}





