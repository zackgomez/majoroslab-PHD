/****************************************************************
 VcfStream.H
 Copyright (C)2022 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_VcfStream_H
#define INCL_VcfStream_H
#include <iostream>
#include "Variant.H"
#include "BOOM/VcfReader.H"
#include "BOOM/Interval.H"
using namespace std;
using namespace BOOM;

class VcfStream {
  VcfReader reader;
  BOOM::VariantAndGenotypes buffer;
  bool buffered;
public:
  VcfStream(const String &filename);
  virtual ~VcfStream();
  void getVariants(const Interval &,Vector<::Variant> &into);
  void close();
};

#endif

