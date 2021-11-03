/****************************************************************
 SamReader.H
 Copyright (C)2020 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_SamReader_H
#define INCL_SamReader_H
#include <iostream>
#include <utility>
#include "SamRecord.H"
#include "BOOM/File.H"
#include "BOOM/String.H"
#include "BOOM/Regex.H"
#include "BOOM/Set.H"
#include "SamRecord.H"
using namespace std;
using namespace BOOM;

class SamReader {
public:
  SamReader();
  SamReader(const String &filename); // sam file can be gzipped!
  virtual ~SamReader();
  virtual SamRecord *nextSequence() { return nextRecord(); }
  virtual SamRecord *nextRecord();
  virtual SamRecord *nextSeqAndText(String &line);
  virtual void close();
protected:
  Regex gzRegex;
  Set<char> headerChars;
  File *fh;
  Vector<String> headerLines;
};

#endif

