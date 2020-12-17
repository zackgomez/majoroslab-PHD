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
  SamReader(const String &filename);
  SamRecord *nextSequence() { return nextRecord(); }
  SamRecord *nextRecord();
  SamRecord *nextSeqAndText(String &line);
  void close();
private:
  Set<char> headerChars;
  File fh;
  Vector<String> headerLines;
};

#endif

