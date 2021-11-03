/****************************************************************
 SamTabix.H
 Copyright (C)2021 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_SamTabix_H
#define INCL_SamTabix_H
#include <iostream>
#include "BOOM/Pipe.H"
#include "SamReader.H"
using namespace std;
using namespace BOOM;

class SamTabix : public SamReader {
public:
  SamTabix(const String &tabix,const String &samFile,const String &substrate,
	   int begin,int end);
protected:
  void openPipe(const String &tabix,const String &samFile,
		const String &substrate,int begin,int end);
};

#endif

