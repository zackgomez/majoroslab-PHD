/****************************************************************
 ReadPairManager.C
 Copyright (C)2015 William H. Majoros (martiandna@gmail.com).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "ReadPairManager.H"
using namespace std;
using namespace BOOM;

ReadPairManager::ReadPairManager()
{
  // ctor
}



bool ReadPairManager::Register(ReadVariants *r)
{
  const String &id=r->getID();
  if(!idToRead.isDefined(id)) {
    idToRead[id]=r;
    r->setFirstOfPair(true);
    return true;
  }
  else {
    ReadVariants *first=idToRead[id]; // ### this is causing the problem
    first->setMate(r);
    r->setFirstOfPair(false);
    return false;
  }
}


