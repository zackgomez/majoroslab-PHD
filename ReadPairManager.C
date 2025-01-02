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
}

void ReadPairManager::Register(ReadVariants *r)
{
  const String &id = r->getID();
  if (!idToRead.isDefined(id))
  {
    idToRead[id] = r;
    r->setFirstOfPair(true);
  }
  else
  {
    ReadVariants *first = idToRead[id];

    if (first->getInterval().overlaps(r->getInterval()))
    {
      first->setMate(r);
      r->setFirstOfPair(false);
    }
  }
}
