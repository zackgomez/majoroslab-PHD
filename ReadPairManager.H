/****************************************************************
 ReadPairManager.H
 Copyright (C)2021 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_ReadPairManager_H
#define INCL_ReadPairManager_H
#include <iostream>
#include "BOOM/Map.H"
#include "ReadVariants.H"
using namespace std;
using namespace BOOM;

class ReadPairManager
{
public:
  ReadPairManager();
  void Register(ReadVariants *);

protected:
  Map<String, ReadVariants *> idToRead;
};

#endif
