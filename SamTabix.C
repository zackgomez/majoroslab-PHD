/****************************************************************
 SamTabix.C
 Copyright (C)2021 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "SamTabix.H"
using namespace std;
using namespace BOOM;

SamTabix::SamTabix(const String &tabix,const String &samFile,
		   const String &substrate,int begin,int end)
  : SamReader()
{
  // ctor

  openPipe(tabix,samFile,substrate,begin,end);
}



void SamTabix::openPipe(const String &tabix,const String &samFile,
			const String &substrate,int begin,int end)
{
  const String cmd=tabix+" "+samFile+" "+substrate+":"
    +String(begin)+"-"+String(end);
  fh=new Pipe(cmd,"r");
}



