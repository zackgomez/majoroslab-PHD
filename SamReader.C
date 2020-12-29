/****************************************************************
 SamReader.C
 Copyright (C)2020 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "SamReader.H"
#include "BOOM/Pipe.H"
using namespace std;
using namespace BOOM;

SamReader::SamReader(const String &filename)
  : gzRegex("\\.gz$")
{
  // ctor

  //if(Regex::search("\\\\.gz$",filename))
  //  throw "Gzipped files are not supported in SamReader at this time";
  //fh.open(filename,"r");
  fh=gzRegex.search(filename) ?
    new GunzipPipe(filename) : new File(filename);
  headerChars+=('@'); 
  headerChars+=('[');
}



SamReader::~SamReader()
{
  delete fh;
}



SamRecord *SamReader::nextRecord()
{
  String line;
  SamRecord *rec=nextSeqAndText(line);
  return rec;
}



SamRecord *SamReader::nextSeqAndText(String &line)
{
  line=fh->getline();
  if(line=="") return NULL;
  while(line!="" && headerChars.isMember(line[0])) {
    if(line[0]=="@") headerLines.push_back(line);
    line=fh->getline();
  }
  if(line=="") return NULL;
  Vector<String> fields;
  line.getFields(fields,"\t");
  if(fields.size()<11) 
    throw RootException("can't parse sam line: "+line);
  String ID=fields[0], flagStr=fields[1], refName=fields[2],
    refPosStr=fields[3], cigar=fields[5], seq=fields[9];
  if(cigar=="*") cigar="";
  //ID,flags,refName,refPos,mapQual,cigar,rnext,pnext,templateLen,seq,qual)=fields[:11]
  const int refPos=int(refPosStr)-1; // convert 1-based to 0-based
  const unsigned int flags=int(flagStr);
  //const CigarString CIGAR=(cigar);
  Vector<String> tags;
  for(int i=11 ; i<fields.size() ; ++i) tags.push_back(fields[i]);
  SamRecord *rec=new SamRecord(ID,refName,refPos,cigar,seq,flags,tags);
  return rec;
}



void SamReader::close()
{
  fh->close();
}


/*
# M03884:303:000000000-C4RM6:1:1101:1776:15706    99      chrX:31786371-31797409  6687    44      150M    =       6813    271     ATACTATTGCTGCGGTAATAACTGTAACTGCAGTTACTATTTAGTGATTTGTATGTAGATGTAGATGTAGTCTATGTCAGACACTATGCTGAGCATTTTATGGTTGCTATGTACTGATACATACAGAAACAAGAGGTACGTTCTTTTACA  BBBBFFFFFFFGGGGGEFGGFGHFHFFFHHHFFHHHFHFHHHGFHEDGGHFHBGFHGBDHFHFFFHHHHFHHHHHGHGFFBGGGHFHFFHHFFFFHHHHGHGFHHGFHGHHHGFHFFHHFHHFFGFFFFGGEHFFEHHFGHHHGHHHHFB  AS:i:300        XN:i:0  
 */
