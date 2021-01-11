/****************************************************************
 SamRecord.C
 Copyright (C)2020 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "SamRecord.H"
using namespace std;
using namespace BOOM;

SamRecord::SamRecord(const String &ID,const String &refName,int refPos,
		     String cigar,String seq,unsigned int flags,
		     const Vector<String> &tags,const String &qual)
  : ID(ID), refName(refName), refPos(refPos), CIGAR(cigar), seq(seq),
    flags(flags), tags(tags), qual(qual)
{
  // ctor
}



const String &SamRecord::getQualityScores() const
{
  return qual;
}



const String &SamRecord::getID() const
{
  return ID;
}



const SamCigarString &SamRecord::getCigar() const
{
  return CIGAR;
}



const String &SamRecord::getSequence() const
{
  return seq;
}


const String &SamRecord::getRefName() const
{
  return refName;
}



int SamRecord::getRefPos() const
{
  return refPos;
}



String SamRecord::getTag(const String &label) const
{
  Regex regex("^([^:]+):[^:]+:(\\S+)");
  for(Vector<String>::const_iterator cur=tags.begin(), end=tags.end() ;
      cur!=end ; ++cur) {
    const String &tag=*cur;
    if(!regex.match(tag)) throw String("Can't parse SAM tag: ")+tag;
    if(regex[1]==label) return regex[2];
    return "";
  }
}



const Vector<String> &SamRecord::getTags() const
{
  return tags;
}



void SamRecord::parseMDtag(Vector<String> &fields) const
{
  fields.clear();
  String md=getTag("MD");
  Regex rex1("^(\\d+)(.*)"), rex2("^([ACGT])(.*)"), rex3("^(\\^[ACGT]+)(.*)");
  while(md.length()>0) {
    if(rex1.search(md)) {
      fields.push_back(rex1[1]);
      md=rex1[2]; 
    }
    else if(rex2.search(md)) {
      fields.push_back(rex2[1]);
      md=rex2[2]; 
    }
    else if(rex3.search(md)) {
      fields.push_back(rex3[1]);
      md=rex3[2]; 
    }
    else throw String("Can't parse MD tag: ")+md;
  }
}



bool SamRecord::flag_hasMultipleSegments() const
{
  return bool(flags & 0x1);
}



bool SamRecord::flag_properlyAligned() const
{
  return bool(flags & 0x2);
}



bool SamRecord::flag_unmapped() const
{
  return bool(flags & 0x4);
}



bool SamRecord::flag_nextSegmentUnmapped() const
{
  return bool(flags & 0x8);
}



bool SamRecord::flag_revComp() const
{
  return bool(flags & 0x10);
}



bool SamRecord::flag_nextSegmentRevComp() const
{
  return bool(flags & 0x20);
}



bool SamRecord::flag_firstOfPair() const
{
  return bool(flags & 0x40);
}



bool SamRecord::flag_secondOfPair() const
{
  return bool(flags & 0x80);
}



bool SamRecord::flag_secondaryAlignment() const
{
  return bool(flags & 0x100);
}



bool SamRecord::flag_failedFilters() const
{
  return bool(flags & 0x200);
}



bool SamRecord::flag_PCRduplicate() const
{
  return bool(flags & 0x400);
}



bool SamRecord::flag_supplAlignment() const
{
  return bool(flags & 0x800);
}



/*
 FLAGS CODES:

   0x1 template having multiple segments in sequencing
   0x2 each segment properly aligned according to the aligner
 > 0x4 segment unmapped
 > 0x8 next segment in the template unmapped
 > 0x10 SEQ being reverse complemented
   0x20 SEQ of the next segment in the template being reverse complemented
 > 0x40 the first segment in the template
 > 0x80 the last segment in the template
   0x100 secondary alignment
   0x200 not passing filters, such as platform/vendor quality controls
 > 0x400 PCR or optical duplicate
   0x800 supplementary alignment
*/
