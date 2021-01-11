/****************************************************************
 SamRecord.H
 Copyright (C)2020 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#ifndef INCL_SamRecord_H
#define INCL_SamRecord_H
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/SamCigarString.H"
using namespace std;
using namespace BOOM;


class SamRecord {
public:
  SamRecord(const String &ID,const String &refName,int refPos,String cigar,
	    String seq,unsigned int flags,const Vector<String> &tags,
	    const String &qual);
  const String &getID() const;
  const SamCigarString &getCigar() const;
  const String &getSequence() const;
  const String &getQualityScores() const;
  const String &getRefName() const;
  int getRefPos() const;
  const Vector<String> &getTags() const;
  void parseMDtag(Vector<String> &into) const;
  String getTag(const String &tag) const; // returns third field, e.g. "122G25" in MD:Z:122G25
  bool flag_hasMultipleSegments() const;
  bool flag_properlyAligned() const;
  bool flag_unmapped() const;
  bool flag_nextSegmentUnmapped() const;
  bool flag_revComp() const;
  bool flag_nextSegmentRevComp() const;
  bool flag_firstOfPair() const;
  bool flag_secondOfPair() const;
  bool flag_secondaryAlignment() const;
  bool flag_failedFilters() const;
  bool flag_PCRduplicate() const;
  bool flag_supplAlignment() const;
private:
  String ID;
  String refName; // name of reference sequence the read aligns to
  int refPos; // position in reference where alignment begins
  SamCigarString CIGAR; // CigarString
  String seq; // read sequence
  String qual; // quality score encoding
  unsigned int flags; // bitfield
  Vector<String> tags; // array of tags at end of record (MD:Z:122G25, NM:i:1, etc.)
};


/*
EXAMPLE SAM RECORD:

M03884:303:000000000-C4RM6:1:1101:1776:15706    99      chrX:31786371-31797409  6687    44      150M    =       6813    271     ATACTATTGCTGCGGTAATAACTGTAACTGCAGTTACTATTTAGTGATTTGTATGTAGATGTAGATGTAGTCTATGTCAGACACTATGCTGAGCATTTTATGGTTGCTATGTACTGATACATACAGAAACAAGAGGTACGTTCTTTTACA  BBBBFFFFFFFGGGGGEFGGFGHFHFFFHHHFFHHHFHFHHHGFHEDGGHFHBGFHGBDHFHFFFHHHHFHHHHHGHGFFBGGGHFHFFHHFFFFHHHHGHGFHHGFHGHHHGFHFFHHFHHFFGFFFFGGEHFFEHHFGHHHGHHHHFB  AS:i:300        XN:i:0  

 */

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

#endif

