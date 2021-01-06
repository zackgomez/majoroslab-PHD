/****************************************************************
 phd.C : Piled Higher & Deeper

 Copyright (C)2021 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/GffReader.H"
#include "BOOM/Regex.H"
#include "BOOM/Interval.H"
#include "BOOM/Vector.H"
#include "BOOM/VectorSorter.H"
#include "BOOM/Pipe.H"
#include "BOOM/Array2D.H"
#include "SamReader.H"
using namespace std;
using namespace BOOM;

struct Variant {
  String ID;
  int pos;
  char ref, alt;
  int genotype[2];
  Array2D<int> edges; // Edges to next variant; 0=ref, 1=alt
  Variant() {}
  Variant(String ID,int pos,char ref,char alt,int g[2])
    : ID(ID), pos(pos), ref(ref), alt(alt) 
  { genotype[0]=g[0]; genotype[1]=g[1]; }
};

class Application {
  Regex pseudogeneRegex;
  Regex chrRegex;
  String vcfFile;
  void processExons(Vector<GffFeature*> &exons,
		    Vector<Interval> &intervals,Vector<Variant> &,
		    String &substrate);
  void deleteExons(Vector<GffFeature*> &exons);
  void getIntervals(Vector<GffFeature*> &exons,Vector<Interval> &into);
  void getVariants(const String &substrate,
		   const Vector<Interval> &intervals,
		   Vector<Variant> &lines);
  bool parseVariant(const String &line,String &ID,int &pos,
		    char &cRef,char &cAlt,int *genotype);
  void filter(Vector<Variant> &,const Vector<Interval> &exons);
  bool find(const Variant &,const Vector<Interval> &exons);
  void processSam(SamReader &,Vector<Variant> &,Vector<Interval> &exons,
		  int geneBegin,int geneEnd,const String &substrate);
  int getLastPos(const Vector<Variant> &);
  int getLastPos(const Vector<Interval> &exons);
  void getGeneLimits(const Vector<GffFeature*> &exons,int &begin,
		     int &end);
  void addEdges(const SamRecord *read,Vector<Variant> &graph);

public:
  Application();
  int main(int argc,char *argv[]);
};


int main(int argc,char *argv[])
{
  try {
    Application app;
    return app.main(argc,argv);
  }
  catch(const char *p) { cerr << p << endl; }
  catch(const string &msg) { cerr << msg.c_str() << endl; }
  catch(const RootException &e)
    {cerr << e.getMessage() << endl; }
  catch(const exception &e)
    {cerr << "STL exception caught in main:\n" << e.what() << endl;}
  catch(...) { cerr << "Unknown exception caught in main" << endl; }
  return -1;
}



Application::Application()
  : pseudogeneRegex("pseudogene"),
    chrRegex("chr")
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=3)
    throw String("phd <indexed.vcf.gz> <sorted.gff> <in.sam>");
  vcfFile=cmd.arg(0);
  const String gffFile=cmd.arg(1);
  const String samFile=cmd.arg(2);

  // Open input files
  GffReader gff(gffFile);
  SamReader sam(samFile);

  // Process the GFF file line-by-line
  String currentGene;
  GffFeature *buffer=NULL;
  Vector<GffFeature*> exons;
  String chrom;
  while(true) {
    GffFeature *feature;
    if(buffer) { feature=buffer; buffer=NULL; }
    else { 
      feature=gff.nextFeature(); if(feature==NULL) break; 
      const String &geneID=feature->lookupExtra("gene_id");
      if(currentGene!="" && geneID!=currentGene) {
	buffer=feature;
	if(exons.size()==0) continue;
	Vector<Variant> variants;
	Vector<Interval> exonIntervals;
	int geneBegin, geneEnd;
	getGeneLimits(exons,geneBegin,geneEnd);
	String substrate;
	processExons(exons,exonIntervals,variants,substrate);
	if(substrate!=chrom)
	  { cout<<"CHROM "<<substrate<<endl; chrom=substrate; }
	deleteExons(exons);
	if(variants.size()==0) continue;
	processSam(sam,variants,exonIntervals,geneBegin,geneEnd,
		   substrate);
	//cout<<"============================================"<<endl;
	continue;
      }
    }
    currentGene=feature->lookupExtra("gene_id");
    if(feature->getFeatureType()!="exon" ||
       pseudogeneRegex.search(feature->lookupExtra("gene_type")) ||
       feature->lookupExtra("level")!="1")
      //feature->lookupExtra("transcript_support_level")!="1")
      { delete feature; continue; }
    exons.push_back(feature);
  }
  return 0;
}



void Application::getGeneLimits(const Vector<GffFeature*> &exons,
				int &begin,int &end)
{
  //cout<<exons.size()<<" exons"<<endl;
  begin=-1; end=-1;
  for(Vector<GffFeature*>::const_iterator cur=exons.begin(),
	End=exons.end() ; cur!=End ; ++cur) {
    const GffFeature *exon=*cur;
    int b=exon->getBegin(), e=exon->getEnd();
    if(begin==-1) { begin=b; end=e; continue; }
    if(b<begin) begin=b;
    if(e>end) end=e;
  }
}



int Application::getLastPos(const Vector<Variant> &variants)
{
  int pos=-1;
  for(Vector<Variant>::const_iterator cur=variants.begin(), end=
	variants.end() ; cur!=end ; ++cur) {
    const Variant &v=*cur;
    if(v.pos>pos) pos=v.pos;
  }
  return pos;
}



int Application::getLastPos(const Vector<Interval> &exons)
{
  int pos=-1;
  for(Vector<Interval>::const_iterator cur=exons.begin(), end=
	exons.end() ; cur!=end ; ++cur) {
    const Interval &v=*cur;
    if(v.getEnd()>pos) pos=v.getEnd();
  }
  return pos;
}


void Application::processSam(SamReader &sam,Vector<Variant> &variants,
			     Vector<Interval> &exons,int geneBegin,
			     int geneEnd,const String &substrate)
{
  static SamRecord *buffer=NULL;
  int kept=0;
  while(true) {
    SamRecord *rec;
    if(buffer) { rec=buffer; buffer=NULL; }
    else rec=sam.nextRecord();
    if(!rec) break;
    if(rec->flag_unmapped()) { delete rec; continue; }
    String readSubstrate=rec->getRefName();
    readSubstrate=chrRegex.substitute(readSubstrate,"");
    if(readSubstrate<substrate) { delete rec; continue; } // ### ???
    if(rec->getRefPos()+rec->getSequence().getLength()<geneBegin)
      { delete rec; continue; }
    if(rec->getRefPos()>geneEnd) { buffer=rec; break;}
    addEdges(rec,variants);
    delete rec;
    ++kept;
  }
  cout<<"KEPT READS: "<<kept<<endl;
}



void Application::processExons(Vector<GffFeature*> &exons,
			       Vector<Interval> &intervals,
			       Vector<Variant> &variants,
			       String &substrate)
{
  if(exons.size()==0) return;
  GffFeature &exon=*exons[0];
  substrate=exon.getSubstrate();
  substrate=chrRegex.substitute(substrate,"");
  getIntervals(exons,intervals);
  //cout<<"SUBSTRATE: "<<substrate<<endl;
  /*for(Vector<Interval>::iterator cur=intervals.begin(), 
	end=intervals.end() ; cur!=end ; ++cur)
	cout<<*cur<<endl;*/
  getVariants(substrate,intervals,variants);
  //cout<<"VARIANTS:"<<endl;
  for(Vector<Variant>::iterator cur=variants.begin(), end=variants.end() ;
      cur!=end ; ++cur) {
    Variant v=*cur;
    cout<<v.ID<<"\t"<<v.pos<<"\t"<<v.ref<<"\t"<<v.alt<<"\t"
	<<v.genotype[0]<<"|"<<v.genotype[1]<<endl;
  }
}



void Application::getIntervals(Vector<GffFeature*> &exons,
			       Vector<Interval> &into)
{
  for(Vector<GffFeature*>::iterator cur=exons.begin(), end=exons.end() ;
      cur!=end ; ++cur) {
    const GffFeature *feature=*cur;
    into.push_back(Interval(feature->getBegin(),feature->getEnd()));
  }

  IntervalComparator cmp;
  VectorSorter<Interval> sorter(into,cmp);
  sorter.sortAscendInPlace();
  Interval::Union(into);
}



void Application::deleteExons(Vector<GffFeature*> &exons)
{
  for(Vector<GffFeature*>::iterator cur=exons.begin(), end=exons.end() ;
      cur!=end ; ++cur)
    delete *cur;
  exons.clear();
}



void Application::getVariants(const String &substrate,
			      const Vector<Interval> &intervals,
			      Vector<Variant> &variants)
{
  for(Vector<Interval>::const_iterator cur=intervals.begin(), 
	end=intervals.end() ; cur!=end ; ++cur) {
    Interval interval=*cur;
    const String cmd=String("tabix ")+vcfFile+" "+substrate+":"
      +String(interval.getBegin())+"-"+String(interval.getEnd());
    Pipe pipe(cmd,"r");
    while(!pipe.eof()) {
      const String line=pipe.getline();
      if(line.length()>0 && line[0]=='#') continue;
      String ID; int pos; char ref,alt; int genotype[2];
      if(!parseVariant(line,ID,pos,ref,alt,genotype)) continue;
      Variant v(ID,pos,ref,alt,genotype);
      variants.push_back(v);
    }
    pipe.close();
  }
  filter(variants,intervals);
}



bool Application::parseVariant(const String &line,String &ID,int &pos,
			       char &cRef,char &cAlt,int *genotype)
{
  Vector<String> fields;
  line.getFields(fields,"\t");
  if(fields.size()<9 || fields[6]!="PASS" || fields[8]!="GT") 
    return false;
  String sGenotype=fields[9];

  // ### This is not fully generic:
  if(sGenotype=="0|1") { genotype[0]=0; genotype[1]=1; }
  else if(sGenotype=="1|0") { genotype[0]=1; genotype[1]=0; }
  else return false;

  pos=fields[1].asInt()-1; // Convert 1-based coord to 0-based
  ID=fields[2]; 
  String ref=fields[3], alt=fields[4];
  if(ref.length()!=1 || alt.length()!=1) return false;
  cRef=ref[0]; cAlt=alt[0];
  if(cRef!='A' && cRef!='C' && cRef!='G' && cRef!='T') return false;
  if(cAlt!='A' && cAlt!='C' && cAlt!='G' && cAlt!='T') return false;
  return true;
}



void Application::filter(Vector<Variant> &variants,
			 const Vector<Interval> &exons)
{
  int n=variants.size();
  for(int i=0 ; i<n ; ++i) {
    if(!find(variants[i],exons)) { variants.cut(i); --i; --n; }
  }
}



bool Application::find(const Variant &v,const Vector<Interval> &exons)
{
  for(Vector<Interval>::const_iterator cur=exons.begin(), end=exons.end();
      cur!=end ; ++cur)
    if((*cur).contains(v.pos)) return true;
  return false;
}



void Application::addEdges(const SamRecord *read,
			   Vector<Variant> &graph)
{
  const CigarString &cigar=read->getCigar();
  CigarAlignment &alignment=*cigar.getAlignment();
  const int L=alignment.length();
  int firstRefPos=alignment[0], lastRefPos=alignment[L-1]; // ### wrong
  const int offset=read->getRefPos();
  Interval refInterval(offset+firstRefPos,offset+lastRefPos+1);
  bool containsVariant=false;
  for(Vector<Variant>::iterator cur=graph.begin(), end=graph.end() ; 
      cur!=end ; ++cur) {
    if(refInterval.contains((*cur).pos)) { containsVariant=true; break; }
  }
  if(containsVariant) {
    CigarAlignment &inverse=*alignment.invert(lastRefPos-firstRefPos+1);
    for(Vector<Variant>::iterator cur=graph.begin(), end=graph.end() ; 
	cur!=end ; ++cur) {
      Variant &v=*cur;
      if(refInterval.contains(v.pos)) { 
	const int readPos=inverse[v.pos-offset];
	cout<<read->getID()<<" CONTAINS VARIANT AT POS "<<readPos<<endl;
	if(readPos<0) cout<<"\tUNMAPPED"<<endl;
	else {
	  const String &seq=read->getSequence();
	  if(readPos>seq.length()) {
	    cout<<"COORD ERROR: "<<readPos<<">"<<seq.length()<<endl;
	    throw "COORD ERROR";
	  }
	  cout<<"\t"<<v.ID<<" ALLELE IN READ: "<<seq[readPos]
	      <<" REV="<<read->flag_revComp()<<endl;
	}
      }
    }
    delete &inverse;
  }

    // ### Have to adress soft masks in CigarAlignment (not implemented)

  delete &alignment;
}

