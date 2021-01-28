/****************************************************************
 phd.C : Piled Higher & Deeper

 Copyright (C)2021 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/Constants.H"
#include "BOOM/CommandLine.H"
#include "BOOM/GffReader.H"
#include "BOOM/Regex.H"
#include "BOOM/Interval.H"
#include "BOOM/Vector.H"
#include "BOOM/VectorSorter.H"
#include "BOOM/Pipe.H"
#include "BOOM/Array2D.H"
#include "BOOM/IlluminaQual.H"
#include "BOOM/SumLogProbs.H"
#include "SamReader.H"
#include "VariantGraph.H"
#include "VariantInRead.H"
using namespace std;
using namespace BOOM;

const float MIN_PROB_CORRECT=0.8;

class Application {
  IlluminaQual illumina;
  char MIN_QUAL; // This is the encoding charater, NOT the integer value!
  Regex pseudogeneRegex;
  Regex chrRegex;
  String vcfFile;
  int readsSeen, readsDiscarded, readsUnmapped, readsWrongChrom;
  int numConcordant, numNonzero;
  void processExons(Vector<GffFeature*> &exons,
		    Vector<Interval> &intervals,VariantGraph &,
		    String &substrate);
  void deleteExons(Vector<GffFeature*> &exons);
  void getIntervals(Vector<GffFeature*> &exons,Vector<Interval> &into);
  void getVariants(const String &substrate,
		   const Vector<Interval> &intervals,VariantGraph &);
  bool parseVariant(const String &line,String &ID,int &pos,
		    char &cRef,char &cAlt,int *genotype);
  void filter(VariantGraph &,const Vector<Interval> &exons);
  bool find(const Variant &,const Vector<Interval> &exons);
  void processSam(SamReader &,VariantGraph &,Vector<Interval> &exons,
		  int geneBegin,int geneEnd,const String &substrate);
  int getLastPos(const VariantGraph &);
  int getLastPos(const Vector<Interval> &exons);
  void getGeneLimits(const Vector<GffFeature*> &exons,int &begin,int &end);
  void addEdges(const SamRecord *read,VariantGraph &graph);
  void installEdges(Vector<VariantInRead> &,const String &readID,
		    const String &qualities);
  void findVariantsInRead(VariantGraph &,const SamRecord *,
			  CigarAlignment &,Vector<VariantInRead> &);
  void processGraph(VariantGraph &);
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
  : pseudogeneRegex("pseudogene"), chrRegex("chr"),
    readsSeen(0), readsDiscarded(0), readsUnmapped(0), readsWrongChrom(0),
    numConcordant(0), numNonzero(0)
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=4)
    throw String("phd <indexed.vcf.gz> <sorted.gff> <in.sam> <min-qual>");
  vcfFile=cmd.arg(0);
  const String gffFile=cmd.arg(1);
  const String samFile=cmd.arg(2);
  const int minQual=cmd.arg(3).asInt();
  MIN_QUAL=illumina.phredToChar(minQual);

  // Open input files
  GffReader gff(gffFile);
  SamReader sam(samFile);

  // Load GFF
  //cout<<"Loading GFF..."<<endl;
  //Vector<GffGene> *genes=gff.loadGenes();
  //cout<<"Done."<<endl;

  // Process the GFF file line-by-line
  String currentGene;
  GffFeature *buffer=NULL;
  Vector<GffFeature*> exons;
  String chrom;
  int prevGeneBegin=-1, prevGeneEnd=-1;
  while(true) {
    GffFeature *feature;
    if(buffer) { feature=buffer; buffer=NULL; }
    else { 
      feature=gff.nextFeature(); if(feature==NULL) break; 
      const String &geneID=feature->lookupExtra("gene_id");
      if(currentGene!="" && geneID!=currentGene) {
	buffer=feature;
	if(exons.size()==0) continue;
	VariantGraph variants;
	Vector<Interval> exonIntervals;
	int geneBegin, geneEnd;
	getGeneLimits(exons,geneBegin,geneEnd);
	//if(geneBegin<prevGeneEnd)
	if(geneEnd<prevGeneBegin)
	  throw RootException(String("GTF is not sorted: (")+prevGeneBegin+
			      ","+prevGeneEnd+") overlaps ("+geneBegin+","+
			      geneEnd+")");
	prevGeneBegin=geneBegin; prevGeneEnd=geneEnd;
	String substrate;
	processExons(exons,exonIntervals,variants,substrate);
	if(substrate!=chrom)
	  { cout<<"CHROM "<<substrate<<endl; chrom=substrate; }
	deleteExons(exons);
	if(variants.size()==0) continue;
	processSam(sam,variants,exonIntervals,geneBegin,geneEnd,
		   substrate);
	processGraph(variants);
	continue;
      }
    }
    currentGene=feature->lookupExtra("gene_id");
    /*if(feature->getFeatureType()!="exon" ||
       pseudogeneRegex.search(feature->lookupExtra("gene_type")) ||
       feature->lookupExtra("level")!="1")*/
    if(feature->getFeatureType()!="exon" ||
       pseudogeneRegex.search(feature->lookupExtra("gene_type")))
      //feature->lookupExtra("transcript_support_level")!="1")
      { delete feature; continue; }
    exons.push_back(feature);
  }
  cout<<readsSeen<<" reads seen, "<<readsDiscarded<<" discarded"<<endl;
  cout<<numConcordant<<" concordant edges out of "<<numNonzero<<" nonzero"
      <<endl;
  return 0;
}



void Application::getGeneLimits(const Vector<GffFeature*> &exons,
				int &begin,int &end)
{
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



int Application::getLastPos(const VariantGraph &graph)
{
  int pos=-1;
  Vector<Variant> &variants=graph.getVariants();
  for(Vector<Variant>::const_iterator cur=variants.begin(), end=
	variants.end() ; cur!=end ; ++cur) {
    const Variant &v=*cur;
    if(v.getPos()>pos) pos=v.getPos();
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


void Application::processSam(SamReader &sam,VariantGraph &graph,
			     Vector<Interval> &exons,int geneBegin,
			     int geneEnd,const String &substrate)
{
  static SamRecord *buffer=NULL;
  while(true) {
    SamRecord *rec;
    if(buffer) { rec=buffer; buffer=NULL; }
    else { rec=sam.nextRecord(); if(rec) ++readsSeen; }
    if(!rec) break;
    if(rec->flag_unmapped()) { delete rec; ++readsUnmapped; continue; }
    String readSubstrate=rec->getRefName();
    readSubstrate=chrRegex.substitute(readSubstrate,"");
    if(readSubstrate<substrate) 
      { delete rec; ++readsWrongChrom; continue; } // ### ???

    // This line is not quite correct:
    if(rec->getRefPos()+rec->getSequence().getLength()<geneBegin)
      { delete rec; ++readsDiscarded; continue; }

    if(rec->getRefPos()>geneEnd) { buffer=rec; break;}
    addEdges(rec,graph);
    delete rec;
  }
}



void Application::processExons(Vector<GffFeature*> &exons,
			       Vector<Interval> &intervals,
			       VariantGraph &graph,
			       String &substrate)
{
  if(exons.size()==0) return;
  GffFeature &exon=*exons[0];
  substrate=exon.getSubstrate();
  substrate=chrRegex.substitute(substrate,"");
  getIntervals(exons,intervals);
  getVariants(substrate,intervals,graph);
  /*
  for(Vector<Variant>::iterator cur=variants.begin(), end=variants.end() ;
      cur!=end ; ++cur) {
    Variant v=*cur;
    //cout<<v.ID<<"\t"<<v.pos<<"\t"<<v.ref<<"\t"<<v.alt<<"\t"
    //	<<v.genotype[0]<<"|"<<v.genotype[1]<<endl;
    }*/
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
			      VariantGraph &graph)
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
      graph.getVariants().push_back(v);
    }
    pipe.close();
  }
  filter(graph,intervals);
}



bool Application::parseVariant(const String &line,String &ID,int &pos,
			       char &cRef,char &cAlt,int *genotype)
{
  Vector<String> fields;
  line.getFields(fields,"\t");
  if(fields.size()<9 || fields[6]!="PASS" || fields[8]!="GT") 
    return false;
  String sGenotype=fields[9];

  // ### This is not fully generic, could be improved:
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



void Application::filter(VariantGraph &graph,
			 const Vector<Interval> &exons)
{
  Vector<Variant> &variants=graph.getVariants();
  int n=variants.size();
  for(int i=0 ; i<n ; ++i) {
    if(!find(variants[i],exons)) { variants.cut(i); --i; --n; }
  }
}



bool Application::find(const Variant &v,const Vector<Interval> &exons)
{
  for(Vector<Interval>::const_iterator cur=exons.begin(), end=exons.end();
      cur!=end ; ++cur)
    if((*cur).contains(v.getPos())) return true;
  return false;
}



void Application::addEdges(const SamRecord *read,
			   VariantGraph &graph)
{
  const CigarString &cigar=read->getCigar();
  CigarAlignment &alignment=*cigar.getAlignment();
  const int numNodes=graph.size();
  Vector<VariantInRead> variantsInRead;
  findVariantsInRead(graph,read,alignment,variantsInRead);
  //if(variantsInRead.size()>1) cout<<"FOUND "<<variantsInRead.size()<<" VARIANTS IN READ"<<endl;
  installEdges(variantsInRead,read->getID(),read->getQualityScores());

  // ### Have to adress soft masks in CigarAlignment (not implemented)
  delete &alignment;
}



void Application::findVariantsInRead(VariantGraph &graph,
				     const SamRecord *read,
				     CigarAlignment &alignment,
				     Vector<VariantInRead> &variants)
{
  const int L=alignment.length();
  const int offset=read->getRefPos();
  const String &seq=read->getSequence();
  const String &qual=read->getQualityScores();
  for(int readPos=0 ; readPos<L ; ++readPos) {
    const int refPos=alignment[readPos]+offset;
    if(refPos==CIGAR_UNDEFINED) continue;
    for(Vector<Variant>::iterator cur=graph.begin(), end=graph.end() ;
	cur!=end ; ++cur) {
      Variant &v=*cur;
      if(v.getPos()==refPos) {
	//	const int q=int(qual[readPos]);
	if(qual[readPos]<MIN_QUAL) continue;
	const char c=seq[readPos];
	SNP_ALLELE allele;
	if(c==v.getRef()) allele=REF;
	else if(c==v.getAlt()) allele=ALT;
	else {
	  const float pError=illumina.charToErrorProb(qual[readPos]);
	  cout<<read->getID()<<" ALLELE MISMATCH P(error)="<<pError
	      <<"="<<illumina.charToPhred(qual[readPos])<<"="
	      <<qual[readPos]<<" "<<c<<" NOT "<<v.getRef()
	      <<" NOR "<<v.getAlt()<<" VAR="<<v.getID()
	      <<" READ POS="<<readPos<<endl;
	  continue;
	}
	//cout<<"NOMISMATCH\t"<<illumina.charToErrorProb(qual[readPos])
	//    <<endl;
	variants.push_back(VariantInRead(v,readPos,allele));
      }
    }
  }
}



void Application::installEdges(Vector<VariantInRead> &variants,
			       const String &readID,
			       const String &qualities)
{
  const int N=variants.size();
  for(int i=0 ; i<N-1 ; ++i) {
    VariantInRead &thisVar=variants[i], &nextVar=variants[i+1];
    ++thisVar.v.getEdges()[thisVar.allele][nextVar.allele];
    Vector<pair<float,float> > &cell=
      thisVar.v.getProbCorrect()[thisVar.allele][nextVar.allele];
    const float p1=1-illumina.charToErrorProb(qualities[thisVar.pos]);
    const float p2=1-illumina.charToErrorProb(qualities[nextVar.pos]);
    if(!isFinite(p1) || !isFinite(p2))
      cout<<"PHRED: "<<qualities[thisVar.pos]<<"="
	  <<illumina.charToErrorProb(qualities[thisVar.pos])
	  <<" "<<qualities[nextVar.pos]
	  <<" "<<illumina.charToErrorProb(qualities[nextVar.pos])<<endl;
    cell.push_back(pair<float,float>(p1,p2));
  }
}



void Application::processGraph(VariantGraph &G)
{
  const int N=G.size();
  int totalEdges=0;
  for (int i=0 ; i<N-1 ; ++i) {
    const Variant &v=G[i];
    for(int j=0 ; j<2 ; ++j)
      for(int k=0 ; k<2 ; ++k) 
	totalEdges+=v.getEdges()[j][k];
    if(v.nonzero()) ++numNonzero;
    if(v.concordant()) ++numConcordant;
  }
  if(totalEdges<1) return;

  for (int i=0 ; i<N-1 ; ++i) {
    const Variant &v=G[i];
    if(v.concordant() || !v.nonzero()) continue;
    //if(!v.nonzero()) continue;
    cout<<"GRAPH:"<<endl;
    cout<<v.getID()<<"\t"<<v.getEdges()<<endl;
    const float prob=v.probInPhase(illumina);
    cout<<"IN-PHASE: "<<prob<<"  ANTI-PHASED: "<<1-prob<<endl;
    if(prob<MIN_PROB_CORRECT && 1-prob<MIN_PROB_CORRECT) 
      throw String(v.getID())+" UNRESOLVED";
  }
}


