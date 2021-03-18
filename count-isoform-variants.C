/****************************************************************
 count-isoform-variants.C

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
#include "BOOM/SumLogProbs.H"
#include "Variant.H"
using namespace std;
using namespace BOOM;

class Application {
  Regex pseudogeneRegex;
  Regex chrRegex;
  String vcfFile;
  bool parseVariant(const String &line,String &ID,int &pos,char &cRef,char &cAlt,
		    int *genotype);
  void getVariants(const String &substrate,const Vector<Interval> &,
		   Vector<Variant> &);
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
  : pseudogeneRegex("pseudogene"), chrRegex("chr")
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=2)
    throw String("count-isoform-variants <indexed.vcf.gz> <sorted.gff>");
  vcfFile=cmd.arg(0);
  const String gffFile=cmd.arg(1);

  // Load GFF
  cout<<"Loading GFF..."<<endl;
  GffReader gff(gffFile);
  Vector<GffGene> &genes=*gff.loadGenes();
  cout<<"Done loading GFF"<<endl;

  // Process each gene
  for(Vector<GffGene>::iterator cur=genes.begin(), End=genes.end() ; 
      cur!=End ; ++cur) {
    GffGene &gene=*cur;
    const int begin=gene.getBegin(), end=gene.getEnd();
    const String ID=gene.getID(), substrate=gene.getSubstrate();
    String vcfSubstrate=String("chr")+substrate;
    const int numTranscripts=gene.numTranscripts();
    cout<<"processing gene "<<ID<<": "<<numTranscripts<<" transcripts "<<" on "
	<<vcfSubstrate<<endl;
    for(int i=0 ; i<numTranscripts ; ++i) {
      GffTranscript &transcript=gene.getIthTranscript(i);
      Vector<GffExon*> rawExons;
      transcript.getRawExons(rawExons);
      Vector<Interval> exonIntervals;
      GffExon::getIntervals(rawExons,exonIntervals);
      Vector<Variant> variants;
      getVariants(vcfSubstrate,exonIntervals,variants);
      cout<<variants.size()<<" variants found in this transcript's exons"<<endl;
      transcript.deleteExons(rawExons);
    }
  }

  delete &genes;
  return 0;
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



