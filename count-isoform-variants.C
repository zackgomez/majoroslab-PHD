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
#include "BOOM/Set.H"
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
		   Set<int> &);
  int countIsoSpecVariants(Array1D< Set<int> > &variants,int &totalVariants);
  bool isInAll(int variantPos,Array1D< Set<int> > &variants);
  void Application::debug(Vector<Interval> &exons,Set<int> &);
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
  Vector<GffGene> &genes=*GffReader::loadGenes(gffFile);
  cout<<genes.size()<<" genes loaded"<<endl;

  // Process each gene
  int genesWithVariants=0, genesWithIsoSpecVariants=0, allVariants=0,
    allIsoSpec=0;
  for(Vector<GffGene>::iterator cur=genes.begin(), End=genes.end() ; 
      cur!=End ; ++cur) {
    GffGene &gene=*cur;
    const int begin=gene.getBegin(), end=gene.getEnd();
    const String ID=gene.getID(), substrate=gene.getSubstrate();
    String vcfSubstrate=chrRegex.substitute(substrate,"");
    const int numTranscripts=gene.numTranscripts();
    cout<<"processing gene "<<ID<<": "<<numTranscripts<<" transcripts "<<" on "
	<<vcfSubstrate<<endl;
    Array1D< Set<int> > variants(numTranscripts);
    for(int i=0 ; i<numTranscripts ; ++i) {
      GffTranscript &transcript=gene.getIthTranscript(i);
      Vector<GffExon*> rawExons;
      transcript.getRawExons(rawExons);
      Vector<Interval> exonIntervals;
      GffExon::getIntervals(rawExons,exonIntervals);
      getVariants(vcfSubstrate,exonIntervals,variants[i]);
      //debug(exonIntervals,variants[i]);
      transcript.deleteExons(rawExons);
    }
    int totalVariants;
    int numIsoSpec=countIsoSpecVariants(variants,totalVariants);
    cout<<"XXX "<<numIsoSpec<<" isoform-specific variants in this gene, out of "
	<<totalVariants<<endl;
    if(numIsoSpec>0) ++genesWithIsoSpecVariants;
    if(totalVariants>0) ++genesWithVariants;
    cout<<genesWithIsoSpecVariants<<" genes had isoform-specific variants, out of "
	<<genesWithVariants<<" genes having exonic het sites = "
	<<float(genesWithIsoSpecVariants)/float(genesWithVariants)<<endl;
    allVariants+=totalVariants;
    allIsoSpec+=numIsoSpec;
    cout<<allIsoSpec<<" isoform-specific variants out of "<<allVariants<<" total = "
	<<float(allIsoSpec)/float(allVariants)<<endl;
  }
  cout<<genesWithIsoSpecVariants<<" genes had isoform-specific variants, out of "
      <<genesWithVariants<<" genes having exonic het sites"<<endl;

  delete &genes;
  return 0;
}



void Application::debug(Vector<Interval> &exons,Set<int> &variants)
{
  cout<<"  TRANSCRIPT:"<<endl<<"\t";
  for(Vector<Interval>::iterator cur=exons.begin(), end=exons.end() ; cur!=end ; 
      ++cur) {
    Interval &interval=*cur;
    cout<<interval<<" ";
  }
  cout<<"\n  VARIANTS IN THAT TRANSCRIPT:"<<endl<<"\t"<<variants<<endl;
}



int Application::countIsoSpecVariants(Array1D< Set<int> > &variants,
				      int &totalVariants)
{
  int num=0;
  Set<int> all;
  for(int i=0 ; i<variants.size() ; ++i) all+=variants[i];
  totalVariants=all.size();
  for(Set<int>::iterator cur=all.begin(), end=all.end() ; cur!=end ; ++cur) {
    const int x=*cur;
    if(!isInAll(x,variants)) ++num;
  }
  return num;
}



bool Application::isInAll(int variantPos,Array1D< Set<int> > &variants)
{
  for(int i=0 ; i<variants.size() ; ++i)
    if(!variants[i].isMember(variantPos)) return false;
  return true;
}



void Application::getVariants(const String &substrate,
			      const Vector<Interval> &intervals,
			      Set<int> &variants)
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
      //Variant v(ID,pos,ref,alt,genotype);
      variants.insert(pos);
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



