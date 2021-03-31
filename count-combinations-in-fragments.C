/****************************************************************
 count-combinations-in-fragments.C
 Copyright (C)2021 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Pipe.H"
#include "BOOM/Random.H"
#include "BOOM/Interval.H"
#include "BOOM/Regex.H"
#include "BOOM/Map.H"
#include "BOOM/Array1D.H"
#include "BOOM/Array2D.H"
using namespace std;
using namespace BOOM;

struct Variant {
  String ID;
  String chrom;
  int pos;
  Array1D<int> genotype;
  Variant(const String &ID,const String &chrom,int pos,const Array1D<int> &genotype)
  : ID(ID), chrom(chrom), pos(pos), genotype(genotype) {}
};

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
  bool parseVariant(const String &sampleID,const String &line,String &ID,int &pos,
		    Array1D<int> &genotype);
  void getVariants(const String &sampleID,const Interval &,
		   Vector<Variant> &variants);
  Interval pickInterval(int intervalLen,int chromLen);
  int countHets(const String &sampleID,const Interval &,int whichHap0or1);
  int countHetsRead(const String &sampleID,const Interval &);
  void countHets(const String &sampleID,const Interval &fragment,
		 const Interval &read1,const Interval &read2,int &fragHets,
		 int &readHets,int whichHap0or1);
  void parseHeader(const String &filename);
  void addToPairCount(const Variant &v1,const Variant &v2,int whichHap0or1);
private:
  Regex genotypeRegex;
  String chromName;
  String vcfFile;
  Map<String,int> sampleColumns; // # column indices in VCF file for samples
  Map<String,Array2D<int>*> allelePairCounts; // key = position x position
};


int main(int argc,char *argv[])
{
  try {
    Application app;
    return app.main(argc,argv);
  }
  catch(const char *p) { cerr << p << endl; }
  catch(const string &msg) { cerr << msg.c_str() << endl; }
  catch(const exception &e)
    {cerr << "STL exception caught in main:\n" << e.what() << endl;}
  catch(...) { cerr << "Unknown exception caught in main" << endl; }
  return -1;
}



Application::Application()
  : genotypeRegex("(.+)\\\\|(.+)")
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=7)
    throw String("count-combinations-in-fragments <indexed-vcf.gz> <chrom> <chrom-len> <fragment-size> <read-length> <num-fragments> <sample-IDs>");
  vcfFile=cmd.arg(0);
  chromName=cmd.arg(1);
  const int chromLen=cmd.arg(2).asInt();
  const int fragmentLen=cmd.arg(3).asInt();
  const int readLen=cmd.arg(4).asInt();
  const int numFragments=cmd.arg(5).asInt();
  const String &sampleIDs=cmd.arg(6);

  // Parse sampleIDs
  Vector<String> samples;
  sampleIDs.getFields(samples,",");
  const int numSamples=samples.size();
  parseHeader(vcfFile);

  // Simulate fragments
  for(int i=0 ; i<numFragments ; ++i) {
    // Generate fragment and read intervals
    String sample=samples[RandomNumber(numSamples)];
    Interval fragment=pickInterval(fragmentLen,chromLen);
    Interval read1=Interval(fragment.getBegin(),
			    fragment.getBegin()+readLen);
    Interval read2=Interval(fragment.getEnd()-readLen,
			    fragment.getEnd());
    int whichHap=RandomNumber(2);

    // Count het sites in fragment and reads
    int fragHets, readHets;
    countHets(sample,fragment,read1,read2,fragHets,readHets,whichHap);
    cout<<fragHets<<"\t"<<readHets<<endl;
  }

  return 0;
}



void Application::parseHeader(const String &filename)
{
  GunzipPipe pipe(filename);
  while(!pipe.eof()) {
    String line=pipe.getline();
    if(line.length()==0) continue;
    if(line[0]!='#') continue;
    Vector<String> fields;
    line.getFields(fields,"\t");
    const int n=fields.size();
    if(n==0) continue;
    if(fields[0]!="#CHROM") continue;
    for(int i=9 ; i<n ; ++i) {
      const String &sample=fields[i];
      sampleColumns[sample]=i;
    }
    pipe.close();
    return;
  }
}



Interval Application::pickInterval(int len,int chromLen)
{
  const int begin=RandomNumber(chromLen-len);
  return Interval(begin,begin+len);
}




void Application::getVariants(const String &sampleID,const Interval &interval,
			      Vector<Variant> &variants)
{
  const String cmd=String("tabix ")+vcfFile+" "+chromName+":"
    +String(interval.getBegin())+"-"+String(interval.getEnd());
  Pipe pipe(cmd,"r");
  while(!pipe.eof()) {
    const String line=pipe.getline();
    if(line.length()>0 && line[0]=='#') continue;
    String ID; int pos; Array1D<int> genotype(2);
    if(!parseVariant(sampleID,line,ID,pos,genotype)) continue;
    Variant v(ID,chromName,pos,genotype);
    variants.push_back(v);
  }
  pipe.close();
}



void Application::addToPairCount(const Variant &v1,const Variant &v2,int whichHap)
{
  String key=String(v1.pos)+" "+String(v2.pos);
  if(!allelePairCounts.isDefined(key)) {
    Array2D<int> &a=allelePairCounts[key]=new Array2D<int>(2,2);
    a.setAllTo(0);
  }
  Array2D<int> &a=*allelePairCounts[key];
  ++a[v1.genotype[whichHap]][v2.genotype[whichHap]];
}



bool Application::parseVariant(const String &sampleID,const String &line,
			       String &ID,int &pos,Array1D<int> &genotype)
{
  Vector<String> fields;
  line.getFields(fields,"\t");
  if(fields.size()<9 || fields[6]!="PASS" || fields[8]!="GT") 
    return false;
  const int index=sampleColumns[sampleID];
  String sGenotype=fields[index];

  if(!genotypeRegex.match(sGenotype)) 
    throw sGenotype+" : can't parse genotype";
  const String gt1=genotypeRegex[2], gt2=genotypeRegex[2];
  if(gt1.length()!=1 || gt2.length()!=1) return false;
  const bool het=gt1!=gt2;
  if(!het) return false;
  genotype[0]=gt1.asInt(); genotype[1]=gt2.asInt();
  if(genotype[0]>1 || genotype[1]>1) return false; // bi-allelic variants only

  pos=fields[1].asInt()-1; // Convert 1-based coord to 0-based
  ID=fields[2]; 
  return true;
}



int Application::countHets(const String &sampleID,const Interval &interval,
			   int whichHap)
{
  Vector<Variant> variants;
  getVariants(sampleID,interval,variants);
  //addToPairCount(const Variant &v1,const Variant &v2,int whichHap)
  return variants.size();
}



int Application::countHetsRead(const String &sampleID,const Interval &interval)
{
  Vector<Variant> variants;
  getVariants(sampleID,interval,variants);
  return variants.size();
}



void Application::countHets(const String &sampleID,const Interval &fragment,
			    const Interval &read1,const Interval &read2,
			    int &fragHets,int &readHets,int whichHap)
{
  fragHets=countHets(sampleID,fragment,whichHap);
  if(!read1.overlaps(read2)) 
    readHets=countHetsRead(sampleID,read1)+countHetsRead(sampleID,read2);
  else 
    readHets=countHetsRead(sampleID,Interval(read1.getBegin(),read2.getEnd()));
}



