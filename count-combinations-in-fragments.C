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
using namespace std;
using namespace BOOM;

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
  bool parseVariant(const String &sampleID,const String &line,String &ID,int &pos);
  void getVariants(const String &sampleID,const Interval &,Set<int> &variants);
  Interval pickInterval(int intervalLen,int chromLen);
  int countHets(const String &sampleID,const Interval &);
  void countHets(const String &sampleID,const Interval &fragment,
		 const Interval &read1,const Interval &read2,int &fragHets,
		 int &readHets);
  void parseHeader(const String &filename);
private:
  Regex genotypeRegex;
  String chromName;
  String vcfFile;
  Map<String,int> sampleColumns; // # column indices in VCF file for samples
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
  const int numSamples=sampleIDs.size();
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

    // Count het sites in fragment and reads
    int fragHets, readHets;
    countHets(sample,fragment,read1,read2,fragHets,readHets);
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
    if(line[0]!="#") continue;
    Vector<String> fields;
    line.getFields(fields,"\t");
    const int n=fields.size();
    if(n==0) continue;
    if(fields[0]!="#CHROM") continue;
    for(int i=9 ; i<n ; ++i) {
      const String &sample=fields[i];
      sampleColumns[sample]=i;
    }
  }
}



Interval Application::pickInterval(int len,int chromLen)
{
  const int begin=RandomNumber(chromLen-len);
  return Interval(begin,begin+len);
}




void Application::getVariants(const String &sampleID,const Interval &interval,
			      Set<int> &variants)
{
  const String cmd=String("tabix ")+vcfFile+" "+chromName+":"
    +String(interval.getBegin())+"-"+String(interval.getEnd());
  Pipe pipe(cmd,"r");
  while(!pipe.eof()) {
    const String line=pipe.getline();
    if(line.length()>0 && line[0]=='#') continue;
    String ID; int pos;
    if(!parseVariant(sampleID,line,ID,pos)) continue;
    variants.insert(pos);
  }
  pipe.close();
}



bool Application::parseVariant(const String &sampleID,const String &line,
			       String &ID,int &pos)
{
  Vector<String> fields;
  line.getFields(fields,"\t");
  if(fields.size()<9 || fields[6]!="PASS" || fields[8]!="GT") 
    return false;
  const int index=sampleColumns[sampleID];
  String sGenotype=fields[index];

  if(!genotypeRegex.match(sGenotype)) 
    throw sGenotype+" : can't parse genotype";
  const bool het=(genotypeRegex[1]!=genotypeRegex[2]);
  //cout<<sGenotype<<" = "<<(het ? "het"  : "hom")<<endl;
  if(!het) return false;

  pos=fields[1].asInt()-1; // Convert 1-based coord to 0-based
  ID=fields[2]; 
  return true;
}



int Application::countHets(const String &sampleID,const Interval &interval)
{
  Set<int> variants;
  getVariants(sampleID,interval,variants);
  return variants.size();
}



void Application::countHets(const String &sampleID,const Interval &fragment,
			    const Interval &read1,
			    const Interval &read2,int &fragHets,int &readHets)
{
  fragHets=countHets(sampleID,fragment);
  if(!read1.overlaps(read2)) 
    readHets=countHets(sampleID,read1)+countHets(sampleID,read2);
  else 
    readHets=countHets(sampleID,Interval(read1.getBegin(),read2.getEnd()));
}



