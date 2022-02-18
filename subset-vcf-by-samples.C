/****************************************************************
 subset-vcf-by-sample.C
 Copyright (C)2022 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/VcfReader.H"
#include "BOOM/Pipe.H"
#include "BOOM/Regex.H"
using namespace std;
using namespace BOOM;

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
private:
  Regex gzRegex;
  int chromBegin, chromEnd;
  void getSampleIndices(const Vector<String> &wantIDs,
			const Vector<String> &sampleIDs,Vector<int> &indices);
  void emitHeaderLines(const Vector<String> &lines,File &);
  void emitChromLine(const Vector<String> &ids,File &);
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
  : gzRegex("gz$"), chromBegin(0), chromEnd(-1)
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"c:");
  if(cmd.numArgs()!=3)
    throw String("subset-vcf-by-sample [-c begin:end] <in.vcf> <sampleIDs> <out.vcf>");
  const String infile=cmd.arg(0);
  const String samples=cmd.arg(1);
  const String outfile=cmd.arg(2);
  Vector<String> wantIDs;
  samples.getFields(wantIDs,",");

  if(cmd.option('c')) {
    const String optparms=cmd.optParm('c');
    Regex r("(\\S+):(\\S+)");
    if(!r.match(optparms)) throw String("Invalid -c option");
    chromBegin=r[1].asInt();
    chromEnd=r[2].asInt();
    if(chromBegin<0 || chromBegin>=chromEnd) 
      throw String("Invalid coordinates in -c option");
  }

  // Open output file
  File &file=gzRegex.search(outfile) ? *new BGzipPipe(outfile)
    : *new File(outfile,"w");

  // Process input file
  VcfReader reader(infile);
  const Vector<String> &sampleIDs=reader.getSampleIDs();
  Vector<int> wantIndices;
  getSampleIndices(wantIDs,sampleIDs,wantIndices);
  if(wantIndices.size()==0) throw String("Can't find samples in VCF file");
  const Vector<String> &headerLines=reader.getHeaderLines();
  emitHeaderLines(headerLines,file);
  const String &chromLine=reader.getChromLine();
  emitChromLine(wantIDs,file);
  Variant variant; Vector<Genotype> genotypes;
  while(reader.nextVariant(variant,genotypes)) {
    const int pos=variant.getPos();
    if(pos<chromBegin || pos>chromEnd) continue;
    file.print(variant.getText());
    for(Vector<int>::iterator cur=wantIndices.begin(), end=wantIndices.end() ;
	cur!=end ; ++cur)
      file.print("\t"+genotypes[*cur].getText());
    file.print("\n");
  }
  reader.close();
  delete &file;

  return 0;
}



void Application::getSampleIndices(const Vector<String> &wantIDs,
				   const Vector<String> &sampleIDs,
				   Vector<int> &indices)
{
  Set<String> want;
  for(Vector<String>::const_iterator cur=wantIDs.begin(), end=wantIDs.end() ;
      cur!=end ; ++cur) want.insert(*cur);
  int index=0;
  for(Vector<String>::const_iterator cur=sampleIDs.begin(), end=sampleIDs.end() ;
      cur!=end ; ++cur) {
    if(want.isMember(*cur)) indices.push_back(index);
    ++index;
  }
}



void Application::emitChromLine(const Vector<String> &ids,File &file)
{
  file.print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
  for(Vector<String>::const_iterator cur=ids.begin(), end=ids.end() ; cur!=end ; 
      ++cur)
    file.print("\t"+*cur);
  file.print("\n");
}



void Application::emitHeaderLines(const Vector<String> &lines,File &f)
{
  for(Vector<String>::const_iterator cur=lines.begin(), end=lines.end() ;
      cur!=end ; ++cur)
    f.print(*cur+"\n");
}



