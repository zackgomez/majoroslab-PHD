/****************************************************************
 split-vcf.C
 Copyright (C)2022 William H. Majoros (bmajoros@alumni.duke.edu).
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Pipe.H"
#include "BOOM/Regex.H"
#include "BOOM/Map.H"
#include "BOOM/VcfReader.H"
using namespace std;
using namespace BOOM;

class Application {
  Regex gzRegex;
  Map<String,File*> fileHandles;
  void closeOutputs();
  void getSampleIndices(const Vector<String> &wantIDs,
			const Vector<String> &sampleIDs,Vector<int> &indices);
  void emitHeaderLines(const Vector<String> &lines,File &);
  void emitChromLine(const Vector<String> &ids,File &);
  bool atLeastOneHet(const Vector<Genotype> &,const Vector<int> &indices);
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
  catch(const exception &e)
    {cerr << "STL exception caught in main:\n" << e.what() << endl;}
  catch(...) { cerr << "Unknown exception caught in main" << endl; }
  return -1;
}



Application::Application()
  : gzRegex("gz$")
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=3)
    throw String("split-vcf <in.vcf[.gz]> <sample-IDs> <out-dir>");
  const String infile=cmd.arg(0);
  const String samples=cmd.arg(1);
  const String outDir=cmd.arg(2);
  Vector<String> wantIDs;
  samples.getFields(wantIDs,",");

  // Process input VCF file
  VcfReader reader(infile);
  const Vector<String> &sampleIDs=reader.getSampleIDs();
  Vector<int> wantIndices;
  getSampleIndices(wantIDs,sampleIDs,wantIndices);
  if(wantIndices.size()==0) throw String("Can't find samples in VCF file");
  Variant variant; Vector<Genotype> genotypes;
  while(reader.nextVariant(variant,genotypes)) {
    if(variant.numAlleles()!=2) continue;
    if(!atLeastOneHet(genotypes,wantIndices)) continue;
    const String allele1=variant.getAllele(0), allele2=variant.getAllele(1);
    if(allele1.length()!=1 || allele2.length()!=1 || allele1==allele2)
      continue;
    const String outfile=outDir+"/"+variant.getChr()+".vcf.gz";
    File *fp=NULL;
    if(fileHandles.isDefined(outfile)) fp=fileHandles[outfile];
    else {
      fp=fileHandles[outfile]=new BGzipPipe(outfile);
      emitHeaderLines(reader.getHeaderLines(),*fp);
      emitChromLine(wantIDs,*fp); }
    fp->print(variant.getText());
    for(Vector<int>::iterator cur=wantIndices.begin(), end=wantIndices.end() ;
	cur!=end ; ++cur) fp->print("\t"+genotypes[*cur].getText());
    fp->print("\n");
  }
  reader.close();
  
  // Close files
  closeOutputs();

  return 0;
}



bool Application::atLeastOneHet(const Vector<Genotype> &genotypes,
				const Vector<int> &indices)
{
  for(Vector<int>::const_iterator cur=indices.begin(), end=indices.end() ;
      cur!=end ; ++cur)
    if(genotypes[*cur].isHet()) return true;
  return false;
}



void Application::closeOutputs()
{
  Set<File*> files;
  fileHandles.getValues(files);
  for(Set<File*>::iterator cur=files.begin(), end=files.end() ; cur!=end ;
      ++cur) {
    File *fp=*cur;
    fp->close();
    delete fp;
  }
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


