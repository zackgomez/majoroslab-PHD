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
using namespace std;
using namespace BOOM;

class Application {
  Regex gzRegex;
  Map<String,File*> fileHandles;
  void closeOutputs();
public:
  Application();
  int main(int argc,char *argv[]);
  void process(File &,const String &outDir);
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
  : gzRegex("*\\.gz$")
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=2)
    throw String("split-vcf <in.vcf[.gz]> <out-dir>");
  const String infile=cmd.arg(0);
  const String outDir=cmd.arg(1);

  // Process input VCF file
  File &file=gzRegex.match(infile) ?
    *new File(infile) : *new GunzipPipe(infile);
  process(file,outDir);

  // Close output files
  closeOutputs();
  
  delete &file;
  return 0;
}



void Application::process(File &file,const String &outDir)
{
  while(!file.eof()) {
    const String line=file.getline();
    if(line.length()==0 || line[0]=='#') continue;
    Vector<String> fields;
    line.getFields(fields,"\t");
    if(fields.size()<9) continue;
    const String chr=fields[0];
    const int pos=fields[1].asInt();
    const String id=fields[2];
    const String ref=fields[3];
    const String alt=fields[4];
    const String pass=fields[6];
    if(ref.length()>1 || alt.length()>1) continue;
    if(pass!="PASS") continue;
    const String outfile=outDir+"/"+chr+".vcf.gz";
    if(!fileHandles.isDefined(outfile))
      fileHandles[outfile]=new GzipPipe(outfile);
    File &fh=*fileHandles[outfile];
    line=chr+'\t'+pos+'\t'+id+'\t'+ref+'\t'+alt+'\t';
    fh.print(line+"\n");
  }
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

