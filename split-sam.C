/****************************************************************
 split-sam.C
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
#include "SamReader.H"
using namespace std;
using namespace BOOM;

class Application {
  Regex gzRegex;
  Map<String,File*> fileHandles;
  void closeOutputs();
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
  if(cmd.numArgs()!=2)
    throw String("split-sam <in.sam[.gz]> <out-dir>");
  const String infile=cmd.arg(0);
  const String outDir=cmd.arg(1);

  // Process input file
  SamReader reader(infile);
  while(true) {
    String text;
    SamRecord *rec=reader.nextSeqAndText(text);
    if(rec==NULL) break;
    const String outfile=outDir+"/"+rec->getRefName()+".sam.gz";
    File *fp=NULL;
    if(fileHandles.isDefined(outfile)) fp=fileHandles[outfile];
    else fp=fileHandles[outfile]=new GzipPipe(outfile);
    fp->print(text+"\n");
    delete rec;
  }
  reader.close();
  
  // Close files
  closeOutputs();

  return 0;
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



