/****************************************************************
 phd.C
 Copyright (C)2020 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/GffReader.H"
#include "BOOM/Regex.H"
#include "SamReader.H"
using namespace std;
using namespace BOOM;

class Application {
  Regex pseudogeneRegex;
  Regex chrRegex;
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
  const String vcfFile=cmd.arg(0);
  const String gffFile=cmd.arg(1);
  const String samFile=cmd.arg(2);

  GffReader gff(gffFile);
  while(true) {
    GffFeature *feature=gff.nextFeature();
    if(feature==NULL) break;
    if(feature->getFeatureType()!="exon" ||
       pseudogeneRegex.search(feature->lookupExtra("gene_type")) ||
       feature->lookupExtra("transcript_support_level")!="1")
      { delete feature; continue; }
    const String &transcriptSupport=
      feature->lookupExtra("transcript_support_level");
    cout<<"SUPPORT="<<transcriptSupport<<endl;
    cout<<"FEATURE="<<feature->getFeatureType()<<endl;
    String substrate=feature->getSubstrate();
    substrate=chrRegex.substitute(substrate,"");
    cout<<"SUBSTRATE="<<substrate<<endl;
    cout<<"GENETYPE="<<feature->lookupExtra("gene_type")<<endl;
    cout<<"gene_id="<<feature->lookupExtra("gene_id")<<endl;
    cout<<"=============================================="<<endl;
    delete feature;
  }

  return 0;

  // Process sam file
  SamReader sam(samFile);
  while(true) {
    SamRecord *rec=sam.nextRecord();
    if(!rec) break;
    if(rec->flag_unmapped()) { delete rec; continue; }
    /*cout<<rec->getID()<<"\t"<<rec->getRefName()<<"\t"<<rec->getRefPos()
	<<"\tunmapped="<<rec->flag_unmapped()<<" rev="
	<<rec->flag_revComp()<<"\t"<<rec->getCigar()<<endl;*/
    delete rec;
  }

  return 0;
}

