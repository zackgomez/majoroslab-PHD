/****************************************************************
 Copyright (C)2021 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include "BOOM/String.H"
#include "BOOM/Regex.H"
#include "BOOM/CommandLine.H"
#include "BOOM/GffReader.H"
#include "BOOM/VectorSorter.H"
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
  : pseudogeneRegex("pseudogene"), chrRegex("chr")
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=1)
    throw String("sort-gff <in.gff>");
  const String filename=cmd.arg(0);

  // Open input file
  GffReader gff(filename);
  cerr<<"Loading GFF..."<<endl;
  Vector<GffGene> &genes=*gff.loadGenes();

  cerr<<"Sorting..."<<endl;
  GffGeneComparator cmp;
  VectorSorter<GffGene> sorter(genes,cmp);
  //sorter.sortAscendInPlace();
  Vector<int> indices;
  sorter.sortAscendByIndex(indices);

  cerr<<"Writing output..."<<endl;
  int n=genes.size();
  cerr<<n<<" genes"<<endl;
  for(int i=0 ; i<n ; ++i) {
    GffGene &gene=genes[indices[i]];
    if(i+1<n) {
      GffGene &nextGene=genes[indices[i+1]];
      Interval int1(gene.getBegin(),gene.getEnd());
      Interval int2(nextGene.getBegin(),nextGene.getEnd());
      if(int1.overlaps(int2) && int1.getLength()<int2.getLength())
	continue;
    }
    int numTrans=gene.numTranscripts();
    for(int j=0 ; j<numTrans ; ++j) {
      GffTranscript &trans=gene.getIthTranscript(j);
      trans.toGff(cout);
    }
  }
  cerr<<"Done."<<endl;
  return 0;
}
