/****************************************************************
 count-combinations-in-fragments.C
 Copyright (C)2021 William H. Majoros (bmajoros@alumni.duke.edu)
 This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
 License (GPL) version 3, as described at www.opensource.org.
 ****************************************************************/
#include <iostream>
#include <fstream>
#include "BOOM/String.H"
#include "BOOM/CommandLine.H"
#include "BOOM/Pipe.H"
#include "BOOM/Random.H"
#include "BOOM/Interval.H"
#include "BOOM/Regex.H"
#include "BOOM/Map.H"
#include "BOOM/Array1D.H"
#include "BOOM/Array2D.H"
#include "BOOM/VcfReader.H"
#include "BOOM/VectorSorter.H"
using namespace std;
using namespace BOOM;

class VcfIndex {
public:
  VcfIndex(const Vector<VariantAndGenotypes> &,int binSize,int minCoord,
	   int maxCoord);
  void lookup(int coord,int fragLen,Vector<VariantAndGenotypes> &into);
protected:
  void initialize();
  const int minCoord, maxCoord, numBins, binSize;
  Array1D<int> pointers; // indices into the array of variants, or -1 if NULL
  const Vector<VariantAndGenotypes> &allVariants;
};

struct SimpleVariant {
  String ID;
  String chrom;
  int pos;
  Array1D<int> genotype;
  SimpleVariant(const String &ID,const String &chrom,int pos,const Array1D<int> &genotype)
  : ID(ID), chrom(chrom), pos(pos), genotype(genotype) {}
};

class Application {
public:
  Application();
  int main(int argc,char *argv[]);
  //bool parseVariant(const String &sampleID,const String &line,String &ID,int &pos,
  //		    Array1D<int> &genotype);
  //void getVariants(const String &sampleID,const Interval &,
  //		   Vector<SimpleVariant> &variants);
  void variantsInInterval(const Interval &,Vector<VariantAndGenotypes> &into);
  void variantsInInterval(const String &sampleID,const Interval &,
			  Vector<SimpleVariant> &into);
  Interval pickInterval(int intervalLen,int chromLen);
  Interval pickInterval(int intervalLen,int chrBegin,int chrEnd);
  int countHets(const String &sampleID,const Interval &,int whichHap0or1);
  int countHetsRead(const String &sampleID,const Interval &);
  void countHets(const String &sampleID,const Interval &fragment,
		 const Interval &read1,const Interval &read2,int &fragHets,
		 int &readHets,int whichHap0or1);
  void parseHeader(const String &filename);
  void addToPairCount(const SimpleVariant &v1,const SimpleVariant &v2,int whichHap0or1);
  void getComboStats(int &allFourHaps,int &totalPairs);
  void loadVCF(const String &filename);
  bool isSorted(const Vector<VariantAndGenotypes> &);
  void sort(Vector<VariantAndGenotypes> &);
private:
  Regex genotypeRegex;
  String chromName;
  String vcfFile;
  Map<String,int> sampleColumns; // column indices in VCF file for samples
  Map<String,Array2D<int>*> allelePairCounts; // key = position x position
  Vector<VariantAndGenotypes> allVariants;
  VcfIndex *vcfIndex;
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



VcfIndex::VcfIndex(const Vector<VariantAndGenotypes> &variants,int binSize,
		   int minCoord,int maxCoord)
  : allVariants(variants), binSize(binSize), minCoord(minCoord), maxCoord(maxCoord), 
    numBins(maxCoord/binSize+1)
{
  initialize();
}



void VcfIndex::lookup(int begin,int fragLen,Vector<VariantAndGenotypes> &into)
{
  if(begin<minCoord) throw "Coord too small in VcfIndex"; // ### DEBUGGING
  if(begin>maxCoord) throw "Coord too large in VcfIndex"; // ### DEBUGGING
  const int bin=(begin-minCoord)/binSize;
  const int index=pointers[bin];
  if(index<0) return;
  const int end=begin+fragLen, numVar=allVariants.size();
  int i=index;
  while(i<numVar && allVariants[i].variant.getPos()<begin) ++i;
  for(; i<numVar ; ++i) {
    const VariantAndGenotypes &v=allVariants[i];
    const int pos=v.variant.getPos();
    if(pos>=end) break;
    into.push_back(v);
  }
}



void VcfIndex::initialize()
{
  const int numVar=allVariants.size();
  int j=0;
  pointers.resize(numBins);
  for(int i=0 ; i<numBins ; ++i) {
    int begin=minCoord+i*binSize, end=minCoord+(i+1)*binSize;
    while(j<numVar && allVariants[j].variant.getPos()<begin) ++j;
    if(j<numVar && allVariants[j].variant.getPos()<end) pointers[i]=j++; 
    else pointers[i]=-1;
  }
  //cout<<pointers<<endl;
}



Application::Application()
  : genotypeRegex("(.+)\\\\|(.+)"), vcfIndex(NULL)
{
  // ctor
}



int Application::main(int argc,char *argv[])
{
  randomize();

  // Process command line
  CommandLine cmd(argc,argv,"");
  if(cmd.numArgs()!=11)
    throw String("count-combinations-in-fragments <indexed-vcf.gz> <chrom> <chrom-len> <fragment-size> <read-length> <num-fragments> <sample-IDs> <chr-begin> <chr-end> <coutfile> <VCF-index-bin-size>");
  vcfFile=cmd.arg(0);
  chromName=cmd.arg(1);
  const int chromLen=cmd.arg(2).asInt();
  const int fragmentLen=cmd.arg(3).asInt();
  const int readLen=cmd.arg(4).asInt();
  const int numFragments=cmd.arg(5).asInt();
  const String &sampleIDs=cmd.arg(6);
  const int chrBegin=cmd.arg(7).asInt();
  const int chrEnd=cmd.arg(8).asInt();
  const String outfile=cmd.arg(9);
  const int vcfBinSize=cmd.arg(10).asInt();

  // Parse sampleIDs
  Vector<String> samples;
  sampleIDs.getFields(samples,",");
  const int numSamples=samples.size();
  parseHeader(vcfFile);

  // Load VCF file
  cout<<"loading VCF file"<<endl;
  loadVCF(vcfFile);

  // Index the VCF file for fast lookup
  cout<<"Indexing VCF"<<endl;
  vcfIndex=new VcfIndex(allVariants,vcfBinSize,chrBegin,chrEnd);
  cout<<"done indexing"<<endl;

  // Simulate fragments
  ofstream os(outfile.c_str());
  for(int i=0 ; i<numFragments ; ++i) {
    // Generate fragment and read intervals
    String sample=samples[RandomNumber(numSamples)];
    //Interval fragment=pickInterval(fragmentLen,chromLen);
    Interval fragment=pickInterval(fragmentLen,chrBegin,chrEnd);
    Interval read1=Interval(fragment.getBegin(),
			    fragment.getBegin()+readLen);
    Interval read2=Interval(fragment.getEnd()-readLen,
			    fragment.getEnd());
    int whichHap=RandomNumber(2);

    // Count het sites in fragment and reads
    int fragHets, readHets;
    countHets(sample,fragment,read1,read2,fragHets,readHets,whichHap);
    os<<fragHets<<"\t"<<readHets<<endl;
  }

  // Collect statistics about combinations of alleles
  int allFourHaps, totalPairs;
  getComboStats(allFourHaps,totalPairs);
  cout<<"Observed all combinations in "<<allFourHaps<<" out of "<<totalPairs
      <<" pairs"<<endl;

  return 0;
}



bool Application::isSorted(const Vector<VariantAndGenotypes> &v)
{
  Vector<VariantAndGenotypes> &variants=const_cast<Vector<VariantAndGenotypes>&>(v);
  VariantAndGenotypesComparator cmp;
  const int n=variants.size();
  for(int i=0 ; i+1<n ; ++i)
    if(cmp.greater(variants[i],variants[i+1])) return false;
  return true;
}



void Application::sort(Vector<VariantAndGenotypes> &variants)
{
  cerr<<"sorting"<<endl;
  VariantAndGenotypesComparator cmp;
  VectorSorter<VariantAndGenotypes> sorter(variants,cmp);
  sorter.sortAscendInPlace();
  cerr<<"done sorting"<<endl;
}



void Application::variantsInInterval(const String &sampleID,const Interval &interval,
				     Vector<SimpleVariant> &into)
{
  Vector<VariantAndGenotypes> vars;
  vcfIndex->lookup(interval.getBegin(),interval.getLength(),vars);
  for(Vector<VariantAndGenotypes>::iterator cur=vars.begin(), end=vars.end() ; 
      cur!=end ; ++cur) {
    const VariantAndGenotypes &vg=*cur;
    const int index=sampleColumns[sampleID];
    Genotype g=vg.genotypes[index];
    if(!g.isHet()) continue;
    Array1D<int> geno(2); geno[0]=g[0]; geno[1]=g[1];
    SimpleVariant v(vg.variant.getID(),vg.variant.getChr(),vg.variant.getPos(),
		    geno);
    into.push_back(v);
  }
}



/*void Application::variantsInInterval(const String &sampleID,const Interval &interval,
				     Vector<SimpleVariant> &into)
{
  Vector<VariantAndGenotypes> vars;
  variantsInInterval(interval,vars);
  for(Vector<VariantAndGenotypes>::iterator cur=vars.begin(), end=vars.end() ; 
      cur!=end ; ++cur) {
    const VariantAndGenotypes &vg=*cur;
    const int index=sampleColumns[sampleID];
    Genotype g=vg.genotypes[index];
    if(!g.isHet()) continue;
    Array1D<int> geno(2); geno[0]=g[0]; geno[1]=g[1];
    SimpleVariant v(vg.variant.getID(),vg.variant.getChr(),vg.variant.getPos(),
		    geno);
    into.push_back(v);
  }
  }*/



void Application::variantsInInterval(const Interval &interval,
				     Vector<VariantAndGenotypes> &into)
{
  const int n=allVariants.size();
  int begin=0, end=n, want=interval.getBegin();
  while(begin<end) {
    const int mid=(begin+end)/2;
    const int midVal=allVariants[mid].variant.getPos();
    if(mid<midVal) end=mid;
    else begin=mid;
  }
  if(want<allVariants[begin].variant.getPos()) throw "XXX1";
  int i;
  for(i=begin ; i<n && allVariants[i].variant.getPos()<want ; ++i);
  end=interval.getEnd();
  for(; i<n && allVariants[i].variant.getPos()<end ; ++i) 
    into.push_back(allVariants[i]);
}



void Application::loadVCF(const String &filename)
{
  VcfReader reader(filename);
  VariantAndGenotypes vg;
  Vector<Genotype> &G=vg.genotypes;
  while(reader.nextVariant(vg)) {
    bool skip=false, seen0=false, seen1=false;
    for(Vector<Genotype>::iterator cur=G.begin(), end=G.end() ; cur!=end ; ++cur) {
      Genotype &gt=*cur;
      const int N=gt.numAlleles();
      if(N!=2) { skip=true; break; }
      int gt0=gt[0], gt1=gt[1];
      if(gt0>1 || gt1>1 ) { skip=true; break; }
      if(gt0==0 || gt1==0) seen0=true;
      if(gt0==1 || gt1==1) seen1=true;
    }
    if(skip || !seen0 || !seen1) continue;
    allVariants.push_back(vg);
  }
  if(!isSorted(allVariants)) sort(allVariants);
}



void Application::getComboStats(int &allFourHaps,int &totalPairs)
{
  totalPairs=0;
  allFourHaps=0;
  for(Map<String,Array2D<int>*>::iterator cur=allelePairCounts.begin(),
	end=allelePairCounts.end() ; cur!=end ; ++cur) {
    pair<String,Array2D<int>*> p=*cur;
    const Array2D<int> &m=*p.second;
    //cout<<p.first<<"\t"<<m[0][0]<<"\t"<<m[0][1]<<"\t"<<m[1][0]<<"\t"<<m[1][1]<<endl;
    int nonzero=0;
    for(int i=0 ; i<2 ; ++i) 
      for(int j=0 ; j<2 ; ++j)
	if(m[i][j]>0) ++nonzero;
    if(nonzero==4) ++allFourHaps;
    ++totalPairs;
  }
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
      sampleColumns[sample]=i-9;
      //cout<<sample<<" assigned col "<<i<<endl;
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



Interval Application::pickInterval(int len,int chrBegin,int chrEnd)
{
  const int lastPos=chrEnd-len;
  if(lastPos<chrBegin) 
    throw String("chromosome region is smaller than fragment length");
  const int begin=chrBegin+RandomNumber(lastPos-chrBegin);
  return Interval(begin,begin+len);
}




/*void Application::getVariants(const String &sampleID,const Interval &interval,
			      Vector<SimpleVariant> &variants)
{
  const String cmd=String("tabix ")+vcfFile+" "+chromName+":"
    +String(interval.getBegin())+"-"+String(interval.getEnd());
  Pipe pipe(cmd,"r");
  while(!pipe.eof()) {
    const String line=pipe.getline();
    if(line.length()>0 && line[0]=='#') continue;
    String ID; int pos; Array1D<int> genotype(2);
    if(!parseVariant(sampleID,line,ID,pos,genotype)) continue;
    SimpleVariant v(ID,chromName,pos,genotype);
    variants.push_back(v);
  }
  pipe.close();
  }*/



void Application::addToPairCount(const SimpleVariant &v1,const SimpleVariant &v2,int whichHap)
{
  String key=String(v1.pos)+" "+String(v2.pos);
  if(!allelePairCounts.isDefined(key)) {
    Array2D<int> *a=new Array2D<int>(2,2);
    allelePairCounts[key]=a;
    a->setAllTo(0);
  }
  Array2D<int> &a=*allelePairCounts[key];
  if(v1.genotype[whichHap]>1) cout<<"v1 "<<whichHap<<" "<<v1.genotype[whichHap]<<endl;
  if(v2.genotype[whichHap]>1) cout<<"v2 "<<whichHap<<" "<<v2.genotype[whichHap]<<endl;
  ++a[v1.genotype[whichHap]][v2.genotype[whichHap]];
}



/*bool Application::parseVariant(const String &sampleID,const String &line,
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
  const String gt1=genotypeRegex[1], gt2=genotypeRegex[2];
  if(gt1.length()!=1 || gt2.length()!=1) return false;
  const bool het=gt1!=gt2;
  if(!het) return false;
  genotype[0]=gt1.asInt(); genotype[1]=gt2.asInt();
  if(genotype[0]>1 || genotype[1]>1) return false; // bi-allelic variants only

  pos=fields[1].asInt()-1; // Convert 1-based coord to 0-based
  ID=fields[2]; 
  return true;
  }*/



int Application::countHets(const String &sampleID,const Interval &interval,
			   int whichHap)
{
  Vector<SimpleVariant> variants;
  //getVariants(sampleID,interval,variants);
  variantsInInterval(sampleID,interval,variants);
  const int n=variants.size();
  for(int i=0 ; i<n-1 ; ++i)
    addToPairCount(variants[i],variants[i+1],whichHap);
  return variants.size();
}



int Application::countHetsRead(const String &sampleID,const Interval &interval)
{
  Vector<SimpleVariant> variants;
  //getVariants(sampleID,interval,variants);
  variantsInInterval(sampleID,interval,variants);
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



