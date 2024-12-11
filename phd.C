/****************************************************************
 phd.C : Piled Higher & Deeper

 Copyright (C)2022 William H. Majoros (bmajoros@alumni.duke.edu)
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
#include "BOOM/IlluminaQual.H"
#include "BOOM/SumLogProbs.H"
#include "BOOM/ConfigFile.H"
#include "SamReader.H"
// #include "SamTabix.H"
#include "VariantGraph.H"
#include "VariantInRead.H"
#include "VcfStream.H"
using namespace std;
using namespace BOOM;

bool DEBUG = false;

/****************************************************************
 TO DO:
 * Some reads extend past annotated transcript; variants outside the
   annotated transcript are discarded, but they could be used
 ****************************************************************/

struct GeneRegion
{
  String geneID;
  String substrate;
  Interval interval;
};

class Application
{
  VcfStream *vcfStream;
  // String tabix;
  float MIN_PROB_CORRECT;
  bool DEDUPLICATE;
  IlluminaQual illumina;
  char MIN_QUAL; // This is the encoding charater, NOT the integer value!
  mutable Regex pseudogeneRegex;
  mutable Regex chrRegex;
  String vcfFile;
  int readsSeen, readsDiscarded, readsUnmapped, readsWrongChrom;
  int numConcordant, numNonzero, duplicatesRemoved;
  Set<int> seenPositions;
  Set<String> exonTypes;
  bool parseVariant(const String &line, String &ID, int &pos,
                    char &cRef, char &cAlt, int *genotype);
  void processSam(SamReader &, VariantGraph &,
                  int geneBegin, int geneEnd, const String &substrate);
  bool addEdges(const SamRecord *read, VariantGraph &graph);
  void installEdges(ReadVariants &, const String &readID,
                    const String &qualities, VariantGraph &);
  void findVariantsInRead(VariantGraph &, const SamRecord *,
                          CigarAlignment &, ReadVariants &,
                          const String &qualities);
  void processGraph(VariantGraph &, const String &geneID);

  std::vector<GeneRegion> readGeneRegions(GffReader &gff) const;

public:
  Application();
  int main(int argc, char *argv[]);
};

int main(int argc, char *argv[])
{
  try
  {
    Application app;
    return app.main(argc, argv);
  }
  catch (const char *p)
  {
    cerr << p << endl;
  }
  catch (const string &msg)
  {
    cerr << msg.c_str() << endl;
  }
  catch (const RootException &e)
  {
    cerr << e.getMessage() << endl;
  }
  catch (const exception &e)
  {
    cerr << "STL exception caught in main:\n"
         << e.what() << endl;
  }
  catch (...)
  {
    cerr << "Unknown exception caught in main" << endl;
  }
  return -1;
}

Application::Application()
    : pseudogeneRegex("pseudogene"), chrRegex("chr"),
      readsSeen(0), readsDiscarded(0), readsUnmapped(0), readsWrongChrom(0),
      numConcordant(0), numNonzero(0), duplicatesRemoved(0),
      vcfStream(NULL)
{
  // ctor

  exonTypes += "exon";
  exonTypes += "EXON";
  exonTypes += "utr";
  exonTypes += "UTR";
  exonTypes += "cds";
  exonTypes += "CDS";
  exonTypes += "five_prime_UTR";
  exonTypes += "three_prime_UTR";
  exonTypes += "single-exon";
  exonTypes += "initial-exon";
  exonTypes += "internal-exon";
  exonTypes += "final-exon";
}

int Application::main(int argc, char *argv[])
{
  // Process command line
  CommandLine cmd(argc, argv, "");
  if (cmd.numArgs() != 4)
    throw String("phd <config> <indexed.vcf.gz> <sorted.gff> <in.sam>");
  const String configFile = cmd.arg(0);
  vcfFile = cmd.arg(1);
  const String gffFile = cmd.arg(2);
  const String samFile = cmd.arg(3);

  // Load parameters from config file
  ConfigFile config(configFile);
  // tabix=config.lookupOrDie("tabix");
  MIN_PROB_CORRECT = config.getFloatOrDie("min-phasing-probability");
  DEDUPLICATE = config.getBoolOrDie("deduplicate");
  MIN_QUAL = illumina.phredToChar(config.getIntOrDie("min-base-quality"));

  // Open input files
  GffReader gff(gffFile);
  vcfStream = new VcfStream(vcfFile);
  SamReader samReader(samFile);

  std::vector<GeneRegion> geneRegions = readGeneRegions(gff);

  for (const auto &region : geneRegions)
  {
    VariantGraph graph;
    vcfStream->getVariants(region.interval, graph.getVariants());

    if (graph.size() == 0)
      continue;

    processSam(samReader, graph, region.interval.getBegin(), region.interval.getEnd(),
               region.substrate);
    processGraph(graph, region.geneID);
  }

  cerr << readsSeen << " reads seen, " << readsDiscarded << " discarded, "
       << readsUnmapped << " unmapped" << endl;
  cerr << numConcordant << " concordant edges out of " << numNonzero << " nonzero"
       << endl;
  cerr << duplicatesRemoved << " duplicate reads removed" << endl;
  return 0;
}

void Application::processSam(SamReader &sam, VariantGraph &graph, int geneBegin,
                             int geneEnd, const String &substrate)
{
  // XXX(zack): does this need to be static?  If so figure out a better way to do this.
  static SamRecord *buffer = NULL;
  // cout<<substrate<<":"<<geneBegin<<"-"<<geneEnd<<endl;
  int debug = 0;
  while (true)
  {
    SamRecord *rec;
    if (buffer)
    {
      rec = buffer;
      buffer = NULL;
    }
    else
    {
      rec = sam.nextRecord();
      if (rec)
        ++readsSeen;
    }
    if (!rec)
      break;
    if (rec->flag_unmapped())
    {
      cout << "unmapped" << endl;
      delete rec;
      ++readsUnmapped;
      continue;
    }
    const int refPos = rec->getRefPos();
    ++debug;

    // if(seenPositions.isMember(refPos)) { delete rec; continue; }
    // seenPositions+=refPos;
    if (DEDUPLICATE && rec->flag_PCRduplicate())
    {
      // cout<<"duplicate"<<endl;
      delete rec;
      ++duplicatesRemoved;
      continue;
    }
    String readSubstrate = rec->getRefName();
    readSubstrate = chrRegex.substitute(readSubstrate, "");
    if (readSubstrate < substrate)
    {
      delete rec;
      // cout << "wrong chrom" << endl;
      ++readsWrongChrom;
      continue;
    } // ### ???

    if (refPos + rec->getCigar().genomicSpan() < geneBegin)
    {
      // cout<<"discarding "<<rec->getID()<<" @ "<<refPos<<"+"<<rec->getSequence().getLength()<<" < "<<geneBegin<<endl;
      delete rec;
      ++readsDiscarded;
      continue;
    }

    if (refPos > geneEnd)
    {
      // cout << "after gene" << endl;
      buffer = rec;
      break;
    }
    const bool containsVariants = addEdges(rec, graph);
    // if(!containsVariants) cout<<"no variants"<<endl;
    // else cout<<"CONTAINS VARIANTS"<<endl;
    delete rec;
  }
  // cout<<debug<<" records retrieved by tabix"<<endl;
}

bool Application::parseVariant(const String &line, String &ID, int &pos,
                               char &cRef, char &cAlt, int *genotype)
{
  Vector<String> fields;
  line.getFields(fields, "\t");
  if (fields.size() < 9 || fields[6] != "PASS" || fields[8] != "GT")
    return false;
  String sGenotype = fields[9];

  // ### This is not fully generic, could be improved:
  if (sGenotype == "0|1")
  {
    genotype[0] = 0;
    genotype[1] = 1;
  }
  else if (sGenotype == "1|0")
  {
    genotype[0] = 1;
    genotype[1] = 0;
  }
  else
    return false;

  // XXX(zack): check out this indexing stuff
  pos = fields[1].asInt() - 1; // Convert 1-based coord to 0-based
  ID = fields[2];
  String ref = fields[3], alt = fields[4];
  if (ref.length() != 1 || alt.length() != 1)
    return false;
  cRef = ref[0];
  cAlt = alt[0];
  if (cRef != 'A' && cRef != 'C' && cRef != 'G' && cRef != 'T')
    return false;
  if (cAlt != 'A' && cAlt != 'C' && cAlt != 'G' && cAlt != 'T')
    return false;
  return true;
}

bool Application::addEdges(const SamRecord *read,
                           VariantGraph &graph)
{
  const CigarString &cigar = read->getCigar();
  CigarAlignment &alignment = *cigar.getAlignment();
  const int numNodes = graph.size();
  ReadVariants readVariants(read->getID());
  findVariantsInRead(graph, read, alignment, readVariants,
                     read->getQualityScores());
  bool ret = false;
  if (readVariants.size() > 0)
  {
    installEdges(readVariants, read->getID(), read->getQualityScores(),
                 graph);
    ret = true;
  }
  // delete &alignment;
  return ret; // = whether variants were found in the read
}

void Application::findVariantsInRead(VariantGraph &graph,
                                     const SamRecord *read,
                                     CigarAlignment &alignment,
                                     ReadVariants &variants,
                                     const String &qualities)
{
  const int L = alignment.length();
  const int offset = read->getRefPos();
  const String &seq = read->getSequence();
  const String &qual = read->getQualityScores();
  for (int readPos = 0; readPos < L; ++readPos)
  {
    const int refPos = alignment[readPos] + offset;
    if (refPos == CIGAR_UNDEFINED)
      continue;
    for (Vector<::Variant>::iterator cur = graph.begin(), end = graph.end();
         cur != end; ++cur)
    {
      ::Variant &v = *cur;
      if (v.getPos() == refPos)
      {
        if (qual[readPos] < MIN_QUAL)
          continue;
        const char c = seq[readPos];
        Allele allele;
        if (c == v.getRef())
          allele = REF;
        else if (c == v.getAlt())
          allele = ALT;
        else if (DEBUG)
        {
          const float pError = illumina.charToErrorProb(qual[readPos]);
          cout << read->getID() << " ALLELE MISMATCH P(error)=" << pError
               << "=" << illumina.charToPhred(qual[readPos]) << "="
               << qual[readPos] << " " << c << " NOT " << v.getRef()
               << " NOR " << v.getAlt() << " VAR=" << v.getID()
               << " READ POS=" << readPos << endl;
          continue;
        }
        else
        {
          continue;
        }
        const float p = 1 - illumina.charToErrorProb(qualities[readPos]);
        variants.push_back(VariantInRead(v, readPos, allele, p));
      }
    }
  }
}

void Application::installEdges(ReadVariants &read, const String &readID,
                               const String &qualities, VariantGraph &G)
{
  G.addRead(read);
  const int N = read.size();
  for (int i = 0; i < N - 1; ++i)
  {
    VariantInRead &thisVar = read[i], &nextVar = read[i + 1];
    thisVar.v->getEdges()[thisVar.allele][nextVar.allele] += 1;
    Vector<pair<float, float>> &cell =
        thisVar.v->getProbCorrect()[thisVar.allele][nextVar.allele];
    // const float p1=1-illumina.charToErrorProb(qualities[thisVar.pos]);
    // const float p2=1-illumina.charToErrorProb(qualities[nextVar.pos]);
    const float p1 = thisVar.probCorrect, p2 = nextVar.probCorrect;
    if (!isFinite(p1) || !isFinite(p2))
      cout << "PHRED: " << qualities[thisVar.pos] << "="
           << illumina.charToErrorProb(qualities[thisVar.pos])
           << " " << qualities[nextVar.pos]
           << " " << illumina.charToErrorProb(qualities[nextVar.pos]) << endl;
    cell.push_back(pair<float, float>(p1, p2));
  }
}

void Application::processGraph(VariantGraph &G, const String &geneID)
{
  // Count concordant/discordant edges
  const int N = G.size();
  int totalEdges = 0;
  for (int i = 0; i < N - 1; ++i)
  {
    const ::Variant &v = G[i];
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k)
        totalEdges += v.getEdges()[j][k];
    if (v.nonzero())
      ++numNonzero;
    if (v.concordant())
      ++numConcordant;
  }
  // if(totalEdges<1) return;

  // Phase the graph using the reads
  G.phase(illumina, MIN_PROB_CORRECT);

  // Get the connected components
  Vector<ConnectedComponent> components;
  G.getComponents(components, illumina, MIN_PROB_CORRECT);
  if (components.size() == 0)
    return;

  // cout<<"Graph has "<<G.size()<<" variants before consistency checking"<<endl;

  // Discard reads inconsistent with chosen phase
  Vector<ReadVariants *> &reads = G.getReads(), filtered;
  for (Vector<ReadVariants *>::iterator cur = reads.begin(), end = reads.end();
       cur != end; ++cur)
  {
    ReadVariants *read = *cur;
    if (read->size() == 0)
      continue;
    if (read->consistentWithPhase())
      filtered.push_back(read);
  }
  G.getReads() = filtered;

  // Assign maternal/paternal haplotypes to connected components
  G.phaseComponents(components);

  // Assign reads to haplotypes
  G.assignReads();

  // Print out variant IDs and read counts for all components
  int compNum = 1;
  for (Vector<ConnectedComponent>::iterator cur = components.begin(),
                                            end = components.end();
       cur != end; ++cur)
  {
    ConnectedComponent &comp = *cur;
    ::Variant &v = comp[0];
    const int ref = v.getCount(REF), alt = v.getCount(ALT);
    if (ref + alt == 0)
      continue;
    cout << geneID << "\t" << comp[0].getID() << "\t" << ref << "\t" << alt << "\t";
    for (int i = 0; i < comp.size(); ++i)
    {
      cout << comp[i].getID();
      if (i + 1 < comp.size())
        switch (comp[i].getPhase())
        {
        case UNPHASED:
          throw "unphased link in graph!";
        case IN_PHASE:
          cout << "-";
          break;
        case ANTI_PHASED:
          cout << "\\";
          break;
        default:
          throw String("unknown VariantPhase value ") + comp[i].getPhase();
        }
    }
    cout << endl;
    ++compNum;
  }
}

std::vector<GeneRegion> Application::readGeneRegions(GffReader &gff) const
{

  std::vector<GeneRegion> geneRegions;
  String currentGeneID;
  String currentSubstrate;
  Interval accumulatedRegion;
  while (GffFeature *feature = gff.nextFeature())
  {
    // Skip non-exon features
    if (!exonTypes.isMember(feature->getFeatureType()) || pseudogeneRegex.search(feature->lookupExtra("gene_type")))
    {
      delete feature;
      continue;
    }

    String geneID = feature->lookupExtra("gene_id");
    // first region
    if (currentGeneID == "")
    {
      currentGeneID = geneID;
      currentSubstrate = chrRegex.substitute(feature->getSubstrate(), "");
      accumulatedRegion = Interval(feature->getBegin(), feature->getEnd());
    }
    // accumulating region
    else if (geneID == currentGeneID)
    {
      accumulatedRegion.setBegin(std::min(accumulatedRegion.getBegin(), feature->getBegin()));
      accumulatedRegion.setEnd(std::max(accumulatedRegion.getEnd(), feature->getEnd()));
    }
    // new region
    else
    {
      assert(accumulatedRegion.length() > 0);
      geneRegions.push_back(GeneRegion{currentGeneID, currentSubstrate, accumulatedRegion});
      accumulatedRegion = Interval(feature->getBegin(), feature->getEnd());
      currentGeneID = geneID;
    }
  }
  if (accumulatedRegion.length() > 0)
  {
    geneRegions.push_back(GeneRegion{currentGeneID, currentSubstrate, accumulatedRegion});
  }

  // assert sorted and non-overlapping
  for (size_t i = 1; i < geneRegions.size(); ++i)
  {
    if (geneRegions[i].interval.getEnd() < geneRegions[i - 1].interval.getBegin())
    {
      throw RootException(String("GTF is not sorted: (") + geneRegions[i - 1].interval.getBegin() + "," + geneRegions[i - 1].interval.getEnd() + ") overlaps (" + geneRegions[i].interval.getBegin() + "," + geneRegions[i].interval.getEnd() + ")");
    }
  }

  return geneRegions;
}
