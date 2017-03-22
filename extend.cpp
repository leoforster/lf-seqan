#include <fstream>
#include <iostream>
#include <assert.h>  

#include <seqan/align.h>
#include <seqan/align_extend.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

using namespace seqan;

//todo: seed chaining 
//http://docs.seqan.de/seqan/2.0.0/demos/seeds/seeds_chaining.cpp
//http://seqan.readthedocs.io/en/master/Tutorial/Algorithms/SeedExtension.html
//Merge() ? see http://docs.seqan.de/seqan/2.1.0/class_SeedSet.html
//handling failed seeds?


int main(int argc, char const ** argv)
{
  int i;
  StringSet<CharString> ids;
  StringSet<CharString> seqs; //DnaString ?

  //check files exist
  for (i = 1; i < 3; i++)
  {
    if (FILE *file = std::fopen(argv[i], "r")) {
      std::fclose(file);
    } else {
      std::cerr << "couldn't open file " << argv[i] << "\n";
      return 1;
    }
  }
  
  const char* seqFile = argv[1];
  const char* seedFile = argv[2];
  //Score<int, Simple> scoring(0, 1, 1);
  Score<int, Simple> scoring(1, -1, -1);
  //Score<int> scoring(2, -1, -2);
  
  //parse sequences
  SeqFileIn seqFileIn(seqFile);
  readRecords(ids, seqs, seqFileIn);
  
  //typedef Iterator<StringSet<CharString>, Standard>::Type TIterator;
  //for (TIterator it = begin(seqs, Standard()); it != end(seqs, Standard()); ++it)
      //std::cout << position(it, seqs) << " " << *it << std::endl;
  
  //parse and extend seeds
  std::ifstream infile(seedFile);
  std::string fail, strand, sid, qid, spos, qpos, len;
  while (infile.good()) //last newline gives repeated alignment
  {
    std::getline(infile, fail, ',');
    std::getline(infile, strand, ',');
    std::getline(infile, sid, ',');     
    std::getline(infile, spos, ',');
    std::getline(infile, qid, ','); 
    std::getline(infile, qpos, ',');
    std::getline(infile, len);
    
    //Seed<Simple> seed(stoi(spos), stoi(qpos), stoi(len)); 
    //consider seedSet -- what about sid and qid though?
    String<char> s = seqs[stoi(sid)]; //CharString fails on compilation
    String<char> q = seqs[stoi(qid)];
    
    //cant use edit distance with gapped xdrop?
    //extendSeed(seed, s, q, EXTEND_BOTH, MatchExtend());
    //extendSeed(seed, s, q, EXTEND_BOTH, scoring, 2, GappedXDrop());
    //extendSeed(seed, s, q, EXTEND_BOTH, scoring, 6, UnGappedXDrop());
    unsigned int sstart = stoi(spos);
    unsigned int send = stoi(spos) + stoi(len);
    unsigned int qstart = stoi(qpos);
    unsigned int qend = stoi(qpos) + stoi(len);
    
    Align<Infix<CharString const>::Type> align;
    AlignmentStats stats;
    resize(rows(align), 2);
    
    assignSource(row(align, 0), infix(s, sstart, send));
    assignSource(row(align, 1), infix(q, qstart, qend));
    //assignSource(row(align, 0), infix(s, beginPositionH(seed), endPositionH(seed)));
    //assignSource(row(align, 1), infix(q, beginPositionV(seed), endPositionV(seed)));
    Tuple<unsigned, 4> positions = { {sstart, qstart, send, qend} };
    extendAlignment(align, s, q, positions, EXTEND_BOTH, 2, scoring);
    
    //globalAlignment(align, scoring);
    computeAlignmentStats(stats, align, scoring);

    //int slen = endPositionH(seed) - stoi(spos);
    //int qlen = endPositionV(seed) - stoi(qpos);
    //int spos_i = beginPositionH(seed);
    //int qpos_i = beginPositionV(seed);
    int slen = clippedEndPosition(row(align, 0)) - clippedBeginPosition(row(align, 0)) - 1;
    int qlen = clippedEndPosition(row(align, 1)) - clippedBeginPosition(row(align, 1)) - 1;
    int spos_i = clippedBeginPosition(row(align, 0));
    int qpos_i = clippedBeginPosition(row(align, 1));
    int score = stats.alignmentScore;
    int edist = stats.numPositiveScores; //stats.numNegativeScores?
    int ident = stats.alignmentIdentity; //float vs int? decimal places?

    std::cout << slen << " "
              << sid << " "
              << spos_i << " "
              << strand << " "
              << qlen << " "
              << qid << " "
              << qpos_i << " "
              << score << " "
              << edist << " "
              << ident << "\n";
  }

  return 0;
}























