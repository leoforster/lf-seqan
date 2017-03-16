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

int main(int argc, char const ** argv)
{
  int i, j;
  
  if (argc <= 3)
    return 1;

  int seqcount = atoi(argv[1]);
  int seedcount = atoi(argv[2]);
  
  assert(argc >= 2 + seqcount + seedcount * 5);
  
  StringSet<String<char> > seqs;
  String<char> seq;
  for (i = 3; i < 3 + seqcount; ++i)
  {
    seq = argv[i];
    appendValue(seqs, seq);
  }
  
  //Score<int, Simple> scoringScheme(1, -1, -1);
  Score<int, Simple> scoring(0, 1, 1);
  for (j = i; j < i + seedcount * 5; j = j + 5)
  {
    int sid = atoi(argv[j]);
    int qid = atoi(argv[j + 1]);
    int spos = atoi(argv[j + 2]);
    int qpos = atoi(argv[j + 3]);
    int len = atoi(argv[j + 4]);
    
    String<char> s = seqs[sid];
    String<char> q = seqs[qid];
    
    //std::cout << s << std::endl;
    //std::cout << q << std::endl;
    
    Seed<Simple> seed(spos, qpos, len);
    //extendSeed(seed, s, q, EXTEND_BOTH, MatchExtend());
    extendSeed(seed, s, q, EXTEND_BOTH, scoring, 6, GappedXDrop());
    
    //Align<Infix<CharString const>::Type> align;
    Align<CharString> align;
    AlignmentStats stats;
    resize(rows(align), 2);
    
    assignSource(row(align, 0), infix(s, beginPositionH(seed), endPositionH(seed)));
    assignSource(row(align, 1), infix(q, beginPositionV(seed), endPositionV(seed)));
    
    globalAlignment(align, scoring);
    computeAlignmentStats(stats, align, scoring);

    int slen = endPositionH(seed) - spos;
    int qlen = endPositionV(seed) - qpos;
    spos = beginPositionH(seed);
    qpos = beginPositionV(seed);
    char strand = '?'; //
    int score = stats.alignmentScore;
    int edist = stats.numPositiveScores; //stats.numNegativeScores?
    int ident = stats.alignmentIdentity; //float vs int? decimal places?

    //std::cout << "len " << slen << " "
              //<< "sid " << sid << " "
              //<< "spos " << spos << " "
              //<< "strand" << strand << " "
              //<< "qlen " << qlen << " "
              //<< "qid " << qid << " "
              //<< "qpos " << qpos << " "
              //<< "score " << score << " "
              //<< "edist " << edist << " "
              //<< "ident " << ident << "\n";
    std::cout << slen << " "
              << sid << " "
              << spos << " "
              << strand << " "
              << qlen << " "
              << qid << " "
              << qpos << " "
              << score << " "
              << edist << " "
              << ident << "\n";
  }
  return 0;
}























