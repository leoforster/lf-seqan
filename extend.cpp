#include <fstream>
#include <iostream>

#include <seqan/basic.h>
#include <seqan/align.h>
#include <seqan/align_extend.h>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/modifier.h>

using namespace seqan;

//todo: seed chaining 
//http://docs.seqan.de/seqan/2.0.0/demos/seeds/seeds_chaining.cpp
//http://seqan.readthedocs.io/en/master/Tutorial/Algorithms/SeedExtension.html
//Merge() ? see http://docs.seqan.de/seqan/2.1.0/class_SeedSet.html

//todo: lengths of extended seeds in each sequence are identical
//how to count gaps to get correct length?

void alignment_based_extension(int sid, int qid, std::string strand, int len, 
                               String<char> s, String<char> q, Score<int> scoring,
                               unsigned sstart, unsigned send, 
                               unsigned qstart, unsigned qend)
{
  Align<Infix<CharString const>::Type> align;
  AlignmentStats stats;
  resize(rows(align), 2);
  
  assignSource(row(align, 0), infix(s, sstart, send));
  assignSource(row(align, 1), infix(q, qstart, qend));
  Tuple<unsigned, 4> positions = { {sstart, qstart, send, qend} };
  extendAlignment(align, s, q, positions, EXTEND_BOTH, 2, scoring);
  
  //globalAlignment(align, scoring);
  computeAlignmentStats(stats, align, scoring);

  int new_slen = clippedEndPosition(row(align, 0)) - 
                 clippedBeginPosition(row(align, 0));
  int new_qlen = clippedEndPosition(row(align, 1)) - 
                 clippedBeginPosition(row(align, 1));
  int new_spos = clippedBeginPosition(row(align, 0));
  //int new_qpos = beginPosition(q) + 
                 //(endPosition(q) - clippedBeginPosition(row(align, 1)));
  int new_qpos = clippedBeginPosition(row(align, 1));
  int score = stats.alignmentScore;
  int edist = stats.numMismatches + stats.numInsertions + stats.numDeletions;
  float ident = stats.alignmentIdentity;

  std::cout << new_slen << " "
            << sid << " "
            << new_spos << " "
            << strand << " "
            << new_qlen << " "
            << qid << " "
            << new_qpos << " "
            << score << " "
            << edist << " "
            << std::setprecision(4)
            << ident << " "
            << len << " "
            << sstart << " "
            << qstart << "\n";
  //std::cout << align << std::endl;
}

void seed_based_extension(int sid, int qid, std::string strand, int len, 
                         String<char> s, String<char> q, Score<int> scoring,
                         unsigned sstart, unsigned send, 
                         unsigned qstart, unsigned qend)
{
  Seed<Simple> seed(sstart, qstart, len);
  extendSeed(seed, s, q, EXTEND_BOTH, scoring, 2, GappedXDrop());
  
  int new_slen = endPositionH(seed) - beginPositionH(seed);
  int new_spos = beginPositionH(seed);
  int new_qlen = endPositionV(seed) - beginPositionV(seed);
  int new_qpos = beginPositionV(seed);
  int score = 0;
  int edist = 0;
  float ident = 0.0;
  
  std::cout << new_slen << " "
            << sid << " "
            << new_spos << " "
            << strand << " "
            << new_qlen << " "
            << qid << " "
            << new_qpos << " "
            << score << " "
            << edist << " "
            << std::setprecision(4)
            << ident << " "
            << len << " "
            << sstart << " "
            << qstart << "\n";
}

int main(int argc, char const ** argv)
{
  int i = 0, j = 0, seedcount = 0;
  double starttime;
  StringSet<CharString> ids;
  StringSet<CharString> seqs, quer; //DnaString ?
  Score<int> scoring(2, -1, -2);
  
  //assert(argc >= 4);
  int filecount = atoi(argv[1]);
  int method = atoi(argv[2]);
  if (filecount == 0) 
  {
    std::cerr << "couldn't convert \"" << argv[1] << "\"\n";
    return 1;
  }
  else 
  {
    filecount == 1 ? j = 4 : j = 5;
  }

  //check files exist
  for (i = 3; i < j; i++)
  {
    if (FILE *file = std::fopen(argv[i], "r")) 
    {
      std::fclose(file);
    } else 
    {
      std::cerr << "couldn't open file " << argv[i] << "\n";
      return 1;
    }
  }
  
  const char* seedFile = argv[3];
  const char* seqFile = argv[4];
  const char* qerFile = NULL;
  if (filecount == 2)
    qerFile = argv[5];

  //parse sequences
  SeqFileIn seqFileIn(seqFile);
  readRecords(ids, seqs, seqFileIn);
  if (filecount > 1)
  {
    SeqFileIn seqFileIn(qerFile);
    readRecords(ids, quer, seqFileIn);
  }
  
  //parse and extend seeds
  std::ifstream infile(seedFile);
  std::string fail, strand, sid, qid, spos, qpos, len;
  
  starttime = sysTime();
  
  while (infile.good()) //last newline gives repeated alignment
  {
    seedcount++;
    std::getline(infile, fail, ',');
    std::getline(infile, strand, ',');
    std::getline(infile, sid, ',');     
    std::getline(infile, spos, ',');
    std::getline(infile, qid, ','); 
    std::getline(infile, qpos, ',');
    std::getline(infile, len);
    
    String<char> s, q;
    s = seqs[stoi(sid)];
    filecount == 1 ? q = seqs[stoi(qid)] : q = quer[stoi(qid)];
    toUpper(s);
    toUpper(q);
    
    //unsure
    if (strand == "P")
      reverseComplement(q);

    unsigned sstart = stoi(spos);
    unsigned send = stoi(spos) + stoi(len);
    unsigned qstart, qend;
    //if (strand == "P")
    //{
      //qstart = beginPosition(q) + (endPosition(q) - stoi(qpos));
      //qend = beginPosition(q) + (endPosition(q) - stoi(qpos) + stoi(len));
    //} else
    //{
      qstart = stoi(qpos);
      qend = stoi(qpos) + stoi(len);
    //}
    
    if (method == 0)
      alignment_based_extension(stoi(sid), stoi(qid), strand, stoi(len), s, q, 
                                scoring, sstart, send, qstart, qend);
    else
      seed_based_extension(stoi(sid), stoi(qid), strand, stoi(len), s, q, 
                           scoring, sstart, send, qstart, qend);
  }

  std::cout << "# ... xdrop extension of " << seedcount << " seeds in " 
            << sysTime() - starttime << " seconds." << std::endl;

  return 0;
}























