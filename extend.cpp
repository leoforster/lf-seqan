#include <fstream>
#include <iostream>
#include <time.h>

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
//handling failed seeds?

int main(int argc, char const ** argv)
{
  int i = 0, j = 0, seedcount = 0;
  struct timespec start, finish;
  double elapsed;
  StringSet<CharString> ids;
  StringSet<CharString> seqs, quer; //DnaString ?
  Score<int> scoring(2, -1, -2);
  
  //assert(argc >= 4);
  int filecount = atoi(argv[1]);
  if (filecount == 0) 
  {
    std::cerr << "couldn't convert \"" << argv[1] << "\"\n";
    return 1;
  }
  else 
  {
    filecount == 1 ? j = 3 : j = 4;
  }

  //check files exist
  for (i = 2; i < j; i++)
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
  
  const char* seedFile = argv[2];
  const char* seqFile = argv[3];
  const char* qerFile = NULL;
  if (filecount == 2)
    qerFile = argv[4];

  //parse sequences
  SeqFileIn seqFileIn(seqFile);
  readRecords(ids, seqs, seqFileIn);
  if (filecount > 1)
  {
    SeqFileIn seqFileIn(qerFile);
    readRecords(ids, quer, seqFileIn);
  }
  
  //typedef Iterator<StringSet<CharString>, Standard>::Type TIterator;
  //for (TIterator it = begin(quer, Standard()); it != end(quer, Standard()); ++it)
    //std::cout << position(it, seqs) << " " << *it << std::endl;
  
  //parse and extend seeds
  std::ifstream infile(seedFile);
  std::string fail, strand, sid, qid, spos, qpos, len;
  
  clock_gettime(CLOCK_MONOTONIC, &start);
  
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
    
    //Seed<Simple> seed(stoi(spos), stoi(qpos), stoi(len)); 
    //extendSeed(seed, s, q, EXTEND_BOTH, MatchExtend());
    //extendSeed(seed, s, q, EXTEND_BOTH, scoring, 2, GappedXDrop());
    //extendSeed(seed, s, q, EXTEND_BOTH, scoring, 6, UnGappedXDrop());
    String<char> s, q;
    s = seqs[stoi(sid)];
    filecount == 1 ? q = seqs[stoi(qid)] : q = quer[stoi(qid)];
    toUpper(s);
    toUpper(q);
    
    //unsure
    if (strand == "P")
      reverseComplement(q);

    unsigned int sstart = stoi(spos);
    unsigned int send = stoi(spos) + stoi(len);
    unsigned int qstart, qend;
    //if (strand == "P")
    //{
      //qstart = beginPosition(q) + (endPosition(q) - stoi(qpos));
      //qend = beginPosition(q) + (endPosition(q) - stoi(qpos) + stoi(len));
    //} else
    //{
      qstart = stoi(qpos);
      qend = stoi(qpos) + stoi(len);
    //}
    
    Align<Infix<CharString const>::Type> align;
    AlignmentStats stats;
    resize(rows(align), 2);
    
    assignSource(row(align, 0), infix(s, sstart, send));
    assignSource(row(align, 1), infix(q, qstart, qend));
    Tuple<unsigned, 4> positions = { {sstart, qstart, send, qend} };
    extendAlignment(align, s, q, positions, EXTEND_BOTH, 2, scoring);
    
    //globalAlignment(align, scoring);
    computeAlignmentStats(stats, align, scoring);

    int new_slen = clippedEndPosition(row(align, 0)) - clippedBeginPosition(row(align, 0));
    int new_qlen = clippedEndPosition(row(align, 1)) - clippedBeginPosition(row(align, 1));
    int new_spos = clippedBeginPosition(row(align, 0));
    //int new_qpos = beginPosition(q) + (endPosition(q) - clippedBeginPosition(row(align, 1)));
    int new_qpos = clippedBeginPosition(row(align, 1));
    int score = stats.alignmentScore;
    int edist = stats.numMismatches + stats.numInsertions + stats.numDeletions;
    float ident = stats.alignmentIdentity; //float vs int? decimal places?

    //std::cout << "# Fields: s.len, s.seqnum, s.start, strand, q.len, q.seqnum, "
                  //" q.start, score, editdist, identity, seedlen, s.seedstart, "
                  //"q.seedstart" << std::endl;

    std::cout << new_slen << " " //lengths are off due to gaps?
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
              << spos << " "
              << qpos << "\n";
    //std::cout << align << std::endl;
  }
  
  clock_gettime(CLOCK_MONOTONIC, &finish);
  
  elapsed = (finish.tv_sec - start.tv_sec);
  elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
  std::cout << "# ... xdrop extension of " << seedcount << " seeds in " 
            << elapsed << " seconds." << std::endl;

  return 0;
}























