#include <iostream>

#include <seqan/align.h>
#include <seqan/align_extend.h>
#include <seqan/sequence.h>

using namespace seqan;

int main()
{
  Align<Infix<CharString const>::Type> align;
  resize(rows(align), 2);
  int score = 0;
  
  String<char> a = "SEEDabXcdXefXXX";
  String<char> b = "SEEDabYcdefYYYY";
  
  int left_a = 0; 
  int left_b = 0;
  int length = 4;
  //int right_a = left_a + length;
  //int right_b = left_b + length;
  
  Seed<Simple> seed1(0, 0, 4);          //left=0; length=4
  Score<> scoring(1, -1, -1);

  assignSource(row(align, 0), infix(a, left_a, left_a + length));
  assignSource(row(align, 1), infix(b, left_b, left_b + length));
  score = globalAlignment(align, scoring);
  
  std::cout << "Initial alignment (score == " << score << ")\n\n"
            << align;
  
  //extendSeed(seed, a, b, 1, MatchExtend());
  //extendSeed(seed1, a, b, EXTEND_BOTH, MatchExtend());
  //extendSeed(seed1, a, b, EXTEND_RIGHT, MatchExtend());
  //extendSeed(seed1, a, b, EXTEND_BOTH, scoring, 2, UnGappedXDrop());
  extendSeed(seed1, a, b, EXTEND_BOTH, scoring, 2, GappedXDrop());
  
  length = endPositionH(seed1) - left_a;
  assignSource(row(align, 0), infix(a, left_a, left_a + length));
  assignSource(row(align, 1), infix(b, left_b, left_b + length));
  score = globalAlignment(align, scoring);

  std::cout << "Resulting alignment (score == " << score << ")\n\n"
            << align;

  //Seed<Simple> seed2(0, 0, 4);          //left=0; length=4
  //Score<> scoring(1, -1, -1);
  //extendSeed(seed2, a, b, EXTEND_BOTH, scoring, 2, UnGappedXDrop());
  //std::cout << "endPositionH(seed2) = " << endPositionH(seed2) << std::endl;
  //std::cout << "endPositionV(seed2) = " << endPositionV(seed2) << std::endl;

  //Seed<Simple> seed3(0, 0, 4);          //left=0; length=4
  //extendSeed(seed3, a, b, EXTEND_BOTH, scoring, 2, GappedXDrop());
  //std::cout << "endPositionH(seed3) = " << endPositionH(seed3) << std::endl;
  //std::cout << "endPositionV(seed3) = " << endPositionV(seed3) << std::endl;

  return 0;
}
