#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

using namespace seqan;

int main(int argc, char ** argv)
{
  typedef StringSet<IupacString> ReferenceSet;
  const char *program = argv[0];
  const char *inputfile = argv[1];

  ReferenceSet referenceSet;
  SeqFileIn seqFileIn;
  CharString header;
  IupacString seq;

  if (!open(seqFileIn, inputfile))
  {
    std::cerr << "ERROR: " << program << ": Could not open the file "
              << inputfile << std::endl;
    return -1;
  }
  while (!atEnd(seqFileIn))
  {
    try
    {
      readRecord(header, seq, seqFileIn);
    }
    catch (Exception const & e)
    {
      std::cerr << "ERROR: " << e.what() << std::endl;
      return -1;
    }
    appendValue(referenceSet, seq);
  }
  return 0;
}
