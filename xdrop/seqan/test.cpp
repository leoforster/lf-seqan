#include <fstream>
#include <iostream>

#include <seqan/stream.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

using namespace seqan;

int main(int argc, char ** argv)
{
  //const char *program = argv[0];
  const char *inputfile = argv[1];
  
  std::cout << inputfile << std::endl;
  
  std::fstream in(inputfile);
  //RecordReader<std::fstream, SinglePass<> > reader(in);
  //CharString id;
  Dna5String seq;
  
  //if (!open(seqFileIn, inputfile))
  //{
    //std::cerr << "ERROR: " << program << ": Could not open the file "
              //<< inputfile << std::endl;
    //return -1;
  //}
  if(in.is_open())
  {
    while (std::getline(in, seq))
    {
      std::cout << seq << '\n';
    }
    in.close();
  }
  return 0;
}
