
//Listing 1: Global alignment computation.
//----------------------------------------

//#include <iostream>
//#include <seqan/stream.h>
//#include <seqan/align.h>

//using namespace seqan;
//int main()
//{
  //StringSet<DnaString> stringSet;

  //appendValue(stringSet, "AGTTTAATCA");
  //appendValue(stringSet, "AGTATACGA");

  ////doesnt compile
  ////Align<DnaString, AnchorGaps> align(stringSet);
  ////compiles
  //Align<DnaString> align(stringSet);

  //int score = globalAlignment(align, EditDistanceScore());
  //std::cout << "Score :" << score << std::endl;
  //std::cout << align << std::endl;

  //return 0;
//}



//Listing 2: Searching with FM-indices
//----------------------------------------

//#include <iostream>
//#include <seqan/find.h>
//#include <seqan/index.h>

//using namespace seqan;

//int main()
//{
  //CharString hstck = "I spy with my little eye something that is yellow";
  //Index<CharString, FMIndex<> > index(hstck);
  //Finder<Index<CharString, FMIndex<> > > finder(hstck);

  //while (find(finder, "y"))
      //std::cout << "Hit at position: " << position(finder) << std::endl;

  //clear(finder);

  //while (find(finder, "t"))
      //std::cout << "Hit at position: " << position(finder) << std::endl;
//}
//doesnt compile (but probably should): 
//example.cpp:(.text+0x25): undefined reference to `aio_error'
//example.cpp:(.text+0x9e): undefined reference to `aio_suspend'
//example.cpp:(.text+0xb4): undefined reference to `aio_return'
//In function `seqan::BufferHandler<seqan::Pool<seqan::Triple<unsigned long, ...



//Listing 3: Reading and writing BAM files
//----------------------------------------
//#include <iostream>
//#include <seqan/stream.h>
//#include <seqan/bam_io.h>

//using namespace seqan;

//int main()
//{
  //BamFileIn bamFileIn;
  ////doesnt compile
  ////if (!open(bamFileIn, example.bam))
  ////compiles
  //if (!open(bamFileIn, "example.bam"))
  //{
    //std::cerr << "Can’t open the file." << std::endl;
    //return 1;
  //}
  
  //// Open output stream to stdcout .
  //BamFileOut bamFileOut(bamFileIn);
  //open(bamFileOut, std::cout, Sam());
  
  //// Copy header.
  //BamHeader header;
  //readHeader(header, bamFileIn);
  //writeHeader(bamFileOut, header);

  //// Copy records.
  //BamAlignmentRecord record;
  //while (!atEnd(bamFileIn))
  //{
      //readRecord(record, bamFileIn);
      //writeRecord(bamFileOut, record);
  //}

  //return 0;
//}














