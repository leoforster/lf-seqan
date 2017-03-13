CXX?=g++
SYSTEM?=$(shell uname -s)

ifeq ($(SYSTEM),Linux)
  RT=-lrt
  GTINSTALL?=/usr/local
  GTLIB=-L${GTINSTALL}/lib
else
  GTINSTALL=/usr/local
  # delete the following line if facebook vector implementation folly/FBvector
  # is not available
  # FOLLY=-DWITHFOLLY -lfolly
endif

SEQANPATH?=${PACKAGES}/seqan

# make sure that $GTINSTALL/lib is in LD_LIBRARY_PATH (for Linux) and
# DYLD_LIBRARY_PATH (for Mac)
GTINCLUDE=-I${GTINSTALL}/include/genometools
SEQANINCLUDE+=-I${SEQANPATH}/include

# do not include cairo including headers - needed by libgenometools
MY_CXXFLAGS=-Wall -Werror -Wpedantic -O3 -fno-strict-aliasing -std=c++14
MY_GT_CXXFLAGS=${MY_CXXFLAGS} -DWITHOUT_CAIRO -DSEQAN_HAS_ZLIB -DSEQAN_HAS_BZIP2
#-flto
# -ldl is needed by libgenometools
# -lrt -lpthread is needed by SeqAn

MY_LDFLAGS=-ldl -lz -lbz2 ${RT} -lpthread ${FOLLY} ${LDFLAGS}
EXECS=extend.x 

# the following does not work due to change of interface: blastread.x blastwrite.x

all: ${EXECS}

%.o:%.cpp
	${CXX} ${MY_CXXFLAGS} ${SEQANINCLUDE} -c -o $@ $^

%.x:%.cpp
	${CXX} ${MY_CXXFLAGS} ${SEQANINCLUDE} -o $@ $^ ${MY_LDFLAGS}

# for Mac need to install genometools.dylib at /usr/local/lib/

clean:
	${RM} *.[ox]
