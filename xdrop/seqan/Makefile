CXX?=g++
SYSTEM?=$(shell uname -s)

SEQANPATH?=${PACKAGES}/seqan

SEQANINCLUDE+=-I${SEQANPATH}/include

# do not include cairo including headers - needed by libgenometools
MY_CXXFLAGS=-Wall -Werror -Wpedantic -O3 -fno-strict-aliasing -std=c++14

MY_LDFLAGS=-ldl -lz -lbz2 ${RT} -lpthread ${FOLLY} ${LDFLAGS}
EXECS=example.x \
	extend.x \
	test.x

# the following does not work due to change of interface: blastread.x blastwrite.x

all: ${EXECS}

%.o:%.cpp
	${CXX} ${MY_CXXFLAGS} ${SEQANINCLUDE} -c -o $@ $^

%.x:%.cpp
	${CXX} ${MY_CXXFLAGS} ${SEQANINCLUDE} -o $@ $^ ${MY_LDFLAGS}

# for Mac need to install genometools.dylib at /usr/local/lib/

clean:
	${RM} *.[ox]
