CXX=g++
RM=rm -rf

#include directories seperated by a space
INC=./utils/include/
INC_PARAMS=$(foreach d, $(INC), -I$d)

#lib directories seperated by a space
LIB=./utils/lib/
LIB_PARAMS=$(foreach d, $(LIB), -L$d)

#paramaters as suggested by seqan
CXXLAGS=$(INC_PARAMS) -std=c++14 -O3 -DNDEBUG -DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_BZIP2=1 -DSEQAN_DISABLE_VERSION_CHECK=YES -W -Wall -pedantic
#Parameters for debug purposes. Comment the above line and out the below one to compile in debug mode
#CXXLAGS=$(INC_PARAMS) -std=c++14 -g -O0 -DSEQAN_ENABLE_DEBUG=1 -DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_BZIP2=1 -DSEQAN_DISABLE_VERSION_CHECK=YES -W -Wall -pedantic

LDFLAGS=
LDLIBS=$(LIB_PARAMS) -lz -lbz2

ifeq ($(OS),Windows_NT)
    CXXLAGS += -D WIN32
    ifeq ($(PROCESSOR_ARCHITEW6432),AMD64)
        CXXLAGS=/W2 /wd4996 -D_CRT_SECURE_NO_WARNINGS
    else
        ifeq ($(PROCESSOR_ARCHITECTURE),AMD64)
            CXXLAGS=/W2 /wd4996 -D_CRT_SECURE_NO_WARNINGS
        endif
        ifeq ($(PROCESSOR_ARCHITECTURE),x86)
            CXXLAGS=/W2 /wd4996 -D_CRT_SECURE_NO_WARNINGS
        endif
    endif
else
    UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
        CXXLAGS += -fopenmp
        LDLIBS += -lpthread -lrt -lgomp
    endif
    ifeq ($(UNAME_S),Darwin)
        CXXLAGS +=
    endif
endif


#can add many src files if it is necessary in future (by putting spaces)
SRC= ./src/HMMDecoder.cpp ./src/HMMTrainer.cpp ./src/HMMGraph.cpp ./src/Polisher.cpp ./src/main.cpp 
OBJS=$(subst .cpp,.o,$(SRC))

default: all

all: $(OBJS)
	mkdir -p ./bin/
	$(CXX) $(LDFLAGS) -o ./bin/apollo $(OBJS) $(LDLIBS)

./src/main.o: ./src/Polisher.o ./src/main.cpp ./src/CommandLineParser.h
	$(CXX) $(CXXLAGS) -c -o ./src/main.o ./src/main.cpp

./src/Polisher.o: ./src/HMMDecoder.o ./src/HMMTrainer.o ./src/HMMGraph.o ./src/Polisher.cpp ./src/Polisher.h
	$(CXX) $(CXXLAGS) -c -o ./src/Polisher.o ./src/Polisher.cpp

./src/HMMDecoder.o: ./src/HMMGraph.o ./src/HMMDecoder.cpp ./src/HMMDecoder.h
	$(CXX) $(CXXLAGS) -c -o ./src/HMMDecoder.o ./src/HMMDecoder.cpp

./src/HMMTrainer.o: ./src/HMMGraph.o ./src/HMMTrainer.cpp ./src/HMMTrainer.h
	$(CXX) $(CXXLAGS) -c -o ./src/HMMTrainer.o ./src/HMMTrainer.cpp

./src/HMMGraph.o: ./src/HMMGraph.cpp ./src/HMMGraph.h ./src/HMMCommons.h ./utils/lib/libz.a ./utils/lib/libbz2.a ./utils/include/seqan/basic.h
	$(CXX) $(CXXLAGS) -c -o ./src/HMMGraph.o ./src/HMMGraph.cpp
	
./utils/lib/libz.a: ./utils/zlib-1.2.11.tar.gz
	mkdir -p ./utils/include/
	mkdir -p ./utils/lib/
	tar -xf ./utils/zlib-1.2.11.tar.gz -C ./utils/
	cd ./utils/zlib-1.2.11/; prefix=./ ./configure; make install prefix=../
	rm -rf ./utils/zlib-1.2.11/

./utils/lib/libbz2.a: ./utils/bzip2-1.0.6.tar.gz
	mkdir -p ./utils/include/
	mkdir -p ./utils/lib/
	tar -xf ./utils/bzip2-1.0.6.tar.gz -C ./utils/
	cd ./utils/bzip2-1.0.6/; make install PREFIX=../
	rm -rf ./utils/bzip2-1.0.6/

./utils/include/seqan/basic.h: ./utils/seqan.tar.gz
	rm -rf ./utils/include/seqan/
	mkdir -p ./utils/include/
	mkdir -p ./utils/lib/
	tar -xzf ./utils/seqan.tar.gz -C ./utils/
	mv ./utils/seqan ./utils/include/
# 	rm -rf ./utils/seqan-library-2.4.0

clean:
	$(RM) $(OBJS) ./bin/ ./utils/include/ ./utils/lib/ ./utils/share/ ./utils/man/ ./utils/bin/
