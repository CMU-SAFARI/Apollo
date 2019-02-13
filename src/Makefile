CXX=g++
RM=rm -rf

#include directories seperated by a space
INC=../utils/include/
INC_PARAMS=$(foreach d, $(INC), -I$d)

#lib directories seperated by a space
LIB=../utils/lib/
LIB_PARAMS=$(foreach d, $(LIB), -L$d)

#paramaters as suggested by seqan
CXXLAGS=$(INC_PARAMS) -std=c++14 -O3 -DNDEBUG -DSEQAN_HAS_ZLIB=1 -DSEQAN_HAS_BZIP2=1 -W -Wall -pedantic
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
SRC=main.cpp
OBJS=$(subst .cpp,.o,$(SRC))

default: all

all: $(OBJS)
	mkdir -p ../bin/
	$(CXX) $(LDFLAGS) -o ../bin/apollo $(OBJS) $(LDLIBS)

main.o: main.cpp CommandLineOptions.h ../utils/lib/libz.a ../utils/lib/libbz2.a ../utils/include/seqan/basic.h
	$(CXX) $(CXXLAGS) -c -o main.o main.cpp $(LDLIBS)

../utils/lib/libz.a: ../utils/zlib-1.2.11.tar.gz
	mkdir -p ../utils/include/
	mkdir -p ../utils/lib/
	tar -xf ../utils/zlib-1.2.11.tar.gz -C ../utils/
	cd ../utils/zlib-1.2.11/; prefix=../ ./configure; make install prefix=../
	rm -rf ../utils/zlib-1.2.11/

../utils/lib/libbz2.a: ../utils/bzip2-1.0.6.tar.gz
	mkdir -p ../utils/include/
	mkdir -p ../utils/lib/
	tar -xf ../utils/bzip2-1.0.6.tar.gz -C ../utils/
	cd ../utils/bzip2-1.0.6/; make install PREFIX=../
	rm -rf ../utils/bzip2-1.0.6/

../utils/include/seqan/basic.h: ../utils/seqan-library-2.4.0.tar.xz
	rm -rf utils/include/seqan/
	mkdir -p ../utils/include/
	mkdir -p ../utils/lib/
	tar -xf ../utils/seqan-library-2.4.0.tar.xz -C ../utils/
	mv ../utils/seqan-library-2.4.0/include/* ../utils/include/
	rm -rf ../utils/seqan-library-2.4.0

clean:
	$(RM) $(OBJS) ../bin/* ../utils/include/* ../utils/lib/* ../utils/share/ ../utils/man/ ../utils/bin/
