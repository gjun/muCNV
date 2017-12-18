#
# Generic Makefile for muCNV
#
#   Author: Goo Jun (goo.jun@uth.tmc.edu)

TARGET=  bin/muCNV
SRCS := $(wildcard muCNV/*.cpp)
OBJS := $(addprefix obj/,$(notdir $(SRCS:.cpp=.o)))
CFLAGS= -g -Wall -O2 -fPIC -std=c++0x
DFLAGS= -D_FILE_OFFSET_BITS=64
CC= gcc
CXX= g++ 
INCLUDES= -I./tclap-1.2.1/include -I./htslib

# .SUFFIXES:.cpp .o

#.cpp.o:
#	$(CXX) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

obj/%.o: muCNV/%.cpp
	$(CXX) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all: $(TARGET) 

$(TARGET) : $(OBJS)
	$(CXX) -o $@ $(OBJS) -lhts -lm -lz

clean :
	-rm -f $(OBJS) $(TARGET) *~

cleandepend:
	makedepend -- $(DFLAGS) --

depend:
	makedepend -- $(DFLAGS) -- $(SRC) >/dev/null 2>&1

# DO NOT DELETE THIS LINE -- make depend depends on it
