#
# Generic Makefile for muCNV
#
#   Author: Goo Jun (goo.jun@uth.tmc.edu)

TARGET=  bin/muCNV
SRCS := $(wildcard muCNV/*.cpp)
OBJS := $(addprefix obj/,$(notdir $(SRCS:.cpp=.o)))
# OBJS= Main.o Error.o File.o deletions.o cluster.o common.o duplications.o
CFLAGS= -g -Wall -O2 -fPIC
DFLAGS= -D_FILE_OFFSET_BITS=64
CC= gcc
CXX= g++ 
INCLUDES=-I./tabix -I./tclap-1.2.1/include

# .SUFFIXES:.cpp .o

#.cpp.o:
#	$(CXX) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

obj/%.o: muCNV/%.cpp
	$(CXX) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all: $(TARGET) 

$(TARGET) : $(OBJS)
	$(CXX) -o $@ $(OBJS) -L./lib -ltabix -lm -lz

clean :
	-rm -f $(OBJS) $(TARGET) *~

cleandepend:
	makedepend -- $(DFLAGS) --

depend:
	makedepend -- $(DFLAGS) -- $(SRC) >/dev/null 2>&1

# DO NOT DELETE THIS LINE -- make depend depends on it
