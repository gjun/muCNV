#
# Generic Makefile for muCNV
#
#   Author: Goo Jun (goo.jun@uth.tmc.edu)

TARGET=  muCNV
SRCS= Main.cpp Error.cpp File.cpp deletions.cpp  cluster.cpp common.cpp duplications.cpp
OBJS= Main.o Error.o File.o deletions.o cluster.o common.o duplications.o
CFLAGS= -g -Wall -O2 -fPIC
DFLAGS= -D_FILE_OFFSET_BITS=64
CC= gcc
CXX= g++ 
INCLUDES=-I./tabix -I./tclap-1.2.1/include
LIBRARY= tabix/libtabix.a

.SUFFIXES:.cpp .o

.cpp.o:
	$(CXX) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all: $(TARGET) 

$(TARGET) : $(LIBRARY) $(OBJS)
	$(CXX) -o $@ $(OBJS) -L./tabix -ltabix -lm -lz -lssl


$(LIBRARY) : 
	make -C tabix lib

clean :
	-rm -f *.o $(TARGET) *~

cleandepend:
	makedepend -- $(DFLAGS) --

depend:
	makedepend -- $(DFLAGS) -- $(SRC) >/dev/null 2>&1

# DO NOT DELETE THIS LINE -- make depend depends on it
