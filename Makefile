
CC=gcc
CXX=g++
RM=rm -f
#CPPFLAGS=-O3 -W -Wall
CPPFLAGS=-O3 -W -Wall
LDFLAGS=-g
LDLIBS=-L.

SRCS=gibbs.cpp
OBJS = $(SRCS:.cpp=.o)

# Run-time variables
EM = 4 # EM iteration
GS = 50 # Gibbs iteration
BE = 4 # Beta
LD = 0.5 # Lambda
ST = 125 # Tile width
HW = 0 # Gibbs method to use
RU = 1 # Number of runs 
LO = 0 # Try for different Gibbs iterations or not.

# CPU:
# ./gibbs 4 50 4 0.5 125 0 1 0
# ./gibbs 4 50 4 0.5 125 0 5 1

# FPGA:
# ./gibbs 4 50 4 0.5 24 1 1 0 
# ./gibbs 4 50 4 0.5 24 1 5 1

all: gibbs

run: gibbs
	./gibbs $(EM) $(GS) $(BE) $(LD) $(ST) $(HW) $(RU) $(LO)

gibbs: gibbs.o
	$(CXX) -o $@ gibbs.o $(CPPFLAGS)

gibbs.o: gibbs.cpp
	$(CXX) -c $< -o $@ $(CPPFLAGS) -std=c++11

clean:
	$(RM) $(OBJS) gibbs err
