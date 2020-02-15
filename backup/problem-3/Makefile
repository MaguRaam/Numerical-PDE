CXX = g++
CXXFLAGS = -std=c++17 -g
LIBRARY_PATH  = ../library/lib
INCLUDE_PATH  = ../library/include
SRC = $(PWD)/src
INCLUDE = $(PWD)/include




main: objdir obj/helmholtz.o $(SRC)/main.cc
	$(CXX) $(CXXFLAGS) -o main $(SRC)/main.cc obj/helmholtz.o -I$(INCLUDE) -I$(INCLUDE_PATH) -L$(LIBRARY_PATH) -llibrary 

run: main
	./main

objdir:
	mkdir -p obj

obj/helmholtz.o: $(SRC)/helmholtz.cc
	$(CXX) $(CXXFLAGS) -c $(SRC)/helmholtz.cc -o obj/helmholtz.o -I$(INCLUDE) -I$(INCLUDE_PATH) -L$(LIBRARY_PATH) -llibrary


clean:
	rm -r obj
	rm plot/*.dat

#export LD_LIBRARY_PATH=/../lib
