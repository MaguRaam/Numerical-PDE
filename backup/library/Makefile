CXX = g++
CXXFLAGS = -std=c++17 -g
LIBRARY	= $(PWD)
INCLUDE = $(LIBRARY)/include
SRC = $(LIBRARY)/src

all: liblibrary.so

liblibrary.so: libdir objdir  obj/poisson.o obj/field.o
	$(CXX) -shared $(CXXFLAGS) -o lib/liblibrary.so obj/poisson.o obj/field.o

libdir:
	mkdir -p lib

objdir:
	mkdir -p obj



obj/poisson.o: $(SRC)/poisson.cc
	$(CXX) -fPIC $(CXXFlags) -o obj/poisson.o -I$(INCLUDE) -c $(SRC)/poisson.cc

obj/field.o: $(SRC)/field.cc
	$(CXX) -fPIC $(CXXFlags) -o obj/field.o -I$(INCLUDE) -c $(SRC)/field.cc

clean:
	rm -r obj
	rm -r lib
