CXX = g++
CXXFLAGS = -std=c++17 -g
LIBRARY	= $(PWD)
INCLUDE = $(LIBRARY)/include
SRC = $(LIBRARY)/src

all: liblibrary.so

liblibrary.so: libdir objdir obj/field.o obj/fourier.o
	$(CXX) -shared $(CXXFLAGS) -o lib/liblibrary.so obj/field.o obj/fourier.o

libdir:
	mkdir -p lib

objdir:
	mkdir -p obj

obj/field.o: $(SRC)/field.cc
	$(CXX) -fPIC $(CXXFlags) -o obj/field.o -I$(INCLUDE) -c $(SRC)/field.cc

obj/fourier.o: $(SRC)/fourier.cc
	$(CXX) -fPIC $(CXXFlags) -o obj/fourier.o -I$(INCLUDE) -c $(SRC)/fourier.cc

clean:
	rm -r obj
	rm -r lib 
	
