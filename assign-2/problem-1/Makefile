CXX = g++
CXXFLAGS = -std=c++17 -g
LIBRARY_PATH  = ../library/lib
INCLUDE_PATH  = ../library/include

all: main

main: main.cc
	$(CXX) $(CXXFLAGS) -o main main.cc

run: main
	./main

clean:
	rm main
	rm plot/*.dat
