CXX = g++
CXXFLAGS = -std=c++17 -g
LIBRARY_PATH  = ../library/lib
INCLUDE_PATH  = ../library/include

all: main

main: main.cc $(LIBRARY_PATH)/liblibrary.so
	$(CXX) $(CXXFLAGS) -o main main.cc -I$(INCLUDE_PATH) -L$(LIBRARY_PATH) -llibrary

run: main
	./main

clean:
	rm main
	rm plot/*.dat
