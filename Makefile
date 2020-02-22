CXX=g++
CXXFLAGS=-std=c++11 -Wall -pedantic
LIBS=

task1: compile1
	out/task1.out
task1_name = src/task1.cpp
compile1: $(task1_name)
	$(CXX) $(CXXFLAGS) $(task1_name) -o out/task1.out
clean:
	rm -rf out/*
