CXX=g++
CXXFLAGS=-std=c++11 -Wall -pedantic -g -Ofast
LIBS=
NAME1=task1

task1: compile1
	out/$(NAME1).out
compile1: src/$(NAME1).cpp
	$(CXX) $(CXXFLAGS) src/$(NAME1).cpp -o out/$(NAME1).out
clean:
	rm -rf out/*
