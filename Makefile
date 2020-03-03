CXX=g++
CXXFLAGS=-std=c++11 -Wall -pedantic
LIBS=
NAME1=task1

task1: compile1
	out/$(NAME1).out
compile1: src/$(NAME1).cpp
	$(CXX) $(CXXFLAGS) -Ofast src/$(NAME1).cpp -o out/$(NAME1).out
debug1: src/$(NAME1).cpp
	$(CXX) $(CXXFLAGS) -g src/$(NAME1).cpp -o out/$(NAME1).out
tests1: compile1
	for file in test/*.txt; do \
	    echo "$$file\n"; \
	    out/$(NAME1).out < "$$file"; \
	    echo "--------------------"; \
	done
clean:
	rm -rf out/*
