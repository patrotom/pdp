CXX=g++
CXXFLAGS=-std=c++11 -Wall -pedantic
LIBS=
NAME1=task1
NAME2=task2

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

task2: compile2
	out/$(NAME2).out
compile2: src/$(NAME2).cpp
	$(CXX) $(CXXFLAGS) -Ofast -fopenmp src/$(NAME2).cpp -o out/$(NAME2).out
debug2: src/$(NAME2).cpp
	$(CXX) $(CXXFLAGS) -g -fopenmp src/$(NAME2).cpp -o out/$(NAME2).out
tests2: compile2
	for file in test/*.txt; do \
	    echo "$$file\n"; \
	    out/$(NAME2).out 15 < "$$file"; \
	    echo "--------------------"; \
	done

clean:
	rm -rf out/*
