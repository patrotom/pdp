CXX=g++
CXXFLAGS=-std=c++11 -Wall -pedantic
CXXSET=$(CXX) $(CXXFLAGS)
COMPILE_SEQ=$(CXXSET) -Ofast
COMPILE_PAR=$(CXXSET) -Ofast -fopenmp
DEBUG_SEQ=$(CXXSET) -g
DEBUG_PAR=$(CXXSET) -g -fopenmp
NAME1=task1
NAME2=task2
NAME3=task3

define tests
	for file in test/*.txt; do \
	    echo "$$file\n"; \
	    out/$(1).out $(2) < "$$file"; \
	    echo "--------------------"; \
	done
endef

task1: compile1
	out/$(NAME1).out
compile1:
	$(COMPILE_SEQ) src/$(NAME1).cpp -o out/$(NAME1).out
debug1:
	$(DEBUG_SEQ) src/$(NAME1).cpp -o out/$(NAME1).out
tests1: compile1
	$(call tests,$(NAME1))

task2: compile2
	out/$(NAME2).out
compile2:
	$(COMPILE_PAR) src/$(NAME2).cpp -o out/$(NAME2).out
debug2:
	$(DEBUG_PAR) src/$(NAME2).cpp -o out/$(NAME2).out
tests2: compile2
	$(call tests,$(NAME2),15)

task3: compile3
	out/$(NAME3).out
compile3:
	$(COMPILE_PAR) src/$(NAME3).cpp -o out/$(NAME3).out
debug3:
	$(DEBUG_PAR) src/$(NAME3).cpp -o out/$(NAME3).out
tests3: compile3
	$(call tests,$(NAME3),1000)

clean:
	rm -rf out/*
