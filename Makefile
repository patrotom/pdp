CXX=g++
CXXFLAGS=-std=c++11 -Wall -pedantic
CXXSET=$(CXX) $(CXXFLAGS)
COMPILE_SEQ=$(CXXSET) -Ofast
COMPILE_PAR=$(CXXSET) -Ofast -fopenmp
COMPILE_MPI=OMPI_CXX=g++ mpic++ -fopenmp -Ofast $(CXXFLAGS)
DEBUG_SEQ=$(CXXSET) -g
DEBUG_PAR=$(CXXSET) -g -fopenmp
NAME1=task1
NAME2=task2
NAME3=task3
NAME4=task4
NP=4

define tests
	for file in test/*.txt; do \
	    echo "$$file\n"; \
		if [ $(1) = "task4" ]; then \
			mpirun -np $(2) out/$(1).out $(3) $(4) "$$file"; \
		else \
	    	out/$(1).out $(2) < "$$file"; \
		fi; \
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
	out/$(NAME2).out 15
compile2:
	$(COMPILE_PAR) src/$(NAME2).cpp -o out/$(NAME2).out
debug2:
	$(DEBUG_PAR) src/$(NAME2).cpp -o out/$(NAME2).out
tests2: compile2
	$(call tests,$(NAME2),15)

task3: compile3
	out/$(NAME3).out 1000
compile3:
	$(COMPILE_PAR) src/$(NAME3).cpp -o out/$(NAME3).out
debug3:
	$(DEBUG_PAR) src/$(NAME3).cpp -o out/$(NAME3).out
tests3: compile3
	$(call tests,$(NAME3),1000)

compile4:
	$(COMPILE_MPI) src/$(NAME4).cpp -o out/$(NAME4).out
tests4: compile4
	$(call tests,$(NAME4),$(NP),15,1000)

clean:
	rm -rf out/*
