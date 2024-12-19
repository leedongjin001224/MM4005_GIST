BoundaryProblem.exec : BoundaryProblem.o
	mkdir results; mpicxx BoundaryProblem.o -o BoundaryProblem.exec -L /home/sprng5/lib/ -l sprng; mpiexec -n 10 ./BoundaryProblem.exec; julia plotting.jl
BoundaryProblem.o : BoundaryProblem.cpp
	mpicxx BoundaryProblem.cpp -c -L /home/sprng5/lib/ -l sprng
clean:
	rm BoundaryProblem.o BoundaryProblem.exec; rm -r results