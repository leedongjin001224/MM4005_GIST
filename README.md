When you run codes here, you can ignore below messages.
"sh: latex: command not found; latex: failed to create a dvi file"
I tried to use latx in plots, but those messages would be appear if it is not installed.

Here are codes on two problems.

First problem is about the Buffon's needle problem.
Enter "julia Buffon.jl" and then enter two parameters to run it; it takes around 25 minutes.
Parmeters are about the distance between two neighbor lines and the length of the needle.
This code uses two distinct methods to approximate pi.
Answers on Buffon's needle problem would be printed on command line, with each approximated pi.
Two plots "Wallis_product.pdf" and "Hit-and-Miss_approximation.pdf" would be saved in current directory.
Plots show error per iteration in log-scale with 95% confidence range of the slope.
Note that Wallis product is deterministic method, which computes product series.

Second problem is about the probabilistic potential theory.
Enter "make" (enter "make clean; make" if exec file already exists) to run.
You should enter random number generator type and number of bins (try 0, 1000 each for test).
This executes two codes and might takes less than 10 minutes.
First c++ code for simulation would save results in directory ./results.
Then plotting.jl would be run if you entered "make" before.
This plotting code would save 6 plots, each names with "{bin*num}*{some plot name}.pdf" format where bin_num is the number of bins you entered at first.

Checking points

1. MPI code in "BoundaryProblem.cpp", "Makefile" in current directory
2. SPRNG in "BoundaryProblem.cpp", quasi-random number in "Buffon.jl"
3. Visualization plots would be saved in current directory if codes were done.
4. Uploaded on github.
