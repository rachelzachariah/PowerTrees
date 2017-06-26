-----------
BASIC USAGE
-----------

The program random_script.sage will generate random power trees and solve them using our algorithm. Its core usage is

sage random_script.sage -n [Number of vertices] -md [Maximum degree of the tree] -reps [Number of iterations]

The program will then generate random trees of the given number of vertices and will continue to generate them until it finds one with degree at most the bound given. It then generates random susceptances and solves using our algorithm. This entire process will be repeated for the given number of times. It will then print how long (on average) the solving process takes.

Example usage:

sage random_script.sage -n 6 -md 5 -reps 10
Average time for 6 nodes : 0.0281903902054

The second line here is the output of the program.

IMPORTANT NOTE: Python's timeit function actually runs multiple trials in order to get the best estimate of how long the function call takes. This may result in the program taking longer to run than its output suggests.


--------------
ADVANCED USAGE
--------------

We may want to compare our algorithm to others. In order to do so, I've included functionality that will write the equations that were generated to text files.

sage random_script -n N -md M -reps R -w filename

will perform the above process. The only difference is that for 0 <= i <= R-1, it will write the equations of the power system to the file filename_i_eqs.txt". For example:

sage random_script.sage -n 6 -md 5 -reps 2 -w example
Average time for 6 nodes : 0.0281903902054

Now if we look at our files, we will see the files
example_0_eqs.txt
example_1_eqs.txt

We can now apply other solvers. The equations are written in a format that phc can accept. If you call

phc -b example_0_eqs.txt example_0_roots.txt

this will generate a file example_0_roots.txt that contains all the roots of the system. At the bottom the file will contain information about how long it took to run.