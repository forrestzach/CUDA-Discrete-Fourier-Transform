CPSC 4780
Final Project
Forrest Zach (fzach)
Developed and tested on Palmetto cluster machines.


To compile this project simply utilize the following command:
nvcc Parallel_DFT.cu

To run this project use the following template:
./a.out <Problem Size> <OPTIONAL input flag>

-Recommneded problem size can be anything above 4000, expected execution time can 
    be found in the writeup
-Input flag can be one of two options
    "1" which produces random numbers for input into the algorithm, this is what
        the project uses by default if no input flag is used
    "2" which utilizes a constant real input of "1.0" for the real numbers and "0.0"
        for the complex numbers.

(Note, running with >10000 problem size becomes very cumbersome for the CPU portion
    and could take a very long time)
Some example command of running the code include:

Run program with problem size N = 8000 and random input 
./a.out 8000 1

Run program with problem size N = 10000 and constant real input
./a.out 10000 2


Expected Output:
There should be two banks of data which appear, first entry and last 20 entries of
    each output data.
When using input flag "1" (random data) there will be all random numbers throughout 
    the two banks of data.
When using input flag "2" (constant real data) there will be all zeroes throughout 
    the two banks of with the exception of the first entry ([0]) which should be
    approximately equal to the problem size N.

The way to verify that the GPU algorithm is working is to compare the entries between
    the bank of CPU and GPU data. The two banks should have matching data to a high
    degree of accuracy (-1e5).
An easy way to compare this data is note the side and index which you want to
    compare, say Xre[9999], for a problem size N = 10,000 and see that it is the same
    in both the CPU and GPU output
For input flag "2" the only one which will be the same is the first entry, as
    mentioned previously. This is due to the Fourier Transform finding that there is 
    N appearances of the input value "1.0" (oversimplification).

I have chosen to present the data in this matter so that it is easy for one to see
    that the output is indeed the same rather than a function which simply spits out
    a yes/no to the user.

At the bottom of the output there is a readout of the processing time for both the 
    CPU and GPU along with showing the ratio of speedup between the two.


Example Output:

[fzach@node1271 finalProj]$ nvcc Parallel_DFT.cu 
[fzach@node1271 finalProj]$ ./a.out 10000 1
DFTW calculation with N = 10000 
grid 20 block 512

Begin CPU Data!
Xre[0] = -53.424818, Xim[0] = 85.833911 
Xre[9980] = 31.476049, Xim[9980] = -111.836622 
Xre[9981] = 6.935445, Xim[9981] = 1.573554 
Xre[9982] = -3.240630, Xim[9982] = -35.670220 
Xre[9983] = -44.496288, Xim[9983] = -43.555883 
Xre[9984] = 54.418956, Xim[9984] = -38.290842 
Xre[9985] = 31.736564, Xim[9985] = 29.535148 
Xre[9986] = -12.403891, Xim[9986] = 51.126122 
Xre[9987] = 48.997278, Xim[9987] = 1.508871 
Xre[9988] = -212.104894, Xim[9988] = -16.705939 
Xre[9989] = 104.856748, Xim[9989] = 36.759333 
Xre[9990] = -27.285625, Xim[9990] = 84.328730 
Xre[9991] = -52.730343, Xim[9991] = -12.172082 
Xre[9992] = -87.746472, Xim[9992] = 17.416480 
Xre[9993] = 9.592329, Xim[9993] = -58.375589 
Xre[9994] = 59.717731, Xim[9994] = -25.352853 
Xre[9995] = -85.793740, Xim[9995] = 71.113798 
Xre[9996] = 57.221915, Xim[9996] = -3.842724 
Xre[9997] = -38.378324, Xim[9997] = -15.558056 
Xre[9998] = -176.220503, Xim[9998] = -22.933954 
Xre[9999] = -10.719111, Xim[9999] = 108.829138 

Begin GPU Data!
Xre[0] = -53.424818, Xim[0] = 85.833911 
Xre[9980] = 31.476049, Xim[9980] = -111.836622 
Xre[9981] = 6.935445, Xim[9981] = 1.573554 
Xre[9982] = -3.240630, Xim[9982] = -35.670220 
Xre[9983] = -44.496288, Xim[9983] = -43.555883 
Xre[9984] = 54.418956, Xim[9984] = -38.290842 
Xre[9985] = 31.736564, Xim[9985] = 29.535148 
Xre[9986] = -12.403891, Xim[9986] = 51.126122 
Xre[9987] = 48.997278, Xim[9987] = 1.508871 
Xre[9988] = -212.104894, Xim[9988] = -16.705939 
Xre[9989] = 104.856748, Xim[9989] = 36.759333 
Xre[9990] = -27.285625, Xim[9990] = 84.328730 
Xre[9991] = -52.730343, Xim[9991] = -12.172082 
Xre[9992] = -87.746472, Xim[9992] = 17.416480 
Xre[9993] = 9.592329, Xim[9993] = -58.375589 
Xre[9994] = 59.717731, Xim[9994] = -25.352853 
Xre[9995] = -85.793740, Xim[9995] = 71.113798 
Xre[9996] = 57.221915, Xim[9996] = -3.842724 
Xre[9997] = -38.378324, Xim[9997] = -15.558056 
Xre[9998] = -176.220503, Xim[9998] = -22.933954 
Xre[9999] = -10.719111, Xim[9999] = 108.829138 

Total CPU processing time(s): 20.061130
Total GPU processing time(s): 0.029815
GPU Speedup over CPU: 672.843310