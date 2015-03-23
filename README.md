#Reduced-Rank Randomized SVD package (RRSVD-Package) 
This repository contains the implementation of the method described in the paper 'Fast Time-Evolving Block-Decimation algorithm through Reduced-Rank Radomized Singular Value Decomposition'

[![Linux Build Status](https://travis-ci.org/qubit-ulm/ebs.svg?branch=master)](https://travis-ci.org/qubit-ulm/ebs)
[![Windows Build Status](https://ci.appveyor.com/api/projects/status/of2bssgx58e8qyi3?svg=true)](https://ci.appveyor.com/project/jrosskopf/ebs)

# Introduction 
From the paper
# Usage
The RRSVD routine comes in three versions: 

*CBLAS-LAPACKE-based (BL-RRSVD)
*Intel MKL-based (MKL-RRSVD)
*Nvidia CUDA-based (CUDA-RRSVD). 

The interface to the library is uniform for all the different implementations. For every implementation, we provide single/double and real/complex RRSVD. 

The RRSVD-Package contains the three implementations of RRSVD. For the sake of modularity we provide them in three different folders. In each folder we put the library source files (RRSVD folder), some examplifying C++ code showing how to use the library (Examples), a Fortran wrapper that allows to call the RRSVD routine directly from Fortran code (FortranWrapper) and some sample matrices to use for testing (Data). All the folders and subfolders have README files illustrating their content and use. Where necessary, a Makefile is provided as well.

## General remarks
The RRSVD routine, for the time being, works only on square (n x n) and tall (m x n, m>n) matrices. In order to handle wide matrices (m x n, n>m) matrices, you must transpose (or adjoin) them in your own code and then manipulate the outpot provided by RRSVD as needed.
The RRSVD code, and most of the example programs, are written in C++11 and use some of the new libraries that come with this version of the language. 
# Building

## Requirements
RRSVD-Package has the following requirements on Unix systems

+BL-RRSVD
gcc >= 4.7
CBLAS and LAPAcke libraries (link)

+MKL-RRSVD
icc >= 14.0
Intel MKL library (Version)

+CUDA-RRSVD
CUDA >= 5.0
CUDA capable device with compute capability >= 2.0
CULA (Version and link)
gcc >=4.7

#Warning
While the building of the MKL- and CUDA-RRSVD is rather standard, the creation of the BL-RRSVD library requires a little bit of attention. In the RRSVD_BL_Release folder we include a precompiled version of LAPacke and CBLAS. The libraries have been compiled with the GNU gcc 4.7 compiler and should work on most Linux distributions with gcc 4.7 installed. Before compiling this libraries we have set complex types to be represented as a struct, with the real and imaginary part accessible via the field names .re and .im. This is one of the options provided by the LAPacke installation. If you recompile the LAPacke and CBLAS libraries pay attention to set the right complex numbers representation in the appropriate configuration files.



## Build Process
RRSVD-Package uses GNU-Make as a build system and allows several flexible build configuration options.
Refer to the README files in the folder and subfolders of the version you need to use to have further information about how to compile and use the code.

First clone the git repository using the command line git client or your favorite gui interface. By default this will create a directory `ebs` in the current working directory. You should change to this repository root, after the clone operation has finished.

$ git clone --recursive https://github.com/qubit-ulm/ebs.git 
$ cd ebs 

Then, make a build directory.  The directory can have any name, not just 'build', but 'build' is sufficient.

$ mkdir build                                     
$ cd build

The next step is to run CMake to configure the project.  Running CMake is the equivalent to running `./configure` with autotools. In this step CMake will also download the boost library, extract and bootstrap it. Boost is a dependency of EBS and is used for different purposes. All necessary parts of boost are statically linked, so a boost installation is not necessary on the computer.

$ cmake ..

Once CMake has finished, the process building the executable depends on the opteration system you are on. On Unix systems the `make` program will start the build.

$ make

On Windows systems you should run from the Visual Studio Command Prompt `msbuild`

$ msbuild ebs.sln

This will build the denoising as well as the graph processing part of EBS.

If the build fails and you cannot figure out why, please file a ticket on the github page the EBS developers will quickly help you figure it out.


# Attribution

Please cite `, (2015)
<http://arxiv.org/abs/1202.3665>`_ if you find this code useful in your
research and add your project or publication to `the users list
<https://github.com/qubit-ulm/ebs/blob/master/docs/users.md>`_.

The BibTeX entry for our paper is::

@article{ebs,
author = {},
title = {An Energy Based Scheme for Reconstructing Piecewise Constant Signals underlying Molecular Machines},
journal = {Biophysical Journal},
year = 2015,
volume = ,
pages = {},
eprint = {},
doi = {}
}

# License
The EBS project is licensed to you under the Apache License, Version 2.0 (the "License"); you may not use the code except in compliance with the License. You may obtain a copy of the License at

[http://www.apache.org/licenses/LICENSE-2.0](http://www.apache.org/licenses/LICENSE-2.0)

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language governing permissions and limitations under the License.

## Required libraries

While CBLAS and LAPacke are included in the distribution folder (and usable under the conditions LICENCE), the MKL and CUDA-enabled libraries require the download of the libraries from the

Each program has documentation which details how to use it, what each of the parameters are, and how to use them:

    $ ./bin/lambdaopt --help

    If one wants to see more details of the inner workings of a program, adding the debug flag might increase the verbosity of the output:

        $ ./bin/lambdaopt --debug

        The [Matrix Market File format](http://math.nist.gov/MatrixMarket/formats.html) as the default input output format. The format ASCII based, allows comment lines, which begin with a percent sign. We use the "array" format for general dense vectors. Details how to handle this format in python or matlab can be found in the particular demos.

        If you find a bug in EBS, have a problem using it or have a question about the method in general, feel free to open a github issue. The development team will try to answer the problem or fix the issue in a timely fashion. 

## Determining the regularization parameter lambda
The first step in a typical processing chain of a piecewise constant singnal is to determine a reasonable choice for the regularization pparameter lambda. Typically this requires a lot of twiddling. The program `lambdaopt` implements the heuristic we proposed in the paper to chose this parameter automatically.

The usage is as simple as:

    $ ./bin/lambdaopt $NOISY_DATA

    `$NOISY_DATA` is the path to the matrix market formated file containing the noisy input vector. The sole output of this command is a floating point value of the optimal lambda.

## Denoising the dataset
Having the optimal lambda at hand, the next step in the processing chain is solving the TVDN on the noisy input file:

    $ ./bin/denoising --lambda $LAMBDA_OPT $NOISY_DATA > $DENOISED_DATA

    `$NOISY_DATA` is again the path to the noisy input vector. `$LAMBDA_OPT` contains the value for the regularization parameter. Typically this is the output of the `lambdaopt` program. The solution of the TVDN optimization problem is written to stdout by default, but a file destination can be chosen by supplying the `--output` parameter. In the example above the output to stdout is redirected to a file `$DENOISED_DATA`.

## Clustering to a set of predefined levels 
The prerequisit for clustering is having a noise free signal, which we assume in a matrix market file `$DENOISED_DATA`. The task for this step as described in the paper is to cluster or assign a level from a predefined set to each sample in the noise free vector. So as a second prerequisite we require a vector containing the level set. 
One option would be to generate a problem specific set of levels which incorporates prior knowledge about the steps to expect in the signal. Another option is to simply build a equidistant grid between the minimal and maximal value in `$DENOISED_DATA`. This is exactly what the `level_generator` program is good for:

    $ ./bin/level_generator --level-distance $DISTANCE $DENOISED_DATA > $LEVEL_DATA

    The output, which is written to stdout by default, contains a vector with a grid with spacing of $DISTANCE. In the example the output is redirected to a file $LEVEL_DATA.

    Now everything is at hand to start the clustering process:

        $ ./bin/graph_processing --input $DENOISED_DATA --levels $LEVEL_DATA > $CLUSTERED_DATA  

        Where in the example the `$DENOISED_DATA` is the output of `denoising` and the `$LEVEL_DATA` is a vector containing the level set in matrix market format. The output is written to stdout by default and contains a vector of the same length as `$DENOISED_DATA`. 

        The energy, which is minimized via graph cuts by the `graph_processing` program, consists of three components: A data term, a smoothing term and a step height prior term. The relative weight, between this terms, can be adjusted by setting the parameters `--rho-d` (default value 100), `--rho-s` (default value 10) and `--rho-p` (which is disabled by default). We found the default values to work well on our test datasets. But the values may need to be tweaked due to the
        specific problem.

        To turn on the step height prior, the parameter `--rho-p` has to be chosen > 0. Futher the parameter `--prior-distance` has to be set to the distance of two adjacent steps, which should NOT be penalized.


# Demos
The source package contains demo code which demonstrates the usage of the above described programs form either [Matlab](http://www.mathworks.com/products/matlab/) or [Python](http://www.python.org). The Matlab demo is in the `matlab/demo.m` file. In this file the simulation code is used to create new test data for each run. The Python one in `python/demo.py`. Here we use pre-generated test data from the `noisy_data.mm` file in the same directory. The result of each demo run should be a plot
which should look like the picture below:

![Output plot of the demos](https://github.com/qubit-ulm/ebs/blob/master/demo_plot_full.png)
![Takeout of the demo plot](https://github.com/qubit-ulm/ebs/blob/master/demo_plot_outtake.png)


