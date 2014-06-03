FFVC-MINI
=========

* version: 1.0.0 (based on FFVC 0.9.2)
* date: 2014/06/03
* contact: miniapp@riken.jp


About FFVC-MINI and FFVC
-------------------------

The mini application FFVC-MINI is based on an FFVC simulation program
which solves the 3D unsteady thermal flow of the incompressible fluid,
being developed by Institute of Industrial Science, the University of Tokyo.
FFVC discretizes and solves the 3D incompressible Navier-Stokes equation
in the even-spaced grid point orthogonal space using finite volume method.
Refer to the user guide "ffvc_ug.pdf" for details.

Contact point of the FFVC full version is: Prof. Kenji Ono <keno@riken.jp>


Compilation
-----------

To build and run this program, C++ and Fortran90 compilers that support OpenMP,
MPI library and GNU Make are required.

 1. Obtain the package and extract its contents.

 2. Go to the src directory and edit the file "make_setting" according to
    your compilers and libraries.
    There are several example "make_setting.*" files:
    * make_setting.fx10 : For K and FX10-based systems
    * make_setting.gcc  : For the GCC compilers
    * make_setting.pgi  : For using PGI compilers

 3. Run make command in the src directory. 
    After successful compilation, an executable file named `ffvc_mini` is 
    generated in the bin directory.


Testing
-------

A test shell script is provided in the test directory.
To run the test interactively, simply run "./go.sh" in the directory.
Or run "make test" in the src directory.
This test script runs the program `ffvc_mini` with 8 MPI processes and 
2 OpenMP threads per process, and compares the computed results with reference data.


What FFVC-MINI does - brief explanation
---------------------------------------

FFVC-MINI is a subset of FFVC full application.
FFVC-MINI includes the minimum program modules that can carry out
the intrinsic performance monitor test (PMT) run.
The intrinsic PMT run computes the cavity flow in the 3D cuboid domain. 
The computation is done in the performance measurement mode, 
performing the fixed number of iterations irrespective of the convergence.
A detailed explanation of the intrinsic PMT
analysis is given in the FFVC user guide page 18.

The full version of FFVC has a rich set of user defined boundary
conditions, including velocity flux, heat flux, temperature, etc.
Outer, i.e. external, boundary conditions can be defined on the
bounding planes of the computed domain.
Local, i.e. internal, boundary conditions can be defined in the
other internal domain.

FFVC-MINI is specifically targeted to the intrinsic PMT analysis
for the cavity flow in the 3D cuboid domain. 
The essential program structure is maintained as is, but the
functions and routines that are not used for this analysis are all
cut off from the package.
The nested loop structure of the full version is as below:

    Time step loop  {

        - compute convection term, diffusion term
        - compute Poisson source term
  
        V-P iteration loop {
            Poisson iteration loop {
                - SOR computation
            }
  
            - Velocity update
            /*
            - Poisson source term update for pressure loss
              (removed in the FFVC-MINI)
            */
        }

    }

This structure is maintained in FFVC-MINI as well.
**Poisson iteration loop** solves pressure Poisson equation.
Outer **V-P iteration loop** handles the pressure loss boundary conditions
expressed as the function of the velocity. 
Since FFVC-MINI does not cover the pressure loss model,
it is enough to conduct only one cycle of V-P iteration.
Users, however, can explicitly specify the number of iterations
for each loops.
Typical numbers of iterations for the FFVC full version range
from 5 to 100 for V-P loop, and from 20 to 1000 for Poisson loop.

The following schemes are used in FFVC-MINI:

 - 1st order explicit Euler scheme for time integration.
 - 3rd order MUSCL scheme for convection term.
 - Strided memory access dual-colored SOR scheme for iterative solver.


How to run the program
----------------------

### Command line arguments ####

There is no input file for FFVC-MINI.
All the necessary parameters should be given as the command line options.
The list of command line options are as below:

    --scale=str     choose "strong" or "weak" (mandatory)
    --size=int      number of cells per domain edge (mandatory)
    --division=str  domain partitioning in "LxMxN" (default 1x1x1)
    --step=int      computing time steps (default 20)
    --dt=float      time step value in CFL unit (default 0.2)
    --comm=str      "sync" or "async" communication mode (default async)
    --p_itr=int     max. number of Poisson iterations (default 30)
    --vp_itr=int    max. number of V-P iterations (default 20)
    --practical     flag to enable the practical computation mode (see remark below)
    --output_interval=int
                    interval of results output (default 0 == no output)
    --help          print help message and exit


#### Sizing options ####

The `--scale=strong` option corresponds to a strong scaling run.
The `--scale=weak` option corresponds to a weak scaling run.
Depending on this option, the meaning of the size option varies.

The command line below corresponds to a strong scaling run,
modeling the total space of 128x128x128 cells, dividing it into
2x2x2 domains (=8 domains).

    $ ffvc_mini --scale=strong --size=128 --division=2x2x2

The command line below corresponds to a weak scaling run,
modeling a domain of 128x128x128 cells, running 2x2x2 domains,
i.e. 256x256x256 cells in total space.

    $ ffvc_mini --scale=weak --size=128 --division=2x2x2

#### Practical computation flag ####

By default (without `-practical` option), the program runs in
the performance measurement mode.
In this mode, both of the Poisson loop
and V-P loop are iterated to their maximum numbers.
If `-practical` option is given, the program is executed under a condition that reflect more realistic fluid characteristics. Also, the iteration count is not pre-fixed, but is dynamically determined based on the convergence criteria.

### Output file ####

The file "history_base.txt" will contain the convergence status,
statistics of the physical variables, and elapsed time for each
time step.

If the integer value of `--output_interval` option is larger than 0,
the instantaneous pressure and velocity values will be output
to files. The files are created for each of the time steps and
for each of the processes with the following naming convention.

  * prs_.dfi   ... pressure index file (text)
  * vel_.dfi   ... velocity index file (text)
  * prs_TTTTTTTTTT_idNNNNNN.sph  ... pressure values (binary)
  * vel_TTTTTTTTTT_idNNNNNN.sph  ... velocity values (binary)

where TTTTTTTTTT is a 10 digit number indicating the time step,
and NNNNNN is a 6 digit number indicating the process number.


Target exa-scale problem setting
--------------------------------

The current model setting targeted at exa-scale simulations is as follows:

- 1mm cell size
- 10^11 cells
- 3-second real time with 10^6 simulation time steps

The target performance of this simulation is to compute the above model within 3 hours.


Modification from the original FFVC full application
-----------------------------------------------------

* Reduced the size of the code by limiting the functionality:

    - Physical model is limited to the intrinsic PMT analysis model
      (3D cavity flow performance measurement)

    - The CFD solver algorithm is limited to FS_EE_EE
      (Fractional Step, Euler Explicit)

    - Iterative solver library is limited to SOR2SMA 
      (Strided memory access dual-colored SOR scheme)

    - Monitoring and sampling feature is removed except for history_base.txt

* Changed the parameter input method from configuration file to command line options.

* Other code clean-up.


Further detail information
------------------------

### Hot spot ###

The "hot spots" are in the SOR computation in the Poisson iteration loop.

  * psor2sma_core_     SOR core part (N_vp*N_p*2 times)
  * poi_residual_      residual computation (N_vp*N_p times)
  * sma_comm_, sma_comm_wait_    ghost cell communication (N_vp*N_p*2 times)

The values in the parentheses show the execution counts, 
where N_vp and N_p are the numbers of V-P and Poisson iterations, respectively. 

The realistic fluid simulation may finish in much fewer iterations.
In such a case, following parts can also be hot spots.

  * update_vec_    velocity update (N_vp) outside Poisson loop
  * pvec_muscl_    convection term (1) outside the V-P loop


### Bit-array based compressed coefficients ###

In FFVC the information of the cell, 
including medium property, boundary condition and neighboring cell information,
is encoded into 32-bit integer arrays as bit flags.
The following arrays are used in the program.

 * d_bcd  medium information
 * d_bcp  pressure related information
 * d_bcv  velocity related information

The coefficient matrix of the linear equations is decoded from these
arrays when needed at run time.


### Single precision and double precision  ###

In the default setting, FFVC-MINI uses single-precision floating point calculation,
which is often sufficiently precise for incompressible fluid simulations.
If the double-precision computation is necessary, the program can be built
as double-precision version as below:

 1. Edit file "src/make_setting" 
    * enable `CXXFLAGS = -D_REAL_IS_DOUBLE_`
    * set the appropriate compiler option which makes default REAL 8 bytes long in `F90FLAGS`
     (`-r8` for Intel compiler, `-CcdRR8` for FX10/K, etc )

 2.  Do "make clean; make" in the src directory

