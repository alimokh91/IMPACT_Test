***
# IMPACT
## by Rolf Henniger, Institute of Fluid Dynamics, ETH Zurich (henniger@ifd.mavt.ethz.ch)
### Mai 2005 - Dec 2011
***



--- general remarks on the fluid solver "IMPACT" ---



Scope:
------
The code is designed for numerical simulations of incompressible flows (according to the Navier-Stokes
equations) along with concentration fields and discrete particles on massively parallel shared-memory
computers. By default the Boussinesq approximation is applied, but non-Boussinesq simulations are also
feasible.
Memory bandwidth and parallel scalability are today by far the most limiting factors for the speed
of PDE solvers. Hence, the highest priority for the code development was the optimization/minimization of
accesses to the main memories as well as the minimization of data transfers between network nodes. To this
end, the simulation code is limited to Cartesian grids and uses (compact) high-order finite differences in
space and a three-step Runge-Kutta time integration. Velocity, pressure and concentrations are found
iteratively in each sub-time step using Krylov subspace primary solvers and geometric multigrid
preconditioning. The design of the simulation code is best suited for so-called "torus" or "mesh" network
computers in which the dimensionality of the network grid is at least the same as the dimensionality of the
domain decomposition.



Prerequisites:
--------------
- MPI-2 for inter-processor and inter-node communication (must be compiled for FORTRAN. MPI-1 works only
  without parallel hdf5)
- HDF5 1.6 or HDF5 1.8 for I/O (must be compiled for FORTRAN and for parallel applications. If you are not
  sure about your hdf5 installation check "[path-to-HDF5]/lib/libhdf5.settings".)
- Lapack for the initialization



Usage of the source files:
--------------------------
The usage and circulation of the simulation code is regulated in file "Licence.pdf". If you are not authorized
to deal with the code and if you wish to use it, please contact me (henniger@ifd.mavt.ethz.ch).

Generally, the development of the source code is still in motion and thus not "finished" yet which makes
sharing the code difficult. The aim for the future code development is to find a more or less general
framework of modules, subroutines and functions, however, the development of code details will probably never
end. The reasons for this are mainly the changing computer architectures and user demands. Also the
improvement of the existing code in terms of speed, readability and comprehensibility never stops.

To this end, the idea was to split the source files in files which contain the kernel (called "mod-files",
i.e. all files named "mod_..." plus the files "impact.f90", "alloc.f90", "sub_other.f90" and "compile") and
others which contain only user-specific code or input (called "usr-files", i.e. all files named "usr_..." plus
the file "Makefile"). Therefore, you should modify only the usr-files because fundamental code changes affect
mostly the kernel whereas the framework of the usr-files is much more stable. This convention is quite
important since it allows much easier updates of the base code by just replacing the mod-files. Hence, if you
miss an important kernel feature you should contact me such that we can work out an appropriate modification
of the mod-files (instead of customizing the mod-files by yourself). I think that this is currently the best
way to go until the basic code structure has more settled.

The usr-files contain different subroutines and functions. Most of them are closely connected to the kernel
subroutines which is correspondingly annotated in brackets ("basic subroutine" or "basic function"). Other
subroutines or functions serve just as examples which is also marked accordingly. In either case, you need to
specify the contents of the basic routines according to the instructions and explanations at the beginning of
the usr-files and inside the routines. Generally, you can add further content to these routines and you can
also add further functions and subroutines to the usr-files. Again, it is important that you don't customize
the mod-files if you want to benefit from future code updates.

The usr-files contain a simple example of a flow simulation which can be used as a basis for your own
configurations. The intention was to make learning the current programming conventions and the code handling
easier. Only the file "usr_stats.f90" is relatively lengthy and comprehensive (actually it is close to my own
production code). It contains a lot of basic subroutine calls to mod-files as well as some basic MPI
programming such that it may be a good template for programming parallel jobs. Generally, you can call all
functions or subroutines of the mod-files. However, you should be aware that they may change in future
releases.

By default, the code is set to a 5th/6th-order explicit or 10th-order compact spatial discretization of the
inner field (the convergence order on and near the boundaries goes down to 4th order). It is possible to
change the discretization with relatively little effort, however, such modifications are currently not very
intuitive and the possibility to change the discretization has also currently no priority in the code
development.



Compilation:
------------
To compile the code you should use the file "Makefile". It contains some compilation targets and macros for
different computers which I use, however, you will need to modify the macros for your own purposes. The
Makefile includes the more general file "compile" which should not be edited. More information about compiling
the code is given in the Makefile.

Note that there are some initialization subroutines in the kernel which may cause underflow floating point
exceptions at run-time. If such numbers are treated as zero, everything is fine and you can switch off
corresponding tripping compiler flags. Generally, it is beneficial to perform array-bound checks from time to
time in order to find possible bugs (note that not every code line has been tested by now ...). The
corresponding compiler flag is usually "-C".

For fast executables you have to set appropriate optimization flags. Particularly, ensure that vectorization
is enabled which is to my experience the most important measure for performance enhancements. Unrolling and
inlining are also helpful but secondary compared to vectorization.

You can compile the code for static as well as for dynamic memory allocation. Static allocation typically
yields somewhat faster executables (improvements of up to 10%), however, this option requires the adjustment
of the mod-file "mod_dims.f90" (which is problematic as stated before). This compilation option is explained
in the Makefile more in detail.



Running the executable:
-----------------------
As mentioned before, the design of the simulation code is best suited for torus or mesh network computers.
Maximum performance is obtained by mapping the mesh of sub-domains/processor-blocks to the network mesh. The
performance can be further enhanced if the torus/mesh network is combined with a so-called (fat-)tree network
to yield almost perfect weak parallel scalability (i.e. the execution time remains constant when the number of
grid points and the number of processors are increased simultaneously).

For the performance on multi-core processors it is crucial that only adjacent processor-blocks are assigned to
each processor/node in order to minimize the data transfers across the network adapters.

When you run a simulation for the first time, please check the file "log_iterations.txt" regularly during time
integration (the output of this file must be enabled, cf. "usr_config.f90"). The simulation is set up
correctly only if each of the (cascaded) iterative solver(s) converges within a "small" number of iterations
(such as one to ten).

There exist currently two interfaces to communicate with the executable during time integration at run-time.
First, the executable writes a small text file called "send_signal.txt" at the beginning of the time
integration which is read again after each full time step. It currently allows an instantaneous termination of
the time integration or a reduction of the maximum number of time steps or termination time. Restart files are
written if specified correspondingly. Second, you can provide a file termed "queue.txt" which contains the
number of seconds (integer) after which the job is terminated (note that the compiler and the operating system
must support the FORTRAN function "alarm").



Concluding remarks:
-------------------
If you observe a bug in the code (especially in the mod-files) please contact me as soon as possible since
other users will run the code with the same bug. Generally, any input or comment on the entire package
(including code examples, Makefile, this readme file etc.) is highly appreciated. Even notifications about
typos can be helpful.
