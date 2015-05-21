
DIRECTORIES
===========
  bin/         -- c++ binaries
  src/         -- c++ source code
  matlab/      -- matlab code
  cartfd/      -- cache of cartesian coords for galewsky plots
  data/        -- input files for C++ code
  resdata/     -- solutions from C++
  nodesets/    -- preprocessed node sets
  orgnodesets/ -- original unprocessed node sets
  tree/        -- cached nearest neighbor information for node sets


SETUP
=====

BASH COMMANDS
  # checkout rbf-sw, superglue, and mpi-superglue (read-only)
  git clone https://github.com/tillenius/rbf-sw.git
  git clone https://github.com/tillenius/superglue.git
  git clone https://github.com/tillenius/mpi-superglue.git

  # or (with write-permission)
  git clone git@github.com/tillenius/rbf-sw.git
  git clone git@github.com:tillenius/superglue.git
  git clone git@github.com:tillenius/mpi-superglue.git

  # create links to superglue and mpi-superglue
  ( cd rbf-sw ; ln -s ../superglue/include/ superglue )
  ( cd rbf-sw ; ln -s ../mpi-superglue/include/ mpi-superglue )

  # make some directories
  ( cd rbf-sw ; mkdir -p nodesets tree data resdata cartfd )


QUICK GETTING STARTED
=====================

There is a small test that can be run to verify that everything works,
to show an example of how to run a simulation from beginning to end,
and to make sure nothing the latest changes didn't break anything.

BASH COMMANDS
  ./test.sh

DETAILS

The test script also takes an argument, a string of numbers that
represent the different steps below.

Example: ./test.sh 0    -- only compiles
Example: ./test.sh 1    -- only preprocess nodesets
Example: ./test.sh 0123 -- compile, preprocess, generate data,
                           run simulation, but don't plot results

There is also a ./testmpi.sh to test with MPI. It will start 4
MPI-processes on the same machine, each using all cores, and
thus runs slow, since the cores are oversubscribed.


WORKFLOW
========

1) Generate data in MATLAB
  1.1) mt_preprocess(filename, fd)
  1.2) mt_save(par)

2) Run simulation
  ./run <INPUT_FILE_DIR> <CHUNK_SIZE> <TIME_STEP> <END_TIME>

The parameters for the MATLAB commands are as follows:

  fd = stencil size
  par = parameters:
    par.test = 'galew' or 'tc5'
    par.n = number of node points
    par.fd = stencil size
    par.ep = shape parameter (epsilon)
    par.order = order of hyperviscosity
    par.gamma_c = coefficient for hyperviscosity

Example (generate data in MATLAB):

  mt_preprocess('..\orgnodesets\icos655362.mat', 31)
  par = struct('test', 'galew', 'n', 655362, 'fd', 31, 'ep', 40, 'order', 4, 'gamma_c', -0.1);
  mt_save(par)

Example (run):
  ./run galew-655362-31-ep40-o4-gc-0.1 5120 5 500


DETAILED WORKFLOW
=================

0) Compile
----------

BASH COMMANDS
  make



1) Preprocess nodeset
---------------------

MATLAB COMMANDS
  mt_preprocess('orgnodesets/x764128.mat', 31);

INPUT
  Original nodeset, stencil size

OUTPUT
  Creates 'nodesets/764128_31.mat' -- permuted nodeset
  Creates 'tree/764128_31.mat' -- nearest neighbor information

NOTES
 - Will convert node points to spherical coordinates and back with radius = 1.
   (This was done in the original code, so I guess it is desired.)
 - The ordering of the nodes depend on the stencil size.
 - The nearest neighbor information is calculated once and cached for performance reasons.


2) Generate input files for C++-code
------------------------------------

MATLAB COMMANDS
  clear;
  par = struct('test', 'galew', 'n', 764128, 'fd', 31, 'ep', 40, 'order', 4, 'gamma_c', -0.1);
  mt_save(par);

INPUT
  Parameters
  Reads 'nodesets/764128_31.mat'

OUTPUT
  Generates 'data/galew_764128-31-ep40-o4-gc-0.1/*'

NOTES
 - Calls galew_setup()
 - Calls galew_computeInitialCondition()
 -- galew_computeInitialCondition() needs int10001.mat
 - Calls mt_rbfmatrix_fd_hyper()
 - Compiles and uses mex_save to store matrices in C++-readable format.


3) Run C++ code
---------------

BASH COMMANDS
  ./run galew-764128-31-ep40-o4-gc-0.1 <CHUNK_SIZE> <TIME_STEP> <END_TIME>

  # CHUNK_SIZE = block size. For 8 nodes, 15 threads/node, use 764128/15/8 ~ 6368
  # TIME_STEP = in seconds
  # END_TIME = in seconds. For Galewsky, use 518400 for 6 days

INPUT
  Reads 'data/galew-764128-31-ep40-o4-gc-0.1/*'

OUTPUT
  Generates 'resdata/galew-764128-31-ep40-o4-gc-0.1/result.txt' -- final solution
  Generates 'resdata/galew-764128-31-ep40-o4-gc-0.1/results-t###.txt' -- solution at time ###


4) Plot results in MATLAB
-------------------------

MATLAB COMMANDS
  clear;
  par = struct('test', 'galew', 'n', 764128, 'fd', 31, 'ep', 40, 'order', 4, 'gamma_c', -0.1);
  dt = 5;
  saveplot(par, dt);

INPUT
  Parameters + timestep

OUTPUT
  Creates 'resdata/galew-764128-31-ep40-o4-gc-0.1-dt5.pdf'
  Creates 'resdata/galew-764128-31-ep40-o4-gc-0.1-dt5.png'

NOTES
  Uses "H = getsol(par, dt)" to read the solution into MATLAB
  Uses "mt_plot(par, H)" to plot the solution


UTILITIES
=========
  getnodes(par)      -- get processed nodeset for parameters 'par'
  gettree(par)       -- get neighbor information for parameters 'par'
  getsol(par, dt)    -- load solution into MATLAB
  getstart(par)      -- get initial solution
  loaddata(filename) -- load matrix from C++ into MATLAB

  mt_plot(par, H)    -- 2D plot of solution 
                        Caches cartesian coordinates for galewsky solution in 'cartfd/'
  mt_plot3(par, H)   -- 3D scatterplot solution 
                        Caches cartesian coordinates for galewsky solution in 'cartfd/'

  mt_plottest        -- plot the solutions from ./test.sh 
