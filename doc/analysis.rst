
.. _label-analysis:

Analysis of 1CPN Simulations
=============================

Sorry! This part of the documentation is still under construction. 

Take a look in the  the `src/analysis` directory to see some sample analysis scripts. If you know a bit of C++ you should be able to copy one of these `.cpp` files and modify it to perform your analysis of interest.

You'll also want to checkout `src/include/trajectory_iterator.h`. This is a library we've put together that makes parsing of LAMMPS files a bit easier. Almost all of the analysis files in `src/analysis/` include `trajectory_iterator.h` and show the basic commands on how to use it.

