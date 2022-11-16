

====================================
Generalised Laplace analytic example
====================================

This example solves the generalised Laplace equation in two dimensions
in a rectangular domain and compares it to an analytic solution.

Building the example
====================

The fortran version of the example can be configured and built with CMake::

  git clone https://github.com/OpenCMISS-Examples/laplace_equation
  mkdir laplace_equation-build
  cd laplace_equation-build
  cmake -DOpenCMISSLibs_DIR=/path/to/opencmisslib/install ../laplace_equation
  make

This will create the example executable "laplace_equation" in ./src/fortran/ directory.

Running the example
===================

Fortran version::

  cd ./src/fortran/
  ./laplace_equation NUMBER_X NUMBER_Y INTERPOLATION_TYPE

where NUMBER_X is the number of elements in the X direction, NUMBER_Y
is the number of elements in the Y direction and INTERPOLATION_TYPE is
the interpolation type to use.

Verifying the example
=====================

Results can be visualised by running `visualise.cmgui <./src/fortran/visualise.cmgui>`_ with the `Cmgui visualiser <http://physiomeproject.org/software/opencmiss/cmgui/download>`_.

