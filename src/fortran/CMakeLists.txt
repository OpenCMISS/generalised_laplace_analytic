# Define the executable and adding the source files required to build the application.
# When we are defining fortran source files if we name the source files with the extension
# .F90 (using a capital *F* and not a lowercase one) CMake will automatically run the *C*
# preprocessor on the file for us, saving us from specifying that through other means.
# This is only necessary if the file requires the the *C* preprocessor.
add_executable(XXXXXXXX ZZZZZZZZ.F90)

# Set the required libraries for this executable.  Here we let CMake know which libraries
# are required to link the application.  Now because the *opencmisslibs* is known to CMake 
# as a target, there is nothing else we need to do.  The *opencmisslibs* target contains
# enough information for CMake to properly build and link the application with the OpenCMISS
# libraries.
#
# While https://github.com/OpenCMISS/iron/issues/88 is not fixed we also need to specify that
# MPI should also be linked into the application.  Once this issue is resolved we will no longer
# be required to add this as a link library.  This only applies to examples making use of Iron, 
# if the example is only using the Zinc library then MPI is not required at all.
target_link_libraries(XXXXXXXX PUBLIC opencmisslibs)
