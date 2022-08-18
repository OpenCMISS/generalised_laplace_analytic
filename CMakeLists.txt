# Specify the minimum version of CMake, OpenCMISS itself requires
# at least version 3.4 of CMake so that will also constrian the 
# minimum version we are able to set here.  We can't continue if this
# is not the case so we may as well stop right here.
cmake_minimum_required(VERSION 3.4 FATAL_ERROR)

# Declare the project name and version number and specify the languages
# used.  We must specify the *C* language irrespective of whether we use 
# it or not as it is required by MPI.
project(XXXXXXXX VERSION 1.0.0 LANGUAGES C Fortran)

# Get CMake to find the OpenCMISS libraries.  Because the OpenCMISS libraries
# are not usually available in the system directories we will have to 
# specify where OpenCMISS libraries can be found.  We can do this through the
# command line by setting the argument *OpenCMISSLibs_DIR* to the location
# of the OpenCMISS libraries install directory.
find_package(OpenCMISSLibs 1.3.0 REQUIRED CONFIG)

# Add the subdirectory for further CMakeLists.txt files that define any
# executables that are to be built.  Python examples do not have a configure
# and build phase so CMake has no work to do thus we will not see any mention
# of Python in the CMake files.
add_subdirectory(src/fortran)
