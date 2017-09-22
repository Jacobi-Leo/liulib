cmake_minimum_required(VERSION 3.5)

project(matrix)
enable_language(Fortran)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -llapack -lblas")
set(CMAKE_Fortran_COMPILER_ID "Intel")
set(CMAKE_BUILD_TYPE "Debug")

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_COMPILER "gfortran")
  set(CMAKE_C_COMPILER "gcc")
  set(CMAKE_CXX_COMPILER "g++")
  # set(dialect "-ffree-form -std=f2008 -fimplicit-none")
  set(dialect "-std=f2008")
  set(bounds "-fbounds-check")
  set(warning "-Wall -Wextra")
  set(precision "-fdefault-real-8")
  set(performance "-Ofast")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_COMPILER "ifort")
  set(CMAKE_C_COMPILER "icc")
  set(CMAKE_CXX_COMPILER "icpc")
  # set(dialect "-stand f08 -free -implicitnone")
  set(dialect "-stand f08")
  set(bounds "-check bounds")
  set(warning "-warn all")
  set(precision "-r8")
  set(performance "-Ofast")
endif()
# if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
#   set(CMAKE_Fortran_COMPILER "pgfortran")
#   set(CMAKE_C_COMPILER "pgcc")
#   set(CMAKE_CXX_COMPILER "pgc++")
#   set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
#   set(bounds "-C")
#   set(precision "-r8")
#   set(performance "-fast -Mipa=fast,inline")
# endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${dialect} ${bounds} ${precision}")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} ${dialect} ${precision} ${performance}")

#
# Compile.
#
file(GLOB_RECURSE sources src/*.f src/*.f90 src/*.h)
add_executable(prog ${sources})

# This one is old, the newer version is in project GridAveragedFVM.
