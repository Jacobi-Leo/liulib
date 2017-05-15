cmake_minimum_required(VERSION 3.5)

project(matrix)
enable_language(Fortran)

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -llapack")
set(CMAKE_Fortran_COMPILER_ID "Intel")

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_COMPILER "gfortran")
  set(dialect "-ffree-form -std=f2008 -fimplicit-none")
  set(bounds "-fbounds-check")
  set(warning "-Wall -Wextra")
  set(precision "-fdefault-real-8")
  set(performance "-Ofast")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_COMPILER "ifort")
  set(dialect "-stand f08 -free -implicitnone")
  set(bounds "-check bounds")
  set(warning "-warn")
  set(precision "-r8")
  set(performance "-Ofast")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
  set(CMAKE_Fortran_COMPILER "pgfortran")
  set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
  set(bounds "-C")
  set(precision "-r8")
  set(performance "-fast -Mipa=fast,inline")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${dialect} ${bounds} ${precision}")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect} ${precision} ${performance}")

#
# Compile.
#
file(GLOB_RECURSE sources  src/test.f90 src/denseMatrix.f90 src/constant.f90 src/*.h)
add_executable(prog ${sources})
