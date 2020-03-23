cmake_minimum_required(VERSION 3.8)
project( transient-heat
         DESCRIPTION "Transient Heat Conduction using MFEM + SUNDIALS + HYPRE"
         LANGUAGES C)

option(ENABLE_MFEM "Enable MFEM" ON)

option(MFEM_DIR "absolute/path/to/mfem@develop/install" [.])
option(SUNDIALS_DIR "absolute/path/to/sundials@develop/install" [.])
option(HYPRE_DIR "absolute/path/to/hypre@2.18.2+superlu-dist+shared/install" [.])

find_package(MPI REQUIRED)
set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

find_library(ZLIB_LIBRARY z
    HINTS ${ZLIB_LIBRARY_DIR} $ENV{METIS_LIBRARY_DIR}
    PATH_SUFFIXES lib
    NO_DEFAULT_PATH
    DOC "The zlib library.")
find_library(ZLIB_LIBRARY z)

find_library(METIS_LIBRARY metis
    HINTS ${METIS_DIR} $ENV{METIS_DIR}
    PATH_SUFFIXES lib
    NO_DEFAULT_PATH
    DOC "The Metis library.")
find_library(METIS_LIBRARY metis)

# Find MFEM
# Find the MFEM header files
find_path(MFEM_INCLUDE_DIRS mfem.hpp
  HINTS ${MFEM_DIR} $ENV{MFEM_DIR}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "Directory with MFEM header.")
find_path(MFEM_INCLUDE_DIRS mfem.hpp)

# Find the MFEM library
find_library(MFEM_LIBRARY mfem
  HINTS ${MFEM_DIR} $ENV{MFEM_DIR}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  DOC "The MFEM library.")
find_library(MFEM_LIBRARY mfem)

find_package(PETSC REQUIRED)

find_package(SUNDIALS REQUIRED)

# Find SUNDIALS
# Find the SUNDIALS header files
find_path(SUNDIALS_INCLUDE_DIRS sundials_version.h
  HINTS ${SUNDIALS_DIR} $ENV{SUNDIALS_DIR}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "Directory with SUNDIALS header.")
find_path(SUNDIALS_INCLUDE_DIRS sundials_version.h)

# Find the SUNDIALS library
find_library(SUNDIALS_LIBRARY sundials_sunmatrixdense
  HINTS ${SUNDIALS_DIR} $ENV{SUNDIALS_DIR}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  DOC "The SUNDIALS library.")
find_library(SUNDIALS_LIBRARY sundials_sunmatrixdense)


# Find HYPRE
# Find the HYPRE header files
find_path(HYPRE_INCLUDE_DIRS HYPRE.h
  HINTS ${HYPRE_DIR} $ENV{HYPRE_DIR}
  PATH_SUFFIXES include
  NO_DEFAULT_PATH
  DOC "Directory with HYPRE header.")
find_path(HYPRE_INCLUDE_DIRS HYPRE.h)

# Find the HYPRE library
find_library(HYPRE_LIBRARY HYPRE
  HINTS ${HYPRE_DIR} $ENV{HYPRE_DIR}
  PATH_SUFFIXES lib
  NO_DEFAULT_PATH
  DOC "The HYPRE library.")
find_library(HYPRE_LIBRARY HYPRE)

add_executable(transient-heat transient-heat.cpp)
target_link_libraries( transient-heat PRIVATE
    ${MFEM_LIBRARY} SUNDIALS::cvode  SUNDIALS::arkode SUNDIALS::kinsol SUNDIALS::nvecparallel SUNDIALS::nvecparhyp
    ${HYPRE_LIBRARY} PETSC::ALL ${METIS_LIBRARY} ${ZLIB_LIBRARY} ${MPI_CXX_LIBRARIES})
target_include_directories(transient-heat PRIVATE
  ${MFEM_INCLUDE_DIRS} ${HYPRE_INCLUDE_DIRS} ${MPI_C_INCLUDE_PATH})

install(TARGETS transient-heat RUNTIME DESTINATION bin)
