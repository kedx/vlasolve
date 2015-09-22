project(VLASOLVE)

set(USE_LONG_INT_DEF false)
set(USE_SIMPLE_PRECISION_DEF false)
set(USE_TR1_DEF true)
set(USE_PTHREADS_DEF true)
set(USE_OPENMP_DEF true)
set(EINSPLINE_DIR_DEF "/home/thierry/libs/")
set(FFTW3_DIR_DEF "/home/thierry/libs/")

SET( CMAKE_CXX_FLAGS_CUSTOM_DEF "-g -O3 -msse2 -mavx -DNDEBUG" CACHE STRING
    "default flags for custom builds"
    FORCE )

#Default build type: Release, Debug or Custom 
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif()

#set (CMAKE_C_COMPILER "gcc-4.4")
#set (CMAKE_CXX_COMPILER "g++-4.4")

set (VLASOLVE_VERSION_MAJOR 0)
set (VLASOLVE_VERSION_MINOR 9)
set (VLASOLVE_VERSION_BUILD 0)

########## NO MODIFICATION NEEDED BELOW THIS LINE ##########

#CMAKE_MINIMUM_REQUIRED(VERSION 2.6.3)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

add_definitions(
  -DVER_MAJOR=${VLASOLVE_VERSION_MAJOR} 
  -DVER_MINOR=${VLASOLVE_VERSION_MINOR} 
  -DVER_BUILD=${VLASOLVE_VERSION_BUILD})

#add_definitions(-Wall)
set(CMAKE_CXX_FLAGS_CUSTOM ${CMAKE_CXX_FLAGS_CUSTOM_DEF} CACHE binary "custom build flags")
set(USE_LONG_INT ${USE_LONG_INT_DEF} CACHE binary "use long integers")
set(USE_SIMPLE_PRECISION ${USE_SIMPLE_PRECISION_DEF} CACHE binary "use simple precision floats")
set(USE_TR1 ${USE_TR1_DEF} CACHE binary "use TR1 C++")
set(USE_PTHREADS ${USE_PTHREADS_DEF} CACHE binary "use pthreads")
set(USE_OPENMP ${USE_OPENMP_DEF} CACHE binary "use OpenMP")
set(EINSPLINE_DIR ${EINSPLINE_DIR_DEF} CACHE PATH "path to Einspline")
set(FFTW3_DIR ${FFTW3_DIR_DEF} CACHE PATH "path to fftw3")
set_property(CACHE USE_PTHREADS PROPERTY STRINGS true false)
set_property(CACHE USE_OPENMP PROPERTY STRINGS true false)

IF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
SET(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}" CACHE PATH "VLASOLVE install prefix" FORCE
)
ENDIF(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

SET(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS  true)
LIST(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/modules")
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6.3)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6.3)
  endif()
endif()

#SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-msse2 -mavx -g -Wa,-adhls -save-temps")
#SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Q -v -DFORTIFY_SOURCE=2 ")
#SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} " -msse2 -mavx")

if (USE_LONG_INT)
  add_definitions(-DUSELONGINT)
endif()

if (USE_SIMPLE_PRECISIONT)
  add_definitions(-DUSESIMPLEPRECISION)
endif()

if (USE_TR1)
  FIND_PACKAGE(TR1)
endif()

if (USE_PTHREADS)
  add_definitions(-DUSE_PTHREADS)
  link_libraries(-lpthread)
endif()

if (USE_OPENMP)
  add_definitions(-fopenmp -DUSE_OPENMP)
  link_libraries(-fopenmp)
endif()


if (TR1_UNORDERED_MAP_FOUND)
  add_definitions(-DHAVE_TR1)
  if (TR1_HEADER_PREFIX) 
    add_definitions(-DHAVE_TR1_HEADER_PREFIX)
  endif()
  if (TR1_NAMESPACE_PREFIX) 
    add_definitions(-DHAVE_TR1_NAMESPACE_PREFIX)
  endif()
endif()

set(NO_EINSPLINE 1)
if(NOT NO_EINSPLINE)
  FIND_PACKAGE(Einspline)
endif()

#set(NO_FFTW3 1)
if(NOT NO_FFTW3)
  FIND_PACKAGE(FFTW3)
endif()

if(NOT NO_MPI)
  FIND_PACKAGE(MPI)
endif()

if(EINSPLINE_FOUND)
  add_definitions(-DHAVE_EINSPLINE)
  include_directories(${EINSPLINE_INCLUDE_DIR})
  link_libraries(${EINSPLINE_LIBRARIES})
endif()

if(FFTW_FOUND)
  include_directories(${FFTW_INCLUDE_DIR} )
  link_libraries(${FFTW_LIBRARIES} )
if(FFTW_THREADS_FOUND)
  add_definitions(-DHAVE_FFTW3)
  include_directories(${FFTW_THREADS_INCLUDE_DIR})
  link_libraries(${FFTW_THREADS_LIBRARIES})
endif()
endif()

if(MPI_CXX_FOUND)
  INCLUDE(CMakeForceCompiler)
  CMAKE_FORCE_CXX_COMPILER(${MPI_CXX_COMPILER} "MPI C++ Compiler")
  add_definitions(-DHAVE_MPI)
  include_directories(${MPI_CXX_INCLUDE_PATH})
  link_libraries(${MPI_CXX_LIBRARIES})
  SET(CMAKE_EXE_LINKER_FLAGS ${CMAKE_EXE_LINKER_FLAGS} ${MPI_CXX_LINK_FLAG})
  SET(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} ${MPI_CXX_COMPILE_FLAGS})
endif()

MACRO(LIST_LENGTH var)
  SET(entries)
  FOREACH(e ${ARGN})
    SET(entries "${entries}.")
  ENDFOREACH(e)
  STRING(LENGTH ${entries} ${var})
ENDMACRO(LIST_LENGTH)

if (CMAKE_COMPILER_IS_GNUCXX)
  execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion
    OUTPUT_VARIABLE GCXX_VERSION)
  string(REGEX MATCHALL "[0-9]+" GCXX_VERSION_COMPONENTS ${GCXX_VERSION})

  LIST_LENGTH(length ${GCXX_VERSION_COMPONENTS})
  
  list(GET GCXX_VERSION_COMPONENTS 0 GCXX_MAJOR)

  if(length GREATER 1)
    list(GET GCXX_VERSION_COMPONENTS 1 GCXX_MINOR)
  else()
    SET(GCXX_MINOR 0)
  endif()
  
  if(length GREATER 2)
    list(GET GCXX_VERSION_COMPONENTS 2 GCXX_BUILD)
  else()
    SET(GCXX_BUILD 0)
  endif()

  message(STATUS "GNU compiler found : g++ v${GCXX_MAJOR}.${GCXX_MINOR}.${GCXX_BUILD}")
endif()

SET(TARGETS_LIST "")
SET(NOT_BUILT_TARGETS_LIST "")
SET(DISABLED_FEATURES_LIST "")

add_subdirectory ("src")

# build a CPack driven installer package
include (InstallRequiredSystemLibraries)
set (CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
set (CPACK_PACKAGE_VERSION_MAJOR "${VLASOLVE_VERSION_MAJOR}")
set (CPACK_PACKAGE_VERSION_MINOR "${VLASOLVE_VERSION_MINOR}")
set (CPACK_PACKAGE_VERSION_PATCH "${VLASOLVE_VERSION_BUILD}")

include (CPack)

message (STATUS "")
message (STATUS "-----------------------------------")
message (STATUS "*** REPORT for Vlasolve v${VLASOLVE_VERSION_MAJOR}.${VLASOLVE_VERSION_MINOR}.${VLASOLVE_VERSION_BUILD} ***")
message (STATUS "-----------------------------------")
message (STATUS "")

if (NOT USE_PTHREADS)
  message(STATUS "   * Threads (pthread) are DISABLED (enable with '-DUSE_PTHREADS=true').") 
else()
   message(STATUS "   * Threads (pthread) are ENABLED (disable with '-DUSE_PTHREADS=false').") 
endif()

if (NOT USE_OPENMP)
  message(STATUS "   * OpenMP support is DISABLED (enable with '-DUSE_OPENMP=true').") 
else()
  message(STATUS "   * OpenMP support is ENABLED (enable with '-DUSE_OPENMP=false').") 
endif()


if (USE_LONG_INT)
  message(STATUS "   * Long integers are ENABLED (disable with '-DUSE_LONG_INT=false')")
  message(STATUS "     It is recommanded to disable unless you are using VERY large datasets.")
else()
  LIST(APPEND DISABLED_FEATURES_LIST all:long_int)
  message(STATUS "   * Long integers are DISABLED (enable with '-DUSE_LONG_INT=true')")
  message(STATUS "     Enable only if you are using VERY large datasets.")
endif()

if (NOT USE_SIMPLE_PRECISION)
  message(STATUS "   * Double precision ENABLED (disable with '-DUSE_SIMPLE_PRECISION=true')")
else()
  LIST(APPEND DISABLED_FEATURES_LIST all:simple_precision)
  message(STATUS "   * Double precision computations DISABLED (enable with '-DUSE_SIMPLE_PRECISION=false')")
endif()

if (USE_TR1)
  if (NOT TR1_UNORDERED_MAP_FOUND)
    LIST(APPEND DISABLED_FEATURES_LIST all:tr1)
    message(STATUS "   * TR1 C++ instructions were not found.")
  else()
    message(STATUS "   * TR1 C++ instructions are ENABLED (disable with '-DUSE_TR1=false').")
  endif()
else()
  LIST(APPEND DISABLED_FEATURES_LIST all:tr1)
  message(STATUS "   * TR1 C++ instructions are DISABLED (enable with '-DUSE_TR1=true').")
endif()

#if (EINSPLINE_FOUND)
#  message(STATUS "   * Einspline was found.")
#else()
#  message(STATUS "   * Einspline was NOT found.")
#endif()

if (FFTW_FOUND) 
  if (FFTW_THREADS_FOUND)
    message(STATUS "   * FFTW3-THREADS was found.")
  else()
    message(STATUS "   * FFTW3-THREADS was NOT found (but FFTW3 was found).")
  endif()
else()
  message(STATUS "   * FFTW3(double) was NOT found.")
endif()

if(MPI_CXX_FOUND)
  message(STATUS "   * MPI was found (using ${MPI_CXX_COMPILER}).")
  message(STATUS "     Disable MPI support with '-DNO_MPI=true'.")
  message(STATUS "     You can also force a specific MPI compiler with:")
  message(STATUS "      -DMPI_CXX_COMPILER=<NAME> (e.g with <NAME> = mpicc)" )
else()
  if(NO_MPI)
    message(STATUS "   * MPI support DISABLED (enable with '-DNO_MPI=false').")
  else()
    message(STATUS "   * MPI support was NOT found.")
  endif()
endif()

message(STATUS "")
message(STATUS " The following binaries will be built:\n     ${TARGETS_LIST}")
if (DEFINED NOT_BUILT_TARGETS_LIST)
  message(STATUS " The following binaries will NOT be built:\n     ${NOT_BUILT_TARGETS_LIST}")
endif()
if (DEFINED DISABLED_FEATURES_LIST)
  message(STATUS " The following features are DISABLED:\n     ${DISABLED_FEATURES_LIST}")
endif()

message(STATUS " ")
message(STATUS " Compiler : ${CMAKE_CXX_COMPILER}")
message(STATUS " Mode : ${CMAKE_BUILD_TYPE}")

message(STATUS "  change compiler with -DCMAKE_CXX_COMPILER=<NAME>.")
message(STATUS "  set compile mode with '-DCMAKE_BUILD_TYPE={Debug,Release,Custom}',")
message(STATUS "  custom flags are set with '-DCMAKE_CXX_FLAGS_CUSTOM=...' ")
message(STATUS " ")
message(STATUS " Vlasolve will be installed in '${CMAKE_INSTALL_PREFIX}'")
message(STATUS " The default install path can be changed with -DCMAKE_INSTALL_PREFIX=<PATH>")


message (STATUS "")
message (STATUS "-----------------------------------")
message (STATUS "***        END OF REPORT        ***")
message (STATUS "-----------------------------------")
message (STATUS "")

