#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <string>
#include <iostream>
#include <fstream>

#include "openmp.h"

#define GLOBAL_DEFINITION
#include "global.h"

#undef GLOBAL_DEFINITION
#include "setup.hxx"

#include "paramsParser.hxx"

int main (int argc, char **argv)
{
  // global variables setup
  num_omp_threads=1;
  verbose=1;
  debug=0;

  mpiComT mpiCom(argc,argv);
  paramsParser params(argc,argv);

  verbose=params.get<>("verbose",paramsParser::defaultCategory(),verbose);
  if (verbose>1) params.report();

  debug=params.get<>("debug",paramsParser::defaultCategory(),debug);
  
#ifdef USE_OPENMP
  int nt_default=1;
#pragma omp parallel
  nt_default=omp_get_num_threads();
  num_omp_threads = params.get<>("threads_per_node",paramsParser::defaultCategory(),nt_default);
  omp_set_num_threads(num_omp_threads);
#endif
  
  std::string reload = params.get<std::string>("reload",paramsParser::defaultCategory(),"");

  solverT solver(mpiCom);
  solver.init(params,reload); 

  params.reportUnused();

  solver.run();
  
  return 0;

}
