#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <string>
#include <iostream>
#include <fstream>
//#include <boost/program_options.hpp>

#include "setup.hxx"
#include "initGen_all.hxx"

#define GLOBAL_DEFINITION
#include "global.h"

using namespace std;
//namespace po = boost::program_options;

int main (int argc, char **argv)
{
  glob_num_threads=1;
  glob_num_omp_threads=1;
  verbose=0;
  debug=0;
  /*
  // OPTIONS LIST GOES HERE
  // allowed only in command line
  po::options_description generic("Command line only options");  
  generic.add_options()
    ("help,h", "produce help message")
    ("config-file,c", po::value< string >(), "configuration file")
    ;

  // allowed in cmd line and config file 
  po::options_description config("Configuration");
  config.add_options()
    ("verbose,v", po::value<int>()->implicit_value(1), "set verbose level")
#if defined (USE_OPENMP) || defined (USE_THREADS)
    ("num-threads,t", po::value<int>(), "the initial number of threads")
#endif
    ("init-file,i", po::value<string>(), "initial conditions filename")
    ;

  // hidden options
  po::options_description hidden("Hidden");
  hidden.add_options();
    

  po::positional_options_description p;
  p.add("config-file", -1);
  
  // END OF OPTIONS LIST
  po::options_description cmdline_options("Command line only options");
  cmdline_options.add(generic).add(config).add(hidden);
  po::options_description config_file_options("General options");
  config_file_options.add(config).add(hidden);
  po::options_description visible_options("Options");
  visible_options.add(generic).add(config);
  
  
  po::variables_map vm;
  try {
    po::store(po::command_line_parser(argc, argv).options(cmdline_options).
	      positional(p).run(), vm);  
    if (vm.count("config-file"))
      {
	ifstream in(vm["config-file"].as<string>().c_str());
	po::store(po::parse_config_file(in, config_file_options),vm);
      }
  }
  
  catch(std::exception &err) {
    cout << "\n" << err.what() << "\n\n";
    cout << visible_options << "\n";
    exit(-1);
    //cout<< err << "\n";
  }
  catch(...) { 
    cerr << "An unknown exception occurs while parsing command line arguments." << "\n"; 
    exit(-1); 
  }
  po::notify(vm);
  
  if (vm.count("help")) {
    cout << visible_options << "\n";
    return 1;
  }
  
  if (vm.count("num-threads")) {
#ifdef USE_OPENMP
    omp_set_num_threads(vm["num-threads"].as<int>());
#pragma omp parallel
    glob_num_omp_threads=omp_get_num_threads();
#pragma omp end parallel
#else
    fprintf(stderr,"WARNING: option '-nThreads' ignored.\n");
    fprintf(stderr,"   Add command 'set (USE_THREADS 1)' in CMakeList.txt and recompile\n");
    fprintf(stderr,"   in order to enable openMP/pthreads support.\n");
#endif
  }
#ifdef USE_THREADS
#ifdef USE_OPENMP
  glob_num_threads=glob_num_omp_threads;
#else
  glob_num_threads=vm["num-threads"].as<int>();
#endif
#endif
  */

  //initT *init=new UniformDensitySphere<gridT>();
  
  //initT *init=new Plummer<gridT>();

  std::string iniTypeStr("UDF");
  if (argc>1) iniTypeStr = std::string(argv[1]);

  initGenT *init = initGenAllT::get(iniTypeStr);
  solverT solver(init);
  solver.solve(); 
}

/*
initialConditions<grid_traits<1,2,double>,UniformDensitySphere>

template <class T> 
initialConditions<grid_traits<1,2,T>,UniformDensitySphere>
{
  UniformDensitySphere()
}
*/
