#ifdef GLOBAL_DEFINITION
#define GLOBAL
#else 
#define GLOBAL extern
#endif

GLOBAL int verbose;
GLOBAL int debug;
//GLOBAL int glob_num_threads;
GLOBAL int num_omp_threads;

#undef GLOBAL
