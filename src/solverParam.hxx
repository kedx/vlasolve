#ifndef __VLASOV_SOLVER_PARAMS_HXX__
#define __VLASOV_SOLVER_PARAMS_HXX__

#include "myIO.hxx"

struct solverParam {
  double G;
  double T0;
  double Tmax;
  double dt;
  int nSpecies;

  long snapshotEvery;
  long restartEvery;
  long statsEvery;

  void write(FILE *f)
    {
      fwrite(&G,sizeof(double),1,f);
      fwrite(&T0,sizeof(double),1,f);
      fwrite(&Tmax,sizeof(double),1,f);
      fwrite(&dt,sizeof(double),1,f);
      fwrite(&nSpecies,sizeof(int),1,f);

      fwrite(&snapshotEvery,sizeof(long),1,f);
      fwrite(&restartEvery,sizeof(long),1,f);    
      fwrite(&statsEvery,sizeof(long),1,f);    
    }

  void read(FILE *f, bool swap)
    {
      myIO::fread(&G,sizeof(double),1,f,swap);
      myIO::fread(&T0,sizeof(double),1,f,swap);
      myIO::fread(&Tmax,sizeof(double),1,f,swap);
      myIO::fread(&dt,sizeof(double),1,f,swap);
      myIO::fread(&nSpecies,sizeof(int),1,f,swap);

      myIO::fread(&snapshotEvery,sizeof(long),1,f,swap);
      myIO::fread(&restartEvery,sizeof(long),1,f,swap);
      myIO::fread(&statsEvery,sizeof(long),1,f,swap);
    }
};

#endif
