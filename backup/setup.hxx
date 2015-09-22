#ifndef __SETUP_HXX__
#define __SETUP_HXX__

#include "initGen_all.hxx"
#include "vlasovSolver.hxx"
#include "dimTraits.hxx"
#include "grid.hxx"
#include "interpol.hxx"

typedef double dataT;
typedef dimTraits<1,2> dimTr; 
typedef regularGrid<dimTr,dataT> gridT;
typedef interpol<interpType::SPLINE,dataT> interpolT;

/* NOTHING TO CHANGE BELOW THIS */
typedef initGenAll<gridT> initGenAllT;
typedef initGenAllT::initGenT initGenT;
typedef vlSolver<gridT,interpolT> solverT;

#endif
