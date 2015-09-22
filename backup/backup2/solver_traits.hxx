#ifndef __SOLVER_TRAITS_HXX__
#define __SOLVER_TRAITS_HXX__

#include "all_grids.hxx"
#include "initGen_all.hxx"
#include "grid_handler.hxx"
#include "paramsParser.hxx"
#include "mpiCommunication.hxx"
#include "slicer21.hxx"
#include "vlasovSolver21_spline.hxx"
#include "valLocationType.hxx"

namespace vlSolver {

  enum type {spline21};
  
  template <int P_dims, int U_dims, typename dataType,type solverType>
  struct traits;

  template <typename dataType>
  struct traits<1,2,dataType,spline21>
  {  
    //typedef void ERROR_SOLVER_TRAITS_ARE_INVALID;
    typedef dimTraits<1,2> dimTr;
    typedef dataType dataT;
    
    typedef grid21<dataT,valLocationVal::VERTEX> gridT;
    typedef slicer21<gridT> slicerT;
    typedef gridHandler<gridT,slicerT> gridHandlerT;
    typedef vlasovSolver21_spline<gridHandlerT> solverT;
  };

}
#endif
