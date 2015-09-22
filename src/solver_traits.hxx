#ifndef __SOLVER_TRAITS_HXX__
#define __SOLVER_TRAITS_HXX__

#include "all_grids.hxx"
#include "initGen_all.hxx"
#include "grid_handler.hxx"
#include "paramsParser.hxx"
#include "mpiCommunication.hxx"

#include "slicer21.hxx"
#include "vlasovSolver21_spline.hxx"
#include "slicerNN.hxx"
#include "vlasovSolverNN_spline.hxx"

#include "valLocation_traits.hxx"

namespace vlSolver {

  enum type {spline};
  
  template <int P_dims, int U_dims, typename dataType,type solverType>
  struct traits;

  template <typename dataType>
  struct traits<1,2,dataType,spline>
  {  
    //typedef void ERROR_SOLVER_TRAITS_ARE_INVALID;
    typedef dimTraits<1,2> dimTr;
    typedef dataType dataT;
    typedef valLocation_traits<valLocationVal::VERTEX,valLocationVal::VERTEX,valLocationVal::CELL> valLocationTr;
    typedef grid21<dataT,valLocationTr> gridT;
    typedef slicer21<gridT> slicerT;
    typedef gridHandler<gridT,slicerT> gridHandlerT;
    typedef vlasovSolver21_spline<gridHandlerT> solverT;
  };

#ifdef HAVE_FFTW3
  template <typename dataType>
  struct traits<1,1,dataType,spline>
  {  
    //typedef void ERROR_SOLVER_TRAITS_ARE_INVALID;
    typedef dimTraits<1,1> dimTr;
    typedef dataType dataT;
    typedef valLocation_traits<valLocationVal::VERTEX,valLocationVal::VERTEX> valLocationTr;
    typedef gridNN<1,dataT,valLocationTr> gridT;
    typedef slicerNN<gridT> slicerT;
    typedef gridHandler<gridT,slicerT> gridHandlerT;
    typedef vlasovSolverNN_spline<gridHandlerT> solverT;
  };


  template <typename dataType>
  struct traits<2,2,dataType,spline>
  {  
    //typedef void ERROR_SOLVER_TRAITS_ARE_INVALID;
    typedef dimTraits<2,2> dimTr;
    typedef dataType dataT;
    typedef valLocation_traits<valLocationVal::VERTEX,valLocationVal::VERTEX> valLocationTr;
    typedef gridNN<2,dataT,valLocationTr> gridT;
    typedef slicerNN<gridT> slicerT;
    typedef gridHandler<gridT,slicerT> gridHandlerT;
    typedef vlasovSolverNN_spline<gridHandlerT> solverT;
  };

  template <typename dataType>
  struct traits<3,3,dataType,spline>
  {  
    //typedef void ERROR_SOLVER_TRAITS_ARE_INVALID;
    typedef dimTraits<3,3> dimTr;
    typedef dataType dataT;
    typedef valLocation_traits<valLocationVal::VERTEX,valLocationVal::VERTEX> valLocationTr;
    typedef gridNN<3,dataT,valLocationTr> gridT;
    typedef slicerNN<gridT> slicerT;
    typedef gridHandler<gridT,slicerT> gridHandlerT;
    typedef vlasovSolverNN_spline<gridHandlerT> solverT;
  };
#endif

}
#endif
