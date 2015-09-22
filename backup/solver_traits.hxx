#ifndef __SOLVER_TRAITS_HXX__
#define __SOLVER_TRAITS_HXX__

#include "all_grids.hxx"
#include "initGen_all.hxx"
#include "grid_handler.hxx"
#include "paramsParser.hxx"
#include "mpiCom.hxx"
#include "slicer21Traits.hxx"
#include "vlasovSolver21_spline.hxx"

namespace vlSolver {

  enum type {spline21};
    
  template <int P_dims, int U_dims, typename dataType, type solverType>
  struct traits
  {
   
  };

  template <1,2,typename dataType,spline21>
  struct traits
  {  
    typedef dimtraits<1,2> dimTr;
    typedef dataType dataT;

    typedef mpiComT mpiT;
    typedef paramsParserT parserT;
    typedef scale<dataT> scaleT;

    typedef grid21<dataT,scaleT> gridT;
    typedef slicer21Traits<gridT> slicerT;
    typedef initGenAll<gridT> initGenAllT;
    typedef gridHandler<gridT,mpiT,slicerT> gridHandlerT;

    typedef vlasovSolver21_spline<gridHandlerT,initGenAllT> solverT;
  } 

}
#endif
