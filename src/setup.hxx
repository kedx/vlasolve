#ifndef __SETUP_HXX__
#define __SETUP_HXX__

#include "solver_traits.hxx"

// setup the solver here
typedef vlSolver::traits<1,2,FLOAT,vlSolver::spline> solver_traits;
//typedef vlSolver::traits<1,1,FLOAT,vlSolver::spline> solver_traits;
//typedef vlSolver::traits<2,2,FLOAT,vlSolver::spline> solver_traits;
//typedef vlSolver::traits<3,3,FLOAT,vlSolver::spline> solver_traits;


/* NOTHING TO CHANGE BELOW THIS */
// this is a cheap static assert, to check solver_traits validity ;)
// typedef typename solver_traits::ERROR_SOLVER_TRAITS_ARE_INVALID static_check;

typedef solver_traits::solverT solverT;
typedef solverT::mpiComT mpiComT;

#endif
