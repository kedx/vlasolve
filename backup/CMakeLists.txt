
SET( VLASOLVE_SRC
  vlasolve.cpp   
  initGen.hxx 
  initGen_all.hxx 
  initGen_UniformDensitySphere.hxx 
  initGen_Plummer.hxx 
  initGen_Hernquist.hxx   
  interpol.hxx 
  interpol_DUMMY.hxx 
  interpol_spline.hxx 
  interpol_PFC2.hxx 
  setup.hxx 
  scale.hxx
  grid.hxx 
  vlasovSolver.hxx
  solverParams.hxx 
  snapshots.hxx
  NDfield.c 
  NDfield.h 
  myendian.c 
  myendian.h 
  types.h 
  )

add_executable(vlasolve ${VLASOLVE_SRC})
  
SET(TARGETS_LIST_LOC ${TARGETS_LIST_LOC} vlasolve)
install (TARGETS vlasolve DESTINATION bin)




SET(TARGETS_LIST ${TARGETS_LIST} ${TARGETS_LIST_LOC} PARENT_SCOPE)
SET(NOT_BUILT_TARGETS_LIST ${NOT_BUILT_TARGETS_LIST} ${NOT_BUILT_TARGETS_LIST_LOC} PARENT_SCOPE)

SET(DISABLED_FEATURES_LIST ${DISABLED_FEATURES_LIST} ${DISABLED_FEATURES_LIST_LOC} PARENT_SCOPE)
