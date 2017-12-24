
MODULE PROBLEM_ROUTINES

  USE OpenCMISS
  USE OpenCMISS_Iron
  USE PARAMETERS
  USE OpenCMISS_Variables
  USE REGION_MESH
  USE DECOMPOSITION 
  USE FIELDS_EQUATIONS_SET 
  USE CELLML
  
#ifndef NOMPIMOD
  USE MPI
#endif

#ifdef WIN32
  USE IFQWIN
#endif

  IMPLICIT NONE

#ifdef NOMPIMOD
#include "mpif.h"
#endif

CONTAINS

SUBROUTINE CreateProblem()

  !--------------------------------------------------------------------------------------------------------------------------------
  !Define the problem
  CALL cmfe_Problem_Initialise(Problem,Err)

  IF (ModelType == 0) THEN    ! 3a, "MultiPhysStrain", old tomo mechanics
    IF (SplittingType == 1) THEN    ! Strang splitting
      CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_MULTI_PHYSICS_CLASS, &
        & CMFE_PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE,CMFE_PROBLEM_STRANG_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE],Problem,Err)
    ELSE      ! Godunov splitting
      CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_MULTI_PHYSICS_CLASS, &
        & CMFE_PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE,CMFE_PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ELASTICITY_SUBTYPE],Problem,Err)
    ENDIF
  ELSEIF (ModelType == 1) THEN ! 3, , "MultiPhysStrain", numerically more stable
    IF (SplittingType == 1) THEN    ! Strang splitting
      CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_MULTI_PHYSICS_CLASS, &
        & CMFE_PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE,CMFE_PROBLEM_STRANG_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE],Problem,Err)
    ELSE      ! Godunov splitting
      CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_MULTI_PHYSICS_CLASS, &
        & CMFE_PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE,CMFE_PROBLEM_GUDUNOV_MONODOMAIN_1D3D_ACTIVE_STRAIN_SUBTYPE],Problem,Err)
      ENDIF
  ELSEIF (ModelType == 2) THEN ! 4, "Titin"
    IF (SplittingType == 1) THEN    ! Strang splitting
      CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_MULTI_PHYSICS_CLASS, &
       & CMFE_PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE,CMFE_PROBLEM_STRANG_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE],Problem,Err)
    ELSE      ! Godunov splitting
      CALL cmfe_Problem_CreateStart(ProblemUserNumber,[CMFE_PROBLEM_MULTI_PHYSICS_CLASS, &
       & CMFE_PROBLEM_BIOELECTRIC_FINITE_ELASTICITY_TYPE,CMFE_PROBLEM_GUDUNOV_MONODOMAIN_ELASTICITY_W_TITIN_SUBTYPE],Problem,Err)
    ENDIF
  ENDIF

  CALL cmfe_Problem_CreateFinish(Problem,Err)

END SUBROUTINE CreateProblem


END MODULE PROBLEM_ROUTINES