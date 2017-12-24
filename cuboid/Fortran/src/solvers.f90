
MODULE SOLVERS

  USE OpenCMISS
  USE OpenCMISS_Iron
  USE PARAMETERS
  USE OpenCMISS_Variables
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

SUBROUTINE CreateSolvers()
  TYPE(cmfe_SolverType) :: linearSolver

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem solvers
  ! create parabolic and ode solver, create subsolver structure
  CALL cmfe_Problem_SolversCreateStart(Problem,Err)

  !Create the DAE solver
  CALL cmfe_Solver_Initialise(SolverDAE,Err)
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopMonodomainNumber,CMFE_CONTROL_LOOP_NODE], &
   & SolverDAEIndex,SolverDAE,Err)
  CALL cmfe_Solver_DAETimeStepSet(SolverDAE,ODETimeStep,Err)  ! #timestepset
  
  ! might be handy some day:
  ! CALL cmfe_Solver_DAETimesSet(SolverDAE,?0.0_CMISSRP?,?0.001_CMISSRP?,Err)
  
  SELECT CASE(ODESolverId) !use bdf instead of default-explicitEuler
    CASE(2) ! BDF
      CALL cmfe_Solver_DAESolverTypeSet(SolverDAE,CMFE_SOLVER_DAE_BDF,Err)
      CALL cmfe_Solver_DAEbdfSetTolerance(SolverDAE,0.0000001_CMISSRP,0.0000001_CMISSRP,err) !default values were both: 1.0E-7
    CASE(3) ! GL, not stable yet
      CALL cmfe_Solver_DAESolverTypeSet(SolverDAE,CMFE_SOLVER_DAE_GL,Err)
    CASE(4) ! Crank-Nicolson, not stable yet
      CALL cmfe_Solver_DAESolverTypeSet(SolverDAE,CMFE_SOLVER_DAE_CRANK_NICOLSON,Err)
    CASE(5) ! improved Euler (Heun)
      CALL cmfe_Solver_DAEEulerSolverTypeSet(SolverDAE,CMFE_SOLVER_DAE_EULER_IMPROVED,Err)
    CASE DEFAULT
      IF (ComputationalNodeNumber == 0) THEN
        PRINT *, ""
        PRINT *, "Warning: For the DAE Problem (0D) the standard explicit Euler method is used, which is slow. " // NEW_LINE('A') &
          & // "          Consider to use another DAE solver via the command line option 'ODESolverId'."
        PRINT *, ""
      ENDIF
  END SELECT
  !-------------------------------------------------------------------------------------------  
  !Set the Number of ODE time steps. CARE: This makes cmfe_Solver_DAETimeStepSet() obsolete!
  IF(ODESolverId==1 .AND. OdeNSteps/=-1) THEN
    CALL cmfe_Solver_DAEEulerForwardSetNSteps(SolverDAE,OdeNSteps,Err)
  ELSEIF(ODESolverId==5 .AND. OdeNSteps/=-1) THEN
    CALL cmfe_Solver_DAEEulerImprovedSetNSteps(SolverDAE,OdeNSteps,Err)
  END IF
  
  !> \todo or not-todo... solve the CellML equations on the GPU for efficiency (later)
  !CALL cmfe_Solver_DAESolverTypeSet(SolverDAE,CMFE_SOLVER_DAE_EXTERNAL,Err)

  IF (DebuggingOutput) THEN
    CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_SOLVER_OUTPUT,Err)
    !CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
    !CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_TIMING_OUTPUT,Err)
    !CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_SOLVER_OUTPUT,Err)
    !CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_MATRIX_OUTPUT,Err)
  ELSE
    CALL cmfe_Solver_OutputTypeSet(SolverDAE,CMFE_SOLVER_NO_OUTPUT,Err)
  ENDIF

  !-----------------------------------------------------------------------------------------------------
  !Create the parabolic solver
  CALL cmfe_Solver_Initialise(SolverParabolic,Err)
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopMonodomainNumber,CMFE_CONTROL_LOOP_NODE], &
   & SolverParabolicIndex,SolverParabolic,Err)
  

  IF (SplittingType == 1) THEN    ! Strang splitting      
    CALL cmfe_Solver_DynamicSchemeSet(SolverParabolic,CMFE_SOLVER_DYNAMIC_CRANK_NICOLSON_SCHEME,Err)
  ELSE    ! Godunov splitting
    CALL cmfe_Solver_DynamicSchemeSet(SolverParabolic,CMFE_SOLVER_DYNAMIC_BACKWARD_EULER_SCHEME,Err)
  ENDIF
  
  ! data structure is as follows:
  ! SolverParabolic
  !   LinearSolver
  !      IterativeSolver
  !   LinkedSolver(1) -> LinearSolver
  ! Retrieve linear solver
  NULLIFY(linearSolver%solver)
  CALL cmfe_Solver_DynamicLinearSolverGet(SolverParabolic, linearSolver, Err)

  ! direct:
  ! CALL cmfe_Solver_LinearTypeSet(linearSolver, CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE, Err)
  ! CALL cmfe_Solver_LinearDirectTypeSet(linearSolver, CMFE_SOLVER_ITERATIVE_CONJUGATE_GRADIENT, Err)
  ! CMFE_SOLVER_DIRECT_LU
  ! CMFE_SOLVER_DIRECT_CHOLESKY
  ! CMFE_SOLVER_DIRECT_SVD
  !
  ! iterative:
  ! CALL cmfe_Solver_LinearTypeSet(linearSolver, CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE, Err)
  ! CALL cmfe_Solver_LinearIterativeTypeSet(linearSolver, CMFE_SOLVER_ITERATIVE_CONJUGATE_GRADIENT, Err)
  ! CALL cmfe_Solver_LinearIterativePreconditionerTypeSet(linearSolver, CMFE_SOLVER_ITERATIVE_INCOMPLETE_LU_PRECONDITIONER, Err)
  ! solver:
  ! CMFE_SOLVER_ITERATIVE_GMRES
  ! CMFE_SOLVER_ITERATIVE_CONJUGATE_GRADIENT
  ! CMFE_SOLVER_ITERATIVE_CONJGRAD_SQUARED
  ! preconditioner:
  ! CMFE_SOLVER_ITERATIVE_NO_PRECONDITIONER
  ! CMFE_SOLVER_ITERATIVE_JACOBI_PRECONDITIONER
  ! CMFE_SOLVER_ITERATIVE_BLOCK_JACOBI_PRECONDITIONER
  ! CMFE_SOLVER_ITERATIVE_SOR_PRECONDITIONER
  ! CMFE_SOLVER_ITERATIVE_INCOMPLETE_CHOLESKY_PRECONDITIONER
  ! CMFE_SOLVER_ITERATIVE_INCOMPLETE_LU_PRECONDITIONER
  ! CMFE_SOLVER_ITERATIVE_ADDITIVE_SCHWARZ_PRECONDITIONER
  
  !MonodomainSolverId, MonodomainPreconditionerId, ODESolverId
  IF (MonodomainSolverId <= 1) THEN    ! direct solver
    
    CALL cmfe_Solver_LinearTypeSet(linearSolver, CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE, Err)
  
    SELECT CASE(MonodomainSolverId)
      CASE(1)
        CALL cmfe_Solver_LinearDirectTypeSet(linearSolver, CMFE_SOLVER_DIRECT_LU, Err)
      CASE DEFAULT
        PRINT*, "Warning: Wrong MonodomainSolverId ",MonodomainSolverId
        CALL cmfe_Solver_LinearDirectTypeSet(linearSolver, CMFE_SOLVER_DIRECT_LU, Err)
    END SELECT
  ELSE ! iterative solver
    
    CALL cmfe_Solver_LinearTypeSet(linearSolver, CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE, Err)
    
    SELECT CASE(MonodomainSolverId)
      CASE(2)
        CALL cmfe_Solver_LinearIterativeTypeSet(linearSolver, CMFE_SOLVER_ITERATIVE_GMRES, Err)
      CASE(3)
        CALL cmfe_Solver_LinearIterativeTypeSet(linearSolver, CMFE_SOLVER_ITERATIVE_CONJUGATE_GRADIENT, Err)
      CASE(4)
        CALL cmfe_Solver_LinearIterativeTypeSet(linearSolver, CMFE_SOLVER_ITERATIVE_CONJGRAD_SQUARED, Err)
      CASE DEFAULT
        PRINT*, "Warning: Wrong MonodomainSolverId ",MonodomainSolverId
        CALL cmfe_Solver_LinearIterativeTypeSet(linearSolver, CMFE_SOLVER_ITERATIVE_GMRES, Err)
    END SELECT
    
    SELECT CASE(MonodomainPreconditionerId)
      CASE(1)
        CALL cmfe_Solver_LinearIterativePreconditionerTypeSet(linearSolver, CMFE_SOLVER_ITERATIVE_NO_PRECONDITIONER, Err)
      CASE(2)
        CALL cmfe_Solver_LinearIterativePreconditionerTypeSet(linearSolver, CMFE_SOLVER_ITERATIVE_JACOBI_PRECONDITIONER, Err)
      CASE(3)
        CALL cmfe_Solver_LinearIterativePreconditionerTypeSet(linearSolver, CMFE_SOLVER_ITERATIVE_BLOCK_JACOBI_PRECONDITIONER, Err)
      CASE(4)
        CALL cmfe_Solver_LinearIterativePreconditionerTypeSet(linearSolver, CMFE_SOLVER_ITERATIVE_SOR_PRECONDITIONER, Err)
      CASE(5)
        CALL cmfe_Solver_LinearIterativePreconditionerTypeSet(linearSolver, &
         & CMFE_SOLVER_ITERATIVE_INCOMPLETE_CHOLESKY_PRECONDITIONER, Err)
      CASE(6)
        CALL cmfe_Solver_LinearIterativePreconditionerTypeSet(linearSolver, CMFE_SOLVER_ITERATIVE_INCOMPLETE_LU_PRECONDITIONER, Err)
      CASE(7)
        CALL cmfe_Solver_LinearIterativePreconditionerTypeSet(linearSolver, & 
         & CMFE_SOLVER_ITERATIVE_ADDITIVE_SCHWARZ_PRECONDITIONER, Err)
      CASE DEFAULT
        PRINT *, "Warning: Wrong MonodomainPreconditionerId, ", MonodomainPreconditionerId
        CALL cmfe_Solver_LinearIterativePreconditionerTypeSet(linearSolver, CMFE_SOLVER_ITERATIVE_NO_PRECONDITIONER, Err)
    END SELECT
    
  ENDIF
  
  !IF (DEBUGGING_PROBLEM_OUTPUT) THEN
  !  PRINT*, ""
  !  PRINT*, "before cmfe_Solver_LinearTypeSet"
  !  CALL cmfe_PrintSolver(SolverParabolic, 5, 10, Err)
  !ENDIF
  
  ! Recreate linearSolver subtype as direct solver (instead of iterative solver)
  !CALL cmfe_Solver_LinearTypeSet(linearSolver, CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE, Err)
  
  !IF (DEBUGGING_PROBLEM_OUTPUT) THEN
  !  PRINT*, ""
  !  PRINT*, "========================================================================="
  !  PRINT*, "After cmfe_Solver_LinearTypeSet"
  !  CALL cmfe_PrintSolver(SolverParabolic, 5, 10, Err)
  !ENDIF
  
  !SOLVER_LINEAR_DIRECT_TYPE_SET
  
  !CALL cmfe_SOLVER_DYNAMICLINEARITYTYPESET(SolverParabolic, CMFE_SOLVER_DYNAMIC_LINEAR, Err)
  !CALL cmfe_Solver_NonlinearTypeSet(SolverParabolic,CMFE_SOLVER_NONLINEAR_NEWTON,Err)

  IF (DebuggingOutput) THEN
    CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_NO_OUTPUT,Err)
    !CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
    !CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_TIMING_OUTPUT,Err)
    !CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_SOLVER_OUTPUT,Err)
    !CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_MATRIX_OUTPUT,Err)
  ELSE
    CALL cmfe_Solver_OutputTypeSet(SolverParabolic,CMFE_SOLVER_NO_OUTPUT,Err)
  ENDIF

  !Create the Finte Elasticity solver
  CALL cmfe_Solver_Initialise(SolverFE,Err)
  CALL cmfe_Solver_Initialise(LinearSolverFE,Err)
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopElasticityNumber,CMFE_CONTROL_LOOP_NODE], &
   & SolverFEIndex,SolverFE,Err)

  ! only enable output for specified ComputeNode
  IF (DebuggingOutput) THEN
    IF (ComputationalNodeNumber == 0) THEN
      !CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_NO_OUTPUT,Err)
      !CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_PROGRESS_OUTPUT,Err)
      !CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_TIMING_OUTPUT,Err)
      CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_SOLVER_OUTPUT,Err)
      !CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_MATRIX_OUTPUT,Err)
    ELSE
      CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_NO_OUTPUT,Err)
    ENDIF
  ELSE
    CALL cmfe_Solver_OutputTypeSet(SolverFE,CMFE_SOLVER_NO_OUTPUT,Err)
  ENDIF

  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(SolverFE,CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED,Err)
!  CALL cmfe_Solver_NewtonJacobianCalculationTypeSet(SolverFE,CMFE_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED,Err) ! 721 steps
  CALL cmfe_Solver_NewtonMaximumIterationsSet(SolverFE,NewtonMaximumNumberOfIterations,Err)
  CALL cmfe_Solver_NewtonAbsoluteToleranceSet(SolverFE,NewtonTolerance,Err)
  CALL cmfe_Solver_NewtonSolutionToleranceSet(SolverFE,NewtonTolerance,Err)
  CALL cmfe_Solver_NewtonRelativeToleranceSet(SolverFE,NewtonTolerance,Err)
  CALL cmfe_Solver_NewtonLinearSolverGet(SolverFE,LinearSolverFE,Err)
  CALL cmfe_Solver_LinearTypeSet(LinearSolverFE,CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE,Err)
!  CALL cmfe_Solver_LinearTypeSet(LinearSolverFE,CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE,Err)

  CALL cmfe_Problem_SolversCreateFinish(Problem,Err)
  
  

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem solver CellML equations
  CALL cmfe_Problem_CellMLEquationsCreateStart(Problem,Err)

  CALL cmfe_Solver_Initialise(SolverDAE,Err)
  CALL cmfe_CellMLEquations_Initialise(CellMLEquations,Err)
  
  !                           problem, control loop indetifiers (in),                    solverindex (in), solver (out)
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopMonodomainNumber,CMFE_CONTROL_LOOP_NODE],SolverDAEIndex,SolverDAE,Err)
  CALL cmfe_Solver_CellMLEquationsGet(SolverDAE,CellMLEquations,Err)
  
  !                                   in              in                out
  CALL cmfe_CellMLEquations_CellMLAdd(CellMLEquations,CellMLEnvironment,CellMLIndex,Err)
    
  CALL cmfe_Problem_CellMLEquationsCreateFinish(Problem,Err)

  !--------------------------------------------------------------------------------------------------------------------------------
  !Create the problem solver equations
  CALL cmfe_Problem_SolverEquationsCreateStart(Problem,Err)

  !Create the problem solver parabolic equations (Monodomain)
  CALL cmfe_Solver_Initialise(SolverParabolic,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquationsM,Err)
  
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopMonodomainNumber,CMFE_CONTROL_LOOP_NODE], &
   & SolverParabolicIndex,SolverParabolic,Err)
   
   
  CALL cmfe_Solver_SolverEquationsGet(SolverParabolic,SolverEquationsM,Err)
  
  
  !PRINT*, ""
  !PRINT*, "After cmfe_Solver_SolverEquationsGet"
  !CALL cmfe_PrintSolverEquationsM(SolverEquationsM, 5, 10, Err)
  
  
  
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsM,CMFE_SOLVER_SPARSE_MATRICES,Err)
  
  
  !PRINT*, ""
  !PRINT*, "After cmfe_SolverEquations_SparsityTypeSet"
  !CALL cmfe_PrintSolverEquationsM(SolverEquationsM, 5, 10, Err)
  
  !CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsM,CMFE_SOLVER_FULL_MATRICES,Err)
  
  
  ! TYPE(cmfe_SolverEquationsType), INTENT(IN) :: solverEquations !<The solver equations to add the equations set for.
  ! TYPE(cmfe_EquationsSetType), INTENT(IN) :: equationsSet !<The equations set to add.
  ! INTEGER(INTG), INTENT(OUT) :: equationsSetIndex !<On return, the index of the added equations set in the solver equations
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquationsM,EquationsSetM,EquationsSetIndexM,Err)
  
  !PRINT*, ""
  !PRINT*, "cmfe_PrintSolverEquationsType:"
  !CALL cmfe_PrintSolverEquations(SolverEquationsM, 5, 10, Err)
  !CALL cmfe_PrintSolverEquationsM(SolverEquationsM, Err)
  !STOP

  !Create the problem solver Finite Elasticity equations
  CALL cmfe_Solver_Initialise(SolverFE,Err)
  CALL cmfe_SolverEquations_Initialise(SolverEquationsFE,Err)
  CALL cmfe_Problem_SolverGet(Problem,[ControlLoopElasticityNumber,CMFE_CONTROL_LOOP_NODE],SolverFEIndex,SolverFE,Err)
  CALL cmfe_Solver_SolverEquationsGet(SolverFE,SolverEquationsFE,Err)
  CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsFE,CMFE_SOLVER_SPARSE_MATRICES,Err)
  !CALL cmfe_SolverEquations_SparsityTypeSet(SolverEquationsFE,CMFE_SOLVER_FULL_MATRICES,Err)
  CALL cmfe_SolverEquations_EquationsSetAdd(SolverEquationsFE,EquationsSetFE,EquationsSetIndexFE,Err)

  CALL cmfe_Problem_SolverEquationsCreateFinish(Problem,Err)

  IF (DEBUGGING_PROBLEM_OUTPUT) THEN
    PRINT*, "==========================="
    PRINT*, "SolverFE"
    CALL cmfe_PrintSolver(SolverFE, 6, 3, Err)
  ENDIF
  
  !SolverFE%solver%NONLINEAR_SOLVER%NEWTON_SOLVER%LINEAR_SOLVER%OUTPUT_TYPE = SOLVER_MATRIX_OUTPUT
  
  !STOP
  
  IF (DEBUGGING_PROBLEM_OUTPUT) THEN
    PRINT*, ""
    PRINT*, "after_cmfe_Problem_SolverEquationsCreateFinish"
    CALL cmfe_PrintSolver(SolverParabolic, 5, 10, Err)
  ENDIF
  
END SUBROUTINE CreateSolvers

END MODULE SOLVERS