
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! SETTING DEFAULT VALUES OF DATA STRUCTURES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  ! default parameters for the solver block
  all_Solver(1)%NewtonRelativeTolerance(1)          = "0"
  all_Solver(1)%NewtonAbsoluteTolerance(1)          = "0"
  all_Solver(1)%IterativeSolverAbsoluteTolerance(1) = "0"
  all_Solver(1)%IterativeSolverRelativeTolerance(1) = "0"
  all_Solver(1)%IterativeSolverAbsoluteTolerance(1) = "0"
  all_Solver(1)%JacobianType(1)                     = "SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED"
  all_Solver(1)%LibraryType(1)                      = "SOLVER_MUMPS_LIBRARY"
  all_Solver(1)%JacobianType(1)                     = "SOLVER_ITERATIVE_NO_PRECONDITIONER"

