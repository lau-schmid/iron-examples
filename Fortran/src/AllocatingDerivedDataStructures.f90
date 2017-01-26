
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING HEADER FILE PERFOMS THE FOLLOWINGS TASKS. ... !!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 1) PARSE THROUGH THE INPUT FILE TO COUNT AND STORE BASIS, REGIONS ETC. DEFINED... !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ...BY USER IN THE INPUT FILE.  ...!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 2) BASED ON THE COUNTS IT ALLOCATES THE DERIVED DATA STRUCTURE. !!!!!!!!!!!!!!!!!!!!!


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FOLLOWING VARIABLES STORE  NUMBER OF TIMES BASIS , PARTS , REGIONS ETC.  .....!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! .... ARE DEFINED BY THE USER. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  NumberOfEquationsSet= 0
  NumberOfProblem = 0
  NumberOfBasis = 0
  NumberOfSolver = 0
  NumberOfMesh = 0
  NumberOfBoundaryCondition = 0
  NumberOfBasis = 0
  NumberOfMaterialField = 0
  NumberOfEquation = 0
  NumberOfGeometricField = 0
  NumberOfControlLoop = 0
  NumberOfDependentField = 0
  NumberOfRegion = 0
  NumberOfCoordinateSystem = 0
  NumberOfFiberField       = 0
  NumberOfPressureBasis    = 0
  NumberOfField    = 0

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!! PARSE INPUT FILE TO STORE INFORMATION. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(12,file=InputFile ,status="old")
  read(12,'(A)') rdline

  do while (trim(rdline).NE."STOP_PARSING")

       read(12,'(A)') rdline

       CALL searching(rdline,"START_EQUATIONS_SET",NumberOfEquationsSet)
       call searching(rdline,"START_PROBLEM",NumberOfProblem)
       call searching(rdline,"START_BOUNDARY_CONDITIONS",NumberOfBoundaryCondition)
       call searching(rdline,"START_BASIS",NumberOfBasis)
       call searching(rdline,"PRESSURE_BASIS",NumberOfPressureBasis)
       call searching(rdline,"START_MESH",NumberOfMesh)
       call searching(rdline,"START_MATERIAL_FIELD",NumberOfMaterialField)
       call searching(rdline,"START_COORDINATE_SYSTEM",NumberOfCoordinateSystem)
       call searching(rdline,"START_REGION",NumberOfRegion)
       call searching(rdline,"START_GEOMETRIC_FIELD",NumberOfGeometricField)
       call searching(rdline,"START_DEPENDENT_FIELD",NumberOfDependentField)
       call searching(rdline,"START_SOLVER_SETTINGS",NumberOfSolver)
       call searching(rdline,"START_CONTROL_LOOP",NumberOfControlLoop)
       call searching(rdline,"START_FIBER_FIELD",NumberOfFiberField)
       call searching(rdline,"START_PRESSURE_BASIS",NumberOfPressureBasis)
       call searching(rdline,"START_FIELD",NumberOfField)

  enddo
  NumberOfDecomposition = NumberOfMesh               !! This is wrong , though works for this study. Gonna fix it later.
  NumberOfGeneratedMesh = NumberOfMesh    	   !! As per my understanding these two always seem  to be equal.


  !!!!!!!!!!!!!!!!!		ALLOCATE DATA STRUCTURES BASED ON THE PARAMETERS INITIALIZED ABOVE.  !!!!!!!!!!!!!!!!!

  allocate(all_Basis%Basis(NumberOfBasis))
  allocate(all_PressureBasis%PressureBasis(NumberOfPressureBasis))
  allocate(all_BoundaryConditions%BoundaryConditions(NumberOfBoundaryCondition))
  allocate(all_CoordinateSystem%CoordinateSystem(NumberOfCoordinateSystem))
  allocate(all_WorldCoordinateSystem%WorldCoordinateSystem(NumberOfCoordinateSystem))
  allocate(all_Mesh%Mesh(NumberOfMesh))
  allocate(all_Decomposition%Decomposition(NumberOfdecomposition))
  allocate(all_Equations%Equations(NumberOfEquationsSet))
  allocate(all_EquationsSet%EquationsSet(NumberOfEquationsSet))
  allocate(all_GeometricField%GeometricField(NumberOfGeometricField))
  allocate(all_FibreField%FibreField(NumberOfFiberField))
  allocate(all_MaterialField%MaterialField(NumberOfMaterialField))
  allocate(all_EquationsSetField%EquationsSetField(NumberOfEquationsSet))
  allocate(all_Fields%Fields(NumberOfField))
  allocate(all_Problem%Problem(NumberOfProblem))
  allocate(all_Region%Region(NumberOfRegion))
  allocate(all_WorldRegion%WorldRegion(NumberOfWorldRegion))
  allocate(all_Solver%Solver(NumberOfSolver))
  allocate(all_LinearSolver%LinearSolver(NumberOfSolver))
  allocate(all_SolverEquations%SolverEquations(NumberOfSolver))
  allocate(all_ControlLoop%ControlLoop(1))     !! Hard coded , gonna fix it later.
  allocate(all_GeneratedMesh%GeneratedMesh(NumberOfGeneratedMesh))
  allocate(all_DependentField%DependentField(NumberOfDependentField))

