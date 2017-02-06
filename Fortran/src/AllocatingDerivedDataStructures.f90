
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
  NumberOfSourceField    = 0
  NumberOfOutput   = 0
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!! PARSE INPUT FILE TO STORE INFORMATION. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  OPEN(12,file=InputFile ,status="old")
  READ(12,'(A)') rdline

  DO WHILE (TRIM(rdline).NE."STOP_PARSING")

       READ(12,'(A)') rdline
       CALL searching(rdline,"START_EQUATIONS_SET",NumberOfEquationsSet)
       CALL searching(rdline,"START_PROBLEM",NumberOfProblem)
       CALL searching(rdline,"START_BOUNDARY_CONDITIONS",NumberOfBoundaryCondition)
       CALL searching(rdline,"START_BASIS",NumberOfBasis)
       CALL searching(rdline,"PRESSURE_BASIS",NumberOfPressureBasis)
       CALL searching(rdline,"START_MESH",NumberOfMesh)
       CALL searching(rdline,"START_MATERIAL_FIELD",NumberOfMaterialField)
       CALL searching(rdline,"START_COORDINATE_SYSTEM",NumberOfCoordinateSystem)
       CALL searching(rdline,"START_REGION",NumberOfRegion)
       CALL searching(rdline,"START_GEOMETRIC_FIELD",NumberOfGeometricField)
       CALL searching(rdline,"START_DEPENDENT_FIELD",NumberOfDependentField)
       CALL searching(rdline,"START_SOLVER_SETTINGS",NumberOfSolver)
       CALL searching(rdline,"START_CONTROL_LOOP",NumberOfControlLoop)
       CALL searching(rdline,"START_FIBER_FIELD",NumberOfFiberField)
       CALL searching(rdline,"START_PRESSURE_BASIS",NumberOfPressureBasis)
       CALL searching(rdline,"START_SOURCE_FIELD",NumberOfSourceField)
       CALL searching(rdline,"START_OUTPUT",NumberOfOutput)
  END DO
  NumberOfDecomposition = NumberOfMesh               !! PLEASE IGNORE IT, I WILL FIX IT AFTER OUR MEETING
  NumberOfGeneratedMesh = NumberOfMesh    	     !! PLEASE IGNORE IT, I WILL FIX IT AFTER OUR MEETING


  !!!!!!!!!!!!!!!!!		ALLOCATE DATA STRUCTURES BASED ON THE PARAMETERS INITIALIZED ABOVE.  !!!!!!!!!!!!!!!!!

  ALLOCATE(all_Basis%Basis(NumberOfBasis))
  ALLOCATE(all_PressureBasis%PressureBasis(NumberOfPressureBasis))
  ALLOCATE(all_BoundaryConditions%BoundaryConditions(NumberOfBoundaryCondition))
  ALLOCATE(all_CoordinateSystem%CoordinateSystem(NumberOfCoordinateSystem))
  ALLOCATE(all_WorldCoordinateSystem%WorldCoordinateSystem(NumberOfCoordinateSystem))
  ALLOCATE(all_Mesh%Mesh(NumberOfMesh))
  ALLOCATE(all_Decomposition%Decomposition(NumberOfdecomposition))
  ALLOCATE(all_Equations%Equations(NumberOfEquationsSet))
  ALLOCATE(all_EquationsSet%EquationsSet(NumberOfEquationsSet))
  ALLOCATE(all_GeometricField%GeometricField(NumberOfGeometricField))
  ALLOCATE(all_FibreField%FibreField(NumberOfFiberField))
  ALLOCATE(all_MaterialField%MaterialField(NumberOfMaterialField))
  ALLOCATE(all_EquationsSetField%EquationsSetField(NumberOfEquationsSet))
  ALLOCATE(all_Fields%Fields(NumberOfSourceField))
  ALLOCATE(all_Problem%Problem(NumberOfProblem))
  ALLOCATE(all_Region%Region(NumberOfRegion))
  ALLOCATE(all_WorldRegion%WorldRegion(NumberOfWorldRegion))
  ALLOCATE(all_Solver%Solver(NumberOfSolver))
  ALLOCATE(all_LinearSolver%LinearSolver(NumberOfSolver))
  ALLOCATE(all_SolverEquations%SolverEquations(NumberOfSolver))
  ALLOCATE(all_ControlLoop%ControlLoop(1))     !! PLEASE IGNORE IT, I WILL FIX IT AFTER OUR MEETING
  ALLOCATE(all_GeneratedMesh%GeneratedMesh(NumberOfGeneratedMesh))
  ALLOCATE(all_DependentField%DependentField(NumberOfDependentField))
  ALLOCATE(all_SourceField%SourceField(NumberOfSourceField))

