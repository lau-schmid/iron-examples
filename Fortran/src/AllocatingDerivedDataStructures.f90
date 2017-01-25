
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING HEADER FILE PERFOMS THE FOLLOWINGS TASKS. ... !!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 1) PARSE THROUGH THE INPUT FILE TO COUNT AND STORE BASIS, REGIONS ETC. DEFINED... !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ...BY USER IN THE INPUT FILE.  ...!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 2) BASED ON THE COUNTS IT ALLOCATES THE DERIVED DATA STRUCTURE. !!!!!!!!!!!!!!!!!!!!!

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! FOLLOWING VARIABLES STORE  NUMBER OF TIMES BASIS , PARTS , REGIONS ETC.  .....!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! .... ARE DEFINED BY THE USER. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  num_of_EquationsSet= 0
  num_of_Problem = 0
  num_of_Basis = 0
  num_of_Solver = 0
  num_of_Mesh = 0
  num_of_BoundaryCondition = 0
  num_of_Basis = 0
  num_of_PressureBasis = 0
  num_of_MaterialField = 0
  num_of_Equation = 0
  num_of_GeometricField = 0
  num_of_ControlLoop = 0
  num_of_DependentField = 0
  num_of_Region = 0
  num_of_CoordinateSystem = 0
  num_of_FiberField       = 0
  num_of_PressureBasis    = 0
  num_of_Field    = 0
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!! PARSE INPUT FILE TO STORE INFORMATION. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  open(12,file=InputFile ,status="old")
  read(12,'(A)') rdline

  do while (trim(rdline).NE."STOP_PARSING")

       read(12,'(A)') rdline

       call searching(rdline,"START_EQUATIONS_SET",num_of_EquationsSet)
       call searching(rdline,"START_PROBLEM",num_of_Problem)
       call searching(rdline,"START_BOUNDARY_CONDITIONS",num_of_BoundaryCondition)
       call searching(rdline,"START_BASIS",num_of_Basis)
       call searching(rdline,"PRESSURE_BASIS",num_of_PressureBasis)
       call searching(rdline,"START_MESH",num_of_Mesh)
       call searching(rdline,"START_MATERIAL_FIELD",num_of_MaterialField)
       call searching(rdline,"START_COORDINATE_SYSTEM",num_of_CoordinateSystem)
       call searching(rdline,"START_REGION",num_of_Region)
       call searching(rdline,"START_GEOMETRIC_FIELD",num_of_GeometricField)
       call searching(rdline,"START_DEPENDENT_FIELD",num_of_DependentField)
       call searching(rdline,"START_SOLVER_SETTINGS",num_of_Solver)
       call searching(rdline,"START_CONTROL_LOOP",num_of_ControlLoop)
       call searching(rdline,"START_FIBER_FIELD",num_of_FiberField)
       call searching(rdline,"START_PRESSURE_BASIS",num_of_PressureBasis)
       call searching(rdline,"START_FIELD",num_of_Field)

  enddo
  num_of_Decomposition = num_of_Mesh               !! This is wrong , though works for this study. Gonna fix it later.                
  num_of_GeneratedMesh = num_of_Mesh    	   !! As per my understanding these two always seem  to be equal.


  !!!!!!!!!!!!!!!!!		ALLOCATE DATA STRUCTURES BASED ON THE PARAMETERS INITIALIZED ABOVE.  !!!!!!!!!!!!!!!!!

  allocate(all_Basis%Basis(num_of_Basis))
  allocate(all_PressureBasis%PressureBasis(num_of_PressureBasis))
  allocate(all_BoundaryConditions%BoundaryConditions(num_of_BoundaryCondition))
  allocate(all_CoordinateSystem%CoordinateSystem(num_of_CoordinateSystem))
  allocate(all_WorldCoordinateSystem%WorldCoordinateSystem(num_of_CoordinateSystem))
  allocate(all_Mesh%Mesh(num_of_Mesh))
  allocate(all_Decomposition%Decomposition(num_of_decomposition))
  allocate(all_Equations%Equations(num_of_EquationsSet))
  allocate(all_EquationsSet%EquationsSet(num_of_EquationsSet))
  allocate(all_GeometricField%GeometricField(num_of_GeometricField))
  allocate(all_FibreField%FibreField(num_of_FiberField))
  allocate(all_MaterialField%MaterialField(num_of_MaterialField))
  allocate(all_EquationsSetField%EquationsSetField(num_of_EquationsSet))
  allocate(all_Fields%Fields(num_of_Field))
  allocate(all_Problem%Problem(num_of_Problem))
  allocate(all_Region%Region(num_of_Region))
  allocate(all_WorldRegion%WorldRegion(num_of_WorldRegion))  
  allocate(all_Solver%Solver(num_of_Solver))
  allocate(all_LinearSolver%LinearSolver(num_of_Solver))
  allocate(all_SolverEquations%SolverEquations(num_of_Solver))
  allocate(all_ControlLoop%ControlLoop(1))     !! Hard coded , gonna fix it later. 
  allocate(all_GeneratedMesh%GeneratedMesh(num_of_GeneratedMesh))
  allocate(all_DependentField%DependentField(num_of_DependentField))




