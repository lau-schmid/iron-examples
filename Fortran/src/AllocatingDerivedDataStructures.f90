
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


  open(12,file=fileplace,status="old")
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
enddo

  num_of_WorldCoordinateSystem = num_of_CoordinateSystem
  num_of_Decomposition = num_of_Mesh
  num_of_GeneratedMesh =  num_of_Mesh
  


!!!!!!!!!!!!!!!!!!!		ALLOCATE DATA STRUCTURES 		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

  allocate(all_EquationsSetField%EquationsSetField(1))
  allocate(all_Fields%Fields(1))
  allocate(all_Problem%Problem(1))
  allocate(all_Region%Region(1))
  allocate(all_WorldRegion%WorldRegion(1))
  allocate(all_Solver%Solver(1))
  allocate(all_LinearSolver%LinearSolver(1))
  allocate(all_SolverEquations%SolverEquations(1))
  allocate(all_ControlLoop%ControlLoop(1))
  allocate(all_GeneratedMesh%GeneratedMesh(1))
  allocate(all_DependentField%DependentField(num_of_DependentField))




