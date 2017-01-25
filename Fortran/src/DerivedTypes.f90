
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! _______________________ HEADER COMMENTS ____________________!!!!!!!!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING HEADER FILE INITILIZES THE POINTER TYPE DERIVED DATA STRCUTRES....!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! .... WITH DERIVED TYPES. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!________________________________________________________!!!!!!!!!!!!!!!!!!!!!!!!
  


  !!!!!!!!!!!!!!!!!!!!!!!!! DEFINING DERIVED TYPES AND ALLOCATING THIER RESPECTIVE DATA STRUCTURES !!!!!!!!!!!!!!!!!!!

  type basis_array

	type(cmfe_BasisType) 	         , allocatable :: Basis(:),PressureBasis(:)

  end type basis_array

  type boundary_conditions_array

	TYPE(cmfe_BoundaryConditionsType), allocatable :: BoundaryConditions(:)

  end type  boundary_conditions_array

  type coordinate_system_array

	type(cmfe_CoordinateSystemType)  , allocatable :: CoordinateSystem(:),  & 
 							  WorldCoordinateSystem(:)
  end type coordinate_system_array

  type mesh_array

	type(cmfe_MeshType)              , allocatable :: Mesh(:)

  end type  mesh_array

  type decomposition_array

        type(cmfe_DecompositionType)     , allocatable :: Decomposition(:)

  end type  decomposition_array

  type equations_array

        type(cmfe_EquationsType)         , allocatable :: Equations(:)

  end type equations_array

  type equations_set_array

        type(cmfe_EquationsSetType)      , allocatable :: EquationsSet(:)

  end type equations_set_array

  type field_type_array

        type(cmfe_FieldType)             , allocatable ::  GeometricField(:),FibreField(:), & 
							   MaterialField(:), EquationsSetField(:)
  end type  field_type_array

  type fields_type_array

        type(cmfe_FieldsType)            , allocatable ::  Fields(:)

  end type  fields_type_array

  type problem_type

	type(cmfe_ProblemType)           , allocatable :: Problem(:)

  end type  problem_type

  type region_type

	TYPE(cmfe_RegionType)            , allocatable :: Region(:),WorldRegion(:)

  end type region_type

  type solver_type

	TYPE(cmfe_SolverType)            , allocatable :: Solver(:),LinearSolver(:)
	
  end type solver_type


  type solvers_equations_type

        TYPE(cmfe_SolverEquationsType)   , allocatable :: SolverEquations(:)
	
  end type solvers_equations_type

  type control_loop_type

	TYPE(cmfe_ControlLoopType)       , allocatable :: ControlLoop(:)

  end type control_loop_type

  type generate_mesh_type

	TYPE(cmfe_GeneratedMeshType)     , allocatable :: GeneratedMesh(:)

  end type generate_mesh_type

  type dependent_field_type

 	TYPE(cmfe_FieldType)             , allocatable :: DependentField(:)

  end type dependent_field_type

 !!!!!! INITIALIZING POINTER THAT POINT AT THE DATA STRUCTURES OF DERIVED TYPES !!!!!!!!!!!!!!!!

 TYPE(basis_array)          	      :: all_Basis, all_PressureBasis
 TYPE(boundary_conditions_array)      :: all_BoundaryConditions
 TYPE(coordinate_system_array)        :: all_CoordinateSystem
 TYPE(coordinate_system_array)        :: all_WorldCoordinateSystem
 TYPE(mesh_array)                     :: all_Mesh
 TYPE(decomposition_array) 	      :: all_Decomposition
 TYPE(equations_array) 	              :: all_Equations
 TYPE(equations_set_array)            :: all_EquationsSet
 TYPE(field_type_array)     	      :: all_GeometricField
 TYPE(field_type_array)     	      :: all_FibreField
 TYPE(field_type_array)     	      :: all_MaterialField
 TYPE(field_type_array)     	      :: all_EquationsSetField
 TYPE(fields_type_array)	      :: all_Fields
 TYPE(problem_type)		      :: all_problem
 TYPE(region_type)		      :: all_Region
 TYPE(region_type)		      :: all_WorldRegion
 TYPE(solver_type)		      :: all_Solver
 TYPE(solver_type)		      :: all_LinearSolver
 TYPE(solvers_equations_type)	      :: all_SolverEquations
 TYPE(control_loop_type)	      :: all_ControlLoop
 TYPE(generate_mesh_type)             :: all_GeneratedMesh
 TYPE(dependent_field_type) 	      :: all_DependentField
 TYPE(cmfe_FieldType) 		      :: SourceField
  
