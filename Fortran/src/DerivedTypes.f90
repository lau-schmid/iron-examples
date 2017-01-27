
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! _______________________ HEADER COMMENTS ____________________!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING HEADER FILE INITILIZES THE POINTER TYPE DERIVED DATA STRCUTRES....!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! .... WITH DERIVED TYPES. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!________________________________________________________!!!!!!!!!!!!!!!!!!!!!!!!



  !!!!!!!!!!!!!!!!!!!!!!!!! DEFINING DERIVED TYPES AND ALLOCATING THIER RESPECTIVE DATA STRUCTURES !!!!!!!!!!!!!!!!!!!

  TYPE basis_array

    TYPE(cmfe_BasisType) 	         , allocatable :: Basis(:),PressureBasis(:)

  END TYPE basis_array

  TYPE boundary_conditions_array

    TYPE(cmfe_BoundaryConditionsType)    , allocatable :: BoundaryConditions(:)

  END TYPE  boundary_conditions_array

  TYPE coordinate_system_array

    TYPE(cmfe_CoordinateSystemType)      , allocatable :: CoordinateSystem(:), WorldCoordinateSystem(:)

  END TYPE coordinate_system_array

  TYPE mesh_array

    TYPE(cmfe_MeshType)                 , allocatable :: Mesh(:)

  END TYPE  mesh_array

  TYPE decomposition_array

    TYPE(cmfe_DecompositionType)       , allocatable :: Decomposition(:)

  END TYPE  decomposition_array

  TYPE equations_array

    TYPE(cmfe_EquationsType)           , allocatable :: Equations(:)

  END TYPE equations_array

  TYPE equations_set_array

    TYPE(cmfe_EquationsSetType)        , allocatable :: EquationsSet(:)

  END TYPE equations_set_array

  TYPE field_type_array

    TYPE(cmfe_FieldType)            , allocatable :: GeometricField(:),FibreField(:),MaterialField(:),EquationsSetField(:)

  END TYPE  field_type_array

  TYPE fields_type_array

    TYPE(cmfe_FieldsType)           , allocatable ::  Fields(:)

  END TYPE  fields_type_array

  TYPE problem_type

    TYPE(cmfe_ProblemType)           , allocatable :: Problem(:)

  END TYPE  problem_type

  TYPE region_type

    TYPE(cmfe_RegionType)            , allocatable :: Region(:),WorldRegion(:)

  END TYPE region_type

  TYPE solver_type

    TYPE(cmfe_SolverType)            , allocatable :: Solver(:),LinearSolver(:)

  END TYPE solver_type


  TYPE solvers_equations_type

    TYPE(cmfe_SolverEquationsType)   , allocatable :: SolverEquations(:)

  END TYPE solvers_equations_type

  TYPE control_loop_type

    TYPE(cmfe_ControlLoopType)       , allocatable :: ControlLoop(:)

  END TYPE control_loop_type

  TYPE generate_mesh_type

    TYPE(cmfe_GeneratedMeshType)     , allocatable :: GeneratedMesh(:)

  END TYPE generate_mesh_type

  TYPE depENDent_field_type

    TYPE(cmfe_FieldType)             , allocatable :: DepENDentField(:)

  END TYPE depENDent_field_type

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
 TYPE(depENDent_field_type) 	      :: all_DepENDentField
 TYPE(cmfe_FieldTYPE) 		      :: SourceField

