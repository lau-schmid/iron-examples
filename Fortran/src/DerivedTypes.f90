
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! _______________________ HEADER COMMENTS ____________________!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING HEADER FILE INITILIZES THE POINTER TYPE DERIVED DATA STRCUTRES....!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! .... WITH DERIVED TYPES. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!________________________________________________________!!!!!!!!!!!!!!!!!!!!!!!!



  !!!!!!!!!!!!!!!!!!!!!!!!! DEFINING DERIVED TYPES AND ALLOCATING THIER RESPECTIVE DATA STRUCTURES !!!!!!!!!!!!!!!!!!!

  TYPE basis_array

    TYPE(cmfe_BasisType) 	         , ALLOCATABLE :: Basis(:),PressureBasis(:)

  END TYPE basis_array

  TYPE boundary_conditions_array

    TYPE(cmfe_BoundaryConditionsType)    , ALLOCATABLE :: BoundaryConditions(:)

  END TYPE  boundary_conditions_array

  TYPE coordinate_system_array

    TYPE(cmfe_CoordinateSystemType)      , ALLOCATABLE :: CoordinateSystem(:), WorldCoordinateSystem(:)

  END TYPE coordinate_system_array

  TYPE mesh_array

    TYPE(cmfe_MeshType)                  , ALLOCATABLE :: Mesh(:)

  END TYPE  mesh_array

  TYPE decomposition_array

    TYPE(cmfe_DecompositionType)         , ALLOCATABLE :: Decomposition(:)

  END TYPE  decomposition_array

  TYPE equations_array

    TYPE(cmfe_EquationsType)             , ALLOCATABLE :: Equations(:)

  END TYPE equations_array

  TYPE equations_set_array

    TYPE(cmfe_EquationsSetType)          , ALLOCATABLE :: EquationsSet(:)

  END TYPE equations_set_array

  TYPE field_type_array

    TYPE(cmfe_FieldType)                 , ALLOCATABLE :: GeometricField(:),FibreField(:),MaterialField(:),EquationsSetField(:)

  END TYPE  field_type_array

  TYPE fields_type_array

    TYPE(cmfe_FieldsType)                , ALLOCATABLE ::  Fields(:)

  END TYPE  fields_type_array

  TYPE problem_type

    TYPE(cmfe_ProblemType)               , ALLOCATABLE :: Problem(:)

  END TYPE  problem_type

  TYPE region_type

    TYPE(cmfe_RegionType)                , ALLOCATABLE :: Region(:),WorldRegion(:)

  END TYPE region_type

  TYPE solver_type

    TYPE(cmfe_SolverType)                , ALLOCATABLE :: Solver(:),LinearSolver(:)

  END TYPE solver_type


  TYPE solvers_equations_type

    TYPE(cmfe_SolverEquationsType)       , ALLOCATABLE :: SolverEquations(:)

  END TYPE solvers_equations_type

  TYPE control_loop_type

    TYPE(cmfe_ControlLoopType)           , ALLOCATABLE :: ControlLoop(:)

  END TYPE control_loop_type

  TYPE generate_mesh_type

    TYPE(cmfe_GeneratedMeshType)         , ALLOCATABLE :: GeneratedMesh(:)

  END TYPE generate_mesh_type

  TYPE Dependent_field_type

    TYPE(cmfe_FieldType)                 , ALLOCATABLE :: DependentField(:),SourceField(:)

  END TYPE Dependent_field_type

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
 TYPE(Dependent_field_type) 	      :: all_DependentField
 TYPE(Dependent_field_type) 	      :: all_SourceField

