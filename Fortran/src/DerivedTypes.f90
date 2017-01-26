
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! _______________________ HEADER COMMENTS ____________________!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING HEADER FILE INITILIZES THE POINTER TYPE DERIVED DATA STRCUTRES....!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! .... WITH DERIVED TYPES. !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!________________________________________________________!!!!!!!!!!!!!!!!!!!!!!!!



  !!!!!!!!!!!!!!!!!!!!!!!!! DEFINING DERIVED TYPES AND ALLOCATING THIER RESPECTIVE DATA STRUCTURES !!!!!!!!!!!!!!!!!!!

  TYPE basis_array

    TYPE(cmfe_BasisType) 	         , allocatable :: Basis(:),PressureBasis(:)

  end TYPE basis_array

  TYPE boundary_conditions_array

    TYPE(cmfe_BoundaryConditionsType)    , allocatable :: BoundaryConditions(:)

  end TYPE  boundary_conditions_array

  TYPE coordinate_system_array

    TYPE(cmfe_CoordinateSystemType)      , allocatable :: CoordinateSystem(:), WorldCoordinateSystem(:)

  end TYPE coordinate_system_array

  TYPE mesh_array

    TYPE(cmfe_MeshType)                 , allocatable :: Mesh(:)

  end TYPE  mesh_array

  TYPE decomposition_array

    TYPE(cmfe_DecompositionType)       , allocatable :: Decomposition(:)

  end TYPE  decomposition_array

  TYPE equations_array

    TYPE(cmfe_EquationsType)           , allocatable :: Equations(:)

  end TYPE equations_array

  TYPE equations_set_array

    TYPE(cmfe_EquationsSetType)        , allocatable :: EquationsSet(:)

  end TYPE equations_set_array

  TYPE field_type_array

    TYPE(cmfe_FieldType)            , allocatable :: GeometricField(:),FibreField(:),MaterialField(:),EquationsSetField(:)

  end TYPE  field_type_array

  TYPE fields_type_array

    TYPE(cmfe_FieldsType)           , allocatable ::  Fields(:)

  end TYPE  fields_type_array

  TYPE problem_type

    TYPE(cmfe_ProblemType)           , allocatable :: Problem(:)

  end TYPE  problem_type

  TYPE region_type

    TYPE(cmfe_RegionType)            , allocatable :: Region(:),WorldRegion(:)

  end TYPE region_type

  TYPE solver_type

    TYPE(cmfe_SolverType)            , allocatable :: Solver(:),LinearSolver(:)

  end TYPE solver_type


  TYPE solvers_equations_type

    TYPE(cmfe_SolverEquationsType)   , allocatable :: SolverEquations(:)

  end TYPE solvers_equations_type

  TYPE control_loop_type

    TYPE(cmfe_ControlLoopType)       , allocatable :: ControlLoop(:)

  end TYPE control_loop_type

  TYPE generate_mesh_type

    TYPE(cmfe_GeneratedMeshType)     , allocatable :: GeneratedMesh(:)

  end TYPE generate_mesh_type

  TYPE dependent_field_type

    TYPE(cmfe_FieldType)             , allocatable :: DependentField(:)

  end TYPE dependent_field_type

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
 TYPE(cmfe_FieldTYPE) 		      :: SourceField

