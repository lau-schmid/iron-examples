




                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        !!!!!!!!!                  LAPLACE CASE STUDY              !!!!!!!!!!!!!!!!!!!make

                        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! EQUATION_SET block defines the governing equation!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! It creates an object of derived type CMFE_EquationsSetType !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_EQUATIONS_SET

MODEL_ID
SOLID

CLASS
EQUATIONS_SET_CLASSICAL_FIELD_CLASS

TYPE
EQUATIONS_SET_LAPLACE_EQUATION_TYPE


SUBTYPE
EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE

EQUATIONS_OUTPUT_TYPE_SET
EQUATIONS_NO_OUTPUT

END_EQUATIONS_SET


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! PROBLEM block defines the governing equation!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! It creates an object of derived type CMFE_ProblemType !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


START_PROBLEM

PROBLEM_ID
SOLID

CLASS
PROBLEM_CLASSICAL_FIELD_CLASS

TYPE
PROBLEM_LAPLACE_EQUATION_TYPE

SUBTYPE
PROBLEM_STANDARD_LAPLACE_SUBTYPE

END_PROBLEM


!!!!!!!!!!!!! Define Region !!!!!!!!!!!!!!!!!!!!!!!
START_REGION

REGION_ID
SOLID

REGION_LABEL_SET
REGION

END_REGION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! Creates a region to which geometric entities like ...!!!!!!!!!!!!!!!
!!!!!!!!! Cooridnate System, Mesh are attached !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!It creates an object of derived type CMFE_RegionType!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_GEOMETRIC_FIELD

GEOMETRIC_FIELD_ID
SOLID

GEOMETRIC_FIELD_LABEL_SET
GEOMETRY

!!!! rest functionailites will be added later

END_GEOMETRIC_FIELD


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! The following block defines an incremental loop!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!It creates an object of derived type CMFE_ControlLoopType!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_CONTROL_LOOP

CONTROL_LOOP_ID
SOLID

TYPE
PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE

INCREMENTS
1

END_CONTROL_LOOP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! The following block defines a solver !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!It creates an object of derived type CMFE_SolverType !!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_SOLVER_SETTINGS

SOLVER_ID
SOLID

LINEAR_SOLVER_TYPE                     !! type of linear solver , direct or iterative
SOLVER_LINEAR_DIRECT_SOLVE_TYPE

SOLVER_LIBRARY_TYPE
SOLVER_MUMPS_LIBRARY

!!!!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FEATURES/KEYWORDS SHOULD BE INCLUDED FOR THE DIRECT SOLVERS !!!!!!!!!!!!!!!!!!

EQUATION_DIRECT_SOLVER_TYPE
SOLVER_DIRECT_LU

!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FEATURES/KEYWORDS SHOULD BE INCLUDED FOR THE ITERATIVE SOLVERS !!!!!!!!!!!!!!!!!!

EQUATION_ITERATIVE_SOLVER_TYPE
SOLVER_ITERATIVE_GMRES

LINEAR_ITERATIVE_RELATIVE_TOLERANCE    !! User can can use both relative and absolute tolerance or ignore one of them
1E-6

LINEAR_ITERATIVE_ABSOLUTE_TOLERANCE    !! User can can use both relative and absolute tolerance or ignore one of them
1E-12

LINEAR_ITERATIVE_MAXIMUM_ITERATIVE     !! Maximum iterations for the linear iterative solvers
100000

LINEAR_ITERATIVE_PRECONDITIONER_TYPE   !! THis parameter can be ignored by the user. The default type is no preconditioner
SOLVER_ITERATIVE_JACOBI_PRECONDITIONER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

NONLINEAR_SOLVER                       !! for a quasi static problem this parameter is STATIC_NONLINEAR
Dynamic_nonlinear                      !! "_nonlinear" refers to the fact that the equilibrium is achieved for each load/time..
                                       !! .. increment using newton raphson iterations

NEWTON_SOLVER_MAXIMUM_ITERATIONS       !! Maximum number of Newton Raphson iterations
1000

NEWTON_SOLVER_ABSOLUTE_TOLERANCE       !! User can can use both relative and absolute tolerance or ignore one of them
1e-6

NEWTON_SOLVER_RELATIVE_TOLERANCE       !! User can can use both relative and absolute tolerance or ignore one of them
1e-6

JACOBIAN_TYPE                          !! default one is the one caluclated analytically using double derivation of the energy
SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED !! .. function

SOLVER_OUTPUT_TYPE_SET                 !! display solver output on screen m informations such as residual , total number of
SOLVER_PROGRESS_OUTPUT                 !! solver iterations

END_SOLVER_SETTINGS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! The following block defines a coordinate system !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!It creates an object of derived type CMFE_CoordinateSystemType!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_COORDINATE_SYSTEM

COORDINATE_SYSTEM_ID
FLUID

COORDINATE_SYSTEM_TYPE                 !! type cooridnate system
COORDINATE_RECTANGULAR_CARTESIAN_TYPE

COORDINATE_SYSTEM_DIMENSION            !! dimensions of the coordinate system
COORDINATE_SYSTEM_3D

END_COORDINATE_SYSTEM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! The following block defines a Dependent Field !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!It creates an object of derived type CMFE_FieldType!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_DEPENDENT_FIELD

DEPENDENT_FIELD_ID
SOLID


NUM_OF_COMPONENTS
1

INITIAL_VALUE_OF_STATE_VARIABLE
0.5

DEPENDENT_FIELD_LABEL_SET
DEPENDENT

END_DEPENDENT_FIELD





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! The following block defines a generated mesh !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!It creates an object of derived type CMFE_GeneratedMeshType!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_GENERATED_MESH

GENERATED_MESH_ID
SOLID


GENERATED_MESH_TYPE
GENERATED_REGULAR_MESH_TYPE

GENERATED_MESH_GEOMETRIC_EXTENTS
1,	1, 1

GENERATED_MESH_ORIGIN_SET
0,0,0

GENERATED_MESH_NUMBER_OF_ELEMENTS
2,  2, 2


END_GENERATED_MESH


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! The following block defines a  mesh !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!It creates an object of derived type CMFE_MeshType !!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_MESH

MESH_ID
SOLID

END_MESH


!!!!!!!!!!!!!!!!!!!! DEFINE BASIS !!!!!!!!!!!!!!!!!!!!!!!!!!!!
START_BASIS

BASIS_ID
SOLID

BASIS_INTERPOLATION_TYPE
LINEAR_LAGRANGE_INTERPOLATION

NumberOfGaussXi
2,2,2

END_BASIS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! The following block defines boundary condition!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!It creates an object of derived type CMFE_BoundaryCOnditionType!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


START_BOUNDARY_CONDITIONS

BOUNDARY_CONDITIONS_ID
SOLID

DIRICHILET  !!!! [ TYPE , LOCATION , COMPONENTS , VALUE ]
SURFACE,	MESH_REGULAR_LEFT_SURFACE,	100,		1.0 ,  BOUNDARY_CONDITION_FIXED
SURFACE,	MESH_REGULAR_RIGHT_SURFACE,	100,		0.0 ,  BOUNDARY_CONDITION_FIXED



END_BOUNDARY_CONDITIONS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! FOLLOWING BLOCK DEFINES EXPORT PREFIX !!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_OUTPUT

OUTPUT_ID
SOLID

NODE_EXPORT
NODE_DATA_LAPLACE

ELEMENT_EXPORT
ELEMENT_DATA_LAPLACE


END_OUTPUT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! The following block defines boundary condition!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!It creates an object of derived type CMFE_DecompositionType!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_DECOMPOSITION

DECOMPOSITION_ID
SOLID

CALCULATE_ELEMENT_FACES
FALSE

END_DECOMPOSITION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!! FIELDS BLOCK initilize object of a data type CMFE_FIELDS !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_FIELDS

FIELDS_ID
SOLID

END_FIELDS



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!! DEFINING Functions that are prescribed in BC block !!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

START_FUNCTION

FUNCTION_ID
QUADRATIC_1

FUNCTION_CONSTANTS                               !!! A,B,C,E,F,G,H,I,J in equation ....
0,-4,0,0,0,0,0,8,0,0
                                                 !!! f(x,y,z) = A*x^2 + B*y^2 + C*z^2 + D*x*y + E*x*z + F*y*z ...
                                                 !!! ... + G*x + H*y + I*z + J
                                                 !!! The given constants correpond to a parabolic profile of ...
                                                 !! .. shape V = 4*V_max * (y/h)(1-(y/h)) , where V_max = 4 units and h is 2 units

END_FUNCTION


STOP_PARSING

