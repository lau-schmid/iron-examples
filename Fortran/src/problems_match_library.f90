
FUNCTION match_problem(type_string) RESULT (type_number)

    IMPLICIT NONE
    character(*), intent(in) :: type_string
    integer                  :: type_number

    ! class
    IF (type_string == "PROBLEM_ELASTICITY_CLASS") THEN
      type_number = CMFE_PROBLEM_ELASTICITY_CLASS
    ELSEIF (type_string == "PROBLEM_FLUID_MECHANICS_CLASS") THEN
      type_number = CMFE_PROBLEM_FLUID_MECHANICS_CLASS
    ELSEIF (type_string == "PROBLEM_ELECTROMAGNETICS_CLASS") THEN
      type_number = CMFE_PROBLEM_ELECTROMAGNETICS_CLASS
    ELSEIF (type_string == "CLASSICAL_FIELD_CLASS") THEN
      type_number = CMFE_PROBLEM_CLASSICAL_FIELD_CLASS
    ELSEIF (type_string == "PROBLEM_BIOELECTRICS_CLASS") THEN
      type_number = CMFE_PROBLEM_BIOELECTRICS_CLASS
    ELSEIF (type_string == "PROBLEM_MODAL_CLASS") THEN
      type_number = CMFE_PROBLEM_MODAL_CLASS
    ELSEIF (type_string == "PROBLEM_FITTING_CLASS") THEN
      type_number = CMFE_PROBLEM_FITTING_CLASS
    ELSEIF (type_string == "PROBLEM_OPTIMISATION_CLASS") THEN
      type_number = CMFE_PROBLEM_OPTIMISATION_CLASS
    ELSEIF (type_string == "PROBLEM_MULTI_PHYSICS_CLASS") THEN
      type_number = CMFE_PROBLEM_MULTI_PHYSICS_CLASS
    ELSEIF (type_string == "PROBLEM_MULTI_PHYSICS_CLASS") THEN
      type_number = CMFE_PROBLEM_MULTI_PHYSICS_CLASS
    ELSEIF (type_string == "PROBLEM_CLASSICAL_FIELD_CLASS") THEN
      type_number = CMFE_PROBLEM_CLASSICAL_FIELD_CLASS
    ELSEIF (type_string == "PROBLEM_NO_CLASS") THEN
      type_number = CMFE_PROBLEM_NO_CLASS
    ELSEIF (type_string == "PROBLEM_FINITE_ELASTICITY_TYPE") THEN
      type_number = CMFE_PROBLEM_FINITE_ELASTICITY_TYPE
    ELSEIF (type_string == "PROBLEM_LAPLACE_EQUATION_TYPE") THEN
      type_number = CMFE_PROBLEM_LAPLACE_EQUATION_TYPE
    ELSEIF (type_string == "PROBLEM_STANDARD_LAPLACE_SUBTYPE") THEN
      type_number =  CMFE_PROBLEM_STANDARD_LAPLACE_SUBTYPE
    ! TODO add more types.. TODO
    ! subtype
    ELSEIF (type_string == "U_VARIABLE") THEN

      type_number = CMFE_FIELD_U_VARIABLE_TYPE
    ELSEIF (type_string == "SEPARATED") THEN
      type_number = CMFE_FIELD_SEPARATED_COMPONENT_DOF_ORDER
    ELSEIF (type_string == "PROBLEM_NO_SUBTYPE") THEN
      type_number = CMFE_PROBLEM_NO_SUBTYPE

    ! TODO add more subtypes.. TODO
    ELSE
      CALL handle_error("Invalid string "// type_string)
    ENDIF

END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!.DEFINES MESH TYPE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION match_generated_mesh(type_string) RESULT (type_number)
    IMPLICIT NONE
    character(*), intent(in) :: type_string
    integer                  :: type_number
    ! number of dimensions
    IF (type_string == "GENERATED_MESH_1D") THEN
      type_number = 1
    ELSEIF (type_string == "GENERATED_MESH_2D") THEN
      type_number = 2
    ELSEIF (type_string == "GENERATED_MESH_3D") THEN
      type_number = 3
    ! generated mesh type
    ELSEIF (type_string == "REGULAR_MESH_TYPE") THEN
      type_number = CMFE_GENERATED_MESH_REGULAR_MESH_TYPE
    ELSEIF (type_string == "POLAR_MESH_TYPE") THEN
      type_number = CMFE_GENERATED_MESH_POLAR_MESH_TYPE
    ELSEIF (type_string == "TREE_MESH_TYPE") THEN
      type_number = CMFE_GENERATED_MESH_FRACTAL_TREE_MESH_TYPE
    ELSEIF (type_string == "CYLINDER_MESH_TYPE") THEN
      type_number = CMFE_GENERATED_MESH_CYLINDER_MESH_TYPE
    ELSEIF (type_string == "ELLIPSOID_MESH_TYPE") THEN
      type_number = CMFE_GENERATED_MESH_ELLIPSOID_MESH_TYPE
    ELSE
      CALL handle_error("Invalid string "//type_string)

    ENDIF
END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!.REMOVES INTERPOLATION TYPE FOR STATE VARIABLE!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION match_basis(type_string) RESULT (type_number)

    IMPLICIT NONE
    character(*), intent(in) :: type_string
    integer                  :: type_number
    SELECT CASE (type_string)
    case("LINEAR_LAGRANGE_INTERPOLATION")
    type_number=CMFE_BASIS_LINEAR_LAGRANGE_INTERPOLATION
    case("QUADRATIC_LAGRANGE_INTERPOLATION")
    type_number=CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
    case("CUBIC_LAGRANGE_INTERPOLATION ")
    type_number=CMFE_BASIS_QUADRATIC_LAGRANGE_INTERPOLATION
    case("QUADRATIC1_HERMITE_INTERPOLATION")
    type_number=CMFE_BASIS_QUADRATIC1_HERMITE_INTERPOLATION
    case("CUBIC_HERMITE_INTERPOLATION")
    type_number=CMFE_BASIS_CUBIC_HERMITE_INTERPOLATION
    case("QUADRATIC2_HERMITE_INTERPOLATION")
    type_number=CMFE_BASIS_QUADRATIC2_HERMITE_INTERPOLATION
    case("LINEAR_SIMPLEX_INTERPOLATION")
    type_number=CMFE_BASIS_LINEAR_SIMPLEX_INTERPOLATION
    case("QUADRATIC_SIMPLEX_INTERPOLATION")
    type_number=CMFE_BASIS_QUADRATIC_SIMPLEX_INTERPOLATION
    case("CUBIC_SIMPLEX_INTERPOLATION")
    type_number=CMFE_BASIS_CUBIC_SIMPLEX_INTERPOLATION
    END SELECT

END FUNCTION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!.DEFINES COORDINATE SYSTEM !!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION match_coordinate_system(type_string) RESULT (type_number)

    IMPLICIT NONE
    character(*), intent(in) :: type_string
    integer                  :: type_number

    IF (type_string == "COORDINATE_RECTANGULAR_CARTESIAN_TYPE") THEN
      type_number = CMFE_COORDINATE_RECTANGULAR_CARTESIAN_TYPE
    ELSEIF (type_string == "COORDINATE_CYLINDRICAL_POLAR_TYPE") THEN
      type_number = CMFE_COORDINATE_CYLINDRICAL_POLAR_TYPE
    ELSEIF (type_string == "COORDINATE_SPHERICAL_POLAR_TYPE") THEN
      type_number = CMFE_COORDINATE_SPHERICAL_POLAR_TYPE
    ELSEIF (type_string == "COORDINATE_PROLATE_SPHEROIDAL_TYPE") THEN
      type_number = CMFE_COORDINATE_PROLATE_SPHEROIDAL_TYPE
    ELSEIF (type_string == "COORDINATE_OBLATE_SPHEROIDAL_TYPE") THEN
      type_number = CMFE_COORDINATE_OBLATE_SPHEROIDAL_TYPE
    ELSEIF (type_string == "COORDINATE_SYSTEM_1D") THEN
      type_number = 1
    ELSEIF (type_string == "COORDINATE_SYSTEM_2D") THEN
      type_number = 2
    ELSEIF (type_string == "COORDINATE_SYSTEM_3D") THEN
      type_number = 3
    ELSE

      CALL handle_error("Invalid string "//type_string)

    ENDIF

END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!.REMOVES INTERPOLATION TYPE FOR STATE VARIABLE!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION match_equations_set(type_string) RESULT (type_number)

    IMPLICIT NONE
    character(*), intent(in) :: type_string
    integer                   :: type_number
    ! sparsity type
    IF (type_string == "EQUATIONS_SPARSE_MATRICES") THEN
      type_number = CMFE_EQUATIONS_SPARSE_MATRICES
    ELSEIF (type_string == "EQUATIONS_FULL_MATRICES") THEN
      type_number = CMFE_EQUATIONS_FULL_MATRICES
    ELSEIF (type_string == "SEPARATED") THEN
      type_number = CMFE_FIELD_SEPARATED_COMPONENT_DOF_ORDER
    ! output type
    ELSEIF (type_string == "NO_OUTPUT") THEN
      type_number = CMFE_EQUATIONS_NO_OUTPUT
    ELSEIF (type_string == "EQUATIONS_TIMING_OUTPUT") THEN
      type_number = CMFE_EQUATIONS_TIMING_OUTPUT
    ELSEIF (type_string == "EQUATIONS_MATRIX_OUTPUT") THEN
      type_number = CMFE_EQUATIONS_MATRIX_OUTPUT
    ELSEIF (type_string == "EQUATIONS_ELEMENT_MATRIX_OUTPUT") THEN
      type_number = CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT
    ELSEIF (type_string == "U_VARIABLE") THEN
      type_number = CMFE_FIELD_U_VARIABLE_TYPE
    ! class
    ELSEIF (type_string == "EQUATIONS_SET_ELASTICITY_CLASS") THEN
      type_number = CMFE_EQUATIONS_SET_ELASTICITY_CLASS
    ELSEIF (type_string == "EQUATIONS_SET_CLASSICAL_FIELD_CLASS") THEN
      type_number = CMFE_EQUATIONS_SET_CLASSICAL_FIELD_CLASS
    ! TODO insert more classes.. TODO
    ! type
    ELSEIF (type_string == "EQUATIONS_SET_FINITE_ELASTICITY_TYPE") THEN
      type_number = CMFE_EQUATIONS_SET_FINITE_ELASTICITY_TYPE

    ELSEIF (type_string == "EQUATIONS_SET_LAPLACE_EQUATION_TYPE") THEN
      type_number = CMFE_EQUATIONS_SET_LAPLACE_EQUATION_TYPE
    ! subtype
    ELSEIF (type_string == "EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE") THEN
      type_number = CMFE_EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE
    ELSEIF (type_string == "EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE") THEN
      type_number = CMFE_EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE
    ELSEIF (type_string == "EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE") THEN
      type_number = CMFE_EQUATIONS_SET_STANDARD_LAPLACE_SUBTYPE
    ELSEIF (type_string == "GUESS_TYPE") THEN
      type_number = CMFE_FIELD_VALUES_SET_TYPE
    ELSEIF (type_string == "BOUNDARY_CONDITIONS_SET") THEN
     type_number =  CMFE_FIELD_BOUNDARY_CONDITIONS_SET_TYPE
    ELSEIF (type_string == "UX_VALUE") THEN
     type_number =  1.0
    ELSEIF (type_string == "UY_VALUE") THEN
     type_number =  2.0
    ELSEIF (type_string == "UZ_VALUE") THEN
     type_number =  3.0
    ! TODO insert more subtypes.. TODO
    ELSE
      CALL handle_error("Invalid string "//type_string)

    ENDIF

END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!.READS KEYWORDDS TO DEFINE the OUTPUT!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!OF THE ANALYSIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


FUNCTION output_type(type_string) RESULT (type_number)

    IMPLICIT NONE
    character(*), intent(in) :: type_string
    integer                  :: type_number
    ! sparsity type
    IF (type_string == "EQUATIONS_SPARSE_MATRICES") THEN
      type_number = CMFE_EQUATIONS_SPARSE_MATRICES
    ELSEIF (type_string == "EQUATIONS_FULL_MATRICES") THEN
      type_number = CMFE_EQUATIONS_FULL_MATRICES
    ELSEIF (type_string == "SEPARATED") THEN
      type_number = CMFE_FIELD_SEPARATED_COMPONENT_DOF_ORDER
    ! output type
    ELSEIF (type_string == "EQUATIONS_NO_OUTPUT") THEN
      type_number = CMFE_EQUATIONS_NO_OUTPUT
    ELSEIF (type_string == "EQUATIONS_TIMING_OUTPUT") THEN
      type_number = CMFE_EQUATIONS_TIMING_OUTPUT
    ELSEIF (type_string == "EQUATIONS_MATRIX_OUTPUT") THEN
      type_number = CMFE_EQUATIONS_MATRIX_OUTPUT
    ELSEIF (type_string == "EQUATIONS_ELEMENT_MATRIX_OUTPUT") THEN
      type_number = CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT
    ELSE

      CALL handle_error("Invalid output string  "//type_string)

    ENDIF
END FUNCTION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!.READS KEYWORDDS TO DEFINE the depENDent field!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION match_depENDent_field(type_string) RESULT (type_number)
    character(*), intent(in) :: type_string
    integer                  :: type_number
    IF (type_string == "FIELD_U_VARIABLE_TYPE") THEN
      type_number = CMFE_FIELD_U_VARIABLE_TYPE
    END IF

END FUNCTION match_depENDent_field


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!.READS KEYWORDDS TO DEFINE BCs!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION bc_def(type_string) RESULT (type_number)
    character(*), intent(in) :: type_string
    integer                :: type_number

     IF  (type_string == "RIGHT") THEN
       type_number = CMFE_GENERATED_MESH_REGULAR_RIGHT_SURFACE
     ELSEIF (type_string == "LEFT") THEN
       type_number = CMFE_GENERATED_MESH_REGULAR_LEFT_SURFACE
     ELSEIF (type_string == "FRONT") THEN
       type_number = CMFE_GENERATED_MESH_REGULAR_FRONT_SURFACE
     ELSEIF (type_string == "BOTTOM") THEN
       type_number = CMFE_GENERATED_MESH_REGULAR_BOTTOM_SURFACE
     ELSEIF (type_string == "TOP") THEN
       type_number = CMFE_GENERATED_MESH_REGULAR_TOP_SURFACE
     ELSE

       CALL handle_error("Invalid string "//type_string)

     END IF
END FUNCTION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!.READS KEYWORDDS TO DEFINE SOLVER!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION solver_def(type_string) RESULT (type_number)
    IMPLICIT NONE
    character(*), intent(in) :: type_string
    integer                  :: type_number
    IF (type_string == "SPARSE_MATRIX") THEN
      type_number = CMFE_EQUATIONS_SPARSE_MATRICES
    ELSEIF (type_string == "FULL_MATRICES") THEN
      type_number = CMFE_EQUATIONS_FULL_MATRICES
    ELSEIF (type_string == "SEPARATED") THEN
      type_number = CMFE_FIELD_SEPARATED_COMPONENT_DOF_ORDER
    ELSEIF (type_string == "NO_OUTPUT") THEN
      type_number = CMFE_EQUATIONS_NO_OUTPUT
    ELSEIF (type_string == "EQUATIONS_TIMING_OUTPUT") THEN
      type_number = CMFE_EQUATIONS_TIMING_OUTPUT
    ELSEIF (type_string == "EQUATIONS_MATRIX_OUTPUT") THEN
      type_number = CMFE_EQUATIONS_MATRIX_OUTPUT
    ELSEIF (type_string == "EQUATIONS_ELEMENT_MATRIX_OUTPUT") THEN
      type_number = CMFE_EQUATIONS_ELEMENT_MATRIX_OUTPUT
    ELSEIF (type_string == "SOLVER_LINEAR_DIRECT_SOLVE_TYPE") THEN
      type_number = CMFE_SOLVER_LINEAR_DIRECT_SOLVE_TYPE
    ELSEIF (type_string == "ITERATIVE") THEN
      type_number = CMFE_SOLVER_LINEAR_ITERATIVE_SOLVE_TYPE
    ELSEIF (type_string == "INCREMENTAL") THEN
      type_number = CMFE_PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE
    ELSEIF (type_string == "SOLVER_NEWTON_JACOBIAN_FD_CALCULATED") THEN
      type_number = CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED
    ELSEIF (type_string == "SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED") THEN
      type_number = CMFE_SOLVER_NEWTON_JACOBIAN_EQUATIONS_CALCULATED
    ELSEIF (type_string == "SOLVER_NEWTON_JACOBIAN_FD_CALCULATED") THEN
      type_number = CMFE_SOLVER_NEWTON_JACOBIAN_FD_CALCULATED
    ELSE
      CALL handle_error("Invalid string "//type_string)

    ENDIF

END FUNCTION


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!.READS KEYWORDDS TO DEFINE SOLVER!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION control_loop_def(type_string) RESULT (type_number)
    IMPLICIT NONE
    character(*), intent(in) :: type_string
    integer                  :: type_number
    IF (type_string == "CONTROL_LOAD_INCREMENT_LOOP_TYPE") THEN
      type_number = CMFE_PROBLEM_CONTROL_LOAD_INCREMENT_LOOP_TYPE
    !!! add more options later
    ELSE
      CALL handle_error("Invalid string "//type_string)
    ENDIF

END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION !!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!.READS KEYWORDDS TO DEFINE MATERIALS!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION material_parameters(type_string) RESULT (type_number)

    IMPLICIT NONE
    character(*), intent(in) :: type_string
    integer                 :: type_number
    IF (type_string == "EQUATIONS_SET_TRANSVERSE_ISOTROPIC_ACTIVE_SUBTYPE") THEN
      type_number = 5
    ELSE IF (type_string == "EQUATIONS_SET_MOONEY_RIVLIN_SUBTYPE") THEN
       type_number =2
    ELSEIF (type_string == "EQUATIONS_SET_THREE_DIMENSIONAL_SUBTYPE")  THEN
      type_number = 2
    ELSE
      CALL handle_error("Invalid string "//type_string)

    ENDIF

END FUNCTION material_parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!! THE FOLLOWING FUNCTION  !!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!..used for error handling !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE HANDLE_ERROR(ERROR_STRING)

    CHARACTER(LEN=*), INTENT(IN) :: ERROR_STRING
    WRITE(*,'(">>ERROR: ",A)') ERROR_STRING(1:LEN_TRIM(ERROR_STRING))
    STOP

END SUBROUTINE HANDLE_ERROR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

