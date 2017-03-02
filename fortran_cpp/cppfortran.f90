! fortran file

PROGRAM TEST
!USE Function

IMPLICIT NONE

!EXTERNAL Function

INTERFACE
  SUBROUTINE CPPFunction()
  END SUBROUTINE CPPFunction
END INTERFACE


  PRINT *, "Hello in fortran file."
  CALL CPPFunction()

  PRINT*, "End fortran routine."


END PROGRAM TEST
