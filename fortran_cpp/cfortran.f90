! fortran file

PROGRAM TEST
!USE Function

IMPLICIT NONE

!EXTERNAL Function

INTERFACE
  SUBROUTINE CFunction()
  END SUBROUTINE CFunction
END INTERFACE


  PRINT *, "Hello in fortran file."
  CALL CFunction()

  PRINT*, "End fortran routine."


END PROGRAM TEST
