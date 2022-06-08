!MB function adapted from Elmer/Ice tutorial material available at http://elmerice.elmerfem.org/courses-tutorials by EMY

FUNCTION MB(  Model, Node, InputArray) RESULT(b_l)
  ! provides you with most Elmer functionality
  USE DefUtils
  ! saves you from stupid errors
  IMPLICIT NONE
  ! the external variables
  !----------------------------------------------------------------------------
  TYPE(Model_t) :: Model     ! the access point to everything about the model
  INTEGER :: Node            ! the current Node number
  REAL(KIND=dp) :: InputArray(3) ! Contains the arguments passed to the function
  REAL(KIND=dp) :: b_l      ! the result
  REAL(KIND=dp) :: coordxs, coordys
  !----------------------------------------------------------------------------
  ! internal variables
  !----------------------------------------------------------------------------
  REAL(KIND=dp) :: res_factor, b, m, time, inittime
  LOGICAL :: FirstTime=.TRUE.
 !  Remember this value
  SAVE FirstTime, inittime

  ! lets hard-code our values (if we have time we can later make them being read from SIF)
! 106SS = 1.4 & 105SS = 1.0   &   506 = 1.1
  res_factor = 1.8_dp
  m = 0.00015_dp
  b = 1.8_dp
  
  ! copy input (should match the arguments!)
  
  !time = InputArray(1)
  Coordxs = InputArray(2)
  Coordys = InputArray(3)
  !WRITE (Message, '(A,E10.2,A,E10.2)')  "time=", time, "Coordxs=", Coordxs, "Coordys=", Coordys
  CALL INFO("MB", Message, Level=9)

  ! store the initial time, to be sure to have relative times
  IF (FirstTime) THEN
    inittime = time
     FirstTime = .FALSE.
  END IF


  ! get the time dependent beta values
  IF (Coordxs > 15000.0) THEN
     b_l = (-Coordxs*m) + b
  ELSE IF (Coordys < 14000.0) THEN
     b_l = (-Coordxs*m) + b
  ELSE 
    IF (((-Coordxs*m) + b) > 0) THEN
      b_l = ((-Coordxs*m) + b)*res_factor
    ELSE
      b_l = (-Coordxs*m) + b
    END IF
  END IF

  ! set the x and y dependent beta values


  RETURN

END FUNCTION MB
