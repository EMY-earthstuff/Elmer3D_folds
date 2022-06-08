!sliding function adapted from Elmer/Ice tutorial material available at http://elmerice.elmerfem.org/courses-tutorials by EMY
!includes some pretty crude surge cycle forcing for 3x surges
FUNCTION beta_val(  Model, Node, InputArray) RESULT(beta)
  ! provides you with most Elmer functionality
  USE DefUtils
  ! saves you from stupid errors
  IMPLICIT NONE
  ! the external variables
  !----------------------------------------------------------------------------
  TYPE(Model_t) :: Model     ! the access point to everything about the model
  INTEGER :: Node            ! the current Node number
  REAL(KIND=dp) :: InputArray(3) ! Contains the arguments passed to the function
  REAL(KIND=dp) :: beta      ! the result
  REAL(KIND=dp) :: coordxs, coordys
  !----------------------------------------------------------------------------
  ! internal variables
  !----------------------------------------------------------------------------
  REAL(KIND=dp) :: Bq_background, Bq, Bsmax, m_surge, b_surge, m_q, b_q, &
       inittime, time, elevation, cutoff
  LOGICAL :: FirstTime=.TRUE.
  ! Remember this value
  SAVE FirstTime, inittime


  ! lets hard-code our values (if we have time we can later make them being read from SIF)
  Bq_background = 0.1_dp
  Bq = 0.01_dp
  Bsmax = 0.000001_dp
  m_surge = -49.95_dp
  b_surge = 1148.95_dp
  m_q = -0.495_dp
  b_q = 10.495_dp
  
  ! copy input (should match the arguments!)
  time = InputArray(1)
  Coordxs = InputArray(2)
  Coordys = InputArray(3)
  !WRITE (Message, '(A,E10.2,A,E10.2)')  "time=", time, "Coordxs=", Coordxs, "Coordys=", Coordys
  CALL INFO("beta_val", Message, Level=9)

  ! store the initial time, to be sure to have relative times
  IF (FirstTime) THEN
     inittime = time
     FirstTime = .FALSE.
  END IF


  ! get the time dependent beta values
  IF (Coordxs > 15000.0) THEN
     beta = Bq
  ELSE IF (Coordys < 14000.0) THEN
     beta = Bq
  ELSE 
    IF (time < 20.0) THEN
      beta = Bq_background
    ELSE IF (20.0 <= time .AND. time < 41.0) THEN
      beta = (time*m_q + 20.395) * Bq
    ELSE IF (41.0 <= time .AND. time <= 43.0) THEN
      beta = (time*m_surge + 2147.95) * Bsmax
    ELSE IF (43.0 < time .AND. time < 64.0) THEN
      beta = Bq_background
    ELSE IF (64.0 <= time .AND. time < 85.0) THEN
      beta = (time*m_q + 42.175) * Bq
    ELSE IF (85.0 <= time .AND. time <= 87.0) THEN
      beta = (time*m_surge + 4345.75) * Bsmax
    ELSE IF (87.0 < time .AND. time < 108.0) THEN
      beta = Bq_background
    ELSE IF (108.0 <= time .AND. time < 129.0) THEN
      beta = (time*m_q + 63.955) * Bq
    ELSE IF (129.0 <= time .AND. time <= 131.0) THEN
      beta = (time*m_surge + 6543.55) * Bsmax
    ELSE 
      beta = Bq_background

    END IF
  END IF

  ! set the x and y dependent beta values


  RETURN

END FUNCTION beta_val
