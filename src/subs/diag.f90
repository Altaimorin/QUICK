SUBROUTINE DIAG(NDIM,LWMAX,A,EVAL1,EVEC1,IERROR)

  IMPLICIT doUBLE PRECISION (A-H,O-Z)
  DIMENSION A(NDIM,NDIM),EVAL1(NDIM)
  DIMENSION EVEC1(NDIM,NDIM)

  integer :: LDA, LWORK, LIWORK, IWORK(LWMAX)
  DIMENSION WORK(LWMAX)
  LDA = NDIM
  EVEC1 = A


  if(NDIM == 1)then
     EVAL1(1) = A(1,1)
     EVEC1(1,1) = 1.0D0
     RETURN
  endif

  ! Query the optimal workspace
  LWORK = -1
  LIWORK = -1
  !print *, "start first call"
  call dsyevd('Vectors', 'Upper', NDIM, EVEC1, LDA, EVAL1, WORK, LWORK,IWORK, LIWORK, IERROR)
  LWORK = min(LWMAX, int(WORK(1)))
  LIWORK = min(LWMAX, IWORK(1))

  ! Solve eigenproblem
  !print *, "start second call"
  !print *, "LWORK = ", LWORK
  !print *, "WORK(1) = ", WORK(1)
  !print *, "LIWORK = ", LIWORK
  !print *, "IWORK(1) = ", IWORK(1)

  call dsyevd('Vectors', 'Upper', NDIM, EVEC1, LDA, EVAL1, WORK, LWORK, IWORK, LIWORK, IERROR)

  if (INFO.gt.0) then
        write(*,*) "The algorithm failed to compute eigenvalues."
        stop
    end if

    RETURN
end SUBROUTINE DIAG
