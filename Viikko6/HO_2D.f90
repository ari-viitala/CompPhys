PROGRAM Ho_2d
!
! Solves interacting bosons (or two opposite spin fermions)
! in 2D harmonic oscillator using VMC

  INTEGER, PARAMETER :: in=5,out=6
  INTEGER, PARAMETER :: dp = KIND(1.d0)
  INTEGER :: Npart=2, d_space=2
  INTEGER :: Nw, iw, ii, Nsteps(2)
  REAL(dp), ALLOCATABLE :: r(:,:,:)
  REAL(dp) :: a, step_size, Coulomb=1.d0

  CALL RANDOM_SEED()

  WRITE(out,*) ' Give Coulomb strength (try 1)'
  READ(in,*) Coulomb

  WRITE(out,*) ' Give variational parameter (around 0.5)'
  READ(in,*) a

  step_size=1.d0

  WRITE(out,*) ' How many walkers? (around 10)'
  READ(in,*) Nw
  ALLOCATE(r(Nw,Npart,d_space))

  CALL RANDOM_NUMBER(r)
  r=r-.5d0

  DO ii=1, 100
    CALL metropolis(5,.TRUE.) ! Thermalization
  END DO

  WRITE(out,*) ' Give number of metropolis steps (around 1000)'
  READ(in,*) Nsteps(1)

  if (Nsteps(1)<0) call para_scan()

  CALL metropolis(Nsteps(1))

  WRITE(out,*) ' Give number of optimization steps (around 1000)'
  READ(in,*) Nsteps(2)

  CALL opt(Nsteps(2))

  WRITE(out,*) 'Doing again Metropolis:'
  CALL metropolis(Nsteps(1))

CONTAINS
  FUNCTION el(r)
! Function for local energy
    REAL(dp) :: r(:,:), el, t, h=1.d-4, d(d_space)
    INTEGER :: ipart, jpart, id
    el=.5d0*SUM(r**2)
    DO ipart=1, Npart
      DO jpart=ipart+1, Npart
        el=el+Coulomb/SQRT(SUM( (r(ipart,:)-r(jpart,:))**2 ))
      END DO
    END DO
    t=0.d0
    d=0.d0
    DO ipart=1, Npart
      DO id=1, d_space
        d(id)=h
! 3-point formula for laplacian:
!        t=t+1.d0-.5d0*(psi_ratio(r,d,ipart)+psi_ratio(r,-d,ipart))
! 5-point formula for laplacian:
        t=t+(15.d0-8.d0*(psi_ratio(r,d,ipart)+psi_ratio(r,-d,ipart))&
             &+.5d0*(psi_ratio(r,2.d0*d,ipart)+psi_ratio(r,-2.d0*d,ipart)))&
             &/12.d0
        d(id)=0.d0
      END DO
    END DO
    t=t/h**2
    el=el+t
  END FUNCTION el

  FUNCTION psi_ratio(r,d,ipart)
! Wave function ratio, particle ipart is moved by d
    REAL(dp) :: r(:,:), d(d_space), psi_ratio
    INTEGER :: ipart, jpart
    psi_ratio=EXP(-.5d0*(SUM( (r(ipart,:)+d)**2)-SUM(r(ipart,:)**2)))
    DO jpart=1, Npart
      IF (jpart==ipart) CYCLE
      psi_ratio=psi_ratio&
          &*jas(SQRT(SUM((r(ipart,:)+d-r(jpart,:))**2)),a)&
          &/jas(SQRT(SUM((r(ipart,:)-r(jpart,:))**2)),a)
    END DO
  END FUNCTION psi_ratio

  FUNCTION psi(r,aa)
! Actual wave function
    REAL(dp) :: r(:,:), psi, aa
    INTEGER :: ipart, jpart
    psi=EXP(-.5d0*(SUM( r(:,:)**2)))
    DO ipart=1, Npart
      DO jpart=ipart+1, Npart
        psi=psi*jas(SQRT(SUM((r(ipart,:)-r(jpart,:))**2)),aa)
      END DO
    END DO
  END FUNCTION psi

  FUNCTION jas(r,aa)
    REAL(dp) :: jas, r, aa
    jas=1.d0+aa*r
  END FUNCTION jas

  SUBROUTINE metropolis(steps,in_quiet)
    ! This does most of work. Samples psi**2 using Metropolis.
    INTEGER :: i, steps, iw, iacc, irej, ipart
    REAL(dp) :: delta(d_space), xt, wft, r_acc
    REAL(dp) :: e(Nw,2), e_l
    LOGICAL, OPTIONAL :: in_quiet
    LOGICAL :: quiet
    quiet=.FALSE.
    IF (PRESENT(in_quiet)) quiet=in_quiet
    iacc=0
    irej=0
    e=0.d0
    DO i=1, steps ! Metropolis loop
       DO iw=1, Nw ! Loop over walkers
          DO ipart=1, Npart ! Loop over particles
             CALL RANDOM_NUMBER(delta)
             delta=step_size*(delta-0.5d0)*2.d0
             wft=psi_ratio(r(iw,:,:),delta,ipart)
             CALL RANDOM_NUMBER(r_acc)
             IF ( wft**2 > r_acc) THEN
                ! Move accepted
                r(iw,ipart,:)=r(iw,ipart,:)+delta
                iacc=iacc+1
             ELSE
                irej=irej+1
             END IF
             IF (.NOT. quiet) THEN
                ! Let's measure energy:
                e_l=el(r(iw,:,:))
                e(iw,:)=e(iw,:)+((/e_l, e_l**2/)-e(iw,:))/REAL(i)
             END IF
          END DO
       END DO
    END DO
    step_size=step_size-(0.5d0-REAL(iacc)/REAL(iacc+irej))*step_size
    IF (.NOT. quiet) THEN
       WRITE(out,*) 'acc ', REAL(iacc)/REAL(iacc+irej), ' step to ', step_size
       WRITE(out,*) 'Energy is ', SUM(e(:,1))/REAL(nw)

       WRITE(out,*) ' +/-',&
            & SUM(ABS(SUM(e(:,1))/REAL(nw)-e(:,1)))&
            &/REAL(nw*SQRT(MAX(nw-1.d0,1.d0)))

       WRITE(out,*) 'variance  ', (ABS(SUM(e(:,2))/REAL(nw)&
            &-(SUM(e(:,1))/REAL(nw))**2))

       WRITE(out,*) 'aEdEv ', a, SUM(e(:,1))/REAL(nw),&
            & SUM(ABS(SUM(e(:,1))/REAL(nw)-e(:,1)))&
            &/REAL(nw*SQRT(MAX(nw-1.d0,1.d0))), (ABS(SUM(e(:,2))/REAL(nw)&
            &-(SUM(e(:,1))/REAL(nw))**2))

    END IF
  END SUBROUTINE metropolis

  SUBROUTINE opt(steps)
! This does the optimization, see lecture notes.
    INTEGER :: steps, i, iw, ii
    REAL(dp) :: grad, elocal, dpsi, gterms(3), alpha=.5d0, gpre, h=1.d-6

    OPEN(unit=1122,file="i_a_grad.txt")

    WRITE(*,*) ' Give optimization parameter (try 0.5)'
    READ(*,*) alpha

    gpre=0.d0
    ii=1
    DO i=1, steps
      CALL metropolis(10,.TRUE.)
      gterms=0.d0
      DO iw=1, Nw
        elocal=el(r(iw,:,:))
        dpsi=(psi(r(iw,:,:),a+h)/psi(r(iw,:,:),a)-1.d0)/h
        gterms(1)=gterms(1)+elocal*dpsi
        gterms(2)=gterms(2)+dpsi
        gterms(3)=gterms(3)+elocal
      END DO
      gterms=gterms/REAL(Nw)
      grad=(gterms(1)-gterms(2)*gterms(3))
      IF (grad*gpre<0.d0) THEN ! we damp only when g changes sign
        ii=ii+1
      END IF
      gpre=grad
      WRITE(1122,*) i, a, grad, ii
      a=a-alpha*grad/REAL(ii)**0.51d0
      IF (MOD(i*10,steps)==0) WRITE(out,*) i, ' grad and a ', grad, a
    END DO
    WRITE(*,*) ' Final parameter ', a
  END SUBROUTINE opt
  subroutine para_scan()
    real(dp) :: lim(2)
    integer :: i, N
    write(*,*) 'min and max parameter'
    read(*,*) lim
    write(*,*) 'number of values'
    read(*,*) N
    do i=1, N
       a=lim(1)+real(i-1)/real(N-1)*(lim(2)-lim(1))
       CALL metropolis(abs(Nsteps(1)))
    end do
    stop
  end subroutine para_scan
END PROGRAM Ho_2d
