! * * * * * * * * * * * * * * * * * * * * * * * * *
! --- DRIVER FOR E-SERK ON COMBUSTION 2D PROBLEM
! * * * * * * * * * * * * * * * * * * * * * * * * *
       include 'ESERK5.f'
c       include '../ESERK_withoutRK34.f'
!
! no implicit variables
!
      implicit none
!
!
! program variables
!
      integer          i, j, ns, neqn, idid, iwork(10), spcrad(2)
      parameter (ns=599,neqn=ns*ns)
      double precision y(neqn), error, truey(neqn), ans, xx, yy
      double precision t, tend, r, h, rtol, atol, work(1+5*neqn)
      external         fcombustion
!
! measure CPU time
!
      real etime          ! Declare the type of etime()
      real elapsed(2)     ! For receiving user and system time
      integer*4 now(3)

! ----- file for solution -----
      open(8,file='sol_serk5.out')
      rewind 8
      open(9,file='sol_serk5_xy.out')
      rewind 9

! ----- initial and end point of integration -----
      t=0.0d0
      tend=1.48d0
! ----- initial values -----
      ans=ns+1.0d0
      do j=1,ns
         do i=1,ns
            y((j-1)*ns+i) = 1.0d0
         end do
      end do
! ----- required tolerance -----
      rtol=1.0d-3
      atol=rtol
! ----- initial step size -----
      h=rtol
! ----- spcrad(1) = 0 ->  routine rho is provided to calculate the spectralradius -----
! ----- spcrad(1) = 1 ->  calculate the spectralradius internally using a nonlinear power method -----
! ----- spcrad(2) = 0 ->  The Jacobian may not be constant. -----
! ----- spcrad(2) = 1 ->  The Jacobian is constant. -----
      spcrad(1) = 1
      spcrad(2) = 0
c-----------------------------
c  Initialize the integration.
c-----------------------------
      idid = 0
      call itime(now)    ! now(1)=hour, (2)=minute, (3)=second
      write(6,*) 'Integration of the 2-dim Combustion problem'
      call ESERK(neqn,t,tend,h,y,fcombustion,atol,spcrad,iwork,
     *           work,idid)
      call itime(now)    ! now(1)=hour, (2)=minute, (3)=second

      write(*,*) 'Time needed', etime(elapsed)
C --- PRINT SOLUTION
      do j=1,neqn
        write (8,*) y(j)
      end do
      
c     print solution only at the points x=y
c      write (9,*) 'r    ', 'u(x,y, tend)'

      do j=1,ns
        xx=(j-1)/ans
        yy=(j-1)/ans
        r=sqrt(xx*xx+yy*yy)
        i=(j-1)*ns+j
c        write(*,*) xx, yy, j, i, r
        write (9,*) r, y(i)
      end do
c------------------------------------------------------
c  Done.  Report some statistics.
c------------------------------------------------------
! ----- print statistics -----
      write(6,*) 'Solution is tabulated in file sol_serk_xy.out'
      write(6,*) 'The value of IDID is',idid
!      write(6,*) 'Max estimation of the spectral radius=',iwork(11)
!      write(6,*) 'Min estimation of the spectral radius=',iwork(12)
      write(6,*) 'Max number of stages used=',iwork(10)
!      write(6,*) 'Number of f eval. for the spectr. radius=',iwork(9)
      write (6,91) iwork(5),iwork(7)+iwork(8),iwork(7),iwork(8)
 91   format(' Number of f evaluations=',i7,' steps=',i6,
     *        ' accpt=',i6,' rejct=',i3)
c------------------------------------------------------
c  Done.  Compute the error and report some statistics.
c------------------------------------------------------

! ----- file with exact solution -----
      open(20,file='sol_exact.out')

      read(20,*) truey
      error = 0.0d0


      do j=1,neqn
        error = max(error,abs(y(j) - truey(j)))
      end do

      close(8)
      close(20)

      write(*,'(/a,e10.4,a,f4.2,a,e10.4)') ' With tol = ',rtol,
     &  ', the maximum error at tend = ',tend,' was ',error


      close(8)
      close(9)
!--------------------------------------------------------
!     End of main program
!--------------------------------------------------------
      end


      double precision function rho(neqn,t,y)
c--------------------------
c  This is a dummy routine.
c--------------------------
      integer          neqn
      double precision t,y(neqn)

      alpha = 0.25d0
      nssq = neqn/2
      rho = 8.0d0*nssq*alpha + 2.0d0


      return
      end

