c * * * * * * * * * * * * * * * * * * * * * * * * *
c    Driver for SERK at Brusselator-2dim problem
c * * * * * * * * * * * * * * * * * * * * * * * * *
c
c    This driver shows how to use SERK2. It solves a
c    system of ODEs resulting from the 2-dimensional space 
c    discretization of the Brusselator equations (u=u(x,y,t),v=v(x,y,t)):
c     
c    u_t=1+u^2*v-4.4*u+0.1*(u_{xx}+u_{yy})+f(x,y,t)
c    v_t=3.4*u-u^2*v+0.1*(v_{xx}+v_{yy})     for t>=0, 0<= x <= 1, 0<= y <= 1
c
c    with initial conditions 
c
c    u(x,y,0)=22*y*(1-y)^{3/2}  v(x,y,0)=27*x*(1-x)^{3/2}
c
c    and periodic boundary conditions
c
c    u(x+1,y,t)=u(x,y,t),  v(x,y+1,t)=v(x,y,t).
c
c    The function f is defined by (inhomogeneity)
c
c    f(x,y,t)=5 if (x-0.3)^2+(y-0.6)^2<= 0.1^2 and t>=1.1
c            =0 else
c
c    We discretize the space variables with
c    x_i=i/(N+1), y_i=i/(N+1) for i=0,1,...,N,
c    with N=160. We obtain a system of 51200
c    equations. The spectral radius of the Jacobian 
c    can be estimated with the Gershgorin theorem   
c    (409602 is an estimation for it). Thus we
c    provide an external function RHO, giving 
c    the spectral radius of the Jacobian matrix.
c    As output point we choose t_out=1.5
c
c--------------------------------------------------------
      include 'ESERK5.f'

      implicit none
      integer i, j, k, ns, nssq, neqn, idid, iwork(10), spcrad(2)
      parameter (ns=400,neqn=ns*ns*2)
      double precision y(neqn), ans, xx, yy, rtol, atol, h, t, tend
      double precision work(1+5*neqn), truey(neqn), error
      external  fbrus
!
! measure CPU time
!
      real etime          ! Declare the type of etime()
      real elapsed(2)     ! For receiving user and system time
      integer*4 now(3)

! ----- file for solution -----
      open(8,file='sol_E_SERK5.out')
      rewind 8
! ----- dimensions -----
      nssq=ns*ns
! ----- initial and end point of integration -----
      t=0.0d0
      tend=1.0d0
! ----- initial values -----
      ans=ns
      do j=1,ns
        yy=(j-1)/ans
        do i=1,ns
          y(((j-1)*ns+i)*2-1)=22.d0*yy*(1.d0-yy)**(1.5d0)
        end do
      end do
      do i=1,ns
        xx=(i-1)/ans
        do j=1,ns
           y(((j-1)*ns+i)*2)=27.d0*xx*(1.d0-xx)**(1.5d0)
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
      spcrad(2) = 1
! ----- integration -----
      call itime(now)    ! now(1)=hour, (2)=minute, (3)=second
      write(6,*) 'Integration of the 2-dim Brusselator problem'
! ----- call of the subroutine SERK -----
      call ESERK(neqn,t,tend,h,y,fbrus,atol,spcrad,iwork,work,idid)
      call itime(now)    ! now(1)=hour, (2)=minute, (3)=second
      write(*,*) 'Time needed', etime(elapsed)
! ----- print solution -----
      !do j=1,neqn,11
      !  write (8,101) y(j)
      !end do
      !101   format(1X,F22.16)
! ----- print statistics -----
      write(6,*) 'Solution is tabulated in file sol_E_SERK.out'
      write(6,*) 'The value of IDID is',idid
!      write(6,*) 'Max estimation of the spectral radius=',iwork(11)
!      write(6,*) 'Min estimation of the spectral radius=',iwork(12)
      write(6,*) 'Max number of stages used=',iwork(10)
!      write(6,*) 'Number of f eval. for the spectr. radius=',iwork(9)
      write (6,91) iwork(5),iwork(7)+iwork(8),iwork(7),iwork(8)
 91   format(' Number of f evaluations=',i7,' steps=',i6,
     *        ' accpt=',i6,' rejct=',i3)
       
       write(8,*) ((y(((j-1)*ns+i)*2-1),i=1,ns),j=1,ns)
       write(8,*) ((y(((j-1)*ns+i)*2),  i=1,ns),j=1,ns)

			 write(6,*) 'Solution is tabulated in file sol.dat'
			 
c------------------------------------------------------
c  Done.  Compute the error and report some statistics.
c------------------------------------------------------
      open(20,file='sol_E_SERK5exact.out')
      read(20,*) truey
      error = 0.0d0
      k=1
      do j=1,ns
         do i=1,ns
            error = max(error,abs(y(((j-1)*ns+i)*2-1) - truey(k)))
           error = max(error,abs(y(((j-1)*ns+i)*2) - truey(ns*ns+k)))
            k=k+1;
         end do
      end do
      write(*,'(/a,e10.4,a,f4.2,a,e10.4)') ' With tol = ',rtol,
     &  ', the maximum error at tend = ',tend,' was ',error
      !close(8)
      close(20)
      !error = 0.0d0

      !do j=1,neqn
      !  error = max(error,abs(y(j) - truey(j)))
      !end do

      close(8)

      !write(*,'(/a,d8.1,a,f4.1,a,e8.2)') ' With rtol = atol =',rtol,
      !&  ', the maximum error at tend =',tend,' was ',error

!--------------------------------------------------------
!     End of main program
!--------------------------------------------------------
      end

!--------------------------------------------------------
!     The subroutine RHO gives an estimation of the spectral
!     radius of the Jacobian matrix of the problem. This
!     is a bound for the whole interval and thus RHO is called
!     once.
!--------------------------------------------------------
      double precision function rho(neqn,t,y)
      integer neqn
      double precision t, y(neqn), alpha
      integer nssq

      alpha = 0.25d0
      nssq = neqn/2
      rho = 8.0d0*nssq*alpha + 2.0d0

      return
      end 

