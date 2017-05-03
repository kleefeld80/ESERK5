c--------------------------------------------------------
c     The subroutine FCOMBUSTION compute the value of f(x,y) and
c     has to be declared as external.
c--------------------------------------------------------
      subroutine fcombustion(neqn,x,y,f)
      implicit none
      integer neqn, ns, nsnsm1, i, ix, iy
      double precision y(neqn), f(neqn), x, ans, xx, yy
      double precision uxx, uyy, uij
      double precision d, alpha, delta, R, oneoverh2, par
c
c include global variables and their definitions
c

      ns = sqrt(DBLE(neqn))
      nsnsm1=ns*(ns-1)
      oneoverh2 = (ns+1.0d0)*(ns+1.0d0)
c      write(*,*) ns
c ----- constants for inhomogenity -----
c      ans=ns
      d=2.5d0
      alpha=1.0d0
      delta=20.0d0
      R=5.0d0
      par=1.0d0
c ----- big loop -----
      do i=1,neqn
c         iy=(i-1)/ns+1
c         ix=i-(iy-1)*ns
c         yy=iy/ans
c         xx=ix/ans
c         write(*,*) 'ix=', ix, 'iy=', iy
c ----- left and right borders -----  x=0 or x=1
         if(mod(i,ns).eq.1) then                                        ! left border: hom. Neumann  u_0 = (4*u_1-u_2)/3
            uxx=( -2*y(i)+2*y(i+1) )/3.0d0
c            write(*,*) 'left x=0'
         else if(mod(i,ns).eq.0) then                                   ! right border: Dirchlet = 1
            uxx=y(i-1)-2*y(i)+1.0d0
c            write(*,*) 'right x=1'
         else                                                           ! inside: centered diff.
            uxx=y(i-1)-2*y(i)+y(i+1)
c            write(*,*) 'inside x'
         end if
c ----- lower and upper borders -----   y=0 or y=1
         if(i.le.ns) then                                               ! lower border: hom. Neumann  u_0 = (4*u_1-u_2)/3
            uyy=( -2*y(i)+2*y(i+ns) )/3.0d0
c            write(*,*) 'lower y=0'
         else if(i.gt.nsnsm1) then                                      ! upper border: Dirchlet = 1
            uyy=y(i-ns)-2*y(i)+1.0d0
c            write(*,*) 'upper y=1'
         else                                                           ! inside: centered diff.
            uyy=y(i-ns)-2*y(i)+y(i+ns)
c            write(*,*) 'inside y'
         end if
c         pause
c ----- the derivative -----
         uij=y(i)
         f(i)=d*oneoverh2*(uxx+uyy)+R/(alpha*delta)*(1.0d0+alpha-uij)
     *         *exp(delta*(1.0d0-par/uij))
c         write(*,*) x, f(i)
      end do
c      pause
      return
      end
