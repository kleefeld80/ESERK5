c--------------------------------------------------------
c     The subroutine FBRUS compute the value of f(x,y) and
c     has to be declared as external.
c--------------------------------------------------------
      subroutine fbrus(neqn,x,y,f)
c ----- brusselator with diffusion in 2 dim. space -----
      implicit none
      integer neqn, ns, nssq, nsnsm1, i, ix, iy
      double precision y(neqn), f(neqn), x, ans, xx, yy
      double precision uleft, vleft, uright, vright, ulow, vlow
      double precision uup, vup, uij, vij, brussa, brussb, alf
c
c include global variables and their definitions
c

      nssq = neqn/2
      ns = sqrt(DBLE(nssq))
      nsnsm1=ns*(ns-1)
    
c ----- constants for inhomogenity -----
      
      brussa=1.3d0
      brussb=2.0d6
      alf=1.d-1
      
c      write(*,*) 'time=', x
c ----- big loop -----
      do i=1,nssq
c ----- left neighbour -----
         if(mod(i,ns).eq.1)then
            uleft=y((i+ns-1)*2-1)
            vleft=y((i+ns-1)*2)
         else
            uleft=y((i-1)*2-1)
            vleft=y((i-1)*2)
         end if
c ----- right neighbour -----
         if(mod(i,ns).eq.0)then
            uright=y((i-ns+1)*2-1)
            vright=y((i-ns+1)*2)
         else
            uright=y((i+1)*2-1)
            vright=y((i+1)*2)
         end if
c ----- lower neighbour -----
         if(i.le.ns)then
            ulow=y((i+nsnsm1)*2-1)
            vlow=y((i+nsnsm1)*2)
         else
            ulow=y((i-ns)*2-1)
            vlow=y((i-ns)*2)
         end if
c ----- upper neighbour -----
         if(i.gt.nsnsm1)then
            uup=y((i-nsnsm1)*2-1)
            vup=y((i-nsnsm1)*2)
         else
            uup=y((i+ns)*2-1)
            vup=y((i+ns)*2)
         end if
c ----- the derivative -----
         uij=y(i*2-1)
         vij=y(i*2)
         f(i*2-1)=alf*nssq*(uleft+uright+ulow+uup-4.d0*uij)+
     *            brussa+uij*uij*vij-(brussb+1.d0)*uij
         f(i*2)=alf*nssq*(vleft+vright+vlow+vup-4.d0*vij)+
     *            brussb*uij - uij*uij*vij
      end do
      
      return
      end

