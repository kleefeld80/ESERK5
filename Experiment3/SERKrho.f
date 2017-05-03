      subroutine SERKrho(neqn,t,f,yn,fn,iwork,hmax,work,sprad,idid)
c---------------------------------------------------------------
c  SERKRHO attempts to compute a close upper bound, SPRAD, on
c  the spectral radius of the Jacobian matrix using a nonlinear
c  power method.  A convergence failure is reported by IDID = 6.
c---------------------------------------------------------------
      integer          neqn, idid, iwork(10)
      double precision t,yn(neqn),fn(neqn),sprad
      double precision hmax,work(*)
      external         f
c
      integer          itmax
      parameter       (itmax=100)
      integer          i,iter,index,ptr5
      double precision uround,sqrtu,ynrm,sigma,sigmal,v(neqn),fv(neqn),
     &                 dynrm,dfnrm,vnrm,small
      integer          nsteps, nfesig
      parameter       (uround=2.22d-16)
c
      sqrtu = sqrt(uround)
      nfesig = iwork(9)
c------------------------------------------------------------
c  sprad smaller than small = 1/hmax are not
c  interesting because they do not constrain the step size.
c------------------------------------------------------------
      small = 1.0d0/hmax
c---------------------------------------------------------
c  The initial slope is used as guess when nsteps = 0 and
c  thereafter the last computed eigenvector.  Some care
c  is needed to deal with special cases. Approximations to
c  the eigenvector are normalized so that their Euclidean
c  norm has the constant value dynrm.
c---------------------------------------------------------
      ptr5 = 4*neqn
c     total steps = accepted steps + rejected steps = iwork(7) + iwork(8)
      nsteps = iwork(7) + iwork(8)
      if(nsteps .eq. 0) then
        do 10 i = 1,neqn
          v(i) = fn(i)
10      continue
      else
        do 20 i = 1,neqn
          v(i) = work(ptr5+i)
20      continue
      endif
      ynrm = 0.0d0
      vnrm = 0.0d0
      do 30 i = 1,neqn
        ynrm = ynrm + yn(i)**2
        vnrm = vnrm + v(i)**2
30    continue
      ynrm = sqrt(ynrm)
      vnrm = sqrt(vnrm)
      if(ynrm .ne. 0.0d0 .and. vnrm .ne. 0.0d0) then
        dynrm = ynrm*sqrtu
        do 40 i = 1,neqn
          v(i) = yn(i) + v(i)*(dynrm/vnrm)
40      continue
      elseif(ynrm .ne. 0.0d0) then
        dynrm = ynrm*sqrtu
        do 50 i = 1, neqn
          v(i) = yn(i) + yn(i)*sqrtu
50      continue
      elseif(vnrm .ne. 0.0d0) then
        dynrm = uround
        do 60 i = 1,neqn
          v(i) = v(i)*(dynrm/vnrm)
60      continue
      else
        dynrm = uround
        do 70 i = 1,neqn
          v(i) = dynrm
70      continue
      endif
c--------------------------------------------
c  Now iterate with a nonlinear power method.
c--------------------------------------------
      sigma = 0.0d0
      do 110 iter = 1, itmax
        call f(neqn,t,v,fv)
        nfesig = nfesig + 1
        dfnrm = 0.0d0
        do 80 i = 1, neqn
          dfnrm = dfnrm + (fv(i) - fn(i))**2
80      continue
        dfnrm = sqrt(dfnrm)
        sigmal = sigma
        sigma = dfnrm/dynrm
c----------------------------------------------------------
c  sprad is a little bigger than the estimate sigma of the
c  spectral radius, so is more likely to be an upper bound.
c----------------------------------------------------------
        sprad = 1.2d0*sigma
        if(iter .ge. 2 .and.
     &     abs(sigma - sigmal) .le. max(sigma,small)*0.01d0) then
          do 90 i = 1,neqn
            work(i) = yn(i)
            work(neqn+i) = fn(i)
            work(2*neqn+i) = v(i)
            work(3*neqn+i) = fv(i)
            work(ptr5+i) = v(i) - yn(i)
90        continue
c     number of evaluations of F used to estimate the spectral radius
          iwork(9) = nfesig
          return
        endif
c--------------------------------------
c  The next v(*) is the change in f
c  scaled so that norm(v - yn) = dynrm.
c--------------------------------------
        if(dfnrm .ne. 0.0d0) then
          do 100 i = 1,neqn
            v(i) = yn(i) + (fv(i) - fn(i))*(dynrm/dfnrm)
100       continue
        else
c-------------------------------------------------------
c  The new v(*) degenerated to yn(*)--"randomly" perturb
c  current approximation to the eigenvector by changing
c  the sign of one component.
c-------------------------------------------------------
          index = 1 + mod(iter,neqn)
          v(index) = yn(index) - (v(index) - yn(index))
        endif
110   continue
c     number of evaluations of F used to estimate the spectral radius
      iwork(9) = nfesig
c-------------------------------------------
c  Set flag to report a convergence failure.
c-------------------------------------------
      idid = 6
      return
      end
