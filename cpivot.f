      subroutine cpivot(nbeg,nend,ndim,a,b,ierr)
c
c-----CAMx v4.02 030709
c
c     CPIVOT solves the equation, using column pivot method when number
c     of variables are greater than 3, or analytically otherwise
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        none
c
c     Input arguments:
c        nbeg                first variable to be solved
c        nend                last variable to be solved
c        ndim                the size of the matrix
c        a                   coefficient matrix
c        b                   right hand side
c
c     Output arguments:
c        b                   solution
c        ierr                non-zero if determinant zero
c
c     Routines called:
c        none
c
c     Called by:
c        TRAP
c
      dimension a(ndim,ndim),b(ndim)
      data eps/1.0e-10/
c
c-----Entry point
c
      ierr=0
      if (nend-nbeg+1.eq.2) then
        det = a(nbeg,nbeg)*a(nend,nend)-a(nend,nbeg)*a(nbeg,nend)
        if (det.eq.0.0) then
           ierr=-1
           return
        endif
        bnbeg = b(nbeg)*a(nend,nend)-b(nend)*a(nbeg,nend)
        bnend = a(nbeg,nbeg)*b(nend)-a(nend,nbeg)*b(nbeg)
        b(nbeg) = bnbeg/det
        b(nend) = bnend/det
        return
      endif
c
      if (nend-nbeg+1.eq.3) then
        det   = a(nbeg,nbeg)    *a(nbeg+1,nbeg+1)*a(nend,nend)
     &        + a(nbeg+1,nbeg)  *a(nend,nbeg+1)  *a(nbeg,nend)
     &        + a(nend,nbeg)    *a(nbeg+1,nend)  *a(nbeg,nbeg+1)
     &        - a(nend,nbeg)    *a(nbeg+1,nbeg+1)*a(nbeg,nend)
     &        - a(nend,nbeg+1)  *a(nbeg+1,nend)  *a(nbeg,nbeg)
     &        - a(nend,nend)    *a(nbeg+1,nbeg)  *a(nbeg,nbeg+1)
        if (det.eq.0.0) then
           ierr=-2
           return
        endif
        bnbeg = b(nbeg)         *a(nbeg+1,nbeg+1)*a(nend,nend)
     &        + b(nbeg+1)       *a(nend,nbeg+1)  *a(nbeg,nend)
     &        + b(nend)         *a(nbeg+1,nend)  *a(nbeg,nbeg+1)
     &        - b(nend)         *a(nbeg+1,nbeg+1)*a(nbeg,nend)
     &        - a(nend,nbeg+1)  *a(nbeg+1,nend)  *b(nbeg)  
     &        - a(nend,nend)    *b(nbeg+1)       *a(nbeg,nbeg+1)
        bnbeg1= a(nbeg,nbeg)    *b(nbeg+1)       *a(nend,nend)
     &        + a(nbeg+1,nbeg)  *b(nend)         *a(nbeg,nend)
     &        + a(nend,nbeg)    *a(nbeg+1,nend)  *b(nbeg)  
     &        - a(nend,nbeg)    *b(nbeg+1)       *a(nbeg,nend)
     &        - b(nend)         *a(nbeg+1,nend)  *a(nbeg,nbeg)
     &        - a(nend,nend)    *a(nbeg+1,nbeg)  *b(nbeg)
        bnend = a(nbeg,nbeg)    *a(nbeg+1,nbeg+1)*b(nend)  
     &        + a(nbeg+1,nbeg)  *a(nend,nbeg+1)  *b(nbeg)
     &        + a(nend,nbeg)    *b(nbeg+1)       *a(nbeg,nbeg+1)
     &        - a(nend,nbeg)    *a(nbeg+1,nbeg+1)*b(nbeg)
     &        - a(nend,nbeg+1)  *b(nbeg+1)       *a(nbeg,nbeg)
     &        - b(nend)         *a(nbeg+1,nbeg)  *a(nbeg,nbeg+1)
 
        b(nbeg) = bnbeg/det
        b(nbeg+1) = bnbeg1/det
        b(nend) = bnend/det
c
        return
      endif
c
      do 100 icurrnt=nbeg,nend-1
c
c-----Find the maximum element in the column
c
         alarge=0.0
         do i=icurrnt,nend
           if (alarge.lt.abs(a(i,icurrnt))) then
             alarge=abs(a(i,icurrnt))
             m=i
           endif
         enddo
c
c-----Switch the rows
c
         do j=icurrnt,nend
           temp=a(icurrnt,j)
           a(icurrnt,j)=a(m,j)
           a(m,j)=temp
         enddo
         temp=b(icurrnt)
         b(icurrnt)=b(m)
         b(m)=temp
c
c-----Check for singularity
c
         if (abs(a(icurrnt,icurrnt)).le.eps) then
           ierr=icurrnt
           return
         endif
c
c-----Gaussian reduction
c
        b(icurrnt)=b(icurrnt)/a(icurrnt,icurrnt)
        do j=nend,icurrnt,-1
          a(icurrnt,j)=a(icurrnt,j)/a(icurrnt,icurrnt)
        enddo
        do k=icurrnt+1,nend
          b(k)=b(k)-b(icurrnt)*a(k,icurrnt)
          do j=nend,icurrnt,-1
            a(k,j)=a(k,j)-a(icurrnt,j)*a(k,icurrnt)
          enddo
        enddo
c
  100 continue
c
c-----Back substitution
c
      b(nend) = b(nend)/a(nend,nend)
      do i=nend-1,nbeg,-1
        do j=i+1,nend
          b(i)=b(i)-a(i,j)*b(j)
        enddo
      enddo
c
      return
      end
