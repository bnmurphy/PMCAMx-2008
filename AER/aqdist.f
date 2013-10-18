c   program aqdist

      subroutine aqdist(fdistx, fdistx2)
c
c-----PMCAMx v3.01 020531
c
c --- Calculate the factors for allocating aerosol production
c     in the aqueous chemistry to the PMCAMx size sections.
c     A coarse mode (2.5 um) and a fine mode (0.4 um) are
c     formed.
c
c     Copyright 2002
c     ENVIRON International Corporation
c
c     Modifications:
c        none
c
c     Input arguments:
c
c        nsect  -  number of size sections
c        cut    -  cut points for the size sections
c
c     Output arguments:
c
c        fdist  -  fractions for bulk AQ chemistry
c        fdist2 -  fractions for size resolved AQ chemistry
c
c     Routines Called:
c
c        AREAG
c
c     Called by:
c    
c      include 'PARACHEM.INC'
      include 'dynamic.inc'
      include 'aerpar.inc'
      include 'droppar.inc'
      include 'dropcom.inc'
c      include 'discf.inc'
      
c      parameter (nsect = 10)     ! delete for PMCAMx

      double precision cut(nsect+1), fdistx(nsect), fdistx2(nsect)
      real*8 cuts(nsect+1)
      double precision x(nsect+1),sumx,sumx1,sumx2
      double precision dfin, dcrs, sfin, scrs
      double precision fdistx2a(nsect), fdistx2b(nsect)
      integer i,actsect,noactsect

      actsect = 0
      noactsect = 0

      do i = 1,nsect
      if (daer(i) .gt. dactiv) then
      actsect = actsect +1  ! activated sections
      else
      noactsect = noactsect+1  !interstitial aerosol
      endif
      enddo
      
      do i = 1,nsect+1
      cut(i) = SNGL(dsecf(i))   ! cut = the cut points for sections
      enddo

c --- These data define the fine and coarse modes
c     of the aqueous phase aerosol production.

      data dfin /0.4d0/
      data sfin /1.8d0/
      data dcrs /2.5d0/
      data scrs /2.15d0/

c --- Entry point

      do i=1,nsect+1
        x(i) = cut(i)
      enddo
c      x(1) = 1.0d-6
c      x(nsect+1) = 100.d0

c --- calculate the fraction of Guassian in each section
c     working in log space for a log-normal distribution

      do i = 1,nsect
      fdistx2(i) = 0.d0
      fdistx(i) = 0.d0
      sumx = 0.d0
      enddo
      
      do i=noactsect+1,nsect
        ffin = areag(dlog10(x(i)),dlog10(x(i+1)),
     &                             dlog10(dfin),dlog10(sfin))
        fcrs = areag(dlog10(x(i)),dlog10(x(i+1)),
     &                             dlog10(dcrs),dlog10(scrs))
        fdistx(i) = ffin + fcrs
	fdistx2a(i) = ffin
	fdistx2b(i) = fcrs
	
      enddo 
      
c Normalize so fdist sums to 1      
      do i =1,nsect
      sumx = sumx + fdistx(i)
      enddo
      
      do i =1,nsect
      fdistx(i) = fdistx(i)/sumx
      enddo
      
      sumx1 = 0.d0
      sumx2 = 0.d0
      do i = isect,nsect
      if(daer(i) .gt. dactiv .and. daer(i) .lt. dsep) then
      sumx1 = sumx1+fdistx2a(i)
      else
      sumx2 = sumx2+fdistx2b(i)
      endif
      enddo
      
      do i = isect,nsect
      if(daer(i) .gt. dactiv .and. daer(i) .lt. dsep) then
      fdistx2(i) = fdistx2a(i)/sumx1
      else
      fdistx2(i) = fdistx2b(i)/sumx2
      endif
      enddo                 

      return
      end

