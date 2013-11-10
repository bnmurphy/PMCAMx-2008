      subroutine linint(dx,yx,m,df,yf,n,ns,diamb,diame,bound,secsize,ax,af)
c
c This routine regards dx as (LOG10) midpoint of each moving section
c
c   MODIFICATIONS
c     04/15/2002 (tmg); modified algorithm for linear interpolation
c                       use known section width to interpolate mass into 
c                       fixed sections rather than using logarithmic mean
c                       between moving sectional diameters
c
c...USE:
c     THIS SUBROUTINE PERFORMS LOG10 LINEAR INTERPOLATION
c     OF THE MOVING SECTIONAL DISTRIBUTION FUNCTION yx(m)
c     WITH SECTIONAL BOUNDARIES dx(m)
c     TO THE FIXED SECTIONAL DISTRIBUTION (STEP FUNCTION)
c     yf(n) WITH SECTIONAL BOUNDARIES df(n+1)
c      
c...INPUT
c     dx(m) : DIAMETERS FOR INPUT yx MOVING SECTIONAL DISTRIBUTION
c     yx(m) : VALUES OF THE FUNCTION AT dx(m) (UNITS ARE MASS CONCENTRATION)
c     m : NUMBER OF MOVING SECTIONS IN THE INPUT DISTRIBUTION
c
c     df(n+1) : DIAMETERS SECTIONAL DISTRIBUTION FOR OUTPUT yf DISTRIBUTION
c     n : NUMBER OF SECTIONS IN THE OUTPUT DISTRIBUTION
c
c...OUTPUT
c     yf(n) : LINEAR INTERPOLATED VALUES (UNITS ARE MASS CONCENTRATION)
c
c
      !BNM - Commented out these declarations (2-23-11)
      !implicit double precision (a-h,o-z)
      !dimension dx(m+1),yx(m*ns),df(n+1),yf(n*ns)
      !dimension diamb(m), diame(m)
      !dimension bound(n+1)
      !dimension ax(ns), af(ns)
      !dimension secsize(m)
      !real*8 tempod
      !integer ifsec,isec,iupper

      !Replace with these - BNM
      implicit none
      integer i,m,ns,n,k,kk,nsp
      real*8 bound(n+1),df(n+1),dx(m+1),secsize(m)
      real*8 diamb(m),diame(m),ax(ns),af(ns)
      real*8 yx(m*ns),yf(n*ns),tempod
      integer ifsec,isec,iupper
      integer igrdchm, ichm, jchm, kchm, iout, idiag
      !End BNM

      print *,'linint: ns=',ns,'   nsp=',nsp
c     set fixed section boundaries
      bound(1)=log10(df(1))
      do i=2,n+1
         bound(i)=log10(df(i))
         secsize(i-1)=bound(i)-bound(i-1)
      enddo

c     The section width stays the same as it moves, so we can calculate
c     the new lower and upper boundaries using the moving sectional 
c     diameter as the midpoint
      do i=1,m
         diamb(i)=log10(dx(i))-secsize(i)/2.0d0
         diame(i)=log10(dx(i))+secsize(i)/2.0d0
      enddo

      do k=1,ns
       ax(k)=0.0
       do i=1,m
        ax(k)=ax(k)+yx((i-1)*ns+k)
       enddo
      enddo

      do k=1,ns
       do i=1,n
        yf((i-1)*ns+k)=0.0
       enddo
      enddo

c     sort by minimum diameter
      do i=1,m-1
         if(diamb(i).gt.diamb(i+1)) then
            tempod=diamb(i)
            diamb(i)=diamb(i+1)
            diamb(i+1)=tempod
            tempod=diame(i)
            diame(i)=diame(i+1)
            diame(i+1)=tempod
            tempod=secsize(i)
            secsize(i)=secsize(i+1)
            secsize(i+1)=tempod
            !BNM Changed all occurences nsp -> ns (2-23-11)
            !nsp = ns
            print *,'linint: nsp=',nsp
            do kk=1,ns
               tempod = yx((i-1)*ns+kk)
               yx((i-1)*ns+kk) = yx(i*ns+kk)
               yx(i*ns+kk) = tempod
            enddo
            !BNM Changed nsp -> ns (2-23-11)
         endif
      enddo

c   Start with first moving (isec) and fixed (ifsec) sections
      isec=1
      ifsec=1

c   Partition mass into fixed sections
 20   if (diamb(isec).lt.bound(ifsec+1)) then
c   Add all mass if whole moving section is within fixed section
       if (diame(isec).le.bound(ifsec+1)) then
        do k=1,ns
         yf((ifsec-1)*ns+k)=yf((ifsec-1)*ns+k)+yx((isec-1)*ns+k)
        enddo
        isec=isec+1
        if (isec.le.m) goto 20
c   Otherwise just add portion of mass within fixed section and
c   then calculate the upper sections separately
       else
        do k=1,ns
         yf((ifsec-1)*ns+k)=yf((ifsec-1)*ns+k)+yx((isec-1)*ns+k)
     &             *(bound(ifsec+1)-diamb(isec))/secsize(isec)
        enddo
c     Add remaining mass in upper sections, taking into account that
c     upper boundary is not necessarily in the next fixed section 
        iupper=ifsec+2
 30     if (diame(isec).gt.bound(iupper)) then
         if (iupper.le.n) then
          do k=1,ns
           yf((iupper-2)*ns+k)=yf((iupper-2)*ns+k)+yx((isec-1)*ns+k)
     &               *(bound(iupper)-bound(iupper-1))/secsize(isec)
          enddo
          iupper=iupper+1
          goto 30
         else
          do k=1,ns
           yf((iupper-2)*ns+k)=yf((iupper-2)*ns+k)+yx((isec-1)*ns+k)
     &               *(diame(isec)-bound(iupper-1))/secsize(isec)
          enddo
         endif
        else
         do k=1,ns
          yf((iupper-2)*ns+k)=yf((iupper-2)*ns+k)+yx((isec-1)*ns+k)
     &              *(diame(isec)-bound(iupper-1))/secsize(isec)
         enddo
        endif
        isec=isec+1
        if (isec.le.m) goto 20
       endif
      else
       ifsec=ifsec+1
       if (ifsec.lt.n) goto 20 
      endif

c     Put remaining mass in largest section
c
 40   if (isec.le.m) then
       do k=1,ns
        yf((n-1)*ns+k)=yf((n-1)*ns+k)+yx((isec-1)*ns+k)
       enddo
       isec=isec+1
       goto 40
      endif

      do k=1,ns
       af(k)=0.0
       do i=1,n
        af(k)=af(k)+yf((i-1)*ns+k)
       enddo
      enddo

c     Check for mass conservation error
      do k=1,ns
       if (abs(ax(k)-af(k))/ax(k).gt.1.0d-6.and.ax(k).gt.1.0d-10) then
          call get_param(igrdchm,ichm,jchm,kchm,iout,idiag)
          write(iout,'(//,A)')'ERROR in LININT: mass conservation'
          write(iout,*)' ax(k)-af(k)=',ax(k)-af(k),' ax(k)=',ax(k)
          write(iout,*)' igrd,i,j,k: ',igrdchm,ichm,jchm,kchm
          call camxerr()
       endif
      enddo
c
      return
      end
