      subroutine kphoto(iozon,ialb,ihaze,hght,zenith,fcloud,cldtrns,
     &                  ldark,iabov)
c 
c-----CAMx v4.02 030709
c 
c     KPHOTO adjusts all photolysis rate constants for height, ozone
c     column, surface albedo, haze turbidity, and zenith angle
c                           
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c           
c     Modifications: 
c        4/2/03        Removed option for UAM-V type cloud adjustment
c 
c     Input arguments: 
c        iozon               ozone column index
c        ialb                surface albedo index
c        ihaze               haze turbidity index
c        hght                height AGL (km)
c        zenith              solar zenith angle (deg)
c        fcloud              cloud coverage
c        cldtrns             enery transmission coefficient
c        ldark               darkness flag (T=dark)
c        iabov               Flag to note above cloud (1=above)
c 
c     Output arguments: 
c        none
c 
c     Routines Called: 
c        none
c 
c     Called by: 
c        CHEMDRIV
c        RADDRIVR
c
      include "camx.prm"
      include "chmstry.com"
      logical ldark
c
      data deg2rad /0.01745329/
c
c-----Entry point
c
c-----If dark, zero all the photolysis reaction rate constants
c
      if (ldark) then
        do irxn = 1,nphot1
          rk(idphot1(irxn)) = 0.
        enddo
        do irxn = 1,nphot2
          rk(idphot2(irxn)) = 0.
        enddo
        goto 999
      endif       
c
c-----Locate the point in height axis
c
      ihght = 1
      do i = 1,NHGHT-1
        if (hght.ge.htint(i)) ihght = i
      enddo
      tmp = (hght - htint(ihght))/(htint(ihght+1) - htint(ihght))
      w1 = 1. - tmp
c
c-----If height greater than uppermost height in the table
c
      if (tmp.gt.1.) then
        ihght = NHGHT - 1
        w1 = 0.
      endif
c
c-----Locate the point in zenith angle axis
c
      izen = 1
      do j = 1,NZEN-1
        if (zenith.ge.zenint(j)) izen = j
      enddo
      tmp = (zenith - zenint(izen))/(zenint(izen+1) - zenint(izen))
      w2 = 1. - tmp
c
c-----If zenith is greater than the uppermost zenith in the table
c
      if (tmp.gt.1.) then
        izen = NZEN-1
        w2 = 0.
      endif
c
c-----Set the primary photolysis rates using interpolation
c     in height and zenith angle
c
      do irxn = 1,nphot1
        prkn11 = prkn(izen,  irxn,ihght,  ihaze,ialb,iozon)
        prkn12 = prkn(izen,  irxn,ihght+1,ihaze,ialb,iozon)
        prkn21 = prkn(izen+1,irxn,ihght,  ihaze,ialb,iozon)
        prkn22 = prkn(izen+1,irxn,ihght+1,ihaze,ialb,iozon)
        tmp1 = w1*prkn11 + (1.-w1)*prkn12
        tmp2 = w1*prkn21 + (1.-w1)*prkn22
        rk(idphot1(irxn)) = w2*tmp1 + (1.-w2)*tmp2
      enddo
c
c-----Cloud coverage adjustment
c
      zenang = amin1(zenith,60.0)
      zenang = deg2rad*zenang
      if (iabov.eq.1) then
        cldrat = 1. + (1. - cldtrns)*cos(zenang)
      else
        cldrat = 1.6*cldtrns*cos(zenang)
      endif
      coefcld = 1. + fcloud*(cldrat - 1.)
c
      do irxn = 1,nphot1
        rk(idphot1(irxn)) = rk(idphot1(irxn))*coefcld
      enddo
c
c-----Set secondary photolysis rates as ratios to primary rates
c
      do irxn = 1,nphot2
        rk(idphot2(irxn)) = rk(idphot3(irxn))*phtscl(irxn)
      enddo
c
 999  return
      end
