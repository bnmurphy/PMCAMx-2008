      subroutine wrtmass(igrd,idate,time,ind)
c
c-----CAMx v4.02 030709
c
c     WRTMASS writes out mass and fluxes at specified time,
c     and zeros out accumulation arrays
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c          
c     Modifications:
c        5/16/00   Modified mass output to facilitate import to a spreadsheet
c
c     Input arguments:
c        igrd                grid index
c        time                simulation time (HHMM)
c        idate               simulation date (YYJJJ)
c        ind                 code to select whether to output or initialize
c
c     Output arguments:
c        none
c
c     Subroutines called:
c        none
c
c     Called by:
c        CAMx
c
      include "camx.prm"
      include "camxfld.com"
      include "chmstry.com"
      include "filunit.com"
      include "flags.com"
c
      real*8 xincrmnt(MXSPEC*14),fluxin,fluxout
      real*8 pigdmp(MXSPEC)
      character*10 nammass(25)
      data nammass/'Start Mass','Final Mass','Surf Emiss','Pnt Emiss ',
     &             'North IN  ','North OUT ','South IN  ','South OUT ',
     &             'East IN   ','East OUT  ','West IN   ','West OUT  ',
     &             'Top IN    ','Top OUT   ','Dry Dep   ','Cloud Dep ',
     &             'BlwCld Dep','Chemistry ','Just Chem ','Aer Proc. ',
     &             'Nest Chnge','PiG Change','Net Change','Residual  ',
     &             'Mass Error'/
      data cf /1.e-6/
c
c-----Entry point
c

c-----Compute residual
c
      if (ind.eq.1) then
        do i=1,nspec
          pigdmp(i) = 0.
        enddo
        if( ipigflg .EQ. GRESPIG ) then
            pigdmp(kno)   = pigdump(1,igrd)
            pigdmp(kno2)  = pigdump(2,igrd)
            pigdmp(khno3) = pigdump(3,igrd)
            pigdmp(ko3)   = pigdump(4,igrd)
        elseif( ipigflg .EQ. IRONPIG ) then
          do isp=1,nspec
               pigdmp(isp)=pigdump(isp,igrd)
          enddo
        endif

        do ilay = 1,14     !<--BNM added for layers
          do l = 1,nspec
            fluxin = 0.
            do i = 1,9,2
              fluxin = fluxin + fluxes((i-1)*nspec*14 + (l-1)*14 + ilay,igrd) 
            enddo
            fluxout = 0.
            do i = 2,10,2
              fluxout = fluxout + fluxes(ilay + (l-1)*14 + (i-1)*nspec*14,igrd) 
            enddo
            fluxout = fluxout + fluxes(ilay+(l-1)*14+(11-1)*nspec*14,igrd) !Dry Dep
            fluxout = fluxout + fluxes(ilay+(l-1)*14+(12-1)*nspec*14,igrd) !InCloud Wet Dep
            fluxout = fluxout + fluxes(ilay+(l-1)*14+(13-1)*nspec*14,igrd) !Below Cloud Wet Dep

            xincrmnt(ilay+(l-1)*14) = fluxin + fluxout + armass(ilay+(l-1)*14,igrd) + 
     &                  ptmass(ilay+(l-1)*14,igrd) + xmschem(ilay+(l-1)*14,igrd) + 
     &                  xmsfin(ilay+(l-1)*14,igrd) + pigdmp(l)
            resid(ilay + (l-1)*14,igrd) = xmsold(ilay+(l-1)*14,igrd) + 
     &                 xincrmnt(ilay+(l-1)*14) - subxmass(ilay+(l-1)*14,igrd)
          enddo
c
c-----Write mass and fluxes at this date/hour
c
          imass = 8 + ilay
          write(imass,*)
          write(imass,'(a10,3x,a5,2a10,25(3x,a10))') 
     &       'Species','Grid','Date','Time',(nammass(i),i=1,25)
          do l = 1,nspec
            indx = ilay + (l-1)*14
            lar = 1
            if (ilay.gt.1) lar = 0
            denom = max(xmsold(indx,igrd), subxmass(indx,igrd),
     &            armass(indx,igrd), ptmass(indx,igrd),
     &            abs(xmschem(indx,igrd)), abs(xmsfin(indx,igrd)),
     &            abs(pigdmp(l)))
            do j = 1,13
              denom = max(abs(real(fluxes(indx+(j-1)*nspec*14,igrd))),denom)
            enddo
            write(imass,'(3x,a10,i5,i10.5,f10.0,25(1pe13.4))') 
     &        spname(l),igrd,idate,time,
     &        xmsold(indx,igrd)*cf,  subxmass(indx,igrd)*cf,    armass(l,igrd)*cf*lar,
     &        ptmass(indx,igrd)*cf,  (fluxes(indx+(i-1)*nspec*14,igrd)*cf,i=1,13),
     &        xmschem(indx,igrd)*cf, xmsjustchem(indx,igrd)*cf, xmspart(indx,igrd)*cf,
     &        xmsfin(indx,igrd)*cf,  pigdmp(l)*cf,
     &        xincrmnt(indx)*cf,     resid(indx,igrd)*cf,    abs(resid(indx,igrd)/denom)
          enddo
c
c-----Move current mass array to old mass array
c
          do l = 1,nspec
            indx = ilay + (l-1)*14
            xmsold(indx,igrd) = subxmass(indx,igrd)
          enddo

        enddo  !<---ilay BNM
      endif
c
c-----Zeros the mass and fluxes
c
      do ilay = 1,14
        do l = 1,nspec
          indx = ilay + (l-1)*14
          armass(indx,igrd) = 0.
          ptmass(indx,igrd) = 0.
          do i=1,13
            fluxes(indx+(i-1)*nspec*14,igrd) = 0.
          enddo
          xmschem(indx,igrd) = 0.
          xmsjustchem(indx,igrd) = 0.      !<- BNM 6/2/09
          xmspart(indx,igrd) = 0.          !<- BNM 6/2/09
          xmsfin(indx,igrd) = 0.
          pigdmp(l) = 0.
        enddo
      enddo        !<---ilay BNM
c
        do l=1,4
          pigdump(l,igrd) = 0.
        enddo
c
      return
      end
