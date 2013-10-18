      subroutine gresdriv(igrd,iptr2d,ncol,nrow,nlay,itzon,
     &                    dt,dx,dy,mapscl,height,rkv,tempk,press,water,
     &                    windu,windv,cldtrns,fcloud,cellat,cellon,
     &                    conc,pigdump,ipsa3d,ipa_cel)
c
c-----CAMx v4.03 031205
c
c     GRESDRIV is the driver program for the GREASD PiG submodel; performs
c     puff growth, chemistry, entrainment and mass dumping
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        5/17/00   Expanded PIGDUMP array to include ozone increment
c        9/27/00   Added local adjustments to photolysis rate constants
c       10/29/01   Map scale factor added to PiG emissions increment
c       01/22/02   Added code to only do OSAT PiG if not using RTRAC
c        7/05/02   Changed name of routine to accommadate IRON-PiG
c        4/2/03    Removed option for UAM-V type cloud adjustment
c
c     Input arguments:
c        igrd                grid index 
c        iptr2d              pointers into vectors for 2-D fields
c        ncol                number of columns
c        nrow                number of rows
c        nlay                number of layers
c        itzon               time zone 
c        dt                  time step (s)
c        dx                  cell size in x-direction (m)
c        dy                  cell size in y-direction (m)
c        height              gridded layer height (m)
c        rkv                 gridded vertical diffusivity (m2/s)
c        tempk               gridded temperature (K)
c        press               gridded pressure (mb)
c        water               gridded water vapor concentration (ppm)
c        windu               gridded x-component windspeed (m/s)
c        windv               gridded y-component windspeed (m/s)
c        cldtrns             gridded energy transmission coefficient (fraction)
c        fcloud              gridded cloud coverage (fraction)
c        cellat              gridded cell centroid latitude (deg)
c        cellon              gridded cell centroid longitude (deg)
c        conc                species concentration in the grid (umol/m3)
c        ipsa3d              pointer into source apportionment vectors
c        ipa_cel             gridded array to identify if cell is
c                            in a IPRM sub-domain
c
c     Output arguments:
c        conc                species concentration in the grid (umol/m3)
c        pigdump             PiG dumped mass (umol)
c
c     Subroutines Called:
c        GRESGROW 
c        KTHERM
c        GETZNTH
c        KPHOTO
c        GRESCHEM
c
c     Called by:
c        EMISTRNS
c
      include "camx.prm"
      include "camx.com"
      include "filunit.com"
      include "pigsty.com"
      include "chmstry.com"
      include "ahomap.com"
c
c======================== Process Analysis Begin ====================================
c
      include "procan.com"
c
      integer ipa_cel(ncol,nrow,nlay)
c
c========================= Process Analysis End =====================================
c
c
c======================== Source Apportion Begin =======================
c
      include 'tracer.com'
c
c========================= Source Apportion End ========================
c
      real*8 pigdump(MXSPEC,MXGRID)
      logical ldark
      dimension height(ncol,nrow,nlay),rkv(ncol,nrow,nlay),
     &          tempk(ncol,nrow,nlay),water(ncol,nrow,nlay),
     &          press(ncol,nrow,nlay),windu(ncol,nrow,nlay),
     &          windv(ncol,nrow,nlay),cldtrns(ncol,nrow,nlay),
     &          fcloud(ncol,nrow,nlay),cellat(ncol,nrow),
     &          cellon(ncol,nrow),conc(ncol,nrow,nlay,nspec),dx(nrow)
      real mapscl(ncol,nrow)
      dimension concamb(3),conpig(4),dumpmass(3),pigrk(9)
c
      data pi/3.1415926/
c
c-----Entry point
c
      do 50 n = 1,npig
c
c-----Skip puffs that are not in the current grid
c
        if (ingrd(n).ne.igrd) goto 50
c
c-----Perform growth/chemistry once for new puffs over time = agepig
c     (skip growth/chemsitry for new puffs in subsequent multiple
c     fine grid steps) 
c
        if (lnewg(n)) then
          delt = agepig(n) 
        elseif (lnewt(n)) then
          goto 50
c
c-----Perform growth/chemistry for old puffs over time = delta t
c     for current grid
c
        else
          delt = dt 
        endif 
c
c-----Locate the pig in the grid
c
        i = iipig(n)
        j = jjpig(n)
        do k = 1,nlay
          if (height(i,j,k).gt.zpig(n)) goto 15
        enddo
        k = nlay
  15    continue
c
c-----Watch the pig grow in size, digest ambient or dump excess
c
        depth = height(i,j,k)
        if (k.gt.1) depth = height(i,j,k) - height(i,j,k-1)
        axiszmx = depth
        axisymx = amin1(dy,dx(j))
        rkpigv = rkv(i,j,k)
        if (k.gt.1) rkpigv = (rkv(i,j,k) + rkv(i,j,k-1))/2.
        concamb(1) = conc(i,j,k,kno)
        concamb(2) = conc(i,j,k,kno2)
        concamb(3) = conc(i,j,k,ko3)
        uu = 0.5*(windu(i,j,k) + windu(i+1,j,k))
        vv = 0.5*(windv(i,j,k) + windv(i,j+1,k))
c
        call gresgrow(n,delt,rkpigv,axiszmx,axisymx,concamb(3),
     &               dumpmass)
c
        volpuff = xlength(n)*axisy(n)*axisz(n)*pi
        volcell = depth*dx(j)*dy/(mapscl(i,j)**2)
        do l = 1,4
          conpig(l) = puffmass(l,n)/(volpuff + 1.e-6)
          if (conpig(l).lt.0.0) then
            write(iout,'(//,a)') 'ERROR in GRESDRIV:'
            write(iout,*) 'Negative concentration in pig'
            write(iout,*) 'Puff#,Grid#,length,axisx/y/z:'
            write(iout,*) n,ingrd(n),xlength(n),axisz(n),axisy(n),
     &                    sigz(n)
            write(iout,*) 'Conc and mass for four species:'
            write(iout,'(2e12.4)') (conpig(l1),puffmass(l1,n),l1=1,4)
            call camxerr()
          endif
        enddo
c
c-----Skip dumping if it's a new puff
c
        if (lnewg(n)) goto 40
c
c-----Slaughter the pig if its NOx concentration leads to a grid NOx
c     increment of less than 0.1 ppb
c
        cnox = 1000.*(puffmass(1,n) + puffmass(2,n))/
     &       (volcell*densfac*(273./tempk(i,j,k))*(press(i,j,k)/1013.))
        if (cnox.le.0.1) then
          ingrd(n) = 0
          do l = 1,3
            dumpmass(l) = dumpmass(l) + puffmass(l,n)
          enddo
        endif
c
c-----Calculate dumped mass and update ambient no and no2 concentrations
c
        pigdump(1,igrd) = pigdump(1,igrd) + dumpmass(1)
        pigdump(2,igrd) = pigdump(2,igrd) + dumpmass(2)
        pigdump(3,igrd) = pigdump(3,igrd) + dumpmass(3)
        conc(i,j,k,kno)  = conc(i,j,k,kno) + dumpmass(1)/volcell
        conc(i,j,k,kno2) = conc(i,j,k,kno2) + dumpmass(2)/volcell
        conc(i,j,k,khno3) = conc(i,j,k,khno3) + dumpmass(3)/volcell
c
c======================== Process Analysis Begin ====================================
c
        if( lipr .AND. ipa_cel(i,j,k) .GT. 0 ) then
           ipa_idx = ipa_cel(i,j,k)
c
c-----NO change from Pig
c
           cipr(IPR_PIGEMIS, ipa_idx, kno) =
     &           cipr(IPR_PIGEMIS, ipa_idx, kno) + dumpmass(1)/volcell
c
c-----NO2 change from Pig
c
           cipr(IPR_PIGEMIS, ipa_idx, kno2) =
     &           cipr(IPR_PIGEMIS, ipa_idx, kno2) + dumpmass(2)/volcell
c
c-----HNO3 change from Pig
c
           cipr(IPR_PIGEMIS, ipa_idx, khno3) =
     &           cipr(IPR_PIGEMIS, ipa_idx, khno3) + dumpmass(3)/volcell
        endif
c
c========================= Process Analysis End =====================================
c
c
c======================== Source Apportion Begin =======================
c
c  --- if doing source apportionment, put NOx mass into tracer arrays ----
c
        if( ltrace .AND. tectyp .NE. RTRAC ) then
            dumpnox = dumpmass(1) + dumpmass(2)
            dumpsa = dumpnox / volcell
            call pigsa(ncol,nrow,nlay,ntotsp,ptconc(ipsa3d),
     &                                               i,j,k,n,dumpsa)
        endif
c
c========================= Source Apportion End ========================
c
c-----Skip chemistry if its a slaughtered puff
c
        if (ingrd(n).eq.0) goto 50
c 
c-----Compute rate constants for given PiG position
c 
 40     continue
        call ktherm(tempk(i,j,k),press(i,j,k)) 
        ij = i + (j-1)*ncol
        iozon = icdozn(iptr2d-1+ij)
        ihaze = icdhaz(iptr2d-1+ij)
        ialb  = icdalb(iptr2d-1+ij)
        hght = height(i,j,k)/2000.
        if (k.gt.1) hght = (height(i,j,k) + height(i,j,k-1))/2000.
        if (cldtrns(i,j,k).ne.1.) then
          iabov = 0 
          ctrns = cldtrns(i,j,k)
          fcld = fcloud(i,j,k)
        else
          iabov = 1
          ctrns = cldtrns(i,j,1)
          fcld = fcloud(i,j,1)
        endif
        call getznth(cellat(i,j),cellon(i,j),time,date,itzon,
     &               zenith,ldark)
        call kphoto(iozon,ialb,ihaze,hght,zenith,fcld,
     &              ctrns,ldark,iabov)
c
c-----Perform chemistry
c
        convfac = densfac*(273./tempk(i,j,k))*(press(i,j,k)/1013.)
        do l = 1,9
          pigrk(l) = rk(ipigrxn(l))
        enddo
        call greschem(delt,ldark,bdnl(ko3),bdnl(kno),water(i,j,k),
     &               pigrk,convfac,conpig,o3dif) 
        do l = 1,4 
          puffmass(l,n) = conpig(l)*volpuff 
        enddo 
        pigdump(4,igrd) = pigdump(4,igrd) + o3dif*volpuff
        conc(i,j,k,ko3) = conc(i,j,k,ko3) + o3dif*volpuff/volcell
        conc(i,j,k,ko3) = amax1(conc(i,j,k,ko3),convfac*bdnl(ko3))
c
c======================== Process Analysis Begin ====================================
c
        if( lipr .AND. ipa_cel(i,j,k) .GT. 0 ) then
           ipa_idx = ipa_cel(i,j,k) 
c
c-----Ozone change from Pig
c
           cipr(IPR_PIGEMIS, ipa_idx, ko3) =
     &          cipr(IPR_PIGEMIS, ipa_idx, ko3) + o3dif*volpuff/volcell
           cipr(IPR_FINAL, ipa_idx, kno) = conc(i,j,k,kno)
           cipr(IPR_FINAL, ipa_idx, kno2) = conc(i,j,k,kno2)
           cipr(IPR_FINAL, ipa_idx, khno3) = conc(i,j,k,khno3)
           cipr(IPR_FINAL, ipa_idx, ko3) = conc(i,j,k,ko3)
c
c-----Ozone change from boundary condition
c
           if( conc(i,j,k,ko3) .LT. convfac*bdnl(ko3) ) then
               cipr(IPR_PIGEMIS, ipa_idx, ko3) =
     &            cipr(IPR_PIGEMIS, ipa_idx, ko3) + 
     &                            (convfac*bdnl(ko3)-conc(i,j,k,ko3))
           endif
        endif
c
c========================= Process Analysis End =====================================
c
c
c======================== Source Apportion Begin =======================
c
c  --- if doing source apportionment, put O3 mass into tracer arrays ----
c
        if( ltrace .AND. tectyp .NE. RTRAC ) then
            dumpsa = o3dif * volpuff / volcell
            zero = 0.
            call adjstsa(igrd,i,j,k,dumpsa,zero,zero)
        endif
c
c========================= Source Apportion End ========================
c
c-----Set LNEWG to false
c
        lnewg(n) = .false.
c
  50  continue
c
      return
      end
