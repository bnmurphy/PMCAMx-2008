      subroutine emiss(igrid,kno,kno2,nspec,narspc,nptspc,
     &                 larmap,lptmap,nsrc,idsrc,isrc,jsrc,ncol,
     &                 nrow,nlay,deltat,dx,dy,mapscl,
     &                 height,depth,windu,windv,tempk,press,aremis,
     &                 pttrace,armass,ptmass,conc,ipa_cel)
c
c-----CAMx v4.02 030709
c
c     EMISS updates species concentrations due to emissions
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c          
c     Modifications: 
c        10/29/01   Map scale factor added to emissions increment
c
c     Input arguments:
c        igrid               grid numbe
c        kno                 index of NO in species arrays
c        kno2                index of NO2 in species arrays
c        nspec               number of species
c        narspc              number of area source species
c        nptspc              number of point source species
c        larmap              area source species map
c        lptmap              point source species map
c        nsrc                number of point sources
c        idsrc               point source ID map
c        isrc                grid column index for point sources
c        jsrc                grid row index for point sources
c        ncol                number of columns
c        nrow                number of rows
c        nlay                number of layers
c        deltat              time step size (s)
c        dx                  cell size in x-direction (m)
c        dy                  cell size in y-direction (m)
c        mapscl              map scale factor
c        height              layer height (m)
c        depth               layer depth (m)
c        windu               wind speed in x-direction (m/s)
c        windv               wind speed in y-direction (m/s)
c        tempk               air temprature (K)
c        press               air pressure (mb)
c        aremis              area source strength (mol/s)
c        pttrace             point source strength (mol/s)
c        armass              mass from area source emission (umol)
c        ptmass              mass from point source emission (umol)
c        conc                species concentrations (umol/m3)
c        ipa_cel             gridded array to identify if cell is
c                            in a IPRM sub-domain
c
c     Output arguments:
c        conc                species concentrations (umol/m3)
c        armass              mass from area source emission (umol)
c        ptmass              mass from point source emission (umol)
c
c     Routines called:
c        PLUMERIS
c
c     Called by:
c        EMISTRNS
c
      include  "camx.prm"
      include  "ptemiss.com"
      include  "bndary.com"
      include  "flags.com"
c
c
c======================== Process Analysis Begin ====================================
c
      include "procan.com"
      include "tracer.com"
c
      integer ipa_cel(ncol,nrow,nlay)
c
c========================= Process Analysis End =====================================
c
      real*8 armass(nspec),ptmass(nspec),dmass
      dimension larmap(narspc),lptmap(nptspc),idsrc(nsrc),isrc(nsrc),
     &          jsrc(nsrc),pttrace(MXPTSRC,nspec),
     &          conc(ncol,nrow,nlay,nspec),aremis(ncol,nrow,narspc),
     &          dx(nrow)
      real      mapscl(ncol,nrow)
      dimension height(ncol,nrow,nlay),depth(ncol,nrow,nlay),
     &          windu(ncol,nrow,nlay),windv(ncol,nrow,nlay),
     &          tempk(ncol,nrow,nlay),press(ncol,nrow,nlay)
      dimension hght1d(MXLAYA),wind1d(MXLAYA),tempk1d(MXLAYA),
     1          dtdz1d(MXLAYA)
c
      data gamma,p0 /0.286,1000./
c
c-----Entry point
c
c-----Update concentration due to area source
c
      if (larsrc) then
        do 10 lar = 1,narspc
          l = larmap(lar)
          if (l.eq.0) goto 10
          do 20 i = 2,ncol-1
c
c-----Check for boundary cell
c
            j1 = 2 
            j2 = nrow-1 
            if (igrid.eq.1) then 
              if (jbeg(i).eq.-999) goto 20 
              j1 = jbeg(i) 
              j2 = jend(i) 
            endif 
            do j = j1,j2
              vol = dx(j)*dy*depth(i,j,1)/(mapscl(i,j)**2)
              dmass = aremis(i,j,lar)*deltat*1e6
              dconc = REAL(dmass)/vol
              armass(l) = armass(l) + dmass
              conc(i,j,1,l) = conc(i,j,1,l) + dconc
c
c======================== Process Analysis Begin ====================================
c
c----- if surface layer in this column is in a sub-domain
c      then track the area emissions in this grid cell ---
c
              if( .NOT. ltrace .AND. lipr .AND. 
     &                                ipa_cel(i,j,1) .GT. 0) then
                    ipa_idx = ipa_cel(i,j,1)
                    cipr(IPR_AEMIS, ipa_idx, l) =
     &                           cipr(IPR_AEMIS, ipa_idx, l) + dconc
              endif
c
c========================= Process Analysis End =====================================
c
            enddo
 20       continue
 10     continue
      endif
c
c-----Update concentration due to point sources
c
      if (lptsrc) then
        do 50 lsrc = 1,nsrc
          n = idsrc(lsrc)
          i = isrc(lsrc)
          j = jsrc(lsrc)
c
c-----If effph is negative, override PLUMERIS
c
          if (effph(n) .lt. 0.) then
            zstk = abs(effph(n))
            goto 14
          endif
          do k = 1,nlay
            hght1d(k) = height(i,j,k)
            tempk1d(k) = tempk(i,j,k)
            w2 = windu(i,j,k)*windu(i,j,k) + windv(i,j,k)*windv(i,j,k)
            wind1d(k) = amax1(sqrt(w2),0.1)
            if (k.lt.nlay) then
              dz = height(i,j,k+1)/2. 
              if (k.gt.1) dz = (height(i,j,k+1) - height(i,j,k-1))/2. 
              dtheta = (tempk(i,j,k+1)*(p0/press(i,j,k+1))**gamma - 
     &                  tempk(i,j,k)*(p0/press(i,j,k))**gamma) 
              dtdz1d(k) = dtheta/dz 
            else
              dtdz1d(k) = dtdz1d(k-1)
            endif 
          enddo
c
c-----Calculate plume rise
c
          dstkabs = abs(dstk(n))
          call plumeris(nlay,hght1d,tempk1d,dtdz1d,wind1d,hstk(n),
     &                  dstkabs,tstk(n),vstk(n),zstk)
c
  14      continue
          do k = 1,nlay
            if (height(i,j,k).gt.zstk) goto 15
          enddo
          k = nlay
  15      continue
          do 40 lpt = 1,nptspc
            l = lptmap(lpt)
            if( ipigflg .EQ. GRESPIG .AND. lpiglet(n) .AND. 
     &                              (l.eq.kno .OR. l .EQ. kno2)) goto 40
            if( ipigflg .EQ. IRONPIG .AND. lpiglet(n) ) goto 40
            if( l .EQ. 0 ) goto 40
            vol = dx(j)*dy*depth(i,j,k)/(mapscl(i,j)**2)
            dmass = pttrace(n,lpt)*deltat*1e6
            dconc = REAL(dmass)/vol
            ptmass(l) = ptmass(l) + dmass
            conc(i,j,k,l) = conc(i,j,k,l) + dconc
c
c======================== Process Analysis Begin ====================================
c
c   --- if layer containing plume for this cell is in a sub-domain
c       then track the point source emissions for this cell ----
c
            if( .NOT. ltrace .AND. lipr .AND. 
     &                                  ipa_cel(i,j,k) .GT. 0 ) then
              ipa_idx = ipa_cel(i,j,k)
              cipr(IPR_PTEMIS, ipa_idx, l) =
     &                           cipr(IPR_PTEMIS, ipa_idx, l) + dconc
c
            endif
c
c========================= Process Analysis End =====================================
c
 40       continue
 50     continue
      endif
c
      return
      end
