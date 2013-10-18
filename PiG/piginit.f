      subroutine piginit(dt,pttrace,height,windu,windv,tempk,press)
c
c-----CAMx v4.02 030709
c
c     PIGINIT initializes PiG puffs from flagged sources
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c          
c     Modifications: 
c        07/05/02    Added code for IRON-PiG
c
c     Input arguments:
c        dt                  time step size (s)
c        pttrace             point source emissions (mol/s)
c        height              layer height (m)
c        windu               wind speed in x-direction (m/s)
c        windv               wind speed in y-direction (m/s)
c        tempk               air temprature (deg.K)
c        press               air pressure (mb)
c
c     Output arguments:
c        none
c
c     Subroutines Called:
c        PLUMERIS
c
c     Called By:
c        CAMx
c
      include  "camx.prm"
      include  "grid.com"
      include  "ptemiss.com"
      include  "pigsty.com"
      include  "chmstry.com"
      include  "flags.com"
      include  "filunit.com"
c
      dimension height(MXVEC3D),windu(MXVEC3D),
     &          windv(MXVEC3D),tempk(MXVEC3D),
     &          press(MXVEC3D),pttrace(MXPTSRC,MXSPEC)
      dimension hght1d(MXLAYA),wind1d(MXLAYA),tempk1d(MXLAYA),
     &          dtdz1d(MXLAYA),zmid(MXLAYA),press1d(MXLAYA)
c
      real*4 frctr(MXRECTR),frctrin(MXRECTR)
c
      data pi/3.1415926/
      data cpcv/1.41/
      data zscore/1.5/
      data gamma /0.286/
      data p0 /1000./
c
c-----Entry point
c
c-----Divide initial mass among reactors
c
c.........................for future reference........................
c                 if lrctin is true, reactor input fractions (frctrin)
c                 must be specified here (or somewhere).
c                 lrctin is set in the .prm block
c......................................................................
c
c-----If fractions not specified, partition according to erf function
c     correct for stopping at 2 sigma
c
      if( ipigflg .EQ. IRONPIG ) then
         if( lrctin ) then
           frsum=0.
           do nr=1,nreactr
             frctr(nr)=frctrin(nr)
             frsum=frsum + frctr(nr)
           enddo
             if( abs(frsum-1.00) .gt. .001) then
                write (iout,'(//,a)') 'ERROR in PIGINIT:'
                write (iout,*) 'Reactor Fractions Do Not Sum To 1.00'
                call camxerr()
             endif
         else
           frsum=0.
           do nr=1,nreactr
              t1=nr*2/float(nreactr)
              t2=(nr-1)*2/float(nreactr)
              frctr(nr)= erf(t1)-erf(t2)
              frsum=frsum + frctr(nr)
           enddo
c
c-----Correct for truncation of erf to conserve mass
c
           do nr=1,nreactr
              frctr(nr)=(1.000/frsum)*frctr(nr)
           enddo
         endif
      endif
c
c-----Get species maps for NO and NO2
c
      do lpt = 1, nptspc
        if (lptmap(lpt).eq.kno) lpt1 = lpt
        if (lptmap(lpt).eq.kno2) lpt2 = lpt
      enddo
c
c-----Loop over elevated point sources in Coarse Grid
c
      m1 = 1
      do 50 lsrc = 1,nosrc(1)
        n = idsrc(lsrc,1)
        if (.not.lpiglet(n)) goto 50
        if (ipigflg .EQ. GRESPIG .AND. pttrace(n,lpt1).le.0. 
     &                      .and. pttrace(n,lpt2).le.0.) goto 50
c
c-----Found a PiG point source with non-zero emissions; load local 1-D
c     vectors from finest grid
c
        i = isrc(lsrc,1)
        j = jsrc(lsrc,1)
        igrd0 = 1
        do ip = 1,ngrid
          do ic = 1,nchdrn(ip)
            igrd = idchdrn(ic,ip)
            ig = mapgrd(igrd)
            if (i.ge.inst1(ig) .and. i.le.inst2(ig) .and.             
     &          j.ge.jnst1(ig) .and. j.le.jnst2(ig)) 
     &      igrd0 = igrd
          enddo
        enddo
        if (igrd0.eq.1) then
          ii = i
          jj = j
        else
          xtmp = xstk(n,1) - (inst1(igrd0) - 1)*delx
          ytmp = ystk(n,1) - (jnst1(igrd0) - 1)*dely
          ii = 2 + INT(xtmp/delx*meshold(igrd0))
          jj = 2 + INT(ytmp/dely*meshold(igrd0))
        endif
        do k = 1,nlay(igrd0)
          n3d = ii + ncol(igrd0)*(jj - 1) + 
     &          ncol(igrd0)*nrow(igrd0)*(k - 1)
          hght1d(k) = height(iptr3d(igrd0)-1+n3d)
          tempk1d(k) = tempk(iptr3d(igrd0)-1+n3d)
          press1d(k) = press(iptr3d(igrd0)-1+n3d)
          w2 = windu(iptr3d(igrd0)-1+n3d)*windu(iptr3d(igrd0)-1+n3d) +
     &         windv(iptr3d(igrd0)-1+n3d)*windv(iptr3d(igrd0)-1+n3d)
          wind1d(k) = amax1(sqrt(w2),0.1) 
        enddo
        do k = 1,nlay(igrd0)
          zmid(k) = hght1d(k)/2.
          if (k.gt.1) zmid(k) = (hght1d(k) + hght1d(k-1))/2.
          if (k.lt.nlay(igrd0)) then 
            dz = hght1d(k+1)/2.  
            if (k.gt.1) dz = (hght1d(k+1) - hght1d(k-1))/2.  
            dtheta = (tempk1d(k+1)*(p0/press1d(k+1))**gamma -  
     &                tempk1d(k)*(p0/press1d(k))**gamma)  
            dtdz1d(k) = dtheta/dz  
          else
            dtdz1d(k) = dtdz1d(k-1)  
          endif  
        enddo
c
c-----Calculate plume rise
c
        dstktmp = abs(dstk(n))
        call plumeris(nlay(igrd0),hght1d,tempk1d,dtdz1d,wind1d,hstk(n), 
     &                dstktmp,tstk(n),vstk(n),zstk) 
        do kk = 1,nlay(igrd0)
          if (hght1d(kk).gt.zstk) goto 15
        enddo
        kk = nlay(igrd0)
        if (zstk.gt.hght1d(kk)) zstk = hght1d(kk)
  15    continue
c
c-----Calculate pressure at stack top
c
        do k0 = 1,nlay(igrd0)
          if (hght1d(k0).gt.hstk(n)) goto 10
        enddo
        k0 = nlay(igrd0)
  10    continue
        dz = hstk(n) - zmid(k0)
        if ((dz.gt.0. .and. k0.eq.nlay(igrd0)) .or.
     &      (dz.le.0. .and. k0.gt.1)) then
          dzmid = zmid(k0) - zmid(k0-1)
          pk0 = press1d(k0) + (dz/dzmid)*(press1d(k0) - press1d(k0-1))
        else
          dzmid = zmid(k0+1) - zmid(k0)
          pk0 = press1d(k0) + (dz/dzmid)*(press1d(k0+1) - press1d(k0))
        endif
c
c-----Calculate pressure, temperature, and potential temperature at top
c     of rise
c
        dz = zstk - zmid(kk)
        if ((dz.gt.0. .and. kk.eq.nlay(igrd0)) .or.
     &      (dz.le.0. .and. kk.gt.1)) then
          dzmid = zmid(kk) - zmid(kk-1)
          pk = press1d(kk) + (dz/dzmid)*(press1d(kk) - press1d(kk-1))
          tk = tempk1d(kk) + (dz/dzmid)*(tempk1d(kk) - tempk1d(kk-1))
        else
          dzmid = zmid(kk+1) - zmid(kk)
          pk = press1d(kk) + (dz/dzmid)*(press1d(kk+1) - press1d(kk))
          tk = tempk1d(kk) + (dz/dzmid)*(tempk1d(kk+1) - tempk1d(kk))
        endif
        thetak = tk*(p0/pk)**gamma
c
c-----Divide the emissions into several puffs if the length exceeds the 
c     user-defined maximum 
c
        npuff = int(wind1d(kk)*dt/xlmax - 0.1) + 1
        dtpuff = dt/float(npuff)
        dtime = dt
        do 40 ipuff = 1,npuff
          dtime = dtime - dtpuff
          dttmp = dtime + .5*dtpuff
c
c-----Initial volume; flowrate*timestep
c
          vol0 = pi*dstk(n)*dstk(n)*vstk(n)*dtpuff/4.
c
c-----Assume plume rise is an adiabatic process and ideal gas law
c     also applies:   p*v**(cp/cv) = constant
c     Assume 100% entrainment factor above adiabatic expansion
c
          adfac = (pk0/pk)**(1./cpcv)
          volpuff = vol0*(1. + adfac)
c
c-----Looking for a spot for the pig; add to the first empty location
c
          do m = m1, npig
            if (ingrd(m).eq.0) then
              m0 = m
              m1 = m0 + 1
              goto 20
            endif
          enddo
          npig = npig + 1
          m0 = npig
          m1 = npig + 1
  20      continue
          if (m0.gt.MXPIG) then
            write(iout,'(//,a)') 'ERROR in PIGINIT:'
            write(iout,*) 'PiG number exceeds maximum of ',MXPIG
            write(iout,*) 'Increase MXPIG and recompile'
            call camxerr()
          endif
c
c-----Load puff variables
c
          lnewt(m0) = .true.
          lnewg(m0) = .true.
          idpig(m0) = n
          ingrd(m0) = igrd0
          iipig(m0) = ii
          jjpig(m0) = jj
          xpig(m0) = xstk(n,1)
          ypig(m0) = ystk(n,1)
          zpig(m0) = zstk
          agepig(m0) = dttmp
          thetapig(m0) = thetak
c
c-----Calculate 3-D puff dimensions and fill with concentrations
c
          if( ipigflg .EQ. GRESPIG ) then
             xlength(m0) = amax1(100.,wind1d(kk)*dtpuff)
             axisz(m0) = sqrt(volpuff/xlength(m0)/pi)
             axisy(m0) = axisz(m0)
             sigz(m0) = axisz(m0)/zscore
             sigy(m0) = axisy(m0)/zscore
             fmspig(m0) = 1. - exp(-zscore*zscore/2.)
c
             puffmass(1,m0) = pttrace(n,lpt1)*dtpuff*1e6
             puffmass(2,m0) = pttrace(n,lpt2)*dtpuff*1e6
             puffmass(3,m0) = 0.
             puffmass(4,m0) = 0.
          else if( ipigflg .EQ. IRONPIG ) then
             xlength(m0) = amax1(100.,wind1d(kk)*dtpuff)
             axisz(m0) = sqrt(volpuff/xlength(m0))
             axisy(m0) = axisz(m0)
             sigz(m0) = axisz(m0)/2.
             sigy(m0) = axisy(m0)/2.
             fmspig(m0) = 1.
             pufftop(m0) = zpig(m0) + sigz(m0)
             puffbot(m0) = zpig(m0) - sigz(m0)
             if (puffbot(m0).lt.0.) then
               axisz(m0) = axisz(m0) + puffbot(m0)
               puffbot(m0) = 0.
             endif
             if (pufftop(m0).gt.hght1d(nlay(igrd0))) then
               axisz(m0) = axisz(m0) - 
     &                     (pufftop(m0) - hght1d(nlay(igrd0)))
               pufftop(m0) = hght1d(nlay(igrd0))
             endif
c
             do nr= 1, nreactr
c
c-----Initialize to zero
c
               do is=1,nspec
                  puffpol(is,nr,m0)= 0.
               enddo
               do lpt = 1,nptspc
                  is = lptmap(lpt)
                  if( is .GT. 0 ) puffpol(is,nr,m0) = pttrace(n,lpt) *
     &                                          dtpuff * 1e6 * frctr(nr)
               enddo
             enddo
c
c-----Initialize puff radicals
c
             do nr=1,nreactr
               do is=1,MXRADCL
                  cncradp(m0,nr,is)=1.e-9
               enddo
             enddo
          endif
  40    continue
  50  continue
c
      return
      end
