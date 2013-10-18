      subroutine drydeprt(igrid,ncol,nrow,nlay,nrtsp,itzon,tsurf,cellat,
     &                  cellon,pwc,cwc,height,press,windu,windv,fcloud,
     &                  water,fsurf,tempk,vdep)
c
c-----CAMx v4.02 030709
c 
c     This version of is for the RTRAC species.
c     DRYDEPRT is the driver for the calculation of gridded dry deposition 
c     velocities for a given grid. Deposition velocities are calculated for
c     each gas species and for each aerosol size bin, weighted by the
c     fractional land use specified for each cell.
c 
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c           
c     Modifications:
c        5/20/03   Taken from DRYDEP for regular model
c
c     Input arguments:
c        igrid               grid index 
c        ncol                number of columns 
c        nrow                number of rows 
c        nlay                number of layers 
c        nrtsp               number of species
c        itzon               time zone
c        tsurf               surface temperature field (K)
c        cellat              cell centroid latitude (deg)
c        cellon              cell centroid longitude (deg)
c        pwc                 precipitation water content (g/m3)
c        cwc                 cloud water content (g/m3)
c        height              layer interface height field (m)
c        press               layer pressure field (mb)
c        windu               layer U-component wind field (m/s)
c        windv               layer V-component wind field (m/s)
c        fcloud              fractional cloud cover field (fraction)
c        water               layer water vapor field (ppm)
c        fsurf               fractional landuse cover field (fraction)
c        tempk               layer temperature field (K)
c             
c     Output arguments: 
c        vdep                species-dependent deposition velocity field (m/s)
c             
c     Routines called: 
c        CALDATE
c        GETZNTH
c        MICROMET
c        VD_GAS
c        VD_AER
c             
c     Called by: 
c        CAMx 
c 
      include 'camx.prm'
      include 'camx.com'
      include 'bndary.com'
      include 'deposit.com'
      include 'chmstry.com'
      include 'filunit.com'
      include 'tracer.com'
      include 'rtracchm.com'
c
      dimension tsurf(ncol,nrow),cellat(ncol,nrow),cellon(ncol,nrow)
      dimension height(ncol,nrow,nlay),press(ncol,nrow,nlay),
     &          windu(ncol,nrow,nlay),windv(ncol,nrow,nlay),
     &          fcloud(ncol,nrow,nlay),water(ncol,nrow,nlay),
     &          fsurf(ncol,nrow,NLU),tempk(ncol,nrow,nlay),
     &          pwc(ncol,nrow,nlay),cwc(ncol,nrow,nlay)
     
      dimension vdep(ncol,nrow,nrtsp)
      real lv
      logical ldark
c
      data eps/0.622/, e0/6.11/, lv/2.5e6/, rv/461./, pi/3.1415927/
c
c-----Entry point
c
c-----Determine season index
c
      idate = date
      call caldate(idate)
      month = (idate - 10000*int(idate/10000.))/100
      if (month.le.2 .or. month.eq.12) then  
        isesn = 4 
      elseif (month.ge.10) then 
        isesn = 3 
      elseif (month.eq.9) then 
        isesn = 2 
      elseif (month.ge.6) then 
        isesn = 1
      else 
        isesn = 5 
      endif 
c
c-----Loop over rows and columns
c
      do 30 j = 2,nrow-1 
        i1 = 2 
        i2 = ncol-1 
        if (igrid.eq.1) then 
          if (ibeg(j).eq.-999) goto 30 
          i1 = ibeg(j) 
          i2 = iend(j) 
        endif 
        do 20 i = i1,i2 
c
c-----Load local met variables
c
          totland = 0.
          deltaz = height(i,j,1)/2.
          temp0 = tsurf(i,j) - 273.15
          prss0 = press(i,j,1) - 
     &            2.*deltaz*(press(i,j,2) - press(i,j,1))/height(i,j,2)
          ucomp = (windu(i,j,1) + windu(i-1,j,1))/2.
          vcomp = (windv(i,j,1) + windv(i,j-1,1))/2.
          wind = sqrt(ucomp**2 + vcomp**2)
          wind = amax1(0.1,wind)
c
c-----Calculate solar flux (Holstag and Van Ulden, 1983; as implemented
c     in CALMET)
c
          call getznth(cellat(i,j),cellon(i,j),time,date,itzon,
     &                 zenith,ldark)
          sinelv = cos(zenith*pi/180.)
          solflux = (990.*sinelv - 30.)*(1. - 0.75*fcloud(i,j,1)**3.4) 
          solflux = amax1(0.,solflux)
c
c-----Determine surface wetness
c
          iwet = 0
          if (pwc(i,j,1).gt.cwmin) then
            iwet = 2
          elseif (cwc(i,j,1).gt.cwmin) then
            iwet = 1
          else
            qwatr = 1.e-6*water(i,j,1)*18./28.8
            ev = qwatr*prss0/(qwatr + eps) 
            es = e0*exp((lv/rv)*(1./273. - 1./tsurf(i,j))) 
            rh = amin1(1.,ev/es)
            dew = (1. - rh)*(wind + 0.6) 
            if (dew.lt.0.19) iwet = 1
          endif
c
c-----Loop over land use; surface roughness for water is dependent on
c     wind speed
c
          do l = 1,nrtrac
            vdep(i,j,l) = 0.
          enddo
c
          do 10 m = 1,NLU
            if (fsurf(i,j,m).lt.0.01) goto 10
            totland = totland + fsurf(i,j,m)
            z0 = z0lu(m)
            if (m.eq.7) z0 = 2.0e-6*wind**2.5 
c 
c-----Get surface layer micrometeorological parameters for this cell and
c     landuse type
c 
            if (prss0.lt.0) then
              write(iout,'(//,a)') 'ERROR in DRYDEPRT:'
              write(iout,*) 'Invalid pressure value'
              write(iout,*) 'Cell   Height  Deltaz'
              write(iout,*) i,j,height(i,j,1),deltaz
              call camxerr()
            endif
            call micromet(tempk(i,j,1),tsurf(i,j),press(i,j,1),prss0,
     &                    deltaz,wind,z0,ustar,psih) 
c
c-----Loop over GAS species, and calculate deposition velocity for this cell,
c     landuse, and current species
c
            henso2 = henso20*exp(tfactso2*(1/298. - 1./tsurf(i,j)))
            do 40 l = 1,nrtgas
              if (rthlaw(l).gt.1.e-6) then
                iflgso2 = 0
                iflgo3 = 0
                henry = rthlaw(l)*exp(rttfact(l)*(1/298.-1./tsurf(i,j)))
                call vd_gas(m,istress(m,isesn),iwet,iflgso2,iflgo3,z0,
     &                     deltaz,psih,ustar,rtdrate(l),henry,henso2,
     &                     rtreact(l),rtscale(l),temp0,solflux,
     &                     rj(m,isesn),rlu(m,isesn),rac(m,isesn),
     &                     rlcs(m,isesn),rlco(m,isesn),rgss(m,isesn),
     &                     rgso(m,isesn),vd)
              else
                vd = 0.
              endif
              vdep(i,j,l) = vdep(i,j,l) + vd*fsurf(i,j,m)
 40         continue
c
c-----Loop over AEROSOL size bins, and calculate deposition velocity for
c     this cell and landuse
c
            if (nrtaero .gt. 0) then
              do 50 l = nrtgas+1,nrtrac
                diam = sqrt(rtlcut(l)*rtucut(l))*1.e-6
                call vd_aer(z0,deltaz,psih,ustar,diam,rtdens(l),
     &                      tsurf(i,j),vd)
                vdep(i,j,l) = vdep(i,j,l) + vd*fsurf(i,j,m)
 50           continue
            endif

 10       continue
c
c-----Calculate final landuse-weighted deposition velocities
c
          do l = 1,nrtrac
            vdep(i,j,l) = vdep(i,j,l)/totland
          enddo
c
 20     continue
 30   continue
c
      return
      end
