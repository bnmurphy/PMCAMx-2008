      subroutine updtmet(igrd,ncol,nrow,nlay,ngas,densfac,deltat,phpt,
     &                   height,depth,pppt,press,pupt,windu,pvpt,windv,
     &                   pspt,tsurf,ptpt,tempk,conc)
c 
c-----CAMx v4.02 030709
c  
c     UPDTMET updates the time-varying vertical layer structure for the
c     current time step and current grid, based upon linear interpolation
c     between last input time and next input time.  Also performs
c     a similar interpolation for meteorological fields, and rescales
c     the boundary concentrations for the changing met.
c 
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c 
c     Modifications:  
c        none  
c  
c     Input arguments:  
c        igrd                grid index 
c        ncol                number of columns
c        nrow                number of rows
c        nlay                number of layers
c        ngas                number of gas species
c        densfac             density conversion factor (mol/m3)
c        deltat              timestep (s)
c        phpt                time-rate change of layer interface height (m/s)
c        height              layer interface height (m)
c        pppt                time-rate change of pressure (mb/s)
c        press               presssure (mb)
c        pupt                time-rate change of u-component wind (m/s2)
c        windu               u-component wind field (m/s)
c        pvpt                time-rate change of v-component wind (m/s2)
c        windv               v-component wind field (m/s)
c        pspt                time-rate change of surface temperature (K/s)
c        tsurf               surface temperature (K)
c        ptpt                time-rate change of 3-D temperature (K/s)
c        tempk               3-D temperature (K)
c        conc                concentration arrray (umol/m3)
c             
c     Output arguments:  
c        height              layer interface height (m)
c        depth               layer depth (m) 
c        press               presssure (mb)
c        windu               u-component wind field (m/s)
c        windv               v-component wind field (m/s)
c        tsurf               surface temperature (K)
c        tempk               3-D temperature (K)
c        conc                concentration arrray (umol/m3)
c             
c     Routines Called:  
c        none 
c             
c     Called by:  
c        EMISTRNS  
c 
      include 'camx.prm'
      include 'bndary.com'
c
      real height(ncol,nrow,nlay), phpt(ncol,nrow,nlay),
     &     depth(ncol,nrow,nlay), press(ncol,nrow,nlay), 
     &     pppt(ncol,nrow,nlay), windu(ncol,nrow,nlay), 
     &     pupt(ncol,nrow,nlay), windv(ncol,nrow,nlay),
     &     pvpt(ncol,nrow,nlay), tempk(ncol,nrow,nlay),
     &     ptpt(ncol,nrow,nlay), tsurf(ncol,nrow), pspt(ncol,nrow)
      real conc(ncol,nrow,nlay,ngas)
c
c-----Entry point 
c
c-----Convert coarse grid boundary concs from umol/m3 to ppm
c
      if (igrd.eq.1) then
        do 30 l = 1,ngas
          do 20 k = 1,nlay
            do 10 j = 1,nrow
              if (ibeg(j).eq.-999) goto 10 
              i = ibeg(j) - 1 
              convfac = densfac*273./tempk(i,j,k)*press(i,j,k)/1013.
              conc(i,j,k,l) = conc(i,j,k,l)/convfac 
              i = iend(j) + 1 
              convfac = densfac*273./tempk(i,j,k)*press(i,j,k)/1013.
              conc(i,j,k,l) = conc(i,j,k,l)/convfac 
 10         continue
            do 15 i = 1,ncol
              if (jbeg(i).eq.-999) goto 15
              j = jbeg(i) - 1
              convfac = densfac*273./tempk(i,j,k)*press(i,j,k)/1013.
              conc(i,j,k,l) = conc(i,j,k,l)/convfac 
              j = jend(i) + 1
              convfac = densfac*273./tempk(i,j,k)*press(i,j,k)/1013.
              conc(i,j,k,l) = conc(i,j,k,l)/convfac 
 15         continue
 20       continue
 30     continue
      endif
c
c-----Update vertical layer structure forward 1 timestep
c
      do 60 k = 1,nlay
        do 50 j = 1,nrow
          do 40 i = 1,ncol
            height(i,j,k) = height(i,j,k) + deltat*phpt(i,j,k)
            depth(i,j,k) = height(i,j,k)
            if (k.gt.1) depth(i,j,k) = height(i,j,k) - height(i,j,k-1)
c
c-----Update other parameters forward 1 timestep
c
            press(i,j,k) = press(i,j,k) + deltat*pppt(i,j,k)
            if (i.lt.ncol .and. j.lt.nrow) then
              windu(i,j,k) = windu(i,j,k) + deltat*pupt(i,j,k)
              windv(i,j,k) = windv(i,j,k) + deltat*pvpt(i,j,k)
            endif
            tempk(i,j,k) = tempk(i,j,k) + deltat*ptpt(i,j,k)
            if (k.eq.1) tsurf(i,j) = tsurf(i,j) + deltat*pspt(i,j)
c
 40       continue
 50     continue
 60   continue
c
c-----Convert coarse grid boundary concs from ppm to umol/m3
c 
      if (igrd.eq.1) then 
        do 90 l = 1,ngas 
          do 80 k = 1,nlay 
            do 70 j = 1,nrow 
              if (ibeg(j).eq.-999) goto 70 
              i = ibeg(j) - 1 
              convfac = densfac*273./tempk(i,j,k)*press(i,j,k)/1013. 
              conc(i,j,k,l) = conc(i,j,k,l)*convfac 
              i = iend(j) + 1 
              convfac = densfac*273./tempk(i,j,k)*press(i,j,k)/1013. 
              conc(i,j,k,l) = conc(i,j,k,l)*convfac 
 70         continue 
            do 75 i = 1,ncol 
              if (jbeg(i).eq.-999) goto 75 
              j = jbeg(i) - 1 
              convfac = densfac*273./tempk(i,j,k)*press(i,j,k)/1013. 
              conc(i,j,k,l) = conc(i,j,k,l)*convfac 
              j = jend(i) + 1 
              convfac = densfac*273./tempk(i,j,k)*press(i,j,k)/1013. 
              conc(i,j,k,l) = conc(i,j,k,l)*convfac 
 75         continue 
 80       continue 
 90     continue 
      endif
c
      return
      end
