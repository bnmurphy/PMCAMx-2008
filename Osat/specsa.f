c**** SPECSA
c
      subroutine specsa(idate,begtim,jdate,endtim)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine sets up the species names and pointers into the species
c   for all of the tracer species.  Pointers will be set up for both the
c   concentration array and the emissions array.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Argument descriptions:
c      Inputs:
c        idate  I   date of the beginning of the simulation (YYJJJ)
c        begtim R   hour of the begining of simulation
c        jdate  I   date of the ending of the simulation (YYJJJ)
c        endtim R   hour of the endng of simulation
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     05/26/96   --gwilson--    Original development
c     12/12/97   --gwilson--    Fixed bug in initializing the timing
c                               tracers
c     11/06/01   --cemery--     Input dates are now Julian
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'tracer.com'
      include 'filunit.com'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer   idate
      integer   jdate
      real      begtim
      real      endtim     
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*10 name 
      integer      ibegdt, ienddt
      integer      ncount, ioff, idtnow, nhours
      integer      i, j, l
      real         timnow, btim, etim
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      ibegdt = idate
      btim = begtim/100.
      ienddt = jdate
      etim = endtim/100.
c
c  --- calculate the beginning of the various tracer types ---
c      there will be (ngroup+1) if there is an extra group for the 
c      "leftover" group  ----
c
      if( lbndry ) then
          nbdic = 6 
      else
          nbdic = 2 
      endif
      if( ngroup .EQ. 0 ) then
          ncount = nregin
      else
          if( leftovr ) then
             ncount = (ngroup + 1) * nregin
          else
             ncount = ngroup * nregin
          endif
      endif
c
      iptnox = 1
      iemnox = iptnox + nbdic
      iptvoc = iptnox + ncount + nbdic
      iemvoc = iptvoc + nbdic
      ipto3n = iptvoc + ncount + nbdic
      ipto3v = ipto3n + ncount + nbdic
      ipttim = ipto3v + ncount + nbdic
      iemtim = ipttim
c
c   --- calculate the number of tracers at first time step 
c       (1 timing release is updated later) ---
c
      nreles = 0
      ntotsp = ipttim-1
      nsaspc = ipttim-1
      if( ntrtim .EQ. 0 ) npttim = 0
c
c   --- set the names for the initial condition tracers ---
c
      if( ipto3v .GT. MXTRSP ) goto 7000
      ptname(iptnox) = 'NOX000IC  '
      ptname(iptvoc) = 'VOC000IC  '
      ptname(ipto3n) = 'O3N000IC  '
      ptname(ipto3v) = 'O3V000IC  '
c
c   --- if stratifying by boundary there will be 5 boundary condition
c       tracers, otherwise there will be only one ---
c
      if( lbndry ) then
          if( ipto3v+IDXBTP .GT. MXTRSP ) goto 7000
          ptname(iptnox + IDXBNT) = 'NOXNTHBC  '
          ptname(iptvoc + IDXBNT) = 'VOCNTHBC  '
          ptname(ipto3n + IDXBNT) = 'O3NNTHBC  '
          ptname(ipto3v + IDXBNT) = 'O3VNTHBC  '
c
          ptname(iptnox + IDXBES) = 'NOXESTBC  '
          ptname(iptvoc + IDXBES) = 'VOCESTBC  '
          ptname(ipto3n + IDXBES) = 'O3NESTBC  '
          ptname(ipto3v + IDXBES) = 'O3VESTBC  '
c
          ptname(iptnox + IDXBST) = 'NOXSTHBC  '
          ptname(iptvoc + IDXBST) = 'VOCSTHBC  '
          ptname(ipto3n + IDXBST) = 'O3NSTHBC  '
          ptname(ipto3v + IDXBST) = 'O3VSTHBC  '
c
          ptname(iptnox + IDXBWS) = 'NOXWSTBC  '
          ptname(iptvoc + IDXBWS) = 'VOCWSTBC  '
          ptname(ipto3n + IDXBWS) = 'O3NWSTBC  '
          ptname(ipto3v + IDXBWS) = 'O3VWSTBC  '
c
          ptname(iptnox + IDXBTP) = 'NOXTOPBC  '
          ptname(iptvoc + IDXBTP) = 'VOCTOPBC  '
          ptname(ipto3n + IDXBTP) = 'O3NTOPBC  '
          ptname(ipto3v + IDXBTP) = 'O3VTOPBC  '
      else
          if( ipto3v+1 .GT. MXTRSP ) goto 7000
          ptname(iptnox+1) = 'NOX000BC  '
          ptname(iptvoc+1) = 'VOC000BC  '
          ptname(ipto3n+1) = 'O3N000BC  '
          ptname(ipto3v+1) = 'O3V000BC  '
      endif
c
c  --- construct the tracer names and put into names array ---
c
      if( ngroup .EQ. 0 ) then
          ioff = nbdic
          do 10 i=1,nregin
             if( ipto3v+ioff .GT. MXTRSP ) goto 7000
             write(name,'(A,I3.3,I3.3)') 'NOX',1,i
             ptname(iptnox+ioff) = name
             write(name,'(A,I3.3,I3.3)') 'VOC',1,i
             ptname(iptvoc+ioff) = name
             write(name,'(A,I3.3,I3.3)') 'O3N',1,i
             ptname(ipto3n+ioff) = name
             write(name,'(A,I3.3,I3.3)') 'O3V',1,i
             ptname(ipto3v+ioff) = name
             ioff = ioff + 1
   10     continue
      else
          ioff = nbdic 
          if( leftovr ) then
             ncount = ngroup + 1
          else
             ncount = ngroup 
          endif
          do 20 j=1,ncount
             do 30 i=1,nregin
                if( ipto3v+ioff .GT. MXTRSP ) goto 7000
                write(name,'(A,I3.3,I3.3)') 'NOX',j,i
                ptname(iptnox+ioff) = name
                write(name,'(A,I3.3,I3.3)') 'VOC',j,i
                ptname(iptvoc+ioff) = name
                write(name,'(A,I3.3,I3.3)') 'O3N',j,i
                ptname(ipto3n+ioff) = name
                write(name,'(A,I3.3,I3.3)') 'O3V',j,i
                ptname(ipto3v+ioff) = name
                ioff = ioff + 1
   30        continue
   20     continue
      endif
c
c  --- calculate the number of timing tracers there will be and put
c      the names into the names array ---
c
      if( ntrtim .GT. 0 ) then
        if( etim .EQ. 0. ) then
            etim = 24.
            ienddt = ienddt - 1
        endif 
        timnow = btim
        idtnow = ibegdt
        nhours = (ienddt-ibegdt)*24 + INT( etim - btim ) 
        npttim = 1
        do 40 i=1,nhours
           if( MOD( INT(timnow), 24/ntrtim ) .EQ. 0 .OR. i .EQ. 1) then
              do 50 j=1,nregin
                  write(name,'(A,I3.3,I2.2,I3.3)') 'I',MOD(idtnow,1000),
     &                                                  INT(timnow),j
                  if( ntotsp+1 .GT. MXTRSP ) goto 7000
                  ptname(ntotsp+1) = name
                  write(name,'(A,I3.3,I2.2,I3.3)') 'D',MOD(idtnow,1000),
     &                                                  INT(timnow),j
                  if( ntotsp+2 .GT. MXTRSP ) goto 7000
                  ptname(ntotsp+2) = name
                  npttim = npttim + 2
                  ntotsp = ntotsp + 2
   50         continue
           endif
           timnow = timnow + 1.0
           if( timnow .EQ. 24.0 ) then
               timnow = 0.
               idtnow = idtnow + 1
           endif
   40   continue
      endif
c
c   --- check for array overflow ---
c
      if( ntotsp .GT. MXTRSP ) goto 7000
c
c  --- initialize all of the tracers concs to zero to start off ---
c
      do 60 i=1,MXSA3D
        ptconc(i) = 0.
   60 continue
      do 80 l=1,MXTRSP
         ptloft(l) = 0.
         lsamap(l) = l
         do 90 i=1,MXRECP
            conrcp(l,i) = 0.        
   90    continue
   80 continue
      do 11 i=1,MXRECP
         rcpnox(i) = 0.        
         rcpvoc(i) = 0.        
         rcpo3(i) = 0.        
   11 continue
c
c  --- set the flag for outputting the species to average file
c      to true automatically for tracer species ---
c
      do i=1,MXTRSP
         loutsa(i) = .TRUE.
      enddo
c      
c  --- return to calling routine ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue 
      write(iout,'(//,a)') 'ERROR in SPECSA:'
      write(iout,'(/,1X,A,I10)') 
     &          'Number of tracer species exceeds max: ',MXTRSP
      write(iout,'(1X,A)') 'Increase parameter: MXTRSP'
      call camxerr()
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
