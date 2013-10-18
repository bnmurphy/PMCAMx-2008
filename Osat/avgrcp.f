c*** AVGRCP
c
      subroutine avgrcp(idate,etime,nox,noy,nspsa,nspc,saconc,conc)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c   Description:
c     This routine writes the average concentrations to the 
c     receptor average file.  
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c   Argument descriptions:
c     Inputs:
c        idate   I  ending date of this time interval
c        etime   R  ending time of this time interval
c        nox     I  number of cells in the domain
c        noy     I  number of cells in the domain
c        nspsa   I  number of regular model species
c        nspc    I  number of regular model species
c        saconc  R  tracer concentrations
c        conc    R  regualar model concentrations
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c    12/11/2000 --gwilson-- Now handles the case where peak cell could
c                           not be found.  This could happen if not
c                           outputting ozone as an average species.
c     04/02/03  --gwilson-- Added grid number to recptors defined by
c                           cell index
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'camx.com'
      include 'tracer.com'
      include 'filunit.com'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer   idate 
      real      etime 
      integer   nox
      integer   noy
      integer   nspsa 
      integer   nspc
      real      saconc(nox,noy,nspsa)
      real      conc(nox,noy,nspc)
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer   idatbg, idatnd, ircp, i, l
      real      timbg, timnd, sumnox, sumvoc, sumo3
      real      modnox, modvoc, modo3
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  ---- set the date and time ---
c
      idatnd = idate
      timnd = ANINT(etime)/100.
      if( timnd .EQ. 0. ) then
         timnd = 24.0
         idatnd = idatnd - 1
      endif
      idatbg = idatnd
      if( dtout .GE. 60.0 ) then
          timbg = ANINT( 1000*(timnd - ANINT(dtout)/60.) )/1000.
      else
          timbg = ANINT( 1000*(timnd - ANINT(dtout)/100.) )/1000.
      endif
      if(timbg .LT. 0. ) then
        timbg = timbg + 24.
        idatbg = idatbg - 1
      endif
c
c  ---  write the header for this hour ---
c
      write(IOWAVG,9000,ERR=7000)
      write(IOWAVG,9000,ERR=7000) 'Data for Period',
     &                                    idatbg, timbg, idatnd, timnd
c
c  --- output the hourly peak receptor, only if peak
c      cell has been found ---
c
      if( ipekcl(1) .GT. 0 .AND. ipekcl(2) .GT. 0 ) then
        write(IOWAVG,9001,ERR=7000) 'Receptor',1
        sumnox = 0.
        do 10 i=iptnox,iptvoc-1
            sumnox = sumnox + saconc(ipekcl(1),ipekcl(2),i)
   10   continue
        write(IOWAVG,9002,ERR=7000) (saconc(ipekcl(1),ipekcl(2),i),
     &                                            i=iptnox,iptvoc-1)
        sumvoc = 0.
        do 20 i=iptvoc,ipto3n-1
            sumvoc = sumvoc + saconc(ipekcl(1),ipekcl(2),i)
   20   continue
        write(IOWAVG,9002,ERR=7000) (saconc(ipekcl(1),ipekcl(2),i),
     &                                            i=iptvoc,ipto3n-1)
        sumo3 = 0.
        do 30 i=ipto3n,ipttim-1
            sumo3 = sumo3 + saconc(ipekcl(1),ipekcl(2),i)
   30   continue
        write(IOWAVG,9002,ERR=7000) (saconc(ipekcl(1),ipekcl(2),i),
     &                                               i=ipto3n,ipto3v-1)
        write(IOWAVG,9002,ERR=7000) (saconc(ipekcl(1),ipekcl(2),i),
     &                                               i=ipto3v,ipttim-1)
      endif
c
      if( ntrtim .GT. 0 ) then
         write(IOWAVG,9002,ERR=7000) (saconc(ipekcl(1),ipekcl(2),i),
     &               i=ipttim,nsaspc-1,2),(BNDLPT,i=nsaspc+1,npttim-1,2)
         write(IOWAVG,9002,ERR=7000) (saconc(ipekcl(1),ipekcl(2),i),
     &                 i=ipttim+1,nsaspc,2),(BNDLPT,i=nsaspc+2,npttim,2)
      endif
      write(IOWAVG,9002,ERR=7000) sumnox, sumvoc, sumo3
      modnox = 0.
      modvoc = 0.
      modo3 = 0.
      do 40 i=1,nspc
         if( lnoxsp(i) ) then
            modnox = modnox + conc(ipekcl(1),ipekcl(2),i)
         else if( lvocsp(i) ) then
            modvoc = modvoc + conc(ipekcl(1),ipekcl(2),i) * crbnum(i)
         else if( lo3sp(i) ) then
            modo3 = modo3 + conc(ipekcl(1),ipekcl(2),i)
         endif
   40 continue
      write(IOWAVG,9002,ERR=7000) modnox, modvoc, modo3
c
c  --- loop over all receptors, except the hourly peak (receptor 1)  ---
c
      do 50 ircp=2,nrecep
c
c  --- write the header for this receptor ----
c
          write(IOWAVG,9001,ERR=7000) 'Receptor',ircp
          sumnox = 0.
          do 60 i=iptnox,iptvoc-1
              sumnox = sumnox + conrcp(i,ircp)
   60     continue
          write(IOWAVG,9002,ERR=7000) (conrcp(i,ircp),i=iptnox,iptvoc-1)
          sumvoc = 0.
          do 70 i=iptvoc,ipto3n-1
              sumvoc = sumvoc + conrcp(i,ircp)
   70     continue
          write(IOWAVG,9002,ERR=7000) (conrcp(i,ircp),i=iptvoc,ipto3n-1)
          sumo3 = 0.
          do 80 i=ipto3n,ipttim-1
              sumo3 = sumo3 + conrcp(i,ircp)
   80     continue
          write(IOWAVG,9002,ERR=7000) (conrcp(i,ircp),i=ipto3n,ipto3v-1)
          write(IOWAVG,9002,ERR=7000) (conrcp(i,ircp),i=ipto3v,ipttim-1)
c
          if( ntrtim .GT. 0 ) then
              write(IOWAVG,9002,ERR=7000) (conrcp(i,ircp),
     &               i=ipttim,nsaspc-1,2),(BNDLPT,i=nsaspc+1,npttim-1,2)
              write(IOWAVG,9002,ERR=7000) (conrcp(i,ircp),
     &                 i=ipttim+1,nsaspc,2),(BNDLPT,i=nsaspc+2,npttim,2)
          endif
          write(IOWAVG,9002,ERR=7000) sumnox, sumvoc, sumo3
          write(IOWAVG,9002,ERR=7000) rcpnox(ircp), rcpvoc(ircp), 
     &                                                      rcpo3(ircp)
c
c  --- next receptor ----
c
   50 continue
      call flush(IOWAVG)
c
c  ---- re-initialize the averages to zero ---
c
      do 90 l=1,MXTRSP
         do 11 i=1,MXRECP
            conrcp(l,i) = 0.        
   11    continue
   90 continue
      do 21 i=1,MXRECP
         rcpnox(i) = 0.        
         rcpvoc(i) = 0.        
         rcpo3(i) = 0.        
   21 continue
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
      write(iout,'(//,a)') 'ERROR in AVGRCP:'
      write(iout,9010,ERR=9999) 'Writing Tracer average file at: ',
     &                           idatbg,timbg,idatnd,timnd
      call camxerr()
c
c-----------------------------------------------------------------------
c    Format statements:
c-----------------------------------------------------------------------
c
 9000 format(A,',',2(I10.5,',',F10.2,','))
 9001 format(A10,',',I5,',')
 9002 format(500(:,1p1e10.4,',',1X))
 9010 format(/,1X,A,',',2(I10.5,',',F10.2,','))
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
