c*** AVGRCPDDM
c
      subroutine avgrcpddm(idate,etime,nox,noy,nspc,conc)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c   Description:
c     This routine writes the average sensitivities to the 
c     DDM receptor file.  
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
c        nspc    I  number of regular model species
c        conc    R  regualar model concentrations
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c    07/18/01  --gyarwood--  Original development from AVGRCP for OSAT
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
      integer   nspc
      real      conc(nox,noy,nspc)
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer      idatbg, idatnd, ircp, i, l
      real         timbg, timnd
      real         modnox, modvoc, modo3
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c --- set the date and time ---
c
      idatnd = idate
      timnd = ANINT(etime)/100.
      if( timnd .EQ. 0. ) then
         timnd = 24.0
         idatnd = idatnd - 1
      endif
      idatbg = idatnd
      if( dtout .GE. 60. ) then
          timbg = AINT( 1000*(timnd - ANINT(dtout)/60.) )/1000.
      else
          timbg = AINT( 1000*(timnd - ANINT(dtout)/100.) )/1000.
      endif
      if(timbg .LT. -0.01 ) then
        timbg = timbg + 24.
        idatbg = idatbg - 1
      endif
c
c ---  write the header for this hour ---
c
      write(IOWAVG,9000,ERR=7000)
      write(IOWAVG,9000,ERR=7000) 'Data for Period',
     &                                    idatbg, timbg, idatnd, timnd
c
c --- output the hourly peak receptor, only if peak
c     cell has been found ---
c
      write(IOWAVG,9001,ERR=7000) 1
      if( ipekcl(1) .GT. 0 .AND. ipekcl(2) .GT. 0 ) then
        do i = 1, ntotsp
          if(loutsa(i)) write(IOWAVG,9002,ERR=7000) conrcp(i,1)
        enddo
      else
        do i = 1, ntotsp
          if(loutsa(i)) write(IOWAVG,9002,ERR=7000) 
     &                        -999.
        enddo
      endif
      modnox = 0.
      modvoc = 0.
      modo3 = 0.
      if( ipekcl(1) .GT. 0 .AND. ipekcl(2) .GT. 0 ) then
        do i=1,nspc
          if( lnoxsp(i) ) then
            modnox = modnox + conc(ipekcl(1),ipekcl(2),i)
          else if( lvocsp(i) ) then
            modvoc = modvoc + conc(ipekcl(1),ipekcl(2),i) * crbnum(i)
          else if( lo3sp(i) ) then
            modo3 = modo3 + conc(ipekcl(1),ipekcl(2),i)
          endif
        enddo
      endif
      write(IOWAVG,9002,ERR=7000) modnox, modvoc, modo3
      write(IOWAVG,*,ERR=7000)
c
c --- loop over all receptors, except the hourly peak (receptor 1)  ---
c
      do ircp=2,nrecep
        write(IOWAVG,9001,ERR=7000) ircp
        do i = 1, ntotsp
          if(loutsa(i)) write(IOWAVG,9002,ERR=7000) 
     &                        conrcp(i,ircp)
        enddo
        write(IOWAVG,9002,ERR=7000) rcpnox(ircp), rcpvoc(ircp), 
     &                                                   rcpo3(ircp)
        write(IOWAVG,*,ERR=7000)
      enddo
      call flush(IOWAVG)
c
c --- re-initialize the averages to zero ---
c
      do l=1,MXTRSP
        do i=1,MXRECP
           conrcp(l,i) = 0.        
        enddo
      enddo
      do i=1,MXRECP
        rcpnox(i) = 0.        
        rcpvoc(i) = 0.        
        rcpo3(i) = 0.        
      enddo
c
c --- return to calling routine ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in AVGRCPDDM:'
      write(iout,9010,ERR=9999) 'Writing DDM average file at: ',
     &                           idatbg,timbg,idatnd,timnd
      call camxerr()
c
c-----------------------------------------------------------------------
c    Format statements:
c-----------------------------------------------------------------------
c
 9000 format(A,',',2(I10,',',F10.2,','))
 9001 format(I5,',',$)
 9002 format(:,1p1e11.4,',',$)
 9010 format(/,1X,A,',',2(I10,',',F10.2,','))
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
