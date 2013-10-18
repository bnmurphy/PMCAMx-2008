      subroutine rdpthdr(xloc,yloc,begtim,begdate,endtim,enddate)
c 
c-----CAMx v4.02 030709
c 
c     RDPTHDR reads the header of binary point source emissions file,
c     initializes time-invariant point source variables, and maps the point
c     source species list to the internal CAMx species list
c                           
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c           
c     Modifications: 
c        10/24/01  Removed BSWAP and converted integer strings to character*4
c        11/06/01  Input dates are now Julian
c        02/09/02  Added code to handle end of year dates
c 
c     Input arguments: 
c        begtim              model start time (HHMM) 
c        begdate             model start date (YYJJJ) 
c        endtim              model end time (HHMM) 
c        enddate             model end date (YYJJJ)
c             
c     Output arguments: 
c        xloc                stack x-coordinates
c        yloc                stack y-coordinates
c             
c     Routines Called: 
c        none
c             
c     Called by: 
c        PNTPREP
c 
      include 'camx.prm'
      include 'grid.com'
      include 'chmstry.com'
      include 'ptemiss.com'
      include 'filunit.com'
      include 'flags.com'
      include 'bndary.com'
c
c======================== Source Apportion Begin =======================
c
      include 'tracer.com'
c
c========================= Source Apportion End ========================
c

      character*4 ifile(10),note(60),ptspec(10,MXSPEC) 
      integer begdate,enddate
      dimension xloc(MXPTSRC),yloc(MXPTSRC)
      character*10 ptfil,infil,ptspc 
c
      data ptfil /'PTSOURCE  '/
c
c-----Entry point
c             
c-----Read 1st PT header record and check inputs 
c             
      rewind(iptem)
      read(iptem) ifile,note,nseg,nptspc,idat1,tim1,idat2,tim2
      if(INT(tim2) .EQ. 24) then
        tim2 = 0.
        idat2 = idat2 + 1
        if( MOD(idat2,1000) .GT. 365 ) then
           if( MOD(INT(idat2/1000),4) .EQ. 0 ) then
              if( MOD(idat2,1000) .EQ. 367 ) 
     &                       idat2 = (INT(idat2/1000)+1)*1000 + 1
           else
              idat2 = (INT(idat2/1000)+1)*1000 + 1
           endif
        endif
      endif
      write(infil,'(10a1)') (ifile(n),n=1,10) 
      if (infil.ne.ptfil) then 
        write(iout,'(//,a)') 'ERROR in RDPTHDR:'
        write(iout,*)'PT input file is not labelled PTSOURCE' 
        call camxerr()
      endif   
      if (nptspc.gt.MXSPEC) then 
        write(iout,'(//,a)') 'ERROR in RDPTHDR:'
        write(iout,*)'Number of species on PT file > MXSPEC' 
        write(iout,*) nptspc,MXSPEC 
        call camxerr()
      endif
      tim1 = 100.*tim1
      tim2 = 100.*tim2
      if (idat1.gt.begdate) then 
        write(iout,'(//,a)') 'WARNING in RDPTHDR:'
        write(iout,*)'PT start date > simulation start date' 
        write(iout,'(a,i10.5,a,i10.5)')
     &        ' PT file: ',idat1,' Sim start: ',begdate 
        if (.not.le1day) then
          write(iout,*) 'Day-specific emissions supplied: CAMx Stopping'
          call camxerr()
        endif
      elseif (idat1.eq.begdate .and. tim1.gt.begtim) then 
        write(iout,'(//,a)') 'WARNING in RDPTHDR:'
        write(iout,*)'PT start time > simulation start time' 
        write(iout,*)'PT file: ',tim1,' Sim start: ',begtim 
        if (.not.le1day) then
          write(iout,*) 'Day-specific emissions supplied: CAMx Stopping'
          call camxerr()
        endif
      elseif (idat2.lt.enddate) then 
        write(iout,'(//,a)') 'WARNING in RDPTHDR:'
        write(iout,*)'PT end date < simulation end date' 
        write(iout,'(a,i10.5,a,i10.5)')
     &        'PT file: ',idat2,' Sim end: ',enddate 
        if (.not.le1day) then
          write(iout,*) 'Day-specific emissions supplied: CAMx Stopping'
          call camxerr()
        endif
      elseif (idat2.eq.enddate .and. tim2.lt.endtim) then 
        write(iout,'(//,a)') 'WARNING in RDPTHDR:'
        write(iout,*)'PT end time < simulation end time' 
        write(iout,*)'PT file: ',tim2,' Sim end: ',endtim 
        if (.not.le1day) then
          write(iout,*) 'Day-specific emissions supplied: CAMx Stopping'
          call camxerr()
        endif
      endif 
c 
c-----Read 2nd PT header record and check inputs 
c 
      read(iptem) orgx,orgy,izone,utmx,utmy,dx,dy,nx,ny,nz 
      if (nx.ne.ncol(1) .or. ny.ne.nrow(1)) then 
        write(iout,'(//,a)') 'WARNING in RDPTHDR:'
        write(iout,*)'PT grid size not equal to model grid size' 
        write(iout,*)'PT file: ',nx,ny,nz, 
     &               ' model: ',ncol(1),nrow(1),nlay(1) 
      endif 
c 
c-----Read 3rd & 4th PT header 
c 
      read(iptem) (idum,idum,idum,idum,n=1,nseg) 
      read(iptem) ((ptspec(n,l),n=1,10),l=1,nptspc) 
c 
c-----Map PT species to model species 
c 
      do 15 lpt = 1,nptspc 
        lptmap(lpt) = 0
        write(ptspc,'(10a1)') (ptspec(n,lpt),n=1,10) 
        if (ptspc.eq.'HNO2      ') ptspc = 'HONO      '
        if (ptspc.eq.'HCHO      ' .and. kHCHO.eq.nspec+1)
     &                                        ptspc = 'FORM      '
        do 20 l = 1,nspec 
          if (ptspc.eq.spname(l)) then 
            lptmap(lpt) = l
            write(idiag,'(2(a,i5,2x,a))')
     &                   'Point source species ',lpt,ptspc, 
     &                   ' mapped to model species ',l,spname(l) 
            goto 15 
          endif 
 20     continue 
        write(idiag,*)'PT species: ',ptspc,' not modeled'
 15   continue 
c
c-----Read time invariant data (assume nseg=1)
c     check number of sources against max 
c
      read(iptem) idum,nptsrc
      if (nptsrc.gt.MXPTSRC) then
        write(iout,'(//,a)') 'ERROR in RDPTHDR:'
        write(iout,*)'Number of point sources > MXPTSRC'
        write(iout,*) nptsrc,MXPTSRC
        call camxerr()
      endif
      read(iptem) (xloc(n),yloc(n),hstk(n),dstk(n),tstk(n),vstk(n),
     &             n=1,nptsrc)
c
      do 30 n = 1,nptsrc
        if (dstk(n).lt.0.) then
          lpiglet(n) = .true.
c         dstk(n) = -dstk(n)
        else
          lpiglet(n) = .false.
        endif
        vstk(n) = vstk(n)/3600.
c
c======================== Source Apportion Begin =======================
c
        xlocpt(n) = xloc(n)
        ylocpt(n) = yloc(n)
c
c========================= Source Apportion End ========================
c
  30  continue
      write(idiag,*)
c
      return
      end
