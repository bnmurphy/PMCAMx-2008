      subroutine bndprep(begtim,begdate,endtim,enddate)
c 
c-----CAMx v4.02 030709
c 
c     BNDPREP reads the boundary conditions file header and sets
c     mapping variables to allow for irregular coarse grid boundaries.
c     It also prepares mapping variables that map the BC species list
c     to the internal CAMx species list
c                           
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c           
c     Modifications: 
c        1/20/99   Grid cell size from file should be meters for all cartesian
c                  projections (UTM, LCP, PSP)
c        10/24/01  Removed BSWAP and converted integer strings to character*4
c        10/31/01  Revised diagnostic checks on ibeg and jbeg for irregular
c                  boundaries
c        11/06/01  Input dates are now Julian
c 
c     Input arguments: 
c        begtim              model start time (HHMM)
c        begdate             model start date (YYJJJ)
c        endtim              model end time (HHMM)
c        enddate             model end date (YYJJJ)
c             
c     Output arguments: 
c        none
c             
c     Routines Called: 
c        none
c             
c     Called by: 
c        STARTUP 
c
      include 'camx.prm'
      include 'filunit.com'
      include 'bndary.com'
      include 'grid.com'
      include 'chmstry.com'
      include 'flags.com'
c
      character*4 ifile(10),note(60),bcspec(10,MXSPEC)
      integer begdate,enddate
      dimension indx(4,MX1D)
      character*10 bcfil,infil,bcspc
c
      data bcfil /'BOUNDARY  '/
c
c-----Entry point
c
c-----Read 1st BC header record and check inputs
c
      rewind(ibc)
      read(ibc) ifile,note,nseg,nbcspc,idat1,tim1,idat2,tim2
      if(INT(tim2) .EQ. 24) then
        idat2 = idat2 + 1
        tim2 = 0.
        if( MOD(idat2,1000) .GT. 365 ) then
           if( MOD(INT(idat2/1000),4) .EQ. 0 ) then
              if( MOD(idat2,1000) .EQ. 367 )
     &                     idat2 = (INT(idat2/1000)+1)*1000 + 1
           else
              idat2 = (INT(idat2/1000)+1)*1000 + 1
           endif
        endif
      endif
      write(infil,'(10a1)') (ifile(n),n=1,10)
      if (infil.ne.bcfil) then
        write(iout,'(//,a)') 'ERROR in BNDPREP:'
        write(iout,*)'BC input file is not labelled BOUNDARY'
        call camxerr()
      endif
      if (nbcspc.gt.MXSPEC) then
        write(iout,'(//,a)') 'ERROR in BNDPREP:'
        write(iout,*)'Number of species on BC file > MXSPEC'
        write(iout,*) nbcspc,MXSPEC
        call camxerr()
      endif
      tim1 = 100.*tim1
      tim2 = 100.*tim2
      if (idat1.gt.begdate) then
        write(iout,'(//,a)') 'ERROR in BNDPREP:'
        write(iout,*)'BC start date > simulation start date'
        write(iout,'(a,i10.5,a,i10.5)')
     &        'BC file: ',idat1,' Sim start: ',begdate
c        call camxerr()
      elseif (idat1.eq.begdate .and. tim1.gt.begtim) then
        write(iout,'(//,a)') 'ERROR in BNDPREP:'
        write(iout,*)'BC start time > simulation start time'
        write(iout,*)'BC file: ',tim1,' Sim start: ',begtim
c        call camxerr()
      elseif (idat2.lt.enddate) then
        write(iout,'(//,a)') 'ERROR in BNDPREP:'
        write(iout,*)'BC end date < simulation end date'
        write(iout,'(a,i10.5,a,i10.5)')
     &        'BC file: ',idat2,' Sim end: ',enddate
c        call camxerr()
      elseif (idat2.eq.enddate .and. tim2.lt.endtim) then
        write(iout,'(//,a)') 'ERROR in BNDPREP:'
        write(iout,*)'BC end time < simulation end time'
        write(iout,*)'BC file: ',tim2,' Sim end: ',endtim
c        call camxerr()
      endif
c
c-----Read 2nd BC header record and check inputs
c
      read(ibc) orgx,orgy,izone,utmx,utmy,dx,dy,nx,ny,nz
      if (.NOT.llatlon) then
        dx = dx/1000.
        dy = dy/1000.
      endif
      if (dx.ne.delx .or. dy.ne.dely) then
        write(iout,'(//,a)') 'WARNING in BNDPREP:'
        write(iout,*)'BC cell size not equal to model cell size'
        write(iout,*)'BC file: ',dx,dy,' model: ',delx,dely
        write(iout,*)
      elseif (nx.ne.ncol(1) .or. ny.ne.nrow(1) .or. nz.ne.nlay(1)) then
        write(iout,'(//,a)') 'ERROR in BNDPREP:'
        write(iout,*)'BC grid size not equal to model grid size'
        write(iout,*)'BC file: ',nx,ny,nz,
     &               ' model: ',ncol(1),nrow(1),nlay(1)
        write(iout,*)
        call camxerr()
      endif 
c
c-----Read 3rd & 4th BC header 
c
      read(ibc) (idum,idum,idum,idum,n=1,nseg)
      read(ibc) ((bcspec(n,l),n=1,10),l=1,nbcspc)
c
c-----Map BC species to model species
c
      do 20 l = 1,nspec
        lbcmap(l) = 0
        do 15 lbc = 1,nbcspc
          write(bcspc,'(10a1)') (bcspec(n,lbc),n=1,10)
          if (bcspc.eq.'HNO2      ') bcspc = 'HONO      '
          if (bcspc.eq.'HCHO      ' .and. kHCHO.eq.nspec+1)
     &                                        bcspc = 'FORM      '
          if (bcspc.eq.spname(l)) then
            lbcmap(l) = lbc
            write(idiag,'(2(a,i5,2x,a))')'Boundary species ',lbc,bcspc,
     &                   ' mapped to model species ',l,spname(l)
            goto 20
          endif
 15     continue
        write(idiag,*)'Did not find species: ',spname(l),' on BC file'
 20   continue
c
c-----Read time invariant data (assume nseg=1)
c     assign coarse grid starting/ending indices for each row/column
c
      do 25 nedge = 1,4
        read(ibc) isegm,iedge,ncells,((indx(n,i),n=1,4),i=1,ncells)
        do 30 i = 1,ncells
          if (nedge.eq.1) then 
            if (indx(1,i).eq.0) then
              ibeg(i) = -999
            else
              ibeg(i) = indx(1,i)
            endif
          elseif (nedge.eq.2) then
            if (indx(1,i).eq.0) then
              iend(i) = -999
            else
              iend(i) = indx(1,i)
            endif
          elseif (nedge.eq.3) then
            if (indx(1,i).eq.0) then
              jbeg(i) = -999
            else
              jbeg(i) = indx(1,i)
            endif
          else
            if (indx(1,i).eq.0) then
              jend(i) = -999
            else
              jend(i) = indx(1,i)
            endif
          endif 
  30    continue
  25  continue
c
c-----Ensure at the very least that rows 1 & ny, and columns 1 & nx
c     are boundary cells
c
      ibeg(1) = -999
      iend(1) = -999
      ibeg(nrow(1)) = -999
      iend(nrow(1)) = -999
      jbeg(1) = -999
      jend(1) = -999
      jbeg(ncol(1)) = -999
      jend(ncol(1)) = -999
      do i = 2,ncol(1)-1
        if (jbeg(i).le.1 .and. jbeg(i).ne.-999) then
          write(idiag,*)'WARNING: First computational cell in col',i,
     &                  ' is',jbeg(i),' setting to',2
          jbeg(i) = 2
        endif
        if (jend(i).ge.nrow(1)) then
          write(idiag,*)'WARNING:  Last computational cell in col',i,
     &                  ' is',jend(i),' setting to',nrow(1)-1
          jend(i) = nrow(1)-1
        endif
      enddo
      do j = 2,nrow(1)-1
        if (ibeg(j).le.1 .and. ibeg(j).ne.-999) then
          write(idiag,*)'WARNING: First computational cell in row',j,
     &                  ' is',ibeg(j),' setting to',2
          ibeg(j) = 2
        endif
        if (iend(j).ge.ncol(1)) then
          write(idiag,*)'WARNING:  Last computational cell in row',j,
     &                  ' is',ibeg(j),' setting to',ncol(1)-1
          iend(j) = ncol(1)-1
        endif
      enddo
c
c-----Echo coarse grid boundary definition, and check that if nested grids
c     are specified, their boundaries are not outside the coarse
c     grid region
c
      write(idiag,*)
      write(idiag,*)'Coarse grid computational domain definition'
      do 50 j = 1,nrow(1)
        write(idiag,*)'row',j,' spans cell',ibeg(j),' through',iend(j)
        do 40 n = 2,ngrid
          if (j.lt.jnst1(n) .or. j.gt.jnst2(n)) goto 40
          if (ibeg(j).eq.-999 .or. inst1(n).lt.ibeg(j) .or.
     &        inst2(n).gt.iend(j)) then
            write(iout,'(//,a)') 'ERROR in BNDPREP:'
            write(iout,*)'Nested grid extends into coarse grid',
     &                   ' boundary region -- Nest = ',n
            write(iout,*)'j,ibeg,iend,inst1,inst2: ',
     &                    j,ibeg(j),iend(j),inst1(n),inst2(n)
            call camxerr()
          endif
 40     continue
 50   continue
      write(idiag,*)
      do 70 i = 1,ncol(1)
        write(idiag,*)'col',i,' spans cell',jbeg(i),' through',jend(i)
        do 60 n = 2,ngrid
          if (i.lt.inst1(n) .or. i.gt.inst2(n)) goto 60
          if (jbeg(i).eq.-999 .or. jnst1(n).lt.jbeg(i) .or.
     &        jnst2(n).gt.jend(i)) then
            write(iout,'(//,a)') 'ERROR in BNDPREP:'
            write(iout,*)'Nested grid extends into coarse grid',
     &                   ' boundary region -- Nest = ',n
            write(iout,*)'i,jbeg,jend,jnst1,jnst2: ',
     &                    i,jbeg(i),jend(i),jnst1(n),jnst2(n)
            call camxerr()
          endif
 60     continue
 70   continue
      write(idiag,*)
c 
      return
      end
