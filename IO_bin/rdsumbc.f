      subroutine rdsumbc(idtnow,timnow,jdlast,ttlast,ncols,nrows,nlays,
     &                   contop,consum)
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c
c-----CAMx v4.03 031205
c
c     11/20/03   --gwilson--    Fixed bug when using the species map
c                               to get the index of modeled species
c
      include 'camx.prm'
      include 'camx.com'
      include 'filunit.com'
      include 'chmstry.com'
      include 'bndary.com'
      include 'tracer.com'
c
      character*4 iname(10)
      real      height(MXCOLA*MXROWA*MXLAYA),depth(MXCOLA*MXROWA*MXLAYA)
      real      consum(MXSPEC,0:IDXBTP), contop(MXSPEC)
      real      conwst(MXLAYA,MXROWA), conest(MXLAYA,MXROWA)
      real      consth(MXLAYA,MXCOLA), connth(MXLAYA,MXCOLA)
      logical   lfound
c
      ibgdhp = 0
      btimhp = 0 
      ibegbc = 0
      btimbc = 0.
      iendbc = 0
      etimbc = 0.
c
c   --- read the ZP file to get the layer heights ----
c
      rewind(ihtp(1))
  333 continue
c
      call getdepth(ncols,nrows,nlays,ly2k,ibgdhp,idtnow,btimhp,timnow,
     &              ihtp(1),height,depth)
c
c   --- initialize the array to zero ---
c
      do 60 i=1,MXSPEC
        consum(i,0) = 0.
        consum(i,1) = 0.
        consum(i,2) = 0.
        consum(i,3) = 0.
        consum(i,4) = 0.
        consum(i,5) = 0.
   60 continue
c
c   --- read boundary conditions for this hour ----
c
      read(ibc,ERR=7002) ibegbc, btimbc, iendbc, etimbc
      lfound = .TRUE.
      if( ibegbc .EQ. idtnow .AND. btimbc .GT. timnow ) lfound = .FALSE.
      if( ibegbc .GT. idtnow ) lfound = .FALSE.
      if( etimbc .EQ. 0. ) then
         iendbc = iendbc - 1
         etimbc = 24.0
      endif
      do 70 ispc=1,nbcspc
          read(ibc,ERR=7002,END=7003) iseg, (iname(i),i=1,10), iedge,
     &                         ((conwst(izcl,j),izcl=1,nlays),j=1,nrows)
          read(ibc,ERR=7002,END=7003) iseg, (iname(i),i=1,10), iedge,
     &                         ((conest(izcl,j),izcl=1,nlays),j=1,nrows)
          read(ibc,ERR=7002,END=7003) iseg, (iname(i),i=1,10), iedge,
     &                         ((consth(izcl,i),izcl=1,nlays),i=1,ncols)
          read(ibc,ERR=7002,END=7003) iseg, (iname(i),i=1,10), iedge,
     &                         ((connth(izcl,i),izcl=1,nlays),i=1,ncols)
c
c   --- if record does not span this hour, skip it ----
c
          if( .NOT. lfound ) goto 70
c
c   --- if the species is a not modelled or not a VOC species skip it ---
c
          idx = 0
          do i=1,nspec
            if( lbcmap(i) .EQ. ispc ) idx = i
          enddo
          if( idx .LE. 0 ) goto 70
          if( .NOT. lvocsp(idx) ) goto 70
          call sumwt4(idx,ncols,nrows,nlays,depth,conwst,conest,
     &                consth,connth,consum)
c
c   --- if record does not span this hour, skip it ----
c
c   --- next species ---
c
   70 continue
c
c   --- add the top concentrations for this hour ---
c
      do 51 ispc=1,nspec
         if( .NOT. lvocsp(ispc) ) goto 51
         do 61 j=1,nrows
            if( ibeg(j) .GT. 0 ) then
                do 71 i=ibeg(j),iend(j)
                    consum(ispc,IDXBTP) = consum(ispc,IDXBTP) + 
     &                                                   contop(ispc)
                    consum(ispc,0) = consum(ispc,0) + contop(ispc)
   71            continue
             endif
   61    continue
c
c   --- get next species in the TOPCONC file ---
c
   51 continue
c
c   --- chack date and time, if it is still in the episode
c       go back to read the next hour ----
c
      if( lfound ) then
         timnow = timnow + 1.0
         if( timnow .GT. 24.0 ) then
             idtnow = idtnow + 1
             timnow = 0.0
             if( MOD(idtnow,1000) .GT. 365 ) then
                if( MOD(INT(idtnow/1000),4) .EQ. 0 ) then
                   if( MOD(idtnow,1000) .EQ. 367 )
     &                     idtnow = (INT(idtnow/1000)+1)*1000 + 1
                else
                   idtnow = (INT(idtnow/1000)+1)*1000 + 1
                endif
             endif
         endif
      endif
      if( idtnow .GT. jdlast .OR. 
     &          (idtnow .EQ. jdlast .AND. timnow .GT. ttlast) ) goto 444
c
c  --- check if we need to back up the boundary conditions file, that is
c      does the record just read contain this hour ----
c
      if( idtnow .LT. iendbc .OR. 
     &           (idtnow .EQ. iendbc .AND.  timnow .LE. etimbc ) ) then
         do 81 i=1,nbcspc*4
            backspace(ibc)
   81    continue
         backspace(ibc)
      endif
c
c  --- process next hour ----
c
      goto 333
c
c  --- all concentrations are summed, calculate the weghted fraction ----
c
  444 continue
c
c   ---- rewind files to be used again by the regular model ---
c
      rewind(ibc)
      rewind(itopc)
      rewind(ihtp(1))
      return
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in RDSUMBC:'
      write(iout,'(/,1X,A)') 'Reading boundary conditions. '
      call camxerr()
c
 7003 continue
      write(iout,'(//,a)') 'ERROR in RDSUMBC:'
      write(iout,'(/,1X,2A)') 'Premature end-of-file reading ',
     &                                         'boundary conditions.'
      call camxerr()
      end
