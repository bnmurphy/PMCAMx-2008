c**** SPECDDM
c
      subroutine specddm( )
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine sets up the species names and pointers into the species
c   for all of the DDM species.  Pointers will be set up for both the
c   concentration array and the emissions array.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     03/23/99   --gwilson--    Original development
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'filunit.com'
      include 'chmstry.com'
      include 'tracer.com'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*10 caffec, cinflu, srcnam
      character*3  edgnam(5)
      integer      ispc, iaffec, iddm, nedge, i
      logical      lout 
c
c-----------------------------------------------------------------------
c    Data statements:
c-----------------------------------------------------------------------
c
      data nedge /5/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- load the names of the boundary edges ---
c
      edgnam(IDXBWS) = 'WST'
      edgnam(IDXBES) = 'EST'
      edgnam(IDXBST) = 'STH'
      edgnam(IDXBNT) = 'NTH'
      edgnam(IDXBTP) = 'TOP'
c
c  --- echo the table of species to the output file ---
c
      write(idiag,*) ' '
      write(idiag,*) ' Affected   Influencing   Source',
     &                   '                       Long            Short' 
      write(idiag,*) ' Species      Species      Type ',
     &                   '      Group   Region   Name            Name' 
      write(idiag,*) ('-',i=1,79)
c
c  --- calculate the number of DDM species per family ---
c
      if( lbndry ) then
          nbdic = 5 * nbcddm
      else
          nbdic = nbcddm
      endif
      if( ngroup .EQ. 0 ) then
          nddmsp = (nicddm + nbdic) + nemddm * nregin
      else
          nddmsp = (nicddm + nbdic) + nemddm * nregin * ngroup
      endif
c
c  --- loop over the modeled species and setup the names ---
c
      do ispc=1,ngas
c
c  --- set the pointer into the array for this family ---
c
          iptddm(ispc) = (ispc-1)*nddmsp + 1
          iaffec = ispc
          caffec = spname(ispc)
c
c  --- set the flag for determining if this species should
c      be output to average file ---
c
          lout = .FALSE.
          do i=1,navspc
            if( lavmap(i) .EQ. ispc ) lout = .TRUE.
          enddo
c
c  --- intialize the index into the array ---
c
          idxddm = iptddm(ispc) - 1
c
c  --- loop over all of the DDM initial condition speicies ---
c
          do iddm = 1,nicddm
             cinflu = icddmsp(iddm)
             srcnam = 'IC'
c
c  --- call routine to fill the family of species names ---
c
             call filspddm(iaffec,caffec,cinflu,
     &                                    srcnam,idxddm,0,0,lout)
          enddo
c
c  --- loop over all of the DDM boundary condition speicies ---
c
          do iddm = 1,nbcddm
             cinflu = bcddmsp(iddm)
c
c  --- if the stratify boundary is off, just fill the one species name ---
c
             if( .NOT. lbndry ) then
                srcnam = 'BCALL'
c
c  --- call routine to fill the family of species names ---
c
                call filspddm(iaffec,caffec,cinflu,
     &                                    srcnam,idxddm,0,0,lout)
c
c  --- otherwise, loop over all of the boundary edges ---
c
             else
                do i=1,nedge
                   srcnam = 'BC'//edgnam(i)(1:3)
c
c  --- call routine to fill the family of species names ---
c
                   call filspddm(iaffec,caffec,cinflu,
     &                                   srcnam,idxddm,0,0,lout)
                enddo
             endif
          enddo
c
c  --- loop over all of the DDM emissions condition speicies ---
c
          do iddm = 1,nemddm
             cinflu = emddmsp(iddm)
             srcnam = 'EM'
c
c  --- call routine to fill the family of species names ---
c
              call filspddm(iaffec,caffec,cinflu,srcnam,idxddm,
     &                                    MAX(1,ngroup),nregin,lout)
          enddo
c
c  --- get then next affect species ---
c
      enddo
      ntotsp = idxddm
      ipttim = ntotsp + 1
      nsaspc = ntotsp
      write(idiag,*) ('-',i=1,79)
      write(idiag,*) 
c
c  --- return to calling routine ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
