      subroutine depprep(endtim,enddate)
c 
c-----CAMx v4.02 030709
c 
c     DEPPREP generates new deposition output species names and writes 
c     headers to new DEPOSITION output files.
c 
c     Copyright 2003
c     ENVIRON International Corporation
c           
c     Modifications: 
c        none
c 
c     Input arguments: 
c        endtim              model end time (HHMM)
c        enddate             model end date (YYJJJ)
c             
c     Output arguments: 
c        none
c             
c     Routines Called: 
c        ISTRLN
c             
c     Called by: 
c        STARTUP 
c
      include 'camx.prm'
      include 'camx.com'
      include 'camxfld.com'
      include 'filunit.com'
      include 'grid.com'
      include 'chmstry.com'
      include 'flags.com'
c
      character*4 ifile(10),note(60),dpspec(10,4*MXSPEC)
      integer enddate
      character*10 avfil,tmpnam
      integer ifgptr(MXGRID),ifglvl(MXGRID)
c
      data avfil /'AVERAGE   '/
      data nseg,izero,ione /1,0,1/
      data zero /0./
c
c-----Entry point
c
      idat1 = begdate
      idat2 = enddate
      tim1 = begtim/100.
      tim2 = endtim/100.
      if (.NOT.llatlon) then
        orgx = 1000.*xorg
        orgy = 1000.*yorg
        dx = 1000.*delx
        dy = 1000.*dely
        izone = 0
        if (lutm) izone = iuzon
      else
        orgx = xorg
        orgy = yorg
        dx = delx
        dy = dely
        izone = 0
      endif
      read(runmsg(1:60),'(60a1)') (note(n),n=1,60)
      read(avfil,'(10a1)') (ifile(n),n=1,10)
      do l = 1,navspc
        tmpnam = spname(lavmap(l))
        ll = istrln(tmpnam)
        tmpnam(ll+1:ll+3) = '_DV' 
        depsp(l) = tmpnam
        read(tmpnam,'(10a1)') (dpspec(n,l),n=1,10)
        tmpnam(ll+1:ll+3) = '_DD' 
        depsp(navspc+l) = tmpnam
        read(tmpnam,'(10a1)') (dpspec(n,navspc+l),n=1,10)
        tmpnam(ll+1:ll+3) = '_WD' 
        depsp(2*navspc+l) = tmpnam
        read(tmpnam,'(10a1)') (dpspec(n,2*navspc+l),n=1,10)
        tmpnam(ll+1:ll+3) = '_LC' 
        depsp(3*navspc+l) = tmpnam
        read(tmpnam,'(10a1)') (dpspec(n,3*navspc+l),n=1,10)
      enddo
c
c-----Coarse grid header
c
      rewind(idep)
      write(idep) ifile,note,nseg,4*navspc,idat1,tim1,idat2,tim2
      write(idep) zero,zero,izone,orgx,orgy,dx,dy,ncol(1),nrow(1),ione,
     &            izero,izero,zero,zero,zero
      write(idep) izero,izero,ncol(1),nrow(1)
      write(idep) ((dpspec(n,l),n=1,10),l=1,4*navspc)
      if( ngrid .EQ. 1 ) goto 9999
c
c-----Fine grid header
c
      do i = 1,ngrid 
        ifglvl(i) = 0 
      enddo 
      do i = 1,ngrid 
        do j = 1,nchdrn(i) 
          idch = idchdrn(j,i) 
          ifgptr(idch) = i - 1 
          ifglvl(idch) = ifglvl(idch) + 1 
          do k = 1,nchdrn(idch) 
            ifglvl(idchdrn(k,idch)) = ifglvl(idchdrn(k,idch)) + 1 
          enddo 
        enddo 
      enddo
      rewind(ifdep)
      write(ifdep) runmsg 
      write(ifdep) ngrid-1, 4*navspc 
      write(ifdep) (depsp(l),l=1,4*navspc) 
      do igrd = 2,ngrid
        write(ifdep) inst1(igrd),jnst1(igrd),inst2(igrd), 
     &               jnst2(igrd),meshold(igrd),meshold(igrd), 
     &               ncol(igrd),nrow(igrd),ione, 
     &               ifgptr(igrd),ifglvl(igrd) 
      enddo
c
 9999 continue
      return
      end
