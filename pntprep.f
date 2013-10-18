      subroutine pntprep(begtim,begdate,endtim,enddate)
c 
c-----CAMx v4.02 030709
c 
c     PNTPREP reads the header of binary point source emissions file,
c     initializes time-invariant point source variables, and maps the point
c     source species list to the internal CAMx species list
c                           
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c           
c     Modifications: 
c        11/06/01  Input dates are now Julian
c        05/01/03  Fixed bug that lost sources in first modeled cell
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
c        RDPTHDR
c        GEOPSP 
c             
c     Called by: 
c        STARTUP 
c 
      include 'camx.prm'
      include 'grid.com'
      include 'chmstry.com'
      include 'ptemiss.com'
      include 'filunit.com'
      include 'flags.com'
      include 'bndary.com'
c
      integer begdate,enddate
c
      dimension xloc(MXPTSRC),yloc(MXPTSRC),indxpt(MXPTSRC)
c
c-----Entry point
c
      call rdpthdr(xloc,yloc,begtim,begdate,endtim,enddate)
c
c-----Convert location to grid coordinates, convert vexit from m/hr to m/s 
c     and make sure stack diameter is positive (PiG flag is negative diameter)
c  
c-----Coarse grid
c
      if (llatlon) then
        do n = 1,nptsrc
          xstk(n,1) = xloc(n) - xorg
          ystk(n,1) = yloc(n) - yorg
        enddo
      else
        do n = 1,nptsrc
          xstk(n,1) = xloc(n)/1000. - xorg
          ystk(n,1) = yloc(n)/1000. - yorg
        enddo
      endif
c
      do 30 n=1,nptsrc
        ii = 1 + INT(xstk(n,1)/delx)
        jj = 1 + INT(ystk(n,1)/dely)
        isrc(n,1) = ii
        jsrc(n,1) = jj
        indxpt(n) = 1
c
        if (abs(dstk(n)).lt.0.01) then
          write(idiag,'(a,i5,2x,f7.3,a)') 'Zero stack diameter ',
     &           n,dstk(n),'  --- not modeled'
          indxpt(n) = 0
        endif
c
c-----Note sources outside domain
c
        if (ii.lt.1 .or. ii.gt.ncol(1) .or. jj.lt.1 .or.
     &      jj.gt.nrow(1)) then
          indxpt(n) = 0
          write(idiag,'(a,i5,2x,2i7,a)')'Source outside domain: ',
     &           n,ii,jj,'  --- not modeled'
        elseif (ibeg(jj).eq.-999 .or. jbeg(ii).eq.-999 .or.
     &          ii.lt.ibeg(jj) .or. ii.gt.iend(jj) .or.
     &          jj.lt.jbeg(ii) .or. jj.gt.jend(ii)) then
          indxpt(n) = 0
          write(idiag,'(a,i5,2x,2i7,a)')'Source in boundary cells: ',
     &           n,ii,jj,'  --- not modeled'
        endif
 30   continue
c
c-----Map source to coarse grid
c
      nsrc = 0
      do n=1,nptsrc
        if (indxpt(n).eq.1) then
          nsrc = nsrc + 1
          idsrc(nsrc,1) = n
          isrc(nsrc,1) = isrc(n,1) 
          jsrc(nsrc,1) = jsrc(n,1) 
        endif
      enddo
      nosrc(1) = nsrc
c
c-----Fine grids
c
      do igrd = 2,ngrid
        nsrc = 0
        do n = 1,nptsrc
          if (indxpt(n).eq.1) then
            xstk(n,igrd) = xstk(n,1) - (inst1(igrd) - 1)*delx
            ystk(n,igrd) = ystk(n,1) - (jnst1(igrd) - 1)*dely
            ii = 2 + INT(xstk(n,igrd)/delx*FLOAT( meshold(igrd) ) )
            jj = 2 + INT(ystk(n,igrd)/dely*FLOAT( meshold(igrd) ) )
c
            if (ii.gt.1 .and. ii.lt.ncol(igrd) .and.
     &          jj.gt.1 .and. jj.lt.nrow(igrd)) then
              nsrc = nsrc + 1
              idsrc(nsrc,igrd) = n
              isrc(nsrc,igrd) = ii
              jsrc(nsrc,igrd) = jj
            endif
          endif
          nosrc(igrd) = nsrc
        enddo
      enddo
c
      write(idiag,*)
      do igrd = 1,ngrid
        write(idiag,*) 'Point sources in grid #',igrd,' = ',nosrc(igrd)
c       write(idiag,'(4i5)') (n,idsrc(n,igrd),isrc(n,igrd),jsrc(n,igrd),
c    &                          n=1,nosrc(igrd))
      enddo
      write(idiag,*)
c
      return
      end
