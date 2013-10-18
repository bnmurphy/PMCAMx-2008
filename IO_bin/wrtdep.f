      subroutine wrtdep(tim2,idat2,iunit,nox,noy,nsptmp,nspdry,
     &                  vdep,depfld)
c
c-----CAMx v4.02 030709
c 
c     WRTDEP writes deposition fields.
c 
c     Copyright 2003
c     ENVIRON International Corporation
c           
c     Modifications: 
c        none
c 
c     Input arguments:
c        tim2                output time (HHMM)
c        idat2               output date (YYJJJ)
c        iunit               output unit
c        nox                 number of cells in x-direction
c        noy                 number of cells in y-direction
c        nsptmp              number of dep field species
c        nspdry              number of dry dep species
c        vdep                dry deposition velocities
c        depfld              deposition field to output
c             
c     Output arguments: 
c        none
c             
c     Routines Called: 
c        none
c             
c     Called by: 
c        CAMx
c 
      include 'camx.prm'
      include 'camx.com'
      include 'grid.com'
      include 'chmstry.com'
c
      character*4 ispec(10,4*MXSPEC)
      real depfld(nox,noy,nsptmp),vdep(nox,noy,nspdry)
c
      data nseg /1/
c
c-----Entry point
c
c-----Determine time/date range
c
      idat1 = idat2
      etim = ANINT(tim2)/100.
      if( dtout .GE. 60.0 ) then
          btim = ANINT( 1000*(etim - ANINT(dtout)/60.) )/1000.
      else
          btim = ANINT( 1000*(etim - ANINT(dtout)/100.) )/1000.
      endif
      if (btim.lt.0.) then 
        btim = btim + 24. 
        idat1 = idat1 - 1 
      endif 
      idat3 = idat2
      etim3 = etim + 0.1
      if (etim3.gt.24.) then 
        etim3 = 24. - etim3 
        idat3 = idat3 + 1 
        if( MOD(idat3,1000) .GT. 365 ) then
            if( MOD(INT(idat3/1000),4) .EQ. 0 ) then
               if( MOD(idat3,1000) .EQ. 367 )
     &                     idat3 = (INT(idat3/1000)+1)*1000 + 1
            else
               idat3 = (INT(idat3/1000)+1)*1000 + 1
            endif
         endif
      endif 
c
c-----Write gridded deposition field
c
      do l = 1,4*navspc
        read(depsp(l),'(10a1)') (ispec(n,l),n=1,10)
      enddo
c      write(iunit) idat1,btim,idat2,etim
      do l = 1,navspc
        ll = lavmap(l)
c        write(iunit) nseg,(ispec(n,l),n=1,10),
c     &               ((vdep(i,j,ll),i=1,nox),j=1,noy)
      enddo
      do l = 1,nsptmp
        ll = l + navspc
c        write(iunit) nseg,(ispec(n,ll),n=1,10),
c     &               ((depfld(i,j,l),i=1,nox),j=1,noy)
      enddo
c
      return
      end
