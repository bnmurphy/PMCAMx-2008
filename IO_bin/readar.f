      subroutine readar(igrd,ncol,nrow,iarem,iout,aremis,dum1,dum2)
c
c-----CAMx v4.02 030709
c
c     READAR reads the time-variant records of the area source file for
c     the given grid and cycles through to current time/date to load 
c     area emission rates
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c          
c     Modifications:
c        10/20/00  Added check for negative emission rates
c        05/01/03  Time span of emissions must now match emiss update interval
c
c     Input arguments:
c        igrd                grid index
c        ncol                number of columns
c        nrow                number of rows
c        iarem               file unit for area emissions file
c        iout                file unit number for output file
c
c     Output arguments:
c        aremis              area emissions rate (mole/s)
c        dum1/dum2           dummy real variables
c
c     Routines Called:
c        none
c
c     Called by:
c        CAMx
c
      include 'camx.prm'
      include 'camx.com'
      include 'chmstry.com'
      include 'flags.com'
c
      character*4 arspec(10)
      real aremis(ncol,nrow,narspc)
c
c-----Entry point
c
      if (iarem.eq.0) goto 999
      kount = 1
 100  read(iarem,end=900) idat1,tim1,idat2,tim2 
      ichktm1 = NINT( 1000*(tim1) ) 
      ichktm2 = NINT( 1000*(tim2) ) 
      if( NINT(tim2) .EQ. 0 ) ichktm2 = 24000
      ichkems = NINT( 1000*(dtems/60.) )
      if( (ichktm2 - ichktm1) .NE. ichkems ) then
          write(iout,'(//,a)')'ERROR in READAR:'
          write(iout,*) 'Time interval in surface emissions file does'
          write(iout,*)  ' not match emissions update time interval.'
          write(iout,*) '   Beginning Date/Time (HHMM): ',idat1,tim1
          write(iout,*) '   Ending Date/Time    (HHMM): ',idat2,tim2
          write(iout,*) '   Emiss Input interval (min): ',dtems
          call camxerr()
      endif
      if( NINT(tim2) .EQ. 0) then
        tim2 = 24.
        idat2 = idat2 - 1
      endif
      tim1 = 100.*tim1 
      tim2 = 100.*tim2 
      do ll = 1,narspc-4
        read(iarem) idum,(arspec(i),i=1,10),
     &      ((aremis(i,j,ll),i=1,ncol),j=1,nrow)
      enddo
      write(iout,'(a40,2(f7.0,i8.5),a,i3)')
     &      'Read area source file at ',tim1,idat1,tim2,idat2,
     &      ' grid',igrd

c
c-----Check times only if LE1DAY = T, otherwise check both time and date
c
      if (le1day) then
        if (tim1.le.time .and. tim2.gt.time) goto 200
        if (tim1.gt.time) goto 900
      else
        if ((idat1.lt.date .or. (idat1.eq.date .and. tim1.le.time))
     &    .and. (idat2.gt.date .or. (idat2.eq.date .and. tim2.gt.time)))
     &    goto 200
      endif
      goto 100

c 
c-----Convert emission rates from moles/(dtems-hours) to moles/s for gas
c     or g/(dtems-hours) to g/s for aero species
c

c-----Split OLE2 into 4 VOCs--------------------------------------

 200  do j = 1,nrow
c        do i = 1,ncol
c          if (idmech.eq.6) then
c            aremis(i,j,narspc-3) = 0.16*aremis(i,j,21)
c            aremis(i,j,narspc-2) = 0.10*aremis(i,j,21)
c            aremis(i,j,narspc-1) = 0.38*aremis(i,j,21)
c            aremis(i,j,narspc) = 0.0*aremis(i,j,21)
c            aremis(i,j,21) = 0.36*aremis(i,j,21)
c          elseif (idmech.eq.5) then
c            aremis(i,j,narspc-3) = 0.16*aremis(i,j,21)
c            aremis(i,j,narspc-2) = 0.10*aremis(i,j,21)
c            aremis(i,j,narspc-1) = 0.38*aremis(i,j,21)
c            aremis(i,j,narspc) = 0.0*aremis(i,j,21)
c            aremis(i,j,21) = 0.36*aremis(i,j,21)
c          endif
c        enddo
      enddo

c-----------------------------------------------------------------

      do 10 l = 1,narspc
        do j = 1,nrow
          do i = 1,ncol
            if (aremis(i,j,l).lt.0.) then
              write(iout,'(//,a)') 'ERROR in READAR:'
              write(iout,'(a,i3)') 'Negative emissions for grid:',igrd
              write(iout,'(a,3i3)') '(i,j,l): ',i,j,l
              call camxerr()
            endif
            aremis(i,j,l) = aremis(i,j,l)/(60.*dtems) 
          enddo 
        enddo
 10   continue 
      goto 999
c
c-----End of file reached; if 1-day emissions requested, rewind and read 
c     through header once more.  Otherwise, report error and stop
c
 900  continue
      if (le1day) then
        if (kount.ge.2) then
          write(iout,'(//,a)') 'ERROR in READAR:'
          write(iout,*)'Cannot match model time with area source time'
          call camxerr()
        endif 
        rewind(iarem)
        read(iarem) idum 
        read(iarem) dum  
        read(iarem) idum  
        read(iarem) idum 
        kount = kount + 1
        goto 100
      else
        write(iout,'(//,a)') 'ERROR in READAR:'
        write(iout,*)'End of area source file reached'
        call camxerr()
      endif
c
 999  return
      end
