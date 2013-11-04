       PROGRAM avgrads
c############################################################################
c#
c#    avgrads reads in radical concentration output from PMCAMx, averages
c#      the concentrations as perscribed by the user, and outputs them in
c#      an array readable by Matlab, Python, etc
c#
c#
c############################################################################
c#
c#    Written: Ben Murphy 03-26-09
c#
c############################################################################

      implicit none

      integer   ndate, nlay, ncol, nrow, mspc, nhr, nspc, nyr
      integer   ncalc, nrad, nt, nsite
      integer    ifin, ifout    !File unit markers

      parameter (ndate = 31)
      parameter (nyr   = 17)
      parameter (nlay  = 1)
      parameter (ncol  = 97)
      parameter (nrow  = 90)
      parameter (nrad  = 1 )
      parameter (nhr   = 24)
      parameter (nt    = 1000)
      parameter (nspc  = 3 )
      parameter (ncalc = 3 )
      parameter (nsite = 9 )
      parameter (ifin  = 30 )
      parameter (ifout = 31 )

      real     dt(nt)
      real     conc(ncol,nrow,nt)
      real     avgconc(ncol,nrow,ndate,nhr)
      real     a

      integer      t, tcount        !Looping Variables
      integer      irad, ilay, icol, irow, idate, it, isite  !Looping Variables
      integer      time1, time2, tcnt, iyr
      character*3  cdate(ndate), csite(nsite)
      character*4  crad(nrad), cyr(nyr), year
      character*2  clay(nlay), chour(nhr)  !String Variables for filenames
      character*199 scenario
      character*199 fname, fname2, indir, outdir

      integer lcalc(ncalc), f
      integer colsite(nsite), rowsite(nsite)

c      data cyr /'1991','1992','1993','1994','1995','1996','1997','1998','1999','2000'
c     &         ,'2051','2052','2053','2054','2055','2056','2057','2058','2059','2060'/
c      data cyr /'1996','2056','2060'/ !remove '182' to run
      data cyr /'1991','1992','1993','1994','1995','1997','1998','1999','2000'
     &         ,'2051','2052','2053','2054','2055','2057','2058','2059'/
c      data cyr /'2057'/ !1996 has trouble

      data cdate /'182','183','184','185','186','187','188','189','190','191','192',
     &            '193','194','195','196','197','198','199','200','201','202','203',
     &            '204','205','206','207','208','209','210','211','212'/

c      data cdate /'193','194','195','196','197','198','199','200','201',
c     &          '202','203','204','205','206','207','208','209'/
c      data cdate /'001','002'/ !  '003','004','005','006','007','008','009','010',
c     &            '011','012','013','014','015','016','017','018','019','020',
c     &            '021','022','023','024','025'/ !,'026','027','028','209','030',
c     &            '031'/

      data chour /'01','02','03','04','05','06','07','08','09','10','11','12',
     &            '13','14','15','16','17','18','19','20','21','22','23','24'/

c      data clay /'01','02','03','04','05','06','07','08','09','10','11','12',
c     &           '13','14'/
      data clay /'01'/ 

c      data crad /'OH','HO2','NO3','OHl'/
      data crad /'OH'/ ! change nrad if switched

      data csite /'PIT','ATL','CHI','NYC','TOT','DFN','LMI','NYO','LER'/
      data colsite /65,58,47,79,1,70, 47, 83, 58/
      data rowsite /51,28,52,55,1,38, 56, 55, 53/
c      data csite /'LMI','NYO','LER'/
c      data colsite /47, 83, 58/
c      data rowsite /56, 55, 53/

      do iyr = 1,nyr
      year = cyr(iyr)

c
c   Input/Output Structure
c
c      print *,'Enter the folder name (eg. SAPRC.July/): '
c      read(*,*), fname

       if (year .eq. '1994') then
       fname = 'PMCAMx2008_MGN3_'//year
       else
       fname = 'PMCAMx2008_MGN_'//year
       end if

c      print *,'Enter the scenario name (eg. semivol.fxEmisBC.060308): '
c      read(*,*), scenario
c      scenario = 'OH.111610'
      indir    = '/home/mday/PMCAMx/met_hd/Output/'// fname(1:INDEX(fname,' ')-1) //'/'
      outdir   = '/home/mday/PMCAMx/met_hd/processed_output/'//year
c     &            //fname(1:INDEX(fname,' ')-1)//'/'
c     &            //scenario

c
c   User-Defined Calculation Options
c
      print *,'Enter 1 to print daily avg output files (for maps): '
c      read(*,*), lcalc(1)
      lcalc(1) = 1
      print *,'Enter 1 to print hourly avg diurnal output files (for maps): '
c      read(*,*), lcalc(2)
      lcalc(2) = 0
      print *,'Enter 1 to print diurnal time series at specific sites: '
c      read(*,*), lcalc(3)
      lcalc(3) = 1

      do irad = 1,nrad
        print *,'Processing Radical: ',crad(irad)
        do ilay = 1,nlay
          print *,'Oh layer: ',clay(ilay)

c
c  Loop over days
c
          do idate = 1,ndate

c
c   Open Input Files
c
c           fname2 = indir(1:INDEX(indir,' ')-1)//scenario(1:INDEX(scenario,' ')-1)//
c     &       '/4rpos.baseE.'//cdate(idate)//'.'//
c     &       scenario(1:INDEX(scenario,' ')-1)//'.avrg'//clay(ilay)//'.'//crad(irad)

            fname2 = indir(1:INDEX(indir,' ')-1)//
     &      '/4rpos.baseE.'//cdate(idate)//'.'//
     &      fname(1:INDEX(fname,' ')-1)//'.avrg'//clay(ilay)//'.'//crad(irad)
            open(unit=ifin, file=fname2, form='UNFORMATTED',status='old')

       print *,'Retrieving Data from file: ',fname2

c
c   Initialize Concentration Array
c
        do icol = 1,ncol
          do irow = 1,nrow
            do it = 1,nt
              if (it.le.24) avgconc(icol,irow,idate,it) = 0
              conc(icol,irow,it) = 0.
            enddo
          end do
        end do
      do it = 1,nt
        dt(it) = 0.0
      enddo
      tcnt = 0

c
c   Read in all radical concentrations
c
      print *,'Reading in from input file'

      it = 1
      do while (1.eq.1)
             read(ifin,end=990) dt(it), ((conc(icol,irow,it),icol=1,ncol),irow=1,nrow)
        it = it + 1
      enddo
 990      continue
      tcnt = it-1
      close(ifin)

      print *,'Input file closed'

c
c   Perform hourly averages, calculations, etc on radical concentrations
c
      tcount = 0
      time2  = 1
      do it = 1,nt
c        !Bin time intervals to appropriate hour of the day
c        ! 0015 -> 1 hr; 1051 -> 11 hr; 2312.65 -> 24 hr
        time1 = aint(dt(it) / 100) + 1
c        print *,'dt = ',dt(it)
c        print *,'time1 = ',time1
c        !Count the number of intervals that occur within an hour
c        !and average by each hourly segment
        if (time1.gt.tcnt) then
            goto 991
        elseif (time1.eq.time2) then
          tcount = tcount + 1
          do icol = 1,ncol
            do irow = 1,nrow
              avgconc(icol,irow,idate,time1) = avgconc(icol,irow,idate,time1) +
     &            conc(icol,irow,it)
            enddo
          enddo
        elseif (time1.eq.time2+1) then
          do icol = 1,ncol
            do irow = 1,nrow
            avgconc(icol,irow,idate,time2) = avgconc(icol,irow,idate,time2)/tcount
                avgconc(icol,irow,idate,time1) = avgconc(icol,irow,idate,time1) +
     &              conc(icol,irow,it)
            enddo
          enddo 
          time2 = time2 + 1
          tcount = 1
        endif
      enddo  !it
 991    continue

        print *,'avgconc = ',avgconc(65,55,idate,4)
      enddo       !idate

c
c   Write to new output files
c
      print *,'Writing to output files for layer: ',clay(ilay),'  day: ','  radical: ',crad(irad)
      if (lcalc(1)) then         !Daily Map Files
        print *,'Calc = 1. Daily Map Files'
        fname2 = outdir(1:INDEX(outdir,' ')-1)//'/'//
     &            crad(irad)(1:INDEX(crad(irad),' ')-1)//'.daily'
c     &           //clay(ilay)
              open(unit=ifout, file=fname2, form='FORMATTED')
c        write(ifout, 46), 'Date', 'X','Y','CONC'
      
        do idate = 1,ndate
          do icol = 1,ncol
            do irow = 1,nrow
          a = 0
          do it = 1,nhr
               a = a + avgconc(icol,irow,idate,it)/nhr
          enddo
          !Write daily avg concentration and transform from ppm -> ppt
           write(ifout,47) cdate(idate), icol ,irow, a*10**6
        enddo
          enddo
            enddo

        close(ifout)
      endif
      if (lcalc(2)) then   !Diurnal Avg Map Files
        print *,'Calc = 2. Diurnal Avg Maps'
        fname2 = outdir(1:INDEX(outdir,' ')-1)//'/'//crad(irad)(1:INDEX(crad(irad),' ')-1)//
     &            '.diurn.'//clay(ilay)
              open(unit=ifout, file=fname2, form='FORMATTED')
c        write(ifout, 46), 'Hour', 'X','Y','CONC'
      
        do it = 1,nhr
          do icol = 1,ncol
            do irow = 1,nrow
          a = 0
              do idate = 1,ndate
            a = a + avgconc(icol,irow,idate,it)/ndate
          enddo
c            !Write daily avg concentration and transform from ppm -> ppt
             write(ifout,47) chour(it), icol ,irow, a*10**6
        enddo
          enddo
        enddo

        close(ifout)
      endif
      if (lcalc(3)) then   !Diurnal Conc'ns at specific sites
        print *,'Calc = 3. Specific site diurnal paterns'
        do isite = 1,nsite
          fname2 = outdir(1:INDEX(outdir,' ')-1)//'/'//crad(irad)(1:INDEX(crad(irad),' ')-1)//
     &            '.'//csite(isite)//'.'//clay(ilay)
                open(unit=ifout, file=fname2, form='FORMATTED')
c          write(ifout, 48), 'Date','Hour', 'X','Y','CONC'
      
          do idate = 1,ndate
          do it = 1,nhr
        a = 0
        if (csite(isite).eq.'TOT') then
            do icol = 1,ncol
              do irow = 1,nrow
            a = a + avgconc(icol,irow,idate,it)/nrow/ncol
              enddo
            enddo
        else
            a = avgconc(colsite(isite),rowsite(isite),idate,it)
        endif
c        !Write daily avg concentration and transform from ppm -> ppt
             write(ifout, 49) cdate(idate),chour(it), colsite(isite),rowsite(isite), a*10**6
          enddo
          enddo
        enddo !isite
    
        close(ifout)

      endif  !Calculation Conditional

      enddo       !ilay
      enddo    !irad

      print *,'All Days Done!!!!'

 46   format(A4,2x,A1,3x,A1,3x,A4)
 47   format(A3,3x,i2,2x,i2,2x,e12.5)
 48   format(A4,2x,A4,2x,A1,3x,A1,3x,A4)
 49   format(A3,3x,A3,3x,i2,2x,i2,2x,e12.5)

      end do
      end
