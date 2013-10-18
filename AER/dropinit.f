      subroutine dropinit
      include 'aerpar.inc'
      include 'droppar.inc'
      include 'dropcom.inc'
      include 'section.inc'
      include 'camx_aero.inc'
c
c     CALCULATION OF AEROSOL DIAMETERS
c
      if (nsect .eq. 1) then
        daer(1)=dmin * 1.e-6
        daer(2)=dmax * 1.e-6
      else
        do i=1, nsect
          daer(i)=dsecf_c(i)*1.e-6
        enddo
      endif
c
c     CALCULATION OF FULL CENTER GRID
c
      do i=1, nsect
      dd(i) = daer(i)
      enddo
c
c     CALCULATION OF FULL GRID 
c
      do i=2, nsp_aq
      diameter(i) = SQRT(dd(i-1)*dd(i))
      enddo
      diameter(1)=dd(1)*dd(1)/diameter(2)
      diameter(nbounds)= dd(nsp_aq)*dd(nsp_aq)/diameter(nsp_aq)
c
c     DOCUMENT GRIDS USED
c
      f=1.e6          ! change units from m to um
c     
c     LOADING OF MOLECULAR WEIGHTS      
c
      wso2 = 64.
      wh2o2 = 34.
      whcho = 30.
      whcooh = 46.
      wnh3 = 17.
      whno3 = 63.
      whcl = 36.5
c
c      MOLECULAR WEIGHTS
c       
       wmol(1)= 81.0e0
       wmol(2)= 96.0e0
       wmol(3)= 47.0e0
       wmol(4)= 62.0e0
       wmol(5)= 62.0e0
       wmol(6)= 34.0e0
       wmol(7)= 48.0e0  ! was previously 60.0e0
       wmol(8)= 46.0e0
       wmol(9)= 30.0e0
       wmol(10)=46.0e0
       wmol(11)=48.0e0
       wmol(12)=121.0e0
       wmol(13)=76.0e0
       wmol(14)=48.0e0
       wmol(15)=35.5e0
       wmol(16)=17.0e0
       wmol(17)=33.0e0
       wmol(18)=62.0e0
       wmol(19)=18.0e0
       wmol(20)=47.0e0
       wmol(21)=32.0e0
       wmol(22)=35.5e0
       wmol(23)=52.50e0
       wmol(24)=96.0e0
       wmol(25)=112.0e0
       wmol(26)=113.0e0
       wmol(27)=111.0e0
       wmol(28)=60.00e0
       wmol(29)=18.0e0
c
       amol(1)= 55.85e0 
       amol(2)= 55.0e0 
       amol(3)= 23.0e0
c       
       gmol(1)=64.0
       gmol(2)=98.08 
       gmol(3)=47.02 
       gmol(4)=63.02 
       gmol(5)=44.01
       gmol(6)=34.02 
       gmol(7)=30.03 
       gmol(8)=46.00 
       gmol(9)=30.01 
       gmol(10)=46.01
       gmol(11)=48.00
       gmol(12)=121.05 
       gmol(13)=76.00 
       gmol(14)=48.00 
       gmol(15)=36.50
       gmol(16)=17.00 
       gmol(17)=33.01 
       gmol(18)=62.01 
       gmol(19)=17.00 
       gmol(20)=47.00
       gmol(21)=32.00 
       gmol(22)=18.00
      
      return
      end
