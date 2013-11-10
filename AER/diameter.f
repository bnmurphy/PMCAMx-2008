c   1/25/2002 (tmg): separated into two functions, one for dry basis (ddiameter)
c                    and one for wet basis (wdiameter)
C ----------------------------------------------------------------------------
C   Diameter-calculates dry sectional diameters
      subroutine ddiameter(q)
      include 'dynamic.inc'
      real*8 q(ntotal)

        do i=1,nsec
          qt(i)=0.0d0
        enddo
        do i=1,nsec     ! Calculate total mass in each section
	  do ki=2,nsp ! tmg (01/25/02) make dry diameter
             qt(i)=qt(i)+q((i-1)*nsp+ki)  
          end do
	  if (imovsec .eq. 1) then
	  if(qn(i).gt.0.and.qt(i).gt.tinys) then ! avoid empty sections
 	     dsec(i)=(qt(i)/qn(i))**(0.33333)
cgy	     if(dsec(i).lt.0.005) write(90,*) 'at t= ',t1,
cgy     &         ' Diameter of section ',i,
cgy     &         ' is smaller than 0.005 um at Dp = ',dsec(i)
cgy	     if(dsec(i).gt.20.0) write(90,*) 'at t= ',t1,
cgy     &         ' Diameter of section ',i,
cgy     &         ' is larger than 20.0 um at Dp = ',dsec(i)
	  endif
	  endif
	end do

      return
      end



C ----------------------------------------------------------------------------
C   Diameter-calculates wet sectional diameters
      subroutine wdiameter(q)
      include 'dynamic.inc'
      real*8 q(ntotal)

        do i=1,nsec
          qt(i)=0.0d0
        enddo
        do i=1,nsec     ! Calculate total mass in each section
	  do ki=1,nsp ! tmg (01/25/02) make dry diameter
             qt(i)=qt(i)+q((i-1)*nsp+ki)  
          end do
	  if (imovsec .eq. 1) then
	  if(qn(i).gt.0.and.qt(i).gt.tinys) then ! avoid empty sections
 	     dsec(i)=(qt(i)/qn(i))**(0.33333)
cgy	     if(dsec(i).lt.0.005) write(90,*) 'at t= ',t1,
cgy     &         ' Diameter of section ',i,
cgy     &         ' is smaller than 0.005 um at Dp = ',dsec(i)
cgy	     if(dsec(i).gt.20.0) write(90,*) 'at t= ',t1,
cgy     &         ' Diameter of section ',i,
cgy     &         ' is larger than 20.0 um at Dp = ',dsec(i)
	  endif
	  endif
	end do

      return
      end
