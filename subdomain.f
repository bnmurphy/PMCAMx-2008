      subroutine subdomain(subd)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Subdomain returns an array compatible with the Eastern US CAMx grid domain
c	that is reduced in size to supress the affect of boundary conditions during
c	a budget analysis
c
c     Ben Murphy 9-23-09
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      integer subd(97,90)	!Columns,Rows or x,y
      integer i,j,x,y,a


C  This is the blank subdomain. Invoke either
C  the US or European versions to output mass 
C  fluxes that are less influenced by the boundary
C  conditions

C     INITIALIZE SUBD
      do i = 1,97
	do j = 1,90
	  subd(i,j) = 1 
	enddo
      enddo

      end
