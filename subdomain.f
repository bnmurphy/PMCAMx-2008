      subroutine subdomain(subd)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Subdomain returns an array compatible with the Eastern US CAMx grid domain
c       that is reduced in size to supress the affect of boundary conditions during
c       a budget analysis
c
c     Ben Murphy 9-23-09
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none

      integer subd(150,162)    !Columns,Rows or x,y
      integer i,j,x,y,a

C     INITIALIZE SUBD
      do i = 1,150
        do j = 1,162
          subd(i,j) = 1
        enddo
      enddo

C     WESTERN BOUNDARY
      do j = 78,30,-1
        subd(8,j) = 1
      enddo
      do j = 30,9,-1
	x = 8-(j-30)
	y = j
	subd(x,y) = 1
      enddo

C     NORTHERN BOUNDARY
      do i = 8,90
	subd(i,78) = 1
      enddo

C     ATLANTIC BOUNDARY
      do j = 66,78
	subd(90,j) = 1
      enddo
      do j = 65,20,-1
	x = 90+int(j/3-65/3)
	y = j
	subd(x,y) = 1
      enddo
      do j = 3,20
	subd(75,j) = 1
      enddo

C     CARIBBEAN BOUNDARY
      do i = 30,60
	subd(i,9) = 1
      enddo
      do j = 3,9
	subd(60,j) = 1
      enddo
      do i = 60,75
	subd(i,3) = 1
      enddo

C     FILL IN DOMAIN WITH 1's SO MASSES CAN BE SCALED
      do i = 8,90
	do j = 64,78
	  subd(i,j) = 1
	enddo
      enddo
      do i = 8,75
	do j = 30,63
	  subd(i,j) = 1
	enddo
      enddo
      do i = 30,60
	do j = 9,29
	  subd(i,j) = 1
	enddo
      enddo
      do i = 60,75
	do j =3,29
	  subd(i,j) = 1
	enddo
      enddo
      do i = 8,29
	do j = 30,31+7-i,-1
	  subd(i,j) = 1
	enddo
      enddo
      do j = 20,65
	a = int(j/3-65/3)
	do i = 75,90+a
	  subd(i,j) = 1
	enddo
      enddo
      
C     PRINT SUBDOMAIN TO CONFIRM
c      print '(90(97I1/))', ((subd(i,j),i=1,97),j=90,1,-1)

      end
