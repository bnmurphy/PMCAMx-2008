      subroutine chemrxn(igrd)
c
c-----CAMx v4.02 030709
c
c     CHEMRXN passes arrays from common decks to CHEMDRIV
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c            
c     Modifications:  
c        none
c
c     Input arguments:
c        igrd                grid index
c
c     Output arguments:
c        none
c
c     Routines called:
c        CHEMDRIV
c        MASSUM
c
c     Called by:
c        CAMx
c        NESTING
c
      include 'camx.prm'
      include 'camx.com'
      include 'camxfld.com'
      include 'grid.com'
      include 'flags.com'
      include 'filunit.com'
      include 'chmstry.com'

      real masschem(nspec*nlay(igrd)), masspart(nspec*nlay(igrd))
c
c======================== Process Analysis Begin ====================================
c
      include "procan.com"
c
c========================= Process Analysis End =====================================
c
      real chemtime_0, chemtime, aertime

c
c-----Entry point
c
      call massum(igrd,nspec,ncol(igrd),nrow(igrd),nlay(igrd),
     &            deltax(1,igrd),deltay(igrd),depth(iptr3d(igrd)),
     &            conc(iptr4d(igrd)),xmstmp(1,igrd), subxmstmp(1,igrd))
c
      if (lchem) then
        write(*,'(a20,$)') 'chemdriv ......'
        write(iout,'(a20,$)') 'chemdriv ......'

 	do l = 1,nspec*nlay(igrd)
	  masschem(l) = 0
	  masspart(l) = 0
	enddo
 
        chemtime_0 = 0
        chemtime = 0
        aertime = 0

        call chemdriv(igrd,ncol(igrd),nrow(igrd),nlay(igrd),
     &                deltat(igrd),itzon,idfin(iptr2d(igrd)),
     &                fcloud(iptr3d(igrd)),cldtrns(iptr3d(igrd)),
     &                water(iptr3d(igrd)),tempk(iptr3d(igrd)),
     &                press(iptr3d(igrd)),height(iptr3d(igrd)),
     &                cwc(iptr3d(igrd)),
     &                conc(iptr4d(igrd)),cncrad(iptrad(igrd)),
     &                cellat(iptr2d(igrd)),cellon(iptr2d(igrd)),
     &                ldark(iptr2d(igrd)),l3davg,
     $                iptr2d(igrd),
     &                ipsa3d(igrd),ipacl_3d(iptr3d(igrd)),
     &                masschem,masspart, deltax(1,igrd),      !<- BNM 6/30/09
     &                deltay(igrd), depth(iptr3d(igrd)),      !<- BNM 6/2/09
     &                chemtime_0, chemtime, aertime)
        !tcpu = dtime(tarray2)
        !write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
        !write(iout,'(a,f10.3)') '   CPU = ', tarray2(1)
        write(*,'(a,f10.3)') '   CPU = ', chemtime_0+chemtime+aertime
        write(iout,'(a,f10.3)') '   CPU = ', chemtime_0+chemtime+aertime
        write(*,'(a,f10.3)') '   Chemistry = ', chemtime
        write(iout,'(a,f10.3)') '   Chemistry = ', chemtime
        write(*,'(a,f10.3)') '   Aerosols = ', aertime
        write(iout,'(a,f10.3)') '   Aerosols = ', aertime


        call flush(6)
        call flush(iout)
      endif
c
      call massum(igrd,nspec,ncol(igrd),nrow(igrd),nlay(igrd),
     &            deltax(1,igrd),deltay(igrd),depth(iptr3d(igrd)),
     &            conc(iptr4d(igrd)),xmass(1,igrd), subxmass(1,igrd))


C    SUM UP BUDGET ANALYSIS ARRAYS
C	XMSCHEM - total chemistry + partitioning
C	XMSJUSTCHEM - just change due to chemistry
C	XMSPART - change due to partitioning
      do l = 1,nspec*nlay(igrd) 
        xmschem(l,igrd) = xmschem(l,igrd) + subxmass(l,igrd) - 
     &                    subxmstmp(l,igrd)
	xmsjustchem(l,igrd) = xmsjustchem(l,igrd) + masschem(l)
	xmspart(l,igrd) = xmspart(l,igrd) + masspart(l)

      enddo
c
      return
      end
