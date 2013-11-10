c     Modifications
c       02/10/2002 (tmg); change basis from wet to dry for calculation of
c                         the number of particles
c       04/15/2002 (tmg); call linint once instead of once for each species
c
      subroutine newdist(t,q)
c
      include 'dynamic.inc'
      real*8 q(ntotal)
      real*8 y0(nsec*nsp),y1(nsec*nsp)
      real*8 t,tempod,chmin,chmax
      real*8 diame(nsec),diamb(nsec),bound(nsec+1),secsize(nsec)
      real*8 ax(nsp),af(nsp)
      real*8 zz                                                   ! bkoo
c
c     dsec crossing
c     - With the new linint.f, this sorting is redundant, so removed - bkoo (08/25/03)

      chmax=0.0d0
      chmin=2.0d0
      do k=1,nsec
         qt(k)=0.0d0
         do kk=1,nsp
            qt(k)=qt(k)+q((k-1)*nsp+kk)
         enddo
cbk         if(qt(k).gt.1.0d-9) then
         if (.false.) then
            i=(k-1)*nsp
            zz=q(i+kna)/emw(kna)+q(i+knh4)/emw(knh4)              ! bkoo
            if(zz.eq.0.0) then                                    ! bkoo
             call get_param(igrdchm,ichm,jchm,kchm,iout,idiag)    ! bkoo
             write(iout,'(//,A)')'ERROR in NEWDIST: zz=0'         ! bkoo
             write(iout,*)' qt(',k,')=',qt(k)                     ! bkoo
             write(iout,*)' igrd,i,j,k: ',igrdchm,ichm,jchm,kchm  ! bkoo
             call camxerr()                                       ! bkoo
            else                                                  ! bkoo
            y0(k)=(2.0d0*q(i+kso4)/emw(kso4)+
     &            q(i+kno3)/emw(kno3)+q(i+kcl)/emw(kcl))/
     &            (q(i+kna)/emw(kna)+q(i+knh4)/emw(knh4))
            chmax=max(chmax,y0(k))
            chmin=min(chmin,y0(k))
            endif                                                 ! bkoo
         endif
      enddo

c
c     LINEAR INTERPOLATION ROUTINE
c
      do i=1,nsp
         do k=1,nsec
            y0((k-1)*nsp+i)=q((k-1)*nsp+i)
         enddo
      enddo
      call linint(dsec,y0,nsec,dsecf,y1,nsec,nsp,diamb,diame,bound,secsize,ax,af)
      do i=1,nsp
         do k=1,nsec
            q((k-1)*nsp+i)=y1((k-1)*nsp+i)
         enddo
      enddo

c      pause

c
c     check charge balances for each section
c
      do k=1,nsec
         qt(k)=0.0d0
         do kk=2,nsp ! tmg (02/19/02)
            qt(k)=qt(k)+q((k-1)*nsp+kk)
         enddo
cbk         if(qt(k).gt.1.0d-9) then
         if (.false.) then
            i=(k-1)*nsp
            zz=q(i+kna)/emw(kna)+q(i+knh4)/emw(knh4)              ! bkoo
            if(zz.eq.0.0) then                                    ! bkoo
             call get_param(igrdchm,ichm,jchm,kchm,iout,idiag)    ! bkoo
             write(iout,'(//,A)')'ERROR in NEWDIST: zz=0'         ! bkoo
             write(iout,*)' qt(',k,')=',qt(k)                     ! bkoo
             write(iout,*)' igrd,i,j,k: ',igrdchm,ichm,jchm,kchm  ! bkoo
             call camxerr()                                       ! bkoo
            else                                                  ! bkoo
            y1(k)=(2.0d0*q(i+kso4)/emw(kso4)+
     &            q(i+kno3)/emw(kno3)+q(i+kcl)/emw(kcl))/
     &            (q(i+kna)/emw(kna)+q(i+knh4)/emw(knh4))
cgy Fred Lurmann's instructions to changed tolerance, 19 Jan 2001 
cgy            if((y1(k).lt.(chmin-1.0d-8)).or.
cgy     &	      (y1(k).gt.(chmax+1.0d-8))) then
            if((y1(k).lt.(chmin-1.0d-6)).or.
     &         (y1(k).gt.(chmax+1.0d-6))) then
               call get_param(igrdchm,ichm,jchm,kchm,iout,idiag)
               write(iout,'(//,A)')'ERROR in NEWDIST: charge balance'
               write(iout,*)' ch.bal.(',k,')=',y1(k),chmin,chmax
cgy               write(iout,*) (chmin-1.0d-8)-y1(k),y1(k)-(chmax+1.0d-8)
               write(iout,*) (chmin-1.0d-6)-y1(k),y1(k)-(chmax+1.0d-6)
               write(iout,*) (dsec(m),m=1,nsec)
               write(iout,*) (dsecf(m),m=1,nsec+1)
               write(iout,*)' igrd,i,j,k: ',igrdchm,ichm,jchm,kchm
               call camxerr()
            endif 
            endif                                                 ! bkoo
         endif
      enddo
c
c     reset moving sectional diameters
c
      do i=1,nsec
         dsec(i)=sqrt(dsecf(i)*dsecf(i+1))
      enddo
c
c     recalculate number of particles in each section
c
      do i=1,nsec
         qn(i)=qt(i)/dsec(i)**3 ! calculate 0th moment(number of particles)
      enddo                     ! if density = 1g/cm3 units are particles/cm3
c
      return
      end
