c last modified by bkoo (Apr, 2000): fix other negative concentrations
c
cgy supress error messages to unit 90 (6/12/00)
c
c-----------------------------------------------------------------------
c  NEGCHK deals with negative concentrations
c  adds negative concentrations to the gas phase
c  adds species to maintain charge balance
c  this routine only needs to apply to evaporating species NH4+, NO3-, and Cl-
      subroutine negchk(t,q,nsection)

      include 'dynamic.inc'
      real*8 q(ntotal)

c     corrected by bkoo (02/09/01, 05/08/01)
      nsecxx = nsection
      ntotalxx = nsection * nsp + ngas

      prod=rgas*temp/(1.01325d5*pres) ! pres (bkoo, 06/09/00)
      ng=nsecxx*nsp

      do i=1,ntotalxx
        if(q(i).lt.1.d-30.and.q(i).gt.-1.d-30) q(i)=0.0D0
      enddo
      
      do isec=1,nsecxx
        indx=(isec-1)*nsp
c        if(q(indx+knh4).lt.-tinys)then        ! for NH4+
        if(q(indx+knh4).lt.0.d0) then         ! for NH4+
          if(q(indx+kno3).gt.q(indx+kcl))then ! maintain charge 
            q(indx+kno3)=q(indx+kno3)-q(indx+knh4)*emw(kno3)/emw(knh4)
            q(ng+ihno3)=q(ng+ihno3)+q(indx+knh4)/emw(knh4)*prod
          else ! add to larger of either NO3- or Cl-
            q(indx+kcl)=q(indx+kcl)-q(indx+knh4)*emw(kcl)/emw(knh4)
            q(ng+ihcl)=q(ng+ihcl)+q(indx+knh4)/emw(knh4)*prod
          endif                               ! maintain mass balance
          q(ng+inh3)=q(ng+inh3)+q(indx+knh4)*prod/gmw(Inh3)
          q(indx+knh4)=0.0d0
        endif
c        if(q(indx+kno3).lt.-tinys)then        ! for NO3-
        if(q(indx+kno3).lt.0.d0) then         ! for NO3-
          if(q(indx+knh4).gt.tinys.and.q(ng+inh3).gt.-q(indx+kno3)
     &     /emw(kno3)*prod)then ! maintain charge if NH4 and NH3 add NH4
            q(indx+knh4)=q(indx+knh4)-q(indx+kno3)*emw(knh4)/emw(kno3)
            q(ng+inh3)=q(ng+inh3)+q(indx+kno3)/emw(kno3)*prod
          elseif(q(indx+kcl).gt.-q(indx+kno3))then ! loose Cl-
            q(indx+kcl)=q(indx+kcl)+q(indx+kno3)*emw(kcl)/emw(kno3)
            q(ng+ihcl)=q(ng+ihcl)-q(indx+kno3)/emw(kno3)*prod
          elseif(q(ng+inh3).gt.-q(indx+kno3)/emw(kno3)*prod) then
c     if there is no Cl- to loose and there is NH3 add NH4
            q(indx+knh4)=q(indx+knh4)-q(indx+kno3)*emw(knh4)/emw(kno3)
            q(ng+inh3)=q(ng+inh3)+q(indx+kno3)/emw(kno3)*prod
          endif                               ! maintain mass balance
          q(ng+ihno3)=q(ng+ihno3)+q(indx+kno3)*prod/gmw(ihno3)
          q(indx+kno3)=0.0d0
        endif
c        if(q(indx+kcl).lt.-tinys)then         ! for Cl-
        if(q(indx+kcl).lt.0.d0) then          ! for Cl-
         if(q(indx+knh4).gt.tinys.and.q(ng+inh3).gt.-q(indx+kCL)
     &      /emw(kcl)*prod)then ! maintain charge if NH4 and NH3 add NH4
            q(indx+knh4)=q(indx+knh4)-q(indx+kcl)*emw(knh4)/emw(kcl)
            q(ng+inh3)=q(ng+inh3)+q(indx+kcl)/emw(kcl)*prod
          elseif(q(indx+kno3).gt.-q(indx+kcl))then ! loose NO3-
            q(indx+kno3)=q(indx+kno3)+q(indx+kcl)*emw(kno3)/emw(kcl)
            q(ng+ihno3)=q(ng+ihno3)-q(indx+kcl)/emw(kcl)*prod
          elseif(q(ng+inh3).gt.-q(indx+kcl)/emw(kcl)*prod) then
c     if there is no NO3- to loose and there is NH3 add NH4
            q(indx+knh4)=q(indx+knh4)-q(indx+kcl)*emw(knh4)/emw(kcl)
            q(ng+inh3)=q(ng+inh3)+q(indx+kcl)/emw(kcl)*prod
          endif                               ! maintain mass balance
          q(ng+ihcl)=q(ng+ihcl)+q(indx+kcl)*prod/gmw(ihcl)
          q(indx+kcl)=0.0d0
        endif
c     check negative concentrations for aerosols (Apr, 2000)
        do isp=1,nsp
          if(q(indx+isp).lt.0.0d0) then
c     subtract from gas to compensate the negative aerosol (bkoo, 10/07/03)
            idg = 0
            if(isp.eq.KSO4) idg = IH2SO4
            if(isp.eq.KNO3) idg = IHNO3
            if(isp.eq.KNH4) idg = INH3
            if(isp.eq.KCL)  idg = IHCL
            if(isp.ge.KAA41 .and. isp.lt.KAA41+ngas-4)
     &                      idg = ICA41+isp-KAA41
            if(idg.gt.0) then
          if(q(ng+idg)+q(indx+isp)*prod/gmw(idg).lt.-0.0001)            ! bkoo_dbg
     &    write(*,222)q(indx+isp),q(ng+idg)+q(indx+isp)*prod/gmw(idg),  ! bkoo_dbg
     &                isp,isec                                          ! bkoo_dbg
 222      format('NEG-AER: q=',G,' dgas=',G,' isp=',I3,' isec=',I3)     ! bkoo_dbg
              q(ng+idg)=dmax1(q(ng+idg)+q(indx+isp)*prod/gmw(idg)
     &                        ,1.d-12)
            endif
c
            q(indx+isp)=0.0d0
          endif
        enddo
      enddo
c     check negative concentrations for gases (Apr, 2000)
      do isp=1,ngas
        if(q(ng+isp).lt.0.0d0) then
        if(q(ng+isp).lt.-0.001) then                                    ! bkoo_dbg
        call get_param(igrdchm,ichm,jchm,kchm,iout,idiag)               ! bkoo_dbg
        write(*,223)q(ng+isp),isp,ichm,jchm,kchm                        ! bkoo_dbg
 223    format('NEG-GAS: q=',G,' isp,i,j,k=',4I4)                       ! bkoo_dbg
        endif                                                           ! bkoo_dbg
          q(ng+isp)=0.0d0
        endif
      enddo
      return
      end

