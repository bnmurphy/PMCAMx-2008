subroutine xyadvec(igrd,xyordr,ncol,nrow,nlay,nspc,nsen,nadv, deltat,&
    dx,dy,windu,windv,depth,mapscl,conc,fluxes,sens,tarray2,isaptr,ipa_cel)

!include "camx.prm"
!include "bndary.com"
!include "chmstry.com"
!include "filunit.com"
!include "flags.com"
!include "tracer.com"
!include "procan.com"
    
INTEGER, DIMENSION(ncol,nrow,nlay) :: ipa_cel

integer xyordr
integer nadv(nlay)

dimension conc(ncol,nrow,nlay,nspc) !what is this?
dimension sens(ncol,nrow,nlay,nsen) !what is this?
real windu(ncol,nrow,nlay),windv(ncol,nrow,nlay), depth(ncol,nrow,nlay), mapscl(ncol,nrow), dx(nrow)
real c1d(MX1D),v1d(MX1D),flxarr(MX1D),m1d(MX1D),saflux(MX1D)
real fluxo3(MX1D), fluxvoc(MX1D), fluxnox(MX1D)
real cnco3(MX1D), cncvoc(MX1D), cncnox(MX1D)
real c1d0(MX1D),fpc(MX1D),fmc(MX1D)
real sen1d(MX1D,MXTRSP)
!real*8 fluxes(nspc,11),flux1,flux2
!real*8 fluxtmp(MXSPEC,8,MXLAYA)
REAL fluxes(nspc,11),flux1,flux2
REAL, DIMENSION (MXSPEC,8,MXLAYA) :: fluxtmp = 0.
dimension tarray2(2) !what is this?
REAL, DIMENSION (MX1D):: fc1, fc2

INTEGER :: i1, i2 !some looping variables that were previously undefined

if( xyordr .eq. 0 ) goto 200 !WTF?

!-----Advection in x-direction

100  CONTINUE !WTF?
!@todo
!$omp master
!write(*,'(a20,$)') 'x advection ......'
!write(iout,'(a20,$)') 'x advection ......'
tcpu = dtime(tarray2)

!@todo
!$omp end master
!$omp parallel default(shared)
!$omp&  private(i,j,k,i1,i2,ispc,l,c1d,v1d,m1d,nn,iddm,c1d0,ioff,
!$omp          sen1d,flxarr,flux1,flux2,saflux,fpc,fmc,fc1,fc2,
!$omp&          ipa_idx,fluxnox,fluxvoc,fluxo3,cncnox,cncvoc,cnco3,
!$omp&          istep,dtuse)

!$omp do schedule(dynamic)

DO 20 k = 1,nlay
    dtuse = deltat/nadv(k)
    DO 21 istep = 1,nadv(k)
        DO 10 j = 2,nrow-1
            i1 = 1
            i2 = ncol
            IF (igrd.eq.1) THEN
                IF (ibeg(j).eq.-999) goto 10
                i1 = ibeg(j) - 1
                i2 = iend(j) + 1
            END IF
            l = 0
            DO i = i1,i2-1
                l = l + 1
                v1d(l) = 2.*windu(i,j,k)/(mapscl(i+1,j) + mapscl(i,j))
                m1d(l) = mapscl(i,j)*mapscl(i,j)
            END DO
            l = l + 1
            v1d(l) = windu(i2,j,k)
            m1d(l) = mapscl(i2,j)*mapscl(i2,j)
            nn = i2 - i1 + 1
            ! Source Apportion Begin
            ! initialize the tracer fluxes/concs to zero
            DO i=1,MX1D
                fluxo3(i) = 0.
                fluxvoc(i) = 0.
                fluxnox(i) = 0.
                cnco3(i) = 0.
                cncvoc(i) = 0.
                cncnox(i) = 0.
            END DO
            ! Source Apportion End
            
            DO 22 ispc = 1,nspc
                l = 0
                DO i = i1,i2
                    l = l + 1
                    c1d(l) = conc(i,j,k,ispc)*dy*depth(i,j,k)
                END DO
                IF (igrd.eq.1 .and. v1d(1).lt.0.) c1d(1) = c1d(2) !@todo check if that syntax is correct
                IF (igrd.eq.1 .and. v1d(nn-1).gt.0.) c1d(nn) = c1d(nn-1)!@todo check if that syntax is correct
                
                !DDM Begin
                IF( lddm ) THEN
                    DO iddm=1,nddmsp
                        l = 0
                        DO i = i1,i2
                            l = l + 1
                            c1d0(l) = c1d(l)
                            ioff = iptddm(ispc)+iddm-1
                            sen1d(l,iddm) = sens(i,j,k,ioff)*dy*depth(i,j,k)
                        END DO

                        IF (igrd.eq.1 .and. v1d(1).lt.0.) sen1d(1,iddm) = sen1d(2,iddm)
                        IF (igrd.eq.1 .and. v1d(nn-1).gt.0.) sen1d(nn,iddm) = sen1d(nn-1,iddm)
                    END DO
                END IF                
                !DDM End
                
                IF (iadvct.eq.2) THEN
                    call hadvbot(nn,dtuse,dx(j),c1d,v1d,m1d,flxarr,flux1,flux2,saflux,fpc,fmc,fc1,fc2)
                ELSE IF( iadvct .eq. 3) THEN
                    call hadvppm(nn,dtuse,dx(j),c1d,v1d,m1d,flxarr,flux1,flux2,saflux,fc1,fc2)
                END IF
                
                !DDM Begin
                IF( lddm .AND. iadvct .eq. 2 ) THEN
                    call bottddm(nn,dtuse,dx(j),sen1d,nddmsp,c1d0,c1d,fpc,fmc,v1d,m1d)
                END IF
                !DDM End 
                !WTF is this DDM ENd end Begin about? 
                
                !Process Analysis Begin
                !are such comments really necessary?

                !@todo WTF?, what are those recursive function calls?
                IF( lipr ) THEN
                    !-----Change from X-direction horizontal advection
                    l = 1
                    DO i=i1+1,i2-1
                        l = l+1
                       IF( ipa_cel(i,j,k) .GT. 0 ) THEN
                          ipa_idx = ipa_cel(i,j,k)
                          !-----Flux at west boundary
                          cipr(IPR_WADV, ipa_idx, ispc) = cipr(IPR_WADV, ipa_idx, ispc) + fc1(l)/dy/depth(i,j,k)          
                          !-----Flux at east boundary
                          cipr(IPR_EADV, ipa_idx, ispc) = cipr(IPR_EADV, ipa_idx, ispc) + fc2(l)/dy/depth(i,j,k)
                          !-----Average volume
                          cipr(IPR_VOL, ipa_idx, ispc) = cipr(IPR_VOL, ipa_idx, ispc) + dx(j)*dy*depth(i,j,k)
                          npastep(ipa_idx,ispc) = npastep(ipa_idx,ispc) + 1
                          END IF
                    END DO
                END IF
                !Process Analysis End
                
                !Source Apportion Begin
                IF( ltrace .AND. tectyp .NE. RTRAC ) THEN
                    l = 0
                    DO i = i1,i2
                        l = l + 1
                        IF( lo3sp(ispc) ) THEN
                            fluxo3(l) = fluxo3(l) + saflux(l)
                            cnco3(l) = cnco3(l) + conc(i,j,k,ispc)
                        ELSE IF( lnoxsp(ispc) ) THEN
                            fluxnox(l) = fluxnox(l) + saflux(l)
                            cncnox(l) = cncnox(l) + conc(i,j,k,ispc)
                        ELSE IF ( lvocsp(ispc) ) THEN
                            fluxvoc(l) = fluxvoc(l) + saflux(l) * crbnum(ispc)
                            cncvoc(l) = cncvoc(l) + conc(i,j,k,ispc) * crbnum(ispc)
                        END IF
                    END DO
                END IF
                !Source Apportion End

                l = 1
                DO i = i1+1,i2-1
                    l = l + 1
                    conc(i,j,k,ispc) = c1d(l)/dy/depth(i,j,k)
                END DO
                
                !DDM Begin
                IF( lddm ) THEN
                    DO iddm=1,nddmsp
                        l = 1
                        DO i = i1+1,i2-1
                            l = l + 1
                            ioff = iptddm(ispc)+iddm-1
                            sens(i,j,k,ioff) = sen1d(l,iddm)/dy/depth(i,j,k)
                        END DO
                    END DO
                END IF
                !DDM End

                !-----Sum up fluxes in east and west sides

                IF (flux1.gt.0.0) THEN
                    fluxtmp(ispc,7,k)  = fluxtmp(ispc,7,k)  + flux1*dtuse
                ELSE
                    fluxtmp(ispc,8,k)  = fluxtmp(ispc,8,k)  + flux1*dtuse
                END IF
                
                IF (flux2.lt.0.0) THEN
                    fluxtmp(ispc,5,k)  = fluxtmp(ispc,5,k)  - flux2*dtuse
                ELSE
                    fluxtmp(ispc,6,k)  = fluxtmp(ispc,6,k)  - flux2*dtuse
                END IF
                22 CONTINUE !@todo, what does this mean?
                
                !Source Apportion Begin
                !call routine to update the tracer concentrations based on the calculated fluxes
                IF( ltrace .AND. tectyp .NE. RTRAC ) THEN
                    call xfluxsa(ncol,nrow,nlay,ntotsp, ptconc(isaptr),i1+1,i2-1,j,k,dy,depth,mapscl,&
                                 fluxo3,cnco3,fluxnox,cncnox,fluxvoc,cncvoc)
                END IF
                
                !Source Apportion End 
                !  --- next row, non-parallelized loop
                10  CONTINUE !@todo, what is this?
                !  --- next layer ---
                21  CONTINUE !@todo, what is this?
                20  CONTINUE !@todo, what is this?

!@todo, shouldn't there be something like this here?
!            END DO
!        END DO
!    END DO
!END DO

!  --- end of parallelized loop ---, @todo, WTF?

!$omp end parallel
!$omp master
                
tcpu = dtime(tarray2)
write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
write(iout,'(a,f10.3)') '   CPU = ', tarray2(1)
!$omp end master

if (xyordr.eq.0) goto 300

!-----Advection in y-direction

200  CONTINUE !

!$omp master
!write(*,'(a20,$)') 'y advection ......'
!write(iout,'(a20,$)') 'y advection ......'
!$omp end master

!@todo
!$omp parallel default(shared)
!$omp&  private(i,j,k,j1,j2,ispc,l,c1d,v1d,m1d,nn,iddm,c1d0,ioff,
!$omp&          sen1d,flxarr,flux1,flux2,saflux,fpc,fmc,fc1,fc2,
!$omp&          ipa_idx,fluxnox,fluxvoc,fluxo3,cncnox,cncvoc,cnco3,
!$omp&          istep,dtuse)

!$omp do schedule(dynamic)

do 40 k = 1,nlay
dtuse = deltat/nadv(k)
do 41 istep = 1,nadv(k)

do 30 i = 2,ncol-1
j1 = 1
j2 = nrow
if (igrd.eq.1) then
if (jbeg(i).eq.-999) goto 30
j1 = jbeg(i) - 1
j2 = jend(i) + 1
endif
l = 0
do j = j1,j2-1
l = l + 1
v1d(l) = 2.*windv(i,j,k)/(mapscl(i,j+1) + mapscl(i,j))
m1d(l) = mapscl(i,j)*mapscl(i,j)
enddo
l = l + 1
v1d(l) = windv(i,j2,k)
m1d(l) = mapscl(i,j2)*mapscl(i,j2)
nn = j2 - j1 + 1

!======================== Source Apportion Begin =======================

 ! --- initialize the tracer fluxes/concs to zero ----

do j=1,MX1D
fluxo3(j) = 0.
fluxvoc(j) = 0.
fluxnox(j) = 0.
cnco3(j) = 0.
cncvoc(j) = 0.
cncnox(j) = 0.
enddo

!========================= Source Apportion End ========================

do 42 ispc = 1,nspc

l = 0
do j = j1,j2
l = l + 1
c1d(l) = conc(i,j,k,ispc)*dx(j)*depth(i,j,k)
enddo
if (igrd.eq.1 .and. v1d(1).lt.0.) c1d(1) = c1d(2)
if (igrd.eq.1 .and. v1d(nn-1).gt.0.) c1d(nn) = c1d(nn-1)

!======================== DDM Begin =======================

if( lddm ) then
do iddm=1,nddmsp
 l = 0
 do j = j1,j2
   l = l + 1
   c1d0(l) = c1d(l)
   ioff = iptddm(ispc)+iddm-1
   sen1d(l,iddm) = sens(i,j,k,ioff)*dx(j)*depth(i,j,k)
 enddo
 if (igrd.eq.1 .and. v1d(1).lt.0.) sen1d(1,iddm) = sen1d(2,iddm)
 if (igrd.eq.1 .and. v1d(nn-1).gt.0.) sen1d(nn,iddm) = sen1d(nn-1,iddm)
enddo
endif

!======================== DDM End =======================

IF (iadvct.eq.2) THEN
    call hadvbot(nn,dtuse,dy,c1d,v1d,m1d,flxarr,flux1,flux2,saflux,fpc,fmc,fc1,fc2)
ELSE IF ( iadvct .eq. 3) THEN
    call hadvppm(nn,dtuse,dy,c1d,v1d,m1d,flxarr,flux1, flux2,saflux,fc1,fc2)
END IF

!======================== DDM Begin =======================


IF( lddm .AND. iadvct .eq. 2 ) THEN
    call bottddm(nn,dtuse,dy,sen1d,nddmsp,c1d0,c1d,fpc,fmc,v1d,m1d)
END IF

!======================== DDM End =======================
!======================== Process Analysis Begin ====================================

if( lipr ) then

!-----Change from Y-direction horizontal advection
 
l = 1
do j=j1+1,j2-1
  l = l+1

  if( ipa_cel(i,j,k) .GT. 0 ) then
     ipa_idx = ipa_cel(i,j,k)

!-----Flux at south boundary

     cipr(IPR_SADV, ipa_idx, ispc) = cipr(IPR_SADV, ipa_idx, ispc) + fc1(l)/dx(j)/depth(i,j,k)

!-----Flux at north boundary

     cipr(IPR_NADV, ipa_idx, ispc) = cipr(IPR_NADV, ipa_idx, ispc) + fc2(l)/dx(j)/depth(i,j,k)

  endif
enddo
endif

!========================= Process Analysis End =====================================


!======================== Source Apportion Begin =======================

if( ltrace .AND. tectyp .NE. RTRAC ) then
l = 0
do j = j1,j2
 l = l + 1
 if( lo3sp(ispc) ) then
     fluxo3(l) = fluxo3(l) + saflux(l)
     cnco3(l) = cnco3(l) + conc(i,j,k,ispc)
 else if( lnoxsp(ispc) ) then
     fluxnox(l) = fluxnox(l) + saflux(l)
     cncnox(l) = cncnox(l) + conc(i,j,k,ispc)
 else if( lvocsp(ispc) ) then
     fluxvoc(l) = fluxvoc(l) + saflux(l) * crbnum(ispc)
     cncvoc(l) = cncvoc(l) + conc(i,j,k,ispc) * crbnum(ispc)
 endif
enddo
endif

!======================== Source Apportion End =======================

l = 1
do j = j1+1,j2-1
l = l+1
conc(i,j,k,ispc) = c1d(l)/dx(j)/depth(i,j,k)
enddo

!======================== DDM Begin =======================

if( lddm ) then
do iddm=1,nddmsp
 l = 1
 do j = j1+1,j2-1
   l = l + 1
   ioff = iptddm(ispc)+iddm-1
   sens(i,j,k,ioff) = sen1d(l,iddm)/dx(j)/depth(i,j,k)
 enddo
enddo
endif

!======================== DDM End =======================


!-----Sum up fluxes in north and south sides

if (flux1.gt.0.0) then
fluxtmp(ispc,3,k)  = fluxtmp(ispc,3,k)  + flux1*dtuse
else
fluxtmp(ispc,4,k)  = fluxtmp(ispc,4,k)  + flux1*dtuse
endif
if (flux2.lt.0.0) then
fluxtmp(ispc,1,k)  = fluxtmp(ispc,1,k)  - flux2*dtuse
else
fluxtmp(ispc,2,k)  = fluxtmp(ispc,2,k)  - flux2*dtuse
endif
42      continue

!======================== Source Apportion Begin =======================

!  --- call routine to update the tracer concentrations
!      based on the calculated fluxes ---

IF( ltrace .AND. tectyp .NE. RTRAC ) THEN
    call yfluxsa(ncol,nrow,nlay,ntotsp, &
    ptconc(isaptr),i,j1+1,j2-1,k,dx,depth,mapscl,&
    fluxo3,cnco3,fluxnox,cncnox,fluxvoc,cncvoc)
endif

!======================== Source Apportion End =======================

30    continue

!  --- next layer, end of parallelized loop ---

41  continue
40  continue

!$omp end parallel

!$omp master
tcpu = dtime(tarray2)
write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
write(iout,'(a,f10.3)') '   CPU = ', tarray2(1)
!$omp end master

if (xyordr.eq.0) goto 100

300  continue

!  ---- put fluxes into global array ----

do i=1,nspc
do j=1,8
do k=1,nlay
fluxes(i,j) = fluxes(i,j) + fluxtmp(i,j,k)
enddo
enddo
enddo

call flush(6)
call flush(iout)
return
end


