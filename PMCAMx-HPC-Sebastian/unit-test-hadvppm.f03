PROGRAM unit_test_hadvppm
    IMPLICIT NONE
    INTEGER,PARAMETER :: NN = 90
    REAL,DIMENSION(NN) :: c1d, v1d, cinit, m1d, flxarr
    REAL :: dt , dx, flux1, flux2
    INTEGER :: ii = 0
    
    !Set up input scalars, vectors, and arrays
    DO ii = 1,NN
        c1d(ii) = 1.0 !Concentration (micrograms / m3)
        cinit(ii) = c1d(ii)
        v1d(ii) = 0.5 !Wind Speed (m/s)
        m1d(ii) = 1.0 !Map Scale Factor (important for projecting coordinates onto grid, set 1.0 for this exercise)
        flxarr(ii) = 0.0 !Output Vector
    END DO

    dt = 900 !Time Step (seconds)
    dx = 36000.0 !Grid size (m)
    flux1 = 0.0 !Output - flux at start boundary
    flux2 = 0.0 !Output - flux at end boundary

    !CALL 1-D ADVECTION SUBROUTINE!
    PRINT *,'Calling hadvppm...'
    CALL hadvppm(NN,dt,dx,c1d,v1d,m1d,flxarr,flux1, flux2)
    !call hadvppm(m2,dtuse,dy,c1d,v1d,a1d,av1d,flxarr,num1d)
    flxarr(NN) = 0.0 !The flux is BETWEEN grid cells, so for the lastcell it equals 0.

!    !Output all for the data to a text file
!    print *,'Writing output data...'
!    open(unit=20,file='Flux_Out.dat')
!
!    write(20,'(A4,2x,A9,2x,A9,2x,A4)'), 'Cell','Start Con','End Con','Flux'
!    do ii = 1,NN
!    write(20,'(I2,4x,E9.4,2x,E9.4,2x,E9.4)'), ii, cinit(ii), c1d(ii), flxarr(ii)
!    enddo
END PROGRAM unit_test_hadvppm
