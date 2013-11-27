program UnitTestHadvppm
    !IMPLICIT NONE
    integer :: ii, nn, status, iEnd
    character(len=10) :: arg
    real :: dt,dx,sigma
    real, allocatable, dimension(:) :: con,conInit,vel,area,areav,flxarr
    
    !read input arguments from command line, 
    !iEnd: number of iterations, nn length of arrays
    call get_command_argument(1,arg)
    read(arg,*) iEnd
    call get_command_argument(2,arg)
    read(arg,*) nn
    
    !allocate memory for dynamic arrays
    allocate(con(nn), stat = status)
    allocate(conInit(nn), stat = status)
    allocate(vel(nn), stat = status)
    allocate(area(nn), stat = status)
    allocate(areav(nn), stat = status)
    allocate(flxarr(nn), stat = status)
    
    
    !build test Gaussian test function  
    dt = 1.0
    dx = 36000
    sigma = 1.0
    
    do ii=1,nn
        con(ii) = EXP((-(ii-10) * (ii-10)/2*sigma*sigma))
        vel(ii) = 150
        area(ii) = 1.
        areav(ii) = 1.
        flxarr(ii) = 0.
    end do
    
    conInit = con
    
    !repeatedly call the advection routine
    !do ii = 1,iEnd
     !   call hadvppm(nn,dt,dx,con,vel,area,areav,flxarr)
    !end do
    
    !Output all for the data to a text file
!    print *,'Writing output data...'
!    open(unit=20,file='Flux_Out.dat')
!
!    
!    do ii = 1,nn
!        write(20,'(I9,4x,E9.4,2x,E9.4,2x,E9.4)') ii, conInit(ii), con(ii), flxarr(ii)
!    enddo
    
    !deallocate memory
    deallocate(con, stat = status)
    deallocate(conInit, stat = status)
    deallocate(vel, stat = status)
    deallocate(area, stat = status)
    deallocate(areav, stat = status)
    deallocate(flxarr, stat = status)
    
    
end program UnitTestHadvppm
