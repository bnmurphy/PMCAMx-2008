program UnitTestHadvppm
    !IMPLICIT NONE
    integer :: ii, nn, mynn
    character(len=10) :: arg
    real :: dt,dx,sigma
    real, dimension(nn) :: con,vel,area,areav,flxarr
    
    call get_command_argument(1,arg)
    read(arg,*) ii
    call get_command_argument(2,arg)
    read(arg,*) nn
    
    mynn = nn
    
    
    dt = 1.0
    dx = 360000
    sigma = 0.1
    
    do ii=1,nn
        con(ii) = EXP((-(ii-9) * (ii-9)/2*sigma*sigma))
        vel(ii) = 150
        area(ii) = 1.
        areav(ii) = 1.
        flxarr(ii) = 0.
    end do
    
    do ii = 1,100000
        call hadvppm(nn,dt,dx,con,vel,area,areav,flxarr,mynn)
    end do
    
end program UnitTestHadvppm
