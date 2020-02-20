subroutine checkboundary(DEM, nr, nc, boundary)
    ! This is the code to identify boundaries of the DEM, created by Kwanghun
    ! Choi (2015).

    ! use external module to use nan.
    !use, intrinsic::iso_fortran_env
    !use, intrinsic::ieee_arithmetic
    implicit none
    ! The "nr", and the "nc" are the number of row and column of the DEM.
    integer, intent( in ) :: nr, nc
    ! The "DEM" is the original DEM before sink filling.
    double precision, dimension( nr, nc ) :: DEM
    ! The "buffer" is the DEM with buffer to identify the boundary of the DEM.
    double precision, dimension( 0 : (nr + 1), 0 : (nc + 1) ) :: buffer
    double precision, dimension( 3, 3 ) :: buffer_m3
    ! The "boundary" is the boundary of the DEM"
    double precision, dimension( nr, nc ) :: boundary
    ! Dummy variables and the minimum and maximum row and column number.
    integer :: i, j
    ! Parameters for the NaN value.
    double precision :: zero, NaN
    !    NaN = transfer(z'7ff8000000000000', 1.0d0)
    
    zero = 0.d0
    NaN = 0.d0 / zero
    !NaN = ieee_value(NaN, ieee_quiet_nan)

    where(DEM .lt. -99999) DEM = NaN
    buffer = NaN; buffer(1:nr, 1:nc) = DEM
    boundary = NaN
    ! If adjacent of target cell has at least one NaN and one valid value
    ! (i.e., not isolated), we define the cell as boundary.
    ! That condition is same as number of adjacent cells not having zero or eight NaNs.
    do i = 1, nr, 1
        do j = 1, nc, 1
            buffer_m3 = buffer( i-1:i+1, j-1:j+1 )
            !if ( ( .not. ieee_is_nan( buffer( i, j ) ) ) .and. &
            !    ( count( ieee_is_nan( buffer( i-1:i+1, j-1:j+1 ) ) ) .gt. 0 ) ) then
            if ( ( buffer( i, j ) .eq. buffer( i, j ) ) .and. &
                ( count( buffer_m3 .ne. buffer_m3 ) .gt. 0 ) ) then
                boundary(i, j) = 1.0d0
            else
                boundary(i, j) = NaN
            end if
        end do
    end do
end subroutine checkboundary



