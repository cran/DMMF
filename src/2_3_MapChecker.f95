subroutine mapchecker(DEM, nr, nc, boundary, sink, stand)

    ! use external module to use nan.
    !use, intrinsic::iso_fortran_env
    !use, intrinsic::ieee_arithmetic
    implicit none
    ! The "nr", and the "nc" are the number of row and column of the DEM.
    integer, intent(in) :: nr, nc
    ! The "DEM" is the original DEM before sink filling.
    ! The "DEM_t" is the temporary DEM for sink filling algorithm.
    ! The "sink" is the sink of the DEM"
    ! The "DEM_nosink" is the sink filled DEM.
    double precision, dimension( nr, nc ) :: DEM
    double precision, dimension( nr, nc ) :: boundary, sink, stand
    double precision, dimension( 0:(nr+1), 0:(nc+1) ) :: DEM_b
    ! Dummy variables and the minimum and maximum row and column number.
    integer :: i, j
    ! Variables of locations for maximum value and minumum value.
    !integer, dimension(2) :: min_loc
    ! Parameters for the NaN value.
    double precision :: zero, NaN
    !    NaN = transfer(z'7ff8000000000000', 1.0d0)
    zero = 0.d0
    NaN = 0.d0 / zero
    !NaN = ieee_value(NaN, ieee_quiet_nan)

    ! Define NaN values
    where( DEM .lt. -99999 ) DEM = NaN
    ! Initialize
    DEM_b = NaN; boundary = NaN; sink = NaN; stand = NaN
    ! Create DEM with buffer cells with NaN
    DEM_b(1:nr, 1:nc) = DEM
    
    !forall (i=1:nr, j=1:nc, ( .not. ieee_is_nan( DEM_b( i, j ) ) ) )
    forall (i=1:nr, j=1:nc, DEM_b( i, j ) .eq. DEM_b( i, j ) )
        boundary(i, j) =&
            dble( count( DEM_b( i-1:i+1, j-1:j+1 ) .ne.&
            DEM_b( i-1:i+1, j-1:j+1 ) ) )
        stand(i, j) =&
            dble( count( DEM_b(i-1:i+1, j-1:j+1) .lt. DEM(i, j) ) )
        sink(i, j) =&
            dble( count( DEM_b(i-1:i+1, j-1:j+1) .le. DEM(i, j) )  - 1 )
    end forall

    where(boundary .gt. 0.0d0) boundary = 1.0d0
    where(boundary .eq. 0.0d0) boundary = NaN
    where(sink .gt. 0.0d0) sink = NaN
    where(stand .gt. 0.0d0) stand = NaN
    where(sink .eq. 0.0d0) sink = 1.0d0
    where(stand .eq. 0.0d0) stand = 1.0d0

end subroutine mapchecker
