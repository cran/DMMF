subroutine sinkfill(DEM, nr, nc, res, boundary, min_angle, DEM_nosink, partition)
    ! This is the Fortran code to produce sink free DEM using an algorithm
    ! proposed by Wang & Liu (2006), and code is created by Kwanghun Choi (2015). 
    ! References
    ! [1] Wang, L. and Liu, H. (2006). An efficient method for identifying and
    ! filling surface depressions in digital elevation models for hydrologic
    ! analysis and modelling. International Journal of Geographical Information
    ! Science, 20(2):193–213.
    ! [2] Volker Wichmann (2007) Module Fill Sinks (Wang & Liu). SAGA-GIS Module 
    ! Library Documentation (v2.1.3) [ cited 2015. 08. 15 ], Available from:
    ! http://www.saga-gis.org/saga_module_doc/2.1.3/ta_preprocessor_4.html
    
    ! use external module to use nan.
    !use, intrinsic::iso_fortran_env
    !use, intrinsic::ieee_arithmetic
    implicit none
    ! The "nr", and the "nc" are the number of row and column of the DEM.
    integer, intent(in) :: nr, nc
    ! The "res" is the reolution of the DEM.
    ! The "min_angle" is the minimum angle between the cell and the adjacent
    ! cells.
    double precision, intent(in) :: res, min_angle
    ! The "DEM" is the original DEM before sink filling.
    ! The "DEM_t" is the temporary DEM for sink filling algorithm.
    ! The "boundary" is the boundary of the DEM"
    ! The "DEM_nosink" is the sink filled DEM.
    double precision, dimension( nr, nc ) :: DEM, DEM_t, boundary, DEM_nosink, partition
    ! min_diff matrix
    double precision, dimension( -1:1, -1:1 ) :: min_diff, min_block
    double precision, dimension( -1:1, -1:1 ) :: bnd_block, prt_block, DEM_block, out_block
    ! Dummy variables and the minimum and maximum row and column number.
    integer :: r, c, mnr, mxr, mnc, mxc, t
    double precision :: mda, mdd
    ! Variables of locations for maximum value and minumum value.
    integer, dimension(2) :: min_loc
    ! Parameters for the NaN value.
    double precision :: zero, NaN
    !    NaN = transfer(z'7ff8000000000000', 1.0d0)

    zero = 0.d0
    NaN = 0.d0 / zero
    !NaN = ieee_value(NaN, ieee_quiet_nan)

    where( DEM .lt. -99999 ) DEM = NaN
    DEM_t = DEM; DEM_nosink = DEM
    partition = NaN; t = 0
    where( boundary .lt. -99999 ) boundary = NaN

    mda = res * dtan(min_angle)               ! min difference for adjoining cells
    mdd = res * dtan(min_angle) * sqrt(2.0d0) ! min difference at diagonal cells  
    min_diff = reshape((/ mdd, mda, mdd, mda, 0.0d0, mda, mdd, mda, mdd /), shape(min_diff)) 

    do while( any( boundary .eq. boundary ) )

        min_loc = minloc( DEM_t, mask = ( boundary .eq. boundary ) )
        r = min_loc(1)
        c = min_loc(2)
        ! Assing partition if the target cell doesn't have it.
        if( partition( r, c ) .ne. partition( r, c ) ) then
            t = t + 1
            partition( r, c ) = t
        end if
        ! Check the location of target cell and its surrounding cells.
        mnr = max0( 1, r - 1 )
        mxr = min0( nr, r + 1 )
        mnc = max0( 1, c - 1 )
        mxc = min0( nc, c + 1 )
        ! Initialized blocks and assign values from original maps.
        ! Initialization
        min_block = NaN; bnd_block = NaN; prt_block = NaN; DEM_block = NaN
        ! Assign values
        bnd_block( (mnr-r):(mxr-r), (mnc-c):(mxc-c) ) = boundary( mnr:mxr, mnc:mxc )
        prt_block( (mnr-r):(mxr-r), (mnc-c):(mxc-c) ) = partition( mnr:mxr, mnc:mxc )
        DEM_block( (mnr-r):(mxr-r), (mnc-c):(mxc-c) ) = DEM_t( mnr:mxr, mnc:mxc )
        out_block( (mnr-r):(mxr-r), (mnc-c):(mxc-c) ) = DEM_nosink( mnr:mxr, mnc:mxc )
        min_block( (mnr-r):(mxr-r), (mnc-c):(mxc-c) ) = &
            DEM_t( r, c ) + min_diff( (mnr-r):(mxr-r), (mnc-c):(mxc-c) )
        ! Assign new boundary and perform sink fill algorithm
        !where( ( ( .not. ieee_is_nan(DEM_block ) ) .and. ieee_is_nan( bnd_block ) )&
        !        .and. (DEM_block .le. min_block) ) DEM_block = min_block
        where( (DEM_block .eq. DEM_block) .and. (DEM_block .le. min_block)&
                .and. (bnd_block .ne. bnd_block) ) DEM_block = min_block
        !where( ( ieee_is_nan( prt_block ) .and. ieee_is_nan( bnd_block) ) &
        !        .and. (.not. ieee_is_nan( out_block ) ) ) prt_block = partition( r, c )
        where( (prt_block .ne. prt_block) .and. (out_block .eq. out_block)&
                .and. (bnd_block .ne. bnd_block) ) prt_block = partition( r, c )
        !where( .not. ieee_is_nan( DEM_block ) )
        where( DEM_block .eq. DEM_block ) 
            bnd_block = 1.0d0
            out_block = DEM_block
        end where

        DEM_t( mnr:mxr, mnc:mxc )       = DEM_block( (mnr-r):(mxr-r), (mnc-c):(mxc-c) )
        DEM_nosink( mnr:mxr, mnc:mxc )  = out_block( (mnr-r):(mxr-r), (mnc-c):(mxc-c) )
        DEM_t( r, c ) = NaN
        boundary( mnr:mxr, mnc:mxc )    = bnd_block( (mnr-r):(mxr-r), (mnc-c):(mxc-c) )
        boundary( r, c ) = NaN
        partition( mnr:mxr, mnc:mxc )   = prt_block( (mnr-r):(mxr-r), (mnc-c):(mxc-c) )
    end do
end subroutine sinkfill

