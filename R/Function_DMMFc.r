DMMFc <-
    function(DEM, R, RI, ET, P_c, P_z, P_s, theta_init, theta_sat, theta_fc, 
             SD, K, P_I, n_s, CC, GC, IMP, PH, D, NV, d_a = 0.005, 
             DK_c = 0.1, DK_z = 0.5, DK_s = 0.3, 
             DR_c = 1.0, DR_z = 1.6, DR_s = 1.5,
             SoilMap, LULCMap,
             Breaking, Init_point, Sinks, R_Type = 0, slpMode = 2, ALL = TRUE)
    {
        # Input variable check
        # 1. Raster maps  
        if ( is( DEM, "RasterLayer" ) & ( all.equal(xres( DEM ), yres( DEM ) ) ) ){
            # Extract basic infomation of a DEM
            res     <- xres( DEM )
            xmn     <- extent(DEM)[1]
            xmx     <- extent(DEM)[2]
            ymn     <- extent(DEM)[3]
            ymx     <- extent(DEM)[4]
            crs     <- proj4string(DEM)

            # Calculate the valid cells in DEM (Cells that does not have "NA").
            vc <- cellStats(!is.na(DEM), sum)

            # Convert format of a DEM from raster to matrix for FORTRAN subroutine.  
            DEM[is.na(DEM)] <- -999999; DEM_m   <- as.matrix( DEM )
            mask_r <- DEM * 0; mask_m <- as.matrix( mask_r )
        }else{
            stop("Missing or invalid type of DEM")
        }
        
        # Convert format of a Sinks from raster to matrix for FORTRAN subroutine.  
        if( missing( Sinks ) )
            {
                Sinks_m <- mask_m - 999999
            }else{
                vc <- vc - cellStats(!is.na(Sinks), sum)
                Sinks[is.na(Sinks)] <- -999999 
                Sinks_m   <- as.matrix( Sinks )
            }

        # Convert format of a SoilMap from raster to matrix for FORTRAN subroutine.  
        if( missing( SoilMap ) )
            {
                stop("Missing or invalid type of soil map")
            }else{
                SoilMap[is.na(SoilMap)] <- -999999 
                SoilMap_m   <- as.matrix( SoilMap )
            }

        # Convert format of a LULCMap from raster to matrix for FORTRAN subroutine.  
        if( missing( LULCMap ) )
            {
                stop("Missing or invalid type of LULC map")
            }else{
                LULCMap[is.na(LULCMap)] <- -999999 
                LULCMap_m   <- as.matrix( LULCMap )
            }

        # 2. Vectors
        # Meteorogical data
        # Calculate the total days of interest (event or events)
        days <- nlayers( R )

        # 3. R, RI, and ET Arrays
        R_m <- as.array( R )
        RI_m <- as.array( RI )
        ET_m <- as.array( ET )

        # 4. Breaking points and initial soil water content

        # If the model run with only one period.
        if( missing( Breaking ) )
        {
            warning( "Breaking point is not specified, consider the case as one event", immediate = TRUE )
            Breaking <- c( rep( 0, (days - 1) ), 1 )
        }else{
            Breaking[ Breaking > 0 ] <- 1 
        }
        # Number of results of days.
        N_out <- as.integer( sum(Breaking) )

        # The number of output layers which is same for Breaking point.
        mask_out <- as.array( brick( replicate( N_out, mask_r ) ) )
        # Because the program can calculate with initial soil water content (Theta_init) for every events)
        # 8. Initial soil water content and point
        # If the model run with only one initial soil water condition.
        if( missing( Init_point ) )
        {
            warning( "No Init_point is specified, model assumes there is only one measured initial soil water content at the first day of simulation", immediate = TRUE )
            Init_point <- c( 1, rep( 0, (days - 1) ) )
        }else{
            Init_point[ Init_point > 0 ] <- 1 
        }
        Init_days <- sum(Init_point)


        MMF_result <- .Fortran( C_dmmfc, DEM = DEM_m, SoilMap = SoilMap_m, LULCMap = LULCMap_m, 
                               nr = nrow(DEM), nc = ncol(DEM), res = res, 
                               option = as.integer(slpMode), days = as.integer(days), 
                               R = R_m, RI = RI_m, R_Type = as.integer(R_Type), 
                               ET = ET_m, P_c = P_c, P_z = P_z, P_s = P_s, 
                               theta_init = unname(as.matrix(theta_init)), theta_sat = theta_sat, 
                               theta_fc = theta_fc, SD = SD, K = K, 
                               P_I = unname(as.matrix(P_I)), n_s = unname(as.matrix(n_s)), 
                               d_a = unname(as.matrix(d_a)), CC = unname(as.matrix(CC)), 
                               GC = unname(as.matrix(GC)), IMP = unname(as.matrix(IMP)), 
                               PH = unname(as.matrix(PH)), D = unname(as.matrix(D)), 
                               NV = unname(as.matrix(NV)), 
                               DK_c = DK_c, DK_z = DK_z, DK_s = DK_s, 
                               DR_c = DR_c, DR_z = DR_z, DR_s = DR_s, 
                               Breaking = as.integer(Breaking), 
                               N_out = as.integer(N_out), 
                               vc = as.integer(vc), Init_days = as.integer(Init_days),
                               Init_point = as.integer(Init_point), sinks = Sinks_m,
                               A = mask_m, Rf_r = mask_out, SW_c_r = mask_out, 
                               theta_r_r = mask_out, TC_r = mask_out,
                               Q_in_r = mask_out, Q_out_r = mask_out, 
                               IF_in_r = mask_out, IF_out_r = mask_out, 
                               SS_c_r = mask_out, SS_z_r = mask_out, SS_s_r = mask_out,
                               G_c_r = mask_out, G_z_r = mask_out, G_s_r = mask_out,
                               SL_c_in_r = mask_out, SL_z_in_r = mask_out, 
                               SL_s_in_r = mask_out, SL_in_r = mask_out,
                               SL_c_out_r = mask_out, SL_z_out_r = mask_out, 
                               SL_s_out_r = mask_out, SL_out_r = mask_out, 
                               NAOK = T)

        # Convert result to visible raster brick
        if(ALL == FALSE){
            # 1. Area of each element (cell)
            A         <- raster(MMF_result$A, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            # 2. Inflow of surface runoff from upslope elements (cells)
            Q_in        <- brick(MMF_result$Q_in_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(Q_in) <- paste("Q_in", c(1:length(names(Q_in))), sep="_") 
            # 3. outflow of surface runoff from elements (cells)
            Q_out        <- brick(MMF_result$Q_out_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(Q_out) <- paste("Q_out", c(1:length(names(Q_out))), sep="_") 
            # 4. Inflow of interflow from upslope elements (cells)
            IF_in        <- brick(MMF_result$IF_in_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(IF_in) <- paste("IF_in", c(1:length(names(IF_in))), sep="_") 
            # 5. outflow of interflow from elements (cells)
            IF_out        <- brick(MMF_result$IF_out_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(IF_out) <- paste("IF_out", c(1:length(names(IF_out))), sep="_") 
            # 6. inflow of soil particles from elements (cells)
            SL_in        <- brick(MMF_result$SL_in_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(SL_in) <- paste("SL_in", c(1:length(names(SL_in))), sep="_") 
            # 7. outflow of soil particles from elements (cells)
            SL_out        <- brick(MMF_result$SL_out_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(SL_out) <- paste("SL_out", c(1:length(names(SL_out))), sep="_") 
            # 8. Outputs
            output <- list(A, Q_out, Q_in, IF_out, IF_in, SL_out, SL_in) 
            names(output) <- c( "A", "Q_out", "Q_in", "IF_out", "IF_in", "SL_out", "SL_in")
        }else{
            # 1. Area of each element (cell)
            A         <- raster(MMF_result$A, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            # 2. Effective rainfall
            Rf        <- brick(MMF_result$Rf_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(Rf) <- paste("Rf", c(1:length(names(Rf))), sep="_") 
            # 3. Soil water infiltration capacity  
            SW_c        <- brick(MMF_result$SW_c_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(SW_c) <- paste("SW_c", c(1:length(names(SW_c))), sep="_") 
            # 4. Simulated initial soil water content. 
            theta_r        <- brick(MMF_result$theta_r_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(theta_r) <- paste("theta_r", c(1:length(names(theta_r))), sep="_") 
            # 5. Transport capacity of surface runoff
            TC        <- brick(MMF_result$TC_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(TC) <- paste("TC", c(1:length(names(TC))), sep="_") 
            # 6. Inflow of surface runoff from upslope elements (cells)
            Q_in        <- brick(MMF_result$Q_in_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(Q_in) <- paste("Q_in", c(1:length(names(Q_in))), sep="_") 
            # 7. outflow of surface runoff from elements (cells)
            Q_out        <- brick(MMF_result$Q_out_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(Q_out) <- paste("Q_out", c(1:length(names(Q_out))), sep="_") 
            # 8. Inflow of interflow from upslope elements (cells)
            IF_in        <- brick(MMF_result$IF_in_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(IF_in) <- paste("IF_in", c(1:length(names(IF_in))), sep="_") 
            # 9. outflow of interflow from elements (cells)
            IF_out        <- brick(MMF_result$IF_out_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(IF_out) <- paste("IF_out", c(1:length(names(IF_out))), sep="_") 
            # 10. Suspended sediments of each particle size class
            SS_c        <- brick(MMF_result$SS_c_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(SS_c) <- paste("SS_c", c(1:length(names(SS_c))), sep="_") 
            SS_z        <- brick(MMF_result$SS_z_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(SS_z) <- paste("SS_z", c(1:length(names(SS_z))), sep="_") 
            SS_s        <- brick(MMF_result$SS_s_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(SS_s) <- paste("SS_s", c(1:length(names(SS_s))), sep="_") 
            # 11. available sediments for transport of each particle size class
            G_c        <- brick(MMF_result$G_c_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(G_c) <- paste("G_c", c(1:length(names(G_c))), sep="_") 
            G_z        <- brick(MMF_result$G_z_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(G_z) <- paste("G_z", c(1:length(names(G_z))), sep="_") 
            G_s        <- brick(MMF_result$G_s_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(G_s) <- paste("G_s", c(1:length(names(G_s))), sep="_") 
            # 12. sediment inputs from upslope elements (cells) of each particle size class
            SL_c_in        <- brick(MMF_result$SL_c_in_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(SL_c_in) <- paste("SL_c_in", c(1:length(names(SL_c_in))), sep="_") 
            SL_z_in        <- brick(MMF_result$SL_z_in_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(SL_z_in) <- paste("SL_z_in", c(1:length(names(SL_z_in))), sep="_") 
            SL_s_in        <- brick(MMF_result$SL_s_in_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(SL_s_in) <- paste("SL_s_in", c(1:length(names(SL_s_in))), sep="_") 
            SL_in        <- brick(MMF_result$SL_in_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(SL_in) <- paste("SL_in", c(1:length(names(SL_in))), sep="_") 
            # 13. sediment outputs from elements (cells) of each particle size class
            SL_c_out        <- brick(MMF_result$SL_c_out_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(SL_c_out) <- paste("SL_c_out", c(1:length(names(SL_c_out))), sep="_") 
            SL_z_out        <- brick(MMF_result$SL_z_out_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(SL_z_out) <- paste("SL_z_out", c(1:length(names(SL_z_out))), sep="_") 
            SL_s_out        <- brick(MMF_result$SL_s_out_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(SL_s_out) <- paste("SL_s_out", c(1:length(names(SL_s_out))), sep="_") 
            SL_out        <- brick(MMF_result$SL_out_r, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx, crs = crs)
            names(SL_out) <- paste("SL_out", c(1:length(names(SL_out))), sep="_") 
            # 14. Outputs
            output <- list( A, Rf, SW_c, theta_r, TC, Q_in, Q_out, IF_in, IF_out, SS_c, SS_z, SS_s, G_c, G_z, G_s, SL_c_in, SL_z_in, SL_s_in, SL_in, SL_c_out, SL_z_out, SL_s_out, SL_out)
            names(output) <- c( "A", "Rf", "SW_c", "theta_r", "TC", "Q_in", "Q_out", "IF_in", "IF_out", "SS_c", "SS_z", "SS_s", "G_c", "G_z", "G_s", "SL_c_in", "SL_z_in", "SL_s_in", "SL_in", "SL_c_out", "SL_z_out", "SL_s_out", "SL_out")
        }
        return(output)
    }
