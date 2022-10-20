!> \file mo_wqm_readcatchdata.f90

!> \brief Reads water quality configure, input files and allocating WQ global variables.

!> \details This module is for water quality modeling.\n
!> 

!> \authors Xiaoqiang Yang
!> \date Sep 2021

MODULE mo_wqm_readcatchdata

  USE mo_kind, ONLY: i4, dp

  
  IMPLICIT NONE

  PUBLIC :: wqm_readcatchdata    

  contains
  ! 
  ! ------------------------------------------------------------------

  !     NAME
  !         wqm_readcatchdata

  !     PURPOSE
  !>        \brief Reads water quality information from input files.

  !>        \details Reads water quality information from input files.  \n
  !>        Four input files are located in "/input/water_quality/".

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Xiaoqiang Yang
  !>        \date Sep 2021



  subroutine wqm_readcatchdata()

    use mo_read_spatial_data,   only: read_header_ascii,     &
                                      read_spatial_data_ascii
    use mo_append,              only: append
    use mo_message,             only: message
    use mo_string_utils,        only: num2str
    use mo_wqm_global_variables, only: &
         dir_catchinfo, fn_catchDomain, model_type, m_Level,resolution, &
         r_Level, rout_resolution,Rout_mask, simPer, &
         Model_mask,L1_nCells, L1_area,  L1_nLinks,       & ! area in square-meter
         dir_gauges, fn_gauge_loc, basin_gauge, L1_gaugeLoc,    &
         nHydroGauges,nWQGauges, nInflowGauges, & 
         fdir_type, fn_flowDir, L1_fDir, fn_spatialorganization, &
         L1_fromN,L1_toN, OutNode, L1_rOrder, L1_netPerm, &
         fn_channel_mask, fn_channel_length, fn_channel_width,   &
         L1_chanmask,L1_length,L1_width, R_areaCell, &
         geoProperty, &
         nLanduseTotal, nLUperGrid, fn_LanduseTypes,fn_flanduse,L1_fLAI, L1_fSealed,  &
         Sealed_LanduseNrs, Sealed_LanduseIDs, &
         nSoilTotal, nSoilperGrid, fn_SoilTypes, fn_fsoil, L1_fsoil, &
         nRotationTotal, nRotationperGrid, fn_RotationTypes, fn_frotation, L1_frotation,  &
         nSoilLayers, fn_soilDepth, L1_soildepth, &
         !semi-distributed
         L1_nSubs, L1_nHRUs, SubID, SubArea,SubHRUs, CropRotationID, &
         fn_classProperties, &
         fn_wiltpoint, fn_satMoist, L1_wiltingpoint, L1_saturatedSM,    &  
         fn_initial_condition,fn_cropinfo,fn_lut_rotation,fn_lut_monLAI, fn_globalRadiation, &
         fn_init_concIN,fn_init_concON,  & 
         nCroptation, &
         num_crops,       & !number of crop types in total
         rotation,        & !rotation infos (type variables)
         cropdata,        & !crop infos (type variables)
         init_type,    &
         init_concIN,  & !
         init_concON,  & !
         init_humusN,  & !
         init_fastN,   & !
         hnhalf,       & !
         LAIUnitList, & !
         LAILUT,      & !
         GR_file_exist, & !
         global_radiation, &  !
         norlai_daily,     & !
         nor_globalradi  
    use mo_wqm_constants,       only: nodata_i4,nodata_dp, maxnCropsInRotation
    use mo_julian,              only: julday
    use mo_wqm_shadingeffect,    only: read_shading_data   !read and normalise data for shading calculation
    use mo_spatial_organize,     only: set_topologic_network, get_connection_sequences,match_RtoM

    implicit none
    !---
    !local
    integer(i4)      :: i,j,ii,jj
    character(256)   :: fName, dummy, line
    integer(i4)      :: funit  
    integer(i4), dimension(:,:), allocatable   :: data_i4_2d
    real(dp), dimension(:,:), allocatable      :: data_dp_2d
    logical, dimension(:,:), allocatable       :: mask_2d, mask_global, mask_global_r
    logical, dimension(:), allocatable         :: mask_tmp
    logical                                    :: file_exist
    !integer(i4)      :: nLinks
    integer(i4), dimension(:), allocatable   :: fromNode,toNode,Rsequence,netPerm
    real(dp),    dimension(:), allocatable   :: slength,swidth
    integer(i4), dimension(:), allocatable   :: tmp_ids
    real(dp),    dimension(:), allocatable   :: tmp_areaS
    real(dp),    dimension(:), allocatable :: file_initINconc, file_initONconc

    integer(i4)      ::  nClasses  !for initial value read-in
	
    !variables for global radiation
    !## added by yangx Nov 2016 ##
     
    integer(i4), dimension(3)    :: start_period, end_period, start_period_file, end_period_file !calendar date
    integer(i4), dimension(3)    :: file_time
    real(dp)                     :: nodata_file      ! no data value of data
    integer(i4)                  :: start_jul,end_jul,start_jul_file,end_jul_file !julian date  
    real(dp), dimension(:), allocatable   :: gr_file  ! raw global radiation data in the file
    logical, dimension(:),  allocatable   :: gr_mask   ! mask where measured data are missing in the raw data file
    real(dp), dimension(:), allocatable :: gr_data     !dim1=1, dim2= period length
    real(dp)                     :: avg_gr  !average gr value through all measured data, also as a default value
                                            !for global variable (global_radiation)  	
    integer(i4)                  :: length_file, length_period ! number of days in file and period respectivly
    integer(i4)                  :: idx_st_period, idx_en_period, idx_st_file, idx_en_file
                                    !index for merge file and period
    integer(i4)                  :: nLAI,idgauges
    real(dp)                     :: cellFactorR


    if (model_type == 1_i4) then
    !*****************************************************************
	!READ SPATIAL DATA OF DISTRIBUTED CATCHMENT INFORMATION
	!*****************************************************************	
    !******************************************************************************************************
    !check headers
        funit = 111
        fName = trim(adjustl(dir_catchinfo)) // trim(adjustl(fn_catchDomain))
        call read_header_ascii(trim(fName), funit,   &
            m_Level%ncols,     m_Level%nrows, m_Level%xllcorner, &
            m_Level%yllcorner, m_Level%cellsize, m_Level%nodata_value)

        ! check for L0 and L1 scale consistency
        if (abs(resolution - m_Level%cellsize) > 0.0001) then
          call message()
          call message('***ERROR: model resolution should be matched with the input data resolution')
          call message('check wqm_config.nml: *resolution* and the file *fn_catchDomain*')
          stop
        end if
    !**********************************************************************************************
    !modeling domain
        fName = trim(adjustl(dir_catchinfo)) // trim(adjustl(fn_catchDomain))
        call read_spatial_data_ascii(trim(fName), funit, &
            m_Level%ncols,     m_Level%nrows, m_Level%xllcorner,&
            m_Level%yllcorner, m_Level%cellsize, data_dp_2d, mask_global)

        call message('    Reading catchment domain: ', trim(fn_catchDomain),' ...')
       !
       ! create overall model mask and save indices
        allocate(Model_mask(size(mask_global, dim=1), size(mask_global, dim=2)))
        Model_mask = mask_global
        !number of cells within the domain
        L1_nCells = count(mask_global)
        !the data
        data_dp_2d = data_dp_2d * real(resolution,dp) * real(resolution,dp)  ![m2]
        data_dp_2d = merge(data_dp_2d, nodata_dp, mask_global)
        call append( L1_area, pack(data_dp_2d, mask_global) )
       ! deallocate arrays
        deallocate(data_dp_2d)
    !**********************************************************************************************
    !flow direction for grid connection
    !check if modeling and routing resolutions are the same or not
        if (abs(rout_resolution-resolution) <0.001_dp ) then
            r_Level = m_Level
            if (.not. allocated(R_areaCell)) allocate(R_areaCell(size(L1_area)))
            R_areaCell = L1_area
        else
            cellFactorR = rout_resolution/resolution
            if ( nint(cellFactorR,i4) < 1_i4) then
                call message()
                call message('***ERROR: rout_resolution should be divided by resolution with no remainder!')
                stop
            end if
            fName = trim(adjustl(dir_catchinfo)) // trim(adjustl(fn_flowDir) )
            funit = 112
            !read headers
            call read_header_ascii(trim(fName), funit,   &
            r_Level%ncols,     r_Level%nrows, r_Level%xllcorner, &
            r_Level%yllcorner, r_Level%cellsize, r_Level%nodata_value)
            if (abs(rout_resolution - r_Level%cellsize) > 0.0001) then
              call message()
              call message('***ERROR: routing resolution should be matched with the input data resolution')
              call message('check wqm_config.nml: *rout_resolution* and the file *fn_flowDir*')
              stop
            end if
        end if
        fName = trim(adjustl(dir_catchinfo)) // trim(adjustl(fn_flowDir) )
        funit = 112
        call read_spatial_data_ascii(trim(fName), funit, &
            r_Level%ncols,     r_Level%nrows, r_Level%xllcorner,&
            r_Level%yllcorner, r_Level%cellsize, data_i4_2d, mask_global_r)
        call message('    Reading flow direction map...')
       ! create overall model mask and save indices
        allocate(Rout_mask(size(mask_global_r, dim=1), size(mask_global_r, dim=2)))
        Rout_mask = mask_global_r
        !number of links within the domain
        L1_nLinks = count(mask_global_r)
        data_i4_2d = merge(data_i4_2d, nodata_i4, mask_global_r )
        call append( L1_fDir, pack(data_i4_2d, mask_global_r) )



		!calculate the connection seqence
        call set_topologic_network(fdir_type)  !get from and to node info
        call get_connection_sequences()  !get the routing sequences
        deallocate(data_i4_2d) 
		
        !following mHM-mRM implementations to match the modeling and routing resolutions
        !implementation not finished yet!!!
        if (rout_resolution > resolution) then
            call match_RtoM() !match the routing grid ids and modeling grid ids
        end if
    !**********************************************************************************************
    !channel mask + length and width information
        fName = trim(adjustl(dir_catchinfo)) // trim(adjustl(fn_channel_mask))
        funit = 113
        call read_spatial_data_ascii(trim(fName), funit, &
            r_Level%ncols,     r_Level%nrows, r_Level%xllcorner,&
            r_Level%yllcorner, r_Level%cellsize, data_i4_2d, mask_2d)

        call message('    Reading channel information ...')
       !check if all channel mask is within the model domain
        if (any(mask_2d .and. (.not. mask_global_r))) then
           call message()
           call message('***ERROR: channel outside of the routing domain')
           stop
        end if
        call append( L1_chanmask, pack(mask_2d, mask_global_r) )
       !deallocate arrays
        deallocate(data_i4_2d, mask_2d)
     !channel length
        fName = trim(adjustl(dir_catchinfo)) // trim(adjustl(fn_channel_length))
        funit = 114
        call read_spatial_data_ascii(trim(fName), funit, &
            r_Level%ncols,     r_Level%nrows, r_Level%xllcorner,&
            r_Level%yllcorner, r_Level%cellsize, data_dp_2d, mask_2d)

        data_dp_2d = merge(data_dp_2d, nodata_dp, mask_global_r)
        call append( L1_length, pack(data_dp_2d, mask_global_r) )
       ! deallocate arrays
        deallocate(data_dp_2d,mask_2d)
     !channel width
        fName = trim(adjustl(dir_catchinfo)) // trim(adjustl(fn_channel_width))
        funit = 115
        INQUIRE(file = fName, exist = file_exist)
        if (file_exist) then
          call read_spatial_data_ascii(trim(fName), funit, &
            r_Level%ncols,     r_Level%nrows, r_Level%xllcorner,&
            r_Level%yllcorner, r_Level%cellsize, data_dp_2d, mask_2d)

          data_dp_2d = merge(data_dp_2d, nodata_dp, mask_global_r)
          call append( L1_width, pack(data_dp_2d, mask_global_r) )
         ! deallocate arrays
          deallocate(data_dp_2d,mask_2d)
        else
          allocate(data_dp_2d(size(mask_global_r, dim=1), size(mask_global_r, dim=2)))
          data_dp_2d= nodata_dp
          call append( L1_width, pack(data_dp_2d, mask_global_r) )
          call message('***Attention: channel width is not provided, but calculated using hard-coded equation')
        end if
    !**********************************************************************************************       
    !INITIAL IN/ON concentration maps
        fName = trim(adjustl(dir_catchinfo)) // trim(adjustl(fn_init_concIN))
        funit = 131
        call read_spatial_data_ascii(trim(fName), funit, &
            m_Level%ncols,     m_Level%nrows, m_Level%xllcorner,&
            m_Level%yllcorner, m_Level%cellsize, data_dp_2d, mask_2d)
       !check if all channel mask is within the model domain
        if (.not. all(mask_2d .eqv. mask_global)) then
           call message()
           call message('***ERROR: IN initial value map dose not match the catchment domain')
           stop
        end if
        data_dp_2d = merge(data_dp_2d, nodata_dp, mask_global)
        call append(init_concIN, pack(data_dp_2d, mask_global) )
       ! deallocate arrays
        deallocate(data_dp_2d,mask_2d)

        fName = trim(adjustl(dir_catchinfo)) // trim(adjustl(fn_init_concON))
        funit = 132
        call read_spatial_data_ascii(trim(fName), funit, &
            m_Level%ncols,     m_Level%nrows, m_Level%xllcorner,&
            m_Level%yllcorner, m_Level%cellsize, data_dp_2d, mask_2d)
       !check if all channel mask is within the model domain
        if (.not. all(mask_2d .eqv. mask_global)) then
           call message()
           call message('***ERROR: ON initial value map dose not match the catchment domain')
           stop
        end if
        data_dp_2d = merge(data_dp_2d, nodata_dp, mask_global)
        call append(init_concON, pack(data_dp_2d, mask_global) )
       ! deallocate arrays
        deallocate(data_dp_2d,mask_2d)

    !********************************************************************************************** 
    !LAND USE map and areal share of each grid (if nLUperGrid > 1)
        !allocate
        allocate(L1_fLAI(L1_nCells,nLanduseTotal))
        L1_fLAI = 0.0_dp
        if (nLUperGrid == 1_i4) then
          !read
          fName = trim(adjustl(dir_catchinfo)) // trim(adjustl(fn_LanduseTypes))
          funit = 116
          INQUIRE(file = fName, exist = file_exist)
          if (.not. file_exist) then
            call message()
            call message('***ERROR: since nLUperGrid =1, land use map (fn_LanduseTypes) should be provided!')
            stop
          end if
          call read_spatial_data_ascii(trim(fName), funit, &
            m_Level%ncols,     m_Level%nrows, m_Level%xllcorner,&
            m_Level%yllcorner, m_Level%cellsize, data_i4_2d, mask_2d)

          call message('    Reading land-use information(fn_LanduseTypes)...')       
           !check if mask matches
          if (.not. all(mask_2d .eqv. mask_global)) then
             call message()
             call message('***ERROR: land use map dose not match catchment domain')
          stop
          end if

          data_i4_2d = merge(data_i4_2d, nodata_i4, mask_2d )
          allocate(mask_tmp(L1_nCells))
          call append(tmp_ids, pack(data_i4_2d, mask_2d))
          do ii = 1, nLanduseTotal
            mask_tmp(:) = .FALSE.
            where(tmp_ids == ii) mask_tmp = .TRUE.
            L1_fLAI(:,ii) = merge(1.0_dp,0.0_dp,mask_tmp)

          end do
          deallocate(data_i4_2d,mask_2d,mask_tmp,tmp_ids)
        else
		  !read maps of areal share for each grid
          if (size(fn_flanduse) .ne. nLUperGrid) then
            call message()
            call message('***ERROR: since nLUperGrid > 1, each areal share maps should be provided!')
            stop
          end if
         call message('    Reading land-use information(fn_flanduse) ...') 
          do ii = 1, nLUperGrid
            fName = trim(adjustl(dir_catchinfo)) // trim(adjustl(fn_flanduse(ii)))
            funit = 126
            call read_spatial_data_ascii(trim(fName), funit, &
              m_Level%ncols,     m_Level%nrows, m_Level%xllcorner,&
              m_Level%yllcorner, m_Level%cellsize, data_dp_2d, mask_2d)
      
           !check if mask matches
            if (.not. all(mask_2d .eqv. mask_global)) then
             call message()
             call message('***ERROR: land use map dose not match catchment domain')
             stop
            end if
            data_dp_2d = merge(data_dp_2d, nodata_dp, mask_2d )
            call append(tmp_areaS,pack(data_dp_2d, mask_2d))
            L1_fLAI(:,ii) =  tmp_areaS(:)      
            deallocate(data_dp_2d,mask_2d,tmp_areaS)  
          end do
        end if
    !Sealed storage and direct runoff

        allocate(L1_fSealed(L1_nCells))
        L1_fSealed = 0.0_dp
        if (Sealed_LanduseNrs >= 1_i4) then
            do jj = 1, Sealed_LanduseNrs
                sloop: do ii =1,nLanduseTotal
                    if (ii == Sealed_LanduseIDs(jj)) then
                        L1_fSealed(:) = L1_fSealed(:) + L1_fLAI(:,ii)
                        exit sloop
                    end if
                end do sloop
            end do
        end if
    !********************************************************************************************** 
    !SOIL TYPE map and areal share of each grid (if nSoilperGrid > 1)
        !allocate
        allocate(L1_fsoil(L1_nCells,nSoilTotal))
        L1_fsoil = 0.0_dp

        if (nSoilperGrid == 1_i4) then
          fName = trim(adjustl(dir_catchinfo)) // trim(adjustl(fn_SoilTypes))
          funit = 117
          INQUIRE(file = fName, exist = file_exist)
          if (.not. file_exist) then
            call message()
            call message('***ERROR: since nSoilperGrid =1, soil map (fn_SoilTypes) should be provided!')
            stop
          end if
          call read_spatial_data_ascii(trim(fName), funit, &
            m_Level%ncols,     m_Level%nrows, m_Level%xllcorner,&
            m_Level%yllcorner, m_Level%cellsize, data_i4_2d, mask_2d)

          call message('    Reading soil information(fn_SoilTypes)...')       
           !check if mask matches
          if (.not. all(mask_2d .eqv. mask_global)) then
             call message()
             call message('***ERROR: soil map dose not match catchment domain')
          stop
          end if

          data_i4_2d = merge(data_i4_2d, nodata_i4, mask_2d )
          allocate(mask_tmp(L1_nCells))
          call append(tmp_ids, pack(data_i4_2d, mask_2d))
          do ii = 1, nSoilTotal
            mask_tmp(:) = .FALSE.
            where(tmp_ids == ii) mask_tmp = .TRUE.
            L1_fsoil(:,ii) = merge(1.0_dp,0.0_dp,mask_tmp)  
          end do
          deallocate(data_i4_2d,mask_2d,tmp_ids,mask_tmp)
        else
		  !read maps of areal share for each grid
          if (size(fn_fsoil) .ne. nSoilperGrid) then
            call message()
            call message('***ERROR: since nSoilperGrid > 1, each areal share maps should be provided!')
            stop
          end if
          call message('    Reading soil information(fn_fsoil) ...')  
          do ii = 1, nSoilperGrid
            fName = trim(adjustl(dir_catchinfo)) // trim(adjustl(fn_fsoil(ii)))
            funit = 127
            call read_spatial_data_ascii(trim(fName), funit, &
              m_Level%ncols,     m_Level%nrows, m_Level%xllcorner,&
              m_Level%yllcorner, m_Level%cellsize, data_dp_2d, mask_2d)

     
           !check if mask matches
            if (.not. all(mask_2d .eqv. mask_global)) then
             call message()
             call message('***ERROR: soil map dose not match catchment domain')
             stop
            end if
            data_dp_2d = merge(data_dp_2d, nodata_dp, mask_2d )
            call append(tmp_areaS,pack(data_dp_2d, mask_2d))
            L1_fsoil(:,ii) = tmp_areaS  
            deallocate(data_dp_2d,mask_2d, tmp_areaS)  
          end do
        end if

    !********************************************************************************************** 
    !SOIL DEPTH and Parameters FOR EACH LAYER
    ! the soil depth here is the lowerboundary depth
        !allocate
        allocate(L1_soildepth(L1_nCells,nSoilLayers))
        allocate(L1_wiltingpoint(L1_nCells,nSoilLayers))
        allocate(L1_saturatedSM(L1_nCells,nSoilLayers))
        L1_soildepth = 0.0_dp
        L1_wiltingpoint = 0.0_dp
        L1_saturatedSM = 0.0_dp
        if (size(fn_soilDepth) .ne. nSoilLayers) then
          call message()
          call message('***ERROR: soil depth map of each layer should be provided!')
          stop
        end if
        do ii = 1, nSoilLayers
          ! soil depth
            fName = trim(adjustl(dir_catchinfo)) // trim(adjustl(fn_soilDepth(ii)))
            funit = 119
            call read_spatial_data_ascii(trim(fName), funit, &
              m_Level%ncols,     m_Level%nrows, m_Level%xllcorner,&
              m_Level%yllcorner, m_Level%cellsize, data_dp_2d, mask_2d)

            !call message('    Reading soil depth information(fn_soilDepth) ...')       
           !check if mask matches
            if (.not. all(mask_2d .eqv. mask_global)) then
             call message()
             call message('***ERROR: soil depth  map dose not match catchment domain')
             stop
            end if
            data_dp_2d = merge(data_dp_2d, nodata_dp, mask_2d )
            call append(tmp_areaS,pack(data_dp_2d, mask_2d))
            L1_soildepth(:,ii) = tmp_areaS    
            deallocate(data_dp_2d,mask_2d,tmp_areaS)
          !soil wilting point
            fName = trim(adjustl(dir_catchinfo)) // trim(adjustl(fn_wiltpoint(ii)))
            funit = 129
            call read_spatial_data_ascii(trim(fName), funit, &
              m_Level%ncols,     m_Level%nrows, m_Level%xllcorner,&
              m_Level%yllcorner, m_Level%cellsize, data_dp_2d, mask_2d)

            !call message('    Reading soil information(fn_wiltpoint) ...')       
           !check if mask matches
            if (.not. all(mask_2d .eqv. mask_global)) then
             call message()
             call message('***ERROR: map dose not match catchment domain')
             stop
            end if
            data_dp_2d = merge(data_dp_2d, nodata_dp, mask_2d )
            call append(tmp_areaS,pack(data_dp_2d, mask_2d))
            L1_wiltingpoint(:,ii) = tmp_areaS * L1_soildepth(:,ii)  !from volumetric content[-] to [m]
            deallocate(data_dp_2d,mask_2d,tmp_areaS)
          !soil saturated soil moisture
            fName = trim(adjustl(dir_catchinfo)) // trim(adjustl(fn_satMoist(ii)))
            funit = 139
            call read_spatial_data_ascii(trim(fName), funit, &
              m_Level%ncols,     m_Level%nrows, m_Level%xllcorner,&
              m_Level%yllcorner, m_Level%cellsize, data_dp_2d, mask_2d)

            !call message('    Reading soil information(fn_satMoist) ...')       
           !check if mask matches
            if (.not. all(mask_2d .eqv. mask_global)) then
             call message()
             call message('***ERROR: map dose not match catchment domain')
             stop
            end if
            data_dp_2d = merge(data_dp_2d, nodata_dp, mask_2d )
            call append(tmp_areaS,pack(data_dp_2d, mask_2d))
            L1_saturatedSM(:,ii) = tmp_areaS * L1_soildepth(:,ii)  !from volumetric content[-] to [m]
            deallocate(data_dp_2d,mask_2d,tmp_areaS)
        end do

    !********************************************************************************************** 
    !CROP ROTATION map and areal share of each grid (if nSoilperGrid > 1)
        !allocate
        allocate(L1_frotation(L1_nCells,nRotationTotal))
        L1_frotation = 0.0_dp

        if (nRotationperGrid == 1_i4) then
          fName = trim(adjustl(dir_catchinfo)) // trim(adjustl(fn_RotationTypes))
          funit = 118
          INQUIRE(file = fName, exist = file_exist)
          if (.not. file_exist) then
            call message()
            call message('***ERROR: since nRotationperGrid =1, crop rotation map(fn_RotationTypes)should be provided!')
            stop
          end if
          call read_spatial_data_ascii(trim(fName), funit, &
            m_Level%ncols,     m_Level%nrows, m_Level%xllcorner,&
            m_Level%yllcorner, m_Level%cellsize, data_i4_2d, mask_2d)

          call message('    Reading crop rotation information(fn_RotationTypes)...')       
           !check if mask matches
          if (.not. all(mask_2d .eqv. mask_global)) then
             call message()
             call message('***ERROR: crop rotation map dose not match catchment domain')
          stop
          end if

          data_i4_2d = merge(data_i4_2d, nodata_i4, mask_2d )
          allocate(mask_tmp(L1_nCells))
          call append(tmp_ids, pack(data_i4_2d, mask_2d))
          do ii = 1, nRotationTotal
            mask_tmp(:) = .FALSE.
            where(tmp_ids == ii) mask_tmp = .TRUE.
            L1_frotation(:,ii) = merge(1.0_dp,0.0_dp,mask_tmp)
          end do
          deallocate(data_i4_2d,mask_2d,tmp_ids,mask_tmp)
        else
		  !read maps of areal share for each grid
          if (size(fn_frotation) .ne. nRotationperGrid) then
            call message()
            call message('***ERROR: since nRotationperGrid > 1, each areal share maps should be provided!')
            stop
          end if
          call message('    Reading crop rotation information(fn_frotation) ...')
          do ii = 1, nRotationperGrid
            fName = trim(adjustl(dir_catchinfo)) // trim(adjustl(fn_frotation(ii)))
            funit = 128
            call read_spatial_data_ascii(trim(fName), funit, &
              m_Level%ncols,     m_Level%nrows, m_Level%xllcorner,&
              m_Level%yllcorner, m_Level%cellsize, data_dp_2d, mask_2d)

           !check if mask matches
            if (.not. all(mask_2d .eqv. mask_global)) then
             call message()
             call message('***ERROR: crop rotation map dose not match catchment domain')
             stop
            end if
            data_dp_2d = merge(data_dp_2d, nodata_dp, mask_2d )
            call append(tmp_areaS,pack(data_dp_2d, mask_2d))
            L1_frotation(:,ii) = tmp_areaS     
            deallocate(data_dp_2d,mask_2d, tmp_areaS)  
          end do
        end if



    !********************************************************************************************** 
    !gauging network information -->routing resolution
        fName = trim(adjustl(dir_gauges)) // trim(adjustl(fn_gauge_loc))
        funit = 311
        call read_spatial_data_ascii(trim(fName), funit, &
            r_Level%ncols,     r_Level%nrows, r_Level%xllcorner,&
            r_Level%yllcorner, r_Level%cellsize, data_i4_2d, mask_2d)

        call message('    Reading gauging network information ...')
        data_i4_2d = merge(data_i4_2d,  nodata_i4, mask_2d)
        call append(L1_gaugeLoc, pack(data_i4_2d, mask_global_r))
        !check input gauge id locations and if id matches the provided data
        do idgauges = 1, nHydroGauges
            if (.not. any(data_i4_2d .EQ. basin_gauge%gaugeid(idgauges))) then
                call message()
                call message('***ERROR: gauge_id "', trim(adjustl(num2str(basin_gauge%gaugeid(idgauges)))), &
                     '" not found in ' )
                call message('Gauge location input file: ', &
                     trim(adjustl(dir_gauges))//trim(adjustl(fn_gauge_loc)))
                stop
            end if
            !get the right index of gauging ids
            do i = 1, size(L1_gaugeLoc)
                if (L1_gaugeLoc(i) .eq. basin_gauge%gaugeid(idgauges)) then
                  basin_gauge%gaugeid_loc(idgauges) = i
                end if
            end do
        end do
        do idgauges = 1, nWQGauges
            if (.not. any(data_i4_2d .EQ. basin_gauge%wqm_gaugeid(idgauges))) then
                call message()
                call message('***ERROR: wqm_gauge_id "', trim(adjustl(num2str(basin_gauge%wqm_gaugeid(idgauges)))), &
                     '" not found in ' )
                call message('Gauge location input file: ', &
                     trim(adjustl(dir_gauges))//trim(adjustl(fn_gauge_loc)))
                stop
            end if
            !get the right index of gauging ids
            do i = 1, size(L1_gaugeLoc)
                if (L1_gaugeLoc(i) .eq. basin_gauge%wqm_gaugeid(idgauges)) then
                  basin_gauge%wqm_gaugeid_loc(idgauges) = i
                end if
            end do
        end do
        do idgauges = 1, nInflowGauges
            if (.not. any(data_i4_2d .EQ. basin_gauge%inflowgaugeid(idgauges))) then
                call message()
                call message('***ERROR: inflowgaugeid "', trim(adjustl(num2str(basin_gauge%inflowgaugeid(idgauges)))), &
                     '" not found in ' )
                call message('Gauge location input file: ', &
                     trim(adjustl(dir_gauges))//trim(adjustl(fn_gauge_loc)))
                stop
            end if
            !get the right index of gauging ids
            do i = 1, size(L1_gaugeLoc)
                if (L1_gaugeLoc(i) .eq. basin_gauge%inflowgaugeid(idgauges)) then
                  basin_gauge%inflowgaugeid_loc(idgauges) = i
                end if
            end do
        end do

       !deallocate arrays
        deallocate(data_i4_2d, mask_2d)

    !-----------------------------------------------------------------------------------
    !***********************************************************************************
    !Semi-distributed model structure
    !***********************************************************************************
    else if (model_type == 2) then
    !modeling domain
        fName = trim(adjustl(dir_catchinfo)) // trim(adjustl(fn_catchDomain) )
        funit =111
        open(funit, file= fName, action='read')
        read(funit, *) dummy, L1_nSubs
        read(funit, *) dummy, L1_nHRUs   !used for the total HRU numbers defined
        read(funit, *) dummy              !title line 
        L1_nCells =  L1_nHRUs * L1_nSubs !total number of terrestrial interations
		!allocate
        if (.not. allocated(Model_mask) ) allocate(Model_mask(L1_nHRUs, L1_nSubs) )
        if (.not. allocated(SubID) )    allocate(SubID(L1_nSubs) )
        if (.not. allocated(SubArea) )  allocate(SubArea(L1_nSubs) )
        if (.not. allocated(SubHRUs) )  allocate(SubHRUs(L1_nHRUs,L1_nSubs) )
        if (.not. allocated(file_initINconc))   allocate(file_initINconc(L1_nSubs))
        if (.not. allocated(file_initINconc))   allocate(file_initONconc(L1_nSubs))
        if (.not. allocated(CropRotationID) )   allocate(CropRotationID(L1_nSubs))
        !assign global variables L1_**
        if (.not. allocated(init_concIN) )  allocate(init_concIN(L1_nCells) )
        if (.not. allocated(init_concON) )  allocate(init_concON(L1_nCells) )
        if (.not. allocated(L1_area))       allocate(L1_area(L1_nCells))
        L1_area = 0.0_dp
        init_concIN = 0.0_dp
        init_concON = 0.0_dp
        if (.not. allocated(L1_frotation))  allocate(L1_frotation(L1_nCells,nRotationTotal))
        L1_frotation = 0.0_dp
        do i =1, L1_nSubs
            read(funit, *) SubID(i), CropRotationID(i), SubArea(i), (SubHRUs(j,i),j=1, L1_nHRUs),&
                           file_initINconc(i),file_initONconc(i)
          do j=1, L1_nHRUs
            L1_area(L1_nHRUs*(i-1)+j) = SubArea(i) * (SubHRUs(j,i))  ![m2]
            init_concIN(L1_nHRUs*(i-1)+j) = file_initINconc(i)
            init_concON(L1_nHRUs*(i-1)+j) = file_initONconc(i)
            inloop: do ii = 1, nRotationTotal
                if (ii == CropRotationID(i)) then
                  L1_frotation(L1_nHRUs*(i-1)+j,ii) = 1.0_dp
                  exit inloop
                end if
            end do inloop
          end do
        end do
        close(funit)
    !river network routing order
          fName = trim(adjustl(dir_catchinfo)) // trim(adjustl(fn_spatialorganization) )
          funit =112
          call message(' Reading spatial organization file...')
          open(funit, file= fName, action='read')
          read(funit, *) dummy, L1_nLinks
          !check nLinks should be L1_nSubs
          !note the outlet should be the case of fromNode == toNode
          if (L1_nLinks .ne. L1_nSubs) then
           call message()
           call message('***ERROR: number of connection-links dose not match catchment domain')
           stop            
          end if
          allocate(L1_chanmask(L1_nLinks))
          L1_chanmask(:) = .TRUE.
          allocate(R_areaCell(L1_nLinks))
          R_areaCell(:) = SubArea(:) !for each subcatchment
          allocate(fromNode(L1_nLinks))
          allocate(toNode(L1_nLinks))
          allocate(Rsequence(L1_nLinks))
          allocate(netPerm(L1_nLinks))
          allocate(slength(L1_nLinks))
          allocate(swidth(L1_nLinks))
          read(funit, *) dummy
          do i =1, L1_nLinks
            read(funit, *) dummy, fromNode(i),toNode(i),Rsequence(i),slength(i), swidth(i)
            !link id, fromcell id, tocell id, sequence id
            !the outlet node ids
            if (fromNode(i) == toNode(i)) then
               OutNode = fromNode(i)
            end if
          end do
          close(funit)
          !assign information to global variables
          call append(L1_fromN,fromNode(:))       
          call append(L1_toN,toNode(:))
          call append(L1_rOrder,Rsequence(:))
          do i=1,L1_nLinks
            netPerm(Rsequence(i)) = i
          end do
          call append(L1_netPerm,netPerm(:))
         !channel morphology info
          call append(L1_length,slength(:) )
          call append(L1_width,swidth(:) )  
    !read in geo properties of each HRU (land use type, soil type and soil properties)
        fName = trim(adjustl(dir_catchinfo)) // trim(adjustl(fn_classProperties) )
        funit =113
        allocate(geoProperty%hruids(L1_nHRUs) )
        allocate(geoProperty%landuseids(L1_nHRUs) )
        allocate(geoProperty%soilids(L1_nHRUs) )
        allocate(geoProperty%soildepth(L1_nHRUs, nSoilLayers))
        allocate(geoProperty%wiltingpoint(L1_nHRUs, nSoilLayers))
        allocate(geoProperty%saturatedSM(L1_nHRUs, nSoilLayers))

        open(funit, file= fName, action='read')
        read(funit, *) dummy
        read(funit, *) dummy
        read(funit, *) dummy
        do i =1, L1_nHRUs
            read(funit, *) geoProperty%hruids(i), geoProperty%landuseids(i), geoProperty%soilids(i), &
            (geoProperty%soildepth(i, j),j=1, nSoilLayers), (geoProperty%wiltingpoint(i, j),j=1, nSoilLayers), &
            (geoProperty%saturatedSM(i, j),j=1, nSoilLayers)
        end do
        close(funit)
        !assign global variable 
        allocate(L1_soildepth(L1_nCells,nSoilLayers))
        allocate(L1_wiltingpoint(L1_nCells,nSoilLayers))
        allocate(L1_saturatedSM(L1_nCells,nSoilLayers))
        L1_soildepth = 0.0_dp
        L1_wiltingpoint = 0.0_dp
        L1_saturatedSM = 0.0_dp 
        !assign to global variables as the same as the grid structure
        allocate(L1_fLAI(L1_nCells,nLanduseTotal))
        L1_fLAI = 0.0_dp
        allocate(L1_fsoil(L1_nCells,nSoilTotal))
        L1_fsoil = 0.0_dp

        do i = 1,L1_nCells, L1_nHRUs
          do j = 1, L1_nHRUs
            inloop2: do ii = 1, nLanduseTotal
                if (ii == geoProperty%landuseids(j)) then
                  L1_fLAI(i+j-1,ii) = 1.0_dp
                  exit inloop2
                end if
            end do inloop2
            inloop3: do ii = 1, nSoilTotal
                if (ii == geoProperty%soilids(j)) then
                  L1_fsoil(i+j-1,ii) = 1.0_dp
                  exit inloop3
                end if
            end do inloop3
            !
            L1_soildepth(i+j-1,:) = geoProperty%soildepth(j,:)
            L1_wiltingpoint(i+j-1,:) = geoProperty%wiltingpoint(j,:) * geoProperty%soildepth(j,:) !content[-] to [m]
            L1_saturatedSM(i+j-1,:) = geoProperty%saturatedSM(j,:) * geoProperty%soildepth(j,:)   !content[-] to [m]
          end do
        end do

    !Sealed storage and direct runoff
        allocate(L1_fSealed(L1_nCells))
        L1_fSealed = 0.0_dp
        if (Sealed_LanduseNrs >= 1_i4) then
            do jj = 1, Sealed_LanduseNrs
                sloop1: do ii =1,nLanduseTotal
                    if (ii == Sealed_LanduseIDs(jj)) then
                        L1_fSealed(:) = L1_fSealed(:) + L1_fLAI(:,ii)
                        exit sloop1
                    end if
                end do sloop1
            end do
        end if

    !gauge locations : the ids in config.nml is the subcatchmend ID
    !so gaugeid_loc = gauge_id
        basin_gauge%gaugeid_loc(:) = basin_gauge%gaugeid(:)
        basin_gauge%wqm_gaugeid_loc(:) = basin_gauge%wqm_gaugeid(:)
        basin_gauge%inflowgaugeid_loc(:) = basin_gauge%inflowgaugeid(:)

    else
        call message('***ERROR: please specify the "model_type" as either 1 or 2! ')
        stop
    end if  ! model_type 

    !*****************
	!READ INITIA VALUES OF LAND-USE DENPENDENT VARIABLES
	!*****************
    fName = trim(adjustl(dir_catchinfo)) // trim(adjustl(fn_initial_condition) )
    funit =211
    open(funit, file= fName, action='read')
    read(funit, *) dummy, init_type
    read(funit, *) dummy, nClasses
    read(funit, *) dummy
    !check category dependency
    select case (init_type)
    case(1)  ! landuse dependent
        if (nClasses .ne. nLanduseTotal) then
            call message()
            call message('***ERROR: initial values are assigned as landuse dependent, but do not match land use information')
            stop
        end if
    case(2)   !soil-type dependent
        if (nClasses .ne. nSoilTotal) then
            call message()
            call message('***ERROR: initial values are assigned as soil type dependent, but do not match soil type information')
            stop
        end if
    end select
    !allocate
    if (.not. allocated(init_humusN) ) allocate(init_humusN(nClasses) )
    if (.not. allocated(init_fastN) )  allocate(init_fastN(nClasses) )
    if (.not. allocated(hnhalf) )      allocate(hnhalf(nClasses) )
    do i =1, nClasses
       read(funit, *) dummy, init_fastN(i), init_humusN(i), hnhalf(i)
    end do
    close(funit)


    !*****************
    !READ CROPDATA
    !*****************
    fName = trim(adjustl(dir_catchinfo)) // trim(adjustl(fn_cropinfo) )
    funit = 212
    open(funit, file = fName, action='read' )

    !read header
    read(funit, *) dummy, num_crops
    read(funit, *) dummy
    allocate(cropdata(num_crops) )
	
    do i =1, num_crops
       read(funit, *) cropdata(i)%cropname, cropdata(i)%cropid, cropdata(i)%frtn1, cropdata(i)%frtday1, cropdata(i)%frtdown1,&
          cropdata(i)%frtn2, cropdata(i)%frtday2, cropdata(i)%frtdown2, cropdata(i)%mann1, cropdata(i)%manday1, &
          cropdata(i)%mandown1, cropdata(i)%mann2, cropdata(i)%manday2, cropdata(i)%mandown2, cropdata(i)%manfIN,&			   
          cropdata(i)%frtperiod, cropdata(i)%resn, cropdata(i)%resday, cropdata(i)%resdown, cropdata(i)%resfast,   & 
          cropdata(i)%resperiod,cropdata(i)%up1, cropdata(i)%up2, cropdata(i)%up3, cropdata(i)%uppsoil, cropdata(i)%deepsoil,&  
          cropdata(i)%plantd,cropdata(i)%emergd, cropdata(i)%harvestd, cropdata(i)%ccrop, &
          cropdata(i)%ccplantd, cropdata(i)%cchavestd

    !applied amount convert units: from kg/ha to g/m2
    cropdata(i)%frtn1 = cropdata(i)%frtn1 * 0.1_dp 
    cropdata(i)%frtn2 = cropdata(i)%frtn2 * 0.1_dp
    cropdata(i)%mann1 = cropdata(i)%mann1 * 0.1_dp 
    cropdata(i)%mann2 = cropdata(i)%mann2 * 0.1_dp
    cropdata(i)%resn  = cropdata(i)%resn * 0.1_dp
  
    end do  
    close(funit)    

    !*****************
    !READ MONTHLY LAI INFORMATION 
    !*****************
    fName = trim(adjustl(dir_catchinfo)) // trim(adjustl(fn_lut_monLAI) )
    funit = 213
    open(funit, file = fName, action='read' )

    !read header
    read(funit, *) dummy, nLAI
    read(funit, *) dummy
    allocate(LAIUnitList(nLAI) )
    allocate(LAILUT(nLAI, int(12,i4)))
 
    do i =1, nLAI
       read(funit, *) LAIUnitList(i), dummy, (LAILUT(i,j),j=1,int(12,i4))
    end do
    close(funit)
	
    !*****************
    !READ ROTATION INFORMATION 
    !*****************
    fName = trim(adjustl(dir_catchinfo)) // trim(adjustl(fn_lut_rotation) )
    funit = 214
    open(funit, file = fName, action='read' )

    !read header
    read(funit, *) dummy, nCroptation
    read(funit, *) dummy
    if (nCroptation .ne. nRotationTotal) then
        call message()
        call message('***ERROR: nRotationTotal does not match the provided information')
        stop
    end if
    allocate(rotation%id(nCroptation) )
    allocate(rotation%ncrops(nCroptation) )
    allocate(rotation%crop(nCroptation, maxnCropsInRotation) ) !maximum 10 crop in a rotation sequence
    do i =1, nCroptation
       read(funit, *) rotation%id(i),dummy, rotation%ncrops(i), &
          (rotation%crop(i, j),j=1, rotation%ncrops(i)) 
    end do
    close(funit)
	

	!##################################################################
	!## BASIN AVERAGE GLOBAL RADIATION (for in-stream shading effect)##
	!##################################################################
    fName = trim(adjustl(dir_catchinfo )) // trim(adjustl(fn_globalRadiation) )
	!to check if the "actual_cropNuptk.txt" exists or not
    INQUIRE(file = fName, exist = file_exist)
    GR_file_exist = file_exist  !update the variable value
    if (GR_file_exist) then
      fName = trim(adjustl(dir_catchinfo )) // trim(adjustl(fn_globalRadiation) )
      funit = 215

      start_period = (/simPer%yStart, simPer%mStart, simPer%dStart/)
      end_period   = (/simPer%yEnd,   simPer%mEnd,   simPer%dEnd  /) 
      open(funit, file = fName, action = 'read')
        ! read header
        read(funit,'(a256)') dummy
        read(funit,*)        dummy, nodata_file
        read(funit,*)        dummy, (start_period_file(i), i = 1, 3)
        read(funit,*)        dummy, (end_period_file(i),   i = 1, 3)
        read(funit,'(a256)') line     
        dummy = dummy//''   ! only to avoid warning 
	    !
        start_jul = julday(start_period(3),start_period(2),start_period(1))
        end_jul   = julday(end_period(3),end_period(2),end_period(1)) 
        start_jul_file = julday(start_period_file(3),start_period_file(2),start_period_file(1))
        end_jul_file   = julday(end_period_file(3),end_period_file(2),end_period_file(1)) 
	  !
        if ((start_jul < start_jul_file) .or. (end_jul > end_jul_file)) then
           print*, 'simulation period is not fully covered by global radiation data,' 
           print*, 'average value is given by default!!'
        end if
	    !allocate period 
        allocate(gr_file(end_jul_file - start_jul_file +1_i4))
        allocate(gr_mask(end_jul_file - start_jul_file +1_i4))
        allocate(gr_data(end_jul - start_jul + 1_i4))

	    !default value
        gr_file = nodata_file
        gr_mask = .TRUE.

	  
        do i = 1, (end_jul_file - start_jul_file + 1_i4)
           read(funit, *) (file_time(j), j = 1,3), gr_file(i)
        end do 
        where (abs(gr_file - nodata_file) .lt. tiny(1.0_dp)) 
           gr_mask = .FALSE.
        end where
        avg_gr = sum(gr_file, gr_mask)/ max(1,count(gr_mask))
        !set average global radiation value as default
        gr_data(:) = avg_gr
	  
        length_file   = (end_jul_file   - start_jul_file   + 1 )
        length_period = (end_jul - start_jul + 1 )
        !Merge file data to global variable
        !       |---------------------------------|    FILE
        !                  |--------------|            PERIOD
        if (( start_jul .ge. start_jul_file ) .and. ( end_jul .le. end_jul_file )) then
           idx_st_period = 1
           idx_en_period = length_period 
           idx_st_file   = start_jul - start_jul_file + 1  
           idx_en_file   = idx_st_file + length_period     - 1 
        end if

        !                  |--------------|            FILE
        !       |---------------------------------|    PERIOD
        if (( start_jul .lt. start_jul_file ) .and. ( end_jul .gt. end_jul_file )) then
           idx_st_period = start_jul_file - start_jul + 1  
           idx_en_period = idx_st_period + length_file     - 1 
           idx_st_file   = 1
           idx_en_file   = length_file      
        end if

        !  |--------------|                            FILE
        !       |---------------------------------|    PERIOD
        if (( start_jul .ge. start_jul_file ) .and. ( end_jul .gt. end_jul_file )) then
           idx_st_period = 1
           idx_en_period =  end_jul_file     - start_jul + 1  
           idx_st_file   =  start_jul - start_jul_file   + 1  
           idx_en_file   = length_file   
        end if

        !                          |--------------|    FILE
        !  |---------------------------------|         PERIOD
        if (( start_jul.lt. start_jul_file ) .and. ( end_jul .le. end_jul_file )) then
           idx_st_period =  start_jul_file - start_jul + 1  
           idx_en_period = length_period
           idx_st_file   = 1
           idx_en_file   =  end_jul - start_jul_file   + 1 
        end if
	    !update global variable when measured data are available
        gr_data(idx_st_period:idx_en_period) = gr_file(idx_st_file: idx_en_file)
        call append(global_radiation, gr_data)
	  
        deallocate(gr_file)
        deallocate(gr_mask)
        deallocate(gr_data)
      end if
    
    !############################################
    !##Riparian zone shading effect coefficient##
	!############################################
	!to get the coefficient of riparian zone shading effect: rz_coeff
    if (GR_file_exist) then
      call read_shading_data(global_radiation, LAILUT, norlai_daily, nor_globalradi)
    end if
			  
  end subroutine wqm_readcatchdata

END MODULE mo_wqm_readcatchdata  
