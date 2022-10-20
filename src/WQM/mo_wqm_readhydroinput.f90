!> \file mo_wqm_readhydroinput.f90

!> \brief Reads hyrologic information as the input of WQM.

!> \details This module is to read

!> \modified from mHM inplementations of reading netcdf and binary climate forcing

!> \authors Xiaoqiang Yang
!> \date Sep 2021


module mo_wqm_readhydroinput
  implicit none
  public :: wqm_readhydroinput
  !
contains


  ! ------------------------------------------------------------------

  !     NAME
  !         wqm_readhydroinput

  !     PURPOSE
  !>        \brief Reads forcing input in NetCDF file format.

  !>        \details Reads netCDF forcing files.  \n
  !>        First, the dimensions given are cross-checked with header.txt information. Second, the data of the
  !>        specified period are read from the specified directory.
  !>        If the optional lower and/or upper bound for the data values is given, the read data are checked for validity.
  !>        The program is stopped if any value lies out of range.\n
  !>        If the optinal argument nocheck is true, the data are not checked for coverage with the input mask.
  !>        Additionally in this case an mask of vild data points can be received from the routine in maskout.

 
  !     RESTRICTIONS
  !>        \note Files have to be called like defined in mo_files. Furthermore the variable names have to be called
  !>              like they are defined in the declaration of this subroutine. The NetCDF file has to have 3 dimensions:
  !>              1. x, 2. y, 3. t. It is expected that the variables (especially)within the NetCDF files contain an
  !>              unit attribute. The timestep has to be equidistant


  subroutine wqm_readhydroinput()

    use mo_kind,             only: i4, dp
    use mo_julian,           only: caldat, julday
    use mo_message,          only: message

    use mo_string_utils,      only: num2str
    use mo_utils,             only: eq, ne

    use mo_wqm_global_variables,    only: &
         simPer, nTstepDay, &
         dir_hydroInput, inputFormat,inputhydro_filename,  &
         inputFormat_rout, inputhydro_rout_filename, inputfluxUT,input_rout_fluxUT,ETfluxUT, &
	!input variable names
         nSoilLayers, nEvapLayer, nTranspLayer, nSubsurfaceflow, nSubsrfToChannel, &
         L1_saturatedSM,L1_wiltingpoint,L1_soildepth, &
         vn_precipitation,vn_airTemperature, vn_WS_interception,vn_evapCanopy, vn_throughfall, &
         vn_infiltration, vn_WS_snowpack, vn_snowmelt, vn_WS_ponding,&
         vn_WS_ponding_ver,    &
         vn_ETevap, vn_ETtransp, &
         vn_surfaceflow,vn_WS_surfaceflow, SrfOVFmixratio,SrfOVFmixratio_vn,L1_SrfmixRatio,  &
         vn_subsurfaceflow,vn_WS_subsurfaceflow, &
         vn_baseflow, vn_WS_baseflow, &
         vn_directflow,vn_Esealed,vn_WS_directflow, &
       !hydro structure
         vn_WS_interception_ver, vn_WS_snowpack_ver, vn_WS_surfaceflow_ver, &
         vn_WS_subsurfaceflow_ver, vn_WS_baseflow_ver, &
         vn_WS_surfaceflow_lat,  &
         vn_WS_subsurfaceflow_lat,vn_WS_baseflow_lat,   &
         vn_WS_soil,vn_SatSM, vn_WP,vn_WS_soil_ver, &
         vn_reinfiltration,vn_reinfil_soil, vn_soilTemperature,&
       !channel
         vn_chan_store, vn_chan_in, vn_chan_out, &
         vn_srftochan,vn_baseftochan,  &
         vn_subsrftochan,  &
    !variable name for computing		 
        !process related, dim1= ncells, dim2 = ntimestep 
        processCase,L1_rainfall_direct, L1_temp, L1_throughfall, L1_snowmelt, L1_evapCanopy, &  
        L1_infiltrate, L1_surfaceRunoff, L1_baseflow, &
        L1_reinfiltSrf, L1_reinfiltSoil,&
        L1_WSsealed, L1_Esealed, L1_directRunoff, &
        !dim1= ncells, dim2 = ntimestep, dim3 = number of subsurface flow components
        L1_subsurfaceRunoff,  &
        !ReservoirStore type
        L1_WSintercept, L1_WSsnowpack,L1_WSponding,L1_WSsurface, L1_WSbase,  L1_WSsubsurface, L1_WSsoil, &
        !dim1= ncells, dim2 = ntimestep, dim3 = number of involved soil layers
        L1_evapSoil,L1_transp, L1_soiltemp,  &
        !channel, dim1 = ncells, dim2 = ntimestep
        L1_WSchannel, L1_WSchannel_upin, L1_WSchannel_downout, &
        !dim1 = dim1 = ncells, dim2 = ntimestep
        L1_srfTOchan, L1_basefTOchan,   &
        !dim1 = dim1 = ncells, dim2 = ntimestep, dim3 = number OF SUBSURFACE runoffs
        L1_subsrfTOchan


    implicit none

    ! local variables
    !real(dp), dimension(:,:), allocatable  :: hydrodata !directly from input
    real(dp), dimension(:,:), allocatable    :: hydrodata_packed !packed to 2-D
    character(256)                         :: varName   ! name of variable
    character(256)            :: fdir, fName        ! name of dir/file
    integer(i4)               :: iver
    logical                   :: isnc
    integer(i4)               :: simPerL
    !integer(i4)               :: datatype

    ! simulation period length
    simPerL = (simPer%julEnd - simPer%julStart + 1 ) * nTstepDay


    if (inputFormat .eq. "nc") then
        fName = trim(dir_hydroInput) // trim(inputhydro_filename) // '.nc'
        isnc = .TRUE.
    else if (inputFormat .eq. "bin") then
        fdir = dir_hydroInput
        isnc = .FALSE.
    else
        call message()
        call message('***ERROR: the hydrological input files should be in either netCDF or binary format!')
        stop
    end if

      !if soil moisture constants provided by external hydro model vn_SatSM, vn_WP
      !this differs from below flux inputs, since they won't change with time steps
        if ((nSoilLayers >= 1) .and. (any(vn_SatSM .ne. "none"))) then
            do iver =1, nSoilLayers
                varName = vn_SatSM(iver)
                if (isnc) then
                    call read_hydro_nc(fName, varName, 1, hydrodata_packed) !take the first-step value is enough
                else
                    call read_hydro_bin(fdir, varName, 1, hydrodata_packed)
                end if
                L1_saturatedSM(:,iver) = hydrodata_packed(:,1) * L1_soildepth(:,iver) 
                ! free memory immediately
                deallocate(hydrodata_packed)
            end do
        end if
        !the externally supplied wilting point info
        if ((nSoilLayers >= 1) .and. (any(vn_WP .ne. "none"))) then
            do iver =1, nSoilLayers
                varName = vn_WP(iver)
                if (isnc) then
                    call read_hydro_nc(fName, varName, 1, hydrodata_packed) !take the first-step value is enough
                else
                    call read_hydro_bin(fdir, varName, 1, hydrodata_packed)
                end if
                L1_wiltingpoint(:,iver) = hydrodata_packed(:,1) * L1_soildepth(:,iver) 
                ! free memory immediately
                deallocate(hydrodata_packed)
            end do
        end if

        !the externally supplied surface flow conc mixing ratio info
        if (processCase(7) == 1 ) then
          L1_SrfmixRatio(:) = SrfOVFmixratio ! the unique value specified from the configuration file
          !if it's externally supplied by the hydro model, then above value will be over-written
          if (SrfOVFmixratio_vn .ne. "none") then
            varName = SrfOVFmixratio_vn
            if (isnc) then
                call read_hydro_nc(fName, varName, 1, hydrodata_packed) !take the first-step value is enough
            else
                call read_hydro_bin(fdir, varName, 1, hydrodata_packed)
            end if
            L1_SrfmixRatio(:) = hydrodata_packed(:,1)
            ! free memory immediately
            deallocate(hydrodata_packed)
          end if
        end if
      !Start input flux time series
        !read data for each provided variable names
        if (vn_precipitation .ne. "none") then
            varName = vn_precipitation
            if (isnc) then
                call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
            else
                call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
            end if
            L1_rainfall_direct = hydrodata_packed * inputfluxUT
            ! free memory immediately
            deallocate(hydrodata_packed)
        end if

        if (vn_airTemperature .ne. "none") then
            varName = vn_airTemperature
            if (isnc) then
                call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
            else
                call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
            end if
            L1_temp = hydrodata_packed
            ! free memory immediately
            deallocate(hydrodata_packed)
        end if

        if (vn_evapCanopy .ne. "none") then
            varName = vn_evapCanopy
            if (isnc) then
                call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
            else
                call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
            end if 
            L1_evapCanopy = hydrodata_packed * ETfluxUT
            ! free memory immediately
            deallocate(hydrodata_packed)
        end if

        if (vn_throughfall .ne. "none") then
            varName = vn_throughfall
            if (isnc) then
                call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
            else
                call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
            end if
            L1_throughfall = hydrodata_packed * inputfluxUT
            ! free memory immediately
            deallocate(hydrodata_packed)
        end if

        if (vn_snowmelt .ne. "none") then
            varName = vn_snowmelt
            if (isnc) then
                call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
            else
                call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
            end if
            L1_snowmelt = hydrodata_packed * inputfluxUT
            ! free memory immediately
            deallocate(hydrodata_packed)
        end if

        if (vn_infiltration .ne. "none") then
            varName = vn_infiltration
            if (isnc) then
                call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
            else
                call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
            end if
            L1_infiltrate = hydrodata_packed * inputfluxUT
            ! free memory immediately
            deallocate(hydrodata_packed)
        end if

        if (vn_surfaceflow .ne. "none") then
            varName = vn_surfaceflow
            if (isnc) then
                call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
            else
                call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
            end if
            L1_surfaceRunoff = hydrodata_packed * inputfluxUT
            ! free memory immediately
            deallocate(hydrodata_packed)
        end if

        if (vn_baseflow .ne. "none") then
            varName = vn_baseflow
            if (isnc) then
                call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
            else
                call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
            end if
            L1_baseflow = hydrodata_packed * inputfluxUT
            ! free memory immediately
            deallocate(hydrodata_packed)
        end if
        !soil evap and transpiration are allocated for all soil layers, but only flux values of involved layers
        !are assigned according to the hydro input, the rest is 0.0_dp by default 
        if (nEvapLayer > 0)  then
            do iver =1, nEvapLayer
                varName = vn_ETevap(iver)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_evapSoil(:,:,iver) = hydrodata_packed * ETfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)
            end do
        end if
        if (nTranspLayer > 0) then
            do iver =1, nTranspLayer
                varName = vn_ETtransp(iver)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_transp(:,:,iver) = hydrodata_packed * ETfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)
            end do
        end if

      !if re-infiltration allowed
        if (vn_reinfiltration .ne. "none") then
            varName = vn_reinfiltration
            if (isnc) then
                call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
            else
                call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
            end if
            L1_reinfiltSrf = hydrodata_packed * inputfluxUT
            ! free memory immediately
            deallocate(hydrodata_packed)
        end if

      !if soil temperature is provided by hydro inputs
        if (vn_soilTemperature .ne. "none") then
            varName = vn_soilTemperature
            if (isnc) then
                call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
            else
                call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
            end if
            L1_soiltemp = hydrodata_packed
            ! free memory immediately
            deallocate(hydrodata_packed)
        end if

        if ((nSoilLayers >= 1) .and. (processCase(6) == 1)) then
            do iver =1, nSoilLayers
                if (vn_reinfil_soil(iver) .ne. "none") then !the last layer might not have reinfil flux
                  varName = vn_reinfil_soil(iver)
                  if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                  else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                  end if
                  L1_reinfiltSoil(:,:,iver) = hydrodata_packed * inputfluxUT
                  ! free memory immediately
                  deallocate(hydrodata_packed)
                end if
            end do
        end if

        if ((nSubsurfaceflow >= 1) .and. any(vn_subsurfaceflow .ne. "none")) then
            do iver =1, nSubsurfaceflow
                varName = vn_subsurfaceflow(iver)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_subsurfaceRunoff(:,:,iver) = hydrodata_packed * inputfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)
            end do
        end if
    !Direct runoff from sealed storage vn_WS_directflow
        if (vn_WS_directflow .ne. "none") then
            varName = vn_WS_directflow
            if (isnc) then
                call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
            else
                call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
            end if
            L1_WSsealed = hydrodata_packed * inputfluxUT
            ! free memory immediately
            deallocate(hydrodata_packed)            
        end if
        if (vn_Esealed .ne. "none") then
            varName = vn_Esealed
            if (isnc) then
                call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
            else
                call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
            end if
            L1_Esealed = hydrodata_packed * inputfluxUT
            ! free memory immediately
            deallocate(hydrodata_packed)
        end if
        if (vn_directflow .ne. "none") then
            varName = vn_directflow
            if (isnc) then
                call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
            else
                call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
            end if
            L1_directRunoff = hydrodata_packed * inputfluxUT
            ! free memory immediately
            deallocate(hydrodata_packed)
        end if

    !vn_WS_interception
        if (vn_WS_interception .ne. "none") then
            varName = vn_WS_interception
            if (isnc) then
                call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
            else
                call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
            end if
            L1_WSintercept%reservoir = hydrodata_packed * inputfluxUT
            ! free memory immediately
            deallocate(hydrodata_packed)
        end if
        if (any(vn_WS_interception_ver .ne. "none")) then
            !four vertical directions
            if (vn_WS_interception_ver(1) .ne. "none")  then
                varName = vn_WS_interception_ver(1)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_WSintercept%up_in = hydrodata_packed * inputfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)
            end if
            if (vn_WS_interception_ver(2) .ne. "none")  then
                varName = vn_WS_interception_ver(2)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_WSintercept%down_out = hydrodata_packed * inputfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)
            end if
            if (vn_WS_interception_ver(3) .ne. "none")  then
                varName = vn_WS_interception_ver(3)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_WSintercept%up_out = hydrodata_packed * inputfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)
            end if
            ! if (vn_WS_interception_ver(4) .ne. "none")  then
                ! varName = vn_WS_interception_ver(4)
                ! if (isnc) then
                    ! call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                ! else
                    ! call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                ! end if
                ! L1_WSintercept%down_in = hydrodata_packed
                ! ! free memory immediately
                ! deallocate(hydrodata_packed)
            ! end if
        end if

        ! if (any(vn_WS_interception_lat .ne. "none")) then
            ! !four vertical directions
            ! if (vn_WS_interception_lat(1) .ne. "none")  then
                ! varName = vn_WS_interception_lat(1)
                ! if (isnc) then
                    ! call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                ! else
                    ! call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                ! end if
                ! L1_WSintercept%lateral_in = hydrodata_packed
                ! ! free memory immediately
                ! deallocate(hydrodata_packed)
            ! end if
            ! if (vn_WS_interception_lat(2) .ne. "none")  then
                ! varName = vn_WS_interception_lat(2)
                ! if (isnc) then
                    ! call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                ! else
                    ! call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                ! end if
                ! L1_WSintercept%lateral_out = hydrodata_packed
                ! ! free memory immediately
                ! deallocate(hydrodata_packed)
            ! end if
        ! end if
    !vn_WS_snowpack

        if (vn_WS_snowpack .ne. "none") then
            varName = vn_WS_snowpack
            if (isnc) then
                call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
            else
                call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
            end if
            L1_WSsnowpack%reservoir = hydrodata_packed * inputfluxUT
            ! free memory immediately
            deallocate(hydrodata_packed)
        end if

        if (any(vn_WS_snowpack_ver .ne. "none")) then
            !four vertical directions
            if (vn_WS_snowpack_ver(1) .ne. "none")  then
                varName = vn_WS_snowpack_ver(1)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_WSsnowpack%up_in = hydrodata_packed * inputfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)
            end if
            if (vn_WS_snowpack_ver(2) .ne. "none")  then
                varName = vn_WS_snowpack_ver(2)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_WSsnowpack%down_out = hydrodata_packed * inputfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)
            end if
            if (vn_WS_snowpack_ver(3) .ne. "none")  then
                varName = vn_WS_snowpack_ver(3)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_WSsnowpack%up_out = hydrodata_packed * inputfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)
            end if
            ! if (vn_WS_snowpack_ver(4) .ne. "none")  then
                ! varName = vn_WS_snowpack_ver(4)
                ! if (isnc) then
                    ! call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                ! else
                    ! call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                ! end if
                ! L1_WSsnowpack%down_in = hydrodata_packed
                ! ! free memory immediately
                ! deallocate(hydrodata_packed)
            ! end if
        end if

        ! if (any(vn_WS_snowpack_lat .ne. "none")) then
            ! !four vertical directions
            ! if (vn_WS_snowpack_lat(1) .ne. "none")  then
                ! varName = vn_WS_snowpack_lat(1)
                ! if (isnc) then
                    ! call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                ! else
                    ! call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                ! end if
                ! L1_WSsnowpack%lateral_in = hydrodata_packed
                ! ! free memory immediately
                ! deallocate(hydrodata_packed)
            ! end if
            ! if (vn_WS_snowpack_lat(2) .ne. "none")  then
                ! varName = vn_WS_snowpack_lat(2)
                ! if (isnc) then
                    ! call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                ! else
                    ! call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                ! end if
                ! L1_WSsnowpack%lateral_out = hydrodata_packed
                ! ! free memory immediately
                ! deallocate(hydrodata_packed)
            ! end if
        ! end if
    !vn_WS_ponding
        if (vn_WS_ponding .ne. "none") then
            varName = vn_WS_ponding
            if (isnc) then
                call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
            else
                call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
            end if
            L1_WSponding%reservoir = hydrodata_packed * inputfluxUT
            ! free memory immediately
            deallocate(hydrodata_packed)
        end if
        if (any(vn_WS_ponding_ver .ne. "none")) then
            !four vertical directions
            if (vn_WS_ponding_ver(1) .ne. "none")  then
                varName = vn_WS_ponding_ver(1)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_WSponding%up_in = hydrodata_packed * inputfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)
            end if
            if (vn_WS_ponding_ver(2) .ne. "none")  then
                varName = vn_WS_ponding_ver(2)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_WSponding%down_out = hydrodata_packed * inputfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)
            end if
            if (vn_WS_ponding_ver(3) .ne. "none")  then
                varName = vn_WS_ponding_ver(3)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_WSponding%up_out = hydrodata_packed * inputfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)
            end if
            if (vn_WS_ponding_ver(4) .ne. "none")  then
                varName = vn_WS_ponding_ver(4)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_WSponding%down_in = hydrodata_packed * inputfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)
            end if
        end if

        ! if (any(vn_WS_snowpack_lat .ne. "none")) then
            ! !four vertical directions
            ! if (vn_WS_snowpack_lat(1) .ne. "none")  then
                ! varName = vn_WS_snowpack_lat(1)
                ! if (isnc) then
                    ! call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                ! else
                    ! call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                ! end if
                ! L1_WSponding%lateral_in = hydrodata_packed
                ! ! free memory immediately
                ! deallocate(hydrodata_packed)
            ! end if
            ! if (vn_WS_snowpack_lat(2) .ne. "none")  then
                ! varName = vn_WS_snowpack_lat(2)
                ! if (isnc) then
                    ! call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                ! else
                    ! call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                ! end if
                ! L1_WSponding%lateral_out = hydrodata_packed
                ! ! free memory immediately
                ! deallocate(hydrodata_packed)
            ! end if
        ! end if
    !vn_WS_surfaceflow
        if (vn_WS_surfaceflow .ne. "none") then
            varName = vn_WS_surfaceflow
            if (isnc) then
                call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
            else
                call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
            end if
            L1_WSsurface%reservoir = hydrodata_packed * inputfluxUT
            ! free memory immediately
            deallocate(hydrodata_packed)
        end if
        if (any(vn_WS_surfaceflow_ver .ne. "none")) then
            !four vertical directions
            if (vn_WS_surfaceflow_ver(1) .ne. "none")  then
                varName = vn_WS_surfaceflow_ver(1)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_WSsurface%up_in = hydrodata_packed * inputfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)
            end if
            if (vn_WS_surfaceflow_ver(2) .ne. "none")  then
                varName = vn_WS_surfaceflow_ver(2)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_WSsurface%down_out = hydrodata_packed * inputfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)
            end if
            if (vn_WS_surfaceflow_ver(3) .ne. "none")  then
                varName = vn_WS_surfaceflow_ver(3)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_WSsurface%up_out = hydrodata_packed * inputfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)
            end if
            if (vn_WS_surfaceflow_ver(4) .ne. "none")  then
                varName = vn_WS_surfaceflow_ver(4)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_WSsurface%down_in = hydrodata_packed * inputfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)
            end if
        end if
        if (any(vn_WS_surfaceflow_lat .ne. "none")) then
            !four vertical directions
            if (vn_WS_surfaceflow_lat(1) .ne. "none")  then
                varName = vn_WS_surfaceflow_lat(1)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_WSsurface%lateral_in = hydrodata_packed * inputfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)
            end if
            if (vn_WS_surfaceflow_lat(2) .ne. "none")  then
                varName = vn_WS_surfaceflow_lat(2)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_WSsurface%lateral_out = hydrodata_packed * inputfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)
            end if
        end if
    !vn_WS_baseflow
        if (vn_WS_baseflow .ne. "none") then
            varName = vn_WS_baseflow
            if (isnc) then
                call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
            else
                call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
            end if
            L1_WSbase%reservoir = hydrodata_packed * inputfluxUT
            ! free memory immediately
            deallocate(hydrodata_packed) 
        end if
        if (any(vn_WS_baseflow_ver .ne. "none")) then
            !four vertical directions
            if (vn_WS_baseflow_ver(1) .ne. "none")  then
                varName = vn_WS_baseflow_ver(1)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_WSbase%up_in = hydrodata_packed * inputfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)
            end if
            if (vn_WS_baseflow_ver(2) .ne. "none") then
                varName = vn_WS_baseflow_ver(2)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_WSbase%down_out = hydrodata_packed * inputfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)
            end if
            if (vn_WS_baseflow_ver(3) .ne. "none") then
                varName = vn_WS_baseflow_ver(3)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_WSbase%up_out = hydrodata_packed * inputfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)
            end if
            if (vn_WS_baseflow_ver(4) .ne. "none") then
                varName = vn_WS_baseflow_ver(4)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_WSbase%down_in = hydrodata_packed * inputfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)
            end if
        end if

        if (any(vn_WS_baseflow_lat .ne. "none")) then
            !four vertical directions
            if (vn_WS_baseflow_lat(1) .ne. "none")  then
                varName = vn_WS_baseflow_lat(1)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_WSbase%lateral_in = hydrodata_packed * inputfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)
            end if
            if (vn_WS_baseflow_lat(2) .ne. "none") then
                varName = vn_WS_baseflow_lat(2)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_WSbase%lateral_out = hydrodata_packed * inputfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)
            end if
        end if

    !vn_WS_subsurfaceflow
      if ((nSubsurfaceflow >= 1) .and. any(vn_WS_subsurfaceflow .ne. "none")) then
            do iver =1, nSubsurfaceflow
                varName = vn_WS_subsurfaceflow(iver)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_WSsubsurface(iver)%reservoir = hydrodata_packed * inputfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)

                if (any(vn_WS_subsurfaceflow_ver(iver,:) .ne. "none")) then
                    !four vertical directions
                    if (vn_WS_subsurfaceflow_ver(iver,1)  .ne. "none")  then
                        varName = vn_WS_subsurfaceflow_ver(iver,1)
                        if (isnc) then
                            call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                        else
                            call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                        end if
                        L1_WSsubsurface(iver)%up_in = hydrodata_packed * inputfluxUT
                        ! free memory immediately
                        deallocate(hydrodata_packed)
                    end if

                    if (vn_WS_subsurfaceflow_ver(iver,2)  .ne. "none") then
                        varName = vn_WS_subsurfaceflow_ver(iver,2)
                        if (isnc) then
                            call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                        else
                            call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                        end if
                        L1_WSsubsurface(iver)%down_out = hydrodata_packed * inputfluxUT
                        ! free memory immediately
                        deallocate(hydrodata_packed)
                    end if
                    if (vn_WS_subsurfaceflow_ver(iver,3)  .ne. "none") then
                        varName = vn_WS_subsurfaceflow_ver(iver,3)
                        if (isnc) then
                            call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                        else
                            call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                        end if
                        L1_WSsubsurface(iver)%up_out = hydrodata_packed * inputfluxUT
                        ! free memory immediately
                        deallocate(hydrodata_packed)
                    end if
                    if (vn_WS_subsurfaceflow_ver(iver,4)  .ne. "none")  then
                        varName = vn_WS_subsurfaceflow_ver(iver,4)
                        if (isnc) then
                            call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                        else
                            call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                        end if
                        L1_WSsubsurface(iver)%down_in = hydrodata_packed * inputfluxUT
                        ! free memory immediately
                        deallocate(hydrodata_packed)
                    end if
                end if

                if (any(vn_WS_subsurfaceflow_lat(iver,:) .ne. "none")) then
                    !two lateral directions
                    if (vn_WS_subsurfaceflow_lat(iver,1) .ne. "none")  then
                        varName = vn_WS_subsurfaceflow_lat(iver,1)
                        if (isnc) then
                            call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                        else
                            call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                        end if
                        L1_WSsubsurface(iver)%lateral_in = hydrodata_packed * inputfluxUT
                        ! free memory immediately
                        deallocate(hydrodata_packed)
                    end if
                    if (vn_WS_subsurfaceflow_lat(iver,2) .ne. "none")  then
                        varName = vn_WS_subsurfaceflow_lat(iver,2)
                        if (isnc) then
                            call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                        else
                            call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                        end if
                        L1_WSsubsurface(iver)%lateral_out = hydrodata_packed * inputfluxUT
                        ! free memory immediately
                        deallocate(hydrodata_packed)
                    end if
                end if
            end do
        end if

    !vn_WS_soil

      if ((nSoilLayers >= 1) .and. any(vn_WS_soil .ne. "none")) then
            do iver =1, nSoilLayers
                varName = vn_WS_soil(iver)
                if (isnc) then
                    call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                end if
                L1_WSsoil(iver)%reservoir = hydrodata_packed * inputfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)

                if (any(vn_WS_soil_ver(iver,:) .ne. "none")) then
                    !four vertical directions
                    if (vn_WS_soil_ver(iver,1)  .ne. "none")  then
                        varName = vn_WS_soil_ver(iver,1)
                        if (isnc) then
                            call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                        else
                            call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                        end if
                        L1_WSsoil(iver)%up_in = hydrodata_packed * inputfluxUT
                        ! free memory immediately
                        deallocate(hydrodata_packed)
                    end if

                    if (vn_WS_soil_ver(iver,2)  .ne. "none") then
                        varName = vn_WS_soil_ver(iver,2)
                        if (isnc) then
                            call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                        else
                            call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                        end if
                        L1_WSsoil(iver)%down_out = hydrodata_packed * inputfluxUT
                        ! free memory immediately
                        deallocate(hydrodata_packed)
                    end if

                    if (vn_WS_soil_ver(iver,3)  .ne. "none") then
                        varName = vn_WS_soil_ver(iver,3)
                        if (isnc) then
                            call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                        else
                            call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                        end if
                        L1_WSsoil(iver)%up_out = hydrodata_packed * inputfluxUT
                        ! free memory immediately
                        deallocate(hydrodata_packed)
                    end if
                    if (vn_WS_soil_ver(iver,4)  .ne. "none")  then
                        varName = vn_WS_soil_ver(iver,4)
                        if (isnc) then
                            call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                        else
                            call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                        end if
                        L1_WSsoil(iver)%down_in = hydrodata_packed * inputfluxUT
                        ! free memory immediately
                        deallocate(hydrodata_packed)
                    end if
                end if
                ! if (any(vn_WS_soil_lat(iver,:) .ne. "none")) then
                    ! !two lateral directions
                    ! if (vn_WS_soil_lat(iver,1) .ne. "none")  then
                        ! varName = vn_WS_soil_lat(iver,1)
                        ! if (isnc) then
                            ! call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                        ! else
                            ! call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                        ! end if
                        ! L1_WSsoil(iver)%lateral_in = hydrodata_packed
                        ! ! free memory immediately
                        ! deallocate(hydrodata_packed)
                    ! end if
                    ! if (vn_WS_soil_lat(iver,2) .ne. "none")  then
                        ! varName = vn_WS_soil_lat(iver,2)
                        ! if (isnc) then
                            ! call read_hydro_nc(fName, varName, simPerL, hydrodata_packed) 
                        ! else
                            ! call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed)
                        ! end if
                        ! L1_WSsoil(iver)%lateral_out = hydrodata_packed
                        ! ! free memory immediately
                        ! deallocate(hydrodata_packed)
                    ! end if
                ! end if
            end do
        end if
    !channel
    if (inputFormat_rout .eq. "nc") then
        fName = trim(dir_hydroInput) // trim(inputhydro_rout_filename) // '.nc'
        isnc = .TRUE.
    else if (inputFormat_rout .eq. "bin") then
        fdir = dir_hydroInput
        isnc = .FALSE.
    else
        call message()
        call message('***ERROR: the hydrological input files should be in either netCDF or binary format!')
        stop
    end if

        if (vn_chan_store .ne. "none") then
            varName = vn_chan_store
            if (isnc) then
                call read_hydro_rout_nc(fName, varName, simPerL, hydrodata_packed) 
            else
                call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed,datatype=2_i4)
            end if
            L1_WSchannel = hydrodata_packed * input_rout_fluxUT
            ! free memory immediately
            deallocate(hydrodata_packed)
        end if

        if (vn_chan_in .ne. "none") then
            varName = vn_chan_in
            if (isnc) then
                call read_hydro_rout_nc(fName, varName, simPerL, hydrodata_packed) 
            else
                call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed,datatype=2_i4)
            end if
            L1_WSchannel_upin = hydrodata_packed * input_rout_fluxUT
            ! free memory immediately
            deallocate(hydrodata_packed)
        end if

        if (vn_chan_out .ne. "none") then
            varName = vn_chan_out
            if (isnc) then
                call read_hydro_rout_nc(fName, varName, simPerL, hydrodata_packed) 
            else
                call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed,datatype=2_i4)
            end if
            L1_WSchannel_downout = hydrodata_packed * input_rout_fluxUT
            ! free memory immediately
            deallocate(hydrodata_packed)
        end if

        if (vn_srftochan .ne. "none") then
            varName = vn_srftochan
            if (isnc) then
                call read_hydro_rout_nc(fName, varName, simPerL, hydrodata_packed) 
            else
                call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed,datatype=2_i4)
            end if
            L1_srfTOchan = hydrodata_packed * inputfluxUT
            ! free memory immediately
            deallocate(hydrodata_packed)
        end if
 
        if (vn_baseftochan .ne. "none") then
            varName = vn_baseftochan
            if (isnc) then
                call read_hydro_rout_nc(fName, varName, simPerL, hydrodata_packed) 
            else
                call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed,datatype=2_i4)
            end if
            L1_basefTOchan = hydrodata_packed * inputfluxUT
            ! free memory immediately
            deallocate(hydrodata_packed)
        end if

        if ((nSubsrfToChannel >= 1) .and. any(vn_subsrftochan .ne. "none")) then
            do iver =1, nSubsrfToChannel
                varName = vn_subsrftochan(iver)
                if (isnc) then
                    call read_hydro_rout_nc(fName, varName, simPerL, hydrodata_packed) 
                else
                    call read_hydro_bin(fdir, varName, simPerL, hydrodata_packed,datatype=2_i4)
                end if
                L1_subsrfTOchan(:,:,iver) = hydrodata_packed * inputfluxUT
                ! free memory immediately
                deallocate(hydrodata_packed)
            end do
        end if

  end subroutine wqm_readhydroinput

  ! ------------------------------------------------------------------

  !     NAME
  !         read_hydro_nc

  !     PURPOSE
  !>        \brief Reads forcing input in NetCDF file format.

  subroutine read_hydro_nc(fName, varName, simPerL, hydrodata_packed, check)

    use mo_kind,             only: i4, dp
    use mo_message,          only: message
    use mo_utils,            only: eq
    use mo_string_utils,      only: num2str
    use mo_wqm_global_variables,    only: model_type, &
         m_Level,Model_mask,  inputcheck,L1_nCells, &
         L1_nSubs, L1_nHRUs
    use mo_wqm_constants,       only: nodata_dp
    use mo_wqm_netcdfinoutputs,  only: ncfile_getDim, ncfile_read
    implicit none

    character(len=*),                                  intent(in)  :: fName     !folder/file name where data are stored
    character(len=*),                                  intent(in)  :: varName   !variable name
    integer(i4),                                       intent(in)  :: simPerL   !length of simulation
    real(dp), dimension(:,:), allocatable,             intent(out) :: hydrodata_packed !data
    logical,                                 optional, intent(in)  :: check    ! .TRUE. if check for nodata values deactivated
    !                                                                            ! default = .FALSE. - check is done	 

    !local 
    real(dp), dimension(:,:,:), allocatable :: hydrodata !data
    integer(i4), dimension(5) :: dimen        ! dimension for NetCDF file
    !integer(i4)               :: dttype     ! datatype of attribute
    !character(256)            :: AttValues    ! netcdf attribute values
    !real(dp)                  :: nodata_value ! data nodata value
    integer(i4)               :: ncells

    integer(i4)               :: t
    inputcheck = .FALSE. !global variables default .FALSE. since it might take long time

    if (present(check)) inputcheck = check

    ! get dimensions -dim1=ncols,dim2=nrows,dim3=timesteps because netCDF is row-major ordered

    dimen = ncfile_getDim(trim(fName), trim(varName))
    !fully distributed
    if (model_type == 1) then
        if ( (dimen(1) .ne. m_Level%ncols) .or. (dimen(2) .ne. m_Level%nrows) ) then
          stop '***ERROR: hydro input (nc) not match model domain dimensions'
        end if
        ncells = L1_nCells
    !semi-distributed
    else if (model_type == 2) then
        if ( (dimen(1) .ne.L1_nHRUs) .or. (dimen(2) .ne.  L1_nSubs) ) then
          stop '***ERROR: hydro input (nc) not match model domain dimensions'
        end if
        ncells = L1_nCells
    end if
    ! determine no data value
    !call Get_NcVarAtt(fName, varName, '_FillValue', AttValues, dtype=dttype)
    ! convert to number
    !read(AttValues, *) nodata_value

    ! Check if time steps in file cover simulation period
    if (simPerL > dimen(3)) then
       call message('***ERROR: time series length of input data: ', trim(varName), &
            '          is not matching modelling period.')
       stop
    end if

    ! alloc and read
    allocate(hydrodata(dimen(1), dimen(2), simPerL))
    call ncfile_read(trim(fName), trim(varName), hydrodata, &
              nstart = (/ 1_i4, 1_i4, 1_i4 /), &
              ncount = (/ dimen(1), dimen(2), simPerL /) )
       
    ! start checking values
    ! if (inputcheck .and. model_type == 1) then
        ! do t = 1, dimen(3)
       ! ! neglect checking for naodata values if optional nocheck is given
          ! if (any(eq(hydrodata(:,:,t),nodata_value) .and. (Model_mask))) then
             ! call message('***ERROR: Hydrological input (nc): nodata value within basin ')
             ! call message('          boundary in variable: ', trim(varName))
             ! call message('          at timestep         : ', trim(num2str(t)))
             ! stop
          ! end if
        ! end do
    ! end if

    if (model_type ==1) then
      allocate(hydrodata_packed(ncells, simPerL))
      hydrodata_packed = nodata_dp
      do t = 1, simPerL
        hydrodata_packed(:,t) = pack( hydrodata(:,:,t), MASK=Model_mask(:,:) ) 
      end do
    else if (model_type ==2) then
      !ncalc = L1_nHRUs*L1_nSubs
      allocate(hydrodata_packed(ncells, simPerL))
      hydrodata_packed = nodata_dp
      do t = 1, simPerL
        hydrodata_packed(:,t) = reshape(hydrodata(:,:,t),(/ncells/)) !loop first HRUs then SUBs
      end do
    end if

  end subroutine read_hydro_nc

  ! ------------------------------------------------------------------

  !     NAME
  !         read_hydro_rout_nc

  !     PURPOSE
  !>        \brief Reads forcing input in NetCDF file format.

  subroutine read_hydro_rout_nc(fName, varName, simPerL, hydrodata_packed, check)

    use mo_kind,             only: i4, dp
    use mo_message,          only: message
    use mo_utils,            only: eq
    use mo_string_utils,      only: num2str
    use mo_wqm_global_variables,    only: model_type, &
         r_Level,Rout_mask,  inputcheck,L1_nLinks, &
         L1_nSubs
    use mo_wqm_constants,       only: nodata_dp
    use mo_wqm_netcdfinoutputs,  only: ncfile_getDim, ncfile_read
    implicit none

    character(len=*),                                  intent(in)  :: fName     !folder/file name where data are stored
    character(len=*),                                  intent(in)  :: varName   !variable name
    integer(i4),                                       intent(in)  :: simPerL   !length of simulation
    real(dp), dimension(:,:), allocatable,             intent(out) :: hydrodata_packed !data
    logical,                                 optional, intent(in)  :: check    ! .TRUE. if check for nodata values deactivated
    !                                                                            ! default = .FALSE. - check is done	 

    !local 
    real(dp), dimension(:,:,:), allocatable :: hydrodata !data
    real(dp), dimension(:,:), allocatable :: hydrodata2D !data
    integer(i4), dimension(5) :: dimen        ! dimension for NetCDF file
    ! integer(i4)               :: dttype     ! datatype of attribute
    ! character(256)            :: AttValues    ! netcdf attribute values
    ! real(dp)                  :: nodata_value ! data nodata value
    integer(i4)               :: ncells
    integer(i4)               :: t
    
    inputcheck = .FALSE. !global variables default .FALSE. since it might take long time

    if (present(check)) inputcheck = check

    ! get dimensions -dim1=ncols,dim2=nrows,dim3=timesteps because netCDF is row-major ordered
    dimen = ncfile_getDim(trim(fName), trim(varName))
  !fully distributed
  if (model_type == 1) then
    if ( (dimen(1) .ne. r_Level%ncols) .or. (dimen(2) .ne. r_Level%nrows) ) then
        stop '***ERROR: hydro input (nc) not match routing domain dimensions'
    end if
    ncells = L1_nLinks
    ! ! determine no data value
    ! call Get_NcVarAtt(fName, varName, '_FillValue', AttValues, dtype=dttype)
    ! ! convert to number
    ! read(AttValues, *) nodata_value
    ! Check if time steps in file cover simulation period
    if (simPerL > dimen(3)) then
       call message('***ERROR: time series length of input data: ', trim(varName), &
            '          is not matching modelling period.')
       stop
    end if

    ! alloc and read
    allocate(hydrodata(dimen(1), dimen(2), simPerL))
    call ncfile_read(trim(fName), trim(varName), hydrodata, &
              nstart = (/ 1_i4, 1_i4, 1_i4 /), &
              ncount = (/ dimen(1), dimen(2), simPerL /) )       
    ! ! start checking values
    ! if (inputcheck .and. model_type == 1) then
        ! do t = 1, dimen(3)
       ! ! neglect checking for naodata values if optional nocheck is given
          ! if (any(eq(hydrodata(:,:,t),nodata_value) .and. (Rout_mask))) then
             ! call message('***ERROR: Hydrological input (nc): nodata value within basin ')
             ! call message('          boundary in variable: ', trim(varName))
             ! call message('          at timestep         : ', trim(num2str(t)))
             ! stop
          ! end if
        ! end do
    ! end if
    !prepare output
    allocate(hydrodata_packed(ncells, simPerL))
    hydrodata_packed = nodata_dp
    do t = 1, simPerL
      hydrodata_packed(:,t) = pack( hydrodata(:,:,t), MASK=Rout_mask(:,:) ) 
    end do

  !semi-distributed
  else if (model_type == 2) then
    if ((dimen(1) .ne. L1_nSubs) ) then
        stop '***ERROR: hydro input (nc) not match routing domain dimensions'
    end if
    ncells = L1_nSubs
    ! ! determine no data value
    ! call Get_NcVarAtt(fName, varName, '_FillValue', AttValues, dtype=dttype)
    ! ! convert to number
    ! read(AttValues, *) nodata_value
    ! Check if time steps in file cover simulation period
    if (simPerL > dimen(2)) then
       call message('***ERROR: time series length of input data: ', trim(varName), &
            '          is not matching modelling period.')
       stop
    end if

    ! alloc and read
    allocate(hydrodata2D(dimen(1), simPerL))
    call ncfile_read(trim(fName), trim(varName), hydrodata2D, &
              nstart = (/ 1_i4, 1_i4 /), &
              ncount = (/ dimen(1), simPerL/) )
    !ncalc = L1_nHRUs*L1_nSubs
    allocate(hydrodata_packed(ncells, simPerL))
    hydrodata_packed = nodata_dp
    hydrodata_packed(:,:) = hydrodata2D(:,:)

  end if

  end subroutine read_hydro_rout_nc
  ! ------------------------------------------------------------------

  !     NAME
  !         read_hydro_bin

  !> 

  !     PURPOSE
  !>        \brief Reads forcing input in binary file format (double precision).
  !>        (1) the each input variable should be stored in a separated file 
  !>            named as varName.bin
  !         (2) values should be organized as 2-D: (number of grid/calculate unit,simulation length)
  !             corresponding to the grid ids/unit ids from the model domain file  


  subroutine read_hydro_bin(fdir, varName, simPerL, hydrodata, datatype)

    use mo_kind,             only: i4, dp
    use mo_message,          only: message

     use mo_wqm_global_variables,    only: &
         L1_nCells, L1_nLinks

    implicit none

    character(len=*),                                  intent(in)  :: fdir     !folder/file name where data are stored
    character(len=*),                                  intent(in)  :: varName   !variable name
    integer(i4),                                       intent(in)  :: simPerL   !length of simulation
    real(dp), dimension(:,:), allocatable,             intent(inout) :: hydrodata!data
    integer(i4),                  optional,            intent(in)  :: datatype ! 1 for modeling level, 2 for routing level

    !local 
    character(256)            :: fName
    integer(i4)               :: ii, checking, ncells,idatatype
    real(dp), dimension(:) ,allocatable   :: tmpdata

    if (present(datatype)) then
      idatatype = datatype
    else
      idatatype = 1_i4
    end if
    if (idatatype ==1_i4) then
      ncells = L1_nCells
    else if (idatatype ==2_i4) then
      ncells = L1_nLinks
    end if

    fName = trim(fdir) // trim(varName) // '.bin'
    open(unit=89, file=trim(fName), &
        form='unformatted', access='direct', recl=8*ncells)  !**double precision**

    allocate(tmpdata(ncells))
    allocate(hydrodata(ncells,simPerL))

    ii = 1_i4
    do while(.TRUE.)
        read(89,rec=ii, iostat = checking) tmpdata
        if (checking < 0) then
            if (ii .ne. simPerL) then
                call message('***ERROR: time series length of input data: ', trim(varName), &
                         '          is not matching modelling period.')
                stop
            end if
            exit
        end if
        if (ii <= simPerL) then
            hydrodata(:,ii) = tmpdata
        end if
        ii = ii + 1
    end do

    close(89)

  end subroutine read_hydro_bin



end module mo_wqm_readhydroinput

