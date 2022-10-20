!> \file mo_wqm_readconfig.f90

!> \brief Reads water quality configure, input files and allocating WQ global variables.

!> \details This module is for water quality modeling.\n
!> 

!> \authors Xiaoqiang Yang
!> \date Sep 2021

MODULE mo_wqm_readconfig


  USE mo_kind, ONLY: i4, sp, dp

  
  IMPLICIT NONE

  PUBLIC :: wqm_readconfig               ! read water qualtiy configure 


  ! ------------------------------------------------------------------

 CONTAINS
  !     NAME
  !         wqm_readconfig

  !     PURPOSE
  !>        \brief .

  !>        \details 

  !     CALLING SEQUENCE
  !         



  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Xiaoqiang Yang
  !>        \date Sep 2021

 subroutine wqm_readconfig()

  use mo_file,                   only: file_namelist_wqm, unamelist, &
                                       file_param_wqm, unamelist_param, &
                                       file_stateoutput_wqm,ustatOutput
  use mo_nml,                    only: open_nml, position_nml, close_nml
  use mo_julian,                 only: dec2date, date2dec
  use mo_wqm_global_variables,   only: hydroModel, hydromodel_info,  &
         timestep, nTstepDay, resolution, rout_resolution, model_type, &   !main configuration from wqm_config.nmml
         warmingDays, evalPer, simPer,   &
       !directories for WQM and coupling
         dir_catchinfo, dir_gauges, dir_hydroInput, dir_Out, &
       !file names of catchment information
         fn_catchDomain, fdir_type,fn_flowDir,fn_spatialorganization, fn_channel_mask, fn_channel_length, fn_channel_width, &
         fn_cropinfo, fn_lut_rotation, fn_lut_monLAI, &
         fn_initial_condition,fn_globalRadiation, const_crain,const_csnowmelt, &
         fn_init_concIN,fn_init_concON, &
       !geophysical related information and file names
         nLanduseTotal, fn_LanduseTypes, nLUperGrid, fn_flanduse, &
         nSoilTotal,fn_SoilTypes, nSoilperGrid, fn_fsoil, &
         nRotationTotal,fn_RotationTypes, nRotationperGrid, fn_frotation, &
         nSoilLayers, fn_soilDepth, fn_wiltpoint, fn_satMoist,  &
         fn_classProperties, &
       !gauge network
         fn_gauge_loc, basin_gauge, nHydroGauges,nWQGauges, nInflowGauges,  &
         hydro_eval_gaugenr,hydro_evalgaugeIDs,wqm_eval_gaugenr,wqm_evalgaugeIDs, &
       !process definition
         processCase,vn_precipitation,vn_airTemperature,vn_WS_interception,vn_evapCanopy, vn_throughfall, &
         vn_WS_snowpack, vn_snowmelt,vn_infiltration,vn_WS_ponding, &
         ETfluxUT,nEvapLayer, vn_ETevap, nTranspLayer, vn_ETtransp, &
         vn_surfaceflow,vn_WS_surfaceflow,vn_fromWS_surfaceflow, SrfOVFmixratio, SrfOVFmixratio_vn, &
         nSubsurfaceflow,vn_subsurfaceflow,vn_WS_subsurfaceflow,vn_fromWS_subsurfaceflow, &
         vn_baseflow, vn_WS_baseflow,ratio_RetentionWS, vn_fromWS_baseflow, &
         vn_directflow,vn_Esealed,vn_WS_directflow, Sealed_LanduseNrs,Sealed_LanduseIDs,&
       !hydro structure
         inputFormat, inputcheck, inputhydro_filename, inputfluxUT, vn_soilTemperature, &
         vn_WS_interception_ver,vn_WS_snowpack_ver, vn_WS_ponding_ver,vn_WS_surfaceflow_ver, &
         vn_WS_subsurfaceflow_ver, vn_WS_baseflow_ver, &
         !vn_WS_interception_lat,vn_WS_snowpack_lat,vn_WS_ponding_lat,   &
         vn_WS_surfaceflow_lat,vn_WS_subsurfaceflow_lat,vn_WS_baseflow_lat,   &
         vn_WS_soil,vn_SatSM, vn_WP, vn_WS_soil_ver, &   !vn_WS_soil_lat,
         vn_reinfiltration,vn_reinfil_soil, & !if re-infiltration is allowed
         all_vn_WS, &
       !channel
         inputFormat_rout, inputhydro_rout_filename,input_rout_fluxUT,  &
         vn_chan_store, vn_chan_in, vn_chan_out, &
         vn_srftochan,vn_baseftochan, &
         nSubsrfToChannel, vn_subsrftochan,  & 
         !wqm parameters
         nTotalParam, wqm_param,   &
         !type variable
         wqm_parameterization, period, &
		 !nc output
         outputFlxState_wqm, timeStep_model_outputs_wqm, nOUTstate_wqm
  !mhm optimization implementations
  use mo_opti_variables,         only: &
         optimize,optimize_restart, opti_method, opti_function, &
         nIterations, seed, dds_r, sa_temp, sce_ngs, sce_npg, sce_nps, mcmc_opti, mcmc_error_params

  use mo_wqm_constants,          only: nodata_i4, maxnLandUse,maxnSoilType,maxnRotation,maxnSoilLayer, &
                                       maxnhydroGauges, maxnWQGauges, maxnInflowGauges, &
                                       maxnSubSrfstore,maxnChannelSeepageConponent,maxnNtransfors
  use mo_string_utils,           only: num2str
  use mo_message,                only: message 
  use mo_append,                 only: append
	   
  !local
    type(period)                                      :: sim_Per, eval_Per  !evaluation period
    !file names
    character(256), dimension(maxnLandUse)            :: fn_arealshare_LU       ! filename land use for each type
    character(256), dimension(maxnSoilType)           :: fn_arealshare_soil     ! filename soil type
    character(256), dimension(maxnRotation)           :: fn_arealshare_rotation ! filename rotation type
    character(256), dimension(maxnSoilLayer)          :: fn_soil_Depth          ! filename soil depth for each layer (lowerdepth - upperdepth)
    character(256), dimension(maxnSoilLayer)          :: fn_soil_wiltpoint      ! filename soil parameter wilting point
    character(256), dimension(maxnSoilLayer)          :: fn_soil_satMoist       ! filename soil parameter saturation soil moisture 
    !gauge network
    character(256), dimension(maxnhydroGauges)              :: hydro_gauge_filename
    character(256), dimension(maxnWQGauges)                 :: wqm_evalgauge_filename
    character(256), dimension(maxnInflowGauges)             :: InflowGauge_filename
    character(256), dimension(maxnInflowGauges)             :: wqm_inflowgauge_filename
    integer(i4), dimension(maxnhydroGauges)            :: gauge_id,hydro_eval_gaugeIDs
    integer(i4), dimension(maxnWQGauges)               :: wqm_gauge_id,wqm_eval_gaugeIDs
    integer(i4), dimension(maxnInflowGauges)           :: InflowGauge_id
    logical,     dimension(maxnInflowGauges)           :: InflowisHead
 
    !variable names
    character(256), dimension(maxnSoilLayer)          :: vn_evap
    character(256), dimension(maxnSoilLayer)          :: vn_transp
    character(256), dimension(maxnSoilLayer)          :: vn_reinfilsoil

    character(256), dimension(maxnSubSrfstore)        :: vn_subsurface_flow
    character(256), dimension(maxnSubSrfstore)        :: vn_WS_subsurface_flow
    character(256), dimension(maxnSubSrfstore)        :: vn_fromWS_subsurface_flow


    character(256), dimension(maxnSubSrfstore,4)        :: vn_WS_subsurface_flow_ver
    character(256), dimension(maxnSubSrfstore,2)        :: vn_WS_subsurface_flow_lat
    !subsurface seepage to channel should not be more than subsurface storage
    character(256), dimension(maxnChannelSeepageConponent)        :: vn_sub_srftochan
    !character(256), dimension(maxnChannelSeepageConponent)        :: vn_fromWS_sub_srftochan

    character(256), dimension(maxnSoilLayer)            ::vn_WSsoil
    character(256), dimension(maxnSoilLayer)            ::vn_Satsoilmoist
    character(256), dimension(maxnSoilLayer)            ::vn_Wiltpnt
    character(256), dimension(maxnSoilLayer,4)          ::vn_WSsoil_ver
    !character(256), dimension(maxnSoilLayer,2)          ::vn_WSsoil_lat

    integer(i4),    dimension(maxnLandUse)              :: Sealed_Landuse_IDs

    character(256)    :: vn_channelstorage
    character(256)    :: vn_upstreaminflow
    character(256)    :: vn_downstreamoutflow

    logical, dimension(99)             :: output_FlxState_wqm

  !parameter name
    logical        :: fexist
    type(wqm_parameterization), dimension(maxnNtransfors)   :: Nparameters

  !loc variables
    integer(i4)       :: idx,idgrp
    real(dp)       :: jday_frac

  !define namelist
  !namelist directories for water quality model
    namelist /mainconfig/ hydroModel, hydromodel_info, timestep, resolution, rout_resolution, model_type
  ! namelist for time settings
    namelist /time_periods/ sim_Per, eval_Per
    namelist /directories/ dir_catchinfo, dir_gauges, dir_hydroInput, dir_Out
  !file names for basic catchment information
    namelist /catchinfo/ fn_catchDomain,fdir_type, fn_flowDir, fn_spatialorganization, &
                         fn_channel_mask, fn_channel_length, fn_channel_width, &
                         fn_cropinfo, fn_lut_rotation, fn_lut_monLAI, fn_initial_condition, & !crop related information
                         fn_init_concIN,fn_init_concON, fn_globalRadiation,const_crain,const_csnowmelt, & !constant concentration in precipitation and snowmelt
                         nLanduseTotal,nLUperGrid,fn_LanduseTypes, fn_arealshare_LU, & !landuse 
                         nSoilTotal,nSoilperGrid,fn_SoilTypes,  fn_arealshare_soil, & !soil type
                         nSoilLayers, fn_soil_Depth, & ! soil layers/depth
                         fn_soil_wiltpoint,fn_soil_satMoist, & !soil parameters for each grid and layer
                         nRotationTotal,nRotationperGrid, fn_RotationTypes, fn_arealshare_rotation, &
                         fn_classProperties
    namelist /gaugenetwork/ fn_gauge_loc, nHydroGauges, gauge_id, hydro_gauge_filename,&
                            hydro_eval_gaugenr,hydro_eval_gaugeIDs, & !hydrological obs
                            nWQGauges, wqm_gauge_id, wqm_evalgauge_filename,&
                            wqm_eval_gaugenr, wqm_eval_gaugeIDs, & !WQ obs
                            nInflowGauges, InflowGauge_id,InflowisHead, InflowGauge_filename, wqm_inflowgauge_filename
    namelist /processes/ processCase, vn_precipitation,vn_airTemperature,vn_WS_interception,vn_evapCanopy,&
                         vn_throughfall,vn_WS_snowpack,vn_snowmelt,vn_infiltration,vn_WS_ponding,   &
                         ETfluxUT,nEvapLayer, vn_evap, nTranspLayer,vn_transp,&
                         vn_surfaceflow,vn_WS_surfaceflow,vn_fromWS_surfaceflow, SrfOVFmixratio, SrfOVFmixratio_vn, &
                         nSubsurfaceflow,vn_subsurface_flow,vn_WS_subsurface_flow, &
                         vn_fromWS_subsurface_flow,vn_baseflow, vn_WS_baseflow,ratio_RetentionWS,vn_fromWS_baseflow, &
                         vn_directflow,vn_Esealed,vn_WS_directflow,Sealed_LanduseNrs, Sealed_Landuse_IDs
    namelist /hydrology_structure/ inputFormat,inputcheck, inputhydro_filename, inputfluxUT, &
                                   vn_soilTemperature,vn_WS_interception_ver,vn_WS_snowpack_ver, vn_WS_ponding_ver, &
                                   !vn_WS_interception_lat,vn_WS_snowpack_lat, vn_WS_ponding_lat, &
                                   vn_WSsoil,vn_Satsoilmoist, vn_Wiltpnt, &
                                   vn_WSsoil_ver, vn_reinfiltration,vn_reinfilsoil, &
                                   vn_WS_surfaceflow_ver,vn_WS_subsurface_flow_ver, vn_WS_baseflow_ver, &
                                   vn_WS_surfaceflow_lat,vn_WS_subsurface_flow_lat,vn_WS_baseflow_lat

    namelist /channel_routing/ inputFormat_rout, inputhydro_rout_filename,input_rout_fluxUT, &
                               vn_channelstorage, vn_upstreaminflow, vn_downstreamoutflow, &
                               vn_srftochan,vn_baseftochan, &
                               nSubsrfToChannel, vn_sub_srftochan 
  !optimization options from mHM implementations
    namelist /optimization/ optimize,optimize_restart, opti_method, opti_function, &
                            nIterations, seed, dds_r, sa_temp, sce_ngs, sce_npg, sce_nps, mcmc_opti, mcmc_error_params
  !parameter namelist: wqm_parameter.nml
  !Nitrogen sub-model
    namelist /parameters/ nTotalParam, Nparameters
                          ! instream_denitrification,autotrophicuptk, netprimaryprod, &
                          ! soil_denitrification,degradationN, mineralisationN,&
                          ! dissolutionN

    namelist /statevarsoutput/timeStep_model_outputs_wqm,nOUTstate_wqm, output_FlxState_wqm 

  !=================================================
  !INITIALIZE VARIABLE NAMES
  !=================================================
    vn_precipitation="none"
    vn_airTemperature="none"
    vn_WS_interception="none"
    vn_evapCanopy="none"
    vn_throughfall="none"
    vn_WS_snowpack="none"
    vn_snowmelt="none"
    vn_infiltration="none"
    vn_WS_ponding="none"
    vn_evap="none"
    vn_transp="none"
    vn_soilTemperature="none"
    vn_surfaceflow="none"
    vn_WS_surfaceflow="none"
    vn_fromWS_surfaceflow="none"
    SrfOVFmixratio_vn = "none"
    vn_subsurface_flow="none"
    vn_WS_subsurface_flow="none"
    vn_fromWS_subsurface_flow="none"
    vn_baseflow="none"
    vn_WS_baseflow="none"
    vn_fromWS_baseflow="none"
    vn_WS_interception_ver="none"
    vn_WS_snowpack_ver="none"
    vn_WS_ponding_ver="none"
    vn_WSsoil="none"
    vn_Satsoilmoist="none"
    vn_Wiltpnt="none"
    vn_WSsoil_ver="none"
    vn_reinfiltration="none"
    vn_reinfilsoil="none"
    vn_WS_surfaceflow_ver="none"
    vn_WS_subsurface_flow_ver="none"
    vn_WS_baseflow_ver="none"
    vn_WS_surfaceflow_lat="none"
    vn_WS_subsurface_flow_lat="none"
    vn_WS_baseflow_lat="none"
    vn_channelstorage="none"
    vn_upstreaminflow="none"
    vn_downstreamoutflow="none"
    vn_srftochan ="none"
    vn_baseftochan ="none"
    vn_chan_store ="none"
    vn_chan_in ="none"
    vn_chan_out ="none"
    vn_sub_srftochan ="none"
  !=================================================
  !READ namelist "wqm_config.nml"
  !=================================================
    call open_nml(file_namelist_wqm, unamelist, quiet =.true.)
    !------------------------
    !read basic config info
    !------------------------
    call position_nml('mainconfig', unamelist)
    read(unamelist, nml=mainconfig)
    nTstepDay = 24_i4/timestep  !number of steps per day

    if (abs(rout_resolution - resolution) > 0.001_dp) then
       call message()
       call message('WARMING: modeling and routing resolution do not match')
       call message('terrestrial lateral flux exchange is therefore deactivated!') 
    end if
    !------------------------
    !  read simulation time periods incl. warming days
    !------------------------
    call position_nml('time_periods', unamelist)
    read(unamelist, nml=time_periods)
    evalPer = eval_Per
    simPer = sim_Per
    ! get julian day of evaluation period
    jday_frac = date2dec(dd=evalPer%dStart, mm=evalPer%mStart, yy=evalPer%yStart)
    evalPer%julStart = nint(jday_frac)
    jday_frac = date2dec(dd=evalPer%dEnd, mm=evalPer%mEnd, yy=evalPer%yEnd)
    evalPer%julEnd  = nint(jday_frac, i4 )

    jday_frac = date2dec(dd=simPer%dStart, mm=simPer%mStart, yy=simPer%yStart)
    simPer%julStart = nint(jday_frac)
    jday_frac = date2dec(dd=simPer%dEnd, mm=simPer%mEnd, yy=simPer%yEnd)
    simPer%julEnd  = nint(jday_frac, i4 )
    !get the warming days
    warmingDays = evalPer%julStart - simPer%julStart
    if (warmingDays == 0_i4) then
       call message()
       call message('WARMING: warming-up period not set, be cautious!')
    ! else if (warmingDays < 0_i4) then
       ! call message()
       ! call message('***ERROR: evaluation period is not allowed to be started earlier than simulation period')
       ! stop
    end if

    !------------------------
    !  read directories
    !------------------------
    call position_nml('directories', unamelist)
    read(unamelist, nml=directories)
    !------------------------
    !  read catchment information
    !------------------------
    call position_nml('catchinfo', unamelist)
    read(unamelist, nml=catchinfo)
    ! allocate 
    if (nLUperGrid > 1_i4) then
        allocate(fn_flanduse(nLUperGrid))
        fn_flanduse(:) =  fn_arealshare_LU(1:nLUperGrid)
    end if
    if (nSoilperGrid > 1_i4) then
        allocate(fn_fsoil(nSoilperGrid))
        fn_fsoil(:) = fn_arealshare_soil(1:nSoilperGrid)
    end if
    if (nRotationperGrid > 1_i4) then
        allocate(fn_frotation(nRotationperGrid))
        fn_frotation(:) = fn_arealshare_rotation(1:nRotationperGrid)
    end if

    allocate(fn_soilDepth(nSoilLayers))
    allocate(fn_wiltpoint(nSoilLayers))
    allocate(fn_satMoist(nSoilLayers))
    fn_soilDepth(1:nSoilLayers) = fn_soil_Depth(1:nSoilLayers)
    fn_wiltpoint(:) = fn_soil_wiltpoint(1:nSoilLayers)
    fn_satMoist(:) = fn_soil_satMoist(1:nSoilLayers)
    !------------------------
    !  read gauge network
    !------------------------
    call position_nml('gaugenetwork', unamelist)
    read(unamelist, nml=gaugenetwork)
    !allocation
    allocate(basin_gauge%fname_dis(max(1, nHydroGauges)))
    allocate(basin_gauge%fname_wqm(max(1, nWQGauges)))
    allocate(basin_gauge%gaugeid(max(1, nHydroGauges)))
    allocate(basin_gauge%gaugeid_loc(max(1, nHydroGauges)))
    allocate(basin_gauge%wqm_gaugeid(max(1, nWQGauges)))
    allocate(basin_gauge%wqm_gaugeid_loc(max(1, nWQGauges)))
    allocate(basin_gauge%Inflowfname_dis(max(1, nInflowGauges)))
    allocate(basin_gauge%Inflowfname_wqm(max(1, nInflowGauges)))
    allocate(basin_gauge%inflowgaugeid(max(1, nInflowGauges)))
    allocate(basin_gauge%inflow_ishead(max(1, nInflowGauges)))
    allocate(basin_gauge%inflowgaugeid_loc(max(1, nInflowGauges)))
    allocate(hydro_evalgaugeIDs(max(1,hydro_eval_gaugenr)))
    allocate(wqm_evalgaugeIDs(max(1,wqm_eval_gaugenr)))
    basin_gauge%gaugeid = nodata_i4
    basin_gauge%wqm_gaugeid = nodata_i4
    basin_gauge%inflowgaugeid = nodata_i4
    basin_gauge%gaugeid_loc = nodata_i4
    basin_gauge%wqm_gaugeid_loc = nodata_i4
    basin_gauge%inflowgaugeid_loc = nodata_i4
    basin_gauge%fname_dis =num2str(nodata_i4)
    basin_gauge%fname_wqm =num2str(nodata_i4)
    basin_gauge%Inflowfname_dis =num2str(nodata_i4)
    basin_gauge%Inflowfname_wqm =num2str(nodata_i4)
    hydro_evalgaugeIDs = nodata_i4
    wqm_evalgaugeIDs = nodata_i4

    basin_gauge%fname_dis(1:max(1, nHydroGauges)) = hydro_gauge_filename(1:max(1, nHydroGauges))
    basin_gauge%fname_wqm(1:max(1, nWQGauges)) = wqm_evalgauge_filename(1:max(1, nWQGauges))
    basin_gauge%Inflowfname_dis(1:max(1, nInflowGauges)) = InflowGauge_filename(1:max(1, nInflowGauges))
    basin_gauge%Inflowfname_wqm(1:max(1, nInflowGauges)) = wqm_inflowgauge_filename(1:max(1, nInflowGauges))

    basin_gauge%gaugeid(:) = gauge_id(1:max(1, nHydroGauges))
    basin_gauge%wqm_gaugeid(:) = wqm_gauge_id(1:max(1, nWQGauges))
    basin_gauge%gaugeid(:) = gauge_id(1:max(1, nHydroGauges))
    basin_gauge%inflowgaugeid(:) = InflowGauge_id(1:max(1, nInflowGauges))
    basin_gauge%inflow_ishead(:) = InflowisHead(1:max(1, nInflowGauges))
    hydro_evalgaugeIDs(:) = hydro_eval_gaugeIDs(1:max(1,hydro_eval_gaugenr))
    wqm_evalgaugeIDs(:) = wqm_eval_gaugeIDs(1:max(1,wqm_eval_gaugenr))

    !------------------------
    !  process and variable names
    !------------------------
    call position_nml('processes', unamelist)
    read(unamelist, nml=processes)
    !check mataches between processes and variables
    if (processCase(1) == 1_i4) then
      if ((vn_throughfall .eq. "none") .or. (vn_WS_interception .eq. "none")) then
        call message()
        call message('***ERROR: canopy interception activated, but storage or flux throughfall not provided')
        stop
      end if
      ! collect WS variable names
      call append(all_vn_WS, vn_WS_interception)
    end if


    if (processCase(2) == 1_i4) then
      if((vn_snowmelt .eq. "none") .or.(vn_WS_snowpack .eq. "none")) then
        call message()
        call message('***ERROR: snowpack and snowmelt activated, but storage or flux snowmelt not provided')
        stop
      end if
      ! collect WS variable names
      call append(all_vn_WS, vn_WS_snowpack)
    end if


    if (processCase(3) == 1_i4) then
      if ((vn_infiltration .eq. "none") .or. (vn_WS_ponding .eq. "none")) then
        call message()
        call message('***ERROR: surface ponding and infiltration activated, but not provided')
        stop
      end if
      ! collect WS variable names
      call append(all_vn_WS, vn_WS_ponding)
    end if


    if ((processCase(5) == 1_i4) .and. &
        all(vn_transp .eq. "none"))then
        call message()
        call message('***ERROR: Overall ET selected, but fluxes from soil layers not provided')
        stop
    else if ((processCase(5) == 2_i4)) then 
        if (all(vn_transp .eq. "none") .or. &
          all(vn_evap .eq. "none")) then
            call message()
            call message('***ERROR: Separated ET selected, but evap or transp fluxes from soil layers not provided')
            stop
        end if
        allocate(vn_ETevap(nSoilLayers))
        vn_ETevap(:) = "none"  !set default variable name as "none" for all soil layers, if not specified in configuration
        vn_ETevap(1:max(1,nEvapLayer)) = vn_evap(1:max(1,nEvapLayer))
    end if
    allocate(vn_ETtransp(nSoilLayers))
    vn_ETtransp(:) = "none" !set default variable name as "none" for all soil layers, if not specified in configuration
    vn_ETtransp(1:max(1,nTranspLayer)) = vn_transp(1:max(1,nTranspLayer))


    if (processCase(7) == 1_i4)  then  
       if ((vn_surfaceflow .eq. "none") .or. &
           (vn_WS_surfaceflow .eq. "none") .or. &
           (vn_fromWS_surfaceflow .eq. "none")) then
          call message()
          call message('***ERROR: surface flow activated, but surfaceflow or surface Reservoir not provided')
          stop
       end if
       ! collect WS variable names
       call append(all_vn_WS, vn_WS_surfaceflow)
    end if

    if (processCase(8) == 1_i4) then
       if (all(vn_subsurface_flow .eq. "none") .or. &
           all(vn_WS_subsurface_flow .eq. "none") .or. &
           all(vn_fromWS_subsurface_flow .eq. "none")) then
          call message()
          call message('***ERROR: subsurface flow activated, but subsurfaceflow or subsurface Reservoir not provided')
          stop
       end if
       allocate(vn_subsurfaceflow(max(1,nSubsurfaceflow)))
       allocate(vn_WS_subsurfaceflow(max(1,nSubsurfaceflow)))
       allocate(vn_fromWS_subsurfaceflow(max(1,nSubsurfaceflow)))
       vn_subsurfaceflow = "none"
       vn_WS_subsurfaceflow ="none"
       vn_fromWS_subsurfaceflow="none"
       vn_subsurfaceflow(:) = vn_subsurface_flow(1:max(1,nSubsurfaceflow))
       vn_WS_subsurfaceflow(:) = vn_WS_subsurface_flow(1:max(1,nSubsurfaceflow))
       vn_fromWS_subsurfaceflow(:) = vn_fromWS_subsurface_flow(1:max(1,nSubsurfaceflow))
       ! collect WS variable names
       call append(all_vn_WS, vn_WS_subsurfaceflow)
    end if


    if (processCase(9) == 1_i4) then
      if  ((vn_baseflow .eq. "none") .or. &
           (vn_WS_baseflow .eq. "none") .or. &
           (vn_fromWS_baseflow .eq. "none")) then
          call message()
          call message('***ERROR: baseflow activated, but baseflow or baseflowReservoir not provided')
          stop
      end if
      ! collect WS variable names
      call append(all_vn_WS, vn_WS_baseflow)
    end if

    !  Surface impermeable storage and direct runoff generation 
    if (processCase(11) == 1_i4) then
       if ((vn_directflow .eq. "none") .or. &
           (vn_WS_directflow .eq. "none") .or. &
           (Sealed_LanduseNrs < 1_i4)) then
          call message()
          call message('***ERROR: direct runoff process activated, but variables not provided')
          stop
       end if
       allocate(Sealed_LanduseIDs(Sealed_LanduseNrs))
       Sealed_LanduseIDs(:) = Sealed_Landuse_IDs(1:Sealed_LanduseNrs)

       ! collect WS variable names
       call append(all_vn_WS, vn_WS_directflow)
    end if

    !------------------------
    !  hydrological structure 
    !------------------------
    call position_nml('hydrology_structure', unamelist)
    read(unamelist, nml=hydrology_structure)

    allocate(vn_WS_subsurfaceflow_ver(max(1,nSubsurfaceflow),4))
    allocate(vn_WS_subsurfaceflow_lat(max(1,nSubsurfaceflow),2))
    allocate(vn_WS_soil(nSoilLayers))
    allocate(vn_SatSM(nSoilLayers))
    allocate(vn_WP(nSoilLayers))
    allocate(vn_WS_soil_ver(nSoilLayers,4))
    vn_WS_subsurfaceflow_ver = "none"
    vn_WS_subsurfaceflow_lat= "none"
    vn_WS_soil = "none"
    vn_WS_soil_ver= "none"
    vn_SatSM = "none"
    vn_WP = "none"

    ! !check if "up_in" is provided for subsurface storage
    ! if (any(vn_WS_subsurface_flow_ver(:,1) .eq. "none")) then
        ! call message()
        ! call message('***ERROR: Input of subsurface water reservoirs up_in should be provided')
        ! stop
    ! end if
    ! if (vn_WS_baseflow_ver(1) .eq. "none") then
        ! call message()
        ! call message('***ERROR: Input of baseflow water reservoir up_in should be provided')
        ! stop
    ! end if

    vn_WS_subsurfaceflow_ver(:,:) = vn_WS_subsurface_flow_ver(1:max(1,nSubsurfaceflow),1:4)
    vn_WS_subsurfaceflow_lat(:,:) = vn_WS_subsurface_flow_lat(1:max(1,nSubsurfaceflow),1:2)
    vn_WS_soil(:) = vn_WSsoil(1:nSoilLayers)
    vn_WS_soil_ver(:,:) = vn_WSsoil_ver(1:nSoilLayers,1:4)
    vn_SatSM(:) = vn_Satsoilmoist(1:nSoilLayers)
    vn_WP(:) = vn_Wiltpnt(1:nSoilLayers)
    !re-infiltration
    if (processCase(6) == 1_i4) then
       if ((vn_reinfiltration .eq. "none") .or. &
           all(vn_reinfilsoil .eq. "none") ) then
          call message()
          call message('***ERROR: reinfiltration activated, but fluxes not provided')
          stop
       end if
       !corresponding to each soil layers
       allocate(vn_reinfil_soil(nSoilLayers))
       vn_reinfil_soil = "none"
       vn_reinfil_soil(1:nSoilLayers) = vn_reinfilsoil(1:nSoilLayers)
    end if


    ! collect WS variable names
    call append(all_vn_WS, vn_WS_soil)

    !------------------------
    !  channel routing structure 
    !------------------------
    call position_nml('channel_routing', unamelist)
    read(unamelist, nml=channel_routing)
    if (processCase(10) == 1_i4) then
        vn_chan_store = vn_channelstorage
        vn_chan_in = vn_upstreaminflow
        vn_chan_out  =  vn_downstreamoutflow
        if (nSubsrfToChannel /= nSubsurfaceflow) then
          call message('WARMING: are you sure:')
          call message('the number of subsurface seepage is not equal to that of subsurface runoff?')
        end if
        allocate(vn_subsrftochan(max(1,nSubsrfToChannel)))
        vn_subsrftochan = "none"
        vn_subsrftochan(:) = vn_sub_srftochan(1:max(1,nSubsrfToChannel))
    end if



    call close_nml(unamelist)

  !===============================================================
  ! Settings for Optimization directly from the mHM implementations
  !===============================================================
    call open_nml(file_namelist_wqm, unamelist, quiet=.true.)
    ! namelist for Optimization settings
    call position_nml('optimization', unamelist)
    read(unamelist, nml=optimization)

    ! check and set default values
    if (nIterations .le. 0_i4) then
       call message('Number of iterations for Optimization (nIterations) must be greater than zero')
       stop
    end if
    if (dds_r .lt. 0.0_dp .or. dds_r .gt. 1.0_dp) then
       call message('dds_r must be between 0.0 and 1.0')
       stop
    end if
    if (sce_ngs .lt. 1_i4) then
       call message ('number of complexes in SCE (sce_ngs) must be at least 1')
       stop
    end if
    ! ! number of points in each complex: default = 2n+1
    ! if (sce_npg .lt. 0_i4) then
       ! n_true_pars = count(nint(global_parameters(:,4)) .eq. 1)
       ! sce_npg = 2 * n_true_pars   + 1_i4! 1111+11_i4+
       ! !***********************
	   ! ! 11 only temporally for sce_nitrate
       ! !***********************
    ! end if
    ! ! number of points in each sub-complex: default = n+1
    ! if (sce_nps .lt. 0_i4) then
       ! n_true_pars = count(nint(global_parameters(:,4)) .eq. 1)
       ! sce_nps = n_true_pars  + 1_i4!1111+11_i4+
    ! end if
    if (sce_npg .lt. sce_nps) then
       call message ('number of points per complex (sce_npg) must be greater or')
       call message ('equal number of points per sub-complex (sce_nps)')
       stop
    end if

    


  !=================================================
  !READ parameter file "wqm_parameter.nml"
  !=================================================
  inquire(file = file_param_wqm, exist = fexist)
    if (fexist) then
      call open_nml(file_param_wqm, unamelist_param, quiet =.true.)
    
      call position_nml('parameters', unamelist_param)
    read(unamelist_param, nml=parameters)
    else
      call message('')
      call message('***ERROR: parameter namelist not exist!')
      stop
    end if

    !check if value is within the parameter ranges
    do idx =1, nTotalParam
      do idgrp = 1, Nparameters(idx)%ngroups
        if ((Nparameters(idx)%pvalue(idgrp) .lt. Nparameters(idx)%lowbound(idgrp)) .or. &
          (Nparameters(idx)%pvalue(idgrp) .gt. Nparameters(idx)%upbound(idgrp))) then
            call message('')
            call message('***ERROR: parameter value out of bounds, please check:')
            call message('parameter:', trim(Nparameters(idx)%param_name))
            stop
        end if
      end do
    end do
    
    allocate(wqm_param(nTotalParam))
    wqm_param(:) = Nparameters(1:nTotalParam)

    call close_nml(unamelist_param)
  
  !======================================================
  !READ parameter namelist "wqm_outputs.nml"
  !======================================================  
  inquire(file = file_stateoutput_wqm, exist = fexist)
  if (fexist) then
    call open_nml(file_stateoutput_wqm, ustatOutput, quiet =.true.)  
    call position_nml ('statevarsoutput', ustatOutput)
    read(ustatOutput, nml = statevarsoutput)
  end if
  allocate(outputFlxState_wqm(nOUTstate_wqm) )  
  outputFlxState_wqm(1:nOUTstate_wqm) = output_FlxState_wqm(1:nOUTstate_wqm)
  
  if (any(outputFlxState_wqm)) then  
     call message( '' )
     call message( 'Following state variables of water quality model will be written:' )
     if (outputFlxState_wqm(1)) then
        call message( 'Soil moisture concentration in each layer     [mg/l]')
     end if  
     if (outputFlxState_wqm(2)) then
        call message( 'direct runoff concentration (L1_cdirectRunoff)     [mg/l]')
     end if 
     if (outputFlxState_wqm(3)) then
        call message( 'Surface runoff concentration (L1_csurfaceRunoff)     [mg/l]')
     end if 
     if (outputFlxState_wqm(4)) then
        call message( 'Subsurface runoff concentration for each component defined (L1_csubsurfaceRunoff)[mg/l]')
     end if   
     if (outputFlxState_wqm(5)) then
        call message( 'Baseflow concentration (L1_cbaseflow)     [mg/l]')
     end if 
     if (outputFlxState_wqm(6)) then
        call message( 'Soil N stock in each layer (L1_dissolvedIN) [g/m^2/setted timestep]')
     end if
     if (outputFlxState_wqm(7)) then
        call message( 'Total uptake amount in terrestrial phase (L1_soilUptakeN) [g/m^2/setted timestep]')
     end if
     if (outputFlxState_wqm(8)) then
        call message( 'Total denitrification amount (L1_soilDenitri) [g/m^2/setted timestep]')
     end if
     if (outputFlxState_wqm(9)) then
        call message( 'Total mineralisation amount (L1_soilMineralN) [g/m^2/setted timestep]')
     end if
     if (outputFlxState_wqm(10)) then
        call message( 'Total frtman applied IN amount (L1_soilINfrtmanapp) [g/m^2/setted timestep]')
     end if
     if (outputFlxState_wqm(11)) then
        call message( 'concentration in direct flow to channel  (L1_concdirTOchan) [mg/l]')
     end if
     if (outputFlxState_wqm(12)) then
        call message( 'concentration in Surface flow to channel  (L1_concsrfTOchan) [mg/l]')
     end if
     if (outputFlxState_wqm(13)) then
        call message( 'concentration in Subsurface flow to channel for each component defined (L1_concsubsrfTOchan)[mg/l]')
     end if
     if (outputFlxState_wqm(14)) then
        call message( 'concentration in baseflow to channel   (L1_concbasefTOchan)[mg/l]')
     end if
     if (outputFlxState_wqm(15)) then
        call message( 'concentration in total flow to channel  (L1_concOUT)[mg/l]')
     end if
     if (outputFlxState_wqm(16)) then
        call message( 'concentration of each stream reach (L1_concTIN) [mg/l]')
     end if
     if (outputFlxState_wqm(17)) then
        call message( 'Total terrestrial export (L1_terXpt) [mg/m^2/setted timestep]')
     end if
     if (outputFlxState_wqm(18)) then
        call message( 'Instream denitrification amount (L1_aquaticDenitri) [kg/setted timestep]')

     end if
     if (outputFlxState_wqm(19)) then
        call message( 'Instream assimilatory uptake amount (L1_aquaticAssimil) [kg/setted timestep]')
     end if

  end if
  call close_nml(ustatOutput)
 end subroutine wqm_readconfig

END MODULE mo_wqm_readconfig
