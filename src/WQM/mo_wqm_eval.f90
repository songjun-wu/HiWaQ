!> \file mo_wqm_eval.f90

!> \brief Runs

!> \details  

!> \authors Xiaoqiang Yang 
!> \date Sep 2021


MODULE mo_wqm_eval

  USE mo_kind, ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: wqm_eval

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !      NAME
  !          wqm_eval

  !>        \brief

  !>        \details 
  !     RETURN
  !         None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         None

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Xiaoqiang Yang
  !>        \date Sep 2021

  SUBROUTINE wqm_eval(runoff,nutrient)

    use mo_message,              only : message
    use mo_utils,                only : ge
    use mo_string_utils,         only : num2str  
    use mo_julian,               only : caldat, julday,dec2date,ndays
    use mo_wqm_constants,        only : nodata_dp,nsubstances
    use mo_wqm_global_variables, only: &
         simPer, nTstepDay, timeStep,warmingDays, &
         nLanduseTotal, nSoilTotal, nSubsurfaceflow, &
         L1_nCells, L1_nLinks,&
         !nitrate submodel parameters		 
         L1_rdegradN, L1_rmineralN, L1_rdissolN,           & ! INOUT NE1 nitrate submodel parameters
         L1_rdeniSoil,                                     & ! INOUT NE1 nitrate submodel parameters
         L1_rdeniAqtc,L1_ratuptkN, L1_rpprodN,            & ! INOUT NE1 nitrate submodel parameters
         L1_fLAI, L1_fsoil,                               & ! 
         nSoilLayers,L1_soildepth,                        & !
         L1_humusN, L1_fastN, L1_dissolvedIN,L1_dissolvedON, & ! INOUT NS Four different Nitrate pools in each soil layer
         L1_cbaseflow, L1_cpercolate,     & ! INOUT NS conc. in each soil layer,dim2=soillayers,dim3=substances 
         year_start, day_preyr,                           &
         L1_netPerm, L1_fromN,L1_toN, L1_rOrder,          &
         const_crain,const_csnowmelt,                     &
         L1_csnowmelt,L1_cthroughfall,                    &
         L1_snowmelt,L1_throughfall, L1_infiltrate,       & !
         L1_evapSoil,L1_transp,                           &
         L1_surfaceRunoff,L1_subsurfaceRunoff,  &
         L1_csubsurfaceRunoff,L1_baseflow, vn_WS_surfaceflow_lat,  &
         !direct runoff from sealed areas
         L1_fSealed, L1_WSsealed, L1_Esealed, L1_directRunoff, prevstep_WSsealed, &
         L1_cWSsealed,L1_cdirectRunoff, &
         prevstep_WSponding, L1_csurfaceRunoff,L1_cinfilsurface,   &
         prevstep_WSsoil, L1_cinfilSoil,L1_creturnflow, L1_creinfiltSoil,           &
         prevstep_WSsurface,prevstep_WSsubstorage,                 &
         prevstep_WSbase,prevstep_baseflow,prevstep_percol,prevstep_baseinflow,prevstep_baseoutflow,&
         prevstep_cpercol,prevstep_cbaseinflow,prevstep_cbaseout,     &
         stat_LATsurface, stat_LATsubstorage,stat_LATbase, &
         !agri-management &soil N transformations
         L1_frotation,nRotationTotal,L1_temp, L1_SrfmixRatio,&
         L1_Puptake, L1_soilINfrtmanapp, L1_wiltingpoint, L1_saturatedSM,&
         vn_soilTemperature,L1_soiltemp, stat_soiltemp, L1_WSsnowpack,L1_soilMineralN, &
         L1_soilUptakeN,L1_soilDenitri,   &
         !
         L1_reinfiltSrf,L1_reinfiltSoil, &
         L1_WSsoil,L1_cWSsoil,    &
         L1_WSponding,L1_WSsurface, L1_WSsubsurface,L1_WSbase,L1_pcbase,L1_cWSbase,  &
         !channel
         L1_chanmask,nSubsrfToChannel, L1_totalRunoffTochan, L1_WSchannel_downout,&
         nLink_yravg_q, L1_width, L1_length, L1_avgtemp, &
         basin_gauge,   &
         R_areaCell, L1_qOUT, L1_concOUT, L1_terXpt,& !terrestrial export from local grid
         L1_dirTOchan,L1_concdirTOchan,L1_srfTOchan,L1_basefTOchan,L1_subsrfTOchan, &
         L1_concsrfTOchan, L1_concbasefTOchan, L1_concsubsrfTOchan, &
         L1_qTIN,L1_concTIN,L1_loadTIN,&
         L1_WSchannel, nLink_criverbox,nLink_rivert_avg10,nLink_rivert_avg20,nLink_rivertemp,  &
         L1_aquaticDenitri, L1_aquaticAssimil,&
         nor_globalradi, L1_rzcoeff,L1_flight, &			   
    ! output
         nHydroGauges,nWQGauges, &
         Hydro_DIS,WQM_nutrient, &
         timeStep_model_outputs_wqm, outputFlxState_wqm

    use mo_wqm_parameterization,     only: wqm_parameterization
    use mo_water_quality,            only: wqm_calc_initialise,&
                                           sealed_storage_directrunoff, &
                                           vertical_phase_dynamic, &
                                           terrestial_verticalconc_calculation, &
                                           terrestial_lateralconc_calculation, &
                                           agri_management, &
                                           soil_nutrient_transformation, &
                                           soil_plant_uptake,  &
                                           soil_denitrification,  &
                                           terrestial_lateralconc_accumulation, &
                                           terTochan_calculation, &
                                           instream_nutrient_processes, &
                                           mix_conc

    use mo_wqm_initialise,           only: wqm_variables_default_init
    use mo_wqm_write,                only: wqm_write_output_fluxes
    use mo_opti_variables,           only: optimize
    !use mo_wqm_restart,              only: wqm_read_restart_states
 
    implicit none

    real(dp), dimension(:,:),   allocatable, optional, intent(out) :: runoff        ! dim1=time dim2=gauge

    real(dp), dimension(:,:,:), allocatable, optional, intent(out) :: nutrient    ! dim1=time dim2=gauge dim3=substances

    ! -------------------------------------
    ! local variables
    ! counters and indexes
    integer(i4)                               :: ic,isl, k, ii, tt, hh      ! Counters
    integer(i4)                               :: iNode, tNode
    !integer(i4)                               :: nLink
    integer(i4)                               :: nTimeSteps
    real(dp)                                  :: newTime
    real(dp)                                  :: sec_TS   !seconds per  timestep(that should be in hours)
    integer(i4)                               :: year_counter     ! for yearly output
    integer(i4)                               :: average_counter  ! for averaging output

    integer(i4)                               :: writeout_counter ! write out time step
    integer(i4)                               :: day, month, year, hour,no_day,no_year
    logical                                   :: isday    !indicate day or night    
    !state-variables
    real(dp), dimension(nSoilLayers)              :: tmp_soilmoist
    real(dp), dimension(nSoilLayers,nsubstances)  :: tmp_concsoil   !to temporally store soilmoist conc
    real(dp)                                  :: depth   !empirical equation based water depth for in-stream process
    real(dp)                                  :: benthic_area !benthic area based on length and width
    real(dp)                                  :: newbox   !

    ! LAI options
    integer(i4)                               :: day_counter
    integer(i4)                               :: month_counter

    
	! for water quality model
    integer(i4)                               :: cg
    !type(OutputDatasetwqm)                    :: ncn       !for state outputs
	

    !-------------------------------------------------------------------
    ! Initalize State variables either to the default value or
    ! from the restart_files.
    ! All variables to be initalized had been allocated to the required
    ! space before this point (see, mo_startup: initialise)
    !-------------------------------------------------------------------
    !if (.NOT. read_restart ) then
       ! as default values,
       !**********************************
	   ! For water quality model, initialise all the variables for continuous running(e.g., optimization)
       !**********************************      
        call wqm_variables_default_init()
    !else
       ! read from restart files
       ! call wqm_read_restart_states(dirRestartIn )
       ! call message('  Reading WQM restart file: ',trim(dirRestartIn) &
	   !	  // 'WQM_restart.nc')
    !end if

    !water quality parameterization
    call wqm_parameterization(L1_nCells,L1_nLinks, nLanduseTotal,nSoilTotal, L1_fLAI,L1_fsoil, &
               L1_rdegradN, L1_rmineralN, L1_rdissolN,L1_rdeniSoil, &
               L1_rdeniAqtc,L1_ratuptkN, L1_rpprodN)

    ! calculate NtimeSteps for this basin
    nTimeSteps = ( simPer%julEnd - simPer%julStart + 1 ) * nTstepDay

    ! -0.5 is due to the fact that dec2date routine
    !   changes the day at 12:00 in NOON
    newTime = real(simPer%julStart,dp)

    ! Loop over time
    average_counter  = 0
    writeout_counter = 0
    hour = -timeStep
    sec_TS = real(timeStep, dp) * 3600_dp 

    do tt = 1, nTimeSteps

        hour = mod(hour+timeStep, 24)

        call caldat(int(newTime), yy=year, mm=month, dd=day)

        ! initalise counters
        if ( tt .EQ. 1 ) then
            day_counter   = day
            month_counter = month
            year_counter  = year
            !initial values for water quality model
            !if (.not. read_restart) then
              call wqm_calc_initialise(sec_TS,nSoilLayers, L1_soildepth, nLanduseTotal, L1_fLAI, &
                 nSoilTotal, L1_fsoil,nSubsurfaceflow, L1_humusN, L1_fastN, L1_dissolvedIN, L1_dissolvedON, &
                 L1_cdirectRunoff, &
                 L1_cbaseflow, L1_cpercolate,nLink_yravg_q,nLink_criverbox,nLink_rivertemp,&
                 nLink_rivert_avg10,nLink_rivert_avg20,stat_soiltemp)
            !!restart not yet implemented...
            !else
            !   call wqm_read_restart_states(iBasin)
            !end if

           !special initial values needed
           call dec2date(newTime-0.5_dp, yy=year_start)

        end if

        !get date of global calendar from Julian day
        call dec2date(newTime-0.5_dp, yy=year, mm=month, dd=day, hh=hour)
        isday = ( hour .gt. 6 ) .AND. ( hour .le. 18 )
	
        !get the number of days in current year and number of years
        no_day = ndays(day,month,year) - ndays(1,1,year) + 1
        no_year = year - year_start + 1
        !initial value for the days of previous year
        if (no_year == 1)  day_preyr =365_i4
        !update at the end of year
        if ((day == 31) .and. (month == 12)) then
        day_preyr = no_day   
        end if


      !******************************
      !call N calcualtion subroutines
      !******************************
        !Note mixing in interception and soilpack is ignored for N dynamics
        !constant values are assigned as input from the configure file
        do ic = 1, nsubstances
            !L1_crainfall_direct(:,ic) = const_crain(ic)  ![mg/l] conc.
            L1_csnowmelt(:,ic) = const_csnowmelt(ic)  ![mg/l] conc.
            !throughfall is given as the same as precipitation
            L1_cthroughfall(:,ic) = const_crain(ic)  ![mg/l] conc.
        end do
        !loop over grid cells/HRUs*SUBs

		!linkage between grid cells/subcatchments and channel routing
        !mostly vertical exchange
        do k = 1, L1_nCells
            !sealed storage and direct runoff generation
            if (L1_fSealed(k) > 1.0E-4_dp) then
                call sealed_storage_directrunoff(L1_snowmelt(k,tt),L1_throughfall(k,tt), &
                  L1_csnowmelt(k,:),L1_cthroughfall(k,:),L1_WSsealed(k,tt), L1_Esealed(k,tt), &
                  prevstep_WSsealed(k),L1_directRunoff(k,tt), &
                  L1_cWSsealed(k,:),L1_cdirectRunoff(k,:))
            end if
            !vertical dynamics --surface ponding and soil layer storage
            call vertical_phase_dynamic(k,tt,nSoilLayers, &
              L1_snowmelt(k,tt),L1_throughfall(k,tt), L1_infiltrate(k,tt),&
              L1_csnowmelt(k,:),L1_cthroughfall(k,:), &
              prevstep_WSponding(k),L1_cinfilsurface(k,:),&
              L1_evapSoil(k,tt,:),L1_transp(k,tt,:), &
              prevstep_WSsoil(k,:),L1_cinfilSoil(k,:,:))

            !update conc of conceptual storages and runoff components
            call terrestial_verticalconc_calculation(k,tt,nSubsurfaceflow, &
                  prevstep_WSsurface(k), L1_surfaceRunoff(k,tt),L1_SrfmixRatio(k),L1_csurfaceRunoff(k,:), &
                  L1_WSsoil(1)%reservoir(k,tt),L1_soildepth(k,1),&
                  prevstep_WSsubstorage(k,:),L1_subsurfaceRunoff(k,tt,:),L1_csubsurfaceRunoff(k,:,:),&
                  prevstep_WSbase(k),L1_baseflow(k,tt),prevstep_baseflow(k),prevstep_percol(k),&
                  prevstep_cpercol(k,:), prevstep_cbaseout(k,:), &
                  L1_cbaseflow(k,:), L1_cpercolate(k,:))

            !agricultural management--nutrient input from fertilizer and manure application (crop data)
            !AND potential plant uptake (potential_uptake)
            do isl=1, nSoilLayers
                tmp_soilmoist(isl) = prevstep_WSsoil(k,isl)
                tmp_concsoil(isl,:) = L1_cWSsoil(isl)%reservoir(k, :)
            end do


          
            call agri_management(k,timeStep, no_day, no_year, day_preyr,nRotationTotal, &
              L1_frotation(k,:),tmp_soilmoist(:), tmp_concsoil(:,:), L1_temp(k, tt),&
              L1_fastN(k,:), L1_humusN(k,:),   &
              L1_Puptake(k,:), L1_soilINfrtmanapp(k))
            !checking flux conc     
            !update dissolved ON and IN pools		  
            L1_dissolvedIN(k,:) = tmp_soilmoist(:) * tmp_concsoil(:,1)
            L1_dissolvedON(k,:) = tmp_soilmoist(:) * tmp_concsoil(:,2)

	        !transformation between differnt pools AND DENITRIFICATION
            call soil_nutrient_transformation(k, timeStep, nSoilLayers, L1_soildepth(k,:), L1_temp(k, tt), &
                   vn_soilTemperature,stat_soiltemp(k),L1_soiltemp(k,tt),L1_WSsnowpack%reservoir(k, tt), L1_wiltingpoint(k,:), &
                   L1_saturatedSM(k,:), tmp_soilmoist(:), tmp_concsoil(:,:), L1_fastN(k,:), &
                   L1_humusN(k,:),L1_dissolvedIN(k,:), L1_dissolvedON(k,:), &
                   L1_rdegradN(k), L1_rmineralN(k), L1_rdissolN(k), L1_soilMineralN(k))
            !soil denitrification process
            call soil_denitrification(k,timeStep, nSoilLayers, L1_saturatedSM(k,:), L1_soiltemp(k,tt), &
                   tmp_soilmoist(:), tmp_concsoil(:,:), L1_dissolvedIN(k,:), &
                   L1_rdeniSoil(k), L1_soilDenitri(k))
   
            !plant uptake and denitrification in soil phase 
            call soil_plant_uptake(k, L1_wiltingpoint(k,:), tmp_soilmoist(:), tmp_concsoil(:,:), &
                   L1_Puptake(k,:), L1_soilUptakeN(k))
            L1_dissolvedIN(k,:) = tmp_soilmoist(:) * tmp_concsoil(:,1)

            !update soil moisture conc.	to global variable		  
            do isl=1, nSoilLayers
                L1_cWSsoil(isl)%reservoir(k, :) = tmp_concsoil(isl,:) 
            end do

        end do

        !initialize some locally used variables
        L1_qTIN = 0.0_dp
        L1_loadTIN = 0.0_dp
        newbox = 0.0_dp
        !cell spatial linkage and channel routing
        do k = 1, L1_nLinks
            !based on the spatial organizations
            ii = L1_netPerm(k)
            iNode = L1_fromN(ii) !iNode is exactly the same as grid ids
            tNode = L1_toN(ii)
            ! catchment spatial organization
            if (L1_nLinks == L1_nCells) then  !lateral exchange might be considered
                if (any(vn_WS_surfaceflow_lat .ne. "none" )) then
                !potentially lateral exchange & update conc of runoff and storage
                  call terrestial_lateralconc_calculation(iNode,tNode,tt,nSubsurfaceflow,nSoilLayers, &
                    prevstep_WSsurface(iNode), stat_LATsurface(tNode), L1_SrfmixRatio(k),L1_csurfaceRunoff(iNode,:), &
                    prevstep_WSsubstorage(iNode,:),stat_LATsubstorage(tNode,:),L1_csubsurfaceRunoff(iNode,:,:),&
                    prevstep_WSbase(iNode),stat_LATbase(tNode),L1_baseflow(iNode,tt),prevstep_baseflow(iNode),&
                    prevstep_percol(iNode),prevstep_baseinflow(iNode), prevstep_baseoutflow(iNode),& 
                    prevstep_cpercol(iNode,:), prevstep_cbaseinflow(iNode,:),prevstep_cbaseout(iNode,:), &
                    L1_cbaseflow(iNode,:), L1_cpercolate(iNode,:),prevstep_WSsoil(iNode,:),&
                    prevstep_WSponding(iNode),L1_reinfiltSrf(iNode,tt),L1_reinfiltSoil(iNode,tt,:), &
                    L1_creinfiltSoil(iNode,:,:),L1_creturnflow(iNode,:,:) )
                end if
                if (L1_chanmask(iNode)) then
                !channel seepage generation (and additional inflows, if provided) mixing and conc update
                   call terTochan_calculation(iNode,tt,sec_TS, nSubsrfToChannel,nsubstances,R_areaCell(iNode),&
                   L1_fSealed(iNode), L1_directRunoff(iNode,tt),L1_cdirectRunoff(iNode,:), &
                   L1_srfTOchan(iNode,tt),L1_basefTOchan(iNode,tt),L1_subsrfTOchan(iNode,tt,:), &
                   L1_concsrfTOchan(iNode,:), L1_concbasefTOchan(iNode,:), &
                   L1_concsubsrfTOchan(iNode,:,:), &
                   L1_totalRunoffTochan(iNode),L1_qOUT(iNode), L1_concOUT(iNode,:),L1_terXpt(iNode,:))
                   L1_avgtemp(iNode) = L1_temp(iNode,tt)
                   !just to update global variables for direct runoff seepage
                   L1_dirTOchan(iNode,tt) = L1_directRunoff(iNode,tt)
                   L1_concdirTOchan(iNode,:) = L1_cdirectRunoff(iNode,:)
                end if

            else  ! semi-distributed or rout_resolution > resolution ==> ter. lateral exchange is not allowed!
                !accumulate modeling resoultion (ncells) to terrestial&in-stream routing resolution (nlinks)
                !runoff components become channel seepage directly
                call terrestial_lateralconc_accumulation(iNode,tt,sec_TS,nsubstances,nSubsurfaceflow, &
                  L1_fSealed(:), L1_directRunoff(:,tt),L1_cdirectRunoff(:,:), & !direct runoff from sealed storage
                  L1_surfaceRunoff(:,tt),L1_subsurfaceRunoff(:,tt,:),L1_baseflow(:,tt), &  !runoff components
                  L1_csurfaceRunoff(:,:),L1_csubsurfaceRunoff(:,:,:),L1_cbaseflow(:,:), &  !runoff conc for each components
                  L1_dirTOchan(iNode,tt),L1_srfTOchan(iNode,tt),L1_basefTOchan(iNode,tt),L1_subsrfTOchan(iNode,tt,:), &
                  L1_concdirTOchan(iNode,:), L1_concsrfTOchan(iNode,:), L1_concbasefTOchan(iNode,:),&
                  L1_concsubsrfTOchan(iNode,:,:),  & !corresponding aggregated conc
                  L1_totalRunoffTochan(iNode),L1_qOUT(iNode), L1_concOUT(iNode,:),L1_terXpt(iNode,:), &
                  L1_temp(:,tt),L1_avgtemp(iNode))
            end if

         	!channel routing & instream processes
            if (L1_chanmask(iNode)) then
              if (iNode /= tNode) then
               ! average channel discharge (for empirial channel width and depth calculation) [m^3/s]
                nLink_yravg_q(iNode) = nLink_yravg_q(iNode) + (L1_WSchannel_downout(iNode,tt) - nLink_yravg_q(iNode)) / (90.0_dp )  
               ! empirial equations for channel width and depth
                if (all((L1_width(:) - nodata_dp) > -0.00001_dp)) then
                    L1_width(iNode) = 5.40_dp *nLink_yravg_q(iNode) **(0.50_dp)  ![m] from M. Rode (2016), EST
                end if
                depth = 0.27_dp *nLink_yravg_q(iNode) **(0.39_dp)  ! [m] from J.A. Moody(2002), Earth Surface Processes and Landforms
                benthic_area = L1_width(iNode) * L1_length(iNode) ![m^2]
               !riverstore
               !total input at "iNode" (which is the same as)
                L1_qTIN(iNode) = L1_qTIN(iNode) + L1_qOUT(iNode) ![m^3/s]
               !conc. of total infow of node(at iNode)
                L1_concTIN(iNode,:) = (L1_loadTIN(iNode,:) + &
                                       L1_qOUT(iNode) * L1_concOUT(iNode,:) * sec_TS) /&
                                      (L1_qTIN(iNode) * sec_TS)

                newbox = L1_WSchannel(iNode,tt) + L1_qTIN(iNode) * sec_TS  ![m^3]

               !update conc. and water volumn of reach iNode
                nLink_criverbox(iNode,:) = (nLink_criverbox(iNode,:) * L1_WSchannel(iNode,tt) + L1_loadTIN(iNode,:) + &
                   L1_qOUT(iNode) * L1_concOUT(iNode,:) * sec_TS) / newbox          

                !instream denitrification and primary production
                !water temperature is assumed as 20-day's average of air temperature
                nLink_rivertemp(iNode) = nLink_rivertemp(iNode) + (L1_avgtemp(iNode) - nLink_rivertemp(iNode)) / 20.0_dp

                !10- and 20-day moving mean temperature of river water
                nLink_rivert_avg10(iNode) = nLink_rivert_avg10(iNode) + &
                                            (nLink_rivertemp(iNode) - nLink_rivert_avg10(iNode)) / (10.0_dp )
                nLink_rivert_avg20(iNode) = nLink_rivert_avg20(iNode) + &
                                            (nLink_rivertemp(iNode) - nLink_rivert_avg20(iNode)) / (20.0_dp )
                !nLink_riverbox(iNode) = newbox

                call instream_nutrient_processes(timeStep,no_day,iNode, nLink_rivertemp(iNode), nLink_rivert_avg10(iNode),&
                   nLink_rivert_avg20(iNode),benthic_area, depth, newbox, nLink_criverbox(iNode,:), &
                   L1_aquaticDenitri(iNode), L1_aquaticAssimil(iNode),L1_rdeniAqtc(iNode),L1_ratuptkN(iNode),&
                   L1_rpprodN(iNode),nor_globalradi(tt), L1_rzcoeff(iNode),L1_flight(iNode))

                !update water volumn of reach(link) 
                !nLink_riverbox(iNode) = nLink_riverbox(iNode) - L1_WSchannel_downout(iNode,tt) * sec_TS
               !calculate the load from each upstream reach and store as input of the to_node (tNode) 
                L1_loadTIN(tNode,:) = L1_loadTIN(tNode,:) + &
                     nLink_criverbox(iNode,:) * L1_WSchannel_downout(iNode,tt) * sec_TS   ![mg/l*m3/s*s] == [g]
               !accumulate flow to the 'to node' of the current link, as the upstream inflow of 'to node'
                L1_qTIN(tNode) = L1_qTIN(tNode) + L1_WSchannel_downout(iNode,tt)

              else
                !*************************************
                ! at the outlet of catchment
                !*************************************
                ! "tNode" == "iNode"		
                L1_qTIN(tNode) = L1_qTIN(tNode) + L1_qOUT(tNode)
                L1_concTIN(tNode,:) = (L1_loadTIN(tNode,:) +  &
                       L1_qOUT(tNode)* L1_concOUT(tNode,:) * sec_TS) / &
                       (L1_qTIN(tNode) * sec_TS)
              end if
            end if              
        end do

    !store previous step fluxes/statues
        prevstep_WSponding(:) = L1_WSponding%reservoir(:,tt)
        do hh =1, nSoilLayers
            prevstep_WSsoil(:,hh) = L1_WSsoil(hh)%reservoir(:,tt)
        end do
        prevstep_WSsurface(:) = L1_WSsurface%reservoir(:,tt)
        do hh =1, nSubsurfaceflow
            prevstep_WSsubstorage(:,hh) = L1_WSsubsurface(hh)%reservoir(:,tt)
        end do
        prevstep_WSbase(:) = L1_WSbase%reservoir(:,tt)
        prevstep_baseflow(:) = L1_baseflow(:,tt)
        prevstep_percol(:) = L1_WSbase%up_in(:,tt) 
        prevstep_baseinflow(:) = L1_WSbase%lateral_in(:,tt)
        prevstep_baseoutflow(:) = L1_WSbase%lateral_out(:,tt)
        prevstep_cpercol(:,:) = L1_pcbase%reservoir(:,:)
        prevstep_cbaseinflow(:,:) = L1_cWSbase%lateral_in(:,:)
        prevstep_cbaseout(:,:) = L1_cbaseflow(:,:)
        stat_LATsurface = 0.0_dp
        stat_LATsubstorage = 0.0_dp
        stat_LATbase = 0.0_dp

    !storing the discharge at evaluation gauges
        do cg=1, nHydroGauges
            Hydro_DIS(tt, basin_gauge%gaugeid(cg)) = L1_WSchannel_downout(basin_gauge%gaugeid_loc(cg),tt)
                                                     !L1_qTIN(basin_gauge%gaugeid_loc(cg)) !
        end do
          !storing calcualted concentration at evaluation gauges
        do cg=1, nWQGauges
            !1for IN, 2 for ON
            WQM_nutrient(tt, basin_gauge%wqm_gaugeid(cg), 1:2) = L1_concTIN(basin_gauge%wqm_gaugeid_loc(cg),1:2)
            !3 for TN (=1+2)                
            WQM_nutrient(tt,basin_gauge%wqm_gaugeid(cg),3)=sum(WQM_nutrient(tt, basin_gauge%wqm_gaugeid(cg), 1:2))
        end do

          ! update the counters
        if (day_counter   .NE. day  ) day_counter   = day
        if (month_counter .NE. month) month_counter = month
        if (year_counter  .NE. year)  year_counter  = year

          ! increment of timestep
        newTime = julday(day,month,year) + real(hour+timeStep,dp)/24._dp

        if (.not. optimize) then
            if (any (outputFlxState_wqm) .and. (tt-warmingDays*nTstepDay) > 0_i4) then! then
                call wqm_write_output_fluxes(&
                    ! output specification
                    timeStep_model_outputs_wqm, &
                    ! time specification. newTime has been updated to next step
                    warmingDays, newTime, nTimeSteps, nTStepDay, &
                    tt, &
                    ! parse present date to water quality states writer
                    day_counter, month_counter, year_counter)
            end if
        end if

    end do !<< TIME STEPS LOOP
    ! =========================================================================
	! SET NUTRIENT OUTPUT VARIABLE
    ! =========================================================================
    if (present(runoff))   runoff = Hydro_DIS
    if (present(nutrient)) nutrient = WQM_nutrient


  end SUBROUTINE wqm_eval


END MODULE mo_wqm_eval
