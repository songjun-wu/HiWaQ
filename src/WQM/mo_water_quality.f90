!> \file mo_water_quality.f90

!> \brief All main processes of water quality model.

!> \authors Xiaoqiang Yang
!> \date Sep 2021

MODULE mo_water_quality


  USE mo_kind, ONLY: i4, sp, dp

  IMPLICIT NONE


  PUBLIC :: wqm_calc_initialise  !initialize for calculation if restart is not activated
  PUBLIC :: sealed_storage_directrunoff !direct runoff from impermeable storage (if activated)
  PUBLIC :: vertical_phase_dynamic    !water and N mixing in different storages in the terrestrial soil phase
  PUBLIC :: terrestial_verticalconc_calculation !update conc due to vertical water dynamics
  PUBLIC :: agri_management !consider agricultural management
  PUBLIC :: soil_nutrient_transformation !soil N biogeochemical transformations
  PUBLIC :: soil_plant_uptake !the actual plant/crop N uptake given soil water/N availability
  PUBLIC :: soil_denitrification !
  PUBLIC :: terrestial_lateralconc_calculation !calculate concentration in lateral runoff components
  PUBLIC :: terrestial_lateralconc_accumulation !accumulated terrestial values for drainage networking and stream routing, for subbasin and mhm-structure sum-ups
  PUBLIC :: terTochan_calculation !calculate concentration of flow generation where channel mask is true
  PUBLIC :: instream_nutrient_processes !instream N transformation processes
  PUBLIC :: mix_conc

CONTAINS

  ! ------------------------------------------------------------------
  !     NAME
  !         wqm_calc_initialise

  !     PURPOSE
  !>        \brief initialisation of water quality model. 

  !     LITERATURE
  !         HYPE model

  !     HISTORY
  !>        \authors Xiaoqiang Yang 
  !>        \date Sep 2021

  subroutine wqm_calc_initialise(sec_TS,nsoillayer, soil_depth, nlandusetype,fLAI,  &
      nsoiltypes, fsoil, nSubsrflow, humus_N, fast_N, diss_IN, diss_ON, cdirectRunoff,cbaseflow, cpercol, &
      yravg_q,criverbox,rivertemp,rivert_avg10,rivert_avg20,stat_soiltemp)

  use mo_wqm_global_variables,   only: L1_cWSsoil,L1_cWSsubsurface, L1_cWSbase, L1_cWSsealed, & !storagereservoir type
       L1_WSsoil,L1_WSbase, L1_WSsurface,L1_WSintercept,L1_WSsubsurface, L1_baseflow,&
       init_type, init_humusN, init_fastN, hnhalf, init_concIN, init_concON, &
       !read in from file: initial_value.txt
       !readin init_humusN&init_fastN unit[mg/m^3] 
       L1_pcbase,  const_crain, &
       prevstep_WSintercept,prevstep_WSsurface,prevstep_WSbase,prevstep_WSsubstorage,prevstep_WSsoil, &
       prevstep_baseflow, prevstep_percol,prevstep_baseinflow,prevstep_baseoutflow, &
       prevstep_cpercol,prevstep_cbaseinflow, prevstep_cbaseout, &
       L1_WSchannel_downout
  implicit none

  real(dp),                           intent(in)    :: sec_TS  ![s]
  integer(i4),                        intent(in)    :: nsoillayer
  real(dp), dimension(:,:),           intent(in)    :: soil_depth   ![m]
  integer(i4),                        intent(in)    :: nlandusetype
  real(dp), dimension(:,:),           intent(in)    :: fLAI
  integer(i4),                        intent(in)    :: nsoiltypes
  real(dp), dimension(:,:),           intent(in)    :: fsoil
  integer(i4),                        intent(in)    :: nSubsrflow
  real(dp), dimension(:,:),           intent(inout) :: humus_N         ![g/m^2]
  real(dp), dimension(:,:),           intent(inout) :: fast_N          ![g/m^2]
  real(dp), dimension(:,:),           intent(inout) :: diss_IN         ![g/m^2]
  real(dp), dimension(:,:),           intent(inout) :: diss_ON         ![g/m^2]
  real(dp), dimension(:,:),           intent(inout) :: cdirectRunoff   ![mg/l]
  real(dp), dimension(:,:),           intent(inout) :: cbaseflow       ![mg/l]
  real(dp), dimension(:,:),           intent(inout) :: cpercol         ![mg/l]
  real(dp), dimension(:),             intent(inout) :: yravg_q
  real(dp), dimension(:,:),           intent(inout) :: criverbox
  real(dp), dimension(:),             intent(inout) :: rivertemp
  real(dp), dimension(:),             intent(inout) :: rivert_avg10
  real(dp), dimension(:),             intent(inout) :: rivert_avg20
  real(dp), dimension(:),             intent(inout) :: stat_soiltemp
  !local
  integer(i4)     :: j,i
  real(dp), dimension(size(soil_depth,dim=1))  :: init_humusN0, init_fastN0 
  real(dp) , dimension(size(soil_depth,dim=1)) :: hnhalf0,half_N
  real(dp), dimension(size(soil_depth,dim=1),size(soil_depth,dim=2))  :: horizon_depth1   !temporal variable for soil layer depth
  !
  horizon_depth1(:,:) = soil_depth(:,:)    ![m]         

  init_humusN0(:) = 0.0_dp              ![g/m^2]
  init_fastN0(:) = 0.0_dp               ![g/m^2]
  hnhalf0(:) = 0.0_dp
  if (init_type == 1_i4) then
    do j=1, nlandusetype
      init_humusN0(:) = init_humusN0(:) + fLAI(:, j) * init_humusN(j) * horizon_depth1(:,1)
      init_fastN0(:)  = init_fastN0(:) + fLAI(:, j) * init_fastN(j) * horizon_depth1(:,1)
      hnhalf0(:) = hnhalf0(:) + fLAI(:,j) * hnhalf(j)
    end do
  else if (init_type == 2_i4) then
    do j=1, nsoiltypes
      init_humusN0(:) = init_humusN0(:) + fsoil(:, j) * init_humusN(j) * horizon_depth1(:,1)
      init_fastN0(:)  = init_fastN0(:) + fsoil(:, j) * init_fastN(j) * horizon_depth1(:,1)
      hnhalf0(:) = hnhalf0(:) + fsoil(:,j) * hnhalf(j)
    end do
  end if
  fast_N(:,1) = init_fastN0(:)
  humus_N(:,1) = init_humusN0(:)
  half_N(:) = log(2.0_dp)/ hnhalf0(:)
  if (nsoillayer > 1_i4) then  
    humus_N(:,2) = init_humusN0(:) * exp(-half_N(:) * horizon_depth1(:,2)/2.0_dp)
    fast_N(:,2) = init_fastN0(:)
    L1_cWSsoil(2)%reservoir(:,1) = init_concIN(:)
    L1_cWSsoil(2)%reservoir(:,2) = init_concON(:)
    diss_IN(:,2) = init_concIN(:) * L1_WSsoil(2)%reservoir(:,1) !here dim2 =1 means geting soil water at the first time step
    diss_ON(:,2) = init_concON(:) * L1_WSsoil(2)%reservoir(:,1)
  end if
  if (nsoillayer > 2_i4) then
    humus_N(:,3) = init_humusN0(:) * exp(-half_N(:) * (horizon_depth1(:,3) + horizon_depth1(:,2) - horizon_depth1(:,1))/2.0_dp)
    fast_N(:,3) = init_fastN0(:)
  end if
  !only the last soil layer conc. should be initially given, due to the plausibly longer residence time
  L1_cWSsoil(nsoillayer)%reservoir(:,1) = init_concIN(:)
  L1_cWSsoil(nsoillayer)%reservoir(:,2) = init_concON(:)
  !L1_cWSsoil(nsoillayer)%lateral_in(:,1) = init_concIN(:)
  !L1_cWSsoil(nsoillayer)%lateral_in(:,2) = init_concON(:)
  diss_IN(:,nsoillayer) = init_concIN(:) * L1_WSsoil(nsoillayer)%reservoir(:,1) !here dim2 =1 means geting soil water at the first time step
  diss_ON(:,nsoillayer) = init_concON(:) * L1_WSsoil(nsoillayer)%reservoir(:,1)
  do i =1, nSubsrflow
  L1_cWSsubsurface(i)%reservoir(:,1) = init_concIN(:)
  L1_cWSsubsurface(i)%reservoir(:,2) = init_concON(:)
  L1_cWSsubsurface(i)%lateral_in(:,1) = init_concIN(:)
  L1_cWSsubsurface(i)%lateral_in(:,2) = init_concON(:)
  end do

  L1_cWSsealed(:,1) = const_crain(1)
  L1_cWSsealed(:,2) = const_crain(2)
  cdirectRunoff(:,1) = const_crain(1)
  cdirectRunoff(:,2) = const_crain(2)
  
  cbaseflow(:,1) = init_concIN(:)
  cbaseflow(:,2) = init_concON(:)
  cpercol = cbaseflow
  !initialise "prevstep_" variables
  prevstep_WSintercept(:) = L1_WSintercept%reservoir(:,1)   !here dim2 =1 means geting soil water at the first time step          
  prevstep_WSsurface(:) = L1_WSsurface%reservoir(:,1)
  prevstep_WSbase(:) = L1_WSbase%reservoir(:,1) 
  do i =1, nSubsrflow
    prevstep_WSsubstorage(:,i) = L1_WSsubsurface(i)%reservoir(:,1)
  end do         
  do i =1, nsoillayer
    prevstep_WSsoil(:,i) = L1_WSsoil(i)%reservoir(:,1)
  end do    
  !reachtemp1 = 0.0_dp
  !initial setting for groundwater conc. calculating
  prevstep_baseflow(:) = L1_baseflow(:,1)
  prevstep_percol(:) = L1_WSbase%up_in(:,1)
  prevstep_baseinflow(:) = L1_WSbase%lateral_in(:,1)
  prevstep_baseoutflow(:) = L1_WSbase%lateral_out(:,1)
  prevstep_cpercol = L1_pcbase%reservoir
  prevstep_cbaseinflow = cbaseflow
  prevstep_cbaseout = cbaseflow
  !some other initial values
  yravg_q(:) = L1_WSchannel_downout(:,1)
  criverbox(:,1) = 1.0_dp  !mg/l
  criverbox(:,2) = 0.0_dp
  rivertemp(:) = 10.0_dp   !water temperature 10 degree celsius
  rivert_avg10 = 10.0_dp   !water temperature 10 degree celsius
  rivert_avg20 = 10.0_dp   !water temperature 10 degree celsius 
  stat_soiltemp= 10.0_dp   !initial soil temperature, if not externally supplied
  
  end subroutine wqm_calc_initialise

  ! ------------------------------------------------------------------

  !     NAME
  !         sealed_storage_directrunoff

  !     PURPOSE
  !>        \brief sealed_storage_directrunoff. 

  !>        \details   

  !     HISTORY
  !>        \author Xiaoqiang Yang 
  !>        \date May 2022
 
  subroutine sealed_storage_directrunoff(snowmelt,throughfall,&
         csnowmelt,cthroughfall,sealed_store, sealed_aet, &
         prev_WSsealed,directrunoff, csealedstore, cdirectRunoff)

    implicit none

    !integer(i4),                      intent(in) :: k
    !integer(i4),                      intent(in) :: tt
    real(dp),                         intent(in) :: snowmelt
    real(dp),                         intent(in) :: throughfall
    real(dp), dimension(:),           intent(in) :: csnowmelt
    real(dp), dimension(:),           intent(in) :: cthroughfall
    real(dp),                         intent(in)    :: sealed_store
    real(dp),                         intent(in)    :: sealed_aet
    real(dp),                         intent(inout)    :: prev_WSsealed
    real(dp),                         intent(in)    :: directrunoff
    real(dp),   dimension(:),         intent(inout)   :: csealedstore
    real(dp),   dimension(:),         intent(out)   :: cdirectRunoff
    !local
    real(dp)         :: sealtmp
    real(dp), dimension(size(csealedstore)) :: ltmp,ctmp

    cdirectRunoff(:) = 0.0_dp
    sealtmp = prev_WSsealed+snowmelt+throughfall
    !mixing
    if ( sealtmp > 1.0E-10_dp) then ![m]
       ltmp(:) = prev_WSsealed * csealedstore(:) +  &
                 snowmelt * csnowmelt(:) + throughfall * cthroughfall(:)
       ctmp(:) = ltmp(:) / sealtmp
    else
       ctmp(:) = 0.0_dp
    end if
    !update concentration in storage and output fluxes
    csealedstore(:) = ctmp(:)
    cdirectRunoff(:) = ctmp(:)
    prev_WSsealed = sealtmp - sealed_aet - directrunoff

  end subroutine sealed_storage_directrunoff

  ! ------------------------------------------------------------------

  !     NAME
  !         vertical_phase_dynamic

  !     PURPOSE
  !>        \brief vertical_phase_dynamic. 

  !>        \details   
  !         

  !     LITERATURE
  !         HYPE model
  !         INCA model

  !     HISTORY
  !>        \author Xiaoqiang Yang 
  !>        \date Sep 2021
  
  subroutine vertical_phase_dynamic(k,tt,nsoillayers,  &
               snowmelt,throughfall,srfinfil, &
               csnowmelt,cthroughfall, &
               prev_WSponding, csrfinfil, &
               evapsoil, transplant, &
               prev_WSsoil, cinfiltration)

    use mo_wqm_global_variables, only: &
          L1_cWSponding,L1_cWSsoil, L1_WSsoil,&
          const_crain


    implicit none
    integer(i4),                      intent(in) :: k
    !integer(i4),                      intent(in) :: tNode
    integer(i4),                      intent(in) :: tt
    integer(i4),                      intent(in) :: nsoillayers
    !real(dp),                         intent(in) :: rainfall_direct
    real(dp),                         intent(in) :: snowmelt
    real(dp),                         intent(in) :: throughfall
    real(dp),                         intent(in) :: srfinfil
    !real(dp), dimension(:),           intent(in) :: crainfall_direct
    real(dp), dimension(:),           intent(in) :: csnowmelt
    real(dp), dimension(:),           intent(in) :: cthroughfall
    real(dp),                         intent(inout) :: prev_WSponding
    real(dp), dimension(:),           intent(out):: csrfinfil
    real(dp), dimension(:),           intent(in) :: evapsoil !soil layers
    real(dp), dimension(:),           intent(in) :: transplant !soil layers
    real(dp), dimension(:),           intent(inout) :: prev_WSsoil
    !real(dp), dimension(:),           intent(inout) :: st_LATsoil
    real(dp), dimension(:,:),         intent(out) :: cinfiltration
    !real(dp), dimension(:,:),         intent(out) :: creturnflow

 
    !local
    real(dp)                                      ::tmp    !temporal var for water
    real(dp), dimension(size(cinfiltration))              ::ctmp,ltmp   !temporal var for conc./load

    integer(i4)                                   ::hh
   !all flux unit is [m]!!
   !********************************
   !above ground
   !********************************
    csrfinfil(:) = 0.0_dp
    ctmp(:) = 0.0_dp
    !csrfrunoff(:) = 0.0_dp
    tmp = snowmelt+throughfall
    ctmp(:) = (snowmelt * csnowmelt(:) + throughfall * cthroughfall(:))/tmp
    !mixing
    call mix_conc(prev_WSponding, L1_cWSponding%reservoir(k,:),  &
                  tmp,ctmp)
    !update concentration in storage and output fluxes
    L1_cWSponding%down_out(k,:) = L1_cWSponding%reservoir(k,:)
    prev_WSponding = prev_WSponding+tmp - srfinfil
    !update the up_in concentration of first soil layer
    L1_cWSsoil(1)%up_in(k,:) = L1_cWSponding%down_out(k,:)
    csrfinfil(:) = L1_cWSponding%down_out(k,:)


   !********************************
   !within soil layers
   !********************************	
    !every step initialise conc. in infiled water of lower two layers 
    cinfiltration(:,:) = 0.0_dp

    !---conc. in soil layers
    do hh =1, nsoillayers
       !considering the mixing of storage and the flux from the upper layer and lateral input
        call mix_conc(prev_WSsoil(hh),L1_cWSsoil(hh)%reservoir(k,:), &
                    L1_WSsoil(hh)%up_in(k,tt),L1_cWSsoil(hh)%up_in(k,:))

       !update conc in out-going fluxes
        L1_cWSsoil(hh)%down_out(k,:) = L1_cWSsoil(hh)%reservoir(k,:)
       !update soil storage
       prev_WSsoil(hh) = prev_WSsoil(hh) + L1_WSsoil(hh)%up_in(k,tt)
       !update up_in of next layer
        if (hh < nsoillayers) then 
            L1_cWSsoil(hh+1)%up_in(k,:) = L1_cWSsoil(hh)%down_out(k,:)
        end if
       !---conc. after soil ET		
       !considering enrichment effect due to evapotranspiration for each layer
        L1_cWSsoil(hh)%reservoir(k,:) = L1_cWSsoil(hh)%reservoir(k,:) *prev_WSsoil(hh)/ &
            (prev_WSsoil(hh)-evapsoil(hh)-transplant(hh))
        prev_WSsoil(hh) = prev_WSsoil(hh)-evapsoil(hh)-transplant(hh)
        prev_WSsoil(hh) = prev_WSsoil(hh)-L1_WSsoil(hh)%down_out(k,tt)
        !to avoid numerical error
        if (prev_WSsoil(hh) < 1.0E-10_dp) then
            L1_cWSsoil(hh)%reservoir(k,:) = 0.0_dp
        end if
        cinfiltration(hh,:) = L1_cWSsoil(hh)%reservoir(k,:)



    end do


  end subroutine vertical_phase_dynamic

  ! ------------------------------------------------------------------

  !     NAME
  !         terrestial_verticalconc_calculation

  !     PURPOSE
  !>        \brief terrestial_verticalconc_calculation. 

  !>        \details   

  !     HISTORY
  !>        \author Xiaoqiang Yang 
  !>        \date Sep 2021
 
  subroutine terrestial_verticalconc_calculation(k,tt,nSubsrflow,  &
         prev_WSsurface,srfrunoff, SrfmixSoilratio, csrfrunoff,  & 
         soilMoisttop, soildepthtop, &	 
         prev_WSsubstorage,subsrfRunoff, csubsrfRunoff, &
         prev_WSbase, baseflow, prev_baseflow,prev_percol, &
         prev_cpercol, prev_cbaseout, cbaseflow, cpercol)

    use mo_wqm_global_variables, only: &
          L1_WSsurface,L1_cWSsurface,L1_cWSponding, vn_fromWS_surfaceflow,L1_cWSsoil,&
          L1_WSsubsurface,L1_cWSsubsurface,L1_WSbase, &
          L1_pcsurface, L1_pcsubsurface,L1_pcbase, ratio_RetentionWS,&
          vn_WS_surfaceflow_ver, &
          vn_WS_subsurfaceflow_ver,processCase

    implicit none

    integer(i4),                      intent(in) :: k
    !integer(i4),                      intent(in) :: tNode
    integer(i4),                      intent(in) :: tt
    integer(i4),                      intent(in) :: nSubsrflow
    real(dp),                         intent(inout) :: prev_WSsurface
    real(dp),                         intent(in)    :: srfrunoff
    real(dp),                         intent(in)    :: SrfmixSoilratio
    real(dp),   dimension(:),         intent(out)   :: csrfrunoff
    real(dp),                         intent(in)    :: soilMoisttop
    real(dp),                         intent(in)    :: soildepthtop
    real(dp),   dimension(:),         intent(inout) :: prev_WSsubstorage
    real(dp),   dimension(:),         intent(in)    :: subsrfRunoff
    real(dp),   dimension(:,:),       intent(inout) :: csubsrfRunoff
    real(dp),                         intent(inout) :: prev_WSbase
    real(dp),                         intent(in)    :: baseflow
    real(dp),                         intent(inout)    :: prev_baseflow     !previous step baseflow related flow&conc
    real(dp),                         intent(inout)    :: prev_percol       !previous step baseflow related flow&conc
    real(dp), dimension(:),           intent(inout) :: prev_cpercol      !previous step baseflow related flow&conc
    real(dp), dimension(:),           intent(inout) :: prev_cbaseout     !previous step baseflow related flow&conc
    real(dp), dimension(:),           intent(inout) :: cbaseflow
    real(dp), dimension(:),           intent(inout) :: cpercol
 
    !local
    integer(i4)         :: ii
    !variables for the numerical solving method
    real(dp), dimension(size(cbaseflow))    :: k1, k2, k3, k4
    real(dp), dimension(size(cbaseflow))    :: y1, y2, y3
    real(dp)   :: RTWS,old_wout,new_wout  !retention storage averaged between current and previous step
    real(dp),dimension(size(cbaseflow))    :: old_influx,old_outflux,new_influx

    ! for each conceptual storage, 
    ! firstly to update the reservoir conc from the source water storage
    ! secondly to mix with lateral inputs and update conc in source water storage
    !---------------
    if (processCase(7) == 1_i4) then
    ! storage for surfaceflow generation
      !recharge the storage since there are "up_in" fluxes
      if (vn_WS_surfaceflow_ver(1) .ne. "none") then 
        !mix due to up_in fluxes
        call mix_conc(prev_WSsurface,L1_cWSsurface%reservoir(k,:), &
            L1_WSsurface%up_in(k,tt),L1_pcsurface%reservoir(k,:))
        !update storage due to up_in
        prev_WSsurface = prev_WSsurface + L1_WSsurface%up_in(k,tt)
      ! otherwise, replace the conc. of this conceptual storage directly by the source storage conc.
      else
        ! directly replace the conc with the source storage conc
        L1_cWSsurface%reservoir(k,:) = L1_pcsurface%reservoir(k,:)
      end if
      !specifically for the surface flow concentration, there might be substaintial mixing with top-layer soil water 
      !during the terrestrial surface transport/routing
      csrfrunoff(:) =  L1_cWSsurface%reservoir(k,:) * (1- SrfmixSoilratio) + &
                       L1_cWSsoil(1)%reservoir(k,:) * SrfmixSoilratio 
    end if
    !update storage due to surface flow generation
    !prev_WSsurface = prev_WSsurface - srfrunoff  
    !---------------    
    ! storage(s) for subsurface flow generation (can be multiple)
    do ii =1, nSubsrflow
        if (vn_WS_subsurfaceflow_ver(ii,1) .ne. "none") then  !the input flux of subsurface storage is specificed
            !mix due to up_in fluxes
            call mix_conc(prev_WSsubstorage(ii),L1_cWSsubsurface(ii)%reservoir(k,:), &
                    L1_WSsubsurface(ii)%up_in(k,tt),L1_pcsubsurface(ii)%Pnter_Rsvr%reservoir(k,:))           
            !update storage with up_in
            !prev_WSsubstorage(ii) = prev_WSsubstorage(ii) + L1_WSsubsurface(ii)%up_in(k,tt)
        else  !otherwise, directly updated by pointing to the source storage
            L1_cWSsubsurface(ii)%reservoir(k,:) = L1_pcsubsurface(ii)%Pnter_Rsvr%reservoir(k,:)
        end if
        !subsurface runoff
        csubsrfRunoff(ii,:) = L1_cWSsubsurface(ii)%reservoir(k,:)
        !update due to flow generation
        !prev_WSsubstorage(ii) = = prev_WSsubstorage(ii) - subsrfRunoff(ii)
    end do

    if (processCase(9) == 1_i4) then 
	!saturated zone --baseflow
	!equation introduced from INCA model, Wade et al., 2002 HESS
    !Also using the fourth-order Runge-Kutta technique to solve the differential equation
    old_influx = 0.0_dp
    old_outflux = 0.0_dp
    new_influx = 0.0_dp
    RTWS = 0.0_dp
    old_wout = 0.0_dp
    new_wout = 0.0_dp
    !
    old_influx(:) = prev_cpercol(:)*prev_percol  ![mg/l * m] == [g/m^2]
    old_wout = prev_baseflow
    old_outflux(:) = old_wout * prev_cbaseout(:)
    new_influx(:) = L1_pcbase%reservoir(k,:)*L1_WSbase%up_in(k,tt) 
    RTWS = (L1_WSbase%reservoir(k,tt) +prev_WSbase) * ratio_RetentionWS*0.5_dp
    new_wout = baseflow 
	!new lateral in prev_baseflow
    k1(:) = (old_influx(:) - old_outflux(:)) / RTWS
    y1(:) = prev_cbaseout(:) + k1(:) * 0.5_dp
    k2(:) = (0.5*(new_influx(:)+old_influx(:))-y1(:)*0.5*(old_wout+new_wout)) / RTWS
    y2(:) = prev_cbaseout(:) + k2(:) * 0.5_dp
    k3(:) = (0.5*(new_influx(:)+old_influx(:))-y2(:)*0.5*(old_wout+new_wout)) / RTWS
    y3(:) = prev_cbaseout(:) + k3(:)
    k4(:) = (new_influx(:) - y3(:)*new_wout) / RTWS
    cbaseflow(:) = prev_cbaseout(:) + (k1(:)+2*k2(:)+2*k3(:)+k4(:))/6.0_dp
    cpercol(:) = L1_pcbase%reservoir(k,:)
    end if


  end subroutine terrestial_verticalconc_calculation
  ! ------------------------------------------------------------------

  !     NAME
  !         agri_management

  !     PURPOSE
  !>        \brief agricultural management. 

  !>        \details calculates nutrient(Nitrogen at the moment) input from agricultural management (e.g. fertilizer 
  !>        and manure application) and from residuals by havesting or ploughing. Also, calcualtes potential nutrient 
  !>        uptake during the growing season (from plant date to havest date).
  !>        Notes: catch crops has not been implemented yet....   

  !     LITERATURE
  !         HYPE model

  !     HISTORY
  !>        \author Xiaoqiang Yang 
  !>        \date Jun 2016
  !>    Modified
  !>        X. Yang Jun 2017 enabled different time-step (hourly)
  !>        X. Yang Aug 2019 modified potential crop uptake calculation
  
  subroutine agri_management(k,TS, noday, noyear, days_prev , nCroptation, &
       frac_rotation, soilmoist, concsoil,temp, &
       fast_N, humus_N, potential_uptake,infrtmanapp)
  
  use mo_wqm_global_variables,   only: &
       cropdata, rotation              !&
	   
  implicit none
  integer(i4),               intent(in)    :: k
  integer(i4),               intent(in)    :: TS  !timestep [h]
  integer(i4),               intent(in)    :: noday, noyear,days_prev
  integer(i4),               intent(in)    :: nCroptation
  real(dp), dimension(:),    intent(in)    :: frac_rotation
  !real(dp), dimension(:),    intent(in)    :: fracRoots
  real(dp), dimension(:),    intent(in)    :: soilmoist        ![m]
  real(dp), dimension(:,:),  intent(inout) :: concsoil
  real(dp),                  intent(in)    :: temp                            !temperature
  real(dp), dimension(:),    intent(inout) :: fast_N, humus_N                 ![g/m^2]
  real(dp), dimension(:),    intent(inout) :: potential_uptake  
  real(dp),                  intent(inout) :: infrtmanapp                     ![g/m^2]


  !local variables
  real(dp),dimension(2,2)   :: frtman_nadd, res_nadd   ! amount of sources added to nitrogen pools [kg/km^2 or mg/m^2]
  integer(i4)       :: j, num_crp
  real(dp)          :: DT  !number of steps in one day
  integer(i4)       :: ir,idxr,jc, jc_prev, jc_next    !rotation index/crop index in a specific rotation
  integer(i4)       :: remidr, remidr_prev, remidr_next    !reminder after mod function to identify cropid in a rotation
  real(dp)          :: uptk_help, uptake_N  
  real(dp)          :: f_temp              ! temperature factor of soil uptake
  !
  
  DT = 24.0_dp / TS
  frtman_nadd = 0.0_dp 
  res_nadd    = 0.0_dp
  potential_uptake = 0.0_dp

  do ir = 1, nCroptation
     idxr= rotation%id(ir)
     if (frac_rotation(idxr) > 0.0_dp ) then 
        ! abbreviation of crop rotation variables
        num_crp = rotation%ncrops(ir)
        ! identify correct crop in specific year
        remidr = mod(noyear, num_crp)
        remidr_prev = mod(noyear-1, num_crp)
        remidr_next = mod(noyear+1, num_crp)  ! the year after, **modified 2019-07-17**
        if (remidr == 0_i4)   remidr = num_crp
        if (remidr_prev == 0_i4)   remidr_prev = num_crp
        if (remidr_next == 0_i4)   remidr_next = num_crp
        jc = rotation%crop(ir,remidr) 
        jc_prev = rotation%crop(ir,remidr_prev)
        jc_next = rotation%crop(ir,remidr_next)  !**modified 2019-07-17**
        !***get crop_id 'jc' from rotation_info.txt and then get corresponding crop management information from cropdata.txt***
        ! fertiliser application
        if ((noday >= cropdata(jc)%frtday1 .and. noday < cropdata(jc)%frtday1 + cropdata(jc)%frtperiod) .or. &
            (noday < (cropdata(jc_prev)%frtday1 + cropdata(jc_prev)%frtperiod - days_prev))) then
        frtman_nadd(1,1) = frtman_nadd(1,1) + frac_rotation(idxr) * cropdata(jc)%frtn1 * (1.0 - cropdata(jc)%frtdown1) &
                                            / (cropdata(jc)%frtperiod * DT)
        frtman_nadd(2,1) = frtman_nadd(2,1) + frac_rotation(idxr) * cropdata(jc)%frtn1 * cropdata(jc)%frtdown1  &
                                            / (cropdata(jc)%frtperiod * DT)
        end if
        !in case fert applied twice a year
        if (cropdata(jc)%frtday2 .ne. 0_i4) then
        if (noday >= cropdata(jc)%frtday2 .and. noday < (cropdata(jc)%frtday2 + cropdata(jc)%frtperiod)) then
        frtman_nadd(1,1) = frtman_nadd(1,1) + frac_rotation(idxr) * cropdata(jc)%frtn2 * (1.0 - cropdata(jc)%frtdown2) &
                                            / (cropdata(jc)%frtperiod * DT)
        frtman_nadd(2,1) = frtman_nadd(2,1) + frac_rotation(idxr) * cropdata(jc)%frtn2 * cropdata(jc)%frtdown2  &
                                            / (cropdata(jc)%frtperiod * DT)
        end if
        end if
        !in case fert applied twice in previous year and the second application period across the calendar year
        if (cropdata(jc_prev)%frtday2 .ne. 0_i4) then
        if (noday < (cropdata(jc_prev)%frtday2 + cropdata(jc_prev)%frtperiod - days_prev)) then
        frtman_nadd(1,1) = frtman_nadd(1,1) + frac_rotation(idxr) * cropdata(jc)%frtn2 * (1.0 - cropdata(jc)%frtdown2) &
                                            / (cropdata(jc)%frtperiod * DT)
        frtman_nadd(2,1) = frtman_nadd(2,1) + frac_rotation(idxr) * cropdata(jc)%frtn2 * cropdata(jc)%frtdown2  &
                                            / (cropdata(jc)%frtperiod * DT)
        end if
        end if
		!manure application
        if ((noday >= cropdata(jc)%manday1 .and. noday < cropdata(jc)%manday1 + cropdata(jc)%frtperiod) .or. &
            (noday < (cropdata(jc_prev)%manday1 + cropdata(jc_prev)%frtperiod - days_prev))) then
        frtman_nadd(1,1) = frtman_nadd(1,1) + frac_rotation(idxr) * cropdata(jc)%mann1 * (1.0 - cropdata(jc)%mandown1) &
                                            * cropdata(jc)%manfIN / (cropdata(jc)%frtperiod * DT)
        frtman_nadd(2,1) = frtman_nadd(2,1) + frac_rotation(idxr) * cropdata(jc)%mann1 * cropdata(jc)%mandown1  &
                                            * cropdata(jc)%manfIN / (cropdata(jc)%frtperiod * DT)
        frtman_nadd(1,2) = frtman_nadd(1,2) + frac_rotation(idxr) * cropdata(jc)%mann1 * (1.0 - cropdata(jc)%mandown1) &
                                            * (1.0 - cropdata(jc)%manfIN) / (cropdata(jc)%frtperiod * DT)
        frtman_nadd(2,2) = frtman_nadd(2,2) + frac_rotation(idxr) * cropdata(jc)%mann1 * cropdata(jc)%mandown1  &
                                            * (1.0 - cropdata(jc)%manfIN) / (cropdata(jc)%frtperiod * DT)
        end if
		!in case manure applied two times a year
        if (cropdata(jc)%manday2 .ne. 0_i4) then
        if (noday >= cropdata(jc)%manday2 .and. noday < (cropdata(jc)%manday2 + cropdata(jc)%frtperiod)) then
        frtman_nadd(1,1) = frtman_nadd(1,1) + frac_rotation(idxr) * cropdata(jc)%mann2 * (1.0 - cropdata(jc)%mandown2) &
                                            * cropdata(jc)%manfIN / (cropdata(jc)%frtperiod * DT)
        frtman_nadd(2,1) = frtman_nadd(2,1) + frac_rotation(idxr) * cropdata(jc)%mann2 * cropdata(jc)%mandown2  &
                                            * cropdata(jc)%manfIN / (cropdata(jc)%frtperiod * DT)
        frtman_nadd(1,2) = frtman_nadd(1,2) + frac_rotation(idxr) * cropdata(jc)%mann2 * (1.0 - cropdata(jc)%mandown2) &
                                            * (1.0 - cropdata(jc)%manfIN) / (cropdata(jc)%frtperiod * DT)
        frtman_nadd(2,2) = frtman_nadd(2,2) + frac_rotation(idxr) * cropdata(jc)%mann2 * cropdata(jc)%mandown2  &
                                            * (1.0 - cropdata(jc)%manfIN) / (cropdata(jc)%frtperiod * DT)
        end if
        end if
        !in case manure applied twice in previous year and the second application period across the calendar year	
        if (cropdata(jc_prev)%manday2 .ne. 0_i4) then
        if (noday < (cropdata(jc_prev)%manday2 + cropdata(jc_prev)%frtperiod - days_prev))  then
        frtman_nadd(1,1) = frtman_nadd(1,1) + frac_rotation(idxr) * cropdata(jc)%mann2 * (1.0 - cropdata(jc)%mandown2) &
                                            * cropdata(jc)%manfIN / (cropdata(jc)%frtperiod * DT)
        frtman_nadd(2,1) = frtman_nadd(2,1) + frac_rotation(idxr) * cropdata(jc)%mann2 * cropdata(jc)%mandown2  &
                                            * cropdata(jc)%manfIN / (cropdata(jc)%frtperiod * DT)
        frtman_nadd(1,2) = frtman_nadd(1,2) + frac_rotation(idxr) * cropdata(jc)%mann2 * (1.0 - cropdata(jc)%mandown2) &
                                            * (1.0 - cropdata(jc)%manfIN) / (cropdata(jc)%frtperiod * DT)
        frtman_nadd(2,2) = frtman_nadd(2,2) + frac_rotation(idxr) * cropdata(jc)%mann2 * cropdata(jc)%mandown2  &
                                            * (1.0 - cropdata(jc)%manfIN) / (cropdata(jc)%frtperiod * DT)
        end if
        end if
        !residuals
        if (cropdata(jc)%resday == 0_i4) cropdata(jc)%resperiod = 365_i4
        if ((noday >= cropdata(jc)%resday .and. noday < cropdata(jc)%resday + cropdata(jc)%resperiod) .or. &
            (noday < (cropdata(jc_prev)%resday + cropdata(jc_prev)%resperiod - days_prev))) then
        res_nadd(1,1) = res_nadd(1,1) + frac_rotation(idxr) * cropdata(jc)%resn * (1.0 - cropdata(jc)%resdown) &
                                      * cropdata(jc)%resfast / (cropdata(jc)%resperiod * DT)
        res_nadd(1,2) = res_nadd(1,2) + frac_rotation(idxr) * cropdata(jc)%resn * (1.0 - cropdata(jc)%resdown) &
                                      * (1.0 - cropdata(jc)%resfast) / (cropdata(jc)%resperiod * DT)
        res_nadd(2,1) = res_nadd(2,1) + frac_rotation(idxr) * cropdata(jc)%resn *  cropdata(jc)%resdown &
                                      * cropdata(jc)%resfast / (cropdata(jc)%resperiod * DT)
        res_nadd(2,2) = res_nadd(2,2) + frac_rotation(idxr) * cropdata(jc)%resn *  cropdata(jc)%resdown &
                                      * (1.0 - cropdata(jc)%resfast) / (cropdata(jc)%resperiod * DT)
        end if
		
		

		!*************	 
        !calculate Nitrate potential uptake because of plant/crop growth 
		!Based on the three-parameter logistic growth function
        !*************  
        uptake_N =0.0_dp
        if (temp < 5.0_dp) then
           f_temp = 0.0_dp
        else
           f_temp = min(1.0_dp, (temp-5.0)/20.0)
        end if
        !!normally plant in Spring and harvest in Autumn
        if (cropdata(jc)%plantd <= cropdata(jc)%harvestd) then    
           if (noday >= cropdata(jc)%plantd .and. noday < cropdata(jc)%harvestd) then
           uptk_help = (cropdata(jc)%up1 - cropdata(jc)%up2) * exp(-cropdata(jc)%up3*(noday - cropdata(jc)%plantd))
           uptake_N = cropdata(jc)%up1 *cropdata(jc)%up2*cropdata(jc)%up3 * uptk_help /(cropdata(jc)%up2+uptk_help)/ &
                     (cropdata(jc)%up2+ uptk_help)
           end if
        end if
        !Winter crop, plant in Autumn of previous year and harvest in the current year, in other words, plant date is latter than harvest date
        !in "cropdata.txt" we still use the day number of a year
        !in this case,the emergence date (emergd >0 ) can be specified, indicating the starting of growing period in spring
		!alternatively (emergd ==0), the parameters (e.g., up2) of logistic growth function should be ajusted
        if (cropdata(jc)%plantd > cropdata(jc)%harvestd) then
            if (cropdata(jc)%emergd > 1) then
               if (noday >= cropdata(jc)%emergd .and. noday < cropdata(jc)%harvestd ) then 
                  uptk_help = (cropdata(jc)%up1 - cropdata(jc)%up2) * exp(-cropdata(jc)%up3*(noday - cropdata(jc)%emergd)) 
                  uptake_N = cropdata(jc)%up1 *cropdata(jc)%up2*cropdata(jc)%up3 * uptk_help/ &
                            (cropdata(jc)%up2+uptk_help)/(cropdata(jc)%up2+ uptk_help)           
               elseif (noday < cropdata(jc)%emergd) then
                  uptk_help = (cropdata(jc)%up1 - cropdata(jc)%up2) * exp(-cropdata(jc)%up3*(noday + 365 - cropdata(jc)%plantd))
                  uptake_N = f_temp* cropdata(jc)%up1 *cropdata(jc)%up2*cropdata(jc)%up3 * uptk_help/&
                             (cropdata(jc)%up2+uptk_help)/(cropdata(jc)%up2+ uptk_help)
               else
                  uptake_N = 0.0_dp  
               end if
            else
                if (noday < cropdata(jc)%harvestd) then
                  uptk_help = (cropdata(jc)%up1 - cropdata(jc)%up2) * exp(-cropdata(jc)%up3*(noday - cropdata(jc)%plantd + 365))   
                  uptake_N = cropdata(jc)%up1 *cropdata(jc)%up2*cropdata(jc)%up3 * uptk_help/ &
                             (cropdata(jc)%up2+uptk_help)/(cropdata(jc)%up2+ uptk_help)
                else
                  uptake_N = 0.0_dp
                end if
            end if
        end if
        ! if the next year cropping winter crops, means the next year's crop will be planted in current years' winter
        if (cropdata(jc_next)%plantd > cropdata(jc_next)%harvestd) then  
           if (noday >= cropdata(jc_next)%plantd ) then
           uptk_help = (cropdata(jc_next)%up1 - cropdata(jc_next)%up2) &
                       * exp(-cropdata(jc_next)%up3*(noday - cropdata(jc_next)%plantd))
           uptake_N = f_temp*cropdata(jc_next)%up1 *cropdata(jc_next)%up2*cropdata(jc_next)%up3* &
                         uptk_help /(cropdata(jc)%up2+uptk_help)/(cropdata(jc)%up2+ uptk_help)
           end if
        end if
        potential_uptake(1) =potential_uptake(1)+frac_rotation(idxr)*uptake_N*  &
                             cropdata(jc)%uppsoil* (1- cropdata(jc)%deepsoil)
        potential_uptake(2) =potential_uptake(2) + frac_rotation(idxr) * uptake_N * &
                             (1-cropdata(jc)%uppsoil)*(1-cropdata(jc)%deepsoil)
        !yangx 2021-08 allowing N uptake occurs in all three layers, in accordance to root fractions
        ! added layer 3		
        potential_uptake(3) =potential_uptake(3) + frac_rotation(idxr) * uptake_N * cropdata(jc)%deepsoil!(1.0 - cropdata(jc)%uppsoil)
		
     end if  !end of (if frac_rotation>0)
  end do     !end of rotation loop

  

  !add to corresponding pools
  do j= 1, 2  !upper two layers
  !fertilizer and manure: orgN goes to fastN pool, inorgN goes to soil water(if > 0) or fastN pool 
  if (frtman_nadd(j,1) > 0.0_dp) then
     if (soilmoist(j) > 1.0E-10_dp) then
        call add_source_water(soilmoist(j), concsoil(j,1), frtman_nadd(j,1))
     else
        fast_N(j) = fast_N(j) + frtman_nadd(j,1)
     end if 
  end if
  if (frtman_nadd(j,2) > 0.0_dp) then
     fast_N(j) = fast_N(j) + frtman_nadd(j,2)
  end if
  !residuals: go to humusN and fastN pools
  if (sum(res_nadd(j,:)) > 0.0_dp ) then
     fast_N(j) = fast_N(j) + res_nadd(j,1)
     humus_N(j)= humus_N(j) + res_nadd(j,2)
  end if

  end do 
  !potential uptake
  do j= 1, 3  !all three layers --yangx 2021-08
  if (potential_uptake(j) > 0.0_dp) then
     potential_uptake(j) = potential_uptake(j) / DT    !g/m^2/d 



  end if  


 
  end do 
  infrtmanapp =  sum(frtman_nadd(:,1))
 
  end subroutine agri_management   

  ! ------------------------------------------------------------------

  !     NAME
  !         soil_nutrient_transformation

  !     PURPOSE
  !>        \brief Nutrient transformation between different pools in soil phase. 

  !>        \details calculates the transformation between four different pools in every soil layer, including 
  !>        degradation, minearalisation and dissolution.         

  !     LITERATURE
  !         HYPE model

  !     HISTORY
  !>        \author Xiaoqiang Yang 
  !>        \date Jun 2016
  !>    Modified
  !>        X. Yang Jun 2017 enabled different time-step (hourly)
  
  subroutine soil_nutrient_transformation(k,TS, nsoillayers, soil_depth, airtemp, vn_soiltemp,st_soiltemp, soiltemp, snowdp, &
         wp, sat, soilmoist, concsoil, fast_N, humus_N, diss_IN, diss_ON, &
         pardegradN, parmineraN, pardissolN, mineralN_soil)

    implicit none
    integer(i4),                   intent(in)    :: k
    integer(i4),                   intent(in)   :: TS
    integer(i4),                   intent(in)   :: nsoillayers     !
    real(dp), dimension(:),        intent(in)   :: soil_depth ![m]
    real(dp),                      intent(in)   :: airtemp       !
    character(256),                intent(in)   :: vn_soiltemp
    real(dp),                      intent(inout):: st_soiltemp
    real(dp),                      intent(inout):: soiltemp      !
    real(dp),                      intent(in)   :: snowdp        ![m]
    real(dp), dimension(:),        intent(in)   :: wp            !
    real(dp), dimension(:),        intent(in)   :: sat           !
    real(dp), dimension(:),        intent(in)   :: soilmoist     ![m]
    real(dp), dimension(:,:),      intent(inout):: concsoil      ![mg/l]
    real(dp), dimension(:),        intent(inout):: fast_N        ![g/m^2]
    real(dp), dimension(:),        intent(inout):: humus_N       ![g/m^2]
    real(dp), dimension(:),        intent(inout):: diss_IN       ![g/m^2]
    real(dp), dimension(:),        intent(inout):: diss_ON       ![g/m^2]
    real(dp),                      intent(in)   :: pardegradN    ![/d]
    real(dp),                      intent(in)   :: parmineraN    ![/d]
    real(dp),                      intent(in)   :: pardissolN    ![/d]
    real(dp),                      intent(out)  :: mineralN_soil ![g/m^2]
    !local variables
    real(dp)    :: DT
    real(dp)    :: fct_temp     ! temperature factor
    real(dp)    :: fct_sm       ! soil moisture factor
    real(dp)    :: degradN, mineraN, mineraN2,dissolN!transformated amount in each process
    integer(i4) :: j
    real(dp), dimension(nsoillayers)   :: mineralN_lyr
	

    DT = 24.0_dp / TS
    !if not externally provided by hydro inputs
    if (vn_soiltemp .eq. "none") then
    !calcualte soil temperature from HYPE functions
    !preserve previous step soil temperature
       call calculate_soil_temperature(airtemp, st_soiltemp, snowdp)
       soiltemp = st_soiltemp 
    end if
    !temperature factor
    fct_temp = tempfactor(soiltemp) 
    do j =1, nsoillayers
       degradN = 0.0_dp
       mineraN = 0.0_dp 
       mineraN2 = 0.0_dp
       dissolN = 0.0_dp

    !soil moisture factor
    fct_sm = moistfactor(soilmoist(j), wp(j), sat(j), soil_depth(j)) 
    !degradation: from humusN pool to fastN pool
       degradN =  pardegradN * fct_temp * fct_sm * humus_N(j) / DT 
       humus_N(j) = humus_N(j) - degradN
       fast_N(j) = fast_N(j) + degradN    
    !mineralisation: from fastN pool to dissolved IN
       mineraN =  parmineraN * fct_temp * fct_sm * fast_N(j) / DT
       fast_N(j)  = fast_N(j) - mineraN
       diss_IN(j) = diss_IN(j) + mineraN
	!update concentration of IN and ON
       concsoil(j,1) = diss_IN(j) / soilmoist(j)
       concsoil(j,2) = diss_ON(j) / soilmoist(j)
 
    !mineralisation: from dissolved ON to dissolved IN
       mineraN2 =  parmineraN * fct_temp * fct_sm * concsoil(j,2) / DT 
       concsoil(j,1) = concsoil(j,1) + mineraN2
       concsoil(j,2) = concsoil(j,2) - mineraN2

    !total mineralistaion amount
       mineralN_lyr(j) = mineraN + mineraN2*soilmoist(j) 
	!update dissolved IN and ON pools
       diss_IN(j) = concsoil(j,1) * soilmoist(j)
       diss_ON(j) = concsoil(j,2) * soilmoist(j)

    !dissolution: from fastN pool to dissolved ON
       ! fast_N(j) = fast_N(j) + diss_ON(j)
       ! newconc = fast_N(j) / (soilmoist(j) + pardissolN * soil_depth(j)/1000.0_dp)
       ! concsoil(j,2) = concsoil(j,2) + (newconc - concsoil(j,2))*(1.0_dp - exp(-0.01_dp / DT))  !1% perday
       dissolN = pardissolN * fct_temp * fct_sm * fast_N(j)
    !update fastN and dissolved ON pools
       diss_ON(j) = diss_ON(j) + dissolN
       fast_N(j) = fast_N(j) - dissolN
      
    end do
       mineralN_soil = sum(mineralN_lyr(:))




  end subroutine soil_nutrient_transformation

  ! ------------------------------------------------------------------

  !     NAME
  !         soil_plant_uptake

  !     PURPOSE
  !>        \brief plant uptake of Nitrate in soil phase. 

  !>        \details calculates nitrate uptake amount in soil water, considering potential uptake(crop dependent),
  !>         concentration, and soil moisture factors.               

  !     LITERATURE
  !         HYPE model

  !     HISTORY
  !>        \author Xiaoqiang Yang 
  !>        \date Jun 2016
  
  subroutine soil_plant_uptake(k,wp, soilmoist, concsoil, potential_uptake, soiluptkN)
  
  implicit none
    integer(i4),                   intent(in)    :: k
  real(dp), dimension(:),         intent(in)     :: wp
  real(dp), dimension(:),         intent(in)     :: soilmoist  ![m]
  real(dp), dimension(:,:),       intent(inout)  :: concsoil   ![mg/l]
  real(dp), dimension(:),       intent(in)     :: potential_uptake  ![g/m^2]
  real(dp),                       intent(out)    :: soiluptkN
  !local
  integer(i4)         :: j
  real(dp), dimension(3)    :: max_uptk, lyr_uptk
  real(dp), dimension(3)    :: diss_INt
  !hype originally only happens in upper two layers
  !yangx 2021-08 allowed three layers

  do j = 1, 3 
     if (soilmoist(j) > 1.0E-10_dp) then
     diss_INt(j) = soilmoist(j) * concsoil(j,1)
     !in case max_uptk less than zero
     if (soilmoist(j) > wp(j)) then
     max_uptk(j) = (soilmoist(j)- wp(j))/ soilmoist(j)
     else
     max_uptk(j) = 0.0_dp
     end if
     lyr_uptk(j) = min(potential_uptake(j), max_uptk(j)*diss_INt(j) )
    
	 !update dissolved inorganic pool and conc. in soil water
     diss_INt(j) = diss_INt(j) - lyr_uptk(j)
     concsoil(j,1) = diss_INt(j) / soilmoist(j)
     else
     lyr_uptk(j) =0.0_dp
     concsoil(j,1) = 0.0_dp
     diss_INt(j) = 0.0_dp
     end if

  end do  
	 
  soiluptkN = sum(lyr_uptk(:))


  end subroutine soil_plant_uptake

  ! ------------------------------------------------------------------

  !     NAME
  !         soil_denitrification

  !     PURPOSE
  !>        \brief Nitrate denitrification in soil water 

  !>        \details calculates nitrate denitrification in soil water, considering temperature, concentration, 
  !>        and soil moisture factors.         

  !     LITERATURE
  !         HYPE model

  !     HISTORY
  !>        \author Xiaoqiang Yang 
  !>        \date Jun 2016
  !>    Modified
  !>        X. Yang Jun 2017 enabled different time-step (hourly)

  subroutine soil_denitrification(k, TS, nsoillayers,sat, soiltemp, soilmoist, concsoil, &
                                  diss_IN, pardenirate, soildenitri)

  implicit none
    integer(i4),                   intent(in)    :: k
  integer(i4),                    intent(in)     :: TS
  integer(i4),                    intent(in)     :: nsoillayers
  real(dp), dimension(:),         intent(in)     :: sat
  real(dp),                       intent(in)     :: soiltemp
  real(dp), dimension(:),         intent(in)     :: soilmoist
  real(dp), dimension(:,:),       intent(inout)  :: concsoil
  real(dp), dimension(:),         intent(inout)  :: diss_IN
  real(dp),                       intent(in)     :: pardenirate
  real(dp),                       intent(out)    :: soildenitri
  !local
  real(dp), dimension(nsoillayers)  :: lyr_deni
  real(dp)    :: fct_temp, fct_sm, fct_conc  !terms for different factors
  integer(i4) :: j
  real(dp)    :: DT
  
  
  DT = 24.0_dp / TS
  ! factors for denitrification
  fct_temp = tempfactor(soiltemp)
  do j=1, nsoillayers
     if (soilmoist(j) > 1.0E-10_dp) then
     fct_sm = deni_moistfactor(soilmoist(j),sat(j) )
     fct_conc = concfactor(concsoil(j,1))
     !denitrification
     lyr_deni(j) = pardenirate * diss_IN(j) * fct_temp * fct_sm * fct_conc / DT
	 !update dissolved inorgainc pool
     diss_IN(j) = diss_IN(j) - lyr_deni(j)
     concsoil(j,1) = diss_IN(j) / soilmoist(j)
     else
     lyr_deni(j) = 0.0_dp
     concsoil(j,1) = 0.0_dp
     diss_IN(j) = 0.0_dp
     end if

  end do
  
  soildenitri = sum(lyr_deni(:) )
    
  end subroutine soil_denitrification


  ! ------------------------------------------------------------------

  !     NAME
  !         terrestial_lateralconc_calculation

  !     PURPOSE
  !>        \brief terrestial_lateralconc_calculation. 

  !>        \details   

  !     HISTORY
  !>        \author Xiaoqiang Yang 
  !>        \date Sep 2021
 
  subroutine terrestial_lateralconc_calculation(iNode,tNode,tt,nSubsrflow,nsoillayers,  &
         prev_WSsurface, st_LATsurface, SrfmixSoilratio, csrfrunoff,  &  
         prev_WSsubstorage,st_LATsubstorage,csubsrfRunoff, &
         prev_WSbase,st_LATbase, baseflow, prev_baseflow,prev_percol, prev_baseinflow,prev_baseoutflow, &
         prev_cpercol, prev_cbaseinflow, prev_cbaseout, cbaseflow, cpercol, &
         prev_WSsoil,prev_WSponding,reinfiltSrf,reinfiltSoil,creinfiltSoil,creturnflow)

    use mo_wqm_global_variables, only: processCase,  &
          L1_WSsurface, L1_cWSsurface, L1_pcsubsurface,&
          L1_WSsubsurface,L1_cWSsubsurface,L1_WSbase, L1_cWSbase, &
          L1_cWSsoil,L1_WSsoil,L1_cWSponding,L1_WSponding,&
          L1_pcbase, ratio_RetentionWS,&
          vn_WS_subsurfaceflow_lat,vn_WS_baseflow_lat, &
          vn_WS_soil_ver,vn_WS_ponding_ver,vn_WS_subsurfaceflow_ver

    implicit none

    integer(i4),                      intent(in) :: iNode
    integer(i4),                      intent(in) :: tNode
    integer(i4),                      intent(in) :: tt
    integer(i4),                      intent(in) :: nSubsrflow
    integer(i4),                      intent(in) :: nsoillayers
    real(dp),                         intent(inout) :: prev_WSsurface
    real(dp),                         intent(inout) :: st_LATsurface  !stat variable, summing up all lateral input of tNode grid cell
    real(dp),                         intent(in)    :: SrfmixSoilratio
    real(dp),   dimension(:),         intent(out)   :: csrfrunoff
    real(dp),   dimension(:),         intent(inout) :: prev_WSsubstorage
    real(dp),   dimension(:),         intent(inout) :: st_LATsubstorage
    real(dp),   dimension(:,:),       intent(inout) :: csubsrfRunoff
    real(dp),                         intent(inout) :: prev_WSbase
    real(dp),                         intent(inout) :: st_LATbase
    real(dp),                         intent(in)    :: baseflow
    real(dp),                         intent(inout)    :: prev_baseflow     !previous step baseflow related flow&conc
    real(dp),                         intent(inout)    :: prev_percol       !previous step baseflow related flow&conc
    real(dp),                         intent(inout)    :: prev_baseinflow   !previous step baseflow related flow&conc
    real(dp),                         intent(inout)    :: prev_baseoutflow  !previous step baseflow related flow&conc
    real(dp), dimension(:),           intent(inout) :: prev_cpercol      !previous step baseflow related flow&conc
    real(dp), dimension(:),           intent(inout) :: prev_cbaseinflow  !previous step baseflow related flow&conc
    real(dp), dimension(:),           intent(inout) :: prev_cbaseout     !previous step baseflow related flow&conc
    real(dp), dimension(:),           intent(inout) :: cbaseflow
    real(dp), dimension(:),           intent(inout) :: cpercol
    real(dp),   dimension(:),         intent(inout) :: prev_WSsoil
    real(dp),                         intent(inout) :: prev_WSponding
    real(dp),                         intent(in)    :: reinfiltSrf
    real(dp), dimension(:),           intent(in)    :: reinfiltSoil
    real(dp), dimension(:,:),         intent(inout) :: creinfiltSoil
	
    real(dp), dimension(:,:),         intent(inout) :: creturnflow
    !local
    integer(i4)         :: ii,hh
    !variables for the numerical solving method
    real(dp), dimension(size(cbaseflow))    :: k1, k2, k3, k4
    real(dp), dimension(size(cbaseflow))    :: y1, y2, y3
    real(dp)   :: RTWS,old_wout,new_wout  !retention storage averaged between current and previous step
    real(dp),dimension(size(cbaseflow))    :: old_influx,old_outflux,new_influx
    real(dp), dimension(size(cbaseflow))  :: creinfiltIN
    real(dp)                              :: reinfiltIN

    ! for each conceptual storage, 
    ! firstly to update the reservoir conc from the source water storage
    ! secondly to mix with lateral inputs and update conc in source water storage
    !---------------
    ! storage for surfaceflow generation
    !consider later fluxes   
        !mix due to lateral in
    call mix_conc(prev_WSsurface,L1_cWSsurface%reservoir(iNode,:), &
            L1_WSsurface%lateral_in(iNode,tt),L1_cWSsurface%lateral_in(iNode,:))
    !update lateral out
    L1_cWSsurface%lateral_out(iNode,:) = L1_cWSsurface%reservoir(iNode,:)
    !update the input concentration of "tNode"
    call mix_conc(st_LATsurface, L1_cWSsurface%lateral_in(tNode,:), &
        L1_WSsurface%lateral_out(iNode,tt), L1_cWSsurface%lateral_out(iNode,:))
    !
    st_LATsurface = st_LATsurface + L1_WSsurface%lateral_out(iNode,tt)
    !update final surface flow conc
    !csrfrunoff(:) =  L1_cWSsurface%reservoir(iNode,:)
    !specifically for the surface flow concentration, there might be substaintial mixing with top-layer soil water 
    !during the terrestrial surface transport/routing
    csrfrunoff(:) =  L1_cWSsurface%reservoir(iNode,:) * (1- SrfmixSoilratio) + &
                       L1_cWSsoil(1)%reservoir(iNode,:) * SrfmixSoilratio 
    !if reinfiltration is allowed update again
    !ponding == surface storage
    !update soil water concentration
    if ( processCase(6) ==1_i4 ) then
        !surface reinf--
        creinfiltSoil(:,:) = 0.0_dp
        creinfiltIN(:) = L1_cWSsurface%reservoir(iNode,:)
        reinfiltIN = reinfiltSrf
                !L1_(k,tt),L1_reinfiltSoil(k,tt,:), &
        
        !first layer
        do hh =1, nsoillayers
          !update soil water conc due to reinfiltration
          call mix_conc(prev_WSsoil(hh),L1_cWSsoil(hh)%reservoir(iNode,:), &
                    reinfiltIN,creinfiltIN(:))

          prev_WSsoil(hh) = prev_WSsoil(hh) + reinfiltIN - reinfiltSoil(hh)
		  
		  !update reinfilt for next layer
          creinfiltSoil(hh,:) =  L1_cWSsoil(hh)%reservoir(iNode,:)
          reinfiltIN = reinfiltSoil(hh)
          creinfiltIN(:) = creinfiltSoil(hh,:) 
          !to avoid numerical error
          if (prev_WSsoil(hh) < 1.0E-10_dp) then
            L1_cWSsoil(hh)%reservoir(iNode,:) = 0.0_dp
          end if
          creinfiltSoil(hh,:) = L1_cWSsoil(hh)%reservoir(iNode,:)

         ! if (iNode == 422) then
          ! write(23,*) "reinfil",hh, reinfiltIN*creinfiltIN(1),L1_cWSsoil(hh)%reservoir(iNode,1)

         ! end if
		  
        end do
       !!Replace previous "terrestial_verticalconc_calculation"::up_in due to the potential reinfiltration
       !update the conceptual storage pointer
        do ii =1, nSubsrflow
          if (vn_WS_subsurfaceflow_ver(ii,1) .ne. "none") then  !the input flux of subsurface storage is specified
            
            call mix_conc(prev_WSsubstorage(ii),L1_cWSsubsurface(ii)%reservoir(iNode,:), &
                    L1_WSsubsurface(ii)%up_in(iNode,tt),L1_pcsubsurface(ii)%Pnter_Rsvr%reservoir(iNode,:))           
          else  !otherwise, directly updated by pointing to the source storage
            L1_cWSsubsurface(ii)%reservoir(iNode,:) = L1_pcsubsurface(ii)%Pnter_Rsvr%reservoir(iNode,:)
          end if
        end do
    end if
    !---------------
    ! storage(s) for subsurface flow generation (can be multiple)
    do ii =1, nSubsrflow
       !considering the mixing of storage and the lateral flux
        if (any(vn_WS_subsurfaceflow_lat(ii,:) .ne. "none")) then
            !update storage with "up_in" first before "lateral_in"
            prev_WSsubstorage(ii) = prev_WSsubstorage(ii) + L1_WSsubsurface(ii)%up_in(iNode,tt) ! if up_in == "none" then L1_WSsubsurface(ii)%up_in(iNode,tt) = 0.0
            !mix due to lateral in
            call mix_conc(prev_WSsubstorage(ii),L1_cWSsubsurface(ii)%reservoir(iNode,:), &
                    L1_WSsubsurface(ii)%lateral_in(iNode,tt),L1_cWSsubsurface(ii)%lateral_in(iNode,:))
            !update the later output conc.
            L1_cWSsubsurface(ii)%lateral_out(iNode,:) = L1_cWSsubsurface(ii)%reservoir(iNode,:)
            !update the input concentration of "tNode"
            call mix_conc(st_LATsubstorage(ii), L1_cWSsubsurface(ii)%lateral_in(tNode,:), &
                L1_WSsubsurface(ii)%lateral_out(iNode,tt), L1_cWSsubsurface(ii)%lateral_out(iNode,:))
            st_LATsubstorage(ii) = st_LATsubstorage(ii) + L1_WSsubsurface(ii)%lateral_out(iNode,tt)
            !if up_in == none, indicating the source storage concentration should also be updated -- 2022-05
            if (vn_WS_subsurfaceflow_ver(ii,1) == "none") then
              L1_pcsubsurface(ii)%Pnter_Rsvr%reservoir(iNode,:) = L1_cWSsubsurface(ii)%reservoir(iNode,:)
            end if
        end if
        !subsurface runoff
        csubsrfRunoff(ii,:) = L1_cWSsubsurface(ii)%reservoir(iNode,:)
    end do


	!saturated zone --baseflow
	!equation introduced from INCA model, Wade et al., 2002 HESS
    !Also using the fourth-order Runge-Kutta technique to solve the differential equation
    old_influx = 0.0_dp
    old_outflux = 0.0_dp
    new_influx = 0.0_dp
    RTWS = 0.0_dp
    old_wout = 0.0_dp
    new_wout = 0.0_dp
 
    if (any(vn_WS_baseflow_lat(:) .ne. "none")) then
    !replace baseflow concentration calc with up_in percolation and lateral_in
        old_influx(:) = prev_cpercol(:)*prev_percol + prev_baseinflow*prev_cbaseinflow(:)
        old_wout = prev_baseflow + prev_baseoutflow
        old_outflux(:) = old_wout * prev_cbaseout(:)
        new_influx(:) = L1_pcbase%reservoir(iNode,:)*L1_WSbase%up_in(iNode,tt) + &
                    L1_cWSbase%lateral_in(iNode,:)*L1_WSbase%lateral_in(iNode,tt)
        RTWS = (L1_WSbase%reservoir(iNode,tt) +prev_WSbase)* ratio_RetentionWS *0.5_dp
        new_wout = baseflow + L1_WSbase%lateral_out(iNode,tt) 
	!new lateral in prev_baseflow
        k1(:) = (old_influx(:) - old_outflux(:)) / RTWS
        y1(:) = prev_cbaseout(:) + k1(:) * 0.5_dp
        k2(:) = (0.5*(new_influx(:)+old_influx(:))-y1(:)*0.5*(old_wout+new_wout)) / RTWS
        y2(:) = prev_cbaseout(:) + k2(:) * 0.5_dp
        k3(:) = (0.5*(new_influx(:)+old_influx(:))-y2(:)*0.5*(old_wout+new_wout)) / RTWS
        y3(:) = prev_cbaseout(:) + k3(:)
        k4(:) = (new_influx(:) - y3(:)*new_wout) / RTWS
        cbaseflow(:) = prev_cbaseout(:) + (k1(:)+2*k2(:)+2*k3(:)+k4(:))/6.0_dp
        cpercol(:) = L1_pcbase%reservoir(iNode,:)
    !update the input base-lateral flow concentration of "tNode"
        L1_cWSbase%lateral_out(iNode,:) = cbaseflow(:)
        call mix_conc(st_LATbase, L1_cWSbase%lateral_in(tNode,:), &
            L1_WSbase%lateral_out(iNode,tt), L1_cWSbase%lateral_out(iNode,:))
        st_LATbase = st_LATbase + L1_WSbase%lateral_out(iNode,tt)
    end if 
 

    !*******************************************************
    !mixing once more should be done if return flow occurs
    !*******************************************************
    ! the last soil layer "down_in" should be from storage reservoir of free water, or zero if not configured
    creturnflow(:,:) = 0.0_dp
    do hh = nsoillayers, 1, -1
        if (vn_WS_soil_ver(hh,4) .ne. "none") then !if "down_in" variable name of current soil layer is specified
            call mix_conc(prev_WSsoil(hh),L1_cWSsoil(hh)%reservoir(iNode,:), &
              L1_WSsoil(hh)%down_in(iNode,tt),L1_cWSsoil(hh)%down_in(iNode,:))
            !update conc in out-going return flow  
            L1_cWSsoil(hh)%up_out(iNode,:) = L1_cWSsoil(hh)%reservoir(iNode,:) 
            !prev_WSsoil(hh) = prev_WSsoil(hh) + L1_WSsoil(hh)%down_in(iNode,tt)
            if (hh > 1_i4) then   !update down_in of upper layer
                L1_cWSsoil(hh-1)%down_in(iNode,:) = L1_cWSsoil(hh)%up_out(iNode,:)
            else ! the first layer "up_out" is surface "down_in", just store it there for next judgement          
                L1_cWSponding%down_in(iNode,:) = L1_cWSsoil(hh)%up_out(iNode,:)
            end if
        end if
        creturnflow(hh,:) = L1_cWSsoil(hh)%reservoir(iNode,:)
    end do
    !surface ponding storage mixing once more should be done if return flow occurs
    if (vn_WS_ponding_ver(4) .ne. "none") then !avoid extremely low values
        call mix_conc(prev_WSponding, L1_cWSponding%reservoir(iNode,:),&
           L1_WSponding%down_in(iNode,tt),L1_cWSponding%down_in(iNode,:))       
    end if

  end subroutine terrestial_lateralconc_calculation




  ! ------------------------------------------------------------------

  !     NAME
  !         terTochan_calculation

  !     PURPOSE
  !>        \brief terTochan_calculation. 

  !>        \details   

  !     HISTORY
  !>        \author Xiaoqiang Yang 
  !>        \date Sep 2021
 
  subroutine terTochan_calculation(iNode,tt,sec_TS, nSubsrtochan,nsubstances, Cell_area,&
         fsealed, directTochan, concdirectTochan, & 
         srfTOchan,basefTOchan,subsrfTOchan, concsrfTOchan, concbasefTOchan,concsubsrfTOchan,&
         totalrunoff, qOUT, concOUT, loadexport)

    use mo_wqm_global_variables, only: &
          nInflowGauges, basin_gauge,L1_cWSsurface,L1_cWSbase,L1_cWSsubsurface

    implicit none

    integer(i4),                      intent(in) :: iNode
    integer(i4),                      intent(in) :: tt
    real(dp),                         intent(in) :: sec_TS        ![s]
    integer(i4),                      intent(in) :: nSubsrtochan
    integer(i4),                      intent(in) :: nsubstances
    real(dp),                         intent(in) :: Cell_area     ![m2]
    real(dp),                         intent(in) :: fsealed ![areal proportion]
    real(dp),                         intent(in) :: directTochan     ![m]
    real(dp),   dimension(:),         intent(inout) :: concdirectTochan   ![mg/l]
    real(dp),                         intent(in) :: srfTOchan     ![m]
    real(dp),                         intent(in) :: basefTOchan
    real(dp),   dimension(:),         intent(in) :: subsrfTOchan
    real(dp),   dimension(:),         intent(inout) :: concsrfTOchan   ![mg/l]
    real(dp),   dimension(:),         intent(inout) :: concbasefTOchan
    real(dp),   dimension(:,:),       intent(inout) :: concsubsrfTOchan
    real(dp),                         intent(inout) :: totalrunoff  !total runoff in[m/timestep]
    real(dp),                         intent(inout) :: qOUT     !terrestrial total output in [m^3/s]
    real(dp),   dimension(:),         intent(inout) :: concOUT  !concentration of total output runoff [mg/l]
    real(dp),   dimension(:),         intent(inout) :: loadexport !terrestrial export load [g/m^2/timestep]
    !real(dp),                         intent(in) :: prev_WSsurface
    !local
    integer(i4)         :: ii, nn
    real(dp),   dimension(nsubstances) :: tmp_locsubsrfinflux,tmp_loctotal_influx
    real(dp)                           :: tmp_locsubsrf

    tmp_locsubsrfinflux = 0.0_dp
    tmp_loctotal_influx = 0.0_dp
    tmp_locsubsrf = 0.0_dp
    concOUT =0.0_dp
    qOUT = 0.0_dp
    totalrunoff = 0.0_dp
    loadexport = 0.0_dp
    !runoff conponents and channel seepage match each other by default
    !channel surface flow generation
    concsrfTOchan(:) = L1_cWSsurface%reservoir(iNode,:) !L1_pcsrfTOchan%reservoir(iNode,:)
    !channel baseflow generation
    concbasefTOchan(:) = L1_cWSbase%lateral_out(iNode,:) !L1_pcbasefTOchan%reservoir(iNode,:)
    ! for channel each subsurface flow generation 
    do ii =1, nSubsrtochan
        concsubsrfTOchan(ii,:) = L1_cWSsubsurface(ii)%reservoir(iNode,:) !L1_pcsubsrfTOchan(ii)%Pnter_Rsvr%reservoir(iNode,:)
        tmp_locsubsrfinflux(:) = tmp_locsubsrfinflux(:) + subsrfTOchan(ii) * concsubsrfTOchan(ii,:)
        tmp_locsubsrf = tmp_locsubsrf + subsrfTOchan(ii) 
    end do
 
    !only channel connected cells consider the additional 
    !the total input flux from local terrestrial export  [mg/l * m] == [g/m^2]
    tmp_loctotal_influx(:) = (1-fsealed)*(tmp_locsubsrfinflux(:) + &   !subsurface
                             srfTOchan * concsrfTOchan(:) + & !surface
                             basefTOchan * concbasefTOchan(:)) + &  !baseflow
                             fsealed * directTochan * concdirectTochan
    totalrunoff = (1-fsealed)*(tmp_locsubsrf+srfTOchan+basefTOchan) + &
                   fsealed*directTochan            !first in [m]
    concOUT(:) =  tmp_loctotal_influx(:)/ totalrunoff  ![mg/l]
    !print*,totalrunoff,Cell_area ,sec_TS
    qOUT = totalrunoff * Cell_area / sec_TS  ![m] -> [m^3/s]

    loadexport(:) = tmp_loctotal_influx(:)  !total export [mg/l * m, the same as other budget outputs]


    !if there are addtional inflow (e.g., point source inputs)
    if (nInflowGauges .gt. 0_i4 ) then 
      !inflow gauge ids
      do nn =1, nInflowGauges 
        if (basin_gauge%inflowgaugeid_loc(nn) == iNode) then
            if (basin_gauge%inflow_ishead(nn)) then
              concOUT(:) = basin_gauge%InflowGaugeConc(tt,nn, 1:2) 
            else
              concOUT(:) = (concOUT(:)* qOUT+ &   !terrestrial
                basin_gauge%InflowGaugeConc(tt,nn, 1:2) *basin_gauge%inflowQ(tt,nn) )/ &  !additional input
                (qOUT+basin_gauge%inflowQ(tt,nn))   
            end if
        end if
      end do 
    end if


  end subroutine terTochan_calculation

  ! ------------------------------------------------------------------

  !     NAME
  !         terrestial_lateralconc_accumulation

  !     PURPOSE
  !>        \brief terrestial_lateralconc_accumulation. 

  !>        \details   

  !     HISTORY
  !>        \author Xiaoqiang Yang 
  !>        \date Sep 2021
 
  subroutine terrestial_lateralconc_accumulation(iNode,tt,sec_TS,nsubstances, nSubsrflow,&
      fsealed, directrunoff, cdirectrunoff, &  
      srfrunoff, subsrfRunoff,baseflow, csrfrunoff,csubsrfRunoff,cbaseflow,&
      dirTOchan, srfTOchan,basefTOchan,subsrfTOchan,&
      concdirTochan,concsrfTOchan, concbasefTOchan,concsubsrfTOchan, &
      totalrunoff, qOUT, concOUT, loadexport, airtemp, avgairtemp)

  use mo_wqm_global_variables,    only: model_type,RIdonM,L1_area, R_areaCell, L1_nHRUs,&
        nInflowGauges, basin_gauge,vn_srftochan,vn_baseftochan,vn_subsrftochan,const_crain
	  
  implicit none
  integer(i4),                    intent(in)     :: iNode       ! number of cells at L11
  integer(i4),                    intent(in)     :: tt           ! [h]time account
  real(dp),                       intent(in)     :: sec_TS           ! [s] time step in second
  integer(i4),                    intent(in)     :: nsubstances  ! muber of substances involved
  integer(i4),                    intent(in)     :: nSubsrflow   ! number of subsurface flow components
  !runoff and conc at nCells level
  real(dp), dimension(:),         intent(in)     :: fsealed      ![areal proportion]
  real(dp), dimension(:),         intent(in)     :: directrunoff     ![m]
  real(dp), dimension(:,:),       intent(inout)  :: cdirectrunoff  ![mg/l]
  real(dp), dimension(:),         intent(in)     :: srfrunoff    ! surface runoff dim = nCells
  real(dp), dimension(:,:),       intent(in)     :: subsrfRunoff ! subsurface runoff dim1 = nCells, dim2=nSubsrflow
  real(dp), dimension(:),         intent(in)     :: baseflow     ! baseflow dim = nCells
  real(dp), dimension(:,:),       intent(in)     :: csrfrunoff   ! conc surface runoff dim1=nCells, dim2=nsubstances
  real(dp), dimension(:,:,:),     intent(in)     :: csubsrfRunoff! conc subsurface runoff dim1 = nCells, dim2=nSubsrflow, dim3= nsubstances
  real(dp), dimension(:,:),       intent(in)     :: cbaseflow    ! conc baseflow  dim1=nCells, dim2=nsubstances
  !acc runoff at nLinks level
  real(dp),                       intent(inout)     :: dirTOchan       !direct runoff seepage
  real(dp),                       intent(inout)     :: srfTOchan       !surface seepage
  real(dp),                       intent(inout)     :: basefTOchan     !baseflow seepage
  real(dp), dimension(:),         intent(inout)     :: subsrfTOchan    !subsurface seepage  dim=nSubsrflow
  real(dp), dimension(:),         intent(inout)     :: concdirTochan   !conc surfaceseepage dim = nsubstances
  real(dp), dimension(:),         intent(inout)     :: concsrfTOchan   !conc surfaceseepage dim = nsubstances
  real(dp), dimension(:),         intent(inout)     :: concbasefTOchan !conc baseflowseepage dim = nsubstances
  real(dp), dimension(:,:),       intent(inout)     :: concsubsrfTOchan !conc surfaceseepage dim1 = nSubsrflow , dim2 = nsubstances
  
  real(dp),                         intent(out) :: totalrunoff ! total runoff in m
  real(dp),                         intent(out) :: qOUT     !terrestrial total output in [m^3/s]
  real(dp),   dimension(:),         intent(out) :: concOUT  !concentration of total output runoff [mg/l]
  real(dp),   dimension(:),         intent(out) :: loadexport !terrestrial export load [g/m^2/timestep]
  !temperature average
  real(dp),  dimension(:),         intent(in)     :: airtemp
  real(dp),                        intent(out)    :: avgairtemp


  !local
  !real(dp)                           :: tmp_dirTochan
  real(dp),   dimension(nsubstances) :: tmp_locsubsrfinflux,tmp_loctotal_influx!, tmp_concdirTochan
  integer(i4)                     :: ii, mm, nn,ij
  integer(i4)                     :: sloc, eloc
  real(dp)                        :: fsealed_rout

  avgairtemp = 0.0_dp
  concOUT =0.0_dp
  qOUT = 0.0_dp
  totalrunoff = 0.0_dp
  loadexport = 0.0_dp
  tmp_locsubsrfinflux = 0.0_dp
  tmp_loctotal_influx = 0.0_dp
  !tmp_dirTochan = 0.0_dp 
  !tmp_concdirTochan = 0.0_dp
  fsealed_rout = 0.0_dp


  !
    if (model_type ==1) then
      !sealed area at routing level & accumulate direct runoff
      fsealed_rout = sum(fsealed(:) * L1_area(:) , RIdonM(:) .eq. iNode )/ R_areaCell(iNode) 
      if (fsealed_rout > 1.0E-4_dp) then
        dirTOchan = sum(directrunoff(:) * L1_area(:) , RIdonM(:) .eq. iNode )/ R_areaCell(iNode) 
      end if
      !runoff seepage if not specified from hydroinput
      if (vn_srftochan .eq. "none") then
        srfTOchan = sum(srfrunoff(:) * L1_area(:) , RIdonM(:) .eq. iNode )/ R_areaCell(iNode) 
      end if
	
      if (vn_baseftochan .eq. "none") then
        basefTOchan = sum(baseflow(:) * L1_area(:), RIdonM(:) .eq. iNode ) / R_areaCell(iNode)
      end if

      if (any(vn_subsrftochan(:) .eq. "none")) then
        do nn = 1, nSubsrflow
          subsrfTOchan(nn) = sum(subsrfRunoff(:,nn) * L1_area(:), RIdonM(:) .eq. iNode)/ R_areaCell(iNode)
        end do
      end if
      !conc of the runoff seepage
      do mm = 1, nsubstances
        if (fsealed_rout > 1.0E-4_dp) then
          if (dirTOchan > 0.0_dp) then
            concdirTochan(mm) = sum(directrunoff(:)* L1_area(:)*cdirectrunoff(:,mm), RIdonM(:) .eq. iNode )/  &
                                      (dirTOchan*R_areaCell(iNode))
          else
            concdirTochan(mm) = const_crain(mm)
          end if
        end if        
        if (srfTOchan > 0.0_dp) then
          concsrfTOchan(mm) =  sum(srfrunoff(:) * L1_area(:)*csrfrunoff(:,mm) , RIdonM(:) .eq. iNode )/ &
                                     (srfTOchan*R_areaCell(iNode)) 
        else
          concsrfTOchan(mm) = const_crain(mm)
        end if
        if (basefTOchan > 0.0_dp) then
          concbasefTOchan(mm) =  sum(baseflow(:) * L1_area(:)*cbaseflow(:,mm) , RIdonM(:) .eq. iNode )/ &
                                     (basefTOchan*R_areaCell(iNode))
        else
          concbasefTOchan(mm) = const_crain(mm)
        end if
        do nn = 1,nSubsrflow
          if (subsrfTOchan(nn) > 0.0_dp) then
            concsubsrfTOchan(nn,mm) = sum(subsrfRunoff(:,nn) * L1_area(:)*csubsrfRunoff(:,nn,mm),RIdonM(:) .eq. iNode )/ &
                            (subsrfTOchan(nn)*R_areaCell(iNode)) 
          else
            concsubsrfTOchan(nn,mm) = const_crain(mm)
          end if
          tmp_locsubsrfinflux(mm) = tmp_locsubsrfinflux(mm) + concsubsrfTOchan(nn,mm) * subsrfTOchan(nn)
        end do
      end do
      !averaged air temperature
      avgairtemp = sum(airtemp(:) * L1_area(:) / R_areaCell(iNode) , RIdonM(:) .eq. iNode )
	!model_type ==2
    else !semi-distributed structure -- nLinks == nSubs
      sloc = L1_nHRUs*(iNode-1)+1
      eloc = L1_nHRUs*iNode
      !sealed area at routing level & accumulate direct runoff
      fsealed_rout = sum(fsealed(sloc:eloc) * L1_area(sloc:eloc)) / R_areaCell(iNode)
      if (fsealed_rout > 1.0E-4_dp) then
        dirTOchan = sum(directrunoff(sloc:eloc) * L1_area(sloc:eloc)) / R_areaCell(iNode)
      end if
      !runoff seepage if not specified from hydroinput
      if (vn_srftochan .eq. "none") then
        srfTOchan = sum(srfrunoff(sloc:eloc) * L1_area(sloc:eloc)) / R_areaCell(iNode)
      end if

      if (vn_baseftochan .eq. "none") then
        basefTOchan = sum(baseflow(sloc:eloc) * L1_area(sloc:eloc)) / R_areaCell(iNode)
      end if

      if (any(vn_subsrftochan(:) .eq. "none")) then
        do nn = 1, nSubsrflow
          subsrfTOchan(nn) = sum(subsrfRunoff(sloc:eloc,nn) * L1_area(sloc:eloc)) / R_areaCell(iNode)
        end do
      end if
      !conc of the runoff seepage
      do mm = 1, nsubstances
        if (fsealed_rout > 1.0E-4_dp) then
          if (dirTOchan > 0.0_dp) then
            concdirTochan(mm) = sum(directrunoff(sloc:eloc) * L1_area(sloc:eloc)*cdirectrunoff(sloc:eloc,mm)) /  &
                                      (dirTOchan*R_areaCell(iNode))
          else
            concdirTochan(mm) = const_crain(mm)
          end if
        end if
        if (srfTOchan > 0.0_dp) then
          concsrfTOchan(mm) =  sum(srfrunoff(sloc:eloc)*L1_area(sloc:eloc)*csrfrunoff(sloc:eloc,mm)) / &
                                     (srfTOchan*R_areaCell(iNode)) 
        else
          concsrfTOchan(mm) = const_crain(mm)
        end if
        if (basefTOchan > 0.0_dp) then
        concbasefTOchan(mm) =  sum(baseflow(sloc:eloc)*L1_area(sloc:eloc)*cbaseflow(sloc:eloc,mm)) / &
                                     (basefTOchan*R_areaCell(iNode))
        else
          concbasefTOchan(mm) = const_crain(mm)
        end if
        do nn = 1,nSubsrflow
          if (subsrfTOchan(nn) > 0.0_dp) then
            concsubsrfTOchan(nn,mm) = sum(subsrfRunoff(sloc:eloc,nn) * L1_area(sloc:eloc)*csubsrfRunoff(sloc:eloc,nn,mm))/&
                            (subsrfTOchan(nn)*R_areaCell(iNode))
          else
            concsubsrfTOchan(nn,mm) = const_crain(mm)
          end if
          tmp_locsubsrfinflux(mm) = tmp_locsubsrfinflux(mm) + concsubsrfTOchan(nn,mm) * subsrfTOchan(nn)
        end do
      end do
      !averaged air temperature
      avgairtemp = sum(airtemp(sloc:eloc) * L1_area(sloc:eloc))/ R_areaCell(iNode)
    end if

    !the total input flux from local terrestrial export [g/m^2]
    tmp_loctotal_influx(:) = (1-fsealed_rout)*(tmp_locsubsrfinflux(:) + &   !subsurface
                             srfTOchan * concsrfTOchan(:) + & !surface
                             basefTOchan * concbasefTOchan(:)) + &  !baseflow
                             fsealed_rout * dirTOchan * concdirTochan(:)
    totalrunoff = (1-fsealed_rout)*(sum(subsrfTOchan(:))+srfTOchan+basefTOchan) + & !first in [m]
                   fsealed_rout * dirTOchan
    if (totalrunoff>0.0_dp) then
      concOUT(:) =  tmp_loctotal_influx(:)/ totalrunoff  ![mg/l]
    end if
	
    qOUT = totalrunoff * R_areaCell(iNode) / sec_TS  ![m] -> [m^3/s]

    loadexport(:) = tmp_loctotal_influx(:)  !total export [mg/l * m, the same as other budget outputs]

    !if there are addtional inflow (e.g., point source inputs)
    if (nInflowGauges .gt. 0_i4 ) then 
      !inflow gauge ids
      do ii =1, nInflowGauges 
        if (basin_gauge%inflowgaugeid_loc(ii) == iNode) then
            if (basin_gauge%inflow_ishead(ii)) then
              concOUT(:) = basin_gauge%InflowGaugeConc(tt,ii, 1:2) 
            else
              concOUT(:) = (concOUT(:)* qOUT+ &   !terrestrial
                basin_gauge%InflowGaugeConc(tt,ii, 1:2) *basin_gauge%inflowQ(tt,ii) )/ &  !additional input
                (qOUT+basin_gauge%inflowQ(tt,ii))   
            end if
        end if
      end do 
    end if


  end subroutine terrestial_lateralconc_accumulation


  ! ------------------------------------------------------------------

  !     NAME
  !         instream_nutrient_processes

  !     PURPOSE
  !>        \brief In-stream processes in a reach(link). 

  !>        \details Calculates in-stream processes in a sepcific reach(link), including denitrification and assimilatory.

  !     LITERATURE
  !         HYPE model
  !         Yang et al., 2019 Water Res. 

  !     HISTORY
  !>        \author Xiaoqiang Yang
  !>        \date Sep 2016
  !>    Modified
  !>        X. Yang Jun 2017 enabled different time-step (hourly)
  !>        X. Yang     2019 HF data-based autotrophic N uptake description

  subroutine instream_nutrient_processes(TS,no_day,i, rivertemp11, rt_avg10, rt_avg20, benthic_area, depth, &
        riverbox, criverbox, aqdenitri, aqassimil, pardeniratew,parautouptkrate, parprodrate,nor_gr, rz_coeff, f_light)

  use mo_wqm_shadingeffect,      only: rz_shading_coeff
  use mo_wqm_global_variables,   only: GR_file_exist
  implicit none
  integer(i4),                    intent(in)     :: TS
  integer(i4),                    intent(in)     :: no_day     ! for light calculate -- shading effect
  integer(i4),                    intent(in)     :: i
  real(dp),                       intent(in)     :: rivertemp11
  real(dp),                       intent(in)     :: rt_avg10
  real(dp),                       intent(in)     :: rt_avg20
  real(dp),                       intent(in)     :: benthic_area  ![m^2]
  real(dp),                       intent(in)     :: depth         ![m]
  real(dp),                       intent(inout)  :: riverbox      ![m^3]
  real(dp), dimension(:),         intent(inout)  :: criverbox     ![mg/l] 
  real(dp),                       intent(out)  :: aqdenitri       ![kg/timestep]
  real(dp),                       intent(out)  :: aqassimil       ![kg/timestep]
  real(dp),                       intent(in)     :: pardeniratew
  real(dp),                       intent(in)     :: parautouptkrate
  real(dp),                       intent(in)     :: parprodrate  
  real(dp),                       intent(in)     :: nor_gr   ! for shading effect
  real(dp),                       intent(inout)  :: rz_coeff ! coefficient for each reach
  real(dp),                       intent(inout)  :: f_light  ! 5 days' moving average of rz_coeff

  !local
  real(dp)                 :: INpool, ONpool, tp_ave,aqassimil0
  real(dp)                 :: f_temp, f_conc, f_tp, f_temp1, f_temp2
  real(dp)                 :: f_tempm, f_netall
  real(dp)              :: DT
  
  !HYPE parameter
  real(dp), parameter      :: tpmean = 0.05_dp
  ! constant variables
  real(dp), parameter      :: halfsatINwater = 1.5  ! mg/l
  real(dp), parameter      :: maxdenitriwater = 0.5_dp !
  real(dp), parameter      :: maxprodwater = 0.5_dp  ! 
  !real(dp), parameter      :: maxdegradwater = 0.5_dp  
  real(dp), parameter      :: halfsatIPwater = 0.05_dp  
  real(dp), parameter      :: activedepth = 1.0_dp  
  
  DT = 24.0_dp / TS

 
  INpool = riverbox * criverbox(1) / 1000.0_dp   ! kg
  f_temp = tempfactor(rivertemp11)
  f_conc = criverbox(1) / (criverbox(1) + halfsatINwater )

  !denitrification amount
  aqdenitri = pardeniratew * f_temp * f_conc * benthic_area / DT
  aqdenitri = min(maxdenitriwater*INpool, aqdenitri) 

  !update pool and conc.
  INpool = INpool - aqdenitri  ![kg]
  if (riverbox >1.0E-10_dp) then
    criverbox(1) = INpool / riverbox * 1000.0_dp   !mg/l
  else
    criverbox(1) = 0.0_dp
  end if
  !primary production and mineralisation (inverse processes)
  ! for temperature effects on mineralisation
  if (rivertemp11 > 0.0_dp) then
     f_temp1 = rivertemp11 /20.0_dp
  else
     f_temp1 = 0.0_dp
  end if  
  f_temp2 = (rt_avg10 - rt_avg20) / 5.0_dp
  !the final empirical temperature effect, the same as HYPE implementation
  f_tempm = f_temp1 * f_temp2
  !get the updated stream ON pool
  ONpool = riverbox * criverbox(2) /1000.0_dp   !kg
  if (GR_file_exist) then
    !##############################################
    !##new method considering light availability ##
    !##############################################
    call rz_shading_coeff(i, no_day, nor_gr, rz_coeff)
    f_light = f_light + (rz_coeff-f_light) / 5.0_dp

    aqassimil =0.0_dp
    aqassimil0 = 0.0_dp
    
    if (rivertemp11 >= 0.0_dp) then
       aqassimil0 =  parautouptkrate* f_light * activedepth * benthic_area * 1.0E-6 /DT   ![kg]
	   ! here the parautouptkrate (potential N autotrophic uptake = 300 mg N m-2 d-1) 
       aqassimil0 = min(0.9_dp *INpool*activedepth, aqassimil0 ) 
       !net uptake via primary production	   
       !part of the gross uptake will be returned via respiration controlled by a parameter and
       !additional temperature effect (f_tempm approx. in [-0.2,0.4]) is taken from original HYPE method
       f_netall = parprodrate*f_tempm
       if (f_netall > 0.0_dp) then
          if (f_netall >= 1.0_dp) f_netall = 1.0_dp !to ensure the overall coeff <= 1.0
          aqassimil = f_netall*aqassimil0
          aqassimil = min(maxprodwater*INpool*activedepth, aqassimil)
       else
          if (f_netall <= -1.0_dp) f_netall = -1.0_dp
          aqassimil = f_netall*aqassimil0
          aqassimil = max(-maxprodwater*ONpool*activedepth, aqassimil)
       end if
       ! !mineralization (respiration) assumed part of autotrophic assimilation
       ! aqassimil = parprodrate*aqassimil0
    end if
  else
    !##method from HYPE##    
    tp_ave = tpmean   ! since phosphorus is not implemented at the monment..
    f_tp = tp_ave / (tp_ave + halfsatIPwater)
    aqassimil =0.0_dp
    if (rivertemp11 >= 0.0_dp) then
       aqassimil = parprodrate * f_tempm * f_tp * activedepth * depth * benthic_area/ DT !
       if (aqassimil > 0.0_dp) then
          aqassimil = min(maxprodwater*INpool*activedepth, aqassimil ) 
       else
          aqassimil = max(-maxprodwater*ONpool*activedepth, aqassimil) 
       end if
    end if
  end if

  
  
  INpool = INpool - aqassimil
  ONpool = ONpool + aqassimil
  if (riverbox >1.0E-10_dp) then
    criverbox(1) = INpool / riverbox * 1000.0_dp   !mg/l
    criverbox(2) = ONpool / riverbox * 1000.0_dp   !mg/l
  else
    criverbox(1) = 0.0_dp
    criverbox(2) = 0.0_dp
  end if 


  end subroutine instream_nutrient_processes 
 !--------------------------------------------------------------------------------------
 
!local subroutines and functions
!--------------------------------------------------------------------------  
  subroutine mix_conc(ws,cws,infil,cinfil, laterin, claterin)

    implicit none
    real(dp),               intent(in)      ::ws    ! soil moisture
    real(dp), dimension(:), intent(inout)   ::cws   ! conc. in the soil layer
    real(dp),               intent(in)      ::infil ! input water from the upper soil layer
    real(dp), dimension(:), intent(in)     ::cinfil ! conc. in the input soil water from upper layer
    real(dp),              optional, intent(in) ::laterin !lateral input,optional
    real(dp), dimension(:),optional, intent(in) ::claterin !conc. in the lateral input, optional	

    if (present(laterin)) then
      if ((ws + infil+laterin) > 1.0E-10_dp) then
        cws(:) = (ws * cws(:) + infil * cinfil(:) + laterin * claterin(:)) &
                 / (ws + infil+ laterin)
      end if
    else
      if ((ws + infil) > 1.0E-10_dp) then
        cws(:) = (ws * cws(:) + infil * cinfil(:)) / (ws + infil)
      end if
    end if
  end subroutine mix_conc
  !----------------------------------------------------------------------------
  subroutine add_source_water(vol, conc, source)
     implicit none
     real(dp),              intent(in)     :: vol          ![m]water volumn
     real(dp),              intent(inout)  :: conc         ![mg/l]concentration in water
     real(dp),              intent(in)     :: source       ![g/m^2]amount of source  
     !local
     !vol*conc-->[m]*[mg/l]=[g/m^2]
     !to aviod numerical error
     if (vol > 1.0E-10_dp) then
        conc = (conc * vol + source) / vol
     end if
     
  end subroutine add_source_water
  !-----------------------------------------------------------------------------
  subroutine production_pool(pool, source)
  
    implicit none
    real(dp),              intent(inout)  :: pool
    real(dp),              intent(in)     :: source

    pool = pool + source
  
  end subroutine production_pool  
  !------------------------------------------------------------------------------
  subroutine calculate_soil_temperature(airt, soilt, snowdp)  
    implicit none
    real(dp),              intent(in)   ::airt     ! air temperature
    real(dp),              intent(inout):: soilt   ! soil temperature
    real(dp),              intent(in)   :: snowdp  ! snow depth [m]
    !local
    !parameter variables, constant
    real(dp), parameter  :: temp_deep = 5.0_dp     ! temperature of deep soil
    real(dp), parameter  :: soil_mem = 30.0_dp     ! soil memory [days]
    real(dp), parameter  :: sp_frost = 10.0_dp     ! -----
    real(dp), parameter  :: w_deep = 0.001_dp  ! weight of deep soil temperature

    real(dp)  :: w_air  ! weight of air temperature

    w_air = 1.0_dp / (soil_mem + sp_frost * snowdp)
    soilt = soilt * (1.0_dp-w_air-w_deep) + airt * w_air + temp_deep * w_deep
	
  end subroutine calculate_soil_temperature
  !-------------------------------------------------------------------------------
  real(kind=dp) function tempfactor(c_temp)
  
  implicit none
  real(dp),               intent(in)  :: c_temp     ! temperature
  real(dp)                            :: f_temp     ! function value
  
  f_temp = 2**((c_temp - 20.0_dp) / 10.0_dp)
  if (c_temp < 5.0_dp) f_temp = f_temp * (c_temp / 5.0_dp)
  if (c_temp < 0.0_dp) f_temp = 0.0_dp
  tempfactor = f_temp
    
 
  end function tempfactor
  !-------------------------------------------------------------------------------
  real(kind=dp) function moistfactor(sm, wp, sat, thickness)
  implicit none
  !values of a sepcific soil layer
  real(dp),       intent(in) :: sm   !soil moisture
  real(dp),       intent(in) :: wp          !wilting point
  real(dp),       intent(in) :: sat         !saturated water content
  real(dp),       intent(in) :: thickness       ! soil thickness
  !local
  ! parameter variables constant
  real(dp), parameter        :: smf_satact = 0.6_dp  !
  real(dp), parameter        :: smf_upp    = 0.12_dp  !  
  real(dp), parameter        :: smf_low    = 0.08_dp  !  
  real(dp), parameter        :: smf_pow    = 1.0_dp  !

  real(dp)    :: f_sm  
  if (sm >= sat) then
     f_sm = smf_satact  
  elseif (sm <= wp) then
     f_sm = 0.0_dp
  else
     f_sm = min(1.0_dp, (1.0_dp-smf_satact)*((sat-sm)/ ((smf_upp/100.0_dp)*thickness))**smf_pow + smf_satact)
     f_sm = min(f_sm,((sm-wp)/((smf_low/100.0_dp)*thickness))**smf_pow)
  end if

  moistfactor = f_sm

  end function moistfactor
  !-----------------------------------------------------------------------------------
  real(kind =dp) function deni_moistfactor(sm, sat)
  implicit none 
  real(dp),          intent(in)  :: sm
  real(dp),          intent(in)  :: sat
  !constant parameter  
  real(dp), parameter         :: fsm_denilimit = 0.7_dp !2021-08 yangx changed based on the updated HYPE implementation
  real(dp), parameter         :: fsm_denipow   = 2.5_dp
  !local
  real(dp)  :: f_denism
  
  f_denism = 0.0_dp
  if (sm/sat > fsm_denilimit) then
     f_denism = (((sm/sat) - fsm_denilimit)/ (1.0-fsm_denilimit)) ** fsm_denipow  
  end if
  deni_moistfactor = f_denism

  end function deni_moistfactor
  
  !-----------------------------------------------------------------------------------
  real(kind=dp) function concfactor(conc)
  implicit none
  real(dp),      intent(in) :: conc
  !constant parameter
  real(dp), parameter         :: halfsatINsoil = 1.0_dp !-yangx 2021-07 changed from original 1 to 10!
  
  concfactor = conc / (conc + halfsatINsoil)

  end function concfactor  

END MODULE mo_water_quality
