!> \file mo_wqm_initialise.f90

!> \brief  allocating WQ global variables that has not yet allocated during the read in.

!> \details This module is: firstly to read water qualtiy input data, including initial 
!>          pool values, cropdata, agricultural management and crop rotation classes 
!>          and spatial distribution; secondly to allocate global variables for water quality modeling.\n
!> 

!> \authors Xiaoqiang Yang
!> \date Sep 2016

MODULE mo_wqm_initialise


  USE mo_kind, ONLY: i4, sp, dp

  
  IMPLICIT NONE

  PUBLIC :: wqm_variables_initalloc      ! initial and allocate variables of water quality model
  PUBLIC :: wqm_variables_default_init   ! default initial value for all variables especially for continuous running 

  ! ------------------------------------------------------------------

CONTAINS

 ! ------------------------------------------------------------------

  !     NAME
  !         wqm_variables_initalloc

  !     PURPOSE
  !>        \brief Allocates initial space for global variables of water quality model.

  !>        \details Allocates initial space for global variables of water quality model.
  !>        

  !     CALLING SEQUENCE
  !         

  !     INTENT(IN)
  !>        \param[in] 

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         
  !         

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Xiaoqiang Yang
  !>        \date Jun 2016

  subroutine wqm_variables_initalloc()

    use mo_append,              only: append
    use mo_wqm_paste,           only: append_3d
    use mo_string_utils,        only: num2str
    use mo_message,             only: message 
    use mo_wqm_constants,       only: nsubstances
    use mo_wqm_global_variables, only: &
     !basic catch info
        L1_nCells, L1_nLinks, nSoilLayers,  &
        simPer,nTstepDay,  &
     !for hydro inputs
        !process related, dim1= ncells, dim2 = ntimestep
        vn_WS_surfaceflow_ver, vn_WS_subsurfaceflow_ver, vn_WS_baseflow_ver,  &
        L1_rainfall_direct, L1_temp,L1_throughfall, L1_snowmelt, L1_infiltrate,  &  
        L1_evapCanopy, L1_surfaceRunoff, L1_baseflow, &
        !dim1= ncells, dim2 = ntimestep, dim3 = number of subsurface flow components
        nSubsurfaceflow, L1_subsurfaceRunoff,  &
        !ReservoirStore type
        L1_WSintercept, L1_WSsnowpack,L1_WSponding,&
        L1_WSsurface, L1_WSbase,  L1_WSsubsurface, L1_WSsoil, &
        L1_WSsealed, L1_Esealed, L1_directRunoff,&
        !dim1= ncells, dim2 = ntimestep, dim3 = number of soil layers
        L1_evapSoil,L1_transp,   &
        L1_reinfiltSoil, L1_reinfiltSrf, &
        !channel, dim1 = ncells, dim2 = ntimestep
        L1_WSchannel, L1_WSchannel_upin, L1_WSchannel_downout, &
        !dim1 = dim1 = ncells, dim2 = ntimestep
        L1_dirTOchan, L1_srfTOchan, L1_basefTOchan,    &
        !dim1 = dim1 = ncells, dim2 = ntimestep, dim3 = number OF SUBSURFACE runoffs
        nSubsrfToChannel,  & !total number of subsurface seepage generation
        L1_subsrfTOchan,L1_totalRunoffTochan, &
     !for WQ corresponding
     ! each variable has one new dimension for nsubstances
     ! no dimension for timesteps since all of them will be upddated for each time step
        L1_crainfall_direct, L1_cthroughfall, L1_csnowmelt, & 
        L1_csurfaceRunoff,L1_cinfilsurface, L1_cinfilSoil,L1_creinfiltSoil,L1_creturnflow, &
        L1_csubsurfaceRunoff, L1_cbaseflow, L1_cdirectRunoff, &
        !ReservoirStore type
        L1_cWSintercept, L1_cWSsnowpack,L1_cWSponding, L1_cWSsealed,&
        L1_cWSsurface, L1_cWSsubsurface, L1_cWSbase, L1_cWSsoil, &
        !channel      
        L1_concdirTOchan,L1_concsrfTOchan, L1_concbasefTOchan,    &
        !dim1 = dim1 = ncells, dim2 = number OF SUBSURFACE runoffs, dim3=nsubstances
        L1_concsubsrfTOchan, &
        vn_WS_interception, vn_WS_snowpack, vn_WS_surfaceflow,vn_WS_soil,vn_WS_subsurfaceflow, &
        vn_fromWS_surfaceflow,vn_fromWS_subsurfaceflow, vn_fromWS_baseflow,& 
        L1_pcsurface,&
        L1_pcsubsurface,L1_pcbase, & !pointer to link reservoirs
        prevstep_WSintercept,prevstep_WSponding,prevstep_WSsurface,prevstep_WSbase,prevstep_WSsubstorage, &
        prevstep_WSsoil,stat_LATsurface,stat_LATbase,stat_LATsubstorage,stat_LATsoil, &
        prevstep_baseflow,  prevstep_percol,   prevstep_baseinflow,prevstep_baseoutflow,  prevstep_WSsealed,  &
       !WQ calculation
        L1_SrfmixRatio,  &
        L1_Puptake,L1_humusN, L1_fastN, L1_dissolvedIN,              & ! INOUT NS Four different Nitrate pools in each soil layer
        L1_dissolvedON,   & ! INOUT NS conc. in each soil layer, dim2=soillayers, dim3=substances  
        L1_soiltemp, & ! INOUT NS soil temperature, variables for calculating baseflow conc.                                                & ! INOUT NX dim2=substances 
        L1_rdegradN, L1_rmineralN, L1_rdissolN,                   & ! INOUT NE1 nitrate submodel parameters
        L1_rdeniSoil, L1_rdeniAqtc,L1_ratuptkN, L1_rpprodN, &     ! INOUT NE1 nitrate submodel parameters  
        L1_cpercolate,L1_cbaseflow,   &
        prevstep_cpercol,prevstep_cbaseinflow,prevstep_cbaseout, &
        L1_soilUptakeN,L1_soilDenitri, L1_soilINfrtmanapp, L1_soilMineralN, &
        L1_avgtemp,stat_soiltemp, &
        nLink_yravg_q, L1_qOUT, L1_concOUT, L1_terXpt, & !terrestrial export from local grid
        L1_qTIN,L1_concTIN,L1_loadTIN,&
        nLink_riverbox, nLink_criverbox,nLink_rivert_avg10,nLink_rivert_avg20,nLink_rivertemp,  &
        L1_aquaticDenitri, L1_aquaticAssimil,&
        L1_rzcoeff,L1_flight,&
        out_csoilMoist,out_cdirectrunoff,out_csurfaceRunoff,out_csubsurfaceRunoff,out_cbaseflow,out_soilstoreN, out_soilUptakeN, &   
        out_soilDenitri, out_soilMineralN, out_soilINfrtmanapp,out_concdirTOchan,out_concsrfTOchan, &  
        out_concsubsrfTOchan, out_concbasefTOchan,out_concOUT,out_concTIN , &   
        out_terXpt, out_aquaticDenitri,out_aquaticAssimil    


    implicit none
    !

    !local
    integer(i4)    :: i,j
    real(dp), dimension(:),     allocatable   :: dummy_vector
    real(dp), dimension(:,:),   allocatable   :: dummy_matrix
    !real(dp), dimension(:,:,:), allocatable   :: dummy_matrix_conc
    real(dp), dimension(:,:,:), allocatable   :: dummy_matrix_3d
  

    integer(i4)    :: totaltimestep
    
!hydro inputs
   !dim1 = number if cells, dim2=number of timesteps
    totaltimestep = (simPer%julEnd - simPer%julStart + 1 ) * nTstepDay
    allocate( dummy_matrix (L1_nCells, totaltimestep))
    dummy_matrix(:,:) = 0.0_dp
    call append(L1_rainfall_direct,dummy_matrix)
    call append(L1_temp,dummy_matrix)
    call append(L1_throughfall,dummy_matrix)
    call append(L1_snowmelt,dummy_matrix)
    call append(L1_infiltrate, dummy_matrix)
    call append(L1_reinfiltSrf, dummy_matrix)
    call append(L1_evapCanopy,dummy_matrix)
    call append(L1_soiltemp, dummy_matrix )   
    call append(L1_WSsealed,dummy_matrix)
    call append(L1_Esealed,dummy_matrix)

    !the type ReservoirStore
    !primary vertical fluxes
    call append(L1_WSintercept%reservoir,dummy_matrix)
    call append(L1_WSintercept%up_in,dummy_matrix)
    call append(L1_WSintercept%down_out,dummy_matrix)
    call append(L1_WSintercept%up_out,dummy_matrix)
    !call append(L1_WSintercept%down_in,dummy_matrix)
    !call append(L1_WSintercept%lateral_in,dummy_matrix)
    !call append(L1_WSintercept%lateral_out,dummy_matrix)

    call append(L1_WSsnowpack%reservoir,dummy_matrix)
    call append(L1_WSsnowpack%up_in,dummy_matrix)
    call append(L1_WSsnowpack%down_out,dummy_matrix)
    call append(L1_WSsnowpack%up_out,dummy_matrix)
    !call append(L1_WSsnowpack%down_in,dummy_matrix)
    !call append(L1_WSsnowpack%lateral_in,dummy_matrix)
    !call append(L1_WSsnowpack%lateral_out,dummy_matrix)

    call append(L1_WSponding%reservoir,dummy_matrix)
    call append(L1_WSponding%up_in,dummy_matrix)
    call append(L1_WSponding%down_out,dummy_matrix)
    call append(L1_WSponding%up_out,dummy_matrix)
    call append(L1_WSponding%down_in,dummy_matrix)
    !call append(L1_WSponding%lateral_in,dummy_matrix)
    !call append(L1_WSponding%lateral_out,dummy_matrix)

    allocate(L1_WSsoil(nSoilLayers))
    do i = 1, nSoilLayers
        call append(L1_WSsoil(i)%reservoir,dummy_matrix)
        call append(L1_WSsoil(i)%up_in,dummy_matrix)
        call append(L1_WSsoil(i)%down_out,dummy_matrix)
        call append(L1_WSsoil(i)%up_out,dummy_matrix)
        call append(L1_WSsoil(i)%down_in,dummy_matrix)
        !call append(L1_WSsoil(i)%lateral_in,dummy_matrix)
        !call append(L1_WSsoil(i)%lateral_out,dummy_matrix)
    end do
   !processes and conceptual storages
    call append(L1_surfaceRunoff,dummy_matrix)
    call append(L1_baseflow,dummy_matrix)
    call append(L1_directRunoff,dummy_matrix)
    !conceptual storages primary lateral, with possible vertical exchange 
    call append(L1_WSsurface%reservoir,dummy_matrix)
    call append(L1_WSsurface%lateral_in,dummy_matrix)
    call append(L1_WSsurface%lateral_out,dummy_matrix)
    if (any(vn_WS_surfaceflow_ver .ne. "none")) then
        call append(L1_WSsurface%up_in,dummy_matrix)
        call append(L1_WSsurface%down_out,dummy_matrix)
        call append(L1_WSsurface%up_out,dummy_matrix)
        call append(L1_WSsurface%down_in,dummy_matrix)
    end if
    call append(L1_WSbase%reservoir,dummy_matrix)
    call append(L1_WSbase%lateral_in,dummy_matrix)
    call append(L1_WSbase%lateral_out,dummy_matrix)
    if (any(vn_WS_baseflow_ver .ne. "none")) then
        call append(L1_WSbase%up_in,dummy_matrix)
        call append(L1_WSbase%down_out,dummy_matrix)
        call append(L1_WSbase%up_out,dummy_matrix)
        call append(L1_WSbase%down_in,dummy_matrix)
    end if

    allocate(L1_WSsubsurface(nSubsurfaceflow))
    do i = 1, nSubsurfaceflow
        call append(L1_WSsubsurface(i)%reservoir,dummy_matrix)
        call append(L1_WSsubsurface(i)%lateral_in,dummy_matrix)
        call append(L1_WSsubsurface(i)%lateral_out,dummy_matrix)
        if (any(vn_WS_subsurfaceflow_ver(i,:) .ne. "none")) then
            call append(L1_WSsubsurface(i)%up_in,dummy_matrix)
            call append(L1_WSsubsurface(i)%down_out,dummy_matrix)
            call append(L1_WSsubsurface(i)%up_out,dummy_matrix)
            call append(L1_WSsubsurface(i)%down_in,dummy_matrix)
        end if
    end do
    deallocate( dummy_matrix)
    !subsurface fluxes and storages 
    allocate( dummy_matrix_3d (L1_nCells, totaltimestep,nSubsurfaceflow))
    dummy_matrix_3d(:,:,:) = 0.0_dp
    call append_3d(L1_subsurfaceRunoff,dummy_matrix_3d)
    deallocate( dummy_matrix_3d)
    !ET related
    !dim1= ncells, dim2 = ntimestep, dim3 = number of involved soil layers
    allocate( dummy_matrix_3d (L1_nCells, totaltimestep,nSoilLayers))
    dummy_matrix_3d(:,:,:) = 0.0_dp
    call append_3d(L1_evapSoil,dummy_matrix_3d)
    call append_3d(L1_transp, dummy_matrix_3d)
    call append_3d(L1_reinfiltSoil, dummy_matrix_3d)
    deallocate( dummy_matrix_3d)

    !lateral "L1_nLinks"----------------
    allocate( dummy_matrix (L1_nLinks, totaltimestep))
    dummy_matrix(:,:) = 0.0_dp
   !channel
    call append(L1_WSchannel,dummy_matrix)
    call append(L1_WSchannel_upin,dummy_matrix)
    call append(L1_WSchannel_downout,dummy_matrix)
    call append(L1_dirTOchan,dummy_matrix)
    call append(L1_srfTOchan,dummy_matrix)
    call append(L1_basefTOchan,dummy_matrix)
    deallocate( dummy_matrix)

   !channel subsurface seepage
    !dim1 = dim1 = ncells, dim2 = ntimestep, dim3 = number OF SUBSURFACE runoffs
    allocate( dummy_matrix_3d (L1_nLinks, totaltimestep,nSubsrfToChannel))
    dummy_matrix_3d(:,:,:) = 0.0_dp
    call append_3d(L1_subsrfTOchan, dummy_matrix_3d)
    deallocate( dummy_matrix_3d)

!WQ corresponding
    !dim1 = number if cells, dim2=number of substances
    allocate( dummy_matrix (L1_nCells, nsubstances))
    dummy_matrix(:,:) = 0.0_dp
    call append(L1_crainfall_direct,dummy_matrix)
    call append(L1_cthroughfall,dummy_matrix)
    call append(L1_csnowmelt,dummy_matrix)
    call append(L1_cinfilsurface, dummy_matrix )
    call append(L1_cWSsealed,dummy_matrix)

	
    !the type ReservoirStore
    call append(L1_cWSintercept%reservoir,dummy_matrix)
    call append(L1_cWSintercept%up_in,dummy_matrix)
    call append(L1_cWSintercept%down_out,dummy_matrix)
    call append(L1_cWSintercept%up_out,dummy_matrix)
    call append(L1_cWSintercept%down_in,dummy_matrix)
    call append(L1_cWSintercept%lateral_in,dummy_matrix)
    call append(L1_cWSintercept%lateral_out,dummy_matrix)

    call append(L1_cWSsnowpack%reservoir,dummy_matrix)
    call append(L1_cWSsnowpack%up_in,dummy_matrix)
    call append(L1_cWSsnowpack%down_out,dummy_matrix)
    call append(L1_cWSsnowpack%up_out,dummy_matrix)
    call append(L1_cWSsnowpack%down_in,dummy_matrix)
    call append(L1_cWSsnowpack%lateral_in,dummy_matrix)
    call append(L1_cWSsnowpack%lateral_out,dummy_matrix)
    call append(L1_cWSponding%reservoir,dummy_matrix)
    call append(L1_cWSponding%up_in,dummy_matrix)
    call append(L1_cWSponding%down_out,dummy_matrix)
    call append(L1_cWSponding%up_out,dummy_matrix)
    call append(L1_cWSponding%down_in,dummy_matrix)
    call append(L1_cWSponding%lateral_in,dummy_matrix)
    call append(L1_cWSponding%lateral_out,dummy_matrix)

    allocate(L1_cWSsoil(nSoilLayers))
    do i = 1, nSoilLayers
        call append(L1_cWSsoil(i)%reservoir,dummy_matrix)
        call append(L1_cWSsoil(i)%up_in,dummy_matrix)
        call append(L1_cWSsoil(i)%down_out,dummy_matrix)
        call append(L1_cWSsoil(i)%up_out,dummy_matrix)
        call append(L1_cWSsoil(i)%down_in,dummy_matrix)
        call append(L1_cWSsoil(i)%lateral_in,dummy_matrix)
        call append(L1_cWSsoil(i)%lateral_out,dummy_matrix)
    end do
	!processes and conceptual storages
    call append(L1_cpercolate, dummy_matrix )
    call append(L1_csurfaceRunoff,dummy_matrix)
    call append(L1_cbaseflow,dummy_matrix)
    call append(L1_cdirectRunoff,dummy_matrix)
    call append(prevstep_cpercol,dummy_matrix )
    call append(prevstep_cbaseinflow,dummy_matrix )
    call append(prevstep_cbaseout, dummy_matrix )
    call append(out_cdirectrunoff,dummy_matrix)
    call append(out_csurfaceRunoff,dummy_matrix)
    call append(out_cbaseflow,dummy_matrix)

    call append(L1_cWSsurface%reservoir,dummy_matrix)
    call append(L1_cWSsurface%lateral_in,dummy_matrix)
    call append(L1_cWSsurface%lateral_out,dummy_matrix)
    if (any(vn_WS_surfaceflow_ver .ne. "none")) then
        call append(L1_cWSsurface%up_in,dummy_matrix)
        call append(L1_cWSsurface%down_out,dummy_matrix)
        call append(L1_cWSsurface%up_out,dummy_matrix)
        call append(L1_cWSsurface%down_in,dummy_matrix)
    end if
    call append(L1_cWSbase%reservoir,dummy_matrix)
    call append(L1_cWSbase%lateral_in,dummy_matrix)
    call append(L1_cWSbase%lateral_out,dummy_matrix)
    if (any(vn_WS_baseflow_ver .ne. "none")) then
        call append(L1_cWSbase%up_in,dummy_matrix)
        call append(L1_cWSbase%down_out,dummy_matrix)
        call append(L1_cWSbase%up_out,dummy_matrix)
        call append(L1_cWSbase%down_in,dummy_matrix)
    end if
    allocate(L1_cWSsubsurface(nSubsurfaceflow))
    do i = 1, nSubsurfaceflow
        call append(L1_cWSsubsurface(i)%reservoir,dummy_matrix)
        call append(L1_cWSsubsurface(i)%lateral_in,dummy_matrix)
        call append(L1_cWSsubsurface(i)%lateral_out,dummy_matrix)
        if (any(vn_WS_subsurfaceflow_ver(i,:) .ne. "none")) then
            call append(L1_cWSsubsurface(i)%up_in,dummy_matrix)
            call append(L1_cWSsubsurface(i)%down_out,dummy_matrix)
            call append(L1_cWSsubsurface(i)%up_out,dummy_matrix)
            call append(L1_cWSsubsurface(i)%down_in,dummy_matrix)
        end if
    end do
    deallocate( dummy_matrix)

    !subsurface conc fluxes and storages 
    !dim1 = number if cells,dim2=number OF SUBSURFACE storages, dim3=number of substances
    allocate( dummy_matrix_3d (L1_nCells,nSubsurfaceflow, nsubstances))
    dummy_matrix_3d(:,:,:) = 0.0_dp
    call append_3d(L1_csubsurfaceRunoff,dummy_matrix_3d)
    call append_3d(out_csubsurfaceRunoff,dummy_matrix_3d)
    deallocate( dummy_matrix_3d)

!dim1 = number if cells, dim2=number of substances
    allocate( dummy_matrix (L1_nLinks, nsubstances))
    dummy_matrix(:,:) = 0.0_dp
    !channel concentration in runoff generation
    call append(L1_concdirTOchan,dummy_matrix)
    call append(L1_concsrfTOchan,dummy_matrix)
    call append(L1_concbasefTOchan,dummy_matrix)
    call append(out_concdirTOchan,dummy_matrix)
    call append(out_concsrfTOchan,dummy_matrix)
    call append(out_concbasefTOchan,dummy_matrix)
    deallocate( dummy_matrix)

    !conc channel subsurface seepage
    !dim1 = dim1 = ncells,  dim2=number OF SUBSURFACE runoffs, dim3 =  nsubstances
    allocate( dummy_matrix_3d (L1_nLinks,nSubsrfToChannel ,nsubstances))
    dummy_matrix_3d(:,:,:) = 0.0_dp
    call append_3d(L1_concsubsrfTOchan, dummy_matrix_3d)
    call append_3d(out_concsubsrfTOchan, dummy_matrix_3d)
    deallocate( dummy_matrix_3d)

   !detect the source storage from where conc should be taken
   !pointers for assigning source WS storages
   !fromWS_surfaceflow
    if (vn_fromWS_surfaceflow .eq. vn_WS_interception) then
      L1_pcsurface=>L1_cWSintercept
    else if (vn_fromWS_surfaceflow .eq. vn_WS_snowpack) then
      L1_pcsurface=>L1_cWSsnowpack
    else if (vn_fromWS_surfaceflow .eq. vn_WS_surfaceflow) then
      L1_pcsurface=>L1_cWSponding
    else  if (vn_fromWS_surfaceflow .eq. vn_WS_soil(1)) then
      L1_pcsurface=>L1_cWSsoil(1)
    else !if not specified as above, default N concentration from the first soil layer
      L1_pcsurface=>L1_cWSsoil(1)
      print*, 'Note: surface flow source stroage not specified, the first-layer soil water concentration is taken by default!'
    end if
   !fromWS_baseflow
    outerbs: do while (.true.)
    do j = 1, nSoilLayers
      if(vn_WS_soil(j) .eq. vn_fromWS_baseflow) then
        L1_pcbase=>L1_cWSsoil(j)
        exit outerbs
      end if
    end do
    do j = 1, nSubsurfaceflow
      if(vn_WS_subsurfaceflow(j) .eq. vn_fromWS_baseflow) then
        L1_pcbase=>L1_cWSsubsurface(j)
        exit outerbs
      end if
      if (j ==nSubsurfaceflow) then
        call message('***ERROR: no source storage assiged to baseflow storage')
        stop
      end if
    end do
    end do outerbs

   !fromWS_subsurfaceflow
    allocate(L1_pcsubsurface(nSubsurfaceflow))
    outer: do i =1,nSubsurfaceflow            
        do j = 1, nSoilLayers 
          if (vn_WS_soil(j) .eq. vn_fromWS_subsurfaceflow(i)) then
            L1_pcsubsurface(i)%Pnter_Rsvr=>L1_cWSsoil(j)
            exit outer
          end if
        end do
        do j = 1, nSubsurfaceflow
          if (vn_WS_subsurfaceflow(j) .eq. vn_fromWS_subsurfaceflow(i)) then
            L1_pcsubsurface(i)%Pnter_Rsvr=>L1_cWSsubsurface(j)
            if (i == j) then
              call message('***ERROR: subsurface storage sourced to itself')
              stop 
            end if
            exit outer
          end if
        end do
    end do outer 

   ! !fromWS_srftochan
    ! if (vn_fromWS_srftochan .eq. vn_WS_interception) then
      ! L1_pcsrfTOchan=>L1_cWSintercept
    ! else if (vn_fromWS_srftochan .eq. vn_WS_snowpack) then
      ! L1_pcsrfTOchan=>L1_cWSsnowpack
    ! else if (vn_fromWS_srftochan .eq. vn_WS_surfaceflow) then
      ! L1_pcsrfTOchan=>L1_cWSsurface
    ! else
      ! outersrfchn: do while (.true.)
        ! do j = 1, nSoilLayers 
          ! if (vn_WS_soil(j) .eq. vn_fromWS_srftochan) then
            ! L1_pcsrfTOchan=>L1_cWSsoil(j)
            ! exit outersrfchn
          ! end if
        ! end do
        ! do j = 1, nSubsurfaceflow
          ! if (vn_WS_subsurfaceflow(j) .eq. vn_fromWS_srftochan) then
            ! L1_pcsrfTOchan=>L1_cWSsubsurface(j)
            ! exit outersrfchn
          ! end if
          ! if (j ==nSubsurfaceflow) then
            ! call message('***ERROR: no source storage assiged to srfTOchn storage')
            ! stop
          ! end if
        ! end do
      ! end do outersrfchn
    ! end if
   ! !fromWS_baseftochan
    ! outerbsfchn: do while (.true.)
    ! do j = 1, nSoilLayers
      ! if(vn_WS_soil(j) .eq. vn_fromWS_baseftochan) then
        ! L1_pcbasefTOchan=>L1_cWSsoil(j)
        ! exit outerbsfchn
      ! end if
    ! end do
    ! do j = 1, nSubsurfaceflow
      ! if(vn_WS_subsurfaceflow(j) .eq. vn_fromWS_baseftochan) then
        ! L1_pcbasefTOchan=>L1_cWSsubsurface(j)
        ! exit outerbsfchn
      ! end if
      ! if (j ==nSubsurfaceflow) then
        ! call message('***ERROR: no source storage assiged to bsfTOchn storage')
        ! stop
      ! end if
    ! end do
    ! end do outerbsfchn
   ! !specifically for arrays

   ! !fromWS_subsrfTOchan
    ! allocate(L1_pcsubsrfTOchan(nSubsrfToChannel))
    ! outersubsrfchn: do i =1,nSubsrfToChannel           
        ! do j = 1, nSoilLayers 
            ! if (vn_WS_soil(j) .eq. vn_fromWS_subsrftochan(i)) then
              ! L1_pcsubsrfTOchan(i)%Pnter_Rsvr=>L1_cWSsoil(j)
              ! exit outersubsrfchn
            ! end if
        ! end do
        ! do j = 1, nSubsurfaceflow
            ! if (vn_WS_subsurfaceflow(j) .eq. vn_fromWS_subsrftochan(i)) then
              ! L1_pcsubsrfTOchan(i)%Pnter_Rsvr=>L1_cWSsubsurface(j)
              ! exit outersubsrfchn
            ! end if
            ! if (j ==nSubsurfaceflow) then
              ! call message('***ERROR: no source storage assiged to subsrfTOchn storage')
              ! stop
            ! end if
        ! end do
    ! end do outersubsrfchn  

    allocate( dummy_vector(L1_nCells))
    allocate( dummy_matrix(L1_nCells, nSoilLayers))
    dummy_vector(:) = 0.0_dp
    dummy_matrix(:,:) = 0.0_dp
    !parameters
    call append(L1_SrfmixRatio, dummy_vector ) !
    call append(L1_rdegradN, dummy_vector )
    call append(L1_rmineralN, dummy_vector )
    call append(L1_rdissolN, dummy_vector )
    call append(L1_rdeniSoil, dummy_vector )
    call append(prevstep_WSponding, dummy_vector)
    call append(prevstep_WSintercept, dummy_vector)
    call append(prevstep_WSsealed, dummy_vector)
    call append(stat_soiltemp, dummy_vector)
 
    call append(L1_soilUptakeN, dummy_vector)
    call append(L1_soilDenitri, dummy_vector)
    call append(L1_soilMineralN, dummy_vector)
    call append(L1_soilINfrtmanapp, dummy_vector)
    call append(out_soilUptakeN, dummy_vector)
    call append(out_soilDenitri, dummy_vector)
    call append(out_soilMineralN, dummy_vector)
    call append(out_soilINfrtmanapp, dummy_vector)


    call append(prevstep_WSsurface,dummy_vector)
    call append(stat_LATsurface,dummy_vector)
    call append(prevstep_WSbase,dummy_vector)
    call append(stat_LATbase,dummy_vector)
    call append(prevstep_baseflow,dummy_vector)
    call append(prevstep_percol,dummy_vector)
    call append(prevstep_baseinflow,dummy_vector)
    call append(prevstep_baseoutflow,dummy_vector)
    deallocate( dummy_vector)
    deallocate( dummy_matrix)
	!nLinks
    allocate( dummy_vector(L1_nLinks))
    allocate( dummy_matrix(L1_nLinks, nsubstances))
    dummy_vector(:) = 0.0_dp
    dummy_matrix(:,:) = 0.0_dp
    call append(nLink_rivertemp, dummy_vector)!store 20day's average air temperature (for instream calcualting)
    call append(nLink_yravg_q, dummy_vector)
    call append(L1_qOUT, dummy_vector)
    call append(L1_qTIN, dummy_vector)
    call append(nLink_riverbox, dummy_vector)
    call append(nLink_rivert_avg10, dummy_vector)
    call append(nLink_rivert_avg20, dummy_vector)
    call append(L1_aquaticDenitri, dummy_vector)
    call append(L1_aquaticAssimil, dummy_vector)
    call append(L1_rzcoeff, dummy_vector)
    call append(L1_flight, dummy_vector)
    call append(L1_avgtemp, dummy_vector)
    call append(L1_totalRunoffTochan,dummy_vector)
    call append(L1_rdeniAqtc, dummy_vector )
    call append(L1_ratuptkN, dummy_vector )
    call append(L1_rpprodN, dummy_vector )

    call append(nLink_criverbox, dummy_matrix )
    call append(L1_concOUT, dummy_matrix )
    call append(L1_terXpt, dummy_matrix )
    call append(L1_loadTIN, dummy_matrix )
    call append(L1_concTIN, dummy_matrix )

    call append(out_concOUT, dummy_matrix )
    call append(out_terXpt, dummy_matrix )
    call append(out_concTIN, dummy_matrix )
    call append(out_aquaticDenitri, dummy_vector)
    call append(out_aquaticAssimil, dummy_vector)

    deallocate( dummy_vector)
    deallocate( dummy_matrix)

    allocate( dummy_matrix (L1_nCells, nSubsurfaceflow))
    dummy_matrix(:,:) = 0.0_dp
    call append(prevstep_WSsubstorage, dummy_matrix)
    call append(stat_LATsubstorage, dummy_matrix)
    deallocate( dummy_matrix)
    allocate( dummy_matrix (L1_nCells, nSoilLayers))
    dummy_matrix(:,:) = 0.0_dp
    call append(prevstep_WSsoil, dummy_matrix)
    call append(stat_LATsoil, dummy_matrix)
    deallocate( dummy_matrix)

    allocate( dummy_matrix_3d(L1_nCells, nSoilLayers, nsubstances))
    dummy_matrix_3d(:,:,:) = 0.0_dp
    call append_3d(L1_cinfilSoil, dummy_matrix_3d )
    call append_3d(L1_creturnflow, dummy_matrix_3d )
    call append_3d(L1_creinfiltSoil, dummy_matrix_3d )
    call append_3d(out_csoilMoist, dummy_matrix_3d )
    deallocate( dummy_matrix_3d)


    ! nitrate pools
    allocate( dummy_matrix (L1_nCells, nSoilLayers))
    dummy_matrix(:,:) = 0.0_dp
    call append(L1_Puptake, dummy_matrix )
    call append(L1_humusN, dummy_matrix )
    call append(L1_fastN, dummy_matrix )
    call append(L1_dissolvedIN, dummy_matrix )
    call append(L1_dissolvedON, dummy_matrix )
    call append(out_soilstoreN, dummy_matrix)
    deallocate( dummy_matrix)

    ! free space
    if ( allocated( dummy_vector          ) ) deallocate( dummy_vector          )
    if ( allocated( dummy_matrix          ) ) deallocate( dummy_matrix          )
    if ( allocated( dummy_matrix_3d     ) )   deallocate( dummy_matrix_3d     )

   
	
  end subroutine wqm_variables_initalloc
 ! ------------------------------------------------------------------
 ! ------------------------------------------------------------------

  !     NAME
  !         wqm_variables_default_init

  !     PURPOSE
  !>        \brief initialise global variables of water quality model for continuous running (e.g. calibration, sensitivity and uncertainty analysis, ...).

  !>        \details Allocates initial space for global variables of water quality model.
  !>        

  !     CALLING SEQUENCE
  !         

  !     INTENT(IN)
  !>        \param[in] 

  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None

  !     INTENT(IN), OPTIONAL
  !         None

  !     INTENT(INOUT), OPTIONAL
  !         None

  !     INTENT(OUT), OPTIONAL
  !         None

  !     RETURN
  !>        None

  !     RESTRICTIONS
  !         None

  !     EXAMPLE
  !         
  !         

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Xiaoqiang Yang
  !>        \date Jun 2016

  subroutine wqm_variables_default_init()

    use mo_wqm_global_variables, only: &
        !dim1= ncells, dim2 = ntimestep, dim3 = number of subsurface flow components
        nSubsurfaceflow,nSoilLayers, &
     !for WQ corresponding
     ! each variable has one new dimension for nsubstances
     ! no dimension for timesteps since all of them will be upddated for each time step
        L1_crainfall_direct, L1_cthroughfall, L1_csnowmelt, & 
        L1_csurfaceRunoff, L1_cbaseflow, L1_cWSsealed, L1_cdirectRunoff,&
        L1_csubsurfaceRunoff, L1_cinfilsurface, L1_cinfilSoil,L1_creinfiltSoil,L1_creturnflow, &
        !ReservoirStore type
        L1_cWSintercept, L1_cWSsnowpack,L1_cWSponding, &
        L1_cWSsurface, L1_cWSsubsurface, L1_cWSbase, L1_cWSsoil,&
        vn_WS_surfaceflow_ver,vn_WS_subsurfaceflow_ver, vn_WS_baseflow_ver,   &
        !channel      
        L1_concdirTOchan,L1_concsrfTOchan, L1_concbasefTOchan,    &
        !dim1 = dim1 = ncells, dim2 = ntimestep, dim3 = number OF SUBSURFACE runoffs, dim4=nsubstances
        L1_concsubsrfTOchan, &        
        prevstep_WSintercept,prevstep_WSsurface,prevstep_WSbase,prevstep_WSsubstorage,prevstep_WSsoil, &
        stat_LATsurface,stat_LATbase,stat_LATsubstorage,stat_LATsoil, &
        prevstep_baseflow,  prevstep_percol,   prevstep_baseinflow,prevstep_baseoutflow,  prevstep_WSsealed,&
       !WQ calculation
        L1_SrfmixRatio, &
        L1_Puptake,L1_humusN, L1_fastN, L1_dissolvedIN,              & ! INOUT NS Four different Nitrate pools in each soil layer
        L1_dissolvedON,   & ! INOUT NS conc. in each soil layer, dim2=soillayers, dim3=substances  
        L1_soiltemp, & ! INOUT NS soil temperature, variables for calculating baseflow conc.                                                & ! INOUT NX dim2=substances 
        L1_rdegradN, L1_rmineralN, L1_rdissolN,                   & ! INOUT NE1 nitrate submodel parameters
        L1_rdeniSoil, L1_rdeniAqtc,L1_ratuptkN, L1_rpprodN, &     ! INOUT NE1 nitrate submodel parameters  
        L1_cpercolate,L1_cbaseflow,   &
        prevstep_cpercol,prevstep_cbaseinflow,prevstep_cbaseout, &
        L1_soilUptakeN,L1_soilDenitri, L1_soilINfrtmanapp, L1_soilMineralN, &
        L1_totalRunoffTochan,L1_avgtemp, &
        nLink_yravg_q, L1_qOUT, L1_concOUT, L1_terXpt, & !terrestrial export from local grid
        L1_qTIN,L1_concTIN,L1_loadTIN,&
        nLink_riverbox, nLink_criverbox,nLink_rivert_avg10,nLink_rivert_avg20,nLink_rivertemp,  &
        L1_aquaticDenitri, L1_aquaticAssimil,&
        L1_rzcoeff,L1_flight

    use mo_wqm_constants,  only: P1_InitStateFluxes
		 
    implicit none
    integer(i4)    :: i
    L1_crainfall_direct = P1_InitStateFluxes
    L1_cthroughfall = P1_InitStateFluxes
    L1_csnowmelt = P1_InitStateFluxes
    L1_csurfaceRunoff = P1_InitStateFluxes
    L1_cbaseflow = P1_InitStateFluxes
    L1_cinfilsurface= P1_InitStateFluxes
    L1_cpercolate= P1_InitStateFluxes
    prevstep_cpercol= P1_InitStateFluxes
    prevstep_cbaseinflow= P1_InitStateFluxes
    prevstep_cbaseout= P1_InitStateFluxes
    L1_csubsurfaceRunoff = P1_InitStateFluxes
    L1_cWSsealed = P1_InitStateFluxes
    L1_cdirectRunoff = P1_InitStateFluxes

    !channel
    L1_concdirTOchan = P1_InitStateFluxes
    L1_concsrfTOchan = P1_InitStateFluxes
    L1_concbasefTOchan = P1_InitStateFluxes
    L1_concsubsrfTOchan = P1_InitStateFluxes
    !ReservoirStore type
    L1_cWSintercept%reservoir = P1_InitStateFluxes
    L1_cWSintercept%up_in = P1_InitStateFluxes
    L1_cWSintercept%down_out = P1_InitStateFluxes
    L1_cWSintercept%up_out = P1_InitStateFluxes
    L1_cWSintercept%down_in = P1_InitStateFluxes
    L1_cWSintercept%lateral_in = P1_InitStateFluxes
    L1_cWSintercept%lateral_out = P1_InitStateFluxes
 
    L1_cWSsnowpack%reservoir = P1_InitStateFluxes
    L1_cWSsnowpack%up_in = P1_InitStateFluxes
    L1_cWSsnowpack%down_out = P1_InitStateFluxes
    L1_cWSsnowpack%up_out = P1_InitStateFluxes
    L1_cWSsnowpack%down_in = P1_InitStateFluxes
    L1_cWSsnowpack%lateral_in = P1_InitStateFluxes
    L1_cWSsnowpack%lateral_out = P1_InitStateFluxes

    L1_cWSponding%reservoir = P1_InitStateFluxes
    L1_cWSponding%up_in = P1_InitStateFluxes
    L1_cWSponding%down_out = P1_InitStateFluxes
    L1_cWSponding%up_out = P1_InitStateFluxes
    L1_cWSponding%down_in = P1_InitStateFluxes
    L1_cWSponding%lateral_in = P1_InitStateFluxes
    L1_cWSponding%lateral_out = P1_InitStateFluxes

    L1_cWSsurface%reservoir = P1_InitStateFluxes
    if (any(vn_WS_surfaceflow_ver .ne. "none")) then
      L1_cWSsurface%up_in = P1_InitStateFluxes
      L1_cWSsurface%down_out = P1_InitStateFluxes
      L1_cWSsurface%up_out = P1_InitStateFluxes
      L1_cWSsurface%down_in = P1_InitStateFluxes
    end if
    L1_cWSsurface%lateral_in = P1_InitStateFluxes
    L1_cWSsurface%lateral_out = P1_InitStateFluxes

    L1_cWSbase%reservoir = P1_InitStateFluxes
    if (any(vn_WS_baseflow_ver .ne. "none")) then
      L1_cWSbase%up_in = P1_InitStateFluxes
      L1_cWSbase%down_out = P1_InitStateFluxes
      L1_cWSbase%up_out = P1_InitStateFluxes
      L1_cWSbase%down_in = P1_InitStateFluxes
    end if
    L1_cWSbase%lateral_in = P1_InitStateFluxes
    L1_cWSbase%lateral_out = P1_InitStateFluxes
    do i = 1, nSubsurfaceflow
        L1_cWSsubsurface(i)%reservoir = P1_InitStateFluxes
        if (any(vn_WS_subsurfaceflow_ver(i,:) .ne. "none")) then
          L1_cWSsubsurface(i)%up_in = P1_InitStateFluxes
          L1_cWSsubsurface(i)%down_out = P1_InitStateFluxes
          L1_cWSsubsurface(i)%up_out = P1_InitStateFluxes
          L1_cWSsubsurface(i)%down_in = P1_InitStateFluxes
        end if
        L1_cWSsubsurface(i)%lateral_in = P1_InitStateFluxes
        L1_cWSsubsurface(i)%lateral_out = P1_InitStateFluxes
    end do
    do i = 1, nSoilLayers
        L1_cWSsoil(i)%reservoir = P1_InitStateFluxes
        L1_cWSsoil(i)%up_in = P1_InitStateFluxes
        L1_cWSsoil(i)%down_out = P1_InitStateFluxes
        L1_cWSsoil(i)%up_out = P1_InitStateFluxes
        L1_cWSsoil(i)%down_in = P1_InitStateFluxes
        L1_cWSsoil(i)%lateral_in = P1_InitStateFluxes
        L1_cWSsoil(i)%lateral_out = P1_InitStateFluxes
    end do
    prevstep_WSsealed = P1_InitStateFluxes 
    prevstep_WSintercept = P1_InitStateFluxes     ! sealed storage in last step
    prevstep_WSsurface =P1_InitStateFluxes       ! unsaturated storage in last step
    prevstep_WSbase=P1_InitStateFluxes
    stat_LATsurface =P1_InitStateFluxes         ! saturated storage in last step
    stat_LATbase=P1_InitStateFluxes
    prevstep_WSsubstorage =P1_InitStateFluxes   
    stat_LATsubstorage =P1_InitStateFluxes   
    prevstep_WSsoil=P1_InitStateFluxes   
    stat_LATsoil =P1_InitStateFluxes   
    L1_cinfilSoil =P1_InitStateFluxes
    L1_creinfiltSoil =P1_InitStateFluxes
    L1_creturnflow =P1_InitStateFluxes   
    prevstep_percol = P1_InitStateFluxes
    prevstep_baseflow = P1_InitStateFluxes
    prevstep_baseinflow = P1_InitStateFluxes          ! soil moisture in last step
    prevstep_baseoutflow= P1_InitStateFluxes   


    nLink_rivertemp= P1_InitStateFluxes
    nLink_yravg_q= P1_InitStateFluxes
    L1_qOUT= P1_InitStateFluxes
    L1_qTIN= P1_InitStateFluxes
    nLink_riverbox= P1_InitStateFluxes
    nLink_rivert_avg10= P1_InitStateFluxes
    nLink_rivert_avg20= P1_InitStateFluxes
    L1_aquaticDenitri= P1_InitStateFluxes
    L1_aquaticAssimil= P1_InitStateFluxes
    L1_rzcoeff= P1_InitStateFluxes
    L1_flight= P1_InitStateFluxes
    L1_avgtemp= P1_InitStateFluxes

    nLink_criverbox= P1_InitStateFluxes
    L1_concOUT= P1_InitStateFluxes
    L1_terXpt= P1_InitStateFluxes
    L1_loadTIN= P1_InitStateFluxes
    L1_concTIN= P1_InitStateFluxes

    L1_SrfmixRatio = P1_InitStateFluxes

    !nitrate submodel parameters
    L1_rdegradN = P1_InitStateFluxes
    L1_rmineralN = P1_InitStateFluxes
    L1_rdissolN = P1_InitStateFluxes
    L1_rdeniSoil =P1_InitStateFluxes
    L1_rdeniAqtc=P1_InitStateFluxes
    L1_ratuptkN=P1_InitStateFluxes
    L1_rpprodN=P1_InitStateFluxes
    L1_totalRunoffTochan=P1_InitStateFluxes
	
    L1_soiltemp = P1_InitStateFluxes
    L1_soilUptakeN = P1_InitStateFluxes
    L1_soilDenitri = P1_InitStateFluxes
    L1_soilMineralN = P1_InitStateFluxes
    L1_soilINfrtmanapp = P1_InitStateFluxes
    ! nitrate pools
    L1_Puptake = P1_InitStateFluxes
    L1_humusN = P1_InitStateFluxes
    L1_fastN =P1_InitStateFluxes
    L1_dissolvedIN = P1_InitStateFluxes
    L1_dissolvedON = P1_InitStateFluxes
	
  end subroutine wqm_variables_default_init
 
END MODULE mo_wqm_initialise
