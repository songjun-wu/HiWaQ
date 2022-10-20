!> \file mo_wqm_write.f90

!> \brief write out water quality model information

!> \details 

!> \author Xiaoqiang Yang
!> \date  Sep 2017

MODULE mo_wqm_write


  USE mo_kind, ONLY: i4, sp, dp

  IMPLICIT NONE


  PUBLIC :: wqm_write   ! Constant Pi in single precision
  PUBLIC :: wqm_write_output_fluxes   !State variables and Fluxes
  PUBLIC :: write_nml_par !wirte out the optimized parameter file (namelist)


CONTAINS


  subroutine write_nml_par()

    USE mo_wqm_global_variables,  ONLY : wqm_param,nTotalParam,dir_Out, &
         wqm_parameterization
    use mo_nml,                   only: open_nml,position_nml
    use mo_file,                  only: file_param_wqm, unamelist_param
    use mo_string_utils,           only: num2str
    use mo_message,                only: message
    type(wqm_parameterization), dimension(nTotalParam)   :: Nparameters

    character(256) :: fName
    integer(i4)    :: err
    !read in the par nml first
    namelist /parameters/ nTotalParam, Nparameters

    call open_nml(file_param_wqm, unamelist_param, quiet =.true.)
    
    call position_nml('parameters', unamelist_param)
    read(unamelist_param, nml=parameters)
    
    Nparameters(1:nTotalParam) = wqm_param(1:nTotalParam)
    !
    fName = trim(adjustl(dir_Out)) // trim(adjustl("optimized_parameter.nml"))
    open(89, file=trim(fName), status='unknown', action='write', iostat=err)
    if( err .ne. 0 ) then
       call message ('  IOError while opening ',trim(fName))
       call message ('  Error-Code ', num2str(err))
       stop
    end if

    write(89,nml=parameters)

  end subroutine write_nml_par

  ! ------------------------------------------------------------------

  !     NAME
  !         wqm_write_output_fluxes

  !     PURPOSE
  !>        \brief write fluxes to netcdf output files
  !
  !>        \details This subroutine creates a netcdf data set
  !>        for writing water quality model related state variables and fluxes for different time averages.
  !
  !     INTENT(IN)
  !         None
  !
  !     INTENT(INOUT)
  !         None

  !     INTENT(OUT)
  !         None
  !
  !     INTENT(IN), OPTIONAL
  !         None
  !
  !     INTENT(INOUT), OPTIONAL
  !         None
  !
  !     INTENT(OUT), OPTIONAL
  !         None
  !
  !     RETURN
  !         None
  !
  !     RESTRICTIONS
  !         None
  !
  !     EXAMPLE
  !       None
  !
  !     LITERATURE
  !       None

  !     HISTORY
  !>        \author Xiaoqiang Yang
  !>        \date Aug 2021


  subroutine wqm_write_output_fluxes( &
       timeStep_model_outputs, & ! timestep of model outputs
       warmingDays, & ! number of warming days
       newTime, & ! julian date of next time step
       nTimeSteps, & ! number of total timesteps
       nTStepDay, & ! number of timesteps per day
       tt, & ! current model timestep
       day, & ! current day of the year
       month, & ! current month of the year
       year) ! current year

    use mo_kind, only: i4, dp
    use mo_julian, only: caldat
    use mo_wqm_global_variables, only: L1_cWSsoil, L1_nCells,nSoilLayers, nsubstances,&
          L1_cdirectRunoff,L1_csurfaceRunoff, L1_csubsurfaceRunoff, &
          L1_cbaseflow, L1_concdirTOchan, L1_concsrfTOchan, L1_concsubsrfTOchan,L1_concbasefTOchan,&
          L1_concOUT,L1_concTIN, L1_terXpt,L1_dissolvedIN,  &
          L1_soilUptakeN, L1_soilDenitri, L1_soilMineralN,L1_soilINfrtmanapp,   &
          L1_aquaticDenitri,L1_aquaticAssimil, &
          !variables for ncoutput
          out_counter, nouter, &
          out_csoilMoist            , &  !soil moisture conc
          out_cdirectrunoff,        & !direct runoff conc
          out_csurfaceRunoff    , &  !surface flow conc
          out_csubsurfaceRunoff  , &  !subsurface flow conc
          out_cbaseflow          , &  !baseflow conc
          out_soilstoreN,     &       !soil N stock at each layer
          out_soilUptakeN  , &        !soil plant/crop actual N uptake
          out_soilDenitri  , &        !soil denitrification
          out_soilMineralN , &        !soil mineralization
          out_soilINfrtmanapp, &      !N external input (fertilizer, residues)
          out_concdirTOchan      , &  !channel seepage direct runoff from sealed areas
          out_concsrfTOchan      , &  !channel seepage surface
          out_concsubsrfTOchan   , &  !channel seepage subsurface
          out_concbasefTOchan    , &  !channel seepage baseflow
          out_concOUT            , &  !total runoff generation to channel
          out_concTIN         , &    !stream water Nitrate-N concentration
          out_terXpt          , &    !terrestrial export load
          out_aquaticDenitri  , &    !in-stream denitrification
          out_aquaticAssimil         !in-stream assimilatory uptake          
    use mo_message,                only: message
    use mo_string_utils,           only: num2str
    use mo_wqm_netcdfinoutputs, only: ncfile_create,ncfile_update

    implicit none
    ! input variables

    integer(i4), intent(in) :: timeStep_model_outputs
    integer(i4), intent(in) :: warmingDays
    real(dp), intent(in) :: newTime
    integer(i4), intent(in) :: nTimeSteps
    integer(i4), intent(in) :: nTStepDay
    integer(i4), intent(in) :: tt
    integer(i4), intent(in) :: day
    integer(i4), intent(in) :: month
    integer(i4), intent(in) :: year
 

    ! local variables
    integer(i4)          :: tIndex_out, tTotal_out
    logical              :: writeout
    integer(i4)          :: new_year ! year of next timestep (newTime)
    integer(i4)          :: new_month ! month of next timestep (newTime)
    integer(i4)          :: new_day ! day of next timestep (newTime)
    integer(i4)          :: isoily
    real(dp), dimension(L1_nCells,nSoilLayers,nsubstances)   :: csoilMoist

    !

    call caldat(int(newTime), yy=new_year, mm=new_month, dd=new_day)

    tIndex_out = (tt-warmingDays*nTStepDay) ! tt if write out of warming period
    tTotal_out = (nTimeSteps-warmingDays*nTStepDay) !
    ! create output dataset	   

    if ( tIndex_out == 1_i4 ) then
       call ncfile_create()
       out_counter = 1_i4
       nouter = 0_i4
    end if

    do isoily=1,nSoilLayers
        csoilMoist(:,isoily,1) = L1_cWSsoil(isoily)%reservoir(:, 1) !IN
        csoilMoist(:,isoily,2) = L1_cWSsoil(isoily)%reservoir(:, 2) !ON
    end do
    ! update Dataset
    out_csoilMoist = out_csoilMoist + csoilMoist   !soil moisture conc
    out_cdirectrunoff = out_cdirectrunoff+ L1_cdirectRunoff
    out_csurfaceRunoff = out_csurfaceRunoff + L1_csurfaceRunoff    !surface flow conc
    out_csubsurfaceRunoff = out_csubsurfaceRunoff + L1_csubsurfaceRunoff    !subsurface flow conc
    out_cbaseflow = out_cbaseflow + L1_cbaseflow !baseflow conc
    out_soilstoreN = out_soilstoreN + L1_dissolvedIN
    out_soilUptakeN = out_soilUptakeN + L1_soilUptakeN        !soil plant/crop actual N uptake
    out_soilDenitri = out_soilDenitri + L1_soilDenitri         !soil denitrification
    out_soilMineralN = out_soilMineralN + L1_soilMineralN        !soil mineralization
    out_soilINfrtmanapp = out_soilINfrtmanapp+ L1_soilINfrtmanapp     !N external input (fertilizer, residues)
    out_concdirTOchan = out_concdirTOchan + L1_concdirTOchan
    out_concsrfTOchan = out_concsrfTOchan+L1_concsrfTOchan        !channel seepage surface
    out_concsubsrfTOchan = out_concsubsrfTOchan+L1_concsubsrfTOchan     !channel seepage subsurface
    out_concbasefTOchan = out_concbasefTOchan+ L1_concbasefTOchan     !channel seepage baseflow
    out_concOUT = out_concOUT+L1_concOUT             !total runoff generation to channel
    out_concTIN = out_concTIN + L1_concTIN             !stream water Nitrate-N concentration
    out_terXpt = out_terXpt + L1_terXpt             !terrestrial export load
    out_aquaticDenitri = out_aquaticDenitri+ L1_aquaticDenitri      !in-stream denitrification
    out_aquaticAssimil = out_aquaticAssimil+L1_aquaticAssimil      !in-stream assimilatory uptake
    out_counter = out_counter + 1_i4
                
    ! determine write flag
    writeout = .false.

    if (timeStep_model_outputs .gt. 0) then
      if ((mod(tIndex_out, timeStep_model_outputs) .eq. 0) .or. (tt .eq. nTimeSteps)) then 
         writeout = .true.
         nouter = nouter + 1_i4
      end if
    else
      select case(timeStep_model_outputs)
      case(0) ! only at last time step
        if (tt .eq. nTimeSteps) then
          writeout = .true.
          nouter = nouter + 1_i4
        end if
      case(-1) ! daily
        if (((tIndex_out .gt. 1) .and. (day .ne. new_day)) .or. (tt .eq. nTimeSteps))  then
          writeout = .true.
          nouter = nouter + 1_i4
        end if
      case(-2) ! monthly
        if (((tIndex_out .gt. 1) .and. (month .ne. new_month)) .or. (tt .eq. nTimeSteps)) then
          writeout = .true.
          nouter = nouter + 1_i4
        end if
      case(-3) ! yearly
        if (((tIndex_out .gt. 1) .and. (year .ne. new_year)) .or. (tt .eq. nTimeSteps)) then
          writeout = .true.
          nouter = nouter + 1_i4
        end if
      case default ! no output at all
        continue
      end select
    endif



    ! write data
    if (writeout) then
      call ncfile_update(nouter,tIndex_out,tTotal_out)
      out_counter = 1_i4
      out_csoilMoist = 0.0_dp
      out_cdirectrunoff = 0.0_dp
      out_csurfaceRunoff = 0.0_dp
      out_csubsurfaceRunoff = 0.0_dp
      out_cbaseflow = 0.0_dp
      out_soilstoreN = 0.0_dp
      out_soilUptakeN = 0.0_dp
      out_soilDenitri = 0.0_dp
      out_soilMineralN = 0.0_dp
      out_soilINfrtmanapp = 0.0_dp
      out_concdirTOchan = 0.0_dp
      out_concsrfTOchan = 0.0_dp
      out_concsubsrfTOchan = 0.0_dp
      out_concbasefTOchan = 0.0_dp
      out_concOUT = 0.0_dp
      out_concTIN = 0.0_dp
      out_terXpt = 0.0_dp
      out_aquaticDenitri = 0.0_dp
      out_aquaticAssimil = 0.0_dp
    end if

    
  end subroutine wqm_write_output_fluxes
  
  
  
  ! ------------------------------------------------------------------

  !     NAME
  !         wqm_write

  !     PURPOSE
  !>        \brief 

  !>        \details 

  !     CALLING SEQUENCE
  !         None

  !     INTENT(IN)
  !>        None

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
  !         None
  !         

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Xiaoqiang Yang - Modified from the corresponding nHM routines
  !>        \date Sep 2016 

  subroutine wqm_write()

    use mo_wqm_global_variables,    only: &
          WQM_nutrient,Hydro_DIS, basin_gauge,nWQGauges, nHydroGauges, &
          evalPer,simPer,nTstepDay,warmingDays
    use mo_wqm_constants,           only: maxnWQinputCols  ! 3 is fixed for N simulation

    implicit none
    integer(i4)     :: ii,gg
    integer(i4)     :: iDay, tt, iS, iE,intertt
    integer(i4)     :: nTimeSteps
    real(dp), dimension(:,:,:), allocatable   :: d_ConcMod
    real(dp), dimension(:,:), allocatable     :: d_QMod

    !nutrient
    ii = evalPer%julEnd - evalPer%julStart + 1
    allocate( d_ConcMod(ii, nWQGauges, maxnWQinputCols) )
    allocate( d_QMod(ii, nHydroGauges) )
    d_ConcMod = 0.0_dp
    !discharge
    d_QMod = 0.0_dp
    !time step till the end of the evaluation period, including the warming period
    nTimeSteps = ( evalPer%julEnd - simPer%julStart + 1 ) * nTstepDay  
    iDay = 0

    print*, size(d_ConcMod)
    print*, size(WQM_nutrient)

    ! loop over timesteps
    do tt = warmingDays*nTstepDay+1, nTimeSteps, nTstepDay
          iS = tt
          iE = tt + nTstepDay - 1
          iDay = iDay + 1
          ! over gauges WQ
          do gg = 1, nWQGauges
             do intertt=iS, iE
             d_ConcMod(iDay, basin_gauge%wqm_gaugeid(gg),:) = d_ConcMod(iDay, basin_gauge%wqm_gaugeid(gg),:) &
                  + WQM_nutrient(intertt, basin_gauge%wqm_gaugeid(gg),:)
             end do
             d_ConcMod(iDay, basin_gauge%wqm_gaugeid(gg),:)=d_ConcMod(iDay, basin_gauge%wqm_gaugeid(gg),:) &
                  / real(nTstepDay,dp)
          end do
          ! over gauges Hydro
          do gg = 1, nHydroGauges
             do intertt=iS, iE
             d_QMod(iDay, basin_gauge%gaugeid(gg)) = d_QMod(iDay, basin_gauge%gaugeid(gg)) &
                  + Hydro_DIS(intertt, basin_gauge%gaugeid(gg))
             end do
             d_QMod(iDay, basin_gauge%gaugeid(gg))=d_QMod(iDay, basin_gauge%gaugeid(gg)) &
                  / real(nTstepDay,dp)
          end do
    end do

    ! write in an ASCII file          ! OBS[nModeling_days X nGauges_total] , SIM[nModeling_days X nGauges_total] 
    call write_daily_obs_sim_discharge( basin_gauge%Q(:,:), d_Qmod(:,:) )
    call write_daily_obs_sim_conc( basin_gauge%GaugeConc(:,:,1), d_ConcMod(:,:,1) )


    ! free space
    deallocate(d_ConcMod)   
    deallocate(d_Qmod) 
	
  end subroutine wqm_write
 ! ------------------------------------------------------------------
  ! ------------------------------------------------------------------

  !     NAME
  !         write_daily_obs_sim_conc

  !     PURPOSE
  !>        \brief 

  !>        \details 
  !         
  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Xiaoqiang Yang - Modified from the corresponding nHM routines
  !>        \date Sep 2016 
 
  subroutine write_daily_obs_sim_discharge( Qobs, Qsim )

  use mo_errormeasures,          only: kge, nse, sse, rmse, bias
  use mo_utils,                  only: ge
  use mo_wqm_global_variables,   only: &
       dir_Out, basin_gauge, evalPer, nHydroGauges,hydro_eval_gaugenr,hydro_evalgaugeIDs

  use mo_file,   only: &
       file_daily_discharge, udaily_dis 
  use mo_julian,                 only: dec2date  
  use mo_message,                only: message
  use mo_string_utils,           only: num2str
  implicit none
  real(dp), dimension(:,:), intent(in)  :: Qobs
  real(dp), dimension(:,:), intent(in)  :: Qsim
  !local
    character(256) :: fName, formHeader, formData
    integer(i4) :: gg,pgg, tt, err
    integer(i4) :: day, month, year
    real(dp) :: newTime


    ! check the existance of file
    fName = trim(adjustl(dir_Out)) // trim(adjustl(file_daily_discharge))
       open(udaily_dis, file=trim(fName), status='unknown', action='write', iostat=err)
       if( err .ne. 0 ) then
          call message ('  IOError while opening ',trim(fName))
          call message ('  Error-Code ', num2str(err))
          stop
       end if
       ! header
       write(formHeader, *) '( 4a8, ' , hydro_eval_gaugenr,'(2X, a8, i7.7, 2X, a8, i7.7) )' 
       write(udaily_dis, formHeader) 'No', 'Day', 'Mon', 'Year', &
            ( 'Qobs_', basin_gauge%gaugeid(hydro_evalgaugeIDs(gg)), &
            'Qsim_', basin_gauge%gaugeid(hydro_evalgaugeIDs(gg)), gg=1, hydro_eval_gaugenr )

       ! form data
       write(formData, *) '( 4I8, ' , hydro_eval_gaugenr,'(2X,   f15.7 , 2X,  f15.7  ) )' 

       ! write data
       newTime  = real(evalPer%julStart,dp) - 0.5_dp

       do tt = 1, (evalPer%julEnd - evalPer%julStart + 1)          
          call dec2date(newTime, yy=year, mm=month, dd=day)
          !currently, only write out to **.out file 
          write(udaily_dis, formData) tt, day, month, year, ( Qobs(tt,hydro_evalgaugeIDs(gg)), Qsim(tt,hydro_evalgaugeIDs(gg)),&
                                                              gg=1, hydro_eval_gaugenr)
          newTime = newTime + 1.0_dp
       end do

       ! close file
       close(udaily_dis)
       ! ======================================================================
       ! screen output
       ! ======================================================================
       call message()
       call message('  OUTPUT: saved daily discharge file to ')
       call message(' ',trim(fname))
! statistical calculation is different from discharge, regarding to those missing value(nodata_dp)
! Here needs future work!!
       do gg=1, hydro_eval_gaugenr
          pgg = hydro_evalgaugeIDs(gg)
          if (count(ge(Qobs(:,pgg), 0.0_dp)) > 1 )  then
             call message('    KGE of daily discharge (gauge #',trim(adjustl(num2str(pgg))),'): ', &
                  trim(adjustl(num2str(kge(Qobs(:,pgg), Qsim(:,pgg), mask=(ge(Qobs(:,pgg), 0.0_dp)))))) )
             call message('    NSE of daily discharge (gauge #',trim(adjustl(num2str(pgg))),'): ', &
                  trim(adjustl(num2str(nse(Qobs(:,pgg), Qsim(:,pgg), mask=(ge(Qobs(:,pgg), 0.0_dp)))))) )
             call message('    PBIAS(%) of daily discharge (gauge #',trim(adjustl(num2str(pgg))),'): ', &
                  trim(adjustl(num2str(bias(Qobs(:,pgg), Qsim(:,pgg), mask=(ge(Qobs(:,pgg), 0.0_dp)))))) )
          end if
       end do

  end subroutine write_daily_obs_sim_discharge 
  ! ------------------------------------------------------------------

  !     NAME
  !         write_daily_obs_sim_conc

  !     PURPOSE
  !>        \brief 

  !>        \details 
  !         

  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Xiaoqiang Yang - Modified from the corresponding nHM routines
  !>        \date Sep 2016 
 
  subroutine write_daily_obs_sim_conc( Concobs, Concsim )

  use mo_errormeasures,          only: kge, nse, sse, rmse, bias
  use mo_utils,                  only: ge
  use mo_wqm_global_variables,   only: &
       dir_Out, basin_gauge, nWQGauges, evalPer,wqm_evalgaugeIDs,wqm_eval_gaugenr

  use mo_file,   only: &
       file_daily_conc, udaily_conc 
  use mo_julian,                 only: dec2date  
  use mo_message,                only: message
  use mo_string_utils,           only: num2str
  implicit none
  real(dp), dimension(:,:), intent(in)  :: Concobs
  real(dp), dimension(:,:), intent(in)  :: Concsim
  !local
    character(256) :: fName, formHeader, formData
    integer(i4) ::  gg,pgg, tt, err
    integer(i4) :: day, month, year
    real(dp) :: newTime

    ! check the existance of file
    fName = trim(adjustl(dir_Out)) // trim(adjustl(file_daily_conc))
       open(udaily_conc, file=trim(fName), status='unknown', action='write', iostat=err)
       if( err .ne. 0 ) then
          call message ('  IOError while opening ',trim(fName))
          call message ('  Error-Code ', num2str(err))
          stop
       end if

       ! header
       write(formHeader, *) '( 4a8, ' , wqm_eval_gaugenr,'(2X, a8, i7.7, 2X, a8, i7.7) )' 
       write(udaily_conc, formHeader) 'No', 'Day', 'Mon', 'Year', &
            ( 'INConcobs_', basin_gauge%gaugeid(wqm_evalgaugeIDs(gg)), &
            'INConcsim_', basin_gauge%gaugeid(wqm_evalgaugeIDs(gg)), gg=1, wqm_eval_gaugenr )

       ! form data
       write(formData, *) '( 4I8, ' , wqm_eval_gaugenr,'(2X,   f15.7 , 2X,  f15.7  ) )' 

       ! write data
       newTime  = real(evalPer%julStart,dp) - 0.5_dp

       do tt = 1, (evalPer%julEnd - evalPer%julStart + 1)          
          call dec2date(newTime, yy=year, mm=month, dd=day)
          !currently, only inorganic nitrogen(nitrate) is write out to **.out file 
          write(udaily_conc, formData) tt, day, month, year, ( Concobs(tt,wqm_evalgaugeIDs(gg)), Concsim(tt,wqm_evalgaugeIDs(gg)),&
                                                               gg=1, wqm_eval_gaugenr)
          newTime = newTime + 1.0_dp
       end do

       ! close file
       close(udaily_conc)
       ! ======================================================================
       ! screen output
       ! ======================================================================
       call message()
       call message('  OUTPUT: saved daily IN concentration file to ')
       call message(' ',trim(fname))
! statistical calculation is different from discharge, regarding to those missing value(nodata_dp)
! Here needs future work!!
       do gg=1, wqm_eval_gaugenr
          pgg = wqm_evalgaugeIDs(gg)
          if (count(ge(Concobs(:,pgg), 0.0_dp)) > 1 )  then
             call message('    KGE of daily nitrate (gauge #',trim(adjustl(num2str(pgg))),'): ', &
                  trim(adjustl(num2str(kge(Concobs(:,pgg), Concsim(:,pgg), mask=(ge(Concobs(:,pgg), 0.0_dp)))))) )
             call message('    NSE of daily nitrate (gauge #',trim(adjustl(num2str(pgg))),'): ', &
                  trim(adjustl(num2str(nse(Concobs(:,pgg), Concsim(:,pgg), mask=(ge(Concobs(:,pgg), 0.0_dp)))))) )
             call message('    PBIAS(%) of daily nitrate (gauge #',trim(adjustl(num2str(pgg))),'): ', &
                  trim(adjustl(num2str(bias(Concobs(:,pgg), Concsim(:,pgg), mask=(ge(Concobs(:,pgg), 0.0_dp)))))) )
          end if
       end do


  end subroutine write_daily_obs_sim_conc 
  
  
END MODULE mo_wqm_write
