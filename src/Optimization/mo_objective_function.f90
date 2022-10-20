!> \file mo_objective_function.f90

!> \authors Xiaoqiang Yang
!> \date Apr 2022
!> \Modified from mHM implementation
!>          All the objective functions are supposed to be minimized!                                      \n
!>               (1)  SO: Q:        1.0 - NSE                                                              \n
!>               (2)  SO: Q:        1.0 - lnNSE                                                            \n
!>               (3)  SO: Q:        1.0 - 0.5*(NSE+lnNSE)                                                  \n
!>               (4)  SO: Q:       -1.0 * loglikelihood with trend removed from absolute errors and
!>                                        then lag(1)-autocorrelation removed                              \n
!>               (5)  SO: Q:        ((1-NSE)**6+(1-lnNSE)**6)**(1/6)                                       \n
!>               (6)  SO: Q:        SSE                                                                    \n
!>               (7)  SO: Q:       -1.0 * loglikelihood with trend removed from absolute errors            \n
!>               (8)  SO: Q:       -1.0 * loglikelihood with trend removed from the relative errors and
!>                                        then lag(1)-autocorrelation removed                              \n
!>               (9)  SO: Q:        1.0 - KGE (Kling-Gupta efficiency measure)                             \n
!>               (14) SO: Q:        sum[((1.0-KGE_i)/ nGauges)**6]**(1/6) > combination of KGE of 
!>                                        every gauging station based on a power-6 norm                    \n
!>               (16) MO: Q:        1st objective: 1.0 - NSE                                               \n
!>                        Q:        2nd objective: 1.0 - lnNSE                                             \n
!>               (18) MO: Q:        1st objective: 1.0 - lnNSE(Q_highflow)  (95% percentile)               \n
!>                        Q:        2nd objective: 1.0 - lnNSE(Q_lowflow)   (5% of data range)             \n
!>               (19) MO: Q:        1st objective: 1.0 - lnNSE(Q_highflow)  (non-low flow)                 \n
!>                        Q:        2nd objective: 1.0 - lnNSE(Q_lowflow)   (5% of data range)eshold for Q \n
!>               (20) MO: Q:        1st objective: absolute difference in FDC's low-segment volume         \n
!>                        Q:        2nd objective: 1.0 - NSE of discharge of months DJF                    \n

!> \authors Juliane Mai
!> \date Dec 2012


MODULE mo_objective_function

  ! This module provides objective functions for optimization of the UFZ CHS mesoscale hydrologic model mHM.

  ! Written  Juliane Mai, Dec 2012
  ! Modified Stephan Thober, Oct 2015 - removed all none runoff objectives,
  !                                     these can be found mo_objective_functions_sm

  USE mo_kind, ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE

  !PUBLIC :: loglikelihood           ! loglikelihood with errormodel including linear trend and lag(1)-correlation
  !PUBLIC :: loglikelihood_stddev    ! loglikelihood where error is computed from difference of obs vs model
  PUBLIC :: objective ! single-objective function wrapper
  !PUBLIC :: multi_objective_runoff  ! multi-objective function wrapper
  PUBLIC :: extract_nitrate          ! extract nitrate time series 
  PUBLIC :: opti_eval

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !      NAME
  !          objective

  !>        \brief Wrapper for objective functions optimizing agains runoff.

  !>        \details The functions selects the objective function case defined in a namelist, 
  !>        i.e. the global variable \e opti\_function.\n
  !>        It return the objective function value for a specific parameter set.

  !     INTENT(IN)
  !>        \param[in] "real(dp) :: parameterset(:)"        1D-array with parameters the model is run with

  !     HISTORY
  !>        \author Juliane Mai
  !>        \date Dec 2012

  FUNCTION objective(parameterset)

    USE mo_opti_variables, ONLY: opti_function

    IMPLICIT NONE

    REAL(dp), DIMENSION(:), INTENT(IN)  :: parameterset
    REAL(dp)                            :: objective

    !write(*,*) 'parameterset: ',parameterset(:)
    select case (opti_function)
    ! case (1)
       ! ! 1.0-nse
       ! objective = objective_nse(parameterset)
    ! case (2)
       ! ! 1.0-lnnse
       ! objective = objective_lnnse(parameterset)
    ! case (3)
       ! ! 1.0-0.5*(nse+lnnse)
       ! objective = objective_equal_nse_lnnse(parameterset)
    ! case (4)
       ! ! -loglikelihood with trend removed from absolute errors and then lag(1)-autocorrelation removed
       ! objective = - loglikelihood_stddev(parameterset, 1.0_dp)
    case (5)
       ! ((1-NSE)**6+(1-lnNSE)**6)**(1/6)
       ! This objective function is modified for calibrating both discharge and nitrate --yangx 2017
       objective = objective_power6_nse_lnnse(parameterset)
    ! case (6)
       ! ! SSE
       ! objective = objective_sse(parameterset)
    ! case (7)
       ! ! -loglikelihood with trend removed from absolute errors
       ! objective = -loglikelihood_trend_no_autocorr(parameterset, 1.0_dp)
    ! case (8)
       ! ! -loglikelihood of approach 2 of Evin et al. (2013),
       ! !  i.e. linear error model with lag(1)-autocorrelation on relative errors
       ! objective = -loglikelihood_evin2013_2(parameterset)
    ! case (9)
       ! ! KGE
       ! objective = objective_kge(parameterset)
    ! case (14)
       ! ! combination of KGE of every gauging station based on a power-6 norm \n
       ! ! sum[((1.0-KGE_i)/ nGauges)**6]**(1/6) 
       ! objective = objective_multiple_gauges_kge_power6(parameterset)
    case default
       stop "Error objective: This opti_function is either not implemented yet or is not a single-objective one."
    end select

  END FUNCTION objective
  ! ------------------------------------------------------------------

  !      NAME
  !          objective_power6_nse_lnnse
  !     INTENT(IN)
  !>        \param[in] "real(dp) :: parameterset(:)"        1D-array with parameters the model is run with

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
  !     HISTORY
  !>        \author Juliane Mai and Matthias Cuntz
  !>        \date March 2014
  !         Modified, Stephan Thober, Jan 2015 - introduced extract_runoff
  !                   Stephan Thober, Aug 2015 - substituted nGaugesTotal variable with size(runoff, dim=2)
  !                                              to not interfere with mRM
  !                   Xiaoqiang Yang, Sep 2016 - Added objective function for nitrate calibration

  FUNCTION objective_power6_nse_lnnse(parameterset)

    use mo_errormeasures,    only: nse, lnnse
    use mo_wqm_global_variables, only: wqm_evalgaugeIDs!,hydro_optigaugeIDs

    implicit none

    real(dp), dimension(:), intent(in) :: parameterset
    real(dp)                           :: objective_power6_nse_lnnse

    ! local
    real(dp), allocatable, dimension(:,:) :: runoff                   ! modelled runoff for a given parameter set
    !                                                                 ! dim1=nTimeSteps, dim2=nGauges
    real(dp), allocatable, dimension(:,:,:) :: nutrient               ! modelled nutrient conc. 
	                                                                  ! dim1= nTimeSteps, dim2=nGauges, dim3= nsubstances
    integer(i4)                           :: gg                       ! gauges counter
    integer(i4)                           :: nGaugesTotal
    real(dp), dimension(:),   allocatable :: runoff_agg               ! aggregated simulated runoff
    real(dp), dimension(:),   allocatable :: runoff_obs               ! measured runoff
    logical,  dimension(:),   allocatable :: runoff_obs_mask          ! mask for measured runoff
    !for nitrate added by yangx 2016-09
    real(dp), dimension(:),   allocatable :: nitrate_agg               ! aggregated simulated nitrate
    real(dp), dimension(:),   allocatable :: nitrate_obs               ! measured nitrate
    logical,  dimension(:),   allocatable :: nitrate_obs_mask          ! mask for measured nitrate	
    real(dp)                              :: qcal,incal
	
    real(dp), parameter :: onesixth = 1.0_dp/6.0_dp

    call opti_eval(parameterset, runoff=runoff, nutrient = nutrient)
    nGaugesTotal = size(wqm_evalgaugeIDs)

    objective_power6_nse_lnnse = 0.0_dp
    do gg=1, nGaugesTotal
       ! extract runoff
       !call extract_runoff( gg, runoff, runoff_agg, runoff_obs, runoff_obs_mask )
	   !added by yangx 2016-09
       call extract_nitrate(wqm_evalgaugeIDs(gg), nutrient(:,:,1), nitrate_agg, nitrate_obs, nitrate_obs_mask)
       ! NSE + lnNSE
       !for discharge
       !qcal=( (1.0_dp-nse(  runoff_obs, runoff_agg, mask=runoff_obs_mask) )**6 + &
       !     (1.0_dp-lnnse(runoff_obs, runoff_agg, mask=runoff_obs_mask) )**6 )**onesixth
       !objective_power6_nse_lnnse = objective_power6_nse_lnnse + qcal
       !for nitrate concentration
       incal = ( (1.0_dp-nse(  nitrate_obs, nitrate_agg, mask=nitrate_obs_mask) )**6 + &
            (1.0_dp-lnnse(nitrate_obs, nitrate_agg, mask=nitrate_obs_mask) )**6 )**onesixth
       objective_power6_nse_lnnse = objective_power6_nse_lnnse + incal
       !aggregate with average value
       !objective_power6_nse_lnnse = objective_power6_nse_lnnse + (0.1_dp * incal+ 0.9_dp * qcal)

    end do
    ! objective function value which will be minimized
    objective_power6_nse_lnnse = objective_power6_nse_lnnse / real(nGaugesTotal,dp)

    write(*,*) 'objective_power6_nse_lnnse = ', objective_power6_nse_lnnse
    ! pause

    !deallocate( runoff_agg, runoff_obs, runoff_obs_mask )
    deallocate( nitrate_agg, nitrate_obs, nitrate_obs_mask)

  END FUNCTION objective_power6_nse_lnnse

  ! ------------------------------------------------------------------

  ! NAME
  !         extract_nitrate


  !     HISTORY
  !>        \author Xiaoqiang Yang 
  !>        \date Sep 2016

  ! ------------------------------------------------------------------
  subroutine extract_nitrate( gaugeId, nitrate, nitrate_agg, nitrate_obs, nitrate_obs_mask )

    use mo_wqm_global_variables, only: nEvalCmeasPerday, basin_gauge,nTstepDay,evalPer, &
        warmingDays
    use mo_message,              only: message
    use mo_utils,                only: ge


    implicit none

    ! input variables
    integer(i4),               intent(in) :: gaugeId      ! current gauge Id
    real(dp),  dimension(:,:), intent(in) :: nitrate       ! simulated nitrate

    ! output variables
    real(dp), dimension(:), allocatable, intent(out) :: nitrate_agg      ! aggregated simulated
    ! nitrate to the resolution
    ! of the measurement
    real(dp), dimension(:), allocatable, intent(out) :: nitrate_obs      ! extracted measured 
    ! nitrate to exactly the
    ! evaluation period
    logical,  dimension(:), allocatable, intent(out) :: nitrate_obs_mask ! mask of no data values
    ! in nitrate_obs

    ! local variables
    integer(i4)                         :: tt      ! timestep counter
    integer(i4)                         :: length  ! length of extracted time series
    integer(i4)                         :: factor  ! between simulated and measured time scale
    integer(i4)                         :: TPD_sim ! simulated Timesteps per Day
    integer(i4)                         :: TPD_obs ! observed Timesteps per Day
    real(dp), dimension(:), allocatable :: dummy

    ! copy time resolution to local variables
    TPD_sim = nTstepDay
    TPD_obs = nEvalCmeasPerday(gaugeId)

    ! check if modelled timestep is an integer multiple of measured timesteps
    if ( modulo( TPD_sim, TPD_obs) .eq. 0 ) then
       factor = TPD_sim / TPD_obs
    else
       call message(' Error: Number of modelled datapoints is no multiple of measured datapoints per day')
       stop
    end if
    ! get length of evaluation period times TPD_obs
    length = ( evalPer%julEnd - evalPer%julStart + 1 ) * TPD_obs
    ! extract measurements
    if ( allocated( nitrate_obs ) ) deallocate( nitrate_obs )
    allocate( nitrate_obs( length ) )
    nitrate_obs = basin_gauge%GaugeConc( 1 : length, gaugeId, 1 )

    ! create mask of observed nitrate
    if ( allocated( nitrate_obs_mask ) ) deallocate( nitrate_obs_mask )
    allocate( nitrate_obs_mask( length ) )
    nitrate_obs_mask = .false.
    forall(tt=1:length) nitrate_obs_mask(tt) = ge( nitrate_obs(tt), 0.000_dp)    

    ! extract and aggregate simulated nitrate
    if ( allocated( nitrate_agg ) ) deallocate( nitrate_agg )
    allocate( nitrate_agg( length ) )
    ! remove warming days
    length = ( evalPer%julEnd - evalPer%julStart + 1 ) * TPD_sim
    allocate( dummy( length ) )
    dummy = nitrate( warmingDays*TPD_sim+1:warmingDays*TPD_sim+length, gaugeId)
    ! aggregate nitrate
    length = ( evalPer%julEnd - evalPer%julStart + 1 ) * TPD_obs
    forall(tt=1:length) nitrate_agg(tt) = sum( dummy( (tt-1)*factor+1: tt*factor ) ) / &
         real(factor,dp)
    ! clean up
    deallocate( dummy )

  end subroutine extract_nitrate

  ! ==================================================================
  ! PRIVATE ROUTINES =================================================
  ! ==================================================================

  ! ------------------------------------------------------------------

  ! NAME
  !         opti_eval
  !     HISTORY
  !>        \author Xiaoqiang Yang
  !>        \date Apr 2022
  !>    Modified from
  !>        \author Stephan Thober
  !>        \date Oct 2015

  ! ------------------------------------------------------------------
  subroutine opti_eval(parameterset, runoff, nutrient)


    USE mo_wqm_eval,            ONLY : wqm_eval
    USE mo_wqm_global_variables,ONLY : wqm_param,nTotalParam

    implicit none

    ! input variables
    real(dp), intent(in) :: parameterset(:)

    ! output variables
    real(dp), allocatable, optional, intent(out) :: runoff(:,:)
    ! added by yangx 2ß16-09
    real(dp), allocatable, optional, intent(out) :: nutrient(:,:,:)

    !local
    integer(i4)   :: it,ip,igrp

    !update the new values from parameterset to the global paraemter: wqm_param
    it = 0_i4
    do ip=1,nTotalParam    
      do igrp = 1, wqm_param(ip)%ngroups
          it = it+1_i4
          wqm_param(ip)%pvalue(igrp) = parameterset(it)!parameter value
      end do
    end do

    if (present(runoff)) then
      call wqm_eval(runoff=runoff, nutrient= nutrient)
    else
      call wqm_eval()
    end if

  end subroutine opti_eval


  ! ! ------------------------------------------------------------------

  ! ! NAME
  ! !         extract_runoff

  ! !>        \brief extracts runoff data from global variables

  ! !>        \details extracts simulated and measured runoff from global variables,
  ! !>                 such that they overlay exactly. For measured runoff, only the runoff
  ! !>                 during the evaluation period are cut, not succeeding nodata values.
  ! !>                 For simulated runoff, warming days as well as succeeding nodata values
  ! !>                 are neglected and the simulated runoff is aggregated to the resolution
  ! !>                 of the observed runoff.\n

  ! !     INTENT(IN)
  ! !>        \param[in] "integer(i4) :: gaugeID"   - ID of the current gauge to process
  ! !>        \param[in] "real(dp)    :: runoff(:)" - simulated runoff at this gauge

  ! !     INTENT(INOUT)
  ! !         None

  ! !     INTENT(OUT)
  ! !>        \param[out] "real(dp)   :: runoff_agg(:)"      - aggregated simulated runoff at this gauge\n
  ! !>        \param[out] "real(dp)   :: runoff_obs(:)"      - extracted observed runoff\n
  ! !>        \param[out] "logical    :: runoff_obs_mask(:)" - masking non-negative values in runoff_obs\n

  ! !     INTENT(IN), OPTIONAL
  ! !         None

  ! !     INTENT(INOUT), OPTIONAL
  ! !         None

  ! !     INTENT(OUT), OPTIONAL
  ! !         None

  ! !     RETURN
  ! !         None

  ! !     RESTRICTIONS
  ! !         None

  ! !     EXAMPLE
  ! !         see use in this module above

  ! !     LITERATURE
  ! !         None

  ! !     HISTORY
  ! !>        \author Stephan Thober
  ! !>        \date Jan 2015

  ! ! ------------------------------------------------------------------
  ! subroutine extract_runoff( gaugeId, runoff, runoff_agg, runoff_obs, runoff_obs_mask )

    ! use mo_mrm_global_variables, only: gauge, nMeasPerDay, evalPer, warmingDays_mrm, nTstepDay
    ! use mo_message,              only: message
    ! use mo_utils,                only: ge

    ! implicit none

    ! ! input variables
    ! integer(i4),               intent(in) :: gaugeId      ! current gauge Id
    ! real(dp),  dimension(:,:), intent(in) :: runoff       ! simulated runoff

    ! ! output variables
    ! real(dp), dimension(:), allocatable, intent(out) :: runoff_agg      ! aggregated simulated
    ! ! runoff to the resolution
    ! ! of the measurement
    ! real(dp), dimension(:), allocatable, intent(out) :: runoff_obs      ! extracted measured 
    ! ! runoff to exactly the
    ! ! evaluation period
    ! logical,  dimension(:), allocatable, intent(out) :: runoff_obs_mask ! mask of no data values
    ! ! in runoff_obs

    ! ! local variables
    ! integer(i4)                         :: iBasin  ! basin id
    ! integer(i4)                         :: tt      ! timestep counter
    ! integer(i4)                         :: length  ! length of extracted time series
    ! integer(i4)                         :: factor  ! between simulated and measured time scale
    ! integer(i4)                         :: TPD_sim ! simulated Timesteps per Day
    ! integer(i4)                         :: TPD_obs ! observed Timesteps per Day
    ! real(dp), dimension(:), allocatable :: dummy

    ! ! copy time resolution to local variables
    ! TPD_sim = nTstepDay
    ! TPD_obs = nMeasPerDay

    ! ! check if modelled timestep is an integer multiple of measured timesteps
    ! if ( modulo( TPD_sim, TPD_obs) .eq. 0 ) then
       ! factor = TPD_sim / TPD_obs
    ! else
       ! call message(' Error: Number of modelled datapoints is no multiple of measured datapoints per day')
       ! stop
    ! end if

    ! ! extract basin Id from gauge Id
    ! iBasin = gauge%basinId( gaugeId )

    ! ! get length of evaluation period times TPD_obs
    ! length = ( evalPer( iBasin )%julEnd - evalPer( iBasin )%julStart + 1 ) * TPD_obs

    ! ! extract measurements
    ! if ( allocated( runoff_obs ) ) deallocate( runoff_obs )
    ! allocate( runoff_obs( length ) )
    ! runoff_obs = gauge%Q( 1 : length, gaugeId )

    ! ! create mask of observed runoff
    ! if ( allocated( runoff_obs_mask ) ) deallocate( runoff_obs_mask )
    ! allocate( runoff_obs_mask( length ) )
    ! runoff_obs_mask = .false.
    ! forall(tt=1:length) runoff_obs_mask(tt) = ge( runoff_obs(tt), 0.0_dp)    

    ! ! extract and aggregate simulated runoff
    ! if ( allocated( runoff_agg ) ) deallocate( runoff_agg )
    ! allocate( runoff_agg( length ) )
    ! ! remove warming days
    ! length = ( evalPer( iBasin )%julEnd - evalPer( iBasin )%julStart + 1 ) * TPD_sim
    ! allocate( dummy( length ) )
    ! dummy = runoff( warmingDays_mrm(iBasin)*TPD_sim + 1:warmingDays_mrm(iBasin)*TPD_sim + length, gaugeId )
    ! ! aggregate runoff
    ! length = ( evalPer( iBasin )%julEnd - evalPer( iBasin )%julStart + 1 ) * TPD_obs
    ! forall(tt=1:length) runoff_agg(tt) = sum( dummy( (tt-1)*factor+1: tt*factor ) ) / &
         ! real(factor,dp)
    ! ! clean up
    ! deallocate( dummy )

  ! end subroutine extract_runoff



END MODULE mo_objective_function
