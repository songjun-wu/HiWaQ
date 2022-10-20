!> \file mo_wqm_readobsdata.f90

!> \brief Reads water quality configure, input files and allocating WQ global variables.

!> \details This module is: firstly to read water qualtiy input data, including initial 
!>          pool values, cropdata, agricultural management and crop rotation classes 
!>          and spatial distribution; secondly to allocate global variables for water quality modeling.\n
!> 

!> \authors Xiaoqiang Yang
!> \date Sep 2016

MODULE mo_wqm_readobsdata


  USE mo_kind, ONLY: i4, sp, dp

  
  IMPLICIT NONE

  PUBLIC :: wqm_readobsdata              ! read observed water quality data (conc. at gauging stations)

  ! ------------------------------------------------------------------

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         wqm_readobsdata

  !     PURPOSE
  !>        \brief .

  !>        \details 

  !     CALLING SEQUENCE
  !         

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


  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Xiaoqiang Yang
  !>        \date Jul 2016

  subroutine wqm_readobsdata()
  

    use mo_wqm_global_variables,    only: &
         evalPer,simPer,nTstepDay,      &
         dir_gauges,basin_gauge,nHydroGauges, nWQGauges, nInflowGauges,   &
         nMeasPerDay, nEvalCmeasPerday,nAddinCmeasPerday, &
         numEvalCol, numAddinCol, &   
         Hydro_DIS,WQM_nutrient
    use mo_wqm_constants,           only: nodata_dp,maxnWQinputCols
    use mo_wqm_readtimeseries,      only: read_timeseries_conc,read_timeseries
    use mo_wqm_paste,               only: paste_conc
    use mo_append,                  only: paste
    use mo_string_utils,           only : num2str
    use mo_message,                only : message
    implicit none
  
    integer(i4)                  :: maxtimestep,i   !d1,d2,d3, 
    integer(i4)                  :: iGauge
    integer(i4), dimension(3)    :: start_tmp, end_tmp
    character(256)               :: fName
    real(dp), dimension(:), allocatable :: data_dp_1d
    logical, dimension(:), allocatable :: mask_1d
    real(dp), dimension(:,:), allocatable  :: data_dp_conc   !dim1= number of observations, dim2= "maxnWQinputCols"
    logical,  dimension(:,:), allocatable  :: mask_conc
	!file heading
!    integer(i4)                            :: numScol        ! number of wq data columns in gauge file
    character(256), dimension(:), allocatable :: head_str    ! name of each wq data column (heading string)
    character(256), dimension(:), allocatable :: Inhead_str    ! name of each wq data column (heading string)
    !allocate
    if (.not. allocated(numEvalCol))      allocate(numEvalCol(nWQGauges) ) 
    if (.not. allocated(nEvalCmeasPerday)) allocate(nEvalCmeasPerday(nWQGauges) )

    !*******************************************************************
    !initialise "Hydro_runoff"
    maxtimestep = (simPer%julEnd - simPer%julStart + 1 ) * nTstepDay
    allocate(Hydro_DIS(maxtimestep,nHydroGauges)) 
    Hydro_DIS = nodata_dp
    !*******************************************************************
    !allocate variable for simulated value (initialise "WQM_nutrient")
    allocate(WQM_nutrient(maxtimestep,nWQGauges,maxnWQinputCols)) 
    WQM_nutrient = nodata_dp

   !----------------------------------------------
    ! get start and end dates -- evalPer is sufficient
    start_tmp = (/evalPer%yStart, evalPer%mStart, evalPer%dStart/)
    end_tmp   = (/evalPer%yEnd,   evalPer%mEnd,   evalPer%dEnd  /)

   !for discharge observations	
    do iGauge = 1, nHydroGauges

       ! evaluation gauge
        fName = trim(adjustl(dir_gauges))//trim(adjustl(basin_gauge%fname_dis(iGauge)))
        call read_timeseries(trim(fName), 89, &
            start_tmp, end_tmp, &
            data_dp_1d, mask=mask_1d, nMeasPerDay=nMeasPerDay)
        data_dp_1d = merge(data_dp_1d, nodata_dp, mask_1d)
        call paste(basin_gauge%Q, data_dp_1d, nodata_dp )
        deallocate (data_dp_1d,mask_1d)
    end do
	
   !for WQ observations	
    do iGauge = 1, nWQGauges

       ! evaluation gauge
        fName = trim(adjustl(dir_gauges))//trim(adjustl(basin_gauge%fname_wqm(iGauge)))
        call read_timeseries_conc(trim(fName), 89, &
            start_tmp, end_tmp, maxnWQinputCols, numEvalCol(iGauge), head_str,&
            data_dp_conc, mask=mask_conc, nCmeasPerday=nEvalCmeasPerday(iGauge))
        data_dp_conc = merge(data_dp_conc, nodata_dp, mask_conc)

        call paste_conc(basin_gauge%GaugeConc, data_dp_conc, nodata_dp )

        deallocate (data_dp_conc,mask_conc)       
    end do

    !*******************************************************************
    ! additional inflow gauge
    ! upstream inflow and sewage plants input
    ! in mhm call InflowGauge%Q has to be initialized -- dummy allocation with period of basin 1 and initialization
    if (nInflowGauges .EQ. 0) then
        allocate( data_dp_1d(simPer%julEnd  - simPer%julStart + 1) ) 
        data_dp_1d = nodata_dp
        call paste(basin_gauge%inflowQ, data_dp_1d, nodata_dp)

        allocate( data_dp_conc(simPer%julEnd  - simPer%julStart + 1 , maxnWQinputCols)) 
        data_dp_conc = nodata_dp
        call paste_conc(basin_gauge%InflowGaugeConc, data_dp_conc, nodata_dp)
    else
       ! allocate
        if (.not. allocated(numAddinCol))   allocate(numAddinCol(nInflowGauges) )
        if (.not. allocated(nAddinCmeasPerday)) allocate(nAddinCmeasPerday(nInflowGauges) )
       !shoul cover all simPer (= Warmingday + evalPer)
       ! get start and end dates
       start_tmp = (/simPer%yStart, simPer%mStart, simPer%dStart/)
       end_tmp   = (/simPer%yEnd,   simPer%mEnd,   simPer%dEnd  /)

        do iGauge = 1, nInflowGauges
       !inflow discharge
          fName = trim(adjustl(dir_gauges))//trim(adjustl(basin_gauge%Inflowfname_dis(iGauge)))
          call read_timeseries(trim(fName), 89, &
               start_tmp, end_tmp, &
               data_dp_1d, mask=mask_1d, nMeasPerDay=nMeasPerDay)
          if ( .NOT. (all(mask_1d)) .and. &
          basin_gauge%inflow_ishead(iGauge)) then
             call message()
             call message('***ERROR: Nodata values in the HEAD inflow gauge time series. File: ', trim(fName))
             call message('          During simulation period from ', num2str(simPer%yStart) &
                  ,' to ', num2str(simPer%yEnd))
             stop
          end if
          data_dp_1d = merge(data_dp_1d, nodata_dp, mask_1d)
          call paste(basin_gauge%inflowQ, data_dp_1d, nodata_dp)
          deallocate (data_dp_1d,mask_1d)

       ! inflow gauge WQM
          fName = trim(adjustl(dir_gauges))//trim(adjustl(basin_gauge%Inflowfname_wqm(iGauge)))
          call read_timeseries_conc(trim(fName), 89, &
               start_tmp, end_tmp, maxnWQinputCols, numAddinCol(iGauge), Inhead_str,&
               data_dp_conc, mask=mask_conc, nCmeasPerday=nAddinCmeasPerday(iGauge) )
          data_dp_conc = merge(data_dp_conc, nodata_dp, mask_conc)

          !!post-process of inflow wq gauge data 
          !only two columns in inflow data file and these can only be IN and TN
          if (numAddinCol(iGauge).eq. 2) then   
             do i=1,size(data_dp_conc,1)
                if (mask_conc(i,1) .and. mask_conc(i,3)) then   !if measured IN and TN both have value
                   data_dp_conc(:,2) = data_dp_conc(:,3) - data_dp_conc(:,1)
                elseif (.not. mask_conc(i,1) .and. mask_conc(i,3)) then ! if only measured TN has value                      
                   data_dp_conc(i,1) = 0.8_dp * data_dp_conc(i,3)
                   data_dp_conc(i,2) = 0.2_dp * data_dp_conc(i,3)
                elseif (mask_conc(i,1) .and. (.not. mask_conc(i,3))) then !if only measured IN has value
                   data_dp_conc(i,2) = 0.25_dp * data_dp_conc(i,1)
                end if    
             end do
          end if
 
          !take care of a specific condition: only "TN" or "IN" was provided in inflow data
          if (numAddinCol(iGauge).eq. 1) then 
             if(any(mask_conc(:,3))) then
               data_dp_conc(:,1) = 0.8_dp * data_dp_conc(:,3)
               data_dp_conc(:,2) = 0.2_dp * data_dp_conc(:,3)
               !!NOTICE INFORMATION
               !print*, 'WARNING: Only total nitrogen data was provided at inflow gauge (', &
               !                basin_gauge%inflowgaugeid,')! '
               !print*, 'WARNING: 80% of TN was given to IN by default'
             else
               data_dp_conc(:,2) = 0.25_dp * data_dp_conc(:,1)
               !print*, 'WARNING: Only Inorganic nitrogen data was provided at inflow gauge (', &
               !                basin_gauge%inflowgaugeid,')! '
               !print*, 'WARNING: ON was given as 25% of IN by default'
             end if
          end if
          if ( any(data_dp_conc(:,1:2)-nodata_dp >= 0.0_dp) .and. &
          basin_gauge%inflow_ishead(iGauge)) then
             call message()
             call message('***ERROR: Nodata values in inflow gauge concentration time series. File: ', trim(fName))
             call message('          During simulation period from ', num2str(simPer%yStart) &
                  ,' to ', num2str(simPer%yEnd))
             stop
          end if
          call paste_conc(basin_gauge%InflowGaugeConc, data_dp_conc, nodata_dp)
          deallocate (data_dp_conc,mask_conc)
        end do
    end if

  end subroutine wqm_readobsdata
 
END MODULE mo_wqm_readobsdata
