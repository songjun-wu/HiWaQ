!> \file wqm_driver.f90
! --------------------------------------------------------------------------
!> \authors   Xiaoqiang Yang (UFZ & IGB)
!  CONTACT    
!> \version   
!> \date      Sep 2021

!  PURPOSE
!>            \brief 

! --------------------------------------------------------------------------

PROGRAM wqm_driver

  USE mo_finish,              ONLY : finish                         ! Finish with style
  USE mo_wqm_global_variables,    ONLY :                         &
       hydroModel,hydromodel_info,   & !some basic information of the coupled hydro model
       !wqm_param,       & !water quality parameters 
       dir_catchinfo,dir_gauges, simPer, nTstepDay, &     ! number of basins, frequency of input read
       dir_Out,dir_hydroInput


  USE mo_kind,                ONLY : i4                         ! number precision
  USE mo_message,             ONLY : message, message_text          ! For print out
  USE mo_string_utils,        ONLY : num2str             ! String magic
  USE mo_timer,               ONLY : timers_init, timer_start, timer_stop, timer_get! Timing of processes                     
  !Water quality model
  USE mo_wqm_readconfig,      ONLY : wqm_readconfig
  USE mo_wqm_readcatchdata,   ONLY : wqm_readcatchdata  
  USE mo_wqm_readobsdata,     ONLY : wqm_readobsdata 
  USE mo_wqm_readhydroinput,  ONLY : wqm_readhydroinput                                   
  USE mo_wqm_initialise,      ONLY : wqm_variables_initalloc
  USE mo_wqm_eval,            ONLY : wqm_eval
  USE mo_wqm_write,           ONLY : wqm_write
  !USE mo_wqm_restart,         ONLY : wqm_write_restart_files

  USE mo_optimization,        ONLY : optimization
  USE mo_opti_variables,      ONLY : optimize

  IMPLICIT NONE

  ! local
  integer(i4), dimension(8)             :: datetime         ! Date and time
              ! Counters
  integer(i4)                           :: iTimer           ! Current timer number
  integer(i4)                           :: nTimeSteps


  ! --------------------------------------------------------------------------
  ! START
  ! --------------------------------------------------------------------------
     call message('')
     call message('********************************************************************')
     call message('       Flexible Water Quality Model HiWaQ - Nitrate                 ')
     call message('               by Xiaoqiang Yang et al                              ')
     call message('  Department of Aquatic Ecosystem Analysis and Management (ASAM),   ')
     call message('      Helmholtz Centre for Environmental Research GmbH - UFZ        ')
     call message('         Department of Ecohydrology and Biogeochemistry,            ')
     call message(' Leibniz Institute of Freshwater Ecology and Inland Fisheries -IGB  ')
     call message('********************************************************************')
     call message('********************************************************************')



  !$OMP PARALLEL
  !$ n_threads = OMP_GET_NUM_THREADS()
  !$OMP END PARALLEL
  !$ call message('Run with OpenMP with ', trim(num2str(n_threads)), ' threads.')

  !--------------------------------
  !read configurations
  !--------------------------------
  call wqm_readconfig()
  call message('')
  call message('!! The coupled model is ', trim(hydroModel), '!!')
  call message('...developed by ', trim(hydromodel_info) )
  call message('')
  call date_and_time(values=datetime)
  message_text = trim(num2str(datetime(3),'(I2.2)'))//"."//trim(num2str(datetime(2),'(I2.2)')) &
       //"."//trim(num2str(datetime(1),'(I4.4)'))//" "//trim(num2str(datetime(5),'(I2.2)')) &
       //":"//trim(num2str(datetime(6),'(I2.2)'))//":"//trim(num2str(datetime(7),'(I2.2)'))
  call message('Start at ', trim(message_text), '.')
  call message('...Read configuration done!')
  !--------------------------------
  call message()
  call message('Input data:')

  call message('Catchment geophysical information:  ',   trim(dir_catchinfo) )
  call message('Gauging network and observations:   ',   trim(dir_gauges))
  call message('Hydrological inputs:                ',   trim(dir_hydroInput))
  call message('Output directory:                   ',   trim(dir_Out))
 
  call message('')


  ! Start timings
  call timers_init

  ! --------------------------------------------------------------------------
  ! READ AND (RE-)INITIALISE
  ! --------------------------------------------------------------------------
  itimer = 1
  call message()
  !--------------------------------
  !read basic catchment input data
  call message('  Read catchment information: ...')
  call timer_start(itimer)
  call wqm_readcatchdata()
  call timer_stop(itimer)
  call message('    in ', trim(num2str(timer_get(itimer),'(F9.3)')), ' seconds.')
  !--------------------------------

  call message('  Initialisation: ...')
  call timer_start(itimer)
  call wqm_variables_initalloc()
  call timer_stop(itimer)
  call message('    in ', trim(num2str(timer_get(itimer),'(F9.3)')), ' seconds.')

  !--------------------------------
  !read observation data
  call message('  Read gauging data: ...')
  call timer_start(itimer)
  call wqm_readobsdata()
  call timer_stop(itimer)
  call message('    in ', trim(num2str(timer_get(itimer),'(F9.3)')), ' seconds.')
  !--------------------------------
  !--------------------------------
  !read hydrological data
  call message('  Read hydrological inputs: ...')
  call timer_start(itimer)
  call wqm_readhydroinput()
  call timer_stop(itimer)
  call message('    in ', trim(num2str(timer_get(itimer),'(F9.3)')), ' seconds.')
  !--------------------------------

  !
  !call write_configfile()

  ! --------------------------------------------------------------------------
  ! RUN OR OPTIMIZE
  ! --------------------------------------------------------------------------
  iTimer = iTimer + 1
  call message()
  if (optimize) then
    call optimization()
  else
     ! --------------------------------------------------------------------------
     ! call **_eval
     ! --------------------------------------------------------------------------
    call message('  Run Water Quality Simulation')
    call timer_start(iTimer)
    call wqm_eval()
    call timer_stop(itimer)
    call message('    in ', trim(num2str(timer_get(itimer),'(F12.3)')), ' seconds.')
    ! ---------------------------------------------------------------------------  
    ! WRITE Discharge and Nutrients at gauging locations
    ! ---------------------------------------------------------------------------
    call wqm_write()

  end if

  
  ! --------------------------------------------------------------------------
  ! FINISH UP
  ! --------------------------------------------------------------------------
  itimer = itimer + 1
  ! call message()
  ! call message('  Write ouput data')
  ! call timer_start(itimer)
  ! ! call write_data()
  ! call timer_stop(itimer)
  ! call message('    in ', trim(num2str(timer_get(itimer),'(F9.3)')), ' seconds.')

  nTimeSteps = (simPer%julEnd - simPer%julStart + 1 ) * nTstepDay
  call date_and_time(values=datetime)
  call message()
  message_text = 'Done '//trim(num2str(nTimeSteps,'(I10)'))//" time steps."
  call message(trim(message_text))
  message_text = trim(num2str(datetime(3),'(I2.2)'))//"."//trim(num2str(datetime(2),'(I2.2)')) &
       //"."//trim(num2str(datetime(1),'(I4.4)'))//" "//trim(num2str(datetime(5),'(I2.2)')) &
       //":"//trim(num2str(datetime(6),'(I2.2)'))//":"//trim(num2str(datetime(7),'(I2.2)'))
  call message('Finished at ', trim(message_text), '.')
  call message()
  call finish(trim(hydroModel)//' + Water Quality','Finished!')

END PROGRAM wqm_driver

