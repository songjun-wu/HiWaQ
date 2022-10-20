!> \file mo_file.f90

!> \brief Provides file names and units for independent WQM-N module

!> \details 

!> \author 
!> \date 

MODULE mo_file

  IMPLICIT NONE
  
  !> Driver file
  CHARACTER(len=*), PARAMETER :: file_main              = 'wqm_driver.f90'              ! Driver
  !> Namelist file name
  CHARACTER(len=*), PARAMETER :: file_namelist_wqm      = 'wqm_config.nml'             ! Namelist
  !> Unit for namelist
  INTEGER,          PARAMETER :: unamelist              = 19                            ! 
  !> Parameter namelists file name
  CHARACTER(len=*), PARAMETER :: file_param_wqm         = 'wqm_parameter.nml'           ! Parameter namelists
  !> Unit for namelist
  INTEGER,          PARAMETER :: unamelist_param        = 89                            ! 
  !> file defining outputs
  CHARACTER(len=*), PARAMETER :: file_stateoutput_wqm    = 'wqm_outputs.nml'             ! output states and fluxes
  !> Unit for file defining mHM's outputs
  INTEGER,          PARAMETER :: ustatOutput             = 84  

  !> file name to store the final daily concentration results.
  character(len=*), parameter :: file_daily_discharge  = 'daily_discharge.out' 
  INTEGER,          PARAMETER :: udaily_dis             = 99  

  !> file name to store the final daily concentration results.
  character(len=*), parameter :: file_daily_conc  = 'daily_concentration.out' 
  INTEGER,          PARAMETER :: udaily_conc             = 25  

  !> Input nCols and nRows of binary meteo files are in header file
  CHARACTER(len=*), PARAMETER :: file_meteo_header      = 'header.txt'                  ! Meteo header
  !> Unit for meteo header file
  INTEGER,          PARAMETER :: umeteo_header          = 50                            ! 
  !> File ending of meteo files
  CHARACTER(len=*), PARAMETER :: file_meteo_binary_end  = '.bin'                        ! Meteo
  !> Unit for meteo files
  INTEGER,          PARAMETER :: umeteo                 = 51                            ! 
  !> Soil database file (iFlag_soilDB = 0) = classical mHM format
  CHARACTER(len=*), PARAMETER :: file_soil_database     = 'soil_classdefinition.txt'    ! Soil data base 
  !>> Soil database file (iFlag_soilDB = 1)
  CHARACTER(len=*), PARAMETER :: file_soil_database_1   = 'soil_classdefinition_iFlag_soilDB_1.txt' 
  !> Unit for soil data base
  INTEGER,          PARAMETER :: usoil_database         = 52                            ! 
  !> DEM input data file
  CHARACTER(len=*), PARAMETER :: file_dem               = 'dem.asc'                     ! DEM
  !> Unit for  DEM input data file
  INTEGER,          PARAMETER :: udem                   = 53                            ! 
  !> slope input data file
  CHARACTER(len=*), PARAMETER :: file_slope             = 'slope.asc'                   ! slope
  !> Unit for  slope input data file
  INTEGER,          PARAMETER :: uslope                 = 54                            ! 
  !> aspect input data file
  CHARACTER(len=*), PARAMETER :: file_aspect            = 'aspect.asc'                  ! aspect
  !> Unit for  aspect input data file
  INTEGER,          PARAMETER :: uaspect                = 55                            ! 
  !> flow accumulation input data file
  CHARACTER(len=*), PARAMETER :: file_facc              = 'facc.asc'                    ! flow accumulation
  !> Unit for  flow accumulation input data file
  INTEGER,          PARAMETER :: ufacc                  = 56                            ! 
  !> flow direction input data file
  CHARACTER(len=*), PARAMETER :: file_fdir              = 'fdir.asc'                    ! flow direction
  !> Unit for  flow direction input data file
  INTEGER,          PARAMETER :: ufdir                  = 57                            ! 
  !> hydrogeological classes input data file
  CHARACTER(len=*), PARAMETER :: file_hydrogeoclass     = 'geology_class.asc'           ! hydrogeological classes
  !> Unit for  hydrogeological classes input data file
  INTEGER,          PARAMETER :: uhydrogeoclass         = 58                            ! 
  !> soil classes input data file
  CHARACTER(len=*), PARAMETER :: file_soilclass         = 'soil_class.asc'              ! soil classes
  !> Unit for  soil classes input data file
  INTEGER,          PARAMETER :: usoilclass             = 59                            ! 
  !> LAI classes input data file
  CHARACTER(len=*), PARAMETER :: file_laiclass          = 'LAI_class.asc'               ! LAI classes
  !> Unit for  LAI input data file
  INTEGER,          PARAMETER :: ulaiclass              = 60                            ! 
  !> Unit for  LCover input data file
  ! Land cover files are explicit given by the user in mhm.nml file
  INTEGER,          PARAMETER :: ulcoverclass           = 61                            ! 

  !> geological formation lookup table file
  CHARACTER(len=*), PARAMETER :: file_geolut            = 'geology_classdefinition.txt' ! geolog. formation lookup table
  !> Unit for geological formation lookup table file
  INTEGER,          PARAMETER :: ugeolut                = 64                            ! 

  !> LAI classes lookup table file
  CHARACTER(len=*), PARAMETER :: file_lailut            = 'LAI_classdefinition.txt'     ! LAI classes lookup table
  !> Unit for LAI classes lookup table file
  INTEGER,          PARAMETER :: ulailut                = 65                            ! 

  !> unit for discharge time series 
  INTEGER,          PARAMETER :: udischarge             = 66                            ! 

  !> file containing mHM configuration
  CHARACTER(len=*), PARAMETER :: file_config            = 'ConfigFile.log'              ! configuration
  !> Unit for file containing mHM configuration
  INTEGER,          PARAMETER :: uconfig                = 68                            ! 

  !> file defining optimization outputs (objective and parameter set)
  CHARACTER(len=*), PARAMETER :: file_opti              = 'FinalParam.out'              ! final parameters & objective
  !> Unit for file optimization outputs (objective and parameter set)
  INTEGER,          PARAMETER :: uopti                  = 72                            ! 
  !> file defining optimization outputs in a namelist format (parameter set)
  CHARACTER(len=*), PARAMETER :: file_opti_nml          = 'FinalParam.nml'              ! final parameters
  !> Unit for file optimization outputs in a namelist format (parameter set)
  INTEGER,          PARAMETER :: uopti_nml              = 73                            ! 

  !> Input nCols and nRows of binary gridded LAI files are in header file
  CHARACTER(len=*), PARAMETER :: file_lai_header        = 'header.txt'                  ! LAI header
  !> Unit for LAI header file
  INTEGER,          PARAMETER :: ulai_header            = 75                            !
  !> Binary file ending of LAI files
  CHARACTER(len=*), PARAMETER :: file_lai_binary_end    = '.bin'                        ! Gridded LAI
  !> Unit for binary LAI files
  INTEGER,          PARAMETER :: ulai                   = 76                            ! 
  !> unit for tws time series 
  INTEGER,          PARAMETER :: utws                   = 77                            ! 

  
END MODULE mo_file
