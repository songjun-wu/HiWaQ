!> \author Xiaoqiang Yang
!> \date Sep 2021

MODULE mo_wqm_constants

  USE mo_kind, ONLY: i4, dp

  IMPLICIT NONE

  PRIVATE

  ! natural
  real(dp),    public, parameter :: H2Odens            = 1000.0_dp ! Density of water (kg/m3)

  ! computational
  integer(i4), public, parameter :: nodata_i4          = -9999_i4  ! [-]     global no data value
  real(dp),    public, parameter :: nodata_dp          = -9999._dp ! [-]     global no data value

  integer(i4), public, parameter :: maxngrps  = 50_i4  !the maximum subgroups for further re-group wqm parameter dependency
  integer(i4), public, parameter :: maxntypes = 50_i4  !the maximum types for each subgroup

  ! hydrologic modeling
  integer(i4), public, parameter :: maxnLandUse        = 50_i4     ! maximum number of allowed Land use types
  integer(i4), public, parameter :: maxnSoilType       = 50_i4     ! maximum number of allowed soil types
  integer(i4), public, parameter :: maxnRotation       = 50_i4     ! maximum number of allowed crop rotation sequences
  integer(i4), public, parameter :: maxnSoilLayer      = 10_i4     ! maximum number of allowed geological classes


  integer(i4), public, parameter :: maxnhydroGauges            = 50_i4     ! maximum number of hydrological gauges
  integer(i4), public, parameter :: maxnWQGauges               = 50_i4     ! maximum number of water quality gauges
  integer(i4), public, parameter :: maxnInflowGauges           = 50_i4     ! maximum number of inflow gauges
  integer(i4), public, parameter :: maxnProcesses              = 20_i4     ! maximum number of hydrological processes related to WQM dynamics
  integer(i4), public, parameter :: maxnSubSrfstore            = 10_i4     ! maximum number of subsurface conceptual storages and subsurface flow 
  integer(i4), public, parameter :: maxnChannelSeepageConponent = 10_i4    ! maximum number of flow components from terrestrial to stream channel
  integer(i4), public, parameter :: maxnNtransfors             = 20_i4     ! maximum number of N biogeochemical transformations  

  integer(i4), public, parameter :: maxnCropsInRotation = 10_i4     ! maximum number of crops in one rotation sequence

  !constant: maximum columns of wq data in gauging station
  ! currently, for nitrogen submodel only IN, ON, TN are available
  integer(i4), parameter, public  :: maxnWQinputCols = 3_i4
  ! number of substances involved, currently only IN and ON
  integer(i4), parameter, public  :: nsubstances = 2_i4 
  
  real(dp),    public, parameter :: P1_InitStateFluxes =    0.00_dp

END MODULE mo_wqm_constants

