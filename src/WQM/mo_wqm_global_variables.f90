!> \file mo_wqm_global_variables.f90

!> \brief Global variables used in reading, writing and startup for water quality modelling.

!> \details

!> \authors Xiaoqiang Yang
!> \date Jun 2016
!>  Modified   X Yang, May 2018   added global variables for riparian zone shading effect


MODULE mo_wqm_global_variables

  
  USE mo_kind,             ONLY: i4, i8, dp
  USE mo_wqm_constants,    ONLY: maxngrps, maxntypes,nsubstances,maxnProcesses, &
                                 maxnSubSrfstore,maxnChannelSeepageConponent


  IMPLICIT NONE
  
  public :: period
  public :: HRUproperty
  public :: wqm_parameterization  
  public :: ReservoirStore
  public :: rotationType
  public :: cropType
  public :: wqm_gaugeinfo
  public :: gridGeoRef
  public :: ncvars

  !********************************************************
  ! READ CONFIGURATION
  !********************************************************

  !main configuration from wqm_config.nmml
    character(256),   public :: hydroModel      ! the name of the coupled hydrologocal model
    character(256),   public :: hydromodel_info  ! brief information of developer and original literature
    integer(i4), public :: timeStep ! [h] simulation time step (= TS) in [h]  
    integer(i4), public :: nTstepDay ! number of steps per day = 24/timestep

    real(dp),   public :: resolution ! [m ]  modeling resolution
    real(dp),   public :: rout_resolution ! [m ] routing resolution  
    integer(i4),   public :: model_type ! 1 for grid-based, 2 for semi-distributed

    integer(i4), public :: warmingDays ! number of days for warm up period

    type period
      integer(i4) :: dStart      ! first day
      integer(i4) :: mStart      ! first month
      integer(i4) :: yStart      ! first year
      integer(i4) :: dEnd        ! last  day
      integer(i4) :: mEnd        ! last  month
      integer(i4) :: yEnd        ! last  year
      integer(i4) :: julStart    ! first julian day 
      integer(i4) :: julEnd      ! last  julian day 
      integer(i4) :: nObs        ! total number of observations
    end type period


    type(period),   public :: warmPer     ! time period for warming
    type(period),   public :: evalPer     ! time period for model evaluation
    type(period),   public :: simPer      ! warmPer + evalPer
 
 !defined HRUs types based on the soil and land use classes
    type HRUproperty
      integer(i4), dimension(:), allocatable    :: hruids           ! HRU id
      integer(i4), dimension(:), allocatable    :: landuseids    ! Landuse id
      integer(i4), dimension(:), allocatable    :: soilids       ! soiltype id
      real(dp),    dimension(:,:), allocatable  :: soildepth     ! soil depth for each layer (dim1 =nHRUs, dim2 = nsoillayers)
      real(dp),    dimension(:,:), allocatable  :: wiltingpoint  ! wilting point
      real(dp),    dimension(:,:), allocatable  :: saturatedSM   ! saturation soil water
    end type HRUproperty
    type (HRUproperty), public    :: geoProperty 

    type ncvars
      integer(i4)         :: varids      ! ID
      character (256) :: var_NAME      ! NAME
      character (256) :: var_UNIT      ! UNIT
    end type ncvars
	
    type(ncvars), dimension(:), allocatable, public :: ter_variables,chan_variables     ! time period for warming


  !directories for WQM and coupling
    character(256),   public :: dir_catchinfo      ! Directory where catchmenet information files are located
    character(256),   public :: dir_gauges         ! Directory where observations at gauges are located
    character(256),   public :: dir_hydroInput     ! Directory where hydrological inputs are located
    character(256),   public :: dir_Out            ! Directory where outputs should be written
  ! file names of catchment information
    character(256),   public :: fn_catchDomain     !catchment domain (basemap for grid, geodata for semi-dis)
    character(256),   public :: fn_flowDir         !flow direction for grid/subcatch connection
    character(256),   public :: fn_spatialorganization !alternatively given directly the spatial organization information
    character(256),   public :: fn_channel_mask
    character(256),   public :: fn_channel_length
    character(256),   public ::fn_channel_width
    character(256),   public ::fn_cropinfo        !all crop related information "cropdata.txt"
    character(256),   public ::fn_lut_rotation    !crop rotation's look-up-table
    character(256),   public ::fn_lut_monLAI
    character(256),   public ::fn_initial_condition
    character(256),   public ::fn_init_concIN
    character(256),   public ::fn_init_concON
    character(256),   public ::fn_globalRadiation
    !dim = number of substances 
    real(dp), dimension(nsubstances), public :: const_crain     !constant given IN and ON concentration (mg/l) in precipitation
    real(dp), dimension(nsubstances), public :: const_csnowmelt  !constant given IN and ON concentration (mg/l) in snowmelt
  !	geophysical related information
    integer(i4), public :: nLanduseTotal
    integer(i4), public :: nLUperGrid
    integer(i4), public :: nSoilTotal
    integer(i4), public :: nSoilperGrid
    integer(i4), public :: nRotationTotal
    integer(i4), public :: nRotationperGrid
    integer(i4), public :: nSoilLayers
  ! file names of related geo data
    character(256), public                                  ::fn_LanduseTypes !land use id map
    character(256), public                                  ::fn_SoilTypes    !soil type id map
    character(256), public                                  ::fn_RotationTypes!crop rotation id map
    !for distributed 
    character(256), dimension(:), allocatable, public       :: fn_flanduse !dim = number of landuse types
    character(256), dimension(:), allocatable, public       :: fn_fsoil    !dim = number of soil types
    character(256), dimension(:), allocatable, public       :: fn_frotation !dim = number of rotation types
    character(256), public                                  :: fn_soil_totalDepth
    character(256), dimension(:), allocatable, public       :: fn_soilDepth !lower depth - upper depth for each layer
    character(256), dimension(:), allocatable, public       :: fn_wiltpoint !filename soil parameter wilting point
    character(256), dimension(:), allocatable, public       :: fn_satMoist  ! filename soil parameter saturation soil moisture
    character(256),                            public       :: fn_classProperties!filename properties for each HRUs (landuse, soil type and soil properties)fn_classProperties
  ! gauging networks and number gauges
    character(256),  public :: fn_gauge_loc  !file name of gauge ids map
    integer(i4),     public :: nHydroGauges
    integer(i4),     public :: nWQGauges
    integer(i4),     public :: nInflowGauges
    integer(i4),     public :: fdir_type
  ! process related variable names
	!process
    integer(i4),     dimension(maxnProcesses), public       :: processCase
    integer(i4),                                     public :: nEvapLayer      !number of soil layer for soil evaporation
    integer(i4),                                     public :: nTranspLayer    !number of soil layer for traspiration
    integer(i4),                                     public :: nSubsurfaceflow !number of subsurface flow components (e.g.,several interflows)
    character(256),                                  public :: vn_precipitation !variable name of flux direct rainfall into surface ponding storage
    character(256),                                  public :: vn_airTemperature   !variable name of air temperature
    character(256),                                  public :: vn_WS_interception  !variable name of storage canopy interception 
    character(256),                                  public :: vn_evapCanopy       !variable name of evaporation from canopy
    character(256),                                  public :: vn_throughfall      !variable name of flux throughfall that goes into ponding storage
    character(256),                                  public :: vn_snowfall         !variable name of flux snowfall that goes into snowpack (water equivalent)
    character(256),                                  public :: vn_WS_snowpack      !variable name of storage snow package 
    character(256),                                  public :: vn_snowmelt         !variable name of flux snowmelt that goes into ponding storage
    character(256),                                  public :: vn_infiltration     !variable name of flux that infiltrated from surface to first soil layer 
    character(256),                                  public :: vn_WS_ponding       !variable name of storage surface ponding water 
    real(dp),                                   public :: ETfluxUT    !unit convert constant for ET fluxes to m/timestep 
    character(256),  dimension(:), allocatable, public :: vn_ETevap   !variable name of fluxes soil evaporation (dim = nEvapLayer)
    character(256),  dimension(:), allocatable, public :: vn_ETtransp !variable name of fluxes transpiration (dim = nTranspLayer)
    character(256),                             public :: vn_reinfiltration   !variable name of surface reinfiltration
    character(256),  dimension(:), allocatable, public :: vn_reinfil_soil     !variable name of reinfiltration to each layer
    character(256),                             public :: vn_surfaceflow      !variable name of flux surface runoff
    character(256),                             public :: vn_WS_surfaceflow   !variable name of storage surface ponding, where surface runoff is generated from
    character(256),                             public :: vn_fromWS_surfaceflow  !variable name of storage from where baseflow storage gets replenished
    real(dp),                                   public :: SrfOVFmixratio      !the weight of top-layer soil water concentration,
                                                                              !if surface flow and top soil water mixing is considered during the terrestrial routing	
    character(256),                             public :: SrfOVFmixratio_vn   !or the ratio information is supplied from the coupled hydro model, which has the higher priority!!
                                                                              !if not specified, the default ratio of soil water is max(0.5,theta_1/satSM_1)	
    character(256),  dimension(:), allocatable, public :: vn_subsurfaceflow    !variable name of fluxes subsurface runoff (dim = nSubsurfaceflow)
    character(256),  dimension(:), allocatable, public :: vn_WS_subsurfaceflow !variable name of storages where subsruface flows are respectively generated from
    character(256),  dimension(:), allocatable, public :: vn_fromWS_subsurfaceflow !variable name of storages from where subsruface storage get replenished
    character(256),                             public :: vn_baseflow         !variable name of flux baseflow runoff
    character(256),                             public :: vn_WS_baseflow      !variable name of storage where baseflow is generated from
    character(256),                             public :: vn_fromWS_baseflow  !variable name of storage from where baseflow storage gets replenished
    character(256),                             public :: vn_directflow         !variable name of flux direct runoff from sealed subareas
    character(256),                             public :: vn_Esealed           !variable name of flux evaporation from sealed subareas
    character(256),                             public :: vn_WS_directflow      !variable name of the sealed storage where direct runoff is generated from
    integer(i4), dimension(:), allocatable,     public :: Sealed_LanduseIDs   !for the areal fraction of the sealed storage (according to specific landuse types)
    integer(i4),                                public :: Sealed_LanduseNrs   !number of landuse types belonging to "impermeable"
  ! structure of hydro conceptualization
    character(256),                             public :: inputFormat         !format of hydrological input file (i.e., "nc" or "bin") the "nc" file is preferred
    logical,                                    public :: inputcheck          !check if inputs cover all model domain, default == .FALSE.
    character(256),                             public :: inputhydro_filename !file name of the hydro information input
    real(dp),                                   public :: inputfluxUT         !flux unit transfer constant. The flux unit of this model is implemented as m/timestep.
    character(256),                             public :: inputFormat_rout    !format of hydrological routing input file (i.e., "nc" or "bin") the "nc" file is preferred
    character(256),                             public :: inputhydro_rout_filename !!file name of the hydro routing information input
    real(dp),                                   public :: input_rout_fluxUT   !channel store and stream flow unit transfer constant, to m3 and m3/s, respectively.
    character(256),                             public :: vn_soilTemperature  !optional external input of soil temperature
    !variable names of fluxes related to each water storage
    !vertically FOUR: in the order of "up_in","down_out","up_out","down_in"
    !laterally TWO: in the order of "lateral in", "lateral out"
    character(256),  dimension(4),                public :: vn_WS_interception_ver!interception
    character(256),  dimension(4),                public :: vn_WS_snowpack_ver!snowpack
    character(256),  dimension(4),                public :: vn_WS_ponding_ver!surface ponding
    character(256),  dimension(4),                public :: vn_WS_surfaceflow_ver      !surfaceflow ponding
    character(256),  dimension(:,:), allocatable, public :: vn_WS_subsurfaceflow_ver   !subsurface storages (dim1=nSubsurfaceflow, dim2 = 4)
    character(256),  dimension(4),                public :: vn_WS_baseflow_ver         !deep stroage where baseflow is generated
    character(256),  dimension(2),                public :: vn_WS_interception_lat     !interception
    character(256),  dimension(2),                public :: vn_WS_snowpack_lat         !snowpack
    character(256),  dimension(2),                public :: vn_WS_ponding_lat          !surface ponding
    character(256),  dimension(2),                public :: vn_WS_surfaceflow_lat      !surfaceflow storage
    character(256),  dimension(:,:), allocatable, public :: vn_WS_subsurfaceflow_lat   !subsurface storages (dim1=nSubsurfaceflow, dim2 = 2)
    character(256),  dimension(2),                public :: vn_WS_baseflow_lat         !deep stroage where baseflow is generated
    character(256),  dimension(:),   allocatable, public :: vn_WS_soil         !soil storage for each soil layer (dim = nSoilLayers)
    character(256),  dimension(:),   allocatable, public :: vn_SatSM           !saturated soil moisture content --soil porosity if supplied by the couple hydro model
    character(256),  dimension(:),   allocatable, public :: vn_WP              !wilting point --if supplied by the coupled hydro model
    character(256),  dimension(:,:), allocatable, public :: vn_WS_soil_ver     !vertical flux exchanges (dim1 = nSoilLayers, dim2 = 4)
    character(256),  dimension(:,:), allocatable, public :: vn_WS_soil_lat     !lateral flux exchanges (dim1 = nSoilLayers, dim2 =2)

    character(256),  dimension(:),   allocatable, public :: all_vn_WS !store all vn_WS_** names

  !channle routing variables
    character(256),                               public :: vn_chan_store    !variable name of storage channel storage
    character(256),                               public :: vn_chan_in       !variable name of flux channl inflow from upstream 
    character(256),                               public :: vn_chan_out      !variable name of flux channel outflow to downstream
    character(256),                               public :: vn_srftochan     !variable name of flux surface seepage to channel
    character(256),                               public :: vn_fromWS_srftochan   !variable name of storage from where srf seepage is generated
    character(256),                               public :: vn_baseftochan   !variable name of flux baseflow seepage to channel
    character(256),                               public :: vn_fromWS_baseftochan !variable name of storage from where baseflow seepage is generated
    integer(i4),                                  public :: nSubsrfToChannel !number of subsurface runoff components from terrestrial to channel
    character(256),  dimension(:),   allocatable, public :: vn_subsrftochan !variable of each runoff component (dim = nSubsrfToChannel)
    character(256),  dimension(:),   allocatable, public :: vn_fromWS_subsrftochan !variable name of storage from where each subsurface flow is generated

    ! number of substances involved, currently only IN and ON
    !integer(i4), parameter, public                    :: nsubstances = 2   !moved to mo_wqm_constants  
    integer(i4), public                               :: nCroptation         ! Number of crop rotation types 
    integer(i4), public                               :: num_crops           ! number of total crop types
    !crop data and land-use fraction

    real(dp), dimension(:,:), allocatable, public     :: L1_frotation        ! area fraction at Level L1
    real(dp), dimension(:,:), allocatable, public     :: L1_fLAI             ! area fraction at Level L1
    real(dp), dimension(:,:), allocatable, public     :: L1_fsoil            ! area fraction at Level L1
    real(dp), dimension(:),   allocatable, public     :: L1_fSealed          ! area fraction at Level L1
    !soil depth for each layer
    real(dp), dimension(:,:), allocatable, public     :: L1_soildepth ! dim1=ncells, dim2 = nlayers
    real(dp), dimension(:,:), allocatable, public     :: L1_wiltingpoint
    real(dp), dimension(:,:), allocatable, public     :: L1_saturatedSM

    !(for crop rotation calculation)
    integer(i4), public                               :: year_start          ! simulation starting year
    integer(i4), public                               :: day_preyr           ! number of days in the previous year.

    type wqm_gaugeinfo
     !dim = number of total stations
      integer(i4), dimension(:), allocatable     :: gaugeid !the id for each station
      integer(i4), dimension(:), allocatable     :: gaugeid_loc !the location within the modeling domain     	  
      integer(i4), dimension(:), allocatable     :: wqm_gaugeid
      integer(i4), dimension(:), allocatable     :: wqm_gaugeid_loc !
      integer(i4), dimension(:), allocatable     :: inflowgaugeid
      logical,     dimension(:), allocatable     :: inflow_ishead   !if model simulations should be replaced by inflow
      integer(i4), dimension(:), allocatable     :: inflowgaugeid_loc
      character(256), dimension(:), allocatable  :: fname_dis        ! file name of measured dis data at evaluation gauging station 
     ! file name of additional inflow wq data(upstream and point source)   
      character(256), dimension(:), allocatable  :: Inflowfname_dis 
     !dim = number of total stations
      character(256), dimension(:), allocatable  :: fname_wqm        ! file name of measured wq data at evaluation gauging station
     ! file name of additional inflow wq data(upstream and point source)   
      character(256), dimension(:), allocatable  :: Inflowfname_wqm  
      real(dp), dimension(:,:), allocatable      :: Q ! [m3 s-1] observed daily mean discharge (simPer)
     !                                                 ! dim1=number observations, dim2=number of gauges
      real(dp), dimension(:,:), allocatable      :: inflowQ ! [m3 s-1] inflow daily mean discharge (simPer)
     !                                                 ! dim1=number observations, dim2=number of gauges
     !dim1= number of observation, dim2=number of total gauges, dim3=number of measured data type
      real(dp), dimension(:,:,:), allocatable    :: GaugeConc        ! measured wq data at evaluation gauging station
      real(dp), dimension(:,:,:), allocatable    :: InflowGaugeConc  ! measured wq data at additional inflow station
    end type wqm_gaugeinfo  
    type(wqm_gaugeinfo), public                            :: basin_gauge        ! basin info for wqm
    ! station information included for optimization/evaluation
    integer(i4), public                               :: hydro_eval_gaugenr
    integer(i4), dimension(:), allocatable,public     :: hydro_evalgaugeIDs
    integer(i4), public                               :: wqm_eval_gaugenr
    integer(i4), dimension(:), allocatable,public     :: wqm_evalgaugeIDs

    !organized parameterization
    type wqm_parameterization
      character(256)                                 :: param_name
      integer(i4)                                    :: category
      integer(i4)                                    :: ngroups
      integer(i4),  dimension(maxngrps,maxntypes)    :: groupids
      real(dp),     dimension(maxngrps)              :: lowbound
      real(dp),     dimension(maxngrps)              :: upbound 
      real(dp),     dimension(maxngrps)              :: pvalue
      integer(i4),  dimension(maxngrps)              :: flag
    end type wqm_parameterization

    !define the global variables 
    type(wqm_parameterization),  dimension(:), allocatable,   public :: wqm_param     ! parameteriazation of wqm
    integer(i4),                                             public :: nTotalParam   !number of total parameters,the dimension of type wqm_param

  !********************************************************
  ! READ DATA
  !********************************************************
  !all L1_$$ related variables should be initially allocated and initilized
   !read catch data  
    type gridGeoRef
      integer(i4)    :: ncols        ! Number of columns
      integer(i4)    :: nrows        ! Number of rows
      real(dp)       :: xllcorner    ! x coordinate of the lowerleft corner
      real(dp)       :: yllcorner    ! y coordinate of the lowerleft corner
      real(dp)       :: cellsize     ! Cellsize x = cellsize y
      real(dp)       :: nodata_value ! Code to define the mask
    end type gridGeoRef
    type (gridGeoRef), public    :: m_Level
    type (gridGeoRef), public    :: r_Level

    logical,        dimension(:,:),     allocatable,   public   :: Model_mask        ! FROM MODEL DOMAIN dim1=ncols, dim2=nrows
    logical,        dimension(:,:),     allocatable,   public   :: Rout_mask         ! from fdir maps dim1 = ncols, dim2 = nrows
    integer(i4),                                     public   :: L1_nCells    !number of cells within modeling domain
    integer(i4),                                     public   :: L1_nLinks    !number of links
    real(dp),       dimension(:),   allocatable,   public     :: L1_area     ! [m2] the areal share of each grid from model domain, multipled with (resolution **2)
                                                                             !considering boundary cells
    real(dp),       dimension(:),   allocatable,   public     :: R_areaCell ![m2] routing grid area if rout_resolution > resolution
    integer(i4),    dimension(:),   allocatable,   public     :: RIdonM   ! modeling grid ids at the routing resolution, for accummulation purpose
    integer(i4),    dimension(:),   allocatable,   public   :: L1_fDir !connection direction for spatial organization, dim = L1_nCells
    integer(i4),    dimension(:),   allocatable,   public   :: L1_fromN   !for each link, "from cell" id
    integer(i4),    dimension(:),   allocatable,   public   :: L1_toN     !for each link, "to cell" id
    integer(i4),    dimension(:),   allocatable,   public   :: L1_rOrder  !for each link, the routing order
    integer(i4),    dimension(:),   allocatable,   public   :: L1_netPerm !for each link, the order/sequence for calculation
    integer(i4),                                   public   :: OutNode

    integer(i4),    dimension(:),     allocatable,   public   :: L1_gaugeLoc !gauge location map

    logical,        dimension(:),     allocatable,   public   :: L1_chanmask        !dim = L1_nCells
    real(dp),       dimension(:),     allocatable,   public   :: L1_width        !dim = L1_nCells
    real(dp),       dimension(:),     allocatable,   public   :: L1_length        !dim = L1_nCells

    real(dp),       dimension(:),     allocatable,   public   :: L1_SrfmixRatio ! the weight of top-layer soil water conc to be mixed, for the surface flow conc

   !read semidistributed model domain information
    integer(i4),                                     public   :: L1_nSubs    !number of sub-catchments
    integer(i4),                                     public   :: L1_nHRUs    !total number of HRUs
    integer(i4),    dimension(:),   allocatable,     public   :: SubID       !subcatch ids
    real(dp),       dimension(:),   allocatable,   public     :: SubArea      !subcatch area [m2]
    real(dp),       dimension(:,:),   allocatable,   public   :: SubHRUs      !areal share of each HRU in each subcatch (dim1=subids,dim2=hrus) 
    integer(i4),    dimension(:),   allocatable,     public   :: CropRotationID !crop rotation IDs for each subcatchment
    type ReservoirStore
      !dim1 = number if cells, dim2=number of  simulation steps/ number of substances
      real(dp), dimension(:,:),  allocatable        :: reservoir  ! 
      real(dp), dimension(:,:),  allocatable        :: up_in       !vertical flux
      real(dp), dimension(:,:),  allocatable        :: down_out    !vertical flux
      real(dp), dimension(:,:),  allocatable        :: up_out      !vertical flux
      real(dp), dimension(:,:),  allocatable        :: down_in     !vertical flux
      real(dp), dimension(:,:),  allocatable        :: lateral_in  !lateral flux
      real(dp), dimension(:,:),  allocatable        :: lateral_out !lateral flux
    end type ReservoirStore
    type(ReservoirStore), public                                :: L1_WSintercept 
    type(ReservoirStore), public                                :: L1_WSsnowpack
    type(ReservoirStore), public                                :: L1_WSponding
    type(ReservoirStore), public                                :: L1_WSsurface 
    type(ReservoirStore), public                                :: L1_WSbase 
    type(ReservoirStore), dimension(:), allocatable, public     :: L1_WSsubsurface !dim=number of subsurface flow components
    type(ReservoirStore), dimension(:), allocatable, public     :: L1_WSsoil !dim = number of soillayer
   !dim1 = number if cells, dim2=number of timesteps
    real(dp), dimension(:,:),   allocatable, public :: L1_rainfall_direct
    real(dp), dimension(:,:),   allocatable, public :: L1_temp
    real(dp), dimension(:,:),   allocatable, public :: L1_throughfall
    real(dp), dimension(:,:),   allocatable, public :: L1_snowmelt
    real(dp), dimension(:,:),   allocatable, public :: L1_infiltrate
    real(dp), dimension(:,:),   allocatable, public :: L1_evapCanopy
    real(dp), dimension(:,:),   allocatable, public :: L1_surfaceRunoff
    real(dp), dimension(:,:),   allocatable, public :: L1_baseflow
    real(dp), dimension(:,:),   allocatable, public :: L1_reinfiltSrf   !surface reinfiltration
    real(dp), dimension(:,:,:), allocatable, public :: L1_reinfiltSoil   !reinfiltration for each soil layer
    real(dp), dimension(:,:),   allocatable, public :: L1_WSsealed
    real(dp), dimension(:,:),   allocatable, public :: L1_Esealed
    real(dp), dimension(:,:),   allocatable, public :: L1_directRunoff

   !dim1 = number if cells, dim2=number of timesteps, dim3 = number of subsurface flow components
    real(dp), dimension(:,:,:),   allocatable, public :: L1_subsurfaceRunoff 
   !dim1= ncells, dim2 = ntimestep, dim3 = number of involved soil layers
    real(dp), dimension(:,:,:),   allocatable, public :: L1_evapSoil
    real(dp), dimension(:,:,:),   allocatable, public :: L1_transp
   !channel, dim1 = number if cells, dim2=number of timesteps
    real(dp), dimension(:,:),   allocatable, public   :: L1_WSchannel
    real(dp), dimension(:,:),   allocatable, public   :: L1_WSchannel_upin
    real(dp), dimension(:,:),   allocatable, public   :: L1_WSchannel_downout
    real(dp), dimension(:,:),   allocatable, public   :: L1_dirTOchan
    real(dp), dimension(:,:),   allocatable, public   :: L1_srfTOchan
    real(dp), dimension(:,:),   allocatable, public   :: L1_basefTOchan
   !channel dim1 = dim1 = ncells, dim2 = ntimestep, dim3 = number OF SUBSURFACE runoffs
    real(dp), dimension(:,:,:),   allocatable, public :: L1_subsrfTOchan
    real(dp), dimension(:),     allocatable, public :: L1_avgtemp !average L1_temp if nLinks < nCells


!for WQ corresponding
     ! each variable has one new dimension for nsubstances
     ! no dimension for timesteps since all of them will be upddated for each time step
	 !dim1 = number if cells,dim2=number of nsubstances
    real(dp), dimension(:,:),   allocatable, public :: L1_crainfall_direct
    real(dp), dimension(:,:),   allocatable, public :: L1_cthroughfall
    real(dp), dimension(:,:),   allocatable, public :: L1_csnowmelt
    real(dp), dimension(:,:),   allocatable, public :: L1_csurfaceRunoff
    real(dp), dimension(:,:),   allocatable, public :: L1_cbaseflow
    real(dp), dimension(:,:),   allocatable, public :: L1_cpercolate     ! conc. in percolated water from unsat storage 	
    real(dp), dimension(:,:),   allocatable, public :: L1_cinfilsurface  !infiltration conc from surface to first soil layer
    real(dp), dimension(:,:),   allocatable, public :: L1_cdirectRunoff
    !the sealed storage only receives precip input and release evapor and direct runoff
    real(dp), dimension(:,:),   allocatable, public :: L1_cWSsealed 
	
   !dim1 = number if cells, dim2=number of nsubstances, dim3 = number of subsurface flow components	
    real(dp), dimension(:,:,:), allocatable, public :: L1_csubsurfaceRunoff
    real(dp), dimension(:,:,:), allocatable, public :: L1_cinfilSoil     ! conc. in infiltrated soil water
                                                                         !  (between each soillayer)
    real(dp), dimension(:,:,:), allocatable, public :: L1_creturnflow     ! conc. in return flow
                                                                         !  (between each soillayer)
    real(dp), dimension(:,:,:), allocatable, public :: L1_creinfiltSoil

    type(ReservoirStore), target,public                                :: L1_cWSintercept 
    type(ReservoirStore), target,public                                :: L1_cWSsnowpack 
    type(ReservoirStore), target,public                                :: L1_cWSponding
    type(ReservoirStore), target,public                                :: L1_cWSsurface 
    type(ReservoirStore), target,public                                :: L1_cWSbase 
    type(ReservoirStore), target,dimension(:), allocatable, public     :: L1_cWSsubsurface !dim=number of subsurface flow components
    type(ReservoirStore), target,dimension(:), allocatable, public     :: L1_cWSsoil !dim = number of soillayer
    !pointer for source input flux and concentration
    type(ReservoirStore), pointer,public                             :: L1_pcsurface !pointer for linking input flux and conc for surfaceflowWS conceptual storag
    type(ReservoirStore), pointer,public                             :: L1_pcbase !pointer for linking input flux and conc for baseflowWS conceptual storage	
    !type(ReservoirStore), pointer,public                             :: L1_pcsrfTOchan    
    !type(ReservoirStore), pointer,public                             :: L1_pcbasefTOchan
    !pointer array
    type pnt_Rsvr
        type(ReservoirStore), pointer :: Pnter_Rsvr
    end type pnt_Rsvr
    type(pnt_Rsvr), dimension(:),allocatable,public          :: L1_pcsubsurface !pointer for linking input flux and conc for subsurface conceptual storage
    !type(pnt_Rsvr), dimension(:),allocatable,public          :: L1_pcsubsrfTOchan
   !channel, dim1 = number if cells, dim2=number of nsubstances    
    real(dp), dimension(:,:),   allocatable, public                  :: L1_concdirTOchan
    real(dp), dimension(:,:),   allocatable, public                  :: L1_concsrfTOchan
    real(dp), dimension(:,:),   allocatable, public                  :: L1_concbasefTOchan
   !channel, dim1 = number if cells, dim2=number of nsubstances, dim3 = number OF SUBSURFACE runoffs
    real(dp), dimension(:,:,:),   allocatable, public                :: L1_concsubsrfTOchan


    real(dp), dimension(:),   allocatable, public :: L1_totalRunoffTochan ![mm]total runoff seepage to channel 
    real(dp), dimension(:),   allocatable, public :: L1_qOUT  ![m3/s] terrestrial flow export, considering also additional inflows
    real(dp), dimension(:,:), allocatable, public :: L1_concOUT ![mg/l]concentration of total ter. export
    real(dp), dimension(:,:), allocatable, public :: L1_terXpt        !total terrestrial export load to channel
    !nitrate parameters
    real(dp), dimension(:),     allocatable, public :: L1_rdegradN       ! degradation rate from humusN pool to fastN pool
    real(dp), dimension(:),     allocatable, public :: L1_rmineralN      ! mineralisation rate from fastN pool 
                                                                         ! to dissolved inorganic pool
    real(dp), dimension(:),     allocatable, public :: L1_rdissolN       ! dissolution rate from fastN pool 
                                                                         ! to dissolved organic pool
    real(dp), dimension(:),     allocatable, public :: L1_rdeniSoil      ! IN denitrification rate in soil water
    real(dp), dimension(:),     allocatable, public :: L1_rdeniAqtc    ! denitrification rate in aquatic system(river)
    real(dp), dimension(:),     allocatable, public :: L1_ratuptkN     ! autotropic N uptake mgNm-2d-1
    real(dp), dimension(:),     allocatable, public :: L1_rpprodN      ! assimilatory rate in aquatic system

	!values of land-use denpent parameters and initial values of different pools
	!map input --dim = ncells  
    real(dp), dimension(:), allocatable,  public   :: init_concIN    ! initial IN concentration in soil water
    real(dp), dimension(:), allocatable,  public   :: init_concON    ! initial ON concentration in soil water
    !dim=nclass
    integer(i4),  public       :: init_type   !inital value dependency: 1 -landuse, 2-soil type should be given in the file
    real(dp), dimension(:), allocatable,  public   :: init_humusN    ! amount of pool
    real(dp), dimension(:), allocatable,  public   :: init_fastN     ! amount of pool
    !real(dp), dimension(:), allocatable,  public   :: upsoil_part    ! fraction of plant uptake in first layer 
    real(dp), dimension(:), allocatable,  public   :: hnhalf         ! half depthe for humusN pool
    integer(i4), dimension(:), allocatable, public    :: Geoform     ! Formation 1 - permeable sedimentary materials
                                                                     !           2 - impermeable bedrock
    !variables as global arguements
    !state variables   dim1=number of cells, dim2=number of soil layers
    real(dp), dimension(:,:),   allocatable, public :: L1_Puptake    ! potential N uptake amount for each soil layer	
    real(dp), dimension(:,:),   allocatable, public :: L1_humusN         ! humus Nitrate pool in each soil layer(organic form)
    real(dp), dimension(:,:),   allocatable, public :: L1_fastN          ! fast Nitrate pool in each soil layer(orgainc form)
    real(dp), dimension(:,:),   allocatable, public :: L1_dissolvedIN    ! dissolved inorganic pool in each soil layer
    real(dp), dimension(:,:),   allocatable, public :: L1_dissolvedON    ! dissolved organic pool in each soil layer

    real(dp), dimension(:,:),     allocatable, public :: L1_soiltemp       ! soil temperature
    real(dp), dimension(:),       allocatable, public :: stat_soiltemp     !state variable for soil temperature if not externally supplied
    !matter fluxes
    real(dp), dimension(:),     allocatable, public :: L1_soilUptakeN    ! soil N uptake amount
    real(dp), dimension(:),     allocatable, public :: L1_soilDenitri    ! soil denitrification amount
    real(dp), dimension(:),     allocatable, public :: L1_soilMineralN    ! soil N mineralisation amount
    real(dp), dimension(:),     allocatable, public :: L1_soilINfrtmanapp ! soil IN surpplied amount

    real(dp), dimension(:),     allocatable, public :: L1_gwresidT       ! groundwater resident time baseflow conc(INCA EQ.)
   
  !state variables to store last step values (previous step) for soil water conc. 
    real(dp), dimension(:),     allocatable, public :: prevstep_WSponding   !surface ponding storage at previous step    
    real(dp), dimension(:),     allocatable, public :: prevstep_WSintercept    !interception storage at previous step    
    real(dp), dimension(:),     allocatable, public :: prevstep_WSsurface    ! surface ponding at previous step   
    real(dp),                                public :: ratio_RetentionWS  !the ratio of retention passive storage to hydrologically active storage in deeper groundwater
    real(dp), dimension(:),     allocatable, public :: prevstep_WSbase    ! DeepGW for baseflow at previous step 	
    real(dp), dimension(:,:),   allocatable, public :: prevstep_WSsubstorage    !subsurface storage for subsurface flow at pervious step  
    real(dp), dimension(:,:),   allocatable, public :: prevstep_WSsoil          ! soil water storage at previous step	
    real(dp), dimension(:),     allocatable, public :: prevstep_WSsealed   !sealed storage at previous step 

    real(dp), dimension(:),     allocatable, public :: stat_LATsurface  !lateral in of previous step
    real(dp), dimension(:),     allocatable, public :: stat_LATbase
    real(dp), dimension(:,:),   allocatable, public :: stat_LATsubstorage
    real(dp), dimension(:,:),   allocatable, public :: stat_LATsoil

    real(dp), dimension(:),     allocatable, public :: prevstep_baseflow          ! baseflow in previous step
    real(dp), dimension(:),     allocatable, public :: prevstep_percol
    real(dp), dimension(:),     allocatable, public :: prevstep_baseinflow
    real(dp), dimension(:),     allocatable, public :: prevstep_baseoutflow
    real(dp), dimension(:,:),   allocatable, public :: prevstep_cpercol
    real(dp), dimension(:,:),   allocatable, public :: prevstep_cbaseinflow
    real(dp), dimension(:,:),   allocatable, public :: prevstep_cbaseout

    integer(i4), dimension(:),     allocatable, public :: LAIUnitList
    real(dp), dimension(:,:),   allocatable, public :: LAILUT

  type rotationType
     ! input data
     integer(i4), dimension(:), allocatable        :: id        !         rotation Id
     integer(i4), dimension(:), allocatable        :: ncrops    !Number of crops in specific rotation type, maximum 10 types
     integer(i4), dimension(:,:), allocatable      :: crop      ! crop_ids in specific rotation type
     !                                                                !            
  end type rotationType
  type(rotationType), public :: rotation  ! rotation infos

  type cropType
     ! input data
     character(256)        :: cropname  ! name decription of crop
     integer(i4)           :: cropid    ! crop id
     ! data of input sources     
     real(dp)              :: frtn1       ! amount of nitrogen in first fertilizer
     integer(i4)           :: frtday1     ! number of days in a year that firstly apply fertilizer
     real(dp)              :: frtdown1    ! fraction of N that goes down to the second layers
     real(dp)              :: frtn2       ! amount of nitrogen in second fertilizer
     integer(i4)           :: frtday2     ! number of days in a year that secondly apply fertilizer
     real(dp)              :: frtdown2    ! fraction of N that goes down to the second layers
     real(dp)              :: mann1       ! amount of nitrogen in first manure
     integer(i4)           :: manday1     ! number of days in a year that firstly apply manure
     real(dp)              :: mandown1    ! fraction of N that goes down to the second layers
     real(dp)              :: mann2       ! amount of nitrogen in second manure
     integer(i4)           :: manday2     ! number of days in a year that secondly apply manure
     real(dp)              :: mandown2    ! fraction of N that goes down to the second layers
     real(dp)              :: manfIN       ! fraction of inorganic N in applied mannure
     integer(i4)           :: frtperiod   ! number of days to spread fert and man
     real(dp)              :: resn      ! amount of nitrogen in residuals
     integer(i4)           :: resday    ! number of days in a year that residuals appear
     real(dp)              :: resdown   ! fraction of N in residuals that goes down to the second layers
     real(dp)              :: resfast   ! fraction of N in residuals that goes into fast N pool
     integer(i4)           :: resperiod ! number of days to spread residuals
     ! uptake parameters
     real(dp)              :: up1
     real(dp)              :: up2
     real(dp)              :: up3
     real(dp)              :: uppsoil   ! fractions of the first layer 
     real(dp)              :: deepsoil  ! fractions from the deeper soil layer (the third layer)
     !here the uptake share of:
     !the first layer f = uppsoil * (1-deepsoil), second layer f = (1-uppsoil) * (1-deepsoil), the third f = deepsoil

     ! farming date
     integer(i4)           :: plantd    ! planting date (sowning date)
     integer(i4)           :: emergd   ! emerging date for winter crops
     integer(i4)           :: harvestd   ! havesting date
     ! catch crop
     integer(i4)           :: ccrop     ! 1= including catch crop; 0= non catch crop
     integer(i4)           :: ccplantd  ! planting date of catch crop
     integer(i4)           :: cchavestd ! havesting date of catch crop (can be more than 365 if it's winter catch crop 
	                                    ! and havest in the spring of next year)
  end type cropType
  type(cropType), dimension(:), allocatable, public      :: cropdata    ! crop managements 
	! -------------------------------------------------------------------

  !VARIABLES FOR IN-STREAM PROCESSES

  real(dp), dimension(:),       allocatable, public      :: nLink_rivertemp     !reach water temperature at L1 level                                                                            
  !state variables for calculating
  real(dp), dimension(:),       allocatable, public      :: nLink_rivert_avg10 ! 10-day's average temperature of river water in each reach 
  real(dp), dimension(:),       allocatable, public      :: nLink_rivert_avg20 ! 20-day's average temperature of river water in each reach    

  real(dp), dimension(:),       allocatable, public      :: nLink_riverbox     ! amount of water stored in each reach
  real(dp), dimension(:,:),     allocatable, public      :: nLink_criverbox    ! conc. of each riverbox
  
  real(dp), dimension(:),       allocatable, public      :: nLink_yravg_q    ! yearly average discharge of each reach [m3/s]
  real(dp), dimension(:),       allocatable, public      :: L1_qTIN          ! total input of q [m3/s] for each reach
  real(dp), dimension(:,:),     allocatable, public      :: L1_concTIN       ! conc. of total inflow in each reach node,
  real(dp), dimension(:,:),     allocatable, public      :: L1_loadTIN       ! [mg/l *m3/s * s]                                                                    
  ! state variables for uptake amount
  real(dp), dimension(:),       allocatable, public      :: L1_aquaticDenitri  ! N denitrified in stream water 
  real(dp), dimension(:),       allocatable, public      :: L1_aquaticAssimil  !N uptake (assimilatory) amount in stream water
  ! global radiation for in-stream Primarily production  
  logical, public                                        :: GR_file_exist = .FALSE.  ! if the "global_radiation.txt" presents or not
  real(dp), dimension(:),     allocatable, public      :: global_radiation 
                                                           ! measured global radiation data, dim1=nbasin,dim2=nosimulatedays
  real(dp), dimension(:),     allocatable, public      :: nor_globalradi   ! normalised global radiation
  real(dp), dimension(:,:),     allocatable, public      :: norlai_daily     ! normalised daily lai, dim1=nLAIclass, dim2=days_yr
 

  real(dp), dimension(:),   allocatable, public          :: L1_rzcoeff  !riparian zone shading coefficient  
                                                           ! dim = number of reaches 
  real(dp), dimension(:),   allocatable, public          :: L1_flight   ! moving average of L11_recoeff, for GPP calculation


  integer(i4), public                            :: nMeasPerDay ! Number of observations per day  
  integer(i4),dimension(:), allocatable, public  :: nEvalCmeasPerday     ! number of measured conc. data per day in files
  integer(i4),dimension(:), allocatable, public  :: nAddinCmeasPerday     ! number of measured conc. data per day in files
  integer(i4),dimension(:), allocatable, public  :: numEvalCol       ! number of wq data columns in evaluation data file  
  integer(i4),dimension(:), allocatable, public  :: numAddinCol      ! number of wq data columns in additional inflow data file


  character(256), dimension(:,:), allocatable, public  :: evalHead_str     ! name of data measured at evaluation station
  character(256), dimension(:,:), allocatable, public  :: inflowHead_str   ! name of data measured at inflow station   
  
  !for model write out
  real(dp), dimension(:,:), allocatable, public          :: Hydro_DIS !dim1 = observation, dim2=gauges
  real(dp), dimension(:,:,:), allocatable  , public      :: WQM_nutrient !dim3=number of measured data type  

  !State Fluxes variable for output
    integer(i4), public                             :: wncid_ter,wncid_chan,wNLONS,wNLATS,wNchans
    integer(i4), public                             :: nouter
    integer(i4), public                             :: out_counter
    real(dp), dimension(:,:,:), allocatable, public :: out_csoilMoist           !soil moisture conc
    real(dp), dimension(:,:),   allocatable, public :: out_cdirectrunoff        !directrunoff flow conc
    real(dp), dimension(:,:),   allocatable, public :: out_csurfaceRunoff       !surface flow conc
    real(dp), dimension(:,:,:), allocatable, public :: out_csubsurfaceRunoff    !subsurface flow conc
    real(dp), dimension(:,:),   allocatable, public :: out_cbaseflow            !baseflow conc
    real(dp), dimension(:,:),   allocatable, public :: out_soilstoreN          !soil stock N dim1 = ncells, dim2 = nsoillayers
    real(dp), dimension(:),     allocatable, public :: out_soilUptakeN         !soil plant/crop actual N uptake
    real(dp), dimension(:),     allocatable, public :: out_soilDenitri         !soil denitrification
    real(dp), dimension(:),     allocatable, public :: out_soilMineralN        !soil mineralization
    real(dp), dimension(:),     allocatable, public :: out_soilINfrtmanapp     !N external input (fertilizer, residues)
    real(dp), dimension(:,:),   allocatable, public :: out_concdirTOchan       !channel seepage surface
    real(dp), dimension(:,:),   allocatable, public :: out_concsrfTOchan       !channel seepage surface
    real(dp), dimension(:,:,:), allocatable, public :: out_concsubsrfTOchan    !channel seepage subsurface
    real(dp), dimension(:,:),   allocatable, public :: out_concbasefTOchan     !channel seepage baseflow
    real(dp), dimension(:,:),   allocatable, public :: out_concOUT             !total runoff generation to channel
    real(dp), dimension(:,:),   allocatable, public :: out_concTIN             !stream water Nitrate-N concentration
    real(dp), dimension(:,:),   allocatable, public :: out_terXpt              !terrestrial export load
    real(dp), dimension(:),     allocatable, public :: out_aquaticDenitri     !in-stream denitrification
    real(dp), dimension(:),     allocatable, public :: out_aquaticAssimil       !in-stream assimilatory uptake  



  integer(i4)                                            :: timeStep_model_outputs_wqm ! timestep for writing model outputs
  integer(i4),      public                               :: nOUTstate_wqm       ! total number of states that can be written out
  logical, dimension(:), allocatable, public             :: outputFlxState_wqm         ! Define model outputs see "wqm_outputs.nml"
  character(256), parameter, public                      :: Version_WQM = '2.0'          ! model version

  ! !directly introduced from mHM implementations
   ! ! -------------------------------------------------------------------
   ! ! OPTIMIZATION
   ! ! -------------------------------------------------------------------
   ! integer(i4), public                              :: opti_method         ! Optimization algorithm:
   ! !                                                                       ! 1 - DDS
   ! !                                                                       ! 2 - Simulated Annealing
   ! !                                                                       ! 3 - SCE
   ! integer(i4), public                              :: opti_function       ! Objective function:
   ! !                                                                       ! 1 - 1.0-NSE
   ! !                                                                       ! 2 - 1.0-lnNSE
   ! !                                                                       ! 3 - 1.0-0.5*(NSE+lnNSE)
   ! logical,     public                              :: optimize            ! Optimization   (.true. ) or
   ! !                                                                       ! Evaluation run (.false.)
   ! logical,     public                              :: optimize_restart    ! Optimization will be restarted from
   ! !                                                                       ! mo_<opti_method>.restart file (.true.)
   ! ! settings for optimization algorithms: 
   ! integer(i8), public                              :: seed                ! seed used for optimization
   ! !                                                                       ! default: -9 --> system time 
   ! integer(i4), public                              :: nIterations         ! number of iterations for optimization
   ! real(dp),    public                              :: dds_r               ! DDS: perturbation rate
   ! !                                                                       !      default: 0.2
   ! real(dp),    public                              :: sa_temp             ! SA:  initial temperature
   ! !                                                                       !      default: -9.0 --> estimated
   ! integer(i4), public                              :: sce_ngs             ! SCE: # of complexes
   ! !                                                                       !      default: 2
   ! integer(i4), public                              :: sce_npg             ! SCE: # of points per complex
   ! !                                                                       !      default: -9 --> 2n+1
   ! integer(i4), public                              :: sce_nps             ! SCE: # of points per subcomplex
   ! !                                                                       !      default: -9 --> n+1
   ! logical,     public                              :: mcmc_opti           ! MCMC: Optimization (.true. ) or
   ! !                                                                       !       Only parameter uncertainty (.false.)
   ! integer(i4), public, parameter                   :: nerror_model = 2    !       # possible parameters in error model
   ! !                                                                       !       e.g. for opti_function=8: 2
   ! real(dp),    public, dimension(nerror_model)     :: mcmc_error_params   !       Parameters of error model if mcmc_opti=.false.
   ! !                                                                       !       e.g. for opti_function=8: 0.01, 0.3  
  
  
END MODULE mo_wqm_global_variables