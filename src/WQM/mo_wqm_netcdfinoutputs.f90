!> \file mo_wqm_netcdfinoutputs.f90

!> \details NetCDF is first initialized and later on variables are put to the NetCDF.

!  HISTORY
!>     \authors Xiaoqiang Yang
!>     \date Apr 2022

module mo_wqm_netcdfinoutputs

  use mo_kind                 , only: i4, dp
  use mo_string_utils         , only: num2str
  use mo_wqm_global_variables     , only: gridGeoRef
  use mo_wqm_constants        , only: nodata_dp,nodata_i4
  use netcdf

  implicit none

  PUBLIC :: ncfile_create
  PUBLIC :: ncfile_update
  PUBLIC :: ncfile_getDim
  PUBLIC :: ncfile_read


  interface ncfile_read
    module procedure ncfile_read_2d, ncfile_read_3d
  end interface ncfile_read

  CONTAINS

  !------------------------------------------------------------------
  !     NAME
  !         ncfile_create
  !
  !     HISTORY
  !>        \author Xiaoqiang Yang
  !>        \date Apr 2022
  
  subroutine ncfile_create()

    use mo_wqm_global_variables,  only : dir_Out, model_type,m_Level,L1_nHRUs, L1_nSubs,&
        L1_chanmask,ter_variables, chan_variables, &
        outputFlxState_wqm, nSoilLayers,nSubsurfaceflow,nSubsrfToChannel, &
        L1_nCells, L1_nLinks
    implicit none

    ! local
    integer(i4)             :: ncid_ter,ncid_chan            !nc file ids
    integer(i4)             :: lat_dimid,lon_dimid,rec_dimid,chan_dimid !dimension ids
    integer(i4)             :: lat_varid,lon_varid,chan_varid !variable ids
    character(256)      :: fname,fname_chan
    character(256)      :: LAT_NAME,LON_NAME,REC_NAME,chanid_NAME
    character(256),parameter  :: UNITS = "units"
    character(256)      :: varUNITS,LAT_UNITS,LON_UNITS,chanid_UNITS
    integer(i4)             :: NLATS,NLONS,Nchans
    integer(i4)             :: nr_tervar, nr_chnvar
    real(dp),dimension(:),allocatable  :: lats, lons,chanids
    integer(i4)             ::dimids(3),dimids_chan(2),nstart(3), ncount(3)
    integer(i4)             :: ii, nn,jj
    real(dp),dimension(:,:,:), allocatable :: fillval


    nr_tervar = count(outputFlxState_wqm(1:10))
    print*, "Total nc output variables: ", nr_tervar
    nr_chnvar = count(outputFlxState_wqm(11:19))
    ! write terrestrial variables
    if (nr_tervar .gt. 0_i4) then
      if (model_type == 1_i4) then 
        LAT_NAME = "northing"
        LON_NAME = "easting"
        NLATS =  m_Level%nrows
        NLONS =  m_Level%ncols
        LAT_UNITS = "degrees_north"
        LON_UNITS = "degrees_east"
        allocate(lats(NLATS))
        allocate(lons(NLONS))
        lats = nodata_dp
        lons = nodata_dp
        lons(1) =  m_Level%xllcorner + 0.5_dp *  m_Level%cellsize
        do jj = 2, NLONS
            lons(jj)   =  lons(jj-1) + m_Level%cellsize
        end do
        lats( m_Level%nrows) =  m_Level%yllcorner + 0.5_dp * m_Level%cellsize
        do jj = NLATS-1,1,-1
            lats(jj)   =  lats(jj+1) + m_Level%cellsize
        end do
      else if (model_type == 2_i4) then
        LAT_NAME = "hrus"
        LON_NAME = "subs"
        NLATS =  L1_nSubs
        NLONS =  L1_nHRUs
        LAT_UNITS = "no_hrus"
        LON_UNITS = "no_subcatchment"
        allocate(lats(NLATS))
        allocate(lons(NLONS))
        lats = nodata_dp
        lons = nodata_dp
        do jj = 1, NLONS
            lons(jj)   =  jj
        end do
        do jj = NLATS,1,-1
            lats(jj)   =  jj
        end do
      else
        print*, "Please specify >model_type< as 1 or 2!"
      end if
      !dimension time
      REC_NAME = "time"

      fname = trim(dir_Out) // 'WQMter_Fluxes_States.nc'
      call check( nf90_create(fname, nf90_clobber, ncid_ter) )

    ! Define the dimensions.
      call check( nf90_def_dim(ncid_ter, LAT_NAME, NLATS, lat_dimid) )
      call check( nf90_def_dim(ncid_ter, LON_NAME, NLONS, lon_dimid) )
      call check( nf90_def_dim(ncid_ter, REC_NAME, NF90_UNLIMITED, rec_dimid) )
    ! Define the coordinate variables. They will hold the coordinate
    ! information, that is, the latitudes and longitudes. A varid is
    ! returned for each.
      call check( nf90_def_var(ncid_ter, LAT_NAME, NF90_DOUBLE, lat_dimid, lat_varid) )
      call check( nf90_def_var(ncid_ter, LON_NAME, NF90_DOUBLE, lon_dimid, lon_varid) )
    ! Assign units attributes to coordinate var data. This attaches a
    ! text attribute to each of the coordinate variables, containing the
    ! units.
      call check( nf90_put_att(ncid_ter, lat_varid, UNITS, LAT_UNITS) )
      call check( nf90_put_att(ncid_ter, lon_varid, UNITS, LON_UNITS) )
    ! The dimids array is used to pass the dimids of the dimensions of
    ! the netCDF variables.
      dimids = (/ lon_dimid, lat_dimid, rec_dimid /)
    !check if output of specific variables is specified 

     if (allocated(ter_variables)) deallocate(ter_variables) 
     allocate(ter_variables(nr_tervar* nSoilLayers))
     ter_variables(:)%varids = nodata_i4
     ter_variables(:)%var_NAME = "varnames"
     ter_variables(:)%var_UNIT = "varunits"
 
     ii = 0
     if (outputFlxState_wqm(1)) then
      do nn = 1, nSoilLayers
        ii = ii+1
        ter_variables(ii)%var_NAME = "SM_conc"//trim(num2str(nn, '(i2.2)'))
        ter_variables(ii)%var_UNIT = "mg l-1"
        ! Define the netCDF variables for the pressure and temperature data.
        call check( nf90_def_var(ncid_ter, ter_variables(ii)%var_NAME, NF90_DOUBLE, dimids, ter_variables(ii)%varids))
        call check( nf90_put_att(ncid_ter, ter_variables(ii)%varids, UNITS, ter_variables(ii)%var_UNIT))
      end do
     end if

     if (outputFlxState_wqm(2)) then
      ii = ii + 1
      ter_variables(ii)%var_NAME = "cdirectRunoff"
      ter_variables(ii)%var_UNIT = "mg l-1"
      ! Define the netCDF variables for the pressure and temperature data.
      call check( nf90_def_var(ncid_ter, ter_variables(ii)%var_NAME, NF90_DOUBLE, dimids, ter_variables(ii)%varids))
      call check( nf90_put_att(ncid_ter, ter_variables(ii)%varids, UNITS, ter_variables(ii)%var_UNIT))
     end if


     if (outputFlxState_wqm(3)) then
      ii = ii + 1
      ter_variables(ii)%var_NAME = "csrfRunoff"
      ter_variables(ii)%var_UNIT = "mg l-1"
      ! Define the netCDF variables for the pressure and temperature data.
      call check( nf90_def_var(ncid_ter, ter_variables(ii)%var_NAME, NF90_DOUBLE, dimids, ter_variables(ii)%varids))
      call check( nf90_put_att(ncid_ter, ter_variables(ii)%varids, UNITS, ter_variables(ii)%var_UNIT))
     end if

     if (outputFlxState_wqm(4)) then
       do nn = 1, nSubsurfaceflow
         ii = ii + 1
         ter_variables(ii)%var_NAME = "csubsrfRunoff"//trim(num2str(nn, '(i2.2)'))
         ter_variables(ii)%var_UNIT = "mg l-1"
         ! Define the netCDF variables for the pressure and temperature data.
         call check( nf90_def_var(ncid_ter, ter_variables(ii)%var_NAME, NF90_DOUBLE, dimids, ter_variables(ii)%varids))
         call check( nf90_put_att(ncid_ter, ter_variables(ii)%varids, UNITS, ter_variables(ii)%var_UNIT))
       end do
     end if



     if (outputFlxState_wqm(5)) then
      ii = ii + 1
      ter_variables(ii)%var_NAME = "cbasefow"
      ter_variables(ii)%var_UNIT = "mg l-1"
      ! Define the netCDF variables for the pressure and temperature data.
      call check( nf90_def_var(ncid_ter, ter_variables(ii)%var_NAME, NF90_DOUBLE, dimids, ter_variables(ii)%varids))
      call check( nf90_put_att(ncid_ter, ter_variables(ii)%varids, UNITS, ter_variables(ii)%var_UNIT))
     end if
	
     if (outputFlxState_wqm(6)) then
       do nn = 1, nSoilLayers
         ii = ii + 1
         ter_variables(ii)%var_NAME = "soilstockN"//trim(num2str(nn, '(i2.2)'))
         ter_variables(ii)%var_UNIT = "g m-2"
         ! Define the netCDF variables for the pressure and temperature data.
         call check( nf90_def_var(ncid_ter, ter_variables(ii)%var_NAME, NF90_DOUBLE, dimids, ter_variables(ii)%varids))
         call check( nf90_put_att(ncid_ter, ter_variables(ii)%varids, UNITS, ter_variables(ii)%var_UNIT))
       end do
     end if

     if (outputFlxState_wqm(7)) then
      ii = ii + 1
      ter_variables(ii)%var_NAME = "soiluptakeN"
      ter_variables(ii)%var_UNIT = "g m-2"
      ! Define the netCDF variables for the pressure and temperature data.
      call check( nf90_def_var(ncid_ter, ter_variables(ii)%var_NAME, NF90_DOUBLE, dimids, ter_variables(ii)%varids))
      call check( nf90_put_att(ncid_ter, ter_variables(ii)%varids, UNITS, ter_variables(ii)%var_UNIT))
     end if
	
     if (outputFlxState_wqm(8)) then
       ii = ii + 1
       ter_variables(ii)%var_NAME ="soildenitri"
       ter_variables(ii)%var_UNIT = "g m-2"
      ! Define the netCDF variables for the pressure and temperature data.
       call check( nf90_def_var(ncid_ter, ter_variables(ii)%var_NAME, NF90_DOUBLE, dimids, ter_variables(ii)%varids))
       call check( nf90_put_att(ncid_ter, ter_variables(ii)%varids, UNITS, ter_variables(ii)%var_UNIT))
     end if



     if (outputFlxState_wqm(9)) then
       ii = ii + 1
       ter_variables(ii)%var_NAME ="soilmineralN"
       ter_variables(ii)%var_UNIT = "g m-2"
      ! Define the netCDF variables for the pressure and temperature data.
       call check( nf90_def_var(ncid_ter, ter_variables(ii)%var_NAME, NF90_DOUBLE, dimids, ter_variables(ii)%varids))
       call check( nf90_put_att(ncid_ter, ter_variables(ii)%varids, UNITS, ter_variables(ii)%var_UNIT))
     end if

     if (outputFlxState_wqm(10)) then
       ii = ii + 1
       ter_variables(ii)%var_NAME ="INfrtmanapp"
       ter_variables(ii)%var_UNIT = "g m-2"
      ! Define the netCDF variables for the pressure and temperature data.
       call check( nf90_def_var(ncid_ter, ter_variables(ii)%var_NAME, NF90_DOUBLE, dimids, ter_variables(ii)%varids))
       call check( nf90_put_att(ncid_ter, ter_variables(ii)%varids, UNITS, ter_variables(ii)%var_UNIT))
     end if
    ! End define mode.
     call check( nf90_enddef(ncid_ter) )
    ! put lon and lat values

    ! Write the coordinate variable data. This will put the latitudes
    ! and longitudes of our data grid into the netCDF file.
     call check( nf90_put_var(ncid_ter, lat_varid, lats) )
     call check( nf90_put_var(ncid_ter, lon_varid, lons) )
    ! ! initial fill values
     ! ncount = (/ NLONS, NLATS, 1 /)
     ! nstart = (/ 1, 1, 1/)
     ! allocate(fillval(NLONS,NLATS,1))
     ! fillval = nodata_dp
     ! do jj=1,ii
       ! call check( nf90_put_var(ncid_ter, ter_variables(jj)%varids, fillval,start = nstart, &
                     ! count = ncount) )
     ! end do

    !close the file. This causes netCDF to flush all buffers 
     call check( nf90_close(ncid_ter) )
    end if


    !for channel related variables
    if (nr_chnvar .gt. 0_i4) then
      chanid_NAME = "chanIDs"
      Nchans = size(L1_chanmask)
      chanid_UNITS = "grid_location"
      !dimension time
      REC_NAME = "time"

      allocate(chanids(Nchans))
      do ii=1,Nchans
        chanids(ii) = ii
      end do
      if (.not.allocated(chan_variables)) allocate(chan_variables(nr_chnvar* nSubsrfToChannel))
      chan_variables(:)%varids = nodata_i4
      chan_variables(:)%var_NAME = "varnames"
      chan_variables(:)%var_UNIT = "varunits"
      fname_chan = trim(dir_Out) // 'WQMchan_Fluxes_States.nc'
      call check( nf90_create(fname_chan, nf90_clobber, ncid_chan) )
    ! Define the dimensions.
      call check( nf90_def_dim(ncid_chan, chanid_NAME, Nchans, chan_dimid) )
      call check( nf90_def_dim(ncid_chan, REC_NAME, NF90_UNLIMITED, rec_dimid) )
    ! Define the coordinate variables. They will hold the coordinate
    ! information, that is, the latitudes and longitudes. A varid is
    ! returned for each.
      call check( nf90_def_var(ncid_chan, chanid_NAME, NF90_DOUBLE, chan_dimid, chan_varid) )
    ! Assign units attributes to coordinate var data. This attaches a
    ! text attribute to each of the coordinate variables, containing the
    ! units.
      call check( nf90_put_att(ncid_chan, chan_varid, UNITS, chanid_UNITS) )

    ! The dimids array is used to pass the dimids of the dimensions of
    ! the netCDF variables.
      dimids_chan = (/ chan_dimid,rec_dimid /)
    !check if output of specific variables is specified 
     ii = 0
    if (outputFlxState_wqm(11)) then
       ii = ii + 1     
       chan_variables(ii)%var_NAME ="cdirTochan"
       chan_variables(ii)%var_UNIT = "mg l-1"
      ! Define the netCDF variables for the pressure and temperature data.
       call check( nf90_def_var(ncid_chan, chan_variables(ii)%var_NAME, NF90_DOUBLE, dimids_chan, chan_variables(ii)%varids))
       call check( nf90_put_att(ncid_chan, chan_variables(ii)%varids, UNITS, chan_variables(ii)%var_UNIT))
    end if
    if (outputFlxState_wqm(12)) then
       ii = ii + 1     
       chan_variables(ii)%var_NAME ="csrfTochan"
       chan_variables(ii)%var_UNIT = "mg l-1"
      ! Define the netCDF variables for the pressure and temperature data.
       call check( nf90_def_var(ncid_chan, chan_variables(ii)%var_NAME, NF90_DOUBLE, dimids_chan, chan_variables(ii)%varids))
       call check( nf90_put_att(ncid_chan, chan_variables(ii)%varids, UNITS, chan_variables(ii)%var_UNIT))
    end if

    if (outputFlxState_wqm(13)) then
       do nn = 1, nSubsrfToChannel
         ii = ii + 1
         chan_variables(ii)%var_NAME ="csubsrfTochan"//trim(num2str(nn, '(i2.2)'))
         chan_variables(ii)%var_UNIT = "mg l-1"
        ! Define the netCDF variables for the pressure and temperature data.
         call check( nf90_def_var(ncid_chan, chan_variables(ii)%var_NAME, NF90_DOUBLE, dimids_chan, chan_variables(ii)%varids))
         call check( nf90_put_att(ncid_chan, chan_variables(ii)%varids, UNITS, chan_variables(ii)%var_UNIT))
       end do
    end if

    if (outputFlxState_wqm(14)) then
       ii = ii + 1
       chan_variables(ii)%var_NAME ="cbaseTochan"
       chan_variables(ii)%var_UNIT = "mg l-1"
      ! Define the netCDF variables for the pressure and temperature data.
       call check( nf90_def_var(ncid_chan, chan_variables(ii)%var_NAME, NF90_DOUBLE, dimids_chan, chan_variables(ii)%varids))
       call check( nf90_put_att(ncid_chan, chan_variables(ii)%varids, UNITS, chan_variables(ii)%var_UNIT))
    end if


    if (outputFlxState_wqm(15)) then
       ii = ii + 1
       chan_variables(ii)%var_NAME ="ctrunoffTochan"
       chan_variables(ii)%var_UNIT = "mg l-1"
      ! Define the netCDF variables for the pressure and temperature data.
       call check( nf90_def_var(ncid_chan, chan_variables(ii)%var_NAME, NF90_DOUBLE, dimids_chan, chan_variables(ii)%varids))
       call check( nf90_put_att(ncid_chan, chan_variables(ii)%varids, UNITS, chan_variables(ii)%var_UNIT))
    end if

    if (outputFlxState_wqm(16)) then
       ii = ii + 1
       chan_variables(ii)%var_NAME ="cstreamwater"
       chan_variables(ii)%var_UNIT = "mg l-1"
      ! Define the netCDF variables for the pressure and temperature data.
       call check( nf90_def_var(ncid_chan, chan_variables(ii)%var_NAME, NF90_DOUBLE, dimids_chan, chan_variables(ii)%varids))
       call check( nf90_put_att(ncid_chan, chan_variables(ii)%varids, UNITS, chan_variables(ii)%var_UNIT))
    end if

    if (outputFlxState_wqm(17)) then
       ii = ii + 1
       chan_variables(ii)%var_NAME ="Terexpt"
       chan_variables(ii)%var_UNIT = "g m-2"
      ! Define the netCDF variables for the pressure and temperature data.
       call check( nf90_def_var(ncid_chan, chan_variables(ii)%var_NAME, NF90_DOUBLE, dimids_chan, chan_variables(ii)%varids))
       call check( nf90_put_att(ncid_chan, chan_variables(ii)%varids, UNITS, chan_variables(ii)%var_UNIT))
    end if

    if (outputFlxState_wqm(18)) then
       ii = ii + 1
       chan_variables(ii)%var_NAME ="waterdenitri"
       chan_variables(ii)%var_UNIT = "kg d-1"
      ! Define the netCDF variables for the pressure and temperature data.
       call check( nf90_def_var(ncid_chan, chan_variables(ii)%var_NAME, NF90_DOUBLE, dimids_chan, chan_variables(ii)%varids))
       call check( nf90_put_att(ncid_chan, chan_variables(ii)%varids, UNITS, chan_variables(ii)%var_UNIT))
    end if
    if (outputFlxState_wqm(19)) then
       ii = ii + 1
       chan_variables(ii)%var_NAME ="waterassimi"
       chan_variables(ii)%var_UNIT = "kg d-1"
      ! Define the netCDF variables for the pressure and temperature data.
       call check( nf90_def_var(ncid_chan, chan_variables(ii)%var_NAME, NF90_DOUBLE, dimids_chan, chan_variables(ii)%varids))
       call check( nf90_put_att(ncid_chan, chan_variables(ii)%varids, UNITS, chan_variables(ii)%var_UNIT))
    end if
    ! End define mode.
    call check( nf90_enddef(ncid_chan) )

    ! Write the coordinate variable data. This will put the chanids of our data grid into the netCDF file.
    call check( nf90_put_var(ncid_chan, chan_varid, chanids) )
    call check( nf90_close(ncid_chan) )
   end if

  end subroutine ncfile_create
!-----------------------------------------------------------------
!update record values/different variables
  !------------------------------------------------------------------
  !     NAME
  !         ncfile_update
  !
  !     HISTORY
  !>        \author Xiaoqiang Yang
  !>        \date Apr 2022
  
  subroutine ncfile_update(step_loc,nt,nTotal)

    use mo_wqm_constants, only: nodata_dp
    use mo_wqm_global_variables,  only : dir_Out, model_type,m_Level,L1_nHRUs, L1_nSubs,&
        L1_chanmask,ter_variables, chan_variables, &
        outputFlxState_wqm, nSoilLayers,nSubsurfaceflow,nSubsrfToChannel, &
        L1_nCells, L1_nLinks, &
        wncid_ter,wncid_chan,wNLONS,wNLATS,wNchans, &
        Model_mask, out_counter, &
          out_csoilMoist            , &  !soil moisture conc
          out_cdirectrunoff      , & !direct runoff conc
          out_csurfaceRunoff    , &  !surface flow conc
          out_csubsurfaceRunoff  , &  !subsurface flow conc
          out_cbaseflow          , &  !baseflow conc
          out_soilstoreN,    &        !soil N stock
          out_soilUptakeN  , &        !soil plant/crop actual N uptake
          out_soilDenitri  , &        !soil denitrification
          out_soilMineralN , &        !soil mineralization
          out_soilINfrtmanapp, &      !N external input (fertilizer, residues)
          out_concdirTOchan,   &      !directrunoff seepage
          out_concsrfTOchan      , &  !channel seepage surface
          out_concsubsrfTOchan   , &  !channel seepage subsurface
          out_concbasefTOchan    , &  !channel seepage baseflow
          out_concOUT            , &  !total runoff generation to channel
          out_concTIN         , &    !stream water Nitrate-N concentration
          out_terXpt          , &    !terrestrial export load
          out_aquaticDenitri  , &    !in-stream denitrification
          out_aquaticAssimil     !in-stream assimilatory uptake

    implicit none

    integer(i4), intent(in) :: step_loc
    integer(i4), intent(in) :: nt
    integer(i4), intent(in) :: nTotal
    ! local
    !integer(i4)             ::            !nc file ids
    !integer(i4)             :: NLATS,NLONS,Nchans
    integer(i4)             :: nr_tervar, nr_chnvar,varids
    real(dp), dimension(L1_nCells)        :: tmpout
    !real(dp), dimension(L1_nLinks)        :: tmpout_chan
    real(dp),dimension(:,:,:),allocatable   :: ncvout
    real(dp),dimension(:,:),allocatable     :: ncvout_chan
    integer(i4) :: nstart(3), ncount(3),nstart_chan(2), ncount_chan(2)
    character(256)      :: fname,fname_chan
    integer(i4)             :: ii, nn,jj

    !if it's the first time call this update write,
    !then, have to open the newly created nc file first

      nr_tervar = count(outputFlxState_wqm(1:10))
      nr_chnvar = count(outputFlxState_wqm(11:19))
    if (step_loc == 1_i4) then      
     ! write terrestrial variables
      if (nr_tervar > 0_i4) then
        if (model_type == 1) then 
          wNLATS =  m_Level%nrows
          wNLONS =  m_Level%ncols
        else
          wNLATS =  L1_nSubs
          wNLONS =  L1_nHRUs
        end if

        fname = trim(dir_Out) // 'WQMter_Fluxes_States.nc'
        ! Open the file. 

        call check( nf90_open(fname, nf90_write, wncid_ter) )



      end if
      !for channel related variables
      if (nr_chnvar > 0_i4) then
        wNchans = size(L1_chanmask)
        fname_chan = trim(dir_Out) // 'WQMchan_Fluxes_States.nc'
        ! Open the file. 
        call check( nf90_open(fname_chan, nf90_write, wncid_chan) )
      end if
    end if

    ! These settings tell netcdf to write one timestep of data.
    if (nr_tervar > 0_i4) then 
	
      ncount = (/ wNLONS, wNLATS, 1 /)
      nstart = (/ 1, 1, step_loc/)

      if (.not. allocated(ncvout)) allocate(ncvout(wNLONS,wNLATS,nr_tervar*nSoilLayers))
      ncvout = nodata_dp

      !check if output of specific variables is specified 
      ii = 0
      if (outputFlxState_wqm(1)) then
        do nn = 1, nSoilLayers
          ii = ii+1
          tmpout = out_csoilMoist(:,nn,1)/out_counter     !mean concentration 
          ncvout(:,:,ii) = unpack(tmpout, Model_mask,nodata_dp)
        end do
      end if

      if (outputFlxState_wqm(2)) then
        tmpout = out_cdirectrunoff(:,1)/out_counter     !mean concentration 
        ii = ii+1
        ncvout(:,:,ii) = unpack(tmpout, Model_mask,nodata_dp)
      end if

      if (outputFlxState_wqm(3)) then
        tmpout = out_csurfaceRunoff(:,1)/out_counter     !mean concentration 
        ii = ii+1
        ncvout(:,:,ii) = unpack(tmpout, Model_mask,nodata_dp)
      end if

      if (outputFlxState_wqm(4)) then
        do nn = 1, nSubsurfaceflow
          tmpout = out_csubsurfaceRunoff(:,nn,1)/out_counter     !mean concentration 
          ii = ii+1
          ncvout(:,:,ii) = unpack(tmpout, Model_mask,nodata_dp)
        end do
      end if

      if (outputFlxState_wqm(5)) then
        tmpout = out_cbaseflow(:,1)/out_counter     !mean concentration 
        ii = ii+1
        ncvout(:,:,ii) = unpack(tmpout, Model_mask,nodata_dp)
      end if

      if (outputFlxState_wqm(6)) then
        do nn = 1, nSoilLayers
          ii = ii+1
          tmpout = out_soilstoreN(:,nn)/out_counter     !mean  
          ncvout(:,:,ii) = unpack(tmpout, Model_mask,nodata_dp)
        end do
      end if
	
      if (outputFlxState_wqm(7)) then
        tmpout = out_soilUptakeN    !accumulated amount
        ii = ii+1
        ncvout(:,:,ii) = unpack(tmpout, Model_mask,nodata_dp)
      end if
	
      if (outputFlxState_wqm(8)) then
        tmpout = out_soilDenitri    !accumulated amount
        ii = ii+1
        ncvout(:,:,ii) = unpack(tmpout, Model_mask,nodata_dp)
      end if

      if (outputFlxState_wqm(9)) then
        tmpout = out_soilMineralN    !accumulated amount
        ii = ii+1
        ncvout(:,:,ii) = unpack(tmpout, Model_mask,nodata_dp)
      end if

      if (outputFlxState_wqm(10)) then
        tmpout = out_soilINfrtmanapp    !accumulated amount
        ii = ii+1
        ncvout(:,:,ii) = unpack(tmpout, Model_mask,nodata_dp)
      end if
     ! Get the varids of the pressure and temperature netCDF variables.
     ! Write the pretend data ncvout.
     ! print*, step_loc, ii
      do jj=1,ii
        call check( nf90_inq_varid(wncid_ter, ter_variables(jj)%var_NAME, ter_variables(jj)%varids) )
        !print*, ncvout(4,35,jj)
        call check( nf90_put_var(wncid_ter, ter_variables(jj)%varids, ncvout(:,:,jj), start = nstart, &
                    count = ncount) )
      end do

      if (nt == nTotal) then
        ! Close the file. This causes netCDF to flush all buffers and make
        ! sure your data are really written to disk.
        call check( nf90_close(wncid_ter) ) 
      end if 

    end if
    !channel internal output
    if (nr_chnvar > 0_i4) then
      ncount_chan = (/ wNchans, 1 /)
      nstart_chan = (/ 1, step_loc/)
      if (.not. allocated(ncvout_chan)) allocate(ncvout_chan(wNchans,nr_chnvar*nSubsrfToChannel))
      ncvout_chan = nodata_dp

      ii = 0_i4

      if (outputFlxState_wqm(11)) then
       ii = ii + 1 
       ncvout_chan(:,ii) = out_concdirTOchan(:,1)/out_counter
      end if

      if (outputFlxState_wqm(12)) then
       ii = ii + 1 
       ncvout_chan(:,ii) = out_concsrfTOchan(:,1)/out_counter
      end if

      if (outputFlxState_wqm(13)) then
        do nn = 1, nSubsrfToChannel
          ii = ii + 1
          ncvout_chan(:,ii) = out_concsubsrfTOchan(:,nn,1)/out_counter
        end do
      end if

      if (outputFlxState_wqm(14)) then
        ii = ii + 1
        ncvout_chan(:,ii) = out_concbasefTOchan(:,1)/out_counter
      end if

      if (outputFlxState_wqm(15)) then
        ii = ii + 1
        ncvout_chan(:,ii) = out_concOUT(:,1)/out_counter
      end if

      if (outputFlxState_wqm(16)) then
        ii = ii + 1
        ncvout_chan(:,ii) = out_concTIN(:,1)/out_counter
      end if

      if (outputFlxState_wqm(17)) then
        ii = ii + 1
        ncvout_chan(:,ii) = out_terXpt(:,1)
      end if

      if (outputFlxState_wqm(18)) then
        ii = ii + 1
        ncvout_chan(:,ii) = out_aquaticDenitri(:)
      end if
      if (outputFlxState_wqm(19)) then
        ii = ii + 1
        ncvout_chan(:,ii) = out_aquaticAssimil(:)
      end if

     ! Get the varids of the pressure and temperature netCDF variables.
     ! Write the pretend data ncvout.

      do jj=1,ii
        call check( nf90_inq_varid(wncid_chan, chan_variables(jj)%var_NAME, chan_variables(jj)%varids) )
        call check( nf90_put_var(wncid_chan, chan_variables(jj)%varids, ncvout_chan(:,jj), start = nstart_chan, &
                    count = ncount_chan) )
      end do

      if (nt == nTotal) then
        ! Close the file. This causes netCDF to flush all buffers and make
        ! sure your data are really written to disk.
        call check( nf90_close(wncid_chan) ) 
      end if 

    end if

  end subroutine ncfile_update

!-----------------------------------------------------------------
!get hydro input nc file dimension information
  !------------------------------------------------------------------
  !     NAME
  !         ncfile_getDim
  !
  !     HISTORY
  !>        \author Xiaoqiang Yang
  !>        \date Apr 2022
  
  function ncfile_getDim(fname,var_name)


    implicit none
    !
    character(len=*),      intent(in)  :: fname
    character(len=*),      intent(in)  :: var_name
    integer(i4), dimension(5)          :: ncfile_getDim
    !

    integer(i4)                           :: ncid    ! id of input stream
    integer(i4)                           :: varid   ! id of variable to be read
    character(NF90_MAX_NAME)              :: dname    ! name of dimension in the file
    integer(i4)                           :: NumDims ! # of dimensions
    integer(i4),dimension(:), allocatable :: DimIds
    integer(i4) :: n,ilength
    !
    ! Open NetCDF filename
    call check(nf90_open(fname, NF90_NOWRITE, ncid))
    !
    call check(nf90_inq_varid(ncid, var_name, varid))
    !
    call check(nf90_inquire_variable(ncid, varid, ndims=NumDims))

    allocate(DimIds(NumDims))
    call check(nf90_inquire_variable(ncid, varid, dimids=DimIds))
    if ( NumDims > size(ncfile_getDim) ) &
       stop 'ERROR*** Dimension size of Variable is greater than 5. Are you sure?'
    ncfile_getDim(:) = -1 
    do n =1, NumDims
       call check(nf90_inquire_dimension(ncid, DimIds(n), name=dname, len=ilength))
       ncfile_getDim(n) = ilength
    end do
    ! close File
    call check(nf90_close(ncid))

  end function ncfile_getDim

!-----------------------------------------------------------------
!read record values/different variables
  !------------------------------------------------------------------
  !     NAME
  !         ncfile_read_3d
  !
  !     HISTORY
  !>        \author Xiaoqiang Yang
  !>        \date Apr 2022
  
  subroutine ncfile_read_3d(fname, var_name, readinData, &
              nstart, ncount)

    implicit none
    character(len=*),                           intent(in)    :: fname
    character(len=*),                           intent(in)    :: var_name ! Variable name
    real(dp),    dimension(:,:,:),              intent(inout) :: readinData    ! array where values should be stored
    integer(i4), dimension(:),                  intent(in) :: nstart
    integer(i4), dimension(:),                  intent(in) :: ncount
    !
    integer(i4) :: ncid, varid
    !
    call check(nf90_open(fname,NF90_NOWRITE, ncid))
    !get the var id
    call check(nf90_inq_varid(ncid, var_name, varid))
    !get the values
    call check(nf90_get_var(ncid, varid, readinData, start = nstart, count = ncount))
    !
    call check(nf90_close(ncid))

  end subroutine ncfile_read_3d

  !------------------------------------------------------------------
  !     NAME
  !         ncfile_read_2d
  !
  !     HISTORY
  !>        \author Xiaoqiang Yang
  !>        \date Apr 2022
  
  subroutine ncfile_read_2d(fname, var_name, readinData, &
              nstart, ncount)

    implicit none
    character(len=*),                           intent(in)    :: fname
    character(len=*),                           intent(in)    :: var_name ! Variable name
    real(dp),    dimension(:,:),                intent(inout) :: readinData    ! array where values should be stored
    integer(i4), dimension(:),                  intent(in) :: nstart
    integer(i4), dimension(:),                  intent(in) :: ncount
    !
    integer(i4) :: ncid, varid
    !
    call check(nf90_open(fname,NF90_NOWRITE, ncid))
    !get the var id
    call check(nf90_inq_varid(ncid, var_name, varid))
    !get the values
    call check(nf90_get_var(ncid, varid, readinData, start = nstart, count = ncount))
    !
    call check(nf90_close(ncid))

  end subroutine ncfile_read_2d


!--------------------------------------------------------------------------------
 subroutine check(status)
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, trim(nf90_strerror(status))
      stop "Stopped"
    end if
 end subroutine check  
!---------------------------------------------------------------------------------

end module mo_wqm_netcdfinoutputs

