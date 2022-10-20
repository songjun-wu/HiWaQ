!> \file mo_spatial_organize.f90

!> \brief Reads water quality configure, input files and allocating WQ global variables.

!> \details This module is for water quality modeling.\n
!> 

!> \authors Xiaoqiang Yang
!> \date Sep 2021

MODULE mo_spatial_organize

  use mo_kind, only: i4, dp
  implicit none
  PUBLIC :: set_topologic_network
  PUBLIC :: get_connection_sequences
  PUBLIC :: match_RtoM

CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         set_topologic_network

  !     PURPOSE
  !>        \brief Set network topology

  !>        \details Set network topology from and to node for all links
  !>                 at level-11 (\ref fig_routing "Routing Network"). \n

  !>                 If a variable is added or removed here, then it also has to 
  !>                 be added or removed in the subroutine L11_config_set in
  !>                 module mo_restart and in the subroutine set_L11_config in module
  !>                 mo_set_netcdf_restart.

  !     LITERATURE
  !         None
  
  !     HISTORY
  !> \authors Xiaoqiang Yang
  !> \date Sep 2021   
  subroutine set_topologic_network(nty)
    use mo_wqm_constants, only: nodata_i4
    use mo_append, only: append
    use mo_wqm_global_variables, only: &
         r_Level, Rout_mask,     &
         L1_nLinks,               &
         L1_fDir,                 &
         L1_fromN,                & ! INOUT: from node 
         L1_toN, OutNode                     ! INOUT: to node

    implicit none
    integer(i4),  intent(in)    :: nty  !identify the flow direction numbering system arcgis or PCRaster map type
    ! local
    !integer(i4)                               :: nNodes
    integer(i4)                               :: nrows, ncols
    logical,     dimension(:,:), allocatable  :: mask_loc, mask_loc2
    integer(i4), dimension(:,:), allocatable  :: tmp_fdir,cellid   
    integer(i4), dimension(:,:), allocatable  :: fDir1
    integer(i4), dimension(:,:), allocatable  :: cellCoor
    integer(i4)                               :: jj, kk, ic, jc ,icc,jcc
    integer(i4)                               :: fn, tn

    integer(i4), dimension(:), allocatable    :: nLinkFromN, nLinkToN 


    nrows = r_Level%nrows
    ncols = r_Level%ncols
    !take care of the Rout_mask (ncols,nrows)
    !and L1_fDir is also ordered as (ncols,nrows) and then packed as one dimension
    allocate ( mask_loc(ncols,nrows) )
    allocate ( tmp_fdir(ncols,nrows))
    mask_loc(:,:) = .FALSE.
    mask_loc(:,:) = Rout_mask
    tmp_fdir(:,:) =  UNPACK( L1_fDir,  mask_loc, nodata_i4 )
    !     Routing network vectors have nNodes size instead of nLinks to
    !     avoid the need of having two extra indices to identify a basin. 
    ! allocate
    allocate ( nLinkFromN ( L1_nLinks    ) )  ! valid from (1 : nLinks)
    allocate ( nLinkToN   ( L1_nLinks    ) )  ! "
    allocate ( fDir1     ( nrows, ncols ) )  ! the noraml geo-type dimension
    allocate ( mask_loc2 ( nrows, ncols ) )  ! the noraml geo-type dimension
    allocate ( cellid ( nrows, ncols ) ) ! the noraml geo-type dimension
    ! initialize
    nLinkFromN(:)   = nodata_i4
    nLinkToN(:)     = nodata_i4
    fDir1(:,:)     = nodata_i4
    fDir1(:,:) =  transpose(tmp_fdir)

    ! ------------------------------------------------------------------
    !  network topology
    ! ------------------------------------------------------------------
    ! allocate
    allocate ( cellCoor(L1_nLinks,2) )

    ! initialize
    cellCoor(:,:) = nodata_i4
    !now mask_loc2(nrows, ncols) which is needed for the actual location and flow direction
    mask_loc2(:,:) = transpose(mask_loc)  
    ! counting valid cells and assign the row and col coordinations to each cell id
    kk = 0
      do icc = 1, nrows
       do jcc = 1, ncols
          if ( mask_loc2(icc,jcc) ) then
             kk = kk+1
             cellid(icc,jcc) = kk
             cellCoor(kk,1) = icc
             cellCoor(kk,2) = jcc

          end if
       end do
      end do

    jj = 0
    do kk = 1, L1_nLinks
       ic = cellCoor(kk,1) !row
       jc = cellCoor(kk,2) !col
       fn = kk

       call moveDownOneCell(fDir1(ic, jc), ic, jc, nty)
       tn = cellid(ic,jc)
       ! commented by yangx 2021-10 
       ! the case of fn == tn identifies the outlet cell,
       if (fn == tn) then
           OutNode = kk 
       else
           jj = jj + 1
           nLinkFromN(jj) = fn
           nLinkToN(jj)   = tn
       end if
    end do

    !--------------------------------------------------------
    ! Start padding up local variables to global variables
    !--------------------------------------------------------

    ! L11 data sets
    call append( L1_fromN, nLinkFromN(:) ) ! sinks are at the end 
    call append( L1_toN,   nLinkToN(:)   )
    ! free space
    deallocate (mask_loc, mask_loc2,cellid,cellCoor, tmp_fdir, fDir1, nLinkFromN, nLinkToN)   

  end subroutine set_topologic_network

  ! ------------------------------------------------------------------

  !     NAME
  !         get_connection_sequences

  !     PURPOSE
  !>        \brief Find routing order, headwater cells and sink

  !>        \details Find routing order, headwater cells and sink. \n
  !>                 If a variable is added or removed here, then it also has to 
  !>                 be added or removed in the subroutine L11_config_set in
  !>                 module mo_restart and in the subroutine set_L11_config in module
  !>                 mo_set_netcdf_restart

  !     INTENT(IN)
  !>        \param[in] "integer(i4)        :: iBasin"             Basin Id         

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
  !         None

  !     RESTRICTIONS
  !         None
  
  !     EXAMPLE
  !         None
  
  !     LITERATURE
  !         None
  
  !     HISTORY
  !>         \authors Xiaoqiang Yang
  !>         \date Sep 2021
  ! ------------------------------------------------------------------

  subroutine get_connection_sequences
    use mo_wqm_constants, only: nodata_i4
    use mo_append, only: append
    use mo_wqm_global_variables, only: &
         L1_nLinks,             &
         L1_fromN,                & ! INOUT: from node 
         L1_toN,                  & ! INOUT: to node
         L1_rOrder,               &
         L1_netPerm,OutNode


    implicit none

    ! local
    !integer(i4)                               :: nNodes
    integer(i4)                               :: nLinks

    integer(i4), dimension(:), allocatable    :: nLinkFromN      ! from node
    integer(i4), dimension(:), allocatable    :: nLinkToN        ! to node
    integer(i4), dimension(:), allocatable    :: nLinkROrder     ! network routing order  
    integer(i4), dimension(:), allocatable    :: netPerm         ! routing order (permutation)
    integer(i4)                               :: ii, jj, kk
    logical                                   :: flag

   
    nLinks  = L1_nLinks - 1_i4
    !  Routing network vectors have nNodes size instead of nLinks to
    !  avoid the need of having two extra indices to identify a basin. 

    ! allocate
    allocate ( nLinkFromN  ( L1_nLinks ) )  ! all vectors valid from (1 : nLinks)
    allocate ( nLinkToN    ( L1_nLinks ) )
    allocate ( nLinkROrder ( L1_nLinks) )
    allocate ( netPerm     ( L1_nLinks ) )
    ! initialize
    nLinkFromN(:)         = nodata_i4
    nLinkToN(:)           = nodata_i4
    nLinkROrder(1:nLinks) = 1
    nLinkROrder(L1_nLinks)   = nodata_i4
    netPerm(:)            = nodata_i4

    ! for a single node model run
    if(L1_nLinks .GT. 1) then
      ! get network vectors
      nLinkFromN(:) = L1_fromN (:)
      nLinkToN(:)   = L1_toN   (:)

      loop1: do ii = 1, nLinks
         loop2: do jj = 1, nLinks
            if ( jj == ii ) cycle loop2
            if ( nLinkFromN(ii) == nLinkToN(jj) ) then
               nLinkROrder(ii) = -9
            end if
            if ( nLinkROrder(ii) == -9 ) cycle loop1
         end do loop2
      end do loop1

      ! counting headwaters
      kk = 0
      do ii = 1, nLinks
         if ( nLinkROrder(ii) == 1) then
            kk = kk + 1
            nLinkROrder(ii) = kk
         end if
      end do

      ! counting downstream
      do while ( minval( nLinkROrder( 1 : nLinks ) ) < 0 )
       !!  print *, count(nLinkROrder .lt. 0), minval(nLinkROrder)
         loop3: do ii = 1, nLinks
            if ( .NOT. nLinkROrder(ii) == -9 ) cycle loop3
            flag = .TRUE.
            loop4: do jj = 1, nLinks
               if ( jj == ii .OR. nLinkFromN(ii)  /=  nLinkToN(jj) ) then
                  cycle loop4
               else if (.NOT. (  nLinkFromN(ii)  == nLinkToN(jj)  .AND. nLinkROrder(jj) > 0 )) then
                  flag = .FALSE.
                  exit loop4
               else
               end if
            end do loop4
            if (flag) then
               kk = kk + 1
               nLinkROrder(ii) = kk
            end if
         end do loop3
      end do


      ! keep routing order
      do ii = 1, nLinks
         netPerm(nLinkROrder(ii)) = ii
      end do
 
      !the outlet grid
      !L1_nCells = nLinks +1
      nLinkROrder(L1_nLinks)   = L1_nLinks
      netPerm(nLinkROrder(L1_nLinks)) = L1_nLinks
      L1_fromN(L1_nLinks) = OutNode
      L1_toN(L1_nLinks) = OutNode
 
      ! end of multi-node network design loop
    end if
   
    !--------------------------------------------------------
    ! Start padding up local variables to global variables
    !--------------------------------------------------------
    ! L11 network data sets 
    call append( L1_rOrder,  nLinkROrder(:) )
    call append( L1_netPerm, netPerm(:)     )

    ! free space
    deallocate (nLinkFromN, nLinkToN, nLinkROrder, netPerm )   

  end subroutine get_connection_sequences


  ! ------------------------------------------------------------------

  !     NAME
  !         match_RtoM

  !     PURPOSE
  !>        \brief matching the grid cell ids at modeling and routing resolutions

  !>        \details 

  !     LITERATURE
  !         mHM/mRM implementations
  
  !     HISTORY
  !> \authors Xiaoqiang Yang
  !> \date Sep 2021   
  subroutine match_RtoM()
    use mo_wqm_constants, only: nodata_i4,nodata_dp
    use mo_append, only: append
    use mo_wqm_global_variables, only: &
         m_Level, Model_mask, RIdonM,&
         r_Level, Rout_mask, L1_area,R_areaCell,    &
         L1_nLinks

    implicit none
    ! local
    !integer(i4)                               :: nNodes
    integer(i4)                               :: Mnrows, Mncols,Rnrows, Rncols
    logical,     dimension(:,:), allocatable  :: mask_loc
    real(dp),    dimension(:,:), allocatable  :: cellarea
    integer(i4), dimension(:,:), allocatable  :: Rcellid,RId_on_M 
    integer(i4)                               ::  kk,icc,jcc
    integer(i4)                              :: il,ir, ju,jd
    real(dp)                                 :: cellFactorR
    !real(dp), dimension(:), allocatable      :: cellareaR


    Rnrows = r_Level%nrows
    Rncols = r_Level%ncols
    Mnrows = m_Level%nrows
    Mncols = m_Level%ncols

    cellFactorR = r_Level%cellsize/m_Level%cellsize
      !take care of the Rout_mask (ncols,nrows)
      !and L1_fDir is also ordered as (ncols,nrows) and then packed as one dimension
      allocate ( mask_loc(Rncols,Rnrows) )
      mask_loc(:,:) = .FALSE.
      mask_loc(:,:) = Rout_mask

      ! allocate

      allocate ( Rcellid(Rncols,Rnrows  ) ) !
      allocate ( RId_on_M(Mncols,Mnrows  ) )
      allocate ( cellarea(Mncols,Mnrows  ))
      allocate ( R_areaCell(L1_nLinks))
      ! initialize
      RId_on_M(:,:)   = nodata_i4
      Rcellid(:,:)    = nodata_i4
      cellarea(:,:)   = nodata_dp
      R_areaCell(:)    = nodata_dp

      cellarea(:,:) = unpack(L1_area, Model_mask, nodata_dp)
      ! counting valid cells and assign the row and col coordinations to each cell id
      kk = 0
      do jcc = 1, Rnrows
         do icc = 1, Rncols
          if ( .not. mask_loc(icc,jcc) ) cycle
          kk = kk+1
          Rcellid(icc,jcc) = kk
          !coord. of four corners 
          il = (icc-1) * nint(cellFactorR,i4) + 1
          ir =     icc * nint(cellFactorR,i4)
          ju = (jcc-1) * nint(cellFactorR,i4) + 1
          jd =     jcc * nint(cellFactorR,i4)
          ! constrain the range of up, down, left, and right boundaries
          if( il < 1   ) il =  1
          if( ir > Mncols ) ir =  Mncols
          if( ju < 1   ) ju =  1
          if( jd > Mnrows ) jd =  Mnrows
          ! Delimitation of level-11 cells on level-1 for L11 resolution lower than L1 resolution
          RId_on_M(il:ir, ju:jd) = Rcellid(icc,jcc)
          ! effective area [m2] & total no. of model cells within a given routing cell
          R_areaCell(kk) = sum( cellarea(il:ir, ju:jd), Model_mask(il:ir, ju:jd) )
          !print*, kk,il,ir, ju,jd
         end do
      end do

      ! Start padding up local variables to global variables
      !call append( R_areaCell, cellareaR(:) ) !
      call append( RIdonM, pack(RId_on_M, Model_mask))

      ! free space
      deallocate (mask_loc,Rcellid, cellarea)   


  end subroutine match_RtoM

  !local subroutine for identifying the downstream location
  !--------------------------------------------------------------------
  subroutine moveDownOneCell(fDir, iRow, jCol, nty)
    use mo_message,             only : message

    implicit none

    integer(i4), intent(IN)     :: fDir
    integer(i4), intent(INOUT)  :: iRow, jCol
    integer(i4), intent(IN)     :: nty

   if (nty == 1) then   !arcgis type
    select case (fDir)
    case(1)   !E
       jCol = jCol + 1
    case(2)   !SE
       iRow = iRow + 1
       jCol = jCol + 1
    case(4)   !S
       iRow = iRow + 1
    case(8)   !SW
       iRow = iRow + 1
       jCol = jCol - 1
    case(16)  !W
       jCol = jCol - 1
    case(32)  !NW
       iRow = iRow - 1
       jCol = jCol - 1
    case(64)  !N
       iRow = iRow - 1
    case(128) !NE
       iRow = iRow - 1
       jCol = jCol + 1
    case default !sink
       ! do nothing
    end select
   else if (nty == 2) then   !PCRater map type
    select case (fDir)
    case(6)   !E
       jCol = jCol + 1
    case(3)   !SE
       iRow = iRow + 1
       jCol = jCol + 1
    case(2)   !S
       iRow = iRow + 1
    case(1)   !SW
       iRow = iRow + 1
       jCol = jCol - 1
    case(4)  !W
       jCol = jCol - 1
    case(7)  !NW
       iRow = iRow - 1
       jCol = jCol - 1
    case(8)  !N
       iRow = iRow - 1
    case(9) !NE
       iRow = iRow - 1
       jCol = jCol + 1
    case default !sink
       ! do nothing
    end select
   else        !the information should be provided in configuration file
    call message()
    call message('***ERROR: the numbering system of flow direction inputs should be given')
    call message('1 for ArcGIS style; 2 for PCRaster style')
    stop
   end if

  end subroutine moveDownOneCell



END MODULE mo_spatial_organize 
