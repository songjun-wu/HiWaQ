!> \file mo_wqm_parameterization.f90

!> \brief Upscale parameters to model levels

!> \details 

!>        \author Xiaoqiang Yang 
!>        \date Sep 2021

MODULE mo_wqm_parameterization

 
  USE mo_kind, ONLY: i4, sp, dp


  IMPLICIT NONE

 
  PUBLIC :: wqm_parameterization  ! parameters of water quality model


CONTAINS

  ! ------------------------------------------------------------------

  !     NAME
  !         wqm_parameterization

  !     PURPOSE
  !>        \brief Water quality parameters 
  !     LITERATURE
  !         None

  !     HISTORY
  !>        \author Xiaoqiang Yang 
  !>        \date Sep 2021

  subroutine wqm_parameterization(ncells, nNodes,nlandusetotal,nsoiltypetotal, fLAI, fsoil, degradN_rate, &
                mineraN_rate, dissolN_rate, sdenitr_rate, adenitr_rate,atruptk_rate, priprod_rate)

  use mo_message,                 only: message
  use mo_wqm_global_variables,    only: nTotalParam, wqm_param,model_type,resolution, rout_resolution, &
        RIdonM,L1_area, R_areaCell,L1_nHRUs

  implicit none
  integer(i4),                       intent(in)     :: ncells      ! number of cells at model level
  integer(i4),                       intent(in)     :: nNodes       ! number of cells at routing level
  integer(i4),                       intent(in)     :: nlandusetotal       ! number of land use types
  integer(i4),                       intent(in)     :: nsoiltypetotal       ! number of soil types
  real(dp), dimension(:,:),          intent(in)     :: fLAI         ! area fraction of each landuse type
  real(dp), dimension(:,:),          intent(in)     :: fsoil        ! area fraction of each landuse type
  real(dp), dimension(:),            intent(out)    :: degradN_rate ! parameter 
  real(dp), dimension(:),            intent(out)    :: mineraN_rate ! parameter   
  real(dp), dimension(:),            intent(out)    :: dissolN_rate ! parameter   
  real(dp), dimension(:),            intent(out)    :: sdenitr_rate ! parameter    
  real(dp), dimension(:),            intent(out)    :: adenitr_rate ! parameter  
  real(dp), dimension(:),            intent(out)    :: atruptk_rate ! parameter   
  real(dp), dimension(:),            intent(out)    :: priprod_rate ! parameter



  !local
  integer(i4)    :: ip, igrp,iids,sloc,eloc,nn,k
  real(dp), dimension(:),  allocatable :: tmp_pvalue     
  real(dp), dimension(:,:),  allocatable ::fLAI_11, fsoil_11  
  
    !loop over no of groups for each parameter --terrestrial grids
    do ip=4,nTotalParam
      allocate(tmp_pvalue(ncells))
      tmp_pvalue = 0.0_dp
      do igrp =1, wqm_param(ip)%ngroups
        if (wqm_param(ip)%category == 1_i4) then !land use dependent
          do iids = 1, nlandusetotal
            if (any(wqm_param(ip)%groupids(igrp,:) .eq. iids)) then
              tmp_pvalue(:) = tmp_pvalue(:) + wqm_param(ip)%pvalue(igrp) * fLAI(:,iids)
            end if
          end do
        else if (wqm_param(ip)%category == 2_i4) then !soil type dependent
          do iids = 1, nsoiltypetotal
            if (any(wqm_param(ip)%groupids(igrp,:) .eq. iids)) then
              tmp_pvalue(:) = tmp_pvalue(:) + wqm_param(ip)%pvalue(igrp) * fsoil(:,iids)
            end if
          end do
        end if
      end do

     !assign parameter values to global variable
      select case (ip) !the same order as in wqm_parameter.nml
      case(4)
        sdenitr_rate(:) = tmp_pvalue(:)
      case(5)
        mineraN_rate(:) = tmp_pvalue(:)
      case(6)
        degradN_rate(:) = tmp_pvalue(:)
      case(7)
        dissolN_rate(:) = tmp_pvalue(:)
      case default
        call message()
        call message('***ERROR: something wrong with N module parameterization')
        call message('please modify mo_wqm_parameterization.f90, if you want to change the parameterization scheme')
        stop
      end select
      deallocate(tmp_pvalue)
    end do
    !parameters at the routing level
    if (.not.allocated(fLAI_11)) allocate(fLAI_11(nNodes,nlandusetotal ))
    if (.not.allocated(fsoil_11)) allocate(fsoil_11(nNodes, nsoiltypetotal))
    !
    if (model_type ==1) then
      if (abs(rout_resolution - resolution) > 0.0001_dp) then
        do k=1,nNodes
          do nn = 1, nlandusetotal
            fLAI_11(k,nn) =  sum(fLAI(:,nn) * L1_area(:) / R_areaCell(k) , RIdonM(:) .eq. k )
          end do
          do nn =1, nsoiltypetotal
            fsoil_11(k,nn) = sum(fsoil(:,nn) * L1_area(:) / R_areaCell(k) , RIdonM(:) .eq. k )
          end do
        end do
      else
        fLAI_11 = fLAI
        fsoil_11 = fsoil
      end if
    else !semi-distributed structure -- nLinks == nSubs
        do k=1,nNodes
          sloc = L1_nHRUs*(k-1)+1
          eloc = L1_nHRUs*k  
          do nn = 1, nlandusetotal
            fLAI_11(k,nn) =  sum(fLAI(sloc:eloc,nn) * L1_area(sloc:eloc))/ R_areaCell(k)
          end do
          do nn =1, nsoiltypetotal
            fsoil_11(k,nn) = sum(fsoil(sloc:eloc,nn) * L1_area(sloc:eloc)) / R_areaCell(k)
          end do
        end do
    end if
    !
    do ip=1,3
      allocate(tmp_pvalue(nNodes))
      tmp_pvalue = 0.0_dp
      do igrp =1, wqm_param(ip)%ngroups
        if (wqm_param(ip)%category == 1_i4) then !land use dependent
          do iids = 1, nlandusetotal
            if (any(wqm_param(ip)%groupids(igrp,:) .eq. iids)) then
              tmp_pvalue(:) = tmp_pvalue(:) + wqm_param(ip)%pvalue(igrp) * fLAI_11(:,iids)
            end if
          end do
        else if (wqm_param(ip)%category == 2_i4) then !soil type dependent
          do iids = 1, nsoiltypetotal
            if (any(wqm_param(ip)%groupids(igrp,:) .eq. iids)) then
              tmp_pvalue(:) = tmp_pvalue(:) + wqm_param(ip)%pvalue(igrp) * fsoil_11(:,iids)
            end if
          end do
        end if
      end do
	  
      select case (ip) !the same order as in wqm_parameter.nml
      case(1)
        adenitr_rate(:) = tmp_pvalue(:)
      case(2)
        atruptk_rate(:) = tmp_pvalue(:)
      case(3)
        priprod_rate(:) = tmp_pvalue(:)
      case default
        call message()
        call message('***ERROR: something wrong with N module parameterization')
        call message('please modify mo_wqm_parameterization.f90, if you want to change the parameterization scheme')
        stop
      end select
      deallocate(tmp_pvalue)
    end do

  end subroutine wqm_parameterization
END MODULE mo_wqm_parameterization
