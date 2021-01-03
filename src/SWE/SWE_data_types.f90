! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE)
MODULE SWE_data_types
#if defined(_SWE_DG)
  use SWE_DG_matrices
   
  implicit none
  
  integer, PARAMETER :: BYTE = selected_int_kind(1)
  integer, PARAMETER :: SHORT = selected_int_kind(4)
  integer, PARAMETER :: GRID_SI = selected_int_kind(8)
  integer, PARAMETER :: GRID_DI = selected_int_kind(16)
  
  integer, PARAMETER :: SR = GRID_SR
  integer, PARAMETER :: SI = GRID_SI
  integer, PARAMETER :: DI = GRID_DI
  
  real (kind = GRID_SR), parameter					:: g = _GRAV_CONSTANT		!< gravitational constant
  
  
  !***********************
  !Entity data
  !***********************
  
  !> state vector of DoFs, either as absoulte values or updates
  
  type t_dof_state
     real (kind = GRID_SR) :: h			!< water change
     real (kind = GRID_SR), dimension(2) :: p   !< momentum change
     
   contains
     
     procedure, pass :: add => dof_state_add
     procedure, pass :: inv => dof_state_inv
     procedure, pass :: scale => dof_state_scale
     
     generic :: operator(+) => add
     generic :: operator(-) => inv
     generic :: operator(*) => scale
  end type t_dof_state
  
  !> cell state vector including bathymetry
  type, extends(t_dof_state) :: t_state
     real (kind = GRID_SR) :: b						!< bathymetry 
   contains
     
     procedure, pass :: add_state => state_add
     generic :: operator(+) => add_state
  end type t_state
  
  !> update vector
  type, extends(t_dof_state) :: t_update
   contains
     procedure, pass :: add_update => update_add
     generic :: operator(+) => add_update
  end type t_update
  
  !> persistent scenario data on a node
  type num_node_data_pers
    integer (kind = BYTE)  :: dummy !< no data
  END type num_node_data_pers

		!> persistent scenario data on an edge
  type num_edge_data_pers
     integer (kind = BYTE), dimension(0) :: dummy !< no data
     ! real(kind=GRID_SR)   ,DIMENSION(_SWE_DG_ORDER+1,4)   :: QP
     ! real(kind=GRID_SR)   ,DIMENSION(2,_SWE_DG_ORDER+1,3) :: FP
     !real(kind = GRID_SR),DIMENSION(_SWE_PATCH_ORDER)    :: H, HU, HV, B
  end type num_edge_data_pers

!> persistent scenario data on a cell
  type num_cell_data_pers
  real (kind = GRID_SR), DIMENSION(_SWE_PATCH_ORDER_SQUARE):: H, HU, HV, B
  type(t_state)        , DIMENSION(_SWE_DG_DOFS)           :: Q   !DG degrees of freedom
  real(kind=GRID_SR)   , DIMENSION(_SWE_DG_DOFS,3)         :: Q_DG_UPDATE !Predictor element update
  real(kind=GRID_SR)   , DIMENSION(3,  _SWE_DG_ORDER+1,4)  :: QP  !Predictor projections on edges
  real(kind=GRID_SR)   , DIMENSION(3,2,_SWE_DG_ORDER+1,3)  :: FP  !Predictor projections on edges
  integer :: troubled

#if defined(_CELL_METRICS)
  real(kind=GRID_SR)   , DIMENSION(_SWE_DG_DOFS)            :: dt    = 0.0
  real(kind=GRID_SR)   , DIMENSION(_SWE_PATCH_ORDER_SQUARE) :: dt_fv = 0.0
  real(kind=GRID_SR)   , DIMENSION(_SWE_DG_DOFS,3)          :: QP_avg = 0.0
#endif
  
#if defined(_DEBUG)
  integer :: debug_flag = 0
#endif                        
  contains
    procedure :: get_dofs_dg => get_dofs_dg
    procedure :: set_dofs_dg => set_dofs_dg
end type num_cell_data_pers

  !> Cell representation on an edge, this would typically be everything required from a cell to compute the flux function on an edge
type num_cell_rep
    !type(t_state), DIMENSION((_SWE_DG_ORDER+1)*3)           :: Q
    real(kind=GRID_SR), DIMENSION(_SWE_DG_ORDER+1,4)        :: QP
    real(kind=GRID_SR), DIMENSION(2,_SWE_DG_ORDER+1,3)      :: FP
    real (kind = GRID_SR), dimension (_SWE_PATCH_ORDER)     :: H, HU, HV, B
    real (kind = GRID_SR), dimension (_DMP_NUM_OBSERVABLES) :: minObservables
    real (kind = GRID_SR), dimension (_DMP_NUM_OBSERVABLES) :: maxObservables
    integer :: troubled
#if defined(_DEBUG)   
    integer :: debug_flag = 0.0_GRID_SR
#endif   
  end type num_cell_rep

  !> Cell update, this would typically be a flux function
  type num_cell_update
    type(t_update), DIMENSION(_SWE_DG_ORDER+1)   :: flux
    real (kind = GRID_SR), DIMENSION(_SWE_PATCH_ORDER) :: H, HU, HV, B !< values of ghost cells
    real (kind = GRID_SR), dimension (_DMP_NUM_OBSERVABLES) :: minObservables
    real (kind = GRID_SR), dimension (_DMP_NUM_OBSERVABLES) :: maxObservables
    integer :: troubled
  end type num_cell_update

		!*************************
		!Temporary per-Entity data
		!*************************

		!> temporary scenario data on a node (deleted after each traversal)
		type num_node_data_temp
			integer (kind = BYTE), dimension(0)										:: dummy					!< no data
		END type num_node_data_temp

		!> temporary scenario data on an edge (deleted after each traversal)
		type num_edge_data_temp
			integer (kind = BYTE), dimension(0)										:: dummy					!< no data
		END type num_edge_data_temp

		!> temporary scenario data on a cell (deleted after each traversal)
		type num_cell_data_temp
		END type num_cell_data_temp

		!***********************
		!Global data
		!***********************

		!> Data type for the scenario configuration
  type num_global_data
     real (kind = GRID_SR)							:: r_time					!< simulation time
     real (kind = GRID_SR)							:: r_dt						!< time step
     real (kind = GRID_SR)							:: r_dt_new					!< new time step for the next iteration
     real (kind = GRID_SR)              :: min_courant
     
     !For Refinment
     real (kind = GRID_SR) :: min_error
     real (kind = GRID_SR) :: max_error
     real (kind = GRID_SR) :: min_error_new
     real (kind = GRID_SR) :: max_error_new

     
#if defined(_ASAGI)
     real (kind = GRID_SR)							:: b_max=TINY(1.0_GRID_SR),b_min=HUGE(1.0_GRID_SR)
     real (kind = GRID_SR)							:: b_max_new=TINY(1.0_GRID_SR),b_min_new=HUGE(1.0_GRID_SR)     

#endif
#     if defined(_XDMF)
        integer (kind = GRID_SI)                               :: xdmf_filter_count_cells !<amount of cells in this section after filtering
#       if defined(_SWE_PATCH)
          integer (kind = GRID_SI)                             :: xdmf_filter_count_patches !<amount of patches in this section after filtering
#       endif
#     endif   
		end type

		contains

    subroutine apply_phi(dg,fv)
      real(kind=GRID_SR),intent(out) :: fv(_SWE_PATCH_ORDER_SQUARE)
      real(kind=GRID_SR),intent(in)  :: dg(_SWE_DG_DOFS)
      real(kind=GRID_SR)             :: int_dg,int_fv
      fv=matmul(phi_hat,dg)

      int_dg = dot_product(weights,dg)
      int_fv = sum(fv) * _REF_TRIANGLE_SIZE
      if(abs(int_fv) > 0.0_GRID_SR) then
         fv = fv * int_dg/int_fv
      end if

    end subroutine apply_phi
    
    subroutine apply_phi_cons(h_dg,hu_dg,hv_dg,h_fv,hu_fv,hv_fv)
      real(kind=GRID_SR),dimension(_SWE_DG_DOFS),intent(in):: h_dg,hu_dg,hv_dg
      real(kind=GRID_SR),dimension(_SWE_PATCH_ORDER_SQUARE),intent(out) :: h_fv,hu_fv,hv_fv

      real(kind=GRID_SR),dimension(_SWE_DG_DOFS) :: u_dg,v_dg
      real(kind=GRID_SR),dimension(_SWE_PATCH_ORDER_SQUARE) :: u_fv,v_fv,h_fv_orig
      
      real(kind = GRID_SR) :: negative_mass,diff_mass
      integer :: num_cells_positive_mass,i,j
      integer :: surface_indeces(_SWE_PATCH_ORDER_SQUARE)
      logical :: surface_mask(_SWE_PATCH_ORDER_SQUARE)
      u_dg = 0.0_GRID_SR
      v_dg = 0.0_GRID_SR
      
      call apply_phi(h_dg,h_fv)
      call apply_phi(hv_dg,hv_fv)
      call apply_phi(hu_dg,hu_fv)
      h_fv_orig = h_fv
      
      num_cells_positive_mass = 0
      do i = 1,_SWE_PATCH_ORDER_SQUARE
         if(h_fv(i) .le. 0.0_GRID_SR) then
            negative_mass = negative_mass - h_fv(i)
            h_fv(i) = 0.0_GRID_SR
         else
            num_cells_positive_mass =  num_cells_positive_mass + 1
         end if
      end do
      
      surface_mask = .True.
      do i = 1,_SWE_PATCH_ORDER_SQUARE
         surface_indeces(i) =  minloc(h_fv,1,surface_mask)
         surface_mask(surface_indeces(i)) =  .False.
      end do
                       
      if (negative_mass > 0.0_GRID_SR) then
         do j = 1,_SWE_PATCH_ORDER_SQUARE
            i = surface_indeces(j)
            if(h_fv(i) > 0.0_GRID_SR) then               
               diff_mass = min(negative_mass/num_cells_positive_mass,h_fv(i))
               
               h_fv(i) = h_fv(i) - diff_mass               
               negative_mass = negative_mass - diff_mass               
               num_cells_positive_mass = num_cells_positive_mass - 1
            end if
         end do
      endif
      
      if(negative_mass > 0.0_GRID_SR) then
         print*,negative_mass
         print*,num_cells_positive_mass
      endif

      where(h_fv_orig > 0.0_GRID_SR)
         hu_fv = hu_fv * (h_fv/h_fv_orig)
         hv_fv = hv_fv * (h_fv/h_fv_orig)
      end where
      
    end subroutine apply_phi_cons

    
    subroutine apply_mue(fv,dg)
      real(kind=GRID_SR),intent(in) :: fv(_SWE_PATCH_ORDER_SQUARE)
      real(kind=GRID_SR),intent(out)  :: dg(_SWE_DG_DOFS)
      real(kind=GRID_SR)             :: q_temp(_SWE_DG_DOFS+1)
      real(kind=GRID_SR)             :: int_dg,int_fv
      dg= matmul(mue_inv,fv)
      int_dg = dot_product(weights,dg)
      int_fv = sum(fv) * _REF_TRIANGLE_SIZE
      if(abs(int_dg) > 0.0_GRID_SR) then
         dg = dg * (int_fv/int_dg)
      end if
    end subroutine apply_mue
    
    subroutine apply_mue_sample(fv,dg)
      real(kind=GRID_SR),intent(in) :: fv(_SWE_PATCH_ORDER_SQUARE)
      real(kind=GRID_SR),intent(out)  :: dg(_SWE_DG_DOFS)
      real(kind=GRID_SR)             :: q_temp(_SWE_DG_DOFS+1)
      real(kind=GRID_SR)             :: fv_max
      real(kind=GRID_SR)             :: fv_min
      real(kind=GRID_SR)             :: int,alpha
      real(kind=GRID_SR)             :: int_dg,int_fv
      !      dg= matmul(sample_fv,fv)

      dg= matmul(mue_inv_slope,fv)
      int_dg = dot_product(weights,dg)
      int_fv = sum(fv) * _REF_TRIANGLE_SIZE
      if(abs(int_dg) > 0.0_GRID_SR) then
         dg = dg * (int_fv/int_dg)
      end if

      ! dg= matmul(mue_inv,fv)
      ! fv_max = maxval(fv)
      ! fv_min = minval(fv)
      ! int = sum(fv) * _REF_TRIANGLE_SIZE
      ! alpha = max( (maxval(dg) - int) / (fv_max - int) , (minval(dg) - int) / (fv_min - int)) 
      ! dg = (dg - int) / alpha + int
      
    end subroutine apply_mue_sample
    
		!adds two state vectors
		elemental function state_add(Q1, Q2)	result(Q_out)
			class (t_state), intent(in)		:: Q1
			type (t_state), intent(in)		:: Q2
			type (t_state)					:: Q_out
			Q_out = t_state(Q1%h + Q2%h, Q1%p + Q2%p, Q1%b + Q2%b)
		end function

		!adds two update vectors
		elemental function update_add(f1, f2)	result(f_out)
			class (t_update), intent(in)		:: f1
			type (t_update), intent(in)		    :: f2
			type (t_update)					    :: f_out

			f_out = t_update(f1%h + f2%h, f1%p + f2%p)
                end function update_add

		!adds two dof state vectors
		elemental function dof_state_add(Q1, Q2)	result(Q_out)
			class (t_dof_state), intent(in)		:: Q1
			type (t_dof_state), intent(in)		:: Q2
			type (t_dof_state)					:: Q_out

			Q_out = t_dof_state(Q1%h + Q2%h, Q1%p + Q2%p)
		end function

		!inverts a dof state vector
		elemental function dof_state_inv(f)	result(f_out)
			class (t_dof_state), intent(in)		:: f
			type (t_dof_state)					:: f_out

			f_out = t_dof_state(-f%h, -f%p)
		end function

		!multiplies a scalar with a dof state vector
		elemental function dof_state_scale(f, s)	result(f_out)
			class (t_dof_state), intent(in)		:: f
			real (kind = GRID_SR), intent(in)		:: s
			type (t_dof_state)					:: f_out

			f_out = t_dof_state(s * f%h, s * f%p)
		end function

		!multiplies a scalar with a dof state vector
                subroutine get_dofs_dg(f,q)
                  integer ::i
                  class (num_cell_data_pers), intent(in)		:: f
                  real (kind = GRID_SR),intent(out)               :: q (size(f%Q,1),3)
                  q(:,1)= f%Q(:)%h
                  q(:,2)= f%Q(:)%p(1)
                  q(:,3)= f%Q(:)%p(2)
                end subroutine get_dofs_dg

		subroutine set_dofs_dg(f,q)
                  class (num_cell_data_pers),intent(inout) 	:: f
                  real (kind = GRID_SR) 		        :: q(size(f%Q,1),3)
                  f%Q(:)%H    = q(:,1)
                  f%Q(:)%p(1) =q(:,2)
                  f%Q(:)%p(2) =q(:,3)
                end subroutine set_dofs_dg
#endif
	END MODULE SWE_data_types
#endif
