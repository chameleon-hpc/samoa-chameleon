! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE)
	MODULE SWE_data_types
#if defined(_SWE_DG)
   use SWE_DG_matrices
   
		implicit none
# else
		implicit none
		PUBLIC

		!data precision
#       if defined(_SINGLE_PRECISION)
            integer, PARAMETER :: GRID_SR = kind(1.0e0)
#       elif defined(_DOUBLE_PRECISION)
            integer, PARAMETER :: GRID_SR = kind(1.0d0)
#       elif defined(_QUAD_PRECISION)
            integer, PARAMETER :: GRID_SR = kind(1.0q0)
#       else
#           error "No floating point precision is chosen!"
#       endif
#endif

		integer, PARAMETER :: BYTE = selected_int_kind(1)
		integer, PARAMETER :: SHORT = selected_int_kind(4)
		integer, PARAMETER :: GRID_SI = selected_int_kind(8)
		integer, PARAMETER :: GRID_DI = selected_int_kind(16)

                integer, PARAMETER :: SR = GRID_SR
                integer, PARAMETER :: SI = GRID_SI
                integer, PARAMETER :: DI = GRID_DI
  
		real (kind = GRID_SR), parameter					:: g = 9.80665_GRID_SR		!< gravitational constant


		!***********************
		!Entity data
		!***********************

  !> state vector of DoFs, either as absoulte values or updates
  
		type t_dof_state
			real (kind = GRID_SR)													:: h						!< water change
			real (kind = GRID_SR), dimension(2)										:: p						!< momentum change

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
			integer (kind = BYTE)															:: dummy					!< no data
		END type num_node_data_pers

		!> persistent scenario data on an edge
		type num_edge_data_pers
			integer (kind = BYTE), dimension(0)											:: dummy					!< no data
                END type num_edge_data_pers

!> persistent scenario data on a cell
type num_cell_data_pers
  real (kind = GRID_SR), DIMENSION(_SWE_PATCH_ORDER_SQUARE) :: H, HU, HV, B
  type(t_state)     , DIMENSION(_SWE_DG_DOFS)     :: Q_DG
  real(kind=GRID_SR), DIMENSION(_SWE_DG_DOFS,3)   :: Q_DG_UPDATE
  real(kind=GRID_SR), DIMENSION(2,_SWE_DG_DOFS,3) :: FP
  real(kind=GRID_SR), DIMENSION(  _SWE_DG_DOFS,4) :: QP
  type(t_state)     , DIMENSION(_SWE_CELL_SIZE)	  :: Q 
  integer :: troubled
#if defined(_DEBUG)
  integer :: debug_flag = 0
#endif                        
  contains
                          procedure :: convert_fv_to_dg => convert_fv_to_dg
                          procedure :: convert_dg_to_fv => convert_dg_to_fv
                          procedure :: convert_dg_to_fv_bathymetry => convert_dg_to_fv_bathymetry
                          procedure :: convert_fv_to_dg_bathymetry => convert_fv_to_dg_bathymetry

                          procedure :: get_dofs_dg => get_dofs_dg
                          procedure :: set_dofs_dg => set_dofs_dg
end type num_cell_data_pers

  !> Cell representation on an edge, this would typically be everything required from a cell to compute the flux function on an edge
  type num_cell_rep
    type(t_state), DIMENSION((_SWE_DG_ORDER+1)*3)           :: Q
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
			integer (kind = BYTE), dimension(0)										:: dummy					!< no data
		END type num_cell_data_temp

		!***********************
		!Global data
		!***********************

		!> Data type for the scenario configuration
  type num_global_data
     real (kind = GRID_SR)							:: r_time					!< simulation time
     real (kind = GRID_SR)							:: r_dt						!< time step
     real (kind = GRID_SR)							:: r_dt_new					!< new time step for the next iteration
     
#if defined(_ASAGI)
#if defined(_SWE_DG)
     real (kind = GRID_SR)							:: b_max=TINY(1.0_GRID_SR),b_min=HUGE(1.0_GRID_SR)
     real (kind = GRID_SR)							:: b_max_new=TINY(1.0_GRID_SR),b_min_new=HUGE(1.0_GRID_SR)     

#endif
#endif
#     if defined(_XDMF)
        integer (kind = GRID_SI)                               :: xdmf_filter_count !<amount of cells in this section after filtering
#     endif   
		end type

		contains
#if defined (_SWE_DG)

                 subroutine convert_fv_to_dg(dofs)
                   class(num_cell_data_pers) :: dofs
                   real(kind=GRID_SR) :: q(_SWE_PATCH_ORDER*_SWE_PATCH_ORDER,3),q_temp(_SWE_DG_DOFS+1,3)
                   real(kind=GRID_SR) :: h_temp(_SWE_DG_DOFS +1), hu_temp(_SWE_DG_DOFS +1), hv_temp(_SWE_DG_DOFS +1)
                   real(kind=GRID_SR) :: q_dg(_SWE_DG_DOFS,3)
                   
                   real(kind=GRID_SR) :: epsilon=0.0_GRID_SR
                   integer :: i,j

!                  PATCHES store h+b DG not

                   !                   q(:,1)=dofs%H-dofs%b
                   
                   q(:,1)=dofs%H
                   q(:,2)=dofs%HU
                   q(:,3)=dofs%HV

                   q_temp(1:_SWE_DG_DOFS,:) = 2.0q0*matmul(transpose(phi),q)

                   q_temp(_SWE_DG_DOFS+1,1) = sum(q(:,1))
                   q_temp(_SWE_DG_DOFS+1,2) = sum(q(:,2))
                   q_temp(_SWE_DG_DOFS+1,3) = sum(q(:,3))

                   q_temp = q_temp /_REF_TRIANGLE_SIZE_INV

                   ! call lusolve(mue_lu,_SWE_DG_DOFS+1,mue_lu_pivot,q_temp(:,1))
                   ! call lusolve(mue_lu,_SWE_DG_DOFS+1,mue_lu_pivot,q_temp(:,2))
                   ! call lusolve(mue_lu,_SWE_DG_DOFS+1,mue_lu_pivot,q_temp(:,3))

                   q_temp(:,1)=matmul(mue_inv,q_temp(:,1))
                   q_temp(:,2)=matmul(mue_inv,q_temp(:,2))
                   q_temp(:,3)=matmul(mue_inv,q_temp(:,3))

!                   print*,q_temp(1:_SWE_DG_DOFS,:)
                   
                   q_temp(1:_SWE_DG_DOFS,1)=q_temp(1:_SWE_DG_DOFS,1) - dofs%Q_DG%b

                   q_dg=q_temp(1:_SWE_DG_DOFS,:)
                
                   
                   call dofs%set_dofs_dg(q_dg)
                   
                 end subroutine convert_fv_to_dg

                 subroutine convert_dg_to_fv(dofs)
                   class(num_cell_data_pers) :: dofs
                   real(kind=GRID_SR)        :: q(_SWE_DG_DOFS,3)
                   real(kind=GRID_SR)        :: fv_temp(_SWE_PATCH_ORDER*_SWE_PATCH_ORDER,3)

                   call dofs%get_dofs_dg(q)

                   !PATCHES store h+b DG not
                   q(:,1)  = q(:,1) + dofs%Q_DG%b

                   fv_temp=matmul(phi,q)*_REF_TRIANGLE_SIZE_INV

                   !                   dofs%H=fv_temp(:,1) + dofs%B
                   dofs%H=fv_temp(:,1)
                   dofs%HU=fv_temp(:,2)
                   dofs%HV=fv_temp(:,3)

                   ! #if defined(_ASAGI)

!                    where(0_GRID_SR)
!                       q_temp(1:_SWE_DG_DOFS,1)=0.0_GRID_SR
!                    end where
! #endif                   



                 end subroutine convert_dg_to_fv

                 subroutine apply_phi(dg,fv)
                   real(kind=GRID_SR),intent(out) :: fv(_SWE_PATCH_ORDER_SQUARE)
                   real(kind=GRID_SR),intent(in)  :: dg(_SWE_DG_DOFS)

                   fv=matmul(phi,dg)*_REF_TRIANGLE_SIZE_INV
                   
                 end subroutine apply_phi

                 subroutine apply_mue(fv,dg)
                   real(kind=GRID_SR),intent(in) :: fv(_SWE_PATCH_ORDER_SQUARE)
                   real(kind=GRID_SR),intent(out)  :: dg(_SWE_DG_DOFS)
                   real(kind=GRID_SR)             :: q_temp(_SWE_DG_DOFS+1)                   

                   q_temp(1:_SWE_DG_DOFS)= 2.0q0*matmul(transpose(phi),fv)
                   
                   q_temp(_SWE_DG_DOFS+1) = sum(fv)
                   
                   q_temp = q_temp /_REF_TRIANGLE_SIZE_INV

                  
                   !                   call lusolve(mue_lu,_SWE_DG_DOFS+1,mue_lu_pivot,q_temp)
                   q_temp=matmul(mue_inv,q_temp)                   

                   dg=q_temp(1:_SWE_DG_DOFS)
                   
                 end subroutine apply_mue



                 subroutine convert_dg_to_fv_bathymetry(dofs)
                   class(num_cell_data_pers) :: dofs
                   real(kind=GRID_SR)        :: q(_SWE_DG_DOFS)
                   real(kind=GRID_SR)        :: fv_temp(_SWE_PATCH_ORDER*_SWE_PATCH_ORDER)

                   q = dofs%Q_DG%B
                   fv_temp=matmul(phi,q)*_REF_TRIANGLE_SIZE_INV
                   dofs%B=fv_temp

                 end subroutine convert_dg_to_fv_bathymetry



                 subroutine convert_fv_to_dg_bathymetry(dofs,normals)
                   class(num_cell_data_pers) :: dofs
                   real(kind=GRID_SR) :: b_temp(_SWE_DG_DOFS +1),b_fv(_SWE_PATCH_ORDER_SQUARE)

                   real(kind=GRID_SR) :: b_x_temp(_SWE_DG_DOFS)
                   real(kind=GRID_SR) :: b_y_temp(_SWE_DG_DOFS)
                   real(kind=GRID_SR),intent(in) :: normals(2,2)
                   real(kind=GRID_SR) :: normals_normed(2,2)
                   integer ::i,j
                   real(kind=GRID_SR) :: epsilon=0.0_GRID_SR

                   normals_normed=normals/NORM2(normals(1,1:2))

                   b_temp(1:_SWE_DG_DOFS)= 2.0q0*matmul(transpose(phi),dofs%b)
                   b_temp(_SWE_DG_DOFS+1) = sum(dofs%b)
                   b_temp = b_temp /_REF_TRIANGLE_SIZE_INV
                   b_temp=matmul(mue_inv,b_temp)                   

                   do i=1,_SWE_DG_DOFS
                      dofs%Q_DG(i)%b=b_temp(i)
                   end do
                 end subroutine convert_fv_to_dg_bathymetry
#endif                   

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

#if defined(_SWE_DG)                
		!multiplies a scalar with a dof state vector
                subroutine get_dofs_dg(f,q)
                  integer ::i
                  class (num_cell_data_pers), intent(in)		:: f
                  real (kind = GRID_SR),intent(out)               :: q (size(f%Q_DG,1),3)
                  q(:,1)= f%Q_DG(:)%h
                  q(:,2)= f%Q_DG(:)%p(1)
                  q(:,3)= f%Q_DG(:)%p(2)
                end subroutine get_dofs_dg

		subroutine set_dofs_dg(f,q)
                  class (num_cell_data_pers),intent(inout) 	:: f
                  real (kind = GRID_SR) 		        :: q(size(f%q_dg,1),3)
                  f%Q_DG(:)%H = q(:,1)
                  f%Q_DG(:)%p(1) =q(:,2)
                  f%Q_DG(:)%p(2) =q(:,3)
                end subroutine set_dofs_dg
#endif
	END MODULE SWE_data_types
#endif
