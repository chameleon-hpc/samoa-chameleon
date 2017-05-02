! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE)
	MODULE SWE_data_types
#if defined(_SWE_DG)
   use SWE_DG_matrices
#define _REF_TRIANGLE_SIZE_INV  (2.0q0 * real(_SWE_PATCH_ORDER_SQUARE,kind=kind(1.0q0)))
   
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
		end type

		!> cell state vector including bathymetry
		type, extends(t_dof_state) :: t_state
     real (kind = GRID_SR) :: b						!< bathymetry

#if defined(_SWE_DG_NODAL)
                        real(kind=GRID_SR) ::       b_x
                        real(kind=GRID_SR) ::       b_y
#endif
                        
            contains

            procedure, pass :: add_state => state_add
            generic :: operator(+) => add_state
         end type t_state

		!> update vector
		type, extends(t_dof_state) :: t_update
			real (kind = GRID_SR)													:: max_wave_speed			!< maximum wave speed required to compute the CFL condition
            contains

            procedure, pass :: add_update => update_add
            generic :: operator(+) => add_update
		end type

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

#		if defined(_SWE_PATCH)
#                       if defined(_SWE_DG)
                        type(t_state), DIMENSION(_SWE_DG_DOFS)                     :: Q_DG
                        type(t_dof_state), DIMENSION(_SWE_DG_DOFS*(_SWE_DG_ORDER+1))   :: Q_DG_P
                        real(kind=GRID_SR) :: dt

#if !defined(_SWE_DG_NODAL)
                        real(kind=GRID_SR) :: b_x(size(st_gl_node_vals,1)), b_y(size(st_gl_node_vals,1))
#endif                       

                        integer :: troubled = HUGE(1)
#                       endif
			type(t_state), DIMENSION(_SWE_CELL_SIZE)			:: Q !TODO: remove this and others Qs --> must handle conflicts with t_gv_Q methods afterwards...
			real (kind = GRID_SR), DIMENSION(_SWE_PATCH_ORDER_SQUARE)	:: H, HU, HV, B !< unknowns + bathymetry in triangular patch
#		else
			type(t_state), DIMENSION(_SWE_CELL_SIZE)				:: Q						!< cell status vector
#		endif

#               if defined(_SWE_DG)
                        contains

                          procedure :: convert_fv_to_dg => convert_fv_to_dg
                          procedure :: convert_dg_to_fv => convert_dg_to_fv
                          procedure :: convert_dg_to_fv_bathymetry => convert_dg_to_fv_bathymetry
                          procedure :: convert_fv_to_dg_bathymetry => convert_fv_to_dg_bathymetry

                          procedure :: get_dofs_dg => get_dofs_dg
                          procedure :: get_dofs_pred => get_dofs_pred
                          procedure :: set_dofs_dg => set_dofs_dg
                          procedure :: set_dofs_pred => set_dofs_pred

                          ! procedure :: vec_to_dofs_dg => vec_to_dofs_dg
                          ! procedure :: vec_to_dofs_dg_p => vec_to_dofs_dg_p
                          ! procedure :: dofs_to_vec_dg  => dofs_to_vec_dg
                          ! procedure :: dofs_to_vec_dg_p => dofs_to_vec_dg_p

#               endif
		END type num_cell_data_pers

		!> Cell representation on an edge, this would typically be everything required from a cell to compute the flux function on an edge
		type num_cell_rep

#if defined(_SWE_PATCH)
     type(t_state), DIMENSION(_SWE_EDGE_SIZE)				:: Q !TODO: remove this and others Qs --> must handle conflicts with t_gv_Q methods afterwards...
			real (kind = GRID_SR), dimension (_SWE_PATCH_ORDER)		:: H, HU, HV, B !< edge stores ghost cells for communication of ghost cells

#if defined(_SWE_DG)
   type(t_state), DIMENSION((_SWE_DG_ORDER+1)*(_SWE_DG_ORDER+1))   :: Q_DG_P
   type(t_state), DIMENSION(_SWE_DG_DOFS)   :: Q_DG
   
!                        real(kind=GRID_SR),Dimension(3) :: min_val,max_val
   integer :: troubled

#endif

#else
!   type(t_state), DIMENSION(1)				:: Q !TODO: remove this and others Qs --> must handle conflicts with t_gv_Q methods afterwards...
   
#endif

                     end type num_cell_rep

		!> Cell update, this would typically be a flux function
		type num_cell_update
#if defined (_SWE_DG)
                        real (kind = GRID_SR), DIMENSION(_SWE_PATCH_ORDER)							:: H, HU, HV, B !< values of ghost cells
!                      type(t_state), DIMENSION((_SWE_DG_ORDER+1)*(_SWE_DG_ORDER+1))   :: Q_DG_P
                        type(t_update), DIMENSION((_SWE_DG_ORDER+1)*(_SWE_DG_ORDER+1))   :: flux

                        !dofs at t_n
                        type(t_state), DIMENSION(_SWE_DG_DOFS)   :: Q_DG

                        integer :: troubled
#else

#if defined (_SWE_PATCH)
                    real (kind = GRID_SR), DIMENSION(_SWE_PATCH_ORDER)							:: H, HU, HV, B !< values of ghost cells
#else                        
                    type(t_update), DIMENSION(_SWE_EDGE_SIZE)									:: flux						!< cell update
#endif
!_SWE_DG
#endif
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

                   call lusolve(mue_lu,_SWE_DG_DOFS+1,mue_lu_pivot,q_temp(:,1))
                   call lusolve(mue_lu,_SWE_DG_DOFS+1,mue_lu_pivot,q_temp(:,2))
                   call lusolve(mue_lu,_SWE_DG_DOFS+1,mue_lu_pivot,q_temp(:,3))
                   
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


                 end subroutine convert_dg_to_fv

                 subroutine apply_phi(dg,fv)
                   real(kind=GRID_SR),intent(out) :: fv(_SWE_PATCH_ORDER_SQUARE)
                   real(kind=GRID_SR),intent(in)  :: dg(_SWE_DG_DOFS)

                   fv=matmul(phi,dg)*_REF_TRIANGLE_SIZE_INV
                   
                 end subroutine apply_phi


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
#if defined(_SWE_DG_NODAL)                   
                   real(kind=GRID_SR) :: b_x_temp(_SWE_DG_DOFS)
                   real(kind=GRID_SR) :: b_y_temp(_SWE_DG_DOFS)
#else
                   real(kind=GRID_SR) :: b_x_temp(size(st_der_x_gl_node_vals,1))
                   real(kind=GRID_SR) :: b_y_temp(size(st_der_y_gl_node_vals,1))

#endif
                   real(kind=GRID_SR),intent(in) :: normals(2,2)
                   real(kind=GRID_SR) :: normals_normed(2,2)
                   integer ::i,j
                   real(kind=GRID_SR) :: epsilon=0.0_GRID_SR


                   normals_normed=normals/NORM2(normals(1,1:2))

                   b_temp(1:_SWE_DG_DOFS)= 2.0q0*matmul(transpose(phi),dofs%b)
                   
                   b_temp(_SWE_DG_DOFS+1) = sum(dofs%b)
                   
                   b_temp = b_temp /_REF_TRIANGLE_SIZE_INV

                  
                   call lusolve(mue_lu,_SWE_DG_DOFS+1,mue_lu_pivot,b_temp)

                   do i=1,_SWE_DG_DOFS
                      dofs%Q_DG(i)%b=b_temp(i)
                   end do

!                   do i=0,_SWE_DG_ORDER
!                      dofs%Q_DG_P(1+_SWE_DG_DOFS*i:_SWE_DG_DOFS*(i+1))%b=dofs%Q_DG(:)%b
!                   end do
#if !defined(_SWE_DG_NODAL)

                   ! b_x_temp=matmul(st_der_x_gl_node_vals,dofs%Q_DG_P(:)%b)
                   ! b_y_temp=matmul(st_der_y_gl_node_vals,dofs%Q_DG_P(:)%b)
                   
                   ! dofs%b_x = normals_normed(1,1) * b_x_temp + normals_normed(1,2) * b_y_temp
                   ! dofs%b_y = normals_normed(2,1) * b_x_temp + normals_normed(2,2) * b_y_temp

#else
                   b_x_temp=matmul(basis_der_x,b_temp(1:_SWE_DG_DOFS))
                   b_y_temp=matmul(basis_der_y,b_temp(1:_SWE_DG_DOFS))

                   dofs%Q_DG(:)%b_x=normals_normed(1,1) * b_x_temp + normals_normed(1,2) * b_y_temp
                   dofs%Q_DG(:)%b_y=normals_normed(2,1) * b_x_temp + normals_normed(2,2) * b_y_temp

                   ! do i=0,_SWE_DG_ORDER
                   !    dofs%Q_DG_P(1+_SWE_DG_DOFS*i:_SWE_DG_DOFS*(i+1))%b_x=dofs%Q_DG(:)%b_x
                   !    dofs%Q_DG_P(1+_SWE_DG_DOFS*i:_SWE_DG_DOFS*(i+1))%b_y=dofs%Q_DG(:)%b_y
                   ! end do
#endif
                 end subroutine convert_fv_to_dg_bathymetry

#endif                   

		!adds two state vectors
		elemental function state_add(Q1, Q2)	result(Q_out)
			class (t_state), intent(in)		:: Q1
			type (t_state), intent(in)		:: Q2
			type (t_state)					:: Q_out
#if !defined(_SWE_DG_NODAL)
			Q_out = t_state(Q1%h + Q2%h, Q1%p + Q2%p, Q1%b + Q2%b)
#else
			Q_out = t_state(Q1%h + Q2%h, Q1%p + Q2%p, Q1%b + Q2%b,Q1%b_x + Q2%b_x,Q1%b_y + Q2%b_y)
#endif
		end function

		!adds two update vectors
		elemental function update_add(f1, f2)	result(f_out)
			class (t_update), intent(in)		:: f1
			type (t_update), intent(in)		    :: f2
			type (t_update)					    :: f_out

			f_out = t_update(f1%h + f2%h, f1%p + f2%p, max_wave_speed = max(f1%max_wave_speed, f2%max_wave_speed))
		end function

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

                subroutine get_dofs_pred(f,q)
                  class (num_cell_data_pers), intent(in)		:: f
                  real (kind = GRID_SR),intent(out)               :: q (size(f%Q_DG_P,1),3)
                  q(:,1)= f%Q_DG_P(:)%h
                  q(:,2)= f%Q_DG_P(:)%p(1)
                  q(:,3)= f%Q_DG_P(:)%p(2)
                end subroutine get_dofs_pred


		subroutine set_dofs_dg(f,q)
                  class (num_cell_data_pers),intent(inout) 	:: f
                  real (kind = GRID_SR) 		        :: q(size(f%q_dg,1),3)
                  f%Q_DG(:)%H = q(:,1)
                  f%Q_DG(:)%p(1) =q(:,2)
                  f%Q_DG(:)%p(2) =q(:,3)
                end subroutine set_dofs_dg

		subroutine set_dofs_pred(f,q)
                  class (num_cell_data_pers),intent(inout)      :: f
                  real (kind = GRID_SR) 		        :: q(size(f%q_dg_p,1),3)
                  f%Q_DG_P(:)%H    = q(:,1)
                  f%Q_DG_P(:)%p(1) = q(:,2)
                  f%Q_DG_P(:)%p(2) = q(:,3)
                end subroutine set_dofs_pred
#endif
	END MODULE SWE_data_types
#endif
