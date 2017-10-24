! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE)
	MODULE SWE_Displace
		use SFC_edge_traversal
		use SWE_euler_timestep
		use SWE_initialize_bathymetry
#       if defined(_SWE_PATCH)
            use SWE_PATCH
#       endif
#       if defined(_SWE_DG)
1           use SWE_DG_solver
#       endif            
		use Samoa_swe

		implicit none

        type num_traversal_data
            integer (kind = GRID_DI)			:: i_refinements_issued
        end type

		type(t_gv_Q)							:: gv_Q
		type(t_lfs_flux)						:: lfs_flux

#		define _GT_NAME							t_swe_displace_traversal

#		define _GT_EDGES

#		define _GT_ELEMENT_OP					element_op

#		define _GT_CELL_TO_EDGE_OP				cell_to_edge_op_dg

!#		define _GT_NODE_MPI_TYPE

#		include "SFC_generic_traversal_ringbuffer.f90"

!         subroutine create_node_mpi_type(mpi_node_type)
!             integer, intent(out)            :: mpi_node_type

! #           if defined(_MPI)
!                 type(t_node_data)                       :: node
!                 integer                                 :: blocklengths(2), types(2), disps(2), type_size, i_error
!                 integer (kind = MPI_ADDRESS_KIND)       :: lb, ub

!                 blocklengths(1) = 1
!                 blocklengths(2) = 1

!                 disps(1) = 0
!                 disps(2) = sizeof(node)

!                 types(1) = MPI_LB
!                 types(2) = MPI_UB

!                 call MPI_Type_struct(2, blocklengths, disps, types, mpi_node_type, i_error); assert_eq(i_error, 0)
!                 call MPI_Type_commit(mpi_node_type, i_error); assert_eq(i_error, 0)

!                 call MPI_Type_size(mpi_node_type, type_size, i_error); assert_eq(i_error, 0)
!                 call MPI_Type_get_extent(mpi_node_type, lb, ub, i_error); assert_eq(i_error, 0)

!                 assert_eq(0, lb)
!                 assert_eq(0, type_size)
!                 assert_eq(sizeof(node), ub)
! #           endif
!         end subroutine

		!******************
		!Geometry operators
		!******************

		subroutine element_op(traversal, section, element)
			type(t_swe_displace_traversal), intent(inout)		:: traversal
			type(t_grid_section), intent(inout)					:: section
			type(t_element_base), intent(inout)					:: element

			type(t_state)			                            :: Q(_SWE_CELL_SIZE)

			call gv_Q%read(element, Q)

			call alpha_volume_op(traversal, section, element, Q)

			call gv_Q%write(element, Q)
		end subroutine

		!*******************************
		!Volume and DoF operators
		!*******************************

		subroutine alpha_volume_op(traversal, section, element, Q)
			type(t_swe_displace_traversal), intent(inout)				:: traversal
			type(t_grid_section), intent(inout)							:: section
			type(t_element_base), intent(inout)						    :: element
			type(t_state), intent(inout)	                            :: Q(:)
            integer (kind = GRID_SI)                                    :: i
#           if defined(_SWE_PATCH)
            real (kind = GRID_SR), dimension(_SWE_PATCH_ORDER_SQUARE)        :: db
                real (kind = GRID_SR), dimension(_SWE_DG_DOFS)        :: db_dg            
                real (kind = GRID_SR), DIMENSION(2)                             :: r_coords     !< cell coords within patch
                integer (kind = GRID_SI)                                        :: j, row, col, cell_id
#           else
                real (kind = GRID_SR)                                       :: db(_SWE_CELL_SIZE)
#           endif

                !evaluate initial function values at dof positions and compute DoFs

                
#  if defined(_ASAGI)

#if defined(_SWE_DG)
                if( element%cell%data_pers%troubled .le. 0 .or. element%cell%data_pers%troubled .ge. 6) then
                   db_dg = -element%cell%data_pers%Q_DG%B + get_bathymetry_at_dg_patch(section, element, section%r_time)

                   if( .not.all(element%cell%data_pers%Q_DG%H +db_dg > cfg%coast_height))then
                      element%cell%data_pers%troubled = 1
                      ! print*,db_dg
                      ! print*,element%cell%data_pers%Q_DG%H
                      ! print*,element%cell%data_pers%Q_DG%p(1)
                      ! print*,element%cell%data_pers%Q_DG%p(2)
                      ! print*,element%cell%data_pers%Q_DG%B
                      
                      call apply_phi(element%cell%data_pers%Q_DG%H+element%cell%data_pers%Q_DG%B+db_dg,element%cell%data_pers%H)
                      call apply_phi(element%cell%data_pers%Q_DG%p(1),element%cell%data_pers%HU)
                      call apply_phi(element%cell%data_pers%Q_DG%p(2),element%cell%data_pers%HV)
                      element%cell%data_pers%B=get_bathymetry_at_patch(section, element, section%r_time)
                      ! print*
                      ! print*,element%cell%data_pers%H
                      ! print*,element%cell%data_pers%HU
                      ! print*,element%cell%data_pers%HV
                      ! print*,element%cell%data_pers%B


                   else
                      element%cell%data_pers%Q_DG%H = element%cell%data_pers%Q_DG%H + db_dg
                      element%cell%data_pers%Q_DG%B = element%cell%data_pers%Q_DG%B + db_dg
                      call bathymetry_derivatives(element%cell%data_pers,ref_plotter_data(abs(element%cell%geometry%i_plotter_type))%jacobian)
                   end if

                else
#endif                

#if defined (_SWE_PATCH)

                   db = -element%cell%data_pers%B + get_bathymetry_at_patch(section, element, section%r_time)
                   ! if(any(abs(db).ge.50))then
                   !    print*,element%cell%data_pers%troubled
                   !    print*,-element%cell%data_pers%B
                   !    print*,get_bathymetry_at_patch(section, element, section%r_time)
                   !    print*,element%cell%data_pers%H
                   !    stop
                   ! end if
                element%cell%data_pers%H = element%cell%data_pers%H + db
                element%cell%data_pers%B = element%cell%data_pers%B + db
#else
                db = -Q%b + get_bathymetry_at_element(section, element, section%r_time)
                Q%h = Q%h + db
                Q%b = Q%b + db
#endif
#           endif
                    
#if defined(_SWE_DG)
             end if
#endif


            !no coarsening while the earthquake takes place
			element%cell%geometry%refinement = max(0, element%cell%geometry%refinement)
   
		end subroutine
	END MODULE
#endif
