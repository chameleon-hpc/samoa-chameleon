! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE)
  MODULE SWE_Adapt
    use SFC_edge_traversal
    use Conformity
    
    use Samoa_swe
    use Tools_noise
    use SWE_initialize_bathymetry
#if defined (_SWE_PATCH)
    use SWE_PATCH
#endif
#if defined(_SWE_DG)
    use SWE_dg_predictor
    use SWE_dg_solver
    use SWE_DG_Limiter    
#endif

    implicit none

    type num_traversal_data
#if defined (_SWE_PATCH)
       logical, DIMENSION(_SWE_PATCH_ORDER_SQUARE)                         :: dry_cell ! used only in coarsen operator
#else
       type(t_state), dimension(_SWE_CELL_SIZE, 2)                         :: Q_in
#endif
    end type num_traversal_data


		type(t_gv_Q)							:: gv_Q
		type(t_lfs_flux)						:: lfs_flux

#		define _GT_NAME					t_swe_adaption_traversal

#		define _GT_EDGES

#               undef _GT_SKELETON_OP  
#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_PRE_TRAVERSAL_OP			pre_traversal_op
#		define _GT_POST_TRAVERSAL_OP			post_traversal_op

#		define _GT_TRANSFER_OP				transfer_op
#		define _GT_REFINE_OP				refine_op
#		define _GT_COARSEN_OP				coarsen_op

#if defined(_SWE_DG)
#		define _GT_CELL_TO_EDGE_OP			cell_to_edge_op_dg
#		define _GT_CELL_UPDATE_OP		        cell_update_op_adapt
#else
#		define _GT_CELL_TO_EDGE_OP			cell_to_edge_op
#endif

!#		define _GT_NODE_MPI_TYPE

!#		define _GT_NODE_WRITE_OP			node_write_op
  
!#		define _GT_EDGE_WRITE_OP			edge_write_op

!#               define _GT_EDGE_MERGE_OP                edge_merge_op_dg

!#               define _GT_CELL_FIRST_TOUCH_OP                  cell_first_touch_op_dg
#		include "SFC_generic_adaptive_traversal.f90"


  subroutine create_node_mpi_type(mpi_node_type)
          integer, intent(out)            :: mpi_node_type

#if defined(_MPI)
          type(t_node_data)                       :: node
          integer                                 :: blocklengths(2), types(2), disps(2), type_size, i_error
          integer (kind = MPI_ADDRESS_KIND)       :: lb, ub
          
          blocklengths(1) = 1
          blocklengths(2) = 1
          
          disps(1) = 0
          disps(2) = sizeof(node)
          
          types(1) = MPI_LB
          types(2) = MPI_UB
          
          call MPI_Type_struct(2, blocklengths, disps, types, mpi_node_type, i_error); assert_eq(i_error, 0)
          call MPI_Type_commit(mpi_node_type, i_error); assert_eq(i_error, 0)
          
          call MPI_Type_size(mpi_node_type, type_size, i_error); assert_eq(i_error, 0)
          call MPI_Type_get_extent(mpi_node_type, lb, ub, i_error); assert_eq(i_error, 0)
          
          assert_eq(0, lb)
          assert_eq(0, type_size)
          assert_eq(sizeof(node), ub)
#endif
        end subroutine


#if defined(_SWE_DG)

        subroutine cell_update_op_adapt(traversal, section, element, update1, update2, update3)
          type(t_swe_adaption_traversal), intent(inout)				:: traversal
          type(t_grid_section), intent(inout)							:: section
          type(t_element_base), intent(inout)						:: element
          type(num_cell_update), intent(inout)						:: update1, update2, update3
          integer :: i
          call updateCellStatus(element%cell%data_pers)
          if(isDG(element%cell%data_pers%troubled)) then
             call dg_predictor(element,section%r_dt)
#if defined (_DEBUG)             
             element%cell%data_pers%debug_flag = -1
          else
             element%cell%data_pers%debug_flag = element%cell%data_pers%troubled
#endif                          
          end if
          
        end subroutine cell_update_op_adapt


        subroutine edge_merge_op_dg(local_edge,neighbor_edge)
          type(t_edge_data), intent(inout)   ::local_edge
          type(t_edge_data), intent(in)      ::neighbor_edge
          type(t_state),Allocatable          ::rep_temp(:)
          integer                            ::i,j
          local_edge%rep%troubled=neighbor_edge%rep%troubled
          local_edge%rep%QP = neighbor_edge%rep%QP
          local_edge%rep%FP = neighbor_edge%rep%FP
        end subroutine edge_merge_op_dg

#endif

        
        subroutine pre_traversal_grid_op(traversal, grid)
          type(t_swe_adaption_traversal), intent(inout)				:: traversal
          type(t_grid), intent(inout)							        :: grid
#if defined(_SWE_DG)

          if (cfg%r_max_time > 0.0_SR) then
             grid%r_dt = min(cfg%r_max_time-grid%r_time, grid%r_dt)
          end if

          !--- Needed for output at distinc timestep ---!
          if(cfg%r_output_time_step > 0) then
             if(mod(grid%r_time,cfg%r_output_time_step) > 0.0_GRID_SR) then
                grid%r_dt = min(cfg%r_output_time_step-mod(grid%r_time,cfg%r_output_time_step), grid%r_dt)
             end if
          end if
          
          ! if (cfg%r_output_time_step > 0.0_SR) then
          !    grid%r_dt = min(cfg%r_output_time_step, grid%r_dt)
          ! end if

#           if defined(_ASAGI)
          !if we are in the earthquake phase, limit the simulation time step by the earthquake time step
          if (grid%r_time < cfg%t_max_eq) then
             grid%r_dt = min(grid%r_dt, cfg%dt_eq)
          end if
#           endif

!_SWE_DG          
#endif
    call scatter(grid%r_dt, grid%sections%elements_alloc%r_dt)


  end subroutine pre_traversal_grid_op

        subroutine pre_traversal_op(traversal, section)
          type(t_swe_adaption_traversal), intent(inout)				:: traversal
          type(t_grid_section), intent(inout)							:: section
        end subroutine pre_traversal_op

        subroutine post_traversal_op(traversal, section)
          type(t_swe_adaption_traversal), intent(inout)				:: traversal
          type(t_grid_section), intent(inout)							:: section
        end subroutine post_traversal_op
        !******************
        !Adaption operators
        !******************

        subroutine transfer_op(traversal, section, src_element, dest_element)
          type(t_swe_adaption_traversal), intent(inout)	                            :: traversal
          type(t_grid_section), intent(inout)											:: section
          type(t_traversal_element), intent(inout)									:: src_element
          type(t_traversal_element), intent(inout)									:: dest_element

#if defined(_SWE_PATCH)
          type(t_state), dimension(_SWE_PATCH_ORDER_SQUARE)						:: Q
          real (kind = GRID_SR), dimension(_SWE_PATCH_ORDER_SQUARE)				:: H, HU, HV, B
          integer :: i
#else
          type(t_state), dimension(_SWE_CELL_SIZE)								:: Q

#endif

#if defined (_SWE_PATCH)
          dest_element%cell%data_pers%H = src_element%cell%data_pers%H
          dest_element%cell%data_pers%HU = src_element%cell%data_pers%HU
          dest_element%cell%data_pers%HV = src_element%cell%data_pers%HV
          dest_element%cell%data_pers%B = src_element%cell%data_pers%B

#else
          call gv_Q%read( src_element%t_element_base, Q)
          call gv_Q%write( dest_element%t_element_base, Q)
#endif
#if defined (_SWE_DG)
          do i=1,size(dest_element%edges)
             dest_element%edges(i)%ptr%rep = src_element%edges(i)%ptr%rep
             dest_element%edges(i)%ptr%data_pers = src_element%edges(i)%ptr%data_pers
          end do

          dest_element%cell%data_pers%QP = src_element%cell%data_pers%QP
          dest_element%cell%data_pers%FP = src_element%cell%data_pers%FP
          dest_element%cell%data_pers%Q = src_element%cell%data_pers%Q
          dest_element%cell%data_pers%troubled=src_element%cell%data_pers%troubled
#endif
        end subroutine transfer_op
        

        subroutine refine_op(traversal, section, src_element, dest_element, refinement_path)
          type(t_swe_adaption_traversal), intent(inout)	                        :: traversal
          type(t_grid_section), intent(inout)					:: section
          type(t_traversal_element), intent(inout)				:: src_element
          type(t_traversal_element), intent(inout)				:: dest_element
          integer, dimension(:), intent(in)					:: refinement_path
          
#if defined(_SWE_PATCH)
                real (kind=GRID_SR), dimension(_SWE_PATCH_ORDER_SQUARE)                     :: H_in, HU_in, HV_in, B_in
                real (kind=GRID_SR), dimension(_SWE_PATCH_ORDER_SQUARE)                     :: H_out, HU_out, HV_out, B_out
                integer                                                                     :: i_plotter_type, j, row, col, cell_id
                integer, DIMENSION(_SWE_PATCH_ORDER_SQUARE,2)                               :: child
                real (kind = GRID_SR), DIMENSION(2)                                         :: r_coords     !< cell coords within patch
                logical, DIMENSION(_SWE_PATCH_ORDER_SQUARE)                                 :: dry_cell_in, dry_cell_out
#else
          type(t_state), dimension(_SWE_CELL_SIZE)                                    :: Q_in
          type(t_state), dimension(_SWE_CELL_SIZE, 2)                                 :: Q_out
          logical                                                                     :: dry_cell

#endif
          integer					                              :: i


#if defined (_SWE_PATCH)
#if defined (_SWE_DG)
          ! print*,"src"
          ! print*,src_element%cell%data_pers%troubled
          ! print*,src_element%cell%data_pers%Q%h
          ! print*,src_element%cell%data_pers%Q%p(1)
          ! print*,src_element%cell%data_pers%Q%p(2)
          ! print*,src_element%cell%data_pers%Q%b          
          
          i_plotter_type = src_element%cell%geometry%i_plotter_type
!          dest_element%cell%data_pers%troubled=src_element%cell%data_pers%troubled
          
          if(.not.isCoast(src_element%cell%data_pers%troubled)) then
             do i=1, size(refinement_path)
                ! decide which child is being computed (first = left to hypot, second = right to hypot.)
                if ( (refinement_path(i) == 1 .and. i_plotter_type>0) .or. (refinement_path(i) == 2 .and. i_plotter_type<0)) then
                   dest_element%cell%data_pers%Q%H=matmul(ref2,src_element%cell%data_pers%Q%H+src_element%cell%data_pers%Q%B)
                   dest_element%cell%data_pers%Q%p(1)=matmul(ref2,src_element%cell%data_pers%Q%p(1))
                   dest_element%cell%data_pers%Q%p(2)=matmul(ref2,src_element%cell%data_pers%Q%p(2))
                else
                   dest_element%cell%data_pers%Q%H=matmul(ref1,src_element%cell%data_pers%Q%H+src_element%cell%data_pers%Q%B)
                   dest_element%cell%data_pers%Q%p(1)=matmul(ref1,src_element%cell%data_pers%Q%p(1))
                   dest_element%cell%data_pers%Q%p(2)=matmul(ref1,src_element%cell%data_pers%Q%p(2))
                end if

                dest_element%cell%data_pers%Q%H=matmul(s_m_inv,dest_element%cell%data_pers%Q%H)
                dest_element%cell%data_pers%Q%p(1)=matmul(s_m_inv,dest_element%cell%data_pers%Q%p(1))
                dest_element%cell%data_pers%Q%p(2)=matmul(s_m_inv,dest_element%cell%data_pers%Q%p(2))
                dest_element%cell%data_pers%Q%B = get_bathymetry_at_dg_patch(section, dest_element%t_element_base, section%r_time)
                !call dest_element%cell%data_pers%convert_fv_to_dg_bathymetry(ref_plotter_data(abs(i_plotter_type))%jacobian)
                dest_element%cell%data_pers%Q%H=dest_element%cell%data_pers%Q%H-dest_element%cell%data_pers%Q%B
             end do
             
          end if
          
          call updateCellStatus(dest_element%cell%data_pers)
          if(isDG(dest_element%cell%data_pers%troubled)) then
             call dg_predictor(dest_element,section%r_dt)
          end if
          

#endif
             
       if(dest_element%cell%data_pers%troubled.ge.1) then

#endif             
          i_plotter_type = src_element%cell%geometry%i_plotter_type

          H_in = src_element%cell%data_pers%H
          HU_in = src_element%cell%data_pers%HU
          HV_in = src_element%cell%data_pers%HV
          B_in = src_element%cell%data_pers%B

          dry_cell_in = .false.
          dry_cell_out = .false.
          
                do i=1, size(refinement_path)
                    ! compute 1st or 2nd child depending on orientation and refinement_path(i)
                    ! and store it in output arrays.
                    ! Store it in input arrays for next refinement.
                    
                    ! decide which child is being computed (first = left to hypot, second = right to hypot.)
                    if ( (refinement_path(i) == 1 .and. i_plotter_type>0) .or. (refinement_path(i) == 2 .and. i_plotter_type<0)) then
                        child = SWE_PATCH_GEOMETRY%first_child
                    else
                        child = SWE_PATCH_GEOMETRY%second_child
                    end if
                    
                    ! the child patch values are always averages of two values from the parent patch
                    do j=1, _SWE_PATCH_ORDER_SQUARE
                        H_out(j) = 0.5*( H_in(child(j,1)) + H_in(child(j,2)) )
                        HU_out(j) = 0.5*( HU_in(child(j,1)) + HU_in(child(j,2)) )
                        HV_out(j) = 0.5*( HV_in(child(j,1)) + HV_in(child(j,2)) )
                        
                        ! find out if resulting fine cell was derived from a dry cell
                        if (i == 1 .and. (H_in(child(j,1)) < B_in(child(j,1)) + cfg%dry_tolerance .or. H_in(child(j,2)) < B_in(child(j,2)) + cfg%dry_tolerance)) then
                            dry_cell_out(j) = .true.
                        end if
                        if (i == 2 .and. (dry_cell_in(child(j,1)) .or. dry_cell_in(child(j,2)))) then
                            dry_cell_out(j) = .true.
                        end if
                    end do
                    
                    ! refined patches will have inverse orientation
                    i_plotter_type = -1 * i_plotter_type
                    
                    ! copy to input arrays for next iteration
                    H_in = H_out
                    HU_in = HU_out
                    HV_in = HV_out
                    dry_cell_in = dry_cell_out
                 end do
                
                ! store results in dest_element
                dest_element%cell%data_pers%H = H_out
                dest_element%cell%data_pers%HU = HU_out
                dest_element%cell%data_pers%HV = HV_out
                
                ! get bathymetry
                dest_element%cell%data_pers%B = get_bathymetry_at_patch(section, dest_element%t_element_base, section%r_time)
                
!!#if defined(_ASAGI)q
          ! if cell was initially dry, we need to check if the fine cells should be initialized with h=0
                where (dry_cell_out(:))
#if defined(_ASAGI)                           
             dest_element%cell%data_pers%H = 0.0_GRID_SR
#else
             dest_element%cell%data_pers%H = dest_element%cell%data_pers%B
#endif                           
                   
             dest_element%cell%data_pers%HU = 0.0_GRID_SR
             dest_element%cell%data_pers%HV = 0.0_GRID_SR
          end where
!!#endif

#else
          call gv_Q%read( src_element%t_element_base, Q_in)

          !convert momentum to velocity
          !Q_in(1)%p = 1.0_GRID_SR / (Q_in(1)%h - Q_in(1)%b) * Q_in(1)%p

          ! find out if resulting fine cell was derived from a dry cell
          if (Q_in(1)%h < Q_in(1)%b + cfg%dry_tolerance) then
             dry_cell = .true.
          else
             dry_cell = .false.
          end if

          do i = 1, size(refinement_path)
             call t_basis_Q_split(Q_in%h, 	Q_out(:, 1)%h, 		Q_out(:, 2)%h)
             call t_basis_Q_split(Q_in%p(1),	Q_out(:, 1)%p(1),	Q_out(:, 2)%p(1))
             call t_basis_Q_split(Q_in%p(2),	Q_out(:, 1)%p(2),	Q_out(:, 2)%p(2))

             Q_in = Q_out(:, refinement_path(i))
          end do

          !convert velocity back to momentum
          !Q_in(1)%p = (Q_in(1)%h - Q_in(1)%b) * Q_in(1)%p

          Q_in%b = get_bathymetry_at_element(section, dest_element%t_element_base, section%r_time)

#if defined(_ASAGI)
          ! if cell was initially dry, we need to check if the fine cells should be initialized with h=0
          if (dry_cell) then
             Q_in%h = max(0.0_GRID_SR, Q_in%b)
             Q_in(1)%p = [0.0_GRID_SR, 0.0_GRID_SR]
          end if
#endif

          call gv_Q%write( dets_element%t_element_base, Q_in)
#endif

#if defined(_SWE_DG)

          end if
!          call dest_element%cell%data_pers%convert_fv_to_dg

!          call dest_element%cell%data_pers%convert_fv_to_dg_bathymetry(ref_plotter_data(abs(i_plotter_type))%jacobian)
          
          !          dest_element%cell%data_pers%troubled=src_element%cell%data_pers%troubled
          ! print*,"dest"
          ! print*,dest_element%cell%data_pers%troubled
          ! print*,dest_element%cell%data_pers%Q%h
          ! print*,dest_element%cell%data_pers%Q%p(1)
          ! print*,dest_element%cell%data_pers%Q%p(2)
          ! print*,dest_element%cell%data_pers%Q%b          



#if defined (_DEBUG)                        
           dest_element%cell%data_pers%debug_flag=-4
#endif           
        end subroutine refine_op

        subroutine coarsen_op(traversal, section, src_element, dest_element, refinement_path)
          type(t_swe_adaption_traversal), intent(inout)	                :: traversal
          type(t_grid_section), intent(inout)													:: section
          type(t_traversal_element), intent(inout)									:: src_element
          type(t_traversal_element), intent(inout)									:: dest_element
          integer, dimension(:), intent(in)											:: refinement_path
          integer                                                                     :: i
#if defined(_SWE_PATCH)
          integer                                                                     :: j,i_plotter_type
          integer, DIMENSION(_SWE_PATCH_ORDER_SQUARE,2)                               :: child
          logical                                                                     :: fv_coarse , dg_coarse , first
#else
          type(t_state), dimension(_SWE_CELL_SIZE)                                :: Q_out
#endif

#if defined(_SWE_PATCH)
#if defined(_SWE_DG)
          

          i_plotter_type=dest_element%cell%geometry%i_plotter_type            
          first = refinement_path(1).eq.1
          
          if(first) then
             if(src_element%cell%data_pers%troubled .le. 0) then
                ! dg = 1 fv = 1
                dg_coarse= .true.
                fv_coarse= .true.
                dest_element%cell%data_pers%troubled = DG
             else
                ! dg = 0 fv = 1
                dg_coarse= .false.
                fv_coarse= .true.
                dest_element%cell%data_pers%troubled = TROUBLED
             end if
          else 
             if(src_element%cell%data_pers%troubled .le. 0) then
                dg_coarse = merge(.true.,.false.,dest_element%cell%data_pers%troubled.le.0)                
                fv_coarse = .not.dg_coarse
                dest_element%cell%data_pers%troubled = DG
             else
                dg_coarse= .false.
                fv_coarse= .true.
                dest_element%cell%data_pers%troubled = TROUBLED
             end if
          end if

          if(dg_coarse) then
                if (first) then
                   dest_element%cell%data_pers%Q%H   =0.0_GRID_SR
                   dest_element%cell%data_pers%Q%p(1)=0.0_GRID_SR
                   dest_element%cell%data_pers%Q%p(2)=0.0_GRID_SR
                   dest_element%cell%data_pers%Q%B = 0.0_GRID_SR

                end if
                
                if (first .and. i_plotter_type > 0 .or. .not.first .and. i_plotter_type < 0)then

                   !preserving constant water height
                   dest_element%cell%data_pers%Q%H   =dest_element%cell%data_pers%Q%H   +matmul(coarsen(:,_SWE_DG_DOFS+1:_SWE_DG_DOFS*2),src_element%cell%data_pers%Q%H &
                                                                                                                                            +  src_element%cell%data_pers%Q%B)
                   dest_element%cell%data_pers%Q%p(1)=dest_element%cell%data_pers%Q%p(1)+matmul(coarsen(:,_SWE_DG_DOFS+1:_SWE_DG_DOFS*2),src_element%cell%data_pers%Q%p(1))
                   dest_element%cell%data_pers%Q%p(2)=dest_element%cell%data_pers%Q%p(2)+matmul(coarsen(:,_SWE_DG_DOFS+1:_SWE_DG_DOFS*2),src_element%cell%data_pers%Q%p(2))
                else
                   !preserving constant water height
                   dest_element%cell%data_pers%Q%H   =dest_element%cell%data_pers%Q%H   +matmul(coarsen(:,1:_SWE_DG_DOFS),src_element%cell%data_pers%Q%H &
                                                                                                                              + src_element%cell%data_pers%Q%B)
                   dest_element%cell%data_pers%Q%p(1)=dest_element%cell%data_pers%Q%p(1)+matmul(coarsen(:,1:_SWE_DG_DOFS),src_element%cell%data_pers%Q%p(1))
                   dest_element%cell%data_pers%Q%p(2)=dest_element%cell%data_pers%Q%p(2)+matmul(coarsen(:,1:_SWE_DG_DOFS),src_element%cell%data_pers%Q%p(2))
                end if

                if(.not.first)then

                   ! call lusolve(s_m_lu, _SWE_DG_DOFS, s_m_lu_pivot,dest_element%cell%data_pers%Q%H)
                   ! call lusolve(s_m_lu, _SWE_DG_DOFS, s_m_lu_pivot,dest_element%cell%data_pers%Q%p(1))
                   ! call lusolve(s_m_lu, _SWE_DG_DOFS, s_m_lu_pivot,dest_element%cell%data_pers%Q%p(2))

                   dest_element%cell%data_pers%Q%H=matmul(s_m_inv,dest_element%cell%data_pers%Q%H)
                   dest_element%cell%data_pers%Q%p(1)=matmul(s_m_inv,dest_element%cell%data_pers%Q%p(1))
                   dest_element%cell%data_pers%Q%p(2)=matmul(s_m_inv,dest_element%cell%data_pers%Q%p(2))
                   !generate bythymetry from initial data
                   
!                   dest_element%cell%data_pers%B = get_bathymetry_at_patch(section, dest_element%t_element_base, section%r_time)
!                   call dest_element%cell%data_pers%convert_fv_to_dg_bathymetry(ref_plotter_data(abs(i_plotter_type))%jacobian)
                   dest_element%cell%data_pers%Q%B = get_bathymetry_at_dg_patch(section, dest_element%t_element_base, section%r_time)



                   dest_element%cell%data_pers%Q%H = dest_element%cell%data_pers%Q%H-dest_element%cell%data_pers%Q%B

                   ! if cell is dry after refinement do fv coarsening
                   if(isWetDryInterface(dest_element%cell%data_pers%Q%H)) then
                      fv_coarse = .true.
                      dest_element%cell%data_pers%troubled=WET_DRY_INTERFACE
                   end if
                   
                end if
             end if
          
             if(fv_coarse) then
                if(src_element%cell%data_pers%troubled.le.0) then
                   call apply_phi(src_element%cell%data_pers%Q%H,src_element%cell%data_pers%H)
                   call apply_phi(src_element%cell%data_pers%Q%p(1),src_element%cell%data_pers%HU)
                   call apply_phi(src_element%cell%data_pers%Q%p(2),src_element%cell%data_pers%HV)
                   call apply_phi(src_element%cell%data_pers%Q%b,src_element%cell%data_pers%B)
                   src_element%cell%data_pers%H=src_element%cell%data_pers%H+src_element%cell%data_pers%B
                end if

#endif
                !IMPORTANT: in the current samoa implementation, this subroutine is always called first with refinement_path=1, and then =2.
                ! The below implementation supposes this. If the samoa core implementation changes, this may become invalid!

                ! find out children order based on geometry
                if ((refinement_path(1) == 1 .and. dest_element%cell%geometry%i_plotter_type>0) .or. (refinement_path(1) == 2 .and. dest_element%cell%geometry%i_plotter_type<0)) then
                    child = SWE_PATCH_GEOMETRY%first_child
                else
                    child = SWE_PATCH_GEOMETRY%second_child
                end if
                
                if (refinement_path(1) == 1) then 
                    traversal%dry_cell = .false.
                end if
                
                associate (data => dest_element%cell%data_pers)
                    ! initialize all values with zero - the final value is an average of 4 children's cells
                    if (refinement_path(1) == 1) then
                        data%H = 0.0_GRID_SR
                        data%HU = 0.0_GRID_SR
                        data%HV = 0.0_GRID_SR
                        data%B = 0.0_GRID_SR
                    end if

                    ! sum all values to their respective cells
                    do i=1, _SWE_PATCH_ORDER_SQUARE
                        do j=1, 2
                            data%H(child(i,j)) = data%H(child(i,j)) + src_element%cell%data_pers%H(i)
                            data%HU(child(i,j)) = data%HU(child(i,j)) + src_element%cell%data_pers%HU(i)
                            data%HV(child(i,j)) = data%HV(child(i,j)) + src_element%cell%data_pers%HV(i)
                            data%B(child(i,j)) = data%B(child(i,j)) + src_element%cell%data_pers%B(i)
                            
                            if (src_element%cell%data_pers%H(i) < src_element%cell%data_pers%B(i) + cfg%dry_tolerance) then
                                traversal%dry_cell(child(i,j)) = .true.
                            end if
                        end do
                    end do

                    ! divide by 4 to compute the average of the 4 values
                    if (refinement_path(1) == 2) then
                        data%H = data%H * 0.25_GRID_SR
                        data%HU = data%HU * 0.25_GRID_SR
                        data%HV = data%HV * 0.25_GRID_SR
                        data%B = data%B * 0.25_GRID_SR
                        
!#if defined(_ASAGI)
                            ! if one of the cells was dry, we need to check if the coarsen cell should be initialized with h=0
                        where (traversal%dry_cell(:))
#if defined(_ASAGI)                           
                           data%H = 0.0_GRID_SR
#else
                           data%H = data%B
#endif                           
                           data%HU = 0.0_GRID_SR
                           data%HV = 0.0_GRID_SR
                        end where
                        !#endif
                        
                     end if

                   end associate
#else

        !state vector
        
        i = refinement_path(1)
        call gv_Q%read( src_element%t_element_base, traversal%Q_in(:, i))
        
        !convert momentum to velocity
        !traversal%Q_in(1, i)%p = 1.0_GRID_SR / (traversal%Q_in(1, i)%h - traversal%Q_in(1, i)%b) * traversal%Q_in(1, i)%p
        
        if (i > 1) then
           call t_basis_Q_merge(traversal%Q_in(:, 1)%h,		traversal%Q_in(:, 2)%h,		Q_out%h)
           call t_basis_Q_merge(traversal%Q_in(:, 1)%p(1),	    traversal%Q_in(:, 2)%p(1),	Q_out%p(1))
           call t_basis_Q_merge(traversal%Q_in(:, 1)%p(2),	    traversal%Q_in(:, 2)%p(2),	Q_out%p(2))
           call t_basis_Q_merge(traversal%Q_in(:, 1)%b,		traversal%Q_in(:, 2)%b,		Q_out%b)
           
           !convert velocity back to momentum
           !Q_out(1)%p = (Q_out(1)%h - Q_out(1)%b) * Q_out(1)%p
           
#if defined(_ASAGI)
           ! if one of the input cells was dry, we need to check if the coarse cell should be initialized with h=0
           if (traversal%Q_in(1, 1)%h < traversal%Q_in(1, 1)%b + cfg%dry_tolerance .or. traversal%Q_in(1, 2)%h < traversal%Q_in(1, 2)%b + cfg%dry_tolerance) then
              Q_out(1)%h = max(0.0_GRID_SR, Q_out(1)%b)
              Q_out(1)%p = [0.0_GRID_SR, 0.0_GRID_SR]
           end if
#endif
           call gv_Q%write( dest_element%t_element_base, Q_out)
        end if
#endif
           
#if defined(_SWE_DG)
     end if

     if(dest_element%cell%data_pers%troubled .le. 0 ) then
        call dg_predictor(dest_element,section%r_dt)
#if defined (_DEBUG)                     
        dest_element%cell%data_pers%debug_flag=-2
     else
        dest_element%cell%data_pers%debug_flag=2
#endif        
     end if   
     
#endif

   end subroutine coarsen_op


      pure subroutine node_write_op(local_node, neighbor_node)
        type(t_node_data), intent(inout)			    :: local_node
        type(t_node_data), intent(in)			    :: neighbor_node
        !do nothing
      end subroutine node_write_op
      
      
      pure subroutine edge_write_op(local_node, neighbor_node)
        type(t_edge_data), intent(inout)			    :: local_node
        type(t_edge_data), intent(in)			    :: neighbor_node
        !do nothing
        local_node%data_pers = neighbor_node%data_pers
        local_node%data_temp = neighbor_node%data_temp
        local_node%rep = neighbor_node%rep        

      end subroutine edge_write_op

    END MODULE
#endif

    
