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
    use SWE_PATCH
    use SWE_dg_predictor
    use SWE_dg_solver
    use SWE_DG_Limiter    

    implicit none

    type num_traversal_data
       logical, DIMENSION(_SWE_PATCH_ORDER_SQUARE) :: dry_cell
    end type num_traversal_data

    type(t_gv_Q)							:: gv_Q
    type(t_lfs_flux)						:: lfs_flux
    
#		define _GT_NAME					t_swe_adaption_traversal
    
#		define _GT_EDGES
#   define _GT_EDGE_MPI_TYPE

#   undef _GT_SKELETON_OP  
#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_PRE_TRAVERSAL_OP			pre_traversal_op
#		define _GT_POST_TRAVERSAL_OP			post_traversal_op

#		define _GT_TRANSFER_OP				transfer_op
#		define _GT_REFINE_OP				refine_op
#		define _GT_COARSEN_OP				coarsen_op
#		define _GT_CELL_TO_EDGE_OP			cell_to_edge_op_dg
#		define _GT_CELL_UPDATE_OP		        cell_update_op_adapt
#		include "SFC_generic_adaptive_traversal.f90"


        subroutine cell_update_op_adapt(traversal, section, element, update1, update2, update3)
          type(t_swe_adaption_traversal), intent(inout)				:: traversal
          type(t_grid_section), intent(inout)							:: section
          type(t_element_base), intent(inout)						:: element
          type(num_cell_update), intent(inout)						:: update1, update2, update3
          integer :: i
          !call updateCellStatus(element%cell%data_pers)
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
        
        subroutine pre_traversal_grid_op(traversal, grid)
          type(t_swe_adaption_traversal), intent(inout)				:: traversal
          type(t_grid), intent(inout)							        :: grid

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
          type(t_state), dimension(_SWE_PATCH_ORDER_SQUARE)						:: Q
          real (kind = GRID_SR), dimension(_SWE_PATCH_ORDER_SQUARE)				:: H, HU, HV, B
          integer :: i

          dest_element%cell%data_pers%H = src_element%cell%data_pers%H
          dest_element%cell%data_pers%HU = src_element%cell%data_pers%HU
          dest_element%cell%data_pers%HV = src_element%cell%data_pers%HV
          dest_element%cell%data_pers%B = src_element%cell%data_pers%B

          do i=1,size(dest_element%edges)
             dest_element%edges(i)%ptr%rep = src_element%edges(i)%ptr%rep
             dest_element%edges(i)%ptr%data_pers = src_element%edges(i)%ptr%data_pers
          end do

          dest_element%cell%data_pers%QFV= src_element%cell%data_pers%QFV
          dest_element%cell%data_pers%QP = src_element%cell%data_pers%QP
          dest_element%cell%data_pers%FP = src_element%cell%data_pers%FP
          dest_element%cell%data_pers%Q  = src_element%cell%data_pers%Q
          
          dest_element%cell%data_pers%troubled=src_element%cell%data_pers%troubled
        end subroutine transfer_op


        subroutine refine_op(traversal, section, src_element, dest_element, refinement_path)
          type(t_swe_adaption_traversal), intent(inout)	           :: traversal
          type(t_grid_section), intent(inout)			   :: section
          type(t_traversal_element), intent(inout)		   :: src_element
          type(t_traversal_element), intent(inout)		   :: dest_element
          integer, dimension(:), intent(in)			   :: refinement_path
          real (kind=GRID_SR), dimension(_SWE_PATCH_ORDER_SQUARE)  :: H_in, HU_in, HV_in, B_in
          real (kind=GRID_SR), dimension(_SWE_PATCH_ORDER_SQUARE)  :: H_out, HU_out, HV_out, B_out
          integer :: i_plotter_type, j, row, col, cell_id
          integer, DIMENSION(_SWE_PATCH_ORDER_SQUARE,2)            :: child
          real (kind = GRID_SR), DIMENSION(2)                      :: r_coords
          logical, DIMENSION(_SWE_PATCH_ORDER_SQUARE)              :: dry_cell_in, dry_cell_out
          integer					           :: i
          
          
          i_plotter_type = src_element%cell%geometry%i_plotter_type
          if(.not.isCoast(src_element%cell%data_pers%troubled)) then
             associate(src => src_element%cell%data_pers, dest => dest_element%cell%data_pers)
               if(size(refinement_path) .ge. 1) then
                  ! decide which child is being computed
                  ! first = left to hypot, second = right to hypot
                  if ( (refinement_path(1) == 1 .and. i_plotter_type > 0) .or.&
                       (refinement_path(1) == 2 .and. i_plotter_type < 0)) then
                     dest%Q%H    = matmul(s_m_inv_ref2,src%Q%H+src%Q%B)
                     dest%Q%p(1) = matmul(s_m_inv_ref2,src%Q%p(1))
                     dest%Q%p(2) = matmul(s_m_inv_ref2,src%Q%p(2))
                  else
                     dest%Q%H    = matmul(s_m_inv_ref1,src%Q%H+src%Q%B)
                     dest%Q%p(1) = matmul(s_m_inv_ref1,src%Q%p(1))
                     dest%Q%p(2) = matmul(s_m_inv_ref1,src%Q%p(2))
                  end if
               end if

               if(size(refinement_path) .ge. 2) then
                  i_plotter_type=i_plotter_type*(-1)
                  if ( (refinement_path(2) == 1 .and. i_plotter_type > 0) .or.&
                       (refinement_path(2) == 2 .and. i_plotter_type < 0)) then
                     dest%Q%H    = matmul(s_m_inv_ref2,dest%Q%H)
                     dest%Q%p(1) = matmul(s_m_inv_ref2,dest%Q%p(1))
                     dest%Q%p(2) = matmul(s_m_inv_ref2,dest%Q%p(2))
                  else
                     dest%Q%H    = matmul(s_m_inv_ref1,dest%Q%H)
                     dest%Q%p(1) = matmul(s_m_inv_ref1,dest%Q%p(1))
                     dest%Q%p(2) = matmul(s_m_inv_ref1,dest%Q%p(2))
                  end if
               end if                  
               
               dest%Q%B = get_bathymetry_at_dg_patch(section, dest_element%t_element_base, section%r_time)
               dest%Q%H   = dest%Q%H-dest%Q%B
               dest%troubled=src%troubled
               call updateCellStatus(dest)
               
               if(isDG(dest%troubled)) then
                  call dg_predictor(dest_element,section%r_dt)
               else
                  call apply_phi(dest%Q(:)%h+dest%Q(:)%b ,dest%h)
                  call apply_phi(dest%Q(:)%p(1)          ,dest%hu)
                  call apply_phi(dest%Q(:)%p(2)          ,dest%hv)
                  call apply_phi(dest%Q(:)%b             ,dest%b)
                  
               end if
             end associate
             
        !Shouldn't happen during runtime   
        else if(isCoast(src_element%cell%data_pers%troubled)) then
          dest_element%cell%data_pers%troubled = src_element%cell%data_pers%troubled
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

          end if
#if defined (_DEBUG)                        
           dest_element%cell%data_pers%debug_flag=-4
#endif           
        end subroutine refine_op

        subroutine coarsen_op(traversal, section, src_element, dest_element, refinement_path)
          type(t_swe_adaption_traversal), intent(inout)  :: traversal
          type(t_grid_section), intent(inout)            :: section
          type(t_traversal_element), intent(inout)       :: src_element
          type(t_traversal_element), intent(inout)       :: dest_element
          integer, dimension(:), intent(in)              :: refinement_path
          integer                                        :: i
          integer                                        :: j,i_plotter_type
          integer, DIMENSION(_SWE_PATCH_ORDER_SQUARE,2)  :: child
          logical                                        :: fv_coarse , dg_coarse , first


          i_plotter_type=dest_element%cell%geometry%i_plotter_type            
          if( refinement_path(1).eq.1 ) then
             dest_element%cell%data_pers%Q%H    = 0.0_GRID_SR
             dest_element%cell%data_pers%Q%p(1) = 0.0_GRID_SR
             dest_element%cell%data_pers%Q%p(2) = 0.0_GRID_SR
             dest_element%cell%data_pers%Q%B    = 0.0_GRID_SR
          endif
          if(.not.isCoast(src_element%cell%data_pers%troubled)) then
             associate(src => src_element%cell%data_pers, dest => dest_element%cell%data_pers)
               if ( (refinement_path(1) == 1 .and. i_plotter_type > 0) .or.&
                    (refinement_path(1) == 2 .and. i_plotter_type < 0)) then
                  !Hpreserving constant water height
                  dest%Q%H    = dest%Q%H    + matmul(s_m_inv_coarsen2,src%Q%H +  src%Q%B)
                  dest%Q%p(1) = dest%Q%p(1) + matmul(s_m_inv_coarsen2,src%Q%p(1))
                  dest%Q%p(2) = dest%Q%p(2) + matmul(s_m_inv_coarsen2,src%Q%p(2))
               else
                  !preserving constant water height
                  dest%Q%H    = dest%Q%H    + matmul(s_m_inv_coarsen1,src%Q%H  + src%Q%B)
                  dest%Q%p(1) = dest%Q%p(1) + matmul(s_m_inv_coarsen1,src%Q%p(1))
                  dest%Q%p(2) = dest%Q%p(2) + matmul(s_m_inv_coarsen1,src%Q%p(2))
               end if
               
               if(refinement_path(1) == 2)then
                  dest%Q%B = &
                       get_bathymetry_at_dg_patch(section, dest_element%t_element_base, section%r_time)

                  dest%Q%H = dest%Q%H-dest%Q%B
                  call updateCellStatus(dest)
                  
                  if(isDG(dest%troubled)) then
                     call dg_predictor(dest_element,section%r_dt)
                  else
                     call apply_phi(dest%Q(:)%h+dest%Q(:)%b ,dest%h)
                     call apply_phi(dest%Q(:)%p(1)          ,dest%hu)
                     call apply_phi(dest%Q(:)%p(2)          ,dest%hv)
                     call apply_phi(dest%Q(:)%b             ,dest%b)
                  end if
               end if
             end associate
          else
             !Shouldnt happen
             print*,"Coarsening coast"
             stop
          end if
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

    
