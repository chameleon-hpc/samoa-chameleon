! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_SWE)
	MODULE SWE_Euler_Timestep
		use SFC_edge_traversal

		use Samoa_swe
		use c_bind_riemannsolvers

#		if defined(_SWE_PATCH)
			use SWE_PATCH
			use SWE_PATCH_Solvers
#               if defined(_SWE_DG)
                        use SWE_DG_Matrices
#		endif
#		endif
#               if defined(_HLLE_FLUX)
                        use SWE_HLLE
#               endif

                        implicit none

        type num_traversal_data
            integer (kind = GRID_DI)			:: i_refinements_issued
        end type num_traversal_data

        interface skeleton_op
            module procedure skeleton_array_op
            module procedure skeleton_scalar_op
        end interface

        interface bnd_skeleton_op
            module procedure bnd_skeleton_array_op
            module procedure bnd_skeleton_scalar_op
        end interface

        PUBLIC cell_to_edge_op
        PUBLIC compute_geoclaw_flux
              
#if defined(_SWE_DG)
#endif
		type(t_gv_Q)							:: gv_Q
		type(t_lfs_flux)						:: lfs_flux

#		define _GT_NAME					t_swe_euler_timestep_traversal

#		define _GT_EDGES

#		define _GT_PRE_TRAVERSAL_OP			pre_traversal_op
#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op
#		define _GT_POST_TRAVERSAL_GRID_OP		post_traversal_grid_op

#		define _GT_CELL_TO_EDGE_OP			cell_to_edge_op

#		define _GT_SKELETON_OP				skeleton_op
#		define _GT_BND_SKELETON_OP			bnd_skeleton_op

#		define _GT_CELL_UPDATE_OP			cell_update_op
#		define _GT_CELL_LAST_TOUCH_OP			cell_last_touch_op
!#		define _GT_CELL_FIRST_TOUCH_OP			cell_first_touch_op
#		define _GT_NODE_WRITE_OP			node_write_op
#		define _GT_EDGE_WRITE_OP			edge_write_op

#		define _GT_NODE_MPI_TYPE

#		include "SFC_generic_traversal_ringbuffer.f90"

        subroutine create_node_mpi_type(mpi_node_type)
            integer, intent(out)            :: mpi_node_type

#           if defined(_MPI)
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
#           endif
        end subroutine

		!*******************************
		!Geometry operators
		!*******************************

        subroutine pre_traversal_grid_op(traversal, grid)
          type(t_swe_euler_timestep_traversal), intent(inout)	:: traversal
          type(t_grid), intent(inout)		                :: grid

          if (cfg%r_max_time > 0.0_SR) then
             grid%r_dt = min(cfg%r_max_time, grid%r_dt)
          end if
          
          if (cfg%r_output_time_step > 0.0_SR) then
             grid%r_dt = min(cfg%r_output_time_step, grid%r_dt)
          end if

#           if defined(_ASAGI)
                !if we are in the earthquake phase, limit the simulation time step by the earthquake time step
                if (grid%r_time < cfg%t_max_eq) then
                    grid%r_dt = min(grid%r_dt, cfg%dt_eq)
                end if
#           endif
                call scatter(grid%r_dt, grid%sections%elements_alloc%r_dt)

		end subroutine

  subroutine post_traversal_grid_op(traversal, grid)
    type(t_swe_euler_timestep_traversal), intent(inout)		:: traversal
    type(t_grid), intent(inout)					:: grid

    grid%r_time = grid%r_time + grid%r_dt
    
    call reduce(traversal%i_refinements_issued, traversal%sections%i_refinements_issued, MPI_SUM, .true.)
    call reduce(grid%r_dt_new, grid%sections%elements_alloc%r_dt_new, MPI_MIN, .true.)
    
    grid%r_dt_new = cfg%courant_number * grid%r_dt_new

    if (rank_MPI == 0) then
       if (cfg%courant_number > grid%r_dt_new / grid%r_dt) then
          _log_write(1, '("WARNING! Time step size was too big. dt (old): ", ES10.3, ", dt (CFL): ", ES10.3, ", maximum courant number: ", F0.3)') grid%r_dt, grid%r_dt_new / cfg%courant_number, grid%r_dt_new / grid%r_dt
       end if
    end if
    
    grid%r_dt = grid%r_dt_new
    grid%r_dt = min(cfg%r_max_time, grid%r_dt)    
    call scatter(grid%r_time, grid%sections%elements_alloc%r_time)
  end subroutine post_traversal_grid_op

  subroutine pre_traversal_op(traversal, section)
    type(t_swe_euler_timestep_traversal), intent(inout)				:: traversal
    type(t_grid_section), intent(inout)						:: section
    
    !this variable will be incremented for each cell with a refinement request
    traversal%i_refinements_issued = 0_GRID_DI
    section%r_dt_new = huge(1.0_SR)
  end subroutine pre_traversal_op

  function cell_to_edge_op(element, edge) result(rep)
    type(t_element_base), intent(in)			:: element
    type(t_edge_data), intent(in)			:: edge
    type(num_cell_rep)					:: rep
    type(t_state), dimension(_SWE_CELL_SIZE)		:: Q
    integer(kind = GRID_SI)				:: i, j, i_edge
    real(kind = GRID_SR), dimension(2, _SWE_EDGE_SIZE)	:: dof_pos
    real(kind = GRID_SR), dimension(2, 3), parameter	:: edge_offsets = reshape([0.0, 0.0, 0.0, 1.0, 1.0, 0.0], [2, 3])
    real(kind = GRID_SR), dimension(2, 3), parameter	:: edge_vectors = reshape([0.0, 1.0, 1.0, -1.0, -1.0, 0.0], [2, 3])
    type(num_cell_data_pers) :: data_temp
    
    
#           if defined(_SWE_PATCH)
                integer                                             :: edge_type !1=left, 2=hypotenuse, 3=right
#if defined (_SWE_DG)           
                data_temp = element%cell%data_pers

!                if(edge%data_pers%troubled) then
                if(element%cell%data_pers%troubled.le.0) then
                   !recover DG from predictor
                   data_temp%Q_DG = element%cell%data_pers%Q_DG_P(1:_SWE_DG_DOFS)
                   call data_temp%convert_dg_to_fv()
                end if
#endif
#           else
                call gv_Q%read(element, Q)
                _log_write(6, '(3X, A)') "swe cell to edge op:"
                _log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "Q in: ", Q
                _log_write(6, '(4X, A, F0.3, 1X, F0.3)') "normal in : ", edge%transform_data%normal
#           endif
                
                
#           if (_SWE_CELL_SIZE > 1)
                i_edge = edge%transform_data%index
                _log_write(6, '(4X, A, I0)') "edge ", i_edge

                forall (i = 1 : _SWE_EDGE_SIZE)
                   dof_pos(:, i) = edge_offsets(:, i_edge) + t_basis_flux_get_dof_coords(i) * edge_vectors(:, i_edge)
                end forall

                call lfs_flux%transform(edge%transform_data, dof_pos(1, :))
                call lfs_flux%transform(edge%transform_data, dof_pos(2, :))

                forall (i = 1 : _SWE_EDGE_SIZE)
                   rep%Q(i)%h = t_basis_Q_eval(dof_pos(:, i), Q%h)
                   rep%Q(i)%p(1) = t_basis_Q_eval(dof_pos(:, i), Q%p(1))
                   rep%Q(i)%p(2) = t_basis_Q_eval(dof_pos(:, i), Q%p(2))
                end forall
#           elif defined(_SWE_PATCH)

                !find out which edge it is comparing its normal with cell normals
                ! obs.: negative plotter_types describe triangles in desired order: left, hypotenuse, right
                associate(cell_edge => ref_plotter_data(- abs(element%cell%geometry%i_plotter_type))%edges, normal => edge%transform_data%normal)
                do i=1,3
                   if ((normal(1) == cell_edge(i)%normal(1) .and. normal(2) == cell_edge(i)%normal(2))   &
                        .or. (normal(1) == -cell_edge(i)%normal(1) .and. normal(2) == -cell_edge(i)%normal(2))) then
                      edge_type = i
                   end if
                end do
              end associate
#if defined(_SWE_DG)
              associate(H => data_temp%H, HU => data_temp%HU, HV => data_temp%HV, B => data_temp%B)
#else
              associate(H => element%cell%data_pers%H, HU => element%cell%data_pers%HU, HV => element%cell%data_pers%HV, B => element%cell%data_pers%B)
#endif
                ! copy boundary values to respective edges
                ! left leg cells go to edge 1
                ! hypotenuse cells go to edge 2
                ! right leg cells go to edge 3

              select case (edge_type)
              case (1) !cells with id i*i+1 (left leg)
                 do i=0, _SWE_PATCH_ORDER - 1
                    rep%H(i+1) = H(i*i + 1)
                    rep%HU(i+1) = HU(i*i + 1)
                    rep%HV(i+1) = HV(i*i + 1)
                    rep%B(i+1) = B(i*i + 1)
                 end do
              case (2) ! hypotenuse
                 do i=1, _SWE_PATCH_ORDER
                    rep%H(i) = H((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*i - 1)
                    rep%HU(i) = HU((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*i - 1)
                    rep%HV(i) = HV((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*i - 1)
                    rep%B(i) = B((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*i - 1)
                 end do
              case (3) !cells with id i*i (right leg)
                 do i=1, _SWE_PATCH_ORDER
                    rep%H(_SWE_PATCH_ORDER + 1 - i) = H(i*i)
                    rep%HU(_SWE_PATCH_ORDER + 1 - i) = HU(i*i)
                    rep%HV(_SWE_PATCH_ORDER + 1 - i) = HV(i*i)
                    rep%B(_SWE_PATCH_ORDER + 1 - i) = B(i*i)
                 end do
              end select
            end associate
#if defined(_SWE_DG)
#endif
#           else
            rep%Q(1) = Q(1)
            _log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "Q out: ", rep%Q
#           endif

          end function

          subroutine skeleton_array_op(traversal, grid, edges, rep1, rep2, update1, update2)
            type(t_swe_euler_timestep_traversal), intent(in)		:: traversal
            type(t_grid_section), intent(in)		                :: grid
            type(t_edge_data), intent(in)			        :: edges(:)
            type(num_cell_rep), intent(in)				:: rep1(:), rep2(:)
            type(num_cell_update), intent(out)				:: update1(:), update2(:)
            integer (kind = GRID_SI)                                    :: i

            do i = 1, size(edges)
                call skeleton_scalar_op(traversal, grid, edges(i), rep1(i), rep2(i), update1(i), update2(i))
            end do
          end subroutine skeleton_array_op

          subroutine skeleton_scalar_op(traversal, grid, edge, rep1, rep2, update1, update2)
            type(t_swe_euler_timestep_traversal), intent(in)		:: traversal
            type(t_grid_section), intent(in)				:: grid
            type(t_edge_data), intent(in)				:: edge
            type(num_cell_rep), intent(in)				:: rep1, rep2
            type(num_cell_update), intent(out)				:: update1, update2
            integer                                                     :: i

            _log_write(6, '(3X, A)') "swe skeleton op:"


#      if defined (_SWE_LF) || defined (_SWE_LF_BATH) || defined (_SWE_LLF) || defined (_SWE_LLF_BATH)
            _log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "Q 1 in: ", rep1%Q
            _log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "Q 2 in: ", rep2%Q
            call compute_lf_flux(edge%transform_data%normal, rep1%Q(1), rep2%Q(1), update1%flux(1), update2%flux(1))
#      elif defined (_SWE_PATCH)
            ! invert values in edges
            ! cells are copied in inverse order because the neighbor 
            ! ghost cells will have a mirrored numbering! See a (poorly-drawn) example:
            !         ___
            !  /\3   1\  |
            ! /  \2   2\ |
            !/____\1   3\|
            do i=1, _SWE_PATCH_ORDER
               update1%H(i) = rep2%H(_SWE_PATCH_ORDER + 1 - i)
               update1%HU(i) = rep2%HU(_SWE_PATCH_ORDER + 1 - i)
               update1%HV(i) = rep2%HV(_SWE_PATCH_ORDER + 1 - i)
               update1%B(i) = rep2%B(_SWE_PATCH_ORDER + 1 - i)


               update2%H(i) = rep1%H(_SWE_PATCH_ORDER + 1 - i)
               update2%HU(i) = rep1%HU(_SWE_PATCH_ORDER + 1 - i)
               update2%HV(i) = rep1%HV(_SWE_PATCH_ORDER + 1 - i)
               update2%B(i) = rep1%B(_SWE_PATCH_ORDER + 1 - i)
            end do
#       else
            call compute_geoclaw_flux(edge%transform_data%normal, rep1%Q(1), rep2%Q(1), update1%flux(1), update2%flux(1))
#   	endif

            _log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "flux 1 out: ", update1%flux
            _log_write(6, '(4X, A, F0.3, 1X, F0.3, 1X, F0.3, 1X, F0.3)') "flux 2 out: ", update2%flux

          end subroutine

          subroutine bnd_skeleton_array_op(traversal, grid, edges, rep, update)
            type(t_swe_euler_timestep_traversal), intent(in)		:: traversal
            type(t_grid_section), intent(in)				:: grid
            type(t_edge_data), intent(in)				:: edges(:)
            type(num_cell_rep), intent(in)				:: rep(:)
            type(num_cell_update), intent(out)				:: update(:)
            integer (kind = GRID_SI)                                    :: i

            do i = 1, size(edges)
               call bnd_skeleton_scalar_op(traversal, grid, edges(i), rep(i), update(i))
            end do
          end subroutine bnd_skeleton_array_op

          subroutine bnd_skeleton_scalar_op(traversal, grid, edge, rep, update)
            type(t_swe_euler_timestep_traversal), intent(in)		:: traversal
            type(t_grid_section), intent(in)			        :: grid
            type(t_edge_data), intent(in)			        :: edge
            type(num_cell_rep), intent(in)				:: rep
            type(num_cell_update), intent(out)				:: update
            type(t_state)						:: bnd_rep
            type(t_update)						:: bnd_flux
            integer 							:: i,j,k,indx,indx_mir
            real(kind=GRID_SR) :: length_flux
            real(kind=GRID_SR) :: normal_normed(2)

            !SLIP: reflect momentum at normal
            !bnd_rep = t_state(rep%Q(1)%h, rep%Q(1)%p - dot_product(rep%Q(1)%p, edge%transform_data%normal) * edge%transform_data%normal, rep%Q(1)%b)

            !NOSLIP: invert momentum (stable)
            !bnd_rep = t_state(rep%Q(1)%h, -rep%Q(1)%p, rep%Q(1)%b)
            normal_normed=(edge%transform_data%normal)/NORM2(edge%transform_data%normal)

            !OUTFLOW: copy values
            bnd_rep = rep%Q(1)


#if defined (_LF_FLUX) || defined (_LF_BATH_FLUX) || defined (_LLF_FLUX) || defined (_LLF_BATH_FLUX)
            call compute_lf_flux(edge%transform_data%normal, rep%Q(1), bnd_rep, update%flux(1), bnd_flux)

#elif defined (_SWE_PATCH)
            ! boundary conditions on ghost cells
            update%H = rep%H
            !                update%HU = rep%HU
            !                update%HV = rep%HV
            update%B = rep%B

            do i=1,_SWE_PATCH_ORDER
               length_flux = (rep%HU(i)* normal_normed(1)+rep%HV(i) * normal_normed(2))
               !reflecting
!               update%HU(i) = rep%HU(i)-2.0_GRID_SR*length_flux*normal_normed(1)
!               update%HV(i) = rep%HV(i)-2.0_GRID_SR*length_flux*normal_normed(2)

               !outflow
               update%HU(i) = rep%HU(i)
               update%HV(i) = rep%HV(i)
            end do

#else

            call compute_geoclaw_flux(edge%transform_data%normal, rep%Q(1), bnd_rep, update%flux(1), bnd_flux)
#endif

          end subroutine bnd_skeleton_scalar_op
          !fv cell update
          subroutine cell_update_op(traversal, section, element, update1, update2, update3)
            type(t_swe_euler_timestep_traversal), intent(inout)		:: traversal
            type(t_grid_section), intent(inout)				:: section
            type(t_element_base), intent(inout)				:: element
            type(num_cell_update), intent(inout)			:: update1, update2, update3
#if defined (_SWE_PATCH)
            integer                                                     :: i, j, ind
            type(num_cell_update)                                       :: tmp !> ghost cells in correct order 
            real(kind = GRID_SR)                                        :: volume, edge_lengths(3), maxWaveSpeed, dQ_max_norm, dt_div_volume

            real(kind = GRID_SR), DIMENSION(_SWE_PATCH_ORDER_SQUARE)     :: dQ_H, dQ_HU, dQ_HV !> deltaQ, used to compute cell updates
            real(kind = GRID_SR), DIMENSION(_SWE_PATCH_SOLVER_CHUNK_SIZE)           :: hL, huL, hvL, bL
            real(kind = GRID_SR), DIMENSION(_SWE_PATCH_SOLVER_CHUNK_SIZE)           :: hR, huR, hvR, bR
            real(kind = GRID_SR), DIMENSION(_SWE_PATCH_SOLVER_CHUNK_SIZE)           :: upd_hL, upd_huL, upd_hvL, upd_hR, upd_huR, upd_hvR
            real(kind = GRID_SR), DIMENSION(_SWE_PATCH_SOLVER_CHUNK_SIZE,2,2)       :: transf
            real(kind= GRID_SR) :: maxWaveSpeeds (_SWE_PATCH_ORDER_SQUARE)

#if !defined (_SWE_USE_PATCH_SOLVER)
            type(t_state), dimension(_SWE_PATCH_SOLVER_CHUNK_SIZE)      :: edges_a, edges_b
            type(t_update)                                              :: update_a, update_b
            real(kind = GRID_SR), dimension(2,3)                        :: normals
#endif
#if defined(_SWE_DG)
            type(t_update),  DIMENSION((_SWE_DG_ORDER+1)**2)		:: flux1, flux2, flux3
            type(t_state), dimension((_SWE_DG_ORDER+1)**2)              :: edge_l, edge_m,edge_r
            real(kind= GRID_SR):: dt,dx,delta(3),max_neighbour(3),min_neighbour(3),data_max_val(3),data_min_val(3)
            integer:: indx,indx_mir,k
            real (kind = GRID_SR), dimension (_SWE_PATCH_ORDER_SQUARE):: H_old, HU_old, HV_old
            logical :: drying,troubled
#endif

            !DIR$ ASSUME_ALIGNED hL: 64
            !DIR$ ASSUME_ALIGNED hR: 64
            !DIR$ ASSUME_ALIGNED huL: 64
            !DIR$ ASSUME_ALIGNED huR: 64
            !DIR$ ASSUME_ALIGNED hvL: 64
            !DIR$ ASSUME_ALIGNED hvR: 64
            !DIR$ ASSUME_ALIGNED bL: 64
            !DIR$ ASSUME_ALIGNED bR: 64
            !DIR$ ASSUME_ALIGNED upd_hL: 64
            !DIR$ ASSUME_ALIGNED upd_hR: 64
            !DIR$ ASSUME_ALIGNED upd_huL: 64
            !DIR$ ASSUME_ALIGNED upd_huR: 64
            !DIR$ ASSUME_ALIGNED upd_hvL: 64
            !DIR$ ASSUME_ALIGNED upd_hvR: 64
            !DIR$ ASSUME_ALIGNED transf: 64

#if !defined(_SWE_USE_PATCH_SOLVER)
            ! using patches, but applying geoclaw solvers on single edges
            ! the normals are only needed in this case.
            ! copy/compute normal vectors
            ! normal for type 2 edges is equal to the 2nd edge's normal
            normals(:,2) = element%edges(2)%ptr%transform_data%normal
            ! normal for types 1 and 3 depend on cell orientation.
            ! notice that normal for type 1 points inwards. That's why it is multiplied by -1.
            if (element%cell%geometry%i_plotter_type < 0) then ! orientation = backward
               normals(:,1) = - element%edges(1)%ptr%transform_data%normal
               normals(:,3) = element%edges(3)%ptr%transform_data%normal
            else ! orientation = forward, so reverse edge order
               normals(:,1) = - element%edges(3)%ptr%transform_data%normal
               normals(:,3) = element%edges(1)%ptr%transform_data%normal
            end if
#endif

            if (element%cell%geometry%i_plotter_type > 0) then ! if orientation = forward, reverse updates
               tmp=update1
               update1=update3
               update3=tmp
            end if

            ! init some variables
            dQ_H = 0.0_GRID_SR
            dQ_HU = 0.0_GRID_SR
            dQ_HV = 0.0_GRID_SR
            maxWaveSpeed = 0.0_GRID_SR

            volume = cfg%scaling * cfg%scaling * element%cell%geometry%get_volume() / (_SWE_PATCH_ORDER_SQUARE)
            dt_div_volume = section%r_dt / volume
            edge_lengths = cfg%scaling * element%cell%geometry%get_edge_sizes() / _SWE_PATCH_ORDER

            associate(data => element%cell%data_pers, geom => SWE_PATCH_geometry)
#if defined(_SWE_DG)
            if(data%troubled.le.0) then
               do i = 1,size(element%edges)
                  if(element%edges(i)%ptr%data_pers%troubled) then
                     element%cell%data_pers%troubled = -4
                  end if
               end do
            end if

            if(data%troubled.ge.1) then
#                               endif                

               ! copy cell values to arrays edges_a and edges_b
               ! obs: cells with id > number of cells are actually ghost cells and come from edges "updates"
               ! see t_SWE_PATCH_geometry for an explanation about ghost cell ordering
               do i=1, _SWE_PATCH_NUM_EDGES, _SWE_PATCH_SOLVER_CHUNK_SIZE ! i -> init of chunk
               ! if this is the last chunk and it is not a full chunk, 
               ! it is necessary to set everything to 0 to avoid using last iteration values
                  if (i + _SWE_PATCH_SOLVER_CHUNK_SIZE -1 > _SWE_PATCH_NUM_EDGES) then
                     hL = 0.0_GRID_SR
                     huL = 0.0_GRID_SR
                     hvL = 0.0_GRID_SR
                     bL = 0.0_GRID_SR
                     hR = 0.0_GRID_SR
                     huR = 0.0_GRID_SR
                     hvR = 0.0_GRID_SR
                     bR = 0.0_GRID_SR
                  end if

                  do j=1, _SWE_PATCH_SOLVER_CHUNK_SIZE ! j -> position inside chunk
                     ind = i + j - 1 ! actual index

                     ! don't go outside array limits
                     if (ind > _SWE_PATCH_NUM_EDGES) then
                        exit
                     end if
                     ! !!!!!!print*,"updates"
                     ! !!!!!!print*,update1%H
                     ! !!!!!!print*,update1%HU
                     ! !!!!!!print*,update1%HV
                     ! !!!!!!print*,update1%B
                     ! !!!!!!print*,update2%H
                     ! !!!!!!print*,update2%HU
                     ! !!!!!!print*,update2%HV
                     ! !!!!!!print*,update2%B
                     ! !!!!!!print*,update3%H
                     ! !!!!!!print*,update3%HU
                     ! !!!!!!print*,update3%HV
                     ! !!!!!!print*,update3%B



                     ! left
                     if (geom%edges_a(ind) <= _SWE_PATCH_ORDER_SQUARE) then
                        hL(j) = data%H(geom%edges_a(ind))
                        huL(j) = data%HU(geom%edges_a(ind))
                        hvL(j) = data%HV(geom%edges_a(ind))
                        bL(j) = data%B(geom%edges_a(ind))
                     else if (geom%edges_a(ind) <= _SWE_PATCH_ORDER_SQUARE + _SWE_PATCH_ORDER) then
                        hL(j) = update1%H(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE)
                        huL(j) = update1%HU(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE)
                        hvL(j) = update1%HV(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE)
                        bL(j) = update1%B(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE)
                     else if (geom%edges_a(ind) <= _SWE_PATCH_ORDER_SQUARE + 2*_SWE_PATCH_ORDER) then
                        hL(j) = update2%H(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                        huL(j) = update2%HU(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                        hvL(j) = update2%HV(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                        bL(j) = update2%B(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                     else 
                        hL(j) = update3%H(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                        huL(j) = update3%HU(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                        hvL(j) = update3%HV(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                        bL(j) = update3%B(geom%edges_a(ind) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                     end if

                     ! right
                     if (geom%edges_b(ind) <= _SWE_PATCH_ORDER_SQUARE) then
                        hR(j) = data%H(geom%edges_b(ind))
                        huR(j) = data%HU(geom%edges_b(ind))
                        hvR(j) = data%HV(geom%edges_b(ind))
                        bR(j) = data%B(geom%edges_b(ind))
                     else if (geom%edges_b(ind) <= _SWE_PATCH_ORDER_SQUARE + _SWE_PATCH_ORDER) then
                        hR(j) = update1%H(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE)
                        huR(j) = update1%HU(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE)
                        hvR(j) = update1%HV(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE)
                        bR(j) = update1%B(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE)
                     else if (geom%edges_b(ind) <= _SWE_PATCH_ORDER_SQUARE + 2*_SWE_PATCH_ORDER) then
                        hR(j) = update2%H(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                        huR(j) = update2%HU(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                        hvR(j) = update2%HV(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                        bR(j) = update2%B(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE - _SWE_PATCH_ORDER)
                     else 
                        hR(j) = update3%H(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                        huR(j) = update3%HU(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                        hvR(j) = update3%HV(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                        bR(j) = update3%B(geom%edges_b(ind) - _SWE_PATCH_ORDER_SQUARE - 2*_SWE_PATCH_ORDER)
                     end if

                     !copy transformation matrices
                     transf(j,:,:) = geom%transform_matrices(geom%edges_orientation(ind),:,:,element%cell%geometry%i_plotter_type)
                  end do

                  ! compute net_updates -> solve Riemann problems within chunk
#                       if defined (_SWE_USE_PATCH_SOLVER)
#                       if defined(_FWAVE_FLUX) || defined(_AUG_RIEMANN_FLUX)
                  call compute_updates_simd(transf, hL, huL, hvL, bL, hR, huR, hvR, bR, upd_hL, upd_huL, upd_hvL, upd_hR, upd_huR, upd_hvR, maxWaveSpeed)
#                       elif defined(_HLLE_FLUX)
                  call compute_updates_hlle_simd(transf, hL, huL, hvL, bL, hR, huR, hvR, bR, upd_hL, upd_huL, upd_hvL, upd_hR, upd_huR, upd_hvR, maxWaveSpeed)
#                       else
                  ! this should never happen -> SCons rules should avoid this before compiling
#                       error "No valid SWE solver for patches/simd implementation has been defined!"
#                       endif
#                       else                    
                  ! using original geoclaw solver
                  maxWaveSpeeds=0
                  do j=1, _SWE_PATCH_SOLVER_CHUNK_SIZE
                     ind = i + j - 1 ! actual index

                     ! don't go outside array limits
                     if (ind > _SWE_PATCH_NUM_EDGES) then
                        exit
                     end if

                     edges_a(j)%h = hL(j)
                     edges_a(j)%p(1) = huL(j)
                     edges_a(j)%p(2) = hvL(j)
                     edges_a(j)%b = bL(j)
                     edges_b(j)%h = hR(j)
                     edges_b(j)%p(1) = huR(j)
                     edges_b(j)%p(2) = hvR(j)
                     edges_b(j)%b = bR(j)
                     call compute_geoclaw_flux(normals(:,geom%edges_orientation(ind)), edges_a(j), edges_b(j), update_a, update_b)
                     upd_hL(j) = update_a%h
                     upd_huL(j) = update_a%p(1)
                     upd_hvL(j) = update_a%p(2)
                     upd_hR(j) = update_b%h
                     upd_huR(j) = update_b%p(1)
                     upd_hvR(j) = update_b%p(2)
                     maxWaveSpeeds(j)=update_a%max_wave_speed
                  end do
#                       endif

                  ! compute dQ
                  do j=1, _SWE_PATCH_SOLVER_CHUNK_SIZE
                     ind = i + j - 1 ! actual index

                     ! don't go outside array limits
                     if (ind > _SWE_PATCH_NUM_EDGES) then
                        exit
                     end if

                     if (geom%edges_a(ind) <= _SWE_PATCH_ORDER_SQUARE) then !ignore ghost cells
                        dQ_H(geom%edges_a(ind)) = dQ_H(geom%edges_a(ind)) + upd_hL(j) * edge_lengths(geom%edges_orientation(ind))
                        dQ_HU(geom%edges_a(ind)) = dQ_HU(geom%edges_a(ind)) + upd_huL(j) * edge_lengths(geom%edges_orientation(ind))
                        dQ_HV(geom%edges_a(ind)) = dQ_HV(geom%edges_a(ind)) + upd_hvL(j) * edge_lengths(geom%edges_orientation(ind))
                     end if
                     if (geom%edges_b(ind) <= _SWE_PATCH_ORDER_SQUARE) then
                        dQ_H(geom%edges_b(ind)) = dQ_H(geom%edges_b(ind)) + upd_hR(j) * edge_lengths(geom%edges_orientation(ind))
                        dQ_HU(geom%edges_b(ind)) = dQ_HU(geom%edges_b(ind)) + upd_huR(j) * edge_lengths(geom%edges_orientation(ind))
                        dQ_HV(geom%edges_b(ind)) = dQ_HV(geom%edges_b(ind)) + upd_hvR(j) * edge_lengths(geom%edges_orientation(ind))
                     end if
                  end do
               end do

               !set refinement condition -> Here I am using the same ones as in the original no-patches implementation, but considering only the max value.
               element%cell%geometry%refinement = 0
               dQ_max_norm = maxval(abs(dQ_H))

               if (element%cell%geometry%i_depth < cfg%i_max_depth .and. dQ_max_norm > 5.0_GRID_SR * cfg%scaling * get_edge_size(cfg%i_max_depth) / _SWE_PATCH_ORDER ) then
                  element%cell%geometry%refinement = 1
                  traversal%i_refinements_issued = traversal%i_refinements_issued + 1_GRID_DI
               else if (element%cell%geometry%i_depth > cfg%i_min_depth .and. dQ_max_norm < 5.0_GRID_SR * cfg%scaling * get_edge_size(cfg%i_max_depth) / (_SWE_PATCH_ORDER * 8.0_SR) ) then
                  element%cell%geometry%refinement = -1
               endif

               dQ_H = dQ_H * (-dt_div_volume)
               dQ_HU = dQ_HU * (-dt_div_volume)
               dQ_HV = dQ_HV * (-dt_div_volume)

               ! if land is flooded, init water height to dry tolerance and
               ! velocity to zero

               where (data%H < data%B + cfg%dry_tolerance .and. dQ_H > 0.0_GRID_SR)
                  data%H = data%B + cfg%dry_tolerance
                  data%HU = 0.0_GRID_SR
                  data%HV = 0.0_GRID_SR
               end where

               ! update unknowns
               data%H = data%H + dQ_H
               data%HU = data%HU + dQ_HU
               data%HV = data%HV + dQ_HV
               

               ! if the water level falls below the dry tolerance, set water level to 0 and velocity to 0
               where (data%H < data%B + cfg%dry_tolerance) 
                  data%H = data%B
                  data%HU = 0.0_GRID_SR
                  data%HV = 0.0_GRID_SR
               end where
               
               maxWaveSpeed=maxval(maxWaveSpeeds)
               section%r_dt_new = min(section%r_dt_new, volume / (edge_lengths(2) * maxWaveSpeed) )
               maxWaveSpeed=0

#if defined (_SWE_DG)

               if(all(data%H - data%B > cfg%dry_tolerance*50.0)) then
                  call data%convert_fv_to_dg()
                  if(all(data%Q_DG%H > cfg%dry_tolerance*50.0)) then
                     data%troubled = -data%troubled
                  else
                     data%troubled = 1
                  end if
               else
                  data%troubled = 1
               end if
            end if

            !consider new dg cells in r_dt
            if(data%troubled .le. 0) then
               do i=1,_SWE_DG_DOFS
                  maxWaveSpeed =  sqrt(g * (data%Q_DG(i)%h)) + maxval(abs(data%Q_DG(i)%p/data%Q_DG(i)%h))
                  section%r_dt_new = min(section%r_dt_new,cfg%scaling*  element%transform_data%custom_data%scaling  / (maxWaveSpeed* (_SWE_DG_ORDER*4.0_GRID_SR +2.0_GRID_SR)))
                  maxWaveSpeed = 0.0
               end do
            end if
#endif
          end associate
#else
          !local variables

          type(t_state)   :: dQ(_SWE_CELL_SIZE)

          call volume_op(element%cell%geometry, traversal%i_refinements_issued, element%cell%geometry%i_depth, &
               element%cell%geometry%refinement, section%r_dt_new, dQ, [update1%flux, update2%flux, update3%flux], section%r_dt)

          !if land is flooded, init water height to dry tolerance and velocity to 0
          if (element%cell%data_pers%Q(1)%h < element%cell%data_pers%Q(1)%b + cfg%dry_tolerance .and. dQ(1)%h > 0.0_GRID_SR) then
             element%cell%data_pers%Q(1)%h = element%cell%data_pers%Q(1)%b + cfg%dry_tolerance
             element%cell%data_pers%Q(1)%p = [0.0_GRID_SR, 0.0_GRID_SR]
             ! ! ! !!!!!!!!!print '("Wetting:", 2(X, F0.0))', cfg%scaling * element%transform_data%custom_data%offset + cfg%offset
          end if

          call gv_Q%add(element, dQ)

          !if the water level falls below the dry tolerance, set water level to 0 and velocity to 0
          if (element%cell%data_pers%Q(1)%h < element%cell%data_pers%Q(1)%b + cfg%dry_tolerance) then
             element%cell%data_pers%Q(1)%h = element%cell%data_pers%Q(1)%b
             element%cell%data_pers%Q(1)%p = [0.0_GRID_SR, 0.0_GRID_SR]
          end if
#          endif
!          element%cell%data_pers%troubled = 1                      
        end subroutine


        subroutine cell_last_touch_op(traversal, section, cell)
          type(t_swe_euler_timestep_traversal), intent(inout)		:: traversal
          type(t_grid_section), intent(inout)				:: section
          type(t_cell_data_ptr), intent(inout)				:: cell
          real(kind = GRID_SR)						:: b_norm

          b_norm = minval(abs(cell%data_pers%Q%h - cell%data_pers%Q%b))

          !refine also on the coasts
          if (cell%geometry%i_depth < cfg%i_max_depth .and. b_norm < 20.0_GRID_SR) then
             !cell%geometry%refinement = 1
             !traversal%i_refinements_issued = traversal%i_refinements_issued + 1_GRID_DI
          else if (b_norm < 100.0_GRID_SR) then
             !cell%geometry%refinement = max(cell%geometry%refinement, 0)
          endif
        end subroutine cell_last_touch_op

        ! subroutine cell_first_touch_op(traversal, section, cell)
        !   type(t_swe_euler_timestep_traversal), intent(inout)		:: traversal
        !   type(t_grid_section), intent(inout)				:: section
        !   type(t_cell_data_ptr), intent(inout)				:: cell
        ! end subroutine cell_first_touch_op


                 !*******************************
                 !Volume and DoF operators
                 !*******************************

                 subroutine volume_op(cell, i_refinements_issued, i_depth, i_refinement, r_dt_new, dQ, fluxes, r_dt)

                         type(fine_triangle), intent(in)				                        :: cell
                         integer (kind = GRID_DI), intent(inout)							    :: i_refinements_issued
                         integer (kind = BYTE), intent(in)							        :: i_depth
                         integer (kind = BYTE), intent(out)							        :: i_refinement
                         real(kind = GRID_SR), intent(inout)								    :: r_dt_new
                         type(t_state), dimension(:), intent(out)						    :: dQ
                         type(t_update), dimension(:), intent(in)						    :: fluxes
                         real(kind = GRID_SR), intent(in)								    :: r_dt

                         real(kind = GRID_SR)											    :: volume, dQ_norm, edge_lengths(3)
                         integer (kind = BYTE)												:: i
                         real (kind = GRID_SR), parameter                                    :: refinement_threshold = 5.0_SR

                         _log_write(6, '(3X, A)') "swe cell update op:"
                         _log_write(6, '(4X, A, 4(X, F0.3))') "edge 1 flux in:", fluxes(1)
                         _log_write(6, '(4X, A, 4(X, F0.3))') "edge 2 flux in:", fluxes(2)
                         _log_write(6, '(4X, A, 4(X, F0.3))') "edge 3 flux in:", fluxes(3)

                         volume = cfg%scaling * cfg%scaling * cell%get_volume()
                         edge_lengths = cfg%scaling * cell%get_edge_sizes()

                         dQ%h = dot_product(edge_lengths, fluxes%h)
                         dQ%p(1) = dot_product(edge_lengths, fluxes%p(1))
                         dQ%p(2) = dot_product(edge_lengths, fluxes%p(2))
                         dQ%b = 0.0_GRID_SR

                         !set refinement condition

                         i_refinement = 0
                         dQ_norm = abs(dQ(1)%h)

                         if (i_depth < cfg%i_max_depth .and. dQ_norm > refinement_threshold * cfg%scaling * get_edge_size(cfg%i_max_depth)) then
                                 i_refinement = 1
                                 i_refinements_issued = i_refinements_issued + 1_GRID_DI
                         else if (i_depth > cfg%i_min_depth .and. dQ_norm < refinement_threshold * cfg%scaling * get_edge_size(cfg%i_max_depth) / 8.0_SR) then
                                 i_refinement = -1
                         endif

             !This will cause a division by zero if the wave speeds are 0.
             !Bue to the min operator, the error will not affect the time step.
             !r_dt_new = min(r_dt_new, volume / dot_product(edge_lengths, fluxes%max_wave_speed))
             r_dt_new = min(r_dt_new, volume / (edge_lengths(2) * maxval(fluxes%max_wave_speed)))


             do i = 1, _SWE_CELL_SIZE
                 dQ(i)%t_dof_state = dQ(i)%t_dof_state * (-r_dt / volume)
             end do

                         _log_write(6, '(4X, A, 4(X, F0.3))') "dQ out: ", dQ
                 end subroutine volume_op

                 !> Lax Friedrichs flux. Depending on compiler flags, the function implements
                 !> the global or local variant with or without bathymetry
                 subroutine compute_lf_flux(normal, QL, QR, fluxL, fluxR)
                         type(t_state), intent(in)							:: QL, QR
                         type(t_update), intent(out) 						:: fluxL, fluxR
                         real(kind = GRID_SR), intent(in)		            :: normal(2)

                         real(kind = GRID_SR)								:: vL, vR, hL, hR, alpha

#           if defined(_LF_BATH_FLUX) || defined(_LLF_BATH_FLUX)
                 if (QL%h - QL%b < cfg%dry_tolerance .or. QR%h - QR%b < cfg%dry_tolerance) then
                     hL = 0.0_SR; hR = 0.0_SR
                     vL = 0.0_SR; vR = 0.0_SR

                     fluxL%max_wave_speed = 0.0_SR; fluxR%max_wave_speed = 0.0_SR
                     fluxL%h = 0.0_SR; fluxR%h = 0.0_SR
                     fluxL%p = 0.0_SR; fluxR%p = 0.0_SR

                     !This boundary treatment assumes a wall condition.
                     !For the mass flux, we choose pR := -hL * vL, vR = 0 (walls are immovable), hence hR must be infinite.
                     !For the momentum flux we choose hR := 0, vR := 0 (implying there is no hydrostatic pressure), bR := bL + hL (there is a wall to the right)

                     if (QL%h - QL%b < cfg%dry_tolerance .and. QR%h - QR%b < cfg%dry_tolerance) then
                     else if (QL%h - QL%b < cfg%dry_tolerance) then
                         hR = max(QR%h - QR%b, 0.0_SR)
                         vR = dot_product(normal, QR%p / (QR%h - QR%b))
                         fluxR%max_wave_speed = sqrt(g * hR) + abs(vR)
                         fluxR%p = -0.5_SR * vR * QR%p - 0.5_GRID_SR * g * hR * hR * normal + 0.5_SR * fluxR%max_wave_speed * QR%p
                     else if (QR%h - QR%b < cfg%dry_tolerance) then
                         hL = max(QL%h - QL%b, 0.0_SR)
                         vL = dot_product(normal, QL%p / (QL%h - QL%b))
                         fluxL%max_wave_speed = sqrt(g * hL) + abs(vL)
                         fluxL%p = 0.5_SR * vL * QL%p + 0.5_GRID_SR * g * hL * hL * normal + 0.5_SR * fluxL%max_wave_speed * QL%p
                     end if

                     return
                 end if

                 hL = max(QL%h - QL%b, 0.0_SR)
                 hR = max(QR%h - QR%b, 0.0_SR)

                 vL = dot_product(normal, QL%p / (QL%h - QL%b))
                 vR = dot_product(normal, QR%p / (QR%h - QR%b))

                 fluxL%max_wave_speed = sqrt(g * hL) + abs(vL)
                 fluxR%max_wave_speed = sqrt(g * hR) + abs(vR)

#               if defined(_LLF_BATH_FLUX)
                     alpha = max(fluxL%max_wave_speed, fluxR%max_wave_speed)
#               else
                     alpha = 100.0_GRID_SR
#               endif

                 !Except for the diffusion term, the mass flux is the standard LF flux
                 fluxL%h = 0.5_GRID_SR * (hL * vL + hR * vR + alpha * (QL%h - QR%h))
                 fluxR%h = -fluxL%h

                 !The base momentum flux is similar to the standard LF flux.
                 fluxL%p = 0.5_GRID_SR * (QL%p * vL + QR%p * vR + 0.5_GRID_SR * g * (hL * hL + hR * hR) * normal + alpha * (QL%p - QR%p))
                 fluxR%p = -fluxL%p

                 !The source term $\Delta x \ \Psi = $-1/2 g \ \frac{1}{2} \ (h_l + h_r) \ (b_r - b_l)$ [LeVeque] is added on both sides with a weight of 1/2.
                 !This factor ensures that the method is well-balanced.
                 fluxL%p = fluxL%p + 0.25_SR * g * (hL + hR) * (QR%b - QL%b) * normal
                 fluxR%p = fluxR%p + 0.25_SR * g * (hL + hR) * (QR%b - QL%b) * normal
#           elif defined(_LF_FLUX) || defined(_LLF_FLUX)
                 real(kind = GRID_SR), parameter					:: b = 0.0_GRID_SR     !default constant bathymetry

                 hL = max(QL%h - b, 0.0_SR)
                 hR = max(QR%h - b, 0.0_SR)

                 vL = dot_product(normal, QL%p / (QL%h - b))
                 vR = dot_product(normal, QR%p / (QR%h - b))

                 fluxL%max_wave_speed = sqrt(g * (QL%h - b)) + abs(vL)
                 fluxR%max_wave_speed = sqrt(g * (QR%h - b)) + abs(vR)

#               if defined(_LLF_FLUX)
                     alpha = max(fluxL%max_wave_speed, fluxR%max_wave_speed)
#               else
                     alpha = 100.0_GRID_SR
#               endif


                 fluxL%h = 0.5_GRID_SR * (vL * hL + vR * hR + alpha * (QL%h - QR%h))
                 fluxR%h = -fluxL%h

                 fluxL%p = 0.5_GRID_SR * (vL * vL * hL + vR * vR * hR + 0.5_GRID_SR * g * (hL * hL + hR * hR) * normal + alpha * (QL%p - QR%p))
                 fluxR%p = -fluxL%p

#           endif
                 end subroutine


                 !> Augmented Riemann solver from geoclaw
                 subroutine compute_geoclaw_flux(normal, QL, QR, fluxL, fluxR)
                         type(t_state), intent(in)           :: QL, QR
                         type(t_update), intent(out)         :: fluxL, fluxR
                         real(kind = GRID_SR), intent(in)    :: normal(2)

                         real(kind = GRID_SR)				:: transform_matrix(2, 2)
                         real(kind = GRID_SR)			    :: net_updatesL(3), net_updatesR(3), max_wave_speed
                         real(kind = GRID_SR)                :: pL(2), pR(2), hL, hR, bL, bR



                         transform_matrix(1, :) = normal
                         transform_matrix(2, :) = [-normal(2), normal(1)]

                         pL = matmul(transform_matrix, QL%p)
                         pR = matmul(transform_matrix, QR%p)

                         hL = QL%h-QL%b
                         hR = QR%h-QR%b

                         bL = QL%b
                         bR = QR%b

#           if defined(_FWAVE_FLUX)
                 call c_bind_geoclaw_solver(GEOCLAW_FWAVE, 1, 3, hL, hR, pL(1), pR(1), pL(2), pR(2), bL, bR, real(cfg%dry_tolerance, GRID_SR), g, net_updatesL, net_updatesR, max_wave_speed)
#           elif defined(_AUG_RIEMANN_FLUX)
                 call c_bind_geoclaw_solver(GEOCLAW_AUG_RIEMANN, 1, 3, hL, hR, pL(1), pR(1), pL(2), pR(2), bL, bR, real(cfg%dry_tolerance, GRID_SR), g, net_updatesL, net_updatesR, max_wave_speed)
#           elif defined(_HLLE_FLUX)
                 call compute_updates_hlle_single(hL, hR, pL(1), pR(1), pL(2), pR(2), bL, bR, net_updatesL, net_updatesR, max_wave_speed)
#           endif

                 fluxL%h = net_updatesL(1)
                 fluxL%p = matmul(net_updatesL(2:3), transform_matrix)
                 fluxL%max_wave_speed = max_wave_speed

                 fluxR%h = net_updatesR(1)
                 fluxR%p = matmul(net_updatesR(2:3), transform_matrix)
                 fluxR%max_wave_speed = max_wave_speed
                 end subroutine

         pure subroutine node_write_op(local_node, neighbor_node)
             type(t_node_data), intent(inout)			    :: local_node
             type(t_node_data), intent(in)				    :: neighbor_node

             !do nothing
         end subroutine


         pure subroutine edge_write_op(local_node, neighbor_node)
             type(t_edge_data), intent(inout)			    :: local_node
             type(t_edge_data), intent(in)				    :: neighbor_node

             !do nothing
         end subroutine
         END MODULE

# if defined(_SWE_DG)

  MODULE SWE_dg_predictor
                 use SFC_edge_traversal

                 use Samoa_swe
                 use c_bind_riemannsolvers

#		if defined(_SWE_PATCH)
                         use SWE_PATCH
                         use SWE_PATCH_Solvers
#               if defined(_SWE_DG)
                         use SWE_DG_Matrices
#		endif
#		endif

#               if defined(_HLLE_FLUX)
                         use SWE_HLLE
#               endif

 implicit none
         type num_traversal_data
         end type

         interface skeleton_op
              module procedure skeleton_array_op
              module procedure skeleton_scalar_op
          end interface

          interface bnd_skeleton_op
              module procedure bnd_skeleton_array_op
              module procedure bnd_skeleton_scalar_op
          end interface


          public dg_predictor

#		define _GT_NAME						t_swe_dg_predictor_traversal
#		define _GT_EDGES

#		define _GT_PRE_TRAVERSAL_OP				pre_traversal_op
#		define _GT_PRE_TRAVERSAL_GRID_OP		        pre_traversal_grid_op


#		define _GT_SKELETON_OP					skeleton_op
#		define _GT_BND_SKELETON_OP				bnd_skeleton_op
#		define _GT_ELEMENT_OP                                   element_op
#		define _GT_CELL_UPDATE_OP				cell_update_op

#		include "SFC_generic_traversal_ringbuffer.f90"


   subroutine pre_traversal_grid_op(traversal, grid)
     type(t_swe_dg_predictor_traversal), intent(inout)		:: traversal
     type(t_grid), intent(inout)							    :: grid
     
     if (cfg%r_max_time > 0.0_SR) then
        grid%r_dt = min(cfg%r_max_time, grid%r_dt)
     end if
     
     if (cfg%r_output_time_step > 0.0_SR) then
        grid%r_dt = min(cfg%r_output_time_step, grid%r_dt)
     end if
     
#           if defined(_ASAGI)
     !if we are in the earthquake phase, limit the simulation time step by the earthquake time step
     if (grid%r_time < cfg%t_max_eq) then
        grid%r_dt = min(grid%r_dt, cfg%dt_eq)
     end if
#           endif
     call scatter(grid%r_dt, grid%sections%elements_alloc%r_dt)
     
   end subroutine pre_traversal_grid_op


   subroutine pre_traversal_op(traversal, section)
     type(t_swe_dg_predictor_traversal), intent(inout)	:: traversal
     type(t_grid_section), intent(inout)		:: section
   end subroutine pre_traversal_op

   subroutine dg_predictor(cell,dt,dx_in)
     class(t_cell_data_ptr), intent(in)		:: cell
     integer                                    :: i,j
     real(kind=GRID_SR) ,intent(in)             :: dt
     real(kind=GRID_SR),intent(in)              :: dx_in
     real(kind=GRID_SR)                         :: dx
     real(kind=GRID_SR)                         :: q_0(_SWE_DG_DOFS,3)
     real(kind=GRID_SR)                         :: q_i(_SWE_DG_DOFS*(_SWE_DG_ORDER+1),3)
     real(kind=GRID_SR)                         :: q_temp(_SWE_DG_DOFS*_SWE_DG_ORDER,3)
     real(kind=GRID_SR),Dimension(_SWE_DG_DOFS*(_SWE_DG_ORDER+1),3)  :: f1_hat,f1_s,f2_hat,f2_s,S_s
     real(kind=GRID_SR),Dimension(_SWE_DG_DOFS*_SWE_DG_ORDER,3)      :: source,volume_flux1,volume_flux2,source_temp,orig_source
     real(kind=GRID_SR) :: epsilon
     real(kind=GRID_SR),Dimension(2,2) :: jacobian,jacobian_inv
#if defined(_SWE_DG_NODAL)
     real(kind=GRID_SR),Dimension(_SWE_DG_DOFS*(_SWE_DG_ORDER+1)) :: H_x,H_y,H_x_temp,H_y_temp,b_x,b_y
     real(kind=GRID_SR),Dimension(_SWE_DG_DOFS*(_SWE_DG_ORDER+1)) :: b
#else
     real(kind=GRID_SR),Dimension(size(st_gl_node_vals,1))   :: b_x,b_y
     real(kind=GRID_SR),Dimension(size(st_gl_node_vals,1),3) :: q_v
     integer :: st_nodes = size(st_gl_node_vals,1)
     integer :: s_nodes = size(s_der_x_gl_node_vals,1)
#endif

     associate(Q_DG =>cell%data_pers%Q_DG ,Q_DG_P =>cell%data_pers%Q_DG_P)

     jacobian=ref_plotter_data(abs(cell%geometry%i_plotter_type))%jacobian
     jacobian= jacobian/NORM2(jacobian(1,:))

     jacobian_inv=ref_plotter_data(abs(cell%geometry%i_plotter_type))%jacobian_inv
     jacobian_inv= jacobian_inv/NORM2(jacobian_inv(1,:))

     call cell%data_pers%dofs_to_vec_dg(q_0)

     q_i(1:_SWE_DG_DOFS,:) = q_0
     b(1:_SWE_DG_DOFS) = Q_DG(:)%B

     do i=1,_SWE_DG_ORDER
        q_i(1+i*_SWE_DG_DOFS:(i+1)*_SWE_DG_DOFS,:) = q_0
        b(1+i*_SWE_DG_DOFS:(i+1)*_SWE_DG_DOFS) = Q_DG(:)%B
     end do
     !!!!!!!print*,q_i
     dx=dx_in*cfg%scaling

     epsilon=1.0_GRID_SR                          

     i=0
     do while(epsilon > 1.0e-13_GRID_SR)
        i=i+1

#if defined(_SWE_DG_NODAL)
        !print*
        !print*,q_0(:,1)
        !print*,q_0(:,2)
        !print*,q_0(:,3)
        !print*

        f1_hat = f1(q_i,(_SWE_DG_ORDER+1)*_SWE_DG_DOFS)
        f2_hat = f2(q_i,(_SWE_DG_ORDER+1)*_SWE_DG_DOFS)

        do j=0,(_SWE_DG_ORDER)
           H_x(_SWE_DG_DOFS * j + 1: _SWE_DG_DOFS * (j+1))= matmul(basis_der_x,q_i(_SWE_DG_DOFS * j + 1: _SWE_DG_DOFS * (j+1),1)+b(_SWE_DG_DOFS * j + 1: _SWE_DG_DOFS * (j+1)))
           H_y(_SWE_DG_DOFS * j + 1: _SWE_DG_DOFS * (j+1))= matmul(basis_der_y,q_i(_SWE_DG_DOFS * j + 1: _SWE_DG_DOFS * (j+1),1)+b(_SWE_DG_DOFS * j + 1: _SWE_DG_DOFS * (j+1)))
        end do

        ! where(q_i(:,1).ne.0 .and. abs(H_x/q_i(:,1))< 1.0e-14)
        !    H_x = 0.0_GRID_SR
        ! end where

        ! where(q_i(:,1).ne.0 .and. abs(H_y/q_i(:,1))< 1.0e-14)
        !    H_y = 0.0_GRID_SR
        ! end where

        !print*,"ders"
        !print*,H_x
        !print*,H_y


        source_temp=0
        source_temp(:,2)=-matmul(st_m,(g*q_i(:,1)*H_x))
        source_temp(:,3)=-matmul(st_m,(g*q_i(:,1)*H_y))

        source_temp(:,2)=source_temp(:,2)+matmul(st_k_x,(0.5_GRID_SR*g*q_i(:,1)**2))
        source_temp(:,3)=source_temp(:,3)+matmul(st_k_y,(0.5_GRID_SR*g*q_i(:,1)**2))

        source=0
        source(:,2) = jacobian(1,1) * source_temp(:,2) + jacobian(1,2) * source_temp(:,3)
        source(:,3) = jacobian(2,1) * source_temp(:,2) + jacobian(2,2) * source_temp(:,3)

        b_x=Q_DG_P(:)%b_x
        b_y=Q_DG_P(:)%b_y
        S_s =S(q_i,b_x,b_y,_SWE_DG_DOFS*(_SWE_DG_ORDER+1))

        orig_source=matmul(st_m,S_s)*dt/dx

#else
        q_v=matmul(st_gl_node_vals,q_i)
        f1_hat=matmul(st_gl_weights,f1(q_v,st_nodes))
        f2_hat=matmul(st_gl_weights,f2(q_v,st_nodes))

        b_x=cell%data_pers%b_x
        b_y=cell%data_pers%b_y

        S_s=matmul(st_gl_weights,S(q_v(:,1),b_x,b_y,st_nodes))

        call lusolve(st_m_lu,_SWE_DG_DOFS*(_SWE_DG_ORDER+1),st_m_lu_pivot,f1_hat(:,1))
        call lusolve(st_m_lu,_SWE_DG_DOFS*(_SWE_DG_ORDER+1),st_m_lu_pivot,f1_hat(:,2))
        call lusolve(st_m_lu,_SWE_DG_DOFS*(_SWE_DG_ORDER+1),st_m_lu_pivot,f1_hat(:,3))
        call lusolve(st_m_lu,_SWE_DG_DOFS*(_SWE_DG_ORDER+1),st_m_lu_pivot,f2_hat(:,1))
        call lusolve(st_m_lu,_SWE_DG_DOFS*(_SWE_DG_ORDER+1),st_m_lu_pivot,f2_hat(:,2))
        call lusolve(st_m_lu,_SWE_DG_DOFS*(_SWE_DG_ORDER+1),st_m_lu_pivot,f2_hat(:,3))

        source = S_s(_SWE_DG_DOFS+1:,:)

#endif

        f1_s=(jacobian_inv(1,1)*f1_hat+jacobian_inv(1,2)*f2_hat)
        f2_s=(jacobian_inv(2,1)*f1_hat+jacobian_inv(2,2)*f2_hat)
        ! !!!!!print*
        ! !!!!!print*,f1_s
        source  = source *dt/dx

        volume_flux1 =matmul(st_k_x,f1_s) *dt/dx
        volume_flux2 =matmul(st_k_y,f2_s) *dt/dx

        !print*
        !print*,abs(cell%geometry%i_plotter_type)
        !print*
        !print*,"orig source"
        !print*,orig_source(:,1)
        !print*,orig_source(:,2)
        !print*,orig_source(:,3)
        !print*
        !print*,"Pred source"
        !print*,source(:,1)
        !print*,source(:,2)
        !print*,source(:,3)
        !print*,"Pred flux"
        !print*, - volume_flux1(:,1) - volume_flux2(:,1)
        !print*, - volume_flux1(:,2) - volume_flux2(:,2)
        !print*, - volume_flux1(:,3) - volume_flux2(:,3)
        !print*
!        source=orig_source
        !print*,"sum"
        !print*, source(:,1) - volume_flux1(:,1) - volume_flux2(:,1)
        !print*, source(:,2) - volume_flux1(:,2) - volume_flux2(:,2)
        !print*, source(:,3) - volume_flux1(:,3) - volume_flux2(:,3)

        q_temp= source &
             - matmul(st_w_k_t_1_0,q_i(1:_SWE_DG_DOFS,:))&
             - volume_flux1 &
             - volume_flux2

        !!!!!!!print*,"q_temp"
        call lusolve(st_w_k_t_1_1_lu,_SWE_DG_ORDER*_SWE_DG_DOFS,st_w_k_t_1_1_lu_pivot ,q_temp(:,1))
        call lusolve(st_w_k_t_1_1_lu,_SWE_DG_ORDER*_SWE_DG_DOFS,st_w_k_t_1_1_lu_pivot ,q_temp(:,2))
        call lusolve(st_w_k_t_1_1_lu,_SWE_DG_ORDER*_SWE_DG_DOFS,st_w_k_t_1_1_lu_pivot ,q_temp(:,3))

        epsilon=0.0_GRID_SR                          
        do j=1,(_SWE_DG_ORDER)*_SWE_DG_DOFS
           if(.not.all(q_i(_SWE_DG_DOFS+j,:).eq.0)) then
              epsilon=max(epsilon,NORM2(q_temp(j,:)-q_i(_SWE_DG_DOFS+j,:))/NORM2(q_i(_SWE_DG_DOFS+j,:)))
           else if(.not.all(q_temp(j,:).eq. 0)) then
              epsilon=max(epsilon,1.0_GRID_SR)
           end if
        end do

        if(i > 100) then                           
           !!!print*,"predictor not converging"
           stop
        end if
        
        q_i(1+_SWE_DG_DOFS:(_SWE_DG_ORDER+1)*_SWE_DG_DOFS,:) = q_temp

        !print*
        !print*,q_i(:,1)
        !print*,q_i(:,2)
        !print*,q_i(:,3)

     end do

     ! q_i(:,2) = q_i(:,2) *q_i(:,1)
     ! ! q_i(:,3) = q_i(:,3) *q_i(:,1)
     ! do i=1,_SWE_DG_DOFS
     !    if(q_0(i,1).ne.0.0_GRID_SR)then
     !    do j=0,_SWE_DG_ORDER
     !       if(abs(q_i(i+j*_SWE_DG_DOFS,1)-q_0(i,1))/q_0(i,1) < 1.0e-14)then
     !          q_i(i+j*_SWE_DG_DOFS,1)=q_0(i,1)
     !       end if
     !    end do
     !    else if(abs(q_i(i+j*(_SWE_DG_DOFS),1)) < 1.0e-14)then
     !          q_i(i+j*_SWE_DG_DOFS,1)=q_0(i,1)
     !       end if
     ! end do

     call cell%data_pers%vec_to_dofs_dg_p(q_i)
   end  associate
 end subroutine dg_predictor

 subroutine element_op(traversal, section, element)
   type(t_swe_dg_predictor_traversal), intent(inout)				:: traversal
   type(t_grid_section), intent(inout)						:: section
   type(t_element_base), intent(inout)						:: element
   integer :: i


   if(element%cell%data_pers%troubled.le.0) then
     !print*,"predict"
      call dg_predictor(element%cell,section%r_dt,element%cell%geometry%get_scaling())
     !print*,"predict done"
   ! else
   !    element%cell%data_pers%Q_DG_P(:)%H=0
   !    element%cell%data_pers%Q_DG_P(:)%p(1)=0
   !    element%cell%data_pers%Q_DG_P(:)%p(2)=0
   !    element%cell%data_pers%Q_DG_P(1:_SWE_DG_DOFS) = element%cell%data_pers%Q_DG
   end if
 end subroutine element_op

 function f1(q,N)
   integer :: N
   real(kind=GRID_SR) ,intent(in) ::q(N,3)
   real(kind=GRID_SR)             ::f1(N,3)
   integer                        ::i

   f1(:,1) = q(:,2)
   do i=1,N
      if(q(i,1).eq.0) then
         f1(i,1) = 0
         f1(i,2) = 0
         f1(i,3) = 0
      else
         f1(i,2) = q(i,2)**2/q(i,1) + 0.5_GRID_SR * g * q(i,1)**2
         f1(i,3) = q(i,2)*q(i,3)/q(i,1)
      end if
   end do

 end function f1

 function f2(q,N)
   integer :: N
   real(kind=GRID_SR) ,intent(in) ::q(N,3)
   real(kind=GRID_SR)             ::f2(N,3)
   integer                        ::i

   f2(:,1) = q(:,3)
   do i=1,N
      if(q(i,1).eq.0) then
         f2(i,1) = 0
         f2(i,2) = 0
         f2(i,3) = 0
      else
         f2(i,2) = q(i,2)*q(i,3)/q(i,1)
         f2(i,3) = q(i,3)**2/q(i,1) + 0.5_GRID_SR * g * q(i,1)**2
      end if
   end do

 end function f2

 !source
 function S(q,b_x,b_y,N)
   integer ::N,i
   real(kind=GRID_SR) ,intent(in) ::q(N,3), b_x(N), b_y(N)
   real(kind=GRID_SR)             ::S(N,3)

   S(:,1) = 0
   do i=1,N
      S(i,2)=-q(i,1)*(b_x(i))*g
      S(i,3)=-q(i,1)*(b_y(i))*g 

   end do

 end function S

 function cell_to_edge_op(element, edge) result(rep)
   type(t_element_base), intent(in)						:: element
   type(t_edge_data), intent(in)						    :: edge
   type(num_cell_rep)										:: rep

 end function cell_to_edge_op

 subroutine skeleton_array_op(traversal, grid, edges, rep1, rep2, update1, update2)
   type(t_swe_dg_predictor_traversal), intent(in)				:: traversal
   type(t_grid_section), intent(in)							    :: grid
   type(t_edge_data), intent(in)								    :: edges(:)
   type(num_cell_rep), intent(in)									:: rep1(:), rep2(:)
   type(num_cell_update), intent(out)								:: update1(:), update2(:)
   integer (kind = GRID_SI)                                        :: i

   do i = 1, size(edges)
      call skeleton_scalar_op(traversal, grid, edges(i), rep1(i), rep2(i), update1(i), update2(i))
   end do

 end subroutine skeleton_array_op

 subroutine skeleton_scalar_op(traversal, grid, edge, rep1, rep2, update1, update2)
   type(t_swe_dg_predictor_traversal), intent(in)				:: traversal
   type(t_grid_section), intent(in)							    :: grid
   type(t_edge_data), intent(in)								    :: edge
   type(num_cell_rep), intent(in)									:: rep1, rep2
   type(num_cell_update), intent(out)								:: update1, update2
   !!!!!!!!!print*,rep1%H

 end subroutine skeleton_scalar_op

 subroutine bnd_skeleton_array_op(traversal, grid, edges, rep, update)
   type(t_swe_dg_predictor_traversal), intent(in)				:: traversal
   type(t_grid_section), intent(in)							    :: grid
   type(t_edge_data), intent(in)								    :: edges(:)
   type(num_cell_rep), intent(in)									:: rep(:)
   type(num_cell_update), intent(out)								:: update(:)
   integer :: i

   do i = 1, size(edges)
      call bnd_skeleton_scalar_op(traversal, grid, edges(i), rep(i), update(i))
   end do


 end subroutine bnd_skeleton_array_op

 subroutine bnd_skeleton_scalar_op(traversal, grid, edge, rep, update)
   type(t_swe_dg_predictor_traversal), intent(in)				:: traversal
   type(t_grid_section), intent(in)							    :: grid
   type(t_edge_data), intent(in)								    :: edge
   type(num_cell_rep), intent(in)									:: rep
   type(num_cell_update), intent(out)								:: update

 end subroutine bnd_skeleton_scalar_op


 subroutine cell_update_op(traversal, section, element, update1, update2, update3)
   type(t_swe_dg_predictor_traversal), intent(inout)				:: traversal
   type(t_grid_section), intent(inout)							:: section
   type(t_element_base), intent(inout)						:: element
   type(num_cell_update), intent(inout)						:: update1, update2, update3
   integer :: i

   element%cell%data_pers%troubled_old = element%cell%data_pers%troubled
 end subroutine cell_update_op

END MODULE SWE_DG_predictor



MODULE SWE_DG_timestep
  use SFC_edge_traversal
  
  use Samoa_swe
  use c_bind_riemannsolvers
  
  use SWE_DG_Matrices
  
#               if defined(_HLLE_FLUX)
  use SWE_HLLE
#               endif
  

 implicit none
 type num_traversal_data
    integer (kind = GRID_DI)			:: i_refinements_issued
 end type num_traversal_data
 
 interface skeleton_op_dg
    module procedure skeleton_array_op_dg
    module procedure skeleton_scalar_op_dg
 end interface skeleton_op_dg
 
 interface bnd_skeleton_op_dg
    module procedure bnd_skeleton_array_op_dg
    module procedure bnd_skeleton_scalar_op_dg
 end interface bnd_skeleton_op_dg
 
 interface edge_first_touch_op_dg
    module procedure edge_first_touch_array_op_dg
    module procedure edge_first_touch_scalar_op_dg
 end interface edge_first_touch_op_dg
 
 
#		define _GT_NAME							t_swe_dg_timestep_traversal
 
#		define _GT_EDGES
 
#		define _GT_PRE_TRAVERSAL_OP				pre_traversal_op_dg
#		define _GT_PRE_TRAVERSAL_GRID_OP		pre_traversal_grid_op_dg

#		define _GT_SKELETON_OP					skeleton_op_dg
#		define _GT_BND_SKELETON_OP				bnd_skeleton_op_dg
#		define _GT_ELEMENT_OP                                   element_op_dg

#		define _GT_CELL_UPDATE_OP				cell_update_op_dg
#               undef _GT_CELL_TO_EDGE_OP
#		define _GT_CELL_TO_EDGE_OP				cell_to_edge_op_dg

#               define _GT_EDGE_FIRST_TOUCH_OP                          edge_first_touch_op_dg

 public cell_to_edge_op_dg

#		include "SFC_generic_traversal_ringbuffer.f90"

 subroutine edge_first_touch_scalar_op_dg(traversal, section, edge)
   type(t_swe_dg_timestep_traversal)  ::traversal
   type(t_grid_section), intent(inout)::section
   type(t_edge_data), intent(inout)   ::edge
   
   edge%data_pers%troubled=.false.
 end subroutine edge_first_touch_scalar_op_dg
 
 subroutine edge_first_touch_array_op_dg(traversal, section, edge)
   type(t_swe_dg_timestep_traversal)  ::traversal
   type(t_grid_section), intent(inout)::section
   type(t_edge_data), intent(inout)   ::edge(:)
   edge%data_pers%troubled=.false.
 end subroutine edge_first_touch_array_op_dg



   subroutine pre_traversal_grid_op_dg(traversal, grid)
     type(t_swe_dg_timestep_traversal), intent(inout)		:: traversal
     type(t_grid), intent(inout)							    :: grid
   end subroutine pre_traversal_grid_op_dg


   subroutine pre_traversal_op_dg(traversal, section)
     type(t_swe_dg_timestep_traversal), intent(inout)				:: traversal
     type(t_grid_section), intent(inout)							:: section
   end subroutine pre_traversal_op_dg

   subroutine element_op_dg(traversal, section, element)
     type(t_swe_dg_timestep_traversal), intent(inout)				:: traversal
     type(t_grid_section), intent(inout)							:: section
     type(t_element_base), intent(inout)						:: element
   end subroutine element_op_dg


   function f1(q,N)
     integer :: N
     real(kind=GRID_SR) ,intent(in) ::q(N,3)
     real(kind=GRID_SR)             ::f1(N,3)
     integer                        ::i

     f1(:,1) = q(:,2)
     do i=1,N
        if(q(i,1).eq.0) then
           f1(i,1) = 0
           f1(i,2) = 0
           f1(i,3) = 0
        else
           f1(i,2) = q(i,2)**2/q(i,1) + 0.5_GRID_SR * g * q(i,1)**2
           f1(i,3) = q(i,2)*q(i,3)/q(i,1)
        end if
     end do

   end function f1

   function f2(q,N)
     integer :: N
     real(kind=GRID_SR) ,intent(in) ::q(N,3)
     real(kind=GRID_SR)             ::f2(N,3)
     integer                        ::i

     f2(:,1) = q(:,3)
     do i=1,N
        if(q(i,1).eq.0) then
           f2(i,1) = 0
           f2(i,2) = 0
           f2(i,3) = 0
        else
           f2(i,2) = q(i,2)*q(i,3)/q(i,1)
           f2(i,3) = q(i,3)**2/q(i,1) + 0.5_GRID_SR * g * q(i,1)**2
        end if
     end do

   end function f2

   function S(q,b_x,b_y,N)
     integer ::N,i
     real(kind=GRID_SR) ,intent(in) ::q(N,3), b_x(N), b_y(N)
     real(kind=GRID_SR)             ::S(N,3)

     S(:,1) = 0
     do i=1,N
        S(i,2)=-q(i,1)*(b_x(i))*g
        S(i,3)=-q(i,1)*(b_y(i))*g 

     end do

   end function S

   function cell_to_edge_op_dg(element, edge) result(rep)
     type(t_element_base), intent(in)						:: element
     type(t_edge_data), intent(in)						    :: edge
     type(num_cell_rep)										:: rep
     integer                                             :: edge_type !1=left, 2=hypotenuse, 3=right
     integer :: i,k,indx,indx_mir
     type(num_cell_rep) ::rep_fv

#if defined(_SWE_DG)		

     associate(cell_edge => ref_plotter_data(- abs(element%cell%geometry%i_plotter_type))%edges, normal => edge%transform_data%normal)
       do i=1,3
          if ((normal(1) == cell_edge(i)%normal(1) .and. normal(2) == cell_edge(i)%normal(2))   &
               .or. (normal(1) == -cell_edge(i)%normal(1) .and. normal(2) == -cell_edge(i)%normal(2))) then
             edge_type = i
          end if
       end do
     end associate

     if(element%cell%data_pers%troubled.le.0) then
     do k=0,_SWE_DG_ORDER
        do i=0,_SWE_DG_ORDER
           indx_mir=i+1+k*(_SWE_DG_ORDER+1)
           select case (edge_type)
           case(1) !left
              indx=1+i+k*_SWE_DG_DOFS
              rep%Q_DG_P(indx_mir)%h = element%cell%data_pers%Q_DG_P(indx)%h
              rep%Q_DG_P(indx_mir)%p = element%cell%data_pers%Q_DG_P(indx)%p
              rep%Q_DG_P(indx_mir)%b = element%cell%data_pers%Q_DG_P(indx)%b
           case(2) !mid
              indx=_SWE_DG_DOFS-(i)*(i+1)/2+k*_SWE_DG_DOFS
              rep%Q_DG_P(indx_mir)%h = element%cell%data_pers%Q_DG_P(indx)%h
              rep%Q_DG_P(indx_mir)%p = element%cell%data_pers%Q_DG_P(indx)%p
              rep%Q_DG_P(indx_mir)%b = element%cell%data_pers%Q_DG_P(indx)%b
           case(3) !right
              indx=_SWE_DG_DOFS-(_SWE_DG_ORDER-i+1)*(_SWE_DG_ORDER-i+2)/2 +1+k*_SWE_DG_DOFS
              rep%Q_DG_P(indx_mir)%h = element%cell%data_pers%Q_DG_P(indx)%h
              rep%Q_DG_P(indx_mir)%p = element%cell%data_pers%Q_DG_P(indx)%p
              rep%Q_DG_P(indx_mir)%b = element%cell%data_pers%Q_DG_P(indx)%b
           end select
        end do
     end do
     end if


!     call element%cell%data_pers%convert_dg_to_fv()

     associate(H => element%cell%data_pers%H, HU => element%cell%data_pers%HU, HV => element%cell%data_pers%HV, B => element%cell%data_pers%B)
     select case (edge_type)
     case (1) !cells with id i*i+1 (left leg)
        do i=0, _SWE_PATCH_ORDER - 1
           rep%H(i+1) = H(i*i + 1)
           rep%HU(i+1) = HU(i*i + 1)
           rep%HV(i+1) = HV(i*i + 1)
           rep%B(i+1) = B(i*i + 1)
        end do
     case (2) ! hypotenuse
        do i=1, _SWE_PATCH_ORDER
           rep%H(i) = H((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*i - 1)
           rep%HU(i) = HU((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*i - 1)
           rep%HV(i) = HV((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*i - 1)
           rep%B(i) = B((_SWE_PATCH_ORDER-1)*(_SWE_PATCH_ORDER-1) + 2*i - 1)
        end do
     case (3) !cells with id i*i (right leg)
        do i=1, _SWE_PATCH_ORDER
           rep%H(_SWE_PATCH_ORDER + 1 - i) = H(i*i)
           rep%HU(_SWE_PATCH_ORDER + 1 - i) = HU(i*i)
           rep%HV(_SWE_PATCH_ORDER + 1 - i) = HV(i*i)
           rep%B(_SWE_PATCH_ORDER + 1 - i) = B(i*i)
        end do
     end select

   end associate


   rep%troubled=element%cell%data_pers%troubled_old  

   rep%max_val(1)=maxval(element%cell%data_pers%H)
   rep%max_val(2)=maxval(element%cell%data_pers%HU)
   rep%max_val(3)=maxval(element%cell%data_pers%HV)

   rep%min_val(1)=minval(element%cell%data_pers%H)
   rep%min_val(2)=minval(element%cell%data_pers%HU)
   rep%min_val(3)=minval(element%cell%data_pers%HV)

#endif				

   end function cell_to_edge_op_dg

   subroutine skeleton_array_op_dg(traversal, grid, edges, rep1, rep2, update1, update2)
     type(t_swe_dg_timestep_traversal), intent(in)	            :: traversal
     type(t_grid_section), intent(in)		            	    :: grid
     type(t_edge_data), intent(in)				    :: edges(:)
     type(num_cell_rep), intent(in)				    :: rep1(:), rep2(:)
     type(num_cell_update), intent(out)				    :: update1(:), update2(:)
     integer (kind = GRID_SI)                                       :: i

     do i = 1, size(edges)
        call skeleton_scalar_op_dg(traversal, grid, edges(i), rep1(i), rep2(i), update1(i), update2(i))
     end do

   end subroutine skeleton_array_op_dg

   subroutine skeleton_scalar_op_dg(traversal, grid, edge, rep1, rep2, update1, update2)
     type(t_swe_dg_timestep_traversal), intent(in)		:: traversal
     type(t_grid_section), intent(in)				:: grid
     type(t_edge_data), intent(in)				:: edge
     type(num_cell_rep), intent(in)				:: rep1, rep2
     type(num_cell_update), intent(out)				:: update1, update2

     update1%Q_DG_P(:)%h   =rep2%Q_DG_P(:)%h
     update1%Q_DG_P(:)%p(1)=rep2%Q_DG_P(:)%p(1)
     update1%Q_DG_P(:)%p(2)=rep2%Q_DG_P(:)%p(2)
     update1%Q_DG_P(:)%b   =rep2%Q_DG_P(:)%b
     update1%troubled      =rep2%troubled
     update1%max_val       =rep2%max_val
     update1%min_val       =rep2%min_val

     update2%Q_DG_P(:)%h   =rep1%Q_DG_P(:)%h
     update2%Q_DG_P(:)%p(1)=rep1%Q_DG_P(:)%p(1)
     update2%Q_DG_P(:)%p(2)=rep1%Q_DG_P(:)%p(2)
     update2%Q_DG_P(:)%b   =rep1%Q_DG_P(:)%b
     update2%troubled      =rep1%troubled
     update2%max_val       =rep1%max_val
     update2%min_val       =rep1%min_val

   end subroutine skeleton_scalar_op_dg

   subroutine bnd_skeleton_array_op_dg(traversal, grid, edges, rep, update)
     type(t_swe_dg_timestep_traversal), intent(in)				:: traversal
     type(t_grid_section), intent(in)						:: grid
     type(t_edge_data), intent(in)						:: edges(:)
     type(num_cell_rep), intent(in)						:: rep(:)
     type(num_cell_update), intent(out)						:: update(:)
     integer :: i
     do i = 1, size(edges)
        call bnd_skeleton_scalar_op_dg(traversal, grid, edges(i), rep(i), update(i))
     end do
   end subroutine bnd_skeleton_array_op_dg

   subroutine bnd_skeleton_scalar_op_dg(traversal, grid, edge, rep, update)
     type(t_swe_dg_timestep_traversal), intent(in)  :: traversal
     type(t_grid_section), intent(in)		    :: grid
     type(t_edge_data), intent(in)		    :: edge
     type(num_cell_rep), intent(in)		    :: rep
     type(num_cell_update), intent(out)		    :: update
     real(kind=GRID_SR)                             :: normal_normed(2)
     real(kind=GRID_SR)                             :: length_flux
     integer                                        :: i

     normal_normed=(edge%transform_data%normal)/NORM2(edge%transform_data%normal)

     update%Q_DG_P(:)%h = rep%Q_DG_P(:)%h
     update%Q_DG_P(:)%b = rep%Q_DG_P(:)%b

     update%max_val=rep%max_val
     update%min_val=rep%min_val

     update%troubled=rep%troubled

     do i=1,(_SWE_DG_ORDER+1)**2
        !                          reflecting boundary
        length_flux = dot_product(rep%Q_DG_P(i)%p, normal_normed)

        update%Q_DG_P(i)%p(1) = rep%Q_DG_P(i)%p(1)-2.0_GRID_SR*length_flux*normal_normed(1)
        update%Q_DG_P(i)%p(2) = rep%Q_DG_P(i)%p(2)-2.0_GRID_SR*length_flux*normal_normed(2)

        !simple outflowing boundary
!        update%Q_DG_P(i)%p(1) =  rep%Q_DG_P(i)%p(1)
!        update%Q_DG_P(i)%p(2) =  rep%Q_DG_P(i)%p(2)
        !                            end if

        !zero vel bnd
        !                          update%Q_DG_P(i)%p(1) =  0
        !                          update%Q_DG_P(i)%p(2) =  0
     end do

   end subroutine bnd_skeleton_scalar_op_dg

   !cup_time
   subroutine cell_update_op_dg(traversal, section, element, update1, update2, update3)
     type(t_swe_dg_timestep_traversal), intent(inout)		:: traversal
     type(t_grid_section), intent(inout)			:: section
     type(t_element_base), intent(inout)			:: element
     type(num_cell_update), intent(inout)			:: update1, update2, update3
     type(num_cell_update) :: tmp
     type(t_update),  DIMENSION((_SWE_DG_ORDER+1)**2)		:: flux1, flux2, flux3
     type(t_state), dimension((_SWE_DG_ORDER+1)**2)             :: edge_l, edge_m,edge_r
     real(kind= GRID_SR) :: dt,dx,delta(3),max_neighbour(3),min_neighbour(3),data_max_val(3),data_min_val(3)
     integer :: indx,indx_mir,k,i
     real (kind=GRID_SR) :: max_wave_speed, wave_speed=0
     real (kind = GRID_SR), dimension (_SWE_PATCH_ORDER_SQUARE) :: H_old, HU_old, HV_old, H, HU, HV, B_old
     logical :: drying,troubled,neighbours_troubled
     integer :: i_depth


     if (element%cell%geometry%i_plotter_type > 0) then ! if orientation = forward, reverse updates
        tmp=update1
        update1=update3
        update3=tmp
     end if

     neighbours_troubled=(update1%troubled.ge.1).or.(update2%troubled.ge.1).or.(update3%troubled.ge.1)

     associate(data => element%cell%data_pers)

     if(neighbours_troubled) then
        data%troubled=merge(data%troubled,2,data%troubled.ge.1)
        call element%cell%data_pers%convert_dg_to_fv_bathymetry()
     end if

     if(data%troubled.le.0) then
        data%troubled=0
        do k=0,_SWE_DG_ORDER
           do i=0,_SWE_DG_ORDER
              indx_mir=i+1+k*(_SWE_DG_ORDER+1)

              indx=i+1+k*_SWE_DG_DOFS
              edge_l(indx_mir)%h = element%cell%data_pers%Q_DG_P(indx)%h
              edge_l(indx_mir)%p = element%cell%data_pers%Q_DG_P(indx)%p
              edge_l(indx_mir)%b = element%cell%data_pers%Q_DG_P(indx)%b
              
              indx=_SWE_DG_DOFS-(_SWE_DG_ORDER-i)*(_SWE_DG_ORDER-i+1)/2+k*_SWE_DG_DOFS
              edge_m(indx_mir)%h = element%cell%data_pers%Q_DG_P(indx)%h
              edge_m(indx_mir)%p = element%cell%data_pers%Q_DG_P(indx)%p
              edge_m(indx_mir)%b = element%cell%data_pers%Q_DG_P(indx)%b
              
              indx=_SWE_DG_DOFS-(_SWE_DG_ORDER-i+1)*(_SWE_DG_ORDER-i+2)/2+1+k*_SWE_DG_DOFS
              edge_r(indx_mir)%h = element%cell%data_pers%Q_DG_P(indx)%h
              edge_r(indx_mir)%p = element%cell%data_pers%Q_DG_P(indx)%p
              edge_r(indx_mir)%b = element%cell%data_pers%Q_DG_P(indx)%b
           end do
        end do
        
       call compute_flux_pred(ref_plotter_data(abs(element%cell%geometry%i_plotter_type))%edges(3)%normal,edge_l,update1%Q_DG_P,flux1)
       call compute_flux_pred(ref_plotter_data(abs(element%cell%geometry%i_plotter_type))%edges(2)%normal,edge_m,update2%Q_DG_P,flux2)
       call compute_flux_pred(ref_plotter_data(abs(element%cell%geometry%i_plotter_type))%edges(1)%normal,edge_r,update3%Q_DG_P,flux3)

        !print*,"updatesl"
        !print*,edge_l(:)%H
        !print*,edge_l(:)%p(1)
        !print*,edge_l(:)%p(2)
        !print*,edge_l(:)%b
        !print*,update1%Q_DG_P(:)%H
        !print*,update1%Q_DG_P(:)%p(1)
        !print*,update1%Q_DG_P(:)%p(2)
        !print*,update1%Q_DG_P(:)%b
        !print*,flux1(:)%H
        !print*,flux1(:)%p(1)
        !print*,flux1(:)%p(2)

        

       H_old =data%H
       B_old =data%B
       HU_old=data%HU
       HV_old=data%HV

       data_max_val(1)=maxval(H_old)
       data_max_val(2)=maxval(HU_old)
       data_max_val(3)=maxval(HV_old)
       
       data_min_val(1)=minval(H_old)
       data_min_val(2)=minval(HU_old)
       data_min_val(3)=minval(HV_old)
       
       max_neighbour=max(update1%max_val,update2%max_val,update3%max_val,data_max_val)
       min_neighbour=min(update1%min_val,update2%min_val,update3%min_val,data_min_val)

       delta(1) = max(0.1_GRID_SR,max_neighbour(1)-min_neighbour(1))
       delta(2) = max(0.1_GRID_SR,max_neighbour(2)-min_neighbour(2))
       delta(3) = max(0.1_GRID_SR,max_neighbour(3)-min_neighbour(3))
       
       delta=delta*1.0e-3_GRID_SR
!       !!!!!!print*,"solver"
       call dg_solver(element,flux1,flux2,flux3,section%r_dt)

       call element%cell%data_pers%convert_dg_to_fv_bathymetry()
       call data%convert_dg_to_fv()

       H =data%H
       HU=data%HU
       HV=data%HV


#if defined(_SWE_DG_LIMITER_UNLIMITED)
       troubled=.false.
#elif defined(_SWE_DG_LIMITER_HEIGHT)
       troubled=.not.all(data%H< max_neighbour(1)+delta(1).and.data%H >min_neighbour(1)-delta(1))
#elif defined(_SWE_DG_LIMITER_ALL)
       troubled=&
            .not.all(data%H< max_neighbour(1)+delta(1).and.data%H >min_neighbour(1)-delta(1)) .or.&
            .not.all(data%HU< max_neighbour(2)+delta(2).and.data%HU>min_neighbour(2)-delta(2)) .or.&
            .not.all(data%HV< max_neighbour(3)+delta(3).and.data%HV>min_neighbour(3)-delta(3))
#endif       
       drying=.not.all(data%H - data%B > cfg%dry_tolerance*50.0).or..not.all(data%Q_DG%H > cfg%dry_tolerance*50.0)

       if(troubled.or.drying) then
          data%H=H_old
          data%HU=HU_old
          data%HV=HV_old
          data%B=B_old
          data%troubled=merge(1,3,drying)
       else

          do i=1,_SWE_DG_DOFS
             wave_speed =  sqrt(g * (data%Q_DG(i)%h)) + maxval(abs(data%Q_DG(i)%p/data%Q_DG(i)%h))
             section%r_dt_new = min(section%r_dt_new,cfg%scaling*  element%transform_data%custom_data%scaling  / (wave_speed* (_SWE_DG_ORDER*4.0_GRID_SR +2.0_GRID_SR)))
             max_wave_speed = max(max_wave_speed,wave_speed)
          end do
          i_depth = element%cell%geometry%i_depth

#define  _SWE_DG_REFINEMENT_COAST_HEIGHT 200
#define  _SWE_DG_REFINEMENT_COAST_DEPTH 12
#define  _SWE_DG_REFINEMENT_WAVE_DEPTH 10
#define  _SWE_DG_REFINEMENT_LAKE_AT_REST_DEPTH 4

          if (minval(H-data%B).le._SWE_DG_REFINEMENT_COAST_HEIGHT*cfg%dry_tolerance .and. i_depth < min(_SWE_DG_REFINEMENT_COAST_DEPTH,cfg%i_max_depth)) then
             element%cell%geometry%refinement = 1
             traversal%i_refinements_issued = traversal%i_refinements_issued + 1_GRID_DI
          else if(max_wave_speed > 1e-15 .and. i_depth < min(_SWE_DG_REFINEMENT_WAVE_DEPTH,cfg%i_max_depth)) then
             element%cell%geometry%refinement = 1
             traversal%i_refinements_issued = traversal%i_refinements_issued + 1_GRID_DI
          else if(max_wave_speed < 1e-15 .and. i_depth > max(_SWE_DG_REFINEMENT_LAKE_AT_REST_DEPTH,cfg%i_min_depth)) then
             element%cell%geometry%refinement = -1
          end if
        end if
    end if

    !if cell is troubled mark all edges as troubled
    if(data%troubled.ge.1) then
       do i=1,size(element%edges)
          element%edges(i)%ptr%data_pers%troubled= .true.
       end do
    end if
  end associate

end subroutine cell_update_op_dg


  subroutine compute_flux_pred(normal,QL, QR,flux)
    
    type(t_state),Dimension((_SWE_DG_ORDER+1)**2), intent(inout)	:: QL
    type(t_state),Dimension((_SWE_DG_ORDER+1)**2), intent(inout)   :: QR
    type(t_update),Dimension((_SWE_DG_ORDER+1)**2), intent(out) :: flux
    type(t_update),Dimension((_SWE_DG_ORDER+1)**2) 		:: flux_r
#if !defined(_SWE_DG_NODAL)
    type(t_update),Dimension(size(bnd_gl_node_vals,1)):: flux_l2
    type(t_update),Dimension(size(bnd_gl_node_vals,1)):: flux_r_l2
    type(t_state),Dimension(size(bnd_gl_node_vals,1)) :: QL_l2
    type(t_state),Dimension(size(bnd_gl_node_vals,1)) :: QR_l2
#endif
    
    real(kind=GRID_SR),Dimension(1,3) :: f1r_s,f1l_s,f2r_s,f2l_s,f,ql_v,qr_v
    real(kind = GRID_SR), intent(in)  :: normal(2)
    real(kind = GRID_SR)	      :: normal_normed(2)
    real(kind = GRID_SR)	      :: epsilon
    
    integer ::i
#if defined(_LLF_FLUX_DG)
    real(kind=GRID_SR)								:: vL, vR, alpha
    real(kind=GRID_SR) :: max_wave_speedR
#endif
    
    normal_normed=normal/NORM2(normal)
    
    flux%h=0
    flux%p(1)=0
    flux%p(2)=0
    flux_r%h=0
    flux_r%p(1)=0
    flux_r%p(2)=0
    
#if !defined(_SWE_DG_NODAL)
    QL_l2(:)%h=matmul(bnd_gl_node_vals,QL%h)
    QL_l2(:)%p(1)=matmul(bnd_gl_node_vals,QL%p(1))
    QL_l2(:)%p(2)=matmul(bnd_gl_node_vals,QL%p(2))
    QL_l2(:)%b=matmul(bnd_gl_node_vals,QL%b)
    
    QR_l2(:)%h=matmul(bnd_gl_node_vals,QR%h)
    QR_l2(:)%p(1)=matmul(bnd_gl_node_vals,QR%p(1))
    QR_l2(:)%p(2)=matmul(bnd_gl_node_vals,QR%p(2))
    QR_l2(:)%b=matmul(bnd_gl_node_vals,QR%b)
#endif
    
#if defined(_LLF_FLUX_DG)
#if defined(_SWE_DG_NODAL)
    alpha=0
    
    do i=1,_SWE_DG_ORDER+1
       
       if(.not.(QL(i)%h < cfg%dry_tolerance) ) then 
          vL = DOT_PRODUCT(normal_normed, QL(i)%p/QL(i)%h)
          flux(i)%max_wave_speed = sqrt(g * (QL(i)%h)) + abs(vL)
       else
          flux(i)%max_wave_speed = 0
       end if
       
       if(.not.(QR(i)%h < cfg%dry_tolerance) ) then 
          vR = DOT_PRODUCT(normal_normed, QR(i)%p/QR(i)%h)
          max_wave_speedR = sqrt(g * (QR(i)%h)) + abs(vR)
       else
          max_wave_speedR = 0
       end if
       
       alpha = max(alpha,flux(i)%max_wave_speed, max_wave_speedR)

    end do

    do i=1,(_SWE_DG_ORDER+1)**2
       ql_v(1,1)= QL(i)%h
       ql_v(1,2:3)= QL(i)%p

       qr_v(1,1)= QR(i)%h
       qr_v(1,2:3)= QR(i)%p

       f1r_s=f1(qr_v,1)
       f2r_s=f2(qr_v,1)
       f1l_s=f1(ql_v,1)
       f2l_s=f2(ql_v,1)
       f=0.5_GRID_SR*((f1r_s-f1l_s) * normal_normed(1) + (f2r_s-f2l_s) *normal_normed(2)+alpha*(ql_v-qr_v))
       flux(i)%h=f(1,1)
       flux(i)%p=f(1,2:3)

       ! end if
    end do
#else
    alpha=0

    do i=1,size(bnd_gl_node_vals,1)

       if(.not.(QL_L2(i)%h < cfg%dry_tolerance) ) then 
          vL = DOT_PRODUCT(normal_normed, QL_L2(i)%p/QL_L2(i)%h)
          flux_L2(i)%max_wave_speed = sqrt(g * (QL_L2(i)%h)) + abs(vL)
       else
          flux_L2(i)%max_wave_speed = 0
       end if

       if(.not.(QR_L2(i)%h < cfg%dry_tolerance) ) then 
          vR = DOT_PRODUCT(normal_normed, QR_L2(i)%p/QR_L2(i)%h)
          max_wave_speedR = sqrt(g * (QR_L2(i)%h)) + abs(vR)
       else
          max_wave_speedR = 0
       end if

       alpha = max(alpha,flux_L2(i)%max_wave_speed, max_wave_speedR)

    end do

    do i=1,size(bnd_gl_node_vals,1)
       ql_v(1,1)= QL_L2(i)%h
       ql_v(1,2:3)= QL_L2(i)%p

       qr_v(1,1)= QR_L2(i)%h
       qr_v(1,2:3)= QR_L2(i)%p

       f1r_s=f1(qr_v,1)
       f2r_s=f2(qr_v,1)
       f1l_s=f1(ql_v,1)
       f2l_s=f2(ql_v,1)

       f=0.5_GRID_SR*((f1r_s-f1l_s) * normal_normed(1) + (f2r_s-f2l_s) *normal_normed(2)+alpha*(ql_v-qr_v))
       flux_l2(i)%h=f(1,1)
       flux_l2(i)%p=f(1,2:3)

       ! end if
    end do
#endif
#elif                 (defined(_FWAVE_FLUX)|defined(_AUG_RIEMANN_FLUX)|defined(_HLLE_FLUX))

#if defined(_SWE_DG_NODAL)

    flux%h=0
    flux%p(1)=0
    flux%p(2)=0
    flux_r%h=0
    flux_r%p(1)=0
    flux_r%p(2)=0

    do i=1,(_SWE_DG_ORDER+1)**2
    !    epsilon=1.0e-14
    !    if(abs(QL(i)%p(1))<epsilon)then
    !       QL(i)%p(1)=0.0_GRID_SR
    !    end if
    !    if(abs(QR(i)%p(1))<epsilon)then
    !       QR(i)%p(1)=0.0_GRID_SR
    !    end if
    !    if(abs(QL(i)%p(2))<epsilon)then
    !       QL(i)%p(2)=0.0_GRID_SR
    !    end if
    !    if(abs(QR(i)%p(2))<epsilon)then
    !       QR(i)%p(2)=0.0_GRID_SR
    !    end if
    !    if(QL(i)%b.ne.0) then
    !       if(abs(QL(i)%b-QR(i)%b)/QL(i)%b < epsilon)then
    !          QL(i)%b=QR(i)%b
    !       end if
    !    else if(abs((QL(i)%b-QR(i)%b)) < epsilon)then
    !          QL(i)%b=QR(i)%b
    !    end if

    !    if(QL(i)%h.ne.0) then
    !       if(abs(QL(i)%h-QR(i)%h)/QL(i)%b < epsilon)then
    !          QL(i)%h=QR(i)%h
    !       end if
    !    else if(abs(QL(i)%h-QR(i)%h) < epsilon)then
    !       QL(i)%h=QR(i)%h
    !    end if



       call compute_geoclaw_flux_pred(normal_normed, QL(i), QR(i), flux(i), flux_r(i))
    end do

#else
    do i=1,size(bnd_gl_node_vals,1)
       call compute_geoclaw_flux_pred(normal_normed, QL_l2(i), QR_l2(i), flux_l2(i), flux_r_l2(i))
    end do
#endif

#endif

#if !defined(_SWE_DG_NODAL)

    flux%h=matmul(bnd_gl_weights,flux_l2%h)
    flux%p(1)=matmul(bnd_gl_weights,flux_l2%p(1))
    flux%p(2)=matmul(bnd_gl_weights,flux_l2%p(2))
#endif
  end subroutine compute_flux_pred

  subroutine dg_solver(element,update1,update2,update3,dt)
    type(t_element_base), intent(in)				:: element
    type(t_update), dimension((_SWE_DG_ORDER+1)**2),intent(in)  :: update1, update2, update3
    real(kind=GRID_SR)                                          :: q(_SWE_DG_DOFS,3),&
         q_p((_SWE_DG_ORDER+1)*_SWE_DG_DOFS,3),q_temp(_SWE_DG_DOFS,3),&
         nF1(_SWE_DG_DOFS*(_SWE_DG_ORDER+1),3),&
         nF2(_SWE_DG_DOFS*(_SWE_DG_ORDER+1),3),&
         nF3(_SWE_DG_DOFS*(_SWE_DG_ORDER+1),3)
    real(kind=GRID_SR)                                          :: dt
    integer :: i,j,k,inx_l,inx_m,inx_r
    real(kind=GRID_SR) :: dx
    real(kind=GRID_SR),Dimension(_SWE_DG_DOFS*(_SWE_DG_ORDER+1),3) :: f1_hat,f1_s,f2_hat,f2_s,S_s
    real(kind=GRID_SR),Dimension(_SWE_DG_DOFS*(_SWE_DG_ORDER+1),3) :: f1_hat_speed,f1_s_speed,f2_hat_speed,f2_s_speed
    real(kind=GRID_SR),Dimension(_SWE_DG_DOFS*(_SWE_DG_ORDER+1),3) :: f1_mom_der_hat,f1_mom_der_s,f2_mom_der_hat,f2_mom_der_s
    real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,3) :: volume_flux1,volume_flux2,bnd_flux_contrib

    real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,3) :: source,bnd_source_contrib,source_temp,bnd_source_contrib_temp,orig_source
    real(kind=GRID_SR),Dimension(_SWE_DG_DOFS,3) :: bnd_flux_l,bnd_flux_m,bnd_flux_r
#if defined(_SWE_DG_NODAL)
    real(kind=GRID_SR),Dimension(_SWE_DG_DOFS*(_SWE_DG_ORDER+1)) :: H_x,H_y,H_x_temp,H_y_temp,b_x,b_y
#else
    real(kind=GRID_SR),Dimension(size(st_gl_node_vals,1)) :: b_x,b_y
    real(kind=GRID_SR),Dimension(size(st_gl_node_vals,1),3) :: q_v
    integer :: st_nodes = size(st_gl_node_vals,1)
    integer :: s_nodes = size(s_der_x_gl_node_vals,1)
#endif

    real(kind=GRID_SR),Dimension(2,2) :: jacobian,jacobian_inv

    !!!!!!!!!print*,"solve"

    jacobian=ref_plotter_data(abs(element%cell%geometry%i_plotter_type))%jacobian
    jacobian= jacobian/NORM2(jacobian(:,1))

    jacobian_inv=ref_plotter_data(abs(element%cell%geometry%i_plotter_type))%jacobian_inv
    jacobian_inv= jacobian_inv/NORM2(jacobian_inv(:,1))

    associate(Q_DG => element%cell%data_pers,Q_DG_P => element%cell%data_pers%Q_DG_P)

    call Q_DG%dofs_to_vec_dg(q)
    call Q_DG%dofs_to_vec_dg_p(q_p)

    dx=element%transform_data%custom_data%scaling*cfg%scaling
    !                        dx=element%transform_data%custom_data%scaling

#if defined(_SWE_DG_NODAL)
    nF1=0
    nF2=0
    nF3=0
    do k=0,_SWE_DG_ORDER
       do i=0,_SWE_DG_ORDER
          j=i+1+k*(_SWE_DG_ORDER+1)
          inx_l=i+1+k*_SWE_DG_DOFS
          inx_m=_SWE_DG_DOFS-(_SWE_DG_ORDER-i)*(_SWE_DG_ORDER-i+1)/2+k*_SWE_DG_DOFS
          inx_r=_SWE_DG_DOFS-(_SWE_DG_ORDER-i+1)*(_SWE_DG_ORDER-i+2)/2+1+k*_SWE_DG_DOFS

          nF1(inx_l,1)=update1(j)%h
          nF1(inx_l,2)=update1(j)%p(1)
          nF1(inx_l,3)=update1(j)%p(2)

          NF2(inx_m,1)=update2(j)%h
          NF2(inx_m,2)=update2(j)%p(1)
          NF2(inx_m,3)=update2(j)%p(2)

          NF3(inx_r,1)=update3(j)%h
          NF3(inx_r,2)=update3(j)%p(1)
          NF3(inx_r,3)=update3(j)%p(2)
       end do
    end do

    bnd_flux_l=matmul(b_m_1,nF1)
    bnd_flux_m=matmul(b_m_2,nF2)
    bnd_flux_r=matmul(b_m_3,nF3)
#else
    bnd_flux_l=0
    bnd_flux_m=0
    bnd_flux_r=0

    do j=0,_SWE_DG_ORDER
       do i=0,_SWE_DG_ORDER
          inx_l=i+1
          inx_m=_SWE_DG_DOFS-(_SWE_DG_ORDER-i)*(_SWE_DG_ORDER-i+1)/2
          inx_r=_SWE_DG_DOFS-(_SWE_DG_ORDER-i+1)*(_SWE_DG_ORDER-i+2)/2+1

          bnd_flux_l(inx_l,1)=bnd_flux_l(inx_l,1)+update1(1+j*(_SWE_DG_ORDER+1)+i)%h
          bnd_flux_m(inx_m,1)=bnd_flux_m(inx_m,1)+update2(1+j*(_SWE_DG_ORDER+1)+i)%h * sqrt(2.0_GRID_SR)
          bnd_flux_r(inx_r,1)=bnd_flux_r(inx_r,1)+update3(1+j*(_SWE_DG_ORDER+1)+i)%h

          bnd_flux_l(inx_l,2)=bnd_flux_l(inx_l,2)+update1(1+j*(_SWE_DG_ORDER+1)+i)%p(1)
          bnd_flux_m(inx_m,2)=bnd_flux_m(inx_m,2)+update2(1+j*(_SWE_DG_ORDER+1)+i)%p(1) * sqrt(2.0_GRID_SR)
          bnd_flux_r(inx_r,2)=bnd_flux_r(inx_r,2)+update3(1+j*(_SWE_DG_ORDER+1)+1)%p(1)

          bnd_flux_l(inx_l,3)=bnd_flux_l(inx_l,3)+update1(1+j*(_SWE_DG_ORDER+1)+1)%p(2)
          bnd_flux_m(inx_m,3)=bnd_flux_m(inx_m,3)+update2(1+j*(_SWE_DG_ORDER+1)+1)%p(2) * sqrt(2.0_GRID_SR)

          bnd_flux_r(inx_r,3)=bnd_flux_r(inx_r,3)+update3(1+j*(_SWE_DG_ORDER+1)+1)%p(2)
       end do
    end do

#endif

#if defined(_SWE_DG_NODAL)

    H_x=matmul(basis_der_st_x,Q_DG_P(:)%h+Q_DG_P(:)%B)
    H_y=matmul(basis_der_st_y,Q_DG_P(:)%h+Q_DG_P(:)%B)


    ! where(q_p(:,1).ne.0 .and. abs(H_x/q_p(:,1))< 1.0e-14)
    !    H_x = 0.0_GRID_SR
    ! end where

    ! where(q_p(:,1).ne.0 .and. abs(H_y/q_p(:,1))< 1.0e-14)
    !    H_y = 0.0_GRID_SR
    ! end where


    source_temp=0
    source_temp(:,2)=-matmul(s_m,(g*Q_DG_P%H*H_x))
    source_temp(:,3)=-matmul(s_m,(g*Q_DG_P%H*H_y))
    
    source_temp(:,2)=source_temp(:,2)-matmul(s_k_x,(0.5_GRID_SR*g*Q_DG_P%H**2))
    source_temp(:,3)=source_temp(:,3)-matmul(s_k_y,(0.5_GRID_SR*g*Q_DG_P%H**2))

    source=0
    source(:,2)=jacobian(1,1) * source_temp(:,2) + jacobian(1,2) * source_temp(:,3)
    source(:,3)=jacobian(2,1) * source_temp(:,2) + jacobian(2,2) * source_temp(:,3)

    bnd_source_contrib_temp=0
    bnd_source_contrib_temp(:,2)=-matmul(b_m_3,(0.5_GRID_SR*g*Q_DG_P%H**2))+matmul(b_m_2,0.5_GRID_SR*g*Q_DG_P%H**2)/sqrt(2.0_GRID_SR)
    bnd_source_contrib_temp(:,3)=-matmul(b_m_1,(0.5_GRID_SR*g*Q_DG_P%H**2))+matmul(b_m_2,0.5_GRID_SR*g*Q_DG_P%H**2)/sqrt(2.0_GRID_SR)

    bnd_source_contrib=0
    bnd_source_contrib(:,2)= bnd_source_contrib_temp(:,2)*jacobian(1,1)  + bnd_source_contrib_temp(:,3)*jacobian(1,2)
    bnd_source_contrib(:,3)= bnd_source_contrib_temp(:,2)*jacobian(2,1)  + bnd_source_contrib_temp(:,3)*jacobian(2,2)  


    b_x=Q_DG_P%b_x
    b_y=Q_DG_P%b_y

    S_s =S(q_p,b_x,b_y,_SWE_DG_DOFS*(_SWE_DG_ORDER+1))

    f1_hat=f1(q_p,_SWE_DG_DOFS*(_SWE_DG_ORDER+1))
    f2_hat=f2(q_p,_SWE_DG_DOFS*(_SWE_DG_ORDER+1))

    orig_source=matmul(s_m,S_s)
#else

    q_v=matmul(st_gl_node_vals,q_p)

    f1_hat=matmul(st_gl_weights,f1(q_v,st_nodes))
    f2_hat=matmul(st_gl_weights,f2(q_v,st_nodes))

    b_x=element%cell%data_pers%b_x
    b_y=element%cell%data_pers%b_y

    S_s =matmul(st_gl_weights,S(q_v(:,1),b_x,b_y,st_nodes))

    source=0

    do i=0,_SWE_DG_ORDER
       source=source+S_s(1+_SWE_DG_DOFS*i:_SWE_DG_DOFS*(1+i),:)
    end do

    call lusolve(st_m_lu,_SWE_DG_DOFS*(_SWE_DG_ORDER+1),st_m_lu_pivot,f1_hat(:,1))
    call lusolve(st_m_lu,_SWE_DG_DOFS*(_SWE_DG_ORDER+1),st_m_lu_pivot,f2_hat(:,1))

    call lusolve(st_m_lu,_SWE_DG_DOFS*(_SWE_DG_ORDER+1),st_m_lu_pivot,f1_hat(:,2))
    call lusolve(st_m_lu,_SWE_DG_DOFS*(_SWE_DG_ORDER+1),st_m_lu_pivot,f2_hat(:,2))

    call lusolve(st_m_lu,_SWE_DG_DOFS*(_SWE_DG_ORDER+1),st_m_lu_pivot,f1_hat(:,3))
    call lusolve(st_m_lu,_SWE_DG_DOFS*(_SWE_DG_ORDER+1),st_m_lu_pivot,f2_hat(:,3))

#endif                        

    f1_s=(jacobian_inv(1,1)*f1_hat+jacobian_inv(1,2)*f2_hat)
    f2_s=(jacobian_inv(2,1)*f1_hat+jacobian_inv(2,2)*f2_hat)

    volume_flux1=matmul(s_k_x,f1_s)
    volume_flux2=matmul(s_k_y,f2_s)

    bnd_flux_contrib=-matmul(b_m_1,f2_s)-matmul(b_m_3,f1_s)+matmul(b_m_2,f1_s+f2_s)/sqrt(2.0_GRID_SR)

    !print*,"type"
    !print*,abs(element%cell%geometry%i_plotter_type)
    !print*,"volflx"
    !print*, volume_flux1(:,1)+volume_flux2(:,1)
    !print*, volume_flux1(:,2)+volume_flux2(:,2)
    !print*, volume_flux1(:,3)+volume_flux2(:,3)
    !print*,"src"
    !print*,source(:,1)
    !print*,source(:,2)
    !print*,source(:,3)
    !print*,"bnd_flux"
    !print*,bnd_flux_contrib(:,1)
    !print*,bnd_flux_contrib(:,2)
    !print*,bnd_flux_contrib(:,3)
    !print*,"bnd_src"
    !print*,bnd_source_contrib(:,1)
    !print*,bnd_source_contrib(:,2)
    !print*,bnd_source_contrib(:,3)

    !print*,"src_comp"
    !print*,source(:,1)-bnd_source_contrib(:,1)
    !print*,source(:,2)-bnd_source_contrib(:,2)
    !print*,source(:,3)-bnd_source_contrib(:,3)

    !print*,"orig_source"
    !print*,orig_source(:,1)
    !print*,orig_source(:,2)
    !print*,orig_source(:,3)
    !print*
    !print*,"sum"
    !print*, -volume_flux1(:,1)-volume_flux2(:,1)-source(:,1) + bnd_flux_contrib(:,1) -bnd_source_contrib(:,1)
    !print*, -volume_flux1(:,2)-volume_flux2(:,2)-source(:,2) + bnd_flux_contrib(:,2) -bnd_source_contrib(:,2)
    !print*, -volume_flux1(:,3)-volume_flux2(:,3)-source(:,3) + bnd_flux_contrib(:,3) -bnd_source_contrib(:,3)
    !print*
    !print*,bnd_flux_l
    !print*,bnd_flux_m
    !print*,bnd_flux_r

    !print*,"Q_DG_P"
    !print*,q_p(:,1)
    !print*,q_p(:,2)
    !print*,q_p(:,3)

    q_temp = ( bnd_flux_l &
         + bnd_flux_m &
         + bnd_flux_r &
         - volume_flux1 &
         - volume_flux2 &
         - source & 
         - bnd_source_contrib &
         + bnd_flux_contrib &
         )* dt/dx  
    !!!!!!!!!print*,"solver"
    !!!!!!!!!print*, "orig ",q
    !!!!!!!!!print*,"orientation",abs(element%cell%geometry%i_plotter_type)
    !!!!!!!!!print*,"jacobian_inv",jacobian_inv
    !print*, "dx ",dx
    !print*, "dt ",dt
    ! !!!!!!print*
    ! !!!!!!print*, "elt flux 1 "
    ! !!!!!!print*,volume_flux1(:,1)
    ! !!!!!!print*,volume_flux1(:,2)
    ! !!!!!!print*,volume_flux1(:,3)
    ! !!!!!!print*
    ! !!!!!!print*, "elt flux 2 "
    ! !!!!!!print*,volume_flux2(:,1)
    ! !!!!!!print*,volume_flux2(:,2)
    ! !!!!!!print*,volume_flux2(:,3)
    ! !!!!!!print*
    ! !!!!!!!!!print*, "elt flux sum"
    ! !!!!!!!!!print*,volume_flux2(:,1)+volume_flux1(:,1)
    ! !!!!!!!!!print*,volume_flux2(:,2)+volume_flux1(:,2)
    ! !!!!!!!!!print*,volume_flux2(:,3)+volume_flux1(:,3)
    ! !!!!!!!!!print*
    ! !!!!!!!!!print*, "bnd_flux_contribution "
    ! !!!!!!!!!print*,bnd_flux_contrib(:,1)
    ! !!!!!!!!!print*,bnd_flux_contrib(:,2)
    ! !!!!!!!!!print*,bnd_flux_contrib(:,3)
    ! !!!!!!!!!print*
    ! !!!!!!!!!print*, "bnd_flux "
    ! !!!!!!!!!print*,bnd_flux_l(:,1)+bnd_flux_r(:,1)+bnd_flux_m(:,1)
    ! !!!!!!!!!print*,bnd_flux_l(:,2)+bnd_flux_r(:,2)+bnd_flux_m(:,2)
    ! !!!!!!!!!print*,bnd_flux_l(:,3)+bnd_flux_r(:,3)+bnd_flux_m(:,3)
    ! !!!!!!!!!print*,                        
    !!!!!!print*, "volume_flux comp"
    !!!!!!print*, volume_flux1(:,1)+volume_flux2(:,1)-bnd_flux_contrib(:,1)
    !!!!!!print*, volume_flux1(:,2)+volume_flux2(:,2)-bnd_flux_contrib(:,2)
    !!!!!!print*, volume_flux1(:,3)+volume_flux2(:,3)-bnd_flux_contrib(:,3)
    !!!!!!print*
    !!!!!!print*, "source"
    !!!!!!print*,source(:,1)
    !!!!!!print*,source(:,2)
    !!!!!!print*,source(:,3)
    !!!!!!print*
    !!!!!!print*, "comp"
    !!!!!!print*,(volume_flux1(:,1)+volume_flux2(:,1)-bnd_flux_contrib(:,1)+source(:,1))*dt/dx
    !!!!!!print*,(volume_flux1(:,2)+volume_flux2(:,2)-bnd_flux_contrib(:,2)+source(:,2))*dt/dx
    !!!!!!print*,(volume_flux1(:,3)+volume_flux2(:,3)-bnd_flux_contrib(:,3)+source(:,3))*dt/dx
    !!!!!!print*
    ! !!!!!!!!!print*, "bnd_flux_1_dof h"
    ! !!!!!!!!!print*,update1%h
    ! !!!!!!!!!print*,update1%p(1)
    ! !!!!!!!!!print*,update1%p(2)
    ! !!!!!!!!!print*
    ! !!!!!!!!!print*, "bnd_flux_2_dof"
    ! !!!!!!!!!print*,update2%h
    ! !!!!!!!!!print*,update2%p(1)
    ! !!!!!!!!!print*,update2%p(2)
    ! !!!!!!!!!print*
    ! !!!!!!!!!print*, "bnd_flux_3_dof"
    ! !!!!!!!!!print*,update3%h
    ! !!!!!!!!!print*,update3%p(1)
    ! !!!!!!!!!print*,update3%p(2)
    ! !!!!!!!!!print*
    ! !!!!!!!!!print*, "volume_flux1 ",volume_flux1
    ! !!!!!!!!!print*,                           
    ! !!!!!!!!!print*, "volume_flux2 ",volume_flux2
    ! !!!!!!!!!print*,
     ! !!!!!!print*, "f1_s "
     ! !!!!!!print*,f1_s(:,1)
     ! !!!!!!print*,f1_s(:,2)
     ! !!!!!!print*,f1_s(:,3)
     ! !!!!!!print*,
     ! !!!!!!print*, "f2_s 2 "
     ! !!!!!!print*,f2_s(:,1)
     ! !!!!!!print*,f2_s(:,2)
     ! !!!!!!print*,f2_s(:,3)
     ! !!!!!!print*,
    ! !!!!!!!!!print*, "f1_hat "
    ! !!!!!!!!!print*,f1_hat(:,1)
    ! !!!!!!!!!print*,f1_hat(:,2)
    ! !!!!!!!!!print*,f1_hat(:,3)
    ! !!!!!!!!!print*,
    ! !!!!!!!!!print*, "f2_hat "
    ! !!!!!!!!!print*,f2_hat(:,1)
    ! !!p!!!!!!!print*,f2_hat(:,2)
    ! !!!!!!!!!print*,f2_hat(:,3)
     ! !!!!!!Print*, "elt_source"
     ! !!!!!!print*,S_s(:,1)
     ! !!!!!!print*,S_s(:,2)
     ! !!!!!!print*,S_s(:,3)
    !  !!!!!!print*
    !  !!!!!!print*,"Bathy x"
    !  !!!!!!print*,b_x
    !  !!!!!!print*,"Bathy y"
    !  !!!!!!print*,b_y
    ! !!!!!!!!!print*,"orig"
    ! !!!!!!!!!print*,q
    ! !!!!!!!!!print*
    ! !!!!!!print*,"Q_DG_P"
    ! !!!!!!print*,q_p(:,1)
    ! !!!!!!print*,q_p(:,2)
    ! !!!!!!print*,q_p(:,3)
    ! !!!!!!print*,Q_DG_P%B
    ! !!!!!!print*
    !!!!!!!!!print*,"pred height"
    !!!!!!!!!print*,element%cell%data_pers%Q_DG_P(:)%h
    !!!!!!!!!print*
    !!!!!!!!!print*,"pred bathy"
    !!!!!!!!!print*,element%cell%data_pers%Q_DG_P(:)%b
    !!!!!!!!!print*
    !!!!!!!!!print*,"q_temp",q_temp
    !!!!!!!!!print*

    
    !print*,"dg temp"
    !print*,q_temp(:,1)
    !print*,q_temp(:,2)
    !print*,q_temp(:,3)
    !print*


    call lusolve(s_m_lu, _SWE_DG_DOFS, s_m_lu_pivot,q_temp(:,1))
    call lusolve(s_m_lu, _SWE_DG_DOFS, s_m_lu_pivot,q_temp(:,2))
    call lusolve(s_m_lu, _SWE_DG_DOFS, s_m_lu_pivot,q_temp(:,3))


    !print*,"dg temp lu"
    !print*,q_temp(:,1)
    !print*,q_temp(:,2)
    !print*,q_temp(:,3)
    !print*

    !print*,"dg orig"
    !print*,q(:,1)
    !print*,q(:,2)
    !print*,q(:,3)
    !print*


    q=q-q_temp

    !print*,"dg res"
    !print*,q(:,1)
    !print*,q(:,2)
    !print*,q(:,3)
    !print*
    call Q_DG%vec_to_dofs_dg(q)

  end associate

end subroutine dg_solver

subroutine compute_geoclaw_flux_pred(normal, QL, QR, fluxL, fluxR)
  type(t_state), intent(in)           :: QL, QR
  type(t_update), intent(out)         :: fluxL, fluxR
  real(kind = GRID_SR), intent(in)    :: normal(2)
  real(kind = GRID_SR)		      :: transform_matrix(2, 2)
  real(kind = GRID_SR)		      :: net_updatesL(3), net_updatesR(3), max_wave_speed
  real(kind = GRID_SR)                :: pL(2), pR(2), hL, hR, bL, bR
  
  transform_matrix(1, :) = normal
  transform_matrix(2, :) = [-normal(2), normal(1)]
  
  pL = matmul(transform_matrix, QL%p)
  pR = matmul(transform_matrix, QR%p)
  
  hL = QL%h
  hR = QR%h
  
  bL = QL%b
  bR = QR%b
  
#           if defined(_FWAVE_FLUX)
  call c_bind_geoclaw_solver(GEOCLAW_FWAVE, 1, 3, hL, hR, pL(1), pR(1), pL(2), pR(2), bL, bR, real(cfg%dry_tolerance, GRID_SR), g, net_updatesL, net_updatesR, max_wave_speed)
#           elif defined(_AUG_RIEMANN_FLUX)
  call c_bind_geoclaw_solver(GEOCLAW_AUG_RIEMANN, 1, 3, hL, hR, pL(1), pR(1), pL(2), pR(2), bL, bR, real(cfg%dry_tolerance, GRID_SR), g, net_updatesL, net_updatesR, max_wave_speed)
#           elif defined(_HLLE_FLUX)
  call compute_updates_hlle_single(hL, hR, pL(1), pR(1), pL(2), pR(2), bL, bR, net_updatesL, net_updatesR, max_wave_speed)
#           endif
  
  fluxL%h = net_updatesL(1)
  fluxL%p = matmul(net_updatesL(2:3), transform_matrix)
  fluxL%max_wave_speed = max_wave_speed
  
  fluxR%h = net_updatesR(1)
  fluxR%p = matmul(net_updatesR(2:3), transform_matrix)
  fluxR%max_wave_speed = max_wave_speed
end subroutine compute_geoclaw_flux_pred

END MODULE SWE_DG_timestep

#endif
#endif
