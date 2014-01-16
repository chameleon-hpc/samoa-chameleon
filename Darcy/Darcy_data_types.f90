! Sam(oa)² - SFCs and Adaptive Meshes for Oceanic And Other Applications
! Copyright (C) 2010 Oliver Meister, Kaveh Rahnema
! This program is licensed under the GPL, for details see the file LICENSE


#include "Compilation_control.f90"

#if defined(_DARCY)
	MODULE Darcy_data_types
		PUBLIC

		!data precision

		integer, PARAMETER :: GRID_SR = selected_real_kind(14,40)
		integer, PARAMETER :: GRID_DR = selected_real_kind(28,80)

		integer, PARAMETER :: GRID_SI = selected_int_kind(8)
		integer, PARAMETER :: GRID_DI = selected_int_kind(16)

		integer, PARAMETER :: GRID_SL = 1

		!*********************************************
		!Persistent Entity data (geometric association)
		!**********************************************

		!> persistent, scenario specific data on a node
		type num_node_data_pers
			real (kind = GRID_SR)   :: p(_DARCY_P_NODE_SIZE)
			real (kind = GRID_SR)   :: A_d(_DARCY_P_NODE_SIZE), d(_DARCY_P_NODE_SIZE), r(_DARCY_P_NODE_SIZE)

			real (kind = GRID_SR)   :: u(2, _DARCY_U_NODE_SIZE)
			real (kind = GRID_SR)   :: saturation(_DARCY_FLOW_NODE_SIZE)    !< water saturation
		END type

		!> persistent, scenario specific data on an edge
		type num_edge_data_pers
			real (kind = GRID_SR)   :: p(_DARCY_P_EDGE_SIZE)
			real (kind = GRID_SR)   :: A_d(_DARCY_P_EDGE_SIZE), d(_DARCY_P_EDGE_SIZE), r(_DARCY_P_EDGE_SIZE)

			real (kind = GRID_SR)   :: u(2, _DARCY_U_EDGE_SIZE)		        !< velocity
			real (kind = GRID_SR)   :: saturation(_DARCY_FLOW_EDGE_SIZE)    !< water saturation

			real (kind = GRID_SR)   :: r_dummy(2)
		END type

		!> persistent, scenario specific data on a cell
		type num_cell_data_pers
			real (kind = GRID_SR)   :: p(_DARCY_P_CELL_SIZE)
			real (kind = GRID_SR)   :: A_d(_DARCY_P_CELL_SIZE), d(_DARCY_P_CELL_SIZE), r(_DARCY_P_CELL_SIZE)

			real (kind = GRID_SR)   :: u(2, _DARCY_U_CELL_SIZE)			    !< velocity
			real (kind = GRID_SR)   :: saturation(_DARCY_FLOW_CELL_SIZE)    !< water saturation

			real (kind = GRID_SR)   :: base_permeability
			real (kind = GRID_SR)   :: permeability
		END type

		!*********************************************
		!Temporary Entity data (geometric association)
		!*********************************************

		!> temporary, scenario specific data on a node (deleted after each traversal)
		type num_node_data_temp
			real (kind = GRID_SR), DIMENSION(_DARCY_P_NODE_SIZE)		:: r
			real (kind = GRID_SR), DIMENSION(_DARCY_P_NODE_SIZE)		:: mat_diagonal

			real (kind = GRID_SR), DIMENSION(_DARCY_FLOW_NODE_SIZE)		:: flux
			real (kind = GRID_SR), DIMENSION(_DARCY_FLOW_NODE_SIZE)		:: volume
			logical 									                :: is_dirichlet_boundary(1)
		END type num_node_data_temp

		!> temporary, scenario specific data on an edge (deleted after each traversal)
		type num_edge_data_temp
			real (kind = GRID_SR), DIMENSION(_DARCY_P_EDGE_SIZE)		:: r
			real (kind = GRID_SR), DIMENSION(_DARCY_P_EDGE_SIZE)		:: mat_diagonal

			real (kind = GRID_SR), DIMENSION(_DARCY_FLOW_EDGE_SIZE)		:: flux
			real (kind = GRID_SR), DIMENSION(_DARCY_FLOW_EDGE_SIZE)		:: volume

			real (kind = GRID_SR)										:: r_dummy
		END type num_edge_data_temp

		!> temporary, scenario specific data on a cell (deleted after each traversal)
		type num_cell_data_temp
			real (kind = GRID_SR), DIMENSION(_DARCY_P_CELL_SIZE)		:: r
			real (kind = GRID_SR), DIMENSION(_DARCY_P_CELL_SIZE)		:: mat_diagonal

			real (kind = GRID_SR), DIMENSION(_DARCY_FLOW_CELL_SIZE)		:: flux
			real (kind = GRID_SR), DIMENSION(_DARCY_FLOW_CELL_SIZE)		:: volume

			real (kind = GRID_SR)										:: r_dummy
		END type num_cell_data_temp

		!*************************
		!Global data
		!*************************

		!> Base data type for the scenario configuration
		type num_global_data
			real (kind = GRID_SR)						:: r_time					!< simulation time
			real (kind = GRID_SR)						:: r_dt						!< time step
			real (kind = GRID_SR)						:: u_max			        !< maximum velocity
		END type
	END MODULE Darcy_data_types
#endif
