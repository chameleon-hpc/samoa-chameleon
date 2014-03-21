Config.o: Tools_log.o M_kracken.o
SFC_data_types.o: Config.o Tools_log.o Heat_Equation/Heat_Eq_data_types.o Darcy/Darcy_data_types.o SWE/SWE_data_types.o Tests/Tests_data_types.o Generic/Generic_data_types.o Flash/FLASH_data_types.o #NUMA/NUMA_data_types.o
SFC_edge_traversal.o: SFC_data_types.o
SFC_node_traversal.o: SFC_edge_traversal.o
SFC_traversal.o: SFC_node_traversal.o Heat_Equation/Heat_Eq.o Darcy/Darcy.o Tests/Tests.o SWE/SWE.o Generic/Generic.o Flash/FLASH.o #NUMA/NUMA.o
SFC_main.o: Config.o SFC_traversal.o
Tools_local_function_space_base.o: SFC_data_types.o
Tools_noise.o: SFC_data_types.o
Tools_stack_base.o: SFC_data_types.o

Conformity/Conformity.o: SFC_edge_traversal.o

geoclaw/c_bind_riemannsolvers.o: geoclaw/riemannsolvers.o geoclaw/riemannsolvers_sp.o

Tests/Tests.o: Tests/Tests_data_types.o Tests/Tests_basis_functions.o Tests/Tests_initialize.o Tests/Tests_node_dummy_traversal.o Tests/Tests_flops_traversal.o Tests/Tests_memory_traversal.o Tests/Tests_consistency_traversal.o Samoa/Samoa.o
Tests/Tests_basis_functions.o: SFC_data_types.o
Tests/Tests_initialize.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o
Tests/Tests_node_dummy_traversal.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o
Tests/Tests_consistency_traversal.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o
Tests/Tests_flops_traversal.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o
Tests/Tests_memory_traversal.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o

Darcy/Darcy.o: Darcy/Darcy_data_types.o Solver/LinearSolver.o Darcy/Darcy_initialize.o Darcy/Darcy_node_dummy.o Darcy/Darcy_output.o Darcy/Darcy_xml_output.o Darcy/Darcy_grad_p.o Darcy/Darcy_laplace_jacobi.o Darcy/Darcy_laplace_cg.o Darcy/Darcy_laplace_pipecg.o Darcy/Darcy_transport_eq.o Darcy/Darcy_permeability.o Darcy/Darcy_adapt.o Samoa/Samoa.o
Darcy/Darcy_local_function_spaces.o: SFC_data_types.o Tools_grid_variable.f90 Tools_local_function_space.f90
Darcy/Darcy_basis.o: SFC_data_types.o Darcy/Darcy_local_function_spaces.o Samoa/Samoa.o
Darcy/Darcy_adapt.o: SFC_generic_adaptive_traversal.f90 Conformity/Conformity.o Darcy/Darcy_local_function_spaces.o Darcy/Darcy_basis.o Samoa/Samoa.o Darcy/Darcy_initialize.o
Darcy/Darcy_node_dummy.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o
Darcy/Darcy_grad_p.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o Darcy/Darcy_local_function_spaces.o Darcy/Darcy_basis.o Samoa/Samoa.o
Darcy/Darcy_initialize.o: SFC_generic_traversal_ringbuffer.f90 Tools_noise.o SFC_node_traversal.o Darcy/Darcy_basis.o Samoa/Samoa.o
Darcy/Darcy_laplace_cg.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o Darcy/Darcy_basis.o Samoa/Samoa.o
Darcy/Darcy_laplace_pipecg.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o Darcy/Darcy_basis.o Samoa/Samoa.o
Darcy/Darcy_laplace_jacobi.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o Darcy/Darcy_basis.o Samoa/Samoa.o
Darcy/Darcy_output.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o Darcy/Darcy_basis.o Samoa/Samoa.o LIB_VTK_IO.o
Darcy/Darcy_xml_output.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o Darcy/Darcy_basis.o Samoa/Samoa.o LIB_VTK_IO.o
Darcy/Darcy_transport_eq.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o Darcy/Darcy_basis.o Samoa/Samoa.o
Darcy/Darcy_permeability.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o Darcy/Darcy_basis.o Samoa/Samoa.o

Heat_Equation/Heat_Eq.o: Heat_Equation/Heat_Eq_data_types.o Heat_Equation/Heat_Eq_basis.o Heat_Equation/Heat_Eq_initialize.o Heat_Equation/Heat_Eq_output.o Heat_Equation/Heat_Eq_xml_output.o Heat_Equation/Heat_Eq_euler_timestep.o Heat_Equation/Heat_Eq_midpoint_timestep.o Heat_Equation/Heat_Eq_heun_timestep.o Heat_Equation/Heat_Eq_adapt.o Samoa/Samoa.o
Heat_Equation/Heat_Eq_local_function_spaces.o: SFC_data_types.o
Heat_Equation/Heat_Eq_basis.o: SFC_data_types.o Heat_Equation/Heat_Eq_local_function_spaces.o Samoa/Samoa.o
Heat_Equation/Heat_Eq_adapt.o: SFC_generic_adaptive_traversal.f90 Conformity/Conformity.o Heat_Equation/Heat_Eq_basis.o Samoa/Samoa.o
Heat_Equation/Heat_Eq_euler_timestep.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o Heat_Equation/Heat_Eq_basis.o Samoa/Samoa.o
Heat_Equation/Heat_Eq_heun_timestep.o: SFC_generic_traversal_ringbuffer.f90 Heat_Equation/Heat_Eq_euler_timestep.o
Heat_Equation/Heat_Eq_initialize.o: SFC_generic_traversal_ringbuffer.f90 Tools_noise.o SFC_node_traversal.o Heat_Equation/Heat_Eq_basis.o Samoa/Samoa.o Heat_Equation/Heat_Eq_local_function_spaces.o
Heat_Equation/Heat_Eq_midpoint_timestep.o: SFC_generic_traversal_ringbuffer.f90 Heat_Equation/Heat_Eq_euler_timestep.o
Heat_Equation/Heat_Eq_output.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o Heat_Equation/Heat_Eq_basis.o Heat_Equation/Heat_Eq_local_function_spaces.o Samoa/Samoa.o LIB_VTK_IO.o
Heat_Equation/Heat_Eq_xml_output.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o Heat_Equation/Heat_Eq_basis.o Heat_Equation/Heat_Eq_local_function_spaces.o Samoa/Samoa.o LIB_VTK_IO.o

SWE/SWE.o: SWE/SWE_data_types.o SWE/SWE_displace.o SWE/SWE_basis.o SWE/SWE_initialize.o SWE/SWE_output.o SWE/SWE_xml_output.o SWE/SWE_euler_timestep.o SWE/SWE_adapt.o Samoa/Samoa.o SWE/SWE_ascii_output.o
SWE/SWE_local_function_spaces.o: SFC_data_types.o
SWE/SWE_basis.o: SFC_data_types.o SWE/SWE_local_function_spaces.o Samoa/Samoa.o
SWE/SWE_adapt.o: SFC_generic_adaptive_traversal.f90 Conformity/Conformity.o SWE/SWE_basis.o Samoa/Samoa.o SWE/SWE_euler_timestep.o SWE/SWE_initialize.o
SWE/SWE_euler_timestep.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o SWE/SWE_basis.o Samoa/Samoa.o geoclaw/c_bind_riemannsolvers.o
SWE/SWE_initialize.o: SFC_generic_traversal_ringbuffer.f90 SWE/SWE_euler_timestep.o Tools_noise.o SFC_node_traversal.o SWE/SWE_basis.o Samoa/Samoa.o SWE/SWE_local_function_spaces.o
SWE/SWE_displace.o: SWE/SWE_initialize.o SFC_generic_traversal_ringbuffer.f90 SWE/SWE_euler_timestep.o Tools_noise.o SFC_node_traversal.o SWE/SWE_basis.o Samoa/Samoa.o SWE/SWE_local_function_spaces.o
SWE/SWE_output.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o SWE/SWE_basis.o SWE/SWE_local_function_spaces.o SWE/SWE_euler_timestep.o Samoa/Samoa.o LIB_VTK_IO.o
SWE/SWE_xml_output.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o SWE/SWE_basis.o SWE/SWE_local_function_spaces.o SWE/SWE_euler_timestep.o Samoa/Samoa.o LIB_VTK_IO.o
SWE/SWE_ascii_output.o: SWE/ascii_output.o SFC_edge_traversal.o LIB_VTK_IO.o SFC_edge_traversal.f90 Samoa/Samoa.o SWE/SWE_euler_timestep.o

#NUMA/NUMA.o: NUMA/NUMA_data_types.o NUMA/NUMA_basis.o NUMA/NUMA_initialize.o NUMA/NUMA_output.o NUMA/NUMA_xml_output.o NUMA/NUMA_euler_timestep.o NUMA/NUMA_adapt.o Samoa/Samoa.o
#NUMA/NUMA_local_function_spaces.o: SFC_data_types.o
#NUMA/NUMA_basis.o: SFC_data_types.o NUMA/NUMA_local_function_spaces.o Samoa/Samoa.o
#NUMA/NUMA_adapt.o: SFC_generic_adaptive_traversal.f90 Conformity/Conformity.o NUMA/NUMA_basis.o Samoa/Samoa.o NUMA/NUMA_euler_timestep.o NUMA/NUMA_initialize.o
#NUMA/NUMA_euler_timestep.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o NUMA/NUMA_basis.o NUMA/NUMA_constants.o Samoa/Samoa.o geoclaw/c_bind_riemannsolvers.o
#NUMA/NUMA_initialize.o: SFC_generic_traversal_ringbuffer.f90 NUMA/NUMA_euler_timestep.o Tools_noise.o SFC_node_traversal.o NUMA/NUMA_basis.o Samoa/Samoa.o NUMA/NUMA_local_function_spaces.o
#NUMA/NUMA_output.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o NUMA/NUMA_basis.o NUMA/NUMA_local_function_spaces.o NUMA/NUMA_euler_timestep.o Samoa/Samoa.o LIB_VTK_IO.o
#NUMA/NUMA_xml_output.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o NUMA/NUMA_basis.o NUMA/NUMA_local_function_spaces.o NUMA/NUMA_euler_timestep.o Samoa/Samoa.o LIB_VTK_IO.o
#NUMA/NUMA_constants.o: NUMA/NUMA_data_types.o


Flash/FLASH.o: Flash/FLASH_data_types.o Flash/FLASH_basis.o Flash/FLASH_initialize.o Flash/FLASH_output.o Flash/FLASH_xml_output.o Flash/FLASH_euler_timestep.o Flash/FLASH_adapt.o Samoa/Samoa.o
Flash/FLASH_local_function_spaces.o: SFC_data_types.o
Flash/FLASH_basis.o: SFC_data_types.o Flash/FLASH_local_function_spaces.o Samoa/Samoa.o
Flash/FLASH_adapt.o: SFC_generic_adaptive_traversal.f90 Conformity/Conformity.o Flash/FLASH_basis.o Samoa/Samoa.o Flash/FLASH_euler_timestep.o Flash/FLASH_initialize.o
Flash/FLASH_euler_timestep.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o Flash/FLASH_basis.o Samoa/Samoa.o Flash/FLASH_dg_element.o
Flash/FLASH_initialize.o: SFC_generic_traversal_ringbuffer.f90 Flash/FLASH_euler_timestep.o Tools_noise.o SFC_node_traversal.o Flash/FLASH_basis.o Samoa/Samoa.o Flash/FLASH_local_function_spaces.o
Flash/FLASH_output.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o Flash/FLASH_basis.o Flash/FLASH_local_function_spaces.o Flash/FLASH_euler_timestep.o Samoa/Samoa.o LIB_VTK_IO.o
Flash/FLASH_xml_output.o: SFC_generic_traversal_ringbuffer.f90 SFC_node_traversal.o Flash/FLASH_basis.o Flash/FLASH_local_function_spaces.o Flash/FLASH_euler_timestep.o Samoa/Samoa.o LIB_VTK_IO.o
Flash/FLASH_dg_element.o: SFC_data_types.o

Generic/Generic.o: Generic/Generic_data_types.o  SFC_node_traversal.o Generic/Generic_initialize.o Generic/Generic_template.o Generic/Generic_adapt_template.o
Generic/Generic_initialize.o: SFC_generic_traversal_ringbuffer.f90 SFC_edge_traversal.o Samoa/Samoa.o
Generic/Generic_template.o: SFC_generic_traversal_ringbuffer.f90 SFC_edge_traversal.o Samoa/Samoa.o
Generic/Generic_adapt_template.o: SFC_generic_adaptive_traversal.f90 Conformity/Conformity.o Samoa/Samoa.o

Samoa/Samoa.o: SFC_data_types.o Samoa/Tools_quadrature_rule_base.o
Samoa/Tools_quadrature_rule_base.o: SFC_data_types.o

Solver/LinearSolver.o: SFC_edge_traversal.o
