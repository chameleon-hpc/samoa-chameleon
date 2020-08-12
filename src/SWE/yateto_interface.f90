module yateto_interface
  use iso_c_binding
  implicit none

  interface
     !subroutine yateto_predictor_execute(P, F, dtdx, i0) bind(C, name="yateto_predictor_execute")
     subroutine yateto_predictor_execute(P, F, Q0, S, dtdx , i0) bind(C, name="yateto_predictor_execute")
       use iso_c_binding
       real(c_double),dimension(_SWE_DG_DOFS,_SWE_DG_ORDER  ,3)   :: P
       real(c_double),dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1,2,3) :: F
       real(c_double),dimension(_SWE_DG_DOFS,3)                   :: Q0
       real(c_double),dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1,2)   :: S
       real(c_double),dimension(1)   :: dtdx    
       integer(c_int),value :: i0
     end subroutine yateto_predictor_execute

     subroutine yateto_compute_source_execute(S2, W)  bind(C, name="yateto_compute_source_execute")
       use iso_c_binding
       real(c_double),dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1,2) :: S2
       real(c_double),dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1)   :: W
     end subroutine yateto_compute_source_execute


     subroutine yateto_volume_execute(U, F, S, i0)  bind(C, name="yateto_volume_execute")
       use iso_c_binding
       real(c_double),dimension(_SWE_DG_DOFS,3)                   :: U
       real(c_double),dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1,2,3) :: F
       real(c_double),dimension(_SWE_DG_DOFS,_SWE_DG_ORDER+1,2)   :: S
       integer(c_int),value :: i0
     end subroutine yateto_volume_execute
  end interface
  
end module yateto_interface