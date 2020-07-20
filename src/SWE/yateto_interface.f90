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
    integer(C_INT),value :: i0
  end subroutine yateto_predictor_execute
  end interface
  
end module yateto_interface
