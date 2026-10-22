module scaling_m
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use kinds_m, only : dp
  use status_m, only : SolverStatusType, clear_status, set_status, STATUS_INVALID_ARGUMENT
  use field_layout_m, only : NDOF_NODE
  use sparse_pattern_m, only : SparseAssemblyContextType
  implicit none
  private

  public :: build_variable_scaling, equilibrate_csr_system, unscale_solution

contains

  subroutine build_variable_scaling(node_count,field_scales,variable_scales,status)
    integer, intent(in) :: node_count
    real(dp), intent(in) :: field_scales(NDOF_NODE)
    real(dp), allocatable, intent(out) :: variable_scales(:)
    type(SolverStatusType), intent(out) :: status
    integer :: equation,field

    call clear_status(status)
    if (node_count < 1 .or. any(field_scales <= 0.0_dp) .or. &
        .not. all(ieee_is_finite(field_scales))) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'Variable field scales must be finite and positive.')
      return
    end if
    allocate(variable_scales(NDOF_NODE*node_count))
    do equation=1,size(variable_scales)
      field=mod(equation-1,NDOF_NODE)+1
      variable_scales(equation)=field_scales(field)
    end do
  end subroutine build_variable_scaling

  subroutine equilibrate_csr_system(context,values,rhs,variable_scales, &
                                    scaled_values,scaled_rhs,residual_scales,status)
    type(SparseAssemblyContextType), intent(in) :: context
    real(dp), intent(in) :: values(:),rhs(:),variable_scales(:)
    real(dp), allocatable, intent(out) :: scaled_values(:),scaled_rhs(:),residual_scales(:)
    type(SolverStatusType), intent(out) :: status
    integer :: row,position,column
    real(dp) :: row_scale

    call clear_status(status)
    if (size(values) /= context%number_of_nonzeros .or. &
        size(rhs) /= context%number_of_equations .or. &
        size(variable_scales) /= context%number_of_equations) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'Scaling array size mismatch.')
      return
    end if
    if (any(variable_scales <= 0.0_dp) .or. .not. all(ieee_is_finite(variable_scales)) .or. &
        .not. all(ieee_is_finite(values)) .or. .not. all(ieee_is_finite(rhs))) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'Scaling input contains non-finite or non-positive data.')
      return
    end if
    allocate(scaled_values(size(values)),scaled_rhs(size(rhs)), &
      residual_scales(context%number_of_equations))
    do row=1,context%number_of_equations
      row_scale=0.0_dp
      do position=context%ia(row),context%ia(row+1)-1
        column=context%ja(position)
        row_scale=max(row_scale,abs(values(position))*variable_scales(column))
      end do
      if (.not. ieee_is_finite(row_scale) .or. row_scale <= tiny(1.0_dp)) then
        call set_status(status,STATUS_INVALID_ARGUMENT,'Cannot equilibrate a zero or non-finite CSR row.')
        return
      end if
      residual_scales(row)=row_scale
      scaled_rhs(row)=rhs(row)/row_scale
      do position=context%ia(row),context%ia(row+1)-1
        column=context%ja(position)
        scaled_values(position)=values(position)*variable_scales(column)/row_scale
      end do
    end do
  end subroutine equilibrate_csr_system

  subroutine unscale_solution(variable_scales,scaled_solution,solution,status)
    real(dp), intent(in) :: variable_scales(:),scaled_solution(:)
    real(dp), intent(out) :: solution(:)
    type(SolverStatusType), intent(out) :: status

    call clear_status(status)
    if (size(variable_scales) /= size(scaled_solution) .or. &
        size(solution) /= size(scaled_solution)) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'Solution scaling array size mismatch.')
      return
    end if
    solution=variable_scales*scaled_solution
  end subroutine unscale_solution

end module scaling_m
