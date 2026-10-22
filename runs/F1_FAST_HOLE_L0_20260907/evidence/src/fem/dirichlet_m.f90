module dirichlet_m
  use kinds_m, only : dp
  use constants_m, only : VERIFICATION_TOL
  use status_m, only : SolverStatusType, clear_status, set_status, &
    STATUS_INPUT_ERROR, STATUS_INVALID_ARGUMENT
  use field_layout_m, only : NDOF_NODE, global_dof
  use model_types_m, only : StateType, LoadStateType
  use sparse_pattern_m, only : SparseAssemblyContextType, find_csr_position
  implicit none
  private

  type, public :: DirichletBCType
    integer :: field_id = 0
    integer, allocatable :: node_ids(:)
    real(dp) :: value_start = 0.0_dp
    real(dp) :: value_end = 0.0_dp
  end type DirichletBCType

  public :: project_dirichlet_state, apply_dirichlet_increment_rows, &
    build_fixed_dof_mask

contains

  subroutine project_dirichlet_state(state, boundary_conditions, load, node_count, status)
    type(StateType), intent(inout) :: state
    type(DirichletBCType), intent(in) :: boundary_conditions(:)
    type(LoadStateType), intent(in) :: load
    integer, intent(in) :: node_count
    type(SolverStatusType), intent(out) :: status
    logical, allocatable :: assigned(:)
    real(dp), allocatable :: assigned_value(:)
    integer :: boundary, entry, dof
    real(dp) :: prescribed

    call clear_status(status)
    if (.not. allocated(state%x) .or. size(state%x) /= NDOF_NODE*node_count) then
      call set_status(status, STATUS_INVALID_ARGUMENT, 'State size is incompatible with Dirichlet projection.')
      return
    end if
    allocate(assigned(size(state%x)),assigned_value(size(state%x)))
    assigned = .false.
    assigned_value = 0.0_dp
    do boundary = 1, size(boundary_conditions)
      call validate_boundary_set(boundary_conditions(boundary),node_count,status)
      if (.not. status%is_ok()) return
      prescribed = boundary_value(boundary_conditions(boundary),load)
      do entry = 1, size(boundary_conditions(boundary)%node_ids)
        dof = global_dof(boundary_conditions(boundary)%node_ids(entry), &
          boundary_conditions(boundary)%field_id)
        if (assigned(dof) .and. abs(assigned_value(dof)-prescribed) > VERIFICATION_TOL) then
          call set_status(status, STATUS_INPUT_ERROR, 'Conflicting Dirichlet values on the same DOF.')
          return
        end if
        assigned(dof) = .true.
        assigned_value(dof) = prescribed
        state%x(dof) = prescribed
      end do
    end do
  end subroutine project_dirichlet_state

  subroutine apply_dirichlet_increment_rows(context, boundary_conditions, residual, &
                                            jacobian_values, status)
    type(SparseAssemblyContextType), intent(in) :: context
    type(DirichletBCType), intent(in) :: boundary_conditions(:)
    real(dp), intent(inout) :: residual(:)
    real(dp), intent(inout) :: jacobian_values(:)
    type(SolverStatusType), intent(out) :: status
    logical, allocatable :: fixed(:)
    integer :: row, position, diagonal_position

    call clear_status(status)
    if (size(residual) /= context%number_of_equations .or. &
        size(jacobian_values) /= context%number_of_nonzeros) then
      call set_status(status, STATUS_INVALID_ARGUMENT, 'Dirichlet row application size mismatch.')
      return
    end if
    call build_fixed_dof_mask(boundary_conditions, context%number_of_equations/NDOF_NODE, &
      fixed, status)
    if (.not. status%is_ok()) return
    do row = 1, context%number_of_equations
      if (.not. fixed(row)) cycle
      residual(row) = 0.0_dp
      do position = context%ia(row), context%ia(row+1)-1
        jacobian_values(position) = 0.0_dp
      end do
      diagonal_position = find_csr_position(context,row,row)
      if (diagonal_position == 0) then
        call set_status(status, STATUS_INVALID_ARGUMENT, 'CSR row lacks a Dirichlet diagonal entry.')
        return
      end if
      jacobian_values(diagonal_position) = 1.0_dp
    end do
  end subroutine apply_dirichlet_increment_rows

  subroutine build_fixed_dof_mask(boundary_conditions, node_count, fixed, status)
    type(DirichletBCType), intent(in) :: boundary_conditions(:)
    integer, intent(in) :: node_count
    logical, allocatable, intent(out) :: fixed(:)
    type(SolverStatusType), intent(out) :: status
    integer :: boundary, entry, dof

    call clear_status(status)
    allocate(fixed(NDOF_NODE*node_count))
    fixed = .false.
    do boundary = 1, size(boundary_conditions)
      call validate_boundary_set(boundary_conditions(boundary),node_count,status)
      if (.not. status%is_ok()) return
      do entry = 1, size(boundary_conditions(boundary)%node_ids)
        dof = global_dof(boundary_conditions(boundary)%node_ids(entry), &
          boundary_conditions(boundary)%field_id)
        fixed(dof) = .true.
      end do
    end do
  end subroutine build_fixed_dof_mask

  subroutine validate_boundary_set(boundary_condition,node_count,status)
    type(DirichletBCType), intent(in) :: boundary_condition
    integer, intent(in) :: node_count
    type(SolverStatusType), intent(out) :: status

    call clear_status(status)
    if (boundary_condition%field_id < 1 .or. boundary_condition%field_id > NDOF_NODE) then
      call set_status(status, STATUS_INPUT_ERROR, 'Dirichlet field identifier is invalid.')
      return
    end if
    if (.not. allocated(boundary_condition%node_ids) .or. &
        size(boundary_condition%node_ids) < 1) then
      call set_status(status, STATUS_INPUT_ERROR, 'Dirichlet node set is empty.')
      return
    end if
    if (any(boundary_condition%node_ids < 1) .or. &
        any(boundary_condition%node_ids > node_count)) then
      call set_status(status, STATUS_INPUT_ERROR, 'Dirichlet node identifier is outside the mesh.')
    end if
  end subroutine validate_boundary_set

  pure real(dp) function boundary_value(boundary_condition,load) result(value)
    type(DirichletBCType), intent(in) :: boundary_condition
    type(LoadStateType), intent(in) :: load
    value = boundary_condition%value_start + load%voltage_factor * &
      (boundary_condition%value_end-boundary_condition%value_start)
  end function boundary_value

end module dirichlet_m
