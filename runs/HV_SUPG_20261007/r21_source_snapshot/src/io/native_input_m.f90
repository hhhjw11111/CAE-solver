module native_input_m
  use kinds_m, only : dp
  use status_m, only : SolverStatusType, clear_status, set_status, STATUS_INPUT_ERROR
  use input_schema_m, only : CaseConfigType, validate_case_config, MODEL_M2
  use field_layout_m, only : NDOF_NODE
  use mesh_types_m, only : MeshType, HEX20_NNODE
  implicit none
  private

  public :: read_native_case

contains

  subroutine read_native_case(case_path, config, mesh, status)
    character(len=*), intent(in) :: case_path
    type(CaseConfigType), intent(out) :: config
    type(MeshType), intent(out) :: mesh
    type(SolverStatusType), intent(out) :: status
    integer :: schema_version, unit_id, io_status,max_newton,newton_stagnation_window
    integer :: continuation_maximum_attempts,continuation_maximum_newton_per_attempt
    integer :: continuation_fast_iterations,continuation_slow_iterations
    logical :: continuation_full_joule_during_voltage
    logical :: transport_supg
    logical :: half_y_symmetry
    logical :: supg_directional_tau
    logical :: global_log_step_limit
    logical :: continuation_secant_predictor
    real(dp) :: supg_scale
    character(len=128) :: case_id
    character(len=64) :: model_level,problem_type,material_profile
    character(len=512) :: mesh_file, mesh_path,restart_nodes_file
    character(len=512) :: io_message
    real(dp) :: target_voltage,target_voltage_factor,target_joule_factor,hole_radius
    real(dp) :: restart_voltage_factor,restart_joule_factor
    real(dp) :: field_scales(NDOF_NODE),residual_tolerance,update_tolerance
    real(dp) :: armijo_c1,backtrack_factor,minimum_step,maximum_log_density_update
    real(dp) :: continuation_voltage_checkpoint_spacing
    real(dp) :: continuation_initial_step,continuation_minimum_step
    real(dp) :: continuation_maximum_step,continuation_cutback_factor
    real(dp) :: continuation_growth_factor
    namelist /case_control/ schema_version, case_id, model_level, mesh_file
    namelist /problem_control/ problem_type,material_profile,target_voltage, &
      target_voltage_factor,target_joule_factor,hole_radius,restart_nodes_file, &
      restart_voltage_factor,restart_joule_factor,half_y_symmetry
    namelist /solver_control/ field_scales,max_newton,newton_stagnation_window,residual_tolerance, &
      continuation_voltage_checkpoint_spacing, &
      update_tolerance,armijo_c1,backtrack_factor,minimum_step, &
      maximum_log_density_update,continuation_initial_step, &
      continuation_full_joule_during_voltage, &
      continuation_minimum_step,continuation_maximum_step, &
      continuation_cutback_factor,continuation_growth_factor, &
      continuation_maximum_attempts,continuation_maximum_newton_per_attempt, &
      continuation_fast_iterations,continuation_slow_iterations,transport_supg,supg_scale,global_log_step_limit, &
      continuation_secant_predictor,supg_directional_tau

    call clear_status(status)
    transport_supg=.false.;supg_scale=1.0_dp
    half_y_symmetry=.false.
    supg_directional_tau=.false.
    global_log_step_limit=.false.
    continuation_secant_predictor=.false.
    schema_version = 0
    case_id = 'UNSET'
    model_level = 'UNSET'
    mesh_file = ''
    problem_type='UNSET'
    material_profile='UNSET'
    restart_nodes_file=''
    target_voltage=0.0_dp
    target_voltage_factor=1.0_dp
    target_joule_factor=1.0_dp
    hole_radius=0.0_dp
    restart_voltage_factor=0.0_dp
    restart_joule_factor=0.0_dp
    field_scales=[1.0e-10_dp,1.0e-10_dp,1.0e-10_dp,10.0_dp, &
      1.0e-1_dp,1.0e-1_dp,10.0_dp]
    max_newton=40
    newton_stagnation_window=0
    continuation_voltage_checkpoint_spacing=0.0_dp
    residual_tolerance=1.0e-8_dp
    update_tolerance=1.0e-8_dp
    armijo_c1=1.0e-4_dp
    backtrack_factor=0.5_dp
    minimum_step=1.0e-8_dp
    maximum_log_density_update=1.0_dp
    continuation_full_joule_during_voltage=.false.
    continuation_initial_step=0.25_dp
    continuation_minimum_step=0.03125_dp
    continuation_maximum_step=0.5_dp
    continuation_cutback_factor=0.5_dp
    continuation_growth_factor=2.0_dp
    continuation_maximum_attempts=100
    continuation_maximum_newton_per_attempt=20
    continuation_fast_iterations=4
    continuation_slow_iterations=10

    open(newunit=unit_id, file=trim(case_path), status='old', action='read', &
      iostat=io_status, iomsg=io_message)
    if (io_status /= 0) then
      call set_status(status, STATUS_INPUT_ERROR, 'Cannot open case file: ' // trim(io_message))
      return
    end if
    read(unit_id, nml=case_control, iostat=io_status, iomsg=io_message)
    if (io_status /= 0) then
      close(unit_id)
      call set_status(status, STATUS_INPUT_ERROR, 'Cannot read case_control namelist: ' // trim(io_message))
      return
    end if
    if (trim(model_level) == MODEL_M2) then
      rewind(unit_id)
      read(unit_id,nml=problem_control,iostat=io_status,iomsg=io_message)
      if (io_status /= 0) then
        close(unit_id)
        call set_status(status,STATUS_INPUT_ERROR, &
          'Cannot read problem_control namelist: '//trim(io_message))
        return
      end if
      rewind(unit_id)
      read(unit_id,nml=solver_control,iostat=io_status,iomsg=io_message)
      if (io_status /= 0) then
        close(unit_id)
        call set_status(status,STATUS_INPUT_ERROR, &
          'Cannot read solver_control namelist: '//trim(io_message))
        return
      end if
    end if
    close(unit_id)

    config%schema_version = schema_version
    config%case_id = trim(case_id)
    config%model_level = trim(model_level)
    config%mesh_file = trim(mesh_file)
    config%case_file = trim(case_path)
    config%problem_type=trim(problem_type)
    config%material_profile=trim(material_profile)
    config%target_voltage=target_voltage
    config%target_voltage_factor=target_voltage_factor
    config%target_joule_factor=target_joule_factor
    config%hole_radius=hole_radius
    if (len_trim(restart_nodes_file) > 0) then
      config%restart_nodes_file=resolve_sibling_path(case_path,restart_nodes_file)
    else
      config%restart_nodes_file=''
    end if
    config%restart_voltage_factor=restart_voltage_factor
    config%restart_joule_factor=restart_joule_factor
    config%field_scales=field_scales
    config%transport_supg=transport_supg;config%supg_scale=supg_scale
    config%half_y_symmetry=half_y_symmetry
    config%supg_directional_tau=supg_directional_tau
    config%nonlinear_options%max_newton=max_newton
    config%nonlinear_options%stagnation_window=newton_stagnation_window
    config%continuation_options%voltage_checkpoint_spacing=continuation_voltage_checkpoint_spacing
    config%nonlinear_options%global_log_step_limit=global_log_step_limit
    config%nonlinear_options%residual_tolerance=residual_tolerance
    config%nonlinear_options%update_tolerance=update_tolerance
    config%nonlinear_options%armijo_c1=armijo_c1
    config%nonlinear_options%backtrack_factor=backtrack_factor
    config%nonlinear_options%minimum_step=minimum_step
    config%nonlinear_options%maximum_log_density_update=maximum_log_density_update
    config%continuation_options%initial_step=continuation_initial_step
    config%continuation_options%secant_predictor=continuation_secant_predictor
    config%continuation_options%full_joule_during_voltage= &
      continuation_full_joule_during_voltage
    config%continuation_options%minimum_step=continuation_minimum_step
    config%continuation_options%maximum_step=continuation_maximum_step
    config%continuation_options%cutback_factor=continuation_cutback_factor
    config%continuation_options%growth_factor=continuation_growth_factor
    config%continuation_options%maximum_attempts=continuation_maximum_attempts
    config%continuation_options%maximum_newton_per_attempt= &
      continuation_maximum_newton_per_attempt
    config%continuation_options%fast_iteration_threshold=continuation_fast_iterations
    config%continuation_options%slow_iteration_threshold=continuation_slow_iterations
    call validate_case_config(config, status)
    if (.not. status%is_ok()) return

    mesh_path = resolve_sibling_path(case_path, mesh_file)
    call read_mesh_file(trim(mesh_path), mesh, status)
  end subroutine read_native_case

  subroutine read_mesh_file(mesh_path, mesh, status)
    character(len=*), intent(in) :: mesh_path
    type(MeshType), intent(out) :: mesh
    type(SolverStatusType), intent(out) :: status
    integer :: unit_id, io_status, node_count, element_count
    integer :: row, identifier, material_identifier
    integer :: nodes(HEX20_NNODE)
    real(dp) :: x, y, z
    character(len=512) :: io_message

    call clear_status(status)
    open(newunit=unit_id, file=trim(mesh_path), status='old', action='read', &
      iostat=io_status, iomsg=io_message)
    if (io_status /= 0) then
      call set_status(status, STATUS_INPUT_ERROR, 'Cannot open mesh file: ' // trim(io_message))
      return
    end if

    read(unit_id,*,iostat=io_status,iomsg=io_message) node_count, element_count
    if (io_status /= 0 .or. node_count < HEX20_NNODE .or. element_count < 1) then
      close(unit_id)
      call set_status(status, STATUS_INPUT_ERROR, 'Invalid mesh header; expected node_count element_count.')
      return
    end if
    allocate(mesh%coordinates(3,node_count))
    allocate(mesh%connectivity(HEX20_NNODE,element_count))
    allocate(mesh%material_id(element_count))

    do row = 1, node_count
      read(unit_id,*,iostat=io_status,iomsg=io_message) identifier, x, y, z
      if (io_status /= 0 .or. identifier /= row) then
        close(unit_id)
        call set_status(status, STATUS_INPUT_ERROR, &
          'Node records must be contiguous and ordered from 1.')
        return
      end if
      mesh%coordinates(:,row) = [x,y,z]
    end do

    do row = 1, element_count
      read(unit_id,*,iostat=io_status,iomsg=io_message) identifier, material_identifier, nodes
      if (io_status /= 0 .or. identifier /= row) then
        close(unit_id)
        call set_status(status, STATUS_INPUT_ERROR, &
          'Element records must be contiguous and ordered from 1.')
        return
      end if
      mesh%material_id(row) = material_identifier
      mesh%connectivity(:,row) = nodes
    end do
    close(unit_id)
  end subroutine read_mesh_file

  function resolve_sibling_path(parent_file, child_file) result(resolved)
    character(len=*), intent(in) :: parent_file, child_file
    character(len=512) :: resolved
    integer :: slash_position

    if (is_absolute_path(child_file)) then
      resolved = trim(child_file)
      return
    end if
    slash_position = max(scan(trim(parent_file), '/', back=.true.), &
                         scan(trim(parent_file), achar(92), back=.true.))
    if (slash_position > 0) then
      resolved = parent_file(:slash_position) // trim(child_file)
    else
      resolved = trim(child_file)
    end if
  end function resolve_sibling_path

  logical function is_absolute_path(path)
    character(len=*), intent(in) :: path
    integer :: length_path

    length_path = len_trim(path)
    is_absolute_path = .false.
    if (length_path >= 2) is_absolute_path = path(2:2) == ':'
    if (length_path >= 1) then
      is_absolute_path = is_absolute_path .or. path(1:1) == '/' .or. path(1:1) == achar(92)
    end if
  end function is_absolute_path

end module native_input_m
