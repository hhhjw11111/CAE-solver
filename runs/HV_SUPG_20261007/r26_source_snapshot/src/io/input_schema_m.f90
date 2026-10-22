module input_schema_m
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use kinds_m, only : dp
  use status_m, only : SolverStatusType, clear_status, set_status, STATUS_INPUT_ERROR
  use field_layout_m, only : NDOF_NODE
  use model_types_m, only : NonlinearOptionsType, ContinuationOptionsType
  implicit none
  private

  integer, parameter, public :: NATIVE_SCHEMA_VERSION = 2
  character(len=*), parameter, public :: MODEL_M0 = 'M0_ANALYTIC'
  character(len=*), parameter, public :: MODEL_M1 = 'M1_TWO_CARRIER_ISOTHERMAL'
  character(len=*), parameter, public :: MODEL_M2 = 'M2_THERMO_TWO_CARRIER'
  character(len=*), parameter, public :: PROBLEM_REFERENCE_BAR = 'REFERENCE_BAR'
  character(len=*), parameter, public :: PROBLEM_PLAIN_PLATE = 'PLAIN_PLATE_CONTROL'
  character(len=*), parameter, public :: PROBLEM_CENTER_HOLE_PLATE = 'CENTER_HOLE_PLATE'
  character(len=*), parameter, public :: MATERIAL_TPS_G08_REFERENCE = 'TPS_G08_REFERENCE'

  type, public :: CaseConfigType
    integer :: schema_version = 0
    character(len=128) :: case_id = 'UNSET'
    character(len=64) :: model_level = 'UNSET'
    character(len=512) :: mesh_file = ''
    character(len=512) :: case_file = ''
    character(len=64) :: problem_type = 'UNSET'
    character(len=64) :: material_profile = 'UNSET'
    character(len=512) :: restart_nodes_file = ''
    real(dp) :: target_voltage = 0.0_dp
    real(dp) :: target_voltage_factor = 1.0_dp
    real(dp) :: target_joule_factor = 1.0_dp
    real(dp) :: hole_radius = 0.0_dp
    real(dp) :: restart_voltage_factor = 0.0_dp
    real(dp) :: restart_joule_factor = 0.0_dp
    real(dp) :: field_scales(NDOF_NODE) = [1.0e-10_dp,1.0e-10_dp,1.0e-10_dp, &
      10.0_dp,1.0e-1_dp,1.0e-1_dp,10.0_dp]
    type(NonlinearOptionsType) :: nonlinear_options
    type(ContinuationOptionsType) :: continuation_options
    logical :: transport_supg = .false.
    logical :: transport_log_weighted = .false.
    integer :: log_weighted_quadrature_order = 3
    logical :: half_y_symmetry = .false.
    logical :: supg_directional_tau = .false.
    real(dp) :: supg_scale = 1.0_dp
  end type CaseConfigType

  public :: validate_case_config

contains

  subroutine validate_case_config(config, status)
    type(CaseConfigType), intent(in) :: config
    type(SolverStatusType), intent(out) :: status

    call clear_status(status)
    if((config%log_weighted_quadrature_order/=3.and.config%log_weighted_quadrature_order/=5).or. &
        (config%log_weighted_quadrature_order/=3.and..not.config%transport_log_weighted)) then
      call set_status(status,STATUS_INPUT_ERROR,'Quadrature order 5 is restricted to log-weighted transport; allowed 3 or 5.')
      return
    end if
    if(config%transport_log_weighted) then
      if(config%transport_supg.or.(trim(config%material_profile)/='ZNO_F1_NT'.and. &
          trim(config%material_profile)/='ZNO_F1_N0')) then
        call set_status(status,STATUS_INPUT_ERROR,'Log-weighted transport requires nonlinear F1 NT/N0 without SUPG.')
        return
      end if
    end if
    if(config%half_y_symmetry) then
      if(trim(config%problem_type)/=PROBLEM_CENTER_HOLE_PLATE.or. &
         (trim(config%material_profile)/='ZNO_F1_NT'.and. &
          trim(config%material_profile)/='ZNO_F1_N0'.and.trim(config%material_profile)/='ZNO_F1_LQ')) then
        call set_status(status,STATUS_INPUT_ERROR,'Y-half symmetry requires the explicitly symmetric F1 ZnO hole profile.')
        return
      end if
    end if
    if (.not.ieee_is_finite(config%supg_scale) .or. config%supg_scale <= 0.0_dp) then
      call set_status(status,STATUS_INPUT_ERROR,'supg_scale must be finite and positive.');return
    end if
    if (config%transport_supg) then
      if (trim(config%material_profile)/='ZNO_F1_NT' .and. &
          trim(config%material_profile)/='ZNO_F1_N0') then
        call set_status(status,STATUS_INPUT_ERROR,'SUPG requires run-f1 with nonlinear ZNO_F1_NT/N0.');return
      end if
    end if
    if (config%schema_version /= NATIVE_SCHEMA_VERSION) then
      call set_status(status, STATUS_INPUT_ERROR, 'Unsupported case schema_version; expected version 2.')
      return
    end if
    if (len_trim(config%case_id) == 0 .or. trim(config%case_id) == 'UNSET') then
      call set_status(status, STATUS_INPUT_ERROR, 'case_id is required.')
      return
    end if
    select case (trim(config%model_level))
    case (MODEL_M0, MODEL_M1, MODEL_M2)
      continue
    case default
      call set_status(status, STATUS_INPUT_ERROR, 'model_level is not one of the frozen M0-M2 levels.')
      return
    end select
    if (len_trim(config%mesh_file) == 0) then
      call set_status(status, STATUS_INPUT_ERROR, 'mesh_file is required.')
      return
    end if
    if (trim(config%model_level) == MODEL_M2) then
      if (trim(config%problem_type) /= PROBLEM_REFERENCE_BAR .and. &
          trim(config%problem_type) /= PROBLEM_PLAIN_PLATE .and. &
          trim(config%problem_type) /= PROBLEM_CENTER_HOLE_PLATE) then
        call set_status(status,STATUS_INPUT_ERROR, &
          'M2 problem_type is not an explicitly supported voltage-driven profile.')
        return
      end if
      if (trim(config%material_profile) /= MATERIAL_TPS_G08_REFERENCE .and. &
          trim(config%material_profile) /= 'ZNO_F1_LQ' .and. &
          trim(config%material_profile) /= 'ZNO_F1_N0' .and. &
          trim(config%material_profile) /= 'ZNO_F1_NT') then
        call set_status(status,STATUS_INPUT_ERROR, &
          'Unsupported material profile; F1 candidate profiles require the run-f1 entry.')
        return
      end if
      if (.not. ieee_is_finite(config%target_voltage) .or. &
          .not. ieee_is_finite(config%hole_radius) .or. config%hole_radius < 0.0_dp .or. &
          config%target_voltage_factor < 0.0_dp .or. config%target_voltage_factor > 1.0_dp .or. &
          config%target_joule_factor < 0.0_dp .or. config%target_joule_factor > 1.0_dp .or. &
          any(config%field_scales <= 0.0_dp) .or. &
          .not. all(ieee_is_finite(config%field_scales))) then
        call set_status(status,STATUS_INPUT_ERROR,'M2 physical targets or field scales are invalid.')
        return
      end if
      if (len_trim(config%restart_nodes_file) > 0) then
        if (config%restart_voltage_factor < 0.0_dp .or. &
            config%restart_voltage_factor > config%target_voltage_factor .or. &
            config%restart_joule_factor < 0.0_dp .or. &
            config%restart_joule_factor > config%target_joule_factor) then
          call set_status(status,STATUS_INPUT_ERROR,'M2 restart load factors are invalid.')
          return
        end if
      else if (config%restart_voltage_factor /= 0.0_dp .or. &
               config%restart_joule_factor /= 0.0_dp) then
        call set_status(status,STATUS_INPUT_ERROR, &
          'M2 restart load factors require restart_nodes_file.')
        return
      end if
      if (trim(config%problem_type) == PROBLEM_CENTER_HOLE_PLATE .and. &
          config%hole_radius <= 0.0_dp) then
        call set_status(status,STATUS_INPUT_ERROR, &
          'CENTER_HOLE_PLATE requires a positive hole_radius.')
        return
      end if
      if (trim(config%problem_type) /= PROBLEM_CENTER_HOLE_PLATE .and. &
          config%hole_radius /= 0.0_dp) then
        call set_status(status,STATUS_INPUT_ERROR, &
          'Only CENTER_HOLE_PLATE may specify a nonzero hole_radius.')
        return
      end if
      if (config%nonlinear_options%max_newton < 1 .or. &
          config%nonlinear_options%residual_tolerance <= 0.0_dp .or. &
          config%nonlinear_options%update_tolerance <= 0.0_dp .or. &
          config%nonlinear_options%armijo_c1 <= 0.0_dp .or. &
          config%nonlinear_options%armijo_c1 >= 0.5_dp .or. &
          config%nonlinear_options%backtrack_factor <= 0.0_dp .or. &
          config%nonlinear_options%backtrack_factor >= 1.0_dp .or. &
          config%nonlinear_options%minimum_step <= 0.0_dp .or. &
          config%nonlinear_options%maximum_log_density_update <= 0.0_dp) then
        call set_status(status,STATUS_INPUT_ERROR,'M2 Newton options are invalid.')
        return
      end if
      if(config%nonlinear_options%trust_region) then
        if (.not.all(ieee_is_finite([config%nonlinear_options%trust_initial_radius, &
            config%nonlinear_options%trust_maximum_radius,config%nonlinear_options%trust_minimum_radius])).or. &
            config%nonlinear_options%trust_minimum_radius<=0.0_dp.or. &
            config%nonlinear_options%trust_initial_radius<config%nonlinear_options%trust_minimum_radius.or. &
            config%nonlinear_options%trust_maximum_radius<config%nonlinear_options%trust_initial_radius) then
          call set_status(status,STATUS_INPUT_ERROR,'M2 Newton trust region radii are invalid.')
          return
        end if
      end if
      if (config%continuation_options%initial_step <= 0.0_dp .or. &
          config%continuation_options%minimum_step <= 0.0_dp .or. &
          config%continuation_options%maximum_step < config%continuation_options%minimum_step .or. &
          config%continuation_options%initial_step > config%continuation_options%maximum_step .or. &
          config%continuation_options%cutback_factor <= 0.0_dp .or. &
          config%continuation_options%cutback_factor >= 1.0_dp .or. &
          config%continuation_options%growth_factor <= 1.0_dp .or. &
          config%continuation_options%maximum_attempts < 1 .or. &
          config%continuation_options%maximum_newton_per_attempt < 1 .or. &
          config%continuation_options%fast_iteration_threshold < 0 .or. &
          config%continuation_options%slow_iteration_threshold <= &
            config%continuation_options%fast_iteration_threshold) then
        call set_status(status,STATUS_INPUT_ERROR,'M2 continuation options are invalid.')
      end if
    end if
  end subroutine validate_case_config

end module input_schema_m
