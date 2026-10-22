module model_types_m
  use kinds_m, only : dp
  use field_layout_m, only : NDOF_NODE
  implicit none
  private

  type, public :: MaterialType
    character(len=64) :: name = 'UNSET'
    logical :: frozen_transport = .false. ! N0/LQ: freeze BOTH mobility and diffusion at T0.
    logical :: linear_carriers = .false. ! LQ relative perturbations, not log densities.
    logical :: transport_supg = .false. ! Opt-in consistent HEX20 residual stabilization.
    logical :: transport_log_weighted = .false. ! Opt-in density-normalized weak transport.
    logical :: heat_source_current_work = .false. ! Opt-in diagnostic Q=E dot (Jn+Jp); not a thermodynamic closure.
    integer :: log_weighted_quadrature_order = 3
    logical :: supg_directional_tau = .false.
    real(dp) :: supg_scale = 1.0_dp
    real(dp) :: elastic(6,6) = 0.0_dp
    real(dp) :: piezoelectric(3,6) = 0.0_dp
    real(dp) :: permittivity(3,3) = 0.0_dp
    real(dp) :: thermal_stress(6) = 0.0_dp
    real(dp) :: pyroelectric(3) = 0.0_dp
    real(dp) :: thermal_conductivity(3,3) = 0.0_dp
    real(dp) :: mobility_n0(3,3) = 0.0_dp
    real(dp) :: mobility_p0(3,3) = 0.0_dp
    real(dp) :: beta_n = 0.0_dp
    real(dp) :: beta_p = 0.0_dp
    real(dp) :: reference_temperature = 300.0_dp
    real(dp) :: intrinsic_concentration = 0.0_dp
    real(dp) :: donor_concentration = 0.0_dp
    real(dp) :: acceptor_concentration = 0.0_dp
    real(dp) :: reference_n = 1.0_dp
    real(dp) :: reference_p = 1.0_dp
  end type MaterialType

  type, public :: StateType
    real(dp), allocatable :: x(:)
  end type StateType

  type, public :: LoadStateType
    real(dp) :: voltage_factor = 0.0_dp
    real(dp) :: joule_factor = 0.0_dp
  end type LoadStateType

  type, public :: NonlinearOptionsType
    logical :: verbose = .false.
    logical :: global_log_step_limit = .false.
    logical :: trust_region = .false. ! Optional dogleg globalization; physical residual unchanged.
    real(dp) :: trust_initial_radius = 0.25_dp ! RMS in existing variable scales.
    real(dp) :: trust_maximum_radius = 4.0_dp
    real(dp) :: trust_minimum_radius = 1.0e-12_dp
    integer :: max_newton = 40
    integer :: stagnation_window = 0
    real(dp) :: residual_tolerance = 1.0e-8_dp
    real(dp) :: update_tolerance = 1.0e-8_dp
    real(dp) :: armijo_c1 = 1.0e-4_dp
    real(dp) :: backtrack_factor = 0.5_dp
    real(dp) :: minimum_step = 1.0e-8_dp
    real(dp) :: maximum_log_density_update = 1.0_dp
  end type NonlinearOptionsType

  type, public :: ContinuationOptionsType
    logical :: full_joule_during_voltage = .false.
    logical :: secant_predictor = .false.
    logical :: arclength_enabled = .false. ! Opt-in pseudo-arclength continuation for fold bifurcations.
    real(dp) :: voltage_checkpoint_spacing = 0.0_dp
    real(dp) :: heat_balance_reference_power = 0.0_dp ! Optional physical scale when QJ is exactly off.
    real(dp) :: initial_step = 0.25_dp
    real(dp) :: minimum_step = 0.03125_dp
    real(dp) :: maximum_step = 0.5_dp
    real(dp) :: cutback_factor = 0.5_dp
    real(dp) :: growth_factor = 2.0_dp
    integer :: maximum_attempts = 100
    integer :: maximum_newton_per_attempt = 20
    integer :: fast_iteration_threshold = 4
    integer :: slow_iteration_threshold = 10
    real(dp) :: arclength_psi_voltage = 1.0_dp ! Arc-length weight: load-factor contribution.
    real(dp) :: arclength_psi_state = 1.0_dp   ! Arc-length weight: state-vector contribution.
    real(dp) :: arclength_fold_tolerance = 1.0e-4_dp ! Scaled RMS width of the tangent-sign bracket.
    integer :: arclength_maximum_refinements = 24
  end type ContinuationOptionsType

  type, public :: DiagnosticsType
    real(dp) :: residual_by_field(NDOF_NODE) = 0.0_dp
    real(dp) :: update_by_field(NDOF_NODE) = 0.0_dp
    real(dp) :: min_temperature = huge(1.0_dp)
    real(dp) :: max_temperature = -huge(1.0_dp)
    real(dp) :: min_n = huge(1.0_dp)
    real(dp) :: max_n = -huge(1.0_dp)
    real(dp) :: min_p = huge(1.0_dp)
    real(dp) :: max_p = -huge(1.0_dp)
    real(dp) :: charge_imbalance = 0.0_dp
    real(dp) :: current_imbalance = 0.0_dp
    real(dp) :: heat_imbalance = 0.0_dp
    real(dp) :: max_peclet_n = 0.0_dp
    real(dp) :: max_peclet_p = 0.0_dp
  end type DiagnosticsType

end module model_types_m
