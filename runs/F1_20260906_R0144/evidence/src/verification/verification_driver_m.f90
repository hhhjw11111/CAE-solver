module verification_driver_m
  use kinds_m, only : dp
  use constants_m, only : VERIFICATION_TOL
  use status_m, only : SolverStatusType, clear_status, set_status, STATUS_VERIFICATION_FAILED
  use field_layout_m, only : NDOF_NODE
  use mesh_types_m, only : MeshType, HEX20_NNODE
  use hex20_shape_m, only : hex20_shape, hex20_geometry, hex20_natural_coordinates
  use hex20_quadrature_m, only : HEX20_NGAUSS, gauss_rule_3x3x3
  use dof_map_m, only : HEX20_NDOF, build_element_dof_map
  use mesh_checks_m, only : validate_mesh
  use material_tangent_check_m, only : verify_material_point_g02
  use element_tangent_check_m, only : verify_element_g03
  use global_assembly_check_m, only : verify_global_g04
  use global_block_jacobian_check_m, only : verify_global_block_jacobian_g06
  use linear_solver_check_m, only : verify_linear_solver_l01
  use analytic_benchmarks_m, only : run_m0_verification
  use newton_solver_check_m, only : verify_newton_solver_g07
  use continuation_solver_check_m, only : verify_continuation_solver_g08
  use center_hole_check_m, only : verify_center_hole_g09
  use mesh_convergence_check_m, only : verify_mesh_convergence_g10
  use matched_control_check_m, only : verify_matched_control_g11
  use radius_matrix_check_m, only : verify_radius_matrix_g12,verify_radius_controls_g12c
  use production_joule_check_m, only : verify_production_joule_a1j
  use a1_joint_check_m, only : verify_a1_joint
  use a2_hole_check_m, only : verify_a2_hole
  implicit none
  private

  public :: run_verification

contains

  subroutine run_verification(verification_id, status, case_path_1, case_path_2,case_path_3,case_path_4)
    character(len=*), intent(in) :: verification_id
    type(SolverStatusType), intent(out) :: status
    character(len=*), intent(in), optional :: case_path_1,case_path_2,case_path_3,case_path_4
    character(len=512) :: g10_paths(3)

    select case (uppercase(trim(verification_id)))
    case ('V01')
      call verify_v01(status)
    case ('G02')
      call verify_material_point_g02(status)
    case ('F1M')
      call verify_material_point_g02(status,'NT')
      if (.not. status%is_ok()) return
      call verify_material_point_g02(status,'N0')
      if (.not. status%is_ok()) return
      call verify_material_point_g02(status,'LQ')
    case ('G03')
      call verify_element_g03(status)
    case ('G04')
      call verify_global_g04(status)
    case ('G06')
      call verify_global_block_jacobian_g06(status)
    case ('G07')
      call verify_newton_solver_g07(status)
    case ('G08')
      call verify_continuation_solver_g08(status)
    case ('A1J')
      call verify_production_joule_a1j(status)
    case ('A1')
      if (present(case_path_1)) then
        call verify_a1_joint(trim(case_path_1),status)
      else
        call verify_a1_joint('',status)
      end if
    case ('A2')
      if (present(case_path_1)) then
        call verify_a2_hole(trim(case_path_1),status)
      else
        call verify_a2_hole('',status)
      end if
    case ('G12C')
      call verify_radius_controls_g12c(status)
    case ('G09')
      if (.not. present(case_path_1) .or. .not. present(case_path_2)) then
        call set_status(status,STATUS_VERIFICATION_FAILED, &
          'G09 requires the frozen H00 and H02 case paths.')
      else
        call verify_center_hole_g09(case_path_1,case_path_2,status)
      end if
    case ('G10')
      if (.not. present(case_path_1) .or. .not. present(case_path_2) .or. &
          .not. present(case_path_3)) then
        call set_status(status,STATUS_VERIFICATION_FAILED, &
          'G10 requires the frozen L0, L1, and L2 case paths.')
      else
        g10_paths=[character(len=512) :: trim(case_path_1),trim(case_path_2), &
          trim(case_path_3)]
        call verify_mesh_convergence_g10(g10_paths,status)
      end if
    case ('G11')
      if (.not. present(case_path_1) .or. .not. present(case_path_2)) then
        call set_status(status,STATUS_VERIFICATION_FAILED, &
          'G11 requires the frozen matched H00-L2 and H02-L2 case paths.')
      else
        call verify_matched_control_g11(case_path_1,case_path_2,status)
      end if
    case ('G12')
      if (.not. present(case_path_1) .or. .not. present(case_path_2) .or. &
          .not. present(case_path_3) .or. .not. present(case_path_4)) then
        call set_status(status,STATUS_VERIFICATION_FAILED, &
          'G12 requires the frozen H00, R060, R090, and R120 case paths.')
      else
        g10_paths=[character(len=512) :: trim(case_path_2),trim(case_path_3), &
          trim(case_path_4)]
        call verify_radius_matrix_g12(trim(case_path_1),g10_paths,status)
      end if
    case ('L01')
      call verify_linear_solver_l01(status)
    case ('V02','V03','V04','V05','V06','G05')
      call run_m0_verification(verification_id,status)
    case default
      call set_status(status, STATUS_VERIFICATION_FAILED, &
        'Unknown verification identifier: ' // trim(verification_id))
    end select
  end subroutine run_verification

  subroutine verify_v01(status)
    type(SolverStatusType), intent(out) :: status
    type(MeshType) :: mesh
    type(SolverStatusType) :: local_status
    real(dp) :: natural_nodes(3,HEX20_NNODE), physical_nodes(3,HEX20_NNODE)
    real(dp) :: points(3,HEX20_NGAUSS), weights(HEX20_NGAUSS)
    real(dp) :: shape(HEX20_NNODE), derivative_natural(3,HEX20_NNODE)
    real(dp) :: derivative_global(3,HEX20_NNODE), determinant
    real(dp) :: expected_shape(HEX20_NNODE), mapped_point(3), expected_point(3)
    real(dp) :: reproduced_gradient(3), expected_gradient(3)
    real(dp) :: max_kronecker_error, max_partition_error, max_derivative_sum
    real(dp) :: max_mapping_error, max_gradient_error, max_determinant_error
    real(dp) :: quadrature_weight_error, minimum_determinant, total_volume, volume_error
    integer :: point, node, coordinate_index
    integer :: connectivity(HEX20_NNODE), element_dof(HEX20_NDOF)
    logical :: passed

    call clear_status(status)
    call hex20_natural_coordinates(natural_nodes)
    physical_nodes = 0.5_dp * (natural_nodes + 1.0_dp)
    call gauss_rule_3x3x3(points, weights)

    max_kronecker_error = 0.0_dp
    do node = 1, HEX20_NNODE
      call hex20_shape(natural_nodes(1,node), natural_nodes(2,node), &
        natural_nodes(3,node), shape, derivative_natural)
      expected_shape = 0.0_dp
      expected_shape(node) = 1.0_dp
      max_kronecker_error = max(max_kronecker_error, maxval(abs(shape-expected_shape)))
    end do

    max_partition_error = 0.0_dp
    max_derivative_sum = 0.0_dp
    max_mapping_error = 0.0_dp
    max_gradient_error = 0.0_dp
    max_determinant_error = 0.0_dp
    do point = 1, HEX20_NGAUSS
      call hex20_geometry(physical_nodes, points(1,point), points(2,point), &
        points(3,point), shape, derivative_global, determinant, local_status)
      if (.not. local_status%is_ok()) then
        call set_status(status, STATUS_VERIFICATION_FAILED, trim(local_status%message))
        return
      end if
      call hex20_shape(points(1,point), points(2,point), points(3,point), &
        shape, derivative_natural)
      max_partition_error = max(max_partition_error, abs(sum(shape)-1.0_dp))
      max_derivative_sum = max(max_derivative_sum, maxval(abs(sum(derivative_natural,dim=2))))
      mapped_point = matmul(physical_nodes, shape)
      expected_point = 0.5_dp * (points(:,point) + 1.0_dp)
      max_mapping_error = max(max_mapping_error, maxval(abs(mapped_point-expected_point)))
      max_determinant_error = max(max_determinant_error, abs(determinant-0.125_dp))
      do coordinate_index = 1, 3
        reproduced_gradient = matmul(derivative_global, physical_nodes(coordinate_index,:))
        expected_gradient = 0.0_dp
        expected_gradient(coordinate_index) = 1.0_dp
        max_gradient_error = max(max_gradient_error, &
          maxval(abs(reproduced_gradient-expected_gradient)))
      end do
    end do
    quadrature_weight_error = abs(sum(weights)-8.0_dp)

    allocate(mesh%coordinates(3,HEX20_NNODE))
    allocate(mesh%connectivity(HEX20_NNODE,1))
    allocate(mesh%material_id(1))
    mesh%coordinates = physical_nodes
    connectivity = [(node,node=1,HEX20_NNODE)]
    mesh%connectivity(:,1) = connectivity
    mesh%material_id(1) = 1
    call validate_mesh(mesh, minimum_determinant, total_volume, local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status, STATUS_VERIFICATION_FAILED, trim(local_status%message))
      return
    end if
    volume_error = abs(total_volume-1.0_dp)

    call build_element_dof_map(connectivity, element_dof, local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status, STATUS_VERIFICATION_FAILED, trim(local_status%message))
      return
    end if

    passed = max_kronecker_error <= VERIFICATION_TOL .and. &
             max_partition_error <= VERIFICATION_TOL .and. &
             max_derivative_sum <= VERIFICATION_TOL .and. &
             max_mapping_error <= VERIFICATION_TOL .and. &
             max_gradient_error <= VERIFICATION_TOL .and. &
             max_determinant_error <= VERIFICATION_TOL .and. &
             quadrature_weight_error <= VERIFICATION_TOL .and. &
             volume_error <= VERIFICATION_TOL .and. &
             abs(minimum_determinant-0.125_dp) <= VERIFICATION_TOL .and. &
             element_dof(1) == 1 .and. element_dof(NDOF_NODE) == NDOF_NODE .and. &
             element_dof(NDOF_NODE+1) == NDOF_NODE+1 .and. &
             element_dof(HEX20_NDOF) == HEX20_NDOF

    write(*,'(a)') 'verification_id=V01'
    call print_metric('max_kronecker_error', max_kronecker_error)
    call print_metric('max_partition_error', max_partition_error)
    call print_metric('max_natural_derivative_sum', max_derivative_sum)
    call print_metric('max_linear_mapping_error', max_mapping_error)
    call print_metric('max_linear_gradient_error', max_gradient_error)
    call print_metric('max_jacobian_determinant_error', max_determinant_error)
    call print_metric('quadrature_weight_error', quadrature_weight_error)
    call print_metric('unit_cube_volume_error', volume_error)
    call print_metric('minimum_jacobian_determinant', minimum_determinant)
    write(*,'(a,i0)') 'hex20_local_dof_count=', HEX20_NDOF
    if (passed) then
      write(*,'(a)') 'verification_result=PASS'
    else
      write(*,'(a)') 'verification_result=FAIL'
      call set_status(status, STATUS_VERIFICATION_FAILED, 'V01 metric exceeded the frozen tolerance.')
    end if
  end subroutine verify_v01

  subroutine print_metric(name, value)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: value
    write(*,'(a,es24.16e3)') trim(name) // '=', value
  end subroutine print_metric

  pure function uppercase(text) result(converted)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: converted
    integer :: index, code
    do index = 1, len(text)
      code = iachar(text(index:index))
      if (code >= iachar('a') .and. code <= iachar('z')) then
        converted(index:index) = achar(code-32)
      else
        converted(index:index) = text(index:index)
      end if
    end do
  end function uppercase

end module verification_driver_m
