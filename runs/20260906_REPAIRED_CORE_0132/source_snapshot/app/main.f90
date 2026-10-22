program TPSolverV2
  use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
  use kinds_m, only : dp
  use status_m, only : SolverStatusType
  use build_info_m, only : print_build_info
  use input_schema_m, only : CaseConfigType
  use native_input_m, only : read_native_case
  use mesh_types_m, only : MeshType
  use mesh_checks_m, only : validate_mesh
  use verification_driver_m, only : run_verification
  use production_runner_m, only : run_production_case
  implicit none

  character(len=512) :: command, argument, result_argument, verification_case_1, &
    verification_case_2,verification_case_3,verification_case_4
  type(SolverStatusType) :: status
  integer :: argument_count

  argument_count = command_argument_count()
  if (argument_count < 1) then
    call print_usage(error_unit)
    stop 2
  end if
  call get_command_argument(1, command)

  select case (trim(command))
  case ('version')
    call print_build_info(output_unit)
  case ('inspect')
    if (argument_count < 2) then
      call print_usage(error_unit)
      stop 2
    end if
    call get_command_argument(2, argument)
    call inspect_case(trim(argument), status)
    call finish_with_status(status)
  case ('verify')
    if (argument_count < 2) then
      call print_usage(error_unit)
      stop 2
    end if
    call get_command_argument(2, argument)
    if (trim(argument) == 'G09' .or. trim(argument) == 'g09' .or. &
        trim(argument) == 'G11' .or. trim(argument) == 'g11') then
      if (argument_count /= 4) then
        call print_usage(error_unit)
        stop 2
      end if
      call get_command_argument(3,verification_case_1)
      call get_command_argument(4,verification_case_2)
      call run_verification(trim(argument),status,trim(verification_case_1), &
        trim(verification_case_2))
    else if (trim(argument) == 'G10' .or. trim(argument) == 'g10') then
      if (argument_count /= 5) then
        call print_usage(error_unit)
        stop 2
      end if
      call get_command_argument(3,verification_case_1)
      call get_command_argument(4,verification_case_2)
      call get_command_argument(5,verification_case_3)
      call run_verification(trim(argument),status,trim(verification_case_1), &
        trim(verification_case_2),trim(verification_case_3))
    else if (trim(argument) == 'G12' .or. trim(argument) == 'g12') then
      if (argument_count /= 6) then
        call print_usage(error_unit)
        stop 2
      end if
      call get_command_argument(3,verification_case_1)
      call get_command_argument(4,verification_case_2)
      call get_command_argument(5,verification_case_3)
      call get_command_argument(6,verification_case_4)
      call run_verification(trim(argument),status,trim(verification_case_1), &
        trim(verification_case_2),trim(verification_case_3),trim(verification_case_4))
    else
      call run_verification(trim(argument), status)
    end if
    call finish_with_status(status)
  case ('run')
    if (argument_count /= 3) then
      call print_usage(error_unit)
      stop 2
    end if
    call get_command_argument(2,argument)
    call get_command_argument(3,result_argument)
    call run_production_case(trim(argument),trim(result_argument),status)
    call finish_with_status(status)
  case ('jacobian')
    write(error_unit,'(a)') &
      'jacobian is intentionally unavailable; use verify G06 for the frozen Jacobian audit.'
    stop 3
  case default
    write(error_unit,'(a)') 'Unknown command: ' // trim(command)
    call print_usage(error_unit)
    stop 2
  end select

contains

  subroutine inspect_case(case_path, local_status)
    character(len=*), intent(in) :: case_path
    type(SolverStatusType), intent(out) :: local_status
    type(CaseConfigType) :: config
    type(MeshType) :: mesh
    real(dp) :: minimum_determinant, total_volume

    call read_native_case(case_path, config, mesh, local_status)
    if (.not. local_status%is_ok()) return
    call validate_mesh(mesh, minimum_determinant, total_volume, local_status)
    if (.not. local_status%is_ok()) return

    write(output_unit,'(a)') 'inspect_result=PASS'
    write(output_unit,'(a)') 'case_id=' // trim(config%case_id)
    write(output_unit,'(a,i0)') 'schema_version=', config%schema_version
    write(output_unit,'(a)') 'model_level=' // trim(config%model_level)
    write(output_unit,'(a,i0)') 'node_count=', mesh%number_of_nodes()
    write(output_unit,'(a,i0)') 'element_count=', mesh%number_of_elements()
    write(output_unit,'(a,es24.16e3)') 'minimum_jacobian_determinant=', minimum_determinant
    write(output_unit,'(a,es24.16e3)') 'integrated_volume=', total_volume
  end subroutine inspect_case

  subroutine finish_with_status(local_status)
    type(SolverStatusType), intent(in) :: local_status
    if (.not. local_status%is_ok()) then
      write(error_unit,'(a,i0)') 'status_code=', local_status%code
      write(error_unit,'(a)') 'status_message=' // trim(local_status%message)
      stop 1
    end if
  end subroutine finish_with_status

  subroutine print_usage(unit_id)
    integer, intent(in) :: unit_id
    write(unit_id,'(a)') 'Usage:'
    write(unit_id,'(a)') '  TPSolverV2.exe version'
    write(unit_id,'(a)') '  TPSolverV2.exe inspect <case.nml>'
    write(unit_id,'(a)') '  TPSolverV2.exe verify V01'
    write(unit_id,'(a)') '  TPSolverV2.exe verify G02'
    write(unit_id,'(a)') '  TPSolverV2.exe verify G03'
    write(unit_id,'(a)') '  TPSolverV2.exe verify G04'
    write(unit_id,'(a)') '  TPSolverV2.exe verify L01'
    write(unit_id,'(a)') '  TPSolverV2.exe verify V02|V03|V04|V05|V06'
    write(unit_id,'(a)') '  TPSolverV2.exe verify G05'
    write(unit_id,'(a)') '  TPSolverV2.exe verify G06'
    write(unit_id,'(a)') '  TPSolverV2.exe verify G07'
    write(unit_id,'(a)') '  TPSolverV2.exe verify G08'
    write(unit_id,'(a)') '  TPSolverV2.exe verify G09 <H00-case.nml> <H02-case.nml>'
    write(unit_id,'(a)') '  TPSolverV2.exe verify G10 <L0-case.nml> <L1-case.nml> <L2-case.nml>'
    write(unit_id,'(a)') '  TPSolverV2.exe verify G11 <H00-L2-case.nml> <H02-L2-case.nml>'
    write(unit_id,'(a)') '  TPSolverV2.exe verify G12 <H00-case.nml> <R060-case.nml> <R090-case.nml> <R120-case.nml>'
    write(unit_id,'(a)') '  TPSolverV2.exe run <case.nml> <result.dat>'
  end subroutine print_usage

end program TPSolverV2
