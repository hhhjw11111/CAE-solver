module a1_joint_check_m
  use kinds_m, only : dp
  use constants_m, only : ELEMENTARY_CHARGE,BOLTZMANN_CONSTANT
  use status_m, only : SolverStatusType,clear_status,set_status,STATUS_VERIFICATION_FAILED
  use field_layout_m, only : F_U1,F_PHI,F_ETA_N,F_ETA_P,F_THETA,global_dof
  use mesh_types_m, only : MeshType
  use pardiso_solver_m, only : PardisoDiagnosticsType
  use m0_benchmark_m, only : M0CoefficientsType,M0_A1,solve_m0_benchmark, &
    compute_m0_field_error,compute_m0_max_current,exact_m0_field
  implicit none
  private
  public :: verify_a1_joint
contains
  subroutine verify_a1_joint(output_directory,status)
    character(len=*), intent(in) :: output_directory
    type(SolverStatusType), intent(out) :: status
    integer, parameter :: counts(5)=[1,2,4,8,16]
    integer, parameter :: fields(5)=[F_U1,F_PHI,F_ETA_N,F_ETA_P,F_THETA]
    character(len=8), parameter :: names(5)=[character(len=8)::'u','phi','eta_n','eta_p','theta']
    type(M0CoefficientsType) :: c
    type(MeshType) :: mesh
    type(PardisoDiagnosticsType) :: diagnostics
    real(dp), allocatable :: solution(:)
    real(dp) :: l2(5,5),h1(5,5),orders(2),jn,jp,maximum_current,maximum_eta
    integer :: level,field,unit_id,io,node,maximum_pivots
    real(dp) :: exact(5),derivative,x
    character(len=64) :: filename,key
    logical :: passed

    call clear_status(status)
    c%thermal_expansion=1.0e-5_dp
    c%diffusion_n=BOLTZMANN_CONSTANT*300.0_dp*c%mobility_n/ELEMENTARY_CHARGE
    c%diffusion_p=BOLTZMANN_CONSTANT*300.0_dp*c%mobility_p/ELEMENTARY_CHARGE
    maximum_pivots=0;maximum_current=0.0_dp;maximum_eta=0.0_dp
    write(*,'(a)') 'verification_id=A1_F1_joint_two_carrier'
    write(*,'(a)') 'model=M0_linear_perturbations_joule_off'
    write(*,'(a)') 'error_quadrature=5x5x5'
    do level=1,5
      call solve_m0_benchmark(counts(level),M0_A1,c,mesh,solution,diagnostics,status)
      if (.not. status%is_ok()) return
      maximum_pivots=max(maximum_pivots,diagnostics%perturbed_pivots)
      do field=1,5
        call compute_m0_field_error(mesh,solution,M0_A1,fields(field),c,l2(field,level),h1(field,level),status)
        if (.not. status%is_ok()) return
        write(key,'("a1_ne",i0,"_",a,"_l2")') counts(level),trim(names(field))
        call metric(trim(key),l2(field,level))
        write(key,'("a1_ne",i0,"_",a,"_h1")') counts(level),trim(names(field))
        call metric(trim(key),h1(field,level))
      end do
      call compute_m0_max_current(mesh,solution,M0_A1,c,jn,status,F_ETA_N)
      if (.not. status%is_ok()) return
      call compute_m0_max_current(mesh,solution,M0_A1,c,jp,status,F_ETA_P)
      if (.not. status%is_ok()) return
      maximum_current=max(maximum_current,jn,jp)
      do node=1,mesh%number_of_nodes()
        maximum_eta=max(maximum_eta,abs(solution(global_dof(node,F_ETA_N))), &
          abs(solution(global_dof(node,F_ETA_P))))
      end do
      if (len_trim(output_directory) > 0) then
        write(filename,'("A1_ne",i0,"_nodes.csv")') counts(level)
        open(newunit=unit_id,file=trim(output_directory)//'/'//trim(filename), &
          status='new',action='write',iostat=io)
        if (io /= 0) then
          call set_status(status,STATUS_VERIFICATION_FAILED,'Cannot create new A1 nodal output.');return
        end if
        write(unit_id,'(a)') 'x_m,y_m,z_m,u_m,phi_V,eta_n,eta_p,theta_K,'// &
          'u_exact_m,phi_exact_V,eta_n_exact,eta_p_exact,theta_exact_K'
        do node=1,mesh%number_of_nodes()
          x=mesh%coordinates(1,node)
          do field=1,5
            call exact_m0_field(M0_A1,fields(field),c,x,exact(field),derivative,status)
            if (.not. status%is_ok()) then
              close(unit_id);return
            end if
          end do
          write(unit_id,'(*(es24.16e3,:,","))') mesh%coordinates(:,node), &
            (solution(global_dof(node,fields(field))),field=1,5),exact
        end do
        close(unit_id)
      end if
    end do
    passed=maximum_pivots == 0 .and. maximum_current < 1.0e-10_dp .and. maximum_eta <= 0.05_dp
    do field=2,4
      orders=[log(l2(field,4)/l2(field,5)),log(h1(field,4)/h1(field,5))]/log(2.0_dp)
      call metric('a1_'//trim(names(field))//'_l2_order',orders(1))
      call metric('a1_'//trim(names(field))//'_h1_order',orders(2))
      passed=passed .and. orders(1)>2.7_dp .and. orders(2)>1.8_dp .and. &
        l2(field,5)<2.0e-5_dp .and. h1(field,5)<2.0e-3_dp
    end do
    passed=passed .and. l2(1,5)<2.0e-5_dp .and. h1(1,5)<2.0e-3_dp .and. &
      maxval(l2(5,:))<1.0e-10_dp .and. maxval(h1(5,:))<1.0e-10_dp
    call metric('a1_maximum_normalized_current',maximum_current)
    call metric('a1_maximum_eta',maximum_eta)
    write(*,'(a,i0)') 'a1_maximum_perturbed_pivots=',maximum_pivots
    if (passed) then
      write(*,'(a)') 'verification_result=PASS'
    else
      write(*,'(a)') 'verification_result=FAIL'
      call set_status(status,STATUS_VERIFICATION_FAILED,'A1 frozen joint analytic criteria not satisfied.')
    end if
  end subroutine verify_a1_joint

  subroutine metric(name,value)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: value
    write(*,'(a,es24.16e3)') trim(name)//'=',value
  end subroutine metric
end module a1_joint_check_m
