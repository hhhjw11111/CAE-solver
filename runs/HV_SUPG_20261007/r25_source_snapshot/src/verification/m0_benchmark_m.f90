module m0_benchmark_m
  use kinds_m, only : dp
  use constants_m, only : ELEMENTARY_CHARGE
  use status_m, only : SolverStatusType, clear_status, set_status, &
    STATUS_INVALID_ARGUMENT, STATUS_VERIFICATION_FAILED
  use field_layout_m, only : NDOF_NODE, F_U1, F_PHI, F_ETA_N, F_ETA_P, F_THETA, &
    global_dof, local_dof
  use mesh_types_m, only : MeshType, HEX20_NNODE
  use hex20_shape_m, only : hex20_natural_coordinates, hex20_shape, hex20_geometry
  use hex20_quadrature_m, only : HEX20_NGAUSS, gauss_rule_3x3x3
  use sparse_pattern_m, only : SparseAssemblyContextType, initialize_sparse_context, &
    find_csr_position
  use scaling_m, only : build_variable_scaling, equilibrate_csr_system, unscale_solution
  use pardiso_solver_m, only : PardisoDiagnosticsType, solve_csr_pardiso
  implicit none
  private

  integer, parameter, public :: M0_V02 = 2
  integer, parameter, public :: M0_V03 = 3
  integer, parameter, public :: M0_V04 = 4
  integer, parameter, public :: M0_V05_LINEAR = 51
  integer, parameter, public :: M0_V05_UNIFORM = 52
  integer, parameter, public :: M0_V06 = 6
  integer, parameter, public :: M0_A1 = 7

  type, public :: M0CoefficientsType
    real(dp) :: length = 1.0e-6_dp
    real(dp) :: width = 1.0e-6_dp
    real(dp) :: elastic_c = 1.0e11_dp
    real(dp) :: piezo_e = 1.0_dp
    real(dp) :: permittivity = 1.0e-8_dp
    real(dp) :: reference_n = 1.0e21_dp
    real(dp) :: reference_p = 8.0e20_dp
    real(dp) :: mobility_n = 5.0e-2_dp
    real(dp) :: diffusion_n = 1.2926e-3_dp
    real(dp) :: mobility_p = 2.0e-2_dp
    real(dp) :: diffusion_p = 5.1704e-4_dp
    real(dp) :: thermal_conductivity = 1.0e2_dp
    real(dp) :: thermal_expansion = 1.0e-7_dp
    real(dp) :: pyroelectric = 1.0e-8_dp
    real(dp) :: voltage = 1.0e-3_dp
    real(dp) :: delta_temperature = 1.0_dp
    real(dp) :: uniform_heat_source = 8.0e14_dp
  end type M0CoefficientsType

  public :: solve_m0_benchmark, compute_m0_field_error, compute_m0_max_current
  public :: compute_m0_uniform_heat_balance, exact_m0_field

contains

  subroutine solve_m0_benchmark(element_count,benchmark,coefficients,mesh,solution, &
                                diagnostics,status)
    integer, intent(in) :: element_count,benchmark
    type(M0CoefficientsType), intent(in) :: coefficients
    type(MeshType), intent(out) :: mesh
    real(dp), allocatable, intent(out) :: solution(:)
    type(PardisoDiagnosticsType), intent(out) :: diagnostics
    type(SolverStatusType), intent(out) :: status
    type(SparseAssemblyContextType) :: context
    type(SolverStatusType) :: local_status
    real(dp), allocatable :: values(:),rhs(:),variable_scales(:),residual_scales(:)
    real(dp), allocatable :: scaled_values(:),scaled_rhs(:),scaled_solution(:)
    real(dp) :: field_scales(NDOF_NODE)

    call clear_status(status)
    if (element_count < 1 .or. .not. valid_benchmark(benchmark)) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'Invalid M0 benchmark request.')
      return
    end if
    call validate_coefficients(coefficients,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,local_status%code,trim(local_status%message))
      return
    end if
    call build_bar_mesh(element_count,coefficients%length,coefficients%width,mesh)
    call initialize_sparse_context(mesh,context,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message))
      return
    end if
    call assemble_m0_system(mesh,context,benchmark,coefficients,values,rhs,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message))
      return
    end if
    call apply_m0_boundary_conditions(mesh,context,benchmark,coefficients,values,rhs, &
      local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message))
      return
    end if

    call benchmark_field_scales(benchmark,coefficients,field_scales)
    call build_variable_scaling(mesh%number_of_nodes(),field_scales,variable_scales,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message))
      return
    end if
    call equilibrate_csr_system(context,values,rhs,variable_scales,scaled_values, &
      scaled_rhs,residual_scales,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message))
      return
    end if
    allocate(scaled_solution(context%number_of_equations),solution(context%number_of_equations))
    call solve_csr_pardiso(context%number_of_equations,context%ia,context%ja, &
      scaled_values,scaled_rhs,scaled_solution,diagnostics,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message))
      return
    end if
    call unscale_solution(variable_scales,scaled_solution,solution,local_status)
    if (.not. local_status%is_ok()) then
      call set_status(status,STATUS_VERIFICATION_FAILED,trim(local_status%message))
    end if
  end subroutine solve_m0_benchmark

  subroutine assemble_m0_system(mesh,context,benchmark,coefficients,values,rhs,status)
    type(MeshType), intent(in) :: mesh
    type(SparseAssemblyContextType), intent(in) :: context
    integer, intent(in) :: benchmark
    type(M0CoefficientsType), intent(in) :: coefficients
    real(dp), allocatable, intent(out) :: values(:),rhs(:)
    type(SolverStatusType), intent(out) :: status
    real(dp) :: points(3,HEX20_NGAUSS),weights(HEX20_NGAUSS)
    real(dp) :: shape(HEX20_NNODE),gradient(3,HEX20_NNODE),determinant,weight
    real(dp) :: gradient_product,gradient_shape,shape_product,charge_coefficient
    real(dp) :: coordinates(3,HEX20_NNODE)
    integer :: element,point,row_node,column_node,node,field,row
    type(SolverStatusType) :: local_status

    call clear_status(status)
    allocate(values(context%number_of_nonzeros),rhs(context%number_of_equations))
    values=0.0_dp
    rhs=0.0_dp
    call gauss_rule_3x3x3(points,weights)
    do element=1,mesh%number_of_elements()
      coordinates=mesh%coordinates(:,mesh%connectivity(:,element))
      do point=1,HEX20_NGAUSS
        call hex20_geometry(coordinates,points(1,point),points(2,point),points(3,point), &
          shape,gradient,determinant,local_status)
        if (.not. local_status%is_ok()) then
          call set_status(status,STATUS_VERIFICATION_FAILED,'M0 element geometry failed.')
          return
        end if
        weight=weights(point)*determinant
        do row_node=1,HEX20_NNODE
          if (benchmark == M0_V05_UNIFORM) then
            row=context%element_dof(local_dof(row_node,F_THETA),element)
            rhs(row)=rhs(row)+weight*shape(row_node)*coefficients%uniform_heat_source
          end if
          do column_node=1,HEX20_NNODE
            gradient_product=weight*dot_product(gradient(:,row_node),gradient(:,column_node))
            gradient_shape=weight*gradient(1,row_node)*shape(column_node)
            shape_product=weight*shape(row_node)*shape(column_node)
            select case(benchmark)
            case(M0_V02)
              charge_coefficient=ELEMENTARY_CHARGE*coefficients%reference_n* &
                coefficients%mobility_n/coefficients%diffusion_n
              call add_entry(context,element,row_node,F_U1,column_node,F_U1, &
                coefficients%elastic_c*gradient_product,values)
              call add_entry(context,element,row_node,F_U1,column_node,F_PHI, &
                coefficients%piezo_e*gradient_product,values)
              call add_entry(context,element,row_node,F_PHI,column_node,F_U1, &
                coefficients%piezo_e*gradient_product,values)
              call add_entry(context,element,row_node,F_PHI,column_node,F_PHI, &
                -coefficients%permittivity*gradient_product- &
                charge_coefficient*shape_product,values)
            case(M0_V03)
              call add_entry(context,element,row_node,F_ETA_N,column_node,F_PHI, &
                -ELEMENTARY_CHARGE*coefficients%reference_n* &
                coefficients%mobility_n*gradient_product,values)
              call add_entry(context,element,row_node,F_ETA_N,column_node,F_ETA_N, &
                ELEMENTARY_CHARGE*coefficients%reference_n* &
                coefficients%diffusion_n*gradient_product,values)
            case(M0_V04)
              call add_entry(context,element,row_node,F_ETA_P,column_node,F_PHI, &
                -ELEMENTARY_CHARGE*coefficients%reference_p* &
                coefficients%mobility_p*gradient_product,values)
              call add_entry(context,element,row_node,F_ETA_P,column_node,F_ETA_P, &
                -ELEMENTARY_CHARGE*coefficients%reference_p* &
                coefficients%diffusion_p*gradient_product,values)
            case(M0_V05_LINEAR,M0_V05_UNIFORM)
              call add_entry(context,element,row_node,F_THETA,column_node,F_THETA, &
                coefficients%thermal_conductivity*gradient_product,values)
            case(M0_V06,M0_A1)
              call add_entry(context,element,row_node,F_U1,column_node,F_U1, &
                coefficients%elastic_c*gradient_product,values)
              call add_entry(context,element,row_node,F_U1,column_node,F_PHI, &
                coefficients%piezo_e*gradient_product,values)
              call add_entry(context,element,row_node,F_U1,column_node,F_THETA, &
                -coefficients%elastic_c*coefficients%thermal_expansion* &
                gradient_shape,values)
              call add_entry(context,element,row_node,F_PHI,column_node,F_U1, &
                coefficients%piezo_e*gradient_product,values)
              call add_entry(context,element,row_node,F_PHI,column_node,F_PHI, &
                -coefficients%permittivity*gradient_product,values)
              call add_entry(context,element,row_node,F_PHI,column_node,F_ETA_N, &
                -ELEMENTARY_CHARGE*coefficients%reference_n*shape_product,values)
              if (benchmark == M0_A1) then
                call add_entry(context,element,row_node,F_PHI,column_node,F_ETA_P, &
                  ELEMENTARY_CHARGE*coefficients%reference_p*shape_product,values)
                call add_entry(context,element,row_node,F_ETA_P,column_node,F_PHI, &
                  -ELEMENTARY_CHARGE*coefficients%reference_p* &
                  coefficients%mobility_p*gradient_product,values)
                call add_entry(context,element,row_node,F_ETA_P,column_node,F_ETA_P, &
                  -ELEMENTARY_CHARGE*coefficients%reference_p* &
                  coefficients%diffusion_p*gradient_product,values)
              end if
              call add_entry(context,element,row_node,F_PHI,column_node,F_THETA, &
                coefficients%pyroelectric*gradient_shape,values)
              call add_entry(context,element,row_node,F_ETA_N,column_node,F_PHI, &
                -ELEMENTARY_CHARGE*coefficients%reference_n* &
                coefficients%mobility_n*gradient_product,values)
              call add_entry(context,element,row_node,F_ETA_N,column_node,F_ETA_N, &
                ELEMENTARY_CHARGE*coefficients%reference_n* &
                coefficients%diffusion_n*gradient_product,values)
              call add_entry(context,element,row_node,F_THETA,column_node,F_THETA, &
                coefficients%thermal_conductivity*gradient_product,values)
            end select
          end do
        end do
      end do
    end do

    do node=1,mesh%number_of_nodes()
      do field=1,NDOF_NODE
        if (active_field(benchmark,field)) cycle
        row=global_dof(node,field)
        values(find_csr_position(context,row,row))=1.0_dp
      end do
    end do
  end subroutine assemble_m0_system

  subroutine add_entry(context,element,row_node,row_field,column_node,column_field, &
                       contribution,values)
    type(SparseAssemblyContextType), intent(in) :: context
    integer, intent(in) :: element,row_node,row_field,column_node,column_field
    real(dp), intent(in) :: contribution
    real(dp), intent(inout) :: values(:)
    integer :: global_row,offset,position

    global_row=context%element_dof(local_dof(row_node,row_field),element)
    offset=context%element_neighbor_offset(row_node,column_node,element)
    position=context%ia(global_row)+offset+column_field-1
    values(position)=values(position)+contribution
  end subroutine add_entry

  subroutine apply_m0_boundary_conditions(mesh,context,benchmark,coefficients,values,rhs,status)
    type(MeshType), intent(in) :: mesh
    type(SparseAssemblyContextType), intent(in) :: context
    integer, intent(in) :: benchmark
    type(M0CoefficientsType), intent(in) :: coefficients
    real(dp), intent(inout) :: values(:),rhs(:)
    type(SolverStatusType), intent(out) :: status
    integer :: node,left_count,right_count
    real(dp) :: tolerance,x

    call clear_status(status)
    tolerance=1.0e-10_dp*coefficients%length
    left_count=0
    right_count=0
    do node=1,mesh%number_of_nodes()
      x=mesh%coordinates(1,node)
      if (abs(x) <= tolerance) left_count=left_count+1
      if (abs(x-coefficients%length) <= tolerance) right_count=right_count+1
    end do
    if (left_count == 0 .or. right_count == 0) then
      call set_status(status,STATUS_VERIFICATION_FAILED,'M0 bar endpoint face was not found.')
      return
    end if

    select case(benchmark)
    case(M0_V02)
      do node=1,mesh%number_of_nodes()
        x=mesh%coordinates(1,node)
        if (abs(x) <= tolerance) then
          call set_total_dirichlet(context,global_dof(node,F_U1),0.0_dp,values,rhs,status)
          if (.not. status%is_ok()) return
          call set_total_dirichlet(context,global_dof(node,F_PHI),0.0_dp,values,rhs,status)
          if (.not. status%is_ok()) return
        else if (abs(x-coefficients%length) <= tolerance) then
          call set_total_dirichlet(context,global_dof(node,F_PHI), &
            coefficients%voltage,values,rhs,status)
          if (.not. status%is_ok()) return
        end if
      end do
    case(M0_V03)
      do node=1,mesh%number_of_nodes()
        x=mesh%coordinates(1,node)
        call set_total_dirichlet(context,global_dof(node,F_PHI), &
          coefficients%voltage*x/coefficients%length,values,rhs,status)
        if (.not. status%is_ok()) return
        if (abs(x) <= tolerance) then
          call set_total_dirichlet(context,global_dof(node,F_ETA_N),0.0_dp,values,rhs,status)
          if (.not. status%is_ok()) return
        end if
      end do
    case(M0_V04)
      do node=1,mesh%number_of_nodes()
        x=mesh%coordinates(1,node)
        call set_total_dirichlet(context,global_dof(node,F_PHI), &
          coefficients%voltage*x/coefficients%length,values,rhs,status)
        if (.not. status%is_ok()) return
        if (abs(x) <= tolerance) then
          call set_total_dirichlet(context,global_dof(node,F_ETA_P),0.0_dp,values,rhs,status)
          if (.not. status%is_ok()) return
        end if
      end do
    case(M0_V05_LINEAR)
      do node=1,mesh%number_of_nodes()
        x=mesh%coordinates(1,node)
        if (abs(x) <= tolerance) then
          call set_total_dirichlet(context,global_dof(node,F_THETA),0.0_dp,values,rhs,status)
        else if (abs(x-coefficients%length) <= tolerance) then
          call set_total_dirichlet(context,global_dof(node,F_THETA), &
            coefficients%delta_temperature,values,rhs,status)
        end if
        if (.not. status%is_ok()) return
      end do
    case(M0_V05_UNIFORM)
      do node=1,mesh%number_of_nodes()
        x=mesh%coordinates(1,node)
        if (abs(x) <= tolerance .or. abs(x-coefficients%length) <= tolerance) then
          call set_total_dirichlet(context,global_dof(node,F_THETA),0.0_dp,values,rhs,status)
          if (.not. status%is_ok()) return
        end if
      end do
    case(M0_V06,M0_A1)
      do node=1,mesh%number_of_nodes()
        x=mesh%coordinates(1,node)
        if (abs(x) <= tolerance) then
          call set_total_dirichlet(context,global_dof(node,F_U1),0.0_dp,values,rhs,status)
          if (.not. status%is_ok()) return
          call set_total_dirichlet(context,global_dof(node,F_PHI),0.0_dp,values,rhs,status)
          if (.not. status%is_ok()) return
          call set_total_dirichlet(context,global_dof(node,F_ETA_N),0.0_dp,values,rhs,status)
          if (.not. status%is_ok()) return
          if (benchmark == M0_A1) then
            call set_total_dirichlet(context,global_dof(node,F_ETA_P),0.0_dp,values,rhs,status)
            if (.not. status%is_ok()) return
          end if
          call set_total_dirichlet(context,global_dof(node,F_THETA),0.0_dp,values,rhs,status)
          if (.not. status%is_ok()) return
        else if (abs(x-coefficients%length) <= tolerance) then
          call set_total_dirichlet(context,global_dof(node,F_PHI),0.0_dp,values,rhs,status)
          if (.not. status%is_ok()) return
          call set_total_dirichlet(context,global_dof(node,F_THETA), &
            coefficients%delta_temperature,values,rhs,status)
          if (.not. status%is_ok()) return
        end if
      end do
    end select
  end subroutine apply_m0_boundary_conditions

  subroutine set_total_dirichlet(context,row,prescribed_value,values,rhs,status)
    type(SparseAssemblyContextType), intent(in) :: context
    integer, intent(in) :: row
    real(dp), intent(in) :: prescribed_value
    real(dp), intent(inout) :: values(:),rhs(:)
    type(SolverStatusType), intent(out) :: status
    integer :: position,diagonal

    call clear_status(status)
    if (row < 1 .or. row > context%number_of_equations) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'M0 Dirichlet row is invalid.')
      return
    end if
    do position=context%ia(row),context%ia(row+1)-1
      values(position)=0.0_dp
    end do
    diagonal=find_csr_position(context,row,row)
    if (diagonal == 0) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'M0 Dirichlet row has no diagonal.')
      return
    end if
    values(diagonal)=1.0_dp
    rhs(row)=prescribed_value
  end subroutine set_total_dirichlet

  subroutine compute_m0_field_error(mesh,solution,benchmark,field,coefficients, &
                                    relative_l2,relative_h1,status)
    type(MeshType), intent(in) :: mesh
    real(dp), intent(in) :: solution(:)
    integer, intent(in) :: benchmark,field
    type(M0CoefficientsType), intent(in) :: coefficients
    real(dp), intent(out) :: relative_l2,relative_h1
    type(SolverStatusType), intent(out) :: status
    real(dp), parameter :: gp(5)=[-0.9061798459386640_dp,-0.5384693101056831_dp, &
      0.0_dp,0.5384693101056831_dp,0.9061798459386640_dp]
    real(dp), parameter :: gw(5)=[0.2369268850561891_dp,0.4786286704993665_dp, &
      0.5688888888888889_dp,0.4786286704993665_dp,0.2369268850561891_dp]
    real(dp) :: shape(HEX20_NNODE),gradient(3,HEX20_NNODE),determinant,weight
    real(dp) :: coordinates(3,HEX20_NNODE),nodal_values(HEX20_NNODE)
    real(dp) :: numerical_value,numerical_gradient,exact_value,exact_gradient,x
    real(dp) :: error_l2,error_h1,reference_l2,reference_h1
    integer :: element,i,j,k,node
    type(SolverStatusType) :: local_status

    call clear_status(status)
    if (size(solution) /= NDOF_NODE*mesh%number_of_nodes()) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'M0 error vector size mismatch.')
      return
    end if
    error_l2=0.0_dp; error_h1=0.0_dp
    reference_l2=0.0_dp; reference_h1=0.0_dp
    do element=1,mesh%number_of_elements()
      coordinates=mesh%coordinates(:,mesh%connectivity(:,element))
      do node=1,HEX20_NNODE
        nodal_values(node)=solution(global_dof(mesh%connectivity(node,element),field))
      end do
      do i=1,5
        do j=1,5
          do k=1,5
        call hex20_geometry(coordinates,gp(i),gp(j),gp(k), &
          shape,gradient,determinant,local_status)
        if (.not. local_status%is_ok()) then
          call set_status(status,STATUS_VERIFICATION_FAILED,'M0 error integration geometry failed.')
          return
        end if
        weight=gw(i)*gw(j)*gw(k)*determinant
        x=dot_product(shape,coordinates(1,:))
        numerical_value=dot_product(shape,nodal_values)
        numerical_gradient=dot_product(gradient(1,:),nodal_values)
        call exact_m0_field(benchmark,field,coefficients,x,exact_value,exact_gradient,status)
        if (.not. status%is_ok()) return
        error_l2=error_l2+weight*(numerical_value-exact_value)**2
        error_h1=error_h1+weight*(numerical_gradient-exact_gradient)**2
        reference_l2=reference_l2+weight*exact_value**2
        reference_h1=reference_h1+weight*exact_gradient**2
          end do
        end do
      end do
    end do
    relative_l2=sqrt(error_l2)/max(sqrt(reference_l2),tiny(1.0_dp))
    relative_h1=sqrt(error_h1)/max(sqrt(reference_h1),tiny(1.0_dp))
  end subroutine compute_m0_field_error

  subroutine compute_m0_max_current(mesh,solution,benchmark,coefficients,normalized_current,status,carrier_field)
    type(MeshType), intent(in) :: mesh
    real(dp), intent(in) :: solution(:)
    integer, intent(in) :: benchmark
    type(M0CoefficientsType), intent(in) :: coefficients
    real(dp), intent(out) :: normalized_current
    type(SolverStatusType), intent(out) :: status
    integer, intent(in), optional :: carrier_field
    real(dp) :: points(3,HEX20_NGAUSS),weights(HEX20_NGAUSS)
    real(dp) :: shape(HEX20_NNODE),gradient(3,HEX20_NNODE),determinant
    real(dp) :: coordinates(3,HEX20_NNODE),phi_nodes(HEX20_NNODE)
    real(dp) :: eta_nodes(HEX20_NNODE),phi_gradient,eta_gradient,current,current_scale
    integer :: element,point,node,eta_field
    type(SolverStatusType) :: local_status

    call clear_status(status)
    if (benchmark == M0_V03 .or. benchmark == M0_V06 .or. benchmark == M0_A1) then
      eta_field=F_ETA_N
      current_scale=ELEMENTARY_CHARGE*coefficients%reference_n* &
        coefficients%mobility_n*max(abs(coefficients%voltage)/coefficients%length,1.0_dp)
    else if (benchmark == M0_V04) then
      eta_field=F_ETA_P
      current_scale=ELEMENTARY_CHARGE*coefficients%reference_p* &
        coefficients%mobility_p*max(abs(coefficients%voltage)/coefficients%length,1.0_dp)
    else
      call set_status(status,STATUS_INVALID_ARGUMENT,'Current is undefined for this M0 benchmark.')
      return
    end if
    if (present(carrier_field)) then
      if (carrier_field /= F_ETA_N .and. carrier_field /= F_ETA_P) then
        call set_status(status,STATUS_INVALID_ARGUMENT,'Invalid current carrier field.')
        return
      end if
      eta_field=carrier_field
      if (eta_field == F_ETA_P) current_scale=ELEMENTARY_CHARGE*coefficients%reference_p* &
        coefficients%mobility_p*max(abs(coefficients%voltage)/coefficients%length,1.0_dp)
    end if
    call gauss_rule_3x3x3(points,weights)
    normalized_current=0.0_dp
    do element=1,mesh%number_of_elements()
      coordinates=mesh%coordinates(:,mesh%connectivity(:,element))
      do node=1,HEX20_NNODE
        phi_nodes(node)=solution(global_dof(mesh%connectivity(node,element),F_PHI))
        eta_nodes(node)=solution(global_dof(mesh%connectivity(node,element),eta_field))
      end do
      do point=1,HEX20_NGAUSS
        call hex20_geometry(coordinates,points(1,point),points(2,point),points(3,point), &
          shape,gradient,determinant,local_status)
        if (.not. local_status%is_ok()) then
          call set_status(status,STATUS_VERIFICATION_FAILED,'M0 current geometry failed.')
          return
        end if
        phi_gradient=dot_product(gradient(1,:),phi_nodes)
        eta_gradient=dot_product(gradient(1,:),eta_nodes)
        if (eta_field == F_ETA_P) then
          current=ELEMENTARY_CHARGE*coefficients%reference_p* &
            (-coefficients%mobility_p*phi_gradient-coefficients%diffusion_p*eta_gradient)
        else
          current=ELEMENTARY_CHARGE*coefficients%reference_n* &
            (-coefficients%mobility_n*phi_gradient+coefficients%diffusion_n*eta_gradient)
        end if
        normalized_current=max(normalized_current,abs(current)/current_scale)
      end do
    end do
  end subroutine compute_m0_max_current

  subroutine compute_m0_uniform_heat_balance(mesh,solution,coefficients,relative_balance,status)
    type(MeshType), intent(in) :: mesh
    real(dp), intent(in) :: solution(:)
    type(M0CoefficientsType), intent(in) :: coefficients
    real(dp), intent(out) :: relative_balance
    type(SolverStatusType), intent(out) :: status
    real(dp), parameter :: gp(3)=[-sqrt(3.0_dp/5.0_dp),0.0_dp,sqrt(3.0_dp/5.0_dp)]
    real(dp), parameter :: gw(3)=[5.0_dp/9.0_dp,8.0_dp/9.0_dp,5.0_dp/9.0_dp]
    real(dp) :: left_flux,right_flux,source_integral

    call clear_status(status)
    call boundary_heat_flux(mesh,solution,coefficients,1,-1.0_dp,gp,gw,left_flux,status)
    if (.not. status%is_ok()) return
    call boundary_heat_flux(mesh,solution,coefficients,mesh%number_of_elements(), &
      1.0_dp,gp,gw,right_flux,status)
    if (.not. status%is_ok()) return
    source_integral=coefficients%uniform_heat_source*coefficients%length* &
      coefficients%width**2
    relative_balance=abs(left_flux+right_flux-source_integral)/abs(source_integral)
  end subroutine compute_m0_uniform_heat_balance

  subroutine boundary_heat_flux(mesh,solution,coefficients,element,xi,gp,gw,flux,status)
    type(MeshType), intent(in) :: mesh
    real(dp), intent(in) :: solution(:)
    type(M0CoefficientsType), intent(in) :: coefficients
    integer, intent(in) :: element
    real(dp), intent(in) :: xi,gp(3),gw(3)
    real(dp), intent(out) :: flux
    type(SolverStatusType), intent(out) :: status
    real(dp) :: coordinates(3,HEX20_NNODE),theta_nodes(HEX20_NNODE)
    real(dp) :: shape(HEX20_NNODE),derivative_natural(3,HEX20_NNODE)
    real(dp) :: gradient(3,HEX20_NNODE),determinant,tangent_eta(3),tangent_zeta(3)
    real(dp) :: area_jacobian,gradient_x,outward_normal
    integer :: i,j,node
    type(SolverStatusType) :: local_status

    call clear_status(status)
    coordinates=mesh%coordinates(:,mesh%connectivity(:,element))
    do node=1,HEX20_NNODE
      theta_nodes(node)=solution(global_dof(mesh%connectivity(node,element),F_THETA))
    end do
    outward_normal=merge(-1.0_dp,1.0_dp,xi < 0.0_dp)
    flux=0.0_dp
    do i=1,3
      do j=1,3
        call hex20_shape(xi,gp(i),gp(j),shape,derivative_natural)
        call hex20_geometry(coordinates,xi,gp(i),gp(j),shape,gradient,determinant,local_status)
        if (.not. local_status%is_ok()) then
          call set_status(status,STATUS_VERIFICATION_FAILED,'M0 boundary heat geometry failed.')
          return
        end if
        tangent_eta=matmul(coordinates,derivative_natural(2,:))
        tangent_zeta=matmul(coordinates,derivative_natural(3,:))
        area_jacobian=vector_norm(cross_product(tangent_eta,tangent_zeta))
        gradient_x=dot_product(gradient(1,:),theta_nodes)
        flux=flux+gw(i)*gw(j)*area_jacobian* &
          (-coefficients%thermal_conductivity*gradient_x*outward_normal)
      end do
    end do
  end subroutine boundary_heat_flux

  subroutine exact_m0_field(benchmark,field,coefficients,x,value,derivative,status)
    integer, intent(in) :: benchmark,field
    type(M0CoefficientsType), intent(in) :: coefficients
    real(dp), intent(in) :: x
    real(dp), intent(out) :: value,derivative
    type(SolverStatusType), intent(out) :: status
    real(dp) :: bar_kappa,gamma_squared,gamma,phi,phi_derivative,s,ratio

    call clear_status(status)
    value=0.0_dp
    derivative=0.0_dp
    bar_kappa=coefficients%permittivity+coefficients%piezo_e**2/coefficients%elastic_c
    select case(benchmark)
    case(M0_V02)
      gamma_squared=ELEMENTARY_CHARGE*coefficients%reference_n*coefficients%mobility_n/ &
        (coefficients%diffusion_n*bar_kappa)
      gamma=sqrt(gamma_squared)
      phi=coefficients%voltage*sinh(gamma*x)/sinh(gamma*coefficients%length)
      phi_derivative=coefficients%voltage*gamma*cosh(gamma*x)/ &
        sinh(gamma*coefficients%length)
      if (field == F_PHI) then
        value=phi; derivative=phi_derivative
      else if (field == F_U1) then
        value=-coefficients%piezo_e*phi/coefficients%elastic_c
        derivative=-coefficients%piezo_e*phi_derivative/coefficients%elastic_c
      end if
    case(M0_V03)
      phi=coefficients%voltage*x/coefficients%length
      phi_derivative=coefficients%voltage/coefficients%length
      ratio=coefficients%mobility_n/coefficients%diffusion_n
      if (field == F_PHI) then
        value=phi; derivative=phi_derivative
      else if (field == F_ETA_N) then
        value=ratio*phi; derivative=ratio*phi_derivative
      end if
    case(M0_V04)
      phi=coefficients%voltage*x/coefficients%length
      phi_derivative=coefficients%voltage/coefficients%length
      ratio=coefficients%mobility_p/coefficients%diffusion_p
      if (field == F_PHI) then
        value=phi; derivative=phi_derivative
      else if (field == F_ETA_P) then
        value=-ratio*phi; derivative=-ratio*phi_derivative
      end if
    case(M0_V05_LINEAR)
      if (field == F_THETA) then
        value=coefficients%delta_temperature*x/coefficients%length
        derivative=coefficients%delta_temperature/coefficients%length
      end if
    case(M0_V05_UNIFORM)
      if (field == F_THETA) then
        value=coefficients%uniform_heat_source*x*(coefficients%length-x)/ &
          (2.0_dp*coefficients%thermal_conductivity)
        derivative=coefficients%uniform_heat_source*(coefficients%length-2.0_dp*x)/ &
          (2.0_dp*coefficients%thermal_conductivity)
      end if
    case(M0_V06,M0_A1)
      gamma_squared=ELEMENTARY_CHARGE*coefficients%reference_n*coefficients%mobility_n/ &
        (coefficients%diffusion_n*bar_kappa)
      if (benchmark == M0_A1) gamma_squared=gamma_squared+ &
        ELEMENTARY_CHARGE*coefficients%reference_p*coefficients%mobility_p/ &
        (coefficients%diffusion_p*bar_kappa)
      gamma=sqrt(gamma_squared)
      s=(coefficients%pyroelectric+coefficients%piezo_e* &
        coefficients%thermal_expansion)*coefficients%delta_temperature/ &
        (bar_kappa*coefficients%length)
      phi=s/gamma_squared*(cosh(gamma*(x-0.5_dp*coefficients%length))/ &
        cosh(0.5_dp*gamma*coefficients%length)-1.0_dp)
      phi_derivative=s/gamma*sinh(gamma*(x-0.5_dp*coefficients%length))/ &
        cosh(0.5_dp*gamma*coefficients%length)
      ratio=coefficients%mobility_n/coefficients%diffusion_n
      select case(field)
      case(F_U1)
        value=coefficients%thermal_expansion*coefficients%delta_temperature*x**2/ &
          (2.0_dp*coefficients%length)-coefficients%piezo_e*phi/coefficients%elastic_c
        derivative=coefficients%thermal_expansion*coefficients%delta_temperature*x/ &
          coefficients%length-coefficients%piezo_e*phi_derivative/coefficients%elastic_c
      case(F_PHI)
        value=phi; derivative=phi_derivative
      case(F_ETA_N)
        value=ratio*phi; derivative=ratio*phi_derivative
      case(F_ETA_P)
        if (benchmark == M0_A1) then
          ratio=coefficients%mobility_p/coefficients%diffusion_p
          value=-ratio*phi; derivative=-ratio*phi_derivative
        end if
      case(F_THETA)
        value=coefficients%delta_temperature*x/coefficients%length
        derivative=coefficients%delta_temperature/coefficients%length
      end select
    case default
      call set_status(status,STATUS_INVALID_ARGUMENT,'Exact M0 field is undefined.')
    end select
  end subroutine exact_m0_field

  subroutine build_bar_mesh(element_count,length,width,mesh)
    integer, intent(in) :: element_count
    real(dp), intent(in) :: length,width
    type(MeshType), intent(out) :: mesh
    real(dp) :: natural_nodes(3,HEX20_NNODE),coordinates(3,HEX20_NNODE*element_count)
    real(dp) :: target(3),element_length,tolerance
    integer :: connectivity(HEX20_NNODE,element_count)
    integer :: element,local_node,node,node_count

    call hex20_natural_coordinates(natural_nodes)
    coordinates=0.0_dp
    connectivity=0
    node_count=0
    element_length=length/real(element_count,dp)
    tolerance=1.0e-12_dp*max(length,width)
    do element=1,element_count
      do local_node=1,HEX20_NNODE
        target(1)=real(element-1,dp)*element_length+ &
          0.5_dp*element_length*(natural_nodes(1,local_node)+1.0_dp)
        target(2)=0.5_dp*width*(natural_nodes(2,local_node)+1.0_dp)
        target(3)=0.5_dp*width*(natural_nodes(3,local_node)+1.0_dp)
        node=locate_node(coordinates,node_count,target,tolerance)
        if (node == 0) then
          node_count=node_count+1
          coordinates(:,node_count)=target
          node=node_count
        end if
        connectivity(local_node,element)=node
      end do
    end do
    allocate(mesh%coordinates(3,node_count),mesh%connectivity(HEX20_NNODE,element_count), &
      mesh%material_id(element_count))
    mesh%coordinates=coordinates(:,1:node_count)
    mesh%connectivity=connectivity
    mesh%material_id=1
  end subroutine build_bar_mesh

  integer function locate_node(coordinates,node_count,target,tolerance) result(node_id)
    real(dp), intent(in) :: coordinates(:,:),target(3),tolerance
    integer, intent(in) :: node_count
    integer :: node
    node_id=0
    do node=1,node_count
      if (maxval(abs(coordinates(:,node)-target)) <= tolerance) then
        node_id=node
        return
      end if
    end do
  end function locate_node

  subroutine benchmark_field_scales(benchmark,coefficients,scales)
    integer, intent(in) :: benchmark
    type(M0CoefficientsType), intent(in) :: coefficients
    real(dp), intent(out) :: scales(NDOF_NODE)
    real(dp) :: phi_scale,eta_scale,temperature_scale

    phi_scale=max(abs(coefficients%voltage),1.0e-4_dp)
    if (benchmark == M0_V06) phi_scale=max(phi_scale,1.0e-4_dp)
    eta_scale=max(coefficients%mobility_n/coefficients%diffusion_n*phi_scale,1.0e-3_dp)
    temperature_scale=max(abs(coefficients%delta_temperature), &
      abs(coefficients%uniform_heat_source*coefficients%length**2/ &
      (8.0_dp*coefficients%thermal_conductivity)),1.0_dp)
    scales=[1.0e-13_dp,1.0e-13_dp,1.0e-13_dp,phi_scale,eta_scale,eta_scale, &
      temperature_scale]
  end subroutine benchmark_field_scales

  logical function active_field(benchmark,field)
    integer, intent(in) :: benchmark,field
    select case(benchmark)
    case(M0_V02)
      active_field=field == F_U1 .or. field == F_PHI
    case(M0_V03)
      active_field=field == F_ETA_N
    case(M0_V04)
      active_field=field == F_ETA_P
    case(M0_V05_LINEAR,M0_V05_UNIFORM)
      active_field=field == F_THETA
    case(M0_V06)
      active_field=field == F_U1 .or. field == F_PHI .or. &
        field == F_ETA_N .or. field == F_THETA
    case(M0_A1)
      active_field=field == F_U1 .or. field == F_PHI .or. &
        field == F_ETA_N .or. field == F_ETA_P .or. field == F_THETA
    case default
      active_field=.false.
    end select
  end function active_field

  logical function valid_benchmark(benchmark)
    integer, intent(in) :: benchmark
    valid_benchmark=benchmark == M0_V02 .or. benchmark == M0_V03 .or. &
      benchmark == M0_V04 .or. benchmark == M0_V05_LINEAR .or. &
      benchmark == M0_V05_UNIFORM .or. benchmark == M0_V06 .or. benchmark == M0_A1
  end function valid_benchmark

  subroutine validate_coefficients(coefficients,status)
    type(M0CoefficientsType), intent(in) :: coefficients
    type(SolverStatusType), intent(out) :: status
    call clear_status(status)
    if (coefficients%length <= 0.0_dp .or. coefficients%width <= 0.0_dp .or. &
        coefficients%elastic_c <= 0.0_dp .or. coefficients%permittivity <= 0.0_dp .or. &
        coefficients%reference_n <= 0.0_dp .or. coefficients%reference_p <= 0.0_dp .or. &
        coefficients%mobility_n <= 0.0_dp .or. coefficients%mobility_p <= 0.0_dp .or. &
        coefficients%diffusion_n <= 0.0_dp .or. coefficients%diffusion_p <= 0.0_dp .or. &
        coefficients%thermal_conductivity <= 0.0_dp) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'M0 coefficients must be positive where required.')
    end if
  end subroutine validate_coefficients

  pure function cross_product(a,b) result(c)
    real(dp), intent(in) :: a(3),b(3)
    real(dp) :: c(3)
    c(1)=a(2)*b(3)-a(3)*b(2)
    c(2)=a(3)*b(1)-a(1)*b(3)
    c(3)=a(1)*b(2)-a(2)*b(1)
  end function cross_product

  pure real(dp) function vector_norm(vector) result(value)
    real(dp), intent(in) :: vector(:)
    value=sqrt(sum(vector*vector))
  end function vector_norm

end module m0_benchmark_m
