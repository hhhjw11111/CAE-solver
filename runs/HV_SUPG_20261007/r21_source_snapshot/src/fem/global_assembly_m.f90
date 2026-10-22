module global_assembly_m
  use kinds_m, only : dp
  use status_m, only : SolverStatusType, clear_status, set_status, &
    STATUS_INPUT_ERROR, STATUS_INVALID_ARGUMENT
  use field_layout_m, only : NDOF_NODE, local_dof
  use mesh_types_m, only : MeshType, HEX20_NNODE
  use dof_map_m, only : HEX20_NDOF
  use model_types_m, only : MaterialType, StateType, LoadStateType, DiagnosticsType
  use element_kernel_m, only : ElementRequestType, ElementDiagnosticsType, evaluate_hex20_element
  use sparse_pattern_m, only : SparseAssemblyContextType
  use dirichlet_m, only : DirichletBCType, apply_dirichlet_increment_rows
  implicit none
  private

  type, public :: GlobalAssemblyRequestType
    logical :: need_jacobian = .true.
  end type GlobalAssemblyRequestType

  type, public :: GlobalSystemType
    real(dp), allocatable :: residual(:)
    real(dp), allocatable :: jacobian_values(:)
  end type GlobalSystemType

  public :: assemble_global_system

contains

  subroutine assemble_global_system(mesh, materials, state, load, context, &
                                    boundary_conditions, request, system, diagnostics, status)
    type(MeshType), intent(in) :: mesh
    type(MaterialType), intent(in) :: materials(:)
    type(StateType), intent(in) :: state
    type(LoadStateType), intent(in) :: load
    type(SparseAssemblyContextType), intent(in) :: context
    type(DirichletBCType), intent(in), optional :: boundary_conditions(:)
    type(GlobalAssemblyRequestType), intent(in) :: request
    type(GlobalSystemType), intent(inout) :: system
    type(DiagnosticsType), intent(out) :: diagnostics
    type(SolverStatusType), intent(out) :: status
    real(dp) :: coordinates_e(3,HEX20_NNODE), state_e(HEX20_NDOF)
    real(dp) :: residual_e(HEX20_NDOF), jacobian_e(HEX20_NDOF,HEX20_NDOF)
    type(ElementRequestType) :: element_request
    type(ElementDiagnosticsType) :: element_diagnostics
    type(SolverStatusType) :: local_status
    integer :: element, material_id, local_row_node, local_column_node
    integer :: row_field, column_field, local_row, local_column, global_row
    integer :: csr_position, offset
    integer, parameter :: batch_size=32
    integer :: slot
    logical :: batch_transport
    real(dp), allocatable :: batch_residual(:,:),batch_jacobian(:,:,:)
    type(ElementDiagnosticsType), allocatable :: batch_diagnostics(:)
    type(SolverStatusType), allocatable :: batch_status(:)

    call clear_status(status)
    diagnostics = DiagnosticsType()
    if (context%number_of_equations /= NDOF_NODE*mesh%number_of_nodes()) then
      call set_status(status, STATUS_INVALID_ARGUMENT, 'Sparse context does not match the mesh.')
      return
    end if
    if (.not. allocated(state%x) .or. size(state%x) /= context%number_of_equations) then
      call set_status(status, STATUS_INVALID_ARGUMENT, 'Global state size does not match the CSR context.')
      return
    end if
    if (size(materials) < 1 .or. .not. allocated(mesh%material_id)) then
      call set_status(status, STATUS_INPUT_ERROR, 'Materials or element material identifiers are missing.')
      return
    end if
    call ensure_system_storage(context,system)
    system%residual = 0.0_dp
    system%jacobian_values = 0.0_dp
    element_request%need_jacobian = request%need_jacobian
    batch_transport=any(materials%transport_supg)
    if(batch_transport) then
      allocate(batch_residual(HEX20_NDOF,batch_size),batch_jacobian(HEX20_NDOF,HEX20_NDOF,batch_size), &
        batch_diagnostics(batch_size),batch_status(batch_size))
    end if

    do element = 1, mesh%number_of_elements()
      material_id = mesh%material_id(element)
      if (material_id < 1 .or. material_id > size(materials)) then
        call set_status(status, STATUS_INPUT_ERROR, 'Element material identifier is invalid.')
        return
      end if
      if(batch_transport) then
        slot=mod(element-1,batch_size)+1
        if(slot==1) call evaluate_transport_batch(mesh,materials,state,load,context,element_request, &
          element,min(batch_size,mesh%number_of_elements()-element+1), &
          batch_residual,batch_jacobian,batch_diagnostics,batch_status)
        residual_e=batch_residual(:,slot);jacobian_e=batch_jacobian(:,:,slot)
        element_diagnostics=batch_diagnostics(slot);local_status=batch_status(slot)
      else
        coordinates_e = mesh%coordinates(:,mesh%connectivity(:,element))
        state_e = state%x(context%element_dof(:,element))
        call evaluate_hex20_element(materials(material_id),coordinates_e,state_e,load, &
          element_request,residual_e,jacobian_e,element_diagnostics,local_status)
      end if
      if (.not. local_status%is_ok()) then
        call set_status(status,local_status%code,'Element '//integer_to_string(element)// &
          ' assembly failed: '//trim(local_status%message))
        return
      end if
      system%residual(context%element_dof(:,element)) = &
        system%residual(context%element_dof(:,element)) + residual_e

      if (request%need_jacobian) then
        do local_row_node = 1, HEX20_NNODE
          do row_field = 1, NDOF_NODE
            local_row = local_dof(local_row_node,row_field)
            global_row = context%element_dof(local_row,element)
            do local_column_node = 1, HEX20_NNODE
              offset = context%element_neighbor_offset(local_row_node,local_column_node,element)
              csr_position = context%ia(global_row) + offset
              do column_field = 1, NDOF_NODE
                local_column = local_dof(local_column_node,column_field)
                system%jacobian_values(csr_position+column_field-1) = &
                  system%jacobian_values(csr_position+column_field-1) + &
                  jacobian_e(local_row,local_column)
              end do
            end do
          end do
        end do
      end if
      diagnostics%min_temperature = min(diagnostics%min_temperature, &
        element_diagnostics%min_temperature)
      diagnostics%max_temperature = max(diagnostics%max_temperature, &
        element_diagnostics%max_temperature)
      diagnostics%min_n = min(diagnostics%min_n,element_diagnostics%min_n)
      diagnostics%max_n = max(diagnostics%max_n,element_diagnostics%max_n)
      diagnostics%min_p = min(diagnostics%min_p,element_diagnostics%min_p)
      diagnostics%max_p = max(diagnostics%max_p,element_diagnostics%max_p)
    end do

    if (present(boundary_conditions)) then
      call apply_dirichlet_increment_rows(context,boundary_conditions,system%residual, &
        system%jacobian_values,local_status)
      if (.not. local_status%is_ok()) then
        call set_status(status,local_status%code,trim(local_status%message))
        return
      end if
    end if
    call compute_field_residual_norms(system%residual,diagnostics%residual_by_field)
  end subroutine assemble_global_system

  subroutine evaluate_transport_batch(mesh,materials,state,load,context,request,first,count,r,k,d,s)
    type(MeshType), intent(in) :: mesh
    type(MaterialType), intent(in) :: materials(:)
    type(StateType), intent(in) :: state
    type(LoadStateType), intent(in) :: load
    type(SparseAssemblyContextType), intent(in) :: context
    type(ElementRequestType), intent(in) :: request
    integer, intent(in) :: first,count
    real(dp), intent(out) :: r(:,:),k(:,:,:)
    type(ElementDiagnosticsType), intent(out) :: d(:)
    type(SolverStatusType), intent(out) :: s(:)
    real(dp) :: coords(3,HEX20_NNODE),x(HEX20_NDOF)
    integer :: i,e,mat
    ! Independent element kernels only. The caller scatters in original element
    ! order, so shared-node sums have identical ordering at every thread count.
    !$omp parallel do default(shared) private(i,e,mat,coords,x) schedule(static)
    do i=1,count
      e=first+i-1;mat=mesh%material_id(e)
      if(mat<1 .or.mat>size(materials)) then
        r(:,i)=0.0_dp;k(:,:,i)=0.0_dp;d(i)=ElementDiagnosticsType()
        call set_status(s(i),STATUS_INPUT_ERROR,'Element material identifier is invalid.')
      else
        coords=mesh%coordinates(:,mesh%connectivity(:,e));x=state%x(context%element_dof(:,e))
        call evaluate_hex20_element(materials(mat),coords,x,load,request,r(:,i),k(:,:,i),d(i),s(i))
      end if
    end do
    !$omp end parallel do
  end subroutine evaluate_transport_batch

  subroutine ensure_system_storage(context,system)
    type(SparseAssemblyContextType), intent(in) :: context
    type(GlobalSystemType), intent(inout) :: system

    if (allocated(system%residual)) then
      if (size(system%residual) /= context%number_of_equations) deallocate(system%residual)
    end if
    if (.not. allocated(system%residual)) allocate(system%residual(context%number_of_equations))
    if (allocated(system%jacobian_values)) then
      if (size(system%jacobian_values) /= context%number_of_nonzeros) &
        deallocate(system%jacobian_values)
    end if
    if (.not. allocated(system%jacobian_values)) &
      allocate(system%jacobian_values(context%number_of_nonzeros))
  end subroutine ensure_system_storage

  subroutine compute_field_residual_norms(residual,norms)
    real(dp), intent(in) :: residual(:)
    real(dp), intent(out) :: norms(NDOF_NODE)
    integer :: equation, field

    norms = 0.0_dp
    do equation = 1, size(residual)
      field = mod(equation-1,NDOF_NODE)+1
      norms(field) = norms(field)+residual(equation)*residual(equation)
    end do
    norms = sqrt(norms)
  end subroutine compute_field_residual_norms

  function integer_to_string(value) result(text)
    integer, intent(in) :: value
    character(len=32) :: text
    write(text,'(i0)') value
  end function integer_to_string

end module global_assembly_m
