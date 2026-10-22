module sparse_pattern_m
  use kinds_m, only : dp, i8
  use status_m, only : SolverStatusType, clear_status, set_status, &
    STATUS_INPUT_ERROR, STATUS_INVALID_ARGUMENT
  use field_layout_m, only : NDOF_NODE, global_dof
  use mesh_types_m, only : MeshType, HEX20_NNODE
  use dof_map_m, only : HEX20_NDOF, build_element_dof_map
  implicit none
  private

  type, public :: SparseAssemblyContextType
    integer :: number_of_equations = 0
    integer :: number_of_nonzeros = 0
    integer, allocatable :: ia(:)
    integer, allocatable :: ja(:)
    integer, allocatable :: element_dof(:,:)             ! (140,nelem)
    integer, allocatable :: element_neighbor_offset(:,:,:) ! (20,20,nelem)
  end type SparseAssemblyContextType

  public :: initialize_sparse_context, csr_matvec, find_csr_position

contains

  subroutine initialize_sparse_context(mesh, context, status)
    type(MeshType), intent(in) :: mesh
    type(SparseAssemblyContextType), intent(out) :: context
    type(SolverStatusType), intent(out) :: status
    integer(i8), allocatable :: pair_keys(:)
    integer, allocatable :: node_neighbor_count(:), node_ia(:), node_ja(:), next_position(:)
    integer :: node_count, element_count, pair_count, unique_count
    integer :: element, local_row, local_column, position, row_node, column_node
    integer :: row_field, column_field, global_row, node_position, offset
    integer(i8) :: key, previous_key
    type(SolverStatusType) :: local_status

    call clear_status(status)
    node_count = mesh%number_of_nodes()
    element_count = mesh%number_of_elements()
    if (node_count < HEX20_NNODE .or. element_count < 1) then
      call set_status(status, STATUS_INPUT_ERROR, 'Cannot build CSR pattern from an empty HEX20 mesh.')
      return
    end if
    if (.not. allocated(mesh%connectivity)) then
      call set_status(status, STATUS_INPUT_ERROR, 'Mesh connectivity is not allocated.')
      return
    end if

    pair_count = element_count * HEX20_NNODE * HEX20_NNODE
    allocate(pair_keys(pair_count))
    position = 0
    do element = 1, element_count
      do local_row = 1, HEX20_NNODE
        row_node = mesh%connectivity(local_row,element)
        do local_column = 1, HEX20_NNODE
          column_node = mesh%connectivity(local_column,element)
          position = position + 1
          pair_keys(position) = int(row_node-1,i8)*int(node_count,i8) + int(column_node,i8)
        end do
      end do
    end do
    call merge_sort_i8(pair_keys)

    unique_count = 0
    previous_key = -1_i8
    do position = 1, pair_count
      if (pair_keys(position) /= previous_key) then
        unique_count = unique_count + 1
        previous_key = pair_keys(position)
      end if
    end do
    allocate(node_neighbor_count(node_count), node_ia(node_count+1), node_ja(unique_count))
    node_neighbor_count = 0
    previous_key = -1_i8
    do position = 1, pair_count
      key = pair_keys(position)
      if (key == previous_key) cycle
      row_node = int((key-1_i8)/int(node_count,i8)) + 1
      node_neighbor_count(row_node) = node_neighbor_count(row_node) + 1
      previous_key = key
    end do
    node_ia(1) = 1
    do row_node = 1, node_count
      node_ia(row_node+1) = node_ia(row_node) + node_neighbor_count(row_node)
    end do
    allocate(next_position(node_count))
    next_position = node_ia(1:node_count)
    previous_key = -1_i8
    do position = 1, pair_count
      key = pair_keys(position)
      if (key == previous_key) cycle
      row_node = int((key-1_i8)/int(node_count,i8)) + 1
      column_node = int(mod(key-1_i8,int(node_count,i8))) + 1
      node_ja(next_position(row_node)) = column_node
      next_position(row_node) = next_position(row_node) + 1
      previous_key = key
    end do

    context%number_of_equations = NDOF_NODE * node_count
    context%number_of_nonzeros = unique_count * NDOF_NODE * NDOF_NODE
    allocate(context%ia(context%number_of_equations+1))
    allocate(context%ja(context%number_of_nonzeros))
    position = 1
    do row_node = 1, node_count
      do row_field = 1, NDOF_NODE
        global_row = global_dof(row_node,row_field)
        context%ia(global_row) = position
        do node_position = node_ia(row_node), node_ia(row_node+1)-1
          column_node = node_ja(node_position)
          do column_field = 1, NDOF_NODE
            context%ja(position) = global_dof(column_node,column_field)
            position = position + 1
          end do
        end do
      end do
    end do
    context%ia(context%number_of_equations+1) = position
    if (position-1 /= context%number_of_nonzeros) then
      call set_status(status, STATUS_INVALID_ARGUMENT, 'Internal CSR nonzero count mismatch.')
      return
    end if

    allocate(context%element_dof(HEX20_NDOF,element_count))
    allocate(context%element_neighbor_offset(HEX20_NNODE,HEX20_NNODE,element_count))
    do element = 1, element_count
      call build_element_dof_map(mesh%connectivity(:,element), &
        context%element_dof(:,element), local_status)
      if (.not. local_status%is_ok()) then
        call set_status(status, local_status%code, trim(local_status%message))
        return
      end if
      do local_row = 1, HEX20_NNODE
        row_node = mesh%connectivity(local_row,element)
        do local_column = 1, HEX20_NNODE
          column_node = mesh%connectivity(local_column,element)
          node_position = find_in_sorted_row(node_ja, node_ia(row_node), &
            node_ia(row_node+1)-1, column_node)
          if (node_position == 0) then
            call set_status(status, STATUS_INVALID_ARGUMENT, &
              'Element node pair was not found in the CSR node graph.')
            return
          end if
          offset = (node_position-node_ia(row_node))*NDOF_NODE
          context%element_neighbor_offset(local_row,local_column,element) = offset
        end do
      end do
    end do
  end subroutine initialize_sparse_context

  subroutine csr_matvec(context, values, vector, product, status)
    type(SparseAssemblyContextType), intent(in) :: context
    real(dp), intent(in) :: values(:), vector(:)
    real(dp), intent(out) :: product(:)
    type(SolverStatusType), intent(out) :: status
    integer :: row, position

    call clear_status(status)
    if (size(values) /= context%number_of_nonzeros .or. &
        size(vector) /= context%number_of_equations .or. &
        size(product) /= context%number_of_equations) then
      product = 0.0_dp
      call set_status(status, STATUS_INVALID_ARGUMENT, 'CSR matvec array size mismatch.')
      return
    end if
    product = 0.0_dp
    do row = 1, context%number_of_equations
      do position = context%ia(row), context%ia(row+1)-1
        product(row) = product(row) + values(position)*vector(context%ja(position))
      end do
    end do
  end subroutine csr_matvec

  integer function find_csr_position(context, row, column) result(position)
    type(SparseAssemblyContextType), intent(in) :: context
    integer, intent(in) :: row, column
    integer :: low, high, middle

    position = 0
    if (row < 1 .or. row > context%number_of_equations) return
    low = context%ia(row)
    high = context%ia(row+1)-1
    do while (low <= high)
      middle = (low+high)/2
      if (context%ja(middle) == column) then
        position = middle
        return
      else if (context%ja(middle) < column) then
        low = middle+1
      else
        high = middle-1
      end if
    end do
  end function find_csr_position

  integer function find_in_sorted_row(values, first, last, target) result(position)
    integer, intent(in) :: values(:), first, last, target
    integer :: low, high, middle

    position = 0
    low = first
    high = last
    do while (low <= high)
      middle = (low+high)/2
      if (values(middle) == target) then
        position = middle
        return
      else if (values(middle) < target) then
        low = middle+1
      else
        high = middle-1
      end if
    end do
  end function find_in_sorted_row

  subroutine merge_sort_i8(values)
    integer(i8), intent(inout) :: values(:)
    integer(i8), allocatable :: buffer(:)
    integer :: width, left, middle, right, i, j, output_position, count

    count = size(values)
    if (count <= 1) return
    allocate(buffer(count))
    width = 1
    do while (width < count)
      left = 1
      do while (left <= count)
        middle = min(left+width-1,count)
        right = min(left+2*width-1,count)
        if (middle >= right) then
          buffer(left:right) = values(left:right)
        else
          i = left
          j = middle+1
          do output_position = left, right
            if (i > middle) then
              buffer(output_position) = values(j)
              j = j+1
            else if (j > right) then
              buffer(output_position) = values(i)
              i = i+1
            else if (values(i) <= values(j)) then
              buffer(output_position) = values(i)
              i = i+1
            else
              buffer(output_position) = values(j)
              j = j+1
            end if
          end do
        end if
        left = left+2*width
      end do
      values = buffer
      if (width > count/2) exit
      width = 2*width
    end do
  end subroutine merge_sort_i8

end module sparse_pattern_m
