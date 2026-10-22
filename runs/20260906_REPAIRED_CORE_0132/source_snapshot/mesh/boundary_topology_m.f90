module boundary_topology_m
  use kinds_m, only : dp
  use mesh_types_m, only : MeshType, HEX20_NNODE
  use status_m, only : SolverStatusType, clear_status, set_status, &
    STATUS_INPUT_ERROR, STATUS_GEOMETRY_ERROR
  implicit none
  private

  integer, parameter :: FACE_COUNT=6
  integer, parameter :: FACE_NNODE=8
  integer, parameter :: FACE_LOCAL(FACE_NNODE,FACE_COUNT)=reshape([ &
    1,2,3,4,9,10,11,12, &
    5,6,7,8,13,14,15,16, &
    1,2,6,5,9,18,13,17, &
    2,3,7,6,10,19,14,18, &
    3,4,8,7,11,20,15,19, &
    4,1,5,8,12,17,16,20],[FACE_NNODE,FACE_COUNT])

  type, public :: BoundaryTopologyAuditType
    integer :: boundary_face_count=0
    integer :: internal_face_count=0
    integer :: nonmanifold_face_count=0
    integer :: left_face_count=0
    integer :: right_face_count=0
    integer :: front_face_count=0
    integer :: back_face_count=0
    integer :: bottom_face_count=0
    integer :: top_face_count=0
    integer :: hole_face_count=0
    integer :: unknown_face_count=0
    integer :: left_node_count=0
    integer :: right_node_count=0
    integer :: hole_node_count=0
    real(dp) :: bounds_min(3)=0.0_dp
    real(dp) :: bounds_max(3)=0.0_dp
    real(dp) :: center(2)=0.0_dp
    real(dp) :: classification_tolerance=0.0_dp
    real(dp) :: minimum_radial_distance=huge(1.0_dp)
    real(dp) :: maximum_hole_radius_error=0.0_dp
  end type BoundaryTopologyAuditType

  public :: audit_boundary_topology

contains

  subroutine audit_boundary_topology(mesh,hole_radius,audit,status)
    type(MeshType), intent(in) :: mesh
    real(dp), intent(in) :: hole_radius
    type(BoundaryTopologyAuditType), intent(out) :: audit
    type(SolverStatusType), intent(out) :: status
    integer, allocatable :: face_nodes(:,:),sorted_corners(:,:),multiplicity(:)
    logical, allocatable :: left_nodes(:),right_nodes(:),hole_nodes(:)
    integer :: maximum_faces,unique_faces,element,local_face,face,match,node
    integer :: candidate(FACE_NNODE),corners(4)
    real(dp) :: scale,radius

    call clear_status(status)
    audit=BoundaryTopologyAuditType()
    if (.not. allocated(mesh%coordinates) .or. .not. allocated(mesh%connectivity)) then
      call set_status(status,STATUS_INPUT_ERROR,'Boundary topology received an empty mesh.')
      return
    end if
    if (size(mesh%connectivity,1) /= HEX20_NNODE .or. hole_radius < 0.0_dp) then
      call set_status(status,STATUS_INPUT_ERROR,'Boundary topology input is invalid.')
      return
    end if

    audit%bounds_min=minval(mesh%coordinates,dim=2)
    audit%bounds_max=maxval(mesh%coordinates,dim=2)
    audit%center=0.5_dp*(audit%bounds_min(1:2)+audit%bounds_max(1:2))
    scale=maxval(audit%bounds_max-audit%bounds_min)
    audit%classification_tolerance=max(1.0e-9_dp*scale,1.0e-18_dp)
    do node=1,mesh%number_of_nodes()
      radius=hypot(mesh%coordinates(1,node)-audit%center(1), &
        mesh%coordinates(2,node)-audit%center(2))
      audit%minimum_radial_distance=min(audit%minimum_radial_distance,radius)
    end do

    maximum_faces=FACE_COUNT*mesh%number_of_elements()
    allocate(face_nodes(FACE_NNODE,maximum_faces),sorted_corners(4,maximum_faces), &
      multiplicity(maximum_faces))
    face_nodes=0; sorted_corners=0; multiplicity=0; unique_faces=0
    do element=1,mesh%number_of_elements()
      do local_face=1,FACE_COUNT
        candidate=mesh%connectivity(FACE_LOCAL(:,local_face),element)
        corners=candidate(1:4)
        call sort_four(corners)
        match=0
        do face=1,unique_faces
          if (all(sorted_corners(:,face) == corners)) then
            match=face
            exit
          end if
        end do
        if (match == 0) then
          unique_faces=unique_faces+1
          sorted_corners(:,unique_faces)=corners
          face_nodes(:,unique_faces)=candidate
          multiplicity(unique_faces)=1
        else
          multiplicity(match)=multiplicity(match)+1
        end if
      end do
    end do

    allocate(left_nodes(mesh%number_of_nodes()),right_nodes(mesh%number_of_nodes()), &
      hole_nodes(mesh%number_of_nodes()))
    left_nodes=.false.; right_nodes=.false.; hole_nodes=.false.
    do face=1,unique_faces
      select case(multiplicity(face))
      case(1)
        audit%boundary_face_count=audit%boundary_face_count+1
        call classify_boundary_face(mesh,face_nodes(:,face),hole_radius,audit, &
          left_nodes,right_nodes,hole_nodes)
      case(2)
        audit%internal_face_count=audit%internal_face_count+1
      case default
        audit%nonmanifold_face_count=audit%nonmanifold_face_count+1
      end select
    end do
    audit%left_node_count=count(left_nodes)
    audit%right_node_count=count(right_nodes)
    audit%hole_node_count=count(hole_nodes)
    if (audit%nonmanifold_face_count /= 0) then
      call set_status(status,STATUS_GEOMETRY_ERROR, &
        'Mesh contains a face shared by more than two HEX20 elements.')
    end if
  end subroutine audit_boundary_topology

  subroutine classify_boundary_face(mesh,nodes,hole_radius,audit,left_nodes, &
                                    right_nodes,hole_nodes)
    type(MeshType), intent(in) :: mesh
    integer, intent(in) :: nodes(FACE_NNODE)
    real(dp), intent(in) :: hole_radius
    type(BoundaryTopologyAuditType), intent(inout) :: audit
    logical, intent(inout) :: left_nodes(:),right_nodes(:),hole_nodes(:)
    real(dp) :: coordinates(3,FACE_NNODE),radii(FACE_NNODE),tolerance
    integer :: entry

    coordinates=mesh%coordinates(:,nodes)
    tolerance=audit%classification_tolerance
    do entry=1,FACE_NNODE
      radii(entry)=hypot(coordinates(1,entry)-audit%center(1), &
        coordinates(2,entry)-audit%center(2))
    end do
    if (all(abs(coordinates(1,:)-audit%bounds_min(1)) <= tolerance)) then
      audit%left_face_count=audit%left_face_count+1
      left_nodes(nodes)=.true.
    else if (all(abs(coordinates(1,:)-audit%bounds_max(1)) <= tolerance)) then
      audit%right_face_count=audit%right_face_count+1
      right_nodes(nodes)=.true.
    else if (all(abs(coordinates(2,:)-audit%bounds_min(2)) <= tolerance)) then
      audit%front_face_count=audit%front_face_count+1
    else if (all(abs(coordinates(2,:)-audit%bounds_max(2)) <= tolerance)) then
      audit%back_face_count=audit%back_face_count+1
    else if (all(abs(coordinates(3,:)-audit%bounds_min(3)) <= tolerance)) then
      audit%bottom_face_count=audit%bottom_face_count+1
    else if (all(abs(coordinates(3,:)-audit%bounds_max(3)) <= tolerance)) then
      audit%top_face_count=audit%top_face_count+1
    else if (hole_radius > 0.0_dp .and. &
             all(abs(radii-hole_radius) <= tolerance)) then
      audit%hole_face_count=audit%hole_face_count+1
      hole_nodes(nodes)=.true.
      audit%maximum_hole_radius_error=max(audit%maximum_hole_radius_error, &
        maxval(abs(radii-hole_radius)))
    else
      audit%unknown_face_count=audit%unknown_face_count+1
    end if
  end subroutine classify_boundary_face

  subroutine sort_four(values)
    integer, intent(inout) :: values(4)
    integer :: left,right,temporary
    do left=1,3
      do right=left+1,4
        if (values(right) < values(left)) then
          temporary=values(left); values(left)=values(right); values(right)=temporary
        end if
      end do
    end do
  end subroutine sort_four

end module boundary_topology_m
