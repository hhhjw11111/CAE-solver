module independent_flux_audit_m
  use kinds_m, only : dp
  use constants_m, only : ELEMENTARY_CHARGE
  use field_layout_m, only : NDOF_NODE, F_PHI, F_ETA_N, F_ETA_P, F_THETA, global_dof
  use mesh_types_m, only : MeshType, HEX20_NNODE, BoundaryFaceType
  use model_types_m, only : MaterialType, StateType, LoadStateType
  use status_m, only : SolverStatusType, clear_status, set_status, &
    STATUS_INVALID_ARGUMENT, STATUS_GEOMETRY_ERROR
  use hex20_shape_m, only : hex20_shape, hex20_geometry, hex20_natural_coordinates
  use material_point_m, only : MaterialPointStateType, MaterialPointResponseType, &
    MaterialPointTangentType, evaluate_material_point
  implicit none
  private

  ! Boundary sets: 1 = minimum-x contact, 2 = maximum-x contact, 3 = other exterior.
  ! Local face numbering here is x-,x+,y-,y+,z-,z+, not boundary_topology_m numbering.
  type, public :: IndependentFluxReportType
    integer :: exterior_faces=0,internal_faces=0
    real(dp) :: volume=0.0_dp,area=0.0_dp,area_normal(3)=0.0_dp
    real(dp) :: electron_out(3)=0.0_dp,hole_out(3)=0.0_dp,heat_out(3)=0.0_dp
    real(dp) :: electric_displacement_out=0.0_dp,absolute_displacement_out=0.0_dp
    real(dp) :: charge=0.0_dp,absolute_charge=0.0_dp,recombination_charge=0.0_dp
    real(dp) :: joule=0.0_dp,terminal_power=0.0_dp,total_current_work=0.0_dp
    real(dp) :: diffusion_work=0.0_dp,absolute_diffusion_work=0.0_dp
    ! Measured one-sided FE flux on the artificial symmetry cut, not forced zero.
    real(dp) :: symmetry_electron_out=0.0_dp,symmetry_hole_out=0.0_dp,symmetry_heat_out=0.0_dp
    real(dp) :: symmetry_terminal_power=0.0_dp
  end type IndependentFluxReportType

  public :: integrate_independent_fluxes, independent_gauss_rule

contains

  subroutine integrate_independent_fluxes(mesh,materials,state,load,order,report,status,half_y)
    type(MeshType), intent(in) :: mesh
    type(MaterialType), intent(in) :: materials(:)
    type(StateType), intent(in) :: state
    type(LoadStateType), intent(in) :: load
    integer, intent(in) :: order
    type(IndependentFluxReportType), intent(out) :: report
    type(SolverStatusType), intent(out) :: status
    logical, optional, intent(in) :: half_y
    type(BoundaryFaceType), allocatable :: faces(:)
    type(MaterialPointResponseType) :: response
    type(MaterialPointStateType) :: point_state
    real(dp), allocatable :: gp(:),gw(:)
    real(dp) :: coords(3,HEX20_NNODE),shape(HEX20_NNODE),gradient(3,HEX20_NNODE)
    real(dp) :: dn(3,HEX20_NNODE),jacobian(3,3),rst(3),area_vector(3)
    real(dp) :: det,weight,phi,current_out,displacement_out,diffusion(3)
    integer :: element,i,j,k,face,axis,axis_a,axis_b,set_id,mat
    logical :: symmetry_face

    call clear_status(status)
    report=IndependentFluxReportType()
    if (.not. allocated(state%x) .or. mesh%number_of_elements() < 1) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'Flux audit needs a nonempty mesh and state.')
      return
    end if
    if (size(state%x) /= NDOF_NODE*mesh%number_of_nodes() .or. &
        any(mesh%material_id < 1) .or. any(mesh%material_id > size(materials))) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'Flux audit state/material dimensions are invalid.')
      return
    end if
    call independent_gauss_rule(order,gp,gw,status)
    if (.not. status%is_ok()) return
    call exterior_faces(mesh,faces,report%internal_faces,status)
    if (.not. status%is_ok()) return
    report%exterior_faces=size(faces)

    ! No assembled residuals, Dirichlet reactions or production integration weights.
    do element=1,mesh%number_of_elements()
      coords=mesh%coordinates(:,mesh%connectivity(:,element))
      mat=mesh%material_id(element)
      do k=1,order
        do j=1,order
          do i=1,order
            call hex20_geometry(coords,gp(i),gp(j),gp(k),shape,gradient,det,status)
            if (.not. status%is_ok()) then
              call set_status(status,STATUS_GEOMETRY_ERROR, &
                'Flux audit volume geometry failed at element '//trim(integer_to_string(element))// &
                ', quadrature point '//trim(integer_to_string(i))//','//trim(integer_to_string(j))//','// &
                trim(integer_to_string(k))//'.')
              return
            end if
            call sample_state(mesh,element,materials(mat),state,load,shape,gradient, &
              point_state,response,phi,status)
            if (.not. status%is_ok()) return
            weight=gw(i)*gw(j)*gw(k)*det
            report%volume=report%volume+weight
            report%joule=report%joule+weight*response%joule_heat
            report%charge=report%charge+weight*response%charge_density
            report%absolute_charge=report%absolute_charge+weight*abs(response%charge_density)
            report%recombination_charge=report%recombination_charge+ &
              weight*ELEMENTARY_CHARGE*response%recombination
            diffusion=response%current_n+response%current_p-ELEMENTARY_CHARGE*( &
              response%n*matmul(response%mobility_n,point_state%electric_field)+ &
              response%p*matmul(response%mobility_p,point_state%electric_field))
            report%total_current_work=report%total_current_work+weight* &
              dot_product(response%current_n+response%current_p,point_state%electric_field)
            report%diffusion_work=report%diffusion_work+weight* &
              dot_product(diffusion,point_state%electric_field)
            report%absolute_diffusion_work=report%absolute_diffusion_work+weight* &
              abs(dot_product(diffusion,point_state%electric_field))
          end do
        end do
      end do
    end do

    do face=1,size(faces)
      element=faces(face)%owner_element
      coords=mesh%coordinates(:,mesh%connectivity(:,element))
      mat=mesh%material_id(element)
      axis=(faces(face)%local_face+1)/2
      axis_a=mod(axis,3)+1
      axis_b=mod(axis+1,3)+1
      set_id=faces(face)%boundary_set
      symmetry_face=.false.
      if(present(half_y)) then
        if(half_y) symmetry_face=all(abs(mesh%coordinates(2,faces(face)%node_ids))<1.0e-15_dp)
      end if
      rst=0.0_dp
      rst(axis)=real(2*mod(faces(face)%local_face+1,2)-1,dp)
      do j=1,order
        do i=1,order
          rst(axis_a)=gp(i); rst(axis_b)=gp(j)
          call hex20_shape(rst(1),rst(2),rst(3),shape,dn)
          jacobian=matmul(coords,transpose(dn))
          area_vector=rst(axis)*cross_product(jacobian(:,axis_a),jacobian(:,axis_b))
          call hex20_geometry(coords,rst(1),rst(2),rst(3),shape,gradient,det,status)
          if (.not. status%is_ok()) then
            ! Invalid geometry is a hard failure. Moving the sampling point
            ! would change the boundary flux and conceal an inverted layer.
            call set_status(status,STATUS_GEOMETRY_ERROR, &
              'Flux audit face geometry failed at element '//trim(integer_to_string(element))// &
              ', local face '//trim(integer_to_string(faces(face)%local_face))// &
              ', quadrature point '//trim(integer_to_string(i))//','//trim(integer_to_string(j))//'.')
            return
          end if
          call sample_state(mesh,element,materials(mat),state,load,shape,gradient, &
            point_state,response,phi,status)
          if (.not. status%is_ok()) return
          area_vector=area_vector*gw(i)*gw(j)
          report%area=report%area+sqrt(sum(area_vector**2))
          report%area_normal=report%area_normal+area_vector
          report%electron_out(set_id)=report%electron_out(set_id)+ &
            dot_product(response%current_n,area_vector)
          report%hole_out(set_id)=report%hole_out(set_id)+ &
            dot_product(response%current_p,area_vector)
          report%heat_out(set_id)=report%heat_out(set_id)+ &
            dot_product(response%heat_flux,area_vector)
          displacement_out=dot_product(response%electric_displacement,area_vector)
          report%electric_displacement_out=report%electric_displacement_out+displacement_out
          report%absolute_displacement_out=report%absolute_displacement_out+abs(displacement_out)
          current_out=dot_product(response%current_n+response%current_p,area_vector)
          report%terminal_power=report%terminal_power-phi*current_out
          if(symmetry_face) then
            report%symmetry_electron_out=report%symmetry_electron_out+dot_product(response%current_n,area_vector)
            report%symmetry_hole_out=report%symmetry_hole_out+dot_product(response%current_p,area_vector)
            report%symmetry_heat_out=report%symmetry_heat_out+dot_product(response%heat_flux,area_vector)
            report%symmetry_terminal_power=report%symmetry_terminal_power-phi*current_out
          end if
        end do
      end do
    end do
  end subroutine integrate_independent_fluxes

  subroutine sample_state(mesh,element,material,state,load,shape,gradient,point,response,phi,status)
    type(MeshType), intent(in) :: mesh
    integer, intent(in) :: element
    type(MaterialType), intent(in) :: material
    type(StateType), intent(in) :: state
    type(LoadStateType), intent(in) :: load
    real(dp), intent(in) :: shape(HEX20_NNODE),gradient(3,HEX20_NNODE)
    type(MaterialPointStateType), intent(out) :: point
    type(MaterialPointResponseType), intent(out) :: response
    real(dp), intent(out) :: phi
    type(SolverStatusType), intent(out) :: status
    type(MaterialPointTangentType) :: tangent
    real(dp) :: nodal(NDOF_NODE,HEX20_NNODE),du(3,3)
    integer :: node,field

    do node=1,HEX20_NNODE
      do field=1,NDOF_NODE
        nodal(field,node)=state%x(global_dof(mesh%connectivity(node,element),field))
      end do
    end do
    point=MaterialPointStateType()
    du=matmul(nodal(1:3,:),transpose(gradient))
    point%strain=[du(1,1),du(2,2),du(3,3),du(2,3)+du(3,2), &
      du(3,1)+du(1,3),du(1,2)+du(2,1)]
    phi=dot_product(shape,nodal(F_PHI,:))
    point%electric_field=-matmul(gradient,nodal(F_PHI,:))
    point%temperature=material%reference_temperature+dot_product(shape,nodal(F_THETA,:))
    point%temperature_gradient=matmul(gradient,nodal(F_THETA,:))
    point%eta_n=dot_product(shape,nodal(F_ETA_N,:))
    point%eta_p=dot_product(shape,nodal(F_ETA_P,:))
    point%eta_n_gradient=matmul(gradient,nodal(F_ETA_N,:))
    point%eta_p_gradient=matmul(gradient,nodal(F_ETA_P,:))
    call evaluate_material_point(material,point,load,response,tangent,status)
  end subroutine sample_state

  subroutine exterior_faces(mesh,faces,internal_count,status)
    type(MeshType), intent(in) :: mesh
    type(BoundaryFaceType), allocatable, intent(out) :: faces(:)
    integer, intent(out) :: internal_count
    type(SolverStatusType), intent(out) :: status
    type(BoundaryFaceType), allocatable :: candidates(:)
    integer, allocatable :: keys(:,:),corners(:,:),multiplicity(:)
    integer :: nodes(8),key(8),corner(4),element,face,axis,unique_count,match,entry,i,j,temp
    real(dp) :: natural(3,HEX20_NNODE),sign_value,x_min,x_max,tol
    logical :: mask(HEX20_NNODE)

    call clear_status(status)
    internal_count=0; unique_count=0
    call hex20_natural_coordinates(natural)
    allocate(candidates(6*mesh%number_of_elements()),keys(8,6*mesh%number_of_elements()), &
      corners(4,6*mesh%number_of_elements()),multiplicity(6*mesh%number_of_elements()))
    multiplicity=0
    do element=1,mesh%number_of_elements()
      do face=1,6
        axis=(face+1)/2
        sign_value=real(2*mod(face+1,2)-1,dp)
        mask=abs(natural(axis,:)-sign_value) < 0.5_dp
        nodes=pack(mesh%connectivity(:,element),mask)
        key=nodes
        corner=pack(mesh%connectivity(1:8,element),mask(1:8))
        do i=1,7
          do j=i+1,8
            if (key(j) < key(i)) then
              temp=key(j); key(j)=key(i); key(i)=temp
            end if
          end do
        end do
        do i=1,3
          do j=i+1,4
            if (corner(j) < corner(i)) then
              temp=corner(j); corner(j)=corner(i); corner(i)=temp
            end if
          end do
        end do
        match=0
        do entry=1,unique_count
          if (all(corners(:,entry) == corner)) then
            if (any(keys(:,entry) /= key)) then
              call set_status(status,STATUS_GEOMETRY_ERROR,'Flux audit found nonconforming midside face nodes.')
              return
            end if
            match=entry
            exit
          end if
        end do
        if (match == 0) then
          unique_count=unique_count+1
          keys(:,unique_count)=key; corners(:,unique_count)=corner
          candidates(unique_count)%owner_element=element
          candidates(unique_count)%local_face=face
          candidates(unique_count)%node_ids=nodes
          multiplicity(unique_count)=1
        else
          multiplicity(match)=multiplicity(match)+1
          if (multiplicity(match) > 2) then
            call set_status(status,STATUS_GEOMETRY_ERROR,'Flux audit found a nonmanifold face.')
            return
          end if
        end if
      end do
    end do
    internal_count=count(multiplicity(1:unique_count) == 2)
    allocate(faces(count(multiplicity(1:unique_count) == 1)))
    x_min=minval(mesh%coordinates(1,:)); x_max=maxval(mesh%coordinates(1,:))
    tol=max(1.0e-10_dp*(x_max-x_min),1.0e-18_dp)
    face=0
    do entry=1,unique_count
      if (multiplicity(entry) /= 1) cycle
      face=face+1
      faces(face)=candidates(entry)
      faces(face)%boundary_set=3
      if (all(abs(mesh%coordinates(1,faces(face)%node_ids)-x_min) < tol)) then
        faces(face)%boundary_set=1
      else if (all(abs(mesh%coordinates(1,faces(face)%node_ids)-x_max) < tol)) then
        faces(face)%boundary_set=2
      end if
    end do
  end subroutine exterior_faces

  subroutine independent_gauss_rule(order,points,weights,status)
    integer, intent(in) :: order
    real(dp), allocatable, intent(out) :: points(:),weights(:)
    type(SolverStatusType), intent(out) :: status
    call clear_status(status)
    select case(order)
    case(4)
      points=[-0.8611363115940526_dp,-0.3399810435848563_dp, &
        0.3399810435848563_dp,0.8611363115940526_dp]
      weights=[0.3478548451374539_dp,0.6521451548625461_dp, &
        0.6521451548625461_dp,0.3478548451374539_dp]
    case(5)
      points=[-0.9061798459386640_dp,-0.5384693101056831_dp,0.0_dp, &
        0.5384693101056831_dp,0.9061798459386640_dp]
      weights=[0.2369268850561891_dp,0.4786286704993665_dp,0.5688888888888889_dp, &
        0.4786286704993665_dp,0.2369268850561891_dp]
    case default
      call set_status(status,STATUS_INVALID_ARGUMENT,'Independent audit supports quadrature order 4 or 5.')
    end select
  end subroutine independent_gauss_rule

  pure function cross_product(a,b) result(value)
    real(dp), intent(in) :: a(3),b(3)
    real(dp) :: value(3)
    value=[a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)]
  end function cross_product

  function integer_to_string(value) result(text)
    integer, intent(in) :: value
    character(len=32) :: text
    write(text,'(i0)') value
  end function integer_to_string
end module independent_flux_audit_m
