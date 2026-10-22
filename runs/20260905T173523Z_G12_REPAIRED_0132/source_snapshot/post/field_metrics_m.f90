module field_metrics_m
  use kinds_m, only : dp
  use status_m, only : SolverStatusType, clear_status, set_status, &
    STATUS_INVALID_ARGUMENT, STATUS_GEOMETRY_ERROR
  use field_layout_m, only : NDOF_NODE,F_U1,F_U2,F_U3,local_dof
  use mesh_types_m, only : MeshType, HEX20_NNODE
  use model_types_m, only : MaterialType, StateType, LoadStateType
  use dof_map_m, only : HEX20_NDOF, build_element_dof_map
  use hex20_shape_m, only : hex20_geometry
  use element_kernel_m, only : build_operators, unpack_point_state
  use material_point_m, only : MATERIAL_INPUT_SIZE, MATERIAL_RESPONSE_SIZE, &
    MaterialPointStateType, MaterialPointResponseType, MaterialPointTangentType, &
    evaluate_material_point
  implicit none
  private

  integer, parameter :: NQXY=5,NQZ=3,NQ=NQXY*NQXY*NQZ

  type, public :: HoleFieldMetricsType
    real(dp) :: volume=0.0_dp
    real(dp) :: annulus_volume=0.0_dp
    real(dp) :: total_joule_power=0.0_dp
    real(dp) :: annulus_joule_mean=0.0_dp
    real(dp) :: annulus_joule_p95=0.0_dp
    real(dp) :: annulus_electric_field_mean=0.0_dp
    real(dp) :: annulus_n_relative_mean=0.0_dp
    real(dp) :: annulus_p_relative_mean=0.0_dp
    real(dp) :: annulus_temperature_rise_mean=0.0_dp
    real(dp) :: annulus_displacement_mean=0.0_dp
  end type HoleFieldMetricsType

  public :: evaluate_hole_field_metrics

contains

  subroutine evaluate_hole_field_metrics(mesh,materials,state,load,hole_radius, &
                                         annulus_width,metrics,status,annulus_subdivisions)
    type(MeshType), intent(in) :: mesh
    type(MaterialType), intent(in) :: materials(:)
    type(StateType), intent(in) :: state
    type(LoadStateType), intent(in) :: load
    real(dp), intent(in) :: hole_radius,annulus_width
    type(HoleFieldMetricsType), intent(out) :: metrics
    type(SolverStatusType), intent(out) :: status
    integer, intent(in), optional :: annulus_subdivisions
    real(dp) :: points(3,NQ),weights(NQ),coordinates_e(3,HEX20_NNODE)
    real(dp) :: state_e(HEX20_NDOF),shape(HEX20_NNODE)
    real(dp) :: derivative_global(3,HEX20_NNODE),determinant,integration_weight
    real(dp) :: state_operator(MATERIAL_INPUT_SIZE,HEX20_NDOF)
    real(dp) :: test_operator(HEX20_NDOF,MATERIAL_RESPONSE_SIZE)
    real(dp) :: independent_state(MATERIAL_INPUT_SIZE),position(3),center(2),radius
    real(dp) :: displacement(3)
    real(dp), allocatable :: annulus_joule(:),annulus_weights(:)
    integer :: element,point,material_id,sample_count,node
    integer :: subdivision_count,subdivision_r,subdivision_s
    real(dp) :: mapped_r,mapped_s,subdivision_center_r,subdivision_center_s
    real(dp) :: subdivision_scale
    integer :: element_dof(HEX20_NDOF)
    type(MaterialPointStateType) :: point_state
    type(MaterialPointResponseType) :: point_response
    type(MaterialPointTangentType) :: point_tangent
    type(SolverStatusType) :: local_status

    call clear_status(status)
    metrics=HoleFieldMetricsType()
    subdivision_count=1
    if (present(annulus_subdivisions)) subdivision_count=annulus_subdivisions
    if (hole_radius <= 0.0_dp .or. annulus_width <= 0.0_dp .or. &
        mesh%number_of_elements() < 1 .or. size(materials) < 1 .or. &
        .not. allocated(state%x) .or. subdivision_count < 1) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'Hole field metric input is invalid.')
      return
    end if
    if (size(state%x) /= NDOF_NODE*mesh%number_of_nodes()) then
      call set_status(status,STATUS_INVALID_ARGUMENT,'Hole field metric state size is invalid.')
      return
    end if
    center=0.5_dp*(minval(mesh%coordinates(1:2,:),dim=2)+ &
      maxval(mesh%coordinates(1:2,:),dim=2))
    allocate(annulus_joule(NQ*mesh%number_of_elements()*subdivision_count**2), &
      annulus_weights(NQ*mesh%number_of_elements()*subdivision_count**2))
    annulus_joule=0.0_dp; annulus_weights=0.0_dp; sample_count=0
    call gauss_rule_5x5x3(points,weights)

    do element=1,mesh%number_of_elements()
      material_id=mesh%material_id(element)
      if (material_id < 1 .or. material_id > size(materials)) then
        call set_status(status,STATUS_INVALID_ARGUMENT, &
          'Hole field metric encountered an invalid material identifier.')
        return
      end if
      coordinates_e=mesh%coordinates(:,mesh%connectivity(:,element))
      call build_element_dof_map(mesh%connectivity(:,element),element_dof,local_status)
      if (.not. local_status%is_ok()) then
        call set_status(status,local_status%code,trim(local_status%message)); return
      end if
      state_e=state%x(element_dof)
      subdivision_scale=1.0_dp/real(subdivision_count,dp)
      do subdivision_r=1,subdivision_count
        subdivision_center_r=-1.0_dp+(2.0_dp*real(subdivision_r,dp)-1.0_dp)* &
          subdivision_scale
        do subdivision_s=1,subdivision_count
          subdivision_center_s=-1.0_dp+(2.0_dp*real(subdivision_s,dp)-1.0_dp)* &
            subdivision_scale
          do point=1,NQ
            mapped_r=subdivision_center_r+subdivision_scale*points(1,point)
            mapped_s=subdivision_center_s+subdivision_scale*points(2,point)
            call hex20_geometry(coordinates_e,mapped_r,mapped_s, &
              points(3,point),shape,derivative_global,determinant,local_status)
            if (.not. local_status%is_ok()) then
              call set_status(status,STATUS_GEOMETRY_ERROR, &
                'Independent post-processing encountered invalid element geometry.')
              return
            end if
            integration_weight=weights(point)*determinant*subdivision_scale**2
            call build_operators(shape,derivative_global,state_operator,test_operator)
            independent_state=matmul(state_operator,state_e)
            call unpack_point_state(materials(material_id),independent_state,point_state)
            call evaluate_material_point(materials(material_id),point_state,load, &
              point_response,point_tangent,local_status)
            if (.not. local_status%is_ok()) then
              call set_status(status,local_status%code,trim(local_status%message)); return
            end if
            position=matmul(coordinates_e,shape)
            displacement=0.0_dp
            do node=1,HEX20_NNODE
              displacement(1)=displacement(1)+shape(node)*state_e(local_dof(node,F_U1))
              displacement(2)=displacement(2)+shape(node)*state_e(local_dof(node,F_U2))
              displacement(3)=displacement(3)+shape(node)*state_e(local_dof(node,F_U3))
            end do
            radius=hypot(position(1)-center(1),position(2)-center(2))
            metrics%volume=metrics%volume+integration_weight
            metrics%total_joule_power=metrics%total_joule_power+ &
              integration_weight*point_response%joule_heat
            if (radius >= hole_radius .and. radius <= hole_radius+annulus_width) then
              sample_count=sample_count+1
              annulus_joule(sample_count)=point_response%joule_heat
              annulus_weights(sample_count)=integration_weight
              metrics%annulus_volume=metrics%annulus_volume+integration_weight
              metrics%annulus_joule_mean=metrics%annulus_joule_mean+ &
                integration_weight*point_response%joule_heat
              metrics%annulus_electric_field_mean=metrics%annulus_electric_field_mean+ &
                integration_weight*sqrt(sum(point_state%electric_field**2))
              metrics%annulus_n_relative_mean=metrics%annulus_n_relative_mean+ &
                integration_weight*abs(point_response%n/materials(material_id)%reference_n-1.0_dp)
              metrics%annulus_p_relative_mean=metrics%annulus_p_relative_mean+ &
                integration_weight*abs(point_response%p/materials(material_id)%reference_p-1.0_dp)
              metrics%annulus_temperature_rise_mean=metrics%annulus_temperature_rise_mean+ &
                integration_weight*abs(point_state%temperature- &
                  materials(material_id)%reference_temperature)
              metrics%annulus_displacement_mean=metrics%annulus_displacement_mean+ &
                integration_weight*sqrt(sum(displacement**2))
            end if
          end do
        end do
      end do
    end do
    if (sample_count < 1 .or. metrics%annulus_volume <= tiny(1.0_dp)) then
      call set_status(status,STATUS_INVALID_ARGUMENT, &
        'No independent quadrature sample lies in the frozen hole annulus.')
      return
    end if
    metrics%annulus_joule_mean=metrics%annulus_joule_mean/metrics%annulus_volume
    metrics%annulus_electric_field_mean= &
      metrics%annulus_electric_field_mean/metrics%annulus_volume
    metrics%annulus_n_relative_mean= &
      metrics%annulus_n_relative_mean/metrics%annulus_volume
    metrics%annulus_p_relative_mean= &
      metrics%annulus_p_relative_mean/metrics%annulus_volume
    metrics%annulus_temperature_rise_mean= &
      metrics%annulus_temperature_rise_mean/metrics%annulus_volume
    metrics%annulus_displacement_mean= &
      metrics%annulus_displacement_mean/metrics%annulus_volume
    metrics%annulus_joule_p95=weighted_percentile(annulus_joule(1:sample_count), &
      annulus_weights(1:sample_count),0.95_dp)
  end subroutine evaluate_hole_field_metrics

  subroutine gauss_rule_5x5x3(points,weights)
    real(dp), intent(out) :: points(3,NQ),weights(NQ)
    real(dp), parameter :: xy(5)=[-0.9061798459386640_dp,-0.5384693101056831_dp, &
      0.0_dp,0.5384693101056831_dp,0.9061798459386640_dp]
    real(dp), parameter :: wx(5)=[0.2369268850561891_dp,0.4786286704993665_dp, &
      0.5688888888888889_dp,0.4786286704993665_dp,0.2369268850561891_dp]
    real(dp), parameter :: z(3)=[-sqrt(3.0_dp/5.0_dp),0.0_dp,sqrt(3.0_dp/5.0_dp)]
    real(dp), parameter :: wz(3)=[5.0_dp/9.0_dp,8.0_dp/9.0_dp,5.0_dp/9.0_dp]
    integer :: i,j,k,point_index

    point_index=0
    do i=1,5
      do j=1,5
        do k=1,3
          point_index=point_index+1
          points(1,point_index)=xy(i)
          points(2,point_index)=xy(j)
          points(3,point_index)=z(k)
          weights(point_index)=wx(i)*wx(j)*wz(k)
        end do
      end do
    end do
    if (point_index /= NQ) error stop 'Internal 5x5x3 quadrature count mismatch.'
  end subroutine gauss_rule_5x5x3

  function weighted_percentile(values,input_weights,fraction) result(value)
    real(dp), intent(in) :: values(:),input_weights(:),fraction
    real(dp) :: value,cumulative,target
    real(dp), allocatable :: sorted_values(:),sorted_weights(:)
    integer :: index

    allocate(sorted_values(size(values)),sorted_weights(size(input_weights)))
    sorted_values=values; sorted_weights=input_weights
    call quicksort_pairs(sorted_values,sorted_weights,1,size(sorted_values))
    target=fraction*sum(sorted_weights)
    cumulative=0.0_dp; value=sorted_values(size(sorted_values))
    do index=1,size(sorted_values)
      cumulative=cumulative+sorted_weights(index)
      if (cumulative >= target) then
        value=sorted_values(index)
        exit
      end if
    end do
  end function weighted_percentile

  recursive subroutine quicksort_pairs(values,weights,left,right)
    real(dp), intent(inout) :: values(:),weights(:)
    integer, intent(in) :: left,right
    integer :: i,j
    real(dp) :: pivot,temporary

    if (left >= right) return
    i=left; j=right; pivot=values((left+right)/2)
    do
      do
        if (i > right) exit
        if (values(i) >= pivot) exit
        i=i+1
      end do
      do
        if (j < left) exit
        if (values(j) <= pivot) exit
        j=j-1
      end do
      if (i > j) exit
      temporary=values(i); values(i)=values(j); values(j)=temporary
      temporary=weights(i); weights(i)=weights(j); weights(j)=temporary
      i=i+1; j=j-1
    end do
    if (left < j) call quicksort_pairs(values,weights,left,j)
    if (i < right) call quicksort_pairs(values,weights,i,right)
  end subroutine quicksort_pairs

end module field_metrics_m
