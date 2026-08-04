

module object_counts
    use iso_fortran_env
    implicit none
    
    ! raw object counts
    integer(KIND=INT32) :: npoin, nele, nface, nedge
    
    ! abstracted object counts
        ! i_ are internal 
        ! b_ are external
    integer(KIND=INT32) :: i_npoin, b_npoin, i_nedge, b_nedge, i_nface, b_nface, i_nele, b_nele
    
end module object_counts

module raw_data
    use iso_fortran_env
    implicit none
    
    real(KIND=REAL64), allocatable :: coords(:,:)
    
    ! raw connectivity arrays 
    integer(KIND=INT32) :: c_p_sum, f_p_sum, e_p_sum
    integer(KIND=INT32), allocatable :: c_p_obj_relation_array(:), c_p_index_array(:)
    integer(KIND=INT32), allocatable :: f_p_obj_relation_array(:), f_p_index_array(:)
    integer(KIND=INT32), allocatable :: e_p_obj_relation_array(:), e_p_index_array(:)
    
    ! premapped cell to face and face to edge relations
    ! c_f :: cell index to face index
    ! f_e :: face index to edge index
    
    integer(KIND=INT32) :: c_f_sum, f_e_sum
    integer(KIND=INT32), allocatable :: c_f_obj_relation_array(:) , c_f_index_array(:)
    integer(KIND=INT32), allocatable :: f_e_obj_relation_array(:) , f_e_index_array(:)
    
end module raw_data

module object_relation_data
    use iso_fortran_env
    implicit none
    
    ! reversed connectivity arrays 
    ! p_c :: point index to cell index
    ! p_f :: point index to face index
    ! p_e :: point index to edge index
    
    integer(KIND=INT32) :: p_c_sum, p_f_sum, p_e_sum
    integer(KIND=INT32), allocatable :: p_c_obj_relation_array(:) , p_c_index_array(:)
    integer(KIND=INT32), allocatable :: p_f_obj_relation_array(:) , p_f_index_array(:)
    integer(KIND=INT32), allocatable :: p_e_obj_relation_array(:) , p_e_index_array(:)
    
    ! mapped cell to edge relation
    ! c_e :: cell index to edge index
    integer(KIND=INT32) :: c_e_sum
    integer(KIND=INT32), allocatable :: c_e_obj_relation_array(:) , c_e_index_array(:)
    
    ! reversed object mapping arrays 
    ! f_c :: face index to cell index
    ! e_f :: edge index to face index
    ! e_c :: edge index to cell index
    
    integer(KIND=INT32) :: f_c_sum, e_f_sum, e_c_sum
    integer(KIND=INT32), allocatable :: f_c_obj_relation_array(:) , f_c_index_array(:)
    integer(KIND=INT32), allocatable :: e_f_obj_relation_array(:) , e_f_index_array(:)
    integer(KIND=INT32), allocatable :: e_c_obj_relation_array(:) , e_c_index_array(:)
    
end module object_relation_data

module boundary_data
    use iso_fortran_env
    implicit none
    
    ! logical arrays of size npoin, nedge, nface, nele with boundary adjacent objects flagged true
    logical, allocatable :: p_bound_array(:) , e_bound_array(:), f_bound_array(:) , c_bound_array(:)
    
    ! indexing arrays as above, boundaries
    integer(KIND=INT32), allocatable :: p_bound_indexing_array(:) , e_bound_indexing_array(:), f_bound_indexing_array(:) , c_bound_indexing_array(:)
    
    
    ! logical arrays of size npoin, nedge, nface, nele with internal objects flagged true
    ! the logical inverse of x_bound_array
    logical, allocatable :: p_internal_array(:) , e_internal_array(:), f_internal_array(:) , c_internal_array(:)
    
    ! indexing arrays as above, internals
    integer(KIND=INT32), allocatable :: p_internal_indexing_array(:) , e_internal_indexing_array(:), f_internal_indexing_array(:) , c_internal_indexing_array(:)
    
    
    ! index arrays for moving from global object index to boundary index
    integer(KIND=INT32), allocatable :: reversed_p_bound_indexing_array(:) , reversed_e_bound_indexing_array(:), reversed_f_bound_indexing_array(:)
    
    ! bounfary condition flagging array
    integer(KIND=INT32), allocatable :: p_boundary_flags(:)
    
    ! for boundary normal vectors it makes sense to both calculate  normals for faces, and then derive point normals
    real(KIND=REAL64), allocatable :: p_normal_vectors(:,:), f_normal_vectors(:,:)
    
    ! logical array identifying feature points and edges
    logical, allocatable :: feature_points(:), feature_edges(:)
    
    ! feature point indexing array
    integer(KIND=INT32), allocatable :: n_fb_point_index(:), n_fb_point_count(:), fbp_f_index(:), fbp_f_obj_relation_array(:,:)
    
    ! number of feature points
    integer(KIND=INT32)  :: n_fb_points, fbp_f_sum
    
    ! boundary flag regions
    integer(KIND=INT32), allocatable :: f_bound_flags(:)
    
    ! the angle we class feature edges with
    real(KIND=REAL64)   :: splitting_angle = 15.0
    
end module boundary_data

module projection_data
    use iso_fortran_env
    implicit none
    
    real(KIND=REAL64), allocatable :: sn(:,:), sb(:,:), sbb(:,:), vol(:)
    
end module projection_data

module centroid_data
    use iso_fortran_env
    implicit none
    
    real(KIND=REAL64), allocatable :: e_centroid(:,:), f_centroid(:,:), c_centroid(:,:)
    
end module centroid_data

module io_data
    use iso_fortran_env
    implicit none
    
    character(:), allocatable :: raw_data_path
    
    integer(KIND=INT32) :: output_requests 
    integer(KIND=INT32),parameter :: no_output=0, binary_internal=1 ,ascii_internal=2, vtu_ascii=4, vtu_binary_appended=8
    
    integer(KIND=INT32),parameter :: real_length = 8
    
    logical :: manual_geom
    
end module io_data
