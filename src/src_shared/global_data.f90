

module object_counts
    use iso_fortran_env
    implicit none
    
    ! raw object counts
    integer(KIND=INT32) :: npoin, nele, nface, nedge
    
    ! abstracted object counts
        ! i_ are internal 
        ! e_ are external
    integer(KIND=INT32) :: i_nedge, e_nedge, i_nface, e_nface, i_npoin, e_npoin
    
end module object_counts

module raw_data
    use iso_fortran_env
    implicit none
    
    real(KIND=REAL64), allocatable :: coords(:,:)
    
    ! raw connectivity arrays 
    integer(KIND=INT32) :: c_p_sum, f_p_sum
    integer(KIND=INT32), allocatable :: c_p_obj_relation_array(:), c_p_index_array(:)
    integer(KIND=INT32), allocatable :: f_p_obj_relation_array(:), f_p_index_array(:), 
    integer(KIND=INT32), allocatable :: e_p_obj_relation_array(:)
    
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
    
end module object_relation_data


module io_data
    use iso_fortran_env
    implicit none
    
    character(:), allocatable :: raw_data_path
    
    integer(KIND=INT32) :: output_requests 
    integer(KIND=INT32),parameter :: no_output=0, binary_internal=1 ,ascii_internal=2, vtu_ascii=4, vtu_binary_appended=8
    
    
    
end module io_data
