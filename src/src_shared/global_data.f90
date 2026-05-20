

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
    
    integer(KIND=INT32) :: element_con_sum, face_con_sum
    integer(KIND=INT32), allocatable :: element_con_index(:), element_connectivity(:)
    integer(KIND=INT32), allocatable :: face_con_index(:), face_connectivity(:)
    integer(KIND=INT32), allocatable :: edge_connectivity(:)
    
    
end module raw_data



module object_relation_data
    use iso_fortran_env
    implicit none
    
    
    
end module object_relation_data


module io_data
    use iso_fortran_env
    implicit none
    
    character(:), allocatable :: raw_data_path
    
    integer(KIND=INT32) :: output_requests 
    integer(KIND=INT32),parameter :: no_output=0, binary_internal=1 ,ascii_internal=2, vtu_ascii=4, vtu_binary_appended=8
    
    
    
end module io_data
