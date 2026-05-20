

module object_counts
    use iso_fortran_env
    implicit none
    
    
    
end module object_counts

module raw_data
    use iso_fortran_env
    implicit none
    
    
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
    integer(KIND=INT32),parameter :: binary_internal=1 ,ascii_internal=2, vtu_ascii=4, vtu_binary_appended=8
    
    
    
end module io_data
