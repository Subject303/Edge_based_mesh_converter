

program converter_serial
    use iso_fortran_env
    
    use ...
    use ...
    use ...
    use ...
    
    implicit none
    
    ! start of program
    
    
    call read_data_serial
    
    
    call data_preprocessing
    
    ! start of data processing
    
    call internal_edge_processing
    call external_edge_processing
    call external_face_processing
    
    ! start of outputting
    
    select case (output_requests)
        
        case(output_requests .gt. binary_internal)
            call output_binary_internal
            
        case(output_requests .gt. ascii_internal)
            call output_ascii_internal
            
        case(output_requests .gt. vtu_ascii)
            call output_vtu_ascii
            
        case(output_requests .gt. vtu_binary_appended)
            call output_vtu_binary_appended
            
    end select
    
    ! end of program
    
    contains 
    
        subroutine data_preprocessing
        
            call ...
            call ...
            call ...
            
        end subroutine data_preprocessing

    !end contains
    
end program converter_serial
