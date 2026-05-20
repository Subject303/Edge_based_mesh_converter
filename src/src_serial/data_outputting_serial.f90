    
module data_outputting_serial
    
    
    contains
    
    subroutine data_outputting
        
        
        
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
        
    end subroutine data_outputting
    
end module data_outputting_serial
