    
module data_outputting_serial
    
    
    contains
    
    subroutine data_outputting
        
        
        
        
        if ((output_requests - binary_internal) .ge. 0) then
                call output_binary_internal
        endif
        
        if ((output_requests - ascii_internal) .ge. 0) then
                call output_ascii_internal
        endif
        
        if ((output_requests - vtu_ascii) .ge. 0) then
                call output_vtu_ascii
        endif
        
        if ((output_requests - vtu_binary_appended) .ge. 0) then
                call output_vtu_binary_appended
        endif
        
        
    end subroutine data_outputting
    
end module data_outputting_serial
