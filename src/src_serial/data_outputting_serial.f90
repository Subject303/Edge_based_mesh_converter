    
module data_outputting_serial
    use io_data
    use object_counts
    use raw_data
    use object_relation_data
    use boundary_data
    use projection_data
    use outputting_routines
    
    contains
    
    subroutine data_outputting
        
        
        implicit none
        
        
        
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
