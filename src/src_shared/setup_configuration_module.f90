

module setup_configuration_module
    
    use io_data
    
    implicit none
    
    contains
    
    subroutine setup_configuration
        implicit none
        
        character(len = :), allocatable :: current_line, &
        & prerocessed_data_file_str, binary_internal_str, ascii_internal_str, vtu_ascii_str, vtu_binary_appended_str, temp_str
        integer(KIND=int32) :: indexer
        !character :: current_line*(*)
        
        output_requests = 0
        
        prerocessed_data_file_str = 'prerocessed_data_file='
        binary_internal_str = 'binary_internal='
        ascii_internal_str = 'ascii_internal='
        vtu_ascii_str = 'vtu_ascii='
        vtu_binary_appended_str = 'vtu_binary_appended='
        
        open(10,file='config.cfg',access='sequential',action='read',status='old')
        
        do 
            read(10,*) current_line
            
            current_line = trim(adjustL(current_line))
            
            if (len(current_line).gt.0) then
                indexer = 1+index(current_line,'=')
                temp_str = current_line(1:indexer)
                
                select case(temp_str)
                
                    case('prerocessed_data_file=')
                        raw_data_path = current_line(len('prerocessed_data_file='):)
                        
                    case('binary_internal=')
                        indexer = index(current_line,'true')
                        if (indexer.ne.0) output_requests = output_requests + binary_internal
                        
                    case('ascii_internal=')
                        indexer = index(current_line,'true')
                        if (indexer.ne.0) output_requests = output_requests + ascii_internal
                        
                    case('vtu_ascii=')
                        indexer = index(current_line,'true')
                        if (indexer.ne.0) output_requests = output_requests + vtu_ascii
                        
                    case('vtu_binary_appended=')
                        indexer = index(current_line,'true')
                        if (indexer.ne.0) output_requests = output_requests + vtu_binary_appended
                        
                        ! last line in config file
                        
                        exit
                        
                end select
                
            endif
        enddo
        
        if (.not.allocated(raw_data_path)) print*, 'no raw data file found', raw_data_path
        
        if (output_requests .eq. no_output) print*, 'no output formats specified'
        
        close(10)
        
    end subroutine setup_configuration
    
    ! end contains
    
end module setup_configuration_module
