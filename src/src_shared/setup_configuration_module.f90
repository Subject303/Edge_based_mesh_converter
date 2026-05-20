

module setup_configuration_module
    
    use iso_fortran_env
    
    use io_data
    
    implicit none
    
    contains
    
    subroutine setup_configuration
        
        
        character(:), allocatable :: current_line&
        & prerocessed_data_file_str, binary_internal_str, ascii_internal_str, vtu_ascii_str, vtu_binary_appended_str
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
            
            if (len(current_line)) then
                select case(current_line)
                
                    case(current_line(1:len(prerocessed_data_file_str)) .eq. prerocessed_data_file_str)
                        raw_data_path = current_line(len(prerocessed_data_file_str):)
                        
                    case(current_line(1:len(binary_internal_str)) .eq. binary_internal_str)
                        if ( trim(adjustL(current_line((1+len(binary_internal_str)):len(current_line)))) .eq. 'true' ) &
                        & output_requests = output_requests + binary_internal
                        
                    case(current_line(1:len(ascii_internal_str)) .eq. ascii_internal_str)
                        if ( trim(adjustL(current_line((1+len(binary_internal_str)):len(current_line)))) .eq. 'true' ) &
                        & output_requests = output_requests + ascii_internal
                        
                    case(current_line(1:len(vtu_ascii_str)) .eq. vtu_ascii_str)
                        if ( trim(adjustL(current_line((1+len(binary_internal_str)):len(current_line)))) .eq. 'true' ) &
                        & output_requests = output_requests + vtu_ascii
                        
                    case(current_line(1:len(vtu_binary_appended_str)) .eq. vtu_binary_appended_str)
                        if ( trim(adjustL(current_line((1+len(binary_internal_str)):len(current_line)))) .eq. 'true' ) &
                        & output_requests = output_requests + vtu_binary_appended
                        
                end select
            endif
        enddo
        
        
        
    end subroutine setup_configuration
    
    ! end contains
    
end module setup_configuration_module
