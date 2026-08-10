

module setup_configuration_module
    
    use io_data
    use boundary_data
    
    implicit none
    
    contains
    
    subroutine setup_configuration
        implicit none
        
        character(len = :), allocatable :: current_line, temp_str
        integer(KIND=int32) :: indexer
        logical :: file_end
        
        output_requests = 0
        manual_geom = .false.
        
        open(10,file='./config.cfg',access='sequential',action='read',status='old')
        
        file_end = .false.
        
        do while(.not.file_end)
            allocate(character(len=512) :: current_line)
            
            read(10,'(A)') current_line
            
            current_line = trim(adjustL(current_line))
            
            if (len(current_line).gt.0) then
                indexer = index(current_line,'=')
                temp_str = current_line(1:indexer)
                
                select case(temp_str)
                    case('feature_angle=')
                        read(current_line(indexer+1:), '(f12.9)') splitting_angle
                        
                    case('manual_geom=')
                        indexer = index(current_line,'true')
                        if (indexer.ne.0) then
                            manual_geom = .true.
                        endif
                        
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
                        file_end = .true.
                        
                end select
            endif
            deallocate(current_line)
        enddo
        
        call GET_COMMAND_ARGUMENT (1 , length=indexer)
        
        if (indexer .eq. 0) then
            print*, 'no raw data file found ', raw_data_path, ' somthing is broken '
            print*, 'stopping'
            stop
        endif
        
        allocate(character(indexer)::raw_data_path)
        
        call GET_COMMAND_ARGUMENT (1 , value=raw_data_path)
        
        ! 28 chars in ./preprocessed_mesh_folder/
        ! 23 in .preprocessed_mesh_file
        binary_intern_data_path = raw_data_path(28:(len(raw_data_path)-23))//".internalmesh"
        print*, binary_intern_data_path
        
        
        if (output_requests .eq. no_output) print*, 'no output formats specified'
        
        close(10)
        
    end subroutine setup_configuration
    
    ! end contains
    
end module setup_configuration_module
