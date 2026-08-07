

program converter_serial
    
    use iso_fortran_env
    
    use setup_configuration_module
    
    use io_data
    
    use read_data_serial
    use data_processing_serial
    use data_outputting_serial
    
    use manual_geom_generation_module
    
    implicit none
    
    ! start of program
    
    print*, 'STARTING MAIN PROCESSOR SCRIPT'
    
    call setup_configuration
    
    print*, 'FINISHED READING CONFIG FILE'
    
    if (manual_geom) then
        call manual_geom_generation
    else
        call read_data
        
        print*, 'FINISHED READING PREPROCESSED MESH FILE'
        
    endif
    
    call data_preprocessing
    
    print*, 'FINISHED CONNECTIVITY PREPROCCESSING'
    
    call data_processing
    
    print*, 'FINISHED DATA PROCESSING'
    
    call data_outputting
    
    print*, 'FINISHED DATA OUTPUT'
    
    ! end of program
    
    print*, 'EXITTING MAIN PROCESSOR SCRIPT'
    
end program converter_serial
