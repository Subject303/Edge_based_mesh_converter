

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
    
    call setup_configuration
    
    if (manual_geom) then
        call manual_geom_generation
    else
        call read_data
    endif
    
    call data_preprocessing
    
    call data_processing
    
    call data_outputting
    
    ! end of program
    
    
end program converter_serial
