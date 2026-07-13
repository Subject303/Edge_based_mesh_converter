module manual_geom_generation_module
    
    use iso_fortran_env
    
    use setup_configuration_module
    
    use read_data_serial
    use data_processing_serial
    use data_outputting_serial
    
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
    
    
end module manual_geom_generation_module
