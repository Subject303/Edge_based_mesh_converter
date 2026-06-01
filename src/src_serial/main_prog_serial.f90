

program converter_serial
    
    use iso_fortran_env
    
    use setup_configuration_module
    
    use read_data_serial
    use data_processing_serial
    use data_outputting_serial
    
    implicit none
    
    ! start of program
    
    call setup_configuration
    
    call read_data
    
    call data_preprocessing
    
    call data_processing
    
    call data_outputting
    
    ! end of program
    
    
end program converter_serial
