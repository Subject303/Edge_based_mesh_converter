

program converter_serial
    use iso_fortran_env
    
    use read_data_serial
    use data_preprocessing_serial
    use data_processing_serial
    use data_outputting_serial
    
    implicit none
    
    ! start of program
    
    call read_data
    
    call data_preprocessing
    
    call data_preprocessing
    
    call data_outputting
    
    ! end of program
    
    
end program converter_serial
