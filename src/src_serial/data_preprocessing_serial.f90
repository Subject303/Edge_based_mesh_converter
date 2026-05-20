

module data_preprocessing_serial
    
    use preprocessor_routine_module
    
    implicit none
    
    contains
    
    subroutine data_preprocessing
        
        ! from inputs we already have 
        ! e_p, f_p, c_p (connectivity arrays)
        ! c_e, f_e      (cell to face to edge relations)
        
        ! so four more object mapping arrays are needed to provide full, immediate access to all connected objects
        
        call c_e_preprocess
        
        call p_c_preprocess
        call p_f_preprocess
        call p_e_preprocess
        
    end subroutine data_preprocessing
    
end module data_preprocessing_serial
