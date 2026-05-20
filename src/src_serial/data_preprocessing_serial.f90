

module data_preprocessing_serial
    
    use preprocessor_routine_module
    
    implicit none
    
    contains
    
    subroutine data_preprocessing
    
        call p_c_preprocess
        call p_f_preprocess
        call p_e_preprocess
        
    end subroutine data_preprocessing
    
end module data_preprocessing_serial
