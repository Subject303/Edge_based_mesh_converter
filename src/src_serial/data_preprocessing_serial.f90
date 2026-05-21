

module data_preprocessing_serial
    
    use preprocessor_routine_module
    
    implicit none
    
    contains
    
    subroutine data_preprocessing
        
        ! from inputs we already have 
        ! e_p, f_p, c_p (connectivity arrays)
        ! c_f, f_e      (cell to face to edge relations)
        
        ! so four more object mapping arrays are needed to provide full, immediate access to all connected objects
        
        call c_e_preprocess
        
        call p_c_preprocess
        call p_f_preprocess
        call p_e_preprocess

        ! I now need the last three directions,
        ! the inverse of the e f c relations
        ! e_c, e_f, e_c 

        call e_c_preprocess
        call e_f_preprocess
        call e_c_preprocess ! remember this specific one must be synced in the mpi implementation

        ! we will then have 

        ! connectivities and inverses
        ! p_e, p_f, p_c
        ! e_p, f_p, c_p

        ! object relations and inverses

        ! c_f, f_e, c_e
        ! f_c, e_f, e_c

    end subroutine data_preprocessing
    
end module data_preprocessing_serial
