

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


        ! once this is sorted, we now need to break out the data we need wrt to boundaries
        ! we can identify internal and external faces by whether they are connected to 1, or 2 unique cells this gives us sbb faces
        ! we can generate normal vectors with vector math in the faces, and we can ensure they're external facing by comparing with the edge barycentre
        ! we then can search for edges within boundary faces, this gives us sb edges
        ! we can flag boundaries by comparing connected boundary face normal angles
        ! we can subtract external edges from all edges for only internal edges

        call identify_external_faces(...)
        call calculate_normal_vectors(...)
        
        call identify_boundary_edges(...)
        call boundary_angle_feature_flagging(...)


    end subroutine data_preprocessing
    
end module data_preprocessing_serial
