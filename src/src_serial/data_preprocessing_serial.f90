

module data_processing_serial
    
    use preprocessor_routine_module
    use boundary_routine_module
    use object_counts
    use volume_processing
    
    implicit none
    
    contains
    
    subroutine data_preprocessing
        
        implicit none
        
        print*, 'beginning calculating object centroids'
        print*, ' '
        
        call calc_e_centroid
        call calc_f_centroid
        call calc_c_centroid
        
        print*, ' '
        print*, 'finished calculating object centroids'
        print*, ' '
        print*, 'beginning object relation preprocessing'
        print*, ' '
        
        ! from inputs we already have 
        ! e_p, f_p, c_p (connectivity arrays)
        ! c_f, f_e      (cell to face to edge relations)
        
        ! so four more object mapping arrays are needed to provide full, immediate access to all connected objects
        
        call c_e_preprocess ! conjoined from c_f and f_e
        
        print*, 'c > e relation array sums: '
        print*, 'c_e_sum'
        print*,  c_e_sum
        
        call p_c_preprocess ! don't think p_c is ever actually needed but weh
        call p_f_preprocess
        call p_e_preprocess
        
        ! I now need the last three directions,
        ! the inverse of the e f c relations
        ! e_c, e_f, e_c 

        call f_c_preprocess
        call e_f_preprocess
        call e_c_preprocess ! remember this specific one must be synced in the mpi implementation
        
        print*, 'Inverted relation array sums: '
        print*, 'p_c_sum,  p_f_sum,  p_e_sum,  f_c_sum,  e_f_sum,  e_c_sum'
        print*,  p_c_sum,  p_f_sum,  p_e_sum,  f_c_sum,  e_f_sum,  e_c_sum
        
        ! we will then have 

        ! connectivities and inverses
        ! p_e, p_f, p_c
        ! e_p, f_p, c_p

        ! object relations and inverses

        ! c_f, f_e, c_e
        ! f_c, e_f, e_c

        print*, ' '
        print*, 'finished object relation preprocessing'
        print*, ' '
        print*, 'beginning boundary identification'
        print*, ' '

        ! once this is sorted, we now need to break out the data we need wrt to boundaries
        ! we can identify internal and external faces by whether they are connected to 1, or 2 unique cells this gives us sbb faces
        ! we can generate normal vectors with vector math in the faces, and we can ensure they're external facing by comparing with the edge barycentre
        ! we then can search for edges within boundary faces, this gives us sb edges
        ! we can subtract external edges from all edges for only internal edges, giving us sn edges
        ! we can flag boundaries by comparing connected boundary face normal angles
        
        call identify_boundary_faces
        call identify_boundary_edges
        call identify_boundary_points
        
        print*, ' '
        print*, 'finished boundary identification'
        print*, ' '
        print*, 'Boundary objects total as follow: '
        print*, 'internal-boundry edges, internal-boundry faces, internal-boundry cells'
        print*, i_npoin, b_npoin, i_nedge, b_nedge, i_nface, b_nface
        print*, 'Total object sums follow: '
        print *, npoin, ' points, ', nedge, ' edges, ', nface, ' faces '
        print*, ' '
        print*, 'starting calculating boundary normal vectors'
        print*, ' '
        
        call calculate_normal_vectors
        
        print*, ' '
        print*, 'finished calculating boundary normal vectors'
        print*, ' '
        
    end subroutine data_preprocessing
    
    subroutine data_processing
        implicit none
                
        print*, ' '
        print*, 'Moving to primary volume processing. '
        print*, ' '
        
        call volume_alloc
        
        
        
    end subroutine data_processing
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
end module data_processing_serial
