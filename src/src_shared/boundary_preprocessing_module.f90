


module preprocessor_routine_module
    use iso_fortran_env
    use object_counts
    use raw_data
    use object_relation_data
    
    implicit none
    
    contains
    
    !!!!!!!!!!!!!!!!!! boundary identification routines
    
    subroutine identify_boundary_faces
        integer(KIND=INT32) :: f
        
        ! this one is easy because all faces are connected to between 1 and 2 cells
        ! all faces connected to 2 cells are internal, and all connected to 1 cell are boundaries
        
        allocate(f_bound_array(nface) ,f_internal_array(nface))
        
        f_bound_array = .false.
        f_internal_array = .false.
        
        do f=1,nface
            if (f_c_index_array(f) .eq. 1) then
            
                f_bound_array(f) = .true.
                
            elseif (f_c_index_array(f) .eq. 2) 
            ! this elseif is technically overkill and has performance overhead
            ! but belt and braces keeps my arse off the frontpages
            
                f_internal_array = .true.
                
            else
                print*, 'face id : ', f, ' , is connected to more than 2 cells or is not connected at all'
                print*, 'impressive considering faces are usually 2 dimensional'
                stop
                
            endif
        enddo
        
    end subroutine identify_boundary_faces
    
    subroutine identify_boundary_edges
        integer(KIND=INT32) :: f, fe_start, fe_end
        
        ! this one is still very simple, we're searching boundary faces, and flagging connected edges
        
        allocate(e_bound_array(nedge) ,e_internal_array(nedge))
        
        e_bound_array = .false.
        e_internal_array = .false.
        
        do f=1, nface
            
            fe_start = (1 + f_e_index_array(f-1))
            fe_end   = f_e_index_array(f)
            
            if (f_bound_array(f))     e_bound_array(   f_e_obj_relation_array(fe_start:fe_end)) = .true.
            if (f_internal_array(f))  e_internal_array(f_e_obj_relation_array(fe_start:fe_end)) = .true.
            
        enddo
        
    end subroutine identify_boundary_edges
    
    subroutine identify_boundary_points
        integer(KIND=INT32) :: f, fp_start, fp_end
        
        ! again same thing but points
        
        allocate(p_bound_array(nedge) ,p_internal_array(nedge))
        
        p_bound_array = .false.
        p_internal_array = .false.
        
        do f=1, nface
            
            fp_start = (1 + f_p_index_array(f-1))
            fp_end   = f_p_index_array(f)
            
            if (f_bound_array(f))     p_bound_array(   f_p_obj_relation_array(fp_start:fp_end)) = .true.
            if (f_internal_array(f))  p_internal_array(f_p_obj_relation_array(fp_start:fp_end)) = .true.
            
        enddo
    end subroutine identify_boundary_points
    
    !!!!!!!!!!!!!!!!!! boundary flagging & normal vector routines
    
    subroutine calculate_normal_vectors
    
    end subroutine calculate_normal_vectors
    
    subroutine boundary_angle_feature_flagging
    
    end subroutine boundary_angle_feature_flagging
    
    ! end contains    
    
end module preprocessor_routine_module
