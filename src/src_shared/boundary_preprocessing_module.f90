


module preprocessor_routine_module
    use iso_fortran_env
    use object_counts
    use raw_data
    use object_relation_data
    
    implicit none
    
    contains
    
    !!!!!!!!!!!!!!!!!! boundary identification routines
    
    subroutine identify_boundary_faces
        integer(KIND=INT32) :: f, i, b
        
        ! this one is easy because all faces are connected to between 1 and 2 cells
        ! all faces connected to 2 cells are internal, and all connected to 1 cell are boundaries
        
        allocate(f_bound_array(nface) ,f_internal_array(nface), f_bound_indexing_array(nface), f_internal_indexing_array(nface))
        
        f_bound_array = .false.
        f_internal_array = .false.
        f_bound_indexing_array    = (/1:nface/)
        allocate(f_internal_indexing_array,source=f_bound_indexing_array)
        
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
        
        f_bound_indexing_array   =pack(f_bound_indexing_array   ,f_bound_array)
        f_internal_indexing_array=pack(f_internal_indexing_array,f_internal_array)
        b_nface = size(f_bound_indexing_array)
        i_nface = size(f_internal_indexing_array)
        
    end subroutine identify_boundary_faces
    
    subroutine identify_boundary_edges
        integer(KIND=INT32) :: f, fe_start, fe_end
        
        ! this one is still very simple, we're searching boundary faces, and flagging connected edges
        
        allocate(e_bound_array(nedge) ,e_internal_array(nedge), e_bound_indexing_array(nedge))
        
        e_bound_array = .false.
        e_internal_array = .false.
        e_bound_indexing_array    = (/1:nedge/)
        allocate(e_internal_indexing_array,source=e_bound_indexing_array)
        
        do f=1, nface
            
            fe_start = (1 + f_e_index_array(f-1))
            fe_end   = f_e_index_array(f)
            
            if (f_bound_array(f))then
                e_bound_array(   f_e_obj_relation_array(fe_start:fe_end)) = .true.
            elseif (f_internal_array(f))
                e_internal_array(f_e_obj_relation_array(fe_start:fe_end)) = .true.
            else
                print*, 'face id : ', f, ' , is not flagged internal or external'
            endif
            
        enddo
        
        e_bound_indexing_array   =pack(e_bound_indexing_array   ,e_bound_array)
        e_internal_indexing_array=pack(e_internal_indexing_array,e_internal_array)
        b_nedge = size(e_bound_indexing_array)
        i_nedge = size(e_internal_indexing_array)
        
    end subroutine identify_boundary_edges
    
    subroutine identify_boundary_points
        integer(KIND=INT32) :: f, fp_start, fp_end
        
        ! again same thing but points
        
        allocate(p_bound_array(npoin) ,p_internal_array(npoin), e_bound_indexing_array(npoin))
        
        p_bound_array = .false.
        p_internal_array = .false.
        p_bound_indexing_array    = (/1:npoin/)
        allocate(p_internal_indexing_array,source=p_bound_indexing_array)
        
        do f=1, nface
            
            fp_start = (1 + f_p_index_array(f-1))
            fp_end   = f_p_index_array(f)
            
            if (f_bound_array(f))then
                p_bound_array(   f_p_obj_relation_array(fp_start:fp_end)) = .true.
            elseif (f_internal_array(f))
                p_internal_array(f_p_obj_relation_array(fp_start:fp_end)) = .true.
            else
                print*, 'face id : ', f, ' , is not flagged internal or external'
            endif
            
        enddo
        
        p_bound_indexing_array   =pack(p_bound_indexing_array   ,p_bound_array)
        p_internal_indexing_array=pack(p_internal_indexing_array,p_internal_array)
        b_npoin = size(p_bound_indexing_array)
        i_npoin = size(p_internal_indexing_array)
        
    end subroutine identify_boundary_points
    
    !!!!!!!!!!!!!!!!!! boundary flagging & normal vector routines
    
    subroutine calculate_normal_vectors
        integer(KIND=INT32) :: f, bf, fp_start, p1, p2, p3, magnitude, c
        real(KIND=REAL32)   :: v1(3), v2(3)
        
        allocate(f_normal_vectors(b_nface,3))
        
        ! so first we loop boundary faces and calc normal vectors from the normal to the cross plane p1 > p2 and p2 > p3
        ! we can assume all faces have at least 3 points (as 2 pointed objects are edges and 1 pointed edges are points) so we only need the start index 
        
        do bf=1,b_nface
            
            f = f_bound_indexing_array(bf)
            
            fp_start = (1 + f_p_index_array(f-1))
            
            p1 = f_p_obj_relation_array(fp_start)
            p2 = f_p_obj_relation_array(fp_start+1)
            p3 = f_p_obj_relation_array(fp_start+2)
            
            v1(:) = coords(p1,:) - coords(p2,:)
            v2(:) = coords(p1,:) - coords(p3,:)
            
            f_normal_vectors(bf,1) = v1(2)*v2(3) - v1(3)*v2(2)
            f_normal_vectors(bf,2) = v1(3)*v2(1) - v1(1)*v2(3)
            f_normal_vectors(bf,3) = v1(1)*v2(2) - v1(2)*v2(1)
            
            magnitude = f_normal_vectors(bf,1)*f_normal_vectors(bf,1) + f_normal_vectors(bf,2)*f_normal_vectors(bf,2) + f_normal_vectors(bf,3)*f_normal_vectors(bf,3)
            
            f_normal_vectors(bf,:)  = f_normal_vectors(bf,:) / magnitude
            
            c = f_c_obj_relation_array(f_c_index_array(f))
            
            ! I want barycentres already
            
        endif
        
    end subroutine calculate_normal_vectors
    
    subroutine boundary_angle_feature_flagging
    
    end subroutine boundary_angle_feature_flagging
    
    ! end contains    
    
end module preprocessor_routine_module
