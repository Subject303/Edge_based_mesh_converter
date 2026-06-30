


module boundary_routine_module
    use iso_fortran_env
    use object_counts
    use raw_data
    use object_relation_data
    use boundary_data
    use quicksort_module
    use utils, ONLY : alignment
    implicit none
    
    contains
    
    !!!!!!!!!!!!!!!!!! boundary identification routines
    
    subroutine identify_boundary_faces
        implicit none
        integer(KIND=INT32) :: f, i, b, nfaces
        
        ! this one is easy because all faces are connected to between 1 and 2 cells
        ! all faces connected to 2 cells are internal, and all connected to 1 cell are boundaries
        
        allocate(f_bound_array(nface) ,f_internal_array(nface), f_bound_indexing_array(nface), reversed_f_bound_indexing_array(nface))
        
        f_bound_array = .false.
        f_internal_array = .false.
        f_bound_indexing_array    = (/(f, f=1,nface)/)
        allocate(f_internal_indexing_array,source=f_bound_indexing_array)
        
        do f=1,nface
        
            nfaces = f_c_index_array(f) - f_c_index_array(f-1)
            
            if (nfaces .eq. 1) then
            
                f_bound_array(f) = .true.
                
            elseif (nfaces .eq. 2) then
            ! this elseif is technically overkill and has performance overhead
            ! but belt and braces keeps my arse off the frontpages
            
                f_internal_array(f) = .true.
                
            else
                print*, 'ERROR ERROR'
                print*, 'face id : ', f, ' , is connected to more than 2 cells or is not connected at all'
                print*, 'impressive considering faces are usually 2 dimensional'
                stop
                
            endif
        enddo
        
        reversed_f_bound_indexing_array = 0
        do i=1,nface
            if (f_internal_array(i)) reversed_f_bound_indexing_array(i) = reversed_f_bound_indexing_array(i) + 1
        enddo
        do i=2,nface
            reversed_f_bound_indexing_array(i) = reversed_f_bound_indexing_array(i-1) + reversed_f_bound_indexing_array(i)
        enddo
        reversed_f_bound_indexing_array = f_bound_indexing_array - reversed_f_bound_indexing_array
        
        f_bound_indexing_array   =pack(f_bound_indexing_array   ,f_bound_array)
        f_internal_indexing_array=pack(f_internal_indexing_array,f_internal_array)
        b_nface = size(f_bound_indexing_array)
        i_nface = size(f_internal_indexing_array)
        
        if (b_nface+i_nface .ne. nface)then
            print*, 'ERROR ERROR'
            print*, 'sum of internal and external faces is different to the total count of faces.'
            stop
        endif
        
    end subroutine identify_boundary_faces
    
    subroutine identify_boundary_edges
        implicit none
        integer(KIND=INT32) :: f, fe_start, fe_end, i
        
        
        ! this one is still very simple, we're searching boundary faces, and flagging connected edges
        
        allocate(e_bound_array(nedge) ,e_internal_array(nedge), e_bound_indexing_array(nedge), reversed_e_bound_indexing_array(nedge))
        
        e_bound_array = .false.
        e_internal_array = .false.
        e_bound_indexing_array    = (/(f, f=1,nedge)/)
        allocate(e_internal_indexing_array,source=e_bound_indexing_array)
        
        do f=1, nface
            
            fe_start = (1 + f_e_index_array(f-1))
            fe_end   = f_e_index_array(f)
            
            if (f_bound_array(f))then
                
                e_bound_array(   f_e_obj_relation_array(fe_start:fe_end)) = .true.
                
            elseif (f_internal_array(f)) then
                ! purely a convinient place to track if the internal array has been flagged right
            else
                print*, 'face id : ', f, ' , is not flagged internal or external'
            endif
            
        enddo
        
        e_internal_array = .not. e_bound_array
        
        reversed_e_bound_indexing_array = 0
        do i=1,nedge
            if (e_internal_array(i)) reversed_e_bound_indexing_array(i) = reversed_e_bound_indexing_array(i) + 1
        enddo
        do i=2,nedge
            reversed_e_bound_indexing_array(i) = reversed_e_bound_indexing_array(i-1) + reversed_e_bound_indexing_array(i)
        enddo
        reversed_e_bound_indexing_array = e_bound_indexing_array - reversed_e_bound_indexing_array
        
        e_bound_indexing_array   =pack(e_bound_indexing_array   ,e_bound_array)
        e_internal_indexing_array=pack(e_internal_indexing_array,e_internal_array)
        b_nedge = size(e_bound_indexing_array)
        i_nedge = size(e_internal_indexing_array)
        
        if (b_nedge+i_nedge .ne. nedge)then
            print*, 'ERROR ERROR'
            print*, 'sum of internal and external edges is different to the total count of edges.'
            stop
        endif
        
    end subroutine identify_boundary_edges
    
    subroutine identify_boundary_points
        implicit none
        integer(KIND=INT32) :: f, fp_start, fp_end, i
        
        
        ! again same thing but points
        
        allocate(p_bound_array(npoin) ,p_internal_array(npoin), p_bound_indexing_array(npoin),reversed_p_bound_indexing_array(npoin))
        
        p_bound_array = .false.
        p_internal_array = .false.
        p_bound_indexing_array    = (/(f, f=1,npoin)/)
        allocate(p_internal_indexing_array,source=p_bound_indexing_array)
        
        do f=1, nface
            
            fp_start = (1 + f_p_index_array(f-1))
            fp_end   = f_p_index_array(f)
            
            if (f_bound_array(f))then
                p_bound_array(   f_p_obj_relation_array(fp_start:fp_end)) = .true.
            elseif (f_internal_array(f)) then
                ! again purely a convinient place to track if the internal array has been flagged right
            else
                print*, 'face id : ', f, ' , is not flagged internal or external'
            endif
            
        enddo
        
        p_internal_array = .not. p_bound_array
        
        reversed_p_bound_indexing_array = 0
        do i=1,npoin
            if (p_internal_array(i)) reversed_p_bound_indexing_array(i) = reversed_p_bound_indexing_array(i) + 1
        enddo
        do i=2,npoin
            reversed_p_bound_indexing_array(i) = reversed_p_bound_indexing_array(i-1) + reversed_p_bound_indexing_array(i)
        enddo
        reversed_p_bound_indexing_array = p_bound_indexing_array - reversed_p_bound_indexing_array
        
        p_bound_indexing_array   =pack(p_bound_indexing_array   ,p_bound_array)
        p_internal_indexing_array=pack(p_internal_indexing_array,p_internal_array)
        b_npoin = size(p_bound_indexing_array)
        i_npoin = size(p_internal_indexing_array)
        
        if (b_npoin+i_npoin .ne. npoin)then
            print*, 'ERROR ERROR'
            print*, 'sum of internal and external points is different to the total count of points.'
            stop
        endif
        
    end subroutine identify_boundary_points
    
    !!!!!!!!!!!!!!!!!! boundary flagging & normal vector routines
    
    subroutine calculate_normal_vectors
        use centroid_data
        implicit none
        integer(KIND=INT32) :: f, bf, fp_start, p1, p2, p3, p, bp, pf_start, pf_end, number_of_points, i, j
        real(KIND=REAL32)   :: v1(3), v2(3), sumV, magnitude, angle
        
        
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
            
            f_normal_vectors(bf,:)  = f_normal_vectors(bf,:) / sqrt(magnitude)
            
            
            ! we're taking the scalar product of the normal vector and the vector of the adjacent cell to the face to guarantee the normal faces outward.
            v1(:) = (c_centroid(f_c_obj_relation_array(f_c_index_array(f)),:) - f_centroid(f,:))
            v2(:) = f_normal_vectors(bf,:)
            !angle = (v1(1)*f_normal_vectors(bf,1) + v1(2)*f_normal_vectors(bf,2) + v1(3)*f_normal_vectors(bf,3))/ sqrt(v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3))
            angle = alignment(v1,v2)
            
            if (angle .lt. 0) f_normal_vectors(bf,:) = -f_normal_vectors(bf,:)
            
        enddo
        
        ! because I've already split my feature edges in my preprocessor I can then get point normal vectors by just averaging connected face normals
        
        allocate(p_normal_vectors(b_npoin,3))
        p_normal_vectors = 0.
        
        do bp=1,b_npoin
        
            p = p_bound_indexing_array(bp)
            
            pf_start = (1 + p_f_index_array(p-1))
            pf_end   = p_f_index_array(p)
            number_of_points = 1 + pf_end - pf_start
            
            do i=1,3
                
                do f=pf_start,pf_end
                    if (f_bound_array(p_f_obj_relation_array(f))) then
                        p_normal_vectors(bp,i) = p_normal_vectors(bp,i) + f_normal_vectors(reversed_f_bound_indexing_array(p_f_obj_relation_array(f)),i)
                    endif
                enddo
                
            enddo
            
            ! this is pretty elegant
            ! I'm dickriding myself pretty hard over this
            
        enddo
        
        do p=1,b_npoin
            p_normal_vectors(p,:) = p_normal_vectors(p,:) / sqrt(p_normal_vectors(p,1)**2 + p_normal_vectors(p,2)**2 + p_normal_vectors(p,3)**2)
        enddo
        
        
    end subroutine calculate_normal_vectors
    
    subroutine boundary_angle_feature_flagging
        implicit none
        integer(KIND=INT32) :: i, j, be, e, f1, f2
        real(KIND=REAL32)   :: nv1(3), nv2(3), angle
        
        allocate(feature_edges(b_nedge))
        feature_edges = .false.
        
        do be=1,b_nedge
            e = e_bound_indexing_array(be)
            
            j = e_f_index_array(e-1)
            do 
                j=j+1
                if (f_bound_array(e_f_obj_relation_array(j))) exit
            enddo
            f1 = e_f_obj_relation_array(j)
            
            do 
                j=j+1
                if (f_bound_array(e_f_obj_relation_array(j))) exit
            enddo
            f2 = e_f_obj_relation_array(j)
			
            if (j.gt.e_f_index_array(e))then
                print*, ' non-feature boundary edge has more than 2 faces '
                stop
            endif
            
            f1 = reversed_f_bound_indexing_array(f1)
            f2 = reversed_f_bound_indexing_array(f2)
            
            nv1(:)=f_normal_vectors(f1,:)
            nv2(:)=f_normal_vectors(f2,:)
            
            angle = acosd(alignment(nv1, nv2))
            
            if (angle.gt.15.) then
                call split_feature_edges(e,be, f1, f2)
            endif
            
        enddo
        
    end subroutine boundary_angle_feature_flagging
    
    subroutine split_feature_edges(e,be,f1,f2)
        implicit none
        integer(KIND=INT32) :: i, j, be, e, p1, p2, f1, f2
        
        
        feature_edges(be) = .true.
        
        e = e_bound_indexing_array(be)
            
        p1 = e_p_obj_relation_array(e_p_index_array(e)  )
        p2 = e_p_obj_relation_array(e_p_index_array(e)-1)
        
        
        
    end subroutine split_feature_edges
    
    subroutine dummy_feature_flagging
        implicit none
        
        
        
        allocate(feature_edges(b_nedge),feature_points(b_npoin))
        feature_edges  = .false.
        feature_points = .false.
        
        
        
    end subroutine dummy_feature_flagging
    
    ! end contains    
    
end module boundary_routine_module
