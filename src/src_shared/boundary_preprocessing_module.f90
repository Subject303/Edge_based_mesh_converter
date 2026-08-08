


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
        real(KIND=REAL64)   :: v1(3), v2(3), sumV, magnitude, angle
        
        
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
        integer(KIND=INT32) :: i, j, be, e, f1, f2, bp1, bp2
        real(KIND=REAL64)   :: nv1(3), nv2(3), angle
        
        allocate(feature_edges(b_nedge), feature_points(b_npoin))
        feature_edges = .false.
        feature_points = .false.
        
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
            
            angle = abs(acosd(alignment(nv1, nv2)))
            
            if (angle.gt.splitting_angle) then
            
                feature_edges(be) = .true.
                
                bp1 = reversed_p_bound_indexing_array(e_p_obj_relation_array(e_p_index_array(e)  ))
                bp2 = reversed_p_bound_indexing_array(e_p_obj_relation_array(e_p_index_array(e)-1))
                
                feature_points(bp1) = .true.
                feature_points(bp2) = .true.
                
                !call split_feature_edges(e,be, f1, f2)
            endif
            
        enddo
        
    end subroutine boundary_angle_feature_flagging
    
    
    subroutine boundary_region_face_flagging
        implicit none
        integer(KIND=INT32) :: i, be, bf, current_flag, ff_index, fe_index, ff_index_old, fe_index_old, f, e, fe_stt, fe_end, ef_stt, ef_end, fe, ef
        integer(KIND=INT32), allocatable :: flagged_faces(:), flagged_edges(:)
        logical :: all_edges_feature
        
        
        
        
        
        allocate(f_boundary_flags(b_nface), e_boundary_flags(b_nedge))
        
        f_boundary_flags = 999
        e_boundary_flags = 999
        
        ! so we have a list of flags, 999 for non feature edges, 1 for feature edges
        do be=1,b_nedge
            if (feature_edges(be)) e_boundary_flags(be) = 1000
        enddo
        
        bf=1
        current_flag = -1
        
        allocate(flagged_faces(b_nface), flagged_edges(b_nedge))
        
        flagged_faces = -999
        flagged_edges = -999
        
        ff_index = 1
        fe_index = 0
        
        ff_index_old = 1
        fe_index_old = 0
        
        flagged_faces(ff_index) = bf
        f_boundary_flags(bf) = current_flag
        
        
        do 
            fe_index_old = fe_index + 1
            
            !if (ff_index_old .ge. 11367) 
!             print*, 'FLAG 1', current_flag, fe_index_old, fe_index, ff_index_old, ff_index
            
            ! loop over all flagged faces we havn't already looped over
            do i=ff_index_old,ff_index
                bf = flagged_faces(i)
                
                f = f_bound_indexing_array(bf)
                
                fe_stt = f_e_index_array(f-1)+1
                fe_end = f_e_index_array(f)
                
                all_edges_feature = .true.
                
                ! loop over all flagged edges connected to flagged faces
                do fe=fe_stt,fe_end
                    be = reversed_e_bound_indexing_array(f_e_obj_relation_array(fe))
                    
                    ! if the flagged edge is not a feature and is not already flagged then flag it
                    !if ((e_boundary_flags(be).ne.1000).and.(e_boundary_flags(be).ne.current_flag)) then
                    
                    ! if the flagged edge is not already flagged
                    if (e_boundary_flags(be).eq.999) then
                        ! if any edges are non feature we know we can keep going
                        all_edges_feature = .false.
                        
                        fe_index = fe_index + 1
                        
                        flagged_edges(fe_index) = be
                        e_boundary_flags(be) = current_flag
                        
                        ! fe_index is the number of flagged edges
                    endif
                enddo
            enddo
            
            !if (ff_index_old .ge. 11367) then
!                 print*, ''
!                 print*, ''
!                 print*, flagged_faces
!                 print*, f_boundary_flags
!                 print*, ''
!                 print*, flagged_edges
!                 print*, e_boundary_flags
!                 print*, ''
!                 print*, ''
            !endif
            
            ! do somthing with exit conditions here
            ! all_edges_feature .true. stuff
            
            ff_index_old = ff_index + 1
            
            ! this is the condition to create a new boundary region
            if (all_edges_feature) then
                
                
                ! first check that we're not just already done.
                if (ff_index .eq. b_nface)then
                    exit
                endif
                
                ! otherwise we need to start a new flag and a new unseeded face
                current_flag = current_flag - 1
                
                do i=1,b_nface
                    if (f_boundary_flags(i).eq.999) then
                        bf = i
                        
                        ff_index = ff_index + 1
                        
                        flagged_faces(ff_index) = bf
                        f_boundary_flags(bf) = current_flag
                        
                        exit
                    endif
                enddo
                
!                 print*, 'NEW FLAG ', current_flag, fe_index_old, fe_index, ff_index_old, ff_index
                
                ! skip the next loop because there shouldnt be any connected faces
                cycle
                
            endif
            
            
            !if (ff_index_old .ge. 11367) 
!             print*, 'FLAG 2', current_flag, fe_index_old, fe_index, ff_index_old, ff_index
            
            ! now we update the flagged faces to restart the loop 
            
            ! so again loop over all flagged edges we havn't already looped over
            do i=fe_index_old,fe_index
                be = flagged_edges(i)
                
                e = e_bound_indexing_array(be)
                
                ef_stt = e_f_index_array(e-1)+1
                ef_end = e_f_index_array(e)
                
                ! then loop over connected faces
                do ef=ef_stt,ef_end
                    f = e_f_obj_relation_array(ef)
                    
                    ! ignore internal faces
                    if (f_internal_array(f)) cycle
                    
                    bf = reversed_f_bound_indexing_array(f)
                    
                    ! if the face is not already flagged
                    if (f_boundary_flags(bf).eq.999) then
                        ff_index = ff_index + 1
                        
                        ! and add all the boundary faces
                        flagged_faces(ff_index) = bf
                        f_boundary_flags(bf) = current_flag
                    endif
                    
                enddo
                
            enddo
            
            
        enddo
        
        do be=1,b_nedge
            
            if (e_boundary_flags(be) .eq. 1000) cycle
            
            current_flag = -999
            
            e = e_bound_indexing_array(be)
            
            ef_stt = e_f_index_array(e-1)+1
            ef_end = e_f_index_array(e)
            
            do ef=ef_stt,ef_end
                f = e_f_obj_relation_array(ef)
                
                ! ignore internal faces
                if (f_internal_array(f)) cycle
                
                bf = reversed_f_bound_indexing_array(f)
                
                if (current_flag .ne. -999) then
                    ! if the face is not already flagged
                    if (f_boundary_flags(bf).ne.current_flag) then
                        
                        if (f_boundary_flags(bf).gt.current_flag) then
                            do i=1, b_nface
                                if (f_boundary_flags(i) .eq. current_flag) f_boundary_flags(i) = f_boundary_flags(bf)
                            enddo
                            
                            do i=1, b_nedge
                                if (e_boundary_flags(i) .eq. current_flag) e_boundary_flags(i) = f_boundary_flags(bf)
                            enddo
                        elseif (f_boundary_flags(bf).lt.current_flag) then
                            do i=1, b_nface
                                if (f_boundary_flags(i) .eq. f_boundary_flags(bf)) f_boundary_flags(i) = current_flag
                            enddo
                            
                            do i=1, b_nedge
                                if (e_boundary_flags(i) .eq. f_boundary_flags(bf)) e_boundary_flags(i) = current_flag
                            enddo
                        endif
                            
                            
                        exit
                    endif
                endif
                
                current_flag = f_boundary_flags(bf)
                
            enddo
        enddo
        
        
    end subroutine boundary_region_face_flagging
     
    subroutine split_corners
        implicit none
        integer(KIND=INT32) :: i, j, f, p, be, bf, bp, pf, pf_stt, pf_end, bp1, bp2, e, new_b_npoin, num_feature_points, num_feat_projections, new_bp, flag
        integer(KIND=INT32),allocatable :: number_of_projections(:), temp_p_bound_indexing_array(:)
        logical, allocatable :: temp_feature_points(:)
        
        ! ok so
        ! realistically here I want to reconstruct my boundry points into to lists, non feature points and feature points
        ! at the minute I have bp indexes and p indexes
        ! this means I can scan connected things to p from bp, and I can have multiple bp occupy the same p index
        ! so ideally I can just add however many projections a p has as extra bp 
        ! this will break the reverse list so a condition of the reversed list is non feature points
        ! as long as this is the last step in processing that should be fine
        !
!         reversed_p_bound_indexing_array
!         p_bound_indexing_array
!         p_internal_indexing_array
!         b_npoin
!         i_npoin
!         p_boundary_flags
!         feature_points
        
        ! first step i think is to count the number of boundaries each point will have
        ! we can do this by adding the number of feature edges connected to a point, all non feature points have a count of 1
        
        allocate(number_of_projections(b_npoin))
        number_of_projections = 0
        
        ! we make the number of projections negative to play with the sorting algo later dw
        do be=1,b_nedge
            
            if (feature_edges(be))then
                
                e = e_bound_indexing_array(be)
            
                bp1 = reversed_p_bound_indexing_array(e_p_obj_relation_array(e_p_index_array(e)  ))
                bp2 = reversed_p_bound_indexing_array(e_p_obj_relation_array(e_p_index_array(e)-1))
                
                number_of_projections(bp1) = number_of_projections(bp1) - 1
                number_of_projections(bp2) = number_of_projections(bp2) - 1
            endif
            
        enddo
        
        do bp=1,b_npoin
            if (number_of_projections(bp).eq.0) number_of_projections(bp) = 1
        enddo
        
        ! this means lines of feature edges will split into two projections, the meeting point of 3 edges becomes 3 bounds ect
        ! feature edges that lead into a smooth surface will only have 1 projection on the last point and 2 on the next
        
        ! fix the positivity of the counts
        number_of_projections(:) = abs(number_of_projections(:))
        
        new_b_npoin = sum(number_of_projections)
        allocate(temp_p_bound_indexing_array,source=p_bound_indexing_array)
        
        ! this will sort the p indexing array by number of bounds
        ! largest number of bounds first
        call quicksort(number_of_projections,1,b_npoin,temp_p_bound_indexing_array)
        
        feature_points = .false.
        do bp=1,b_npoin
            if (number_of_projections(bp) .eq. 1) then
                feature_points(bp) = .false.
            else
                feature_points(bp) = .true.
            endif
        enddo
        
        ! this finds the last feature point in the index so we can seperate these out
        num_feature_points = 0
        do bp=1,b_npoin-1
            if (number_of_projections(bp+1).eq.1) num_feature_points = bp
        enddo
        
        ! so now we can loop from 1:num_feature_points or num_feature_points+1:b_npoin
        ! if we want feat or non feat
        
        num_feat_projections = sum(number_of_projections(1:num_feature_points))
        
        deallocate(p_bound_indexing_array)
        allocate(p_bound_indexing_array(new_b_npoin), p_boundary_flags(new_b_npoin))
        
        allocate(temp_feature_points,source=feature_points)
        deallocate(feature_points, p_normal_vectors)
        allocate(feature_points(new_b_npoin), p_normal_vectors(new_b_npoin,3))
        feature_points = .false.
        p_normal_vectors = 0.0
        
        new_bp = 0
        do bp=1,b_npoin
            
            flag = 999
            
            p = temp_p_bound_indexing_array(bp)
                
            pf_stt = p_f_index_array(p-1)+1
            pf_end = p_f_index_array(p)
            
            pf = pf_stt
            
            do i=1,number_of_projections(bp)
                new_bp = new_bp + 1
                
                feature_points(new_bp) = temp_feature_points(bp)
                
                p_bound_indexing_array(new_bp) = p
                
                ! this is so we only look forward into the face array
                do pf=pf,pf_end
                    
                    f = p_f_obj_relation_array(pf)
                    
                    if (f_internal_array(f)) cycle
                    
                    bf = reversed_f_bound_indexing_array(f)
                    
                    if (f_boundary_flags(bf) .ne. flag)then
                        flag = f_boundary_flags(bf)
                        exit
                    endif
                    
                    
                enddo
                
                p_boundary_flags(new_bp) = flag
                
                print*,p_boundary_flags(new_bp)
                
            enddo
            
        enddo
        
        b_npoin = new_b_npoin
        
!         reversed_p_bound_indexing_array
!         p_bound_indexing_array          ! done
!         p_internal_indexing_array
!         b_npoin                         ! done
!         i_npoin
!         p_boundary_flags                ! done
!         feature_points                  ! done
        
        ! I dont think any of the other arrays are necessarily in use post here anyway
        ! so I'm going to be lazy
        ! wish me luck
        
        ! forgot normals :(
        

        do bp=1,b_npoin
            
            flag = p_boundary_flags(bp)
            
            p = p_bound_indexing_array(bp)
            
            pf_stt = p_f_index_array(p-1)+1
            pf_end = p_f_index_array(p)
            do pf=pf_stt,pf_end
            
                f = p_f_obj_relation_array(pf)
                
                if (f_internal_array(f)) cycle
                
                bf = reversed_f_bound_indexing_array(f)
                
                if (flag .eq. f_boundary_flags(bf)) then
                    p_normal_vectors(bp,:) = p_normal_vectors(bp,:) + f_normal_vectors(bf,:)
                endif
            enddo
        enddo
        
        do bp=1,b_npoin
            p_normal_vectors(bp,:) = p_normal_vectors(bp,:) / sqrt(p_normal_vectors(bp,1)**2 + p_normal_vectors(bp,2)**2 + p_normal_vectors(bp,3)**2)
        enddo
        
        
        
        
    end subroutine split_corners

    subroutine dummy_feature_flagging
        implicit none
        
        
        allocate(feature_edges(b_nedge),feature_points(b_npoin))
        feature_edges  = .false.
        feature_points = .false.
        
        allocate(p_boundary_flags(b_npoin))
        
        
        
        
        p_boundary_flags = -999
        
    end subroutine dummy_feature_flagging
    
    ! end contains    
    
end module boundary_routine_module
