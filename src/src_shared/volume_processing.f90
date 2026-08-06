

module volume_processing
    use iso_fortran_env
    use object_counts
    use raw_data
    use object_relation_data
    use centroid_data
    use boundary_data
    use quicksort_module
    use projection_data
    use utils, ONLY : alignment, planar_alignment
    ! need everything p much here.
    implicit none
        
    contains
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine volume_alloc
        implicit none
        allocate(sn(i_nedge,3),sb(b_nedge,3),sbb(b_npoin,3))
        
        sn  = 0.0
        sb  = 0.0
        sbb = 0.0
        
        
    end subroutine volume_alloc
        
    subroutine internal_edge_volume_processing
        implicit none
        integer(KIND=INT32) :: ie, e, i, j, i1, i2, centroid_array_count, centroid_array_count_old, ef, ec, ef_start, ec_start, ef_end, ec_end
        integer(KIND=INT32) :: cell_count, face_count, current_face, cell_1, cell_2, prev_cell
        integer(KIND=INT32),allocatable :: centroid_index_array(:)
        real(KIND=REAL64),allocatable   :: centroid_array(:,:)
        real(KIND=REAL64)               :: in_progress_projection(3), in_progress_centroid(3), direction_array(3), c1(3), c2(3), angle, vol1, vol2
        logical, allocatable :: non_viable_faces(:)
        
        ! it's upsetting this is the easiest of the three jobs I gotta do
        
        ! so, step wise here's the process per edge
        ! 1:
        !   find all connected cells and faces
        ! 2:
        !   order an array of cells and faces
        !   hopefully this is not too dificult
        ! 3:
        !   call the processer subroutine and sum the projections and volumes
        
        !centroid_array_count_old = -1
        allocate(centroid_index_array(0),non_viable_faces(0),centroid_array(0,0))
        
        do ie = 1, i_nedge
            e = e_internal_indexing_array(ie)
            ! this is a loop of all internal edges.
            
!             i1 = e_p_obj_relation_array(e_p_index_array(e)-1)
!             i2 = e_p_obj_relation_array(e_p_index_array(e))
            ! i1 and i2 are the constituent points of edge e
            
            call centroid_assembler(e, 1, centroid_index_array)
            
            
            
!             ef_start = e_f_index_array(e-1)
!             ec_start = e_c_index_array(e-1)
!             ef_end   = e_f_index_array(e)
!             ec_end   = e_c_index_array(e)
!             
!             cell_count = ec_end-ec_start
!             face_count = ef_end-ef_start
!             
!             centroid_array_count = 1 + cell_count + face_count
!             
!             ! sum of number of faces, number of edges plus 1 for the duplicate starting edge
!             
!             if (centroid_array_count_old .ne. centroid_array_count) then
!                 deallocate(centroid_index_array,non_viable_faces,centroid_array)
!                 allocate(centroid_index_array(centroid_array_count), non_viable_faces(face_count+1), centroid_array(centroid_array_count,3))
!             endif
!             
!             non_viable_faces = .false.
!             
!             current_face = e_f_obj_relation_array(ef_start+1)
!             cell_1 = f_c_obj_relation_array(f_c_index_array(current_face)    )
!             cell_2 = f_c_obj_relation_array(f_c_index_array(current_face-1)+1)
!             
!             centroid_index_array(1) = cell_1
!             centroid_index_array(2) = current_face
!             centroid_index_array(3) = cell_2
!             
!             centroid_array(1,:) = c_centroid(cell_1,:)
!             centroid_array(2,:) = f_centroid(current_face,:)
!             centroid_array(3,:) = c_centroid(cell_2,:)
!                     
!             prev_cell = cell_2
!             
!             !print*, 'number of faces ', face_count, ' number of cells ', cell_count, ' count ', centroid_array_count
!             
!             i=4
!             
!             ef = 1
!             
!             !print*, ef, current_face, cell_1, cell_2 ,prev_cell, 'prev_cell'
!             
!             do
!                 
!                 if (ef.eq.face_count) then
!                     ef = 1
!                 else
!                     ef=ef+1
!                 endif
! 
!                 
!                 if (non_viable_faces(ef)) cycle
!                 
!                 current_face = e_f_obj_relation_array(ef_start + ef)
!                 cell_1 = f_c_obj_relation_array(f_c_index_array(current_face)    )
!                 cell_2 = f_c_obj_relation_array(f_c_index_array(current_face-1)+1)
!                 
!                 if (prev_cell .eq. cell_1) then
!             
!                     centroid_index_array(i)   = current_face
!                     centroid_index_array(i+1) = cell_2
!                     
!                     centroid_array(i,:)   = f_centroid(current_face,:)
!                     centroid_array(i+1,:) = c_centroid(cell_2,:)
!                     
!                     non_viable_faces(ef) = .true.
!                     
!                     prev_cell = centroid_index_array(i+1)
!                     i=i+2 
!                     
!                     !print*, 'i-1 = cell_1, i+1 = cell_2'
!                     
!                     if (i.eq.centroid_array_count+1) exit
!                     
!                 elseif (prev_cell .eq. cell_2) then
!                     
!                     centroid_index_array(i)   = current_face
!                     centroid_index_array(i+1) = cell_1
!                     
!                     centroid_array(i,:)   = f_centroid(current_face,:)
!                     centroid_array(i+1,:) = c_centroid(cell_1,:)
!                     
!                     non_viable_faces(ef) = .true.
!                     
!                     prev_cell = centroid_index_array(i+1)
!                     i=i+2 
!                     
!                     !print*, 'swapped'
!                     
!                     if (i.eq.centroid_array_count+1) exit
!                     
!                 else
!                     !print*, 'none, looping'
!                 endif
!                 
!                 !print*, ef, current_face, cell_1, cell_2 , prev_cell
!                 
!             enddo
!             
!             ! this might be exactly negative, alternatives are 1, and cell_1
!             centroid_index_array(centroid_array_count) = centroid_index_array(1)
!             centroid_array(centroid_array_count,:) = c_centroid(centroid_index_array(1),:)
! 			
!             centroid_array_count_old = centroid_array_count
!             
!             
!             !if (i.ne.centroid_array_count+1) print*, 'centroid array count broken'
!             !if (centroid_index_array(3).ne.centroid_index_array(centroid_array_count)) print*, 'internal centroid array start and end wrong'
!             
!             !print*, centroid_index_array
!             
!             in_progress_projection(:) = 0.0
!             in_progress_centroid = e_centroid(e,:)
!             
!             c1(:) = coords(i1,:)
!             c2(:) = coords(i2,:)
!             
!             call centroid_array_routine(in_progress_projection,  in_progress_centroid, centroid_array_count, centroid_array)
!             
!             direction_array(:) = c1(:) - c2(:)
!             
!             angle = alignment(in_progress_projection, direction_array)
!             
!             if (angle .lt. 0.) then
!                 sn(ie,:) = in_progress_projection(:)
! 			else
!                 sn(ie,:) = -in_progress_projection
! 			endif
!             
! !             print*,'aa'
! !             !print*, centroid_array_count, SIZE(centroid_array,1)
! !             print*, ie, ' :: ', i1-1, i2-1
! !             print*, in_progress_projection
! !             print*, sn(ie,:)
! !             print*, angle
! !             print*,'aa'
! !             print*,in_progress_centroid
! !             print*,'aa'
! !             do i=1, centroid_array_count
! !                 print*,i,' :: ', centroid_array(i,:)
! !             enddo
! !             print*,'aa'
!             
        enddo
        
        !print*,'aaaa'
        !do i=1,i_nedge
        !    print*, i, sqrt(sn(i,1)*sn(i,1) + sn(i,2)*sn(i,2) + sn(i,3)*sn(i,3)), sn(i,:)
        !enddo
        !print*,'aaaa'
        
        print*, 'backstop i_edge volume procs'
        stop
        
    end subroutine internal_edge_volume_processing
    
        
    subroutine boundary_edge_volume_processing
        implicit none
        integer(KIND=INT32) :: be, e, i, j, i1, i2, centroid_array_count, centroid_array_count_old, ef, ec, ef_start, ec_start, ef_end, ec_end
        integer(KIND=INT32) :: cell_count, face_count, current_face, cell_1, cell_2, prev_cell
        integer(KIND=INT32),allocatable :: centroid_index_array(:)
        real(KIND=REAL64),allocatable   :: centroid_array(:,:)
        real(KIND=REAL64)               :: in_progress_projection(3), in_progress_centroid(3), direction_array(3), c1(3), c2(3), angle, vol1, vol2
        logical, allocatable :: non_viable_faces(:)
        
!         do e=1,nedge
!             i1       = e_p_index_array(e-1)	+ 1
! 			ef_start = e_f_index_array(e-1)	+ 1
!             ec_start = e_c_index_array(e-1)	+ 1
!             i2       = e_p_index_array(e) 
!             ef_end   = e_f_index_array(e)
!             ec_end   = e_c_index_array(e)
!             print*, 'bound ', e_bound_array(e), e_centroid(e,:)
!             print*, 'point'
!             print*, e_p_obj_relation_array(i1:i2)
!             print*, p_bound_array(e_p_obj_relation_array(i1:i2))
!             print*, 'edge'
!             print*, e_f_obj_relation_array(ef_start:ef_end)
!             print*, f_bound_array(e_f_obj_relation_array(ef_start:ef_end))
!             print*, 'face'
!             print*, e_c_obj_relation_array(ec_start:ec_end)
!             
! 	    enddo
        
        
        ! it's upsetting this is the easiest of the three jobs I gotta do
        
        ! so, step wise here's the process per edge
        ! 1:
        !   find all connected cells and faces
        ! 2:
        !   order an array of cells and faces
        !   hopefully this is not too dificult
        ! 3:
        !   call the processer subroutine and sum the projections and volumes
        
        centroid_array_count_old = -1
        allocate(centroid_index_array(0),non_viable_faces(0),centroid_array(0,0))
        
        do be = 1, b_nedge
            e = e_bound_indexing_array(be)
            ! this is a loop of all boundary edges.
            
            i1 = e_p_obj_relation_array(e_p_index_array(e)-1)
            i2 = e_p_obj_relation_array(e_p_index_array(e))
            ! i1 and i2 are the constituent points of edge e
            
            ef_start = e_f_index_array(e-1)
            ec_start = e_c_index_array(e-1)
            ef_end   = e_f_index_array(e)
            ec_end   = e_c_index_array(e)
            
            cell_count = ec_end-ec_start
            face_count = ef_end-ef_start
            
            centroid_array_count = cell_count + face_count
            
            ! sum of number of faces, number of edges plus 1 for the duplicate starting edge
            
            if (centroid_array_count_old .ne. centroid_array_count) then
                deallocate(centroid_index_array,non_viable_faces,centroid_array)
                allocate(centroid_index_array(centroid_array_count), non_viable_faces(face_count), centroid_array(centroid_array_count,3))
            endif
            
            ef = 1
            
            non_viable_faces = .false.
            
            do ! we must start on a boundary face
                current_face = e_f_obj_relation_array(ef_start+ef)
!                 print*, current_face, f_c_obj_relation_array(f_c_index_array(current_face)    ), f_c_obj_relation_array(f_c_index_array(current_face-1)+1), f_bound_array(current_face), ef, ef_start+ef, ef_end
                if (f_bound_array(current_face)) exit
                ef = ef + 1
            enddo

            non_viable_faces(ef) = .true.
            
            cell_1 = f_c_obj_relation_array(f_c_index_array(current_face)    )
            
            centroid_index_array = -1
            centroid_index_array(1) = current_face
            centroid_index_array(2) = cell_1
            
            centroid_array(1,:) = f_centroid(current_face,:)
            centroid_array(2,:) = c_centroid(cell_1,:)
                    
            prev_cell = cell_1
            cell_2 = -1
            
            !print*, 'number of faces ', face_count, ' number of cells ', cell_count, ' count ', centroid_array_count
            
            i=3
            
            !print*, ef, current_face, cell_1, cell_2 ,prev_cell, 'prev_cell'
!             print*, 'aa ', centroid_index_array
            do
                
                if (ef.eq.face_count) then
                    ef = 1
                else
                    ef=ef+1
                endif

                
                if (non_viable_faces(ef)) cycle
                
                current_face = e_f_obj_relation_array(ef_start + ef)
                cell_1 = f_c_obj_relation_array(f_c_index_array(current_face)    )
                cell_2 = f_c_obj_relation_array(f_c_index_array(current_face-1)+1)
                
                if (prev_cell .eq. cell_1) then
            
                    centroid_index_array(i)   = current_face
                    centroid_array(i,:)   = f_centroid(current_face,:)
                    
!                     print*, current_face, cell_1, cell_2, f_bound_array(current_face), ef, ef_start+ef, ef_end
                    if (f_bound_array(current_face)) exit
                    
                    centroid_index_array(i+1) = cell_2
                    centroid_array(i+1,:) = c_centroid(cell_2,:)
                    
                    non_viable_faces(ef) = .true.
                    
                    prev_cell = centroid_index_array(i+1)
                    i=i+2 
                    
                    !print*, 'i-1 = cell_1, i+1 = cell_2'
                    
                elseif (prev_cell .eq. cell_2) then
                    
                    centroid_index_array(i)   = current_face
                    centroid_array(i,:)   = f_centroid(current_face,:)
                    
!                     print*, current_face, cell_1, cell_2, f_bound_array(current_face), ef, ef_start+ef, ef_end
                    if (f_bound_array(current_face)) exit
                    
                    centroid_index_array(i+1) = cell_1
                    centroid_array(i+1,:) = c_centroid(cell_1,:)
                    
                    non_viable_faces(ef) = .true.
                    
                    prev_cell = centroid_index_array(i+1)
                    i=i+2 
                    
                    !print*, 'swapped'
                    
                else
                    !print*, 'none, looping'
                endif
                
                !print*, ef, current_face, cell_1, cell_2 , prev_cell
!                 print*, 'aa ', centroid_index_array
                
            enddo
			
!             print*, 'aa ', centroid_index_array
            
            
            centroid_array_count_old = centroid_array_count
            
            
            !if (i.ne.centroid_array_count+1) print*, 'centroid array count broken'
            !if (centroid_index_array(1).ne.centroid_index_array(centroid_array_count)) print*, 'internal centroid array start and end wrong'
            
            
            in_progress_projection(:) = 0.0
            in_progress_centroid = e_centroid(e,:)
            
            c1(:) = coords(i1,:)
            c2(:) = coords(i2,:)
            
!             print*, be, (c1(1)+c2(1))/2, (c1(2)+c2(2))/2, (c1(3)+c2(3))/2, in_progress_centroid
            
            call centroid_array_routine(in_progress_projection,  in_progress_centroid, centroid_array_count, centroid_array)
            
            
            direction_array(:) = c1(:) - c2(:)
            
            angle = alignment(in_progress_projection, direction_array)
            
            if (angle .lt. 0.) then
                sb(be,:) = in_progress_projection(:)
			else
                sb(be,:) = -in_progress_projection
			endif
            
!             print*,'aa'
!             !print*, centroid_array_count, SIZE(centroid_array,1)
!             print*, be, ' :: ', i1-1, i2-1, ' :: ', centroid_array_count
!             print*, in_progress_projection
!             print*, sb(be,:)
!             print*, angle
!             print*,'aa'
!             print*,in_progress_centroid
!             print*,'aa'
!             do i=1, centroid_array_count
!                 print*,i,' :: ', centroid_array(i,:)
!             enddo
!             print*,'aa'
            
        enddo
        
        !do i=1,b_nedge
        !    print*, i, sqrt(sb(i,1)*sb(i,1) + sb(i,2)*sb(i,2) + sb(i,3)*sb(i,3)), sb(i,:)
        !enddo
        !print*,'aaaa'
        
    end subroutine boundary_edge_volume_processing
    
    
    subroutine boundary_face_volume_processing
        implicit none
        integer(KIND=INT32) :: bp, p, i, j, i1, i2, centroid_array_count, centroid_array_count_real, centroid_array_count_old, pe, pe_start, pf_start, pe_end, pf_end
        integer(KIND=INT32) :: edge_count, face_count, current_edge, face_1, face_2, prev_face,be
        integer(KIND=INT32) :: k, pf, edge_count_real, face_count_real, flag, remainder_count, current_face
        integer(KIND=INT32),allocatable :: centroid_index_array(:)
        real(KIND=REAL64),allocatable   :: centroid_array(:,:)
        real(KIND=REAL64)               :: in_progress_projection(3), in_progress_centroid(3), direction_array(3), c1(3), c2(3), angle
        logical, allocatable :: non_viable_edges(:)
        
        ! it's upsetting this is the easiest of the three jobs I gotta do
        
        ! so, step wise here's the process per edge
        ! 1:
        !   find all connected cells and faces
        ! 2:
        !   order an array of cells and faces
        !   hopefully this is not too dificult
        ! 3:
        !   call the processer subroutine and sum the projections and volumes
        
        centroid_array_count_old = -1
        allocate(centroid_index_array(0),non_viable_edges(0),centroid_array(0,0))
        
        do bp = 1, b_npoin
            p = p_bound_indexing_array(bp)
            ! this is a loop of all boundary edges.
            
            
            i1 = p
            i2 = p
            ! i1 and i2 are the constituent points of edge e
            
            
            
            if (feature_points(bp)) then
                
                ! FEATURE POINTS
                ! IE CORNERS
                
                ! so remember in boundary preprocessing, I've flagged faces and points with boundary flags already
                ! this means we just need to produce the edge, face array comprised only of the correct boundary flag 
                ! AND remembering to start on the correct two feature edges
                
                flag = p_boundary_flags(bp)
                
                pe_start = p_e_index_array(p-1)
                pf_start = p_f_index_array(p-1)
                pe_end   = p_e_index_array(p)
                pf_end   = p_f_index_array(p)
                
                edge_count = pe_end-pe_start
                
                edge_count_real = 2 ! remember that feature edges wont be counted here
                do pe=(pe_start+1),pe_end
                    current_edge = p_e_obj_relation_array(pe)
                    if (e_internal_array(current_edge)) cycle
                    if (e_boundary_flags(reversed_e_bound_indexing_array(current_edge)) .eq. flag) edge_count_real = edge_count_real + 1
                enddo
                
                face_count_real = 0
                do pf=(pf_start+1),pf_end
                    current_face = p_f_obj_relation_array(pf)
                    if (f_internal_array(current_face)) cycle
                    if (f_boundary_flags(reversed_f_bound_indexing_array(current_face)) .eq. flag) face_count_real = face_count_real + 1
                enddo
                
                centroid_array_count = edge_count_real + face_count_real
                
                ! ok so here we're going to pick a random edge in the list that does have the correct flag, and search connections forward until we hit a feature edge
                ! then count how many short we are and shuffle the array forward and search backwards
                
                if (centroid_array_count_old .ne. centroid_array_count) then
                    deallocate(centroid_index_array,non_viable_edges,centroid_array)
                    allocate(centroid_index_array(centroid_array_count), non_viable_edges(edge_count), centroid_array(centroid_array_count,3))
                endif
                
                non_viable_edges = .false.
                
                ! flag all wrong flag edges
                do pe=1, edge_count
                    current_edge  = p_e_obj_relation_array(pe_start+pe)
                    if (e_internal_array(current_edge)) then
                        non_viable_edges(pe) = .true.
                        cycle
                    endif
                    
                    be = reversed_e_bound_indexing_array(current_edge)
                    if ((e_boundary_flags(be).ne.flag) .or. (.not.feature_edges(be)) ) non_viable_edges(pe) = .true.
                    
                enddo
                
                ! select a correctly flagged edge
                do pe=1,edge_count
                    if (non_viable_edges(pe)) cycle
                    
                    current_edge  = p_e_obj_relation_array(pe_start+pe)
                    be = reversed_e_bound_indexing_array(current_edge)
                    
                    if (e_boundary_flags(be).eq.flag) exit
                enddo
                
                ! we dont need to double back on edges here
                non_viable_edges(pe) = .true.
                
                ! find two correctly flagged faces
                j = e_f_index_array(current_edge-1)
                do 
                    j=j+1
                    if (f_bound_array(e_f_obj_relation_array(j))) then
                        if (f_boundary_flags(reversed_f_bound_indexing_array(e_f_obj_relation_array(j))) .eq. flag ) exit
                    endif
                enddo
                face_1 = e_f_obj_relation_array(j)
                do 
                    j=j+1
                    if (f_bound_array(e_f_obj_relation_array(j))) then
                        if (f_boundary_flags(reversed_f_bound_indexing_array(e_f_obj_relation_array(j))) .eq. flag ) exit
                    endif
                enddo
                face_2 = e_f_obj_relation_array(j)
                
                if (j.gt.e_f_index_array(current_edge))then
                    print*, ' non-feature boundary edge has more than 2 faces '
                    stop
                endif
                
                
                
                ! ok so we set positions 2,3,4 as this edge because we know 100% that we cant start on a flagged edge
                centroid_index_array(2) = face_1
                centroid_index_array(3) = current_edge
                centroid_index_array(4) = face_2
            
                centroid_array_count_real = 4
                
                centroid_array(2,:) = f_centroid(face_1,:)
                centroid_array(3,:) = e_centroid(current_edge,:)
                centroid_array(4,:) = f_centroid(face_2,:)
                
                !print*, bp , face_1, current_edge, face_2
                
                prev_face = face_2
                
                ! i is the forward index
                ! k is the backwards index
                i=5
                k=1
                
                do
                
                    if (pe.eq.edge_count) then
                        pe = 1
                    else
                        pe=pe+1
                    endif
                    

!                     print*, pe, non_viable_edges

                    ! skip anything inviable
                    if (non_viable_edges(pe)) cycle
                
                    current_edge = p_e_obj_relation_array(pe_start+pe)
                    be = reversed_e_bound_indexing_array(current_edge)
                    
                    if (.not.feature_edges(be)) then
                        ! find two faces and proceed as normal
                        j = e_f_index_array(current_edge-1)
                        do 
                            j=j+1
                            if (f_bound_array(e_f_obj_relation_array(j))) then
                                if (f_boundary_flags(reversed_f_bound_indexing_array(e_f_obj_relation_array(j))) .eq. flag ) exit
                            endif
                        enddo
                        face_1 = e_f_obj_relation_array(j)
                        do 
                            j=j+1
                            if (f_bound_array(e_f_obj_relation_array(j))) then
                                if (f_boundary_flags(reversed_f_bound_indexing_array(e_f_obj_relation_array(j))) .eq. flag ) exit
                            endif
                        enddo
                        face_2 = e_f_obj_relation_array(j)
                        if (j.gt.e_f_index_array(current_edge))then
                            print*, ' non-feature boundary edge has more than 2 faces '
                            stop
                        endif
                        
                        if (prev_face .eq. face_1) then
                            
                            if (i+2 .gt. centroid_array_count) then
                                
                                do j=k,i
                                    centroid_index_array(j-2) = centroid_index_array(j)
                                    centroid_array(j-2,:)     = centroid_array(j,:)
                                enddo
                                
                                i = i - 2
                                
                            endif
                            
                            centroid_index_array(i)   = current_edge
                            centroid_index_array(i+1) = face_2
                        
                            centroid_array(i,:)   = e_centroid(current_edge,:)
                            centroid_array(i+1,:) = f_centroid(face_2,:)
                        
                            centroid_array_count_real = centroid_array_count_real + 2
                            
                            non_viable_edges(pe) = .true.
                        
                            prev_face = centroid_index_array(i+1)
                            i=i+2 
                        
                            print*, 'i-1 = face_1, i+1 = face_2',flag,e_boundary_flags(reversed_e_bound_indexing_array(current_edge)),feature_edges(reversed_e_bound_indexing_array(current_edge))!, i, centroid_array_count+1
                        
                        elseif (prev_face .eq. face_2) then
                        
                            if (k-2 .lt. 1) then
                                
                                do j=k,i,-1
                                    centroid_index_array(j+2) = centroid_index_array(j)
                                    centroid_array(j+2,:)     = centroid_array(j,:)
                                enddo
                                
                                k = k + 2
                                
                            endif
                        
                            centroid_index_array(k)   = current_edge
                            centroid_index_array(k-1) = face_1
                        
                            centroid_array(k,:)   = e_centroid(current_edge,:)
                            centroid_array(k-1,:) = f_centroid(face_1,:)
                        
                            centroid_array_count_real = centroid_array_count_real + 2
                            
                            non_viable_edges(pe) = .true.
                        
                            prev_face = centroid_index_array(k-1)
                            k=k-2 
                        
                            print*, 'swapped',flag,e_boundary_flags(reversed_e_bound_indexing_array(current_edge)),feature_edges(reversed_e_bound_indexing_array(current_edge))!, i, centroid_array_count+1
                        
                        else
                            print*, 'none, looping',flag,e_boundary_flags(reversed_e_bound_indexing_array(current_edge)),feature_edges(reversed_e_bound_indexing_array(current_edge))!, i, centroid_array_count+1, edge_count, face_count
                            
                        endif
                        
                        
                    else
                        ! only find a single face
                        ! and make sure its flagged right
                        j = e_f_index_array(current_edge-1)
                        do 
                            j=j+1
                            if (f_bound_array(e_f_obj_relation_array(j))) then
                                if (f_boundary_flags(reversed_f_bound_indexing_array(e_f_obj_relation_array(j))) .eq. flag ) exit
                            endif
                        enddo
                        face_1 = e_f_obj_relation_array(j)
                        
                        
                        
                    endif
                    
                    
                
!                     print*, pe, current_edge, face_1, face_2 , prev_face
                    
                    
                    print*, ' '
                    print*, pe, current_edge, face_1, face_2 , prev_face
                    print*, ' '
                    print*, 'centroid array ', centroid_index_array
                    print*, ' '
                    print*, ' '
                    print*, ' '
                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                enddo
                
                print*, 'centroid array end', centroid_index_array
                
            else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! NON FEATURE POINTS
                ! IE NOT CORNERS
                
                pe_start = p_e_index_array(p-1)
                pf_start = p_f_index_array(p-1)
                pe_end   = p_e_index_array(p)
                pf_end   = p_f_index_array(p)
                
                edge_count = pe_end-pe_start
                face_count = pf_end-pf_start
                
                centroid_array_count = 1 + edge_count + face_count
            
                ! sum of number of faces, number of edges plus 1 for the duplicate starting edge
            
                if (centroid_array_count_old .ne. centroid_array_count) then
                    deallocate(centroid_index_array,non_viable_edges,centroid_array)
                    allocate(centroid_index_array(centroid_array_count), non_viable_edges(edge_count), centroid_array(centroid_array_count,3))
                endif
            
                pe = 1
                
                non_viable_edges = .false.
            
                do ! we must only use boundary edges
                    current_edge = p_e_obj_relation_array(pe_start+pe)
                    if (e_bound_array(current_edge)) exit
                    non_viable_edges(pe) = .true.
                    
                    if (pe.eq.edge_count) then
                        pe = 1
                    else
                        pe=pe+1
                    endif
                enddo
                
                j = e_f_index_array(current_edge-1)
                do 
                    j=j+1
                    if (f_bound_array(e_f_obj_relation_array(j))) exit
                enddo
                face_1 = e_f_obj_relation_array(j)
                do 
                    j=j+1
                    if (f_bound_array(e_f_obj_relation_array(j))) exit
                enddo
                face_2 = e_f_obj_relation_array(j)
                
                if (j.gt.e_f_index_array(current_edge))then
                    print*, ' non-feature boundary edge has more than 2 faces '
                    stop
                endif
                
                centroid_index_array(1) = face_1
                centroid_index_array(2) = current_edge
                centroid_index_array(3) = face_2
            
                centroid_array_count_real = 3
                
                centroid_array(1,:) = f_centroid(face_1,:)
                centroid_array(2,:) = e_centroid(current_edge,:)
                centroid_array(3,:) = f_centroid(face_2,:)
                
                prev_face = face_2
            
                !print*, 'number of faces ', face_count, ' number of cells ', cell_count, ' count ', centroid_array_count
            
                i=4
            
                !print*, pe, current_edge, face_1, face_2 , prev_face
                
                do
                
                    if (pe.eq.edge_count) then
                        pe = 1
                    else
                        pe=pe+1
                    endif

                    !print*, pe, non_viable_edges
                    if (non_viable_edges(pe)) cycle
                
                    do ! we must only use boundary edges
                        current_edge = p_e_obj_relation_array(pe_start+pe)
                        if (e_bound_array(current_edge)) exit
                        non_viable_edges(pe) = .true.
                        
                        if (pe.eq.edge_count) then
                            pe = 1
                        else
                            pe=pe+1
                        endif
                        
                    enddo
                    
                    j = e_f_index_array(current_edge-1)
                    do 
                        j=j+1
                        if (f_bound_array(e_f_obj_relation_array(j))) exit
                    enddo
                    face_1 = e_f_obj_relation_array(j)
                    do 
                        j=j+1
                        if (f_bound_array(e_f_obj_relation_array(j))) exit
                    enddo
                    face_2 = e_f_obj_relation_array(j)
                    
                    if (j.gt.e_f_index_array(current_edge))then
                        print*, ' non-feature boundary edge has more than 2 faces '
                        stop
                    endif
                        
                    if (prev_face .eq. face_1) then
            
                        centroid_index_array(i)   = current_edge
                        centroid_index_array(i+1) = face_2
                    
                        centroid_array(i,:)   = e_centroid(current_edge,:)
                        centroid_array(i+1,:) = f_centroid(face_2,:)
                    
                        centroid_array_count_real = centroid_array_count_real + 2
                        
                        non_viable_edges(pe) = .true.
                    
                        prev_face = centroid_index_array(i+1)
                        i=i+2 
                    
                        !print*, 'i-1 = face_1, i+1 = face_2', i, centroid_array_count+1
                    
                        if (face_2 .eq. centroid_index_array(1)) exit
                    
                    elseif (prev_face .eq. face_2) then
                    
                        centroid_index_array(i)   = current_edge
                        centroid_index_array(i+1) = face_1
                    
                        centroid_array(i,:)   = e_centroid(current_edge,:)
                        centroid_array(i+1,:) = f_centroid(face_1,:)
                    
                        centroid_array_count_real = centroid_array_count_real + 2
                        
                        non_viable_edges(pe) = .true.
                    
                        prev_face = centroid_index_array(i+1)
                        i=i+2 
                    
                        !print*, 'swapped', i, centroid_array_count+1
                    
                        if (face_1 .eq. centroid_index_array(1)) exit
                    
                    else
                        !print*, 'none, looping', i, centroid_array_count+1, edge_count, face_count
                        
                    endif
                
                    !print*, pe, current_edge, face_1, face_2 , prev_face
                    
                    !print*, 'centroid array ', centroid_index_array
                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                enddo
            
               ! print*, 'centroid array ', centroid_index_array(1:centroid_array_count_real)
                
            endif
            centroid_array_count_old = centroid_array_count
            
            
            !if (i.ne.centroid_array_count+1) print*, 'centroid array count broken'
            !if (centroid_index_array(1).ne.centroid_index_array(centroid_array_count)) print*, 'internal centroid array start and end wrong'
            
            
            in_progress_projection(:) = 0.0
            in_progress_centroid = coords(p,:)
            ! by making i1, i2, i3 all p, the volume change should be zero
            
            
            c1(:) = coords(i1,:)
            c2(:) = coords(i2,:)
            
            call centroid_array_routine(in_progress_projection,  in_progress_centroid, centroid_array_count_real, centroid_array)
            
            direction_array(:) = p_normal_vectors(bp,:)
            
            angle = alignment(in_progress_projection, direction_array)
            
            if (angle .lt. 0.) then
                in_progress_projection(:) = in_progress_projection(:)
            else
                in_progress_projection(:) = -in_progress_projection(:)
            endif
            
            sbb(bp,:) = in_progress_projection
            
        enddo
        
        !do i=1,b_npoin
        !    print*, i, sqrt(sbb(i,1)*sbb(i,1) + sbb(i,2)*sbb(i,2) + sbb(i,3)*sbb(i,3)), sbb(i,:)
        !enddo
        !print*, 'aaa'
        
        print*, 'BACKSTOP IN BOUND FACE VOLUME PROCESSING'
        stop
        
    end subroutine boundary_face_volume_processing

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine centroid_assembler(obj, obj_type, centroid_obj_array)
        implicit none
        integer(KIND=INT32) :: obj, obj_type
        integer(KIND=INT32) :: mm, tt(2), i, j, k, m_i, m_stt, m_end, t_stt, t_end, fwd_i, bck_i
        integer(KIND=INT32),parameter :: internal_edge=1, boundary_edge=2, featre_point=3, non_feature_point=4
        integer(KIND=INT32) :: main_count, tert_count, centroid_array_count
        integer(KIND=INT32),allocatable :: main_obj(:), tert_obj(:), centroid_obj_array(:)
        !real(KIND=REAL64)  ,allocatable :: centroid_array(:,3)
        logical, allocatable :: viable_mains(:)
        logical ::  e_state, state
        
        deallocate(centroid_obj_array)
        
        select case(obj_type)
            case(internal_edge)
                
                ! faces are main
                ! cells are tertiary
                m_stt   = e_f_index_array(obj-1)
                t_stt   = e_c_index_array(obj-1)
                m_end   = e_f_index_array(obj)
                t_end   = e_c_index_array(obj)
                
            main_count = m_end-m_stt
            tert_count = t_end-t_stt
            
            ! for i_edge the centroid will consist of all adjacent centroids plus one repeat
            centroid_array_count = 1 + main_count + tert_count
            
            ! inital edge setup doesnt matter much for i_edge
            ! just pick a random edge, my end condition is whenever the first cell is equal to the last cell
            m_i = 1
            
            case(boundary_edge)
                
                
                
                
            case(non_feature_point)
                
                
                
                
            case(featre_point)
                
                
                
                
        end select
        
        ! allocate my temp arrays
        allocate(viable_mains(main_count), centroid_obj_array(centroid_array_count))
        
        ! all main objects are viable to begin with.
        viable_mains = .true.
        
        call obj_select(m_i, m_stt, mm, tt, obj_type)
        
        ! we just assign our values into the first available slots
        centroid_obj_array(1) = tt(1)
        centroid_obj_array(2) = mm
        centroid_obj_array(3) = tt(2)
        
        ! and initialise our indexers
        fwd_i = 3
        bck_i = 1
        
        do
            ! now we loop about our centroids
            if (m_i.eq.main_count) then
                m_i = 1
            else
                m_i=m_i+1
            endif
            
            ! if we know the centroid isnt viable we skip it
            if (.not.viable_mains(m_i)) cycle
            
            ! obviously we get our objects
            call obj_select(m_i, m_stt, mm, tt, obj_type)
            
!                     print*, fwd_i, bck_i
            
            call centroid_swapper(mm, tt, fwd_i, bck_i, centroid_obj_array, state)
            
!                     print*, tt(1), mm, tt(2)
!                     print*, centroid_obj_array
!                     print*, fwd_i, bck_i, state
!                     print*, ''
!                     print*, ''
            
            if (state) then
                
                viable_mains(m_i) = .false.
                
                call obj_endconditions(e_state, fwd_i, bck_i, centroid_obj_array, centroid_array_count, obj_type)
                ! our end condition is the first edge and last edge being the same
                if (e_state) exit
            endif
            
            
        enddo
        
        print*, centroid_obj_array
        
    end subroutine centroid_assembler
    
    subroutine obj_select(m_i, m_stt, mm, tt, obj_type)
        implicit none
        integer(KIND=INT32) :: obj_type
        integer(KIND=INT32) :: m_i, m_stt, mm, tt(2)
        integer(KIND=INT32),parameter :: internal_edge=1, boundary_edge=2, featre_point=3, non_feature_point=4
        
        
        select case(obj_type)
            case(internal_edge)
                mm    = e_f_obj_relation_array(m_stt + m_i)
                tt(1) = f_c_obj_relation_array(f_c_index_array(mm)    )
                tt(2) = f_c_obj_relation_array(f_c_index_array(mm-1)+1)
                
            case(boundary_edge)
                
                
                
                
            case(non_feature_point)
                
                
                
                
            case(featre_point)
            
        end select
            
    end subroutine obj_select
    
    subroutine obj_endconditions(e_state, fwd_i, bck_i, centroid_obj_array, centroid_array_count, obj_type)
        implicit none
        integer(KIND=INT32) :: fwd_i, bck_i, centroid_array_count
        integer(KIND=INT32) :: centroid_obj_array(:)
        integer(KIND=INT32),parameter :: internal_edge=1, boundary_edge=2, featre_point=3, non_feature_point=4
        logical ::  e_state
        
        e_state = .false.
        select case(obj_type)
            
            case(internal_edge)
            
                if (centroid_obj_array(1) .eq. centroid_obj_array(centroid_array_count)) e_state = .true.
                
            case(boundary_edge)
                
                
                
                
            case(non_feature_point)
                
                
                
                
            case(featre_point)
            
        end select
            
    end subroutine obj_endconditions
    
    subroutine centroid_swapper(mm, tt, fwd_i, bck_i, c_array, state)
        implicit none
        integer(KIND=INT32) :: mm, tt(2), jj(2), fwd_i, bck_i
        integer(KIND=INT32) :: c_array(:)
        logical ::  state
        
        state = .false.
        
        ! jj are the end values
        jj(1) = c_array(fwd_i)
        jj(2) = c_array(bck_i)
        
        ! there are four worlds 
        ! two tertiary objects and two end objects
        ! so four conditions
        ! imagine tt(1) is the left value of the couple tt(1), mm, tt(2)
        ! if tt(1) = jj(1) then our centroid array should look like jj, mm, tt(2) No swap go forward
        ! if tt(2) = jj(1) then our centroid array should look like jj, mm, tt(1) swap tt go forward
        
        ! if tt(1) = jj(2) then our centroid array should look like tt(2), mm, jj swap tt go backward
        ! if tt(2) = jj(2) then our centroid array should look like tt(1), mm, jj no swap go backward
        if    (tt(1).eq.jj(1)) then
            
            call centroid_shuffle(c_array, fwd_i, bck_i,  2)
            
            c_array(fwd_i+1)   = mm
            c_array(fwd_i+2) = tt(2)
            
            state = .true.
            fwd_i = fwd_i + 2
            
        elseif(tt(2).eq.jj(1)) then
            
            call centroid_shuffle(c_array, fwd_i, bck_i,  2)
            
            c_array(fwd_i+1)   = mm
            c_array(fwd_i+2) = tt(1)
            
            state = .true.
            fwd_i = fwd_i + 2
            
        elseif(tt(1).eq.jj(2)) then
            
            call centroid_shuffle(c_array, fwd_i, bck_i, -2)
            
            c_array(bck_i-1)   = mm
            c_array(bck_i-2) = tt(2)
            
            state = .true.
            bck_i = bck_i - 2
            
        elseif(tt(2).eq.jj(2)) then
            
            call centroid_shuffle(c_array, fwd_i, bck_i, -2)
            
            c_array(bck_i-1)   = mm
            c_array(bck_i-2) = tt(1)
            
            state = .true.
            bck_i = bck_i - 2
            
        endif
        
    end subroutine centroid_swapper
    
    subroutine centroid_shuffle(c_array, i, j, change)
        implicit none
        integer(KIND=INT32) :: change, i, j, c_size, k
        integer(KIND=INT32) :: c_array(:)
        
        c_size = size(c_array)
        
        if (change.gt.0)then
        
            if (i + change .gt. c_size) then
                do k=j, i
                    c_array(k-change) = c_array(k)
                enddo
                
                i = i - change
                j = j - change
                
            endif
            
        elseif (change.lt.0)then
        
            if (j + change .lt. 0) then
                do k=i, j, -1
                    c_array(k-change) = c_array(k)
                enddo
                
                i = i - change
                j = j - change
                
            endif
            
        endif
        
    end subroutine centroid_shuffle
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine projection_test
        implicit none
        integer(KIND=INT32)           :: i, i1, i2
        real(KIND=REAL64)			  :: total_volume
        real(KIND=REAL64),allocatable :: projection_mag(:), tot(:,:)
        
        
        allocate(tot(npoin,3))
        tot = 0.0
        do i=1,i_nedge
            i1 = e_p_obj_relation_array(e_p_index_array(e_internal_indexing_array(i))-1)
            i2 = e_p_obj_relation_array(e_p_index_array(e_internal_indexing_array(i)))
            
            tot(i1,:) = tot(i1,:) + sn(i,:)
            tot(i2,:) = tot(i2,:) - sn(i,:)
        enddo
        do i=1,b_nedge
            i1 = e_p_obj_relation_array(e_p_index_array(e_bound_indexing_array(i))-1)
            i2 = e_p_obj_relation_array(e_p_index_array(e_bound_indexing_array(i)))
            
            tot(i1,:) = tot(i1,:) + sb(i,:)
            tot(i2,:) = tot(i2,:) - sb(i,:)
        enddo
        do i=1,b_npoin
            i1 = p_bound_indexing_array(i)
            
            tot(i1,:) = tot(i1,:) + sbb(i,:)
        enddo
        
        print*, ' testing projections '
        
        allocate(projection_mag(i_nedge))
        projection_mag(:) = sqrt(sn(:,1)*sn(:,1) + sn(:,2)*sn(:,2) + sn(:,3)*sn(:,3))
        !print*, projection_mag
        i1=1
        i2=1
        do i=1,i_nedge
            if (projection_mag(i).lt.projection_mag(i1)) i1 = i
            if (projection_mag(i).gt.projection_mag(i2)) i2 = i
            
        enddo
        print'(a,4e14.4e3,i,a,4e14.4e3,i)', 'lowest sn  ', projection_mag(i1), sn(i1,:), i1, ' highest sn  ', projection_mag(i2), sn(i2,:), i2
        print*, ' '
        
        deallocate(projection_mag)
        allocate(projection_mag(b_nedge))
        projection_mag(:) = sqrt(sb(:,1)*sb(:,1) + sb(:,2)*sb(:,2) + sb(:,3)*sb(:,3))
        !print*, projection_mag
        i1=1
        i2=1
        do i=1,b_nedge
            if (projection_mag(i).lt.projection_mag(i1)) i1 = i
            if (projection_mag(i).gt.projection_mag(i2)) i2 = i
            
        enddo
        print'(a,4e14.4e3,i,a,4e14.4e3,i)', 'lowest sb  ', projection_mag(i1), sb(i1,:), i1, ' highest sb  ', projection_mag(i2), sb(i2,:), i2
        print*, ' '
        
        deallocate(projection_mag)
        allocate(projection_mag(b_npoin))
        projection_mag(:) = sqrt(sbb(:,1)*sbb(:,1) + sbb(:,2)*sbb(:,2) + sbb(:,3)*sbb(:,3))
        !print*, projection_mag
        i1=1
        i2=1
        do i=1,b_npoin
            if (projection_mag(i).lt.projection_mag(i1)) i1 = i
            if (projection_mag(i).gt.projection_mag(i2)) i2 = i
            
        enddo
        print'(a,4e14.4e3,i,a,4e14.4e3,i)', 'lowest sbb ', projection_mag(i1), sbb(i1,:), i1, ' highest sbb ', projection_mag(i2), sbb(i2,:), i2
        print*, ' '
        
        deallocate(projection_mag)
        allocate(projection_mag(npoin))
        projection_mag(:) = sqrt(tot(:,1)*tot(:,1) + tot(:,2)*tot(:,2) + tot(:,3)*tot(:,3))
        !print*, projection_mag
        i1=1
        i2=1
        do i=1,npoin
            if (projection_mag(i).lt.projection_mag(i1)) i1 = i
            if (projection_mag(i).gt.projection_mag(i2)) i2 = i
            
        enddo
        print'(a,4e14.4e3,i,a,4e14.4e3,i)', 'lowest tot ', projection_mag(i1), tot(i1,:), i1, ' highest tot ', projection_mag(i2), tot(i2,:), i2
        print*, ' '
        
        total_volume = 0.0
        i1=1
        i2=1
        do i=1,npoin
            if (vol(i).lt.vol(i1)) i1 = i
            if (vol(i).gt.vol(i1)) i2 = i
            total_volume = total_volume + vol(i)
        enddo
        
        print'(a,e14.4e3,i,a,e14.4e3,i,a,f)', 'lowest vol ', vol(i1), i1, ' highest vol ', vol(i2), i2, ' total volume ', total_volume
        print*, ' '
        
        print*, MAXVAL(coords(:,1)),MINVAL(coords(:,1)),MAXVAL(coords(:,2)),MINVAL(coords(:,2)),MAXVAL(coords(:,3)),MINVAL(coords(:,3))
        
    end subroutine projection_test
    
    subroutine centroid_array_routine(obj_projection, obj_centroid, centroid_array_count, centroid_array)
        implicit none
        integer(KIND=INT32)        :: i, centroid_array_count
        real(KIND=REAL64)          :: obj_projection(3), centroid_array(:,:), obj_centroid(3), baryobj_1_centroid(3), baryobj_2_centroid(3)
        

        do i=2,centroid_array_count
            baryobj_1_centroid = centroid_array(i-1,:)
            baryobj_2_centroid = centroid_array(i,:) 
            call cvolume(obj_projection, obj_centroid, baryobj_1_centroid, baryobj_2_centroid)
        enddo

    end subroutine centroid_array_routine
    
    
!     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine volume_generation
        implicit none
        integer(KIND=INT32)        :: i, i1, i2, e, ie, be
        real(KIND=REAL64)          :: vol1, vol2, c1(3), c2(3), c3(3), sxx(3)
        
        allocate(vol(npoin))
        vol = 0.0
        
        do ie = 1, i_nedge
            e = e_internal_indexing_array(ie)
            ! this is a loop of all internal edges.
            
            i1 = e_p_obj_relation_array(e_p_index_array(e)-1)
            i2 = e_p_obj_relation_array(e_p_index_array(e))
            
            c1 = coords(i1,:)
            c2 = coords(i2,:)
            c3 = e_centroid(e,:)
            
            sxx(:) = sn(ie,:)
            
            call volume_from_proj(sxx, c1, c2, c3, vol1, vol2)
            
            vol(i1)  = vol(i1) + vol1
            vol(i2)  = vol(i2) + vol2
            
        enddo
        
        do be = 1, b_nedge
            e = e_bound_indexing_array(be)
            ! this is a loop of all boundary edges.
            
            i1 = e_p_obj_relation_array(e_p_index_array(e)-1)
            i2 = e_p_obj_relation_array(e_p_index_array(e))
            
            c1 = coords(i1,:)
            c2 = coords(i2,:)
            c3 = e_centroid(e,:)
            
            sxx(:) = sb(be,:)
            
            call volume_from_proj(sxx, c1, c2, c3, vol1, vol2)
            
            vol(i1)  = vol(i1) + vol1
            vol(i2)  = vol(i2) + vol2
            
        enddo
        
    end subroutine volume_generation
    
    subroutine volume_from_proj(sx, i1, i2, i3, vol1, vol2)
        implicit none
        real(KIND=REAL64)          :: vol1, vol2
        real(KIND=REAL64)          :: sx(3), i1(3), i2(3), i3(3), v31(3), v32(3)
        
        v31    = i3 - i1
        v32    = i3 - i2
        vol1   = abs((sx(1)*v31(1)) + (sx(2)*v31(2)) + (sx(3)*v31(3)) )/3
        vol2   = abs((sx(1)*v32(1)) + (sx(2)*v32(2)) + (sx(3)*v32(3)) )/3
        
    end subroutine volume_from_proj
    
    
    subroutine cvolume(sx, i3, i4, i5)
        implicit none
        real(KIND=REAL64),dimension(3) :: sx, i3, i4, i5
        real(KIND=REAL64),dimension(3) :: v34,v35
        real(KIND=REAL64),dimension(3) :: c3435
        
        integer(KIND=INT32) :: i
        

        do i=1,3
            v34(i) = i3(i) - i4(i) ! vector from 3 to 4
            v35(i) = i3(i) - i5(i) ! vector from 3 to 5
        enddo
        
        ! Volume=∥a×b∥ ∥c∥ |cosϕ|=|(a×b)⋅c|. volume of parallelepiped of sides abc
        ! Volume = Volume/6 think a parallelepiped is 6 tets
        !https://mathinsight.org/scalar_triple_product
        ! a = v14 b = v15
        
        c3435(1) = (v34(2) * v35(3) - v34(3) * v35(2))/2
        c3435(2) = (v34(3) * v35(1) - v34(1) * v35(3))/2
        c3435(3) = (v34(1) * v35(2) - v34(2) * v35(1))/2

        sx(1) = sx(1) + c3435(1) ! projection in the xx axis
        sx(2) = sx(2) + c3435(2) ! projection in the yy axis
        sx(3) = sx(3) + c3435(3) ! projection in the zz axis
        
        
    !                               i4--------------i5
    !                                \  >>>>>>>   /
    !                                 \         /
    !                                  \      / 
    !                                   \   /
    !           i1 -------------------- i3 -------------------------i2
    end subroutine cvolume
    
    ! end contains
    
end module volume_processing
