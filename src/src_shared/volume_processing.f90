

module volume_processing
    use iso_fortran_env
    use object_counts
    use raw_data
    use object_relation_data
    use centroid_data
    use boundary_data
    use quicksort_module
    use projection_data
    ! need everything p much here.
    implicit none
    
        
    contains
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine volume_alloc
        implicit none
        allocate(vol(npoin))
        
    end subroutine volume_alloc
        
    subroutine internal_edge_volume_processing
        implicit none
        integer(KIND=INT32) :: ie, e, i, j, i1, i2, centroid_array_count, centroid_array_count_old, ef, ec, ef_start, ec_start, ef_end, ec_end
        integer(KIND=INT32) :: cell_count, face_count, current_face, cell_1, cell_2, prev_cell
        integer(KIND=INT32),allocatable :: centroid_index_array(:)
        real(KIND=REAL32),allocatable   :: centroid_array(:,:)
        real(KIND=REAL32)               :: in_progress_projection(3), in_progress_centroid(3)
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
        
        allocate(sn(i_nedge,3))
        
        centroid_array_count_old = -1
        allocate(centroid_index_array(0),non_viable_faces(0),centroid_array(0,0))
        
        do ie = 1, i_nedge
            e = e_internal_indexing_array(ie)
            ! this is a loop of all internal edges.
            
            i1 = e_p_obj_relation_array(e_p_index_array(e)-1)
            i2 = e_p_obj_relation_array(e_p_index_array(e))
            ! i1 and i2 are the constituent points of edge e
            
            ef_start = e_f_index_array(e-1)
            ec_start = e_c_index_array(e-1)
            ef_end   = e_f_index_array(e)
            ec_end   = e_c_index_array(e)
            
            cell_count = ec_end-ec_start
            face_count = ef_end-ef_start
            
            centroid_array_count = 1 + cell_count + face_count
            
            ! sum of number of faces, number of edges plus 1 for the duplicate starting edge
            
            if (centroid_array_count_old .ne. centroid_array_count) then
                deallocate(centroid_index_array,non_viable_faces,centroid_array)
                allocate(centroid_index_array(centroid_array_count), non_viable_faces(face_count+1), centroid_array(centroid_array_count,3))
            endif
            
            non_viable_faces = .false.
            
            current_face = e_f_obj_relation_array(ef_start+1)
            cell_1 = f_c_obj_relation_array(f_c_index_array(current_face)    )
            cell_2 = f_c_obj_relation_array(f_c_index_array(current_face-1)+1)
            
            centroid_index_array(1) = cell_1
            centroid_index_array(2) = current_face
            centroid_index_array(3) = cell_2
            
            centroid_array(1,:) = c_centroid(cell_1,:)
            centroid_array(2,:) = f_centroid(current_face,:)
            centroid_array(3,:) = c_centroid(cell_2,:)
                    
            prev_cell = cell_2
            
            !print*, 'number of faces ', face_count, ' number of cells ', cell_count, ' count ', centroid_array_count
            
            i=4
            
            ef = 1
            
            !print*, ef, current_face, cell_1, cell_2 ,prev_cell, 'prev_cell'
            
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
                    centroid_index_array(i+1) = cell_2
                    
                    centroid_array(i,:)   = f_centroid(current_face,:)
                    centroid_array(i+1,:) = c_centroid(cell_2,:)
                    
                    non_viable_faces(ef) = .true.
                    
                    prev_cell = centroid_index_array(i+1)
                    i=i+2 
                    
                    !print*, 'i-1 = cell_1, i+1 = cell_2'
                    
                    if (i.eq.centroid_array_count+1) exit
                    
                elseif (prev_cell .eq. cell_2) then
                    
                    centroid_index_array(i)   = current_face
                    centroid_index_array(i+1) = cell_1
                    
                    centroid_array(i,:)   = f_centroid(current_face,:)
                    centroid_array(i+1,:) = c_centroid(cell_1,:)
                    
                    non_viable_faces(ef) = .true.
                    
                    prev_cell = centroid_index_array(i+1)
                    i=i+2 
                    
                    !print*, 'swapped'
                    
                    if (i.eq.centroid_array_count+1) exit
                    
                else
                    !print*, 'none, looping'
                endif
                
                !print*, ef, current_face, cell_1, cell_2 , prev_cell
                
            enddo
            
            centroid_index_array(centroid_array_count) = centroid_index_array(3)
            ! this might be exactly negative, alternatives are 1, and cell_1
            centroid_array(centroid_array_count,:) = c_centroid(cell_2,:)
            
            
            centroid_array_count_old = centroid_array_count
            
            
            !if (i.ne.centroid_array_count+1) print*, 'centroid array count broken'
            if (centroid_index_array(1).ne.centroid_index_array(centroid_array_count)) print*, 'internal centroid array start and end wrong'
            
            !print*, centroid_index_array
            
            in_progress_projection(:) = 0.0
            in_progress_centroid = e_centroid(e,:)
            
            call centroid_array_routine(in_progress_projection, in_progress_centroid, centroid_array_count, centroid_array, i1, i2)
            
            sn(e,:) = in_progress_projection
            
        enddo
        
        print*,'aaaa'
        do i=1,i_nedge
            print*, i, sqrt(sn(i,1)*sn(i,1) + sn(i,2)*sn(i,2) + sn(i,3)*sn(i,3)), sn(i,:)
        enddo
        print*,'aaaa'
        
    end subroutine internal_edge_volume_processing
    
        
    subroutine boundary_edge_volume_processing
        implicit none
        integer(KIND=INT32) :: be, e, i, j, i1, i2, centroid_array_count, centroid_array_count_old, ef, ec, ef_start, ec_start, ef_end, ec_end
        integer(KIND=INT32) :: cell_count, face_count, current_face, cell_1, cell_2, prev_cell
        integer(KIND=INT32),allocatable :: centroid_index_array(:)
        real(KIND=REAL32),allocatable   :: centroid_array(:,:)
        real(KIND=REAL32)               :: in_progress_projection(3), in_progress_centroid(3)
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
        
        allocate(sb(b_nedge,3))
        
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
                if (f_bound_array(current_face)) exit
                ef = ef + 1
            enddo

            non_viable_faces(ef) = .true.
            
            cell_1 = f_c_obj_relation_array(f_c_index_array(current_face)    )
            
            centroid_index_array(1) = current_face
            centroid_index_array(2) = cell_1
            
            centroid_array(1,:) = f_centroid(current_face,:)
            centroid_array(2,:) = c_centroid(cell_1,:)
                    
            prev_cell = cell_1
            
            !print*, 'number of faces ', face_count, ' number of cells ', cell_count, ' count ', centroid_array_count
            
            i=3
            
            !print*, ef, current_face, cell_1, cell_2 ,prev_cell, 'prev_cell'
            
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
                
            enddo
            
            
            centroid_array_count_old = centroid_array_count
            
            
            !if (i.ne.centroid_array_count+1) print*, 'centroid array count broken'
            !if (centroid_index_array(1).ne.centroid_index_array(centroid_array_count)) print*, 'internal centroid array start and end wrong'
            
            
            in_progress_projection(:) = 0.0
            in_progress_centroid = e_centroid(e,:)
            
            
            call centroid_array_routine(in_progress_projection, in_progress_centroid, centroid_array_count, centroid_array, i1, i2)
            
            sb(be,:) = in_progress_projection
            
        enddo
        
        do i=1,b_nedge
            print*, i, sqrt(sb(i,1)*sb(i,1) + sb(i,2)*sb(i,2) + sb(i,3)*sb(i,3)), sb(i,:)
        enddo
        print*,'aaaa'
        
    end subroutine boundary_edge_volume_processing
    
    
    subroutine boundary_face_volume_processing
        implicit none
        integer(KIND=INT32) :: bp, p, i, j, i1, i2, centroid_array_count, centroid_array_count_old, pe, pe_start, pf_start, pe_end, pf_end
        integer(KIND=INT32) :: edge_count, face_count, current_edge, face_1, face_2, prev_face
        integer(KIND=INT32),allocatable :: centroid_index_array(:)
        real(KIND=REAL32),allocatable   :: centroid_array(:,:)
        real(KIND=REAL32)               :: in_progress_projection(3), in_progress_centroid(3)
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
        
        allocate(sbb(b_npoin,3))
        
        centroid_array_count_old = -1
        allocate(centroid_index_array(0),non_viable_edges(0),centroid_array(0,0))
        
        do bp = 1, b_npoin
            p = p_bound_indexing_array(bp)
            ! this is a loop of all boundary edges.
            
            print*, bp, p, 'inside loop'
            
            i1 = p
            i2 = p
            ! i1 and i2 are the constituent points of edge e
            
            pe_start = p_e_index_array(p-1)
            pf_start = p_f_index_array(p-1)
            pe_end   = p_e_index_array(p)
            pf_end   = p_f_index_array(p)
            
            edge_count = pe_end-pe_start
            face_count = pf_end-pf_start
            
            if (feature_points(bp)) then
                
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! FEATURE POINTS
                ! IE CORNERS
                
                print*, bp, p, 'inside feature statement'
                
                centroid_array_count = edge_count + face_count
            
                ! sum of number of faces, number of edges plus 1 for the duplicate starting edge
            
                if (centroid_array_count_old .ne. centroid_array_count) then
                    deallocate(centroid_index_array,non_viable_edges,centroid_array)
                    allocate(centroid_index_array(centroid_array_count), non_viable_edges(edge_count), centroid_array(centroid_array_count,3))
                endif
            
                pe = 1
            
                non_viable_edges = .false.
                
                do ! we must start on a feature edge
                    current_edge = p_e_obj_relation_array(pe_start+pe)
                    if (e_internal_array(current_edge)) non_viable_edges(pe) = .true.
                    if (feature_edges(reversed_e_bound_indexing_array(current_edge))) exit
                    pe = pe + 1
                enddo
                
                non_viable_edges(pe) = .true.
                
                j = e_f_index_array(current_edge-1)
                do 
                    j=j+1
                    if (f_bound_array(e_f_obj_relation_array(j))) exit
                enddo
                face_1 = e_f_obj_relation_array(j)
                face_2 = -1
                
                if (j.gt.e_f_index_array(current_edge)) print*, ' feature boundary edge has more than 1 faces '
                
                centroid_index_array(1) = current_edge
                centroid_index_array(2) = face_1
            
                centroid_array(1,:) = e_centroid(current_edge,:)
                centroid_array(2,:) = f_centroid(face_1,:)
                    
                prev_face = face_1
            
                !print*, 'number of faces ', face_count, ' number of cells ', cell_count, ' count ', centroid_array_count
            
                i=3
            
                !print*, ef, current_face, cell_1, cell_2 ,prev_cell, 'prev_cell'
            
                do
                
                    if (pe.eq.edge_count) then
                        pe = 1
                    else
                        pe=pe+1
                    endif

                
                    if (non_viable_edges(pe)) cycle
                
                    do ! we must only use boundary edges
                        current_edge = p_e_obj_relation_array(pe_start+pe)
                        if (e_bound_array(current_edge)) exit
                        non_viable_edges(pe) = .true.
                        pe = pe + 1
                    enddo
                    
                    j = e_f_index_array(current_edge-1)
                    do 
                        j=j+1
                        if (f_bound_array(e_f_obj_relation_array(j))) exit
                    enddo
                    face_1 = e_f_obj_relation_array(j)
                    face_2 = -1
                    
                    if (j.gt.e_f_index_array(current_edge)) print*, ' feature boundary edge has more than 1 faces '
                    
                    if (prev_face .eq. face_1) then
            
                        centroid_index_array(i)   = current_edge
                        centroid_array(i,:)   = e_centroid(current_edge,:)
                    
                        if (feature_edges(reversed_e_bound_indexing_array(current_edge))) exit  
                    
                        centroid_index_array(i+1) = face_2
                        centroid_array(i+1,:) = f_centroid(face_2,:)
                    
                        non_viable_edges(pe) = .true.
                    
                        prev_face = centroid_index_array(i+1)
                        i=i+2 
                    
                        print*, 'i-1 = face_1, i+1 = face_2'
                    
                    elseif (prev_face .eq. face_2) then
                    
                        centroid_index_array(i)   = current_edge
                        centroid_array(i,:)   = e_centroid(current_edge,:)
                    
                        if (feature_edges(reversed_e_bound_indexing_array(current_edge))) exit
                    
                        centroid_index_array(i+1) = face_1
                        centroid_array(i+1,:) = f_centroid(face_1,:)
                    
                        non_viable_edges(pe) = .true.
                    
                        prev_face = centroid_index_array(i+1)
                        i=i+2 
                    
                        print*, 'swapped'
                    
                    else
                        !print*, 'none, looping'
                    endif
                
                    !print*, ef, current_edge, face_1, face_2 , prev_face
                
                enddo
            else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                ! NON FEATURE POINTS
                ! IE NOT CORNERS
                
                centroid_array_count = 1 + edge_count + face_count
            
                ! sum of number of faces, number of edges plus 1 for the duplicate starting edge
            
                if (centroid_array_count_old .ne. centroid_array_count) then
                    deallocate(centroid_index_array,non_viable_edges,centroid_array)
                    allocate(centroid_index_array(centroid_array_count), non_viable_edges(edge_count), centroid_array(centroid_array_count,3))
                else
                    centroid_index_array = 0
                    centroid_array = 0.0
                endif
            
                pe = 1
                
                non_viable_edges = .false.
            
                do ! we must only use boundary edges
                    current_edge = p_e_obj_relation_array(pe_start+pe)
                    if (e_bound_array(current_edge)) exit
                    non_viable_edges(pe) = .true.
                    pe = pe + 1
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
                
                if (j.gt.e_f_index_array(current_edge)) print*, ' non-feature boundary edge has more than 2 faces '
            
                centroid_index_array(1) = face_1
                centroid_index_array(2) = current_edge
                centroid_index_array(3) = face_2
            
                centroid_array(1,:) = f_centroid(face_1,:)
                centroid_array(2,:) = e_centroid(current_edge,:)
                centroid_array(3,:) = f_centroid(face_2,:)
                
                prev_face = face_2
            
                !print*, 'number of faces ', face_count, ' number of cells ', cell_count, ' count ', centroid_array_count
            
                i=4
            
                print*, pe, current_edge, face_1, face_2 , prev_face
                
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
                        pe = pe + 1
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
                    
                    if (j.gt.e_f_index_array(current_edge)) print*, ' non-feature boundary edge has more than 2 faces '
                
                    if (prev_face .eq. face_1) then
            
                        centroid_index_array(i)   = current_edge
                        centroid_index_array(i+1) = face_2
                    
                        centroid_array(i,:)   = e_centroid(current_edge,:)
                        centroid_array(i+1,:) = f_centroid(face_2,:)
                    
                        non_viable_edges(pe) = .true.
                    
                        prev_face = centroid_index_array(i+1)
                        i=i+2 
                    
                        print*, 'i-1 = face_1, i+1 = face_2', i, centroid_array_count+1
                    
                        if (face_2 .eq. centroid_index_array(1)) exit
                    
                    elseif (prev_face .eq. face_2) then
                    
                        centroid_index_array(i)   = current_edge
                        centroid_index_array(i+1) = face_1
                    
                        centroid_array(i,:)   = e_centroid(current_edge,:)
                        centroid_array(i+1,:) = f_centroid(face_1,:)
                    
                        non_viable_edges(pe) = .true.
                    
                        prev_face = centroid_index_array(i+1)
                        i=i+2 
                    
                        print*, 'swapped', i, centroid_array_count+1
                    
                        if (face_1 .eq. centroid_index_array(1)) exit
                    
                    else
                        print*, 'none, looping', i, centroid_array_count+1, edge_count, face_count
                        
                    endif
                
                    print*, pe, current_edge, face_1, face_2 , prev_face
                    
                    !print*, 'centroid array ', centroid_index_array
                    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                enddo
            
                print*, 'centroid array ', centroid_index_array
                
            endif
            centroid_array_count_old = centroid_array_count
            
            
            !if (i.ne.centroid_array_count+1) print*, 'centroid array count broken'
            !if (centroid_index_array(1).ne.centroid_index_array(centroid_array_count)) print*, 'internal centroid array start and end wrong'
            
            
            in_progress_projection(:) = 0.0
            in_progress_centroid = coords(p,:)
            ! by making i1, i2, i3 all p, the volume change should be zero
            
            call centroid_array_routine(in_progress_projection, in_progress_centroid, centroid_array_count, centroid_array, i1, i2)
            
            sbb(bp,:) = in_progress_projection
            
        enddo
        
        do i=1,b_npoin
            print*, i, sqrt(sbb(i,1)*sbb(i,1) + sbb(i,2)*sbb(i,2) + sbb(i,3)*sbb(i,3)), sbb(i,:)
        enddo
        
    end subroutine boundary_face_volume_processing
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine centroid_array_routine(obj_projection, obj_centroid, centroid_array_count, centroid_array, i1, i2)
        implicit none
        integer(KIND=INT32)          :: i, centroid_array_count, i1
        integer(KIND=INT32),optional :: i2
        real(KIND=REAL32) :: obj_projection(3), centroid_array(:,:), obj_centroid(3), baryobj_1_centroid(3), baryobj_2_centroid(3)
        
        do i=2,centroid_array_count
            baryobj_1_centroid = centroid_array(i-1,:)
            baryobj_2_centroid = centroid_array(i,:) 
            call cvolume(obj_projection, vol(i1), vol(i2), i1, i2, obj_centroid, baryobj_1_centroid, baryobj_2_centroid)
        enddo
        
        
    end subroutine
    
    subroutine cvolume(sn, vol1, vol2, i1, i2, i3, i4, i5)
        implicit none
        real(KIND=REAL32) :: vol1,vol2
        real(KIND=REAL32),dimension(3) :: v13,v14,v15,v32,v42,v52,c1415,c4252
        real(KIND=REAL32),dimension(3) :: sn, i3, i4, i5
        real(KIND=REAL32),dimension(3) :: v34,v35
        integer(KIND=INT32) :: i, i1,i2
        
        do i=1,3
            v13(i) = coords(i1,i) - i3(i) ! vector from 1 to 3
            v14(i) = coords(i1,i) - i4(i) ! vector from 1 to 4
            v15(i) = coords(i1,i) - i5(i) ! vector from 1 to 5
            v32(i) = i3(i) - coords(i2,i) ! vector from 3 to 2
            v42(i) = i4(i) - coords(i2,i) ! vector from 4 to 2
            v52(i) = i5(i) - coords(i2,i) ! vector from 5 to 2
        enddo
        ! Volume=∥a×b∥ ∥c∥ |cosϕ|=|(a×b)⋅c|. volume of parallelepiped of sides abc
        ! Volume = Volume/6 think a parallelepiped is 6 tets
        !https://mathinsight.org/scalar_triple_product
        ! a = v14 b = v15
        c1415(1) = v14(2) * v15(3) - v14(3) * v15(2)
        c1415(2) = v14(3) * v15(1) - v14(1) * v15(3)
        c1415(3) = v14(1) * v15(2) - v14(2) * v15(1)
        c4252(1) = v42(2) * v52(3) - v42(3) * v52(2)
        c4252(2) = v42(3) * v52(1) - v42(1) * v52(3)
        c4252(3) = v42(1) * v52(2) - v42(2) * v52(1)
        
        vol1 = vol1 + abs((c1415(1)*v13(1) ) - (c1415(2)*v13(2) ) + (c1415(3)*v13(3) ) )/6
        vol2 = vol2 + abs((c1415(1)*v32(1) ) - (c1415(2)*v32(2) ) + (c1415(3)*v32(3) ) )/6

        do i=1,3
            v34(i) = i3(i) - i4(i) ! vector from 3 to 4
            v35(i) = i3(i) - i5(i) ! vector from 3 to 5
        enddo

        sn(1) = sn(1) + ((v34(2) * v35(3) - v34(3) * v35(2))/2) ! projection in the xx axis
        sn(2) = sn(2) + ((v34(3) * v35(1) - v34(1) * v35(3))/2) ! projection in the yy axis
        sn(3) = sn(3) + ((v34(1) * v35(2) - v34(2) * v35(1))/2) ! projection in the zz axis
        
    !                               i4--------------i5
    !                                \  >>>>>>>   /
    !                                 \         /
    !                                  \      / 
    !                                   \   /
    !           i1 -------------------- i3 -------------------------i2
    end subroutine cvolume
    
    ! end contains
    
end module volume_processing
