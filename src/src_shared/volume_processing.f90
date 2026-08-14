

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
    integer(KIND=INT32),parameter :: internal_edge=1, boundary_edge=2, featre_point=3, non_feature_point=4
    integer(KIND=INT32) :: obj_type, flag
    
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
        
        obj_type = internal_edge
        
        do ie = 1, i_nedge
            e = e_internal_indexing_array(ie)
            ! this is a loop of all internal edges.
            
            i1 = e_p_obj_relation_array(e_p_index_array(e)-1)
            i2 = e_p_obj_relation_array(e_p_index_array(e))
            ! i1 and i2 are the constituent points of edge e
            
            call centroid_assembler(e, centroid_array)
            
            in_progress_projection(:) = 0.0
            in_progress_centroid = e_centroid(e,:)
            
            c1(:) = coords(i1,:)
            c2(:) = coords(i2,:)
            
            call centroid_array_routine(in_progress_projection,  in_progress_centroid, centroid_array)
            
            direction_array(:) = c1(:) - c2(:)
            
            angle = alignment(in_progress_projection, direction_array)
            
            if (angle .le. 0.) then
                sn(ie,:) = in_progress_projection(:)
			else
                sn(ie,:) = -in_progress_projection
			endif
			
! 			print*, ie, sn(ie,:), angle
            
        enddo
        
    end subroutine internal_edge_volume_processing
    
        
    subroutine boundary_edge_volume_processing
        implicit none
        integer(KIND=INT32) :: be, e, i, j, i1, i2, centroid_array_count, centroid_array_count_old, ef, ec, ef_start, ec_start, ef_end, ec_end
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
        
        centroid_array_count_old = -1
        allocate(centroid_index_array(0),non_viable_faces(0),centroid_array(0,0))
        
        obj_type = boundary_edge
        
        do be = 1, b_nedge
            e = e_bound_indexing_array(be)
            ! this is a loop of all boundary edges.
            
            i1 = e_p_obj_relation_array(e_p_index_array(e)-1)
            i2 = e_p_obj_relation_array(e_p_index_array(e))
            ! i1 and i2 are the constituent points of edge e
            
            call centroid_assembler(e, centroid_array)
            
            in_progress_projection(:) = 0.0
            in_progress_centroid = e_centroid(e,:)
            
            c1(:) = coords(i1,:)
            c2(:) = coords(i2,:)
            
            call centroid_array_routine(in_progress_projection,  in_progress_centroid, centroid_array)
            
            direction_array(:) = c1(:) - c2(:)
            
            angle = alignment(in_progress_projection, direction_array)
            
            if (angle .le. 0.) then
                sb(be,:) = in_progress_projection(:)
			else
                sb(be,:) = -in_progress_projection
			endif
			
! 			print*, be, sb(be,:), angle
            
        enddo
        
    end subroutine boundary_edge_volume_processing
    
    
    subroutine boundary_face_volume_processing
        implicit none
        integer(KIND=INT32) :: bp, p, i, j, i1, i2, centroid_array_count, centroid_array_count_real, centroid_array_count_old, pe, pe_start, pf_start, pe_end, pf_end
        integer(KIND=INT32) :: edge_count, face_count, current_edge, face_1, face_2, prev_face,be
        integer(KIND=INT32) :: k, pf, edge_count_real, face_count_real, remainder_count, current_face
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
            
            
            if (feature_points(bp)) then
                ! corners
                obj_type = featre_point
                
                i = bp ! we need bp inside of the feature edges and we cant modify bp in the loop
                call centroid_assembler(i, centroid_array)
                
            else
                ! not corners
                obj_type = non_feature_point
                
                call centroid_assembler(p, centroid_array)
                
            endif
            
            in_progress_projection(:) = 0.0
            in_progress_centroid = coords(p,:)
            
            call centroid_array_routine(in_progress_projection,  in_progress_centroid, centroid_array)
            
            direction_array(:) = p_normal_vectors(bp,:)
            
            angle = alignment(in_progress_projection, direction_array)
            
            if (angle .le. 0.) then
                sbb(bp,:) = in_progress_projection(:)
            else
                sbb(bp,:) = -in_progress_projection(:)
            endif
            
! 			print*, bp, sbb(bp,:), angle
            
        enddo
        
    end subroutine boundary_face_volume_processing

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine centroid_assembler(obj, centroid_array)
        implicit none
        integer(KIND=INT32) :: obj
        integer(KIND=INT32) :: mm, tt(2), i, j, k, m_i, m_stt, m_end, t_stt, t_end, fwd_i, bck_i
        integer(KIND=INT32) :: main_count, tert_count, centroid_array_count
        integer(KIND=INT32),allocatable :: main_obj(:), tert_obj(:), centroid_obj_array(:)
        real(KIND=REAL64)  ,allocatable :: centroid_array(:,:)
        logical, allocatable :: viable_mains(:)
        logical ::  e_state, state
        
        select case(obj_type)
            case(internal_edge)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
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
                
                ! allocate my temp arrays
                allocate(viable_mains(main_count), centroid_obj_array(centroid_array_count))
                centroid_obj_array = -1
                
                ! all main objects are viable to begin with.
                viable_mains = .true.
                
                ! inital edge setup doesnt matter much for i_edge
                ! just pick a random edge, my end condition is whenever the first cell is equal to the last cell
                m_i = 1
            
        
                call obj_select(m_i, m_stt, mm, tt)
                
                ! we just assign our values into the first available slots
                centroid_obj_array(1) = tt(1)
                centroid_obj_array(2) = mm
                centroid_obj_array(3) = tt(2)
                
                ! and initialise our indexers
                fwd_i = 3
                bck_i = 1
                
            case(boundary_edge)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                ! faces are main
                ! cells are tertiary
                m_stt   = e_f_index_array(obj-1)
                t_stt   = e_c_index_array(obj-1)
                m_end   = e_f_index_array(obj)
                t_end   = e_c_index_array(obj)
                
                main_count = m_end-m_stt
                tert_count = t_end-t_stt
                
                ! for b_edge the centroid will consist of all adjacent centroids
                centroid_array_count = 2 + main_count + tert_count
                ! but for convinience we accept that the last entry can be dropped so subroutines can be simpler
                
                ! allocate my temp arrays
                allocate(viable_mains(main_count), centroid_obj_array(centroid_array_count))
                centroid_obj_array = -1
                
                ! all main objects are viable to begin with.
                viable_mains = .true.
                
                ! we must start on a boundary face
                m_i = 1
                do 
                    mm = e_f_obj_relation_array(m_stt+m_i)
                    if (f_bound_array(mm)) exit
                    m_i = m_i + 1
                enddo
                
                call obj_select(m_i, m_stt, mm, tt)
                
!                 if ((obj.eq. 4).or.(obj.eq.10)) then
!                     print*, obj, mm, tt(1), tt(2)
!                     print*,e_f_obj_relation_array(m_stt:m_end)
!                     print*,f_bound_array(e_f_obj_relation_array(m_stt:m_end))
!                     print*,e_c_obj_relation_array(t_stt:t_end)
!                 endif
                
                ! this differs from internal because our first val will be a face and only connects to one cell
                centroid_obj_array(1) = -1
                centroid_obj_array(2) = mm
                centroid_obj_array(3) = tt(1)
                ! tt(1) and tt(2) should be identical here 
                
                ! and initialise our indexers
                fwd_i = 3
                bck_i = 1
                
                ! other than the current edge
                viable_mains(m_i) = .false.
                
            case(non_feature_point)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
                ! edges are main
                ! faces are tertiary
                m_stt   = p_e_index_array(obj-1)
                t_stt   = p_f_index_array(obj-1)
                m_end   = p_e_index_array(obj)
                t_end   = p_f_index_array(obj)
                
                main_count = m_end-m_stt
                tert_count = t_end-t_stt
                
                ! allocate my temp arrays
                allocate(viable_mains(main_count))
                
                ! here we flag only boundary edges as viable
                viable_mains = .false.
                
                centroid_array_count = 1 + main_count + tert_count
                
                ! we also remove the internal centroids from the centroid array count
                m_i = 1
                do m_i=1, main_count
                    mm = p_e_obj_relation_array(m_stt+m_i)
                    if (e_bound_array(mm)) then
                        viable_mains(m_i) = .true.
                    else
                        centroid_array_count = centroid_array_count - 1
                    endif
                enddo
                
                ! and the same for internal faces
                m_i = 1
                do m_i=1, tert_count
                    mm = p_f_obj_relation_array(t_stt+m_i)
                    if (.not. f_bound_array(mm)) centroid_array_count = centroid_array_count - 1
                enddo
                
                allocate(centroid_obj_array(centroid_array_count))
                centroid_obj_array = -1
                
                m_i = 1
                ! then select a viable edge (ie a boundary edge)
                do 
                    if (viable_mains(m_i)) exit
                    m_i = m_i + 1
                enddo
                
                call obj_select(m_i, m_stt, mm, tt)
                
                ! we just assign our values into the first available slots
                centroid_obj_array(1) = tt(1)
                centroid_obj_array(2) = mm
                centroid_obj_array(3) = tt(2)
                
                ! and initialise our indexers
                fwd_i = 3
                bck_i = 1
                
            case(featre_point)!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                
                ! we need bp only for the flag
                flag = p_boundary_flags(obj)
                
                ! and p for everything else
                obj = p_bound_indexing_array(obj)
                
                ! edges are main
                ! faces are tertiary
                m_stt   = p_e_index_array(obj-1)
                t_stt   = p_f_index_array(obj-1)
                m_end   = p_e_index_array(obj)
                t_end   = p_f_index_array(obj)
                
                main_count = m_end-m_stt
                tert_count = t_end-t_stt
                
                ! allocate my temp arrays
                allocate(viable_mains(main_count))
                
                ! here we flag only boundary edges as viable
                ! unlike with non feature points however, we must also only flag feature edges and edges with the correct boundary flag.
                viable_mains = .false.
                
                !centroid_array_count = 2 + main_count + tert_count
                centroid_array_count = 2
                ! we also remove the wrong flagged centroids from the centroid array count
                m_i = 1
                do m_i=1, main_count
                    mm = p_e_obj_relation_array(m_stt+m_i)
                    if (e_bound_array(mm)) then
!                         print*, m_i, mm, e_boundary_flags(reversed_e_bound_indexing_array(mm)), viable_mains(m_i)
                        
                        if     (e_boundary_flags(reversed_e_bound_indexing_array(mm)) .eq. flag) then
                            viable_mains(m_i) = .true.
                            centroid_array_count = centroid_array_count + 1
                        elseif (e_boundary_flags(reversed_e_bound_indexing_array(mm)) .eq. 1000) then
                            ! now we only want the feature edges that are adjacent to at least one correctly flagged face
                            
                            ! so we pull our boundary faces
                            i = e_f_index_array(mm-1)
                            do 
                                i=i+1
                                if (f_bound_array(e_f_obj_relation_array(i))) exit
                            enddo
                            tt(1) = e_f_obj_relation_array(i)
                            do 
                                i=i+1
                                if (f_bound_array(e_f_obj_relation_array(i))) exit
                            enddo
                            tt(2) = e_f_obj_relation_array(i)
                            
                            ! and then we check if either of them are flagged correctly
                            
                            if ((f_boundary_flags(reversed_f_bound_indexing_array(tt(1))) .eq. flag ) .or. (f_boundary_flags(reversed_f_bound_indexing_array(tt(2))) .eq. flag )) then
                                viable_mains(m_i) = .true.
                                centroid_array_count = centroid_array_count + 1
                            endif
                            
                        else
                            !centroid_array_count = centroid_array_count - 1
                        endif
                        
!                         print*, m_i, mm, e_boundary_flags(reversed_e_bound_indexing_array(mm)), viable_mains(m_i)
!                         
!                         print*, ' ' 
                    else
                        !centroid_array_count = centroid_array_count - 1
                    endif
                    !print*, m_i, mm, e_boundary_flags(reversed_e_bound_indexing_array(mm))
                enddo
                
                ! and the same for internal faces but we dont need to worry about features
                m_i = 1
                do m_i=1, tert_count
                    mm = p_f_obj_relation_array(t_stt+m_i)
                    if (f_bound_array(mm)) then
                        ! no wrong flags
                        if (f_boundary_flags(reversed_f_bound_indexing_array(mm)) .eq. flag) then
                            centroid_array_count = centroid_array_count + 1
                        endif
                    ! no internals
                    else
                        
                    endif
                enddo
                
                allocate(centroid_obj_array(centroid_array_count))
                centroid_obj_array = -1
                
                ! starting arrangements here are more difficult, 
                ! we need a feature edge and then we need to make sure we select only the correctly flagged face
                
                ! so we make sure we start on a feature edge, preferably one that is bounding the region
                
                do m_i=1, main_count
                    if (viable_mains(m_i)) then
                        call obj_select(m_i, m_stt, mm, tt)
                        if (e_boundary_flags(reversed_e_bound_indexing_array(mm)) .eq. 1000) then
                            if (tt(1) .eq. tt(2)) exit
                        endif
                    endif
                enddo
                viable_mains(m_i) = .false.
                
                ! we just assign our values into the first available slots
                centroid_obj_array(1) = -1
                centroid_obj_array(2) = mm
                centroid_obj_array(3) = tt(2)
                
                ! and initialise our indexers
                fwd_i = 3
                bck_i = 1
                
                
            
        end select!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
        
        
!         if (obj_type.eq.featre_point) then
!             print*, 'new point', obj, flag, coords(obj,:)
!             print*, tt(1), mm, tt(2)
!             print*, centroid_obj_array
!         endif
        if (obj_type.eq.boundary_edge) then
            if ((obj.eq. 4).or.(obj.eq.10)) then
                print*, flag, mm, tt(1), tt(2), f_bound_array(mm)
                print*, centroid_obj_array
                print*, ' '
                print*, ' '
            endif
        endif
        
        
        do k=1,1000 ! this k limit is here to stop infinite loops, it needs to be set to an unreasonable number to not trip incorrectly
            
            ! now we loop about our centroids
            if (m_i.eq.main_count) then
                m_i = 1
            else
                m_i=m_i+1
            endif
            
            ! if we know the centroid isnt viable we skip it
            if (.not.viable_mains(m_i)) cycle
            
            ! obviously we get our objects
            call obj_select(m_i, m_stt, mm, tt)
            
            call centroid_swapper(mm, tt, fwd_i, bck_i, centroid_obj_array, state)
            
            if (obj_type.eq.boundary_edge) then
                if ((obj.eq. 4).or.(obj.eq.10)) then
                    print*, flag, mm, tt(1), tt(2), f_bound_array(mm)
                    print*, centroid_obj_array
                    print*, ' '
                    print*, ' '
                endif
            endif
            
            if (state) then
                
                viable_mains(m_i) = .false.
                
                call obj_endconditions(e_state, fwd_i, bck_i, centroid_obj_array, centroid_array_count)
                ! our end condition is the first edge and last edge being the same
                if (e_state) exit
            endif
            
            
        enddo
        
        if (k.eq.1001) then
        
!             if (obj_type.eq.boundary_edge)then
!                 print*, 'new edge', obj
!                 print*, centroid_obj_array
!             endif    
            
            print*, 'weve hit 1000 loops in the centroid assembler so the end conditions are probably munted'
            
            
            return
        endif
        
!         if (obj_type.eq.featre_point) then
!             print*, centroid_obj_array
!         endif
        
        call centroid_float_assembler(centroid_obj_array, fwd_i, bck_i, centroid_array, centroid_array_count)
        
        
    end subroutine centroid_assembler
    
    subroutine obj_select(m_i, m_stt, mm, tt)
        implicit none
        integer(KIND=INT32) :: i, m_i, m_stt, mm, tt(2)
        
        
        select case(obj_type)
            case(internal_edge)
            
                ! faces are main
                ! cells are tertiary
                mm    = e_f_obj_relation_array(m_stt + m_i)
                tt(1) = f_c_obj_relation_array(f_c_index_array(mm)    )
                tt(2) = f_c_obj_relation_array(f_c_index_array(mm-1)+1)
                
            case(boundary_edge)
                
                ! faces are main
                ! cells are tertiary
                mm    = e_f_obj_relation_array(m_stt + m_i)
                tt(1) = f_c_obj_relation_array(f_c_index_array(mm)    )
                tt(2) = f_c_obj_relation_array(f_c_index_array(mm-1)+1)
                
            case(non_feature_point)
                
                ! edges are main
                ! faces are tertiary
                mm = p_e_obj_relation_array(m_stt+m_i)
                
                i = e_f_index_array(mm-1)
                do 
                    i=i+1
                    if (f_bound_array(e_f_obj_relation_array(i))) exit
                enddo
                tt(1) = e_f_obj_relation_array(i)
                do 
                    i=i+1
                    if (f_bound_array(e_f_obj_relation_array(i))) exit
                enddo
                tt(2) = e_f_obj_relation_array(i)
                
                if (i.gt.e_f_index_array(mm))then
                    print*, ' non-feature boundary edge has more than 2 faces '
                    stop
                endif
                
                
            case(featre_point)
                
                mm = p_e_obj_relation_array(m_stt+m_i)
                
                if (e_boundary_flags(reversed_e_bound_indexing_array(mm)) .eq. 1000) then
                    
                    ! our feature edge will only have one correctly flagged face
                    i = e_f_index_array(mm-1)
                    do 
                        i=i+1
                        if (f_bound_array(e_f_obj_relation_array(i))) exit
                    enddo
                    tt(1) = e_f_obj_relation_array(i)
                    do 
                        i=i+1
                        if (f_bound_array(e_f_obj_relation_array(i))) exit
                    enddo
                    tt(2) = e_f_obj_relation_array(i)
                    
!                     print*, tt(1), f_boundary_flags(reversed_f_bound_indexing_array(tt(1))), tt(2), f_boundary_flags(reversed_f_bound_indexing_array(tt(2)))
                    
                    if (f_boundary_flags(reversed_f_bound_indexing_array(tt(1))) .ne. flag ) then
                        tt(1) = tt(2)
                    endif
                    if (f_boundary_flags(reversed_f_bound_indexing_array(tt(2))) .ne. flag ) then
                        tt(2) = tt(1)
                    endif
                    
                else
                    ! else proceed as a normal edge would in the non feature point list
                    ! because the only normal edges left should be edges where both faces are flagged right
                    i = e_f_index_array(mm-1)
                    do 
                        i=i+1
                        if (f_bound_array(e_f_obj_relation_array(i))) exit
                    enddo
                    tt(1) = e_f_obj_relation_array(i)
                    do 
                        i=i+1
                        if (f_bound_array(e_f_obj_relation_array(i))) exit
                    enddo
                    tt(2) = e_f_obj_relation_array(i)
                    
                    if (i.gt.e_f_index_array(mm))then
                        print*, ' non-feature boundary edge has more than 2 faces '
                        stop
                    endif
                    
                endif
                
                
        end select
            
    end subroutine obj_select
    
    subroutine centroid_float_assembler(centroid_obj_array, fwd_i, bck_i, centroid_array, centroid_array_count)
        implicit none
        integer(KIND=INT32) :: i, centroid_obj_array(:), centroid_array_count, fwd_i, bck_i
        real(KIND=REAL64)  ,allocatable :: centroid_array(:,:)
        
        deallocate(centroid_array)
        
        select case(obj_type)
            case(internal_edge)
                
                allocate(centroid_array(centroid_array_count,3))
                
                ! start on a tertiary
                
                ! faces are main
                ! cells are tertiary
                do i=1,centroid_array_count,2
                    centroid_array(i,:) = c_centroid(centroid_obj_array(i),:)
                enddo
                do i=2,centroid_array_count-1,2
                    centroid_array(i,:) = f_centroid(centroid_obj_array(i),:)
                enddo
                
            case(boundary_edge)
                
                ! we allowed ourselves one extra space on our allocations to make the swapper happy, 
                ! so we drop that
                allocate(centroid_array(centroid_array_count-2,3))
                
                ! start on a main
                
!                 faces are main
!                 cells are tertiary
                do i=2,centroid_array_count-1,2
                    centroid_array(i,:) = f_centroid(centroid_obj_array(i),:)
                enddo
                do i=3,centroid_array_count-2,2
                    centroid_array(i,:) = c_centroid(centroid_obj_array(i),:)
                enddo
                
                
            case(non_feature_point)
                
                allocate(centroid_array(centroid_array_count,3))
                
                ! start on a tertiary
                
!                 edges are main
!                 faces are tertiary
                do i=1,centroid_array_count,2
                    centroid_array(i,:) = f_centroid(centroid_obj_array(i),:)
                enddo
                do i=2,centroid_array_count-1,2
                    centroid_array(i,:) = e_centroid(centroid_obj_array(i),:)
                enddo
                
                
            case(featre_point)
                
                ! again we allowed ourselves one extra space on our allocations to make the swapper happy, 
                ! so we drop that
                ! and also a space at the beginning that we drop as well
                allocate(centroid_array(centroid_array_count-2,3))
            
                ! start on a main
                
!                 edges are main
!                 faces are tertiary
                do i=2,centroid_array_count-1,2
                    centroid_array(i-1,:) = e_centroid(centroid_obj_array(i),:)
                enddo
                do i=3,centroid_array_count-2,2
                    centroid_array(i-1,:) = f_centroid(centroid_obj_array(i),:)
                enddo
                
        end select
        
    end subroutine centroid_float_assembler
    
    subroutine obj_endconditions(e_state, fwd_i, bck_i, centroid_obj_array, centroid_array_count)
        implicit none
        integer(KIND=INT32) :: fwd_i, bck_i, centroid_array_count, i
        integer(KIND=INT32) :: centroid_obj_array(:)
        logical ::  e_state
        
        e_state = .false.
        select case(obj_type)
            
            case(internal_edge)
            
                if (centroid_obj_array(1) .eq. centroid_obj_array(centroid_array_count)) e_state = .true.
                
            case(boundary_edge)
                
                ! these two conditions are effectivly identical in function because if 
                ! the last and 2nd to last values are the same it implies a face is only connected to one cell
                ! and is therefore a boundary face
                
                if (centroid_obj_array(centroid_array_count-1) .ne. -1) then
                    if (f_bound_array(centroid_obj_array(centroid_array_count-1))) e_state = .true.
                endif
                
                ! this causes problems however because it doesnt check if the equal centroids are -1
                
!                 if (centroid_obj_array(centroid_array_count) .eq. centroid_obj_array(centroid_array_count-2)) e_state = .true.
                
            case(non_feature_point)
                
                if (centroid_obj_array(1) .eq. centroid_obj_array(centroid_array_count)) e_state = .true.
                
!                 if (centroid_obj_array(fwd_i) .eq. centroid_obj_array(bck_i)) e_state = .true.
                
            case(featre_point)
                
!                 print*, centroid_array_count
!                 print*, centroid_obj_array
                
                if ((centroid_obj_array(2) .ne. -1) .and. (centroid_obj_array(centroid_array_count-1) .ne. -1)) then
                    if (feature_edges(reversed_e_bound_indexing_array(centroid_obj_array(2))) .and. feature_edges(reversed_e_bound_indexing_array(centroid_obj_array(centroid_array_count-1)))) e_state = .true.
                    !if (centroid_obj_array(centroid_array_count) .eq. centroid_obj_array(centroid_array_count-2)) e_state = .true.
                endif
                
                
                if ((centroid_obj_array(3) .ne. -1) .and. (centroid_obj_array(centroid_array_count-1) .ne. -1)) then
                    if (centroid_obj_array(3) .eq. centroid_obj_array(centroid_array_count-1)) e_state = .true.
                    !if (centroid_obj_array(centroid_array_count) .eq. centroid_obj_array(centroid_array_count-2)) e_state = .true.
                endif
                
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
    
    subroutine centroid_array_routine(obj_projection, obj_centroid, centroid_array)
        implicit none
        integer(KIND=INT32)        :: i, centroid_array_count
        real(KIND=REAL64)          :: obj_projection(3), centroid_array(:,:), obj_centroid(3), baryobj_1_centroid(3), baryobj_2_centroid(3)
        
        centroid_array_count = size(centroid_array,1)

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
