

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
        integer(KIND=INT32) :: cell_count, face_count, current_face, current_cell, cell_1, cell_2, prev_cell
        integer(KIND=INT32),allocatable :: centroid_index_array(:)
        real(KIND=REAL32),allocatable   :: centroid_array(:)
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
        allocate(centroid_index_array(0),non_viable_faces(0))
        
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
                deallocate(centroid_index_array,non_viable_faces)
                allocate(centroid_index_array(centroid_array_count), non_viable_faces(face_count+1))
            endif
            
            current_face = e_f_obj_relation_array(ef_start+1)
            cell_1 = f_c_obj_relation_array(f_c_index_array(current_face)    )
            cell_2 = f_c_obj_relation_array(f_c_index_array(current_face-1)+1)
            
            centroid_index_array(1) = cell_1
            centroid_index_array(2) = current_face
            centroid_index_array(3) = cell_2
            prev_cell = cell_1
            
            print*, 'number of faces ', face_count, ' number of cells ', cell_count, ' count ', centroid_array_count
            
            i=2
            non_viable_faces = .false.
            
            ef = 1
            
            do
                
                if (ef.eq.face_count) then
                    ef = 1
                else
                    ef=ef+1
                endif
                
                print*, current_face, cell_1, cell_2 , prev_face

                
                if (non_viable_faces(ef)) cycle
                
                current_face = e_f_obj_relation_array(ef_start + ef)
                cell_1 = f_c_obj_relation_array(f_c_index_array(current_face)    )
                cell_2 = f_c_obj_relation_array(f_c_index_array(current_face-1)+1)
                
                if (prev_cell .eq. cell_1) then
            
                    centroid_index_array(i)   = current_face
                    centroid_index_array(i+1) = cell_2
                    non_viable_faces(ef) = .true.
                    
                    prev_cell = centroid_index_array(i-1)
                    i=i+2 
                    
                    if (i.eq.centroid_array_count+1) exit
                    
                elseif (prev_cell .eq. cell_2) then
                    
                    centroid_index_array(i)   = current_face
                    centroid_index_array(i+1) = cell_1
                    non_viable_faces(ef) = .true.
                    
                    prev_cell = centroid_index_array(i-1)
                    i=i+2 
                    
                    if (i.eq.centroid_array_count+1) exit
                    
                else
                endif
                
            enddo
            
            centroid_index_array(centroid_array_count) = centroid_index_array(1)
            
            
            centroid_array_count_old = centroid_array_count
            
            
            if (i.ne.centroid_array_count+1) print*, 'centroid array count broken'
            if (centroid_index_array(1).ne.centroid_index_array(centroid_array_count)) print*, 'internal centroid array start and end wrong'
            
            print*, centroid_index_array
            
            
            
            !call centroid_array_routine(sn(e,:), e_centroid(e,:), centroid_array_count, centroid_array, i1, i2)
            
        enddo
        
        
        
    end subroutine internal_edge_volume_processing
    
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine centroid_array_routine(obj_projection, obj_centroid, centroid_array_count, centroid_array, i1, i2)
        implicit none
        integer(KIND=INT32)          :: i, centroid_array_count, i1
        integer(KIND=INT32),optional :: i2
        real(KIND=REAL32) :: obj_projection(3), centroid_array(:,:), obj_centroid(3)
        
        do i=2,centroid_array_count
            call cvolume(obj_projection, vol(i1), vol(i2), i1, i2, obj_centroid, centroid_array(i-1,:), centroid_array(i,:))
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
        sn(2) = sn(1) + ((v34(3) * v35(1) - v34(1) * v35(3))/2) ! projection in the yy axis
        sn(3) = sn(1) + ((v34(1) * v35(2) - v34(2) * v35(1))/2) ! projection in the zz axis
        
    !                               i4--------------i5
    !                                \  >>>>>>>   /
    !                                 \         /
    !                                  \      / 
    !                                   \   /
    !           i1 -------------------- i3 -------------------------i2
    end subroutine cvolume
    
    ! end contains
    
end module volume_processing
