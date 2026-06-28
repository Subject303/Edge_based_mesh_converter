


module read_data_serial
    
    use iso_fortran_env
    
    use object_counts
    use raw_data
    use io_data
    
    implicit none
    
    contains
    
    subroutine read_data
        implicit none
        
        print*, ' '
        print *, 'Reading_raw_data'
        
        open(10,file=raw_data_path,access='stream',action='read',status='old')
        
        read(10) npoin, nedge, nface, nele
        print *, 'mesh contains ', npoin, ' points, ', nele, ' elements, ', nface, ' faces, ', nedge, ' edges '
        print*, ' '
        
        allocate( coords(npoin,3), c_p_index_array(0:nele), f_p_index_array(0:nface), e_p_index_array(0:nedge) )
        
        Print*, 'reading coordinates '

        read(10) coords(:,1)
        read(10) coords(:,2)
        read(10) coords(:,3)
        
        Print*, 'reading element connectivity index array '
        
        read(10) c_p_sum
        read(10) c_p_index_array(:)
        
        Print*, 'reading face connectivity index array '
        read(10) f_p_sum
        read(10) f_p_index_array(:)
        
        allocate ( c_p_obj_relation_array(c_p_sum), f_p_obj_relation_array(f_p_sum), e_p_obj_relation_array(2*nedge) )
        
        Print*, 'reading element connectivity array '
        read(10) c_p_obj_relation_array(:)
        
        Print*, 'reading face connectivity array '
        read(10) f_p_obj_relation_array(:)
        
        Print*, 'reading edge connectivity array '
        read(10) e_p_obj_relation_array(:)
        
        allocate ( c_f_index_array(0:nele), f_e_index_array(0:nface)  )
        e_p_sum = nedge*2
        
        print*, 'reading cell to face index array '
        read(10) c_f_sum
        read(10) c_f_index_array(:)
        
        print*, 'reading face to edge index array '
        read(10) f_e_sum
        read(10) f_e_index_array(:)
        ! note: an assumption here is that I've broken feature edges but not points into two for boundaries in the preprocessor
        ! this is so I can alternate boundary faces then edges to grow my boundary regions,
        ! but still have it reference the same point in all other point relation arrays
        !
        ! IE.: feature edges that are split will have two different edge IDs despite having the same points inside, e1 = (p1,p2) & e2 = (p1,p2)
        !
        ! might have to filter out stuff for kbface?
        
        allocate ( c_f_obj_relation_array(c_f_sum), f_e_obj_relation_array(f_e_sum) )
        
        print*, 'reading cell to face relation array '
        read(10) c_f_obj_relation_array(:)
        
        print*, 'reading face to edge relation array '
        read(10) f_e_obj_relation_array(:)
        
        c_f_obj_relation_array = c_f_obj_relation_array - 1
        f_e_obj_relation_array = f_e_obj_relation_array - 1
        
        print*,c_f_sum, SIZE(c_f_obj_relation_array), nele, SIZE(c_f_index_array)
        print*, 'aa'
        print*,c_f_index_array
        print*, 'aa'
        print*,c_f_obj_relation_array
        print*, 'aa'        
        print*,f_e_sum, SIZE(f_e_obj_relation_array), nface, SIZE(f_e_index_array)
        print*, 'aa'
        print*,f_e_index_array
        print*, 'aa'
        print*,f_e_obj_relation_array
        print*, 'aa'
        
        c_f_obj_relation_array = c_f_obj_relation_array + 1
        f_e_obj_relation_array = f_e_obj_relation_array + 1
        
        close(10)
        
        print*, 'Relation array sums: '
        print*, 'c_p_sum,  f_p_sum,  e_p_sum,  c_f_sum,  f_e_sum'
        print*,  c_p_sum,f_p_sum,e_p_sum,c_f_sum,f_e_sum
        
        print*, 'finished reading raw data'
        
        print*, 'total memory size of raw data :: ', ( (8*3*npoin) + (4*( 4 + nele + nface + c_p_sum + f_p_sum + (2*nedge) + 2 + nele + nface + c_f_sum + f_e_sum + 4)))  , ' bytes'
        print*, ' '
        
    end subroutine read_data
    
end module read_data_serial
