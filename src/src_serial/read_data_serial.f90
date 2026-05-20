


module read_data_serial
    
    use iso_fortran_env
    
    use object_counts
    use raw_data
    use io_data
    
    implicit none
    
    contains
    
    subroutine read_data
        
        print*, ' '
        print *, 'Reading_raw_data'
        
        open(10,file=raw_data_path,access='stream',action='read',status='old')
        
        read(10) npoin, nele, nface, nedge
        print *, 'mesh contains ', npoin, ' points, ', nele, ' elements, ', nface, ' faces, ', nedge, ' edges '
        print*, ' '
        
        allocate( coords(npoin,3), c_p_index_array(nele), f_p_index_array(nface) )
        
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
        
        close(10)
        
        print*, 'total memory size of raw data :: ', ( (8*3*npoin) + (4*(4+nele+nface+c_p_sum+f_p_sum+(2*nedge)))) ) , ' bytes'
        print*, ' '
        
    end subroutine read_data
    
end module read_data_serial
