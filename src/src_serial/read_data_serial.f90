


module read_data_serial
    
    use iso_fortran_env
    
    use object_counts
    use raw_data
    use io_data
    
    implicit none
    
    contains
    
    subroutine read_data
        
        print *, 'Reading_raw_data'
        
        open(10,file=raw_data_path,access='stream',action='read',status='old')
        
        read(10) npoin, nele, nface, nedge
        print *, 'mesh contains ', npoin, ' points, ', nele, ' elements, ', nface, ' faces, ', nedge, ' edges '
        
        allocate( coords(npoin,3), element_con_index(nele), face_con_index(nface) )
        
        Print*, 'reading coordinates '

        read(10) coords(:,1)
        read(10) coords(:,2)
        read(10) coords(:,3)
        
        Print*, 'reading element connectivity index array '
        
        read(10) element_con_sum
        read(10) element_con_index(:)
        
        Print*, 'reading face connectivity index array '
        
        read(10) face_con_sum
        read(10) face_con_index(:)
        
        allocate ( element_connectivity(element_con_sum), face_connectivity(face_con_sum), edge_connectivity(2*nedge) )
        
        Print*, 'reading element connectivity array '
        read(10) element_connectivity(:)
        
        Print*, 'reading face connectivity array '
        read(10) face_connectivity(:)
        
        Print*, 'reading edge connectivity array '
        read(10) edge_connectivity(:)
        
    end subroutine read_data
    
end module read_data_serial
