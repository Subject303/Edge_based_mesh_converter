    
module data_outputting_serial
    use io_data
    use object_counts
    use raw_data
    use object_relation_data
    use boundary_data
    use projection_data
    
    contains
    
    subroutine data_outputting
        
        
        implicit none
        
        
        
        if ((output_requests - binary_internal) .ge. 0) then
			call output_binary_internal
        endif
        
        if ((output_requests - ascii_internal) .ge. 0) then
			call output_ascii_internal
        endif
        
        if ((output_requests - vtu_ascii) .ge. 0) then
			call output_vtu_ascii
        endif
        
        if ((output_requests - vtu_binary_appended) .ge. 0) then
			call output_vtu_binary_appended
        endif
        
        
    end subroutine data_outputting
    
    subroutine output_binary_internal
        implicit none
        
    end subroutine output_binary_internal
    
    subroutine output_ascii_internal
        implicit none
        
    end subroutine output_ascii_internal
    
    subroutine output_vtu_ascii
        implicit none
        
    end subroutine output_vtu_ascii
    
    subroutine output_vtu_binary_appended
        implicit none
		integer(KIND=INT32)             :: i, c, f, e, offset_count, cf, fe, p, num_e, written_p
		integer(KIND=INT8 ),allocatable :: types(:)
        integer(KIND=INT64)			    :: offset, face_offset
        character(:), allocatable		:: outstring
        character(16)                   :: offset_text, TXTnele, TXTnpoin
        logical, allocatable		    :: viable_edge(:)
        
        print*, NEW_LINE(''),NEW_LINE(''),NEW_LINE(''),NEW_LINE(''),NEW_LINE('')
        print*, 'WRITING VTU BINARY FILE'
        print*, NEW_LINE(''),NEW_LINE(''),NEW_LINE(''),NEW_LINE(''),NEW_LINE('')
        
        write(TXTnpoin, '(I16)') npoin
        write(TXTnele,  '(I16)') nele
        
		face_offset = 0
		do c=1,nele
            face_offset = face_offset + 1! + c_f_index_array(c) - c_f_index_array(c-1)
            do cf=(1+c_f_index_array(c-1)),c_f_index_array(c)
                f = c_f_obj_relation_array(cf)
                face_offset = face_offset + 1 + (f_p_index_array(f)-f_p_index_array(f-1))
            enddo
        enddo
        face_offset = face_offset * 4
        
            outstring = '&
&<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" header_type="UInt64">'//NEW_LINE('')//'&
&<UnstructuredGrid>'//NEW_LINE('')//'&
&<Piece NumberOfPoints="'//TXTnpoin//'" NumberOfCells="'//TXTnele//'">'//NEW_LINE('')//'&
&<Points>'//NEW_LINE('')
            offset = 0 ;write(offset_text, '(i16)')  offset
    outstring = outstring//'<DataArray type="Float32" NumberOfComponents="3"  Name="Points" format="appended" offset="'//offset_text//'" ></DataArray>'//NEW_LINE('')//'&
&</Points>'//NEW_LINE('')//'&
&<Cells>'//NEW_LINE('')
        offset = offset + 8 + npoin*3*4
        write(offset_text, '(i16)')  offset
    outstring = outstring//'<DataArray type="Int32" Name="offsets" format="appended" offset="'//offset_text//'" ></DataArray>'//NEW_LINE('')
        offset = offset + 8 + nele*4
        write(offset_text, '(i16)')  offset
    outstring = outstring//'<DataArray type="UInt8" Name="types" format="appended" offset="'//offset_text//'" ></DataArray>'//NEW_LINE('')
        offset = offset + 8 + nele 
        write(offset_text, '(i16)')  offset
    outstring = outstring//'<DataArray type="Int32" Name="connectivity" format="appended" offset="'//offset_text//'" ></DataArray>'//NEW_LINE('')
        offset = offset + 8 + c_p_sum*4
        write(offset_text, '(i16)')  offset
    outstring = outstring//'<DataArray type="Int32" Name="faces" format="appended" offset="'//offset_text//'" ></DataArray>'//NEW_LINE('')
        offset = offset + 8 + face_offset
        write(offset_text, '(i16)')  offset
    outstring = outstring//'<DataArray type="Int32" Name="faceoffsets" format="appended" offset="'//offset_text//'" ></DataArray>'//NEW_LINE('')//'&
&</Cells>'//NEW_LINE('')//'&
&<PointData>'//NEW_LINE('')
        outstring = outstring//'</PointData>'//NEW_LINE('')//'&
&</Piece>'//NEW_LINE('')//'&
&</UnstructuredGrid>'//NEW_LINE('')//'&
&<AppendedData encoding="raw">'//&
&'_'
        
        open(1,file='vtu_mesh.vtu',access='stream',status='replace',convert='little_endian')
        
        write(1) outstring
        
        offset = 8+npoin*3*4
        write(1) offset
        do i=1,npoin
            write(1) coords(i,:)
		enddo
        
        offset = 8+nele*4
        write(1) offset, c_p_index_array(1:nele)
    
        allocate(types(nele))
        types = 42
        offset = 8+nele
        write(1) offset, types(:)
        
        offset = 8+c_p_sum*4
        write(1) offset, c_p_obj_relation_array(:)-1
    
        write(1) 8+face_offset
        do c=1,nele
            ! number of faces
            write(1) (c_f_index_array(c)-c_f_index_array(c-1)) 
            do cf=(1+c_f_index_array(c-1)),c_f_index_array(c)
                f = c_f_obj_relation_array(cf)
                ! number of points in face, points in face
                write(1) (f_p_index_array(f)-f_p_index_array(f-1)) , f_p_obj_relation_array((f_p_index_array(f-1)+1):f_p_index_array(f))-1
                
                enddo
            enddo
        enddo
    
        offset = 8+nele*4
        write(1) offset
        offset_count = 0
        do c=1,nele
            offset_count =  offset_count + 1 + c_f_index_array(c) - c_f_index_array(c-1)
            do cf=(1+c_f_index_array(c-1)),c_f_index_array(c)
                f = c_f_obj_relation_array(cf)
                offset_count = offset_count + (f_p_index_array(f)-f_p_index_array(f-1))
            enddo
            write(1) offset_count
        enddo
    
        outstring = NEW_LINE('')//'&
&</AppendedData>'//NEW_LINE('')//'&
&</VTKFile>'
    
        write(1) outstring
		
        close(1)
        
        print*, NEW_LINE(''),NEW_LINE(''),NEW_LINE(''),NEW_LINE(''),NEW_LINE('')
        print*, 'FINISHED WRITING VTU BINARY FILE'
        print*, NEW_LINE(''),NEW_LINE(''),NEW_LINE(''),NEW_LINE(''),NEW_LINE('')
        
    end subroutine output_vtu_binary_appended
    
end module data_outputting_serial
