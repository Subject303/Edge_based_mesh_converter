    
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
		integer(KIND=INT32)             :: c, f, offset_count, cf, fp1, fp2
		integer(KIND=INT32 )            :: elesize(nele)
        integer(KIND=INT64)			    :: offset, face_offset
        
        
        open(2,file='internal_mesh.internalmesh',access='stream',status='replace',convert='little_endian')
        
        write(2) npoin, i_nedge, b_nedge, b_npoin, nele
        
        ! double check this should be i1 i1 i1 i2 i2 i2 or i1 i2 i1 i2 i1 i2
        write(2) e_p_obj_relation_array(e_p_index_array(e_internal_indexing_array(:))  )
        write(2) e_p_obj_relation_array(e_p_index_array(e_internal_indexing_array(:))-1)
        
!         print*, e_p_obj_relation_array(e_p_index_array(e_internal_indexing_array(:))  ), e_p_obj_relation_array(e_p_index_array(e_internal_indexing_array(:))-1)
        
        write(2) e_p_obj_relation_array(e_p_index_array(e_bound_indexing_array(:))  )
        write(2) e_p_obj_relation_array(e_p_index_array(e_bound_indexing_array(:))-1)
        
!         print*, e_p_obj_relation_array(e_p_index_array(e_bound_indexing_array(:))  ), e_p_obj_relation_array(e_p_index_array(e_bound_indexing_array(:))-1)
        
        write(2) p_bound_indexing_array
        write(2) p_boundary_flags
        
!         print*, p_bound_indexing_array
        
        write(2) coords(:,1)
        write(2) coords(:,2)
        write(2) coords(:,3)
        write(2) vol
        write(2) sn(:,1)
        write(2) sn(:,2)
        write(2) sn(:,3)
        write(2) sb(:,1)
        write(2) sb(:,2)
        write(2) sb(:,3)
        write(2) sbb(:,1)
        write(2) sbb(:,2)
        write(2) sbb(:,3)
        write(2) p_normal_vectors(:,1)
        write(2) p_normal_vectors(:,2)
        write(2) p_normal_vectors(:,3)
        
        do c=1,nele
            
            elesize(c) = c_p_index_array(c) - c_p_index_array(c-1)
            
        enddo
        
        print*, elesize
        
!         ! offsets
!         write(2) c_p_index_array(1:nele)
! !         ! types
! !         write(2) types(:)
!         ! connectivity
!         write(2) c_p_obj_relation_array(:)-1
!         ! faces
!         do c=1,nele
!             ! number of faces
!             write(2) (c_f_index_array(c)-c_f_index_array(c-1)) 
!             do cf=(1+c_f_index_array(c-1)),c_f_index_array(c)
!                 f = c_f_obj_relation_array(cf)
!                 fp1 = f_p_index_array(f-1)
!                 fp2 = f_p_index_array(f)
!                 ! number of points in face, points in face
!                 write(2) (fp2-fp1) , f_p_obj_relation_array((fp1+1):fp2)-1
!                 
!             enddo
!         enddo
!         ! faceoffsets
!         offset_count = 0
!         do c=1,nele
!             offset_count =  offset_count + 1 + c_f_index_array(c) - c_f_index_array(c-1)
!             do cf=(1+c_f_index_array(c-1)),c_f_index_array(c)
!                 f = c_f_obj_relation_array(cf)
!                 offset_count = offset_count + (f_p_index_array(f)-f_p_index_array(f-1))
!             enddo
!             write(2) offset_count
!         enddo
        
        
        
        close(2)
        
    end subroutine output_binary_internal
    
    subroutine output_ascii_internal
        implicit none
        
    end subroutine output_ascii_internal
    
    subroutine output_vtu_ascii
        implicit none
        
    end subroutine output_vtu_ascii
    
    subroutine output_vtu_binary_appended
        implicit none
		integer(KIND=INT32)             :: i, i1, i2, c, f, e, offset_count, cf, fe, p, fp1, fp2
		integer(KIND=INT8 ),allocatable :: types(:)
        integer(KIND=INT64)			    :: offset, face_offset
        real(KIND=REAL64),dimension(npoin,3) :: outvar, tot_proj
        character(:), allocatable		:: outstring
        character(16)                   :: offset_text, TXTnele, TXTnpoin
        
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
    outstring = outstring//'<DataArray type="Float64" NumberOfComponents="3"  Name="Points" format="appended" offset="'//offset_text//'" ></DataArray>'//NEW_LINE('')//'&
&</Points>'//NEW_LINE('')//'&
&<Cells>'//NEW_LINE('')
        offset = offset + 8 + npoin*3*real_length
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

        offset = offset + 8 + nele*4
        write(offset_text, '(i16)')  offset
    outstring = outstring//'<DataArray type="Float64" NumberOfComponents="3" Name="sn" format="appended" offset="'//offset_text//'" ></DataArray>'//NEW_LINE('')
        offset = offset + 8 + npoin*3*real_length
        write(offset_text, '(i16)')  offset
    outstring = outstring//'<DataArray type="Float64" NumberOfComponents="3" Name="sb" format="appended" offset="'//offset_text//'" ></DataArray>'//NEW_LINE('')
        offset = offset + 8 + npoin*3*real_length
        write(offset_text, '(i16)')  offset
    outstring = outstring//'<DataArray type="Float64" NumberOfComponents="3" Name="sbb" format="appended" offset="'//offset_text//'" ></DataArray>'//NEW_LINE('')
        offset = offset + 8 + npoin*3*real_length
        write(offset_text, '(i16)')  offset
    outstring = outstring//'<DataArray type="Float64" NumberOfComponents="3" Name="projections" format="appended" offset="'//offset_text//'" ></DataArray>'//NEW_LINE('')

        offset = offset + 8 + npoin*3*real_length
        write(offset_text, '(i16)')  offset
    outstring = outstring//'<DataArray type="Float64" Name="volume" format="appended" offset="'//offset_text//'" ></DataArray>'//NEW_LINE('')
    
        offset = offset + 8 + npoin*real_length
        write(offset_text, '(i16)')  offset
    outstring = outstring//'<DataArray type="Float64" NumberOfComponents="3" Name="normals" format="appended" offset="'//offset_text//'" ></DataArray>'//NEW_LINE('')
    
        outstring = outstring//'</PointData>'//NEW_LINE('')//'&
&</Piece>'//NEW_LINE('')//'&
&</UnstructuredGrid>'//NEW_LINE('')//'&
&<AppendedData encoding="raw">'//&
&'_'
        
        open(1,file='vtu_mesh.vtu',access='stream',status='replace',convert='little_endian')
        
        write(1) outstring
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !offsets
        offset = 8+npoin*3*real_length
        write(1) offset
        do i=1,npoin
            write(1) coords(i,:)
		enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !types
        offset = 8+nele*4
        write(1) offset, c_p_index_array(1:nele)
    
        allocate(types(nele))
        types = 42
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !connectivity
        offset = 8+nele
        write(1) offset, types(:)
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !faces
        offset = 8+c_p_sum*4
        write(1) offset, c_p_obj_relation_array(:)-1
    
        write(1) 8+face_offset
        do c=1,nele
            ! number of faces
            write(1) (c_f_index_array(c)-c_f_index_array(c-1)) 
            do cf=(1+c_f_index_array(c-1)),c_f_index_array(c)
                f = c_f_obj_relation_array(cf)
                fp1 = f_p_index_array(f-1)
                fp2 = f_p_index_array(f)
                ! number of points in face, points in face
                write(1) (fp2-fp1) , f_p_obj_relation_array((fp1+1):fp2)-1
                
            enddo
        enddo
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !faceoffsets
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
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !sn
        offset = 8+npoin*3*real_length
        write(1) offset
        
        tot_proj = 0.0 
        
        outvar = 0.0
        do i=1,i_nedge
        
            i1 = e_p_obj_relation_array(e_p_index_array(e_internal_indexing_array(i))-1)
            i2 = e_p_obj_relation_array(e_p_index_array(e_internal_indexing_array(i)))
            
            outvar(i1,:) = outvar(i1,:) + sn(i,:)
            outvar(i2,:) = outvar(i2,:) - sn(i,:)
            
        enddo
        
        do i=1,npoin
            write(1) outvar(i,:)
		enddo
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !sb
        offset = 8+npoin*3*real_length
        write(1) offset
        
        tot_proj = tot_proj + outvar
        outvar = 0.0
        do i=1,b_nedge
            i1 = e_p_obj_relation_array(e_p_index_array(e_bound_indexing_array(i))-1)
            i2 = e_p_obj_relation_array(e_p_index_array(e_bound_indexing_array(i)))
            
            outvar(i1,:) = outvar(i1,:) + sb(i,:)
            outvar(i2,:) = outvar(i2,:) - sb(i,:)
            
        enddo
        
        do i=1,npoin
            write(1) outvar(i,:)
		enddo
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !sbb
        offset = 8+npoin*3*real_length
        write(1) offset
        
        tot_proj = tot_proj + outvar
        outvar = 0.0
        do i=1,b_npoin
            i1 = p_bound_indexing_array(i)
            
            outvar(i1,:) = outvar(i1,:) + sbb(i,:)
            
        enddo
        
        tot_proj = tot_proj + outvar
        
        do i=1,npoin
            write(1) outvar(i,:)
		enddo
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !total projections
        offset = 8+npoin*3*real_length
        write(1) offset
        
        do i=1,npoin
            write(1) tot_proj(i,:)
		enddo
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !volumes
        offset = 8+npoin*real_length
        write(1) offset
        
        write(1) vol
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !normals
        offset = 8+npoin*3*real_length
        write(1) offset
        
        outvar = 0.0
        do i=1,b_npoin
            i1 = p_bound_indexing_array(i)
            outvar(i1,:) = outvar(i1,:) + p_normal_vectors(i,:)
        enddo
        
        do i=1,npoin
            write(1) outvar(i,:)
		enddo
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        outstring = NEW_LINE('')//'&
&</AppendedData>'//NEW_LINE('')//'&
&</VTKFile>'
    
        write(1) outstring
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        close(1)
        
        print*, NEW_LINE(''),NEW_LINE(''),NEW_LINE(''),NEW_LINE(''),NEW_LINE('')
        print*, 'FINISHED WRITING VTU BINARY FILE'
        print*, NEW_LINE(''),NEW_LINE(''),NEW_LINE(''),NEW_LINE(''),NEW_LINE('')
        
    end subroutine output_vtu_binary_appended
    
end module data_outputting_serial
