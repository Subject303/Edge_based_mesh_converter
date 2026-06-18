

module volume_processing
    use iso_fortran_env
    use object_counts
    use raw_data
    use object_relation_data
    use centroid_data
    use boundary_data
    use quicksort_module
    ! need everything p much here.
    implicit none
    
        
    contains
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine internal_edge_volume_processing
        implicit none
        integer(KIND=INT32) :: e, i
        ! it's upsetting this is the easiest of the three jobs I gotta do
        
        ! so, step wise here's the process per edge
        ! 1:
        !   find all connected cells and faces
        ! 2:
        !   order an array of cells and faces
        !   hopefully this is not too dificult
        ! 3:
        !   call the processer subroutine and sum the projections and volumes
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    end subroutine internal_edge_volume_processing
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine cvolume(sn,vol1,vol2,i1,i2,i3,i4,i5)
        use shared
        real(KIND=REAL32) :: vol1,vol2
        real(KIND=REAL32),dimension(3) :: v13,v14,v15,v32,v42,v52,c1415,c4252
        real(KIND=REAL32),dimension(3) :: sn,i3,i4,i5
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
        
        vol1 = abs((c1415(1)*v13(1) ) - (c1415(2)*v13(2) ) + (c1415(3)*v13(3) ) )/6
        vol2 = abs((c1415(1)*v32(1) ) - (c1415(2)*v32(2) ) + (c1415(3)*v32(3) ) )/6

        do i=1,3
            v34(i) = i3(i) - i4(i) ! vector from 3 to 4
            v35(i) = i3(i) - i5(i) ! vector from 3 to 5
        enddo

        sn(1) = ((v34(2) * v35(3) - v34(3) * v35(2))/2) ! projection in the xx axis
        sn(2) = ((v34(3) * v35(1) - v34(1) * v35(3))/2) ! projection in the yy axis
        sn(3) = ((v34(1) * v35(2) - v34(2) * v35(1))/2) ! projection in the zz axis
        
        i=1
    !                               i4--------------i5
    !                                \  >>>>>>>   /
    !                                 \         /
    !                                  \      / 
    !                                   \   /
    !           i1 -------------------- i3 -------------------------i2
    end subroutine
    
    ! end contains
    
end module volume_processing
