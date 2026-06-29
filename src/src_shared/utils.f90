

module utils
    use iso_fortran_env
    implicit none
    
	private alignment_single, alignment_double
    
    interface alignment

        module procedure &
        & alignment_single, &
        & alignment_double
        
    end interface alignment
    
    interface planar_alignment

        module procedure &
        & planar_alignment_single, &
        & planar_alignment_double
        
    end interface planar_alignment
    
    contains
    
    subroutine sort_and_flag_duplicates(array)
        use quicksort_module
        implicit none
        integer(KIND=INT32) :: array(:), start_index, end_index, indexer
        
        
        !call quicksort(temparray,1,size(temparray)) ! must be start_index because indexing arrays start at zero
        
        call quicksort(array,1 ,size(array) ) ! must be start_index because indexing arrays start at zero
        
        do indexer=2,size(array)
            if (array(indexer-1) .eq. array(indexer)) array(indexer-1) = -1
        enddo
        
    end subroutine sort_and_flag_duplicates
    
    subroutine remove_flagged_duplicates(array)
        implicit none
        integer(KIND=INT32),allocatable :: array(:)
        logical,allocatable :: mask(:)
        
        allocate(mask(size(array)))
        
        mask=(/ (array(:).ge.0) /)
        
        array=pack(array,mask)
        
    end subroutine remove_flagged_duplicates
    
    real(KIND=REAL32) function alignment_single(v1,v2)
        implicit none
        real(KIND=REAL32) :: v1(3), v2(3)
        
        alignment_single = (v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3))/ (sqrt(v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3)) * sqrt(v2(1)*v2(1) + v2(2)*v2(2) + v2(3)*v2(3)))
        
    end function alignment_single
    
    real(KIND=REAL64) function alignment_double(v1,v2)
        implicit none
        real(KIND=REAL64) :: v1(3), v2(3)
        
        alignment_double = (v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3))/ (sqrt(v1(1)*v1(1) + v1(2)*v1(2) + v1(3)*v1(3)) * sqrt(v2(1)*v2(1) + v2(2)*v2(2) + v2(3)*v2(3)))
        
    end function alignment_double
    
    real(KIND=REAL32) function planar_alignment_single(v1,v2,v3)
        implicit none
        real(KIND=REAL32) :: v1(3), v2(3), v3(3), mag(3), det, dot
        
        mag(1) = sqrt(v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3))
        v1(:) = v1(:) / mag(1)
        
        mag(2) = sqrt(v2(1)*v2(1)+v2(2)*v2(2)+v2(3)*v2(3))
        v2(:) = v2(:) / mag(2)
        
        mag(3) = sqrt(v3(1)*v3(1)+v3(2)*v3(2)+v3(3)*v3(3))
        v3(:) = v3(:) / mag(3)
        
        dot = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
        det = v1(1)*v2(2)*v3(3) + v2(1)*v3(2)*v1(3) + v3(1)*v1(2)*v2(3) - v1(3)*v2(2)*v3(1) - v2(3)*v3(2)*v1(1) - v3(3)*v1(2)*v2(1)
                
        planar_alignment_single = ATAN(det,dot)
        
    end function planar_alignment_single
    
    real(KIND=REAL64) function planar_alignment_double(v1,v2,v3)
        implicit none
        real(KIND=REAL64) :: v1(3), v2(3), v3(3), mag(3), det, dot
        
        mag(1) = sqrt(v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3))
        v1(:) = v1(:) / mag(1)
        
        mag(2) = sqrt(v2(1)*v2(1)+v2(2)*v2(2)+v2(3)*v2(3))
        v2(:) = v2(:) / mag(2)
        
        mag(3) = sqrt(v3(1)*v3(1)+v3(2)*v3(2)+v3(3)*v3(3))
        v3(:) = v3(:) / mag(3)
        
        dot = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
        det = v1(1)*v2(2)*v3(3) + v2(1)*v3(2)*v1(3) + v3(1)*v1(2)*v2(3) - v1(3)*v2(2)*v3(1) - v2(3)*v3(2)*v1(1) - v3(3)*v1(2)*v2(1)
                
        planar_alignment_double = ATAN(det,dot)
        
    end function planar_alignment_double
    
    !end contains
    
end module
