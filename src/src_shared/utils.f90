

module utils
    use iso_fortran_env
    implicit none
    
	private alignment_single, alignment_double
    
    interface alignment

        module procedure &
        & alignment_single, &
        & alignment_double
        
    end interface alignment
    
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
    
    !end contains
    
end module
