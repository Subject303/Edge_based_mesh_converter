

module utils
    use iso_fortran_env
    implicit none
    
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
        print*, mask
        
        array=pack(array,mask)
        
        print*, array
        
    end subroutine remove_flagged_duplicates
    
    !end contains
    
end module
