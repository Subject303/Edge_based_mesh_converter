

module utils
    use iso_fortran_env
    implicit none
    
    contains
    
    subroutine sort_and_flag_duplicates(array, start_index, end_index)
        use quicksort_module
        implicit none
        integer(KIND=INT32) :: array(:), start_index, end_index, indexer
        integer(KIND=INT32),allocatable :: temparray(:)
        
        allocate(temparray,source=array(start_index:end_index))
        
        print*, start_index, end_index, size(temparray)
        
        call quicksort(temparray,1,size(temparray)) ! must be start_index because indexing arrays start at zero
        
        
        indexer = start_index + 1
        
        do 
            if (array(indexer-1) .eq. array(indexer)) array(indexer-1) = -1
            
            
            indexer = indexer + 1
            if (indexer.gt.end_index) exit
        enddo
        
    end subroutine sort_and_flag_duplicates
    
    subroutine remove_flagged_duplicates(array)
        implicit none
        integer(KIND=INT32) :: array(:)
        logical,allocatable :: mask(:)
        
        allocate(mask(size(array)))
        
        mask=(/ (array(:).ge.0) /)
        array=pack(array,mask)
        
    end subroutine remove_flagged_duplicates
    
    !end contains
    
end module
