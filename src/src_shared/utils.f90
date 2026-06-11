

module utils
    use iso_fortran_env
    implicit none
    
    contains
    
!     subroutine sort_and_flag_duplicates(array, start_index, end_index)
    subroutine sort_and_flag_duplicates(array)!, start_index, end_index)
        use quicksort_module
        implicit none
        integer(KIND=INT32),allocatable :: array(:)
        integer(KIND=INT32) :: start_index, end_index, indexer
        !integer(KIND=INT32),allocatable :: temparray(:)
        
        !allocate(temparray,source=array(start_index:end_index))
        
        !call quicksort(temparray,1,size(temparray)) ! must be start_index because indexing arrays start at zero
        
        start_index = 1
        end_index = size(array)
        
        call quicksort(array,start_index,end_index) ! must be start_index because indexing arrays start at zero
        
        indexer = start_index + 1
        
        do 
            if (array(indexer-1) .eq. array(indexer)) array(indexer-1) = -1
            
            
            indexer = indexer + 1
            if (indexer.gt.end_index) exit
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
