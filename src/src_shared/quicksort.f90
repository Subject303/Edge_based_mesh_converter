
module quicksort_module
    use iso_fortran_env
    implicit none 
    
    private
    
    public quicksort, quicksort_with_indexer

    interface quicksort

        module procedure &
        & quicksort_unindexed, &
        & quicksort_arr  , quicksort_matrix
        
    end interface quicksort

    contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PRIMARY QUICKSORT ROUTINES

    recursive subroutine quicksort_unindexed(array, first_unbound , last_unbound)
        
        implicit none
        class(*)                     :: array(:), first_unbound , last_unbound
        integer(KIND=int32)          :: first, last, i, j

        first = type_bounding(first_unbound)
        last  = type_bounding(last_unbound)
        
        i = first
        j = last
        
        select type (a)
            type is (integer(KIND=int16))
                call qsort_swap_i16(array, i, j)
                if (first < i-1) call qsort_i16(array, first, i-1)
                if (j+1 < last)  call qsort_i16(array, j+1, last)
                
            type is (integer(KIND=int32))
                call qsort_swap_i32(array, i, j)
                if (first < i-1) call qsort_i32(array, first, i-1)
                if (j+1 < last)  call qsort_i32(array, j+1, last)
                
            type is (integer(KIND=int64))
                call qsort_swap_i64(array, i, j)
                if (first < i-1) call qsort_i64(array, first, i-1)
                if (j+1 < last)  call qsort_i64(array, j+1, last)
                
            type is (real(KIND=real32))
                call qsort_swap_r32(array, i, j)
                if (first < i-1) call qsort_r32(array, first, i-1)
                if (j+1 < last)  call qsort_r32(array, j+1, last)
                
            type is (real(KIND=real64))
                call qsort_swap_r64(array, i, j)
                if (first < i-1) call qsort_r64(array, first, i-1)
                if (j+1 < last)  call qsort_r64(array, j+1, last)
                
            class default
                print*,"MISSING TYPING DEFINITON IN quicksort_unindexed"
                STOP
        end select
        
    end subroutine quicksort_unindexed
    
    recursive subroutine quicksort_with_indexer(array, first_unbound , last_unbound, Sorted_index)
        
        implicit none
        class(*)                     :: array(:), first_unbound , last_unbound
        integer(KIND=int32)          :: first, last, i, j, indexer
        integer(KIND=int32),allocatable :: Sorted_index

        first = type_bounding(first_unbound)
        last  = type_bounding(last_unbound)
        
        i = first
        j = last
        
        Sorted_index = (/first:last/)
        
        select type (a)
            type is (integer(KIND=int16))
                call qsort_swap_i16(array, i, j)
                indexer         = Sorted_index(i)
                Sorted_index(i) = Sorted_index(j)
                Sorted_index(j) = indexer
                if (first < i-1) call qsort_i16(array, first, i-1)
                if (j+1 < last)  call qsort_i16(array, j+1, last)
                
            type is (integer(KIND=int32))
                call qsort_swap_i32(array, i, j)
                indexer         = Sorted_index(i)
                Sorted_index(i) = Sorted_index(j)
                Sorted_index(j) = indexer
                if (first < i-1) call qsort_i32(array, first, i-1)
                if (j+1 < last)  call qsort_i32(array, j+1, last)
                
            type is (integer(KIND=int64))
                call qsort_swap_i64(array, i, j)
                indexer         = Sorted_index(i)
                Sorted_index(i) = Sorted_index(j)
                Sorted_index(j) = indexer
                if (first < i-1) call qsort_i64(array, first, i-1)
                if (j+1 < last)  call qsort_i64(array, j+1, last)
                
            type is (real(KIND=real32))
                call qsort_swap_r32(array, i, j)
                indexer         = Sorted_index(i)
                Sorted_index(i) = Sorted_index(j)
                Sorted_index(j) = indexer
                if (first < i-1) call qsort_r32(array, first, i-1)
                if (j+1 < last)  call qsort_r32(array, j+1, last)
                
            type is (real(KIND=real64))
                call qsort_swap_r64(array, i, j)
                indexer         = Sorted_index(i)
                Sorted_index(i) = Sorted_index(j)
                Sorted_index(j) = indexer
                if (first < i-1) call qsort_r64(array, first, i-1)
                if (j+1 < last)  call qsort_r64(array, j+1, last)
                
            class default
                print*,"MISSING TYPING DEFINITON IN quicksort_with_indexer"
                STOP
        end select
        
    end subroutine quicksort_with_indexer
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TERTIARY QUICKSORT ROUTINES


    recursive subroutine quicksort_arr(array, first_unbound , last_unbound, secondary_array&
    &, third_array, fourth_array, fifth_array, sixth_array, seventh_array, eighth_array, ninth_array, tenth_array)
    
        
        implicit none
        class(*)                        :: array(:), secondary_array(:), first_unbound , last_unbound
        class(*),optional               :: third_array(:), fourth_array(:), fifth_array(:), sixth_array(:), seventh_array(:), eighth_array(:), ninth_array(:), tenth_array(:)
        class(*),allocatable            :: tertiary_array_temp(:)
        integer(KIND=int32),allocatable :: Sorted_index
        
        call quicksort_with_indexer(array, first_unbound , last_unbound, Sorted_index)
        
        allocate(tertiary_array_temp,source=secondary_array)
        secondary_array(:) = tertiary_array_temp(Sorted_index)
        
        ! this is probably huge overkill, but it means that I could in theory dump an entire tensor in this and it'll handle it
        
        if (present(third_array)) then
            deallocate(tertiary_array_temp)
            allocate(tertiary_array_temp,source=third_array)
            third_array(:) = tertiary_array_temp(Sorted_index)
        endif 
        
        if (present(fourth_array)) then
            deallocate(tertiary_array_temp)
            allocate(tertiary_array_temp,source=fourth_array)
            fourth_array(:) = tertiary_array_temp(Sorted_index)
        endif 
        
        if (present(fifth_array)) then
            deallocate(tertiary_array_temp)
            allocate(tertiary_array_temp,source=fifth_array)
            fifth_array(:) = tertiary_array_temp(Sorted_index)
        endif 
        
        if (present(sixth_array)) then
            deallocate(tertiary_array_temp)
            allocate(tertiary_array_temp,source=sixth_array)
            sixth_array(:) = tertiary_array_temp(Sorted_index)
        endif 
        
        if (present(seventh_array)) then
            deallocate(tertiary_array_temp)
            allocate(tertiary_array_temp,source=seventh_array)
            seventh_array(:) = tertiary_array_temp(Sorted_index)
        endif 
        
        if (present(eighth_array)) then
            deallocate(tertiary_array_temp)
            allocate(tertiary_array_temp,source=eighth_array)
            eighth_array(:) = tertiary_array_temp(Sorted_index)
        endif 
        
        if (present(ninth_array)) then
            deallocate(tertiary_array_temp)
            allocate(tertiary_array_temp,source=ninth_array)
            ninth_array(:) = tertiary_array_temp(Sorted_index)
        endif 
        
        if (present(tenth_array)) then
            deallocate(tertiary_array_temp)
            allocate(tertiary_array_temp,source=tenth_array)
            tenth_array(:) = tertiary_array_temp(Sorted_index)
        endif 
        
    end subroutine quicksort_arr
    
    recursive subroutine quicksort_matrix(array, first_unbound , last_unbound, matrix_one, matrix_two, matrix_three)
    
        
        implicit none
        class(*)                        :: array(:), matrix_one(:,:), first_unbound , last_unbound
        class(*),optional               :: matrix_two(:,:), matrix_three(:,:)
        class(*),allocatable            :: tertiary_array_temp(:)
        integer(KIND=int32),allocatable :: Sorted_index
        integer(KIND=int32)             :: i
        
        call quicksort_with_indexer(array, first_unbound , last_unbound, Sorted_index)
    
        do i=1,size(matrix_one,2)
            allocate(tertiary_array_temp,source=matrix_one(:,i))
            matrix_one(:,i) = tertiary_array_temp(Sorted_index)
        enddo
        
        if (present(matrix_two)) then
            do i=1,size(matrix_two,2)
                allocate(tertiary_array_temp,source=matrix_two(:,i))
                matrix_two(:,i) = tertiary_array_temp(Sorted_index)
            enddo
        endif
        
        if (present(matrix_three)) then
            do i=1,size(matrix_three,2)
                allocate(tertiary_array_temp,source=matrix_three(:,i))
                matrix_three(:,i) = tertiary_array_temp(Sorted_index)
            enddo
        endif
    
    end subroutine quicksort_matrix
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TYPE-UNBOUND FUNCTION FOR CONVERTING TO INT32
    
    integer(KIND=int32) function type_bounding(input)
        class(*) :: input
        
        implicit none
    
        select type (a)
            type is (integer(KIND=int16))
                type_bounding = input
            type is (integer(KIND=int32))
                type_bounding = input
            type is (integer(KIND=int64))
                type_bounding = input
            type is (real(KIND=real32))
                type_bounding = input
            type is (real(KIND=real64))
                type_bounding = input
            class default
                print*,"MISSING TYPING DEFINITON IN bound_typing"
                STOP
        end select
    end function type_bounding
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TYPE-BOUND QSORT ROUTINES
    
    recursive subroutine qsort_i16(array, first , last)
        
        implicit none
        integer(KIND=int16) :: array(:)
        real(kind=real32)   pivot
        integer(KIND=int32) first, last
        integer(KIND=int32) i, j
        
        i = first
        j = last
        
        call qsort_swap_i16(array, i, j)
        
        if (first < i-1) call qsort_i16(array, first, i-1)
        if (j+1 < last)  call qsort_i16(array, j+1, last)
    end subroutine qsort_i16
    
    recursive subroutine qsort_i32(array, first , last)
        
        implicit none
        integer(KIND=int32) :: array(:)
        real(kind=real32)   pivot
        integer(KIND=int32) first, last
        integer(KIND=int32) i, j
        
        i = first
        j = last
        
        call qsort_swap_i32(array, i, j)
        
        if (first < i-1) call qsort_i32(array, first, i-1)
        if (j+1 < last)  call qsort_i32(array, j+1, last)
    end subroutine qsort_i32
    
    recursive subroutine qsort_i64(array, first , last)
        
        implicit none
        integer(KIND=int64) :: array(:)
        real(kind=real32)   pivot
        integer(KIND=int32) first, last
        integer(KIND=int32) i, j
        
        i = first
        j = last
        
        call qsort_swap_i64(array, i, j)
        
        if (first < i-1) call qsort_i64(array, first, i-1)
        if (j+1 < last)  call qsort_i64(array, j+1, last)
    end subroutine qsort_i64
    
    recursive subroutine qsort_r32(array, first , last)
        
        implicit none
        real(KIND=real32) :: array(:)
        real(kind=real32)   pivot
        integer(KIND=int32) first, last
        integer(KIND=int32) i, j
        
        i = first
        j = last
        
        call qsort_swap_r32(array, i, j)
        
        if (first < i-1) call qsort_r32(array, first, i-1)
        if (j+1 < last)  call qsort_r32(array, j+1, last)
    end subroutine qsort_r32
    
    recursive subroutine qsort_r64(array, first , last)
        
        implicit none
        real(KIND=real64) :: array(:)
        real(kind=real32)   pivot
        integer(KIND=int32) first, last
        integer(KIND=int32) i, j
        
        i = first
        j = last
        
        call qsort_swap_r64(array, i, j)
        
        if (first < i-1) call qsort_r64(array, first, i-1)
        if (j+1 < last)  call qsort_r64(array, j+1, last)
    end subroutine qsort_r64
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TYPE-BOUND SWAPPING ROUTINES

    subroutine qsort_swap_i16(array,i,j)
        
        implicit none
        integer(KIND=int32) :: i, j
        real(kind=real32)   :: pivot
        integer(KIND=int16) ::  array(:)
        integer(KIND=int16) ::  swap_variable
        pivot = array( (i+j) / 2 )
        do
            do while (array(i) < x)
                i=i+1
            end do
            do while (x < array(j))
                j=j-1
            end do
            if (i >= j) exit
            
            swap_variable = array(i1)
            array(i1) = array(i2)
            array(i2) = swap_variable

            i=i+1
            j=j-1
        end do
    end subroutine qsort_swap_i16

    subroutine qsort_swap_i32(array,i,j)
        
        implicit none
        integer(KIND=int32) :: i, j
        real(kind=real32)   :: pivot
        integer(KIND=int32) ::  array(:)
        integer(KIND=int32) ::  swap_variable
        pivot = array( (i+j) / 2 )
        do
            do while (array(i) < x)
                i=i+1
            end do
            do while (x < array(j))
                j=j-1
            end do
            if (i >= j) exit
            
            swap_variable = array(i1)
            array(i1) = array(i2)
            array(i2) = swap_variable

            i=i+1
            j=j-1
        end do
    end subroutine qsort_swap_i32

    subroutine qsort_swap_i64(array,i,j)
        
        implicit none
        integer(KIND=int32) :: i, j
        real(kind=real32)   :: pivot
        integer(KIND=int64) ::  array(:)
        integer(KIND=int64) ::  swap_variable
        pivot = array( (i+j) / 2 )
        do
            do while (array(i) < x)
                i=i+1
            end do
            do while (x < array(j))
                j=j-1
            end do
            if (i >= j) exit
            
            swap_variable = array(i1)
            array(i1) = array(i2)
            array(i2) = swap_variable

            i=i+1
            j=j-1
        end do
    end subroutine qsort_swap_i64

    subroutine qsort_swap_r32(array,i,j)
        
        implicit none
        integer(KIND=int32) :: i, j
        real(kind=real32)   :: pivot
        real(KIND=real32) ::  array(:)
        real(KIND=real32) ::  swap_variable
        pivot = array( (i+j) / 2 )
        do
            do while (array(i) < x)
                i=i+1
            end do
            do while (x < array(j))
                j=j-1
            end do
            if (i >= j) exit
            
            swap_variable = array(i1)
            array(i1) = array(i2)
            array(i2) = swap_variable

            i=i+1
            j=j-1
        end do
    end subroutine qsort_swap_r32

    subroutine qsort_swap_r64(array,i,j)
        
        implicit none
        integer(KIND=int32) :: i, j
        real(kind=real32)   :: pivot
        real(KIND=real64)  ::  array(:)
        real(KIND=real64)  ::  swap_variable
        pivot = array( (i+j) / 2 )
        do
            do while (array(i) < x)
                i=i+1
            end do
            do while (x < array(j))
                j=j-1
            end do
            if (i >= j) exit
            
            swap_variable = array(i1)
            array(i1) = array(i2)
            array(i2) = swap_variable

            i=i+1
            j=j-1
        end do
    end subroutine qsort_swap_i64

    ! end contains
    
end module quicksort_module
