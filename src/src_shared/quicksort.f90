
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
        
        select type (array)
            type is (integer(KIND=int16))
                call qsort_i6(array, i, j)
                
            type is (integer(KIND=int32))
                call qsort_i32(array, i, j)
                
            type is (integer(KIND=int64))
                call qsort_i64(array, i, j)
                
            type is (real(KIND=real32))
                call qsort_r32(array, i, j)
                
            type is (real(KIND=real64))
                call qsort_r64(array, i, j)
                
            class default
                print*,"MISSING TYPING DEFINITON IN quicksort_unindexed"
                STOP
        end select
        
!         select type (array)
!             type is (integer(KIND=int16))
!                 call qsort_swap_i6(array, i, j)
!                 if (first < i-1) call qsort_i6(array, first, i-1)
!                 if (j+1 < last)  call qsort_i6(array, j+1, last)
!                 
!             type is (integer(KIND=int32))
!                 call qsort_swap_i32(array, i, j)
!                 if (first < i-1) call qsort_i32(array, first, i-1)
!                 if (j+1 < last)  call qsort_i32(array, j+1, last)
!                 
!             type is (integer(KIND=int64))
!                 call qsort_swap_i64(array, i, j)
!                 if (first < i-1) call qsort_i64(array, first, i-1)
!                 if (j+1 < last)  call qsort_i64(array, j+1, last)
!                 
!             type is (real(KIND=real32))
!                 call qsort_swap_r32(array, i, j)
!                 if (first < i-1) call qsort_r32(array, first, i-1)
!                 if (j+1 < last)  call qsort_r32(array, j+1, last)
!                 
!             type is (real(KIND=real64))
!                 call qsort_swap_r64(array, i, j)
!                 if (first < i-1) call qsort_r64(array, first, i-1)
!                 if (j+1 < last)  call qsort_r64(array, j+1, last)
!                 
!             class default
!                 print*,"MISSING TYPING DEFINITON IN quicksort_unindexed"
!                 STOP
!         end select
        
    end subroutine quicksort_unindexed
    
    recursive subroutine quicksort_with_indexer(array, first_unbound , last_unbound, Sorted_index)
        
        implicit none
        class(*)                     :: array(:), first_unbound , last_unbound
        integer(KIND=int32)          :: first, last, i, j, indexer
        integer(KIND=int32),allocatable :: Sorted_index(:)

        first = type_bounding(first_unbound)
        last  = type_bounding(last_unbound)
        
        i = first
        j = last
        
        Sorted_index = (/(i,i=first,last)/)
        
        select type (array)
            type is (integer(KIND=int16))
                call qsort_i6(array, i, j, Sorted_index)
                
            type is (integer(KIND=int32))                
                print*, 'qsort indexer 1'
                print*, array
                print*, Sorted_index
                call qsort_i32(array, i, j, Sorted_index)
                print*, 'qsort indexer 2'
                print*, array
                print*, Sorted_index
                
            type is (integer(KIND=int64))
                call qsort_i64(array, i, j, Sorted_index)
                
            type is (real(KIND=real32))
                call qsort_r32(array, i, j, Sorted_index)
                
            type is (real(KIND=real64))
                call qsort_r64(array, i, j, Sorted_index)
                
            class default
                print*,"MISSING TYPING DEFINITON IN quicksort_with_indexer"
                STOP
        end select
        
    end subroutine quicksort_with_indexer
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TERTIARY QUICKSORT ROUTINES


    recursive subroutine quicksort_arr(array, first_unbound , last_unbound, secondary_array&
    &, third_array, fourth_array, fifth_array, sixth_array, seventh_array, eighth_array, ninth_array, tenth_array)
    
        
        implicit none
        class(*)                        :: array(:), first_unbound , last_unbound
        class(*)                        :: secondary_array(:)
        class(*),optional               :: third_array(:), fourth_array(:), fifth_array(:), sixth_array(:), seventh_array(:), eighth_array(:), ninth_array(:), tenth_array(:)
        integer(KIND=int32),allocatable :: Sorted_index(:)
        
        
        call quicksort_with_indexer(array, first_unbound , last_unbound, Sorted_index)
        
        call tertiary_array_swapper(secondary_array, Sorted_index)
        
        ! this is probably huge overkill, but it means that I could in theory dump an entire tensor in this and it'll handle it
        
        if (present(third_array))   call tertiary_array_swapper(third_array,  Sorted_index)
        
        if (present(fourth_array))  call tertiary_array_swapper(fourth_array, Sorted_index)
        
        if (present(fifth_array))   call tertiary_array_swapper(fifth_array,  Sorted_index)
        
        if (present(sixth_array))   call tertiary_array_swapper(sixth_array,  Sorted_index)
        
        if (present(seventh_array)) call tertiary_array_swapper(seventh_array,Sorted_index)
        
        if (present(eighth_array))  call tertiary_array_swapper(eighth_array, Sorted_index)
        
        if (present(ninth_array))   call tertiary_array_swapper(ninth_array,  Sorted_index)
        
        if (present(tenth_array))   call tertiary_array_swapper(tenth_array,  Sorted_index) 
        
    end subroutine quicksort_arr
    
    recursive subroutine quicksort_matrix(array, first_unbound , last_unbound, matrix_one, matrix_two, matrix_three)
    
        
        implicit none
        class(*)                        :: array(:), first_unbound , last_unbound
        class(*)            :: matrix_one(:,:)
        class(*),optional   :: matrix_two(:,:), matrix_three(:,:)
        integer(KIND=int32),allocatable :: Sorted_index(:)
        integer(KIND=int32)             :: i
        
        call quicksort_with_indexer(array, first_unbound , last_unbound, Sorted_index)
    
        do i=1,size(matrix_one,2)
            call tertiary_array_swapper(matrix_one(:,i),Sorted_index)
        enddo
        
        if (present(matrix_two)) then
            do i=1,size(matrix_two,2)
                call tertiary_array_swapper(matrix_one(:,i),Sorted_index)
            enddo
        endif
        
        if (present(matrix_three)) then
            do i=1,size(matrix_three,2)
                call tertiary_array_swapper(matrix_one(:,i),Sorted_index)
            enddo
        endif
        
    end subroutine quicksort_matrix
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TYPE-UNBOUND FUNCTION FOR CONVERTING TO INT32
    
    integer(KIND=int32) function type_bounding(input)
        implicit none
        class(*) :: input
    
        select type (input)
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
    
    recursive subroutine qsort_i6(array, first , last, sorted_index)
        
        implicit none
        integer(KIND=int16) :: array(:)
        integer(KIND=int32) first, last
        integer(KIND=int32) i, j, indexer
        integer(KIND=int32), optional :: sorted_index(:)
        
        i = first
        j = last
        
        call qsort_swap_i6(array, i, j)
        
        if (present(sorted_index)) then
            call indexer_swap(sorted_index,first,last,i,j)
            if (first < i-1) call qsort_i6(array, first, i-1, sorted_index)
            if (j+1 < last)  call qsort_i6(array, j+1, last,  sorted_index)
        else
            if (first < i-1) call qsort_i6(array, first, i-1)
            if (j+1 < last)  call qsort_i6(array, j+1, last)
        endif
            
    end subroutine qsort_i6
    
    recursive subroutine qsort_i32(array, first , last, sorted_index)
        
        implicit none
        integer(KIND=int32) :: array(:)
        integer(KIND=int32) first, last
        integer(KIND=int32) i, j, indexer
        integer(KIND=int32), optional :: sorted_index(:)
        
        i = first
        j = last
        
        call qsort_swap_i32(array, i, j)
        
        if (present(sorted_index)) then
            call indexer_swap(sorted_index,first,last,i,j)
            if (first < i-1) call qsort_i32(array, first, i-1, sorted_index)
            if (j+1 < last)  call qsort_i32(array, j+1, last,  sorted_index)
        else
            if (first < i-1) call qsort_i32(array, first, i-1)
            if (j+1 < last)  call qsort_i32(array, j+1, last)
        endif
    end subroutine qsort_i32
    
    recursive subroutine qsort_i64(array, first , last, sorted_index)
        
        implicit none
        integer(KIND=int64) :: array(:)
        integer(KIND=int32) first, last
        integer(KIND=int32) i, j, indexer
        integer(KIND=int32), optional :: sorted_index(:)
        
        i = first
        j = last
        
        call qsort_swap_i64(array, i, j)
        
        if (present(sorted_index)) then
            call indexer_swap(sorted_index,first,last,i,j)
            if (first < i-1) call qsort_i64(array, first, i-1, sorted_index)
            if (j+1 < last)  call qsort_i64(array, j+1, last,  sorted_index)
        else
            if (first < i-1) call qsort_i64(array, first, i-1)
            if (j+1 < last)  call qsort_i64(array, j+1, last)
        endif
    end subroutine qsort_i64
    
    recursive subroutine qsort_r32(array, first , last, sorted_index)
        
        implicit none
        real(KIND=real32) :: array(:)
        integer(KIND=int32) first, last
        integer(KIND=int32) i, j, indexer
        integer(KIND=int32), optional :: sorted_index(:)
        
        i = first
        j = last
        
        call qsort_swap_r32(array, i, j)
        
        if (present(sorted_index)) then
            call indexer_swap(sorted_index,first,last,i,j)
            if (first < i-1) call qsort_r32(array, first, i-1, sorted_index)
            if (j+1 < last)  call qsort_r32(array, j+1, last,  sorted_index)
        else
            if (first < i-1) call qsort_r32(array, first, i-1)
            if (j+1 < last)  call qsort_r32(array, j+1, last)
        endif
    end subroutine qsort_r32
    
    recursive subroutine qsort_r64(array, first , last, sorted_index)
        
        implicit none
        real(KIND=real64) :: array(:)
        integer(KIND=int32) first, last
        integer(KIND=int32) i, j, indexer
        integer(KIND=int32), optional :: sorted_index(:)
        
        i = first
        j = last
        
        call qsort_swap_r64(array, i, j)
        
        if (present(sorted_index)) then
            call indexer_swap(sorted_index,first,last,i,j)
            if (first < i-1) call qsort_r64(array, first, i-1, sorted_index)
            if (j+1 < last)  call qsort_r64(array, j+1, last,  sorted_index)
        else
            if (first < i-1) call qsort_r64(array, first, i-1)
            if (j+1 < last)  call qsort_r64(array, j+1, last)
        endif
    end subroutine qsort_r64
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TYPE-BOUND SWAPPING ROUTINES

    subroutine qsort_swap_i6(array,i,j)
        
        implicit none
        integer(KIND=int32) :: i, j
        real(kind=real32)   :: pivot
        integer(KIND=int16) ::  array(:)
        integer(KIND=int16) ::  swap_variable
        pivot = array( (i+j) / 2 )
        do
            do while (array(i) < pivot)
                i=i+1
            end do
            do while (pivot < array(j))
                j=j-1
            end do
            if (i >= j) exit
            
            swap_variable = array(i)
            array(i) = array(j)
            array(j) = swap_variable

            i=i+1
            j=j-1
        end do
    end subroutine qsort_swap_i6

    subroutine qsort_swap_i32(array,i,j)
        
        implicit none
        integer(KIND=int32) :: i, j
        real(kind=real32)   :: pivot
        integer(KIND=int32) ::  array(:)
        integer(KIND=int32) ::  swap_variable
        pivot = array( (i+j) / 2 )
        do
            do while (array(i) < pivot)
                i=i+1
            end do
            do while (pivot < array(j))
                j=j-1
            end do
            if (i >= j) exit
            
            swap_variable = array(i)
            array(i) = array(j)
            array(j) = swap_variable

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
            do while (array(i) < pivot)
                i=i+1
            end do
            do while (pivot < array(j))
                j=j-1
            end do
            if (i >= j) exit
            
            swap_variable = array(i)
            array(i) = array(j)
            array(j) = swap_variable

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
            do while (array(i) < pivot)
                i=i+1
            end do
            do while (pivot < array(j))
                j=j-1
            end do
            if (i >= j) exit
            
            swap_variable = array(i)
            array(i) = array(j)
            array(j) = swap_variable

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
            do while (array(i) < pivot)
                i=i+1
            end do
            do while (pivot < array(j))
                j=j-1
            end do
            if (i >= j) exit
            
            swap_variable = array(i)
            array(i) = array(j)
            array(j) = swap_variable

            i=i+1
            j=j-1
        end do
    end subroutine qsort_swap_r64
    
    subroutine indexer_swap(array, old_i, old_j, i, j)
        
        implicit none
        integer(KIND=int32) :: i, j, old_i, old_j
        integer(KIND=int32)  ::  array(:)
        integer(KIND=int32)  ::  swap_variable

! c         swap_variable = array(i-1)
! c         array(i-1) = array(j+1)
! c         array(j+1) = swap_variable
        
!         swap_variable = array(i)
        array(i-1) = array(old_i)
        array(j+1) = array(old_j)

    end subroutine indexer_swap

    subroutine tertiary_array_swapper(tertiary_array_temp, Sorted_index)
        implicit none
        class(*) :: tertiary_array_temp(:)
        integer(KIND=int16),allocatable :: temp_arrayi2(:)
        integer(KIND=int32),allocatable :: temp_arrayi4(:)
        integer(KIND=int64),allocatable :: temp_arrayi8(:)
        real(KIND=real32),allocatable   :: temp_arrayr4(:)
        real(KIND=real64),allocatable   :: temp_arrayr8(:)
        integer(KIND=int32), intent(in)  :: Sorted_index(:)
        integer(KIND=int32) :: i
        select type (tertiary_array_temp)
            type is (integer(KIND=int16))
                allocate(temp_arrayi2, source=tertiary_array_temp)
                tertiary_array_temp(:) = temp_arrayi2(Sorted_index(:))
            type is (integer(KIND=int32))
                allocate(temp_arrayi4, source=tertiary_array_temp)
!                 print*,temp_arrayi4
!                 print*,tertiary_array_temp
!                 print*,Sorted_index
!                 print*,'AAA'
!                 do i=1,size(tertiary_array_temp)
!                     temp_arrayi4(i) = tertiary_array_temp(Sorted_index(i))
!                 enddo
!                 tertiary_array_temp = temp_arrayi4
                tertiary_array_temp(:) = temp_arrayi4(Sorted_index(:))
!                 print*,temp_arrayi4
!                 print*,tertiary_array_temp
!                 print*,Sorted_index
                
            type is (integer(KIND=int64))
                allocate(temp_arrayi8, source=tertiary_array_temp)
                tertiary_array_temp(:) = temp_arrayi8(Sorted_index(:))
            type is (real(KIND=real32))
                allocate(temp_arrayr4, source=tertiary_array_temp)
                tertiary_array_temp(:) = temp_arrayr4(Sorted_index(:))
            type is (real(KIND=real64))
                allocate(temp_arrayr8, source=tertiary_array_temp)
                tertiary_array_temp(:) = temp_arrayr8(Sorted_index(:))
            class default
                print*,"MISSING TYPING DEFINITON IN tertiary_array_swapper"
                STOP
        end select
    
    end subroutine tertiary_array_swapper
    
    ! end contains
    
end module quicksort_module
