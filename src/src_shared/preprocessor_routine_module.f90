


module preprocessor_routine_module
    use iso_fortran_env
    use object_counts
    use raw_data
    use object_relation_data
    use centroid_data
    use quicksort_module
    implicit none
    
    contains
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!  centroid routines below ::
    
!     subroutine qsort_test
!         implicit none
!         
!         integer(KIND=INT32),allocatable :: temp(:), temp2(:)
!         print*, 'beginning calculating edge centroids'
!         
!         temp = (/3,5,4,3,2,1/)
!         
!         temp2 = (/77,2,4,6,8,10/)
!         print*, 'aaa'
!         print*, temp
!         print*, temp2
!         print*, 'aaa'
!         call quicksort(temp, 1, 6, temp2)
!         print*, 'aaa'
!         print*, temp
!         print*, temp2
!         print*, 'aaa'
!         
!         stop
! 
!         
!     end subroutine qsort_test
    
    subroutine calc_e_centroid
        implicit none
        integer(KIND=INT32) :: e, i, ep_start, ep_end
        
        print*, 'beginning calculating edge centroids'
        
        allocate(e_centroid(nedge,3))
        
        e_p_sum = nedge*2
        e_p_index_array = (/(e,e=0,nedge*2,2)/) 
        
        do e=1,nedge
            ep_start = (1 + e_p_index_array(e-1))
            ep_end   = e_p_index_array(e)
            
            do i=1,3
                e_centroid(e,i) = sum(coords(e_p_obj_relation_array(ep_start:ep_end),i)) / 2
            enddo
        enddo
        
        print*, 'finished calculating edge centroids'
        
    end subroutine calc_e_centroid
    
    subroutine calc_f_centroid
        implicit none
        integer(KIND=INT32) :: f, i, fp_start, fp_end, number_of_points
        
        print*, 'beginning calculating face centroids'
        
        allocate(f_centroid(nface,3))

        do f=1,nface
            fp_start = (1 + f_p_index_array(f-1))
            fp_end   = f_p_index_array(f)
            number_of_points = 1 + fp_end - fp_start
            
            do i=1,3
                f_centroid(f,i) = sum(coords(f_p_obj_relation_array(fp_start:fp_end),i)) / number_of_points
            enddo
            
        enddo
        
        print*, 'finished calculating face centroids'
        
    end subroutine calc_f_centroid    
    
    subroutine calc_c_centroid
        implicit none
        integer(KIND=INT32) :: c, i, cp_start, cp_end, number_of_points
        
        
        print*, 'beginning calculating cell centroids'
        
        allocate(c_centroid(nele,3))
        
        do c=1,nele
            cp_start = (1 + c_p_index_array(c-1))
            cp_end   = c_p_index_array(c)
            number_of_points = 1 + cp_end - cp_start
            
            do i=1,3
                c_centroid(c,i) = sum(coords(c_p_obj_relation_array(cp_start:cp_end),i)) / number_of_points
            enddo
            
        enddo
        
        print*, 'finished calculating cell centroids'
        
    end subroutine calc_c_centroid
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!  cell edge routine below ::
    
    
    subroutine c_e_preprocess
        use utils
        implicit none
        integer(KIND=INT32) :: c, cf_index, cf_start_index, cf_end_index, f, fe_index, fe_start_index, fe_end_index, ce_index, n_con_edges, ce_start_index, ce_end_index, backward_index_duplicates(0:nele), i
        
        
        print *, 'beginning processing cell edge relation array'
        
!         print *, c_f_index_array
!         print *, 'AAAAAAAAAAAAAA1'
!         print *, c_f_obj_relation_array
!         print *, 'AAAAAAAAAAAAAA'
!         print *, f_e_index_array
        
        
        c_e_sum = 0
        
        allocate(c_e_index_array(0:nele))
        c_e_index_array(0) = 0
        
        
        !print*, size(c_f_obj_relation_array), c_f_sum, maxval(c_f_index_array), size(c_f_index_array), nele
        
        do c=1,nele
        
            cf_start_index = 1 + c_f_index_array(c-1)
            cf_end_index   = c_f_index_array(c)
            
            !print*, 'CELL ' ,c,c_f_index_array(c-1),c_f_index_array(c)
            
            n_con_edges = 0
            do cf_index = cf_start_index, cf_end_index
                
                if (cf_index .gt. 3111791) print*, cf_index, cf_start_index, cf_end_index
                
                f = c_f_obj_relation_array(cf_index)
                
                n_con_edges =  n_con_edges + f_e_index_array(f) - f_e_index_array(f-1) ! all index arrays start at 0
            enddo
            
            c_e_index_array(c) = n_con_edges
        enddo
        
        do c=1,nele
            c_e_index_array(c) = c_e_index_array(c-1) + c_e_index_array(c)
        enddo
        
        c_e_sum = c_e_index_array(nele)
        
        allocate( c_e_obj_relation_array(c_e_sum) )
        
        ! finished establishing the bounds of c_e
        
        ce_index = 0
        
        do c=1,nele
            n_con_edges = c_e_index_array(c) - c_e_index_array(c-1)
            
            
            cf_start_index = 1 + c_f_index_array(c-1)
            cf_end_index   = c_f_index_array(c)
            
            do cf_index=cf_start_index, cf_end_index
            
                f = c_f_obj_relation_array(cf_index)
                
                fe_start_index = 1 + f_e_index_array(f-1)
                fe_end_index   = f_e_index_array(f)
                
                do fe_index=fe_start_index, fe_end_index
                
                    ce_index = ce_index + 1
                    
                    c_e_obj_relation_array(ce_index) = f_e_obj_relation_array(fe_index)
                    
                enddo
            enddo
        enddo
        
        backward_index_duplicates = 0
        
        do c=1,nele
            ce_start_index = 1 + c_e_index_array(c-1)
            ce_end_index   = c_e_index_array(c)
            
            call sort_and_flag_duplicates(c_e_obj_relation_array(ce_start_index:ce_end_index))
            
            
            do i=ce_start_index,ce_end_index
                if (c_e_obj_relation_array(i) .lt. 0) backward_index_duplicates(c) = backward_index_duplicates(c) + 1 
            enddo
            
        enddo
        
        do c=1,nele
            backward_index_duplicates(c) = backward_index_duplicates(c-1) + backward_index_duplicates(c)
        enddo
        
        c_e_index_array(:) = c_e_index_array(:) - backward_index_duplicates(:)
        c_e_sum = c_e_index_array(nele)
        
        call remove_flagged_duplicates(c_e_obj_relation_array)
        
        print *, 'finished processing cell edge relation array'
        
    end subroutine c_e_preprocess
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!  obj relation inversion routines below ::
    
    !!!!!!!!!!! raw connectivity inversions ::
    
    subroutine p_c_preprocess
    
        print *, 'beginning processing cell connectivity array inversion'
        call obj_relation_inverter(nele, npoin, c_p_sum, c_p_index_array, c_p_obj_relation_array, p_c_sum, p_c_index_array, p_c_obj_relation_array)
        print *, 'finished processing cell connectivity array inversion'
        
    end subroutine p_c_preprocess
    
    
    subroutine p_f_preprocess
    
        print *, 'beginning processing face connectivity array inversion'
        call obj_relation_inverter(nface, npoin, f_p_sum, f_p_index_array, f_p_obj_relation_array, p_f_sum, p_f_index_array, p_f_obj_relation_array)
        print *, 'finished processing face connectivity array inversion'
        
    end subroutine p_f_preprocess
    
    subroutine p_e_preprocess
    
        print *, 'beginning processing edge connectivity array inversion'
        call obj_relation_inverter(nedge, npoin, e_p_sum, e_p_index_array, e_p_obj_relation_array, p_e_sum, p_e_index_array, p_e_obj_relation_array)
        print *, 'finished processing edge connectivity array inversion'
        
    end subroutine p_e_preprocess
    
    !!!!!!!!!!! e f c obj relation inversions ::
    
    subroutine f_c_preprocess
        
        print *, 'beginning processing face cell relation array inversion'
        call obj_relation_inverter(nele, nface, c_f_sum, c_f_index_array, c_f_obj_relation_array, f_c_sum, f_c_index_array, f_c_obj_relation_array)
        print *, 'finished processing face cell relation array inversion'
        
    end subroutine f_c_preprocess
    
    subroutine e_f_preprocess
    
        print *, 'beginning processing edge face relation array inversion'
        call obj_relation_inverter(nface, nedge, f_e_sum, f_e_index_array, f_e_obj_relation_array, e_f_sum, e_f_index_array, e_f_obj_relation_array)
        print *, 'finished processing edge face relation array inversion'
        
    end subroutine e_f_preprocess
    
    subroutine e_c_preprocess
        
        print *, 'beginning processing edge cell relation array inversion'
        call obj_relation_inverter(nele, nedge, c_e_sum, c_e_index_array, c_e_obj_relation_array, e_c_sum, e_c_index_array, e_c_obj_relation_array)
        print *, 'finished processing edge cell relation array inversion'
        
    end subroutine e_c_preprocess
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!! inverter routine below ::
    
    subroutine obj_relation_inverter(forward_leading_obj_count, backward_leading_obj_count, forward_sum, forward_index, forward_obj_relation_array, backward_sum, backward_index, backward_obj_relation_array)
        use utils
        implicit none
    
        integer(KIND=INT32),allocatable :: forward_index(:), backward_index(:), forward_obj_relation_array(:), backward_obj_relation_array(:)
        integer(KIND=INT32),allocatable :: y_index(:), backward_index_duplicates(:)
        integer(KIND=INT32) :: x, y, x_count, total_x_count, forward_sum, backward_sum, xy_start_index, xy_end_index, forward_leading_obj_count, backward_leading_obj_count, i
        
        
        allocate(backward_obj_relation_array(forward_sum), backward_index(0:backward_leading_obj_count), backward_index_duplicates(0:backward_leading_obj_count) )
        allocate(y_index,source=forward_obj_relation_array)
        
        !do i=1,forward_leading_obj_count
        !    print*,forward_obj_relation_array( (forward_index(i-1)+1) : forward_index(i) )
        !enddo
        
        !print*,forward_obj_relation_array
        
        do x=1, forward_leading_obj_count
            xy_start_index = 1 + forward_index(x-1)
            xy_end_index   = forward_index(x)

            !print*, xy_start_index, xy_end_index, x, forward_leading_obj_count, backward_leading_obj_count
            
            backward_obj_relation_array(xy_start_index:xy_end_index) = x

        enddo

        call quicksort(y_index , 1 , forward_sum , backward_obj_relation_array)
        
        !print*,y_index
        !print*,backward_obj_relation_array
        
		
        backward_index = 0
        backward_index_duplicates = 0
        total_x_count = 0
        x_count=0
        do y=1, backward_leading_obj_count
        
            ! count the number of cells adjacent to each point
			
            !print*, y , y_index(x_count + 1)
            
            do while (y .eq. y_index(x_count + 1))
                x_count = x_count + 1
                if (x_count .eq. forward_sum) exit
            enddo
            
            ! adjust x_count
			
            backward_index(y) = x_count
            
            ! check if there exists any possible duplicates
            if (backward_index(y-1)+1 .ne. x_count) then
            
                ! flag duplicates
                
                call sort_and_flag_duplicates(backward_obj_relation_array((backward_index(y-1)+1):x_count))
                
                do i=(backward_index(y-1)+1),x_count
                    if (backward_obj_relation_array(i) .lt. 0) backward_index_duplicates(y) = backward_index_duplicates(y) + 1 
                enddo
                
                
            endif
            
                
            ! I'ma be real here
            ! I'm pretty smug about this loop
            ! update: less so now

        enddo
        
        do i=1,backward_leading_obj_count
            backward_index_duplicates(i) = backward_index_duplicates(i-1) + backward_index_duplicates(i)
        enddo
        
        backward_index(:) = backward_index(:) - backward_index_duplicates(:)

        call remove_flagged_duplicates(backward_obj_relation_array)

        backward_sum = backward_index(backward_leading_obj_count)
        
        ! based on ::
!         integer(KIND=INT32),allocatable :: p_index, p_c_index_array_duplicates
!         integer(KIND=INT32) :: c, p, c_count, total_c_count
! 
!         allocate(p_c_obj_relation_array(c_p_sum), p_c_index_array(0:npoin), p_c_index_array_duplicates(0:npoin) )
!         allocate(p_index,source=c_p_obj_relation_array)
! 
!         do c=1,nele
!             cp_start_index = 1 + c_p_index_array(c-1)
!             cp_end_index   = c_p_index_array(c)
! 
!             p_c_obj_relation_array(cp_start_index:cp_end_index) = c
! 
!         enddo
!         
!         call quicksort(p_index , 1 , c_p_sum , p_c_obj_relation_array)
! 
!         p_c_index_array(0) = 0
!         p_c_index_array_duplicates = 0
!         total_c_count = 0
!         c_count=0
!         do p=1,npoin
!         
!             ! c_count the number of cells adjacent to each point
! 
!             do while (p .eq. p_index(c_count))
!                 c_count = c_count + 1
!             enddo
! 
!             ! flag duplicates
! 
!             call sort_and_flag_duplicates(p_c_obj_relation_array, p_c_index_array(p-1)+1, c_count)
! 
!             ! adjust c_count
! 
!             p_c_index_array(p) = c_count
!             
!             
!             p_c_index_array_duplicates = c_count - 1 ! starts looping at beginning of a point and counts up
!             do while (p_c_obj_relation_array(p_c_index_array(p)-p_c_index_array_duplicates(p)) .lt. 0)
!                 p_c_index_array_duplicates(p) = p_c_index_array_duplicates(p) - 1
!             enddo
! 
!             ! I'ma be real here
!             ! I'm pretty smug about this loop
!             ! update: less so now
! 
!         enddo
! 
!         p_c_index_array(:) = p_c_index_array(:) - p_c_index_array_duplicates(:)
! 
!         call remove_flagged_duplicates(p_c_obj_relation_array)
! 
!         p_c_sum = size(p_c_obj_relation_array)
        
        
    end subroutine obj_relation_inverter
    
    
    
    
    
    ! end contains    
    
end module preprocessor_routine_module
