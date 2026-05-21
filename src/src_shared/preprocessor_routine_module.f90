


module preprocessor_routine_module
    use iso_fortran_env
    use raw_data
    use object_relation_data
    
    implicit none
    
    contains
    
    subroutine c_e_preprocess
        integer(KIND=INT32) :: c, cf_index, cf_start_index, cf_end_index, f, fe_index, fe_start_index, fe_end_index, ce_index, n_con_edges
        integer(KIND=INT32) :: c_e_sum
        
        c_start_index = 1
        c_e_sum = 0
        
        allocate(c_e_index_array(0:nele))
        c_e_index_array(0) = 0
        
        do c=1,nele
        
            cf_start_index = 1 + c_f_index_array(c-1)
            cf_end_index   = c_f_index_array(c)
            
            n_con_edges = 0
            do cf_index = cf_start_index, cf_end_index
                f = c_f_obj_relation_array(cf_index)
                
                n_con_edges =  n_con_edges + f_e_index_array(f) - f_e_index_array(f-1) ! all index arrays start at 0
            enddo
            
            c_e_index_array(c) = n_con_edges
        enddo
        
        c_e_sum = sum(c_e_index_array)
        
        do c=1,nele
            c_e_index_array(c) = c_e_index_array(c-1) + c_e_index_array(c)
        enddo
        
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
        
        do c=1,nele
            ce_start_index = 1 + c_e_index_array(c-1)
            ce_end_index   = c_e_index_array(c)
            
            call sort_and_flag_duplicates(c_e_obj_relation_array, ce_start_index, ce_end_index)
            
        enddo
        
        call remove_flagged_duplicates(c_e_obj_relation_array)
        
    end subroutine c_e_preprocess
    
    
    subroutine p_c_preprocess
        integer(KIND=INT32),allocatable :: point_index, p_c_index_array_duplicates
        integer(KIND=INT32) :: c, p, count, total_count, index

        allocate(p_c_obj_relation_array(c_p_sum), p_c_index_array(0:npoin), p_c_index_array_duplicates(0:npoin) )
        allocate(point_index,source=c_p_obj_relation_array)

        do c=1,nele
            cp_start_index = 1 + c_p_index_array(c-1)
            cp_end_index   = c_p_index_array(c)

            p_c_obj_relation_array(cp_start_index:cp_end_index) = c

        enddo
        
        call quicksort(point_index , 1 , c_p_sum , p_c_obj_relation_array)

        p_c_index_array(0) = 0
        p_c_index_array_duplicates = 0
        total_count = 0
        count=0
        index=1
        do p=1,npoin
        
            ! count the number of cells adjacent to each point

            do while (p .eq. point_index(count))
                count = count + 1
            enddo

            ! flag duplicates

            call sort_and_flag_duplicates(p_c_obj_relation_array, p_c_index_array(p-1)+1, count)

            ! adjust count

            p_c_index_array(p) = count

            do while (p_c_obj_relation_array(p_c_index_array(p)-p_c_index_array_duplicates(p)) .lt. 0)
                p_c_index_array_duplicates(p) = p_c_index_array_duplicates(p) + 1
            enddo

            ! I'ma be real here
            ! I'm pretty smug about this loop
            ! update: less so now

        enddo

        p_c_index_array(:) = p_c_index_array(:) - p_c_index_array_duplicates(:)

        call remove_flagged_duplicates(p_c_obj_relation_array)

        p_c_sum = size(p_c_obj_relation_array)

    end subroutine p_c_preprocess
    
    
    
    
    subroutine p_f_preprocess
    
    end subroutine p_f_preprocess
    
    
    
    
    
    subroutine p_e_preprocess
    
    end subroutine p_e_preprocess
    
    
    ! end contains    
    
end module preprocessor_routine_module
