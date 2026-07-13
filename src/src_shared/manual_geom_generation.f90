module manual_geom_generation_module
    
    use iso_fortran_env
    
    use object_counts
    use raw_data
    
    implicit none
    
    contains
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine manual_geom_generation
        implicit none
        
        call t_birch_slender(1.0, 3.0, -3.0, 15.0, 0.0, 100, 100, 100)
        
        stop
        
    end subroutine manual_geom_generation
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine t_birch_slender (inner_r, outer_r, start_length, end_length, nose, n_lengthways ,n_highways, n_rotational)
        implicit none
        integer(KIND=INT32)            :: i, n, n_lengthways, n_highways, n_rotational
        real(KIND=REAL32)              :: inner_r, outer_r, start_length, end_length, nose, single_coord(2), lambda, dx, dydx
        real(KIND=REAL32),allocatable  :: temp_coords(:,:,:)
        
        allocate(temp_coords(n_lengthways,2,n_highways))
        
        temp_coords = 0.0
        
        dx   = (end_length - start_length)/(n_lengthways-1)
        dydx = (outer_r - inner_r)        /(n_lengthways-1)
        
        do i=1,n_lengthways
            
            lambda = start_length + (i-1)*dx
            
            temp_coords(i,:,1) = paramaterised_line(lambda, 3, nose)
            
            temp_coords(i,1,:) = temp_coords(i,1,1)
            
            temp_coords(i,2,n_highways) = inner_r + (i-1)*dydx
            
        end do
        
!         do i=1,n_lengthways
!             
!             dx = (temp_coords(i,2,n_highways) - temp_coords(i,2,1))/n_highways
!                 
!             do n=2,(n_highways-1)
!                 
!                 temp_coords(i,2,n) = dx*n
!                 
!             enddo
!             
!         end do
        
        open(12,file="planar_mesh_coords.csv",access='sequential',action='write',status='replace')
        
        do i=1,n_lengthways
            
            do n=1,n_highways
                
                write(12,'(f,a,f)') temp_coords(i,1,n) , ',', temp_coords(i,2,n)
                
            enddo
            
        enddo
        
        close(12)
        
    end subroutine t_birch_slender
    
    function paramaterised_line(lambda, funct, beta)
        implicit none
        real(KIND=REAL32), dimension(2) :: paramaterised_line
        real(KIND=REAL32)             :: lambda, D
        real(KIND=REAL32),optional    :: beta
        integer(KIND=INT32)           :: funct
        integer(KIND=INT32),parameter :: straight_line=1, angled_line=2, trever_birch_proj=3
        
        select case (funct)
            case(straight_line)
            
                paramaterised_line(1) = lambda
                paramaterised_line(2) = beta
                
                return
                
            case(angled_line)
            
                paramaterised_line(1) = lambda
                paramaterised_line(2) = lambda / beta
                
                return
                
            case(trever_birch_proj)
                
                D  = 1.0
                
                if (lambda .le. beta) then
                
                    paramaterised_line(1) = lambda
                    paramaterised_line(2) = 0.0
                    return
                endif
                if ((lambda - beta) .ge. 10.0) then
                    
                    paramaterised_line(1) = lambda
                    paramaterised_line(2) = 0.0
                    return
                endif
                
                if ((lambda - beta) .gt. 3.0) then
                    
                    paramaterised_line(1) = lambda
                    paramaterised_line(2) = D / 2
                    return
                endif
                
                if (lambda .gt. beta) then
                    
                    !
                    ! r/D = -0.002615(x/D)^3 - 0.039867(x/D)^2 + 0.30984(x/D)
                    !
                    
                    D = lambda/D
                    
                    paramaterised_line(1) = lambda
                    paramaterised_line(2) = -0.002615*(D**3) - 0.039867*(D**2) + 0.30984*(D)
                    return
                    
                endif
                
        end select
        
        
        
        print*,  'missing function definition in mesh generator'
    end function
    
    !end contains
    
    
end module manual_geom_generation_module
