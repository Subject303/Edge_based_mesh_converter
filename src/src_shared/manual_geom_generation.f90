module manual_geom_generation_module
    
    use iso_fortran_env
    
    use object_counts
    use raw_data
    
    implicit none
    
    contains
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine manual_geom_generation
        implicit none
        
        call t_birch_slender(1.0, 2.0, -3.0, 15.0, 0.0, 1, 2, 4, 2)
        
        stop
        
    end subroutine manual_geom_generation
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    subroutine t_birch_slender (inner_r, outer_r, start_length, end_length, nose, scale_var ,n_highways, n_rotational, n_highways_wake)
        implicit none
        integer(KIND=INT32)            :: i, j, n, n_lengthways, n_highways, n_rotational, n_highways_wake, wake_length, wake_start, wake_end, wake_i, scale_var, plane_npoin
        real(KIND=REAL32)              :: inner_r, outer_r, start_length, end_length, nose, single_coord(2), lambda, dx, dydx
        real(KIND=REAL32),allocatable  :: temp_coords(:,:,:), temp_wake_coords(:,:,:)
        real(KIND=REAL64)              :: pi, dOO
        
        n_lengthways = 1 + (end_length - start_length) * scale_var
        
        allocate(temp_coords(n_lengthways,2,n_highways))
        
        temp_coords = 0.0
        
        dx   = (end_length - start_length)/(n_lengthways-1)
        dydx = (outer_r - inner_r)        /(n_lengthways-1)
        
        do i=1,n_lengthways
            
            lambda = start_length + (i-1)*dx
            
            temp_coords(i,:,1) = paramaterised_line(lambda, 3, nose)
            
            temp_coords(i,1,:) = temp_coords(i,1,1)
            
            if ((lambda - nose) .gt. 10.) then
                wake_end = wake_end + 1
                temp_coords(i,2,n_highways) = temp_coords(wake_start,2,n_highways)
            else
                temp_coords(i,2,n_highways) = inner_r + (i-1)*dydx
                wake_start= i + 1
                wake_end  = i
            endif
            
        end do
        
        do i=1,n_lengthways
            
            dx = (temp_coords(i,2,n_highways) - temp_coords(i,2,1))/(n_highways-1)
                
            do n=1,n_highways
                
                temp_coords(i,2,n) = temp_coords(i,2,1) + dx*(n-1)
                
            enddo
            
        end do
        
        wake_length = 1 + wake_end - wake_start
        allocate(temp_wake_coords(wake_length,2,n_highways_wake))
        
        temp_wake_coords(:,2,1) = 0.0
        
        wake_i = 0
        
        do i=wake_start,wake_end
            
            wake_i = wake_i + 1
            
            dx = temp_coords(i,2,1)/(n_highways_wake-1)
            
            temp_wake_coords(wake_i,1,:) = temp_coords(i,1,1)
            
            do n=2,n_highways_wake
                
                temp_wake_coords(wake_i,2,n) = temp_wake_coords(wake_i,2,1) + dx*(n-1)
                
            enddo
            
        end do
        
        open(12,file="planar_mesh_coords.csv",access='sequential',action='write',status='replace')
        
        do i=1,n_lengthways
            
            do n=1,n_highways
                
                !print*, temp_coords(i,1,n) , ',', temp_coords(i,2,n)
                
                write(12,'(f,a,f)') temp_coords(i,1,n) , ',', temp_coords(i,2,n)
                
            enddo
            
        enddo        
        
        do i=1,wake_length
            
            do n=1,n_highways_wake
                
                !print*, temp_wake_coords(i,1,n) , ',', temp_wake_coords(i,2,n)
                
                write(12,'(f,a,f)') temp_wake_coords(i,1,n) , ',', temp_wake_coords(i,2,n)
                
            enddo
            
        enddo
        
        close(12)
        
        pi = 3.14159265359
        
        dOO = 2 * pi / n_rotational
        
        print*, pi, dOO, sind(dOO), cosd(dOO), tand(dOO)
        
        npoin = ((n_lengthways*n_highways) + (wake_length*n_highways_wake)) * n_rotational
        
        allocate(coords(npoin,3))
        
        coords = 0.0
        
        j=0
        do i=1,n_lengthways
            do n=1,n_highways
                j=j+1
                coords(j,1) = temp_coords(i,1,n)
                coords(j,2) = temp_coords(i,2,n)
            enddo
        enddo  
        do i=1,wake_length
            do n=1,n_highways_wake
                j=j+1
                coords(j,1) = temp_wake_coords(i,1,n)
                coords(j,2) = temp_wake_coords(i,2,n)
            enddo
        enddo  
        
        plane_npoin = j
        
        coords(1:plane_npoin,3) = 0.0
        
        do n=2, n_rotational
            
            i = (plane_npoin * (n-1)) + 1
            j = (plane_npoin *  n   )
            
            coords(i:j,1) =                coords(1:plane_npoin,1)
            coords(i:j,2) = cos(doo*(n-1))*coords(1:plane_npoin,2)! - sin(doo*n)*coords(1:plane_npoin,3)
            coords(i:j,3) = sin(doo*(n-1))*coords(1:plane_npoin,2)! + cos(doo*n)*coords(1:plane_npoin,3)
            
        enddo
        
        do i=1,npoin
            print*, coords(i,:)
        enddo
        
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
                
!                 if ((lambda - beta) .ge. 10.0) then
!                     
!                     paramaterised_line(1) = lambda
!                     paramaterised_line(2) = 0.0
!                     return
!                 endif
                
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
