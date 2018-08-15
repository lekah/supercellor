! MODULE constants
!--------------------------------------------------------------!
! Numerical constants and constants for variable declarations. !
!--------------------------------------------------------------!
! IMPLICIT NONE
! integer, parameter :: sp = selected_real_kind(6,37) ! single precision
! INTEGER,PARAMETER :: dp=selected_real_kind(15,307) ! double precision
!
! real(sp) :: r_sp = 1.0
! REAL(dp) :: r_dp = 1.0_dp
! For now commenting out constants module, this just doesn't properly work with f2py,
! have to declare REAL*8 everywhere
! END MODULE constants


 MODULE utils
!--------------------------!
! Miscellaneous utilities. !
!--------------------------!
 USE constants
 USE, intrinsic :: iso_fortran_env
 IMPLICIT NONE
 
 CONTAINS

 INTEGER FUNCTION gcd(int_1,int_2)
!----------------------------------------------------------------------!
! Calculate the greatest common divisor of two positive integers using !
! Euclid's algorithm.                                                  !
!----------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: int_1,int_2
 INTEGER :: a,b,temp

 a=int_1 ; b=int_2

 if(a<b)then
  temp=a
  a=b
  b=temp
 endif ! a<b

 do
  temp=mod(a,b)
  if(temp==0)exit
  a=b
  b=temp
 enddo
 gcd=b
 END FUNCTION gcd


FUNCTION cross(a, b)
  !
  ! Calculates the cross product of 2 3d-vectors
  !
  IMPLICIT NONE
  REAL*8, DIMENSION(3) :: cross
  REAL*8, DIMENSION(3), INTENT(IN) :: a, b
  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)
END FUNCTION cross

FUNCTION diag_vol(cell, radius)
  !
  ! Calculates the diagonal cell and gives the number of repetitions necessary
  !
  IMPLICIT NONE
  INTEGER  :: diag_vol
  REAL*8, DIMENSION(3,3), INTENT(IN) :: cell
  REAL*8, INTENT(IN) :: radius

  REAL*8, DIMENSION(3) :: cross_vector

  INTEGER :: i
  diag_vol = 1
  DO i=1, 3
    cross_vector = cross(cell(MOD(i,3)+1,1:3), cell(MOD(i+1,3)+1,1:3))
    cross_vector(:) = cross_vector(:) / SQRT(SUM(cross_vector(:)*cross_vector(:)))
    diag_vol = diag_vol * ceiling(radius / abs(sum(cross_vector(:)*cell(i,1:3))))
  ENDDO

END FUNCTION diag_vol




 SUBROUTINE gauss_reduce(a,b)
 IMPLICIT NONE
 REAL*8 :: a(3), b(3), temp_vec(3), a_len, b_len, len
 REAL*8,PARAMETER :: tol_zero=1.d-8
 INTEGER :: r
 
 iter: DO
 a_len = sum(a * a)
 b_len = sum(b * b)
 
 IF (a_len - b_len > tol_zero) THEN
    temp_vec = b
    len = b_len
    !
    b = a 
    b_len = a_len
    !
    a = temp_vec
    a_len = len 
 END IF
 r = nint(sum(b(1:3) * a(1:3)) /a_len)
 temp_vec = b - r * a
 len = sum(temp_vec * temp_vec)
 IF (a_len - len > tol_zero) THEN
    b = a
    a = temp_vec
    CYCLE
 END IF
 EXIT 
 END DO iter
  
 END SUBROUTINE gauss_reduce
 
 SUBROUTINE reduce_vecs(vecs)
 !-----------------------------------------------------------------------------!
 ! Given n vectors a(i) that form a basis for a lattice L in n dimensions, the ! 
 ! subroutine returns the shortest possible basis vectors for L, that is:      !
 !                                                                             !
 ! - a(1) is the shortest non-zero vector in L                                 !
 ! - for i>1, a(i) is the shortest possible vector in L such that a(i)>=a(i-1) !
 !  and the set of vectors a(1) to a(i) are linearly independent               !
 !-----------------------------------------------------------------------------!
 ! From http://www.csie.nuk.edu.tw/~cychen/Lattices/A%203-Dimensional%20Lattice%20Reduction%20Algorithm.pdf
  
 IMPLICIT NONE
 REAL*8,INTENT(inout) :: vecs(3,3)
 REAL*8,PARAMETER :: tol_zero=1.d-8
 INTEGER :: longest, i, j, x1, x2 
 REAL*8 :: temp_vec(3), best_vec(3), maxlen, minlen, nlen, len2(3), &
 &prod12,prod13,prod23, denom, y1, y2

iter: DO
! Determine which of the three input vectors is the longest
 maxlen=0
 DO i=1,3
  nlen=vecs(i,1)**2+vecs(i,2)**2+vecs(i,3)**2
! Test nlen>maxlen within some tolerance to avoid floating point problems
  IF(nlen-maxlen>tol_zero)THEN
   maxlen=nlen
   longest=i
  ENDIF
 ENDDO ! i
 
 ! Put the longest vector as the third one 
 temp_vec(1:3) = vecs(3,1:3)
 vecs(3,1:3) = vecs(longest,1:3)
 vecs(longest,1:3) = temp_vec(1:3)
 
 ! Gauss-reduce the first two vectors
 call gauss_reduce(vecs(1,1:3),vecs(2,1:3))
 
 DO i=1,3
   len2(i) = sum(vecs(i,1:3) * vecs(i,1:3))
 END DO 
 prod12 = sum(vecs(1,1:3) * vecs(2,1:3))
 prod13 = sum(vecs(1,1:3) * vecs(3,1:3))
 prod23 = sum(vecs(2,1:3) * vecs(3,1:3))
 
 denom = 1.0d0 - prod12**2 / len2(1) /len2(2)
 y1 = - (prod23 - prod12 * prod13/len2(1))/len2(2)/denom
 y2 = - (prod13 - prod12 * prod23/len2(2))/len2(1)/denom
 
 ! Initialize to a large number the minimum length
 minlen = 100 * maxlen
 DO i=-1,1
    x1 = nint(y1+i)
    DO j=-1,1
       x2 = nint(y2+j)
       temp_vec = vecs(3,1:3) + x2 * vecs(2,1:3) + x1 * vecs(1,1:3)
       nlen = sum(temp_vec * temp_vec)
       IF (minlen - nlen > tol_zero) THEN
          minlen = nlen
          best_vec = temp_vec
       END IF 
    END DO ! j
 END DO ! i
 
 IF (maxlen - minlen > tol_zero) THEN
    vecs(3,1:3) = best_vec(1:3)
    CYCLE
 END IF
 EXIT
 END DO iter

 END SUBROUTINE reduce_vecs


 FUNCTION determinant33(A)
!-----------------------------------------------------!
! Given a 3x3 matrix A, this function returns det(A). !
!-----------------------------------------------------!
 IMPLICIT NONE
 REAL*8,INTENT(in) :: A(3,3)
 REAL*8 :: determinant33

 determinant33=A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))&
  &+A(1,2)*(A(3,1)*A(2,3)-A(2,1)*A(3,3))&
  &+A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))

 END FUNCTION determinant33

 FUNCTION determinant33_int(A)
    IMPLICIT NONE
    INTEGER,  INTENT(in) :: A(3,3)
    INTEGER :: determinant33_int

    determinant33_int = A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3)) &
        +A(1,2)*(A(3,1)*A(2,3)-A(2,1)*A(3,3)) &
        +A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))

 END FUNCTION determinant33_int



 SUBROUTINE inv33(v,inv)
!-----------------------!
! Inverts 3x3 matrices. !
!-----------------------!
 IMPLICIT NONE
 REAL*8,INTENT(in) :: v(3,3)
 REAL*8,INTENT(out) :: inv(3,3)
 REAL*8 :: d
 d=determinant33(v)
 if(d==0.d0)then
  write(*,*)'Trying to invert a singular determinant.'
  stop
 endif ! d
 d=1.d0/d
 inv(1,1)=(v(2,2)*v(3,3)-v(2,3)*v(3,2))*d
 inv(1,2)=(v(3,2)*v(1,3)-v(1,2)*v(3,3))*d
 inv(1,3)=(v(1,2)*v(2,3)-v(1,3)*v(2,2))*d
 inv(2,1)=(v(3,1)*v(2,3)-v(2,1)*v(3,3))*d
 inv(2,2)=(v(1,1)*v(3,3)-v(3,1)*v(1,3))*d
 inv(2,3)=(v(2,1)*v(1,3)-v(1,1)*v(2,3))*d
 inv(3,1)=(v(2,1)*v(3,2)-v(2,2)*v(3,1))*d
 inv(3,2)=(v(3,1)*v(1,2)-v(1,1)*v(3,2))*d
 inv(3,3)=(v(1,1)*v(2,2)-v(1,2)*v(2,1))*d

 END SUBROUTINE inv33

 END MODULE utils


 SUBROUTINE optimal_supercell_hnf(prim_latt_vecs, radius, verbosity, scaling_matrix, supercell_latt_vecs)
!---------------------!
!  OPTIMAL_SUPERCELL  !
!---------------------!
 USE utils
 use, intrinsic :: iso_fortran_env
 IMPLICIT NONE
 REAL*8, INTENT(IN), DIMENSION(3,3) :: prim_latt_vecs
 
 INTEGER, INTENT(OUT), DIMENSION(3,3) :: scaling_matrix
 REAL*8, INTENT(OUT), DIMENSION(3,3) :: supercell_latt_vecs
 REAL*8, INTENT(IN) :: radius
 INTEGER, INTENT(IN) :: verbosity ! verbosity settings 0-3

 REAL*8,PARAMETER :: tol=1.d-8
 INTEGER :: i,j,k,s11,s12,s13,s22,s23,s33, & !abs_norm, best_abs_norm, &
 &quotient,min_super_size, max_super_size, super_size,hnf(3,3) !,best_supercell(3,3) 
 REAL*8 :: rec_latt_vecs(3,3), temp_latt_vecs(3,3), cross_vector(3)
 REAL*8 :: cell_volume, abs_best_min_image_distance, & !min_image_distance, &
        min_image_distance_this_hnf, minimum_image_distances(3)
 LOGICAL :: found 
 print*, radius
 print*, scaling_matrix
 print*, prim_latt_vecs
! Get the primitive cell lattice vectors and the target radius

 
 call inv33(prim_latt_vecs,rec_latt_vecs)
 rec_latt_vecs=transpose(rec_latt_vecs)

 cell_volume = determinant33(prim_latt_vecs)
 
 ! The minimal value for the supercell size is set by requiring that the volume
 ! of the supercell is at least equal to a cube that contains the target sphere
 min_super_size = int(radius**3/cell_volume)
 max_super_size = diag_vol(prim_latt_vecs, radius)
 IF ( verbosity > 0 ) THEN
    write(*,*) "Unit cell volume: ", cell_volume
    write(*,*) "Minimal supercell size: ", min_super_size
    write(*,*) "Maximal supercell size: ", max_super_size
 ENDIF

 ! Initialize to zero the best minimum image distance
 abs_best_min_image_distance = 0.0d0
 
 found = .false.
 
 volume_loop: do super_size = min_super_size, max_super_size
   s11_loop: do s11=1,super_size
    if ( .not. mod(super_size,s11) == 0 ) cycle
    quotient=super_size/s11
    s22_loop: do s22=1,quotient
     if(.not.mod(quotient,s22)==0) cycle
     s33=quotient/s22
     ! IF (s33>1) CYCLE ! uncomment this line for 2D
     s12_loop: do s12=0,s22-1
      s13_loop: do s13=0,s33-1
       s23_loop: do s23=0,s33-1
        ! construct the supercell matrix
        hnf(1:3,1:3)=0
        hnf(1,1)=s11
        hnf(1,2)=s12
        hnf(1,3)=s13
        hnf(2,2)=s22
        hnf(2,3)=s23
        hnf(3,3)=s33
        do k=1,3
          do j=1,3
           temp_latt_vecs(k,j)=sum(dble(hnf(k,1:3))*prim_latt_vecs(1:3,j))
          enddo ! j
        enddo ! k
        IF ( verbosity > 1 ) THEN
            print*, 'HNF:'
            do k=1,3
                print*, hnf(k, :)
            end do
            print*, 'temp latt:'
            do k=1,3
                print*, temp_latt_vecs(k, :)
            end do
        ENDIF
        
        call reduce_vecs(temp_latt_vecs) ! comment this line for 2D

        do k=1,3
         do j=1,3
          hnf(k,j)=nint(sum(temp_latt_vecs(k,1:3)*rec_latt_vecs(j,1:3)))  
         enddo ! j
        enddo ! k

        IF ( verbosity > 1 ) THEN
            print*, 'reduced:'
            do k=1,3
                print*, temp_latt_vecs(k, :)
            end do
            print*, 'Reduced scaling matrix:'
            do k=1,3
                print*, hnf(k, :)
            end do
        ENDIF
        
        ! After the reduction, the minimum image distance is simply
        ! the length of the first lattice vector
        ! LK: No, I disagree there, this is only the minimum image distance for 
        ! an orthorhombic system. One has to take the angle into account!
        ! min_image_distance = sqrt(sum(temp_latt_vecs(1,1:3)*temp_latt_vecs(1,1:3)))
        ! print*, 'Min image distance:', min_image_distance
        !min_image_distance_this_hnf = 0.0d0
        DO k=1, 3
            cross_vector = cross(temp_latt_vecs(MOD(k,3)+1,1:3), temp_latt_vecs(MOD(k+1,3)+1,1:3))
            cross_vector(:) = cross_vector(:) / SQRT(SUM(cross_vector(:)*cross_vector(:)))
            minimum_image_distances(k) = abs(sum(cross_vector(:)*temp_latt_vecs(k,1:3)))
            IF ( verbosity > 1) THEN
                print*, 'Min image distance2:', minimum_image_distances(k)
            ENDIF

            IF ( minimum_image_distances(k) < radius ) THEN
                IF ( verbosity > 1 ) print*, 'found too small mim image distance'
                cycle s23_loop
            ENDIF
        END DO
        ! If I got here, it means I did not break minimum image distance criterion
        ! and that I have found a valid supercell
        found = .true.
        min_image_distance_this_hnf = minval(minimum_image_distances(:))
        if (min_image_distance_this_hnf > abs_best_min_image_distance) then
           abs_best_min_image_distance = min_image_distance_this_hnf
           ! best_supercell = hnf
           scaling_matrix(1:3,1:3) = hnf(1:3,1:3)
           supercell_latt_vecs(1:3,1:3) = temp_latt_vecs(1:3,1:3)
           !best_abs_norm = 0
           !do j = 1, 3
           !   do k = 1, 3
           !      best_abs_norm = best_abs_norm + abs(best_supercell(k,j))
           !   end do
           !end do
           IF ( verbosity > 0 ) print*, 'New best minimum image distance:', abs_best_min_image_distance
!~         elseif (abs(min_image_distance-best_min_image_distance)<tol) then
!~            abs_norm = 0 
!~            do j = 1, 3
!~               do k = 1, 3
!~                  abs_norm = abs_norm + abs(hnf(k,j))
!~               end do
!~            end do
!~            if (abs_norm < best_abs_norm) then
!~                best_supercell = hnf
!~                best_abs_norm = abs_norm
!~            end if 
        end if
        
       enddo s23_loop
      enddo s13_loop
     enddo s12_loop
    enddo s22_loop
   enddo s11_loop
   if ( found ) THEN 
    if (verbosity > 0 ) print*, "Exiting, I found a good volume", super_size
    exit volume_loop
   endif
 enddo volume_loop

 if(.not.found)then
  write(*,*)'Unable to find an optimal supercell.'
  stop
 endif 
 if (verbosity > 0 ) THEN
    write(*,*) "Best minimum image distance: ", abs_best_min_image_distance
    write(*,*) "Optimal supercell: "

     do i = 1, 3
        write(*,*) scaling_matrix(i,1:3)
     end do
 endif
 END SUBROUTINE optimal_supercell_hnf


PROGRAM optimal_supercell
 USE utils
 USE, intrinsic :: iso_fortran_env
 implicit none
 INTEGER ierr
 REAL*8, dimension(3,3) :: cell
 INTEGER, dimension(3,3) ::  scaling_matrix
 REAL*8 :: radius
 INTEGER :: k
 open(unit=11, file='cell_and_radius.dat', status='old', iostat=ierr)
 if(ierr/=0)then
  write(*,*)'Problem opening cell_and_radius.dat file.'
  stop
 endif ! ierr
 read(11,*,iostat=ierr) cell(1,1:3)
 if(ierr==0) read(11,*,iostat=ierr) cell(2,1:3)
 if(ierr==0) read(11,*,iostat=ierr) cell(3,1:3)
 if(ierr==0) read(11,*,iostat=ierr) radius
 if(ierr/=0) then
  write(*,*)'Problem reading cell_and_radius.dat file.'
  stop
 endif ! ierr
 close(11)
 call optimal_supercell_hnf(cell, radius, scaling_matrix)
 open(unit=12,file='supercell.dat',status='replace',iostat=ierr)
 if(ierr/=0)then
    write(*,*)'Problem opening supercell.dat file.'
    stop
 endif ! ierr    

 do k=1, 3
    write(12,*) scaling_matrix(k,1:3)
    !print*, scaling_matrix(k)
 end do
 close(12)
END PROGRAM optimal_supercell
