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
 ! USE constants
 ! USE, intrinsic :: iso_fortran_env
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


 FUNCTION cross_normed(a, b)
 
  !
  ! Calculates the cross product of 2 3d-vectors
  ! and devides by the norm
  !
  IMPLICIT NONE
  REAL*8, DIMENSION(3) :: cross_normed
  REAL*8 :: norm
  REAL*8, DIMENSION(3), INTENT(IN) :: a, b
  cross_normed = cross(a, b)
  norm = sqrt(sum(cross_normed(1:3)*cross_normed(1:3)))
  cross_normed(:) = cross_normed(:) / norm


END FUNCTION cross_normed


 FUNCTION dot(a,b)
  ! Calculates the dot-product of 2 3d-vectors
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: a(3), b(3)
  REAL*8 :: dot
  dot = SUM(a(1:3) * b(1:3))
 END FUNCTION dot

 FUNCTION angle(a,b)
  ! Calculates the dot-product of 2 3d-vectors and divides by the norm
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: a(3), b(3)
  REAL*8 :: angle
  angle = SUM(a(1:3) * b(1:3)) / sqrt(sum(a(1:3) * a(1:3))) / sqrt(sum(b(1:3) * b(1:3)))
 END FUNCTION angle


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
 y1 = - (prod13 - prod12 * prod23/len2(2))/len2(1)/denom
 y2 = - (prod23 - prod12 * prod13/len2(1))/len2(2)/denom
 
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


 SUBROUTINE fort_optimal_supercell_bec(norms_of_sorted_Gr_r2, sorted_Gc_r2, &
    sorted_Gr_r2, r_outer_, v_diag, r_inner, verbosity, N, R_best, C_best)
  USE utils
  use, intrinsic :: iso_fortran_env
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: N
  REAL*8, INTENT(IN)  :: norms_of_sorted_Gr_r2(N)
  INTEGER, INTENT(IN) :: sorted_Gc_r2(N,3)
  REAL*8, INTENT(IN)  :: sorted_Gr_r2(N,3)
  REAL*8, INTENT(IN)  :: r_inner
  REAL*8,  INTENT(IN) :: r_outer_
  INTEGER, INTENT(IN) :: v_diag
  INTEGER, INTENT(IN) :: verbosity


  ! variable used here
  INTEGER  :: volume, min_volume=0


  INTEGER :: i1, i2, i3 ! counting vars
  REAL*8  :: norm1, norm2, norm3 ! Here I load the norms, for convenience
  REAL*8  :: d(3) ! for interplane-distancess
  REAL*8  :: r_outer
  REAL*8  :: vector1(3), vector2(3), vector3(3) ! Here I load vectors, convencience
  REAL*8  :: cross13(3), cross23(3), cross12(3)
  REAL*8  :: min_inter_face_dist, max_min_inter_face_dist=0.0D0



  REAL*8, PARAMETER :: EPSILON=1.0D-6


  INTEGER :: R(3,3) ! the R matrix
  INTEGER,  INTENT(OUT) :: R_best(3,3) ! the R matrix
  REAL*8, INTENT(OUT):: C_best(3,3) ! the C matrix = R*cell


  R_best(:,:) = 0
  C_best(:,:) = 0.0D0
  min_volume  = huge(min_volume)
  r_outer = r_outer_


  DO i1=1,N
    
    norm1 = norms_of_sorted_Gr_r2(i1)
    IF ( norm1 > (r_outer - EPSILON) ) EXIT
    vector1(1:3) = sorted_Gr_r2(i1,1:3)
    if (verbosity > 1) print*, '  setting vector1', vector1
    DO i2=i1+1, N
      norm2 = norms_of_sorted_Gr_r2(i2)
      IF ( norm2 > (r_outer - EPSILON) ) EXIT
      vector2(1:3) = sorted_Gr_r2(i2,1:3)
      if (verbosity > 1) print*, '    setting vector2', vector2
      IF (ABS(dot(vector1, vector2) / (norm1 * norm2)) > 0.5D0 - EPSILON) THEN
        IF (verbosity > 1) print*, '   -> Angle 12 < 60, continue'
        CYCLE
      ENDIF

      DO  i3=i2+1,N
        norm3 = norms_of_sorted_Gr_r2(i3)
        IF ( norm3 > (r_outer - EPSILON) ) EXIT
        vector3(1:3) = sorted_Gr_r2(i3,1:3)
        if (verbosity > 1) print*, '      setting vector3', vector3
        IF (ABS(dot(vector1, vector3) / (norm1 * norm3)) > 0.5D0 - EPSILON) THEN
          IF (verbosity > 1) print*, '   -> Angle 13 < 60, continue'
          CYCLE
        ENDIF
        IF (ABS(dot(vector2, vector3) / (norm2 * norm3)) > 0.5D0 - EPSILON) THEN
          IF (verbosity > 1) print*, '   -> Angle 23 < 60, continue'
          CYCLE
        ENDIF
        ! checking intersections of each plane
        cross23 = cross_normed(vector2, vector3)
        d(1) = ABS(dot(cross23, vector1))
        if ( d(1) < r_inner + EPSILON) THEN
          if (verbosity > 1) print*, '     -> d1', r_inner,' < r_inner, continue'
          cycle
        endif
        cross13 = cross_normed(vector1, vector3)
        d(2) = ABS(dot(cross13, vector2))
        if ( d(2) < r_inner + EPSILON) THEN
          if (verbosity > 1) print*, '     -> d2', r_inner,' < r_inner, continue'
          cycle
        endif
        cross12 = cross_normed(vector1, vector2)
        d(3) = ABS(dot(cross12, vector3))
        if ( d(3) < r_inner + EPSILON) THEN
          if (verbosity > 1) print*, '     -> d3', r_inner,' < r_inner, continue'
          cycle
        endif
        ! Now I load into R
        R(1, :) = sorted_Gc_r2(i1, :)
        R(2, :) = sorted_Gc_r2(i2, :)
        R(3, :) = sorted_Gc_r2(i3, :)
        volume = abs(determinant33_int(R))
        min_inter_face_dist = minval(d)
        if ((volume < min_volume ) .or. ( (volume .eq. min_volume) .and. &
                        (min_inter_face_dist > max_min_inter_face_dist) )) THEN
            min_volume = volume
            r_outer = DBLE(min_volume) / r_inner**2
            max_min_inter_face_dist = min_inter_face_dist

            if ( verbosity > 0) THEN
              print*, "New optimal supercell",volume, max_min_inter_face_dist
            endif
            R_best(:,:) = R(:,:)
            C_best(1,:) = vector1(:)
            C_best(2,:) = vector2(:)
            C_best(3,:) = vector3(:)
        endif
      ENDDO
    ENDDO
  ENDDO
 END SUBROUTINE fort_optimal_supercell_bec

 SUBROUTINE fort_optimal_supercell_hnf(prim_latt_vecs, radius, verbosity, scaling_matrix, supercell_latt_vecs)
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
 INTEGER :: counter, i,j,k,s11,s12,s13,s22,s23,s33, & !abs_norm, best_abs_norm, &
 &quotient,min_super_size, max_super_size, super_size,hnf(3,3) !,best_supercell(3,3) 
 REAL*8 :: rec_latt_vecs(3,3), temp_latt_vecs(3,3)
 REAL*8 :: cell_volume, abs_best_min_image_distance, & !min_image_distance, &
        min_image_distance, best_min_image_distance, best_abs_norm, abs_norm
        
 LOGICAL :: found 

 call inv33(prim_latt_vecs,rec_latt_vecs)
 rec_latt_vecs=transpose(rec_latt_vecs)

 cell_volume = determinant33(prim_latt_vecs)
 
 ! The minimal value for the supercell size is set by requiring that the volume
 ! of the supercell is at least equal to a cube that contains the target sphere
 min_super_size = 1 ! int(radius**3/cell_volume)
 max_super_size = diag_vol(prim_latt_vecs, radius)
 IF ( verbosity > 0 ) THEN
    write(*,*) "Unit cell volume: ", cell_volume
    write(*,*) "Minimal supercell size: ", min_super_size
    write(*,*) "Maximal supercell size: ", max_super_size
    print*, "The primitive cell:"
    do i=1,3
      write(*,*) prim_latt_vecs(i, 1:3)
    enddo
 ENDIF

 ! Initialize to zero the best_min_image_distance to required one
 found = .false.
 best_min_image_distance = radius
 counter = 0
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
        counter = counter + 1
        do k=1,3
          do j=1,3
           temp_latt_vecs(k,j)=sum(dble(hnf(k,1:3))*prim_latt_vecs(1:3,j))
          enddo ! j
        enddo ! k
        IF ( verbosity > 1 ) THEN
            print*, '-------------------------------------------'
            print*, counter
            print*, 'HNF:'
            do k=1,3
                print*, hnf(k, :)
            end do
            print*, 'temp latt:'
            do k=1,3
                print*, temp_latt_vecs(k, :)
            end do
        ENDIF
        !
        call reduce_vecs(temp_latt_vecs) ! comment this line for 2D

        do k=1,3
         do j=1,3
          hnf(k,j)=nint(sum(temp_latt_vecs(k,1:3)*rec_latt_vecs(j,1:3)))  
         enddo ! j
        enddo ! k

        IF ( verbosity > 1 ) THEN
            print*, 'Reduced scaling matrix:'
            do k=1,3
                print*, hnf(k, :)
            end do
            print*, 'reduced:'
            do k=1,3
                print*, temp_latt_vecs(k, :)
            end do
            print*, 'shortest'
            do k=1,3
                print*, sqrt(sum(temp_latt_vecs(k,1:3)*temp_latt_vecs(k,1:3)))
            enddo

        ENDIF
        ! After the reduction, the minimum image distance is simply
        ! the length of the first lattice vector
        min_image_distance = sqrt(sum(temp_latt_vecs(1,1:3)*temp_latt_vecs(1,1:3)))
        if (min_image_distance > best_min_image_distance) then
            ! This first case, I have a min_image_distance that's better than the
            ! set one. I found a valid solution
            found = .true.
            scaling_matrix(1:3,1:3) = hnf(1:3,1:3)
            supercell_latt_vecs(1:3,1:3) = temp_latt_vecs(1:3,1:3)

            ! However, it might be that there's an even better solution for this volume
            ! So, setting best_min_image_distance to current min_image_distance
            best_min_image_distance = min_image_distance

            ! I also calculate the best_abs_norm, in case of equal min_image_distance
            ! I will decide based on the norm.
            best_abs_norm = 0
            do j = 1, 3
                do k = 1, 3
                    best_abs_norm = best_abs_norm + abs(scaling_matrix(k,j))
                end do
            end do
        elseif (abs(min_image_distance-best_min_image_distance)<tol) then
            abs_norm = 0 
            do j = 1, 3
                do k = 1, 3
                    abs_norm = abs_norm + abs(hnf(k,j))
                end do
            end do
            if (abs_norm < best_abs_norm) then
                scaling_matrix(1:3,1:3) = hnf(1:3,1:3)
                best_abs_norm = abs_norm
            end if
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
 END SUBROUTINE fort_optimal_supercell_hnf


PROGRAM optimal_supercell
 USE utils
 USE, intrinsic :: iso_fortran_env
 implicit none
 INTEGER ierr
 REAL*8, dimension(3,3) :: cell, supercell
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
 call fort_optimal_supercell_hnf(cell, radius, 1, scaling_matrix, supercell)
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
