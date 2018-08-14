 MODULE constants
!--------------------------------------------------------------!
! Numerical constants and constants for variable declarations. !
!--------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,PARAMETER :: dp=kind(1.d0)
 END MODULE constants


 MODULE utils
!--------------------------!
! Miscellaneous utilities. !
!--------------------------!
 USE constants
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



 SUBROUTINE gauss_reduce(a,b)
 IMPLICIT NONE
 REAL(dp) :: a(3), b(3), temp_vec(3), a_len, b_len, len
 REAL(dp),PARAMETER :: tol_zero=1.d-8
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
 REAL(dp),INTENT(inout) :: vecs(3,3)
 REAL(dp),PARAMETER :: tol_zero=1.d-8
 INTEGER :: longest, i, j, x1, x2 
 REAL(dp) :: temp_vec(3), best_vec(3), maxlen, minlen, nlen, len2(3), &
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
 REAL(dp),INTENT(in) :: A(3,3)
 REAL(dp) :: determinant33

 determinant33=A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))&
  &+A(1,2)*(A(3,1)*A(2,3)-A(2,1)*A(3,3))&
  &+A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))

 END FUNCTION determinant33


 SUBROUTINE supercells_generator(num_pcells,num_hnf,hnf)
!------------------------------------------------------------------------------!
! Generate all unique supercells that contain a given number of primitive unit !
! cells. See 'Hart and Forcade, Phys. Rev. B 77, 224115 (2008)' for details of !
! the algorithm.                                                               !
!------------------------------------------------------------------------------!
 IMPLICIT NONE
 INTEGER,INTENT(in) :: num_pcells
 INTEGER,INTENT(out) :: num_hnf
 INTEGER,POINTER :: hnf(:,:,:)
 INTEGER :: a,b,c,d,e,f,ialloc,count_hnf,quotient

 count_hnf=0

 do a=1,num_pcells 
  if(.not.mod(num_pcells,a)==0)cycle
  quotient=num_pcells/a
  do c=1,quotient  
   if(.not.mod(quotient,c)==0)cycle
   f=quotient/c
   count_hnf=count_hnf+c*f**2
  enddo ! c
 enddo ! a

 num_hnf=count_hnf
 count_hnf=0

 allocate(hnf(3,3,num_hnf),stat=ialloc)
 if(ialloc/=0)then
  write(*,*)'Problem allocating hnf array in supercells_generator.'
  stop
 endif

 hnf(1:3,1:3,1:num_hnf)=0

 do a=1,num_pcells 
  if(.not.mod(num_pcells,a)==0)cycle
  quotient=num_pcells/a
  do c=1,quotient  
   if(.not.mod(quotient,c)==0)cycle
   f=quotient/c
   do b=0,c-1
    do d=0,f-1
     do e=0,f-1
      count_hnf=count_hnf+1
      hnf(1,1,count_hnf)=a
      hnf(1,2,count_hnf)=b
      hnf(2,2,count_hnf)=c
      hnf(1,3,count_hnf)=d
      hnf(2,3,count_hnf)=e
      hnf(3,3,count_hnf)=f
     enddo ! e
    enddo ! d
   enddo ! b
  enddo ! c
 enddo ! a

 if(count_hnf/=num_hnf)then
  write(*,*)'Did not generate all HNF matrices.'
  stop
 endif 
 
 END SUBROUTINE supercells_generator




 SUBROUTINE inv33(v,inv)
!-----------------------!
! Inverts 3x3 matrices. !
!-----------------------!
 IMPLICIT NONE
 REAL(dp),INTENT(in) :: v(3,3)
 REAL(dp),INTENT(out) :: inv(3,3)
 REAL(dp) :: d

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


 PROGRAM optimal_supercell
!---------------------!
!  OPTIMAL_SUPERCELL  !
!---------------------!
 USE utils
 IMPLICIT NONE
 REAL(dp),PARAMETER :: tol=1.d-8
 
 INTEGER :: i,j,k,ierr,s11,s12,s13,s22,s23,s33, abs_norm, best_abs_norm, &
 &quotient,min_super_size,super_size,hnf(3,3),best_supercell(3,3) 
 REAL(dp) :: prim_latt_vecs(3,3), rec_latt_vecs(3,3), temp_latt_vecs(3,3), cross_vector(3)
 real(dp) :: radius, cell_volume, min_image_distance, best_min_image_distance, min_image_distance2
 LOGICAL :: found 
 
! Get the primitive cell lattice vectors and the target radius
 open(unit=11,file='cell_and_radius.dat',status='old',iostat=ierr)
 if(ierr/=0)then
  write(*,*)'Problem opening cell_and_radius.dat file.'
  stop
 endif ! ierr
 read(11,*,iostat=ierr)prim_latt_vecs(1,1:3)
 if(ierr==0)read(11,*,iostat=ierr)prim_latt_vecs(2,1:3)
 if(ierr==0)read(11,*,iostat=ierr)prim_latt_vecs(3,1:3)
 if(ierr==0)read(11,*,iostat=ierr)radius
 if(ierr/=0)then
  write(*,*)'Problem reading cell_and_radius.dat file.'
  stop
 endif ! ierr
 close(11)
 
 call inv33(prim_latt_vecs,rec_latt_vecs)
 rec_latt_vecs=transpose(rec_latt_vecs)

 cell_volume = determinant33(prim_latt_vecs)
 
 ! The minimal value for the supercell size is set by requiring that the volume
 ! of the supercell is at least equal to a cube that contains the target sphere
 min_super_size = int(radius**3/cell_volume)

 write(*,*) "Unit cell volume: ", cell_volume
 ! write(*,*) "Minimal cubic volume: ", (2*radius)**3
 write(*,*) "Minimal supercell size: ", min_super_size
 
 ! Initialize to zero the best minimum image distance
 best_min_image_distance = 0.0d0
 
 found = .false.
 
 do super_size = min_super_size,min_super_size + 2
   do s11=1,super_size
    if ( .not. mod(super_size,s11) == 0 ) cycle
    quotient=super_size/s11
    do s22=1,quotient
     if(.not.mod(quotient,s22)==0) cycle
     s33=quotient/s22
     ! IF (s33>1) CYCLE ! uncomment this line for 2D
     do s12=0,s22-1
      do s13=0,s33-1
       do s23=0,s33-1
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
        print*, 'HNF:'
        do k=1,3
            print*, hnf(k, :)
        end do
        print*, 'temp latt:'
        do k=1,3
            print*, temp_latt_vecs(k, :)
        end do
        
        call reduce_vecs(temp_latt_vecs) ! comment this line for 2D
        print*, 'reduced:'
        do k=1,3
            print*, temp_latt_vecs(k, :)
        end do
        
        !call gauss_reduce(temp_latt_vecs(1,1:3),temp_latt_vecs(2,1:3)) ! uncomment this line for 2D
        
        do k=1,3
         do j=1,3
          hnf(k,j)=nint(sum(temp_latt_vecs(k,1:3)*rec_latt_vecs(j,1:3)))  
         enddo ! j
        enddo ! k
        
        ! After the reduction, the minimum image distance is simply
        ! the length of the first lattice vector
        ! LK: No, I disagree there, this is only the minimum image distance for 
        ! an orthorhombic system. One has to take the angle into account!
        min_image_distance = sqrt(sum(temp_latt_vecs(1,1:3)*temp_latt_vecs(1,1:3)))
        print*, 'Min image distance:', min_image_distance
        DO k=1, 3
            cross_vector = cross(temp_latt_vecs(MOD(k,3)+1,1:3), temp_latt_vecs(MOD(k+1,3)+1,1:3))
            cross_vector(:) = cross_vector(:) / SQRT(SUM(cross_vector(:)*cross_vector(:)))
            min_image_distance2 = abs(sum(cross_vector(:)*temp_latt_vecs(k,1:3)))
            print*, 'Min image distance2:', min_image_distance2
        END DO
        
        if (min_image_distance > best_min_image_distance) then
           best_min_image_distance = min_image_distance
           best_supercell = hnf
           best_abs_norm = 0
           do j = 1, 3
              do k = 1, 3
                 best_abs_norm = best_abs_norm + abs(best_supercell(k,j))
              end do
           end do
           write(*,*) 'Minimum image distance:', min_image_distance
        elseif (abs(min_image_distance-best_min_image_distance)<tol) then
           abs_norm = 0 
           do j = 1, 3
              do k = 1, 3
                 abs_norm = abs_norm + abs(hnf(k,j))
              end do
           end do
           if (abs_norm < best_abs_norm) then
               best_supercell = hnf
               best_abs_norm = abs_norm
           end if 
        end if
        
       enddo ! s23
      enddo ! s13
     enddo ! s12
    enddo ! s22
   enddo ! s11
   
   if (best_min_image_distance .ge. radius) then
      found = .true.
      exit
   end if
 enddo ! super_size


 if(.not.found)then
  write(*,*)'Unable to find an optimal supercell.'
  stop
 endif 
 
 write(*,*) "Best minimum image distance: ", best_min_image_distance
 write(*,*) "Optimal supercell: "
 do i = 1, 3
    write(*,*) best_supercell(i,1:3)
 end do
 
 open(unit=12,file='supercell.dat',status='replace',iostat=ierr)
 if(ierr/=0)then
    write(*,*)'Problem opening supercell.dat file.'
    stop
 endif ! ierr
 write(12,*) best_supercell(1,1:3)
 write(12,*) best_supercell(2,1:3)
 write(12,*) best_supercell(3,1:3)
 close(12)

 END PROGRAM optimal_supercell
