module Simple_Sorting
   
   integer, private, parameter :: wp=KIND(0.0) ! Single or double precision?
   
contains

SUBROUTINE HeapSort ( array )
 !
 !  Purpose:
 !     Subroutine to sort a real array. This routine uses the
 !     heapsort technique.  It was based upon the examples
 !     found in NUMERICAL RECIPES, by Press, Flannery, Teukolsky,
 !     and Vetterling. 
 !
 IMPLICIT NONE

 ! Declare local parameters
 INTEGER, PARAMETER :: SGL = wp ! Precision

 ! Declare calling arguments
 REAL(KIND=SGL), DIMENSION(:), INTENT(INOUT) :: array
 ! Array to sort

 ! List of local variables:
 INTEGER :: n
 INTEGER :: i                ! Index of a heap node
 INTEGER :: u_heap           ! Highest heap index
 INTEGER :: j                ! Index of a child of node i
 INTEGER :: l_heap           ! Lowest heap index
 REAL(KIND=SGL) :: temp      ! Temp variable for value at node i
 
 n=SIZE(array)

 ! If l_heap<=i<=u_heap, i is a heap node. ! If 2*i<=u_heap, 2*i is 
 ! a child of i and if 2*i+1<=u_heap, 2*i+1 is a child of i.
 ! If i/=l_heap and j is a child of i, array(i)>=array(j).
 ! For example, l_heap = 1 and U_heap = 14, the heap has the form
 !                               1
 !               |---------------|---------------|
 !               2                               3
 !       |-------|-------|               |-------|-------|
 !       4               5               6               7
 !  |----|----|     |----|----|     |----|----|     |----|   
 !  8         9    10        11    12        13    14      

 ! Initialize the heap (no node has a child)
    l_heap  = n / 2 + 1
    u_heap = n

    mainloop: DO
    ! In the first n/2 iterations, the stack is extended by 
    ! decrementing l_heap.
    ! In the following n-1 iterations, the largest element is in array(1),
    ! which is swapped with array(u_heap), then u_heap is decremented. 
       IF ( l_heap > 1 ) THEN
          l_heap    = l_heap - 1
          temp = array(l_heap)
       ELSE
          temp      = array(u_heap)
          array(u_heap) = array(1)
          u_heap        = u_heap - 1
          IF ( u_heap == 1 ) THEN
             !
             !            All done.  Store final value.
             !
             array(1) = temp
             EXIT mainloop
             !
          END IF
       END IF
       ! temp now holds array(l_heap). 
       i = l_heap
       j = 2*i

       DO WHILE (j<=u_heap)
          IF ( j < u_heap ) THEN
            ! Make j be the child of i with the greater value
             IF ( array(j) < array(j+1) ) j = j + 1
          END IF
          IF ( temp < array(j) ) THEN
            ! Swap i with j so that the value at i is greater than
            ! the value at a child of i. temp now holds the value at j.
             array(i) = array(j)
            ! Move to the child j.
             i = j
             j = 2*i
          ELSE
            ! No need for more swaps. 
             j = u_heap + 1
          END IF
          ! Note that the number of iterations of this DO WHILE 
          ! loop cannot exceed log_2(n). 
       END DO
       array(i) = temp
    END DO mainloop

END SUBROUTINE

SUBROUTINE QuickSort(list)

! Quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! Modified by Alan Miller to include an associated integer array which gives
! the positions of the elements in the original order.

IMPLICIT NONE
REAL(wp), DIMENSION (:), INTENT(IN OUT)  :: list

! Local variable
!INTEGER, DIMENSION (SIZE(list))  :: order
INTEGER :: i

!DO i = 1, SIZE(list)
!  order(i) = i
!END DO

CALL quick_sort_1(1, SIZE(list))

CONTAINS

RECURSIVE SUBROUTINE quick_sort_1(left_end, right_end)

INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j, itemp
REAL(wp)                :: reference, temp
INTEGER, PARAMETER  :: max_simple_sort_size = 6

IF (right_end < left_end + max_simple_sort_size) THEN
  ! Use interchange sort for small lists
  CALL interchange_sort(left_end, right_end)

ELSE
  ! Use partition ("quick") sort
  reference = list((left_end + right_end)/2)
  i = left_end - 1; j = right_end + 1

  DO
    ! Scan list from left end until element >= reference is found
    DO
      i = i + 1
      IF (list(i) >= reference) EXIT
    END DO
    ! Scan list from right end until element <= reference is found
    DO
      j = j - 1
      IF (list(j) <= reference) EXIT
    END DO


    IF (i < j) THEN
      ! Swap two out-of-order elements
      temp = list(i); list(i) = list(j); list(j) = temp
      !itemp = order(i); order(i) = order(j); order(j) = itemp
    ELSE IF (i == j) THEN
      i = i + 1
      EXIT
    ELSE
      EXIT
    END IF
  END DO

  IF (left_end < j) CALL quick_sort_1(left_end, j)
  IF (i < right_end) CALL quick_sort_1(i, right_end)
END IF

END SUBROUTINE quick_sort_1


SUBROUTINE interchange_sort(left_end, right_end)

INTEGER, INTENT(IN) :: left_end, right_end

!     Local variables
INTEGER             :: i, j, itemp
REAL(wp)               :: temp

DO i = left_end, right_end - 1
  DO j = i+1, right_end
    IF (list(i) > list(j)) THEN
      temp = list(i); list(i) = list(j); list(j) = temp
      !itemp = order(i); order(i) = order(j); order(j) = itemp
    END IF
  END DO
END DO

END SUBROUTINE interchange_sort

END SUBROUTINE
end module
