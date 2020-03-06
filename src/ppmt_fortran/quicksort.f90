! Module for sorting various arrays
! Author: John Manchuk, 2010
! QUICKSORT modified from:
!   Algorithms and Data Structures in F and Fortran
!   Robin A. Vowels
!   Unicomp, 1998
!
! An interface is provided so that several types of data can be
! sorted by calling the same function name.  Options include:
! 1. for A, B, ... real or real*8 or integer, with [] optional inputs
!
!   call qsortem(A,left,right,[B,C,...H])
!
! 2. for A a multidimensional real array and B a one dimensional
!    real array.  If B is present, all columns of A are sorted
!    based on the values in B, otherwise all columns of A are sorted
!    based on the specified column of A.
!
!   call qsortem(left,right,A,column,[B])
!
module quicksort
    implicit  none
    
    !Sort routines
    PUBLIC  :: qsortem, &   !Sort a set of variables based on another
               quick_sort, &   !Sort an individual variable and return a sorted index
               quick_sort_int

    interface QSORTEM
        module procedure SORTEM_REAL, SORTEM_REAL8, SORTEM_INT
        module procedure ARRAY_SORTEM_REAL, ARRAY_SORTEM_INT
        module procedure ARRAY_SORTEM_INT_2, ARRAY_SORTEM_REAL_2
    end interface QSORTEM
    
    private
    
    integer :: ALLOC_TEST
    
    !These pointers are used for a recursive array sort
    integer, pointer :: pA(:,:) => null()
    integer, pointer :: pT(:) => null()
    integer, pointer :: pI(:) => null()
    real*8, pointer :: pAr(:,:) => null()
    real*8, pointer :: pTr(:) => null()
    
contains

    !Sort a set of real arrays in terms of A.
    subroutine SORTEM_REAL8 ( A, left, right, B, C, D, E, F, G, H)
        real*8, dimension ( : ), intent (in out) :: A
        real*8, dimension ( : ), optional        :: B,C,D,E,F,G,H
        real*8,    allocatable, dimension ( : )  :: TEMP
        integer, allocatable, dimension ( : )  :: IDS
        integer, intent (in) :: left, right
        integer :: i, n
        
        !Setup indexing arrays
        n = right - left + 1
        allocate(IDS(left:right),TEMP(left:right), stat = ALLOC_TEST)
        if(ALLOC_TEST /= 0)then
            write(*,'(/,a,/)') 'ERROR: ALLOCATION FAILED IN SORTEM_REAL8'
            stop
        endif
        
        do i = left, right
            IDS(i) = i
        enddo
        
        !Sort A and indices
        call QUICK_SORT (A(left:right), IDS, 1, n)
        
        !If any of the other arrays are present, sort them by A
        if(present(B))then
            TEMP = B(left:right) ; B(left:right) = TEMP(IDS(:))
            if(present(C))then
                TEMP = C(left:right) ; C(left:right) = TEMP(IDS(:))
                if(present(D))then
                    TEMP = D(left:right) ; D(left:right) = TEMP(IDS(:))
                    if(present(E))then
                        TEMP = E(left:right) ; E(left:right) = TEMP(IDS(:))
                        if(present(F))then
                            TEMP = F(left:right) ; F(left:right) = TEMP(IDS(:))
                            if(present(G))then
                                TEMP = G(left:right) ; G(left:right) = TEMP(IDS(:))
                                if(present(H))then
                                    TEMP = H(left:right) ; H(left:right) = TEMP(IDS(:))
                                endif
                            endif
                        endif
                    endif
                endif
            endif
        endif
        deallocate(IDS,TEMP)
    end subroutine SORTEM_REAL8
    subroutine SORTEM_REAL ( A, left, right, B, C, D, E, F, G, H)
        real, dimension ( : ), intent (in out) :: A
        real, dimension ( : ), optional        :: B,C,D,E,F,G,H
        real,    allocatable, dimension ( : )  :: TEMP
        real*8,     allocatable, dimension ( : )  :: A_
        integer, allocatable, dimension ( : )  :: IDS
        integer, intent (in) :: left, right
        integer :: i, n
               
        !Setup indexing arrays
        n = right - left + 1
        allocate(A_(left:right),IDS(left:right),TEMP(left:right), stat = ALLOC_TEST)
        if(ALLOC_TEST /= 0)then
            write(*,'(/,a,/)') 'ERROR: ALLOCATION FAILED IN SORTEM_REAL'
            stop
        endif
        
        do i = left, right
            IDS(i) = i
            A_(i) = real(A(i),8)
        enddo
        
        !Sort A and indices
        call QUICK_SORT (A_(left:right), IDS, 1, n)
        A(left:right) = real(A_(left:right))
        
        !If any of the other arrays are present, sort them by A
        if(present(B))then
            TEMP = B(left:right) ; B(left:right) = TEMP(IDS(:))
            if(present(C))then
                TEMP = C(left:right) ; C(left:right) = TEMP(IDS(:))
                if(present(D))then
                    TEMP = D(left:right) ; D(left:right) = TEMP(IDS(:))
                    if(present(E))then
                        TEMP = E(left:right) ; E(left:right) = TEMP(IDS(:))
                        if(present(F))then
                            TEMP = F(left:right) ; F(left:right) = TEMP(IDS(:))
                            if(present(G))then
                                TEMP = G(left:right) ; G(left:right) = TEMP(IDS(:))
                                if(present(H))then
                                    TEMP = H(left:right) ; H(left:right) = TEMP(IDS(:))
                                endif
                            endif
                        endif
                    endif
                endif
            endif
        endif
        deallocate(IDS,TEMP,A_)
    end subroutine SORTEM_REAL
    subroutine SORTEM_INT ( A, left, right, B, C, D, E, F, G, H) !integer version
        integer, dimension ( : ), intent (in out) :: A
        integer, dimension ( : ), optional        :: B,C,D,E,F,G,H
        integer,    allocatable, dimension ( : )  :: TEMP
        real*8,     allocatable, dimension ( : )  :: A_
        integer, allocatable, dimension ( : )  :: IDS
        integer, intent (in) :: left, right
        integer :: i, n
               
        !Setup indexing arrays
        n = right - left + 1
        allocate(A_(left:right),IDS(left:right),TEMP(left:right), stat = ALLOC_TEST)
        if(ALLOC_TEST /= 0)then
            write(*,'(/,a,/)') 'ERROR: ALLOCATION FAILED IN SORTEM_INT'
            stop
        endif
        
        do i = left, right
            IDS(i) = i
        enddo
        
        do i = left, right
            A_(i) = real(A(i),8)
        enddo
        
        !Sort A and indices
        call QUICK_SORT (A_(left:right), IDS, 1, n)
        do i = left, right
            A(i) = int(A_(i))
        enddo
        
        !If any of the other arrays are present, sort them by A
        if(present(B))then
            TEMP = B(left:right) ; B(left:right) = TEMP(IDS(:))
            if(present(C))then
                TEMP = C(left:right) ; C(left:right) = TEMP(IDS(:))
                if(present(D))then
                    TEMP = D(left:right) ; D(left:right) = TEMP(IDS(:))
                    if(present(E))then
                        TEMP = E(left:right) ; E(left:right) = TEMP(IDS(:))
                        if(present(F))then
                            TEMP = F(left:right) ; F(left:right) = TEMP(IDS(:))
                            if(present(G))then
                                TEMP = G(left:right) ; G(left:right) = TEMP(IDS(:))
                                if(present(H))then
                                    TEMP = H(left:right) ; H(left:right) = TEMP(IDS(:))
                                endif
                            endif
                        endif
                    endif
                endif
            endif
        endif
        deallocate(IDS,TEMP,A_)
    end subroutine SORTEM_INT
    
    !Similar to above, but for sorting multidimensional arrays
    !Sort matrix B based on column col of B, or by values in array A
    subroutine ARRAY_SORTEM_REAL ( left, right, B, col, A)
        real*8, dimension (:,:), intent (in out) :: B
        integer, intent (in) :: col
        real*8, dimension ( : ), optional        :: A
        real*8,    allocatable, dimension ( : )  :: TEMP
        integer, allocatable, dimension ( : )  :: IDS
        integer, intent (in) :: left, right
        integer :: i, n
        logical :: rowmajor
        
        !Setup indexing arrays
        n = right - left + 1
        allocate(IDS(left:right),TEMP(left:right), stat = ALLOC_TEST)
        if(ALLOC_TEST /= 0)then
            write(*,'(/,a,/)') 'ERROR: ALLOCATION FAILED IN ARRAY_SORTEM_REAL'
            stop
        endif
        do i = left, right
            IDS(i) = i
        enddo
        
        !Determine array order
        rowmajor = .true.
        if(ubound(B,1) < ubound(B,2))then
            rowmajor = .false.
        endif
        
        !Determine sorted order of elements (in A or B)
        if(present(A))then
            TEMP = A(left:right)
        else
            if(rowmajor)then
                TEMP = B(left:right,col)
            else
                TEMP = B(col,left:right)
            endif
        endif
        
        call QUICK_SORT (TEMP, IDS, 1, n)
        
        !Re-order B and possibly A
        if(present(A))then
            A(left:right) = TEMP
        endif
        
        !Determine array order
        if(.not.rowmajor)then
            do i = 1, ubound(B,1)
                TEMP = B(i,left:right)
                B(i,left:right) = TEMP(IDS(:))
            enddo
        else
            do i = 1, ubound(B,2)
                TEMP = B(left:right,i)
                B(left:right,i) = TEMP(IDS(:))
            enddo
        endif

        deallocate(IDS,TEMP)
    end subroutine ARRAY_SORTEM_REAL
    subroutine ARRAY_SORTEM_INT ( left, right, B, col, A)
        integer, dimension (:,:), intent (in out) :: B
        integer, intent (in) :: col
        integer, dimension ( : ), optional :: A
        real*8, allocatable, dimension ( : ) :: TEMP
        integer, allocatable, dimension ( : ) :: ITEMP
        integer, allocatable, dimension ( : ) :: IDS
        integer, intent (in) :: left, right
        integer :: i, n
        logical :: rowmajor
        
        !Setup indexing arrays
        n = right - left + 1
        allocate(IDS(left:right),TEMP(left:right),ITEMP(left:right), stat = ALLOC_TEST)
        if(ALLOC_TEST /= 0)then
            write(*,'(/,a,/)') 'ERROR: ALLOCATION FAILED IN ARRAY_SORTEM_INT'
            stop
        endif
        do i = left, right
            IDS(i) = i
        enddo
        
        !Determine array order
        rowmajor = .true.
        if(ubound(B,1) < ubound(B,2))then
            rowmajor = .false.
        endif
        
        !Determine sorted order of elements (in A or B)
        if(present(A))then
            TEMP = real(A(left:right),8)
        else
            if(rowmajor)then
                TEMP = real(B(left:right,col),8)
            else
                TEMP = real(B(col,left:right),8)
            endif
        endif
        
        call QUICK_SORT (TEMP, IDS, 1, n)
        
        !Re-order B and possibly A
        if(present(A))then
            A(left:right) = floor(TEMP+0.1)
        endif
        
        !Determine array order
        if(.not.rowmajor)then
            do i = 1, ubound(B,1)
                ITEMP = B(i,left:right)
                B(i,left:right) = ITEMP(IDS(:))
            enddo
        else
            do i = 1, ubound(B,2)
                ITEMP = B(left:right,i)
                B(left:right,i) = ITEMP(IDS(:))
            enddo
        endif

        deallocate(IDS,TEMP,ITEMP)
    end subroutine ARRAY_SORTEM_INT
    
    subroutine ARRAY_SORTEM_INT_2 ( left, right, B, col)
        integer, target, dimension (:,:), intent (in out) :: B
        integer, intent (in) :: col(:)
        integer, allocatable, target, dimension ( : )  :: TEMP
        integer, allocatable, target, dimension ( : )  :: IDS
        integer, intent (in) :: left, right
        integer :: i, n, ncol, icol, t, j, i0, i1
        logical :: is_sorted(ubound(col,1))
        
        !Setup indexing arrays
        ncol = ubound(col,1)
        n = right - left + 1
        allocate(IDS(left:right),TEMP(left:right), stat = ALLOC_TEST)
        if(ALLOC_TEST /= 0)then
            write(*,'(/,a,/)') 'ERROR: ALLOCATION FAILED IN ARRAY_SORTEM_INT_2'
            stop
        endif
        do i = left, right
            IDS(i) = i
        enddo
        
        pA => B
        pT => TEMP
        pI => IDS
        
        call array_sortem_int_2_r ( left, right, col, 1 )
        
        deallocate(IDS,TEMP)
    end subroutine ARRAY_SORTEM_INT_2
    !Recursive sort for array_sortem_int_2, which sets global pointers pA,pT,pI
    recursive subroutine array_sortem_int_2_r ( left, right, cols, icol )
        integer, intent (in) :: left, right !left and right bounds
        integer, intent (in) :: cols(:) !columns being sorted
        integer, intent (in) :: icol !current column index
        integer :: i, j, i0, i1
        
        !Do nothing for 1 element
        if(right == left) return
        
        !Sort the current column and reorder the others
        pT(left:right) = pA(cols(icol),left:right)
        call quick_sort_int (pT, pI, left, right)
        RLOOP : do i = 1, ubound(pA,1)
            !do not re-order a column that is already sorted
            do j = 1, icol-1
                if(i == cols(j)) cycle RLOOP
            enddo
            pT(left:right) = pA(i,left:right)
            pA(i,left:right) = pT(pI(left:right))
        enddo RLOOP
        
        !Compute bounds based on the most recently sorted column
        i0 = left
        j = cols(icol)
        do while (i0 < right)
            i1 = i0
            pI(i0) = i0
            do while ( pA(j,i1) == pA(j,i0) )
                i1 = i1 + 1
                if(i1 > right) exit
                pI(i1) = i1
            enddo
            i1 = i1 - 1
            
            !Determine sorted order of elements
            if(icol < ubound(cols,1))then
                call array_sortem_int_2_r ( i0, i1, cols, icol + 1 )
            endif
            
            !Next subset
            i0 = i1 + 1
        enddo

        return
    end subroutine array_sortem_int_2_r
    
    
    
    
    
    subroutine ARRAY_SORTEM_REAL_2 ( left, right, B, col)
        real*8, target, dimension (:,:), intent (in out) :: B
        integer, intent (in) :: col(:)
        real*8, allocatable, target, dimension ( : )  :: TEMP
        integer, allocatable, target, dimension ( : )  :: IDS
        integer, intent (in) :: left, right
        integer :: i, n, ncol, icol, t, j, i0, i1
        logical :: is_sorted(ubound(col,1))
        
        !Setup indexing arrays
        ncol = ubound(col,1)
        n = right - left + 1
        allocate(IDS(left:right),TEMP(left:right), stat = ALLOC_TEST)
        if(ALLOC_TEST /= 0)then
            write(*,'(/,a,/)') 'ERROR: ALLOCATION FAILED IN ARRAY_SORTEM_INT_2'
            stop
        endif
        do i = left, right
            IDS(i) = i
        enddo
        
        pAr => B
        pTr => TEMP
        pI => IDS
        
        call array_sortem_real_2_r ( left, right, col, 1 )
        
        nullify(pTr)
        nullify(pI)
        nullify(pAr)
        
        deallocate(IDS,TEMP)
    end subroutine ARRAY_SORTEM_REAL_2
    !Recursive sort for array_sortem_int_2, which sets global pointers pA,pT,pI
    recursive subroutine array_sortem_real_2_r ( left, right, cols, icol )
        integer, intent (in) :: left, right !left and right bounds
        integer, intent (in) :: cols(:) !columns being sorted
        integer, intent (in) :: icol !current column index
        integer :: i, j, i0, i1
        
        !Do nothing for 1 element
        if(right == left) return
        
        !Sort the current column and reorder the others
        !The following throws a stack overflow if right-left+1 is big enough
!        pTr(left:right) = pAr(cols(icol),left:right)
        do i = left, right
            pTr(i) = pAr(cols(icol),i)
        enddo
        
        call quick_sort (pTr, pI, left, right)
        RLOOP : do i = 1, ubound(pAr,1)
            !do not re-order a column that is already sorted
            do j = 1, icol-1
                if(i == cols(j)) cycle RLOOP
            enddo
!            pTr(left:right) = pAr(i,left:right)
!            pAr(i,left:right) = pTr(pI(left:right))
            do j = left, right
                pTr(j) = pAr(i,j)
            enddo
            do j = left, right
                pAr(i,j) = pTr(pI(j))
            enddo
        enddo RLOOP
        
        !Compute bounds based on the most recently sorted column
        i0 = left
        j = cols(icol)
        do while (i0 < right)
            i1 = i0
            pI(i0) = i0
            do while ( pAr(j,i1) == pAr(j,i0) )
                i1 = i1 + 1
                if(i1 > right) exit
                pI(i1) = i1
            enddo
            i1 = i1 - 1
            
            !Determine sorted order of elements
            if(icol < ubound(cols,1))then
                call array_sortem_real_2_r ( i0, i1, cols, icol + 1 )
            endif
            
            !Next subset
            i0 = i1 + 1
        enddo

        return
    end subroutine array_sortem_real_2_r

    !! This is a professional version of Quicksort.
    recursive subroutine quick_sort ( A, IDS, LLeft, RRight )
       !! This subroutine implements an improved Quicksort, based on the algorithm by
       !! C. A. R. Hoare.
       !! INCOMING: A      = an array whose elements are to be sorted;
       !!           LLeft  = a pointer to the left-hand end of the partition to be sorted;
       !!           RRight = a pointer to the right-hand end of the partition to be sorted;
       !! OUTGOING: A      = an array whose elements A(LLeft) through A(RRight) are sorted.
       !!           IDS    = an array of indices holding original positions
       real*8, intent (inout) :: A( : )
       integer, intent (in out) :: IDS( : )
       integer, intent (in) :: LLeft, RRight
       integer :: Mid, Left, Right
       real*8    :: Ref, Temp
       integer :: iRef,iTemp, L, R, La
       integer :: No_Swaps, J, Rml
       integer, parameter :: One = 1

       Left = LLeft
       Right = RRight
       if (Right <= Left) then                        !! The partition is empty or contains one element.
          return                                      !! There's no work to do.
       end  if

       !! Section to select the median of three elements, and to move it to the left.
       Mid = (Left + Right)/2
       if (A(Mid)  > A(Right))  then
          call SWAP (A(Mid), A(Right), IDS(Mid), IDS(Right))
       endif

       if (Left+1 == Right) then                      !! There are 2 elements in the partition,
          return                                      !! & they are now in sort.
       endif
       if (A(Left) > A(Mid))  then
          call SWAP (A(Left), A(Mid), IDS(Left), IDS(Mid))
       endif
       if (A(Mid)  > A(Right))  then
          call SWAP (A(Mid), A(Right), IDS(Mid), IDS(Right))
       endif
       if (Left+ 2 == Right) then                     !! There are 3 elements in the partition,
          return                                      !! & they are now in sort.
       endif
       if (A(Mid) == A(Right)) then                   !! Some elements are equal!
          Ref = A(Left)                               !! Forces the left partition to omit equal elements.
         iRef = IDS(Left)
       else
          call SWAP (A(Left), A(Mid), IDS(Left), IDS(Mid))
          Ref = A(Left)                               !! Select the Reference Element.
         iRef = IDS(Left)
       endif

       L = Left
       R = Right + 1

       !! Partition the elements into three groups.
       No_Swaps = 0
       do
          if (L >= R) then
             exit
          endif
          do L = L + 1, R - 1                         !! Scan from the left for an element
             if (A(L) > Ref) then                     !! larger than the Reference.
                exit
             endif
          enddo

          do                                          !! Scan from the right for an element
             R = R - 1
             if (A(R) <= Ref) then                    !! less than or equal to the Reference.
                exit
             endif
          enddo

          if (L < R) then                             !! Swap two elements that are in the wrong partitions.
             Temp = A(R) ; iTemp  = IDS(R)
             A(R) = A(L) ; IDS(R) = IDS(L)
             A(L) = Temp ; IDS(L) = iTemp
             No_Swaps = No_Swaps + 1                  !! Count each swap as we go.
          endif
       enddo
                                                      !! Partitioning is complete.
       if (Left < R) then                             !! Swap the Reference Element into its final position R in the array.
          A(Left) = A(R) ; IDS(Left) = IDS(R)
          A(R)    = Ref  ; IDS(R)    = iRef
       endif
       !! At this point, A(R) is in its correct position in the list.  Elements A(Left) to A(R-1)
       !! are less than or equal to A(R), and elements A(R+1) to A(Right) are greater then A(R).

       !! Section to find out why no swaps were performed.
       if (No_Swaps == 0) then                        !! Something funny happened: not one element was moved. Investigate cause.
          INCREASING: do
             !! Look for any pre-existing order.
             do J = Left, Right-1
                if (A(J) > A(J+1))  then
                   exit  INCREASING
                endif
             enddo
             return !! The elements are already in order.
          enddo INCREASING

       !! Section to take a strong hand when the maximum number of elements is swapped.
       !! It's possible that the elements were in decreasing order.  Check it.
       elseif (No_Swaps+1 == (Right-Left)/2) then
                                                      !! All possible pairs were swapped. Perhaps the elements were in reverse
       DECREASING: do                                 !! order?  Find out why.
             Rml = Right - Left
             if (iand(Rml, One) /= 0) then            !! A partition containing an even number f elements was disarranged during partitioning.
                if (Left < R-1)   then
                   call SWAP (A(Left), A(R-1), IDS(Left), IDS(R-1))        !! Restore order.
                endif
             endif
             do J = Left, Right-1                     !! Check that the elements are sorted.
                if (A(J) > A(J+1)) then
                   exit DECREASING
                endif
             enddo
             return                                   !! The entire sub-list is sorted.
          enddo DECREASING
       endif

       do La = R-1, Left+1, -1                        !! Pass over any elements that are equal to Ref.
          if (A(La) /= Ref) then
             exit
          endif
       enddo
        !! At this point, elements A(La+1) through A(R) are equal.
        !! A(L) is in its correct position too, even
        !! if it is not equal to A(Mid)!
        !! But if Left=La already, the partition
        !! was just lop-sided.

       if (Left < La) then
          call QUICK_SORT (A, IDS, Left, La)               !! Partition the left segment.
       endif
                                                      !! The element at R is in its correct position.
       if (R+1 < Right) then
          call QUICK_SORT (A, IDS, R+1, Right)             !! Partition the right segment.
       endif

    end subroutine QUICK_SORT

    !! This subroutine swaps the element Left_Element with Right_Element.
    subroutine swap ( Left_Element, Right_Element, Left_ID, Right_ID)
       real*8,    intent (in out)     :: Left_Element, Right_Element
       integer, intent (in out)     :: Left_ID, Right_ID
       real*8                         :: Temp
       integer                      :: TempID

       Temp     = Left_Element        ; TempID = Left_ID
       Left_Element   = Right_Element ; Left_ID = Right_ID
       Right_Element  = Temp          ; Right_ID = TempID
    end subroutine swap
    
    !! This is a professional version of Quicksort.
    recursive subroutine quick_sort_int ( A, IDS, LLeft, RRight )
       !! This subroutine implements an improved Quicksort, based on the algorithm by
       !! C. A. R. Hoare.
       !! INCOMING: A      = an array whose elements are to be sorted;
       !!           LLeft  = a pointer to the left-hand end of the partition to be sorted;
       !!           RRight = a pointer to the right-hand end of the partition to be sorted;
       !! OUTGOING: A      = an array whose elements A(LLeft) through A(RRight) are sorted.
       !!           IDS    = an array of indices holding original positions
       integer, intent (inout) :: A( : )
       integer, intent (in out) :: IDS( : )
       integer, intent (in) :: LLeft, RRight
       integer :: Mid, Left, Right
       integer :: Ref, Temp
       integer :: iRef,iTemp, L, R, La
       integer :: No_Swaps, J, Rml
       integer, parameter :: One = 1

       Left = LLeft
       Right = RRight
       if (Right <= Left) then                        !! The partition is empty or contains one element.
          return                                      !! There's no work to do.
       end  if

       !! Section to select the median of three elements, and to move it to the left.
       Mid = (Left + Right)/2
       if (A(Mid)  > A(Right))  then
          call SWAPI (A(Mid), A(Right), IDS(Mid), IDS(Right))
       endif

       if (Left+1 == Right) then                      !! There are 2 elements in the partition,
          return                                      !! & they are now in sort.
       endif
       if (A(Left) > A(Mid))  then
          call SWAPI (A(Left), A(Mid), IDS(Left), IDS(Mid))
       endif
       if (A(Mid)  > A(Right))  then
          call SWAPI (A(Mid), A(Right), IDS(Mid), IDS(Right))
       endif
       if (Left+ 2 == Right) then                     !! There are 3 elements in the partition,
          return                                      !! & they are now in sort.
       endif
       if (A(Mid) == A(Right)) then                   !! Some elements are equal!
          Ref = A(Left)                               !! Forces the left partition to omit equal elements.
         iRef = IDS(Left)
       else
          call SWAPI (A(Left), A(Mid), IDS(Left), IDS(Mid))
          Ref = A(Left)                               !! Select the Reference Element.
         iRef = IDS(Left)
       endif

       L = Left
       R = Right + 1

       !! Partition the elements into three groups.
       No_Swaps = 0
       do
          if (L >= R) then
             exit
          endif
          do L = L + 1, R - 1                         !! Scan from the left for an element
             if (A(L) > Ref) then                     !! larger than the Reference.
                exit
             endif
          enddo

          do                                          !! Scan from the right for an element
             R = R - 1
             if (A(R) <= Ref) then                    !! less than or equal to the Reference.
                exit
             endif
          enddo

          if (L < R) then                             !! Swap two elements that are in the wrong partitions.
             Temp = A(R) ; iTemp  = IDS(R)
             A(R) = A(L) ; IDS(R) = IDS(L)
             A(L) = Temp ; IDS(L) = iTemp
             No_Swaps = No_Swaps + 1                  !! Count each swap as we go.
          endif
       enddo
                                                      !! Partitioning is complete.
       if (Left < R) then                             !! Swap the Reference Element into its final position R in the array.
          A(Left) = A(R) ; IDS(Left) = IDS(R)
          A(R)    = Ref  ; IDS(R)    = iRef
       endif
       !! At this point, A(R) is in its correct position in the list.  Elements A(Left) to A(R-1)
       !! are less than or equal to A(R), and elements A(R+1) to A(Right) are greater then A(R).

       !! Section to find out why no swaps were performed.
       if (No_Swaps == 0) then                        !! Something funny happened: not one element was moved. Investigate cause.
          INCREASING: do
             !! Look for any pre-existing order.
             do J = Left, Right-1
                if (A(J) > A(J+1))  then
                   exit  INCREASING
                endif
             enddo
             return !! The elements are already in order.
          enddo INCREASING

       !! Section to take a strong hand when the maximum number of elements is swapped.
       !! It's possible that the elements were in decreasing order.  Check it.
       elseif (No_Swaps+1 == (Right-Left)/2) then
                                                      !! All possible pairs were swapped. Perhaps the elements were in reverse
       DECREASING: do                                 !! order?  Find out why.
             Rml = Right - Left
             if (iand(Rml, One) /= 0) then            !! A partition containing an even number f elements was disarranged during partitioning.
                if (Left < R-1)   then
                   call SWAPI (A(Left), A(R-1), IDS(Left), IDS(R-1))        !! Restore order.
                endif
             endif
             do J = Left, Right-1                     !! Check that the elements are sorted.
                if (A(J) > A(J+1)) then
                   exit DECREASING
                endif
             enddo
             return                                   !! The entire sub-list is sorted.
          enddo DECREASING
       endif

       do La = R-1, Left+1, -1                        !! Pass over any elements that are equal to Ref.
          if (A(La) /= Ref) then
             exit
          endif
       enddo
        !! At this point, elements A(La+1) through A(R) are equal.
        !! A(L) is in its correct position too, even
        !! if it is not equal to A(Mid)!
        !! But if Left=La already, the partition
        !! was just lop-sided.

       if (Left < La) then
          call QUICK_SORT_INT (A, IDS, Left, La)               !! Partition the left segment.
       endif
                                                      !! The element at R is in its correct position.
       if (R+1 < Right) then
          call QUICK_SORT_INT (A, IDS, R+1, Right)             !! Partition the right segment.
       endif

    end subroutine QUICK_SORT_INT

    !! This subroutine swaps the element Left_Element with Right_Element.
    subroutine swapi ( Left_Element, Right_Element, Left_ID, Right_ID)
       integer, intent (in out)     :: Left_Element, Right_Element
       integer, intent (in out)     :: Left_ID, Right_ID
       integer                      :: Temp
       integer                      :: TempID

       Temp     = Left_Element        ; TempID = Left_ID
       Left_Element   = Right_Element ; Left_ID = Right_ID
       Right_Element  = Temp          ; Right_ID = TempID
    end subroutine swapi

end module quicksort