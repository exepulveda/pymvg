    module ppmt_mod
      implicit none
    !--------------------------------------------------------------------
    ! Modules containing subroutines for the PPMT program.
    ! R. Barnett - 2015
    !--------------------------------------------------------------------

    contains


    subroutine nscore( z, wt, y)
    !--------------------------------------------------------------------
    ! Subroutine for performing a normal score transformation on multiple
    ! input variables.
    ! R. Barnett - 2015
    ! Variables:
    !   z(ndata,nvar) - input data to be transformed
    !   wt(ndata) - input weight attributed to each record
    !   y(ndata,nvar) - output normal scored data
    !--------------------------------------------------------------------
        use quicksort
        use acorn_rng
        use normaldist, only : snorm_inv
        implicit none
        real*8, intent(in) :: z(:,:), wt(:)
        real*8, intent(out) :: y(:,:)
        integer i, j, k, ierr, ndata, nvar, test
        real*8 twt, cp, oldcp, wtfac, w
        real*8, allocatable ::  wt_ns(:), oindex(:), z_ns(:), y1(:,:)
        real*8, parameter :: MACHP = sqrt(epsilon(real(1,8)))
        integer, parameter :: seed = 69069
        logical :: isacorninit = .false.
        type(acorn_type) :: rng

        ! Allocate working arrays
        ndata = ubound( y, 1 )
        nvar = ubound( y, 2 )
        test = 0
        allocate( wt_ns(ndata), oindex(ndata), z_ns(ndata) , y1(ndata,nvar), stat = test)
        if(test /= 0)then; write(*,'(/,a,/)') 'ERROR: ALLOCATION FAILED';stop;endif
        ! Initialize RNG
        if(.not.isacorninit)then
            call acorninit(rng,seed)
            isacorninit = .true.
        endif
        twt = sum( wt )
        wtfac = 1.d0/twt
        ! Proceed through all variables
        do j = 1, nvar
            wt_ns = wt
            ! Perform random despiking
            do i = 1, ndata
              z_ns(i) = z(i,j) + acorni(rng)*MACHP
              oindex(i) = i
            enddo
            ! Sort data by value:
            call qsortem( z_ns, 1, ndata, wt_ns, oindex)
            ! Compute the cumulative probabilities
            oldcp = 0.d0
            cp    = 0.d0
            do i=1,ndata
                cp = cp + wtfac*wt_ns(i)
                wt_ns(i) = (cp + oldcp)/2.d0
                y1(i,j) = snorm_inv( wt_ns(i) )
                oldcp = cp
            end do
            ! Reassign value to entry index
            y( idnint(oindex) , j ) = y1(:,j)
        enddo
        deallocate( wt_ns, oindex, z_ns, y1 )
    end subroutine nscore

    subroutine sphere( y, x, S )
    !--------------------------------------------------------------------
    ! Subroutine for performing a data sphereing transformation.
    ! R. Barnett - 2015
    ! Variables:
    !   y(ndata,nvar) - input data to be transformed
    !   x(ndata,nvar) - output sphere data
    !   S(ndata,nvar) - output sphereing matrix
    !--------------------------------------------------------------------
        implicit none
        real*8, intent(in) :: y(:,:)
        real*8, intent(out) :: x(:,:),S(:,:)
        integer :: i, ndata, nvar, test
        real*8, allocatable, dimension(:) :: diag
        real*8, allocatable, dimension(:,:) :: C, U, D, Dinv, identity

        ! Allocate working arrays
        ndata = ubound( y, 1 )
        nvar = ubound( y, 2 )
        test = 0
        allocate( C(nvar,nvar), U(nvar,nvar), D(nvar,nvar), Dinv(nvar,nvar), &
        & diag(nvar), identity(nvar,nvar), stat = test)
        if(test /= 0)then; write(*,'(/,a,/)') 'ERROR: ALLOCATION FAILED';stop;endif
        ! Calculate the covariance matrix
        call covariance( y, C )
        !print *,"C",C
        ! Perform eigenvalue/eigenvector decomposition of the covariance matrix
        call eig( C, diag, U, nvar )
        !print *,"diag",diag
        !print *,"U",U
        D = 0.d0
        do i = 1, nvar
            D(i,i) = diag(i)
        enddo
        ! Find the inverse of the eigenvalue matrix
        Dinv = D
        call invert(Dinv,nvar)
        !print *,"Dinv",Dinv
        identity = matmul(D,Dinv)
        if (sum(identity)-real(nvar) .gt. 1.d-5)then
            write(*,*) 'Error determining covariance inverse!'
            write(*,*)
            stop
        endif
        Dinv = sqrt(Dinv)
        !print *,"Dinv sqrt",Dinv
        ! Calculate the S sphering matrix
        S = matmul( U, matmul(Dinv,transpose(U)) ) ! Spectral decomposition
        !print *,"S",S

        ! S = matmul( U, Dinv )  ! Dimension reduction
        ! Transform the variables
        x = matmul( y , S )
        deallocate( C, U, D, Dinv )
    end subroutine sphere

    subroutine decorrelate( y, C, x, U )
    !--------------------------------------------------------------------
    ! Subroutine for decorrelating data, given an input covariance matrix
    ! Originally intended for spatially decorrelating data.
    ! R. Barnett - 2015
    ! Variables:
    !   y(ndata,nvar) - input data to be transformed
    !   C(nvar,nvar) - input covariance matrix
    !   x(ndata,nvar) - output transformed data
    !   U(nvar,nvar) - output rotation matrix
    !--------------------------------------------------------------------
        implicit none
        real*8, intent(in) :: y(:,:)
        real*8, intent(in) :: C(:,:)
        real*8, intent(out) :: x(:,:), U(:,:)
        real*8, allocatable, dimension(:) :: diag
        integer ndata, nvar, test

        ! Allocate working arrays
        ndata = ubound( y, 1 )
        nvar = ubound( y, 2 )
        test = 0
        allocate( diag(nvar),  stat = test)
        if(test /= 0)then; write(*,'(/,a,/)') 'ERROR: ALLOCATION FAILED';stop;endif
        ! Perform eigenvalue/eigenvector decomposition of the covariance matrix
        call eig(C,diag,U,nvar)
        ! Transform the variables
        x = matmul( y , U )
        deallocate( diag )
    end subroutine decorrelate

    subroutine covariance( z, C )
    !--------------------------------------------------------------------
    ! Subroutine for calculating data covariance.
    ! R. Barnett - 2015
    ! Variables:
    !   z(ndata,nvar) - input data matrix
    !   C(nvar,nvar) - output covariance matrix
    !--------------------------------------------------------------------
        implicit none
        integer :: nvar, ndata, test, i, j
        real*8, intent(in), dimension(:,:) :: z
        real*8, intent(out), dimension(:,:) :: C
        real*8, allocatable, dimension(:) :: mean

        ndata = size( z, 1 )
        nvar = size( z, 2 )
        allocate( mean(nvar), stat = test)
        if(test /= 0)then; write(*,'(/,a,/)') 'ERROR: ALLOCATION FAILED';stop;endif
        do i = 1,nvar
            mean(i) = sum( z(:,i) ) / dble(ndata)
        enddo
        do i = 1,nvar
            do j = 1,nvar
                C(i,j) = sum( (z(:,i)-mean(i))*(z(:,j)-mean(j)) ) / dble(ndata)
            enddo
        enddo
        deallocate( mean )
    end subroutine covariance

    subroutine correlation( z, C )
    !--------------------------------------------------------------------
    ! Subroutine for calculating data correlation (pearson).
    ! R. Barnett - 2015
    ! Variables:
    !   z(ndata,nvar) - input data matrix
    !   C(nvar,nvar) - output correlation matrix
    !--------------------------------------------------------------------
        implicit none
        integer :: nvar, ndata, test, i, j
        real*8, intent(in), dimension(:,:) :: z
        real*8, intent(out), dimension(:,:) :: C
        real*8, allocatable, dimension(:) :: mean, std

        ndata = size( z, 1 )
        nvar = size( z, 2 )
        allocate( mean(nvar), std(nvar), stat = test)
        if(test /= 0)then; write(*,'(/,a,/)') 'ERROR: ALLOCATION FAILED';stop;endif
        do i = 1,nvar
            mean(i) = sum( z(:,i) ) / dble(ndata)
            std(i) = sqrt( sum( (z(:,i) - mean(i))**2 )/dble(ndata) )
        enddo
        do i = 1,nvar
            do j = 1,nvar
                C(i,j) = sum( (z(:,i)-mean(i))*(z(:,j)-mean(j)) ) / dble(ndata)
                C(i,j) = C(i,j) / (std(i)*std(j));
            enddo
        enddo
        deallocate( mean, std )
    end subroutine correlation

    subroutine orthogs(A)
    !--------------------------------------------------------------------
    ! Subroutine for calculating a matrix of orthogonal columns
    ! using Gram Schmidt algorithm (first column left intact)
    ! R. Barnett - 2015
    ! Variables:
    !   A(n,n) - input/output matrix.
    !--------------------------------------------------------------------
        implicit none
        real*8, intent(inout) :: A(:,:)
        integer n, i, j, test
        real*8, allocatable, dimension(:,:) :: v, vt, vo

        n = size( A, 1 )
        allocate( v(n,n), vt(1,n), vo(n,1), stat = test)
        if(test /= 0)then; write(*,'(/,a,/)') 'ERROR: ALLOCATION FAILED';stop;endif

        v = 0.d0
        v(:,1) = A(:,1)/norm2(A(:,1))
        do j = 2, n
            v(:,j) = A(:,j)
            do i = 1, j-1
                vo(:,1) = v(:,i)
                vt(1,:) = v(:,i)
                v(:,j) = v(:,j) - matmul( vo, matmul(vt,A(:,j)) )
            enddo
            v(:,j) = v(:,j)/norm2(v(:,j))
        enddo
        A = v
        deallocate( v, vt, vo )
    endsubroutine orthogs

    subroutine projection( c, y, z, r )
    !--------------------------------------------------------------------
    ! Subroutine for calculating the projection of multivariate data
    ! on a given vector
    ! Optionally - calculates the r Gaussian test statistics
    ! R. Barnett - 2015
    ! Variables:
    !   c(1,nvar) - input vector matrix.
    !   y(ndata,nvar) - input data matrix.
    !   z(ndata,1) - output projection
    !   r(ndata,1) - output test statistics
    !--------------------------------------------------------------------
        use normaldist, ONLY : snorm_cdf
        ! project the data based on linear combination c
        implicit none
        real*8, intent(in) :: y(:,:), c(:,:)
        real*8, intent(out) :: z(:,:)
        real*8, optional, intent(out):: r(:,:)
        integer ndata, nvar, i
        ndata = size( y, 1 )
        nvar = size( y, 2 )
        if( ubound(c,1) == nvar )then
            z = matmul( y, c )
        elseif( ubound(c,2) == nvar )then
            z = matmul( y, transpose(c) )
        else
            write(*,*) ' Error - dimensions of data and vector do not align!'
            stop
        endif
        if(present(r))then
            do i=1,ndata
                r(i,1) = 2.d0 * snorm_cdf ( dble( z(i,1) ) ) - 1.d0
            enddo
        endif
    end subroutine projection

    real*8 function projindex( c, y, w, l ) result ( ix )
    !--------------------------------------------------------------------
    ! Subroutine for calculating the projection index of multivariate
    ! data on a vector
    ! R. Barnett - 2015
    ! Variables:
    !   c(1,nvar) - input vector matrix.
    !   y(ndata,nvar) - input data matrix.
    !   w(ndata) - weight of each data.
    !   l - number of legrendre polynomial expansions
    !   ix - output projection index for this vector and data
    !--------------------------------------------------------------------
        implicit none
        real*8,intent(in) :: y(:,:), c(:,:), w(:)
        integer, intent(in) :: l
        integer ndata, nvar, i, j, test
        real*8, allocatable, dimension(:,:) ::  p, z, r
        real*8, allocatable, dimension(:) ::  E, TT

        ndata = size( y, 1 )
        nvar = size( y, 2 )
        allocate( p(ndata,l+1), E(l+1), z(ndata,1), r(ndata,1), stat = test)
        if(test /= 0)then; write(*,'(/,a,/)') 'ERROR: ALLOCATION FAILED';stop;endif

        if (.false.) then
          print *, "Inputs display:"

          print *,'C:'
          do i=1,nvar
            print *,i,'C',C(1,i)
          end do

          print *,'Y:'
          do i=1,ndata
            do j=1,nvar
              print *,i,j,'Y',Y(i,j)
            end do
          end do

          print *,'W:'
          do i=1,ndata
            print *,i,'W',W(i)
          end do

          print *,'L:'
          print *,l

          stop
        end if

        call projection( c, y, z, r )

        !print *,'projection computing'
        !print *,'Z/R'
        !do i=1,ndata
        !  print *,i,',',z(i,1),',',r(i,1)
        !end do

        p(:,1) = 1 !Legendre 0
        p(:,2) = r(:,1) !Legendre 1
        do i = 3, l+1
            j = i-1
            p(:,i) = (( 2.d0*dble(j) - 1.d0 ) * r(:,1) * p(:,i-1) &
            &           - ( dble(j)-1.d0)*p(:,i-2) )/dble(j) * w(:) !Legendre >=2
        enddo

        !print *, "LEGENDRE:"
        !do i=1,ndata
        !  write(*,'(2000f22.6)')( p(i,j),j=1,l+1)
        !enddo

        !print *, "E:"
        do i = 1, l+1
            E(i) = ( sum( p(:,i) ) / sum(w) )**2.d0 !E do equation 23
        !    print *,E(i)
        enddo

        !print *, "IX:"
        ix = 0.d0
        do i = 2, l+1
            ix = ix + (2.d0*dble(i) + 1.d0)/2.d0 * E(i) ! I (Projection index)
        !    print *,i,ix
        enddo



        !stop


        deallocate( p, E, z, r )
    end function projindex

    subroutine gradient( c, y, w, l, E )
    !--------------------------------------------------------------------
    ! Subroutine for calculating the gradient of the projection index
    ! based on derivatives of the legendre polynomials
    ! R. Barnett - 2015
    ! Variables:
    !   c(1,nvar) - input vector matrix.
    !   y(ndata,nvar) - input data matrix.
    !   w(ndata) - weight of each data.
    !   l - number of legrendre polynomial expansions
    !   ix - output projection index for this vector and data
    !--------------------------------------------------------------------
        implicit none
        real*8, intent(in) :: y(:,:), c(:,:), w(:)
        integer, intent(in) :: l
        real*8, intent(out) :: E(:,:)
        integer ndata, nvar, i, j, k, test
        real*8, allocatable, dimension(:,:) ::  p, pd, z, r
        real*8, parameter :: pi = 3.141592653589793238462643383279502884197

        ndata = size( y, 1 )
        nvar = size( y, 2 )
        allocate( p(ndata,l+1), pd(ndata,l), z(ndata,1), r(ndata,1), stat = test)
        if(test /= 0)then; write(*,'(/,a,/)') 'ERROR: ALLOCATION FAILED';stop;endif

        call projection( c, y, z, r)
        ! Compute P
        p(:,1) = 1 !Legendre 0
        p(:,2) = r(:,1) !Legendre 1
        do i = 3, l+1
            j = i-1
            p(:,i) = (( 2.d0*dble(j) - 1.d0 ) * r(:,1) * p(:,i-1) &
            &           - ( dble(j)-1.d0)*p(:,i-2) )/dble(j) * w(:) !Legendre >=2
        enddo
        ! Compute Pd (derivative)
        pd(:,1) = 1 !Legendre 1 derivative (eq12 - pg252)
        do i = 2,l
             pd(:,i) = (r(:,1)*pd(:,i-1) + i*p(:,i-1)) * w(:); !Legendre >=1
        enddo
        E = 0.d0
        do k = 1, nvar
            do i = 1, l
                E(1,k) = ( 2.d0*dble(i)+1 ) * sum( p(:,i)/sum(w) ) *&
                &( sum( pd(:,i) * exp(-z(:,1)**2/2.d0) * ( y(:,k)-c(1,k)*z(:,1) ) ) / sum(w) ) !equation 11
            enddo
        enddo
        E = E*2 / (sqrt(2.d0*pi))
        deallocate( p, pd, z, r )
    end subroutine gradient

    end module ppmt_mod
