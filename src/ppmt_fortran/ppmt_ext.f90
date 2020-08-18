
!***********************************************************************************!
!                                                                                   !
!  PROGRAM: PPMT                                                                    !
!                                                                                   !
!  PURPOSE:  Projection Pursuit Multivariate Transform                              !
!                                                                                   !
!  AUTHOR: Ryan Barnett and John Manchuk 2015                                       !
!                                                                                   !
!  PAPER: Barnett, R., Manchuk, J., & Deutsch, C. (2012). Projection Pursuit        !
!  Multivariate Transform. CCG Annual Report 14, Paper103.                          !
!                                                                                   !
!  ACKNOWLEDGMENTS: Based largely on the algorithm originally outlined in the paper:!
!  Friedman, J. (1987). Exploratory Projection Pursuit. Journal of the American     !
!  Statistical Association, vol.82, pp.249-266.                                     !
!                                                                                   !
!  Version 3.1 performs an optional spatial decorrelation (MAF) transform following !
!              projection pursuit                                                   !
!***********************************************************************************!

module ppmt_ext
implicit none

contains

function pp_index(direction, data, w, lo) result ( ix )
  use ppmt_mod
  implicit none

  real*8, intent(in) :: direction(:)       !covariates to fit and transform
  real*8, intent(in) :: data(:,:)       !covariates to fit and transform
  real*8, intent(in) :: w(:)     !the sphering matrix
  integer, intent(in) :: lo
  real*8 :: ix
  !locals
  real*8 :: c(1,size(direction))
  !integer i

  !print *, size(direction)
  !print *, size(data,1)
  !print *, size(data,2)

  !do i=1,size(data,2)
  !  print *,sum(data(:,i)) / size(data,1)
  !end do

  !stop
  c(1,:) = direction
  ix = projindex(c,data,w,lo)
end function

subroutine ppmt_fit(data,min_index,maxiters,iterations, &
    y,                 & !data transformed
    snorm_out,         & !normal score used for Gaussianisation
    z_first,             & !first Gaussianisation before iteration
    sph_out,           & !sphering matrix
    sph_inv_out,       & !inverse of sphering matrix
    U_out,             & !Orthogonal rotation matrices
    projections_out,   & !projections
    zd,                & !final normal score
    seed,              & !seed number
    trace)
!f2py intent(inout) y
!f2py intent(inout) sph
!f2py intent(inout) sph_inv
!f2py intent(inout) directions
!f2py intent(inout) projections
!f2py intent(inout) normal_scores
  use normaldist
  use ppmt_mod
  use quicksort

  implicit none
  real*8, intent(in)    :: data(:,:)
  real*8, intent(in)    :: min_index
  integer, intent(in)   :: maxiters
  integer, intent(out)  :: iterations         !number of final iterations
  real*8, intent(inout) :: y(:,:)
  real*8, intent(inout) :: snorm_out(:)
  real*8, intent(inout) :: z_first(:,:)
  real*8, intent(inout) :: sph_out(:,:)
  real*8, intent(inout) :: sph_inv_out(:,:)
  real*8, intent(inout) :: U_out(:,:,:)
  real*8, intent(inout) :: projections_out(:,:)
  real*8, intent(inout) :: zd(:,:) 
  integer, intent(in)   :: seed         !seed number
  integer, intent(in)   :: trace         !verbosity

  !locals
  real*8, allocatable :: y1(:,:)
  real*8, allocatable :: x(:,:)
  real*8, allocatable :: w(:,:)
  !real*8, allocatable :: Am(:,:)
  real*8, allocatable :: R(:,:)
  real*8, allocatable :: ei(:,:)
  real*8, allocatable :: a(:,:)
  real*8, allocatable :: atest1(:,:)
  real*8, allocatable :: atest2(:,:)
  real*8, allocatable :: id(:,:)
  real*8, allocatable :: gradtest(:,:)
  real*8, allocatable :: U(:,:)
  real*8, allocatable :: Ux(:,:)
  real*8, allocatable :: indexi(:)
  real*8, allocatable :: r8temp2(:,:)
  real*8, allocatable :: eit(:,:)
  real*8, allocatable :: pr(:,:)
  real*8, allocatable :: eye(:,:)
  !locals not allocatable
  integer test
  real*8              :: c(1,1)
  real*8              :: c1, c2
  integer exc
  integer i,j,inner_iter, iter1, iters_type
  real*8 maxpix, maxpix0,pix,pixm,pixp,s,step
  integer ndata, nvar,nt
  logical useweight

  useweight = .false.

  ndata = size(data,1)
  nvar = size(data,2)

  if (trace>=10) then
    print *,'ndata',ndata
    print *,'nvar',nvar
  end if

  ! Allocate and populate some working arrays
  allocate(y1(ndata,nvar), x(ndata,nvar), w(ndata,2),stat = test )

  if(test /= 0)then; write(*,'(/,a,/)') 'ERROR: ALLOCATION FAILED';stop;endif

  z_first = data
  !
  ! Normal score preprocessing (always equally weighted)
  !
  if (trace >= 1) then
    write(*,*)
    write(*,'(a,a)') ' normal scoring the data'
  end if
  w(:,1) = 1 ! No normal score weights are applied prior to the PPMT (applied after if necessary)
  call nscore( z_first , w(:,1), y )
  y1 = y ! Store the univariate Gaussian data for writing out later
  ! These standard normal Gaussian values will be referenced later for Gaussianizing
  snorm_out = y(:,1)
  call qsortem( snorm_out , 1 , ndata )
  do i = 1, nvar
    call qsortem( z_first(:,i) , 1 , ndata )
  enddo
  !print *,"snorm min.max",snorm_out(1),snorm_out(ndata)
  !write(ltrn) snorm
  !write(ltrn) z !ndata x nvar
  !deallocate(z)
  !
  ! Sphere the data
  !
  if (trace >= 1) then
    write(*,*)
    write(*,'(a,a)') ' sphereing the data'
  end if
  call sphere( y  , x , sph_inv_out )
  !print *,'sphering matrix',sph
  sph_out = sph_inv_out
  call invert ( sph_inv_out , nvar )
  ! Write out the inverse of the Sphere matrix (needed for back-transform rather than Sphere)
  !write(ltrn) sph
  !print *,'sphering inverse matrix',sph
  !stop
  !
  ! Initialize looping variables
  !
  allocate(R(nvar,nvar),ei(nvar,nvar),a(1,nvar),atest1(1,nvar),atest2(1,nvar),id(nvar,nvar),&
  & gradtest(1,nvar),U(nvar,nvar),Ux(ndata,nvar),indexi(ndata),r8temp2(ndata,1),&
  & eit(nvar,1),pr(ndata,1), eye(1,nvar), stat = test)
  if(test /= 0)then; write(*,'(/,a,/)') 'ERROR: ALLOCATION FAILED';stop;endif
  u=0.d0;id=0.d0;ei=0.d0;
  do i = 1, nvar
      id(i,i) = 1
      ei(i,i) = 1
  enddo
  ! A random basis will be used intially for R
  ! pseudo-random number generator is good enough for this purpose
  call srand(seed)
  do i=1,nvar
      do j=1,nvar
          R(i,j)=rand()
      enddo
  enddo
  call orthogs(R)
  do i = 1, ndata
      indexi(i) = i
  enddo

  !
  ! Main projection pursuit loop
  !
  write(*,*)
  write(*,'(a,a)') ' performing projection pursuit iterations'
  exc = 0
  iter1 = 0

  !ispatial = 0
  iters_type = 0
  do nt=1,maxiters
      if (trace >= 1) write(*,'(a,i6)') 'performing projection pursuit at iteration:',nt
      !ispatial = ispatial + 1

      ! First determine the maximum direction along the coordinate axes
      ! Since X is sphered, the coordinate axes are the principal components
      ! of the original data Y
      maxpix = 0.d0
      do i = 1,nvar
          eye(1,:) = ei(i,:)
          !do j=1,5
          !  print *,j,'XX ',x(j,1),x(j,2)
          !end do

          pix = projindex( eye, x, w(:,1), 10 )
          if( abs(pix) > maxpix )then
              a(1,:) = ei(i,:)
              maxpix = abs(pix)
          endif
          !print *, 'axis',i,'eye',eye,'pi',pix
      enddo
      ! Now perform a coarse stepping around the most interesting coordinate
      ! axes according to Freidman '85 (pg256)
      if (trace >= 5) write(*,'(a,f8.6)') ' maximum projection index stage 1:',  maxpix

      maxpix0 = 0.d0
      inner_iter = 0;
      do while ( (maxpix-maxpix0) >= 1.0d-20 .and. inner_iter < 100 )
          maxpix0 = maxpix
          do j = 1,nvar
              eit(:,1) = ei(j,:)
              c = (1.d0+ matmul( a,eit ))**0.5d0
              c1 = c(1,1)
              atest1(1,:) = 1/sqrt(2.0)*( a(1,:)+ei(j,:) ) / c1
              pixp = projindex( atest1, x, w(:,1), 10)
              !print *, 'coarse opt::inner_iter',inner_iter,j,'ei',ei(j,:),'c1',c1,'a',atest1(1,:),'pi',pixp
              c = (1.d0- matmul( a,eit ))**0.5d0
              c2 = c(1,1)
              atest2(1,:) = 1/sqrt(2.0)*( a(1,:)-ei(j,:) )/ c2
              pixm = projindex( atest2, x, w(:,1), 10 )
              !print *, 'coarse opt::inner_iter',inner_iter,j,'ei',ei(j,:),'c2',c2,'a',atest2(1,:),'pi',pixm

              if( isnan( pixp ) .or. any(isnan( atest1 )) ) pixp = 0.d0
              if( isnan( pixm ) .or. any(isnan( atest2 )) ) pixm = 0.d0
              if( pixp == 0.d0 .and. pixm == 0.d0 ) cycle
              if( pixp> pixm )then
                  pix = pixp;s = 1.d0 !ones(1,length(a));
              else
                  pix = pixm;s = -1.d0 !*ones(1,length(a));
              endif
              if( pix > maxpix )then
                  a(1,:) = 1/sqrt(2.d0)*(a(1,:)+s*ei(j,:)) / (1+s*a(1,j))**0.5
                  maxpix = pix
              endif

              if (trace >= 5) write(*,'(a,i5,i5,f12.6)') ' maximum projection index stage 1: ',inner_iter,j, maxpix

          enddo
          inner_iter = inner_iter + 1
      enddo
      !stop "AAA"
      if (trace >= 5) write(*,'(a,f8.6)') 'maximum projection index stage 2:',  maxpix

      ! Now apply gradient descent to maximize on a smaller scale
      step = 0.1d0
      inner_iter = 0
      maxpix0 = 0.d0
      gradtest = 0.d0
      do while ( (maxpix-maxpix0) >= 1.0d-20 .and. inner_iter < 100 )
          maxpix0 = maxpix
          !Apply gradient descent
          call gradient( a, x, w(:,1), 10, gradtest )
          atest1=a+step*gradtest
          ! Update objective
          pix = projindex( atest1, x, w(:,1), 10 )
          if( pix >  maxpix )then
              a = a + step*gradtest
              maxpix = pix
              if (trace >= 5) write(*,'(a,f14.12,i5)') 'new maximum projection index reached:',  maxpix,inner_iter
          endif
          inner_iter = inner_iter + 1
      enddo
      !On occasion, projection of full zeros passes optimization
      !The following hotfix fixes the anomaly... must be better understood
      if(maxpix > 1)then
          exc = exc + 1
          if( exc > nvar ) exc = 1
          a = 0.d0
          a(1,exc) = 1.d0
      else
          if (trace >= 5) write(*,'(a,f8.6)') 'at the end of iteration, maximum projection index:',  maxpix
      endif

      if (trace >= 5) write(*,'(a,i4,f8.6,200(f8.2))') 'END ITERATION:',nt,maxpix,(a(1,j), j=1,nvar)

      !print *,'pi',maxpix


      !stop "FIN ITER 1"

      call projection( a, x, pr )

      !Gaussianize along the projection
      U = R
      U(:,1) = a(1,:)
      call orthogs(U)!orthogonolize with Gram-Schmidt
      Ux = matmul(x,U)
      r8temp2(:,1) = indexi
      call qsortem( pr(:,1), 1, ndata, r8temp2(:,1) )
      Ux(int( r8temp2(:,1) ),1) = snorm_out !nscore along projection
      x = matmul(Ux,transpose(U))   !unproject

      !store the outputs
      U_out(:,:,nt) = U           !orthogonal basis
      projections_out(:,nt) = pr(:,1) !projection


      !Write information to the transformation table
      !EXE itertype = 54321 ! Code for projection pursuit
      !write(ltrn) itertype
      !write(ltrn) U
      !write(ltrn) pr(:,1)

      ! Stopping criteria
      if( maxpix <= min_index)then
        if (trace >= 5) write(*,*) 'PI tolerance reached'
        exit
      endif
  enddo

  iterations = min(nt,maxiters)
  !
  ! Iterations complete, report stopping criteria
  !
  if( maxpix <= min_index) then
      write(*,*)
      write(*,'(a,i4,a,f0.4)') 'PP finished based on targeted Gaussianity index at iteration:', &
        iterations,'. pi=',maxpix
  else
      write(*,*)
      write(*,'(a,i4,a,f0.4)') 'PP finished based on max iters.', &
        iterations,'. pi=',maxpix
  endif

  zd = x
  !
  ! Normalize the margins - apply declustering weights if necessary
  !
  !EXE itertype = 6789 ! Code for weighted transform
  !write(ltrn) itertype
  if( useweight )then
      call nscore( x , w(:,2) , y )
  else
      call nscore( x , w(:,1) , y )
  endif

  deallocate(y1,x,w,R,ei,a,atest1,atest2,id,&
    gradtest,U,Ux,indexi,r8temp2,&
    eit,pr,eye)
end subroutine

subroutine ppmt_transform(data,& !data to transform
    iterations, &
    snorm,             & !normal score used for Gaussianisation
    zd,                & !first Gaussianisation before iteration
    sph,               & !sphering matrix
    pr_input,          & !original projections
    Us,                & !Orthogonal rotation matrices
    y,                 & !data transformed
    trace)
!f2py intent(inout) y
  use normaldist
  use ppmt_mod
  use quicksort
  use backtrmod
  implicit none
  real*8, intent(in)    :: data(:,:)
  real*8, intent(inout) :: y(:,:)
  integer, intent(in)   :: iterations         ! the number of iterations of the pp loop
  real*8, intent(in)    :: snorm(:) !the inverse of the sphering matrix
  real*8, intent(in)    :: pr_input(:,:) !projections
  real*8, intent(in)    :: zd(:,:)       !first nscore
  real*8, intent(in)    :: sph(:,:) !the inverse of the sphering matrix
  real*8, intent(in)    :: Us(:,:,:) !the inverse of the sphering matrix
  integer, intent(in)   :: trace         !verbosity

  !locals
  real*8, allocatable :: x(:,:)
  real*8, allocatable :: w(:)
  real*8, allocatable :: a(:,:)
  real*8, allocatable :: U(:,:)
  real*8, allocatable :: Ux(:,:)
  real*8, allocatable :: indexi(:)
  real*8, allocatable :: r8temp2(:,:)
  real*8, allocatable :: pr(:,:)
  real*8, allocatable :: p_trans(:)
  !locals not allocatable
  integer test
  integer i,j
  integer ndata, nvar,nt,utail,ltail
  logical useweight
  real*8 :: ymin,ymax,utpar,ltpar

  useweight = .false.

  ndata = size(data,1)
  nvar = size(data,2)

  if (trace>=10) then
    print *,'ndata',ndata
    print *,'nvar',nvar
  end if

  ! Allocate and populate some working arrays
  allocate(x(ndata,nvar),w(ndata),stat = test)

  if(test /= 0)then; write(*,'(/,a,/)') 'ERROR: ALLOCATION FAILED';stop;endif

  !
  ! Step 1) Apply normal score preprocessing
  !
  if (trace >= 1) then
    write(*,*)
    write(*,'(a,a)') ' normal scoring the data'
  end if
  w = 1 ! No normal score weights are applied prior to the PPMT (applied after if necessary)

  ymin = -6.0
  ymax =  6.0

  utail = 1
  utpar = 1.95
  ltail = 1
  ltpar = 1.95

  do i=1,nvar
    do j = 1, ndata
        print *, data(j,i),minval(zd(:,i)),maxval(zd(:,i)),ymin,ymax
      y(j,i) = forward(data(j,i),zd(:,i),snorm)
    enddo
  enddo

  print *, sum(y(:,1))
  print *, sum(y(:,2))
  print *, sum(snorm)
  !stop

  !
  ! Step 2) Sphere the data
  !
  if (trace >= 1) then
    write(*,*)
    write(*,'(a,a)') ' sphereing the data'
  end if

  x = matmul(y,sph)

  !x now holds the current transformed data after sphering

  !
  ! Initialize looping variables
  !
  allocate(a(1,nvar),p_trans(ndata), &
  & U(nvar,nvar),Ux(ndata,nvar),indexi(ndata),r8temp2(ndata,1),&
  & pr(ndata,1), stat = test)
  if(test /= 0)then; write(*,'(/,a,/)') 'ERROR: ALLOCATION FAILED';stop;endif

  do i = 1, ndata
      indexi(i) = i
  enddo

  !
  ! Step 3) Main projection pursuit loop
  !
  write(*,*)
  write(*,'(a,i6)') 'performing projection pursuit iterations:',iterations

  do nt=1,iterations
      if (trace >= 1) write(*,'(a,i6)') 'performing projection pursuit at iteration:',nt

      ! Step 3a) Project onto direction

      a(1,:) = Us(:,1,nt) !extract direction a from saved Us

      call projection( a, x, pr )

      ! Step 3b) Gaussianize along the projection
      U = Us(:,:,nt)
      Ux = matmul(x,U)
      r8temp2(:,1) = indexi
      call qsortem( pr(:,1), 1, ndata, r8temp2(:,1) )

      do i=1,ndata
        p_trans(i) = forward(pr(i,1), pr_input(:,nt), snorm)
      enddo


      Ux(int( r8temp2(:,1) ),1) = p_trans !nscore along projection
      x = matmul(Ux,transpose(U))   !unproject
  enddo

  !
  ! Step 4) Normalize the margins - apply declustering weights if necessary
  !
  call nscore(x,w,y)

  deallocate(x,w,a,U,Ux,indexi,r8temp2,pr,p_trans)

end subroutine

subroutine ppmt_back_transform(data, back_data, & !data to back transform and back transformed
    iterations,        &    !iterations made
    snorm,             & !normal score used for Gaussianisation
    zd,                & !first Gaussianisation before iteration
    sph_inv,           & !inverse of sphering matrix
    Us,                & !Orthogonal rotation matrices
    pr,                & !projections
    zd2_in,            & !final normal score
    yd2_in,            & !final normal score
    trace)
  use normaldist
  use ppmt_mod
  use quicksort
  use backtrmod
  implicit none
  !f2py intent(inout) back_data


  real*8, intent(in) :: data(:,:)       !covariates to fit and transform
  real*8, intent(inout) :: back_data(:,:)       !covariates to fit and transform
  integer, intent(in) :: iterations
  real*8, intent(in) :: snorm(:) !the inverse of the sphering matrix
  real*8, intent(in) :: zd(:,:)       !covariates to fit and transform
  real*8, intent(in) :: sph_inv(:,:) !the inverse of the sphering matrix
  real*8, intent(in) :: Us(:,:,:) !the inverse of the sphering matrix
  real*8, intent(in) :: pr(:,:) !the inverse of the sphering matrix
  real*8, intent(in) :: zd2_in(:,:) !the inverse of the sphering matrix
  real*8, intent(in) :: yd2_in(:,:) !the inverse of the sphering matrix
  integer, intent(in)   :: trace         !verbosity

  !locals
  real*8, allocatable    :: xt(:,:)       !covariates to fit and transform
  real*8, allocatable :: U(:,:)
  real*8, allocatable :: Ux(:,:)
  real*8, allocatable :: zd2(:,:)
  real*8, allocatable :: yd2(:,:)
  !locals not allocatable
  integer test
  real*8              :: ymin,ymax,zmin,zmax
  !real*8              :: c1, c2
  !integer exc
  integer i,j, np, np0, utail,ltail
  real*8 utpar,ltpar,pdelta
  integer ndata, nvar
  logical useweight

  np = size(data,1)
  nvar = size(data,2)
  ndata = size(snorm)
  np0 = np
  print *, 'np',np
  print *, 'nvar',nvar

  print *, 'iterations',iterations

  useweight = .false.


  !Determine extrapolation values
  ymin = snorm(1)  !snorm_inv ( dble(1/dble(nsim)) )
  ymax = snorm(ndata)   !snorm_inv ( (dble(nsim)-1.0)/dble(nsim) )

  print *, 'ymin',ymin
  print *, 'ymax',ymax


  allocate( xt(np,nvar), Ux(np,nvar),U(nvar,nvar),zd2(ndata,nvar),yd2(ndata,nvar), stat=test )

  zd2 = zd2_in
  yd2 = yd2_in
  do i = 1,nvar
      call qsortem(yd2(:,i),1,ndata,zd2(:,i))
  enddo


  !allocate( data_(nvari,nxyz), ixp(nxyz), stat=test ) ! This doesn't change in the iterations
  if(test /= 0)then; write(*,'(/,a,/)') 'ERROR: ALLOCATION FAILED';stop;endif

  utail = 2
  utpar = 1.95
  ltail = 2
  ltpar = 1.95
  pdelta = 1.d0 /dble(ndata) - 1.d0/dble(np)

  xt = data
  !Back-nscore
  do i=1,nvar
      zmin = minval( zd2(:,i) )
      zmax = maxval( zd2(:,i) )
      do j = 1, np
          back_data(j,i) = backtr( xt(j,i), zd2(:,i), yd2(:,i), zmin, zmax, ltail, ltpar, utail, utpar )
      enddo
  enddo

  ! Proceed backward through the iterations
  do i = iterations, 1, -1
    ! Grab the correct vector
    U = Us(:,:,i)
    ! Project the data along the vector
    Ux = matmul(xt,U)
    ! Back-transform the projection
    xt(:,1) = Ux(:,1)
    zmin = minval( pr(:,i) )
    zmax = maxval( pr(:,i) )
    do j = 1, np
        Ux(j,1) = backtr( xt(j,1), pr(:,i), snorm, zmin, zmax, ltail, ltpar, utail, utpar )
    enddo
    ! Unproject - could do this in one step but likely to cause allocation errors
    do j = 1, np
        back_data(j,:) = matmul( Ux(j,:), transpose(U))
    enddo
    xt = back_data
  enddo
  !Back-sphere (again - do one by one
  do i = 1, np
      back_data(i,:) = matmul( xt(i,:), sph_inv)
  enddo
  !Back-nscore
  xt = back_data
  do i=1,nvar
      zmin = minval( zd(:,i) )
      zmax = maxval( zd(:,i) )
      do j = 1, np
          back_data(j,i) = backtr( xt(j,i), zd(:,i), snorm, zmin, zmax, ltail, ltpar, utail, utpar )
      enddo
  enddo

  deallocate(xt,Ux,U,zd2)


end subroutine

subroutine ppmt_extract_param_dim(transfile, &
    nvar,      &
    ndata,     &
    iterations)
  implicit none

  character(len=*), intent(in) :: transfile
  integer, intent(out)   :: nvar,ndata
  integer, intent(out)   :: iterations

  !locals
  integer lin,len,nvar_tmp,ndata_tmp,iters,iters0,iters1,itertype,i
  character(len=1024) :: str
  real*8, allocatable :: snorm(:),zd(:,:),sph_inv(:,:), Us(:,:),pr(:),Cs(:,:)
  integer :: test
  !begin

  open (unit=lin, file = transfile, status = 'old', access = 'stream' )

  read( unit=lin,ERR = 97) len
  read( lin) str(1:len)
  if( str(1:len) /= 'PPMT transformation table file' ) then
     write(*,*) trim(str(1:len))
     go to 97
  end if

  read( lin, err = 97 ) nvar

  read( lin, err = 97 ) ndata

  allocate(snorm(ndata),zd(ndata,nvar),sph_inv(nvar,nvar), Us(nvar,nvar),pr(ndata),Cs(nvar,nvar), stat=test )
  if(test /= 0)then; write(*,'(/,a,/)') 'ERROR: ALLOCATION FAILED';goto 97;endif

  !Read in the normal score preprocessing info
  read( lin, err = 97 ) snorm
  read( lin, err = 97 ) zd
  !Read in the inverse sphere matrix
  read( lin, err = 97 ) sph_inv

  iters = 0
  iters0 = 0 ! Counter for spatial decorrelation
  iters1 = 0 ! Counter for projection pursuit
  !spatialiter = .false.
  do
      if( iters > 500 )then
          write(*,*)  ' PPMT back-transform only handles 500 iterations for now....'
          stop
      endif
      read( lin, err = 97 ) itertype
      if( itertype == 12345 )then
          iters0 = iters0 + 1 ! Spatial decorrelation
          read( lin, err = 97 ) Cs
      elseif( itertype == 54321 )then
          iters1 = iters1 + 1 ! Projection pursuit
          read( lin, err = 97 ) Us
          read( lin, err = 97 ) pr
      elseif( itertype == 6789 )then
          exit ! Final marginal normal score
      else
          go to 97 ! Shouldn't get here - error
      endif
      iters = iters + 1
  end do

  iterations = iters

  goto 99

97 iterations = 0
   ndata = 0
   nvar = 0

99 close(lin)
  deallocate(snorm,zd,sph_inv,Us,pr,Cs)
end subroutine

subroutine ppmt_extract_param(transfile, &
    nvar,      &
    ndata,     &
    max_iterations,&
    iterations,&
    snorm,     & !normal score used for Gaussianisation
    zd,        & !first Gaussianisation before iteration
    sph_inv,   & !inverse of sphering matrix
    Us,        & !Orthogonal rotation matrices
    pr,        & !projections
    zd2,       & !final normal score
    yd2,       & !final normal score
    trace)
!f2py intent(inout) snorm
!f2py intent(inout) zd
!f2py intent(inout) sph_inv
!f2py intent(inout) Us
!f2py intent(inout) pr
!f2py intent(inout) zd2
!f2py intent(inout) yd2
  use quicksort

  implicit none

  character(len=*), intent(in) :: transfile
  integer, intent(in)   :: nvar,ndata
  integer, intent(in)   :: max_iterations
  integer, intent(out)   :: iterations
  real*8, intent(inout) :: snorm(:)
  real*8, intent(inout) :: zd(:,:)
  real*8, intent(inout) :: sph_inv(:,:)
  real*8, intent(inout) :: Us(:,:,:)
  real*8, intent(inout) :: pr(:,:)
  real*8, intent(inout) :: zd2(:,:) !the inverse of the sphering matrix
  real*8, intent(inout) :: yd2(:,:) !the inverse of the sphering matrix
  integer, intent(in)   :: trace         !verbosity

  !locals
  integer lin,len,nvar_tmp,ndata_tmp,iters,iters0,iters1,itertype,i
  character(len=1024) :: str
  real*8, allocatable :: Cs(:,:)
  integer :: test
  !begin

  iterations = 10

  !
  ! Load the transformation file information
  !
  !lin = file_open( transfile, .false., 'stream' )

  write(*,*) 'transfile: ' //  trim(transfile)
  write(*,*) 'nvar: ',nvar
  write(*,*) 'ndata: ',ndata
  write(*,*) 'max_iterations: ',max_iterations

  open (unit=lin, file = transfile, status = 'old', access = 'stream' )

  read( unit=lin,ERR = 97) len
  read( lin) str(1:len)
  if( str(1:len) /= 'PPMT transformation table file' ) then
     write(*,*) trim(str(1:len))
     go to 97
  end if

  read( lin, err = 97 ) nvar_tmp

  if (nvar /= nvar_tmp) then
    write(*,*) 'nvar_tmp: ',nvar_tmp
  end if

  read( lin, err = 97 ) ndata_tmp

  if (ndata /= ndata_tmp) then
    write(*,*) 'ndata_tmp: ',ndata_tmp
  end if


  iters = max_iterations
  allocate(Cs(nvar,nvar), stat=test )

  !Read in the normal score preprocessing info
  read( lin, err = 97 ) snorm
  read( lin, err = 97 ) zd
  !Read in the inverse sphere matrix
  read( lin, err = 97 ) sph_inv

  print *, 'stage 2: iterating:: size(snorm)',size(snorm)


  iters = 0
  iters0 = 0 ! Counter for spatial decorrelation
  iters1 = 0 ! Counter for projection pursuit
  !spatialiter = .false.
  do
      if( iters > 500 )then
          write(*,*)  ' PPMT back-transform only handles 500 iterations for now....'
          stop
      endif
      read( lin, err = 97 ) itertype
      if( itertype == 12345 )then
          iters0 = iters0 + 1 ! Spatial decorrelation
          read( lin, err = 97 ) Cs
      elseif( itertype == 54321 )then
          iters1 = iters1 + 1 ! Projection pursuit
          read( lin, err = 97 ) Us(:,:,iters1)
          read( lin, err = 97 ) pr(:,iters1)
      elseif( itertype == 6789 )then
          exit ! Final marginal normal score
      else
          go to 97 ! Shouldn't get here - error
      endif
      iters = iters + 1
  end do
  read(lin, err = 97 ) zd2
  read(lin, err = 97 ) yd2

  do i = 1,nvar
      call qsortem( yd2(:,i), 1, ndata, zd2(:,i)  )
  enddo


  iterations = iters

  goto 99

97 iterations = 0

99 close(lin)

  deallocate(Cs)
end subroutine

end module
