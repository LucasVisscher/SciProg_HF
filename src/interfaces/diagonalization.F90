module Diagonalization

implicit none

private

public :: diagonalize, solve_genev

contains

   subroutine solve_genev (matrix,metric,eigenvectors,eigenvalues)

   ! Solve generalized eigenvalue problem using Lowdins transformation to orthonormal basis

   real*8, intent(in)               :: matrix(:,:)
   real*8, intent(in)               :: metric(:,:)
   real*8, intent(out)              :: eigenvalues(:)
   real*8, intent(out)              :: eigenvectors(:,:)
   real*8, allocatable              :: u(:,:),v(:,:),d(:),o(:,:)
   integer i, n ! counter, matrix dimension

   n = size(matrix,1)
   ! Error checking
   if (size(metric,1) /= n .or. size(metric,2) /= n) then
      print*," inconsistent metric for generalize_diagonalize"
      stop "error in diagonalize routine"
   end if

   allocate (u(n,n))  ! a scratch matrix used for various purposes
   allocate (v(n,n))  ! the Lowdin transformation to the orthonormal basis
   allocate (o(n,n))  ! the matrix transformed to the orthonormal basis
   allocate (d(n))    ! eigenvalues of the metric

   ! Diagonalize the metric (typically this is the overlap) matrix
   call diagonalize (metric,u,d)

   ! Form the Lowdin transformation matrix v = u d^{-1/2} u^t
   ! Note that for simplicty we do NOT check for small eigenvalues
   ! In real QC codes one reduces the basis set size when these are encountered (linear dependency)
   forall (i=1:n) u(:,i) = u(:,i) * d(i)**(-0.25D0)
   ! With the scaling above we have v = u u^t which can be obtained by a single matrix multiplication
   call dgemm ('n','t',n,n,n,1.D0,u,n,u,n,0.D0,v,n)
   
   ! transform matrix to the orthogonal basis: o = v^t matrix v
   call dgemm ('n','n',n,n,n,1.D0,matrix,n,v,n,0.D0,u,n)
   call dgemm ('t','n',n,n,n,1.D0,v,n,u,n,0.D0,o,n)

   ! diagonalize the matrix that is now a standard eigenvalue problem, re-using u again to store the eigenvectors
   call diagonalize (o,u,eigenvalues)

   ! transform the eigenvectors to the original non-orthogonal basis
   call dgemm ('n','n',n,n,n,1.D0,v,n,u,n,0.D0,eigenvectors,n)

   deallocate (u,v,o,d)

   end

   subroutine diagonalize (matrix,eigenvectors,eigenvalues)

   real*8, intent(in)               :: matrix(:,:)
   real*8, intent(out)              :: eigenvalues(:)
   real*8, intent(out), optional    :: eigenvectors(:,:)
   real*8, allocatable              :: a(:,:),work(:)
   character, parameter             :: uplo='L'
   character                        :: jobz

   integer                          :: n, lda, info, lwork
   integer                          :: ILAENV

!  Determine the calculation mode
   if (present(eigenvectors)) then
      jobz = 'V'
   else
      jobz = 'N'
   end if

!  Error checking
   if (size(matrix,1) /= size(matrix,2)) then
      print*," diagonalize assumes square matrices"
      stop "error in diagonalize routine"
   end if
   if (size(matrix,1) > size(eigenvalues)) then
      print*," dimension of eigenvalue array too small "
      stop "error in diagonalize routine"
   end if
   if (present(eigenvectors)) then
      if (size(eigenvectors,1) /= size(eigenvectors,2)) then
         print*," diagonalize assumes square matrices"
         stop "error in diagonalize routine"
      end if
      if (size(matrix,1) /= size(eigenvectors,1)) then
         print*," diagonalize assumes same size array for matrix and eigenvectors "
         stop "error in diagonalize routine"
      end if
   end if

!  Initialize and diagonalize using lapack's dsyev routine
   n = size(matrix,1)
   lda = n
   lwork = n * ( 2 + ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 ))
   allocate (work(lwork))
   allocate (a(n,n))
   a = matrix

   call DSYEV( JOBZ, UPLO, N, A, LDA, eigenvalues, WORK, LWORK, INFO )

   if (present(eigenvectors)) eigenvectors = a

   end subroutine

end module
