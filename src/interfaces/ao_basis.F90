       module ao_basis

        implicit none
        private

        integer, parameter, private:: REALD=8 ! default real for this module

        public add_shell_to_basis, clear_basis 

!       Basis function info:
        type, public:: basis_func_info_t
         integer:: orb_momentum=-1       ! orbital momentum: (0,1,2,3,...)
         integer:: atom_number=-1        ! atom sequential number [1..MAX] (if centered on a nucleus of some atom)
         integer:: atom_element=0        ! atomic element number (if centered on a nucleus of some atom)
         integer:: n_primitives=0        ! number of primitives in the contraction
         integer:: n_contracted=0        ! number of contracted that are made
         real(REALD):: coord(1:3)        ! coordinate of the basis function center (typically a nucleus)
         real(REALD), allocatable :: exponent(:) ! exponents (size is n_primitives)
         real(REALD), allocatable :: coeff(:) !contraction coefficients (size is n_primitives * n_contracted)
        end type basis_func_info_t
!       Basis set info:
        type, public:: basis_set_info_t
         integer:: nshells=0             ! number of shells in the basis set
         integer:: nao=0                 ! number of basis functions
         integer:: basis_angular=1       ! 1=cartesian, 2=spherical
         type(basis_func_info_t), pointer:: gtos(:)
        end type basis_set_info_t

        contains

        subroutine clear_gto(gto)
          type(basis_func_info_t)  :: gto
          if (allocated(gto%exponent)) deallocate(gto%exponent)
          if (allocated(gto%coeff))    deallocate(gto%coeff)
        end subroutine

        subroutine clear_basis(ao_basis)
          type(basis_set_info_t)  :: ao_basis
          integer                 :: ishell
          do ishell = 1, ao_basis%nshells
             call clear_gto(ao_basis%gtos(ishell))
          enddo
          ao_basis%nshells = 0
          ao_basis%nao     = 0
        end subroutine

        subroutine add_shell_to_basis(ao_basis,angular,coord,alpha,exponents,coefficients)
          type(basis_set_info_t)               :: ao_basis
          type(basis_func_info_t), allocatable :: gtos(:)
          integer, intent(in)                  :: angular
          real(8), intent(in)                  :: coord(3)
          real(8), intent(in), optional        :: alpha ! used for uncontracted basis sets
          real(8), intent(in), optional        :: exponents(:),coefficients(:,:) ! for contracted sets

          integer  :: i, ishell, nshells, nao, nprim, ncont
          real(8), allocatable :: exponent(:), coeff(:)
   
          if (present(alpha)) then
             allocate (exponent(1))
             exponent(1) = alpha
          elseif (present(exponents)) then
             allocate (exponent(size(exponents)))
             exponent = exponents
          else
             stop "no exponent(s) given in add_shell_to_basis"
          end if

          ! first make a copy of the original basis
          nshells = ao_basis%nshells
          nao     = ao_basis%nao
          if (allocated(gtos)) deallocate(gtos)
          allocate(gtos(nshells+1))
          do ishell = 1, nshells
             nprim = ao_basis%gtos(ishell)%n_primitives
             ncont = ao_basis%gtos(ishell)%n_contracted
             allocate(gtos(ishell)%exponent(nprim))
             allocate(gtos(ishell)%coeff(nprim*ncont))
             gtos(ishell)= ao_basis%gtos(ishell)
          end do

          ! prepare the new shell (define coefficients in right format, normalize the functions)
          nprim = size(exponent)
          if (present(coefficients)) then
             ! contracted block of functions, flatten the coefficient matrix to 1D
             if (size(coeff,1) == nprim) then
                ncont = size(coefficients,2)
                allocate (coeff(ncont*nprim))
                coeff = reshape(coefficients,(/nprim*ncont/))
             else
                stop "inconsistency in contraction coefficients detected"
             endif
          else
            ncont = nprim
            ! uncontracted block of functions, the coefficient matrix is a unit matrix
            allocate (coeff(ncont*nprim))
            coeff = 0.D0
            forall (i=1:ncont) coeff(nprim*(i-1)+i) = 1.D0
          end if

          ! copy all of this into the gtos array
          gtos(nshells+1)%orb_momentum = angular
          gtos(nshells+1)%n_primitives = nprim
          gtos(nshells+1)%n_contracted = ncont
          gtos(nshells+1)%coord = coord
          allocate(gtos(nshells+1)%exponent(nprim))
          allocate(gtos(nshells+1)%coeff(nprim*ncont))
          gtos(nshells+1)%exponent = exponent
          gtos(nshells+1)%coeff    = coeff

          ! copy into the basis
          call clear_basis(ao_basis)
          allocate(ao_basis%gtos(nshells+1))
          do ishell = 1, nshells+1
             nprim = gtos(ishell)%n_primitives
             ncont = gtos(ishell)%n_contracted
             allocate(ao_basis%gtos(ishell)%exponent(nprim))
             allocate(ao_basis%gtos(ishell)%coeff(nprim*ncont))
             ao_basis%gtos(ishell)= gtos(ishell)
          end do
          ao_basis%nshells = nshells+1
          ao_basis%nao     = nao + n_ang(angular) * ncont
    
        end subroutine

        integer function n_ang(angular)
           ! function to give the number of gtos in a block (s:1 function, p:3 functions, d:6 functions, etc.)
           integer, intent(in) :: angular
           n_ang = (angular+1)*(angular+2)/2
         end function

       end module ao_basis

