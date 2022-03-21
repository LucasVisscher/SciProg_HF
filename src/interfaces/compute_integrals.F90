module compute_integrals
 implicit none

 private

 public compute_1e_integrals, compute_2e_integrals, generate_2int, n_ang

 contains

    subroutine compute_1e_integrals (property,ao_basis_bra,ao_basis_ket,ao_integrals,molecule)

     use molecular_structure
     use ao_basis
     use gen1int
     implicit none

     type(basis_set_info_t),intent(in) :: ao_basis_bra, ao_basis_ket
     real(8), intent(out)              :: ao_integrals(ao_basis_bra%nao,ao_basis_ket%nao)
     character(3), intent(in)          :: property
     type(molecular_structure_t), intent(in), optional :: molecule

     type(one_prop_t) prop_operator
     integer info_prop
     integer num_prim_bra, num_prim_ket, num_contr_bra, num_contr_ket, num_gto_bra, num_gto_ket, num_prop
     real(8), allocatable :: contr_ints(:,:,:,:,:)
   
     integer :: ierr, ksh, lsh, koff, loff
     integer :: k, l, kcon, lcon, kang, lang

     ! Initialize for the right kind of property
     num_prop = 1
     select case (property)
     case ("OVL")
       call OnePropCreate(prop_name=INT_OVERLAP,  &
                          one_prop=prop_operator, &
                          info_prop=info_prop)
     case ("KIN")
       call OnePropCreate(prop_name=INT_KIN_ENERGY,  &
                          one_prop=prop_operator, &
                          info_prop=info_prop)
     case ("POT")
       if (present(molecule)) then
       call OnePropCreate(prop_name=INT_POT_ENERGY,  &
                          one_prop=prop_operator, &
                          charge_nuclei=molecule%charge, &
                          coord_nuclei=molecule%coord,   &
                          info_prop=info_prop)
       else
          print*, " Computing potential energy elements requires specification of the molecule !"
          stop " necessary argument missing in compute_1e_integrals"
       end if
      end select
  
      loff = 0
      do lsh = 1, ao_basis_ket%nshells
          num_prim_ket   = size(ao_basis_ket%gtos(lsh)%exponent)
          num_contr_ket  = size(ao_basis_ket%gtos(lsh)%coeff) / num_prim_ket
          num_gto_ket    = n_ang(ao_basis_ket%gtos(lsh)%orb_momentum)
          koff = 0
          do ksh = 1, ao_basis_bra%nshells
             num_prim_bra   = size(ao_basis_bra%gtos(ksh)%exponent)
             num_contr_bra  = size(ao_basis_bra%gtos(ksh)%coeff) / num_prim_bra
             num_gto_bra    = n_ang(ao_basis_bra%gtos(ksh)%orb_momentum)
             ! allocate the contracted integrals to be provided by OnePropGetIntegral
             allocate(contr_ints(num_gto_bra,num_contr_bra, &
                      num_gto_ket,num_contr_ket, &
                      num_prop), stat=ierr)
             call OnePropGetIntegral(idx_bra=1,                                 &
                             coord_bra=ao_basis_bra%gtos(ksh)%coord,            &
                             angular_bra=ao_basis_bra%gtos(ksh)%orb_momentum,   &
                             num_prim_bra=num_prim_bra,                         &
                             exponent_bra=ao_basis_bra%gtos(ksh)%exponent,      &
                             num_contr_bra=num_contr_bra,                       &
                             contr_coef_bra=ao_basis_bra%gtos(ksh)%coeff,       &
                             idx_ket=2,                                         &
                             coord_ket=ao_basis_ket%gtos(lsh)%coord,            &
                             angular_ket=ao_basis_ket%gtos(lsh)%orb_momentum,   &
                             num_prim_ket=num_prim_ket,                         &
                             exponent_ket=ao_basis_ket%gtos(lsh)%exponent,      &
                             num_contr_ket=num_contr_ket,                       &
                             contr_coef_ket=ao_basis_ket%gtos(lsh)%coeff,       &
                             spher_gto=.false., one_prop=prop_operator,         &
                             num_gto_bra=num_gto_bra, num_gto_ket=num_gto_ket,  &
                             num_opt=num_prop, contr_ints=contr_ints)
             ! put these integrals at the right place in the AO matrix
             l = 0
             do lcon = 1, num_contr_ket
               do lang = 1, num_gto_ket
                 l = l + 1
                 k = 0
                 do kcon = 1,num_contr_bra
                   do kang = 1, num_gto_bra
                      k = k + 1
                      ao_integrals(k+koff,l+loff) = contr_ints(kang,kcon,lang,lcon,1)
                  end do
                end do
              end do
             end do
             ! clean up and go to the next block
             deallocate(contr_ints)
             koff = koff + num_contr_bra * num_gto_bra
          end do
          loff = loff + num_contr_ket * num_gto_ket
        end do

   end subroutine compute_1e_integrals

   subroutine generate_2int (ao_basis,ao_integrals)

     ! generate all two-electron integrals for an uncontracted basis

     use ao_basis
     implicit none

     type(basis_set_info_t),intent(in) :: ao_basis
     real(8), intent(out)              :: ao_integrals(:,:,:,:)
     real(8), allocatable              :: ints(:,:,:,:)

     integer :: ish, jsh, ksh, lsh, i, j, k, l, ioff, joff, koff, loff
     integer :: iang, jang, kang, lang, nang_i, nang_j, nang_k, nang_l
     integer :: npri_i, npri_j, npri_k, npri_l

     call interest_initialize()

     loff = 0
     do lsh = 1, ao_basis%nshells
      nang_l = n_ang(ao_basis%gtos(lsh)%orb_momentum)
      npri_l = size(ao_basis%gtos(lsh)%exponent)
      do l = 1, npri_l
       koff = 0
       do ksh = 1, ao_basis%nshells
        nang_k = n_ang(ao_basis%gtos(ksh)%orb_momentum)
        npri_k = size(ao_basis%gtos(ksh)%exponent)
        do k = 1, npri_k
         joff = 0
         do jsh = 1, ao_basis%nshells
          nang_j = n_ang(ao_basis%gtos(jsh)%orb_momentum)
          npri_j = size(ao_basis%gtos(jsh)%exponent)
          do j = 1, npri_j
           ioff = 0
           do ish = 1, ao_basis%nshells
            nang_i = n_ang(ao_basis%gtos(ish)%orb_momentum)
            npri_i = size(ao_basis%gtos(ish)%exponent)
            allocate (ints(nang_i,nang_j,nang_k,nang_l))
            do i = 1, npri_i
               call  compute_2e_integrals(ints, &
                     ao_basis%gtos(ish)%orb_momentum, &
                     ao_basis%gtos(jsh)%orb_momentum, &
                     ao_basis%gtos(ksh)%orb_momentum, &
                     ao_basis%gtos(lsh)%orb_momentum, &
                     ao_basis%gtos(ish)%exponent(i), &
                     ao_basis%gtos(jsh)%exponent(j), &
                     ao_basis%gtos(ksh)%exponent(k), &
                     ao_basis%gtos(lsh)%exponent(l), &
                     ao_basis%gtos(ish)%coord, &
                     ao_basis%gtos(jsh)%coord, &
                     ao_basis%gtos(ksh)%coord, &
                     ao_basis%gtos(lsh)%coord)
               ! put these integrals at the right place in the AO matrix to be consistent with gen1int
               ! we use offsets that indicate the position of a certain shell in the full list of functions
               do lang = 1, nang_l
                do kang = 1, nang_k
                 do jang = 1, nang_j
                  do iang = 1, nang_i
                     ao_integrals(ioff+iang,joff+jang,koff+kang,loff+lang) = ints(iang,jang,kang,lang)
                  end do
                 end do
                end do
               end do

             ioff = ioff + nang_i
            end do ! i
            deallocate (ints)
           end do ! ish
           joff = joff + nang_j
          end do ! j
         end do ! jsh
         koff = koff + nang_k
        end do ! k
       end do ! ksh
       loff = loff + nang_l
      end do ! l 
     end do ! lsh
     
   end subroutine

   subroutine compute_2e_integrals (ints, &
                                    angular_bra_1,angular_ket_1,angular_bra_2,angular_ket_2, &
                                    exponent_bra_1,exponent_ket_1,exponent_bra_2,exponent_ket_2, &
                                    coord_bra_1,coord_ket_1,coord_bra_2,coord_ket_2)
     use module_interest_eri
     real(8), intent(out) :: ints(:,:,:,:)
     integer, intent(in)  :: angular_bra_1,angular_ket_1,angular_bra_2,angular_ket_2
     real(8), intent(in)  :: exponent_bra_1,exponent_ket_1,exponent_bra_2,exponent_ket_2
     real(8), dimension(3), intent(in) :: coord_bra_1,coord_ket_1,coord_bra_2,coord_ket_2

     ! local variables used to convert gen1int-style arguments to interest-style arguments
     real(8) :: ci, cj, ck, cl, fijkl
     real(8) :: ei, ej, ek, el
     real(8) :: xi, yi, zi, xj, yj, zj, xk, yk, zk, xl, yl, zl
     integer :: li, lj, lk, ll, ni, nj, nk, nl
     real(8), allocatable :: gout(:)

     ! copy and restructure the data to fit interest
     ni = size(ints,1)
     nj = size(ints,2)
     nk = size(ints,3)
     nl = size(ints,4)
     allocate (gout(ni*nj*nk*nl))
     fijkl = 1.0D0
     ci = 1.D0
     cj = 1.D0
     ck = 1.D0
     cl = 1.D0
     li = angular_bra_1 + 1
     lj = angular_ket_1 + 1
     lk = angular_bra_2 + 1
     ll = angular_ket_2 + 1
     ei = exponent_bra_1
     ej = exponent_ket_1
     ek = exponent_bra_2
     el = exponent_ket_2
     xi = coord_bra_1(1)
     yi = coord_bra_1(2)
     zi = coord_bra_1(3)
     xj = coord_ket_1(1)
     yj = coord_ket_1(2)
     zj = coord_ket_1(3)
     xk = coord_bra_2(1)
     yk = coord_bra_2(2)
     zk = coord_bra_2(3)
     xl = coord_ket_2(1)
     yl = coord_ket_2(2)
     zl = coord_ket_2(3)
     call interest_eri('llll',fijkl,gout,&
                       lk,ek,xk,yk,zk,ck,&
                       ll,el,xl,yl,zl,cl,&
                       li,ei,xi,yi,zi,ci,&
                       lj,ej,xj,yj,zj,cj )
     ! copy integrals into the 4-dimensional output array
     ints = reshape(gout(1:ni*nj*nk*nl),(/ni,nj,nk,nl/))
     deallocate(gout)

   end subroutine

   integer function n_ang(angular)
     ! function to give the number of gtos in a block (s:1 function, p:3 functions, d:6 functions, etc.)
     integer, intent(in) :: angular
     n_ang = (angular+1)*(angular+2)/2
   end function


end module
