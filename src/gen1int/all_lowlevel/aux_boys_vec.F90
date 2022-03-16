!!  gen1int: compute one-electron integrals using rotational London atomic-orbitals
!!  Copyright 2009-2012 Bin Gao, and Andreas Thorvaldsen
!! 
!!  gen1int is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU Lesser General Public License as published by
!!  the Free Software Foundation, either version 3 of the License, or
!!  (at your option) any later version.
!! 
!!  gen1int is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!!  GNU Lesser General Public License for more details.
!! 
!!  You should have received a copy of the GNU Lesser General Public License
!!  along with gen1int. If not, see <http://www.gnu.org/licenses/>.
!!
!!  This file computes the Boys function.
!!
!!  2009-06-27, Bin Gao:
!!  * first version

#include "stdout.h"

  !> \brief computes the Boys function for a given argument up to some order
  !> \detail computes the Boys function using (1) Taylor expansion and downward
  !>         recurrence relation for small argument, and (2) modified asymptotic
  !>         series and upward recurrence relation for large argument
  !> \author Bin Gao
  !> \date 2009-06-27
  !> \param min_order_boys is the minimum order of Boys functions
  !> \param max_order_boys is the maximum order of Boys functions
  !> \param arg_boys is the given argument of Boys functions
  !> \return val_boys contains the values of Boys functions
  !> \note For Boys function \f$F_n(T)\f$, according to our tests, the current asymptotic
  !>       series expansion \fn(dboys_asymp) can achieve \f$10^{-12}\f$ relative
  !>       accuracy after \f$T>25\f$ (But for quite smaller \f$n\f$ if \f$T\f$ is not large enough),
  !>       we may need to provide a new asymptotic series expansion later ...;
  !>       the modified asymptotic series expansion DO NOT have a stable upward
  !>       recurrence relations when \f$n\ge2.7T-7.7\f$ (\f$T>12\f$) according to
  !>       our tests (negtive values may occur); but for the time being, it seems to be enough;
  !>       last but not least, for the power series expansion, we have used downward recurrence
  !>       relation, which seems to be very stable according to our tests (\f$0\le T\le 400\f$,
  !>       \f$n=10,1000\f$), even the value of largest order is not correct (for \f$T<24\f$
  !>       it accurate enough for \f$n=1000\f$; and it is accurate enough for \f$T<28\f$, \f$n=100\f$;
  !>       according to our tests, for larger \f$T\f$, it may become reasonable after \f$n<2T\f$)
  subroutine aux_boys_vec(min_order_boys, max_order_boys, arg_boys, val_boys)
    use xkind
    implicit none
    integer, intent(in) :: min_order_boys
    integer, intent(in) :: max_order_boys
    real(REALK), intent(in) :: arg_boys
    real(REALK), intent(out) :: val_boys(min_order_boys:max_order_boys)
!f2py intent(in) :: min_order_boys
!f2py intent(in) :: max_order_boys
!f2py intent(in) :: arg_boys
!f2py intent(out) :: val_boys
!f2py depend(min_order_boys) :: val_boys
!f2py depend(max_order_boys) :: val_boys
#include "private/pi.h"
#include "private/boys_power.h"
#include "private/max_gen_order.h"
#include "private/tab_boys.h"
    logical gen_tab_boys                            !if generating the tabulated Boys function during runtime
    integer max_order_tab                           !maximum order of the tabulated Boys function
    real(REALK), allocatable :: tab_boys(:,:)       !tabulated Boys function during runtime
    integer, parameter :: NUM_TAYLOR_TERM = 6       !number of terms of Taylor expansion for
                                                    !the tabulated Boys functions
    logical use_asym                                !indicates using the four term modified asymptotic series
    real(REALK), parameter :: COEF_ASYM(0:3) =   &  !four term modified asymptotic series
      (/0.4999489092_REALK, -0.2473631686_REALK, &
        0.321180909_REALK, -0.3811559346_REALK/)
    integer iarg              !incremental recorder for the argument
    integer iord              !incremental recorder for the orders
    real(REALK) curr_arg      !\f$T^*-T\f$
    real(REALK) twice_arg     !twice of the argument
    real(REALK) max_val_boys  !\f$F_J(T)\f$
    real(REALK) curr_exp_arg  !\f$exp(-T)\f$
    real(REALK) curr_iarg     !\f$1/T\f$
    real(REALK) half_iarg     !\f$1/(2T)\f$
    real(REALK) div_down      !divisor in the downward recurrence relations
    real(REALK), allocatable :: val_low(:)
                              !lower orders Boys functions (temporarily used)
    integer ierr              !error information
#if defined(XTIME)
    real(REALK) curr_time     !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! checks if the order of pre-tabulated Boys functions is enough
    gen_tab_boys = NUM_TAYLOR_TERM+max_order_boys>MAX_GEN_ORDER
    ! \f$0\le T<12\f$, using Taylor expansion and downward recurrence relation
    if (arg_boys<=MAX_ARG_TAB) then
      ! the nearest \f$T^*\f$, notice that the first argument of tabulated Boys function is 0!
      iarg = nint(arg_boys/INTERV_TAB)
      ! computes \f$T^*-T\f$
      curr_arg = real(iarg,REALK)*INTERV_TAB-arg_boys
      ! \f$\frac{(T*-T)^k}{k!}=1\f$, for \f$k=0\f$
      curr_exp_arg = 1.0_REALK
      ! generates the tabulated Boys functions during runtime
      if (gen_tab_boys) then
        ! sets the maximum order of tabulated Boys functions
        max_order_tab = NUM_TAYLOR_TERM+max_order_boys
        allocate(tab_boys(max_order_boys:max_order_tab,0:NSTEPS_TAB), stat=ierr)
        if (ierr/=0)                                                     &
          call error_stop("aux_boys_vec", "failed to allocate tab_boys", &
                          (max_order_tab-max_order_boys+1)*(NSTEPS_TAB+1))
        ! computes the tabulated Boys functions
        call boys_power(max_order_boys, max_order_tab, MIN_ARG_TAB, &
                        INTERV_TAB, NSTEPS_TAB, tab_boys)
        ! the first Taylor's term
        val_boys(max_order_boys) = tab_boys(max_order_boys,iarg)
        ! the left Taylor's terms
        do iord = 1, NUM_TAYLOR_TERM
          curr_exp_arg = curr_exp_arg*curr_arg/real(iord,REALK)
          val_boys(max_order_boys) = val_boys(max_order_boys) &
            + curr_exp_arg*tab_boys(iord+max_order_boys,iarg)
        end do
        deallocate(tab_boys)
      ! uses the pre-tabulated Boys functions
      else
        ! the first Taylor's term
        val_boys(max_order_boys) = PRE_TAB_BOYS(iarg*(MAX_GEN_ORDER+1)+max_order_boys+1)
        ! the left Taylor's terms
        do iord = 1, NUM_TAYLOR_TERM
          curr_exp_arg = curr_exp_arg*curr_arg/real(iord,REALK)
          val_boys(max_order_boys) = val_boys(max_order_boys) &
            + curr_exp_arg*PRE_TAB_BOYS(iarg*(MAX_GEN_ORDER+1)+iord+max_order_boys+1)
        end do
      end if
      curr_exp_arg = exp(-arg_boys)
      twice_arg = arg_boys+arg_boys
      div_down = real(max_order_boys+max_order_boys+1,REALK)
      ! using downward recurrence relations for other \f$F_j(T)\f$
      do iord = max_order_boys-1, min_order_boys, -1
        div_down = div_down-2.0_REALK
        val_boys(iord) = (twice_arg*val_boys(iord+1)+curr_exp_arg)/div_down
      end do
    ! \f$T>12\f$, using modified asymptotic series expansion and upward recurrence relation
    else
!FIXME: this modified asymptotic series expansion seems to be bad when
!FIXME: \f$n>4T/3+16\f$ (\f$T>12\f$)
!     radovan: the next expression can overflow the integer
!     if (max_order_boys>=4*floor(arg_boys/3.0_REALK)+16) then
!     more robust:
      if ((max_order_boys-16) >= 4.0_REALK*arg_boys/3.0_REALK) then
        write(STDOUT,100) "maximum order of Boys functions:", max_order_boys
        write(STDOUT,110) "argument of Boys functions:", arg_boys
        write(STDOUT,100) "the upward recurrence relation is unstable!"
        call error_stop("aux_boys_vec", "please write to authors", max_order_boys)
      end if
      ! allocates memory for the lower orders Boys functions (temporarily used)
      allocate(val_low(0:min_order_boys), stat=ierr)
      if (ierr/=0)                                                    &
        call error_stop("aux_boys_vec", "failed to allocate val_low", &
                        min_order_boys+1)
      ! for \f$T\le2J+36\f$, we use four asymptotic terms, and exact upward recurrence relation
      if (gen_tab_boys) then
        use_asym = arg_boys<=real(2*max_order_boys+36,REALK)
      else
        use_asym = arg_boys<=real(2*(MAX_GEN_ORDER-NUM_TAYLOR_TERM)+36,REALK)
      end if
      if (use_asym) then
        curr_iarg = 1.0_REALK/arg_boys
        curr_exp_arg = exp(-arg_boys)
        ! \f$g\f$
        max_val_boys = COEF_ASYM(0)+curr_iarg*(COEF_ASYM(1) &
                     + curr_iarg*(COEF_ASYM(2)+curr_iarg*COEF_ASYM(3)))
        ! \f$F_0(T)=\sqrt{\pi/T}/2-\exp(-T)*g/T\f$
        val_low(0) = sqrt(PI*curr_iarg)*0.5_REALK &
                   - curr_exp_arg*max_val_boys*curr_iarg
        ! \f$1/(2T)\f$
        half_iarg = curr_iarg*0.5_REALK
        ! \f$\exp(-T)/(2T)\f$
        curr_exp_arg = half_iarg*curr_exp_arg
        ! upward recurrence relation for other \f$F_j(T)\f$
        do iord = 1, min_order_boys
          val_low(iord) = half_iarg*val_low(iord-1)-curr_exp_arg
          half_iarg = curr_iarg+half_iarg
        end do
        val_boys(min_order_boys) = val_low(min_order_boys)
        do iord = min_order_boys+1, max_order_boys
          val_boys(iord) = half_iarg*val_boys(iord-1)-curr_exp_arg
          half_iarg = curr_iarg+half_iarg
        end do
      ! for \f$T>2J+36\f$, we use |f$F_0(T)=\sqrt(\pi/T)/2\f$, and approximate
      ! upward recurrence relation \f$F_{j+1}(T)=(2T)^{-1}(2j+1)F_j(T)\f$
      else
        curr_iarg = 1.0_REALK/arg_boys
        val_low(0) = sqrt(PI*curr_iarg)*0.5_REALK
        ! \f$1/(2T)\f$
        half_iarg = curr_iarg*0.5_REALK
        ! upward approximate recurrence relation for other \f$F_j(T)\f$
        do iord = 1, min_order_boys
          val_low(iord) = half_iarg*val_low(iord-1)
          half_iarg = curr_iarg+half_iarg
        end do
        val_boys(min_order_boys) = val_low(min_order_boys)
        do iord = min_order_boys+1, max_order_boys
          val_boys(iord) = half_iarg*val_boys(iord-1)
          half_iarg = curr_iarg+half_iarg
        end do
      end if
      deallocate(val_low)
    end if
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "aux_boys_vec", STDOUT)
#endif
    return
100 format('aux_boys_vec>> ',A,I6)
110 format('aux_boys_vec>> ',A,F14.8)
  end subroutine aux_boys_vec

  !> \brief power series expansion of Boys function (for small argument)
  !> \detail the power series expansion can be found, for example, in:
  !>         V. R. Saunders. An introduction to molecular integral evaluation.
  !>         In G.H.F. Diercksen, B.T. Sutcliffe, and A. Veillard, editors,
  !>         Computational Techniques in Quantum Chemistry and Molecular Physics,
  !>         Eq. (39), page 347, 1975.
  !> \author Bin Gao
  !> \date 2009-06-27
  !> \param strt_order is the start order of Boys function
  !> \param end_order is the end order of Boys function
  !> \param strt_arg is the start argument of Boys function
  !> \param step_arg is the step of the argument
  !> \param num_arg is the number of arguments
  !> \return val_boys contains the values of Boys function
  subroutine boys_power(strt_order, end_order, strt_arg, step_arg, num_arg, val_boys)
    use xkind
    implicit none
    integer, intent(in) :: strt_order
    integer, intent(in) :: end_order
    real(REALK), intent(in) :: strt_arg
    real(REALK), intent(in) :: step_arg
    integer, intent(in) :: num_arg
    real(REALK), intent(out) :: val_boys(strt_order:end_order,0:num_arg)
    real(REALK), parameter :: MAX_ARG_POW = 30.0_REALK
                              !maximum argument for power series expansion
!FIXME how to determine?? which depends on the order & argument
!FIXME 200 terms are enough for arguments smaller than, like 20.0
    integer, parameter :: MAX_NTERM_POW = 200
                              !maximum number of terms for the power series expansion
    real(REALK), parameter :: CUT_OFF_POW = 10.0_REALK**(-18)
                              !cut off for the terms in power series expansion
    real(REALK) curr_arg      !argument of Boys function at each step
    real(REALK) twice_arg     !twice of the argument
    integer iord              !incremental recorder for the orders
    integer iarg              !incremental recorder for the arguments
    real(REALK) val_power     !value of power series expansion's term
    real(REALK) div_power     !divisor in the term of power series expansion
    logical not_converged     !indicates the power series expansion is not converged
    integer ipower            !incremental recorder for the summation
    real(REALK) curr_exp_arg  !\f$exp(-T)\f$ with \f$T\f$ being the argument at each step
#if defined(XTIME)
    real(REALK) curr_time     !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! the minimum argument
    curr_arg = min(strt_arg+step_arg*real(num_arg,REALK), strt_arg)
    ! the maximum argument
    twice_arg = max(strt_arg+step_arg*real(num_arg,REALK), strt_arg)
    ! checks if the value of largest order is reasonable
    if (twice_arg>23.0_REALK .and. end_order>=nint(curr_arg+curr_arg)) then
      write(STDOUT,100) "minimum argument of Boys function:", curr_arg
      write(STDOUT,110) "maximum order of Boys function:", end_order
      write(STDOUT,100) "warning! values of the largest order may not be accurate or correct!"
    end if
    ! loops over arguments
    do iarg = 0, num_arg
      ! for accuracy, this is much better than adding \var(step_arg) to \var(curr_arg) at each time
      curr_arg = strt_arg+step_arg*real(iarg,REALK)
      ! if the argument exceeds \var(MAX_ARG_POW), we change to asymptotic series expansion
      if (curr_arg>MAX_ARG_POW) then
        write(STDOUT,100) "argument of Boys function:", curr_arg
        write(STDOUT,100) "it is too large, we change to asymptotic series expansion ..."
        call boys_asymp(strt_order, end_order, curr_arg, step_arg, 0, val_boys(:,iarg))
      else
        twice_arg = curr_arg+curr_arg
        ! we use power series expansion for the largest order,
        ! and downward recurrence relations for others
        div_power = real(end_order+end_order+1,REALK)
        val_power = 1.0_REALK/div_power
        ! the first term
        val_boys(end_order,iarg) = val_power
        ! loops over power series expansion
        not_converged = .true.
        do ipower = 1, MAX_NTERM_POW
          div_power = div_power+2.0_REALK
          val_power = val_power*twice_arg/div_power
          val_boys(end_order,iarg) = val_boys(end_order,iarg)+val_power
          if (val_power<=CUT_OFF_POW) then
            not_converged = .false.
            exit
          end if
        end do
        if (not_converged) then
          write(STDOUT,110) "power series expansion is not converged after step:", MAX_NTERM_POW
          write(STDOUT,100) "with argument:", curr_arg
          write(STDOUT,110) "and order:", end_order
          call error_stop("boys_power",                                        &
                          "try to increase MAX_NTERM_POW in aux_boys_vec.F90", &
                          end_order)
        end if
        curr_exp_arg = exp(-curr_arg)
        val_boys(end_order,iarg) = curr_exp_arg*val_boys(end_order,iarg)
        ! using downward recurrence relations for other orders
        div_power = real(end_order+end_order+1,REALK)
        do iord = end_order-1, strt_order, -1
          div_power = div_power-2.0_REALK
          val_boys(iord,iarg) = (twice_arg*val_boys(iord+1,iarg)+curr_exp_arg)/div_power
        end do
      end if
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "boys_power", STDOUT)
#endif
    return
100 format("boys_power>> ",A,F14.8)
110 format("boys_power>> ",A,I6)
  end subroutine boys_power

  !> \brief asymptotic series expansion of Boys funciton (for large argument!!)
  !> \detail We use \f$F_n(T) = \frac{(2n-1)!!}{2^{n+1}}\sqrt{\frac{\pi}{T^{2n+1}}}\f$,
  !>         see, for example, Trygve Helgaker, Poul Jorgensen, Jeppe Olsen,
  !>         Molecular Electronic Structure Theory, Eq. (9.8.9), page 365.
  !> \author Bin Gao
  !> \date 2009-06-27
  !> \param strt_order is the start order of Boys function
  !> \param end_order is the end order of Boys function
  !> \param strt_arg is the start argument of Boys function
  !> \param step_arg is the step of the argument
  !> \param num_arg is the number of arguments
  !> \return val_boys contains the values of Boys function
  subroutine boys_asymp(strt_order, end_order, strt_arg, step_arg, num_arg, val_boys)
    use xkind
    implicit none
    integer, intent(in) :: strt_order
    integer, intent(in) :: end_order
    real(REALK), intent(in) :: strt_arg
    real(REALK), intent(in) :: step_arg
    integer, intent(in) :: num_arg
    real(REALK), intent(out) :: val_boys(strt_order:end_order,0:num_arg)
#include "private/pi.h"
    real(REALK), parameter :: MIN_ARG_ASYM = 30.0_REALK
                           !minimum argument for asymptotic series expansion
    real(REALK) curr_arg   !argument of Boys function at each step
    real(REALK) curr_iarg  !\f$\frac{1}{T}\f$ with \f$T\f$ being the current argument
    integer iord           !incremental recorder for the orders
    integer iarg           !incremental recorder for the arguments
    real(REALK) div_asymp  !dividend in the asymptotic term
    integer iasymp         !incremental recorder for the asymptotic series expansion
#if defined(XTIME)
    real(REALK) curr_time  !current CPU time
    ! sets current CPU time
    call xtimer_set(curr_time)
#endif
    ! the minimum argument
    curr_arg = min(strt_arg+step_arg*real(num_arg,REALK), strt_arg)
    ! checks if the minimum argument is large enough
    if (curr_arg<=MIN_ARG_ASYM) then
      write(STDOUT,100) "minimum argument of Boys function:", curr_arg
      write(STDOUT,100) "minimum argument for asymptotic series expansion:", MIN_ARG_ASYM
      call error_stop("boys_asymp", &
                      "too small argument for asymptotic series expansion", -1)
    end if
    ! loops over number of arguments
    do iarg = 0, num_arg
      curr_arg = strt_arg+step_arg*real(iarg,REALK)
      curr_iarg = 1.0_REALK/curr_arg
      div_asymp = -curr_iarg*0.5_REALK
      ! \f$\sqrt{\frac{\pi}{T}}/2\f$
      val_boys(strt_order,iarg) = sqrt(PI/curr_arg)*0.5_REALK
      ! the smallest order
      do iasymp = 1, strt_order
        ! adds \f$\frac{1}{T}\f$
        div_asymp = div_asymp+curr_iarg
        ! multiplied by \f$\frac{2n-1}{2T}\f$
        val_boys(strt_order,iarg) = val_boys(strt_order,iarg)*div_asymp
      end do
      ! the left orders
      do iord = strt_order, end_order-1
        ! adds \f$\frac{1}{T}\f$
        div_asymp = div_asymp+curr_iarg
        val_boys(iord+1,iarg) = val_boys(iord,iarg)*div_asymp
      end do
    end do
#if defined(XTIME)
    ! prints the CPU elapsed time
    call xtimer_view(curr_time, "boys_asymp", STDOUT)
#endif
    return
100 format("boys_asymp>> ",A,F14.8)
  end subroutine boys_asymp
