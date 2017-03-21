!< Define the Riemann solver for ideal compressible fluid for FORESEER library.

module foreseer_riemann_solver_compressible_object
!< Define the Riemann solver for ideal compressible fluid for FORESEER library.

use foreseer_conservative_compressible, only : conservative_compressible
use foreseer_conservative_object, only : conservative_object
use foreseer_eos_compressible, only : eos_compressible
use foreseer_eos_object, only : eos_object
use foreseer_riemann_solver_object, only : riemann_solver_object
use penf, only : R8P, str, ZeroR8
use vecfor, only : vector

implicit none
private
public :: riemann_solver_compressible_object

type, extends(riemann_solver_object), abstract :: riemann_solver_compressible_object
   !< Riemann solver for ideal compressible fluid object class.
   !<
   !< The ideal compressible fluid generates a 3-waves pattern: 2 genuinely non-linear acoustic waves and 1
   !< linear-degener contact discontinuity.
   !<
   !<```
   !<   t ^                                .
   !<     |     S1 _       S2 _            .      _ S=u23   _ S3       _ S4
   !<     |       |\_        |\_           .      /|      __/|       __/|
   !<     |          \__        \_     U2  .     /  U3  _/       ___/
   !<     |             \___      \_       .    /    __/     ___/
   !<     |                 \____   \_     .   /   _/    ___/
   !<     |                      \___ \_   .  / __/  ___/
   !<     |          UL=U1           \__\_ . /_/____/          UR=U4
   !<     |                              \\.///
   !<  ---+--------------------------------o--------------------------------->
   !<     |                                xo                                x
   !<```
   type(eos_compressible) :: eos_1      !< Equation of state 1.
   real(R8P)              :: r_1=0._R8P !< Density of state 1.
   real(R8P)              :: u_1=0._R8P !< Velocity (normal) of state 1.
   real(R8P)              :: p_1=0._R8P !< Pressure of state 1.
   real(R8P)              :: a_1=0._R8P !< Speed of sound of state 1.
   type(eos_compressible) :: eos_4      !< Equation of state 4.
   real(R8P)              :: r_4=0._R8P !< Density of state 4.
   real(R8P)              :: u_4=0._R8P !< Velocity (normal) of state 4.
   real(R8P)              :: p_4=0._R8P !< Pressure of state 4.
   real(R8P)              :: a_4=0._R8P !< Speed of sound of state 4.
   real(R8P)              :: u23=0._R8P !< Velocity (normal) of intermediate states.
   real(R8P)              :: p23=0._R8P !< Pressure of intermediate states.
   real(R8P)              :: r_2=0._R8P !< Density of state 2.
   real(R8P)              :: a_2=0._R8P !< Speed of sound of state 2.
   real(R8P)              :: r_3=0._R8P !< Density of state 3.
   real(R8P)              :: a_3=0._R8P !< Speed of sound of state 3.
   real(R8P)              :: s_1=0._R8P !< Left-front of left wave.
   real(R8P)              :: s_2=0._R8P !< Right-front of left wave.
   real(R8P)              :: s_3=0._R8P !< Left-front of right wave.
   real(R8P)              :: s_4=0._R8P !< Right-front of right wave.
   contains
      ! deferred methods
      procedure, pass(self) :: description !< Return pretty-printed object description.
      ! public methods
      procedure, pass(self) :: compute_fluxes            !< Compute fluxes at interface `x=xo`.
      procedure, pass(self) :: compute_states23_from_u23 !< Compute interstates 2 and 3 given veloctiy `S=u23`.
      procedure, pass(self) :: set_states14              !< Set states 1 and 4.
endtype riemann_solver_compressible_object

contains
   ! deferred methods
   pure function description(self, prefix) result(desc)
   !< Return a pretty-formatted object description.
   class(riemann_solver_compressible_object), intent(in)           :: self             !< Solver.
   character(*),                              intent(in), optional :: prefix           !< Prefixing string.
   character(len=:), allocatable                                   :: prefix_          !< Prefixing string, local variable.
   character(len=:), allocatable                                   :: desc             !< Description.
   character(len=1), parameter                                     :: NL=new_line('a') !< New line character.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   desc = ''
   desc = desc//prefix_//'r_1 = '//trim(str(n=self%r_1))//NL
   desc = desc//prefix_//'u_1 = '//trim(str(n=self%u_1))//NL
   desc = desc//prefix_//'p_1 = '//trim(str(n=self%p_1))//NL
   desc = desc//prefix_//'a_1 = '//trim(str(n=self%a_1))//NL
   desc = desc//prefix_//'r_4 = '//trim(str(n=self%r_4))//NL
   desc = desc//prefix_//'u_4 = '//trim(str(n=self%u_4))//NL
   desc = desc//prefix_//'p_4 = '//trim(str(n=self%p_4))//NL
   desc = desc//prefix_//'a_4 = '//trim(str(n=self%a_4))//NL
   desc = desc//prefix_//'u23 = '//trim(str(n=self%u23))//NL
   desc = desc//prefix_//'p23 = '//trim(str(n=self%p23))//NL
   desc = desc//prefix_//'r_2 = '//trim(str(n=self%r_2))//NL
   desc = desc//prefix_//'a_2 = '//trim(str(n=self%a_2))//NL
   desc = desc//prefix_//'r_3 = '//trim(str(n=self%r_3))//NL
   desc = desc//prefix_//'a_3 = '//trim(str(n=self%a_3))//NL
   desc = desc//prefix_//'s_1 = '//trim(str(n=self%s_1))//NL
   desc = desc//prefix_//'s_2 = '//trim(str(n=self%s_2))//NL
   desc = desc//prefix_//'s_3 = '//trim(str(n=self%s_3))//NL
   desc = desc//prefix_//'s_4 = '//trim(str(n=self%s_4))
   endfunction description

   ! public methods
   elemental subroutine compute_fluxes(self, eos_left, eos_right, normal, fluxes)
   !< Compute fluxes at interface `x=xo`.
   !<
   !< Sampling the pattern, the interface states are computed.
   class(riemann_solver_compressible_object), intent(in)    :: self      !< Solver.
   class(eos_object),                         intent(in)    :: eos_left  !< Equation of state for left state.
   class(eos_object),                         intent(in)    :: eos_right !< Equation of state for right state.
   type(vector),                              intent(in)    :: normal    !< Normal (versor) of face where fluxes are given.
   class(conservative_object),                intent(inout) :: fluxes    !< Fluxes at interface `x=xo`.
   real(R8P)                                                :: a         !< Speed of sound at interface `x=xo`.
   real(R8P)                                                :: p         !< Pressure at interface `x=xo`.
   real(R8P)                                                :: r         !< Desnity at interface `x=xo`.

   call fluxes%destroy
   associate(s1=>self%s_1, s2=>self%s_2, u23=>self%u23, s3=>self%s_3, s4=>self%s_4, &
             g1=>self%eos_1%g(), d1=>self%eos_1%delta(), e1=>self%eos_1%eta(),      &
             g4=>self%eos_4%g(), d4=>self%eos_4%delta(), e4=>self%eos_4%eta(),      &
             p1=>self%p_1, r1=>self%r_1, u1=>self%u_1, a1=>self%a_1,                &
             p4=>self%p_4, r4=>self%r_4, u4=>self%u_4, a4=>self%a_4,                &
             p23=>self%p23, r2=>self%r_2, r3=>self%r_3)
      select type(fluxes)
      class is(conservative_compressible)
         select case(minloc([-s1, s1 * s2, s2 * u23, u23 * s3, s3 * s4, s4], dim=1))
         case(1) ! left supersonic
            call fluxes%compute_fluxes_from_primitive(eos=eos_left, p=p1, r=r1, u=u1, normal=normal)
         case(2) ! left transonic
            a = (a1 + u1 * d1) / (1._R8P + d1)
            p = p1* (a / a1)**e1
            r = eos_left%density(pressure=p, speed_of_sound=a)
            call fluxes%compute_fluxes_from_primitive(eos=eos_left, p=p, r=r, u=a, normal=normal)
         case(3) ! left subsonic
            call fluxes%compute_fluxes_from_primitive(eos=eos_left, p=p23, r=r2, u=u23, normal=normal)
         case(4) ! right subsonic
            call fluxes%compute_fluxes_from_primitive(eos=eos_right, p=p23, r=r3, u=u23, normal=normal)
         case(5) ! right transonic
            a = (a4 - u4 * d4) / (1._R8P + d4)
            p = p4 * (a / a4) ** e4
            r = eos_right%density(pressure=p, speed_of_sound=a)
            call fluxes%compute_fluxes_from_primitive(eos=eos_right, p=p, r=r, u=-a, normal=normal)
         case(6) ! right supersonic
            call fluxes%compute_fluxes_from_primitive(eos=eos_right, p=p4, r=r4, u=u4, normal=normal)
         endselect
      endselect
   endassociate
   endsubroutine compute_fluxes

   elemental subroutine compute_states23_from_u23(self, p_2, p_3)
   !< Compute interstates 2 and 3 given (an approximation of) veloctiy `S=u23`.
   !<
   !< @Note the pressure of interstates, that should be equal, are returned into separate arguments for allowing a *convergence
   !< checking*, namely if the approximation of `u23` is exact the output is `p_2=p_3=p23`.
   class(riemann_solver_compressible_object), intent(inout) :: self !< Solver.
   real(R8P),                                 intent(out)   :: p_2  !< Pressure of state 2.
   real(R8P),                                 intent(out)   :: p_3  !< Pressure of state 3.

   ! left wave
   if (abs(self%u23 - self%u_1) <= ZeroR8) then
       call compute_post_rarefaction(eos=self%eos_1, sgn=-1._R8P,                        &
                                     u0=self%u_1, p0=self%p_1, a0=self%a_1, ux=self%u23, &
                                     rx=self%r_2, px=     p_2, ax=self%a_2, s0=self%S_1, sx=self%S_2)
   else
      if (self%u23 < self%u_1) then
         call compute_post_shock(eos=self%eos_1, sgn=-1._R8P,                        &
                                 u0=self%u_1, p0=self%p_1, a0=self%a_1, ux=self%u23, &
                                 rx=self%r_2, px=     p_2, ax=self%a_2, ss=self%S_1)
         self%S_2 = self%S_1
      else
         call compute_post_rarefaction(eos=self%eos_1, sgn=-1._R8P,                        &
                                       u0=self%u_1, p0=self%p_1, a0=self%a_1, ux=self%u23, &
                                       rx=self%r_2, px=     p_2, ax=self%a_2, s0=self%S_1, sx=self%S_2)
      endif
   endif
   ! right wave
   if (abs(self%u23 - self%u_4) <= ZeroR8) then
       call compute_post_rarefaction(eos=self%eos_4, sgn=1._R8P,                         &
                                     u0=self%u_4, p0=self%p_4, a0=self%a_4, ux=self%u23, &
                                     rx=self%r_3, px=     p_3, ax=self%a_3, s0=self%S_4, sx=self%S_3)
   else
      if (self%u23 > self%u_4) then
         call compute_post_shock(eos=self%eos_4, sgn=1._R8P,                         &
                                 u0=self%u_4, p0=self%p_4, a0=self%a_4, ux=self%u23, &
                                 rx=self%r_3, px=     p_3, ax=self%a_3, ss=self%S_4)
         self%S_3 = self%S_4
      else
         call compute_post_rarefaction(eos=self%eos_4, sgn=1._R8P,                         &
                                       u0=self%u_4, p0=self%p_4, a0=self%a_4, ux=self%u23, &
                                       rx=self%r_3, px=     p_3, ax=self%a_3, s0=self%S_4, sx=self%S_3)
      endif
   endif
  endsubroutine compute_states23_from_u23

   elemental subroutine set_states14(self, eos_left, state_left, eos_right, state_right, normal)
   !< Set states 1 and 4.
   class(riemann_solver_compressible_object), intent(inout) :: self        !< Solver.
   class(eos_object),                         intent(in)    :: eos_left    !< Equation of state for left state.
   class(conservative_object),                intent(in)    :: state_left  !< Left Riemann state.
   class(eos_object),                         intent(in)    :: eos_right   !< Equation of state for right state.
   class(conservative_object),                intent(in)    :: state_right !< Right Riemann state.
   type(vector),                              intent(in)    :: normal      !< Normal (versor) of face where fluxes are given.

   select type(eos_left)
   class is(eos_compressible)
      self%eos_1 = eos_left
   endselect
   select type(state_left)
   class is(conservative_compressible)
      self%u_1  = state_left%velocity().dot.normal
      self%p_1  = state_left%pressure(eos=eos_left)
      self%r_1  = state_left%density
      self%a_1  = eos_left%speed_of_sound(density=state_left%density, pressure=self%p_1)
   endselect
   select type(eos_right)
   class is(eos_compressible)
      self%eos_4 = eos_right
   endselect
   select type(state_right)
   class is(conservative_compressible)
      self%u_4 = state_right%velocity().dot.normal
      self%p_4 = state_right%pressure(eos=eos_right)
      self%r_4 = state_right%density
      self%a_4 = eos_right%speed_of_sound(density=state_right%density, pressure=self%p_4)
   endselect
   endsubroutine set_states14

   ! non TBP
   elemental subroutine compute_post_rarefaction(eos, sgn, u0, p0, a0, ux, rx, px, ax, s0, sx)
   !< Compute an unknown state `x` from a known state `0` when the two states are separated by a rarefaction, given the velocity
   !< `ux`.
   !<
   !< The `sgn` dummy argument indicates if the rarefaction propagates on `u-a (sgn=-1)` or `u+a (sgn=1)`.
   class(eos_object), intent(in)  :: eos        !< Equation of state.
   real(R8P),         intent(in)  :: sgn        !< Sign for distinguishing *left* (-1) from *right* (1) wave.
   real(R8P),         intent(in)  :: u0, p0, a0 !< Known state (speed, pressure and speed of sound).
   real(R8P),         intent(in)  :: ux         !< Known speed of unknown state.
   real(R8P),         intent(out) :: rx, px, ax !< Unknown pressure and density.
   real(R8P),         intent(out) :: s0, sx     !< Wave speeds (head and back fronts).

   ax = a0 + sgn * eos%delta() * (ux - u0)          ! unknown speed of sound
   px = p0 * ((ax / a0) ** (eos%eta()))             ! unknown pressure
   rx = eos%density(pressure=px, speed_of_sound=ax) ! unknown density
   s0 = u0 + sgn * a0                               ! left wave speed
   sx = ux + sgn * ax                               ! right wave speed
   endsubroutine compute_post_rarefaction

   elemental subroutine compute_post_shock(eos, sgn, u0, p0, a0, ux, rx, px, ax, ss)
   !< Computing an unknown state `x` from a known state `0` when the two states are separated by a shock, given the velocity
   !< `ux`.
   !<
   !< The `sgn` dummy argument indicates if the shock propagates on `u-a (sgn=-1)` or `u+a (sgn=1)`.
   class(eos_object), intent(in)  :: eos        !< Equation of state.
   real(R8P),         intent(in)  :: sgn        !< Sign for distinguishing *left* (-1) from *right* (1) wave.
   real(R8P),         intent(in)  :: u0, p0, a0 !< Known state (speed, pressure and speed of sound).
   real(R8P),         intent(in)  :: ux         !< Unknown speed.
   real(R8P),         intent(out) :: rx, px, ax !< Unknown state (density, pressure and speed of sound).
   real(R8P),         intent(out) :: ss         !< Shock wave speed.
   real(R8P)                      :: M0         !< Relative Mach number of known state.
   real(R8P)                      :: x          !< Dummy variable.

   x   = 0.25_R8P * eos%gp1() * (ux - u0) / a0                     ! dummy variable
   M0  = x + sgn * sqrt(1.0_R8P + x * x)                           ! relative Mach number of known state
   x   = 1._R8P + 2._R8P * eos%g() * (M0 * M0 - 1._R8P) /eos%gp1() ! dummy variable (pressure ratio px/p0)

   ax = a0 * sqrt((eos%gp1() + eos%gm1() * x)/(eos%gp1() + eos%gm1() / x)) ! unknown speed of sound
   px = p0 * x                                                             ! unknown pressure
   rx = eos%density(pressure=px, speed_of_sound=ax)                        ! unknown density
   ss = u0 + a0 * M0                                                       ! shock wave speed
   endsubroutine compute_post_shock
endmodule foreseer_riemann_solver_compressible_object
