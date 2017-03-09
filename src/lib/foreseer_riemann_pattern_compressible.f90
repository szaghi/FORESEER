!< Define the Riemann (waves) pattern for ideal compressible fluid for FORESEER library.

module foreseer_riemann_pattern_compressible
!< Define the Riemann (waves) pattern for ideal compressible fluid for FORESEER library.

use foreseer_conservative_compressible, only : conservative_compressible
use foreseer_conservative_object, only : conservative_object
use foreseer_eos_object, only : eos_object
use foreseer_riemann_pattern_object, only : riemann_pattern_object
use penf, only : R8P
use vecfor, only : vector

implicit none
private
public :: riemann_pattern_compressible

type, extends(riemann_pattern_object) :: riemann_pattern_compressible
   !< Riemann (waves) pattern for ideal compressible fluid object class.
   real(R8P) :: u_1=0._R8P !< Velocity (normal) of state 1.
   real(R8P) :: u_4=0._R8P !< Velocity (normal) of state 4.
   real(R8P) :: p_1=0._R8P !< Pressure of state 1.
   real(R8P) :: p_4=0._R8P !< Pressure of state 4.
   real(R8P) :: r_1=0._R8P !< Density of state 1.
   real(R8P) :: r_4=0._R8P !< Density of state 4.
   real(R8P) :: a_1=0._R8P !< Speed of sound of state 1.
   real(R8P) :: a_4=0._R8P !< Speed of sound of state 4.
   real(R8P) :: u23=0._R8P !< Velocity (normal) of intermediate states.
   real(R8P) :: p23=0._R8P !< Pressure of intermediate states.
   real(R8P) :: r_2=0._R8P !< Density of state 2.
   real(R8P) :: r_3=0._R8P !< Density of state 3.
   real(R8P) :: s_1=0._R8P !< Left-front of left wave.
   real(R8P) :: s_2=0._R8P !< Right-front of left wave.
   real(R8P) :: s_3=0._R8P !< Left-front of right wave.
   real(R8P) :: s_4=0._R8P !< Right-front of right wave.
   contains
      ! deferred methods
      procedure, pass(self) :: compute !< Compute pattern given left and right states.
      ! private methods
      procedure, pass(self) :: compute_from_u23 !< Compute waves pattern given the speed u23 of intermediates states.
endtype riemann_pattern_compressible

contains
   ! deferred methods
   elemental subroutine compute(self, eos_left, state_left, eos_right, state_right, normal)
   !< Compute pattern given left and right states.
   class(riemann_pattern_compressible), intent(inout) :: self         !< Riemann pattern.
   class(eos_object),                   intent(in)    :: eos_left     !< Equation of state for left state.
   class(conservative_object),          intent(in)    :: state_left   !< Left Riemann state.
   class(eos_object),                   intent(in)    :: eos_right    !< Equation of state for right state.
   class(conservative_object),          intent(in)    :: state_right  !< Right Riemann state.
   type(vector),                        intent(in)    :: normal       !< Normal (versor) of face where fluxes are given.

   select type(state_left)
   class is(conservative_compressible)
      self%u_1 = state_left%velocity().dot.normal
      self%p_1 = state_left%pressure(eos=eos_left)
      self%r_1 = state_left%density
      self%a_1 = eos_left%speed_of_sound(density=state_left%density, pressure=self%p_1)
   endselect
   select type(state_right)
   class is(conservative_compressible)
      self%u_4 = state_right%velocity().dot.normal
      self%p_4 = state_right%pressure(eos=eos_right)
      self%r_4 = state_right%density
      self%a_4 = eos_right%speed_of_sound(density=state_right%density, pressure=self%p_4)
   endselect
   self%u23 = 0.5_R8P * (self%u_1 + self%u_4) - 2.0_R8P * (self%p_4 - self%p_1) / ((self%r_1 + self%r_4) * (self%a_1 + self%a_4))
   call self%compute_from_u23(eos_left=eos_left, eos_right=eos_right)
   endsubroutine compute

   ! private methods
   elemental subroutine compute_from_u23(self, eos_left, eos_right)
   !< Compute waves pattern given the speed u23 of intermediates states.
   class(riemann_pattern_compressible), intent(inout) :: self         !< Riemann pattern.
   class(eos_object),                   intent(in)    :: eos_left     !< Equation of state for left state.
   class(eos_object),                   intent(in)    :: eos_right    !< Equation of state for right state.
   real(R8P)                                          :: x            !< Dummy variable.

   ! compute left state
   if (self%u23 < self%u_1) then ! shock
     x  = 0.25_R8P * (eos_left%gam() + 1._R8P) * (self%u23 - self%u_1) / self%a_1
     self%s_1 = self%u_1 + self%a_1 * (x - sqrt(1.0_R8P + x * x))
   else
     ! rarefaction
     self%s_1 = self%u_1 - self%a_1
   endif
   ! computing right state
   if (self%u23 > self%u_4) then
     ! shock
     x  = 0.25_R8P * (eos_right%gam() + 1._R8P) * (self%u23 - self%u_4) / self%a_4
     self%s_4 = self%u_4 + self%a_4 * (x + sqrt(1.0_R8P + x * x))
   else
     ! rarefaction
     self%s_4 = self%u_4  + self%a_4
   endif
   endsubroutine compute_from_u23

endmodule foreseer_riemann_pattern_compressible
