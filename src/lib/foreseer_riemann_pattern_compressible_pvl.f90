!< Define the Riemann (waves) pattern for ideal compressible PVL-approximated fluid for FORESEER library.

module foreseer_riemann_pattern_compressible_pvl
!< Define the Riemann (waves) pattern for ideal compressible PVL-approximated fluid for FORESEER library.

use foreseer_conservative_compressible, only : conservative_compressible
use foreseer_conservative_object, only : conservative_object
use foreseer_eos_object, only : eos_object
use foreseer_riemann_pattern_object, only : riemann_pattern_object
use foreseer_riemann_pattern_compressible_object, only : riemann_pattern_compressible_object
use penf, only : R8P
use vecfor, only : vector

implicit none
private
public :: riemann_pattern_compressible_pvl

type, extends(riemann_pattern_compressible_object) :: riemann_pattern_compressible_pvl
   !< Riemann (waves) pattern for ideal compressible PVL-approximated fluid object class.
   procedure(compute_interface), pointer :: compute_ => compute_pvl_u23 !< Compute whole pattern, actual algorithm.
   contains
      ! deferred methods
      procedure, pass(self) :: compute !< Compute whole pattern.
      ! public methods
      procedure, pass(self) :: initialize !< Initialize pattern.
      ! private methods
      procedure, pass(self), private :: compute_p23        !< Compute interstates pressure.
      procedure, pass(self), private :: compute_pvl_u23    !< Compute whole pattern by `u23` algorithm.
      procedure, pass(self), private :: compute_pvl_up23   !< Compute whole pattern by `up23` algorithm.
      procedure, pass(self), private :: compute_u23        !< Compute interstates velocity.
      procedure, pass(self), private :: compute_up23       !< Compute interstates velocity and pressure.
      procedure, pass(self), private :: compute_waves_u23  !< Compute waves speed by `u23` algorithm.
      procedure, pass(self), private :: compute_waves_up23 !< Compute waves speed by `up23` algorithm.
endtype riemann_pattern_compressible_pvl

abstract interface
   pure subroutine compute_interface(self, eos_left, state_left, eos_right, state_right, normal)
   !< Compute whole pattern.
   import :: conservative_object, eos_object, riemann_pattern_compressible_pvl, vector
   class(riemann_pattern_compressible_pvl), intent(inout) :: self        !< Riemann pattern.
   class(eos_object),                       intent(in)    :: eos_left    !< Equation of state for left state.
   class(conservative_object),              intent(in)    :: state_left  !< Left Riemann state.
   class(eos_object),                       intent(in)    :: eos_right   !< Equation of state for right state.
   class(conservative_object),              intent(in)    :: state_right !< Right Riemann state.
   type(vector),                            intent(in)    :: normal      !< Normal (versor) of face where fluxes are given.
   endsubroutine compute_interface
endinterface

contains
   ! deferred methods
   pure subroutine compute(self, eos_left, state_left, eos_right, state_right, normal)
   !< Compute whole pattern.
   class(riemann_pattern_compressible_pvl), intent(inout) :: self        !< Riemann pattern.
   class(eos_object),                       intent(in)    :: eos_left    !< Equation of state for left state.
   class(conservative_object),              intent(in)    :: state_left  !< Left Riemann state.
   class(eos_object),                       intent(in)    :: eos_right   !< Equation of state for right state.
   class(conservative_object),              intent(in)    :: state_right !< Right Riemann state.
   type(vector),                            intent(in)    :: normal      !< Normal (versor) of face where fluxes are given.

   call self%compute_(eos_left=eos_left, state_left=state_left, eos_right=eos_right, state_right=state_right, normal=normal)
   endsubroutine compute

   subroutine initialize(self, config)
   !< Initialize pattern.
   class(riemann_pattern_compressible_pvl), intent(inout)        :: self    !< Riemann pattern.
   character(len=*),                        intent(in), optional :: config  !< Configuration for pattern computation.
   character(len=:), allocatable                                 :: config_ !< Configuration for pattern computation, local var.

   self%compute_ => compute_pvl_u23
   config_ = '' ; if (present(config)) config_ = config
   select case(config_)
   case('u23')
     self%compute_ => compute_pvl_u23
   case('up23')
     self%compute_ => compute_pvl_up23
   case('upr23')
     self%compute_ => compute_pvl_up23
   endselect
   endsubroutine initialize

   ! private methods
   elemental subroutine compute_p23(self)
   !< Compute interstates pressure.
   class(riemann_pattern_compressible_pvl), intent(inout) :: self !< Riemann pattern.

   self%p23 = 0.5_R8P * ((self%p_1 + self%p_4) - 0.25_R8P * (self%u_4 - self%u_1) * (self%r_1 + self%r_4) * (self%a_1 + self%a_4))
   endsubroutine compute_p23

   pure subroutine compute_pvl_u23(self, eos_left, state_left, eos_right, state_right, normal)
   !< Compute whole pattern by `u23` algorithm.
   !<
   !< Use Primitive Variables Linearization algorithm by means of only `u23` approximation.
   class(riemann_pattern_compressible_pvl), intent(inout) :: self        !< Riemann pattern.
   class(eos_object),                       intent(in)    :: eos_left    !< Equation of state for left state.
   class(conservative_object),              intent(in)    :: state_left  !< Left Riemann state.
   class(eos_object),                       intent(in)    :: eos_right   !< Equation of state for right state.
   class(conservative_object),              intent(in)    :: state_right !< Right Riemann state.
   type(vector),                            intent(in)    :: normal      !< Normal (versor) of face where fluxes are given.

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
   call self%compute_u23
   call self%compute_waves_u23(eos_left=eos_left, eos_right=eos_right)
   endsubroutine compute_pvl_u23

   pure subroutine compute_pvl_up23(self, eos_left, state_left, eos_right, state_right, normal)
   !< Compute whole pattern `up23` algorithm.
   !<
   !< Use Primitive Variables Linearization algorithm by means of only `up23` approximation.
   class(riemann_pattern_compressible_pvl), intent(inout) :: self        !< Riemann pattern.
   class(eos_object),                       intent(in)    :: eos_left    !< Equation of state for left state.
   class(conservative_object),              intent(in)    :: state_left  !< Left Riemann state.
   class(eos_object),                       intent(in)    :: eos_right   !< Equation of state for right state.
   class(conservative_object),              intent(in)    :: state_right !< Right Riemann state.
   type(vector),                            intent(in)    :: normal      !< Normal (versor) of face where fluxes are given.

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
   call self%compute_up23
   call self%compute_waves_up23(eos_left=eos_left, eos_right=eos_right)
   endsubroutine compute_pvl_up23

   elemental subroutine compute_u23(self)
   !< Compute interstates velocity.
   class(riemann_pattern_compressible_pvl), intent(inout) :: self !< Riemann pattern.

   self%u23 = 0.5_R8P * (self%u_1 + self%u_4) - 2.0_R8P * (self%p_4 - self%p_1) / ((self%r_1 + self%r_4) * (self%a_1 + self%a_4))
   endsubroutine compute_u23

   elemental subroutine compute_up23(self)
   !< Compute interstates velocity and pressure.
   class(riemann_pattern_compressible_pvl), intent(inout) :: self !< Riemann pattern.
   real(R8P)                                              :: ram  !< Mean value of `r * a`.

   ram = 0.25_R8P * (self%r_1 + self%r_4) * (self%a_1 + self%a_4)
   self%u23 = 0.5_R8P * ((self%u_1 + self%u_4) - (self%p_4 - self%p_1) / ram)
   self%p23 = 0.5_R8P * ((self%p_1 + self%p_4) - (self%u_4 - self%u_1) * ram)
   endsubroutine compute_up23

   elemental subroutine compute_waves_u23(self, eos_left, eos_right)
   !< Compute waves speed `u23` algorithm.
   !<
   !< Use Primitive Variables Linearization algorithm by means of only `u23` approximation.
   class(riemann_pattern_compressible_pvl), intent(inout) :: self      !< Riemann pattern.
   class(eos_object),                       intent(in)    :: eos_left  !< Equation of state for left state.
   class(eos_object),                       intent(in)    :: eos_right !< Equation of state for right state.
   real(R8P)                                              :: x         !< Dummy variable.

   ! compute left state
   if (self%u23 < self%u_1) then ! shock
     x  = 0.25_R8P * (eos_left%gam() + 1._R8P) * (self%u23 - self%u_1) / self%a_1
     self%s_1 = self%u_1 + self%a_1 * (x - sqrt(1.0_R8P + x * x))
   else ! rarefaction
     self%s_1 = self%u_1 - self%a_1
   endif
   ! compute right state
   if (self%u23 > self%u_4) then ! shock
     x  = 0.25_R8P * (eos_right%gam() + 1._R8P) * (self%u23 - self%u_4) / self%a_4
     self%s_4 = self%u_4 + self%a_4 * (x + sqrt(1.0_R8P + x * x))
   else ! rarefaction
     self%s_4 = self%u_4  + self%a_4
   endif
   endsubroutine compute_waves_u23

   elemental subroutine compute_waves_up23(self, eos_left, eos_right)
   !< Compute waves speed `u23` algorithm.
   !<
   !< Use Primitive Variables Linearization algorithm by means of only `up23` approximation.
   class(riemann_pattern_compressible_pvl), intent(inout) :: self      !< Riemann pattern.
   class(eos_object),                       intent(in)    :: eos_left  !< Equation of state for left state.
   class(eos_object),                       intent(in)    :: eos_right !< Equation of state for right state.

   ! compute left state
   if (self%u23 < self%u_1) then ! shock
      self%s_1 = self%u_1 - self%a_1 * sqrt(1._R8P + 0.5_R8P * (eos_left%gam() + 1._R8P) / eos_left%gam() * &
                                            (self%p23 / self%p_1 - 1._R8P))
   else ! rarefaction
      self%s_1 = self%u_1 - self%a_1
   endif
   ! compute right state
   if (self%u23 > self%u_4) then ! shock
      self%s_4 = self%u_4 + self%a_4 * sqrt(1._R8P + 0.5_R8P * (eos_right%gam() + 1._R8P)/ eos_right%gam() * &
                                            (self%p23/self%p_4 - 1._R8P))
   else ! rarefaction
      self%s_4 = self%u_4  + self%a_4
   endif
   endsubroutine compute_waves_up23
endmodule foreseer_riemann_pattern_compressible_pvl
