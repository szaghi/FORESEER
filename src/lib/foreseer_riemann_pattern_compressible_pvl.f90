!< Define the compressible Riemann (states) PVL pattern for FORESEER library.

module foreseer_riemann_pattern_compressible_pvl
!< Define the compressible Riemann (states) PVL pattern for FORESEER library.

use flow_conservative_object, only : conservative_object
use flow_eos_object, only : eos_object
use foreseer_riemann_pattern_compressible_object, only : riemann_pattern_compressible_object
use penf, only : R8P
use vecfor, only : vector

implicit none
private
public :: riemann_pattern_compressible_pvl

type, extends(riemann_pattern_compressible_object) :: riemann_pattern_compressible_pvl
   !< Compressible Riemann (states) PVL pattern object class.
   procedure(compute_waves_interface), pointer :: compute_waves_ => compute_waves_up23 !< Compute waves speed.
   contains
      ! deferred methods
      procedure, pass(self) :: compute_waves !< Compute waves speed.
      ! private methods
      procedure, pass(self), private :: compute_u23        !< Compute interstates velocity.
      procedure, pass(self), private :: compute_up23       !< Compute interstates velocity and pressure.
      procedure, pass(self), private :: compute_waves_u23  !< Compute waves speed by `u23` algorithm.
      procedure, pass(self), private :: compute_waves_up23 !< Compute waves speed by `up23` algorithm.
endtype riemann_pattern_compressible_pvl

abstract interface
   pure subroutine compute_waves_interface(self)
   !< Compute waves pattern.
   import :: riemann_pattern_compressible_pvl
   class(riemann_pattern_compressible_pvl), intent(inout) :: self !< Riemann (states) pattern solution.
   endsubroutine compute_waves_interface
endinterface

contains
   ! deferred methods
   pure subroutine compute_waves(self)
   !< Compute waves speed.
   class(riemann_pattern_compressible_pvl), intent(inout) :: self !< Riemann (states) pattern solution.

   call self%compute_waves_
   endsubroutine compute_waves

   ! public methods
   elemental subroutine compute_u23(self)
   !< Compute interstates velocity.
   class(riemann_pattern_compressible_pvl), intent(inout) :: self !< Riemann (states) pattern solution.

   self%u23 = 0.5_R8P * (self%u_1 + self%u_4) - 2.0_R8P * (self%p_4 - self%p_1) / ((self%r_1 + self%r_4) * (self%a_1 + self%a_4))
   endsubroutine compute_u23

   elemental subroutine compute_up23(self)
   !< Compute interstates velocity and pressure.
   class(riemann_pattern_compressible_pvl), intent(inout) :: self !< Riemann (states) pattern solution.
   real(R8P)                                              :: ram  !< Mean value of `r * a`.

   ram = 0.25_R8P * (self%r_1 + self%r_4) * (self%a_1 + self%a_4)
   self%u23 = 0.5_R8P * ((self%u_1 + self%u_4) - (self%p_4 - self%p_1) / ram)
   self%p23 = 0.5_R8P * ((self%p_1 + self%p_4) - (self%u_4 - self%u_1) * ram)
   endsubroutine compute_up23

   pure subroutine compute_waves_u23(self)
   !< Compute waves speed `u23` algorithm.
   !<
   !< Use Primitive Variables Linearization algorithm by means of only `u23` approximation.
   class(riemann_pattern_compressible_pvl), intent(inout) :: self !< Riemann (states) pattern solution.
   real(R8P)                                              :: x    !< Dummy variable.

   call self%compute_u23

   ! compute left state
   if (self%u23 < self%u_1) then ! shock
     x  = 0.25_R8P * (self%eos_1%g() + 1._R8P) * (self%u23 - self%u_1) / self%a_1
     self%s_1 = self%u_1 + self%a_1 * (x - sqrt(1.0_R8P + x * x))
   else ! rarefaction
     self%s_1 = self%u_1 - self%a_1
   endif
   ! compute right state
   if (self%u23 > self%u_4) then ! shock
     x  = 0.25_R8P * (self%eos_4%g() + 1._R8P) * (self%u23 - self%u_4) / self%a_4
     self%s_4 = self%u_4 + self%a_4 * (x + sqrt(1.0_R8P + x * x))
   else ! rarefaction
     self%s_4 = self%u_4 + self%a_4
   endif
   endsubroutine compute_waves_u23

   pure subroutine compute_waves_up23(self)
   !< Compute waves speed `u23` algorithm.
   !<
   !< Use Primitive Variables Linearization algorithm by means of only `up23` approximation.
   class(riemann_pattern_compressible_pvl), intent(inout) :: self !< Riemann (states) pattern solution.

   call self%compute_up23

   ! compute left state
   if (self%u23 < self%u_1) then ! shock
      self%s_1 = self%u_1 - self%a_1 * sqrt(1._R8P + 0.5_R8P * (self%eos_1%g() + 1._R8P) / &
                                            self%eos_1%g() * (self%p23 / self%p_1 - 1._R8P))
   else ! rarefaction
      self%s_1 = self%u_1 - self%a_1
   endif
   ! compute right state
   if (self%u23 > self%u_4) then ! shock
      self%s_4 = self%u_4 + self%a_4 * sqrt(1._R8P + 0.5_R8P * (self%eos_4%g() + 1._R8P) / &
                                            self%eos_4%g() * (self%p23 / self%p_4 - 1._R8P))
   else ! rarefaction
      self%s_4 = self%u_4  + self%a_4
   endif
   endsubroutine compute_waves_up23
endmodule foreseer_riemann_pattern_compressible_pvl
