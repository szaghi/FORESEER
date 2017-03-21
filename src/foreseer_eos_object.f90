!< Define the abstract equation of state (EOS) for FORESEER library.

module foreseer_eos_object
!< Define the abstract equation of state (EOS) for FORESEER library.

use penf, only : R8P

implicit none
private
public :: eos_object

type, abstract :: eos_object
   !< Equation of State (EOS) object class.
   contains
      ! deferred methods
      procedure(scalar_interface),         pass(self), deferred :: cp             !< Return specific heat at constant pressure.
      procedure(scalar_interface),         pass(self), deferred :: cv             !< Return specific heat at constant volume.
      procedure(density_interface),        pass(self), deferred :: density        !< Return density.
      procedure(scalar_interface),         pass(self), deferred :: delta          !< Return `(gamma - 1) / 2`.
      procedure(description_interface),    pass(self), deferred :: description    !< Return pretty-printed object description.
      procedure(energy_interface),         pass(self), deferred :: energy         !< Return specific internal energy.
      procedure(assignment_interface),     pass(lhs),  deferred :: eos_assign_eos !< Operator `=`.
      procedure(scalar_interface),         pass(self), deferred :: eta            !< Return `2 * gamma / (gamma - 1)`.
      procedure(scalar_interface),         pass(self), deferred :: g              !< Return specific heats ratio `gamma=cp/cv`.
      procedure(scalar_interface),         pass(self), deferred :: gm1            !< Return `gamma - 1`.
      procedure(scalar_interface),         pass(self), deferred :: gp1            !< Return `gamma + 1`.
      procedure(pressure_interface),       pass(self), deferred :: pressure       !< Return pressure.
      procedure(scalar_interface),         pass(self), deferred :: R              !< Return fluid constant `R=cp-cv`.
      procedure(speed_of_sound_interface), pass(self), deferred :: speed_of_sound !< Return speed of sound.
      procedure(temperature_interface),    pass(self), deferred :: temperature    !< Return temperature.
      ! operators
      generic :: assignment(=) => eos_assign_eos !< Overload `=`.
endtype eos_object

abstract interface
   !< Abstract interfaces of deferred methods of [[eos_object]].
   pure subroutine assignment_interface(lhs, rhs)
   !< Operator `=`.
   import :: eos_object
   class(eos_object), intent(inout) :: lhs !< Left hand side.
   class(eos_object), intent(in)    :: rhs !< Right hand side.
   endsubroutine assignment_interface

   elemental function density_interface(self, energy, pressure, speed_of_sound, temperature) result(density_)
   !< Return density.
   import :: eos_object, R8P
   class(eos_object), intent(in)           :: self           !< Equation of state.
   real(R8P),         intent(in), optional :: energy         !< Specific internal energy value.
   real(R8P),         intent(in), optional :: pressure       !< Pressure value.
   real(R8P),         intent(in), optional :: speed_of_sound !< Speed of sound value.
   real(R8P),         intent(in), optional :: temperature    !< Temperature value.
   real(R8P)                               :: density_       !< Density value.
   endfunction density_interface

   pure function description_interface(self, prefix) result(desc)
   !< Return a pretty-formatted object description.
   import :: eos_object
   class(eos_object), intent(in)           :: self   !< Equation of state.
   character(*),      intent(in), optional :: prefix !< Prefixing string.
   character(len=:), allocatable           :: desc   !< Description.
   endfunction description_interface

   elemental function energy_interface(self, density, pressure, temperature) result(energy_)
   !< Return specific internal energy.
   import :: eos_object, R8P
   class(eos_object), intent(in)           :: self        !< Equation of state.
   real(R8P),         intent(in), optional :: density     !< Density value.
   real(R8P),         intent(in), optional :: pressure    !< Pressure value.
   real(R8P),         intent(in), optional :: temperature !< Temperature value.
   real(R8P)                               :: energy_     !< Energy value.
   endfunction energy_interface

   elemental function pressure_interface(self, density, energy, temperature) result(pressure_)
   !< Return pressure.
   import :: eos_object, R8P
   class(eos_object), intent(in)           :: self        !< Equation of state.
   real(R8P),         intent(in), optional :: density     !< Density value.
   real(R8P),         intent(in), optional :: energy      !< Specific internal energy value.
   real(R8P),         intent(in), optional :: temperature !< Temperature value.
   real(R8P)                               :: pressure_   !< Pressure value.
   endfunction pressure_interface

   elemental function scalar_interface(self) result(scalar_)
   !< Return a scalar real value by only `self` data.
   import :: eos_object, R8P
   class(eos_object), intent(in) :: self    !< Equation of state.
   real(R8P)                     :: scalar_ !< Scalar value.
   endfunction scalar_interface

   elemental function speed_of_sound_interface(self, density, pressure) result(speed_of_sound_)
   !< Return speed of sound.
   import :: eos_object, R8P
   class(eos_object), intent(in) :: self            !< Equation of state.
   real(R8P),         intent(in) :: density         !< Density value.
   real(R8P),         intent(in) :: pressure        !< Pressure value.
   real(R8P)                     :: speed_of_sound_ !< Speed of sound value.
   endfunction speed_of_sound_interface

   elemental function temperature_interface(self, density, energy, pressure) result(temperature_)
   !< Return temperature.
   import :: eos_object, R8P
   class(eos_object), intent(in)           :: self         !< Equation of state.
   real(R8P),         intent(in), optional :: density      !< Density value.
   real(R8P),         intent(in), optional :: energy       !< Specific internal energy value.
   real(R8P),         intent(in), optional :: pressure     !< Pressure value.
   real(R8P)                               :: temperature_ !< Temperature value.
   endfunction temperature_interface
endinterface
endmodule foreseer_eos_object
