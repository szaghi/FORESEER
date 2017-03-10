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
      procedure(cp_interface),             pass(self), deferred :: cp             !< Return specific heat at constant pressure.
      procedure(cv_interface),             pass(self), deferred :: cv             !< Return specific heat at constant volume.
      procedure(density_interface),        pass(self), deferred :: density        !< Return density.
      procedure(energy_interface),         pass(self), deferred :: energy         !< Return specific internal energy.
      procedure(gam_interface),            pass(self), deferred :: gam            !< Return specific heats ratio `gamma=cp/cv`.
      procedure(pressure_interface),       pass(self), deferred :: pressure       !< Return pressure.
      procedure(R_interface),              pass(self), deferred :: R              !< Return fluid constant `R=cp-cv`.
      procedure(speed_of_sound_interface), pass(self), deferred :: speed_of_sound !< Return speed of sound.
      procedure(temperature_interface),    pass(self), deferred :: temperature    !< Return temperature.
endtype eos_object

abstract interface
   !< Abstract interfaces of deferred methods of [[eos_object]].
   elemental function cp_interface(self) result(cp_)
   !< Return specific heat at constant pressure.
   import :: eos_object, R8P
   class(eos_object), intent(in) :: self !< Equation of state.
   real(R8P)                     :: cp_  !< `cp` value.
   endfunction cp_interface

   elemental function cv_interface(self) result(cv_)
   !< Return specific heat at constant volume.
   import :: eos_object, R8P
   class(eos_object), intent(in) :: self !< Equation of state.
   real(R8P)                     :: cv_  !< `cv` value.
   endfunction cv_interface

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

   elemental function energy_interface(self, density, pressure, temperature) result(energy_)
   !< Return specific internal energy.
   import :: eos_object, R8P
   class(eos_object), intent(in)           :: self        !< Equation of state.
   real(R8P),         intent(in), optional :: density     !< Density value.
   real(R8P),         intent(in), optional :: pressure    !< Pressure value.
   real(R8P),         intent(in), optional :: temperature !< Temperature value.
   real(R8P)                               :: energy_     !< Energy value.
   endfunction energy_interface

   elemental function gam_interface(self) result(gam_)
   !< Return specific heats ratio `gamma=cp/cv`.
   import :: eos_object, R8P
   class(eos_object), intent(in) :: self !< Equation of state.
   real(R8P)                     :: gam_ !< Specific heats ratio value.
   endfunction gam_interface

   elemental function pressure_interface(self, density, energy, temperature) result(pressure_)
   !< Return pressure.
   import :: eos_object, R8P
   class(eos_object), intent(in)           :: self        !< Equation of state.
   real(R8P),         intent(in), optional :: density     !< Density value.
   real(R8P),         intent(in), optional :: energy      !< Specific internal energy value.
   real(R8P),         intent(in), optional :: temperature !< Temperature value.
   real(R8P)                               :: pressure_   !< Pressure value.
   endfunction pressure_interface

   elemental function R_interface(self) result(R_)
   !< Return fluid constant `R=cp-cv`.
   import :: eos_object, R8P
   class(eos_object), intent(in) :: self !< Equation of state.
   real(R8P)                     :: R_   !< Fluid constant value.
   endfunction R_interface

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
