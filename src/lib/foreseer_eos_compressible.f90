!< Define the equation of state (EOS) of ideal compressible fluid for FORESEER library.

module foreseer_eos_compressible
!< Define the equation of state (EOS) of ideal compressible fluid for FORESEER library.

use foreseer_eos_object, only : eos_object
use penf, only : R8P

implicit none
private
public :: eos_compressible

type, extends(eos_object) :: eos_compressible
   !< Equation of state (EOS) of ideal compressible object class.
   private
   real(R8P) :: cp_=0._R8P !< Specific heat at constant pressure `cp`.
   real(R8P) :: cv_=0._R8P !< Specific heat at constant volume `cv`.
   contains
      ! deferred methods
      procedure, pass(self) :: cp          !< Return specific heat at constant pressure `cp` value.
      procedure, pass(self) :: cv          !< Return specific heat at constant volume `cv` value.
      procedure, pass(self) :: density     !< Return density value.
      procedure, pass(self) :: energy      !< Return specific internal energy value.
      procedure, pass(self) :: gam         !< Return specific heats ratio `gamma=cp/cv` value.
      procedure, pass(self) :: pressure    !< Return pressure value.
      procedure, pass(self) :: R           !< Return fluid constant `R=cp-cv` value.
      procedure, pass(self) :: temperature !< Return temperature value.
endtype eos_compressible

interface eos_compressible
   !< Overload [[eos_compressible]] name with its constructor.
   module procedure eos_compressible_instance
endinterface

contains
   ! deferred methods
   elemental function cp(self) result(cp_)
   !< Return specific heat at constant pressure `cp` value.
   class(eos_compressible), intent(in) :: self !< Equation of state.
   real(R8P)                           :: cp_  !< `cp` value.

   cp_ = self%cp_
   endfunction cp

   elemental function cv(self) result(cv_)
   !< Return specific heat at constant volume `cv` value.
   class(eos_compressible), intent(in) :: self !< Equation of state.
   real(R8P)                           :: cv_  !< `cv` value.

   cv_ = self%cv_
   endfunction cv

   elemental function density(self, energy, pressure, temperature) result(density_)
   !< Return density value.
   class(eos_compressible), intent(in)           :: self        !< Equation of state.
   real(R8P),               intent(in), optional :: energy      !< Specific internal energy value.
   real(R8P),               intent(in), optional :: pressure    !< Pressure value.
   real(R8P),               intent(in), optional :: temperature !< Temperature value.
   real(R8P)                                     :: density_    !< Density value.

   density_ = 0._R8P
   if (present(energy).and.present(pressure)) then
      density_ = pressure / ((self%gam() - 1._R8P) * energy)
   elseif (present(pressure).and.present(temperature)) then
      density_ = pressure / (self%R() * temperature)
   endif
   endfunction density

   elemental function energy(self, density, pressure, temperature) result(energy_)
   !< Return specific internal energy value.
   class(eos_compressible), intent(in)           :: self        !< Equation of state.
   real(R8P),               intent(in), optional :: density     !< Density value.
   real(R8P),               intent(in), optional :: pressure    !< Pressure value.
   real(R8P),               intent(in), optional :: temperature !< Temperature value.
   real(R8P)                                     :: energy_     !< Energy value.

   energy_ = 0._R8P
   if (present(density).and.present(pressure)) then
      energy_ = pressure / ((self%gam() - 1._R8P) * density)
   elseif (present(temperature)) then
      energy_ = self%cv() * temperature
   endif
   endfunction energy

   elemental function gam(self) result(gam_)
   !< Return specific heats ratio `gamma=cp/cv` value.
   class(eos_compressible), intent(in) :: self !< Equation of state.
   real(R8P)                           :: gam_ !< Specific heats ratio value.

   gam_ =  self%cp_ / self%cv_
   endfunction gam

   elemental function pressure(self, density, energy, temperature) result(pressure_)
   !< Return pressure value.
   class(eos_compressible), intent(in)           :: self        !< Equation of state.
   real(R8P),               intent(in), optional :: density     !< Density value.
   real(R8P),               intent(in), optional :: energy      !< Specific internal energy value.
   real(R8P),               intent(in), optional :: temperature !< Temperature value.
   real(R8P)                                     :: pressure_   !< Pressure value.

   pressure_ = 0._R8P
   if (present(density).and.present(energy)) then
      pressure_ = density * (self%gam() - 1._R8P) * energy
   elseif (present(density).and.present(temperature)) then
      pressure_ = density * self%R() * temperature
   endif
   endfunction pressure

   elemental function R(self) result(R_)
   !< Return fluid constant `R=cp-cv` value.
   class(eos_compressible), intent(in) :: self !< Equation of state.
   real(R8P)                           :: R_   !< Fluid constant value.

   R_ =  self%cp_ - self%cv_
   endfunction R

   elemental function temperature(self, density, energy, pressure) result(temperature_)
   !< Return temperature value.
   class(eos_compressible), intent(in)           :: self         !< Equation of state.
   real(R8P),               intent(in), optional :: density      !< Density value.
   real(R8P),               intent(in), optional :: energy       !< Specific internal energy value.
   real(R8P),               intent(in), optional :: pressure     !< Pressure value.
   real(R8P)                                     :: temperature_ !< Temperature value.

   temperature_ = 0._R8P
   if (present(density).and.present(pressure)) then
      temperature_ = pressure / (self%R() * density)
   elseif (present(energy)) then
      temperature_ = energy / self%cv()
   endif
   endfunction temperature

   ! private non TBP
   elemental function eos_compressible_instance(cp, cv, gam, R) result(instance)
   !< Return and instance of [[eos_compressible]].
   !<
   !< @note This procedure is used for overloading [[eos_compressible]] name.
   real(R8P), intent(in), optional :: cp       !< Specific heat at constant pressure `cp` value.
   real(R8P), intent(in), optional :: cv       !< Specific heat at constant volume `cv` value.
   real(R8P), intent(in), optional :: gam      !< Specific heats ratio `gamma=cp/cv` value.
   real(R8P), intent(in), optional :: R        !< Fluid constant `R=cp-cv` value.
   type(eos_compressible)          :: instance !< Instance of [[eos_compressible]].

   if (present(cp).and.present(cv)) then
      instance%cp_ = cp
      instance%cv_ = cv
   elseif (present(gam).and.present(R)) then
      instance%cv_ = R/(gam - 1._R8P)
      instance%cp_ = gam * instance%cv_
   elseif (present(gam).and.present(cp)) then
      instance%cp_ = cp
      instance%cv_ = cp / gam
   elseif (present(gam).and.present(cv)) then
      instance%cp_ = gam * cv
      instance%cv_ = cv
   elseif (present(R).and.present(cp)) then
      instance%cp_ = cp
      instance%cv_ = cp - R
   elseif (present(R).and.present(cv)) then
      instance%cp_ = cv + R
      instance%cv_ = cv
   endif
   endfunction eos_compressible_instance
endmodule foreseer_eos_compressible
