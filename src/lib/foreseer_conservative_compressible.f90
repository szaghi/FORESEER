!< Define the abstract conservative compressible state of a Riemann Problem for FORESEER library.

module foreseer_conservative_compressible
!< Define the abstract conservative compressible state of a Riemann Problem for FORESEER library.

use foreseer_conservative_object, only : conservative_object
use penf, only : R8P
use vecfor, only : vector

implicit none
private
public :: conservative_compressible

type, extends(conservative_object) :: conservative_compressible
   !< Convervative compressible object class.
   real(R8P)    :: cp=0._R8P      !< Specific heat at constant pressure.
   real(R8P)    :: cv=0._R8P      !< Specific heat at constant volume.
   real(R8P)    :: density=0._R8P !< Density, `rho`.
   type(vector) :: momentum       !< Momentum, `rho * v`, `rho` being the density and `v` the velocity vector.
   real(R8P)    :: energy=0._R8P  !< Energy, `rho * E`, `rho` being the density and `E` the specific energy.
   contains
      ! public methods
      procedure, pass(self) :: pressure             !< Return pressure value.
      procedure, pass(self) :: specific_heats_ratio !< Return specific heats ratio `cp/cv`.
      procedure, pass(self) :: velocity             !< Return velocity vector.
      ! deferred methods
      procedure, pass(self) :: array              !< Return serialized array of conservative.
      procedure, pass(self) :: compute_fluxes     !< Compute conservative fluxes.
      procedure, pass(self) :: destroy            !< Destroy conservative.
      procedure, pass(self) :: initialize         !< Initialize conservative.
      procedure, pass(lhs)  :: cons_assign_cons   !< Operator `=`.
      procedure, pass(lhs)  :: cons_multiply_cons !< Operator `*`.
      procedure, pass(rhs)  :: real_multiply_cons !< Operator `real * cons`.
      procedure, pass(lhs)  :: add                !< Operator `+`.
      procedure, pass(lhs)  :: sub                !< Operator `-`.
endtype conservative_compressible

interface conservative_compressible
   !< Overload [[conservative_compressible]] name with its constructor.
   module procedure conservative_compressible_instance
endinterface

contains
   ! public methods
   elemental function pressure(self) result(pressure_)
   !< Return pressure value.
   class(conservative_compressible), intent(in) :: self      !< Conservative.
   real(R8P)                                    :: pressure_ !< Pressure value.
   type(vector)                                 :: velocity_ !< Velocity vector.

   velocity_ = self%velocity()
   pressure_ = (self%specific_heats_ratio() - 1._R8P) * (self%energy - 0.5_R8P * self%density * velocity_%sq_norm())
   endfunction pressure

   elemental function specific_heats_ratio(self) result(gamma_)
   !< Return specific heats ratio `cv/cv`.
   class(conservative_compressible), intent(in) :: self   !< Conservative.
   real(R8P)                                    :: gamma_ !< Specific heats ratio `cp/cv`.

   gamma_ = self%cp / self%cv
   endfunction specific_heats_ratio

   elemental function velocity(self) result(velocity_)
   !< Return velocity vector.
   class(conservative_compressible), intent(in) :: self      !< Conservative.
   type(vector)                                 :: velocity_ !< Velocity vector.

   velocity_ = self%momentum / self%density
   endfunction velocity

   ! deferred methods
   pure function array(self) result(array_)
   !< Return serialized array of conservative.
   class(conservative_compressible), intent(in) :: self      !< Conservative.
   real(R8P), allocatable                       :: array_(:) !< Serialized array of conservative.

   allocate(array_(1:7))
   array_(1) = self%cp
   array_(2) = self%cv
   array_(3) = self%density
   array_(4) = self%momentum%x
   array_(5) = self%momentum%y
   array_(6) = self%momentum%z
   array_(7) = self%energy
   endfunction array

   subroutine compute_fluxes(self, normal, fluxes)
   !< Compute conservative fluxes.
   class(conservative_compressible), intent(in)  :: self             !< Conservative.
   type(vector),                     intent(in)  :: normal           !< Normal (versor) of face where fluxes are given.
   class(conservative_object),       intent(out) :: fluxes           !< Conservative fluxes.
   real(R8P)                                     :: pressure_        !< Pressure value.
   type(vector)                                  :: velocity_        !< Velocity vector.
   real(R8P)                                     :: velocity_normal_ !< Velocity component parallel to given normal.

   select type(fluxes)
   class is(conservative_compressible)
      pressure_ = self%pressure()
      velocity_ = self%velocity()
      velocity_normal_ = velocity_.dot.normal
      fluxes%cp = self%cp
      fluxes%cv = self%cv
      fluxes%density = self%momentum.dot.normal
      fluxes%momentum = self%density * velocity_ * velocity_normal_ + pressure_ * normal
      fluxes%energy = (self%energy + pressure_) * velocity_normal_
   endselect
   endsubroutine compute_fluxes

   elemental subroutine destroy(self)
   !< Destroy conservative.
   class(conservative_compressible), intent(inout) :: self  !< Conservative.
   type(conservative_compressible)                 :: fresh !< Fresh instance of conservative object.

   self = fresh
   endsubroutine destroy

   subroutine initialize(self, initial_state)
   !< Initialize conservative.
   class(conservative_compressible), intent(inout)        :: self          !< Conservative.
   class(conservative_object),       intent(in), optional :: initial_state !< Initial state.

   if (present(initial_state)) then
      select type(initial_state)
      class is(conservative_compressible)
         self = initial_state
      endselect
   else
      call self%destroy
   endif
   endsubroutine initialize

   ! operators
   pure subroutine cons_assign_cons(lhs, rhs)
   !< Operator `=`.
   class(conservative_compressible), intent(inout) :: lhs !< Left hand side.
   class(conservative_object),       intent(in)    :: rhs !< Right hand side.

   select type(rhs)
   class is (conservative_compressible)
      lhs%cp       = rhs%cp
      lhs%cv       = rhs%cv
      lhs%density  = rhs%density
      lhs%momentum = rhs%momentum
      lhs%energy   = rhs%energy
   endselect
   endsubroutine cons_assign_cons

   function real_multiply_cons(lhs, rhs) result(operator_result)
   !< Operator `real * cons`.
   !<
   !< @note Specific heats `cp, cv` are not affected by this operator: they are copied by right hand side.
   real(R8P),                        intent(in) :: lhs             !< Left hand side.
   class(conservative_compressible), intent(in) :: rhs             !< Right hand side.
   class(conservative_object), allocatable      :: operator_result !< Operator result.

   allocate(conservative_compressible :: operator_result)
   select type(operator_result)
   class is(conservative_compressible)
      operator_result%cp       = rhs%cp
      operator_result%cv       = rhs%cv
      operator_result%density  = lhs * rhs%density
      operator_result%momentum = lhs * rhs%momentum
      operator_result%energy   = lhs * rhs%energy
   endselect
   endfunction real_multiply_cons

   function cons_multiply_cons(lhs, rhs) result(operator_result)
   !< Operator `*`.
   !<
   !< @note Specific heats `cp, cv` are not affected by this operator: they are copied by left hand side.
   class(conservative_compressible), intent(in) :: lhs             !< Left hand side.
   class(conservative_object),       intent(in) :: rhs             !< Right hand side.
   class(conservative_object), allocatable      :: operator_result !< Operator result.

   allocate(conservative_compressible :: operator_result)
   select type(operator_result)
   class is(conservative_compressible)
      operator_result = lhs
      select type(rhs)
      class is (conservative_compressible)
         operator_result%density  = lhs%density  * rhs%density
         operator_result%momentum = lhs%momentum * rhs%momentum
         operator_result%energy   = lhs%energy   * rhs%energy
      endselect
   endselect
   endfunction cons_multiply_cons

   function add(lhs, rhs) result(operator_result)
   !< Operator `+`.
   !<
   !< @note Specific heats `cp, cv` are not affected by this operator: they are copied by left hand side.
   class(conservative_compressible), intent(in) :: lhs             !< Left hand side.
   class(conservative_object),       intent(in) :: rhs             !< Right hand side.
   class(conservative_object), allocatable      :: operator_result !< Operator result.

   allocate(conservative_compressible :: operator_result)
   select type(operator_result)
   class is(conservative_compressible)
      operator_result = lhs
      select type(rhs)
      class is (conservative_compressible)
         operator_result%density  = lhs%density  + rhs%density
         operator_result%momentum = lhs%momentum + rhs%momentum
         operator_result%energy   = lhs%energy   + rhs%energy
      endselect
   endselect
   endfunction add

   function sub(lhs, rhs) result(operator_result)
   !< Operator `+`.
   !<
   !< @note Specific heats `cp, cv` are not affected by this operator: they are copied by left hand side.
   class(conservative_compressible), intent(in) :: lhs             !< Left hand side.
   class(conservative_object),       intent(in) :: rhs             !< Right hand side.
   class(conservative_object), allocatable      :: operator_result !< Operator result.

   allocate(conservative_compressible :: operator_result)
   select type(operator_result)
   class is(conservative_compressible)
      operator_result = lhs
      select type(rhs)
      class is (conservative_compressible)
         operator_result%density  = lhs%density  - rhs%density
         operator_result%momentum = lhs%momentum - rhs%momentum
         operator_result%energy   = lhs%energy   - rhs%energy
      endselect
   endselect
   endfunction sub

   ! private non TBP
   elemental function conservative_compressible_instance(cp, cv, density, momentum, energy) result(instance)
   !< Return and instance of [[conservative_compressible]].
   !<
   !< @note This procedure is used for overloading [[conservative_compressible]] name.
   real(R8P),    intent(in)           :: cp       !< Specific heat at constant pressure.
   real(R8P),    intent(in)           :: cv       !< Specific heat at constant volume.
   real(R8P),    intent(in), optional :: density  !< Density, `rho`.
   type(vector), intent(in), optional :: momentum !< Momentum, `rho * v`, `rho` being the density and `v` the velocity vector.
   real(R8P),    intent(in), optional :: energy   !< Energy, `rho * E`, `rho` being the density and `E` the specific energy.
   type(conservative_compressible)    :: instance !< Instance of [[conservative_compressible]].

   instance%cp = cp
   instance%cv = cv
   if (present(density)) instance%density = density
   if (present(momentum)) instance%momentum = momentum
   if (present(energy)) instance%energy = energy
   endfunction conservative_compressible_instance
endmodule foreseer_conservative_compressible
