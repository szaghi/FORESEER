!< Define the abstract conservative state of a Riemann Problem for FORESEER library.

module foreseer_conservative_object
!< Define the abstract conservative state of a Riemann Problem for FORESEER library.

use penf, only : R8P
use vecfor, only : vector

implicit none
private
public :: conservative_object

type, abstract :: conservative_object
  !< Convervative object class.
  contains
     ! deferred methods
     procedure(array_interface),      pass(self), deferred :: array              !< Return serialized array of conservative.
     procedure(destroy_interface),    pass(self), deferred :: destroy            !< Destroy conservative.
     procedure(initialize_interface), pass(self), deferred :: initialize         !< Initialize conservative.
     procedure(fluxes_interface),     pass(self), deferred :: fluxes             !< Return conservative fluxes.
     procedure(assignment_interface), pass(lhs),  deferred :: cons_assign_cons   !< Operator `=`.
     procedure(symmetric_operator),   pass(lhs),  deferred :: cons_multiply_cons !< Operator `*`.
     procedure(real_operator_cons),   pass(rhs),  deferred :: real_multiply_cons !< Operator `real * cons`.
     procedure(symmetric_operator),   pass(lhs),  deferred :: add                !< Operator `+`.
     procedure(symmetric_operator),   pass(lhs),  deferred :: sub                !< Operator `-`.
     ! operators
     generic :: assignment(=) => cons_assign_cons                     !< Overload `=`.
     generic :: operator(+) => add                                    !< Overload `+`.
     generic :: operator(-) => sub                                    !< Overload `-`.
     generic :: operator(*) => cons_multiply_cons, real_multiply_cons
endtype conservative_object

abstract interface
   !< Abstract interfaces of deferred methods of [[conservative_object]].
   pure function array_interface(self) result(array_)
   !< Return serialized array of conservative.
   import :: conservative_object, R8P
   class(conservative_object), intent(in) :: self      !< Conservative.
   real(R8P), allocatable                 :: array_(:) !< Serialized array of conservative.
   endfunction array_interface

   elemental subroutine destroy_interface(self)
   !< Destroy conservative.
   import :: conservative_object
   class(conservative_object), intent(inout) :: self !< Conservative.
   endsubroutine destroy_interface

   subroutine initialize_interface(self, initial_state)
   !< Initialize conservative.
   import :: conservative_object
   class(conservative_object),           intent(inout) :: self          !< Conservative.
   class(conservative_object), optional, intent(in)    :: initial_state !< Initial state.
   endsubroutine initialize_interface

   subroutine fluxes_interface(self, normal, conservative_fluxes)
   !< Return conservative fluxes.
   import :: conservative_object, vector
   class(conservative_object), intent(in)  :: self                !< Conservative.
   type(vector),               intent(in)  :: normal              !< Normal (versor) of face where fluxes are given.
   class(conservative_object), intent(out) :: conservative_fluxes !< Conservative fluxes.
   endsubroutine fluxes_interface

   pure subroutine assignment_interface(lhs, rhs)
   !< Operator `=`.
   import :: conservative_object
   class(conservative_object), intent(inout) :: lhs !< Left hand side.
   class(conservative_object), intent(in)    :: rhs !< Right hand side.
   endsubroutine assignment_interface

   function real_operator_cons(lhs, rhs) result(operator_result)
   !< Operator `real * cons`.
   import :: conservative_object, R8P
   real(R8P),                  intent(in)  :: lhs             !< Left hand side.
   class(conservative_object), intent(in)  :: rhs             !< Right hand side.
   class(conservative_object), allocatable :: operator_result !< Operator result.
   endfunction real_operator_cons

   function symmetric_operator(lhs, rhs) result(operator_result)
   !< Symmetric operator `cons.op.cons`.
   import :: conservative_object
   class(conservative_object), intent(in)  :: lhs             !< Left hand side.
   class(conservative_object), intent(in)  :: rhs             !< Right hand side.
   class(conservative_object), allocatable :: operator_result !< Operator result.
   endfunction symmetric_operator
endinterface
endmodule foreseer_conservative_object
