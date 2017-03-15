!< Define the abstract primitive state of a Riemann Problem for FORESEER library.

module foreseer_primitive_object
!< Define the abstract primitive state of a Riemann Problem for FORESEER library.

use foreseer_eos_object, only : eos_object
use penf, only : R8P
use vecfor, only : vector

implicit none
private
public :: primitive_object

type, abstract :: primitive_object
   !< Convervative object class.
   contains
      ! deferred methods
      procedure(array_interface),       pass(self), deferred :: array              !< Return serialized array of primitive.
      procedure(description_interface), pass(self), deferred :: description        !< Return pretty-printed object description.
      procedure(destroy_interface),     pass(self), deferred :: destroy            !< Destroy primitive.
      procedure(energy_interface),      pass(self), deferred :: energy             !< Return energy value.
      procedure(initialize_interface),  pass(self), deferred :: initialize         !< Initialize primitive.
      procedure(momentum_interface),    pass(self), deferred :: momentum           !< Return momentum vector.
      procedure(assignment_interface),  pass(lhs),  deferred :: prim_assign_prim   !< Operator `=`.
      procedure(prim_operator_real),    pass(lhs),  deferred :: prim_divide_real   !< Operator `prim / real`.
      procedure(prim_operator_real),    pass(lhs),  deferred :: prim_multiply_real !< Operator `prim * real`.
      procedure(symmetric_operator),    pass(lhs),  deferred :: prim_multiply_prim !< Operator `*`.
      procedure(real_operator_prim),    pass(rhs),  deferred :: real_multiply_prim !< Operator `real * prim`.
      procedure(symmetric_operator),    pass(lhs),  deferred :: add                !< Operator `+`.
      procedure(unary_operator),        pass(self), deferred :: positive           !< Unary operator `+ prim`.
      procedure(symmetric_operator),    pass(lhs),  deferred :: sub                !< Operator `-`.
      procedure(unary_operator),        pass(self), deferred :: negative           !< Unary operator `- prim`.
      ! operators
      generic :: assignment(=) => prim_assign_prim                                         !< Overload `=`.
      generic :: operator(+) => add, positive                                              !< Overload `+`.
      generic :: operator(-) => sub, negative                                              !< Overload `-`.
      generic :: operator(*) => prim_multiply_prim, prim_multiply_real, real_multiply_prim !< Overload `*`.
      generic :: operator(/) => prim_divide_real                                           !< Overload `/`.
endtype primitive_object

abstract interface
   !< Abstract interfaces of deferred methods of [[primitive_object]].
   pure function array_interface(self) result(array_)
   !< Return serialized array of primitive.
   import :: primitive_object, R8P
   class(primitive_object), intent(in) :: self      !< Primitive.
   real(R8P), allocatable              :: array_(:) !< Serialized array of primitive.
   endfunction array_interface

   pure function description_interface(self, prefix) result(desc)
   !< Return a pretty-formatted object description.
   import :: primitive_object
   class(primitive_object), intent(in)           :: self   !< Primitive.
   character(*),            intent(in), optional :: prefix !< Prefixing string.
   character(len=:), allocatable                 :: desc   !< Description.
   endfunction description_interface

   elemental subroutine destroy_interface(self)
   !< Destroy primitive.
   import :: primitive_object
   class(primitive_object), intent(inout) :: self !< Primitive.
   endsubroutine destroy_interface

   elemental function energy_interface(self, eos) result(energy_)
   !< Return energy value.
   import :: primitive_object, eos_object, R8P
   class(primitive_object), intent(in) :: self    !< Primitive.
   class(eos_object),       intent(in) :: eos     !< Equation of state.
   real(R8P)                           :: energy_ !< Energy value.
   endfunction energy_interface

   subroutine initialize_interface(self, initial_state)
   !< Initialize primitive.
   import :: primitive_object
   class(primitive_object),           intent(inout) :: self          !< Primitive.
   class(primitive_object), optional, intent(in)    :: initial_state !< Initial state.
   endsubroutine initialize_interface

   elemental function momentum_interface(self) result(momentum_)
   !< Return momentum vector.
   import :: primitive_object, vector
   class(primitive_object), intent(in) :: self      !< Primitive.
   type(vector)                        :: momentum_ !< Momentum vector.
   endfunction momentum_interface

   pure subroutine assignment_interface(lhs, rhs)
   !< Operator `=`.
   import :: primitive_object
   class(primitive_object), intent(inout) :: lhs !< Left hand side.
   class(primitive_object), intent(in)    :: rhs !< Right hand side.
   endsubroutine assignment_interface

   function prim_operator_real(lhs, rhs) result(operator_result)
   !< Operator `prim.op.real`.
   import :: primitive_object, R8P
   class(primitive_object), intent(in)  :: lhs             !< Left hand side.
   real(R8P),               intent(in)  :: rhs             !< Right hand side.
   class(primitive_object), allocatable :: operator_result !< Operator result.
   endfunction prim_operator_real

   function real_operator_prim(lhs, rhs) result(operator_result)
   !< Operator `real * prim`.
   import :: primitive_object, R8P
   real(R8P),               intent(in)  :: lhs             !< Left hand side.
   class(primitive_object), intent(in)  :: rhs             !< Right hand side.
   class(primitive_object), allocatable :: operator_result !< Operator result.
   endfunction real_operator_prim

   function symmetric_operator(lhs, rhs) result(operator_result)
   !< Symmetric operator `prim.op.prim`.
   import :: primitive_object
   class(primitive_object), intent(in)  :: lhs             !< Left hand side.
   class(primitive_object), intent(in)  :: rhs             !< Right hand side.
   class(primitive_object), allocatable :: operator_result !< Operator result.
   endfunction symmetric_operator

   function unary_operator(self) result(operator_result)
   !< Unary operator `.op.prim`.
   import :: primitive_object
   class(primitive_object), intent(in)  :: self            !< Primitive.
   class(primitive_object), allocatable :: operator_result !< Operator result.
   endfunction unary_operator
endinterface
endmodule foreseer_primitive_object
