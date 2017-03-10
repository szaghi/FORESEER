!< Define the abstract Riemann (waves) pattern object for FORESEER library.

module foreseer_riemann_pattern_object
!< Define the abstract Riemann (waves) pattern object for FORESEER library.

use foreseer_conservative_object, only : conservative_object
use foreseer_eos_object, only : eos_object
use penf, only : R8P
use vecfor, only : vector

implicit none
private
public :: riemann_pattern_object

type, abstract :: riemann_pattern_object
   !< Riemann (waves) pattern object class.
   !<
   !< The pattern arises after the *breaking* of the initial discontinuity of the Riemann Problem.
   contains
      ! deferred methods
      procedure(compute_interface),       pass(self), deferred :: compute        !< Compute whole pattern.
      procedure(description_interface),   pass(self), deferred :: description    !< Return pretty-printed object description.
      procedure(assignment_interface),    pass(lhs),  deferred :: pat_assign_pat !< Operator `=`.
      procedure(waves_extrema_interface), pass(self), deferred :: waves_extrema  !< Return waves speed extrema.
      ! operators
      generic :: assignment(=) => pat_assign_pat !< Overload `=`.
endtype riemann_pattern_object

abstract interface
   !< Abstract interfaces of deferred methods of [[riemann_pattern_object]].
   pure subroutine compute_interface(self, eos_left, state_left, eos_right, state_right, normal)
   !< Compute whole pattern.
   !<
   !< Depending on the concrete extension, this procedure computes all the interstates and waves speed.
   import :: conservative_object, eos_object, riemann_pattern_object, vector
   class(riemann_pattern_object), intent(inout) :: self        !< Riemann pattern.
   class(eos_object),             intent(in)    :: eos_left    !< Equation of state for left state.
   class(conservative_object),    intent(in)    :: state_left  !< Left Riemann state.
   class(eos_object),             intent(in)    :: eos_right   !< Equation of state for right state.
   class(conservative_object),    intent(in)    :: state_right !< Right Riemann state.
   type(vector),                  intent(in)    :: normal      !< Normal (versor) of face where fluxes are given.
   endsubroutine compute_interface

   pure function description_interface(self, prefix) result(desc)
   !< Return a pretty-formatted object description.
   import :: riemann_pattern_object
   class(riemann_pattern_object), intent(in)           :: self   !< Riemann pattern.
   character(*),                  intent(in), optional :: prefix !< Prefixing string.
   character(len=:), allocatable                       :: desc   !< Description.
   endfunction description_interface

   elemental subroutine assignment_interface(lhs, rhs)
   !< Operator `=`.
   import :: riemann_pattern_object
   class(riemann_pattern_object), intent(inout) :: lhs !< Left hand side.
   class(riemann_pattern_object), intent(in)    :: rhs !< Right hand side.
   endsubroutine assignment_interface

   pure function waves_extrema_interface(self) result(waves)
   !< Return waves speed extrema.
   import :: riemann_pattern_object, R8P
   class(riemann_pattern_object), intent(in) :: self       !< Riemann pattern.
   real(R8P)                                 :: waves(1:2) !< Waves speed extrema.
   endfunction waves_extrema_interface
endinterface
endmodule foreseer_riemann_pattern_object
