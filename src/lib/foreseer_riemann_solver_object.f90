!< Define the abstract Riemann solver of FORESEER library.

module foreseer_riemann_solver_object
!< Define the abstract Riemann solver of FORESEER library.

use foreseer_conservative_object, only : conservative_object
use foreseer_eos_object, only : eos_object
use foreseer_riemann_pattern_object, only : riemann_pattern_object
use vecfor, only : vector

implicit none
private
public :: riemann_solver_object

type, abstract :: riemann_solver_object
   !< Abstract Riemann Solver.
   contains
      ! public deferred methods
      procedure(initialize_interface), pass(self), deferred :: initialize !< Initialize solver.
      procedure(solve_interface),      pass(self), deferred :: solve      !< Solve Riemann Problem.
endtype riemann_solver_object

abstract interface
   !< Abstract interfaces of [[riemann_solver_object]] deferred methods.
   subroutine initialize_interface(self)
   !< Initialize solver.
   import :: riemann_solver_object
   class(riemann_solver_object), intent(inout) :: self !< Solver.
   endsubroutine initialize_interface

   subroutine solve_interface(self, eos_left, state_left, eos_right, state_right, normal, fluxes, pattern)
   !< Solve Riemann Problem.
   import :: conservative_object, eos_object, riemann_pattern_object, riemann_solver_object, vector
   class(riemann_solver_object),  intent(inout)         :: self        !< Solver.
   class(eos_object),             intent(in)            :: eos_left    !< Equation of state for left state.
   class(conservative_object),    intent(in)            :: state_left  !< Left Riemann state.
   class(eos_object),             intent(in)            :: eos_right   !< Equation of state for right state.
   class(conservative_object),    intent(in)            :: state_right !< Right Riemann state.
   type(vector),                  intent(in)            :: normal      !< Normal (versor) of face where fluxes are given.
   class(conservative_object),    intent(out)           :: fluxes      !< Fluxes of the Riemann Problem solution.
   class(riemann_pattern_object), intent(out), optional :: pattern     !< Riemann pattern.
   endsubroutine solve_interface
endinterface
endmodule foreseer_riemann_solver_object
