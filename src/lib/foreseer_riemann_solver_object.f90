!< Define the abstract Riemann solver of FORESEER library.

module foreseer_riemann_solver_object
!< Define the abstract Riemann solver of FORESEER library.

use foreseer_conservative_object, only : conservative_object
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

   subroutine solve_interface(self, state_left, state_right, normal, fluxes)
   !< Solve Riemann Problem.
   import :: conservative_object, riemann_solver_object, vector
   class(riemann_solver_object), intent(inout) :: self        !< Solver.
   class(conservative_object),   intent(in)    :: state_left  !< Left Riemann state.
   class(conservative_object),   intent(in)    :: state_right !< Right Riemann state.
   type(vector),                 intent(in)    :: normal      !< Normal (versor) of face where fluxes are given.
   class(conservative_object),   intent(out)   :: fluxes      !< Fluxes of the Riemann Problem solution.
   endsubroutine solve_interface
endinterface
endmodule foreseer_riemann_solver_object
