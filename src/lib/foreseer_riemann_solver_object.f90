!< Define the abstract Riemann solver of FORESEER library.

module foreseer_riemann_solver_object
!< Define the abstract Riemann solver of FORESEER library.

use penf, only : I4P, R8P

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

   subroutine solve_interface(self, state_left, state_right, flux)
   !< Solve Riemann Problem.
   import :: R8P, riemann_solver_object
   class(riemann_solver_object), intent(inout) :: self            !< Solver.
   real(R8P),                    intent(in)    :: state_left(1:)  !< Left Riemann state.
   real(R8P),                    intent(in)    :: state_right(1:) !< Right Riemann state.
   real(R8P),                    intent(out)   :: flux(1:)        !< Fluxes of the Riemann Problem solution.
   endsubroutine solve_interface
endinterface
endmodule foreseer_riemann_solver_object
