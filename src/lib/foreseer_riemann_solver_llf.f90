!< Define the Local Lax-Friedrichs (known also as Rusanov) Riemann solver of FORESEER library.

module foreseer_riemann_solver_llf
!< Define the Local Lax-Friedrichs (known also as Rusanov) Riemann solver of FORESEER library.

use foreseer_riemann_solver_object, only : riemann_solver_object
use penf, only : I4P, R8P

implicit none
private
public :: riemann_solver_llf

type, extends(riemann_solver_object) :: riemann_solver_llf
   !< Abstract Riemann Solver.
   contains
      ! public deferred methods
      procedure, pass(self) :: initialize !< Initialize solver.
      procedure, pass(self) :: solve      !< Solve Riemann Problem.
endtype riemann_solver_llf

contains
   ! public deferred methods
   subroutine initialize(self)
   !< Initialize solver.
   class(riemann_solver_llf), intent(inout) :: self !< Solver.
   endsubroutine initialize

   subroutine solve(self, state_left, state_right, flux)
   !< Solve Riemann Problem.
   !<
   !< Approximate Riemann Solver based on (local) Lax-Friedrichs (known also as Rusanov) algorithm.
   class(riemann_solver_llf), intent(inout) :: self                       !< Solver.
   real(R8P),                 intent(in)    :: state_left(1:)             !< Left Riemann state.
   real(R8P),                 intent(in)    :: state_right(1:)            !< Right Riemann state.
   real(R8P),                 intent(out)   :: flux(1:)                   !< Fluxes of the Riemann Problem solution.
   real(R8P)                                :: waves14(1:2)               !< Waves speeds 1 and 4.
   real(R8P)                                :: flux1(1:size(flux, dim=1)) !< Fluxes of state 1.
   real(R8P)                                :: flux4(1:size(flux, dim=1)) !< Fluxes of state 4.
   real(R8P)                                :: lmax                       !< Maximum wave speed estimation.

   call compute_waves(state_left=state_left, state_right=state_right, waves14=waves14)
   lmax = max(abs(waves14(1)), abs(waves14(2)))
   call compute_flux(state=state_left,  flux=flux1)
   call compute_flux(state=state_right, flux=flux4)
   flux = 0.5_R8P * (flux1 + flux4 - lmax * (state_right - state_left))
   endsubroutine solve

   ! non TBP
   subroutine compute_flux(state, flux)
   !< Compute flux given a state.
   real(R8P), intent(in)  :: state(1:) !< A Riemann state.
   real(R8P), intent(out) :: flux(1:)  !< Fluxes of the Riemann state..
   endsubroutine compute_flux

   subroutine compute_waves(state_left, state_right, waves14)
   !< Compute waves pattern.
   real(R8P), intent(in)  :: state_left(1:)  !< Left Riemann state.
   real(R8P), intent(in)  :: state_right(1:) !< Right Riemann state.
   real(R8P), intent(out) :: waves14(1:)     !< Waves speeds 1 and 4.
   endsubroutine compute_waves
endmodule foreseer_riemann_solver_llf
