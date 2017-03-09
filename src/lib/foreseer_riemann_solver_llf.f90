!< Define the Local Lax-Friedrichs (known also as Rusanov) Riemann solver of FORESEER library.

module foreseer_riemann_solver_llf
!< Define the Local Lax-Friedrichs (known also as Rusanov) Riemann solver of FORESEER library.

use foreseer_conservative_compressible, only : conservative_compressible
use foreseer_conservative_object, only : conservative_object
use foreseer_eos_object, only : eos_object
use foreseer_riemann_pattern_compressible, only : riemann_pattern_compressible
use foreseer_riemann_pattern_object, only : riemann_pattern_object
use foreseer_riemann_solver_object, only : riemann_solver_object
use penf, only : I4P, R8P
use vecfor, only : vector

implicit none
private
public :: riemann_solver_llf

type, extends(riemann_solver_object) :: riemann_solver_llf
   !< Local Lax-Friedrichs (known also as Rusanov) Riemann Solver.
   !<
   !< @note This is the implemention for [[conservative_compressible]] Riemann states.
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

   subroutine solve(self, eos_left, state_left, eos_right, state_right, normal, fluxes, pattern)
   !< Solve Riemann Problem.
   !<
   !< Approximate Riemann Solver based on (local) Lax-Friedrichs (known also as Rusanov) algorithm.
   class(riemann_solver_llf),     intent(inout)         :: self         !< Solver.
   class(eos_object),             intent(in)            :: eos_left     !< Equation of state for left state.
   class(conservative_object),    intent(in)            :: state_left   !< Left Riemann state.
   class(eos_object),             intent(in)            :: eos_right    !< Equation of state for right state.
   class(conservative_object),    intent(in)            :: state_right  !< Right Riemann state.
   type(vector),                  intent(in)            :: normal       !< Normal (versor) of face where fluxes are given.
   class(conservative_object),    intent(out)           :: fluxes       !< Fluxes of the Riemann Problem solution.
   class(riemann_pattern_object), intent(out), optional :: pattern      !< Riemann pattern.
   type(conservative_compressible)                      :: fluxes_left  !< Fluxes of left state.
   type(conservative_compressible)                      :: fluxes_right !< Fluxes of right state.
   type(riemann_pattern_compressible)                   :: pattern_     !< Riemann pattern.
   real(R8P)                                            :: lmax         !< Maximum wave speed estimation.

   call pattern_%compute(eos_left=eos_left, state_left=state_left, eos_right=eos_right, state_right=state_right, normal=normal)
   lmax = max(abs(pattern_%s_1), abs(pattern_%s_4))
   call state_left%compute_fluxes(eos=eos_left, normal=normal, fluxes=fluxes_left)
   call state_right%compute_fluxes(eos=eos_right, normal=normal, fluxes=fluxes_right)
   select type(fluxes)
   type is(conservative_compressible)
      fluxes = 0.5_R8P * (fluxes_left + fluxes_right - (lmax * (state_right - state_left)))
   endselect
   if (present(pattern)) then
      select type(pattern)
      type is(riemann_pattern_compressible)
         pattern = pattern_
      endselect
   endif
   endsubroutine solve
endmodule foreseer_riemann_solver_llf
