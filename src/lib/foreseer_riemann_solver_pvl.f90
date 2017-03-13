!< Define the Primitive Variables Linearization based Riemann solver of FORESEER library.

module foreseer_riemann_solver_pvl
!< Define the Primitive Variables Linearization based Riemann solver of FORESEER library.

use foreseer_conservative_compressible, only : conservative_compressible
use foreseer_conservative_object, only : conservative_object
use foreseer_eos_object, only : eos_object
use foreseer_riemann_pattern_compressible_pvl, only : riemann_pattern_compressible_pvl
use foreseer_riemann_pattern_object, only : riemann_pattern_object
use foreseer_riemann_solver_object, only : riemann_solver_object
use penf, only : I4P, R8P
use vecfor, only : vector

implicit none
private
public :: riemann_solver_pvl

type, extends(riemann_solver_object) :: riemann_solver_pvl
   !< Primitive Variables Linearization based Riemann solver.
   !<
   !< @note This is the implemention for [[conservative_compressible]] Riemann states.
   type(riemann_pattern_compressible_pvl) :: pattern !< Riemann pattern.
   contains
      ! public deferred methods
      procedure, pass(self) :: initialize !< Initialize solver.
      procedure, pass(self) :: solve      !< Solve Riemann Problem.
endtype riemann_solver_pvl

contains
   ! public deferred methods
   subroutine initialize(self)
   !< Initialize solver.
   class(riemann_solver_pvl), intent(inout) :: self !< Solver.

   call self%pattern%initialize(config='upr23')
   endsubroutine initialize

   subroutine solve(self, eos_left, state_left, eos_right, state_right, normal, fluxes, pattern)
   !< Solve Riemann Problem.
   !<
   !< Approximate Riemann Solver based on (local) Lax-Friedrichs (known also as Rusanov) algorithm.
   class(riemann_solver_pvl),     intent(inout)         :: self         !< Solver.
   class(eos_object),             intent(in)            :: eos_left     !< Equation of state for left state.
   class(conservative_object),    intent(in)            :: state_left   !< Left Riemann state.
   class(eos_object),             intent(in)            :: eos_right    !< Equation of state for right state.
   class(conservative_object),    intent(in)            :: state_right  !< Right Riemann state.
   type(vector),                  intent(in)            :: normal       !< Normal (versor) of face where fluxes are given.
   class(conservative_object),    intent(inout)         :: fluxes       !< Fluxes of the Riemann Problem solution.
   class(riemann_pattern_object), intent(out), optional :: pattern      !< Riemann pattern.

   call self%pattern%compute(eos_left=eos_left, state_left=state_left, eos_right=eos_right, state_right=state_right, normal=normal)
   call self%pattern%compute_fluxes(eos_left=eos_left, eos_right=eos_right, normal=normal, fluxes=fluxes)
   if (present(pattern)) then
      select type(pattern)
      type is(riemann_pattern_compressible_pvl)
         pattern = self%pattern
      endselect
   endif
   endsubroutine solve
endmodule foreseer_riemann_solver_pvl
