!< Define the HLLC Riemann solver of FORESEER library.

module foreseer_riemann_solver_compressible_hllc
!< Define the HLLC Riemann solver of FORESEER library.

use foreseer_conservative_compressible, only : conservative_compressible, conservative_compressible_pointer
use foreseer_conservative_object, only : conservative_object
use foreseer_eos_object, only : eos_object
use foreseer_riemann_solver_compressible_object, only : riemann_solver_compressible_object
use foreseer_riemann_solver_compressible_pvl, only : riemann_solver_compressible_pvl
use foreseer_riemann_solver_object, only : riemann_solver_object
use penf, only : R8P
use vecfor, only : vector

implicit none
private
public :: riemann_solver_compressible_hllc

type, extends(riemann_solver_compressible_object) :: riemann_solver_compressible_hllc
   !< HLLC (Harten, Lax, Van Leer, Toro) Riemann Solver.
   !<
   !< @note This is the implemention for [[conservative_compressible]] Riemann states.
   type(riemann_solver_compressible_pvl) :: solver_pvl !< PVL Riemann solver.
   contains
      ! public deferred methods
      procedure, pass(self) :: compute_waves !< Compute waves pattern.
      procedure, pass(self) :: initialize    !< Initialize solver.
      procedure, pass(self) :: solve         !< Solve Riemann Problem.
endtype riemann_solver_compressible_hllc

contains
   ! public deferred methods
   pure subroutine compute_waves(self, eos_left, state_left, eos_right, state_right, normal, waves)
   !< Compute waves pattern.
   !<
   !< The PVL approximation is based on a 3 waves pattern where the acoustic waves are reduced to a single linear wave instead
   !< of a non linear fan one.
   class(riemann_solver_compressible_hllc), intent(inout) :: self        !< Solver.
   class(eos_object),                       intent(in)    :: eos_left    !< Equation of state for left state.
   class(conservative_object),              intent(in)    :: state_left  !< Left Riemann state.
   class(eos_object),                       intent(in)    :: eos_right   !< Equation of state for right state.
   class(conservative_object),              intent(in)    :: state_right !< Right Riemann state.
   type(vector),                            intent(in)    :: normal      !< Normal (versor) of face where fluxes are given.
   real(R8P),                               intent(out)   :: waves(1:)   !< Waves pattern.

   call self%solver_pvl%compute_waves(eos_left=eos_left, state_left=state_left, eos_right=eos_right, state_right=state_right, &
                                      normal=normal, waves=waves)
   endsubroutine compute_waves

   subroutine initialize(self, config)
   !< Initialize solver.
   class(riemann_solver_compressible_hllc), intent(inout)        :: self    !< Solver.
   character(len=*),                        intent(in), optional :: config  !< Configuration for solver algorithm.
   character(len=:), allocatable                                 :: config_ !< Configuration for solver algorithm, local variable.

   config_ = 'up23' ; if (present(config)) config_ = config
   call self%solver_pvl%initialize(config=config_)
   endsubroutine initialize

   subroutine solve(self, eos_left, state_left, eos_right, state_right, normal, fluxes)
   !< Solve Riemann Problem.
   !<
   !< Approximate Riemann Solver based on (local) Lax-Friedrichs (known also as Rusanov) algorithm.
   class(riemann_solver_compressible_hllc), intent(inout) :: self         !< Solver.
   class(eos_object),                       intent(in)    :: eos_left     !< Equation of state for left state.
   class(conservative_object),              intent(in)    :: state_left   !< Left Riemann state.
   class(eos_object),                       intent(in)    :: eos_right    !< Equation of state for right state.
   class(conservative_object),              intent(in)    :: state_right  !< Right Riemann state.
   type(vector),                            intent(in)    :: normal       !< Normal (versor) of face where fluxes are given.
   class(conservative_object),              intent(inout) :: fluxes       !< Fluxes of the Riemann Problem solution.
   type(conservative_compressible)                        :: state23      !< Intermediate states.
   type(conservative_compressible), pointer               :: state_left_  !< Left Riemann state, local variable.
   type(conservative_compressible), pointer               :: state_right_ !< Right Riemann state, local variable.
   real(R8P)                                              :: waves(1:5)   !< Waves speed pattern.
   real(R8P)                                              :: u23          !< Maximum wave speed estimation.

   state_left_ => conservative_compressible_pointer(to=state_left)
   state_right_ => conservative_compressible_pointer(to=state_right)
   call self%solver_pvl%compute_waves(eos_left=eos_left, state_left=state_left, eos_right=eos_right, state_right=state_right, &
                                      normal=normal, waves=waves)
   associate(r_1=>self%solver_pvl%r_1, u_1=>self%solver_pvl%u_1, p_1=>self%solver_pvl%p_1, g_1=>self%solver_pvl%eos_1%g(), &
             r_4=>self%solver_pvl%r_4, u_4=>self%solver_pvl%u_4, p_4=>self%solver_pvl%p_4, g_4=>self%solver_pvl%eos_4%g(), &
             s_1=>self%solver_pvl%s_1, s_4=>self%solver_pvl%s_4,                                                           &
             E_1=>state_left_%energy/state_left_%density, E_4=>state_right_%energy/state_right_%density)
      u23 = (r_4 * u_4 * (s_4 - u_4) - r_1 * u_1 * (s_1 - u_1) + p_1 - p_4) / &
            (r_4 * (s_4 - u_4) - r_1 * (s_1 - u_1))
      select case(minloc([-s_1, s_1 * u23, u23 * s_4, s_4], dim=1))
      case(1)
         call state_left%compute_fluxes(eos=eos_left, normal=normal, fluxes=fluxes)
      case(2)
         call state_left%compute_fluxes(eos=eos_left, normal=normal, fluxes=fluxes)
         state23%density  = r_1 * (s_1 - u_1) / (s_1 - u23)
         state23%momentum = state23%density * u23 * normal
         state23%energy   = state23%density * (E_1 + (u23 - u_1) * (u23 + p_1 / (r_1 * (s_1 - u_1))))
         select type(fluxes)
         type is(conservative_compressible)
            fluxes = fluxes + s_1 * (state23 - state_left_)
         endselect
      case(3)
         call state_right%compute_fluxes(eos=eos_right, normal=normal, fluxes=fluxes)
         state23%density  = r_4 * (s_4 - u_4) / (s_4 - u23)
         state23%momentum = state23%density * u23 * normal
         state23%energy   = state23%density * (E_4 + (u23 - u_4) * (u23 + p_4 / (r_4 * (s_4 - u_4))))
         select type(fluxes)
         type is(conservative_compressible)
            fluxes = fluxes + s_4 * (state23 - state_right_)
         endselect
      case(4)
         call state_right%compute_fluxes(eos=eos_right, normal=normal, fluxes=fluxes)
      endselect
   endassociate

   endsubroutine solve
endmodule foreseer_riemann_solver_compressible_hllc
