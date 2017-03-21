!< Define the exact (Newton-iterative) Riemann solver of FORESEER library.

module foreseer_riemann_solver_compressible_exact
!< Define the exact (Newton-iterative) Riemann solver of FORESEER library.

use foreseer_conservative_compressible, only : conservative_compressible, conservative_compressible_pointer
use foreseer_conservative_object, only : conservative_object
use foreseer_eos_object, only : eos_object
use foreseer_riemann_solver_compressible_object, only : riemann_solver_compressible_object
use foreseer_riemann_solver_compressible_pvl, only : riemann_solver_compressible_pvl
use foreseer_riemann_solver_object, only : riemann_solver_object
use penf, only : cton, R8P
use vecfor, only : vector

implicit none
private
public :: riemann_solver_compressible_exact

type, extends(riemann_solver_compressible_object) :: riemann_solver_compressible_exact
   !< Exact (Newton-iterative) Riemann Solver.
   !<
   !< @note This is the implemention for [[conservative_compressible]] Riemann states.
   real(R8P) :: tolerance=1.e-10_R8P !< Tolerance on Newton convergence.
   contains
      ! public deferred methods
      procedure, pass(self) :: compute_waves !< Compute waves pattern.
      procedure, pass(self) :: initialize    !< Initialize solver.
      procedure, pass(self) :: solve         !< Solve Riemann Problem.
endtype riemann_solver_compressible_exact

contains
   ! public deferred methods
   pure subroutine compute_waves(self, eos_left, state_left, eos_right, state_right, normal, waves)
   !< Compute waves pattern.
   class(riemann_solver_compressible_exact), intent(inout) :: self        !< Solver.
   class(eos_object),                        intent(in)    :: eos_left    !< Equation of state for left state.
   class(conservative_object),               intent(in)    :: state_left  !< Left Riemann state.
   class(eos_object),                        intent(in)    :: eos_right   !< Equation of state for right state.
   class(conservative_object),               intent(in)    :: state_right !< Right Riemann state.
   type(vector),                             intent(in)    :: normal      !< Normal (versor) of face where fluxes are given.
   real(R8P),                                intent(out)   :: waves(1:)   !< Waves pattern.
   type(conservative_compressible)                         :: fluxes      !< Fluxes of the Riemann Problem solution.

   call self%solve(eos_left=eos_left, state_left=state_left, eos_right=eos_right, state_right=state_right, normal=normal, &
                   fluxes=fluxes)
   waves(1) = self%S_1
   waves(2) = self%S_2
   waves(3) = self%u23
   waves(4) = self%S_3
   waves(5) = self%S_4
   endsubroutine compute_waves

   subroutine initialize(self, config)
   !< Initialize solver.
   class(riemann_solver_compressible_exact), intent(inout)        :: self    !< Solver.
   character(len=*),                         intent(in), optional :: config  !< Configuration for solver algorithm.
   character(len=:), allocatable                                  :: config_ !< Configuration for solver algorithm, local variable.

   config_ = '1.e-10' ; if (present(config)) config_ = config
   self%tolerance = cton(config_, knd=1._R8P)
   endsubroutine initialize

   pure subroutine solve(self, eos_left, state_left, eos_right, state_right, normal, fluxes)
   !< Solve Riemann Problem.
   !<
   !< Approximate Riemann Solver based on (local) Lax-Friedrichs (known also as Rusanov) algorithm.
   class(riemann_solver_compressible_exact), intent(inout) :: self            !< Solver.
   class(eos_object),                        intent(in)    :: eos_left        !< Equation of state for left state.
   class(conservative_object),               intent(in)    :: state_left      !< Left Riemann state.
   class(eos_object),                        intent(in)    :: eos_right       !< Equation of state for right state.
   class(conservative_object),               intent(in)    :: state_right     !< Right Riemann state.
   type(vector),                             intent(in)    :: normal          !< Normal (versor) of face where fluxes are given.
   class(conservative_object),               intent(inout) :: fluxes          !< Fluxes of the Riemann Problem solution.
   real(R8P)                                               :: dum, alfa, beta !< Dummies coefficients.
   real(R8P)                                               :: p_2, p_3        !< Pessure of state 2 and 3.
   real(R8P)                                               :: dp2, dp3        !< Derivate of pessure (dp/du) of state 2 and 3.

   call self%set_states14(eos_left=eos_left, state_left=state_left, eos_right=eos_right, state_right=state_right, normal=normal)

   ! initiale u23 speed
   if (self%p_1 < self%p_4) then
     dum  = 0.5_R8P * self%eos_4%gm1() / self%eos_4%g() ! (gamma - 1) / (gamma * 2)
   else
     dum  = 0.5_R8P * self%eos_1%gm1() / self%eos_1%g() ! (gamma - 1) / (gamma * 2)
   endif
   alfa = (self%p_1 / self%p_4) ** dum
   beta = alfa * self%eos_1%delta() / self%a_1 + self%eos_4%delta()/ self%a_4
   self%u23 = (alfa - 1.0_R8P) / beta +         &
              0.5_R8P * (self%u_1 + self%u_4) + &
              0.5_R8P * (self%u_1 - self%u_4) * (alfa * self%eos_1%delta() / self%a_1 - self%eos_4%delta()/ self%a_4) / beta

   Newton: do
      call self%compute_states23_from_u23(p_2=p_2, p_3=p_3)
      ! evaluate the Newton-Rapson convergence
      if (abs(1.0_R8P - (p_2 / p_3)) >= self%tolerance) then
         dp2 = -1._R8P * self%eos_1%g() * p_2 / self%a_2
         dp3 =  1._R8P * self%eos_4%g() * p_3 / self%a_3
         self%u23 = self%u23 - ((p_2 - p_3) / (dp2-dp3))
      else
        self%p23 = p_2 ! p_2 ~= p_3
        exit Newton
      endif
   enddo Newton

   call self%compute_fluxes(eos_left=eos_left, eos_right=eos_right, normal=normal, fluxes=fluxes)
   endsubroutine solve
endmodule foreseer_riemann_solver_compressible_exact
