!< Define the Primitive Variables Linearization based Riemann solver of FORESEER library.

module foreseer_riemann_solver_compressible_pvl
!< Define the Primitive Variables Linearization based Riemann solver of FORESEER library.

use foreseer_conservative_compressible, only : conservative_compressible
use foreseer_conservative_object, only : conservative_object
use foreseer_eos_object, only : eos_object
use foreseer_riemann_solver_compressible_object, only : riemann_solver_compressible_object
use foreseer_riemann_solver_object, only : riemann_solver_object
use penf, only : I4P, R8P
use vecfor, only : vector

implicit none
private
public :: riemann_solver_compressible_pvl

type, extends(riemann_solver_compressible_object) :: riemann_solver_compressible_pvl
   !< Primitive Variables Linearization based Riemann solver.
   !<
   !< @note This is the implemention for [[conservative_compressible]] Riemann states.
   procedure(compute_waves_interface), pointer :: compute_waves_ => compute_waves_up23 !< Compute waves pattern
   procedure(solve_interface),         pointer :: solve_ => solve_up23                 !< Solve Riemann problem.
   contains
      ! public deferred methods
      procedure, pass(self) :: compute_waves !< Compute waves pattern.
      procedure, pass(self) :: initialize    !< Initialize solver.
      procedure, pass(self) :: solve         !< Solve Riemann Problem.
      ! private methods
      procedure, pass(self), private :: compute_u23        !< Compute interstates velocity.
      procedure, pass(self), private :: compute_up23       !< Compute interstates velocity and pressure.
      procedure, pass(self), private :: compute_waves_u23  !< Compute waves speed by `u23` algorithm.
      procedure, pass(self), private :: compute_waves_up23 !< Compute waves speed by `up23` algorithm.
      procedure, pass(self), private :: solve_u23          !< Compute whole pattern by `u23` algorithm.
      procedure, pass(self), private :: solve_up23         !< Compute whole pattern by `up23` algorithm.
endtype riemann_solver_compressible_pvl

abstract interface
   pure subroutine compute_waves_interface(self, eos_left, state_left, eos_right, state_right, normal, waves)
   !< Compute waves pattern.
   import :: conservative_object, eos_object, riemann_solver_compressible_pvl, R8P, vector
   class(riemann_solver_compressible_pvl), intent(inout) :: self        !< Solver.
   class(eos_object),                      intent(in)    :: eos_left    !< Equation of state for left state.
   class(conservative_object),             intent(in)    :: state_left  !< Left Riemann state.
   class(eos_object),                      intent(in)    :: eos_right   !< Equation of state for right state.
   class(conservative_object),             intent(in)    :: state_right !< Right Riemann state.
   type(vector),                           intent(in)    :: normal      !< Normal (versor) of face where fluxes are given.
   real(R8P),                              intent(out)   :: waves(1:)   !< Waves pattern.
   endsubroutine compute_waves_interface

   pure subroutine solve_interface(self, eos_left, state_left, eos_right, state_right, normal, fluxes)
   !< Solve Riemann problem.
   import :: conservative_object, eos_object, riemann_solver_compressible_pvl, vector
   class(riemann_solver_compressible_pvl), intent(inout) :: self        !< Solver.
   class(eos_object),                      intent(in)    :: eos_left    !< Equation of state for left state.
   class(conservative_object),             intent(in)    :: state_left  !< Left Riemann state.
   class(eos_object),                      intent(in)    :: eos_right   !< Equation of state for right state.
   class(conservative_object),             intent(in)    :: state_right !< Right Riemann state.
   type(vector),                           intent(in)    :: normal      !< Normal (versor) of face where fluxes are given.
   class(conservative_object),             intent(inout) :: fluxes      !< Fluxes of the Riemann Problem solution.
   endsubroutine solve_interface
endinterface

contains
   ! public deferred methods
   pure subroutine compute_waves(self, eos_left, state_left, eos_right, state_right, normal, waves)
   !< Compute waves pattern.
   !<
   !< The PVL approximation is based on a 3 waves pattern where the acoustic waves are reduced to a single linear wave instead
   !< of a non linear fan one.
   class(riemann_solver_compressible_pvl), intent(inout) :: self        !< Solver.
   class(eos_object),                      intent(in)    :: eos_left    !< Equation of state for left state.
   class(conservative_object),             intent(in)    :: state_left  !< Left Riemann state.
   class(eos_object),                      intent(in)    :: eos_right   !< Equation of state for right state.
   class(conservative_object),             intent(in)    :: state_right !< Right Riemann state.
   type(vector),                           intent(in)    :: normal      !< Normal (versor) of face where fluxes are given.
   real(R8P),                              intent(out)   :: waves(1:)   !< Waves pattern.

   call self%compute_waves_(eos_left=eos_left, state_left=state_left, eos_right=eos_right, state_right=state_right, &
                            normal=normal, waves=waves)
   endsubroutine compute_waves

   pure subroutine solve(self, eos_left, state_left, eos_right, state_right, normal, fluxes)
   !< Solve Riemann Problem.
   class(riemann_solver_compressible_pvl), intent(inout) :: self        !< Solver.
   class(eos_object),                      intent(in)    :: eos_left    !< Equation of state for left state.
   class(conservative_object),             intent(in)    :: state_left  !< Left Riemann state.
   class(eos_object),                      intent(in)    :: eos_right   !< Equation of state for right state.
   class(conservative_object),             intent(in)    :: state_right !< Right Riemann state.
   type(vector),                           intent(in)    :: normal      !< Normal (versor) of face where fluxes are given.
   class(conservative_object),             intent(inout) :: fluxes      !< Fluxes of the Riemann Problem solution.

   call self%solve_(eos_left=eos_left, state_left=state_left, eos_right=eos_right, state_right=state_right, &
                    normal=normal, fluxes=fluxes)
   endsubroutine solve

   subroutine initialize(self, config)
   !< Initialize solver.
   class(riemann_solver_compressible_pvl), intent(inout)        :: self    !< Solver.
   character(len=*),                       intent(in), optional :: config  !< Configuration for solver algorithm.
   character(len=:), allocatable                                :: config_ !< Configuration for solver algorithm, local var.

   self%compute_waves_ => compute_waves_u23
   self%solve_ => solve_u23
   config_ = '' ; if (present(config)) config_ = config
   select case(config_)
   case('u23')
     self%compute_waves_ => compute_waves_u23
     self%solve_ => solve_u23
   case('up23')
     self%compute_waves_ => compute_waves_up23
     self%solve_ => solve_up23
   case('upr23')
     self%compute_waves_ => compute_waves_up23
     self%solve_ => solve_up23
   endselect
   endsubroutine initialize

   ! private methods
   elemental subroutine compute_u23(self)
   !< Compute interstates velocity.
   class(riemann_solver_compressible_pvl), intent(inout) :: self !< Solver.

   self%u23 = 0.5_R8P * (self%u_1 + self%u_4) - 2.0_R8P * (self%p_4 - self%p_1) / ((self%r_1 + self%r_4) * (self%a_1 + self%a_4))
   endsubroutine compute_u23

   elemental subroutine compute_up23(self)
   !< Compute interstates velocity and pressure.
   class(riemann_solver_compressible_pvl), intent(inout) :: self !< Solver.
   real(R8P)                                             :: ram  !< Mean value of `r * a`.

   ram = 0.25_R8P * (self%r_1 + self%r_4) * (self%a_1 + self%a_4)
   self%u23 = 0.5_R8P * ((self%u_1 + self%u_4) - (self%p_4 - self%p_1) / ram)
   self%p23 = 0.5_R8P * ((self%p_1 + self%p_4) - (self%u_4 - self%u_1) * ram)
   endsubroutine compute_up23

   pure subroutine compute_waves_u23(self, eos_left, state_left, eos_right, state_right, normal, waves)
   !< Compute waves speed `u23` algorithm.
   !<
   !< Use Primitive Variables Linearization algorithm by means of only `u23` approximation.
   class(riemann_solver_compressible_pvl), intent(inout) :: self        !< Solver.
   class(eos_object),                      intent(in)    :: eos_left    !< Equation of state for left state.
   class(conservative_object),             intent(in)    :: state_left  !< Left Riemann state.
   class(eos_object),                      intent(in)    :: eos_right   !< Equation of state for right state.
   class(conservative_object),             intent(in)    :: state_right !< Right Riemann state.
   type(vector),                           intent(in)    :: normal      !< Normal (versor) of face where fluxes are given.
   real(R8P),                              intent(out)   :: waves(1:)   !< Waves pattern.
   real(R8P)                                             :: x           !< Dummy variable.

   call self%set_states14(eos_left=eos_left, state_left=state_left, eos_right=eos_right, state_right=state_right, normal=normal)

   call self%compute_u23

   ! compute left state
   if (self%u23 < self%u_1) then ! shock
     x  = 0.25_R8P * (self%eos_1%g() + 1._R8P) * (self%u23 - self%u_1) / self%a_1
     self%s_1 = self%u_1 + self%a_1 * (x - sqrt(1.0_R8P + x * x))
   else ! rarefaction
     self%s_1 = self%u_1 - self%a_1
   endif
   ! compute right state
   if (self%u23 > self%u_4) then ! shock
     x  = 0.25_R8P * (self%eos_4%g() + 1._R8P) * (self%u23 - self%u_4) / self%a_4
     self%s_4 = self%u_4 + self%a_4 * (x + sqrt(1.0_R8P + x * x))
   else ! rarefaction
     self%s_4 = self%u_4  + self%a_4
   endif

   waves(1) = self%s_1
   waves(2) = self%s_1
   waves(3) = self%u23
   waves(4) = self%s_4
   waves(5) = self%s_4
   endsubroutine compute_waves_u23

   pure subroutine compute_waves_up23(self, eos_left, state_left, eos_right, state_right, normal, waves)
   !< Compute waves speed `u23` algorithm.
   !<
   !< Use Primitive Variables Linearization algorithm by means of only `up23` approximation.
   class(riemann_solver_compressible_pvl), intent(inout) :: self        !< Solver.
   class(eos_object),                      intent(in)    :: eos_left    !< Equation of state for left state.
   class(conservative_object),             intent(in)    :: state_left  !< Left Riemann state.
   class(eos_object),                      intent(in)    :: eos_right   !< Equation of state for right state.
   class(conservative_object),             intent(in)    :: state_right !< Right Riemann state.
   type(vector),                           intent(in)    :: normal      !< Normal (versor) of face where fluxes are given.
   real(R8P),                              intent(out)   :: waves(1:)   !< Waves pattern.

   call self%set_states14(eos_left=eos_left, state_left=state_left, eos_right=eos_right, state_right=state_right, normal=normal)

   call self%compute_up23

   ! compute left state
   if (self%u23 < self%u_1) then ! shock
      self%s_1 = self%u_1 - self%a_1 * sqrt(1._R8P + 0.5_R8P * (self%eos_1%g() + 1._R8P) / &
                                            self%eos_1%g() * (self%p23 / self%p_1 - 1._R8P))
   else ! rarefaction
      self%s_1 = self%u_1 - self%a_1
   endif
   ! compute right state
   if (self%u23 > self%u_4) then ! shock
      self%s_4 = self%u_4 + self%a_4 * sqrt(1._R8P + 0.5_R8P * (self%eos_4%g() + 1._R8P) / &
                                            self%eos_4%g() * (self%p23 / self%p_4 - 1._R8P))
   else ! rarefaction
      self%s_4 = self%u_4  + self%a_4
   endif

   waves(1) = self%s_1
   waves(2) = self%s_1
   waves(3) = self%u23
   waves(4) = self%s_4
   waves(5) = self%s_4
   endsubroutine compute_waves_up23

   pure subroutine solve_u23(self, eos_left, state_left, eos_right, state_right, normal, fluxes)
   !< Solve Riemann problem by `u23` algorithm.
   !<
   !< Use Primitive Variables Linearization algorithm by means of only `u23` approximation.
   class(riemann_solver_compressible_pvl), intent(inout) :: self        !< Solver.
   class(eos_object),                      intent(in)    :: eos_left    !< Equation of state for left state.
   class(conservative_object),             intent(in)    :: state_left  !< Left Riemann state.
   class(eos_object),                      intent(in)    :: eos_right   !< Equation of state for right state.
   class(conservative_object),             intent(in)    :: state_right !< Right Riemann state.
   type(vector),                           intent(in)    :: normal      !< Normal (versor) of face where fluxes are given.
   class(conservative_object),             intent(inout) :: fluxes      !< Fluxes of the Riemann Problem solution.
   real(R8P)                                             :: waves(1:5)  !< Waves pattern.

   call self%compute_waves_u23(eos_left=eos_left, state_left=state_left, eos_right=eos_right, state_right=state_right, &
                               normal=normal, waves=waves)
   call self%compute_fluxes(eos_left=eos_left, eos_right=eos_right, normal=normal, fluxes=fluxes)
   endsubroutine solve_u23

   pure subroutine solve_up23(self, eos_left, state_left, eos_right, state_right, normal, fluxes)
   !< Solve Riemann problem `up23` algorithm.
   !<
   !< Use Primitive Variables Linearization algorithm by means of only `up23` approximation.
   class(riemann_solver_compressible_pvl), intent(inout) :: self        !< Riemann pattern.
   class(eos_object),                      intent(in)    :: eos_left    !< Equation of state for left state.
   class(conservative_object),             intent(in)    :: state_left  !< Left Riemann state.
   class(eos_object),                      intent(in)    :: eos_right   !< Equation of state for right state.
   class(conservative_object),             intent(in)    :: state_right !< Right Riemann state.
   type(vector),                           intent(in)    :: normal      !< Normal (versor) of face where fluxes are given.
   class(conservative_object),             intent(inout) :: fluxes      !< Fluxes of the Riemann Problem solution.
   real(R8P)                                             :: waves(1:5)  !< Waves pattern.

   call self%compute_waves_up23(eos_left=eos_left, state_left=state_left, eos_right=eos_right, state_right=state_right, &
                                normal=normal, waves=waves)
   call self%compute_fluxes(eos_left=eos_left, eos_right=eos_right, normal=normal, fluxes=fluxes)
   endsubroutine solve_up23
endmodule foreseer_riemann_solver_compressible_pvl
