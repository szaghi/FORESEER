!< Define the Roe (with the Harten-Hyman entropy fix) Riemann solver of FORESEER library.

module foreseer_riemann_solver_compressible_roe
!< Define the Roe (with the Harten-Hyman entropy fix) Riemann solver of FORESEER library.

use flow_conservative_compressible, only : conservative_compressible, conservative_compressible_pointer
use flow_conservative_object, only : conservative_object
use flow_eos_compressible, only : eos_compressible
use flow_eos_object, only : eos_object
use foreseer_riemann_pattern_compressible_pvl, only : riemann_pattern_compressible_pvl
use foreseer_riemann_solver_object, only : riemann_solver_object
use penf, only : R8P
use vecfor, only : vector

implicit none
private
public :: riemann_solver_compressible_roe

type, extends(riemann_solver_object) :: riemann_solver_compressible_roe
   !< Roe (with the Harten-Hyman entropy fix) Riemann Solver.
   !<
   !< @note This is the implemention for [[conservative_compressible]] Riemann states.
   contains
      ! public deferred methods
      procedure, pass(self) :: initialize !< Initialize solver.
      procedure, pass(self) :: solve      !< Solve Riemann Problem.
      ! private methods
      procedure, nopass, private :: compute_roe_state !< Compute intermediate state.
endtype riemann_solver_compressible_roe

contains
   ! public deferred methods
   subroutine initialize(self, config)
   !< Initialize solver.
   class(riemann_solver_compressible_roe), intent(inout)        :: self    !< Solver.
   character(len=*),                       intent(in), optional :: config  !< Configuration for solver algorithm.
   character(len=:), allocatable                                :: config_ !< Configuration for solver algorithm, local variable.

   config_ = 'up23' ; if (present(config)) config_ = config
   ! call self%solver_pvl%initialize(config=config_)
   endsubroutine initialize

   subroutine compute_roe_state(eos_left, state_left, eos_right, state_right, r_d, r_u, r_e, r_a)
   !< Evaluate the intermediate state from the known states U1,U4 using the Roe linearization.
   class(eos_object),                      intent(in)    :: eos_left     !< Equation of state for left state.
   class(conservative_object),             intent(in)    :: state_left   !< Left Riemann state.
   class(eos_object),                      intent(in)    :: eos_right    !< Equation of state for right state.
   class(conservative_object),             intent(in)    :: state_right  !< Right Riemann state.
   real(R8P),                              intent(out)   :: r_d          !< Roe intermediate state density.
   type(vector),                           intent(out)   :: r_u          !< Roe intermediate state velocity vector..
   real(R8P),                              intent(out)   :: r_e          !< Roe intermediate state enthalpy.
   real(R8P),                              intent(out)   :: r_a          !< Roe intermediate state sound speed.
   real(R8P)                                             :: x, omx       !< x = sqrt(r1)/(sqrt(r1)+sqrt(r4)),  omx = 1-x
   real(R8P)                                             :: cp, cv       !< Roe intermediate state Cp and Cv.
   type(conservative_compressible), pointer              :: state_left_  !< Left Riemann state, local variable.
   type(conservative_compressible), pointer              :: state_right_ !< Right Riemann state, local variable.

   state_left_ => conservative_compressible_pointer(to=state_left)
   state_right_ => conservative_compressible_pointer(to=state_right)
   x = sqrt(state_left_%density)/(sqrt(state_left_%density)+sqrt(state_right_%density)) ; omx = 1._R8P - x

   r_d   = sqrt(state_left_%density*state_right_%density)
   r_u%x = state_left_%momentum%x / state_left_%density * x + state_right_%momentum%x / state_right_%density * omx
   r_u%y = state_left_%momentum%y / state_left_%density * x + state_right_%momentum%y / state_right_%density * omx
   r_u%z = state_left_%momentum%z / state_left_%density * x + state_right_%momentum%z / state_right_%density * omx
   r_e   = (state_left_%energy  + state_left_%pressure(eos_left))   / state_left_%density  * x + &
           (state_right_%energy + state_right_%pressure(eos_right)) / state_right_%density * omx
   cp    = eos_left%cp() * x + eos_right%cp() * omx
   cv    = eos_left%cv() * x + eos_right%cv() * omx
   r_a   = sqrt((cp/cv - 1._R8P)*(r_e - 0.5_R8P * r_u%sq_norm()))

   endsubroutine compute_roe_state

   subroutine solve(self, eos_left, state_left, eos_right, state_right, normal, fluxes)
   !< Solve Riemann Problem.
   !<
   !< Approximate Riemann Solver based on Roe (with Harten-Hyman entropy fix) algorithm.
   class(riemann_solver_compressible_roe), intent(in)    :: self         !< Solver.
   class(eos_object),                      intent(in)    :: eos_left     !< Equation of state for left state.
   class(conservative_object),             intent(in)    :: state_left   !< Left Riemann state.
   class(eos_object),                      intent(in)    :: eos_right    !< Equation of state for right state.
   class(conservative_object),             intent(in)    :: state_right  !< Right Riemann state.
   type(vector),                           intent(in)    :: normal       !< Normal (versor) of face where fluxes are given.
   class(conservative_object),             intent(inout) :: fluxes       !< Fluxes of the Riemann Problem solution.
   type(conservative_compressible), pointer              :: state_left_  !< Left Riemann state, local variable.
   type(conservative_compressible), pointer              :: state_right_ !< Right Riemann state, local variable.
   type(conservative_compressible)                       :: fluxes_left  !< Fluxes of left state.
   type(conservative_compressible)                       :: fluxes_right !< Fluxes of right state.
   type(riemann_pattern_compressible_pvl)                :: pattern      !< Riemann (states) PVL pattern solution.
   type(conservative_compressible)                       :: state23      !< Intermediate states.
   real(R8P)                                             :: r_d          !< Roe intermediate state density.
   type(vector)                                          :: r_u          !< Roe intermediate state velocity vector..
   real(R8P)                                             :: r_e          !< Roe intermediate state enthalpy.
   real(R8P)                                             :: r_a          !< Roe intermediate state sound speed.
   real(R8P)                                             :: waves(1:5)   !< Waves speed pattern.
   real(R8P)                                             :: lmax         !< Maximum wave speed estimation.
   type(vector)                                          :: vec_a        !< Vector of sound speeds, local variable.
   type(vector)                                          :: vec_r        !< Vector of densities, local variable.
   type(vector)                                          :: vec_m        !< Vector of momentums, local variable.
   type(vector)                                          :: vec_e        !< Vector of energies, local variable.
   real(R8P)                                             :: Dr           !< Density difference  Dr = r4-r1.
   real(R8P)                                             :: Du           !< Velocity difference Du = u4-u1.
   real(R8P)                                             :: Dp           !< Pressure difference Dp = p4-p1.
   real(R8P)                                             :: aa1,aa2,aa3  !< Wawes amplitudes Roe's estimation.
   real(R8P)                                             :: ll1,ll2,ll3  !< Wawes speeds Roe's estimation.
   real(R8P)                                             :: ls1,    ls3  !< Wawes speeds Roe's estimation with entropy fix of Harten-Hyman.
   real(R8P)                                             :: p_2, p_3     !< Pessure of state 2 and 3.
   real(R8P)                                             :: u23          !< Maximum wave speed estimation.

   call pattern%initialize(eos_left=eos_left, state_left=state_left, eos_right=eos_right, state_right=state_right, normal=normal)
   call pattern%compute_waves
   associate(r_1=>pattern%r_1, u_1=>pattern%u_1, p_1=>pattern%p_1, s_1=>pattern%s_1, &
             r_4=>pattern%r_4, u_4=>pattern%u_4, p_4=>pattern%p_4, s_4=>pattern%s_4, &
             s_2=>pattern%s_2, s_3=>pattern%s_3, u23=>pattern%u23)
     call pattern%compute_states23_from_u23(p_2=p_2, p_3=p_3)

     select case(minloc([-s_1, s_1 * s_2, s_2 * u23, u23 * s_3, s_3 * s_4, s_4], dim=1))
     case(1)
       call state_left%compute_fluxes(eos=eos_left, normal=normal, fluxes=fluxes_left)
     case(2)
       call compute_roe_state(eos_left =eos_left,  state_left =state_left,  &
                              eos_right=eos_right, state_right=state_right, &
                              r_d=r_d, r_u=r_u, r_e=r_e, r_a=r_a)
       Du  = u_4 - u_1
       Dp  = p_4 - p_1
       aa1 = 0.5_R8P*(Dp - r_d * r_a * Du) / (r_a * r_a)
       ll1 = (r_u .dot. normal) - r_a
       ls1 = s_1 * (s_2 - ll1) / (s_2 - s_1)
       call state_left%compute_fluxes(eos=eos_left, normal=normal, fluxes=fluxes_left)
       state23%density  = aa1
       state23%momentum = aa1 * ll1 * normal
       state23%energy   = aa1 * (r_e - (r_u .dot. normal) * r_a)
       select type(fluxes)
        type is(conservative_compressible)
#ifdef __GFORTRAN__
           fluxes = fluxes + ls1 * state23
#else
            error stop 'error: Intel fortran still does not support abstract math!'
#endif
       endselect
     case(3,4)
       call compute_roe_state(eos_left =eos_left,  state_left =state_left,  &
                              eos_right=eos_right, state_right=state_right, &
                              r_d=r_d, r_u=r_u, r_e=r_e, r_a=r_a)
       Dr  = r_4 - r_1
       Du  = u_4 - u_1
       Dp  = p_4 - p_1
       aa1 = 0.5_R8P * (Dp - r_d * r_a * Du) / (r_a * r_a)
       aa2 = Dr - Dp / (r_a * r_a)
       aa3 = 0.5_R8P * (Dp + r_d * r_a * Du) / (r_a * r_a)
       ll1 = (r_u .dot. normal) - r_a
       ll2 = (r_u .dot. normal)
       ll3 = (r_u .dot. normal) + r_a
       call state_left%compute_fluxes (eos=eos_left,  normal=normal, fluxes=fluxes_left)
       call state_right%compute_fluxes(eos=eos_right, normal=normal, fluxes=fluxes_right)
       vec_a%x = abs(ll1)                              ; vec_a%y = abs(ll2)                 ; vec_a%z = abs(ll3)
       vec_r%x = aa1                                   ; vec_r%y = aa2                      ; vec_r%z = aa3
       vec_m%x = aa1 * ll1                             ; vec_m%y = aa2 * ll2                ; vec_m%z = aa3 * ll3
       vec_e%x = aa1 * (r_e - (r_u .dot. normal) * r_a); vec_e%y = aa2 * 0.5_R8P * ll2 * ll2
       vec_e%x = aa3 * (r_e - (r_u .dot. normal) * r_a)
       state23%density  =  vec_a .dot. vec_r
       state23%momentum = (vec_a .dot. vec_m) * normal
       state23%energy   =  vec_a .dot. vec_e
       select type(fluxes)
        type is(conservative_compressible)
#ifdef __GFORTRAN__
           fluxes = 0.5_R8P * (fluxes_left + fluxes_right - state23)
#else
            error stop 'error: Intel fortran still does not support abstract math!'
#endif
       endselect
     case(5)
       call compute_roe_state(eos_left =eos_left,  state_left =state_left,  &
                              eos_right=eos_right, state_right=state_right, &
                              r_d=r_d, r_u=r_u, r_e=r_e, r_a=r_a)
       Du  = u_4 - u_1
       Dp  = p_4 - p_1
       aa3 = 0.5_R8P*(Dp - r_d * r_a * Du) / (r_a * r_a)
       ll3 = (r_u .dot. normal) - r_a
       ls3 = s_4 * (ll3 - s_3) / (s_4 - s_3)
       call state_right%compute_fluxes(eos=eos_right, normal=normal, fluxes=fluxes_right)
       state23%density  = aa3
       state23%momentum = aa3 * ll3 * normal
       state23%energy   = aa3 * (r_e - (r_u .dot. normal) * r_a)
       select type(fluxes)
        type is(conservative_compressible)
#ifdef __GFORTRAN__
           fluxes = fluxes + ls3 * state23
#else
            error stop 'error: Intel fortran still does not support abstract math!'
#endif
       endselect
     case(6)
       call state_right%compute_fluxes(eos=eos_right, normal=normal, fluxes=fluxes_right)
     endselect
   endassociate

   endsubroutine solve
endmodule foreseer_riemann_solver_compressible_roe
