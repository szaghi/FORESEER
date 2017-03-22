!< Define the Roe (with the Harten-Hyman entropy fix) Riemann solver of FORESEER library.

module foreseer_riemann_solver_compressible_roe
!< Define the Roe (with the Harten-Hyman entropy fix) Riemann solver of FORESEER library.

use foreseer_conservative_compressible, only : conservative_compressible, conservative_compressible_pointer
use foreseer_conservative_object, only : conservative_object
use foreseer_eos_compressible, only : eos_compressible, eos_compressible_pointer
use foreseer_eos_object, only : eos_object
use foreseer_riemann_solver_compressible_object, only : riemann_solver_compressible_object
use foreseer_riemann_solver_compressible_pvl, only : riemann_solver_compressible_pvl
use foreseer_riemann_solver_object, only : riemann_solver_object
use penf, only : R8P
use vecfor, only : vector

implicit none
private
public :: riemann_solver_compressible_roe

type, extends(riemann_solver_compressible_object) :: riemann_solver_compressible_roe
   !< Roe (with the Harten-Hyman entropy fix) Riemann Solver.
   !<
   !< @note This is the implemention for [[conservative_compressible]] Riemann states.
   type(riemann_solver_compressible_pvl) :: solver_pvl !< PVL Riemann solver.
   contains
      ! public deferred methods
      procedure, pass(self) :: compute_waves       !< Compute waves pattern.
      procedure, pass(self) :: initialize          !< Initialize solver.
      procedure, pass(self) :: solve               !< Solve Riemann Problem.
      procedure, nopass     :: compute_roe_state   !< Compute intermediate state.
endtype riemann_solver_compressible_roe

contains
   ! public deferred methods
   pure subroutine compute_waves(self, eos_left, state_left, eos_right, state_right, normal, waves)
   !< Compute waves pattern.
   !<
   !< The PVL approximation is based on a 3 waves pattern where the acoustic waves are reduced to a single linear wave instead
   !< of a non linear fan one.
   class(riemann_solver_compressible_roe), intent(inout) :: self        !< Solver.
   class(eos_object),                      intent(in)    :: eos_left    !< Equation of state for left state.
   class(conservative_object),             intent(in)    :: state_left  !< Left Riemann state.
   class(eos_object),                      intent(in)    :: eos_right   !< Equation of state for right state.
   class(conservative_object),             intent(in)    :: state_right !< Right Riemann state.
   type(vector),                           intent(in)    :: normal      !< Normal (versor) of face where fluxes are given.
   real(R8P),                              intent(out)   :: waves(1:)   !< Waves pattern.

   call self%solver_pvl%compute_waves(eos_left=eos_left, state_left=state_left, eos_right=eos_right, state_right=state_right, &
                                      normal=normal, waves=waves)
   endsubroutine compute_waves

   subroutine initialize(self, config)
   !< Initialize solver.
   class(riemann_solver_compressible_roe), intent(inout)        :: self   !< Solver.
   character(len=*),                       intent(in), optional :: config !< Configuration for solver algorithm.

   call self%solver_pvl%initialize(config='u23')
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
   type(eos_compressible),          pointer              :: eos_left_    !< Left Riemann state, local variable.
   type(eos_compressible),          pointer              :: eos_right_   !< Right Riemann state, local variable.

   state_left_ => conservative_compressible_pointer(to=state_left)
   state_right_ => conservative_compressible_pointer(to=state_right)
   eos_left_ => eos_compressible_pointer(to=eos_left)
   eos_right_ => eos_compressible_pointer(to=eos_right)
   x = sqrt(state_left_%density)/(sqrt(state_left_%density)+sqrt(state_right_%density)) ; omx = 1._R8P - x

   r_d   = sqrt(state_left_%density*state_right_%density)
   r_u%x = state_left_%momentum%x / state_left_%density * x + state_right_%momentum%x / state_right_%density * omx
   r_u%y = state_left_%momentum%y / state_left_%density * x + state_right_%momentum%y / state_right_%density * omx
   r_u%z = state_left_%momentum%z / state_left_%density * x + state_right_%momentum%z / state_right_%density * omx
   r_e   = (state_left_%energy  + state_left_%pressure(eos_left))   / state_left_%density  * x + &
           (state_right_%energy + state_right_%pressure(eos_right)) / state_right_%density * omx
   cp    = eos_left_%cp_ * x + eos_right_%cp_ * omx
   cv    = eos_left_%cv_ * x + eos_right_%cv_ * omx
   r_a   = sqrt((cp/cv - 1._R8P)*(r_e - 0.5_R8P * r_u%sq_norm()))

   endsubroutine compute_roe_state

   subroutine solve(self, eos_left, state_left, eos_right, state_right, normal, fluxes)
   !< Solve Riemann Problem.
   !<
   !< Approximate Riemann Solver based on Roe (with Harten-Hyman entropy fix) algorithm.
   class(riemann_solver_compressible_roe), intent(inout) :: self         !< Solver.
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

   call self%set_states14(eos_left=eos_left, state_left=state_left, eos_right=eos_right, state_right=state_right, normal=normal)
   call self%solver_pvl%compute_waves(eos_left=eos_left, state_left=state_left, eos_right=eos_right, state_right=state_right, &
                                      normal=normal, waves=waves)
   associate(r_1=>self%solver_pvl%r_1, u_1=>self%solver_pvl%u_1, p_1=>self%solver_pvl%p_1, &
             r_4=>self%solver_pvl%r_4, u_4=>self%solver_pvl%u_4, p_4=>self%solver_pvl%p_4, &
             s_1=>self%s_1, s_2=>self%s_2, s_3=>self%s_3, s_4=>self%s_4, u23=>self%u23)
     s_1 = waves(1)
     u23 = waves(3)
     s_4 = waves(5)
     call self%compute_states23_from_u23(p_2=p_2, p_3=p_3)

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
           fluxes = fluxes + ls1 * state23
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
       ll1 = (r_u .dot. normal) - r_d
       ll2 = (r_u .dot. normal)
       ll3 = (r_u .dot. normal) + r_d
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
           fluxes = 0.5_R8P * (fluxes_left + fluxes_right - state23)
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
           fluxes = fluxes + ls3 * state23
       endselect
     case(6)
       call state_right%compute_fluxes(eos=eos_right, normal=normal, fluxes=fluxes_right)
     endselect
   endassociate

   endsubroutine solve
endmodule foreseer_riemann_solver_compressible_roe
