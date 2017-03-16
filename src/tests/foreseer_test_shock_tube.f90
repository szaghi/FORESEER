!< FORESEER test: shock tube tester, 1D Euler equation.

module foreseer_euler_1d
!< Definition of Euler 1D class for FORESEER test.

use foreseer, only : conservative_object, conservative_compressible, primitive_compressible,         &
                     conservative_to_primitive_compressible, primitive_to_conservative_compressible, &
                     eos_object, eos_compressible,                                                   &
                     riemann_solver_object, riemann_solver_compressible_llf, riemann_solver_compressible_pvl
use penf, only : I4P, R8P
use foodie, only : integrand
use vecfor, only : ex, vector
use wenoof, only : interpolator_object, wenoof_create
! use wenoof, only : weno_factory, weno_constructor_upwind, weno_interpolator, weno_interpolator_upwind

implicit none
private
public :: euler_1d

type, extends(integrand) :: euler_1d
   !< Euler 1D PDEs system field.
   !<
   !< It is a FOODIE integrand class concrete extension.
   !<
   !<### 1D Euler PDEs system
   !< The 1D Euler PDEs system considered is a non linear, hyperbolic (inviscid) system of conservation laws for compressible gas
   !< dynamics, that reads as
   !<$$
   !<\begin{matrix}
   !<U_t = R(U)  \Leftrightarrow U_t = F(U)_x \\
   !<U = \begin{bmatrix}
   !<\rho \\
   !<\rho u \\
   !<\rho E
   !<\end{bmatrix}\;\;\;
   !<F(U) = \begin{bmatrix}
   !<\rho u \\
   !<\rho u^2 + p \\
   !<\rho u H
   !<\end{bmatrix}
   !<\end{matrix}
   !<$$
   !< where \(\rho\) is the density, \(u\) is the velocity, \(p\) the pressure, \(E\) the total internal specific energy and \(H\)
   !< the total specific enthalpy. The PDEs system must completed with the proper initial and boundary conditions. Moreover, an
   !< ideal (thermally and calorically perfect) gas is considered
   !<$$
   !<\begin{matrix}
   !<R = c_p - c_v \\
   !<\gamma = \frac{c_p}{c_v}\\
   !<e = c_v T \\
   !<h = c_p T
   !<\end{matrix}
   !<$$
   !< where *R* is the gas constant, \(c_p\,c_v\) are the specific heats at constant pressure and volume (respectively), *e* is the
   !< internal energy, *h* is the internal enthalpy and *T* is the temperature. The following addition equations of state hold:
   !<$$
   !<\begin{matrix}
   !<T = \frac{p}{\rho R} \\
   !<E = \rho e + \frac{1}{2} \rho u^2 \\
   !<H = \rho h + \frac{1}{2} \rho u^2 \\
   !<a = \sqrt{\frac{\gamma p}{\rho}}
   !<\end{matrix}
   !<$$
   !<
   !<#### Numerical grid organization
   !< The finite volume, Godunov's like approach is employed. The conservative variables (and the primitive ones) are co-located at
   !< the cell center. The cell and (inter)faces numeration is as follow.
   !<```
   !<                cell            (inter)faces
   !<                 |                   |
   !<                 v                   v
   !<     |-------|-------|-.....-|-------|-------|-------|-------|-.....-|-------|-------|-------|-.....-|-------|-------|
   !<     | 1-Ng  | 2-Ng  | ..... |  -1   |   0   |   1   |  2    | ..... |  Ni   | Ni+1  | Ni+1  | ..... |Ni+Ng-1| Ni+Ng |
   !<     |-------|-------|-.....-|-------|-------|-------|-------|-.....-|-------|-------|-------|-.....-|-------|-------|
   !<    0-Ng                             -1      0       1       2      Ni-1     Ni                                    Ni+Ng
   !<```
   !< Where *Ni* are the finite volumes (cells) used for discretizing the domain and *Ng* are the ghost cells used for imposing the
   !< left and right boundary conditions (for a total of *2Ng* cells).
   integer(I4P)                                 :: ord=0        !< Space accuracy formal order.
   integer(I4P)                                 :: Ni=0         !< Space dimension.
   integer(I4P)                                 :: Ng=0         !< Number of ghost cells for boundary conditions handling.
   real(R8P)                                    :: Dx=0._R8P    !< Space step.
   type(eos_compressible)                       :: eos          !< Equation of state.
   type(conservative_compressible), allocatable :: U(:)         !< Integrand (state) variables, whole physical domain.
   character(:),                    allocatable :: BC_L         !< Left boundary condition type.
   character(:),                    allocatable :: BC_R         !< Right boundary condition type.
   class(interpolator_object),      allocatable :: interpolator !< WENO interpolator.
   contains
      ! auxiliary methods
      procedure, pass(self) :: initialize       !< Initialize field.
      procedure, pass(self) :: destroy          !< Destroy field.
      procedure, pass(self) :: output           !< Extract Euler field.
      procedure, pass(self) :: dt => compute_dt !< Compute the current time step, by means of CFL condition.
      ! ADT integrand deferred methods
      procedure, pass(self) :: t => dEuler_dt                                       !< Time derivative, residuals function.
      procedure, pass(lhs)  :: local_error => euler_local_error                     !< Operator `||euler-euler||`.
      procedure, pass(lhs)  :: integrand_multiply_integrand => euler_multiply_euler !< Operator `*`.
      procedure, pass(lhs)  :: integrand_multiply_real => euler_multiply_real       !< Operator `euler * real`.
      procedure, pass(rhs)  :: real_multiply_integrand => real_multiply_euler       !< Operator `real * euler`.
      procedure, pass(lhs)  :: add => add_euler                                     !< Operator `+`.
      procedure, pass(lhs)  :: sub => sub_euler                                     !< Operator `-`.
      procedure, pass(lhs)  :: assign_integrand => euler_assign_euler               !< Operator `=`.
      procedure, pass(lhs)  :: assign_real => euler_assign_real                     !< Operator `euler = real`.
      ! private methods
      procedure, pass(self), private :: impose_boundary_conditions    !< Impose boundary conditions.
      procedure, pass(self), private :: reconstruct_interfaces_states !< Reconstruct interfaces states.
      ! procedure, pass(self), private :: riemann_solver                !< Solve the Riemann Problem at cell interfaces.
endtype euler_1d

contains
   ! auxiliary methods
   subroutine initialize(self, Ni, Dx, BC_L, BC_R, initial_state, eos, ord)
   !< Initialize field.
   class(euler_1d),              intent(inout)        :: self              !< Euler field.
   integer(I4P),                 intent(in)           :: Ni                !< Space dimension.
   real(R8P),                    intent(in)           :: Dx                !< Space step.
   character(*),                 intent(in)           :: BC_L              !< Left boundary condition type.
   character(*),                 intent(in)           :: BC_R              !< Right boundary condition type.
   type(primitive_compressible), intent(in)           :: initial_state(1:) !< Initial state of primitive variables.
   type(eos_compressible),       intent(in)           :: eos               !< Equation of state.
   integer(I4P),                 intent(in), optional :: ord               !< Space accuracy formal order.
   integer(I4P)                                       :: i                 !< Space counter.

   call self%destroy
   self%ord = 1 ; if (present(ord)) self%ord = ord
   self%Ng = (self%ord + 1) / 2
   if (self%ord>1) then
      call wenoof_create(interpolator_type='reconstructor-JS', &
                         S=self%Ng,                            &
                         interpolator=self%interpolator)
   endif
   self%Ni = Ni
   self%Dx = Dx
   if (allocated(self%U)) deallocate(self%U) ; allocate(self%U(1-self%Ng:self%Ni+self%Ng))
   do i=1, Ni
      self%U(i) = primitive_to_conservative_compressible(primitive=initial_state(i), eos=eos)
   enddo
   self%eos = eos
   self%BC_L = BC_L
   self%BC_R = BC_R
   endsubroutine initialize

   pure subroutine destroy(self)
   !< Destroy field.
   class(euler_1d), intent(inout) :: self !< Euler field.

   self%ord = 0
   self%Ni = 0
   self%Ng = 0
   self%Dx = 0._R8P
   if (allocated(self%U)) deallocate(self%U)
   if (allocated(self%BC_L)) deallocate(self%BC_L)
   if (allocated(self%BC_R)) deallocate(self%BC_R)
   if (allocated(self%interpolator)) deallocate(self%interpolator)
   endsubroutine destroy

   pure function output(self, is_primitive) result(state)
   !< Output the Euler field state.
   class(euler_1d), intent(in)           :: self          !< Euler field.
   logical,         intent(in), optional :: is_primitive  !< Output in primitive variables.
   real(R8P), allocatable                :: state(:,:)    !< Euler state vector.
   logical                               :: is_primitive_ !< Output in primitive variables, local variable.
   type(primitive_compressible)          :: primitive     !< Primitive state.
   integer(I4P)                          :: i             !< Counter.

   is_primitive_ = .false. ; if (present(is_primitive)) is_primitive_ = is_primitive
   allocate(state(1:5, 1:self%Ni))
   if (is_primitive_) then
      do i=1, self%Ni
         primitive = conservative_to_primitive_compressible(conservative=self%U(i), eos=self%eos)
      enddo
   else
      do i=1, self%Ni
         state(:, i) = self%U(i)%array()
      enddo
   endif
   endfunction output

   pure function compute_dt(self, Nmax, Tmax, t, CFL) result(Dt)
   !< Compute the current time step by means of CFL condition.
   class(euler_1d), intent(in) :: self !< Euler field.
   integer(I4P),    intent(in) :: Nmax !< Maximun number of iterates.
   real(R8P),       intent(in) :: Tmax !< Maximum time (ignored if Nmax>0).
   real(R8P),       intent(in) :: t    !< Time.
   real(R8P),       intent(in) :: CFL  !< CFL value.
   real(R8P)                   :: Dt   !< Time step.
   type(vector)                :: u    !< Velocity vector.
   real(R8P)                   :: a    !< Speed of sound.
   real(R8P)                   :: vmax !< Maximum propagation speed of signals.
   integer(I4P)                :: i    !< Counter.

   associate(Ni=>self%Ni, Dx=>self%Dx)
      vmax = 0._R8P
      do i=1, Ni
         u = self%U(i)%velocity()
         a = self%eos%speed_of_sound(density=self%U(i)%density, pressure=self%U(i)%pressure(eos=self%eos))
         vmax = max(vmax, u%normL2() + a)
      enddo
      Dt = Dx * CFL / vmax
      if (Nmax <= 0) then
         if ((t + Dt) > Tmax) Dt = Tmax - t
      endif
   endassociate
   endfunction compute_dt

   ! ADT integrand deferred methods
   function dEuler_dt(self, t) result(dState_dt)
   !< Time derivative of Euler field, the residuals function.
   class(euler_1d), intent(in)           :: self                         !< Euler field.
   real(R8P),       intent(in), optional :: t                            !< Time.
   class(integrand), allocatable         :: dState_dt                    !< Euler field time derivative.
   type(primitive_compressible)          :: P(1-self%Ng:self%Ni+self%Ng) !< Primitive variables.
   type(primitive_compressible)          :: PR(1:2,0:self%Ni+1)          !< Reconstructed primitive variables.
   type(conservative_compressible)       :: UR(1:2,0:self%Ni+1)          !< Reconstructed conservative variables.
   type(conservative_compressible)       :: F(0:self%Ni)                 !< Fluxes of conservative variables.
   integer(I4P)                          :: i                            !< Counter.
   type(riemann_solver_compressible_llf) :: riemann_solver               !< Riemann solver.

   do i=1, self%Ni
      P(i) = conservative_to_primitive_compressible(conservative=self%U(i), eos=self%eos)
   enddo
   call self%impose_boundary_conditions(P=P)
   call self%reconstruct_interfaces_states(primitive=P, r_primitive=PR)
   do i=0, self%Ni+1
      UR(1, i) = primitive_to_conservative_compressible(primitive=PR(1, i), eos=self%eos)
      UR(2, i) = primitive_to_conservative_compressible(primitive=PR(2, i), eos=self%eos)
   enddo
   do i=0, self%Ni
      call riemann_solver%solve(eos_left=self%eos,  state_left=UR( 2, i  ), &
                                eos_right=self%eos, state_right=UR(1, i+1), normal=ex, fluxes=F(i))
   enddo
   allocate(euler_1d :: dState_dt)
   select type(dState_dt)
   class is(euler_1d)
      dState_dt = self
      do i=1, self%Ni
          dState_dt%U(i) = (F(i - 1) - F(i)) / self%Dx
      enddo
   endselect
   endfunction dEuler_dt

  function euler_local_error(lhs, rhs) result(error)
  !< Estimate local truncation error between 2 euler approximations.
  !<
  !< The estimation is done by norm L2 of U:
  !<
  !< $$ error = \sqrt{ \sum_i{\sum_i{ \frac{(lhs\%U_i - rhs\%U_i)^2}{lhs\%U_i^2} }} } $$
  class(euler_1d),  intent(in) :: lhs      !< Left hand side.
  class(integrand), intent(in) :: rhs      !< Right hand side.
  real(R8P)                    :: error    !< Error estimation.
  real(R8P), allocatable       :: U_lhs(:) !< Serialized conservative variables.
  real(R8P), allocatable       :: U_rhs(:) !< Serialized conservative variables.
  integer(I4P)                 :: i        !< Space counter.
  integer(I4P)                 :: v        !< Variables counter.

  select type(rhs)
  class is (euler_1d)
     error = 0._R8P
     do i=1, lhs%Ni
       U_lhs = lhs%U(i)%array()
       U_rhs = rhs%U(i)%array()
       do v=1, size(U_lhs, dim=1)
         error = error + (U_lhs(v) - U_rhs(v)) ** 2 / U_lhs(v) ** 2
       enddo
     enddo
     error = sqrt(error)
  endselect
  endfunction euler_local_error

  function euler_multiply_euler(lhs, rhs) result(opr)
  !< Multiply an Euler field by another one.
  class(euler_1d),  intent(in)  :: lhs !< Left hand side.
  class(integrand), intent(in)  :: rhs !< Right hand side.
  class(integrand), allocatable :: opr !< Operator result.
  integer(I4P)                  :: i   !< Counter.

  allocate(euler_1d :: opr)
  select type(opr)
  class is(euler_1d)
     opr = lhs
  endselect
  select type(opr)
  class is(euler_1d)
     select type(rhs)
     class is (euler_1d)
        do i=1, lhs%Ni
           opr%U(i) = lhs%U(i) * rhs%U(i)
        enddo
     endselect
  endselect
  endfunction euler_multiply_euler

  function euler_multiply_real(lhs, rhs) result(opr)
  !< Multiply an Euler field by a real scalar.
  class(euler_1d), intent(in)   :: lhs !< Left hand side.
  real(R8P),       intent(in)   :: rhs !< Right hand side.
  class(integrand), allocatable :: opr !< Operator result.
  integer(I4P)                  :: i   !< Counter.

  allocate(euler_1d :: opr)
  select type(opr)
  class is(euler_1d)
     opr = lhs
  endselect
  select type(opr)
  class is(euler_1d)
     do i=1, lhs%Ni
        opr%U(i) = rhs * lhs%U(i)
     enddo
  endselect
  endfunction euler_multiply_real

  function real_multiply_euler(lhs, rhs) result(opr)
  !< Multiply a real scalar by an Euler field.
  real(R8P),       intent(in)   :: lhs !< Left hand side.
  class(euler_1d), intent(in)   :: rhs !< Right hand side.
  class(integrand), allocatable :: opr !< Operator result.
  integer(I4P)                  :: i   !< Counter.

  allocate(euler_1d :: opr)
  select type(opr)
  class is(euler_1d)
     opr = rhs
  endselect
  select type(opr)
  class is(euler_1d)
     do i=1, rhs%Ni
        opr%U(i) = lhs * rhs%U(i)
     enddo
  endselect
  endfunction real_multiply_euler

  function add_euler(lhs, rhs) result(opr)
  !< Add two Euler fields.
  class(euler_1d),  intent(in)  :: lhs !< Left hand side.
  class(integrand), intent(in)  :: rhs !< Right hand side.
  class(integrand), allocatable :: opr !< Operator result.
  integer(I4P)                  :: i   !< Counter.

  allocate (euler_1d :: opr)
  select type(opr)
  class is(euler_1d)
     opr = lhs
  endselect
  select type(opr)
  class is(euler_1d)
     select type(rhs)
     class is (euler_1d)
        do i=1, lhs%Ni
           opr%U(i) = lhs%U(i) + rhs%U(i)
        enddo
     endselect
  endselect
  endfunction add_euler

  function sub_euler(lhs, rhs) result(opr)
  !< Subtract two Euler fields.
  class(euler_1d),  intent(in)  :: lhs !< Left hand side.
  class(integrand), intent(in)  :: rhs !< Right hand side.
  class(integrand), allocatable :: opr !< Operator result.
  integer(I4P)                  :: i   !< Counter.

  allocate (euler_1d :: opr)
  select type(opr)
  class is(euler_1d)
     opr = lhs
  endselect
  select type(opr)
  class is(euler_1d)
     select type(rhs)
     class is (euler_1d)
        do i=1, lhs%Ni
           opr%U(i) = lhs%U(i) - rhs%U(i)
        enddo
     endselect
  endselect
  endfunction sub_euler

  subroutine euler_assign_euler(lhs, rhs)
  !< Assign one Euler field to another.
  class(euler_1d),  intent(inout) :: lhs !< Left hand side.
  class(integrand), intent(in)    :: rhs !< Right hand side.
  integer(I4P)                    :: i   !< Counter.

  select type(rhs)
  class is(euler_1d)
     lhs%ord  = rhs%ord
     lhs%Ni   = rhs%Ni
     lhs%Ng   = rhs%Ng
     lhs%Dx   = rhs%Dx
     lhs%eos  = rhs%eos
     if (allocated(rhs%U)) then
        if (allocated(lhs%U)) deallocate(lhs%U) ; allocate(lhs%U(1:lhs%Ni))
        select type(rhs)
        class is(euler_1d)
           if (allocated(rhs%U)) then
              do i=1, lhs%Ni
                 lhs%U(i) = rhs%U(i)
              enddo
           endif
        endselect
     endif
     if (allocated(rhs%BC_L)) lhs%BC_L = rhs%BC_L
     if (allocated(rhs%BC_R)) lhs%BC_R = rhs%BC_R
     if (allocated(rhs%interpolator)) then
        if (allocated(lhs%interpolator)) deallocate(lhs%interpolator)
        allocate(lhs%interpolator, source=rhs%interpolator)
     endif
  endselect
  endsubroutine euler_assign_euler

  subroutine euler_assign_real(lhs, rhs)
  !< Assign one real to an Euler field.
  class(euler_1d), intent(inout) :: lhs !< Left hand side.
  real(R8P),       intent(in)    :: rhs !< Right hand side.
  integer(I4P)                   :: i   !< Counter.

  if (allocated(lhs%U)) then
     do i=1, lhs%Ni
        lhs%U(i)%density = rhs
        lhs%U(i)%momentum = rhs
        lhs%U(i)%energy = rhs
     enddo
  endif
  endsubroutine euler_assign_real

   ! private methods
   pure subroutine impose_boundary_conditions(self, P)
   !< Impose boundary conditions.
   !<
   !< The boundary conditions are imposed on the primitive variables by means of the ghost cells approach.
   class(euler_1d),              intent(in)    :: self          !< Euler field.
   type(primitive_compressible), intent(inout) :: P(1-self%Ng:) !< Primitive variables.
   integer(I4P)                                :: i             !< Space counter.

   select case(trim(adjustl(self%BC_L)))
      case('TRA') ! trasmissive (non reflective) BC
         do i=1-self%Ng, 0
            P(i) = P(-i+1)
         enddo
      case('REF') ! reflective BC
         do i=1-self%Ng, 0
            P(i)%density  =   P(-i+1)%density
            P(i)%velocity = - P(-i+1)%velocity
            P(i)%pressure =   P(-i+1)%pressure
         enddo
   endselect

   select case(trim(adjustl(self%BC_R)))
      case('TRA') ! trasmissive (non reflective) BC
         do i=self%Ni+1, self%Ni+self%Ng
            P(i) = P(self%Ni-(i-self%Ni-1))
         enddo
      case('REF') ! reflective BC
         do i=self%Ni+1, self%Ni+self%Ng
            P(i)%density  =   P(self%Ni-(i-self%Ni-1))%density
            P(i)%velocity = - P(self%Ni-(i-self%Ni-1))%velocity
            P(i)%pressure =   P(self%Ni-(i-self%Ni-1))%pressure
         enddo
   endselect
   endsubroutine impose_boundary_conditions

  subroutine reconstruct_interfaces_states(self, primitive, r_primitive)
  !< Reconstruct the interfaces states (into primitive variables formulation) by the requested order of accuracy.
  class(euler_1d),              intent(in)    :: self                                 !< Euler field.
  type(primitive_compressible), intent(in)    :: primitive(1-self%Ng:self%Ni+self%Ng) !< Primitive variables.
  type(primitive_compressible), intent(inout) :: r_primitive(1:2, 0:self%Ni+1)        !< Reconstructed primitive variables.
  type(primitive_compressible)                :: Pm(1:2)                              !< Mean of primitive variables.
  real(R8P)                                   :: LPm(1:3, 1:3, 1:2)                   !< Mean left eigenvectors matrix.
  real(R8P)                                   :: RPm(1:3, 1:3, 1:2)                   !< Mean right eigenvectors matrix.
  real(R8P)                                   :: C(1:2, 1-self%Ng:-1+self%Ng, 1:3)    !< Pseudo characteristic variables.
  real(R8P)                                   :: CR(1:2, 1:3)                         !< Pseudo characteristic reconst. vars.
  real(R8P)                                   :: buffer(1:3)                          !< Dummy buffer.
  integer(I4P)                                :: i                                    !< Counter.
  integer(I4P)                                :: j                                    !< Counter.
  integer(I4P)                                :: f                                    !< Counter.
  integer(I4P)                                :: v                                    !< Counter.
  class(interpolator_object), allocatable     :: interpolator                         !< WENO interpolator.

  select case(self%ord)
  case(1) ! 1st order piecewise constant reconstruction
     do i=0, self%Ni+1
        r_primitive(1, i) = primitive(i)
        r_primitive(2, i) = r_primitive(1, i)
     enddo
  case(3, 5, 7) ! 3rd, 5th or 7th order WENO reconstruction
     call wenoof_create(interpolator_type='reconstructor-JS', &
                        S=self%Ng,                            &
                        interpolator=interpolator)
     do i=0, self%Ni+1
        ! trasform primitive variables to pseudo charteristic ones
        do f=1, 2
           Pm(f) = 0.5_R8P * (primitive(i+f-2) + primitive(i+f-1))
        enddo
        do f=1, 2
           LPm(:, :, f) = Pm(f)%left_eigenvectors(eos=self%eos)
           RPm(:, :, f) = Pm(f)%right_eigenvectors(eos=self%eos)
        enddo
        do j=i+1-self%Ng, i-1+self%Ng
           do f=1, 2
              do v=1, 3
                 C(f, j-i, v) = dot_product(LPm(v, :, f), [primitive(j)%density, primitive(j)%velocity%x, primitive(j)%pressure])
              enddo
           enddo
        enddo
        ! compute WENO reconstruction of pseudo charteristic variables
        do v=1, 3
           call interpolator%interpolate(stencil=C(:, :, v), interpolation=CR(:, v))
        enddo
        ! trasform back reconstructed pseudo charteristic variables to primitive ones
        do f=1, 2
           do v=1, 3
              buffer(v) = dot_product(RPm(v, :, f), CR(f, :))
           enddo
           r_primitive(f, i)%density  = buffer(1)
           r_primitive(f, i)%velocity = buffer(2) * ex
           r_primitive(f, i)%pressure = buffer(3)
        enddo
     enddo
  endselect
  endsubroutine reconstruct_interfaces_states
endmodule foreseer_euler_1d

program foreseer_test_shock_tube
!< FORESEER test: shock tube tester, 1D Euler equation.

use flap, only : command_line_interface
use foodie, only : tvd_runge_kutta_integrator
use foreseer, only : conservative_compressible, primitive_compressible,         &
                     conservative_to_primitive_compressible, primitive_to_conservative_compressible, &
                     eos_compressible
use foreseer_euler_1d, only : euler_1d
use penf, only : FR8P, I4P, R8P, str
! use pyplot_module, only :  pyplot
use vecfor, only : ex, vector

implicit none
type(command_line_interface)     :: cli                   !< Command line interface handler.
type(tvd_runge_kutta_integrator) :: rk_integrator         !< Runge-Kutta integrator.
integer, parameter               :: rk_stages=1           !< Runge-Kutta stages number.
type(euler_1d)                   :: rk_stage(1:rk_stages) !< Runge-Kutta stages.
real(R8P)                        :: dt                    !< Time step.
real(R8P)                        :: t                     !< Time.
integer(I4P)                     :: steps                 !< Time steps counter.
type(euler_1d)                   :: domain                !< Domain of Euler equations.
real(R8P), parameter             :: CFL=0.7_R8P           !< CFL value.
character(3)                     :: BC_L                  !< Left boundary condition type.
character(3)                     :: BC_R                  !< Right boundary condition type.
integer(I4P)                     :: Ni                    !< Number of grid cells.
real(R8P)                        :: Dx                    !< Space step discretization.
real(R8P), allocatable           :: x(:)                  !< Cell center x-abscissa values.
integer(I4P)                     :: error                 !< Error handler.
integer(I4P)                     :: steps_max             !< Maximum number of time steps.
logical                          :: plots                 !< Flag for activating plots saving.
logical                          :: results               !< Flag for activating results saving.
logical                          :: time_serie            !< Flag for activating time serie-results saving.
logical                          :: verbose               !< Flag for activating more verbose output.
! real(R_P), allocatable           :: final_state(:,:)      !< Final state.

call cli%init(progname    = 'foreseer_test_shock_tube',                             &
              authors     = 'S. Zaghi',                                             &
              license     = 'GNU GPLv3',                                            &
              description = 'FORESEER test: shock tube tester, 1D Euler equation.', &
              examples    = ["foreseer_test_shock_tube --results  ",                &
                             "foreseer_test_shock_tube -r -t -v -p",                &
                             "foreseer_test_shock_tube            ",                &
                             "foreseer_test_shock_tube --plots -r "])
call cli%add(switch='--Ni', help='Number finite volumes used', required=.false., act='store', def='100', error=error)
call cli%add(switch='--steps', help='Number time steps performed', required=.false., act='store', def='60', error=error)
call cli%add(switch='--results', switch_ab='-r', help='Save results', required=.false., act='store_true', def='.false.', &
             error=error)
call cli%add(switch='--plots', switch_ab='-p', help='Save plots of results', required=.false., act='store_true', def='.false.', &
             error=error)
call cli%add(switch='--tserie', switch_ab='-t', help='Save time-serie-results', required=.false., act='store_true', def='.false.', &
             error=error)
call cli%add(switch='--verbose', help='Verbose output', required=.false., act='store_true', def='.false.', error=error)
call cli%parse(error=error)
call cli%get(switch='--Ni', val=Ni, error=error) ; if (error/=0) stop
call cli%get(switch='--steps', val=steps_max, error=error) ; if (error/=0) stop
call cli%get(switch='-r', val=results, error=error) ; if (error/=0) stop
call cli%get(switch='-p', val=plots, error=error) ; if (error/=0) stop
call cli%get(switch='-t', val=time_serie, error=error) ; if (error/=0) stop
call cli%get(switch='--verbose', val=verbose, error=error) ; if (error/=0) stop

call initialize
call save_time_serie(filename='euler_1D_integration-tvdrk-'//trim(str(n=rk_stages, no_sign=.true.))//'-time_serie.dat', t=t)
do steps=1, steps_max
   if (verbose) print "(A)", '    Time step: '//str(n=dt)//', Time: '//str(n=t)
   dt = domain%dt(Nmax=steps_max, Tmax=0._R8P, t=t, CFL=CFL)
   call rk_integrator%integrate(U=domain, stage=rk_stage, dt=dt, t=t)
   t = t + dt
   call save_time_serie(t=t)
enddo
if (verbose) print "(A)", '    Time step: '//str(n=dt)//', Time: '//str(n=t)
call save_time_serie(t=t, finish=.true.)

contains
   subroutine initialize()
   !< Initialize the test.
   type(primitive_compressible), allocatable :: initial_state(:) !< Initial state of primitive variables.
   integer(I4P)                              :: i                !< Space counter.

   call rk_integrator%init(stages=rk_stages)
   t = 0._R8P
   allocate(x(1:Ni))
   allocate(initial_state(1:Ni))
   Dx = 1._R8P / Ni
   ! Sod's problem
   BC_L = 'TRA'
   BC_R = 'TRA'
   do i=1, Ni / 2
      x(i) = Dx * i - 0.5_R8P * Dx
      initial_state(i)%density  = 1._R8P
      initial_state(i)%velocity = 0._R8P
      initial_state(i)%pressure = 1._R8P
   enddo
   do i=Ni / 2 + 1, Ni
      x(i) = Dx * i - 0.5_R8P * Dx
      initial_state(i)%density  = 0.125_R8P
      initial_state(i)%velocity = 0._R8P
      initial_state(i)%pressure = 0.1_R8P
   enddo
   call domain%initialize(Ni=Ni, Dx=Dx, BC_L=BC_L, BC_R=BC_R, &
                          initial_state=initial_state, eos=eos_compressible(cp=1040.004_R8P, cv=742.86_R8P), ord=1)
   endsubroutine initialize

   subroutine save_time_serie(filename, finish, t)
   !< Save time-serie results.
   character(*), intent(in), optional :: filename  !< Output filename.
   logical,      intent(in), optional :: finish    !< Flag for triggering the file closing.
   real(R8P),    intent(in)           :: t         !< Current integration time.
   integer(I4P), save                 :: tsfile    !< File unit for saving time serie results.
   type(primitive_compressible)       :: primitive !< Primitive variables.
   integer(I4P)                       :: i         !< Counter.

   if (time_serie) then
      if (present(filename)) then
         open(newunit=tsfile, file=filename)
      endif
      write(tsfile, '(A)')'VARIABLES = "x" "rho" "u" "p"'
      write(tsfile, '(A)')'ZONE T="'//str(n=t)//'"'
      do i=1, Ni
         primitive = conservative_to_primitive_compressible(conservative=domain%U(i), eos=domain%eos)
         write(tsfile, '(4'//'('//FR8P//',1X))')x(i), primitive%density, primitive%velocity%x, primitive%pressure
      enddo
      if (present(finish)) then
         if (finish) close(tsfile)
      endif
   endif
   endsubroutine save_time_serie
endprogram foreseer_test_shock_tube
