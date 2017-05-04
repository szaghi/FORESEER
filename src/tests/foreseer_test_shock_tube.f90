!< FORESEER test: shock tube tester, 1D Euler equation.

module foreseer_euler_1d
!< Definition of Euler 1D class for FORESEER test.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use foreseer, only : conservative_object, conservative_compressible, primitive_compressible,         &
                     conservative_to_primitive_compressible, primitive_to_conservative_compressible, &
                     eos_object, eos_compressible,                                                   &
                     riemann_solver_object, riemann_solver_compressible_exact,                       &
                     riemann_solver_compressible_hllc, riemann_solver_compressible_llf,              &
                     riemann_solver_compressible_pvl, riemann_solver_compressible_roe
use penf, only : I4P, R_P
use foodie, only : integrand_object
use vecfor, only : ex, vector
use wenoof, only : interpolator_object, wenoof_create

implicit none
private
public :: euler_1d

type, extends(integrand_object) :: euler_1d
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
   integer(I4P)                                 :: Ni=0                                  !< Space dimension.
   integer(I4P)                                 :: Nc=0                                  !< Number of conservative variables.
   integer(I4P)                                 :: No=0                                  !< Dimension of ODE system.
   integer(I4P)                                 :: Ng=0                                  !< Ghost cells number.
   real(R_P)                                    :: Dx=0._R_P                             !< Space step.
   type(eos_compressible)                       :: eos                                   !< Equation of state.
   integer(I4P)                                 :: weno_order=0                          !< WENO reconstruction order.
   character(:),                    allocatable :: weno_scheme                           !< WENO scheme.
   character(:),                    allocatable :: weno_variables                        !< WENO reconstruction variables.
   character(:),                    allocatable :: riemann_solver_scheme                 !< Riemann solver scheme.
   type(conservative_compressible), allocatable :: U(:)                                  !< Integrand (state) variables.
   character(:),                    allocatable :: BC_L                                  !< Left boundary condition type.
   character(:),                    allocatable :: BC_R                                  !< Right boundary condition type.
   class(interpolator_object),      allocatable :: interpolator                          !< WENO interpolator.
   class(riemann_solver_object),    allocatable :: riemann_solver                        !< Riemann solver.
   procedure(reconstruct_interfaces_), pointer  :: reconstruct_interfaces=>&
                                                   reconstruct_interfaces_characteristic !< Reconstruct interface states.
   contains
      ! auxiliary methods
      procedure, pass(self) :: initialize       !< Initialize field.
      procedure, pass(self) :: destroy          !< Destroy field.
      procedure, pass(self) :: dt => compute_dt !< Compute the current time step, by means of CFL condition.
      ! ADT integrand deferred methods
      procedure, pass(self) :: t => dEuler_dt !< Time derivative, residuals function.
      ! operators
      procedure, pass(lhs) :: local_error => euler_local_error !< Operator `||euler-euler||`.
      ! +
      procedure, pass(lhs) :: integrand_add_integrand => euler_add_euler !< `+` operator.
      procedure, pass(lhs) :: integrand_add_real => euler_add_real       !< `+ real` operator.
      procedure, pass(rhs) :: real_add_integrand => real_add_euler       !< `real +` operator.
      ! *
      procedure, pass(lhs) :: integrand_multiply_integrand => euler_multiply_euler         !< `*` operator.
      procedure, pass(lhs) :: integrand_multiply_real => euler_multiply_real               !< `* real` operator.
      procedure, pass(rhs) :: real_multiply_integrand => real_multiply_euler               !< `real *` operator.
      procedure, pass(lhs) :: integrand_multiply_real_scalar => euler_multiply_real_scalar !< `* real_scalar` operator.
      procedure, pass(rhs) :: real_scalar_multiply_integrand => real_scalar_multiply_euler !< `real_scalar *` operator.
      ! -
      procedure, pass(lhs) :: integrand_sub_integrand => euler_sub_euler !< `-` operator.
      procedure, pass(lhs) :: integrand_sub_real => euler_sub_real       !< `- real` operator.
      procedure, pass(rhs) :: real_sub_integrand => real_sub_euler       !< `real -` operator.
      ! =
      procedure, pass(lhs) :: assign_integrand => euler_assign_euler !< `=` operator.
      procedure, pass(lhs) :: assign_real => euler_assign_real       !< `= real` operator.
      ! private methods
      procedure, pass(self), private :: impose_boundary_conditions            !< Impose boundary conditions.
      procedure, pass(self), private :: reconstruct_interfaces_characteristic !< Reconstruct (charc.) interface states.
      procedure, pass(self), private :: reconstruct_interfaces_conservative   !< Reconstruct (cons.) interface states.
      procedure, pass(self), private :: reconstruct_interfaces_primitive      !< Reconstruct (prim.) interface states.
endtype euler_1d

abstract interface
   !< Abstract interfaces of [[euler_1d]] pointer methods.
   subroutine reconstruct_interfaces_(self, conservative, r_conservative)
   !< Reconstruct interface states.
   import :: conservative_compressible, euler_1d, primitive_compressible
   class(euler_1d),                 intent(in)    :: self                     !< Euler field.
   type(conservative_compressible), intent(in)    :: conservative(1-self%Ng:) !< Conservative variables.
   type(conservative_compressible), intent(inout) :: r_conservative(1:, 0:)   !< Reconstructed conservative variables.
   endsubroutine reconstruct_interfaces_
endinterface

contains
   ! auxiliary methods
   subroutine initialize(self, Ni, Dx, BC_L, BC_R, initial_state, eos, weno_scheme, weno_order, weno_variables, &
                         riemann_solver_scheme)
   !< Initialize field.
   class(euler_1d),              intent(inout)        :: self                   !< Euler field.
   integer(I4P),                 intent(in)           :: Ni                     !< Space dimension.
   real(R_P),                    intent(in)           :: Dx                     !< Space step.
   character(*),                 intent(in)           :: BC_L                   !< Left boundary condition type.
   character(*),                 intent(in)           :: BC_R                   !< Right boundary condition type.
   type(primitive_compressible), intent(in)           :: initial_state(1:)      !< Initial state of primitive variables.
   type(eos_compressible),       intent(in)           :: eos                    !< Equation of state.
   character(*),                 intent(in), optional :: weno_scheme            !< WENO scheme.
   integer(I4P),                 intent(in), optional :: weno_order             !< WENO reconstruction order.
   character(*),                 intent(in), optional :: weno_variables         !< Variables on which WENO reconstruction is done.
   character(*),                 intent(in), optional :: riemann_solver_scheme  !< Riemann solver scheme.
   integer(I4P)                                       :: i                      !< Space couner.

   call self%destroy

   self%weno_scheme = 'reconstructor-JS' ; if (present(weno_scheme)) self%weno_scheme = weno_scheme
   self%weno_order = 1 ; if (present(weno_order)) self%weno_order = weno_order
   self%weno_variables = 'characteristic' ; if (present(weno_variables)) self%weno_variables = trim(adjustl(weno_variables))

   self%Ng = (self%weno_order + 1) / 2
   if (self%weno_order>1) call wenoof_create(interpolator_type=self%weno_scheme, S=self%Ng, interpolator=self%interpolator)

   self%riemann_solver_scheme = 'llf'
   if (present(riemann_solver_scheme)) self%riemann_solver_scheme = trim(adjustl(riemann_solver_scheme))

   self%Ni = Ni
   self%Dx = Dx

   self%eos = eos
   if (allocated(self%U)) deallocate(self%U) ; allocate(self%U(1-self%Ng:self%Ni+self%Ng))
   do i=1, Ni
      self%U(i) = primitive_to_conservative_compressible(primitive=initial_state(i), eos=eos)
   enddo

   self%Nc = size(self%U(1)%array(), dim=1)
   self%No = self%Ni * self%Nc

   self%BC_L = BC_L
   self%BC_R = BC_R

   select case(self%weno_variables)
   case('characteristic')
      self%reconstruct_interfaces => reconstruct_interfaces_characteristic
   case('conservative')
      self%reconstruct_interfaces => reconstruct_interfaces_conservative
   case('primitive')
      self%reconstruct_interfaces => reconstruct_interfaces_primitive
   case default
      write(stderr, '(A)') 'error: WENO reconstruction variables set "'//self%weno_variables//'" unknown!'
      stop
   endselect

   select case(self%riemann_solver_scheme)
   case('exact')
      allocate(riemann_solver_compressible_exact :: self%riemann_solver)
   case('hllc')
      allocate(riemann_solver_compressible_hllc :: self%riemann_solver)
   case('llf')
      allocate(riemann_solver_compressible_llf :: self%riemann_solver)
   case('pvl')
      allocate(riemann_solver_compressible_pvl :: self%riemann_solver)
   case('roe')
      allocate(riemann_solver_compressible_roe :: self%riemann_solver)
   case default
      write(stderr, '(A)') 'error: Riemann Solver scheme "'//self%riemann_solver_scheme//'" unknown!'
      stop
   endselect
   endsubroutine initialize

   pure subroutine destroy(self)
   !< Destroy field.
   class(euler_1d), intent(inout) :: self !< Euler field.

   self%Ni = 0
   self%Nc = 0
   self%No = 0
   self%Ng = 0
   self%Dx = 0._R_P
   self%weno_order = 0
   if (allocated(self%weno_scheme)) deallocate(self%weno_scheme)
   if (allocated(self%weno_variables)) deallocate(self%weno_variables)
   if (allocated(self%riemann_solver_scheme)) deallocate(self%riemann_solver_scheme)
   if (allocated(self%U)) deallocate(self%U)
   if (allocated(self%BC_L)) deallocate(self%BC_L)
   if (allocated(self%BC_R)) deallocate(self%BC_R)
   if (allocated(self%interpolator)) deallocate(self%interpolator)
   if (allocated(self%riemann_solver)) deallocate(self%riemann_solver)
   self%reconstruct_interfaces => reconstruct_interfaces_characteristic
   endsubroutine destroy

   pure function compute_dt(self, steps_max, t_max, t, CFL) result(Dt)
   !< Compute the current time step by means of CFL condition.
   class(euler_1d), intent(in) :: self      !< Euler field.
   integer(I4P),    intent(in) :: steps_max !< Maximun number of time steps.
   real(R_P),       intent(in) :: t_max     !< Maximum integration time.
   real(R_P),       intent(in) :: t         !< Time.
   real(R_P),       intent(in) :: CFL       !< CFL value.
   real(R_P)                   :: Dt        !< Time step.
   type(vector)                :: u         !< Velocity vector.
   real(R_P)                   :: a         !< Speed of sound.
   real(R_P)                   :: vmax      !< Maximum propagation speed of signals.
   integer(I4P)                :: i         !< Counter.

   associate(Ni=>self%Ni, Dx=>self%Dx)
      vmax = 0._R_P
      do i=1, Ni
         u = self%U(i)%velocity()
         a = self%eos%speed_of_sound(density=self%U(i)%density, pressure=self%U(i)%pressure(eos=self%eos))
         vmax = max(vmax, u%normL2() + a)
      enddo
      Dt = Dx * CFL / vmax
      if (steps_max <= 0 .and. t_max > 0._R_P) then
         if ((t + Dt) > t_max) Dt = t_max - t
      endif
   endassociate
   endfunction compute_dt

   ! integrand_object deferred methods
   function dEuler_dt(self, t) result(dState_dt)
   !< Time derivative of Euler field, the residuals function.
   class(euler_1d), intent(in)           :: self                         !< Euler field.
   real(R_P),       intent(in), optional :: t                            !< Time.
   real(R_P), allocatable                :: dState_dt(:)                 !< Euler field time derivative.
   type(conservative_compressible)       :: U(1-self%Ng:self%Ni+self%Ng) !< Conservative variables.
   type(conservative_compressible)       :: UR(1:2,0:self%Ni+1)          !< Reconstructed conservative variables.
   type(conservative_compressible)       :: F(0:self%Ni)                 !< Fluxes of conservative variables.
   integer(I4P)                          :: i                            !< Counter.

   do i=1, self%Ni
      U(i) = self%U(i)
   enddo
   call self%impose_boundary_conditions(U=U)
   call self%reconstruct_interfaces(conservative=U, r_conservative=UR)
   do i=0, self%Ni
      call self%riemann_solver%solve(eos_left=self%eos,  state_left=UR( 2, i  ), &
                                     eos_right=self%eos, state_right=UR(1, i+1), normal=ex, fluxes=F(i))
   enddo
   allocate(dState_dt(1:self%No))
   do i=1, self%Ni
      dState_dt(i + (i-1) * (self%Nc-1) : i + i * (self%Nc-1)) = (F(i - 1) - F(i)) / self%Dx
   enddo
   endfunction dEuler_dt

   ! operators
   function euler_local_error(lhs, rhs) result(error)
   !< Estimate local truncation error between 2 euler approximations.
   class(euler_1d),         intent(in) :: lhs   !< Left hand side.
   class(integrand_object), intent(in) :: rhs   !< Right hand side.
   real(R_P)                           :: error !< Error estimation.

   error stop 'error: local error is not definite for Euler field!'
   endfunction euler_local_error

   ! +
   pure function euler_add_euler(lhs, rhs) result(opr)
   !< `+` operator.
   class(euler_1d),         intent(in) :: lhs    !< Left hand side.
   class(integrand_object), intent(in) :: rhs    !< Right hand side.
   real(R_P), allocatable              :: opr(:) !< Operator result.
   integer(I4P)                        :: i      !< Counter.

   allocate(opr(1:lhs%No))
   select type(rhs)
   class is (euler_1d)
      do i=1, lhs%Ni
         opr(i + (i-1) * (lhs%Nc-1) : i + i * (lhs%Nc-1)) = lhs%U(i) + rhs%U(i)
      enddo
   endselect
   endfunction euler_add_euler

   pure function euler_add_real(lhs, rhs) result(opr)
   !< `+ real` operator.
   class(euler_1d),  intent(in) :: lhs     !< Left hand side.
   real(R_P),        intent(in) :: rhs(1:) !< Right hand side.
   real(R_P), allocatable       :: opr(:)  !< Operator result.
   integer(I4P)                 :: i       !< Counter.

   allocate(opr(1:lhs%No))
   do i=1, lhs%Ni
      opr(i + (i-1) * (lhs%Nc-1) : i + i * (lhs%Nc-1)) = lhs%U(i) + rhs(i + (i-1) * (lhs%Nc-1) : i + i * (lhs%Nc-1))
   enddo
   endfunction euler_add_real

   pure function real_add_euler(lhs, rhs) result(opr)
   !< `real +` operator.
   real(R_P),        intent(in) :: lhs(1:) !< Left hand side.
   class(euler_1d),  intent(in) :: rhs     !< Right hand side.
   real(R_P), allocatable       :: opr(:)  !< Operator result.
   integer(I4P)                 :: i       !< Counter.

   allocate(opr(1:rhs%No))
   do i=1, rhs%Ni
      opr(i + (i-1) * (rhs%Nc-1) : i + i * (rhs%Nc-1)) = lhs(i + (i-1) * (rhs%Nc-1) : i + i * (rhs%Nc-1)) + rhs%U(i)
   enddo
   endfunction real_add_euler

   ! *
   pure function euler_multiply_euler(lhs, rhs) result(opr)
   !< `*` operator.
   class(euler_1d),         intent(in) :: lhs    !< Left hand side.
   class(integrand_object), intent(in) :: rhs    !< Right hand side.
   real(R_P), allocatable              :: opr(:) !< Operator result.
   integer(I4P)                        :: i      !< Counter.

   allocate(opr(1:lhs%No))
   select type(rhs)
   class is (euler_1d)
      do i=1, lhs%Ni
         opr(i + (i-1) * (lhs%Nc-1) : i + i * (lhs%Nc-1)) = lhs%U(i) * rhs%U(i)
      enddo
   endselect
   endfunction euler_multiply_euler

   pure function euler_multiply_real(lhs, rhs) result(opr)
   !< `* real` operator.
   class(euler_1d), intent(in) :: lhs     !< Left hand side.
   real(R_P),       intent(in) :: rhs(1:) !< Right hand side.
   real(R_P), allocatable      :: opr(:)  !< Operator result.
   integer(I4P)                :: i       !< Counter.

   allocate(opr(1:lhs%No))
   do i=1, lhs%Ni
      opr(i + (i-1) * (lhs%Nc-1) : i + i * (lhs%Nc-1)) = lhs%U(i) * rhs(i + (i-1) * (lhs%Nc-1) : i + i * (lhs%Nc-1))
   enddo
   endfunction euler_multiply_real

   pure function real_multiply_euler(lhs, rhs) result(opr)
   !< `real *` operator.
   real(R_P),       intent(in) :: lhs(1:) !< Left hand side.
   class(euler_1d), intent(in) :: rhs     !< Right hand side.
   real(R_P), allocatable      :: opr(:)  !< Operator result.
   integer(I4P)                :: i       !< Counter.

   allocate(opr(1:rhs%No))
   do i=1, rhs%Ni
      opr(i + (i-1) * (rhs%Nc-1) : i + i * (rhs%Nc-1)) = lhs(i + (i-1) * (rhs%Nc-1) : i + i * (rhs%Nc-1)) * rhs%U(i)
   enddo
   endfunction real_multiply_euler

   pure function euler_multiply_real_scalar(lhs, rhs) result(opr)
   !< `* real_scalar` operator.
   class(euler_1d), intent(in) :: lhs    !< Left hand side.
   real(R_P),       intent(in) :: rhs    !< Right hand side.
   real(R_P), allocatable      :: opr(:) !< Operator result.
   integer(I4P)                :: i      !< Counter.

   allocate(opr(1:lhs%No))
   do i=1, lhs%Ni
      opr(i + (i-1) * (lhs%Nc-1) : i + i * (lhs%Nc-1)) = lhs%U(i) * rhs
   enddo
   endfunction euler_multiply_real_scalar

   pure function real_scalar_multiply_euler(lhs, rhs) result(opr)
   !< `real_scalar *` operator.
   real(R_P),       intent(in) :: lhs    !< Left hand side.
   class(euler_1d), intent(in) :: rhs    !< Right hand side.
   real(R_P), allocatable      :: opr(:) !< Operator result.
   integer(I4P)                :: i      !< Counter.

   allocate(opr(1:rhs%No))
   do i=1, rhs%Ni
      opr(i + (i-1) * (rhs%Nc-1) : i + i * (rhs%Nc-1)) = lhs * rhs%U(i)
   enddo
   endfunction real_scalar_multiply_euler

   ! -
   pure function euler_sub_euler(lhs, rhs) result(opr)
   !< `-` operator.
   class(euler_1d),         intent(in) :: lhs    !< Left hand side.
   class(integrand_object), intent(in) :: rhs    !< Right hand side.
   real(R_P), allocatable              :: opr(:) !< Operator result.
   integer(I4P)                        :: i      !< Counter.

   allocate(opr(1:lhs%No))
   select type(rhs)
   class is (euler_1d)
      do i=1, lhs%Ni
         opr(i + (i-1) * (lhs%Nc-1) : i + i * (lhs%Nc-1)) = lhs%U(i) - rhs%U(i)
      enddo
   endselect
   endfunction euler_sub_euler

   pure function euler_sub_real(lhs, rhs) result(opr)
   !< `- real` operator.
   class(euler_1d),  intent(in) :: lhs     !< Left hand side.
   real(R_P),        intent(in) :: rhs(1:) !< Right hand side.
   real(R_P), allocatable       :: opr(:)  !< Operator result.
   integer(I4P)                 :: i       !< Counter.

   allocate(opr(1:lhs%No))
   do i=1, lhs%Ni
      opr(i + (i-1) * (lhs%Nc-1) : i + i * (lhs%Nc-1)) = lhs%U(i) - rhs(i + (i-1) * (lhs%Nc-1) : i + i * (lhs%Nc-1))
   enddo
   endfunction euler_sub_real

   pure function real_sub_euler(lhs, rhs) result(opr)
   !< `real -` operator.
   real(R_P),        intent(in) :: lhs(1:) !< Left hand side.
   class(euler_1d),  intent(in) :: rhs     !< Right hand side.
   real(R_P), allocatable       :: opr(:)  !< Operator result.
   integer(I4P)                 :: i       !< Counter.

   allocate(opr(1:rhs%No))
   do i=1, rhs%Ni
      opr(i + (i-1) * (rhs%Nc-1) : i + i * (rhs%Nc-1)) = lhs(i + (i-1) * (rhs%Nc-1) : i + i * (rhs%Nc-1)) - rhs%U(i)
   enddo
   endfunction real_sub_euler

   ! =
   pure subroutine euler_assign_euler(lhs, rhs)
   !< `=` operator.
   class(euler_1d),         intent(inout) :: lhs !< Left hand side.
   class(integrand_object), intent(in)    :: rhs !< Right hand side.
   integer(I4P)                           :: i   !< Counter.

   select type(rhs)
   class is(euler_1d)
      lhs%Ni         = rhs%Ni
      lhs%Nc         = rhs%Nc
      lhs%No         = rhs%No
      lhs%Ng         = rhs%Ng
      lhs%Dx         = rhs%Dx
      lhs%eos        = rhs%eos
      lhs%weno_order = rhs%weno_order
      if (allocated(rhs%weno_scheme          )) lhs%weno_scheme           = rhs%weno_scheme
      if (allocated(rhs%weno_variables       )) lhs%weno_variables        = rhs%weno_variables
      if (allocated(rhs%riemann_solver_scheme)) lhs%riemann_solver_scheme = rhs%riemann_solver_scheme
      if (allocated(rhs%U)) then
         if (allocated(lhs%U)) then
            if (size(lhs%U) /= size(rhs%U)) deallocate(lhs%U)
         endif
         if (.not.allocated(lhs%U)) allocate(lhs%U(1:lhs%Ni))
         do i=1, lhs%Ni
            lhs%U(i) = rhs%U(i)
         enddo
      endif
      if (allocated(rhs%BC_L)) lhs%BC_L = rhs%BC_L
      if (allocated(rhs%BC_R)) lhs%BC_R = rhs%BC_R
      if (allocated(rhs%interpolator)) then
         if (allocated(lhs%interpolator)) deallocate(lhs%interpolator)
         allocate(lhs%interpolator, mold=rhs%interpolator)
         lhs%interpolator = rhs%interpolator
      endif
      if (associated(rhs%reconstruct_interfaces)) lhs%reconstruct_interfaces => rhs%reconstruct_interfaces
      if (allocated(rhs%riemann_solver)) then
         if (allocated(lhs%riemann_solver)) deallocate(lhs%riemann_solver)
         allocate(lhs%riemann_solver, mold=rhs%riemann_solver)
         lhs%riemann_solver = rhs%riemann_solver
      endif
      if (associated(rhs%reconstruct_interfaces)) lhs%reconstruct_interfaces => rhs%reconstruct_interfaces
   endselect
   endsubroutine euler_assign_euler

   pure subroutine euler_assign_real(lhs, rhs)
   !< `= real` operator.
   class(euler_1d), intent(inout) :: lhs     !< Left hand side.
   real(R_P),       intent(in)    :: rhs(1:) !< Right hand side.
   integer(I4P)                   :: i       !< Counter.

   do i=1, lhs%Ni
      lhs%U(i) = rhs(i + (i-1) * (lhs%Nc-1) : i + i * (lhs%Nc-1))
   enddo
   endsubroutine euler_assign_real

   ! private methods
   pure subroutine impose_boundary_conditions(self, U)
   !< Impose boundary conditions.
   !<
   !< The boundary conditions are imposed on the primitive variables by means of the ghost cells approach.
   class(euler_1d),                 intent(in)    :: self          !< Euler field.
   type(conservative_compressible), intent(inout) :: U(1-self%Ng:) !< Conservative variables.
   integer(I4P)                                   :: i             !< Space counter.

   select case(trim(adjustl(self%BC_L)))
      case('TRA') ! trasmissive (non reflective) BC
         do i=1-self%Ng, 0
            U(i) = U(-i+1)
         enddo
      case('REF') ! reflective BC
         do i=1-self%Ng, 0
            U(i)%density  =   U(-i+1)%density
            U(i)%momentum = - U(-i+1)%momentum
            U(i)%energy   =   U(-i+1)%energy
         enddo
   endselect

   select case(trim(adjustl(self%BC_R)))
      case('TRA') ! trasmissive (non reflective) BC
         do i=self%Ni+1, self%Ni+self%Ng
            U(i) = U(self%Ni-(i-self%Ni-1))
         enddo
      case('REF') ! reflective BC
         do i=self%Ni+1, self%Ni+self%Ng
            U(i)%density  =   U(self%Ni-(i-self%Ni-1))%density
            U(i)%momentum = - U(self%Ni-(i-self%Ni-1))%momentum
            U(i)%energy   =   U(self%Ni-(i-self%Ni-1))%energy
         enddo
   endselect
   endsubroutine impose_boundary_conditions

   subroutine reconstruct_interfaces_characteristic(self, conservative, r_conservative)
   !< Reconstruct interfaces states.
   !<
   !< The reconstruction is done in pseudo characteristic variables.
   class(euler_1d),                 intent(in)    :: self                                 !< Euler field.
   type(conservative_compressible), intent(in)    :: conservative(1-self%Ng:)             !< Conservative variables.
   type(conservative_compressible), intent(inout) :: r_conservative(1:, 0:)               !< Reconstructed conservative vars.
   type(primitive_compressible)                   :: primitive(1-self%Ng:self%Ni+self%Ng) !< Primitive variables.
   type(primitive_compressible)                   :: r_primitive(1:2, 0:self%Ni+1)        !< Reconstructed primitive variables.
   type(primitive_compressible)                   :: Pm(1:2)                              !< Mean of primitive variables.
   real(R_P)                                      :: LPm(1:3, 1:3, 1:2)                   !< Mean left eigenvectors matrix.
   real(R_P)                                      :: RPm(1:3, 1:3, 1:2)                   !< Mean right eigenvectors matrix.
   real(R_P)                                      :: C(1:2, 1-self%Ng:-1+self%Ng, 1:3)    !< Pseudo characteristic variables.
   real(R_P)                                      :: CR(1:2, 1:3)                         !< Pseudo characteristic reconst.
   real(R_P)                                      :: buffer(1:3)                          !< Dummy buffer.
   integer(I4P)                                   :: i                                    !< Counter.
   integer(I4P)                                   :: j                                    !< Counter.
   integer(I4P)                                   :: f                                    !< Counter.
   integer(I4P)                                   :: v                                    !< Counter.

   select case(self%weno_order)
   case(1) ! 1st order piecewise constant reconstruction
      do i=0, self%Ni+1
         r_conservative(1, i) = conservative(i)
         r_conservative(2, i) = r_conservative(1, i)
      enddo
   case(3, 5, 7, 9, 11, 13, 15, 17) ! 3rd-17th order WENO reconstruction
      do i=1-self%Ng, self%Ni+self%Ng
         primitive(i) = conservative_to_primitive_compressible(conservative=conservative(i), eos=self%eos)
      enddo
      do i=0, self%Ni+1
         ! compute pseudo charteristic variables
         do f=1, 2
            Pm(f) = 0.5_R_P * (primitive(i+f-2) + primitive(i+f-1))
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
            call self%interpolator%interpolate(stencil=C(:, :, v), interpolation=CR(:, v))
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
      do i=0, self%Ni+1
         r_conservative(1, i) = primitive_to_conservative_compressible(primitive=r_primitive(1, i), eos=self%eos)
         r_conservative(2, i) = primitive_to_conservative_compressible(primitive=r_primitive(2, i), eos=self%eos)
      enddo
   endselect
   endsubroutine reconstruct_interfaces_characteristic

   subroutine reconstruct_interfaces_conservative(self, conservative, r_conservative)
   !< Reconstruct interfaces states.
   !<
   !< The reconstruction is done in conservative variables.
   class(euler_1d),                 intent(in)    :: self                                 !< Euler field.
   type(conservative_compressible), intent(in)    :: conservative(1-self%Ng:)             !< Conservative variables.
   type(conservative_compressible), intent(inout) :: r_conservative(1:, 0:)               !< Reconstructed conservative vars.
   real(R_P), allocatable                         :: U(:)                                 !< Serialized conservative variables.
   real(R_P)                                      :: C(1:2, 1-self%Ng:-1+self%Ng, 1:3)    !< Pseudo characteristic variables.
   real(R_P)                                      :: CR(1:2, 1:3)                         !< Pseudo characteristic reconst.
   integer(I4P)                                   :: i                                    !< Counter.
   integer(I4P)                                   :: j                                    !< Counter.
   integer(I4P)                                   :: f                                    !< Counter.
   integer(I4P)                                   :: v                                    !< Counter.

   select case(self%weno_order)
   case(1) ! 1st order piecewise constant reconstruction
      do i=0, self%Ni+1
         r_conservative(1, i) = conservative(i)
         r_conservative(2, i) = r_conservative(1, i)
      enddo
   case(3, 5, 7, 9, 11, 13, 15, 17) ! 3rd-17th order WENO reconstruction
      do i=0, self%Ni+1
         do j=i+1-self%Ng, i-1+self%Ng
             U = conservative(j)%array()
            do f=1, 2
               C(f, j-i, 1) = U(1)
               C(f, j-i, 2) = U(2)
               C(f, j-i, 3) = U(5)
            enddo
         enddo
         do v=1, 3
            call self%interpolator%interpolate(stencil=C(:, :, v), interpolation=CR(:, v))
         enddo
         do f=1, 2
            r_conservative(f, i)%density  = CR(f, 1)
            r_conservative(f, i)%momentum = CR(f, 2) * ex
            r_conservative(f, i)%energy   = CR(f, 3)
         enddo
      enddo
   endselect
   endsubroutine reconstruct_interfaces_conservative

   subroutine reconstruct_interfaces_primitive(self, conservative, r_conservative)
   !< Reconstruct interfaces states.
   !<
   !< The reconstruction is done in primitive variables.
   class(euler_1d),                 intent(in)    :: self                                 !< Euler field.
   type(conservative_compressible), intent(in)    :: conservative(1-self%Ng:)             !< Conservative variables.
   type(conservative_compressible), intent(inout) :: r_conservative(1:, 0:)               !< Reconstructed conservative vars.
   type(primitive_compressible)                   :: primitive(1-self%Ng:self%Ni+self%Ng) !< Primitive variables.
   type(primitive_compressible)                   :: r_primitive(1:2, 0:self%Ni+1)        !< Reconstructed primitive variables.
   real(R_P), allocatable                         :: P(:)                                 !< Serialized primitive variables.
   real(R_P)                                      :: C(1:2, 1-self%Ng:-1+self%Ng, 1:3)    !< Pseudo characteristic variables.
   real(R_P)                                      :: CR(1:2, 1:3)                         !< Pseudo characteristic reconst.
   integer(I4P)                                   :: i                                    !< Counter.
   integer(I4P)                                   :: j                                    !< Counter.
   integer(I4P)                                   :: f                                    !< Counter.
   integer(I4P)                                   :: v                                    !< Counter.

   select case(self%weno_order)
   case(1) ! 1st order piecewise constant reconstruction
      do i=0, self%Ni+1
         r_conservative(1, i) = conservative(i)
         r_conservative(2, i) = r_conservative(1, i)
      enddo
   case(3, 5, 7, 9, 11, 13, 15, 17) ! 3rd-17th order WENO reconstruction
      do i=1-self%Ng, self%Ni+self%Ng
         primitive(i) = conservative_to_primitive_compressible(conservative=conservative(i), eos=self%eos)
      enddo
      do i=0, self%Ni+1
         do j=i+1-self%Ng, i-1+self%Ng
             P = primitive(j)%array()
            do f=1, 2
               C(f, j-i, 1) = P(1)
               C(f, j-i, 2) = P(2)
               C(f, j-i, 3) = P(5)
            enddo
         enddo
         do v=1, 3
            call self%interpolator%interpolate(stencil=C(:, :, v), interpolation=CR(:, v))
         enddo
         do f=1, 2
            r_primitive(f, i)%density  = CR(f, 1)
            r_primitive(f, i)%velocity = CR(f, 2) * ex
            r_primitive(f, i)%pressure = CR(f, 3)
         enddo
      enddo
      do i=0, self%Ni+1
         r_conservative(1, i) = primitive_to_conservative_compressible(primitive=r_primitive(1, i), eos=self%eos)
         r_conservative(2, i) = primitive_to_conservative_compressible(primitive=r_primitive(2, i), eos=self%eos)
      enddo
   endselect
   endsubroutine reconstruct_interfaces_primitive
endmodule foreseer_euler_1d

program foreseer_test_shock_tube
!< FORESEER test: shock tube tester, 1D Euler equation.

use flap, only : command_line_interface
use foodie, only : integrator_runge_kutta_ssp
use foreseer, only : conservative_compressible, primitive_compressible,         &
                     conservative_to_primitive_compressible, primitive_to_conservative_compressible, &
                     eos_compressible
use foreseer_euler_1d, only : euler_1d
use penf, only : cton, FR_P, I4P, R_P, str
use vecfor, only : ex, vector

implicit none
character(99), allocatable       :: riemann_solver_schemes(:)  !< Riemann Problem solver scheme(s).
character(99), allocatable       :: riemann_problems(:)        !< Riemann problems.
character(99)                    :: weno_scheme                !< WENO scheme.
integer(I4P)                     :: weno_order                 !< WENO accuracy order.
real(R_P)                        :: weno_eps                   !< WENO epsilon paramter.
character(99)                    :: weno_variables             !< WENO variables.
character(99)                    :: rk_scheme                  !< Runge-Kutta scheme.
type(integrator_runge_kutta_ssp) :: rk_integrator              !< Runge-Kutta integrator.
type(euler_1d), allocatable      :: rk_stage(:)                !< Runge-Kutta stages.
real(R_P)                        :: dt                         !< Time step.
real(R_P)                        :: t                          !< Time.
integer(I4P)                     :: step                       !< Time steps counter.
type(euler_1d)                   :: domain                     !< Domain of Euler equations.
real(R_P)                        :: CFL                        !< CFL value.
character(3)                     :: BC_L                       !< Left boundary condition type.
character(3)                     :: BC_R                       !< Right boundary condition type.
integer(I4P)                     :: Ni                         !< Number of grid cells.
real(R_P)                        :: Dx                         !< Space step discretization.
real(R_P), allocatable           :: x(:)                       !< Cell center x-abscissa values.
integer(I4P)                     :: steps_max                  !< Maximum number of time steps.
real(R_P)                        :: t_max                      !< Maximum integration time.
logical                          :: results                    !< Flag for activating results saving.
logical                          :: time_serie                 !< Flag for activating time serie-results saving.
logical                          :: verbose                    !< Flag for activating more verbose output.
integer(I4P)                     :: s                          !< Schemes counter.
integer(I4P)                     :: p                          !< Problems counter.
real(R_P), parameter             :: pi = 4._R_P * atan(1._R_P) !< Pi greek.

call parse_command_line_interface
do s=1, size(riemann_solver_schemes, dim=1)
   do p=1, size(riemann_problems, dim=1)
      if (verbose) print "(A)", 'Use Riemann Problem solver "'//trim(adjustl(riemann_solver_schemes(s)))//'"'
      call initialize(riemann_solver_scheme=riemann_solver_schemes(s), riemann_problem=riemann_problems(p))
      call save_time_serie(filename='foreseer_test_shock_tube-'//                                                 &
                                    trim(adjustl(riemann_problems(p)))//'-'//                                     &
                                    trim(adjustl(weno_scheme))//'-'//trim(str(weno_order, no_sign=.true.))//'-'// &
                                    trim(adjustl(weno_variables))//'-'//                                          &
                                    trim(adjustl(rk_scheme))//'-'//                                               &
                                    trim(adjustl(riemann_solver_schemes(s)))//'-'//                               &
                                    'Ni_'//trim(str(Ni, no_sign=.true.))//'.dat', t=t)
      step = 0
      time_loop: do
         step = step + 1
         dt = domain%dt(steps_max=steps_max, t_max=t_max, t=t, CFL=CFL)
         call rk_integrator%integrate(U=domain, stage=rk_stage, dt=dt, t=t)
         t = t + dt
         if (verbose) print "(A)", 'step = '//str(n=step)//', time step = '//str(n=dt)//', time = '//str(n=t)
         if ((t == t_max).or.(step == steps_max)) exit time_loop
         call save_time_serie(t=t)
      enddo time_loop
      call save_time_serie(t=t, finish=.true.)
   enddo
enddo

contains
   subroutine initialize(riemann_solver_scheme, riemann_problem)
   !< Initialize the test.
   character(*), intent(in)                  :: riemann_solver_scheme !< Riemann Problem solver scheme.
   character(*), intent(in)                  :: riemann_problem       !< Riemann problem.
   type(primitive_compressible), allocatable :: initial_state(:)      !< Initial state of primitive variables.
   integer(I4P)                              :: i                     !< Space counter.

   call rk_integrator%initialize(scheme=rk_scheme)
   if (allocated(rk_stage)) deallocate(rk_stage) ; allocate(rk_stage(1:rk_integrator%stages))
   t = 0._R_P
   if (allocated(x)) deallocate(x) ; allocate(x(1:Ni))
   if (allocated(initial_state)) deallocate(initial_state) ; allocate(initial_state(1:Ni))
   Dx = 1._R_P / Ni
   select case(trim(adjustl(riemann_problem)))
   case('sod')
      call riemann_problem_sod_initial_state(initial_state=initial_state)
   case('lax')
      call riemann_problem_lax_initial_state(initial_state=initial_state)
   case('shu-osher')
      call riemann_problem_shu_osher_initial_state(initial_state=initial_state)
   case('123')
      call riemann_problem_123_initial_state(initial_state=initial_state)
   case('woodward-colella')
      call riemann_problem_woodward_colella_initial_state(initial_state=initial_state)
   case('SS1')
      call riemann_problem_ss1_initial_state(initial_state=initial_state)
   case('SS2')
      call riemann_problem_ss2_initial_state(initial_state=initial_state)
   case('SS3')
      call riemann_problem_ss3_initial_state(initial_state=initial_state)
   case('SS4')
      call riemann_problem_ss4_initial_state(initial_state=initial_state)
   case('SS5')
      call riemann_problem_ss5_initial_state(initial_state=initial_state)
   endselect
   call domain%initialize(Ni=Ni, Dx=Dx,                                         &
                          BC_L=BC_L, BC_R=BC_R,                                 &
                          initial_state=initial_state,                          &
                          eos=eos_compressible(cp=1040.004_R_P, cv=742.86_R_P), &
                          weno_scheme=weno_scheme,                              &
                          weno_order=weno_order,                                &
                          weno_variables=weno_variables,                        &
                          riemann_solver_scheme=riemann_solver_scheme)
   endsubroutine initialize

   subroutine parse_command_line_interface()
   !< Parse Command Line Interface (CLI).
   type(command_line_interface)  :: cli                   !< Command line interface handler.
   character(99)                 :: riemann_solver_scheme !< Riemann Problem solver scheme.
   character(99)                 :: riemann_problem       !< Riemann problem.
   integer(I4P)                  :: error                 !< Error handler.
   character(len=:), allocatable :: buffer                !< String buffer.

   call cli%init(description = 'FORESEER test: shock tube tester, 1D Euler equations', &
                 examples    = ["foreseer_test_shock_tube         ",                   &
                                "foreseer_test_shock_tube --tserie"])

   call cli%add(switch='--p', help='Riemann problem', required=.false., act='store', def='sod', &
                choices='all,sod,lax,shu-osher,123,woodward-colella,SS1,SS2,SS3,SS4,SS5')
   call cli%add(switch='--riemann', help='Riemann Problem solver', required=.false., act='store', def='all', &
                choices='all,exact,hllc,llf,pvl,roe')
   call cli%add(switch='--Ni', help='Number finite volumes used', required=.false., act='store', def='100')
   call cli%add(switch='--steps', help='Number time steps performed', required=.false., act='store', def='60')
   call cli%add(switch='--t-max', help='Maximum integration time', required=.false., act='store', def='0.')
   call cli%add(switch='--weno-scheme', help='WENO scheme', required=.false., act='store', def='reconstructor-JS', &
                choices='reconstructor-JS,reconstructor-M-JS,reconstructor-M-Z,reconstructor-Z')
   call cli%add(switch='--weno-order', help='WENO order', required=.false., act='store', def='1')
   call cli%add(switch='--weno-eps', help='WENO epsilon parameter', required=.false., act='store', def='0.000001')
   call cli%add(switch='--weno-variables', help='WENO variables', required=.false., act='store', def='characteristic', &
                choices='characteristic,conservative,primitive')
   call cli%add(switch='--rk-scheme', help='Runge-Kutta intergation scheme', required=.false., act='store',                     &
                def='runge_kutta_ssp_stages_1_order_1',                                                                         &
                choices='runge_kutta_ssp_stages_1_order_1,runge_kutta_ssp_stages_2_order_2,runge_kutta_ssp_stages_3_order_3,'// &
                        'runge_kutta_ssp_stages_5_order_4')
   call cli%add(switch='--cfl', help='CFL value', required=.false., act='store', def='0.7')
   call cli%add(switch='--tserie', switch_ab='-t', help='Save time-serie-result', required=.false., act='store_true', def='.false.')
   call cli%add(switch='--verbose', help='Verbose output', required=.false., act='store_true', def='.false.')

   call cli%parse(error=error)

   call cli%get(switch='--p',              val=riemann_problem,       error=error) ; if (error/=0) stop
   call cli%get(switch='--riemann',        val=riemann_solver_scheme, error=error) ; if (error/=0) stop
   call cli%get(switch='--Ni',             val=Ni,                    error=error) ; if (error/=0) stop
   call cli%get(switch='--steps',          val=steps_max,             error=error) ; if (error/=0) stop
   call cli%get(switch='--t-max',          val=t_max,                 error=error) ; if (error/=0) stop
   call cli%get(switch='--weno-scheme',    val=weno_scheme,           error=error) ; if (error/=0) stop
   call cli%get(switch='--weno-order',     val=weno_order,            error=error) ; if (error/=0) stop
   call cli%get(switch='--weno-eps',       val=weno_eps,              error=error) ; if (error/=0) stop
   call cli%get(switch='--weno-variables', val=weno_variables,        error=error) ; if (error/=0) stop
   call cli%get(switch='--rk-scheme',      val=rk_scheme,             error=error) ; if (error/=0) stop
   call cli%get(switch='--cfl',            val=CFL,                   error=error) ; if (error/=0) stop
   call cli%get(switch='--tserie',         val=time_serie,            error=error) ; if (error/=0) stop
   call cli%get(switch='--verbose',        val=verbose,               error=error) ; if (error/=0) stop

   if (t_max > 0._R_P) steps_max = 0

   if (trim(adjustl(riemann_problem))=='all') then
      riemann_problems = ['sod             ', &
                          'lax             ', &
                          'shu-osher       ', &
                          '123             ', &
                          'woodward-colella', &
                          'SS1             ', &
                          'SS2             ', &
                          'SS3             ', &
                          'SS4             ', &
                          'SS5             ']
   else
      riemann_problems = [trim(adjustl(riemann_problem))]
   endif

   if (trim(adjustl(riemann_solver_scheme))=='all') then
      riemann_solver_schemes = ['exact', &
                                'hllc ', &
                                'llf  ', &
                                'pvl  ', &
                                'roe  ']
   else
      riemann_solver_schemes = [trim(adjustl(riemann_solver_scheme))]
   endif
   endsubroutine parse_command_line_interface

   subroutine riemann_problem_sod_initial_state(initial_state)
   type(primitive_compressible), intent(inout) :: initial_state(1:) !< Initial state of primitive variables.
   integer(I4P)                                :: i                 !< Space counter.

   BC_L = 'TRA'
   BC_R = 'TRA'
   do i=1, Ni / 2
      x(i) = Dx * i - 0.5_R_P * Dx
      initial_state(i)%density  = 1._R_P
      initial_state(i)%velocity = 0._R_P
      initial_state(i)%pressure = 1._R_P
   enddo
   do i=Ni / 2 + 1, Ni
      x(i) = Dx * i - 0.5_R_P * Dx
      initial_state(i)%density  = 0.125_R_P
      initial_state(i)%velocity = 0._R_P
      initial_state(i)%pressure = 0.1_R_P
   enddo
   endsubroutine riemann_problem_sod_initial_state

   subroutine riemann_problem_lax_initial_state(initial_state)
   type(primitive_compressible), intent(inout) :: initial_state(1:) !< Initial state of primitive variables.
   integer(I4P)                                :: i                 !< Space counter.

   BC_L = 'TRA'
   BC_R = 'TRA'
   do i=1, Ni / 2
      x(i) = Dx * i - 0.5_R_P * Dx
      initial_state(i)%density  = 0.445_R_P
      initial_state(i)%velocity = 0.698_R_P
      initial_state(i)%pressure = 3.528_R_P
   enddo
   do i=Ni / 2 + 1, Ni
      x(i) = Dx * i - 0.5_R_P * Dx
      initial_state(i)%density  = 0.5_R_P
      initial_state(i)%velocity = 0._R_P
      initial_state(i)%pressure = 0.571_R_P
   enddo
   endsubroutine riemann_problem_lax_initial_state

   subroutine riemann_problem_shu_osher_initial_state(initial_state)
   type(primitive_compressible), intent(inout) :: initial_state(1:) !< Initial state of primitive variables.
   integer(I4P)                                :: i                 !< Space counter.

   BC_L = 'TRA'
   BC_R = 'TRA'
   do i=1, Ni
      x(i) = Dx * i - 0.5_R_P * Dx
      if (x(i)<1._R_P/8._R_P) then
        initial_state(i)%density  = 3.857143_R_P
        initial_state(i)%velocity = 2.629369_R_P
        initial_state(i)%pressure = 10.3333_R_P
      else
      initial_state(i)%density  = 1._R_P + 0.2_R_P * sin(pi * x(i) / 2._R_P)
      initial_state(i)%velocity = 0._R_P
      initial_state(i)%pressure = 1._R_P
      endif
   enddo
   endsubroutine riemann_problem_shu_osher_initial_state

   subroutine riemann_problem_123_initial_state(initial_state)
   type(primitive_compressible), intent(inout) :: initial_state(1:) !< Initial state of primitive variables.
   integer(I4P)                                :: i                 !< Space counter.

   BC_L = 'TRA'
   BC_R = 'TRA'
   do i=1, Ni / 2
      x(i) = Dx * i - 0.5_R_P * Dx
      initial_state(i)%density  = 1._R_P
      initial_state(i)%velocity = -2._R_P
      initial_state(i)%pressure = 0.4_R_P
   enddo
   do i=Ni / 2 + 1, Ni
      x(i) = Dx * i - 0.5_R_P * Dx
      initial_state(i)%density  = 1._R_P
      initial_state(i)%velocity = 2._R_P
      initial_state(i)%pressure = 0.4_R_P
   enddo
   endsubroutine riemann_problem_123_initial_state

   subroutine riemann_problem_woodward_colella_initial_state(initial_state)
   type(primitive_compressible), intent(inout) :: initial_state(1:) !< Initial state of primitive variables.
   integer(I4P)                                :: i                 !< Space counter.

   BC_L = 'REF'
   BC_R = 'REF'
   do i=1, Ni
      x(i) = Dx * i - 0.5_R_P * Dx
      initial_state(i)%density  = 1._R_P
      initial_state(i)%velocity = 0._R_P
      if (x(i)<=0.1_R_P) then
        initial_state(i)%pressure = 1000._R_P
      elseif (x(i)>0.9_R_P) then
        initial_state(i)%pressure = 100._R_P
      else
        initial_state(i)%pressure = 0.01_R_P
      endif
   enddo
   endsubroutine riemann_problem_woodward_colella_initial_state

   subroutine riemann_problem_ss1_initial_state(initial_state)
   type(primitive_compressible), intent(inout) :: initial_state(1:) !< Initial state of primitive variables.
   integer(I4P)                                :: i                 !< Space counter.

   BC_L = 'TRA'
   BC_R = 'TRA'
   do i=1, Ni / 2
      x(i) = Dx * i - 0.5_R_P * Dx
      initial_state(i)%density  = 1._R_P
      initial_state(i)%velocity = 0._R_P
      initial_state(i)%pressure = 1._R_P
   enddo
   do i=Ni / 2 + 1, Ni
      x(i) = Dx * i - 0.5_R_P * Dx
      initial_state(i)%density  = 1._R_P
      initial_state(i)%velocity = 0._R_P
      initial_state(i)%pressure = 0.1_R_P
   enddo
   endsubroutine riemann_problem_ss1_initial_state

   subroutine riemann_problem_ss2_initial_state(initial_state)
   type(primitive_compressible), intent(inout) :: initial_state(1:) !< Initial state of primitive variables.
   integer(I4P)                                :: i                 !< Space counter.

   BC_L = 'TRA'
   BC_R = 'TRA'
   do i=1, Ni / 2
      x(i) = Dx * i - 0.5_R_P * Dx
      initial_state(i)%density  = 1._R_P
      initial_state(i)%velocity = 0._R_P
      initial_state(i)%pressure = 10._R_P
   enddo
   do i=Ni / 2 + 1, Ni
      x(i) = Dx * i - 0.5_R_P * Dx
      initial_state(i)%density  = 1._R_P
      initial_state(i)%velocity = 0._R_P
      initial_state(i)%pressure = 0.1_R_P
   enddo
   endsubroutine riemann_problem_ss2_initial_state

   subroutine riemann_problem_ss3_initial_state(initial_state)
   type(primitive_compressible), intent(inout) :: initial_state(1:) !< Initial state of primitive variables.
   integer(I4P)                                :: i                 !< Space counter.

   BC_L = 'TRA'
   BC_R = 'TRA'
   do i=1, Ni / 2
      x(i) = Dx * i - 0.5_R_P * Dx
      initial_state(i)%density  = 1._R_P
      initial_state(i)%velocity = 0._R_P
      initial_state(i)%pressure = 100._R_P
   enddo
   do i=Ni / 2 + 1, Ni
      x(i) = Dx * i - 0.5_R_P * Dx
      initial_state(i)%density  = 1._R_P
      initial_state(i)%velocity = 0._R_P
      initial_state(i)%pressure = 0.1_R_P
   enddo
   endsubroutine riemann_problem_ss3_initial_state

   subroutine riemann_problem_ss4_initial_state(initial_state)
   type(primitive_compressible), intent(inout) :: initial_state(1:) !< Initial state of primitive variables.
   integer(I4P)                                :: i                 !< Space counter.

   BC_L = 'TRA'
   BC_R = 'TRA'
   do i=1, Ni / 2
      x(i) = Dx * i - 0.5_R_P * Dx
      initial_state(i)%density  = 1._R_P
      initial_state(i)%velocity = 0._R_P
      initial_state(i)%pressure = 1000._R_P
   enddo
   do i=Ni / 2 + 1, Ni
      x(i) = Dx * i - 0.5_R_P * Dx
      initial_state(i)%density  = 1._R_P
      initial_state(i)%velocity = 0._R_P
      initial_state(i)%pressure = 0.1_R_P
   enddo
   endsubroutine riemann_problem_ss4_initial_state

   subroutine riemann_problem_ss5_initial_state(initial_state)
   type(primitive_compressible), intent(inout) :: initial_state(1:) !< Initial state of primitive variables.
   integer(I4P)                                :: i                 !< Space counter.

   BC_L = 'TRA'
   BC_R = 'TRA'
   do i=1, Ni / 2
      x(i) = Dx * i - 0.5_R_P * Dx
      initial_state(i)%density  = 1._R_P
      initial_state(i)%velocity = 0._R_P
      initial_state(i)%pressure = 10000._R_P
   enddo
   do i=Ni / 2 + 1, Ni
      x(i) = Dx * i - 0.5_R_P * Dx
      initial_state(i)%density  = 1._R_P
      initial_state(i)%velocity = 0._R_P
      initial_state(i)%pressure = 0.1_R_P
   enddo
   endsubroutine riemann_problem_ss5_initial_state

   subroutine save_time_serie(t, filename, finish)
   !< Save time-serie results.
   real(R_P),    intent(in)           :: t         !< Current integration time.
   character(*), intent(in), optional :: filename  !< Output filename.
   logical,      intent(in), optional :: finish    !< Flag for triggering the file closing.
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
         write(tsfile, '(4'//'('//FR_P//',1X))')x(i), primitive%density, primitive%velocity%x, primitive%pressure
      enddo
      if (present(finish)) then
         if (finish) close(tsfile)
      endif
   endif
   endsubroutine save_time_serie
endprogram foreseer_test_shock_tube
