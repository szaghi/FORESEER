!< Define the abstract primitive compressible state of a Riemann Problem for FORESEER library.

module foreseer_primitive_compressible
!< Define the abstract primitive compressible state of a Riemann Problem for FORESEER library.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use foreseer_primitive_object, only : primitive_object
use foreseer_eos_object, only : eos_object
use penf, only : R8P, str
use vecfor, only : vector

implicit none
private
public :: primitive_compressible
public :: primitive_compressible_pointer

type, extends(primitive_object) :: primitive_compressible
   !< Convervative compressible object class.
   real(R8P)    :: density=0._R8P  !< Density, `rho`.
   type(vector) :: velocity        !< Velocity, `v`.
   real(R8P)    :: pressure=0._R8P !< Pressure, `p`.
   contains
      ! public methods
      procedure, pass(self) :: left_eigenvectors  !< Return the left eigenvectors matrix `L` as `dF/dP = A = R ^ L`.
      procedure, pass(self) :: right_eigenvectors !< Return the right eigenvectors matrix `R` as `dF/dP = A = R ^ L`.
      ! deferred methods
      procedure, pass(self) :: array              !< Return serialized array of primitive.
      procedure, pass(self) :: description        !< Return pretty-printed object description.
      procedure, pass(self) :: destroy            !< Destroy primitive.
      procedure, pass(self) :: energy             !< Return energy value.
      procedure, pass(self) :: initialize         !< Initialize primitive.
      procedure, pass(self) :: momentum           !< Return momentum vector.
      procedure, pass(lhs)  :: prim_assign_prim   !< Operator `=`.
      procedure, pass(lhs)  :: prim_divide_real   !< Operator `prim / real`.
      procedure, pass(lhs)  :: prim_multiply_real !< Operator `prim * real`.
      procedure, pass(lhs)  :: prim_multiply_prim !< Operator `*`.
      procedure, pass(rhs)  :: real_multiply_prim !< Operator `real * prim`.
      procedure, pass(lhs)  :: add                !< Operator `+`.
      procedure, pass(self) :: positive           !< Unary operator `+ prim`.
      procedure, pass(lhs)  :: sub                !< Operator `-`.
      procedure, pass(self) :: negative           !< Unary operator `- prim`.
endtype primitive_compressible

interface primitive_compressible
   !< Overload [[primitive_compressible]] name with its constructor.
   module procedure primitive_compressible_instance
endinterface

contains
   ! public non TBP
   function primitive_compressible_pointer(to, error_message) result(pointer_)
   !< Return [[primitive_compressible]] pointer associated to [[primitive_object]] or its extensions until
   !< [[primitive_compressible]] included.
   !<
   !< @note A type-guard check is performed and error stop is raised if necessary.
   class(primitive_object), intent(in), target   :: to            !< Target of associate.
   character(*),            intent(in), optional :: error_message !< Auxiliary error message.
   class(primitive_compressible), pointer        :: pointer_      !< Associated pointer.

   select type(to)
   type is(primitive_compressible)
      pointer_ => to
   class default
      write(stderr, '(A)') 'error: cast primitive_object to primitive_compressible failed!'
      if (present(error_message)) write(stderr, '(A)') error_message
      stop
   endselect
   endfunction primitive_compressible_pointer

   ! public methods
   pure function left_eigenvectors(self, eos) result(eig)
   !< Return the left eigenvectors matrix `L` as `dF/dP = A = R ^ L`.
   class(primitive_compressible), intent(in) :: self          !< Primitive.
   class(eos_object),             intent(in) :: eos           !< Equation of state.
   real(R8P)                                 :: eig(1:3, 1:3) !< Eigenvectors.
   real(R8P)                                 :: gp            !< `g*p`.
   real(R8P)                                 :: gp_a          !< `g*p/a`.

   gp = eos%gam() * self%pressure
   gp_a = gp / eos%speed_of_sound(density=self%density, pressure=self%pressure)
   eig(1, 1) = 0._R8P            ; eig(1, 2) = -gp_a  ; eig(1, 3) =  1._R8P
   eig(2, 1) = gp / self%density ; eig(2, 2) = 0._R8P ; eig(2, 3) = -1._R8P
   eig(3, 1) = 0._R8P            ; eig(3, 2) =  gp_a  ; eig(3, 3) =  1._R8P
   endfunction left_eigenvectors

   pure function right_eigenvectors(self, eos) result(eig)
   !< Return the right eigenvectors matrix `R` as `dF/dP = A = R ^ L`.
   class(primitive_compressible), intent(in) :: self          !< Primitive.
   class(eos_object),             intent(in) :: eos           !< Equation of state.
   real(R8P)                                 :: eig(1:3, 1:3) !< Eigenvectors.
   real(R8P)                                 :: gp            !< `g*p`.
   real(R8P)                                 :: gp_inv        !< `1/(g*p)`.
   real(R8P)                                 :: a             !< Speed of sound, `sqrt(g*p/r)`.

   gp = eos%gam() * self%pressure
   gp_inv = 1._R8P / gp
   a = eos%speed_of_sound(density=self%density, pressure=self%pressure)
   eig(1, 1) =  0.5_R8P * self%density * gp_inv ; eig(1, 2) = self%density * gp_inv  ; eig(1, 3) =  eig(1, 1)
   eig(2, 1) = -0.5_R8P * a * gp_inv            ; eig(2, 2) = 0._R8P                 ; eig(2, 3) = -eig(2, 1)
   eig(3, 1) =  0.5_R8P                         ; eig(3, 2) = 0._R8P                 ; eig(3, 3) =  eig(3, 1)
   endfunction right_eigenvectors

   ! deferred methods
   pure function array(self) result(array_)
   !< Return serialized array of primitive.
   class(primitive_compressible), intent(in) :: self      !< Primitive.
   real(R8P), allocatable                    :: array_(:) !< Serialized array of primitive.

   allocate(array_(1:5))
   array_(1) = self%density
   array_(2) = self%velocity%x
   array_(3) = self%velocity%y
   array_(4) = self%velocity%z
   array_(5) = self%pressure
   endfunction array

   pure function description(self, prefix) result(desc)
   !< Return a pretty-formatted object description.
   class(primitive_compressible), intent(in)           :: self             !< Primitive.
   character(*),                  intent(in), optional :: prefix           !< Prefixing string.
   character(len=:), allocatable                       :: prefix_          !< Prefixing string, local variable.
   character(len=:), allocatable                       :: desc             !< Description.
   character(len=1), parameter                         :: NL=new_line('a') !< New line character.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   desc = ''
   desc = desc//prefix_//'density  = '//trim(str(n=self%density))//NL
   desc = desc//prefix_//'velocity = '//trim(str(n=[self%velocity%x, self%velocity%y, self%velocity%z]))//NL
   desc = desc//prefix_//'pressure = '//trim(str(n=self%pressure))
   endfunction description

   elemental subroutine destroy(self)
   !< Destroy primitive.
   class(primitive_compressible), intent(inout) :: self  !< Primitive.
   type(primitive_compressible)                 :: fresh !< Fresh instance of primitive object.

   self = fresh
   endsubroutine destroy

   elemental function energy(self, eos) result(energy_)
   !< Return energy value.
   class(primitive_compressible), intent(in) :: self    !< Primitive.
   class(eos_object),             intent(in) :: eos     !< Equation of state.
   real(R8P)                                 :: energy_ !< Energy value.

   energy_ = self%pressure / (eos%gam() - 1._R8P) + 0.5_R8P * self%density * self%velocity%sq_norm()
   endfunction energy

   subroutine initialize(self, initial_state)
   !< Initialize primitive.
   class(primitive_compressible), intent(inout)        :: self          !< Primitive.
   class(primitive_object),       intent(in), optional :: initial_state !< Initial state.

   if (present(initial_state)) then
      select type(initial_state)
      class is(primitive_compressible)
         self = initial_state
      endselect
   else
      call self%destroy
   endif
   endsubroutine initialize

   elemental function momentum(self) result(momentum_)
   !< Return momentum vector.
   class(primitive_compressible), intent(in) :: self      !< Primitive.
   type(vector)                              :: momentum_ !< Momentum vector.

   momentum_ = self%density * self%velocity
   endfunction momentum

   ! operators
   pure subroutine prim_assign_prim(lhs, rhs)
   !< Operator `=`.
   class(primitive_compressible), intent(inout) :: lhs !< Left hand side.
   class(primitive_object),       intent(in)    :: rhs !< Right hand side.

   select type(rhs)
   class is (primitive_compressible)
      lhs%density  = rhs%density
      lhs%velocity = rhs%velocity
      lhs%pressure = rhs%pressure
   endselect
   endsubroutine prim_assign_prim

   function prim_divide_real(lhs, rhs) result(operator_result)
   !< Operator `prim / real`.
   class(primitive_compressible), intent(in) :: lhs             !< Left hand side.
   real(R8P),                     intent(in) :: rhs             !< Right hand side.
   class(primitive_object), allocatable      :: operator_result !< Operator result.

   allocate(primitive_compressible :: operator_result)
   select type(operator_result)
   class is(primitive_compressible)
      operator_result%density  = lhs%density  / rhs
      operator_result%velocity = lhs%velocity / rhs
      operator_result%pressure = lhs%pressure / rhs
   endselect
   endfunction prim_divide_real

   function prim_multiply_real(lhs, rhs) result(operator_result)
   !< Operator `prim * real`.
   class(primitive_compressible), intent(in) :: lhs             !< Left hand side.
   real(R8P),                     intent(in) :: rhs             !< Right hand side.
   class(primitive_object), allocatable      :: operator_result !< Operator result.

   allocate(primitive_compressible :: operator_result)
   select type(operator_result)
   class is(primitive_compressible)
      operator_result%density  = lhs%density  * rhs
      operator_result%velocity = lhs%velocity * rhs
      operator_result%pressure = lhs%pressure * rhs
   endselect
   endfunction prim_multiply_real

   function real_multiply_prim(lhs, rhs) result(operator_result)
   !< Operator `real * prim`.
   real(R8P),                     intent(in) :: lhs             !< Left hand side.
   class(primitive_compressible), intent(in) :: rhs             !< Right hand side.
   class(primitive_object), allocatable      :: operator_result !< Operator result.

   allocate(primitive_compressible :: operator_result)
   select type(operator_result)
   class is(primitive_compressible)
      operator_result%density  = lhs * rhs%density
      operator_result%velocity = lhs * rhs%velocity
      operator_result%pressure = lhs * rhs%pressure
   endselect
   endfunction real_multiply_prim

   function prim_multiply_prim(lhs, rhs) result(operator_result)
   !< Operator `*`.
   class(primitive_compressible), intent(in) :: lhs             !< Left hand side.
   class(primitive_object),       intent(in) :: rhs             !< Right hand side.
   class(primitive_object), allocatable      :: operator_result !< Operator result.

   allocate(primitive_compressible :: operator_result)
   select type(operator_result)
   class is(primitive_compressible)
      operator_result = lhs
      select type(rhs)
      class is (primitive_compressible)
         operator_result%density  = lhs%density  * rhs%density
         operator_result%velocity = lhs%velocity * rhs%velocity
         operator_result%pressure = lhs%pressure * rhs%pressure
      endselect
   endselect
   endfunction prim_multiply_prim

   function add(lhs, rhs) result(operator_result)
   !< Operator `+`.
   class(primitive_compressible), intent(in) :: lhs             !< Left hand side.
   class(primitive_object),       intent(in) :: rhs             !< Right hand side.
   class(primitive_object), allocatable      :: operator_result !< Operator result.

   allocate(primitive_compressible :: operator_result)
   select type(operator_result)
   class is(primitive_compressible)
      operator_result = lhs
      select type(rhs)
      class is (primitive_compressible)
         operator_result%density  = lhs%density  + rhs%density
         operator_result%velocity = lhs%velocity + rhs%velocity
         operator_result%pressure = lhs%pressure + rhs%pressure
      endselect
   endselect
   endfunction add

   function positive(self) result(operator_result)
   !< Unary operator `+ prim`.
   class(primitive_compressible), intent(in) :: self            !< Primitive.
   class(primitive_object), allocatable      :: operator_result !< Operator result.

   allocate(primitive_compressible :: operator_result)
   select type(operator_result)
   class is(primitive_compressible)
      operator_result%density  = + self%density
      operator_result%velocity = + self%velocity
      operator_result%pressure = + self%pressure
   endselect
   endfunction positive

   function sub(lhs, rhs) result(operator_result)
   !< Operator `+`.
   class(primitive_compressible), intent(in) :: lhs             !< Left hand side.
   class(primitive_object),       intent(in) :: rhs             !< Right hand side.
   class(primitive_object), allocatable      :: operator_result !< Operator result.

   allocate(primitive_compressible :: operator_result)
   select type(operator_result)
   class is(primitive_compressible)
      operator_result = lhs
      select type(rhs)
      class is (primitive_compressible)
         operator_result%density  = lhs%density  - rhs%density
         operator_result%velocity = lhs%velocity - rhs%velocity
         operator_result%pressure = lhs%pressure - rhs%pressure
      endselect
   endselect
   endfunction sub

   function negative(self) result(operator_result)
   !< Unary operator `- prim`.
   class(primitive_compressible), intent(in) :: self            !< Primitive.
   class(primitive_object), allocatable      :: operator_result !< Operator result.

   allocate(primitive_compressible :: operator_result)
   select type(operator_result)
   class is(primitive_compressible)
      operator_result%density  = - self%density
      operator_result%velocity = - self%velocity
      operator_result%pressure = - self%pressure
   endselect
   endfunction negative

   ! private non TBP
   elemental function primitive_compressible_instance(density, velocity, pressure) result(instance)
   !< Return and instance of [[primitive_compressible]].
   !<
   !< @note This procedure is used for overloading [[primitive_compressible]] name.
   real(R8P),    intent(in), optional :: density  !< Density, `rho`.
   type(vector), intent(in), optional :: velocity !< Velocity, `v`.
   real(R8P),    intent(in), optional :: pressure !< Pressure, `p`.
   type(primitive_compressible)       :: instance !< Instance of [[primitive_compressible]].

   if (present(density)) instance%density = density
   if (present(velocity)) instance%velocity = velocity
   if (present(pressure)) instance%pressure = pressure
   endfunction primitive_compressible_instance
endmodule foreseer_primitive_compressible
