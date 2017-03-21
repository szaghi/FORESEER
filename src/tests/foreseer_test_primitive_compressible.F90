!< FORESEER test: primitive compressible class test.

program foreseer_test_primitive_compressible
!< FORESEER test: primitive compressible class test.

use foreseer, only : eos_compressible, primitive_compressible, primitive_compressible_pointer
use penf, only : R8P, ZeroR8
use vecfor, only : ex, vector

implicit none
type(eos_compressible)                :: eos                  !< An equation of state.
type(primitive_compressible)          :: p                    !< A primitive compressible instance.
type(primitive_compressible)          :: another_p            !< A primitive compressible instance.
type(primitive_compressible), pointer :: p_pointer            !< A primitive compressible pointer.
type(vector)                          :: momentum             !< Momentum vector.
real(R8P), allocatable                :: p_serialized(:)      !< Primitive variable serialized.
real(R8P), allocatable                :: l_eigenvectors(:,:)  !< Left eigenvectors matrix.
real(R8P), allocatable                :: r_eigenvectors(:,:)  !< Right eigenvectors matrix.
real(R8P), allocatable                :: identity(:,:)        !< Identity tensor.
#ifdef __GFORTRAN__
logical                               :: are_tests_passed(16) !< List of passed tests.
#else
logical                               :: are_tests_passed(8)  !< List of passed tests.
#endif

are_tests_passed = .false.

call p%initialize

are_tests_passed(1) = (p%density  == 0._R8P).and. &
                      (p%velocity == 0._R8P).and. &
                      (p%pressure == 0._R8P)
print "(A,L1)", 'primitive = 0, is right? ', are_tests_passed(1)

eos = eos_compressible(cp=1040.004_R8P, cv=742.86_R8P)
p = primitive_compressible(density=1._R8P, pressure=1._R8P)

are_tests_passed(2) = (p%energy(eos=eos) >= 2.5_R8P - ZeroR8).and.(p%energy(eos=eos) <= 2.5_R8P + ZeroR8)
print "(A,L1)", 'p%energy() = 2.5, is right? ', are_tests_passed(2)

momentum = p%momentum()
are_tests_passed(3) = (p%momentum() >= 0._R8P - ZeroR8).and.(p%momentum() <= 0._R8P + ZeroR8)
print "(A,L1)", 'p%momentum() = 0, is right? ', are_tests_passed(3)

p_serialized = p%array()
are_tests_passed(4) = (size(p_serialized, dim=1) == 5).and.            &
                       (p_serialized(1)          == 1._R8P).and.       &
                       (p_serialized(2)          == 0._R8P).and.       &
                       (p_serialized(3)          == 0._R8P).and.       &
                       (p_serialized(4)          == 0._R8P).and.       &
                       (p_serialized(5)          == 1._R8P)
print "(A,L1)", 'p => serialized, is done right? ', are_tests_passed(4)

call p%destroy
are_tests_passed(5) = (p%density  == 0._R8P).and. &
                      (p%velocity == 0._R8P).and. &
                      (p%pressure == 0._R8P)
print "(A,L1)", 'p destroyed, is right? ', are_tests_passed(5)

p = primitive_compressible(density=1._R8P, velocity=ex, pressure=1._R8P)

another_p = p
are_tests_passed(6) = (another_p%density  >= 1._R8P - ZeroR8).and.(another_p%density  <= 1._R8P + ZeroR8).and. &
                      (another_p%velocity >= 1._R8P - ZeroR8).and.(another_p%velocity <= 1._R8P + ZeroR8).and. &
                      (another_p%pressure >= 1._R8P - ZeroR8).and.(another_p%pressure <= 1._R8P + ZeroR8)
print "(A,L1)", 'another_p = p, is done right? ', are_tests_passed(6)

p_pointer => primitive_compressible_pointer(to=p)
are_tests_passed(7) = (p_pointer%density  >= 1._R8P - ZeroR8).and.(p_pointer%density  <= 1._R8P + ZeroR8).and. &
                      (p_pointer%velocity >= 1._R8P - ZeroR8).and.(p_pointer%velocity <= 1._R8P + ZeroR8).and. &
                      (p_pointer%pressure >= 1._R8P - ZeroR8).and.(p_pointer%pressure <= 1._R8P + ZeroR8)
print "(A,L1)", 'p => p, is done right? ', are_tests_passed(7)

call another_p%initialize(initial_state=p)
are_tests_passed(8) = (another_p%density  >= 1._R8P - ZeroR8).and.(another_p%density  <= 1._R8P + ZeroR8).and. &
                      (another_p%velocity >= 1._R8P - ZeroR8).and.(another_p%velocity <= 1._R8P + ZeroR8).and. &
                      (another_p%pressure >= 1._R8P - ZeroR8).and.(another_p%pressure <= 1._R8P + ZeroR8)
print "(A,L1)", 'antoher_p == p, is right? ', are_tests_passed(8)

#ifdef __GFORTRAN__
p = 2._R8P * p
are_tests_passed(9) = (p%density  >= 2._R8P - ZeroR8).and.(p%density  <= 2._R8P + ZeroR8).and. &
                      (p%velocity >= 2._R8P - ZeroR8).and.(p%velocity <= 2._R8P + ZeroR8).and. &
                      (p%pressure >= 2._R8P - ZeroR8).and.(p%pressure <= 2._R8P + ZeroR8)
print "(A,L1)", '2 * p, is done right? ', are_tests_passed(9)

p = p * p
are_tests_passed(10) = (p%density  >= 4._R8P - ZeroR8).and.(p%density  <= 4._R8P + ZeroR8).and. &
                       (p%velocity >= 4._R8P - ZeroR8).and.(p%velocity <= 4._R8P + ZeroR8).and. &
                       (p%pressure >= 4._R8P - ZeroR8).and.(p%pressure <= 4._R8P + ZeroR8)
print "(A,L1)", 'p * p, is done right? ', are_tests_passed(10)

p = p + p
are_tests_passed(11) = (p%density  >= 8._R8P - ZeroR8).and.(p%density  <= 8._R8P + ZeroR8).and. &
                       (p%velocity >= 8._R8P - ZeroR8).and.(p%velocity <= 8._R8P + ZeroR8).and. &
                       (p%pressure >= 8._R8P - ZeroR8).and.(p%pressure <= 8._R8P + ZeroR8)
print "(A,L1)", 'p + p, is done right? ', are_tests_passed(11)

p = p - p
are_tests_passed(12) = (p%density  >= 0._R8P - ZeroR8).and.(p%density  <= 0._R8P + ZeroR8).and. &
                       (p%velocity >= 0._R8P - ZeroR8).and.(p%velocity <= 0._R8P + ZeroR8).and. &
                       (p%pressure >= 0._R8P - ZeroR8).and.(p%pressure <= 0._R8P + ZeroR8)
print "(A,L1)", 'p - p, is done right? ', are_tests_passed(12)

p = primitive_compressible(density=1._R8P, pressure=1._R8P)

another_p = - p
are_tests_passed(13) = (another_p%density  >= -1._R8P - ZeroR8).and.(another_p%density  <= -1._R8P + ZeroR8).and. &
                       (another_p%pressure >= -1._R8P - ZeroR8).and.(another_p%pressure <= -1._R8P + ZeroR8)
print "(A,L1)", 'another_p = - p, is done right? ', are_tests_passed(13)

another_p = + p
are_tests_passed(14) = (another_p%density  >= 1._R8P - ZeroR8).and.(another_p%density  <= 1._R8P + ZeroR8).and. &
                       (another_p%pressure >= 1._R8P - ZeroR8).and.(another_p%pressure <= 1._R8P + ZeroR8)
print "(A,L1)", 'another_p = + p, is done right? ', are_tests_passed(14)

p = p * 2._R8P
are_tests_passed(15) = (p%density  >= 2._R8P - ZeroR8).and.(p%density  <= 2._R8P + ZeroR8).and. &
                       (p%pressure >= 2._R8P - ZeroR8).and.(p%pressure <= 2._R8P + ZeroR8)
print "(A,L1)", 'p * 2, is done right? ', are_tests_passed(15)

p = p / 2._R8P
are_tests_passed(16) = (p%density  >= 1._R8P - ZeroR8).and.(p%density  <= 1._R8P + ZeroR8).and. &
                       (p%pressure >= 1._R8P - ZeroR8).and.(p%pressure <= 1._R8P + ZeroR8)
print "(A,L1)", 'p / 2, is done right? ', are_tests_passed(16)
#endif

p = primitive_compressible(density=1._R8P, velocity=ex, pressure=1._R8P)

print "(A)", "Test pretty printing"
print "(A)", p_pointer%description()

print "(A)", "Test eigenvector computing"
l_eigenvectors = p%left_eigenvectors(eos=eos)
r_eigenvectors = p%right_eigenvectors(eos=eos)
identity = matmul(l_eigenvectors, r_eigenvectors)
print "(A)", "Left"
print "(3F7.3)", l_eigenvectors(1,:)
print "(3F7.3)", l_eigenvectors(2,:)
print "(3F7.3)", l_eigenvectors(3,:)
print "(A)", "Right"
print "(3F7.3)", r_eigenvectors(1,:)
print "(3F7.3)", r_eigenvectors(2,:)
print "(3F7.3)", r_eigenvectors(3,:)
print "(A)", "L * R"
print "(3F7.3)", identity(1,:)
print "(3F7.3)", identity(2,:)
print "(3F7.3)", identity(3,:)

print "(A,L1)", new_line('a')//'Are all tests passed? ', all(are_tests_passed)
endprogram foreseer_test_primitive_compressible
