!< FORESEER test: compressible variables transformations test.

program foreseer_test_compressible_transformations
!< FORESEER test: compressible variables transformations test.

use foreseer, only : eos_compressible, conservative_compressible, primitive_compressible, &
                     conservative_to_primitive_compressible, primitive_to_conservative_compressible
use penf, only : R8P, ZeroR8

implicit none
type(eos_compressible)          :: eos                 !< An equation of state.
type(primitive_compressible)    :: p                   !< A primitive compressible instance.
type(conservative_compressible) :: u                   !< A conservative compressible instance.
logical                         :: are_tests_passed(2) !< List of passed tests.

are_tests_passed = .false.

eos = eos_compressible(cp=1040.004_R8P, cv=742.86_R8P)
u = conservative_compressible(density=1._R8P, energy=2.5_R8P)
p = conservative_to_primitive_compressible(conservative=u, eos=eos)
are_tests_passed(1) = (p%density  >= 1._R8P - ZeroR8).and.(p%density  <= 1._R8P + ZeroR8).and. &
                      (p%velocity >= 0._R8P - ZeroR8).and.(p%velocity <= 0._R8P + ZeroR8).and. &
                      (p%pressure >= 1._R8P - ZeroR8).and.(p%pressure <= 1._R8P + ZeroR8)
print "(A,L1)", 'u to p, is done right? ', are_tests_passed(1)

p = primitive_compressible(density=1._R8P, pressure=1._R8P)
u = primitive_to_conservative_compressible(primitive=p, eos=eos)
are_tests_passed(2) = (u%density  >= 1._R8P  - ZeroR8).and.(u%density  <= 1._R8P  + ZeroR8).and. &
                      (u%momentum >= 0._R8P  - ZeroR8).and.(u%momentum <= 0._R8P  + ZeroR8).and. &
                      (u%energy   >= 2.5_R8P - ZeroR8).and.(u%energy   <= 2.5_R8P + ZeroR8)
print "(A,L1)", 'p to u, is done right? ', are_tests_passed(2)

print "(A,L1)", new_line('a')//'Are all tests passed? ', all(are_tests_passed)
endprogram foreseer_test_compressible_transformations
