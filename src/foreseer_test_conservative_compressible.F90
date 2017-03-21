!< FORESEER test: conservative compressible class test.

program foreseer_test_conservative_compressible
!< FORESEER test: conservative compressible class test.

use foreseer, only : eos_compressible, conservative_compressible, conservative_compressible_pointer
use penf, only : R8P, ZeroR8
use vecfor, only : ex, vector

implicit none
type(eos_compressible)                   :: eos                  !< An equation of state.
type(conservative_compressible)          :: u                    !< A conservative compressible instance.
type(conservative_compressible)          :: another_u            !< Another conservative compressible instance.
type(conservative_compressible), pointer :: u_pointer            !< A conservative compressible pointer.
type(conservative_compressible)          :: f                    !< Conservative fluxes.
type(vector)                             :: velocity             !< Velocity vector.
real(R8P), allocatable                   :: u_serialized(:)      !< Conservative variable serialized.
#ifdef __GFORTRAN__
logical                                  :: are_tests_passed(17) !< List of passed tests.
#else
logical                                  :: are_tests_passed(9)  !< List of passed tests.
#endif

are_tests_passed = .false.

call u%initialize

are_tests_passed(1) = (u%density  == 0._R8P).and. &
                      (u%momentum == 0._R8P).and. &
                      (u%energy   == 0._R8P)
print "(A,L1)", 'conservative = 0, is right? ', are_tests_passed(1)

eos = eos_compressible(cp=1040.004_R8P, cv=742.86_R8P)
u = conservative_compressible(density=1._R8P, energy=2.5_R8P)

are_tests_passed(2) = (u%pressure(eos=eos) >= 1._R8P - ZeroR8).and.(u%pressure(eos=eos) <= 1._R8P + ZeroR8)
print "(A,L1)", 'u%pressure() = 1, is right? ', are_tests_passed(2)

velocity = u%velocity()
are_tests_passed(3) = (u%velocity() >= 0._R8P - ZeroR8).and.(u%velocity() <= 0._R8P + ZeroR8)
print "(A,L1)", 'u%velocity() = 0, is right? ', are_tests_passed(3)

u_serialized = u%array()
are_tests_passed(4) = (size(u_serialized, dim=1) == 5).and.            &
                       (u_serialized(1)          == 1._R8P).and.       &
                       (u_serialized(2)          == 0._R8P).and.       &
                       (u_serialized(3)          == 0._R8P).and.       &
                       (u_serialized(4)          == 0._R8P).and.       &
                       (u_serialized(5)          == 2.5_R8P)
print "(A,L1)", 'u => serialized, is done right? ', are_tests_passed(4)

call u%destroy
are_tests_passed(5) = (u%density  == 0._R8P).and. &
                      (u%momentum == 0._R8P).and. &
                      (u%energy   == 0._R8P)
print "(A,L1)", 'u destroyed, is right? ', are_tests_passed(5)

u = conservative_compressible(density=1._R8P, momentum=ex, energy=2.5_R8P)

call u%compute_fluxes(eos=eos, normal=ex, fluxes=f)
are_tests_passed(6) = (f%density  >= 1._R8P -  ZeroR8).and.(f%density  <= 1._R8P  + ZeroR8).and. &
                      (f%momentum >= 1.8_R8P - ZeroR8).and.(f%momentum <= 1.8_R8P + ZeroR8).and. &
                      (f%energy   >= 3.3_R8P - ZeroR8).and.(f%energy   <= 3.3_R8P + ZeroR8)
print "(A,L1)", 'compute fluxes along X, is done right? ', are_tests_passed(6)

u = f
are_tests_passed(7) = (u%density  >= 1._R8P -  ZeroR8).and.(u%density  <= 1._R8P  + ZeroR8).and. &
                      (u%momentum >= 1.8_R8P - ZeroR8).and.(u%momentum <= 1.8_R8P + ZeroR8).and. &
                      (u%energy   >= 3.3_R8P - ZeroR8).and.(u%energy   <= 3.3_R8P + ZeroR8)
print "(A,L1)", 'u = f, is done right? ', are_tests_passed(7)

u = conservative_compressible(density=1._R8P, momentum=ex, energy=2.5_R8P)

u_pointer => conservative_compressible_pointer(to=u)
are_tests_passed(8) = (u_pointer%density  >= 1._R8P -  ZeroR8).and.(u_pointer%density  <= 1._R8P  + ZeroR8).and. &
                      (u_pointer%momentum >= 1._R8P -  ZeroR8).and.(u_pointer%momentum <= 1._R8P  + ZeroR8).and. &
                      (u_pointer%energy   >= 2.5_R8P - ZeroR8).and.(u_pointer%energy   <= 2.5_R8P + ZeroR8)
print "(A,L1)", 'u => u, is done right? ', are_tests_passed(8)

call another_u%initialize(initial_state=u)
are_tests_passed(9) = (another_u%density  >= 1._R8P -  ZeroR8).and.(another_u%density  <= 1._R8P  + ZeroR8).and. &
                      (another_u%momentum >= 1._R8P -  ZeroR8).and.(another_u%momentum <= 1._R8P  + ZeroR8).and. &
                      (another_u%energy   >= 2.5_R8P - ZeroR8).and.(another_u%energy   <= 2.5_R8P + ZeroR8)
print "(A,L1)", 'antoher_u == u, is right? ', are_tests_passed(9)

#ifdef __GFORTRAN__
u = conservative_compressible(density=1._R8P, momentum=ex, energy=2.5_R8P)
u = 2._R8P * u
are_tests_passed(10) = (u%density  >= 2._R8P - ZeroR8).and.(u%density  <= 2._R8P + ZeroR8).and. &
                       (u%momentum >= 2._R8P - ZeroR8).and.(u%momentum <= 2._R8P + ZeroR8).and. &
                       (u%energy   >= 5._R8P - ZeroR8).and.(u%energy   <= 5._R8P + ZeroR8)
print "(A,L1)", '2 * u, is done right? ', are_tests_passed(10)

u = u * u
are_tests_passed(11) = (u%density  >= 4._R8P -  ZeroR8).and.(u%density  <= 4._R8P  + ZeroR8).and. &
                       (u%momentum >= 4._R8P -  ZeroR8).and.(u%momentum <= 4._R8P  + ZeroR8).and. &
                       (u%energy   >= 25._R8P - ZeroR8).and.(u%energy   <= 25._R8P + ZeroR8)
print "(A,L1)", 'u * u, is done right? ', are_tests_passed(11)

u = u + u
are_tests_passed(12) = (u%density  >= 8._R8P -  ZeroR8).and.(u%density  <= 8._R8P  + ZeroR8).and. &
                       (u%momentum >= 8._R8P -  ZeroR8).and.(u%momentum <= 8._R8P  + ZeroR8).and. &
                       (u%energy   >= 50._R8P - ZeroR8).and.(u%energy   <= 50._R8P + ZeroR8)
print "(A,L1)", 'u + u, is done right? ', are_tests_passed(12)

u = u - u
are_tests_passed(13) = (u%density  >= 0._R8P - ZeroR8).and.(u%density  <= 0._R8P + ZeroR8).and. &
                       (u%momentum >= 0._R8P - ZeroR8).and.(u%momentum <= 0._R8P + ZeroR8).and. &
                       (u%energy   >= 0._R8P - ZeroR8).and.(u%energy   <= 0._R8P + ZeroR8)
print "(A,L1)", 'u - u, is done right? ', are_tests_passed(13)

u = conservative_compressible(density=1._R8P, energy=2.5_R8P)

u = - u
are_tests_passed(14) = (u%density >= -1._R8P  - ZeroR8).and.(u%density <= -1._R8P  + ZeroR8).and. &
                       (u%energy  >= -2.5_R8P - ZeroR8).and.(u%energy  <= -2.5_R8P + ZeroR8)
print "(A,L1)", 'u = - u, is done right? ', are_tests_passed(14)

u = + u
are_tests_passed(15) = (u%density >= -1._R8P  - ZeroR8).and.(u%density <= -1._R8P  + ZeroR8).and. &
                       (u%energy  >= -2.5_R8P - ZeroR8).and.(u%energy  <= -2.5_R8P + ZeroR8)
print "(A,L1)", 'u = + u, is done right? ', are_tests_passed(15)

u = u * 2._R8P
are_tests_passed(16) = (u%density >= -2._R8P - ZeroR8).and.(u%density <= -2._R8P + ZeroR8).and. &
                       (u%energy  >= -5._R8P - ZeroR8).and.(u%energy  <= -5._R8P + ZeroR8)
print "(A,L1)", 'u * 2, is done right? ', are_tests_passed(16)

u = u / 2._R8P
are_tests_passed(17) = (u%density >= -1._R8P  - ZeroR8).and.(u%density <= -1._R8P  + ZeroR8).and. &
                       (u%energy  >= -2.5_R8P - ZeroR8).and.(u%energy  <= -2.5_R8P + ZeroR8)
print "(A,L1)", 'u / 2, is done right? ', are_tests_passed(17)
#endif

print "(A,L1)", new_line('a')//'Are all tests passed? ', all(are_tests_passed)
endprogram foreseer_test_conservative_compressible
