!< FORESEER test: conservative compressible class test.

program foreseer_test_conservative_compressible
!< FORESEER test: conservative compressible class test.

use foreseer, only : conservative_compressible
use penf, only : R8P, ZeroR8
use vecfor, only : ex, vector

implicit none
type(conservative_compressible)          :: u                    !< A conservative compressible instance.
type(conservative_compressible), pointer :: u_pointer            !< A conservative compressible pointer.
type(conservative_compressible)          :: f                    !< Conservative fluxes.
type(vector)                             :: velocity             !< Velocity vector.
real(R8P), allocatable                   :: u_serialized(:)      !< Conservative variable serialized.
logical                                  :: are_tests_passed(13) !< List of passed tests.

are_tests_passed = .false.

call u%initialize

are_tests_passed(1) = (u%cp       == 0._R8P).and. &
                      (u%cv       == 0._R8P).and. &
                      (u%density  == 0._R8P).and. &
                      (u%momentum == 0._R8P).and. &
                      (u%energy   == 0._R8P)
print "(A,L1)", 'u%cp = 0, is right? ', are_tests_passed(1)

u = conservative_compressible(cp=1040.004_R8P, cv=742.86_R8P, density=1._R8P, energy=2.5_R8P)

are_tests_passed(2) = (u%specific_heats_ratio() >= 1.4_R8P - ZeroR8).and.(u%specific_heats_ratio() <= 1.4_R8P + ZeroR8)
print "(A,L1)", 'u%specific_heats_ratio() = 1.4, is right? ', are_tests_passed(2)

are_tests_passed(3) = (u%pressure() >= 1._R8P - ZeroR8).and.(u%pressure() <= 1._R8P + ZeroR8)
print "(A,L1)", 'u%pressure() = 1, is right? ', are_tests_passed(3)

velocity = u%velocity()
are_tests_passed(4) = (u%velocity() >= 0._R8P - ZeroR8).and.(u%velocity() <= 0._R8P + ZeroR8)
print "(A,L1)", 'u%velocity() = 0, is right? ', are_tests_passed(4)

u_serialized = u%array()
are_tests_passed(5) = (size(u_serialized, dim=1) == 7).and.            &
                       (u_serialized(1)          == 1040.004_R8P).and. &
                       (u_serialized(2)          == 742.86_R8P).and.   &
                       (u_serialized(3)          == 1._R8P).and.       &
                       (u_serialized(4)          == 0._R8P).and.       &
                       (u_serialized(5)          == 0._R8P).and.       &
                       (u_serialized(6)          == 0._R8P).and.       &
                       (u_serialized(7)          == 2.5_R8P)
print "(A,L1)", 'u => serialized, is done right? ', are_tests_passed(5)

call u%destroy
are_tests_passed(6) = (u%cp       == 0._R8P).and. &
                      (u%cv       == 0._R8P).and. &
                      (u%density  == 0._R8P).and. &
                      (u%momentum == 0._R8P).and. &
                      (u%energy   == 0._R8P)
print "(A,L1)", 'u destroyed, is right? ', are_tests_passed(6)

u = conservative_compressible(cp=1040.004_R8P, cv=742.86_R8P, density=1._R8P, momentum=ex, energy=2.5_R8P)

call u%compute_fluxes(normal=ex, fluxes=f)
are_tests_passed(7) = (f%cp       >= 1040.004_R8P - ZeroR8).and.(f%cp       <= 1040.004_R8P + ZeroR8).and. &
                      (f%cv       >= 742.86_R8P -   ZeroR8).and.(f%cv       <= 742.86_R8P   + ZeroR8).and. &
                      (f%density  >= 1._R8P -       ZeroR8).and.(f%density  <= 1._R8P       + ZeroR8).and. &
                      (f%momentum >= 1.8_R8P -      ZeroR8).and.(f%momentum <= 1.8_R8P      + ZeroR8).and. &
                      (f%energy   >= 3.3_R8P -      ZeroR8).and.(f%energy   <= 3.3_R8P      + ZeroR8)
print "(A,L1)", 'compute fluxes along X, is done right? ', are_tests_passed(7)

u = f
are_tests_passed(8) = (u%cp       >= 1040.004_R8P - ZeroR8).and.(u%cp       <= 1040.004_R8P + ZeroR8).and. &
                      (u%cv       >= 742.86_R8P -   ZeroR8).and.(u%cv       <= 742.86_R8P   + ZeroR8).and. &
                      (u%density  >= 1._R8P -       ZeroR8).and.(u%density  <= 1._R8P       + ZeroR8).and. &
                      (u%momentum >= 1.8_R8P -      ZeroR8).and.(u%momentum <= 1.8_R8P      + ZeroR8).and. &
                      (u%energy   >= 3.3_R8P -      ZeroR8).and.(u%energy   <= 3.3_R8P      + ZeroR8)
print "(A,L1)", 'u = f, is done right? ', are_tests_passed(8)

u = conservative_compressible(cp=1040.004_R8P, cv=742.86_R8P, density=1._R8P, momentum=ex, energy=2.5_R8P)

u = 2._R8P * u
are_tests_passed(9) = (u%cp       >= 1040.004_R8P - ZeroR8).and.(u%cp       <= 1040.004_R8P + ZeroR8).and. &
                      (u%cv       >= 742.86_R8P -   ZeroR8).and.(u%cv       <= 742.86_R8P   + ZeroR8).and. &
                      (u%density  >= 2._R8P -       ZeroR8).and.(u%density  <= 2._R8P       + ZeroR8).and. &
                      (u%momentum >= 2._R8P -       ZeroR8).and.(u%momentum <= 2._R8P       + ZeroR8).and. &
                      (u%energy   >= 5._R8P -       ZeroR8).and.(u%energy   <= 5._R8P       + ZeroR8)
print "(A,L1)", '2 * u, is done right? ', are_tests_passed(9)

u = u * u
are_tests_passed(10) = (u%cp       >= 1040.004_R8P - ZeroR8).and.(u%cp       <= 1040.004_R8P + ZeroR8).and. &
                       (u%cv       >= 742.86_R8P -   ZeroR8).and.(u%cv       <= 742.86_R8P   + ZeroR8).and. &
                       (u%density  >= 4._R8P -       ZeroR8).and.(u%density  <= 4._R8P       + ZeroR8).and. &
                       (u%momentum >= 4._R8P -       ZeroR8).and.(u%momentum <= 4._R8P       + ZeroR8).and. &
                       (u%energy   >= 25._R8P -      ZeroR8).and.(u%energy   <= 25._R8P      + ZeroR8)
print "(A,L1)", 'u * u, is done right? ', are_tests_passed(10)

u = u + u
are_tests_passed(11) = (u%cp       >= 1040.004_R8P - ZeroR8).and.(u%cp       <= 1040.004_R8P + ZeroR8).and. &
                       (u%cv       >= 742.86_R8P -   ZeroR8).and.(u%cv       <= 742.86_R8P   + ZeroR8).and. &
                       (u%density  >= 8._R8P -       ZeroR8).and.(u%density  <= 8._R8P       + ZeroR8).and. &
                       (u%momentum >= 8._R8P -       ZeroR8).and.(u%momentum <= 8._R8P       + ZeroR8).and. &
                       (u%energy   >= 50._R8P -      ZeroR8).and.(u%energy   <= 50._R8P      + ZeroR8)
print "(A,L1)", 'u + u, is done right? ', are_tests_passed(11)

u = u - u
are_tests_passed(12) = (u%cp       >= 1040.004_R8P - ZeroR8).and.(u%cp       <= 1040.004_R8P + ZeroR8).and. &
                       (u%cv       >= 742.86_R8P -   ZeroR8).and.(u%cv       <= 742.86_R8P   + ZeroR8).and. &
                       (u%density  >= 0._R8P -       ZeroR8).and.(u%density  <= 0._R8P       + ZeroR8).and. &
                       (u%momentum >= 0._R8P -       ZeroR8).and.(u%momentum <= 0._R8P       + ZeroR8).and. &
                       (u%energy   >= 0._R8P -       ZeroR8).and.(u%energy   <= 0._R8P       + ZeroR8)
print "(A,L1)", 'u - u, is done right? ', are_tests_passed(12)

u = conservative_compressible(cp=1040.004_R8P, cv=742.86_R8P, density=1._R8P, momentum=ex, energy=2.5_R8P)

u_pointer => u_pointer%associate_guarded(to=u)
are_tests_passed(13) = (u%cp       >= 1040.004_R8P - ZeroR8).and.(u%cp       <= 1040.004_R8P + ZeroR8).and. &
                       (u%cv       >= 742.86_R8P -   ZeroR8).and.(u%cv       <= 742.86_R8P   + ZeroR8).and. &
                       (u%density  >= 1._R8P -       ZeroR8).and.(u%density  <= 1._R8P       + ZeroR8).and. &
                       (u%momentum >= 1._R8P -       ZeroR8).and.(u%momentum <= 1._R8P       + ZeroR8).and. &
                       (u%energy   >= 2.5_R8P -      ZeroR8).and.(u%energy   <= 2.5_R8P      + ZeroR8)
print "(A,L1)", 'u => u, is done right? ', are_tests_passed(13)

print "(A,L1)", new_line('a')//'Are all tests passed? ', all(are_tests_passed)
endprogram foreseer_test_conservative_compressible
