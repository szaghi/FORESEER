!< FORESEER test: equation of state compressible class test.

program foreseer_test_eos_compressible
!< FORESEER test: equation of state compressible class test.

use foreseer, only : eos_compressible
use penf, only : R8P, ZeroR8

implicit none
type(eos_compressible) :: eos                  !< An equation of state.
logical                :: are_tests_passed(14) !< List of passed tests.

are_tests_passed = .false.

eos = eos_compressible(cp=1040.004_R8P, cv=742.86_R8P)
are_tests_passed(1) = (eos%g() >= 1.4_R8P - ZeroR8).and.(eos%g() <= 1.4_R8P + ZeroR8)
print "(A,L1)", 'eos%g() = 1.4, is right? ', are_tests_passed(1)

eos = eos_compressible(cp=1040.004_R8P, gam=1.4_R8P)
are_tests_passed(2) = (eos%cv() >= 742.86_R8P - ZeroR8).and.(eos%cv() <= 742.86_R8P + ZeroR8)
print "(A,L1)", 'eos%cv() = 742.86, is right? ', are_tests_passed(2)

eos = eos_compressible(cv=742.86_R8P, gam=1.4_R8P)
are_tests_passed(3) = (eos%R() >= 297.144_R8P - 1000*ZeroR8).and.(eos%R() <= 297.144_R8P + 1000*ZeroR8)
print "(A,L1)", 'eos%R() = 297.144, is right? ', are_tests_passed(3)

eos = eos_compressible(cv=742.86_R8P, R=297.144_R8P)
are_tests_passed(4) = (eos%cp() >= 1040.004_R8P - ZeroR8).and.(eos%cp() <= 1040.004_R8P + ZeroR8)
print "(A,L1)", 'eos%cp() = 1040.004, is right? ', are_tests_passed(4)

eos = eos_compressible(cp=1040.004_R8P, R=297.144_R8P)
are_tests_passed(5) = (eos%cv() >= 742.86_R8P - 1000*ZeroR8).and.(eos%cv() <= 742.86_R8P + 1000*ZeroR8)
print "(A,L1)", 'eos%cv() = 742.86, is right? ', are_tests_passed(5)

eos = eos_compressible(gam=1.4_R8P, R=297.144_R8P)
are_tests_passed(6) = (eos%cp() >= 1040.004_R8P - 1000*ZeroR8).and.(eos%cp() <= 1040.004_R8P + 1000*ZeroR8)
print "(A,L1)", 'eos%cp() = 1040.004, is right? ', are_tests_passed(6)

eos = eos_compressible(cp=1040.004_R8P, cv=742.86_R8P)
are_tests_passed(7) = (eos%pressure(density=1._R8P, energy=1._R8P) >= 0.4_R8P - ZeroR8).and.&
                      (eos%pressure(density=1._R8P, energy=1._R8P) <= 0.4_R8P + ZeroR8)
print "(A,L1)", 'eos%pressure(density=1, energy=1) = 0.4, is right? ', are_tests_passed(7)

are_tests_passed(8) = (eos%pressure(density=1._R8P, temperature=1._R8P) >= 297.144_R8P - 1000*ZeroR8).and.&
                      (eos%pressure(density=1._R8P, temperature=1._R8P) <= 297.144_R8P + 1000*ZeroR8)
print "(A,L1)", 'eos%pressure(density=1, temperature=1) = 297.144, is right? ', are_tests_passed(8)

are_tests_passed(9) = (eos%density(pressure=1._R8P, temperature=1._R8P) >= 1._R8P/297.144_R8P - 1000*ZeroR8).and.&
                      (eos%density(pressure=1._R8P, temperature=1._R8P) <= 1._R8P/297.144_R8P + 1000*ZeroR8)
print "(A,L1)", 'eos%density(pressure=1, temperature=1) = 1/297.144, is right? ', are_tests_passed(9)

are_tests_passed(10) = (eos%density(pressure=1._R8P, energy=1._R8P) >= 1._R8P/0.4_R8P - ZeroR8).and.&
                       (eos%density(pressure=1._R8P, energy=1._R8P) <= 1._R8P/0.4_R8P + ZeroR8)
print "(A,L1)", 'eos%density(pressure=1, energy=1) = 1/0.4, is right? ', are_tests_passed(10)

are_tests_passed(11) = (eos%energy(density=1._R8P, pressure=1._R8P) >= 1._R8P/0.4_R8P - ZeroR8).and.&
                       (eos%energy(density=1._R8P, pressure=1._R8P) <= 1._R8P/0.4_R8P + ZeroR8)
print "(A,L1)", 'eos%energy(density=1, pressure=1) = 1/0.4, is right? ', are_tests_passed(11)

are_tests_passed(12) = (eos%energy(temperature=1._R8P) >= 742.86_R8P - ZeroR8).and.&
                       (eos%energy(temperature=1._R8P) <= 742.86_R8P + ZeroR8)
print "(A,L1)", 'eos%energy(temperature=1) = 742.86, is right? ', are_tests_passed(12)

are_tests_passed(13) = (eos%temperature(density=1._R8P, pressure=1._R8P) >= 1._R8P/297.144_R8P - ZeroR8).and.&
                       (eos%temperature(density=1._R8P, pressure=1._R8P) <= 1._R8P/297.144_R8P + ZeroR8)
print "(A,L1)", 'eos%temperature(density=1, pressure=1) = 1/297.144, is right? ', are_tests_passed(13)

are_tests_passed(14) = (eos%temperature(energy=1._R8P) >= 1._R8P/742.86_R8P - ZeroR8).and.&
                       (eos%temperature(energy=1._R8P) <= 1._R8P/742.86_R8P + ZeroR8)
print "(A,L1)", 'eos%temperature(energy=1) = 1/742.86, is right? ', are_tests_passed(14)

print "(A,L1)", new_line('a')//'Are all tests passed? ', all(are_tests_passed)
endprogram foreseer_test_eos_compressible
