!< FORESEER test: Riemann solver LLF class test.

program foreseer_test_riemann_solver_llf
!< FORESEER test: Riemann solver LLF class test.

use foreseer, only : eos_compressible, conservative_compressible, riemann_pattern_compressible, riemann_solver_llf
use penf, only : R8P, str
use vecfor, only : ex

implicit none
type(eos_compressible)             :: eos                 !< The equation of state.
type(conservative_compressible)    :: state_left          !< Left state.
type(conservative_compressible)    :: state_right         !< Right state.
type(conservative_compressible)    :: fluxes              !< Conservative fluxes.
type(riemann_solver_llf)           :: riemann_solver      !< Riemann solver.
type(riemann_pattern_compressible) :: riemann_pattern     !< Riemann pattern.
logical                            :: are_tests_passed(1) !< List of passed tests.

are_tests_passed = .true.

eos = eos_compressible(cp=1040.004_R8P, cv=742.86_R8P)
state_left  = conservative_compressible(density=1._R8P,    energy=1._R8P   *eos%energy(density=1._R8P,    pressure=1._R8P) )
state_right = conservative_compressible(density=0.125_R8P, energy=0.125_R8P*eos%energy(density=0.125_R8P, pressure=0.1_R8P))

call riemann_solver%solve(eos_left=eos, state_left=state_left, eos_right=eos, state_right=state_right, normal=ex, &
                          fluxes=fluxes, pattern=riemann_pattern)
print '(A)', 'Fluxes at interface:'
print '(A)', fluxes%description(prefix='  ')
print '(A)', 'Riemann (waves) pattern approximation:'
print '(A)', riemann_pattern%description(prefix='  ')
print '(A)', 'Exact intemediate states:'
print '(A)', '  r_2 = '//str(n=0.426319003105163574_R8P)
print '(A)', '  r_3 = '//str(n=0.265574008226394653_R8P)
print '(A)', '  p23 = '//str(n=0.303130000829696655_R8P)
print '(A)', '  u23 = '//str(n=0.927452981472015381_R8P)

print "(A,L1)", new_line('a')//'Are all tests passed? ', all(are_tests_passed)
endprogram foreseer_test_riemann_solver_llf
