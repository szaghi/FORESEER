!< FORESEER test: Riemann solver compressible Roe class test.

program foreseer_test_riemann_solver_compressible_roe
!< FORESEER test: Riemann solver compressible Roe class test.

use foreseer, only : eos_compressible, conservative_compressible, riemann_solver_compressible_roe
use penf, only : R8P, str
use vecfor, only : ex

implicit none
type(eos_compressible)                :: eos                          !< The equation of state.
type(conservative_compressible)       :: state_left                   !< Left state.
type(conservative_compressible)       :: state_right                  !< Right state.
type(conservative_compressible)       :: fluxes                       !< Conservative fluxes.
type(riemann_solver_compressible_roe) :: riemann_solver               !< Riemann solver.
real(R8P)                             :: waves(1:5)                   !< Waves pattern.
real(R8P), parameter                  :: r_2=0.426319003105163574_R8P !< Exact value of density in state 2.
real(R8P), parameter                  :: r_3=0.265574008226394653_R8P !< Exact value of density in state 3.
real(R8P), parameter                  :: p23=0.303130000829696655_R8P !< Exact value of pressure in states 2 and 3.
real(R8P), parameter                  :: u23=0.927452981472015381_R8P !< Exact value of velocity in states 2 and 3.
logical                               :: are_tests_passed(1)          !< List of passed tests.

are_tests_passed = .true.

eos = eos_compressible(cp=1040.004_R8P, cv=742.86_R8P)
state_left  = conservative_compressible(density=1._R8P,    energy=1._R8P   *eos%energy(density=1._R8P,    pressure=1._R8P) )
state_right = conservative_compressible(density=0.125_R8P, energy=0.125_R8P*eos%energy(density=0.125_R8P, pressure=0.1_R8P))

print '(A)', 'Test solution with "u23" algorithm:'
call riemann_solver%initialize(config='u23')
call riemann_solver%solve(eos_left=eos, state_left=state_left, eos_right=eos, state_right=state_right, normal=ex, fluxes=fluxes)
print '(A)', 'Fluxes at interface:'
print '(A)', fluxes%description(prefix='  ')
call fluxes%compute_fluxes_from_primitive(eos=eos, p=p23, r=r_2, u=u23, normal=ex)
print '(A)', 'Exact fluxes at interface:'
print '(A)', fluxes%description(prefix='  ')
print '(A)', 'Exact intemediate states:'
print '(A)', '  r_2 = '//str(n=r_2)
print '(A)', '  r_3 = '//str(n=r_3)
print '(A)', '  p23 = '//str(n=p23)
print '(A)', '  u23 = '//str(n=u23)
print '(A)', 'Test solution with "up23" algorithm:'
call riemann_solver%initialize(config='up23')
call riemann_solver%solve(eos_left=eos, state_left=state_left, eos_right=eos, state_right=state_right, normal=ex, fluxes=fluxes)
print '(A)', 'Fluxes at interface:'
print '(A)', fluxes%description(prefix='  ')
call fluxes%compute_fluxes_from_primitive(eos=eos, p=p23, r=r_2, u=u23, normal=ex)
print '(A)', 'Exact fluxes at interface:'
print '(A)', fluxes%description(prefix='  ')
print '(A)', 'Exact intemediate states:'
print '(A)', '  r_2 = '//str(n=r_2)
print '(A)', '  r_3 = '//str(n=r_3)
print '(A)', '  p23 = '//str(n=p23)
print '(A)', '  u23 = '//str(n=u23)

call riemann_solver%compute_waves(eos_left=eos, state_left=state_left, eos_right=eos, state_right=state_right, normal=ex, &
                                  waves=waves)
print '(A)', 'Waves pattern:'
print '(A)', '  S1 = '//str(n=waves(1))
print '(A)', '  S2 = '//str(n=waves(2))
print '(A)', '  S  = '//str(n=waves(3))
print '(A)', '  S3 = '//str(n=waves(4))
print '(A)', '  S4 = '//str(n=waves(5))

print "(A,L1)", new_line('a')//'Are all tests passed? ', all(are_tests_passed)
endprogram foreseer_test_riemann_solver_compressible_roe
