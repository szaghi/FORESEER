!< FORESEER, FOrtran RiEmann SolvErs EnviRonment.
module foreseer
!< FORESEER, FOrtran RiEmann SolvErs EnviRonment.

use, intrinsic :: iso_fortran_env, only : stderr=>error_unit
use flow_compressible_transformations, only : conservative_to_primitive_compressible, primitive_to_conservative_compressible
use flow_conservative_compressible, only : conservative_compressible, conservative_compressible_pointer
use flow_conservative_object, only : conservative_object
use flow_eos_compressible, only : eos_compressible
use flow_eos_object, only : eos_object
use flow_primitive_compressible, only : primitive_compressible, primitive_compressible_pointer
use flow_primitive_object, only : primitive_object
use foreseer_riemann_pattern_compressible_object, only : riemann_pattern_compressible_object
use foreseer_riemann_pattern_compressible_pvl, only : riemann_pattern_compressible_pvl
use foreseer_riemann_pattern_object, only : riemann_pattern_object
use foreseer_riemann_solver_compressible_exact, only : riemann_solver_compressible_exact, riemann_solver_compressible_exact_id
use foreseer_riemann_solver_compressible_hllc, only : riemann_solver_compressible_hllc, riemann_solver_compressible_hllc_id
use foreseer_riemann_solver_compressible_llf, only : riemann_solver_compressible_llf, riemann_solver_compressible_llf_id
use foreseer_riemann_solver_compressible_pvl, only : riemann_solver_compressible_pvl, riemann_solver_compressible_pvl_id
use foreseer_riemann_solver_compressible_roe, only : riemann_solver_compressible_roe, riemann_solver_compressible_roe_id
use foreseer_riemann_solver_object, only : riemann_solver_object

implicit none
private
public :: conservative_to_primitive_compressible, primitive_to_conservative_compressible
public :: conservative_compressible, conservative_compressible_pointer
public :: conservative_object
public :: eos_compressible
public :: eos_object
public :: primitive_compressible, primitive_compressible_pointer
public :: primitive_object
public :: riemann_pattern_compressible_object
public :: riemann_pattern_compressible_pvl
public :: riemann_pattern_object
public :: riemann_solver_compressible_exact, riemann_solver_compressible_exact_id
public :: riemann_solver_compressible_hllc, riemann_solver_compressible_hllc_id
public :: riemann_solver_compressible_llf, riemann_solver_compressible_llf_id
public :: riemann_solver_compressible_pvl, riemann_solver_compressible_pvl_id
public :: riemann_solver_compressible_roe, riemann_solver_compressible_roe_id
public :: riemann_solver_object
public :: foreseer_factory

contains
   subroutine foreseer_factory(riemann_solver_scheme, riemann_solver, config)
   !< FORESEER factory, create concrete instance of Riemann Solver.
   character(len=*),                          intent(in)           :: riemann_solver_scheme !< Riemann Solver scheme ID.
   class(riemann_solver_object), allocatable, intent(inout)        :: riemann_solver        !< Riemann Solver.
   character(len=*),                          intent(in), optional :: config                !< Configuration for solver algorithm.

   if (allocated(riemann_solver)) then
      call riemann_solver%destroy
      deallocate(riemann_solver)
   endif
   select case(trim(adjustl(riemann_solver_scheme)))
   case(riemann_solver_compressible_exact_id)
      allocate(riemann_solver_compressible_exact :: riemann_solver)
   case(riemann_solver_compressible_hllc_id)
      allocate(riemann_solver_compressible_hllc :: riemann_solver)
   case(riemann_solver_compressible_llf_id)
      allocate(riemann_solver_compressible_llf :: riemann_solver)
   case(riemann_solver_compressible_pvl_id)
      allocate(riemann_solver_compressible_pvl :: riemann_solver)
   case(riemann_solver_compressible_roe_id)
      allocate(riemann_solver_compressible_roe :: riemann_solver)
   case default
      write(stderr, '(A)') 'error: Riemann Solver scheme "'//riemann_solver_scheme//'" is unknown!'
      write(stderr, '(A)') 'valid schemes are:'
      write(stderr, '(A)') riemann_solver_compressible_exact_id
      write(stderr, '(A)') riemann_solver_compressible_hllc_id
      write(stderr, '(A)') riemann_solver_compressible_llf_id
      write(stderr, '(A)') riemann_solver_compressible_pvl_id
      write(stderr, '(A)') riemann_solver_compressible_roe_id
      stop
   endselect
   call riemann_solver%initialize(config=config)
   endsubroutine foreseer_factory
endmodule foreseer
