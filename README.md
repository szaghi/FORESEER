<a name="top"></a>

# FORESEER [![GitHub tag](https://img.shields.io/github/tag/szaghi/FORESEER.svg)]() [![Join the chat at https://gitter.im/szaghi/FORESEER](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/szaghi/FORESEER?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

[![License](https://img.shields.io/badge/license-GNU%20GeneraL%20Public%20License%20v3,%20GPLv3-blue.svg)]()
[![License](https://img.shields.io/badge/license-BSD2-red.svg)]()
[![License](https://img.shields.io/badge/license-BSD3-red.svg)]()
[![License](https://img.shields.io/badge/license-MIT-red.svg)]()

[![Status](https://img.shields.io/badge/status-unstable-red.svg)]()
[![Build Status](https://travis-ci.org/szaghi/FORESEER.svg?branch=master)](https://travis-ci.org/szaghi/FORESEER)
[![Coverage Status](https://img.shields.io/codecov/c/github/szaghi/FORESEER.svg)](http://codecov.io/github/szaghi/FORESEER?branch=master)

### FORESEER, FOrtran RiEmann SolvErs EnviRonment

A KISS pure Fortran Library providing a large set of Riemann Solvers:

- FORESEER is a pure Fortran (KISS) library for solving Riemann problems (1D);
- FORESEER is Fortran 2008+ standard compliant;
- FORESEER is OOP designed;
- FORESEER is TDD developed;
- FORESEER is a Free, Open Source Project.

#### A taste of FORESEER

```fortran
use foreseer
type(eos_compressible)                 :: eos            ! Equation of state.
type(conservative_compressible)        :: state_left     ! Left state.
type(conservative_compressible)        :: state_right    ! Right state.
type(conservative_compressible)        :: fluxes         ! Conservative fluxes.
class(riemann_solver_compressible_llf) :: riemann_solver ! Riemman Problem solver.

! air
eos = eos_compressible(cp=1040.004_R8P, cv=742.86_R8P)

! SOD's Riemann Problem
state_left  = conservative_compressible(density=1._R8P,                          &
                                        energy=1._R8P*eos%energy(density=1._R8P, &
                                                                 pressure=1._R8P))
state_right = conservative_compressible(density=0.125_R8P,                             &
                                        energy=0.125_R8P*eos%energy(density=0.125_R8P, &
                                                                    pressure=0.1_R8P))

! solve Riemann Problem
call riemann_solver%solve(eos_left=eos, state_left=state_left,               &
                          eos_right=eos, state_right=state_right, normal=ex, &
                          fluxes=fluxes)

! print results
print '(A)', 'Fluxes at interface:'
print '(A)', fluxes%description(prefix='  ')
```

#### Issues

[![GitHub issues](https://img.shields.io/github/issues/szaghi/FORESEER.svg)]()
[![Ready in backlog](https://badge.waffle.io/szaghi/FORESEER.png?label=ready&title=Ready)](https://waffle.io/szaghi/FORESEER)
[![In Progress](https://badge.waffle.io/szaghi/FORESEER.png?label=in%20progress&title=In%20Progress)](https://waffle.io/szaghi/FORESEER)
[![Open bugs](https://badge.waffle.io/szaghi/FORESEER.png?label=bug&title=Open%20Bugs)](https://waffle.io/szaghi/FORESEER)

#### Compiler Support

[![Compiler](https://img.shields.io/badge/GNU-v6.3.1+-brightgreen.svg)]()
[![Compiler](https://img.shields.io/badge/Intel-v17.0.1+-brightgreen.svg)]()
[![Compiler](https://img.shields.io/badge/IBM%20XL-not%20tested-yellow.svg)]()
[![Compiler](https://img.shields.io/badge/g95-not%20tested-yellow.svg)]()
[![Compiler](https://img.shields.io/badge/NAG-not%20tested-yellow.svg)]()
[![Compiler](https://img.shields.io/badge/PGI-not%20tested-yellow.svg)]()

---

[What is FORESEER?](#what-is-foreseer) | [Main features](#main-features) | [Copyrights](#copyrights) | [Download](#download) | [Compilation](#compilation) | [Documentation](#documentation) | [References](#references)

---

## What is FORESEER?

> **FORESEER** is a modern Fortran library providing a large set of Riemann Solvers.

To be completed.

### How to use

To be written.

Go to [Top](#top)

## Main features

To be written.

Any feature request is welcome.

Go to [Top](#top)

## Copyrights

FORESEER is a Free and Open Source Software (FOSS), it is distributed under a **very permissive** multi-licensing system: selectable licenses are [GPLv3](http://www.gnu.org/licenses/gpl-3.0.html), [BSD2-Clause](http://opensource.org/licenses/BSD-2-Clause), [BSD3-Clause](http://opensource.org/licenses/BSD-3-Clause) and [MIT](http://opensource.org/licenses/MIT), feel free to select the license that best matches your workflow.

> Anyone is interest to use, to develop or to contribute to FORESEER is welcome.

More details can be found on [wiki](https://github.com/szaghi/FORESEER/wiki/Copyrights).

Go to [Top](#top)

## Download

To be written.

Go to [Top](#top)

## Compilation

To be written.

## Documentation

Besides this README file the FORESEER documentation is contained into its own [wiki](https://github.com/szaghi/FORESEER/wiki). Detailed documentation of the API is contained into the [GitHub Pages](http://szaghi.github.io/FORESEER/index.html) that can also be created locally by means of [ford tool](https://github.com/cmacmackin/ford).

Go to [Top](#top)
