# Brownies

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://edwinb-ai.github.io/Brownies.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://edwinb-ai.github.io/Brownies.jl/dev)
[![CI](https://github.com/edwinb-ai/Brownies.jl/workflows/CI/badge.svg)](https://github.com/edwinb-ai/Brownies.jl/actions?query=workflow%3ACI)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/edwinb-ai/Brownies.jl?svg=true)](https://ci.appveyor.com/project/edwinb-ai/Brownies-jl)
[![Codecov](https://codecov.io/gh/edwinb-ai/Brownies.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/edwinb-ai/Brownies.jl)
[![Coverage Status](https://coveralls.io/repos/github/edwinb-ai/Brownies.jl/badge.svg?branch=master)](https://coveralls.io/github/edwinb-ai/Brownies.jl?branch=master)

## Library
Brownian dynamics library for molecular dynamics simulations using the classical
Ermak-McCammon algorithm.

## Observables
The current implementation can handle the following physical observables:

- _Radial distribution function_ to study the structure of the fluid.
- _Mean square displacement_ to study the dynamics of the molecules in the fluid following a probe particle.
- _Structure factor_ to further study the structure of the fluid when being probed by a particle beam, using scattering properties.

