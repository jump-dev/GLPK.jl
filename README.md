Julia GLPK module
=================


| **Documentation**                                                               | **Build Status**                                                                                    |
|:-------------------------------------------------------------------------------:|:---------------------------------------------------------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-latest-img]][docs-latest-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] [![][coveralls-img]][coveralls-url] |


GLPK.jl is a wrapper for the [GNU Linear Programming Kit library](http://www.gnu.org/software/glpk).
It makes it possible to access nearly all of GLPK functionality from within Julia programs.

See also the [GLPKMathProgInterface.jl](https://github.com/JuliaOpt/GLPKMathProgInterface.jl) package for using it with
[MathProgBase.jl](https://github.com/JuliaOpt/MathProgBase.jl) and [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl).

This package is part of [the JuliaOpt project](http://www.juliaopt.org/).

## Installation

The package is registered in the [General registry](https://github.com/JuliaRegistries/General/) and so can be installed with `Pkg.add`.

```julia
julia> import Pkg

julia> Pkg.add("GLPK")
```

GLPK.jl will use [BinaryProvider.jl](https://github.com/JuliaPackaging/BinaryProvider.jl) to automatically install the GLPK binaries with [GMP](https://gmplib.org) support.

## Custom Installation

To install custom built GLPK binaries set the environmental variable `JULIA_GLPK_LIBRARY_PATH` and call `import Pkg; Pkg.build("GLPK")`. For instance, if the libraries are installed in `/opt/lib` just call
```julia
ENV["JULIA_GLPK_LIBRARY_PATH"]="/opt/lib"
Pkg.build("GLPK")
```

If you do not want BinaryProvider to download the default binaries on install set  `JULIA_GLPK_LIBRARY_PATH`  before calling `import Pkg; Pkg.add("GLPK")`.

To switch back to the default binaries clear `JULIA_GLPK_LIBRARY_PATH` and call `import Pkg; Pkg.build("GLPK")`.

## Documentation

- [**STABLE**][docs-stable-url] &mdash; **most recently tagged version of the documentation.**
- [**LATEST**][docs-latest-url] &mdash; *in-development version of the documentation.*

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://gplkjl.readthedocs.org/en/latest/glpk.html

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://gplkjl.readthedocs.org/en/stable/glpk.html

[travis-img]: https://api.travis-ci.org/JuliaOpt/GLPK.jl.svg?branch=master
[travis-url]: https://travis-ci.org/JuliaOpt/GLPK.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/4t5e2dir3gp7fb6h?svg=true
[appveyor-url]: https://ci.appveyor.com/project/JuliaOpt/glpk-jl

[coveralls-img]: https://img.shields.io/coveralls/JuliaOpt/GLPK.jl.svg
[coveralls-url]: https://coveralls.io/r/JuliaOpt/GLPK.jl

[pkg-0.6-img]: http://pkg.julialang.org/badges/GLPK_0.6.svg
[pkg-0.6-url]: http://pkg.julialang.org/?pkg=GLPK
[pkg-0.7-img]: http://pkg.julialang.org/badges/GLPK_0.7.svg
[pkg-0.7-url]: http://pkg.julialang.org/?pkg=GLPK
