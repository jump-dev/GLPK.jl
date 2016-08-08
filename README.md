Julia GLPK module
=================


| **Documentation**                                                               | **PackageEvaluator**                                            | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:---------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-latest-img]][docs-latest-url] | [![][pkg-0.4-img]][pkg-0.4-url] [![][pkg-0.5-img]][pkg-0.5-url] | [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] [![][coveralls-img]][coveralls-url] |


GLPK.jl is a wrapper for the [GNU Linear Programming Kit library](http://www.gnu.org/software/glpk).
It makes it possible to access nearly all of GLPK functionality from within Julia programs.

See also the [GLPKMathProgInterface.jl](https://github.com/JuliaOpt/GLPKMathProgInterface.jl) package for using it with
[MathProgBase.jl](https://github.com/JuliaOpt/MathProgBase.jl) and [JuMP.jl](https://github.com/JuliaOpt/JuMP.jl).

This package is part of [the JuliaOpt project](http://www.juliaopt.org/).

## Installation

The package is registered in `METADATA.jl` and so can be installed with `Pkg.add`.

```
julia> Pkg.add("GLPK")
```

In case `Pkg.add("GLPK")` gives you an error on Linux, you may need to install the GMP library headers.
For example, on Ubuntu/Debian and similar, give the following command from a terminal:

```
$ sudo apt-get install libgmp-dev
```

After that, restart the installation of the package with:

```
julia> Pkg.build("GLPK")
```


## Documentation

- [**STABLE**][docs-stable-url] &mdash; **most recently tagged version of the documentation.**
- [**LATEST**][docs-latest-url] &mdash; *in-development version of the documentation.*

## Project Status

The package is tested against Julia `0.4`, `0.5` and *current* `0.6-dev` on Linux, OS X, and Windows.

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: https://gplkjl.readthedocs.org/en/latest/glpk.html

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://gplkjl.readthedocs.org/en/stable/glpk.html

[travis-img]: https://api.travis-ci.org/JuliaOpt/GLPK.jl.svg?branch=master
[travis-url]: https://travis-ci.org/JuliaOpt/GLPK.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/9yaxhsw24br4c0ux/branch/master?svg=true
[appveyor-url]: https://ci.appveyor.com/project/tkelman/glpk-jl/branch/master

[coveralls-img]: https://img.shields.io/coveralls/JuliaOpt/GLPK.jl.svg
[coveralls-url]: https://coveralls.io/r/JuliaOpt/GLPK.jl

[pkg-0.4-img]: http://pkg.julialang.org/badges/GLPK_0.4.svg
[pkg-0.4-url]: http://pkg.julialang.org/?pkg=GLPK
[pkg-0.5-img]: http://pkg.julialang.org/badges/GLPK_0.5.svg
[pkg-0.5-url]: http://pkg.julialang.org/?pkg=GLPK
