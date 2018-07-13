using BinDeps
using Compat.Libdl

@BinDeps.setup

include("verreq.jl")

glpkname = "glpk-$glpkdefver"
glpkwinname = "glpk-$glpkwinver"
glpkdllname = "glpk_$(replace(glpkdefver, ".", "_"))"
glpkwindllname = "glpk_$(replace(glpkwinver, ".", "_"))"

function glpkvalidate(name, handle)
    ver_str = unsafe_string(ccall(Libdl.dlsym(handle, :glp_version), Ptr{UInt8}, ()))
    ver = VersionNumber(ver_str)
    glpkminver <= ver <= glpkmaxver
end
glpkdep = library_dependency("libglpk", aliases = [glpkdllname,glpkwindllname],
                             validate = glpkvalidate)

# Build from sources (used by Linux, BSD)
if VERSION >= v"0.7-"
    bindir = Base.Sys.BINDIR
else
    bindir = JULIA_HOME
end
julia_usrdir = normpath("$bindir/../") # This is a stopgap, we need a better builtin solution to get the included libraries
libdirs = String["$(julia_usrdir)/lib"]
includedirs = String["$(julia_usrdir)/include"]

provides(Sources, Dict(URI("http://ftp.gnu.org/gnu/glpk/$glpkname.tar.gz") => glpkdep), os = :Unix)
provides(BuildProcess, Dict(
    Autotools(libtarget = joinpath("src", ".libs", "libglpk.la"),
              configure_options = AbstractString["--with-gmp"],
              lib_dirs = libdirs,
              include_dirs = includedirs) => glpkdep
    ), os = :Unix)


# Homebrew (OS X section)
if is_apple()
    using Homebrew
    provides(Homebrew.HB, "staticfloat/juliadeps/glpk461", glpkdep, os = :Darwin)
end

# Windows
provides(Binaries, URI("https://bintray.com/artifact/download/tkelman/generic/win$glpkwinname.zip"),
         glpkdep, unpacked_dir="$glpkwinname/w$(Sys.WORD_SIZE)", os = :Windows)

@BinDeps.install Dict(:libglpk => :libglpk)
