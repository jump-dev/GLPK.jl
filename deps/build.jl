using BinDeps

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
glpkdep = library_dependency("libglpk", aliases = ["libglpk-40"], # it is called libglpk-40 on the glpk-devel WinRPM package
                             validate = glpkvalidate)

# Build from sources (used by Linux, BSD)
julia_usrdir = normpath("$JULIA_HOME/../") # This is a stopgap, we need a better builtin solution to get the included libraries
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
@static if is_windows()
    # GMP
    using WinRPM
    WinRPM.install("libgmp10", yes = true)
    push!(WinRPM.sources, "https://cache.julialang.org/http://download.opensuse.org/repositories/home:/blegat:/branches:/windows:/mingw:/win32/openSUSE_13.2")
    push!(WinRPM.sources, "https://cache.julialang.org/http://download.opensuse.org/repositories/home:/blegat:/branches:/windows:/mingw:/win64/openSUSE_13.2")
    WinRPM.update()
    provides(WinRPM.RPM, "glpk-devel", [glpkdep], os = :Windows)
end

@BinDeps.install Dict(:libglpk => :libglpk)
