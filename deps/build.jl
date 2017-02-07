using BinDeps

@BinDeps.setup

include("verreq.jl")

glpkname = "glpk-$glpkdefver"
glpkdllname = "glpk_$(replace(glpkdefver, ".", "_"))"

function glpkvalidate(name, handle)
    ver_str = unsafe_string(ccall(Libdl.dlsym(handle, :glp_version), Ptr{UInt8}, ()))
    ver = VersionNumber(ver_str)
    glpkminver <= ver <= glpkmaxver
end
glpkdep = library_dependency("libglpk", aliases = [glpkdllname],
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
    if Homebrew.installed("glpk") # remove old conflicting version
        Homebrew.rm("glpk")
    end
    provides(Homebrew.HB, "glpk452", glpkdep, os = :Darwin)
end

# Windows
provides(Binaries, URI("https://bintray.com/artifact/download/tkelman/generic/win$glpkname.zip"),
         glpkdep, unpacked_dir="$glpkname/w$(Sys.WORD_SIZE)", os = :Windows)

@BinDeps.install Dict(:libglpk => :libglpk)
