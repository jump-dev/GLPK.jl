using BinDeps

using Compat
using Compat: unsafe_string, is_apple

@BinDeps.setup

include("verreq.jl")

glpkname = "glpk-$glpkdefver"
glpkdllname = "glpk_$(replace(glpkdefver, ".", "_"))"

const _dlsym = (VERSION >= v"0.4.0-dev+3844" ? Libdl.dlsym : dlsym)

if VERSION >= v"0.4-dev"
    function string_version(str)
        VersionNumber(str)
    end
else
    function string_version(str)
        major_ver, minor_ver = match(r"(\d+)\.(\d+)", str).captures
        # No need for @compat since this is only on <= 0.3
        VersionNumber(parseint(major_ver), parseint(minor_ver))
    end
end

function glpkvalidate(name, handle)
    ver_str = unsafe_string(ccall(_dlsym(handle, :glp_version), Ptr{UInt8}, ()))
    ver = string_version(ver_str)
    glpkminver <= ver <= glpkmaxver
end
glpkdep = library_dependency("libglpk", aliases = [glpkdllname],
                             validate = glpkvalidate)

# Build from sources (used by Linux, BSD)
julia_usrdir = normpath("$JULIA_HOME/../") # This is a stopgap, we need a better builtin solution to get the included libraries
libdirs = AbstractString["$(julia_usrdir)/lib"]
includedirs = AbstractString["$(julia_usrdir)/include"]

@compat provides(Sources, Dict(URI("http://ftp.gnu.org/gnu/glpk/$glpkname.tar.gz") => glpkdep), os = :Unix)
@compat provides(BuildProcess, Dict(
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

@compat @BinDeps.install Dict(:libglpk => :libglpk)
