using BinDeps

@BinDeps.setup

glpkvers = "4.52"
glpkname = "glpk-$glpkvers"
glpkdllname = "glpk_$(replace(glpkvers, ".", "_"))"
glpkdllname32 = "$(glpkdllname)_32bit"
glpkdllname64 = "$(glpkdllname)_64bit"

glpkvalidate(name, handle) = (bytestring(ccall(dlsym(handle, :glp_version), Ptr{Uint8}, ())) == glpkvers)
glpkdep = library_dependency("libglpk", aliases = (WORD_SIZE == 32 ? [glpkdllname32] : [glpkdllname64]), validate = glpkvalidate)

# Build from sources (used by Linux, BSD)
julia_usrdir = normpath("$JULIA_HOME/../") # This is a stopgap, we need a better builtin solution to get the included libraries
libdirs = String["$(julia_usrdir)/lib"]
includedirs = String["$(julia_usrdir)/include"]

provides(Sources, {URI("http://ftp.gnu.org/gnu/glpk/$glpkname.tar.gz") => glpkdep}, os = :Unix)
provides(BuildProcess, {
    Autotools(libtarget = joinpath("src", ".libs", "libglpk.la"),
              configure_options = String["--with-gmp"],
              lib_dirs = libdirs,
              include_dirs = includedirs) => glpkdep
    }, os = :Unix)


# Homebrew (OS X section)
@osx_only begin
    using Homebrew
    provides(Homebrew.HB, "glpk", glpkdep, os = :Darwin)
end

# Windows
glpklibdir = BinDeps.libdir(glpkdep)
glpksrcdir = BinDeps.srcdir(glpkdep)
glpkdownloadsdir = BinDeps.downloadsdir(glpkdep)
glpkdlldir32 = joinpath(glpksrcdir, glpkname, "w32")
glpkdlldir64 = joinpath(glpksrcdir, glpkname, "w64")
glpkdest32 = joinpath(glpklibdir, "$glpkdllname32.dll")
glpkdest64 = joinpath(glpklibdir, "$glpkdllname64.dll")
provides(BuildProcess,
    (@build_steps begin
        FileDownloader("http://downloads.sourceforge.net/project/winglpk/winglpk/GLPK-$glpkvers/win$glpkname.zip",
                       joinpath(glpkdownloadsdir, "win$glpkname.zip"))
        CreateDirectory(glpksrcdir, true)
        FileUnpacker(joinpath(glpkdownloadsdir, "win$glpkname.zip"),
                     glpksrcdir, glpkdlldir32)
        CreateDirectory(glpklibdir, true)
        @build_steps begin
            ChangeDirectory(glpkdlldir32)
            FileRule(glpkdest32, @build_steps begin
                `cp $(glpkdllname).dll $(glpkdest32)`
            end)
        end
        @build_steps begin
            ChangeDirectory(glpkdlldir64)
            FileRule(glpkdest64, @build_steps begin
                `cp $(glpkdllname).dll $(glpkdest64)`
            end)
        end
    end), glpkdep, os = :Windows)

@windows_only push!(BinDeps.defaults, BuildProcess)

@BinDeps.install

@windows_only pop!(BinDeps.defaults)
