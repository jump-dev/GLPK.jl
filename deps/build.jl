using BinDeps

@BinDeps.setup

glpkvers = "4.52"
glpkname = "glpk-$glpkvers"
glpkdllname = "glpk_$(replace(glpkvers, ".", "_"))"
glpkdllnamew(w) = "$(glpkdllname)_$(w)bit"

glpkvalidate(name, handle) = (bytestring(ccall(dlsym(handle, :glp_version), Ptr{Uint8}, ())) == glpkvers)
glpkdep = library_dependency("libglpk", aliases = [glpkdllnamew(WORD_SIZE)], validate = glpkvalidate)

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
glpkdlldirw(w) = joinpath(glpksrcdir, glpkname, "w$w")
glpkdestw(w) = joinpath(glpklibdir, "$(glpkdllnamew(w)).dll")
provides(BuildProcess,
    (@build_steps begin
        FileDownloader("http://downloads.sourceforge.net/project/winglpk/winglpk/GLPK-$glpkvers/win$glpkname.zip",
                       joinpath(glpkdownloadsdir, "win$glpkname.zip"))
        CreateDirectory(glpksrcdir, true)
        FileUnpacker(joinpath(glpkdownloadsdir, "win$glpkname.zip"),
                     glpksrcdir, glpkdlldirw(32))
        CreateDirectory(glpklibdir, true)
        for w in (32,64)
            @build_steps begin
                ChangeDirectory(glpkdlldirw(w))
                FileRule(glpkdestw(w), @build_steps begin
                    `cp $(glpkdllname).dll $(glpkdestw(w))`
                end)
            end
        end
    end), glpkdep, os = :Windows)

@windows_only push!(BinDeps.defaults, BuildProcess)

@BinDeps.install

@windows_only pop!(BinDeps.defaults)
