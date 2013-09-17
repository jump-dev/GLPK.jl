using BinDeps

@BinDeps.setup

glpkvers = "4.52"
glpkname = "glpk-$glpkvers"
glpkdllname = "glpk_$(replace(glpkvers, ".", "_"))"

glpkvalidate(name, handle) = (bytestring(ccall(dlsym(handle, :glp_version), Ptr{Uint8}, ())) == glpkvers)
glpkdep = library_dependency("libglpk", aliases = [glpkdllname], validate = glpkvalidate)

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
    if( Pkg.installed("Homebrew") === nothing )
        error("Homebrew package not installed, please run Pkg.add(\"Homebrew\")")
    else
        using Homebrew
        provides( Homebrew.HB, "glpk", glpkdep, os = :Darwin )
    end
end

# Windows
glpklibdir = BinDeps.libdir(glpkdep)
glpksrcdir = BinDeps.srcdir(glpkdep)
glpkdownloadsdir = BinDeps.downloadsdir(glpkdep)
println(glpksrcdir)
glpkdlldir = joinpath(glpksrcdir, glpkname, "w$WORD_SIZE")
provides(SimpleBuild,
    (@build_steps begin
        FileDownloader("http://downloads.sourceforge.net/project/winglpk/winglpk/GLPK-$glpkvers/win$glpkname.zip",
                       joinpath(glpkdownloadsdir, "win$glpkname.zip"))
        CreateDirectory(glpksrcdir, true)
        FileUnpacker(joinpath(glpkdownloadsdir, "win$glpkname.zip"),
                     glpksrcdir, glpkdlldir)
        CreateDirectory(glpklibdir, true)
        @build_steps begin
            ChangeDirectory(glpkdlldir)
            FileRule(joinpath(glpklibdir, "$glpkdllname.dll"), @build_steps begin
                `cp $glpkdllname.dll $glpklibdir`
            end)
        end
    end), glpkdep, os = :Windows)

@BinDeps.install
