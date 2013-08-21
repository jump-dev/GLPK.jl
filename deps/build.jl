using BinDeps

@BinDeps.setup

glpkvers = "4.48"
glpkname = "glpk-$glpkvers"

glpkdep = library_dependency("libglpk", validate = ((name,handle)->(bytestring(ccall(dlsym(handle, :glp_version), Ptr{Uint8}, ())) == glpkvers)))

provides(Sources, {URI("http://ftp.gnu.org/gnu/glpk/$glpkname.tar.gz") => glpkdep}, os = :Unix)
provides(Sources, {URI("http://downloads.sourceforge.net/project/winglpk/winglpk/GLPK-$glpkvers/win$glpkname.zip") => glpkdep}, os = :Windows)

provides(Homebrew, {"https://raw.github.com/Homebrew/homebrew-science/master/glpk.rb" => glpkdep})

julia_usrdir = normpath(JULIA_HOME*"/../") # This is a stopgap, we need a better builtin solution to get the included libraries
libdirs = String["$(julia_usrdir)/lib"]
includedirs = String["$(julia_usrdir)/include"]

provides(BuildProcess, {
    Autotools(libtarget = joinpath("src", ".libs", "libglpk.la"), 
              # configure_options = String["--with-gmp", "--enable-dl"],
              configure_options = String["--with-gmp"],
              lib_dirs = libdirs, 
              include_dirs = includedirs) => glpkdep
    }, os = :Unix)

@BinDeps.install
