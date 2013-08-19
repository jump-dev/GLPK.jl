using BinDeps

@BinDeps.setup

glpkvers = "4.48"
glpkname = "glpk-$glpkvers"

glpkdep = library_dependency("libglpk", validate = ((name,handle)->(bytestring(ccall(dlsym(handle, :glp_version), Ptr{Uint8}, ())) == glpkvers)))

provides(Sources, {URI("http://ftp.gnu.org/gnu/glpk/$glpkname.tar.gz") => glpkdep}, os = :Unix)
provides(Sources, {URI("http://downloads.sourceforge.net/project/winglpk/winglpk/GLPK-$glpkvers/win$glpkname.zip") => glpkdep}, os = :Windows)

provides(Homebrew, {"libglpk" => glpkdep})

provides(BuildProcess, {
    Autotools(libtarget = joinpath("src", ".libs", "libglpk.la"), configure_options = String["--with-gmp", "--enable-dl"]) => glpkdep
    }, os = :Unix)

@BinDeps.install
