require("BinDeps")

glpkvers = "4.47"
glpkname = "glpk-$glpkvers"
glpkarchive = "$glpkname.tar.gz"
glpkprefix = joinpath(Pkg.dir(), "GLPK", "deps", "usr")
glpkincdir = joinpath(glpkprefix, "include")

CC = get(ENV, "CC", "cc")

tagfile = "installed_vers"

if !isfile(tagfile) || readchomp(tagfile) != glpkvers
    GLPK_CONST = "0x[0-9a-fA-F]+|[-+]?\\s*[0-9]+"
    PERL_CMD = "/^\\s*#define\\s+GLP_(\\w*)\\s*\\(?($GLPK_CONST)\\)?\\s*$$/ and print \"const $$1 = int32($$2)\""

    glpkheader = joinpath(glpkincdir, "glpk.h")

    if !isfile("$glpkarchive")
        run(download_cmd("http://ftp.gnu.org/gnu/glpk/$glpkarchive", glpkarchive))
    end
    run(unpack_cmd(glpkarchive, "."))
    cd("$glpkname") do
        run(`./configure --prefix=$glpkprefix`)
        run(`make install`)
    end

    glpk_h_jl = joinpath(Pkg.dir(), "GLPK", "src", "glpk_h.jl")
    run(`$CC -E -dM $glpkheader` | `perl -nle $PERL_CMD` | (`sort` > glpk_h_jl))

    run(`echo $glpkvers` > tagfile)
end
