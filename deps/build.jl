require("BinDeps")

glpkvers = "4.47"
glpkname = "glpk-$glpkvers"
glpkarchive = "$glpkname.tar.gz"
glpkprefix = joinpath(Pkg.dir(), "GLPK", "deps", "usr")

tagfile = "installed_vers"

if !isfile(tagfile) || readchomp(tagfile) != glpkvers
    if !isfile("$glpkarchive")
        run(download_cmd("http://ftp.gnu.org/gnu/glpk/$glpkarchive", glpkarchive))
    end
    run(unpack_cmd(glpkarchive, "."))
    cd("$glpkname") do
        run(`./configure --prefix=$glpkprefix --with-gmp --enable-dl`)
        run(`make`)
        run(`make check`)
        run(`make install`)
    end

    run(`echo $glpkvers` > tagfile)
end
