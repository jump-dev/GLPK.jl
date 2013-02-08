require("BinDeps")

glpkvers = "4.48"
glpkname = "glpk-$glpkvers"
@unix_only glpkarchive = "$glpkname.tar.gz"
@windows_only glpkarghive = "win$(glpkname).zip"
glpkprefix = Pkg.dir("GLPK", "deps", "usr")

tagfile = "installed_vers"

if !isfile(tagfile) || readchomp(tagfile) != glpkvers
    @unix_only begin
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
    end

    @windows_only begin
        error("sorry, Windows is currently unsupported by the GLPK module")
        if !isfile("$glpkarchive")
            run(download_cmd("http://ftp.gnu.org/gnu/glpk/$glpkarchive", glpkarchive))
            run(download_cmd("http://downloads.sourceforge.net/project/winglpk/winglpk/GLPK-$glpkvers/$glpkarchive", glpkarchive))
        end
        run(`7z x $glpkarchive -y`)
        # TODO move binaries to usr directory
    end

    run(`echo $glpkvers` > tagfile)
end
