# This script is used to auto-generate the file GLPK_constants.jl
# It is intended for use by the Package maintaners whenever the GLPK
# version is updated.

require("BinDeps")

glpkvers = "4.47"
glpkname = "glpk-$glpkvers"
glpkarchive = "$glpkname.tar.gz"

CC = get(ENV, "CC", "cc")

GLPK_CONST = "0x[0-9a-fA-F]+|[-+]?\\s*[0-9]+"
PERL_CMD = "/^\\s*#define\\s+GLP_(\\w*)\\s*\\(?($GLPK_CONST)\\)?\\s*$$/ and print \"const $$1 = int32($$2)\""

if !isfile("$glpkarchive")
    run(download_cmd("http://ftp.gnu.org/gnu/glpk/$glpkarchive", glpkarchive))
end
run(unpack_cmd(glpkarchive, "."))

glpkheader = joinpath(glpkname, "src", "glpk.h")

glpk_h_jl = joinpath("..", "src", "GLPK_constants.jl")
run(`$CC -E -dM $glpkheader` | `perl -nle $PERL_CMD` | (`sort` > glpk_h_jl))

run(`rm -fr $glpkname`)
rm(glpkarchive)
