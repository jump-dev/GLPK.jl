# This script is used to auto-generate the file GLPK_constants.jl
# It is intended for use by the Package maintaners whenever the GLPK
# version is updated.

using BinDeps

include("../deps/verreq.jl")

function with_temp_dir(fn::Function)
    tmpdir = mktempdir()
    try
        cd(fn, tmpdir)
    finally
        rm(tmpdir, recursive=true)
    end
end

function parse_header(iname, oname)
    CC = get(ENV, "CC", "cc")
    glpk_consts = Tuple{AbstractString, Cint}[]
    major_ver = 0
    minor_ver = 0
    for line in eachline(`$CC -E -dM $iname`)
        m = match(r"^#define\s+GLP_(\w+)\s+([^\s]+)", line)
        if m != nothing
            name, value = m.captures
            value = parse(Cint, value)
            if name == "MAJOR_VERSION"
                major_ver = value
            elseif name == "MINOR_VERSION"
                minor_ver = value
            else
                push!(glpk_consts, (name, value))
            end
        end
    end
    check_glpk_version(major_ver, minor_ver)
    sort!(glpk_consts, by=(key)->key[1])
    open(oname, "w") do hdl
        for (name::AbstractString, value::Cint) in glpk_consts
            write(hdl, "const $name = convert(Cint, $value)\n")
        end
    end
end

function get_header()
    if length(ARGS) >= 1 && isfile(ARGS[1])
        return ARGS[1]
    end
    glpkname = "glpk-$glpkdefver"
    glpkarchive = "$glpkname.tar.gz"

    run(download_cmd("http://ftp.gnu.org/gnu/glpk/$glpkarchive", glpkarchive))
    run(`tar xzf $glpkarchive`)

    joinpath(glpkname, "src", "glpk.h")
end

with_temp_dir() do
    glpkheader = get_header()

    glpk_h_jl = joinpath(dirname(@__FILE__), "..", "src", "GLPK_constants.jl")
    parse_header(glpkheader, glpk_h_jl)
end
