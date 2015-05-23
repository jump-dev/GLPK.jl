# Version requirements for libglpk

const glpkminver = v"4.52" # Minimum requirement
const glpkmaxver = v"4.55" # Maximum version known to work
const glpkdefver = "4.52" # Default version

function check_glpk_version(major_ver, minor_ver)
    if !(glpkminver <= VersionNumber(major_ver, minor_ver) <= glpkmaxver)
        error("""
              GLPK: libglpk version not supported:
              Requires: >=$glpkminver, <=$glpkmaxver
              Version found: $major_ver.$minor_ver
              """)
    end
end
