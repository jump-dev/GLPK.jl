# Automatically generated using Clang.jl


const GLP_MAJOR_VERSION = 4
const GLP_MINOR_VERSION = 64
const GLP_MIN = 1
const GLP_MAX = 2
const GLP_CV = 1
const GLP_IV = 2
const GLP_BV = 3
const GLP_FR = 1
const GLP_LO = 2
const GLP_UP = 3
const GLP_DB = 4
const GLP_FX = 5
const GLP_BS = 1
const GLP_NL = 2
const GLP_NU = 3
const GLP_NF = 4
const GLP_NS = 5
const GLP_SF_GM = 0x01
const GLP_SF_EQ = 0x10
const GLP_SF_2N = 0x20
const GLP_SF_SKIP = 0x40
const GLP_SF_AUTO = 0x80
const GLP_SOL = 1
const GLP_IPT = 2
const GLP_MIP = 3
const GLP_UNDEF = 1
const GLP_FEAS = 2
const GLP_INFEAS = 3
const GLP_NOFEAS = 4
const GLP_OPT = 5
const GLP_UNBND = 6
const GLP_BF_LUF = 0x00
const GLP_BF_BTF = 0x10
const GLP_BF_FT = 0x01
const GLP_BF_BG = 0x02
const GLP_BF_GR = 0x03
const GLP_MSG_OFF = 0
const GLP_MSG_ERR = 1
const GLP_MSG_ON = 2
const GLP_MSG_ALL = 3
const GLP_MSG_DBG = 4
const GLP_PRIMAL = 1
const GLP_DUALP = 2
const GLP_DUAL = 3
const GLP_PT_STD = 0x11
const GLP_PT_PSE = 0x22
const GLP_RT_STD = 0x11
const GLP_RT_HAR = 0x22
const GLP_RT_FLIP = 0x33
const GLP_USE_AT = 1
const GLP_USE_NT = 2
const GLP_ORD_NONE = 0
const GLP_ORD_QMD = 1
const GLP_ORD_AMD = 2
const GLP_ORD_SYMAMD = 3
const GLP_BR_FFV = 1
const GLP_BR_LFV = 2
const GLP_BR_MFV = 3
const GLP_BR_DTH = 4
const GLP_BR_PCH = 5
const GLP_BT_DFS = 1
const GLP_BT_BFS = 2
const GLP_BT_BLB = 3
const GLP_BT_BPH = 4
const GLP_PP_NONE = 0
const GLP_PP_ROOT = 1
const GLP_PP_ALL = 2
const GLP_RF_REG = 0
const GLP_RF_LAZY = 1
const GLP_RF_CUT = 2
const GLP_RF_GMI = 1
const GLP_RF_MIR = 2
const GLP_RF_COV = 3
const GLP_RF_CLQ = 4
const GLP_ON = 1
const GLP_OFF = 0
const GLP_IROWGEN = 0x01
const GLP_IBINGO = 0x02
const GLP_IHEUR = 0x03
const GLP_ICUTGEN = 0x04
const GLP_IBRANCH = 0x05
const GLP_ISELECT = 0x06
const GLP_IPREPRO = 0x07
const GLP_NO_BRNCH = 0
const GLP_DN_BRNCH = 1
const GLP_UP_BRNCH = 2
const GLP_EBADB = 0x01
const GLP_ESING = 0x02
const GLP_ECOND = 0x03
const GLP_EBOUND = 0x04
const GLP_EFAIL = 0x05
const GLP_EOBJLL = 0x06
const GLP_EOBJUL = 0x07
const GLP_EITLIM = 0x08
const GLP_ETMLIM = 0x09
const GLP_ENOPFS = 0x0a
const GLP_ENODFS = 0x0b
const GLP_EROOT = 0x0c
const GLP_ESTOP = 0x0d
const GLP_EMIPGAP = 0x0e
const GLP_ENOFEAS = Float32(0x00)
const GLP_ENOCVG = 0x10
const GLP_EINSTAB = 0x11
const GLP_EDATA = 0x12
const GLP_ERANGE = 0x13
const GLP_KKT_PE = 1
const GLP_KKT_PB = 2
const GLP_KKT_DE = 3
const GLP_KKT_DB = 4
const GLP_KKT_CS = 5
const GLP_MPS_DECK = 1
const GLP_MPS_FILE = 2

# Skipping MacroDefinition: glp_error glp_error_ ( __FILE__ , __LINE__ )
# Skipping MacroDefinition: glp_assert ( expr ) ( ( void ) ( ( expr ) || ( glp_assert_ ( # expr , __FILE__ , __LINE__ ) , 1 ) ) )
# Skipping MacroDefinition: glp_malloc ( size ) glp_alloc ( 1 , size )
# Skipping MacroDefinition: glp_calloc ( n , size ) glp_alloc ( n , size )

const GLP_ASN_MIN = 1
const GLP_ASN_MAX = 2
const GLP_ASN_MMP = 3
const glp_prob = Cvoid

mutable struct glp_bfcp
    msg_lev::Cint
    type::Cint
    lu_size::Cint
    piv_tol::Cdouble
    piv_lim::Cint
    suhl::Cint
    eps_tol::Cdouble
    max_gro::Cdouble
    nfs_max::Cint
    upd_tol::Cdouble
    nrs_max::Cint
    rs_size::Cint
    foo_bar::NTuple{38, Cdouble}
end

Base.cconvert(::Type{Ptr{glp_bfcp}}, x::glp_bfcp) = x
function Base.unsafe_convert(::Type{Ptr{glp_bfcp}}, x::glp_bfcp)
    return convert(Ptr{glp_bfcp}, pointer_from_objref(x))
end

mutable struct glp_smcp
    msg_lev::Cint
    meth::Cint
    pricing::Cint
    r_test::Cint
    tol_bnd::Cdouble
    tol_dj::Cdouble
    tol_piv::Cdouble
    obj_ll::Cdouble
    obj_ul::Cdouble
    it_lim::Cint
    tm_lim::Cint
    out_frq::Cint
    out_dly::Cint
    presolve::Cint
    excl::Cint
    shift::Cint
    aorn::Cint
    foo_bar::NTuple{33, Cdouble}
    glp_smcp() = new()
end

Base.cconvert(::Type{Ptr{glp_smcp}}, x::glp_smcp) = x
function Base.unsafe_convert(::Type{Ptr{glp_smcp}}, x::glp_smcp)
    return convert(Ptr{glp_smcp}, pointer_from_objref(x))
end

mutable struct glp_iptcp
    msg_lev::Cint
    ord_alg::Cint
    foo_bar::NTuple{48, Cdouble}
    glp_iptcp() = new()
end

Base.cconvert(::Type{Ptr{glp_iptcp}}, x::glp_iptcp) = x
function Base.unsafe_convert(::Type{Ptr{glp_iptcp}}, x::glp_iptcp)
    return convert(Ptr{glp_iptcp}, pointer_from_objref(x))
end

const glp_tree = Cvoid

mutable struct glp_iocp
    msg_lev::Cint
    br_tech::Cint
    bt_tech::Cint
    tol_int::Cdouble
    tol_obj::Cdouble
    tm_lim::Cint
    out_frq::Cint
    out_dly::Cint
    cb_func::Ptr{Cvoid}
    cb_info::Ptr{Cvoid}
    cb_size::Cint
    pp_tech::Cint
    mip_gap::Cdouble
    mir_cuts::Cint
    gmi_cuts::Cint
    cov_cuts::Cint
    clq_cuts::Cint
    presolve::Cint
    binarize::Cint
    fp_heur::Cint
    ps_heur::Cint
    ps_tm_lim::Cint
    sr_heur::Cint
    use_sol::Cint
    save_sol::Cstring
    alien::Cint
    flip::Cint
    foo_bar::NTuple{23, Cdouble}
    glp_iocp() = new()
end

Base.cconvert(::Type{Ptr{glp_iocp}}, x::glp_iocp) = x
function Base.unsafe_convert(::Type{Ptr{glp_iocp}}, x::glp_iocp)
    return convert(Ptr{glp_iocp}, pointer_from_objref(x))
end

mutable struct glp_attr
    level::Cint
    origin::Cint
    klass::Cint
    foo_bar::NTuple{7, Cdouble}
end

Base.cconvert(::Type{Ptr{glp_attr}}, x::glp_attr) = x
function Base.unsafe_convert(::Type{Ptr{glp_attr}}, x::glp_attr)
    return convert(Ptr{glp_attr}, pointer_from_objref(x))
end

mutable struct glp_mpscp
    blank::Cint
    obj_name::Cstring
    tol_mps::Cdouble
    foo_bar::NTuple{17, Cdouble}
    glp_mpscp() = new()
end

Base.cconvert(::Type{Ptr{glp_mpscp}}, x::glp_mpscp) = x
function Base.unsafe_convert(::Type{Ptr{glp_mpscp}}, x::glp_mpscp)
    return convert(Ptr{glp_mpscp}, pointer_from_objref(x))
end

mutable struct glp_cpxcp
    foo_bar::NTuple{20, Cdouble}
    glp_cpxcp() = new()
end

Base.cconvert(::Type{Ptr{glp_cpxcp}}, x::glp_cpxcp) = x
function Base.unsafe_convert(::Type{Ptr{glp_cpxcp}}, x::glp_cpxcp)
    return convert(Ptr{glp_cpxcp}, pointer_from_objref(x))
end

const glp_tran = Cvoid
const glp_errfunc = Ptr{Cvoid}

mutable struct glp_arc
    tail::Ptr
    head::Ptr
    data::Ptr{Cvoid}
    temp::Ptr{Cvoid}
    t_prev::Ptr{glp_arc}
    t_next::Ptr{glp_arc}
    h_prev::Ptr{glp_arc}
    h_next::Ptr{glp_arc}
end

Base.cconvert(::Type{Ptr{glp_arc}}, x::glp_arc) = x
function Base.unsafe_convert(::Type{Ptr{glp_arc}}, x::glp_arc)
    return convert(Ptr{glp_arc}, pointer_from_objref(x))
end

mutable struct glp_vertex
    i::Cint
    name::Cstring
    entry::Ptr{Cvoid}
    data::Ptr{Cvoid}
    temp::Ptr{Cvoid}
    in::Ptr{glp_arc}
    out::Ptr{glp_arc}
end

Base.cconvert(::Type{Ptr{glp_vertex}}, x::glp_vertex) = x
function Base.unsafe_convert(::Type{Ptr{glp_vertex}}, x::glp_vertex)
    return convert(Ptr{glp_vertex}, pointer_from_objref(x))
end

mutable struct glp_graph
    pool::Ptr{Cvoid}
    name::Cstring
    nv_max::Cint
    nv::Cint
    na::Cint
    v::Ptr{Ptr{glp_vertex}}
    index::Ptr{Cvoid}
    v_size::Cint
    a_size::Cint
end

Base.cconvert(::Type{Ptr{glp_graph}}, x::glp_graph) = x
function Base.unsafe_convert(::Type{Ptr{glp_graph}}, x::glp_graph)
    return convert(Ptr{glp_graph}, pointer_from_objref(x))
end
