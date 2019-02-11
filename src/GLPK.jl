###
### GLPK API Wrapper
###

__precompile__()

module GLPK

using Compat
using Compat.SparseArrays

## Exports
#{{{
export
    # Types
    Param,
    GLPKError,
    GLPKFatalError,
    Prob,
    SimplexParam,
    InteriorParam,
    IntoptParam,
    BasisFactParam,
    MathProgWorkspace,

    # Methods
    version,
    set_prob_name,
    set_obj_name,
    set_row_name,
    set_col_name,
    set_obj_dir,
    add_rows,
    add_cols,
    set_row_bnds,
    set_col_bnds,
    set_obj_coef,
    set_mat_row,
    set_mat_col,
    load_matrix,
    check_dup,
    sort_matrix,
    del_rows,
    del_cols,
    copy_prob,
    erase_prob,
    get_prob_name,
    get_obj_name,
    get_obj_dir,
    get_num_rows,
    get_num_cols,
    get_row_name,
    get_col_name,
    get_row_type,
    get_row_lb,
    get_row_ub,
    get_col_type,
    get_col_lb,
    get_col_ub,
    get_obj_coef,
    get_num_nz,
    get_mat_row,
    get_mat_col,
    create_index,
    find_row,
    find_col,
    delete_index,
    set_rii,
    set_sjj,
    get_rii,
    get_sjj,
    scale_prob,
    unscale_prob,
    set_row_stat,
    set_col_stat,
    std_basis,
    adv_basis,
    cpx_basis,
    simplex,
    exact,
    init_smcp,
    get_status,
    get_prim_stat,
    get_dual_stat,
    get_obj_val,
    get_row_stat,
    get_row_prim,
    get_row_dual,
    get_col_stat,
    get_col_prim,
    get_col_dual,
    get_unbnd_ray,
    interior,
    init_iptcp,
    ipt_status,
    ipt_obj_val,
    ipt_row_prim,
    ipt_row_dual,
    ipt_col_prim,
    ipt_col_dual,
    set_col_kind,
    get_col_kind,
    get_num_int,
    get_num_bin,
    intopt,
    init_iocp,
    mip_status,
    mip_obj_val,
    mip_row_val,
    mip_col_val,
    check_kkt,
    read_mps,
    write_mps,
    read_lp,
    write_lp,
    read_prob,
    write_prob,
    mpl_read_model,
    mpl_read_data,
    mpl_generate,
    mpl_build_prob,
    mpl_postsolve,
    print_sol,
    read_sol,
    write_sol,
    print_ipt,
    read_ipt,
    write_ipt,
    print_mip,
    read_mip,
    write_mip,
    print_ranges,
    print_ranges,
    bf_exists,
    factorize,
    bf_updated,
    get_bfcp,
    set_bfcp,
    get_bhead,
    get_row_bind,
    get_col_bind,
    ftran,
    btran,
    warm_up,
    eval_tab_row,
    eval_tab_col,
    transform_row,
    transform_col,
    prim_rtest,
    dual_rtest,
    analyze_bound,
    analyze_coef,
    ios_reason,
    ios_get_prob,
    ios_row_attr,
    ios_mip_gap,
    ios_node_data,
    ios_select_node,
    ios_heur_sol,
    ios_can_branch,
    ios_branch_upon,
    ios_terminate,
    ios_tree_size,
    ios_curr_node,
    ios_next_node,
    ios_prev_node,
    ios_up_node,
    ios_node_level,
    ios_node_bound,
    ios_best_node,
    ios_pool_size,
    ios_add_row,
    ios_del_row,
    ios_clear_pool,
    init_env,
    free_env,
    term_out,
    open_tee,
    close_tee,
    malloc,
    calloc,
    free,
    mem_usage,
    mem_limit,
    read_cnfsat,
    check_cnfsat,
    write_cnfsat,
    minisat1,
    intfeas1
#}}}

import Base.setindex!, Base.getindex

## Shared library interface setup
#{{{
if isfile(joinpath(dirname(@__FILE__),"..","deps","deps.jl"))
    include("../deps/deps.jl")
else
    error("GLPK not properly installed. Please run Pkg.build(\"GLPK\")")
end

include("GLPK_constants.jl")

# General recoverable exception: all GLPK functions
# throw this in case of recoverable errors
mutable struct GLPKError <: Exception
    msg::AbstractString
end

# Fatal exception: when this is thrown, all GLPK
# objects are no longer valid
mutable struct GLPKFatalError <: Exception
    msg::AbstractString
end

# Error hook, used to catch internal errors when calling
# GLPK functions
function _err_hook(info::Ptr{Cvoid})
    ccall((:glp_error_hook, libglpk), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), C_NULL, C_NULL)
    ccall((:glp_free_env, libglpk), Cvoid, ())
    _del_all_objs()
    throw(GLPKFatalError("GLPK call failed. All GLPK objects you defined so far are now invalidated."))
end

macro glpk_ccall(f, args...)
    quote
        ccall((:glp_error_hook, libglpk), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), @cfunction(_err_hook, Cvoid, (Ptr{Cvoid},)), C_NULL)
        ret = ccall(($"glp_$f", libglpk), $(map(esc,args)...))
        ccall((:glp_error_hook, libglpk), Cvoid, (Ptr{Cvoid}, Ptr{Cvoid}), C_NULL, C_NULL)
        ret
    end
end

# We need to define GLPK.version as first thing
# in order to perform a sanity check
# (since we import structs from the header,
# we must ensure that the binary is the correct
# one)
function version()
    vstr = unsafe_string(@glpk_ccall version Ptr{Cchar} ())
    return tuple(map(x->parse(Int, x), split(vstr, '.'))...)
end

include("../deps/verreq.jl")

function __init__()
    major_ver, minor_ver = version()
    check_glpk_version(major_ver, minor_ver)
end

#}}}

## Preliminary definitions
#{{{

# General structure for the parameters types

abstract type Param end

function setindex!(param::T, val, field_name::AbstractString) where T<:Param
    s = Symbol(field_name)
    i = findfirst(x->x==s, fieldnames(T))
    i > 0 || error("Parameter type $T has no field $field_name")
    t = T.types[i]
    setfield!(param, s, convert(t, val))
end

function getindex(param::Param, field_name::AbstractString)
    getfield(param, Symbol(field_name))
end

# We define some types which allow to pass optional agruments
# to the function.
# In this framework, optional arguments can be passed either
# as an empty vector [] or as the 'nothing' constant

const VecOrNothing = Union{AbstractVector,Nothing}
cint_vec(a::Nothing) = Cint[]
cint_vec(a::AbstractVector) = convert(Vector{Cint}, a)

cdouble_vec(a::Nothing) = Cdouble[]
cdouble_vec(a::AbstractVector) = convert(Vector{Cdouble}, a)

vecornothing_length(a::Nothing) = 0
vecornothing_length(a) = length(a)

#}}}


## Main types definitions
#{{{
# All structs in original glpk are wrapped up in
# composite types, which initialize and destroy themselves
# as needed, and expose pointers when asked to by
# ccall's.
#
# Therefore, the original C glp API
#
#  int glp_simplex(prob * lp, glp_smpc * param)
#
# becomes
#
#  GLPK.simplex(lp::GLPK.Prob, param::GLPK.SimplexParam)
#
#
# The map between names is as follows:
#
# +-------------+--------------------------+
# |  C          |  Julia                   |
# +-------------+--------------------------+
# |  glp_prob   |  GLPK.Prob               |
# |  glp_smcp   |  GLPK.SimplexParam       |
# |  glp_iptcp  |  GLPK.InteriorParam      |
# |  glp_iocp   |  GLPK.IntoptParam        |
# |  glp_bfcp   |  GLPK.BasisFactParam     |
# |  glp_tran   |  GLPK.MathProgWorkspace  |
# |  glp_attr   |  GLPK.Attr               |
# +-------------+--------------------------+
#
# In order to get/set the value of a Param field, no
# special syntax is required:
#
#   lps_opts = GLPK.SimplexParam()
#   lps_opts.msg_lev = GLPK.MSG_ERR
#   lps_opts.presolve = GLPK.ON
#
# However, as a special case, the "type" field in
# glp_bfcp is renamed as "bftype" in GLPK.BasisFactParam:
#
#   bf_opts = GLPK.BasisFactParam()
#   bf_opts.bftype = ...
#

mutable struct Prob
    p::Ptr{Cvoid}
    function Prob(p::Ptr{Cvoid})
        create = (p == C_NULL)
        if create
            p = @glpk_ccall create_prob Ptr{Cvoid} ()
        end
        prob = new(p)
        _add_obj(prob)
        if create
            if VERSION >= v"0.7-"
                finalizer(delete_prob, prob)
            else
                finalizer(prob, delete_prob)
            end
        else
            if VERSION >= v"0.7-"
                finalizer(_del_obj, prob)
            else
                finalizer(prob, _del_obj)
            end
        end
        return prob
    end
end

Prob() = Prob(C_NULL)

function delete_prob(prob::Prob)
    prob.p == C_NULL && return
    if jl_obj_is_valid(prob)
        @glpk_ccall delete_prob Cvoid (Ptr{Cvoid},) prob.p
        _del_obj(prob)
    end
    prob.p = C_NULL
    return
end

mutable struct SimplexParam <: Param
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
    _reserved01::Cdouble
    _reserved02::Cdouble
    _reserved03::Cdouble
    _reserved04::Cdouble
    _reserved05::Cdouble
    _reserved06::Cdouble
    _reserved07::Cdouble
    _reserved08::Cdouble
    _reserved09::Cdouble
    _reserved10::Cdouble
    _reserved11::Cdouble
    _reserved12::Cdouble
    _reserved13::Cdouble
    _reserved14::Cdouble
    _reserved15::Cdouble
    _reserved16::Cdouble
    _reserved17::Cdouble
    _reserved18::Cdouble
    _reserved19::Cdouble
    _reserved20::Cdouble
    _reserved21::Cdouble
    _reserved22::Cdouble
    _reserved23::Cdouble
    _reserved24::Cdouble
    _reserved25::Cdouble
    _reserved26::Cdouble
    _reserved27::Cdouble
    _reserved28::Cdouble
    _reserved29::Cdouble
    _reserved30::Cdouble
    _reserved31::Cdouble
    _reserved32::Cdouble
    _reserved33::Cdouble
    _reserved34::Cdouble
    _reserved35::Cdouble
    _reserved36::Cdouble

    function SimplexParam()
        p = new()
        @glpk_ccall "init_smcp" Cvoid (Ref{SimplexParam},) p
        return p
    end
end

mutable struct InteriorParam <: Param
    msg_lev::Cint
    ord_alg::Cint
    _reserved01::Cdouble
    _reserved02::Cdouble
    _reserved03::Cdouble
    _reserved04::Cdouble
    _reserved05::Cdouble
    _reserved06::Cdouble
    _reserved07::Cdouble
    _reserved08::Cdouble
    _reserved09::Cdouble
    _reserved10::Cdouble
    _reserved11::Cdouble
    _reserved12::Cdouble
    _reserved13::Cdouble
    _reserved14::Cdouble
    _reserved15::Cdouble
    _reserved16::Cdouble
    _reserved17::Cdouble
    _reserved18::Cdouble
    _reserved19::Cdouble
    _reserved20::Cdouble
    _reserved21::Cdouble
    _reserved22::Cdouble
    _reserved23::Cdouble
    _reserved24::Cdouble
    _reserved25::Cdouble
    _reserved26::Cdouble
    _reserved27::Cdouble
    _reserved28::Cdouble
    _reserved29::Cdouble
    _reserved30::Cdouble
    _reserved31::Cdouble
    _reserved32::Cdouble
    _reserved33::Cdouble
    _reserved34::Cdouble
    _reserved35::Cdouble
    _reserved36::Cdouble
    _reserved37::Cdouble
    _reserved38::Cdouble
    _reserved39::Cdouble
    _reserved40::Cdouble
    _reserved41::Cdouble
    _reserved42::Cdouble
    _reserved43::Cdouble
    _reserved44::Cdouble
    _reserved45::Cdouble
    _reserved46::Cdouble
    _reserved47::Cdouble
    _reserved48::Cdouble

    function InteriorParam()
        p = new()
        @glpk_ccall "init_iptcp" Cvoid (Ref{InteriorParam},) p
        return p
    end
end

mutable struct IntoptParam <: Param
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
    use_sol::Cint
    save_sol::Ptr{Cchar}
    alien::Cint
    _reserved01::Cdouble
    _reserved02::Cdouble
    _reserved03::Cdouble
    _reserved04::Cdouble
    _reserved05::Cdouble
    _reserved06::Cdouble
    _reserved07::Cdouble
    _reserved08::Cdouble
    _reserved09::Cdouble
    _reserved10::Cdouble
    _reserved11::Cdouble
    _reserved12::Cdouble
    _reserved13::Cdouble
    _reserved14::Cdouble
    _reserved15::Cdouble
    _reserved16::Cdouble
    _reserved17::Cdouble
    _reserved18::Cdouble
    _reserved19::Cdouble
    _reserved20::Cdouble
    _reserved21::Cdouble
    _reserved22::Cdouble
    _reserved23::Cdouble
    _reserved24::Cdouble
    _reserved25::Cdouble

    function IntoptParam()
        p = new()
        @glpk_ccall "init_iocp" Cvoid (Ref{IntoptParam},) p
        return p
    end
end

mutable struct BasisFactParam <: Param
    msg_lev::Cint
    bftype::Cint # NOTE: changed type->bftype
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
    _reserved01::Cdouble
    _reserved02::Cdouble
    _reserved03::Cdouble
    _reserved04::Cdouble
    _reserved05::Cdouble
    _reserved06::Cdouble
    _reserved07::Cdouble
    _reserved08::Cdouble
    _reserved09::Cdouble
    _reserved10::Cdouble
    _reserved11::Cdouble
    _reserved12::Cdouble
    _reserved13::Cdouble
    _reserved14::Cdouble
    _reserved15::Cdouble
    _reserved16::Cdouble
    _reserved17::Cdouble
    _reserved18::Cdouble
    _reserved19::Cdouble
    _reserved20::Cdouble
    _reserved21::Cdouble
    _reserved22::Cdouble
    _reserved23::Cdouble
    _reserved24::Cdouble
    _reserved25::Cdouble
    _reserved26::Cdouble
    _reserved27::Cdouble
    _reserved28::Cdouble
    _reserved29::Cdouble
    _reserved30::Cdouble
    _reserved31::Cdouble
    _reserved32::Cdouble
    _reserved33::Cdouble
    _reserved34::Cdouble
    _reserved35::Cdouble
    _reserved36::Cdouble
    _reserved37::Cdouble
    _reserved38::Cdouble

    function BasisFactParam()
        return new(0, 0, 0, 0.0, 0, 0, 0.0, 0.0, 0, 0.0, 0, 0)
    end
end

mutable struct MathProgWorkspace
    p::Ptr{Cvoid}
    function MathProgWorkspace()
        tran = @glpk_ccall mpl_alloc_wksp Ptr{Cvoid} ()
        wksp = new(tran)
        _add_obj(wksp)
        if VERSION >= v"0.7-"
            finalizer(GLPK.mpl_free_wksp, wksp)
        else
            finalizer(wksp, GLPK.mpl_free_wksp)
        end
        return wksp
    end
end

function mpl_free_wksp(tran::MathProgWorkspace)
    tran.p == C_NULL && return
    if jl_obj_is_valid(tran)
        @glpk_ccall mpl_free_wksp Cvoid (Ptr{Cvoid},) tran.p
        _del_obj(tran)
    end
    tran.p = C_NULL
    return
end

mutable struct Attr
    level::Cint
    origin::Cint
    klass::Cint
    _reserved01::Cdouble
    _reserved02::Cdouble
    _reserved03::Cdouble
    _reserved04::Cdouble
    _reserved05::Cdouble
    _reserved06::Cdouble
    _reserved07::Cdouble
    function Attr()
        new(0, 0, 0)
    end
end
#}}}

## GLP functions
#{{{
# The API interface is as close as possible to the original
# one.
# The general translation rules are:
#
#  * function names translate like this: glp_func -> GLPK.func
#  * constant names translate like this: GLPK_CONST -> GLPK.CONST
#  * whenever the C library accepts NULL as argument,
#    the Julia one will accept the nothing constant.
#  * vectors do not need to have an extra element at the
#    beginning to accomodate to the 1-based GLP indexing.
#  * most functions will accept any kind of vectors as inputs,
#    provided they can be converted to be C-compatible
#    (i.e. to either Cint or Cdouble elements).
#  * the exceptions to the above are those functions which write
#    their output in a vector, in which case the vector type must
#    be strictly C-compatible.
#  * all char[] strings become Strings, both in inputs and in output.
#
#  A single exception to the strict compatibility is GLPK.version(),
#  which returns a tuple of integers in the form (major, minor)
#  rather than a string.

include("GLPK_checks.jl")

function set_prob_name(prob::Prob, name::Union{AbstractString,Nothing})
    @check! _prob(prob)
    name == nothing && (name = "")
    @check _string_length(name, 255)
    @glpk_ccall set_prob_name Cvoid (Ptr{Cvoid}, Ptr{Cchar}) prob.p string(name)
end

function set_obj_name(prob::Prob, name::Union{AbstractString,Nothing})
    @check! _prob(prob)
    name == nothing && (name = "")
    @check _string_length(name, 255)
    @glpk_ccall set_obj_name Cvoid (Ptr{Cvoid}, Ptr{Cchar}) prob.p string(name)
end

function set_row_name(prob::Prob, row::Integer, name::Union{AbstractString,Nothing})
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    name == nothing && (name = "")
    @check _string_length(name, 255)
    @glpk_ccall set_row_name Cvoid (Ptr{Cvoid}, Cint, Ptr{Cchar}) prob.p row string(name)
end

function set_col_name(prob::Prob, col::Integer, name::Union{AbstractString,Nothing})
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    name == nothing && (name = "")
    @check _string_length(name, 255)
    @glpk_ccall set_col_name Cvoid (Ptr{Cvoid}, Cint, Ptr{Cchar}) prob.p col string(name)
end

function set_obj_dir(prob::Prob, dir::Integer)
    @check! _prob(prob)
    @check _obj_dir_is_valid(dir)
    @glpk_ccall set_obj_dir Cvoid (Ptr{Cvoid}, Cint) prob.p dir
end

function add_rows(prob::Prob, rows::Integer)
    @check! _prob(prob)
    @glpk_ccall add_rows Cint (Ptr{Cvoid}, Cint) prob.p rows
end

function add_cols(prob::Prob, cols::Integer)
    @check! _prob(prob)
    @glpk_ccall add_cols Cint (Ptr{Cvoid}, Cint) prob.p cols
end

function set_row_bnds(prob::Prob, row::Integer, bounds_type::Integer, lb::Real, ub::Real)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    @check _bounds_type_is_valid(bounds_type)
    @check _bounds_are_valid(bounds_type, lb, ub, "constraint")
    @glpk_ccall set_row_bnds Cvoid (Ptr{Cvoid}, Cint, Cint, Cdouble, Cdouble) prob.p row bounds_type lb ub
end

function set_col_bnds(prob::Prob, col::Integer, bounds_type::Integer, lb::Real, ub::Real)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @check _bounds_type_is_valid(bounds_type)
    @check _bounds_are_valid(bounds_type, lb, ub, "variable")
    @glpk_ccall set_col_bnds Cvoid (Ptr{Cvoid}, Cint, Cint, Cdouble, Cdouble) prob.p col bounds_type lb ub
end

function set_obj_coef(prob::Prob, col::Integer, coef::Real)
    @check! _prob(prob)
    @check _col_is_valid_w0(prob, col)
    @glpk_ccall set_obj_coef Cvoid (Ptr{Cvoid}, Cint, Cdouble) prob.p col coef
end

function set_mat_row(prob::Prob, row::Integer, len::Integer, ind::VecOrNothing, val::VecOrNothing)
    @check! _prob(prob)
    @check! _vectors_size(len, ind, val)
    @check _row_is_valid(prob, row)
    if len > 0
        ind32 = cint_vec(ind)
        val64 = cdouble_vec(val)
        off32 = sizeof(Cint)
        off64 = sizeof(Cdouble)
        ind32p = pointer(ind32) - off32
        val64p = pointer(val64) - off64
        @check! _cols_ids_size(prob, 0, len, ind32)
        @check _cols_ids_content(prob, len, ind32)
    else
        ind32p = C_NULL
        val64p = C_NULL
    end

    @glpk_ccall set_mat_row Cvoid (Ptr{Cvoid}, Cint, Cint, Ptr{Cint}, Ptr{Cdouble}) prob.p row len ind32p val64p
end

function set_mat_row(prob::Prob, row::Integer, ind::VecOrNothing, val::VecOrNothing)
    @check! _vectors_all_same_size(ind, val)
    l = vecornothing_length(ind)
    set_mat_row(prob, row, l, ind, val)
end


function set_mat_col(prob::Prob, col::Integer, len::Integer, ind::VecOrNothing, val::VecOrNothing)
    @check! _prob(prob)
    @check! _vectors_size(len, ind, val)
    @check _col_is_valid(prob, col)
    if len > 0
        ind32 = cint_vec(ind)
        val64 = cdouble_vec(val)
        off32 = sizeof(Cint)
        off64 = sizeof(Cdouble)
        ind32p = pointer(ind32) - off32
        val64p = pointer(val64) - off64
        @check! _rows_ids_size(prob, 0, len, ind32)
        @check _rows_ids_content(prob, len, ind32)
    else
        ind32p = C_NULL
        val64p = C_NULL
    end

    @glpk_ccall set_mat_col Cvoid (Ptr{Cvoid}, Cint, Cint, Ptr{Cint}, Ptr{Cdouble}) prob.p col len ind32p val64p
end

function set_mat_col(prob::Prob, col::Integer, ind::VecOrNothing, val::VecOrNothing)
    @check! _vectors_all_same_size(ind, val)
    l = vecornothing_length(ind)
    set_mat_col(prob, col, l, ind, val)
end

function load_matrix(prob::Prob, numel::Integer, ia::VecOrNothing, ja::VecOrNothing, ar::VecOrNothing)
    @check! _prob(prob)
    @check! _vectors_size(numel, ia, ja, ar)
    numel == 0 && return
    ia32 = cint_vec(ia)
    ja32 = cint_vec(ja)
    ar64 = cdouble_vec(ar)
    @check _indices_vectors_dup(prob, numel, ia32, ja32)

    off32 = sizeof(Cint)
    off64 = sizeof(Cdouble)
    ia32p = pointer(ia32) - off32
    ja32p = pointer(ja32) - off32
    ar64p = pointer(ar64) - off64

    @glpk_ccall load_matrix Cvoid (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}) prob.p numel ia32p ja32p ar64p
end

function load_matrix(prob::Prob, ia::VecOrNothing, ja::VecOrNothing, ar::VecOrNothing)
    @check! _vectors_all_same_size(ia, ja, ar)
    l = vecornothing_length(ar)
    load_matrix(prob, l, ia, ja, ar)
end

function load_matrix(prob::Prob, a::Compat.SparseArrays.AbstractSparseMatrix)
    (ia, ja, ar) = findnz(a)
    load_matrix(prob, ia, ja, ar)
end

function check_dup(rows::Integer, cols::Integer, numel::Integer, ia::VecOrNothing, ja::VecOrNothing)
    @check _rows_and_cols(rows, cols)
    @check! _vectors_size(numel, ia, ja)
    ia32 = cint_vec(ia)
    ja32 = cint_vec(ja)

    off32 = sizeof(Cint)
    ia32p = pointer(ia32) - off32
    ja32p = pointer(ja32) - off32

    @glpk_ccall check_dup Cint (Cint, Cint, Cint, Ptr{Cint}, Ptr{Cint}) rows cols numel ia32p ja32p
end

function check_dup(rows::Integer, cols::Integer, ia::VecOrNothing, ja::VecOrNothing)
    @check! _vectors_all_same_size(ia, ja)
    l = vecornothing_length(ia)
    check_dup(rows, cols, l, ia, ja)
end

function sort_matrix(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall sort_matrix Cvoid (Ptr{Cvoid},) prob.p
end

function del_rows(prob::Prob, num_rows::Integer, rows_ids::VecOrNothing)
    @check! _prob(prob)
    rows_ids32 = cint_vec(rows_ids)
    @check! _rows_ids_size(prob, 1, num_rows, rows_ids32)
    @check _rows_ids_content(prob, num_rows, rows_ids32)

    off32 = sizeof(Cint)
    rows_ids32p = pointer(rows_ids32) - off32
    @glpk_ccall del_rows Cvoid (Ptr{Cvoid}, Cint, Ptr{Cint}) prob.p num_rows rows_ids32p
end
del_rows(prob::Prob, rows_ids::VecOrNothing) =
    del_rows(prob, length(rows_ids), rows_ids)

function del_cols(prob::Prob, num_cols::Integer, cols_ids::VecOrNothing)
    @check! _prob(prob)
    cols_ids32 = cint_vec(cols_ids)
    @check! _cols_ids_size(prob, 1, num_cols, cols_ids32)
    @check _cols_ids_content(prob, num_cols, cols_ids32)

    off32 = sizeof(Cint)
    cols_ids32p = pointer(cols_ids32) - off32
    @glpk_ccall del_cols Cvoid (Ptr{Cvoid}, Cint, Ptr{Cint}) prob.p num_cols cols_ids32p
end
del_cols(prob::Prob, cols_ids::VecOrNothing) =
    del_cols(prob, length(cols_ids), cols_ids)

function copy_prob(prob_dest::Prob, prob::Prob, copy_names::Integer)
    @check! _prob(prob)
    @check _copy_names_flag(copy_names)
    @glpk_ccall copy_prob Cvoid (Ptr{Cvoid}, Ptr{Cvoid}, Cint) prob_dest.p prob.p copy_names
end

function erase_prob(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall erase_prob Cvoid (Ptr{Cvoid},) prob.p
end

function get_prob_name(prob::Prob)
    @check! _prob(prob)
    name_cstr = @glpk_ccall get_prob_name Ptr{Cchar} (Ptr{Cvoid},) prob.p
    name_cstr == C_NULL && return ""
    return unsafe_string(name_cstr)
end

function get_obj_name(prob::Prob)
    @check! _prob(prob)
    name_cstr = @glpk_ccall get_obj_name Ptr{Cchar} (Ptr{Cvoid},) prob.p
    name_cstr == C_NULL && return ""
    return unsafe_string(name_cstr)
end

function get_obj_dir(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall get_obj_dir Cint (Ptr{Cvoid},) prob.p
end

function get_num_rows(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall get_num_rows Cint (Ptr{Cvoid},) prob.p
end

function get_num_cols(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall get_num_cols Cint (Ptr{Cvoid},) prob.p
end

function get_row_name(prob::Prob, row::Integer)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    name_cstr = @glpk_ccall get_row_name Ptr{Cchar} (Ptr{Cvoid}, Cint) prob.p row
    name_cstr == C_NULL && return ""
    return unsafe_string(name_cstr)
end

function get_col_name(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    name_cstr = @glpk_ccall get_col_name Ptr{Cchar} (Ptr{Cvoid}, Cint) prob.p col
    name_cstr == C_NULL && return ""
    return unsafe_string(name_cstr)
end

function get_row_type(prob::Prob, row::Integer)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    @glpk_ccall get_row_type Cint (Ptr{Cvoid}, Cint) prob.p row
end

function get_row_lb(prob::Prob, row::Integer)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    @glpk_ccall get_row_lb Cdouble (Ptr{Cvoid}, Cint) prob.p row
end

function get_row_ub(prob::Prob, row::Integer)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    @glpk_ccall get_row_ub Cdouble (Ptr{Cvoid}, Cint) prob.p row
end

function get_col_type(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @glpk_ccall get_col_type Cint (Ptr{Cvoid}, Cint) prob.p col
end

function get_col_lb(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @glpk_ccall get_col_lb Cdouble (Ptr{Cvoid}, Cint) prob.p col
end

function get_col_ub(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @glpk_ccall get_col_ub Cdouble (Ptr{Cvoid}, Cint) prob.p col
end

function get_obj_coef(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid_w0(prob, col)
    @glpk_ccall get_obj_coef Cdouble (Ptr{Cvoid}, Cint) prob.p col
end

function get_num_nz(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall get_num_nz Cint (Ptr{Cvoid},) prob.p
end

function get_mat_row(prob::Prob, row::Integer, ind::Union{Vector{Cint},Nothing}, val::Union{Vector{Cdouble},Nothing})
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    numel = @glpk_ccall get_mat_row Cint (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}) prob.p row C_NULL C_NULL
    numel == 0 && return 0

    if ind != nothing
        @check! _vectors_size(numel, ind)
        off32 = sizeof(Cint)
        ind32p = pointer(ind) - off32
    else
        ind32p = C_NULL
    end
    if val != nothing
        @check! _vectors_size(numel, val)
        off64 = sizeof(Cdouble)
        val64p = pointer(val) - off64
    else
        val64p = C_NULL
    end
    @glpk_ccall get_mat_row Cint (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}) prob.p row ind32p val64p
end

function get_mat_row(prob::Prob, row::Integer)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    numel = @glpk_ccall get_mat_row Cint (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}) prob.p row C_NULL C_NULL
    numel == 0 && return (Cint[], Cdouble[])

    ind = Array{Cint}(undef, numel)
    val = Array{Cdouble}(undef, numel)

    off32 = sizeof(Cint)
    ind32p = pointer(ind) - off32
    off64 = sizeof(Cdouble)
    val64p = pointer(val) - off64
    @glpk_ccall get_mat_row Cint (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}) prob.p row ind32p val64p
    return ind, val
end

function get_mat_col(prob::Prob, col::Integer, ind::Union{Vector{Cint},Nothing}, val::Union{Vector{Cdouble},Nothing})
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    numel = @glpk_ccall get_mat_col Cint (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}) prob.p col C_NULL C_NULL
    numel == 0 && return 0

    if ind != nothing
        @check! _vectors_size(numel, ind)
        off32 = sizeof(Cint)
        ind32p = pointer(ind) - off32
    else
        ind32p = C_NULL
    end
    if val != nothing
        @check! _vectors_size(numel, val)
        off64 = sizeof(Cdouble)
        val64p = pointer(val) - off64
    else
        val64p = C_NULL
    end
    @glpk_ccall get_mat_col Cint (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}) prob.p col ind32p val64p
end

function get_mat_col(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    numel = @glpk_ccall get_mat_col Cint (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}) prob.p col C_NULL C_NULL
    numel == 0 && return (Cint[], Cdouble[])

    ind = Array{Cint}(undef, numel)
    val = Array{Cdouble}(undef, numel)

    off32 = sizeof(Cint)
    ind32p = pointer(ind) - off32
    off64 = sizeof(Cdouble)
    val64p = pointer(val) - off64
    @glpk_ccall get_mat_col Cint (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}) prob.p col ind32p val64p
    return ind, val
end

function create_index(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall create_index Cvoid (Ptr{Cvoid},) prob.p
end

function find_row(prob::Prob, name::AbstractString)
    @check! _prob(prob)
    @glpk_ccall find_row Cint (Ptr{Cvoid}, Ptr{Cchar}) prob.p string(name)
end

function find_col(prob::Prob, name::AbstractString)
    @check! _prob(prob)
    @glpk_ccall find_col Cint (Ptr{Cvoid}, Ptr{Cchar}) prob.p string(name)
end

function delete_index(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall delete_index Cvoid (Ptr{Cvoid},) prob.p
end

function set_rii(prob::Prob, row::Integer, rii::Real)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    @glpk_ccall set_rii Cvoid (Ptr{Cvoid}, Cint, Cdouble) prob.p row rii
end

function set_sjj(prob::Prob, col::Integer, sjj::Real)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @glpk_ccall set_sjj Cvoid (Ptr{Cvoid}, Cint, Cdouble) prob.p col sjj
end

function get_rii(prob::Prob, row::Integer)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    @glpk_ccall get_rii Cdouble (Ptr{Cvoid}, Cint) prob.p row
end

function get_sjj(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @glpk_ccall get_sjj Cdouble (Ptr{Cvoid}, Cint) prob.p col
end

function scale_prob(prob::Prob, flags::Integer)
    @check! _prob(prob)
    @check _scale_flags(flags)
    @glpk_ccall scale_prob Cvoid (Ptr{Cvoid}, Cint) prob.p flags
end

function unscale_prob(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall unscale_prob Cvoid (Ptr{Cvoid},) prob.p
end

function set_row_stat(prob::Prob, row::Integer, stat::Integer)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    @check _stat_is_valid(stat)
    @glpk_ccall set_row_stat Cvoid (Ptr{Cvoid}, Cint, Cint) prob.p row stat
end

function set_col_stat(prob::Prob, col::Integer, stat::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @check _stat_is_valid(stat)
    @glpk_ccall set_col_stat Cvoid (Ptr{Cvoid}, Cint, Cint) prob.p col stat
end

function std_basis(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall std_basis Cvoid (Ptr{Cvoid},) prob.p
end

function adv_basis(prob::Prob, flags::Integer)
    @check! _prob(prob)
    @check _adv_basis_flags(flags)
    @glpk_ccall adv_basis Cvoid (Ptr{Cvoid}, Cint) prob.p flags
end
adv_basis(prob::Prob) = adv_basis(prob, 0)

function cpx_basis(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall cpx_basis Cvoid (Ptr{Cvoid},) prob.p
end

function simplex(prob::Prob, param::SimplexParam)
    @check! _prob(prob)
    return @glpk_ccall simplex Cint (Ptr{Cvoid}, Ref{SimplexParam}) prob.p param
end
function simplex(prob::Prob, param::Nothing)
    @check! _prob(prob)
    return @glpk_ccall simplex Cint (Ptr{Cvoid}, Ptr{Cvoid}) prob.p C_NULL
end
simplex(prob::Prob) = simplex(prob, nothing)

function exact(prob::Prob, param::SimplexParam)
    @check! _prob(prob)
    return @glpk_ccall exact Cint (Ptr{Cvoid}, Ref{SimplexParam}) prob.p param
end
function exact(prob::Prob, param::Nothing)
    @check! _prob(prob)
    return @glpk_ccall exact Cint (Ptr{Cvoid}, Ptr{Cvoid}) prob.p C_NULL
end
exact(prob::Prob) = exact(prob, nothing)

function init_smcp(param::SimplexParam)
    @glpk_ccall init_smcp Cint (Ref{SimplexParam},) param
end

function get_status(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall get_status Cint (Ptr{Cvoid},) prob.p
end

function get_prim_stat(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall get_prim_stat Cint (Ptr{Cvoid},) prob.p
end

function get_dual_stat(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall get_dual_stat Cint (Ptr{Cvoid},) prob.p
end

function get_obj_val(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall get_obj_val Cdouble (Ptr{Cvoid},) prob.p
end

function get_row_stat(prob::Prob, row::Integer)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    @glpk_ccall get_row_stat Cint (Ptr{Cvoid}, Cint) prob.p row
end

function get_row_prim(prob::Prob, row::Integer)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    @glpk_ccall get_row_prim Cdouble (Ptr{Cvoid}, Cint) prob.p row
end

function get_row_dual(prob::Prob, row::Integer)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    @glpk_ccall get_row_dual Cdouble (Ptr{Cvoid}, Cint) prob.p row
end

function get_col_stat(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @glpk_ccall get_col_stat Cint (Ptr{Cvoid}, Cint) prob.p col
end

function get_col_prim(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @glpk_ccall get_col_prim Cdouble (Ptr{Cvoid}, Cint) prob.p col
end

function get_col_dual(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @glpk_ccall get_col_dual Cdouble (Ptr{Cvoid}, Cint) prob.p col
end

function get_unbnd_ray(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall get_unbnd_ray Cint (Ptr{Cvoid},) prob.p
end

function interior(prob::Prob, param::InteriorParam)
    @check! _prob(prob)
    return @glpk_ccall interior Cint (Ptr{Cvoid}, Ref{InteriorParam}) prob.p param
end
function interior(prob::Prob, param::Nothing)
    @check! _prob(prob)
    return @glpk_ccall interior Cint (Ptr{Cvoid}, Ptr{Cvoid}) prob.p C_NULL
end
interior(prob::Prob) = interior(prob, nothing)

function init_iptcp(param::InteriorParam)
    @glpk_ccall init_iptcp Cint (Ref{InteriorParam},) param
end

function ipt_status(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall ipt_status Cint (Ptr{Cvoid},) prob.p
end

function ipt_obj_val(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall ipt_obj_val Cdouble (Ptr{Cvoid},) prob.p
end

function ipt_row_prim(prob::Prob, row::Integer)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    @glpk_ccall ipt_row_prim Cdouble (Ptr{Cvoid}, Cint) prob.p row
end

function ipt_row_dual(prob::Prob, row::Integer)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    @glpk_ccall ipt_row_dual Cdouble (Ptr{Cvoid}, Cint) prob.p row
end

function ipt_col_prim(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @glpk_ccall ipt_col_prim Cdouble (Ptr{Cvoid}, Cint) prob.p col
end

function ipt_col_dual(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @glpk_ccall ipt_col_dual Cdouble (Ptr{Cvoid}, Cint) prob.p col
end

function set_col_kind(prob::Prob, col::Integer, kind::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @check _kind_is_valid(kind)
    @glpk_ccall set_col_kind Cvoid (Ptr{Cvoid}, Cint, Cint) prob.p col kind
end

function get_col_kind(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @glpk_ccall get_col_kind Cint (Ptr{Cvoid}, Cint) prob.p col
end

function get_num_int(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall get_num_int Cint (Ptr{Cvoid},) prob.p
end

function get_num_bin(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall get_num_bin Cint (Ptr{Cvoid},) prob.p
end

function intopt(prob::Prob, param::IntoptParam)
    @check! _prob(prob)
    return @glpk_ccall intopt Cint (Ptr{Cvoid}, Ref{IntoptParam}) prob.p param
end
function intopt(prob::Prob, param::Nothing)
    @check! _prob(prob)
    return @glpk_ccall intopt Cint (Ptr{Cvoid}, Ptr{Cvoid}) prob.p C_NULL
end
intopt(prob::Prob) = intopt(prob, nothing)

function init_iocp(param::IntoptParam)
    @glpk_ccall init_iocp Cint (Ref{IntoptParam},) param
end

function mip_status(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall mip_status Cint (Ptr{Cvoid},) prob.p
end

function mip_obj_val(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall mip_obj_val Cdouble (Ptr{Cvoid},) prob.p
end

function mip_row_val(prob::Prob, row::Integer)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    @glpk_ccall mip_row_val Cdouble (Ptr{Cvoid}, Cint) prob.p row
end

function mip_col_val(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @glpk_ccall mip_col_val Cdouble (Ptr{Cvoid}, Cint) prob.p col
end

function check_kkt(prob::Prob, sol, cond, ae_max, ae_ind, re_max, re_ind)
    error("unsupported. Use GLPK.check_kkt(prob, sol, cond) instead.")
end

function check_kkt(prob::Prob, sol::Integer, cond::Integer)
    @check! _prob(prob)
    @check _sol_param(sol)
    @check _kkt_cond_param(cond, sol)

    ae_max = Ref{Cdouble}()
    ae_ind = Ref{Cint}()
    re_max = Ref{Cdouble}()
    re_ind = Ref{Cint}()

    @glpk_ccall check_kkt Cvoid (Ptr{Cvoid}, Cint, Cint, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}) prob.p sol cond ae_max ae_ind re_max re_ind

    return ae_max[], ae_ind[], re_max[], re_ind[]
end

function read_mps(prob::Prob, format::Integer, param, filename::AbstractString)
    @check! _prob(prob)
    @check _mps_format(format)
    if param == nothing
        param = C_NULL
    else
        @check _mps_param(param)
    end

    @check _file_is_readable(filename)
    ret = @glpk_ccall read_mps Cint (Ptr{Cvoid}, Cint, Ptr{Cvoid}, Ptr{Cchar}) prob.p format param string(filename)
    @check! _succeeded(ret, "read_mps")
    return ret
end

read_mps(prob::Prob, format::Integer, filename::AbstractString) =
    read_mps(prob, format, C_NULL, filename)

function write_mps(prob::Prob, format::Integer, param, filename::AbstractString)
    @check! _prob(prob)
    @check _mps_format(format)
    if param == nothing
        param = C_NULL
    else
        @check _mps_param(param)
    end
    @check _file_is_writable(filename)
    ret = @glpk_ccall write_mps Cint (Ptr{Cvoid}, Cint, Ptr{Cvoid}, Ptr{Cchar}) prob.p format param string(filename)
    @check! _succeeded(ret, "write_mps")
    return ret
end

write_mps(prob::Prob, format::Integer, filename::AbstractString) =
    write_mps(prob, format, C_NULL, filename)

function read_lp(prob::Prob, param, filename::AbstractString)
    @check! _prob(prob)
    @check _lp_param(param)
    @check _file_is_readable(filename)
    ret = @glpk_ccall read_lp Cint (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cchar}) prob.p param string(filename)
    @check! _succeeded(ret, "read_lp")
    return ret
end

read_lp(prob::Prob, filename::AbstractString) =
    read_lp(prob, C_NULL, filename)

function write_lp(prob::Prob, param, filename::AbstractString)
    @check! _prob(prob)
    @check _lp_param(param)
    @check _file_is_writable(filename)
    ret = @glpk_ccall write_lp Cint (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cchar}) prob.p param string(filename)
    @check! _succeeded(ret, "write_lp")
    return ret
end

write_lp(prob::Prob, filename::AbstractString) =
    write_lp(prob, C_NULL, filename)

function read_prob(prob::Prob, flags::Integer, filename::AbstractString)
    @check! _prob(prob)
    @check _read_prob_flags(flags)
    @check _file_is_readable(filename)
    ret = @glpk_ccall read_prob Cint (Ptr{Cvoid}, Cint, Ptr{Cchar}) prob.p flags string(filename)
    @check! _succeeded(ret, "read_prob")
    return ret
end

read_prob(prob::Prob, filename::AbstractString) =
    read_prob(prob, 0, filename)

function write_prob(prob::Prob, flags::Integer, filename::AbstractString)
    @check! _prob(prob)
    @check _write_prob_flags(flags)
    @check _file_is_writable(filename)
    ret = @glpk_ccall write_prob Cint (Ptr{Cvoid}, Cint, Ptr{Cchar}) prob.p flags string(filename)
    @check! _succeeded(ret, "write_prob")
    return ret
end

write_prob(prob::Prob, filename::AbstractString) =
    write_prob(prob, 0, filename)

function mpl_read_model(tran::MathProgWorkspace, filename::AbstractString, skip::Integer)
    @check! _mpl_workspace(tran)
    @check _file_is_readable(filename)
    ret = @glpk_ccall mpl_read_model Cint (Ptr{Cvoid}, Ptr{Cchar}, Cint) tran.p string(filename) skip
    @check! _succeeded(ret, "mpl_read_model")
    return ret
end

function mpl_read_data(tran::MathProgWorkspace, filename::AbstractString)
    @check! _mpl_workspace(tran)
    @check _file_is_readable(filename)
    ret = @glpk_ccall mpl_read_data Cint (Ptr{Cvoid}, Ptr{Cchar}) tran.p string(filename)
    @check! _succeeded(ret, "mpl_read_data")
    return ret
end

function mpl_generate(tran::MathProgWorkspace, filename::Union{AbstractString,Nothing})
    @check! _mpl_workspace(tran)
    if filename == nothing
        cfilename = C_NULL
    else
        @check _file_is_writable(filename)
        cfilename = string(filename)
    end
    ret = @glpk_ccall mpl_generate Cint (Ptr{Cvoid}, Ptr{Cchar}) tran.p cfilename
    @check! _succeeded(ret, "mpl_generate")
    return ret

end
mpl_generate(tran::MathProgWorkspace) = mpl_generate(tran, nothing)

function mpl_build_prob(tran::MathProgWorkspace, prob::Prob)
    @check! _mpl_workspace(tran)
    @check! _prob(prob)
    @glpk_ccall mpl_build_prob Cvoid (Ptr{Cvoid}, Ptr{Cvoid}) tran.p prob.p
end

function mpl_postsolve(tran::MathProgWorkspace, prob::Prob, sol::Integer)
    @check! _mpl_workspace(tran)
    @check! _prob(prob)
    @check _sol_param(sol)
    ret = @glpk_ccall mpl_postsolve Cint (Ptr{Cvoid}, Ptr{Cvoid}, Cint) tran.p prob.p sol
    @check! _succeeded(ret, "mpl_postsolve")
    return ret
end

function print_sol(prob::Prob, filename::AbstractString)
    @check! _prob(prob)
    @check _file_is_writable(filename)
    ret = @glpk_ccall print_sol Cint (Ptr{Cvoid}, Ptr{Cchar}) prob.p string(filename)
    @check! _succeeded(ret, "print_sol")
    return ret
end

function read_sol(prob::Prob, filename::AbstractString)
    @check! _prob(prob)
    @check _file_is_readable(filename)
    ret = @glpk_ccall read_sol Cint (Ptr{Cvoid}, Ptr{Cchar}) prob.p string(filename)
    @check! _succeeded(ret, "read_sol")
    return ret
end

function write_sol(prob::Prob, filename::AbstractString)
    @check! _prob(prob)
    @check _file_is_writable(filename)
    ret = @glpk_ccall write_sol Cint (Ptr{Cvoid}, Ptr{Cchar}) prob.p string(filename)
    @check! _succeeded(ret, "write_sol")
    return ret
end

function print_ipt(prob::Prob, filename::AbstractString)
    @check! _prob(prob)
    @check _file_is_writable(filename)
    ret = @glpk_ccall print_ipt Cint (Ptr{Cvoid}, Ptr{Cchar}) prob.p string(filename)
    @check! _succeeded(ret, "print_ipt")
    return ret
end

function read_ipt(prob::Prob, filename::AbstractString)
    @check! _prob(prob)
    @check _file_is_readable(filename)
    ret = @glpk_ccall read_ipt Cint (Ptr{Cvoid}, Ptr{Cchar}) prob.p string(filename)
    @check! _succeeded(ret, "read_ipt")
    return ret
end

function write_ipt(prob::Prob, filename::AbstractString)
    @check! _prob(prob)
    @check _file_is_writable(filename)
    ret = @glpk_ccall write_ipt Cint (Ptr{Cvoid}, Ptr{Cchar}) prob.p string(filename)
    @check! _succeeded(ret, "write_ipt")
    return ret
end

function print_mip(prob::Prob, filename::AbstractString)
    @check! _prob(prob)
    @check _file_is_writable(filename)
    ret = @glpk_ccall print_mip Cint (Ptr{Cvoid}, Ptr{Cchar}) prob.p string(filename)
    @check! _succeeded(ret, "print_mip")
    return ret
end

function read_mip(prob::Prob, filename::AbstractString)
    @check! _prob(prob)
    @check _file_is_readable(filename)
    ret = @glpk_ccall read_mip Cint (Ptr{Cvoid}, Ptr{Cchar}) prob.p string(filename)
    @check! _succeeded(ret, "read_mip")
    return ret
end

function write_mip(prob::Prob, filename::AbstractString)
    @check! _prob(prob)
    @check _file_is_writable(filename)
    ret = @glpk_ccall write_mip Cint (Ptr{Cvoid}, Ptr{Cchar}) prob.p string(filename)
    @check! _succeeded(ret, "write_mip")
    return ret
end

function print_ranges(prob::Prob, len::Integer, list::VecOrNothing, flags::Integer, filename::AbstractString)
    @check! _prob(prob)
    @check! _vectors_size(len, list)
    @check _status_is_optimal(prob)
    @check _bf_exists(prob)
    @check _print_ranges_flags(flags)
    @check _file_is_writable(filename)

    if len > 0
        list32 = cint_vec(list)
        @check _list_ids(prob, len, list32)

        off32 = sizeof(Cint)
        list32p = pointer(list32) - off32
    else
        list32p = C_NULL
    end

    @glpk_ccall print_ranges Cint (Ptr{Cvoid}, Cint, Ptr{Cint}, Cint, Ptr{Cchar}) prob.p len list32p flags string(filename)
end

print_ranges(prob::Prob, list::VecOrNothing, flags::Integer, filename::AbstractString) =
    print_ranges(prob, length(list), list, flags, filename)

print_ranges(prob::Prob, len::Integer, list::VecOrNothing, filename::AbstractString) =
    print_ranges(prob, len, list, 0, filename)

print_ranges(prob::Prob, list::VecOrNothing, filename::AbstractString) =
    print_ranges(prob, length(list), list, 0, filename)

print_ranges(prob::Prob, filename::AbstractString) =
    print_ranges(prob, 0, nothing, 0, filename)

function bf_exists(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall bf_exists Cint (Ptr{Cvoid},) prob.p
end

function factorize(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall factorize Cint (Ptr{Cvoid},) prob.p
end

function bf_updated(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall bf_updated Cint (Ptr{Cvoid},) prob.p
end

function get_bfcp(prob::Prob, param::BasisFactParam)
    @check! _prob(prob)
    @glpk_ccall get_bfcp Cvoid (Ptr{Cvoid}, Ref{BasisFactParam}) prob.p param
end

function set_bfcp(prob::Prob, param::BasisFactParam)
    @check! _prob(prob)
    return @glpk_ccall set_bfcp Cvoid (Ptr{Cvoid}, Ref{BasisFactParam}) prob.p param
end
function set_bfcp(prob::Prob, param::Nothing)
    @check! _prob(prob)
    return @glpk_ccall set_bfcp Cvoid (Ptr{Cvoid}, Ptr{Cvoid}) prob.p C_NULL
end
set_bfcp(prob::Prob) = set_bfcp(prob, nothing)

function get_bhead(prob::Prob, k::Integer)
    @check! _prob(prob)
    @check _bf_exists(prob)
    @check _row_is_valid(prob, k)
    @glpk_ccall get_bhead Cint (Ptr{Cvoid}, Cint) prob.p k
end

function get_row_bind(prob::Prob, row::Integer)
    @check! _prob(prob)
    @check _bf_exists(prob)
    @check _row_is_valid(prob, row)
    @glpk_ccall get_row_bind Cint (Ptr{Cvoid}, Cint) prob.p row
end

function get_col_bind(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _bf_exists(prob)
    @check _col_is_valid(prob, col)
    @glpk_ccall get_col_bind Cint (Ptr{Cvoid}, Cint) prob.p col
end

function ftran(prob::Prob, x::Vector{Cdouble})
    @check! _prob(prob)
    rows = @glpk_ccall get_num_rows Cint (Ptr{Cvoid},) prob.p
    @check! _vectors_size(rows, x)
    off64 = sizeof(Cdouble)
    x64p = pointer(x) - off64
    @glpk_ccall ftran Cvoid (Ptr{Cvoid}, Ptr{Cdouble}) prob.p x64p
end

function btran(prob::Prob, x::Vector{Cdouble})
    @check! _prob(prob)
    rows = @glpk_ccall get_num_rows Cint (Ptr{Cvoid},) prob.p
    @check! _vectors_size(rows, x)
    off64 = sizeof(Cdouble)
    x64p = pointer(x) - off64
    @glpk_ccall btran Cvoid (Ptr{Cvoid}, Ptr{Cdouble}) prob.p x64p
end

function warm_up(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall warm_up Cint (Ptr{Cvoid},) prob.p
end

function eval_tab_row(prob::Prob, k::Integer, ind::Vector{Cint}, val::Vector{Cdouble})
    @check! _prob(prob)
    @check _bf_exists(prob)
    rows = @glpk_ccall get_num_rows Cint (Ptr{Cvoid},) prob.p
    cols = @glpk_ccall get_num_cols Cint (Ptr{Cvoid},) prob.p
    k_max = rows + cols
    @check _rowcol_is_valid(k, k_max)
    @check _var_is_basic(prob, k)

    resize!(ind, k_max)
    resize!(val, k_max)

    off32 = sizeof(Cint)
    off64 = sizeof(Cdouble)
    ind32p = pointer(ind) - off32
    val64p = pointer(val) - off64

    len = @glpk_ccall eval_tab_row Cint (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}) prob.p k ind32p val64p

    splice!(ind, len+1:length(ind))
    splice!(val, len+1:length(val))

    return len
end

function eval_tab_row(prob::Prob, k::Integer)
    @check! _prob(prob)
    @check _bf_exists(prob)
    rows = @glpk_ccall get_num_rows Cint (Ptr{Cvoid},) prob.p
    cols = @glpk_ccall get_num_cols Cint (Ptr{Cvoid},) prob.p
    k_max = rows + cols
    @check _rowcol_is_valid(k, k_max)
    @check _var_is_basic(prob, k)

    ind = Array{Cint}(undef, k_max)
    val = Array{Cdouble}(undef, k_max)

    off32 = sizeof(Cint)
    off64 = sizeof(Cdouble)
    ind32p = pointer(ind) - off32
    val64p = pointer(val) - off64

    len = @glpk_ccall eval_tab_row Cint (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}) prob.p k ind32p val64p

    splice!(ind, len+1:length(ind))
    splice!(val, len+1:length(val))

    return ind, val
end

function eval_tab_col(prob::Prob, k::Integer, ind::Vector{Cint}, val::Vector{Cdouble})
    @check! _prob(prob)
    @check _bf_exists(prob)
    rows = @glpk_ccall get_num_rows Cint (Ptr{Cvoid},) prob.p
    cols = @glpk_ccall get_num_cols Cint (Ptr{Cvoid},) prob.p
    k_max = rows + cols
    @check _rowcol_is_valid(k, k_max)
    @check _var_is_non_basic(prob, k)

    resize!(ind, k_max)
    resize!(val, k_max)

    off32 = sizeof(Cint)
    off64 = sizeof(Cdouble)
    ind32p = pointer(ind) - off32
    val64p = pointer(val) - off64

    len = @glpk_ccall eval_tab_col Cint (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}) prob.p k ind32p val64p

    splice!(ind, len+1:length(ind))
    splice!(val, len+1:length(val))

    return len
end

function eval_tab_col(prob::Prob, k::Integer)
    @check! _prob(prob)
    @check _bf_exists(prob)
    rows = @glpk_ccall get_num_rows Cint (Ptr{Cvoid},) prob.p
    cols = @glpk_ccall get_num_cols Cint (Ptr{Cvoid},) prob.p
    k_max = rows + cols
    @check _rowcol_is_valid(k, k_max)
    @check _var_is_non_basic(prob, k)

    ind = Array{Cint}(undef, k_max)
    val = Array{Cdouble}(undef, k_max)

    off32 = sizeof(Cint)
    off64 = sizeof(Cdouble)
    ind32p = pointer(ind) - off32
    val64p = pointer(val) - off64

    len = @glpk_ccall eval_tab_col Cint (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}) prob.p k ind32p val64p

    splice!(ind, len+1:length(ind))
    splice!(val, len+1:length(val))

    return ind, val
end

function transform_row(prob::Prob, len::Integer, ind::Vector{Cint}, val::Vector{Cdouble})
    @check! _prob(prob)
    @check _bf_exists(prob)
    @check _col_is_valid(prob, len)
    @check! _vectors_size(len, ind, val)

    cols = @glpk_ccall get_num_cols Cint (Ptr{Cvoid},) prob.p

    resize!(ind, cols)
    resize!(val, cols)

    off32 = sizeof(Cint)
    off64 = sizeof(Cdouble)
    ind32p = pointer(ind) - off32
    val64p = pointer(val) - off64

    len1 = @glpk_ccall transform_row Cint (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}) prob.p len ind32p val64p

    splice!(ind, len1+1:length(ind))
    splice!(val, len1+1:length(val))

    return len1
end

function transform_row(prob::Prob, ind::Vector{Cint}, val::Vector{Cdouble})
    @check! _vectors_all_same_size(ind, val)
    transform_row(prob, length(ind), ind, val)
end

function transform_col(prob::Prob, len::Integer, ind::Vector{Cint}, val::Vector{Cdouble})
    @check! _prob(prob)
    @check _bf_exists(prob)
    @check _row_is_valid(prob, len)
    @check! _vectors_size(len, ind, val)

    rows = @glpk_ccall get_num_rows Cint (Ptr{Cvoid},) prob.p

    resize!(ind, rows)
    resize!(val, rows)

    off32 = sizeof(Cint)
    off64 = sizeof(Cdouble)
    ind32p = pointer(ind) - off32
    val64p = pointer(val) - off64

    len1 = @glpk_ccall transform_col Cint (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}) prob.p len ind32p val64p

    splice!(ind, len1+1:length(ind))
    splice!(val, len1+1:length(val))

    return len1
end

function transform_col(prob::Prob, ind::Vector{Cint}, val::Vector{Cdouble})
    @check! _vectors_all_same_size(ind, val)
    transform_col(prob, length(ind), ind, val)
end

function prim_rtest(prob::Prob, len::Integer, ind::VecOrNothing, val::VecOrNothing, dir::Integer, eps::Real)
    @check! _prob(prob)
    @check _bf_exists(prob)
    @check _is_prim_feasible(prob)
    @check _row_is_valid(prob, len)
    @check! _vectors_size(len, ind, val)
    @check _dir_is_valid(dir)
    @check _eps_is_valid(eps)

    for i = 1:len
        @check _var_is_basic(prob, ind[i])
    end

    ind32 = cint_vec(ind)
    val64 = cdouble_vec(val)
    off32 = sizeof(Cint)
    off64 = sizeof(Cdouble)
    ind32p = pointer(ind32) - off32
    val64p = pointer(val64) - off64

    piv = @glpk_ccall prim_rtest Cint (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}, Cint, Cdouble) prob.p len ind32p val64p dir eps
    return piv
end

function prim_rtest(prob::Prob, ind::VecOrNothing, val::VecOrNothing, dir::Integer, eps::Real)
    @check! _vectors_all_same_size(ind, val)
    prim_rtest(prob, length(ind), ind, val, dir, eps)
end

function dual_rtest(prob::Prob, len::Integer, ind::VecOrNothing, val::VecOrNothing, dir::Integer, eps::Real)
    @check! _prob(prob)
    @check _bf_exists(prob)
    @check _is_dual_feasible(prob)
    @check _col_is_valid(prob, len)
    @check! _vectors_size(len, ind, val)
    @check _dir_is_valid(dir)
    @check _eps_is_valid(eps)

    for i = 1:len
        @check _var_is_non_basic(prob, ind[i])
    end

    ind32 = cint_vec(ind)
    val64 = cdouble_vec(val)
    off32 = sizeof(Cint)
    off64 = sizeof(Cdouble)
    ind32p = pointer(ind32) - off32
    val64p = pointer(val64) - off64

    piv = @glpk_ccall dual_rtest Cint (Ptr{Cvoid}, Cint, Ptr{Cint}, Ptr{Cdouble}, Cint, Cdouble) prob.p len ind32p val64p dir eps
    return piv
end

function dual_rtest(prob::Prob, ind::VecOrNothing, val::VecOrNothing, dir::Integer, eps::Real)
    @check! _vectors_all_same_size(ind, val)
    dual_rtest(prob, length(ind), ind, val, dir, eps)
end

analyze_bound(prob::Prob, k, limit1, var1, limit2, var2) =
    error("unsupported. Use GLPK.analyze_bound(prob, k) instead.")

function analyze_bound(prob::Prob, k::Integer)
    @check! _prob(prob)
    @check _bf_exists(prob)
    @check _rowcol_is_valid(prob, k)
    @check _var_is_non_basic(prob, k)

    limit1 = Ref{Cdouble}()
    var1 = Ref{Cint}()
    limit2 = Ref{Cdouble}()
    var2 = Ref{Cint}()

    @glpk_ccall analyze_bound Cvoid (Ptr{Cvoid}, Cint, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}) prob.p k limit1 var1 limit2 var2

    return limit1[], var1[], limit2[], var2[]
end

analyze_coef(prob::Prob, k, coef1, var1, value1, coef2, var2, value2) =
    error("unsupported. Use GLPK.analyze_coef(prob, k) instead.")

function analyze_coef(prob::Prob, k::Integer)
    @check! _prob(prob)
    @check _bf_exists(prob)
    @check _rowcol_is_valid(prob, k)
    @check _var_is_basic(prob, k)

    coef1 = Ref{Cdouble}()
    var1 = Ref{Cint}()
    value1 = Ref{Cdouble}()
    coef2 = Ref{Cdouble}()
    var2 = Ref{Cint}()
    value2 = Ref{Cdouble}()

    @glpk_ccall analyze_coef Cvoid (Ptr{Cvoid}, Cint, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}) prob.p k coef1 var1 value1 coef2 var2 value2

    return coef1[], var1[], value1[], coef2[], var2[], value2[]
end

function ios_reason(tree::Ptr{Cvoid})
    @check! _tree(tree)
    @glpk_ccall ios_reason Cint (Ptr{Cvoid},) tree
end

function ios_get_prob(tree::Ptr{Cvoid})
    @check! _tree(tree)
    p = @glpk_ccall ios_get_prob Ptr{Cvoid} (Ptr{Cvoid},) tree
    return Prob(p)
end

function ios_row_attr(tree::Ptr{Cvoid}, row::Integer, attr::Attr)
    @check! _tree(tree)
    prob = ios_get_prob(tree)
    @check _row_is_valid(prob, row)
    @glpk_ccall ios_row_attr Cvoid (Ptr{Cvoid}, Cint, Ref{Attr}) tree row attr
end

function ios_row_attr(tree::Ptr{Cvoid}, row::Integer)
    @check! _tree(tree)
    prob = ios_get_prob(tree)
    @check _row_is_valid(prob, row)
    attr = Attr()
    @glpk_ccall ios_row_attr Cvoid (Ptr{Cvoid}, Cint, Ref{Attr}) tree row attr
    return attr
end

function ios_mip_gap(tree::Ptr{Cvoid})
    @check! _tree(tree)
    @glpk_ccall ios_mip_gap Cdouble (Ptr{Cvoid},) tree
end

function ios_node_data(tree::Ptr{Cvoid}, p::Integer)
    @check! _tree(tree)
    @check _ios_node_is_valid(tree, p)
    @glpk_ccall ios_node_data Ptr{Cvoid} (Ptr{Cvoid}, Cint) tree p
end

function ios_select_node(tree::Ptr{Cvoid}, p::Integer)
    @check! _tree(tree)
    @check _reason(tree, [GLPK.ISELECT])
    @check _ios_node_is_active(tree, p)
    @glpk_ccall ios_select_node Cvoid (Ptr{Cvoid}, Cint) tree p
end

function ios_heur_sol(tree::Ptr{Cvoid}, x::VecOrNothing)
    @check! _tree(tree)
    prob = ios_get_prob(tree)
    cols = @glpk_ccall get_num_cols Cint (Ptr{Cvoid},) prob.p
    @check! _vectors_size(cols, x)
    x64 = cdouble_vec(x)
    off64 = sizeof(Cdouble)
    x64p = pointer(x64) - off64

    @glpk_ccall ios_heur_sol Cint (Ptr{Cvoid}, Ptr{Cdouble}) tree x64p
end

function ios_can_branch(tree::Ptr{Cvoid}, col::Integer)
    @check! _tree(tree)
    prob = ios_get_prob(tree)
    @check _col_is_valid(prob, col)
    @glpk_ccall ios_can_branch Cint (Ptr{Cvoid}, Cint) tree col
end

function ios_branch_upon(tree::Ptr{Cvoid}, col::Integer, sel::Integer)
    @check! _tree(tree)
    prob = ios_get_prob(tree)
    @check _col_is_valid(prob, col)
    @check _can_branch(tree, col)
    @check _sel_is_valid(sel)
    @glpk_ccall ios_branch_upon Cvoid (Ptr{Cvoid}, Cint, Cint) tree col sel
end

function ios_terminate(tree::Ptr{Cvoid})
    @check! _tree(tree)
    @glpk_ccall ios_terminate Cvoid (Ptr{Cvoid},) tree
end

function ios_tree_size(tree::Ptr{Cvoid})
    @check! _tree(tree)
    a_cnt = Ref{Cint}()
    n_cnt = Ref{Cint}()
    t_cnt = Ref{Cint}()
    @glpk_ccall ios_tree_size Cvoid (Ptr{Cvoid}, Ptr{Cint}, Ptr{Cint}, Ptr{Cint}) tree a_cnt n_cnt t_cnt

    return a_cnt[], n_cnt[], t_cnt[]
end

ios_tree_size(tree::Ptr{Cvoid}, a_cnt, n_cnt, t_cnt) =
    error("unsupported. Use GLPK.ios_tree_size(tree) instead.")

function ios_curr_node(tree::Ptr{Cvoid})
    @check! _tree(tree)
    @glpk_ccall ios_curr_node Cint (Ptr{Cvoid},) tree
end

function ios_next_node(tree::Ptr{Cvoid}, p::Integer)
    @check! _tree(tree)
    p != 0 && @check _ios_node_is_active(tree, p)
    @glpk_ccall ios_next_node Cint (Ptr{Cvoid}, Cint) tree p
end

function ios_prev_node(tree::Ptr{Cvoid}, p::Integer)
    @check! _tree(tree)
    p != 0 && @check _ios_node_is_active(tree, p)
    @glpk_ccall ios_prev_node Cint (Ptr{Cvoid}, Cint) tree p
end

function ios_up_node(tree::Ptr{Cvoid}, p::Integer)
    @check! _tree(tree)
    @check _ios_node_is_valid(tree, p)
    @glpk_ccall ios_up_node Cint (Ptr{Cvoid}, Cint) tree p
end

function ios_node_level(tree::Ptr{Cvoid}, p::Integer)
    @check! _tree(tree)
    @check _ios_node_is_valid(tree, p)
    @glpk_ccall ios_node_level Cint (Ptr{Cvoid}, Cint) tree p
end

function ios_node_bound(tree::Ptr{Cvoid}, p::Integer)
    @check! _tree(tree)
    @check _ios_node_is_valid(tree, p)
    @glpk_ccall ios_node_bound Cdouble (Ptr{Cvoid}, Cint) tree p
end

function ios_best_node(tree::Ptr{Cvoid})
    @check! _tree(tree)
    @glpk_ccall ios_best_node Cint (Ptr{Cvoid},) tree
end

function ios_pool_size(tree::Ptr{Cvoid})
    @check! _tree(tree)
    @check _reason(tree, [GLPK.ICUTGEN])
    @glpk_ccall ios_pool_size Cint (Ptr{Cvoid},) tree
end

function ios_add_row(tree::Ptr{Cvoid}, name::Union{AbstractString,Nothing}, klass::Integer, flags::Integer,
                     len::Integer, ind::VecOrNothing, val::VecOrNothing, constr_type::Integer, rhs::Real)

    @check! _tree(tree)
    @check _reason(tree, [GLPK.ICUTGEN])
    name == nothing && (name = "")
    @check _string_length(name, 255)
    @check _klass_is_valid(klass)
    @check _ios_add_row_flags(flags)
    @check! _vectors_size(len, ind, val)
    if len > 0
        ind32 = cint_vec(ind)
        val64 = cdouble_vec(val)
        off32 = sizeof(Cint)
        off64 = sizeof(Cdouble)
        ind32p = pointer(ind32) - off32
        val64p = pointer(val64) - off64
        prob = ios_get_prob(tree)
        @check! _cols_ids_size(prob, 0, len, ind32)
        @check _cols_ids_content(prob, len, ind32)
    else
        ind32p = C_NULL
        val64p = C_NULL
    end
    @check _constr_type_is_valid(constr_type)

    @glpk_ccall ios_add_row Cint (Ptr{Cvoid}, Ptr{Cchar}, Cint, Cint, Cint, Ptr{Cint}, Ptr{Cdouble}, Cint, Cdouble) tree string(name) klass flags len ind32p val64p constr_type rhs
end

ios_add_row(tree::Ptr{Cvoid}, klass::Integer, flags::Integer, len::Integer, ind::VecOrNothing, val::VecOrNothing,
            constr_type::Integer, rhs::Real) =
    ios_add_row(tree, nothing, klass, flags, len, ind, val, constr_type, rhs)

function ios_add_row(tree::Ptr{Cvoid}, name::Union{AbstractString,Nothing}, klass::Integer, flags::Integer,
                     ind::VecOrNothing, val::VecOrNothing, constr_type::Integer, rhs::Real)
    @check! _vectors_all_same_size(ind, val)
    ios_add_row(tree, name, klass, flags, length(ind), ind, val, constr_type, rhs)
end

ios_add_row(tree::Ptr{Cvoid}, name::Union{AbstractString,Nothing}, klass::Integer, ind::VecOrNothing, val::VecOrNothing,
            constr_type::Integer, rhs::Real) =
    ios_add_row(tree, name, klass, 0, length(ind), ind, val, constr_type, rhs)

ios_add_row(tree::Ptr{Cvoid}, klass::Integer, ind::VecOrNothing, val::VecOrNothing, constr_type::Integer, rhs::Real) =
    ios_add_row(tree, nothing, klass, 0, length(ind), ind, val, constr_type, rhs)

function ios_del_row(tree::Ptr{Cvoid}, row::Integer)
    @check! _tree(tree)
    @check _reason(tree, [GLPK.ICUTGEN])
    @check _ios_row_is_valid(tree, row)
    @glpk_ccall ios_del_row Cvoid (Ptr{Cvoid}, Cint) tree row
end

function ios_clear_pool(tree::Ptr{Cvoid})
    @check! _tree(tree)
    @check _reason(tree, [GLPK.ICUTGEN])
    @glpk_ccall ios_clear_pool Cvoid (Ptr{Cvoid},) tree
end

function init_env()
    ret = @glpk_ccall init_env Cint ()
    @check! _init_env_succeeded(ret)
    return ret
end

function free_env()
    ret = @glpk_ccall free_env Cint ()
    _del_all_objs()
    return ret
end

function term_out(flag::Integer)
    @check _term_out_flag(flag)
    @glpk_ccall term_out Cint (Cint,) flag
end

function open_tee(filename::AbstractString)
    ret = @glpk_ccall open_tee Cint (Ptr{Cchar},) string(filename)
    @check! _open_tee_succeeded(ret)
    return ret
end

function close_tee()
    @glpk_ccall close_tee Cint ()
end

function malloc(size::Integer)
    @check _alloc_size(size)
    @glpk_ccall alloc Ptr{Cvoid} (Cint, Cint) 1 size
end

function calloc(n::Integer, size::Integer)
    @check _alloc_size(n)
    @check _alloc_size(size)
    @glpk_ccall alloc Ptr{Cvoid} (Cint, Cint) n size
end

function free(ptr::Ptr)
    @check _pointer_is_valid(ptr)
    @glpk_ccall free Cvoid (Ptr{Cvoid},) ptr
end

mem_usage(count, cpeak, total, tpeak) =
    error("unsupported. Use GLPK.mem_usage() instead.")

function mem_usage()
    data32 = Array{Cint}(undef, 2)
    data32_p = pointer(data32)
    off32 = sizeof(Cint)
    count_p = data32_p
    cpeak_p = data32_p + off32

    data64 = Array{Clong}(undef, 2)
    data64_p = pointer(data64)
    off64 = sizeof(Clong)

    total_p = data64_p
    tpeak_p = data64_p + off64

    @glpk_ccall mem_usage Cvoid (Ptr{Cint}, Ptr{Cint}, Ptr{Clong}, Ptr{Clong}) count_p cpeak_p total_p tpeak_p

    count = data32[1]
    cpeak = data32[2]
    total = data64[1]
    tpeak = data64[2]

    return count, cpeak, total, tpeak
end

function mem_limit(limit::Integer)
    @glpk_ccall mem_limit Cvoid (Cint,) limit
end

function read_cnfsat(prob::Prob, filename::AbstractString)
    @check! _prob(prob)
    @check _file_is_readable(filename)
    ret = @glpk_ccall read_cnfsat Cint (Ptr{Cvoid}, Ptr{Cchar}) prob.p string(filename)
    @check! _succeeded(ret, "read_cnfsat")
    return ret
end

function check_cnfsat(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall check_cnfsat Cint (Ptr{Cvoid},) prob.p
end

function write_cnfsat(prob::Prob, filename::AbstractString)
    @check! _prob(prob)
    @check _file_is_writable(filename)
    ret = @glpk_ccall write_cnfsat Cint (Ptr{Cvoid}, Ptr{Cchar}) prob.p string(filename)
    @check! _succeeded(ret, "write_cnfsat")
    return ret
end

function minisat1(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall minisat1 Cint (Ptr{Cvoid},) prob.p
end

function intfeas1(prob::Prob, use_bound::Integer, obj_bound::Integer)
    @check! _prob(prob)
    # TODO : more checks:
    #   1) columns must be GLPK.BV od GLPK.FX
    #   2) constraints and objj coeffs must be integer
    @glpk_ccall intfeas1 Cint (Ptr{Cvoid}, Cint, Cint) prob.p use_bound obj_bound
end


# FUNCTIONS NOT WRAPPED:
#
# glp_printf      (uses varargs)
# glp_vprintf     (uses varargs)
# glp_term_hook
# glp_error       (implemented as macro)
# glp_assert      (implemented as macro)
# glp_error_hook  (reserved for internal use by the wrapper)
#

#}}}

include("MOI_wrapper.jl")

end # module
