###
### GLPK API Wrapper
###

module GLPK

## Exports
#{{{
export
    # Types
    Param,
    GLPKError,
    Prob,
    SimplexParam,
    InteriorParam,
    IntoptParam,
    BasisFactParam,
    Data,
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
    time,
    difftime,
    sdf_open_file,
    sdf_read_int,
    sdf_read_num,
    sdf_read_item,
    sdf_read_text,
    sdf_line,
    sdf_close_file,
    read_cnfsat,
    check_cnfsat,
    write_cnfsat,
    minisat1,
    intfeas1
#}}}

import Base.pointer, Base.assign, Base.ref

## Shared library interface setup
#{{{
include("GLPK_constants.jl")

include(joinpath(Pkg.dir(),"GLPK","deps","ext.jl"))

# General recoverable exception: all GLPK functions
# throw this in case of recoverable errors
type GLPKError <: Exception
    msg::String
end

# Fatal exception: when this is thrown, all GLPK
# objects are no longer valid
type GLPKFatalError <: Exception
    msg::String
end

# Ugly way to deal with C's jmp_buf: just provides enough storage
# (should NEVER be accessed directly)
type JmpBuf
    _pad01::Float64
    _pad02::Float64
    _pad03::Float64
    _pad04::Float64
    _pad05::Float64
    _pad06::Float64
    _pad07::Float64
    _pad08::Float64
    _pad09::Float64
    _pad10::Float64
    _pad11::Float64
    _pad12::Float64
    _pad13::Float64
    _pad14::Float64
    _pad15::Float64
    _pad16::Float64
    _pad17::Float64
    _pad18::Float64
    _pad19::Float64
    _pad20::Float64
    _pad21::Float64
    _pad22::Float64
    _pad23::Float64
    _pad24::Float64
    _pad25::Float64
    _pad26::Float64
    _pad27::Float64
    _pad28::Float64
    _pad29::Float64
    _pad30::Float64
    _pad31::Float64
    _pad32::Float64
    JmpBuf() = new()
end

# Error hook, used to catch internal errors when calling
# GLPK functions
function _err_hook(info::Ptr{Void})
    ccall(:longjmp, Void, (Ptr{JmpBuf}, Int32), info, 1)
end

macro glpk_ccall(f, args...)
    quote
        let jb = JmpBuf()
            if ccall(:setjmp, Int, (Ptr{JmpBuf},), &jb) == 0
                # install a safeguard hook so that errors don't crash Julia
                ccall((:glp_error_hook, _jl_libGLPK), Void, (Ptr{Void}, Ptr{Void}), cfunction(_err_hook, Void, (Ptr{Void},)), &jb)
                # actual call to the function
                ret = ccall(($"glp_$(f)", _jl_libGLPK), $(args...))
                # the call went fine: uninstall the hook
                ccall((:glp_error_hook, _jl_libGLPK), Void, (Ptr{Void}, Ptr{Void}), C_NULL, C_NULL)
                # return the result of the call
                ret
            else
                # something went wrong
                ccall((:glp_free_env, _jl_libGLPK), Void, ())
                throw(GLPKFatalError("GLPK call failed. All GLPK objects you defined so far are now invalidated."))
            end
        end
    end
end

# We need to define GLPK.version as first thing
# in order to perform a sanity check
# (since we import structs from the header,
# we must ensure that the binary is the correct
# one)
function version()
    csp = @glpk_ccall version Ptr{Uint8} ()
    str = Array(Uint8, 100)
    strp = pointer(str)
    k = 0
    for i = 1 : 100
        ccall(:memcpy, Ptr{Void}, (Ptr{Void}, Ptr{Void}, Int), strp, csp, sizeof(Uint8))
        if str[i] == '\0'
            k = i
            break
        end
        strp += sizeof(Uint8)
        csp += sizeof(Uint8)
    end
    if k == 0
        throw(GLPKError("error reading version"))
    end
    vstr = ASCIIString(str[1:k - 1])
    return tuple(map(x->int32(parse_int(x)), split(vstr, '.'))...)
end

if version() != (MAJOR_VERSION, MINOR_VERSION)
    bv = version()
    hv = (MAJOR_VERSION, MINOR_VERSION)
    error("GLPK: mismatched versions: header=$(hv[1]).$(hv[2]) binary=$(bv[1]).$(bv[2])")
end
#}}}

## Preliminary definitions
#{{{

# General structure for the parameters types

abstract Param

function assign{T}(param::Param, val::T, field_name::String)
    param.(symbol(field_name)) = val
end

function ref(param::Param, field_name::String)
    param.(symbol(field_name))
end

# We define some types which allow to pass optional agruments
# to the function.
# In this framework, optional arguments can be passed either
# as an empty vector [] or as the 'nothing' constant

typealias VecOrNothing Union(Vector, Nothing)
function convert_vecornothing{T}(::Type{T}, a::VecOrNothing)
    if isequal(a, nothing) || isa(a, Array{None})
        return T[]
    elseif T <: Integer
        if !(eltype(a) <: Integer)
            throw(GLPKError("integer-valued array required, or [] or nothing"))
        end
    elseif T <: Real
        if !(eltype(a) <: Real)
            throw(GLPKError("real-valued array required, or [] or nothing"))
        end
    end
    convert(Array{T}, a)
end
vecornothing_length(a::VecOrNothing) = is(a, nothing) ? 0 : length(a)


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
# |  glp_data   |  GLPK.Data               |
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

type Prob
    p::Ptr{Void}
    function Prob()
        p = @glpk_ccall create_prob Ptr{Void} ()
        prob = new(p)
        finalizer(prob, delete_prob)
        return prob
    end
    function Prob(p::Ptr{Void})
        new(p)
    end
end

function delete_prob(prob::Prob)
    if prob.p == C_NULL
        return
    end
    @glpk_ccall delete_prob Void (Ptr{Void},) prob.p
    prob.p = C_NULL
    return
end

type SimplexParam <: Param
    msg_lev::Int32
    meth::Int32
    pricing::Int32
    r_test::Int32
    tol_bnd::Float64
    tol_dj::Float64
    tol_piv::Float64
    obj_ll::Float64
    obj_ul::Float64
    it_lim::Int32
    tm_lim::Int32
    out_frq::Int32
    out_dly::Int32
    presolve::Int32
    _reserved01::Float64
    _reserved02::Float64
    _reserved03::Float64
    _reserved04::Float64
    _reserved05::Float64
    _reserved06::Float64
    _reserved07::Float64
    _reserved08::Float64
    _reserved09::Float64
    _reserved10::Float64
    _reserved11::Float64
    _reserved12::Float64
    _reserved13::Float64
    _reserved14::Float64
    _reserved15::Float64
    _reserved16::Float64
    _reserved17::Float64
    _reserved18::Float64
    _reserved19::Float64
    _reserved20::Float64
    _reserved21::Float64
    _reserved22::Float64
    _reserved23::Float64
    _reserved24::Float64
    _reserved25::Float64
    _reserved26::Float64
    _reserved27::Float64
    _reserved28::Float64
    _reserved29::Float64
    _reserved30::Float64
    _reserved31::Float64
    _reserved32::Float64
    _reserved33::Float64
    _reserved34::Float64
    _reserved35::Float64
    _reserved36::Float64

    function SimplexParam()
        p = new()
        @glpk_ccall "init_smcp" Void (Ptr{SimplexParam},) &p
        return p
    end
end

type InteriorParam <: Param
    msg_lev::Int32
    ord_alg::Int32
    _reserved01::Float64
    _reserved02::Float64
    _reserved03::Float64
    _reserved04::Float64
    _reserved05::Float64
    _reserved06::Float64
    _reserved07::Float64
    _reserved08::Float64
    _reserved09::Float64
    _reserved10::Float64
    _reserved11::Float64
    _reserved12::Float64
    _reserved13::Float64
    _reserved14::Float64
    _reserved15::Float64
    _reserved16::Float64
    _reserved17::Float64
    _reserved18::Float64
    _reserved19::Float64
    _reserved20::Float64
    _reserved21::Float64
    _reserved22::Float64
    _reserved23::Float64
    _reserved24::Float64
    _reserved25::Float64
    _reserved26::Float64
    _reserved27::Float64
    _reserved28::Float64
    _reserved29::Float64
    _reserved30::Float64
    _reserved31::Float64
    _reserved32::Float64
    _reserved33::Float64
    _reserved34::Float64
    _reserved35::Float64
    _reserved36::Float64
    _reserved37::Float64
    _reserved38::Float64
    _reserved39::Float64
    _reserved40::Float64
    _reserved41::Float64
    _reserved42::Float64
    _reserved43::Float64
    _reserved44::Float64
    _reserved45::Float64
    _reserved46::Float64
    _reserved47::Float64
    _reserved48::Float64

    function InteriorParam()
        p = new()
        @glpk_ccall "init_iptcp" Void (Ptr{InteriorParam},) &p
        return p
    end
end

type IntoptParam <: Param
    msg_lev::Int32
    br_tech::Int32
    bt_tech::Int32
    tol_int::Float64
    tol_obj::Float64
    tm_lim::Int32
    out_frq::Int32
    out_dly::Int32
    cb_func::Ptr{Void}
    cb_info::Ptr{Void}
    cb_size::Int32
    pp_tech::Int32
    mip_gap::Float64
    mir_cuts::Int32
    gmi_cuts::Int32
    cov_cuts::Int32
    clq_cuts::Int32
    presolve::Int32
    binarize::Int32
    fp_heur::Int32
    alien::Int32
    _reserved01::Float64
    _reserved02::Float64
    _reserved03::Float64
    _reserved04::Float64
    _reserved05::Float64
    _reserved06::Float64
    _reserved07::Float64
    _reserved08::Float64
    _reserved09::Float64
    _reserved10::Float64
    _reserved11::Float64
    _reserved12::Float64
    _reserved13::Float64
    _reserved14::Float64
    _reserved15::Float64
    _reserved16::Float64
    _reserved17::Float64
    _reserved18::Float64
    _reserved19::Float64
    _reserved20::Float64
    _reserved21::Float64
    _reserved22::Float64
    _reserved23::Float64
    _reserved24::Float64
    _reserved25::Float64
    _reserved26::Float64
    _reserved27::Float64
    _reserved28::Float64
    _reserved29::Float64

    function IntoptParam()
        p = new()
        @glpk_ccall "init_iocp" Void (Ptr{IntoptParam},) &p
        return p
    end
end

type BasisFactParam <: Param
    msg_lev::Int32
    bftype::Int32 # NOTE: changed type->bftype
    lu_size::Int32
    piv_tol::Float64
    piv_lim::Int32
    suhl::Int32
    eps_tol::Float64
    max_gro::Float64
    nfs_max::Int32
    upd_tol::Float64
    nrs_max::Int32
    rs_size::Int32
    _reserved01::Float64
    _reserved02::Float64
    _reserved03::Float64
    _reserved04::Float64
    _reserved05::Float64
    _reserved06::Float64
    _reserved07::Float64
    _reserved08::Float64
    _reserved09::Float64
    _reserved10::Float64
    _reserved11::Float64
    _reserved12::Float64
    _reserved13::Float64
    _reserved14::Float64
    _reserved15::Float64
    _reserved16::Float64
    _reserved17::Float64
    _reserved18::Float64
    _reserved19::Float64
    _reserved20::Float64
    _reserved21::Float64
    _reserved22::Float64
    _reserved23::Float64
    _reserved24::Float64
    _reserved25::Float64
    _reserved26::Float64
    _reserved27::Float64
    _reserved28::Float64
    _reserved29::Float64
    _reserved30::Float64
    _reserved31::Float64
    _reserved32::Float64
    _reserved33::Float64
    _reserved34::Float64
    _reserved35::Float64
    _reserved36::Float64
    _reserved37::Float64
    _reserved38::Float64

    function BasisFactParam()
        return new(0, 0, 0, 0.0, 0, 0, 0.0, 0.0, 0, 0.0, 0, 0)
    end
end

type Data
    p::Ptr{Void}
end

pointer(data::Data) = data.p

type MathProgWorkspace
    p::Ptr{Void}
    function MathProgWorkspace()
        tran = @glpk_ccall mpl_alloc_wksp Ptr{Void} ()
        wksp = new(tran)
        finalizer(wksp, GLPK.mpl_free_wksp)
        return wksp
    end
end

function mpl_free_wksp(tran::MathProgWorkspace)
    if tran.p == C_NULL
        return
    end
    @glpk_ccall mpl_free_wksp Void (Ptr{Void},) tran.p
    tran.p = C_NULL
    return
end

type Attr
    level::Int32
    origin::Int32
    klass::Int32
    _reserved01::Float64
    _reserved02::Float64
    _reserved03::Float64
    _reserved04::Float64
    _reserved05::Float64
    _reserved06::Float64
    _reserved07::Float64
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
#  * function names tranlsate like this: glp_func -> GLPK.func
#  * constant names tranlsate like this: GLPK_CONST -> GLPK.CONST
#  * whenever the C library accepts NULL as argument,
#    the Julia one will accept the nothing constant.
#  * vectors do not need to have an extra element at the
#    beginning to accomodate to the 1-based GLP indexing.
#  * most functions will accept any kind of vectors as inputs,
#    provided they can be converted to be C-compatible
#    (i.e. to either Int32 (int) or Float64 (double) elements).
#  * the exceptions to the above are those functions which write
#    their output in a vector, in which case the vector type must
#    be strictly C-compatible.
#  * all char[] strings become Strings, both in inputs and in output.
#
#  A single exception to the strict compatibility is GLPK.version(),
#  which returns a tuple of integers in the form (major, minor)
#  rather than a string.

include("GLPK_checks.jl")

function set_prob_name(prob::Prob, name::Union(String,Nothing))
    @check! _prob(prob)
    if is(name, nothing)
        name = ""
    end
    @check _string_length(name, 255)
    @glpk_ccall set_prob_name Void (Ptr{Void}, Ptr{Uint8}) prob.p bytestring(name)
end

function set_obj_name(prob::Prob, name::Union(String,Nothing))
    @check! _prob(prob)
    if is(name, nothing)
        name = ""
    end
    @check _string_length(name, 255)
    @glpk_ccall set_obj_name Void (Ptr{Void}, Ptr{Uint8}) prob.p bytestring(name)
end

function set_row_name(prob::Prob, row::Integer, name::Union(String,Nothing))
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    if is(name, nothing)
        name = ""
    end
    @check _string_length(name, 255)
    @glpk_ccall set_row_name Void (Ptr{Void}, Int32, Ptr{Uint8}) prob.p row bytestring(name)
end

function set_col_name(prob::Prob, col::Integer, name::Union(String,Nothing))
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    if is(name, nothing)
        name = ""
    end
    @check _string_length(name, 255)
    @glpk_ccall set_col_name Void (Ptr{Void}, Int32, Ptr{Uint8}) prob.p col bytestring(name)
end

function set_obj_dir(prob::Prob, dir::Integer)
    @check! _prob(prob)
    @check _obj_dir_is_valid(dir)
    @glpk_ccall set_obj_dir Void (Ptr{Void}, Int32) prob.p dir
end

function add_rows(prob::Prob, rows::Integer)
    @check! _prob(prob)
    @glpk_ccall add_rows Int32 (Ptr{Void}, Int32) prob.p rows
end

function add_cols(prob::Prob, cols::Integer)
    @check! _prob(prob)
    @glpk_ccall add_cols Int32 (Ptr{Void}, Int32) prob.p cols
end

function set_row_bnds(prob::Prob, row::Integer, bounds_type::Integer, lb::Real, ub::Real)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    @check _bounds_type_is_valid(bounds_type)
    @check _bounds_are_valid(bounds_type, lb, ub)
    @glpk_ccall set_row_bnds Void (Ptr{Void}, Int32, Int32, Float64, Float64) prob.p row bounds_type lb ub
end

function set_col_bnds(prob::Prob, col::Integer, bounds_type::Integer, lb::Real, ub::Real)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @check _bounds_type_is_valid(bounds_type)
    @check _bounds_are_valid(bounds_type, lb, ub)
    @glpk_ccall set_col_bnds Void (Ptr{Void}, Int32, Int32, Float64, Float64) prob.p col bounds_type lb ub
end

function set_obj_coef(prob::Prob, col::Integer, coef::Real)
    @check! _prob(prob)
    @check _col_is_valid_w0(prob, col)
    @glpk_ccall set_obj_coef Void (Ptr{Void}, Int32, Float64) prob.p col coef
end

function set_mat_row{Ti<:Integer, Tv<:Real}(prob::Prob, row::Integer, len::Integer, ind::Vector{Ti}, val::Vector{Tv})
    @check! _prob(prob)
    @check! _vectors_size(len, ind, val)
    @check _row_is_valid(prob, row)
    if len > 0
        ind32 = int32(ind)
        val64 = float64(val)
        off32 = sizeof(Int32)
        off64 = sizeof(Float64)
        ind32p = pointer(ind32) - off32
        val64p = pointer(val64) - off64
        @check! _cols_ids_size(prob, 0, len, ind32)
        @check _cols_ids_content(prob, len, ind32)
    else
        ind32p = C_NULL
        val64p = C_NULL
    end

    @glpk_ccall set_mat_row Void (Ptr{Void}, Int32, Int32, Ptr{Int32}, Ptr{Float64}) prob.p row len ind32p val64p
end
function set_mat_row(prob::Prob, row::Integer, len::Integer, ind::VecOrNothing, val::VecOrNothing)
    ind = convert_vecornothing(Int32, ind)
    val = convert_vecornothing(Float64, ar)
    set_mat_row(prob, row, len, ind, val)
end
function set_mat_row(prob::Prob, row::Integer, ind::VecOrNothing, val::VecOrNothing)
    @check! _vectors_all_same_size(ind, val)
    l = vecornothing_length(ind)
    set_mat_row(prob, row, l, ind, val)
end


function set_mat_col{Ti<:Integer, Tv<:Real}(prob::Prob, col::Integer, len::Integer, ind::Vector{Ti}, val::Vector{Tv})
    @check! _prob(prob)
    @check! _vectors_size(len, ind, val)
    @check _col_is_valid(prob, col)
    if len > 0
        ind32 = int32(ind)
        val64 = float64(val)
        off32 = sizeof(Int32)
        off64 = sizeof(Float64)
        ind32p = pointer(ind32) - off32
        val64p = pointer(val64) - off64
        @check! _rows_ids_size(prob, 0, len, ind32)
        @check _rows_ids_content(prob, len, ind32)
    else
        ind32p = C_NULL
        val64p = C_NULL
    end

    @glpk_ccall set_mat_col Void (Ptr{Void}, Int32, Int32, Ptr{Int32}, Ptr{Float64}) prob.p col len ind32p val64p
end
function set_mat_col(prob::Prob, col::Integer, len::Integer, ind::VecOrNothing, val::VecOrNothing)
    ind = convert_vecornothing(Int32, ind)
    val = convert_vecornothing(Float64, ar)
    set_mat_col(prob, col, len, ind, val)
end
function set_mat_col(prob::Prob, col::Integer, ind::VecOrNothing, val::VecOrNothing)
    @check! _vectors_all_same_size(ind, val)
    l = vecornothing_length(ind)
    set_mat_col(prob, col, l, ind, val)
end

function load_matrix{Ti<:Integer, Tv<:Real}(prob::Prob, numel::Integer, ia::Vector{Ti}, ja::Vector{Ti}, ar::Vector{Tv})
    @check! _prob(prob)
    @check! _vectors_size(numel, ia, ja, ar)
    if numel == 0
        return
    end
    ia32 = int32(ia)
    ja32 = int32(ja)
    ar64 = float64(ar)
    @check _indices_vectors_dup(prob, numel, ia32, ja32)

    off32 = sizeof(Int32)
    off64 = sizeof(Float64)
    ia32p = pointer(ia32) - off32
    ja32p = pointer(ja32) - off32
    ar64p = pointer(ar64) - off64

    @glpk_ccall load_matrix Void (Ptr{Void}, Int32, Ptr{Int32}, Ptr{Int32}, Ptr{Float64}) prob.p numel ia32p ja32p ar64p
end

function load_matrix(prob::Prob, numel::Integer, ia::VecOrNothing, ja::VecOrNothing, ar::VecOrNothing)
    ia = convert_vecornothing(Int32, ia)
    ja = convert_vecornothing(Int32, ja)
    ar = convert_vecornothing(Float64, ar)
    load_matrix(prob, numel, ia, ja, ar)
end

function load_matrix(prob::Prob, ia::VecOrNothing, ja::VecOrNothing, ar::VecOrNothing)
    @check! _vectors_all_same_size(ia, ja, ar)
    l = vecornothing_length(ar)
    load_matrix(prob, l, ia, ja, ar)
end

function load_matrix{Ti<:Integer, Tv<:Real}(prob::Prob, a::SparseMatrixCSC{Tv, Ti})
    (ia, ja, ar) = findn_nzs(a)
    load_matrix(prob, ia, ja, ar)
end

function check_dup{Ti<:Integer}(rows::Integer, cols::Integer, numel::Integer, ia::Vector{Ti}, ja::Vector{Ti})
    @check _rows_and_cols(rows, cols)
    @check! _vectors_size(numel, ia, ja)
    ia32 = int32(ia)
    ja32 = int32(ja)

    off32 = sizeof(Int32)
    ia32p = pointer(ia32) - off32
    ja32p = pointer(ja32) - off32

    @glpk_ccall check_dup Int32 (Int32, Int32, Int32, Ptr{Int32}, Ptr{Int32}) rows cols numel ia32p ja32p
end

function check_dup(rows::Integer, cols::Integer, numel::Integer, ia::VecOrNothing, ja::VecOrNothing)
    ia = convert_vecornothing(Int32, ia)
    ja = convert_vecornothing(Int32, ja)
    check_dup(rows, cols, numel, ia, ja)
end

function check_dup(rows::Integer, cols::Integer, ia::VecOrNothing, ja::VecOrNothing)
    @check! _vectors_all_same_size(ia, ja)
    l = vecornothing_length(ia)
    check_dup(rows, cols, l, ia, ja)
end

function sort_matrix(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall sort_matrix Void (Ptr{Void},) prob.p
end

function del_rows{Ti<:Integer}(prob::Prob, num_rows::Integer, rows_ids::AbstractVector{Ti})
    @check! _prob(prob)
    rows_ids32 = int32(rows_ids)
    @check! _rows_ids_size(prob, 1, num_rows, rows_ids32)
    @check _rows_ids_content(prob, num_rows, rows_ids32)

    off32 = sizeof(Int32)
    rows_ids32p = pointer(rows_ids32) - off32
    @glpk_ccall del_rows Void (Ptr{Void}, Int32, Ptr{Int32}) prob.p num_rows rows_ids32p
end
del_rows{Ti<:Integer}(prob::Prob, rows_ids::AbstractVector{Ti}) =
    del_rows(prob, length(rows_ids), rows_ids)

function del_cols{Ti<:Integer}(prob::Prob, num_cols::Integer, cols_ids::AbstractVector{Ti})
    @check! _prob(prob)
    cols_ids32 = int32(cols_ids)
    @check! _cols_ids_size(prob, 1, num_cols, cols_ids32)
    @check _cols_ids_content(prob, num_cols, cols_ids32)

    off32 = sizeof(Int32)
    cols_ids32p = pointer(cols_ids32) - off32
    @glpk_ccall del_cols Void (Ptr{Void}, Int32, Ptr{Int32}) prob.p num_cols cols_ids32p
end
del_cols{Ti<:Integer}(prob::Prob, cols_ids::AbstractVector{Ti}) =
    del_cols(prob, length(cols_ids), cols_ids)

function copy_prob(prob_dest::Prob, prob::Prob, copy_names::Integer)
    @check! _prob(prob)
    @check _copy_names_flag(copy_names)
    @glpk_ccall copy_prob Void (Ptr{Void}, Ptr{Void}, Int32) prob_dest.p prob.p copy_names
end

function erase_prob(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall erase_prob Void (Ptr{Void},) prob.p
end

function get_prob_name(prob::Prob)
    @check! _prob(prob)
    name_cstr = @glpk_ccall get_prob_name Ptr{Uint8} (Ptr{Void},) prob.p
    if name_cstr == C_NULL
        return ""
    else
        return bytestring(name_cstr)
    end
end

function get_obj_name(prob::Prob)
    @check! _prob(prob)
    name_cstr = @glpk_ccall get_obj_name Ptr{Uint8} (Ptr{Void},) prob.p
    if name_cstr == C_NULL
        return ""
    else
        return bytestring(name_cstr)
    end
end

function get_obj_dir(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall get_obj_dir Int32 (Ptr{Void},) prob.p
end

function get_num_rows(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall get_num_rows Int32 (Ptr{Void},) prob.p
end

function get_num_cols(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall get_num_cols Int32 (Ptr{Void},) prob.p
end

function get_row_name(prob::Prob, row::Integer)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    name_cstr = @glpk_ccall get_row_name Ptr{Uint8} (Ptr{Void}, Int32) prob.p row
    if name_cstr == C_NULL
        return ""
    else
        return bytestring(name_cstr)
    end
end

function get_col_name(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    name_cstr = @glpk_ccall get_col_name Ptr{Uint8} (Ptr{Void}, Int32) prob.p col
    if name_cstr == C_NULL
        return ""
    else
        return bytestring(name_cstr)
    end
end

function get_row_type(prob::Prob, row::Integer)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    @glpk_ccall get_row_type Int32 (Ptr{Void}, Int32) prob.p row
end

function get_row_lb(prob::Prob, row::Integer)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    @glpk_ccall get_row_lb Float64 (Ptr{Void}, Int32) prob.p row
end

function get_row_ub(prob::Prob, row::Integer)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    @glpk_ccall get_row_ub Float64 (Ptr{Void}, Int32) prob.p row
end

function get_col_type(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @glpk_ccall get_col_type Int32 (Ptr{Void}, Int32) prob.p col
end

function get_col_lb(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @glpk_ccall get_col_lb Float64 (Ptr{Void}, Int32) prob.p col
end

function get_col_ub(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @glpk_ccall get_col_ub Float64 (Ptr{Void}, Int32) prob.p col
end

function get_obj_coef(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid_w0(prob, col)
    @glpk_ccall get_obj_coef Float64 (Ptr{Void}, Int32) prob.p col
end

function get_num_nz(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall get_num_nz Int32 (Ptr{Void},) prob.p
end

function get_mat_row(prob::Prob, row::Integer, ind::Union(Vector{Int32},Nothing), val::Union(Vector{Float64},Nothing))
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    numel = @glpk_ccall get_mat_row Int32 (Ptr{Void}, Int32, Ptr{Int32}, Ptr{Float64}) prob.p row C_NULL C_NULL
    if numel == 0
        return 0
    end
    if !isequal(ind, nothing)
        @check! _vectors_size(numel, ind)
        off32 = sizeof(Int32)
        ind32p = pointer(ind) - off32
    else
        ind32p = C_NULL
    end
    if !isequal(val, nothing)
        @check! _vectors_size(numel, val)
        off64 = sizeof(Float64)
        val64p = pointer(val) - off64
    else
        val64p = C_NULL
    end
    @glpk_ccall get_mat_row Int32 (Ptr{Void}, Int32, Ptr{Int32}, Ptr{Float64}) prob.p row ind32p val64p
end

function get_mat_row(prob::Prob, row::Integer)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    numel = @glpk_ccall get_mat_row Int32 (Ptr{Void}, Int32, Ptr{Int32}, Ptr{Float64}) prob.p row C_NULL C_NULL
    if numel == 0
        return (Int32[], Float64[])
    end
    ind = Array(Int32, numel)
    val = Array(Float64, numel)

    off32 = sizeof(Int32)
    ind32p = pointer(ind) - off32
    off64 = sizeof(Float64)
    val64p = pointer(val) - off64
    @glpk_ccall get_mat_row Int32 (Ptr{Void}, Int32, Ptr{Int32}, Ptr{Float64}) prob.p row ind32p val64p
    return ind, val
end

function get_mat_col(prob::Prob, col::Integer, ind::Union(Vector{Int32},Nothing), val::Union(Vector{Float64},Nothing))
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    numel = @glpk_ccall get_mat_col Int32 (Ptr{Void}, Int32, Ptr{Int32}, Ptr{Float64}) prob.p col C_NULL C_NULL
    if numel == 0
        return 0
    end
    if !isequal(ind, nothing)
        @check! _vectors_size(numel, ind)
        off32 = sizeof(Int32)
        ind32p = pointer(ind) - off32
    else
        ind32p = C_NULL
    end
    if !isequal(val, nothing)
        @check! _vectors_size(numel, val)
        off64 = sizeof(Float64)
        val64p = pointer(val) - off64
    else
        val64p = C_NULL
    end
    @glpk_ccall get_mat_col Int32 (Ptr{Void}, Int32, Ptr{Int32}, Ptr{Float64}) prob.p col ind32p val64p
end

function get_mat_col(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    numel = @glpk_ccall get_mat_col Int32 (Ptr{Void}, Int32, Ptr{Int32}, Ptr{Float64}) prob.p col C_NULL C_NULL
    if numel == 0
        return (Int32[], Float64[])
    end
    ind = Array(Int32, numel)
    val = Array(Float64, numel)

    off32 = sizeof(Int32)
    ind32p = pointer(ind) - off32
    off64 = sizeof(Float64)
    val64p = pointer(val) - off64
    @glpk_ccall get_mat_col Int32 (Ptr{Void}, Int32, Ptr{Int32}, Ptr{Float64}) prob.p col ind32p val64p
    return ind, val
end

function create_index(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall create_index Void (Ptr{Void},) prob.p
end

function find_row(prob::Prob, name::String)
    @check! _prob(prob)
    @glpk_ccall find_row Int32 (Ptr{Void}, Ptr{Uint8}) prob.p bytestring(name)
end

function find_col(prob::Prob, name::String)
    @check! _prob(prob)
    @glpk_ccall find_col Int32 (Ptr{Void}, Ptr{Uint8}) prob.p bytestring(name)
end

function delete_index(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall delete_index Void (Ptr{Void},) prob.p
end

function set_rii(prob::Prob, row::Integer, rii::Real)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    @glpk_ccall set_rii Void (Ptr{Void}, Int32, Float64) prob.p row rii
end

function set_sjj(prob::Prob, col::Integer, sjj::Real)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @glpk_ccall set_sjj Void (Ptr{Void}, Int32, Float64) prob.p col sjj
end

function get_rii(prob::Prob, row::Integer)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    @glpk_ccall get_rii Float64 (Ptr{Void}, Int32) prob.p row
end

function get_sjj(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @glpk_ccall get_sjj Float64 (Ptr{Void}, Int32) prob.p col
end

function scale_prob(prob::Prob, flags::Integer)
    @check! _prob(prob)
    @check _scale_flags(flags)
    @glpk_ccall scale_prob Void (Ptr{Void}, Int32) prob.p flags
end

function unscale_prob(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall unscale_prob Void (Ptr{Void},) prob.p
end

function set_row_stat(prob::Prob, row::Integer, stat::Integer)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    @check _stat_is_valid(stat)
    @glpk_ccall set_row_stat Void (Ptr{Void}, Int32, Int32) prob.p row stat
end

function set_col_stat(prob::Prob, col::Integer, stat::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @check _stat_is_valid(stat)
    @glpk_ccall set_col_stat Void (Ptr{Void}, Int32, Int32) prob.p col stat
end

function std_basis(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall std_basis Void (Ptr{Void},) prob.p
end

function adv_basis(prob::Prob, flags::Integer)
    @check! _prob(prob)
    @check _adv_basis_flags(flags)
    @glpk_ccall adv_basis Void (Ptr{Void}, Int32) prob.p flags
end
adv_basis(prob::Prob) = adv_basis(prob, 0)

function cpx_basis(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall cpx_basis Void (Ptr{Void},) prob.p
end

function simplex{Tp<:Union(SimplexParam, Nothing)}(prob::Prob, param::Tp)
    @check! _prob(prob)
    if param == nothing
        return @glpk_ccall simplex Int32 (Ptr{Void}, Ptr{Void}) prob.p C_NULL
    else
        return @glpk_ccall simplex Int32 (Ptr{Void}, Ptr{SimplexParam}) prob.p &param
    end
end

simplex(prob::Prob) =
    simplex(prob, nothing)

function exact{Tp<:Union(SimplexParam, Nothing)}(prob::Prob, param::Tp)
    @check! _prob(prob)
    if param == nothing
        return @glpk_ccall exact Int32 (Ptr{Void}, Ptr{Void}) prob.p C_NULL
    else
        return @glpk_ccall exact Int32 (Ptr{Void}, Ptr{SimplexParam}) prob.p &param
    end
end

exact(prob::Prob) =
    exact(prob, nothing)

function init_smcp(param::SimplexParam)
    @glpk_ccall init_smcp Int32 (Ptr{SimplexParam},) &param
end

function get_status(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall get_status Int32 (Ptr{Void},) prob.p
end

function get_prim_stat(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall get_prim_stat Int32 (Ptr{Void},) prob.p
end

function get_dual_stat(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall get_dual_stat Int32 (Ptr{Void},) prob.p
end

function get_obj_val(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall get_obj_val Float64 (Ptr{Void},) prob.p
end

function get_row_stat(prob::Prob, row::Integer)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    @glpk_ccall get_row_stat Int32 (Ptr{Void}, Int32) prob.p row
end

function get_row_prim(prob::Prob, row::Integer)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    @glpk_ccall get_row_prim Float64 (Ptr{Void}, Int32) prob.p row
end

function get_row_dual(prob::Prob, row::Integer)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    @glpk_ccall get_row_dual Float64 (Ptr{Void}, Int32) prob.p row
end

function get_col_stat(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @glpk_ccall get_col_stat Int32 (Ptr{Void}, Int32) prob.p col
end

function get_col_prim(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @glpk_ccall get_col_prim Float64 (Ptr{Void}, Int32) prob.p col
end

function get_col_dual(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @glpk_ccall get_col_dual Float64 (Ptr{Void}, Int32) prob.p col
end

function get_unbnd_ray(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall get_unbnd_ray Int32 (Ptr{Void},) prob.p
end

function interior{Tp<:Union(InteriorParam, Nothing)}(prob::Prob, param::Tp)
    @check! _prob(prob)
    if param == nothing
        return @glpk_ccall interior Int32 (Ptr{Void}, Ptr{Void}) prob.p C_NULL
    else
        return @glpk_ccall interior Int32 (Ptr{Void}, Ptr{InteriorParam}) prob.p &param
    end
end

interior(prob::Prob) = interior(prob, nothing)

function init_iptcp(param::InteriorParam)
    @glpk_ccall init_iptcp Int32 (Ptr{InteriorParam},) &param
end

function ipt_status(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall ipt_status Int32 (Ptr{Void},) prob.p
end

function ipt_obj_val(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall ipt_obj_val Float64 (Ptr{Void},) prob.p
end

function ipt_row_prim(prob::Prob, row::Integer)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    @glpk_ccall ipt_row_prim Float64 (Ptr{Void}, Int32) prob.p row
end

function ipt_row_dual(prob::Prob, row::Integer)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    @glpk_ccall ipt_row_dual Float64 (Ptr{Void}, Int32) prob.p row
end

function ipt_col_prim(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @glpk_ccall ipt_col_prim Float64 (Ptr{Void}, Int32) prob.p col
end

function ipt_col_dual(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @glpk_ccall ipt_col_dual Float64 (Ptr{Void}, Int32) prob.p col
end

function set_col_kind(prob::Prob, col::Integer, kind::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @check _kind_is_valid(kind)
    @glpk_ccall set_col_kind Void (Ptr{Void}, Int32, Int32) prob.p col kind
end

function get_col_kind(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @glpk_ccall get_col_kind Int32 (Ptr{Void}, Int32) prob.p col
end

function get_num_int(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall get_num_int Int32 (Ptr{Void},) prob.p
end

function get_num_bin(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall get_num_bin Int32 (Ptr{Void},) prob.p
end

function intopt{Tp<:Union(IntoptParam, Nothing)}(prob::Prob, param::Tp)
    @check! _prob(prob)
    if param == nothing
        return @glpk_ccall intopt Int32 (Ptr{Void}, Ptr{Void}) prob.p C_NULL
    else
        return @glpk_ccall intopt Int32 (Ptr{Void}, Ptr{IntoptParam}) prob.p &param
    end
end

intopt(prob::Prob) = intopt(prob, nothing)

function init_iocp(param::IntoptParam)
    @glpk_ccall init_iocp Int32 (Ptr{IntoptParam},) &param
end

function mip_status(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall mip_status Int32 (Ptr{Void},) prob.p
end

function mip_obj_val(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall mip_obj_val Float64 (Ptr{Void},) prob.p
end

function mip_row_val(prob::Prob, row::Integer)
    @check! _prob(prob)
    @check _row_is_valid(prob, row)
    @glpk_ccall mip_row_val Float64 (Ptr{Void}, Int32) prob.p row
end

function mip_col_val(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _col_is_valid(prob, col)
    @glpk_ccall mip_col_val Float64 (Ptr{Void}, Int32) prob.p col
end

#TODO
#function lpx_check_kkt(prob::Prob, scaled::Integer, kkt)

function read_mps(prob::Prob, format::Integer, param, filename::String)
    @check! _prob(prob)
    @check _mps_format(format)
    if is(param, nothing)
        param = C_NULL
    else
        @check _mps_param(param)
    end

    @check _file_is_readable(filename)
    ret = @glpk_ccall read_mps Int32 (Ptr{Void}, Int32, Ptr{Void}, Ptr{Uint8}) prob.p format param bytestring(filename)
    @check! _succeeded(ret, "read_mps")
    return ret
end

read_mps(prob::Prob, format::Integer, filename::String) =
    read_mps(prob, format, C_NULL, filename)

function write_mps(prob::Prob, format::Integer, param, filename::String)
    @check! _prob(prob)
    @check _mps_format(format)
    if is(param, nothing)
        param = C_NULL
    else
        @check _mps_param(param)
    end
    @check _file_is_writable(filename)
    ret = @glpk_ccall write_mps Int32 (Ptr{Void}, Int32, Ptr{Void}, Ptr{Uint8}) prob.p format param bytestring(filename)
    @check! _succeeded(ret, "write_mps")
    return ret
end

write_mps(prob::Prob, format::Integer, filename::String) =
    write_mps(prob, format, C_NULL, filename)

function read_lp(prob::Prob, param, filename::String)
    @check! _prob(prob)
    @check _lp_param(param)
    @check _file_is_readable(filename)
    ret = @glpk_ccall read_lp Int32 (Ptr{Void}, Ptr{Void}, Ptr{Uint8}) prob.p param bytestring(filename)
    @check! _succeeded(ret, "read_lp")
    return ret
end

read_lp(prob::Prob, filename::String) =
    read_lp(prob, C_NULL, filename)

function write_lp(prob::Prob, param, filename::String)
    @check! _prob(prob)
    @check _lp_param(param)
    @check _file_is_writable(filename)
    ret = @glpk_ccall write_lp Int32 (Ptr{Void}, Ptr{Void}, Ptr{Uint8}) prob.p param bytestring(filename)
    @check! _succeeded(ret, "write_lp")
    return ret
end

write_lp(prob::Prob, filename::String) =
    write_lp(prob, C_NULL, filename)

function read_prob(prob::Prob, flags::Integer, filename::String)
    @check! _prob(prob)
    @check _read_prob_flags(flags)
    @check _file_is_readable(filename)
    ret = @glpk_ccall read_prob Int32 (Ptr{Void}, Int32, Ptr{Uint8}) prob.p flags bytestring(filename)
    @check! _succeeded(ret, "read_prob")
    return ret
end

read_prob(prob::Prob, filename::String) =
    read_prob(prob, 0, filename)

function write_prob(prob::Prob, flags::Integer, filename::String)
    @check! _prob(prob)
    @check _write_prob_flags(flags)
    @check _file_is_writable(filename)
    ret = @glpk_ccall write_prob Int32 (Ptr{Void}, Int32, Ptr{Uint8}) prob.p flags bytestring(filename)
    @check! _succeeded(ret, "write_prob")
    return ret
end

write_prob(prob::Prob, filename::String) =
    write_prob(prob, 0, filename)

function mpl_read_model(tran::MathProgWorkspace, filename::String, skip::Integer)
    @check! _mpl_workspace(tran)
    @check _file_is_readable(filename)
    ret = @glpk_ccall mpl_read_model Int32 (Ptr{Void}, Ptr{Uint8}, Int32) tran.p bytestring(filename) skip
    @check! _succeeded(ret, "mpl_read_model")
    return ret
end

function mpl_read_data(tran::MathProgWorkspace, filename::String)
    @check! _mpl_workspace(tran)
    @check _file_is_readable(filename)
    ret = @glpk_ccall mpl_read_data Int32 (Ptr{Void}, Ptr{Uint8}) tran.p bytestring(filename)
    @check! _succeeded(ret, "mpl_read_data")
    return ret
end

function mpl_generate(tran::MathProgWorkspace, filename::Union(String, Nothing))
    @check! _mpl_workspace(tran)
    if is(filename, nothing)
        cfilename = C_NULL
    else
        @check _file_is_writable(filename)
        cfilename = bytestring(filename)
    end
    ret = @glpk_ccall mpl_generate Int32 (Ptr{Void}, Ptr{Uint8}) tran.p cfilename
    @check! _succeeded(ret, "mpl_generate")
    return ret

end
mpl_generate(tran::MathProgWorkspace) = mpl_generate(tran, nothing)

function mpl_build_prob(tran::MathProgWorkspace, prob::Prob)
    @check! _mpl_workspace(tran)
    @check! _prob(prob)
    @glpk_ccall mpl_build_prob Void (Ptr{Void}, Ptr{Void}) tran.p prob.p
end

function mpl_postsolve(tran::MathProgWorkspace, prob::Prob, sol::Integer)
    @check! _mpl_workspace(tran)
    @check! _prob(prob)
    @check _mpl_postsolve_param(sol)
    ret = @glpk_ccall mpl_postsolve Int32 (Ptr{Void}, Ptr{Void}, Int32) tran.p prob.p sol
    @check! _succeeded(ret, "mpl_postsolve")
    return ret
end

function print_sol(prob::Prob, filename::String)
    @check! _prob(prob)
    @check _file_is_writable(filename)
    ret = @glpk_ccall print_sol Int32 (Ptr{Void}, Ptr{Uint8}) prob.p bytestring(filename)
    @check! _succeeded(ret, "print_sol")
    return ret
end

function read_sol(prob::Prob, filename::String)
    @check! _prob(prob)
    @check _file_is_readable(filename)
    ret = @glpk_ccall read_sol Int32 (Ptr{Void}, Ptr{Uint8}) prob.p bytestring(filename)
    @check! _succeeded(ret, "read_sol")
    return ret
end

function write_sol(prob::Prob, filename::String)
    @check! _prob(prob)
    @check _file_is_writable(filename)
    ret = @glpk_ccall write_sol Int32 (Ptr{Void}, Ptr{Uint8}) prob.p bytestring(filename)
    @check! _succeeded(ret, "write_sol")
    return ret
end

function print_ipt(prob::Prob, filename::String)
    @check! _prob(prob)
    @check _file_is_writable(filename)
    ret = @glpk_ccall print_ipt Int32 (Ptr{Void}, Ptr{Uint8}) prob.p bytestring(filename)
    @check! _succeeded(ret, "print_ipt")
    return ret
end

function read_ipt(prob::Prob, filename::String)
    @check! _prob(prob)
    @check _file_is_readable(filename)
    ret = @glpk_ccall read_ipt Int32 (Ptr{Void}, Ptr{Uint8}) prob.p bytestring(filename)
    @check! _succeeded(ret, "read_ipt")
    return ret
end

function write_ipt(prob::Prob, filename::String)
    @check! _prob(prob)
    @check _file_is_writable(filename)
    ret = @glpk_ccall write_ipt Int32 (Ptr{Void}, Ptr{Uint8}) prob.p bytestring(filename)
    @check! _succeeded(ret, "write_ipt")
    return ret
end

function print_mip(prob::Prob, filename::String)
    @check! _prob(prob)
    @check _file_is_writable(filename)
    ret = @glpk_ccall print_mip Int32 (Ptr{Void}, Ptr{Uint8}) prob.p bytestring(filename)
    @check! _succeeded(ret, "print_mip")
    return ret
end

function read_mip(prob::Prob, filename::String)
    @check! _prob(prob)
    @check _file_is_readable(filename)
    ret = @glpk_ccall read_mip Int32 (Ptr{Void}, Ptr{Uint8}) prob.p bytestring(filename)
    @check! _succeeded(ret, "read_mip")
    return ret
end

function write_mip(prob::Prob, filename::String)
    @check! _prob(prob)
    @check _file_is_writable(filename)
    ret = @glpk_ccall write_mip Int32 (Ptr{Void}, Ptr{Uint8}) prob.p bytestring(filename)
    @check! _succeeded(ret, "write_mip")
    return ret
end

function print_ranges{Ti<:Integer}(prob::Prob, len::Integer, list::Vector{Ti}, flags::Integer, filename::String)
    @check! _prob(prob)
    @check! _vectors_size(len, list)
    @check _status_is_optimal(prob)
    @check _bf_exists(prob)
    @check _print_ranges_flags(flags)
    @check _file_is_writable(filename)

    if len > 0
        list32 = int32(list)
        @check _list_ids(prob, len, list32)

        off32 = sizeof(Int32)
        list32p = pointer(list32) - off32
    else
        list32p = C_NULL
    end

    @glpk_ccall print_ranges Int32 (Ptr{Void}, Int32, Ptr{Int32}, Int32, Ptr{Uint8}) prob.p len list32p flags bytestring(filename)
end

print_ranges{Ti<:Integer}(prob::Prob, list::Vector{Ti}, flags::Integer, filename::String) =
    print_ranges(prob, length(list), list, flags, filename)

print_ranges{Ti<:Integer}(prob::Prob, len::Integer, list::Vector{Ti}, filename::String) = 
    print_ranges(prob, len, list, 0, filename)

print_ranges{Ti<:Integer}(prob::Prob, list::Vector{Ti}, filename::String) =
    print_ranges(prob, length(list), list, 0, filename)

function print_ranges(prob::Prob, len::Integer, list::VecOrNothing, flags::Integer, filename::String)
    list = convert_vecornothing(Int32, list)
    print_ranges(prob, len, list, flags, filename)
end

print_ranges(prob::Prob, list::VecOrNothing, flags::Integer, filename::String) =
    print_ranges(prob, vecornothing_length(list), list, flags, filename)

print_ranges(prob::Prob, len::Integer, list::VecOrNothing, filename::String) =
    print_ranges(prob, len, list, 0, filename)

print_ranges(prob::Prob, list::VecOrNothing, filename::String) =
    print_ranges(prob, vecornothing_length(list), list, 0, filename)

print_ranges(prob::Prob, filename::String) =
    print_ranges(prob, 0, nothing, 0, filename)

function bf_exists(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall bf_exists Int32 (Ptr{Void},) prob.p
end

function factorize(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall factorize Int32 (Ptr{Void},) prob.p
end

function bf_updated(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall bf_updated Int32 (Ptr{Void},) prob.p
end

function get_bfcp(prob::Prob, param::BasisFactParam)
    @check! _prob(prob)
    @glpk_ccall get_bfcp Void (Ptr{Void}, Ptr{BasisFactParam}) prob.p &param
end

function set_bfcp(prob::Prob, param::Union(BasisFactParam,Nothing))
    @check! _prob(prob)
    if is(param, nothing)
        return @glpk_ccall set_bfcp Void (Ptr{Void}, Ptr{Void}) prob.p C_NULL
    else
        return @glpk_ccall set_bfcp Void (Ptr{Void}, Ptr{BasisFactParam}) prob.p &param
    end
end
set_bfcp(prob::Prob) = set_bfcp(prob, nothing)

function get_bhead(prob::Prob, k::Integer)
    @check! _prob(prob)
    @check _bf_exists(prob)
    @check _row_is_valid(prob, k)
    @glpk_ccall get_bhead Int32 (Ptr{Void}, Int32) prob.p k
end

function get_row_bind(prob::Prob, row::Integer)
    @check! _prob(prob)
    @check _bf_exists(prob)
    @check _row_is_valid(prob, row)
    @glpk_ccall get_row_bind Int32 (Ptr{Void}, Int32) prob.p row
end

function get_col_bind(prob::Prob, col::Integer)
    @check! _prob(prob)
    @check _bf_exists(prob)
    @check _col_is_valid(prob, col)
    @glpk_ccall get_col_bind Int32 (Ptr{Void}, Int32) prob.p col
end

function ftran(prob::Prob, x::Vector{Float64})
    @check! _prob(prob)
    rows = @glpk_ccall get_num_rows Int32 (Ptr{Void},) prob.p
    @check! _vectors_size(rows, x)
    off64 = sizeof(Float64)
    x64p = pointer(x) - off64
    @glpk_ccall ftran Void (Ptr{Void}, Ptr{Float64}) prob.p x64p
end

function btran(prob::Prob, x::Vector{Float64})
    @check! _prob(prob)
    rows = @glpk_ccall get_num_rows Int32 (Ptr{Void},) prob.p
    @check! _vectors_size(rows, x)
    off64 = sizeof(Float64)
    x64p = pointer(x) - off64
    @glpk_ccall btran Void (Ptr{Void}, Ptr{Float64}) prob.p x64p
end

function warm_up(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall warm_up Int32 (Ptr{Void},) prob.p
end

function eval_tab_row(prob::Prob, k::Integer, ind::Vector{Int32}, val::Vector{Float64})
    @check! _prob(prob)
    @check _bf_exists(prob)
    rows = @glpk_ccall get_num_rows Int32 (Ptr{Void},) prob.p
    cols = @glpk_ccall get_num_cols Int32 (Ptr{Void},) prob.p
    k_max = rows + cols
    @check _rowcol_is_valid(k, k_max)
    @check _var_is_basic(prob, k)

    grow!(ind, k_max - length(ind))
    grow!(val, k_max - length(val))

    off32 = sizeof(Int32)
    off64 = sizeof(Float64)
    ind32p = pointer(ind) - off32
    val64p = pointer(val) - off64

    len = @glpk_ccall eval_tab_row Int32 (Ptr{Void}, Int32, Ptr{Int32}, Ptr{Float64}) prob.p k ind32p val64p

    delete!(ind, len+1:length(ind))
    delete!(val, len+1:length(val))

    return len
end

function eval_tab_row(prob::Prob, k::Integer)
    @check! _prob(prob)
    @check _bf_exists(prob)
    rows = @glpk_ccall get_num_rows Int32 (Ptr{Void},) prob.p
    cols = @glpk_ccall get_num_cols Int32 (Ptr{Void},) prob.p
    k_max = rows + cols
    @check _rowcol_is_valid(k, k_max)
    @check _var_is_basic(prob, k)

    ind = Array(Int32, k_max)
    val = Array(Float64, k_max)

    off32 = sizeof(Int32)
    off64 = sizeof(Float64)
    ind32p = pointer(ind) - off32
    val64p = pointer(val) - off64

    len = @glpk_ccall eval_tab_row Int32 (Ptr{Void}, Int32, Ptr{Int32}, Ptr{Float64}) prob.p k ind32p val64p

    delete!(ind, len+1:length(ind))
    delete!(val, len+1:length(val))

    return ind, val
end

function eval_tab_col(prob::Prob, k::Integer, ind::Vector{Int32}, val::Vector{Float64})
    @check! _prob(prob)
    @check _bf_exists(prob)
    rows = @glpk_ccall get_num_rows Int32 (Ptr{Void},) prob.p
    cols = @glpk_ccall get_num_cols Int32 (Ptr{Void},) prob.p
    k_max = rows + cols
    @check _rowcol_is_valid(k, k_max)
    @check _var_is_non_basic(prob, k)

    grow!(ind, k_max - length(ind))
    grow!(val, k_max - length(val))

    off32 = sizeof(Int32)
    off64 = sizeof(Float64)
    ind32p = pointer(ind) - off32
    val64p = pointer(val) - off64

    len = @glpk_ccall eval_tab_col Int32 (Ptr{Void}, Int32, Ptr{Int32}, Ptr{Float64}) prob.p k ind32p val64p

    delete!(ind, len+1:length(ind))
    delete!(val, len+1:length(val))

    return len
end

function eval_tab_col(prob::Prob, k::Integer)
    @check! _prob(prob)
    @check _bf_exists(prob)
    rows = @glpk_ccall get_num_rows Int32 (Ptr{Void},) prob.p
    cols = @glpk_ccall get_num_cols Int32 (Ptr{Void},) prob.p
    k_max = rows + cols
    @check _rowcol_is_valid(k, k_max)
    @check _var_is_non_basic(prob, k)

    ind = Array(Int32, k_max)
    val = Array(Float64, k_max)

    off32 = sizeof(Int32)
    off64 = sizeof(Float64)
    ind32p = pointer(ind) - off32
    val64p = pointer(val) - off64

    len = @glpk_ccall eval_tab_col Int32 (Ptr{Void}, Int32, Ptr{Int32}, Ptr{Float64}) prob.p k ind32p val64p

    delete!(ind, len+1:length(ind))
    delete!(val, len+1:length(val))

    return ind, val
end

function transform_row(prob::Prob, len::Integer, ind::Vector{Int32}, val::Vector{Float64})
    @check! _prob(prob)
    @check _bf_exists(prob)
    @check _col_is_valid(prob, len)
    @check! _vectors_size(len, ind, val)

    cols = @glpk_ccall get_num_cols Int32 (Ptr{Void},) prob.p

    grow!(ind, cols - length(ind))
    grow!(val, cols - length(val))

    off32 = sizeof(Int32)
    off64 = sizeof(Float64)
    ind32p = pointer(ind) - off32
    val64p = pointer(val) - off64

    len1 = @glpk_ccall transform_row Int32 (Ptr{Void}, Int32, Ptr{Int32}, Ptr{Float64}) prob.p len ind32p val64p

    delete!(ind, len1+1:length(ind))
    delete!(val, len1+1:length(val))

    return len1
end

function transform_row(prob::Prob, ind::Vector{Int32}, val::Vector{Float64})
    @check! _vectors_all_same_size(ind, val)
    transform_row(prob, length(ind), ind, val)
end

function transform_col(prob::Prob, len::Integer, ind::Vector{Int32}, val::Vector{Float64})
    @check! _prob(prob)
    @check _bf_exists(prob)
    @check _row_is_valid(prob, len)
    @check! _vectors_size(len, ind, val)

    rows = @glpk_ccall get_num_rows Int32 (Ptr{Void},) prob.p

    grow!(ind, rows - length(ind))
    grow!(val, rows - length(val))

    off32 = sizeof(Int32)
    off64 = sizeof(Float64)
    ind32p = pointer(ind) - off32
    val64p = pointer(val) - off64

    len1 = @glpk_ccall transform_col Int32 (Ptr{Void}, Int32, Ptr{Int32}, Ptr{Float64}) prob.p len ind32p val64p

    delete!(ind, len1+1:length(ind))
    delete!(val, len1+1:length(val))

    return len1
end

function transform_col(prob::Prob, ind::Vector{Int32}, val::Vector{Float64})
    @check! _vectors_all_same_size(ind, val)
    transform_col(prob, length(ind), ind, val)
end

function prim_rtest{Ti<:Integer, Tv<:Real}(prob::Prob, len::Integer, ind::Vector{Ti}, val::Vector{Tv}, dir::Integer, eps::Real)
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

    ind32 = int32(ind)
    val64 = float64(val)
    off32 = sizeof(Int32)
    off64 = sizeof(Float64)
    ind32p = pointer(ind32) - off32
    val64p = pointer(val64) - off64

    piv = @glpk_ccall prim_rtest Int32 (Ptr{Void}, Int32, Ptr{Int32}, Ptr{Float64}, Int32, Float64) prob.p len ind32p val64p dir eps
    return piv
end

function prim_rtest{Ti<:Integer, Tv<:Real}(prob::Prob, ind::Vector{Ti}, val::Vector{Tv}, dir::Integer, eps::Real)
    @check! _vectors_all_same_size(ind, val)
    prim_rtest(prob, length(ind), ind, val, dir, eps)
end

function dual_rtest{Ti<:Integer, Tv<:Real}(prob::Prob, len::Integer, ind::Vector{Ti}, val::Vector{Tv}, dir::Integer, eps::Real)
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

    ind32 = int32(ind)
    val64 = float64(val)
    off32 = sizeof(Int32)
    off64 = sizeof(Float64)
    ind32p = pointer(ind32) - off32
    val64p = pointer(val64) - off64

    piv = @glpk_ccall dual_rtest Int32 (Ptr{Void}, Int32, Ptr{Int32}, Ptr{Float64}, Int32, Float64) prob.p len ind32p val64p dir eps
    return piv
end

function dual_rtest{Ti<:Integer, Tv<:Real}(prob::Prob, ind::Vector{Ti}, val::Vector{Tv}, dir::Integer, eps::Real)
    @check! _vectors_all_same_size(ind, val)
    dual_rtest(prob, length(ind), ind, val, dir, eps)
end

function analyze_bound(prob::Prob, k, limit1, var1, limit2, var2)
    error("unsupported. Use GLPK.analyze_bound(prob, k) instead.")
end

function analyze_bound(prob::Prob, k::Integer)
    @check! _prob(prob)
    @check _bf_exists(prob)
    @check _rowcol_is_valid(prob, k)
    @check _var_is_non_basic(prob, k)

    limit1 = Array(Float64, 1)
    var1 = Array(Int32, 1)
    limit2 = Array(Float64, 1)
    var2 = Array(Int32, 1)

    @glpk_ccall analyze_bound Void (Ptr{Void}, Int32, Ptr{Float64}, Ptr{Int32}, Ptr{Float64}, Ptr{Int32}) prob.p k pointer(limit1) pointer(var1) pointer(limit2) pointer(var2)

    return limit1[1], var1[1], limit2[1], var2[1]
end

function analyze_coef(prob::Prob, k, coef1, var1, value1, coef2, var2, value2)
    error("unsupported. Use GLPK.analyze_coef(prob, k) instead.")
end

function analyze_coef(prob::Prob, k::Integer)
    @check! _prob(prob)
    @check _bf_exists(prob)
    @check _rowcol_is_valid(prob, k)
    @check _var_is_basic(prob, k)

    coef1 = Array(Float64, 1)
    var1 = Array(Int32, 1)
    value1 = Array(Float64, 1)
    coef2 = Array(Float64, 1)
    var2 = Array(Int32, 1)
    value2 = Array(Float64, 1)

    @glpk_ccall analyze_coef Void (Ptr{Void}, Int32, Ptr{Float64}, Ptr{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Int32}, Ptr{Float64}) prob.p k pointer(coef1) pointer(var1) pointer(value1) pointer(coef2) pointer(var2) pointer(value2)

    return coef1[1], var1[1], value1[1], coef2[1], var2[1], value2[1] 
end

function ios_reason(tree::Ptr{Void})
    @check! _tree(tree)
    @glpk_ccall ios_reason Int32 (Ptr{Void},) tree
end

function ios_get_prob(tree::Ptr{Void})
    @check! _tree(tree)
    p = @glpk_ccall ios_get_prob Ptr{Void} (Ptr{Void},) tree
    return Prob(p)
end

function ios_row_attr(tree::Ptr{Void}, row::Integer, attr::Attr)
    @check! _tree(tree)
    prob = ios_get_prob(tree)
    @check _row_is_valid(prob, row)
    @glpk_ccall ios_row_attr Void (Ptr{Void}, Int32, Ptr{Attr}) tree row &attr
end

function ios_row_attr(tree::Ptr{Void}, row::Integer)
    @check! _tree(tree)
    prob = ios_get_prob(tree)
    @check _row_is_valid(prob, row)
    attr = Attr()
    @glpk_ccall ios_row_attr Void (Ptr{Void}, Int32, Ptr{Attr}) tree row &attr
    return attr
end

function ios_mip_gap(tree::Ptr{Void})
    @check! _tree(tree)
    @glpk_ccall ios_mip_gap Float64 (Ptr{Void},) tree
end

function ios_node_data(tree::Ptr{Void}, p::Integer)
    @check! _tree(tree)
    @check _ios_node_is_valid(tree, p)
    @glpk_ccall ios_node_data Ptr{Void} (Ptr{Void}, Int32) tree p
end

function ios_select_node(tree::Ptr{Void}, p::Integer)
    @check! _tree(tree)
    @check _reason(tree, [GLPK.ISELECT])
    @check _ios_node_is_active(tree, p)
    @glpk_ccall ios_select_node Void (Ptr{Void}, Int32) tree p
end

function ios_heur_sol{Tv<:Real}(tree::Ptr{Void}, x::Vector{Tv})
    @check! _tree(tree)
    prob = ios_get_prob(tree)
    cols = @glpk_ccall get_num_cols Int32 (Ptr{Void},) prob.p
    @check! _vectors_size(cols, x)
    x64 = float64(x)
    off64 = sizeof(Float64)
    x64p = pointer(x64) - off64

    @glpk_ccall ios_heur_sol Int32 (Ptr{Void}, Ptr{Float64}) tree x64p
end

function ios_can_branch(tree::Ptr{Void}, col::Integer)
    @check! _tree(tree)
    prob = ios_get_prob(tree)
    @check _col_is_valid(prob, col)
    @glpk_ccall ios_can_branch Int32 (Ptr{Void}, Int32) tree col
end

function ios_branch_upon(tree::Ptr{Void}, col::Integer, sel::Integer)
    @check! _tree(tree)
    prob = ios_get_prob(tree)
    @check _col_is_valid(prob, col)
    @check _can_branch(tree, col)
    @check _sel_is_valid(sel)
    @glpk_ccall ios_branch_upon Void (Ptr{Void}, Int32, Int32) tree col sel
end

function ios_terminate(tree::Ptr{Void})
    @check! _tree(tree)
    @glpk_ccall ios_terminate Void (Ptr{Void},) tree
end

function ios_tree_size(tree::Ptr{Void})
    @check! _tree(tree)
    a_cnt = Array(Int32, 1)
    n_cnt = Array(Int32, 1)
    t_cnt = Array(Int32, 1)
    @glpk_ccall ios_tree_size Void (Ptr{Void}, Ptr{Int32}, Ptr{Int32}, Ptr{Int32}) tree pointer(a_cnt) pointer(n_cnt) pointer(t_cnt)

    return a_cnt[1], n_cnt[1], t_cnt[1]
end

function ios_tree_size(tree::Ptr{Void}, a_cnt, n_cnt, t_cnt)
    error("unsupported. Use GLPK.ios_tree_size(tree) instead.")
end

function ios_curr_node(tree::Ptr{Void})
    @check! _tree(tree)
    @glpk_ccall ios_curr_node Int32 (Ptr{Void},) tree
end

function ios_next_node(tree::Ptr{Void}, p::Integer)
    @check! _tree(tree)
    if p != 0
        @check _ios_node_is_active(tree, p)
    end
    @glpk_ccall ios_next_node Int32 (Ptr{Void}, Int32) tree p
end

function ios_prev_node(tree::Ptr{Void}, p::Integer)
    @check! _tree(tree)
    if p != 0
        @check _ios_node_is_active(tree, p)
    end
    @glpk_ccall ios_prev_node Int32 (Ptr{Void}, Int32) tree p
end

function ios_up_node(tree::Ptr{Void}, p::Integer)
    @check! _tree(tree)
    @check _ios_node_is_valid(tree, p)
    @glpk_ccall ios_up_node Int32 (Ptr{Void}, Int32) tree p
end

function ios_node_level(tree::Ptr{Void}, p::Integer)
    @check! _tree(tree)
    @check _ios_node_is_valid(tree, p)
    @glpk_ccall ios_node_level Int32 (Ptr{Void}, Int32) tree p
end

function ios_node_bound(tree::Ptr{Void}, p::Integer)
    @check! _tree(tree)
    @check _ios_node_is_valid(tree, p)
    @glpk_ccall ios_node_bound Float64 (Ptr{Void}, Int32) tree p
end

function ios_best_node(tree::Ptr{Void})
    @check! _tree(tree)
    @glpk_ccall ios_best_node Int32 (Ptr{Void},) tree
end

function ios_pool_size(tree::Ptr{Void})
    @check! _tree(tree)
    @check _reason(tree, [GLPK.ICUTGEN])
    @glpk_ccall ios_pool_size Int32 (Ptr{Void},) tree
end

function ios_add_row{Ti<:Integer, Tv<:Real}(tree::Ptr{Void}, name::Union(String,Nothing),
                    klass::Integer, flags::Integer, len::Integer, ind::Vector{Ti}, val::Vector{Tv},
                    constr_type::Integer, rhs::Float64)

    @check! _tree(tree)
    @check _reason(tree, [GLPK.ICUTGEN])
    if name === nothing
        name = ""
    end
    @check _string_length(name, 255)
    @check _klass_is_valid(klass)
    @check _ios_add_row_flags(flags)
    @check! _vectors_size(len, ind, val)
    if len > 0
        ind32 = int32(ind)
        val64 = float64(val)
        off32 = sizeof(Int32)
        off64 = sizeof(Float64)
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

    @glpk_ccall ios_add_row Int32 (Ptr{Void}, Ptr{Uint8}, Int32, Int32, Int32, Ptr{Int32}, Ptr{Float64}, Int32, Float64) tree bytestring(name) klass flags len ind32p val64p constr_type rhs
end

function ios_add_row{Ti<:Integer, Tv<:Real}(tree::Ptr{Void}, name::Union(String,Nothing),
                    klass::Integer, flags::Integer, ind::Vector{Ti}, val::Vector{Tv},
                    constr_type::Integer, rhs::Float64)
    @check! _vectors_all_same_size(ind, val)
    ios_add_row(tree, name, klass, flags, length(ind), ind, val, constr_type, rhs)
end

ios_add_row{Ti<:Integer, Tv<:Real}(tree::Ptr{Void}, name::Union(String,Nothing),
           klass::Integer, ind::Vector{Ti}, val::Vector{Tv}, constr_type::Integer,
           rhs::Float64) =
    ios_add_row(tree, name, klass, 0, length(ind), ind, val, constr_type, rhs)

function ios_del_row(tree::Ptr{Void}, row::Integer)
    @check! _tree(tree)
    @check _reason(tree, [GLPK.ICUTGEN])
    @check _ios_row_is_valid(tree, row)
    @glpk_ccall ios_del_row Void (Ptr{Void}, Int32) tree row
end

function ios_clear_pool(tree::Ptr{Void})
    @check! _tree(tree)
    @check _reason(tree, [GLPK.ICUTGEN])
    @glpk_ccall ios_clear_pool Void (Ptr{Void},) tree
end

function init_env()
    ret = @glpk_ccall init_env Int32 ()
    @check! _init_env_succeeded(ret)
    return ret
end

function free_env()
    ret = @glpk_ccall free_env Int32 ()
end

function term_out(flag::Integer)
    @check _term_out_flag(flag)
    @glpk_ccall term_out Int32 (Int32,) flag
end

function open_tee(filename::String)
    ret = @glpk_ccall open_tee Int32 (Ptr{Uint8},) bytestring(filename)
    @check! _open_tee_succeeded(ret)
    return ret
end

function close_tee()
    @glpk_ccall close_tee Int32 ()
end

function malloc(size::Integer)
    @check _alloc_size(size)
    @glpk_ccall malloc Ptr{Void} (Int32,) size
end

function calloc(n::Integer, size::Integer)
    @check _alloc_size(n)
    @check _alloc_size(size)
    @glpk_ccall calloc Ptr{Void} (Int32, Int32) n size
end

function free(ptr::Ptr)
    @check _pointer_is_valid(ptr)
    @glpk_ccall free Void (Ptr{Void},) ptr
end

function mem_usage(count, cpeak, total, tpeak)
    error("unsupported. Use GLPK.mem_usage() instead.")
end

function mem_usage()
    data32 = Array(Int32, 2)
    data32_p = pointer(data32)
    off32 = sizeof(Int32)
    count_p = data32_p
    cpeak_p = data32_p + off32

    data64 = Array(Int64, 2)
    data64_p = pointer(data64)
    off64 = sizeof(Int64)

    total_p = data64_p
    tpeak_p = data64_p + off64

    @glpk_ccall mem_usage Void (Ptr{Int32}, Ptr{Int32}, Ptr{Int64}, Ptr{Int64}) count_p cpeak_p total_p tpeak_p 

    count = data32[1]
    cpeak = data32[2]
    total = data64[1]
    tpeak = data64[2]

    return count, cpeak, total, tpeak
end

function mem_limit(limit::Integer)
    @glpk_ccall mem_limit Void (Int32,) limit
end

function time()
    @glpk_ccall time Int64 ()
end

function difftime(t1::Integer, t0::Integer)
    @glpk_ccall difftime Float64 (Int64, Int64) t1 t0
end

function sdf_open_file(filename::String)
    @check _file_is_readable(filename)
    data_p = @glpk_ccall sdf_open_file Ptr{Void} (Ptr{Uint8},) bytestring(filename)
    @check! _sdf_file_opened(data_p)
    return Data(data_p)
end

function sdf_read_int(data::Data)
    @check! _data(data)
    @glpk_ccall sdf_read_int Int32 (Ptr{Void},) pointer(data)
end

function sdf_read_num(data::Data)
    @check! _data(data)
    @glpk_ccall sdf_read_num Float64 (Ptr{Void},) pointer(data)
end

function sdf_read_item(data::Data)
    @check! _data(data)
    item_cstr = @glpk_ccall sdf_read_item Ptr{Uint8} (Ptr{Void},) pointer(data)
    return bytestring(item_cstr)
end

function sdf_read_text(data::Data)
    @check! _data(data)
    text_cstr = @glpk_ccall sdf_read_text Ptr{Uint8} (Ptr{Void},) pointer(data)
    return bytestring(text_cstr)
end

function sdf_line(data::Data)
    @check! _data(data)
    @glpk_ccall sdf_line Int32 (Ptr{Void},) pointer(data)
end

function sdf_close_file(data::Data)
    @check! _data(data)
    @glpk_ccall sdf_close_file Void (Ptr{Void},) pointer(data)
end

function read_cnfsat(prob::Prob, filename::String)
    @check! _prob(prob)
    @check _file_is_readable(filename)
    ret = @glpk_ccall read_cnfsat Int32 (Ptr{Void}, Ptr{Uint8}) prob.p bytestring(filename)
    @check! _succeeded(ret, "read_cnfsat")
    return ret
end

function check_cnfsat(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall check_cnfsat Int32 (Ptr{Void},) prob.p
end

function write_cnfsat(prob::Prob, filename::String)
    @check! _prob(prob)
    @check _file_is_writable(filename)
    ret = @glpk_ccall write_cnfsat Int32 (Ptr{Void}, Ptr{Uint8}) prob.p bytestring(filename)
    @check! _succeeded(ret, "write_cnfsat")
    return ret
end

function minisat1(prob::Prob)
    @check! _prob(prob)
    @glpk_ccall minisat1 Int32 (Ptr{Void},) prob.p
end

function intfeas1(prob::Prob, use_bound::Integer, obj_bound::Integer)
    @check! _prob(prob)
    # TODO : more checks:
    #   1) columns must be GLPK.BV od GLPK.FX
    #   2) constraints and objj coeffs must be integer
    @glpk_ccall intfeas1 Int32 (Ptr{Void}, Int32, Int32) prob.p use_bound obj_bound
end


# FUNCTIONS NOT WRAPPED:
#
# 1) Printout functions [because they use varargs/va_list,
#    or need callbacks, or are implemented as macros]:
#
#    glp_printf
#    glp_vprintf
#    glp_term_hook
#    glp_error
#    glp_assert
#    glp_error_hook
#
# 2) Plain data file reading functions [because they either
#    use longjmp or varargs]:
#
#    glp_sdf_set_jump
#    glp_sdf_error
#    glp_sdf_warning
#
# 3) Additional functions [because wat?]:
#
#    lpx_check_kkt

#}}}

end # module
