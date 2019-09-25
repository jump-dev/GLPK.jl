## Check functions for internal use

# Functions which perform all sorts of sanity checks on input parameters and
# throw exceptions in case of errors. Ideally, it should never be possible to
# pass an invalid parameter to the underlying GLPK API when preemptive checks
# are active.
#
# All functions start with an underscore and are supposed to be called in one
# of these 2 forms:
#
#  @check _func
#  @check! _func
#
# where the first one is used for tests which can be disabled via
# `jl_set_preemptive_check(false)`, since they will be performed internally by
# GLPK anyway (but will likely produce a fatal error if they fail, while Julia
# handles things more gracefully).

if VERSION < v"0.7-"
    objectid(x) = object_id(x)
end

let valid_objs = Dict{UInt, Bool}()
    global jl_obj_is_valid, _add_obj, _del_obj, _del_all_objs
    function jl_obj_is_valid(x)
        haskey(valid_objs, objectid(x))
    end
    function _add_obj(x)
        valid_objs[objectid(x)] = true
    end
    function _del_obj(x)
        delete!(valid_objs, objectid(x))
    end
    function _del_all_objs()
        empty!(valid_objs)
    end
end

let PREEMPTIVE_CHECK = true
    global jl_get_preemptive_check,
           jl_set_preemptive_check
    function jl_get_preemptive_check()
        return PREEMPTIVE_CHECK
    end
    function jl_set_preemptive_check(pc::Bool)
        PREEMPTIVE_CHECK = pc
    end
end

macro check(expr)
    quote
        if jl_get_preemptive_check()
            $(esc(expr))
        else
            true
        end
    end
end

macro check!(expr)
    quote
        $(esc(expr))
    end
end

function _prob(prob::Prob)
    if prob.p == C_NULL || !jl_obj_is_valid(prob)
        throw(GLPKError("invalid GLPK.Prob"))
    end
    return true
end

function _string_length(s::AbstractString, maxl::Integer)
    l = length(s)
    if l > maxl
        throw(GLPKError("invalid string length $l (must be <= $maxl)"))
    end
    return true
end

function _row_is_valid(prob::Prob, row::Integer)
    rows = @glpk_ccall get_num_rows Cint (Ptr{Cvoid},) prob.p
    if (row < 1 || row > rows)
        throw(GLPKError("invalid row $row (must be 1 <= row <= $rows)"))
    end
    return true
end

function _col_is_valid(prob::Prob, col::Integer)
    cols = @glpk_ccall get_num_cols Cint (Ptr{Cvoid},) prob.p
    if (col < 1 || col > cols)
        throw(GLPKError("invalid col $col (must be 1 <= col <= $cols)"))
    end
    return true
end

function _col_is_valid_w0(prob::Prob, col::Integer)
    cols = @glpk_ccall get_num_cols Cint (Ptr{Cvoid},) prob.p
    if (col < 0 || col > cols)
        throw(GLPKError("invalid col $col (must be 0 <= col <= $cols)"))
    end
    return true
end

function _obj_dir_is_valid(dir::Integer)
    if !(dir == MIN || dir == MAX)
        throw(GLPKError("invalid obj_dir $dir (use MIN or MAX)"))
    end
    return true
end

function _bounds_type_is_valid(bounds_type::Integer)
    if !(bounds_type == FR ||
         bounds_type == LO ||
         bounds_type == UP ||
         bounds_type == DB ||
         bounds_type == FX)
        throw(GLPKError("invalid bounds_type $bounds_type (allowed values: GLPK.FR, GLPK.LO, GLPK.UP, GLPK.DB, GLPK.FX)"))
    end
    return true
end

function _bounds_are_valid(bounds_type::Integer, lb::Real, ub::Real, name::String)
    if bounds_type == DB && lb > ub
        throw(GLPKError("invalid bounds for double-bounded $name: $lb > $ub"))
    elseif bounds_type == FX && lb != ub
        throw(GLPKError("invalid bounds for fixed $name: $lb != $ub"))
    end
    return true
end

function _vectors_size(numel::Integer, vecs...)
    if numel < 0
        throw(GLPKError("invalid numer of elements: $numel"))
    end
    if numel > 0
        for v = vecs
            if isempty(v)
                throw(GLPKError("number of elements is $numel but vector is empty or nothing"))
            elseif length(v) < numel
                throw(GLPKError("wrong vector size: $(length(v)) (numel declared as $numel)"))
            end
        end
    end
    return true
end

function _vectors_all_same_size(vec0::VecOrNothing, vecs::VecOrNothing...)
    l0 = vecornothing_length(vec0)
    for v in vecs
        l = vecornothing_length(v)
        if l != l0
            throw(GLPKError("incosistent vector lengths: $l0 and $l"))
        end
    end
    return true
end

function _indices_vectors_dup(prob::Prob, numel::Integer, ia::Vector{Cint}, ja::Vector{Cint})
    rows = @glpk_ccall get_num_rows Cint (Ptr{Cvoid},) prob.p
    cols = @glpk_ccall get_num_cols Cint (Ptr{Cvoid},) prob.p
    #numel = length(ia)

    off32 = sizeof(Cint)
    iap = pointer(ia) - off32
    jap = pointer(ja) - off32

    k = @glpk_ccall check_dup Cint (Cint, Cint, Cint, Ptr{Cint}, Ptr{Cint}) rows cols numel iap jap
    if k < 0
        throw(GLPKError("indices out of bounds: $(ia[-k]),$(ja[-k]) (bounds are (1,1) <= (ia,ja) <= ($rows,$cols))"))
    elseif k > 0
        throw(GLPKError("duplicate index entry: $(ia[k]),$(ja[k])"))
    end
    return true
end

function _rows_and_cols(rows::Integer, cols::Integer)
    if rows < 0
        throw(GLPKError("rows < 0 : $rows"))
    end
    if cols < 0
        throw(GLPKError("cols < 0 : $rows"))
    end
end

function _rows_ids_size(prob::Prob, min_size::Integer, num_rows::Integer, rows_ids::Vector{Cint})
    rows = @glpk_ccall get_num_rows Cint (Ptr{Cvoid},) prob.p
    if num_rows < min_size || num_rows > rows
        throw(GLPKError("invalid vector size: $num_rows (min=$min_size max=$rows)"))
    end
    if num_rows == 0
        return true
    end
    if length(rows_ids) < num_rows
        throw(GLPKError("invalid vector size: declared>=$num_rows actual=$(length(rows_ids))"))
    end
    return true
end

function _rows_ids_content(prob::Prob, num_rows::Integer, rows_ids::Vector{Cint})
    rows = @glpk_ccall get_num_rows Cint (Ptr{Cvoid},) prob.p
    ind_set = BitSet()
    union!(ind_set, rows_ids[1 : num_rows])
    if minimum(ind_set) < 1 || maximum(ind_set) > rows
        throw(GLPKError("index out of bounds (min=1 max=$rows)"))
    elseif length(ind_set) != num_rows
        throw(GLPKError("one or more duplicate index(es) found"))
    end
    return true
end

function _cols_ids_size(prob::Prob, min_size::Integer, num_cols::Integer, cols_ids::Vector{Cint})
    cols = @glpk_ccall get_num_cols Cint (Ptr{Cvoid},) prob.p
    if num_cols < min_size || num_cols > cols
        throw(GLPKError("invalid vector size: $num_cols (min=$min_size max=$cols)"))
    end
    if num_cols == 0
        return 0
    end
    if length(cols_ids) < num_cols
        throw(GLPKError("invalid vector size: declared>=$num_cols actual=$(length(cols_ids))"))
    end
    return true
end

function _cols_ids_content(prob::Prob, num_cols::Integer, cols_ids::Vector{Cint})
    cols = @glpk_ccall get_num_cols Cint (Ptr{Cvoid},) prob.p
    ind_set = BitSet()
    union!(ind_set, cols_ids[1 : num_cols])
    if minimum(ind_set) < 1 || maximum(ind_set) > cols
        throw(GLPKError("index out of bounds (min=1 max=$cols)"))
    elseif length(ind_set) != num_cols
        throw(GLPKError("one or more duplicate index(es) found"))
    end
    return true
end

function _list_ids(prob::Prob, len::Integer, list_ids::Vector{Cint})
    if len == 0
        return true
    end
    rows = @glpk_ccall get_num_rows Cint (Ptr{Cvoid},) prob.p
    cols = @glpk_ccall get_num_cols Cint (Ptr{Cvoid},) prob.p
    # note1 the documentation does not mention forbidding duplicates in this case
    # note2 the size should already be checked as this function is only called
    #       by GLPK.print_ranges
    #if len < 0 #|| len > rows + cols
    ##throw(GLPKError("invalid vector size: $len (min=0 max=$(rows + cols))"))
    #throw(GLPKError("invalid vector size: $len < 0"))
    #end
    #if length(list_ids) < len
    #throw(GLPKError("invalid vector size: declared>=$len actual=$(length(list_ids))"))
    #end
    if minimum(list_ids[1:len]) < 1 || maximum(list_ids[1:len]) > rows + cols
        throw(GLPKError("index out of bounds (min=1 max=$(rows + cols))"))
    end
    return true
end

function _status_is_optimal(prob::Prob)
    ret = @glpk_ccall get_status Cint (Ptr{Cvoid},) prob.p
    if ret == OPT
        throw(GLPKError("current basic solution is not optimal"))
    end
    return true
end

function _bf_exists(prob::Prob)
    ret = @glpk_ccall bf_exists Cint (Ptr{Cvoid},) prob.p
    if ret == 0
        throw(GLPKError("no bf solution found (use GLPK.factorize)"))
    end
    return true
end

function _var_is_basic(prob::Prob, ind::Integer)
    rows = @glpk_ccall get_num_rows Cint (Ptr{Cvoid},) prob.p
    if ind <= rows
        j = @glpk_ccall get_row_stat Cint (Ptr{Cvoid}, Cint) prob.p ind
        if j != BS
            throw(GLPKError("variable $ind is non-basic"))
        end
    else
        j = @glpk_ccall get_col_stat Cint (Ptr{Cvoid}, Cint) prob.p ind-rows
        if j != BS
            throw(GLPKError("variable $ind is non-basic"))
        end
    end
end

function _var_is_non_basic(prob::Prob, ind::Integer)
    rows = @glpk_ccall get_num_rows Cint (Ptr{Cvoid},) prob.p
    if ind <= rows
        j = @glpk_ccall get_row_stat Cint (Ptr{Cvoid}, Cint) prob.p ind
        if j == BS
            throw(GLPKError("variable $ind is basic"))
        end
    else
        j = @glpk_ccall get_col_stat Cint (Ptr{Cvoid}, Cint) prob.p ind-rows
        if j == BS
            throw(GLPKError("variable $ind is basic"))
        end
    end
end

function _is_prim_feasible(prob::Prob)
    if FEAS != @glpk_ccall get_prim_stat Cint (Ptr{Cvoid},) prob.p
        throw(GLPKError("problem is not primal feasible"))
    end
    return true
end

function _is_dual_feasible(prob::Prob)
    if FEAS != @glpk_ccall get_dual_stat Cint (Ptr{Cvoid},) prob.p
        throw(GLPKError("problem is not dual feasible"))
    end
    return true
end

function _copy_names_flag(names::Integer)
    if names != ON && names != OFF
        throw(GLPKError("invalid copy_names flag $names (use GLPK.ON or GLPK.OFF)"))
    end
    return true
end

function _scale_flags(flags::Integer)
    all = (SF_GM | SF_EQ | SF_2N | SF_SKIP)
    if (flags | all) != all && flags != SF_AUTO
        throw(GLPKError("invalid scale flags $flags"))
    end
    return true
end

function _stat_is_valid(stat::Integer)
    if (stat != BS &&
        stat != NL &&
        stat != NU &&
        stat != NF &&
        stat != NS)
        throw(GLPKError("invalid status $stat (use GLPK.BS or GLPK.NL or GLPK.NU or GLPK.NF or GLPK.NS)"))
    end
end

function _adv_basis_flags(flags::Integer)
    if flags != 0
        throw(GLPKError("adv_basis flags must be set to 0 (found $flags instead)"))
    end
    return true
end

function _kind_is_valid(kind::Integer)
    if (kind != CV &&
        kind != IV &&
        kind != BV)
        throw(GLPKError("invalid kind $kind (use GLPK.CV or GLPK.IV or GLPK.BV)"))
    end
    return true
end

function _file_is_readable(filename::AbstractString)
    try
        f = open(filename, "r")
        close(f)
    catch err
        throw(GLPKError("file $filename not readable"))
    end
    return true
end

function _file_is_writable(filename::AbstractString)
    try
        f = open(filename, "w")
        close(f)
    catch err
        throw(GLPKError("file $filename not writable"))
    end
    return true
end

function _succeeded(ret::Integer, msg::AbstractString)
    if ret != 0
        throw(GLPKError(msg * " failed"))
    end
    return true
end

function _kkt_cond_param(cond::Integer, sol::Integer)
    if !(cond in [KKT_PE, KKT_PB, KKT_DE, KKT_DB])
        throw(GLPKError("invalid cond parameter $cond (use GLPK.KKT_PE, GLPK.KKT_PB, GLPK.KKT_DE or GLPK.KKT_DB)"))
    end
    if sol != IPT && cond in [KKT_DE, KKT_DB]
        throw(GLPKError("invalid cond parameter $cond used with sol parameter $sol (use GLPK.KKT_PE or GLPK.KKT_PB)"))
    end
    return true
end

function _mps_format(format::Integer)
    if (format != MPS_DECK &&
        format != MPS_FILE)
        throw(GLPKError("invalid MPS format $format (use GLPK.MPS_DECK or GLPK.MPS_FILE)"))
    end
    return true
end

function _mps_param(param)
    if param != C_NULL
        throw(GLPKError("MPS param must be C_NULL"))
    end
    return true
end

function _lp_param(param)
    if param != C_NULL
        throw(GLPKError("LP param must be C_NULL"))
    end
    return true
end

function _read_prob_flags(flags::Integer)
    if flags != 0
        throw(GLPKError("read_prob flags must be 0"))
    end
    return true
end

function _write_prob_flags(flags::Integer)
    if flags != 0
        throw(GLPKError("write_prob flags must be 0"))
    end
    return true
end

function _print_ranges_flags(flags::Integer)
    if flags != 0
        throw(GLPKError("print_ranges flags must be set to 0 (found $flags instead)"))
    end
    return true
end

function _mpl_workspace(tran::MathProgWorkspace)
    if tran.p == C_NULL || !jl_obj_is_valid(tran)
        throw(GLPKError("invalid GLPK.MathProgWorkspace"))
    end
    return true
end

function _sol_param(sol::Integer)
    if !(sol == SOL || sol == IPT || sol == MIP)
        throw(GLPKError("invalid parameter sol $sol (use GLPK.SOL, GLPK.IPT or GLPK.MIP)"))
    end
    return true
end

function _rowcol_is_valid(k::Integer, k_max::Integer)
    if !(1 <= k <= k_max)
        throw(GLPKError("index out of bounds: $k (bounds are 1 <= k <= $k_max"))
    end
    return true
end

function _rowcol_is_valid(prob::Prob, k::Integer)
    rows = @glpk_ccall get_num_rows Cint (Ptr{Cvoid},) prob.p
    cols = @glpk_ccall get_num_cols Cint (Ptr{Cvoid},) prob.p

    k_max = rows + cols

    _rowcol_is_valid(k, k_max)
end


function _dir_is_valid(dir::Integer)
    if !(dir == 1 || dir == -1)
        throw(GLPKError("invalid direction $dir (must be 1 or -1)"))
    end
    return true
end

function _eps_is_valid(eps::Real)
    if (eps < 0)
        throw(GLPKError("invalid eps $eps (must be >= 0)"))
    end
    return true
end

function _tree(tree::Ptr{Cvoid})
    if tree == C_NULL
        throw(GLPKError("invalid tree pointer"))
    end
    return true
end

function _reason(tree::Ptr{Cvoid}, allowed::Vector{Cint})
    reason = @glpk_ccall ios_reason Cint (Ptr{Cvoid},) tree
    if !(reason in allowed)
        throw(GLPKError("callback operation not allowed at current stage"))
    end
    return true
end

function _can_branch(tree::Ptr{Cvoid}, col::Integer)
    if 0 == @glpk_ccall ios_can_branch Cint (Ptr{Cvoid}, Cint) tree col
        throw(GLPKError("column $col cannot branch"))
    end
    return true
end

function _ios_node_is_valid(tree::Ptr{Cvoid}, node::Integer)
    valid_nodes = BitSet()
    n = 0
    while true
        n = @glpk_ccall ios_next_node Cint (Ptr{Cvoid}, Cint) tree n
        if n == 0
            break
        end
        if node == n
            return true
        end
        p = n
        while !(p in valid_nodes)
            push!(valid_nodes, p)
            p = @glpk_ccall ios_up_node Cint (Ptr{Cvoid}, Cint) tree p
            if p == 0
                break
            end
            if node == p
                return true
            end
        end
    end
    throw(GLPKError("invalid node $node"))
end

function _ios_node_is_active(tree::Ptr{Cvoid}, node::Integer)
    c = 0
    while node > c
        c = @glpk_ccall ios_next_node Cint (Ptr{Cvoid}, Cint) tree c
        if c == 0
            break
        end
    end
    if c != node
        throw(GLPKError("node $node is not in the active list"))
    end
    return true
end

function _sel_is_valid(sel::Integer)
    if !(sel == DN_BRNCH || sel == UP_BRNCH || sel == NO_BRNCH)
        throw(GLPKError("invalid select flag: $sel (allowed values: GLPK.DN_BRNCH, GLPK.UP_BRNCH, GLPK.NO_BRNCH)"))
    end
end

function _klass_is_valid(klass::Integer)
    if klass != 0 && !(101 <= klass <= 200)
        throw(GLPKError("invalis klass $klass (should be 0 or 101 <= klass <= 200)"))
    end
    return true
end

function _ios_add_row_flags(flags::Integer)
    if flags != 0
        throw(GLPKError("ios_add_row flags must be 0"))
    end
    return true
end

function _constr_type_is_valid(constr_type::Integer)
    if !(constr_type == LO || constr_type == UP)
        throw(GLPKError("invalid constr_type $constr_type (allowed values: GLPK.LO, GLPK.UP)"))
    end
    return true
end

function _ios_row_is_valid(tree::Ptr{Cvoid}, row::Integer)
    size = @glpk_ccall ios_pool_size Cint (Ptr{Cvoid},) tree
    if !(1 <= row <= size)
        throw(GLPKError("invalid ios row $row (must be 1 <= row <= $size)"))
    end
    return true
end

function _init_env_succeeded(ret::Integer)
    if !(0 <= ret <= 1)
        throw(GLPKError("initialization failed"))
    end
    return true
end

function _term_out_flag(flag::Integer)
    if !(flag == ON || flag == OFF)
        throw(GLPKError("invalid flag $flag (use GLPK.ON or GLPK.OFF)"))
    end
    return true
end

function _open_tee_succeeded(ret::Integer)
    if !(0 <= ret <= 1)
        throw(GLPKError("GLPK.open_tee failed"))
    end
    return true
end

function _alloc_size(n::Integer)
    if n <= 0
        throw(GLPKError("invalid alloc size $n"))
    end
    return true
end

function _pointer_is_valid(ptr::Ptr)
    if ptr == C_NULL
        throw(GLPKError("invalid pointer"))
    end
    return true
end
