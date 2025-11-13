module FancyArrows

import ..Makie
import ..Makie: BezierPath, Point2, Point2d, MoveTo, LineTo, CurveTo, ClosePath
import ..Makie: Observable, ComputeGraph, add_input!, add_constant!, update!, register_computation!, attributes
import ..Makie.ComputePipeline: alias!, get_observable!
import ..Makie: translate, scale, rotate 
import ..Makie: PlotSpec, PlotList, Poly, Lines
import ..Makie: to_color, miter_angle_to_distance, miter_distance_to_angle
import ..Makie: init_arrow_style, plotspecs, shrinksize
import LinearAlgebra: qr

#=
An `AbstractDecoration` object has a `rotation` that describes how its oriented 
with respect to the positive x-axis and a `connector` (2d vector like) 
where the axis line should attach to.
=#
"Abstract supertype for `LineAxis` decoration objects that are similar to plot objects."
abstract type AbstractDecoration end

"Return the current `connector` point."
connector_vector(::AbstractDecoration)=nothing::Vector{Float32}
"Set target rotation with respect to positive x-axis."
set_deco_rotation!(deco::AbstractDecoration, rot; kwargs...)=deco
"Set connector to match target end point."
match_deco_connector!(deco::AbstractDecoration, target)=deco
"Return a list of plotting instruction fed to `plotlist!(parent, ...)`."
deco_plotspecs(::AbstractDecoration; kwargs...) = PlotSpec[]

## compatibility with annotations
shrinksize(::AbstractDecoration) = 0f0
function plotspecs(deco::AbstractDecoration, pos; rotation, kwargs...)
    set_deco_rotation!(deco, rotation; check_reversed=true)
    match_deco_connector!(deco, pos)
    return deco_plotspecs(deco; space=:pixel, kwargs...)
end

"Helper type implementing the `AbstractDecoration` interface without actually producing plots."
struct NoDeco <: AbstractDecoration
    connector::Vector{Float32}
end

NoDeco() = NoDeco(zeros(Float32, 2))
connector_vector(deco::NoDeco) = deco.connector

function match_deco_connector!(deco::NoDeco, target)
    deco.connector .= target
    return deco
end

"""
    BezigonDeco(attributes :: ComputeGraph)

Decoration based on some `path::BezierPath` that is stored in `attributes`
and can be plotted stroked and/or filled."""
struct BezigonDeco{attrs_Type} <: AbstractDecoration
    attributes :: attrs_Type
end

attributes(deco::BezigonDeco)=deco.attributes
connector_vector(deco::BezigonDeco) = attributes(deco)[:connector][]::Vector{Float32}
shrinksize(deco::BezigonDeco) = attributes(deco)[:shrinksize][]

function set_deco_rotation!(deco::BezigonDeco, rotation; check_reversed=true)
    attrs = attributes(deco)
    if check_reversed
        rev = attrs[:reversed][]
        rotation = rev ? rotation - π : rotation
    end
    update!(attributes(deco); rotation)
    return deco
end

function match_deco_connector!(deco::BezigonDeco, target!)
    update!(attributes(deco); target!)
    return deco
end

function deco_plotspecs(deco::BezigonDeco; space=:data, make_static=true, kwargs...)
    graph = attributes(deco)
    ## NOTE `PlotSpec` currently does not like `Observable` args.
    ##      Neither are `Computed` kwargs supported.
    ##      With `annotation`, there is also issues if Observables are used.
    ##      Hence, we read the values and return static specs.
    ##      In context, we have to ensure that everything is up-to-date when calling `deco_plotspecs`.
    sval = make_static ? Makie.to_value : identity
    stroke_plot = PlotSpec(
        Lines, sval(graph[:path]); 
        visible = sval(graph[:stroke_visible]), 
        linewidth = sval(graph[:strokewidth_fin]),
        color = sval(graph[:strokecolor_fin]), 
        linecap = sval(graph[:linecap]), 
        joinstyle = sval(graph[:joinstyle]),
        miter_limit = sval(get(graph, :miter_limit_fin, graph[:miter_limit])),
        space
    )
    fill_plot = PlotSpec(
        Poly, sval(graph[:path]);
        visible = sval(graph[:fill_visible]), 
        strokewidth = 0, 
        color = sval(graph[:fillcolor_fin]),
        space
    )
    specs = PlotSpec[fill_plot, stroke_plot]
    return specs
end

abstract type AbstractDecoSpecification end

## compatibility with annotations
struct NoDecoSpec <: AbstractDecoSpecification end
convert_to_deco_spec(spec::AbstractDecoSpecification)=spec
convert_to_deco_spec(spec)=NoDecoSpec()
init_arrow_style(::NoDecoSpec; kwargs...) = NoDeco()

abstract type AbstractBezigonDecoSpecification <: AbstractDecoSpecification end

function init_arrow_style(spec::T; kwargs...) where T<:AbstractBezigonDecoSpecification
    graph = _common_bezigon_setup(spec; kwargs...)
    graph = _register_path_computations!(graph, T)
    return BezigonDeco(graph)
end

struct ComputerModernRightarrowSpecification{attrs_Type} <: AbstractBezigonDecoSpecification
    attrs :: attrs_Type
end

function ComputerModernRightarrowSpecification(;
    length = (1.9f0, 2.2f0),
    width_ = (0f0, 2.096744f0),
    line_width = (0f0, 1f0),
    linecap = :round,
    joinstyle = :round,
    is_filled = false,
    kwargs...
    )
    attrs = _setup_base_bezigon_attributes(;
        length, width_, line_width, linecap, joinstyle, is_filled, kwargs...)
    return ComputerModernRightarrowSpecification(attrs)
end

function ComputerModernRightarrowTip(; kwargs...)
    ComputerModernRightarrowSpecification(; reversed=false, kwargs...)
end
function ComputerModernRightarrowTail(; kwargs...)
    ComputerModernRightarrowSpecification(; reversed=true, kwargs...)
end

function _register_path_computations!(graph, ::Type{<:ComputerModernRightarrowSpecification})
    graph = _register_appearance_attributes!(graph)
    cx1 = -.81731f0; cy1 = .2f0
    cx2 = -.41019f0; cy2 = .05833333f0
    map!(
        graph, 
        [:arrow_length, :arrow_width, :strokewidth_fin, :stroke_bgon, :reversed, :joinstyle, :miter_limit], 
        [:inner_length, :inner_width, :line_end, :back_end, :tip_end, :miter_limit_fin]
    ) do arr_length, arr_width, _sw, stroke_bgon, rev, join, mlimit
        sw = stroke_bgon ? _sw : 0f0
        inner_length = arr_length - sw
        inner_width = arr_width - sw

        tan_psi_tip = cy2 * inner_width / ((1-cx2) * inner_length)
        sin_psi_tip_inv = sqrt( 1/tan_psi_tip^2 + 1)
        miter_half_len = sin_psi_tip_inv * sw / 2
        tip_end = if join == :round
            sw / 2
        elseif join == :miter
            miter_half_len
        else# if join == :bevel
            1 / sin_psi_tip_inv * sw
        end
        miter_limit = if join == :miter
            _mlimit = miter_distance_to_angle(miter_half_len / sw) / 2
            max(mlimit, _mlimit)
        else
            mlimit
        end
        back_end = inner_length - sw / 2
        line_end = - sw / 2
        return (inner_length, inner_width, line_end, back_end, tip_end, miter_limit)
    end
    
    map!(graph, [:inner_length, :inner_width], :bezier_points) do L, W
        p0 = Point2(-L, W/2)
        p1 = Point2(cx1 * L, cy1 * W)
        p2 = Point2(cx2 * L, cy2 * W)
        p3 = Point2(0, 0)
        return (p0, p1, p2, p3)
    end
    
    map!(
        graph, :bezier_points, :path0
    ) do (P0, P1, P2, P3) 
        swap(p) = Point2d(p[1], -p[2])
        _P2 = swap(P2)
        _P1 = swap(P1)
        _P0 = swap(P0)
        BezierPath([
            MoveTo(P0),             # top left
            CurveTo(P1,P2,P3),      # tip
            CurveTo(_P2,_P1,_P0)    # bottom left
        ])
    end
    
    graph = _register_bezigon_drawing_cmds!(graph)

    return graph
end

struct LatexSpecification{attrs_Type} <: AbstractBezigonDecoSpecification 
    attrs :: attrs_Type
end

function LatexSpecification(;
    length = (3f0, 4.8f0),
    width_ = (0f0, 0.75f0),
    line_width = (0f0, 1f0),
    kwargs... 
)
    attrs = _setup_base_bezigon_attributes(;
        length, width_, line_width, kwargs...)
    return LatexSpecification(attrs)
end

function LatexTip(; kwargs...)
    return LatexSpecification(; reversed = false, kwargs...)
end

function LatexTail(; kwargs...)
    return LatexSpecification(; reversed = true, kwargs...)
end

const LATEX_ARROW_BEZIER_CONSTANTS = (;
    cx1 = .877192f0, cy1 = .077922f0,
    cx2 = .337381f0, cy2 = .519480f0,
)

function _register_path_computations!(graph, ::Type{<:LatexSpecification})
    global LATEX_ARROW_BEZIER_CONSTANTS
    ## cap actual strokewidth to be at most fifth of arrow length    
    map!(graph, [:arrow_length, :strokewidth], :strokewidth_fin) do arr_length, sw
        min(sw, arr_length/5)
    end
    
    graph = _register_appearance_attributes!(graph)
    
    cx1, cy1, cx2, cy2 = LATEX_ARROW_BEZIER_CONSTANTS
    map!(
        graph, 
        [:arrow_length, :arrow_width, :strokewidth_parent, :strokewidth_fin, :stroke_bgon, :reversed, :joinstyle, :miter_limit], 
        [:line_end, :back_end, :tip_end, :bezier_points, :miter_limit_fin]
    ) do arr_length, arr_width, sw_par, _sw, stroke_bgon, rev, join, mlimit
        sw = stroke_bgon ? _sw : 0f0
        arr_halfwidth = arr_width / 2
        if sw > 0   ## if there is no stroke, we don't need these computations
            ## compute miter length
            ## formula: `miter_len = sw / sin(ψ)`
            ##          `ψ` is the angle between bottom leg and hypotenuse of right-angled triangle
            ##          assume bottom leg has length `λ` and other leg has length `ω`
            ##          Then `1 / sin(ψ) = sqrt( λ^2 / ω^2 + 1)`.
            ##          Looking at the drawing code the tangents at the tip form a suitable triangle
            ##          with `λ = (1 - cx1) * arr_length` (or rather `inner_length`) and `ω = cy1 * arr_width / 2`
            tan_psi_tip = cy1 * arr_halfwidth / ((1-cx1) * arr_length)
            psi_tip = atan(tan_psi_tip)     # half angle at tip
            sin_psi_tip = sin(psi_tip)      # 1/sqrt(9 * L^2 / H^2 + 1) in LaTeX
            miter_half_len = ((1 / sin_psi_tip) * sw) / 2

            harpoon_extra_len = (sw / 2) / tan_psi_tip
            ## for inner length, substract half miter_len (front), and half strokewidth (back)
            back_end = -sw / 2
            ## we want non-trivial inner_length, i.e., arr_length + back_end - miter_half_len >= inner_length_min
            inner_length = arr_length + back_end
            inner_length_min = eps(Float32)
            if inner_length - miter_half_len < inner_length_min
                miter_half_len = inner_length - inner_length_min
            end

            ## (vertical) back miter
            ## φ = 2 * ψ
            tan_phi_tail = (arr_length * cx2) / ((1 - cy2) * arr_halfwidth)
            psi_tail = atan(tan_phi_tail) / 2
            sin_psi_tail = sin(psi_tail)
            bmiter_half_len = ((1 / sin_psi_tail) * sw) / 2
            inner_halfwidth = arr_halfwidth
            inner_halfwidth_min = inner_length_min
            if inner_halfwidth - sqrt(bmiter_half_len^2 - (sw/2)^2) < inner_halfwidth_min
                bmiter_half_len = sqrt((inner_halfwidth - inner_halfwidth_min)^2 + (sw/2)^2)
            end

            ## with Cairo, there are issues with very small `miter_limit` (angle) values
            ## we compute the miter limit that would give `miter_half_len` and use it instead
            if join == :miter
                _mlimit = min(
                    miter_distance_to_angle(miter_half_len / sw),
                    miter_distance_to_angle(bmiter_half_len / sw),
                ) / 2   ## divide by 2 to prevent rounding errors
                mlimit, _mlimit
                miter_limit = max(mlimit, _mlimit)
            else
                miter_limit = mlimit
            end
            miter_half_len_max = miter_angle_to_distance(miter_limit) * sw
            miter_half_len = min(miter_half_len, miter_half_len_max)
            bmiter_half_len = min(bmiter_half_len, miter_half_len_max)

            inner_length -= miter_half_len
            inner_halfwidth -= sqrt(bmiter_half_len^2 - (sw/2)^2)

            line_end = rev ? inner_length - sw_par / 2 : 0f0 |> f32
            tip_end = if join == :round
                inner_length + sw / 2
            elseif join == :miter
                inner_length + miter_half_len
            else# join == :bevel
                inner_length + sin_psi_tip * sw
            end
            
            ## modify control points to match scaled geometry
            bezier_points = _latex_bezier_points(arr_length, arr_halfwidth, inner_length, inner_halfwidth)
        else
            line_end = rev ? arr_length - sw_par / 2 : 0f0
            tip_end = arr_length
            back_end = 0f0
            bezier_points = _latex_default_bezier_points(arr_length, arr_halfwidth)
            miter_limit = mlimit
        end
        return (line_end, back_end, tip_end, bezier_points, miter_limit)
    end
    
    map!(
        graph, 
        :bezier_points,
        :path0
    ) do (P0, P1, P2, P3) 
        swap(p) = Point2d(p[1], -p[2])
        _P3 = swap(P3)
        _P2 = swap(P2)
        _P1 = swap(P1)
        BezierPath([
            MoveTo(P0),             # tip
            CurveTo(P1,P2,P3),      # top left
            LineTo(_P3),            # bottom left
            CurveTo(_P2,_P1,P0)     # tip
        ])
    end
    
    graph = _register_bezigon_drawing_cmds!(graph)

    return graph
end

function _latex_default_bezier_points(arrow_length, arrow_halfwidth)
    global LATEX_ARROW_BEZIER_CONSTANTS
    cx1, cy1, cx2, cy2 = LATEX_ARROW_BEZIER_CONSTANTS
    l = arrow_length
    h = arrow_halfwidth
    p0 = [l, 0f0]
    p1 = [cx1 * l, cy1 * h]
    p2 = [cx2 * l, cy2 * h]
    p3 = [0, h]
    return [Point2(p0), Point2(p1), Point2(p2), Point2(p3)]
 end

function _latex_bezier_points(
    arrow_length, arrow_halfwidth,
    inner_length, inner_halfwidth;
    kwargs...,
)
    p0, p1, p2, p3 = _latex_default_bezier_points(arrow_length, arrow_halfwidth)   
    
    L = inner_length
    H = inner_halfwidth
    P0 = [L; 0f0]
    P3 = [0f0; H]
    return _match_bezier_curveto_tangents(
        p0, p1, p2, p3, P0, P3; kwargs...
    )
end

function _match_bezier_curveto_tangents(
    # target curve from p0 to p3 with controls p1 and p2
    p0, p1, p2, p3,
    # new endpoints
    P0, P3;
    # curve parameter values where tangents should match
    eq = [0f0, 1f0], ls = Float32[1/6, 2/6, 3/6, 4/6, 5/6]
)
    dt_b(t) = (1 - t)^2 .* (p1 .- p0) .+ 2 * (1-t) * t .* (p2 .- p1) .+ t^2 .* (p3 .- p2)   # * 3

    # (1 - t)^2 * (x1 - P0[1]) + 2 * (1-t) * t * (x2 - x1) +  t^2 * (P3[1] - x2)
    # (1-t)^2*x1 - (1-t)^2*P0[1] + 2*(1-t)*t*x2 - 2*(1-t)*t*x1 + t^2*P3[1] - t^2*x2
    # ((1-t)^2 - 2*(1-t)*t)*x1 + (2*(1-t)*t - t^2)*x2 + (t^2*P3[1] - (1-t)^2*P0[1])
    # (similar for y1 & y2)
    dt_B1(t) = (1-t)^2 - 2*(1-t)*t
    dt_B2(t) = 2*(1-t)*t - t^2
    dt_Bx(t) = t^2*P3[1] - (1-t)^2*P0[1]
    dt_By(t) = t^2*P3[2] - (1-t)^2*P0[2]
    at = vcat(eq, ls)
    N = length(at)
    A = zeros(Float32, N, 2)
    b = zeros(Float32, N, 2)
    for (i, t) in enumerate(at)
        A[i, 1] = dt_B1(t) 
        A[i, 2] = dt_B2(t) 
        b[i, :] .= dt_b(t)
        b[i, 1] -= dt_Bx(t)
        b[i, 2] -= dt_By(t)
    end
    if length(eq) < 2
        β = A \ b
    else
        # https://en.wikipedia.org/wiki/Ordinary_least_squares#Constrained_estimation
        Qt = A[1:length(eq), :]
        Q = transpose(Qt)
        c = b[1:2, :]
        X = A
        y = b
        XtX = qr(X'X)
        Xy = X'y
        α = XtX \ Xy
        tmp0 = Qt * α - c
        tmp1 = Qt * (XtX \ Q)
        tmp2 = Q * (tmp1 \ tmp0)
        tmp3 = XtX \ tmp2
        β = α - tmp3
    end 
    P1 = β[1, :]
    P2 = β[2, :]
    return [Point2(P0), Point2(P1), Point2(P2), Point2(P3)]
end

#=
function _latex_bezier_points(
    arrow_length, arrow_halfwidth,
    inner_length, inner_halfwidth, _sw;
    u1 = 1/3, u2 = 2/3, max_iter=50, tol=sqrt(eps(Float32))
)
    # https://perso.liris.cnrs.fr/victor.ostromoukhov/publications/pdf/ICCG93_Hermite.pdf
    global LATEX_ARROW_BEZIER_CONSTANTS
    cx1, cy1, cx2, cy2 = LATEX_ARROW_BEZIER_CONSTANTS
    sw = _sw/2
    l = arrow_length
    h = arrow_halfwidth
    p0 = [l, 0f0]
    p1 = [cx1 * l, cy1 * h]
    p2 = [cx2 * l, cy2 * h]
    p3 = [0, h]
    
    b(t) = (1-t)^3 .* p0 .+ 3*(1-t)^2*t .* p1 .+ 3*(1-t)*t^2 .* p2 .+ t^3 .* p3             # generator
    dt_b(t) = 3 .* ((1-t)^2 .* (p1 .- p0) .+ 2*(1-t)*t .* (p2 .- p1) .+ t^2 .* (p3 .- p2))  # tangent
    nb(t) = let bt=dt_b(t); n = [-bt[2]; bt[1]]; n ./ sqrt(sum(abs2, n)) end                # normal
   
    @show dt_b(u1), nb(u1)
    @show Gd1 = b(u1) .+ sw .* nb(u1)
    @show Gd2 = b(u2) .+ sw .* nb(u2)
    
    L = inner_length
    H = inner_halfwidth
    P0 = [L; 0f0]
    P3 = [0f0; H]

    @show G0 = dt_b(0)
    @show G1 = dt_b(1)
    ## tangent at t=0 is 3 * (P1 - P0), should match 3 * λ0 * G0, i.e. P1 = λ0 .* G0 .+ P0
    ## tangent at t=1 is 3 * (P3 - P2), should match 3 * λ1 * G1, i.e. P2 = P3 .- λ1 .* G1

    #B(t, λ0, λ1) = (1-t)^3 .* P0 .+ 3*(1-t)^2*t .* (λ0 .* G0 .+ P0) .+ 3*(1-t)*t^2 .* (P3 - λ1 .* G1) .+ t^3 .* P3 # offset
    a0(t) = 3 * (1-t)^2 * t .* G0
    a1(t) = -3 * (1-t) * t^2 .* G1
    rhs(t) = ((1-t)^3 + 3*(1-t)^2*t) .* P0 .+ (t^3 + 3*(1-t)*t^2) .* P3
    B(t, λ0, λ1) = λ0 .* a0(t) .+ λ1 .* a1(t) .+ rhs(t)

    c = (sqrt(l^2 + h^2) - sqrt(L^2 + H^2)) / 2
    t1 = u1 - c
    t2 = u2 + c
 
    A = hcat(a0(t1), a1(t1))
    r = Gd1 .- rhs(t1)
    Λ = _Λ1 = A \ r

    A = hcat(a0(t2), a1(t2))
    r = Gd2 .- rhs(t2)
    Λ = _Λ2 = A \ r

    _t1 = t1
    _t2 = t2

    t1 = t1 + 1/20
    t2 = t2 - 1/20
    
    err1 = abs(t1 - _t1)
    err2 = abs(t2 - _t2)
    for i = 1:max_iter

        if err1 != 0
            A = hcat(a0(t1), a1(t1))
            r = Gd1 .- rhs(t1)
            Λ = A \ r
            @show A * Λ .- r
            @show B(t1, Λ[1], Λ[2]) .- Gd1

            _Q1 = B(t2, _Λ1[1], _Λ1[2])
            Q1 = B(t2, Λ[1], Λ[2])
            ## on line from `_Q1` to `Q1`, find nearest point to `Gd2`
            v = Q1 .- _Q1
            u = _Q1 .- Gd2
            tmid = min(1f0, max(0f0, -v'u / v'v))
            #tmid = -v'u / v'v
            __t1 = t1
            t1 = tmid * (t1 - _t1) + _t1
            _t1 = __t1
            _Λ1 = Λ
            @show err1 = abs(t1 - _t1)
        end
        if err2 != 0
            A = hcat(a0(t2), a1(t2))
            r = Gd2 .- rhs(t2)
            @show Λ = A \ r
            @show B(t2, Λ[1], Λ[2]) .- Gd2
            
            _Q2 = B(t1, _Λ2[1], _Λ2[2])
            Q2 = B(t1, Λ[1], Λ[2])
            v = Q2 .- _Q2
            u = _Q2 .- Gd1
            tmid = min(1f0, max(0f0, -v'u / v'v))
            #tmid = -v'u / v'v
            __t2 = t2
            t2 = tmid * (t2 - _t2) + _t2
            _t2 = __t2
            _Λ2 = Λ
            @show err2 = abs(t2 - _t2)
        end

        @show err1 + err2
        if err1 + err2 <= tol
            break
        end        
    end
    @show t1, t2, Λ
    P1 = P0 .+ Λ[1] .* G0
    P2 = P3 .- Λ[2] .* G1

    bpath(t, p0, p1, p2, p3) = (1-t)^3 .* p0 .+ 3*(1-t)^2*t .* p1 .+ 3*(1-t)*t^2 .* p2 .+ t^3 .* p3
    dt_bpath(t, p0, p1, p2, p3) = 3 .* ((1-t)^2 .* (p1 .- p0) .+ 2*(1-t)*t .* (p2 .- p1) .+ t^2 .* (p3 .- p2))
    @show p0, p1, p2, p3
    @show P0, P1, P2, P3
    @show dt_bpath(0, P0, P1, P2, P3) ./ dt_bpath(0, p0, p1, p2, p3)
    @show bpath(t1, P0, P1, P2, P3) .- Gd1
    @show bpath(t2, P0, P1, P2, P3) .- Gd2
    return [Point2d(p) for p in [P0, P1, P2, P3]]
end
=#

function _setup_base_bezigon_attributes(;
    ## size attributes,
    length = nothing,
    width = nothing,
    width_ = nothing,
    angle = nothing,
    angle_ = nothing,
    line_width = nothing,
    line_width_ = nothing,
    ## geometric appearance
    miter_limit::Number=eps(Float32),
    reversed::Bool=false,
    swapped::Bool=false,
    scale::Number=1f0,
    align::Number=0f0,
    linecap::Symbol=:butt,
    joinstyle::Symbol=:miter,
    ## coloring
    color=nothing,
    fillcolor=nothing,
    strokecolor=nothing,
    is_filled::Bool=true,
    is_stroked::Bool=true,
    kwargs...
)
    graph = ComputeGraph()
    ## add size attributes and automatically convert to tuples
    add_input!(make_spec2, graph, :length, length)
    add_input!(make_spec2, graph, :width, width)
    add_input!(make_spec2, graph, :width_, width_)
    add_input!(make_spec3, graph, :angle, angle)
    add_input!(make_spec1, graph, :angle_, angle_)
    add_input!(make_spec2, graph, :line_width, line_width)
    add_input!(make_spec2, graph, :line_width_, line_width_)

    ## geometric attributes
    add_input!((k, t) -> min(1f0, max(0f0, f32(t))), graph, :align, align)
    
    miter_lim_min = eps(Float32) 
    miter_lim_max = π/2 - miter_lim_min
    add_input!((k,v) -> min(miter_lim_max, max(miter_lim_min, f32(v))), graph, :miter_limit, miter_limit)
    
    add_input!(f32, graph, :scale, scale)
    add_input!(boolean, graph, :reversed, reversed)
    add_input!(boolean, graph, :swapped, swapped)

    ## register global `scaling` vector
    map!(graph, [:scale, :reversed, :swapped], :scaling) do fac, reversed, swapped
        scaling = fill(fac, 2)
        if reversed
            scaling[1] *= -1
        end
        if swapped
            scaling[2] *= -1
        end
        return scaling
    end
    
    add_input!(graph, :linecap, convert(Symbol, linecap))
    add_input!(graph, :joinstyle, convert(Symbol, joinstyle))
 
    ## drawing attributes
    add_input!((k, v) -> convert(Bool, v), graph, :is_filled, is_filled)
    add_input!((k, v) -> convert(Bool, v), graph, :is_stroked, is_stroked)
    add_input!(makeref, graph, :color, color)               # makeref to allow for `nothing` value
    add_input!(makeref, graph, :fillcolor, fillcolor)
    add_input!(makeref, graph, :strokecolor, strokecolor)

    return graph
end

function _common_bezigon_setup(
    spec; 
    ## parent attributes, could be Observables
    color=nothing, linewidth=1f0, visible=true, kwargs...
)
    graph = _connected_child_graph(spec.attrs)
    add_input!(graph, :color_parent, color)
    add_input!(graph, :strokewidth_parent, linewidth)
    add_input!(graph, :visible, visible)
    graph = _register_common_computations!(graph)
    return graph
end

function _connected_child_graph(attrs)
    ## TODO we don't really need a new graph, do we?
    graph = ComputeGraph()
    for symb in (
        :length, :width, :width_, :angle, :angle_, :line_width, :line_width_,
        :align, :miter_limit, :reversed, :swapped, :scaling, :linecap, :joinstyle, 
        :is_filled, :is_stroked, :color, :strokecolor, :fillcolor
    )
        add_input!(graph, symb, attrs[symb])
    end
    return graph
end

function _register_common_computations!(graph)
    ## register `arrowsize` and `arrow_length`
    map!(
        graph, [:strokewidth_parent, :length, :angle, :angle_, :width, :width_], 
        [:arrow_length, :arrow_width]
    ) do sw, len, ang, ang_, wdth, wdth_
        arrow_width = NaN32
        if !isnan(len[1])
            L, l = len
            arrow_length = L + l * sw
        elseif !isnan(ang[1])
            θ, L, l = ang
            arrow_length = (L + l * sw) * cos(θ)
            arrow_width = (2 * sin(θ/2)) * arrow_length
        else
            error("Cannot determine arrow length. Give valid `length` or `angle` specification.")
        end
        if isnan(arrow_width)
            if !isnan(wdth[1])
                L, l = wdth
                arrow_width = L + l * sw
            elseif !isnan(wdth_[1])
                L, l = wdth_
                arrow_width = L + l * arrow_length
            elseif !isnan(ang_)
                arrow_width = tan(ang_ / 2) * arrow_length
            else
                error("Cannot determine arrow width. Give valid `width`, `width_`, `angle`, or `angle_` specification.")
            end
        end
        return arrow_length, arrow_width
    end

    map!(
        graph, [:line_width, :line_width_, :strokewidth_parent, :arrow_length], :strokewidth
    ) do lw, lw_, sw, len
        if !isnan(lw[1]) 
            fw = lw[1] + lw[2] * sw
        elseif !isnan(lw_[2])
            fw = lw_[1] + lw_[2] * len
        else
            @debug "Cannot determine arrow strokewidth, using parent strokewidth."
            fw = sw
        end
        return f32(fw)
    end

    return graph
end

function _register_appearance_attributes!(graph)

    if !haskey(graph, :strokewidth_fin)
        alias!(graph, :strokewidth, :strokewidth_fin)
    end

    ## get strokecolor, try specific color first, then general color
    map!(graph, [:strokecolor, :color], :strokecolor_user) do strokecolor, color
        isnothing(strokecolor) && return to_color(color)
        return to_color(strokecolor)
    end

    ## get fillcolor, try specific color first, then general color
    map!(graph, [:fillcolor, :color], :fillcolor_user) do fillcolor, color
        isnothing(fillcolor) && return to_color(color)
        return to_color(fillcolor)
    end

    ## final colors can inherit from parent
    map!(
        graph, [:fillcolor_user, :color_parent], :fillcolor_fin
    ) do color_user, color_parent
        isnothing(color_user) && return to_color(color_parent)
        return to_color(color_user)
    end
    map!(
        graph, [:strokecolor_user, :color_parent], :strokecolor_fin
    ) do color_user, color_parent
        isnothing(color_user) && return to_color(color_parent)
        return to_color(color_user)
    end

    ## Bool `stroke_bgon` to indicate whether or not outlines are drawn
    map!(
        graph, [:strokecolor_user, :strokewidth_fin, :is_stroked], :stroke_bgon
    ) do color, strokewidth, is_stroked
        !is_stroked && return false
        iszero(strokewidth) && return false
        color === :transparent && return false
        (!isnothing(color) && iszero(Makie.alpha(to_color(color)))) && return false
        return true
    end

    ## Bool `fill_bgon` to indicate whether or not bezigon is filled
    map!(graph, [:fillcolor_user, :is_filled], :fill_bgon) do color, is_filled
        !is_filled && return is_filled
        color === :transparent && return false
        (!isnothing(color) && iszero(Makie.alpha(to_color(color)))) && return false
        return true
    end

    map!(graph, [:visible, :fill_bgon], :fill_visible) do vis, flag
        return vis && flag
    end

    map!(graph, [:visible, :stroke_bgon], :stroke_visible) do vis, flag
        return vis && flag
    end

    return graph
end

function _register_bezigon_drawing_cmds!(graph)    
    add_input!(v32, graph, :target!, zeros(Float32, 2))
    add_input!(f32, graph, :rotation, 0f0)

    map!(
        graph, 
        [:target!, :line_end, :tip_end, :align, :scaling, :rotation], 
        [:shift, :connector, :shrinksize]
    ) do target!, line_end, tip_end, align, size, ang
        #= 
        * `target!` = where arrow should point modulo alignment; 
        * compute a `shift` vector for bezigon `path0` accordingly
        * `connector`: line end vector after shifting
        * `align==0` => tip end vector equals `target!`, exact pointing
        * `align==1` => line end vector equals `target!`, overshooting
        =#

        anchor = _rotate_2d(_scale_2d([(1-align) * tip_end + align * line_end, 0f0], size), ang)
        connector = _rotate_2d(_scale_2d([line_end, 0f0], size), ang)
        shift = target! .- anchor 
        connector .+= shift
        shrinksize = sqrt(sum(abs2, target! .- connector))
        return (shift, v32(connector), shrinksize)
    end

    ## `:connector` is set above, only `:path` left to do
    map!(_scale_rotate_translate_poly, graph, [:path0, :scaling, :rotation, :shift], :path)

    return graph
end

## geometry helpers
function _scale_rotate_translate_2d(vec, size, rotation, shift)
    v = _scale_2d(vec, size)
    x, y = _rotate_2d(v, rotation)
    z = [x, y]
    return z .+ shift
end
_scale_2d(vec, size) = size .* vec
function _rotate_2d(v, rotation)
    _x, _y = v
    x = _x * cos(rotation) - _y * sin(rotation)
    y = _y * cos(rotation) + _x * sin(rotation)
    return [x; y]
end

function _scale_rotate_translate_poly(poly, size, rotation, shift)
    if size isa AbstractVector
        size = convert(Point2d, size)
    end
    p = scale(poly, size)
    p = rotate(p, rotation)
    return translate(p, convert(Point2d, shift))
end

## simple helpers
f32(t) = convert(Float32, t)
f32(k, t) = f32(t)
v32(t) = convert(Vector{Float32}, t)
v32(k, t) = v32(t)
boolean(k, t) = boolean(t)
boolean(t) = convert(Bool, t)
makeref(t, T::Type)=Ref{T}(t)
makeref(k, t) = makeref(t, Any)

validate_tup2(spec::Number, symb) = (f32(spec), 0f0)
validate_tup2(spec::Tuple{<:Number, <:Number}, symb) = f32.(spec)
validate_tup2(spec, symb) = validate_tup(spec, symb, Union{Number, Tuple{<:Number, <:Number}})

function validate_tup(spec, symb, T=Union{})
    error("Argument `:$(symb)` requires value of type `$(T)`.")
end

make_spec1(k, v::Nothing) = NaN32
make_spec1(k, v::Number) = f32(v)
make_spec1(k, v) = validate_tup(v, k, Number)

make_spec2(k, v::Nothing) = (NaN32, NaN32)
make_spec2(k, v) = validate_tup2(v, k)

make_spec3(k, v::Nothing) = (NaN32, NaN32, NaN32)
make_spec3(k, v::Tuple{<:Number, <:Number}) = (f32.(v)..., 0f0)
make_spec3(k, v::Tuple{<:Number, <:Number, <:Number}) = f32.(v)
make_spec3(k, v) = validate_tup(v, k, Union{Tuple{<:Number, <:Number}, Tuple{<:Number, <:Number, <:Number}})


function ParsedSpec(
    spec_str::AbstractString;
    kwargs...
)
    if spec_str == ">"
        return ComputerModernRightarrowTip(; kwargs...)
    elseif spec_str == "<"
        return ComputerModernRightarrowTail(; kwargs...)
    end
    @warn "Could not parse arrow specification string \"$(spec_str)\"."
    return NoDecoSpec()
end

end#module FancyArrows