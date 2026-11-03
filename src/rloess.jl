# One-predictor translation of R's `loess(y ~ x, weights=w, degree=2)` with
# the default interpolated surface (netlib dloess in R's stats/src/loessf.f:
# lowesd, lowesb with ehg126/ehg124/ehg129/ehg139/ehg127, lowese with ehg128)
# and of DADA2's `loessErrfun`. The k-d tree, the local fits at its vertices
# (value and slope) and the cubic Hermite interpolation follow the Fortran;
# linear algebra uses LAPACK instead of LINPACK, so results agree to rounding.

"Kd-tree cells of a one-dimensional interpolated loess surface."
struct LoessSurface
    vertices::Vector{Float64}
    values::Vector{Float64}   # fitted value at each vertex
    slopes::Vector{Float64}   # fitted derivative at each vertex
    split::Vector{Int}        # 0 for a leaf, otherwise 1 (the only predictor)
    xi::Vector{Float64}       # split value
    lo::Vector{Int}
    hi::Vector{Int}
    lower::Vector{Int}        # vertex indices of each cell
    upper::Vector{Int}
    range::Tuple{Float64,Float64}
end

"Upstream `ehg127`: local quadratic fit at `q`; returns value and slope."
function loess_vertex_fit(q::Float64, x, y, w, nf::Int, span::Float64)
    n = length(x)
    dist = [(x[i]-q)^2 for i in 1:n]
    psi = partialsortperm(dist, 1:nf)
    rho = dist[psi[nf]] * max(1.0, span)
    rho > 0 || throw(ArgumentError("zero-width loess neighbourhood; make span bigger"))
    # All-zero weights only warn upstream; the pseudoinverse then gives zero.
    weights = [sqrt(w[i] * (1 - sqrt(dist[i]/rho)^3)^3) for i in psi]
    b = Matrix{Float64}(undef, nf, 3)
    for (r, i) in enumerate(psi)
        b[r,1] = weights[r]
        b[r,2] = weights[r]*(x[i]-q)
        b[r,3] = weights[r]*(x[i]-q)^2
    end
    eta = [weights[r]*y[i] for (r, i) in enumerate(psi)]
    colnorm = ones(3)
    for j in 1:3
        scale = sqrt(sum(abs2, @view b[:,j]))
        if scale > 0
            b[:,j] ./= scale; colnorm[j] = scale
        end
    end
    # QR without pivoting, then the SVD of R: a truncated pseudoinverse.
    k = min(nf, 3)
    f = qr(b)
    qty = f.Q' * eta
    r = zeros(3, 3)
    r[1:k, :] .= f.R[1:k, :]
    s = svd(r)
    tol = s.S[1] * (100 * eps())
    qty3 = [i <= length(qty) ? qty[i] : 0.0 for i in 1:3]
    gamma = [s.S[j] > tol ? dot(@view(s.U[:,j]), qty3) / s.S[j] : 0.0 for j in 1:3]
    coefficients = [sum(s.V[j,i]/colnorm[j] * gamma[i] for i in 1:3) for j in 1:3]
    coefficients[1], coefficients[2]
end
loess_vertex_fit(q, x, y, w, nf, span) = loess_vertex_fit(Float64(q), x, y, w, nf, Float64(span))

"Build the interpolated loess surface (upstream `lowesb`, one predictor, degree 2)."
function loess_surface(x::AbstractVector{<:Real}, y::AbstractVector{<:Real}, w::AbstractVector{<:Real};
                       span::Real=0.75, cell::Real=0.2)
    n = length(x)
    length(y) == n == length(w) && n > 0 || throw(ArgumentError("Invalid loess data"))
    x = Float64.(x); y = Float64.(y); w = Float64.(w)
    span = Float64(span)
    nf = min(n, floor(Int, n*span + 1e-5))
    nf > 0 || throw(ArgumentError("span is too small"))
    fc = floor(Int, n * (span*Float64(cell)))
    nvmax = max(200, n); ncmax = nvmax
    # ehg126: bounding box, expanded by 0.5% of the range.
    alpha, beta = extrema(x)
    mu = 0.005 * max(beta - alpha, 1e-10*max(abs(alpha), abs(beta)) + 1e-30)
    vertices = [alpha - mu, beta + mu]
    # ehg124: cells hold contiguous ranges of the sorted predictor.
    xs = sort(x)
    lo = [1]; hi = [n]; split = Int[]; xis = Float64[]; lower = [1]; upper = [2]
    p = 1
    while p <= length(lo)
        l = lo[p]; u = hi[p]
        leaf = u - l + 1 <= fc || vertices[upper[p]] - vertices[lower[p]] <= 0 ||
               ncmax < length(lo) + 2 || nvmax < length(vertices) + 1
        m = 0
        if !leaf
            m = (l + u) ÷ 2
            # All ties go with the upper child: look for the nearest gap.
            offset = 0
            while l <= m + offset < u
                if xs[m+offset] == xs[m+offset+1]
                    offset = -offset
                    offset >= 0 && (offset += 1)
                else
                    m += offset
                    break
                end
            end
            leaf = vertices[lower[p]] == xs[m] || vertices[upper[p]] == xs[m]
        end
        if leaf
            push!(split, 0); push!(xis, 0.0)
        else
            push!(split, 1); push!(xis, xs[m])
            # ehg125: reuse an identical vertex.
            v = findfirst(==(xs[m]), vertices)
            v === nothing && (push!(vertices, xs[m]); v = length(vertices))
            push!(lo, l, m+1); push!(hi, m, u)
            push!(lower, lower[p], v); push!(upper, v, upper[p])
            lo[p], hi[p] = length(lo) - 1, length(lo)
        end
        p += 1
    end
    fits = [loess_vertex_fit(v, x, y, w, nf, span) for v in vertices]
    LoessSurface(vertices, first.(fits), last.(fits), split, xis, lo, hi, lower, upper, (alpha, beta))
end

"Upstream `ehg128` in one dimension; NaN outside the range of the data, as `predict`."
function (s::LoessSurface)(z::Real)
    z = Float64(z)
    s.range[1] <= z <= s.range[2] || return NaN
    j = 1
    while s.split[j] != 0
        j = z <= s.xi[j] ? s.lo[j] : s.hi[j]
    end
    ll = s.lower[j]; ur = s.upper[j]
    width = s.vertices[ur] - s.vertices[ll]
    h = (z - s.vertices[ll]) / width
    phi0 = (1-h)^2*(1+2*h); phi1 = h^2*(3-2*h)
    psi0 = h*(1-h)^2; psi1 = h^2*(h-1)
    phi0*s.values[ll] + phi1*s.values[ur] + (psi0*s.slopes[ll] + psi1*s.slopes[ur]) * width
end

# Sum of doubles in x87 extended precision, then rounded, as R's colSums.
# MPFR rounds each addition into the 64-bit destination, independently of the
# caller's BigFloat precision (DD8).
function long_double_sum(values)
    acc = BigFloat(0; precision=64)
    for v in values
        ccall((:mpfr_add_d, Base.MPFR.libmpfr), Cint, (Ref{BigFloat}, Ref{BigFloat}, Cdouble, Cint),
              acc, acc, Float64(v), 0)  # MPFR_RNDN
    end
    Float64(acc)
end

"""
DADA2 `loessErrfun`: per substitution, R's interpolated quadratic LOESS of
`log10((errors+1)/exposure)` on quality, weighted by exposure, predicted at
every quality, extended flat beyond the predicted range, clamped to
[1e-7, 0.25]; self-transitions are one minus the substitutions. A source base
without any exposure keeps the rates of `fallback`; without a fallback this
is an error, as in R. Sources observed at fewer than `minpoints` qualities take
their rows from `sparse` when given.
"""
function r_loess_errors(trans::Matrix{Int}, fallback::Union{Nothing,ErrorModel}; span=0.75,
                        sparse::Union{Nothing,ErrorModel}=nothing, minpoints::Int=4)
    nq = size(trans, 2)
    rates = fallback === nothing ? zeros(16, nq) : copy(fallback.rates)
    qualities = collect(0.0:nq-1)
    for from in 1:4
        rows = 4(from-1)+1:4from
        totals = vec(sum(@view(trans[rows,:]); dims=1))
        observed = findall(>(0), totals)
        if isempty(observed)
            fallback === nothing &&
                throw(ArgumentError("Error rates could not be estimated: no exposure of base $("ACGT"[from])"))
            continue
        end
        if sparse !== nothing && length(observed) < minpoints
            rates[rows, :] .= sparse.rates[rows, :]
            continue
        end
        alternatives = [4(from-1)+to for to in 1:4 if to != from]
        for row in alternatives
            response = [log10((trans[row,q]+1)/totals[q]) for q in observed]
            surface = loess_surface(qualities[observed], response, totals[observed]; span)
            pred = surface.(qualities)
            defined = findall(!isnan, pred)
            pred[1:first(defined)-1] .= pred[first(defined)]
            pred[last(defined)+1:end] .= pred[last(defined)]
            rates[row,:] .= clamp.(exp10.(pred), 1e-7, 0.25)
        end
        for q in 1:nq
            rates[4(from-1)+from, q] = 1 - long_double_sum(rates[row,q] for row in alternatives)
        end
    end
    ErrorModel(rates)
end
