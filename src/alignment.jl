struct Alignment
    reference::String
    query::String
    score::Int
    mapping::Vector{Int} # query position → reference position; zero means insertion
    substitutions::Int
    indels::Int
end

"Banded global Needleman–Wunsch; linear gap cost. O(L × band) trace storage."
function alignSequences(reference::DNAInput, query::DNAInput;
                        band=16, match=5, mismatch=-4, gap=-8)
    band >= 0 && match > 0 && mismatch <= 0 && gap < 0 || throw(ArgumentError("Invalid alignment parameters"))
    a = codeunits(validated_dna(reference)); b = codeunits(validated_dna(query))
    m = length(a); n = length(b)
    width = max(Int(band), abs(m-n))
    sentinel = typemin(Int) ÷ 4
    prev = fill(sentinel, n+1); curr = fill(sentinel, n+1)
    trace = zeros(UInt8, 2width+1, m+1)
    boundary_bound = sentinel
    upperbound(i,j,s) = s + min(m-i,n-j)*match + abs((m-i)-(n-j))*gap
    for j in 0:min(n, width)
        prev[j+1] = j*gap
        j > 0 && (trace[j+width+1,1] = 3)
        j == width && (boundary_bound = max(boundary_bound,upperbound(0,j,prev[j+1])))
    end
    for i in 1:m
        lo = max(0, i-width); hi = min(n, i+width)
        fill!(@view(curr[max(1,lo):min(n+1,hi+2)]), sentinel)
        if lo == 0
            curr[1] = i*gap
            trace[width+1-i,i+1] = 2
            i == width && (boundary_bound = max(boundary_bound,upperbound(i,0,curr[1])))
        end
        for j in max(1,lo):hi
            diag = prev[j] + (a[i] == b[j] ? match : mismatch)
            up = prev[j+1] + gap; left = curr[j] + gap
            best, dir = diag, UInt8(1)
            up > best && ((best, dir) = (up, UInt8(2)))
            left > best && ((best, dir) = (left, UInt8(3)))
            curr[j+1] = best
            trace[j-i+width+1,i+1] = dir
            abs(j-i) == width && (boundary_bound = max(boundary_bound,upperbound(i,j,best)))
        end
        prev, curr = curr, prev
    end
    i = m; j = n
    map = zeros(Int, n); aa = UInt8[]; bb = UInt8[]
    subs = indels = 0; touched = false
    while i > 0 || j > 0
        abs(j-i) == width && i > 0 && j > 0 && (touched = true)
        dir = trace[j-i+width+1,i+1]
        if dir == 1
            push!(aa, a[i]); push!(bb, b[j]); map[j] = i
            subs += a[i] != b[j]; i -= 1; j -= 1
        elseif dir == 2
            push!(aa, a[i]); push!(bb, 0x2d); i -= 1; indels += 1
        elseif dir == 3
            push!(aa, 0x2d); push!(bb, b[j]); j -= 1; indels += 1
        else
            error("Alignment traceback escaped band")
        end
    end
    # Any better path escaping the band must cross a boundary whose optimistic
    # remaining score exceeds the current result. Expand until certified exact.
    if (touched || boundary_bound > prev[n+1]) && width < max(m,n)
        return alignSequences(reference, query; band=min(max(m,n), max(1,2width)), match, mismatch, gap)
    end
    Alignment(String(reverse(aa)), String(reverse(bb)), prev[n+1], map, subs, indels)
end

"Sorted packed k-mer multiset. O(L) storage instead of 4^k per read."
function kmers(seq::DNAInput, k::Int=5)
    1 <= k <= 15 || throw(ArgumentError("k must be in 1:15"))
    result = UInt32[]; bits = UInt32(0); mask = (UInt32(1) << (2k))-1; run = 0
    for b in codeunits(seq)
        if !(b in DNA)
            bits = 0; run = 0
            continue
        end
        bits = ((bits << 2) | UInt32(baseindex(b)-1)) & mask
        run += 1
        run >= k && push!(result, bits)
    end
    sort!(result)
end
function kmer_distance(a::Vector{UInt32}, b::Vector{UInt32})
    min(length(a),length(b)) == 0 && return 0.0
    i = j = 1; common = 0
    while i <= length(a) && j <= length(b)
        if a[i] == b[j]
            common += 1; i += 1; j += 1
        elseif a[i] < b[j]
            i += 1
        else
            j += 1
        end
    end
    1 - common / min(length(a),length(b))
end
