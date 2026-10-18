using DADA2, Dates, TOML

const DEFAULTS = Dict(
    "forward"=>"", "reverse"=>"", "samples"=>"", "sample"=>"sample1", "out"=>"",
    "trunc-f"=>"0", "trunc-r"=>"0", "trim-left-f"=>"0", "trim-left-r"=>"0",
    "max-ee-f"=>"2", "max-ee-r"=>"2", "trunc-q"=>"2", "min-len"=>"20",
    "min-overlap"=>"12", "min-asv-length"=>"1", "max-asv-length"=>"0",
    "max-iter"=>"10", "quality-mode"=>"raw", "chimera-method"=>"consensus")

function usage()
    println("""
    Analysis of real FASTQ files with DADA2.jl (Phred+33; plain or gzip).

    julia --project=. --threads=4 examples/real_data.jl --forward F.fastq.gz --out results
    julia --project=. --threads=4 examples/real_data.jl --forward F.fastq.gz --reverse R.fastq.gz --out results --trunc-f 240 --trunc-r 160
    julia --project=. --threads=4 examples/real_data.jl --samples samples.tsv --out results

    Input: --forward FILE [--reverse FILE] [--sample NAME], or --samples TSV.
    The TSV has the header sample<TAB>forward<TAB>reverse; leave reverse empty for single-end.
    Relative paths in the TSV are resolved against the TSV directory.
    --out DIR                 New or empty directory (required).
    --trunc-f N --trunc-r N    Truncation lengths; 0 = variable (default).
    --trim-left-f N --trim-left-r N   Leading bases to remove (default 0).
    --max-ee-f X --max-ee-r X  Maximum expected errors (default 2).
    --trunc-q N               Truncate at the first quality <= N (default 2).
    --min-len N               Minimum length after filtering (default 20).
    --min-overlap N           Minimum paired-end overlap (default 12).
    --min-asv-length N        Minimum ASV length (default 1).
    --max-asv-length N        Maximum ASV length; 0 = no limit.
    --max-iter N              Maximum learning iterations (default 10).
    --quality-mode raw|mean   Original or mean qualities (default raw).
    --chimera-method consensus|pooled   Chimera method (default consensus).
    --help                    Show this message.
    """)
end

function options(args)
    config = copy(DEFAULTS)
    iseven(length(args)) || error("Every option requires a value. Use --help.")
    seen = Set{String}()
    for i in 1:2:length(args)
        startswith(args[i],"--") || error("Invalid option: $(args[i])")
        key = args[i][3:end]
        haskey(config,key) || error("Unknown option: --$key")
        key in seen && error("Repeated option: --$key")
        push!(seen,key); config[key] = args[i+1]
    end
    isempty(config["out"]) && error("Specify --out DIR.")
    for key in ("trunc-f","trunc-r","trim-left-f","trim-left-r","trunc-q",
                "min-len","min-overlap","min-asv-length","max-asv-length","max-iter")
        value = tryparse(Int,config[key])
        value !== nothing && value >= 0 || error("--$key requires a non-negative integer.")
    end
    for key in ("min-len","min-overlap","min-asv-length","max-iter")
        parse(Int,config[key]) > 0 || error("--$key must be positive.")
    end
    parse(Int,config["trunc-q"]) <= 93 || error("--trunc-q must be between 0 and 93.")
    upper = parse(Int,config["max-asv-length"])
    upper == 0 || upper >= parse(Int,config["min-asv-length"]) || error("Inconsistent ASV length limits.")
    for key in ("max-ee-f","max-ee-r")
        value = tryparse(Float64,config[key])
        value !== nothing && isfinite(value) && value >= 0 || error("--$key requires a finite non-negative number.")
    end
    config["quality-mode"] in ("raw","mean") || error("--quality-mode must be raw or mean.")
    config["chimera-method"] in ("consensus","pooled") || error("--chimera-method must be consensus or pooled.")
    config
end

function samples(config)
    rows = NamedTuple{(:sample,:forward,:reverse),Tuple{String,String,String}}[]
    if !isempty(config["samples"])
        isempty(config["forward"]) && isempty(config["reverse"]) || error("Use either --samples or --forward/--reverse.")
        manifest = abspath(config["samples"])
        lines = readlines(manifest)
        !isempty(lines) && split(replace(lines[1], '\ufeff'=>""),'\t';keepempty=true) ==
            ["sample","forward","reverse"] || error("Required TSV header: sample, forward, reverse.")
        for (i,line) in enumerate(lines[2:end])
            isempty(strip(line)) && continue
            fields = String.(split(line,'\t';keepempty=true))
            length(fields) == 3 || error("Three columns required on line $(i+1).")
            name,f,r = fields
            push!(rows,(sample=name,forward=abspath(joinpath(dirname(manifest),f)),
                        reverse=isempty(r) ? "" : abspath(joinpath(dirname(manifest),r))))
        end
    else
        !isempty(config["forward"]) || error("Specify --forward FILE or --samples TSV.")
        push!(rows,(sample=config["sample"],forward=abspath(config["forward"]),
                    reverse=isempty(config["reverse"]) ? "" : abspath(config["reverse"])))
    end
    !isempty(rows) || error("No samples specified.")
    length(unique(row.sample for row in rows)) == length(rows) || error("Duplicate sample names.")
    paired = !isempty(rows[1].reverse)
    for row in rows
        occursin(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$",row.sample) || error("Invalid sample name: $(row.sample). Use letters, digits, _, - or a dot.")
        (!isempty(row.reverse)) == paired || error("All samples must use the same single/paired-end mode.")
        isfile(row.forward) || error("FASTQ not found: $(row.forward)")
        if paired
            isfile(row.reverse) || error("FASTQ not found: $(row.reverse)")
            realpath(row.forward) != realpath(row.reverse) || error("Forward and reverse must be different files.")
        end
    end
    # File names must also be unique on case-insensitive file systems.
    length(unique(lowercase(row.sample) for row in rows)) == length(rows) || error("Sample names differ only by letter case.")
    rows,paired
end

function filteroptions(config,direction)
    FilterOptions(truncLen=parse(Int,config["trunc-$direction"]),
        trimLeft=parse(Int,config["trim-left-$direction"]),maxEE=parse(Float64,config["max-ee-$direction"]),
        truncQ=parse(Int,config["trunc-q"]),minLen=parse(Int,config["min-len"]))
end

function writecounts(path,table,ids)
    open(path,"w") do io
        println(io,join(["sample";ids],'\t'))
        for i in eachindex(table.samples)
            println(io,join([table.samples[i];string.(table.counts[i,:])],'\t'))
        end
    end
end

function writelearning(out,direction,learned)
    diag = errorDiagnostics(learned)
    open(joinpath(out,"errors_$direction.tsv"),"w") do io
        println(io,join(["transition";["Q$q" for q in diag.qualities]],'\t'))
        for from in 1:4, to in 1:4
            row = 4(from-1)+to
            println(io,join(["$("ACGT"[from])2$("ACGT"[to])";string.(learned.errors.rates[row,:])],'\t'))
        end
    end
    open(joinpath(out,"error_support_$direction.tsv"),"w") do io
        println(io,"source\tquality\texposures\tobserved_substitutions\texpected_substitutions\tsupported")
        for b in 1:4, q in eachindex(diag.qualities)
            println(io,join(("ACGT"[b],diag.qualities[q],diag.exposures[b,q],diag.observed_substitutions[b,q],
                diag.expected_substitutions[b,q],diag.supported[b,q]),'\t'))
        end
    end
    open(joinpath(out,"learning_$direction.tsv"),"w") do io
        println(io,"iteration\tmax_log_change")
        for (i,delta) in enumerate(learned.history)
            println(io,i,'\t',delta)
        end
    end
    learned.converged || @warn "Error learning did not converge" direction iterations=learned.iterations
end

function main(args)
    if isempty(args) || args == ["--help"]
        usage(); return
    end
    config = options(args); input,paired = samples(config)
    out = abspath(config["out"])
    ispath(out) && (!isdir(out) || !isempty(readdir(out))) && error("The output directory must be new or empty: $out")
    mkpath(joinpath(out,"filtered"))
    fpaths = [joinpath(out,"filtered",row.sample*"_F.fastq.gz") for row in input]
    rpaths = paired ? [joinpath(out,"filtered",row.sample*"_R.fastq.gz") for row in input] : String[]
    stats = FilterStats[]
    println("Filtering $(length(input)) samples...")
    open(joinpath(out,"filtering.tsv"),"w") do io
        println(io,"sample\treads_in\treads_filtered")
        for (i,row) in enumerate(input)
            stat = paired ? filterAndTrim(row.forward,row.reverse,fpaths[i],rpaths[i];
                optionsF=filteroptions(config,"f"),optionsR=filteroptions(config,"r")) :
                filterAndTrim(row.forward,fpaths[i];options=filteroptions(config,"f"))
            push!(stats,stat)
            println(io,join((row.sample,stat.reads_in,stat.reads_out),'\t')); flush(io)
            println("  $(row.sample): $(stat.reads_in) → $(stat.reads_out)")
        end
    end
    active = findall(s -> s.reads_out > 0,stats)
    isempty(active) && error("No reads retained. Check the parameters and filtering.tsv.")
    println("Shared error learning ($(config["quality-mode"]))...")
    ef = learnErrors(fpaths[active];qualityMode=Symbol(config["quality-mode"]),maxIter=parse(Int,config["max-iter"]))
    writelearning(out,"forward",ef)
    er = paired ? learnErrors(rpaths[active];qualityMode=Symbol(config["quality-mode"]),maxIter=parse(Int,config["max-iter"])) : nothing
    paired && writelearning(out,"reverse",er)
    results = []; tracking = []
    println("ASV inference",paired ? " and pair merging..." : "...")
    for (i,row) in enumerate(input)
        if stats[i].reads_out == 0
            push!(results,MergeResult(String[],Int[],0,0))
            push!(tracking,(forward=0,reverse=paired ? 0 : "NA",merged=paired ? 0 : "NA",converged="NA"))
            continue
        end
        f = dada(derepFastq(fpaths[i];keepmap=paired),ef)
        if paired
            r = dada(derepFastq(rpaths[i]),er)
            merged = mergePairs(f,r;minOverlap=parse(Int,config["min-overlap"]))
            push!(results,merged)
            push!(tracking,(forward=sum(f.counts),reverse=sum(r.counts),merged=merged.accepted,converged=f.converged && r.converged))
        else
            push!(results,f)
            push!(tracking,(forward=sum(f.counts),reverse="NA",merged="NA",converged=f.converged))
        end
        tracking[end].converged || @warn "ASV inference did not converge" sample=row.sample
    end
    alltable = makeSequenceTable(results;samples=[row.sample for row in input])
    ids = ["ASV$i" for i in eachindex(alltable.sequences)]
    writeFasta(joinpath(out,"asvs_all.fasta"),alltable.sequences;ids)
    writecounts(joinpath(out,"asv_counts_all.tsv"),alltable,ids)
    minlen = parse(Int,config["min-asv-length"]); maxlen = parse(Int,config["max-asv-length"])
    keep = findall(s -> length(s) >= minlen && (maxlen == 0 || length(s) <= maxlen),alltable.sequences)
    lengthtable = SequenceTable(alltable.samples,alltable.sequences[keep],alltable.counts[:,keep])
    clean = removeBimeraDenovo(lengthtable;method=Symbol(config["chimera-method"]))
    idmap = Dict(zip(alltable.sequences,ids)); cleanids = [idmap[s] for s in clean.sequences]
    writeFasta(joinpath(out,"asvs.fasta"),clean.sequences;ids=cleanids)
    writecounts(joinpath(out,"asv_counts.tsv"),clean,cleanids)
    kept = Set(clean.sequences); lengthkept = Set(lengthtable.sequences)
    open(joinpath(out,"asv_summary.tsv"),"w") do io
        println(io,"ASV\tlength\ttotal_reads\tstatus")
        for (j,s) in enumerate(alltable.sequences)
            status = s in kept ? "kept" : s in lengthkept ? "chimera" : "length"
            println(io,join((ids[j],length(s),sum(alltable.counts[:,j]),status),'\t'))
        end
    end
    open(joinpath(out,"read_tracking.tsv"),"w") do io
        println(io,"sample\treads_in\treads_filtered\tdenoised_forward\tdenoised_reverse\tmerged\tafter_length_filter\tnonchimeric\tfinal_asvs\tinference_converged")
        for (i,row) in enumerate(input)
            t = tracking[i]
            println(io,join((row.sample,stats[i].reads_in,stats[i].reads_out,t.forward,t.reverse,t.merged,
                sum(lengthtable.counts[i,:]),sum(clean.counts[i,:]),count(>(0),clean.counts[i,:]),t.converged),'\t'))
        end
    end
    metadata = Dict("package_version"=>string(pkgversion(DADA2)),"julia_version"=>string(VERSION),
        "threads"=>Threads.nthreads(),"completed_utc"=>string(now(UTC)),"mode"=>paired ? "paired" : "single",
        "arguments"=>args,"options"=>config,
        "learning_forward_converged"=>ef.converged,"learning_forward_iterations"=>ef.iterations,
        "samples"=>[Dict("sample"=>row.sample,"forward"=>row.forward,"reverse"=>row.reverse) for row in input])
    if paired
        metadata["learning_reverse_converged"] = er.converged
        metadata["learning_reverse_iterations"] = er.iterations
    end
    open(joinpath(out,"run.toml"),"w") do io
        TOML.print(io,metadata;sorted=true)
    end
    println("Completed: $(length(clean.sequences)) ASVs, $(sum(clean.counts)) final reads. Results: $out")
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        main(ARGS)
    catch err
        println(stderr,"Analysis stopped: ",sprint(showerror,err))
        exit(1)
    end
end
