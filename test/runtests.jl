using Test, DADA2, BioSequences, Random

function fixture(sequences, counts; q=40)
    Derep(sequences, counts, [fill(Float64(q),length(s)) for s in sequences])
end
function dense_score(a,b; match=5,mismatch=-4,gap=-8)
    scores = zeros(Int,length(a)+1,length(b)+1)
    scores[:,1] .= (0:length(a)).*gap
    scores[1,:] .= (0:length(b)).*gap
    for i in 1:length(a), j in 1:length(b)
        scores[i+1,j+1] = max(scores[i,j]+(a[i]==b[j] ? match : mismatch),scores[i,j+1]+gap,scores[i+1,j]+gap)
    end
    scores[end,end]
end

@testset "FASTQ streaming and filtering" begin
    mktempdir() do dir
        raw = joinpath(dir,"reads.fastq")
        write(raw,"@r1/1\nACGTACGT\n+\nIIIIIIII\n@r2/1\nACGTACGT\n+\nHHHHHHHH\n@r3/1\nACGTNCGT\n+\nIIIIIIII\n")
        out = joinpath(dir,"filtered.fastq.gz")
        stats = filterAndTrim(raw,out;options=FilterOptions(minLen=4))
        @test (stats.reads_in,stats.reads_out) == (3,2)
        derep = derepFastq(out)
        @test derep.sequences == ["ACGTACGT"]
        @test derep.counts == [2]
        @test derep.qualities == [fill(39.5,8)]
        @test derep.map == [1,1]
        @test_throws ArgumentError filterAndTrim(raw,raw)
        @test_throws ArgumentError filterAndTrim(raw,out)
        @test_throws ArgumentError derepFastq(raw)
        rev = joinpath(dir,"reverse.fastq")
        write(rev,replace(read(raw,String),"/1"=>"/2"))
        of = joinpath(dir,"f.fastq"); or = joinpath(dir,"r.fastq")
        stats = filterAndTrim(raw,rev,of,or;optionsF=FilterOptions(minLen=4),optionsR=FilterOptions(minLen=4))
        @test stats.reads_out == 2
        @test derepFastq(of).map == derepFastq(or).map
        write(rev,"@wrong/2\nACGT\n+\nIIII\n")
        @test_throws ArgumentError filterAndTrim(raw,rev,joinpath(dir,"badf"),joinpath(dir,"badr"))
        @test !isfile(joinpath(dir,"badf")) && !isfile(joinpath(dir,"badr"))
        write(joinpath(dir,"truncated"),"@read\nACGT\n+\n")
        @test_throws ArgumentError derepFastq(joinpath(dir,"truncated"))
        write(joinpath(dir,"crlf.fastq"),"@r1 x/1\r\nacgt\r\n+\r\nII#!\r\n")
        reader = FastqReader(joinpath(dir,"crlf.fastq"))
        records = try collect(reader) finally close(reader) end
        @test (records[1].id,records[1].sequence,records[1].quality) == ("r1 x/1","ACGT",UInt8[40,40,2,0])
        write(joinpath(dir,"plus.fastq"),"@r1\nACGT\n+r2\nIIII\n")
        @test_throws ArgumentError derepFastq(joinpath(dir,"plus.fastq"))
        @test_throws ArgumentError FastqRecord("r","AC\xff",[40,40,40])
        write(joinpath(dir,"spaces.fa"),">s1 Bacteria;\nAC GT \nacgn\n")
        @test readFasta(joinpath(dir,"spaces.fa")) == (ids=["s1 Bacteria;"],sequences=["ACGTACGN"])
        qraw = joinpath(dir,"low.fastq")
        write(qraw,"@low\nACGTACGT\n+\nIIII!III\n")
        stats = filterAndTrim(qraw,joinpath(dir,"qtrim");options=FilterOptions(minLen=4,trimLeft=1))
        @test stats.reads_out == 0
    end
    @test expectedErrors(fill(UInt8(20),100)) ≈ 1.0
    @test reverseComplement("ACGTN") == "NACGT"
    @test_throws ArgumentError FastqRecord("r","AC",[40])
    @test_throws ArgumentError FastqRecord("r","AX",[40,40])
end

@testset "Quality trimming must preserve truncLen" begin
    mktempdir() do dir
        raw = joinpath(dir,"reads.fastq")
        # At truncLen, quality <= truncQ rejects; immediately after it is allowed.
        write(raw,"@early/1\nACGTACGT\n+\nIII#IIII\n" *
                  "@boundary/1\nACGTACGT\n+\nIIIII#II\n" *
                  "@after/1\nACGTACGT\n+\nIIIIII#I\n" *
                  "@good/1\nACGTACGT\n+\nIIIIIIII\n" *
                  "@short/1\nACGTA\n+\nIIIII\n")
        out = joinpath(dir,"fixed.fastq")
        opts = FilterOptions(truncLen=6,trimLeft=2,minLen=1,maxEE=Inf)
        stats = filterAndTrim(raw,out;options=opts)
        @test (stats.reads_in,stats.reads_out) == (5,2)
        reader = FastqReader(out)
        records = try collect(reader) finally close(reader) end
        @test getfield.(records,:id) == ["after/1","good/1"]
        @test getfield.(records,:sequence) == ["GTAC","GTAC"]
        @test all(r -> r.quality == fill(UInt8(40),4),records)
        free = joinpath(dir,"variable.fastq")
        stats = filterAndTrim(raw,free;options=FilterOptions(minLen=1,maxEE=Inf))
        @test stats.reads_out == 5
        reader = FastqReader(free)
        records = try collect(reader) finally close(reader) end
        @test length.(getfield.(records,:sequence)) == [3,5,6,8,5]
        reverse = joinpath(dir,"reverse.fastq")
        write(reverse,join(["@$(id)/2\nACGTACGT\n+\nIIIIIIII\n"
                           for id in ("early","boundary","after","good","short")]))
        for (index,f,r) in ((1,raw,reverse),(2,reverse,raw))
            of = joinpath(dir,"paired$(index)F.fastq")
            or = joinpath(dir,"paired$(index)R.fastq")
            stats = filterAndTrim(f,r,of,or;optionsF=opts,optionsR=opts)
            @test (stats.reads_in,stats.reads_out) == (5,2)
            for path in (of,or)
                reader = FastqReader(path)
                records = try collect(reader) finally close(reader) end
                @test DADA2.pair_id.(getfield.(records,:id)) == ["after","good"]
                @test all(r -> length(r.sequence) == 4,records)
            end
        end
    end
end

@testset "Alignment compared with independent dense DP" begin
    rng = MersenneTwister(12)
    for _ in 1:120
        a = String(rand(rng,UInt8.(['A','C','G','T']),rand(rng,1:45)))
        b = String(rand(rng,UInt8.(['A','C','G','T']),rand(rng,1:45)))
        al = alignSequences(a,b;band=3)
        @test al.score == dense_score(a,b)
        @test replace(al.reference,"-"=>"") == a
        @test replace(al.query,"-"=>"") == b
        @test count(==(0),al.mapping) <= al.indels
    end
    @test alignSequences("ACGT","ACGGT").indels == 1
    @test alignSequences("","AC").score == -16
    @test alignSequences("AC","").score == -16
    @test DADA2.kmer_distance(DADA2.kmers("AAAAAA"),DADA2.kmers("AAAAAA")) == 0
    @test DADA2.kmer_distance(DADA2.kmers("AAAAAA"),DADA2.kmers("CCCCCC")) == 1
end

@testset "Poisson tail and inference" begin
    @test abundancePvalue(1,1e-12) ≈ 1
    @test abundancePvalue(2,0.1) ≈ (1-exp(-0.1)*1.1)/(1-exp(-0.1)) rtol=1e-12
    @test abundancePvalue(2,1e-12) ≈ 5e-13 rtol=1e-10
    @test abundancePvalue(2,0.0) == 0
    @test abundancePvalue(1,0.0) == 1
    @test abundancePvalue(1,1.0;prior=true) ≈ 1-exp(-1)
    @test_throws ArgumentError abundancePvalue(0,1)
    parent = "ACGTACGTACGTACGTACGTACGTACGTACGT"
    variant = "ACGTACGTACGTTCGTACGTACGTACGTACGT"
    error = "ACGTACGTACGTACGTACGTACGTACGAACGT"
    d = fixture([parent,variant,error],[1000,100,1])
    result = dada(d)
    @test result.sequences == [parent,variant]
    @test result.counts == [1001,100]
    @test result.assignments == [1,2,1]
    @test sum(result.transitions) == sum(d.counts)*length(parent)
    @test result.converged
    @test sum(result.counts) == sum(d.counts)
    @test dada(d;threaded=false).assignments == result.assignments
    @test dada(d;maxClusters=1).converged == false
    @test isempty(dada(fixture(String[],Int[])).sequences)
    @test_throws ArgumentError dada(d,phredErrors(maxQ=30))
    @test_throws ArgumentError dada(Derep(["ACGT"],[1],[fill(38.5,4)]),phredErrors(maxQ=38))
    @test_throws ArgumentError ErrorModel(ones(16,41))
    learned = learnErrors([d];maxIter=3)
    @test learned.iterations <= 3
    @test all(isfinite,learned.errors.rates)
    @test all(q -> all(b -> sum(learned.errors.rates[4b+1:4b+4,q]) ≈ 1,0:3),1:94)
    @test sum(dada(d,learned).counts) == 1101
    @test_throws ArgumentError learnErrors(Derep[])
    zero_rates = copy(phredErrors().rates)
    zero_rates[1,1] = 0; zero_rates[2,1] = 0.5
    @test all(isfinite,learnErrors(d;initial=ErrorModel(zero_rates),maxIter=2).history)
    residual = fixture([repeat("A",32),repeat("C",32)],[100,1])
    strict = dada(residual)
    @test strict.rejected == 1
    @test strict.assignments == [1,0]
    @test dada(residual;omegaC=0).counts == [101]
    rng = MersenneTwister(3)
    seqs = [String(rand(rng,UInt8.(['A','C','G','T']),50)) for _ in 1:70]
    many = fixture(seqs,fill(2,70))
    @test dada(many;threaded=true).counts == dada(many;threaded=false).counts
    @test dada([d,d];threaded=true)[2].assignments == result.assignments
end

include("error_learning.jl")

@testset "Pair merging and community tables" begin
    f = "AAAACCCCGGGG"; r = reverseComplement("CCCCGGGGTTTT")
    fd = dereplicate([FastqRecord("r$i",f,fill(40,12)) for i in 1:3])
    rd = dereplicate([FastqRecord("r$i",r,fill(40,12)) for i in 1:3])
    merged = mergePairs(dada(fd),dada(rd);minOverlap=8)
    @test merged.sequences == ["AAAACCCCGGGGTTTT"]
    @test merged.counts == [3]
    @test (merged.accepted,merged.rejected) == (3,0)
    @test mergePairs(dada(fd),dada(rd);minOverlap=13).rejected == 3
    table = makeSequenceTable([merged,merged];samples=["A","B"])
    @test table.counts == reshape([3,3],2,1)
    p1 = "AAAAAAAAAAAAAAAA"; p2 = "CCCCCCCCCCCCCCCC"; chim = "AAAAAAAACCCCCCCC"
    @test isBimeraDenovo(chim,[p1,p2])
    @test !isBimeraDenovo(chim,[p1])
    @test !isBimeraDenovo(p1,[p1,p2])
    table = SequenceTable(["a","b"],[p1,p2,chim],[100 80 5; 100 80 4])
    clean = removeBimeraDenovo(table)
    @test clean.sequences == [p1,p2]
    @test removeBimeraDenovo(table;method=:pooled).counts == clean.counts
    rng = MersenneTwister(13)
    for _ in 1:40
        parents = [String(rand(rng,UInt8.(['A','C','G','T']),24)) for _ in 1:8]
        splitpos = rand(rng,8:16)
        query = parents[1][1:splitpos]*parents[2][splitpos+1:end]
        expected = any(i != j && parents[i] != query && parents[j] != query &&
            parents[i][1:cut] == query[1:cut] && parents[j][cut+1:end] == query[cut+1:end]
            for i in eachindex(parents),j in eachindex(parents),cut in 8:16)
        @test isBimeraDenovo(query,parents) == expected
    end
end

@testset "Taxonomy and FASTA" begin
    mktempdir() do dir
        ref = joinpath(dir,"tax.fa.gz")
        writeFasta(ref,[repeat("A",60),repeat("C",60)];ids=["Bacteria;Alpha;GenusA;","Bacteria;Beta;GenusB;"])
        @test readFasta(ref).sequences == [repeat("A",60),repeat("C",60)]
        biofasta = readFasta(ref;bio=true)
        @test biofasta.sequences isa Vector{LongDNA{4}}
        @test String.(biofasta.sequences) == [repeat("A",60),repeat("C",60)]
        taxonomy = assignTaxonomy([repeat("A",30),repeat("C",30)],ref;nBoot=20)
        @test taxonomy.taxa == ["Bacteria" "Alpha" "GenusA"; "Bacteria" "Beta" "GenusB"]
        @test all(==(100),taxonomy.bootstrap)
        @test taxonomy.bootstrap == assignTaxonomy([repeat("A",30),repeat("C",30)],ref;nBoot=20).bootstrap
        species = joinpath(dir,"species.fa")
        writeFasta(species,["AACCGGTT","AACCGGTT"];ids=["ref1 Genus first","ref2 Genus second"])
        @test ismissing(addSpecies(["AACCGGTT"],species)[1])
        @test addSpecies(["AACCGGTT"],species;allowMultiple=true) == ["Genus first/Genus second"]
        model = trainTaxonomy(ref)
        @test assignTaxonomy([repeat("A",30)],model;nBoot=20).taxa[1,3] == "GenusA"
        ambiguous = joinpath(dir,"ambiguous.fa")
        writeFasta(ambiguous,[repeat("A",20)*"NRY"*repeat("A",20)];ids=["Bacteria;Alpha;"])
        @test readFasta(ambiguous).sequences[1] == repeat("A",20)*"NRY"*repeat("A",20)
        @test assignTaxonomy([repeat("A",12)*"N"],ambiguous;nBoot=5).taxa[1,2] == "Alpha"
        @test all(ismissing,assignTaxonomy(["NNNNNNNNN"],model).taxa)
        @test_throws ArgumentError writeFasta(joinpath(dir,"empty.fa"),[""])
        write(joinpath(dir,"bad.fa"),">empty\n>next\nACGT\n")
        @test_throws ArgumentError readFasta(joinpath(dir,"bad.fa"))
        fragments = joinpath(dir,"fragments.fa")
        writeFasta(fragments,["AAAAGGCTTAGCGTCCCC","AAAAGGCTTAGAGTCCCC"];ids=["a Genus one","b Genus two"])
        smodel = trainSpecies(fragments)
        @test addSpecies(["GGCTTAGCGT"],smodel;tryRC=false) == ["Genus one"]
        @test addSpecies([reverseComplement("GGCTTAGCGT")],smodel) == ["Genus one"]
        @test ismissing(addSpecies(["TTTTTTTT"],smodel;tryRC=false)[1])
        @test ismissing(addSpecies(["AAA"],smodel;tryRC=false)[1])
        @test_throws ArgumentError addSpecies([""],smodel)
        bioquery=[LongDNA{4}("GGCTTAGCGT")]
        @test addSpecies(bioquery,smodel;tryRC=false)==["Genus one"]
        @test alignSequences(dna"AACCGGTT",dna"AACCGGTA").substitutions==1
        @test FastqRecord("bio",dna"ACGT",fill(40,4)).sequence=="ACGT"
        bioout=joinpath(dir,"bio.fa")
        writeFasta(bioout,[dna"ACGN"];ids=["bio"])
        @test readFasta(bioout;bio=true).sequences==[dna"ACGN"]
    end
    @test reverseComplement("ACGTRYSWKMBDHVN") == "NBDHVKMWSRYACGT"
    @test reverseComplement(dna"ACGTRYSWKMBDHVN") == dna"NBDHVKMWSRYACGT"
    @test isempty(DADA2.kmers("AAAANAAAA",5))
end

@testset "Bundled DADA2 paired-end fixture" begin
    include("upstream_smoke.jl")
end
