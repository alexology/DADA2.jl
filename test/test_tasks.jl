# Gate sample access inside the worker, without changing inference methods.
struct GatedSamples <: AbstractVector{Derep}
    values::Vector{Derep}
    started::Channel{Int}
    release::Vector{Channel{Nothing}}
    workers::Vector{Union{Nothing,Task}}
end
Base.size(x::GatedSamples) = size(x.values)
function Base.getindex(x::GatedSamples, i::Int)
    x.workers[i] = current_task()
    put!(x.started, i)
    take!(x.release[i])
    x.values[i]
end

@testset "Inference and raw recount finish sibling tasks after failure" begin
    if Threads.nthreads() == 1
        @test_skip "requires two threads"
    else
        mktempdir() do dir
            bad = fixture(["ACGT"], [2]; q = 40)
            good = fixture(["TGCA"], [2]; q = 10)
            records = [bad, good]
            results = dada(records; threaded = false)
            path = joinpath(dir, "good.fastq")
            write(path, "@one\nTGCA\n+\n++++\n@two\nTGCA\n+\n++++\n")
            sources = [DADA2.ErrorTrainingSource(joinpath(dir, "missing.fastq"), 2, 8),
                       DADA2.ErrorTrainingSource(path, 2, 8)]
            for raw in (false, true)
                gated = GatedSamples(records, Channel{Int}(2),
                                     [Channel{Nothing}(1) for _ in 1:2],
                                     Union{Nothing,Task}[nothing, nothing])
                caches = DADA2.InferenceCache.(records)
                caller = @async try
                    if raw
                        DADA2.raw_transitions!(zeros(Int, 16, 94), gated, results, sources, caches)
                    else
                        dada(gated, phredErrors(maxQ = 20); threaded = true)
                    end
                    nothing
                catch e
                    e
                end
                try
                    for _ in 1:2
                        @test timedwait(() -> isready(gated.started), 10) == :ok
                        take!(gated.started)
                    end
                    put!(gated.release[1], nothing)
                    @test timedwait(() -> istaskdone(gated.workers[1]), 10) == :ok
                    @test timedwait(() -> istaskdone(caller), 0.1) == :timed_out
                finally
                    # Release both workers even when a preceding assertion fails.
                    for channel in gated.release
                        isready(channel) || put!(channel, nothing)
                    end
                end
                @test timedwait(() -> istaskdone(caller), 10) == :ok
                @test fetch(caller) isa TaskFailedException
                @test all(istaskdone, gated.workers)
            end
        end
    end
end
