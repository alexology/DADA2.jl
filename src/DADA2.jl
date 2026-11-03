module DADA2

using CodecZlib
import BioSequences
import FASTX
using Distributions: Poisson, ccdf, logccdf
using Random
using LinearAlgebra: pinv, qr, svd, dot

export FastqRecord, FastqReader, Derep, ErrorModel, ErrorLearningResult, DadaResult, SequenceTable,
       FilterOptions, FilterStats, MergeResult, TaxonomyResult,
       filterAndTrim, derepFastq, dereplicate, dada, learnErrors, phredErrors, errorDiagnostics,
       mergePairs, makeSequenceTable, removeBimeraDenovo, isBimeraDenovo, isBimeraDenovoTable, isBimera,
       reverseComplement, expectedErrors, alignSequences,
       assignTaxonomy, trainTaxonomy, TaxonomyModel, addSpecies, trainSpecies,
       SpeciesModel, readFasta, writeFasta, abundancePvalue, isPhiX

include("types.jl")
include("fastq.jl")
include("alignment.jl")
include("core.jl")
include("rloess.jl")
include("inference.jl")
include("community.jl")
include("phix.jl")
include("merge.jl")
include("chimeras.jl")
include("taxonomy.jl")

end
