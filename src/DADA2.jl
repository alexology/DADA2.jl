module DADA2

using CodecZlib
import BioSequences
import FASTX
using Distributions: Poisson, logccdf
using Random
using LinearAlgebra: pinv

export FastqRecord, FastqReader, Derep, ErrorModel, ErrorLearningResult, DadaResult, SequenceTable,
       FilterOptions, FilterStats, MergeResult, TaxonomyResult,
       filterAndTrim, derepFastq, dereplicate, dada, learnErrors, phredErrors, errorDiagnostics,
       mergePairs, makeSequenceTable, removeBimeraDenovo, isBimeraDenovo,
       reverseComplement, expectedErrors, alignSequences,
       assignTaxonomy, trainTaxonomy, TaxonomyModel, addSpecies, trainSpecies,
       SpeciesModel, readFasta, writeFasta, abundancePvalue

include("types.jl")
include("fastq.jl")
include("alignment.jl")
include("inference.jl")
include("community.jl")
include("taxonomy.jl")

end
