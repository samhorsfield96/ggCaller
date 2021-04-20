#imports
import argparse
import sys
import json
from Bio.Seq import Seq
import ggCaller_cpp
from ggCaller.graph_traversal import *
from multiprocessing import Pool
from functools import partial
from memory_profiler import profile


def get_options():
    description = 'Generates ORFs from a Bifrost graph.'
    parser = argparse.ArgumentParser(description=description,
                                     prog='ggcaller')

    IO = parser.add_argument_group('Input/options.out')
    IO.add_argument('--graph',
                    default=None,
                    help='Bifrost GFA file generated by Bifrost build. ')
    IO.add_argument('--colours',
                    default=None,
                    help='Bifrost colours file generated by Bifrost build.  ')
    IO.add_argument('--not-ref',
                    action="store_false",
                    default=True,
                    help='If using existing graph, was not graph built exclusively with assembled genomes.  '
                         '[Default = False] ')
    IO.add_argument('--refs',
                    default=None,
                    help='List of reference genomes (one file path per line). ')
    IO.add_argument('--reads',
                    default=None,
                    help='List of reference genomes (one file path per line). ')
    IO.add_argument('--codons',
                    default=None,
                    help='JSON file containing start and stop codon sequences. ')
    IO.add_argument('--kmer',
                    type=int,
                    default=31,
                    help='K-mer size used in Bifrost build (bp). '
                         '[Default = 31] ')
    IO.add_argument('--path',
                    type=int,
                    default=10000,
                    help='Maximum path length during traversal (bp). '
                         '[Default = 10000] ')
    IO.add_argument('--orf',
                    type=int,
                    default=90,
                    help='Minimum ORF length to return (bp). '
                         '[Default = 90] ')
    IO.add_argument('--maxoverlap',
                    type=int,
                    default=60,
                    help='Maximum overlap allowed between overlapping ORFs. '
                         '[Default = 60] ')
    IO.add_argument('--min-path-score',
                    type=int,
                    default=100,
                    help='Minimum total score for a path of ORFs to be returned. '
                         '[Default = 100] ')
    IO.add_argument('--min-orf-score',
                    type=int,
                    default=100,
                    help='Minimum individual score for an ORF to be returned. '
                         '[Default = 100] ')
    IO.add_argument('--no-filter',
                    action="store_true",
                    default=False,
                    help='Do not filter ORF calls using Balrog. Will return all ORF calls. '
                         '[Default = False] ')
    IO.add_argument('--no-write-idx',
                    action="store_false",
                    default=True,
                    help='Do not write FMIndexes to file. '
                         '[Default = False] ')
    IO.add_argument('--no-write-graph',
                    action="store_false",
                    default=True,
                    help='Do not write Bifrost GFA and colours to file. '
                         '[Default = False] ')
    IO.add_argument('--repeat',
                    action="store_true",
                    default=False,
                    help='Enable traversal of nodes mulitple times. '
                         '[Default = False] ')
    IO.add_argument('--threads',
                    type=int,
                    default=1,
                    help='Number of threads for FMIndexing '
                         '[Default = 1] ')
    IO.add_argument('--out',
                    default='calls.fasta',
                    help='options.out FASTA file containing ORF sequences. ')
    return parser.parse_args()

def main():
    options = get_options()

    # parse command line arguments
    options.out = options.out

    # define start/stop codons
    if options.codons != None:
        with open(options.codons, "r") as json_file:
            try:
                data = json.load(json_file)
                start_codons = data["codons"]["start"]
                stop_codons_for = data["codons"]["stop"]
                stop_codons_rev = [str((Seq(i)).reverse_complement()) for i in stop_codons_for]
            except:
                print("Please specify codons in the format shown in codons.json.")
                sys.exit(1)
    else:
        start_codons = ["ATG", "GTG", "TTG"]
        stop_codons_for = ["TAA", "TGA", "TAG"]
        stop_codons_rev = ["TTA", "TCA", "CTA"]

    # if build graph specified, build graph and then call ORFs
    if (options.graph != None) and (options.colours != None) and (options.refs == None) and (options.reads == None):
        graph_tuple = ggCaller_cpp.index_existing(options.graph, options.colours, stop_codons_for, stop_codons_rev,
                                                  options.threads, True)
    # if refs file specified for building
    elif (options.graph == None) and (options.colours == None) and (options.refs != None) and (options.reads == None):
        graph_tuple = ggCaller_cpp.index_build(options.refs, options.kmer, stop_codons_for, stop_codons_rev,
                                               options.threads, True, options.no_write_graph)
    # if reads file specified for building
    elif (options.graph == None) and (options.colours == None) and (options.refs == None) and (options.reads != None):
        graph_tuple = ggCaller_cpp.index_build(options.reads, options.kmer, stop_codons_for, stop_codons_rev,
                                               options.threads, False, options.no_write_graph)
    # if both reads and refs file specified for building
    elif (options.graph == None) and (options.colours == None) and (options.refs != None) and (options.reads != None):
        graph_tuple = ggCaller_cpp.index_build(options.refs, options.kmer, stop_codons_for, stop_codons_rev,
                                               options.threads, False, options.no_write_graph, options.reads)
    else:
        print("Error: incorrect number of input files specified. Please only specify the below combinations:\n"
              "- Bifrost GFA and Bifrost colours file\n"
              "- List of reference files\n"
              "- List of read files\n"
              "- A list of reference files and a list of read files.")
        sys.exit(1)

    # unpack ORF pair into overlap dictionary and list for gene scoring
    graph_vector, node_colour_vector, input_colours, nb_colours, overlap = graph_tuple

    if not options.no_filter:
        print("Loading gene scoring models...")
        model, model_tis, aa_kmer_set = load_models(options.threads)
    else:
        model, model_tis, aa_kmer_set = None, None, None

    # run run_calculate_ORFs with multithreading
    # create list for high scoring ORFs to return
    true_genes = {}
    print("Generating high scoring ORF calls...")

    # Turn gt threading off
    if gt.openmp_enabled():
        gt.openmp_set_num_threads(1)

    with Pool(processes=options.threads) as executor:
        for colour_ID, col_true_genes in executor.map(
                partial(run_calculate_ORFs, graph_vector=graph_vector, repeat=options.repeat, overlap=overlap,
                        max_path_length=options.path,
                        is_ref=options.not_ref, no_filter=options.no_filter, stop_codons_for=stop_codons_for,
                        start_codons=start_codons,
                        min_ORF_length=options.orf,
                        max_ORF_overlap=options.maxoverlap, minimum_ORF_score=options.min_orf_score,
                        minimum_path_score=options.min_path_score, write_idx=options.no_write_idx,
                        input_colours=input_colours, nb_colours=nb_colours, model=model, model_tis=model_tis,
                        aa_kmer_set=aa_kmer_set),
                enumerate(node_colour_vector)):
            # iterate over entries in col_true_genes to generate the sequences
            for ORFNodeVector in col_true_genes:
                gene = generate_seq(graph_vector, ORFNodeVector[0], ORFNodeVector[1], overlap)
                if gene not in true_genes:
                    # create tuple to hold ORF sequence, colours and graph traversal information
                    empty_colours_list = ["0"] * nb_colours
                    true_genes[gene] = (empty_colours_list, ORFNodeVector)
                # update colours with current colour_ID
                true_genes[gene][0][colour_ID] = "1"

    print("Generating fasta file of gene calls...")
    # print output to file
    ORF_count = 1
    with open(options.out, "w") as f:
        for gene, info_pair in true_genes.items():
            colour_str = "".join(info_pair[0])
            f.write(">" + str(ORF_count) + "_" + str(colour_str) + "\n" + gene + "\n")
            ORF_count += 1

    print("Finished.")

    sys.exit(0)
