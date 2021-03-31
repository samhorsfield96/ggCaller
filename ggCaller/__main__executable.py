# imports
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

    IO = parser.add_argument_group('Input/output')
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
                    default=31,
                    help='K-mer size used in Bifrost build (bp). '
                         '[Default = 31] ')
    IO.add_argument('--path',
                    default=10000,
                    help='Maximum path length during traversal (bp). '
                         '[Default = 10000] ')
    IO.add_argument('--orf',
                    default=90,
                    help='Minimum ORF length to return (bp). '
                         '[Default = 90] ')
    IO.add_argument('--maxoverlap',
                    default=60,
                    help='Maximum overlap allowed between overlapping ORFs. '
                         '[Default = 60] ')
    IO.add_argument('--min-path-score',
                    default=100,
                    help='Minimum total score for a path of ORFs to be returned. '
                         '[Default = 100] ')
    IO.add_argument('--min-orf-score',
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
                    default=1,
                    help='Number of threads for FMIndexing '
                         '[Default = 1] ')
    IO.add_argument('--out',
                    default='calls.fasta',
                    help='Output FASTA file containing ORF sequences. ')
    return parser.parse_args()


@profile
def main():
    # options = get_options()
    #
    # # parse command line arguments
    # graph_file = options.graph
    # colours_file = options.colours
    # is_ref = bool(options.not_ref)
    # refs_file = options.refs
    # reads_file = options.reads
    # codons = options.codons
    # ksize = int(options.kmer)
    # max_path_length = int(options.path)
    # min_ORF_length = int(options.orf)
    # max_ORF_overlap = int(options.maxoverlap)
    # minimum_path_score = int(options.min_path_score)
    # minimum_ORF_score = int(options.min_orf_score)
    # no_filter = bool(options.no_filter)
    # write_idx = bool(options.no_write_idx)
    # write_graph = bool(options.no_write_graph)
    # repeat = bool(options.repeat)
    # num_threads = int(options.threads)
    # output = options.out
    #
    # # define start/stop codons
    # if codons != None:
    #     with open(codons, "r") as json_file:
    #         try:
    #             data = json.load(json_file)
    #             start_codons = data["codons"]["start"]
    #             stop_codon_for = data["codons"]["stop"]
    #             stop_codon_rev = [str((Seq(i)).reverse_complement()) for i in stop_codon_for]
    #         except:
    #             print("Please specify codons in the format shown in codons.json.")
    #             sys.exit(1)
    # else:
    #     start_codons = ["ATG", "GTG", "TTG"]
    #     stop_codon_for = ["TAA", "TGA", "TAG"]
    #     stop_codon_rev = ["TTA", "TCA", "CTA"]

    start_codons = ["ATG", "GTG", "TTG"]
    stop_codon_for = ["TAA", "TGA", "TAG"]
    stop_codon_rev = ["TTA", "TCA", "CTA"]

    output = "/mnt/c/Users/sth19/PycharmProjects/Genome_Graph_project/ggCaller/group3_capsular_fa_list_integer_paths.fasta"
    # set mimimum path score
    minimum_path_score = 100
    minimum_ORF_score = 100
    no_filter = False
    repeat = False
    max_path_length = 10000
    is_ref = True
    min_ORF_length = 90
    max_ORF_overlap = 60
    write_idx = True

    num_threads = 4

    graph_tuple = ggCaller_cpp.index_existing(
        "/mnt/c/Users/sth19/PycharmProjects/Genome_Graph_project/ggCaller/data/group3_capsular_fa_list.gfa",
        "/mnt/c/Users/sth19/PycharmProjects/Genome_Graph_project/ggCaller/data/group3_capsular_fa_list.bfg_colors",
        start_codons,
        stop_codon_for, stop_codon_rev, num_threads, True)

    # # if build graph specified, build graph and then call ORFs
    # if (graph_file != None) and (colours_file != None) and (refs_file == None) and (reads_file == None):
    #     called_ORF_tuple = ggCaller_cpp.call_genes_existing(graph_file, colours_file, start_codons, stop_codon_for,
    #                                                         stop_codon_rev, num_threads, is_ref, write_idx, repeat,
    #                                                         no_filter, max_path_length, min_ORF_length, max_ORF_overlap)
    # # if refs file specified for building
    # elif (graph_file == None) and (colours_file == None) and (refs_file != None) and (reads_file == None):
    #     called_ORF_tuple = ggCaller_cpp.call_genes_build(refs_file, ksize, start_codons, stop_codon_for, stop_codon_rev,
    #                                                      num_threads, True, write_idx, repeat, write_graph, no_filter,
    #                                                      max_path_length, min_ORF_length, max_ORF_overlap)
    # # if reads file specified for building
    # elif (graph_file == None) and (colours_file == None) and (refs_file == None) and (reads_file != None):
    #     called_ORF_tuple = ggCaller_cpp.call_genes_build(reads_file, ksize, start_codons, stop_codon_for,
    #                                                      stop_codon_rev, num_threads, False, write_idx, repeat,
    #                                                      write_graph, no_filter, max_path_length, min_ORF_length,
    #                                                      max_ORF_overlap)
    # # if both reads and refs file specified for building
    # elif (graph_file == None) and (colours_file == None) and (refs_file != None) and (reads_file != None):
    #     called_ORF_tuple = ggCaller_cpp.call_genes_build(refs_file, ksize, start_codons, stop_codon_for, stop_codon_rev,
    #                                                      num_threads, False, write_idx, repeat, write_graph, no_filter,
    #                                                      max_path_length, min_ORF_length, max_ORF_overlap, reads_file)
    # else:
    #     print("Error: incorrect number of input files specified. Please only specify the below combinations:\n"
    #           "- Bifrost GFA and Bifrost colours file\n"
    #           "- List of reference files\n"
    #           "- List of read files\n"
    #           "- A list of reference files and a list of read files.")
    #     sys.exit(1)

    # unpack ORF pair into overlap dictionary and list for gene scoring
    graph_vector, node_colour_vector, input_colours, nb_colours, overlap = graph_tuple

    # create list for high scoring ORFs to return
    true_genes = {}

    # calculate ORFs within graph
    for colour_ID, node_set in enumerate(node_colour_vector):
        ORF_overlap_dict, ORF_vector = ggCaller_cpp.calculate_ORFs(graph_vector, colour_ID, node_set, repeat,
                                                                   overlap, max_path_length, is_ref, no_filter,
                                                                   stop_codon_for, start_codons, min_ORF_length,
                                                                   max_ORF_overlap, write_idx, input_colours[colour_ID])

        # calculate scores for genes
        ORF_score_dict = score_genes(ORF_vector, graph_vector, minimum_ORF_score, overlap, num_threads)

        # determine
        high_scoring_ORFs = call_true_genes(ORF_score_dict, ORF_overlap_dict, minimum_path_score)

    # if no filter specified, generate ORF sequences and combine colours of identical ORFS
    # if no_filter == True:
    #     for colour, gene_set in ORF_colour_ID_map.items():
    #         for ORF_ID in gene_set.items():
    #             ORFNodeVector = full_ORF_dict[ORF_ID]
    #             gene = generate_seq(unitig_map, ORFNodeVector[0], ORFNodeVector[1], ORFNodeVector[2], overlap)
    #             if gene not in true_genes:
    #                 # create string of zeros, make nth colour 1
    #                 true_genes[gene] = ["0"] * nb_colours
    #                 true_genes[gene][colour] = "1"
    #             else:
    #                 # make nth colour 1
    #                 true_genes[gene][colour] = "1"
    # else:
    #     # generate gene scores using Balrog
    #     full_ORF_dict, ORF_score_dict = score_genes(full_ORF_dict, minimum_ORF_score, unitig_map, overlap, num_threads)
    #
    #     print("Generating highest scoring gene paths...")
    #
    #     for colour_ORF_tuple in ORF_colour_ID_map.items():
    #         colour, high_scoring_ORFs = call_true_genes(colour_ORF_tuple, minimum_path_score, ORF_score_dict,
    #                                                     ORF_overlap_dict)
    #         # merge any matching genes which had differing TIS but were called together
    #         for ORF_ID in high_scoring_ORFs:
    #             # parse out gene string from full_ORF_dict, generate sequence
    #             ORFNodeVector = full_ORF_dict[ORF_ID]
    #             gene = generate_seq(unitig_map, ORFNodeVector[0], ORFNodeVector[1], overlap)
    #             if gene not in true_genes:
    #                 # create string of zeros, make nth colour 1
    #                 true_genes[gene] = ["0"] * nb_colours
    #                 true_genes[gene][colour] = "1"
    #             else:
    #                 # make nth colour 1
    #                 true_genes[gene][colour] = "1"

    # with Pool(processes=num_threads) as pool:
    #     for colour, high_scoring_ORFs in pool.map(
    #             partial(call_true_genes, minimum_path_score=minimum_path_score,
    #                     ORF_score_dict=ORF_score_dict, ORF_overlap_dict=ORF_overlap_dict),
    #             ORF_colour_ID_map.items()):
    #         # remove TIS, and merge any matching genes which had differing TIS but were called together
    #         for ORF_ID in high_scoring_ORFs:
    #             # parse out gene string from full_ORF_dict
    #             gene = str(full_ORF_dict[ORF_ID][0])[16:]
    #             if gene not in true_genes:
    #                 # create string of zeros, make nth colour 1
    #                 true_genes[gene] = ["0"] * nb_colours
    #                 true_genes[gene][colour] = "1"
    #             else:
    #                 # make nth colour 1
    #                 true_genes[gene][colour] = "1"

    # print("Generating fasta file of gene calls...")
    # # print output to file
    # ORF_count = 1
    # with open(output, "w") as f:
    #     for gene, colour in true_genes.items():
    #         colour_str = "".join(colour)
    #         f.write(">" + str(ORF_count) + "_" + str(colour_str) + "\n" + gene + "\n")
    #         ORF_count += 1

    print("Finished.")

    sys.exit(0)


if __name__ == '__main__':
    main()
