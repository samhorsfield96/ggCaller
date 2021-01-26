#imports
import argparse
import sys
import json
from Bio.Seq import Seq
import ggCaller_cpp
#from balrog.main import *

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
    IO.add_argument('--no-write-idx',
                    action="store_false",
                    default=True,
                    help='Do not write FMIndexes to file. '
                         '[Default = True] ')
    IO.add_argument('--no-write-graph',
                    action="store_false",
                    default=True,
                    help='Do not write Bifrost GFA and colours to file. '
                         '[Default = True] ')
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

def main():
    options = get_options()

    # parse command line arguments
    graph_file = options.graph
    colours_file = options.colours
    is_ref = bool(options.not_ref)
    refs_file = options.refs
    reads_file = options.reads
    codons = options.codons
    ksize = int(options.kmer)
    max_path_length = int(options.path)
    min_ORF_length = int(options.orf)
    write_idx = bool(options.no_write_idx)
    write_graph = bool(options.no_write_graph)
    repeat = bool(options.repeat)
    num_threads = int(options.threads)
    output = options.out

    # define start/stop codons
    if codons != None:
        with open(codons, "r") as json_file:
            data = json.load(json_file)
            start_codons = data["codons"]["start"]
            stop_codon_for = data["codons"]["stop"]
            stop_codon_rev = [str((Seq(i)).reverse_complement()) for i in stop_codon_for]
    else:
        start_codons = ["ATG", "GTG", "TTG"]
        stop_codon_for = ["TAA", "TGA", "TAG"]
        stop_codon_rev = ["TTA", "TCA", "CTA"]

    # if build graph specified, build graph and then call ORFs
    if (graph_file != None) and (colours_file != None) and (refs_file == None) and (reads_file == None):
        ggCaller_cpp.call_genes_existing(graph_file, colours_file, output, start_codons, stop_codon_for, stop_codon_rev, num_threads, is_ref, write_idx, repeat, max_path_length, min_ORF_length)
    # if refs file specified for building
    elif (graph_file == None) and (colours_file == None) and (refs_file != None) and (reads_file == None):
        ggCaller_cpp.call_genes_build(refs_file, ksize, output, start_codons, stop_codon_for, stop_codon_rev, num_threads, True, write_idx, repeat, write_graph, max_path_length, min_ORF_length)
    # if reads file specified for building
    elif (graph_file == None) and (colours_file == None) and (refs_file == None) and (reads_file != None):
        ggCaller_cpp.call_genes_build(reads_file, ksize, output, start_codons, stop_codon_for, stop_codon_rev, num_threads, False, write_idx, repeat, write_graph, max_path_length, min_ORF_length)
    # if both reads and refs file specified for building
    elif (graph_file == None) and (colours_file == None) and (refs_file != None) and (reads_file != None):
        ggCaller_cpp.call_genes_build(refs_file, ksize, output, start_codons, stop_codon_for, stop_codon_rev, num_threads, False, write_idx, repeat, write_graph, max_path_length, min_ORF_length, reads_file)
    else:
        print("Error: incorrect number of input files specified. Please only specify the below combinations:\n"
              "- Bifrost GFA and Bifrost colours file\n"
              "- List of reference files\n"
              "- List of read files\n"
              "- A list of reference files and a list of read files.")

    sys.exit(0)

if __name__ == '__main__':
    main()
