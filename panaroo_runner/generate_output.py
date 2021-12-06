from functools import partial
import shutil
from ggCaller.shared_memory import *
import ggCaller_cpp
from ggCaller import __version__
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from Bio import AlignIO
import itertools as iter
from intbitset import intbitset
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from scipy.optimize import curve_fit
from scipy.interpolate import make_interp_spline, BSpline
from random import shuffle
from uncertainties import ufloat
from .generate_alignments import *


def back_translate_dir(high_scoring_ORFs, isolate_names, annotation_dir, overlap, shd_arr_tup, pool):
    # iterate over all files in annotation directory multithreaded
    pool.map(partial(back_translate, isolate_names=isolate_names, shd_arr_tup=shd_arr_tup,
                     high_scoring_ORFs=high_scoring_ORFs, overlap=overlap, annotation_dir=annotation_dir),
             os.listdir(annotation_dir))

    return True


def back_translate(file, annotation_dir, shd_arr_tup, high_scoring_ORFs, isolate_names, overlap):
    # load shared memory items
    existing_shm = shared_memory.SharedMemory(name=shd_arr_tup.name)
    shd_arr = np.ndarray(shd_arr_tup.shape, dtype=shd_arr_tup.dtype, buffer=existing_shm.buf)

    file = os.path.join(annotation_dir, file)

    output_sequences = []
    with open(file, "r") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            record_ID = record.id
            protein = str(record.seq)
            mem = int(record_ID.split("_")[0])
            ORF_ID = int(record_ID.split("_")[-1])

            # parse DNA sequence
            ORFNodeVector = high_scoring_ORFs[mem][ORF_ID]
            if ORF_ID < 0:
                dna = ORFNodeVector[5]
            else:
                dna = shd_arr[0].generate_sequence(ORFNodeVector[0], ORFNodeVector[1], overlap)

            # # determine if Ns have been added into alignment. If so, need to account for unaligned stop codon
            # if dna[-1] == "N":
            #     dna += "---"
            #
            # # add on stop codon if not present in aa sequence
            # if protein[-1] != "*":
            #     protein += "*"

            # back translate sequence
            aligned_dna = ""
            dna_idx = 0
            for aa_idx in range(len(protein)):
                if protein[aa_idx] == "-":
                    aligned_dna += "---"
                else:
                    to_add = dna[dna_idx: dna_idx + 3]
                    # remove any Ns present in DNA sequence
                    if "N" in to_add:
                        to_add.replace("N", "-")
                    aligned_dna += to_add
                    dna_idx += 3

            id = isolate_names[mem] + "_" + str(ORF_ID).zfill(5)
            output_sequences.append(SeqRecord(Seq(aligned_dna), id=id, description=""))

    # overwrite existing alignment file
    output_sequences = (x for x in output_sequences)
    SeqIO.write(output_sequences, file, 'fasta')

    return


def print_ORF_calls(high_scoring_ORFs, outfile, input_colours, overlap, DBG, truncation_threshold=0, G=None):
    isolate_names = [
        os.path.splitext(os.path.basename(x))[0] for x in input_colours
    ]
    if G == None:
        with open(outfile, "w") as f:
            for colour, gene_dict in high_scoring_ORFs.items():
                for ORF_ID, ORFNodeVector in gene_dict.items():
                    gene = DBG.generate_sequence(ORFNodeVector[0], ORFNodeVector[1], overlap)
                    f.write(">" + isolate_names[colour] + "_" + str(ORF_ID).zfill(5) + "\n" + gene + "\n")
    else:
        with open(outfile, "w") as f:
            for node in G.nodes():
                node_annotation = G.nodes[node]["description"]
                centroid_sequence_ids = set(G.nodes[node]["centroid"])
                length_centroid = G.nodes[node]['longCentroidID'][0]
                for i in range(0, len(G.nodes[node]["centroid"])):
                    colour = int(G.nodes[node]["centroid"][i].split("_")[0])
                    ORF_ID = int(G.nodes[node]["centroid"][i].split("_")[-1])
                    ORF_info = high_scoring_ORFs[colour][ORF_ID]
                    gene = G.nodes[node]["dna"][i]
                    ORF_len = ORF_info[2]
                    gene_annotation = node_annotation
                    if ORF_len < (length_centroid * truncation_threshold) or (
                            ORF_ID < 0 and (ORF_info[3] is True
                                            or ORF_info[
                                                2] % 3 != 0)):
                        gene_annotation += "(potential psuedogene)"
                    if gene_annotation != "":
                        f.write(">" + isolate_names[colour] + "_" + str(ORF_ID).zfill(
                            5) + " " + gene_annotation + "\n" + gene + "\n")
                    else:
                        f.write(">" + isolate_names[colour] + "_" + str(ORF_ID).zfill(5) + "\n" + gene + "\n")
                for sid in G.nodes[node]['seqIDs']:
                    if sid not in centroid_sequence_ids:
                        colour = int(sid.split("_")[0])
                        ORF_ID = int(sid.split("_")[-1])
                        ORF_info = high_scoring_ORFs[colour][ORF_ID]
                        gene = DBG.generate_sequence(ORF_info[0], ORF_info[1], overlap)
                        ORF_len = ORF_info[2]
                        gene_annotation = node_annotation
                        if ORF_len < (length_centroid * truncation_threshold) or (
                                ORF_ID < 0 and (ORF_info[3] is True
                                                or ORF_info[
                                                    2] % 3 != 0)):
                            gene_annotation += "(potential psuedogene)"
                        if gene_annotation != "":
                            f.write(">" + isolate_names[colour] + "_" + str(ORF_ID).zfill(
                                5) + " " + gene_annotation + "\n" + gene + "\n")
                        else:
                            f.write(">" + isolate_names[colour] + "_" + str(ORF_ID).zfill(5) + "\n" + gene + "\n")

    return


def generate_GFF(input_colours, isolate_names, contig_annotation, output_dir):
    # create directory for gffs
    GFF_dir = os.path.join(output_dir, "GFF")
    if not os.path.exists(GFF_dir):
        os.mkdir(GFF_dir)
    # make sure trailing forward slash is present
    GFF_dir = os.path.join(GFF_dir, "")

    # iterate over colours in contig_annotation, writing output file
    for colour in range(len(input_colours)):
        # iterate over the entries in input_colors
        with open(input_colours[colour]) as handle:
            gff_record_list = []
            record_id = 1
            for record in SeqIO.parse(handle, "fasta"):
                # sort records based on first position
                contig_annotation[colour][record_id].sort(key=lambda x: x[1][0][1][0])
                gff_record = record
                gff_record.features = []
                entry_ID = 1
                for entry in contig_annotation[colour][record_id]:
                    qualifiers = {
                        "source": "ggCaller:" + __version__,
                        "ID": isolate_names[colour] + "_" + str(entry[0]).zfill(5),
                        "inference": entry[2][0],
                        "score": entry[2][2],
                        "annotation": [entry[2][3].replace(",", " ").replace("|", "-")
                                           .replace("=", "-").replace("(", "").replace(")", "")],
                    }
                    feature = SeqFeature(
                        FeatureLocation(entry[1][0][1][0], entry[1][0][1][1]),
                        type="gene", strand=int(entry[1][1]), qualifiers=qualifiers
                    )
                    gff_record.features.append(feature)
                    entry_ID += 1
                gff_record_list.append(gff_record)
                record_id += 1

        # write to GFF file
        gff_record_list = (x for x in gff_record_list)
        outfile = GFF_dir + isolate_names[colour] + ".gff"

        with open(outfile, "w") as out_handle:
            GFF.write(gff_record_list, out_handle)

    return

def output_aa_sequence(node_pair):
    # unpack node_pair
    node = node_pair[1]

    ref_output_sequences = []

    # iterate over centroids to generate fasta files
    for i in range(0, len(node["centroid"])):
        name = str(node_pair[0]) + ";" + node["centroid"][i]
        ref_output_sequences.append(SeqRecord(Seq(node["protein"][i]), id=name, description=""))

    return ref_output_sequences


def output_alignment_sequence(node_pair, temp_directory, outdir, shd_arr_tup, high_scoring_ORFs, overlap,
                              ref_aln, ignore_pseduogenes, truncation_threshold):
    # load shared memory items
    existing_shm = shared_memory.SharedMemory(name=shd_arr_tup.name)
    shd_arr = np.ndarray(shd_arr_tup.shape, dtype=shd_arr_tup.dtype, buffer=existing_shm.buf)

    # unpack node_pair
    node = node_pair[1]

    # Counter for the number of sequences to
    seq_no = 0

    # get outname for reference alignment file
    ref_outname = None

    # determine sequences to aligned
    if ignore_pseduogenes:
        length_centroid = node['longCentroidID'][0]
        # identify pseudogenes based on truncation threshold, and if annotated as
        # having premature stop or not multiple of 3 long
        sequence_ids = [x for x in node["seqIDs"] if
                        not ((high_scoring_ORFs[int(x.split("_")[0])][int(x.split("_")[-1])][2]
                              < length_centroid * truncation_threshold) or (int(x.split("_")[-1]) < 0 and
                                                                            (high_scoring_ORFs[int(x.split("_")[0])][
                                                                                 int(x.split("_")[-1])][3] is True or
                                                                             high_scoring_ORFs[int(x.split("_")[0])][
                                                                                 int(x.split("_")[-1])][2] % 3 != 0
                                                                             )))]
        centroid_sequence_ids = set(
            [x for x in node["centroid"] if
             not ((high_scoring_ORFs[int(x.split("_")[0])][int(x.split("_")[-1])][2]
                   < length_centroid * truncation_threshold) or (int(x.split("_")[-1]) < 0 and
                                                                 (high_scoring_ORFs[int(x.split("_")[0])][
                                                                      int(x.split("_")[-1])][3] is True or
                                                                  high_scoring_ORFs[int(x.split("_")[0])][
                                                                      int(x.split("_")[-1])][2] % 3 != 0
                                                                  )))])
    else:
        sequence_ids = node["seqIDs"]
        centroid_sequence_ids = set(node["centroid"])

    # if reference-guided alignment being done, separate centroids and other sequences
    if ref_aln:
        ref_output_sequences = [SeqRecord(Seq(node["protein"][i]), id=node["centroid"][i], description="") for i in
                                range(len(node["centroid"])) if node["centroid"][i] in centroid_sequence_ids]
        centroid_no = len(ref_output_sequences)
        # Put gene of interest sequences in a generator, with corrected isolate names
        ref_output_sequences_gen = (x for x in ref_output_sequences)
        if len(sequence_ids) == centroid_no and centroid_no == 1:
            # If only one sequence, output it to aligned directory and break
            # if no other sequences, then just output with no alignment
            ref_outname = outdir + "/aligned_gene_sequences/" + node["name"] + ".aln.fas"
            if len(ref_outname) >= 248:
                ref_outname = ref_outname[:248] + ".fasta"
            SeqIO.write(ref_output_sequences_gen, ref_outname, 'fasta')
            return None, None
        else:
            # if centroid is on it's own, give name aln for aligned, otherwise ref
            if centroid_no > 1:
                ref_outname = temp_directory + node["name"] + "_ref.fasta"
            else:
                ref_outname = temp_directory + node["name"] + "_ref.aln.fas"
        if len(ref_outname) >= 248:
            ref_outname = ref_outname[:248] + ".fasta"
        SeqIO.write(ref_output_sequences_gen, ref_outname, 'fasta')

    # Look for gene sequences among all genes (from memory)
    output_sequences = []
    for seq in sequence_ids:
        # check if reference true and seqID is centroid. Is so, pass
        if ref_aln and seq in centroid_sequence_ids:
            continue
        member = int(seq.split('_')[0])
        ORF_ID = int(seq.split('_')[-1])
        # generate protein sequence. If refound, add found protein sequence from find_missing.py
        ORFNodeVector = high_scoring_ORFs[member][ORF_ID]
        if ORF_ID < 0:
            protein = ORFNodeVector[4]
        else:
            protein = str(Seq(shd_arr[0].generate_sequence(ORFNodeVector[0], ORFNodeVector[1], overlap)).translate())

        output_sequences.append(
            SeqRecord(Seq(protein), id=seq, description=""))
        seq_no += 1
    # Put gene of interest sequences in a generator, with corrected isolate names
    output_sequences = (x for x in output_sequences)
    # set filename to gene name, if more than one sequence to be aligned
    if (not ref_aln and seq_no > 1) or (ref_aln and seq_no > 0):
        outname = temp_directory + node["name"] + ".fasta"
    else:
        # if number sequences is 0, do not output
        if seq_no > 0:
            # If only one sequence, output it to aligned directory and break
            outname = outdir + "/aligned_gene_sequences/" + node["name"] + ".aln.fas"
            if len(outname) >= 248:
                outname = outname[:248] + ".fasta"
            SeqIO.write(output_sequences, outname, 'fasta')
        return None, ref_outname

    # Write them to disk
    if len(outname) >= 248:
        outname = outname[:248] + ".fasta"
    SeqIO.write(output_sequences, outname, 'fasta')
    return outname, ref_outname


def power_law(x, a, b):
    return a * np.power(x, b)


def generate_summary_graphs(output_dir, gene_frequencies, cluster_sizes, genes_per_isolate, noSamples):
    # write gene frequency histogram
    gene_frequencies = np.array(gene_frequencies)
    plt.hist(gene_frequencies, bins=50, range=(0, 100),
             color='#0504aa', alpha=0.7, rwidth=1.0, edgecolor='black')
    plt.grid(axis='y', alpha=0.75)
    plt.xlim(xmin=0, xmax=100)
    plt.xlabel('Proportion of genomes (%)')
    plt.ylabel('Frequency')
    plt.savefig(output_dir + "gene_frequency.png", format="png")
    plt.clf()

    # write cluster frequency histogram
    cluster_sizes = np.array(cluster_sizes)
    cluster_range = cluster_sizes.max() - cluster_sizes.min()
    cluster_range = 1 if cluster_range <= 0 else cluster_range
    bins = cluster_range if cluster_range < 50 else 50
    plt.hist(cluster_sizes, bins=bins,
             color='#0504aa', alpha=0.7, rwidth=1.0, edgecolor='black')
    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('No. genes per cluster')
    plt.ylabel('Frequency')
    plt.savefig(output_dir + "cluster_size.png", format="png")
    plt.clf()

    # generate rarefaction curve
    rarefaction_list = np.empty(0)
    genome_list = []

    # shuffle genes_per_isolate to generate many samples of genomes
    shuffle_iterations = 50
    for i in range(shuffle_iterations):
        temp_rarefaction_list = []
        prev_set = intbitset([])
        shuffle(genes_per_isolate)
        for mem in range(noSamples):
            genome_list.append(mem + 1)
            new_genes = genes_per_isolate[mem].difference(prev_set)
            temp_rarefaction_list.append(len(new_genes))
            prev_set.update(new_genes)
        temp_rarefaction_list = np.cumsum(temp_rarefaction_list)
        rarefaction_list = np.append(rarefaction_list, [temp_rarefaction_list])

    # model power-law
    genome_list = np.array(genome_list)
    df = pd.DataFrame({'x': genome_list, 'y': rarefaction_list})
    fig, ax = plt.subplots()

    # plot scatter with jitter
    sns.regplot(data=df, x="x", y="y", x_jitter=0.2, fit_reg=False, scatter=True,
                color='silver', ax=ax, scatter_kws={'s': 10, 'alpha': 1, 'edgecolors': 'k'})

    # fit power-law regression, determine confidence intervals
    pars, cov = curve_fit(f=power_law, xdata=genome_list,
                          ydata=rarefaction_list)
    sigma_ab = np.sqrt(np.diagonal(cov))

    b = ufloat(pars[1], sigma_ab[1])

    # determine if pangenome open is upper confidence interval is above 0 for b (aka gamma)
    if (pars[1] + sigma_ab[1]) < 0:
        pangenome_openess = "closed"
    elif (pars[1] - sigma_ab[1]) > 0:
        pangenome_openess = "open"
    else:
        pangenome_openess = "non-significant"
    text_res = r'$\gamma$' + " = {}\nPangenome: {}".format(b, pangenome_openess)

    # generate power-law spline
    genome_list = np.sort(np.unique(genome_list))
    spl = make_interp_spline(genome_list, power_law(genome_list, *pars), k=3)
    xnew = np.linspace(1, genome_list.max(), 100)
    power_smooth = spl(xnew)

    # fit may produce inf values in covariance matrix. If occurs don't plot confidence interval splines.
    try:
        upper_spl = make_interp_spline(genome_list, power_law(genome_list, *(pars + sigma_ab)), k=3)
        lower_spl = make_interp_spline(genome_list, power_law(genome_list, *(pars - sigma_ab)), k=3)
        bound_upper = upper_spl(xnew)
        bound_lower = lower_spl(xnew)
        plt.fill_between(xnew, bound_lower, bound_upper, color='red', alpha=0.15)
    except ValueError:
        print("Fitting of rarefaction curve failed. Power-law fit may look strange.")
        pass

    plt.plot(xnew, power_smooth, linestyle='--', linewidth=2, color='red')

    # set axis variables
    plt.ylim(ymin=0)
    plt.xlim(xmin=0)
    plt.xlabel('Number of genomes sampled')
    plt.ylabel('Cumulative number of genes discovered')
    plt.text(0.7, 0.2, text_res, horizontalalignment='center', verticalalignment='center', transform=ax.transAxes)
    plt.savefig(output_dir + "rarefaction_curve.png", format="png")
    plt.clf()

    return


def generate_nwk_tree(matrix_in, threads, isolate_names, output_dir, alignment):
    if alignment is True:
        # determine distance matrix from gene_presence/absence
        distance_mat = np.array(ggCaller_cpp.get_distances_align(matrix_in, threads)).reshape(len(isolate_names),
                                                                                              len(isolate_names))
        phylip_name = output_dir + 'core_gene_alignment.phylip'
        tree_filename = output_dir + 'core_tree_NJ.nwk'

    else:
        # determine distance matrix from core gene alignment
        distance_mat = np.array(ggCaller_cpp.get_distances_pa(matrix_in, threads)).reshape(len(isolate_names),
                                                                                           len(isolate_names))
        phylip_name = output_dir + 'pangenome_gene_presence_absence.phylip'
        tree_filename = output_dir + 'pangenome_NJ.nwk'

    # generate phylip matrix
    with open(phylip_name, 'w') as pFile:
        pFile.write(str(len(isolate_names)) + "\n")
        for coreDist, iso in zip(distance_mat, isolate_names):
            pFile.write(iso)
            pFile.write(' ' + ' '.join(map(str, coreDist)))
            pFile.write("\n")

    # run rapidnj
    command = ["rapidnj", phylip_name, "-n", "-i", "pd", "-o", "t", "-x",
               tree_filename, "-c", str(threads)]

    result = subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if result.returncode != 0:
        raise Exception("RapidNJ failed to run on file: " + phylip_name)
        sys.exit(1)

    # remove phylip file
    os.remove(phylip_name)

    return


def generate_roary_gene_presence_absence(G, mems_to_isolates, orig_ids,
                                         ids_len_stop, output_dir, threads):
    # hold gene proportions for gene frequency histogram
    gene_frequencies = []

    # hold genes per cluster for cluster frequency histogram
    cluster_sizes = []

    # arange isolates
    isolates = []
    mems_to_index = {}
    for i, mem in enumerate(mems_to_isolates):
        isolates.append(mems_to_isolates[mem])
        mems_to_index[str(mem)] = i

    noSamples = len(isolates)
    # Layout categories
    noCore = 0
    noSoftCore = 0
    noShell = 0
    noCloud = 0
    total_genes = 0

    # hold set of nodes found in each colour for rarefaction curve
    genes_per_isolate = [intbitset([]) for _ in range(noSamples)]

    # hold list of nodes found in each colour for distance calculation
    isolate_gene_list = [[False] * G.number_of_nodes() for _ in range(noSamples)]

    # generate file
    with open(output_dir + "gene_presence_absence_roary.csv", 'w') as roary_csv_outfile, \
            open(output_dir + "gene_presence_absence.csv", 'w') as csv_outfile, \
            open(output_dir + "gene_presence_absence.Rtab", 'w') as Rtab_outfile:
        header = [
                     "Gene", "Non-unique Gene name", "Annotation", "No. isolates",
                     "No. sequences", "Avg sequences per isolate", "Genome Fragment",
                     "Order within Fragment", "Accessory Fragment",
                     "Accessory Order with Fragment", "QC", "Min group size nuc",
                     "Max group size nuc", "Avg group size nuc"
                 ] + isolates
        roary_csv_outfile.write(",".join(header) + "\n")
        csv_outfile.write(",".join(header[:3] + isolates) + "\n")
        Rtab_outfile.write("\t".join((["Gene"] + isolates)) + "\n")

        # Iterate through components writing out to file
        used_gene_names = set([""])
        unique_id_count = 0
        frag = 0
        entry_list = []
        entry_ext_list = []
        pres_abs_list = []
        entry_sizes = []
        entry_count = 0
        for component in nx.connected_components(G):
            frag += 1
            count = 0
            for node in component:
                count += 1
                len_mode = max(G.nodes[node]['lengths'],
                               key=G.nodes[node]['lengths'].count)
                name = '~~~'.join([
                    gn for gn in G.nodes[node]['annotation'].strip().strip(
                        ';').split(';') if gn != ''
                ])
                name = ''.join(e for e in name
                               if e.isalnum() or e in ["_", "~", "|"])
                name = name.replace("|", "-")
                if name not in used_gene_names:
                    entry = [name]
                    used_gene_names.add(name)
                    G.nodes[node]['name'] = name
                else:
                    G.nodes[node]['name'] = "group_" + str(unique_id_count)
                    entry = [G.nodes[node]['name']]
                    unique_id_count += 1
                if G.nodes[node]['annotation'] != '':
                    entry.append(G.nodes[node]['annotation'])
                    entry.append(G.nodes[node]['description'].replace(',', ' '))
                else:
                    entry.append("hypothetical protein")
                    entry.append("NA")
                entry.append(G.nodes[node]['size'])
                entry.append(len(G.nodes[node]['seqIDs']))
                entry.append((1.0 * len(G.nodes[node]['seqIDs'])) /
                             G.nodes[node]['size'])
                entry.append(frag)
                entry.append(count)
                entry += ["", "", ""]
                entry.append(np.min(G.nodes[node]['lengths']))
                entry.append(np.max(G.nodes[node]['lengths']))
                entry.append(np.mean(G.nodes[node]['lengths']))
                pres_abs = [""] * len(isolates)
                pres_abs_ext = [""] * len(isolates)
                entry_size = 0
                for seq in G.nodes[node]['seqIDs']:
                    sample_id = mems_to_index["_".join(seq.split("_")[:-2])]
                    if pres_abs[
                        sample_id] == "":  # ensures we only take the first one
                        if seq in orig_ids:
                            pres_abs[sample_id] = orig_ids[seq]
                            pres_abs_ext[sample_id] = orig_ids[seq]
                        else:
                            pres_abs[sample_id] = seq
                            pres_abs_ext[sample_id] = seq
                        entry_size += 1
                    else:
                        # this is similar to PIRATE output
                        if seq in orig_ids:
                            pres_abs[sample_id] += ";" + orig_ids[seq]
                            pres_abs_ext[sample_id] += ";" + orig_ids[seq]
                        else:
                            pres_abs[sample_id] += ";" + seq
                            pres_abs_ext[sample_id] += ";" + seq
                    if (abs(ids_len_stop[seq][0] - len_mode) /
                        len_mode) > (0.05 * len_mode):
                        pres_abs_ext[sample_id] += "_len"
                    if ids_len_stop[seq][1]:
                        pres_abs_ext[sample_id] += "_stop"

                entry += pres_abs
                entry_list.append(entry)
                entry_ext_list.append(entry[:3] + pres_abs_ext)
                pres_abs_list.append(pres_abs)
                entry_sizes.append((entry_size, entry_count))
                entry_count += 1

                # determine which genomes genes clusters are found in
                for mem in G.nodes[node]['members']:
                    genes_per_isolate[mem].add(node)
                    isolate_gene_list[mem][total_genes] = True

                # determine gene presence/absence
                num_isolates = G.nodes[node]['size']
                proportion_present = float(num_isolates) / noSamples * 100.0
                gene_frequencies.append(proportion_present)
                cluster_sizes.append(entry_size)
                if proportion_present >= 99:
                    noCore += 1
                elif proportion_present >= 95:
                    noSoftCore += 1
                elif proportion_present >= 15:
                    noShell += 1
                else:
                    noCloud += 1
                total_genes += 1

        # sort so that the most common genes are first (as in roary)
        entry_sizes = sorted(entry_sizes, reverse=True)
        for s, i in entry_sizes:
            roary_csv_outfile.write(",".join([str(e)
                                              for e in entry_list[i]]) + "\n")
            csv_outfile.write(",".join([str(e)
                                        for e in entry_ext_list[i]]) + "\n")
            Rtab_outfile.write(entry_list[i][0] + "\t")
            Rtab_outfile.write("\t".join(
                (["0" if e == "" else "1" for e in pres_abs_list[i]])) + "\n")

    # write summary output
    with open(output_dir + "summary_statistics.txt", 'w') as outfile:
        output = ("Core genes\t(99% <= strains <= 100%)\t" + str(noCore) +
                  "\n" + "Soft core genes\t(95% <= strains < 99%)\t" +
                  str(noSoftCore) + "\n" +
                  "Shell genes\t(15% <= strains < 95%)\t" + str(noShell) +
                  "\n" + "Cloud genes\t(0% <= strains < 15%)\t" +
                  str(noCloud) + "\n" +
                  "Total genes\t(0% <= strains <= 100%)\t" + str(total_genes))
        outfile.write(output)

    # generate summary graphs
    generate_summary_graphs(output_dir, gene_frequencies, cluster_sizes, genes_per_isolate, noSamples)

    # generate nwk tree from pangenome if genes detected
    isolate_names = G.graph['isolateNames']
    if G.number_of_nodes() > 0:
        generate_nwk_tree(isolate_gene_list, threads, isolate_names, output_dir, False)

    return G


def generate_pan_genome_reference(G, output_dir, split_paralogs=False):
    # need to treat paralogs differently?
    centroids = set()
    records = []

    for node in G.nodes():
        if not split_paralogs and G.nodes[node]['centroid'][0] in centroids:
            continue
        records.append(
            SeqRecord(Seq(max(G.nodes[node]['dna'], key=lambda x: len(x))),
                      id=G.nodes[node]['name'],
                      description=""))
        for centroid in G.nodes[node]['centroid']:
            centroids.add(centroid)

    with open(output_dir + "pan_genome_reference.fa", 'w') as outfile:
        SeqIO.write(records, outfile, "fasta")

    return


def generate_common_struct_presence_absence(G,
                                            output_dir,
                                            mems_to_isolates,
                                            min_variant_support=2):
    # arange isolates
    isolates = []
    members = []
    for mem in mems_to_isolates:
        isolates.append(mems_to_isolates[mem])
        members.append(mem)

    struct_variants = {}
    for node in G.nodes():
        if G.degree[node] < 3: continue  # skip as linear
        for path in iter.combinations(G.edges(node), 2):
            in_both = (G[path[0][0]][path[0][1]]['members']
                       & G[path[1][0]][path[1][1]]['members'])
            if len(in_both) >= min_variant_support:
                struct_variants[(path[0][0], path[0][1], path[1][1])] = in_both

    header = []
    for variant in struct_variants:
        header.append("-".join([
            G.nodes[variant[1]]['name'], G.nodes[variant[0]]['name'],
            G.nodes[variant[2]]['name']
        ]))

    with open(output_dir + "struct_presence_absence.Rtab",
              'w') as Rtab_outfile:
        Rtab_outfile.write("\t".join((["Gene"] + isolates)) + "\n")
        for h, variant in zip(header, struct_variants):
            variant_calls = [h]
            for member in members:
                if member in struct_variants[variant]:
                    variant_calls.append("1")
                else:
                    variant_calls.append("0")
            Rtab_outfile.write("\t".join(variant_calls) + "\n")

    return


def generate_pan_genome_alignment(G, temp_dir, output_dir, threads,
                                  isolate_names, shd_arr_tup, high_scoring_ORFs, overlap, pool, ref_aln,
                                  call_variants, verbose, ignore_pseduogenes, truncation_threshold):
    unaligned_sequence_files = []
    unaligned_reference_files = []
    # Make a folder for the output alignments, clear if present and remake
    try:
        os.mkdir(output_dir + "aligned_gene_sequences")
    except FileExistsError:
        shutil.rmtree(output_dir + "aligned_gene_sequences")
        os.mkdir(output_dir + "aligned_gene_sequences")

    # Multithread writing gene sequences to disk (temp directory) so aligners can find them
    for outname, ref_outname in pool.map(partial(output_alignment_sequence,
                                                 temp_directory=temp_dir, outdir=output_dir, shd_arr_tup=shd_arr_tup,
                                                 high_scoring_ORFs=high_scoring_ORFs, overlap=overlap, ref_aln=ref_aln,
                                                 ignore_pseduogenes=ignore_pseduogenes,
                                                 truncation_threshold=truncation_threshold),
                                         G.nodes(data=True)):
        unaligned_sequence_files.append(outname)
        unaligned_reference_files.append(ref_outname)
    if ref_aln:
        # centroid files with paired sequence files
        ref_seq_pairs = [
            (temp_dir + unaligned_reference_files[i].split("/")[-1].split("_ref.")[0] + "_ref.aln.fas",
             unaligned_sequence_files[i])
            for i in range(0, len(unaligned_reference_files)) if unaligned_sequence_files[i] is not None
                                                                 and unaligned_reference_files[i] is not None]
        # reference files with no paired sequence files
        ref_seq_singles = [
            temp_dir + unaligned_reference_files[i].split("/")[-1].split("_ref.")[0] + "_ref.aln.fas"
            for i in range(0, len(unaligned_reference_files))
            if unaligned_sequence_files[i] is None and unaligned_reference_files[i] is not None]

    # remove single sequence files
    unaligned_sequence_files = filter(None, unaligned_sequence_files)
    unaligned_reference_files = filter(None, unaligned_reference_files)

    # create set of files not to be used for variant calling
    if call_variants:
        no_vc_set = set(os.listdir(output_dir + "aligned_gene_sequences"))

    # Get Biopython command calls for each output gene sequences
    # check if ref_alignment being done
    if ref_aln:
        # conduct MSA on reference files, first using standard MSA
        commands = [
            get_alignment_commands(fastafile, None, "def")
            for fastafile in unaligned_reference_files if "_ref.aln.fas" not in fastafile
        ]
        if verbose: print("Aligning centroids...")
        multi_align_sequences(commands, temp_dir, threads, "def", not verbose)

        # move any centroid alignments which do not have associated sequence files
        for file in ref_seq_singles:
            os.rename(file, output_dir + "aligned_gene_sequences/" + file.split("/")[-1].split("_ref.")[0] + '.aln.fas')

        # repeat with reference-guided alignment
        commands = [
            get_alignment_commands(fastapair, output_dir, "ref")
            for fastapair in ref_seq_pairs
        ]
        if verbose: print("Aligning remaining sequences...")
        multi_align_sequences(commands, output_dir + "aligned_gene_sequences/",
                              threads, "ref", not verbose)
    else:
        commands = [
            get_alignment_commands(fastafile, output_dir, "def")
            for fastafile in unaligned_sequence_files
        ]
        # Run these commands in a multi-threaded way
        multi_align_sequences(commands, output_dir + "aligned_gene_sequences/",
                              threads, "def", not verbose)

    # back translate sequences
    back_translate_dir(high_scoring_ORFs, isolate_names, output_dir + "aligned_gene_sequences/", overlap, shd_arr_tup,
                       pool)

    # call variants using snp-sites
    if call_variants:
        try:
            os.mkdir(output_dir + "VCF")
        except FileExistsError:
            None
        run_snpsites_dir(output_dir + "aligned_gene_sequences", output_dir + "VCF", no_vc_set, pool)

    return


def get_unannotated_nodes(G):
    # Get the unannotated nodes based on bitscore
    unannotated_nodes = []
    for node in G.nodes(data=True):
        if float(G.nodes[node[0]]["bitscore"]) == 0:
            unannotated_nodes.append(node)
    return unannotated_nodes


def get_core_gene_nodes(G, threshold, num_isolates):
    # Get the core genes based on percent threshold
    core_nodes = []
    for node in G.nodes(data=True):
        if float(G.nodes[node[0]]["size"]) / float(num_isolates) > threshold:
            core_nodes.append(node)
    return core_nodes


def check_rapidnj_install():
    command = ["rapidnj", "-h"]

    p = str(
        subprocess.run(command,
                       stdout=subprocess.PIPE,
                       stderr=subprocess.PIPE))

    present = False

    find_ver = re.search(r'Rapid neighbour-joining', p)
    if find_ver != None:
        present = True

    if present == False:
        sys.stderr.write("Need rapidnj to be installed and available in PATH. " + "\n")
        sys.exit(1)

    return present


def concatenate_core_genome_alignments(core_names, output_dir, isolate_names, threads):
    alignments_dir = output_dir + "/aligned_gene_sequences/"
    # Open up each alignment that is associated with a core node
    alignment_filenames = os.listdir(alignments_dir)
    core_filenames = [
        x for x in alignment_filenames if x.split('.')[0] in core_names
    ]
    # Read in all these alignments
    gene_alignments = []
    isolates = set()
    for filename in core_filenames:
        gene_name = os.path.splitext(os.path.basename(filename))[0]
        alignment = AlignIO.read(alignments_dir + filename, 'fasta')
        gene_dict = {}
        for record in alignment:
            genome_id = record.id.rsplit("_", 1)[0]
            if genome_id in gene_dict:
                if str(record.seq).count("-") < str(
                        gene_dict[genome_id][1]).count("-"):
                    gene_dict[genome_id] = (record.id, record.seq)
            else:
                gene_dict[genome_id] = (record.id, record.seq)
            gene_length = len(record.seq)
            isolates.add(genome_id)
        gene_alignments.append((gene_name, gene_dict, gene_length))
    # Combine them
    isolate_aln = []
    for iso in isolates:
        seq = ""
        for gene in gene_alignments:
            if iso in gene[1]:
                seq += gene[1][iso][1]
            else:
                seq += "-" * gene[2]
        isolate_aln.append(SeqRecord(seq, id=iso, description=""))

    # generate nwk tree from core genome alignment if alignments present
    if isolate_aln:
        generate_nwk_tree([str(record.seq) for record in isolate_aln], threads, isolate_names, output_dir, True)

    # Write out the two output files
    SeqIO.write(isolate_aln, output_dir + 'core_gene_alignment.aln', 'fasta')

    write_alignment_header(gene_alignments, output_dir)
    return core_filenames


def generate_core_genome_alignment(G, temp_dir, output_dir, threads,
                                   isolate_names, threshold, num_isolates, shd_arr_tup, high_scoring_ORFs,
                                   overlap, pool, ref_aln, call_variants, verbose,
                                   ignore_pseduogenes, truncation_threshold):
    unaligned_sequence_files = []
    unaligned_reference_files = []
    # Make a folder for the output alignments, clear if present and remake
    try:
        os.mkdir(output_dir + "aligned_gene_sequences")
    except FileExistsError:
        shutil.rmtree(output_dir + "aligned_gene_sequences")
        os.mkdir(output_dir + "aligned_gene_sequences")

    # Get core nodes
    core_genes = get_core_gene_nodes(G, threshold, num_isolates)
    core_gene_names = [G.nodes[x[0]]["name"] for x in core_genes]

    # Multithread writing gene sequences to disk (temp directory) so aligners can find them
    for outname, ref_outname in pool.map(partial(output_alignment_sequence,
                                                 temp_directory=temp_dir, outdir=output_dir, shd_arr_tup=shd_arr_tup,
                                                 high_scoring_ORFs=high_scoring_ORFs, overlap=overlap, ref_aln=ref_aln,
                                                 ignore_pseduogenes=ignore_pseduogenes,
                                                 truncation_threshold=truncation_threshold),
                                         core_genes):
        unaligned_sequence_files.append(outname)
        unaligned_reference_files.append(ref_outname)

    if ref_aln:
        # centroid files with paired sequence files
        ref_seq_pairs = [
            (temp_dir + unaligned_reference_files[i].split("/")[-1].split("_ref.")[0] + "_ref.aln.fas",
             unaligned_sequence_files[i])
            for i in range(0, len(unaligned_reference_files)) if unaligned_sequence_files[i] is not None
                                                                 and unaligned_reference_files[i] is not None]
        # reference files with no paired sequence files
        ref_seq_singles = [
            temp_dir + unaligned_reference_files[i].split("/")[-1].split("_ref.")[0] + "_ref.aln.fas"
            for i in range(0, len(unaligned_reference_files))
            if unaligned_sequence_files[i] is None and unaligned_reference_files[i] is not None]

    # remove single sequence files
    unaligned_sequence_files = filter(None, unaligned_sequence_files)
    unaligned_reference_files = filter(None, unaligned_reference_files)

    # create set of files not to be used for variant calling
    if call_variants:
        no_vc_set = set(os.listdir(output_dir + "aligned_gene_sequences/"))

    # Get Biopython command calls for each output gene sequences
    # check if ref_alignment being done
    if ref_aln:
        # conduct MSA on reference files, first using standard MSA
        commands = [
            get_alignment_commands(fastafile, None, "def")
            for fastafile in unaligned_reference_files if "_ref.aln.fas" not in fastafile
        ]
        if verbose: print("Aligning centroids...")
        multi_align_sequences(commands, temp_dir, threads, "def", not verbose)

        # move any centroid alignments which do not have associated sequence files
        for file in ref_seq_singles:
            os.rename(file, output_dir + "aligned_gene_sequences/" + file.split("/")[-1].split("_ref.")[0] + '.aln.fas')

        # repeat with reference-guided alignment
        commands = [
            get_alignment_commands(fastapair, output_dir, "ref")
            for fastapair in ref_seq_pairs
        ]
        if verbose: print("Aligning remaining sequences...")
        multi_align_sequences(commands, output_dir + "aligned_gene_sequences/",
                              threads, "ref", not verbose)
    else:
        # Get alignment commands
        commands = [
            get_alignment_commands(fastafile, output_dir, "def")
            for fastafile in unaligned_sequence_files
        ]
        # Run alignment commands
        multi_align_sequences(commands, output_dir + "aligned_gene_sequences/",
                              threads, "def", not verbose)
    # back translate sequences
    back_translate_dir(high_scoring_ORFs, isolate_names, output_dir + "aligned_gene_sequences/", overlap, shd_arr_tup,
                       pool)

    # call variants using snp-sites
    if call_variants:
        try:
            os.mkdir(output_dir + "VCF")
        except FileExistsError:
            None
        run_snpsites_dir(output_dir + "aligned_gene_sequences", output_dir + "VCF", no_vc_set, pool)

    # Concatenate them together to produce the two output files
    concatenate_core_genome_alignments(core_gene_names, output_dir, isolate_names, threads)

    return


def generate_summary_stats(output_dir):
    with open(output_dir + "gene_presence_absence_roary.csv", 'r') as inhandle:
        gene_presence_absence = inhandle.read().splitlines()[1:]
    noSamples = len(gene_presence_absence[0].split(',')) - 14
    # Layout categories
    noCore = 0
    noSoftCore = 0
    noShell = 0
    noCloud = 0
    total_genes = 0
    # Iterate through GPA and summarise
    for gene in gene_presence_absence:
        proportion_present = float(gene.split(',')[3]) / noSamples * 100.0
        if proportion_present >= 99:
            noCore += 1
        elif proportion_present >= 95:
            noSoftCore += 1
        elif proportion_present >= 15:
            noShell += 1
        else:
            noCloud += 1
        total_genes += 1

    # write output
    with open(output_dir + "summary_statistics.txt", 'w') as outfile:
        output = ("Core genes\t(99% <= strains <= 100%)\t" + str(noCore) +
                  "\n" + "Soft core genes\t(95% <= strains < 99%)\t" +
                  str(noSoftCore) + "\n" +
                  "Shell genes\t(15% <= strains < 95%)\t" + str(noShell) +
                  "\n" + "Cloud genes\t(0% <= strains < 15%)\t" +
                  str(noCloud) + "\n" +
                  "Total genes\t(0% <= strains <= 100%)\t" + str(total_genes))
        outfile.write(output)

    return True
