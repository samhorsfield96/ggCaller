Tutorial
==================================

Here we'll walk through a typical run of ggCaller, including both :ref:`Gene-calling` and :ref:`Querying`.

Example results can be found `here <https://figshare.com/articles/dataset/ggCaller_example_results/21831204>`_.

.. important::
    Results will be consistent, but may not exactly match between your run
    and the example. This is due to the greedy clustering
    algorithm used by ggCaller, which can cause small differences
    in genes counts.

Installation and setup
----------------------

Follow the guide in :doc:`installation` for downloading and installing ggCaller.

Working Dataset
---------------

We'll use a dataset from `Bentley et al. (2006) <https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.0020031>`_.
This dataset contains 91 sequences pneumococcal capsular polysaccharide synthetic (CPS) loci. These sequences are structurally diverse,
but are only ~20,000 bp in length, so can be analysed quickly (~5-10 minutes) on a standard laptop or desktop.

Download the files from `here <https://figshare.com/articles/dataset/Bentley_et_al_2006_CPS_sequences/21829038>`_ and unzip::

    tar xvf Bentley_et_al_2006_CPS_sequences.tar.bz2

We will also provide our own custom annotation database for DIAMOND. These will be the manually curated protein sequences
from Bentley et al. Download from `here <https://figshare.com/articles/dataset/Bentley_et_al_2006_CPS_protein_sequences/21829071>`_ and unzip::

    tar xvf Bentley_et_al_2006_CPS_protein_sequences.tar.bz2

Gene-calling
------------

First generate an input file for ggCaller. This must be a file containing paths (absolute recommended) to all sequences to be analysed.
We recommend running the below command within the unzipped to generate this file::

    cd Bentley_et_al_2006_CPS_sequences
    ls -d -1 $PWD/*.fa > input.txt
    cd ..

``input.txt`` will now contain absolute paths to all ``.fa`` files in the directory ``Bentley_et_al_2006_CPS_sequences``.

Now we will run ggCaller specifying the below settings:

- Sensitive DIAMOND annotation using a custom database, and HMMER3 using the default database
- Pangenome-wide alignment using default MAFFT
- Saved intermediate datastructures, enabling sequence querying

To do this using 4 threads, run::

    ggcaller --refs Bentley_et_al_2006_CPS_sequences/input.txt --annotation ultrasensitive --diamonddb Bentley_et_al_2006_CPS_protein_sequences.faa --aligner def --alignment pan --save --out ggc_Bentley_et_al_CPS --threads 4

You will find the following files in the output directory ``ggc_Bentley_et_al_CPS``:

- ``cluster_size.png``: a frequency distribution of clusters by the number of genes found within them

.. image:: images/cluster_size_example.png
   :alt:  A frequency distribution of clusters by size
   :align: center

- ``gene_frequency.png``: a frequency distribution of clusters by proportion of dataset

.. image:: images/gene_frequency_example.png
   :alt:  A frequency distribution of clusters by proportion of dataset
   :align: center

- ``rarefaction_curve.png``: rarefaction curve, describes the number of new genes discovered with random addition of a single genome. Also includes power-law fit for determination of pangenome openness, based on `Tettelin et al. (2005) <https://www.pnas.org/doi/10.1073/pnas.0506758102>`_.

.. image:: images/rarefaction_curve_example.png
   :alt:  A gene rarefaction curve
   :align: center

- ``core_gene_alignment.aln``: concatenated core genome alignment
- ``core_alignment_header.embl``: core genome alignment in EMBL format
- ``core_tree_NJ.nwk``: Neighbour joining tree from core genome alignment generated by `RapidNJ <https://birc.au.dk/software/rapidnj/>`_. This can be visualised in `Microreact <https://microreact.org/>`_

.. image:: images/core_genome_tree.png
   :alt:  A core genome neighbour joining tree
   :align: center

- ``pangenome_NJ.nwk``: Neighbour joining tree from gene presence/absence matrix generated by RapidNJ (can also be visualised in Microreact).

.. image:: images/pangenome_tree.png
   :alt:  A gene presence/absence neighbour joining tree
   :align: center

- ``pan_genome_reference.fa``: contains centroids for each cluster in FASTA format
- ``gene_calls.faa`` and ``gene_calls.ffn``: gene predictions with annotations in amino-acid and nucleotide FASTA formats
- ``pre_filt_graph.gml`` and ``final_graph.gml``: gene graphs pre- and post-quality control with Panaroo
- ``gene_presence_absence*``: gene presence absence files in three formats; Roary-CSV, CSV and Rtab
- ``struct_presence_absence.Rtab``: structural variant presnce/absence matrix
- ``summary_statistics.txt``: summary of gene frequencies based on Roary
- ``VCF``: directory containing VCF files for each cluster generated by `SNP-SITES <https://github.com/sanger-pathogens/snp-sites>`_
- ``aligned_gene_sequences``: directory of alignment files for each cluster in FASTA format
- ``GFF``: directory of GFF files for each sample in GFF3 format
- ``ggc_data``: intermediate datastructures written to disk, required for querying.
- ``ORF_dir``: intermediate datastructures written to disk, containing gene predictions.
- ``Path_dir``: intermediate datastructures written to disk, containing genome paths through the DBG.

Querying the graph
------------------

We can now query the graph. To do so, run::

    ggcaller --query CPS_queries.fasta --graph Bentley_et_al_2006_CPS_sequences/input.gfa --colours Bentley_et_al_2006_CPS_sequences/input.color.bfg --prev_run ggc_Bentley_et_al_CPS --out ggc_Bentley_et_al_CPS --threads 4

Results will be saved in ``ggc_Bentley_et_al_CPS/matched_queries.fasta``.

Details on the output can be found in :ref:`Interpreting results`.

From ``matched_queries.fasta``, we can see that all the genes queried were identified in the graph.

As we searched for specific gene variants, this search was too stringent to return orthologues in other genomes.

.. important::
    We recommend searching for partial gene sequences,
    or lowering ``--query-id`` to return more distantly related sequences.