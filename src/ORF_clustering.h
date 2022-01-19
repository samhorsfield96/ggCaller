//
// Created by sth19 on 01/09/2021.
//

#ifndef GGCALLER_ORF_CLUSTERING_H
#define GGCALLER_ORF_CLUSTERING_H

#include "unitigDict.h"
#include "gene_overlap.h"
#include "translation.h"

ORFGroupTuple group_ORFs(const ColourORFMap& colour_ORF_map,
                         const std::vector<Kmer>& head_kmer_arr);

ORFClusterMap produce_clusters(const ColourORFMap& colour_ORF_map,
                               const ColoredCDBG<MyUnitigMap>& ccdbg,
                               const std::vector<Kmer>& head_kmer_arr,
                               const size_t& DBG_overlap,
                               const ORFMatrixVector& ORF_mat_vector,
                               const std::vector<std::unordered_set<size_t>>& ORF_group_vector,
                               std::vector<std::vector<std::pair<size_t, size_t>>>& centroid_vector,
                               const double& id_cutoff,
                               const double& len_diff_cutoff);

//std::string generate_sequence_private(const std::vector<int>& nodelist,
//                                      const std::vector<indexPair>& node_coords,
//                                      const size_t& overlap,
//                                      const GraphVector& _GraphVector);

double align_seqs(const ORFNodeVector& ORF1_info,
                  const ORFNodeVector& ORF2_info,
                  const ColoredCDBG<MyUnitigMap>& ccdbg,
                  const std::vector<Kmer>& head_kmer_arr,
                  const size_t& DBG_overlap);

#endif //GGCALLER_ORF_CLUSTERING_H
