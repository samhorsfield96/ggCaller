//
// Created by sth19 on 01/09/2021.
//

#ifndef GGCALLER_ORF_CLUSTERING_H
#define GGCALLER_ORF_CLUSTERING_H

#include "unitigDict.h"
#include "gene_overlap.h"
#include "translation.h"

void assign_centroids(const ColoredCDBG<MyUnitigMap>& ccdbg,
                      const std::vector<Kmer>& head_kmer_arr,
                      const size_t& overlap, 
                      const ORFNodeVector& ORF_info,
                      std::vector<std::tuple<int, int, size_t, size_t, std::shared_ptr<std::string>>>& centroid_vector);

ORFGroupPair group_ORFs(const std::map<size_t, std::string>& ORF_file_paths,
                        const ColoredCDBG<MyUnitigMap>& ccdbg,
                        const std::vector<Kmer>& head_kmer_arr,
                        const size_t& overlap,
                        tbb::concurrent_unordered_map<std::string, ORFNodeVector> centroid_map);

std::pair<ORFClusterMap, robin_hood::unordered_map<std::string, std::string>> produce_clusters(const std::map<size_t, std::string>& ORF_file_paths,
                               const ColoredCDBG<MyUnitigMap>& ccdbg,
                               const std::vector<Kmer>& head_kmer_arr,
                               const size_t& overlap,
                               ORFGroupPair& ORF_group_pair,
                               const double& id_cutoff,
                               const double& len_diff_cutoff);

double align_seqs(const std::string& ORF1_aa,
                  const std::string& ORF2_aa);

#endif //GGCALLER_ORF_CLUSTERING_H
