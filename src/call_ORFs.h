#ifndef BIFROST_API_GGCALLER_H
#define BIFROST_API_GGCALLER_H

#include "unitigDict.h"
#include "match_string.h"
#include "indexing.h"

// call_ORFs
ORFNodeMap generate_ORFs(const GraphVector& graph_vector,
                         const std::vector<std::string>& stop_codons,
                         const std::vector<std::string>& start_codons,
                         const std::vector<int>& unitig_path,
                         const int& overlap,
                         const size_t min_len,
                         const bool is_ref,
                         const fm_index_coll& fm_idx);

std::tuple<std::string, std::vector<int>, std::vector<indexPair>> calculate_coords(const std::pair<std::size_t, std::size_t>& codon_pair,
                                                                                    const std::vector<int>& nodelist,
                                                                                    const std::vector<std::vector<size_t>>& node_ranges,
                                                                                    const int& overlap);

std::pair<ORFVector, NodeStrandMap> call_ORFs(const AllPaths& all_paths,
                                             const GraphVector& graph_vector,
                                             const std::vector<std::string>& stop_codons_for,
                                             const std::vector<std::string>& start_codons_for,
                                             const int overlap,
                                             const size_t min_ORF_length,
                                             const bool is_ref,
                                             const fm_index_coll& fm_idx);

ORFVector sort_ORF_indexes(ORFNodeMap& ORF_node_map);

NodeStrandMap calculate_pos_strand(const ORFNodeMap& ORF_node_map);

#endif //BIFROST_API_GGCALLER_H
