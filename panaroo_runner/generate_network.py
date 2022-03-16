from collections import defaultdict
import networkx as nx
from intbitset import intbitset


def generate_network(high_scoring_ORFs, high_scoring_ORF_edges,
                     cluster_id_list, cluster_dict):
    # associate sequences with their clusters
    seq_to_cluster = {}
    seqid_to_centroid = {}
    cluster_centroids = {}
    cluster_members = defaultdict(list)
    cluster_centroid_data = {}
    # iterate over cluster_dict, parsing all sequences to clusters
    cluster_id = 0
    for centroid, ORF_list in cluster_dict.items():
        # generate a panaroo sequence ID for current centroid
        centroid_genome_id = cluster_id_list[centroid][0]
        centroid_local_id = cluster_id_list[centroid][1]

        # access ORF information for centroid from high_scoring_ORFs, ensuring cluster is present
        if centroid_local_id not in high_scoring_ORFs[centroid_genome_id]:
            # if centroid not present, go through ORF_list and find next largest ORF present
            current_length = 0
            for ORF_id in ORF_list:
                genome_id = cluster_id_list[ORF_id][0]
                local_id = cluster_id_list[ORF_id][1]
                if local_id in high_scoring_ORFs[genome_id]:
                    ORFNodeVector = high_scoring_ORFs[genome_id][local_id]
                    # if longer, assign as centroid
                    if ORFNodeVector[2] > current_length:
                        current_length = ORFNodeVector[2]
                        centroid_genome_id = genome_id
                        centroid_local_id = local_id
            # if no new centroids assigned, pass
            if current_length == 0:
                continue

        # access ORF information for centroid from high_scoring_ORFs
        ORFNodeVector = high_scoring_ORFs[centroid_genome_id][centroid_local_id]
        pan_centroid_ID = str(centroid_genome_id) + "_0_" + str(centroid_local_id)

        # add information to cluster_centroid_data
        cluster_centroid_data[cluster_id] = {
            'ORF_info': ORFNodeVector,
            'annotation': '',
            'description': '',
        }

        # append centroid to cluster
        cluster_centroids[cluster_id] = pan_centroid_ID

        # for each ORF, add ORF_id and cluster_id to respective dictionaries
        for ORF_id in ORF_list:
            # generate a panaroo sequence ID for current ORF
            genome_id = cluster_id_list[ORF_id][0]
            local_id = cluster_id_list[ORF_id][1]

            # check if ORF is present in high_scoring_ORFs
            if local_id not in high_scoring_ORFs[genome_id]:
                continue

            pan_ORF_id = str(genome_id) + "_0_" + str(local_id)

            # index sequences to clusters
            seq_to_cluster[pan_ORF_id] = cluster_id
            cluster_members[cluster_id].append(pan_ORF_id)

        cluster_id += 1

    # determine paralogs if required
    paralogs = set()
    for clust in cluster_members:
        # determine if paralogs present by checking colours for each entry in cluster_id_list
        genomes = [int(s.split("_")[0]) for s in cluster_members[clust]]
        if len(genomes) != len(set(genomes)):
            paralogs.add(clust)

    # build graph using adjacency information and optionally split paralogs
    G = nx.Graph()
    centroid_context = defaultdict(list)
    n_nodes = len(cluster_members)

    # iterating over each cluster, adding edges between clusters that each have a connection
    for current_cluster, ORF_members in cluster_members.items():
        # check if cluster is currently in graph
        add_cluster = True if not G.has_node(current_cluster) else False

        # check if current cluster contains paralogs
        has_paralogs = True if current_cluster in paralogs else False

        for ORF_id in ORF_members:
            # parse genome id and local ORF id
            parsed_id = ORF_id.split("_")
            genome_id = int(parsed_id[0])
            local_id = int(parsed_id[-1])

            # map ORF ID to centroid ID
            seqid_to_centroid[ORF_id] = cluster_centroids[current_cluster]

            # parse neighbour information for current ORF if present, otherwise make empty edge_set
            if local_id in high_scoring_ORF_edges[genome_id]:
                edge_set = high_scoring_ORF_edges[genome_id][local_id]
                has_end = False
            else:
                edge_set = set()
                # check if current ORF is end of contig
                has_end = True

            # initialise cluster to add
            cluster_to_add = current_cluster

            if add_cluster:
                if has_paralogs:
                    # create a new paralog
                    n_nodes += 1
                    cluster_to_add = n_nodes
                    centroid_context[
                        cluster_centroids[current_cluster]].append(
                        [cluster_to_add, genome_id])
                G.add_node(
                    cluster_to_add,
                    size=1,
                    centroid=[cluster_centroids[current_cluster]],
                    maxLenId=0,
                    members=intbitset([genome_id]),
                    seqIDs=set([ORF_id]),
                    hasEnd=has_end,
                    ORF_info=[(cluster_centroid_data[current_cluster]['ORF_info'][0],
                               cluster_centroid_data[current_cluster]['ORF_info'][1])],
                    annotation=cluster_centroid_data[current_cluster]
                    ['annotation'],
                    bitscore=0,
                    description=cluster_centroid_data[current_cluster]
                    ['description'],
                    lengths=[
                        cluster_centroid_data[current_cluster]['ORF_info'][2]
                    ],
                    paralog=has_paralogs,
                    mergedDNA=False)
                # if has_paralogs == true, then need to add new node for each sequence
                # otherwise can stop adding clusters
                if not has_paralogs:
                    add_cluster = False

            else:
                # check if ORF_id already added to the cluster
                if ORF_id not in G.nodes[current_cluster]['seqIDs']:
                    G.nodes[current_cluster]['size'] += 1
                    G.nodes[current_cluster]['members'].add(genome_id)
                    G.nodes[current_cluster]['seqIDs'].add(ORF_id)
                    if G.nodes[current_cluster]['hasEnd'] == False: G.nodes[current_cluster]['hasEnd'] = has_end

            # iterate through edge_set, adding nodes and then adding edges if required
            for neighbour in edge_set:
                # generate panaroo neighbour id
                pan_neigbour_id = str(genome_id) + "_0_" + str(neighbour)

                neighbour_cluster = seq_to_cluster[pan_neigbour_id]

                # check that neighbour exists, if not add to graph
                add_neighbour = True if not G.has_node(neighbour_cluster) else False

                # check if current cluster contains paralogs
                neighbour_has_paralogs = True if neighbour_cluster in paralogs else False

                # initialise cluster to add
                neighbour_cluster_to_add = neighbour_cluster

                # parse neighbour information for current ORF
                if neighbour in high_scoring_ORF_edges[genome_id]:
                    neighbour_has_end = False
                else:
                    neighbour_has_end = True

                # add neighbour if not present in graph already
                if add_neighbour:
                    if neighbour_has_paralogs:
                        # create a new paralog
                        n_nodes += 1
                        neighbour_cluster_to_add = n_nodes
                        centroid_context[
                            cluster_centroids[neighbour_cluster]].append(
                            [neighbour_cluster_to_add, genome_id])
                    G.add_node(
                        neighbour_cluster_to_add,
                        size=1,
                        centroid=[cluster_centroids[neighbour_cluster]],
                        maxLenId=0,
                        members=intbitset([genome_id]),
                        seqIDs=set([pan_neigbour_id]),
                        hasEnd=neighbour_has_end,
                        ORF_info=[(cluster_centroid_data[neighbour_cluster]['ORF_info'][0],
                                   cluster_centroid_data[neighbour_cluster]['ORF_info'][1])],
                        annotation=cluster_centroid_data[neighbour_cluster]
                        ['annotation'],
                        bitscore=0,
                        description=cluster_centroid_data[neighbour_cluster]
                        ['description'],
                        lengths=[
                            cluster_centroid_data[neighbour_cluster]['ORF_info'][2]
                        ],
                        paralog=neighbour_has_paralogs,
                        mergedDNA=False)
                else:
                    # check if pan_neigbour_id already added to the cluster
                    if pan_neigbour_id not in G.nodes[neighbour_cluster]['seqIDs']:
                        G.nodes[neighbour_cluster]['size'] += 1
                        G.nodes[neighbour_cluster]['members'].add(genome_id)
                        G.nodes[neighbour_cluster]['seqIDs'].add(pan_neigbour_id)
                        if G.nodes[neighbour_cluster]['hasEnd'] == False:
                            G.nodes[neighbour_cluster]['hasEnd'] = neighbour_has_end

                # add edge between current ORF and neighbour
                if G.has_edge(cluster_to_add, neighbour_cluster_to_add):
                    G[cluster_to_add][neighbour_cluster_to_add]['size'] += 1
                    G[cluster_to_add][neighbour_cluster_to_add]['members'].add(genome_id)
                else:
                    G.add_edge(cluster_to_add,
                               neighbour_cluster_to_add,
                               size=1,
                               members=intbitset([genome_id]))

                # currently don't add edge between current cluster paralog and neighbour paralog.

    return G, centroid_context, seqid_to_centroid