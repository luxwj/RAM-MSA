#ifndef ANYTIMEASTAR_TPP
#define ANYTIMEASTAR_TPP

#include "AnytimeAstar.hpp"

// return false if all cur_crd[d] < seq_lens[d]
template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
bool AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::is_crd_out_of_bound(const NodeCoord &crd, const int dim) {
    bool res = true;
    for (int d = 0; d < dim; ++d)
        if (crd[d] < seq_lens[d]) {
            res = false;
            break;
        }
    return res;
}


// pre-load the match / mismatch scores of each sequence pair
template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
void AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::get_match_mismatch_scores(STYPE *mm_scores, const NodeCoord &cur_crd, const int dim, const bool reverse_seq) {
    NodeCoord all_dir_nbr_crd (cur_crd); // "+= 1" on each dim
    for (int idx = 0; idx < dim; ++idx) all_dir_nbr_crd[idx] += 1;

    int symbols[dim];        // symbols of all_dir_nbr_crd
    for (int i = 0; i < dim; ++i) symbols[i] = -1;

    for (int seq_idx = 0; seq_idx < dim; ++seq_idx) {
        if (all_dir_nbr_crd[seq_idx] > seq_lens[seq_idx])
            continue;

        if (reverse_seq) {  // read the symbol from the end of the sequence
            // crd = 0, symbol = seq_lens[seq_idx]
            // crd = 1, symbol = seq_lens[seq_idx] - 1
            symbols[seq_idx] = sequences[seq_idx][seq_lens[seq_idx] - all_dir_nbr_crd[seq_idx]];     // -1 is necessary
        } else {            // read the symbol normally
            // crd = 0, symbol = -1
            // crd = 1, symbol = 0
            symbols[seq_idx] = sequences[seq_idx][all_dir_nbr_crd[seq_idx] - 1];     // -1 is necessary
        }
    }

    // load the match / mismatch scores
    for (int seq_i = 0; seq_i < dim - 1; ++seq_i) 
        for (int seq_j = seq_i + 1; seq_j < dim; ++seq_j) 
            if (symbols[seq_i] != -1 && symbols[seq_j] != -1)
                mm_scores[seq_i * NN + seq_j] = blosum_table.get_score_char(symbols[seq_i], symbols[seq_j]);

}


template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
bool AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::closed_list_duplication_detection(bool &use_reexp_fscore_instead, const NodeCoord &nbr_crd, const GapLen &nbr_tb_info, STYPE cur_gscore, STYPE &cur_fscore) {
    auto closed_iter = closed_list_gap_len.find(nbr_crd);
    if (closed_iter != closed_list_gap_len.end()) {
        auto inner_iter = closed_iter->second.find(nbr_tb_info);
        if (inner_iter != closed_iter->second.end()) {
            // inner_iter->second = std::pair{gscore, reexp_fscore}
            if (cost_instead_of_score && (cur_gscore > inner_iter->second.first) || 
                !cost_instead_of_score && (cur_gscore < inner_iter->second.first)) {
                workload_recorder.node_pruned_cnt_gscore += 1;
                return true;   
            } else if (cur_gscore == inner_iter->second.first) {
                // memory bound
                // if reexp-fscore is not the initial value, use its reexp-fscore instead if it's found in the Closed List
                if (inner_iter->second.second != init_reexp_fscore) {
                    cur_fscore = inner_iter->second.second;
                    use_reexp_fscore_instead = true;
                } else {
                    // not a reexp node
                    workload_recorder.node_pruned_cnt_gscore += 1;
                    return true;
                }
            } else {
                /* --- Shouldn't erase any node because we might need them in backtracking suboptimal alignments --- */
                
                // // erase the node in closed list
                // closed_iter->second.erase(inner_iter);
                // if (closed_iter->second.empty())
                //     closed_list_gap_len.erase(closed_iter);
            }
        }
    }
    return false;
}


template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
void AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::open_list_insertion(STYPE cur_gscore, STYPE cur_fscore, std::vector<OLCoord*>* target_index_by_coord,
    const NodeCoord &nbr_crd, const GapLen &nbr_tb_info) {

    typename AstarMultiIndex<NN>::CoordTBInfo nbr_crd_tb_info = {nbr_crd, nbr_tb_info};

    if (sum_of_crd_as_level) {
        int nbr_bin = 0;
        for (int i = 0; i < NN; ++i) nbr_bin += nbr_crd[i];
        if (nbr_bin >= bin_cnt) nbr_bin = bin_cnt - 1;      // mainly for the recursive MSA because bin_cnt is set to 1 in lower dimensions 
        // nbr_crd cannot be in other bins because the sum of coordiantes is the same
        auto iterFound = (*target_index_by_coord)[nbr_bin]->find(nbr_crd_tb_info);
        if (iterFound != (*target_index_by_coord)[nbr_bin]->end()) {
            // Updated! Prev: compare gscores and overwrite gscore; Curr: compare fscores and overwrite all node infos
            if (cost_instead_of_score && (cur_fscore >= iterFound->fscore) || 
                !cost_instead_of_score && (cur_fscore <= iterFound->fscore)) {
                workload_recorder.node_pruned_cnt_gscore += 1;      // prune this node
            } else {    // update gscore, fscore, and tbinfo
                (*target_index_by_coord)[nbr_bin]->modify(iterFound, typename AstarMultiIndex<NN>::ChangeGscoreFscoreTBinfoAffineGap(cur_gscore, cur_fscore, nbr_tb_info));
                workload_recorder.node_pruned_cnt_gscore += 1;      // prune this node
            }
        } else {    // insert a new node
            typename AstarMultiIndex<NN>::OpenListNodeAffineGap nbr_node = {
                nbr_crd_tb_info,
                cur_fscore,
                cur_gscore
            };
            (*target_index_by_coord)[nbr_bin]->insert(nbr_node);
            workload_recorder.open_list_cnt += 1;
            workload_recorder.cur_open_list_cnt += 1;
            bin_op_cnt[nbr_bin].first += 1;     // insertion + 1
        }
    } else {
        // Detect duplicated nodes in OpenList
        int nbr_bin = (int)floor((cur_gscore - bin_score_offset) / bin_score_thres);

        if (nbr_bin < 0) nbr_bin = 0;
        if (nbr_bin >= bin_cnt) nbr_bin = bin_cnt - 1;

        bool found_in_higher_bin = false;    // First find the nbr_node in higher bins
        if (cost_instead_of_score) {        // from nbr_bin - 1 to 0
            for (int tmp_bin_idx = nbr_bin - 1; tmp_bin_idx >= 0; --tmp_bin_idx) {
                auto iterFound = (*target_index_by_coord)[tmp_bin_idx]->find(nbr_crd_tb_info);
                if (iterFound != (*target_index_by_coord)[tmp_bin_idx]->end()) {     // cur_gscore cannot be larger than the previous one
                    workload_recorder.node_pruned_cnt_gscore += 1;
                    found_in_higher_bin = true;
                    break;
                }
            }
        } else {                            // from nbr_bin + 1 to bin_cnt - 1
            for (int tmp_bin_idx = nbr_bin + 1; tmp_bin_idx < bin_cnt; ++tmp_bin_idx) {
                auto iterFound = (*target_index_by_coord)[tmp_bin_idx]->find(nbr_crd_tb_info);
                if (iterFound != (*target_index_by_coord)[tmp_bin_idx]->end()) {     // cur_gscore cannot be larger than the previous one
                    workload_recorder.node_pruned_cnt_gscore += 1;
                    found_in_higher_bin = true;
                    break;
                }
            }
        }
        if (found_in_higher_bin == false) {  // If the node is not in bin[nbr_bin + 1], search bin[nbr_bin]
            auto iterFound = (*target_index_by_coord)[nbr_bin]->find(nbr_crd_tb_info);
            if (iterFound != (*target_index_by_coord)[nbr_bin]->end()) {
                // Updated! Prev: compare gscores and overwrite gscore; Curr: compare fscores and overwrite all node infos
                if (cost_instead_of_score && (cur_fscore >= iterFound->fscore) || 
                    !cost_instead_of_score && (cur_fscore <= iterFound->fscore)) {
                    workload_recorder.node_pruned_cnt_gscore += 1;      // prune this node
                } else {    // update gscore, fscore, and tbinfo
                    (*target_index_by_coord)[nbr_bin]->modify(iterFound, typename AstarMultiIndex<NN>::ChangeGscoreFscoreTBinfoAffineGap(cur_gscore, cur_fscore, nbr_tb_info));
                    workload_recorder.node_pruned_cnt_gscore += 1;      // prune this node
                }
            } else {    // insert a new node
                typename AstarMultiIndex<NN>::OpenListNodeAffineGap nbr_node = {
                    nbr_crd_tb_info,
                    cur_fscore,
                    cur_gscore
                };
                (*target_index_by_coord)[nbr_bin]->insert(nbr_node);
                workload_recorder.open_list_cnt += 1;
                workload_recorder.cur_open_list_cnt += 1;
                bin_op_cnt[nbr_bin].first += 1;     // insertion + 1
            }
        }
        // prune this crd in lower bin
        if (cost_instead_of_score) {
            if (nbr_bin + 1 < bin_cnt) {
                auto iterFound = (*target_index_by_coord)[nbr_bin + 1]->find(nbr_crd_tb_info);
                if (iterFound != (*target_index_by_coord)[nbr_bin + 1]->end()) {
                    workload_recorder.node_pruned_cnt_gscore += 1;      // prune this node
                    (*target_index_by_coord)[nbr_bin + 1]->erase(iterFound);
                    workload_recorder.cur_open_list_cnt -= 1;
                    bin_op_cnt[nbr_bin + 1].second += 1;
                }
            }
        } else {
            if (nbr_bin > 0) {
                auto iterFound = (*target_index_by_coord)[nbr_bin - 1]->find(nbr_crd_tb_info);
                if (iterFound != (*target_index_by_coord)[nbr_bin - 1]->end()) {
                    workload_recorder.node_pruned_cnt_gscore += 1;      // prune this node
                    (*target_index_by_coord)[nbr_bin - 1]->erase(iterFound);
                    workload_recorder.cur_open_list_cnt -= 1;
                    bin_op_cnt[nbr_bin - 1].second += 1;
                }
            }
        }
    }
    

}


template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
void AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::preprocessing_MSA_with_PSA_and_greedy() {
    Timer MSA_timer;
    MSA_timer.start();
    pairwise_heuristic->compute_all_PSA();        // pre-compute all pairwise alignments
    workload_recorder.all_PSA_time = MSA_timer.elapsed(false);

    printf("all-PSA done!\n");

    /* Initialize the tmp_msa_gscore with DFS */
    bool greedy_reverse_seq = false;
    compute_greedy(nn, greedy_reverse_seq);
    workload_recorder.greedy_MSA_time = MSA_timer.elapsed(false);
    printf("Greedy MSA done!\n");
    /* End of Initialize the tmp_msa_gscore with DFS */
}

template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
void AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::compute_anytime_Astar() {
    const int dim = nn;
    NodeCoord cur_crd;
    for (int i = 0; i < NN; ++i) cur_crd[i] = 0;

    GapLen cur_tb_info;
    for (int i = 0; i < NN; ++i) cur_tb_info[i] = 0;
    typename AstarMultiIndex<NN>::CoordTBInfo cur_crd_tb_info = {cur_crd, cur_tb_info};

    // Move the previous closed list results into recursive_closed_list
    recursive_closed_list_gap_len.clear();
    for (auto closed_iter = closed_list_gap_len.begin(); closed_iter != closed_list_gap_len.end(); ++closed_iter) {
        recursive_closed_list_gap_len[closed_iter->first] = closed_iter->second;
    }
    closed_list_gap_len.clear();
    
    STYPE highest_score = 0;

    CRDTYPE pairwise_heuristic_crd[NN];
    for (int i = 0; i < dim; ++i) pairwise_heuristic_crd[i] = seq_lens[i] - cur_crd[i];
    bool access_reversed_mats = true;
    DIRTYPE pairwise_gap_len[NN];
    for (int i = 0; i < dim; ++i) pairwise_gap_len[i] = cur_tb_info[i];
    highest_score = pairwise_heuristic->get_sum_of_scores(dim, access_reversed_mats, pairwise_heuristic_crd, pairwise_gap_len);

    bin_score_thres = (highest_score - bin_score_offset) / bin_cnt;

    printf("h-score of the source node = %f, bin_score_offset = %g, bin_score_thres = %g\n", highest_score, bin_score_offset, bin_score_thres);

    typename AstarMultiIndex<NN>::OpenListNodeAffineGap source_node = {
        cur_crd_tb_info,
        0 + highest_score,
        0
    };
    index_by_fscore_arr[0]->insert(source_node);         // bin[0] because gscore = 0

    workload_recorder.open_list_cnt += 1;
    workload_recorder.cur_open_list_cnt += 1;

    bin_op_cnt[0].first = 1;
    
    bool is_opt = false;        // if the current alignment result is optimal
    while (is_opt == false && time_limit_timer->check_time_limit() == false && !open_list_arr_empty()) {
        for (int bin_level = 0; bin_level < bin_cnt; ++bin_level) {     // perform Anytime Column Search iterations
            int beam = 0;
            while (is_opt == false && time_limit_timer->check_time_limit() == false && !open_list_arr_empty(bin_level) && beam < beam_width) {
                print_progress_per_1M_iter();
                workload_recorder.iter_cnt += 1;
                bool specific_bin = true;   // expand a node in bin[bin_level]
                bool reverse_seq = false;   // we are not computing recursive MSA
                is_opt = expand_node(dim, specific_bin, bin_level, reverse_seq);
                beam += 1;
            }
        }
        for (int astar_iter = 0; astar_iter < astar_iter_cnt; ++astar_iter) {       // perform normal A* Search iterations
            if (is_opt || time_limit_timer->check_time_limit() || open_list_arr_empty()) break;
            print_progress_per_1M_iter();
            workload_recorder.iter_cnt += 1;
            bool specific_bin = false;      // expand a node in bin_arr by searching one with the largest fscore
            bool reverse_seq = false;       // we are not computing recursive MSA
            is_opt = expand_node(dim, specific_bin, 0, reverse_seq);
        }        
        // print the iter count when memory bound is reached for the first time
        if (memory_bound_trigger == true && memory_bound_first_iter == 0) 
            memory_bound_first_iter = workload_recorder.iter_cnt;

        // prune nodes in the closed list 
        if (get_closed_list_size() >= max_closed_size) {
            bool access_reversed_mats = true;
            printf("The closed list consumed all available memory!\n");
            int prunable_cnt = prune_closed_list(access_reversed_mats, dim);
            if (prunable_cnt == 0) break;
        }
    }

    print_workload_report_after_search(NN);
}


// If the current alignment result is optimal, return true. Otherwise, return false
// reverse_seq: read the symbol in each sequences from the end instead of the beginning to generate the recursive closed list
// specific_bin = false: find the node with the largest fscore in all bins

template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
bool AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::expand_node(int dim, bool specific_bin, int bin_level, bool reverse_seq) {
    STYPE ext_pen = blosum_table.get_ext_penalty(), open_pen = blosum_table.get_open_penalty();

    std::vector<OLFscore*>* target_index_by_fscore = &index_by_fscore_arr;
    std::vector<OLCoord*>* target_index_by_coord = &index_by_coord_arr;

    bool access_reversed_mats = !reverse_seq;

    // search the node with the highest fscore and update the bin_level
    if (bin_cnt > 1 && specific_bin == false) {
        // update the bin_level with the bin_idx containing the highest fscore
        get_highest_fscore(bin_level);
    }
    
    typename AstarMultiIndex<NN>::CoordTBInfo cur_crd_tb_info;
    STYPE local_gscore, local_fscore;

    auto cur_node = (*target_index_by_fscore)[bin_level]->begin();
    cur_crd_tb_info.crd = cur_node->crd_tb_info.crd;
    cur_crd_tb_info.tb_info = cur_node->crd_tb_info.tb_info;
    local_gscore = cur_node->gscore;
    local_fscore = cur_node->fscore;
    
    (*target_index_by_fscore)[bin_level]->erase(cur_node);
    workload_recorder.cur_open_list_cnt -= 1;
    bin_op_cnt[bin_level].second += 1;
    
    bool target_reached = is_crd_out_of_bound(cur_crd_tb_info.crd, dim);

    // find the cur_node in ClosedList

    // first find the coord
    auto closed_list_gap_len_iter = closed_list_gap_len.find(cur_crd_tb_info.crd);
    bool write_into_closed_list = true;
    // if we found the node, update the gscore
    if (closed_list_gap_len_iter != closed_list_gap_len.end()) {
        // then find the gap_len and update the gscore in ClosedList
        auto gap_len_score_iter = closed_list_gap_len_iter->second.find(cur_crd_tb_info.tb_info);
        if (gap_len_score_iter != closed_list_gap_len_iter->second.end()) {

            // gap_len_score_iter->second.second is reexp_fscore
            if (gap_len_score_iter->second.second == init_reexp_fscore) {
                // Case 1: Not reexpanding the cur node, return when cur_gscore <= closed_gscore
                if (cost_instead_of_score && (local_gscore >= gap_len_score_iter->second.first) || 
                    !cost_instead_of_score && (local_gscore <= gap_len_score_iter->second.first)) {
                    return false;
                }
            } else {
                // Case 2: Reexpanding the cur node, return when cur_gscore < closed_gscore
                if (cost_instead_of_score && (local_gscore > gap_len_score_iter->second.first) || 
                    !cost_instead_of_score && (local_gscore < gap_len_score_iter->second.first)) {
                    return false;
                }
            }
        }
    }

    // if we didn't find the node or cur_node is better, insert it into the closed list
    if (write_into_closed_list) {
        closed_list_gap_len[cur_crd_tb_info.crd][cur_crd_tb_info.tb_info].first = local_gscore;
        closed_list_gap_len[cur_crd_tb_info.crd][cur_crd_tb_info.tb_info].second = init_reexp_fscore;
        workload_recorder.closed_list_cnt += 1;
    }

    if (target_reached) {
        if (cost_instead_of_score && (local_gscore < tmp_msa_gscore) ||
            !cost_instead_of_score && (local_gscore > tmp_msa_gscore)) {
            tmp_msa_gscore = local_gscore;      // update the score lower bound with the beam search result 
            // don't call backtrack when dim < nn, i.e., during recursive MSA with a lower dimension
            if (dim == nn) {
                printf("Before backtracking\n");
                tmp_msa_alignment_length = backtrack_affine(tmp_msa_result, cur_crd_tb_info, false);
                file_writer->write(time_limit_timer->get_time_from_start(), tmp_msa_result, tmp_msa_alignment_length, local_gscore, workload_recorder.iter_cnt);
                printf("Found a better MSA (SP score = %f)!\n", local_gscore);
                
                // write the TC score
                if (blosum_table.is_ref_MSA_set()) {
                    STYPE tc_score = blosum_table.compute_total_column_score(tmp_msa_result, nn, tmp_msa_alignment_length);
                    STYPE sp_accuracy = blosum_table.compute_sum_of_pairs_accuracy(tmp_msa_result, nn, tmp_msa_alignment_length);

                    printf("TC score = %f, SP accuracy = %f\n", tc_score, sp_accuracy);
                    file_writer->write_TC_score(tc_score);
                    file_writer->write_SP_accuracy(sp_accuracy);
                }

                prune_closed_list(access_reversed_mats, dim);
            }
        }

        bool is_opt_ret = true;
        STYPE global_best_result = get_highest_fscore();
        if (cost_instead_of_score && (local_gscore > global_best_result) || 
            !cost_instead_of_score && (local_gscore < global_best_result)) is_opt_ret = false;

        return is_opt_ret;
    }

    DIRTYPE dir_cnt = (1 << dim);       // direction counts

    // pre-load the match / mismatch scores of each sequence pair
    STYPE mm_scores[NN * NN];
    get_match_mismatch_scores(mm_scores, cur_crd_tb_info.crd, dim, reverse_seq);

    for (DIRTYPE dir = 1; dir < dir_cnt; ++dir) {
        NodeCoord nbr_crd (cur_crd_tb_info.crd); // first set nbr_crd = cur_crd, then check if we should "+= 1" to each dim

        bool oob_flag = false;  // node out-of-bound flag
        for (int bit = 0; bit < dim; ++bit) {
            if ((dir>>bit) & 0x01) {
                nbr_crd[bit] += 1;                        // if i == 0b011, then dir_vec = [1, 1, 0]
                if (nbr_crd[bit] > seq_lens[bit]) oob_flag = true;      // don't go beyond the last cell
            }
        }
        if (oob_flag) continue;

        // 1. f-score = cur_g-score + match/mismatch/gap + h-score
        // 2. if f-score >= tmp_msa_gscore, push into open_list
        STYPE cur_gscore = local_gscore;        // cur_gscore = pre_gscore + the score of match/mismatch/gap, need check: use local gscore instead of global gscore 
        STYPE cur_fscore = 0;                   // cur_fscore = cur_gscore + cur_hscore

        // check all pairs of sequences
        for (int seq_i = 0; seq_i < dim - 1; ++seq_i) {
            for (int seq_j = seq_i + 1; seq_j < dim; ++seq_j) {
                // compute and update the score of match/mismatch/gap
                if (nbr_crd[seq_i] != cur_crd_tb_info.crd[seq_i] && nbr_crd[seq_j] != cur_crd_tb_info.crd[seq_j]) {     // check the score table
                    cur_gscore += mm_scores[seq_i * NN + seq_j];
                } else {        // gap(s) in the two sequences
                    // There could be gaps in both sequences, so we need to check the gap length
                    bool gap_in_prev_seq_i = false, gap_in_prev_seq_j = false;
                    
                    int gap_len_seq_i = cur_crd_tb_info.tb_info[seq_i], gap_len_seq_j = cur_crd_tb_info.tb_info[seq_j];
                    if (gap_len_seq_i > gap_len_seq_j) gap_in_prev_seq_i = true;
                    if (gap_len_seq_j > gap_len_seq_i) gap_in_prev_seq_j = true;
                    
                    if (nbr_crd[seq_i] == cur_crd_tb_info.crd[seq_i] && nbr_crd[seq_j] == cur_crd_tb_info.crd[seq_j]) {     // two gaps here
                    } else if (nbr_crd[seq_i] == cur_crd_tb_info.crd[seq_i]) {    // gap in seq i
                        cur_gscore -= (!gap_in_prev_seq_i || cur_crd_tb_info.crd[seq_i] == 0 && cur_crd_tb_info.crd[seq_j] == 0) * open_pen + ext_pen;
                    } else {                                    // gap in seq j
                        cur_gscore -= (!gap_in_prev_seq_j || cur_crd_tb_info.crd[seq_i] == 0 && cur_crd_tb_info.crd[seq_j] == 0) * open_pen + ext_pen;
                    }
                }
            }
        }

        typename AstarMultiIndex<NN>::CoordTBInfo nbr_crd_tb_info = {nbr_crd, cur_crd_tb_info.tb_info};

        // update the tb_info
        for (int seq_idx = 0; seq_idx < dim; ++seq_idx) {
            if ((dir >> seq_idx) & 0x01) 
                nbr_crd_tb_info.tb_info[seq_idx] = 0;       // non-gap in seq[seq_idx]
            else
                nbr_crd_tb_info.tb_info[seq_idx] += 1;
        }

        // find nbr node in ClosedList and compare the gscores
        bool use_reexp_fscore_instead = false;

        bool is_nbr_duplicated = closed_list_duplication_detection(use_reexp_fscore_instead, nbr_crd, nbr_crd_tb_info.tb_info, cur_gscore, cur_fscore);
        if (is_nbr_duplicated) continue;

        if (use_reexp_fscore_instead == false) {
            cur_fscore += cur_gscore;
            // compute the hscore separately
            STYPE cur_hscore = score_upper_bound_recursive(dim, access_reversed_mats, nbr_crd, nbr_crd_tb_info.tb_info);
            cur_fscore += cur_hscore;
        }

        // if cur_fscore is worse than tmp_msa_gscore, don't insert it into the open_list and closed_list
        if (cost_instead_of_score && (cur_fscore > tmp_msa_gscore + __FLT_EPSILON__) || 
            !cost_instead_of_score && (cur_fscore + __FLT_EPSILON__ < tmp_msa_gscore)) {
            workload_recorder.node_pruned_cnt_hscore += 1;
            continue;
        }

        open_list_insertion(cur_gscore, cur_fscore, target_index_by_coord, nbr_crd, nbr_crd_tb_info.tb_info);
    }

    bool is_opt_ret = true;
    STYPE global_best_result = get_highest_fscore();
    if (cost_instead_of_score && (tmp_msa_gscore > global_best_result) ||
        !cost_instead_of_score && (tmp_msa_gscore < global_best_result)) is_opt_ret = false;


    // memory bound
    // check once per 100 iterations
    if (workload_recorder.iter_cnt % 100 == 0 && is_opt_ret == false && enable_memory_bound == true) {

        check_memory_bound_trigger();

        while (memory_bound_trigger && check_memory_thres(false, get_closed_list_size())) {
            /*
                1. Remove the worst node from the lowest non-empty bin. Update its parent's reexp-fscore in Closed List & fscore in Open List
                    1-1. Specifically, if the fscore of the current node is better than the parent's reexp-fscore/fscore, update the score
                    1-2. The source node cannot be removed. Find the second worst node
                2. Insert its parent into the Open List (check duplication)

                Note:
                1. Reset the reexp-fscore when we reexpand a node (if closed_list.find(cur_node) then reset())
                    1-1. Reason: So that its children won't be pruned unexpectedly when trying to access the update the reexp-fscore again.
                2. When computing the fscore of a neighbor, use its reexp-fscore instead if it's found in the Closed List
            */

            int non_empty_bin_idx;
            // if last_node is the source node, remove the second last node instead
            NodeCoord src_crd;
            for (int i = 0; i < NN; ++i) src_crd[i] = 0;
            auto last_node_iter = (*target_index_by_fscore)[0]->begin();        // updated in either of the if-block

            non_empty_bin_idx = 0;
            // find the lowest non-empty bin (Notice: we skip the source node)
            while (non_empty_bin_idx < bin_cnt && (*target_index_by_fscore)[non_empty_bin_idx]->empty())
                non_empty_bin_idx += 1;
            if (non_empty_bin_idx == bin_cnt)     // open list is empty
                return is_opt_ret;

            last_node_iter = prev((*target_index_by_fscore)[non_empty_bin_idx]->end());

            NodeCoord last_node_crd = last_node_iter->get_crd();

            if (last_node_crd == src_crd) {
                // if source is the only node, go to the next bin
                if ((*target_index_by_fscore)[non_empty_bin_idx]->size() == 1) {
                    non_empty_bin_idx += 1;
                    while (non_empty_bin_idx < bin_cnt && (*target_index_by_fscore)[non_empty_bin_idx]->empty())
                        non_empty_bin_idx += 1;
                    // access the last node in non_empty_bin_idx
                    if (non_empty_bin_idx == bin_cnt) {
                        return is_opt_ret;
                    }
                    last_node_iter = prev((*target_index_by_fscore)[non_empty_bin_idx]->end());
                }

                last_node_crd = last_node_iter->get_crd();
            }

            // Calculate the parent's crd and erase the last_node
            STYPE last_node_fscore = last_node_iter->fscore;
            GapLen last_node_tb_info = last_node_iter->get_tb_info();            // a single DIRTYPE variable instead of a vector
            (*target_index_by_fscore)[non_empty_bin_idx]->erase(last_node_iter);
            workload_recorder.cur_open_list_cnt -= 1;
            bin_op_cnt[non_empty_bin_idx].second += 1;     // erasion + 1

            // find a qualified parent in this closed list node and update the reexp_fscore
            // Condition 1: parent.gap_len + delta_gap_len = cur.gap_len
            // Condition 2: parent.gscore + delta_gscore >= cur.gscore
            GapLen parent_gap_len = find_parent_in_closed_list(last_node_crd, last_node_tb_info, last_node_iter->gscore);
            if (parent_gap_len[0] == -1) {
                printf("No qualified parent in Astar memory bound extension!\n");
                break;
            }

            NodeCoord parent_crd = last_node_crd;
            for (int i = 0; i < nn; ++i)
                if (is_gap(last_node_tb_info, i) == false) parent_crd[i] -= 1;

            // closed_list_gap_len[parent_crd][parent_gap_len].second is reexp_fscore
            if (cost_instead_of_score && (last_node_fscore < closed_list_gap_len[parent_crd][parent_gap_len].second) || 
                !cost_instead_of_score && (last_node_fscore > closed_list_gap_len[parent_crd][parent_gap_len].second)) {
                // potential of this node: should contain the info of the REMOVED children
                // When a child is pruned,
                // Case 1: found its parent in closed, child has better fscore
                // Case 1-1: found its parent in open, child has better fscore
                //      Update both
                // Case 1-2: found its parent in open, child has worse fscore 
                //      Update the closed list node
                // Case 1-3: didn't find its parent in open
                //      Update the closed list node and reinsert the parent
                // Case 2: found its parent in closed, child has worse fscore
                // Case 2-1: found its parent in open, child has better fscore
                //      Update the open list node
                // Case 2-2: found its parent in open, child has worse fscore
                //      Do nothing
                // Case 2-3: didn't find its parent in open
                //      Reinsert? (duplication?)
                closed_list_gap_len[parent_crd][parent_gap_len].second = last_node_fscore;
            }

            // Find the parent in Open List and update the fscore (only search the lower bins)
            bool found_parent_in_open_list = false;
            typename AstarMultiIndex<NN>::CoordTBInfo parent_crd_tb_info = {parent_crd, parent_gap_len};

            // score: parent_node <= cur_node, step = -1
            // cost: parent_node <= cur_node, step = -1
            for (int tmp_bin_idx = non_empty_bin_idx; tmp_bin_idx >= 0; --tmp_bin_idx) {
                auto parent_open_list_iter = (*target_index_by_coord)[tmp_bin_idx]->find(parent_crd_tb_info);
                if (parent_open_list_iter != (*target_index_by_coord)[tmp_bin_idx]->end()) {

                    if (cost_instead_of_score && (last_node_fscore < parent_open_list_iter->fscore) || 
                        !cost_instead_of_score && (last_node_fscore > parent_open_list_iter->fscore)) {
                        (*target_index_by_coord)[tmp_bin_idx]->modify(parent_open_list_iter, typename AstarMultiIndex<NN>::ChangeFscoreAffineGap(last_node_fscore));
                    }
                    
                    STYPE closed_list_gscore = closed_list_gap_len[parent_crd][parent_gap_len].first;
                    // compare the gscores
                    if (cost_instead_of_score && (closed_list_gscore >= parent_open_list_iter->gscore) || 
                        !cost_instead_of_score && (closed_list_gscore <= parent_open_list_iter->gscore)) {      // if closed_list_gscore <= gscore, skip; else update gscore
                        // do nothing
                    } else {    // update gscore
                        (*target_index_by_coord)[tmp_bin_idx]->modify(parent_open_list_iter, typename AstarMultiIndex<NN>::ChangeGscoreAffineGap(closed_list_gscore));
                    }

                    found_parent_in_open_list = true;
                    break;
                }

                // the parent can only appear in [non_empty_bin_idx - dim, non_empty_bin_idx - 1]
                if (sum_of_crd_as_level && tmp_bin_idx == non_empty_bin_idx - dim)
                    break;
            }

            // Insert the parent into the Open List
            if (found_parent_in_open_list == false) {
                typename AstarMultiIndex<NN>::OpenListNodeAffineGap reexp_parent_node = {
                    parent_crd_tb_info,
                    closed_list_gap_len[parent_crd][parent_gap_len].second,     // reexp_fscore
                    closed_list_gap_len[parent_crd][parent_gap_len].first       // gscore
                };
                
                // prune the parent node with tmp_msa_gscore
                if (cost_instead_of_score && (reexp_parent_node.fscore > tmp_msa_gscore + __FLT_EPSILON__) || 
                    !cost_instead_of_score && (reexp_parent_node.fscore + __FLT_EPSILON__ < tmp_msa_gscore)) {
                    workload_recorder.node_pruned_cnt_hscore += 1;
                    continue;
                }

                int parent_bin;
                if (sum_of_crd_as_level) {
                    parent_bin = 0;
                    for (int i = 0; i < NN; ++i) parent_bin += parent_crd[i];
                } else {
                    parent_bin = (int)floor((reexp_parent_node.gscore - bin_score_offset) / bin_score_thres);
                    if (parent_bin < 0) parent_bin = 0;
                    if (parent_bin >= bin_cnt) parent_bin = bin_cnt - 1;
                }

                (*target_index_by_fscore)[parent_bin]->insert(reexp_parent_node);
                workload_recorder.open_list_cnt += 1;
                workload_recorder.cur_open_list_cnt += 1;
                bin_op_cnt[parent_bin].first += 1;     // insertion + 1
            }
            

        }
    }
    return is_opt_ret;
}


// initialize the tmp_msa_gscore with DFS
// record the lowest gscore and set the bin_score_offset
template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
void AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::compute_greedy(int dim, bool reverse_seq) {
    STYPE ext_pen = blosum_table.get_ext_penalty(), open_pen = blosum_table.get_open_penalty();
    // only stores the coordinate and gscore, generate the MSA results (tmp_msa_result) during DFS

    // we only use the first dim elements so seems to be safe
    NodeCoord cur_node, next_node;
    for (int i = 0; i < NN; ++i) cur_node[i] = 0;
    for (int i = 0; i < NN; ++i) next_node[i] = 0;

    GapLen cur_tb_info, next_tb_info;       // only used when computing affine gap penalty
    for (int i = 0; i < NN; ++i) cur_tb_info[i] = 0;
    for (int i = 0; i < NN; ++i) next_tb_info[i] = 0;

    STYPE local_gscore = 0, next_gscore;
    int cur_alignment_length = 0;
    STYPE lowest_gscore = INT_MAX;

    bool target_reached = false;
    while (target_reached == false) {
        // traverse all children
        STYPE best_child_score = cost_instead_of_score ? INT_MAX : INT_MIN;

        DIRTYPE dir_cnt = (1 << dim);       // direction counts

        // pre-load the match / mismatch scores of each sequence pair
        STYPE mm_scores[NN * NN];
        get_match_mismatch_scores(mm_scores, cur_node, dim, reverse_seq);
        
        for (DIRTYPE dir = 1; dir < dir_cnt; ++dir) {
            // we only use the first dim elements so seems to be safe
            NodeCoord nbr_crd (cur_node);
            GapLen nbr_tb_info;       // only used when computing affine gap penalty
            for (int i = 0; i < NN; ++i) nbr_tb_info[i] = 0;
            
            bool oob_flag = false;  // node out-of-bound flag
            for (int bit = 0; bit < dim; ++bit) {
                if ((dir>>bit) & 0x01) {
                    nbr_crd[bit] += 1;                        // if i == 0b011, then dir_vec = [1, 1, 0]
                    if (nbr_crd[bit] > seq_lens[bit]) oob_flag = true;      // don't go beyond the last cell
                }
            }
            if (oob_flag) continue;

            STYPE cur_gscore = local_gscore;
            STYPE cur_fscore = 0;

            // check all pairs of sequences
            for (int seq_i = 0; seq_i < dim - 1; ++seq_i) {
                for (int seq_j = seq_i + 1; seq_j < dim; ++seq_j) {
                    // compute and update the score of match/mismatch/gap
                    if (nbr_crd[seq_i] != cur_node[seq_i] && nbr_crd[seq_j] != cur_node[seq_j]) {     // check the score table
                        cur_gscore += mm_scores[seq_i * NN + seq_j];
                    } else {        // gap(s) in the two sequences
                        if (nbr_crd[seq_i] == cur_node[seq_i] && nbr_crd[seq_j] == cur_node[seq_j]) {     // two gaps here
                            // do nothing
                        } else if (nbr_crd[seq_i] == cur_node[seq_i]) {    // gap in seq i
                            if (use_linear_gap_penalty)
                                cur_gscore -= ext_pen;
                            else {
                                bool gap_in_prev_seq_i = cur_tb_info[seq_i] > cur_tb_info[seq_j];
                                cur_gscore -= (!gap_in_prev_seq_i || cur_node[seq_i] == 0 && cur_node[seq_j] == 0) * open_pen + ext_pen;
                            }
                        } else {    // gap in seq j
                            if (use_linear_gap_penalty)
                                cur_gscore -= ext_pen;
                            else {
                                bool gap_in_prev_seq_j = cur_tb_info[seq_i] < cur_tb_info[seq_j];
                                cur_gscore -= (!gap_in_prev_seq_j || cur_node[seq_i] == 0 && cur_node[seq_j] == 0) * open_pen + ext_pen;
                            }
                        }
                    }
                }
            }

            if (cur_gscore < lowest_gscore) lowest_gscore = cur_gscore;

            if (use_linear_gap_penalty == false) {
                // update tb_info (gap lens)
                for (int seq_idx = 0; seq_idx < dim; ++seq_idx) {
                    if ((dir >> seq_idx) & 0x01) 
                        nbr_tb_info[seq_idx] = 0;       // non-gap in seq[seq_idx]
                    else
                        nbr_tb_info[seq_idx] = cur_tb_info[seq_idx] + 1;
                }
            }

            bool access_reversed_mats = !reverse_seq;
            STYPE cur_hscore = 0.0;

            CRDTYPE pairwise_heuristic_crd[NN];
            for (int i = 0; i < dim; ++i) 
                pairwise_heuristic_crd[i] = seq_lens[i] - nbr_crd[i];

            if (use_linear_gap_penalty) {
                cur_hscore = pairwise_heuristic->get_sum_of_scores_linear_gap(dim, access_reversed_mats, pairwise_heuristic_crd);
            } else {
                DIRTYPE pairwise_gap_len[NN];
                for (int i = 0; i < dim; ++i) pairwise_gap_len[i] = nbr_tb_info[i];
                cur_hscore = pairwise_heuristic->get_sum_of_scores(dim, access_reversed_mats, pairwise_heuristic_crd, pairwise_gap_len);
            }

            cur_fscore = cur_gscore + cur_hscore;

            if (cost_instead_of_score && (cur_fscore < best_child_score) || !cost_instead_of_score && (cur_fscore > best_child_score)) {
                best_child_score = cur_fscore;
                next_gscore = cur_gscore;
                next_node = nbr_crd;
                if (use_linear_gap_penalty == false)
                    next_tb_info = nbr_tb_info;
            }
            
        }
        

        // write a new column into tmp_msa_result
        for (int i = 0; i < dim; ++i) {
            bool is_gap = (next_node[i] == cur_node[i]);
            tmp_msa_result[i][cur_alignment_length] = is_gap ? GAP : sequences[i][next_node[i] - 1];
        }
        cur_alignment_length += 1;
        cur_node = next_node;
        local_gscore = next_gscore;
        if (use_linear_gap_penalty == false)
            cur_tb_info = next_tb_info;


        target_reached = is_crd_out_of_bound(cur_node, dim);
    }

    if (cost_instead_of_score && (local_gscore < tmp_msa_gscore) ||
        !cost_instead_of_score && (local_gscore > tmp_msa_gscore)) {
        tmp_msa_gscore = local_gscore;    // update the score lower bound with the beam search result 

        tmp_msa_alignment_length = cur_alignment_length;
        file_writer->write(time_limit_timer->get_time_from_start(), tmp_msa_result, tmp_msa_alignment_length, local_gscore, 0);
            
        // write the TC score
        if (blosum_table.is_ref_MSA_set()) {
            STYPE tc_score = blosum_table.compute_total_column_score(tmp_msa_result, nn, tmp_msa_alignment_length);
            STYPE sp_accuracy = blosum_table.compute_sum_of_pairs_accuracy(tmp_msa_result, nn, tmp_msa_alignment_length);
            printf("TC score = %f, SP accuracy = %f\n", tc_score, sp_accuracy);
            file_writer->write_TC_score(tc_score);
            file_writer->write_SP_accuracy(sp_accuracy);
        }
    }
    printf("In greedy MSA, local_gscore = %f, lowest_gscore = bin_score_offset = %f\n", local_gscore, lowest_gscore);
    if (cost_instead_of_score == false) 
        bin_score_offset = lowest_gscore;

}


// store closed list nodes with gap length
template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
void AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::compute_recursive_Astar(int seq_cnt, bool reverse_seq) {
    // temporarily change the hyper-parameters: only use Astar in recursive MSA
    int prev_beam_width = beam_width, prev_bin_cnt = bin_cnt;
    beam_width = 0; bin_cnt = 1; 

    /* --- Open list of nn-D MSA --- */
    // initialize open_list_arr, deallocation seems unnecessary
    open_list_arr.clear();
    index_by_fscore_arr.clear();
    index_by_coord_arr.clear();

    open_list_arr = std::vector<OLMultiIdx>(bin_cnt);
    for (int i = 0; i < bin_cnt; ++i) {
        index_by_fscore_arr.push_back(&open_list_arr[i].template get<typename AstarMultiIndex<NN>::IndexByFscore>());
        index_by_coord_arr.push_back(&open_list_arr[i].template get<typename AstarMultiIndex<NN>::IndexByCoord>());
    }

    // we don't need the alignment but the closed list nodes
    const int dim = seq_cnt;
    NodeCoord cur_crd;
    for (int i = 0; i < NN; ++i) cur_crd[i] = 0;

    /* --- diff --- */
    GapLen cur_tb_info;
    for (int i = 0; i < NN; ++i) cur_tb_info[i] = 0;
    typename AstarMultiIndex<NN>::CoordTBInfo cur_crd_tb_info = {cur_crd, cur_tb_info};

    // Move the previous closed list results into recursive_closed_list
    recursive_closed_list_gap_len.clear();
    for (auto closed_iter = closed_list_gap_len.begin(); closed_iter != closed_list_gap_len.end(); ++closed_iter) {
        recursive_closed_list_gap_len[closed_iter->first] = closed_iter->second;
    }
    closed_list_gap_len.clear();

    // if reverse_seq == false, access reversed_mats
    bool access_reversed_mats = !reverse_seq;
    STYPE highest_score = score_upper_bound_recursive(dim, access_reversed_mats, cur_crd, cur_tb_info);
    /* --- diff --- */

    bin_score_thres = (highest_score - bin_score_offset) / bin_cnt;

    printf("In compute recursive Astar, dim = %d, h-score of the source node = %f, bin_score_offset = %g, bin_score_thres = %g\n", dim, highest_score, bin_score_offset, bin_score_thres);

    /* --- diff --- */
    typename AstarMultiIndex<NN>::OpenListNodeAffineGap source_node = {
        cur_crd_tb_info,
        0 + highest_score,
        0
    };
    index_by_fscore_arr[0]->insert(source_node);         // bin[0] because gscore = 0
    workload_recorder.open_list_cnt += 1;
    workload_recorder.cur_open_list_cnt += 1;

    bin_op_cnt[0].first = 1;
    
    bool is_opt = false;        // if the current alignment result is optimal
    while (is_opt == false && time_limit_timer->check_time_limit() == false && get_closed_list_size() < max_closed_size && !open_list_arr_empty()) {
        for (int astar_iter = 0; astar_iter < astar_iter_cnt; ++astar_iter) {       // perform normal A* Search iterations
            if (is_opt || time_limit_timer->check_time_limit() || open_list_arr_empty()) break;
            print_progress_per_1M_iter();
            workload_recorder.iter_cnt += 1;
            bool specific_bin = false;      // expand a node in bin_arr by searching one with the largest fscore

            /* --- diff --- */
            is_opt = expand_node(dim, specific_bin, 0, reverse_seq);
        }
    }

    print_workload_report_after_search(dim);

    // restore the hyper-parameters
    beam_width = prev_beam_width; bin_cnt = prev_bin_cnt;
}


// store_gap_len = true
// write the result into msa_result and return alignment length
// print_path == true: only print the backtrack path in terminal; print_path == false: output the alignment
template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
int AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::backtrack_affine(char **&msa_result, typename AstarMultiIndex<NN>::CoordTBInfo crd_tb_info, bool print_path) {
    NodeCoord seq_indices = crd_tb_info.crd;        // indices of the TB_mat, use --idx{i} when accessing sequences[][]

    int res_idx = 0;                                // counter of the alignment length
    bool all_zero_indices = true;
    for (int i = 0; i < nn; ++i)
        if (seq_indices[i] != 0) {
            all_zero_indices = false;
            break;
        }

    GapLen gap_len = crd_tb_info.tb_info;

    if (print_path) {
        printf("Inverted path is {");
        for (int i = 0; i < nn; ++i) {
            printf("%d", seq_indices[i]);
            if (i != nn - 1) printf(", ");
        }
        printf("} -> ");
    }

    while (!all_zero_indices) {    // until all indices are 0
        for (int i = 0; i < nn; ++i) {
            bool is_gap = gap_len[i] > 0;
            gap_len[i] -= 1;        // prepared for updating gap_len
            if (print_path) {
                if (is_gap == false) 
                    --seq_indices[i];
            } else
                msa_result[i][res_idx] = is_gap ? GAP : sequences[i][--seq_indices[i]];
        }
        ++res_idx;

        if (print_path) {
            printf("{");
            for (int i = 0; i < nn; ++i) {
                printf("%d", seq_indices[i]);
                if (i != nn - 1) printf(", ");
            }
            printf("} -> ");
        }

        // check the loop condition
        all_zero_indices = true;
        for (int i = 0; i < nn; ++i)
            if (seq_indices[i] != 0) {
                all_zero_indices = false;
                break;
            }
        // update gap_len by traversing the node in closed list
        STYPE max_gscore = cost_instead_of_score ? INT_MAX : INT_MIN;
        GapLen max_gap_len;
        for (auto iter = closed_list_gap_len[seq_indices].begin(); iter != closed_list_gap_len[seq_indices].end(); ++iter) {
            // parent's gap_len[i] = gap_len[i] - 1 (already minus 1 above)
            // if (gap_len[i] - 1 == -1), parent's gap_len[i] can be arbitrary number (no constraint)

            bool is_qualified = true;
            for (int seq_idx = 0; seq_idx < nn; ++seq_idx) {
                if (gap_len[seq_idx] == -1)
                    continue;
                else if (gap_len[seq_idx] != iter->first[seq_idx]) {
                    is_qualified = false;
                    break;
                }
            }
            if (is_qualified == true) {
                if (cost_instead_of_score && (iter->second.first < max_gscore) ||
                    !cost_instead_of_score && (iter->second.first > max_gscore)) {
                    max_gscore = iter->second.first;
                    max_gap_len = iter->first;
                }
            } 
        }
        gap_len = max_gap_len;
    }

    if (print_path == false) {
        // reverse the results
        for (int i = 0; i < res_idx / 2; ++i) {
            for (int seq = 0; seq < nn; ++seq) {
                char tmp = msa_result[seq][i];
                msa_result[seq][i] = msa_result[seq][res_idx-1-i];
                msa_result[seq][res_idx-1-i] = tmp;
            }
        }
    } else {
        printf("End! Path length = %d\n", res_idx);
    }

    return res_idx;
}

template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
void AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::check_memory_bound_trigger() {
    if (memory_bound_trigger == false) {
        bool actual_memory_usage = true;
        if (check_memory_thres(actual_memory_usage, get_closed_list_size())) {
            // triggered memory bound, set open_size_MB and closed_size_MB
            memory_bound_trigger = true;
            open_size_MB = workload_recorder.cur_open_list_cnt;
            closed_size_MB = get_closed_list_size();
            size_t total_size_MB = open_size_MB * open_node_size + closed_size_MB * closed_node_size;
            // at least one open node (the source) so: total_size_MB - open_node_size
            max_closed_size = (total_size_MB - open_node_size + closed_node_size - 1) / closed_node_size - 1;
            printf("Trigger memory bound, open_size_MB = %ld, closed_size_MB = %ld, max_closed_size = %ld\n", open_size_MB, closed_size_MB, max_closed_size);
        }
    }
}

// set res[0] to -1 when cannot find a qualified parent
template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
typename AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::GapLen AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::find_parent_in_closed_list(const NodeCoord &last_node_crd, const GapLen &last_node_tb_info, STYPE last_node_gscore) {
    STYPE max_gscore = cost_instead_of_score ? INT_MAX : INT_MIN;
    GapLen max_gap_len;
    NodeCoord parent_crd = last_node_crd;
    for (int i = 0; i < NN; ++i)
        if (is_gap(last_node_tb_info, i) == false) parent_crd[i] -= 1;

    STYPE ext_pen = blosum_table.get_ext_penalty(), open_pen = blosum_table.get_open_penalty();

    auto parent_closed_list_iter = closed_list_gap_len.find(parent_crd);
    if (parent_closed_list_iter != closed_list_gap_len.end()) {
        // traverse all directions
        // Condition 1: parent.gap_len + delta_gap_len = cur.gap_len
        // Condition 2: parent.gscore + delta_gscore >= cur.gscore
        for (auto inner_iter = parent_closed_list_iter->second.begin(); inner_iter != parent_closed_list_iter->second.end(); ++inner_iter) {
            // parent's gap_len[i] = gap_len[i] - 1
            // if (gap_len[i] == 0), parent's gap_len[i] can be arbitrary number (no constraint)
            bool is_gap_len_matched = true;
            for (int seq_idx = 0; seq_idx < NN; ++seq_idx) {
                if (last_node_tb_info[seq_idx] == 0)
                    continue;
                else if (last_node_tb_info[seq_idx] - 1 != inner_iter->first[seq_idx]) {
                    is_gap_len_matched = false;
                    break;
                }
            }

            if (is_gap_len_matched == true) {
                STYPE parent_gscore = inner_iter->second.first;
                STYPE updated_parent_gscore = parent_gscore;
                // check all pairs of sequences
                for (int seq_i = 0; seq_i < NN - 1; ++seq_i) {
                    for (int seq_j = seq_i + 1; seq_j < NN; ++seq_j) {
                        // compute and update the score of match/mismatch/gap
                        if (is_gap(last_node_tb_info, seq_i) == false && is_gap(last_node_tb_info, seq_j) == false) {     // check the score table
                            int symbol_i = sequences[seq_i][last_node_crd[seq_i] - 1], symbol_j = sequences[seq_j][last_node_crd[seq_j] - 1];
                            updated_parent_gscore += blosum_table.get_score_char(symbol_i, symbol_j);
                        } else {        // gap(s) in the two sequences
                            // There could be gaps in both sequences, so we need to check the gap length
                            bool gap_in_prev_seq_i = false, gap_in_prev_seq_j = false;
                            
                            int gap_len_seq_i = last_node_tb_info[seq_i], gap_len_seq_j = last_node_tb_info[seq_j];
                            if (gap_len_seq_i > gap_len_seq_j) gap_in_prev_seq_i = true;
                            if (gap_len_seq_j > gap_len_seq_i) gap_in_prev_seq_j = true;
                            
                            if (is_gap(last_node_tb_info, seq_i) && is_gap(last_node_tb_info, seq_j)) {     // two gaps here
                                // do nothing
                            } else if (is_gap(last_node_tb_info, seq_i)) {    // gap in seq i
                                updated_parent_gscore -= (!gap_in_prev_seq_i || last_node_crd[seq_i] == 0 && last_node_crd[seq_j] == 0) * open_pen + ext_pen;
                            } else {                                    // gap in seq j
                                updated_parent_gscore -= (!gap_in_prev_seq_j || last_node_crd[seq_i] == 0 && last_node_crd[seq_j] == 0) * open_pen + ext_pen;
                            }
                        }
                    }
                }

                if (cost_instead_of_score && (updated_parent_gscore > last_node_gscore) || 
                    !cost_instead_of_score && (updated_parent_gscore < last_node_gscore)) {
                        continue;       // condition 2 not fulfilled, skip this direction
                }

                if (cost_instead_of_score && (updated_parent_gscore < max_gscore) ||
                    !cost_instead_of_score && (updated_parent_gscore > max_gscore)) {
                    max_gscore = updated_parent_gscore;
                    max_gap_len = inner_iter->first;
                }
            } 
        }
    } else {
        max_gap_len[0] = -1;
    }
    return max_gap_len;
}


template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
void AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::print_progress_per_1M_iter() {
    if (workload_recorder.iter_cnt % 1000000 == 0) {
        printf("Total iter = %ld M, memory usage = %g bytes, cur_open_list_cnt = %ld, cur_closed_list_cnt = %ld, highest_fscore = %f\n", 
            workload_recorder.iter_cnt / 1000000, (float) get_memory_usage(), workload_recorder.cur_open_list_cnt, get_closed_list_size(), get_highest_fscore());
    }
}

template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
void AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::print_workload_report_after_search(int dim) {
    if (dim < NN) {
        if (get_closed_list_size() < max_closed_size) 
            printf("%dD recursive MSA completed!\n", dim);
        else
            printf("%dD recursive MSA failed! Closed list size > max_closed_size = %ld\n", dim, max_closed_size);
    } else {
        if (memory_bound_trigger) {
            printf("Iter %ld: First time memory bound is reached, open_list_MB = %ld, closed_list_MB = %ld, max_closed_list = %ld\n", 
                memory_bound_first_iter, open_size_MB, closed_size_MB, max_closed_size);
        } else {
            printf("Didn't reach the memory bound!\n");
        }

        if (get_closed_list_size() >= max_closed_size) 
            printf("Maximal closed list size reached! Insufficient memory space\n");
    }
    printf("Open list count = %ld, cur open list count = %ld, closed list count = %ld, cur closed list count = %ld, iteration count = %ld, nodes pruned by gscore = %ld, nodes pruned by hscore = %ld\n", 
        workload_recorder.open_list_cnt, workload_recorder.cur_open_list_cnt, workload_recorder.closed_list_cnt, get_closed_list_size(), workload_recorder.iter_cnt, workload_recorder.node_pruned_cnt_gscore, workload_recorder.node_pruned_cnt_hscore);
    printf("Gap open penalty = %g, Gap extension penalty = %g, Beam search iter = %d, Astar search iter = %d, Bin count = %d, Insertions of each bin = [", blosum_table.get_open_penalty(), blosum_table.get_ext_penalty(), beam_width, astar_iter_cnt, bin_cnt);
    for (int i = 0; i < bin_cnt; ++i) {
        if (i == bin_cnt - 1) printf("%g", bin_op_cnt[i].first);
        else printf("%g, ", bin_op_cnt[i].first);
    }
    printf("]\n");
    printf("Erasions (excluding the pruned nodes) of each bin = [");
    for (int i = 0; i < bin_cnt; ++i) {
        if (i == bin_cnt - 1) printf("%g", bin_op_cnt[i].second);
        else printf("%g, ", bin_op_cnt[i].second);
    }
    printf("]\n");
    // print the current memory consumption
    printf("Memory usage = %g bytes, memory threshold = %g bytes\n", (float) get_memory_usage(), memory_limit_ratio * physical_RAM_size);

}

template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
void AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::write_file_log() {
    std::string line = "Recursive MSA completed.\n";
    file_writer->write(line);

    // if (memory_bound_trigger) {
    //     printf("Iter %ld: First time memory bound is reached, open_list_MB = %ld, closed_list_MB = %ld, max_closed_list = %ld\n", 
    //         memory_bound_first_iter, open_size_MB, closed_size_MB, max_closed_size);
    // } else {
    //     printf("Didn't reach the memory bound!\n");
    // }

    if (memory_bound_trigger) {
        line = "Iter "; line += std::to_string(memory_bound_first_iter); line += ": First time memory bound is reached, open_list_MB = ";
        line += std::to_string(open_size_MB) + ", closed_list_MB = " + std::to_string(closed_size_MB) + ", max_closed_size = " + std::to_string(max_closed_size) + "\n";
        file_writer->write(line);
    } else {
        line = "Didn't reach the memory bound!\n";
        file_writer->write(line);
    }

    if (get_closed_list_size() >= max_closed_size) {
        line = "Maximal closed list size reached! Insufficient memory space\n";
        file_writer->write(line);
    }

    // printf("Open list count = %g, cur open list count = %g, closed list count = %g, cur closed list count = %ld, iteration count = %g, nodes pruned by gscore = %g, nodes pruned by hscore = %g\n", 
    //     workload_recorder.open_list_cnt, workload_recorder.cur_open_list_cnt, workload_recorder.closed_list_cnt, get_closed_list_size(), workload_recorder.iter_cnt, workload_recorder.node_pruned_cnt_gscore, workload_recorder.node_pruned_cnt_hscore);
    line = "Total open list insertions = " + std::to_string(workload_recorder.open_list_cnt);
    line += ", cur open list count = " + std::to_string(workload_recorder.cur_open_list_cnt);
    line += ", total closed list insertions = " + std::to_string(workload_recorder.closed_list_cnt);
    line += ", cur closed list count = " + std::to_string(get_closed_list_size());
    line += ", iteration count = " + std::to_string(workload_recorder.iter_cnt);
    line += ", nodes pruned by gscore = " + std::to_string(workload_recorder.node_pruned_cnt_gscore);
    line += ", nodes pruned by hscore = " + std::to_string(workload_recorder.node_pruned_cnt_hscore) + "\n";
    file_writer->write(line);

    // printf("Gap open penalty = %g, Gap extension penalty = %g, Beam search iter = %d, Astar search iter = %d, Bin count = %d, Insertions of each bin = [", blosum_table.get_open_penalty(), blosum_table.get_ext_penalty(), beam_width, astar_iter_cnt, bin_cnt);
    line = "Gap open penalty = " + std::to_string(blosum_table.get_open_penalty());
    line = "Gap extension penalty = " + std::to_string(blosum_table.get_ext_penalty());
    line = "Beam search iter = " + std::to_string(beam_width);
    line += ", Astar search iter = " + std::to_string(astar_iter_cnt);
    line += ", Bin count = " + std::to_string(bin_cnt) + "\n";
    file_writer->write(line);

    line = "Memory usage = " + std::to_string(get_memory_usage()) + " bytes, memory threshold = " + std::to_string(memory_limit_ratio * physical_RAM_size) + " bytes\n";
    file_writer->write(line);


    // printf("All PSA time = %f s, greedy MSA time = %f s\n", workload_recorder.all_PSA_time, workload_recorder.greedy_MSA_time);
    // printf("Recursive MSA time = [");
    // for (int i = 3; i <= nn; ++i) {
    //     if (i == nn) printf("%dD: %g s", i, workload_recorder.recursive_MSA_time[i]);
    //     else printf("%dD: %g s, ", i, workload_recorder.recursive_MSA_time[i]);
    // }
    // printf("]\n");

    line = "All PSA time = " + std::to_string(workload_recorder.all_PSA_time) + " s, greedy MSA time = " + std::to_string(workload_recorder.greedy_MSA_time) + " s\n";
    file_writer->write(line);

    line = "Recursive MSA time = [";
    for (int i = 3; i <= nn; ++i) {
        line += std::to_string(i) + "D: " + std::to_string(workload_recorder.recursive_MSA_time[i]) + " s";
        if (i != nn) line += ", ";
    }
    line += "]\n";
    file_writer->write(line);


    // total time
    float total_time = time_limit_timer->get_time_from_start();
    line = "Total time: " + std::to_string(total_time) + " s\n";
    file_writer->write(line);
}

#endif  // ANYTIMEASTAR_TPP