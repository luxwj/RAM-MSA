#ifndef ANYTIMEASTARLINEARGAP_HPP
#define ANYTIMEASTARLINEARGAP_HPP

#include "AnytimeAstar.hpp"

#include <ctime>

template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
class AnytimeAstarLinearGapSolver : public AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>{
private:
    bool is_gap(const DIRTYPE &tb_info, int seq_idx) {return ((tb_info >> seq_idx) & 0x01) == 0;};

public:
    using typename AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::NodeCoord;
    using typename AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::WorkloadRecorder;

    struct ClosedListNodeLinearGap {
        STYPE gscore = INT32_MAX;
        STYPE reexp_fscore = 0;        // valid when enable_memory_bound = true, init to INT32_MIN when using scores and INT32_MAX when using costs
        DIRTYPE came_from;
    };

    // VectorIntHash in utils.hpp
    std::unordered_map<NodeCoord, ClosedListNodeLinearGap, VectorIntHash<NN>> closed_list_linear_gap;      // key, value = index, node_info


    size_t get_closed_list_size() override {return closed_list_linear_gap.size();};

    /* --- recursive MSA --- */
    std::unordered_map<NodeCoord, STYPE, VectorIntHash<NN>> recursive_closed_list_linear_gap;
    STYPE score_upper_bound_recursive_linear_gap(int dim, bool access_reversed_mats, const NodeCoord &cur_crd);
    void compute_recursive_Astar_linear_gap(int seq_cnt, bool reverse_seq);
    /* --- recursive MSA --- */

    // prune the unnecessary nodes in the closed list when reaching the memory bound and return the prunable count
    int prune_closed_list_linear_gap(bool access_reversed_mats, int dim);

    // return: if nbr_crd is duplicated
    bool closed_list_duplication_detection_linear_gap(bool &use_reexp_fscore_instead, const NodeCoord &nbr_crd, DIRTYPE dir, STYPE cur_gscore, STYPE &cur_fscore);

    void open_list_insertion_linear_gap(STYPE cur_gscore, STYPE cur_fscore, std::vector<OLCoord*>* target_index_by_coord, const NodeCoord &nbr_crd, DIRTYPE nbr_tb_info);
    

    // If the current alignment result is optimal, return true. Otherwise, return false
    // reverse_seq: read the symbol in each sequences from the end instead of the beginning to generate the recursive closed list
    // specific_bin = false: find the node with the largest fscore in all bins 
    bool expand_node_linear_gap(int dim, bool specific_bin, int bin_level, bool reverse_seq);
    void compute_anytime_Astar_linear_gap();
    int backtrack_linear_gap(char **&msa_result, NodeCoord cur_crd, DIRTYPE cur_path, bool print_path);        // return the alignment length
    
    void AnytimeAstar_MSA_linear_gap(char **&msa_result, int &alignment_len, bool &reached_time_limit_output);       // return the access count and the alignment length

    AnytimeAstarLinearGapSolver(bool _cost_instead_of_score, int sequence_count, int sequence_max_length, char **original_sequences, int *sequence_lengths, const ScoreTable &_blosum_table, STYPE mafft_score, 
        int _bin_cnt, int _beam_width, int _astar_iter_cnt, float _time_limit, float _memory_limit_ratio, std::string anytime_result_dir, bool _enable_recursive_MSA, bool _sum_of_crd_as_level, std::string input_dir) : 
        AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::AnytimeAstarSolver(_cost_instead_of_score, sequence_count, sequence_max_length, original_sequences, sequence_lengths, _blosum_table, mafft_score, 
        _bin_cnt, _beam_width, _astar_iter_cnt, _time_limit, _memory_limit_ratio, anytime_result_dir, _enable_recursive_MSA, _sum_of_crd_as_level, input_dir) {

    }
};

template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
inline void AnytimeAstarLinearGapSolver<OLMultiIdx, OLFscore, OLCoord, NN>::compute_recursive_Astar_linear_gap(int seq_cnt, bool reverse_seq) {
    // temporarily change the hyper-parameters: only use Astar in recursive MSA
    int prev_beam_width = this->beam_width, prev_bin_cnt = this->bin_cnt;
    this->beam_width = 0; this->bin_cnt = 1;

    /* --- Open list of nn-D MSA --- */
    // initialize open_list_arr, deallocation seems unnecessary
    this->open_list_arr.clear();
    this->index_by_fscore_arr.clear();
    this->index_by_coord_arr.clear();

    this->open_list_arr = std::vector<OLMultiIdx>(this->bin_cnt);
    for (int i = 0; i < this->bin_cnt; ++i) {
        this->index_by_fscore_arr.push_back(&this->open_list_arr[i].template get<typename AstarMultiIndex<NN>::IndexByFscore>());
        this->index_by_coord_arr.push_back(&this->open_list_arr[i].template get<typename AstarMultiIndex<NN>::IndexByCoord>());
    }

    const int dim = seq_cnt;
    NodeCoord cur_crd;
    for (int i = 0; i < NN; ++i) cur_crd[i] = 0;

    DIRTYPE cur_path = (1<<dim)-1;     // 1 in each direction

    // Move the previous closed list results into recursive_closed_list_linear_gap
    recursive_closed_list_linear_gap.clear();
    for (auto closed_iter = closed_list_linear_gap.begin(); closed_iter != closed_list_linear_gap.end(); ++closed_iter) {
        recursive_closed_list_linear_gap[closed_iter->first] = closed_iter->second.gscore;
    }
    closed_list_linear_gap.clear();
    
    // if reverse_seq == false, access reversed_mats
    bool access_reversed_mats = !reverse_seq;
    STYPE highest_score = score_upper_bound_recursive_linear_gap(dim, access_reversed_mats, cur_crd);

    this->bin_score_thres = (highest_score - this->bin_score_offset) / this->bin_cnt;

    printf("In compute recursive Astar, dim = %d, highest score (recursive) = %f, bin_score_offset = %g, bin_score_thres = %g\n", dim, highest_score, this->bin_score_offset, this->bin_score_thres);


    typename AstarMultiIndex<NN>::OpenListNodeLinearGap source_node = {
        cur_crd,
        cur_path,
        0 + highest_score,
        0
    };
    this->index_by_fscore_arr[0]->insert(source_node);          // bin[0] because gscore = 0
    this->workload_recorder.open_list_cnt += 1;
    this->workload_recorder.cur_open_list_cnt = 1;              // set to 1 because we cleared the open list

    this->bin_op_cnt[0].first = 1;

    bool is_opt = false;        // if the current alignment result is optimal
    while (is_opt == false && this->time_limit_timer->check_time_limit() == false && get_closed_list_size() < this->max_closed_size && !this->open_list_arr_empty()) {
        for (int astar_iter = 0; astar_iter < this->astar_iter_cnt; ++astar_iter) {       // perform normal A* Search iterations
            if (is_opt || this->time_limit_timer->check_time_limit() || this->open_list_arr_empty()) break;
            this->print_progress_per_1M_iter();
            this->workload_recorder.iter_cnt += 1;
            bool specific_bin = false;      // expand a node in bin_arr by searching one with the largest fscore
            is_opt = expand_node_linear_gap(dim, specific_bin, 0, reverse_seq);
        }
    }

    this->print_workload_report_after_search(dim);

    // restore the hyper-parameters
    this->beam_width = prev_beam_width; this->bin_cnt = prev_bin_cnt;
}


template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
inline int AnytimeAstarLinearGapSolver<OLMultiIdx, OLFscore, OLCoord, NN>::prune_closed_list_linear_gap(bool access_reversed_mats, int dim) {
    printf("Start searching prunable closed list node, closed list size = %ld\n", closed_list_linear_gap.size());
    int prunable_cnt = 0;
    for (auto closed_iter = closed_list_linear_gap.begin(); closed_iter != closed_list_linear_gap.end();) {
        STYPE closed_gscore = closed_iter->second.gscore;
        STYPE closed_hscore = score_upper_bound_recursive_linear_gap(dim, access_reversed_mats, closed_iter->first);

        if (this->cost_instead_of_score && (closed_gscore + closed_hscore > this->tmp_msa_gscore) ||
            !this->cost_instead_of_score && (closed_gscore + closed_hscore < this->tmp_msa_gscore)) {
            closed_iter = closed_list_linear_gap.erase(closed_iter);
            prunable_cnt += 1;
        } else {
            ++closed_iter;
        }
    }
    printf("Total prunable count = %d, closed list size = %ld\n", prunable_cnt, closed_list_linear_gap.size());
    return prunable_cnt;
}

// return: if nbr_crd is duplicated. Erase the node from closed list if the node is worse than nbr_node
template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
inline bool AnytimeAstarLinearGapSolver<OLMultiIdx, OLFscore, OLCoord, NN>::closed_list_duplication_detection_linear_gap(bool &use_reexp_fscore_instead, const NodeCoord &nbr_crd, DIRTYPE dir, STYPE cur_gscore, STYPE &cur_fscore) {
    auto closed_list_iter = closed_list_linear_gap.find(nbr_crd);
    if (closed_list_iter != closed_list_linear_gap.end()) {
        // prune this node after checking the reexp_fscore
        if (this->cost_instead_of_score && (cur_gscore > closed_list_iter->second.gscore) ||
            !this->cost_instead_of_score && (cur_gscore < closed_list_iter->second.gscore)) {
            this->workload_recorder.node_pruned_cnt_gscore += 1;
            return true;
        } else if (cur_gscore == closed_list_iter->second.gscore) {
            // memory bound
            // if reexp-fscore is not the initial value, use its reexp-fscore instead if it's found in the Closed List
            if (closed_list_iter->second.reexp_fscore != this->init_reexp_fscore) {
                cur_fscore = closed_list_iter->second.reexp_fscore;
                use_reexp_fscore_instead = true;
            } else {
                // not a reexp node
                this->workload_recorder.node_pruned_cnt_gscore += 1;
                return true;
            }
        } else {
            /* --- Cannot erase nodes because some suboptimal alignment found before the search terminates may require the nodes --- */
            // Solution: In backtracking, find the node in open list?

            // erase the node in closed list
            closed_list_linear_gap.erase(closed_list_iter);
        }
    }
    return false;
}

template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
inline void AnytimeAstarLinearGapSolver<OLMultiIdx, OLFscore, OLCoord, NN>::open_list_insertion_linear_gap(STYPE cur_gscore, STYPE cur_fscore, std::vector<OLCoord*>* target_index_by_coord,
    const NodeCoord &nbr_crd, DIRTYPE nbr_tb_info) {

    if (this->sum_of_crd_as_level) {
        int nbr_bin = 0;
        for (int i = 0; i < NN; ++i) nbr_bin += nbr_crd[i];
        if (nbr_bin >= this->bin_cnt) nbr_bin = this->bin_cnt - 1;      // mainly for the recursive MSA because bin_cnt is set to 1 in lower dimensions 
        // nbr_crd cannot be in other bins because the sum of coordiantes is the same
        auto iterFound = (*target_index_by_coord)[nbr_bin]->find(nbr_crd);
        if (iterFound != (*target_index_by_coord)[nbr_bin]->end()) {
            // Updated! Prev: compare gscores and overwrite gscore; Curr: compare fscores and overwrite all node infos
            if (this->cost_instead_of_score && (cur_fscore >= iterFound->fscore) || 
                !this->cost_instead_of_score && (cur_fscore <= iterFound->fscore)) {
                this->workload_recorder.node_pruned_cnt_gscore += 1;      // prune this node
            } else {    // update gscore, fscore, and tbinfo
                (*target_index_by_coord)[nbr_bin]->modify(iterFound, typename AstarMultiIndex<NN>::ChangeGscoreFscoreTBinfoLinearGap(cur_gscore, cur_fscore, nbr_tb_info));
                this->workload_recorder.node_pruned_cnt_gscore += 1;      // prune this node
            }
        } else {    // insert a new node
            typename AstarMultiIndex<NN>::OpenListNodeLinearGap nbr_node = {
                nbr_crd,
                nbr_tb_info,
                cur_fscore,
                cur_gscore
            };
            (*target_index_by_coord)[nbr_bin]->insert(nbr_node);
            this->workload_recorder.open_list_cnt += 1;
            this->workload_recorder.cur_open_list_cnt += 1;
            this->bin_op_cnt[nbr_bin].first += 1;     // insertion + 1
        }
    } else {
        // Detect duplicated nodes in OpenList
        int nbr_bin = (int)floor((cur_gscore - this->bin_score_offset) / this->bin_score_thres);

        if (nbr_bin < 0) nbr_bin = 0;
        if (nbr_bin >= this->bin_cnt) nbr_bin = this->bin_cnt - 1;

        bool found_in_higher_bin = false;    // First find the nbr_node in higher bins
        if (this->cost_instead_of_score) {
            for (int tmp_bin_idx = nbr_bin - 1; tmp_bin_idx >= 0; --tmp_bin_idx) {
                auto iterFound = (*target_index_by_coord)[tmp_bin_idx]->find(nbr_crd);
                if (iterFound != (*target_index_by_coord)[tmp_bin_idx]->end()) {     // cur_gscore cannot be larger than the previous one
                    this->workload_recorder.node_pruned_cnt_gscore += 1;
                    found_in_higher_bin = true;
                    break;
                }
            }
        } else {
            for (int tmp_bin_idx = nbr_bin + 1; tmp_bin_idx < this->bin_cnt; ++tmp_bin_idx) {
                auto iterFound = (*target_index_by_coord)[tmp_bin_idx]->find(nbr_crd);
                if (iterFound != (*target_index_by_coord)[tmp_bin_idx]->end()) {     // cur_gscore cannot be larger than the previous one
                    this->workload_recorder.node_pruned_cnt_gscore += 1;
                    found_in_higher_bin = true;
                    break;
                }
            }
        }
        if (found_in_higher_bin == false) {  // If the node is not in bin[nbr_bin + 1], search bin[nbr_bin]
            auto iterFound = (*target_index_by_coord)[nbr_bin]->find(nbr_crd);
            if (iterFound != (*target_index_by_coord)[nbr_bin]->end()) {
                if (this->cost_instead_of_score && (cur_fscore >= iterFound->fscore) || 
                    !this->cost_instead_of_score && (cur_fscore <= iterFound->fscore)) {      // if cur_gscore <= gscore, skip; else update gscore
                    this->workload_recorder.node_pruned_cnt_gscore += 1;      // prune this node
                } else {    // update gscore, fscore, and tbinfo
                    (*target_index_by_coord)[nbr_bin]->modify(iterFound, typename AstarMultiIndex<NN>::ChangeGscoreFscoreTBinfoLinearGap(cur_gscore, cur_fscore, nbr_tb_info));
                    this->workload_recorder.node_pruned_cnt_gscore += 1;      // prune this node
                }
            } else {    // insert a new node
                typename AstarMultiIndex<NN>::OpenListNodeLinearGap nbr_node = {
                    nbr_crd,
                    nbr_tb_info,
                    cur_fscore,
                    cur_gscore
                };
                (*target_index_by_coord)[nbr_bin]->insert(nbr_node);
                this->workload_recorder.open_list_cnt += 1;
                this->workload_recorder.cur_open_list_cnt += 1;
                this->bin_op_cnt[nbr_bin].first += 1;     // insertion + 1
            }
        }
        // prune this crd in lower bin
        if (this->cost_instead_of_score) {
            if (nbr_bin + 1 < this->bin_cnt) {
                auto iterFound = (*target_index_by_coord)[nbr_bin + 1]->find(nbr_crd);
                if (iterFound != (*target_index_by_coord)[nbr_bin + 1]->end()) {
                    this->workload_recorder.node_pruned_cnt_gscore += 1;      // prune this node
                    (*target_index_by_coord)[nbr_bin + 1]->erase(iterFound);
                    this->workload_recorder.cur_open_list_cnt -= 1;
                    this->bin_op_cnt[nbr_bin + 1].second += 1;
                }
            }
        } else {
            if (nbr_bin > 0) {
                auto iterFound = (*target_index_by_coord)[nbr_bin - 1]->find(nbr_crd);
                if (iterFound != (*target_index_by_coord)[nbr_bin - 1]->end()) {
                    this->workload_recorder.node_pruned_cnt_gscore += 1;      // prune this node
                    (*target_index_by_coord)[nbr_bin - 1]->erase(iterFound);
                    this->workload_recorder.cur_open_list_cnt -= 1;
                    this->bin_op_cnt[nbr_bin - 1].second += 1;
                }
            }
        }
    }
}



// If the current alignment result is optimal, return true. Otherwise, return false
// reverse_seq: read the symbol in each sequences from the end instead of the beginning to generate the recursive closed list
// specific_bin = false: find the node with the largest fscore in all bins
template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
inline bool AnytimeAstarLinearGapSolver<OLMultiIdx, OLFscore, OLCoord, NN>::expand_node_linear_gap(int dim, bool specific_bin, int bin_level, bool reverse_seq) {
    STYPE ext_pen = this->blosum_table.get_ext_penalty(), open_pen = this->blosum_table.get_open_penalty();

    std::vector<OLFscore*>* target_index_by_fscore = &this->index_by_fscore_arr;
    std::vector<OLCoord*>* target_index_by_coord = &this->index_by_coord_arr;

    bool access_reversed_mats = !reverse_seq;

    // search the node with the highest fscore and update the bin_level
    if (this->bin_cnt > 1 && specific_bin == false) {
        // update the bin_level with the bin_idx containing the highest fscore
        this->get_highest_fscore(bin_level);
    }

    auto cur_node = (*target_index_by_fscore)[bin_level]->begin();


    NodeCoord cur_crd = cur_node->crd;
    DIRTYPE cur_tb_info = cur_node->tb_info;
    STYPE local_gscore = cur_node->gscore;
    STYPE local_fscore = cur_node->fscore;

    (*target_index_by_fscore)[bin_level]->erase(cur_node);
    this->workload_recorder.cur_open_list_cnt -= 1;
    this->bin_op_cnt[bin_level].second += 1;

    bool target_reached = this->is_crd_out_of_bound(cur_crd, dim);

    // find the cur_node in ClosedList

    // first find the coord
    auto closed_list_iter = closed_list_linear_gap.find(cur_crd);
    bool write_into_closed_list = true;
    // if we found the node, update the gscore and direction
    if (closed_list_iter != closed_list_linear_gap.end()) {
        if (closed_list_iter->second.reexp_fscore == this->init_reexp_fscore) {
            // Case 1: Not reexpanding the cur node, return when cur_gscore <= closed_gscore
            if (this->cost_instead_of_score && (local_gscore >= closed_list_iter->second.gscore) || 
                !this->cost_instead_of_score && (local_gscore <= closed_list_iter->second.gscore)) {
                return false;
            }
        } else {
            // Case 2: Reexpanding the cur node, return when cur_gscore < closed_gscore
            if (this->cost_instead_of_score && (local_gscore > closed_list_iter->second.gscore) || 
                !this->cost_instead_of_score && (local_gscore < closed_list_iter->second.gscore)) {
                return false;
            }
        }
    }

    // if we didn't find the node or cur_node is better, insert it into the closed list
    if (write_into_closed_list) {
        closed_list_linear_gap[cur_crd] = {local_gscore, this->init_reexp_fscore, cur_tb_info};
        this->workload_recorder.closed_list_cnt += 1;
    }

    if (target_reached) {
        if (this->cost_instead_of_score && (local_gscore < this->tmp_msa_gscore) ||
            !this->cost_instead_of_score && (local_gscore > this->tmp_msa_gscore)) {
            this->tmp_msa_gscore = local_gscore;    // update the score lower bound with the beam search result 
            // don't call backtrack when dim < nn, i.e., during recursive MSA with a lower dimension
            if (dim == this->nn) {
                printf("Before backtracking\n");
                this->tmp_msa_alignment_length = backtrack_linear_gap(this->tmp_msa_result, cur_crd, cur_tb_info, false);
                this->file_writer->write(this->time_limit_timer->get_time_from_start(), this->tmp_msa_result, this->tmp_msa_alignment_length, local_gscore, this->workload_recorder.iter_cnt);
                printf("Found a better MSA (score = %f)!\n", local_gscore);

                // write the TC score
                if (this->blosum_table.is_ref_MSA_set()) {
                    STYPE tc_score = this->blosum_table.compute_total_column_score(this->tmp_msa_result, this->nn, this->tmp_msa_alignment_length);
                    STYPE sp_accuracy = this->blosum_table.compute_sum_of_pairs_accuracy(this->tmp_msa_result, this->nn, this->tmp_msa_alignment_length);
                    printf("TC score = %f, SP accuracy = %f\n", tc_score, sp_accuracy);
                    this->file_writer->write_TC_score(tc_score);
                    this->file_writer->write_SP_accuracy(sp_accuracy);
                }

                // find the prunable nodes in the closed list
                prune_closed_list_linear_gap(access_reversed_mats, dim);
            }
        }

        bool is_opt_ret = true;
        STYPE global_best_result = this->get_highest_fscore();
        if (this->cost_instead_of_score && (local_gscore > global_best_result) || 
            !this->cost_instead_of_score && (local_gscore < global_best_result)) is_opt_ret = false;

        return is_opt_ret;
    }
    
    DIRTYPE dir_cnt = (1 << dim);       // direction counts

    // pre-load the match / mismatch scores of each sequence pair
    STYPE mm_scores[NN * NN];
    this->get_match_mismatch_scores(mm_scores, cur_crd, dim, reverse_seq);

    for (DIRTYPE dir = 1; dir < dir_cnt; ++dir) {
        NodeCoord nbr_crd (cur_crd); // first set nbr_crd = cur_crd, then check if we should "+= 1" to each dim

        bool oob_flag = false;  // node out-of-bound flag
        for (int bit = 0; bit < dim; ++bit) {
            if ((dir>>bit) & 0x01) {
                nbr_crd[bit] += 1;                        // if i == 0b011, then dir_vec = [1, 1, 0]
                if (nbr_crd[bit] > this->seq_lens[bit]) oob_flag = true;      // don't go beyond the last cell
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
                if (nbr_crd[seq_i] != cur_crd[seq_i] && nbr_crd[seq_j] != cur_crd[seq_j]) {     // check the score table
                    // cur_gscore += this->blosum_table.get_score_char(symbols[seq_i], symbols[seq_j]);
                    cur_gscore += mm_scores[seq_i * NN + seq_j];
                } else {        // gap(s) in the two sequences
                    if (nbr_crd[seq_i] == cur_crd[seq_i] && nbr_crd[seq_j] == cur_crd[seq_j]) {     // two gaps here
                        // do nothing
                    } else {    // gap in seq i or seq j
                        cur_gscore -= ext_pen;
                    }
                }
            }
        }
        bool use_reexp_fscore_instead = false;

        bool is_nbr_duplicated = closed_list_duplication_detection_linear_gap(use_reexp_fscore_instead, nbr_crd, dir, cur_gscore, cur_fscore);
        if (is_nbr_duplicated) continue;

        if (use_reexp_fscore_instead == false) {
            cur_fscore += cur_gscore;
            // compute the hscore separately
            STYPE cur_hscore = score_upper_bound_recursive_linear_gap(dim, access_reversed_mats, nbr_crd);
            cur_fscore += cur_hscore;
        }
        
        // if cur_fscore < tmp_msa_gscore, don't insert it into the open_list and closed_list
        if (this->cost_instead_of_score && (cur_fscore > this->tmp_msa_gscore + __FLT_EPSILON__) || 
            !this->cost_instead_of_score && (cur_fscore + __FLT_EPSILON__ < this->tmp_msa_gscore)) {
            this->workload_recorder.node_pruned_cnt_hscore += 1;
            continue;
        }

        open_list_insertion_linear_gap(cur_gscore, cur_fscore, target_index_by_coord, nbr_crd, dir);

    }
    
    bool is_opt_ret = true;
    STYPE global_best_result = this->get_highest_fscore();
    if (this->cost_instead_of_score && (this->tmp_msa_gscore > global_best_result) ||
        !this->cost_instead_of_score && (this->tmp_msa_gscore < global_best_result)) is_opt_ret = false;

    // memory bound
    // check once per 100 iterations
    if (this->workload_recorder.iter_cnt % 100 == 0 && is_opt_ret == false && this->enable_memory_bound == true) {

        this->check_memory_bound_trigger();

        while (this->memory_bound_trigger && this->check_memory_thres(false, closed_list_linear_gap.size())) {
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
            while (non_empty_bin_idx < this->bin_cnt && (*target_index_by_fscore)[non_empty_bin_idx]->empty())
                non_empty_bin_idx += 1;
            if (non_empty_bin_idx == this->bin_cnt)     // open list is empty
                return is_opt_ret;

            last_node_iter = prev((*target_index_by_fscore)[non_empty_bin_idx]->end());

            NodeCoord last_node_crd = last_node_iter->get_crd();

            if (last_node_crd == src_crd) {
                // if source is the only node, go to the next bin
                if ((*target_index_by_fscore)[non_empty_bin_idx]->size() == 1) {
                    non_empty_bin_idx += 1;
                    while (non_empty_bin_idx < this->bin_cnt && (*target_index_by_fscore)[non_empty_bin_idx]->empty())
                        non_empty_bin_idx += 1;
                    // access the last node in non_empty_bin_idx
                    if (non_empty_bin_idx == this->bin_cnt) {
                        return is_opt_ret;
                    }
                    last_node_iter = prev((*target_index_by_fscore)[non_empty_bin_idx]->end());
                } else 
                    last_node_iter = prev(last_node_iter);

                last_node_crd = last_node_iter->get_crd();
            }

            // Calculate the parent's crd and erase the last_node
            STYPE last_node_fscore = last_node_iter->fscore;
            DIRTYPE last_node_tb_info = last_node_iter->get_tb_info();            // a single DIRTYPE variable instead of a vector
            (*target_index_by_fscore)[non_empty_bin_idx]->erase(last_node_iter);
            this->workload_recorder.cur_open_list_cnt -= 1;
            this->bin_op_cnt[non_empty_bin_idx].second += 1;     // erasion + 1

            NodeCoord parent_crd = last_node_crd;
            for (int i = 0; i < this->nn; ++i)
                if (is_gap(last_node_tb_info, i) == false) parent_crd[i] -= 1;

            // Find the parent in Closed List and update the reexp_fscore
            auto parent_closed_list_iter = closed_list_linear_gap.find(parent_crd);
            if (parent_closed_list_iter != closed_list_linear_gap.end()) {
                DIRTYPE parent_tb_info = parent_closed_list_iter->second.came_from;
                if (this->cost_instead_of_score && (last_node_fscore < parent_closed_list_iter->second.reexp_fscore) || 
                    !this->cost_instead_of_score && (last_node_fscore > parent_closed_list_iter->second.reexp_fscore)) {
                    // printf("Updated parent in closed list\n");

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
                    parent_closed_list_iter->second.reexp_fscore = last_node_fscore;
                }

                // Find the parent in Open List and update the fscore (only search the lower bins)
                // Collect the tb_info in Closed List first. There is only one direction in linear cases
                bool found_parent_in_open_list = false;
                // score: parent_node <= cur_node, step = -1
                // cost: parent_node <= cur_node, step = -1
                for (int tmp_bin_idx = non_empty_bin_idx; tmp_bin_idx >= 0; --tmp_bin_idx) {
                    auto parent_open_list_iter = (*target_index_by_coord)[tmp_bin_idx]->find(parent_crd);
                    if (parent_open_list_iter != (*target_index_by_coord)[tmp_bin_idx]->end()) {
                        if (this->cost_instead_of_score && (last_node_fscore < parent_open_list_iter->fscore) || 
                            !this->cost_instead_of_score && (last_node_fscore > parent_open_list_iter->fscore)) {
                            (*target_index_by_coord)[tmp_bin_idx]->modify(parent_open_list_iter, typename AstarMultiIndex<NN>::ChangeFscoreLinearGap(last_node_fscore));
                        }
                        
                        STYPE closed_list_gscore = parent_closed_list_iter->second.gscore;
                        // compare the gscores
                        if (this->cost_instead_of_score && (closed_list_gscore >= parent_open_list_iter->gscore) || 
                            !this->cost_instead_of_score && (closed_list_gscore <= parent_open_list_iter->gscore)) {      // if closed_list_gscore <= gscore, skip; else update gscore
                            // do nothing
                        } else {    // update gscore
                            (*target_index_by_coord)[tmp_bin_idx]->modify(parent_open_list_iter, typename AstarMultiIndex<NN>::ChangeGscoreLinearGap(closed_list_gscore));
                        }

                        found_parent_in_open_list = true;
                        break;
                    }

                    // the parent can only appear in [non_empty_bin_idx - dim, non_empty_bin_idx - 1]
                    if (this->sum_of_crd_as_level && tmp_bin_idx == non_empty_bin_idx - dim)
                        break;
                }

                // Insert the parent into the Open List
                if (found_parent_in_open_list == false) {
                    typename AstarMultiIndex<NN>::OpenListNodeLinearGap reexp_parent_node = {
                        parent_crd, 
                        parent_tb_info,
                        parent_closed_list_iter->second.reexp_fscore,
                        parent_closed_list_iter->second.gscore
                    };

                    // prune the parent node with tmp_msa_gscore
                    if (this->cost_instead_of_score && (reexp_parent_node.fscore > this->tmp_msa_gscore + __FLT_EPSILON__) || 
                        !this->cost_instead_of_score && (reexp_parent_node.fscore + __FLT_EPSILON__ < this->tmp_msa_gscore)) {
                        this->workload_recorder.node_pruned_cnt_hscore += 1;
                        continue;
                    }
                    
                    int parent_bin;
                    if (this->sum_of_crd_as_level) {
                        parent_bin = 0;
                        for (int i = 0; i < NN; ++i) parent_bin += parent_crd[i];
                    } else {
                        parent_bin = (int)floor((reexp_parent_node.gscore - this->bin_score_offset) / this->bin_score_thres);
                        if (parent_bin < 0) parent_bin = 0;
                        if (parent_bin >= this->bin_cnt) parent_bin = this->bin_cnt - 1;
                    }

                    (*target_index_by_fscore)[parent_bin]->insert(reexp_parent_node);
                    this->workload_recorder.open_list_cnt += 1;
                    this->workload_recorder.cur_open_list_cnt += 1;
                    this->bin_op_cnt[parent_bin].first += 1;     // insertion + 1
                }
            }
        }
    }

    return is_opt_ret;
}

template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
inline void AnytimeAstarLinearGapSolver<OLMultiIdx, OLFscore, OLCoord, NN>::compute_anytime_Astar_linear_gap() {
    const int dim = this->nn;
    NodeCoord cur_crd;
    for (int i = 0; i < NN; ++i) cur_crd[i] = 0;

    DIRTYPE cur_path = (1<<dim)-1;     // 1 in each direction
    
    // Move the previous closed list results into recursive_closed_list_linear_gap
    recursive_closed_list_linear_gap.clear();
    for (auto closed_iter = closed_list_linear_gap.begin(); closed_iter != closed_list_linear_gap.end(); ++closed_iter) {
        recursive_closed_list_linear_gap[closed_iter->first] = closed_iter->second.gscore;
    }
    closed_list_linear_gap.clear();

    printf("Initialization in compute_anytime_Astar_linear_gap(), closed list size = %d (should be 0)\n", (int) closed_list_linear_gap.size());

    CRDTYPE pairwise_heuristic_crd[NN];
    for (int i = 0; i < dim; ++i) 
        pairwise_heuristic_crd[i] = this->seq_lens[i] - cur_crd[i];
    bool access_reversed_mats = true;
    STYPE highest_score = this->pairwise_heuristic->get_sum_of_scores_linear_gap(dim, access_reversed_mats, pairwise_heuristic_crd);

    this->bin_score_thres = (highest_score - this->bin_score_offset) / this->bin_cnt;

    printf("highest score (recursive) = %f, bin_score_offset = %g, bin_score_thres = %g\n", highest_score, this->bin_score_offset, this->bin_score_thres);

    typename AstarMultiIndex<NN>::OpenListNodeLinearGap source_node = {
        cur_crd, 
        cur_path,
        0 + highest_score,
        0
    };
    this->index_by_fscore_arr[0]->insert(source_node);         // bin[0] because gscore = 0
    this->workload_recorder.open_list_cnt += 1;
    this->workload_recorder.cur_open_list_cnt = 1;              // set to 1 because we cleared the open list

    this->bin_op_cnt[0].first = 1;
    
    bool is_opt = false;        // if the current alignment result is optimal
    while (is_opt == false && this->time_limit_timer->check_time_limit() == false && !this->open_list_arr_empty()) {
        for (int bin_level = 0; bin_level < this->bin_cnt; ++bin_level) {     // perform Anytime Column Search iterations
            int beam = 0;
            while (is_opt == false && this->time_limit_timer->check_time_limit() == false && !this->open_list_arr_empty(bin_level) && beam < this->beam_width) {
                this->print_progress_per_1M_iter();
                this->workload_recorder.iter_cnt += 1;
                bool specific_bin = true;   // expand a node in bin[bin_level]
                bool reverse_seq = false;   // we are not computing recursive MSA
                is_opt = expand_node_linear_gap(dim, specific_bin, bin_level, reverse_seq);
                beam += 1;
            }
        }
        for (int astar_iter = 0; astar_iter < this->astar_iter_cnt; ++astar_iter) {       // perform normal A* Search iterations
            if (is_opt || this->time_limit_timer->check_time_limit() || this->open_list_arr_empty()) break;
            this->print_progress_per_1M_iter();
            this->workload_recorder.iter_cnt += 1;
            bool specific_bin = false;      // expand a node in bin_arr by searching one with the largest fscore
            bool reverse_seq = false;       // we are not computing recursive MSA
            is_opt = expand_node_linear_gap(dim, specific_bin, 0, reverse_seq);
        }
        // print the iter count when memory bound is reached for the first time
        // printf("iter = %g, trigger = %d, first iter = %d\n", this->workload_recorder.iter_cnt, (int)this->memory_bound_trigger, this->memory_bound_first_iter);
        if (this->memory_bound_trigger == true && this->memory_bound_first_iter == 0) 
            this->memory_bound_first_iter = this->workload_recorder.iter_cnt;

        // prune nodes in the closed list 
        if (closed_list_linear_gap.size() >= this->max_closed_size) {
            bool access_reversed_mats = true;
            printf("The closed list consumed all available memory!\n");
            int prunable_cnt = prune_closed_list_linear_gap(access_reversed_mats, dim);
            if (prunable_cnt == 0) break;
        }
    }
    this->print_workload_report_after_search(NN);
}


template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
inline void AnytimeAstarLinearGapSolver<OLMultiIdx, OLFscore, OLCoord, NN>::AnytimeAstar_MSA_linear_gap(char **&msa_result, int &alignment_len, bool &reached_time_limit_output) {
    this->time_limit_timer->start(this->time_limit);

    this->preprocessing_MSA_with_PSA_and_greedy();

    Timer MSA_timer;

    if (this->enable_recursive_MSA) {
        // compute recursive MSA from 3D to (n-1)-D
        STYPE prev_tmp_msa_gscore = this->tmp_msa_gscore;
        for (int recur_idx = 3; recur_idx < this->nn; ++recur_idx) {
            // compute recur_idx dimensional MSA
            bool reverse_seq = (this->nn - recur_idx) & 0x01;         // reverse the sequences when computing (nn-1)-D MSA, (nn-3)-D MSA, ...
            
            this->tmp_msa_gscore = this->cost_instead_of_score ? INT_MAX : INT_MIN;
            this->memory_bound_trigger = false;

            MSA_timer.start();
            compute_recursive_Astar_linear_gap(recur_idx, reverse_seq);
            this->workload_recorder.recursive_MSA_time[recur_idx] = MSA_timer.elapsed(false);
        }
        this->tmp_msa_gscore = prev_tmp_msa_gscore;
    }

    /* --- Open list of nn-D MSA --- */
    // initialize open_list_arr, deallocation seems unnecessary

    this->open_list_arr.clear();
    this->index_by_fscore_arr.clear();
    this->index_by_coord_arr.clear();

    this->open_list_arr = std::vector<OLMultiIdx>(this->bin_cnt);
    for (int i = 0; i < this->bin_cnt; ++i) {
        this->index_by_fscore_arr.push_back(&this->open_list_arr[i].template get<typename AstarMultiIndex<NN>::IndexByFscore>());
        this->index_by_coord_arr.push_back(&this->open_list_arr[i].template get<typename AstarMultiIndex<NN>::IndexByCoord>());
    }

    this->memory_bound_trigger = false;
    
    MSA_timer.start();
    compute_anytime_Astar_linear_gap();
    this->workload_recorder.recursive_MSA_time[this->nn] = MSA_timer.elapsed(false);

    this->print_computation_time_after_search();

    // init the result with gaps, assuming the size of the msa_result is [nn][nn*ll]
    for (int row = 0; row < this->nn; ++row)
        for (int col = 0; col < this->nn * this->ll; ++col) {
            msa_result[row][col] = GAP;
        }

    // backtrack is done whenever the target node is reached in expand_node_linear_gap(), just copy the results from tmp_msa_result to msa_result here
    for (int row = 0; row < this->nn; ++row) {
        for (int col = 0; col < this->tmp_msa_alignment_length; ++col) {
            msa_result[row][col] = this->tmp_msa_result[row][col];
        }
    }
    alignment_len = this->tmp_msa_alignment_length;
    this->reached_time_limit = this->time_limit_timer->check_time_limit();
    reached_time_limit_output = this->reached_time_limit;

    this->write_file_log();

    // deallocate the closed list
    closed_list_linear_gap.clear();

    return;
}

template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
inline STYPE AnytimeAstarLinearGapSolver<OLMultiIdx, OLFscore, OLCoord, NN>::score_upper_bound_recursive_linear_gap(int dim, bool access_reversed_mats, const NodeCoord &cur_crd) {
    STYPE res = 0;        // we first compute the recursive score (without the last sequence), then add the PSA scores of the last sequence

    if (this->enable_recursive_MSA) {
        // first find the crd in the recursive closed list
        // use the first dim - 1 elements because recursive_closed_list_linear_gap contains (dim - 1)-D results
        NodeCoord reversed_crd;
        for (int i = 0; i < NN; ++i) reversed_crd[i] = 0;
        for (int i = 0; i < dim - 1; ++i) reversed_crd[i] = this->seq_lens[i] - cur_crd[i];

        auto recur_closed_list_iter = recursive_closed_list_linear_gap.find(reversed_crd);
        // compare the (dim-1)*(dim-1) distance matrix
        if (recur_closed_list_iter != recursive_closed_list_linear_gap.end()) {
            STYPE recur_gscore = recur_closed_list_iter->second;
            res = recur_gscore;

            STYPE hscore_last_seq = 0;
            // add the PSA scores of the last sequence (deterministic)
            for (int seq_i = 0; seq_i < dim - 1; ++seq_i) {
                int seq_j = dim - 1;
                hscore_last_seq += this->pairwise_heuristic->get_score(access_reversed_mats, seq_i, seq_j, cur_crd[seq_i], cur_crd[seq_j]);
            }
            res += hscore_last_seq;
            return res;
        }
    }

    CRDTYPE pairwise_heuristic_crd[NN];
    for (int i = 0; i < dim; ++i) 
        pairwise_heuristic_crd[i] = this->seq_lens[i] - cur_crd[i];
    res = this->pairwise_heuristic->get_sum_of_scores_linear_gap(dim, access_reversed_mats, pairwise_heuristic_crd);

    return res;
}


// write the result into msa_result and return alignment length
// print_path == true: only print the backtrack path in terminal; print_path == false: output the alignment
template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
inline int AnytimeAstarLinearGapSolver<OLMultiIdx, OLFscore, OLCoord, NN>::backtrack_linear_gap(char **&msa_result, NodeCoord cur_crd, DIRTYPE cur_path, bool print_path) {
    NodeCoord seq_indices = cur_crd;        // indices of the TB_mat, use --idx{i} when accessing sequences[][]

    int res_idx = 0;                                // counter of the alignment length
    bool all_zero_indices = true;
    for (int i = 0; i < this->nn; ++i)
        if (seq_indices[i] != 0) {
            all_zero_indices = false;
            break;
        }
    
    DIRTYPE parent_dir = cur_path;           // a single DIRTYPE variable instead of a vector

    if (print_path) {
        printf("Inverted path is {");
        for (int i = 0; i < this->nn; ++i) {
            printf("%d", seq_indices[i]);
            if (i != this->nn - 1) printf(", ");
        }
        printf("} -> ");
    }

    while (!all_zero_indices) {    // until all indices are 0
        for (int i = 0; i < this->nn; ++i) {
            bool is_gap = ((parent_dir >> i) & 0x01) == 0;
            if (print_path) {
                if (is_gap == false) 
                    --seq_indices[i];
            } else
                msa_result[i][res_idx] = is_gap ? GAP : this->sequences[i][--seq_indices[i]];
        }
        ++res_idx;

        if (print_path) {
            printf("{");
            for (int i = 0; i < this->nn; ++i) {
                printf("%d", seq_indices[i]);
                if (i != this->nn - 1) printf(", ");
            }
            printf("} -> ");
        }

        // check the loop condition
        all_zero_indices = true;
        for (int i = 0; i < this->nn; ++i)
            if (seq_indices[i] != 0) {
                all_zero_indices = false;
                break;
            }

        // update parent_dir
        auto parent_in_closed_iter = closed_list_linear_gap.find(seq_indices);
        if (parent_in_closed_iter == closed_list_linear_gap.end()) {
            // find the node in the open list
            for (int tmp_bin_idx = this->bin_cnt - 1; tmp_bin_idx >= 0; --tmp_bin_idx) {
                auto parent_in_open_iter = this->index_by_coord_arr[tmp_bin_idx]->find(seq_indices);
                if (parent_in_open_iter != this->index_by_coord_arr[tmp_bin_idx]->end()) {
                    parent_dir = parent_in_open_iter->tb_info;
                    break;
                }
            }
        } else {
            parent_dir = parent_in_closed_iter->second.came_from;
        }

        // naive solution: access the node directly. May fail when erasing nodes during search
        // parent_dir = closed_list_linear_gap[seq_indices].came_from;
    }

    if (print_path == false) {
        // reverse the results
        for (int i = 0; i < res_idx / 2; ++i) {
            for (int seq = 0; seq < this->nn; ++seq) {
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

#endif // ANYTIMEASTARLINEARGAP_HPP
