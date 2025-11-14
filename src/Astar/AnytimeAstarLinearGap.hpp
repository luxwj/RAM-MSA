#ifndef ANYTIMEASTARLINEARGAP_HPP
#define ANYTIMEASTARLINEARGAP_HPP

#include "AnytimeAstar.hpp"

template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
class AnytimeAstarLinearGapSolver : public AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>{
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

    /* --- recursive MSA --- */
    std::unordered_map<NodeCoord, STYPE, VectorIntHash<NN>> recursive_closed_list_linear_gap;
    STYPE score_upper_bound_recursive_linear_gap(int dim, bool access_reversed_mats, NodeCoord cur_crd);
    void compute_recursive_Astar_linear_gap(int seq_cnt, bool reverse_seq);
    /* --- recursive MSA --- */

    STYPE score_upper_bound_2D_linear_gap(bool access_reversed_mats, int seq_idx[], int seq_offset[], int cur_seq_lens[]);     // return the score_upper_bound
    
    // initialize the tmp_msa_gscore with DFS
    void compute_greedy_linear_gap(int dim, bool reverse_seq, bool ordered_by_fscore_instead_of_gscore);

    // prune the unnecessary nodes in the closed list when reaching the memory bound and return the prunable count
    int prune_closed_list(bool access_reversed_mats, int dim);

    // If the current alignment result is optimal, return true. Otherwise, return false
    // reverse_seq: read the symbol in each sequences from the end instead of the beginning to generate the recursive closed list
    // specific_bin = false: find the node with the largest fscore in all bins 
    bool expand_node_linear_gap(int dim, bool specific_bin, int bin_level, bool reverse_seq, WorkloadRecorder& workload_recorder, std::vector<std::pair<double, double>>& bin_op_cnt);
    void compute_anytime_Astar_linear_gap();
    // codeupdate
    int backtrack_linear_gap(char **&msa_result, typename AstarMultiIndex<NN>::template CoordTBInfo<LINEAR_TB_TYPE> crd_tb_info, bool print_path);        // return the alignment length
    
    double AnytimeAstar_MSA_linear_gap(char **&msa_result, int &alignment_len, bool &reached_time_limit_output);       // return the access count and the alignment length

    void write_file_log_linear_gap();

    AnytimeAstarLinearGapSolver(bool _cost_instead_of_score, int sequence_count, int sequence_max_length, char **original_sequences, int *sequence_lengths, ScoreTable _blosum_table, STYPE mafft_score, 
        int _bin_cnt, int _beam_width, int _astar_iter_cnt, double _time_limit, double _memory_limit_ratio, std::string anytime_result_dir, bool _enable_recursive_MSA, std::string input_dir) : 
        AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::AnytimeAstarSolver(_cost_instead_of_score, sequence_count, sequence_max_length, original_sequences, sequence_lengths, _blosum_table, mafft_score, 
        _bin_cnt, _beam_width, _astar_iter_cnt, _time_limit, _memory_limit_ratio, anytime_result_dir, _enable_recursive_MSA, input_dir) {
    }
};

template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
inline void AnytimeAstarLinearGapSolver<OLMultiIdx, OLFscore, OLCoord, NN>::compute_recursive_Astar_linear_gap(int seq_cnt, bool reverse_seq) {
    // temporarily change the hyper-parameters: only use Astar in recursive MSA
    int prev_beam_width = this->beam_width, prev_bin_cnt = this->bin_cnt;
    this->beam_width = 0; this->bin_cnt = 1;         // todo: should we use beam search in lower dimension?

    /* --- Open list of nn-D MSA --- */
    // initialize open_list_arr, deallocation seems unnecessary
    this->open_list_arr.clear();
    this->index_by_fscore_arr.clear();        // todo: do we need to deallocate the pointers?
    this->index_by_coord_arr.clear();

    this->open_list_arr = std::vector<OLMultiIdx>(this->bin_cnt);
    for (int i = 0; i < this->bin_cnt; ++i) {
        this->index_by_fscore_arr.push_back(&this->open_list_arr[i].template get<typename AstarMultiIndex<NN>::IndexByFscore>());
        this->index_by_coord_arr.push_back(&this->open_list_arr[i].template get<typename AstarMultiIndex<NN>::IndexByCoord>());
    }

    // we don't need the alignment but the closed list nodes
    const int dim = seq_cnt;
    NodeCoord cur_crd;
    for (int i = 0; i < NN; ++i) cur_crd[i] = 0;

    DIRTYPE cur_path = (1<<dim)-1;     // 1 in each direction

    // codeupdate
    typename AstarMultiIndex<NN>::template CoordTBInfo<LINEAR_TB_TYPE> cur_crd_tb_info = {cur_crd, cur_path};

    // Move the previous closed list results into recursive_closed_list_linear_gap
    recursive_closed_list_linear_gap.clear();
    for (auto closed_iter = closed_list_linear_gap.begin(); closed_iter != closed_list_linear_gap.end(); ++closed_iter) {
        recursive_closed_list_linear_gap[closed_iter->first] = closed_iter->second.gscore;
    }
    closed_list_linear_gap.clear();

    STYPE source_reexp_fscore = this->init_reexp_fscore;

    closed_list_linear_gap[cur_crd_tb_info.crd] = {0.0, source_reexp_fscore, (DIRTYPE)(1<<dim)-1};
    
    // if reverse_seq == false, access reversed_mats
    bool access_reversed_mats = !reverse_seq;
    STYPE highest_score = score_upper_bound_recursive_linear_gap(dim, access_reversed_mats, cur_crd);

    printf("In compute recursive Astar, dim = %d, highest score (recursive) = %f\n", dim, highest_score);

    this->bin_score_thres = highest_score / this->bin_cnt;

    // codeupdate
    typename AstarMultiIndex<NN>::template OpenListNode<LINEAR_TB_TYPE> source_node = {
        cur_crd_tb_info,
        0 + highest_score,
        0
    };
    this->index_by_fscore_arr[0]->insert(source_node);          // bin[0] because gscore = 0
    this->workload_recorder.open_list_cnt += 1;
    this->workload_recorder.cur_open_list_cnt = 1;              // set to 1 because we cleared the open list

    std::vector<std::pair<double, double>> bin_op_cnt(this->bin_cnt, {0, 0});         // count insertions and erasions of each bin
    bin_op_cnt[0].first = 1;

    bool is_opt = false;        // if the current alignment result is optimal
    while (is_opt == false && this->time_limit_timer->check_time_limit() == false && (int) closed_list_linear_gap.size() < this->max_closed_size && !this->open_list_arr_empty()) {
        for (int astar_iter = 0; astar_iter < this->astar_iter_cnt; ++astar_iter) {       // perform normal A* Search iterations
            if (is_opt || this->time_limit_timer->check_time_limit() || this->open_list_arr_empty()) break;
            if ((int)this->workload_recorder.iter_cnt % 1000000 == 0) {
                printf("Total iter = %g, memory usage = %g bytes, cur_open_list_cnt = %d, highest_fscore = %f\n", 
                    this->workload_recorder.iter_cnt, (double) get_memory_usage(), (int) this->workload_recorder.cur_open_list_cnt, this->get_highest_fscore());
            }
            this->workload_recorder.iter_cnt += 1;
            bool specific_bin = false;      // expand a node in bin_arr by searching one with the largest fscore
            is_opt = expand_node_linear_gap(dim, specific_bin, 0, reverse_seq, this->workload_recorder, bin_op_cnt);
        }
    }

    if ((int) closed_list_linear_gap.size() < this->max_closed_size) 
        printf("%dD recursive MSA completed!\n", dim);
    else
        printf("%dD recursive MSA failed! Closed list size > max_closed_size = %ld\n", dim, this->max_closed_size);
    
    printf("Total open list insertions = %g, total closed list insertions = %g, iteration count = %g, nodes pruned by gscore = %g, nodes pruned by hscore = %g\n", 
        this->workload_recorder.open_list_cnt, this->workload_recorder.closed_list_cnt, this->workload_recorder.iter_cnt, this->workload_recorder.node_pruned_cnt_gscore, this->workload_recorder.node_pruned_cnt_hscore);
    printf("Insertions of each bin = [");
    for (int i = 0; i < this->bin_cnt; ++i) {
        if (i == this->bin_cnt - 1) printf("%g", bin_op_cnt[i].first);
        else printf("%g, ", bin_op_cnt[i].first);
    }
    printf("]\n");
    printf("Erasions (excluding the pruned nodes) of each bin = [");
    for (int i = 0; i < this->bin_cnt; ++i) {
        if (i == this->bin_cnt - 1) printf("%g", bin_op_cnt[i].second);
        else printf("%g, ", bin_op_cnt[i].second);
    }
    printf("]\n");
    // print the current memory consumption
    printf("Memory usage = %g bytes, memory threshold = %g bytes\n", (double) get_memory_usage(), 0.8 * this->physical_RAM_size);

    // restore the hyper-parameters
    this->beam_width = prev_beam_width; this->bin_cnt = prev_bin_cnt;
}


// initialize the tmp_msa_gscore with DFS
template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
inline void AnytimeAstarLinearGapSolver<OLMultiIdx, OLFscore, OLCoord, NN>::compute_greedy_linear_gap(int dim, bool reverse_seq, bool ordered_by_fscore_instead_of_gscore) {
    STYPE ext_pen = this->blosum_table.get_ext_penalty(), open_pen = this->blosum_table.get_open_penalty();
    // only stores the coordinate and gscore, generate the MSA results (this->tmp_msa_result) during DFS

    // we only use the first dim elements so seems to be safe
    NodeCoord cur_node, next_node;
    for (int i = 0; i < NN; ++i) cur_node[i] = 0;
    for (int i = 0; i < NN; ++i) next_node[i] = 0;

    STYPE local_gscore = 0, next_gscore;
    int cur_alignment_length = 0;

    bool target_reached = false;
    while (target_reached == false) {
        // traverse over all children
        DIRTYPE best_child_dir;
        STYPE best_child_score = this->cost_instead_of_score ? INT_MAX : INT_MIN;

        DIRTYPE dir_cnt = (1 << dim);       // direction counts
        for (DIRTYPE dir = 1; dir < dir_cnt; ++dir) {
            // we only use the first dim elements so seems to be safe
            NodeCoord nbr_crd;
            for (int i = 0; i < NN; ++i) nbr_crd[i] = 0;

            bool oob_flag = false;  // node out-of-bound flag
            for (int bit = 0; bit < dim; ++bit) {
                if ((dir>>bit) & 0x01) {
                    nbr_crd[bit] = cur_node[bit] + 1;                        // if i == 0b011, then dir_vec = [1, 1, 0]
                    if (nbr_crd[bit] > this->seq_lens[bit]) oob_flag = true;      // don't go beyond the last cell
                } else {
                    nbr_crd[bit] = cur_node[bit];
                }
            }
            if (oob_flag) continue;

            int symbols[dim];        // symbols of nbr_crd
            for (int seq_idx = 0; seq_idx < dim; ++seq_idx) {
                if (nbr_crd[seq_idx] == 0)  symbols[seq_idx] = SPACE;       // there is a space in front of each sequence
                else {
                    if (reverse_seq) {  // read the symbol from the end of the sequence
                        symbols[seq_idx] = this->sequences[seq_idx][this->seq_lens[seq_idx] - nbr_crd[seq_idx]];     // -1 is necessary
                    } else {            // read the symbol normally
                        symbols[seq_idx] = this->sequences[seq_idx][nbr_crd[seq_idx] - 1];     // -1 is necessary
                    }
                }
            }

            STYPE cur_gscore = local_gscore;
            STYPE cur_fscore = 0;               // valid when ordered_by_fscore_instead_of_gscore = true

            // check all pairs of sequences
            for (int seq_i = 0; seq_i < dim - 1; ++seq_i) {
                for (int seq_j = seq_i + 1; seq_j < dim; ++seq_j) {
                    // compute and update the score of match/mismatch/gap
                    if (nbr_crd[seq_i] != cur_node[seq_i] && nbr_crd[seq_j] != cur_node[seq_j]) {     // check the score table
                        cur_gscore += this->blosum_table.get_score_char(symbols[seq_i], symbols[seq_j]);
                    } else {        // gap(s) in the two sequences
                        if (nbr_crd[seq_i] == cur_node[seq_i] && nbr_crd[seq_j] == cur_node[seq_j]) {     // two gaps here
                            // do nothing
                        } else {    // gap in seq i or seq j
                            cur_gscore -= ext_pen;
                        }
                    }
                }
            }

            if (ordered_by_fscore_instead_of_gscore == true) {
                bool access_reversed_mats = !reverse_seq;
                STYPE cur_hscore = 0.0;
                for (int seq_i = 0; seq_i < dim - 1; ++seq_i) {
                    for (int seq_j = seq_i + 1; seq_j < dim; ++seq_j) {
                        int seq_idx[2] = {seq_i, seq_j}, seq_offset[2] = {nbr_crd[seq_i], nbr_crd[seq_j]}, cur_seq_lens[2] = {this->seq_lens[seq_i] - nbr_crd[seq_i], this->seq_lens[seq_j] - nbr_crd[seq_j]};
                        cur_hscore += score_upper_bound_2D_linear_gap(access_reversed_mats, seq_idx, seq_offset, cur_seq_lens);
                    }
                }

                cur_fscore = cur_gscore + cur_hscore;
                if (this->cost_instead_of_score && (cur_fscore < best_child_score) || !this->cost_instead_of_score && (cur_fscore > best_child_score)) {
                    best_child_score = cur_fscore;
                    best_child_dir = dir;
                    next_gscore = cur_gscore;
                    next_node = nbr_crd;
                }
            } else {        // compare using gscore
                if (this->cost_instead_of_score && (cur_gscore < best_child_score) || !this->cost_instead_of_score && (cur_gscore > best_child_score)) {
                    best_child_score = cur_gscore;
                    best_child_dir = dir;
                    next_gscore = cur_gscore;
                    next_node = nbr_crd;
                }
            }
        }

        // write a new column into this->tmp_msa_result
        // todo: check the correctness of the score
        for (int i = 0; i < dim; ++i) {
            bool is_gap = (next_node[i] == cur_node[i]);
            this->tmp_msa_result[i][cur_alignment_length] = is_gap ? GAP : this->sequences[i][next_node[i] - 1];
        }
        cur_alignment_length += 1;
        cur_node = next_node;
        local_gscore = next_gscore;

        target_reached = true;
        for (int d = 0; d < dim; ++d)
            if (cur_node[d] < this->seq_lens[d]) {
                target_reached = false;
                break;
            }
    }

    if (this->cost_instead_of_score && (local_gscore < this->tmp_msa_gscore) ||
        !this->cost_instead_of_score && (local_gscore > this->tmp_msa_gscore)) {
        this->tmp_msa_gscore = local_gscore;    // update the score lower bound with the beam search result 

        this->tmp_msa_alignment_length = cur_alignment_length;
        this->file_writer->write(this->time_limit_timer->get_time_from_start(), this->tmp_msa_result, this->tmp_msa_alignment_length, local_gscore, 0);
    }
    printf("In greedy MSA, ordered_by_fscore_instead_of_gscore = %d, local_gscore = %f\n", (int) ordered_by_fscore_instead_of_gscore, local_gscore);

}

template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
inline int AnytimeAstarLinearGapSolver<OLMultiIdx, OLFscore, OLCoord, NN>::prune_closed_list(bool access_reversed_mats, int dim) {
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

// If the current alignment result is optimal, return true. Otherwise, return false
// reverse_seq: read the symbol in each sequences from the end instead of the beginning to generate the recursive closed list
// specific_bin = false: find the node with the largest fscore in all bins
template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
inline bool AnytimeAstarLinearGapSolver<OLMultiIdx, OLFscore, OLCoord, NN>::expand_node_linear_gap(int dim, bool specific_bin, int bin_level, bool reverse_seq, WorkloadRecorder& workload_recorder, std::vector<std::pair<double, double>>& bin_op_cnt) {
    STYPE ext_pen = this->blosum_table.get_ext_penalty(), open_pen = this->blosum_table.get_open_penalty();

    std::vector<OLFscore*>* target_index_by_fscore = &this->index_by_fscore_arr;
    std::vector<OLCoord*>* target_index_by_coord = &this->index_by_coord_arr;

    bool access_reversed_mats = !reverse_seq;

    // search the node with the highest fscore and update the bin_level
    if (specific_bin == false) {
        // update the bin_level with the bin_idx containing the highest fscore
        this->get_highest_fscore(bin_level);
    }

    // codeupdate
    typename AstarMultiIndex<NN>::template CoordTBInfo<LINEAR_TB_TYPE> cur_crd_tb_info;
    STYPE local_gscore;

    auto cur_node = (*target_index_by_fscore)[bin_level]->begin();
    cur_crd_tb_info.crd = cur_node->crd_tb_info.crd;
    cur_crd_tb_info.tb_info = cur_node->crd_tb_info.tb_info;
    local_gscore = cur_node->gscore;

    STYPE local_fscore = cur_node->fscore;

    (*target_index_by_fscore)[bin_level]->erase((*target_index_by_fscore)[bin_level]->begin());
    this->workload_recorder.cur_open_list_cnt -= 1;
    bin_op_cnt[bin_level].second += 1;

    bool target_reached = true;
    for (int d = 0; d < dim; ++d)
        if (cur_crd_tb_info.crd[d] < this->seq_lens[d]) {
            target_reached = false;
            break;
        }

    // find the cur_node in ClosedList

    // first find the coord
    auto closed_list_iter = closed_list_linear_gap.find(cur_crd_tb_info.crd);
    bool write_into_closed_list = true;
    // if we found the node, update the gscore and direction
    if (closed_list_iter != closed_list_linear_gap.end()) {
        if (this->cost_instead_of_score && (local_gscore > closed_list_iter->second.gscore) || 
            !this->cost_instead_of_score && (local_gscore < closed_list_iter->second.gscore)) {
            return false;
        }
    }

    // if we didn't find the node or cur_node is better, insert it into the closed list
    if (write_into_closed_list) {
        STYPE init_reexp_fscore = this->init_reexp_fscore;      // memory-bound method
        closed_list_linear_gap[cur_crd_tb_info.crd] = {local_gscore, init_reexp_fscore, cur_crd_tb_info.tb_info};
        this->workload_recorder.closed_list_cnt += 1;
    }

    // same as that in affine gap ver.
    if (target_reached) {
        if (this->cost_instead_of_score && (local_gscore < this->tmp_msa_gscore) ||
            !this->cost_instead_of_score && (local_gscore > this->tmp_msa_gscore)) {
            this->tmp_msa_gscore = local_gscore;    // update the score lower bound with the beam search result 
            // don't call backtrack when dim < nn, i.e., during recursive MSA with a lower dimension
            if (dim == this->nn) {
                this->tmp_msa_alignment_length = backtrack_linear_gap(this->tmp_msa_result, cur_crd_tb_info, false);
                this->file_writer->write(this->time_limit_timer->get_time_from_start(), this->tmp_msa_result, this->tmp_msa_alignment_length, local_gscore, this->workload_recorder.iter_cnt);

                // find the prunable nodes in the closed list
                printf("Found a better MSA (score = %f)!\n", local_gscore);
                prune_closed_list(access_reversed_mats, dim);
            }
        }

        bool is_opt_ret = true;
        STYPE global_best_result = this->get_highest_fscore();
        if (this->cost_instead_of_score && (local_gscore > global_best_result) || 
            !this->cost_instead_of_score && (local_gscore < global_best_result)) is_opt_ret = false;

        return is_opt_ret;
    }
    
    DIRTYPE dir_cnt = (1 << dim);       // direction counts
    int newly_inserted_node_cnt = 0;             // count how many nodes to remove when enable_memory_bound

    for (DIRTYPE dir = 1; dir < dir_cnt; ++dir) {
        NodeCoord nbr_crd;
        for (int i = 0; i < NN; ++i) nbr_crd[i] = 0;

        bool oob_flag = false;  // node out-of-bound flag
        for (int bit = 0; bit < dim; ++bit) {
            if ((dir>>bit) & 0x01) {
                nbr_crd[bit] = cur_crd_tb_info.crd[bit] + 1;                        // if i == 0b011, then dir_vec = [1, 1, 0]
                if (nbr_crd[bit] > this->seq_lens[bit]) oob_flag = true;      // don't go beyond the last cell
            } else {
                nbr_crd[bit] = cur_crd_tb_info.crd[bit];
            }
        }
        if (oob_flag) continue;
        int symbols[dim];        // symbols of nbr_crd
        for (int seq_idx = 0; seq_idx < dim; ++seq_idx) {
            if (nbr_crd[seq_idx] == 0)  symbols[seq_idx] = SPACE;       // there is a space in front of each sequence
            else {
                if (reverse_seq) {  // read the symbol from the end of the sequence
                    // crd = 0, symbol = seq_lens[seq_idx]
                    // crd = 1, symbol = seq_lens[seq_idx] - 1
                    symbols[seq_idx] = this->sequences[seq_idx][this->seq_lens[seq_idx] - nbr_crd[seq_idx]];     // -1 is necessary
                } else {            // read the symbol normally
                    // crd = 0, symbol = -1
                    // crd = 1, symbol = 0
                    symbols[seq_idx] = this->sequences[seq_idx][nbr_crd[seq_idx] - 1];     // -1 is necessary
                }
            }
        }

        // 1. f-score = cur_g-score + match/mismatch/gap + h-score
        // 2. if f-score >= tmp_msa_gscore, push into open_list
        STYPE cur_gscore = local_gscore;        // cur_gscore = pre_gscore + the score of match/mismatch/gap
        STYPE cur_fscore = 0;                   // cur_fscore = cur_gscore + cur_hscore

        // check all pairs of sequences
        for (int seq_i = 0; seq_i < dim - 1; ++seq_i) {
            for (int seq_j = seq_i + 1; seq_j < dim; ++seq_j) {
                // compute and update the score of match/mismatch/gap
                if (nbr_crd[seq_i] != cur_crd_tb_info.crd[seq_i] && nbr_crd[seq_j] != cur_crd_tb_info.crd[seq_j]) {     // check the score table
                    cur_gscore += this->blosum_table.get_score_char(symbols[seq_i], symbols[seq_j]);
                } else {        // gap(s) in the two sequences
                    if (nbr_crd[seq_i] == cur_crd_tb_info.crd[seq_i] && nbr_crd[seq_j] == cur_crd_tb_info.crd[seq_j]) {     // two gaps here
                        // do nothing
                    } else {    // gap in seq i or seq j
                        cur_gscore -= ext_pen;
                    }
                }
            }
        }
        // tb_info.size() = 1; tb_info = dir
        // codeupdate
        typename AstarMultiIndex<NN>::template CoordTBInfo<LINEAR_TB_TYPE> nbr_crd_tb_info = {nbr_crd, dir};

        bool use_reexp_fscore_instead = false;

        auto closed_list_iter = closed_list_linear_gap.find(nbr_crd);
        if (closed_list_iter != closed_list_linear_gap.end()) {
            // prune this node after checking the reexp_fscore
            if (this->cost_instead_of_score && (cur_gscore > closed_list_iter->second.gscore) ||
                !this->cost_instead_of_score && (cur_gscore < closed_list_iter->second.gscore)) {
                this->workload_recorder.node_pruned_cnt_gscore += 1;
                continue;
            } else if (cur_gscore == closed_list_iter->second.gscore) {
                // memory bound
                // if reexp-fscore is not the initial value, use its reexp-fscore instead if it's found in the Closed List
                if (this->enable_memory_bound && dir == closed_list_iter->second.came_from && closed_list_iter->second.reexp_fscore != this->init_reexp_fscore) {
                    cur_fscore = closed_list_iter->second.reexp_fscore;
                    use_reexp_fscore_instead = true;
                } else {
                    // not a reexp node
                    this->workload_recorder.node_pruned_cnt_gscore += 1;
                    continue;
                }
            }
        }

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

        
        // Detect duplicated nodes in OpenList
        int nbr_bin = (int)floor(cur_gscore / this->bin_score_thres);

        // printf("computing nbr_bin, score thres = %f, cur_gscore = %f, nbr_bin = %d\n", 
        //     this->bin_score_thres, cur_gscore, nbr_bin);

        if (nbr_bin < 0) nbr_bin = 0;
        if (nbr_bin >= this->bin_cnt) {
            nbr_bin = this->bin_cnt - 1;
        }

        bool found_in_higher_bin = false;    // First find the nbr_node in higher bins
        if (this->cost_instead_of_score) {
            for (int tmp_bin_idx = nbr_bin - 1; tmp_bin_idx >= 0; --tmp_bin_idx) {
                auto iterFound = (*target_index_by_coord)[tmp_bin_idx]->find(nbr_crd_tb_info);
                if (iterFound != (*target_index_by_coord)[tmp_bin_idx]->end()) {     // cur_gscore cannot be larger than the previous one
                    this->workload_recorder.node_pruned_cnt_gscore += 1;
                    found_in_higher_bin = true;
                    break;
                }
            }
        } else {
            for (int tmp_bin_idx = nbr_bin + 1; tmp_bin_idx < this->bin_cnt; ++tmp_bin_idx) {
                auto iterFound = (*target_index_by_coord)[tmp_bin_idx]->find(nbr_crd_tb_info);
                if (iterFound != (*target_index_by_coord)[tmp_bin_idx]->end()) {     // cur_gscore cannot be larger than the previous one
                    this->workload_recorder.node_pruned_cnt_gscore += 1;
                    found_in_higher_bin = true;
                    break;
                }
            }
        }
        if (found_in_higher_bin == false) {  // If the node is not in bin[nbr_bin + 1], search bin[nbr_bin]
            auto iterFound = (*target_index_by_coord)[nbr_bin]->find(nbr_crd_tb_info);
            if (iterFound != (*target_index_by_coord)[nbr_bin]->end()) {
                if (this->cost_instead_of_score && (cur_gscore >= iterFound->gscore) || 
                    !this->cost_instead_of_score && (cur_gscore <= iterFound->gscore)) {      // if cur_gscore <= gscore, skip; else update gscore
                    this->workload_recorder.node_pruned_cnt_gscore += 1;      // prune this node
                } else {    // update gscore
                    // codeupdate
                    (*target_index_by_coord)[nbr_bin]->modify(iterFound, typename AstarMultiIndex<NN>::template ChangeGscore<LINEAR_TB_TYPE>(cur_gscore));
                    this->workload_recorder.node_pruned_cnt_gscore += 1;      // prune this node
                }
            } else {    // insert a new node
                // codeupdate
                typename AstarMultiIndex<NN>::template OpenListNode<LINEAR_TB_TYPE> nbr_node = {
                    nbr_crd_tb_info,
                    cur_fscore,
                    cur_gscore
                };
                (*target_index_by_fscore)[nbr_bin]->insert(nbr_node);

                this->workload_recorder.open_list_cnt += 1;
                this->workload_recorder.cur_open_list_cnt += 1;
                bin_op_cnt[nbr_bin].first += 1;     // insertion + 1
                newly_inserted_node_cnt += 1;                // count how many nodes to remove when the memory usage exceeds the threshold
            }
        }
        // prune this crd in lower bin
        if (this->cost_instead_of_score) {
            if (nbr_bin + 1 < this->bin_cnt) {
                auto iterFound = (*target_index_by_coord)[nbr_bin + 1]->find(nbr_crd_tb_info);
                if (iterFound != (*target_index_by_coord)[nbr_bin + 1]->end()) {
                    this->workload_recorder.node_pruned_cnt_gscore += 1;      // prune this node
                    (*target_index_by_coord)[nbr_bin + 1]->erase(iterFound);
                    this->workload_recorder.cur_open_list_cnt -= 1;
                    bin_op_cnt[nbr_bin + 1].second += 1;
                }
            }
        } else {
            if (nbr_bin > 0) {
                auto iterFound = (*target_index_by_coord)[nbr_bin - 1]->find(nbr_crd_tb_info);
                if (iterFound != (*target_index_by_coord)[nbr_bin - 1]->end()) {
                    this->workload_recorder.node_pruned_cnt_gscore += 1;      // prune this node
                    (*target_index_by_coord)[nbr_bin - 1]->erase(iterFound);
                    this->workload_recorder.cur_open_list_cnt -= 1;
                    bin_op_cnt[nbr_bin - 1].second += 1;
                }
            }
        }
    }

    bool is_opt_ret = true;
    STYPE global_best_result = this->get_highest_fscore();
    if (this->cost_instead_of_score && (this->tmp_msa_gscore > global_best_result) ||
        !this->cost_instead_of_score && (this->tmp_msa_gscore < global_best_result)) is_opt_ret = false;

    // memory bound
    // check once per 100 iterations
    if ((int)this->workload_recorder.iter_cnt % 100 == 0 && is_opt_ret == false && this->enable_memory_bound == true) {

        if (this->memory_bound_trigger == false) {
            bool actual_memory_usage = true;
            if (this->check_memory_thres(actual_memory_usage, closed_list_linear_gap.size())) {
                // triggered memory bound, set open_size_MB and closed_size_MB
                this->memory_bound_trigger = true;
                this->open_size_MB = workload_recorder.cur_open_list_cnt;
                this->closed_size_MB = closed_list_linear_gap.size();
                size_t total_size_MB = this->open_size_MB * this->open_node_size + this->closed_size_MB * this->closed_node_size;
                // at least one open node (the source) so: total_size_MB - this->open_node_size
                this->max_closed_size = (total_size_MB - this->open_node_size + this->closed_node_size - 1) / this->closed_node_size - 1;
                printf("Trigger memory bound, open_size_MB = %ld, closed_size_MB = %ld, max_closed_size = %ld\n", this->open_size_MB, this->closed_size_MB, this->max_closed_size);
            }
        }

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

            if (last_node_iter->crd_tb_info.crd == src_crd) {

                // if source is the only node, go to the next bin
                if ((*target_index_by_fscore)[non_empty_bin_idx]->size() == 1) {
                    non_empty_bin_idx += 1;
                    while (non_empty_bin_idx < this->bin_cnt && (*target_index_by_fscore)[non_empty_bin_idx]->empty())
                        non_empty_bin_idx += 1;
                    // access the last node in non_empty_bin_idx
                    if (non_empty_bin_idx == this->bin_cnt) {
                        return is_opt_ret;
                        // int prunable_cnt = prune_closed_list(access_reversed_mats, dim);
                        // if (prunable_cnt == 0) return is_opt_ret;
                    }
                    last_node_iter = prev((*target_index_by_fscore)[non_empty_bin_idx]->end());
                } else 
                    last_node_iter = prev(last_node_iter);
            }

            // Calculate the parent's crd and erase the last_node
            NodeCoord last_node_crd = last_node_iter->crd_tb_info.crd;
            STYPE last_node_fscore = last_node_iter->fscore;
            DIRTYPE last_node_tb_info = last_node_iter->crd_tb_info.tb_info;            // a single DIRTYPE variable instead of a vector
            (*target_index_by_fscore)[non_empty_bin_idx]->erase(last_node_iter);
            this->workload_recorder.cur_open_list_cnt -= 1;
            bin_op_cnt[non_empty_bin_idx].second += 1;     // erasion + 1

            NodeCoord parent_crd = last_node_crd;
            for (int i = 0; i < this->nn; ++i) {
                bool is_gap = ((last_node_tb_info >> i) & 0x01) == 0;
                if (is_gap == false) parent_crd[i] -= 1;
            }

            // Find the parent in Closed List and update the reexp_fscore
            auto parent_closed_list_iter = closed_list_linear_gap.find(parent_crd);
            if (parent_closed_list_iter != closed_list_linear_gap.end()) {
                DIRTYPE parent_tb_info = parent_closed_list_iter->second.came_from;
                if (this->cost_instead_of_score && (last_node_fscore < parent_closed_list_iter->second.reexp_fscore) || 
                    !this->cost_instead_of_score && (last_node_fscore > parent_closed_list_iter->second.reexp_fscore)) {
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
                // The key is crd_tb_info, so collect the tb_info in Closed List first. There is only one direction in linear cases.
                bool found_parent_in_open_list = false;
                // codeupdate
                typename AstarMultiIndex<NN>::template CoordTBInfo<LINEAR_TB_TYPE> parent_crd_tb_info = {parent_crd, parent_tb_info};
                // score: parent_node <= cur_node, step = -1
                // cost: parent_node <= cur_node, step = -1
                for (int tmp_bin_idx = non_empty_bin_idx; tmp_bin_idx >= 0; --tmp_bin_idx) {
                    auto parent_open_list_iter = (*target_index_by_coord)[tmp_bin_idx]->find(parent_crd_tb_info);
                    if (parent_open_list_iter != (*target_index_by_coord)[tmp_bin_idx]->end()) {

                        if (this->cost_instead_of_score && (last_node_fscore < parent_open_list_iter->fscore) || 
                            !this->cost_instead_of_score && (last_node_fscore > parent_open_list_iter->fscore)) {
                            // codeupdate
                            (*target_index_by_coord)[tmp_bin_idx]->modify(parent_open_list_iter, typename AstarMultiIndex<NN>::template ChangeFscore<LINEAR_TB_TYPE>(last_node_fscore));
                        }
                        
                        STYPE closed_list_gscore = parent_closed_list_iter->second.gscore;
                        // compare the gscores
                        if (this->cost_instead_of_score && (closed_list_gscore >= parent_open_list_iter->gscore) || 
                            !this->cost_instead_of_score && (closed_list_gscore <= parent_open_list_iter->gscore)) {      // if closed_list_gscore <= gscore, skip; else update gscore
                            // do nothing
                        } else {    // update gscore
                            // codeupdate
                            (*target_index_by_coord)[tmp_bin_idx]->modify(parent_open_list_iter, typename AstarMultiIndex<NN>::template ChangeGscore<LINEAR_TB_TYPE>(closed_list_gscore));
                        }

                        found_parent_in_open_list = true;
                        break;
                    }
                }

                // Insert the parent into the Open List
                if (found_parent_in_open_list == false) {
                    // codeupdate
                    typename AstarMultiIndex<NN>::template OpenListNode<LINEAR_TB_TYPE> reexp_parent_node = {
                        parent_crd_tb_info,
                        parent_closed_list_iter->second.reexp_fscore,
                        parent_closed_list_iter->second.gscore
                    };

                    int parent_bin = (int)floor(reexp_parent_node.gscore / this->bin_score_thres);
                    if (parent_bin < 0) parent_bin = 0;
                    if (parent_bin >= this->bin_cnt) parent_bin = this->bin_cnt - 1;

                    // prune the parent node with tmp_msa_gscore
                    if (this->cost_instead_of_score && (reexp_parent_node.fscore > this->tmp_msa_gscore + __FLT_EPSILON__) || 
                        !this->cost_instead_of_score && (reexp_parent_node.fscore + __FLT_EPSILON__ < this->tmp_msa_gscore)) {
                        this->workload_recorder.node_pruned_cnt_hscore += 1;
                        continue;
                    }
        
                    (*target_index_by_fscore)[parent_bin]->insert(reexp_parent_node);
                    this->workload_recorder.open_list_cnt += 1;
                    this->workload_recorder.cur_open_list_cnt += 1;
                    bin_op_cnt[parent_bin].first += 1;     // insertion + 1
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
    // codeupdate
    typename AstarMultiIndex<NN>::template CoordTBInfo<LINEAR_TB_TYPE> cur_crd_tb_info = {cur_crd, cur_path};

    STYPE source_reexp_fscore = this->init_reexp_fscore;
    
    // Move the previous closed list results into recursive_closed_list_linear_gap
    recursive_closed_list_linear_gap.clear();
    for (auto closed_iter = closed_list_linear_gap.begin(); closed_iter != closed_list_linear_gap.end(); ++closed_iter) {
        recursive_closed_list_linear_gap[closed_iter->first] = closed_iter->second.gscore;
    }
    closed_list_linear_gap.clear();

    closed_list_linear_gap[cur_crd_tb_info.crd] = {0.0, source_reexp_fscore, (DIRTYPE)(1<<dim)-1};
    this->workload_recorder.closed_list_cnt += 1;       // total insertion including every dimension
    printf("Initialization in compute_anytime_Astar_linear_gap(), closed list size = %d (should be 1)\n", (int) closed_list_linear_gap.size());
    
    STYPE highest_score = 0;
    for (int seq_i = 0; seq_i < dim - 1; ++seq_i) {
        for (int seq_j = seq_i + 1; seq_j < dim; ++seq_j) {
            // compute PSA score with NW_2D
            int seq_idx[2] = {seq_i, seq_j}, seq_offset[2] = {0, 0}, cur_seq_lens[2] = {this->seq_lens[seq_i], this->seq_lens[seq_j]};
            bool access_reversed_mats = true;
            STYPE score_upper_bound = score_upper_bound_2D_linear_gap(access_reversed_mats, seq_idx, seq_offset, cur_seq_lens);
            highest_score += score_upper_bound;
        }
    }

    printf("h-score of the source node = %f\n", highest_score);

    this->bin_score_thres = highest_score / this->bin_cnt;
    // codeupdate
    typename AstarMultiIndex<NN>::template OpenListNode<LINEAR_TB_TYPE> source_node = {
        cur_crd_tb_info,
        0 + highest_score,
        0
    };
    this->index_by_fscore_arr[0]->insert(source_node);         // bin[0] because gscore = 0
    this->workload_recorder.open_list_cnt += 1;
    this->workload_recorder.cur_open_list_cnt = 1;              // set to 1 because we cleared the open list

    std::vector<std::pair<double, double>> bin_op_cnt(this->bin_cnt, {0, 0});         // count insertions and erasions of each bin
    bin_op_cnt[0].first = 1;
    
    bool is_opt = false;        // if the current alignment result is optimal
    while (is_opt == false && this->time_limit_timer->check_time_limit() == false && !this->open_list_arr_empty()) {
        for (int bin_level = 0; bin_level < this->bin_cnt; ++bin_level) {     // perform Anytime Column Search iterations
            int beam = 0;
            while (is_opt == false && this->time_limit_timer->check_time_limit() == false && !this->open_list_arr_empty(bin_level) && beam < this->beam_width) {
                if ((int)this->workload_recorder.iter_cnt % 1000000 == 0) {
                    printf("Total iter = %g, memory usage = %g bytes, cur_open_list_cnt = %d, cur_closed_list_cnt = %d, highest_fscore = %f\n", 
                        this->workload_recorder.iter_cnt, (double) get_memory_usage(), (int) this->workload_recorder.cur_open_list_cnt, (int) closed_list_linear_gap.size(), this->get_highest_fscore());
                }
                this->workload_recorder.iter_cnt += 1;
                bool specific_bin = true;   // expand a node in bin[bin_level]
                bool reverse_seq = false;   // we are not computing recursive MSA
                is_opt = expand_node_linear_gap(dim, specific_bin, bin_level, reverse_seq, this->workload_recorder, bin_op_cnt);
                beam += 1;
            }
        }
        for (int astar_iter = 0; astar_iter < this->astar_iter_cnt; ++astar_iter) {       // perform normal A* Search iterations
            if (is_opt || this->time_limit_timer->check_time_limit() || this->open_list_arr_empty()) break;
            if ((int)this->workload_recorder.iter_cnt % 1000000 == 0) {
                printf("Total iter = %g, memory usage = %g bytes, cur_open_list_cnt = %d, cur_closed_list_cnt = %d, highest_fscore = %f\n", 
                    this->workload_recorder.iter_cnt, (double) get_memory_usage(), (int) this->workload_recorder.cur_open_list_cnt, (int) closed_list_linear_gap.size(), this->get_highest_fscore());
            }
            this->workload_recorder.iter_cnt += 1;
            bool specific_bin = false;      // expand a node in bin_arr by searching one with the largest fscore
            bool reverse_seq = false;       // we are not computing recursive MSA
            is_opt = expand_node_linear_gap(dim, specific_bin, 0, reverse_seq, this->workload_recorder, bin_op_cnt);
        }
        // print the iter count when memory bound is reached for the first time
        // printf("iter = %g, trigger = %d, first iter = %d\n", this->workload_recorder.iter_cnt, (int)this->memory_bound_trigger, this->memory_bound_first_iter);
        if (this->memory_bound_trigger == true && this->memory_bound_first_iter == -1) 
            this->memory_bound_first_iter = (int)this->workload_recorder.iter_cnt;

        // prune nodes in the closed list 
        if (closed_list_linear_gap.size() >= this->max_closed_size) {
            bool access_reversed_mats = true;
            printf("The closed list consumed all available memory!\n");
            int prunable_cnt = prune_closed_list(access_reversed_mats, dim);
            if (prunable_cnt == 0) break;
        }
    }

    printf("Iter %d: First time memory bound is reached, open_list_MB = %ld, closed_list_MB = %ld, max_closed_list = %ld\n", 
        this->memory_bound_first_iter, this->open_size_MB, this->closed_size_MB, this->max_closed_size);

    if ((int) closed_list_linear_gap.size() >= this->max_closed_size) 
        printf("Maximal closed list size reached! Insufficient memory space\n");

    printf("Open list count = %g, cur open list count = %g, closed list count = %g, cur closed list count = %ld, iteration count = %g, nodes pruned by gscore = %g, nodes pruned by hscore = %g\n", 
        this->workload_recorder.open_list_cnt, this->workload_recorder.cur_open_list_cnt, this->workload_recorder.closed_list_cnt, closed_list_linear_gap.size(), this->workload_recorder.iter_cnt, this->workload_recorder.node_pruned_cnt_gscore, this->workload_recorder.node_pruned_cnt_hscore);
    printf("Beam search iter = %d, Astar search iter = %d, Bin count = %d, Insertions of each bin = [", this->beam_width, this->astar_iter_cnt, this->bin_cnt);
    for (int i = 0; i < this->bin_cnt; ++i) {
        if (i == this->bin_cnt - 1) printf("%g", bin_op_cnt[i].first);
        else printf("%g, ", bin_op_cnt[i].first);
    }
    printf("]\n");
    printf("Erasions (excluding the pruned nodes) of each bin = [");
    for (int i = 0; i < this->bin_cnt; ++i) {
        if (i == this->bin_cnt - 1) printf("%g", bin_op_cnt[i].second);
        else printf("%g, ", bin_op_cnt[i].second);
    }
    printf("]\n");
    // print the current memory consumption
    printf("Memory usage = %g bytes, memory threshold = %g bytes\n", (double) get_memory_usage(), 0.8 * this->physical_RAM_size);
}

template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
inline void AnytimeAstarLinearGapSolver<OLMultiIdx, OLFscore, OLCoord, NN>::write_file_log_linear_gap() {
    std::string line = "Recursive MSA completed.\n";
    this->file_writer->write(line);

    // printf("Iter %d: First time memory bound is reached, open_list_MB = %d, closed_list_MB = %d, max_closed_list = %d\n", 
    //     this->memory_bound_first_iter, this->open_size_MB, this->closed_size_MB, this->max_closed_size);
    line = "Iter "; line += std::to_string(this->memory_bound_first_iter); line += ": First time memory bound is reached, open_list_MB = ";
    line += std::to_string(this->open_size_MB) + ", closed_list_MB = " + std::to_string(this->closed_size_MB) + ", max_closed_size = " + std::to_string(this->max_closed_size) + "\n";
    this->file_writer->write(line);

    if ((int) closed_list_linear_gap.size() >= this->max_closed_size) {
        line = "Maximal closed list size reached! Insufficient memory space\n";
        this->file_writer->write(line);
    }

    // printf("Open list count = %g, cur open list count = %g, closed list count = %g, cur closed list count = %ld, iteration count = %g, nodes pruned by gscore = %g, nodes pruned by hscore = %g\n", 
    //     this->workload_recorder.open_list_cnt, this->workload_recorder.cur_open_list_cnt, this->workload_recorder.closed_list_cnt, closed_list_linear_gap.size(), this->workload_recorder.iter_cnt, this->workload_recorder.node_pruned_cnt_gscore, this->workload_recorder.node_pruned_cnt_hscore);
    line = "Total open list insertions = " + std::to_string(this->workload_recorder.open_list_cnt);
    line += ", cur open list count = " + std::to_string(this->workload_recorder.cur_open_list_cnt);
    line += ", total closed list insertions = " + std::to_string(this->workload_recorder.closed_list_cnt);
    line += ", cur closed list count = " + std::to_string(closed_list_linear_gap.size());
    line += ", iteration count = " + std::to_string(this->workload_recorder.iter_cnt);
    line += ", nodes pruned by gscore = " + std::to_string(this->workload_recorder.node_pruned_cnt_gscore);
    line += ", nodes pruned by hscore = " + std::to_string(this->workload_recorder.node_pruned_cnt_hscore) + "\n";
    this->file_writer->write(line);

    // printf("Beam search iter = %d, Astar search iter = %d, Bin count = %d", this->beam_width, this->astar_iter_cnt, this->bin_cnt);
    line = "Beam search iter = " + std::to_string(this->beam_width);
    line += ", Astar search iter = " + std::to_string(this->astar_iter_cnt);
    line += ", Bin count = " + std::to_string(this->bin_cnt) + "\n";
    this->file_writer->write(line);

    line = "Memory usage = " + std::to_string(get_memory_usage()) + " bytes, memory threshold = " + std::to_string(0.8 * this->physical_RAM_size) + " bytes\n";
    this->file_writer->write(line);

    line = "All PSA time = " + std::to_string(this->workload_recorder.all_PSA_time) + " s, greedy MSA time = " + std::to_string(this->workload_recorder.greedy_MSA_time) + " s\n";
    this->file_writer->write(line);

    line = "Recursive MSA (linear gap) time = [";
    for (int i = 3; i <= this->nn; ++i) {
        line += std::to_string(i) + "D: " + std::to_string(this->workload_recorder.recursive_MSA_time[i]) + " s";
        if (i != this->nn) line += ", ";
    }
    line += "]\n";
    this->file_writer->write(line);


    // total time
    double total_time = this->time_limit_timer->get_time_from_start();
    line = "Total time: " + std::to_string(total_time) + " s\n";
    this->file_writer->write(line);
}


template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
inline double AnytimeAstarLinearGapSolver<OLMultiIdx, OLFscore, OLCoord, NN>::AnytimeAstar_MSA_linear_gap(char **&msa_result, int &alignment_len, bool &reached_time_limit_output) {
    this->time_limit_timer->start(this->time_limit);

    double access_cnt = 0;  // unused

    Timer MSA_timer;
    MSA_timer.start();
    this->pairwise_heuristic->compute_all_PSA();        // pre-compute all pairwise alignments
    this->workload_recorder.all_PSA_time = MSA_timer.elapsed(false);

    printf("All-to-all PSA done!\n");

    bool greedy_reverse_seq = false;
    bool ordered_by_fscore_instead_of_gscore = true;
    compute_greedy_linear_gap(this->nn, greedy_reverse_seq, ordered_by_fscore_instead_of_gscore);
    this->workload_recorder.greedy_MSA_time = MSA_timer.elapsed(false);
    printf("Greedy MSA done!\n");

    bool read_closed_list_from_file = false;        // valid when this->enable_recursive_MSA == true
    double closed_list_write_time = 0.0, closed_list_read_time = 0.0;

    if (this->enable_recursive_MSA) {
        std::string closed_dir = this->sequences_dir;
        closed_dir.append("_recursive_closed_list");
        std::fstream closed_list_file_in;
        if (this->enable_file_storage) {
            closed_list_file_in.open(closed_dir, std::ios_base::in);
            if (closed_list_file_in.is_open()) {
                read_closed_list_from_file = true;
                printf("read_closed_list_from_file = %d\n", read_closed_list_from_file);
            }
        }

        if (read_closed_list_from_file == false) {      // compute recursive MSA from 3D to (n-1)-D
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

            if (this->enable_file_storage) {
                // write the closed list into a file
                std::fstream closed_list_file_out(closed_dir, std::ios::out | std::ios::trunc);
                if (closed_list_file_out.is_open()) {
                    printf("Writing the %dD closed list into file: \"%s\"!\n", (this->nn) - 1, closed_dir.c_str());
                    closed_list_file_out << (this->nn) - 1 << std::endl;      // first write the dimension
                    for (auto closed_iter = closed_list_linear_gap.begin(); closed_iter != closed_list_linear_gap.end(); ++closed_iter) {
                        for (int d = 0; d < (this->nn) - 1; ++d) {
                            closed_list_file_out << closed_iter->first[d] << " ";
                        }
                        closed_list_file_out << closed_iter->second.gscore << " " << closed_iter->second.reexp_fscore << " " << closed_iter->second.came_from << std::endl;
                    }
                    closed_list_file_out.close();
                    closed_list_write_time = MSA_timer.elapsed(false);
                    printf("Successfully written the %dD closed list into file: \"%s\"! Writing time = %g s\n", (this->nn) - 1, closed_dir.c_str(), closed_list_write_time);
                }
            }
        } else if (this->enable_file_storage) {            // read the (nn-1)-D closed list
            std::string line;
            std::getline(closed_list_file_in, line);
            int dim_in_file = stoi(line);
            printf("Reading closed list from the file, dim (nn-1) = %d\n", dim_in_file);
            NodeCoord key;
            for (int i = 0; i < NN; ++i) key[i] = 0;

            ClosedListNodeLinearGap value;

            while (std::getline(closed_list_file_in, line)) {
                // printf("Read a new line: %s\n", line.c_str());
                std::istringstream iss(line);
                for (int d = 0; d < (this->nn) - 1; ++d) {
                    iss >> key[d];
                }
                iss >> value.gscore >> value.reexp_fscore >> value.came_from;
                // printf("Before writing %s\n", line.c_str());
                closed_list_linear_gap[key] = value;
            }
            closed_list_file_in.close();
            printf("Successfully read the %dD closed list from file: %s!\n", (this->nn) - 1, closed_dir.c_str());
            closed_list_read_time = MSA_timer.elapsed(false);
        }
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

    printf("All PSA time = %f s, greedy MSA time = %f s\n", this->workload_recorder.all_PSA_time, this->workload_recorder.greedy_MSA_time);
    printf("Recursive MSA (linear gap) time = [");
    for (int i = 3; i <= this->nn; ++i) {
        if (i == this->nn) printf("%dD: %g s", i, this->workload_recorder.recursive_MSA_time[i]);
        else printf("%dD: %g s, ", i, this->workload_recorder.recursive_MSA_time[i]);
    }
    printf("]\n");

    if (this->enable_file_storage) {
        if (read_closed_list_from_file) printf("Closed list file reading time = %f s\n", closed_list_read_time);
        else                            printf("Closed list file writing time = %f s\n", closed_list_write_time);
    }


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

    write_file_log_linear_gap();

    // deallocate the closed list
    closed_list_linear_gap.clear();

    return access_cnt;
}


template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
inline STYPE AnytimeAstarLinearGapSolver<OLMultiIdx, OLFscore, OLCoord, NN>::score_upper_bound_2D_linear_gap(bool access_reversed_mats, int seq_idx[], int seq_offset[], int cur_seq_lens[]) {
    STYPE score_upper_bound = 0;
    // if both have len = 0, return 0
    if (cur_seq_lens[0] == 0 && cur_seq_lens[1] == 0) return 0.0;
    else if (cur_seq_lens[0] == 0) {
        STYPE pen_e = this->blosum_table.get_ext_penalty();
        score_upper_bound -= pen_e * cur_seq_lens[1];
        return score_upper_bound;
    } else if (cur_seq_lens[1] == 0){
        STYPE pen_e = this->blosum_table.get_ext_penalty();
        score_upper_bound -= pen_e * cur_seq_lens[0];
        return score_upper_bound;
    }

    STYPE DP_val;
    int TB_val;     // unused, because there's no gap open penalty
    this->pairwise_heuristic->get_score(access_reversed_mats, seq_idx[0], seq_idx[1], seq_offset[0], seq_offset[1], DP_val, TB_val);
    score_upper_bound = DP_val;
    return score_upper_bound;
};


template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
inline STYPE AnytimeAstarLinearGapSolver<OLMultiIdx, OLFscore, OLCoord, NN>::score_upper_bound_recursive_linear_gap(int dim, bool access_reversed_mats, NodeCoord cur_crd) {

    STYPE res = 0;        // we first compute the recursive score (without the last sequence), then add the PSA scores of the last sequence
    
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
            int seq_idx[2] = {seq_i, seq_j}, seq_offset[2] = {cur_crd[seq_i], cur_crd[seq_j]}, cur_seq_lens[2] = {this->seq_lens[seq_i] - cur_crd[seq_i], this->seq_lens[seq_j] - cur_crd[seq_j]};
            STYPE cell_val = score_upper_bound_2D_linear_gap(access_reversed_mats, seq_idx, seq_offset, cur_seq_lens);

            hscore_last_seq += cell_val;
        }
        res += hscore_last_seq;
    } else {
        for (int seq_i = 0; seq_i < dim - 1; ++seq_i) {
            for (int seq_j = seq_i + 1; seq_j < dim; ++seq_j) {
                int seq_idx[2] = {seq_i, seq_j}, seq_offset[2] = {cur_crd[seq_i], cur_crd[seq_j]}, cur_seq_lens[2] = {this->seq_lens[seq_i] - cur_crd[seq_i], this->seq_lens[seq_j] - cur_crd[seq_j]};
                STYPE cell_val = score_upper_bound_2D_linear_gap(access_reversed_mats, seq_idx, seq_offset, cur_seq_lens);
                res += cell_val;
            }
        }
    }
    return res;
}

// write the result into msa_result and return alignment length
// print_path == true: only print the backtrack path in terminal; print_path == false: output the alignment
template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
inline int AnytimeAstarLinearGapSolver<OLMultiIdx, OLFscore, OLCoord, NN>::backtrack_linear_gap(char **&msa_result, typename AstarMultiIndex<NN>::template CoordTBInfo<LINEAR_TB_TYPE> crd_tb_info, bool print_path) {   // codeupdate
    NodeCoord seq_indices = crd_tb_info.crd;        // indices of the TB_mat, use --idx{i} when accessing sequences[][]

    int res_idx = 0;                                // counter of the alignment length
    bool all_zero_indices = true;
    for (int i = 0; i < this->nn; ++i)
        if (seq_indices[i] != 0) {
            all_zero_indices = false;
            break;
        }
    
    DIRTYPE parent_dir = crd_tb_info.tb_info;           // a single DIRTYPE variable instead of a vector

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
        parent_dir = closed_list_linear_gap[seq_indices].came_from;
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
