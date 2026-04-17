#ifndef ANYTIMEASTAR_HPP
#define ANYTIMEASTAR_HPP

#include "../parameters.hpp"
#include "../utils.hpp"

#include <vector>
#include <utility>      // pair
#include <queue>
#include <unordered_map>
#include <iostream>
#include <cmath>        // pow, log
#include <algorithm>    // max

#include "PairwiseHeuristic.hpp"        // used when upper bound method = PSA
#include "MultiIndexUtils.hpp"

template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
class AnytimeAstarSolver {
public:
    // store the coordinates as integers
    using NodeCoord = std::array<CRDTYPE, NN>;
    using GapLen = std::array<CRDTYPE, NN>;

protected:
    std::string sequences_dir;              // input file directory, used by recursive MSA to read and write the closed list file

    // open_list_cnt, iter_cnt, node_pruned_cnt_gscore, node_pruned_cnt_hscore
    struct WorkloadRecorder {
        std::vector<float> recursive_MSA_time;     // recur_time[i] represents the time of i-D MSA
        size_t open_list_cnt;               // nodes pushed into the open_list
        size_t cur_open_list_cnt;           // current node count in the open_list
        size_t closed_list_cnt;             // nodes pushed into the closed_list
        size_t iter_cnt;
        size_t node_pruned_cnt_gscore;      // nodes pruned by gscore before inserting into the open_list
        size_t node_pruned_cnt_hscore;      // nodes pruned by hscore before inserting into the open_list
        float all_PSA_time;                // time of computing all PSA results
        float greedy_MSA_time;             // time of computing greedy MSA
    };
    WorkloadRecorder workload_recorder;     // initialized when instantiation

    std::vector<std::pair<float, float>> bin_op_cnt;         // count insertions and erasions of each bin

    // compute the path cost instead of the score
    bool cost_instead_of_score;

    // the node when storing the gap lengths
    // ClosedListGapLen: key = coordinates, value = ClosedListNodeGapLen
    // ClosedListNodeGapLen: key = gap_len, value = {gscore, reexp_fscore}
    using ClosedListNodeGapLen = std::unordered_map<GapLen, std::pair<STYPE, STYPE>, VectorIntHash<NN>>;

    bool use_linear_gap_penalty;        // if gap open penalty == 0, switch to linear gap penalty mode
    const int nn;   // seq_cnt
    const int ll;   // seq_max_len
    int *seq_lens;  // sequence lengths
    char **sequences;           // size = [nn][ll]
    bool enable_recursive_MSA;          // Use the closed list of lower dimensional MSAs to improve the heuristic function
    bool sum_of_crd_as_level;           // In anytime column search, use the sum of coordinates as the level

    char **tmp_NW_sequences;   // array for computing NW PSA, size = [2][ll]
    char **tmp_NW_psa_result;  // NW PSA results, size = [2][2*ll]

    // Pairwise Heuristic, used when upper bound method = PSA
    PSA_heuristic *pairwise_heuristic;

    const ScoreTable &blosum_table;

    char **tmp_msa_result;      // updated when the target node is reached (could be a non-optimal result)
    STYPE tmp_msa_gscore;       // score lower bound updated with the greedy search or beam search result
    int tmp_msa_alignment_length;

    Timer *time_limit_timer;    // Timer begins at the beginning of AnytimeAstar_MSA
    float time_limit;          // When time_limit < 0, there is no time limit
    bool reached_time_limit;    // init to false

    bool enable_memory_bound = true;
    float physical_RAM_size;   // get from sysconf()
    float memory_limit_ratio;        // When memory_limit_ratio < 0 or > 1, there is no memory limit
    size_t open_node_size;         // compare open_node_size * open_cnt + closed_node_size * closed_cnt with memory limit
    size_t closed_node_size;
    // record the iter count when memory bound is reached for the first time
    bool memory_bound_trigger = false;
    size_t memory_bound_first_iter = 0;
    // use the following two variables to verify the memory threshold
    size_t open_size_MB = 0;           // assigned when memory_bound_trigger is set for the first time
    size_t closed_size_MB = 0;
    // updated once open_size_MB and closed_size_MB are initialized
    size_t max_closed_size = INT_MAX;

    // Used to initialize the reexp_fscore in the closed list. Also used to deterine if a node is reexpanded.
    STYPE init_reexp_fscore;

    FileWriter *file_writer;    // output intermediate msa results to a file

    // the above variables are initialized in instantiation

    // Data structures for anytime A*
    int beam_width, astar_iter_cnt;       // hyper-parameters
    int bin_cnt;      // count of open_lists in open_list_arr
    // Initialized after getting the first hscore
    STYPE bin_score_offset = 0;
    STYPE bin_score_thres;
    std::vector<OLMultiIdx> open_list_arr;
    std::vector<OLFscore*> index_by_fscore_arr;
    std::vector<OLCoord*> index_by_coord_arr;

    bool open_list_arr_empty() {
        for (int i = 0; i < bin_cnt; ++i) {
            if (index_by_fscore_arr[i]->empty() == false) return false;
        }
        return true;
    }

    // query a specific bin
    bool open_list_arr_empty(int specific_bin_idx) {
        std::vector<OLFscore*>* target_list = &index_by_fscore_arr;

        if (specific_bin_idx < 0 || specific_bin_idx >= bin_cnt) {
            printf("Invalid specific_bin_idx in open_list_arr_empty()!\n");
            return false;
        }
        return (*target_list)[specific_bin_idx]->empty();
    }

    STYPE get_highest_fscore() {
        std::vector<OLFscore*>* target_list = &index_by_fscore_arr;

        STYPE highest_fscore;
        if (cost_instead_of_score)  highest_fscore = INT_MAX;
        else                        highest_fscore = INT_MIN;

        for (int i = bin_cnt - 1; i >= 0; --i) {
            if ((*target_list)[i]->empty() == false) {
                STYPE cur_fscore = (*target_list)[i]->begin()->fscore;
                bool better_result = cost_instead_of_score ? (cur_fscore < highest_fscore) : (cur_fscore > highest_fscore);
                if (better_result) highest_fscore = cur_fscore;
            }
        }
        return highest_fscore;
    }

    // also outputs target_bin_level
    STYPE get_highest_fscore(int &target_bin_level) {
        std::vector<OLFscore*>* target_list = &index_by_fscore_arr;

        STYPE highest_fscore;
        if (cost_instead_of_score)  highest_fscore = INT_MAX;
        else                        highest_fscore = INT_MIN;

        int max_bin_level = -1;
        for (int i = bin_cnt - 1; i >= 0; --i) {
            if ((*target_list)[i]->empty() == false) {
                STYPE cur_fscore = (*target_list)[i]->begin()->fscore;
                bool better_result = cost_instead_of_score ? (cur_fscore < highest_fscore) : (cur_fscore > highest_fscore);
                if (better_result) {
                    highest_fscore = cur_fscore;
                    max_bin_level = i;
                }
            }
        }
        target_bin_level = max_bin_level;
        if (max_bin_level < 0) printf("Warning! target_bin_level < 0 in get_highest_fscore()\n");
        return highest_fscore;
    }

    // return false if all cur_crd[d] < seq_lens[d]
    bool is_crd_out_of_bound(const NodeCoord &crd, const int dim);

    // pre-load the match / mismatch scores of each sequence pair
    void get_match_mismatch_scores(STYPE *mm_scores, const NodeCoord &cur_crd, const int dim, const bool reverse_seq);

    // return true if 
    // actual_memory_usage == false: the ESTIMATED memory usage exceeds the predefined memory threshold
    // actual_memory_usage == true: actual memory usage > 0.9 * physical RAM size
    // Automatically read the open node cnt from workload_recorder. Receive the closed node cnt from passed parameter.
    bool check_memory_thres(bool actual_memory_usage, size_t closed_node_cnt) {
        bool estimated_res = false, actual_res = false;
        
        if (actual_memory_usage) {
            actual_res = (float) get_memory_usage() > physical_RAM_size * memory_limit_ratio;
        } else {
            float cur_size = workload_recorder.cur_open_list_cnt * open_node_size + closed_node_cnt * closed_node_size;
            float MB_size = (float) open_size_MB * open_node_size + closed_size_MB * closed_node_size;
            estimated_res = cur_size > MB_size;
            // printf("cur_size = %f, MB_size = %f, estimated_res = %d\n", cur_size, MB_size, (int) estimated_res);
        }

        if (estimated_res | actual_res) {
            return true;
        } else
            return false;
    }

    // each dimension of index has a max length of ll+1, while sequences has a size of ll
    // VectorIntHash in utils.hpp
    // key = coord, value = unordered_map<gap_len, gscore>
    std::unordered_map<NodeCoord, ClosedListNodeGapLen, VectorIntHash<NN>> closed_list_gap_len;

    /* --- recursive MSA --- */
    // key = coord, value = unordered_map<gap_len, gscore>
    std::unordered_map<NodeCoord, ClosedListNodeGapLen, VectorIntHash<NN>> recursive_closed_list_gap_len;
    STYPE score_upper_bound_recursive(int dim, bool access_reversed_mats, const NodeCoord &cur_crd, const GapLen &gap_len);
    void compute_recursive_Astar(int seq_cnt, bool reverse_seq);
    /* --- recursive MSA --- */

    // prune the unnecessary nodes in the closed list when reaching the memory bound and return the prunable count
    int prune_closed_list(bool access_reversed_mats, int dim);

    bool closed_list_duplication_detection(bool &use_reexp_fscore_instead, const NodeCoord &nbr_crd, const GapLen &nbr_tb_info, STYPE cur_gscore, STYPE &cur_fscore);

    void open_list_insertion(STYPE cur_gscore, STYPE cur_fscore, std::vector<OLCoord*>* target_index_by_coord, const NodeCoord &nbr_crd, const GapLen &nbr_tb_info);
    

    STYPE score_upper_bound_2D(bool access_reversed_mats, int seq_idx[], int seq_offset[], bool gap_in_prev_seq[], int cur_seq_lens[]);     // return the score_upper_bound
    
    void preprocessing_MSA_with_PSA_and_greedy();

    // initialize the tmp_msa_gscore with DFS
    // Called by both linear and affine versions
    void compute_greedy(int dim, bool reverse_seq);
    
    // If the current alignment result is optimal, return true. Otherwise, return false
    // reverse_seq: read the symbol in each sequences from the end instead of the beginning to generate the recursive closed list
    // specific_bin = false: find the node with the largest fscore in all bins 
    bool expand_node(int dim, bool specific_bin, int bin_level, bool reverse_seq);

    void compute_anytime_Astar();
    int backtrack_affine(char **&msa_result, typename AstarMultiIndex<NN>::CoordTBInfo crd_tb_info, bool print_path);        // return the alignment length
    
    /* --- memory bound extension --- */
    void check_memory_bound_trigger();
    GapLen find_parent_in_closed_list(const NodeCoord &last_node_crd, const GapLen &last_node_tb_info, STYPE last_node_gscore);

    bool is_gap(const GapLen &tb_info, int seq_idx) {return tb_info[seq_idx] > 0;};

    void print_progress_per_1M_iter();
    void print_computation_time_after_search();
    void print_workload_report_after_search(int dim);
    void write_file_log();

public:
    virtual size_t get_closed_list_size() {return closed_list_gap_len.size();};

    // _bin_cnt = count of open_lists in open_list_arr
    AnytimeAstarSolver(bool _cost_instead_of_score, int sequence_count, int sequence_max_length, char **original_sequences, int *sequence_lengths, const ScoreTable &_blosum_table, STYPE mafft_score, 
            int _bin_cnt, int _beam_width, int _astar_iter_cnt, float _time_limit, float _memory_limit_ratio, std::string anytime_result_dir, bool _enable_recursive_MSA, bool _sum_of_crd_as_level, std::string input_dir):
        cost_instead_of_score{_cost_instead_of_score}, nn{sequence_count}, ll{sequence_max_length}, sequences{original_sequences}, blosum_table{_blosum_table}, tmp_msa_gscore{mafft_score}, 
            bin_cnt{_bin_cnt}, beam_width{_beam_width}, astar_iter_cnt{_astar_iter_cnt}, time_limit{_time_limit}, memory_limit_ratio{_memory_limit_ratio}, enable_recursive_MSA{_enable_recursive_MSA}, sum_of_crd_as_level{_sum_of_crd_as_level},
            sequences_dir{input_dir} {

        if (blosum_table.get_open_penalty() == 0) use_linear_gap_penalty = true;
        else use_linear_gap_penalty = false;

        seq_lens = new int[nn];

        if (nn <= 3) enable_recursive_MSA = false;          // recursive MSA is only valid when nn > 3

        // allocate memory for the float-buffered closed lists
        // if we use linear gap penalty, allocate the memory for closed_list_linear_gap in other functions

        for (int i = 0; i < nn; ++i) {
            seq_lens[i] = sequence_lengths[i];
        }

        if (sum_of_crd_as_level) {
            bin_cnt = 1;        // bin_cnt = 1 + sum(seq_lens)
            for (int i = 0; i < nn; ++i) bin_cnt += seq_lens[i];
            printf("bin_cnt = %d\n", bin_cnt);
        }

        tmp_NW_sequences = new char*[2];
        tmp_NW_psa_result = new char*[2];
        tmp_NW_psa_result[0] = new char[2 * ll]; tmp_NW_psa_result[1] = new char[2 * ll];
        pairwise_heuristic = new PSA_heuristic(use_linear_gap_penalty, cost_instead_of_score, nn, ll, sequences, seq_lens, blosum_table);

        // initialize tmp_msa_result
        tmp_msa_result = new char*[nn];
        for (int i = 0; i < nn; ++i)
            tmp_msa_result[i] = new char[nn * ll];

        // init the result with gaps, assuming the size of the msa_result is [nn][nn*ll]
        for (int row = 0; row < nn; ++row)
            for (int col = 0; col < nn * ll; ++col) {
                tmp_msa_result[row][col] = GAP;
            }
        tmp_msa_alignment_length = 0;

        if (memory_limit_ratio < 0.0 || memory_limit_ratio > 1.0) {
            enable_memory_bound = false;
            printf("Memory bound extension disabled!\n");
        } else {
            // long pages = sysconf(_SC_PHYS_PAGES);
            // long page_size = sysconf(_SC_PAGE_SIZE);
            // physical_RAM_size = (float) pages * page_size;

            physical_RAM_size = get_available_memory();

            printf("Physical RAM size = %f bytes\n", physical_RAM_size);

            if (use_linear_gap_penalty) {
                // theoretically
                // open_node_size = 2 * NN + 8 + 4 + 4
                // closed_node_size = ?
                // practically
                open_node_size = 92 + ((nn - 1) / 8) * 16;
                closed_node_size = 60 + ((nn + 3) / 8) * 16;
            } else {
                // theoretically
                // open_node_size = 4 * NN + 4 + 4
                // closed_node_size = ?
                // practically
                open_node_size = 92 + ((nn - 1) / 4) * 16;
                closed_node_size = 268 + ((nn - 1) / 8) * 32;
            }
            init_reexp_fscore = cost_instead_of_score ? INT_MAX : INT_MIN;
        }

        time_limit_timer = new Timer();
        file_writer = new FileWriter(anytime_result_dir, nn, input_dir, sequences, seq_lens, blosum_table.get_open_penalty(), blosum_table.get_ext_penalty());

        workload_recorder = {
            std::vector<float>(nn + 1, 0),     // the time of computing recursive MSA
            0,  // open_list_cnt
            0,  // cur_open_list_cnt
            0,  // closed_list_cnt
            1,  // iter_cnt
            0,  // node_pruned_cnt_gscore
            0,  // node_pruned_cnt_hscore
            0,  // all PSA
            0   // greedy MSA
        };

        for (int i = 0; i < bin_cnt; ++i)
            bin_op_cnt.push_back({0, 0});
    }

    ~AnytimeAstarSolver() {
        delete [] seq_lens;
        delete [] tmp_NW_sequences;
        delete [] tmp_NW_psa_result[0]; delete [] tmp_NW_psa_result[1]; delete [] tmp_NW_psa_result;

        for (int i = 0; i < nn; ++i)
            delete [] tmp_msa_result[i];
        delete [] tmp_msa_result;

        closed_list_gap_len.clear();
        recursive_closed_list_gap_len.clear();

        delete time_limit_timer;
        delete file_writer;

        delete pairwise_heuristic;

        for (int i = 0; i < bin_cnt; ++i) {
            open_list_arr[i].clear();
        }
        open_list_arr.clear();
        index_by_fscore_arr.clear();
        index_by_coord_arr.clear();
    }

    // return the access count and the alignment length
    void AnytimeAstar_MSA(char **&msa_result, int &alignment_len, bool &reached_time_limit_output) {
        time_limit_timer->start(time_limit);

        float access_cnt = 0;  // unused
        preprocessing_MSA_with_PSA_and_greedy();

        Timer MSA_timer;

        if (enable_recursive_MSA) {
            // recursive Astar has its own tmp_msa_gscore
            STYPE prev_tmp_msa_gscore = tmp_msa_gscore;
            for (int recur_idx = 3; recur_idx < nn; ++recur_idx) {
                // compute recur_idx dimensional MSA
                bool reverse_seq = (nn - recur_idx) & 0x01;         // reverse the sequences when computing (nn-1)-D MSA, (nn-3)-D MSA, ...
                
                tmp_msa_gscore = cost_instead_of_score ? INT_MAX : INT_MIN;
                memory_bound_trigger = false;

                MSA_timer.start();
                compute_recursive_Astar(recur_idx, reverse_seq);
                workload_recorder.recursive_MSA_time[recur_idx] = MSA_timer.elapsed(false);
            }
            tmp_msa_gscore = prev_tmp_msa_gscore;
        }


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

        memory_bound_trigger = false;

        MSA_timer.start();
        compute_anytime_Astar();
        workload_recorder.recursive_MSA_time[nn] = MSA_timer.elapsed(false);

        print_computation_time_after_search();

        // init the result with gaps, assuming the size of the msa_result is [nn][nn*ll]
        for (int row = 0; row < nn; ++row)
            for (int col = 0; col < nn * ll; ++col) {
                msa_result[row][col] = GAP;
            }

        // backtrack is done whenever the target node is reached in expand_node(), just copy the results from tmp_msa_result to msa_result here
        for (int row = 0; row < nn; ++row) {
            for (int col = 0; col < tmp_msa_alignment_length; ++col) {
                msa_result[row][col] = tmp_msa_result[row][col];
            }
        }
        alignment_len = tmp_msa_alignment_length;
        reached_time_limit = time_limit_timer->check_time_limit();
        reached_time_limit_output = reached_time_limit;

        write_file_log();

        return;
    }
};

template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
inline int AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::prune_closed_list(bool access_reversed_mats, int dim) {
    printf("Start searching prunable closed list node, closed list size = %ld\n", get_closed_list_size());
    int prunable_cnt = 0;
    for (auto closed_iter = closed_list_gap_len.begin(); closed_iter != closed_list_gap_len.end();) {
        NodeCoord cur_crd = closed_iter->first;
        for (auto inner_iter = closed_iter->second.begin(); inner_iter != closed_iter->second.end(); ) {
            GapLen cur_tb_info = inner_iter->first;
            STYPE closed_gscore = inner_iter->second.first;       // inner_iter->second = std::pair{gscore, reexp_fscore}
            STYPE closed_hscore = score_upper_bound_recursive(dim, access_reversed_mats, cur_crd, cur_tb_info);

            if (cost_instead_of_score && (closed_gscore + closed_hscore > tmp_msa_gscore) ||
                !cost_instead_of_score && (closed_gscore + closed_hscore < tmp_msa_gscore)) {
                inner_iter = closed_iter->second.erase(inner_iter);
                prunable_cnt += 1;
            } else {
                ++inner_iter;
            }
        }

        if (closed_iter->second.empty()) {
            closed_iter = closed_list_gap_len.erase(closed_iter);
        } else {
            ++closed_iter;
        }

    }
    printf("Total prunable count = %d, closed list size = %ld\n", prunable_cnt, get_closed_list_size());
    return prunable_cnt;
}

/**
 * @brief compute the alignment score upper bound with the specific method (Same with that in AstarSolver)
 * 
 * @param seq_idx the index of sequences[][]
 * @param seq_offset the beginning index of sequences[seq_idx][]
 * @param prev_seq_offset similar to traceback flag for resolving gap open penalties when upper_bound_method = PSA
 * @param cur_seq_lens the length of the sequences
 * @return return the score_upper_bound
 */
template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
inline STYPE AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::score_upper_bound_2D(bool access_reversed_mats, int seq_idx[], int seq_offset[], bool gap_in_prev_seq[], int cur_seq_lens[]) {
    STYPE score_upper_bound = 0;
    // if both have len = 0, return 0
    if (cur_seq_lens[0] == 0 && cur_seq_lens[1] == 0) return 0.0;
    else if (cur_seq_lens[0] == 0) {
        STYPE pen_e = blosum_table.get_ext_penalty(), pen_o = blosum_table.get_open_penalty();
        score_upper_bound -= pen_e * cur_seq_lens[1];
        score_upper_bound -= pen_o * (gap_in_prev_seq[0] == false);     // apply gap open penalty only when prev_seq has no gap
        return score_upper_bound;
    } else if (cur_seq_lens[1] == 0){
        STYPE pen_e = blosum_table.get_ext_penalty(), pen_o = blosum_table.get_open_penalty();
        score_upper_bound -= pen_e * cur_seq_lens[0];
        score_upper_bound -= pen_o * (gap_in_prev_seq[1] == false);     // apply gap open penalty only when prev_seq has no gap
        return score_upper_bound;
    }

    score_upper_bound = pairwise_heuristic->get_score(access_reversed_mats, seq_idx[0], seq_idx[1], seq_offset[0], seq_offset[1]);

    // simpler duplicated gap open detection: do "score += open_pen" when there is a gap in the current node
    if (gap_in_prev_seq[0] ^ gap_in_prev_seq[1]) 
        score_upper_bound += blosum_table.get_open_penalty();


    return score_upper_bound;
};



template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
inline STYPE AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::score_upper_bound_recursive(int dim, bool access_reversed_mats, const NodeCoord &cur_crd, const GapLen &gap_len) {
    STYPE open_pen = blosum_table.get_open_penalty();
    STYPE res = 0;        // we first compute the recursive score (without the last sequence), then add the PSA scores of the last sequence
    
    if (enable_recursive_MSA) {
        // first find the crd in the recursive closed list
        NodeCoord reversed_crd;        // use dim - 1 because recursive_closed_list_gap_len contains (dim - 1)-D results
        for (int i = 0; i < NN; ++i) reversed_crd[i] = 0;
        for (int i = 0; i < dim - 1; ++i) reversed_crd[i] = seq_lens[i] - cur_crd[i];

        auto recur_closed_list_iter = recursive_closed_list_gap_len.find(reversed_crd);
        // compare the (dim-1)*(dim-1) distance matrix
        if (recur_closed_list_iter != recursive_closed_list_gap_len.end()) {
            // init the res
            if (cost_instead_of_score)  res = INT_MAX;
            else                        res = INT_MIN;

            // subtract the duplicated gap open penalty using the gap lengths
            for (auto node_iter = recur_closed_list_iter->second.begin(); node_iter != recur_closed_list_iter->second.end(); ++node_iter) {
                GapLen recur_gap_len = node_iter->first;
                STYPE recur_gscore = node_iter->second.first;       // node_iter->second = std::pair{gscore, reexp_fscore}
                // the recursive closed list has no information of seq(dim - 1)
                for (int seq_i = 0; seq_i < dim - 2; ++seq_i) {
                    for (int seq_j = seq_i + 1; seq_j < dim - 1; ++seq_j) {
                        // simpler duplicated gap open detection: do "score += open_pen" when there is a gap in the current node
                        DIRTYPE cur_gap_diff = gap_len[seq_i] - gap_len[seq_j];
                        if (cur_gap_diff != 0) 
                            recur_gscore += open_pen;       // we did -= open_pen twice, so here we += open_pen
                    }
                }
                if (cost_instead_of_score && (recur_gscore < res) ||
                    !cost_instead_of_score && (recur_gscore > res)) 
                    res = recur_gscore;
            }

            STYPE hscore_last_seq = 0;
            // add the PSA scores of the last sequence (deterministic)
            for (int seq_i = 0; seq_i < dim - 1; ++seq_i) {
                int seq_j = dim - 1;
                int seq_idx[2] = {seq_i, seq_j}, seq_offset[2] = {cur_crd[seq_i], cur_crd[seq_j]}, cur_seq_lens[2] = {seq_lens[seq_i] - cur_crd[seq_i], seq_lens[seq_j] - cur_crd[seq_j]};
                bool gap_in_prev_seq[2] = {false, false};
                gap_in_prev_seq[0] = gap_len[seq_i] > gap_len[seq_j];
                gap_in_prev_seq[1] = gap_len[seq_i] < gap_len[seq_j];
                hscore_last_seq += score_upper_bound_2D(access_reversed_mats, seq_idx, seq_offset, gap_in_prev_seq, cur_seq_lens);
            }
            res += hscore_last_seq;
            return res;
        } 
    }
    
    CRDTYPE pairwise_heuristic_crd[NN];
    for (int i = 0; i < dim; ++i) pairwise_heuristic_crd[i] = seq_lens[i] - cur_crd[i];
    DIRTYPE pairwise_gap_len[NN];
    for (int i = 0; i < dim; ++i) pairwise_gap_len[i] = gap_len[i];

    res = pairwise_heuristic->get_sum_of_scores(dim, access_reversed_mats, pairwise_heuristic_crd, pairwise_gap_len);

    return res;
}


template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
inline void AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::print_computation_time_after_search() {
    printf("All PSA time = %f s, greedy MSA time = %f s\n", workload_recorder.all_PSA_time, workload_recorder.greedy_MSA_time);
    printf("Recursive MSA time = [");
    for (int i = 3; i <= nn; ++i) {
        if (i == nn) printf("%dD: %g s", i, workload_recorder.recursive_MSA_time[i]);
        else printf("%dD: %g s, ", i, workload_recorder.recursive_MSA_time[i]);
    }
    printf("]\n");
}

#include "AnytimeAstar.tpp"

#endif // ANYTIMEASTAR_HPP