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
    using NodeCoord = std::array<int, NN>;

protected:
    std::string sequences_dir;              // input file directory, used by recursive MSA to read and write the closed list file

    // open_list_cnt, iter_cnt, node_pruned_cnt_gscore, node_pruned_cnt_hscore
    struct WorkloadRecorder {
        std::vector<double> recursive_MSA_time;     // recur_time[i] represents the time of i-D MSA
        double open_list_cnt;               // nodes pushed into the open_list
        double cur_open_list_cnt;           // current node count in the open_list
        double closed_list_cnt;             // nodes pushed into the closed_list
        double iter_cnt;
        double node_pruned_cnt_gscore;      // nodes pruned by gscore before inserting into the open_list
        double node_pruned_cnt_hscore;      // nodes pruned by hscore before inserting into the open_list
        double all_PSA_time;                // time of computing all PSA results
        double greedy_MSA_time;             // time of computing greedy MSA
    };
    WorkloadRecorder workload_recorder;     // initialized when instantiation

    // compute the path cost instead of the score
    bool cost_instead_of_score;

    // the node when storing the gap lengths
    // ClosedListGapLen: key = coordinates, value = ClosedListNodeGapLen
    // ClosedListNodeGapLen: key = gap_len, value = gscore
    using ClosedListNodeGapLen = std::unordered_map<std::vector<DIRTYPE>, STYPE, VectorIntHash<NN>>;

    bool use_linear_gap_penalty;        // if gap open penalty == 0, switch to linear gap penalty mode
    const int nn;   // seq_cnt
    const int ll;   // seq_max_len
    int *seq_lens;  // sequence lengths
    char **sequences;           // size = [nn][ll]
    bool enable_recursive_MSA;          // Use the closed list of lower dimensional MSAs to improve the heuristic function

    char **tmp_NW_sequences;   // array for computing NW PSA, size = [2][ll]
    char **tmp_NW_psa_result;  // NW PSA results, size = [2][2*ll]

    // Pairwise Heuristic, used when upper bound method = PSA
    PSA_heuristic *pairwise_heuristic;

    ScoreTable blosum_table;

    char **tmp_msa_result;      // updated when the target node is reached (could be a non-optimal result)
    STYPE tmp_msa_gscore;       // score lower bound updated with the greedy search or beam search result
    int tmp_msa_alignment_length;

    Timer *time_limit_timer;    // Timer begins at the beginning of AnytimeAstar_MSA
    double time_limit;          // When time_limit < 0, there is no time limit
    bool reached_time_limit;    // init to false

    bool enable_memory_bound = true;
    double physical_RAM_size;   // get from sysconf()
    double memory_limit_ratio;        // When memory_limit_ratio < 0 or > 1, there is no memory limit
    size_t open_node_size;         // compare open_node_size * open_cnt + closed_node_size * closed_cnt with memory limit
    size_t closed_node_size;
    // record the iter count when memory bound is reached for the first time
    bool memory_bound_trigger = false;
    int memory_bound_first_iter = -1;
    // use the following two variables to verify the memory threshold
    size_t open_size_MB = 0;           // assigned when memory_bound_trigger is set for the first time
    size_t closed_size_MB = 0;
    // updated once open_size_MB and closed_size_MB are initialized
    size_t max_closed_size = INT_MAX;

    bool enable_file_storage = false;       // read/write the (n-1)-D closed list from/into a file when using recursive MSA
    
    // Used to initialize the reexp_fscore in the closed list. Also used to deterine if a node is reexpanded.
    STYPE init_reexp_fscore;

    FileWriter *file_writer;    // output intermediate msa results to a file

    // the above variables are initialized in instantiation

    // Data structures for anytime A*
    int beam_width, astar_iter_cnt;       // hyper-parameters
    int bin_cnt;      // count of open_lists in open_list_arr
    // Initialized after getting the first hscore
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
        if (cost_instead_of_score)  highest_fscore = __FLT_MAX__;
        else                        highest_fscore = __FLT_MIN__;

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
        if (cost_instead_of_score)  highest_fscore = __FLT_MAX__;
        else                        highest_fscore = __FLT_MIN__;

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
        return highest_fscore;
    }

    // return true if 
    // actual_memory_usage == false: the ESTIMATED memory usage exceeds the predefined memory threshold
    // actual_memory_usage == true: actual memory usage > 0.8 * physical RAM size
    // Automatically read the open node cnt from workload_recorder. Receive the closed node cnt from passed parameter.
    bool check_memory_thres(bool actual_memory_usage, size_t closed_node_cnt) {
        bool estimated_res = false, actual_res = false;
        
        if (actual_memory_usage) {
            actual_res = (double) get_memory_usage() > physical_RAM_size * memory_limit_ratio;
        } else {
            double cur_size = workload_recorder.cur_open_list_cnt * open_node_size + closed_node_cnt * closed_node_size;
            double MB_size = (double) open_size_MB * open_node_size + closed_size_MB * closed_node_size;
            estimated_res = cur_size > MB_size;
        }

        if (estimated_res | actual_res) {
            return true;
        } else
            return false;
    }

    // key = coord, value = unordered_map<gap_len, gscore>
    std::unordered_map<NodeCoord, ClosedListNodeGapLen, VectorIntHash<NN>> closed_list_gap_len;

    /* --- recursive MSA --- */
    // key = coord, value = unordered_map<gap_len, gscore>
    std::unordered_map<NodeCoord, ClosedListNodeGapLen, VectorIntHash<NN>>** recursive_closed_lists;
    STYPE score_upper_bound_recursive(int dim, bool access_reversed_mats, NodeCoord cur_crd, std::vector<DIRTYPE> gap_len, std::unordered_map<NodeCoord, ClosedListNodeGapLen, VectorIntHash<NN>>* recur_closed_list);
    void compute_recursive_Astar(int seq_cnt, bool reverse_seq);
    /* --- recursive MSA --- */

    STYPE score_upper_bound_2D(bool access_reversed_mats, int seq_idx[], int seq_offset[], bool gap_in_prev_seq[], int cur_seq_lens[]);     // return the score_upper_bound
    // If the current alignment result is optimal, return true. Otherwise, return false
    // reverse_seq: read the symbol in each sequences from the end instead of the beginning to generate the recursive closed list
    // specific_bin = false: find the node with the largest fscore in all bins 
    bool expand_node(int dim, bool specific_bin, int bin_level, bool reverse_seq, WorkloadRecorder& workload_recorder, std::vector<std::pair<double, double>>& bin_op_cnt);

    void compute_anytime_Astar();
    // codeupdate
    int backtrack_affine(char **&msa_result, typename AstarMultiIndex<NN>::template CoordTBInfo<AFFINE_TB_TYPE> crd_tb_info, bool print_path);        // return the alignment length
public:
    // _bin_cnt = count of open_lists in open_list_arr
    AnytimeAstarSolver(bool _cost_instead_of_score, int sequence_count, int sequence_max_length, char **original_sequences, int *sequence_lengths, ScoreTable _blosum_table, STYPE mafft_score, 
            int _bin_cnt, int _beam_width, int _astar_iter_cnt, double _time_limit, double _memory_limit_ratio, std::string anytime_result_dir, bool _enable_recursive_MSA, std::string input_dir):
        cost_instead_of_score{_cost_instead_of_score}, nn{sequence_count}, ll{sequence_max_length}, sequences{original_sequences}, blosum_table{_blosum_table}, tmp_msa_gscore{mafft_score}, 
            bin_cnt{_bin_cnt}, beam_width{_beam_width}, astar_iter_cnt{_astar_iter_cnt}, time_limit{_time_limit}, memory_limit_ratio{_memory_limit_ratio}, enable_recursive_MSA{_enable_recursive_MSA}, sequences_dir{input_dir} {

        if (blosum_table.get_open_penalty() == 0) use_linear_gap_penalty = true;
        else use_linear_gap_penalty = false;
        
        seq_lens = new int[nn];

        if (nn <= 3) enable_recursive_MSA = false;          // recursive MSA is only valid when nn > 3

        // allocate memory for the double-buffered closed lists
        // if we use linear gap penalty, allocate the memory for closed_list_linear_gap in other functions

        if (use_linear_gap_penalty == false) {
            recursive_closed_lists = new std::unordered_map<NodeCoord, ClosedListNodeGapLen, VectorIntHash<NN>>*[2];
            for (int i = 0; i < 2; ++i)
                recursive_closed_lists[i] = new std::unordered_map<NodeCoord, ClosedListNodeGapLen, VectorIntHash<NN>>;
        }

        for (int i = 0; i < nn; ++i) {
            seq_lens[i] = sequence_lengths[i];
        }
        tmp_NW_sequences = new char*[2];
        tmp_NW_psa_result = new char*[2];
        tmp_NW_psa_result[0] = new char[2 * ll]; tmp_NW_psa_result[1] = new char[2 * ll];
        pairwise_heuristic = new PSA_heuristic(cost_instead_of_score, nn, ll, sequences, seq_lens, blosum_table);

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
        
        if (memory_limit_ratio < 0.0 || memory_limit_ratio > 1.0) enable_memory_bound = false;
        else {
            if (use_linear_gap_penalty == false)
                printf("Warning! memory threshold is calculated using the node size of linear gap penalty. But now we are using affine gap penalty!\n");
            physical_RAM_size = get_available_memory();

            printf("Physical RAM size = %f bytes\n", physical_RAM_size);

            // For closed list, nn should be padded to a multiple of 2
            // collected with ideal test cases:
            open_node_size = 108;
            closed_node_size = 8 * ((nn + 1) / 2 * 2) + 48;         // pad nn to a multiple of 2
            init_reexp_fscore = this->cost_instead_of_score ? INT_MAX : INT_MIN;
        }

        time_limit_timer = new Timer();
        file_writer = new FileWriter(anytime_result_dir, nn, input_dir, sequences, seq_lens);

        workload_recorder = {
            std::vector<double>(nn + 1, 0),     // the time of computing recursive MSA
            0,  // open_list_cnt
            0,  // cur_open_list_cnt
            0,  // closed_list_cnt
            1,  // iter_cnt
            0,  // node_pruned_cnt_gscore
            0,  // node_pruned_cnt_hscore
            0,  // all PSA
            0   // greedy MSA
        };
    }

    ~AnytimeAstarSolver() {
        delete [] seq_lens;
        delete [] tmp_NW_sequences;
        delete [] tmp_NW_psa_result[0]; delete [] tmp_NW_psa_result[1]; delete [] tmp_NW_psa_result;

        for (int i = 0; i < nn; ++i)
            delete [] tmp_msa_result[i];
        delete [] tmp_msa_result;

        closed_list_gap_len.clear();
        if (use_linear_gap_penalty == false) {
            for (int i = 0; i < 2; ++i)
                delete recursive_closed_lists[i];
            delete [] recursive_closed_lists;
        }

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
    double AnytimeAstar_MSA(char **&msa_result, int &alignment_len, bool &reached_time_limit_output) {
        time_limit_timer->start(time_limit);

        double access_cnt = 0;  // unused

        Timer MSA_timer;
        MSA_timer.start();
        pairwise_heuristic->compute_all_PSA();        // pre-compute all pairwise alignments
        workload_recorder.all_PSA_time = MSA_timer.get_time_from_start();

        if (enable_recursive_MSA) {
            // recursive Astar has its own tmp_msa_gscore
            STYPE prev_tmp_msa_gscore = tmp_msa_gscore;
            for (int recur_idx = 3; recur_idx < nn; ++recur_idx) {
                // compute recur_idx dimensional MSA
                bool reverse_seq = (nn - recur_idx) & 0x01;         // reverse the sequences when computing (nn-1)-D MSA, (nn-3)-D MSA, ...
                
                tmp_msa_gscore = cost_instead_of_score ? INT_MAX : INT_MIN;

                // The results are stored in recursive_closed_lists[0]. Swap recursive_closed_lists[0] and [1] at the beginning of compute_recursive_Astar
                MSA_timer.start();
                compute_recursive_Astar(recur_idx, reverse_seq);
                workload_recorder.recursive_MSA_time[recur_idx] = MSA_timer.get_time_from_start();
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

        MSA_timer.start();
        compute_anytime_Astar();
        workload_recorder.recursive_MSA_time[nn] = MSA_timer.get_time_from_start();

        printf("All PSA time = %f s\n", workload_recorder.all_PSA_time);
        printf("Recursive MSA time = [");
        for (int i = 3; i <= nn; ++i) {
            if (i == nn) printf("%dD: %g s", i, workload_recorder.recursive_MSA_time[i]);
            else printf("%dD: %g s, ", i, workload_recorder.recursive_MSA_time[i]);
        }
        printf("]\n");


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
        return access_cnt;
    }
};

// If the current alignment result is optimal, return true. Otherwise, return false
// reverse_seq: read the symbol in each sequences from the end instead of the beginning to generate the recursive closed list
// specific_bin = false: find the node with the largest fscore in all bins

template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
inline bool AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::expand_node(int dim, bool specific_bin, int bin_level, bool reverse_seq, WorkloadRecorder& workload_recorder, std::vector<std::pair<double, double>>& bin_op_cnt) {
    STYPE ext_pen = blosum_table.get_ext_penalty(), open_pen = blosum_table.get_open_penalty();

    std::vector<OLFscore*>* target_index_by_fscore = &index_by_fscore_arr;
    std::vector<OLCoord*>* target_index_by_coord = &index_by_coord_arr;

    std::unordered_map<NodeCoord, ClosedListNodeGapLen, VectorIntHash<NN>>* target_closed_list = nullptr;
    std::unordered_map<NodeCoord, ClosedListNodeGapLen, VectorIntHash<NN>>* recur_closed_list = nullptr;

    if (dim < nn)   target_closed_list = recursive_closed_lists[0];
    else            target_closed_list = &closed_list_gap_len;

    if (dim < nn)   recur_closed_list = recursive_closed_lists[1];
    else            recur_closed_list = recursive_closed_lists[0];
    bool access_reversed_mats = !reverse_seq;

    // search the node with the highest fscore and update the bin_level
    if (specific_bin == false) {
        // update the bin_level with the bin_idx containing the highest fscore
        get_highest_fscore(bin_level);
    }
    // codeupdate
    typename AstarMultiIndex<NN>::template CoordTBInfo<AFFINE_TB_TYPE> cur_crd_tb_info;
    STYPE local_gscore, local_fscore;

    auto cur_node = (*target_index_by_fscore)[bin_level]->begin();
    cur_crd_tb_info.crd = cur_node->crd_tb_info.crd;
    cur_crd_tb_info.tb_info = cur_node->crd_tb_info.tb_info;
    local_gscore = cur_node->gscore;
    local_fscore = cur_node->fscore;
    (*target_index_by_fscore)[bin_level]->erase((*target_index_by_fscore)[bin_level]->begin());

    workload_recorder.cur_open_list_cnt -= 1;
    
    bin_op_cnt[bin_level].second += 1;
    
    bool target_reached = true;
    for (int d = 0; d < dim; ++d)
        if (cur_crd_tb_info.crd[d] < seq_lens[d]) {
            target_reached = false;
            break;
        }

    // find the cur_node in ClosedList

    // first find the coord
    auto closed_list_gap_len_iter = target_closed_list->find(cur_crd_tb_info.crd);
    bool write_into_closed_list = true;
    // if we found the node, update the gscore
    if (closed_list_gap_len_iter != target_closed_list->end()) {
        // then find the gap_len and update the gscore in ClosedList
        auto gap_len_score_iter = closed_list_gap_len_iter->second.find(cur_crd_tb_info.tb_info);
        if (gap_len_score_iter != closed_list_gap_len_iter->second.end()) {
            bool better_result = cost_instead_of_score ? (local_gscore > gap_len_score_iter->second) : (local_gscore < gap_len_score_iter->second);
            if (better_result) {
                // existing node has a better gscore, read from the closed list
                local_gscore = gap_len_score_iter->second;
                write_into_closed_list = false;
            }
        }
    }
    // if we didn't find the node or cur_node is better, insert it into the closed list
    if (write_into_closed_list) {
        (*target_closed_list)[cur_crd_tb_info.crd][cur_crd_tb_info.tb_info] = local_gscore;
        workload_recorder.closed_list_cnt += 1;
    }

    if (target_reached) {
        // don't call backtrack when dim < nn, i.e., during recursive MSA with a lower dimension
        if (dim == nn) {
            if (cost_instead_of_score && (local_gscore < tmp_msa_gscore) ||
                !cost_instead_of_score && (local_gscore > tmp_msa_gscore)) {
                tmp_msa_alignment_length = backtrack_affine(tmp_msa_result, cur_crd_tb_info, false);
                tmp_msa_gscore = local_gscore;      // update the score lower bound with the beam search result 
                file_writer->write(time_limit_timer->get_time_from_start(), tmp_msa_result, tmp_msa_alignment_length, local_gscore, workload_recorder.iter_cnt);
            }
        }
        bool is_opt_ret = true;
        STYPE global_best_result = get_highest_fscore();
        if (cost_instead_of_score && (local_gscore > global_best_result)) is_opt_ret = false;
        else if (!cost_instead_of_score && (local_gscore < global_best_result)) is_opt_ret = false;
        return is_opt_ret;
    }
    
    bool on_demand_completed = false;
    // on_demand_results[seq_i][seq_j] = {gap_in_seq_i, gap_in_seq_j}
    std::vector<std::vector<std::pair<bool, bool>>> on_demand_results(dim, std::vector<std::pair<bool, bool>>(dim, {false, false}));      

    DIRTYPE dir_cnt = (1 << dim);       // direction counts
    for (DIRTYPE dir = 1; dir < dir_cnt; ++dir) {
        NodeCoord nbr_crd;
        for (int i = 0; i < NN; ++i) nbr_crd[i] = 0;

        bool oob_flag = false;  // node out-of-bound flag
        for (int bit = 0; bit < dim; ++bit) {
            if ((dir>>bit) & 0x01) {
                nbr_crd[bit] = cur_crd_tb_info.crd[bit] + 1;                        // if i == 0b011, then dir_vec = [1, 1, 0]
                if (nbr_crd[bit] > seq_lens[bit]) oob_flag = true;      // don't go beyond the last cell
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
                    symbols[seq_idx] = sequences[seq_idx][seq_lens[seq_idx] - nbr_crd[seq_idx]];     // -1 is necessary
                } else {            // read the symbol normally
                    // crd = 0, symbol = -1
                    // crd = 1, symbol = 0
                    symbols[seq_idx] = sequences[seq_idx][nbr_crd[seq_idx] - 1];     // -1 is necessary
                }
            }
        }
        // 1. f-score = cur_g-score + match/mismatch/gap + h-score
        // 2. if f-score >= tmp_msa_gscore, push into open_list
        STYPE cur_gscore = local_gscore;        // cur_gscore = pre_gscore + the score of match/mismatch/gap, need check: use local gscore instead of global gscore 
        STYPE cur_fscore = 0;                   // cur_fscore = cur_gscore + cur_hscore

        // check all pairs of sequences
        for (int seq_i = 0; seq_i < dim - 1; ++seq_i) {
            for (int seq_j = seq_i + 1; seq_j < dim; ++seq_j) {
                bool gap_in_prev_seq[2] = {false, false};      // pass it to hscore function

                // compute and update the score of match/mismatch/gap
                if (nbr_crd[seq_i] != cur_crd_tb_info.crd[seq_i] && nbr_crd[seq_j] != cur_crd_tb_info.crd[seq_j]) {     // check the score table
                    cur_gscore += blosum_table.get_score_char(symbols[seq_i], symbols[seq_j]);
                    // gap_in_prev_seq[0] = false; gap_in_prev_seq[1] = false;      // default values
                } else {        // gap(s) in the two sequences
                    // There could be gaps in both sequences, so we need to check the gap length
                    bool gap_in_prev_seq_i = false, gap_in_prev_seq_j = false;
                    
                    // todo: matrix of on_demand_results is not used, optimize the memory allcoation
                    int gap_len_seq_i = cur_crd_tb_info.tb_info[seq_i], gap_len_seq_j = cur_crd_tb_info.tb_info[seq_j];
                    if (gap_len_seq_i > gap_len_seq_j) gap_in_prev_seq_i = true;
                    if (gap_len_seq_j > gap_len_seq_i) gap_in_prev_seq_j = true;
                    
                    // TODO: update the score computation
                    if (nbr_crd[seq_i] == cur_crd_tb_info.crd[seq_i] && nbr_crd[seq_j] == cur_crd_tb_info.crd[seq_j]) {     // two gaps here
                        gap_in_prev_seq[0] = gap_in_prev_seq_i; gap_in_prev_seq[1] = gap_in_prev_seq_j;
                    } else if (nbr_crd[seq_i] == cur_crd_tb_info.crd[seq_i]) {    // gap in seq i
                        cur_gscore -= (!gap_in_prev_seq_i || cur_crd_tb_info.crd[seq_i] == 0 && cur_crd_tb_info.crd[seq_j] == 0) * open_pen + ext_pen;
                        gap_in_prev_seq[0] = true;
                    } else {                                    // gap in seq j
                        cur_gscore -= (!gap_in_prev_seq_j || cur_crd_tb_info.crd[seq_i] == 0 && cur_crd_tb_info.crd[seq_j] == 0) * open_pen + ext_pen;
                        gap_in_prev_seq[1] = true;
                    }
                }
            }
        }
        cur_fscore += cur_gscore;

        // codeupdate
        typename AstarMultiIndex<NN>::template CoordTBInfo<AFFINE_TB_TYPE> nbr_crd_tb_info = {nbr_crd, std::vector<DIRTYPE>()};

        // find nbr node in ClosedList and compare the gscores
        nbr_crd_tb_info.tb_info.resize(dim);
        for (int seq_idx = 0; seq_idx < dim; ++seq_idx) {
            if ((dir >> seq_idx) & 0x01) 
                nbr_crd_tb_info.tb_info[seq_idx] = 0;       // non-gap in seq[seq_idx]
            else
                nbr_crd_tb_info.tb_info[seq_idx] = cur_crd_tb_info.tb_info[seq_idx] + 1;
        }

        // compute the hscore separately
        STYPE cur_hscore = score_upper_bound_recursive(dim, access_reversed_mats, nbr_crd, nbr_crd_tb_info.tb_info, recur_closed_list);
        cur_fscore += cur_hscore;

        // if cur_fscore is worse than tmp_msa_gscore, don't insert it into the open_list and closed_list
        if (cost_instead_of_score && (cur_fscore > tmp_msa_gscore + __FLT_EPSILON__) || 
            !cost_instead_of_score && (cur_fscore + __FLT_EPSILON__ < tmp_msa_gscore)) {
            workload_recorder.node_pruned_cnt_hscore += 1;
            // printf("cur_fscore - tmp_msa_gscore = %f\n", cur_fscore - tmp_msa_gscore);
            continue;
        }
        
        auto closed_list_gap_len_iter = target_closed_list->find(nbr_crd_tb_info.crd);
        if (closed_list_gap_len_iter != target_closed_list->end()) {
            auto gap_len_score_iter = closed_list_gap_len_iter->second.find(nbr_crd_tb_info.tb_info);
            if (gap_len_score_iter != closed_list_gap_len_iter->second.end()) {
                if (cost_instead_of_score && (cur_gscore >= gap_len_score_iter->second) || 
                    !cost_instead_of_score && (cur_gscore <= gap_len_score_iter->second)) continue;
            }
        }

        // Detect duplicated nodes in OpenList
        int nbr_bin = (int)floor(cur_gscore / bin_score_thres);
        // printf("bin_score_thres = %f, cur_gscore = %f, nbr_bin = %d\n", bin_score_thres, cur_gscore, nbr_bin);
        if (nbr_bin < 0) nbr_bin = 0;
        if (nbr_bin >= bin_cnt)  {
            if (cur_fscore > bin_score_thres * bin_cnt) {
                printf("------Error! Invalid bin number! ");
                printf("nbr_crd = [");
                for (int i = 0; i < dim; ++i) {
                    if (i == dim - 1) printf("%d", nbr_crd_tb_info.crd[i]);
                    else printf("%d, ", nbr_crd_tb_info.crd[i]);
                }
                printf("]. cur_gscore = %g, cur_fscore = %g, bin_score_thres = %g------\n", 
                    cur_gscore, cur_fscore, bin_score_thres);
            }
            nbr_bin = bin_cnt - 1;
        }

        bool found_in_higher_bin = false;    // First find the nbr_node in higher bins
        if (cost_instead_of_score) {    // from nbr_bin - 1 to 0
            for (int tmp_bin_idx = nbr_bin - 1; tmp_bin_idx >= 0; --tmp_bin_idx) {
                auto iterFound = (*target_index_by_coord)[tmp_bin_idx]->find(nbr_crd_tb_info);
                if (iterFound != (*target_index_by_coord)[tmp_bin_idx]->end()) {     // cur_gscore cannot be larger than the previous one
                    workload_recorder.node_pruned_cnt_gscore += 1;
                    found_in_higher_bin = true;
                    break;
                }
            }
        } else {                        // from nbr_bin + 1 to bin_cnt - 1
            for (int tmp_bin_idx = nbr_bin + 1; tmp_bin_idx < bin_cnt; ++tmp_bin_idx) {
                auto iterFound = (*target_index_by_coord)[tmp_bin_idx]->find(nbr_crd_tb_info);
                if (iterFound != (*target_index_by_coord)[tmp_bin_idx]->end()) {     // cur_gscore cannot be larger than the previous one
                    workload_recorder.node_pruned_cnt_gscore += 1;
                    found_in_higher_bin = true;
                    break;
                }
            }
        }
        if (found_in_higher_bin == false) {  // If the node is not in found, search bin[nbr_bin]
            auto iterFound = (*target_index_by_coord)[nbr_bin]->find(nbr_crd_tb_info);
            if (iterFound != (*target_index_by_coord)[nbr_bin]->end()) {
                if (cost_instead_of_score && (cur_gscore >= iterFound->gscore) || 
                    !cost_instead_of_score && (cur_gscore <= iterFound->gscore)) {      // if cur_gscore is worse than existing gscore, skip; else update gscore
                    workload_recorder.node_pruned_cnt_gscore += 1;      // prune this node
                } else {    // update gscore
                    // codeupdate
                    (*target_index_by_coord)[nbr_bin]->modify(iterFound, typename AstarMultiIndex<NN>::template ChangeGscore<AFFINE_TB_TYPE>(cur_gscore));
                    workload_recorder.node_pruned_cnt_gscore += 1;      // prune this node
                }
            } else {    // insert a new node
                // codeupdate
                typename AstarMultiIndex<NN>::template OpenListNode<AFFINE_TB_TYPE> nbr_node = {
                    nbr_crd_tb_info,
                    cur_fscore,
                    cur_gscore
                };
                (*target_index_by_fscore)[nbr_bin]->insert(nbr_node);

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
                }
            }
        } else {
            if (nbr_bin > 0) {
                auto iterFound = (*target_index_by_coord)[nbr_bin - 1]->find(nbr_crd_tb_info);
                if (iterFound != (*target_index_by_coord)[nbr_bin - 1]->end()) {
                    workload_recorder.node_pruned_cnt_gscore += 1;      // prune this node
                    (*target_index_by_coord)[nbr_bin - 1]->erase(iterFound);
                    workload_recorder.cur_open_list_cnt -= 1;
                }
            }
        }
    }

    bool is_opt_ret = true;
    STYPE global_best_result = get_highest_fscore();
    if (cost_instead_of_score && (tmp_msa_gscore > global_best_result)) is_opt_ret = false;
    else if (!cost_instead_of_score && (tmp_msa_gscore < global_best_result)) is_opt_ret = false;

    return is_opt_ret;
}



template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
inline void AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::compute_anytime_Astar() {
    const int dim = nn;
    NodeCoord cur_crd;
    for (int i = 0; i < NN; ++i) cur_crd[i] = 0;

    std::vector<DIRTYPE> cur_path(0);
    cur_path = std::vector<DIRTYPE>(dim, 0);
    // codeupdate
    typename AstarMultiIndex<NN>::template CoordTBInfo<AFFINE_TB_TYPE> cur_crd_tb_info = {cur_crd, cur_path};

    closed_list_gap_len[cur_crd_tb_info.crd][cur_crd_tb_info.tb_info] = 0.0;
    workload_recorder.closed_list_cnt += 1;      // reset to 1
    
    STYPE highest_score = 0;
    for (int seq_i = 0; seq_i < dim - 1; ++seq_i) {
        for (int seq_j = seq_i + 1; seq_j < dim; ++seq_j) {
            // compute PSA score with NW_2D
            int seq_idx[2] = {seq_i, seq_j}, seq_offset[2] = {0, 0}, cur_seq_lens[2] = {seq_lens[seq_i], seq_lens[seq_j]};
            bool gap_in_prev_seq[2] = {false, false};
            bool access_reversed_mats = true;
            STYPE score_upper_bound = score_upper_bound_2D(access_reversed_mats, seq_idx, seq_offset, gap_in_prev_seq, cur_seq_lens);
            // printf("seq_i = %d, seq_j = %d, sccore_upper_bound = %g\n", seq_i, seq_j, score_upper_bound);
            highest_score += score_upper_bound;
        }
    }

    printf("In compute_anytime_Astar, dim = %d, highest score (All-to-all PSA) = %f\n", dim, highest_score);

    bin_score_thres = highest_score / bin_cnt;
    // codeupdate
    typename AstarMultiIndex<NN>::template OpenListNode<AFFINE_TB_TYPE> source_node = {
        cur_crd_tb_info,
        0 + highest_score,
        0
    };
    index_by_fscore_arr[0]->insert(source_node);         // bin[0] because gscore = 0
    workload_recorder.open_list_cnt += 1;
    workload_recorder.cur_open_list_cnt += 1;

    std::vector<std::pair<double, double>> bin_op_cnt(bin_cnt, {0, 0});         // count insertions and erasions of each bin
    bin_op_cnt[0].first = 1;
    
    bool is_opt = false;        // if the current alignment result is optimal
    while (is_opt == false && time_limit_timer->check_time_limit() == false && !open_list_arr_empty()) {
        for (int bin_level = 0; bin_level < bin_cnt; ++bin_level) {     // perform Anytime Column Search iterations
            int beam = 0;
            while (is_opt == false && time_limit_timer->check_time_limit() == false && !open_list_arr_empty(bin_level) && beam < beam_width) {
                // printf("In beam search, bin_level = %d, total iter %g\n", bin_level, iter_cnt);
                if ((int)workload_recorder.iter_cnt % 1000000 == 0) {
                    printf("Total iter = %g, memory usage = %g bytes\n", workload_recorder.iter_cnt, (double) get_memory_usage());
                }
                workload_recorder.iter_cnt += 1;
                bool specific_bin = true;   // expand a node in bin[bin_level]
                bool reverse_seq = false;   // we are not computing recursive MSA
                is_opt = expand_node(dim, specific_bin, bin_level, reverse_seq, workload_recorder, bin_op_cnt);
                beam += 1;
            }
            // filter(V_ext): filter the nodes with lower_bound?
        }
        for (int astar_iter = 0; astar_iter < astar_iter_cnt; ++astar_iter) {       // perform normal A* Search iterations
            if (is_opt || time_limit_timer->check_time_limit() || open_list_arr_empty()) break;
            // printf("In astar search, total iter %g\n", iter_cnt);
            if ((int)workload_recorder.iter_cnt % 1000000 == 0) {
                printf("Total iter = %g, memory usage = %g bytes\n", workload_recorder.iter_cnt, (double) get_memory_usage());
            }
            workload_recorder.iter_cnt += 1;
            bool specific_bin = false;      // expand a node in bin_arr by searching one with the largest fscore
            bool reverse_seq = false;       // we are not computing recursive MSA
            is_opt = expand_node(dim, specific_bin, 0, reverse_seq, workload_recorder, bin_op_cnt);
        }

    }
    printf("Open list count = %g, current open list count = %g, closed list count = %g, iteration count = %g, nodes pruned by gscore = %g, nodes pruned by hscore = %g\n", 
        workload_recorder.open_list_cnt, workload_recorder.cur_open_list_cnt, workload_recorder.closed_list_cnt, workload_recorder.iter_cnt, workload_recorder.node_pruned_cnt_gscore, workload_recorder.node_pruned_cnt_hscore);
    printf("Beam search iter = %d, Astar search iter = %d, Bin count = %d, Insertions of each bin = [", beam_width, astar_iter_cnt, bin_cnt);
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
    printf("Memory usage = %g bytes\n", (double) get_memory_usage());
}




/**
 * @brief compute the alignment score upper bound with the specific method
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

    STYPE DP_val;
    int TB_val;
    pairwise_heuristic->get_score(access_reversed_mats, seq_idx[0], seq_idx[1], seq_offset[0], seq_offset[1], DP_val, TB_val);
    score_upper_bound = DP_val;
    // Deal with the TB_val. If both sides are gaps, and the gaps are in the same direction, score += gap_open_penalty
    bool gap_in_seq_i = ((TB_val >> 0) & 0x01) == 0, gap_in_seq_j = ((TB_val >> 1) & 0x01) == 0;

    if (gap_in_prev_seq[0] && gap_in_seq_i || gap_in_prev_seq[1] && gap_in_seq_j) 
        score_upper_bound += blosum_table.get_open_penalty();
    return score_upper_bound;
};

// store_gap_len = true
// write the result into msa_result and return alignment length
// print_path == true: only print the backtrack path in terminal; print_path == false: output the alignment
template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
inline int AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::backtrack_affine(char **&msa_result, typename AstarMultiIndex<NN>::template CoordTBInfo<AFFINE_TB_TYPE> crd_tb_info, bool print_path) {       // codeupdate
    NodeCoord seq_indices = crd_tb_info.crd;        // indices of the TB_mat, use --idx{i} when accessing sequences[][]

    int res_idx = 0;                                // counter of the alignment length
    bool all_zero_indices = true;
    for (int i = 0; i < nn; ++i)
        if (seq_indices[i] != 0) {
            all_zero_indices = false;
            break;
        }

    std::vector<DIRTYPE> gap_len = crd_tb_info.tb_info;

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
        STYPE max_gscore = INT_MIN;
        std::vector<DIRTYPE> max_gap_len;
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
                if (cost_instead_of_score && (iter->second < max_gscore) ||
                    !cost_instead_of_score && (iter->second > max_gscore)) {
                    max_gscore = iter->second;
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
inline STYPE AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::score_upper_bound_recursive(int dim, bool access_reversed_mats, NodeCoord cur_crd, std::vector<DIRTYPE> gap_len, 
    std::unordered_map<NodeCoord, ClosedListNodeGapLen, VectorIntHash<NN>>* recur_closed_list) {

    STYPE open_pen = blosum_table.get_open_penalty();
    STYPE res = 0;        // we first compute the recursive score (without the last sequence), then add the PSA scores of the last sequence
    
    // first find the crd in the recursive closed list
    NodeCoord reversed_crd;        // use dim - 1 because recur_closed_list contains (dim - 1)-D results
    for (int i = 0; i < NN; ++i) reversed_crd[i] = 0;

    for (int i = 0; i < dim - 1; ++i) reversed_crd[i] = seq_lens[i] - cur_crd[i];

    auto recur_closed_list_iter = recur_closed_list->find(reversed_crd);
    // compare the (dim-1)*(dim-1) distance matrix
    if (recur_closed_list_iter != recur_closed_list->end()) {
        
        // subtract the duplicated gap open penalty using the gap lengths
        for (auto node_iter = recur_closed_list_iter->second.begin(); node_iter != recur_closed_list_iter->second.end(); ++node_iter) {
            std::vector<DIRTYPE> recur_gap_len = node_iter->first;
            STYPE recur_gscore = node_iter->second;
            // the recursive closed list has no information of seq(dim - 1)
            for (int seq_i = 0; seq_i < dim - 2; ++seq_i) {
                for (int seq_j = seq_i + 1; seq_j < dim - 1; ++seq_j) {
                    DIRTYPE cur_gap_diff = gap_len[seq_i] - gap_len[seq_j];
                    DIRTYPE recur_gap_diff = recur_gap_len[seq_i] - recur_gap_len[seq_j];
                    // cur_gap_diff * recur_gap_diff > 0: Two gaps in either seq_i and seq_j
                    if (cur_gap_diff * recur_gap_diff > 0) 
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
    } else {
        for (int seq_i = 0; seq_i < dim - 1; ++seq_i) {
            for (int seq_j = seq_i + 1; seq_j < dim; ++seq_j) {
                // printf("seq_i = %d, seq_j = %d\n", seq_i, seq_j);
                int seq_idx[2] = {seq_i, seq_j}, seq_offset[2] = {cur_crd[seq_i], cur_crd[seq_j]}, cur_seq_lens[2] = {seq_lens[seq_i] - cur_crd[seq_i], seq_lens[seq_j] - cur_crd[seq_j]};
                bool gap_in_prev_seq[2] = {false, false};
                gap_in_prev_seq[0] = gap_len[seq_i] > gap_len[seq_j];
                gap_in_prev_seq[1] = gap_len[seq_i] < gap_len[seq_j];
                res += score_upper_bound_2D(access_reversed_mats, seq_idx, seq_offset, gap_in_prev_seq, cur_seq_lens);
            }
        }
    }

    return res;
}

// store closed list nodes with gap length
template <typename OLMultiIdx, typename OLFscore, typename OLCoord, int NN>
inline void AnytimeAstarSolver<OLMultiIdx, OLFscore, OLCoord, NN>::compute_recursive_Astar(int seq_cnt, bool reverse_seq) {
    // temporarily change the hyper-parameters: only use Astar in recursive MSA
    int prev_beam_width = beam_width, prev_bin_cnt = bin_cnt;
    beam_width = 0; bin_cnt = 1;         // todo: should we use beam search in lower dimension?

    /* --- Open list of nn-D MSA --- */
    // initialize open_list_arr, deallocation seems unnecessary
    open_list_arr.clear();
    index_by_fscore_arr.clear();        // todo: do we need to deallocate the pointers?
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

    std::vector<DIRTYPE> cur_path(0);
    cur_path = std::vector<DIRTYPE>(dim, 0);
    
    // codeupdate
    typename AstarMultiIndex<NN>::template CoordTBInfo<AFFINE_TB_TYPE> cur_crd_tb_info = {cur_crd, cur_path};

    // the previous recursive closed list results are in recursive_closed_lists[0]
    // swap the two pointers and clear recursive_closed_lists[0]
    auto tmp_ptr = recursive_closed_lists[1];
    recursive_closed_lists[1] = recursive_closed_lists[0];
    recursive_closed_lists[0] = tmp_ptr;
    recursive_closed_lists[0]->clear();

    (*recursive_closed_lists[0])[cur_crd_tb_info.crd][cur_crd_tb_info.tb_info] = 0.0;
    
    // codeupdate
    // bool gap_in_prev_seq[dim] = {false};
    bool gap_in_prev_seq[dim];
    for (int i = 0; i < dim; ++i) gap_in_prev_seq[i] = false;

    // if reverse_seq == false, access reversed_mats
    bool access_reversed_mats = !reverse_seq;
    STYPE highest_score = score_upper_bound_recursive(dim, access_reversed_mats, cur_crd, cur_path, recursive_closed_lists[1]);

    printf("In compute recursive Astar, dim = %d, highest score (recursive) = %f\n", dim, highest_score);

    bin_score_thres = highest_score / bin_cnt;
    // codeupdate
    typename AstarMultiIndex<NN>::template OpenListNode<AFFINE_TB_TYPE> source_node = {
        cur_crd_tb_info,
        0 + highest_score,
        0
    };
    index_by_fscore_arr[0]->insert(source_node);         // bin[0] because gscore = 0
    workload_recorder.open_list_cnt += 1;
    workload_recorder.cur_open_list_cnt += 1;

    std::vector<std::pair<double, double>> bin_op_cnt(bin_cnt, {0, 0});         // count insertions and erasions of each bin
    bin_op_cnt[0].first = 1;
    
    bool is_opt = false;        // if the current alignment result is optimal
    while (is_opt == false && time_limit_timer->check_time_limit() == false && !open_list_arr_empty()) {
        for (int astar_iter = 0; astar_iter < astar_iter_cnt; ++astar_iter) {       // perform normal A* Search iterations
            if (is_opt || time_limit_timer->check_time_limit() || open_list_arr_empty()) break;
            if ((int)workload_recorder.iter_cnt % 1000000 == 0) {
                printf("Total iter = %g, memory usage = %g bytes\n", workload_recorder.iter_cnt, (double) get_memory_usage());
            }
            workload_recorder.iter_cnt += 1;
            bool specific_bin = false;      // expand a node in bin_arr by searching one with the largest fscore
            is_opt = expand_node(dim, specific_bin, 0, reverse_seq, workload_recorder, bin_op_cnt);
        }
    }

    printf("%dD recursive MSA completed!\n", dim);
    printf("Total open list insertions = %g, total closed list insertions = %g, iteration count = %g, nodes pruned by gscore = %g, nodes pruned by hscore = %g\n", 
        workload_recorder.open_list_cnt, workload_recorder.closed_list_cnt, workload_recorder.iter_cnt, workload_recorder.node_pruned_cnt_gscore, workload_recorder.node_pruned_cnt_hscore);
    printf("Insertions of each bin = [");
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
    printf("Memory usage = %g bytes\n", (double) get_memory_usage());

    // restore the hyper-parameters
    beam_width = prev_beam_width; bin_cnt = prev_bin_cnt;
}

#endif // ANYTIMEASTAR_HPP