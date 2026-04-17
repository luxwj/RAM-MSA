#ifndef UTILS_HPP
#define UTILS_HPP

#include <unistd.h>     // sysconf
#include <fstream>
#include <iomanip>	    // setprecision
#include <sstream>
#include <algorithm>    // sort
#include <vector>
#include <unordered_set>          // resolve the score table constraints
#include <set>          // resolve the score table constraints
#include <map>          // resolve the score table constraints

#include <chrono>       // timer

#include "parameters.hpp"

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))

template <int NN>
struct VectorIntHash {
    std::size_t operator() (const std::vector<int>& v) const {
        std::size_t seed = 0;
        for (int i : v) {
            seed ^= std::hash<int>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }

    std::size_t operator() (const std::array<DIRTYPE, NN>& v) const {
        std::size_t seed = 0;
        for (DIRTYPE i : v) {
            seed ^= std::hash<DIRTYPE>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }

    std::size_t operator() (const std::array<CRDTYPE, NN>& arr) const {
        std::size_t seed = 0;
        for (CRDTYPE i : arr) {
            seed ^= std::hash<CRDTYPE>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }

    std::size_t operator() (const DIRTYPE& val) const {
        std::size_t seed = 0;
        seed ^= std::hash<DIRTYPE>{}(val) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};

class ScoreTable {
public:
    typedef std::pair<STYPE, std::pair<int, int>> scorepair;   // pair score
    
private:
    STYPE open_penalty;
    STYPE ext_penalty;
    int blosum_map[128];        // map characters onto numbers, blosum_map[A] = 0, blosum_map[B] = 1, ...
    // CHAR_NUM + 1 for symbol 'X', which has the average cost of all other symbol pairs
    STYPE **blosum_scores;

    STYPE lowest_score_in_table;        // the lowest score/cost in the substitution matrix

    /* --- TC score & SP accuracy variables --- */
    struct TCHash {
        std::size_t operator() (const std::vector<int>& v) const {
            std::size_t seed = 0;
            for (int i : v) {
                seed ^= std::hash<int>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };

    std::unordered_set<std::vector<int>, TCHash> us_tc;     // each element is the symbol indices in a column

    std::unordered_set<std::vector<int>, TCHash> us_sp;     // each element is {seq_i, seq_j, symbol_idx_i, symbol_idx_j}

    bool is_ref_msa_set = false;
    int ref_seq_cnt;
    int ref_alignment_len;
    /* --- TC score & SP accuracy variables --- */

    std::vector<scorepair> blosum_scores_sorted;    // {score, {symbol0, symbol1}}
    void load_blosum_scores(const char* file_dir);
    std::map<std::pair<int, int>, std::vector<scorepair>> problematic_pairs;    // {outer pair, {score gain, inner pair}}
    void sort_blosum_scores();
public:
    static bool score_cmp(scorepair p0, scorepair p1) {return p0.first > p1.first;}
    ScoreTable(STYPE _open_penalty, STYPE _ext_penalty, const char* score_table_dir) {
        open_penalty = _open_penalty;
        ext_penalty = _ext_penalty;
        blosum_scores = new STYPE *[CHAR_NUM + 1];
        for (int i = 0; i < CHAR_NUM + 1; ++i)
            blosum_scores[i] = new STYPE[CHAR_NUM + 1];

        load_blosum_scores(score_table_dir);
        sort_blosum_scores();
        printf("ScoreTable init complete!\n");
    }

    ~ScoreTable() {
        for (int i = 0; i < CHAR_NUM + 1; ++i)
            delete [] blosum_scores[i];
        delete [] blosum_scores;
    }

    bool is_ref_MSA_set() const {return is_ref_msa_set;}
    void set_reference_msa(char **ref_msa_results, int ref_seq_cnt, int ref_alignment_len);
    STYPE compute_total_column_score(char **test_msa_results, int test_seq_cnt, int test_alignment_len) const;
    STYPE compute_sum_of_pairs_accuracy(char **test_msa_results, int test_seq_cnt, int test_alignment_len) const;

    STYPE get_open_penalty() const {return open_penalty;}
    STYPE get_ext_penalty() const {return ext_penalty;}
    const int *get_blosum_map() const {return blosum_map;}
    STYPE get_lowest_score_in_table() const {return lowest_score_in_table;}

    const STYPE *const *get_blosum_scores() const {return blosum_scores;}

    std::vector<scorepair> get_blosum_scores_sorted() const {return blosum_scores_sorted;}
    std::map<std::pair<int, int>, std::vector<scorepair>> get_problematic_pairs() const {return problematic_pairs;}

    // input two characters and get the blosum score
    STYPE get_score_char(int c0, int c1) const {
        int idx0 = blosum_map[c0], idx1 = blosum_map[c1];
        if (idx0 < 0 || idx0 > CHAR_NUM || idx1 < 0 || idx1 > CHAR_NUM) 
            printf("Index out of bound in get_score_char()!\n");
        return blosum_scores[blosum_map[c0]][blosum_map[c1]];
    }
    // input two mapped symbols and get the blosum score
    STYPE get_score_char_mapped(int mapped_c0, int mapped_c1) const {
        if (mapped_c0 < 0 || mapped_c0 > CHAR_NUM || mapped_c1 < 0 || mapped_c1 > CHAR_NUM) 
            printf("Index out of bound in get_score_char_mapped()!\n");
        return blosum_scores[mapped_c0][mapped_c1];
    }
    STYPE calc_score_v1(char const* const*msa, const int seq_cnt, const int row_len) const;
    STYPE calc_score(char const* const*msa, const int seq_cnt, const int row_len) const;
};

// build the unordered set for TC scores so there is no need to store the ref msa
inline void ScoreTable::set_reference_msa(char **ref_msa_results, int _ref_seq_cnt, int _ref_alignment_len) {
    us_tc.clear(); us_sp.clear();
    ref_seq_cnt = _ref_seq_cnt; ref_alignment_len = _ref_alignment_len;

    // construct unordered_set for TC score
    // if ref_msa_results[ref_row][ref_col] != GAP, cur_col[ref_row] = ref_indices[ref_row]++;
    std::vector<int> ref_indices(ref_seq_cnt, 0);
    for (int ref_col = 0; ref_col < ref_alignment_len; ++ref_col) {
        std::vector<int> cur_col(ref_seq_cnt);
        for (int ref_row = 0; ref_row < ref_seq_cnt; ++ref_row) {
            int symbol = (int) ref_msa_results[ref_row][ref_col];
            if (symbol != GAP)
                cur_col[ref_row] = ref_indices[ref_row]++;
            else
                cur_col[ref_row] = -1;
        }
        us_tc.insert(cur_col);
    }

    // construct unordered_set for SP accuracy
    std::vector<int> ref_indices_sp(ref_seq_cnt, 0);
    std::vector<bool> is_gap_sp(ref_seq_cnt, false);
    for (int ref_col = 0; ref_col < ref_alignment_len; ++ref_col) {
        // update ref_indices
        for (int ref_row = 0; ref_row < ref_seq_cnt; ++ref_row) {
            int symbol = (int) ref_msa_results[ref_row][ref_col];
            if (symbol != GAP) {
                ref_indices_sp[ref_row] += 1;
                is_gap_sp[ref_row] = false;
            } else
                is_gap_sp[ref_row] = true;
        }
        // insert elements into us
        for (int seq_i = 0; seq_i < ref_seq_cnt - 1; ++seq_i) {
            for (int seq_j = seq_i + 1; seq_j < ref_seq_cnt; ++seq_j) {
                if (is_gap_sp[seq_i] || is_gap_sp[seq_j]) continue;

                std::vector<int> cur_col_sp;
                cur_col_sp.push_back(seq_i); 
                cur_col_sp.push_back(seq_j);
                cur_col_sp.push_back(ref_indices_sp[seq_i]);
                cur_col_sp.push_back(ref_indices_sp[seq_j]);

                us_sp.insert(cur_col_sp);
            }
        }
    }

    is_ref_msa_set = true;
    printf("Set the reference MSA!\n");
}

inline STYPE ScoreTable::compute_total_column_score(char **test_msa_results, int test_seq_cnt, int test_alignment_len) const {
    if (is_ref_msa_set == false) {
        printf("Haven't set the reference MSA in ScoreTable!\n");
        return -1;
    }

    if (ref_seq_cnt != test_seq_cnt) {
        printf("Different ref_seq_cnt and test_seq_cnt in ScoreTable::compute_total_column_score()!\n");
        return -1;
    }
    STYPE score = 0;
    std::vector<int> test_indices(test_seq_cnt, 0);
    for (int test_col = 0; test_col < test_alignment_len; ++test_col) {
        std::vector<int> cur_col(test_seq_cnt);
        for (int test_row = 0; test_row < test_seq_cnt; ++test_row) {
            int symbol = (int) test_msa_results[test_row][test_col];
            if (symbol != GAP)
                cur_col[test_row] = test_indices[test_row]++;
            else
                cur_col[test_row] = -1;
        }
        if (us_tc.find(cur_col) != us_tc.end()) {
            score += 1;
        } else {
            // nothing
        }
    }
    score = score / ref_alignment_len;
    return score;

}

inline STYPE ScoreTable::compute_sum_of_pairs_accuracy(char **test_msa_results, int test_seq_cnt, int test_alignment_len) const {
    if (is_ref_msa_set == false) {
        printf("Haven't set the reference MSA in ScoreTable!\n");
        return -1;
    }

    if (ref_seq_cnt != test_seq_cnt) {
        printf("Different ref_seq_cnt and test_seq_cnt in ScoreTable::compute_sum_of_pairs_accuracy()!\n");
        return -1;
    }
    STYPE score = 0;

    std::vector<int> test_indices(test_seq_cnt, 0);
    std::vector<bool> is_gap(test_seq_cnt, false);
    for (int test_col = 0; test_col < test_alignment_len; ++test_col) {
        // update test_indices
        for (int test_row = 0; test_row < test_seq_cnt; ++test_row) {
            int symbol = (int) test_msa_results[test_row][test_col];
            if (symbol != GAP) {
                test_indices[test_row] += 1;
                is_gap[test_row] = false;
            } else
                is_gap[test_row] = true;
        }

        for (int seq_i = 0; seq_i < test_seq_cnt - 1; ++seq_i) {
            for (int seq_j = seq_i + 1; seq_j < test_seq_cnt; ++seq_j) {
                if (is_gap[seq_i] || is_gap[seq_j]) continue;

                std::vector<int> cur_col;
                cur_col.push_back(seq_i); 
                cur_col.push_back(seq_j);
                cur_col.push_back(test_indices[seq_i]);
                cur_col.push_back(test_indices[seq_j]);

                if (us_sp.find(cur_col) != us_sp.end()) {
                    score += 1;
                } else {
                    // nothing
                }
            }
        }
    }
    score = score / us_sp.size();
    return score;
}


inline void ScoreTable::load_blosum_scores(const char* file_dir) {
    for (int i = 0; i < 128; ++i)
        blosum_map[i] = -1;
	std::ifstream file(file_dir);
	if (file.good()) {
        int line_num = 0;
        STYPE digit_item;
        lowest_score_in_table = 65535;      // sentinel value
        
        STYPE total_score = 0;      // computing the average cost/score of X

        std::string line, item;
        while (getline(file, line)) {
            int col_num = 0;
			std::stringstream ss(line);
            while (ss >> item) {
                if (line_num == 0)
                    blosum_map[(int)item[0]] = col_num;
                else if (col_num > 0) {
                    digit_item = stod(item);
                    blosum_scores[line_num - 1][col_num - 1] = digit_item;
                    total_score += digit_item;      // for symbol X
                    if (digit_item < lowest_score_in_table) lowest_score_in_table = digit_item;
                }
                ++col_num;
            }
            // write X
            if (line_num == 0) blosum_map[(int)'X'] = col_num;
            ++line_num;
        }
        STYPE average_score = total_score / (CHAR_NUM * CHAR_NUM);
        for (int row = 0; row < CHAR_NUM + 1; ++row) blosum_scores[row][CHAR_NUM] = average_score;
        for (int col = 0; col < CHAR_NUM + 1; ++col) blosum_scores[CHAR_NUM][col] = average_score;
        printf("In load_blosum_scores(), total score = %f, average score = score of 'X' = %f\n", total_score, average_score);

	} else {
        printf("In load_blosum_scores(), found issues with the score table file!\n");
    }
}

// Sort the scores in descending order and find the problematic pairs
// problem_map<outer pair, sorted inner pairs>: the inner pairs are sorted based on the score gain
inline void ScoreTable::sort_blosum_scores() {
    blosum_scores_sorted.reserve((CHAR_NUM + 1) * (CHAR_NUM + 1));
    for (int symbol0 = 0; symbol0 < (CHAR_NUM + 1); ++symbol0)
        for (int symbol1 = 0; symbol1 < (CHAR_NUM + 1); ++symbol1) {
            STYPE score = blosum_scores[symbol0][symbol1];
            blosum_scores_sorted.push_back({score, {symbol0, symbol1}});
        }
    std::sort(blosum_scores_sorted.begin(), blosum_scores_sorted.end(), score_cmp);

    // Resolve the score constaints
    // Descending order
    // 1.   a   b   2.  a   a
    //      a   d       b   a
    //      c   b       a   c
    //      c   d       b   c
    // In case 1, score(a, b) + score(c, d) >= score(a, d) + score(c, b)
    // In case 2, score(a, a) + score(b, c) >= score(b, a) + score(a, c)
    
    // How to resolve: Similar to two sum problem, create some sets with four symbols
    // The outer loop traverses from the first element (base pair) in the sorted score table, and the inner loop verifies the remaining elements
    // Assumption: a match has a higher score than any mismatch/gap
    // If so, we dont need to consider the following cases
    // 3.   a   b   4.  b   c
    //      a   a       a   a   (match should has the highest score)
    //      c   b       a   c
    //      c   a       b   a
    // Consequently, the inner loop starts from the end and only inserts the pairs, which have no common symbols with the base pair, into the map

    for (auto outer_iter = blosum_scores_sorted.begin(); outer_iter + 1 != blosum_scores_sorted.end(); ++outer_iter) {
        int symbol0 = outer_iter->second.first, symbol1 = outer_iter->second.second;
        STYPE outer_score = outer_iter->first;
        std::vector<int> outer_pair = {symbol0, symbol1};
        std::map<std::vector<int>, std::pair<STYPE, bool>> score_map;   // key: four symbols, value: {score, once_accessed}
        std::map<std::pair<int, int>, std::vector<std::vector<int>>> pair_set_map;   // map the second and the third pairs to a set of four symbols, like mapping {a, d} to {a, b, c, d}
        // iter from the end to the beginning
        for (auto inner_iter = blosum_scores_sorted.end() - 1; inner_iter != outer_iter; --inner_iter) {
            int symbol2 = inner_iter->second.first, symbol3 = inner_iter->second.second;
            STYPE inner_score = inner_iter->first;
            if (symbol2 != symbol0 && symbol2 != symbol1 && symbol3 != symbol0 && symbol3 != symbol1) {     // no common symbols, insert the set into the score map
                std::vector<int> cur_set = outer_pair;
                cur_set.push_back(symbol2); cur_set.push_back(symbol3);
                score_map[cur_set] = {outer_score + inner_score, false};
                pair_set_map[{symbol0, symbol3}].push_back(cur_set);
                pair_set_map[{symbol2, symbol1}].push_back(cur_set);
            } else {
                std::pair<int, int> inner_pair = {symbol2, symbol3};
                auto found_iter = pair_set_map.find(inner_pair);
                if (found_iter != pair_set_map.end()) {
                    for (auto cur_set : found_iter->second) {
                        STYPE new_score = score_map[cur_set].first - inner_score;
                        if (new_score < 0 && score_map[cur_set].second == true) {
                            // push a new element into the problematic score vector, score gain = -new_score
                            problematic_pairs[{symbol0, symbol1}].push_back({-new_score, {cur_set[2], cur_set[3]}});
                        }
                        score_map[cur_set] = {new_score, true};
                    }
                }
                
            }
        }
    }
    for (auto iter = problematic_pairs.begin(); iter != problematic_pairs.end(); ++iter) {
        std::sort(iter->second.begin(), iter->second.end(), score_cmp);
    }
}

// compute the SP-score of the input MSA, the length of each row should be n*l
inline STYPE ScoreTable::calc_score_v1(char const* const*msa, int seq_cnt, int row_len) const {
    if (msa == nullptr || !msa[0]) {
        printf("Invalid array! Null detected!\n");
        return 0;
    }
    STYPE score = 0;
    for (int i = 0; i < seq_cnt - 1; ++i) {
        for (int j = i + 1; j < seq_cnt; ++j) {
            int gap_len1 = 0, gap_len2 = 0;     // current gap lengths in seq1 and seq2
            // loop over all row_len pairs of symbols
            for (int k = 0; k < row_len; ++k) {
                int c1 = msa[i][k], c2 = msa[j][k];
                if (c1 == GAP && c2 == GAP) {
                    gap_len1 += 1; gap_len2 += 1;
                } else if (c1 == GAP) {
                    // Case 1: gap_len1 > gap_len2
                    // simply gap_len1++
                    // Case 2: gap_len1 = gap_len2
                    // simply gap_len1++
                    // Case 3: gap_len1 < gap_len2 (non-gap in seq2)
                    // calc the gap penalty of seq2 and reset the gap_lens, gap_len1 = 1, gap_len2 = 0
                    
                    if (gap_len1 >= gap_len2) {
                        gap_len1 += 1;
                    } else {
                        score -= open_penalty + gap_len2 * ext_penalty;
                        printf("i = %d, j = %d, k = %d, gap_len1 = %d, gap_len2 = %d, compute gap2, score = %g\n", i, j, k, gap_len1, gap_len2, score);
                        gap_len1 = 1;
                    }
                    gap_len2 = 0;

                } else if (c2 == GAP) {
                    // Case 1: gap_len1 < gap_len2
                    // simply gap_len2++
                    // Case 2: gap_len1 = gap_len2
                    // simply gap_len2++
                    // Case 3: gap_len1 > gap_len2 (non-gap in seq1)
                    // calc the gap penalty of seq1 and reset the gap_lens, gap_len2 = 1, gap_len1 = 0
                    
                    if (gap_len1 <= gap_len2) {
                        gap_len2 += 1;
                    } else {
                        score -= open_penalty + gap_len1 * ext_penalty;
                        printf("i = %d, j = %d, k = %d, gap_len1 = %d, gap_len2 = %d, compute gap1, score = %g\n", i, j, k, gap_len1, gap_len2, score);
                        gap_len2 = 1;
                    }
                    gap_len1 = 0;

                } else {
                    score += blosum_scores[blosum_map[c1]][blosum_map[c2]];
                    // reset gap_len1 and gap_len2
                    int gap_len_max = gap_len1;
                    if (gap_len2 > gap_len1) gap_len_max = gap_len2;
                    if (gap_len_max > 0) {
                        score -= open_penalty + gap_len_max * ext_penalty;
                        printf("i = %d, j = %d, k = %d, gap_len1 = %d, gap_len2 = %d, compute gap%d, score = %g\n", i, j, k, gap_len1, gap_len2, (gap_len1 > gap_len2) ? 1 : 2, score);
                        gap_len1 = 0;
                        gap_len2 = 0;
                    }
                }
            }
            // reset gap_len1 and gap_len2
            int gap_len_max = gap_len1;
            if (gap_len2 > gap_len1) gap_len_max = gap_len2;
            if (gap_len_max > 0) {
                score -= open_penalty + gap_len_max * ext_penalty;
                printf("i = %d, j = %d, last gaps, gap_len1 = %d, gap_len2 = %d, compute gap%d, score = %g\n", i, j, gap_len1, gap_len2, (gap_len1 > gap_len2) ? 1 : 2, score);
                gap_len1 = 0;
                gap_len2 = 0;
            }
        }
    }
    return score;
}


// compute the SP-score of the input MSA, the length of each row should be n*l
inline STYPE ScoreTable::calc_score(char const* const*msa, int seq_cnt, int row_len) const {
    if (msa == nullptr || !msa[0]) {
        printf("Invalid array! Null detected!\n");
        return 0;
    }
    STYPE score = 0;
    for (int i = 0; i < seq_cnt - 1; ++i) {
        for (int j = i + 1; j < seq_cnt; ++j) {
            int gap_len1 = 0, gap_len2 = 0;     // current gap lengths in seq1 and seq2
            // loop over all row_len pairs of symbols
            for (int k = 0; k < row_len; ++k) {
                int c1 = msa[i][k], c2 = msa[j][k];
                if (c1 == GAP && c2 == GAP) {     // both gaps, skip
                    continue;
                } else if (c1 == GAP) {
                    gap_len1 += 1;
                    // reset gap_len2
                    if (gap_len2 > 0) {
                        score -= open_penalty + gap_len2 * ext_penalty;
                        gap_len2 = 0;
                    }
                } else if (c2 == GAP) {
                    gap_len2 += 1;
                    // reset gap_len1
                    if (gap_len1 > 0) {
                        score -= open_penalty + gap_len1 * ext_penalty;
                        gap_len1 = 0;
                    }
                } else {
                    score += blosum_scores[blosum_map[c1]][blosum_map[c2]];
                    // reset gap_len1 and gap_len2
                    if (gap_len1 > 0) {
                        score -= open_penalty + gap_len1 * ext_penalty;
                        gap_len1 = 0;
                    }
                    if (gap_len2 > 0) {
                        score -= open_penalty + gap_len2 * ext_penalty;
                        gap_len2 = 0;
                    }
                }
            }
            // reset gap_len1 and gap_len2
            if (gap_len1 > 0) {
                score -= open_penalty + gap_len1 * ext_penalty;
                gap_len1 = 0;
            }
            if (gap_len2 > 0) {
                score -= open_penalty + gap_len2 * ext_penalty;
                gap_len2 = 0;
            }
        }
    }
    return score;
}

class Timer {
private:
    std::chrono::time_point<std::chrono::system_clock> start_;
    std::chrono::time_point<std::chrono::system_clock> prev_;
    double time_limit;
public:
    Timer(): start_(std::chrono::system_clock::time_point::min()), prev_(std::chrono::system_clock::time_point::min()) {}
    // Timer(double _time_limit): start_{std::chrono::system_clock::time_point::min()}, prev_{std::chrono::system_clock::time_point::min()}, time_limit{_time_limit} {}

    void start() { prev_ = start_ = std::chrono::system_clock::now(); }

    void start(double _time_limit) {
        time_limit = _time_limit;
        prev_ = start_ = std::chrono::system_clock::now(); 
    }

    // return the lap time in s (from n-th elapsed to (n-1)-th elapsed)
    double elapsed(bool print_time) {
        auto now = std::chrono::system_clock::now();
        auto total = std::chrono::duration_cast<std::chrono::microseconds>(now - start_);
        auto lap = std::chrono::duration_cast<std::chrono::microseconds>(now - prev_);
        if (print_time) printf("Total time: %g s. Lap time: %g s.\n", total.count() / 1e6, lap.count() / 1e6);
        prev_ = now;
        return lap.count() / 1e6;
    }

    // time in second
    double get_time_from_start() {
        auto now = std::chrono::system_clock::now();
        auto total = std::chrono::duration_cast<std::chrono::microseconds>(now - start_);
        return total.count() / 1e6;
    }

    bool check_time_limit() {
        if (time_limit < 0) return false;
        auto now = std::chrono::system_clock::now();
        auto total = std::chrono::duration_cast<std::chrono::microseconds>(now - start_);
        return total.count() / 1e6 > time_limit;
    }
};

// get available virtual memory
static size_t get_available_memory() {
    std::ifstream mem_info("/proc/meminfo");
    std::string key;
    size_t value;
    std::string unit;
    size_t mem = 0;

    while (mem_info >> key >> value >> unit) {
        if (key == "MemAvailable:") mem += value;
        if (key == "SwapFree:") mem += value;
    }
    return mem * 1024;  // KB to B
}

// track the memory consumption on Linux (bytes)
static size_t get_memory_usage() {
    std::ifstream statm("/proc/self/statm");
    size_t size, resident, share, text, lib, data, dt;
    statm >> size >> resident >> share >> text >> lib >> data >> dt;
    return resident * sysconf(_SC_PAGESIZE);
}

// remove all empty columns and print the MSA
static void show_result(char **msa, STYPE score, int seq_cnt, int alignment_len) {
    printf("The final MSA result (score = %.4f):\n", score);
    for (int i = 0; i < seq_cnt; ++i) {
        printf("Sequence %d: ", i);
        for (int j = 0; j < alignment_len; ++j) {
            printf("%c", msa[i][j]);
        }
        printf("\n");
    }
}

class FileWriter {
private:
    std::string RES_OUT;
    std::ofstream res_out;
    int result_cnt;
    int seq_cnt;

public:
    // print the file dir and the input sequences in constructor
    FileWriter(std::string _res_out_file_dir, int _seq_cnt, std::string sequences_dir, char **sequences, int *seq_lens, STYPE pen_o, STYPE pen_e) {
        RES_OUT = _res_out_file_dir;
        seq_cnt = _seq_cnt;
        res_out.open(RES_OUT, std::fstream::out | std::fstream::app);
        if (!res_out.is_open()) {
            printf("Failed to open Anytime A* intermediate result file!\n");
        }

        res_out << "\n\nInput file is " << sequences_dir << ", gap open penalty is " << std::fixed << std::setprecision(4) << pen_o << ", gap extension penalty is " << pen_e << std::endl;
        for (int i = 0; i < seq_cnt; ++i) {
            res_out << "Input sequence " << i << ": ";
            res_out.write(sequences[i], seq_lens[i]);
            res_out << std::endl;
        }

        result_cnt = 0;
        res_out << "\nOutput intermediate MSA results!\n";
    }

    ~FileWriter() {
        res_out << "Output complete!\n\n\n";
        res_out.close();
    }

    std::string get_target_file() { return RES_OUT; }

    void write(double msa_time, char** msa_result, int msa_alignment_length, STYPE sp_score, double iter_cnt) {
        res_out << "\nMSA number " << ++result_cnt << ", alignment length = " << msa_alignment_length << 
            ", SP score = " << std::fixed << std::setprecision(4) << sp_score << ", execution time = " << std::fixed << std::setprecision(5) << msa_time << " s, iteration count = " << (int) iter_cnt << "\n";
        for (int i = 0; i < seq_cnt; ++i) {
            res_out << "Sequence " << i << ": ";
            res_out.write(msa_result[i], msa_alignment_length);
            res_out << std::endl;
            // res_out << "Sequence " << i << ": " << msa_result[i] << std::endl;
        }
    }

    void write_TC_score(STYPE tc_score) {
        res_out << "TC score = " << std::fixed << std::setprecision(4) << tc_score << std::endl;
    }

    void write_SP_accuracy(STYPE sp_accuracy) {
        res_out << "SP accuracy = " << std::fixed << std::setprecision(4) << sp_accuracy << std::endl;
    }

    void write(std::string info) {
        res_out << info << std::endl;
    }
};

#endif // UTILS_HPP