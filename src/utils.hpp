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

    std::size_t operator() (const std::vector<DIRTYPE>& v) const {
        std::size_t seed = 0;
        for (DIRTYPE i : v) {
            seed ^= std::hash<DIRTYPE>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }

    std::size_t operator() (const std::array<int, NN>& arr) const {
        std::size_t seed = 0;
        for (int i : arr) {
            seed ^= std::hash<int>{}(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
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
    STYPE blosum_scores[CHAR_NUM + 1][CHAR_NUM + 1];
    void load_blosum_scores(const char* file_dir);
    std::map<std::pair<int, int>, std::vector<scorepair>> problematic_pairs;    // {outer pair, {score gain, inner pair}}
public:
    static bool score_cmp(scorepair p0, scorepair p1) {return p0.first > p1.first;}
    ScoreTable(STYPE _open_penalty, STYPE _ext_penalty, const char* file_dir) {
        open_penalty = _open_penalty;
        ext_penalty = _ext_penalty;
        load_blosum_scores(file_dir);
        printf("ScoreTable init complete!\n");
    }
    STYPE get_open_penalty() {return open_penalty;}
    STYPE get_ext_penalty() {return ext_penalty;}
    const int *get_blosum_map() {return blosum_map;}
    std::map<std::pair<int, int>, std::vector<scorepair>> get_problematic_pairs() {return problematic_pairs;}

    // input two characters and get the blosum score
    STYPE get_score_char(int c0, int c1) {return blosum_scores[blosum_map[c0]][blosum_map[c1]];}
    // input two mapped symbols and get the blosum score
    STYPE get_score_char_mapped(int mapped_c0, int mapped_c1) {return blosum_scores[mapped_c0][mapped_c1];}
    STYPE calc_score(char const* const*msa, const int seq_cnt, const int row_len);
};


inline void ScoreTable::load_blosum_scores(const char* file_dir) {
    for (int i = 0; i < 128; ++i)
        blosum_map[i] = -1;
	std::ifstream file(file_dir);
	if (file.good()) {
        int line_num = 0;
        STYPE digit_item;
        
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


// compute the SP-score of the input MSA, the length of each row should be n*l
inline STYPE ScoreTable::calc_score(char const* const*msa, int seq_cnt, int row_len) {
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
    FileWriter(std::string _res_out_file_dir, int _seq_cnt, std::string sequences_dir, char **sequences, int *seq_lens) {
        RES_OUT = _res_out_file_dir;
        seq_cnt = _seq_cnt;
        res_out.open(RES_OUT, std::fstream::out | std::fstream::app);
        if (!res_out.is_open()) {
            printf("Anytime A* intermediate result file open failed!\n");
        }

        res_out << "\n\nInput file is " << sequences_dir << std::endl;
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

    void write(double msa_time, char** msa_result, int msa_alignment_length, STYPE msa_score, double iter_cnt) {
        res_out << "\nMSA number " << ++result_cnt << ", alignment length = " << msa_alignment_length << 
            ", score = " << std::fixed << std::setprecision(4) << msa_score << ", execution time = " << std::fixed << std::setprecision(5) << msa_time << " s, iteration count = " << (int) iter_cnt << "\n";
        for (int i = 0; i < seq_cnt; ++i) {
            res_out << "Sequence " << i << ": ";
            res_out.write(msa_result[i], msa_alignment_length);
            res_out << std::endl;
            // res_out << "Sequence " << i << ": " << msa_result[i] << std::endl;
        }
    }

    void write(std::string info) {
        res_out << info << std::endl;
    }
};

#endif // UTILS_HPP