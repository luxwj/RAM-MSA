#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <vector>

#include "utils.hpp"        // contains the ScoreTable class
#include "parameters.hpp"
#include "NW/NW.hpp"
#include "Astar/AnytimeAstar.hpp"
#include "Astar/AnytimeAstarLinearGap.hpp"
#include "Astar/AnytimeAstarTemplateHelper.hpp"

const char* blosum_scores_dir = "../score_tables/blosum_double";            // double type
const char* pam250_scores_dir = "../score_tables/pam250_int";               // int type

std::string input_score_table;

std::string sequences_dir = "../data/1ac5.fasta";              // n = 5
// std::string sequences_dir = "../data/2ack.fasta";              // n = 5

const char* anytime_astar_result_dir = "../anytime_results.txt"; // intermediate results

int seq_cnt = 3, seq_max_len = 0;

double open_penalty = 1.53, ext_penalty = 0.0;      // {1.53, 0.00} in MAFFT
bool cost_instead_of_score = false;                 // set to true when using PAM250
double memory_limit_ratio = 0.8;      // default memory limit ratio is 0.8, disable memory-bound strategy when this value is < 0 or > 1
double anytime_time_limit = 10000000;          // time limit in second

char **sequences;
int *sequence_lengths;      // allocated in load_sequences()
char **anytime_Astar_msa_result;

// initialize the results arrays with seq_max_len
void init_res_arr() {
    anytime_Astar_msa_result = new char*[seq_cnt];
    for (int i = 0; i < seq_cnt; ++i)
        anytime_Astar_msa_result[i] = new char[seq_cnt * seq_max_len];

    for (int i = 0; i < seq_cnt; ++i)
        for (int j = 0; j < seq_cnt * seq_max_len; ++j) {
            anytime_Astar_msa_result[i][j] = GAP;
        }
}

void dealloc() {
    for (int i = 0; i < seq_cnt; ++i) {
        delete [] sequences[i];
        delete [] anytime_Astar_msa_result[i];
    }
    delete [] sequences;
    delete [] anytime_Astar_msa_result;
}

/**
 * @brief load a file in FASTA format: Each sequence begins with a '>'
 * New lines are acceptable within a single sequence
 * ';' represents a comment line
 * @param file_dir location of the FASTA file
 * @param target array to store the sequences
 * @param set_parameters determines whether this call sets seq_cnt and seq_max_len
 * @return seq_cnt and seq_max_len
 */
int* load_sequences(std::string file_dir, char **&target, bool set_parameters) {
    static int ret[2] = {0, 0};
    std::vector<int> seq_len;       // store the sequence lengths when loading the sequences
    std::vector<std::string> sequence_lines;
	std::ifstream file(file_dir);
	if (file.good()) {
        int sequence_count = 0, sequence_max_length = 0;
        std::string cur_seq = "", line;
        while (getline(file, line)) {
            if (line[0] == '>') {
                if (cur_seq != "") {
                    sequence_lines.push_back(cur_seq);
                    seq_len.push_back(cur_seq.length());
                    if (cur_seq.length() > sequence_max_length)
                        sequence_max_length = cur_seq.length();
                    cur_seq = "";
                }
                ++sequence_count;
            }
            else if (line.length() != 0 && line[0] != ';') {      // check for non-sequence lines, ';' means a comment, TODO: change it to a robuster condition
                cur_seq += line;
            }
        }
        if (cur_seq != "") {
            // cur_seq += "\0";
            sequence_lines.push_back(cur_seq);
            seq_len.push_back(cur_seq.length());
            if (cur_seq.length() > sequence_max_length)
                sequence_max_length = cur_seq.length();
        }

        if (set_parameters) {
            seq_cnt = sequence_count; seq_max_len = sequence_max_length;
            printf("seq_cnt = %d, seq_max_len = %d\n", seq_cnt, seq_max_len);
            sequence_lengths = new int[seq_cnt];
            for (int ii = 0; ii < seq_cnt; ++ii) sequence_lengths[ii] = seq_len[ii];
        }

        target = new char*[sequence_count];
        // assign values to the target, pad gaps/spaces to the end of each sequence
        for (int i = 0; i < sequence_count; ++i) {
            target[i] = new char[sequence_max_length];
            for (int j = 0; j < sequence_max_length; ++j) {
                if (j < sequence_lines[i].length())
                    target[i][j] = sequence_lines[i][j];
                else {
                    // target[i][j] = SPACE;
                    target[i][j] = GAP;
                }
            }
        }

        for (int i = 0; i < sequence_count; ++i) {
            printf("In load sequeunces, target sequence %d: ", i);
            for (int c = 0; c < sequence_max_length; ++c)
                printf("%c", target[i][c]);
            printf("\n");
        }
        ret[0] = sequence_count; ret[1] = sequence_max_length;
    } else {
        printf("Failed to read the FASTA file!\n");
    }
    return ret;
}

// The computed MSA with the highest SP-score is written into the global variable anytime_Astar_msa_result
void compute_anytime_MSA_Astar(int bin_cnt, int beam_width, int astar_iter_cnt, ScoreTable blosum_table, double time_limit, bool enable_recursive_MSA) {
    const int seq_cnt_c = seq_cnt;

    STYPE heu_score = cost_instead_of_score ? INT_MAX : INT_MIN;

    int alignment_len;
    bool reached_time_limit = false;
    Timer anytime_Astar_timer;
    double computation_time;
    STYPE score = 0;


    if (blosum_table.get_open_penalty() == 0) {
        printf("\n\nBefore Anytime Astar algorithm (linear gap penalty)\n");
        anytime_Astar_timer.start();
        if (cost_instead_of_score) {
            SWITCH_COST_LINEAR_ASTAR(seq_cnt);
        } else {
            SWITCH_SCORE_LINEAR_ASTAR(seq_cnt);
        }
        computation_time = anytime_Astar_timer.elapsed(true);
        printf("Anytime Astar algorithm (linear gap penalty) completed, alignment length = %d, reached_time_limit = %d\n", alignment_len, (int)reached_time_limit);
        score = blosum_table.calc_score(anytime_Astar_msa_result, seq_cnt, alignment_len);
    } else {
        printf("\n\nBefore Anytime Astar algorithm (affine gap penalty)\n");
        anytime_Astar_timer.start();
        if (cost_instead_of_score) {
            SWITCH_COST_AFFINE_ASTAR(seq_cnt);
        } else {
            SWITCH_SCORE_AFFINE_ASTAR(seq_cnt);
        }
        computation_time = anytime_Astar_timer.elapsed(true);
        printf("Anytime Astar algorithm completed (affine gap penalty), alignment length = %d, reached_time_limit = %d\n", alignment_len, (int)reached_time_limit);
        score = blosum_table.calc_score(anytime_Astar_msa_result, seq_cnt, alignment_len);
    }

    show_result(anytime_Astar_msa_result, score, seq_cnt, alignment_len);
}

// argv[0] = executable, -f: input FASTA file, -t: score table, -m: memory limit (GB)
std::map<std::string, std::string> parse_params(int argc, char* argv[]) {
    std::map<std::string, std::string> params;
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-') {
            if (!strcmp(argv[i], "-f")) {
                params["file"] = std::string(argv[++i]);
            } else if (!strcmp(argv[i], "-t")) {
                params["table"] = std::string(argv[++i]);
            } else if (!strcmp(argv[i], "-m")) {
                params["memory"] = std::string(argv[++i]);
            } else if (!strcmp(argv[i], "-bc")) {
                params["bin count"] = std::string(argv[++i]);
            } else if (!strcmp(argv[i], "-bw")) {
                params["beam width"] = std::string(argv[++i]);
            } else if (!strcmp(argv[i], "-non_recur")) {
                params["non recursive"] = "1";
            }
        } else {
            printf("Invalid input parameters!\n-f: input FASTA file, -t: score table, -m: memory limit (GB)\n");
            break;
        }
    }   
    return params;
}

// argv[0] = executable, -f: input FASTA file, -t: score table, -m: memory limit (GB)
int main(int argc, char* argv[]) {
    int bin_cnt = 4, beam_width = 10, astar_iter_cnt = 100;     // RAM-MSA hyper-parameters
    bool enable_recursive_MSA = true;                           // RAM-MSA hyper-parameters

    std::map<std::string, std::string> params = parse_params(argc, argv);
    
    // input FASTA file
    if (params.find("file") != params.end()) {
        sequences_dir = params["file"];
        printf("Read input FASTA file from %s\n", sequences_dir.c_str());
    }

    // default score table = PAM250
    const char *score_table_dir = pam250_scores_dir;
    open_penalty = 0.00;
    ext_penalty = -30;
    cost_instead_of_score = true;

    if (params.find("table") != params.end()) {
        std::string input_score_table = params["table"];
        if (input_score_table == "BLOSUM62") {
            open_penalty = 1.53;
            ext_penalty = 0.00;
            score_table_dir = blosum_scores_dir;
            cost_instead_of_score = false;
        } else if (input_score_table == "PAM250") {
            open_penalty = 0.00;
            ext_penalty = -30;      // we do score -= ext_penalty, so the PENALTY of a gap is -30
            score_table_dir = pam250_scores_dir;
            cost_instead_of_score = true;
        } else {
            printf("Unknown input score table! Default is BLOSUM62!\n");
        }
    }

    // score -= open_penalty + gap_len * ext_penalty;
    ScoreTable score_table(open_penalty, ext_penalty, score_table_dir);
    printf("Loaded score table successfully, cost_instead_of_score = %d\n", (int) cost_instead_of_score);

    // default memory limit ratio is 0.8

    // disable memory-bound on mac OS, codeupdate
#if defined(__APPLE__) && defined(__MACH__)
    params["memory"] = "-1";
#endif

    if (params.find("memory") != params.end()) {
        double input_memory_limit_ratio = stod(params["memory"]);
        memory_limit_ratio = input_memory_limit_ratio;
    }

    printf("Memory limit = %.2f * Available RAM\n", memory_limit_ratio);

    if (params.find("bin count") != params.end()) {
        int input_bin_cnt = stoi(params["bin count"]);
        if (input_bin_cnt < 1) {
            printf("bin count must be positive!\n");
            return -1;
        }
        bin_cnt = input_bin_cnt;
    }

    if (params.find("beam width") != params.end()) {
        int input_beam_width = stoi(params["beam width"]);
        if (input_beam_width < 0) {
            printf("bin_cnt must be non-negative!\n");
            return -1;
        }
        beam_width = input_beam_width;
    }

    printf("Bin count = %d, beam width = %d\n", bin_cnt, beam_width);

    if (params.find("non recursive") != params.end()) {
        enable_recursive_MSA = false;
        printf("Recursive MSA disabled!\n");
    }

    load_sequences(sequences_dir, sequences, true);     // true: set_parameters
    printf("Loaded sequences successfully\n");
    init_res_arr();     // initialize the results arrays with seq_max_len

    compute_anytime_MSA_Astar(bin_cnt, beam_width, astar_iter_cnt, score_table, anytime_time_limit, enable_recursive_MSA);

    dealloc();
}