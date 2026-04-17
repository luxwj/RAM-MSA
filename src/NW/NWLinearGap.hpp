#ifndef NW_LINEAR_GAP_HPP
#define NW_LINEAR_GAP_HPP

#include "../parameters.hpp"
#include "../utils.hpp"

#include <vector>
#include <iostream>
#include <cmath>    // pow
#include <algorithm> // max

// Only computes the DP matrix and perform no traceback
class NWLinearGapSolver {
public:
    typedef STYPE DPCellLinear2D;

private:
    bool cost_instead_of_score;
    int *seq_lens;  // sequence lengths

    // sequences[][i] corresponds to DP_mat[][i+1]
    char const* const*sequences;    
    const ScoreTable &blosum_table;
    // the above variables are initialized in instantiation

    DPCellLinear2D *DP_mat;        // size = (seq_lens[0]+1) * (seq_lens[1]+1)
    
    void init_mat();        // init the first rows in DP_mat
    void compute_NW_2D();
public:
    NWLinearGapSolver(bool _cost_instead_of_score, char **original_sequences, int *sequence_lengths, const ScoreTable &_blosum_table):
        cost_instead_of_score{_cost_instead_of_score}, sequences{original_sequences}, blosum_table{_blosum_table} {
        seq_lens = new int[2];
        for (int i = 0; i < 2; ++i)
            seq_lens[i] = sequence_lengths[i];
    }

    ~NWLinearGapSolver() {
        free(seq_lens);
    }

    // set new sequences and reset variables
    // sequence count must be identical
    void set_sequences(char **original_sequences, int *sequence_lengths) {
        sequences = original_sequences;
        for (int i = 0; i < 2; ++i)
            seq_lens[i] = sequence_lengths[i];
    }

    void NW_MSA(DPCellLinear2D *& _DP_mat);  // output the DP matrix
};

inline void NWLinearGapSolver::compute_NW_2D() {
    STYPE ext_pen = blosum_table.get_ext_penalty();
    DP_mat[0] = 0;

    // first row
    for (int idx1 = 1; idx1 < seq_lens[1] + 1; ++idx1) {
        int idx0 = 0;
        STYPE pen = idx1 * ext_pen;
        DP_mat[idx0 + idx1 * (seq_lens[0] + 1)] = -pen;
    }

    for (int idx0 = 1; idx0 < seq_lens[0] + 1; ++idx0) {
        // first col
        int tmp_idx1 = 0;
        STYPE tmp_pen = idx0 * ext_pen;
        DP_mat[idx0 + tmp_idx1 * (seq_lens[0] + 1)] = -tmp_pen;

        for (int idx1 = 1; idx1 < seq_lens[1] + 1; ++idx1) {
            STYPE score_dir1 = DP_mat[(idx0 - 1) + idx1 * (seq_lens[0] + 1)] - ext_pen;
            STYPE score_dir2 = DP_mat[idx0 + (idx1 - 1) * (seq_lens[0] + 1)] - ext_pen;
            STYPE score_dir3 = DP_mat[(idx0 - 1) + (idx1 - 1) * (seq_lens[0] + 1)];   // match/mismatch
            int c0 = sequences[0][idx0 - 1], c1 = sequences[1][idx1 - 1];
            score_dir3 += blosum_table.get_score_char(c0, c1);

            STYPE max_score;
            if (cost_instead_of_score)  max_score = std::min(score_dir3, std::min(score_dir2, score_dir1));
            else                        max_score = std::max(score_dir3, std::max(score_dir2, score_dir1));

            DP_mat[idx0 + idx1 * (seq_lens[0] + 1)] = max_score;
        }

    }
}

// output the DP matrix
inline void NWLinearGapSolver::NW_MSA(DPCellLinear2D *& _DP_mat) {
    DP_mat = _DP_mat;
    compute_NW_2D();
}

#endif // NW_LINEAR_GAP_HPP