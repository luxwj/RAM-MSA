#ifndef NW_HPP
#define NW_HPP

#include "../parameters.hpp"
#include "../utils.hpp"

#include <vector>
#include <iostream>
#include <cmath>    // pow
#include <algorithm> // max

class NWSolver {
public:
    struct DPCell {
        STYPE score;
        std::vector<int> gap_lens;      // the gap length of each sequence until the current node, 0 if current is a symbol
    };

    // scores[1] = score when direction == 1, scores[3] == score when direction == 3, scores[0] is unused
    struct DPCellAffine2D {
        STYPE scores[4];
    };

private:
    bool cost_instead_of_score;
    const int nn;   // seq_cnt
    int ll;   // seq_max_len
    int *seq_lens;  // sequence lengths
    // the index offsets of each dimension, init as [1, (ll+1), (ll+1)^2, ...]
    int *offsets;

    // sequences[][i] corresponds to DP_mat[][i+1] because DP_mat has an edge length of ll+1
    char const* const*sequences;       
    ScoreTable blosum_table;
    // the above variables are initialized in instantiation
    double access_cnt;

    // These two arrays are init when calling NW_MSA()
    DPCell *DP_mat;        // size = (seq_lens[0]+1) * (seq_lens[1]+1) * ... * (seq_lens[nn-1]+1)
    int *TB_mat;        // traceback matrix, size = (seq_lens[0]+1) * (seq_lens[1]+1) * ... * (seq_lens[nn-1]+1)
    // If the score of cell x comes from a match of sequence i: (TB_mat[x] >> i) & 0x01

    // For example, consider the following alignment and matrix: 
    // S1: A _ _ B C
    // S2: _ A B _ C
    //  DP_mat                  TB_mat
    //          A   B   C               A   B   C
    //      0   -                   11  01
    //  A       |               A       10
    //  B       |   -           B       10  01
    //  C               \\      C               11

    // When computing cell (A, B), we access TB_mat(A, A) and get 10.
    // The result 10 means there is already a gap in seq1. Then, we only apply the ext_penalty to the score.


    // the results of get_DP_mat(a, a, 0) are the same in 2D and 3D tasks
    // dim is the dimension of indices
    template<typename T>
    void set_mat(T *mat, T val, int const *indices, int dim) {
        int m_idx = 0;
        for (int i = 0; i < dim; ++i)
            m_idx += indices[i] * offsets[i];
        mat[m_idx] = val;
    }

    // dim is the dimension of indices
    template<typename T>
    T get_mat(T *mat, int const *indices, int dim) {
        int m_idx = 0;
        for (int i = 0; i < dim; ++i)
            m_idx += indices[i] * offsets[i];
        return mat[m_idx];
    }

    // used in solving 2D problem
    template<typename T>
    void set_mat_offset(T *mat, T val, int idx0, int idx1, int idx0_offset, int idx1_offset) 
        {mat[idx0 * idx0_offset + idx1 * idx1_offset] = val;}
    template<typename T>
    T get_mat_offset(T *mat, int idx0, int idx1, int idx0_offset, int idx1_offset) 
        {return mat[idx0 * idx0_offset + idx1 * idx1_offset];}

    void init_mat();        // init the first rows in DP_mat and TB_mat
    void compute_NW_2D(int seq0, int seq1);
    void compute_NW();
    int backtrack(char **&msa_result);          // return the alignment length
public:
    // In addition to the passed variables, init offsets
    NWSolver(bool _cost_instead_of_score, int sequence_count, int sequence_max_length, char **original_sequences, int *sequence_lengths, ScoreTable _blosum_table):
        cost_instead_of_score{_cost_instead_of_score}, nn{sequence_count}, ll{sequence_max_length}, sequences{original_sequences}, blosum_table{_blosum_table} {
        int offset = 1;
        offsets = new int[nn]; seq_lens = new int[nn];
        for (int i = 0; i < nn; ++i) {
            seq_lens[i] = sequence_lengths[i];
            offsets[i] = offset;
            // offset *= (ll+1);
            offset *= (seq_lens[i] + 1);
        }
    }

    ~NWSolver() {
        delete[] offsets;
        delete[] seq_lens;
    }

    // set new sequences and reset variables
    // sequence count must be identical
    bool set_sequences(int sequence_count, int sequence_max_length, char **original_sequences, int *sequence_lengths) {
        if (sequence_count != nn) {
            printf("Unidentical sequence count!\n");
            return false;
        }
        ll = sequence_max_length; sequences = original_sequences;
        int offset = 1;
        for (int i = 0; i < nn; ++i) {
            seq_lens[i] = sequence_lengths[i];
            offsets[i] = offset;
            offset *= (seq_lens[i] + 1);
        }
        return true;
    }

    double NW_MSA(char **&msa_result, int &alignment_len);    // allocate (ll+1)^nn cells for the DP array. return the workload
    double NW_MSA(DPCell *& _DP_mat, int *& _TB_mat);  // output the DP matrix and traceback matrix. return the workload
};

// init the first cell in DP_mat and TB_mat
inline void NWSolver::init_mat() {
    // codeupdate
    // int zero_vec[nn] = {0};
    int zero_vec[nn];
    for (int i = 0; i < nn; ++i) zero_vec[i] = 0;

    DPCell orig_cell{0, std::vector<int>(nn, 0)};
    set_mat(DP_mat, orig_cell, zero_vec, nn);
    set_mat(TB_mat, 0, zero_vec, nn);               // set TB_mat[0] to all 0, easy to implement backtrack()
}

// seq0, seq1 are the indices of sequences
inline void NWSolver::compute_NW_2D(int seq0, int seq1) {
    STYPE open_pen = blosum_table.get_open_penalty(), ext_pen = blosum_table.get_ext_penalty();

    std::vector<DPCellAffine2D>DP_mat_affine_2D((seq_lens[seq0] + 1) * (seq_lens[seq1] + 1));
    DP_mat_affine_2D[0] = {{0, 0, 0, 0}};     // zeros in all three directions

    // first row
    for (int idx1 = 1; idx1 < seq_lens[seq1] + 1; ++idx1) {
        int idx0 = 0;
        STYPE score_dir2 = DP_mat_affine_2D[idx0 * offsets[seq0] + (idx1 - 1) * offsets[seq1]].scores[2];   // gap in seq0
        access_cnt += 1;

        bool gap_open_seq0 = (idx1 == 1) ? true : false;
        score_dir2 -= gap_open_seq0 * open_pen + ext_pen;

        DP_mat[idx0 * offsets[seq0] + idx1 * offsets[seq1]].score = score_dir2;
        DP_mat_affine_2D[idx0 * offsets[seq0] + idx1 * offsets[seq1]].scores[1] = cost_instead_of_score ? INT16_MAX : INT16_MIN;    // cannot choose this direction
        DP_mat_affine_2D[idx0 * offsets[seq0] + idx1 * offsets[seq1]].scores[2] = score_dir2;
        DP_mat_affine_2D[idx0 * offsets[seq0] + idx1 * offsets[seq1]].scores[3] = cost_instead_of_score ? INT16_MAX : INT16_MIN;    // cannot choose this direction
        // DP_mat_affine_2D[idx0 * offsets[seq0] + idx1 * offsets[seq1]].scores[3] = 0;     // score_dir3 will no be accessed

        int tb_res = (1<<seq1);     // because gap in seq0, result from seq1
        TB_mat[idx0 * offsets[seq0] + idx1 * offsets[seq1]] = tb_res;
        access_cnt += 2;
    }

    for (int idx0 = 1; idx0 < seq_lens[seq0] + 1; ++idx0) {
        // first col
        int tmp_idx1 = 0;
        STYPE tmp_score_dir1 = DP_mat[(idx0 - 1) * offsets[seq0] + tmp_idx1 * offsets[seq1]].score;         // gap in seq1
        access_cnt += 1;

        bool tmp_gap_open_seq1 = (idx0 == 1) ? true : false;
        tmp_score_dir1 -= tmp_gap_open_seq1 * open_pen + ext_pen;

        DP_mat[idx0 * offsets[seq0] + tmp_idx1 * offsets[seq1]].score = tmp_score_dir1;
        DP_mat_affine_2D[idx0 * offsets[seq0] + tmp_idx1 * offsets[seq1]].scores[1] = tmp_score_dir1;
        DP_mat_affine_2D[idx0 * offsets[seq0] + tmp_idx1 * offsets[seq1]].scores[2] = cost_instead_of_score ? INT16_MAX : INT16_MIN;   // cannot choose this direction
        DP_mat_affine_2D[idx0 * offsets[seq0] + tmp_idx1 * offsets[seq1]].scores[3] = cost_instead_of_score ? INT16_MAX : INT16_MIN;   // cannot choose this direction
        // DP_mat_affine_2D[idx0 * offsets[seq0] + tmp_idx1 * offsets[seq1]].scores[3] = 0;     // score_dir3 will no be accessed

        int tmp_tb_res = (1<<seq0);     // because gap in seq1, result from seq0
        TB_mat[idx0 * offsets[seq0] + tmp_idx1 * offsets[seq1]] = tmp_tb_res;
        access_cnt += 2;

        for (int idx1 = 1; idx1 < seq_lens[seq1] + 1; ++idx1) {
            STYPE score_dir11 = DP_mat_affine_2D[(idx0 - 1) * offsets[seq0] + idx1 * offsets[seq1]].scores[1];
            STYPE score_dir12 = DP_mat_affine_2D[(idx0 - 1) * offsets[seq0] + idx1 * offsets[seq1]].scores[2];
            STYPE score_dir13 = DP_mat_affine_2D[(idx0 - 1) * offsets[seq0] + idx1 * offsets[seq1]].scores[3];
            STYPE score_dir1;
            if (cost_instead_of_score)  score_dir1 = std::min(score_dir11 - ext_pen, std::min(score_dir12 - open_pen - ext_pen, score_dir13 - open_pen - ext_pen));
            else                        score_dir1 = std::max(score_dir11 - ext_pen, std::max(score_dir12 - open_pen - ext_pen, score_dir13 - open_pen - ext_pen));

            STYPE score_dir21 = DP_mat_affine_2D[idx0 * offsets[seq0] + (idx1 - 1) * offsets[seq1]].scores[1];
            STYPE score_dir22 = DP_mat_affine_2D[idx0 * offsets[seq0] + (idx1 - 1) * offsets[seq1]].scores[2];
            STYPE score_dir23 = DP_mat_affine_2D[idx0 * offsets[seq0] + (idx1 - 1) * offsets[seq1]].scores[3];
            STYPE score_dir2;
            if (cost_instead_of_score)  score_dir2 = std::min(score_dir22 - ext_pen, std::min(score_dir21 - open_pen - ext_pen, score_dir23 - open_pen - ext_pen));
            else                        score_dir2 = std::max(score_dir22 - ext_pen, std::max(score_dir21 - open_pen - ext_pen, score_dir23 - open_pen - ext_pen));
            
            STYPE score_dir3 = DP_mat[(idx0 - 1) * offsets[seq0] + (idx1 - 1) * offsets[seq1]].score;   // match/mismatch
            int c0 = sequences[seq0][idx0-1], c1 = sequences[seq1][idx1-1];
            score_dir3 += blosum_table.get_score_char(c0, c1);
            access_cnt += 7;

            STYPE max_score;
            if (cost_instead_of_score)  max_score = std::min(score_dir3, std::min(score_dir2, score_dir1));
            else                        max_score = std::max(score_dir3, std::max(score_dir2, score_dir1));

            DP_mat[idx0 * offsets[seq0] + idx1 * offsets[seq1]].score = max_score;
            // we don't care about gap_lens in NW_2D
            // DP_mat[idx0 * offsets[seq0] + idx1 * offsets[seq1]].gap_lens = std::vector<int>(nn, 0);
            DP_mat_affine_2D[idx0 * offsets[seq0] + idx1 * offsets[seq1]].scores[1] = score_dir1;
            DP_mat_affine_2D[idx0 * offsets[seq0] + idx1 * offsets[seq1]].scores[2] = score_dir2;
            DP_mat_affine_2D[idx0 * offsets[seq0] + idx1 * offsets[seq1]].scores[3] = score_dir3;
            access_cnt += 4;

            int tb_res = 0;
            if (max_score == score_dir1) {
                tb_res = (1<<seq0);     // gap in seq1, result from seq0
            } else if (max_score == score_dir2) {
                tb_res = (1<<seq1);     // gap in seq0, result from seq1
            } else {                    // max_score == score_dir3
                tb_res = (1<<seq0) | (1<<seq1);
            }

            TB_mat[idx0 * offsets[seq0] + idx1 * offsets[seq1]] = tb_res;
            access_cnt += 2;
        }

    }
}


// seq[i] are the indices of sequences
// idx_offsets[i] are the offsets for accessing mat
inline void NWSolver::compute_NW() {
    const int dim = nn;
    STYPE open_pen = blosum_table.get_open_penalty(), ext_pen = blosum_table.get_ext_penalty();
    int dir_cnt = pow(2, dim);       // direction counts
    int cur_cell[dim];
    for (int i = 0; i < dim; ++i) cur_cell[i] = 0;
    
    // if cur_cell[i] reaches ii+1, it will become 0 and cur_cell[i+1] += 1
    while (cur_cell[dim - 1] < seq_lens[dim - 1] + 1) {
        int symbols[dim];
        for(int i = 0; i < dim; ++i)
            symbols[i] = cur_cell[i] == 0 ? GAP : sequences[i][cur_cell[i] - 1];

        access_cnt += dim;

        DPCell DPCells[dir_cnt];
        STYPE scores[dir_cnt];      // idx 0 is not used; idx 1101 represents the score from (0, 0, ...,-1, -1, 0, -1) cell
        int tbs[dir_cnt];
        STYPE max_score;
        if (cost_instead_of_score)  max_score = INT32_MAX;
        else                        max_score = INT32_MIN;

        for (int dir = 1; dir < dir_cnt; ++dir) {
            int src_cell[dim];
            bool oob_flag = false;      // out-of-bound flag
            for (int ii = 0; ii < dim; ++ii) {
                src_cell[ii] = cur_cell[ii];
                if ((dir >> ii) & 0x01) src_cell[ii] -= 1;     // if dir == 0b011, then src_cell = [idx0-1, idx1-1, idx2]
                if (src_cell[ii] < 0) oob_flag = true;
            }
            if (oob_flag) continue;

            DPCells[dir] = get_mat(DP_mat, src_cell, dim);
            scores[dir] = DPCells[dir].score;
            if (dir != dir_cnt - 1)   // no need to check TB_mat if there is no gap
                tbs[dir] = get_mat(TB_mat, src_cell, dim);

            if (dir != dir_cnt - 1) access_cnt += 2;
            else access_cnt += 1;

            // update the score, check all pairs of sequences
            for (int i = 0; i < dim - 1; ++i) {
                for (int j = i + 1; j < dim; ++j) {
                    if (cur_cell[i] != src_cell[i] && cur_cell[j] != src_cell[j]) {     // check the score table
                        scores[dir] += blosum_table.get_score_char(symbols[i], symbols[j]);
                    } else if (cur_cell[i] == src_cell[i] && cur_cell[j] == src_cell[j]) {
                        // two gaps so do nothing
                    } else {
                        // use gap_lens in each DPCell to get the gap infos
                        bool gap_in_prev_seq_i = false, gap_in_prev_seq_j = false;
                        if (DPCells[dir].gap_lens[i] > 0 && DPCells[dir].gap_lens[i] > DPCells[dir].gap_lens[j])
                            gap_in_prev_seq_i = true;
                        if (DPCells[dir].gap_lens[j] > 0 && DPCells[dir].gap_lens[j] > DPCells[dir].gap_lens[i])
                            gap_in_prev_seq_j = true;

                        if (cur_cell[i] == src_cell[i]) {    // gap in seq i; TB(0, 0) is 0x00 that cannot be detected by bit operation, so add a condition
                            scores[dir] -= (!gap_in_prev_seq_i || src_cell[i] == 0 && src_cell[j] == 0) * open_pen + ext_pen;
                        } else {                                    // gap in seq j
                            scores[dir] -= (!gap_in_prev_seq_j || src_cell[i] == 0 && src_cell[j] == 0) * open_pen + ext_pen;
                        }
                    }
                    
                    
                }
            }

            bool better_result = cost_instead_of_score ? (scores[dir] < max_score) : (scores[dir] > max_score);
            if (better_result) {
                max_score = scores[dir];
                DPCell cur_DPCell{max_score, std::vector<int>(nn, 0)};
                for (int ii = 0; ii < nn; ++ii) {       // check gap_len of each seq
                    if (cur_cell[ii] != src_cell[ii]) {
                        cur_DPCell.gap_lens[ii] = 0;
                    } else {
                        cur_DPCell.gap_lens[ii] = DPCells[dir].gap_lens[ii] + 1;
                    }
                }
                set_mat(DP_mat, cur_DPCell, cur_cell, dim);
                set_mat(TB_mat, dir, cur_cell, dim);
                access_cnt += 2;
            }
        }

        // update cur_cell
        ++cur_cell[0];
        for (int i = 0; i < dim; ++i) {
            if (cur_cell[i] == seq_lens[i] + 1) {
                if (i == dim - 1) break;
                cur_cell[i] = 0;
                ++cur_cell[i + 1];
            }
        }
    }
}


// write the result into msa_result (begin with the ends of sequence) and return the aligment length
inline int NWSolver::backtrack(char **&msa_result) {
    int seq_indices[nn];
    for (int i = 0; i < nn; ++i)
        seq_indices[i] = seq_lens[i];            // indices of the TB_mat, use --idx{i} when accessing sequences[][]
    int res_idx = 0;
    bool all_zero_indices = false;
    while (!all_zero_indices) {    // until all indices are 0
        int tb = get_mat(TB_mat, seq_indices, nn);
        for (int i = 0; i < nn; ++i)
            msa_result[i][res_idx] = ((tb >> i) & 0x01) ? sequences[i][--seq_indices[i]] : GAP;

        ++res_idx;
        // check the loop condition
        all_zero_indices = true;
        for (int i = 0; i < nn; ++i)
            if (seq_indices[i] != 0) {
                all_zero_indices = false;
                break;
            }
    }
    // reverse the results
    for (int i = 0; i < res_idx / 2; ++i) {
        for (int seq = 0; seq < nn; ++seq) {
            char tmp = msa_result[seq][i];
            msa_result[seq][i] = msa_result[seq][res_idx-1-i];
            msa_result[seq][res_idx-1-i] = tmp;
        }
    }
    return res_idx;
}

inline double NWSolver::NW_MSA(char **&msa_result, int &alignment_len) {
    access_cnt = 0;
    long long DP_mat_size = 1;
    for (int ii = 0; ii < nn; ++ii) {
        DP_mat_size *= (seq_lens[ii] + 1);        // use individual seq_len
    }
    DP_mat = new DPCell[DP_mat_size];
    TB_mat = new int[DP_mat_size];
    init_mat();
    // codeupdate
    // int seqs[nn] = {0};
    int seqs[nn];
    for (int i = 0; i < nn; ++i) seqs[i] = 0;

    for (int i = 1; i < nn; ++i) seqs[i] = i;
    compute_NW();

    // init the result with gaps, assuming the size of the msa_result is [nn][nn*ll]
    for (int row = 0; row < nn; ++row)
        for (int col = 0; col < nn * ll; ++col) {
            msa_result[row][col] = GAP;
        }
    alignment_len = backtrack(msa_result);
    delete[] DP_mat;
    delete[] TB_mat;

    return access_cnt;
}

// output the DP matrix and traceback matrix. return the workload
inline double NWSolver::NW_MSA(DPCell *& _DP_mat, int *& _TB_mat) {
    access_cnt = 0;
    DP_mat = _DP_mat; 
    TB_mat = _TB_mat;
    init_mat();
    if (nn == 2) {
        compute_NW_2D(0, 1);
    } else {
        compute_NW();
    }
    return access_cnt;
}

#endif // NW_HPP