#ifndef PAIRWISE_HEURISTIC_HPP
#define PAIRWISE_HEURISTIC_HPP

#include "../parameters.hpp"
#include "../utils.hpp"

#include "../NW/NW.hpp"     // use NW algorithm to compute 2D scores
/*
    Usage:
    NWSolver NW_solver(seq_cnt, seq_max_len, sequences, blosum_table);
    NW_solver.NW_MSA(NW_msa_result);
    STYPE score = blosum_table.calc_score(NW_msa_result, seq_cnt, seq_cnt * seq_max_len);
*/


#include <cstring>          // memcpy

// A helper used in Astar MSA. It computes all PSAs beforehand, stores the scores in DP matrices, and outputs the score to Astar MSA.
class PSA_heuristic {
private:
    // compute the path cost instead of the score
    bool cost_instead_of_score;

    const int nn;
    const int ll;
    int *seq_lens;  // sequence lengths
    char **sequences;       // size = [nn][ll]
    ScoreTable blosum_table;
    NWSolver *NW_solver;

    // nn * nn matrices in total, index is computed by i * nn + j, 0 <= i < j < nn. Space at some indices are not allocated.
    // Each matrix has a size = (seq_lens[i]+1) * (seq_lens[j]+1)
    NWSolver::DPCell **DP_mats;
    int **TB_mats;              // traceback matrix

    NWSolver::DPCell **reversed_DP_mats;
    int **reversed_TB_mats;     // traceback matrix

    inline int mat_idx(int i, int j) {return (i * nn + j);};

    // because we want the DP scores begin from the bottom-right corner
    void reverse_seqs() {
        for (int idx = 0; idx < nn; ++idx) {
            std::reverse(sequences[idx], sequences[idx] + seq_lens[idx]);
        }
    }

public:
    PSA_heuristic(bool _cost_instead_of_score, int sequence_count, int sequence_max_length, char **original_sequences, int *sequence_lengths, ScoreTable _blosum_table):
        cost_instead_of_score{ _cost_instead_of_score}, nn{sequence_count}, ll{sequence_max_length}, blosum_table{_blosum_table} {
        seq_lens = new int[nn];
        sequences = new char*[nn];
        for (int i = 0; i < nn; ++i) {
            seq_lens[i] = sequence_lengths[i];
            sequences[i] = new char[seq_lens[i]];
            memcpy(sequences[i], original_sequences[i], seq_lens[i] * sizeof(char));
        }
        
        NW_solver = new NWSolver(cost_instead_of_score, 2, ll, original_sequences, seq_lens, blosum_table);  // pass original_sequences and do nothing
        // init DP_mats and TB_mats
        DP_mats = new NWSolver::DPCell*[nn * nn]; TB_mats = new int*[nn * nn];
        for (int i = 0; i < nn - 1; ++i)
            for (int j = i + 1; j < nn; ++j) {
                long long DP_mat_size = (seq_lens[i] + 1) * (seq_lens[j] + 1);
                DP_mats[mat_idx(i, j)] = new NWSolver::DPCell[DP_mat_size];
                TB_mats[mat_idx(i, j)] = new int[DP_mat_size];
            }

        reversed_DP_mats = new NWSolver::DPCell*[nn * nn]; reversed_TB_mats = new int*[nn * nn];
        for (int i = 0; i < nn - 1; ++i)
            for (int j = i + 1; j < nn; ++j) {
                long long DP_mat_size = (seq_lens[i] + 1) * (seq_lens[j] + 1);
                reversed_DP_mats[mat_idx(i, j)] = new NWSolver::DPCell[DP_mat_size];
                reversed_TB_mats[mat_idx(i, j)] = new int[DP_mat_size];
            }
    }

    ~PSA_heuristic() {
        delete [] seq_lens;
        delete NW_solver;

        for (int i = 0; i < nn; ++i)
            delete [] sequences[i];
        delete [] sequences;

        for (int i = 0; i < nn - 1; ++i)
            for (int j = i + 1; j < nn; ++j) {
                delete [] DP_mats[mat_idx(i, j)];
                delete [] TB_mats[mat_idx(i, j)];
                delete [] reversed_DP_mats[mat_idx(i, j)];
                delete [] reversed_TB_mats[mat_idx(i, j)];
            }
        delete [] DP_mats;
        delete [] TB_mats;
        delete [] reversed_DP_mats;
        delete [] reversed_TB_mats;
    }

    // return the workload
    double compute_all_PSA() {
        double access_cnt = 0;
        char* tmp_seqs[2];
        int cur_seq_lens[2];
        // first compute the normal DP matrices, then compute the reversed DP matrices
        for (int i = 0; i < nn - 1; ++i) 
            for (int j = i + 1; j < nn; ++j) {
                tmp_seqs[0] = sequences[i]; tmp_seqs[1] = sequences[j];
                int seq_max_len = seq_lens[i] < seq_lens[j] ? seq_lens[j] : seq_lens[i];
                cur_seq_lens[0] = seq_lens[i]; cur_seq_lens[1] = seq_lens[j];
                NW_solver->set_sequences(2, seq_max_len, tmp_seqs, cur_seq_lens);
                access_cnt += NW_solver->NW_MSA(DP_mats[mat_idx(i, j)], TB_mats[mat_idx(i, j)]);
            }

        // flip the sequences and compute the reversed DP matrices
        reverse_seqs();
        for (int i = 0; i < nn - 1; ++i) 
            for (int j = i + 1; j < nn; ++j) {
                tmp_seqs[0] = sequences[i]; tmp_seqs[1] = sequences[j];
                int seq_max_len = seq_lens[i] < seq_lens[j] ? seq_lens[j] : seq_lens[i];
                cur_seq_lens[0] = seq_lens[i]; cur_seq_lens[1] = seq_lens[j];
                NW_solver->set_sequences(2, seq_max_len, tmp_seqs, cur_seq_lens);
                access_cnt += NW_solver->NW_MSA(reversed_DP_mats[mat_idx(i, j)], reversed_TB_mats[mat_idx(i, j)]);
            }
        return access_cnt;
    }

    // output the corresponding values in the reversed DP and TB matrices
    void get_score(bool access_reversed_mats, int seq_i, int seq_j, int symbol_i, int symbol_j, STYPE &DP_val, int &TB_val) {
        // For a sequence ACGT, the reversed DP matrix is like: 
        //          T   G   C   A
        //      0   1   2   3   4
        //  T   1
        //  G   2
        //  C   3
        //  A   4

        // normal DP matrix is like: 
        //          A   C   G   T
        //      0   1   2   3   4
        //  A   1
        //  C   2
        //  G   3
        //  T   4

        int reversed_idx_i = seq_lens[seq_i] - symbol_i, reversed_idx_j = seq_lens[seq_j] - symbol_j;
        if (access_reversed_mats) {
            // For instance, node (A, G) is symbol (0, 2). It reversed index is (4, 2).
            // int reversed_idx_i = seq_lens[seq_i] - symbol_i, reversed_idx_j = seq_lens[seq_j] - symbol_j;
            // Warning: The offset here is hard-coded (copied from NW.hpp). This should be changed accordingly when NW.hpp changes.
            DP_val = reversed_DP_mats[mat_idx(seq_i, seq_j)][reversed_idx_i + reversed_idx_j * (seq_lens[seq_i] + 1)].score;
            TB_val = reversed_TB_mats[mat_idx(seq_i, seq_j)][reversed_idx_i + reversed_idx_j * (seq_lens[seq_i] + 1)];
        } else {
            // The Astar node is reversed: normal = {ACGT, ACGT}, reversed = {TGCA, TGCA}
            // cur_crd = {TGC, T}, symbol_i = 3, symbol_j = 1. So we return normal DP[1][3] representing {A, ACG}
            // Warning: The offset here is hard-coded (copied from NW.hpp). This should be changed accordingly when NW.hpp changes.
            DP_val = DP_mats[mat_idx(seq_i, seq_j)][reversed_idx_i + reversed_idx_j * (seq_lens[seq_i] + 1)].score;
            TB_val = TB_mats[mat_idx(seq_i, seq_j)][reversed_idx_i + reversed_idx_j * (seq_lens[seq_i] + 1)];
        }
    }

};

#endif // PAIRWISE_HEURISTIC_HPP