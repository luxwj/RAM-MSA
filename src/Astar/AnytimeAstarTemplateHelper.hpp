#ifndef ANYTIME_ASTAR_TEMPLATE_HELPER_H
#define ANYTIME_ASTAR_TEMPLATE_HELPER_H

#define COST_LINEAR_ASTAR(X) \
    AnytimeAstarLinearGapSolver<AstarMultiIndex<X>::OpenListMultiIdx_Cost_Linear, AstarMultiIndex<X>::OpenListByFscore_Cost_Linear, AstarMultiIndex<X>::OpenListByCoord_Cost_Linear, X>  anytime_Astar_linear_gap_solver \
        (cost_instead_of_score, seq_cnt, seq_max_len, sequences, sequence_lengths, blosum_table, heu_score, \
        bin_cnt, beam_width, astar_iter_cnt, time_limit, memory_limit_ratio, anytime_astar_result_dir, enable_recursive_MSA, sequences_dir); \
    anytime_Astar_linear_gap_solver.AnytimeAstar_MSA_linear_gap(anytime_Astar_msa_result, alignment_len, reached_time_limit);

#define SCORE_LINEAR_ASTAR(X) \
    AnytimeAstarLinearGapSolver<AstarMultiIndex<X>::OpenListMultiIdx_Score_Linear, AstarMultiIndex<X>::OpenListByFscore_Score_Linear, AstarMultiIndex<X>::OpenListByCoord_Score_Linear, X>  anytime_Astar_linear_gap_solver \
        (cost_instead_of_score, seq_cnt, seq_max_len, sequences, sequence_lengths, blosum_table, heu_score, \
        bin_cnt, beam_width, astar_iter_cnt, time_limit, memory_limit_ratio, anytime_astar_result_dir, enable_recursive_MSA, sequences_dir); \
    anytime_Astar_linear_gap_solver.AnytimeAstar_MSA_linear_gap(anytime_Astar_msa_result, alignment_len, reached_time_limit);
        
#define COST_AFFINE_ASTAR(X) \
    AnytimeAstarSolver<AstarMultiIndex<X>::OpenListMultiIdx_Cost_Affine, AstarMultiIndex<X>::OpenListByFscore_Cost_Affine, AstarMultiIndex<X>::OpenListByCoord_Cost_Affine, X> anytime_Astar_solver \
        (cost_instead_of_score, seq_cnt, seq_max_len, sequences, sequence_lengths, blosum_table, heu_score, \
        bin_cnt, beam_width, astar_iter_cnt, time_limit, memory_limit_ratio, anytime_astar_result_dir, enable_recursive_MSA, sequences_dir); \
    anytime_Astar_solver.AnytimeAstar_MSA(anytime_Astar_msa_result, alignment_len, reached_time_limit);

#define SCORE_AFFINE_ASTAR(X) \
    AnytimeAstarSolver<AstarMultiIndex<X>::OpenListMultiIdx_Score_Affine, AstarMultiIndex<X>::OpenListByFscore_Score_Affine, AstarMultiIndex<X>::OpenListByCoord_Score_Affine, X> anytime_Astar_solver \
        (cost_instead_of_score, seq_cnt, seq_max_len, sequences, sequence_lengths, blosum_table, heu_score, \
        bin_cnt, beam_width, astar_iter_cnt, time_limit, memory_limit_ratio, anytime_astar_result_dir, enable_recursive_MSA, sequences_dir); \
    anytime_Astar_solver.AnytimeAstar_MSA(anytime_Astar_msa_result, alignment_len, reached_time_limit);

#define SWITCH_COST_LINEAR_ASTAR(X) \
    switch (X) { \
        case 2: {COST_LINEAR_ASTAR(2) break;} \
        case 3: {COST_LINEAR_ASTAR(3) break;} \
        case 4: {COST_LINEAR_ASTAR(4) break;} \
        case 5: {COST_LINEAR_ASTAR(5) break;} \
        case 6: {COST_LINEAR_ASTAR(6) break;} \
        case 7: {COST_LINEAR_ASTAR(7) break;} \
        case 8: {COST_LINEAR_ASTAR(8) break;} \
        case 9: {COST_LINEAR_ASTAR(9) break;} \
        case 10: {COST_LINEAR_ASTAR(10) break;} \
        case 11: {COST_LINEAR_ASTAR(11) break;} \
        case 12: {COST_LINEAR_ASTAR(12) break;} \
        case 13: {COST_LINEAR_ASTAR(13) break;} \
        case 14: {COST_LINEAR_ASTAR(14) break;} \
        case 15: {COST_LINEAR_ASTAR(15) break;} \
        case 16: {COST_LINEAR_ASTAR(16) break;} \
        case 17: {COST_LINEAR_ASTAR(17) break;} \
        case 18: {COST_LINEAR_ASTAR(18) break;} \
        case 19: {COST_LINEAR_ASTAR(19) break;} \
        case 20: {COST_LINEAR_ASTAR(20) break;} \
        default: printf("Invalid seq_cnt for Anytime Astar MSA!\n"); \
    }
    
#define SWITCH_SCORE_LINEAR_ASTAR(X) \
    switch (X) { \
        case 2: {SCORE_LINEAR_ASTAR(2) break;} \
        case 3: {SCORE_LINEAR_ASTAR(3) break;} \
        case 4: {SCORE_LINEAR_ASTAR(4) break;} \
        case 5: {SCORE_LINEAR_ASTAR(5) break;} \
        case 6: {SCORE_LINEAR_ASTAR(6) break;} \
        case 7: {SCORE_LINEAR_ASTAR(7) break;} \
        case 8: {SCORE_LINEAR_ASTAR(8) break;} \
        case 9: {SCORE_LINEAR_ASTAR(9) break;} \
        case 10: {SCORE_LINEAR_ASTAR(10) break;} \
        case 11: {SCORE_LINEAR_ASTAR(11) break;} \
        case 12: {SCORE_LINEAR_ASTAR(12) break;} \
        case 13: {SCORE_LINEAR_ASTAR(13) break;} \
        case 14: {SCORE_LINEAR_ASTAR(14) break;} \
        case 15: {SCORE_LINEAR_ASTAR(15) break;} \
        case 16: {SCORE_LINEAR_ASTAR(16) break;} \
        case 17: {SCORE_LINEAR_ASTAR(17) break;} \
        case 18: {SCORE_LINEAR_ASTAR(18) break;} \
        case 19: {SCORE_LINEAR_ASTAR(19) break;} \
        case 20: {SCORE_LINEAR_ASTAR(20) break;} \
        default: printf("Invalid seq_cnt for Anytime Astar MSA!\n"); \
    }

#define SWITCH_COST_AFFINE_ASTAR(X) \
    switch (X) { \
        case 2: {COST_AFFINE_ASTAR(2) break;} \
        case 3: {COST_AFFINE_ASTAR(3) break;} \
        case 4: {COST_AFFINE_ASTAR(4) break;} \
        case 5: {COST_AFFINE_ASTAR(5) break;} \
        case 6: {COST_AFFINE_ASTAR(6) break;} \
        case 7: {COST_AFFINE_ASTAR(7) break;} \
        case 8: {COST_AFFINE_ASTAR(8) break;} \
        case 9: {COST_AFFINE_ASTAR(9) break;} \
        case 10: {COST_AFFINE_ASTAR(10) break;} \
        case 11: {COST_AFFINE_ASTAR(11) break;} \
        case 12: {COST_AFFINE_ASTAR(12) break;} \
        case 13: {COST_AFFINE_ASTAR(13) break;} \
        case 14: {COST_AFFINE_ASTAR(14) break;} \
        case 15: {COST_AFFINE_ASTAR(15) break;} \
        case 16: {COST_AFFINE_ASTAR(16) break;} \
        case 17: {COST_AFFINE_ASTAR(17) break;} \
        case 18: {COST_AFFINE_ASTAR(18) break;} \
        case 19: {COST_AFFINE_ASTAR(19) break;} \
        case 20: {COST_AFFINE_ASTAR(20) break;} \
        default: printf("Invalid seq_cnt for Anytime Astar MSA!\n"); \
    }

#define SWITCH_SCORE_AFFINE_ASTAR(X) \
    switch (X) { \
        case 2: {SCORE_AFFINE_ASTAR(2) break;} \
        case 3: {SCORE_AFFINE_ASTAR(3) break;} \
        case 4: {SCORE_AFFINE_ASTAR(4) break;} \
        case 5: {SCORE_AFFINE_ASTAR(5) break;} \
        case 6: {SCORE_AFFINE_ASTAR(6) break;} \
        case 7: {SCORE_AFFINE_ASTAR(7) break;} \
        case 8: {SCORE_AFFINE_ASTAR(8) break;} \
        case 9: {SCORE_AFFINE_ASTAR(9) break;} \
        case 10: {SCORE_AFFINE_ASTAR(10) break;} \
        case 11: {SCORE_AFFINE_ASTAR(11) break;} \
        case 12: {SCORE_AFFINE_ASTAR(12) break;} \
        case 13: {SCORE_AFFINE_ASTAR(13) break;} \
        case 14: {SCORE_AFFINE_ASTAR(14) break;} \
        case 15: {SCORE_AFFINE_ASTAR(15) break;} \
        case 16: {SCORE_AFFINE_ASTAR(16) break;} \
        case 17: {SCORE_AFFINE_ASTAR(17) break;} \
        case 18: {SCORE_AFFINE_ASTAR(18) break;} \
        case 19: {SCORE_AFFINE_ASTAR(19) break;} \
        case 20: {SCORE_AFFINE_ASTAR(20) break;} \
        default: printf("Invalid seq_cnt for Anytime Astar MSA!\n"); \
    }

#endif