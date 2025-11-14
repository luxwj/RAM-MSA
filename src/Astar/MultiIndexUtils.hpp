#ifndef MULTIINDEXUTILS_HPP
#define MULTIINDEXUTILS_HPP

#include "../parameters.hpp"
#include "../utils.hpp"
#include <array>    // codeupdate
#include <vector>
#include <boost/multi_index_container.hpp>
#include <boost/multi_index/member.hpp>
#include <boost/multi_index/hashed_index.hpp>
#include <boost/multi_index/ordered_index.hpp>
#include <boost/multi_index/composite_key.hpp>  // tie breaker
#include <functional>   // for std::hash

// An instance of AstarMultiIndex holds an OpenList
template <int NN>
class AstarMultiIndex{
public:
    template <typename TB_TYPE>
    struct CoordTBInfo {         // hash map key
        std::array<int, NN> crd;

        // // When each node includes the gap length, tb[0] = gap_len[0], tb[1] = gap_len[1], ...
        TB_TYPE tb_info;        // trace back info

        bool operator==(const CoordTBInfo& other) const {
            return crd == other.crd && tb_info == other.tb_info;
        }
    };

    struct CoordTBInfoHash {
        std::size_t operator()(const CoordTBInfo<AFFINE_TB_TYPE>& k) const {
            // VectorIntHash in utils.hpp
            std::size_t h1 = VectorIntHash<NN>{}(k.crd);
            std::size_t h2 = VectorIntHash<NN>{}(k.tb_info);
            return h1 ^ (h2 << 1);
        }

        // // only consider the coordinates
        // std::size_t operator()(const CoordTBInfo<LINEAR_TB_TYPE>& k) const {
        //     // VectorIntHash in utils.hpp
        //     std::size_t h1 = VectorIntHash<NN>{}(k.crd);
        //     return h1;
        // }

        // consider both the coordinates and the tb_info
        std::size_t operator()(const CoordTBInfo<LINEAR_TB_TYPE>& k) const {
            // VectorIntHash in utils.hpp
            std::size_t h1 = VectorIntHash<NN>{}(k.crd);
            std::size_t h2 = VectorIntHash<NN>{}(k.tb_info);
            return h1 ^ (h2 << 1);
        }
    };

    template <typename TB_TYPE>
    struct OpenListNode {
        CoordTBInfo<TB_TYPE> crd_tb_info;
        STYPE fscore;
        STYPE gscore;
    };

    /**
     * This one not only changes the gscore but also modifies fscore accordingly.
     * fscore += new_gscore - node.gscore;
     * Usage: 
     * iter = indexByCoord.find(some_idx);
     * indexByCoord.modify(iter, ChangeGscore(new_gscore));
     */
    template <typename TB_TYPE>
    struct ChangeGscore {
        ChangeGscore(const STYPE& new_gscore):new_gscore(new_gscore){}
        
        void operator()(OpenListNode<TB_TYPE>& node) {
            node.fscore += new_gscore - node.gscore;
            node.gscore = new_gscore;
        }

        private:
            STYPE new_gscore;
    };

    /**
     * @brief Update the fscore of a reexpanded parent node when enabling memory bound
     * Usage: 
     * iter = indexByCoord.find(some_idx);
     * indexByCoord.modify(iter, ChangeFscore(new_fscore));
     */
    template <typename TB_TYPE>
    struct ChangeFscore {
        ChangeFscore(const STYPE& new_fscore):new_fscore(new_fscore){}
        
        void operator()(OpenListNode<TB_TYPE>& node) {
            node.fscore = new_fscore;
        }

        private:
            STYPE new_fscore;
    };

    struct IndexByFscore {};
    struct IndexByCoord {}; // coordinates and came from direction

    using OpenListMultiIdx_Score_Affine = boost::multi_index_container<
        OpenListNode<AFFINE_TB_TYPE>,     // the data type stored, AFFINE_TB_TYPE for affine gap penalty
        boost::multi_index::indexed_by<     // list of indexes
            boost::multi_index::ordered_non_unique<
                boost::multi_index::tag<IndexByFscore>,
                boost::multi_index::composite_key<
                    OpenListNode<AFFINE_TB_TYPE>,
                    boost::multi_index::member<OpenListNode<AFFINE_TB_TYPE>, STYPE, &OpenListNode<AFFINE_TB_TYPE>::fscore>,     // primary
                    boost::multi_index::member<OpenListNode<AFFINE_TB_TYPE>, STYPE, &OpenListNode<AFFINE_TB_TYPE>::gscore>      // secondary as tie breaker
                >,
                boost::multi_index::composite_key_compare<
                    std::greater<STYPE>,        // fscore in descending order
                    std::greater<STYPE>         // gscore in descending order
                >
            >,
            boost::multi_index::hashed_non_unique<
                boost::multi_index::tag<IndexByCoord>,
                boost::multi_index::member<OpenListNode<AFFINE_TB_TYPE>, CoordTBInfo<AFFINE_TB_TYPE>, &OpenListNode<AFFINE_TB_TYPE>::crd_tb_info>, // what will be the index's key
                CoordTBInfoHash
            >
        >
    >;

    using OpenListMultiIdx_Score_Linear = boost::multi_index_container<
        OpenListNode<LINEAR_TB_TYPE>,     // the data type stored, LINEAR_TB_TYPE for linear gap penalty
        boost::multi_index::indexed_by<     // list of indexes
            boost::multi_index::ordered_non_unique<
                boost::multi_index::tag<IndexByFscore>,
                boost::multi_index::composite_key<
                    OpenListNode<LINEAR_TB_TYPE>,
                    boost::multi_index::member<OpenListNode<LINEAR_TB_TYPE>, STYPE, &OpenListNode<LINEAR_TB_TYPE>::fscore>,     // primary
                    boost::multi_index::member<OpenListNode<LINEAR_TB_TYPE>, STYPE, &OpenListNode<LINEAR_TB_TYPE>::gscore>      // secondary as tie breaker
                >,
                boost::multi_index::composite_key_compare<
                    std::greater<STYPE>,        // fscore in descending order
                    std::greater<STYPE>         // gscore in descending order
                >
            >,
            boost::multi_index::hashed_non_unique<
                boost::multi_index::tag<IndexByCoord>,
                boost::multi_index::member<OpenListNode<LINEAR_TB_TYPE>, CoordTBInfo<LINEAR_TB_TYPE>, &OpenListNode<LINEAR_TB_TYPE>::crd_tb_info>, // what will be the index's key
                CoordTBInfoHash
            >
        >
    >;

    using OpenListByFscore_Score_Affine = typename boost::multi_index::index<OpenListMultiIdx_Score_Affine, IndexByFscore>::type;
    using OpenListByCoord_Score_Affine = typename boost::multi_index::index<OpenListMultiIdx_Score_Affine, IndexByCoord>::type;

    using OpenListByFscore_Score_Linear = typename boost::multi_index::index<OpenListMultiIdx_Score_Linear, IndexByFscore>::type;
    using OpenListByCoord_Score_Linear = typename boost::multi_index::index<OpenListMultiIdx_Score_Linear, IndexByCoord>::type;
  

    /* --- cost instead of score, fscore and gscore in ascending order --- */
    using OpenListMultiIdx_Cost_Affine = boost::multi_index_container<
        OpenListNode<AFFINE_TB_TYPE>,         // the data type stored, AFFINE_TB_TYPE for affine gap penalty
        boost::multi_index::indexed_by<             // list of indexes
            boost::multi_index::ordered_non_unique<
                boost::multi_index::tag<IndexByFscore>,
                boost::multi_index::composite_key<
                    OpenListNode<AFFINE_TB_TYPE>,
                    boost::multi_index::member<OpenListNode<AFFINE_TB_TYPE>, STYPE, &OpenListNode<AFFINE_TB_TYPE>::fscore>,     // primary
                    boost::multi_index::member<OpenListNode<AFFINE_TB_TYPE>, STYPE, &OpenListNode<AFFINE_TB_TYPE>::gscore>      // secondary as tie breaker
                >,
                boost::multi_index::composite_key_compare<
                    std::less<STYPE>,        // fscore in ascending order
                    std::greater<STYPE>         // updated: gscore in descending order
                >
            >,
            boost::multi_index::hashed_non_unique<
                boost::multi_index::tag<IndexByCoord>,
                boost::multi_index::member<OpenListNode<AFFINE_TB_TYPE>, CoordTBInfo<AFFINE_TB_TYPE>, &OpenListNode<AFFINE_TB_TYPE>::crd_tb_info>, // what will be the index's key
                CoordTBInfoHash
            >
        >
    >;

    using OpenListMultiIdx_Cost_Linear = boost::multi_index_container<
        OpenListNode<LINEAR_TB_TYPE>,         // the data type stored, LINEAR_TB_TYPE for linear gap penalty
        boost::multi_index::indexed_by<             // list of indexes
            boost::multi_index::ordered_non_unique<
                boost::multi_index::tag<IndexByFscore>,
                boost::multi_index::composite_key<
                    OpenListNode<LINEAR_TB_TYPE>,
                    boost::multi_index::member<OpenListNode<LINEAR_TB_TYPE>, STYPE, &OpenListNode<LINEAR_TB_TYPE>::fscore>,     // primary
                    boost::multi_index::member<OpenListNode<LINEAR_TB_TYPE>, STYPE, &OpenListNode<LINEAR_TB_TYPE>::gscore>      // secondary as tie breaker
                >,
                boost::multi_index::composite_key_compare<
                    std::less<STYPE>,        // fscore in ascending order
                    std::greater<STYPE>         // updated: gscore in descending order
                >
            >,
            boost::multi_index::hashed_non_unique<
                boost::multi_index::tag<IndexByCoord>,
                boost::multi_index::member<OpenListNode<LINEAR_TB_TYPE>, CoordTBInfo<LINEAR_TB_TYPE>, &OpenListNode<LINEAR_TB_TYPE>::crd_tb_info>, // what will be the index's key
                CoordTBInfoHash
            >
        >
    >;

    using OpenListByFscore_Cost_Affine = typename boost::multi_index::index<OpenListMultiIdx_Cost_Affine, IndexByFscore>::type;
    using OpenListByCoord_Cost_Affine = typename boost::multi_index::index<OpenListMultiIdx_Cost_Affine, IndexByCoord>::type;

    using OpenListByFscore_Cost_Linear = typename boost::multi_index::index<OpenListMultiIdx_Cost_Linear, IndexByFscore>::type;
    using OpenListByCoord_Cost_Linear = typename boost::multi_index::index<OpenListMultiIdx_Cost_Linear, IndexByCoord>::type;


    /* --- cost instead of score --- */

    /*  usage
        AstarMultiIndex::OpenListMultiIdx_{Cost, Score}_{Linear, Affine} open_list_MI;
        AstarMultiIndex::OpenListByFscore_{Cost, Score}_{Linear, Affine}& indexByFscore = open_list_MI.get<AstarMultiIndex<NN>::IndexByFscore>();
        AstarMultiIndex::OpenListByCoord_{Cost, Score}_{Linear, Affine}& indexByCoord = open_list_MI.get<AstarMultiIndex<NN>::IndexByCoord>();

        insert:
            indexByCoord.insert(new_node);
        erase:
            if (!indexByFscore.empty()) {
                auto cur_node = indexByFscore.begin();
                // compute with cur_node
                indexByFscore.erase(indexByFscore.begin());
            }
        find:
            auto iterFound = indexByCoord.find(query_idx);
            if (iterFound != indexByCoord.end()) {
                found_count += 1;
            }
        modify:
            auto iter = indexByCoord.find(some_idx);
            indexByCoord.modify(iter, ChangeGscore(new_gscore));
    */
};

#endif // MULTIINDEXUTILS