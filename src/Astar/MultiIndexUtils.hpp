#ifndef MULTIINDEXUTILS_HPP
#define MULTIINDEXUTILS_HPP

#include "../parameters.hpp"
#include "../utils.hpp"
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
    // Two types of data structure for OpenList nodes:
    // 1. Each node includes the direction to its parent and grandparent
    // 2. Each node includes the gap length of each sequence
public:

    /* --- Affine --- */
    struct CoordTBInfo {         // hash map key
        std::array<CRDTYPE, NN> crd;

        // Each node includes the gap length, tb[0] = gap_len[0], tb[1] = gap_len[1], ...
        std::array<CRDTYPE, NN> tb_info;        // trace back info

        bool operator==(const CoordTBInfo& other) const {
            return crd == other.crd && tb_info == other.tb_info;
        }
    };

    struct CoordTBInfoHash {
        size_t operator()(const CoordTBInfo& k) const {
            // VectorIntHash in utils.hpp
            size_t h1 = VectorIntHash<NN>{}(k.crd);
            size_t h2 = VectorIntHash<NN>{}(k.tb_info);
            size_t value = h1 ^ (h2 << 1);
            return value;
        }
    };

    struct OpenListNodeAffineGap {
        CoordTBInfo crd_tb_info;
        STYPE fscore;
        STYPE gscore;

        std::array<CRDTYPE, NN> get_crd() const {return crd_tb_info.crd;}
        std::array<CRDTYPE, NN> get_tb_info() const {return crd_tb_info.tb_info;}
    };

    /**
     * This one not only changes the gscore but also modifies fscore accordingly.
     * fscore += new_gscore - node.gscore;
     * Usage: 
     * iter = indexByCoord.find(some_idx);
     * indexByCoord.modify(iter, ChangeGscore(new_gscore));
     */
    struct ChangeGscoreAffineGap {
        ChangeGscoreAffineGap(const STYPE& new_gscore):new_gscore(new_gscore){}
        
        void operator()(OpenListNodeAffineGap& node) {
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
    struct ChangeFscoreAffineGap {
        ChangeFscoreAffineGap(const STYPE& new_fscore):new_fscore(new_fscore){}
        
        void operator()(OpenListNodeAffineGap& node) {
            node.fscore = new_fscore;
        }

        private:
            STYPE new_fscore;
    };

    /**
     * Usage: 
     * iter = indexByCoord.find(some_idx);
     * indexByCoord.modify(iter, ChangeGscoreFscoreTBinfo(new_gscore, new_fscore, new_tbinfo));
     */
    struct ChangeGscoreFscoreTBinfoAffineGap {
        ChangeGscoreFscoreTBinfoAffineGap(const STYPE& new_gscore, const STYPE& new_fscore, const std::array<CRDTYPE, NN>& new_tbinfo):new_gscore(new_gscore), new_fscore(new_fscore), new_tbinfo(new_tbinfo){}
        
        void operator()(OpenListNodeAffineGap& node) {
            node.fscore = new_fscore;
            node.gscore = new_gscore;
            node.crd_tb_info.tb_info = new_tbinfo;
        }

        private:
            STYPE new_gscore;
            STYPE new_fscore;
            std::array<CRDTYPE, NN> new_tbinfo;
    };

    struct IndexByFscore {};
    struct IndexByCoord {}; // coordinates and came from direction

    using OpenListMultiIdx_Score_Affine = boost::multi_index_container<
        OpenListNodeAffineGap,              // the data type stored
        boost::multi_index::indexed_by<     // list of indexes
            boost::multi_index::ordered_non_unique<
                boost::multi_index::tag<IndexByFscore>,
                boost::multi_index::member<
                    OpenListNodeAffineGap, 
                    STYPE, 
                    &OpenListNodeAffineGap::fscore
                >,
                std::greater<STYPE>        // fscore in descending order
            >,
            boost::multi_index::hashed_non_unique<
                boost::multi_index::tag<IndexByCoord>,
                boost::multi_index::member<OpenListNodeAffineGap, CoordTBInfo, &OpenListNodeAffineGap::crd_tb_info>, // what will be the index's key
                CoordTBInfoHash
            >
        >
    >;


    using OpenListByFscore_Score_Affine = typename boost::multi_index::index<OpenListMultiIdx_Score_Affine, IndexByFscore>::type;
    using OpenListByCoord_Score_Affine = typename boost::multi_index::index<OpenListMultiIdx_Score_Affine, IndexByCoord>::type;


    /* --- cost instead of score, fscore in ascending order --- */
    using OpenListMultiIdx_Cost_Affine = boost::multi_index_container<
        OpenListNodeAffineGap,                  // the data type stored
        boost::multi_index::indexed_by<             // list of indexes
            boost::multi_index::ordered_non_unique<
                boost::multi_index::tag<IndexByFscore>,
                boost::multi_index::member<
                    OpenListNodeAffineGap, 
                    STYPE, 
                    &OpenListNodeAffineGap::fscore
                >,
                std::less<STYPE>        // fscore in ascending order
            >,
            boost::multi_index::hashed_non_unique<
                boost::multi_index::tag<IndexByCoord>,
                boost::multi_index::member<OpenListNodeAffineGap, CoordTBInfo, &OpenListNodeAffineGap::crd_tb_info>, // what will be the index's key
                CoordTBInfoHash
            >
        >
    >;

    using OpenListByFscore_Cost_Affine = typename boost::multi_index::index<OpenListMultiIdx_Cost_Affine, IndexByFscore>::type;
    using OpenListByCoord_Cost_Affine = typename boost::multi_index::index<OpenListMultiIdx_Cost_Affine, IndexByCoord>::type;




    /* --- Linear --- */
    struct OpenListNodeLinearGap {
        std::array<CRDTYPE, NN> crd;
        DIRTYPE tb_info;        // trace back info
        STYPE fscore;
        STYPE gscore;

        std::array<CRDTYPE, NN> get_crd() const {return crd;}
        DIRTYPE get_tb_info() const {return tb_info;}
    };


    /**
     * This one not only changes the gscore but also modifies fscore accordingly.
     * fscore += new_gscore - node.gscore;
     * Usage: 
     * iter = indexByCoord.find(some_idx);
     * indexByCoord.modify(iter, ChangeGscore(new_gscore));
     */
    struct ChangeGscoreLinearGap {
        ChangeGscoreLinearGap(const STYPE& new_gscore):new_gscore(new_gscore){}
        
        void operator()(OpenListNodeLinearGap& node) {
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
    struct ChangeFscoreLinearGap {
        ChangeFscoreLinearGap(const STYPE& new_fscore):new_fscore(new_fscore){}
        
        void operator()(OpenListNodeLinearGap& node) {
            node.fscore = new_fscore;
        }

        private:
            STYPE new_fscore;
    };

    /**
     * Usage: 
     * iter = indexByCoord.find(some_idx);
     * indexByCoord.modify(iter, ChangeGscore(new_gscore));
     */
    struct ChangeGscoreFscoreTBinfoLinearGap {
        ChangeGscoreFscoreTBinfoLinearGap(const STYPE& new_gscore, const STYPE& new_fscore, const DIRTYPE& new_tbinfo):new_gscore(new_gscore), new_fscore(new_fscore), new_tbinfo(new_tbinfo){}
        
        void operator()(OpenListNodeLinearGap& node) {
            node.fscore = new_fscore;
            node.gscore = new_gscore;
            node.tb_info = new_tbinfo;
        }

        private:
            STYPE new_gscore;
            STYPE new_fscore;
            DIRTYPE new_tbinfo;
    };

    using OpenListMultiIdx_Score_Linear = boost::multi_index_container<
        OpenListNodeLinearGap,
        boost::multi_index::indexed_by<     // list of indexes
            boost::multi_index::ordered_non_unique<
                boost::multi_index::tag<IndexByFscore>,
                boost::multi_index::member<
                    OpenListNodeLinearGap, 
                    STYPE, 
                    &OpenListNodeLinearGap::fscore
                >,
                std::greater<STYPE>        // fscore in descending order
            >,
            boost::multi_index::hashed_non_unique<
                boost::multi_index::tag<IndexByCoord>,
                boost::multi_index::member<OpenListNodeLinearGap, std::array<CRDTYPE, NN>, &OpenListNodeLinearGap::crd> // what will be the index's key
            >
        >
    >;

    using OpenListByFscore_Score_Linear = typename boost::multi_index::index<OpenListMultiIdx_Score_Linear, IndexByFscore>::type;
    using OpenListByCoord_Score_Linear = typename boost::multi_index::index<OpenListMultiIdx_Score_Linear, IndexByCoord>::type;
  

    /* --- Tie breaker: None --- */
    using OpenListMultiIdx_Cost_Linear = boost::multi_index_container<
        OpenListNodeLinearGap,         // the data type stored, LINEAR_TB_TYPE for linear gap penalty
        boost::multi_index::indexed_by<             // list of indexes
            boost::multi_index::ordered_non_unique<
                boost::multi_index::tag<IndexByFscore>,
                boost::multi_index::member<
                    OpenListNodeLinearGap, 
                    STYPE, 
                    &OpenListNodeLinearGap::fscore
                >,
                std::less<STYPE>        // fscore in ascending order
            >,
            boost::multi_index::hashed_non_unique<
                boost::multi_index::tag<IndexByCoord>,
                boost::multi_index::member<OpenListNodeLinearGap, std::array<CRDTYPE, NN>, &OpenListNodeLinearGap::crd> // what will be the index's key
            >
        >
    >;


    using OpenListByFscore_Cost_Linear = typename boost::multi_index::index<OpenListMultiIdx_Cost_Linear, IndexByFscore>::type;
    using OpenListByCoord_Cost_Linear = typename boost::multi_index::index<OpenListMultiIdx_Cost_Linear, IndexByCoord>::type;


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