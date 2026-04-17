#ifndef PARAMETERS_HPP
#define PARAMETERS_HPP

#define STYPE float    // score type
#define IDXTYPE size_t  // index type
#define CRDTYPE short   // coordinate type
#define DIRTYPE size_t  // direction type, each bit controls one dimension

// used in templates of anytime Astar MSA
#define LINEAR_TB_TYPE DIRTYPE

#define CHAR_NUM 20     // In BLOSUM62_double, exclude B, J, X, Z. But we use CHAR_NUM + 1 for X: average cost
#define GAP 45      // ASCII of '-'

#endif  // PARAMETERS_HPP