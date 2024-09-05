//
// Created by Jannik Hiller on 06/09/2023.
//

#ifndef STORM_SPARSESMGRPATLHELPER_H
#define STORM_SPARSESMGRPATLHELPER_H
//
// Created by Jannik Hiller on 16/08/2023.
//

#include "storm/solver/SolveGoal.h"
#include "storm/storage/BitVector.h"

namespace storm::modelchecker::helper {

template<typename ValueType, typename SolutionType = ValueType>
class SparseSmgRpatlHelper {
   private:
    static storm::storage::BitVector computeStrongAttractors(storm::storage::BitVector const& maximizerCoalition,
                                                             storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
                                                             storm::storage::SparseMatrix<ValueType> const& backwardTransitions,
                                                             storm::storage::BitVector const& targetStateSet, storm::storage::BitVector const& allowedStateSet);

    static storm::storage::BitVector computeWeakAttractors(storm::storage::BitVector const& maximizerCoalition,
                                                           storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
                                                           storm::storage::SparseMatrix<ValueType> const& backwardTransitions,
                                                           storm::storage::BitVector const& targetStateSet, storm::storage::BitVector const& allowedStateSet);

   public:
    static storm::storage::BitVector computeUntilProp1(storm::storage::BitVector const& maximizerCoalition,
                                                       storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
                                                       storm::storage::SparseMatrix<ValueType> const& backwardTransitions,
                                                       storm::storage::BitVector const& phiStates, storm::storage::BitVector const& psiStates);
    static storm::storage::BitVector computeUntilProp0(storm::storage::BitVector const& maximizerCoalition,
                                                       storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
                                                       storm::storage::SparseMatrix<ValueType> const& backwardTransitions,
                                                       storm::storage::BitVector const& phiStates, storm::storage::BitVector const& psiStates);
};

}  // namespace storm::modelchecker::helper

#endif  // STORM_SPARSESMGRPATLHELPER_H
