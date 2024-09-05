//
// Created by janni on 30/07/2024.
//

#include "storm/modelchecker/rpatl/helper/SparseSmgRpatlHelper.h"

#include <boost/format.hpp>

#include "storm/exceptions/NotSupportedException.h"
#include "storm/storage/BitVector.h"
#include "storm/storage/SparseMatrix.h"

namespace storm::modelchecker::helper {

template<typename ValueType, typename SolutionType>
storm::storage::BitVector storm::modelchecker::helper::SparseSmgRpatlHelper<ValueType, SolutionType>::computeStrongAttractors(
    storm::storage::BitVector const& maximizerCoalition, storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
    storm::storage::SparseMatrix<ValueType> const& backwardTransitions, storm::storage::BitVector const& targetStateSet,
    storm::storage::BitVector const& allowedStateSet) {
    if constexpr (std::is_same_v<ValueType, storm::Interval> or std::is_same_v<ValueType, RationalNumber>) {
        STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Intervals and RationalNumbers not supported for SMGs.");
    } else {
        //STORM_LOG_ERROR(boost::format("maximizer: %1%, target: %2%, allowed: %3%") % maximizerCoalition.size() % targetStateSet.size() % allowedStateSet.size());
        auto stateSet{targetStateSet};
        std::vector stack(targetStateSet.begin(), targetStateSet.end());

        while (not stack.empty()) {
            const auto currentState = stack.back();
            stack.pop_back();

            for (const auto& predecessorEntry : backwardTransitions.getRow(currentState)) {
                const auto predecessor = predecessorEntry.getColumn();
                if (not allowedStateSet.get(predecessor))
                    break;
                if (stateSet.get(predecessor))
                    break;

                const auto nonDeterministicChoices = transitionMatrix.getRowGroupIndices(predecessor);

                const auto leadsToCurrentStateSet = [&](const auto& row_index) {
                    // it's enough when one of the successors is in stateSet as we tread the 'nature' states as belonging to the maximizer
                    const auto& transition_row = transitionMatrix.getRow(row_index);
                    return std::any_of(transition_row.begin(), transition_row.end(), [&](auto entry) { return stateSet.get(entry.getColumn()); });
                };

                // if the predecessor belongs to the maximizer one choice leading to currentStateSet is enough
                // if not, all need to lead to stateSet
                const bool predecessorIsValid = maximizerCoalition.get(predecessor)
                                                    ? std::any_of(nonDeterministicChoices.begin(), nonDeterministicChoices.end(), leadsToCurrentStateSet)
                                                    : std::all_of(nonDeterministicChoices.begin(), nonDeterministicChoices.end(), leadsToCurrentStateSet);

                if (predecessorIsValid) {
                    stateSet.set(predecessor, true);
                    stack.push_back(predecessor);
                }
            }
        }
        return stateSet;
    }
}
template<typename ValueType, typename SolutionType>
storm::storage::BitVector storm::modelchecker::helper::SparseSmgRpatlHelper<ValueType, SolutionType>::computeWeakAttractors(
    storm::storage::BitVector const& maximizerCoalition, storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
    storm::storage::SparseMatrix<ValueType> const& backwardTransitions, storm::storage::BitVector const& targetStateSet,
    storm::storage::BitVector const& allowedStateSet) {
    auto currentStateSet{targetStateSet};
    auto nextStateSet = computeStrongAttractors(maximizerCoalition, transitionMatrix, backwardTransitions, targetStateSet, allowedStateSet);

    while (currentStateSet != nextStateSet) {
        currentStateSet = std::move(nextStateSet);
        auto badStateSet =
            computeStrongAttractors(
                ~maximizerCoalition,
                transitionMatrix,
                backwardTransitions,
                ~currentStateSet,
                allowedStateSet & ~targetStateSet);
        nextStateSet = computeStrongAttractors(maximizerCoalition, transitionMatrix, backwardTransitions, targetStateSet, allowedStateSet & ~badStateSet);
    }
    return currentStateSet;
}

template<typename ValueType, typename SolutionType>
storm::storage::BitVector storm::modelchecker::helper::SparseSmgRpatlHelper<ValueType, SolutionType>::computeUntilProp1(
    storm::storage::BitVector const& maximizerCoalition, storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
    storm::storage::SparseMatrix<ValueType> const& backwardTransitions, storm::storage::BitVector const& phiStates,
    storm::storage::BitVector const& psiStates) {
    return computeWeakAttractors(maximizerCoalition, transitionMatrix, backwardTransitions, psiStates, phiStates);
}

template<typename ValueType, typename SolutionType>
storm::storage::BitVector storm::modelchecker::helper::SparseSmgRpatlHelper<ValueType, SolutionType>::computeUntilProp0(
    storm::storage::BitVector const& maximizerCoalition, storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
    storm::storage::SparseMatrix<ValueType> const& backwardTransitions, storm::storage::BitVector const& phiStates,
    storm::storage::BitVector const& psiStates) {
    return ~computeStrongAttractors(maximizerCoalition, transitionMatrix, backwardTransitions, psiStates, phiStates);;
}

template class SparseSmgRpatlHelper<double>;
template class SparseSmgRpatlHelper<storm::Interval, double>;
template class SparseSmgRpatlHelper<storm::RationalNumber>;

}  // namespace storm::modelchecker::helper
