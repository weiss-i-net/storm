//
// Created by Jannik Hiller on 06/09/2023.
//

#ifndef STORM_SPARSESMGRPATLHELPER_H
#define STORM_SPARSESMGRPATLHELPER_H
//
// Created by Jannik Hiller on 16/08/2023.
//

#include <vector>
#include "solver/SolveGoal.h"
#include "storage/BitVector.h"

template<typename ValueType>
storm::storage::BitVector computeStrongAttractors(storm::storage::BitVector const& maximizerCoalition,
                                                  storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
                                                  storm::storage::SparseMatrix<ValueType> const& backwardTransitions,
                                                  storm::storage::BitVector const& targetStates, storm::storage::BitVector const& allowedStates) {
    auto stateSet{targetStates};
    std::vector<uint_fast64_t> stack{targetStates.begin(), targetStates.end()};
    while (not stack.empty()) {
        const auto currentState = stack.back();
        stack.pop_back();

        for (const auto& predecessorEntry : backwardTransitions.getRow(currentState)) {
            const auto predecessor = predecessorEntry.getColumn();
            if (not allowedStates.get(predecessor))
                break;
            if (stateSet.get(predecessor))
                break;

            const auto nonDeterministicChoices = transitionMatrix.getRowGroup(predecessor);
            const auto leadsToCurrentStateSet = [&](const auto& transition_row) {
                // it's enough when one of the successors is in stateSet as we tread the 'nature' states as belonging to the maximizer
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

template<typename ValueType>
auto computeWeakAttractors(storm::storage::BitVector const& maximizerCoalition, storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
                           storm::storage::SparseMatrix<ValueType> const& backwardTransitions, storm::storage::BitVector const& targetStateSet,
                           storm::storage::BitVector const& allowedStateSet) -> storm::storage::BitVector {
    auto currentStateSet{targetStateSet};
    auto nextStateSet = computeStrongAttractors(maximizerCoalition, transitionMatrix, backwardTransitions, targetStateSet, allowedStateSet);

    while (currentStateSet != nextStateSet) {
        currentStateSet = std::move(nextStateSet);
        auto badStateSet =
            computeStrongAttractors(~maximizerCoalition, transitionMatrix, backwardTransitions, ~currentStateSet, allowedStateSet % ~targetStateSet);
        nextStateSet = computeStrongAttractors(maximizerCoalition, transitionMatrix, backwardTransitions, targetStateSet, allowedStateSet % ~badStateSet);
    }
    return currentStateSet;
}

#endif  // STORM_SPARSESMGRPATLHELPER_H
