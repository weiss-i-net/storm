//
// Created by janni on 30/07/2024.
//

#include "storm/modelchecker/rpatl/helper/SparseSmgRpatlHelper.h"

#include <any>
#include <boost/format.hpp>
#include <ranges>

#include "storm/exceptions/NotSupportedException.h"
#include "storm/storage/BitVector.h"
#include "storm/storage/SparseMatrix.h"

namespace storm::modelchecker::helper {

/*template<typename ValueType, typename SolutionType>
storm::storage::BitVector storm::modelchecker::helper::SparseSmgRpatlHelper<ValueType, SolutionType>::computeStrongAttractors(
    storm::storage::BitVector const& maximizerCoalition, storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
    storm::storage::SparseMatrix<ValueType> const& backwardTransitions, storm::storage::BitVector const& targetStateSet,
    storm::storage::BitVector const& allowedStateSet, storm::storage::BitVector const& allowedTransitions) {
    if constexpr (std::is_same_v<ValueType, storm::Interval> or std::is_same_v<ValueType, RationalNumber>) {
        STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Intervals and RationalNumbers not supported for SMGs.");
    } else {
        //STORM_LOG_ERROR(boost::format("maximizer: %1%, target: %2%, allowed: %3%") % maximizerCoalition.size() % targetStateSet.size() %
allowedStateSet.size()); storm::storage::BitVector attractors(transitionMatrix.getColumnCount()); storm::storage::BitVector
usedTransitions(transitionMatrix.getRowCount());

        std::vector stack(targetStateSet.begin(), targetStateSet.end());

        while (not stack.empty()) {
            const auto currentState = stack.back();
            stack.pop_back();

            for (const auto& predecessorEntry : backwardTransitions.getRow(currentState)) {
                const auto predecessor = predecessorEntry.getColumn();

                const auto nonDeterministicChoices = transitionMatrix.getRowGroupIndices(predecessor);

                const auto leadsToCurrentStateSet = [&](const auto& row_index) {
                    // it's enough when one of the successors is in stateSet as we tread the 'nature' states as belonging to the maximizer
                    const auto& transition_row = transitionMatrix.getRow(row_index);
                    return std::any_of(transition_row.begin(), transition_row.end(), [&](auto entry) { return attractors.get(entry.getColumn()); });
                };

                const auto markTransitionsAsUsed = [&usedTransitions](const auto& row_index) {
                    usedTransitions.set(row_index);
                    return row_index;
                };

                // if the predecessor belongs to the maximizer one choice leading to currentStateSet is enough
                // if not, all need to lead to stateSet
                const bool predecessorIsValid = maximizerCoalition.get(predecessor)
                                                    ? std::any_of(nonDeterministicChoices.begin(), nonDeterministicChoices.end(), leadsToCurrentStateSet)
                                                    : std::all_of(nonDeterministicChoices.begin(), nonDeterministicChoices.end(), leadsToCurrentStateSet);

                if (predecessorIsValid) {
                    attractors.set(predecessor, true);
                    stack.push_back(predecessor);
                }
            }
        }
        return attractors;
    }
}*/

template<>
std::pair<storm::storage::BitVector, storm::storage::BitVector>
storm::modelchecker::helper::SparseSmgRpatlHelper<storm::Interval, double>::computeStrongAttractors(
    storm::storage::BitVector const& maximizerCoalition, storm::storage::SparseMatrix<storm::Interval> const& transitionMatrix,
    storm::storage::SparseMatrix<storm::Interval> const& backwardTransitions, storm::storage::BitVector const& targetStateSet,
    storm::storage::BitVector const& allowedStateSet, storm::storage::BitVector const& allowedTransitions) {
    STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "Intervals not supported for SMGs.");
    return {};
}


template<>
std::pair<storm::storage::BitVector, storm::storage::BitVector>
storm::modelchecker::helper::SparseSmgRpatlHelper<storm::RationalNumber>::computeStrongAttractors(
    storm::storage::BitVector const& maximizerCoalition, storm::storage::SparseMatrix<storm::RationalNumber> const& transitionMatrix,
    storm::storage::SparseMatrix<storm::RationalNumber> const& backwardTransitions, storm::storage::BitVector const& targetStateSet,
    storm::storage::BitVector const& allowedStateSet, storm::storage::BitVector const& allowedTransitions) {
    STORM_LOG_THROW(false, storm::exceptions::NotSupportedException, "RationalNumbers not supported for SMGs.");
    return {};
}

template<typename ValueType, typename SolutionType>
std::pair<storm::storage::BitVector, storm::storage::BitVector>
storm::modelchecker::helper::SparseSmgRpatlHelper<ValueType, SolutionType>::computeStrongAttractors(
    storm::storage::BitVector const& maximizerCoalition, storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
    storm::storage::SparseMatrix<ValueType> const& backwardTransitions, storm::storage::BitVector const& targetStateSet,
    storm::storage::BitVector const& allowedStateSet, storm::storage::BitVector const& allowedTransitions) {
    /* The algorithm works by treating random choices as nature states owned by the maximizer coalition and then exploring backwards from the targetStateSet,
     * adding states if they belong to the coalition and have one successor leading to the current attractors or otherwise if all successors lead to the current
     * attractors. We go backwards in layers, with each workingStateSet being a new layer consisting of states that are backwards reachable from the previous.
     */
    storm::storage::BitVector strongAttractors(targetStateSet);
    storm::storage::BitVector transitionsLeadingToAttractors(transitionMatrix.getRowCount());
    storm::storage::BitVector workingStateSet{targetStateSet};

    //STORM_PRINT("matrix: " << transitionMatrix << "\n");

    while (not workingStateSet.empty()) {
        //STORM_PRINT("attractors: " << strongAttractors << "\n");
        //STORM_PRINT("working: " << workingStateSet << "\n");

        auto predecessorStateView =
            workingStateSet |
            std::views::transform([&](auto const index) { return backwardTransitions.getRow(index); })  // get predecessor states for each state
            | std::views::join                                                                          // flatten
            | std::views::transform([](auto const entry) { return entry.getColumn(); })                 // get indexes
            | std::views::filter([&](auto const index) { return not strongAttractors.get(index); });    // already explored

        // we cant use the constructor with iterators because the temlplate is only instantiated for some specific iterators.
        storm::storage::BitVector predecessorStates(transitionMatrix.getColumnCount());
        for (auto const index : predecessorStateView) {
            predecessorStates.set(index);
        };
        //STORM_PRINT("predecessors: " << predecessorStates << "\n");

        auto transitionLeadsToAttractors = [&](auto const index) {
            return std::ranges::any_of(transitionMatrix.getRow(index), [&](auto const entry) { return strongAttractors.get(entry.getColumn()); });
        };
        auto transitionIndexView =
            predecessorStates | std::views::transform([&](auto const index) { return transitionMatrix.getRowGroupIndices(index); })  // get transition indexes
            | std::views::join                                                                                                       // flatten
            | std::views::filter([&](auto const index) { return allowedTransitions.get(index); })  // only use allowed transitions
            | std::views::filter(transitionLeadsToAttractors);                                     // only use transitions that lead into attractors

        for (auto const index : transitionIndexView) {
            transitionsLeadingToAttractors.set(index);
        }
        //STORM_PRINT("transitions: " << transitionsLeadingToAttractors << "\n\n");

        auto isGoodState = [&](auto const stateIndex) {
            auto isUsed = [&](auto const rowIndex) { return transitionsLeadingToAttractors.get(rowIndex); };
            auto outGoingTransitions = transitionMatrix.getRowGroupIndices(stateIndex);
            if (maximizerCoalition.get(stateIndex)) {
                return std::ranges::any_of(outGoingTransitions, isUsed);
            } else {
                return std::ranges::all_of(outGoingTransitions, isUsed);
            }
        };
        auto notAlreadyAttractor = [&](auto const index) {
            return not strongAttractors.get(index);
        };

        predecessorStates &= allowedStateSet;
        workingStateSet.clear();
        for (auto const index : predecessorStates | std::views::filter(isGoodState) | std::views::filter(notAlreadyAttractor)) {
            workingStateSet.set(index);
        }
        strongAttractors |= workingStateSet;
    }

    return std::pair{strongAttractors, transitionsLeadingToAttractors};
}

template<typename ValueType, typename SolutionType>
storm::storage::BitVector storm::modelchecker::helper::SparseSmgRpatlHelper<ValueType, SolutionType>::computeWeakAttractors(
    storm::storage::BitVector const& maximizerCoalition, storm::storage::SparseMatrix<ValueType> const& transitionMatrix,
    storm::storage::SparseMatrix<ValueType> const& backwardTransitions, storm::storage::BitVector const& targetStateSet,
    storm::storage::BitVector const& allowedStateSet) {
    auto currentStateSet{targetStateSet};
    storm::storage::BitVector currentTransitions(transitionMatrix.getRowCount());
    auto [nextStateSet, nextTransitions] =
        computeStrongAttractors(maximizerCoalition, transitionMatrix, backwardTransitions, targetStateSet, allowedStateSet, ~currentTransitions);

    while (currentStateSet != nextStateSet) {
        currentStateSet = std::move(nextStateSet);
        auto [badStateSet, badTransitions] = computeStrongAttractors(~maximizerCoalition, transitionMatrix, backwardTransitions, ~currentStateSet,
                                                                     allowedStateSet & ~targetStateSet, ~currentTransitions);
        std::tie(nextStateSet, nextTransitions) =
            computeStrongAttractors(maximizerCoalition, transitionMatrix, backwardTransitions, targetStateSet, allowedStateSet & ~badStateSet, ~badTransitions);
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
    auto [states, transitions] = computeStrongAttractors(maximizerCoalition, transitionMatrix, backwardTransitions, psiStates, phiStates,
                                                         storm::storage::BitVector(transitionMatrix.getRowCount(), true));
    return ~states;
}

template class SparseSmgRpatlHelper<double>;
template class SparseSmgRpatlHelper<storm::Interval, double>;
template class SparseSmgRpatlHelper<storm::RationalNumber>;

}  // namespace storm::modelchecker::helper
