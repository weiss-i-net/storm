#include "storm/modelchecker/rpatl/SparseSmgRpatlModelChecker.h"

#include <boost/format.hpp>
#include <memory>
#include <vector>

#include "storm/adapters/RationalNumberAdapter.h"
#include "storm/exceptions/InvalidPropertyException.h"
#include "storm/exceptions/NotImplementedException.h"
#include "storm/logic/FragmentSpecification.h"
#include "storm/modelchecker/helper/utility/SetInformationFromCheckTask.h"
#include "storm/modelchecker/results/ExplicitQualitativeCheckResult.h"
#include "storm/modelchecker/results/ExplicitQuantitativeCheckResult.h"
#include "storm/models/sparse/StandardRewardModel.h"
#include "storm/storage/BitVector.h"
#include "storm/utility/FilteredRewardModel.h"
#include "storm/utility/macros.h"

#include "storm/modelchecker/rpatl/helper/SparseSmgRpatlHelper.h"
#include "storm/models/sparse/Smg.h"

#include "storm/adapters/RationalNumberAdapter.h"

namespace storm::modelchecker {
template<typename SparseSmgModelType>
SparseSmgRpatlModelChecker<SparseSmgModelType>::SparseSmgRpatlModelChecker(SparseSmgModelType const& model)
    : SparsePropositionalModelChecker<SparseSmgModelType>(model) {
    // Intentionally left empty.
}

template<typename SparseSmgModelType>
bool SparseSmgRpatlModelChecker<SparseSmgModelType>::canHandleStatic(CheckTask<storm::logic::Formula, ValueType> const& checkTask,
                                                                     bool* requiresSingleInitialState) {
    storm::logic::Formula const& formula = checkTask.getFormula();

    return formula.isInFragment(storm::logic::rpatl());
}

template<typename SparseSmgModelType>
bool SparseSmgRpatlModelChecker<SparseSmgModelType>::canHandle(CheckTask<storm::logic::Formula, ValueType> const& checkTask) const {
    bool requiresSingleInitialState = false;
    if (canHandleStatic(checkTask, &requiresSingleInitialState)) {
        return !requiresSingleInitialState || this->getModel().getInitialStates().getNumberOfSetBits() == 1;
    } else {
        return false;
    }
}

template<typename SparseSmgModelType>
std::unique_ptr<CheckResult> SparseSmgRpatlModelChecker<SparseSmgModelType>::checkGameFormula(
    Environment const& env, CheckTask<storm::logic::GameFormula, ValueType> const& checkTask) {
    storm::logic::GameFormula const& gameFormula = checkTask.getFormula();
    storm::logic::Formula const& subFormula = gameFormula.getSubformula();
    if (subFormula.isOperatorFormula()) {
        return this->checkStateFormula(env, checkTask.substituteFormula(subFormula.asStateFormula()).setPlayerCoalition(gameFormula.getCoalition()));
    }
    STORM_LOG_THROW(false, storm::exceptions::NotImplementedException, "Only game formulas with Operatorformulas as subformula are supported.");
}

template<typename SparseSmgModelType>
std::unique_ptr<CheckResult> SparseSmgRpatlModelChecker<SparseSmgModelType>::computeLongRunAverageProbabilities(
    Environment const&, CheckTask<storm::logic::StateFormula, ValueType> const&) {
    STORM_LOG_THROW(false, storm::exceptions::NotImplementedException, "Not implemented.");
}

template<typename SparseSmgModelType>
std::unique_ptr<CheckResult> SparseSmgRpatlModelChecker<SparseSmgModelType>::computeLongRunAverageRewards(
    Environment const&, CheckTask<storm::logic::LongRunAverageRewardFormula, ValueType> const& checkTask) {
    auto rewardModel = storm::utility::createFilteredRewardModel(this->getModel(), checkTask);
    STORM_LOG_THROW(checkTask.isPlayerCoalitionSet(), storm::exceptions::InvalidPropertyException, "No player coalition was set.");
    auto coalitionStates = this->getModel().computeStatesOfCoalition(checkTask.getPlayerCoalition());
    std::cout << "Found " << coalitionStates.getNumberOfSetBits() << " states in coalition.\n";
    STORM_LOG_THROW(false, storm::exceptions::NotImplementedException, "Not implemented.");
}

template<typename SparseSmgModelType>
std::unique_ptr<CheckResult> SparseSmgRpatlModelChecker<SparseSmgModelType>::computeUntilProbabilities(
    Environment const& env, CheckTask<storm::logic::UntilFormula, ValueType> const& checkTask) {
    const storm::RationalNumber probablilityThreshold = checkTask.getBound().evaluateThresholdAsRational();
    STORM_LOG_THROW(probablilityThreshold == 0 or probablilityThreshold == 1, storm::exceptions::NotSupportedException,
                    boost::format("Quantitative rpatl is not supported. Probability threshold %1% should be either 0.0 or 1.0.") % probablilityThreshold);

    using helper_t = storm::modelchecker::helper::SparseSmgRpatlHelper<ValueType>;
    const auto& model = this->getModel();

    storm::logic::UntilFormula const& pathFormula = checkTask.getFormula();

    const auto leftResult = this->check(env, pathFormula.getLeftSubformula())->asExplicitQualitativeCheckResult();
    const auto rightResult = this->check(env, pathFormula.getRightSubformula())->asExplicitQualitativeCheckResult();

    const storage::BitVector playerCoalitionBitVector = model.computeStatesOfCoalition(checkTask.getPlayerCoalition());

    const auto prob1States = helper_t::computeUntilProp1(playerCoalitionBitVector, model.getTransitionMatrix(), model.getBackwardTransitions(),
                                                         leftResult.getTruthValuesVector(), rightResult.getTruthValuesVector());
    const auto prob0States = helper_t::computeUntilProp0(playerCoalitionBitVector, model.getTransitionMatrix(), model.getBackwardTransitions(),
                                                         leftResult.getTruthValuesVector(), rightResult.getTruthValuesVector());

    std::vector<ValueType> resultValues(model.getNumberOfStates(), 0.5);  // for uncertain states use 0.5 as default value
    for (const auto& index : prob1States) {
        resultValues.at(index) = 1;
    }
    for (const auto& index : prob0States) {
        STORM_LOG_ERROR_COND(resultValues.at(index) != 1, boost::format("State with index %1% has both probability 0 and 1.") % index);
        resultValues.at(index) = 0;
    }

    return std::make_unique<ExplicitQuantitativeCheckResult<SolutionType>>(std::move(resultValues));
}

template class SparseSmgRpatlModelChecker<storm::models::sparse::Smg<double>>;
template class SparseSmgRpatlModelChecker<storm::models::sparse::Smg<storm::RationalNumber>>;

}  // namespace storm::modelchecker
