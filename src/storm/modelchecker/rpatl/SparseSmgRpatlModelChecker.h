#pragma once

#include "storm/modelchecker/propositional/SparsePropositionalModelChecker.h"

namespace storm {
namespace modelchecker {

template<class SparseSmgModelType>
class SparseSmgRpatlModelChecker : public SparsePropositionalModelChecker<SparseSmgModelType> {
   public:
    typedef typename SparseSmgModelType::ValueType ValueType;
    typedef typename SparseSmgModelType::RewardModelType RewardModelType;
    using SolutionType = typename std::conditional<std::is_same_v<ValueType, storm::Interval>, double, ValueType>::type;

    explicit SparseSmgRpatlModelChecker(SparseSmgModelType const& model);

    /*!
     * Returns false, if this task can certainly not be handled by this model checker (independent of the concrete model).
     * @param requiresSingleInitialState if not nullptr, this flag is set to true iff checking this formula requires a model with a single initial state
     */
    static bool canHandleStatic(CheckTask<storm::logic::Formula, ValueType> const& checkTask, bool* requiresSingleInitialState = nullptr);

    // The implemented methods of the AbstractModelChecker interface.
    virtual bool canHandle(CheckTask<storm::logic::Formula, ValueType> const& checkTask) const override;
    virtual std::unique_ptr<CheckResult> checkGameFormula(Environment const& env, CheckTask<storm::logic::GameFormula, ValueType> const& checkTask) override;

    virtual std::unique_ptr<CheckResult> computeLongRunAverageProbabilities(Environment const& env,
                                                                            CheckTask<storm::logic::StateFormula, ValueType> const& checkTask) override;
    virtual std::unique_ptr<CheckResult> computeLongRunAverageRewards(
        Environment const& env, CheckTask<storm::logic::LongRunAverageRewardFormula, ValueType> const& checkTask) override;

    virtual std::unique_ptr<CheckResult> computeUntilProbabilities(Environment const& env,
                                                                   CheckTask<storm::logic::UntilFormula, ValueType> const& checkTask) override;
};
}  // namespace modelchecker
}  // namespace storm
