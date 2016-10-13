#include "src/builder/ExplicitDFTModelBuilderApprox.h"
#include <src/models/sparse/MarkovAutomaton.h>
#include <src/models/sparse/Ctmc.h>
#include <src/utility/constants.h>
#include <src/utility/vector.h>
#include <src/exceptions/UnexpectedException.h>
#include "src/settings/modules/DFTSettings.h"
#include "src/settings/SettingsManager.h"
#include <map>

namespace storm {
    namespace builder {

        template<typename ValueType, typename StateType>
        ExplicitDFTModelBuilderApprox<ValueType, StateType>::ModelComponents::ModelComponents() : transitionMatrix(), stateLabeling(), markovianStates(), exitRates(), choiceLabeling() {
            // Intentionally left empty.
        }

        template<typename ValueType, typename StateType>
        ExplicitDFTModelBuilderApprox<ValueType, StateType>::MatrixBuilder::MatrixBuilder(bool canHaveNondeterminism) : mappingOffset(0), stateRemapping(), currentRowGroup(0), currentRow(0), canHaveNondeterminism((canHaveNondeterminism)) {
            // Create matrix builder
            builder = storm::storage::SparseMatrixBuilder<ValueType>(0, 0, 0, false, canHaveNondeterminism, 0);
        }

        template<typename ValueType, typename StateType>
        ExplicitDFTModelBuilderApprox<ValueType, StateType>::ExplicitDFTModelBuilderApprox(storm::storage::DFT<ValueType> const& dft, storm::storage::DFTIndependentSymmetries const& symmetries, bool enableDC) :
                dft(dft),
                stateGenerationInfo(std::make_shared<storm::storage::DFTStateGenerationInfo>(dft.buildStateGenerationInfo(symmetries))),
                enableDC(enableDC),
                generator(dft, *stateGenerationInfo, enableDC, mergeFailedStates),
                matrixBuilder(!generator.isDeterministicModel()),
                stateStorage(((dft.stateVectorSize() / 64) + 1) * 64),
                statesToExplore(storm::storage::DynamicPriorityQueue<DFTStatePointer, std::vector<DFTStatePointer>, std::function<bool(DFTStatePointer, DFTStatePointer)>>(&storm::builder::compareDepth<ValueType>))
        {
            // stateVectorSize is bound for size of bitvector
            heuristic = storm::settings::getModule<storm::settings::modules::DFTSettings>().getApproximationHeuristic();
        }

        template<typename ValueType, typename StateType>
        void ExplicitDFTModelBuilderApprox<ValueType, StateType>::buildModel(LabelOptions const& labelOpts, bool firstTime, double approximationThreshold) {
            STORM_LOG_TRACE("Generating DFT state space");

            if (firstTime) {
                // Initialize
                modelComponents.markovianStates = storm::storage::BitVector(INITIAL_BITVECTOR_SIZE);

                if(mergeFailedStates) {
                    // Introduce explicit fail state
                    storm::generator::StateBehavior<ValueType, StateType> behavior = generator.createMergeFailedState([this] (DFTStatePointer const& state) {
                        this->failedStateId = newIndex++;
                        matrixBuilder.stateRemapping.push_back(0);
                        return this->failedStateId;
                    } );

                    matrixBuilder.setRemapping(failedStateId);
                    STORM_LOG_ASSERT(!behavior.empty(), "Behavior is empty.");
                    matrixBuilder.newRowGroup();
                    setMarkovian(behavior.begin()->isMarkovian());

                    // Now add self loop.
                    // TODO Matthias: maybe use general method.
                    STORM_LOG_ASSERT(behavior.getNumberOfChoices() == 1, "Wrong number of choices for failed state.");
                    STORM_LOG_ASSERT(behavior.begin()->size() == 1, "Wrong number of transitions for failed state.");
                    std::pair<StateType, ValueType> stateProbabilityPair = *(behavior.begin()->begin());
                    STORM_LOG_ASSERT(stateProbabilityPair.first == failedStateId, "No self loop for failed state.");
                    STORM_LOG_ASSERT(storm::utility::isOne<ValueType>(stateProbabilityPair.second), "Probability for failed state != 1.");
                    matrixBuilder.addTransition(stateProbabilityPair.first, stateProbabilityPair.second);
                    matrixBuilder.finishRow();
                }

                // Build initial state
                this->stateStorage.initialStateIndices = generator.getInitialStates(std::bind(&ExplicitDFTModelBuilderApprox::getOrAddStateIndex, this, std::placeholders::_1));
                STORM_LOG_ASSERT(stateStorage.initialStateIndices.size() == 1, "Only one initial state assumed.");
                initialStateIndex = stateStorage.initialStateIndices[0];
                STORM_LOG_TRACE("Initial state: " << initialStateIndex);

            } else {
                initializeNextIteration();
            }

            exploreStateSpace(approximationThreshold);

            size_t stateSize = stateStorage.getNumberOfStates() + (mergeFailedStates ? 1 : 0);
            modelComponents.markovianStates.resize(stateSize);
            modelComponents.deterministicModel = generator.isDeterministicModel();

            // Replace pseudo states in matrix
            if (!mPseudoStatesMapping.empty()) {
                // TODO Matthias: avoid hack with fixed int type
                std::vector<uint_fast64_t> pseudoStatesVector;
                for (auto const& pseudoStatePair : mPseudoStatesMapping) {
                    pseudoStatesVector.push_back(matrixBuilder.mappingOffset + pseudoStatePair.first);
                }
                STORM_LOG_ASSERT(std::find(pseudoStatesVector.begin(), pseudoStatesVector.end(), 0) == pseudoStatesVector.end(), "Unexplored pseudo state still contained.");
                STORM_LOG_TRACE("Replace pseudo states: " << pseudoStatesVector << ", offset: " << OFFSET_PSEUDO_STATE);
                // TODO Matthias: combine replacement with later one
                matrixBuilder.builder.replaceColumns(pseudoStatesVector, OFFSET_PSEUDO_STATE);
                mPseudoStatesMapping.clear();
            }

            // Fix the entries in the transition matrix according to the mapping of ids to row group indices
            STORM_LOG_ASSERT(matrixBuilder.stateRemapping[initialStateIndex] == initialStateIndex, "Initial state should not be remapped.");
            // TODO Matthias: do not consider all rows?
            STORM_LOG_TRACE("Remap matrix: " << matrixBuilder.stateRemapping << ", offset: " << matrixBuilder.mappingOffset);
            matrixBuilder.remap();

            STORM_LOG_TRACE("State remapping: " << matrixBuilder.stateRemapping);
            STORM_LOG_TRACE("Markovian states: " << modelComponents.markovianStates);
            STORM_LOG_DEBUG("Generated " << stateSize << " states");
            STORM_LOG_DEBUG("Skipped " << skippedStates.size() << " states");
            STORM_LOG_DEBUG("Model is " << (generator.isDeterministicModel() ? "deterministic" : "non-deterministic"));

            // Build transition matrix
            modelComponents.transitionMatrix = matrixBuilder.builder.build(stateSize, stateSize);
            if (stateSize <= 15) {
                STORM_LOG_TRACE("Transition matrix: " << std::endl << modelComponents.transitionMatrix);
            } else {
                STORM_LOG_TRACE("Transition matrix: too big to print");
            }

            buildLabeling(labelOpts);
        }

        template<typename ValueType, typename StateType>
        void ExplicitDFTModelBuilderApprox<ValueType, StateType>::initializeNextIteration() {
            STORM_LOG_TRACE("Refining DFT state space");

            // Initialize matrix builder again
            // TODO Matthias: avoid copy
            std::vector<uint_fast64_t> copyRemapping = matrixBuilder.stateRemapping;
            matrixBuilder = MatrixBuilder(!generator.isDeterministicModel());
            matrixBuilder.stateRemapping = copyRemapping;
            StateType nrStates = modelComponents.transitionMatrix.getRowGroupCount();
            STORM_LOG_ASSERT(nrStates == matrixBuilder.stateRemapping.size(), "No. of states does not coincide with mapping size.");

            // Start by creating a remapping from the old indices to the new indices
            std::vector<StateType> indexRemapping = std::vector<StateType>(nrStates, 0);
            auto iterSkipped = skippedStates.begin();
            size_t skippedBefore = 0;
            for (size_t i = 0; i < indexRemapping.size(); ++i) {
                while (iterSkipped->first <= i) {
                    ++skippedBefore;
                    ++iterSkipped;
                }
                indexRemapping[i] = i - skippedBefore;
            }

            // Set remapping
            size_t nrExpandedStates = nrStates - skippedBefore;
            matrixBuilder.mappingOffset = nrStates;
            STORM_LOG_TRACE("# expanded states: " << nrExpandedStates);
            StateType skippedIndex = nrExpandedStates;
            std::map<StateType, DFTStatePointer> skippedStatesNew;
            for (size_t id = 0; id < matrixBuilder.stateRemapping.size(); ++id) {
                StateType index = matrixBuilder.stateRemapping[id];
                auto itFind = skippedStates.find(index);
                if (itFind != skippedStates.end()) {
                    // Set new mapping for skipped state
                    matrixBuilder.stateRemapping[id] = skippedIndex;
                    skippedStatesNew[skippedIndex] = itFind->second;
                    indexRemapping[index] = skippedIndex;
                    ++skippedIndex;
                } else {
                    // Set new mapping for expanded state
                    matrixBuilder.stateRemapping[id] = indexRemapping[index];
                }
            }
            STORM_LOG_TRACE("New state remapping: " << matrixBuilder.stateRemapping);
            std::stringstream ss;
            ss << "Index remapping:" << std::endl;
            for (auto tmp : indexRemapping) {
                ss << tmp << " ";
            }
            STORM_LOG_TRACE(ss.str());

            // Remap markovian states
            storm::storage::BitVector markovianStatesNew = storm::storage::BitVector(modelComponents.markovianStates.size(), true);
            // Iterate over all not set bits
            modelComponents.markovianStates.complement();
            size_t index = modelComponents.markovianStates.getNextSetIndex(0);
            while (index < modelComponents.markovianStates.size()) {
                markovianStatesNew.set(indexRemapping[index], false);
                index = modelComponents.markovianStates.getNextSetIndex(index+1);
            }
            STORM_LOG_ASSERT(modelComponents.markovianStates.size() - modelComponents.markovianStates.getNumberOfSetBits() == markovianStatesNew.getNumberOfSetBits(), "Remapping of markovian states is wrong.");
            STORM_LOG_ASSERT(markovianStatesNew.size() == nrStates, "No. of states does not coincide with markovian size.");
            modelComponents.markovianStates = markovianStatesNew;

            // Build submatrix for expanded states
            // TODO Matthias: only use row groups when necessary
            for (StateType oldRowGroup = 0; oldRowGroup < modelComponents.transitionMatrix.getRowGroupCount(); ++oldRowGroup) {
                if (indexRemapping[oldRowGroup] < nrExpandedStates) {
                    // State is expanded -> copy to new matrix
                    matrixBuilder.newRowGroup();
                    for (StateType oldRow = modelComponents.transitionMatrix.getRowGroupIndices()[oldRowGroup]; oldRow < modelComponents.transitionMatrix.getRowGroupIndices()[oldRowGroup+1]; ++oldRow) {
                        for (typename storm::storage::SparseMatrix<ValueType>::const_iterator itEntry = modelComponents.transitionMatrix.begin(oldRow); itEntry != modelComponents.transitionMatrix.end(oldRow); ++itEntry) {
                            auto itFind = skippedStates.find(itEntry->getColumn());
                            if (itFind != skippedStates.end()) {
                                // Set id for skipped states as we remap it later
                                matrixBuilder.addTransition(matrixBuilder.mappingOffset + itFind->second->getId(), itEntry->getValue());
                            } else {
                                // Set newly remapped index for expanded states
                                matrixBuilder.addTransition(indexRemapping[itEntry->getColumn()], itEntry->getValue());
                            }
                        }
                        matrixBuilder.finishRow();
                    }
                }
            }

            skippedStates = skippedStatesNew;

            STORM_LOG_ASSERT(matrixBuilder.getCurrentRowGroup() == nrExpandedStates, "Row group size does not match.");

            // Push skipped states to explore queue
            // TODO Matthias: remove
            for (auto const& skippedState : skippedStates) {
                statesToExplore.push(skippedState.second);
            }
            skippedStates.clear();
        }

        template<typename ValueType, typename StateType>
        void ExplicitDFTModelBuilderApprox<ValueType, StateType>::exploreStateSpace(double approximationThreshold) {
            bool explorationFinished = false;
            size_t pseudoStatesToCheck = 0;
            while (!explorationFinished) {
                // Get the first state in the queue
                DFTStatePointer currentState = statesToExplore.top();
                STORM_LOG_ASSERT(stateStorage.stateToId.contains(currentState->status()), "State is not contained in state storage.");
                STORM_LOG_ASSERT(stateStorage.stateToId.getValue(currentState->status()) == currentState->getId(), "Ids of states do not coincide.");
                statesToExplore.pop();

                // Remember that the current row group was actually filled with the transitions of a different state
                matrixBuilder.setRemapping(currentState->getId());

                matrixBuilder.newRowGroup();

                // Try to explore the next state
                generator.load(currentState);

                if (currentState->isSkip(approximationThreshold, heuristic)) {
                    // Skip the current state
                    STORM_LOG_TRACE("Skip expansion of state: " << dft.getStateString(currentState));
                    setMarkovian(true);
                    // Add transition to target state with temporary value 0
                    // TODO Matthias: what to do when there is no unique target state?
                    matrixBuilder.addTransition(failedStateId, storm::utility::zero<ValueType>());
                    // Remember skipped state
                    skippedStates[matrixBuilder.getCurrentRowGroup() - 1] = currentState;
                    matrixBuilder.finishRow();
                } else {
                    // Explore the current state
                    storm::generator::StateBehavior<ValueType, StateType> behavior = generator.expand(std::bind(&ExplicitDFTModelBuilderApprox::getOrAddStateIndex, this, std::placeholders::_1));
                    STORM_LOG_ASSERT(!behavior.empty(), "Behavior is empty.");
                    setMarkovian(behavior.begin()->isMarkovian());

                    // Now add all choices.
                    for (auto const& choice : behavior) {
                        // Add the probabilistic behavior to the matrix.
                        for (auto const& stateProbabilityPair : choice) {
                            if (stateProbabilityPair.first >= OFFSET_PSEUDO_STATE) {
                                // Check that pseudo state and its instantiation do not appear together
                                // TODO Matthias: prove that this is not possible and remove
                                StateType newId = stateProbabilityPair.first - OFFSET_PSEUDO_STATE;
                                STORM_LOG_ASSERT(newId < mPseudoStatesMapping.size(), "Id is not valid.");
                                if (mPseudoStatesMapping[newId].first > 0) {
                                    // State exists already
                                    newId = mPseudoStatesMapping[newId].first;
                                    for (auto itFind = choice.begin(); itFind != choice.end(); ++itFind) {
                                        STORM_LOG_ASSERT(itFind->first != newId, "Pseudo state and instantiation occur together in a distribution.");
                                    }
                                }
                                // Set transtion to pseudo state
                                matrixBuilder.addTransition(stateProbabilityPair.first, stateProbabilityPair.second);
                            } else {
                                // Set transition to state id + offset. This helps in only remapping all previously skipped states.
                                matrixBuilder.addTransition(matrixBuilder.mappingOffset + stateProbabilityPair.first, stateProbabilityPair.second);
                            }
                        }
                        matrixBuilder.finishRow();
                    }
                }

                // Update priority queue
                statesToExplore.fix();

                if (statesToExplore.empty()) {
                    explorationFinished = true;
                    // Before ending the exploration check for pseudo states which are not initialized yet
                    for ( ; pseudoStatesToCheck < mPseudoStatesMapping.size(); ++pseudoStatesToCheck) {
                        std::pair<StateType, storm::storage::BitVector> pseudoStatePair = mPseudoStatesMapping[pseudoStatesToCheck];
                        if (pseudoStatePair.first == 0) {
                            // Create state from pseudo state and explore
                            STORM_LOG_ASSERT(stateStorage.stateToId.contains(pseudoStatePair.second), "Pseudo state not contained.");
                            STORM_LOG_ASSERT(stateStorage.stateToId.getValue(pseudoStatePair.second) >= OFFSET_PSEUDO_STATE, "State is no pseudo state.");
                            STORM_LOG_TRACE("Create pseudo state from bit vector " << pseudoStatePair.second);
                            DFTStatePointer pseudoState = std::make_shared<storm::storage::DFTState<ValueType>>(pseudoStatePair.second, dft, *stateGenerationInfo, newIndex);
                            STORM_LOG_ASSERT(pseudoStatePair.second == pseudoState->status(), "Pseudo states do not coincide.");
                            STORM_LOG_TRACE("Explore pseudo state " << dft.getStateString(pseudoState) << " with id " << pseudoState->getId());

                            getOrAddStateIndex(pseudoState);
                            explorationFinished = false;
                            break;
                        }
                    }
                }

            } // end exploration
        }

        template<typename ValueType, typename StateType>
        void ExplicitDFTModelBuilderApprox<ValueType, StateType>::buildLabeling(LabelOptions const& labelOpts) {
            // Build state labeling
            modelComponents.stateLabeling = storm::models::sparse::StateLabeling(modelComponents.transitionMatrix.getRowGroupCount());
            // Initial state
            modelComponents.stateLabeling.addLabel("init");
            modelComponents.stateLabeling.addLabelToState("init", initialStateIndex);
            // Label all states corresponding to their status (failed, failsafe, failed BE)
            if(labelOpts.buildFailLabel) {
                modelComponents.stateLabeling.addLabel("failed");
            }
            if(labelOpts.buildFailSafeLabel) {
                modelComponents.stateLabeling.addLabel("failsafe");
            }

            // Collect labels for all BE
            std::vector<std::shared_ptr<storage::DFTBE<ValueType>>> basicElements = dft.getBasicElements();
            for (std::shared_ptr<storage::DFTBE<ValueType>> elem : basicElements) {
                if(labelOpts.beLabels.count(elem->name()) > 0) {
                    modelComponents.stateLabeling.addLabel(elem->name() + "_fail");
                }
            }

            // Set labels to states
            if(mergeFailedStates) {
                modelComponents.stateLabeling.addLabelToState("failed", failedStateId);
            }
            for (auto const& stateIdPair : stateStorage.stateToId) {
                storm::storage::BitVector state = stateIdPair.first;
                size_t stateId = stateIdPair.second;
                if (!mergeFailedStates && labelOpts.buildFailLabel && dft.hasFailed(state, *stateGenerationInfo)) {
                    modelComponents.stateLabeling.addLabelToState("failed", stateId);
                }
                if (labelOpts.buildFailSafeLabel && dft.isFailsafe(state, *stateGenerationInfo)) {
                    modelComponents.stateLabeling.addLabelToState("failsafe", stateId);
                };
                // Set fail status for each BE
                for (std::shared_ptr<storage::DFTBE<ValueType>> elem : basicElements) {
                    if (labelOpts.beLabels.count(elem->name()) > 0 && storm::storage::DFTState<ValueType>::hasFailed(state, stateGenerationInfo->getStateIndex(elem->id())) ) {
                        modelComponents.stateLabeling.addLabelToState(elem->name() + "_fail", stateId);
                    }
                }
            }
        }

        template<typename ValueType, typename StateType>
        std::shared_ptr<storm::models::sparse::Model<ValueType>> ExplicitDFTModelBuilderApprox<ValueType, StateType>::getModel() {
            STORM_LOG_ASSERT(skippedStates.size() == 0, "Concrete model has skipped states");
            return createModel(false);
        }

        template<typename ValueType, typename StateType>
        std::shared_ptr<storm::models::sparse::Model<ValueType>> ExplicitDFTModelBuilderApprox<ValueType, StateType>::getModelApproximation(bool lowerBound) {
            // TODO Matthias: handle case with no skipped states
            if (lowerBound) {
                changeMatrixLowerBound(modelComponents.transitionMatrix);
            } else {
                changeMatrixUpperBound(modelComponents.transitionMatrix);
            }
            return createModel(true);
        }

        template<typename ValueType, typename StateType>
        std::shared_ptr<storm::models::sparse::Model<ValueType>> ExplicitDFTModelBuilderApprox<ValueType, StateType>::createModel(bool copy) {
            std::shared_ptr<storm::models::sparse::Model<ValueType>> model;

            if (modelComponents.deterministicModel) {
                // Build CTMC
                if (copy) {
                    model = std::make_shared<storm::models::sparse::Ctmc<ValueType>>(modelComponents.transitionMatrix, modelComponents.stateLabeling);
                } else {
                    model = std::make_shared<storm::models::sparse::Ctmc<ValueType>>(std::move(modelComponents.transitionMatrix), std::move(modelComponents.stateLabeling));
                }
            } else {
                // Build MA
                // Compute exit rates
                // TODO Matthias: avoid computing multiple times
                modelComponents.exitRates = std::vector<ValueType>(modelComponents.markovianStates.size());
                std::vector<typename storm::storage::SparseMatrix<ValueType>::index_type> indices = modelComponents.transitionMatrix.getRowGroupIndices();
                for (StateType stateIndex = 0; stateIndex < modelComponents.markovianStates.size(); ++stateIndex) {
                    if (modelComponents.markovianStates[stateIndex]) {
                        modelComponents.exitRates[stateIndex] = modelComponents.transitionMatrix.getRowSum(indices[stateIndex]);
                    } else {
                        modelComponents.exitRates[stateIndex] = storm::utility::zero<ValueType>();
                    }
                }
                STORM_LOG_TRACE("Exit rates: " << modelComponents.exitRates);

                std::shared_ptr<storm::models::sparse::MarkovAutomaton<ValueType>> ma;
                if (copy) {
                    ma = std::make_shared<storm::models::sparse::MarkovAutomaton<ValueType>>(modelComponents.transitionMatrix, modelComponents.stateLabeling, modelComponents.markovianStates, modelComponents.exitRates);
                } else {
                    ma = std::make_shared<storm::models::sparse::MarkovAutomaton<ValueType>>(std::move(modelComponents.transitionMatrix), std::move(modelComponents.stateLabeling), std::move(modelComponents.markovianStates), std::move(modelComponents.exitRates));
                }
                if (ma->hasOnlyTrivialNondeterminism()) {
                    // Markov automaton can be converted into CTMC
                    // TODO Matthias: change components which were not moved accordingly
                    model = ma->convertToCTMC();
                } else {
                    model = ma;
                }
            }

            STORM_LOG_DEBUG("No. states: " << model->getNumberOfStates());
            STORM_LOG_DEBUG("No. transitions: " << model->getNumberOfTransitions());
            if (model->getNumberOfStates() <= 15) {
                STORM_LOG_TRACE("Transition matrix: " << std::endl << model->getTransitionMatrix());
            } else {
                STORM_LOG_TRACE("Transition matrix: too big to print");
            }
            return model;
        }

        template<typename ValueType, typename StateType>
        void ExplicitDFTModelBuilderApprox<ValueType, StateType>::changeMatrixLowerBound(storm::storage::SparseMatrix<ValueType> & matrix) const {
            // Set lower bound for skipped states
            for (auto it = skippedStates.begin(); it != skippedStates.end(); ++it) {
                auto matrixEntry = matrix.getRow(it->first, 0).begin();
                STORM_LOG_ASSERT(matrixEntry->getColumn() == failedStateId, "Transition has wrong target state.");
                // Get the lower bound by considering the failure of all possible BEs
                ValueType newRate = storm::utility::zero<ValueType>();
                for (size_t index = 0; index < it->second->nrFailableBEs(); ++index) {
                    newRate += it->second->getFailableBERate(index);
                }
                for (size_t index = 0; index < it->second->nrNotFailableBEs(); ++index) {
                    newRate += it->second->getNotFailableBERate(index);
                }
                matrixEntry->setValue(newRate);
            }
        }

        template<typename ValueType, typename StateType>
        void ExplicitDFTModelBuilderApprox<ValueType, StateType>::changeMatrixUpperBound(storm::storage::SparseMatrix<ValueType> & matrix) const {
            // Set uppper bound for skipped states
            for (auto it = skippedStates.begin(); it != skippedStates.end(); ++it) {
                auto matrixEntry = matrix.getRow(it->first, 0).begin();
                STORM_LOG_ASSERT(matrixEntry->getColumn() == failedStateId, "Transition has wrong target state.");
                // Get the upper bound by considering the failure of all BEs
                // The used formula for the rate is 1/( 1/a + 1/b + ...)
                // TODO Matthias: improve by using closed formula for AND of all BEs
                ValueType newRate = storm::utility::zero<ValueType>();
                for (size_t index = 0; index < it->second->nrFailableBEs(); ++index) {
                    newRate += storm::utility::one<ValueType>() / it->second->getFailableBERate(index);
                }
                for (size_t index = 0; index < it->second->nrNotFailableBEs(); ++index) {
                    newRate += storm::utility::one<ValueType>() / it->second->getNotFailableBERate(index);
                }
                newRate = storm::utility::one<ValueType>() / newRate;
                matrixEntry->setValue(newRate);
            }
        }

        template<typename ValueType, typename StateType>
        StateType ExplicitDFTModelBuilderApprox<ValueType, StateType>::getOrAddStateIndex(DFTStatePointer const& state) {
            StateType stateId;
            bool changed = false;

            if (stateGenerationInfo->hasSymmetries()) {
                // Order state by symmetry
                STORM_LOG_TRACE("Check for symmetry: " << dft.getStateString(state));
                changed = state->orderBySymmetry();
                STORM_LOG_TRACE("State " << (changed ? "changed to " : "did not change") << (changed ? dft.getStateString(state) : ""));
            }

            if (stateStorage.stateToId.contains(state->status())) {
                // State already exists
                stateId = stateStorage.stateToId.getValue(state->status());
                STORM_LOG_TRACE("State " << dft.getStateString(state) << " with id " << stateId << " already exists");

                if (!changed && stateId >= OFFSET_PSEUDO_STATE) {
                    // Pseudo state can be created now
                    STORM_LOG_ASSERT(stateId >= OFFSET_PSEUDO_STATE, "State is no pseudo state.");
                    stateId -= OFFSET_PSEUDO_STATE;
                    STORM_LOG_ASSERT(stateId < mPseudoStatesMapping.size(), "Pseudo state not known.");
                    STORM_LOG_ASSERT(mPseudoStatesMapping[stateId].first == 0, "Pseudo state already created.");
                    // Create pseudo state now
                    STORM_LOG_ASSERT(mPseudoStatesMapping[stateId].second == state->status(), "Pseudo states do not coincide.");
                    state->setId(newIndex++);
                    mPseudoStatesMapping[stateId].first = state->getId();
                    stateId = state->getId();
                    stateStorage.stateToId.setOrAdd(state->status(), stateId);
                    STORM_LOG_TRACE("Now create state " << dft.getStateString(state) << " with id " << stateId);
                    statesToExplore.push(state);
                }
            } else {
                // State does not exist yet
                if (changed) {
                    // Remember state for later creation
                    state->setId(mPseudoStatesMapping.size() + OFFSET_PSEUDO_STATE);
                    mPseudoStatesMapping.push_back(std::make_pair(0, state->status()));
                    stateId = stateStorage.stateToId.findOrAdd(state->status(), state->getId());
                    STORM_LOG_ASSERT(stateId == state->getId(), "Ids do not match.");
                    STORM_LOG_TRACE("Remember state for later creation: " << dft.getStateString(state));
                    // Reserve one slot for the coming state in the remapping
                    matrixBuilder.stateRemapping.push_back(0);
                } else {
                    // Create new state
                    state->setId(newIndex++);
                    stateId = stateStorage.stateToId.findOrAdd(state->status(), state->getId());
                    STORM_LOG_ASSERT(stateId == state->getId(), "Ids do not match.");
                    STORM_LOG_TRACE("New state: " << dft.getStateString(state));
                    statesToExplore.push(state);

                    // Reserve one slot for the new state in the remapping
                    matrixBuilder.stateRemapping.push_back(0);
                }
            }
            return stateId;
        }

        template<typename ValueType, typename StateType>
        void ExplicitDFTModelBuilderApprox<ValueType, StateType>::setMarkovian(bool markovian) {
            if (matrixBuilder.getCurrentRowGroup() > modelComponents.markovianStates.size()) {
                // Resize BitVector
                modelComponents.markovianStates.resize(modelComponents.markovianStates.size() + INITIAL_BITVECTOR_SIZE);
            }
            modelComponents.markovianStates.set(matrixBuilder.getCurrentRowGroup() - 1, markovian);
        }


        // Explicitly instantiate the class.
        template class ExplicitDFTModelBuilderApprox<double>;

#ifdef STORM_HAVE_CARL
        template class ExplicitDFTModelBuilderApprox<storm::RationalFunction>;
#endif

    } // namespace builder
} // namespace storm


