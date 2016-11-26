#include "storm/abstraction/MenuGameRefiner.h"

#include "storm/abstraction/AbstractionInformation.h"
#include "storm/abstraction/MenuGameAbstractor.h"

#include "storm/utility/dd.h"

#include "storm/settings/SettingsManager.h"
#include "storm/settings/modules/AbstractionSettings.h"

namespace storm {
    namespace abstraction {
        
        template<storm::dd::DdType Type, typename ValueType>
        MenuGameRefiner<Type, ValueType>::MenuGameRefiner(MenuGameAbstractor<Type, ValueType>& abstractor, std::unique_ptr<storm::solver::SmtSolver>&& smtSolver) : abstractor(abstractor), splitPredicates(storm::settings::getModule<storm::settings::modules::AbstractionSettings>().isSplitPredicatesSet()), splitter(), equivalenceChecker(std::move(smtSolver)) {
            // Intentionally left empty.
        }
        
        template<storm::dd::DdType Type, typename ValueType>
        void MenuGameRefiner<Type, ValueType>::refine(std::vector<storm::expressions::Expression> const& predicates) const {
            abstractor.get().refine(predicates);
        }
        
        template<storm::dd::DdType Type>
        storm::dd::Bdd<Type> pickPivotState(storm::dd::Bdd<Type> const& initialStates, storm::dd::Bdd<Type> const& transitions, std::set<storm::expressions::Variable> const& rowVariables, std::set<storm::expressions::Variable> const& columnVariables, storm::dd::Bdd<Type> const& pivotStates) {
            
            // Perform a BFS and pick the first pivot state we encounter.
            storm::dd::Bdd<Type> pivotState;
            storm::dd::Bdd<Type> frontier = initialStates;
            storm::dd::Bdd<Type> frontierPivotStates = frontier && pivotStates;
            
            uint64_t level = 0;
            bool foundPivotState = !frontierPivotStates.isZero();
            if (foundPivotState) {
                pivotState = frontierPivotStates.existsAbstractRepresentative(rowVariables);
                STORM_LOG_TRACE("Picked pivot state from " << frontierPivotStates.getNonZeroCount() << " candidates on level " << level << ", " << pivotStates.getNonZeroCount() << " candidates in total.");
            } else {
                while (!foundPivotState) {
                    frontier = frontier.relationalProduct(transitions, rowVariables, columnVariables);
                    frontierPivotStates = frontier && pivotStates;
                    
                    if (!frontierPivotStates.isZero()) {
                        STORM_LOG_TRACE("Picked pivot state from " << frontierPivotStates.getNonZeroCount() << " candidates on level " << level << ", " << pivotStates.getNonZeroCount() << " candidates in total.");
                        pivotState = frontierPivotStates.existsAbstractRepresentative(rowVariables);
                        foundPivotState = true;
                    }
                    ++level;
                }
            }
            
            return pivotState;
        }
        
        template <storm::dd::DdType Type, typename ValueType>
        void MenuGameRefiner<Type, ValueType>::refine(storm::dd::Bdd<Type> const& pivotState, storm::dd::Bdd<Type> const& player1Choice, storm::dd::Bdd<Type> const& lowerChoice, storm::dd::Bdd<Type> const& upperChoice) const {
            AbstractionInformation<Type> const& abstractionInformation = abstractor.get().getAbstractionInformation();
            
            // Decode the index of the command chosen by player 1.
            storm::dd::Add<Type, ValueType> player1ChoiceAsAdd = player1Choice.template toAdd<ValueType>();
            auto pl1It = player1ChoiceAsAdd.begin();
            uint_fast64_t player1Index = abstractionInformation.decodePlayer1Choice((*pl1It).first, abstractionInformation.getPlayer1VariableCount());
#ifdef LOCAL_DEBUG
            std::cout << "command index " << commandIndex << std::endl;
            std::cout << program.get() << std::endl;
            
            for (auto stateValue : pivotState.template toAdd<ValueType>()) {
                std::stringstream stateName;
                stateName << "\tpl1_";
                for (auto const& var : currentGame->getRowVariables()) {
                    std::cout << "var " << var.getName() << std::endl;
                    if (stateValue.first.getBooleanValue(var)) {
                        stateName << "1";
                    } else {
                        stateName << "0";
                    }
                }
                std::cout << "pivot is " << stateName.str() << std::endl;
            }
#endif
            // storm::abstraction::prism::AbstractCommand<Type, ValueType>& abstractCommand = modules.front().getCommands()[commandIndex];
            // storm::prism::Command const& concreteCommand = abstractCommand.getConcreteCommand();
#ifdef LOCAL_DEBUG
            player1Choice.template toAdd<ValueType>().exportToDot("pl1choice_ref.dot");
            std::cout << concreteCommand << std::endl;
            
            (currentGame->getTransitionMatrix() * player1Choice.template toAdd<ValueType>()).exportToDot("cuttrans.dot");
#endif
            
            // Check whether there are bottom states in the game and whether one of the choices actually picks the
            // bottom state as the successor.
            bool buttomStateSuccessor = !((abstractionInformation.getBottomStateBdd(false, false) && lowerChoice) || (abstractionInformation.getBottomStateBdd(false, false) && upperChoice)).isZero();
            
            // If one of the choices picks the bottom state, the new predicate is based on the guard of the appropriate
            // command (that is the player 1 choice).
            if (buttomStateSuccessor) {
                STORM_LOG_TRACE("One of the successors is a bottom state, taking a guard as a new predicate.");
                storm::expressions::Expression newPredicate = abstractor.get().getGuard(player1Index);
                STORM_LOG_DEBUG("Derived new predicate: " << newPredicate);
                this->performRefinement({newPredicate});
            } else {
                STORM_LOG_TRACE("No bottom state successor. Deriving a new predicate using weakest precondition.");
                
#ifdef LOCAL_DEBUG
                lowerChoice.template toAdd<ValueType>().exportToDot("lowerchoice_ref.dot");
                upperChoice.template toAdd<ValueType>().exportToDot("upperchoice_ref.dot");
#endif
                
                // Decode both choices to explicit mappings.
#ifdef LOCAL_DEBUG
                std::cout << "lower" << std::endl;
#endif
                std::map<uint_fast64_t, storm::storage::BitVector> lowerChoiceUpdateToSuccessorMapping = abstractionInformation.decodeChoiceToUpdateSuccessorMapping(lowerChoice);
#ifdef LOCAL_DEBUG
                std::cout << "upper" << std::endl;
#endif
                std::map<uint_fast64_t, storm::storage::BitVector> upperChoiceUpdateToSuccessorMapping = abstractionInformation.decodeChoiceToUpdateSuccessorMapping(upperChoice);
                STORM_LOG_ASSERT(lowerChoiceUpdateToSuccessorMapping.size() == upperChoiceUpdateToSuccessorMapping.size(), "Mismatching sizes after decode (" << lowerChoiceUpdateToSuccessorMapping.size() << " vs. " << upperChoiceUpdateToSuccessorMapping.size() << ").");
                
#ifdef LOCAL_DEBUG
                std::cout << "lower" << std::endl;
                for (auto const& entry : lowerChoiceUpdateToSuccessorMapping) {
                    std::cout << entry.first << " -> " << entry.second << std::endl;
                }
                std::cout << "upper" << std::endl;
                for (auto const& entry : upperChoiceUpdateToSuccessorMapping) {
                    std::cout << entry.first << " -> " << entry.second << std::endl;
                }
#endif
                
                // Now go through the mappings and find points of deviation. Currently, we take the first deviation.
                storm::expressions::Expression newPredicate;
                auto lowerIt = lowerChoiceUpdateToSuccessorMapping.begin();
                auto lowerIte = lowerChoiceUpdateToSuccessorMapping.end();
                auto upperIt = upperChoiceUpdateToSuccessorMapping.begin();
                for (; lowerIt != lowerIte; ++lowerIt, ++upperIt) {
                    STORM_LOG_ASSERT(lowerIt->first == upperIt->first, "Update indices mismatch.");
                    uint_fast64_t updateIndex = lowerIt->first;
#ifdef LOCAL_DEBUG
                    std::cout << "update idx " << updateIndex << std::endl;
#endif
                    bool deviates = lowerIt->second != upperIt->second;
                    if (deviates) {
                        for (uint_fast64_t predicateIndex = 0; predicateIndex < lowerIt->second.size(); ++predicateIndex) {
                            if (lowerIt->second.get(predicateIndex) != upperIt->second.get(predicateIndex)) {
                                // Now we know the point of the deviation (command, update, predicate).
                                std::cout << "ref" << std::endl;
                                std::cout << abstractionInformation.getPredicateByIndex(predicateIndex) << std::endl;
                                newPredicate = abstractionInformation.getPredicateByIndex(predicateIndex).substitute(abstractor.get().getVariableUpdates(player1Index, updateIndex)).simplify();
                                break;
                            }
                        }
                    }
                }
                STORM_LOG_ASSERT(newPredicate.isInitialized(), "Could not derive new predicate as there is no deviation.");
                
                STORM_LOG_DEBUG("Derived new predicate: " << newPredicate);
                this->performRefinement({newPredicate});
            }
            
            STORM_LOG_TRACE("Current set of predicates:");
            for (auto const& predicate : abstractionInformation.getPredicates()) {
                STORM_LOG_TRACE(predicate);
            }
        }
        
        template<storm::dd::DdType Type, typename ValueType>
        bool MenuGameRefiner<Type, ValueType>::refine(storm::abstraction::MenuGame<Type, ValueType> const& game, storm::dd::Bdd<Type> const& transitionMatrixBdd, QualitativeResultMinMax<Type> const& qualitativeResult) const {
            STORM_LOG_TRACE("Trying refinement after qualitative check.");
            // Get all relevant strategies.
            storm::dd::Bdd<Type> minPlayer1Strategy = qualitativeResult.prob0Min.getPlayer1Strategy();
            storm::dd::Bdd<Type> minPlayer2Strategy = qualitativeResult.prob0Min.getPlayer2Strategy();
            storm::dd::Bdd<Type> maxPlayer1Strategy = qualitativeResult.prob1Max.getPlayer1Strategy();
            storm::dd::Bdd<Type> maxPlayer2Strategy = qualitativeResult.prob1Max.getPlayer2Strategy();
            
            // Redirect all player 1 choices of the min strategy to that of the max strategy if this leads to a player 2
            // state that is also a prob 0 state.
            minPlayer1Strategy = (maxPlayer1Strategy && qualitativeResult.prob0Min.getPlayer2States()).existsAbstract(game.getPlayer1Variables()).ite(maxPlayer1Strategy, minPlayer1Strategy);
            
            // Build the fragment of transitions that is reachable by both the min and the max strategies.
            storm::dd::Bdd<Type> reachableTransitions = transitionMatrixBdd && (minPlayer1Strategy || minPlayer2Strategy) && maxPlayer1Strategy && maxPlayer2Strategy;
            reachableTransitions = reachableTransitions.existsAbstract(game.getNondeterminismVariables());
            storm::dd::Bdd<Type> pivotStates = storm::utility::dd::computeReachableStates(game.getInitialStates(), reachableTransitions, game.getRowVariables(), game.getColumnVariables());
            
            // Require the pivot state to be a [0, 1] state.
            // TODO: is this restriction necessary or is it already implied?
            //            pivotStates &= prob01.min.first.getPlayer1States() && prob01.max.second.getPlayer1States();
            
            // Then constrain these states by the requirement that for either the lower or upper player 1 choice the player 2 choices must be different and
            // that the difference is not because of a missing strategy in either case.
            
            // Start with constructing the player 2 states that have a prob 0 (min) and prob 1 (max) strategy.
            storm::dd::Bdd<Type> constraint = minPlayer2Strategy.existsAbstract(game.getPlayer2Variables()) && maxPlayer2Strategy.existsAbstract(game.getPlayer2Variables());
            
            // Now construct all player 2 choices that actually exist and differ in the min and max case.
            constraint &= minPlayer2Strategy.exclusiveOr(maxPlayer2Strategy);
            
            // Then restrict the pivot states by requiring existing and different player 2 choices.
            pivotStates &= ((minPlayer1Strategy || maxPlayer1Strategy) && constraint).existsAbstract(game.getNondeterminismVariables());
            
            // We can only refine in case we have a reachable player 1 state with a player 2 successor (under either
            // player 1's min or max strategy) such that from this player 2 state, both prob0 min and prob0 max define
            // strategies and they differ. Hence, it is possible that we arrive at a point where no suitable pivot state
            // is found. In this case, we abort the qualitative refinement here.
            if (pivotStates.isZero()) {
                return false;
            }
            
            STORM_LOG_ASSERT(!pivotStates.isZero(), "Unable to proceed without pivot state candidates.");
            
            // Now that we have the pivot state candidates, we need to pick one.
            storm::dd::Bdd<Type> pivotState = pickPivotState<Type>(game.getInitialStates(), reachableTransitions, game.getRowVariables(), game.getColumnVariables(), pivotStates);
            
            // Compute the lower and the upper choice for the pivot state.
            std::set<storm::expressions::Variable> variablesToAbstract = game.getNondeterminismVariables();
            variablesToAbstract.insert(game.getRowVariables().begin(), game.getRowVariables().end());
            storm::dd::Bdd<Type> lowerChoice = pivotState && game.getExtendedTransitionMatrix().toBdd() && minPlayer1Strategy;
            storm::dd::Bdd<Type> lowerChoice1 = (lowerChoice && minPlayer2Strategy).existsAbstract(variablesToAbstract);
            storm::dd::Bdd<Type> lowerChoice2 = (lowerChoice && maxPlayer2Strategy).existsAbstract(variablesToAbstract);
            
            bool lowerChoicesDifferent = !lowerChoice1.exclusiveOr(lowerChoice2).isZero();
            if (lowerChoicesDifferent) {
                STORM_LOG_TRACE("Refining based on lower choice.");
                auto refinementStart = std::chrono::high_resolution_clock::now();
                
                this->refine(pivotState, (pivotState && minPlayer1Strategy).existsAbstract(game.getRowVariables()), lowerChoice1, lowerChoice2);
                auto refinementEnd = std::chrono::high_resolution_clock::now();
                STORM_LOG_TRACE("Refinement completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(refinementEnd - refinementStart).count() << "ms.");
                return true;
            } else {
                storm::dd::Bdd<Type> upperChoice = pivotState && game.getExtendedTransitionMatrix().toBdd() && maxPlayer1Strategy;
                storm::dd::Bdd<Type> upperChoice1 = (upperChoice && minPlayer2Strategy).existsAbstract(variablesToAbstract);
                storm::dd::Bdd<Type> upperChoice2 = (upperChoice && maxPlayer2Strategy).existsAbstract(variablesToAbstract);
                
                bool upperChoicesDifferent = !upperChoice1.exclusiveOr(upperChoice2).isZero();
                if (upperChoicesDifferent) {
                    STORM_LOG_TRACE("Refining based on upper choice.");
                    auto refinementStart = std::chrono::high_resolution_clock::now();
                    this->refine(pivotState, (pivotState && maxPlayer1Strategy).existsAbstract(game.getRowVariables()), upperChoice1, upperChoice2);
                    auto refinementEnd = std::chrono::high_resolution_clock::now();
                    STORM_LOG_TRACE("Refinement completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(refinementEnd - refinementStart).count() << "ms.");
                    return true;
                } else {
                    STORM_LOG_ASSERT(false, "Did not find choices from which to derive predicates.");
                }
            }
            return false;
        }
        
        template<storm::dd::DdType Type, typename ValueType>
        bool MenuGameRefiner<Type, ValueType>::refine(storm::abstraction::MenuGame<Type, ValueType> const& game, storm::dd::Bdd<Type> const& transitionMatrixBdd, QuantitativeResultMinMax<Type, ValueType> const& quantitativeResult) const {
            STORM_LOG_TRACE("Refining after quantitative check.");
            // Get all relevant strategies.
            storm::dd::Bdd<Type> minPlayer1Strategy = quantitativeResult.min.player1Strategy;
            storm::dd::Bdd<Type> minPlayer2Strategy = quantitativeResult.min.player2Strategy;
            storm::dd::Bdd<Type> maxPlayer1Strategy = quantitativeResult.max.player1Strategy;
            storm::dd::Bdd<Type> maxPlayer2Strategy = quantitativeResult.max.player2Strategy;
            
            // TODO: fix min strategies to take the max strategies if possible.
            
            // Build the fragment of transitions that is reachable by both the min and the max strategies.
            storm::dd::Bdd<Type> reachableTransitions = transitionMatrixBdd && (minPlayer1Strategy || maxPlayer1Strategy) && minPlayer2Strategy && maxPlayer2Strategy;
            reachableTransitions = reachableTransitions.existsAbstract(game.getNondeterminismVariables());
            storm::dd::Bdd<Type> pivotStates = storm::utility::dd::computeReachableStates(game.getInitialStates(), reachableTransitions, game.getRowVariables(), game.getColumnVariables());
            
            // Require the pivot state to be a state with a lower bound strictly smaller than the upper bound.
            pivotStates &= quantitativeResult.min.values.less(quantitativeResult.max.values);
            
            STORM_LOG_ASSERT(!pivotStates.isZero(), "Unable to refine without pivot state candidates.");
            
            // Then constrain these states by the requirement that for either the lower or upper player 1 choice the player 2 choices must be different and
            // that the difference is not because of a missing strategy in either case.
            
            // Start with constructing the player 2 states that have a (min) and a (max) strategy.
            // TODO: necessary?
            storm::dd::Bdd<Type> constraint = minPlayer2Strategy.existsAbstract(game.getPlayer2Variables()) && maxPlayer2Strategy.existsAbstract(game.getPlayer2Variables());
            
            // Now construct all player 2 choices that actually exist and differ in the min and max case.
            constraint &= minPlayer2Strategy.exclusiveOr(maxPlayer2Strategy);
            
            // Then restrict the pivot states by requiring existing and different player 2 choices.
            pivotStates &= ((minPlayer1Strategy || maxPlayer1Strategy) && constraint).existsAbstract(game.getNondeterminismVariables());
            
            STORM_LOG_ASSERT(!pivotStates.isZero(), "Unable to refine without pivot state candidates.");
            
            // Now that we have the pivot state candidates, we need to pick one.
            storm::dd::Bdd<Type> pivotState = pickPivotState<Type>(game.getInitialStates(), reachableTransitions, game.getRowVariables(), game.getColumnVariables(), pivotStates);
            
            // Compute the lower and the upper choice for the pivot state.
            std::set<storm::expressions::Variable> variablesToAbstract = game.getNondeterminismVariables();
            variablesToAbstract.insert(game.getRowVariables().begin(), game.getRowVariables().end());
            storm::dd::Bdd<Type> lowerChoice = pivotState && game.getExtendedTransitionMatrix().toBdd() && minPlayer1Strategy;
            storm::dd::Bdd<Type> lowerChoice1 = (lowerChoice && minPlayer2Strategy).existsAbstract(variablesToAbstract);
            storm::dd::Bdd<Type> lowerChoice2 = (lowerChoice && maxPlayer2Strategy).existsAbstract(variablesToAbstract);
            
            bool lowerChoicesDifferent = !lowerChoice1.exclusiveOr(lowerChoice2).isZero();
            if (lowerChoicesDifferent) {
                STORM_LOG_TRACE("Refining based on lower choice.");
                auto refinementStart = std::chrono::high_resolution_clock::now();
                this->refine(pivotState, (pivotState && minPlayer1Strategy).existsAbstract(game.getRowVariables()), lowerChoice1, lowerChoice2);
                auto refinementEnd = std::chrono::high_resolution_clock::now();
                STORM_LOG_TRACE("Refinement completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(refinementEnd - refinementStart).count() << "ms.");
            } else {
                storm::dd::Bdd<Type> upperChoice = pivotState && game.getExtendedTransitionMatrix().toBdd() && maxPlayer1Strategy;
                storm::dd::Bdd<Type> upperChoice1 = (upperChoice && minPlayer2Strategy).existsAbstract(variablesToAbstract);
                storm::dd::Bdd<Type> upperChoice2 = (upperChoice && maxPlayer2Strategy).existsAbstract(variablesToAbstract);
                
                bool upperChoicesDifferent = !upperChoice1.exclusiveOr(upperChoice2).isZero();
                if (upperChoicesDifferent) {
                    STORM_LOG_TRACE("Refining based on upper choice.");
                    auto refinementStart = std::chrono::high_resolution_clock::now();
                    this->refine(pivotState, (pivotState && maxPlayer1Strategy).existsAbstract(game.getRowVariables()), upperChoice1, upperChoice2);
                    auto refinementEnd = std::chrono::high_resolution_clock::now();
                    STORM_LOG_TRACE("Refinement completed in " << std::chrono::duration_cast<std::chrono::milliseconds>(refinementEnd - refinementStart).count() << "ms.");
                } else {
                    STORM_LOG_ASSERT(false, "Did not find choices from which to derive predicates.");
                }
            }
        }

        
        template<storm::dd::DdType Type, typename ValueType>
        bool MenuGameRefiner<Type, ValueType>::performRefinement(std::vector<storm::expressions::Expression> const& predicates) const {
            if (splitPredicates) {
                std::vector<storm::expressions::Expression> cleanedAtoms;

                for (auto const& predicate : predicates) {
                    AbstractionInformation<Type> const& abstractionInformation = abstractor.get().getAbstractionInformation();

                    // Split the predicates.
                    std::vector<storm::expressions::Expression> atoms = splitter.split(predicate);
                    
                    // Check which of the atoms are redundant in the sense that they are equivalent to a predicate we already have.
                    for (auto const& atom : atoms) {
                        bool addAtom = true;
                        for (auto const& oldPredicate : abstractionInformation.getPredicates()) {
                            if (equivalenceChecker.areEquivalent(atom, oldPredicate)) {
                                addAtom = false;
                                break;
                            }
                        }
                        
                        if (addAtom) {
                            cleanedAtoms.push_back(atom);
                        }
                    }
                }
                
                abstractor.get().refine(cleanedAtoms);
            } else {
                abstractor.get().refine(predicates);
            }
        }
        
        template class MenuGameRefiner<storm::dd::DdType::CUDD, double>;
        template class MenuGameRefiner<storm::dd::DdType::Sylvan, double>;
        
    }
}