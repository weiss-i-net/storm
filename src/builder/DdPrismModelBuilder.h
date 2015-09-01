#ifndef STORM_BUILDER_DDPRISMMODELBUILDER_H_
#define STORM_BUILDER_DDPRISMMODELBUILDER_H_

#include <map>
#include <boost/optional.hpp>

#include "src/logic/Formulas.h"
#include "src/adapters/AddExpressionAdapter.h"
#include "src/utility/macros.h"

namespace storm {
    namespace prism {
        class Program;
        class Module;
        class RewardModel;
        class Update;
        class Command;
    }
    
    namespace dd {
        template<storm::dd::DdType T> class Bdd;
    }
    
    namespace models {
        namespace symbolic {
            template<storm::dd::DdType T> class Model;
            template<storm::dd::DdType T, typename ValueType> class StandardRewardModel;
        }
    }
    
    namespace builder {
        
        template <storm::dd::DdType Type>
        class DdPrismModelBuilder {
        public:
            struct Options {
                /*!
                 * Creates an object representing the default building options.
                 */
                Options();
                
                /*! Creates an object representing the suggested building options assuming that the given formula is the
                 * only one to check. Additional formulas may be preserved by calling <code>preserveFormula</code>.
                 *
                 * @param formula The formula based on which to choose the building options.
                 */
                Options(storm::logic::Formula const& formula);
                
                /*! Creates an object representing the suggested building options assuming that the given formulas are
                 * the only ones to check. Additional formulas may be preserved by calling <code>preserveFormula</code>.
                 *
                 * @param formula Thes formula based on which to choose the building options.
                 */
                Options(std::vector<std::shared_ptr<storm::logic::Formula>> const& formulas);
                
                /*!
                 * Sets the constants definitions from the given string. The string must be of the form 'X=a,Y=b,Z=c',
                 * etc. where X,Y,Z are the variable names and a,b,c are the values of the constants.
                 *
                 * @param program The program managing the constants that shall be defined. Note that the program itself
                 * is not modified whatsoever.
                 * @param constantDefinitionString The string from which to parse the constants' values.
                 */
                void addConstantDefinitionsFromString(storm::prism::Program const& program, std::string const& constantDefinitionString);
                
                /*!
                 * Changes the options in a way that ensures that the given formula can be checked on the model once it
                 * has been built.
                 *
                 * @param formula The formula that is to be ''preserved''.
                 */
                void preserveFormula(storm::logic::Formula const& formula);
                
                /*!
                 * Analyzes the given formula and sets an expression for the states states of the model that can be
                 * treated as terminal states. Note that this may interfere with checking properties different than the
                 * one provided.
                 *
                 * @param formula The formula used to (possibly) derive an expression for the terminal states of the
                 * model.
                 */
                void setTerminalStatesFromFormula(storm::logic::Formula const& formula);
                
                // A flag that indicates whether or not all reward models are to be build.
                bool buildAllRewardModels;
                
                // A list of reward models to be build in case not all reward models are to be build.
                std::set<std::string> rewardModelsToBuild;
                
                // An optional mapping that, if given, contains defining expressions for undefined constants.
                boost::optional<std::map<storm::expressions::Variable, storm::expressions::Expression>> constantDefinitions;
                
                // A flag indicating whether all labels are to be build.
                bool buildAllLabels;
                
                // An optional set of labels that, if given, restricts the labels that are built.
                boost::optional<std::set<std::string>> labelsToBuild;
                
                // An optional set of expressions for which labels need to be built.
                boost::optional<std::vector<storm::expressions::Expression>> expressionLabels;
                
                // An optional expression or label that characterizes the terminal states of the model. If this is set,
                // the outgoing transitions of these states are replaced with a self-loop.
                boost::optional<boost::variant<storm::expressions::Expression, std::string>> terminalStates;
            };
            
            /*!
             * Translates the given program into a symbolic model (i.e. one that stores the transition relation as a
             * decision diagram).
             *
             * @param program The program to translate.
             * @return A pointer to the resulting model.
             */
            static std::shared_ptr<storm::models::symbolic::Model<Type>> translateProgram(storm::prism::Program const& program, Options const& options = Options());
            
        private:
            // This structure can store the decision diagrams representing a particular action.
            struct UpdateDecisionDiagram {
                UpdateDecisionDiagram() : updateDd(), assignedGlobalVariables() {
                    // Intentionally left empty.
                }

                UpdateDecisionDiagram(storm::dd::Add<Type> const& updateDd, std::set<storm::expressions::Variable> const& assignedGlobalVariables) : updateDd(updateDd), assignedGlobalVariables(assignedGlobalVariables) {
                    // Intentionally left empty.
                }
                
                // The DD representing the update behaviour.
                storm::dd::Add<Type> updateDd;
                
                // Keep track of the global variables that were written by this update.
                std::set<storm::expressions::Variable> assignedGlobalVariables;
            };
            
            // This structure can store the decision diagrams representing a particular action.
            struct ActionDecisionDiagram {
                ActionDecisionDiagram() : guardDd(), transitionsDd(), numberOfUsedNondeterminismVariables(0) {
                    // Intentionally left empty.
                }
                
                ActionDecisionDiagram(storm::dd::DdManager<Type> const& manager, std::set<storm::expressions::Variable> const& assignedGlobalVariables = std::set<storm::expressions::Variable>(), uint_fast64_t numberOfUsedNondeterminismVariables = 0) : guardDd(manager.getAddZero()), transitionsDd(manager.getAddZero()), numberOfUsedNondeterminismVariables(numberOfUsedNondeterminismVariables), assignedGlobalVariables(assignedGlobalVariables) {
                    // Intentionally left empty.
                }
                
                ActionDecisionDiagram(storm::dd::Add<Type> guardDd, storm::dd::Add<Type> transitionsDd, std::set<storm::expressions::Variable> const& assignedGlobalVariables = std::set<storm::expressions::Variable>(), uint_fast64_t numberOfUsedNondeterminismVariables = 0) : guardDd(guardDd), transitionsDd(transitionsDd), numberOfUsedNondeterminismVariables(numberOfUsedNondeterminismVariables), assignedGlobalVariables(assignedGlobalVariables) {
                    // Intentionally left empty.
                }
                
                ActionDecisionDiagram(ActionDecisionDiagram const& other) = default;
                ActionDecisionDiagram& operator=(ActionDecisionDiagram const& other) = default;
                
                // The guard of the action.
                storm::dd::Add<Type> guardDd;
                
                // The actual transitions (source and target states).
                storm::dd::Add<Type> transitionsDd;
                
                // The number of variables that are used to encode the nondeterminism.
                uint_fast64_t numberOfUsedNondeterminismVariables;
                
                // Keep track of the global variables that were written by this action.
                std::set<storm::expressions::Variable> assignedGlobalVariables;
            };
            
            // This structure holds all decision diagrams related to a module.
            struct ModuleDecisionDiagram {
                ModuleDecisionDiagram() : independentAction(), synchronizingActionToDecisionDiagramMap(), identity(), numberOfUsedNondeterminismVariables(0) {
                    // Intentionally left empty.
                }
                
                ModuleDecisionDiagram(storm::dd::DdManager<Type> const& manager) : independentAction(manager), synchronizingActionToDecisionDiagramMap(), identity(manager.getAddZero()), numberOfUsedNondeterminismVariables(0) {
                    // Intentionally left empty.
                }

                ModuleDecisionDiagram(ActionDecisionDiagram const& independentAction, std::map<uint_fast64_t, ActionDecisionDiagram> const& synchronizingActionToDecisionDiagramMap, storm::dd::Add<Type> const& identity, uint_fast64_t numberOfUsedNondeterminismVariables = 0) : independentAction(independentAction), synchronizingActionToDecisionDiagramMap(synchronizingActionToDecisionDiagramMap), identity(identity), numberOfUsedNondeterminismVariables(numberOfUsedNondeterminismVariables) {
                    // Intentionally left empty.
                }
                
                ModuleDecisionDiagram(ModuleDecisionDiagram const& other) = default;
                ModuleDecisionDiagram& operator=(ModuleDecisionDiagram const& other) = default;
                
                bool hasSynchronizingAction(uint_fast64_t actionIndex) {
                    return synchronizingActionToDecisionDiagramMap.find(actionIndex) != synchronizingActionToDecisionDiagramMap.end();
                }
                
                // The decision diagram for the independent action.
                ActionDecisionDiagram independentAction;
                
                // A mapping from synchronizing action indices to the decision diagram.
                std::map<uint_fast64_t, ActionDecisionDiagram> synchronizingActionToDecisionDiagramMap;
                
                // A decision diagram that represents the identity of this module.
                storm::dd::Add<Type> identity;
                
                // The number of variables encoding the nondeterminism that were actually used.
                uint_fast64_t numberOfUsedNondeterminismVariables;
            };
            
            /*!
             * Structure to store all information required to generate the model from the program.
             */
            class GenerationInformation;
            
            /*!
             * Structure to store the result of the system creation phase.
             */
            struct SystemResult;
        private:
            static std::set<storm::expressions::Variable> equalizeAssignedGlobalVariables(GenerationInformation const& generationInfo, ActionDecisionDiagram& action1, ActionDecisionDiagram& action2);
            
            static std::set<storm::expressions::Variable> equalizeAssignedGlobalVariables(GenerationInformation const& generationInfo, std::vector<ActionDecisionDiagram>& actionDds);
            
            static storm::dd::Add<Type> encodeChoice(GenerationInformation& generationInfo, uint_fast64_t nondeterminismVariableOffset, uint_fast64_t numberOfBinaryVariables, int_fast64_t value);
            
            static UpdateDecisionDiagram createUpdateDecisionDiagram(GenerationInformation& generationInfo, storm::prism::Module const& module, storm::dd::Add<Type> const& guard, storm::prism::Update const& update);

            static ActionDecisionDiagram createCommandDecisionDiagram(GenerationInformation& generationInfo, storm::prism::Module const& module, storm::prism::Command const& command);

            static ActionDecisionDiagram createActionDecisionDiagram(GenerationInformation& generationInfo, storm::prism::Module const& module, uint_fast64_t synchronizationActionIndex, uint_fast64_t nondeterminismVariableOffset);

            static ActionDecisionDiagram combineCommandsToActionMarkovChain(GenerationInformation& generationInfo, std::vector<ActionDecisionDiagram>& commandDds);

            static ActionDecisionDiagram combineCommandsToActionMDP(GenerationInformation& generationInfo, std::vector<ActionDecisionDiagram>& commandDds, uint_fast64_t nondeterminismVariableOffset);

            static ActionDecisionDiagram combineSynchronizingActions(GenerationInformation const& generationInfo, ActionDecisionDiagram const& action1, ActionDecisionDiagram const& action2);

            static ActionDecisionDiagram combineUnsynchronizedActions(GenerationInformation const& generationInfo, ActionDecisionDiagram& action1, ActionDecisionDiagram& action2, storm::dd::Add<Type> const& identityDd1, storm::dd::Add<Type> const& identityDd2);

            static ModuleDecisionDiagram createModuleDecisionDiagram(GenerationInformation& generationInfo, storm::prism::Module const& module, std::map<uint_fast64_t, uint_fast64_t> const& synchronizingActionToOffsetMap);

            static storm::dd::Add<Type> getSynchronizationDecisionDiagram(GenerationInformation& generationInfo, uint_fast64_t actionIndex = 0);
            
            static storm::dd::Add<Type> createSystemFromModule(GenerationInformation& generationInfo, ModuleDecisionDiagram const& module);
            
            static storm::models::symbolic::StandardRewardModel<Type, double> createRewardModelDecisionDiagrams(GenerationInformation& generationInfo, storm::prism::RewardModel const& rewardModel, ModuleDecisionDiagram const& globalModule, storm::dd::Add<Type> const& transitionMatrix, storm::dd::Add<Type> const& reachableStatesAdd, storm::dd::Add<Type> const& stateActionDd);
            
            static SystemResult createSystemDecisionDiagram(GenerationInformation& generationInfo);
            
            static storm::dd::Bdd<Type> createInitialStatesDecisionDiagram(GenerationInformation& generationInfo);

            static storm::dd::Bdd<Type> computeReachableStates(GenerationInformation& generationInfo, storm::dd::Bdd<Type> const& initialStates, storm::dd::Bdd<Type> const& transitions);
            
        };
        
    } // namespace adapters
} // namespace storm

#endif /* STORM_BUILDER_DDPRISMMODELBUILDER_H_ */