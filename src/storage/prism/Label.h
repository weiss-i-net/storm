#ifndef STORM_STORAGE_PRISM_LABEL_H_
#define STORM_STORAGE_PRISM_LABEL_H_

#include <map>

#include "src/storage/prism/LocatedInformation.h"
#include "src/storage/expressions/Expression.h"

namespace storm {
    namespace prism {
        class Label : public LocatedInformation {
        public:
            /*!
             * Creates a label with the given name and state predicate expression.
             *
             * @param labelName The name of the label.
             * @param statePredicateExpression The predicate that needs to hold before taking a transition with the previously
             * specified name in order to obtain the reward.
             * @param filename The filename in which the transition reward is defined.
             * @param lineNumber The line number in which the transition reward is defined.
             */
            Label(std::string const& labelName, storm::expressions::Expression const& statePredicateExpression, std::string const& filename = "", uint_fast64_t lineNumber = 0);
            
            // Create default implementations of constructors/assignment.
            Label() = default;
            Label(Label const& other) = default;
            Label& operator=(Label const& other)= default;
            Label(Label&& other) = default;
            Label& operator=(Label&& other) = default;
            
            /*!
             * Retrieves the name that is associated with this label.
             *
             * @return The name that is associated with this label.
             */
            std::string const& getLabelName() const;
            
            /*!
             * Retrieves the state predicate expression that is associated with this label.
             *
             * @return The state predicate expression that is associated with this label.
             */
            storm::expressions::Expression const& getStatePredicateExpression() const;
            
            /*!
             * Substitutes all identifiers in the expression of the label according to the given map.
             *
             * @param substitution The substitution to perform.
             * @return The resulting label.
             */
            Label substitute(std::map<std::string, storm::expressions::Expression> const& substitution) const;
            
            friend std::ostream& operator<<(std::ostream& stream, Label const& label);
            
        private:
            // The name of the label.
            std::string labelName;
            
            // A predicate that needs to be satisfied by states for the label to be attached.
            storm::expressions::Expression statePredicateExpression;
        };
    } // namespace prism
} // namespace storm

#endif /* STORM_STORAGE_PRISM_LABEL_H_ */