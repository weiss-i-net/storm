#include "src/storage/prism/menu_games/AbstractionExpressionInformation.h"

#include "src/storage/expressions/ExpressionManager.h"
#include "src/storage/expressions/Expression.h"

namespace storm {
    namespace prism {
        namespace menu_games {
            
            AbstractionExpressionInformation::AbstractionExpressionInformation(storm::expressions::ExpressionManager& expressionManager, std::vector<storm::expressions::Expression> const& predicates, std::set<storm::expressions::Variable> const& variables) : expressionManager(expressionManager), predicates(predicates), variables(variables) {
                // Intentionally left empty.
            }
            
        }
    }
}