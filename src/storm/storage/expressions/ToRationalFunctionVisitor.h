#ifndef STORM_STORAGE_EXPRESSIONS_TORATIONALFUNCTIONVISITOR_H_
#define STORM_STORAGE_EXPRESSIONS_TORATIONALFUNCTIONVISITOR_H_

#include <unordered_map>

#include "storm/adapters/CarlAdapter.h"

#include "storm/storage/expressions/Expression.h"
#include "storm/storage/expressions/Expressions.h"
#include "storm/storage/expressions/ExpressionVisitor.h"
#include "storm/storage/expressions/Variable.h"

namespace storm {
    namespace expressions {

#ifdef STORM_HAVE_CARL
        template<typename RationalFunctionType>
        class ToRationalFunctionVisitor : public ExpressionVisitor {
        public:
            ToRationalFunctionVisitor();
            
            RationalFunctionType toRationalFunction(Expression const& expression);
            
            virtual boost::any visit(IfThenElseExpression const& expression, boost::any const& data) override;
            virtual boost::any visit(BinaryBooleanFunctionExpression const& expression, boost::any const& data) override;
            virtual boost::any visit(BinaryNumericalFunctionExpression const& expression, boost::any const& data) override;
            virtual boost::any visit(BinaryRelationExpression const& expression, boost::any const& data) override;
            virtual boost::any visit(VariableExpression const& expression, boost::any const& data) override;
            virtual boost::any visit(UnaryBooleanFunctionExpression const& expression, boost::any const& data) override;
            virtual boost::any visit(UnaryNumericalFunctionExpression const& expression, boost::any const& data) override;
            virtual boost::any visit(BooleanLiteralExpression const& expression, boost::any const& data) override;
            virtual boost::any visit(IntegerLiteralExpression const& expression, boost::any const& data) override;
            virtual boost::any visit(RationalLiteralExpression const& expression, boost::any const& data) override;
            
            void setMapping(storm::expressions::Variable const& variable, RationalFunctionType const& value);
            
        private:
            template<typename TP = typename RationalFunctionType::PolyType, carl::EnableIf<carl::needs_cache<TP>> = carl::dummy>
            RationalFunctionType convertVariableToPolynomial(carl::Variable const& variable) {
                return RationalFunctionType(typename RationalFunctionType::PolyType(typename RationalFunctionType::PolyType::PolyType(variable), cache));
            }
            
            template<typename TP = typename RationalFunctionType::PolyType, carl::DisableIf<carl::needs_cache<TP>> = carl::dummy>
            RationalFunctionType convertVariableToPolynomial(carl::Variable const& variable) {
                return RationalFunctionType(variable);
            }
            
            // A mapping from our variables to carl's.
            std::unordered_map<storm::expressions::Variable, carl::Variable> variableToVariableMap;
            
            // The cache that is used in case the underlying type needs a cache.
            std::shared_ptr<carl::Cache<carl::PolynomialFactorizationPair<RawPolynomial>>> cache;
            
            // A mapping from variables to their values.
            std::unordered_map<storm::expressions::Variable, RationalFunctionType> valueMapping;
        };
#endif
    }
}

#endif /* STORM_STORAGE_EXPRESSIONS_TORATIONALFUNCTIONVISITOR_H_ */
