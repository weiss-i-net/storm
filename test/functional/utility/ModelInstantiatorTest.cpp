#include "gtest/gtest.h"
#include "storm-config.h"

#ifdef STORM_HAVE_CARL

#include "src/adapters/CarlAdapter.h"
#include<carl/numbers/numbers.h>
#include<carl/core/VariablePool.h>


#include "src/settings/SettingsManager.h"
#include "src/settings/modules/GeneralSettings.h"

#include "utility/storm.h"
#include "utility/ModelInstantiator.h"
#include "src/models/sparse/Model.h"
#include "src/models/sparse/Dtmc.h"
#include "src/models/sparse/Mdp.h"

TEST(ModelInstantiatorTest, BrpProb) {
    carl::VariablePool::getInstance().clear();
    
    std::string programFile = STORM_CPP_TESTS_BASE_PATH "/functional/utility/brp16_2.pm";
    std::string formulaAsString = "P=? [F s=5 ]";
    std::string constantsAsString = ""; //e.g. pL=0.9,TOACK=0.5
    
    // Program and formula
    storm::prism::Program program = storm::parseProgram(programFile);
    program.checkValidity();
    std::vector<std::shared_ptr<storm::logic::Formula>> formulas = storm::parseFormulasForProgram(formulaAsString, program);
    ASSERT_TRUE(formulas.size()==1);
    // Parametric model
    typename storm::builder::ExplicitPrismModelBuilder<storm::RationalFunction>::Options options = storm::builder::ExplicitPrismModelBuilder<storm::RationalFunction>::Options(*formulas[0]);
    options.addConstantDefinitionsFromString(program, constantsAsString); 
    options.preserveFormula(*formulas[0]);
    std::shared_ptr<storm::models::sparse::Dtmc<storm::RationalFunction>> dtmc = storm::builder::ExplicitPrismModelBuilder<storm::RationalFunction>(program, options).translate()->as<storm::models::sparse::Dtmc<storm::RationalFunction>>();
    
    storm::utility::ModelInstantiator<storm::models::sparse::Dtmc<storm::RationalFunction>, storm::models::sparse::Dtmc<double>> modelInstantiator(*dtmc);
    EXPECT_FALSE(dtmc->hasRewardModel());
    
    {
        std::map<storm::Variable, storm::CarlRationalNumber> valuation;
        storm::Variable const& pL = carl::VariablePool::getInstance().findVariableWithName("pL");
        ASSERT_NE(pL, carl::Variable::NO_VARIABLE);
        storm::Variable const& pK = carl::VariablePool::getInstance().findVariableWithName("pK");
        ASSERT_NE(pK, carl::Variable::NO_VARIABLE);
        valuation.insert(std::make_pair(pL,carl::rationalize<storm::CarlRationalNumber>(0.8)));
        valuation.insert(std::make_pair(pK,carl::rationalize<storm::CarlRationalNumber>(0.9)));

        storm::models::sparse::Dtmc<double> const& instantiated(modelInstantiator.instantiate(valuation));

        ASSERT_EQ(dtmc->getTransitionMatrix().getRowGroupIndices(), instantiated.getTransitionMatrix().getRowGroupIndices());
        for(std::size_t rowGroup = 0; rowGroup < dtmc->getTransitionMatrix().getRowGroupCount(); ++rowGroup){
            for(std::size_t row = dtmc->getTransitionMatrix().getRowGroupIndices()[rowGroup]; row < dtmc->getTransitionMatrix().getRowGroupIndices()[rowGroup+1]; ++row){
                auto instantiatedEntry = instantiated.getTransitionMatrix().getRow(row).begin();
                for(auto const& paramEntry : dtmc->getTransitionMatrix().getRow(row)){
                    EXPECT_EQ(paramEntry.getColumn(), instantiatedEntry->getColumn());
                    double evaluatedValue = carl::toDouble(paramEntry.getValue().evaluate(valuation));
                    EXPECT_EQ(evaluatedValue, instantiatedEntry->getValue());
                    ++instantiatedEntry;
                }
                EXPECT_EQ(instantiated.getTransitionMatrix().getRow(row).end(),instantiatedEntry);
            }
        }
        EXPECT_EQ(dtmc->getStateLabeling(), instantiated.getStateLabeling());
        EXPECT_EQ(dtmc->getOptionalChoiceLabeling(), instantiated.getOptionalChoiceLabeling());

        storm::modelchecker::SparseDtmcPrctlModelChecker<storm::models::sparse::Dtmc<double>> modelchecker(instantiated);
        std::unique_ptr<storm::modelchecker::CheckResult> chkResult = modelchecker.check(*formulas[0]);
        storm::modelchecker::ExplicitQuantitativeCheckResult<double>& quantitativeChkResult = chkResult->asExplicitQuantitativeCheckResult<double>();
        EXPECT_NEAR(0.2989278941, quantitativeChkResult[*instantiated.getInitialStates().begin()], storm::settings::generalSettings().getPrecision());
    }
    
    {
        std::map<storm::Variable, storm::CarlRationalNumber> valuation;
        storm::Variable const& pL = carl::VariablePool::getInstance().findVariableWithName("pL");
        ASSERT_NE(pL, carl::Variable::NO_VARIABLE);
        storm::Variable const& pK = carl::VariablePool::getInstance().findVariableWithName("pK");
        ASSERT_NE(pK, carl::Variable::NO_VARIABLE);
        valuation.insert(std::make_pair(pL,carl::rationalize<storm::CarlRationalNumber>(1)));
        valuation.insert(std::make_pair(pK,carl::rationalize<storm::CarlRationalNumber>(1)));

        storm::models::sparse::Dtmc<double> const& instantiated(modelInstantiator.instantiate(valuation));

        ASSERT_EQ(dtmc->getTransitionMatrix().getRowGroupIndices(), instantiated.getTransitionMatrix().getRowGroupIndices());
        for(std::size_t rowGroup = 0; rowGroup < dtmc->getTransitionMatrix().getRowGroupCount(); ++rowGroup){
            for(std::size_t row = dtmc->getTransitionMatrix().getRowGroupIndices()[rowGroup]; row < dtmc->getTransitionMatrix().getRowGroupIndices()[rowGroup+1]; ++row){
                auto instantiatedEntry = instantiated.getTransitionMatrix().getRow(row).begin();
                for(auto const& paramEntry : dtmc->getTransitionMatrix().getRow(row)){
                    EXPECT_EQ(paramEntry.getColumn(), instantiatedEntry->getColumn());
                    double evaluatedValue = carl::toDouble(paramEntry.getValue().evaluate(valuation));
                    EXPECT_EQ(evaluatedValue, instantiatedEntry->getValue());
                    ++instantiatedEntry;
                }
                EXPECT_EQ(instantiated.getTransitionMatrix().getRow(row).end(),instantiatedEntry);
            }
        }
        EXPECT_EQ(dtmc->getStateLabeling(), instantiated.getStateLabeling());
        EXPECT_EQ(dtmc->getOptionalChoiceLabeling(), instantiated.getOptionalChoiceLabeling());

        storm::modelchecker::SparseDtmcPrctlModelChecker<storm::models::sparse::Dtmc<double>> modelchecker(instantiated);
        std::unique_ptr<storm::modelchecker::CheckResult> chkResult = modelchecker.check(*formulas[0]);
        storm::modelchecker::ExplicitQuantitativeCheckResult<double>& quantitativeChkResult = chkResult->asExplicitQuantitativeCheckResult<double>();
        EXPECT_EQ(0.0 , quantitativeChkResult[*instantiated.getInitialStates().begin()]);
    }
    
    {
        std::map<storm::Variable, storm::CarlRationalNumber> valuation;
        storm::Variable const& pL = carl::VariablePool::getInstance().findVariableWithName("pL");
        ASSERT_NE(pL, carl::Variable::NO_VARIABLE);
        storm::Variable const& pK = carl::VariablePool::getInstance().findVariableWithName("pK");
        ASSERT_NE(pK, carl::Variable::NO_VARIABLE);
        valuation.insert(std::make_pair(pL,carl::rationalize<storm::CarlRationalNumber>(1)));
        valuation.insert(std::make_pair(pK,carl::rationalize<storm::CarlRationalNumber>(0.9)));

        storm::models::sparse::Dtmc<double> const& instantiated(modelInstantiator.instantiate(valuation));

        ASSERT_EQ(dtmc->getTransitionMatrix().getRowGroupIndices(), instantiated.getTransitionMatrix().getRowGroupIndices());
        for(std::size_t rowGroup = 0; rowGroup < dtmc->getTransitionMatrix().getRowGroupCount(); ++rowGroup){
            for(std::size_t row = dtmc->getTransitionMatrix().getRowGroupIndices()[rowGroup]; row < dtmc->getTransitionMatrix().getRowGroupIndices()[rowGroup+1]; ++row){
                auto instantiatedEntry = instantiated.getTransitionMatrix().getRow(row).begin();
                for(auto const& paramEntry : dtmc->getTransitionMatrix().getRow(row)){
                    EXPECT_EQ(paramEntry.getColumn(), instantiatedEntry->getColumn());
                    double evaluatedValue = carl::toDouble(paramEntry.getValue().evaluate(valuation));
                    EXPECT_EQ(evaluatedValue, instantiatedEntry->getValue());
                    ++instantiatedEntry;
                }
                EXPECT_EQ(instantiated.getTransitionMatrix().getRow(row).end(),instantiatedEntry);
            }
        }
        EXPECT_EQ(dtmc->getStateLabeling(), instantiated.getStateLabeling());
        EXPECT_EQ(dtmc->getOptionalChoiceLabeling(), instantiated.getOptionalChoiceLabeling());

        storm::modelchecker::SparseDtmcPrctlModelChecker<storm::models::sparse::Dtmc<double>> modelchecker(instantiated);
        std::unique_ptr<storm::modelchecker::CheckResult> chkResult = modelchecker.check(*formulas[0]);
        storm::modelchecker::ExplicitQuantitativeCheckResult<double>& quantitativeChkResult = chkResult->asExplicitQuantitativeCheckResult<double>();
        EXPECT_NEAR(0.01588055832, quantitativeChkResult[*instantiated.getInitialStates().begin()], storm::settings::generalSettings().getPrecision());
    }
}

TEST(ModelInstantiatorTest, Brp_Rew) {
    carl::VariablePool::getInstance().clear();
    
    std::string programFile = STORM_CPP_TESTS_BASE_PATH "/functional/utility/brp16_2.pm";
    std::string formulaAsString = "R=? [F ((s=5) | (s=0&srep=3)) ]";
    std::string constantsAsString = ""; //e.g. pL=0.9,TOACK=0.5
    
    // Program and formula
    storm::prism::Program program = storm::parseProgram(programFile);
    program.checkValidity();
    std::vector<std::shared_ptr<storm::logic::Formula>> formulas = storm::parseFormulasForProgram(formulaAsString, program);
    ASSERT_TRUE(formulas.size()==1);
    // Parametric model
    typename storm::builder::ExplicitPrismModelBuilder<storm::RationalFunction>::Options options = storm::builder::ExplicitPrismModelBuilder<storm::RationalFunction>::Options(*formulas[0]);
    options.addConstantDefinitionsFromString(program, constantsAsString); 
    options.preserveFormula(*formulas[0]);
    std::shared_ptr<storm::models::sparse::Dtmc<storm::RationalFunction>> dtmc = storm::builder::ExplicitPrismModelBuilder<storm::RationalFunction>(program, options).translate()->as<storm::models::sparse::Dtmc<storm::RationalFunction>>();

    storm::utility::ModelInstantiator<storm::models::sparse::Dtmc<storm::RationalFunction>, storm::models::sparse::Dtmc<double>> modelInstantiator(*dtmc);
    
    {
        std::map<storm::Variable, storm::CarlRationalNumber> valuation;
        storm::Variable const& pL = carl::VariablePool::getInstance().findVariableWithName("pL");
        ASSERT_NE(pL, carl::Variable::NO_VARIABLE);
        storm::Variable const& pK = carl::VariablePool::getInstance().findVariableWithName("pK");
        ASSERT_NE(pK, carl::Variable::NO_VARIABLE);
        storm::Variable const& TOMsg = carl::VariablePool::getInstance().findVariableWithName("TOMsg");
        ASSERT_NE(pK, carl::Variable::NO_VARIABLE);
        storm::Variable const& TOAck = carl::VariablePool::getInstance().findVariableWithName("TOAck");
        ASSERT_NE(pK, carl::Variable::NO_VARIABLE);
        valuation.insert(std::make_pair(pL,carl::rationalize<storm::CarlRationalNumber>(0.9)));
        valuation.insert(std::make_pair(pK,carl::rationalize<storm::CarlRationalNumber>(0.3)));
        valuation.insert(std::make_pair(TOMsg,carl::rationalize<storm::CarlRationalNumber>(0.3)));
        valuation.insert(std::make_pair(TOAck,carl::rationalize<storm::CarlRationalNumber>(0.5)));

        storm::models::sparse::Dtmc<double> const& instantiated(modelInstantiator.instantiate(valuation));

        ASSERT_EQ(dtmc->getTransitionMatrix().getRowGroupIndices(), instantiated.getTransitionMatrix().getRowGroupIndices());
        for(std::size_t rowGroup = 0; rowGroup < dtmc->getTransitionMatrix().getRowGroupCount(); ++rowGroup){
            for(std::size_t row = dtmc->getTransitionMatrix().getRowGroupIndices()[rowGroup]; row < dtmc->getTransitionMatrix().getRowGroupIndices()[rowGroup+1]; ++row){
                auto instantiatedEntry = instantiated.getTransitionMatrix().getRow(row).begin();
                for(auto const& paramEntry : dtmc->getTransitionMatrix().getRow(row)){
                    EXPECT_EQ(paramEntry.getColumn(), instantiatedEntry->getColumn());
                    double evaluatedValue = carl::toDouble(paramEntry.getValue().evaluate(valuation));
                    EXPECT_EQ(evaluatedValue, instantiatedEntry->getValue());
                    ++instantiatedEntry;
                }
                EXPECT_EQ(instantiated.getTransitionMatrix().getRow(row).end(),instantiatedEntry);
            }
        }
        ASSERT_TRUE(instantiated.hasUniqueRewardModel());
        EXPECT_FALSE(instantiated.getUniqueRewardModel()->second.hasStateRewards());
        EXPECT_FALSE(instantiated.getUniqueRewardModel()->second.hasTransitionRewards());
        EXPECT_TRUE(instantiated.getUniqueRewardModel()->second.hasStateActionRewards());
        ASSERT_TRUE(dtmc->getUniqueRewardModel()->second.hasStateActionRewards());
        std::size_t stateActionEntries = dtmc->getUniqueRewardModel()->second.getStateActionRewardVector().size();
        ASSERT_EQ(stateActionEntries, instantiated.getUniqueRewardModel()->second.getStateActionRewardVector().size());
        for(std::size_t i =0; i<stateActionEntries; ++i){
            double evaluatedValue = carl::toDouble(dtmc->getUniqueRewardModel()->second.getStateActionRewardVector()[i].evaluate(valuation));
            EXPECT_EQ(evaluatedValue, instantiated.getUniqueRewardModel()->second.getStateActionRewardVector()[i]);
        }
        EXPECT_EQ(dtmc->getStateLabeling(), instantiated.getStateLabeling());
        EXPECT_EQ(dtmc->getOptionalChoiceLabeling(), instantiated.getOptionalChoiceLabeling());

        storm::modelchecker::SparseDtmcPrctlModelChecker<storm::models::sparse::Dtmc<double>> modelchecker(instantiated);
        std::unique_ptr<storm::modelchecker::CheckResult> chkResult = modelchecker.check(*formulas[0]);
        storm::modelchecker::ExplicitQuantitativeCheckResult<double>& quantitativeChkResult = chkResult->asExplicitQuantitativeCheckResult<double>();
        EXPECT_NEAR(1.308324495, quantitativeChkResult[*instantiated.getInitialStates().begin()], storm::settings::generalSettings().getPrecision());
    }
    
}
    

TEST(ModelInstantiatorTest, Consensus) {
    carl::VariablePool::getInstance().clear();
    
    std::string programFile = STORM_CPP_TESTS_BASE_PATH "/functional/utility/coin2_2.pm";
    std::string formulaAsString = "Pmin=? [F \"finished\"&\"all_coins_equal_1\" ]";
    std::string constantsAsString = ""; //e.g. pL=0.9,TOACK=0.5
    
    // Program and formula
    storm::prism::Program program = storm::parseProgram(programFile);
    program.checkValidity();
    std::vector<std::shared_ptr<storm::logic::Formula>> formulas = storm::parseFormulasForProgram(formulaAsString, program);
    ASSERT_TRUE(formulas.size()==1);
    // Parametric model
    typename storm::builder::ExplicitPrismModelBuilder<storm::RationalFunction>::Options options = storm::builder::ExplicitPrismModelBuilder<storm::RationalFunction>::Options(*formulas[0]);
    options.addConstantDefinitionsFromString(program, constantsAsString); 
    options.preserveFormula(*formulas[0]);
    std::shared_ptr<storm::models::sparse::Mdp<storm::RationalFunction>> mdp = storm::builder::ExplicitPrismModelBuilder<storm::RationalFunction>(program, options).translate()->as<storm::models::sparse::Mdp<storm::RationalFunction>>();

    storm::utility::ModelInstantiator<storm::models::sparse::Mdp<storm::RationalFunction>, storm::models::sparse::Mdp<double>> modelInstantiator(*mdp);
    
    std::map<storm::Variable, storm::CarlRationalNumber> valuation;
    storm::Variable const& p1 = carl::VariablePool::getInstance().findVariableWithName("p1");
    ASSERT_NE(p1, carl::Variable::NO_VARIABLE);
    storm::Variable const& p2 = carl::VariablePool::getInstance().findVariableWithName("p2");
    ASSERT_NE(p2, carl::Variable::NO_VARIABLE);
    valuation.insert(std::make_pair(p1,carl::rationalize<storm::CarlRationalNumber>(0.51)));
    valuation.insert(std::make_pair(p2,carl::rationalize<storm::CarlRationalNumber>(0.49)));
    storm::models::sparse::Mdp<double> const& instantiated(modelInstantiator.instantiate(valuation));

    ASSERT_EQ(mdp->getTransitionMatrix().getRowGroupIndices(), instantiated.getTransitionMatrix().getRowGroupIndices());
    for(std::size_t rowGroup = 0; rowGroup < mdp->getTransitionMatrix().getRowGroupCount(); ++rowGroup){
        for(std::size_t row = mdp->getTransitionMatrix().getRowGroupIndices()[rowGroup]; row < mdp->getTransitionMatrix().getRowGroupIndices()[rowGroup+1]; ++row){
            auto instantiatedEntry = instantiated.getTransitionMatrix().getRow(row).begin();
            for(auto const& paramEntry : mdp->getTransitionMatrix().getRow(row)){
                EXPECT_EQ(paramEntry.getColumn(), instantiatedEntry->getColumn());
                double evaluatedValue = carl::toDouble(paramEntry.getValue().evaluate(valuation));
                EXPECT_EQ(evaluatedValue, instantiatedEntry->getValue());
                ++instantiatedEntry;
            }
            EXPECT_EQ(instantiated.getTransitionMatrix().getRow(row).end(),instantiatedEntry);
        }
    }
    EXPECT_EQ(mdp->getStateLabeling(), instantiated.getStateLabeling());
    EXPECT_EQ(mdp->getOptionalChoiceLabeling(), instantiated.getOptionalChoiceLabeling());

    storm::modelchecker::SparseMdpPrctlModelChecker<storm::models::sparse::Mdp<double>> modelchecker(instantiated);
    std::unique_ptr<storm::modelchecker::CheckResult> chkResult = modelchecker.check(*formulas[0]);
    storm::modelchecker::ExplicitQuantitativeCheckResult<double>& quantitativeChkResult = chkResult->asExplicitQuantitativeCheckResult<double>();
    EXPECT_NEAR(0.3526577219, quantitativeChkResult[*instantiated.getInitialStates().begin()], storm::settings::generalSettings().getPrecision());
}

#endif