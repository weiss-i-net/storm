/*
 * CslParser.cpp
 *
 *  Created on: 08.04.2013
 *      Author: Thomas Heinemann
 */

#include "src/parser/CslParser.h"
#include "src/utility/OsDetection.h"
#include "src/utility/constants.h"

// The action class headers.
#include "src/formula/Actions/AbstractAction.h"
#include "src/formula/Actions/BoundAction.h"
#include "src/formula/Actions/InvertAction.h"
#include "src/formula/Actions/FormulaAction.h"
#include "src/formula/Actions/RangeAction.h"
#include "src/formula/Actions/SortAction.h"

// If the parser fails due to ill-formed data, this exception is thrown.
#include "src/exceptions/WrongFormatException.h"

// Used for Boost spirit.
#include <boost/typeof/typeof.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix.hpp>

// Include headers for spirit iterators. Needed for diagnostics and input stream iteration.
#include <boost/spirit/include/classic_position_iterator.hpp>
#include <boost/spirit/include/support_multi_pass.hpp>

// Needed for file IO.
#include <fstream>
#include <iomanip>
#include <map>


// Some typedefs and namespace definitions to reduce code size.
typedef std::string::const_iterator BaseIteratorType;
typedef boost::spirit::classic::position_iterator2<BaseIteratorType> PositionIteratorType;
namespace qi = boost::spirit::qi;
namespace phoenix = boost::phoenix;

namespace storm {
namespace parser {

template<typename Iterator, typename Skipper>
struct CslParser::CslGrammar : qi::grammar<Iterator, storm::property::csl::CslFilter<double>*(), Skipper > {
	CslGrammar() : CslGrammar::base_type(start) {
		//This block contains helper rules that may be used several times
		freeIdentifierName = qi::lexeme[qi::alpha >> *(qi::alnum | qi::char_('_'))];
		comparisonType = (
				(qi::lit(">="))[qi::_val = storm::property::GREATER_EQUAL] |
				(qi::lit(">"))[qi::_val = storm::property::GREATER] |
				(qi::lit("<="))[qi::_val = storm::property::LESS_EQUAL] |
				(qi::lit("<"))[qi::_val = storm::property::LESS]);
		sortingCategory = (
				(qi::lit("index"))[qi::_val = storm::property::action::SortAction<double>::INDEX] |
				(qi::lit("value"))[qi::_val = storm::property::action::SortAction<double>::VALUE]
				);
		//Comment: Empty line or line starting with "//"
		comment = (qi::lit("//") >> *(qi::char_))[qi::_val = nullptr];

		//This block defines rules for parsing state formulas
		stateFormula %= orFormula;
		stateFormula.name("state formula");
		orFormula = andFormula[qi::_val = qi::_1] > *(qi::lit("|") > andFormula)[qi::_val =
				phoenix::new_<storm::property::csl::Or<double>>(qi::_val, qi::_1)];
		orFormula.name("or formula");
		andFormula = notFormula[qi::_val = qi::_1] > *(qi::lit("&") > notFormula)[qi::_val =
				phoenix::new_<storm::property::csl::And<double>>(qi::_val, qi::_1)];
		andFormula.name("and formula");
		notFormula = atomicStateFormula[qi::_val = qi::_1] | (qi::lit("!") > atomicStateFormula)[qi::_val =
				phoenix::new_<storm::property::csl::Not<double>>(qi::_1)];
		notFormula.name("not formula");

		//This block defines rules for "atomic" state formulas
		//(Propositions, probabilistic/reward formulas, and state formulas in brackets)
		atomicStateFormula %= probabilisticBoundOperator | steadyStateBoundOperator | atomicProposition | qi::lit("(") >> stateFormula >> qi::lit(")");
		atomicStateFormula.name("atomic state formula");
		atomicProposition = (freeIdentifierName)[qi::_val =
				phoenix::new_<storm::property::csl::Ap<double>>(qi::_1)];
		atomicProposition.name("atomic proposition");
		probabilisticBoundOperator = (
				(qi::lit("P") >> comparisonType > qi::double_ > qi::lit("[") > pathFormula > qi::lit("]"))[qi::_val =
						phoenix::new_<storm::property::csl::ProbabilisticBoundOperator<double> >(qi::_1, qi::_2, qi::_3)]
				);
		probabilisticBoundOperator.name("probabilistic bound operator");
		steadyStateBoundOperator = (
				(qi::lit("S") >> comparisonType > qi::double_ > qi::lit("[") > stateFormula > qi::lit("]"))[qi::_val =
										phoenix::new_<storm::property::csl::SteadyStateBoundOperator<double> >(qi::_1, qi::_2, qi::_3)]
				);
		steadyStateBoundOperator.name("steady state bound operator");

		//This block defines rules for parsing probabilistic path formulas
		pathFormula = (timeBoundedEventually | eventually | globally | next | timeBoundedUntil | until);
		pathFormula.name("path formula");
		timeBoundedEventually = (
				(qi::lit("F") >> qi::lit("[") > qi::double_ > qi::lit(",") > qi::double_ > qi::lit("]") > stateFormula)[qi::_val =
				phoenix::new_<storm::property::csl::TimeBoundedEventually<double>>(qi::_1, qi::_2, qi::_3)] |
				(qi::lit("F") >> (qi::lit("<=") | qi::lit("<")) > qi::double_ > stateFormula)[qi::_val =
								phoenix::new_<storm::property::csl::TimeBoundedEventually<double>>(0, qi::_1, qi::_2)] |
				(qi::lit("F") >> (qi::lit(">=") | qi::lit(">")) > qi::double_ > stateFormula)[qi::_val =
								phoenix::new_<storm::property::csl::TimeBoundedEventually<double>>(qi::_1, std::numeric_limits<double>::infinity(), qi::_2)]
				);
		timeBoundedEventually.name("time bounded eventually");
		eventually = (qi::lit("F") > stateFormula)[qi::_val =
				phoenix::new_<storm::property::csl::Eventually<double> >(qi::_1)];
		eventually.name("eventually");
		next = (qi::lit("X") > stateFormula)[qi::_val =
				phoenix::new_<storm::property::csl::Next<double> >(qi::_1)];
		next.name("next");
		globally = (qi::lit("G") > stateFormula)[qi::_val =
				phoenix::new_<storm::property::csl::Globally<double> >(qi::_1)];
		globally.name("globally");
		timeBoundedUntil = (
					(stateFormula[qi::_a = phoenix::construct<std::shared_ptr<storm::property::csl::AbstractStateFormula<double>>>(qi::_1)] >> qi::lit("U") >> qi::lit("[") > qi::double_ > qi::lit(",") > qi::double_ > qi::lit("]")  > stateFormula)
					[qi::_val = phoenix::new_<storm::property::csl::TimeBoundedUntil<double>>(qi::_2, qi::_3, phoenix::bind(&storm::property::csl::AbstractStateFormula<double>::clone, phoenix::bind(&std::shared_ptr<storm::property::csl::AbstractStateFormula<double>>::get, qi::_a)), qi::_4)] |
					(stateFormula[qi::_a = phoenix::construct<std::shared_ptr<storm::property::csl::AbstractStateFormula<double>>>(qi::_1)] >> qi::lit("U") >> (qi::lit("<=") | qi::lit("<")) > qi::double_  > stateFormula)
					[qi::_val = phoenix::new_<storm::property::csl::TimeBoundedUntil<double>>(0, qi::_2, phoenix::bind(&storm::property::csl::AbstractStateFormula<double>::clone, phoenix::bind(&std::shared_ptr<storm::property::csl::AbstractStateFormula<double>>::get, qi::_a)), qi::_3)] |
					(stateFormula[qi::_a = phoenix::construct<std::shared_ptr<storm::property::csl::AbstractStateFormula<double>>>(qi::_1)] >> qi::lit("U") >> (qi::lit(">=") | qi::lit(">")) > qi::double_  > stateFormula)
					[qi::_val = phoenix::new_<storm::property::csl::TimeBoundedUntil<double>>(qi::_2, std::numeric_limits<double>::infinity(), phoenix::bind(&storm::property::csl::AbstractStateFormula<double>::clone, phoenix::bind(&std::shared_ptr<storm::property::csl::AbstractStateFormula<double>>::get, qi::_a)), qi::_3)]
				);
		timeBoundedUntil.name("time bounded until");
		until = (stateFormula[qi::_a = phoenix::construct<std::shared_ptr<storm::property::csl::AbstractStateFormula<double>>>(qi::_1)] >> qi::lit("U") > stateFormula)[qi::_val =
				phoenix::new_<storm::property::csl::Until<double>>(phoenix::bind(&storm::property::csl::AbstractStateFormula<double>::clone, phoenix::bind(&std::shared_ptr<storm::property::csl::AbstractStateFormula<double>>::get, qi::_a)), qi::_2)];
		until.name("until formula");

		formula = (pathFormula | stateFormula);
		formula.name("CSL formula");

		//This block defines rules for parsing formulas with noBoundOperators
		noBoundOperator = (probabilisticNoBoundOperator | steadyStateNoBoundOperator);
		noBoundOperator.name("no bound operator");
		probabilisticNoBoundOperator =
				(qi::lit("P") >> qi::lit("min") >> qi::lit("=") >> qi::lit("?") >> qi::lit("[") >> pathFormula >> qi::lit("]"))[qi::_val =
						phoenix::new_<storm::property::csl::CslFilter<double>>(qi::_1, storm::property::MINIMIZE)] |
				(qi::lit("P") >> qi::lit("max") >> qi::lit("=") >> qi::lit("?") >> qi::lit("[") >> pathFormula >> qi::lit("]"))[qi::_val =
						phoenix::new_<storm::property::csl::CslFilter<double>>(qi::_1, storm::property::MAXIMIZE)] |
				(qi::lit("P") >> qi::lit("=") >> qi::lit("?") >> qi::lit("[") >> pathFormula >> qi::lit("]"))[qi::_val =
						phoenix::new_<storm::property::csl::CslFilter<double>>(qi::_1)];
		probabilisticNoBoundOperator.name("probabilistic no bound operator");
		steadyStateNoBoundOperator = (qi::lit("S") >> qi::lit("=") >> qi::lit("?") >> qi::lit("[") >> stateFormula >> qi::lit("]"))[qi::_val =
				phoenix::new_<storm::property::csl::CslFilter<double>>(qi::_1, storm::property::UNDEFINED, true)];
		steadyStateNoBoundOperator.name("steady state no bound operator");

		// This block defines rules for parsing filter actions.
		boundAction = (qi::lit("bound") > qi::lit("(") >> comparisonType >> qi::lit(",") >> qi::double_ >> qi::lit(")"))[qi::_val =
				        phoenix::new_<storm::property::action::BoundAction<double>>(qi::_1, qi::_2)];
		boundAction.name("bound action");

		invertAction = qi::lit("invert")[qi::_val = phoenix::new_<storm::property::action::InvertAction<double>>()];
		invertAction.name("invert action");

		formulaAction = (qi::lit("formula") > qi::lit("(") >> stateFormula >> qi::lit(")"))[qi::_val =
						phoenix::new_<storm::property::action::FormulaAction<double>>(qi::_1)];
		formulaAction.name("formula action");

		rangeAction = (
				(qi::lit("range") >> qi::lit("(") >> qi::uint_ >> qi::lit(",") > qi::uint_ >> qi::lit(")"))[qi::_val =
						phoenix::new_<storm::property::action::RangeAction<double>>(qi::_1, qi::_2)] |
				(qi::lit("range") >> qi::lit("(") >> qi::uint_ >> qi::lit(")"))[qi::_val =
						phoenix::new_<storm::property::action::RangeAction<double>>(qi::_1, qi::_1 + 1)]
				);
		rangeAction.name("range action");

		sortAction = (
				(qi::lit("sort") > qi::lit("(") >> sortingCategory >> qi::lit(")"))[qi::_val =
						phoenix::new_<storm::property::action::SortAction<double>>(qi::_1)] |
				(qi::lit("sort") > qi::lit("(") >> sortingCategory >> qi::lit(", ") >> qi::lit("asc") > qi::lit(")"))[qi::_val =
						phoenix::new_<storm::property::action::SortAction<double>>(qi::_1, true)] |
				(qi::lit("sort") > qi::lit("(") >> sortingCategory >> qi::lit(", ") >> qi::lit("desc") > qi::lit(")"))[qi::_val =
						phoenix::new_<storm::property::action::SortAction<double>>(qi::_1, false)]
				);
		sortAction.name("sort action");

		abstractAction = (boundAction | invertAction | formulaAction | rangeAction | sortAction) >> (qi::eps | qi::lit(";"));
		abstractAction.name("filter action");

		filter = (qi::lit("filter") >> qi::lit("[") >> +abstractAction >> qi::lit("]") >> qi::lit("(") >> formula >> qi::lit(")"))[qi::_val =
					phoenix::new_<storm::property::csl::CslFilter<double>>(qi::_2, qi::_1)] |
				 (noBoundOperator)[qi::_val =
					qi::_1] |
				 (formula)[qi::_val =
					phoenix::new_<storm::property::csl::CslFilter<double>>(qi::_1)];

		filter.name("CSL formula filter");

		start = (((filter) > (comment | qi::eps))[qi::_val = qi::_1] | comment[qi::_val = nullptr] ) > qi::eoi;
		start.name("CSL formula filter start");

	}

	qi::rule<Iterator, storm::property::csl::CslFilter<double>*(), Skipper> start;
	qi::rule<Iterator, storm::property::csl::CslFilter<double>*(), Skipper> filter;

	qi::rule<Iterator, storm::property::csl::CslFilter<double>*(), Skipper> noBoundOperator;
	qi::rule<Iterator, storm::property::csl::CslFilter<double>*(), Skipper> probabilisticNoBoundOperator;
	qi::rule<Iterator, storm::property::csl::CslFilter<double>*(), Skipper> steadyStateNoBoundOperator;

	qi::rule<Iterator, storm::property::action::AbstractAction<double>*(), Skipper> abstractAction;
	qi::rule<Iterator, storm::property::action::BoundAction<double>*(), Skipper> boundAction;
	qi::rule<Iterator, storm::property::action::InvertAction<double>*(), Skipper> invertAction;
	qi::rule<Iterator, storm::property::action::FormulaAction<double>*(), Skipper> formulaAction;
	qi::rule<Iterator, storm::property::action::RangeAction<double>*(), Skipper> rangeAction;
	qi::rule<Iterator, storm::property::action::SortAction<double>*(), Skipper> sortAction;

	qi::rule<Iterator, storm::property::csl::AbstractCslFormula<double>*(), Skipper> formula;
	qi::rule<Iterator, storm::property::csl::AbstractCslFormula<double>*(), Skipper> comment;

	qi::rule<Iterator, storm::property::csl::AbstractStateFormula<double>*(), Skipper> stateFormula;
	qi::rule<Iterator, storm::property::csl::AbstractStateFormula<double>*(), Skipper> atomicStateFormula;

	qi::rule<Iterator, storm::property::csl::AbstractStateFormula<double>*(), Skipper> andFormula;
	qi::rule<Iterator, storm::property::csl::AbstractStateFormula<double>*(), Skipper> atomicProposition;
	qi::rule<Iterator, storm::property::csl::AbstractStateFormula<double>*(), Skipper> orFormula;
	qi::rule<Iterator, storm::property::csl::AbstractStateFormula<double>*(), Skipper> notFormula;
	qi::rule<Iterator, storm::property::csl::ProbabilisticBoundOperator<double>*(), Skipper> probabilisticBoundOperator;
	qi::rule<Iterator, storm::property::csl::SteadyStateBoundOperator<double>*(), Skipper> steadyStateBoundOperator;

	qi::rule<Iterator, storm::property::csl::AbstractPathFormula<double>*(), Skipper> pathFormula;
	qi::rule<Iterator, storm::property::csl::TimeBoundedEventually<double>*(), Skipper> timeBoundedEventually;
	qi::rule<Iterator, storm::property::csl::Eventually<double>*(), Skipper> eventually;
	qi::rule<Iterator, storm::property::csl::Next<double>*(), Skipper> next;
	qi::rule<Iterator, storm::property::csl::Globally<double>*(), Skipper> globally;
	qi::rule<Iterator, storm::property::csl::TimeBoundedUntil<double>*(), qi::locals< std::shared_ptr<storm::property::csl::AbstractStateFormula<double>>>, Skipper> timeBoundedUntil;
	qi::rule<Iterator, storm::property::csl::Until<double>*(), qi::locals< std::shared_ptr<storm::property::csl::AbstractStateFormula<double>>>, Skipper> until;


	qi::rule<Iterator, std::string(), Skipper> freeIdentifierName;
	qi::rule<Iterator, storm::property::ComparisonType(), Skipper> comparisonType;
	qi::rule<Iterator, storm::property::action::SortAction<double>::SortingCategory(), Skipper> sortingCategory;

};

storm::property::csl::CslFilter<double>* CslParser::parseCslFormula(std::string formulaString) {
	// Prepare iterators to input.
	BaseIteratorType stringIteratorBegin = formulaString.begin();
	BaseIteratorType stringIteratorEnd = formulaString.end();
	PositionIteratorType positionIteratorBegin(stringIteratorBegin, stringIteratorEnd, formulaString);
	PositionIteratorType positionIteratorEnd;


	// Prepare resulting intermediate representation of input.
	storm::property::csl::CslFilter<double>* result_pointer = nullptr;

	CslGrammar<PositionIteratorType,  BOOST_TYPEOF(boost::spirit::ascii::space)> grammar;

	// Now, parse the formula from the given string
	try {
		qi::phrase_parse(positionIteratorBegin, positionIteratorEnd, grammar, boost::spirit::ascii::space, result_pointer);
	} catch(const qi::expectation_failure<PositionIteratorType>& e) {
		// If the parser expected content different than the one provided, display information
		// about the location of the error.
		const boost::spirit::classic::file_position_base<std::string>& pos = e.first.get_position();

		// Construct the error message including a caret display of the position in the
		// erroneous line.
		std::stringstream msg;
		msg << pos.file << ", line " << pos.line << ", column " << pos.column
				<< ": parse error: expected " << e.what_ << std::endl << "\t"
				<< e.first.get_currentline() << std::endl << "\t";
		int i = 0;
		for (i = 0; i < pos.column; ++i) {
			msg << "-";
		}
		msg << "^";
		for (; i < 80; ++i) {
			msg << "-";
		}
		msg << std::endl;

		std::cerr << msg.str();

		// Now propagate exception.
		throw storm::exceptions::WrongFormatException() << msg.str();
	}

	// The syntax can be so wrong that no rule can be matched at all
	// In that case, no expectation failure is thrown, but the parser just returns nullptr
	// Then, of course the result is not usable, hence we throw a WrongFormatException, too.
	if (result_pointer == nullptr) {
		throw storm::exceptions::WrongFormatException() << "Syntax error in formula";
	}

	return result_pointer;
}

} /* namespace parser */
} /* namespace storm */
