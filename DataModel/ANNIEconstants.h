#ifndef ANNIECONSTANTS_H
#define ANNIECONSTANTS_H

#include <math.h>
#include <regex>
#include <boost/regex.hpp> // since std::regex doesn't work
#include <boost/regex/pattern_except.hpp>

/// @brief The impedance (in Ohms) to assume when computing charge values for
/// calibrated ADC hits
constexpr double ADC_IMPEDANCE = 50.; // Ohm

/// @brief Multiplying by this constant converts ADC counts to Volts
const double ADC_TO_VOLT = 2.415 / pow(2., 12);

/// @brief The number of nanoseconds per ADC sample
constexpr unsigned int NS_PER_ADC_SAMPLE = 2; // ns

/// @brief A dummy value to use to initialize integers
constexpr int BOGUS_INT = -9999;

/// @brief The format code for a multievent binary BoostStore
constexpr int BOOST_STORE_BINARY_FORMAT = 0;
constexpr int BOOST_STORE_ASCII_FORMAT = 1;
constexpr int BOOST_STORE_MULTIEVENT_FORMAT = 2;

// Used to convert from seconds and nanoseconds to milliseconds when
// creating timestamps in ms while interacting with the Intensity Frontier
// beam database via the BeamChecker and BeamFetcher tools
constexpr unsigned long long THOUSAND = 1000ull;
constexpr unsigned long long MILLION = 1000000ull;
constexpr unsigned long long BILLION = 1000000000ull;

// constants for mapping regular expression error codes to descriptive strings
// ---------------------------------------------------------------------------
const std::map<std::regex_constants::error_type,std::string> regex_err_strings{
	{std::regex_constants::error_collate, "The expression contained an invalid collating element name."},
	{std::regex_constants::error_ctype, "The expression contained an invalid character class name."},
	{std::regex_constants::error_escape, "The expression contained an invalid escaped character, or a trailing escape."},
	{std::regex_constants::error_backref, "The expression contained an invalid back reference."},
	{std::regex_constants::error_brack, "The expression contained mismatched brackets ([ and ])."},
	{std::regex_constants::error_paren, "The expression contained mismatched parentheses (( and ))."},
	{std::regex_constants::error_brace, "The expression contained mismatched braces ({ and })."},
	{std::regex_constants::error_badbrace, "The expression contained an invalid range between braces ({ and })."},
	{std::regex_constants::error_range, "The expression contained an invalid character range."},
	{std::regex_constants::error_space, "There was insufficient memory to convert the expression into a finite state machine."},
	{std::regex_constants::error_badrepeat, "The expression contained a repeat specifier (one of *?+{) that was not preceded by a valid regular expression."},
	{std::regex_constants::error_complexity, "The complexity of an attempted match against a regular expression exceeded a pre-set level."},
	{std::regex_constants::error_stack, "There was insufficient memory to determine whether the regular expression could match the specified character sequence."}
};

// boost version
const std::map<boost::regex_constants::error_type,std::string> bregex_err_strings{
	{boost::regex_constants::error_collate, "The expression contained an invalid collating element name."},
	{boost::regex_constants::error_ctype, "The expression contained an invalid character class name."},
	{boost::regex_constants::error_escape, "The expression contained an invalid escaped character, or a trailing escape."},
	{boost::regex_constants::error_backref, "The expression contained an invalid back reference."},
	{boost::regex_constants::error_brack, "The expression contained mismatched brackets ([ and ])."},
	{boost::regex_constants::error_paren, "The expression contained mismatched parentheses (( and ))."},
	{boost::regex_constants::error_brace, "The expression contained mismatched braces ({ and })."},
	{boost::regex_constants::error_badbrace, "The expression contained an invalid range between braces ({ and })."},
	{boost::regex_constants::error_range, "The expression contained an invalid character range."},
	{boost::regex_constants::error_space, "There was insufficient memory to convert the expression into a finite state machine."},
	{boost::regex_constants::error_badrepeat, "The expression contained a repeat specifier (one of *?+{) that was not preceded by a valid regular expression."},
	{boost::regex_constants::error_complexity, "The complexity of an attempted match against a regular expression exceeded a pre-set level."},
	{boost::regex_constants::error_stack, "There was insufficient memory to determine whether the regular expression could match the specified character sequence."},
	{boost::regex_constants::error_bad_pattern, "Invalid regex pattern."}
};

#endif
