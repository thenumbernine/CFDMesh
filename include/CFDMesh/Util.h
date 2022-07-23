#pragma once

#include <numeric>
#include <functional>
#include <map>
#include <regex>
#include <string>
#include <sstream>
#include <chrono>

template<typename T>
constexpr bool has_fields_v = requires(T const & t) { t.fields; };

//if a typename has 'fields'...
template<typename T> requires has_fields_v<T>
inline std::ostream& operator<<(std::ostream& o, T const & b) {
	o << "[";
	Common::TupleForEach(T::fields, [&o, &b](auto const & x, size_t i) constexpr -> bool {
		auto const & name = std::get<0>(x);
		auto const & field = std::get<1>(x);
		if (i > 0) o << ", ";
		o << name << "=" << b.*field;
		return false;
	});
	o << "]";
	return o;
}

template<typename T>
std::string objectStringFromOStream(T const & x) {
	std::ostringstream ss;
	ss << x;
	return ss.str();
}

template<typename T>
std::ostream& operator<<(std::ostream& o, std::vector<T> const & v) {
	o << "[";
	char const * sep = "";
	for (auto const & x : v) {
		o << sep << x;
		sep = ", ";
	}
	return o << "]";
}

namespace std {

//if a typename has 'fields' then to_string uses objectStringFromOStream
template<typename T> requires has_fields_v<T>
inline std::string to_string(T const & x) {
	return objectStringFromOStream(x);
}

std::string const & to_string(std::string const & s) {
	return s;
}
	
//TODO put this somewhere ...
template<typename T, int n>
std::string to_string(Tensor::Vector<T, n> const & x) {
	return objectStringFromOStream(x);
}

template<typename T>
std::string to_string(std::vector<T> const & x) {
	return objectStringFromOStream(x);
}

}

namespace CFDMesh {

//https://stackoverflow.com/questions/16749069/c-split-string-by-regex
template<typename T>
T split(std::string const & string_to_split, std::string const & regexPattern) {
	std::regex rgx(regexPattern);
	std::sregex_token_iterator iter(string_to_split.begin(),
		string_to_split.end(),
		rgx,
		-1);
	std::sregex_token_iterator end;
	T result;
	for ( ; iter != end; ++iter) {
		result.push_back(*iter);
	}
	return result;
}

template<typename T>
std::string concat(T const & v, std::string const & sep) {
	bool first = true;
	std::string result = "";
	for (auto const & s : v) {
		if (!first) result += sep;
		result += std::to_string(s);
		first = false;
	}
	return result;
}

#warning TODO use Common/Meta.h vectorMap
template<typename From, typename To>
To map(From const & from, std::function<typename To::value_type(typename From::value_type)> f) {
	To to;
	for (auto const & v : from) {
		to.push_back(f(v));
	}
	return to;
}

template<typename T>
typename T::value_type sum(T const & t) {
	return std::accumulate(t.begin(), t.end(), typename T::value_type());
}

inline void timeFunc(std::string name, std::function<void()> cb) {
	std::cout << name << "..." << std::endl;
	std::chrono::high_resolution_clock::time_point startTime = std::chrono::high_resolution_clock::now();
	cb();
	std::chrono::high_resolution_clock::time_point endTime = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> dt = std::chrono::duration_cast<std::chrono::duration<double>>(endTime - startTime);
	std::cout << "...done " << name << " (" << dt.count() << "s)" << std::endl;
}

}
