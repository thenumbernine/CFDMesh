#pragma once

#include <numeric>
#include <functional>
#include <map>
#include <regex>
#include <string>
#include <sstream>
#include <chrono>
#include <experimental/type_traits>	//is_detected_v

template<typename T>
using has_fields_t = decltype(T::fields);

//if a typename has 'fields'...
template<typename T>
inline std::enable_if_t<std::experimental::is_detected_v<has_fields_t, T>, std::ostream&>
operator<<(std::ostream& o, const T& b) {
	o << "[";
	Common::TupleForEach(T::fields, [&o, &b](const auto& x, size_t i) constexpr {
		const auto& name = std::get<0>(x);
		const auto& field = std::get<1>(x);
		if (i > 0) o << ", ";
		o << name << "=" << b.*field;
	});
	o << "]";
	return o;
}

template<typename T>
std::string objectStringFromOStream(const T& x) {
	std::ostringstream ss;
	ss << x;
	return ss.str();
}

template<typename T>
std::ostream& operator<<(std::ostream& o, const std::vector<T>& v) {
	o << "[";
	const char* sep = "";
	for (const auto& x : v) {
		o << sep << x;
		sep = ", ";
	}
	return o << "]";
}

namespace std {

//if a typename has 'fields' then to_string uses objectStringFromOStream
template<typename T>
inline std::enable_if_t<std::experimental::is_detected_v<has_fields_t, T>, std::string>
to_string(const T& x) {
	return objectStringFromOStream(x);
}

const std::string& to_string(const std::string& s) {
	return s;
}
	
//TODO put this somewhere ...
template<typename T, int n>
std::string to_string(const Tensor::Vector<T, n>& x) {
	return objectStringFromOStream(x);
}

template<typename T>
std::string to_string(const std::vector<T>& x) {
	return objectStringFromOStream(x);
}

}

namespace CFDMesh {

//https://stackoverflow.com/questions/16749069/c-split-string-by-regex
template<typename T>
T split(const std::string &string_to_split, const std::string& regexPattern) {
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
std::string concat(const T& v, const std::string& sep) {
	bool first = true;
	std::string result = "";
	for (const auto& s : v) {
		if (!first) result += sep;
		result += std::to_string(s);
		first = false;
	}
	return result;
}

template<typename From, typename To>
To map(const From& from, std::function<typename To::value_type(typename From::value_type)> f) {
	To to;
	for (const auto& v : from) {
		to.push_back(f(v));
	}
	return to;
}

template<typename T>
typename T::value_type sum(const T& t) {
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
