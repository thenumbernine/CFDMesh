#pragma once

#include <numeric>
#include <functional>
#include <map>
#include <regex>
#include <string>
#include <sstream>
#include <experimental/type_traits>	//is_detected_v


template<typename T>
std::ostream& ostreamForFields(std::ostream& a, const T& b) {
	a << "[";
	const char* sep = "";
	Common::TupleForEach(T::fields, [&a, &b, &sep](auto x, size_t i) constexpr {
		auto name = std::get<0>(x);
		auto field = std::get<1>(x);
		auto& value = b.*field;
		a << sep << name << " = " << value;
		sep = ", ";
	});
	return a << "]";
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

//TODO also in GUI.h
template<typename T>
using has_field_t = decltype(std::declval<T&>().fields);

//if a class has 'fields' then automatically use ostreamForFields as its serialization
template<
	typename T,
	std::enable_if_t<std::experimental::is_detected_v<has_field_t, T>, int> = 0
>
std::ostream& operator<<(std::ostream& o, const T& x) {
	return ostreamForFields(o, x);
}

namespace std {

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

}
