#pragma once

#include "CFDMesh/Vector.h"
#include "CFDMesh/Util.h"
#include "Common/crtp_cast.h"
#include <vector>
#include <memory>
#include <string>

namespace CFDMesh {
namespace Equation {


template<typename T, typename Enable = void>
struct FloatTypeForType;

template<typename T>
struct FloatTypeForType<T, typename std::enable_if_t<std::is_floating_point_v<T>>> { 
	using Type = float;
};

template<typename T>
struct FloatTypeForType<Tensor::Vector<T, 3>, typename std::enable_if_t<std::is_floating_point_v<T>>> {
	using Type = float3;
};

template<
	typename real, 
	typename Cons_,
	typename Base
>
struct Equation {
	using real2 = Tensor::Vector<real, 2>;
	using real3 = Tensor::Vector<real, 3>;
	
	using Cons = Cons_;

	struct InitCond {
		virtual ~InitCond() {}
		virtual const char* name() const = 0;
		virtual Cons initCell(const Base* eqn, real3 pos) const = 0;
		virtual void updateGUI() {}
	};

	using DisplayFunc = std::function<float(const Base*, const Cons&)>;

	struct DisplayMethod {
		std::string name;
		DisplayFunc f;
		DisplayMethod(const std::string& name_, DisplayFunc f_) : name(name_), f(f_) {}
	};

	std::vector<std::shared_ptr<InitCond>> initConds;
	std::vector<const char*> initCondNames;

	std::vector<std::shared_ptr<DisplayMethod>> displayMethods;
	std::vector<const char*> displayMethodNames;

	Equation() {
		crtp_cast<Base>(*this).buildInitCondsAndDisplayVars();

		updateNames();
	}

	void addDisplayScalar(const std::string& name, DisplayFunc func) {
		displayMethods.push_back(std::make_shared<DisplayMethod>(name, func));
	}
	
	void addDisplayVector(const std::string& name, std::function<float3(const Base* eqn, const Cons& U)> func) {
		displayMethods.push_back(std::make_shared<DisplayMethod>(name + " len", [func](const Base* eqn, const Cons& U) -> float { return real3::length(func(eqn, U)); }));
		for (int i = 0; i < 3; ++i) {
			displayMethods.push_back(std::make_shared<DisplayMethod>(name + std::to_string(i), [i, func](const Base* eqn, const Cons& U) -> float { return func(eqn, U)(i); }));
		}
	}

	template<typename Type>
	void addDisplayForType(const std::string& name, std::function<typename FloatTypeForType<Type>::Type(const Base*, const Cons&)> func) {
		if constexpr (std::is_same_v<Type, real>) {
			addDisplayScalar(name, func);
		} else if constexpr (std::is_same_v<Type, real3>) {
			addDisplayVector(name, func);
		}
	}

	void updateNames() {
		initCondNames = map<
			decltype(initConds),
			std::vector<const char*>
		>(
			initConds,
			[](std::shared_ptr<InitCond> ic) -> const char* { return ic->name(); }
		);
		
		displayMethodNames = map<
			decltype(displayMethods),
			std::vector<const char*>
		>(
			displayMethods,
			[](const std::shared_ptr<DisplayMethod>& m) -> const char* { return m->name.c_str(); }
		);
	}
};

}
}
