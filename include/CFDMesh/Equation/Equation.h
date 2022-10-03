#pragma once

#include "CFDMesh/Vector.h"
#include "CFDMesh/Mesh/Cell.h"
#include "CFDMesh/Mesh/Face.h"
#include "Common/String.h"
#include "Common/Exception.h"
#include "Common/crtp_cast.h"
#include <vector>
#include <memory>
#include <string>

namespace CFDMesh {
namespace Equation {


// reduce any floating point scalar/vector to 32-bit:

template<typename T>
struct FloatTypeForType;

template<typename T>
requires (std::is_floating_point_v<T>)
struct FloatTypeForType<T> { 
	using Type = float;
};

template<typename T>
requires (std::is_floating_point_v<T>)
struct FloatTypeForType<Tensor::_vec<T, 3>> {
	using Type = float3;
};


template<
	typename Base,
	typename real, 
	typename Cons_,
	typename Prim_
>
struct Equation {
	using real2 = Tensor::_vec<real, 2>;
	using real3 = Tensor::_vec<real, 3>;
	using real3x3 = Tensor::_mat<real, 3, 3>;
	
	using Cons = Cons_;
	using Prim = Prim_;
	
	using Cell = Mesh::Cell<real, Cons>;
	using Face = Mesh::Face<real, Cons>;
	
	enum { numStates = Cons::size };
	enum { numIntStates = Cons::size };

	struct InitCond {
		virtual ~InitCond() {}
		virtual const char* name() const = 0;
		virtual Cons initCell(const Base* eqn, real3 pos) const = 0;
		virtual void updateGUI() {}
	};

	using DisplayFunc = std::function<float(const Base*, const Cell*)>;

	struct DisplayMethod {
		std::string name;
		DisplayFunc f;
		DisplayMethod(const std::string& name_, DisplayFunc f_) : name(name_), f(f_) {}
	};

	std::vector<std::shared_ptr<InitCond>> initConds;
	std::vector<const char*> initCondNames;

	std::vector<std::shared_ptr<DisplayMethod>> displayMethods;

	Equation() {
		crtp_cast<Base>(*this).buildInitCondsAndDisplayVars();

		//cycle through the Cons::fields tuple and add each of these
		//clang++ has a bug?  seems removed 'this' and you get an error, cuz 'this' needs to be captured, but keep 'this' and clang++ gives a warning that 'this' isn't needed ...
		Common::TupleForEach(Cons::fields, [this](auto x, size_t i) constexpr -> bool {
			auto field = std::get<1>(x);
			using FieldType = typename Common::MemberPointer<decltype(field)>::FieldType;
			addDisplayForType<FieldType>(
				std::get<0>(x),
				[field](const Base* eqn, const Cell* c) -> typename FloatTypeForType<FieldType>::Type { 
					const auto& U = c->U;
					return U.*field; 
				}
			);
			return false;
		});

		if constexpr (!std::is_same_v<Cons, Prim>) {
			Common::TupleForEach(Prim::fields, [this](auto x, size_t i) constexpr -> bool {
				auto field = std::get<1>(x);
				using FieldType = typename Common::MemberPointer<decltype(field)>::FieldType;
				addDisplayForType<FieldType>(
					std::get<0>(x),
					[field](const Base* eqn, const Cell* c) -> typename FloatTypeForType<FieldType>::Type { 
						const auto& U = c->U;
						auto W = eqn->primFromCons(U);
						return W.*field; 
					}
				);
				return false;
			});	
		}

		initCondNames = Common::mapElems<
			decltype(initConds),
			std::vector<const char*>
		>(
			initConds,
			[](std::shared_ptr<InitCond> ic) -> const char* { return ic->name(); }
		);
	}

	void addDisplayScalar(const std::string& name, DisplayFunc func) {
		displayMethods.push_back(std::make_shared<DisplayMethod>(name, func));
	}
	
	void addDisplayVector(const std::string& name, std::function<float3(const Base* eqn, const Cell* c)> func) {
		displayMethods.push_back(std::make_shared<DisplayMethod>(name + " len", [func](const Base* eqn, const Cell* c) -> float { return func(eqn, c).length(); }));
		for (int i = 0; i < 3; ++i) {
			displayMethods.push_back(std::make_shared<DisplayMethod>(name + std::to_string(i), [i, func](const Base* eqn, const Cell* c) -> float { return func(eqn, c)(i); }));
		}
	}

	template<typename Type>
	void addDisplayForType(const std::string& name, std::function<typename FloatTypeForType<Type>::Type(const Base*, const Cell*)> func) {
		if constexpr (std::is_same_v<Type, real>) {
			addDisplayScalar(name, func);
		} else if constexpr (std::is_same_v<Type, real3>) {
			addDisplayVector(name, func);
		}
	}

	Cons rotateFrom(Cons U, const real3x3& nb) const {
		Common::TupleForEach(Cons::fields, [&U, &nb](auto x, size_t i) constexpr -> bool {
			auto field = std::get<1>(x);
			using FieldType = typename Common::MemberPointer<decltype(field)>::FieldType;
			if constexpr (std::is_same_v<FieldType, real>) {
				//real = unchanged
			} else if constexpr (std::is_same_v<FieldType, real3>) {
				//real3 = multiply
				//nb is stored transposed, and rotateFrom is an inverse rotation
				//U.*field = CFDMesh::rotateFrom<real3>(U.*field, normal);
				real3 src = U.*field;
				for (int i = 0; i < 3; ++i) {
					real sum = 0;
					for (int j = 0; j < 3; ++j) {
						sum += nb(i,j) * src(j);
					}
					(U.*field)(i) = sum;
				}
			} else {
				throw Common::Exception() << "TODO implement this rotateFrom for this type";
				//static_assert(false, "TODO implement this rotateFrom for this type");
			}
			return false;
		});
		return U;
	}

	Cons rotateTo(Cons U, const real3x3& nb) const {
		Common::TupleForEach(Cons::fields, [&U, &nb](auto x, size_t i) constexpr -> bool {
			auto field = std::get<1>(x);
			using FieldType = typename Common::MemberPointer<decltype(field)>::FieldType;
			if constexpr (std::is_same_v<FieldType, real>) {
				//real = unchanged
			} else if constexpr (std::is_same_v<FieldType, real3>) {
				//real3 = multiply
				//nb is stored transposed, and rotateFrom is a forward rotation
				//U.*field = CFDMesh::rotateTo<real3>(U.*field, normal);
				real3 src = U.*field;
				for (int i = 0; i < 3; ++i) {
					real sum = 0;
					for (int j = 0; j < 3; ++j) {
						sum += nb(j,i) * src(j);
					}
					(U.*field)(i) = sum;
				}
			} else {
				throw Common::Exception() << "TODO implement this rotateFrom for this type";
				//static_assert(false, "TODO implement this rotateFrom for this type");
			}
			return false;
		});
		return U;
	}

	//TODO - restitution isn't the only solution
	// at the moment if an edge is only connected to one face then it will not update
	// and that means if there are constant values on both faces of the other edges ... 
	Cons reflect(Cons U, const real3& normal, real restitution) {
		Common::TupleForEach(Cons::fields, [&U, &normal, restitution](auto x, size_t i) constexpr -> bool {
			auto field = std::get<1>(x);
			using FieldType = typename Common::MemberPointer<decltype(field)>::FieldType;
			if constexpr (std::is_same_v<FieldType, real3>) {
				U.*field = U.*field - normal * ((1 + restitution) * Tensor::dot(normal, U.*field));
			}
			return false;
		});
		return U;
	}
};

}
}
