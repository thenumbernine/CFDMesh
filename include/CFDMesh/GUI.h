#pragma once

#include "Common/Meta.h"
#include "cimgui.h"
#include <experimental/type_traits>	//is_detected_v

namespace CFDMesh {

template<typename T>
using has_field_t = decltype(std::declval<T&>().fields);

template<typename T>
struct UpdateGUI {
	static void exec(T* ptr, std::string prefix) {
		if constexpr (std::experimental::is_detected_v<has_field_t, T>) {
			Common::TupleForEach(T::fields, [ptr, &prefix](auto x, size_t i) constexpr {
				auto name = x.first;
				auto field = x.second;
				using FieldType = typename Common::MemberPointer<decltype(field)>::FieldType;
				UpdateGUI<FieldType>::exec(&(ptr->*field), prefix + name);
			});
		} else if constexpr (std::is_floating_point_v<T>) {
			float f = *ptr;
			if (igInputFloat(prefix.c_str(), &f, .1, 1, "%f", 0)) {
				*ptr = f;
			}
		} else if constexpr (std::is_integral_v<T>) {
			int i = *ptr;
			if (igInputInt(prefix.c_str(), &i, 1, 10, 0)) {
				*ptr = i;
			}
		} else {
			igText(prefix.c_str());
		}
	}
};

template<typename T>
struct UpdateGUI<Tensor::Vector<T, 2>> {
	static void exec(Tensor::Vector<T, 2>* ptr, std::string prefix) {
		UpdateGUI<T>::exec(ptr->v+0, prefix + " 0");
		UpdateGUI<T>::exec(ptr->v+1, prefix + " 1");
	}
};

template<typename T>
struct UpdateGUI<Tensor::Vector<T, 3>> {
	static void exec(Tensor::Vector<T, 3>* ptr, std::string prefix) {
		UpdateGUI<T>::exec(ptr->v+0, prefix + " 0");
		UpdateGUI<T>::exec(ptr->v+1, prefix + " 1");
		UpdateGUI<T>::exec(ptr->v+2, prefix + " 2");
	}
};

template<typename T>
void updateGUI(T* ptr, std::string prefix = {}) {
	UpdateGUI<T>::exec(ptr, prefix);
}

}
