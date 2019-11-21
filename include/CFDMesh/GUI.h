#pragma once

#include "cimgui.h"
#include <experimental/type_traits>	//is_detected_v

namespace CFDMesh {

template<typename T>
using has_field_t = decltype(std::declval<T&>().fields);

template<typename T>
void updateGUI(T* ptr, std::string prefix = {}) {
	if constexpr (std::experimental::is_detected_v<has_field_t, T>) {
		tuple_for_each(T::fields, [ptr, &prefix](auto x, size_t i) constexpr {
			auto name = x.first;
			auto field = x.second;
			using FieldType = typename MemberPointerInfo<decltype(field)>::FieldType;
			updateGUI<FieldType>(&(ptr->*field), prefix + name);
		});
	
	} else if constexpr (std::is_floating_point_v<T>) {
		float f = *ptr;
		if (igInputFloat(prefix.c_str(), &f, .1, 1, "%f", 0)) {
			*ptr = f;
		}
	
	//TODO handle all floatN, doubleN, etc types
	} else if constexpr (std::is_same_v<T, float3>) {
		updateGUI<float>(ptr->v+0, prefix + ".x");
		updateGUI<float>(ptr->v+1, prefix + ".y");
		updateGUI<float>(ptr->v+2, prefix + ".z");
	}
}

}
