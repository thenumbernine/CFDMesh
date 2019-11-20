#pragma once

#include "cimgui.h"
#include <experimental/type_traits>	//is_detected_v

namespace CFDMesh {


template<typename T>
using helper_field_t = decltype(std::declval<T&>().fields);

template<typename T>
void updateGUIForFields(T* ptr, std::string prefix = {}) {
	if constexpr (std::experimental::is_detected_v<helper_field_t, T>) {
		tuple_for_each(T::fields, [ptr, &prefix](auto x, size_t i) constexpr {
			auto name = x.first;
			auto field = x.second;
			using FieldType = typename MemberPointerInfo<decltype(field)>::FieldType;
			updateGUIForFields<FieldType>(&(ptr->*field), prefix + name);
		});
	
	} else if constexpr (std::is_floating_point_v<T>) {
		float f = *ptr;
		if (igInputFloat(prefix.c_str(), &f, .1, 1, "%f", 0)) {
			*ptr = f;
		}
	
	//TODO handle all floatN, doubleN, etc types
	} else if constexpr (std::is_same_v<T, float3>) {
		updateGUIForFields<float>(ptr->v+0, prefix + ".x");
		updateGUIForFields<float>(ptr->v+1, prefix + ".y");
		updateGUIForFields<float>(ptr->v+2, prefix + ".z");
	}
}


}
