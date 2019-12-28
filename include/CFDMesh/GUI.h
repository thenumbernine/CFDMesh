#pragma once

#include "CFDMesh/Util.h"
#include "Common/Meta.h"
#include "cimgui.h"
#include <experimental/type_traits>	//is_detected_v

namespace CFDMesh {

template<typename T>
struct UpdateGUI {
	static void exec(T* ptr, std::string prefix);
};

template<typename T, typename Field>
struct UpdateGUIField {
	//static void exec(T* ptr, std::string prefix, const Field& x) {}
};

template<typename T, typename B>
struct UpdateGUIField<T, std::pair<const char*, B>> {
	static void exec(T* ptr, std::string prefix, const std::pair<const char*, B>& x) {
		// TODO make it a part of the template type declaration
		if constexpr (std::is_member_object_pointer_v<B>) {
			auto name = std::get<0>(x);
			auto field = std::get<1>(x);
			using FieldType = typename Common::MemberPointer<decltype(field)>::FieldType;
			UpdateGUI<FieldType>::exec(&(ptr->*field), prefix + name);
		}
	}
};

template<typename T>
struct UpdateGUIField<T, std::pair<const char*, void (T::*)()>> {
	static void exec(T* ptr, std::string prefix, const std::pair<const char*, void(T::*)()>& x) {
		auto name = std::get<0>(x);
		auto field = std::get<1>(x);
		if (igSmallButton((prefix + name).c_str())) {
			(ptr->*field)();
		}
	}
};

struct GUISeparator {};

template<typename T>
struct UpdateGUIField<T, GUISeparator> {
	static void exec(T* ptr, std::string prefix, const GUISeparator& x) {
		igSeparator();
	}
};

template<typename T>
struct GUICall {
	std::function<void(T*)> f;
	GUICall(std::function<void(T*)> f_) : f(f_) {}
};

//TODO not working
template<typename T>
struct UpdateGUIField<T, GUICall<T>> {
	static void exec(T* ptr, std::string prefix, GUICall<T> x) {
		x.f(ptr);
	}
};

struct GUIReadOnly {};

template<typename T, typename B>
struct UpdateGUIField<T, std::tuple<const char*, B, GUIReadOnly>> {
	static void exec(T* ptr, std::string prefix, const std::tuple<const char*, B, GUIReadOnly>& x) {
		auto name = std::get<0>(x);
		auto field = std::get<1>(x);
		igText((prefix + name + ": " + std::to_string(ptr->*field)).c_str());
	}
};

template<typename T>
void UpdateGUI<T>::exec(T* ptr, std::string prefix) {
	if constexpr (std::experimental::is_detected_v<has_fields_t, T>) {
		igPushIDPtr(ptr);
		Common::TupleForEach(T::fields, [ptr, &prefix](auto x, size_t i) constexpr {
			UpdateGUIField<T, decltype(x)>::exec(ptr, prefix, x);
		});
		igPopID();
	} else if constexpr (std::is_same_v<T, bool>) {
		igCheckbox(prefix.c_str(), ptr);
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

template<>
struct UpdateGUI<bool> {
	static void exec(bool* ptr, std::string prefix) {
		igCheckbox(prefix.c_str(), ptr);
	}
};

template<typename A, typename B>
struct UpdateGUI<std::pair<A, B>> {
	static void exec(std::pair<A, B>* ptr, std::string prefix) {
		UpdateGUI<A>::exec(&ptr->first, prefix + " 0");
		UpdateGUI<B>::exec(&ptr->second, prefix + " 1");
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
