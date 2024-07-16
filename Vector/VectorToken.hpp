#pragma once

namespace twmath{
	struct VectorToken{};
	struct StaticVectorToken : public VectorToken{};
	struct DynamicVectorToken : public VectorToken{};
	
	template<class V> struct value_type{using type = V;};
	template<class V> using value_type_t = typename value_type<V>::type;
	
	
	template<class V> constexpr bool is_vector_v = std::is_base_of_v<V, VectorToken>;
	template<class V> constexpr bool is_static_vector_v = std::is_base_of_v<V, StaticVectorToken>;
	template<class V> constexpr bool is_dynamic_vector_v = std::is_base_of_v<V, DynamicVectorToken>;
}