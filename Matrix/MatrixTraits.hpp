#pragma once

namespace twmath{
	struct MatrixToken{};
	struct StaticMatrixToken : public MatrixToken{};
	struct DynamicMatrixToken : public MatrixToken{};
	
	// already defined in VectorTraits:
	// template<class M> struct value_type{using type = M;};
	// template<class M> using value_type_t = typename value_type<M>::type;
	
	template<class M> constexpr bool is_matrix_v = std::is_base_of_v<twmath::MatrixToken, M>;
	template<class M> constexpr bool is_static_matrix_v = std::is_base_of_v<twmath::StaticMatrixToken, M>;
	template<class M> constexpr bool is_dynamic_matrix_v = std::is_base_of_v<twmath::DynamicMatrixToken, M>;
}