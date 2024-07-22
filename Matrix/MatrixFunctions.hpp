#pragma once


#include "definitions.hpp"
#include "BaseTypeArithmetic.hpp"

#include "MatrixTraits.hpp"
#include "MatrixRef.hpp"
#include "Matrix.hpp"


namespace twmath{
	
	// --------------------- transpose -----------------------
	
	template<class M, TWMATH_ENABLE_IF(is_static_matrix_v<M>)>
	constexpr MatrixRef<value_type_t<M>, M::smem_rows(), M::smem_columns(), !M::stransposed()> transpose(M& m) {
		MatrixRef<value_type_t<M>, M::smem_rows(), M::smem_columns(), !M::stransposed()> result(m.data(), m.row_increment(), m.column_increment());
		return result;
	}
	
	template<class M, TWMATH_ENABLE_IF(is_static_matrix_v<M>)>
	constexpr ConstMatrixRef<value_type_t<M>, M::smem_rows(), M::smem_columns(), !M::stransposed()> transpose(const M& m) {
		ConstMatrixRef<value_type_t<M>, M::smem_rows(), M::smem_columns(), !M::stransposed()> result(m.cdata(), m.row_increment(), m.column_increment());
		return result;
	}

	// --------------------- for_each -----------------------
	
	template<class Ml, class Mr, class F, 
		TWMATH_ENABLE_IF(is_static_matrix_v<Ml> && is_static_matrix_v<Mr> 
		&& Ml::srows() == Mr::srows() && Ml::scolumns() == Mr::scolumns())>
	constexpr Matrix<F, Ml::srows(), Ml::scolumns()> for_each (const Ml& ml, const Mr& mr, F (*f)(const value_type_t<Ml>&, const value_type_t<Mr>&)){
		Matrix<F, Ml::srows(), Ml::scolumns()> result;
		for(size_t row = 0; row != Ml::srows(); ++row){
			for(size_t col = 0; col != Ml::scolumns(); ++col){
				result.at(row, col) = f(ml.at(row, col), mr.at(row, col));
			}
		}
		return result;
	}
	
	template<class Tl, class Mr, class F, TWMATH_ENABLE_IF(!is_matrix_v<Tl> && !is_vector_v<Tl> && is_static_matrix_v<Mr>)>
	constexpr Matrix<F, Mr::srows(), Mr::scolumns()> for_each (const Tl& tl, const Mr& mr, F (*f)(const Tl&, const value_type_t<Mr>&)){
		Matrix<F, Mr::srows(), Mr::scolumns()> result;
		for(size_t row = 0; row != Mr::srows(); ++row){
			for(size_t col = 0; col != Mr::scolumns(); ++col){
				result.at(row, col) = f(tl, mr.at(row, col));
			}
		}
		return result;
	}
	template<class Ml, class Tr, class F, TWMATH_ENABLE_IF(is_static_matrix_v<Ml> && !is_matrix_v<Tr> && !is_vector_v<Tr>)>
	constexpr Matrix<F, Ml::srows(), Ml::scolumns()>  for_each (const Ml& ml, const Tr& tr, F (*f)(const value_type_t<Ml>&, const Tr&)){
		Matrix<F, Ml::srows(), Ml::scolumns()> result;
		for(size_t row = 0; row != Ml::srows(); ++row){
			for(size_t col = 0; col != Ml::scolumns(); ++col){
				result.at(row, col) = f(ml.at(row, col), tr);
			}
		}
		return result;
	}
	template<class M, class F, TWMATH_ENABLE_IF(is_static_matrix_v<M>)>
	constexpr Matrix<F, M::srows(), M::scolumns()> for_each (const M& m, F (*f)(const value_type_t<M>&)){
		Matrix<F, M::srows(), M::scolumns()> result;
		for(size_t row = 0; row != M::srows(); ++row){
			for(size_t col = 0; col != M::scolumns(); ++col){
				result.at(row, col) = f(m.at(row, col));
			}
		}
		return result;
	}
	
	// ======================= comparisons =======================
	// ----------------------- equal / operator== -----------------------
	
	template<class Ml, class Mr, TWMATH_ENABLE_IF((is_matrix_v<Ml> || is_matrix_v<Mr>) && !(is_vector_v<Ml> || is_vector_v<Mr>))> 
	constexpr auto equal (const Ml& lhs, const Mr& rhs){
		return for_each(lhs, rhs, twmath_base::equal<value_type_t<Ml>, value_type_t<Mr>>);
	}
	template<class Ml, class Mr, TWMATH_ENABLE_IF((is_matrix_v<Ml> || is_matrix_v<Mr>) && !(is_vector_v<Ml> || is_vector_v<Mr>))> 
	constexpr auto operator== (const Ml& lhs, const Mr& rhs){
		return equal(lhs, rhs);
	}
	
	// ----------------------- not_equal / operator!= -----------------------
	
	template<class Ml, class Mr, TWMATH_ENABLE_IF((is_matrix_v<Ml> || is_matrix_v<Mr>) && !(is_vector_v<Ml> || is_vector_v<Mr>))> 
	constexpr auto not_equal (const Ml& lhs, const Mr& rhs){
		return for_each(lhs, rhs, twmath_base::not_equal<value_type_t<Ml>, value_type_t<Mr>>);
	}
	template<class Ml, class Mr, TWMATH_ENABLE_IF((is_matrix_v<Ml> || is_matrix_v<Mr>) && !(is_vector_v<Ml> || is_vector_v<Mr>))> 
	constexpr auto operator!= (const Ml& lhs, const Mr& rhs){
		return not_equal(lhs, rhs);
	}
	
	// ----------------------- less / operator< -----------------------
	
	template<class Ml, class Mr, TWMATH_ENABLE_IF((is_matrix_v<Ml> || is_matrix_v<Mr>) && !(is_vector_v<Ml> || is_vector_v<Mr>))> 
	constexpr auto less (const Ml& lhs, const Mr& rhs){
		return for_each(lhs, rhs, twmath_base::less<value_type_t<Ml>, value_type_t<Mr>>);
	}
	template<class Ml, class Mr, TWMATH_ENABLE_IF((is_matrix_v<Ml> || is_matrix_v<Mr>) && !(is_vector_v<Ml> || is_vector_v<Mr>))> 
	constexpr auto operator< (const Ml& lhs, const Mr& rhs){
		return less(lhs, rhs);
	}
	
	// ----------------------- less_equal / operator<= -----------------------
	
	template<class Ml, class Mr, TWMATH_ENABLE_IF((is_matrix_v<Ml> || is_matrix_v<Mr>) && !(is_vector_v<Ml> || is_vector_v<Mr>))> 
	constexpr auto less_equal (const Ml& lhs, const Mr& rhs){
		return for_each(lhs, rhs, twmath_base::less_equal<value_type_t<Ml>, value_type_t<Mr>>);}
	
	template<class Ml, class Mr, TWMATH_ENABLE_IF((is_matrix_v<Ml> || is_matrix_v<Mr>) && !(is_vector_v<Ml> || is_vector_v<Mr>))> 
	constexpr auto operator<= (const Ml& lhs, const Mr& rhs){
		return less_equal(lhs, rhs);
	}
	
	// ----------------------- greater / operator> -----------------------
	
	template<class Ml, class Mr, TWMATH_ENABLE_IF((is_matrix_v<Ml> || is_matrix_v<Mr>) && !(is_vector_v<Ml> || is_vector_v<Mr>))> 
	constexpr auto greater (const Ml& lhs, const Mr& rhs){
		return for_each(lhs, rhs, twmath_base::greater<value_type_t<Ml>, value_type_t<Mr>>);
	}
	template<class Ml, class Mr, TWMATH_ENABLE_IF((is_matrix_v<Ml> || is_matrix_v<Mr>) && !(is_vector_v<Ml> || is_vector_v<Mr>))> 
	constexpr auto operator> (const Ml& lhs, const Mr& rhs){
		return greater(lhs, rhs);
	}
	
	// ----------------------- greater_equal / operator>= -----------------------
	
	template<class Ml, class Mr, TWMATH_ENABLE_IF((is_matrix_v<Ml> || is_matrix_v<Mr>) && !(is_vector_v<Ml> || is_vector_v<Mr>))> 
	constexpr auto greater_equal (const Ml& lhs, const Mr& rhs){
		return for_each(lhs, rhs, twmath_base::greater_equal<value_type_t<Ml>, value_type_t<Mr>>);
	}
	template<class Ml, class Mr, TWMATH_ENABLE_IF((is_matrix_v<Ml> || is_matrix_v<Mr>) && !(is_vector_v<Ml> || is_vector_v<Mr>))> 
	constexpr auto operator>= (const Ml& lhs, const Mr& rhs){
		return greater_equal(lhs, rhs);
	}
	
	// ======================= binary arithmetic operators =======================
	// ----------------------- add / operator+ -----------------------
	
	template<class Ml, class Mr, TWMATH_ENABLE_IF((is_matrix_v<Ml> || is_matrix_v<Mr>) && !(is_vector_v<Ml> || is_vector_v<Mr>))> 
	constexpr auto add (const Ml& lhs, const Mr& rhs){
		return for_each(lhs, rhs, twmath_base::add<value_type_t<Ml>, value_type_t<Mr>>);
	}
	template<class Ml, class Mr, TWMATH_ENABLE_IF((is_matrix_v<Ml> || is_matrix_v<Mr>) && !(is_vector_v<Ml> || is_vector_v<Mr>))> 
	constexpr auto operator+ (const Ml& lhs, const Mr& rhs){
		return add(lhs, rhs);
	}
	
	// ---- add, +: for containers of element 1 ----
	template<class M, class V, TWMATH_ENABLE_IF(is_static_matrix_v<M> && is_static_vector_v<V> && M::ssize() == 1 && V::ssize() == 1)> 
	constexpr auto add (const M& lhs, const V& rhs){return lhs.at(0,0) + rhs.at(0);}
	template<class M, class V, TWMATH_ENABLE_IF(is_static_matrix_v<M> && is_static_vector_v<V> && M::ssize() == 1 && V::ssize() == 1)> 
	constexpr auto add (const V& lhs, const M& rhs){return lhs.at(0) + rhs.at(0,0);}
	template<class M, class V, TWMATH_ENABLE_IF(is_static_matrix_v<M> && is_static_vector_v<V> && M::ssize() == 1 && V::ssize() == 1)> 
	constexpr auto operator+ (const M& lhs, const V& rhs){return lhs.at(0,0) + rhs.at(0);}
	template<class M, class V, TWMATH_ENABLE_IF(is_static_matrix_v<M> && is_static_vector_v<V> && M::ssize() == 1 && V::ssize() == 1)> 
	constexpr auto operator+ (const V& lhs, const M& rhs){return lhs.at(0) + rhs.at(0,0);}
	
	// ----------------------- sub / operator- -----------------------
	
	template<class Ml, class Mr, TWMATH_ENABLE_IF((is_matrix_v<Ml> || is_matrix_v<Mr>) && !(is_vector_v<Ml> || is_vector_v<Mr>))> 
	constexpr auto sub (const Ml& lhs, const Mr& rhs){
		return for_each(lhs, rhs, twmath_base::sub<value_type_t<Ml>, value_type_t<Mr>>);
	}
	template<class Ml, class Mr, TWMATH_ENABLE_IF((is_matrix_v<Ml> || is_matrix_v<Mr>) && !(is_vector_v<Ml> || is_vector_v<Mr>))> 
	constexpr auto operator- (const Ml& lhs, const Mr& rhs){
		return sub(lhs, rhs);
	}
	
	// ---- sub, -: for containers of element 1 ----
	template<class M, class V, TWMATH_ENABLE_IF(is_static_matrix_v<M> && is_static_vector_v<V> && M::ssize() == 1 && V::ssize() == 1)> 
	constexpr auto sub (const M& lhs, const V& rhs){return lhs.at(0,0) - rhs.at(0);}
	template<class M, class V, TWMATH_ENABLE_IF(is_static_matrix_v<M> && is_static_vector_v<V> && M::ssize() == 1 && V::ssize() == 1)> 
	constexpr auto sub (const V& lhs, const M& rhs){return lhs.at(0) - rhs.at(0,0);}
	template<class M, class V, TWMATH_ENABLE_IF(is_static_matrix_v<M> && is_static_vector_v<V> && M::ssize() == 1 && V::ssize() == 1)> 
	constexpr auto operator- (const M& lhs, const V& rhs){return lhs.at(0,0) - rhs.at(0);}
	template<class M, class V, TWMATH_ENABLE_IF(is_static_matrix_v<M> && is_static_vector_v<V> && M::ssize() == 1 && V::ssize() == 1)> 
	constexpr auto operator- (const V& lhs, const M& rhs){return lhs.at(0) - rhs.at(0,0);}
	
	// ----------------------- mul / operator* -----------------------
	
	template<class Ml, class Mr, TWMATH_ENABLE_IF((is_matrix_v<Ml> || is_matrix_v<Mr>) && !(is_vector_v<Ml> || is_vector_v<Mr>))> 
	constexpr auto mul (const Ml& lhs, const Mr& rhs){
		return for_each(lhs, rhs, twmath_base::mul<value_type_t<Ml>, value_type_t<Mr>>);
	}
	template<class Ml, class Mr, TWMATH_ENABLE_IF((is_matrix_v<Ml> != is_matrix_v<Mr>) && !is_vector_v<Ml> && !is_vector_v<Mr>)> 
	constexpr auto operator* (const Ml& lhs, const Mr& rhs){return mul(lhs, rhs);}
	
	// ---- mul, *: for containers of element 1 ----
	template<class M, class V, TWMATH_ENABLE_IF(is_static_matrix_v<M> && is_static_vector_v<V> && M::ssize() == 1 && V::ssize() == 1)> 
	constexpr auto mul (const M& lhs, const V& rhs){return lhs.at(0,0) * rhs.at(0);}
	template<class M, class V, TWMATH_ENABLE_IF(is_static_matrix_v<M> && is_static_vector_v<V> && M::ssize() == 1 && V::ssize() == 1)> 
	constexpr auto mul (const V& lhs, const M& rhs){return lhs.at(0) * rhs.at(0,0);}
	template<class M, class V, TWMATH_ENABLE_IF(is_static_matrix_v<M> && is_static_vector_v<V> && M::ssize() == 1 && V::ssize() == 1)> 
	constexpr auto operator* (const M& lhs, const V& rhs){return lhs.at(0,0) * rhs.at(0);}
	template<class M, class V, TWMATH_ENABLE_IF(is_static_matrix_v<M> && is_static_vector_v<V> && M::ssize() == 1 && V::ssize() == 1)> 
	constexpr auto operator* (const V& lhs, const M& rhs){return lhs.at(0) * rhs.at(0,0);}
	
	// ----------------------- div / operator/ -----------------------
	
	template<class Ml, class Mr, TWMATH_ENABLE_IF((is_matrix_v<Ml> || is_matrix_v<Mr>) && !(is_vector_v<Ml> || is_vector_v<Mr>))> 
	constexpr auto div (const Ml& lhs, const Mr& rhs){
		return for_each(lhs, rhs, twmath_base::div<value_type_t<Ml>, value_type_t<Mr>>);
	}
	template<class Ml, class Mr, TWMATH_ENABLE_IF((is_matrix_v<Ml> != is_matrix_v<Mr>) && !is_vector_v<Ml> && !is_vector_v<Mr>)> 
	constexpr auto operator/ (const Ml& lhs, const Mr& rhs){return div(lhs, rhs);}
	
	// ---- div, /: for containers of element 1 ----
	template<class M, class V, TWMATH_ENABLE_IF(is_static_matrix_v<M> && is_static_vector_v<V> && M::ssize() == 1 && V::ssize() == 1)> 
	constexpr auto div (const M& lhs, const V& rhs){return lhs.at(0,0) / rhs.at(0);}
	template<class M, class V, TWMATH_ENABLE_IF(is_static_matrix_v<M> && is_static_vector_v<V> && M::ssize() == 1 && V::ssize() == 1)> 
	constexpr auto div (const V& lhs, const M& rhs){return lhs.at(0) / rhs.at(0,0);}
	template<class M, class V, TWMATH_ENABLE_IF(is_static_matrix_v<M> && is_static_vector_v<V> && M::ssize() == 1 && V::ssize() == 1)> 
	constexpr auto operator/ (const M& lhs, const V& rhs){return lhs.at(0,0) / rhs.at(0);}
	template<class M, class V, TWMATH_ENABLE_IF(is_static_matrix_v<M> && is_static_vector_v<V> && M::ssize() == 1 && V::ssize() == 1)> 
	constexpr auto operator/ (const V& lhs, const M& rhs){return lhs.at(0) / rhs.at(0,0);}
	
	// ----------------------- negate -----------------------
	
	template<class V, TWMATH_ENABLE_IF(is_matrix_v<V>)> 
	constexpr auto negate (const V& a){return for_each(a, twmath_base::negate<value_type_t<V>>);}
	
	template<class V, TWMATH_ENABLE_IF(is_matrix_v<V>)>
	constexpr auto operator- (const V& a){return for_each(a, twmath_base::negate<value_type_t<V>>);}
	
	// ----------------------- matmul / operator * -----------------------

	template<class Ml, class Mr, TWMATH_ENABLE_IF(is_static_matrix_v<Ml> && is_static_matrix_v<Mr> && Ml::scolumns() == Mr::srows())>
	constexpr Matrix<result_type2((twmath_base::mul<value_type_t<Ml>, value_type_t<Mr>>), value_type_t<Ml>, value_type_t<Mr>), Ml::srows(), Mr::scolumns()>
	matmul (const Ml& ml, const Mr& mr){
		Matrix<result_type2((twmath_base::mul<value_type_t<Ml>, value_type_t<Mr>>), value_type_t<Ml>, value_type_t<Mr>), Ml::srows(), Mr::scolumns()> result;
		for(size_t row=0; row != result.rows(); ++row){
			for(size_t col=0; col != result.columns(); ++col){
				result.at(row, col) = ml.row(row) * mr.column(col);
			}
		}
		return result;
	}
	
	template<class Ml, class Mr, TWMATH_ENABLE_IF(is_matrix_v<Ml> && is_matrix_v<Mr>)>
	constexpr auto operator* (const Ml& ml, const Mr& mr){return matmul(ml, mr);}
	
	// ----------------------- matrix vector mul / operator * -----------------------
	
	template<class Ml, class Vr, TWMATH_ENABLE_IF(is_static_matrix_v<Ml> && is_static_vector_v<Vr> && Ml::scolumns() == Vr::ssize())>
	constexpr Vector<result_type2((twmath_base::mul<value_type_t<Ml>, value_type_t<Vr>>), value_type_t<Ml>, value_type_t<Vr>), Ml::srows()> 
	matmul (const Ml& ml, const Vr& vr){
		Vector<result_type2((twmath_base::mul<value_type_t<Ml>, value_type_t<Vr>>), value_type_t<Ml>, value_type_t<Vr>), Ml::srows()> result;
		for(size_t row = 0; row < Ml::srows(); ++row){
			result.at(row) = ml.row(row) * vr;
		}
		return result;
	}
	template<class Ml, class Vr, TWMATH_ENABLE_IF(is_matrix_v<Ml> && is_vector_v<Vr>)>
	constexpr auto operator* (const Ml& ml, const Vr& vr){return matmul(ml, vr);}
	
	// ----------------------- vector matrix mul / operator * -----------------------
	
	template<class Vl, class Mr, TWMATH_ENABLE_IF(is_static_vector_v<Vl> && is_static_matrix_v<Mr> && Vl::ssize() == Mr::srows())>
	constexpr Vector<result_type2((twmath_base::mul<value_type_t<Vl>, value_type_t<Mr>>), value_type_t<Vl>, value_type_t<Mr>), Mr::scolumns()> 
	matmul (const Vl& vl, const Mr& mr){
		Vector<result_type2((twmath_base::mul<value_type_t<Vl>, value_type_t<Mr>>), value_type_t<Vl>, value_type_t<Mr>), Mr::scolumns()>  result;
		for(size_t col = 0; col < Mr::scolumns(); ++col){
			result.at(col) = vl * mr.column(col);
		}		
		return result;
	}
	template<class Vl, class Mr, TWMATH_ENABLE_IF(is_vector_v<Vl> && is_matrix_v<Mr>)>
	constexpr auto operator* (const Vl& vl, const Mr& mr){return matmul(vl, mr);}
	
	// TODO: matrix inverse
	// TODO: matrix adjunctate
	// TODO: matrix division
	
	
	// ------------------------- Matrix Exponentiation ----------------------------------------
	
	// approximates the exponential of the matrix with the tailor expansion:
	// e^A = I + A + A^2 / 2 + ... + A^n / n!
	// the minimum number of n is 2. If n is set lower than 2, then 2 increments will be calculated regardless
	template<class M, TWMATH_ENABLE_IF(is_static_matrix_v<M> && M::srows() == M::scolumns())>
	constexpr M exp_taylor(const M& A, size_t n = 4){
		M A_pow = A * A;
		size_t factorial = 2;
		M result = unity_like(A) + A + A_pow / factorial;
		for(size_t i = 3; i <= n; ++i){
			A_pow = A_pow * A;
			factorial *= i;
			result = result + A_pow / factorial;
		}
		return result;
	}
	
	// improves the accuracy of the taylor expansion using e^x = (e^(x/s))^s
	// accuracy is increased because the taylor expansion is more precise for smaller values of x.
	// e^x ~ taylor_n(x/(2^s))^(2^s)
	template<class M, TWMATH_ENABLE_IF(is_matrix_v<M>)>
	constexpr M exp_taylor_square_scale(const M& x, size_t taylor_order = 4, size_t scaling = 4){
		M t = exp_taylor(x/(1 << scaling), taylor_order);
		for(size_t i = 0; i < scaling; ++i) t = t * t;
		return t;
	}
	
	// default matrix exponentiation
	template<class M, TWMATH_ENABLE_IF(is_matrix_v<M>)>
	constexpr M mexp(const M& x){return exp_taylor_square_scale(x);}
	
	// ------------------------------------------- generators ---------------------------------------
	
	template<class T, size_t _rows, size_t _columns>
	constexpr Matrix<T, _rows, _columns> zeros_like(const Matrix<T, _rows, _columns>& A){
		return zeros<T, _rows, _columns>();
	}
	
	template<class T, size_t _rows, size_t _columns=_rows>
	constexpr Matrix<T, _rows, _columns> unity(){
		Matrix<T, _rows, _columns> result;
		for(size_t row = 0; row < _rows; ++row){
			for(size_t col = 0; col < _columns; ++col){
				if(row == col){
					result.at(row, col) = T(1);
				}else{
					result.at(row, col) = T(0);
				}
			}
		}
		return result;
	}
	
	template<class T, size_t _rows, size_t _columns>
	constexpr Matrix<T, _rows, _columns> unity_like(const Matrix<T, _rows, _columns>& A){
		return unity<T, _rows, _columns>();
	}
	
	// ------------------------------------------- printing ---------------------------------------
	
	template<class Stream, class M, TWMATH_ENABLE_IF(is_matrix_v<M>)>
	Stream& print_pretty(Stream& stream, const M& m, const char* name = nullptr, const char* indentation=""){
		size_t name_len = 0; 
		if(name!=nullptr) name_len = std::strlen(name);  
		else name_len = 0;
		
		size_t max_elem_size = 0;
		if(m.rows() > 1) for(size_t row=0; row<m.rows(); ++row) for(size_t column=0; column<m.columns(); ++column){
			auto& elem = m.at(row, column);
			std::stringstream s;
			s << elem;
			max_elem_size = (s.str().size() > max_elem_size) ? s.str().size() : max_elem_size;
		}
		
		for(size_t row = 0; row < m.rows(); ++row){
			char open_bracket; 
			char closed_bracket;
			
			if(m.rows() == 1){
				open_bracket = '(';
				closed_bracket = ')';
			}else if(row == 0){
				open_bracket = '/';
				closed_bracket = '\\';
			}else if (row == m.rows()-1){
				open_bracket = '\\';
				closed_bracket = '/';
			}else{
				open_bracket = '|';
				closed_bracket = '|';
			}
			
			// print Matrix
			stream << indentation;
			if (name!=nullptr){
				if(row == m.rows()/2) stream << name << " = ";	
				else for(size_t i = 0; i < name_len + 3; ++i) stream << ' ';
			}
			stream << open_bracket;
			for(size_t column = 0; column < m.columns(); ++column){
				std::stringstream s;
				if(column!=0) stream << ' ';
				s << m.at(row, column);
				if(m.rows() > 1) for(size_t i = 0; i < max_elem_size - s.str().size(); ++i) stream << ' ';
				stream << s.str();
			}
			stream << closed_bracket << '\n';
		}
		
		return stream;
	}

	template<class Stream, class M, TWMATH_ENABLE_IF(is_matrix_v<M>)>
	Stream& print_csv(Stream& stream, const M& m, const char* value_separator=", ", const char* line_separator="\n") {
		for (size_t row = 0; row < m.rows(); row++) {
			for (size_t column = 0; column < m.columns(); column++) {
				if(column == 0) stream << m.at(row, column);
				else stream << value_separator << m.at(row, column);
			}
			stream << line_separator;
		}
		return stream;
	}
	
	template<class Stream, class M, TWMATH_ENABLE_IF(is_matrix_v<M>)>
	Stream& operator<< (Stream& stream, const M& m){
		return print_csv(stream, m);
	}
	
}