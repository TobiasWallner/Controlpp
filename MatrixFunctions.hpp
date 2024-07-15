#pragma once

#include "Matrix.hpp"
#include "MatrixRef.hpp"


namespace twmath{
	
	
	template<class T, size_t rows, size_t columns>
	constexpr MatrixRef<T, rows, columns> transpose(MatrixRef<T, rows, columns>& M){
		return MatrixRef<T, rows, columns>(M.data(), M.row_increment(), M.column_increment(), !M.is_transposed());
	}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns, class F>
	constexpr Matrix<F, _rows, _columns> 
	for_each (const MatrixRef<Tl, _rows, _columns>& lhs, const MatrixRef<Tr, _rows, _columns>& rhs, F (*f)(const Tl&, const Tr&)){
		Matrix<F, _rows, _columns> result;
		auto result_itr = result.begin();
		auto lhs_itr = lhs.begin();
		auto rhs_itr = rhs.begin();
		for(; result_itr != result.end(); ++result_itr, (void)++lhs_itr, (void)++rhs_itr) *result_itr = f(*lhs_itr, *rhs_itr);
		return result;
	}
	template<class Tl, class Tr, size_t _rows, size_t _columns, class F>
	constexpr Matrix<F, _rows, _columns> 
	for_each (const Tl& lhs, const MatrixRef<Tr, _rows, _columns>& rhs, F (*f)(const Tl&, const Tr&)){
		Matrix<F, _rows, _columns> result;
		auto result_itr = result.begin();
		auto rhs_itr = rhs.begin();
		for(; result_itr != result.end(); ++result_itr, (void)++rhs_itr) *result_itr = f(lhs, *rhs_itr);
		return result;
	}
	template<class Tl, class Tr, size_t _rows, size_t _columns, class F>
	constexpr Matrix<F, _rows, _columns> 
	for_each (const MatrixRef<Tl, _rows, _columns>& lhs, const Tr& rhs, F (*f)(const Tl&, const Tr&)){
		Matrix<F, _rows, _columns> result;
		auto result_itr = result.begin();
		auto lhs_itr = lhs.begin();
		for(; result_itr != result.end(); ++result_itr, (void)++lhs_itr) *result_itr = f(*lhs_itr, rhs);
		return result;
	}
	template<class T, size_t _rows, size_t _columns, class F>
	constexpr Matrix<F, _rows, _columns> 
	for_each (const MatrixRef<T, _rows, _columns>& v, F (*f)(const T&)){
		Matrix<F, _rows, _columns> result;
		auto result_itr = result.begin();
		auto v_itr = v.begin();
		for(; result_itr != result.end(); ++result_itr, (void)++v_itr) *result_itr = f(*v_itr);
		return result;
	}
	
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::add<Tl, Tr>), Tl, Tr), _rows, _columns>  
	add (const MatrixRef<Tl, _rows, _columns>& lhs, const MatrixRef<Tr, _rows, _columns>& rhs){return for_each(lhs, rhs, twmath_base::add<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::add<Tl, Tr>), Tl, Tr), _rows, _columns>  
	add (const MatrixRef<Tl, _rows, _columns>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::add<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::add<Tl, Tr>), Tl, Tr), _rows, _columns>  
	add (const Tl& lhs, const MatrixRef<Tr, _rows, _columns>& rhs){return for_each(lhs, rhs, twmath_base::add<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::add<Tl, Tr>), Tl, Tr), _rows, _columns>  
	operator+ (const MatrixRef<Tl, _rows, _columns>& lhs, const MatrixRef<Tr, _rows, _columns>& rhs){return for_each(lhs, rhs, twmath_base::add<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::add<Tl, Tr>), Tl, Tr), _rows, _columns>  
	operator+ (const MatrixRef<Tl, _rows, _columns>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::add<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::add<Tl, Tr>), Tl, Tr), _rows, _columns>  
	operator+ (const Tl& lhs, const MatrixRef<Tr, _rows, _columns>& rhs){return for_each(lhs, rhs, twmath_base::add<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows>
	constexpr Vector<result_type2((twmath_base::add<Tl, Tr>), Tl, Tr), _rows>  
	operator+ (const MatrixRef<Tl, _rows, 1>& lhs, const Vector<Tr, _rows>& rhs){
		Vector<result_type2((twmath_base::add<Tl, Tr>), Tl, Tr), _rows> result; 
		const auto result_end = result.end();
		auto result_itr = result.begin();
		auto lhs_itr = lhs.begin();
		auto rhs_itr = rhs.begin();
		while(result_itr != result_end){
			*(result_itr++) = *(lhs_itr++) + *(rhs_itr++);
		}
		return result;
	}
	
	template<class Tl, class Tr, size_t _columns>
	constexpr Vector<result_type2((twmath_base::add<Tl, Tr>), Tl, Tr), _columns>  
	operator+ (const MatrixRef<Tl, 1, _columns>& lhs, const Vector<Tr, _columns>& rhs){
		Vector<result_type2((twmath_base::add<Tl, Tr>), Tl, Tr), _columns> result; 
		const auto result_end = result.end();
		auto result_itr = result.begin();
		auto lhs_itr = lhs.begin();
		auto rhs_itr = rhs.begin();
		while(result_itr != result_end){
			*(result_itr++) = *(lhs_itr++) + *(rhs_itr++);
		}
		return result;
	}
	
	template<class Tl, class Tr>
	constexpr result_type2((twmath_base::add<Tl, Tr>), Tl, Tr)
	operator+ (const MatrixRef<Tl, 1, 1>& lhs, const Vector<Tr, 1>& rhs){
		return lhs.at(0, 0) + rhs.at(0);
	}
	
	template<class Tl, class Tr, size_t _rows>
	constexpr Vector<result_type2((twmath_base::add<Tl, Tr>), Tl, Tr), _rows>  
	operator+ (const Vector<Tr, _rows>& lhs, const MatrixRef<Tl, _rows, 1>& rhs){
		Vector<result_type2((twmath_base::add<Tl, Tr>), Tl, Tr), _rows> result; 
		const auto result_end = result.end();
		auto result_itr = result.begin();
		auto lhs_itr = lhs.begin();
		auto rhs_itr = rhs.begin();
		while(result_itr != result_end){
			*(result_itr++) = *(lhs_itr++) + *(rhs_itr++);
		}
		return result;
	}
	
	template<class Tl, class Tr, size_t _columns>
	constexpr Vector<result_type2((twmath_base::add<Tl, Tr>), Tl, Tr), _columns>  
	operator+ (const Vector<Tr, _columns>& lhs, const MatrixRef<Tl, 1, _columns>& rhs){
		Vector<result_type2((twmath_base::add<Tl, Tr>), Tl, Tr), _columns> result; 
		const auto result_end = result.end();
		auto result_itr = result.begin();
		auto lhs_itr = lhs.begin();
		auto rhs_itr = rhs.begin();
		while(result_itr != result_end){
			*(result_itr++) = *(lhs_itr++) + *(rhs_itr++);
		}
		return result;
	}
	
	template<class Tl, class Tr>
	constexpr result_type2((twmath_base::add<Tl, Tr>), Tl, Tr)
	operator+ (const Vector<Tr, 1>& lhs, const MatrixRef<Tl, 1, 1>& rhs){
		return lhs.at(0) + rhs.at(0, 0);
	}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::sub<Tl, Tr>), Tl, Tr), _rows, _columns>  
	sub (const MatrixRef<Tl, _rows, _columns>& lhs, const MatrixRef<Tr, _rows, _columns>& rhs){return for_each(lhs, rhs, twmath_base::sub<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::sub<Tl, Tr>), Tl, Tr), _rows, _columns>  
	sub (const MatrixRef<Tl, _rows, _columns>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::sub<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::sub<Tl, Tr>), Tl, Tr), _rows, _columns>  
	sub (const Tl& lhs, const MatrixRef<Tr, _rows, _columns>& rhs){return for_each(lhs, rhs, twmath_base::sub<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::sub<Tl, Tr>), Tl, Tr), _rows, _columns>  
	operator- (const MatrixRef<Tl, _rows, _columns>& lhs, const MatrixRef<Tr, _rows, _columns>& rhs){return for_each(lhs, rhs, twmath_base::sub<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::sub<Tl, Tr>), Tl, Tr), _rows, _columns>  
	operator- (const MatrixRef<Tl, _rows, _columns>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::sub<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::sub<Tl, Tr>), Tl, Tr), _rows, _columns>  
	operator- (const Tl& lhs, const MatrixRef<Tr, _rows, _columns>& rhs){return for_each(lhs, rhs, twmath_base::sub<Tl, Tr>);}
	
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::mul<Tl, Tr>), Tl, Tr), _rows, _columns>  
	mul (const MatrixRef<Tl, _rows, _columns>& lhs, const MatrixRef<Tr, _rows, _columns>& rhs){return for_each(lhs, rhs, twmath_base::mul<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::mul<Tl, Tr>), Tl, Tr), _rows, _columns>  
	mul (const MatrixRef<Tl, _rows, _columns>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::mul<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::mul<Tl, Tr>), Tl, Tr), _rows, _columns> 
	mul (const Tl& lhs, const MatrixRef<Tr, _rows, _columns>& rhs){return for_each(lhs, rhs, twmath_base::mul<Tl, Tr>);}

	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::mul<Tl, Tr>), Tl, Tr), _rows, _columns> 
	operator* (const MatrixRef<Tl, _rows, _columns>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::mul<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::mul<Tl, Tr>), Tl, Tr), _rows, _columns> 
	operator* (const Tl& lhs, const MatrixRef<Tr, _rows, _columns>& rhs){return for_each(lhs, rhs, twmath_base::mul<Tl, Tr>);}


	template<class Tl, class Tr, size_t lcols_rrows, size_t lrows, size_t rcols>
	constexpr Matrix<result_type2((twmath_base::mul<Tl, Tr>), Tl, Tr), rcols, lrows> 
	matmul (const MatrixRef<Tl, lrows, lcols_rrows>& lhs, const MatrixRef<Tr, lcols_rrows, rcols>& rhs){
		Matrix<result_type2((twmath_base::mul<Tl, Tr>), Tl, Tr), rcols, lrows> result;
		
		for(size_t result_row = 0; result_row < lrows; ++result_row){
			for(size_t result_col = 0; result_col < rcols; ++result_col){
				
				result.at(result_row, result_col) = lhs.at(result_row, 0) * rhs.at(0, result_col);
				for(size_t lcol_rrow = 1; lcol_rrow < lcols_rrows; ++lcol_rrow){
					result.at(result_row, result_col) += lhs.at(result_row, lcol_rrow) * rhs.at(lcol_rrow, result_col);
				}
				
			}
			
		}
		return result;
	}
	template<class Tl, class Tr, size_t lcols_rrows, size_t lrows, size_t rcols>
	constexpr Matrix<result_type2((twmath_base::mul<Tl, Tr>), Tl, Tr), rcols, lrows> 
	operator* (const MatrixRef<Tl, lrows, lcols_rrows>& lhs, const MatrixRef<Tr, lcols_rrows, rcols>& rhs){return matmul(lhs, rhs);}
	
	template<class Tl, class Tr, size_t rows, size_t cols>
	constexpr Vector<result_type2((twmath_base::mul<Tl, Tr>), Tl, Tr), rows> mul (const MatrixRef<Tl, rows, cols>& l, const Vector<Tr, cols>& r){
		Vector<result_type2((twmath_base::mul<Tl, Tr>), Tl, Tr), rows> result;
		for(size_t row = 0; row < rows; ++row){
			result.at(row) = l.at(row, 0) * r.at(0);
			for(size_t col = 1; col < cols; ++col){
				result.at(row) += l.at(row, col) * r.at(col);
			}	
		}
		return result;
	}
	template<class Tl, class Tr, size_t rows, size_t cols>
	constexpr Vector<result_type2((twmath_base::mul<Tl, Tr>), Tl, Tr), rows> operator* (const MatrixRef<Tl, rows, cols>& l, const Vector<Tr, cols>& r){
		return mul<Tl, Tr, rows, cols>(l, r);
	}
	
	template<class Tl, class Tr, size_t rows, size_t cols>
	constexpr Vector<result_type2((twmath_base::mul<Tl, Tr>), Tl, Tr), cols> mul (const Vector<Tr, rows>& l, const MatrixRef<Tl, rows, cols>& r){
		Vector<result_type2((twmath_base::mul<Tl, Tr>), Tl, Tr), cols> result;
		for(size_t col = 0; col < cols; ++col){
			result.at(0) += l.at(col) * r.at(0, col);
			for(size_t row = 1; row < rows; ++row){
				result.at(row) += r.at(col) * l.at(row, col);
			}
		}		
		return result;
	}
	template<class Tl, class Tr, size_t rows, size_t cols>
	constexpr Vector<result_type2((twmath_base::mul<Tl, Tr>), Tl, Tr), rows> operator* (const Vector<Tr, rows>& l, const MatrixRef<Tl, rows, cols>& r){
		return mul<Tl, Tr, rows, cols>(l, r);
	}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::div<Tl, Tr>), Tl, Tr), _rows, _columns> 
	div (const MatrixRef<Tl, _rows, _columns>& lhs, const MatrixRef<Tr, _rows, _columns>& rhs){return for_each(lhs, rhs, twmath_base::div<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::div<Tl, Tr>), Tl, Tr), _rows, _columns> 
	div (const MatrixRef<Tl, _rows, _columns>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::div<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::div<Tl, Tr>), Tl, Tr), _rows, _columns>  
	div (const Tl& lhs, const MatrixRef<Tr, _rows, _columns>& rhs){return for_each(lhs, rhs, twmath_base::div<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::div<Tl, Tr>), Tl, Tr), _rows, _columns>   
	operator/ (const MatrixRef<Tl, _rows, _columns>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::div<Tl, Tr>);}
	
	template<class T, size_t _rows, size_t _columns>
	constexpr Matrix<T, _rows, _columns> zeros(){
		Matrix<T, _rows, _columns> result;
		for(const auto& elem : result) elem = T(0);
		return result;
	}
	
	template<class T, size_t _rows, size_t _columns>
	constexpr Matrix<T, _rows, _columns> zeros_like(const MatrixRef<T, _rows, _columns>& A){
		return zeros<T, _rows, _columns>();
	}
	
	template<class T, size_t _rows, size_t _columns>
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
	constexpr Matrix<T, _rows, _columns> unity_like(const MatrixRef<T, _rows, _columns>& A){
		return unity<T, _rows, _columns>();
	}
	
	template<class Stream, class T, size_t _rows, size_t _columns>
	Stream& print_pretty(Stream& stream, const MatrixRef<T, _rows, _columns>& M, const char* name = nullptr, const char* indentation=""){
		size_t name_len = 0; 
		if(name!=nullptr) name_len = std::strlen(name);  
		else name_len = 0;
		
		size_t max_elem_size = 0;
		if(M.rows() > 1) for(const auto& elem : M){
			std::stringstream s;
			s << elem;
			max_elem_size = (s.str().size() > max_elem_size) ? s.str().size() : max_elem_size;
		}
		
		for(size_t row = 0; row < M.rows(); ++row){
			char open_bracket; 
			char closed_bracket;
			
			if(M.rows() == 1){
				open_bracket = '(';
				closed_bracket = ')';
			}else if(row == 0){
				open_bracket = '/';
				closed_bracket = '\\';
			}else if (row == M.rows()-1){
				open_bracket = '\\';
				closed_bracket = '/';
			}else{
				open_bracket = '|';
				closed_bracket = '|';
			}
			
			// print Matrix
			stream << indentation;
			if (name!=nullptr){
				if(row == M.rows()/2) stream << name << " = ";	
				else for(size_t i = 0; i < name_len + 3; ++i) stream << ' ';
			}
			stream << open_bracket;
			for(size_t column = 0; column < M.columns(); ++column){
				std::stringstream s;
				if(column!=0) stream << ' ';
				s << M.at(row, column);
				if(M.rows() > 1) for(size_t i = 0; i < max_elem_size - s.str().size(); ++i) stream << ' ';
				stream << s.str();
			}
			stream << closed_bracket << '\n';
		}
		
		return stream;
	}

	template<class Stream, class T, size_t _rows, size_t _columns>
	Stream& print_csv(Stream& stream, const MatrixRef<T, _rows, _columns>& m, const char* value_separator=", ", const char* line_separator="\n") {
		for (size_t row = 0; row < m.rows(); row++) {
			for (size_t column = 0; column < m.columns(); column++) {
				if(column == 0) stream << m.at(row, column);
				else stream << value_separator << m.at(row, column);
			}
			stream << line_separator;
		}
		return stream;
	}
	
	template<class Stream, class T, size_t _rows, size_t _columns>
	Stream& operator<< (Stream& stream, const MatrixRef<T, _rows, _columns>& m){
		return print_csv(stream, m);
	}
	
	// approximates the exponential of the matrix with the tailor expansion:
	// e^A = I + A + A^2 / 2 + ... + A^n / n!
	// the minimum number of n is 2. If n is set lower than 2, then 2 increments will be calculated regardless
	template<class T, size_t _rows, size_t _columns>
	constexpr Matrix<T, _rows, _columns> exp_taylor(const MatrixRef<T, _rows, _columns>& A, size_t n = 4){
		Matrix<T, _rows, _columns> A_pow = A * A;
		T factorial = T(2);
		Matrix<T, _rows, _columns> result = unity_like(A) + A + A_pow / factorial;
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
	template<class T, size_t _rows, size_t _columns>
	constexpr Matrix<T, _rows, _columns> exp_taylor_square_scale(const MatrixRef<T, _rows, _columns>& x, size_t taylor_order = 4, size_t scaling = 4){
		Matrix<T, _rows, _columns> t = exp_taylor(x/(1 << scaling), taylor_order);
		for(size_t i = 0; i < scaling; ++i) t = t * t;
		return t;
	}
	
	// default matrix exponentiation
	template<class T, size_t _rows, size_t _columns>
	constexpr Matrix<T, _rows, _columns> exp(const MatrixRef<T, _rows, _columns>& x){
		return exp_taylor_square_scale(x);
	}
	
}