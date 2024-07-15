#pragma once
#pragma once

#include <cstddef>
#include <cstdint>
#include <cmath>
#include <utility>
#include <cstring>
#include <type_traits>


#include "definitions.hpp"
#include "BaseTypeArithmetic.hpp"
#include "Vector.hpp"
#include "MatrixRef.hpp"

namespace twmath{
	
	
	
	template<class T, size_t _rows, size_t _columns>
	class Matrix{
	private:
		T _values[_columns * _rows];
		
	public:
	
		constexpr operator MatrixRef<T, _rows, _columns> () {return reference(*this);}
	
		Matrix() = default;
		
		explicit Matrix(const T& v){for(T& elem : *this) elem = v;}
		
		explicit constexpr Matrix(const T (&array)[_rows][_columns]){
			for(size_t row = 0; row < _rows; ++row){
				for(size_t column = 0; column < _columns; ++column){
					this->at(row, column) = array[row][column];
				}
			}
		}
		
		constexpr Matrix& operator=(const T (&array)[_rows][_columns]){
			for(size_t row = 0; row < _rows; ++row){
				for(size_t column = 0; column < _columns; ++column){
					this->at(row, column) = array[row][column];
				}
			}
		}
		
		constexpr Matrix(const Matrix&) = default;
		constexpr Matrix& operator=(const Matrix&) = default;
		
		constexpr size_t columns()const{return _columns;}
		constexpr size_t rows()const{return _rows;}
		
		constexpr size_t row_increment() const {return this->columns();}
		constexpr size_t column_increment() const {return 1;}
		constexpr size_t is_transposed() const {return false;}
		
		constexpr size_t size()const{return this->columns() * this->rows();}
		
		constexpr T* data() {return &(this->_values[0]);}
		constexpr const T* data() const {return &(this->_values[0]);}
		
		constexpr T& at(ptrdiff_t row, ptrdiff_t column){
			row = (row >= 0) ? row : row + this->rows();
			column = (column >= 0) ? column : column + this->columns();
			twmath_assert(column < this->columns(), "Column '" << column << "' is out of bound of '" << this->columns() << "'.");
			twmath_assert(row < this->rows(), "Row '" << row << "' is out of bound of '" << this->rows() << "'.");
			return _values[column + row * _columns];
		} 
		
		constexpr const T& at(ptrdiff_t row, ptrdiff_t column)const{
			row = (row >= 0) ? row : row + this->size();
			column = (column >= 0) ? column : column + this->size();
			twmath_assert(column < this->columns(), "Column '" << column << "' is out of bound of '" << this->columns() << "'.");
			twmath_assert(row < this->rows(), "Row '" << row << "' is out of bound of '" << this->rows() << "'.");
			return _values[column + row * _columns];
		}
	
		constexpr T* begin(){return this->_values;}
		constexpr const T* begin()const{return this->_values;}
		constexpr const T* cbegin()const{return this->_values;}
		
		constexpr T* end(){return this->begin() + this->size();}
		constexpr const T* end()const{return this->begin() + this->size();}
		constexpr const T* cend()const{return this->begin() + this->size();}

		template<size_t sub_rows, size_t sub_columns>
		constexpr MatrixRef<T, sub_rows, sub_columns> sub_matrix(size_t row_offset=0, size_t column_offset=0, size_t row_stride=1, size_t column_stride=1){
			twmath_assert(row_stride!=0, "Matrix Error: row_stride cannot be zero.");
			twmath_assert(column_stride!=0, "Matrix Error: column_stride cannot be zero.");
			T* ptr = &this->at(row_offset, column_offset);
			const size_t sub_column_increment = column_stride;
			const size_t sub_row_increment = row_stride * this->columns();
			MatrixRef<T, sub_rows, sub_columns> result(ptr, sub_row_increment, sub_column_increment);
			return result;
		}
		
		template<size_t rows=_rows>
		constexpr MatrixRef<T, rows, 1> sub_column(size_t column_number, size_t offset=0, size_t stride=1) {
			MatrixRef<T, rows, 1> result = this->sub_matrix<rows, 1>(offset, column_number, 1, stride);
			return result;
		}
		
		template<size_t columns=_columns>
		constexpr MatrixRef<T, 1, columns> sub_row(size_t row_number, size_t offset=0, size_t stride=1) {
			MatrixRef<T, 1, columns> result = this->sub_matrix<1, columns>(row_number, offset, stride, 1);
			return result;
		}
		
		constexpr MatrixRef<T, (_rows < _columns) ? _rows : _columns, 1> diag() {
			MatrixRef<T, (_rows < _columns) ? _rows : _columns, 1> result(this->data(), this->columns() + 1, 1);
			return result;
		}
		
		
		template<class Ty, size_t sub_rows, size_t sub_columns>
		constexpr Matrix& operator=(const Matrix<Ty, sub_rows, sub_columns>& sub_matrix){
			return this->assign(sub_matrix);
		}

		template<size_t sub_size>
		constexpr Vector<T, sub_size> sub_column(ptrdiff_t column_offset, ptrdiff_t row_offset)const{
			column_offset = (column_offset >= 0) ? column_offset : column_offset + this->columns();
			row_offset = (row_offset >= 0) ? row_offset : row_offset + this->rows();
			
			Vector<T, sub_size> result(0);
			auto result_itr = result.begin();
			size_t this_row = row_offset;
			for(; (result_itr != result.end()) && (this_row < this->rows()); (++result_itr), ((void)++this_row)){
				*result_itr = this->at(this_row, column_offset);
			}
			return result;
		}
		
		constexpr Vector<T, _rows> column(ptrdiff_t column_offset)const{
			column_offset = (column_offset >= 0) ? column_offset : column_offset + this->columns();
			Vector<T, _rows> result;
			for(size_t i = 0; i < this->rows(); ++i){result[i] = this->at(i, column_offset);}
			return result;
		}
		
		template<size_t sub_size>
		constexpr Matrix& assign_sub_column(const Vector<T, sub_size>& column, ptrdiff_t row_offset, ptrdiff_t column_offset){
			column_offset = (column_offset >= 0) ? column_offset : column_offset + this->columns();
			row_offset = (row_offset >= 0) ? row_offset : row_offset + this->rows();
			
			size_t result_i = 0;
			size_t this_row = row_offset;
			for(; result_i < sub_size && this_row < _rows; ++result_i, (void)++this_row){
				this->at(this_row, column_offset) = column.at(result_i);
			}
			return *this;
		}
		
		constexpr Matrix& assign_column(const Vector<T, _rows>& column, ptrdiff_t column_offset){
			column_offset = (column_offset >= 0) ? column_offset : column_offset + this->columns();
			
			return this->assign_sub_column<_rows>(column, 0, column_offset);
		}
		
		template<class Iterator>
		constexpr Matrix& assign_column(Iterator first, Iterator last, ptrdiff_t row_offset, ptrdiff_t column_offset){
			row_offset = (row_offset >= 0) ? row_offset : row_offset + this->rows();
			column_offset = (column_offset >= 0) ? column_offset : column_offset + this->columns();
			
			for(size_t row = row_offset; first != last && row < _rows; ++first, (void)++row) this->at(row, column_offset) = *first;
			return *this;
		}
		
		template<class Iterator>
		constexpr Matrix& assign_column(Iterator first, Iterator last, ptrdiff_t column_offset){
			return assign_column(first, last, column_offset, 0);
		}
		
		template<class Iterator>
		constexpr Matrix& assign_column(Iterator first, Iterator last){
			return assign_column(first, last, 0, 0);
		}
		
		template<class Ta, size_t Na>
		constexpr Matrix& assign_column(const Vector<Ta, Na>& v, ptrdiff_t row_offset, ptrdiff_t column_offset){
			return assign_column(v.begin(), v.end(), column_offset, row_offset);
		}
		
		template<class Iterator>
		constexpr Matrix& assign_column(Iterator first, ptrdiff_t row_offset, ptrdiff_t column_offset){
			row_offset = (row_offset >= 0) ? row_offset : row_offset + this->rows();
			column_offset = (column_offset >= 0) ? column_offset : column_offset + this->columns();
			for(size_t row = row_offset; row < _rows; ++first, (void)++row) this->at(row, column_offset) = *first;
			return *this;
		}
		
		template<class Iterator>
		constexpr Matrix& assign_column(Iterator first, ptrdiff_t column_offset){
			return assign_column(first, column_offset, 0);
		}
		
		template<class Iterator>
		constexpr Matrix& assign_column(Iterator first){
			return assign_column(first, 0, 0);
		}
		
		template<size_t sub_size>
		constexpr Vector<T, sub_size> sub_row(ptrdiff_t row_offset, ptrdiff_t column_offset)const{
			row_offset = (row_offset >= 0) ? row_offset : row_offset + this->rows();
			column_offset = (column_offset >= 0) ? column_offset : column_offset + this->columns();
			Vector<T, sub_size> result;
			size_t result_i = 0;
			size_t this_column = column_offset;
			for(; result_i < sub_size && this_column < this->columns(); ++result_i, (void)++this_column){
				result.at(result_i) = this->at(row_offset, this_column);
			}
			return result;
		}

		constexpr Vector<T, _columns> row(ptrdiff_t row_offset)const{
			return sub_row<_columns>(row_offset, 0);
		}
		
		template<size_t sub_size>
		constexpr Matrix& assign_sub_row(const Vector<T, sub_size>& row, ptrdiff_t row_number, ptrdiff_t column_offset){
			row_number = (row_number >= 0) ? row_number : row_number + this->rows();
			column_offset = (column_offset >= 0) ? column_offset : column_offset + this->columns();
			size_t result_i = 0;
			size_t this_column = column_offset;
			for(; result_i < sub_size && this_column < _rows; ++result_i, (void)++this_column){
				this->at(row_number, this_column) = row.at(result_i);
			}
			return *this;
		}
		
		constexpr Matrix& assign_row(const Vector<T, _columns>& row, ptrdiff_t row_number){
			row_number = (row_number >= 0) ? row_number : row_number + this->rows();
			return this->assign_sub_row<_columns>(row, row_number, 0);
		}
		
		template<class Iterator>
		constexpr Matrix& assign_row(Iterator first, Iterator last, ptrdiff_t row_number, ptrdiff_t column_offset){
			row_number = (row_number >= 0) ? row_number : row_number + this->rows();
			column_offset = (column_offset >= 0) ? column_offset : column_offset + this->columns();
			size_t column = column_offset;
			for(;first != last && column < _columns; ++first, (void)++column){
				this->at(row_number, column) = *first;
			}
			return *this;
		}
		
		template<class Iterator>
		constexpr Matrix& assign_row(Iterator first, Iterator last, ptrdiff_t row_offset){
			return this->assign_row(first, last, 0, row_offset);
		}
		
		template<class Iterator>
		constexpr Matrix& assign_row(Iterator first, Iterator last){
			return this->assign_row(first, last, 0, 0);
		}
		
		template<class Iterator>
		constexpr Matrix& assign_row(Iterator first, ptrdiff_t row_offset, ptrdiff_t column_offset){
			row_offset = (row_offset >= 0) ? row_offset : row_offset + this->rows();
			column_offset = (column_offset >= 0) ? column_offset : column_offset + this->columns();
			size_t column = column_offset;
			for(; column < _columns; ++first, (void)++column){
				this->at(row_offset, column) = *first;
			}
			return *this;
		}
		
		template<class Iterator>
		constexpr Matrix& assign_row(Iterator first, ptrdiff_t row_offset){
			row_offset = (row_offset >= 0) ? row_offset : row_offset + this->rows();
			return this->assign_row(first, row_offset, 0);
		}
		
		template<class Iterator>
		constexpr Matrix& assign_row(Iterator first){
			return this->assign_row(first, 0, 0);
		}
		
		template<class Ty, size_t sub_rows, size_t sub_columns>
		constexpr Matrix& assign(const Matrix<Ty, sub_rows, sub_columns>& sub_matrix, ptrdiff_t row_offset=0, ptrdiff_t column_offset=0){
			row_offset = (row_offset >= 0) ? row_offset : row_offset + this->rows();
			column_offset = (column_offset >= 0) ? column_offset : column_offset + this->columns();
			size_t this_row = row_offset;
			size_t sm_row = 0;
			for(; this_row < _rows && sm_row < sub_rows; ++this_row, (void)++sm_row){
				size_t this_column = column_offset;
				size_t sm_column = 0;
				for(; this_column < _columns && sm_column < sub_columns; ++this_column, (void)++sm_column){
					this->at(this_row, this_column) = sub_matrix.at(sm_row, sm_column);
				}
			}
			return *this;
		}

		constexpr Matrix& fill(const T& value){
			for(T& elem : *this) elem = value;
			return *this;
		}
		
		constexpr Matrix& fill(const T& value, ptrdiff_t row_offset, ptrdiff_t column_offset, ptrdiff_t nrows=-1, ptrdiff_t ncolumns=-1){
			row_offset = (row_offset >= 0) ? row_offset : row_offset + this->rows();
			column_offset = (column_offset >= 0) ? column_offset : column_offset + this->columns();
			nrows = (nrows>=0) ? nrows : this->rows()+nrows+1;
			ncolumns = (ncolumns>=0) ? ncolumns : this->rows()+ncolumns+1;
			for(size_t row = row_offset; row != (row_offset + nrows) && row < _rows; ++row){
				for(size_t column = column_offset; column != (column_offset + ncolumns) && column < _columns; ++column){
					this->at(row, column) = value;
				}
			}
			return *this;
		}

		constexpr Matrix& set_zero(){
			return this->fill(T(0));
		}
		constexpr Matrix& set_zero(ptrdiff_t row_offset, ptrdiff_t column_offset, ptrdiff_t nrows=-1, ptrdiff_t ncolumns=-1){
			return this->fill(T(0), row_offset, column_offset, nrows, ncolumns);
		}
		
		constexpr Matrix& set_one(){
			return this->fill(T(1));
		}
		constexpr Matrix& set_one(ptrdiff_t row_offset, ptrdiff_t column_offset, ptrdiff_t nrows=-1, ptrdiff_t ncolumns=-1){
			return this->fill(T(0), row_offset, column_offset, nrows, ncolumns);
		}
		
		constexpr Matrix& set_unity(){
			for(size_t row = 0; row < _rows; ++row){
				for(size_t column = 0; column < _columns; ++column){
					this->at(row, column) = (row == column) ? T(1) : T(0);
				}
			}
			return *this;
		}
		
		constexpr Matrix& set_unity(ptrdiff_t row_offset, ptrdiff_t column_offset){
			row_offset = (row_offset >= 0) ? row_offset : row_offset + this->rows();
			column_offset = (column_offset >= 0) ? column_offset : column_offset + this->columns();

			for(size_t row = 0; (row + row_offset) < this->rows() && (row + column_offset) < this->columns(); ++row) {
				for(size_t column = 0; (column + column_offset) < this->columns(); ++column){
					this->at(row + row_offset, column + column_offset) = (row == column) ? T(1) : T(0);
				}
			}
			return *this;
		}
		
		constexpr Matrix& set_unity(ptrdiff_t row_offset, ptrdiff_t column_offset, ptrdiff_t nrows, ptrdiff_t ncolumns){
			row_offset = (row_offset >= 0) ? row_offset : row_offset + this->rows();
			column_offset = (column_offset >= 0) ? column_offset : column_offset + this->columns();
			nrows = (nrows>=0)?nrows:this->rows()-nrows+1;
			ncolumns = (ncolumns>=0)?ncolumns:this->rows()-ncolumns+1;
			for(size_t row = 0; row < nrows && (row + row_offset) < this->rows(); ++row){
				for(size_t column = 0; column < ncolumns && (column + column_offset) < this->columns(); ++column){
					this->at(row + row_offset, column + column_offset) = (row == column) ? T(1) : T(0);
				}
			}
			return *this;
		}
		
		// conversion to vector
		template <class U = T, TWMATH_ENABLE_IF((std::is_same<T, U>::value && _rows > 1 && _columns == 1))>
		constexpr operator Vector<U, _rows> () const {
			return Vector<U, _rows>(this->begin(), this->end()); 
		}
		
		// conversion to vector
		template <class U = T, TWMATH_ENABLE_IF((std::is_same<T, U>::value && _rows == 1 && _columns > 1))>
		constexpr operator Vector<U, _columns> () const {
			return Vector<U, _columns>(this->begin(), this->end());
		}
		
		// conversion to scalar
		template <class U = T, TWMATH_ENABLE_IF((std::is_same<T, U>::value && _rows == 1 && _columns == 1))>
		constexpr operator U () const {
			return U(*this->begin());
		}
		
	};
	
	template<class T, size_t rows, size_t columns>
	constexpr MatrixRef<T, rows, columns> transpose(Matrix<T, rows, columns>& M){
		return MatrixRef<T, rows, columns>(M.data(), M.row_increment(), M.column_increment(), !M.is_transposed());
	}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns, class F>
	constexpr Matrix<F, _rows, _columns> 
	for_each (const Matrix<Tl, _rows, _columns>& lhs, const Matrix<Tr, _rows, _columns>& rhs, F (*f)(const Tl&, const Tr&)){
		Matrix<F, _rows, _columns> result;
		auto result_itr = result.begin();
		auto lhs_itr = lhs.begin();
		auto rhs_itr = rhs.begin();
		for(; result_itr != result.end(); ++result_itr, (void)++lhs_itr, (void)++rhs_itr) *result_itr = f(*lhs_itr, *rhs_itr);
		return result;
	}
	template<class Tl, class Tr, size_t _rows, size_t _columns, class F>
	constexpr Matrix<F, _rows, _columns> 
	for_each (const Tl& lhs, const Matrix<Tr, _rows, _columns>& rhs, F (*f)(const Tl&, const Tr&)){
		Matrix<F, _rows, _columns> result;
		auto result_itr = result.begin();
		auto rhs_itr = rhs.begin();
		for(; result_itr != result.end(); ++result_itr, (void)++rhs_itr) *result_itr = f(lhs, *rhs_itr);
		return result;
	}
	template<class Tl, class Tr, size_t _rows, size_t _columns, class F>
	constexpr Matrix<F, _rows, _columns> 
	for_each (const Matrix<Tl, _rows, _columns>& lhs, const Tr& rhs, F (*f)(const Tl&, const Tr&)){
		Matrix<F, _rows, _columns> result;
		auto result_itr = result.begin();
		auto lhs_itr = lhs.begin();
		for(; result_itr != result.end(); ++result_itr, (void)++lhs_itr) *result_itr = f(*lhs_itr, rhs);
		return result;
	}
	template<class T, size_t _rows, size_t _columns, class F>
	constexpr Matrix<F, _rows, _columns> 
	for_each (const Matrix<T, _rows, _columns>& v, F (*f)(const T&)){
		Matrix<F, _rows, _columns> result;
		auto result_itr = result.begin();
		auto v_itr = v.begin();
		for(; result_itr != result.end(); ++result_itr, (void)++v_itr) *result_itr = f(*v_itr);
		return result;
	}
	
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::add<Tl, Tr>), Tl, Tr), _rows, _columns>  
	add (const Matrix<Tl, _rows, _columns>& lhs, const Matrix<Tr, _rows, _columns>& rhs){return for_each(lhs, rhs, twmath_base::add<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::add<Tl, Tr>), Tl, Tr), _rows, _columns>  
	add (const Matrix<Tl, _rows, _columns>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::add<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::add<Tl, Tr>), Tl, Tr), _rows, _columns>  
	add (const Tl& lhs, const Matrix<Tr, _rows, _columns>& rhs){return for_each(lhs, rhs, twmath_base::add<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::add<Tl, Tr>), Tl, Tr), _rows, _columns>  
	operator+ (const Matrix<Tl, _rows, _columns>& lhs, const Matrix<Tr, _rows, _columns>& rhs){return for_each(lhs, rhs, twmath_base::add<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::add<Tl, Tr>), Tl, Tr), _rows, _columns>  
	operator+ (const Matrix<Tl, _rows, _columns>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::add<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::add<Tl, Tr>), Tl, Tr), _rows, _columns>  
	operator+ (const Tl& lhs, const Matrix<Tr, _rows, _columns>& rhs){return for_each(lhs, rhs, twmath_base::add<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows>
	constexpr Vector<result_type2((twmath_base::add<Tl, Tr>), Tl, Tr), _rows>  
	operator+ (const Matrix<Tl, _rows, 1>& lhs, const Vector<Tr, _rows>& rhs){
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
	operator+ (const Matrix<Tl, 1, _columns>& lhs, const Vector<Tr, _columns>& rhs){
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
	operator+ (const Matrix<Tl, 1, 1>& lhs, const Vector<Tr, 1>& rhs){
		return lhs.at(0, 0) + rhs.at(0);
	}
	
	template<class Tl, class Tr, size_t _rows>
	constexpr Vector<result_type2((twmath_base::add<Tl, Tr>), Tl, Tr), _rows>  
	operator+ (const Vector<Tr, _rows>& lhs, const Matrix<Tl, _rows, 1>& rhs){
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
	operator+ (const Vector<Tr, _columns>& lhs, const Matrix<Tl, 1, _columns>& rhs){
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
	operator+ (const Vector<Tr, 1>& lhs, const Matrix<Tl, 1, 1>& rhs){
		return lhs.at(0) + rhs.at(0, 0);
	}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::sub<Tl, Tr>), Tl, Tr), _rows, _columns>  
	sub (const Matrix<Tl, _rows, _columns>& lhs, const Matrix<Tr, _rows, _columns>& rhs){return for_each(lhs, rhs, twmath_base::sub<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::sub<Tl, Tr>), Tl, Tr), _rows, _columns>  
	sub (const Matrix<Tl, _rows, _columns>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::sub<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::sub<Tl, Tr>), Tl, Tr), _rows, _columns>  
	sub (const Tl& lhs, const Matrix<Tr, _rows, _columns>& rhs){return for_each(lhs, rhs, twmath_base::sub<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::sub<Tl, Tr>), Tl, Tr), _rows, _columns>  
	operator- (const Matrix<Tl, _rows, _columns>& lhs, const Matrix<Tr, _rows, _columns>& rhs){return for_each(lhs, rhs, twmath_base::sub<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::sub<Tl, Tr>), Tl, Tr), _rows, _columns>  
	operator- (const Matrix<Tl, _rows, _columns>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::sub<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::sub<Tl, Tr>), Tl, Tr), _rows, _columns>  
	operator- (const Tl& lhs, const Matrix<Tr, _rows, _columns>& rhs){return for_each(lhs, rhs, twmath_base::sub<Tl, Tr>);}
	
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::mul<Tl, Tr>), Tl, Tr), _rows, _columns>  
	mul (const Matrix<Tl, _rows, _columns>& lhs, const Matrix<Tr, _rows, _columns>& rhs){return for_each(lhs, rhs, twmath_base::mul<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::mul<Tl, Tr>), Tl, Tr), _rows, _columns>  
	mul (const Matrix<Tl, _rows, _columns>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::mul<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::mul<Tl, Tr>), Tl, Tr), _rows, _columns> 
	mul (const Tl& lhs, const Matrix<Tr, _rows, _columns>& rhs){return for_each(lhs, rhs, twmath_base::mul<Tl, Tr>);}

	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::mul<Tl, Tr>), Tl, Tr), _rows, _columns> 
	operator* (const Matrix<Tl, _rows, _columns>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::mul<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::mul<Tl, Tr>), Tl, Tr), _rows, _columns> 
	operator* (const Tl& lhs, const Matrix<Tr, _rows, _columns>& rhs){return for_each(lhs, rhs, twmath_base::mul<Tl, Tr>);}


	template<class Tl, class Tr, size_t lcols_rrows, size_t lrows, size_t rcols>
	constexpr Matrix<result_type2((twmath_base::mul<Tl, Tr>), Tl, Tr), rcols, lrows> 
	matmul (const Matrix<Tl, lrows, lcols_rrows>& lhs, const Matrix<Tr, lcols_rrows, rcols>& rhs){
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
	operator* (const Matrix<Tl, lrows, lcols_rrows>& lhs, const Matrix<Tr, lcols_rrows, rcols>& rhs){return matmul(lhs, rhs);}
	
	template<class Tl, class Tr, size_t rows, size_t cols>
	constexpr Vector<result_type2((twmath_base::mul<Tl, Tr>), Tl, Tr), rows> mul (const Matrix<Tl, rows, cols>& l, const Vector<Tr, cols>& r){
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
	constexpr Vector<result_type2((twmath_base::mul<Tl, Tr>), Tl, Tr), rows> operator* (const Matrix<Tl, rows, cols>& l, const Vector<Tr, cols>& r){
		return mul<Tl, Tr, rows, cols>(l, r);
	}
	
	template<class Tl, class Tr, size_t rows, size_t cols>
	constexpr Vector<result_type2((twmath_base::mul<Tl, Tr>), Tl, Tr), cols> mul (const Vector<Tr, rows>& l, const Matrix<Tl, rows, cols>& r){
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
	constexpr Vector<result_type2((twmath_base::mul<Tl, Tr>), Tl, Tr), rows> operator* (const Vector<Tr, rows>& l, const Matrix<Tl, rows, cols>& r){
		return mul<Tl, Tr, rows, cols>(l, r);
	}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::div<Tl, Tr>), Tl, Tr), _rows, _columns> 
	div (const Matrix<Tl, _rows, _columns>& lhs, const Matrix<Tr, _rows, _columns>& rhs){return for_each(lhs, rhs, twmath_base::div<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::div<Tl, Tr>), Tl, Tr), _rows, _columns> 
	div (const Matrix<Tl, _rows, _columns>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::div<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::div<Tl, Tr>), Tl, Tr), _rows, _columns>  
	div (const Tl& lhs, const Matrix<Tr, _rows, _columns>& rhs){return for_each(lhs, rhs, twmath_base::div<Tl, Tr>);}
	
	template<class Tl, class Tr, size_t _rows, size_t _columns>
	constexpr Matrix<result_type2((twmath_base::div<Tl, Tr>), Tl, Tr), _rows, _columns>   
	operator/ (const Matrix<Tl, _rows, _columns>& lhs, const Tr& rhs){return for_each(lhs, rhs, twmath_base::div<Tl, Tr>);}
	
	template<class T, size_t _rows, size_t _columns>
	constexpr Matrix<T, _rows, _columns> zeros(){
		Matrix<T, _rows, _columns> result;
		for(const auto& elem : result) elem = T(0);
		return result;
	}
	
	template<class T, size_t _rows, size_t _columns>
	constexpr Matrix<T, _rows, _columns> zeros_like(const Matrix<T, _rows, _columns>& A){
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
	constexpr Matrix<T, _rows, _columns> unity_like(const Matrix<T, _rows, _columns>& A){
		return unity<T, _rows, _columns>();
	}
	
	template<class Stream, class T, size_t _rows, size_t _columns>
	Stream& print_pretty(Stream& stream, const Matrix<T, _rows, _columns>& M, const char* name = nullptr, const char* indentation=""){
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
	Stream& print_csv(Stream& stream, const Matrix<T, _rows, _columns>& m, const char* value_separator=", ", const char* line_separator="\n") {
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
	Stream& operator<< (Stream& stream, const Matrix<T, _rows, _columns>& m){
		return print_csv(stream, m);
	}
	
	// approximates the exponential of the matrix with the tailor expansion:
	// e^A = I + A + A^2 / 2 + ... + A^n / n!
	// the minimum number of n is 2. If n is set lower than 2, then 2 increments will be calculated regardless
	template<class T, size_t _rows, size_t _columns>
	constexpr Matrix<T, _rows, _columns> exp_taylor(const Matrix<T, _rows, _columns>& A, size_t n = 4){
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
	constexpr Matrix<T, _rows, _columns> exp_taylor_square_scale(const Matrix<T, _rows, _columns>& x, size_t taylor_order = 4, size_t scaling = 4){
		Matrix<T, _rows, _columns> t = exp_taylor(x/(1 << scaling), taylor_order);
		for(size_t i = 0; i < scaling; ++i) t = t * t;
		return t;
	}
	
	// default matrix exponentiation
	template<class T, size_t _rows, size_t _columns>
	constexpr Matrix<T, _rows, _columns> exp(const Matrix<T, _rows, _columns>& x){
		return exp_taylor_square_scale(x);
	}
	
}