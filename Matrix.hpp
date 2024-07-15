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
	
}