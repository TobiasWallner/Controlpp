#pragma once

#
#include "MatrixToken.hpp"

namespace twmath{
	template<class T, size_t _rows, size_t _columns>
	class MatrixRef{
	private:
		T* _ptr = nullptr;
		const size_t _row_increment = _columns;
		const size_t _column_increment = 1;
		const bool _transposed = false;
		
	public:
		constexpr MatrixRef(const MatrixRef& other) = default;
		
		constexpr MatrixRef(T* ptr, const size_t row_increment = 1, const size_t column_increment = 1, const bool transposed = false)
			: _ptr(ptr)
			, _row_increment(row_increment)
			, _column_increment(column_increment)
			, _transposed(transposed){
				twmath_assert(row_increment != 0, "Matrix::MatrixReference Error: row_increment cannot be zero.")
				twmath_assert(column_increment != 0, "Matrix::MatrixReference Error: column_increment cannot be zero.")
			}
			
		constexpr T& at(size_t row, size_t column){
			if(!this->_transposed){
				twmath_assert(row < _rows, "Matrix::MatrixReference Error: Out of bound access. Tried to access row " << row << ", but matrix has only " << _rows << "rows.")
				twmath_assert(column < _columns, "Matrix::MatrixReference Error: Out of bound access. Tried to access column " << column << ", but matrix has only " << _columns << "columns.")
				const size_t element_increment = this->_row_increment * row + this->_column_increment * column;
				return *(this->_ptr + element_increment);
			}else{
				twmath_assert(column < _rows, "Matrix::MatrixReference Error: Out of bound access. Tried to access row " << row << ", but matrix has only " << _rows << "rows.")
				twmath_assert(row < _columns, "Matrix::MatrixReference Error: Out of bound access. Tried to access column " << column << ", but matrix has only " << _columns << "columns.")
				const size_t element_increment = this->_row_increment * column + this->_column_increment * row;
				return *(this->_ptr + element_increment);
			}
		}
		
		constexpr const T& at(size_t row, size_t column) const {
			if(!this->_transposed){
				twmath_assert(row < _rows, "Matrix::MatrixReference Error: Out of bound access. Tried to access row " << row << ", but matrix has only " << _rows << "rows.")
				twmath_assert(column < _columns, "Matrix::MatrixReference Error: Out of bound access. Tried to access column " << column << ", but matrix has only " << _columns << "columns.")
				const size_t element_increment = this->_row_increment * row + this->_column_increment * column;
				return *(this->_ptr + element_increment);
			}else{
				twmath_assert(column < _rows, "Matrix::MatrixReference Error: Out of bound access. Tried to access row " << row << ", but matrix has only " << _rows << "rows.")
				twmath_assert(row < _columns, "Matrix::MatrixReference Error: Out of bound access. Tried to access column " << column << ", but matrix has only " << _columns << "columns.")
				const size_t element_increment = this->_row_increment * column + this->_column_increment * row;
				return *(this->_ptr + element_increment);
			}
		}
		
		constexpr size_t row_increment() const {return this->_row_increment;}
		constexpr size_t column_increment() const {return this->_column_increment;}
		
		constexpr bool is_transposed() const {return this->_transposed;}
		
		constexpr size_t columns() const {return (_transposed) ? _rows : _columns;}
		constexpr size_t rows() const {return (_transposed) ? _columns : _rows;}
		size_t size() const {return this->rows() * this->columns();}
		
		constexpr T* data(){return this->_ptr;}
		constexpr const T* data() const {return this->_ptr;}
		
		template<class Tm>
		MatrixRef& assign(const MatrixRef<Tm, _rows, _columns>& M){
			for(size_t r = 0; r < _rows; ++r){
				for(size_t c = 0; c < _columns; ++c){
					this->at(r, c) = M.at(r, c);
				}
			}
		}
		
		template<class Tm>
		MatrixRef& operator=(const MatrixRef<Tm, _rows, _columns>& M){return this->assign(M);}
		
		template<size_t rows, size_t columns, TWMATH_ENABLE_IF(rows != 0 && columns != 0)>
		constexpr MatrixRef<T, rows, columns> sub_matrix(size_t row_offset=0, size_t column_offset=0, size_t row_stride=1, size_t column_stride=1) {
			twmath_assert(row_stride != 0, "Matrix::MatrixReference Error: row_stride cannot be zero.")
			twmath_assert(column_stride != 0, "Matrix::MatrixReference Error: column_stride cannot be zero.")
			auto ptr = &this->at(row_offset, column_offset);
			MatrixRef<T, rows, columns> result(ptr, row_stride * this->_row_increment, column_stride * this->_column_increment);
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
			MatrixRef<T, (_rows < _columns) ? _rows : _columns, 1> result (this->data(), this->_row_increment+1, 1);
			return result;
		}
	};

	
	
}// namespace twmath