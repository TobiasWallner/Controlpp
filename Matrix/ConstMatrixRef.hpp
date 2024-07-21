#pragma once

#include "definitions.hpp"

namespace twmath{
	
	template<class T, size_t _rows, size_t _columns, bool _transposed=false>
	class ConstMatrixRef : public StaticMatrixToken{
	private:
		const T* const _ptr = nullptr;
		const size_t _row_increment = _columns;
		const size_t _column_increment = 1;
		
	public:
		constexpr ConstMatrixRef(const ConstMatrixRef& other) = default;
		
		constexpr ConstMatrixRef(const T* ptr, const size_t row_increment = 1, const size_t column_increment = 1)
			: _ptr(ptr)
			, _row_increment(row_increment)
			, _column_increment(column_increment){
				twmath_assert(row_increment != 0, "Matrix::MatrixReference Error: row_increment cannot be zero.")
				twmath_assert(column_increment != 0, "Matrix::MatrixReference Error: column_increment cannot be zero.")
			}
		
		//------------------- const at ---------------------
		
		constexpr const T& at(size_t row, size_t col) const {
			twmath_assert(col < this->columns(), "Column '" << col << "' is out of bound of '" << this->columns() << "'.");
			twmath_assert(row < this->rows(), "Row '" << row << "' is out of bound of '" << this->rows() << "'.");
			return *(_ptr + ((!_transposed) ? col : row) + ((!_transposed) ? row : col) * _columns);
		}
		
		template<class UInteger, TWMATH_ENABLE_IF(std::is_unsigned_v<UInteger>)>
		constexpr const T& at(UInteger row, UInteger col) const {
			return this->at(static_cast<size_t>(row), static_cast<size_t>(col));
		}
		
		template<class UInteger, class SInteger, TWMATH_ENABLE_IF(std::is_unsigned_v<UInteger>), TWMATH_ENABLE_IF(std::is_signed_v<SInteger>)>
		constexpr const T& at(UInteger row, SInteger col) const {
			return this->at(static_cast<size_t>(row), static_cast<size_t>((col >= 0) ? col : col + this->columns()));
		}
		
		template<class UInteger, class SInteger, TWMATH_ENABLE_IF(std::is_unsigned_v<UInteger>), TWMATH_ENABLE_IF(std::is_signed_v<SInteger>)>
		constexpr const T& at(SInteger row, UInteger col) const {
			return this->at(static_cast<size_t>((row >= 0) ? row : row + this->rows()), static_cast<size_t>(col));
		}
		
		template<class SInteger, TWMATH_ENABLE_IF(std::is_signed_v<SInteger>)>
		constexpr const T& at(SInteger row, SInteger col) const {
			return this->at(static_cast<size_t>((row >= 0) ? row : row + this->rows()), static_cast<size_t>((col >= 0) ? col : col + this->columns()));
		}
		
		// ----------------------- iterators -----------------------------
		// TODO: ?
		
		// ----------------------- size / capacity / co. -----------------------
		
		constexpr size_t row_increment() const {return this->_row_increment;}
		constexpr size_t column_increment() const {return this->_column_increment;}
		
		constexpr size_t columns() const {return (_transposed) ? _rows : _columns;}
		static constexpr scolumns(){return (_transposed) ? _rows : _columns;}
		constexpr mem_columns()const{return _columns;} // returns the number of columns as layed out in memory
		static constexpr smem_columns(){return _columns;} // returns the number of columns as layed out in memory
		
		constexpr size_t rows() const {return (_transposed) ? _columns : _rows;}
		static constexpr size_t srows() {return (_transposed) ? _columns : _rows;}
		constexpr mem_rows()const{return _rows;} // returns the number of rows as layed out in memory
		static constexpr smem_rows(){return _rows;} // returns the number of rows as layed out in memory
		
		constexpr size_t transposed() const {return _transposed;}
		static constexpr size_t stransposed() {return _transposed;}
		
		constexpr size_t size() const {return _rows * _columns;}
		static constexpr size_t ssize() {return _rows * _columns;}
		
		constexpr const T* data() const {return this->_ptr;}
		constexpr const T* cdata() const {return this->_ptr;}
		
		//------------------- sub_matrix ---------------------

		template<size_t sub_rows, size_t sub_columns>
		constexpr ConstMatrixRef<T, sub_rows, sub_columns> csub_matrix(size_t row_offset=0, size_t column_offset=0, size_t row_stride=1, size_t column_stride=1) const {
			twmath_assert(row_stride!=0, "Matrix Error: row_stride cannot be zero.");
			twmath_assert(column_stride!=0, "Matrix Error: column_stride cannot be zero.");
			T* ptr = &this->at(row_offset, column_offset);
			const size_t sub_column_increment = column_stride;
			const size_t sub_row_increment = row_stride * this->columns();
			ConstMatrixRef<T, sub_rows, sub_columns> result(ptr, sub_row_increment, sub_column_increment);
			return result;
		}
		template<size_t sub_rows, size_t sub_columns>
		constexpr ConstMatrixRef<T, sub_rows, sub_columns> sub_matrix(size_t row_offset=0, size_t column_offset=0, size_t row_stride=1, size_t column_stride=1) const {
			return this->csub_matrix(row_offset, column_offset, row_stride, column_stride);
		}

		// ----------------------- column vector slice -----------------------------

		template<size_t rows=_rows>
		constexpr ConstVectorRef<T, rows> ccolumn(size_t column_number=0, size_t offset=0, size_t stride=1) const {
			VectorRef<T, rows> result(&this->at(offset, column_number), stride * this->columns());
			return result;
		}
		template<size_t rows=_rows>
		constexpr ConstVectorRef<T, rows> column(size_t column_number=0, size_t offset=0, size_t stride=1) const {
			return this->ccolumn(column_number, offset, stride);
		}
		
		// ----------------------- row vector slice -----------------------------
		
		template<size_t columns=_columns>
		constexpr ConstVectorRef<T, columns> crow(size_t row_number=0, size_t offset=0, size_t stride=1) const {
			VectorRef<T, columns> result(&this->at(row_number, offset), stride);
			return result;
		}
		template<size_t columns=_columns>
		constexpr ConstVectorRef<T, columns> row(size_t row_number=0, size_t offset=0, size_t stride=1) const {
			return this->crow(row_number, offset, stride);
		}
		constexpr ConstVectorRef<T, _columns> operator[](size_t row_number) const {return this->row(row_number);}
		
		// ----------------------- diag -----------------------------
		
		constexpr ConstVectorRef<T, (_rows < _columns) ? _rows : _columns> cdiag() const {
			VectorRef<T, (_rows < _columns) ? _rows : _columns> result(this->data(), this->columns() + 1);
			return result;
		}
		constexpr VectorRef<T, (_rows < _columns) ? _rows : _columns> diag() const {return this->cdiag();}
	
		// ---------------------- conversion to scalar ---------------------------
		
		template <class U = T, TWMATH_ENABLE_IF((std::is_same_v<T, U> && _rows == 1 && _columns == 1))>
		constexpr operator const U& () const {return this->at(0, 0);}
		
		template <class U = T, TWMATH_ENABLE_IF((std::is_same_v<T, U> && _rows == 1 && _columns == 1))>
		constexpr operator U () const {return U(this->at(0, 0));}
	
	};

	template<class T, size_t rows, size_t columns> struct value_type<ConstMatrixRef<T, rows, columns>>{using type = T;};
	
}// namespace twmath