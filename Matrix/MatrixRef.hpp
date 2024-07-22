#pragma once

#include "definitions.hpp"

#include "MatrixTraits.hpp"
#include "ConstMatrixRef.hpp"

namespace twmath{
	template<class T, size_t _rows, size_t _columns, bool _transposed=false>
	class MatrixRef : public StaticMatrixToken{
	private:
		T* const _ptr = nullptr;
		const size_t _row_increment = _columns;
		const size_t _column_increment = 1;
		
	public:
		constexpr MatrixRef(const MatrixRef& other) = default;
		
		constexpr MatrixRef(T* ptr, const size_t row_increment = 1, const size_t column_increment = 1)
			: _ptr(ptr)
			, _row_increment(row_increment)
			, _column_increment(column_increment){
				twmath_assert(row_increment != 0, "Matrix::MatrixReference Error: row_increment cannot be zero.")
				twmath_assert(column_increment != 0, "Matrix::MatrixReference Error: column_increment cannot be zero.")
			}
			
		
		//------------------- at ---------------------
		
		constexpr T& at(size_t row, size_t col){
			twmath_assert(col < this->columns(), "Column '" << col << "' is out of bound of '" << this->columns() << "'.");
			twmath_assert(row < this->rows(), "Row '" << row << "' is out of bound of '" << this->rows() << "'.");
			return *(_ptr + ((!_transposed) ? col : row) * _column_increment + ((!_transposed) ? row : col) * _row_increment);
		}

		template<class UInteger, TWMATH_ENABLE_IF(std::is_unsigned_v<UInteger>)>
		constexpr T& at(UInteger row, UInteger col){
			return this->at(static_cast<size_t>(row), static_cast<size_t>(col));
		}
		
		template<class UInteger, class SInteger, TWMATH_ENABLE_IF(std::is_unsigned_v<UInteger>), TWMATH_ENABLE_IF(std::is_signed_v<SInteger>)>
		constexpr T& at(UInteger row, SInteger col){
			return this->at(static_cast<size_t>(row), static_cast<size_t>((col >= 0) ? col : col + this->columns()));
		}
		
		template<class UInteger, class SInteger, TWMATH_ENABLE_IF(std::is_unsigned_v<UInteger>), TWMATH_ENABLE_IF(std::is_signed_v<SInteger>)>
		constexpr T& at(SInteger row, UInteger col){
			return this->at(static_cast<size_t>((row >= 0) ? row : row + this->rows()), static_cast<size_t>(col));
		}
		
		template<class SInteger, TWMATH_ENABLE_IF(std::is_signed_v<SInteger>)>
		constexpr T& at(SInteger row, SInteger col){
			return this->at(static_cast<size_t>((row >= 0) ? row : row + this->rows()), static_cast<size_t>((col >= 0) ? col : col + this->columns()));
		}
		
		//------------------- const at ---------------------
		
		constexpr const T& at(size_t row, size_t col) const {
			twmath_assert(col < this->columns(), "Column '" << col << "' is out of bound of '" << this->columns() << "'.");
			twmath_assert(row < this->rows(), "Row '" << row << "' is out of bound of '" << this->rows() << "'.");
			return *(_ptr + ((!_transposed) ? col : row) * _column_increment + ((!_transposed) ? row : col) * _row_increment);
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
		
		//------------------- size / capacity / co. ---------------------
		
		constexpr size_t row_increment() const {return this->_row_increment;}
		constexpr size_t column_increment() const {return this->_column_increment;}
		
		constexpr size_t columns() const {return (_transposed) ? _rows : _columns;}
		static size_t constexpr scolumns(){return (_transposed) ? _rows : _columns;}
		constexpr size_t mem_columns()const{return _columns;} // returns the number of columns as layed out in memory
		static size_t constexpr smem_columns(){return _columns;} // returns the number of columns as layed out in memory
		
		constexpr size_t rows() const {return (_transposed) ? _columns : _rows;}
		static constexpr size_t srows() {return (_transposed) ? _columns : _rows;}
		constexpr size_t mem_rows()const{return _rows;} // returns the number of rows as layed out in memory
		static constexpr size_t smem_rows(){return _rows;} // returns the number of rows as layed out in memory
		
		constexpr size_t transposed() const {return _transposed;}
		static constexpr size_t stransposed() {return _transposed;}
		
		constexpr size_t size() const {return _rows * _columns;}
		static constexpr size_t ssize() {return _rows * _columns;}
		
		//------------------- data ---------------------
		
		constexpr T* data(){return this->_ptr;}
		constexpr const T* data() const {return this->_ptr;}
		constexpr const T* cdata() const {return this->_ptr;}
		
		//------------------- assignment ---------------------
		
		template<class M, TWMATH_ENABLE_IF(is_static_matrix_v<M> && M::srows() == srows() && M::scolumns() == scolumns())>
		MatrixRef& assign(const M& m){
			for(size_t r = 0; r < _rows; ++r){
				for(size_t c = 0; c < _columns; ++c){
					this->at(r, c) = m.at(r, c);
				}
			}
			return *this;
		}
		
		template<class M, TWMATH_ENABLE_IF(is_static_matrix_v<M> && M::srows() == srows() && M::scolumns() == scolumns())>
		MatrixRef& operator=(const M& m){return this->assign(m);}
		
		// ----------------------- iterators -----------------------------
		// TODO: ?
		
		//------------------- sub_matrix ---------------------
		
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
		constexpr VectorRef<T, rows> column(size_t column_number=0, size_t offset=0, size_t stride=1) {
			VectorRef<T, rows> result(&this->at(offset, column_number), stride * this->columns());
			return result;
		}
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
		constexpr VectorRef<T, columns> row(size_t row_number=0, size_t offset=0, size_t stride=1) {
			VectorRef<T, columns> result(&this->at(row_number, offset), stride);
			return result;
		}
		template<size_t columns=_columns>
		constexpr ConstVectorRef<T, columns> crow(size_t row_number=0, size_t offset=0, size_t stride=1) const {
			ConstVectorRef<T, columns> result(&this->at(row_number, offset), stride);
			return result;
		}
		template<size_t columns=_columns>
		constexpr ConstVectorRef<T, columns> row(size_t row_number=0, size_t offset=0, size_t stride=1) const {
			return this->crow(row_number, offset, stride);
		}
		constexpr VectorRef<T, _columns> operator[](size_t row_number){return this->row(row_number);}
		constexpr ConstVectorRef<T, _columns> operator[](size_t row_number) const {return this->row(row_number);}
		
		// ----------------------- diag -----------------------------
		
		constexpr VectorRef<T, (_rows < _columns) ? _rows : _columns> diag() {
			VectorRef<T, (_rows < _columns) ? _rows : _columns> result(this->data(), this->columns() + 1);
			return result;
		}
		constexpr ConstVectorRef<T, (_rows < _columns) ? _rows : _columns> cdiag() const {
			ConstVectorRef<T, (_rows < _columns) ? _rows : _columns> result(this->data(), this->columns() + 1);
			return result;
		}
		constexpr ConstVectorRef<T, (_rows < _columns) ? _rows : _columns> diag() const {return this->cdiag();}
		
		// ----------------------- setters -----------------------------
		
		template<class Ty>
		constexpr MatrixRef& fill(const Ty& value){
			for(size_t row = 0; row < this->rows(); ++row) 
				for(size_t col = 0; col < this->columns(); ++col)
					this->at(row, col) = value;
			return *this;
		}
		
		constexpr MatrixRef& set_zero(){return this->fill(0);}
		constexpr MatrixRef& set_one(){return this->fill(1);}
		constexpr MatrixRef& set_unity(){this->set_zero().diag().set_one(); return *this;}
		
		// ---------------------- conversion to scalar ---------------------------
		template <class U = T, TWMATH_ENABLE_IF((std::is_same_v<T, U> && _rows == 1 && _columns == 1))>
		constexpr operator U& () {return this->at(0, 0);}
		
		template <class U = T, TWMATH_ENABLE_IF((std::is_same_v<T, U> && _rows == 1 && _columns == 1))>
		constexpr operator const U& () const {return this->at(0, 0);}
		
		template <class U = T, TWMATH_ENABLE_IF((std::is_same_v<T, U> && _rows == 1 && _columns == 1))>
		constexpr operator U () const {return U(this->at(0, 0));}
	};
	
	// ----------------------- type trait -----------------------------
	template<class T, size_t rows, size_t columns> struct value_type<MatrixRef<T, rows, columns>>{using type = T;};
	
}// namespace twmath