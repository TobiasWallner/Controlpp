/*
	Author: Tobias Wallner
	tobias.wallner1@gmx.net
	
*/

#include <fstream>
#include <iostream>
#include <algorithm>
#include <utility>
#include <sstream>
#include <cstring>
#include <vector>

#define TEST_CASE(function)										\
	if(function()){ 											\
		std::cout << "[  Ok  ] - " << #function << std::endl;	\
	}else{ 														\
		std::cout << "[Failed] - " << #function << std::endl;	\
	}

#include "Vector.hpp"
#include "Polynom.hpp"
#include "TransferFunction.hpp"
#include "ContinuousTransferFunction.hpp"
#include "DiscreteTransferFunction.hpp"

/* TODO:	* Put 'TransferFunction' into its own header
			* add 'diff' / 'differentiate' and 'integ' / 'integrate' functions
			
			* (Optional: add a dynamic vector that is made per default, with 'N=-1')
				* (Optional: specialise all functions for the dynamic vector or use .size() instead of the template parameter 'N')
				
			* start working on the 'StateSpace' type
			
*/

#include "Matrix.hpp"
#include "MatrixFunctions.hpp"

#define result_type_bop(Tl, op, Tr) decltype(std::declval<Tl>() op std::declval<Tr>())

using namespace twmath;

template<class T, size_t internal_states, size_t inputs, size_t outputs>
struct StateSpace{
	Matrix<T, internal_states, internal_states> _A;
	Matrix<T, internal_states, inputs> _B;
	Matrix<T, outputs, internal_states> _C;
	Matrix<T, outputs, inputs> _D;
	
	constexpr StateSpace() = default;
	
	constexpr static size_t number_of_states = internal_states;
	constexpr static size_t number_of_inputs = inputs;
	constexpr static size_t number_of_outputs = outputs;
	
	constexpr Matrix<T, internal_states, internal_states>& A(){return this->_A;}
	constexpr const Matrix<T, internal_states, internal_states>& A()const{return this->_A;}
	
	constexpr Matrix<T, internal_states, inputs>& B(){return this->_B;}
	constexpr const Matrix<T, internal_states, inputs>& B() const {return this->_B;}
	
	constexpr Matrix<T, outputs, internal_states>& C(){return this->_C;}
	constexpr const Matrix<T, outputs, internal_states>& C() const {return this->_C;}
	
	constexpr Matrix<T, outputs, inputs> D() {return this->_D;}
	constexpr const Matrix<T, outputs, inputs> D() const {return this->_D;}
};

template<class T, size_t num_size, size_t den_size>
constexpr StateSpace<T, den_size-1, 1, 1> state_space(const TransferFunction<T, num_size, den_size>& rp){
	StateSpace<T, den_size-1, 1, 1> result;
	const T a_n = rp.den[rp.den.order()];
	const Polynom<T, den_size> a = -(rp.den/a_n);
	const Polynom<T, num_size> b = rp.num / a_n;
	
	result.A().template sub_matrix<result.A().rows()-1, result.A().columns()-1>(0, 1).set_unity();
	result.A().template column<result.A().rows()-1>(0).set_zero();
	result.A().row(result.A().rows()-1) = a.vec.template sub_vector<a.vec.size()-1>();
	
	result.B().template column<result.B().size()-1>(0).set_zero();
	result.B().at(-1, 0) = T(1);
	
	result.C().template row<b.size()>(0) = b.vec;
	result.C().template row<result.C().columns() - b.size()>(0, b.size()).set_zero();
	
	result.D().at(0, 0) = (b.size() > (den_size-1)) ? b[den_size-1] : T(0);
	
	return result;
}

template<class Stream, class T, size_t states, size_t inputs = 1, size_t outputs = 1>
Stream& print_pretty(Stream& stream, const StateSpace<T, states, inputs, outputs>& lsys, const char* indentation = ""){
	size_t max_elem_sizeAC = 0;
	if(lsys.A().rows() > 1) for(size_t row=0; row<lsys.A().rows(); ++row) for(size_t col=0; col<lsys.A().columns(); ++col){
		auto elem = lsys.A().at(row, col);
		std::stringstream s;
		s << elem;
		max_elem_sizeAC = (s.str().size() > max_elem_sizeAC) ? s.str().size() : max_elem_sizeAC;
	}
	if(lsys.C().rows() > 1) for(size_t row=0; row<lsys.C().rows(); ++row) for(size_t col=0; col<lsys.C().columns(); ++col){
		auto elem = lsys.C().at(row, col);
		std::stringstream s;
		s << elem;
		max_elem_sizeAC = (s.str().size() > max_elem_sizeAC) ? s.str().size() : max_elem_sizeAC;
	}
	
	size_t max_elem_sizeBD = 0;
	if(lsys.B().rows() > 1) for(size_t row=0; row<lsys.B().rows(); ++row) for(size_t col=0; col<lsys.B().columns(); ++col){
		auto elem = lsys.B().at(row, col);
		std::stringstream s;
		s << elem;
		max_elem_sizeBD = (s.str().size() > max_elem_sizeBD) ? s.str().size() : max_elem_sizeBD;
	}
	if(lsys.D().rows() > 1) for(size_t row=0; row<lsys.D().rows(); ++row) for(size_t col=0; col<lsys.D().columns(); ++col){
		auto elem = lsys.D().at(row, col);
		std::stringstream s;
		s << elem;
		max_elem_sizeBD = (s.str().size() > max_elem_sizeBD) ? s.str().size() : max_elem_sizeBD;
	}
	// ---------------- Print A and B matrix ----------------
	for(size_t row = 0; row < lsys.A().rows(); ++row){
		const char* x_dash;
		const char* mul_x;
		const char* mul_u;
		if(row == lsys.A().rows()/2){
			x_dash = "x' = ";
			mul_x = " * x + ";
			mul_u = " * u";
		}else{
			x_dash = "     ";
			mul_x = "       ";
			mul_u = "    ";
		}
		
		char open_bracket; 
		char closed_bracket;
		if(lsys.A().rows() == 1){
			open_bracket = '(';
			closed_bracket = ')';
		}else if(row == 0){
			open_bracket = '/';
			closed_bracket = '\\';
		}else if (row == lsys.A().rows()-1){
			open_bracket = '\\';
			closed_bracket = '/';
		}else{
			open_bracket = '|';
			closed_bracket = '|';
		}
		
		// A Matrix
		stream << indentation << x_dash << open_bracket;
		for(size_t column = 0; column < lsys.A().columns(); ++column){
			std::stringstream s;
			if(column!=0) stream << ' ';
			s << lsys.A().at(row, column);
			if(lsys.A().rows() > 1) for(size_t i = s.str().size(); i < max_elem_sizeAC; ++i) stream << ' ';
			stream << s.str();
		}
		stream << closed_bracket << mul_x;
		
		// B Matrix
		stream << open_bracket;
		for(size_t column = 0; column < lsys.B().columns(); ++column){
			std::stringstream s;
			if(column!=0) stream << ' ';
			s << lsys.B().at(row, column);
			if(lsys.B().rows() > 1) for(size_t i = s.str().size(); i < max_elem_sizeBD; ++i) stream << ' ';
			stream << s.str();
		}
		stream << closed_bracket;
		stream << mul_u << '\n';
	}
	
	// ---------------- Print C and D matrix ----------------
	for(size_t row = 0; row < lsys.C().rows(); ++row){
		
		const char* x_dash;
		const char* mul_x;
		const char* mul_u;
		if(row == lsys.C().rows()/2){
			x_dash = " y = ";
			mul_x = " * x + ";
			mul_u = " * u";
		}else{
			x_dash = "     ";
			mul_x = "       ";
			mul_u = "    ";
		}
		
		char open_bracket; 
		char closed_bracket;
		if(lsys.C().rows() == 1){
			open_bracket = '(';
			closed_bracket = ')';
		}else if(row == 0){
			open_bracket = '/';
			closed_bracket = '\\';
		}else if (row == lsys.C().rows()-1){
			open_bracket = '\\';
			closed_bracket = '/';
		}else{
			open_bracket = '|';
			closed_bracket = '|';
		}
		
		// C Matrix
		stream << indentation << x_dash << open_bracket;
		for(size_t column = 0; column < lsys.C().columns(); ++column){
			std::stringstream s;
			if(column!=0) stream << ' ';
			s << lsys.C().at(row, column);
			for(size_t i = s.str().size(); i < max_elem_sizeAC; ++i) stream << ' ';
			stream << s.str();
		}
		stream << closed_bracket << mul_x;
		
		// D Matrix
		stream << open_bracket;
		for(size_t column = 0; column < lsys.D().columns(); ++column){
			std::stringstream s;
			if(column!=0) stream << ' ';
			s << lsys.D().at(row, column);
			for(size_t i = s.str().size(); i < max_elem_sizeBD; ++i) stream << ' ';
			stream << s.str();
		}
		stream << closed_bracket;
		stream << mul_u << '\n';
	}
	return stream;
}

template<class T, size_t internal_states, size_t inputs = 1, size_t outputs = 1>
struct ContinuousStateSpace{
	StateSpace<T, internal_states, inputs, outputs> ls;
	
	constexpr ContinuousStateSpace() = default;
	constexpr ContinuousStateSpace(const ContinuousStateSpace&) = default;
	
	constexpr ContinuousStateSpace(const StateSpace<T, internal_states, inputs, outputs>& ls) : ls(ls){};
	
	constexpr Matrix<T, internal_states, internal_states>& A(){return this->ls.A();}
	constexpr const Matrix<T, internal_states, internal_states>& A()const{return this->ls.A();}
	
	constexpr Matrix<T, internal_states, inputs>& B(){return this->ls.B();}
	constexpr const Matrix<T, internal_states, inputs>& B() const {return this->ls.B();}
	
	constexpr Matrix<T, outputs, internal_states>& C(){return this->ls.C();}
	constexpr const Matrix<T, outputs, internal_states>& C() const {return this->ls.C();}
	
	constexpr Matrix<T, outputs, inputs> D() {return this->ls.D();}
	constexpr const Matrix<T, outputs, inputs> D() const {return this->ls.D();}
};

template<class T, size_t num_size, size_t den_size>
static constexpr ContinuousStateSpace<T, den_size-1, 1, 1> c_state_space (const ContinuousTransferFunction<T, num_size, den_size>& cls){
	const StateSpace<T, den_size-1, 1, 1> ss = state_space(cls.rp);
	ContinuousStateSpace<T, den_size-1, 1, 1> result(ss);
	return result;
}


template<class Stream, class T, size_t states, size_t inputs = 1, size_t outputs = 1>
Stream& print_pretty(Stream& stream, const ContinuousStateSpace<T, states, inputs, outputs>& cls, const char* indentation = ""){
	return print_pretty(stream, cls.ls, indentation);
}

template<class T, size_t internal_states, size_t inputs = 1, size_t outputs = 1>
struct DiscreteStateSpace{
	StateSpace<T, internal_states, inputs, outputs> ls;
	T ts;
	
	constexpr DiscreteStateSpace() = default;
	
	template<size_t num_size>
	constexpr DiscreteStateSpace(const DiscreteTransferFunction<T, num_size, internal_states>& dls) : ls(dls.rp){}
	
	template<size_t num_size>
	constexpr DiscreteStateSpace(const StateSpace<T, internal_states, inputs, outputs>& ls, const T& ts) : ls(ls), ts(ts){}
	
	constexpr Matrix<T, internal_states, internal_states>& A(){return this->ls.A();}
	constexpr const Matrix<T, internal_states, internal_states>& A()const{return this->ls.A();}
	
	constexpr Matrix<T, internal_states, inputs>& B(){return this->ls.B();}
	constexpr const Matrix<T, internal_states, inputs>& B() const {return this->ls.B();}
	
	constexpr Matrix<T, outputs, internal_states>& C(){return this->ls.C();}
	constexpr const Matrix<T, outputs, internal_states>& C() const {return this->ls.C();}
	
	constexpr Matrix<T, outputs, inputs> D() {return this->ls.D();}
	constexpr const Matrix<T, outputs, inputs> D() const {return this->ls.D();}
};

template<class Tm, class Ts, size_t states>
constexpr DiscreteStateSpace<result_type_bop(Tm, *, Ts), states, 1, 1> discretise(const ContinuousStateSpace<Tm, states, 1, 1>& sys, const Ts& sample_time){
	using ValueType = result_type_bop(Tm, *, Ts);
	DiscreteStateSpace<ValueType, states, 1, 1> result;
	Matrix<ValueType, states+1, states+1> M;

	// preparation
	M.template sub_matrix<states, states>(0, 0) = sys.A();
	M.template column<states>(states) = sys.B().column();
	M.row(states).set_zero();
	
	// calculation
	M = mexp(M * sample_time);
	
	// re-assignment
	result.A() = M.template sub_matrix<states, states>(0, 0);
	result.B().column() = M.template column<states>(states);
	result.C() = sys.C();
	result.D() = sys.D();
	result.ts = sample_time;
	
	return result;
}

template<class Stream, class T, size_t states, size_t inputs = 1, size_t outputs = 1>
Stream& print_pretty(Stream& stream, const DiscreteStateSpace<T, states, inputs, outputs>& dsys, const char* indentation = ""){
	print_pretty(stream, dsys.ls, indentation);
	return stream << indentation << "Ts = " << dsys.ts << '\n';
}

template<class T, size_t internal_states, size_t inputs = 1, size_t outputs = 1>
struct DiscreteLinearSystem{
	DiscreteStateSpace<T, internal_states, inputs, outputs> sys;
	Vector<T, internal_states> x;
	
	Vector<T, outputs> input(const Vector<T, inputs>& u){
		auto y = this->sys.C * this->x + this->sys.D * u;
		this->x = this->sys.A * this->x + this->sys.B * u;
		return y;
	}
	
	void reset() {this->x.set_zero();}
	
	friend Vector<T, outputs> operator >> (const Vector<T, inputs>& u, const DiscreteLinearSystem& df){
		return df.input(u);
	}
};

template<class Ts, class Tv>
class TimeSeries{
	std::vector<Ts> times; // todo: change to my vector once I support dynamic vectors
	std::vector<Tv> values; // todo: change to my vector once I support dynamic vectors

public:

	void push_back(const Ts& time, const Tv& value){
		this->times.push_back(time);
		this->values.push_back(value);
	}
	
	void reserve(size_t n){
		times.reserve(n);
		values.reserve(n);
	}
	
	Ts& time(size_t i) {return this->times[i];}
	Tv& value(size_t i) {return this->values[i];}
	
	const Ts& time(size_t i) const {return this->times[i];}
	const Tv& value(size_t i)  const {return this->values[i];}
	
	size_t size() const {return this->times.size();}
};

template<class Stream, class Ts, class Tv>
Stream& operator << (Stream& stream, const TimeSeries<Ts, Tv>& timeSeries){
	for(size_t i = 0; i < timeSeries.size(); ++i){
		stream << timeSeries.time(i) << ", " << timeSeries.value(i) << '\n';
	}
	return stream; 
}

template<class Ts, class Tv, size_t states, size_t inputs = 1, size_t outputs = 1>
TimeSeries<Ts, Vector<Tv, outputs>> step(const DiscreteLinearSystem<Tv, states, inputs, outputs>& filter, const Ts& simulation_time){
	filter.reset();
	TimeSeries<Ts, Vector<Tv, outputs>> ts;
	ts.reserve(simulation_time / filter.ts + 1);
	const Tv u(1);
	for(Tv t = 0; t <= simulation_time; t += filter.ts){
		ts.times.push_back(t);
		ts.values.push_back(filter.input(u));
	}
	return ts;
}

template<class T, size_t states>
TimeSeries<double, T> step(const DiscreteStateSpace<T, states, 1, 1>& sys, double simulation_time){
	TimeSeries<double, T> ts;
	ts.reserve(simulation_time / sys.ts + 1);
	
	const T u(1);
	Vector<T, states> x(T(0));
	
	for(double t = 0; t <= simulation_time; t += sys.ts){
		const auto v1 = sys.C() * x;
		const auto v2 = sys.D() * u;
		const auto value = v1 + v2; 
		
		const auto v3 = sys.B() * u;
		const auto v4 = sys.A() * x;
		x = v3.column() + v4;
		
		ts.push_back(t, value);
	}
	return ts;
}



#undef result_type_bop

int main(){
	
	
	std::cout << "<Project Name> Tests:" << std::endl;
	std::cout << "---------------" << std::endl;
	
	const auto x = Polynom({0, 1});
	
	std::cout << "x = " << x << std::endl;
	
	auto g = 1 + x + x * x * 2;
	auto h = 3 + 2*x;
	
	std::cout << "g = " << g << std::endl;
	std::cout << "h = " << h << std::endl;
	
	std::cout << "g / h = " << (g / h) << std::endl;
	std::cout << "g % h = " << (g % h) << '\n' << std::endl;
	
	std::cout << "h * (g / h) + (g % h) = " << (h * (g / h) + (g % h)) << std::endl;
	
	const ContinuousTransferFunction<float, 2, 1> s({0.f, 1.f}, 1.f);
	
	auto G = (1.f) / (1 + s + s*s);
	
	print_pretty(std::cout, G, "G", "    ") << std::endl;
		
	auto Gsys = c_state_space(G); // have a look at the conversion of the transfer function into state space ... seems to have wrong number of dimensions (one to many)
	
	print_pretty(std::cout, Gsys, "    ") << std::endl;
	
	const auto expG = mexp(Gsys.A());
	
	print_pretty(std::cout, expG, "mexp(G.A)", "    ") << std::endl;
	
	
	const auto Gsys_d = discretise(Gsys, 0.01);
	
	print_pretty(std::cout, Gsys_d, "    ") << std::endl;
	
	const auto ts = step(Gsys_d, 20);
	
	std::ofstream file("step.csv");
	file << "t, value" << std::endl;
	file << ts << std::endl;
	
	
	Matrix<float, 3, 4> M({
		{1, 2, 3, 4},
		{5, 6, 7, 8},
		{9, 10, 11, 12},
		//{13, 14, 15, 16}
	});
	
	std::cout << "M:\n" << M << std::endl;
	
	std::cout << "transpose(M):\n" << transpose(M) << std::endl;
	
	return 0;
	
}