#pragma once
/*

	TODO: implicit conversion between matrix<n, 1> and matrix<1, n> to vector.
			- vector types are just a specialisatio of the matrix type?
			- duplicate matrix class and add cast?
			- add constructor for vector?

*/

/*
	To disable boundary checks:
		Define: DISABLE_TWMATH_BOUNDARY_CHECKS
		
	To handle errors, define one of the following
		TWMATH_THROW 				(default)
		TWMATH_EXIT 					(calls exit(-1) to exit the program)
		TWMATH_TRAP 					(traps the execution in a while loop)
		TWMATH_CUSTOM_ERROR_HANDLER 	(calls the user defined function: 'void twmath_custom_error_handler(const char* error_message)')
		
	Set the error output stream by defining (only for TWMATH_EXIT, TWMATH_TRAP or TWMATH_RETURN_ZERO):
		TWMATH_CERR	(default is std::cerr)
*/

#define TWMATH_ENABLE_IF(condition) typename std::enable_if_t<(condition), int> = 0

#ifdef DISABLE_TWMATH_BOUNDARY_CHECKS
	#define twmath_assert(condition, message)
#else
	
	#ifndef TWMATH_CERR
		#define TWMATH_CERR std::cout
	#endif

	#if defined(TWMATH_CUSTOM_ERROR_HANDLER)
		void twmath_custom_error_handler(const char* error_message);
		#define twmath_assert(condition, message) if(!(condition)){std::stringstream s; s << message; twmath_custom_error_handler(s.str().c_str());}
	#elif defined(TWMATH_TRAP)
		#define twmath_assert(condition, message) if(!(condition)){TWMATH_CERR << message; while(true){};}
	#elif defined(TWMATH_EXIT)
		#define twmath_assert(condition, message) if(!(condition)){TWMATH_CERR << message; exit(-1);}
	#else
		#include <exception>
		#define twmath_assert(condition, message) if(!(condition)){std::stringstream s; s << message; throw std::runtime_error(s.str().c_str());}
	#endif

#endif

#define result_type2(function, Tl, Tr) decltype(function(std::declval<Tl>(), std::declval<Tr>()))