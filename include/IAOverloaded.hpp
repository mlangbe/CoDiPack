/*
 * Adaptation class 
 * providing classes/methods which have special overloads for intervals.
 *
 *
 * Copyright (C) 2018-2019
 * RHRK - Regional Computing Centre, Kaiserslautern,
 * &
 * Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 *
 * This file is free software: you can redistribute it and/or
 * modify it under the terms of the Apache License Version 2
 * as published by the Apache Software Foundation.
 *
 * it is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * You should have received a copy of the
 * Apache License 2 together with this file.
 * If not, see <https://www.apache.org/licenses/>.
 *
 * Author: Max Langbein (RHRK TU Kaiserslautern)
*/

#ifndef IA_OVERLOADED_HPP
#define IA_OVERLOADED_HPP
#pragma once


//for abs:
#include<cmath>

//for std::max/min
#include<algorithm>

#include<cassert>

/**
* Adaptations
* providing classes/methods which have special overloads for intervals.
* utilities which behave normal for non-interval types, but need overloads for intervals. 
* mostly the boolean types.
* here the "normal" behavior is given.
* <br>
* general contract for the overloads is that f( [x,x] ) = f(x), f( [x,x], [y,y] ) = f(x,y) .
* <br>
* for those functions that may return an interval (here the sqr, sqr3d and firstif functions)
* the following conditions for the overloads must hold (as for all ia function versions):
* <br>
* and  f(x) \in  f( [xa,xb] ) \forall x \in [xa,xb],
* <br>
* and  f(x,y) \in  f( [xa,xb],[ya,yb] ) \forall (x,y) \in [xa,xb] \times [ya,yb]
* <br>
*
* The overloading is done by providing specializations
* of the structs defined in 
* IAOverloaded for the interval types
* for which defaults are given below.
* <br>
* The calling is done using the functions defined below
* having the lowercase name of the structs.
*
* 
*/
namespace ia{


/**
struct containing all overloadable structs
*/
struct IAOverloaded
{	

	/* 
	struct containing all functions returning the base value of the interval type,
	* that is upper, lower, mid, and minAbs.
	* the codi versions will also return the base value, and discard
	* the derivative information
	*/
	template<typename T>
	struct baseValueFunctions;

	
	/*------------------------------------------------------------------------*/

	/* functions converting bool intervals to bool (isCompletelyTrue and isCompletelyFalse)*/
	template<typename BoolType>
	struct iabool;	


	
	/*------------------------------------------------------------------------*/
	/* single-function template classes used to overload the behavior of      */
	/* the function of the same name in lowercase.                            */
	/*------------------------------------------------------------------------*/
	
	/* x^2. */
	template<typename T,typename enable=void>
	struct Sqr;

	/* upper(x)==upper(y)&&lower(x)==lower(y) */
	template<typename T,typename S=T>
	struct SameInterval;

	/* x?a:b */
	template<typename T,typename BoolType>
	struct Firstif;	

	/* join two values to form an interval, or join two intervals */
	template<typename FirstArgument,typename SecondArgument=FirstArgument>
	struct Join;	

	/* claculate the intersection of two intervals */
	template<typename FirstArgument,typename SecondArgument=FirstArgument>
	struct Intersect;	
		
	/** calculate the 2-norm (=length) of a 3d vector of type T*/
	template<typename T,int n=3>
	struct Norm2;

	/** normalize a 3d vector*/
	template<typename T,int n=3>
	struct Normalize;

	/** calculate x * || x ||*/
	template<typename T,int n=3>
	struct MultiplyByNorm;

	/** overloadable min and max functions.*/
	template<typename T,typename S,typename enable=void>
	struct MaxMin;

	
	/** evaluate Unary function with an ia counterpart.
	*Fun is an object which may have alternatives for intervals, and derivatives etc
	*/
	template<typename Fun, typename T,typename ReturnType=T>
	struct EvalUnary;

	/** evaluate Binary function with an ia counterpart.
	*Fun is an object which may have alternatives for intervals, and derivatives etc
	*/
	template<typename Fun, typename FirstArgument,typename SecondArgument=FirstArgument,typename ReturnType=FirstArgument>
	struct EvalBinary;

	/** evaluate vector-valued function
	*Fun is an object which may have alternatives for intervals, and derivatives etc
	*/
	template<typename Fun, typename T,int n=3,typename ReturnType=T>
	struct EvalVector;

};

template<typename Fun,typename T,typename ReturnType>
struct IAOverloaded::EvalUnary{
	/** 
	 * @brief evaluate f(x)
	 * @param[in] x a number
	 * @return value corresponding to f(x) 
	*/
	inline static void _eval(ReturnType&y,const Fun&f,const T&x)
	{
		y=f(x);
	}	
};


template<typename Fun,typename T,int n,typename ReturnType>
struct IAOverloaded::EvalVector{
	
	/** 
	 * @brief evaluate f(x)
	 * @param[in] x a number
	 * @return value corresponding to f(x) 
	*/
	inline static void _eval(ReturnType&y,const Fun&f,const T (&x)[n])
	{
		y=f(x);
	}	
};

template<typename Fun, typename FirstArgument,typename SecondArgument,typename ReturnType>
struct IAOverloaded::EvalBinary{

	/** 
	 * @brief evaluate f(x,y)
	 * @param[in] f a function object
	 * @param[in] x a number
	 * @param[in] y a number
	 * @return value corresponding to f(x) 
	*/
	inline static void _eval(ReturnType&y,const Fun&f,const FirstArgument&a,const SecondArgument&b )
	{
		y=f(a,b);
	}	
};

/** 
 * @brief evaluate f(x),
 * possibly using alternative implementations 
 * for derivatives and intervals.
 * @param[in] x a number
 * @return value corresponding to x*x  . the interval version returns value corresp. to  [ lower(abs(x))^2, upper(abs(x))^2 ].
*/
template<typename Fun,typename T,typename RetType>
inline 
void eval(RetType&y,const Fun&f, const T&x)
{
	IAOverloaded::EvalUnary<Fun,T,RetType>::_eval(y,f,x);
}

/** 
 * @brief evaluate f(x),
 * possibly using alternative implementations 
 * for derivatives and intervals for a certain vector length
 * @param[in] x a number
 * @return value corresponding to x*x  . the interval version returns value corresp. to  [ lower(abs(x))^2, upper(abs(x))^2 ].
*/
template<typename Fun,typename T,int n,typename RetType>
inline
void eval(RetType&y,const Fun&f, const T (&x)[n] )
{
	IAOverloaded::EvalVector<Fun,T,n,RetType>::_eval(y,f,x);
}

/** 
 * @brief evaluate f(x,y),
 * possibly using alternative implementations 
 * for derivatives and intervals.
 * @param[in] x a number
 * @param[in] y a number
 * @return value corresponding to f(x)
*/
template<typename Fun,typename FirstArgument,typename SecondArgument,typename RetType>
inline 
void eval(RetType&ret, const Fun&f, const FirstArgument&x,const SecondArgument&y)
{
	IAOverloaded::EvalBinary<Fun,FirstArgument,SecondArgument,RetType>::_eval(ret,f,x,y);
}


template<typename T,typename enable>
struct IAOverloaded::Sqr{

	/** 
	 * @brief the square of x . The interval overload should return the smallest interval containing all squares of the interval in x.
	 * @param[in] x a number
	 * @return value corresponding to x*x  . the interval version returns value corresp. to  [ lower(abs(x))^2, upper(abs(x))^2 ].
	*/
	inline static T _sqr(const T&x)
	{
		return x*x;
	}	
};


/** 
 * @brief the square of x . The interval overload should return the smallest interval containing all squares of the interval in x.
 * @param[in] x a number
 * @return value corresponding to x*x  . the interval version returns value corresp. to  [ lower(abs(x))^2, upper(abs(x))^2 ].
*/
template<typename T>
inline T sqr(const T&x)
{
	return IAOverloaded::Sqr<T>::_sqr(x);
}



/**
* @tparam T : any type 
* @tparam BoolType : bool or Interval<bool>
*/
template<typename T,typename BoolType>
struct IAOverloaded::Firstif {

	using ReturnType = const T&;

	/** 
	 * @brief return second arg if first is true,  otherwise third. 
	 * corresponds to operator ?: which sadly cannot be IAOverloaded.
	 * the interval overload (with x being an interval )
	 * will return join(a,b) if x==[false,true]
	 *
	 * @param[in] x boolean value
	 * @param[in] a first T argument
	 * @param[in] b 2nd T argument
	 * @return : a if x, else b
	 */
	inline static ReturnType _firstif (BoolType x,const T&a,const T&b )
	{
		return x ? a : b;
	};
};

/** 
 * @brief return second arg if first is true,  otherwise third. 
 * corresponds to operator ?: which sadly cannot be IAOverloaded.
 * the interval overload (with x being an interval )
 * will return join(a,b) if x==[false,true]
 *
 * @param[in] x boolean value
 * @param[in] a first T argument
 * @param[in] b 2nd T argument
 * @return : a if x, else b
 * @tparam T : anything
 * @tparam BoolType : bool or Interval<bool>
*/
template<typename T,typename BoolType>
inline typename IAOverloaded::Firstif<T,BoolType>::ReturnType firstif (const BoolType&x, const T&a, const T&b )
{
	return IAOverloaded::Firstif<T,BoolType>::_firstif(x,a,b);
}


/**
 * @tparam T type of 1st argument 
 * @tparam S type of 2nd argument 
*/
template <typename FirstArgument,typename SecondArgument>
struct IAOverloaded::Join
{
	using ReturnType = FirstArgument;
	
	/**
	 * @brief returns NaN if a!=b, or a  else.
	 * The overload for an interval should return the minimum-sized interval which contains a and b. 
	 * @tparam T number type
	 * @param[in] a first interval
	 * @param[in] b 2nd interval
	 * @return a==b ? a:  NaN . The overload for an interval should return the minimum-sized interval containing a and b if a and b are intervals.
	*/
	inline static ReturnType _join(FirstArgument a, const SecondArgument&b )
	{
		//if(a==b)
		//	return a;

		(void)a;
		(void)b;

		assert( false && "join can only be used on intervals. perhaps a NaN value lead to this path ?" );
		
		return ReturnType();
		//return NaN
		//return (ReturnType)(a/0.0) ;
	}
};

template <typename FirstArgument,typename SecondArgument>
inline typename IAOverloaded::Join<FirstArgument,SecondArgument>::ReturnType join(const FirstArgument& a, const SecondArgument&b )
{
	return IAOverloaded::Join<FirstArgument,SecondArgument>::_join(a,b);
}

/**
 * @tparam T type of 1st argument 
 * @tparam S type of 2nd argument 
*/
template <typename FirstArgument,typename SecondArgument>
struct IAOverloaded::Intersect
{
	using ReturnType = FirstArgument;
	
	/**
	 * @brief returns NaN if a!=b, or a  else.
	 * The overload for an interval should return the minimum-sized interval which contains a and b. 
	 * @tparam T number type
	 * @param[in] a first interval
	 * @param[in] b 2nd interval
	 * @return a==b ? a:  NaN . The overload for an interval should return the minimum-sized interval containing a and b if a and b are intervals.
	*/
	inline static ReturnType _intersect(FirstArgument a, const SecondArgument&b )
	{
		if(a==b)
			return a;

		//assert( ("join can only be used on intervals. perhaps a NaN value lead to this path ?",false) );	
		
		//return NaN
		return (ReturnType)(a/0.0) ;
	}
};

template <typename FirstArgument,typename SecondArgument>
inline typename IAOverloaded::Intersect<FirstArgument,SecondArgument>::ReturnType intersect(const FirstArgument& a, const SecondArgument&b )
{
	return IAOverloaded::Intersect<FirstArgument,SecondArgument>::_intersect(a,b);
}


/**
* class to be specialized for different types
*( currently, there are codi::ActiveReal, and ia::Interval ) 
* which contains methods returning values of the floating point type of T
*
* the codi overload may return the value of the basevalue.
*
*/
template<class T>
struct IAOverloaded::baseValueFunctions{

	/** this should be a basic Float Type. for Interval<double> e.g. double*/
	using ReturnType = T;
	
	/**
	 * @brief returns x. The overload for an interval should return the upper bound of the interval.
	 * the overload for codi should return the upper bound of the primal value.
	 * <br>
	 * Take care that you do not use the result in expressions which which be used for deruvatuive computations,
	 * but only use them for comparisons.
	 * 
	 * @param[in] x interval or number
	 * @return upper bound as number type
	 * @tparam T a number or interval type
	*/
	static const ReturnType& _upper(const T &x)
	{
		return x;
	}
	

	/**
	 * @brief returns x.  The overload for an interval should return the lower bound of the interval.
	 * the overload for codi should return the lower bound of the primal value.
	 * <br>
	 * Take care that you do not use the result in expressions which which be used for deruvatuive computations,
	 * but only use them for comparisons.
	 * 
	 * @param[in] x interval or number
	 * @return lower bound as number type
	 * @tparam T a number or interval type
	*/
	static const ReturnType& _lower(const T &x)
	{
		return x;
	}

	/**
	 * @brief returns x.  The overload for an interval should return the center of the interval
	 * the overload for codi should return the center of thze interval of the primal value.
	 * <br>
	 * Take care that you do not use the result in expressions which which be used for deruvatuive computations,
	 * but only use them for comparisons.
	 * 
	 * @param[in] x interval or number
	 * @return lower bound as number type
	 * @tparam T a number or interval type
	*/
	static const ReturnType& _mid(const T &x)
	{
		return x;
	}
	
	/**
	* @brief returns abs(x). The overload for an interval should  return the lower absolute value, if x does not contain 0, and 0 elsewise
	* @param[in] x value or interval
	* @return value equal to lower(abs(x))
	* @tparam T number  or interval of numbers
	*/
	static ReturnType _minAbs(const T&x)
	{
		using std::abs;
		return abs(x);
	}


};

/**
 * @brief returns x. The overload for an interval should return the upper bound of the interval.
 * the overload for codi should return the upper bound of the primal value.
 * <br>
 * Take care that you do not use the result in expressions which which be used for deruvatuive computations,
 * but only use them for comparisons.
 * 
 * @param[in] x interval or number
 * @return upper bound as number type
 * @tparam T a number or interval type
*/
template<class T>
inline typename IAOverloaded::baseValueFunctions<T>::ReturnType upper(const T &x)
{
	return IAOverloaded::baseValueFunctions<T>::_upper(x);
}

/**
 * @brief returns x.  The overload for an interval should return the lower bound of the interval.
 * the overload for codi should return the lower bound of the primal value.
 * <br>
 * Take care that you do not use the result in expressions which which be used for deruvatuive computations,
 * but only use them for comparisons.
 * 
 * @param[in] x interval or number
 * @return lower bound as number type
 * @tparam T a number or interval type
*/
template<class T>
inline typename IAOverloaded::baseValueFunctions<T>::ReturnType lower(const T &x)
{
	return IAOverloaded::baseValueFunctions<T>::_lower(x);
}

/**
 * @brief returns x.  The overload for an interval should return the center of the interval
 * the overload for codi should return the center of the interval of the primal value.
 * <br>
 * Take care that you do not use the result in expressions which which be used for deruvatuive computations,
 * but only use them for comparisons.
 * 
 * @param[in] x interval or number
 * @return lower bound as number type
 * @tparam T a number or interval type
*/
template<class T>
inline typename IAOverloaded::baseValueFunctions<T>::ReturnType mid(const T &x)
{
	return IAOverloaded::baseValueFunctions<T>::_mid(x);
}

/**
 * @brief returns abs(x). The overload for an interval should  return the lower absolute value, if x does not contain 0, and 0 elsewise
 * @param[in] x value or interval
 * @return value equal to lower(abs(x))
 * @tparam T number  or interval of numbers
*/
template<class T> 
inline typename IAOverloaded::baseValueFunctions<T>::ReturnType minAbs(const T&x)
{
	return IAOverloaded::baseValueFunctions<T>::_minAbs(x);
}



/**
* @brief return true if the two values represent exactly the same interval.
* for non-intervals, just the == operator
* @param[in] x first interval or number
* @param[in] x second interval or  number
* @return value equal to  (upper(x)==upper(y) && lower(x)==lower(y))
* @tparam T a number or interval type
* @tparam S a number or interval type
*/
template<class T,class S>
struct IAOverloaded::SameInterval
{
	/**
	* @brief return true if the two values represent exactly the same interval. 
	* for non-intervals, just the == operator
	* @param[in] x first interval or number
	* @param[in] x second interval or  number
	* @return value equal to  (upper(x)==upper(y) && lower(x)==lower(y))
	* @tparam T a number or interval type
	*/
	inline static bool _sameInterval(const T&x,const S&y)
	{
		return x == y;
	}	
};







/**
* @brief return true if the two values represent exactly the same interval. 
* for non-intervals, just the == operator
* @param[in] x first interval or number
* @param[in] x second interval or  number
* @return value equal to  (upper(x)==upper(y) && lower(x)==lower(y))
* @tparam T a number or interval type
*/
template<class T,class S>
inline bool sameInterval(const T&x,const S&y)
{
	return IAOverloaded::SameInterval<T,S>::_sameInterval(x,y);
}


template<class BoolType>
struct IAOverloaded::iabool
{

	/**
	 * @brief returns !x. the overload for an interval should return true if bool interval only contains false
	 * @param[in] bool or interval of bool
	 * @return value equal to !upper(x)
	 * @tparam T bool or myinterval<bool>
	*/
	inline static bool _isCompletelyFalse(BoolType x)
	{
		return !x;
	}

	/**
	 * @brief returns x. the overload for an interval should returns true if bool interval only contains true
	 * @param[in] bool or interval of bool
	 * @return value eual to lower(x)
	 * @tparam T bool or myinterval<bool>
	*/
	inline static bool _isCompletelyTrue(BoolType x)
	{
		return x;
	}

};



/**
 * @brief returns !x. The overload for an interval should
 * return true if the bool interval only contains false.
 * @param[in] bool or interval of bool
 * @return value eual to lower(x)
 * @tparam T bool or myinterval<bool>
*/
template<typename BoolType>
inline bool isCompletelyFalse(const BoolType& x)
{
	return IAOverloaded::iabool<BoolType>::_isCompletelyFalse(x);
}

/**
 * @brief returns x. The overload for an interval should returns true if bool interval only contains true
 * @param[in] bool or interval of bool
 * @return value eual to lower(x)
 * @tparam T bool or myinterval<bool>
*/
template<typename BoolType>
inline bool isCompletelyTrue(const BoolType& x)
{
	return IAOverloaded::iabool<BoolType>::_isCompletelyTrue(x);
}







/**
 * @brief utility template function return the sum of square-roots of the three values (square of euclid.norm) 
 * @param[in] x array with 3 entries
 * @return ||x||^2
*/
template<typename T>
inline T sqr3d(const T x[3])
{
	return sqr(x[0])+sqr(x[1])+sqr(x[2]);
}


template<int n,typename T>
inline T sqrND(const T x[n])
{
	T ret=sqr(x[0]);
	for(int i=1;i<n;++i)
		ret+=sqr(x[i]);
	return ret;
}

template<int n,typename T>
inline T sqr(const T (&x)[n])
{
	return sqrND<n>(x);
}



template<typename T,int n>
struct IAOverloaded::Norm2
{
	/**
	* @brief return the vector's length.
	* It should be overloaded for codi intervals,
	* as the gradient would not be calculable if the box contains zero.
	* @param[in] vec vector with three entries representing a 3d vector
	* @return value equal to ||vec||_2
	* @tparam T a number or interval type
	*/
	inline static T _norm2(const T vec[n])
	{
		return sqrt(sqrND<n>(vec));
	}
};


template<typename T>
struct IAOverloaded::Norm2<T,3>
{
	/**
	* @brief return the vector's length.
	* It should be overloaded for codi intervals,
	* as the gradient would not be calculable if the box contains zero.
	* @param[in] vec vector with three entries representing a 3d vector
	* @return value equal to ||vec||_2
	* @tparam T a number or interval type
	*/
	inline static T _norm2(const T vec[3])
	{
		return sqrt(sqr3d(vec));
	}	
};


/**
* @brief return the vector's length.
* @param[in] vec vector with three entries representing a 3d vector
* @return value equal to ||vec||_2
* @tparam T a number or interval type
*/
template<class T>
inline T norm2(const T x[3])
{
	return IAOverloaded::Norm2<T,3>::_norm2(x);
}

template<int n,class T>
inline T norm2ND(const T x[n])
{
	return IAOverloaded::Norm2<T,n>::_norm2(x);
}

/** normalize an n-d vector*/
template<typename T,int n>
struct IAOverloaded::Normalize
{
	/**
	 * @brief normalize the 3d vector x, (scale so it has length 1).
	 * the interval overload should return a box which contains all
	 * possible directions which stem from the original box.
	 * It should be overloaded because a box containing zero would
	 * return an infinite box otherwise.
	 *
	 * @param[inout] x the 3d vector
	 * @return the vectors original length
	 * @tparam T a number type , may be an interval , may have derivatives
	 */
	inline static T _normalize(T vec[n])
	{
		T len = Norm2<T,n>::_norm2(vec);
		for(int i=0;i<n;++i)
			vec[i]/=len;
		return len;
	}
};


template<typename T>
struct IAOverloaded::Normalize<T,3>
{
	/**
	 * @brief normalize the 3d vector x, (scale so it has length 1).
	 * the interval overload should return a box which contains all
	 * possible directions which stem from the original box.
	 * It should be overloaded because a box containing zero would
	 * return an infinite box otherwise.
	 *
	 * @param[inout] x the 3d vector
	 * @return the vectors original length
	 * @tparam T a number type , may be an interval , may have derivatives
	 */
	inline static T _normalize(T vec[3])
	{
		T len = norm2(vec);
		vec[0]/=len;
		vec[1]/=len;
		vec[2]/=len;
		return len;
	}	
};

/**
 * @brief normalize the 3d vector x, (scale so it has length 1).
 *
 * @param[inout] x the 3d vector
 * @return the vectors original length
 * @tparam T a number type , may be an interval , may have derivatives
 */
template<int n,class T>
inline T normalizeND(T x[n])
{
	return IAOverloaded::Normalize<T,n>::_normalize(x);
}

template<class T>
inline T normalize(T x[3])
{
	return IAOverloaded::Normalize<T,3>::_normalize(x);
}


/** multiply by norm of an  n-d vector. default implementation.*/
template<typename T,int n>
struct IAOverloaded::MultiplyByNorm
{
	/**
	* @brief nultiplyByNorm vector
	* @param[inout] vec vector with three entries represnting a 3d vector,
	* which will be nultiplyByNormd
	* @return the norm
	*/
	inline static T _multiplyByNorm(T vec[n])
	{
		T len = norm2<T,n>(vec);
		for(int i=0;i<n;++i)
			vec[i]*=len;
		return len;
	}
};

/** multiply by norm of a 3d vector*/
template<typename T>
struct IAOverloaded::MultiplyByNorm<T>
{
	/**
	* @brief nultiplyByNorm vector
	* @param[inout] vec vector with three entries represnting a 3d vector,
	* which will be nultiplyByNormd
	* @return the norm
	*/
	inline static T _multiplyByNorm(T vec[3])
	{
		T len = norm2(vec);
		vec[0]*=len;
		vec[1]*=len;
		vec[2]*=len;
		return len;
	}
};


/**
 * @brief multiply the the vector x by it's length
  * @param[inout] x vector with three entries representing a 3d vector,
  * which will be multiplied by it's length
  * @return the length
 * @tparam T a number type , may be an interval , may have derivatives
 */
template<class T>
inline T multiplyByNorm(T x[3])
{
	return IAOverloaded::MultiplyByNorm<T,3>::_multiplyByNorm(x);
}

/**
 * @brief multiply the the vector x by it's length
  * @param[inout] x vector with three entries representing a 3d vector,
  * which will be multiplied by it's length
  * @return the length
 * @tparam T a number type , may be an interval , may have derivatives
 */
template<int n,class T>
inline T multiplyByNormND(T x[n])
{
	return IAOverloaded::MultiplyByNorm<T,n>::_multiplyByNorm(x);
}


template<typename T,typename S,typename enable=void>
struct IAOverloaded::MaxMin
{
	static auto _max(const T&x,const S&y) ->decltype(std::max(x,y))
	{
		using std::max;
		return max(x,y);
	}

	static auto _min(const T&x,const S&y) ->decltype(std::max(x,y))
	{
		using std::min;
		return min(x,y);
	}
};

template<typename T,typename S>
auto max(const T&a, const S&b) -> decltype( IAOverloaded::MaxMin<T,S>::_max(a,b))
{
	return IAOverloaded::MaxMin<T,S>::_max(a,b);
}

template<typename T,typename S>
auto min(const T&a, const S&b) -> decltype( IAOverloaded::MaxMin<T,S>::_min(a,b))
{
	return IAOverloaded::MaxMin<T,S>::_min(a,b);
}


}//end ns ia


#endif // IA_IAOverloaded_HPP
