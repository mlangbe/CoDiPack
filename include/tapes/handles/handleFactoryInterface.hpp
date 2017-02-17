/*
 * CoDiPack, a Code Differentiation Package
 *
 * Copyright (C) 2015 Chair for Scientific Computing (SciComp), TU Kaiserslautern
 * Homepage: http://www.scicomp.uni-kl.de
 * Contact:  Prof. Nicolas R. Gauger (codi@scicomp.uni-kl.de)
 *
 * Lead developers: Max Sagebaum, Tim Albring (SciComp, TU Kaiserslautern)
 *
 * This file is part of CoDiPack (http://www.scicomp.uni-kl.de/software/codi).
 *
 * CoDiPack is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation, either version 2 of the
 * License, or (at your option) any later version.
 *
 * CoDiPack is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty
 * of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU General Public License for more details.
 * You should have received a copy of the GNU
 * General Public License along with CoDiPack.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 * Authors: Max Sagebaum, Tim Albring, (SciComp, TU Kaiserslautern)
 */

#pragma once

#include "../../configure.h"
#include "../../typeTraits.hpp"

/**
 * @brief Global namespace for CoDiPack - Code Differentiation Package
 */
namespace codi {

  /**
   * @brief The interface for the handle creation of a primal value tape.
   *
   * The interface has two function one for creating the handles and one for evaluating the handles.
   *
   * The handles that are create evaluate the reverse interpretation of an expression.
   *
   * @tparam    HandleType  The type of the handle, that this factory provides.
   * @tparam          Real  A calculation type that supports all mathematical operations.
   * @tparam     IndexType  An integer type that is used to identify the AD objects.
   * @tparam GradientValue  A value type that supports add and scaling operations.
   */
  template<typename HandleType, typename Real, typename IndexType, typename GradientValue=Real>
  struct HandleFactoryInterface {

    /** @brief The passive value of the Real type */
    typedef typename TypeTraits<Real>::PassiveReal PassiveReal;

    /** @brief The handle type, that is provided by this factory.*/
    typedef HandleType Handle;

    /**
     * @brief Create the handle for the given tape and the given expression.
     *
     * @tparam Expr  The expression that performs the evaluation of the reverse AD operations.
     * @tparam Tape  The tape that is performing the reverse AD evaluation.
     */
    template<typename Expr, typename Tape>
    static CODI_INLINE Handle createHandle();

    /**
     * @brief The evaluation of the handle, that was created by this factory.
     *
     * @param[in]           handle  The handle the was generated by this factory and is called with the arguments.
     * @param[in]              adj  The seeding for the expression
     * @param[in]   passiveActives  The number of passive values in the expression call.
     * @param[in,out]     indexPos  The position in the index buffer.
     * @param[in]          indices  The index buffer.
     * @param[in,out]  constantPos  The position in the constant value buffer.
     * @param[in]        constants  The constant value buffer.
     * @param[in,out] primalVector  The vector with the values of the primal variables.
     * @param[in,out]     adjoints  THe vector with the values of the adjoint variables.
     *
     * @tparam Tape  The tape that is performing the reverse AD evaluation.
     */
    template<typename Tape>
    static CODI_INLINE void callHandle(Handle handle, const GradientValue& adj, const StatementInt& passiveActives, size_t& indexPos, IndexType* &indices, size_t& constantPos, PassiveReal* &constants, Real* primalVector, GradientValue* adjoints);
  };
}