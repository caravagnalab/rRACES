/*
 * This file is part of the rRACES (https://github.com/caravagnalab/rRACES/).
 * Copyright (c) 2023-2024 Alberto Casagrande <alberto.casagrande@uniud.it>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __RRACES_LOGICS_IMPL__
#define __RRACES_LOGICS_IMPL__

#include <Rcpp.h>

#include <logics.hpp>

class SpatialSimulation;

namespace Logics
{

class Variable: public RACES::Mutants::Logics::Variable
{
public: 
  Variable(RACES::Mutants::Logics::Variable&& variable);

  void show() const;

  friend class SpatialSimulation;
};

class Expression : public RACES::Mutants::Logics::Expression
{
public:
  Expression(RACES::Mutants::Logics::Expression&& expression);

  void show() const;
};

SEXP sum(const SEXP& lhs, const SEXP& rhs);
SEXP subtract(const SEXP& lhs, const SEXP& rhs);
SEXP multiply(const SEXP& lhs, const SEXP& rhs);

class Relation : public RACES::Mutants::Logics::Relation
{
public:
  Relation(RACES::Mutants::Logics::Relation&& relation);

  void show() const;
};

SEXP gt(const SEXP& lhs, const SEXP& rhs);
SEXP ge(const SEXP& lhs, const SEXP& rhs);

SEXP eq(const SEXP& lhs, const SEXP& rhs);
SEXP ne(const SEXP& lhs, const SEXP& rhs);

SEXP lt(const SEXP& lhs, const SEXP& rhs);
SEXP le(const SEXP& lhs, const SEXP& rhs);

class Formula: public RACES::Mutants::Logics::Formula
{
public:
  Formula(RACES::Mutants::Logics::Formula&& Formula);

  void show() const;
};

SEXP logics_not(const SEXP& subformula);
SEXP logics_and(const SEXP& lhs, const SEXP& rhs);
SEXP logics_or(const SEXP& lhs, const SEXP& rhs);

}

RCPP_EXPOSED_WRAP(Logics::Variable);
RCPP_EXPOSED_WRAP(Logics::Expression);
RCPP_EXPOSED_WRAP(Logics::Relation);
RCPP_EXPOSED_WRAP(Logics::Formula);
RCPP_EXPOSED_AS(Logics::Formula);

#endif // __RRACES_LOGICS_IMPL__
