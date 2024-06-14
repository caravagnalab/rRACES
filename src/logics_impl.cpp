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

#include "logics_impl.hpp"

namespace Logics
{

Variable::Variable(RACES::Mutants::Logics::Variable&& variable):
    RACES::Mutants::Logics::Variable(std::move(variable))
{}

void Variable::show() const
{
  Rcpp::Rcout << *this << std::endl;
}

Expression to_expression(const SEXP& exp)
{
  switch (TYPEOF(exp)) {
    case REALSXP:
    {
      return Expression(Rcpp::as<double>(exp));
    }
    case INTSXP:
    {
      double e_value = Rcpp::as<int>(exp);

      return Expression(e_value);
    }
    case S4SXP:
    {
      Rcpp::S4 s4obj( exp );
      if ( s4obj.is("Rcpp_Variable" )) {
        Rcpp::Environment env( s4obj );
        Rcpp::XPtr<Variable> var_ptr( env.get(".pointer") );

        return Expression(*var_ptr);
      }
      if ( s4obj.is("Rcpp_Expression" )) {
        Rcpp::Environment env( s4obj );
        Rcpp::XPtr<Expression> expr_ptr( env.get(".pointer") );

        return Expression(*expr_ptr);
      }
    }
    default:
      ::Rf_error("Unsupported expression type");
  }
}

Expression::Expression(RACES::Mutants::Logics::Expression&& expression):
    RACES::Mutants::Logics::Expression(std::move(expression))
{}

void Expression::show() const
{
  Rcpp::Rcout << *this << std::endl;
}

SEXP sum(const SEXP& lhs, const SEXP& rhs)
{
  return Rcpp::wrap(Expression(to_expression(lhs) + to_expression(rhs)));
}

SEXP subtract(const SEXP& lhs, const SEXP& rhs)
{
  return Rcpp::wrap(Expression(to_expression(lhs) - to_expression(rhs)));
}

SEXP multiply(const SEXP& lhs, const SEXP& rhs)
{
  return Rcpp::wrap(Expression(to_expression(lhs) * to_expression(rhs)));
}

Relation::Relation(RACES::Mutants::Logics::Relation&& relation):
    RACES::Mutants::Logics::Relation(std::move(relation))
{}

void Relation::show() const
{
  Rcpp::Rcout << static_cast<const RACES::Mutants::Logics::Relation&>(*this)
              << std::endl;
}

SEXP gt(const SEXP& lhs, const SEXP& rhs)
{
  return Rcpp::wrap(Formula(Relation(to_expression(lhs) > to_expression(rhs))));
}

SEXP ge(const SEXP& lhs, const SEXP& rhs)
{
  return Rcpp::wrap(Formula(Relation(to_expression(lhs) >= to_expression(rhs))));
}

SEXP eq(const SEXP& lhs, const SEXP& rhs)
{
  return Rcpp::wrap(Formula(Relation(to_expression(lhs) == to_expression(rhs))));
}

SEXP ne(const SEXP& lhs, const SEXP& rhs)
{
  return Rcpp::wrap(Formula(Relation(to_expression(lhs) != to_expression(rhs))));
}

SEXP le(const SEXP& lhs, const SEXP& rhs)
{
  return Rcpp::wrap(Formula(Relation(to_expression(lhs) <= to_expression(rhs))));
}

SEXP lt(const SEXP& lhs, const SEXP& rhs)
{
  return Rcpp::wrap(Formula(Relation(to_expression(lhs) < to_expression(rhs))));
}

Formula::Formula(RACES::Mutants::Logics::Formula&& formula):
    RACES::Mutants::Logics::Formula(std::move(formula))
{}

void Formula::show() const
{
  Rcpp::Rcout << static_cast<const RACES::Mutants::Logics::Formula&>(*this)
              << std::endl;
}

Formula to_formula(const SEXP& formula)
{
  switch (TYPEOF(formula)) {
    case S4SXP:
    {
      Rcpp::S4 s4obj( formula );
      if ( s4obj.is("Rcpp_Relation" )) {
        Rcpp::Environment env( s4obj );
        Rcpp::XPtr<Relation> rel_ptr( env.get(".pointer") );

        return Formula(*rel_ptr);
      }
      if ( s4obj.is("Rcpp_Formula" )) {
        Rcpp::Environment env( s4obj );
        Rcpp::XPtr<Formula> formula_ptr( env.get(".pointer") );

        return Formula(*formula_ptr);
      }
    }
    default:
      ::Rf_error("Unsupported type");
  }
}

SEXP logics_not(const SEXP& subformula)
{
  return Rcpp::wrap(Formula(!to_formula(subformula)));
}

SEXP logics_and(const SEXP& lhs, const SEXP& rhs)
{ 
  return Rcpp::wrap(Formula(to_formula(lhs) && to_formula(rhs)));
}

SEXP logics_or(const SEXP& lhs, const SEXP& rhs)
{
  return Rcpp::wrap(Formula(to_formula(lhs) || to_formula(rhs)));
}

}