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

#include <vector>
#include <string>

#include <Rcpp.h>

#include "logics_impl.hpp"

using namespace Rcpp;

RCPP_MODULE(Logics){

//' @name Variable
//' @title Represent a simulation quantity
//' @description The objects of this class represent one among the following quantities:
//'
//'   -  the cardinality of a species;
//'   -  the number fired event among deaths, duplications and switches in a species;
//'   -  the elapse simulation time.
//'
//' @seealso `Simulation$var()` to build a variable
  class_<Logics::Variable>("Variable")
    .method("show", &Logics::Variable::show,
            "Show a variable");

//' @name Expression
//' @title Represent an expression of value and variable
//' @description The objects of this class are polynomial expression of variables
//'   and values built by using the binary operators `+`, `-`, and `*` with the 
//'   usual semantics.
//'
//'   An expression is one of the following object:
//'
//'   -  a variable;
//'   -  a numeric value, e.g., `3.4`;
//'   -  the sum of two expressions;
//'   -  the subtraction of two expressions;
//'   -  the multiplication of two expressions.
//'
//' @seealso `Variable`
  class_<Logics::Expression>("Expression")
    .method("show", &Logics::Expression::show,
            "Show an expression");

  function("logics_sum", (SEXP (*)(const SEXP&, const SEXP&))&Logics::sum,
           "Sums two expressions");
  function("logics_subtract", (SEXP (*)(const SEXP&, const SEXP&))&Logics::subtract,
           "Sums two expressions");
  function("logics_multiply", (SEXP (*)(const SEXP&, const SEXP&))&Logics::multiply,
           "Multiplies two expressions");

/*
  class_<Logics::Relation>("Relation")
    .method("show", &Logics::Relation::show,
            "Show a relation");
*/

  function("logics_lt", (SEXP (*)(const SEXP&, const SEXP&))&Logics::lt,
           "Builds a \"less than\" relation");
  function("logics_le", (SEXP (*)(const SEXP&, const SEXP&))&Logics::le,
           "Builds a \"less than or equal to\" relation");
  function("logics_eq", (SEXP (*)(const SEXP&, const SEXP&))&Logics::eq,
           "Builds an \"equal to\" relation");
  function("logics_ne", (SEXP (*)(const SEXP&, const SEXP&))&Logics::ne,
           "Builds an \"not equal to\" relation");
  function("logics_ge", (SEXP (*)(const SEXP&, const SEXP&))&Logics::ge,
           "Builds an \"greater than or equal to\" relation");
  function("logics_gt", (SEXP (*)(const SEXP&, const SEXP&))&Logics::gt,
           "Builds an \"greater than\" relation");

//' @name Formula
//' @title First order formulas about simulations
//' @description The objects of this class describes properties of the 
//'   simulation status.
//'
//'   A formula is:
//'
//'   -  a relation among two expressions (operators `<`, `<=`, `==`,
//'      `!=`, `>=`, `>`);
//'   -  the conjunction of two formulas (operator `&`);
//'   -  the disjunction of two formulas (operator `|`).
//'
//' @seealso `Variable`, `Simulation$var()`, `vignette("tissue_simulation")`
  class_<Logics::Formula>("Formula")
    .method("show", &Logics::Formula::show,
            "Show a formula");

  function("not", (SEXP (*)(const SEXP&))&Logics::logics_not,
           "Builds the negation of a formula");
  function("logics_and", (SEXP (*)(const SEXP& lhs, const SEXP& rhs))&Logics::logics_and,
           "Builds the logical conjunction of two formulas");
  function("logics_or", (SEXP (*)(const SEXP&, const SEXP&))&Logics::logics_or,
           "Builds the logical disjunction of two formulas");
}
