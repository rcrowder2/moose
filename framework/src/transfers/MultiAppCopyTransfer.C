//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

// MOOSE includes
#include "MultiAppCopyTransfer.h"
#include "FEProblemBase.h"
#include "MooseVariableFEBase.h"
#include "MooseMesh.h"
#include "SystemBase.h"

#include "libmesh/system.h"
#include "libmesh/id_types.h"
#include "libmesh/string_to_enum.h"

#include<algorithm>

registerMooseObject("MooseApp", MultiAppCopyTransfer);

defineLegacyParams(MultiAppCopyTransfer);

InputParameters
MultiAppCopyTransfer::validParams()
{
  MooseEnum reduction_types("COPY SUM AVG MIN MAX PROD", "COPY");

  InputParameters params = MultiAppFieldTransfer::validParams();
  params.addRequiredParam<std::vector<AuxVariableName>>(
      "variable", "The auxiliary variable to store the transferred values in.");
  params.addRequiredParam<std::vector<VariableName>>("source_variable",
                                                     "The variable to transfer from.");
  params.addParam<MooseEnum>("reduction",reduction_types, "The type of reduction to perform on the multiapps.");

  params.addClassDescription(
      "Copies variables (nonlinear and auxiliary) between multiapps that have identical meshes.");
  return params;
}

MultiAppCopyTransfer::MultiAppCopyTransfer(const InputParameters & parameters)
  : MultiAppFieldTransfer(parameters),
    _from_var_names(getParam<std::vector<VariableName>>("source_variable")),
    _to_var_names(getParam<std::vector<AuxVariableName>>("variable")),
    _reduction_type(getParam<MooseEnum>("reduction").getEnum<ReductionType>())
{
  /* Right now, most of transfers support one variable only */
  _to_var_name = _to_var_names[0];
  _from_var_name = _from_var_names[0];

  if (_reduction_type != ReductionType::COPY && _current_direction == TO_MULTIAPP)
    paramError("direction", "TO_MULTIAPP is only supported for COPY reductions");
}

void
MultiAppCopyTransfer::execute()
{
  _console << "Beginning MultiAppCopyTransfer " << name() << std::endl;

  if (_current_direction == TO_MULTIAPP)
  {
    FEProblemBase & from_problem = _multi_app->problemBase();
    for (unsigned int i = 0; i < _multi_app->numGlobalApps(); i++)
      if (_multi_app->hasLocalApp(i))
        transfer(_multi_app->appProblemBase(i), from_problem);
  }

  else if (_current_direction == FROM_MULTIAPP)
  {
    FEProblemBase & to_problem = _multi_app->problemBase();
    for (unsigned int i = 0; i < _multi_app->numGlobalApps(); i++)
      if (_multi_app->hasLocalApp(i))
        transfer(to_problem, _multi_app->appProblemBase(i));
  }

  _console << "Finished MultiAppCopyTransfer " << name() << std::endl;
}

void
MultiAppCopyTransfer::transferDofObject(libMesh::DofObject * to_object,
                                         libMesh::DofObject * from_object,
                                         MooseVariableFEBase & to_var,
                                         MooseVariableFEBase & from_var,
                                         NumericVector<Number> & to_solution,
                                         NumericVector<Number> & from_solution)
{
  for (unsigned int vc = 0; vc < to_var.count(); ++vc)
    if (to_object->n_dofs(to_var.sys().number(), to_var.number() + vc) >
        0) // If this variable has dofs at this node
      for (unsigned int comp = 0;
           comp < to_object->n_comp(to_var.sys().number(), to_var.number() + vc);
           ++comp)
      {
        dof_id_type dof = to_object->dof_number(to_var.sys().number(), to_var.number() + vc, comp);
        dof_id_type from_dof =
            from_object->dof_number(from_var.sys().number(), from_var.number() + vc, comp);
        Real from_value = from_solution(from_dof);
        Real to_value = to_solution(dof);

        switch(_reduction_type){
          case ReductionType::COPY:
            to_solution.set(dof, from_value);
            break;

          case ReductionType::SUM:
            to_solution.set(dof, (to_value + from_value));
            break;

          case ReductionType::AVG:
            to_solution.set(dof, to_value + from_value/_multi_app->numGlobalApps());
            break;

          case ReductionType::MIN:
            // will this work if to_value initialized as 0?
            to_solution.set(dof, std::min(to_value, from_value));

          case ReductionType::MAX:
            to_solution.set(dof, std::max(to_value, from_value));
            break;

          case ReductionType::PROD:
            to_solution.set(dof, to_value * from_value);
            break;
        }
      }
}

void
MultiAppCopyTransfer::transfer(FEProblemBase & to_problem, FEProblemBase & from_problem)
{
  // Perform error checking
  if (getToVarNames().size() != getFromVarNames().size())
    mooseError("Number of variables transfered must be same in both systems.");
  for (auto & to_var : getToVarNames())
    checkVariable(to_problem, to_var);
  for (auto & from_var : getFromVarNames())
    checkVariable(from_problem, from_var);

  for (unsigned int v = 0; v < getToVarNames().size(); ++v)
  {
    // Populate the to/from variables needed to perform the transfer
    MooseVariableFEBase & to_var = to_problem.getVariable(
        0, getToVarNames()[v], Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_ANY);
    MeshBase & to_mesh = to_problem.mesh().getMesh();

    MooseVariableFEBase & from_var = from_problem.getVariable(
        0, getFromVarNames()[v], Moose::VarKindType::VAR_ANY, Moose::VarFieldType::VAR_FIELD_ANY);
    MeshBase & from_mesh = from_problem.mesh().getMesh();

    auto & to_solution = isParamValid("to_solution_tag")
                             ? to_var.sys().getVector(
                                   to_problem.getVectorTagID(getParam<TagName>("to_solution_tag")))
                             : to_var.sys().solution();
    auto & from_solution = isParamValid("from_solution_tag")
                               ? from_var.sys().getVector(from_problem.getVectorTagID(
                                     getParam<TagName>("from_solution_tag")))
                               : from_var.sys().solution();

    // Check integrity
    if (to_var.feType() != from_var.feType())
      paramError("variable",
                 "Corresponding 'variable' and 'source_variable' inputs must be the same type "
                 "(order and family): ",
                 libMesh::Utility::enum_to_string<FEFamily>(to_var.feType().family),
                 moose::internal::incompatVarMsg(to_var, from_var));
    if (to_var.fieldType() != from_var.fieldType())
      mooseError(
          "Corresponding transfer variables must be same field type (STANDARD | VECTOR | ARRAY).");
    if (to_var.fieldType() == Moose::VarFieldType::VAR_FIELD_VECTOR)
      mooseError("Unable to transfer vector variables.");
    if (to_var.count() != from_var.count())
      mooseError("Corresponding transfer variables must have same number of components.");

    if ((to_mesh.n_nodes() != from_mesh.n_nodes()) || (to_mesh.n_elem() != from_mesh.n_elem()))
      mooseError("The meshes must be identical to utilize MultiAppCopyTransfer.");

    // Transfer node dofs
    for (const auto & node : as_range(to_mesh.local_nodes_begin(), to_mesh.local_nodes_end()))
      transferDofObject(
          node, from_mesh.node_ptr(node->id()), to_var, from_var, to_solution, from_solution);

    // Transfer elem dofs
    for (auto & to_elem : as_range(to_mesh.local_elements_begin(), to_mesh.local_elements_end()))
    {
      Elem * from_elem = from_mesh.elem_ptr(to_elem->id());
      mooseAssert(to_elem->type() == from_elem->type(), "The elements must be the same type.");
      transferDofObject(to_elem, from_elem, to_var, from_var, to_solution, from_solution);
    }

    to_solution.close();
    to_var.sys().update();
  }
}
