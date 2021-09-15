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
#include "MultiApp.h"

registerMooseObject("MooseApp", MultiAppCopyTransfer);

defineLegacyParams(MultiAppCopyTransfer);

InputParameters
MultiAppCopyTransfer::validParams()
{
  MooseEnum reduction_types("COPY SUM AVG", "COPY");

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

  if (_reduction_type == ReductionType::AVG)
    paramError("reduction", "AVG reduction type is not currently supported");
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
        if (_reduction_type == ReductionType::COPY)
            to_solution.set(dof, from_value);
        else if (_reduction_type == ReductionType::SUM)
            to_solution.set(dof, (to_solution(dof) + from_value)/_multi_app->numGlobalApps());
      }
}
