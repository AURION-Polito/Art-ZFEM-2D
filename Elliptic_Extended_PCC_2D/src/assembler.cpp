// _LICENSE_HEADER_
//
// Copyright (C) 2019 - 2025.
// Terms register on the GPL-3.0 license.
//
// This file can be redistributed and/or modified under the license terms.
//
// See top level LICENSE file for more details.
//
// This file can be used citing references in CITATION.cff file.

#include "assembler.hpp"

#include "Assembler_Utilities.hpp"
#include "EllipticEquation.hpp"
#include "Quadrature_Gauss1D.hpp"
#include "VEM_PCC_Utilities.hpp"

namespace Polydim
{
namespace examples
{
namespace Elliptic_Extended_PCC_2D
{
//***************************************************************************
void Assembler::ComputeStrongTerm(const Gedim::GeometryUtilities& geometry_utilities,
                                  const unsigned int cell2D_index,
                                  const unsigned int cell2D_domain_position,
                                  const Eigen::Matrix3d& domain_rotation,
                                  const Eigen::Vector3d& domain_translation,
                                  const Gedim::MeshMatricesDAO &mesh,
                                  const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                  const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                  const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                  const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data,
                                  const test::I_Test &test,
                                  Elliptic_Extended_PCC_2D_Problem_Data &assembler_data) const
{
    // Assemble strong boundary condition on Cell0Ds
    for (unsigned int v = 0; v < mesh.Cell2DNumberVertices(cell2D_index); ++v)
    {
        const unsigned int cell0D_index = mesh.Cell2DVertex(cell2D_index, v);
        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(0).at(cell0D_index);

        if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong)
            continue;

        const auto coordinates_3D = mesh.Cell0DCoordinates(cell0D_index);
        const auto strong_boundary_values = test.strong_boundary_condition(cell2D_domain_position,
                                                                           boundary_info.Marker,
                                                                           coordinates_3D);

        const auto local_dofs = dofs_data.CellsDOFs.at(0).at(cell0D_index);

        assert(local_dofs.size() == strong_boundary_values.size());

        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
        {
            const auto &local_dof_i = local_dofs.at(loc_i);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong: {
                assembler_data.solutionDirichlet.SetValue(local_dof_i.Global_Index, strong_boundary_values[loc_i]);
            }
            break;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                continue;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }

    // Assemble strong boundary condition on Cell1Ds
    for (unsigned int ed = 0; ed < mesh.Cell2DNumberEdges(cell2D_index); ++ed)
    {
        const unsigned int cell1D_index = mesh.Cell2DEdge(cell2D_index, ed);

        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(1).at(cell1D_index);
        const auto local_dofs = dofs_data.CellsDOFs.at(1).at(cell1D_index);

        if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong ||
            local_dofs.size() == 0)
            continue;

        const auto edge_dofs_coordinates_2D =
            Polydim::PDETools::LocalSpace_PCC_2D::EdgeDofsCoordinates(reference_element_data, local_space_data, ed);
        const auto coordinates_3D = geometry_utilities.RotatePointsFrom2DTo3D(edge_dofs_coordinates_2D,
                                                                              domain_rotation,
                                                                              domain_translation);

        const auto strong_boundary_values = test.strong_boundary_condition(cell2D_domain_position,
                                                                           boundary_info.Marker,
                                                                           coordinates_3D);

        assert(local_dofs.size() == strong_boundary_values.size());

        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
        {
            const auto &local_dof_i = local_dofs.at(loc_i);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong: {
                assembler_data.solutionDirichlet.SetValue(local_dof_i.Global_Index, strong_boundary_values[loc_i]);
            }
            break;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                continue;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }
}
// ***************************************************************************
void Assembler::ComputeWeakTerm(const Gedim::GeometryUtilities& geometry_utilities,
                                const unsigned int cell2DIndex,
                                const unsigned int cell2D_domain_position,
                                const Eigen::Matrix3d& domain_rotation,
                                const Eigen::Vector3d& domain_translation,
                                const Gedim::MeshMatricesDAO &mesh,
                                const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data,
                                const Polydim::examples::Elliptic_Extended_PCC_2D::test::I_Test &test,
                                Elliptic_Extended_PCC_2D_Problem_Data &assembler_data) const
{
    const unsigned numVertices = mesh_geometric_data.Cell2DsVertices.at(cell2DIndex).cols();

    for (unsigned int ed = 0; ed < numVertices; ed++)
    {
        const unsigned int cell1D_index = mesh.Cell2DEdge(cell2DIndex, ed);

        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(1).at(cell1D_index);

        if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Weak)
            continue;

        // compute vem values
        const auto weakReferenceSegment =
            Gedim::Quadrature::Quadrature_Gauss1D::FillPointsAndWeights(2 * reference_element_data.Order);

        const Eigen::VectorXd pointsCurvilinearCoordinates = weakReferenceSegment.Points.row(0);

        // map edge internal quadrature points
        const Eigen::Vector3d &edgeStart = mesh_geometric_data.Cell2DsEdgeDirections.at(cell2DIndex)[ed]
                                               ? mesh_geometric_data.Cell2DsVertices.at(cell2DIndex).col(ed)
                                               : mesh_geometric_data.Cell2DsVertices.at(cell2DIndex).col((ed + 1) % numVertices);

        const Eigen::Vector3d &edgeTangent = mesh_geometric_data.Cell2DsEdgeTangents.at(cell2DIndex).col(ed);
        const double direction = mesh_geometric_data.Cell2DsEdgeDirections.at(cell2DIndex)[ed] ? 1.0 : -1.0;

        const unsigned int numEdgeWeakQuadraturePoints = weakReferenceSegment.Points.cols();
        Eigen::MatrixXd weakQuadraturePoints_2D(3, numEdgeWeakQuadraturePoints);
        for (unsigned int q = 0; q < numEdgeWeakQuadraturePoints; q++)
            weakQuadraturePoints_2D.col(q) = edgeStart + direction * weakReferenceSegment.Points(0, q) * edgeTangent;

        const auto coordinates_3D = geometry_utilities.RotatePointsFrom2DTo3D(weakQuadraturePoints_2D,
                                                                              domain_rotation,
                                                                              domain_translation);

        const double absMapDeterminant = std::abs(mesh_geometric_data.Cell2DsEdgeLengths.at(cell2DIndex)[ed]);
        const Eigen::MatrixXd weakQuadratureWeights = weakReferenceSegment.Weights * absMapDeterminant;

        const Eigen::VectorXd neumannValues = test.weak_boundary_condition(cell2D_domain_position,
                                                                           boundary_info.Marker,
                                                                           coordinates_3D);
        const auto weak_basis_function_values =
            Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsValuesOnEdge(ed, reference_element_data, local_space_data, pointsCurvilinearCoordinates);

        // compute values of Neumann condition
        const Eigen::VectorXd neumannContributions =
            weak_basis_function_values.transpose() * weakQuadratureWeights.asDiagonal() * neumannValues;

        for (unsigned int p = 0; p < 2; ++p)
        {
            const unsigned int cell0D_index = mesh.Cell1DVertex(cell1D_index, p);

            const auto local_dofs = dofs_data.CellsDOFs.at(0).at(cell0D_index);

            for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
            {
                const auto &local_dof_i = local_dofs.at(loc_i);

                switch (local_dof_i.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                    continue;
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF: {
                    assembler_data.rightHandSide.AddValue(local_dof_i.Global_Index, neumannContributions[p]);
                }
                break;
                default:
                    throw std::runtime_error("Unknown DOF Type");
                }
            }
        }

        const auto local_dofs = dofs_data.CellsDOFs.at(1).at(cell1D_index);
        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
        {
            const auto &local_dof_i = local_dofs.at(loc_i);

            const unsigned int localIndex = loc_i;

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                continue;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF: {
                assembler_data.rightHandSide.AddValue(local_dof_i.Global_Index, neumannContributions[localIndex + 2]);
            }
            break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }
}
// ***************************************************************************
Assembler::Elliptic_Extended_PCC_2D_Problem_Data Assembler::Assemble(
    const Polydim::examples::Elliptic_Extended_PCC_2D::Program_configuration &config,
    const PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D_Collection &domains,
    const Gedim::MeshMatricesDAO &mesh,
    const PDETools::Mesh::PDE_Mesh_Utilities::Extended_MeshGeometricData2D &mesh_geometric_data,
    const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
    const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
    const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
    const Polydim::examples::Elliptic_Extended_PCC_2D::test::I_Test &test) const
{
    Gedim::GeometryUtilitiesConfig geometry_utilities_config;
    geometry_utilities_config.Tolerance1D = config.GeometricTolerance1D();
    geometry_utilities_config.Tolerance2D = config.GeometricTolerance2D();
    Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    Elliptic_Extended_PCC_2D_Problem_Data result;

    result.globalMatrixA.SetSize(dofs_data.NumberDOFs, dofs_data.NumberDOFs, Gedim::ISparseArray::SparseArrayTypes::None);
    result.dirichletMatrixA.SetSize(dofs_data.NumberDOFs, dofs_data.NumberStrongs);
    result.rightHandSide.SetSize(dofs_data.NumberDOFs);
    result.solution.SetSize(dofs_data.NumberDOFs);
    result.solutionDirichlet.SetSize(dofs_data.NumberStrongs);

    Polydim::PDETools::Equations::EllipticEquation equation;

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); ++c)
    {
      if (!mesh.Cell2DIsActive(c))
        continue;

      const auto cell2D_domain_position = mesh_geometric_data.cell2Ds_domain_position.at(c);
      const auto& domain_rotation = domains.domains_rotation.at(cell2D_domain_position);
      const auto& domain_translation = domains.domains_translation.at(cell2D_domain_position);

        const auto local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(config.GeometricTolerance1D(),
                                                                                             config.GeometricTolerance2D(),
                                                                                             mesh_geometric_data.mesh_geometric_data,
                                                                                             c,
                                                                                             reference_element_data);

        const auto basis_functions_values =
            Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsValues(reference_element_data, local_space_data, Polydim::VEM::PCC::ProjectionTypes::Pi0k);

        const auto basis_functions_derivative_values =
            Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsDerivativeValues(reference_element_data, local_space_data);

        const auto cell2D_internal_quadrature =
            Polydim::PDETools::LocalSpace_PCC_2D::InternalQuadrature(reference_element_data, local_space_data);

        const auto cell2D_internal_quadrature_3D = geometry_utilities.RotatePointsFrom2DTo3D(cell2D_internal_quadrature.Points,
                                                                                             domain_rotation,
                                                                                             domain_translation);

        const auto diffusion_term_values = test.diffusion_term(cell2D_domain_position,
                                                               cell2D_internal_quadrature_3D);
        const auto source_term_values = test.source_term(cell2D_domain_position,
                                                         cell2D_internal_quadrature_3D);

        const Eigen::MatrixXd local_A = equation.ComputeCellDiffusionMatrix(diffusion_term_values,
                                                                            basis_functions_derivative_values,
                                                                            cell2D_internal_quadrature.Weights);

        Eigen::VectorXd local_rhs =
            equation.ComputeCellForcingTerm(source_term_values, basis_functions_values, cell2D_internal_quadrature.Weights);

        const double k_max = diffusion_term_values.cwiseAbs().maxCoeff();
        const double &diameter = mesh_geometric_data.mesh_geometric_data.Cell2DsDiameters.at(c);

        const Eigen::MatrixXd local_A_stab =
            k_max *
            Polydim::PDETools::LocalSpace_PCC_2D::StabilizationMatrix(reference_element_data, local_space_data);

        const auto &global_dofs = dofs_data.CellsGlobalDOFs[2].at(c);

        assert(Polydim::PDETools::LocalSpace_PCC_2D::Size(reference_element_data, local_space_data) == global_dofs.size());

        Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data =
            {{std::cref(dofs_data)}, {0}, {0}, {0}};

        Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<2>(c,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_A + local_A_stab,
                                                                                          local_rhs,
                                                                                          result.globalMatrixA,
                                                                                          result.dirichletMatrixA,
                                                                                          result.rightHandSide);

        ComputeStrongTerm(geometry_utilities,
                          c,
                          cell2D_domain_position,
                          domain_rotation,
                          domain_translation,
                          mesh, mesh_dofs_info, dofs_data, reference_element_data, local_space_data, test, result);

        ComputeWeakTerm(geometry_utilities,
                        c,
                        cell2D_domain_position,
                        domain_rotation,
                        domain_translation,
                        mesh, mesh_geometric_data.mesh_geometric_data, mesh_dofs_info, dofs_data, reference_element_data, local_space_data, test, result);
    }

    result.rightHandSide.Create();
    result.solutionDirichlet.Create();
    result.globalMatrixA.Create();
    result.dirichletMatrixA.Create();

    if (dofs_data.NumberStrongs > 0)
        result.rightHandSide.SubtractionMultiplication(result.dirichletMatrixA, result.solutionDirichlet);

    return result;
}
// ***************************************************************************
Assembler::Performance_Data Assembler::ComputePerformance(const Polydim::examples::Elliptic_Extended_PCC_2D::Program_configuration &config,
                                                          const Gedim::MeshMatricesDAO &mesh,
                                                          const PDETools::Mesh::PDE_Mesh_Utilities::Extended_MeshGeometricData2D &mesh_geometric_data,
                                                          const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data) const
{
    Assembler::Performance_Data result;
    result.Cell2DsPerformance.resize(mesh.Cell2DTotalNumber());

    // Assemble equation elements
    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
      if (!mesh.Cell2DIsActive(c))
        continue;

        const auto local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(config.GeometricTolerance1D(),
                                                                                             config.GeometricTolerance2D(),
                                                                                             mesh_geometric_data.mesh_geometric_data,
                                                                                             c,
                                                                                             reference_element_data);

        result.Cell2DsPerformance[c] =
            Polydim::PDETools::LocalSpace_PCC_2D::ComputePerformance(reference_element_data, local_space_data);
    }

    return result;
}
// ***************************************************************************
Assembler::PostProcess_Data Assembler::PostProcessSolution(const Polydim::examples::Elliptic_Extended_PCC_2D::Program_configuration &config,
                                                           const PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D_Collection &domains,
                                                           const Gedim::MeshMatricesDAO &mesh,
                                                           const PDETools::Mesh::PDE_Mesh_Utilities::Extended_MeshGeometricData2D &mesh_geometric_data,
                                                           const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                                           const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                                           const Elliptic_Extended_PCC_2D_Problem_Data &assembler_data,
                                                           const Polydim::examples::Elliptic_Extended_PCC_2D::test::I_Test &test) const
{
  Gedim::GeometryUtilitiesConfig geometry_utilities_config;
  geometry_utilities_config.Tolerance1D = config.GeometricTolerance1D();
  geometry_utilities_config.Tolerance2D = config.GeometricTolerance2D();
  Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

    PostProcess_Data result;

    result.residual_norm = 0.0;
    if (dofs_data.NumberDOFs > 0)
    {
        Gedim::Eigen_Array<> residual;
        residual.SetSize(dofs_data.NumberDOFs);
        residual.SumMultiplication(assembler_data.globalMatrixA, assembler_data.solution);
        residual -= assembler_data.rightHandSide;

        result.residual_norm = residual.Norm();
    }

    result.cell0Ds_numeric.setZero(mesh.Cell0DTotalNumber());
    for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); p++)
    {
        const auto local_dofs = dofs_data.CellsDOFs.at(0).at(p);

        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
        {
            const auto &local_dof_i = local_dofs.at(loc_i);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                result.cell0Ds_numeric[p] = assembler_data.solutionDirichlet.GetValue(local_dof_i.Global_Index);
                break;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                result.cell0Ds_numeric[p] = assembler_data.solution.GetValue(local_dof_i.Global_Index);
                break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }

    result.cell0Ds_exact.setZero(mesh.Cell0DTotalNumber());
    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); ++c)
    {
      if (!mesh.Cell2DIsActive(c))
        continue;

      const auto cell2D_domain_position = mesh_geometric_data.cell2Ds_domain_position.at(c);
      const auto& domain_rotation = domains.domains_rotation.at(cell2D_domain_position);
      const auto& domain_translation = domains.domains_translation.at(cell2D_domain_position);

      for (unsigned int v = 0; v < mesh.Cell2DNumberVertices(c); ++v)
      {
        result.cell0Ds_exact[mesh.Cell2DVertex(c, v)] =
            test.exact_solution(cell2D_domain_position,
                                mesh.Cell2DVertexCoordinates(c, v))[0];

      }
    }

    result.cell2Ds_error_L2.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_numeric_norm_L2.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_exact_norm_L2.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_error_H1.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_numeric_norm_H1.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_exact_norm_H1.setZero(mesh.Cell2DTotalNumber());
    result.error_L2 = 0.0;
    result.numeric_norm_L2 = 0.0;
    result.error_H1 = 0.0;
    result.numeric_norm_H1 = 0.0;
    result.mesh_size = 0.0;

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
      if (!mesh.Cell2DIsActive(c))
        continue;

      const auto cell2D_domain_position = mesh_geometric_data.cell2Ds_domain_position.at(c);
      const auto& domain_rotation = domains.domains_rotation.at(cell2D_domain_position);
      const auto& domain_translation = domains.domains_translation.at(cell2D_domain_position);


        const auto local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(config.GeometricTolerance1D(),
                                                                                             config.GeometricTolerance2D(),
                                                                                             mesh_geometric_data.mesh_geometric_data,
                                                                                             c,
                                                                                             reference_element_data);

        const auto basis_functions_values =
            Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsValues(reference_element_data,
                                                                       local_space_data,
                                                                       Polydim::VEM::PCC::ProjectionTypes::Pi0k);

        const auto basis_functions_derivative_values =
            Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsDerivativeValues(reference_element_data, local_space_data);

        const auto cell2D_internal_quadrature =
            Polydim::PDETools::LocalSpace_PCC_2D::InternalQuadrature(reference_element_data, local_space_data);

        const auto cell2D_internal_quadrature_3D = geometry_utilities.RotatePointsFrom2DTo3D(cell2D_internal_quadrature.Points,
                                                                                             domain_rotation,
                                                                                             domain_translation);


        const auto exact_solution_values = test.exact_solution(cell2D_domain_position,
                                                               cell2D_internal_quadrature_3D);
        const auto exact_derivative_solution_values = test.exact_derivative_solution(cell2D_domain_position,
                                                                                     cell2D_internal_quadrature_3D);

        const auto local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(c, dofs_data);
        const Eigen::VectorXd dofs_values =
            PDETools::Assembler_Utilities::global_solution_to_local_solution<2>(c,
                                                                                dofs_data,
                                                                                local_count_dofs.num_total_dofs,
                                                                                local_count_dofs.offsets_DOFs,
                                                                                {0},
                                                                                {0},
                                                                                assembler_data.solution,
                                                                                assembler_data.solutionDirichlet);

        const Eigen::VectorXd local_error_L2 = (basis_functions_values * dofs_values - exact_solution_values).array().square();
        const Eigen::VectorXd local_norm_L2 = (basis_functions_values * dofs_values).array().square();
        const Eigen::VectorXd local_exact_norm_L2 = (exact_solution_values).array().square();

        result.cell2Ds_error_L2[c] = cell2D_internal_quadrature.Weights.transpose() * local_error_L2;
        result.cell2Ds_numeric_norm_L2[c] = cell2D_internal_quadrature.Weights.transpose() * local_norm_L2;
        result.cell2Ds_exact_norm_L2[c] = cell2D_internal_quadrature.Weights.transpose() * local_exact_norm_L2;

        std::vector<Eigen::VectorXd> numerical_derivatives_solution_2D(2);
        numerical_derivatives_solution_2D[0] = basis_functions_derivative_values.at(0) *
                                            dofs_values;
        numerical_derivatives_solution_2D[1] = basis_functions_derivative_values.at(1) *
                                            dofs_values;

        std::vector<Eigen::VectorXd> numerical_derivatives_solution_3D(3);
        numerical_derivatives_solution_3D[0] =
            domain_rotation(0, 0) * numerical_derivatives_solution_2D.at(0) +
            domain_rotation(0, 1) * numerical_derivatives_solution_2D.at(1);
        numerical_derivatives_solution_3D[1] =
            domain_rotation(1, 0) * numerical_derivatives_solution_2D.at(0) +
            domain_rotation(1, 1) * numerical_derivatives_solution_2D.at(1);
        numerical_derivatives_solution_3D[2] =
            domain_rotation(2, 0) * numerical_derivatives_solution_2D.at(0) +
            domain_rotation(2, 1) * numerical_derivatives_solution_2D.at(1);

        const Eigen::VectorXd local_error_H1 =
            (numerical_derivatives_solution_3D.at(0) - exact_derivative_solution_values[0]).array().square() +
            (numerical_derivatives_solution_3D.at(1) - exact_derivative_solution_values[1]).array().square() +
            (numerical_derivatives_solution_3D.at(2) - exact_derivative_solution_values[2]).array().square();

        const Eigen::VectorXd local_norm_H1 = numerical_derivatives_solution_3D.at(0).array().square() +
                                              numerical_derivatives_solution_3D.at(1).array().square() +
                                              numerical_derivatives_solution_3D.at(2).array().square();

        const Eigen::VectorXd local_exact_norm_H1 = exact_derivative_solution_values.at(0).array().square() +
                                                    exact_derivative_solution_values.at(1).array().square() +
                                                    exact_derivative_solution_values.at(2).array().square();

        result.cell2Ds_error_H1[c] = cell2D_internal_quadrature.Weights.transpose() * local_error_H1;
        result.cell2Ds_numeric_norm_H1[c] = cell2D_internal_quadrature.Weights.transpose() * local_norm_H1;
        result.cell2Ds_exact_norm_H1[c] = cell2D_internal_quadrature.Weights.transpose() * local_exact_norm_H1;

        if (mesh_geometric_data.mesh_geometric_data.Cell2DsDiameters.at(c) > result.mesh_size)
            result.mesh_size = mesh_geometric_data.mesh_geometric_data.Cell2DsDiameters.at(c);
    }

    result.error_L2 = std::sqrt(result.cell2Ds_error_L2.sum());
    result.numeric_norm_L2 = std::sqrt(result.cell2Ds_numeric_norm_L2.sum());
    result.exact_norm_L2 = std::sqrt(result.cell2Ds_exact_norm_L2.sum());
    result.error_H1 = std::sqrt(result.cell2Ds_error_H1.sum());
    result.numeric_norm_H1 = std::sqrt(result.cell2Ds_numeric_norm_H1.sum());
    result.exact_norm_H1 = std::sqrt(result.cell2Ds_exact_norm_H1.sum());

    return result;
}
} // namespace Elliptic_Extended_PCC_2D
} // namespace examples
} // namespace Polydim
