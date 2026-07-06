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

namespace Polydim
{
namespace examples
{
namespace Parabolic_PCC_2D
{
//***************************************************************************
void Assembler::ComputeStrongTerm(const double &time_value,
                                  const unsigned int cell2DIndex,
                                  const unsigned int region_id,
                                  const Gedim::MeshMatricesDAO &mesh,
                                  const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                  const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                  const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                  const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data,
                                  const test::I_Test &test,
                                  Gedim::Eigen_Array<> &solution_dirichlet) const
{
    // Assemble strong boundary condition on Cell0Ds
    for (unsigned int v = 0; v < mesh.Cell2DNumberVertices(cell2DIndex); ++v)
    {
        const unsigned int cell0D_index = mesh.Cell2DVertex(cell2DIndex, v);
        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(0).at(cell0D_index);

        if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong)
            continue;

        const auto coordinates = mesh.Cell0DCoordinates(cell0D_index);

        const auto strong_boundary_values = test.strong_boundary_condition(boundary_info.Marker, coordinates, time_value);

        const auto local_dofs = dofs_data.CellsDOFs.at(0).at(cell0D_index);

        assert(local_dofs.size() == strong_boundary_values.size());

        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
        {
            const auto &local_dof_i = local_dofs.at(loc_i);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong: {
                solution_dirichlet.SetValue(local_dof_i.Global_Index, strong_boundary_values[loc_i]);
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
    for (unsigned int ed = 0; ed < mesh.Cell2DNumberEdges(cell2DIndex); ++ed)
    {
        const unsigned int cell1D_index = mesh.Cell2DEdge(cell2DIndex, ed);

        const auto &boundary_info = mesh_dofs_info.CellsBoundaryInfo.at(1).at(cell1D_index);
        const auto local_dofs = dofs_data.CellsDOFs.at(1).at(cell1D_index);

        if (boundary_info.Type != Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo::BoundaryTypes::Strong ||
            local_dofs.size() == 0)
            continue;

        const auto edge_dofs_coordinates =
            Polydim::PDETools::LocalSpace_PCC_2D::EdgeDofsCoordinates(reference_element_data, local_space_data, ed);

        const auto strong_boundary_values = test.strong_boundary_condition(boundary_info.Marker, edge_dofs_coordinates, time_value);

        assert(local_dofs.size() == strong_boundary_values.size());

        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
        {
            const auto &local_dof_i = local_dofs.at(loc_i);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong: {
                solution_dirichlet.SetValue(local_dof_i.Global_Index, strong_boundary_values[loc_i]);
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
void Assembler::ComputeWeakTerm(const double &time_value,
                                const unsigned int cell2DIndex,
                                const unsigned int region_id,
                                const Gedim::MeshMatricesDAO &mesh,
                                const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data,
                                const Polydim::examples::Parabolic_PCC_2D::test::I_Test &test,
                                Gedim::Eigen_Array<> &rightHandSide) const
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
        Eigen::MatrixXd weakQuadraturePoints(3, numEdgeWeakQuadraturePoints);
        for (unsigned int q = 0; q < numEdgeWeakQuadraturePoints; q++)
            weakQuadraturePoints.col(q) = edgeStart + direction * weakReferenceSegment.Points(0, q) * edgeTangent;

        const double absMapDeterminant = std::abs(mesh_geometric_data.Cell2DsEdgeLengths.at(cell2DIndex)[ed]);
        const Eigen::MatrixXd weakQuadratureWeights = weakReferenceSegment.Weights * absMapDeterminant;

        const Eigen::VectorXd neumannValues =
            test.weak_boundary_condition(boundary_info.Marker, region_id, weakQuadraturePoints, time_value);
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
                    rightHandSide.AddValue(local_dof_i.Global_Index, neumannContributions[p]);
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
                rightHandSide.AddValue(local_dof_i.Global_Index, neumannContributions[localIndex + 2]);
            }
            break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }
}
// ***************************************************************************
void Assembler::ComputeInitialCondition(const Polydim::examples::Parabolic_PCC_2D::Program_Configuration &config,
                                        const Gedim::IMeshDAO &mesh,
                                        const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                        const std::vector<unsigned int> &region_id,
                                        const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                        const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                        const test::I_Test &test,
                                        Gedim::Eigen_Array<> &initial_condition,
                                        Gedim::Eigen_Array<> &initial_condition_dirichlet) const
{

    initial_condition.SetSize(dofs_data.NumberDOFs);
    initial_condition_dirichlet.SetSize(dofs_data.NumberStrongs);

    // Assemble equation elements
    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
        // DOFs: vertices
        const Eigen::MatrixXd coordinates = mesh.Cell2DVerticesCoordinates(c);
        const Eigen::VectorXd dofs_vertices = test.initial_solution(coordinates);

        // Assemble local numerical solution
        unsigned int count = 0;
        for (unsigned int p = 0; p < mesh.Cell2DNumberVertices(c); p++)
        {
            const unsigned int cell0D_index = mesh.Cell2DVertex(c, p);

            const auto local_dofs = dofs_data.CellsDOFs.at(0).at(cell0D_index);
            for (unsigned int loc_i = 0; loc_i < local_dofs.size(); loc_i++)
            {
                const auto &local_dof_i = local_dofs.at(loc_i);
                const int global_i = local_dof_i.Global_Index;

                switch (local_dof_i.Type)
                {
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong: {
                    initial_condition_dirichlet.SetValue(global_i, dofs_vertices(count++));
                }
                break;
                case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF: {
                    initial_condition.SetValue(global_i, dofs_vertices(count++));
                }
                break;
                default:
                    throw std::runtime_error("Unknown DOF Type");
                }
            }
        }

        // Assemble strong boundary condition on Cell1Ds
        if (reference_element_data.Order > 1)
        {

            const auto local_space_data =
                Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(config.GeometricTolerance1D(),
                                                                       config.GeometricTolerance2D(),
                                                                       mesh_geometric_data,
                                                                       c,
                                                                       reference_element_data);

            // Assemble strong boundary condition on Cell1Ds
            for (unsigned int ed = 0; ed < mesh.Cell2DNumberEdges(c); ++ed)
            {
                const unsigned int cell1DIndex = mesh.Cell2DEdge(c, ed);

                const auto local_dofs = dofs_data.CellsDOFs.at(1).at(cell1DIndex);

                const auto edge_dofs_coordinates =
                    Polydim::PDETools::LocalSpace_PCC_2D::EdgeDofsCoordinates(reference_element_data, local_space_data, ed);

                const Eigen::VectorXd dofs_edge = test.initial_solution(edge_dofs_coordinates);

                for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
                {
                    const auto &local_dof_i = local_dofs.at(loc_i);
                    const int global_i = local_dof_i.Global_Index;

                    switch (local_dof_i.Type)
                    {
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong: {
                        initial_condition_dirichlet.SetValue(global_i, dofs_edge(loc_i));
                    }
                    break;
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF: {
                        initial_condition.SetValue(global_i, dofs_edge(loc_i));
                    }
                    break;
                    default:
                        throw std::runtime_error("Unknown DOF Type");
                    }
                }
            }

            const auto local_dofs = dofs_data.CellsDOFs.at(2).at(c);

            if (local_dofs.size())
            {
                const auto internal_dofs_coordinates =
                    Polydim::PDETools::LocalSpace_PCC_2D::InternalDofsCoordinates(reference_element_data, local_space_data);

                const Eigen::VectorXd initial_values_at_dofs = test.initial_solution(internal_dofs_coordinates.Points);

                const Eigen::VectorXd dofs_internal =
                    Polydim::PDETools::LocalSpace_PCC_2D::InternalDofs(reference_element_data,
                                                                       local_space_data,
                                                                       initial_values_at_dofs,
                                                                       internal_dofs_coordinates);

                for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
                {
                    const auto &local_dof_i = local_dofs.at(loc_i);
                    const int global_i = local_dof_i.Global_Index;

                    switch (local_dof_i.Type)
                    {
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong: {
                        initial_condition_dirichlet.SetValue(global_i, dofs_internal(loc_i));
                    }
                    break;
                    case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF: {
                        initial_condition.SetValue(global_i, dofs_internal(loc_i));
                    }
                    break;
                    default:
                        throw std::runtime_error("Unknown DOF Type");
                    }
                }
            }
        }
    }

    initial_condition.Create();
    initial_condition_dirichlet.Create();
}
// ***************************************************************************
void Assembler::Assemble(const Polydim::examples::Parabolic_PCC_2D::Program_Configuration &config,
                         const double &time_value,
                         const Gedim::MeshMatricesDAO &mesh,
                         const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                         const std::vector<unsigned int> &region_id,
                         const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                         const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                         const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                         const Polydim::examples::Parabolic_PCC_2D::test::I_Test &test,
                         Gedim::Eigen_SparseArray<> &globalMatrixA,
                         Gedim::Eigen_SparseArray<> &dirichletMatrixA,
                         Gedim::Eigen_Array<> &rightHandSide,
                         Gedim::Eigen_Array<> &solution_dirichlet) const
{
    globalMatrixA.SetSize(dofs_data.NumberDOFs, dofs_data.NumberDOFs, Gedim::ISparseArray::SparseArrayTypes::None);
    dirichletMatrixA.SetSize(dofs_data.NumberDOFs, dofs_data.NumberStrongs);
    rightHandSide.SetSize(dofs_data.NumberDOFs);
    solution_dirichlet.SetSize(dofs_data.NumberStrongs);

    Polydim::PDETools::Equations::EllipticEquation equation;

    Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data =
        {{std::cref(dofs_data)}, {0}, {0}, {0}};

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); ++c)
    {
        const auto local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(config.GeometricTolerance1D(),
                                                                                             config.GeometricTolerance2D(),
                                                                                             mesh_geometric_data,
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

        const auto diffusion_term_values = test.diffusion_term(region_id[c], cell2D_internal_quadrature.Points, time_value);
        const auto source_term_values = test.source_term(region_id[c], cell2D_internal_quadrature.Points, time_value);

        const Eigen::MatrixXd local_A = equation.ComputeCellDiffusionMatrix(diffusion_term_values,
                                                                            basis_functions_derivative_values,
                                                                            cell2D_internal_quadrature.Weights);

        Eigen::VectorXd local_rhs =
            equation.ComputeCellForcingTerm(source_term_values, basis_functions_values, cell2D_internal_quadrature.Weights);

        double k_max = diffusion_term_values[0].cwiseAbs().maxCoeff();
        for (unsigned int d = 1; d < diffusion_term_values.size(); d++)
        {
            const double local_k_max = diffusion_term_values[d].cwiseAbs().maxCoeff();
            if (local_k_max > k_max)
                k_max = local_k_max;
        }

        const Eigen::MatrixXd local_A_stab =
            k_max * Polydim::PDETools::LocalSpace_PCC_2D::StabilizationMatrix(reference_element_data, local_space_data);

        const auto &global_dofs = dofs_data.CellsGlobalDOFs[2].at(c);

        assert(Polydim::PDETools::LocalSpace_PCC_2D::Size(reference_element_data, local_space_data) == global_dofs.size());

        Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<2>(c,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_A + local_A_stab,
                                                                                          local_rhs,
                                                                                          globalMatrixA,
                                                                                          dirichletMatrixA,
                                                                                          rightHandSide);

        ComputeStrongTerm(time_value, c, region_id[c], mesh, mesh_dofs_info, dofs_data, reference_element_data, local_space_data, test, solution_dirichlet);
        ComputeWeakTerm(time_value, c, region_id[c], mesh, mesh_geometric_data, mesh_dofs_info, dofs_data, reference_element_data, local_space_data, test, rightHandSide);
    }

    rightHandSide.Create();
    solution_dirichlet.Create();
    globalMatrixA.Create();
    dirichletMatrixA.Create();

    if (dofs_data.NumberStrongs > 0)
        rightHandSide.SubtractionMultiplication(dirichletMatrixA, solution_dirichlet);
}
// ***************************************************************************
void Assembler::AssembleMassMatrix(const Polydim::examples::Parabolic_PCC_2D::Program_Configuration &config,
                                   const Gedim::MeshMatricesDAO &mesh,
                                   const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                   const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                                   const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                   const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                   Gedim::Eigen_SparseArray<> &globalMatrixM,
                                   Gedim::Eigen_SparseArray<> &dirichletMatrixM) const
{
    globalMatrixM.SetSize(dofs_data.NumberDOFs, dofs_data.NumberDOFs, Gedim::ISparseArray::SparseArrayTypes::None);
    dirichletMatrixM.SetSize(dofs_data.NumberDOFs, dofs_data.NumberStrongs);

    Polydim::PDETools::Assembler_Utilities::local_matrix_to_global_matrix_dofs_data local_matrix_to_global_matrix_dofs_data =
        {{std::cref(dofs_data)}, {0}, {0}, {0}};

    Polydim::PDETools::Equations::EllipticEquation equation;

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); ++c)
    {
        const auto local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(config.GeometricTolerance1D(),
                                                                                             config.GeometricTolerance2D(),
                                                                                             mesh_geometric_data,
                                                                                             c,
                                                                                             reference_element_data);

        const auto basis_functions_values =
            Polydim::PDETools::LocalSpace_PCC_2D::BasisFunctionsValues(reference_element_data,
                                                                       local_space_data,
                                                                       Polydim::VEM::PCC::ProjectionTypes::Pi0k);

        const auto cell2D_internal_quadrature =
            Polydim::PDETools::LocalSpace_PCC_2D::InternalQuadrature(reference_element_data, local_space_data);

        const Eigen::VectorXd reaction_term_values = Eigen::VectorXd::Ones(cell2D_internal_quadrature.Points.cols());
        Eigen::MatrixXd local_C =
            equation.ComputeCellReactionMatrix(reaction_term_values, basis_functions_values, cell2D_internal_quadrature.Weights);

        local_C += Polydim::PDETools::LocalSpace_PCC_2D::StabilizationMatrix(reference_element_data,
                                                                             local_space_data,
                                                                             Polydim::VEM::PCC::ProjectionTypes::Pi0k);

        const auto &global_dofs = dofs_data.CellsGlobalDOFs[2].at(c);

        assert(Polydim::PDETools::LocalSpace_PCC_2D::Size(reference_element_data, local_space_data) == global_dofs.size());
        Polydim::PDETools::Assembler_Utilities::assemble_local_matrix_to_global_matrix<2>(c,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_matrix_to_global_matrix_dofs_data,
                                                                                          local_C,
                                                                                          globalMatrixM,
                                                                                          dirichletMatrixM);
    }

    globalMatrixM.Create();
    dirichletMatrixM.Create();
}
// ***************************************************************************
Assembler::Performance_Data Assembler::ComputePerformance(const Polydim::examples::Parabolic_PCC_2D::Program_Configuration &config,
                                                          const Gedim::MeshMatricesDAO &mesh,
                                                          const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                                          const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data) const
{
    Assembler::Performance_Data result;
    result.Cell2DsPerformance.resize(mesh.Cell2DTotalNumber());

    // Assemble equation elements
    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
        const auto local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(config.GeometricTolerance1D(),
                                                                                             config.GeometricTolerance2D(),
                                                                                             mesh_geometric_data,
                                                                                             c,
                                                                                             reference_element_data);

        result.Cell2DsPerformance[c] =
            Polydim::PDETools::LocalSpace_PCC_2D::ComputePerformance(reference_element_data, local_space_data);
    }

    return result;
}
// ***************************************************************************
Assembler::PostProcess_Data Assembler::PostProcessSolution(const Polydim::examples::Parabolic_PCC_2D::Program_Configuration &config,
                                                           const double &time_value,
                                                           const Gedim::MeshMatricesDAO &mesh,
                                                           const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                                           const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                                           const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                                           const Parabolic_PCC_2D_Problem_Data &assembler_data,
                                                           const Polydim::examples::Parabolic_PCC_2D::test::I_Test &test) const
{
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
    result.cell0Ds_exact.setZero(mesh.Cell0DTotalNumber());

    for (unsigned int p = 0; p < mesh.Cell0DTotalNumber(); p++)
    {
        result.cell0Ds_exact[p] = test.exact_solution(mesh.Cell0DCoordinates(p), time_value)[0];

        const auto local_dofs = dofs_data.CellsDOFs.at(0).at(p);

        for (unsigned int loc_i = 0; loc_i < local_dofs.size(); ++loc_i)
        {
            const auto &local_dof_i = local_dofs.at(loc_i);

            switch (local_dof_i.Type)
            {
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::Strong:
                result.cell0Ds_numeric[p] = assembler_data.solution_dirichlet.GetValue(local_dof_i.Global_Index);
                break;
            case Polydim::PDETools::DOFs::DOFsManager::DOFsData::DOF::Types::DOF:
                result.cell0Ds_numeric[p] = assembler_data.solution.GetValue(local_dof_i.Global_Index);
                break;
            default:
                throw std::runtime_error("Unknown DOF Type");
            }
        }
    }

    result.cell2Ds_error_L2.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_norm_L2.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_error_H1.setZero(mesh.Cell2DTotalNumber());
    result.cell2Ds_norm_H1.setZero(mesh.Cell2DTotalNumber());
    result.mesh_size = 0.0;

    for (unsigned int c = 0; c < mesh.Cell2DTotalNumber(); c++)
    {
        const auto local_space_data = Polydim::PDETools::LocalSpace_PCC_2D::CreateLocalSpace(config.GeometricTolerance1D(),
                                                                                             config.GeometricTolerance2D(),
                                                                                             mesh_geometric_data,
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

        const auto exact_solution_values = test.exact_solution(cell2D_internal_quadrature.Points, time_value);
        const auto exact_derivative_solution_values = test.exact_derivative_solution(cell2D_internal_quadrature.Points, time_value);

        const auto local_count_dofs = Polydim::PDETools::Assembler_Utilities::local_count_dofs<2>(c, dofs_data);
        const Eigen::VectorXd dofs_values =
            PDETools::Assembler_Utilities::global_solution_to_local_solution<2>(c,
                                                                                dofs_data,
                                                                                local_count_dofs.num_total_dofs,
                                                                                local_count_dofs.offsets_DOFs,
                                                                                {0},
                                                                                {0},
                                                                                assembler_data.solution,
                                                                                assembler_data.solution_dirichlet);

        const Eigen::VectorXd local_error_L2 = (basis_functions_values * dofs_values - exact_solution_values).array().square();
        const Eigen::VectorXd local_norm_L2 = (basis_functions_values * dofs_values).array().square();

        result.cell2Ds_error_L2[c] = cell2D_internal_quadrature.Weights.transpose() * local_error_L2;
        result.cell2Ds_norm_L2[c] = cell2D_internal_quadrature.Weights.transpose() * local_norm_L2;

        const Eigen::VectorXd local_error_H1 =
            (basis_functions_derivative_values[0] * dofs_values - exact_derivative_solution_values[0]).array().square() +
            (basis_functions_derivative_values[1] * dofs_values - exact_derivative_solution_values[1]).array().square();

        const Eigen::VectorXd local_norm_H1 = (basis_functions_derivative_values[0] * dofs_values).array().square() +
                                              (basis_functions_derivative_values[1] * dofs_values).array().square();

        result.cell2Ds_error_H1[c] = cell2D_internal_quadrature.Weights.transpose() * local_error_H1;
        result.cell2Ds_norm_H1[c] = cell2D_internal_quadrature.Weights.transpose() * local_norm_H1;

        if (mesh_geometric_data.Cell2DsDiameters.at(c) > result.mesh_size)
            result.mesh_size = mesh_geometric_data.Cell2DsDiameters.at(c);
    }

    result.error_L2 = std::sqrt(result.cell2Ds_error_L2.sum());
    result.norm_L2 = std::sqrt(result.cell2Ds_norm_L2.sum());
    result.error_H1 = std::sqrt(result.cell2Ds_error_H1.sum());
    result.norm_H1 = std::sqrt(result.cell2Ds_norm_H1.sum());

    return result;
}
} // namespace Parabolic_PCC_2D
} // namespace examples
} // namespace Polydim
