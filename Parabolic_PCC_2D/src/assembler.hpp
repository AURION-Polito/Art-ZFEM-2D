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

#ifndef __assembler_H
#define __assembler_H

#include "Assembler_Utilities.hpp"
#include "DOFsManager.hpp"
#include "Eigen_Array.hpp"
#include "Eigen_SparseArray.hpp"
#include "MeshUtilities.hpp"
#include "program_configuration.hpp"
#include "test_definition.hpp"

namespace Polydim
{
namespace examples
{
namespace Parabolic_PCC_2D
{
class Assembler final
{
  public:
    struct Parabolic_PCC_2D_Problem_Data final
    {
        Gedim::Eigen_SparseArray<> globalMatrixM;
        Gedim::Eigen_SparseArray<> dirichletMatrixM;

        Gedim::Eigen_SparseArray<> globalMatrixA;
        Gedim::Eigen_SparseArray<> dirichletMatrixA;
        Gedim::Eigen_Array<> rightHandSide;

        Gedim::Eigen_SparseArray<> initial_globalMatrixA;
        Gedim::Eigen_SparseArray<> initial_dirichletMatrixA;
        Gedim::Eigen_Array<> initial_rightHandSide;

        Gedim::Eigen_Array<> solution;
        Gedim::Eigen_Array<> solution_dirichlet;

        Gedim::Eigen_Array<> initial_solution;
        Gedim::Eigen_Array<> initial_solution_dirichlet;
    };

    struct Performance_Data final
    {
        std::vector<Polydim::PDETools::LocalSpace_PCC_2D::Performance_Data> Cell2DsPerformance;
    };

    struct PostProcess_Data final
    {
        Eigen::VectorXd cell0Ds_numeric;
        Eigen::VectorXd cell0Ds_exact;

        Eigen::VectorXd cell2Ds_error_L2;
        Eigen::VectorXd cell2Ds_norm_L2;
        double error_L2;
        double norm_L2;
        Eigen::VectorXd cell2Ds_error_H1;
        Eigen::VectorXd cell2Ds_norm_H1;
        double error_H1;
        double norm_H1;

        double mesh_size;

        double residual_norm;
        double conditioning;
    };

  private:
    void ComputeStrongTerm(const double &time_value,
                           const unsigned int cell2DIndex,
                           const unsigned int region_id,
                           const Gedim::MeshMatricesDAO &mesh,
                           const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                           const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                           const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                           const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data,
                           const test::I_Test &test,
                           Gedim::Eigen_Array<> &solution_dirichlet) const;

    void ComputeWeakTerm(const double &time_value,
                         const unsigned int cell2DIndex,
                         const unsigned int region_id,
                         const Gedim::MeshMatricesDAO &mesh,
                         const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                         const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                         const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                         const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                         const Polydim::PDETools::LocalSpace_PCC_2D::LocalSpace_Data &local_space_data,
                         const Polydim::examples::Parabolic_PCC_2D::test::I_Test &test,
                         Gedim::Eigen_Array<> &rightHandSide) const;

  public:
    void ComputeInitialCondition(const Polydim::examples::Parabolic_PCC_2D::Program_Configuration &config,
                                 const Gedim::IMeshDAO &mesh,
                                 const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                 const std::vector<unsigned int> &region_id,
                                 const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                 const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                 const test::I_Test &test,
                                 Gedim::Eigen_Array<> &initial_solution,
                                 Gedim::Eigen_Array<> &initial_solution_dirichlet) const;

    void Assemble(const Polydim::examples::Parabolic_PCC_2D::Program_Configuration &config,
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
                  Gedim::Eigen_Array<> &solution_dirichlet) const;

    void AssembleMassMatrix(const Polydim::examples::Parabolic_PCC_2D::Program_Configuration &config,
                            const Gedim::MeshMatricesDAO &mesh,
                            const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                            const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                            const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                            const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                            Gedim::Eigen_SparseArray<> &globalMatrixM,
                            Gedim::Eigen_SparseArray<> &dirichletMatrixM) const;

    Performance_Data ComputePerformance(const Polydim::examples::Parabolic_PCC_2D::Program_Configuration &config,
                                        const Gedim::MeshMatricesDAO &mesh,
                                        const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                        const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data) const;

    PostProcess_Data PostProcessSolution(const Polydim::examples::Parabolic_PCC_2D::Program_Configuration &config,
                                         const double &time_value,
                                         const Gedim::MeshMatricesDAO &mesh,
                                         const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                                         const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                                         const Polydim::PDETools::LocalSpace_PCC_2D::ReferenceElement_Data &reference_element_data,
                                         const Parabolic_PCC_2D_Problem_Data &assembler_data,
                                         const Polydim::examples::Parabolic_PCC_2D::test::I_Test &test) const;
};
} // namespace Parabolic_PCC_2D
} // namespace examples
} // namespace Polydim

#endif
