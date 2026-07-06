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

#ifndef __test_definition_H
#define __test_definition_H

#include "DOFsManager.hpp"
#include "PDE_Mesh_Utilities.hpp"

namespace Polydim
{
namespace examples
{
namespace Parabolic_PCC_2D
{
namespace test
{
enum struct Test_Types
{
    Patch_Test = 1,
    Covergence_in_space = 2,
    Covergence_in_time = 3,
    Benchmark_problem_1 = 4, /// Test 1 in Berrone 2006: "Robust a posteriori error estimates for finite element
                             /// discretizations of the heat equation with discontinuous coefficients
    Benchmark_problem_2 = 5  /// Test 2 in Berrone 2006: "Robust a posteriori error estimates for finite element
                             /// discretizations of the heat equation with discontinuous coefficients
};

struct I_Test
{
    virtual Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Time_Domain_2D domain() const = 0;
    virtual std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info() const = 0;
    virtual std::array<Eigen::VectorXd, 9> diffusion_term(const unsigned int region_id,
                                                          const Eigen::MatrixXd &points,
                                                          const double &time_value) const = 0;
    virtual Eigen::VectorXd source_term(const unsigned int region_id, const Eigen::MatrixXd &points, const double &time_value) const = 0;
    virtual Eigen::VectorXd strong_boundary_condition(const unsigned int marker,
                                                      const Eigen::MatrixXd &points,
                                                      const double &time_value) const = 0;
    virtual Eigen::VectorXd weak_boundary_condition(const unsigned int marker,
                                                    const unsigned int region_id,
                                                    const Eigen::MatrixXd &points,
                                                    const double &time_value) const = 0;
    virtual Eigen::VectorXd exact_solution(const Eigen::MatrixXd &points, const double &time_value) const = 0;
    virtual Eigen::VectorXd initial_solution(const Eigen::MatrixXd &points) const = 0;
    virtual std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd &points,
                                                                     const double &time_value) const = 0;
};
// ***************************************************************************
struct Patch_Test final : public I_Test
{
    static unsigned int space_order;
    static unsigned int time_order;

    Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Time_Domain_2D domain() const
    {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Time_Domain_2D domain;

        domain.spatial_domain.area = 1.0;

        domain.spatial_domain.vertices = Eigen::MatrixXd::Zero(3, 4);
        domain.spatial_domain.vertices.row(0) << 0.0, 1.0, 1.0, 0.0;
        domain.spatial_domain.vertices.row(1) << 0.0, 0.0, 1.0, 1.0;

        domain.spatial_domain.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types::Parallelogram;

        domain.time_domain = {0.0, 1.0};

        return domain;
    }

    std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info() const
    {
        return {{0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}}};
    }

    std::array<Eigen::VectorXd, 9> diffusion_term(const unsigned int region_id, const Eigen::MatrixXd &points, const double &time_value) const
    {
        return {Eigen::VectorXd::Constant(points.cols(), 1.0),
                Eigen::VectorXd::Constant(points.cols(), 0.0),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 0.0),
                Eigen::VectorXd::Constant(points.cols(), 1.0),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 0.0)};
    }

    Eigen::VectorXd source_term(const unsigned int region_id, const Eigen::MatrixXd &points, const double &time_value) const
    {
        Eigen::ArrayXd source_term = Eigen::VectorXd::Constant(points.cols(), 2.0 * space_order * (space_order - 1));
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        const int max_order = space_order - 2;
        for (int i = 0; i < max_order; ++i)
            source_term *= polynomial;

        if (time_order == 1)
            source_term -= 1.0;
        else if (time_order == 2)
            source_term -= 2.0 * time_value;
        else
            throw std::runtime_error("not valid time order");

        return -source_term;
    }

    Eigen::VectorXd strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points, const double &time_value) const
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        return exact_solution(points, time_value);
    }

    Eigen::VectorXd weak_boundary_condition(const unsigned int marker,
                                            const unsigned int region_id,
                                            const Eigen::MatrixXd &points,
                                            const double &time_value) const
    {
        throw std::runtime_error("not valid marker");
    }

    Eigen::VectorXd exact_solution(const Eigen::MatrixXd &points, const double &time_value) const
    {
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        Eigen::ArrayXd result = Eigen::VectorXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < space_order; ++i)
            result *= polynomial;

        if (time_order == 1)
            result += time_value;
        else if (time_order == 2)
            result += time_value * time_value;
        else
            throw std::runtime_error("not valid time order");

        return result;
    }

    Eigen::VectorXd initial_solution(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        Eigen::ArrayXd result = Eigen::VectorXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < space_order; ++i)
            result *= polynomial;

        return result;
    }

    std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd &points, const double &time_value) const
    {
        Eigen::VectorXd derivatives = Eigen::VectorXd::Constant(points.cols(), space_order);
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        const int max_order = space_order - 1;
        for (int i = 0; i < max_order; ++i)
            derivatives.array() *= polynomial;

        return {derivatives, derivatives, Eigen::VectorXd::Zero(points.cols())};
    }
};
// ***************************************************************************

} // namespace test
} // namespace Parabolic_PCC_2D
} // namespace examples
} // namespace Polydim

#endif
