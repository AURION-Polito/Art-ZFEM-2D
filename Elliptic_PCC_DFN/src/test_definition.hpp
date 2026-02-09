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
#include "Extended_PDE_Mesh_Utilities.hpp"
#include <numbers>

namespace Polydim
{
namespace examples
{
namespace Elliptic_PCC_DFN
{
namespace test
{
enum struct Test_Types
{
    Patch_Test = 1,
    Elliptic_Polynomial_Problem = 2, /// Test 1: S. Berrone, G. Teora, F. Vicini, "Improving high-order VEM stability on
    /// badly-shaped elements", doi: https://doi.org/10.1016/j.matcom.2023.10.003.
    Elliptic_Problem = 3,
    Patch_Test_Rotated = 4,
    DFN_Frac_3 = 5
};

struct I_Test
{
    virtual PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D_Collection domains() const = 0;
    virtual std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info() const = 0;
    virtual Eigen::VectorXd diffusion_term(const unsigned int domain_position, const Eigen::MatrixXd &points) const = 0;
    virtual Eigen::VectorXd source_term(const unsigned int domain_position, const Eigen::MatrixXd &points) const = 0;
    virtual Eigen::VectorXd strong_boundary_condition(const unsigned int domain_position,
                                                      const unsigned int marker,
                                                      const Eigen::MatrixXd &points) const = 0;
    virtual Eigen::VectorXd weak_boundary_condition(const unsigned int domain_position,
                                                    const unsigned int marker,
                                                    const Eigen::MatrixXd &points) const = 0;
    virtual Eigen::VectorXd exact_solution(const unsigned int domain_position, const Eigen::MatrixXd &points) const = 0;
    virtual std::array<Eigen::VectorXd, 3> exact_derivative_solution(const unsigned int domain_position,
                                                                     const Eigen::MatrixXd &points) const = 0;
};
// ***************************************************************************
struct Patch_Test final : public I_Test
{
    static unsigned int order;

    PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D_Collection domains() const
    {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain;

        domain.area = 1.0;

        domain.vertices = Eigen::MatrixXd::Zero(3, 4);
        domain.vertices.row(0) << 0.0, 1.0, 1.0, 0.0;
        domain.vertices.row(1) << 0.0, 0.0, 1.0, 1.0;

        domain.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types::Parallelogram;

        return {{domain}, {Eigen::Matrix3d::Identity()}, {Eigen::Vector3d::Zero()}};
    }

    std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info() const
    {
        return {{0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 4}},
                {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}}};
    }

    Eigen::VectorXd diffusion_term(const unsigned int domain_position, const Eigen::MatrixXd &points) const
    {
        if (domain_position != 0)
            throw std::runtime_error("Unknown domain position");

        return Eigen::VectorXd::Constant(points.cols(), 1.0);
    };

    Eigen::VectorXd source_term(const unsigned int domain_position, const Eigen::MatrixXd &points) const
    {
        if (domain_position != 0)
            throw std::runtime_error("Unknown domain position");

        Eigen::VectorXd source_term = Eigen::VectorXd::Constant(points.cols(), 2.0 * order * (order - 1));
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        const int max_order = order - 2;
        for (int i = 0; i < max_order; ++i)
            source_term.array() *= polynomial;

        return -source_term;
    };

    Eigen::VectorXd strong_boundary_condition(const unsigned int domain_position, const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (domain_position != 0)
            throw std::runtime_error("Unknown domain position");

        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        Eigen::VectorXd result = Eigen::VectorXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order; ++i)
            result.array() *= polynomial;

        return result;
    };

    Eigen::VectorXd weak_boundary_condition(const unsigned int domain_position, const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (domain_position != 0)
            throw std::runtime_error("Unknown domain position");

        Eigen::VectorXd derivatives = Eigen::VectorXd::Constant(points.cols(), 1.0);
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        const int max_order = order - 1;
        for (int i = 0; i < max_order; ++i)
            derivatives.array() *= polynomial;

        std::array<Eigen::VectorXd, 3> der = {derivatives, derivatives, Eigen::VectorXd::Zero(points.cols())};

        switch (marker)
        {
        case 2:
            return order * derivatives.array();
        case 4:
            return order * derivatives.array();
        default:
            throw std::runtime_error("not valid marker");
        }

        throw std::runtime_error("Not supported");
    }

    Eigen::VectorXd exact_solution(const unsigned int domain_position, const Eigen::MatrixXd &points) const
    {
        if (domain_position != 0)
            throw std::runtime_error("Unknown domain position");

        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        Eigen::VectorXd result = Eigen::VectorXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order; ++i)
            result.array() *= polynomial;

        return result;
    };

    std::array<Eigen::VectorXd, 3> exact_derivative_solution(const unsigned int domain_position, const Eigen::MatrixXd &points) const
    {
        if (domain_position != 0)
            throw std::runtime_error("Unknown domain position");

        Eigen::VectorXd derivatives = Eigen::VectorXd::Constant(points.cols(), order);
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        const int max_order = order - 1;
        for (int i = 0; i < max_order; ++i)
            derivatives.array() *= polynomial;

        return {derivatives, derivatives, Eigen::VectorXd::Zero(points.cols())};
    }
};
// ***************************************************************************
struct Elliptic_Polynomial_Problem final : public I_Test
{
    PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D_Collection domains() const
    {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain;

        domain.area = 1.0;

        domain.vertices = Eigen::MatrixXd::Zero(3, 4);
        domain.vertices.row(0) << 0.0, 1.0, 1.0, 0.0;
        domain.vertices.row(1) << 0.0, 0.0, 1.0, 1.0;

        domain.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types::Parallelogram;

        return {{domain}, {Eigen::Matrix3d::Identity()}, {Eigen::Vector3d::Zero()}};
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

    Eigen::VectorXd diffusion_term(const unsigned int domain_position, const Eigen::MatrixXd &points) const
    {
        if (domain_position != 0)
            throw std::runtime_error("Unknown domain position");

        const double k = 1.0;
        return Eigen::VectorXd::Constant(points.cols(), k);
    };

    Eigen::VectorXd source_term(const unsigned int domain_position, const Eigen::MatrixXd &points) const
    {
        if (domain_position != 0)
            throw std::runtime_error("Unknown domain position");

        return 32.0 * (points.row(1).array() * (1.0 - points.row(1).array()) +
                       points.row(0).array() * (1.0 - points.row(0).array()));
    };

    Eigen::VectorXd strong_boundary_condition(const unsigned int domain_position, const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (domain_position != 0)
            throw std::runtime_error("Unknown domain position");

        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        return 16.0 * (points.row(1).array() * (1.0 - points.row(1).array()) * points.row(0).array() *
                       (1.0 - points.row(0).array())) +
               1.1;
    };

    Eigen::VectorXd weak_boundary_condition(const unsigned int domain_position, const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (domain_position != 0)
            throw std::runtime_error("Unknown domain position");

        switch (marker)
        {
        case 2: // co-normal derivatives on the right
            return 16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array());
        case 4: // co-normal derivatives on the left
            return -16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array());
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    Eigen::VectorXd exact_solution(const unsigned int domain_position, const Eigen::MatrixXd &points) const
    {
        if (domain_position != 0)
            throw std::runtime_error("Unknown domain position");

        return 16.0 * (points.row(1).array() * (1.0 - points.row(1).array()) * points.row(0).array() *
                       (1.0 - points.row(0).array())) +
               1.1;
    };

    std::array<Eigen::VectorXd, 3> exact_derivative_solution(const unsigned int domain_position, const Eigen::MatrixXd &points) const
    {
        if (domain_position != 0)
            throw std::runtime_error("Unknown domain position");

        return {16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array()),
                16.0 * (1.0 - 2.0 * points.row(1).array()) * points.row(0).array() * (1.0 - points.row(0).array()),
                Eigen::VectorXd::Zero(points.cols())};
    }
};
// ***************************************************************************
struct Elliptic_Problem final : public I_Test
{
    PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D_Collection domains() const
    {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain;

        domain.area = 1.0;

        domain.vertices = Eigen::MatrixXd::Zero(3, 4);
        domain.vertices.row(0) << 0.0, 1.0, 1.0, 0.0;
        domain.vertices.row(1) << 0.0, 0.0, 1.0, 1.0;

        domain.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types::Parallelogram;

        return {{domain}, {Eigen::Matrix3d::Identity()}, {Eigen::Vector3d::Zero()}};
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

    Eigen::VectorXd diffusion_term(const unsigned int domain_position, const Eigen::MatrixXd &points) const
    {
        if (domain_position != 0)
            throw std::runtime_error("Unknown domain position");

        const double k = 2.0;
        return Eigen::VectorXd::Constant(points.cols(), k);
    };

    Eigen::VectorXd source_term(const unsigned int domain_position, const Eigen::MatrixXd &points) const
    {
        if (domain_position != 0)
            throw std::runtime_error("Unknown domain position");

        return 16.0 * std::numbers::pi * std::numbers::pi * sin(2.0 * std::numbers::pi * points.row(0).array()) *
               sin(2.0 * std::numbers::pi * points.row(1).array());
    };

    Eigen::VectorXd strong_boundary_condition(const unsigned int domain_position, const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (domain_position != 0)
            throw std::runtime_error("Unknown domain position");

        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        return exact_solution(domain_position, points);
    };

    Eigen::VectorXd weak_boundary_condition(const unsigned int domain_position, const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (domain_position != 0)
            throw std::runtime_error("Unknown domain position");

        switch (marker)
        {
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    Eigen::VectorXd exact_solution(const unsigned int domain_position, const Eigen::MatrixXd &points) const
    {
        if (domain_position != 0)
            throw std::runtime_error("Unknown domain position");

        return sin(2.0 * std::numbers::pi * points.row(0).array()) * sin(2.0 * std::numbers::pi * points.row(1).array());
    };

    std::array<Eigen::VectorXd, 3> exact_derivative_solution(const unsigned int domain_position, const Eigen::MatrixXd &points) const
    {
        if (domain_position != 0)
            throw std::runtime_error("Unknown domain position");

        return {2.0 * std::numbers::pi * cos(2.0 * std::numbers::pi * points.row(0).array()) *
                    sin(2.0 * std::numbers::pi * points.row(1).array()),
                2.0 * std::numbers::pi * sin(2.0 * std::numbers::pi * points.row(0).array()) *
                    cos(2.0 * std::numbers::pi * points.row(1).array()),
                Eigen::VectorXd::Zero(points.cols())};
    }
};
// ***************************************************************************
struct Patch_Test_Rotated final : public I_Test
{
  private:
    double tol_1D;
    double tol_2D;
    unsigned int order;

  public:
    Patch_Test_Rotated(const double &tol_1D, const double &tol_2D, const unsigned int order)
        : tol_1D(tol_1D), tol_2D(tol_2D), order(order)
    {
    }

    PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D_Collection domains() const
    {
        Gedim::GeometryUtilitiesConfig geometry_utilities_config;
        geometry_utilities_config.Tolerance1D = tol_1D;
        geometry_utilities_config.Tolerance2D = tol_2D;
        Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

        Eigen::MatrixXd vertices_3D(3, 4);
        vertices_3D.row(0) << 0.0, 0.0, 0.0, 0.0;
        vertices_3D.row(1) << 0.0, 1.0, 1.0, 0.0;
        vertices_3D.row(2) << 0.0, 0.0, 1.0, 1.0;

        const auto domain_normal = geometry_utilities.PolygonNormal(vertices_3D);
        const auto domain_translation = geometry_utilities.PolygonTranslation(vertices_3D);
        const auto domain_rotation = geometry_utilities.PolygonRotationMatrix(vertices_3D, domain_normal, domain_translation);

        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain;

        domain.area = 1.0;

        domain.vertices = geometry_utilities.RotatePointsFrom3DTo2D(vertices_3D, domain_rotation.transpose(), domain_translation);

        domain.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types::Parallelogram;

        return {{domain}, {domain_rotation}, {domain_translation}};
    }

    std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info() const
    {
        return {{0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 2}},
                {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 4}},
                {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}}};
    }

    Eigen::VectorXd diffusion_term(const unsigned int domain_position, const Eigen::MatrixXd &points) const
    {
        if (domain_position != 0)
            throw std::runtime_error("Unknown domain position");

        return Eigen::VectorXd::Constant(points.cols(), 1.0);
    };

    Eigen::VectorXd source_term(const unsigned int domain_position, const Eigen::MatrixXd &points) const
    {
        if (domain_position != 0)
            throw std::runtime_error("Unknown domain position");

        if (order <= 1)
            return Eigen::VectorXd::Zero(points.cols());

        Eigen::VectorXd source_term = Eigen::VectorXd::Constant(points.cols(), 2.0 * order * (order - 1));
        const Eigen::ArrayXd polynomial = points.row(1).array() + points.row(2).array() + 0.5;

        const int max_order = order - 2;
        for (int i = 0; i < max_order; ++i)
            source_term.array() *= polynomial;

        return -source_term;
    };

    Eigen::VectorXd strong_boundary_condition(const unsigned int domain_position, const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (domain_position != 0)
            throw std::runtime_error("Unknown domain position");

        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        return exact_solution(domain_position, points);
    };

    Eigen::VectorXd weak_boundary_condition(const unsigned int domain_position, const unsigned int marker, const Eigen::MatrixXd &points) const
    {

        if (domain_position != 0)
            throw std::runtime_error("Unknown domain position");

        const auto derivatives = exact_derivative_solution(domain_position, points);

        switch (marker)
        {
        case 2:
            return derivatives.at(1).array();
        case 4:
            return derivatives.at(2).array();
        default:
            throw std::runtime_error("not valid marker");
        }

        throw std::runtime_error("Not supported");
    }

    Eigen::VectorXd exact_solution(const unsigned int domain_position, const Eigen::MatrixXd &points) const
    {
        if (domain_position != 0)
            throw std::runtime_error("Unknown domain position");

        const Eigen::ArrayXd polynomial = points.row(1).array() + points.row(2).array() + 0.5;

        Eigen::VectorXd result = Eigen::VectorXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order; ++i)
            result.array() *= polynomial;

        return result;
    };

    std::array<Eigen::VectorXd, 3> exact_derivative_solution(const unsigned int domain_position, const Eigen::MatrixXd &points) const
    {
        if (domain_position != 0)
            throw std::runtime_error("Unknown domain position");

        Eigen::VectorXd derivatives = Eigen::VectorXd::Constant(points.cols(), order);
        const Eigen::ArrayXd polynomial = points.row(1).array() + points.row(2).array() + 0.5;

        if (order > 0)
        {
            const int max_order = order - 1;
            for (int i = 0; i < max_order; ++i)
                derivatives.array() *= polynomial;
        }

        return {Eigen::VectorXd::Zero(points.cols()), derivatives, derivatives};
    }
};

// ***************************************************************************
struct DFN_Frac_3 final : public I_Test
{
  private:
    double tol_1D;
    double tol_2D;

  public:
    DFN_Frac_3(const double &tol_1D, const double &tol_2D) : tol_1D(tol_1D), tol_2D(tol_2D)
    {
    }

    PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D_Collection domains() const
    {
        Gedim::GeometryUtilitiesConfig geometry_utilities_config;
        geometry_utilities_config.Tolerance1D = tol_1D;
        geometry_utilities_config.Tolerance2D = tol_2D;
        Gedim::GeometryUtilities geometry_utilities(geometry_utilities_config);

        PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D_Collection domains;
        domains.domains_2D.resize(3);
        domains.domains_rotation.resize(3);
        domains.domains_translation.resize(3);

        std::vector<Eigen::MatrixXd> vertices_3D(3);

        vertices_3D.at(0).resize(3, 4);
        vertices_3D.at(0).col(0) << -1.0000, -1.0000, 0.0000;
        vertices_3D.at(0).col(1) << 0.5000, -1.0000, 0.0000;
        vertices_3D.at(0).col(2) << 0.5000, 1.0000, 0.0000;
        vertices_3D.at(0).col(3) << -1.0000, 1.0000, 0.0000;

        vertices_3D.at(1).resize(3, 4);
        vertices_3D.at(1).col(0) << -1.0000, 0.0000, -1.0000;
        vertices_3D.at(1).col(1) << 0.0000, 0.0000, -1.0000;
        vertices_3D.at(1).col(2) << 0.0000, 0.0000, 1.0000;
        vertices_3D.at(1).col(3) << -1.0000, 0.0000, 1.0000;

        vertices_3D.at(2).resize(3, 4);
        vertices_3D.at(2).col(0) << -0.5, -1.0000, -1.0000;
        vertices_3D.at(2).col(1) << -0.5, 1.0000, -1.0000;
        vertices_3D.at(2).col(2) << -0.5, 1.0000, 1.0000;
        vertices_3D.at(2).col(3) << -0.5, -1.0000, 1.0000;

        for (unsigned int d = 0; d < 3; ++d)
        {
            const auto &domain_vertices_3D = vertices_3D.at(d);

            auto &domain_2D = domains.domains_2D.at(d);
            auto &domain_translation = domains.domains_translation.at(d);
            auto &domain_rotation = domains.domains_rotation.at(d);

            const auto domain_normal = geometry_utilities.PolygonNormal(domain_vertices_3D);
            domain_translation = geometry_utilities.PolygonTranslation(domain_vertices_3D);
            domain_rotation = geometry_utilities.PolygonRotationMatrix(domain_vertices_3D, domain_normal, domain_translation);

            domain_2D.vertices =
                geometry_utilities.RotatePointsFrom3DTo2D(domain_vertices_3D, domain_rotation.transpose(), domain_translation);
            domain_2D.area = geometry_utilities.PolygonArea(domain_2D.vertices);

            domain_2D.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types::Polygon;
        }

        return domains;
    }

    std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info() const
    {
        return {{0, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 0}},
                {1, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {2, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {3, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::None, 1}},
                {4, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {5, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}},
                {6, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 1}},
                {7, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Weak, 1}},
                {8, {Polydim::PDETools::DOFs::DOFsManager::BoundaryTypes::Strong, 1}}};
    }

    Eigen::VectorXd diffusion_term(const unsigned int domain_position, const Eigen::MatrixXd &points) const
    {
        if (domain_position > 2)
            throw std::runtime_error("Unknown domain position");

        return Eigen::VectorXd::Constant(points.cols(), 1.0);
    };

    Eigen::VectorXd source_term(const unsigned int domain_position, const Eigen::MatrixXd &points) const
    {
        if (domain_position > 2)
            throw std::runtime_error("Unknown domain position");

        Eigen::VectorXd result = Eigen::VectorXd::Zero(points.cols());

        const auto &x = points.row(0);
        const auto &y = points.row(1);
        const auto &z = points.row(2);

        switch (domain_position)
        {
        case 0:
            for (unsigned int p = 0; p < points.cols(); p++)
            {
                result[p] =
                    0.1 * (-0.5 - x(p)) *
                        (6.0 * x(p) - 16.0 * y(p) * y(p) - (16.0 * x(p) * x(p) * y(p) * y(p)) / (x(p) * x(p) + y(p) * y(p)) +
                         48.0 * x(p) * y(p) * atan2(y(p), x(p))) +
                    0.1 * (-0.5 - x(p)) *
                        (16.0 * x(p) * x(p) + (16.0 * x(p) * x(p) * y(p) * y(p)) / (x(p) * x(p) + y(p) * y(p)) +
                         48.0 * x(p) * y(p) * atan2(y(p), x(p))) +
                    0.2 * (-3.0 * x(p) * x(p) + 8.0 * x(p) * y(p) * y(p) - 16.0 * x(p) * x(p) * y(p) * atan2(y(p), x(p)) -
                           8.0 * y(p) * (x(p) * x(p) + y(p) * y(p)) * atan2(y(p), x(p)));
            }
            break;
        case 1:
            for (unsigned int p = 0; p < points.cols(); p++)
            {
                result[p] = 0.6 * (-0.5 - x(p)) * x(p) - (3.0 * x(p) * x(p)) / 5.0 -
                            4.8 * M_PI * (-0.5 - x(p)) * x(p) * abs(z(p)) + 4.8 * M_PI * x(p) * x(p) * abs(z(p));
            }
            break;
        case 2:
            for (unsigned int p = 0; p < points.cols(); p++)
            {
                result[p] = 2.0 * (-1.0 + y(p)) * y(p) * (1.0 + y(p)) + 2.0 * (-1.0 + y(p)) * (-1.0 + z(p)) * z(p) +
                            2.0 * y(p) * (-1.0 + z(p)) * z(p) + 2.0 * (1.0 + y(p)) * (-1.0 + z(p)) * z(p);
            }
            break;
        default:
            throw std::runtime_error("Domain index " + std::to_string(domain_position) + " not supported");
        }

        return -1.0 * result;
    };

    Eigen::VectorXd strong_boundary_condition(const unsigned int domain_position, const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (domain_position > 2)
            throw std::runtime_error("Unknown domain position");

        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        return exact_solution(domain_position, points);
    };

    Eigen::VectorXd weak_boundary_condition(const unsigned int domain_position, const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (domain_position > 2)
            throw std::runtime_error("Unknown domain position");

        switch (marker)
        {
        default:
            throw std::runtime_error("not valid marker");
        }

        throw std::runtime_error("Not supported");
    }

    Eigen::VectorXd exact_solution(const unsigned int domain_position, const Eigen::MatrixXd &points) const
    {
        if (domain_position > 2)
            throw std::runtime_error("Unknown domain position");

        Eigen::VectorXd result = Eigen::VectorXd::Zero(points.cols());

        const auto &x = points.row(0);
        const auto &y = points.row(1);
        const auto &z = points.row(2);

        switch (domain_position)
        {
        case 0:
            // ((-10^(-1)).*(1/2 + x)).*(x.^3 + 8.*x.*y.*(x.^2 + y.^2).*atan2(y,x))
            for (unsigned int p = 0; p < points.cols(); p++)
            {
                result[p] = -0.1 * (0.5 + x(p)) *
                            (x(p) * x(p) * x(p) + 8.0 * x(p) * y(p) * (x(p) * x(p) + y(p) * y(p)) * atan2(y(p), x(p)));
            }
            break;
        case 1:
            // (1/10).*(-(1/2) - x).*x.^3 - abs(z).*pi.*(4/5).*(-(1/2) - x).*x.^3
            for (unsigned int p = 0; p < points.cols(); p++)
            {
                result[p] =
                    0.1 * (-0.5 - x(p)) * x(p) * x(p) * x(p) - abs(z(p)) * M_PI * 0.8 * (-0.5 - x(p)) * x(p) * x(p) * x(p);
            }
            break;
        case 2:
            // y.*(y - 1).*(y + 1).*z.*(z - 1)
            for (unsigned int p = 0; p < points.cols(); p++)
            {
                result[p] = y(p) * (y(p) - 1.0) * (y(p) + 1.0) * z(p) * (z(p) - 1.0);
            }
            break;
        default:
            throw std::runtime_error("Domain index " + std::to_string(domain_position) + " not supported");
        }

        return result;
    };

    std::array<Eigen::VectorXd, 3> exact_derivative_solution(const unsigned int domain_position, const Eigen::MatrixXd &points) const
    {
        if (domain_position > 2)
            throw std::runtime_error("Unknown domain position");

        std::array<Eigen::VectorXd, 3> result;

        const auto &x = points.row(0);
        const auto &y = points.row(1);
        const auto &z = points.row(2);

        result.at(0).setZero(points.cols());
        result.at(1).setZero(points.cols());
        result.at(2).setZero(points.cols());

        {
            switch (domain_position)
            {
            case 0:
                for (unsigned int p = 0; p < points.cols(); p++)
                {
                    result.at(0)[p] =
                        0.1 * (-0.5 - x(p)) *
                            (3.0 * x(p) * x(p) - 8.0 * x(p) * y(p) * y(p) + 16.0 * x(p) * x(p) * y(p) * atan2(y(p), x(p)) +
                             8.0 * y(p) * (x(p) * x(p) + y(p) * y(p)) * atan2(y(p), x(p))) +
                        0.1 * (-x(p) * x(p) * x(p) - 8.0 * x(p) * y(p) * (x(p) * x(p) + y(p) * y(p)) * atan2(y(p), x(p)));

                    result.at(1)[p] = 0.1 * (-0.5 - x(p)) *
                                      (8.0 * x(p) * x(p) * y(p) + 16.0 * x(p) * y(p) * y(p) * atan2(y(p), x(p)) +
                                       8.0 * x(p) * (x(p) * x(p) + y(p) * y(p)) * atan2(y(p), x(p)));
                }
                break;
            case 1:
                for (unsigned int p = 0; p < points.cols(); p++)
                {
                    result.at(0)[p] = 0.3 * (-0.5 - x(p)) * x(p) * x(p) - x(p) * x(p) * x(p) / 10.0 -
                                      2.4 * M_PI * (-0.5 - x(p)) * x(p) * x(p) * abs(z(p)) +
                                      0.8 * M_PI * x(p) * x(p) * x(p) * abs(z(p));

                    result.at(2)[p] = -1.0 * (z(p) == 0.0 ? 0.0 : (z(p) > 0.0 ? 1.0 : -1.0)) * 0.8 * M_PI *
                                      (-0.5 - x(p)) * x(p) * x(p) * x(p);
                }
                break;
            case 2:
                for (unsigned int p = 0; p < points.cols(); p++)
                {
                    result.at(1)[p] = (-1.0 + y(p)) * y(p) * (-1.0 + z(p)) * z(p) +
                                      (-1.0 + y(p)) * (1.0 + y(p)) * (-1.0 + z(p)) * z(p) +
                                      y(p) * (1.0 + y(p)) * (-1.0 + z(p)) * z(p);
                    result.at(2)[p] =
                        (-1.0 + y(p)) * y(p) * (1.0 + y(p)) * (-1.0 + z(p)) + (-1.0 + y(p)) * y(p) * (1.0 + y(p)) * z(p);
                }
                break;
            default:
                throw std::runtime_error("Domain index " + std::to_string(domain_position) + " not supported");
            }
        }

        return result;
    }
};

// ***************************************************************************
} // namespace test
} // namespace Elliptic_PCC_DFN
} // namespace examples
} // namespace Polydim

#endif
