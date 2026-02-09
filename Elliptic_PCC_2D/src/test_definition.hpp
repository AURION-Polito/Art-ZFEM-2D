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
namespace Elliptic_PCC_2D
{
namespace test
{
enum struct Test_Types
{
    Patch_Test = 1,
    Elliptic_Polynomial_Problem = 2,
    Elliptic_Problem = 3,
    Computational_Comparison = 4
};

struct I_Test
{
    virtual Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain() const = 0;
    virtual std::map<unsigned int, Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo::BoundaryInfo> boundary_info() const = 0;
    virtual std::array<Eigen::VectorXd, 9> diffusion_term(const Eigen::MatrixXd &points) const = 0;
    virtual Eigen::VectorXd reaction_term(const Eigen::MatrixXd &points) const = 0;
    virtual Eigen::VectorXd source_term(const Eigen::MatrixXd &points) const = 0;
    virtual Eigen::VectorXd strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const = 0;
    virtual Eigen::VectorXd weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const = 0;
    virtual Eigen::VectorXd exact_solution(const Eigen::MatrixXd &points) const = 0;
    virtual std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd &points) const = 0;
};
// ***************************************************************************
struct Elliptic_Problem final : public I_Test
{
    Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain() const
    {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain;

        domain.area = 1.0;

        domain.vertices = Eigen::MatrixXd::Zero(3, 4);
        domain.vertices.row(0) << 0.0, 1.0, 1.0, 0.0;
        domain.vertices.row(1) << 0.0, 0.0, 1.0, 1.0;

        domain.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types::Parallelogram;

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

    std::array<Eigen::VectorXd, 9> diffusion_term(const Eigen::MatrixXd &points) const
    {
        return {1.0 + points.row(1).array() * points.row(1).array(),
                -points.row(0).array() * points.row(1).array(),
                Eigen::VectorXd::Zero(points.cols()),
                -points.row(0).array() * points.row(1).array(),
                1.0 + points.row(0).array() * points.row(0).array(),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Zero(points.cols()),
                Eigen::VectorXd::Constant(points.cols(), 0.0)};
    };

    Eigen::VectorXd reaction_term(const Eigen::MatrixXd &points) const
    {
        return points.row(0).array() * points.row(1).array();
    };

    Eigen::VectorXd source_term(const Eigen::MatrixXd &points) const
    {
        Eigen::ArrayXd x = points.row(0);
        Eigen::ArrayXd y = points.row(1);
        return 4.0 * std::numbers::pi * std::numbers::pi *
                   ((2.0 + x * x + y * y) * sin(2.0 * std::numbers::pi * x) * sin(2.0 * std::numbers::pi * y) +
                    2.0 * x * y * cos(2.0 * std::numbers::pi * x) * cos(2.0 * std::numbers::pi * y)) +
               2.0 * y * std::numbers::pi * sin(2.0 * std::numbers::pi * x) * cos(2.0 * std::numbers::pi * y) +
               2.0 * x * std::numbers::pi * cos(2.0 * std::numbers::pi * x) * sin(2.0 * std::numbers::pi * y) +
               reaction_term(points).array() * exact_solution(points).array();
    };

    Eigen::VectorXd strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        return exact_solution(points);
    };

    Eigen::VectorXd weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        switch (marker)
        {
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    Eigen::VectorXd exact_solution(const Eigen::MatrixXd &points) const
    {
        return sin(2.0 * std::numbers::pi * points.row(0).array()) * sin(2.0 * std::numbers::pi * points.row(1).array());
    };

    std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd &points) const
    {
        return {2.0 * std::numbers::pi * cos(2.0 * std::numbers::pi * points.row(0).array()) *
                    sin(2.0 * std::numbers::pi * points.row(1).array()),
                2.0 * std::numbers::pi * sin(2.0 * std::numbers::pi * points.row(0).array()) *
                    cos(2.0 * std::numbers::pi * points.row(1).array()),
                Eigen::VectorXd::Zero(points.cols())};
    }
};
// ***************************************************************************
struct Patch_Test final : public I_Test
{
    static unsigned int order;

    Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain() const
    {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain;

        domain.area = 1.0;

        domain.vertices = Eigen::MatrixXd::Zero(3, 4);
        domain.vertices.row(0) << 0.0, 1.0, 1.0, 0.0;
        domain.vertices.row(1) << 0.0, 0.0, 1.0, 1.0;

        domain.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types::Parallelogram;

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

    std::array<Eigen::VectorXd, 9> diffusion_term(const Eigen::MatrixXd &points) const
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
    };

    Eigen::VectorXd reaction_term(const Eigen::MatrixXd &points) const
    {
        const double k = 0.0;
        return Eigen::VectorXd::Constant(points.cols(), k);
    };

    Eigen::VectorXd source_term(const Eigen::MatrixXd &points) const
    {
        Eigen::VectorXd source_term = Eigen::VectorXd::Constant(points.cols(), 2.0 * order * (order - 1));
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        const int max_order = order - 2;
        for (int i = 0; i < max_order; ++i)
            source_term.array() *= polynomial;

        return -source_term;
    };

    Eigen::VectorXd strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        Eigen::VectorXd result = Eigen::VectorXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order; ++i)
            result.array() *= polynomial;

        return result;
    };

    Eigen::VectorXd weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {

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

    Eigen::VectorXd exact_solution(const Eigen::MatrixXd &points) const
    {
        const Eigen::ArrayXd polynomial = points.row(0).array() + points.row(1).array() + 0.5;

        Eigen::VectorXd result = Eigen::VectorXd::Constant(points.cols(), 1.0);
        for (int i = 0; i < order; ++i)
            result.array() *= polynomial;

        return result;
    };

    std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd &points) const
    {
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
    Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain() const
    {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain;

        domain.area = 1.0;

        domain.vertices = Eigen::MatrixXd::Zero(3, 4);
        domain.vertices.row(0) << 0.0, 1.0, 1.0, 0.0;
        domain.vertices.row(1) << 0.0, 0.0, 1.0, 1.0;

        domain.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types::Parallelogram;

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
    std::array<Eigen::VectorXd, 9> diffusion_term(const Eigen::MatrixXd &points) const
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
    };

    Eigen::VectorXd reaction_term(const Eigen::MatrixXd &points) const
    {
        const double k = 0.0;
        return Eigen::VectorXd::Constant(points.cols(), k);
    };

    Eigen::VectorXd source_term(const Eigen::MatrixXd &points) const
    {
        return 32.0 * (points.row(1).array() * (1.0 - points.row(1).array()) +
                       points.row(0).array() * (1.0 - points.row(0).array()));
    };

    Eigen::VectorXd strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        return 16.0 * (points.row(1).array() * (1.0 - points.row(1).array()) * points.row(0).array() *
                       (1.0 - points.row(0).array())) +
               1.1;
    };

    Eigen::VectorXd weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
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

    Eigen::VectorXd exact_solution(const Eigen::MatrixXd &points) const
    {
        return 16.0 * (points.row(1).array() * (1.0 - points.row(1).array()) * points.row(0).array() *
                       (1.0 - points.row(0).array())) +
               1.1;
    };

    std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd &points) const
    {
        return {16.0 * (1.0 - 2.0 * points.row(0).array()) * points.row(1).array() * (1.0 - points.row(1).array()),
                16.0 * (1.0 - 2.0 * points.row(1).array()) * points.row(0).array() * (1.0 - points.row(0).array()),
                Eigen::VectorXd::Zero(points.cols())};
    }
};
// ***************************************************************************
struct Computational_Comparison final : public I_Test
{
  private:
    double eps;
    double c;
    double s;
    double power;

  public:
    Computational_Comparison()
    {
        eps = 1.0e-13;
        c = 2;
        s = (0.5 * M_PI) * exp(c);
        power = 20.0;
    }

    Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain() const
    {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D domain;

        domain.area = 1.0;

        domain.vertices = Eigen::MatrixXd::Zero(3, 4);
        domain.vertices.row(0) << 0.0, 1.0, 1.0, 0.0;
        domain.vertices.row(1) << 0.0, 0.0, 1.0, 1.0;

        domain.shape_type = Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D::Domain_Shape_Types::Parallelogram;

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

    std::array<Eigen::VectorXd, 9> diffusion_term(const Eigen::MatrixXd &points) const
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
    };

    Eigen::VectorXd reaction_term(const Eigen::MatrixXd &points) const
    {
        return points.row(0).array() * points.row(1).array();
    };

    Eigen::VectorXd source_term(const Eigen::MatrixXd &points) const
    {
        const Eigen::Index n = points.cols();
        Eigen::VectorXd lap(n);
        lap.setZero();

        for (Eigen::Index i = 0; i < n; ++i)
        {
            const double x = points(0, i);
            const double y = points(1, i);

            const double xb = x + eps;
            const double yb = y + eps;

            if (std::abs(xb) < eps || std::abs(yb) < eps)
                throw std::runtime_error("error on f evaluation, NaN");

            // Helpers for x
            const double Ex = std::exp(-c / xb);
            const double Ax = s * Ex;
            const double Fx = 1.0 - std::cos(Ax);

            // F''(x)
            const double xb2 = xb * xb;
            const double xb4 = xb2 * xb2;
            const double common_x = (c * s * Ex) / xb4; // (cs e^{-c/(x+b)})/(x+b)^4
            const double Fxx = common_x * ((c * s * Ex) * std::cos(Ax) + (c - 2.0 * xb) * std::sin(Ax));

            // Helpers for y
            const double Ey = std::exp(-c / yb);
            const double Ay = s * Ey;
            const double Gy = 1.0 - std::cos(Ay);

            // G''(y)
            const double yb2 = yb * yb;
            const double yb4 = yb2 * yb2;
            const double common_y = (c * s * Ey) / yb4;
            const double Gyy = common_y * ((c * s * Ey) * std::cos(Ay) + (c - 2.0 * yb) * std::sin(Ay));

            // Laplacian: Δu = p [F''(x) G(y) + F(x) G''(y)]
            lap(i) = power * (Fxx * Gy + Fx * Gyy);
        }

        return -lap;
    };

    Eigen::VectorXd strong_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        if (marker != 1)
            throw std::runtime_error("Unknown marker");

        return exact_solution(points);
    };

    Eigen::VectorXd weak_boundary_condition(const unsigned int marker, const Eigen::MatrixXd &points) const
    {
        switch (marker)
        {
        default:
            throw std::runtime_error("Unknown marker");
        }
    }

    Eigen::VectorXd exact_solution(const Eigen::MatrixXd &points) const
    {
        return power * (1.0 - cos(exp(-c / (points.row(0).array() + eps)) * s)) *
               (1.0 - cos(exp(-c / (points.row(1).array() + eps)) * s));
    };

    std::array<Eigen::VectorXd, 3> exact_derivative_solution(const Eigen::MatrixXd &points) const
    {
        return {power * sin(exp(-c / (points.row(0).array() + eps)) * s) *
                    (1.0 - cos(exp(-c / (points.row(1).array() + eps)) * s)) * s *
                    exp(-c / (points.row(0).array() + eps)) * c / (points.row(0).array() + eps).square(),

                power * sin(exp(-c / (points.row(1).array() + eps)) * s) *
                    (1.0 - cos(exp(-c / (points.row(0).array() + eps)) * s)) * s *
                    exp(-c / (points.row(1).array() + eps)) * c / (points.row(1).array() + eps).square(),

                Eigen::VectorXd::Zero(points.cols())};
    }
};
// ***************************************************************************

} // namespace test
} // namespace Elliptic_PCC_2D
} // namespace examples
} // namespace Polydim

#endif
