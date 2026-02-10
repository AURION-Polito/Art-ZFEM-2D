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

#ifndef __program_configuration_H
#define __program_configuration_H

#include "Configurations.hpp"
#include "LocalSpace_PCC_2D.hpp"
#include "PDE_Mesh_Utilities.hpp"
#include "test_definition.hpp"

namespace Polydim
{
namespace examples
{
namespace Elliptic_PCC_2D
{
struct Program_Configuration final
{
    enum struct Solver_Types
    {
      Cholesky = 0,
      PCG = 1
    };

    Program_Configuration()
    {
        Gedim::Configurations::AddProperty("TestType",
                                           static_cast<unsigned int>(Polydim::examples::Elliptic_PCC_2D::test::Test_Types::Patch_Test),
                                           "Test Type 1 - Patch_Test; 2 - Elliptic_Polynomial_Problem; 3 - "
                                           "Elliptic_Problem; 4 - Computational_Comparison (Default: 1)");
        // Export parameters
        Gedim::Configurations::AddProperty("ExportFolder", "./Run", "Folder where to export data (Default: ./Export)");
        // Mesh parameters
        Gedim::Configurations::AddProperty(
            "MeshGenerator",
            static_cast<unsigned int>(Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Triangular),
            "Mesh 2D gereator type, 0 - Triangular; 1 - Minimal; 2 - "
            "Polygonal; 3 - OFF Importer; 4 - CsvImporter (; separator); 5 - Squared (Default: 0)");
        Gedim::Configurations::AddProperty("MeshImportFilePath", "./", "Mesh imported file path (Default: './')");
        Gedim::Configurations::AddProperty("MeshMaxArea", 0.1, "Mesh 2D maximum relative cell area (Default: 0.1)");

        Gedim::Configurations::AddProperty("GeometricTolerance1D", 1.0e-12, "Geometric Tolerance 1D (Default: 1.0e-12)");

        Gedim::Configurations::AddProperty("GeometricTolerance2D", 1.0e-14, "Geometric Tolerance 2D (Default: 1.0e-14)");

        Gedim::Configurations::AddProperty("SubTriangulate", false, "SubTraingulate Mesh (Default: false)");

        // Method parameters
        Gedim::Configurations::AddProperty("MethodType",
                                           static_cast<unsigned int>(PDETools::LocalSpace_PCC_2D::MethodTypes::ZFEM_PCC),
                                           "Method Type, 0 - FEM, 1 - EVem; 2 - EVem_Inertia; 3 - EVem_Ortho "
                                           "; 4 - ZFEM; (Default: 4)");

        Gedim::Configurations::AddProperty("MethodOrder", static_cast<unsigned int>(1), "Method order (Default: 1)");
        Gedim::Configurations::AddProperty("ComputeMethodPerformance", false, "Compute Method Performance (Default: false)");
        Gedim::Configurations::AddProperty("ComputeConditionNumber", false, "Compute Condition Number (Default: false)");
        Gedim::Configurations::AddProperty("PostProcess", true, "Post Process (Default: true)");
        Gedim::Configurations::AddProperty("ExportMatrix", false, "Export matrix (Default: false)");

        Gedim::Configurations::AddProperty("ComputationalTime",
                                           static_cast<unsigned int>(1),
                                           "Number of code execution to compute computational time (Default: 1)");

        Gedim::Configurations::AddProperty("SolverType",
                                           static_cast<unsigned int>(Polydim::examples::Elliptic_PCC_2D::Program_Configuration::Solver_Types::Cholesky),
                                           "Solver Type 0 - Cholesky; 1 - PCG (Default: 0)");
    }

    inline std::string ExportFolder() const
    {
        return Gedim::Configurations::GetPropertyValue<std::string>("ExportFolder");
    }

    inline Polydim::examples::Elliptic_PCC_2D::test::Test_Types TestType() const
    {
        return static_cast<Polydim::examples::Elliptic_PCC_2D::test::Test_Types>(
            Gedim::Configurations::GetPropertyValue<unsigned int>("TestType"));
    }
    inline Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D MeshGenerator() const
    {
        return static_cast<Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D>(
            Gedim::Configurations::GetPropertyValue<unsigned int>("MeshGenerator"));
    }
    inline std::string MeshImportFilePath() const
    {
        return Gedim::Configurations::GetPropertyValue<std::string>("MeshImportFilePath");
    }
    inline double MeshMaxArea() const
    {
        return Gedim::Configurations::GetPropertyValue<double>("MeshMaxArea");
    }
    inline double GeometricTolerance1D() const
    {
        return Gedim::Configurations::GetPropertyValue<double>("GeometricTolerance1D");
    }
    inline double GeometricTolerance2D() const
    {
        return Gedim::Configurations::GetPropertyValue<double>("GeometricTolerance2D");
    }

    inline PDETools::LocalSpace_PCC_2D::MethodTypes MethodType() const
    {
        return static_cast<PDETools::LocalSpace_PCC_2D::MethodTypes>(Gedim::Configurations::GetPropertyValue<unsigned int>("MethodType"));
    }

    inline bool ComputeMethodPerformance() const
    {
        return Gedim::Configurations::GetPropertyValue<bool>("ComputeMethodPerformance");
    }
    inline bool ComputeConditionNumber() const
    {
        return Gedim::Configurations::GetPropertyValue<bool>("ComputeConditionNumber");
    }
    inline bool PostProcess() const
    {
        return Gedim::Configurations::GetPropertyValue<bool>("PostProcess");
    }
    inline bool ExportMatrix() const
    {
        return Gedim::Configurations::GetPropertyValue<bool>("ExportMatrix");
    }
    inline unsigned int MethodOrder() const
    {
        return Gedim::Configurations::GetPropertyValue<unsigned int>("MethodOrder");
    }
    inline unsigned int ComputationalTime() const
    {
        return Gedim::Configurations::GetPropertyValue<unsigned int>("ComputationalTime");
    }

    inline bool SubTriangulate() const
    {
        return Gedim::Configurations::GetPropertyValue<bool>("SubTriangulate");
    }
    inline Polydim::examples::Elliptic_PCC_2D::Program_Configuration::Solver_Types SolverType() const
    {
        return static_cast<Polydim::examples::Elliptic_PCC_2D::Program_Configuration::Solver_Types>(
            Gedim::Configurations::GetPropertyValue<unsigned int>("SolverType"));
    }
};
} // namespace Elliptic_PCC_2D
} // namespace examples
} // namespace Polydim

#endif
