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

#include "Eigen_CholeskySolver.hpp"
#include "Eigen_LUSolver.hpp"
#include "Eigen_PCGSolver.hpp"
#include "MeshMatricesDAO_mesh_connectivity_data.hpp"
#include "VTKUtilities.hpp"
#include "program_utilities.hpp"
#include "test_definition.hpp"

unsigned int Polydim::examples::Elliptic_PCC_2D::test::Patch_Test::order;

int main(int argc, char **argv)
{
    Polydim::examples::Elliptic_PCC_2D::Program_Configuration config;

    if (!Gedim::Output::FileExists("./Parameters.ini"))
        Gedim::Configurations::ExportToIni("./Parameters.ini", false);
    else
        Gedim::Configurations::InitializeFromIni("./Parameters.ini");

    Gedim::Configurations::Initialize(argc, argv);

    /// Create folders
    const std::string exportFolder = config.ExportFolder();
    Gedim::Output::CreateFolder(exportFolder);

    const std::string exportCsvFolder = exportFolder + "/Mesh";
    Gedim::Output::CreateFolder(exportCsvFolder);
    const std::string exportVtuFolder = exportFolder + "/Paraview";
    Gedim::Output::CreateFolder(exportVtuFolder);
    const std::string exportSolutionFolder = exportFolder + "/Solution";
    Gedim::Output::CreateFolder(exportSolutionFolder);

    const std::string logFolder = exportFolder + "/Log";

    Gedim::MeshUtilities meshUtilities;

    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    /// Set Profiler
    Gedim::Profiler::ActivateProfiler = true;

    /// Set Log folder
    Gedim::Output::CreateFolder(logFolder);
    Gedim::LogFile::LogFolder = logFolder;

    /// Export Configuration of the following Run
    Gedim::Configurations::ExportToIni(exportFolder + "/Parameters.ini", false);

    /// Set problem
    Gedim::Output::PrintGenericMessage("SetProblem...", true);
    Gedim::Profiler::StartTime("SetProblem");

    Polydim::examples::Elliptic_PCC_2D::test::Patch_Test::order = config.MethodOrder();

    const auto test = Polydim::examples::Elliptic_PCC_2D::program_utilities::create_test(config);

    const auto domain = test->domain();
    const auto boundary_info = test->boundary_info();

    // export domain
    {
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPolygon(domain.vertices);
        vtkUtilities.Export(exportVtuFolder + "/Domain.vtu");
    }

    Gedim::Profiler::StopTime("SetProblem");
    Gedim::Output::PrintStatusProgram("SetProblem");

    /// Create domain mesh
    Gedim::MeshMatrices meshData;
    Gedim::MeshMatricesDAO mesh(meshData);
    Gedim::MeshUtilities::MeshGeometricData2D meshGeometricData;

    if (!config.SubTriangulate())
    {

        Gedim::Output::PrintGenericMessage("CreateMesh...", true);
        Gedim::Profiler::StartTime("CreateMesh");

        Polydim::examples::Elliptic_PCC_2D::program_utilities::create_domain_mesh(config, domain, mesh);

        Gedim::Profiler::StopTime("CreateMesh");
        Gedim::Output::PrintStatusProgram("CreateMesh");
    }
    else
    {
        Gedim::MeshMatrices polygonal_mesh_data;
        Gedim::MeshMatricesDAO polygonal_mesh(polygonal_mesh_data);

        Gedim::Output::PrintGenericMessage("CreateMesh...", true);
        Gedim::Profiler::StartTime("CreateMesh");

        Polydim::examples::Elliptic_PCC_2D::program_utilities::create_domain_mesh(config, domain, polygonal_mesh);

        Gedim::Profiler::StopTime("CreateMesh");
        Gedim::Output::PrintStatusProgram("CreateMesh");

        // Export the domain mesh
        {
            meshUtilities.ExportMeshToVTU(polygonal_mesh, exportVtuFolder, "Domain_PolygonalMesh");
        }

        Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
        geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
        geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
        Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

        std::vector<Eigen::Vector3d> internal_points(polygonal_mesh.Cell2DTotalNumber());
        for (unsigned int c = 0; c < polygonal_mesh.Cell2DTotalNumber(); c++)
        {
            const Eigen::MatrixXd &cell2DVertices = polygonal_mesh.Cell2DVerticesCoordinates(c);

            const double Cell2DsDiameters = geometryUtilities.PolygonDiameter(cell2DVertices);
            const auto Cell2DsEdgeNormals = geometryUtilities.PolygonEdgeNormals(cell2DVertices);
            double in_radius;
            geometryUtilities.PolygonChebyshevCenter(cell2DVertices, Cell2DsEdgeNormals, internal_points[c], in_radius, Cell2DsDiameters, true);
        }

        Polydim::examples::Elliptic_PCC_2D::program_utilities::make_mesh_triangular_by_internal_point(internal_points, polygonal_mesh);

        const auto filter_data = meshUtilities.FilterActiveMesh(polygonal_mesh);
        meshUtilities.ExtractMesh2D(filter_data.Cell0Ds, filter_data.Cell1Ds, filter_data.Cell2Ds, polygonal_mesh, mesh);
    }

    if (config.MeshGenerator() == Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::OFFImporter)
    {
        const unsigned int num_vertices_domain = domain.vertices.cols();
        for (unsigned int e = 0; e < num_vertices_domain; e++)
        {
            const Eigen::Vector3d start = domain.vertices.col(e);
            const Eigen::Vector3d tangent = domain.vertices.col((e + 1) % num_vertices_domain) - start;
            const double tangent_squared_length = tangent.squaredNorm();
            meshUtilities.SetMeshMarkersOnSegment(geometryUtilities, start, tangent, tangent_squared_length, e + 1 + num_vertices_domain, mesh);
        }

        for (unsigned int v = 0; v < num_vertices_domain; v++)
        {
            const Eigen::Vector3d vert = domain.vertices.col(v);

            for (unsigned int c = 0; c < mesh.Cell0DTotalNumber(); c++)
            {
                if (geometryUtilities.PointsAreCoincident(vert, mesh.Cell0DCoordinates(c)))
                {
                    mesh.Cell0DSetMarker(c, v + 1);
                    break;
                }
            }
        }
    }

    // Export the domain mesh
    {
        meshUtilities.ExportMeshToVTU(mesh, exportVtuFolder, "Domain_Mesh");
    }

    Gedim::Output::PrintGenericMessage("ComputeGeometricProperties...", true);
    Gedim::Profiler::StartTime("ComputeGeometricProperties");

    meshGeometricData =
        Polydim::examples::Elliptic_PCC_2D::program_utilities::create_domain_mesh_geometric_properties(config, mesh);

    Gedim::Profiler::StopTime("ComputeGeometricProperties");
    Gedim::Output::PrintStatusProgram("ComputeGeometricProperties");

    /// Initialize Discrete Space
    Gedim::Output::PrintGenericMessage("CreateDiscreteSpace of order " + std::to_string(config.MethodOrder()) + " and DOFs...", true);
    Gedim::Profiler::StartTime("CreateDiscreteSpace");

    const auto reference_element_data =
        Polydim::PDETools::LocalSpace_PCC_2D::CreateReferenceElement(config.MethodType(), config.MethodOrder());

    Polydim::PDETools::Mesh::MeshMatricesDAO_mesh_connectivity_data mesh_connectivity_data = {mesh};

    Polydim::PDETools::DOFs::DOFsManager dofManager;

    const auto meshDOFsInfo = Polydim::PDETools::LocalSpace_PCC_2D::SetMeshDOFsInfo(reference_element_data, mesh, boundary_info);
    const auto dofs_data = dofManager.CreateDOFs_2D(meshDOFsInfo, mesh_connectivity_data);

    Gedim::Output::PrintGenericMessage("Discrete Space with " + std::to_string(dofs_data.NumberDOFs) + " DOFs and " +
                                           std::to_string(dofs_data.NumberStrongs) + " STRONGs",
                                       true);

    Gedim::Profiler::StopTime("CreateDiscreteSpace");
    Gedim::Output::PrintStatusProgram("CreateDiscreteSpace");

    Polydim::examples::Elliptic_PCC_2D::Assembler assembler;
    Polydim::examples::Elliptic_PCC_2D::Assembler::Elliptic_PCC_2D_Problem_Data assembler_data;

    Gedim::ILinearSolver::SolutionInfo solver_data;
    double time_assembler = 0.0;
    double time_solver = 0.0;
    for (unsigned int i = 0; i < config.ComputationalTime(); i++)
    {
        const auto start_time_assembler = Gedim::Profiler::GetTime();

        assembler_data =
            assembler.Assemble(config, mesh, meshGeometricData, meshDOFsInfo, dofs_data, reference_element_data, *test);

        const auto end_time_assembler = Gedim::Profiler::GetTime();

        time_assembler += Gedim::Profiler::ComputeTime(start_time_assembler, end_time_assembler);

        switch (config.SolverType())
        {
          case Polydim::examples::Elliptic_PCC_2D::Program_Configuration::Solver_Types::PCG:
          {
            const auto start_time_solver = Gedim::Profiler::GetTime();
            if (dofs_data.NumberDOFs > 0)
            {
                 Gedim::Eigen_PCGSolver<> solver;
                 solver.Initialize(assembler_data.globalMatrixA, { dofs_data.NumberDOFs, 1.0e-15 });
                 solver_data = solver.Solve(assembler_data.rightHandSide, assembler_data.solution);
            }
            const auto end_time_solver = Gedim::Profiler::GetTime();
            time_solver += Gedim::Profiler::ComputeTime(start_time_solver, end_time_solver);
          }
            break;
          case Polydim::examples::Elliptic_PCC_2D::Program_Configuration::Solver_Types::Cholesky:
          {
            const auto start_time_solver = Gedim::Profiler::GetTime();
            if (dofs_data.NumberDOFs > 0)
            {
                Gedim::Eigen_CholeskySolver solver;
                solver.Initialize(assembler_data.globalMatrixA);
                solver_data = solver.Solve(assembler_data.rightHandSide, assembler_data.solution);
            }
            const auto end_time_solver = Gedim::Profiler::GetTime();
            time_solver += Gedim::Profiler::ComputeTime(start_time_solver, end_time_solver);
          }
            break;
          default:
            throw std::runtime_error("Unknown solver type");
        }
    }

    time_assembler /= config.ComputationalTime();
    time_solver /= config.ComputationalTime();

    const unsigned int Method_ID = static_cast<unsigned int>(config.MethodType());
    const unsigned int TEST_ID = static_cast<unsigned int>(config.TestType());

    Gedim::Output::PrintGenericMessage("ComputeErrors...", true);
    Gedim::Profiler::StartTime("ComputeErrors");

    const auto post_process_data =
        assembler.PostProcessSolution(config, mesh, meshGeometricData, dofs_data, reference_element_data, assembler_data, *test);

    Gedim::Profiler::StopTime("ComputeErrors");
    Gedim::Output::PrintStatusProgram("ComputeErrors");


    if (config.ExportMatrix())
    {
      Gedim::Output::PrintGenericMessage("ExportMatrix...", true);
      Gedim::Profiler::StartTime("ExportMatrix");

      std::string s;
      std::ostringstream str(s);
      str.precision(6);
      str << std::scientific << post_process_data.mesh_size;

      assembler_data.globalMatrixA.ToBinaryFile(exportSolutionFolder + "/Matrix_" + std::to_string(TEST_ID) + "_" +
                                                std::to_string(Method_ID) + "_" + std::to_string(config.MethodOrder()) +
                                                "_" + str.str() + ".txt");

      Gedim::Profiler::StopTime("ExportMatrix");
      Gedim::Output::PrintStatusProgram("ExportMatrix");
    }

    Gedim::Output::PrintGenericMessage("ExportSolution...", true);
    Gedim::Profiler::StartTime("ExportSolution");

    Polydim::examples::Elliptic_PCC_2D::program_utilities::export_solution(config,
                                                                           mesh,
                                                                           dofs_data,
                                                                           assembler_data,
                                                                           post_process_data,
                                                                           time_assembler,
                                                                           time_solver,
                                                                           solver_data,
                                                                           exportSolutionFolder,
                                                                           exportVtuFolder);
    if (config.PostProcess())
        Polydim::examples::Elliptic_PCC_2D::program_utilities::export_dofs(config,
                                                                           mesh,
                                                                           meshGeometricData,
                                                                           meshDOFsInfo,
                                                                           dofs_data,
                                                                           assembler_data,
                                                                           post_process_data,
                                                                           exportVtuFolder);

    Gedim::Profiler::StopTime("ExportSolution");
    Gedim::Output::PrintStatusProgram("ExportSolution");

    Gedim::Output::PrintGenericMessage("ComputeMethodPerformance...", true);
    Gedim::Profiler::StartTime("ComputeMethodPerformance");

    if (config.ComputeMethodPerformance())
    {
        const auto performance = assembler.ComputePerformance(config, mesh, meshGeometricData, reference_element_data);

        Polydim::examples::Elliptic_PCC_2D::program_utilities::export_performance(config, performance, exportSolutionFolder);
    }

    Gedim::Profiler::StopTime("ComputeMethodPerformance");
    Gedim::Output::PrintStatusProgram("ComputeMethodPerformance");

    return 0;
}
