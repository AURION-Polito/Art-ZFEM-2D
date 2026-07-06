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
#include "MeshMatricesDAO_mesh_connectivity_data.hpp"
#include "VTKUtilities.hpp"
#include "program_utilities.hpp"
#include "test_definition.hpp"

int main(int argc, char **argv)
{
    Polydim::examples::Parabolic_PCC_2D::Program_Configuration config;

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

    const auto test = Polydim::examples::Parabolic_PCC_2D::program_utilities::create_test(config);

    const auto domain = test->domain();
    const auto boundary_info = test->boundary_info();

    // export domain
    {
        Gedim::VTKUtilities vtkUtilities;
        vtkUtilities.AddPolygon(domain.spatial_domain.vertices);
        vtkUtilities.Export(exportVtuFolder + "/Domain.vtu");
    }

    Gedim::Profiler::StopTime("SetProblem");
    Gedim::Output::PrintStatusProgram("SetProblem");

    /// Create domain mesh

    Gedim::Output::PrintGenericMessage("CreateMesh...", true);
    Gedim::Profiler::StartTime("CreateMesh");

    Gedim::MeshMatrices meshData;
    Gedim::MeshMatricesDAO mesh(meshData);
    Gedim::MeshUtilities::MeshGeometricData2D mesh_geometric_data;

    if (!config.SubTriangulate())
    {
        Polydim::examples::Parabolic_PCC_2D::program_utilities::create_domain_mesh(config, domain.spatial_domain, mesh);
    }
    else
    {
        Gedim::MeshMatrices polygonal_mesh_data;
        Gedim::MeshMatricesDAO polygonal_mesh(polygonal_mesh_data);

        Polydim::examples::Parabolic_PCC_2D::program_utilities::create_domain_mesh(config, domain.spatial_domain, polygonal_mesh);

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

        Polydim::examples::Parabolic_PCC_2D::program_utilities::make_mesh_triangular_by_internal_point(internal_points, polygonal_mesh);

        const auto filter_data = meshUtilities.FilterActiveMesh(polygonal_mesh);
        meshUtilities.ExtractMesh2D(filter_data.Cell0Ds, filter_data.Cell1Ds, filter_data.Cell2Ds, polygonal_mesh, mesh);
    }

    std::vector<unsigned int> region_id(mesh.Cell2DTotalNumber(), 0);

    switch (config.MeshGenerator())
    {

    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Triangular:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Minimal:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Polygonal:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::OFFImporter:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Squared:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::RandomDistorted:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::TriangularSimpleImporter:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::StructuredTriangular:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::QuadFromTriangular:
        break;
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::CsvImporter: {
        std::ifstream layer_file(config.MeshImportFilePath() + "/RegionId.csv");
        if (!layer_file.fail())
        {
            unsigned int num_lines;
            layer_file >> num_lines;

            if (num_lines != mesh.Cell2DTotalNumber())
                throw std::runtime_error("not valid region id file");

            unsigned int cell_id = 0;
            char tmp = ';';
            unsigned int marker_layer = 0;
            unsigned int c = 0;
            while (layer_file >> cell_id >> tmp >> marker_layer)
                region_id[c++] = marker_layer;
        }
    }
    break;
    default:
        throw std::runtime_error("not valid mesh generator");
    }

    // Export the domain mesh
    {
        meshUtilities.ExportMeshToVTU(mesh, exportVtuFolder, "Domain_Mesh");
    }

    Gedim::Profiler::StopTime("CreateMesh");
    Gedim::Output::PrintStatusProgram("CreateMesh");

    Gedim::Output::PrintGenericMessage("ComputeGeometricProperties...", true);
    Gedim::Profiler::StartTime("ComputeGeometricProperties");

    mesh_geometric_data =
        Polydim::examples::Parabolic_PCC_2D::program_utilities::create_domain_mesh_geometric_properties(config, mesh);

    Gedim::Profiler::StopTime("ComputeGeometricProperties");
    Gedim::Output::PrintStatusProgram("ComputeGeometricProperties");

    /// Initialize Discrete Space
    Gedim::Output::PrintGenericMessage("CreateDiscreteSpace of order " + std::to_string(config.MethodOrder()) + " and DOFs...", true);
    Gedim::Profiler::StartTime("CreateDiscreteSpace");

    const auto reference_element_data =
        Polydim::PDETools::LocalSpace_PCC_2D::CreateReferenceElement(config.MethodType(), config.MethodOrder());

    Polydim::PDETools::Mesh::MeshMatricesDAO_mesh_connectivity_data mesh_connectivity_data = {mesh};

    Polydim::PDETools::DOFs::DOFsManager dofManager;

    const auto mesh_dofs_info = Polydim::PDETools::LocalSpace_PCC_2D::SetMeshDOFsInfo(reference_element_data, mesh, boundary_info);
    const auto dofs_data = dofManager.CreateDOFs_2D(mesh_dofs_info, mesh_connectivity_data);

    Gedim::Output::PrintGenericMessage("Discrete Space with " + std::to_string(dofs_data.NumberDOFs) + " DOFs and " +
                                           std::to_string(dofs_data.NumberStrongs) + " STRONGs",
                                       true);

    Gedim::Profiler::StopTime("CreateDiscreteSpace");
    Gedim::Output::PrintStatusProgram("CreateDiscreteSpace");

    Polydim::examples::Parabolic_PCC_2D::Assembler assembler;
    Polydim::examples::Parabolic_PCC_2D::Assembler::Parabolic_PCC_2D_Problem_Data assembler_data;

    std::vector<double> time_steps =
        Polydim::examples::Parabolic_PCC_2D::program_utilities::create_time_steps(config, domain.time_domain);
    const double theta_parameter = config.ThetaParameter();

    assembler.ComputeInitialCondition(config,
                                      mesh,
                                      mesh_geometric_data,
                                      region_id,
                                      dofs_data,
                                      reference_element_data,
                                      *test,
                                      assembler_data.initial_solution,
                                      assembler_data.initial_solution_dirichlet);

    assembler.AssembleMassMatrix(config,
                                 mesh,
                                 mesh_geometric_data,
                                 mesh_dofs_info,
                                 dofs_data,
                                 reference_element_data,
                                 assembler_data.globalMatrixM,
                                 assembler_data.dirichletMatrixM);

    assembler.Assemble(config,
                       time_steps[0],
                       mesh,
                       mesh_geometric_data,
                       region_id,
                       mesh_dofs_info,
                       dofs_data,
                       reference_element_data,
                       *test,
                       assembler_data.initial_globalMatrixA,
                       assembler_data.initial_dirichletMatrixA,
                       assembler_data.initial_rightHandSide,
                       assembler_data.initial_solution_dirichlet);

    for (unsigned int t = 1; t < time_steps.size(); t++)
    {
        const double time_value = time_steps[t];
        const double delta_time = time_steps[t] - time_steps[t - 1];

        assembler.Assemble(config,
                           time_value,
                           mesh,
                           mesh_geometric_data,
                           region_id,
                           mesh_dofs_info,
                           dofs_data,
                           reference_element_data,
                           *test,
                           assembler_data.globalMatrixA,
                           assembler_data.dirichletMatrixA,
                           assembler_data.rightHandSide,
                           assembler_data.solution_dirichlet);

        Gedim::Eigen_SparseArray<> global_lhs;
        Gedim::Eigen_Array<> global_rhs;

        global_lhs.SetSize(dofs_data.NumberDOFs, dofs_data.NumberDOFs, Gedim::ISparseArray::SparseArrayTypes::Symmetric);
        global_rhs.SetSize(dofs_data.NumberDOFs);

        global_lhs = theta_parameter * delta_time * assembler_data.globalMatrixA;
        global_lhs += assembler_data.globalMatrixM;

        assembler_data.initial_rightHandSide.SubtractionMultiplication(assembler_data.initial_globalMatrixA,
                                                                       assembler_data.initial_solution);
        assembler_data.initial_rightHandSide *= (1 - theta_parameter);
        global_rhs = assembler_data.rightHandSide;
        global_rhs *= theta_parameter;
        global_rhs += assembler_data.initial_rightHandSide;
        global_rhs *= delta_time;

        global_rhs.SubtractionMultiplication(assembler_data.dirichletMatrixM, assembler_data.solution_dirichlet);
        global_rhs.SumMultiplication(assembler_data.globalMatrixM, assembler_data.initial_solution);
        global_rhs.SumMultiplication(assembler_data.dirichletMatrixM, assembler_data.initial_solution_dirichlet);

        if (dofs_data.NumberDOFs > 0)
        {
            Gedim::Eigen_CholeskySolver solver;
            solver.Initialize(global_lhs);
            solver.Solve(global_rhs, assembler_data.solution);
        }

        Gedim::Output::PrintGenericMessage("ComputeErrors...", true);
        Gedim::Profiler::StartTime("ComputeErrors");

        const auto post_process_data =
            assembler.PostProcessSolution(config, time_value, mesh, mesh_geometric_data, dofs_data, reference_element_data, assembler_data, *test);

        Gedim::Profiler::StopTime("ComputeErrors");
        Gedim::Output::PrintStatusProgram("ComputeErrors");

        Gedim::Output::PrintGenericMessage("ExportSolution...", true);
        Gedim::Profiler::StartTime("ExportSolution");

        Polydim::examples::Parabolic_PCC_2D::program_utilities::export_solution(config,
                                                                                time_value,
                                                                                delta_time,
                                                                                t,
                                                                                theta_parameter,
                                                                                mesh,
                                                                                dofs_data,
                                                                                assembler_data,
                                                                                post_process_data,
                                                                                exportSolutionFolder,
                                                                                exportVtuFolder);

        Gedim::Profiler::StopTime("ExportSolution");
        Gedim::Output::PrintStatusProgram("ExportSolution");

        assembler_data.solution.Copy(assembler_data.initial_solution);
        assembler_data.solution_dirichlet.Copy(assembler_data.initial_solution_dirichlet);

        assembler_data.initial_globalMatrixA = assembler_data.globalMatrixA;
        assembler_data.initial_dirichletMatrixA = assembler_data.dirichletMatrixA;
        assembler_data.initial_rightHandSide = assembler_data.rightHandSide;
    }

    Gedim::Output::PrintGenericMessage("ComputeMethodPerformance...", true);
    Gedim::Profiler::StartTime("ComputeMethodPerformance");

    if (config.ComputeMethodPerformance())
    {
        const auto performance = assembler.ComputePerformance(config, mesh, mesh_geometric_data, reference_element_data);

        Polydim::examples::Parabolic_PCC_2D::program_utilities::export_performance(config, performance, exportSolutionFolder);
    }

    Gedim::Profiler::StopTime("ComputeMethodPerformance");
    Gedim::Output::PrintStatusProgram("ComputeMethodPerformance");

    return 0;
}
