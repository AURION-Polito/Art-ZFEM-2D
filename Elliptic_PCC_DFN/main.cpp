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
#include "program_utilities.hpp"
#include "test_definition.hpp"

unsigned int Polydim::examples::Elliptic_PCC_DFN::test::Patch_Test::order;

int main(int argc, char **argv)
{
    Polydim::examples::Elliptic_PCC_DFN::Program_configuration config;

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

    Polydim::examples::Elliptic_PCC_DFN::test::Patch_Test::order = config.MethodOrder();

    const auto test = Polydim::examples::Elliptic_PCC_DFN::program_utilities::create_test(config);

    const auto domains = test->domains();
    const auto boundary_info = test->boundary_info();

    Polydim::examples::Elliptic_PCC_DFN::program_utilities::export_domains(config, domains, exportVtuFolder);

    Gedim::Profiler::StopTime("SetProblem");
    Gedim::Output::PrintStatusProgram("SetProblem");

    /// Create domain mesh
    Gedim::Output::PrintGenericMessage("CreateMesh...", true);
    Gedim::Profiler::StartTime("CreateMesh");

    Gedim::MeshMatrices meshData;
    Gedim::MeshMatricesDAO mesh(meshData);

    Polydim::examples::Elliptic_PCC_DFN::program_utilities::create_domain_mesh(config, domains, mesh);

    Gedim::Profiler::StopTime("CreateMesh");
    Gedim::Output::PrintStatusProgram("CreateMesh");

    Gedim::Output::PrintGenericMessage("ComputeGeometricProperties...", true);
    Gedim::Profiler::StartTime("ComputeGeometricProperties");

    const auto meshGeometricData =
        Polydim::examples::Elliptic_PCC_DFN::program_utilities::create_domain_mesh_geometric_properties(config, domains, mesh);

    Polydim::examples::Elliptic_PCC_DFN::program_utilities::export_domain_mesh(config, domains, mesh, meshGeometricData, exportVtuFolder);

    Gedim::Profiler::StopTime("ComputeGeometricProperties");
    Gedim::Output::PrintStatusProgram("ComputeGeometricProperties");

    /// Initialize Discrete Space
    Gedim::Output::PrintGenericMessage("CreateDiscreteSpace of order " + std::to_string(config.MethodOrder()) + " and DOFs...", true);
    Gedim::Profiler::StartTime("CreateDiscreteSpace");

    const auto reference_element_data =
        Polydim::PDETools::LocalSpace_PCC_2D::CreateReferenceElement(config.MethodType(), config.MethodOrder());

    Polydim::PDETools::Mesh::MeshMatricesDAO_mesh_connectivity_data mesh_connectivity_data(mesh);

    Polydim::PDETools::DOFs::DOFsManager dofManager;

    const auto meshDOFsInfo = Polydim::PDETools::LocalSpace_PCC_2D::SetMeshDOFsInfo(reference_element_data, mesh, boundary_info);
    const auto dofs_data = dofManager.CreateDOFs_2D(meshDOFsInfo, mesh_connectivity_data);

    Gedim::Output::PrintGenericMessage("Discrete Space with " + std::to_string(dofs_data.NumberDOFs) + " DOFs and " +
                                           std::to_string(dofs_data.NumberStrongs) + " STRONGs",
                                       true);

    Gedim::Profiler::StopTime("CreateDiscreteSpace");
    Gedim::Output::PrintStatusProgram("CreateDiscreteSpace");

    Polydim::examples::Elliptic_PCC_DFN::Assembler assembler;
    Polydim::examples::Elliptic_PCC_DFN::Assembler::Elliptic_PCC_DFN_Problem_Data assembler_data;

    double time_assembler = 0.0;
    double time_solver = 0.0;

    for (unsigned int i = 0; i < config.ComputationalTime(); i++)
    {
        const auto start_time_assembler = Gedim::Profiler::GetTime();
        Polydim::examples::Elliptic_PCC_DFN::Assembler assembler;
        assembler_data =
            assembler.Assemble(config, domains, mesh, meshGeometricData, meshDOFsInfo, dofs_data, reference_element_data, *test);

        const auto end_time_assembler = Gedim::Profiler::GetTime();

        time_assembler += Gedim::Profiler::ComputeTime(start_time_assembler, end_time_assembler);

        const auto start_time_solver = Gedim::Profiler::GetTime();

        if (dofs_data.NumberDOFs > 0)
        {
            Gedim::Eigen_CholeskySolver solver;
            solver.Initialize(assembler_data.globalMatrixA);
            solver.Solve(assembler_data.rightHandSide, assembler_data.solution);
        }

        const auto end_time_solver = Gedim::Profiler::GetTime();

        time_solver += Gedim::Profiler::ComputeTime(start_time_solver, end_time_solver);
    }

    time_assembler /= config.ComputationalTime();
    time_solver /= config.ComputationalTime();

    Gedim::Output::PrintGenericMessage("ComputeErrors...", true);
    Gedim::Profiler::StartTime("ComputeErrors");

    auto post_process_data =
        assembler.PostProcessSolution(config, domains, mesh, meshGeometricData, dofs_data, reference_element_data, assembler_data, *test);

    Gedim::Profiler::StopTime("ComputeErrors");
    Gedim::Output::PrintStatusProgram("ComputeErrors");

    Gedim::Output::PrintGenericMessage("ExportSolution...", true);
    Gedim::Profiler::StartTime("ExportSolution");

    Polydim::examples::Elliptic_PCC_DFN::program_utilities::export_solution(config,
                                                                            domains,
                                                                            mesh,
                                                                            meshGeometricData,
                                                                            dofs_data,
                                                                            assembler_data,
                                                                            post_process_data,
                                                                            time_assembler,
                                                                            time_solver,
                                                                            exportSolutionFolder,
                                                                            exportVtuFolder);

    Polydim::examples::Elliptic_PCC_DFN::program_utilities::export_dofs(config,
                                                                        mesh,
                                                                        meshGeometricData,
                                                                        meshDOFsInfo,
                                                                        dofs_data,
                                                                        reference_element_data,
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

        Polydim::examples::Elliptic_PCC_DFN::program_utilities::export_performance(config, performance, exportSolutionFolder);
    }

    Gedim::Profiler::StopTime("ComputeMethodPerformance");
    Gedim::Output::PrintStatusProgram("ComputeMethodPerformance");

    return 0;
}
