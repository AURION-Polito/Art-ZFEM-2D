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

#include "program_utilities.hpp"

#include "VTKUtilities.hpp"

namespace Polydim
{
namespace examples
{
namespace Elliptic_PCC_2D
{
namespace program_utilities
{
// ***************************************************************************
std::unique_ptr<Polydim::examples::Elliptic_PCC_2D::test::I_Test> create_test(const Polydim::examples::Elliptic_PCC_2D::Program_Configuration &config)
{
    switch (config.TestType())
    {
    case Polydim::examples::Elliptic_PCC_2D::test::Test_Types::Patch_Test:
        return std::make_unique<Polydim::examples::Elliptic_PCC_2D::test::Patch_Test>();
    case Polydim::examples::Elliptic_PCC_2D::test::Test_Types::Elliptic_Polynomial_Problem:
        return std::make_unique<Polydim::examples::Elliptic_PCC_2D::test::Elliptic_Polynomial_Problem>();
    case Polydim::examples::Elliptic_PCC_2D::test::Test_Types::Elliptic_Problem:
        return std::make_unique<Polydim::examples::Elliptic_PCC_2D::test::Elliptic_Problem>();
      case Polydim::examples::Elliptic_PCC_2D::test::Test_Types::Computational_Comparison:
          return std::make_unique<Polydim::examples::Elliptic_PCC_2D::test::Computational_Comparison>();
    default:
        throw std::runtime_error("Test type " + std::to_string((unsigned int)config.TestType()) + " not supported");
    }
}
// ***************************************************************************
void create_domain_mesh(const Polydim::examples::Elliptic_PCC_2D::Program_Configuration &config,
                        const Polydim::PDETools::Mesh::PDE_Mesh_Utilities::PDE_Domain_2D &domain,
                        Gedim::MeshMatricesDAO &mesh)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;

    switch (config.MeshGenerator())
    {
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Triangular:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Minimal:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Polygonal:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::Squared:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::RandomDistorted: {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::create_mesh_2D(geometryUtilities,
                                                                    meshUtilities,
                                                                    config.MeshGenerator(),
                                                                    domain,
                                                                    config.MeshMaxArea(),
                                                                    mesh);
    }
    break;
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::CsvImporter:
    case Polydim::PDETools::Mesh::PDE_Mesh_Utilities::MeshGenerator_Types_2D::OFFImporter: {
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::import_mesh_2D(meshUtilities, config.MeshGenerator(), config.MeshImportFilePath(), mesh);
    }
    break;
    default:
        throw std::runtime_error("MeshGenerator " + std::to_string((unsigned int)config.MeshGenerator()) + " not supported");
    }
    meshUtilities.ComputeCell1DCell2DNeighbours(mesh);
    Gedim::MeshUtilities::CheckMesh2DConfiguration configuration;
    configuration.Cell2D_CheckConvexity = false;
    meshUtilities.CheckMesh2D(configuration, geometryUtilities, mesh);
}
// ***************************************************************************
void make_mesh_triangular_by_internal_point(const std::vector<Eigen::Vector3d> &internal_points, Gedim::IMeshDAO &mesh)
{
    const unsigned int num_cells = mesh.Cell2DTotalNumber();

    unsigned int new_cell0D_id = mesh.Cell0DAppend(num_cells);
    for (unsigned int c = 0; c < num_cells; c++)
    {
        if (!mesh.Cell2DIsActive(c))
            continue;

        const unsigned int num_vertices = mesh.Cell2DNumberVertices(c);

        const auto &cell0Ds_idx = mesh.Cell2DVertices(c);
        const auto &cell1Ds_idx = mesh.Cell2DEdges(c);

        // add new vertex
        mesh.Cell0DInsertCoordinates(new_cell0D_id, internal_points.at(c));
        mesh.Cell0DSetState(new_cell0D_id, true);
        mesh.Cell0DSetMarker(new_cell0D_id, mesh.Cell2DMarker(c));

        // add new edges
        const unsigned int new_cell1D_starting = mesh.Cell1DAppend(num_vertices);
        std::vector<unsigned int> cell2Ds_new_edges(num_vertices);

        for (unsigned int t = 0; t < num_vertices; t++)
        {
            const unsigned int new_cell1D_index = new_cell1D_starting + t;

            cell2Ds_new_edges[t] = new_cell1D_index;
            mesh.Cell1DSetState(new_cell1D_index, true);
            mesh.Cell1DSetMarker(new_cell1D_index, mesh.Cell2DMarker(c));
            mesh.Cell1DInsertExtremes(new_cell1D_index, new_cell0D_id, cell0Ds_idx.at(t));
        }

        // add new faces
        const unsigned int new_cell2D_starting = mesh.Cell2DAppend(num_vertices);
        for (unsigned int t = 0; t < num_vertices; t++)
        {
            const unsigned int new_cell2D_index = new_cell2D_starting + t;

            mesh.Cell2DSetMarker(new_cell2D_index, mesh.Cell2DMarker(c));
            mesh.Cell2DSetState(new_cell2D_index, true);

            Eigen::MatrixXi vertices_edges(2, 3);
            vertices_edges.row(0) << new_cell0D_id, cell0Ds_idx.at(t), cell0Ds_idx.at((t + 1) % num_vertices);
            vertices_edges.row(1) << cell2Ds_new_edges[t], cell1Ds_idx.at(t), cell2Ds_new_edges[(t + 1) % num_vertices];

            mesh.Cell2DAddVerticesAndEdges(new_cell2D_index, vertices_edges);
            mesh.Cell2DInsertUpdatedCell2D(c, new_cell2D_index);
        }

        mesh.Cell2DSetState(c, false);
        new_cell0D_id++;
    }
}
// ***************************************************************************
Gedim::MeshUtilities::MeshGeometricData2D create_domain_mesh_geometric_properties(const Polydim::examples::Elliptic_PCC_2D::Program_Configuration &config,
                                                                                  const Gedim::MeshMatricesDAO &mesh)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    Gedim::MeshUtilities meshUtilities;

    Gedim::MeshUtilities::MeshGeometricData2D geometric_data =
        Polydim::PDETools::Mesh::PDE_Mesh_Utilities::compute_mesh_2D_geometry_data(geometryUtilities, meshUtilities, mesh);

    return geometric_data;
}
// ***************************************************************************
void export_solution(const Polydim::examples::Elliptic_PCC_2D::Program_Configuration &config,
                     const Gedim::MeshMatricesDAO &mesh,
                     const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                     const Polydim::examples::Elliptic_PCC_2D::Assembler::Elliptic_PCC_2D_Problem_Data &assembler_data,
                     const Polydim::examples::Elliptic_PCC_2D::Assembler::PostProcess_Data &post_process_data,
                     const double &time_assembler,
                     const double &time_solver,
                     const std::string &exportSolutionFolder,
                     const std::string &exportVtuFolder)
{
    const unsigned int Method_ID = static_cast<unsigned int>(config.MethodType());
    const unsigned int TEST_ID = static_cast<unsigned int>(config.TestType());

    {
        const char separator = ';';

        std::cout << "ProgramType" << separator;
        std::cout << "MethodType" << separator;
        std::cout << "MethodOrder" << separator;
        std::cout << "Cell2Ds" << separator;
        std::cout << "Dofs" << separator;
        std::cout << "Strongs" << separator;
        std::cout << "h" << separator;
        std::cout << "errorL2" << separator;
        std::cout << "errorH1" << separator;
        std::cout << "normL2" << separator;
        std::cout << "normH1" << separator;
        std::cout << "nnzA" << separator;
        std::cout << "residual" << separator;
        std::cout << "cond" << separator;
        std::cout << "TimeA" << separator;
        std::cout << "TimeS" << std::endl;

        std::cout.precision(2);
        std::cout << std::scientific << TEST_ID << separator;
        std::cout << std::scientific << Method_ID << separator;
        std::cout << std::scientific << config.MethodOrder() << separator;
        std::cout << std::scientific << mesh.Cell2DTotalNumber() << separator;
        std::cout << std::scientific << dofs_data.NumberDOFs << separator;
        std::cout << std::scientific << dofs_data.NumberStrongs << separator;
        std::cout << std::scientific << post_process_data.mesh_size << separator;
        std::cout << std::scientific << post_process_data.error_L2 << separator;
        std::cout << std::scientific << post_process_data.error_H1 << separator;
        std::cout << std::scientific << post_process_data.norm_L2 << separator;
        std::cout << std::scientific << post_process_data.norm_H1 << separator;
        std::cout << std::scientific << assembler_data.globalMatrixA.NonZeros() << separator;
        std::cout << std::scientific << post_process_data.residual_norm << separator;
        std::cout << std::scientific << post_process_data.conditioning << separator;
        std::cout << std::scientific << time_assembler << separator;
        std::cout << std::scientific << time_solver << std::endl;
    }

    {
        const char separator = ';';
        const std::string errorFileName = exportSolutionFolder + "/Errors_" + std::to_string(TEST_ID) + "_" +
                                          std::to_string(Method_ID) + "_" + std::to_string(config.MethodOrder()) + ".csv";
        const bool errorFileExists = Gedim::Output::FileExists(errorFileName);

        std::ofstream errorFile(errorFileName, std::ios_base::app | std::ios_base::out);
        if (!errorFileExists)
        {
            errorFile << "ProgramType" << separator;
            errorFile << "MethodType" << separator;
            errorFile << "MethodOrder" << separator;
            errorFile << "Cell2Ds" << separator;
            errorFile << "Dofs" << separator;
            errorFile << "Strongs" << separator;
            errorFile << "h" << separator;
            errorFile << "errorL2" << separator;
            errorFile << "errorH1" << separator;
            errorFile << "normL2" << separator;
            errorFile << "normH1" << separator;
            errorFile << "nnzA" << separator;
            errorFile << "residual" << separator;
            errorFile << "cond" << separator;
            errorFile << "TimeA" << separator;
            errorFile << "TimeS" << std::endl;
        }

        errorFile.precision(16);
        errorFile << std::scientific << TEST_ID << separator;
        errorFile << std::scientific << Method_ID << separator;
        errorFile << std::scientific << config.MethodOrder() << separator;
        errorFile << std::scientific << mesh.Cell2DTotalNumber() << separator;
        errorFile << std::scientific << dofs_data.NumberDOFs << separator;
        errorFile << std::scientific << dofs_data.NumberStrongs << separator;
        errorFile << std::scientific << post_process_data.mesh_size << separator;
        errorFile << std::scientific << post_process_data.error_L2 << separator;
        errorFile << std::scientific << post_process_data.error_H1 << separator;
        errorFile << std::scientific << post_process_data.norm_L2 << separator;
        errorFile << std::scientific << post_process_data.norm_H1 << separator;
        errorFile << std::scientific << assembler_data.globalMatrixA.NonZeros() << separator;
        errorFile << std::scientific << post_process_data.residual_norm << separator;
        errorFile << std::scientific << post_process_data.conditioning << separator;
        errorFile << std::scientific << time_assembler << separator;
        errorFile << std::scientific << time_solver << std::endl;

        errorFile.close();
    }

    if (config.PostProcess())
    {

        Gedim::VTKUtilities exporter;
        exporter.AddPolygons(mesh.Cell0DsCoordinates(),
                             mesh.Cell2DsVertices(),
                             {{"Numeric",
                               Gedim::VTPProperty::Formats::Points,
                               static_cast<unsigned int>(post_process_data.cell0Ds_numeric.size()),
                               post_process_data.cell0Ds_numeric.data()},
                              {"Exact",
                               Gedim::VTPProperty::Formats::Points,
                               static_cast<unsigned int>(post_process_data.cell0Ds_exact.size()),
                               post_process_data.cell0Ds_exact.data()},
                              {"ErrorL2",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(post_process_data.cell2Ds_error_L2.size()),
                               post_process_data.cell2Ds_error_L2.data()},
                              {"ErrorH1",
                               Gedim::VTPProperty::Formats::Cells,
                               static_cast<unsigned int>(post_process_data.cell2Ds_error_H1.size()),
                               post_process_data.cell2Ds_error_H1.data()}});

        exporter.Export(exportVtuFolder + "/Solution_" + std::to_string(TEST_ID) + "_" + std::to_string(Method_ID) +
                        "_" + std::to_string(config.MethodOrder()) + ".vtu");
    }
}
// ***************************************************************************
void export_dofs(const Polydim::examples::Elliptic_PCC_2D::Program_Configuration &config,
                 const Gedim::MeshMatricesDAO &mesh,
                 const Gedim::MeshUtilities::MeshGeometricData2D &mesh_geometric_data,
                 const Polydim::PDETools::DOFs::DOFsManager::MeshDOFsInfo &mesh_dofs_info,
                 const Polydim::PDETools::DOFs::DOFsManager::DOFsData &dofs_data,
                 const Polydim::examples::Elliptic_PCC_2D::Assembler::Elliptic_PCC_2D_Problem_Data &assembler_data,
                 const Polydim::examples::Elliptic_PCC_2D::Assembler::PostProcess_Data &post_process_data,
                 const std::string &exportVtuFolder)
{
    Gedim::GeometryUtilitiesConfig geometryUtilitiesConfig;
    geometryUtilitiesConfig.Tolerance1D = config.GeometricTolerance1D();
    geometryUtilitiesConfig.Tolerance2D = config.GeometricTolerance2D();
    Gedim::GeometryUtilities geometryUtilities(geometryUtilitiesConfig);

    const unsigned int Method_ID = static_cast<unsigned int>(config.MethodType());
    const unsigned int TEST_ID = static_cast<unsigned int>(config.TestType());
    const std::string file_path = exportVtuFolder + "/dofs_" + std::to_string(TEST_ID) + "_" +
                                  std::to_string(Method_ID) + "_" + std::to_string(config.MethodOrder()) + ".vtu";

    PDETools::LocalSpace_PCC_2D::export_dofs(geometryUtilities,
                                             mesh,
                                             mesh_geometric_data,
                                             mesh_dofs_info,
                                             dofs_data,
                                             assembler_data.rightHandSide,
                                             assembler_data.solution,
                                             assembler_data.solutionDirichlet,
                                             file_path);
}
// ***************************************************************************
void export_performance(const Polydim::examples::Elliptic_PCC_2D::Program_Configuration &config,
                        const Assembler::Performance_Data &performance_data,
                        const std::string &exportFolder)
{

    const char separator = ',';
    std::ofstream exporter;
    const unsigned int Method_ID = static_cast<unsigned int>(config.MethodType());
    const unsigned int TEST_ID = static_cast<unsigned int>(config.TestType());
    exporter.open(exportFolder + "/Cell2Ds_MethodPerformance_" + std::to_string(TEST_ID) + "_" +
                  std::to_string(Method_ID) + "_" + std::to_string(config.MethodOrder()) + ".csv");
    exporter.precision(16);

    if (exporter.fail())
        throw std::runtime_error("Error on mesh cell2Ds file");

    switch (config.MethodType())
    {
    case PDETools::LocalSpace_PCC_2D::MethodTypes::VEM_PCC_Inertia:
    case PDETools::LocalSpace_PCC_2D::MethodTypes::VEM_PCC_Ortho:
    case PDETools::LocalSpace_PCC_2D::MethodTypes::VEM_PCC: {
        exporter << "Cell2D_Index" << separator;
        exporter << "NumQuadPoints_Boundary" << separator;
        exporter << "NumQuadPoints_Internal" << separator;
        exporter << "PiNabla_Cond" << separator;
        exporter << "Pi0k_Cond" << separator;
        exporter << "Pi0km1_Cond" << separator;
        exporter << "PiNabla_Error" << separator;
        exporter << "Pi0k_Error" << separator;
        exporter << "Pi0km1_Error" << separator;
        exporter << "HCD_Error" << separator;
        exporter << "GBD_Error" << separator;
        exporter << "Stab_Error" << std::endl;

        for (unsigned int v = 0; v < performance_data.Cell2DsPerformance.size(); v++)
        {
            const auto &cell2D_performance = performance_data.Cell2DsPerformance[v].performance_data;

            exporter << std::scientific << v << separator;
            exporter << std::scientific << cell2D_performance.NumBoundaryQuadraturePoints << separator;
            exporter << std::scientific << cell2D_performance.NumInternalQuadraturePoints << separator;
            exporter << std::scientific << cell2D_performance.vem_analysis_data.PiNablaConditioning << separator;
            exporter << std::scientific << cell2D_performance.vem_analysis_data.Pi0kConditioning << separator;
            exporter << std::scientific << cell2D_performance.vem_analysis_data.Pi0km1Conditioning << separator;
            exporter << std::scientific << cell2D_performance.vem_analysis_data.ErrorPiNabla << separator;
            exporter << std::scientific << cell2D_performance.vem_analysis_data.ErrorPi0k << separator;
            exporter << std::scientific << cell2D_performance.vem_analysis_data.ErrorPi0km1 << separator;
            exporter << std::scientific << cell2D_performance.vem_analysis_data.ErrorHCD << separator;
            exporter << std::scientific << cell2D_performance.vem_analysis_data.ErrorGBD << separator;
            exporter << std::scientific << cell2D_performance.vem_analysis_data.ErrorStabilization << std::endl;
        }
    }
    break;
    case PDETools::LocalSpace_PCC_2D::MethodTypes::FEM_PCC:
        break;
    case PDETools::LocalSpace_PCC_2D::MethodTypes::ZFEM_PCC: {
        exporter << "Cell2D_Index" << separator;
        exporter << "NumQuadPoints_Internal" << separator;
        exporter << "consistency0" << separator;
        exporter << "consistency1" << std::endl;

        for (unsigned int v = 0; v < performance_data.Cell2DsPerformance.size(); v++)
        {

            const auto &cell2D_performance = performance_data.Cell2DsPerformance[v].performance_data;
            const double consistency_error =
                cell2D_performance.zfem_analysis_data.ConsistencyError.array().sqrt().matrix().maxCoeff();
            const double consistency_error_derivatives =
                (cell2D_performance.zfem_analysis_data.DerivativesConsistencyError[0] +
                 cell2D_performance.zfem_analysis_data.DerivativesConsistencyError[1])
                    .array()
                    .sqrt()
                    .matrix()
                    .maxCoeff();

            exporter << std::scientific << v << separator;
            exporter << std::scientific << cell2D_performance.NumInternalQuadraturePoints << separator;
            exporter << std::scientific << consistency_error << separator;
            exporter << std::scientific << consistency_error_derivatives << std::endl;
        }
    }
    break;
    default:
        throw std::runtime_error("not valid method type");
    }

    exporter.close();
}
// ***************************************************************************
} // namespace program_utilities
} // namespace Elliptic_PCC_2D
} // namespace examples
} // namespace Polydim
