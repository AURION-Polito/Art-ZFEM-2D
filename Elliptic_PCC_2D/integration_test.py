import os
import csv
import matplotlib.pyplot as plt
import numpy as np


def run_program(program_folder,
                program_path,
                run_folder,
                method_type,
                method_order,
                test_type,
                mesh_generator,
                num_ref,
                sub_triangulate,
                compute_conditioning,
                num_code_executions,
                mesh_max_area=0.1,
                mesh_import_path="./",
                post_process=True
                ):
    export_path = os.path.join(program_folder,
                               export_folder,
                               "{0}_TT{1}".format(
                                   run_folder,
                                   test_type),
                               "{0}_TT{1}_VT{2}".format(
                                   run_folder,
                                   test_type,
                                   method_type),
                               "{0}_TT{1}_VT{2}_VO{3}".format(
                                   run_folder,
                                   test_type,
                                   method_type,
                                   method_order))

    program_parameters = "MethodType:uint={0}".format(method_type)
    program_parameters += " MethodOrder:uint={0}".format(method_order)
    program_parameters += " ExportFolder:string={0}".format(export_path)
    program_parameters += " TestType:uint={0}".format(test_type)
    program_parameters += " MeshGenerator:uint={0}".format(mesh_generator)
    program_parameters += " MeshMaxArea:double={0}".format(mesh_max_area)
    program_parameters += " ComputeMethodPerformance:bool={0}".format(0)
    program_parameters += " ComputeConditionNumber:bool={0}".format(compute_conditioning)
    program_parameters += " MeshImportFilePath:string={0}".format(mesh_import_path)
    program_parameters += " SubTriangulate:bool={0}".format(sub_triangulate)
    program_parameters += " ComputationalTime:uint={0}".format(num_code_executions)
    program_parameters += " PostProcess:bool={0}".format(post_process)
    program_parameters += " GeometricTolerance1D:double={0}".format(1.0e-8)
    program_parameters += " GeometricTolerance2D:double={0}".format(1.0e-8)

    output_file = os.path.join(program_folder,
                               "terminal.log")

    run_label = "MethodType {0}".format(method_type)
    run_label += " MethodOrder {0}".format(method_order)
    run_label += " TestType {0}".format(test_type)
    run_label += " MeshGenerator {0}".format(mesh_generator)
    run_label += " NumRefinement {0}".format(num_ref)
    print("Run " + run_label + "...")
    os.system(program_path + " " + program_parameters + " > " + output_file)
    os.system("mv " + output_file + " " + export_path)
    print("Run SUCCESS")

    return export_path


def import_errors(export_path, method_type, method_order, test_type):
    errors_file = os.path.join(export_path,
                               "Solution",
                               "Errors_" + str(test_type) + "_" + str(method_type) + "_" + str(method_order) + ".csv")
    errors = []
    with open(errors_file, newline='') as csvfile:
        file_reader = csv.reader(csvfile, delimiter=';')
        data = list(file_reader)

        counter = 0
        for row in data:
            errors_row = []
            if counter == 0:
                errors_row.append(row[4])  # dofs
                errors_row.append(row[5])  # strongs
                errors_row.append(row[6])  # h
                errors_row.append(row[7])  # error L2
                errors_row.append(row[8])  # error H1
                errors_row.append(row[9])
                errors_row.append(row[10])
                errors_row.append(row[3])  # cell2Ds
                errors_row.append(row[11])  # nnz
                errors_row.append(row[13])  # cond
                errors_row.append(row[14])  # time a
                errors_row.append(row[15])  # time solver
            else:
                errors_row.append(float(row[4]))
                errors_row.append(float(row[5]))
                errors_row.append(float(row[6]))
                errors_row.append(float(row[7]))
                errors_row.append(float(row[8]))
                errors_row.append(float(row[9]))
                errors_row.append(float(row[10]))
                errors_row.append(float(row[3]))
                errors_row.append(float(row[11]))  # nnz
                errors_row.append(float(row[13]))  # cond
                errors_row.append(float(row[14]))  # time a
                errors_row.append(float(row[15]))  # time solver
            errors.append(errors_row)
            counter += 1

    return errors


def test_errors(errors,
                method_order,
                tol):
    num_rows = len(errors)

    if num_rows == 2:
        print("Num. Ref. 1: ", abs(errors[1][1]) / abs(errors[1][3]), abs(errors[1][2]) / abs(errors[1][4]))
    else:
        errors = np.array(errors[1:])
        slope_L2 = np.polyfit(np.log(errors[:, 2]), np.log(errors[:, 3]), 1)[0]
        slope_H1 = np.polyfit(np.log(errors[:, 2]), np.log(errors[:, 4]), 1)[0]
        print("Num. Ref. ", str(num_rows - 1), ": ", slope_L2, slope_H1, errors[:, 0], errors[:, 7])

        return slope_L2, slope_H1


def plot_errors(list_errors, list_errors_fem, list_errors_fem_2, method_order, method_types, plot_err, plot_time, plot_conditioning):
    if plot_err:
        fig, ax = plt.subplots(figsize=(12, 12))

        errors = list_errors_fem[method_order - 1]
        ax.plot(errors[:, 0], errors[:, 3], '-k^', linewidth=2, markersize=14,
                label="FEM")

        for h in range(len(method_types)):

            if method_types[h] == 1:
                errors = list_errors[h]
                ax.plot(errors[:, 0], errors[:, 3], '-ro', linewidth=2, markersize=12,
                        label="VEM")
            elif method_types[h] == 4:
                errors = list_errors[h]
                ax.plot(errors[:, 0], errors[:, 3], '-bs', linewidth=2, markersize=12,
                        label="Z-FEM")
            else:
                raise ValueError("Not valid method type")

        errors = list_errors_fem[method_order - 1]
        ax.plot(errors[:, 0], errors[:, 4], '--k^', linewidth=2, markersize=14,
                label="FEM")

        for h in range(len(method_types)):
            if method_types[h] == 1:
                errors = list_errors[h]
                ax.plot(errors[:, 0], errors[:, 4], '--ro', linewidth=2, markersize=12,
                        label="VEM")
            elif method_types[h] == 4:
                errors = list_errors[h]
                ax.plot(errors[:, 0], errors[:, 4], '--bs', linewidth=2, markersize=12,
                        label="Z-FEM")
            else:
                raise ValueError("Not valid method type")

        # Get handles and labels
        handles, labels = plt.gca().get_legend_handles_labels()

        # Reorder them by row (instead of column)
        ll = [0, 1, 2]
        ordered_handles = [handles[i] for i in ll]
        ordered_labels = [labels[i] for i in ll]

        plt.legend(ordered_handles, ordered_labels, bbox_to_anchor=(0., 1.02, 1.0, 0.2), loc="lower left",
                   mode="expand", borderaxespad=0, ncol=3, fontsize=30)

        plt.xlabel('$N_{\\mathrm{dof}}$', fontsize=30)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.xscale('log')
        plt.yscale('log')
        plt.grid(True, which="both", ls="--")
        plt.savefig(export_folder + "/{}_{}_{}_decay_plot.png".format(test_type, method_order, name_test),
                    bbox_inches='tight', dpi=300)
        plt.show()

    if plot_conditioning:
        fig, ax = plt.subplots(figsize=(12, 12))

        errors = list_errors_fem_2[method_order - 1]
        ax.plot(errors[:, 0], errors[:, 9], '-k^', linewidth=2, markersize=14,
                label="FEM")

        for h in range(len(method_types)):

            if method_types[h] == 1:
                errors = list_errors[h]
                ax.plot(errors[:, 0], errors[:, 9], '-ro', linewidth=2, markersize=12,
                        label="VEM")
            elif method_types[h] == 4:
                errors = list_errors[h]
                ax.plot(errors[:, 0], errors[:, 9], '-bs', linewidth=2, markersize=12,
                        label="Z-FEM")
            else:
                raise ValueError("Not valid method type")

        plt.legend(bbox_to_anchor=(0., 1.02, 1.0, 0.2), loc="lower left",
                   mode="expand", borderaxespad=0, ncol=3, fontsize=30)

        plt.xlabel('$N_{\\mathrm{dof}}$', fontsize=30)
        plt.ylabel('$\\mathrm{cond}(\\mathbf{A})$', fontsize=30)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
        plt.xscale('log')
        plt.yscale('log')
        plt.grid(True, which="both", ls="--")
        plt.savefig(export_folder + "/{}_{}_{}_conditioning.png".format(test_type, method_order, name_test),
                    bbox_inches='tight', dpi=300)
        plt.show()

    if plot_time:
        fig, ax = plt.subplots(figsize=(12, 12))
        errors = list_errors_fem[method_order - 1]
        x = np.arange(len(errors[:, 0]))
        width = 0.22

        methods = ["FEM", "VEM", "Z-FEM"]

        for i, method in enumerate(methods):
            if i == 0:
                errors = list_errors_fem[method_order - 1]
            else:
                errors = list_errors[i - 1]
            plt.bar(x + i * width, errors[:, 10], width,
                    label=f"{method} – Assembler")
            plt.bar(x + i * width, errors[:, 11], width,
                    bottom=errors[:, 10],
                    label=f"{method} – Solver")

        plt.yscale("log")
        plt.xticks(x + width, errors[:, 0])
        plt.xlabel('$N_{\\mathrm{dof}}$', fontsize=30)
        plt.ylabel("Tempo [ms] (log)")
        # plt.title("Confronto tempi (barplot stacked, asse y log)")
        plt.legend(ncol=3, fontsize=8)
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    program_folder = os.path.dirname(os.path.realpath(__file__))
    program_path = os.path.join(".", program_folder, "Elliptic_PCC_2D")

    remove_folder = False

    export_folder = "integration_tests"
    os.system("rm -rf " + os.path.join(program_folder, export_folder))
    tol = 1.0e-12

    plot_err = False
    plot_time = False
    plot_conditioning = False

    on_cluster = False

    print("RUN TESTS...")

    if not on_cluster:
        test_type = 3
        mesh_generator = 0
        method_type = 0
        mesh_max_areas = [0.05, 0.02, 0.01, 0.005, 0.0025]
        method_orders = [1, 2, 3, 4, 5, 6]
        list_errors_fem = []
        vv = 0
        for method_order in method_orders:
            num_ref = 0
            for mesh_max_area in mesh_max_areas:
                export_path = run_program(program_folder,
                                          program_path,
                                          "Run_MG{0}".format(mesh_generator),
                                          method_type,
                                          method_order,
                                          test_type,
                                          mesh_generator,
                                          num_ref,
                                          sub_triangulate=False,
                                          compute_conditioning=True,
                                          num_code_executions=1,
                                          mesh_max_area=mesh_max_area,
                                          mesh_import_path="./", )
                num_ref += 1

            errors = import_errors(export_path, method_type, method_order, test_type)
            test_errors(errors,
                        method_order,
                        tol)
            list_errors_fem.append(np.array(errors[1:]))

            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

        test_type = 3
        mesh_generator = 2
        method_type = 0
        mesh_max_areas = [0.1, 0.05, 0.02, 0.01, 0.005]
        method_orders = [1, 2, 3, 4, 5, 6]
        list_errors_fem_2 = []
        vv = 0
        for method_order in method_orders:
            num_ref = 0
            for mesh_max_area in mesh_max_areas:
                export_path = run_program(program_folder,
                                          program_path,
                                          "Run_MG{0}".format(mesh_generator),
                                          method_type,
                                          method_order,
                                          test_type,
                                          mesh_generator,
                                          num_ref,
                                          sub_triangulate=True,
                                          compute_conditioning=True,
                                          num_code_executions=1,
                                          mesh_max_area=mesh_max_area,
                                          mesh_import_path="./", )
                num_ref += 1

            errors = import_errors(export_path, method_type, method_order, test_type)
            test_errors(errors,
                        method_order,
                        tol)
            list_errors_fem_2.append(np.array(errors[1:]))

            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

        test_type = 3
        mesh_generator = 2
        method_types = [1, 4]
        mesh_max_areas = [0.1, 0.05, 0.02, 0.01, 0.005]
        method_orders = [1, 2, 3, 4, 5, 6]
        name_test = "voro"
        table_order = np.zeros([len(method_types) * 2, len(method_orders)])
        for method_order in method_orders:
            list_errors = []
            vv = 0
            for method_type in method_types:
                num_ref = 0
                for mesh_max_area in mesh_max_areas:
                    export_path = run_program(program_folder,
                                              program_path,
                                              "Run_MG{0}".format(mesh_generator),
                                              method_type,
                                              method_order,
                                              test_type,
                                              mesh_generator,
                                              num_ref,
                                              sub_triangulate=False,
                                              compute_conditioning=True,
                                              num_code_executions=1,
                                              mesh_max_area=mesh_max_area,
                                              mesh_import_path="./")
                    num_ref += 1

                errors = import_errors(export_path, method_type, method_order, test_type)
                slope_L2, slope_H1 = test_errors(errors,
                                                 method_order,
                                                 tol)

                table_order[2 * vv, method_order - 1] = slope_L2
                table_order[2 * vv + 1, method_order - 1] = slope_H1

                vv += 1

                list_errors.append(np.array(errors[1:]))

                if remove_folder:
                    os.system("rm -rf " + os.path.join(program_folder, export_path))

            plot_errors(list_errors, list_errors_fem, list_errors_fem_2, method_order, method_types, plot_err, plot_time,
                        plot_conditioning)

        with np.printoptions(precision=2):
            print(table_order)

        test_type = 3
        mesh_generator = 6
        method_type = 0
        mesh_max_areas = [0.1, 0.05, 0.02, 0.01, 0.005]
        method_orders = [1, 2, 3, 4, 5, 6]
        list_errors_fem_3 = []
        vv = 0
        for method_order in method_orders:
            num_ref = 0
            for mesh_max_area in mesh_max_areas:
                export_path = run_program(program_folder,
                                          program_path,
                                          "Run_MG{0}".format(mesh_generator),
                                          method_type,
                                          method_order,
                                          test_type,
                                          mesh_generator,
                                          num_ref,
                                          sub_triangulate=True,
                                          compute_conditioning=True,
                                          num_code_executions=1,
                                          mesh_max_area=mesh_max_area,
                                          mesh_import_path="./", )
                num_ref += 1

            errors = import_errors(export_path, method_type, method_order, test_type)
            test_errors(errors,
                        method_order,
                        tol)
            list_errors_fem_3.append(np.array(errors[1:]))

            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

        test_type = 3
        mesh_generator = 6
        method_types = [1, 4]
        mesh_max_areas = [0.1, 0.05, 0.02, 0.01, 0.005]
        method_orders = [1, 2, 3, 4, 5, 6]
        name_test = "rdquad"
        table_order = np.zeros([len(method_types) * 2, len(method_orders)])
        for method_order in method_orders:
            list_errors = []
            vv = 0
            for method_type in method_types:
                num_ref = 0
                for mesh_max_area in mesh_max_areas:
                    export_path = run_program(program_folder,
                                              program_path,
                                              "Run_MG{0}".format(mesh_generator),
                                              method_type,
                                              method_order,
                                              test_type,
                                              mesh_generator,
                                              num_ref,
                                              sub_triangulate=False,
                                              compute_conditioning=True,
                                              num_code_executions=5,
                                              mesh_max_area=mesh_max_area,
                                              mesh_import_path="./")
                    num_ref += 1

                errors = import_errors(export_path, method_type, method_order, test_type)

                slope_L2, slope_H1 = test_errors(errors,
                                                 method_order,
                                                 tol)

                table_order[2 * vv, method_order - 1] = slope_L2
                table_order[2 * vv + 1, method_order - 1] = slope_H1

                vv += 1

                list_errors.append(np.array(errors[1:]))

                if remove_folder:
                    os.system("rm -rf " + os.path.join(program_folder, export_path))

            plot_errors(list_errors, list_errors_fem, list_errors_fem_3, method_order, method_types, plot_err, plot_time,
                        plot_conditioning)

        with np.printoptions(precision=2):
            print(table_order)

        test_type = 3
        mesh_generator = 4
        method_type = 0
        mesh_max_areas = [program_folder + "/../../PolyDiM/Mesh/2D/StructuredConcave/StructuredConcave_3x3",
                          program_folder + "/../../PolyDiM/Mesh/2D/StructuredConcave/StructuredConcave_6x6",
                          program_folder + "/../../PolyDiM/Mesh/2D/StructuredConcave/StructuredConcave_9x9",
                          program_folder + "/../../PolyDiM/Mesh/2D/StructuredConcave/StructuredConcave_12x12",
                          program_folder + "/../../PolyDiM/Mesh/2D/StructuredConcave/StructuredConcave_15x15"]
        method_orders = [1, 2, 3, 4, 5, 6]
        list_errors_fem_4 = []
        vv = 0
        for method_order in method_orders:
            num_ref = 0
            for mesh_max_area in mesh_max_areas:
                export_path = run_program(program_folder,
                                          program_path,
                                          "Run_MG{0}".format(mesh_generator),
                                          method_type,
                                          method_order,
                                          test_type,
                                          mesh_generator,
                                          num_ref,
                                          sub_triangulate=True,
                                          compute_conditioning=True,
                                          num_code_executions=1,
                                          mesh_max_area=0.0,
                                          mesh_import_path=mesh_max_area, )
                num_ref += 1

            errors = import_errors(export_path, method_type, method_order, test_type)
            test_errors(errors,
                        method_order,
                        tol)
            list_errors_fem_4.append(np.array(errors[1:]))

            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

        test_type = 3
        mesh_generator = 4
        method_types = [1, 4]
        mesh_max_areas = [program_folder + "/../../PolyDiM/Mesh/2D/StructuredConcave/StructuredConcave_3x3",
                          program_folder + "/../../PolyDiM/Mesh/2D/StructuredConcave/StructuredConcave_6x6",
                          program_folder + "/../../PolyDiM/Mesh/2D/StructuredConcave/StructuredConcave_9x9",
                          program_folder + "/../../PolyDiM/Mesh/2D/StructuredConcave/StructuredConcave_12x12",
                          program_folder + "/../../PolyDiM/Mesh/2D/StructuredConcave/StructuredConcave_15x15"]
        method_orders = [1, 2, 3, 4, 5, 6]
        name_test = "conc_struct"
        table_order = np.zeros([len(method_types) * 2, len(method_orders)])
        for method_order in method_orders:
            list_errors = []
            vv = 0
            for method_type in method_types:
                num_ref = 0
                for mesh_max_area in mesh_max_areas:
                    export_path = run_program(program_folder,
                                              program_path,
                                              "Run_MG{0}".format(mesh_generator),
                                              method_type,
                                              method_order,
                                              test_type,
                                              mesh_generator,
                                              num_ref,
                                              sub_triangulate=False,
                                              compute_conditioning=True,
                                              num_code_executions=1,
                                              mesh_max_area=0.0,
                                              mesh_import_path=mesh_max_area)
                    num_ref += 1

                errors = import_errors(export_path, method_type, method_order, test_type)
                slope_L2, slope_H1 = test_errors(errors,
                                                 method_order,
                                                 tol)

                table_order[2 * vv, method_order - 1] = slope_L2
                table_order[2 * vv + 1, method_order - 1] = slope_H1
                vv += 1
                list_errors.append(np.array(errors[1:]))

                if remove_folder:
                    os.system("rm -rf " + os.path.join(program_folder, export_path))

            plot_errors(list_errors, list_errors_fem, list_errors_fem_4, method_order, method_types, plot_err, plot_time,
                        plot_conditioning)

        with np.printoptions(precision=2):
            print(table_order)

        test_type = 3
        mesh_generator = 0
        method_types = [1, 4]
        mesh_max_areas = [0.05, 0.02, 0.01, 0.005, 0.0025]
        method_orders = [1, 2, 3, 4, 5, 6]
        name_test = "triangle"
        table_order = np.zeros([len(method_types) * 2, len(method_orders)])
        for method_order in method_orders:
            list_errors = []
            vv = 0
            for method_type in method_types:
                num_ref = 0
                for mesh_max_area in mesh_max_areas:
                    export_path = run_program(program_folder,
                                              program_path,
                                              "Run_MG{0}".format(mesh_generator),
                                              method_type,
                                              method_order,
                                              test_type,
                                              mesh_generator,
                                              num_ref,
                                              sub_triangulate=False,
                                              compute_conditioning=True,
                                              num_code_executions=1,
                                              mesh_max_area=mesh_max_area,
                                              mesh_import_path="./")
                    num_ref += 1

                errors = import_errors(export_path, method_type, method_order, test_type)
                slope_L2, slope_H1 = test_errors(errors,
                                                 method_order,
                                                 tol)

                table_order[2 * vv, method_order - 1] = slope_L2
                table_order[2 * vv + 1, method_order - 1] = slope_H1
                vv += 1
                list_errors.append(np.array(errors[1:]))

                if remove_folder:
                    os.system("rm -rf " + os.path.join(program_folder, export_path))

            plot_errors(list_errors, list_errors_fem, method_order, method_types, plot_err, plot_time,
                        plot_conditioning)

        with np.printoptions(precision=2):
            print(table_order)

    if on_cluster:
        test_type = 3
        mesh_generator = 2
        method_type = 0
        mesh_max_areas = [0.005, 0.002, 0.001, 0.0005, 0.0003]
        method_orders = [6]
        list_errors_fem = []
        vv = 0
        for method_order in method_orders:
            num_ref = 0
            for mesh_max_area in mesh_max_areas:
                export_path = run_program(program_folder,
                                          program_path,
                                          "Run_MG{0}".format(mesh_generator),
                                          method_type,
                                          method_order,
                                          test_type,
                                          mesh_generator,
                                          num_ref,
                                          sub_triangulate=True,
                                          compute_conditioning=True,
                                          num_code_executions=5,
                                          mesh_max_area=mesh_max_area,
                                          mesh_import_path="./", )
                num_ref += 1

            errors = import_errors(export_path, method_type, method_order, test_type)

            list_errors_fem.append(np.array(errors[1:]))

            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

        method_types = [1, 4]
        name_test = "voro_time"
        for method_order in method_orders:
            list_errors = []
            vv = 0
            for method_type in method_types:
                num_ref = 0
                for mesh_max_area in mesh_max_areas:
                    export_path = run_program(program_folder,
                                              program_path,
                                              "Run_MG{0}".format(mesh_generator),
                                              method_type,
                                              method_order,
                                              test_type,
                                              mesh_generator,
                                              num_ref,
                                              sub_triangulate=False,
                                              compute_conditioning=True,
                                              num_code_executions=5,
                                              mesh_max_area=mesh_max_area,
                                              mesh_import_path="./")
                    num_ref += 1

                errors = import_errors(export_path, method_type, method_order, test_type)

                vv += 1

                if remove_folder:
                    os.system("rm -rf " + os.path.join(program_folder, export_path))

    if remove_folder:
        os.system("rm -rf " + os.path.join(program_folder, export_folder))

    print("TESTS SUCCESS")