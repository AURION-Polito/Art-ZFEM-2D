import os
import csv
import numpy as np
import matplotlib.pyplot as plt


def run_program(program_folder,
                program_path,
                run_folder,
                method_type,
                method_order,
                test_type,
                mesh_generator,
                num_ref,
                num_code_executions,
                mesh_max_area = 0.1,
                mesh_import_path = "./",
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
    program_parameters += " MeshImportFilePath:string={0}".format(mesh_import_path)
    program_parameters += " ComputationalTime:uint={0}".format(num_code_executions)

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
                errors_row.append(row[7])  # h
                errors_row.append(row[8])  # error L2
                errors_row.append(row[9])  # error H1
                errors_row.append(row[10]) # norm l2 numeric
                errors_row.append(row[11]) # norm h1 numeirc
                errors_row.append(row[3])  # cell2Ds
                errors_row.append(row[14])  # nnz
                errors_row.append(row[16])  # time a
                errors_row.append(row[17])  # time solver
            else:
                errors_row.append(float(row[4])) # dofs
                errors_row.append(float(row[5]))  # strongs
                errors_row.append(float(row[7])) # h
                errors_row.append(float(row[8])) # error L2
                errors_row.append(float(row[9])) # error H1
                errors_row.append(float(row[10])) # norm l2 numeric
                errors_row.append(float(row[11])) # norm h1 numeirc
                errors_row.append(float(row[3]))  # cell2Ds
                errors_row.append(float(row[14]))  # nnz
                errors_row.append(float(row[16]))  # time a
                errors_row.append(float(row[17])) # time solver
            errors.append(errors_row)
            counter += 1

    return errors


def test_errors(errors,
                method_order,
                tol):
    num_rows = len(errors)

    if num_rows == 2:
        print("Num. Ref. 1: ", abs(errors[1][3]) / abs(errors[1][5]), abs(errors[1][4]) / abs(errors[1][6]))
        assert abs(errors[1][3]) < tol * abs(errors[1][5])
        assert abs(errors[1][4]) < tol * abs(errors[1][6])
    else:
        errors = np.array(errors[1:])
        slope_L2 = np.polyfit(np.log(errors[:, 0]), np.log(errors[:, 3]), 1)[0]
        slope_H1 = np.polyfit(np.log(errors[:, 0]), np.log(errors[:, 4]), 1)[0]
        print("Num. Ref. ", str(num_rows-1), ": ", slope_L2*(-2), slope_H1*(-2))
        assert round(slope_L2*(-2)) == round(float(method_order + 1.0))
        assert round(slope_H1*(-2)) == round(float(method_order))

def loglog_slope_triangle(ax, x, y, alpha):


    # Parametri del triangolo
    p = -alpha  # slope desiderata
    x0 = x[-2] # punto di partenza del triangolo (asse x)
    y0 = y[-2] * 0.5  # punto di partenza del triangolo (asse y)
    scale = 4  # quanto grande deve essere il triangolo

    # Costruzione del triangolo
    x1 = x0 * scale
    y1 = y0 * (scale ** p)

    # Disegno del triangolo
    ax.loglog([x0, x1], [y1, y1], 'k-', linewidth=2)  # base orizzontale
    ax.loglog([x0, x0], [y1, y0], 'k-', linewidth=2) # lato verticale (discendente)
    ax.loglog([x0, x1], [y0, y1], 'k-', linewidth=2) # ipotenusa con slope p

    ax.text(x0*0.45, np.sqrt(y0 * y1), "{:<.2f}".format(-p), fontsize=30, ha='left', va='center')



def plot_errors(list_errors, list_errors_fem, method_order, method_types, plot_errors, plot_time):

    if plot_errors:
        fig, ax = plt.subplots(figsize=(12, 12))

        errors = list_errors_fem[method_order - 1]
        ax.plot(errors[:, 0], errors[:, 3], '-k^', linewidth=2, markersize=17,
                label="FEM")


        for h in range(len(method_types)):

            if method_types[h] == 1:
                errors = list_errors[h]
                ax.plot(errors[:, 0], errors[:, 3], '-ro', linewidth=2, markersize=15,
                        label="VEM")
            elif method_types[h] == 4:
                errors = list_errors[h]
                ax.plot(errors[:, 0], errors[:, 3], '-bs', linewidth=2, markersize=15,
                        label="Z-FEM")

                loglog_slope_triangle(ax, errors[:, 0], errors[:, 3], min(2.5, method_order + 1) * 0.5)

            else:
                raise ValueError("Not valid method type")

        errors = list_errors_fem[method_order - 1]
        ax.plot(errors[:, 0], errors[:, 4], '--k^', linewidth=2, markersize=14,
                label="FEM")


        for h in range(len(method_types)):
            if method_types[h] == 1:
                errors = list_errors[h]
                ax.plot(errors[:, 0], errors[:, 4], '--ro', linewidth=2, markersize=15,
                        label="VEM")
            elif method_types[h] == 4:
                errors = list_errors[h]
                ax.plot(errors[:, 0], errors[:, 4], '--bs', linewidth=2, markersize=15,
                        label="Z-FEM")

                loglog_slope_triangle(ax, errors[:, 0], errors[:, 4], min(1.5, method_order) * 0.5)
            else:
                raise ValueError("Not valid method type")

        # Get handles and labels
        handles, labels = plt.gca().get_legend_handles_labels()

        # Reorder them by row (instead of column)
        ll = [0, 1, 2]
        ordered_handles = [handles[i] for i in ll]
        ordered_labels = [labels[i] for i in ll]

        plt.legend(ordered_handles, ordered_labels, bbox_to_anchor=(0., 1.02, 1.0, 0.2), loc="lower left",
                   mode="expand", borderaxespad=0, ncol=3, fontsize=40)

        plt.xlabel('$N_{\\mathrm{dof}}$', fontsize=40)
        plt.xticks(fontsize=40)
        plt.yticks(fontsize=40)
        plt.xscale('log')
        plt.yscale('log')
        plt.grid(True, which="both", ls="--")
        plt.savefig(export_folder + "/{}_{}_{}_decay_plot.png".format(test_type, method_order, name_test),
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
                errors = list_errors[i-1]
            plt.bar(x + i * width, errors[:, 9], width,
                    label=f"{method} – Assembler")
            plt.bar(x + i * width, errors[:, 10], width,
                    bottom=errors[:, 9],
                    label=f"{method} – Solver")

        plt.yscale("log")
        plt.xticks(x + width, errors[:, 0])
        plt.xlabel('$N_{\\mathrm{dof}}$', fontsize=30)
        plt.ylabel("Tempo [ms] (log)")
        #plt.title("Confronto tempi (barplot stacked, asse y log)")
        plt.legend(ncol=3, fontsize=8)
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    program_folder = os.path.dirname(os.path.realpath(__file__))
    program_path = os.path.join(".", program_folder, "Elliptic_PCC_DFN")

    remove_folder = False

    export_folder = "integration_tests"
    os.system("rm -rf " + os.path.join(program_folder, export_folder))
    tol = 1.0e-12

    print("RUN TESTS...")

    test_type = 4
    mesh_generator = 1
    mesh_max_area = 0.0
    method_types = [0, 1, 2, 3]
    method_orders = [1, 2, 3]
    for method_type in method_types:
        for method_order in method_orders:
            export_path = run_program(program_folder,
                                      program_path,
                                      "Run_MG{0}".format(mesh_generator),
                                      method_type,
                                      method_order,
                                      test_type,
                                      mesh_generator,
                                      num_ref = 0,
                                      num_code_executions = 1,
                                      mesh_max_area=mesh_max_area)
            errors = import_errors(export_path, method_type, method_order, test_type)
            test_errors(errors,
                        method_order,
                        tol)

            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 4
    mesh_generator = 0
    mesh_max_area = 0.1
    method_types = [0]
    method_orders = [1, 2, 3]
    for method_type in method_types:
        for method_order in method_orders:
            export_path = run_program(program_folder,
                                      program_path,
                                      "Run_MG{0}".format(mesh_generator),
                                      method_type,
                                      method_order,
                                      test_type,
                                      mesh_generator,
                                      num_ref = 0,
                                      num_code_executions = 1,
                                      mesh_max_area=mesh_max_area)
            errors = import_errors(export_path, method_type, method_order, test_type)
            test_errors(errors,
                        method_order,
                        tol)

            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 5
    mesh_generator = 4
    mesh_import_paths = ["../../Mesh/2D/DFN_Frac_3/Triangle"]
    method_type = 0
    method_orders = np.arange(1, 4)
    list_errors_fem = []
    for method_order in method_orders:
        for mesh_size in [1, 2, 3, 4, 5]:
            mesh_import_path = mesh_import_paths[0] + "/M" + str(mesh_size)
            export_path = run_program(program_folder,
                                      program_path,
                                      "Run_MG{0}".format(mesh_generator),
                                      method_type,
                                      method_order,
                                      test_type,
                                      mesh_generator,
                                      num_ref=mesh_size,
                                      num_code_executions=5,
                                      mesh_import_path=mesh_import_path)

        errors = import_errors(export_path, method_type, method_order, test_type)

        if method_order == 1:
            test_errors(errors,
                        method_order,
                        tol)

        list_errors_fem.append(np.array(errors[1:]))

        if remove_folder:
            os.system("rm -rf " + os.path.join(program_folder, export_path))

    test_type = 5
    mesh_generator = 4
    mesh_import_paths = ["../../Mesh/2D/DFN_Frac_3/Triangle"]
    method_types = [1, 4]
    method_orders = np.arange(1, 4)
    name_test = "triangle_dfn"
    for method_order in method_orders:
        list_errors = []
        for method_type in method_types:
            for mesh_size in [1, 2, 3, 4, 5]:
                mesh_import_path = mesh_import_paths[0] + "/M" + str(mesh_size)
                export_path = run_program(program_folder,
                                          program_path,
                                          "Run_{1}_MG{0}".format(mesh_generator, name_test),
                                          method_type,
                                          method_order,
                                          test_type,
                                          mesh_generator,
                                          num_ref=mesh_size,
                                          num_code_executions=1,
                                          mesh_import_path=mesh_import_path)

            errors = import_errors(export_path, method_type, method_order, test_type)

            list_errors.append(np.array(errors[1:]))

            if method_order == 1:
                test_errors(errors,
                            method_order,
                            tol)

            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

        plot_errors(list_errors, list_errors_fem, method_order, method_types, False, False)

    test_type = 5
    mesh_generator = 4
    mesh_import_paths = ["../../Mesh/2D/DFN_Frac_3/Polygon"]
    method_types = [1, 4]
    method_orders = np.arange(1, 4)
    name_test = "polygonal_dfn"
    for method_order in method_orders:
        list_errors = []
        for method_type in method_types:
            for mesh_size in [1,2,3,4,5]:
                mesh_import_path = mesh_import_paths[0] + "/M" + str(mesh_size)
                export_path = run_program(program_folder,
                                        program_path,
                                        "Run_{1}_MG{0}".format(mesh_generator, name_test),
                                        method_type,
                                        method_order,
                                        test_type,
                                        mesh_generator,
                                        num_ref = mesh_size,
                                        num_code_executions = 1,
                                        mesh_import_path=mesh_import_path)

            errors = import_errors(export_path, method_type, method_order, test_type)

            if method_order == 1:
                test_errors(errors,
                            method_order,
                            tol)

            list_errors.append(np.array(errors[1:]))

            if remove_folder:
                os.system("rm -rf " + os.path.join(program_folder, export_path))

        plot_errors(list_errors, list_errors_fem, method_order, method_types, False, False)


    if remove_folder:
        os.system("rm -rf " + os.path.join(program_folder, export_folder))

    print("TESTS SUCCESS")
