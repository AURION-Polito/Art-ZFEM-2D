import os
import csv
import matplotlib.pyplot as plt
import numpy as np


def run_program(program_folder,
                export_folder,
                program_path,
                run_folder,
                method_type,
                method_order,
                test_type,
                mesh_generator,
                num_ref,
                sub_triangulate,
                theta_parameter,
                time_partition_type,
                num_time_steps,
                mesh_max_area=0.1,
                mesh_import_path="./",
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
    program_parameters += " ComputeMethodPerformance:bool={0}".format(False)
    program_parameters += " MeshImportFilePath:string={0}".format(mesh_import_path)
    program_parameters += " SubTriangulate:bool={0}".format(sub_triangulate)
    program_parameters += " PostProcess:bool={0}".format(False)
    program_parameters += " GeometricTolerance1D:double={0}".format(1.0e-12)
    program_parameters += " GeometricTolerance2D:double={0}".format(1.0e-14)
    program_parameters += " ThetaParameter:double={0}".format(theta_parameter)
    program_parameters += " TimePartition:uint={0}".format(time_partition_type)
    program_parameters += " NumTimeStep:uint={0}".format(num_time_steps)

    output_file = os.path.join(program_folder,
                               "terminal.log")

    run_label = "MethodType {0}".format(method_type)
    run_label += " MethodOrder {0}".format(method_order)
    run_label += " TestType {0}".format(test_type)
    run_label += " MeshGenerator {0}".format(mesh_generator)
    run_label += " NumRefinement {0}".format(num_ref)
    run_label += " ThetaParameter {0}".format(theta_parameter)
    print("Run " + run_label + "...")
    os.system(program_path + " " + program_parameters + " > " + output_file)
    os.system("mv " + output_file + " " + str(export_path))
    print("Run SUCCESS")

    return export_path


def import_errors(export_path, method_type, method_order, test_type, num_time_steps):
    errors_file = os.path.join(export_path,
                               "Solution",
                               "Errors_" + str(test_type) + "_" + str(method_type) + "_" + str(method_order) + "_" + str(num_time_steps) + ".csv")
    errors = []
    with open(errors_file, newline='') as csvfile:
        file_reader = csv.reader(csvfile, delimiter=';')
        data = list(file_reader)

        counter = 0
        for row in data:
            errors_row = []
            if counter == 0:
                errors_row.append(row[2])  # delta time
                errors_row.append(row[7])  # dof
                errors_row.append(row[8])  # strong
                errors_row.append(row[9])  # h
                errors_row.append(row[10])  # error L2
                errors_row.append(row[11])  # error H1
                errors_row.append(row[12]) # norm l2
                errors_row.append(row[13]) # norm h1
            else:
                errors_row.append(float(row[2]))  # delta time
                errors_row.append(float(row[7]))  # dof
                errors_row.append(float(row[8]))  # strong
                errors_row.append(float(row[9]))  # h
                errors_row.append(float(row[10]))  # error L2
                errors_row.append(float(row[11]))  # error H1
                errors_row.append(float(row[12])) # norm l2
                errors_row.append(float(row[13])) # norm h1
            errors.append(errors_row)
            counter += 1

    return errors


def check_errors_space(errors,
                       method_order,
                       _test_type,
                       tol):
    num_rows = len(errors)

    if num_rows == 2:
        print("Num. Ref. 1: ", abs(errors[1][4]) / abs(errors[1][6]), abs(errors[1][5]) / abs(errors[1][7]))

        assert abs(errors[1][4]) / abs(errors[1][6]) < tol
        assert abs(errors[1][5]) / abs(errors[1][7]) < tol

        return None
    else:
        errors = np.array(errors[1:])
        slope_l2 = float(np.polyfit(np.log(errors[:, 3]), np.log(errors[:, 4]), 1)[0])
        slope_h1 = float(np.polyfit(np.log(errors[:, 3]), np.log(errors[:, 5]), 1)[0])
        print("Num. Ref. ", str(num_rows - 1), ": ", slope_l2, slope_h1, errors[:, 0], errors[:, 7])

        assert round(slope_l2) >= round(float(method_order + 1.0))
        assert round(slope_h1) >= round(float(method_order))

        return slope_l2, slope_h1


def loglog_slope_triangle(ax, x, y, alpha):
    # triangle
    p = -alpha  # slope
    x0 = x[-2]  # origin (x-axis)
    y0 = y[-2] * 0.5  # origin (y-axis)
    scale = 2.5  # triangle scaling

    x1 = x0 * scale
    y1 = y0 * (scale ** p)

    # Plot triangle
    ax.loglog([x0, x1], [y1, y1], 'k-', linewidth=2)  # base
    ax.loglog([x0, x0], [y1, y0], 'k-', linewidth=2)  # height
    ax.loglog([x0, x1], [y0, y1], 'k-', linewidth=2)  # hypotenuse with slope p

    ax.text(x0 * 0.6, np.sqrt(y0 * y1), "{:<.2f}".format(-p), fontsize=30, ha='left', va='center')


def plot_errors(export_folder, list_errors, list_errors_fem, method_order, method_types, test_type, name_test, plot_err,
                plot_time, plot_conditioning):
    if plot_err:
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

                loglog_slope_triangle(ax, errors[:, 0], errors[:, 3], (method_order + 1) * 0.5)
            else:
                raise ValueError("Not valid method type")

        errors = list_errors_fem[method_order - 1]
        ax.plot(errors[:, 0], errors[:, 4], '--k^', linewidth=2, markersize=17,
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
                loglog_slope_triangle(ax, errors[:, 0], errors[:, 4], method_order * 0.5)
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

    if plot_conditioning:
        fig, ax = plt.subplots(figsize=(12, 12))

        errors = list_errors_fem[method_order - 1]
        ax.plot(errors[:, 0], errors[:, 9], '-k^', linewidth=2, markersize=17,
                label="FEM")

        for h in range(len(method_types)):

            if method_types[h] == 1:
                errors = list_errors[h]
                ax.plot(errors[:, 0], errors[:, 9], '-ro', linewidth=2, markersize=15,
                        label="VEM")
            elif method_types[h] == 4:
                errors = list_errors[h]
                ax.plot(errors[:, 0], errors[:, 9], '-bs', linewidth=2, markersize=15,
                        label="Z-FEM")
            else:
                raise ValueError("Not valid method type")

        plt.legend(bbox_to_anchor=(0., 1.02, 1.0, 0.2), loc="lower left",
                   mode="expand", borderaxespad=0, ncol=3, fontsize=40)

        plt.xlabel('$N_{\\mathrm{dof}}$', fontsize=40)
        plt.ylabel('$\\mathrm{cond}(\\mathbf{A})$', fontsize=40)
        plt.xticks(fontsize=40)
        plt.yticks(fontsize=40)
        plt.xscale('log')
        plt.yscale('log')
        plt.grid(True, which="both", ls="--")
        plt.savefig(export_folder + "/{}_{}_{}_conditioning.png".format(test_type, method_order, name_test),
                    bbox_inches='tight', dpi=300)
        plt.show()

    if plot_time:
        fig, ax = plt.subplots(figsize=(12, 12))

        errors = list_errors_fem[method_order - 1]
        ax.plot(errors[:, 0], errors[:, 10], ':k^', linewidth=2, markersize=17,
                label="FEM")

        ax.plot(errors[:, 0], errors[:, 11], '--k^', linewidth=2, markersize=17,
                label="FEM")

        ax.plot(errors[:, 0], errors[:, 10] + errors[:, 11], '-k^', linewidth=2, markersize=17,
                label="FEM")

        for h in range(len(method_types)):

            if method_types[h] == 1:
                errors = list_errors[h]
                ax.plot(errors[:, 0], errors[:, 10], ':ro', linewidth=2, markersize=15,
                        label="VEM")

                ax.plot(errors[:, 0], errors[:, 11], '--ro', linewidth=2, markersize=15,
                        label="VEM")

                ax.plot(errors[:, 0], errors[:, 10] + errors[:, 11], '-ro', linewidth=2, markersize=15,
                        label="VEM")
            elif method_types[h] == 4:
                errors = list_errors[h]
                ax.plot(errors[:, 0], errors[:, 10], ':bs', linewidth=2, markersize=15,
                        label="Z-FEM")

                ax.plot(errors[:, 0], errors[:, 11], '--bs', linewidth=2, markersize=15,
                        label="Z-FEM")

                ax.plot(errors[:, 0], errors[:, 10] + errors[:, 11], '-bs', linewidth=2, markersize=15,
                        label="Z-FEM")
            else:
                raise ValueError("Not valid method type")

        # Get handles and labels
        handles, labels = plt.gca().get_legend_handles_labels()

        # Reorder them by row (instead of column)
        ll = [2, 5, 8]
        ordered_handles = [handles[i] for i in ll]
        ordered_labels = [labels[i] for i in ll]

        plt.legend(ordered_handles, ordered_labels, bbox_to_anchor=(0., 1.02, 1.0, 0.2), loc="lower left",
                   mode="expand", borderaxespad=0, ncol=3, fontsize=40)

        plt.xscale("log")
        plt.yscale("log")
        plt.xticks(fontsize=40)
        plt.yticks(fontsize=40)
        plt.xlabel('$N_{\\mathrm{dof}}$', fontsize=40)
        plt.ylabel("Time [s]", fontsize=40)
        plt.grid(True, which="both", ls="--")
        plt.savefig(export_folder + "/{}_{}_{}_time.png".format(test_type, method_order, name_test),
                    bbox_inches='tight', dpi=300)
        plt.show()


def main():
    program_folder = os.path.dirname(os.path.realpath(__file__))
    program_path = os.path.join(".", program_folder, "Parabolic_PCC_2D")

    remove_folder = True

    export_folder = "integration_tests"
    os.system("rm -rf " + os.path.join(program_folder, export_folder))
    tol = 1.0e-8

    print("RUN TESTS...")

    test_type = 1
    mesh_generators = [0, 5]
    method_types = [0, 1, 4]
    mesh_max_areas = [0.1]
    method_orders = [1, 2, 3, 4]
    theta_parameters = [1.0, 0.0, 0.5]
    list_errors_fem = []
    num_time_steps = 1
    for method_type in method_types:
        for method_order in method_orders:
            for theta_parameter in theta_parameters:
                for mesh_generator in mesh_generators:
                    num_ref = 0
                    for mesh_max_area in mesh_max_areas:
                        export_path = run_program(program_folder,
                                                  export_folder,
                                                  program_path,
                                                  "Run_MG{0}".format(mesh_generator),
                                                  method_type,
                                                  method_order,
                                                  test_type,
                                                  mesh_generator,
                                                  num_ref,
                                                  num_time_steps=num_time_steps,
                                                  theta_parameter=theta_parameter,
                                                  time_partition_type=1,
                                                  sub_triangulate=False,
                                                  mesh_max_area=mesh_max_area,
                                                  mesh_import_path="./")
                    num_ref += 1

                    errors = import_errors(export_path, method_type, method_order, test_type, num_time_steps)
                    check_errors_space(errors,
                                       method_order,
                                       test_type,
                                       tol)
                    list_errors_fem.append(np.array(errors[1:]))

                    if remove_folder:
                        os.system("rm -rf " + os.path.join(str(program_folder), str(export_path)))

    test_type = 1
    mesh_generators = [2]
    method_types = [1, 4]
    mesh_max_areas = [0.1]
    method_orders = [1, 2, 3, 4]
    theta_parameters = [1.0, 0.0, 0.5]
    list_errors_fem = []
    num_time_steps = 1
    for method_type in method_types:
        for method_order in method_orders:
            for theta_parameter in theta_parameters:
                for mesh_generator in mesh_generators:
                    num_ref = 0
                    for mesh_max_area in mesh_max_areas:
                        export_path = run_program(program_folder,
                                                  export_folder,
                                                  program_path,
                                                  "Run_MG{0}".format(mesh_generator),
                                                  method_type,
                                                  method_order,
                                                  test_type,
                                                  mesh_generator,
                                                  num_ref,
                                                  num_time_steps=num_time_steps,
                                                  theta_parameter=theta_parameter,
                                                  time_partition_type=1,
                                                  sub_triangulate=False,
                                                  mesh_max_area=mesh_max_area,
                                                  mesh_import_path="./")
                    num_ref += 1

                    errors = import_errors(export_path, method_type, method_order, test_type, num_time_steps)
                    check_errors_space(errors,
                                       method_order,
                                       test_type,
                                       tol)
                    list_errors_fem.append(np.array(errors[1:]))

                    if remove_folder:
                        os.system("rm -rf " + os.path.join(str(program_folder), str(export_path)))

    if remove_folder:
        os.system("rm -rf " + os.path.join(program_folder, export_folder))

    print("TESTS SUCCESS")


if __name__ == "__main__":
    main()