# -*- coding: utf-8 -*-

_author_ = ""
_copyright_ = ""
_license_ = ""

import csv
import numpy as np


class Export:
    @staticmethod
    def evaluate_max_list(parameter_list, time_list):
        """Method for evaluating the max value of a time dependent
        simulation parameter along the time of max

        Args:
            parameter_list (list): time dependent simulation parameter
            time_list (list): time steps for the supplied parameter

        Returns:
            tuple: tuple containing the max value and its time
        """
        max_index = np.argmax(parameter_list)
        max_value = parameter_list[max_index]
        time_of_max = time_list[max_index]
        return max_value, time_of_max

    @staticmethod
    def evaluate_max_variables_list(time, variables):
        max_variable_list = []

        for data_variables in variables:
            max_variable_list.append(Export.evaluate_max_list(data_variables, time))

        return max_variable_list

    @staticmethod
    def raw_simulation_data_export(
        data, filepath: str, header_line, append: bool = False
    ):
        """Adapter method for direct export of simulation solution given
        by solve_ivp

        Args:
            data (list): solution data list to be exported
            filepath (str): filepath to the csv
            header_line (list): list of strings as csv header
            append (bool, optional): boolean option for appending or overwriting
            existing csv. Defaults to False.

        Returns:
            None
        """

        append_boolean = "a" if append else "w"

        with open(filepath, append_boolean, newline="") as file_data:
            solution_writer = csv.writer(file_data)
            solution_writer.writerow(header_line)

            for data_array in zip(
                *data,
            ):
                solution_writer.writerow(data_array)

            return None
