# -*- coding: utf-8 -*-

_author_ = ""
_copyright_ = ""
_license_ = ""

import csv
import numpy as np


class Export:
    @staticmethod
    def evaluate_max_list(list, time_list):
        max_index = np.argmax(list)
        max_value = list[max_index]
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

        append_boolean = "a" if append else "w"

        with open(filepath, append_boolean, newline="") as file_data:
            solution_writer = csv.writer(file_data)
            solution_writer.writerow(header_line)

            for data_array in zip(
                *data,
            ):
                solution_writer.writerow(data_array)

            return None
