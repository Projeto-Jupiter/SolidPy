import sys, ast, csv, json, os
import numpy as np

from solidpy import (
    Bates,
    Star,
    CustomGeometry,
    Motor,
    Propellant,
    Environment,
    Rail,
    BurnSimulation,
    Export,
)
import StructuralCasing

from solidpy.Burn import BurnExport

import matplotlib
matplotlib.use('Qt5Agg')

from PyQt5 import QtGui
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (
    QApplication,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QPushButton,
    QRadioButton,
    QCheckBox,
    QStackedLayout,
    QVBoxLayout,
    QTabWidget,
    QGridLayout,
    QWidget,
    QLineEdit,
    QComboBox,
    QFileDialog,
    QDialog,
    QMessageBox,
)

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.figure import Figure

class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, parent=None, width=5, height=4, dpi=100):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowIcon(QtGui.QIcon('logo.ico'))
        self.setWindowIconText("logo")
        self.setWindowTitle("SolidPy GUI")

        #creating plots
        plot_thrust = MplCanvas(self, width=15, height=8)
        plot_pressure = MplCanvas(self, width=15, height=8)
        plot_kn = MplCanvas(self, width=15, height=8)
        plot_mass = MplCanvas(self, width=15, height=8)
        casing_output_widget = QWidget()
        casing_output_layout = QGridLayout(casing_output_widget)
        plot_bolts_number = MplCanvas(self, width=15, height=6)
        plot_bolts_type = MplCanvas(self, width=15, height=6)
        plot_wall_thickness = MplCanvas(self, width=15, height=6)


        #creating tabs
        information_tabs = QTabWidget()
    
        Grain_Widget = QWidget()
        Motor_Widget = QWidget()
        Bolts_Widget = QWidget()
        Casing_Widget = QWidget()


        #creating tabs layouts

        #grain tab
        Grain_Layout = QGridLayout(Grain_Widget)

        Grain_Layout.addWidget(QLabel('Initial inner diameter(mm):'),0,0,Qt.AlignRight)
        self.initial_inner_diameter_input = QLineEdit('31.92')
        Grain_Layout.addWidget(self.initial_inner_diameter_input,0,1,Qt.AlignLeft)

        Grain_Layout.addWidget(QLabel('Outer Diameter(mm):'),1,0,Qt.AlignRight)
        self.outer_diameter_input = QLineEdit('71.92')
        Grain_Layout.addWidget(self.outer_diameter_input,1,1,Qt.AlignLeft)

        Grain_Layout.addWidget(QLabel('Height(mm):'),2,0,Qt.AlignRight)
        self.height_input = QLineEdit('123')
        Grain_Layout.addWidget(self.height_input,2,1,Qt.AlignLeft)

        Grain_Layout.addWidget(QLabel('Number:'),3,0,Qt.AlignRight)
        self.number_input = QLineEdit('4')
        Grain_Layout.addWidget(self.number_input,3,1,Qt.AlignLeft)

        #create combobox for propellant selection
        Grain_Layout.addWidget(QLabel('Propellant:'),4,0,Qt.AlignRight)
        self.propellant_combobox = QComboBox()
        #read the names of the files in the propellants data folder
        propellants = os.listdir(os.getcwd() + '\GUI\propellants data')
        #add the names to the combobox without the .json extension
        for propellant in propellants:
            self.propellant_combobox.addItem(propellant[:-5])
        Grain_Layout.addWidget(self.propellant_combobox,4,1,Qt.AlignLeft)





        #motor tab
        Motor_Layout = QGridLayout(Motor_Widget)

        Motor_Layout.addWidget(QLabel('Chamber Inner Diameter(mm):'),0,0,Qt.AlignRight)
        self.chamber_inner_diameter_input = QLineEdit('77.92')
        Motor_Layout.addWidget(self.chamber_inner_diameter_input,0,1,Qt.AlignLeft)

        Motor_Layout.addWidget(QLabel('Nozzle Throat Diameter(mm):'),1,0,Qt.AlignRight)
        self.nozzle_throat_diameter_input = QLineEdit('17.5')
        Motor_Layout.addWidget(self.nozzle_throat_diameter_input,1,1,Qt.AlignLeft)

        Motor_Layout.addWidget(QLabel('Nozzle Exit Diameter(mm):'),2,0,Qt.AlignRight)
        self.nozzle_exit_diameter_input = QLineEdit('44.44')
        Motor_Layout.addWidget(self.nozzle_exit_diameter_input,2,1,Qt.AlignLeft)

        Motor_Layout.addWidget(QLabel('Nozzle Angle(deg):'),3,0,Qt.AlignRight)
        self.nozzle_angle_input = QLineEdit('15')
        Motor_Layout.addWidget(self.nozzle_angle_input,3,1,Qt.AlignLeft)

        Motor_Layout.addWidget(QLabel('Chamber Lenght(mm):'),5,0,Qt.AlignRight)
        self.chamber_lenght_input = QLineEdit('600')
        Motor_Layout.addWidget(self.chamber_lenght_input,5,1,Qt.AlignLeft)

        #bolts tab
        Bolts_Layout = QGridLayout(Bolts_Widget)

        Bolts_Layout.addWidget(QLabel('Ultimate Tensile Strenght(MPa):'),0,0,Qt.AlignRight)
        self.ultimate_tensile_strength_input = QLineEdit('515')
        Bolts_Layout.addWidget(self.ultimate_tensile_strength_input,0,1,Qt.AlignLeft)

        Bolts_Layout.addWidget(QLabel('Bolts number:'),7,0,Qt.AlignRight)
        self.bolts_number_input = QLineEdit('12')
        Bolts_Layout.addWidget(self.bolts_number_input,7,1,Qt.AlignLeft)

        Bolts_Layout.addWidget(QLabel('Bolt type:'),8,0,Qt.AlignRight)
        self.bolt_types_widget = QComboBox()
        self.bolt_types_widget.addItems(['M5x0.5', 'M5x0.8', 'M5.5x0.5', 'M6x0.75', 'M6x1', 'M7x0.75', 'M7x1', 'M8x0.75', 'M8x1', 'M8x1.25'])
        Bolts_Layout.addWidget(self.bolt_types_widget,8,1,Qt.AlignLeft)
        self.minor_bolt_diameter_dict = {'M5x0.8':4.134,'M5x0.5':4.459,'M5.5x0.5':4.959,'M6x0.75':5.188,'M6x1':4.917,'M7x0.75':6.188,'M7x1':5.917,'M8x0.75':7.188,'M8x1':6.917,'M8x1.25':6.647} 
        self.major_bolt_diameter_dict = {'M5x0.5':5,'M5x0.8':5,'M5.5x0.5':5.5,'M6x0.75':6,'M6x1':6,'M7x0.75':7,'M7x1':7,'M8x0.75':8,'M8x1':8,'M8x1.25':8}
        self.bolt_types_widget.setCurrentIndex(4)

        #casing input tab
        Casing_Layout = QGridLayout(Casing_Widget)

        Casing_Layout.addWidget(QLabel('Yield Strength(MPa):'),1,0,Qt.AlignRight)
        self.yield_strength_input = QLineEdit('172')
        Casing_Layout.addWidget(self.yield_strength_input,1,1,Qt.AlignLeft)

        Casing_Layout.addWidget(QLabel('Shear Strength(MPa):'),2,0,Qt.AlignRight)
        self.shear_modulus_input = QLineEdit('138')
        Casing_Layout.addWidget(self.shear_modulus_input,2,1,Qt.AlignLeft)

        Casing_Layout.addWidget(QLabel('Bearing Strength(MPa):'),3,0,Qt.AlignRight)
        self.bearing_strength_input = QLineEdit('434')
        Casing_Layout.addWidget(self.bearing_strength_input,3,1,Qt.AlignLeft)

        Casing_Layout.addWidget(QLabel('Casing inner diameter(mm):'),4,0,Qt.AlignRight)
        self.casing_inner_diameter_input = QLineEdit('setted in motor tab')
        self.casing_inner_diameter_input.setDisabled(True)
        Casing_Layout.addWidget(self.casing_inner_diameter_input,4,1,Qt.AlignLeft)
        
        Casing_Layout.addWidget(QLabel('Wall thickness(mm):'),5,0,Qt.AlignRight)
        self.wall_thickness_input = QLineEdit('5.49')
        Casing_Layout.addWidget(self.wall_thickness_input,5,1,Qt.AlignLeft)

        Casing_Layout.addWidget(QLabel('Edge distance(mm):'),6,0,Qt.AlignRight)
        self.edge_distance_input = QLineEdit('15')
        Casing_Layout.addWidget(self.edge_distance_input,6,1,Qt.AlignLeft)


        #casing output tab
        casing_output_layout.addWidget(QLabel('Bolt Shear:'),0,0,Qt.AlignRight)
        self.bolt_shear_output = QLabel('-')
        casing_output_layout.addWidget(self.bolt_shear_output,0,1)
        casing_output_layout.addWidget(QLabel('Bolt Tear Out:'),1,0,Qt.AlignRight)
        self.bolt_tear_out_output = QLabel('-')
        casing_output_layout.addWidget(self.bolt_tear_out_output,1,1)
        casing_output_layout.addWidget(QLabel('Casing Tensile:'),0,2,Qt.AlignRight)
        self.casing_tensile_output = QLabel('-')
        casing_output_layout.addWidget(self.casing_tensile_output,0,3)
        casing_output_layout.addWidget(QLabel('Bearing:'),1,2,Qt.AlignRight)
        self.bearing_output = QLabel('-')
        casing_output_layout.addWidget(self.bearing_output,1,3)
        casing_output_layout.addWidget(QLabel('Casing Circumferential:'),0,4,Qt.AlignRight)
        self.casing_circumferential_output = QLabel('-')
        casing_output_layout.addWidget(self.casing_circumferential_output,0,5)
        casing_output_plots = QHBoxLayout()
        casing_output_layout.addLayout(casing_output_plots,2,0,1,6)
        casing_output_plots.addWidget(plot_bolts_number)
        casing_output_plots.addWidget(plot_bolts_type)
        casing_output_plots.addWidget(plot_wall_thickness)


        #input tabs
        information_tabs.addTab(Grain_Widget, 'Grain')
        information_tabs.addTab(Motor_Widget, 'Motor')
        information_tabs.addTab(Bolts_Widget, 'Bolts')
        information_tabs.addTab(Casing_Widget, 'Casing')

        #output tabs
        plot_tabs = QTabWidget()
        plot_tabs.addTab(plot_thrust,'Thrust')
        plot_tabs.addTab(plot_pressure,'Pressure')
        # plot_tabs.addTab(plot_kn,'Kn')
        # plot_tabs.addTab(plot_mass,'Mass')
        plot_tabs.addTab(casing_output_widget,'Casing Structural Analysis')

        pagelayout = QHBoxLayout()
        input_layout = QVBoxLayout()
        button_layout = QHBoxLayout()
        info_layout = QGridLayout()
        output_layout = QVBoxLayout()

        output_layout.addWidget(plot_tabs)
        output_layout.addLayout(info_layout)
        info_layout.addWidget(QLabel('Total Impulse:'),0,0,Qt.AlignRight)
        info_layout.addWidget(QLabel(''),1,0)
        info_layout.addWidget(QLabel('ISP:'),2,0,Qt.AlignRight)
        info_layout.addWidget(QLabel('Burnout Time:'),0,2,Qt.AlignRight)
        info_layout.addWidget(QLabel('Propellant Mass:'),2,2,Qt.AlignRight)
        total_impulse_output = QLabel('-')
        info_layout.addWidget(total_impulse_output,0,1,Qt.AlignLeft)
        ISP_output = QLabel('-')
        info_layout.addWidget(ISP_output,2,1,Qt.AlignLeft)
        burnout_time_output = QLabel('-')
        info_layout.addWidget(burnout_time_output,0,3,Qt.AlignLeft)
        propellant_mass_output = QLabel('-')
        info_layout.addWidget(propellant_mass_output,2,3,Qt.AlignLeft)
        info_layout.addWidget(QLabel('Average Thrust:'),0,4,Qt.AlignRight)
        info_layout.addWidget(QLabel('Maximum Thrust:'),2,4,Qt.AlignRight)
        info_layout.addWidget(QLabel('Average Pressure:'),0,6,Qt.AlignRight)
        info_layout.addWidget(QLabel('Maximum Pressure:'),2,6,Qt.AlignRight)
        average_thrust_output = QLabel('-')
        info_layout.addWidget(average_thrust_output,0,5,Qt.AlignLeft)
        maximum_thrust_output = QLabel('-')
        info_layout.addWidget(maximum_thrust_output,2,5,Qt.AlignLeft)
        average_pressure_output = QLabel('-')
        info_layout.addWidget(average_pressure_output,0,7,Qt.AlignLeft)
        maximum_pressure_output = QLabel('-')
        info_layout.addWidget(maximum_pressure_output,2,7,Qt.AlignLeft)
        info_layout.addWidget(QLabel('Expansion Ratio:'),0,8,Qt.AlignRight)
        info_layout.addWidget(QLabel('Port/Throat Ratio:'),2,8,Qt.AlignRight)
        expansion_ratio_output = QLabel('-')
        info_layout.addWidget(expansion_ratio_output,0,9,Qt.AlignLeft)
        port_throat_ratio_output = QLabel('-')
        info_layout.addWidget(port_throat_ratio_output,2,9,Qt.AlignLeft)


        def propellant_editor():
            propellant_window = QDialog()
            propellant_window.setWindowTitle("Propellant Editor")
            propellant_window.setFixedSize(400, 300)
            propellant_window_layout = QGridLayout(propellant_window)

            propellant_window_layout.addWidget(QLabel('Specific Heat Ratio:'),0,0,Qt.AlignRight)
            self.specific_heat_ratio_input = QLineEdit('1.1361')
            propellant_window_layout.addWidget(self.specific_heat_ratio_input,0,1,Qt.AlignLeft)

            propellant_window_layout.addWidget(QLabel('Density:'),1,0,Qt.AlignRight)
            self.density_input = QLineEdit('1700')
            propellant_window_layout.addWidget(self.density_input,1,1,Qt.AlignLeft)

            propellant_window_layout.addWidget(QLabel('Products Molecular Mass:'),2,0,Qt.AlignRight)
            self.products_molecular_mass_input = QLineEdit('39.9e-3')
            propellant_window_layout.addWidget(self.products_molecular_mass_input,2,1,Qt.AlignLeft)

            propellant_window_layout.addWidget(QLabel('Combustion Temperature(K):'),3,0,Qt.AlignRight)
            self.combustion_temperature_input = QLineEdit('1600')
            propellant_window_layout.addWidget(self.combustion_temperature_input,3,1,Qt.AlignLeft)

            propellant_window_layout.addWidget(QLabel('Burn Rate a factor:'),4,0,Qt.AlignRight)
            self.burn_rate_a_factor_input = QLineEdit('5.8')
            propellant_window_layout.addWidget(self.burn_rate_a_factor_input,4,1,Qt.AlignLeft)

            propellant_window_layout.addWidget(QLabel('Burn Rate n factor:'),5,0,Qt.AlignRight)
            self.burn_rate_n_factor_input = QLineEdit('0.22')
            propellant_window_layout.addWidget(self.burn_rate_n_factor_input,5,1,Qt.AlignLeft)

            propellant_window_layout.addWidget(QLabel('Propellant Name:'),6,0,Qt.AlignRight)
            self.propellant_name_input = QLineEdit('Projeto Jupiter KNSB')
            propellant_window_layout.addWidget(self.propellant_name_input,6,1,Qt.AlignLeft)

            save_new_propellant_button = QPushButton('Save New Propellant')
            propellant_window_layout.addWidget(save_new_propellant_button,7,0)
            save_new_propellant_button.clicked.connect(save_new_propellant)

            edit_existing_propellant_button = QPushButton('Edit Existing Propellant')
            propellant_window_layout.addWidget(edit_existing_propellant_button,7,1)
            edit_existing_propellant_button.clicked.connect(overwrite_existing_propellant)

            propellant_window.exec()

        #export json with the propellant inputs to the propellants folder with the name of the propellant
        def save_new_propellant():
            propellant_data = {
                'name': self.propellant_name_input.text(),
                'specific_heat_ratio': self.specific_heat_ratio_input.text(),
                'density': self.density_input.text(),
                'products_molecular_mass': self.products_molecular_mass_input.text(),
                'combustion_temperature': self.combustion_temperature_input.text(),
                'burn_rate_a': self.burn_rate_a_factor_input.text(),
                'burn_rate_n': self.burn_rate_n_factor_input.text()
            }

            #export the data to a json file in the propellants data folder
            propellant_data_json = json.dumps(propellant_data)
            propellant_name = self.propellant_name_input.text()
            folder = os.getcwd() + '\GUI\propellants data\ '
            propellant_file = open(folder + propellant_name + '.json', 'w')
            propellant_file.write(propellant_data_json)
            propellant_file.close()

            #warning dialog if the user ytries to save a propellant with the same name as another propellant
            if os.path.exists(folder + propellant_name + '.json'):
                warning_dialog = QMessageBox()
                warning_dialog.setWindowTitle('Warning')
                warning_dialog.setText('A propellant with this name already exists. Try to edit the existing propellant.')
                warning_dialog.exec()

        #button to edit existing propellants
        def overwrite_existing_propellant():
            propellant_data = {
                'name': self.propellant_name_input.text(),
                'specific_heat_ratio': self.specific_heat_ratio_input.text(),
                'density': self.density_input.text(),
                'products_molecular_mass': self.products_molecular_mass_input.text(),
                'combustion_temperature': self.combustion_temperature_input.text(),
                'burn_rate_a': self.burn_rate_a_factor_input.text(),
                'burn_rate_n': self.burn_rate_n_factor_input.text()
            }
           
            propellant_data_json = json.dumps(propellant_data)
            propellant_name = self.propellant_name_input.text()
            folder = os.getcwd() + '\GUI\propellants data\ '

            #warning dialog if the user ytries to overwrite a propellant that doesn't exist
            if not os.path.exists(folder + propellant_name + '.json'):
                warning_dialog = QMessageBox()
                warning_dialog.setWindowTitle('Warning')
                warning_dialog.setText('A propellant with this name does not exist. Try to create a new propellant.')
                warning_dialog.exec()
            else:
                #export the data to a json file in the propellants data folder
                propellant_file = open(folder + propellant_name + '.json', 'w')
                propellant_file.write(propellant_data_json)
                propellant_file.close()

        #refresh the propellant combobox
        def refresh_propellant_combobox():
            self.propellant_combobox.clear()
            folder = os.getcwd() + '\GUI\propellants data'
            for file in os.listdir(folder):
                if file.endswith('.json'):
                    self.propellant_combobox.addItem(file[:-5])

        def run():
            self.casing_inner_diameter_input = self.chamber_inner_diameter_input.text()
            export_button.setEnabled(True)
            #open the propellant json file and load the data
            folder = os.getcwd() + '\GUI\propellants data\ '
            propellant_file = open(folder[:-1] + self.propellant_combobox.currentText() + '.json', 'r')
            propellant_data = json.loads(propellant_file.read())
            propellant_file.close()

            KNSB = Propellant(
            specific_heat_ratio=ast.literal_eval(propellant_data['specific_heat_ratio']),
            density=ast.literal_eval(propellant_data['density']),
            products_molecular_mass=ast.literal_eval(propellant_data['products_molecular_mass']),
            combustion_temperature=ast.literal_eval(propellant_data['combustion_temperature']),
            burn_rate_a=ast.literal_eval(propellant_data['burn_rate_a']),
            burn_rate_n=ast.literal_eval(propellant_data[ 'burn_rate_n'])
            )

            Bates_Grain = Bates(
            outer_radius=ast.literal_eval(self.outer_diameter_input.text()) / 2000,
            inner_radius=ast.literal_eval(self.initial_inner_diameter_input.text()) / 2000,
            ) 

            Leviata = Motor(
            Bates_Grain,
            grain_number=ast.literal_eval(self.number_input.text()),
            chamber_inner_radius=ast.literal_eval(self.chamber_inner_diameter_input.text()) / 2000,
            nozzle_throat_radius=ast.literal_eval(self.nozzle_throat_diameter_input.text()) / 2000,
            nozzle_exit_radius=ast.literal_eval(self.nozzle_exit_diameter_input.text()) / 2000,
            nozzle_angle=ast.literal_eval(self.nozzle_angle_input.text()) * np.pi / 180,
            chamber_length=ast.literal_eval(self.chamber_lenght_input.text()) / 1000,
            )
            
            Simulation = BurnSimulation(Leviata, KNSB)
            SimulationData = BurnExport(Simulation)
            SimulationData.all_info()

            Material = StructuralCasing.Material(
            ast.literal_eval(self.ultimate_tensile_strength_input.text())*1e6,
            ast.literal_eval(self.shear_modulus_input.text())*1e6,
            ast.literal_eval(self.yield_strength_input.text())*1e6,
            ast.literal_eval(self.bearing_strength_input.text())*1e6
            )

            Geometry = StructuralCasing.Geometry(
            casing_inside_diameter = ast.literal_eval(self.casing_inner_diameter_input)/1e3,
            minor_bolt_diameter = self.minor_bolt_diameter_dict[self.bolt_types_widget.currentText()]/1e3,
            bolts_number = ast.literal_eval(self.bolts_number_input.text()),
            wall_thickness = ast.literal_eval(self.wall_thickness_input.text())/1e3,
            edge_distance = ast.literal_eval(self.edge_distance_input.text())/1e3,
            major_bolt_diameter = self.major_bolt_diameter_dict[self.bolt_types_widget.currentText()]/1e3
            )

            Casing = StructuralCasing.StructuralAnalysis(
            SimulationData.max_chamber_pressure[0],
            Material,
            Geometry
            )

            time_array = Simulation.time
            thrust_array = Simulation.thrust
            pressure_array = Simulation.chamber_pressure
            for i,pressure in enumerate(pressure_array):
                pressure_array[i] = pressure*1e-5

            plot_thrust.axes.clear()
            plot_thrust.axes.set_title('Thrust Curve')
            plot_thrust.axes.set_ylabel('Thrust(N)')
            plot_thrust.axes.set_xlabel('Time(s)')
            plot_thrust.axes.plot(time_array, thrust_array)
            plot_thrust.axes.grid()
            plot_thrust.draw()

            plot_pressure.axes.clear()
            plot_pressure.axes.set_title('Pressure Curve')
            plot_pressure.axes.set_ylabel('Pressure(bar)')
            plot_pressure.axes.set_xlabel('Time(s)')
            plot_pressure.axes.plot(time_array, pressure_array)
            plot_pressure.axes.grid()
            plot_pressure.draw()

            total_impulse_output.setText("{:.0f}Ns".format(Simulation.evaluate_total_impulse(thrust_array,time_array)))
            ISP_output.setText("{:.0f}s".format(Simulation.evaluate_specific_impulse(thrust_array,time_array)))
            burnout_time_output.setText("{:.2f}s".format(time_array[-1]))
            propellant_mass_output.setText("{:.2f}kg".format(Simulation.initial_propellant_volume*Simulation.propellant.density))
            average_thrust_output.setText("{:.0f}N".format(Export.positive_mean(thrust_array)))
            maximum_thrust_output.setText("{:.0f}N at {:.2f}s".format(*SimulationData.max_thrust))
            average_pressure_output.setText("{:.2f}bar".format(Export.positive_mean(pressure_array)))
            maximum_pressure_output.setText("{:.2f}bar at {:.2f}s".format(SimulationData.max_chamber_pressure[0]*1e-5,SimulationData.max_chamber_pressure[1]))
            expansion_ratio_output.setText("{:.2f}".format((ast.literal_eval(self.nozzle_exit_diameter_input.text())/ast.literal_eval(self.nozzle_throat_diameter_input.text()))**2))
            port_throat_ratio_output.setText("{:.2f}".format((ast.literal_eval(self.initial_inner_diameter_input.text())/ast.literal_eval(self.nozzle_throat_diameter_input.text()))**2))


            #Structural Analysis
            self.bolt_shear_output.setText("{:.2f}".format(Casing.bolt_shear()))
            self.bolt_tear_out_output.setText("{:.2f}".format(Casing.bolt_tear_out()))
            self.casing_tensile_output.setText("{:.2f}".format(Casing.casing_tensile()))
            self.bearing_output.setText("{:.2f}".format(Casing.bearing()))
            self.casing_circumferential_output.setText("{:.2f}".format(Casing.casing_circumferential_tresca()))

            #safety_factor versus bolts_number plot
            bolt_number_array = [4,8,12,16,20,24]
            safety_factor_bolt_shear_array = []
            safety_factor_bolt_tear_out_array = []
            safety_factor_casing_tensile_array = []
            safety_factor_bearing_array = []
            
            for bolt_number in bolt_number_array:
                Geometry = StructuralCasing.Geometry(
                casing_inside_diameter=ast.literal_eval(self.casing_inner_diameter_input)/1e3,
                minor_bolt_diameter=self.minor_bolt_diameter_dict[self.bolt_types_widget.currentText()]/1e3,
                bolts_number=bolt_number,
                wall_thickness=ast.literal_eval(self.wall_thickness_input.text())/1e3,
                edge_distance=ast.literal_eval(self.edge_distance_input.text())/1e3,
                major_bolt_diameter = self.major_bolt_diameter_dict[self.bolt_types_widget.currentText()]/1e3
                )
                Casing = StructuralCasing.StructuralAnalysis(
                SimulationData.max_chamber_pressure[0],
                Material,
                Geometry
                )
                safety_factor_bolt_shear_array.append(Casing.bolt_shear())
                safety_factor_bolt_tear_out_array.append(Casing.bolt_tear_out())
                safety_factor_casing_tensile_array.append(Casing.casing_tensile())
                safety_factor_bearing_array.append(Casing.bearing())

            plot_bolts_number.axes.clear()
            plot_bolts_number.axes.set_title('Safety Factor vs Bolts Number')
            plot_bolts_number.axes.set_ylabel('Safety Factor')
            plot_bolts_number.axes.set_xlabel('Bolts Number')
            plot_bolts_number.axes.set_xticks(bolt_number_array)
            plot_bolts_number.axes.plot(bolt_number_array, safety_factor_bolt_shear_array, label = 'Bolt Shear')
            plot_bolts_number.axes.plot(bolt_number_array, safety_factor_bolt_tear_out_array, label = 'Bolt Tear Out')
            plot_bolts_number.axes.plot(bolt_number_array, safety_factor_casing_tensile_array, label = 'Casing Tensile')
            plot_bolts_number.axes.plot(bolt_number_array, safety_factor_bearing_array, label = 'Bearing')
            plot_bolts_number.axes.legend()
            plot_bolts_number.axes.grid()
            plot_bolts_number.draw()

            #safety_factor versus metric bolts plot
            x_label_array = list(self.major_bolt_diameter_dict.keys())
            safety_factor_bolt_shear_array = []
            safety_factor_bolt_tear_out_array = []
            safety_factor_casing_tensile_array = []
            safety_factor_bearing_array = []
            
            for element in x_label_array:
                Geometry = StructuralCasing.Geometry(
                casing_inside_diameter=ast.literal_eval(self.casing_inner_diameter_input),
                minor_bolt_diameter=self.minor_bolt_diameter_dict[element],
                bolts_number=ast.literal_eval(self.bolts_number_input.text()),
                wall_thickness=ast.literal_eval(self.wall_thickness_input.text()),
                edge_distance=ast.literal_eval(self.edge_distance_input.text()),
                major_bolt_diameter=self.major_bolt_diameter_dict[element],
                )
                Casing = StructuralCasing.StructuralAnalysis(
                SimulationData.max_chamber_pressure[0],
                Material,
                Geometry
                )
                safety_factor_bolt_shear_array.append(Casing.bolt_shear())
                safety_factor_bolt_tear_out_array.append(Casing.bolt_tear_out())
                safety_factor_casing_tensile_array.append(Casing.casing_tensile())
                safety_factor_bearing_array.append(Casing.bearing())

            plot_bolts_type.axes.clear()
            plot_bolts_type.axes.set_title('Safety Factor vs Bolts Type')
            plot_bolts_type.axes.set_ylabel('Safety Factor')
            plot_bolts_type.axes.set_xticklabels(x_label_array,rotation=45)
            plot_bolts_type.axes.plot(x_label_array, safety_factor_bolt_shear_array, label = 'Bolt Shear')
            plot_bolts_type.axes.plot(x_label_array, safety_factor_bolt_tear_out_array, label = 'Bolt Tear Out')
            plot_bolts_type.axes.plot(x_label_array, safety_factor_casing_tensile_array, label = 'Casing Tensile')
            plot_bolts_type.axes.plot(x_label_array, safety_factor_bearing_array, label = 'Bearing')
            plot_bolts_type.axes.legend()
            plot_bolts_type.axes.grid()
            plot_bolts_type.draw()

            #safety_factor versus wall_thickness plot
            wall_thickness_array = list(range(1,10))
            safety_factor_bolt_tear_out_array = []
            safety_factor_casing_tensile_array = []
            safety_factor_bearing_array = []
            safety_factor_circumferential_tresca_array = []

            for wall_thickness in wall_thickness_array:
                Geometry = StructuralCasing.Geometry(
                    casing_inside_diameter=ast.literal_eval(self.casing_inner_diameter_input)/1e3,
                    minor_bolt_diameter=self.minor_bolt_diameter_dict[self.bolt_types_widget.currentText()]/1e3,
                    bolts_number=ast.literal_eval(self.bolts_number_input.text()),
                    wall_thickness=wall_thickness/1e3,
                    edge_distance=ast.literal_eval(self.edge_distance_input.text())/1e3,
                    major_bolt_diameter=self.major_bolt_diameter_dict[self.bolt_types_widget.currentText()]/1e3,
                )
                Casing = StructuralCasing.StructuralAnalysis(
                    SimulationData.max_chamber_pressure[0],
                    Material,
                    Geometry
                )
                safety_factor_bolt_tear_out_array.append(Casing.bolt_tear_out())
                safety_factor_casing_tensile_array.append(Casing.casing_tensile())
                safety_factor_bearing_array.append(Casing.bearing())
                safety_factor_circumferential_tresca_array.append(Casing.casing_circumferential_tresca())

            plot_wall_thickness.axes.clear()
            plot_wall_thickness.axes.set_title('Safety Factor vs Wall Thickness')
            plot_wall_thickness.axes.set_ylabel('Safety Factor')
            plot_wall_thickness.axes.set_xlabel('Wall Thickness [mm]')
            plot_wall_thickness.axes.set_xticks(wall_thickness_array)
            plot_wall_thickness.axes.plot(wall_thickness_array, safety_factor_bolt_tear_out_array, label = 'Bolt Tear Out')
            plot_wall_thickness.axes.plot(wall_thickness_array, safety_factor_casing_tensile_array, label = 'Casing Tensile')
            plot_wall_thickness.axes.plot(wall_thickness_array, safety_factor_bearing_array, label = 'Bearing')
            plot_wall_thickness.axes.plot(wall_thickness_array, safety_factor_circumferential_tresca_array, label = 'Casing Circumferential Tresca')
            plot_wall_thickness.axes.legend()
            plot_wall_thickness.axes.grid()
            plot_wall_thickness.draw()

            return Simulation

        def export_thrust_curve():
            #Export Thrust Curve
            Simulation = run()
            filename = QFileDialog.getSaveFileName(self, 'Export Thrust Curve', 'Thrust Curve', 'CSV(*.csv)')
            if filename[0]:
                with open(filename[0], 'w', newline='') as file:
                    writer = csv.writer(file)
                    writer.writerow(['Time [s]', 'Thrust [N]'])
                    for i in range(len(Simulation.time)):
                        writer.writerow([Simulation.time[i], Simulation.thrust[i]])



        input_layout.addWidget(information_tabs)
        propellant_editor_button = QPushButton('Propellant Editor')
        propellant_editor_button.clicked.connect(propellant_editor)
        Grain_Layout.addWidget(propellant_editor_button,5,0)
        refresh_propellant_combobox_button = QPushButton('Refresh Propellants')
        refresh_propellant_combobox_button.clicked.connect(refresh_propellant_combobox)
        Grain_Layout.addWidget(refresh_propellant_combobox_button,5,1)
        #organize the grain layout widgets distances




        run_button = QPushButton('Run Simulation')
        run_button.clicked.connect(run)
        button_layout.addWidget(run_button)
        # button_layout.addWidget(QPushButton('Save Simulation'))
        export_button = QPushButton('Export Thrust Curve')
        export_button.setEnabled(False)
        export_button.clicked.connect(export_thrust_curve)
        button_layout.addWidget(export_button)
        input_layout.addLayout(button_layout)
        pagelayout.addLayout(input_layout)
        pagelayout.addLayout(output_layout)
        output_layout.addStretch()



        widget = QWidget()
        widget.setLayout(pagelayout)
        self.setCentralWidget(widget)


app = QApplication(sys.argv)

window = MainWindow()
window.showMaximized()

app.exec()