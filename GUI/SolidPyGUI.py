import sys, ast
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
    Export
)

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

        
        plot_thrust = MplCanvas(self, width=15, height=8)


        plot_pressure = MplCanvas(self, width=15, height=8)


        plot_kn = MplCanvas(self, width=15, height=8)
        plot_kn.axes.set_title('Kn Curve')
        plot_kn.axes.set_ylabel('Kn')
        plot_kn.axes.set_xlabel('Time(s)')

        plot_mass = MplCanvas(self, width=15, height=8)
        plot_mass.axes.set_title('Mass Curve')
        plot_mass.axes.set_ylabel('Propellant Mass(kg)')
        plot_mass.axes.set_xlabel('Time(s)')

        information_tabs = QTabWidget()

        Propellant_Widget = QWidget()
        Grain_Widget = QWidget()
        Motor_Widget = QWidget()
        Structure_Widget = QWidget()

        Propellant_Layout = QGridLayout(Propellant_Widget)

        Propellant_Layout.addWidget(QLabel('Specific Heat Ratio:'),0,0,Qt.AlignRight)
        self.specific_heat_ratio_input = QLineEdit('1.1361')
        Propellant_Layout.addWidget(self.specific_heat_ratio_input,0,1,Qt.AlignLeft)

        Propellant_Layout.addWidget(QLabel('Density:'),1,0,Qt.AlignRight)
        self.density_input = QLineEdit('1700')
        Propellant_Layout.addWidget(self.density_input,1,1,Qt.AlignLeft)

        Propellant_Layout.addWidget(QLabel('Products Molecular Mass:'),2,0,Qt.AlignRight)
        self.products_molecular_mass_input = QLineEdit('39.9e-3')
        Propellant_Layout.addWidget(self.products_molecular_mass_input,2,1,Qt.AlignLeft)

        Propellant_Layout.addWidget(QLabel('Combustion Temperature(K):'),3,0,Qt.AlignRight)
        self.combustion_temperature_input = QLineEdit('1600')
        Propellant_Layout.addWidget(self.combustion_temperature_input,3,1,Qt.AlignLeft)

        Propellant_Layout.addWidget(QLabel('Burn Rate a factor:'),4,0,Qt.AlignRight)
        self.burn_rate_a_factor_input = QLineEdit('5.8')
        Propellant_Layout.addWidget(self.burn_rate_a_factor_input,4,1,Qt.AlignLeft)

        Propellant_Layout.addWidget(QLabel('Burn Rate n factor:'),5,0,Qt.AlignRight)
        self.burn_rate_n_factor_input = QLineEdit('0.22')
        Propellant_Layout.addWidget(self.burn_rate_n_factor_input,5,1,Qt.AlignLeft)

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

        

        Structure_Layout = QVBoxLayout(Structure_Widget)

        information_tabs.addTab(Propellant_Widget, 'Propellant')
        information_tabs.addTab(Grain_Widget, 'Grain')
        information_tabs.addTab(Motor_Widget, 'Motor')
        information_tabs.addTab(Structure_Widget, 'Structure')

        plot_tabs = QTabWidget()
        plot_tabs.addTab(plot_thrust,'Thrust')
        plot_tabs.addTab(plot_pressure,'Pressure')
        plot_tabs.addTab(plot_kn,'Kn')
        plot_tabs.addTab(plot_mass,'Mass')

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

        def run():

            KNSB = Propellant(
            specific_heat_ratio=ast.literal_eval(self.specific_heat_ratio_input.text()),
            density=ast.literal_eval(self.density_input.text()),
            products_molecular_mass=ast.literal_eval(self.products_molecular_mass_input.text()),
            combustion_temperature=ast.literal_eval(self.combustion_temperature_input.text()),
            burn_rate_a=ast.literal_eval(self.burn_rate_a_factor_input.text()),
            burn_rate_n=ast.literal_eval(self.burn_rate_n_factor_input.text())
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



        input_layout.addWidget(information_tabs)
        a = QPushButton('Run Simulation')
        a.clicked.connect(run)
        button_layout.addWidget(a)
        button_layout.addWidget(QPushButton('Save Simulation'))
        button_layout.addWidget(QPushButton('Export Thrust Curve'))
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