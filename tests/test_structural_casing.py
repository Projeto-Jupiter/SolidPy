from solidpy import StructuralCasing
import numpy as np
import matplotlib.pyplot as plt

#material class inputs 
UTS = 12e4
shear_strength = 3e4
YTS = 38e3
BYS =56e3

#geometry class inputs
casing_inside_diameter = 3.624
minor_bolt_diameter = 0.2052
bolts_number = 12
wall_thickness = 0.188
edge_distance = 0.4375
casing_outer_diameter = 4
major_bolt_diameter = 0.25

#MEOP
MEOP = 1400

#defining classes
material = StructuralCasing.Material(UTS,shear_strength,YTS,BYS)
geometry = StructuralCasing.Geometry(casing_inside_diameter,minor_bolt_diameter,bolts_number,wall_thickness,edge_distance,casing_outer_diameter,major_bolt_diameter)
structural_analysis = StructuralCasing.StructuralAnalysis(MEOP,material,geometry)

#testing
print(structural_analysis.bolt_shear())
print(structural_analysis.bolt_tear_out())
print(structural_analysis.casing_tensile())
print(structural_analysis.bearing())
print(structural_analysis.casing_tresca())

#safety_factor versus bolts_number plot
bolts_number_array = np.linspace(4,24,6)
safety_factor_bolt_shear_array = []
safety_factor_bolt_tear_out_array = []
safety_factor_casing_tensile_array = []
safety_factor_bearing_array = []

for i in bolts_number_array:
    geometry = StructuralCasing.Geometry(casing_inside_diameter,minor_bolt_diameter,i,wall_thickness,edge_distance,casing_outer_diameter,major_bolt_diameter)
    structural_analysis = StructuralCasing.StructuralAnalysis(MEOP,material,geometry)
    safety_factor_bolt_shear_array.append(structural_analysis.bolt_shear())
    safety_factor_bolt_tear_out_array.append(structural_analysis.bolt_tear_out())
    safety_factor_casing_tensile_array.append(structural_analysis.casing_tensile())
    safety_factor_bearing_array.append(structural_analysis.bearing())

plt.title('Safety Factor vs Bolts Number')
plt.xlabel('Bolts Number')
plt.xticks(bolts_number_array)
plt.ylabel('Safety Factor')
plt.plot(bolts_number_array,safety_factor_bolt_shear_array,label='bolt shear')
plt.plot(bolts_number_array,safety_factor_bolt_tear_out_array,label='bolt tear out')
plt.plot(bolts_number_array,safety_factor_casing_tensile_array,label='casing tensile')
plt.plot(bolts_number_array,safety_factor_bearing_array,label='bearing')
plt.legend()
plt.show()

#safety_factor versus metric bolts plot
metric_bolts_minor_diameter_dict = {
    'M5x0.5':4.459,
    'M5x0.8':4.134,
    'M5.5x0.5':4.959,
    'M6x0.75':5.188,
    'M6x1':4.917,
    'M7x0.75':6.188,
    'M7x1':5.917,
    'M8x0.75':7.188,
    'M8x1':6.917,
    'M8x1.25':6.647,
}

metric_bolts_major_diameter_dict = {
    'M5x0.5':5,
    'M5x0.8':5,
    'M5.5x0.5':5.5,
    'M6x0.75':6,
    'M6x1':6,
    'M7x0.75':7,
    'M7x1':7,
    'M8x0.75':8,
    'M8x1':8,
    'M8x1.25':8,
}

metric_bolts_minor_diameter_array = []
metric_bolts_major_diameter_array = []
safety_factor_bolt_shear_array = []
safety_factor_bolt_tear_out_array = []
safety_factor_casing_tensile_array = []
safety_factor_bearing_array = []
x_label_array = list(metric_bolts_major_diameter_dict.keys())
print(x_label_array)

for element in x_label_array:
    metric_bolts_minor_diameter_array.append(metric_bolts_minor_diameter_dict[element])
    metric_bolts_major_diameter_array.append(metric_bolts_major_diameter_dict[element])
    geometry = StructuralCasing.Geometry(casing_inside_diameter,metric_bolts_minor_diameter_dict[element],bolts_number,wall_thickness,edge_distance,casing_outer_diameter,metric_bolts_major_diameter_dict[element])
    structural_analysis = StructuralCasing.StructuralAnalysis(MEOP,material,geometry)
    safety_factor_bolt_shear_array.append(structural_analysis.bolt_shear())
    safety_factor_bolt_tear_out_array.append(structural_analysis.bolt_tear_out())
    safety_factor_casing_tensile_array.append(structural_analysis.casing_tensile())
    safety_factor_bearing_array.append(structural_analysis.bearing())

plt.title('Safety Factor vs Metric Bolts')
plt.xlabel('Metric Bolts')
plt.xticks(rotation=45)
plt.ylabel('Safety Factor')
plt.plot(x_label_array,safety_factor_bolt_shear_array,label='bolt shear')
plt.plot(x_label_array,safety_factor_bolt_tear_out_array,label='bolt tear out')
plt.plot(x_label_array,safety_factor_casing_tensile_array,label='casing tensile')
plt.plot(x_label_array,safety_factor_bearing_array,label='bearing')
plt.legend()
plt.show()

#safety_factor versus wall thickness plot
wall_thickness_array = np.linspace(0.039,0.39,100)
safety_factor_bolt_shear_array = []
safety_factor_bolt_tear_out_array = []
safety_factor_casing_tensile_array = []
safety_factor_bearing_array = []
safety_factor_casing_tresca_array = []

for i in wall_thickness_array:
    geometry = StructuralCasing.Geometry(casing_inside_diameter,minor_bolt_diameter,bolts_number,i,edge_distance,casing_outer_diameter,major_bolt_diameter)
    structural_analysis = StructuralCasing.StructuralAnalysis(MEOP,material,geometry)
    safety_factor_bolt_shear_array.append(structural_analysis.bolt_shear())
    safety_factor_bolt_tear_out_array.append(structural_analysis.bolt_tear_out())
    safety_factor_casing_tensile_array.append(structural_analysis.casing_tensile())
    safety_factor_bearing_array.append(structural_analysis.bearing())
    safety_factor_casing_tresca_array.append(structural_analysis.casing_tresca())

plt.title('Safety Factor vs Wall Thickness')
plt.xlabel('Wall Thickness')
plt.ylabel('Safety Factor')
plt.plot(wall_thickness_array,safety_factor_bolt_shear_array,label='bolt shear')
plt.plot(wall_thickness_array,safety_factor_bolt_tear_out_array,label='bolt tear out')
plt.plot(wall_thickness_array,safety_factor_casing_tensile_array,label='casing tensile')
plt.plot(wall_thickness_array,safety_factor_bearing_array,label='bearing')
plt.plot(wall_thickness_array,safety_factor_casing_tresca_array,label='casing tresca')
plt.legend()
plt.show()







