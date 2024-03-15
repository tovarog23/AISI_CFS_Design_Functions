from aisi_design_functions import *

# Axial Examples -------------------------------------------------------------------------------------------------------------------------------
axialResultSingle = getAxialStrength_Single('600S200-97','stud',50,10,29500,11300,1.5,True,1,1,1,10,10,True,4)
print(f'P_allow single: {axialResultSingle}')

axialResultBoxed = getAxialStrength_Boxed('600S200-97','stud',50,4,29500,11300,1.5,True,1,1,11,12)
print(f'P_allow Boxed: {axialResultBoxed}')

axialResultB2B = getAxialStrength_B2B('362S162-43','stud',50,7,29500,11300,1.5,True,1,1,1,7,7,12)
print(f'P_allow Back to Back: {axialResultB2B}')

# Flexural Examples -------------------------------------------------------------------------------------------------------------------------------
flexuralSingleResult = getFlexuralStrength_Single('600T150-43','track',50,10,29500,11300,4,1.5,False,1,1,10,1.14,False)
print(f'M_allow Single Track: {flexuralSingleResult}')

flexuralBoxedResult = getFlexuralStrength_Boxed('362S162-43','stud',14,50,29500,11300,4,1.5,True,1,1.14,12)
print(f'M_allow Boxed: {flexuralBoxedResult}')

flexuralB2BResult = getFlexuralStrength_B2B('362S162-43','stud',11,50,4,29500,11300,4,1.5,False,1,1,11,1.12,12,True,4)
print(f'M_allow Back to Back: {flexuralB2BResult}')

# Globally Braced Flexural Examples ---------------------------------------------------------------------------------

globallyBracedSingle = getFlexuralStrength_Single_GB('600S200-97','stud',50,10,29500,11300,4,1.5,True,1,1,10,1.12)
print(f'Globally Braced M_allow single: {globallyBracedSingle}')

globallyBracedSingleTrack = getFlexuralStrength_Single_GB('600T150-97','track',50,12,29500,11300,4,1.5,False,1,1,12,1.1)
print(f'Globally Braced M_allow single TRACK: {globallyBracedSingleTrack}')

globallyBracedBoxed = getFlexuralStrength_Boxed_GB('600S200-97','stud',10,50,29500,11300,4,1.5,True,1,1.14,12)
print(f'Globally Braced M_allow Boxed: {globallyBracedBoxed}')

globallyBracedBackToBack = getFlexuralStrength_B2B_GB('600S200-97','stud',10,8,50,29500,11300,4,1.5,False,1,1,10,1.14)
print(f'Globally Braced M_allow back to back: {globallyBracedBackToBack}')

# Shear Examples ---------------------------------------------------------------------------------
shearSingleResult = getShearStrength_Single('600S162-54','stud',50,29500,11300,1.5,True)
print(f'V_allow single: {shearSingleResult}')

shearSingleResultTrack = getShearStrength_Single('600T150-54','track',50,29500,11300,1.5,True)
print(f'V_allow single TRACK: {shearSingleResultTrack}')

shearBuiltUpResult = getShearStrength_BuiltUp('362S162-68','stud',50,29500,11300,1.5,False)
print(f'V_allow built up: {shearBuiltUpResult}')

shearBuiltUpResultTrack = getShearStrength_BuiltUp('600T200-97','track',50,29500,11300,1.5,True)
print(f'V_allow built up TRACK: {shearBuiltUpResultTrack}')

# HSS section functions: Examples ---------------------------------------------------------------------------------
axialResultHSS = getHSSAxialStrength('HSS8X4X1/8',46,10,10)
print(f'P_allow HSS: {axialResultHSS}')

HSSFlexural = getHSSFlexuralStrength('HSS8X6X3/16',46,12,'strong')
print(f'HSS M_allow: {HSSFlexural}')

HSSShear = getHSSShearStrength('HSS6X4X1/2',46,'weak')
print(f'HSS V_allow: {HSSShear}')