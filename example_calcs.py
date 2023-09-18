from aisi_design_functions import *
import aisi_design_functions

# Design examples using aisi_design_functions module:

# Axial Examples -------------------------------------------------------------------------------------------------------------------------------
axialResultSingle = getAxialStrength_Single('600S200-97','stud',10,50,4,29500,11300,1.5,True,1,1,1,10,4,4)
print(f'P_allow single: {axialResultSingle}')

axialResultBoxed = getAxialStrength_Boxed('600S200-97', 'stud', 11,50,4,29500,11300,1.5,True,1,1,11,4,4,12)
print(f'P_allow Boxed: {axialResultBoxed}')

axialResultB2B = getAxialStrength_B2B('362S162-43','stud',7,50,7,29500,11300,1.5,True,1,1,1,7,7,7,12)
print(f'P_allow Back to Back: {axialResultB2B}')


# Flexural Examples -------------------------------------------------------------------------------------------------------------------------------
flexuralSingleResult = getFlexuralStrength_Single('600T150-43', 'track', 15,50,0.001,29500,11300,4,1.5,False,1,1,15,15,15,1.14,False,15)
print(f'M_allow Single Track: {flexuralSingleResult}')

flexuralBoxedResult = getFlexuralStrength_Boxed('362S162-43','stud',18,33,29500,11300,4.0,1.5,True,1,1,1.13,12)
print(f'M_allow Boxed: {flexuralBoxedResult}')

flexuralBoxedResultTrack = getFlexuralStrength_Boxed('600T150-43','track',15,50,29500,11300,4.0,1.5,False,1,1,1.13,12)
print(f'M_allow Boxed TRACKS: {flexuralBoxedResultTrack}')

flexuralB2BResult = getFlexuralStrength_B2B('362S162-43','stud',15,50,15,29500,11300,4.0,1.5,True,1,1,15,15,15,1.00,False,15,12)
print(f'M_allow Back to Back: {flexuralB2BResult}')

flexuralB2BResultTrack = getFlexuralStrength_B2B('600T150-97','track',15,50,15,29500,11300,4.0,1.5,False,1,1,15,15,15,1.00,False,15,12)
print(f'M_allow Back to Back Track Section: {flexuralB2BResultTrack}')

# Globally Braced Flexural Examples ---------------------------------------------------------------------------------

globallyBracedSingle = getFlexuralStrength_Single_GB('600S200-97','stud',10,10,10,10,10,50,29500,11300,True, 1.5,4.5,1,1,1.15)
print(f'Globally Braced M_allow single: {globallyBracedSingle}')

globallyBracedSingleTrack = getFlexuralStrength_Single_GB('600T150-97','track',10,10,10,10,10,50,29500,11300,False, 1.5,4.5,1,1,1.15)
print(f'Globally Braced M_allow single TRACK: {globallyBracedSingleTrack}')

globallyBracedBoxed = getFlexuralStrength_Boxed_GB('600S200-97', 'stud',10,50,29500,11300,True,1.5,4.0,1,1,1)
print(f'Globally Braced M_allow Boxed: {globallyBracedBoxed}')

globallyBracedBoxedTrack = getFlexuralStrength_Boxed_GB('600T150-43', 'track',10,50,29500,11300,False,1.5,4.0,1,1,1)
print(f'Globally Braced M_allow Boxed TRACK: {globallyBracedBoxedTrack}')

globallyBracedBackToBack = getFlexuralStrength_B2B_GB('600S200-97','stud', 10,10,50,29500,11300,False,1.5,4.0,10,10,10,1,1,1)
print(f'Globally Braced M_allow back to back: {globallyBracedBackToBack}')

globallyBracedBackToBackTrack = getFlexuralStrength_B2B_GB('600T200-97','track', 10,10,50,29500,11300,False,1.5,4.0,10,10,10,1,1,1)
print(f'Globally Braced M_allow back to back TRACK: {globallyBracedBackToBackTrack}')


# Shear Examples ---------------------------------------------------------------------------------
shearSingleResult = getShearStrength_Single('600S162-54','stud',50,29500,11300,1.5,True)
print(f'V_allow single: {shearSingleResult}')

shearSingleResultTrack = getShearStrength_Single('600T150-54','track',50,29500,11300,1.5,False)
print(f'V_allow single TRACK: {shearSingleResultTrack}')

shearBuiltUpResult = getShearStrength_BuiltUp('362S162-68','stud',50,29500,11300,1.5,False)
print(f'V_allow built up: {shearBuiltUpResult}')

shearBuiltUpResultTrack = getShearStrength_BuiltUp('600T200-97','track',50,29500,11300,1.5,False)
print(f'V_allow built up TRACK: {shearBuiltUpResultTrack}')

# HSS section functions: Examples ---------------------------------------------------------------------------------
axialResultHSS = getHSSAxialStrength('HSS8X4X1/8',46,15,15)
print(f'P_allow HSS: {axialResultHSS}')

HSSFlexural = getHSSFlexuralStrength('HSS8X6X3/16',46,15,'strong')
print(f'HSS M_allow: {HSSFlexural}')

HSSShear = getHSSShearStrength('HSS6X4X1/2', 46, 'weak')
print(f'HSS V_allow: {HSSShear}')




