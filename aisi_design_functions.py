from math import pi, sqrt, sin, cos
import forallpeople
import numpy as np
import csv
import os
from scipy import integrate
import time

''' Python version used: 3.9.13
    
    Codes Used: 
        AISI S100-16
        AISI Cold Formed Steel Design Manual 2017
        AISC 360-10 (14th Edition Manual)
    
    The following script was developed with the objective of providing useful design functions for engineers and professionals involved in the analysis and design of cold formed steel structures, particularly buildings consisting of load bearing wall panels and joist floor systems.
    
    The functions themselves can be imported into other python scripts or be integrated directly to an excel spreadsheet if necessary through the use of xlwings, which allows any user to import user defined functions (UDF).
    
    CSV files were included with adequate shape information (from AISI cold fomed steel manual) for testing and verification of each function under different design conditions.   
    
    The only exception in this script is the addition of functions to calculate the flexural, axial and shear capacity of HSS hot rolled steel sections. The design of this type of shape was included since it is very common to use this type of members as posts, jambs, or headers integrated into typical load bearing cold formed steel wall panels.   
'''

forallpeople.environment('structural', top_level=True)

# GENERAL CALCS:-------------------------------------------------------------------------------------------------

def calcEffectiveSectionModulus(Depth, F_n, procedure_alt, section_type, base, R, t, poiss, E, r, D, d_h, P_width, sectionType):
    centroid = Depth/2
    num_iterations = 10;
    Max_stress = F_n
    results_list_centroids = []
    if procedure_alt == False:
        if section_type == 'track':
            for i in range(num_iterations):
                # COMPRESSION FLANGE
                w = base - (R + t)
                f = Max_stress
                k_flange = 0.43
                F_crl = k_flange*(((pi**2)*E)/(12*(1-poiss**2)))*((t/w)**2)
                lamb_1 = sqrt(f/F_crl)
                rho = (1-(0.22/lamb_1))/(lamb_1) # Eq.1.1-2
                
                if lamb_1 <= 0.673: b = w
                else: b = rho*w 

                # WEB REVISION
                w_web = Depth - 2*(R+t/2)-t
                f_web = Max_stress
                f1_web = (f_web*((centroid)-(t+R)))/(centroid)
                f2_web = -1*((f_web*(Depth - (t+R) - centroid))/centroid)
                psi_web = abs(f2_web/f1_web)

                k_web = 4+2*((1+psi_web)**3)+2*(1+psi_web)
                F_crl_web = k_web*(((pi**2)*E)/(12*(1-poiss**2)))*((t/w_web)**2)
                lamb_web = sqrt(f1_web/F_crl_web)
                rho_web = (1-(0.22/lamb_web))/(lamb_web) # Eq.1.1-2

                if lamb_web <= 0.673: 
                    b_web = w_web
                else: 
                    b_web = rho_web*w_web 
                ho = Depth
                bo = base
                if ho/bo <= 4:
                    b1 = b_web/(3+psi_web)
                    if psi_web > 0.236:
                        b2 = b_web/2
                    else: 
                        b2 = b_web - b1
                else:
                    b1 = b_web/(3+psi_web)
                    b2 = b_web/(1+psi_web)-b1
                # verify if ineffective portion exists for property calculation process
                sum = b1 + b2
                limit = centroid - (t + R)
                if sum <= limit:
                    diff = limit - sum
                else:
                    diff = 0*inch

                # Recompute properties by parts:
                b_neg = -1*(diff) # negative length of negative element, if it exists
                y_b_neg = t/2 + r + b1 - (b_neg/2) # centroidal location of negative component from top fiber

                # generation of table values according to reference example in AISI manual
                '''LENGTHS'''
                top_flange_length = b
                bottom_flange_length = w
                web_length = w_web
                web_neg_length = b_neg

                top_inside_corner = (2*pi*(R + t/2))/4
                bottom_inside_corner = (2*pi*(R + t/2))/4

                Length_vector = np.array([top_flange_length, bottom_flange_length, web_length, web_neg_length, top_inside_corner, bottom_inside_corner])
                
                Length_sum = Length_vector.sum()

                '''Y DISTANCES'''
                top_flange_ytf = t/2
                bottom_flange_ytf = Depth - t/2
                web_ytf = Depth/2
                web_neg_ytf = y_b_neg

                top_inside_corner_ytf = (t + R) - ((2*(R + t/2))/pi); # print(top_inside_corner_ytf)
                bottom_inside_corner_ytf = Depth - ((R + t) - ((2*(R + t/2))/pi)); # print(bottom_inside_corner_ytf)

                ytf_vector = np.array([top_flange_ytf, bottom_flange_ytf, web_ytf, web_neg_ytf, top_inside_corner_ytf, bottom_inside_corner_ytf])
                
                Ly_results = Length_vector*ytf_vector
                
                Ly_total = Ly_results.sum()

                Ly_2_results = Length_vector*(ytf_vector**2)
                Ly_2_total = Ly_2_results.sum()

                '''I'x CALCULATION'''

                top_flange_Ipx = (t**3)/12
                bottom_flange_Ipx = (t**3)/12
                web_Ipx = (web_length**3)/12
                web_neg_Ipx = (web_neg_length**3)/12

                top_inside_corner_Ipx = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)
                bottom_inside_corner_Ipx = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)

                Ipx_vector = np.array([top_flange_Ipx, bottom_flange_Ipx, web_Ipx, web_neg_Ipx, top_inside_corner_Ipx, bottom_inside_corner_Ipx])

                Ipx_total = Ipx_vector.sum()

                y_c = Ly_total/Length_sum

                I_xe = (Ipx_total+Ly_2_total-((y_c**2)*Length_sum))*t
                
                S_e = I_xe/y_c
                
                results_list_centroids.append(y_c)
                
                if y_c > centroid:
                    # centroid is below centerline
                    centroid = y_c
                    
                elif y_c < centroid:
                    # centroid is above centerline (max stress at bottom flange)
                    centroid = y_c
                    ref_stress = (Max_stress*centroid)/(Depth -  centroid)
                    Max_stress = ref_stress
                    
        elif section_type == 'stud':
            for i in range(num_iterations):
                # COMPRESSION FLANGE
                w = base - 2*(R + t)
                f = Max_stress
                S = 1.28*sqrt(E/f)

                if w/t <= 0.328*S:
                    d_small = (D - R - t)
                    R_I = 1
                    b = w
                else: 
                    # calculate effective width of flange
                    I_a = min(399*(t**4)*(((w/t)/S)-0.328)**3,(t**4)*(115*((w/t)/S)+5))
                    d_small = (D - R - t)
                    I_s = ((d_small**3)*t*(sin(0.5*pi))**2)/12
                    R_I = min(1, I_s/I_a)
                    n = max(1/3,0.582-((w/t)/(4*S)))
                    
                    if D/w > 0.25 and D/w <= 0.8: k = min(4,(4.82-((5*D)/w))*(R_I**n)+0.43)
                    elif D/w > 0.8: k = 1.25
                    else: k = min(3.75*(R_I**n) + 0.43,4)
                    
                    F_crl = k*(((pi**2)*E)/(12*(1-poiss**2)))*((t/w)**2)
                    lamb_1 = sqrt(f/F_crl)
                    rho = (1-(0.22/lamb_1))/(lamb_1) # Eq.1.1-2
                    
                    if lamb_1 <= 0.673: b = w
                    else: b = rho*w 
                    
                # STIFFENER LIP REVISION
                '''Max stress in lip'''
                f1_lip = f*((centroid)-(t + R))/(centroid)
                f2_lip = f*((centroid)-D)/(centroid)
                f_lip = f1_lip
                psi_lip = abs(f2_lip/f1_lip)
                k_lip = 0.578/(psi_lip+0.34)
                F_crl_lip = k_lip*(((pi**2)*E)/(12*(1-poiss**2)))*((t/d_small)**2)

                lamb_lip = sqrt(f_lip/F_crl_lip)
                rho_lip = (1-(0.22/lamb_lip))/(lamb_lip) # Eq.1.1-2

                if lamb_lip <= 0.673: 
                    # b_lip = w
                    d_p_s = 1*d_small
                    if w/t <= 0.328*S:
                        d_s = d_p_s 
                    else:
                        d_s = d_p_s*R_I
                else: 
                    # b = rho_lip*w 
                    d_p_s = rho_lip*d_small
                    if w/t <= 0.328*S:
                        d_s = d_p_s 
                    else:
                        d_s = d_p_s*R_I

                # WEB REVISION
                w_web = Depth - 2*(t + R)
                f_web = Max_stress
                f1_web = (f_web*((centroid)-(t+R)))/(centroid)
                f2_web = -1*((f_web*(Depth - (t+R) - centroid))/centroid)
                psi_web = abs(f2_web/f1_web)

                k_web = 4+2*((1+psi_web)**3)+2*(1+psi_web)
                F_crl_web = k_web*(((pi**2)*E)/(12*(1-poiss**2)))*((t/w_web)**2)
                lamb_web = sqrt(f1_web/F_crl_web)
                rho_web = (1-(0.22/lamb_web))/(lamb_web) # Eq.1.1-2

                if lamb_web <= 0.673: 
                    b_web = w_web
                else: 
                    b_web = rho_web*w_web 

                ho = Depth
                bo = base
                if ho/bo <= 4:
                    b1 = b_web/(3+psi_web)
                    if psi_web > 0.236:
                        b2 = b_web/2
                    else: 
                        b2 = b_web - b1
                else:
                    b1 = b_web/(3+psi_web)
                    b2 = b_web/(1+psi_web)-b1
                # verify if ineffective portion exists for property calculation process
                sum = b1 + b2
                limit = centroid - (t + R)
                if sum <= limit:
                    diff = limit - sum
                else:
                    diff = 0*inch

                # Recompute properties by parts:

                b_neg = -1*(diff) # negative length of negative element, if it exists
                y_b_neg = t/2 + r + b1 - (b_neg/2)# centroidal location of negative component from top fiber

                # generation of table values according to reference example in AISI manual
                '''LENGTHS'''
                top_flange_length = b
                bottom_flange_length = w
                web_length = w_web
                web_neg_length = b_neg

                top_inside_corner = (2*pi*(R + t/2))/4
                bottom_inside_corner = (2*pi*(R + t/2))/4

                top_outside_corner = (2*pi*(R + t/2))/4
                bottom_outside_corner = (2*pi*(R + t/2))/4

                top_lip = d_s
                bottom_lip = d_small

                Length_vector = np.array([top_flange_length, bottom_flange_length, web_length, web_neg_length, top_inside_corner, bottom_inside_corner, top_outside_corner, 
                                bottom_outside_corner, top_lip, bottom_lip])
                Length_sum = Length_vector.sum()

                '''Y DISTANCES'''
                top_flange_ytf = t/2
                bottom_flange_ytf = Depth - t/2
                web_ytf = Depth/2
                web_neg_ytf = y_b_neg

                top_inside_corner_ytf = (t + R) - ((2*(R + t/2))/pi); # print(top_inside_corner_ytf)
                bottom_inside_corner_ytf = Depth - ((R + t) - ((2*(R + t/2))/pi)); # print(bottom_inside_corner_ytf)

                top_outside_corner_ytf = (t + R) - ((2*(R + t/2))/pi); # print(top_inside_corner_ytf)
                bottom_outside_corner_ytf = Depth - ((R + t) - ((2*(R + t/2))/pi)); # print(bottom_inside_corner_ytf)

                top_lip_ytf = (R+t)+0.5*(d_s)
                bottom_lip_ytf = Depth - ((R+t)+0.5*(d_small))

                ytf_vector = np.array([top_flange_ytf, bottom_flange_ytf, web_ytf, web_neg_ytf, top_inside_corner_ytf, bottom_inside_corner_ytf, top_outside_corner_ytf,
                                    bottom_outside_corner_ytf, top_lip_ytf, bottom_lip_ytf])

                Ly_results = Length_vector*ytf_vector
                Ly_total = Ly_results.sum()

                Ly_2_results = Length_vector*(ytf_vector**2)
                Ly_2_total = Ly_2_results.sum()

                '''I'x CALCULATION'''

                top_flange_Ipx = (t**3)/12
                bottom_flange_Ipx = (t**3)/12
                web_Ipx = (web_length**3)/12
                web_neg_Ipx = (web_neg_length**3)/12

                top_inside_corner_Ipx = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)
                bottom_inside_corner_Ipx = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)

                top_outside_corner_Ipx = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)
                bottom_outside_corner_Ipx = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)

                top_lip_Ipx = (top_lip**3)/12
                bottom_lip_Ipx = (bottom_lip**3)/12

                Ipx_vector = np.array([top_flange_Ipx, bottom_flange_Ipx, web_Ipx, web_neg_Ipx, top_inside_corner_Ipx, bottom_inside_corner_Ipx, top_outside_corner_Ipx,
                            bottom_outside_corner_Ipx, top_lip_Ipx, bottom_lip_Ipx])

                Ipx_total = Ipx_vector.sum()
                y_c = Ly_total/Length_sum

                I_xe = (Ipx_total+Ly_2_total-((y_c**2)*Length_sum))*t
                
                S_e = I_xe/y_c
                
                results_list_centroids.append(y_c)
                
                if y_c > centroid:
                    # centroid is below centerline
                    centroid = y_c
                elif y_c < centroid:
                    # centroid is above centerline (max stress at bottom flange)
                    centroid = y_c
                    ref_stress = (Max_stress*centroid)/(Depth -  centroid)
                    Max_stress = ref_stress
        
        else:
            return 'ERROR: Incorrect section type value' 
    elif procedure_alt == True:
        if section_type == 'track':
            for i in range(num_iterations):
                # COMPRESSION FLANGE 
                w = base - (R + t)
                f = Max_stress
                k_flange = 0.43
                F_crl = k_flange*(((pi**2)*E)/(12*(1-poiss**2)))*((t/w)**2)
                lamb_1 = sqrt(f/F_crl)
                rho = (1-(0.22/lamb_1))/(lamb_1) # Eq.1.1-2
                
                if lamb_1 <= 0.673: b = w
                else: b = rho*w

                # WEB REVISION - CONSIDERED AS A UNIFORMLY COMPRESSED UNSTIFFENED ELEMENT - PROCESS ACCOUNTS TOP COMP WEB
                ''' ALTERNATIVE PROCESS CONSIDERING HOLES'''
                w_web = (Depth/2) - (d_h/2) - (R + t) 
                k_web = 0.43
                f_web = Max_stress
                f1_web = (f_web*((centroid)-(t+R)))/(centroid)
                F_crl_web = k_web*(((pi**2)*E)/(12*(1-poiss**2)))*((t/w_web)**2)
                lamb_web = sqrt(f1_web/F_crl_web)
                rho_web = (1-(0.22/lamb_web))/(lamb_web) # Eq.1.1-2
                
                if lamb_web <= 0.673: b_web = w_web
                else: b_web = rho_web*w_web 

                # TRACK SECITON - LIPS ARE NOT CHECKED
                
                #HOLE PROPERTIES
                hole_neg = -1*P_width
                y_hole_neg = (Depth/2)
                
                # generation of table values according to reference example in AISI manual
                '''LENGTHS'''
                top_flange_length = b
                bottom_flange_length = w
                top_web_length = b_web
                bottom_web_length = w_web
                web_hole_neg = hole_neg

                top_inside_corner = (2*pi*(R + t/2))/4
                bottom_inside_corner = (2*pi*(R + t/2))/4

                Length_vector = np.array([top_flange_length, bottom_flange_length, top_web_length, bottom_web_length, web_hole_neg, top_inside_corner, bottom_inside_corner])
                
                Length_sum = Length_vector.sum()

                '''Y DISTANCES'''
                top_flange_ytf = t/2
                bottom_flange_ytf = Depth - t/2
                top_web_ytf = (R+t)+0.5*top_web_length
                bottom_web_ytf = (Depth - (R+t)) - 0.5*(w_web)
                web__hole_neg_ytf = y_hole_neg

                top_inside_corner_ytf = (t + R) - ((2*(R + t/2))/pi); # print(top_inside_corner_ytf)
                bottom_inside_corner_ytf = Depth - ((R + t) - ((2*(R + t/2))/pi)); # print(bottom_inside_corner_ytf)

                ytf_vector = np.array([top_flange_ytf, bottom_flange_ytf, top_web_ytf, bottom_web_ytf, web__hole_neg_ytf, top_inside_corner_ytf, bottom_inside_corner_ytf])
                
                Ly_results = Length_vector*ytf_vector
                
                Ly_total = Ly_results.sum()

                Ly_2_results = Length_vector*(ytf_vector**2)
                Ly_2_total = Ly_2_results.sum()

                '''I'x CALCULATION'''

                top_flange_Ipx = (t**3)/12
                bottom_flange_Ipx = (t**3)/12
                top_web_Ipx = (top_web_length**3)/12
                bottom_web_Ipx = (bottom_web_length**3)/12
                web_hole_neg_Ipx = (web_hole_neg**3)/12

                top_inside_corner_Ipx = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)
                bottom_inside_corner_Ipx = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)

                Ipx_vector = np.array([top_flange_Ipx, bottom_flange_Ipx, top_web_Ipx, bottom_web_Ipx, web_hole_neg_Ipx, top_inside_corner_Ipx, bottom_inside_corner_Ipx])

                Ipx_total = Ipx_vector.sum()

                y_c = Ly_total/Length_sum

                I_xe = (Ipx_total+Ly_2_total-((y_c**2)*Length_sum))*t
                
                S_e = I_xe/y_c
                
                results_list_centroids.append(y_c)
                
                if y_c > centroid:
                    centroid = y_c
                    
                elif y_c < centroid:
                    centroid = y_c
                    ref_stress = (Max_stress*centroid)/(Depth -  centroid)
                    Max_stress = ref_stress
 
        elif section_type == 'stud':
            for i in range(num_iterations):
                # COMPRESSION FLANGE
                w = base - 2*(R + t)
                f = Max_stress
                S = 1.28*sqrt(E/f)
                
                if w/t <= 0.328*S:
                    d_small = (D - R - t)
                    R_I = 1
                    b = w 
                else: 
                    # calculate effective width of flange
                    I_a = min(399*(t**4)*(((w/t)/S)-0.328)**3,(t**4)*(115*((w/t)/S)+5))
                    d_small = (D - R - t)
                    I_s = ((d_small**3)*t*(sin(0.5*pi))**2)/12
                    R_I = min(1, I_s/I_a)
                    n = max(1/3,0.582-((w/t)/(4*S)))
                    
                    if D/w > 0.25 and D/w <= 0.8: k = min(4,(4.82-((5*D)/w))*(R_I**n)+0.43)
                    elif D/w > 0.8: k = 1.25
                    else: k = min(3.75*(R_I**n) + 0.43,4)
                    
                    F_crl = k*(((pi**2)*E)/(12*(1-poiss**2)))*((t/w)**2)
                    lamb_1 = sqrt(f/F_crl)
                    rho = (1-(0.22/lamb_1))/(lamb_1) # Eq.1.1-2
                    
                    if lamb_1 <= 0.673: b = w
                    else: b = rho*w 
                    
                # STIFFENER LIP REVISION
                '''Max stress in lip'''
                f1_lip = f*((centroid)-(t + R))/(centroid)
                f2_lip = f*((centroid)-D)/(centroid)
                f_lip = f1_lip
                psi_lip = abs(f2_lip/f1_lip)
                k_lip = 0.578/(psi_lip+0.34)
                F_crl_lip = k_lip*(((pi**2)*E)/(12*(1-poiss**2)))*((t/d_small)**2)

                lamb_lip = sqrt(f_lip/F_crl_lip)
                rho_lip = (1-(0.22/lamb_lip))/(lamb_lip) # Eq.1.1-2

                if lamb_lip <= 0.673: 
                    d_p_s = 1*d_small # Effective width of stiffener, unmodified
                    if w/t <= 0.328*S:
                        d_s = d_p_s
                    else:
                        d_s = d_p_s*R_I
                else: 
                    d_p_s = rho_lip*d_small # Effective width of stiffener modified accordingly
                    if w/t <= 0.328*S:
                        d_s = d_p_s 
                    else:
                        d_s = d_p_s*R_I

                # WEB REVISION - CONSIDERED AS A UNIFORMLY COMPRESSED UNSTIFFENED ELEMENT - PROCESS ACCOUNTS TOP COMP WEB
                ''' ALTERNATIVE PROCESS CONSIDERING HOLES : STUDS'''
                # if sectionType == 'boxed':
                #     w_web = ((Depth/2) - (d_h/2) - (R + t))
                # else: 
                #     w_web = ((Depth/2) - (d_h/2) - (R + t))*2
                
                w_web = ((Depth/2) - (d_h/2) - (R + t))*2
            
                k_web = 0.43
                f_web = Max_stress
                f1_web = (f_web*((centroid)-(t+R)))/(centroid)
                F_crl_web = k_web*(((pi**2)*E)/(12*(1-poiss**2)))*((t/w_web)**2)
                lamb_web = sqrt(f1_web/F_crl_web)
                rho_web = (1-(0.22/lamb_web))/(lamb_web) # Eq.1.1-2
                
                if lamb_web <= 0.673: b_web = w_web
                else: b_web = rho_web*w_web 
                
                #HOLE PROPERTIES
                hole_neg = -1*P_width
                y_hole_neg = (Depth/2)
                
                # generation of table values according to reference example in AISI manual
                '''LENGTHS'''
                top_flange_length = b
                bottom_flange_length = w
                top_web_length = b_web
                bottom_web_length = w_web
                web_hole_neg = hole_neg

                top_inside_corner = (2*pi*(R + t/2))/4
                bottom_inside_corner = (2*pi*(R + t/2))/4

                top_outside_corner = (2*pi*(R + t/2))/4
                bottom_outside_corner = (2*pi*(R + t/2))/4

                top_lip = d_s
                bottom_lip = d_small

                Length_vector = np.array([top_flange_length, bottom_flange_length, top_web_length, bottom_web_length, web_hole_neg, top_inside_corner, bottom_inside_corner, top_outside_corner, 
                                bottom_outside_corner, top_lip, bottom_lip])
                Length_sum = Length_vector.sum()

                '''Y DISTANCES'''
                top_flange_ytf = t/2
                bottom_flange_ytf = Depth - t/2
                top_web_ytf = (R+t)+0.5*top_web_length
                bottom_web_ytf = (Depth - (R+t)) - 0.5*(w_web) 
                web_hole_neg_ytf = y_hole_neg

                top_inside_corner_ytf = (t + R) - ((2*(R + t/2))/pi); # print(top_inside_corner_ytf)
                bottom_inside_corner_ytf = Depth - ((R + t) - ((2*(R + t/2))/pi)); # print(bottom_inside_corner_ytf)

                top_outside_corner_ytf = (t + R) - ((2*(R + t/2))/pi); # print(top_inside_corner_ytf)
                bottom_outside_corner_ytf = Depth - ((R + t) - ((2*(R + t/2))/pi)); # print(bottom_inside_corner_ytf)

                top_lip_ytf = (R+t)+0.5*(d_s)
                bottom_lip_ytf = Depth - ((R+t)+0.5*(d_small))

                ytf_vector = np.array([top_flange_ytf, bottom_flange_ytf, top_web_ytf, bottom_web_ytf, web_hole_neg_ytf, top_inside_corner_ytf, bottom_inside_corner_ytf, top_outside_corner_ytf,
                                    bottom_outside_corner_ytf, top_lip_ytf, bottom_lip_ytf])

                Ly_results = Length_vector*ytf_vector
                Ly_total = Ly_results.sum()

                Ly_2_results = Length_vector*(ytf_vector**2)
                Ly_2_total = Ly_2_results.sum()

                '''I'x CALCULATION'''

                top_flange_Ipx = (t**3)/12
                bottom_flange_Ipx = (t**3)/12
                top_web_Ipx = (top_web_length**3)/12
                bottom_web_Ipx = (bottom_web_length**3)/12
                web_neg_Ipx = (web_hole_neg**3)/12

                top_inside_corner_Ipx = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)
                bottom_inside_corner_Ipx = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)

                top_outside_corner_Ipx = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)
                bottom_outside_corner_Ipx = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)

                top_lip_Ipx = (top_lip**3)/12
                bottom_lip_Ipx = (bottom_lip**3)/12
                
                ho = Depth
                bo = base

                Ipx_vector = np.array([top_flange_Ipx, bottom_flange_Ipx, top_web_Ipx, bottom_web_Ipx, web_neg_Ipx, top_inside_corner_Ipx, bottom_inside_corner_Ipx, top_outside_corner_Ipx,
                            bottom_outside_corner_Ipx, top_lip_Ipx, bottom_lip_Ipx])

                Ipx_total = Ipx_vector.sum()

                y_c = Ly_total/Length_sum

                I_xe = (Ipx_total+Ly_2_total-((y_c**2)*Length_sum))*t
                
                S_e = I_xe/y_c
                
                results_list_centroids.append(y_c)
                
                if y_c > centroid:
                    # centroid is below centerline
                    centroid = y_c
                    
                elif y_c < centroid:
                    # centroid is above centerline (max stress at bottom flange)
                    centroid = y_c
                    ref_stress = (Max_stress*centroid)/(Depth -  centroid)
                    Max_stress = ref_stress
        else: 
            pass   
    
    if sectionType == 'single':
        return S_e, ho, bo # properties will be used for distortional buckling
    else:
        return S_e

def calculateEffectiveSectionModulus_BuiltUp(Depth, F_n, procedure_alt, section_type, base, R, t, poiss, E, r, D, d_h, P_width, AnalysisType, builtUpSectionType):
    # EFFECTIVE SECTION MODULUS FOR STUD SECTIONS ONLY, NO B2B TRACKS ALLOWED
    centroid = Depth/2
    num_iterations = 20;
    Max_stress = F_n
    results_list_centroids = []
    if procedure_alt == False: # no holes considered
        if section_type == 'track':
            return 'ERROR: Track sections are not allowed when using boxed or back to back members'
        elif section_type == 'stud':
            for i in range(num_iterations):
                # COMPRESSION FLANGE
                w = base - 2*(R + t)
                f = Max_stress
                S = 1.28*sqrt(E/f)
                if w/t <= 0.328*S:
                    d_small = (D - R - t)
                    R_I = 1
                    b = w
                else: 
                    # calculate effective width of flange
                    I_a = min(399*(t**4)*(((w/t)/S)-0.328)**3,(t**4)*(115*((w/t)/S)+5))
                    d_small = (D - R - t)
                    I_s = ((d_small**3)*t*(sin(0.5*pi))**2)/12
                    R_I = min(1, I_s/I_a)
                    n = max(1/3,0.582-((w/t)/(4*S)))
                    
                    if D/w > 0.25 and D/w <= 0.8: k = min(4,(4.82-((5*D)/w))*(R_I**n)+0.43)
                    elif D/w > 0.8: k = 1.25
                    else: k = min(3.75*(R_I**n) + 0.43,4)
                    
                    F_crl = k*(((pi**2)*E)/(12*(1-poiss**2)))*((t/w)**2)
                    lamb_1 = sqrt(f/F_crl)
                    rho = (1-(0.22/lamb_1))/(lamb_1) # Eq.1.1-2
                    
                    if lamb_1 <= 0.673: b = w
                    else: b = rho*w 
                # both compression flanges assumed to have equal correction buckling factors  
                    
                # STIFFENER LIP REVISION
                '''Max stress in lip'''
                f1_lip = f*((centroid)-(t + R))/(centroid)
                f2_lip = f*((centroid)-D)/(centroid)
                f_lip = f1_lip
                psi_lip = abs(f2_lip/f1_lip)
                k_lip = 0.578/(psi_lip+0.34)
                F_crl_lip = k_lip*(((pi**2)*E)/(12*(1-poiss**2)))*((t/d_small)**2)
                # print(F_crl_lip)

                lamb_lip = sqrt(f_lip/F_crl_lip)
                rho_lip = (1-(0.22/lamb_lip))/(lamb_lip) # Eq.1.1-2

                if lamb_lip <= 0.673: 
                    # b_lip = w
                    d_p_s = 1*d_small
                    if w/t <= 0.328*S:
                        d_s = d_p_s 
                    else:
                        d_s = d_p_s*R_I
                else: 
                    # b = rho_lip*w 
                    d_p_s = rho_lip*d_small
                    if w/t <= 0.328*S:
                        d_s = d_p_s 
                    else:
                        d_s = d_p_s*R_I
                # Both compression lips assumed to be the same
                # WEB REVISION
                w_web = Depth - 2*(t + R)
                f_web = Max_stress
                f1_web = (f_web*((centroid)-(t+R)))/(centroid)
                f2_web = -1*((f_web*(Depth - (t+R) - centroid))/centroid)
                psi_web = abs(f2_web/f1_web)

                k_web = 4+2*((1+psi_web)**3)+2*(1+psi_web)
                F_crl_web = k_web*(((pi**2)*E)/(12*(1-poiss**2)))*((t/w_web)**2)
                lamb_web = sqrt(f1_web/F_crl_web)
                rho_web = (1-(0.22/lamb_web))/(lamb_web) # Eq.1.1-2

                if lamb_web <= 0.673: 
                    b_web = w_web
                else: 
                    b_web = rho_web*w_web 
                ho = Depth
                bo = base
                if ho/bo <= 4:
                    b1 = b_web/(3+psi_web)
                    if psi_web > 0.236:
                        b2 = b_web/2
                    else: 
                        b2 = b_web - b1
                else:
                    b1 = b_web/(3+psi_web)
                    b2 = b_web/(1+psi_web)-b1
                # verify if ineffective portion exists for property calculation process
                sum = b1 + b2
                limit = centroid - (t + R)
                if sum <= limit:
                    diff = limit - sum
                else:
                    diff = 0*inch

                # Recompute properties by parts:
                b_neg = -1*(diff) # negative length of negative element, if it exists
                y_b_neg = t/2 + r + b1 - (b_neg/2) # centroidal location of negative component from top fiber

                # generation of table values according to reference example in AISI manual
                '''LENGTHS'''
                top_flange_length = b
                top_flange_length_2 = b
                bottom_flange_length = w
                bottom_flange_length_2 = w
                web_length = w_web
                web_length_2 = w_web
                web_neg_length = b_neg # represents ineffective portion, NOT A HOLE IN THE STUD
                web_neg_length_2 = b_neg

                top_inside_corner = (2*pi*(R + t/2))/4
                top_inside_corner_2 = (2*pi*(R + t/2))/4
                bottom_inside_corner = (2*pi*(R + t/2))/4
                bottom_inside_corner_2 = (2*pi*(R + t/2))/4
                
                top_outside_corner = (2*pi*(R + t/2))/4
                top_outside_corner_2 = (2*pi*(R + t/2))/4
                bottom_outside_corner = (2*pi*(R + t/2))/4
                bottom_outside_corner_2 = (2*pi*(R + t/2))/4

                top_lip = d_s
                top_lip_2 = d_s
                bottom_lip = d_small
                bottom_lip_2 = d_small
                
                Length_vector = np.array([top_flange_length, top_flange_length_2, bottom_flange_length, bottom_flange_length_2, web_length, web_length_2,web_neg_length, web_neg_length_2, top_inside_corner, top_inside_corner_2, bottom_inside_corner, bottom_inside_corner_2, top_outside_corner,top_outside_corner_2, bottom_outside_corner, bottom_outside_corner_2, top_lip, top_lip_2, bottom_lip, bottom_lip_2])
                
                Length_sum = Length_vector.sum()

                '''Y DISTANCES'''
                top_flange_ytf = t/2
                top_flange_ytf_2 = t/2
                bottom_flange_ytf = Depth - t/2
                bottom_flange_ytf_2 = Depth - t/2
                web_ytf = Depth/2
                web_ytf_2 = Depth/2
                web_neg_ytf = y_b_neg
                web_neg_ytf_2 = y_b_neg

                top_inside_corner_ytf = (t + R) - ((2*(R + t/2))/pi); 
                top_inside_corner_ytf_2 = (t + R) - ((2*(R + t/2))/pi); 
                bottom_inside_corner_ytf = Depth - ((R + t) - ((2*(R + t/2))/pi)); 
                bottom_inside_corner_ytf_2 = Depth - ((R + t) - ((2*(R + t/2))/pi)); 

                top_outside_corner_ytf = (t + R) - ((2*(R + t/2))/pi); 
                top_outside_corner_ytf_2 = (t + R) - ((2*(R + t/2))/pi); 
                bottom_outside_corner_ytf = Depth - ((R + t) - ((2*(R + t/2))/pi)); 
                bottom_outside_corner_ytf_2 = Depth - ((R + t) - ((2*(R + t/2))/pi)); 

                top_lip_ytf = (R+t)+0.5*(d_s)
                top_lip_ytf_2 = (R+t)+0.5*(d_s)
                bottom_lip_ytf = Depth - ((R+t)+0.5*(d_small))
                bottom_lip_ytf_2 = Depth - ((R+t)+0.5*(d_small))

                ytf_vector = np.array([top_flange_ytf, top_flange_ytf_2, bottom_flange_ytf, bottom_flange_ytf_2, web_ytf, web_ytf_2, web_neg_ytf, web_neg_ytf_2, top_inside_corner_ytf, top_inside_corner_ytf_2, bottom_inside_corner_ytf, bottom_inside_corner_ytf_2, top_outside_corner_ytf, top_outside_corner_ytf_2, bottom_outside_corner_ytf, bottom_outside_corner_ytf_2, top_lip_ytf, top_lip_ytf_2,bottom_lip_ytf, bottom_lip_ytf_2])

                Ly_results = Length_vector*ytf_vector
                Ly_total = Ly_results.sum()

                Ly_2_results = Length_vector*(ytf_vector**2)
                Ly_2_total = Ly_2_results.sum()

                '''I'x CALCULATION'''

                top_flange_Ipx = (t**3)/12
                top_flange_Ipx_2 = (t**3)/12
                bottom_flange_Ipx = (t**3)/12
                bottom_flange_Ipx_2 = (t**3)/12
                web_Ipx = (web_length**3)/12
                web_Ipx_2 = (web_length**3)/12
                web_neg_Ipx = (web_neg_length**3)/12
                web_neg_Ipx_2 = (web_neg_length**3)/12

                top_inside_corner_Ipx = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)
                top_inside_corner_Ipx_2 = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)
                bottom_inside_corner_Ipx = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)
                bottom_inside_corner_Ipx_2 = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)

                top_outside_corner_Ipx = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)
                top_outside_corner_Ipx_2 = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)
                bottom_outside_corner_Ipx = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)
                bottom_outside_corner_Ipx_2 = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)

                top_lip_Ipx = (top_lip**3)/12
                top_lip_Ipx_2 = (top_lip**3)/12
                bottom_lip_Ipx = (bottom_lip**3)/12
                bottom_lip_Ipx_2 = (bottom_lip**3)/12

                Ipx_vector = np.array([top_flange_Ipx, top_flange_Ipx_2, bottom_flange_Ipx, bottom_flange_Ipx_2, web_Ipx, web_Ipx_2, web_neg_Ipx, web_neg_Ipx_2,top_inside_corner_Ipx, top_inside_corner_Ipx_2, bottom_inside_corner_Ipx, bottom_inside_corner_Ipx_2, top_outside_corner_Ipx, top_outside_corner_Ipx_2, bottom_outside_corner_Ipx, bottom_outside_corner_Ipx_2, top_lip_Ipx, top_lip_Ipx_2, bottom_lip_Ipx, bottom_lip_Ipx_2])

                Ipx_total = Ipx_vector.sum()
                y_c = Ly_total/Length_sum

                I_xe = (Ipx_total+Ly_2_total-((y_c**2)*Length_sum))*t
                
                S_e = I_xe/y_c
                
                results_list_centroids.append(y_c)
                
                if y_c > centroid:
                    # centroid is below centerline
                    centroid = y_c
                    
                elif y_c < centroid:
                    # centroid is above centerline (max stress at bottom flange)
                    centroid = y_c
                    ref_stress = (Max_stress*centroid)/(Depth -  centroid)
                    Max_stress = ref_stress
        
        else:   
            return 'ERROR: Incorrect Section, check source code'
    elif procedure_alt == True:
        if section_type == 'track':
            return 'ERROR: Track sections are not allowed when using boxed or back to back members'             
        elif section_type == 'stud':
            for i in range(num_iterations):
                # COMPRESSION FLANGE
                w = base - 2*(R + t)
                f = Max_stress
                S = 1.28*sqrt(E/f)

                if w/t <= 0.328*S:
                    d_small = (D - R - t)
                    R_I = 1
                    b = w
                else: 
                    # calculate effective width of flange
                    I_a = min(399*(t**4)*(((w/t)/S)-0.328)**3,(t**4)*(115*((w/t)/S)+5))
                    d_small = (D - R - t)
                    I_s = ((d_small**3)*t*(sin(0.5*pi))**2)/12
                    R_I = min(1, I_s/I_a)
                    n = max(1/3,0.582-((w/t)/(4*S)))
                    
                    if D/w > 0.25 and D/w <= 0.8: k = min(4,(4.82-((5*D)/w))*(R_I**n)+0.43)
                    elif D/w > 0.8: k = 1.25
                    else: k = min(3.75*(R_I**n) + 0.43,4)
                    
                    F_crl = k*(((pi**2)*E)/(12*(1-poiss**2)))*((t/w)**2)
                    lamb_1 = sqrt(f/F_crl)
                    rho = (1-(0.22/lamb_1))/(lamb_1) # Eq.1.1-2
                    
                    if lamb_1 <= 0.673: b = w
                    else: b = rho*w 
                    
                # STIFFENER LIP REVISION
                '''Max stress in lip'''
                f1_lip = f*((centroid)-(t + R))/(centroid)
                f2_lip = f*((centroid)-D)/(centroid)
                f_lip = f1_lip
                psi_lip = abs(f2_lip/f1_lip)
                k_lip = 0.578/(psi_lip+0.34)
                F_crl_lip = k_lip*(((pi**2)*E)/(12*(1-poiss**2)))*((t/d_small)**2)

                lamb_lip = sqrt(f_lip/F_crl_lip)
                rho_lip = (1-(0.22/lamb_lip))/(lamb_lip) # Eq.1.1-2

                if lamb_lip <= 0.673: 
                    # b_lip = w
                    d_p_s = 1*d_small # Effective width of stiffener, unmodified
                    if w/t <= 0.328*S:
                        d_s = d_p_s
                    else:
                        d_s = d_p_s*R_I
                else: 
                    # b = rho_lip*w 
                    d_p_s = rho_lip*d_small # Effective width of stiffener modified accordingly
                    if w/t <= 0.328*S:
                        d_s = d_p_s 
                    else:
                        d_s = d_p_s*R_I
                        
                ho = Depth
                bo = base
                # WEB REVISION - CONSIDERED AS A UNIFORMLY COMPRESSED UNSTIFFENED ELEMENT - PROCESS ACCOUNTS TOP COMP WEB
                ''' ALTERNATIVE PROCESS CONSIDERING HOLES : STUDS'''
                w_web = ((Depth/2) - (d_h/2) - (R + t))*2
                k_web = 0.43
                f_web = Max_stress
                f1_web = (f_web*((centroid)-(t+R)))/(centroid)
                F_crl_web = k_web*(((pi**2)*E)/(12*(1-poiss**2)))*((t/w_web)**2)
                lamb_web = sqrt(f1_web/F_crl_web)
                rho_web = (1-(0.22/lamb_web))/(lamb_web) # Eq.1.1-2
                
                if lamb_web <= 0.673: b_web = w_web
                else: b_web = rho_web*w_web 
                
                #HOLE PROPERTIES
                hole_neg = -1*P_width
                y_hole_neg = (Depth/2)
                
                # generation of table values according to reference example in AISI manual
                '''LENGTHS'''
                top_flange_length = b
                top_flange_length_2 = b
                bottom_flange_length = w
                bottom_flange_length_2 = w
                top_web_length = b_web
                top_web_length_2 = b_web
                bottom_web_length = w_web
                bottom_web_length_2 = w_web
                web_hole_neg = hole_neg
                web_hole_neg_2 = hole_neg

                top_inside_corner = (2*pi*(R + t/2))/4
                top_inside_corner_2 = (2*pi*(R + t/2))/4
                bottom_inside_corner = (2*pi*(R + t/2))/4
                bottom_inside_corner_2 = (2*pi*(R + t/2))/4

                top_outside_corner = (2*pi*(R + t/2))/4
                top_outside_corner_2 = (2*pi*(R + t/2))/4
                bottom_outside_corner = (2*pi*(R + t/2))/4
                bottom_outside_corner_2 = (2*pi*(R + t/2))/4

                top_lip = d_s
                top_lip_2 = d_s
                bottom_lip = d_small
                bottom_lip_2 = d_small

                Length_vector = np.array([top_flange_length,top_flange_length_2, bottom_flange_length,bottom_flange_length_2, top_web_length,top_web_length_2, bottom_web_length,bottom_web_length_2, web_hole_neg,web_hole_neg_2, top_inside_corner,top_inside_corner_2, bottom_inside_corner,bottom_inside_corner_2, top_outside_corner,top_outside_corner_2, 
                                bottom_outside_corner,bottom_outside_corner_2, top_lip, top_lip_2, bottom_lip,  bottom_lip_2])
                Length_sum = Length_vector.sum()

                '''Y DISTANCES'''
                top_flange_ytf = t/2
                top_flange_ytf_2 = t/2
                bottom_flange_ytf = Depth - t/2
                bottom_flange_ytf_2 = Depth - t/2
                top_web_ytf = (R+t)+0.5*top_web_length
                top_web_ytf_2 = (R+t)+0.5*top_web_length
                bottom_web_ytf = (Depth - (R+t)) - 0.5*(w_web) 
                bottom_web_ytf_2 = (Depth - (R+t)) - 0.5*(w_web) 
                web_hole_neg_ytf = y_hole_neg
                web_hole_neg_ytf_2 = y_hole_neg

                top_inside_corner_ytf = (t + R) - ((2*(R + t/2))/pi); # print(top_inside_corner_ytf)
                top_inside_corner_ytf_2 = (t + R) - ((2*(R + t/2))/pi);
                bottom_inside_corner_ytf = Depth - ((R + t) - ((2*(R + t/2))/pi)); # print(bottom_inside_corner_ytf)
                bottom_inside_corner_ytf_2 = Depth - ((R + t) - ((2*(R + t/2))/pi));

                top_outside_corner_ytf = (t + R) - ((2*(R + t/2))/pi); # print(top_inside_corner_ytf)
                top_outside_corner_ytf_2 = (t + R) - ((2*(R + t/2))/pi);
                bottom_outside_corner_ytf = Depth - ((R + t) - ((2*(R + t/2))/pi)); # print(bottom_inside_corner_ytf)
                bottom_outside_corner_ytf_2 = Depth - ((R + t) - ((2*(R + t/2))/pi));

                top_lip_ytf = (R+t)+0.5*(d_s)
                top_lip_ytf_2 = (R+t)+0.5*(d_s)
                bottom_lip_ytf = Depth - ((R+t)+0.5*(d_small))
                bottom_lip_ytf_2 = Depth - ((R+t)+0.5*(d_small))

                ytf_vector = np.array([top_flange_ytf,top_flange_ytf_2, bottom_flange_ytf,bottom_flange_ytf_2, top_web_ytf,top_web_ytf_2, bottom_web_ytf,bottom_web_ytf_2, web_hole_neg_ytf,web_hole_neg_ytf_2, top_inside_corner_ytf,top_inside_corner_ytf_2, bottom_inside_corner_ytf,bottom_inside_corner_ytf_2, top_outside_corner_ytf,top_outside_corner_ytf_2,
                                    bottom_outside_corner_ytf,bottom_outside_corner_ytf_2, top_lip_ytf,top_lip_ytf_2, bottom_lip_ytf,bottom_lip_ytf_2])

                Ly_results = Length_vector*ytf_vector
                Ly_total = Ly_results.sum()

                Ly_2_results = Length_vector*(ytf_vector**2)
                Ly_2_total = Ly_2_results.sum()

                '''I'x CALCULATION'''

                top_flange_Ipx = (t**3)/12
                top_flange_Ipx_2 = (t**3)/12
                bottom_flange_Ipx = (t**3)/12
                bottom_flange_Ipx_2 = (t**3)/12
                top_web_Ipx = (top_web_length**3)/12
                top_web_Ipx_2 = (top_web_length**3)/12
                bottom_web_Ipx = (bottom_web_length**3)/12
                bottom_web_Ipx_2 = (bottom_web_length**3)/12
                web_neg_Ipx = (web_hole_neg**3)/12
                web_neg_Ipx_2 = (web_hole_neg**3)/12

                top_inside_corner_Ipx = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)
                top_inside_corner_Ipx_2 = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)
                bottom_inside_corner_Ipx = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)
                bottom_inside_corner_Ipx_2 = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)

                top_outside_corner_Ipx = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)
                top_outside_corner_Ipx_2 = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)
                bottom_outside_corner_Ipx = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)
                bottom_outside_corner_Ipx_2 = ((pi/16)*(R+t)**3)-((pi/16)*(R)**3)

                top_lip_Ipx = (top_lip**3)/12
                top_lip_Ipx_2 = (top_lip**3)/12
                bottom_lip_Ipx = (bottom_lip**3)/12
                bottom_lip_Ipx_2 = (bottom_lip**3)/12

                Ipx_vector = np.array([top_flange_Ipx,top_flange_Ipx_2, bottom_flange_Ipx,bottom_flange_Ipx_2, top_web_Ipx,top_web_Ipx_2, bottom_web_Ipx,bottom_web_Ipx_2, web_neg_Ipx,web_neg_Ipx_2, top_inside_corner_Ipx,top_inside_corner_Ipx_2, bottom_inside_corner_Ipx,bottom_inside_corner_Ipx_2, top_outside_corner_Ipx,top_outside_corner_Ipx_2,
                            bottom_outside_corner_Ipx,bottom_outside_corner_Ipx_2, top_lip_Ipx,top_lip_Ipx_2, bottom_lip_Ipx,bottom_lip_Ipx_2])

                Ipx_total = Ipx_vector.sum()
                y_c = Ly_total/Length_sum
                I_xe = (Ipx_total+Ly_2_total-((y_c**2)*Length_sum))*t
                S_e = I_xe/y_c
                results_list_centroids.append(y_c)
                
                if y_c > centroid:
                    # centroid is below centerline
                    centroid = y_c
                    
                elif y_c < centroid:
                    # centroid is above centerline (max stress at bottom flange)
                    centroid = y_c
                    ref_stress = (Max_stress*centroid)/(Depth -  centroid)
                    Max_stress = ref_stress

        else: 
            return 'ERROR: Incorrect section type value' 
              
    if AnalysisType == 'F_n' and builtUpSectionType == 'b2b':
        return S_e, ho, bo
    else: 
        return S_e

def getSingleShear(target_id, section_type, Fy, E, G, P_width, P_cond):
    '''Calculation of CFS single stud shear capacity
    Calculations apply to both stud and track sections. 
    Args:
        target_id (str): Section shape name. 
        section_type (str): 'stud' or 'track'.
        Fy (float): Yield stress (ksi).
        E (float): Elastic modulus (ksi).
        G (float): Shear Modulus (ksi).
        P_width (float): Punchout width, typically 1.5 inches (in).
        P_cond (bool): Punchout condition.

    Returns:
        forallpeople.Physical: Allowable shear strength single stud/track.
    '''
    poiss = 0.3
    Fy = Fy*ksi
    E = E*ksi
    G = G*ksi
    P_width = P_width*inch
    s_section_path = os.path.dirname(__file__) + '\CFS_S_SECTION_DATA.csv'
    track_section_path = os.path.dirname(__file__) + '\CFS_TRACK_SECTION_DATA.csv'
    rawPath_s_section = r'' + s_section_path
    rawPath_track_section = r'' + track_section_path
    
    list_sections = []
    if section_type == 'stud':
        with open(rawPath_s_section, 'r', encoding="utf-8-sig") as csv_file:
            csv_reader = csv.DictReader(csv_file)
            keys = csv_reader.fieldnames
            
            for line in csv_reader:
                list_sections.append(line['ID'])
            
            if target_id not in list_sections:
                return 'ERROR: Section was not found in current reference file. Check your input or add the section manually.'
            else:
                location = list_sections.index(target_id)
                csv_file.seek(0)
                next(csv_reader)
                for count,line in enumerate(csv_reader):
                    if count == location:
                        Depth = float(line[keys[1]])
                        t = float(line[keys[3]])
                        R = float(line[keys[5]])
                        x_m = float(line[keys[14]])
                        m = float(line[keys[15]])
                        J = float(line[keys[16]])
    else:
        with open(rawPath_track_section, 'r', encoding="utf-8-sig") as csv_file:
            csv_reader = csv.DictReader(csv_file)
            keys = csv_reader.fieldnames
            
            for line in csv_reader:
                list_sections.append(line['ID'])
            if target_id not in list_sections:
                return 'ERROR: Section was not found in current reference file. Check your input or add the section manually.'
            else: 
                location = list_sections.index(target_id)
                csv_file.seek(0)
                next(csv_reader)
                for count,line in enumerate(csv_reader):
                    if count == location:
                        Depth = float(line[keys[1]])
                        t = float(line[keys[3]])
                        R = float(line[keys[4]])
                        x_m = float(line[keys[13]])
                        m = float(line[keys[14]])
                        J = float(line[keys[15]])

    if target_id in list_sections:
        Depth = Depth*inch
        t = t*inch
        R = R*inch
        x_m = x_m*inch
        J = J*inch**4
        m = m*inch
    
    # SHEAR STRENGTH CALCULATIONS:
    
    # Shear strength of web elements:   
    h = Depth - 2*(R+t)
    A_w = h * t #Area of web element
    k_v = 5.34 #Unreinforced webs
    
    F_cr = ((pi**2)*E*k_v)/(12*(1-(poiss**2))*((h/t)**2))
    V_cr = A_w*F_cr #Shear Buckling Force
    V_y = 0.6*A_w*Fy
    lamb_v = sqrt(V_y/V_cr)
    
    if lamb_v <= 0.815: V_n = V_y
    elif 0.815 < lamb_v <= 1.227: V_n = 0.815*sqrt(V_cr*V_y)*1*kip
    elif lamb_v > 1.227: V_n = V_cr
    
    if P_cond == True:
        c = h/2 - P_width/2 #For non circular holes
        if c/t >= 54: q_s = 1
        elif 5 <= c/t < 54: q_s = c/(54*t)
        V_n = V_n*q_s

    V_ally = (V_n)/(1.60)
        
    FinalShearStrength = V_ally.to('kip')
    
    return FinalShearStrength

# AXIAL STRENGTH FUNCTIONS:------------------------------------------------------------------------------------
def getAxialStrength_Single(target_id:str, section_type:str, Fy:float, S_bracing:float, E:float, G:float, P_width:float, P_cond:bool, K_x:float, K_y:float, K_t:float, L_x:float, L_t:float, dist_bracing=False, d_bracing=0.01):
    '''Calculation of CFS single stud/track capacity
    Calculations apply to stud and track sections. Holes are supported for both studs and tracks.
    Args:
        target_id (str): Section shape name. 
        section_type (str): 'stud' or 'track'.
        Fy (float): Yield stress (ksi).
        S_bracing (float): Weak axis bracing (ft).
        E (float): Elastic modulus (ksi).
        G (float): Shear Modulus (ksi).
        P_width (float): Punchout width, typically 1.5 inches (in).
        P_cond (bool): Punchout condition.
        K_x (float): Effective length factor x-axis.
        K_y (float): Effective length factor y-axis.
        K_t (float): Effective length factor for twisting/torsional buckling.
        Lx (float): Strong axis unbraced length (ft).
        L_t (float): Twisting unbraced length (ft).
        dist_bracing (bool): Distortional buckling bracing status.
        d_bracing (float): Distortional buckling bracing value if True (ft).
        
        For distortional buckling, it is recommended to use the same value as the axial bracing, unless a conservative approach is preferred.  

    Returns:
        forallpeople.Physical: Allowable axial strength single stud.
    '''
    # Get Section properties from CSV file:
    target_id = str(target_id)
    Fy = Fy*ksi
    S_bracing = S_bracing*ft
    E = E*ksi
    G = G*ksi
    P_width = P_width*inch
    L_x = L_x*ft
    L_y = S_bracing
    L_t = S_bracing
    poiss = 0.3
    if dist_bracing == True:
        L_m = d_bracing
    else:
        L_m = L_x
    list_sections = []
    s_section_path = os.path.dirname(__file__) + '\CFS_S_SECTION_DATA.csv'
    track_section_path = os.path.dirname(__file__) + '\CFS_TRACK_SECTION_DATA.csv'
    rawPath_s_section = r'' + s_section_path
    rawPath_track_section = r'' + track_section_path
    if section_type == 'stud':
        try:
            with open(rawPath_s_section, 'r', encoding="utf-8-sig") as csv_file:
                csv_reader = csv.DictReader(csv_file)
                keys = csv_reader.fieldnames
                
                for line in csv_reader:
                    list_sections.append(line['ID'])
                
                if target_id not in list_sections:
                    return "ERROR: INVALID SECTION"
                else:
                    location = list_sections.index(target_id)
                    csv_file.seek(0)
                    next(csv_reader)
                    for count,line in enumerate(csv_reader):
                        if count == location:
                            Depth = float(line[keys[1]])
                            base = float(line[keys[2]])
                            t = float(line[keys[3]])
                            D = float(line[keys[4]])
                            R = float(line[keys[5]])
                            Area = float(line[keys[6]])
                            Ix = float(line[keys[8]])
                            Iy = float(line[keys[11]])
                            x_m = float(line[keys[14]])
                            m = float(line[keys[15]])
                            J = float(line[keys[16]])
                            Cw = float(line[keys[17]])
                            xo = float(line[keys[20]])
        except FileNotFoundError as error:
            return 'FILE NOT FOUND ' + rawPath_s_section
    else:
        try:
            with open(rawPath_track_section, 'r', encoding="utf-8-sig") as csv_file:
                csv_reader = csv.DictReader(csv_file)
                keys = csv_reader.fieldnames
                for line in csv_reader:
                    list_sections.append(line['ID'])
                if target_id not in list_sections:
                    return "ERROR: INVALID SECTION"
                else: 
                    location = list_sections.index(target_id)
                    csv_file.seek(0)
                    next(csv_reader)
                    for count,line in enumerate(csv_reader):
                        if count == location:
                            Depth = float(line[keys[1]])
                            base = float(line[keys[2]])
                            t = float(line[keys[3]])
                            R = float(line[keys[4]])
                            Area = float(line[keys[5]])
                            Ix = float(line[keys[7]])
                            Iy = float(line[keys[10]])
                            x_m = float(line[keys[13]])
                            m = float(line[keys[14]])
                            J = float(line[keys[15]])
                            Cw = float(line[keys[16]])
                            xo = float(line[keys[19]])
        except FileNotFoundError:
            return 'FILE NOT FOUND ' + rawPath_track_section

    if target_id in list_sections:
        Depth = Depth*inch
        Ag = Area*inch**2
        Ixxg = Ix*inch**4
        Iyyg = Iy*inch**4
        x_o = xo*inch
        Cw = Cw*inch**6
        base = base*inch
        t = t*inch
        R = R*inch
        if section_type == 'stud':
            D = D*inch # this is the small d value - lip dimension
        x_m = x_m*inch
        J = J*inch**4
        r_x = sqrt(Ixxg/Ag)*inch
        r_y = sqrt(Iyyg/Ag)*inch
        r_o = sqrt((r_x**2)+(r_y**2)+(x_o**2))*inch
        m = m*inch
    
    # CHECK FLEXURAL BUCKLING
    x_slender = (K_x*L_x)/(r_x)
    y_slender = (K_y*L_y)/(r_y)
    if x_slender > 200 or y_slender > 200: 
        Slender_status = "NG, TOO SLENDER"
        return Slender_status
    else:
        # PROCEED WITH CALCULATIONS
        Slender_status = "SLENDERNESS OK"
        crit_slender = max(x_slender, y_slender)
        F_cre = (((pi)**2)*E)/(crit_slender**2)
        
        # CHECK TORSIONAL FLEXURAL BUCKLING
        beta = 1 - ((x_o/r_o)**2)
        sigma_ex = (((pi)**2)*E)/(x_slender**2)
        sigma_t = (1/(Ag*r_o**2))*(G*J+(((pi**2)*E*Cw)/((K_t*L_t)**2)))
        F_creftb = (1/(2*beta))*((sigma_ex+sigma_t)-sqrt(((sigma_ex+sigma_t)**2)-4*beta*sigma_ex*sigma_t))
        
        # CONTROLLING BUCKLING MODE REVISION:
        F_cre_crit = min(F_cre,F_creftb)
        lamb_c = sqrt(Fy/F_cre_crit)
        if lamb_c <= 1.5: F_n = (0.658**(lamb_c**2))*Fy
        else: F_n = (0.877/(lamb_c**2))*Fy

        # EFFECTIVE AREA AT Fn CALCULATION:
        if section_type == 'stud':
            w = base - 2*(R+t)
            rel = (w/t)
            d = D-R-t
            f = F_n
            I_s = ((d**3)*t*(sin(0.5*pi))**2)/12
            S = 1.28*sqrt(E/f)
            if rel<=(0.328*S): I_a = 0*inch**4 
            else: I_a = min(399*(t**4)*((rel/S)-0.328)**3,(t**4)*(115*(rel/S)+5))

            if rel<=(0.328*S): R_I = 0
            else: R_I = min((I_s/I_a),1)

            relDw = min(D/w,0.8)
            n = max(1/3,0.582-(rel/(4*S)))

            if relDw > 0.25 and relDw <= 0.8: k = min(4,4.82-((5*D)/w)*(R_I**n)+0.43)
            elif relDw > 0.8: k = 1.25
            else: k = min(3.75*(R_I**n) + 0.43,4)
            F_crl = k*(((pi**2)*E)/(12*(1-poiss**2)))*((t/w)**2)
            lamb_st = sqrt(f/F_crl)
            rho = (1-(0.22/lamb_st))/(lamb_st) # Eq.1.1-2

            if lamb_st <= 0.673: b = w
            else: b = rho*w
            
            # (a) CHECK STIFFENER LIP USING APPENDIX 1 SECTION 1.2.1
            k = 0.43
            w_lip = D - t - R
            F_crl_2 = k*(((pi**2)*E)/(12*(1-poiss**2)))*((t/w_lip)**2)
            lamb_lip = sqrt(f/F_crl_2)
            rho = (1-(0.22/lamb_lip))/(lamb_lip)

            # Use of Equations in section 1.3 of Appendix 1
            dps = d
            if rel<=(0.328*S): d_lip = dps
            else: d_lip = dps*R_I

            A_lips = 2*d_lip*t
            A_w_lips = 2*d*t
            A_w_lips_diff = A_w_lips - A_lips
            
        elif section_type == 'track':
            w = base - (R + t)
            f = F_n
            k_flange = 0.43
            F_crl = k_flange*(((pi**2)*E)/(12*(1-poiss**2)))*((t/w)**2)
            lamb_1 = sqrt(f/F_crl)
            rho = (1-(0.22/lamb_1))/(lamb_1) # Eq.1.1-2
            if lamb_1 <= 0.673: b = w
            else: b = rho*w
            
        else:
            return 'ERROR: Incorrect section type value'

        A_flanges = 2*b*t
        A_wflange = 2*w*t
        A_dif_flange = A_wflange - A_flanges

        # WEB REVISION
        # If punchout is true, treat web as two unstiffened elements, one on each side of the 1.50 inch wide punchout.
        if P_cond == True: w_web = (Depth - 2*(R + t) - P_width)/2; k_web = 0.43
        else: w_web = (Depth - 2*(R + t)); k_web = 4
        F_cr_web = k_web*(((pi**2)*E)/(12*(1-poiss**2)))*((t/w_web)**2)
        lamb_web = sqrt(f/F_cr_web)
        rho_web = (1-(0.22/lamb_web))/lamb_web
        if lamb_web <= 0.673: b_web = w_web
        else: b_web = rho_web*w_web

        if P_cond == True: A_web_ineff = 2 * (w_web - b_web) * t
        else: A_web_ineff = (w_web - b_web) * t


        if section_type == 'stud':
            if P_cond == True: A_effective = Ag-A_dif_flange-A_w_lips_diff-A_web_ineff-P_width*t
            else: A_effective = Ag-A_dif_flange-A_w_lips_diff-A_web_ineff
        
        elif section_type == 'track':
            if P_cond == True: A_effective = Ag-A_dif_flange-A_web_ineff-P_width*t
            else: A_effective = Ag-A_dif_flange-A_web_ineff
        
        A_etotal = A_effective
        
        # DISTORTIONAL BUCKLING:
        if section_type == 'stud':
            d_midline = D - (t/2)
            b_midline = base - t
            h_o = Depth
            A_web_gross = (Depth-2*(R+t))*t
            A_web_net = 2*w_web*t

            if P_cond == True: tr = t*((A_web_net/A_web_gross)**(1/3))
            else: tr = t
            A_net = (Ag-(t*P_width))

            # * Calculation of basic distortional buckling geometric properties
            A_f = (b_midline+d_midline)*tr
            J_f = (1/3)*b_midline*(tr**3)+(1/3)*d_midline*(tr**3)
            I_xf = (tr*((tr**2)*(b_midline**2)+4*b_midline*(d_midline**3)+(tr**2)*b_midline*d_midline+(d_midline**4)))/(12*(b_midline+d_midline))
            I_yf = (tr*((b_midline**4)+4*d_midline*(b_midline**3)))/(12*(b_midline+d_midline))
            I_xyf = (tr*(b_midline**2)*(d_midline**2))/(4*(b_midline+d_midline))
            Cw_f = 0*inch**6
            h_xf = (-1*((b_midline**2)+2*b_midline*d_midline))/(2*(b_midline+d_midline))
            x_of = ((b_midline**2))/(2*(b_midline+d_midline))
            h_yf = (-1*(d_midline**2))/(2*(b_midline+d_midline))
            y_of = h_yf

            L_crd = (((6*(pi**4)*h_o*(1-(poiss**2)))/(tr**3))*(I_xf*((x_of-h_xf)**2)+Cw_f-((I_xyf**2)/I_yf)*((x_of-h_xf)**2)))**(1/4) # Critical unbraced length of distortional buckling
            L_crd = ((L_crd*(1/(1*mm)))/(304.8))*(1*ft) #Ignore this calculation, just a unit calc proof

            # L_m = S_bracing # Distance between discrete restraints that restrict buckling
            L_min = min(L_m,L_crd)
            K_phi_fe = (((pi)/L_min)**4)*(E*I_xf*((x_of-h_xf)**2)+E*Cw_f-E*((I_xyf**2)/I_yf)*((x_of-h_xf)**2))+(((pi)/L_min)**2)*G*J_f
            K_phi_fe = (K_phi_fe*(1/(4448.2216153*N)))*1*kip # Ignore, just a unit conversion

            K_phi_we = (E*(tr**3))/(6*h_o*(1-(poiss**2)))
            K_phi = 0*kip
            K_phi_fg = ((pi/L_min)**2)*(A_f*(((x_of-h_xf)**2)*(((I_xyf)/(I_yf))**2)-2*y_of*(x_of-h_xf)*(I_xyf/I_yf)+(h_xf**2)+(y_of**2))+I_xf+I_yf)
            K_phi_fg = K_phi_fg*(1/(645.16*mm**2))*1*inch**2 # Ignore, unit conversions
            K_phi_wg = ((pi/L_min)**2)*((tr*(h_o**3))/(60))
            K_phi_wg = K_phi_wg*(1/(645.16*mm**2))*1*inch**2 # Ignore, unit conversion
            F_crd = (K_phi_fe+K_phi_we+K_phi)/(K_phi_fg+K_phi_wg)

            P_crd = Ag*F_crd
            P_y = Ag*Fy
            P_ynet = A_net*Fy
            lamb_d = sqrt(P_y/P_crd)
            lamb_d1 = 0.561*(P_ynet/P_y)
            lamb_d2 = 0.561*(14*((P_y/P_ynet)**0.4)-13)
            # * Distortional Buckling Allowable Axial Load with No Holes:

            if lamb_d <= 0.561: P_ndNH = P_y
            else: P_ndNH = (1-0.25*((P_crd/P_y)**0.6))*((P_crd/P_y)**0.6)*P_y

            # * Distortional Buckling Allowable Axial Load WITH Holes:
            P_d2 = (1-0.25*((1/lamb_d2)**1.2))*((1/lamb_d2)**1.2)*P_y

            if lamb_d <= lamb_d2:
                if lamb_d <= lamb_d1: P_ndWH = P_ynet
                else: P_ndWH = P_ynet-((P_ynet-P_d2)/(lamb_d2-lamb_d1))*(lamb_d-lamb_d1)
            elif lamb_d > lamb_d2:
                if lamb_d <= 0.561:
                    P_ndWH = P_y
                else: P_ndWH = (1-0.25*((P_crd/P_y)**0.6))*((P_crd/P_y)**0.6)*P_y
        
            if P_ndWH <= P_ndNH: 
                P_ndfinal = P_ndWH if P_cond == True else P_ndNH
                # COMPUTE Pn:
                P_nl = A_etotal*F_n
                P_n = min(P_nl,P_ndfinal)
                P_a = P_n/1.80
                Final_Pa = P_a.to('kip')
                return Final_Pa
            else: 
                return "ERROR - CAPACITY WITH HOLES IS HIGHER COMPARED TO CAPACITY WITH NO HOLES - CHECK SCRIPT" 
        elif section_type == 'track':
            P_nl = A_etotal*F_n
            P_a = P_nl/1.80
            Final_Pa = P_a.to('kip')
            return Final_Pa
        else:
            return 'ERROR: Incorrect section type value'

def getAxialStrength_Boxed(target_id:str, section_type:str, Fy:float, S_bracing:float, E:float, G:float, P_width:float, P_cond:bool, K_x:float, K_y:float, L_x:float, a:float):
    '''Calculation of CFS built-up boxed member axial capacity.
    Calculations apply ONLY to stud sections. Tracks will be added in the future.
    
    Args:
        target_id (str): Section shape name. 
        section_type (str): 'stud' or 'track'.
        Fy (float): Yield stress (ksi).
        S_bracing (float): Weak axis bracing (ft).
        E (float): Elastic modulus (ksi).
        G (float): Shear Modulus (ksi).
        P_width (float): Punchout width, typically 1.5 inches (in).
        P_cond (bool): Punchout condition.
        K_x (float): Effective length factor x-axis.
        K_y (float): Effective length factor y-axis.
        L_x (float): Strong axis unbraced length (ft).
        a (float): Interconnection spacing (in)
        
        No distortional buckling is considered in boxed stud sections

    Returns:
        forallpeople.Physical: Allowable axial strength boxed studs.
    '''
    # Get Section properties from CSV file:
    target_id = str(target_id)
    Fy = Fy*ksi
    S_bracing = S_bracing*ft
    E = E*ksi
    G = G*ksi
    P_width = P_width*inch
    L_x = L_x*ft
    L_y = S_bracing
    poiss = 0.3
    a = a*inch
    list_sections = []
    s_section_path = os.path.dirname(__file__) + '\CFS_S_SECTION_DATA.csv'
    rawPath_s_section = r'' + s_section_path
    if section_type == 'stud':
        try:
            with open(rawPath_s_section, 'r', encoding="utf-8-sig") as csv_file:
                csv_reader = csv.DictReader(csv_file)
                keys = csv_reader.fieldnames
                
                for line in csv_reader:
                    list_sections.append(line['ID'])
                
                if target_id not in list_sections:
                    return "ERROR: INVALID SECTION"
                else:
                    location = list_sections.index(target_id)
                    csv_file.seek(0)
                    next(csv_reader)
                    for count,line in enumerate(csv_reader):
                        if count == location:
                            Depth = float(line[keys[1]])
                            base = float(line[keys[2]])
                            t = float(line[keys[3]])
                            D = float(line[keys[4]])
                            R = float(line[keys[5]])
                            Area = float(line[keys[6]])
                            Ix = float(line[keys[8]])
                            Iy = float(line[keys[11]])
                            x_m = float(line[keys[14]])
                            m = float(line[keys[15]])
                            J = float(line[keys[16]])
                            Cw = float(line[keys[17]])
        except FileNotFoundError as error:
            return 'FILE NOT FOUND ' + rawPath_s_section
    else:
        return 'ERROR: Incorrect section used'

    if target_id in list_sections:
        Depth = Depth*inch
        Ag = Area*inch**2
        Ixxg = Ix*inch**4
        Iyyg = Iy*inch**4
        Cw = Cw*inch**6
        base = base*inch
        t = t*inch
        R = R*inch
        if section_type == 'stud':
            D = D*inch # this is the small d value - lip dimension
        x_m = x_m*inch
        J = J*inch**4
        r_y = sqrt(Iyyg/Ag)*inch
        m = m*inch

    # BUILT UP SECTION PROPERTIES: 
    Atg = 2*Ag
    Ixxgt = 2*Ixxg
    Iyygt = 2*(Iyyg+Ag*((base-x_m)**2))
    rxt = sqrt(Ixxgt/Atg)*inch
    ryt = sqrt(Iyygt/Atg)*inch
    
    # FLEXURAL BUCKLING STRESS:
    x_slender = (K_x*L_x)/(rxt)
    r_i = r_y
    F_cre = (((pi)**2)*E)/(x_slender**2)
    y_slender = (K_y*L_y)/(ryt)
    y_slendermod = sqrt((y_slender**2)+(a/r_i)**2)
    crit_slender = max(x_slender, y_slendermod)
    ratio_1 = (a/r_i)
    
    if ratio_1 <= (0.5*crit_slender):
        pass
    else: 
        status = 'NG'
        return f'Interconnection spacing is {status}'

    if x_slender > 200 and y_slender > 200:
        return 'Slenderness check failed, too slender'
        
    F_cre_2 = (((pi)**2)*E)/(y_slendermod**2)
    F_cre_crit = min(F_cre, F_cre_2)
    lamb_c = sqrt(Fy/F_cre_crit)
    
    # FN CALCULATION:
    if lamb_c <= 1.5: Fn = (0.658**(lamb_c**2))*Fy
    else: Fn = (0.877/(lamb_c**2))*Fy
    
    # EFFECTIVE AREA AT FN:
    # - Check flange as uniformly compressed element with an edge stiffener    
    w = base - 2*(R+t)
    rel = w/t
    d1 = D-R-t
    f = Fn
    I_s = ((d1**3)*t*(sin(0.5*pi))**2)/12
    S = 1.28*sqrt(E/f)
    if rel <= (0.328*S): I_a = 0*inch**4
    else: I_a = min(399*(t**4)*(((w/t)/S)-0.328)**3,(t**4)*(115*((w/t)/S)+5))
    if rel <= (0.328*S): R_I = 0
    else: R_I = (I_s/I_a)
    if R_I > 1: R_I = 1
    relDw = (D/w)
    if relDw > 0.8: relDw = 0.8
    n = max((1/3), 0.582-((w/t)/(4*S)))
    
    if relDw > 0.25 and relDw <= 0.8: k = min((4.82-((5*D)/w))*((R_I)**n)+0.43,4)
    else: k = min(3.75*(R_I**n)+0.43,4)
    
    F_cr = k*(((pi**2)*E)/(12*(1-(poiss**2))))*((t/w)**2)
    lamb_gen = sqrt(f/F_cr)
    rho = (1-(0.22/lamb_gen))/lamb_gen
    
    if lamb_gen <= 0.673: b = w
    else: b = rho*w
    A_flanges = 2*(b*t)
    A_w_flange = 2*w*t
    A_dif_flange = A_w_flange - A_flanges
    d_s = d1
    if rel <= (0.328*S): d_lip = d_s
    else: d_lip = d_s*R_I
    # Two lips considered for the typical stud section
    A_lips = 2*(d_lip*t)
    A_w_lips = 2*d1*t
    A_w_lips_diff = A_w_lips-A_lips
    
    # - Check web:
    if P_cond == True:
        w_web = (Depth-2*(R+t)-P_width)/2
        k_web = 0.43
    else: 
        w_web = (Depth-2*(R+t))
        k_web = 4.0
    
    F_cr_web = k_web*(((pi**2)*E)/(12*(1-(poiss**2))))*((t/w_web)**2)
    lamb_web = sqrt(f/F_cr_web)
    rho_web = ((1-(0.22/lamb_web))/lamb_web)
    if lamb_web <= 0.673: b_web = w_web
    else: b_web = rho_web*w_web
    if P_cond == True: A_web_ineff = (2*(w_web-b_web))*t
    else: A_web_ineff = (w_web-b_web)*t
    
    if P_cond == True: A_effective = Ag - A_dif_flange - A_w_lips_diff - A_web_ineff - (P_width*t)
    else: A_effective = Ag - A_dif_flange - A_w_lips_diff - A_web_ineff

    A_e_total = 2*A_effective
    P_nl = A_e_total * Fn
    
    # ALLOWABLE AXIAL LOAD:
    P_a = P_nl/1.80
    P_a_final = P_a.to('kip')
    return P_a_final

def getAxialStrength_B2B(target_id:str, section_type:str, Fy:float, S_bracing:float, E:float, G:float, P_width:float, P_cond:bool, K_x:float, K_y:float, K_t:float, L_x:float, L_t:float, a:float):
    '''Calculation of CFS built-up back to back member axial capacity.
    Calculations apply ONLY to stud sections. Tracks will be added in the future.  
    
    Args:
        target_id (str): Section shape name. 
        section_type (str): 'stud' or 'track'.
        Fy (float): Yield stress (ksi).
        S_bracing (float): Weak axis bracing (ft).
        E (float): Elastic modulus (ksi).
        G (float): Shear Modulus (ksi).
        P_width (float): Punchout width, typically 1.5 inches (in).
        P_cond (bool): Punchout condition.
        K_x (float): Effective length factor x-axis.
        K_y (float): Effective length factor y-axis.
        K_t (float): Effective length factor for twisting.
        L_x (float): Strong axis unbraced length (ft).
        L_t (float): Twisting unbraced length (ft).
        a (float): Interconnection spacing (in).

    Returns:
        forallpeople.Physical: Allowable axial strength back to back stud member.
    '''
    # Get Section properties from CSV file:
    target_id = str(target_id)
    Fy = Fy*ksi
    S_bracing = S_bracing*ft
    E = E*ksi
    G = G*ksi
    P_width = P_width*inch
    L_x = L_x*ft
    L_y = S_bracing
    L_t = L_t*ft
    poiss = 0.3
    a = a*inch
    list_sections = []
    s_section_path = os.path.dirname(__file__) + '\CFS_S_SECTION_DATA.csv'
    rawPath_s_section = r'' + s_section_path
    if section_type == 'stud':
        try:
            with open(rawPath_s_section, 'r', encoding="utf-8-sig") as csv_file:
                # print('INSIDE CSV')
                csv_reader = csv.DictReader(csv_file)
                keys = csv_reader.fieldnames
                
                for line in csv_reader:
                    list_sections.append(line['ID'])
                
                if target_id not in list_sections:
                    return "ERROR: INVALID SECTION"
                else:
                    location = list_sections.index(target_id)
                    csv_file.seek(0)
                    next(csv_reader)
                    for count,line in enumerate(csv_reader):
                        if count == location:
                            Depth = float(line[keys[1]])
                            base = float(line[keys[2]])
                            t = float(line[keys[3]])
                            D = float(line[keys[4]])
                            R = float(line[keys[5]])
                            Area = float(line[keys[6]])
                            Ix = float(line[keys[8]])
                            Iy = float(line[keys[11]])
                            x_m = float(line[keys[14]])
                            m = float(line[keys[15]])
                            J = float(line[keys[16]])
                            Cw = float(line[keys[17]])
        except FileNotFoundError as error:
            return 'FILE NOT FOUND ' + rawPath_s_section
    else:
        return 'ERROR: Incorrect section used'

    if target_id in list_sections:
        Depth = Depth*inch
        Ag = Area*inch**2
        Ixxg = Ix*inch**4
        Iyyg = Iy*inch**4
        Cw = Cw*inch**6
        base = base*inch
        t = t*inch
        R = R*inch
        if section_type == 'stud':
            D = D*inch # this is the small d value - lip dimension
        x_m = x_m*inch
        J = J*inch**4
        r_y = sqrt(Iyyg/Ag)*inch
        m = m*inch
    
    # Built Up Section Properties:
    Atg = 2*Ag
    Ixxgt = 2*Ixxg
    Iyygt = 2*(Iyyg+Ag*((x_m)**2))
    Cwt = 2*Cw
    xot = 0*inch
    rxt = sqrt(Ixxgt/Atg)*inch
    ryt = sqrt(Iyygt/Atg)*inch
    rot = sqrt(rxt**2 + ryt**2 + xot**2)*inch
    Jt = 2*J
    
    # FLEXURAL BUCKLING STRESS
    
    x_slender = (K_x*L_x)/(rxt) # X-axis slenderness ratio
    r_i = r_y
    F_cre = (((pi)**2)*E)/(x_slender**2)
    y_slender = (K_y*L_y)/(ryt) # Y-axis slenderness ratio
    y_slendermod = sqrt((y_slender**2)+(a/r_i)**2)
    crit_slender = max(x_slender, y_slendermod)
    ratio_1 = (a/r_i)
    
    if ratio_1 <= (0.5*crit_slender):
        pass
    else: 
        status = 'NG'
        return f'Interconnection spacing is {status}'

    if x_slender > 200 and y_slender > 200:
        return 'Slenderness check failed, too slender' 
    
    F_cre_2 = (((pi)**2)*E)/(y_slendermod**2)
    F_cre_crit = min(F_cre, F_cre_2)
    
    # CHECK TORSIONAL FLEXURAL BUCKLING:
    sigma_t = (1/(Atg*(rot**2)))*(G*Jt+(((pi**2)*E*Cwt)/((K_t*L_t)**2)))
    F_cre_tb = sigma_t
    F_cre_tb = 5000000000*ksi # EXAGGERATED VALUE TO IGNORE THE EFFECTS OF TORSIONAL BUCKLING, MATCHING CFS DESIGNER
    # Fn calculation:
    F_cre_critf = min(F_cre_tb, F_cre_crit)
    lamb_c = sqrt(Fy/F_cre_critf)
    if lamb_c <= 1.5: Fn = (0.658**(lamb_c**2))*Fy
    else: Fn = (0.877/(lamb_c**2))*Fy
    
    # Effective Area at Fn:
    
    w = base - 2*(R+t)
    rel = w/t
    d1 = D-R-t
    f = Fn
    I_s = ((d1**3)*t*(sin(0.5*pi))**2)/12
    S = 1.28*sqrt(E/f)
    if rel <= (0.328*S): I_a = 0*inch**4
    else: I_a = min(399*(t**4)*(((w/t)/S)-0.328)**3,(t**4)*(115*((w/t)/S)+5))
    if rel <= (0.328*S): R_I = 0
    else: R_I = (I_s/I_a)
    if R_I > 1: R_I = 1
    relDw = (D/w)
    if relDw > 0.8: relDw = 0.8
    n = max((1/3), 0.582-((w/t)/(4*S)))
    
    if relDw > 0.25 and relDw <= 0.8: k = min((4.82-((5*D)/w))*((R_I)**n)+0.43,4)
    else: k = min(3.75*(R_I**n)+0.43,4)
    
    F_cr = k*(((pi**2)*E)/(12*(1-(poiss**2))))*((t/w)**2)
    lamb = sqrt(f/F_cr)
    rho = (1-(0.22/lamb))/lamb
    
    if lamb <= 0.673: b = w
    else: b = rho*w
    A_flanges = 2*(b*t)
    A_w_flange = 2*w*t
    A_dif_flange = A_w_flange - A_flanges
    
    k_lip = 0.43
    w_t = (d1/t)
    F_cr_2 = k_lip*(((pi**2)*E)/(12*(1-(poiss**2))))*((1/w_t)**2)
    lamb_2 = sqrt(f/F_cr_2)
    p_lip = (1-(0.22/lamb_2))/lamb_2
    d_s = d1
    if rel <= (0.328*S): d_lip = d_s
    else: d_lip = d_s*R_I
    # Two lips considered for the typical stud section
    A_lips = 2*(d_lip*t)
    A_w_lips = 2*d1*t
    A_w_lips_diff = A_w_lips-A_lips
    
    if P_cond == True:
        w_web = (Depth-2*(R+t)-P_width)/2
        k_web = 0.43
    else: 
        w_web = (Depth-2*(R+t))
        k_web = 4.0
        
    F_cr_web = k_web*(((pi**2)*E)/(12*(1-(poiss**2))))*((t/w_web)**2)
    lamb_web = sqrt(f/F_cr_web)
    rho_web = ((1-(0.22/lamb_web))/lamb_web)
    if lamb_web <= 0.673: b_web = w_web
    else: b_web = rho_web*w_web
    if P_cond == True: A_web_ineff = (2*(w_web-b_web))*t
    else: A_web_ineff = (w_web-b_web)*t
    
    if P_cond == True: A_effective = Ag - A_dif_flange - A_w_lips_diff - A_web_ineff - (P_width*t)
    else: A_effective = Ag - A_dif_flange - A_w_lips_diff - A_web_ineff
    
    # DISTORTIONAL BUCKLING CALCULATIONS:
    d_midline = D - (t/2)
    b_midline = base - t
    h_o = Depth
    A_web_gross = (Depth-2*(R+t))*t
    A_web_net = 2*w_web*t
    
    if P_cond == True: tr = t*((A_web_net/A_web_gross)**(1/3))
    else: tr = t
    A_net = 2*(Ag-(t*P_width))
    
    A_f = (b_midline+d_midline)*tr
    J_f = (1/3)*b_midline*(tr**3)+(1/3)*d_midline*(tr**3)
    I_xf = (tr*((tr**2)*(b_midline**2)+4*b_midline*(d_midline**3)+(tr**2)*b_midline*d_midline+(d_midline**4)))/(12*(b_midline+d_midline))
    I_yf = (tr*((b_midline**4)+4*d_midline*(b_midline**3)))/(12*(b_midline+d_midline))
    I_xyf = (tr*(b_midline**2)*(d_midline**2))/(4*(b_midline+d_midline))
    Cw_f = 0*inch**6
    h_xf = (-1*((b_midline**2)+2*b_midline*d_midline))/(2*(b_midline+d_midline))
    x_of = ((b_midline**2))/(2*(b_midline+d_midline))
    h_yf = (-1*(d_midline**2))/(2*(b_midline+d_midline))
    y_of = h_yf
    
    L_crd = (((6*(pi**4)*h_o*(1-(poiss**2)))/(tr**3))*(I_xf*((x_of-h_xf)**2)+Cw_f-((I_xyf**2)/I_yf)*((x_of-h_xf)**2)))**(1/4) # Critical unbraced length 
    
    L_m = S_bracing # Distance between discrete restraints that restrict buckling
    L_min = min(L_m,L_crd)
    K_phi_fe = (((pi)/L_min)**4)*(E*I_xf*((x_of-h_xf)**2)+E*Cw_f-E*((I_xyf**2)/I_yf)*((x_of-h_xf)**2))+(((pi)/L_min)**2)*G*J_f
    
    K_phi_we = (E*(tr**3))/(6*h_o*(1-(poiss**2)))
    K_phi = 0*kip
    K_phi_fg = ((pi/L_min)**2)*(A_f*(((x_of-h_xf)**2)*(((I_xyf)/(I_yf))**2)-2*y_of*(x_of-h_xf)*(I_xyf/I_yf)+(h_xf**2)+(y_of**2))+I_xf+I_yf)
    K_phi_fg = K_phi_fg*(1/(645.16*mm**2))*1*inch**2 # Ignore, unit conversion
    K_phi_wg = ((pi/L_min)**2)*((tr*(h_o**3))/(60))
    K_phi_wg = K_phi_wg*(1/(645.16*mm**2))*1*inch**2 # Ignore, unit conversion
    F_crd = (K_phi_fe+K_phi_we+K_phi)/(K_phi_fg+K_phi_wg)
    
    P_crd = Atg*F_crd
    P_y = Atg*Fy
    P_ynet = A_net*Fy
    lamb_d = sqrt(P_y/P_crd)
    lamb_d1 = 0.561*(P_ynet/P_y)
    lamb_d2 = 0.561*(14*((P_y/P_ynet)**0.4)-13)
    
    if lamb_d <= 0.561: P_ndNH = P_y
    else: P_ndNH = (1-0.25*((P_crd/P_y)**0.6))*((P_crd/P_y)**0.6)*P_y
    
    P_d2 = (1-0.25*((1/lamb_d2)**1.2))*((1/lamb_d2)**1.2)*P_y
    
    if lamb_d <= lamb_d2:
        if lamb_d <= lamb_d1: P_ndWH = P_ynet
        else: P_ndWH = P_ynet-((P_ynet-P_d2)/(lamb_d2-lamb_d1))*(lamb_d-lamb_d1)
    elif lamb_d > lamb_d2:
        if lamb_d <= 0.561:
            P_ndWH = P_y
        else: P_ndWH = (1-0.25*((P_crd/P_y)**0.6))*((P_crd/P_y)**0.6)*P_y
    
    if P_ndWH <= P_ndNH: pass
    else: return 'ERROR: P_ndWH > P_ndNH'
    
    P_ndfinal = P_ndWH if P_cond == True else P_ndNH
    
    # COMPUTE Pn:
    A_etotal = 2*A_effective
    P_nl = A_etotal*Fn
    P_n = min(P_nl,P_ndfinal)
    P_a = P_n/1.80    
    P_a_final = P_a.to('kip')
    
    return P_a_final

# BENDING STRENGTH FUNCTIONS:------------------------------------------------------------------------------------
def getFlexuralStrength_Single(target_id:str, section_type:str, Fy:float, F_bracing:float, E:float, G:float, P_length:float, P_width:float, P_cond:bool, K_y:float, K_t:float, L_x:float, Cb:float, dist_bracing=False,db_bracing=0.01):
    '''Calculation of CFS single member bending capacity.
    Calculations apply to stud and track sections. Calculations for punched tracks are not available, will be added in the future. 
    
    Args:
        target_id (str): Section shape name. 
        section_type (str): 'stud' or 'track'.
        Fy (float): Yield stress (ksi).
        F_bracing (float): Weak axis flexural bracing (ft).
        E (float): Elastic modulus (ksi).
        G (float): Shear Modulus (ksi).
        P_length (float): Punchout length, typically 4.0 inches (in)
        P_width (float): Punchout width, typically 1.5 inches (in).
        P_cond (bool): Punchout condition.
        K_y (float): Effective length factor y-axis.
        K_t (float): Effective length factor for twisting.
        L_x (float): Total unbraced length strong axis (ft).
        Cb (float): Bending coefficient.
        dist_bracing (bool): Distortional bracing status.
        db_bracing (float): Distortional bracing actual value in case dist_bracing is True (ft).
        
        It is recommended to keep the distortional bracing length equal or less than the value used for flexural weak axis bracing.
        Conservatively, bracing can be assumed to not exist.

        For torsional bracing, the same value used for weak axis bracing is automatically considered. 
    Returns:
        forallpeople.Physical: Allowable bending strength.
    '''
    poiss = 0.3     
    Fy = Fy*ksi
    F_bracing = F_bracing*ft
    E = E*ksi
    G = G*ksi   
    P_length = P_length*inch
    P_width = P_width*inch
    L_x = L_x*ft
    L_y = F_bracing
    L_t = F_bracing
    db_bracing = db_bracing*ft
    if dist_bracing == True:
        L_dist_br = db_bracing
    else:       
        L_dist_br = L_x

    if section_type == 'track':
        if P_cond == True:
            return 'ERROR: Punched tracks are not supported yet.'

    P_corner_radi = 0.25*inch # Punchout corner radii FIXED VALUE
    k_phi = 0*kip # rotational stiffness (0 to assume no stiffness provided - typically considered this way) FIXED VALUE

    s_section_path = os.path.dirname(__file__) + '\CFS_S_SECTION_DATA.csv'
    track_section_path = os.path.dirname(__file__) + '\CFS_TRACK_SECTION_DATA.csv'
    rawPath_s_section = r'' + s_section_path
    rawPath_track_section = r'' + track_section_path
    
    list_sections = []
    if section_type == 'stud':
        with open(rawPath_s_section, 'r', encoding="utf-8-sig") as csv_file:
            csv_reader = csv.DictReader(csv_file)
            keys = csv_reader.fieldnames
            
            for line in csv_reader:
                list_sections.append(line['ID'])
            
            if target_id not in list_sections:
                return 'ERROR: Section was not found in current reference file. Check your input or add the section manually.'
            else:
                location = list_sections.index(target_id)
                csv_file.seek(0)
                next(csv_reader)
                for count,line in enumerate(csv_reader):
                    if count == location:
                        Depth = float(line[keys[1]])
                        base = float(line[keys[2]])
                        t = float(line[keys[3]])
                        D = float(line[keys[4]])
                        R = float(line[keys[5]])
                        Area = float(line[keys[6]])
                        Ix = float(line[keys[8]])
                        Sx = float(line[keys[9]])
                        Iy = float(line[keys[11]])
                        x_m = float(line[keys[14]])
                        m = float(line[keys[15]])
                        J = float(line[keys[16]])
                        Cw = float(line[keys[17]])
                        xo = float(line[keys[20]])
    else:
        with open(rawPath_track_section, 'r', encoding="utf-8-sig") as csv_file:
            csv_reader = csv.DictReader(csv_file)
            keys = csv_reader.fieldnames
            for line in csv_reader:
                list_sections.append(line['ID'])
            if target_id not in list_sections:
                return 'ERROR: Section was not found in current reference file. Check your input or add the section manually.'
            else: 
                location = list_sections.index(target_id)
                csv_file.seek(0)
                next(csv_reader)
                for count,line in enumerate(csv_reader):
                    if count == location:
                        Depth = float(line[keys[1]])
                        base = float(line[keys[2]])
                        t = float(line[keys[3]])
                        R = float(line[keys[4]])
                        Area = float(line[keys[5]])
                        Ix = float(line[keys[7]])
                        Sx = float(line[keys[8]])
                        Iy = float(line[keys[10]])
                        x_m = float(line[keys[13]])
                        m = float(line[keys[14]])
                        J = float(line[keys[15]])
                        Cw = float(line[keys[16]])
                        xo = float(line[keys[19]])

    if target_id in list_sections:
        Depth = Depth*inch
        Ag = Area*inch**2
        Ixxg = Ix*inch**4
        Iyyg = Iy*inch**4
        x_o = xo*inch
        Cw = Cw*inch**6
        base = base*inch
        t = t*inch
        R = R*inch
        if section_type == 'stud':
            D = D*inch # this is the small d value - lip dimension
        x_m = x_m*inch
        J = J*inch**4
        r_x = sqrt(Ixxg/Ag)*inch
        r_y = sqrt(Iyyg/Ag)*inch
        r_o = sqrt((r_x**2)+(r_y**2)+(x_o**2))*inch
        m = m*inch
        S_x = Sx*inch**3
        S_f = S_x # Elastic section modulus of full unreduced section relative to extreme compression fiber
        S_fy = S_x # Elastic section modulus of full unreduced section relative to extreme fiber in first yielding
        r = R + t/2

    # revision of limits for shape with holes:
    procedure_alt = False
    if P_cond == True:
        d_h = P_width
        L_h = P_length
        h = Depth - 2*(R + t)
        clear_dist = 24*inch - P_length
        
        if d_h/h <= 0.7 and h/t <= 200 and clear_dist >= 18*inch and P_corner_radi>= 2*t and d_h<=2.5*inch and L_h <= 4.5*inch and d_h > (9/16)*inch:
            limitations_met = True # check if all limitations are met. If True, Fcre does not need to be modified
            if d_h/h >= 0.38:
                procedure_alt = True # procedure is altered
            else:
                procedure_alt = False # procedure remains the same, no alterations needed, hole ignored
        else:
            return 'ERROR: SECTION DOES NOT COMPLY WITH ALL LIMITATIONS, CANNOT CHECK' 

        # Average Properties:
        hole_area = P_width*t
        hole_length = P_width

        if section_type == 'stud':
            web_L = Depth - 2*(R+t)
            corner_length = (2*pi*(R + t/2))/4
            flange_length = base - 2*(R + t)
            lip_length = D - (R+t)
            Lg = web_L + 4*(corner_length)+2*(flange_length)+2*(lip_length) # segment lenght without holes
            L_net = Lg - hole_length
            A_net = Ag - hole_area
            L = Lg + L_net
            
            # NET INERTIA CALCULATOR STUD SECTION CALCULATED MANUALLY (NOT IN AISI MANUAL)
            I_net_xx = Ixxg - ((t*(P_width**3))/12)
            I_net_yy = Iyyg - (((P_width*(t**3))/12)+P_width*t*((x_m - t/2)**2))
            I_x_avg = (Ixxg*Lg + I_net_xx*L_net)/(L)
            I_y_avg = (Iyyg*Lg + I_net_yy*L_net)/(L)
            
            block_depth = ((Depth/2)-(P_width/2)-t)
            lip_block_depth = D - t
            
            x_m_net = (2*(block_depth)*(t)*(t/2)+2*(t)*(base)*(base/2)+2*lip_block_depth*t*(base-(t/2)))/(2*(block_depth)*(t)+2*t*base+2*t*lip_block_depth)
            
            A_avg = (Ag*Lg + A_net*L_net)/(L)
            x_o_net = (m-(t/2))+x_m_net
            y_o_net = 0*inch
            y_o = 0*inch
            x_o_avg = (x_o*Lg + x_o_net*L_net)/L
            y_o_avg = (y_o*Lg + y_o_net*L_net)/L
            ro_avg = sqrt((x_o_avg**2)+(y_o_avg**2)+((I_x_avg+I_y_avg)/A_avg))*inch
            
    # Initiation of Yielding Strength
    sigma_ey = ((pi**2)*E)/(((K_y*L_y)/(r_y))**2)
    sigma_t = (1/(Ag*r_o**2))*(G*J+(((pi**2)*E*Cw)/((K_t*L_t)**2)))
    if P_cond == True and limitations_met == False: F_cre = ((Cb*ro_avg*Ag)/(S_f))*sqrt(sigma_ey*sigma_t)*ksi
    else: F_cre = ((Cb*r_o*Ag)/(S_f))*sqrt(sigma_ey*sigma_t)*ksi
    
    if F_cre >= (2.78*Fy): F_n = Fy
    elif F_cre > (0.56*Fy) and F_cre < (2.78*Fy): F_n = (10/9)*Fy*(1-((10*Fy)/(36*F_cre)))
    elif F_cre <= (0.56*Fy): F_n = F_cre

    if P_cond == True: Sfy_net = I_net_xx/(Depth/2)
    if P_cond == True: M_y_net = Sfy_net*Fy
    M_y = S_fy*Fy
     
    # Cold work increases are not taken into consideration 
    if P_cond == True: M_ne = min(Sfy_net*F_n, M_y_net) # Punchout considered in section modulus
    else: M_ne = min(S_f*F_n, M_y)

    if P_cond == False:
        if section_type == 'track':
            S_e, ho, bo = calcEffectiveSectionModulus(Depth, F_n, procedure_alt, section_type, base, R, t, poiss, E, r, D=0, d_h=0,P_width=P_width, sectionType='single')   
            S_et, ho, bo = calcEffectiveSectionModulus(Depth, Fy, procedure_alt, section_type, base, R, t, poiss, E, r, D=0, d_h=0, P_width=P_width, sectionType='single')
        else:    
            S_e, ho, bo = calcEffectiveSectionModulus(Depth, F_n, procedure_alt, section_type, base, R, t, poiss, E, r, D, d_h=0,P_width=P_width, sectionType='single') # At max stress
            S_et, ho, bo = calcEffectiveSectionModulus(Depth, Fy, procedure_alt, section_type, base, R, t, poiss, E, r, D, d_h=0, P_width=P_width, sectionType='single') # At yielding point
    else: 
        S_e, ho, bo = calcEffectiveSectionModulus(Depth, F_n, procedure_alt, section_type, base, R, t, poiss, E, r, D, d_h, P_width, sectionType='single') # At max stress
        S_et, ho, bo = calcEffectiveSectionModulus(Depth, Fy, procedure_alt, section_type, base, R, t, poiss, E, r, D, d_h, P_width, sectionType='single') # At yielding point

    Mn_L = S_e*F_n
    Mn_L_Limit = S_et*Fy
    if Mn_L > Mn_L_Limit: Mn_L = Mn_L_Limit

    Mn_final = min(Mn_L, M_ne)
    M_allow = ((Mn_final/1.67)*(1000/12)*(1/(1*kip*inch)))*lb*ft

    # DISTORTIONAL BUCKLING: 
    '''
    USE OF GENERAL APPROACH ACCORDING TO THE COMMENTARY APPENDIX 2 SECTION 2.3.3.3(c) 
    NOTE: 
        - DISTORTIONAL BUCKLING ONLY APPLIES TO SECTIONS WITH LIP STIFFENERS - TRACKS ARE NOT CHECKED FOR THIS LIMIT STATE
        - MODIFIED THICKNESS IS APPLIED IN ORDER TO CONSIDER THE EFFECT OF TYPICAL PUNCHOUTS IN STUD SECTIONS. 
    '''
    if P_cond == True:
        # PUNCHOUTS CONSIDERED - USE OF MODIFIED THICKNESS IN EQ. 2.3.3.3-5 & 2.3.3.3-6 IN AISI SPECIFICATION VOL. 2
        if section_type == 'stud':
            h = ho - t
            b = bo - t
            do = D
            d = do - 0.5*t
            A_f = (b+d)*t
            I_xf = (t*((t**2)*(b**2)+4*b*(d**3)+(t**2)*b*d+(d**4)))/(12*(b+d))
            # I_xf = 5.4651*inch**4 #according to CFS designer, but does not match manual calc
            I_yf = (t*((b**4)+4*d*(b**3)))/(12*(b+d))
            I_xyf = (t*(b**2)*(d**2))/(4*(b+d))
            x_of = (b**2)/(2*(b+d))
            y_of = (-1*(d**2))/(2*(b+d))
            hxf = (-1*((b**2)+2*d*b))/(2*(b+d))
            J_f = (b*(t**3)+d*(t**3))/(3)
            C_wf = 0*inch**6

            L_crd = (((4*(pi**4)*ho*(1-(poiss**2)))/(t**3))*(I_xf*((x_of-hxf)**2)+C_wf-((I_xyf**2)/I_yf)*((x_of-hxf)**2))+(((pi**4)*(ho**4))/720))**0.25

            L = min(L_dist_br, L_crd)

            # TR CALCULATION METHOD ONE - ALTERNATIVE 
            # t_r = t*((1 - (L_h/L_crd))**(1/3))
            # TR CALCULATION METHOD TWO - ALTERNATIVE
            A_web_gross = L_x*(Depth - 2*(R+t))
            A_h = P_width*P_length # Hole area (estimate)
            typ_hole_spacing = 2 # Units in feet
            Lx_feet = L_x.to('ft')
            Lx_feet = Lx_feet*(1/(1*ft))
            Lx_feet = round(Lx_feet,4)
            if Lx_feet % typ_hole_spacing == 0:
                num_holes = Lx_feet/typ_hole_spacing
                num_holes -= 1
            else:
                num_holes = Lx_feet/typ_hole_spacing
                num_holes = num_holes.__floor__()
            A_web_net = A_web_gross - (num_holes*A_h)
            t_r = t*(A_web_net/A_web_gross)**(1/3)
            
            epsilon_web = 2
            
            k_phi_fe = ((pi/L)**4)*(E*I_xf*((x_of-hxf)**2)+E*C_wf-E*((I_xyf**2)/I_yf)*((x_of-hxf)**2))+((pi/L)**2)*G*J_f
            k_phi_we = ((E*(t_r**3))/(12*(1-(poiss**2))))*((3/ho)+((pi/L)**2)*((19*ho)/60)+((pi/L)**4)*((ho**3)/240))
            k_phi_fg = ((pi/L)**2)*(A_f*(((x_of-hxf)**2)*((I_xyf/I_yf)**2)-2*y_of*(x_of-hxf)*(I_xyf/I_yf)+(hxf**2)+(y_of**2))+I_xf+I_yf)
            k_phi_fg = k_phi_fg*(1/(1*mm**2))*0.00155*inch**2 # Manage the units in the factor for consistency in calculations
            k_phi_wg = ((ho*t_r*(pi**2))/(13440))*(((45360*(1-epsilon_web)+62160)*((L/ho)**2)+448*(pi**2)+((ho/L)**2)*(53+3*(1+epsilon_web))*(pi**4))/((pi**4)+28*(pi**2)*((L/ho)**2)+420*((L/ho)**4)))
            
            beta = 1.0
            F_crd = beta*((k_phi_fe+k_phi_we+k_phi)/(k_phi_fg + k_phi_wg))
            M_crd = S_f*F_crd # This value of moment is influenced by holes through the use of the modified thickness
            lamb_d = sqrt(M_y/M_crd) # typically the value of My is altered to include the effect of cold work forming increase for distortional buckling. 
            # In this case it is ignored to stick with the conservative side
            lamb_d1 = 0.673*(M_y_net/M_y)**3
            lamb_d2 = 0.673*(1.7*((M_y/M_y_net)**2.7)-0.7)
            
            if lamb_d <= lamb_d2:
                if lamb_d <= lamb_d1:
                    M_nd = M_y_net
                elif lamb_d1 < lamb_d <= lamb_d2:
                    M_d2 = (1-0.22*(1/lamb_d2))*(1/lamb_d2)*M_y
                    M_nd = M_y_net - ((M_y_net-M_d2)/(lamb_d2 - lamb_d1))*(lamb_d - lamb_d1)
                    M_nd_limit = (1-0.22*((M_crd/M_y)**0.5))*((M_crd/M_y)**0.5)*M_y
                    M_nd = min(M_nd, M_nd_limit)
            else:
                if lamb_d <= 0.673:
                    M_nd = M_y
                else:
                    M_nd = (1-0.22*((M_crd/M_y)**0.5))*((M_crd/M_y)**0.5)*M_y
            M_allow_db = M_nd/1.67
    else:
        # NO PUNCHOUTS CONSIDERED
        if section_type == 'stud':
            h = ho - t
            b = bo - t
            do = D
            d = do - 0.5*t
            A_f = (b+d)*t
            I_xf = (t*((t**2)*(b**2)+4*b*(d**3)+(t**2)*b*d+(d**4)))/(12*(b+d))
            I_yf = (t*((b**4)+4*d*(b**3)))/(12*(b+d))
            I_xyf = (t*(b**2)*(d**2))/(4*(b+d))
            x_of = (b**2)/(2*(b+d))
            y_of = (-1*(d**2))/(2*(b+d))
            hxf = (-1*((b**2)+2*d*b))/(2*(b+d))
            J_f = (b*(t**3)+d*(t**3))/(3)
            C_wf = 0*inch**6

            L_crd = (((4*(pi**4)*ho*(1-(poiss**2)))/(t**3))*(I_xf*((x_of-hxf)**2)+C_wf-((I_xyf**2)/I_yf)*((x_of-hxf)**2))+(((pi**4)*(ho**4))/720))**0.25

            L = min(L_dist_br, L_crd)
                
            epsilon_web = 2
            
            k_phi_fe = ((pi/L)**4)*(E*I_xf*((x_of-hxf)**2)+E*C_wf-E*((I_xyf**2)/I_yf)*((x_of-hxf)**2))+((pi/L)**2)*G*J_f
            k_phi_we = ((E*(t**3))/(12*(1-(poiss**2))))*((3/ho)+((pi/L)**2)*((19*ho)/60)+((pi/L)**4)*((ho**3)/240))
            k_phi_fg = ((pi/L)**2)*(A_f*(((x_of-hxf)**2)*((I_xyf/I_yf)**2)-2*y_of*(x_of-hxf)*(I_xyf/I_yf)+(hxf**2)+(y_of**2))+I_xf+I_yf)
            k_phi_wg = ((ho*t*(pi**2))/(13440))*(((45360*(1-epsilon_web)+62160)*((L/ho)**2)+448*(pi**2)+((ho/L)**2)*(53+3*(1+epsilon_web))*(pi**4))/((pi**4)+28*(pi**2)*((L/ho)**2)+420*((L/ho)**4)))
            
            beta = 1.0
            
            F_crd = beta*((k_phi_fe+k_phi_we+k_phi)/(k_phi_fg + k_phi_wg))
            
            M_crd = S_f*F_crd

            lamb_d = sqrt(M_y/M_crd) # typically the value of My is altered to include the effect of cold work forming increase for distortional buckling. 
            # In this case it is ignored to stick with the conservative side.
            if lamb_d <= 0.673:
                M_nd = M_y
            else:
                M_nd = (1-0.22*((M_crd/M_y)**0.5))*((M_crd/M_y)**0.5)*M_y
                
            M_allow_db = M_nd/1.67
            M_allow_db = (M_allow_db/(12*inch))*1*ft
        
    if section_type == 'stud': M_final_allow = min(M_allow_db, M_allow)
    else: M_final_allow = M_allow
    
    # final result correct units
    
    Final_Moment_Capacity = M_final_allow.to('kipft')    
    
    return Final_Moment_Capacity

def getFlexuralStrength_Boxed(target_id:str, section_type:str, L_stud:float, Fy:float, E:float, G:float, P_length:float, P_width:float, P_cond:bool, K_y:float, Cb:float, a:float):
    '''Calculation of CFS boxed member bending capacity.
    Calculations apply to boxed stud sections only. Other boxed sections are not available.  
    
    Args:
        target_id (str): Section shape name. 
        section_type (str): 'stud' or 'track'.
        L_stud (float): Total element length (ft).
        Fy (float): Yield stress (ksi).
        E (float): Elastic modulus (ksi).
        G (float): Shear Modulus (ksi).
        P_length (float): Punchout length, typically 4.0 inches (in)
        P_width (float): Punchout width, typically 1.5 inches (in).
        P_cond (bool): Punchout condition.
        K_y (float): Effective length factor y-axis.
        Cb (float): Bending coefficient.
        a (float): Interconnection spacing (in).
    
        No distortional buckling is considered for boxed section configurations. 
    Returns:
        forallpeople.Physical: Allowable bending strength.
    '''
    # Inputs: 
    target_id = target_id           
    section_type = section_type
    L_stud = L_stud*ft
    Fy = Fy*ksi
    P_length = P_length*inch
    P_width = P_width*inch
    P_corner_radi = 0.25*inch # Punchout corner radii
    F_bracing = 0.001*inch 
    E = E*ksi
    G = G*ksi
    poiss = 0.3
    L_y = F_bracing
    
    a = a*inch
    s_max = L_stud/6
    if a > s_max:
        return 'ERROR: Interconnection spacing is inadequate'

    s_section_path = os.path.dirname(__file__) + '\CFS_S_SECTION_DATA.csv'
    rawPath_s_section = r'' + s_section_path
    
    list_sections = []
    if section_type == 'stud':
        with open(rawPath_s_section, 'r', encoding="utf-8-sig") as csv_file:
            csv_reader = csv.DictReader(csv_file)
            keys = csv_reader.fieldnames
            
            for line in csv_reader:
                list_sections.append(line['ID'])
            
            if target_id not in list_sections:
                return 'ERROR: Section was not found in current reference file. Check your input or add the section manually.'
            else:
                location = list_sections.index(target_id)
                csv_file.seek(0)
                next(csv_reader)
                for count,line in enumerate(csv_reader):
                    if count == location:
                        Depth = float(line[keys[1]])
                        base = float(line[keys[2]])
                        t = float(line[keys[3]])
                        D = float(line[keys[4]])
                        R = float(line[keys[5]])
                        Area = float(line[keys[6]])
                        Ix = float(line[keys[8]])
                        Iy = float(line[keys[11]])
                        x_m = float(line[keys[14]])
                        m = float(line[keys[15]])
                        J = float(line[keys[16]])
                        Cw = float(line[keys[17]])
    else:
        return 'ERROR: Incorrect section used'
    
    if target_id in list_sections:
        Depth = Depth*inch
        Ag = Area*inch**2
        Ixxg = Ix*inch**4
        Iyyg = Iy*inch**4
        Cw = Cw*inch**6
        base = base*inch
        t = t*inch
        R = R*inch
        if section_type == 'stud':
            D = D*inch # this is the small d value - lip dimension
        x_m = x_m*inch
        J = J*inch**4
        m = m*inch
        r = R + t/2
    # Boxed section overall properties
    Ixxgt = 2*Ixxg
    Iyygt = 2*(Iyyg+Ag*((base - x_m)**2))
    S_ftot = (Ixxgt/(Depth/2))
    S_fy_tot = S_ftot
    # Interconnection spacing revision
    # revision of limits for shape with holes:
    procedure_alt = False
    if P_cond == True:
        d_h = P_width
        L_h = P_length
        h = Depth - 2*(R + t)
        clear_dist = 24*inch - P_length
        
        if d_h/h <= 0.7 and h/t <= 200 and clear_dist >= 18*inch and P_corner_radi>= 2*t and d_h<=2.5*inch and L_h <= 4.5*inch and d_h > (9/16)*inch:
            if d_h/h >= 0.38:
                procedure_alt = True # procedure is altered
            else:
                procedure_alt = False # procedure remains the same, no alterations needed, hole ignored
        else:
            return 'ERROR: SECTION DOES NOT COMPLY WITH ALL LIMITATIONS, CANNOT CHECK' 

        if section_type == 'stud':
            # NET INERTIA CALCULATOR STUD SECTION CALCULATED MANUALLY (NOT IN AISI MANUAL)
            I_net_xxt = Ixxgt - (((t*(P_width**3))/12)*2)

    # YIELDING AND GLOBAL LATERAL TORSIONAL BUCKLING:
    a_J = (base*2) - t
    b_J = (Depth) - t
    J_boxed = (2*(a_J*b_J)**2)/((a_J/t)+(b_J/t))
    L_u = ((0.36*Cb*pi)/(Fy*S_ftot))*(E*G*J_boxed*Iyygt)**(1/2)
    L_u = L_u.to('ft')
    
    if L_stud <= L_u:
        F_n = Fy
    else:
        F_cre = ((Cb*pi)/(K_y*L_y*S_fy_tot))*(E*G*J_boxed*Iyygt)**(1/2)
        if F_cre >= (2.78*Fy): F_n = Fy
        elif F_cre > (0.56*Fy) and F_cre < (2.78*Fy): F_n = (10/9)*Fy*(1-((10*Fy)/(36*F_cre)))
        elif F_cre <= (0.56*Fy): F_n = F_cre
    
    if P_cond == True: Sfy_net_t = I_net_xxt/(Depth/2)
    if P_cond == True: M_y_net = Sfy_net_t*Fy
    M_y = S_fy_tot*Fy 
    if P_cond == True: M_ne = min(Sfy_net_t*F_n, M_y_net)
    else: M_ne = min(S_fy_tot*F_n, M_y)
    
    # LOCAL BUCKLING INTERACTING WITH YIELDING AND GLOBAL BUCKLING:
    if procedure_alt == False:
        S_e = calculateEffectiveSectionModulus_BuiltUp(Depth, F_n, procedure_alt, section_type, base, R, t, poiss, E, r, D, d_h=0, P_width=P_width, AnalysisType='F_n', builtUpSectionType='boxed') # At max. stress
        S_et = calculateEffectiveSectionModulus_BuiltUp(Depth, Fy, procedure_alt, section_type, base, R, t, poiss, E, r, D, d_h=0, P_width=P_width, AnalysisType='Fy', builtUpSectionType='boxed')
    else:
        S_e = calculateEffectiveSectionModulus_BuiltUp(Depth, F_n, procedure_alt, section_type, base, R, t, poiss, E, r, D, d_h, P_width, AnalysisType='F_n', builtUpSectionType='boxed') # At max. stress
        S_et = calculateEffectiveSectionModulus_BuiltUp(Depth, Fy, procedure_alt, section_type, base, R, t, poiss, E, r, D, d_h, P_width, AnalysisType='Fy', builtUpSectionType='boxed')
    
    Mn_L = S_e*F_n
    Mn_L_Limit = S_et*Fy
    if Mn_L > Mn_L_Limit: Mn_L = Mn_L_Limit
    Mn_final = min(Mn_L, M_ne)
    M_allow = Mn_final/1.67
    M_allow_final = M_allow.to('kipft')

    return M_allow_final

def getFlexuralStrength_B2B(target_id:str, section_type:str, L_stud:float, Fy:float, F_bracing:float, E:float, G:float, P_length:float, P_width:float, P_cond:bool, K_y:float, K_t:float, L_x:float, Cb:float, a, dist_bracing, db_bracing=0.01):
    '''Calculation of CFS back to back member bending capacity.
    Calculations apply to back to back stud sections only. Other back to back sections are not available.  
    
    Args:
        target_id (str): Section shape name. 
        section_type (str): 'stud' or 'track'.
        L_stud (float): Total element length for interconnection spacing revision (ft).
        Fy (float): Yield stress (ksi).
        F_bracing (float): Weak axis flexural bracing distance applicable for torsion as well (ft).
        E (float): Elastic modulus (ksi).
        G (float): Shear Modulus (ksi).
        P_length (float): Punchout length, typically 4.0 inches (in)
        P_width (float): Punchout width, typically 1.5 inches (in).
        P_cond (bool): Punchout condition.
        K_y (float): Effective length factor y-axis.
        K_t (float): Effective length factor for twisting.
        L_x (float): Strong axis unbraced length (ft).
        Cb (float): Bending coefficient.
        a (float): Interconnection spacing (in).
        dist_bracing (bool): Distortional bracing status. 
        db_bracing (float): Distortional bracing value if True (ft).
     
    Returns:
        forallpeople.Physical: Allowable bending strength.
    '''
    # Inputs: 
    target_id = target_id
    section_type = section_type
    L_stud = L_stud*ft
    Fy = Fy*ksi
    P_length = P_length*inch
    P_width = P_width*inch
    P_corner_radi = 0.25*inch # Punchout corner radii
    F_bracing = F_bracing*ft # Flexural Bracing spacing (For fully braced condition, use a small number close to 0)
    E = E*ksi
    G = G*ksi
    poiss = 0.3
    L_x = L_x*ft
    L_t = F_bracing
    L_y = F_bracing
    db_bracing = db_bracing*ft
    if dist_bracing == True:
        L_dist_br = db_bracing
    else:       
        L_dist_br = L_x
    
    a = a*inch
    s_max = L_stud/6
    if a > s_max:
        return 'ERROR: Interconnection spacing is inadequate'

    k_phi = 0*kip # rotational stiffness (0 to assume no stiffness provided - typically considered this way)

    s_section_path = os.path.dirname(__file__) + '\CFS_S_SECTION_DATA.csv'
    rawPath_s_section = r'' + s_section_path
    
    list_sections = []
    if section_type == 'stud':
        with open(rawPath_s_section, 'r', encoding="utf-8-sig") as csv_file:
            csv_reader = csv.DictReader(csv_file)
            keys = csv_reader.fieldnames
            
            for line in csv_reader:
                list_sections.append(line['ID'])
            
            if target_id not in list_sections:
                return 'ERROR: Section was not found in current reference file. Check your input or add the section manually.'
            else:
                location = list_sections.index(target_id)
                csv_file.seek(0)
                next(csv_reader)
                for count,line in enumerate(csv_reader):
                    if count == location:
                        Depth = float(line[keys[1]])
                        base = float(line[keys[2]])
                        t = float(line[keys[3]])
                        D = float(line[keys[4]])
                        R = float(line[keys[5]])
                        Area = float(line[keys[6]])
                        Ix = float(line[keys[8]])
                        Sx = float(line[keys[9]])
                        Iy = float(line[keys[11]])
                        x_m = float(line[keys[14]])
                        m = float(line[keys[15]])
                        J = float(line[keys[16]])
                        Cw = float(line[keys[17]])
                        xo = float(line[keys[20]])
    else:
        return 'ERROR: Incorrect section used'

    if target_id in list_sections:
        Depth = Depth*inch
        Ag = Area*inch**2
        Ixxg = Ix*inch**4
        Iyyg = Iy*inch**4
        x_o = xo*inch
        Cw = Cw*inch**6
        base = base*inch
        t = t*inch
        R = R*inch
        if section_type == 'stud':
            D = D*inch # this is the small d value - lip dimension
        x_m = x_m*inch
        J = J*inch**4
        m = m*inch
        S_x = Sx*inch**3
        S_f = S_x # Elastic section modulus of full unreduced section relative to extreme compression fibers
        r = R + t/2
    
    # BUILT UP SECTION PROPERTIES:
    Atg = 2*Ag
    Ixxgt = 2*Ixxg
    Iyygt = 2*(Iyyg+Ag*((x_m)**2))
    Cwt = 2*Cw
    xot = 0*inch
    rxt = sqrt(Ixxgt/Atg)*inch
    ryt = sqrt(Iyygt/Atg)*inch
    rot = sqrt(rxt**2 + ryt**2 + xot**2)*inch
    J_tot = J*2
    S_ftot = (Ixxgt/(Depth/2))
    S_fy_tot = S_ftot
    
    # revision of limits for shape with holes -  single section only:
    procedure_alt = False
    if P_cond == True:
        d_h = P_width
        L_h = P_length
        h = Depth - 2*(R + t)
        clear_dist = 24*inch - P_length
        
        if d_h/h <= 0.7 and h/t <= 200 and clear_dist >= 18*inch and P_corner_radi>= 2*t and d_h<=2.5*inch and L_h <= 4.5*inch and d_h > (9/16)*inch:
            limitations_met = True # check if all limitations are met. If True, Fcre does not need to be modified
            if d_h/h >= 0.38:
                procedure_alt = True # procedure is altered
            else:
                procedure_alt = False # procedure remains the same, no alterations needed, hole ignored
        else:
            return 'ERROR: SECTION DOES NOT COMPLY WITH ALL LIMITATIONS, CANNOT CHECK' 

        # Average Properties:
        hole_area = P_width*t
        hole_length = P_width

        if section_type == 'stud':
            web_L = Depth - 2*(R+t)
            corner_length = (2*pi*(R + t/2))/4
            flange_length = base - 2*(R + t)
            lip_length = D - (R+t)
            Lg = web_L + 4*(corner_length)+2*(flange_length)+2*(lip_length) # segment lenght without holes
            L_net = Lg - hole_length
            A_net = Ag - hole_area
            L = Lg + L_net
            
            # NET INERTIA CALCULATOR STUD SECTION CALCULATED MANUALLY (NOT IN AISI MANUAL)
            I_net_xx = Ixxg - ((t*(P_width**3))/12)
            I_net_yy = Iyyg - (((P_width*(t**3))/12)+P_width*t*((x_m - t/2)**2))
            I_x_avg = (Ixxg*Lg + I_net_xx*L_net)/(L)
            I_y_avg = (Iyyg*Lg + I_net_yy*L_net)/(L)
            block_depth = ((Depth/2)-(P_width/2)-t)
            lip_block_depth = D - t
            x_m_net = (2*(block_depth)*(t)*(t/2)+2*(t)*(base)*(base/2)+2*lip_block_depth*t*(base-(t/2)))/(2*(block_depth)*(t)+2*t*base+2*t*lip_block_depth)
            A_avg = (Ag*Lg + A_net*L_net)/(L)
            x_o_net = (m-(t/2))+x_m_net
            y_o_net = 0*inch
            y_o = 0*inch
            x_o_avg = (x_o*Lg + x_o_net*L_net)/L
            y_o_avg = (y_o*Lg + y_o_net*L_net)/L
            ro_avg = sqrt((x_o_avg**2)+(y_o_avg**2)+((I_x_avg+I_y_avg)/A_avg))*inch

    # YIELDING AND GLOBAL (LATERAL TORSIONAL) BUCKLINGL:
    # Initiation of Yielding Strength
    sigma_ey = ((pi**2)*E)/(((K_y*L_y)/(ryt))**2)
    sigma_t = (1/(Atg*rot**2))*(G*J_tot+(((pi**2)*E*Cwt)/((K_t*L_t)**2)))
    if P_cond == True and limitations_met == False: F_cre = ((Cb*ro_avg*Atg)/(S_ftot))*(sigma_ey*sigma_t)**(1/2)
    else: F_cre = ((Cb*rot*Atg)/(S_ftot))*(sigma_ey*sigma_t)**(1/2)
    
    if F_cre >= (2.78*Fy): F_n = Fy
    elif F_cre > (0.56*Fy) and F_cre < (2.78*Fy): F_n = (10/9)*Fy*(1-((10*Fy)/(36*F_cre)))
    elif F_cre <= (0.56*Fy): F_n = F_cre
    
    if P_cond == True: Sfy_net = (I_net_xx*2)/(Depth/2)
    if P_cond == True: M_y_net = Sfy_net*Fy
    M_y = S_fy_tot*Fy 
    if P_cond == True: M_ne = min(Sfy_net*F_n, M_y_net)
    else: M_ne = min(S_ftot*F_n, M_y)
    
    # LOCAL BUCKLING INTERACTING WITH YIELDING AND GLOBAL BUCKLING:
    # Effective Section Modulus Calculations:
    if P_cond == False:    
        S_e, ho, bo = calculateEffectiveSectionModulus_BuiltUp(Depth, F_n, procedure_alt, section_type, base, R, t, poiss, E, r, D, d_h=0, P_width=P_width, AnalysisType='F_n', builtUpSectionType='b2b')
        S_et = calculateEffectiveSectionModulus_BuiltUp(Depth, Fy, procedure_alt, section_type, base, R, t, poiss, E, r, D, d_h=0, P_width=P_width, AnalysisType='Yielding', builtUpSectionType='b2b')
    else:
        S_e, ho, bo  = calculateEffectiveSectionModulus_BuiltUp(Depth, F_n, procedure_alt, section_type, base, R, t, poiss, E, r, D, d_h, P_width,'F_n','b2b')
        S_et = calculateEffectiveSectionModulus_BuiltUp(Depth, Fy, procedure_alt, section_type, base, R, t, poiss, E, r, D, d_h, P_width, 'Yielding','b2b')
    
    # MOMENT CAPACITY CALCULATION: 
    Mn_L = S_e*F_n
    Mn_L_Limit = S_et*Fy
    if Mn_L > Mn_L_Limit: Mn_L = Mn_L_Limit
    Mn_L = Mn_L.to('kipft')
    
    # DISTORTIONAL BUCKLING:
    M_y_dist = S_f*Fy # Yielding Moment for a single section
    if P_cond == True:
        # PUNCHOUTS CONSIDERED - USE OF MODIFIED THICKNESS IN EQ. 2.3.3.3-5 & 2.3.3.3-6 IN AISI SPECIFICATION VOL. 2
        if section_type == 'stud':
            h = ho - t
            b = bo - t
            do = D
            d = do - 0.5*t
            A_f = (b+d)*t
            I_xf = (t*((t**2)*(b**2)+4*b*(d**3)+(t**2)*b*d+(d**4)))/(12*(b+d))
            I_yf = (t*((b**4)+4*d*(b**3)))/(12*(b+d))
            I_xyf = (t*(b**2)*(d**2))/(4*(b+d))
            x_of = (b**2)/(2*(b+d))
            y_of = (-1*(d**2))/(2*(b+d))
            hxf = (-1*((b**2)+2*d*b))/(2*(b+d))
            J_f = (b*(t**3)+d*(t**3))/(3)
            C_wf = 0*inch**6
            
            L_crd = (((4*(pi**4)*ho*(1-(poiss**2)))/(t**3))*(I_xf*((x_of-hxf)**2)+C_wf-((I_xyf**2)/I_yf)*((x_of-hxf)**2))+(((pi**4)*(ho**4))/720))**0.25

            L = min(L_dist_br, L_crd)
            
            # TR CALCULATION METHOD 3
            t_r = t*((1 - (L_h/L_crd))**(1/3))
            
            #region
            # TR CALCULATION ALTERNATIVE
            # A_web_gross = L_x*(Depth - 2*(R+t))
            # A_h = P_width*P_length # Hole area (estimate)
            # typ_hole_spacing = 2 # Units in feet
            # Lx_feet = L_x.to('ft')
            # Lx_feet = Lx_feet*(1/(1*ft))
            # Lx_feet = round(Lx_feet,4)
            # if Lx_feet % typ_hole_spacing == 0:
            #     num_holes = Lx_feet/typ_hole_spacing
            #     num_holes -= 1
            # else:
            #     num_holes = Lx_feet/typ_hole_spacing
            #     num_holes = num_holes.__floor__()
            # A_web_net = A_web_gross - (num_holes*A_h)
            # t_r = t*(A_web_net/A_web_gross)**(1/3)
            #endregion
            
            epsilon_web = 2
            
            k_phi_fe = ((pi/L)**4)*(E*I_xf*((x_of-hxf)**2)+E*C_wf-E*((I_xyf**2)/I_yf)*((x_of-hxf)**2))+((pi/L)**2)*G*J_f
            k_phi_we = ((E*(t_r**3))/(12*(1-(poiss**2))))*((3/ho)+((pi/L)**2)*((19*ho)/60)+((pi/L)**4)*((ho**3)/240))
            k_phi_fg = ((pi/L)**2)*(A_f*(((x_of-hxf)**2)*((I_xyf/I_yf)**2)-2*y_of*(x_of-hxf)*(I_xyf/I_yf)+(hxf**2)+(y_of**2))+I_xf+I_yf)
            k_phi_fg = k_phi_fg*(1/(1*mm**2))*0.00155*inch**2 # Manage the units in the factor for consistency in calculations
            k_phi_wg = ((ho*t_r*(pi**2))/(13440))*(((45360*(1-epsilon_web)+62160)*((L/ho)**2)+448*(pi**2)+((ho/L)**2)*(53+3*(1+epsilon_web))*(pi**4))/((pi**4)+28*(pi**2)*((L/ho)**2)+420*((L/ho)**4)))
            
            beta = 1.0
            F_crd = beta*((k_phi_fe+k_phi_we+k_phi)/(k_phi_fg + k_phi_wg))
            
            M_crd = S_f*F_crd # This value of moment is influenced by holes through the use of the modified thickness
                
            lamb_d = sqrt(M_y_dist/M_crd) # typically the value of My is altered to include the effect of cold work forming increase for distortional buckling. 
            # In this case it is ignored to stick with the conservative side
            lamb_d1 = 0.673*(M_y_net/M_y_dist)**3
            lamb_d2 = 0.673*(1.7*((M_y_dist/M_y_net)**2.7)-0.7)
            
            if lamb_d <= lamb_d2:
                if lamb_d <= lamb_d1:
                    M_nd = M_y_net
                elif lamb_d1 < lamb_d <= lamb_d2:
                    M_d2 = (1-0.22*(1/lamb_d2))*(1/lamb_d2)*M_y
                    M_nd = M_y_net - ((M_y_net-M_d2)/(lamb_d2 - lamb_d1))*(lamb_d - lamb_d1)
                    M_nd_limit = (1-0.22*((M_crd/M_y_dist)**0.5))*((M_crd/M_y_dist)**0.5)*M_y_dist
                    M_nd = min(M_nd, M_nd_limit)
            else:
                if lamb_d <= 0.673:
                    M_nd = M_y_dist
                else:
                    M_nd = (1-0.22*((M_crd/M_y_dist)**0.5))*((M_crd/M_y_dist)**0.5)*M_y_dist
            M_allow_db = M_nd/1.67
            
        else:
            return 'DISTORTIONAL BUCKLING NOT REVISED, TRACK SECTION USED' 
    else:
        # NO PUNCHOUTS CONSIDERED
        if section_type == 'stud':
            h = ho - t
            b = bo - t
            do = D
            d = do - 0.5*t
            A_f = (b+d)*t
            I_xf = (t*((t**2)*(b**2)+4*b*(d**3)+(t**2)*b*d+(d**4)))/(12*(b+d))
            I_yf = (t*((b**4)+4*d*(b**3)))/(12*(b+d))
            I_xyf = (t*(b**2)*(d**2))/(4*(b+d))
            x_of = (b**2)/(2*(b+d))
            y_of = (-1*(d**2))/(2*(b+d))
            hxf = (-1*((b**2)+2*d*b))/(2*(b+d))
            J_f = (b*(t**3)+d*(t**3))/(3)
            C_wf = 0*inch**6

            L_crd = (((4*(pi**4)*ho*(1-(poiss**2)))/(t**3))*(I_xf*((x_of-hxf)**2)+C_wf-((I_xyf**2)/I_yf)*((x_of-hxf)**2))+(((pi**4)*(ho**4))/720))**0.25
            L = min(L_dist_br, L_crd)
            epsilon_web = 2
            
            k_phi_fe = ((pi/L)**4)*(E*I_xf*((x_of-hxf)**2)+E*C_wf-E*((I_xyf**2)/I_yf)*((x_of-hxf)**2))+((pi/L)**2)*G*J_f
            k_phi_we = ((E*(t**3))/(12*(1-(poiss**2))))*((3/ho)+((pi/L)**2)*((19*ho)/60)+((pi/L)**4)*((ho**3)/240))
            k_phi_fg = ((pi/L)**2)*(A_f*(((x_of-hxf)**2)*((I_xyf/I_yf)**2)-2*y_of*(x_of-hxf)*(I_xyf/I_yf)+(hxf**2)+(y_of**2))+I_xf+I_yf)
            k_phi_wg = ((ho*t*(pi**2))/(13440))*(((45360*(1-epsilon_web)+62160)*((L/ho)**2)+448*(pi**2)+((ho/L)**2)*(53+3*(1+epsilon_web))*(pi**4))/((pi**4)+28*(pi**2)*((L/ho)**2)+420*((L/ho)**4)))
            beta = 1.0
            
            F_crd = beta*((k_phi_fe+k_phi_we+k_phi)/(k_phi_fg + k_phi_wg))
            M_crd = S_f*F_crd
            lamb_d = sqrt(M_y_dist/M_crd) # typically the value of My is altered to include the effect of cold work forming increase for distortional buckling. 
            # In this case it is ignored to stick with the conservative side
            
            if lamb_d <= 0.673:
                M_nd = M_y_dist
            else:
                M_nd = (1-0.22*((M_crd/M_y_dist)**0.5))*((M_crd/M_y_dist)**0.5)*M_y_dist
                
            M_allow_db = M_nd/1.67
            M_allow_db = (M_allow_db/(12*inch))*1*ft
            
        else:
            return 'ERROR: DISTORTIONAL BUCKLING NOT REVISED, TRACK OR UNSUPPORTED SECTION USED'

    M_allow_db_total = M_allow_db*2    

    # FINAL FLEXURAL CAPACITY: 
    Mn_allow_ltb = Mn_L/1.67
    Mn_allow_ltb = Mn_allow_ltb.to('kipft')
    Mn_allow_gtb = M_ne/1.67
    Mn_allow_gtb = Mn_allow_gtb.to('kipft')
    M_allow_final = min(M_allow_db_total, Mn_allow_ltb, Mn_allow_gtb)
    M_allow_final = M_allow_final.to('kipft')

    return M_allow_final

def getFlexuralStrength_Single_GB(target_id:str, section_type:str, Fy:float, F_bracing:float, E:float, G:float, P_length:float, P_width:float, P_cond:bool, K_y:float, K_t:float, L_x:float, Cb:float):
    '''Calculation of globally braced CFS single member bending capacity.
    Calculations apply to stud and track sections. Calculations for punched tracks will be added in the future.
    
    Args:
        target_id (str): Section shape name. 
        section_type (str): 'stud' or 'track'.
        Fy (float): Yield stress (ksi).
        F_bracing (float): Weak axis flexural bracing (ft).
        E (float): Elastic modulus (ksi).
        G (float): Shear Modulus (ksi).
        P_length (float): Punchout length, typically 4.0 inches (in)
        P_width (float): Punchout width, typically 1.5 inches (in).
        P_cond (bool): Punchout condition.
        K_y (float): Effective length factor y-axis.
        K_t (float): Effective length factor for twisting.
        L_x (float): Total unbraced length strong axis (ft).
        Cb (float): Bending coefficient.
        
        It is recommended to keep the distortional bracing length equal or less than the value used for flexural weak axis bracing.
        Conservatively, bracing can be assumed to not exist.

        For torsional bracing, the same value used for weak axis bracing is automatically considered. 
    Returns:
        forallpeople.Physical: Allowable bending strength of globally braced single stud/track section.
    '''
    poiss = 0.3
    Fy = Fy*ksi
    F_bracing = F_bracing*ft
    E = E*ksi
    G = G*ksi
    P_length = P_length*inch
    P_width = P_width*inch
    L_x = L_x*ft
    L_y = F_bracing
    L_t = F_bracing
    
    if section_type == 'track':
        if P_cond == True:
            return 'ERROR: Punched tracks are not supported yet.'

    P_corner_radi = 0.25*inch # Punchout corner radii FIXED VALUE

    s_section_path = os.path.dirname(__file__) + '\CFS_S_SECTION_DATA.csv'
    track_section_path = os.path.dirname(__file__) + '\CFS_TRACK_SECTION_DATA.csv'
    rawPath_s_section = r'' + s_section_path
    rawPath_track_section = r'' + track_section_path
    
    list_sections = []
    if section_type == 'stud':
        with open(rawPath_s_section, 'r', encoding="utf-8-sig") as csv_file:
            csv_reader = csv.DictReader(csv_file)
            keys = csv_reader.fieldnames
            
            for line in csv_reader:
                list_sections.append(line['ID'])
            
            if target_id not in list_sections:
                return 'ERROR: Section was not found in current reference file. Check your input or add the section manually.'
            else:
                location = list_sections.index(target_id)
                csv_file.seek(0)
                next(csv_reader)
                for count,line in enumerate(csv_reader):
                    if count == location:
                        Depth = float(line[keys[1]])
                        base = float(line[keys[2]])
                        t = float(line[keys[3]])
                        D = float(line[keys[4]])
                        R = float(line[keys[5]])
                        Area = float(line[keys[6]])
                        Ix = float(line[keys[8]])
                        Sx = float(line[keys[9]])
                        Iy = float(line[keys[11]])
                        x_m = float(line[keys[14]])
                        m = float(line[keys[15]])
                        J = float(line[keys[16]])
                        Cw = float(line[keys[17]])
                        xo = float(line[keys[20]])
    else:
        with open(rawPath_track_section, 'r', encoding="utf-8-sig") as csv_file:
            csv_reader = csv.DictReader(csv_file)
            keys = csv_reader.fieldnames
            
            for line in csv_reader:
                list_sections.append(line['ID'])
            if target_id not in list_sections:
                return 'ERROR: Section was not found in current reference file. Check your input or add the section manually.'
            else: 
                location = list_sections.index(target_id)
                csv_file.seek(0)
                next(csv_reader)
                for count,line in enumerate(csv_reader):
                    if count == location:
                        Depth = float(line[keys[1]])
                        base = float(line[keys[2]])
                        t = float(line[keys[3]])
                        R = float(line[keys[4]])
                        Area = float(line[keys[5]])
                        Ix = float(line[keys[7]])
                        Sx = float(line[keys[8]])
                        Iy = float(line[keys[10]])
                        x_m = float(line[keys[13]])
                        m = float(line[keys[14]])
                        J = float(line[keys[15]])
                        Cw = float(line[keys[16]])
                        xo = float(line[keys[19]])

    if target_id in list_sections:
        Depth = Depth*inch
        Ag = Area*inch**2
        Ixxg = Ix*inch**4
        Iyyg = Iy*inch**4
        x_o = xo*inch
        Cw = Cw*inch**6
        base = base*inch
        t = t*inch
        R = R*inch
        if section_type == 'stud':
            D = D*inch # this is the small d value - lip dimension
        x_m = x_m*inch
        J = J*inch**4
        r_x = sqrt(Ixxg/Ag)*inch
        r_y = sqrt(Iyyg/Ag)*inch
        r_o = sqrt((r_x**2)+(r_y**2)+(x_o**2))*inch
        m = m*inch
        S_x = Sx*inch**3
        S_f = S_x # Elastic section modulus of full unreduced section relative to extreme compression fiber
        r = R + t/2

    # revision of limits for shape with holes:
    procedure_alt = False
    if P_cond == True:
        d_h = P_width
        L_h = P_length
        h = Depth - 2*(R + t)
        clear_dist = 24*inch - P_length
        
        if d_h/h <= 0.7 and h/t <= 200 and clear_dist >= 18*inch and P_corner_radi>= 2*t and d_h<=2.5*inch and L_h <= 4.5*inch and d_h > (9/16)*inch:
            limitations_met = True # check if all limitations are met. If True, Fcre does not need to be modified
            if d_h/h >= 0.38:
                procedure_alt = True # procedure is altered
            else:
                procedure_alt = False # procedure remains the same, no alterations needed, hole ignored
        else:
            return 'ERROR: SECTION DOES NOT COMPLY WITH ALL LIMITATIONS, CANNOT CHECK' 

        # Average Properties:

        hole_area = P_width*t
        hole_length = P_width

        if section_type == 'stud':
            web_L = Depth - 2*(R+t)
            corner_length = (2*pi*(R + t/2))/4
            flange_length = base - 2*(R + t)
            lip_length = D - (R+t)
            Lg = web_L + 4*(corner_length)+2*(flange_length)+2*(lip_length) # segment lenght without holes
            L_net = Lg - hole_length
            A_net = Ag - hole_area
            L = Lg + L_net
            
            # NET INERTIA CALCULATOR STUD SECTION CALCULATED MANUALLY (NOT IN AISI MANUAL)
            I_net_xx = Ixxg - ((t*(P_width**3))/12)
            I_net_yy = Iyyg - (((P_width*(t**3))/12)+P_width*t*((x_m - t/2)**2))
            I_x_avg = (Ixxg*Lg + I_net_xx*L_net)/(L)
            I_y_avg = (Iyyg*Lg + I_net_yy*L_net)/(L)
            
            block_depth = ((Depth/2)-(P_width/2)-t)
            lip_block_depth = D - t
            
            x_m_net = (2*(block_depth)*(t)*(t/2)+2*(t)*(base)*(base/2)+2*lip_block_depth*t*(base-(t/2)))/(2*(block_depth)*(t)+2*t*base+2*t*lip_block_depth)
            
            A_avg = (Ag*Lg + A_net*L_net)/(L)
            x_o_net = (m-(t/2))+x_m_net
            y_o_net = 0*inch
            y_o = 0*inch
            x_o_avg = (x_o*Lg + x_o_net*L_net)/L
            y_o_avg = (y_o*Lg + y_o_net*L_net)/L
            ro_avg = sqrt((x_o_avg**2)+(y_o_avg**2)+((I_x_avg+I_y_avg)/A_avg))*inch
    
    # Initiation of Yielding Strength
    sigma_ey = ((pi**2)*E)/(((K_y*L_y)/(r_y))**2)
    sigma_t = (1/(Ag*r_o**2))*(G*J+(((pi**2)*E*Cw)/((K_t*L_t)**2)))
    if P_cond == True and limitations_met == False: F_cre = ((Cb*ro_avg*Ag)/(S_f))*sqrt(sigma_ey*sigma_t)*ksi
    else: F_cre = ((Cb*r_o*Ag)/(S_f))*sqrt(sigma_ey*sigma_t)*ksi

    if F_cre >= (2.78*Fy): F_n = Fy
    elif F_cre > (0.56*Fy) and F_cre < (2.78*Fy): F_n = (10/9)*Fy*(1-((10*Fy)/(36*F_cre)))
    elif F_cre <= (0.56*Fy): F_n = F_cre

    sectionType = 'other'
    if P_cond == False:    
        S_e = calcEffectiveSectionModulus(Depth, F_n, procedure_alt, section_type, base, R, t, poiss, E, r, D, d_h=0, P_width=P_width, sectionType=sectionType)
        S_et = calcEffectiveSectionModulus(Depth, Fy, procedure_alt, section_type, base, R, t, poiss, E, r, D, d_h=0, P_width=P_width, sectionType=sectionType)
    else: 
        S_e = calcEffectiveSectionModulus(Depth, F_n, procedure_alt, section_type, base, R, t, poiss, E, r, D, d_h, P_width, sectionType)
        S_et = calcEffectiveSectionModulus(Depth, Fy, procedure_alt, section_type, base, R, t, poiss, E, r, D, d_h, P_width, sectionType)
    M_alo = min(S_e*Fy, S_et*Fy)
    M_alo = M_alo/1.67
    M_alo_final = M_alo.to('kipft')
    return M_alo_final 

def getFlexuralStrength_Boxed_GB(target_id:str, section_type:str, L_stud:float, Fy:float, E:float, G:float, P_length:float, P_width:float, P_cond:bool, K_y:float, Cb:float, a:float):
    '''Calculation of globally braced CFS boxed member bending capacity.
    Calculations apply to boxed stud sections only. Other boxed sections are not available.  
    
    Args:
        target_id (str): Section shape name. 
        section_type (str): 'stud' or 'track'.
        L_stud (float): Total element length (ft).
        Fy (float): Yield stress (ksi).
        E (float): Elastic modulus (ksi).
        G (float): Shear Modulus (ksi).
        P_length (float): Punchout length, typically 4.0 inches (in)
        P_width (float): Punchout width, typically 1.5 inches (in).
        P_cond (bool): Punchout condition.
        K_y (float): Effective length factor y-axis.
        Cb (float): Bending coefficient.
        a (float): Interconnection spacing (in).
    
        No distortional buckling is considered for boxed section configurations. 
    Returns:
        forallpeople.Physical: Allowable bending strength of globally braced boxed stud section.
    '''
    target_id = target_id
    section_type = section_type
    L_stud = L_stud*ft
    Fy = Fy*ksi
    P_length = P_length*inch
    P_width = P_width*inch
    P_corner_radi = 0.25*inch # Punchout corner radii
    F_bracing = 0.001*inch # Flexural Bracing spacing (For fully braced condition, use a small number close to 0)
    E = E*ksi
    G = G*ksi
    poiss = 0.3
    L_y = F_bracing
    
    a = a*inch
    s_max = L_stud/6
    if a > s_max:
        return 'ERROR: Interconnection spacing is inadequate'

    s_section_path = os.path.dirname(__file__) + '\CFS_S_SECTION_DATA.csv'
    rawPath_s_section = r'' + s_section_path
    
    list_sections = []
    if section_type == 'stud':
        with open(rawPath_s_section, 'r', encoding="utf-8-sig") as csv_file:
            csv_reader = csv.DictReader(csv_file)
            keys = csv_reader.fieldnames
            
            for line in csv_reader:
                list_sections.append(line['ID'])
            
            if target_id not in list_sections:
                return 'ERROR: Section was not found in current reference file. Check your input or add the section manually.'
            else:
                location = list_sections.index(target_id)
                csv_file.seek(0)
                next(csv_reader)
                for count,line in enumerate(csv_reader):
                    if count == location:
                        Depth = float(line[keys[1]])
                        base = float(line[keys[2]])
                        t = float(line[keys[3]])
                        D = float(line[keys[4]])
                        R = float(line[keys[5]])
                        Area = float(line[keys[6]])
                        Ix = float(line[keys[8]])
                        Iy = float(line[keys[11]])
                        x_m = float(line[keys[14]])
                        m = float(line[keys[15]])
                        J = float(line[keys[16]])
                        Cw = float(line[keys[17]])
    else:
        return 'ERROR: Incorrect section used.'

    if target_id in list_sections:
        Depth = Depth*inch
        Ag = Area*inch**2
        Ixxg = Ix*inch**4
        Iyyg = Iy*inch**4
        Cw = Cw*inch**6
        base = base*inch
        t = t*inch
        R = R*inch
        if section_type == 'stud':
            D = D*inch # this is the small d value - lip dimension
        x_m = x_m*inch
        J = J*inch**4
        m = m*inch
        r = R + t/2
    # Boxed section properties:
    Ixxgt = 2*Ixxg
    Iyygt = 2*(Iyyg+Ag*((base - x_m)**2))
    S_ftot = (Ixxgt/(Depth/2))
    S_fy_tot = S_ftot
    # revision of limits for shape with holes:
    procedure_alt = False
    if P_cond == True:
        d_h = P_width
        L_h = P_length
        h = Depth - 2*(R + t)
        clear_dist = 24*inch - P_length
        
        if d_h/h <= 0.7 and h/t <= 200 and clear_dist >= 18*inch and P_corner_radi>= 2*t and d_h<=2.5*inch and L_h <= 4.5*inch and d_h > (9/16)*inch:
            if d_h/h >= 0.38:
                procedure_alt = True # procedure is altered
            else:
                procedure_alt = False # procedure remains the same, no alterations needed, hole ignored
        else:
            return 'ERROR: SECTION DOES NOT COMPLY WITH ALL LIMITATIONS, CANNOT CHECK' 
    # YIELDING AND GLOBAL LATERAL TORSIONAL BUCKLING:
    a_J = (base*2) - t
    b_J = (Depth) - t
    J_boxed = (2*(a_J*b_J)**2)/((a_J/t)+(b_J/t))
    L_u = ((0.36*Cb*pi)/(Fy*S_ftot))*(E*G*J_boxed*Iyygt)**(1/2)
    L_u = L_u.to('ft')
    
    if L_stud <= L_u:
        F_n = Fy
    else:
        F_cre = ((Cb*pi)/(K_y*L_y*S_fy_tot))*(E*G*J_boxed*Iyygt)**(1/2)
        if F_cre >= (2.78*Fy): F_n = Fy
        elif F_cre > (0.56*Fy) and F_cre < (2.78*Fy): F_n = (10/9)*Fy*(1-((10*Fy)/(36*F_cre)))
        elif F_cre <= (0.56*Fy): F_n = F_cre
    
    # LOCAL BUCKLING INTERACTING WITH YIELDING AND GLOBAL BUCKLING:
    if procedure_alt == False:
        S_e = calculateEffectiveSectionModulus_BuiltUp(Depth, F_n, procedure_alt, section_type, base, R, t, poiss, E, r, D, d_h=0, P_width=P_width, AnalysisType='F_n', builtUpSectionType='boxed') # At max. stress
        S_et = calculateEffectiveSectionModulus_BuiltUp(Depth, Fy, procedure_alt, section_type, base, R, t, poiss, E, r, D, d_h=0, P_width=P_width, AnalysisType='Fy', builtUpSectionType='boxed')
    else:
        S_e = calculateEffectiveSectionModulus_BuiltUp(Depth, F_n, procedure_alt, section_type, base, R, t, poiss, E, r, D, d_h, P_width, AnalysisType='F_n', builtUpSectionType='boxed') # At max. stress
        S_et = calculateEffectiveSectionModulus_BuiltUp(Depth, Fy, procedure_alt, section_type, base, R, t, poiss, E, r, D, d_h, P_width, AnalysisType='Fy', builtUpSectionType='boxed')
    
    M_alo = min(S_e*Fy, S_et*Fy)
    M_alo = M_alo/1.67
    M_alo_final = M_alo.to('kipft')
    return M_alo_final

def getFlexuralStrength_B2B_GB(target_id:str, section_type:str, L_stud:float, F_bracing:float, Fy:float, E:float, G:float, P_length:float, P_width:float, P_cond:bool, K_y:float, K_t:float, L_x:float, Cb:float):
    '''Calculation of globally braced CFS back to back member bending capacity.
    Calculations apply to back to back stud sections only. Other back to back sections are not available.  
    
    Args:
        target_id (str): Section shape name.
        section_type (str): 'stud' or 'track'.
        L_stud (float): Total element length for interconnection spacing revision (ft).
        F_bracing (float): Weak axis flexural bracing distance applicable for torsion as well (ft).
        Fy (float): Yield stress (ksi).
        E (float): Elastic modulus (ksi).
        G (float): Shear Modulus (ksi).
        P_length (float): Punchout length, typically 4.0 inches (in)
        P_width (float): Punchout width, typically 1.5 inches (in).
        P_cond (bool): Punchout condition.
        K_y (float): Effective length factor y-axis.
        K_t (float): Effective length factor for twisting.
        L_x (float): Strong axis unbraced length (ft).
        Cb (float): Bending coefficient.
     
    Returns:
        forallpeople.Physical: Allowable bending strength of globally braced back to back stud member.
    '''
    # Inputs: 
    target_id = target_id
    section_type = section_type
    L_stud = L_stud*ft
    Fy = Fy*ksi
    P_length = P_length*inch
    P_width = P_width*inch
    P_corner_radi = 0.25*inch # Punchout corner radii
    F_bracing = F_bracing*ft # Flexural Bracing spacing (For fully braced condition, use a small number close to 0)
    E = E*ksi
    G = G*ksi
    poiss = 0.3
    L_x = L_x*ft
    L_t = F_bracing
    L_y = F_bracing

    s_section_path = os.path.dirname(__file__) + '\CFS_S_SECTION_DATA.csv'
    rawPath_s_section = r'' + s_section_path
    
    list_sections = []
    if section_type == 'stud':
        with open(rawPath_s_section, 'r', encoding="utf-8-sig") as csv_file:
            csv_reader = csv.DictReader(csv_file)
            keys = csv_reader.fieldnames
            
            for line in csv_reader:
                list_sections.append(line['ID'])
            
            if target_id not in list_sections:
                return 'ERROR: Section was not found in current reference file. Check your input or add the section manually.'
            else:
                location = list_sections.index(target_id)
                csv_file.seek(0)
                next(csv_reader)
                for count,line in enumerate(csv_reader):
                    if count == location:
                        Depth = float(line[keys[1]])
                        base = float(line[keys[2]])
                        t = float(line[keys[3]])
                        D = float(line[keys[4]])
                        R = float(line[keys[5]])
                        Area = float(line[keys[6]])
                        Ix = float(line[keys[8]])
                        Iy = float(line[keys[11]])
                        x_m = float(line[keys[14]])
                        m = float(line[keys[15]])
                        J = float(line[keys[16]])
                        Cw = float(line[keys[17]])
                        xo = float(line[keys[20]])
    else:
        return 'ERROR: Incorrect section used'

    if target_id in list_sections:
        Depth = Depth*inch
        Ag = Area*inch**2
        Ixxg = Ix*inch**4
        Iyyg = Iy*inch**4
        x_o = xo*inch
        Cw = Cw*inch**6
        base = base*inch
        t = t*inch
        R = R*inch
        if section_type == 'stud':
            D = D*inch # this is the small d value - lip dimension
        x_m = x_m*inch
        J = J*inch**4
        m = m*inch
        r = R + t/2
    
    # BUILT UP SECTION PROPERTIES:
    Atg = 2*Ag
    Ixxgt = 2*Ixxg
    Iyygt = 2*(Iyyg+Ag*((x_m)**2))
    Cwt = 2*Cw
    xot = 0*inch
    rxt = sqrt(Ixxgt/Atg)*inch
    ryt = sqrt(Iyygt/Atg)*inch
    rot = sqrt(rxt**2 + ryt**2 + xot**2)*inch
    J_tot = J*2
    S_ftot = (Ixxgt/(Depth/2))
    
    # revision of limits for shape with holes -  single section only:
    procedure_alt = False
    if P_cond == True:
        d_h = P_width
        L_h = P_length
        h = Depth - 2*(R + t)
        clear_dist = 24*inch - P_length
        
        if d_h/h <= 0.7 and h/t <= 200 and clear_dist >= 18*inch and P_corner_radi>= 2*t and d_h<=2.5*inch and L_h <= 4.5*inch and d_h > (9/16)*inch:
            limitations_met = True # check if all limitations are met. If True, Fcre does not need to be modified
            if d_h/h >= 0.38:
                procedure_alt = True # procedure is altered
            else:
                procedure_alt = False # procedure remains the same, no alterations needed, hole ignored
        else:
            return 'ERROR: SECTION DOES NOT COMPLY WITH ALL LIMITATIONS, CANNOT CHECK'
        
        # Average Properties:
        hole_area = P_width*t
        hole_length = P_width

        if section_type == 'stud':
            web_L = Depth - 2*(R+t)
            corner_length = (2*pi*(R + t/2))/4
            flange_length = base - 2*(R + t)
            lip_length = D - (R+t)
            Lg = web_L + 4*(corner_length)+2*(flange_length)+2*(lip_length) # segment lenght without holes
            L_net = Lg - hole_length
            A_net = Ag - hole_area
            L = Lg + L_net
            
            # NET INERTIA CALCULATOR STUD SECTION CALCULATED MANUALLY (NOT IN AISI MANUAL)
            I_net_xx = Ixxg - ((t*(P_width**3))/12)
            I_net_yy = Iyyg - (((P_width*(t**3))/12)+P_width*t*((x_m - t/2)**2))
            I_x_avg = (Ixxg*Lg + I_net_xx*L_net)/(L)
            I_y_avg = (Iyyg*Lg + I_net_yy*L_net)/(L)
            
            block_depth = ((Depth/2)-(P_width/2)-t)
            lip_block_depth = D - t
            
            x_m_net = (2*(block_depth)*(t)*(t/2)+2*(t)*(base)*(base/2)+2*lip_block_depth*t*(base-(t/2)))/(2*(block_depth)*(t)+2*t*base+2*t*lip_block_depth)
            
            A_avg = (Ag*Lg + A_net*L_net)/(L)
            x_o_net = (m-(t/2))+x_m_net
            y_o_net = 0*inch
            y_o = 0*inch
            x_o_avg = (x_o*Lg + x_o_net*L_net)/L
            y_o_avg = (y_o*Lg + y_o_net*L_net)/L
            ro_avg = sqrt((x_o_avg**2)+(y_o_avg**2)+((I_x_avg+I_y_avg)/A_avg))*inch
    
    # YIELDING AND GLOBAL (LATERAL TORSIONAL) BUCKLINGL:
    
    # Initiation of Yielding Strength
    sigma_ey = ((pi**2)*E)/(((K_y*L_y)/(ryt))**2)
    sigma_t = (1/(Atg*rot**2))*(G*J_tot+(((pi**2)*E*Cwt)/((K_t*L_t)**2)))
    if P_cond == True and limitations_met == False: F_cre = ((Cb*ro_avg*Atg)/(S_ftot))*sqrt(sigma_ey*sigma_t)*ksi
    else: F_cre = ((Cb*rot*Atg)/(S_ftot))*sqrt(sigma_ey*sigma_t)*ksi
    
    if F_cre >= (2.78*Fy): F_n = Fy
    elif F_cre > (0.56*Fy) and F_cre < (2.78*Fy): F_n = (10/9)*Fy*(1-((10*Fy)/(36*F_cre)))
    elif F_cre <= (0.56*Fy): F_n = F_cre
    
    if P_cond == False:
        S_e = calculateEffectiveSectionModulus_BuiltUp(Depth, F_n, procedure_alt, section_type, base, R, t, poiss, E, r, D, d_h=0, P_width=P_width, AnalysisType='F_n', builtUpSectionType='b2b2')
        S_et = calculateEffectiveSectionModulus_BuiltUp(Depth, Fy, procedure_alt, section_type, base, R, t, poiss, E, r, D, d_h=0, P_width=P_width, AnalysisType='Yielding', builtUpSectionType='b2b2')
    else: 
        S_e = calculateEffectiveSectionModulus_BuiltUp(Depth, F_n, procedure_alt, section_type, base, R, t, poiss, E, r, D, d_h, P_width,'F_n', 'b2b2')
        S_et = calculateEffectiveSectionModulus_BuiltUp(Depth, Fy, procedure_alt, section_type, base, R, t, poiss, E, r, D, d_h, P_width,'Yielding', 'b2b2')
    
    print('Se: ', S_e)
    print('S_et: ', S_et)
    
    M_alo = min(S_e*Fy, S_et*Fy)
    M_alo = M_alo/1.67
    M_alo_final = M_alo.to('kipft')
    return M_alo_final

# SHEAR STRENGTH FUNCTIONS: ------------------------------------------------------------------------------------
def getShearStrength_Single(target_id:str, section_type:str, Fy:float, E:float, G:float, P_width:float, P_cond:bool):
    '''Calculation of single stud/track member shear capacity.
    Calculations apply to both stud and track sections, including punched and unpunched.   
    
    Args:
        target_id (str): Section shape name.
        section_type (str): 'stud' or 'track'.
        Fy (float): Yield stress (ksi).
        E (float): Elastic modulus (ksi).
        G (float): Shear Modulus (ksi).
        P_width (float): Punchout width, typically 1.5 inches (in).
        P_cond (bool): Punchout condition.
     
    Returns:
        forallpeople.Physical: Allowable shear strength of single stud/track member.
    '''
    finalShearStrength = getSingleShear(target_id, section_type, Fy, E, G, P_width, P_cond)
    return finalShearStrength

def getShearStrength_BuiltUp(target_id:str, section_type:str, Fy:float, E:float, G:float, P_width:float, P_cond:bool):
    '''Calculation of built up CFS member shear capacity (boxed or back to back)
    Calculations apply to both stud and track sections, including punched and unpunched.
    
    Args:
        target_id (str): Section shape name.
        section_type (str): 'stud' or 'track'.
        Fy (float): Yield stress (ksi).
        E (float): Elastic modulus (ksi).
        G (float): Shear Modulus (ksi).
        P_width (float): Punchout width, typically 1.5 inches (in).
        P_cond (bool): Punchout condition.
    
    Returns:
        forallpeople.Physical: Allowable shear strength of built up CFS member (boxed or back to back).
    
    '''
    # Applicable for both boxed and back to back built up sections
    finalShearStrength = getSingleShear(target_id, section_type, Fy, E, G, P_width, P_cond)*2
    return finalShearStrength

# HSS SECTION DESIGN LIMIT STATES: -----------------------------------------------------------------------------
def getHSSAxialStrength(hss_section:str, Fy:float, L_x:float, L_y:float):
    '''Calculation of HSS section axial capacity.
    
    Args:
        hss_section (str): HSS section shape name.
        Fy (float): Yield stress (ksi).
        L_x (float): Strong axis unbraced length (ft).
        L_y (float): Weak axis unbraced length (ft).
    
    Returns:
        forallpeople.Physical: Allowable axial strength of HSS section.
    
    '''
    # Get relevant section data: 
    E = 29000*ksi
    K_x = 1
    K_y = 1
    L_x = L_x*ft
    L_y = L_y*ft
    Fy = Fy*ksi
    hss_section_path = os.path.dirname(__file__) + '\HSS_SECTION_DATA_JAMBS.csv'
    rawPath_hss_section = r'' + hss_section_path
    list_sections = []
    try:
        with open(rawPath_hss_section, 'r', encoding="utf-8-sig") as csv_file:
            csv_reader = csv.DictReader(csv_file)
            keys = csv_reader.fieldnames
            for line in csv_reader:
                list_sections.append(line['SECTION'])
            if hss_section not in list_sections:
                return 'ERROR: Section not found in csv file'
            else:
                sectionLocation = list_sections.index(hss_section)
                csv_file.seek(0)
                next(csv_reader)
                for count, line in enumerate(csv_reader):
                    if count == sectionLocation:
                        A_g = float(line[keys[1]])
                        I_x = float(line[keys[2]])
                        r_x = float(line[keys[3]])
                        I_y = float(line[keys[4]])
                        r_y = float(line[keys[5]])
                        Height = float(line[keys[6]])
                        Width = float(line[keys[7]])
                        b_t_rel = float(line[keys[8]])
                        h_t_rel = float(line[keys[9]])
                        t = float(line[keys[10]])
    except FileNotFoundError as error:
        return 'ERROR: File not found - ' + error
    if hss_section in list_sections:
        A_g = A_g*inch**2
        I_x = I_x*inch**4
        r_x = r_x*inch
        I_y = I_y*inch**4
        r_y = r_y*inch
        Height = Height*inch
        Width = Width*inch
        t = t*inch
    
    # Start axial calculations
    slenderness_x = (K_x*L_x)/(r_x)
    slenderness_y = (K_y*L_y)/(r_y)
    if slenderness_x >200 or slenderness_y >200:
        return 'ERROR: Section is too slender'
    # Obtain total axial capacity (compression)      
    lambda_r = 1.4*sqrt(E/Fy)
    if b_t_rel > lambda_r:
        flange_class = 'slender'
    else: 
        flange_class = 'non-slender'
    if h_t_rel > lambda_r:
        web_class = 'slender'
    else: 
        web_class = 'non-slender'       
    if flange_class=='slender' or web_class=='slender':
        shape_class = 'slender'
    else: 
        shape_class = 'non-slender'    
    critical_slenderness = max(slenderness_x, slenderness_y)
    F_e = ((pi**2)*E)/((critical_slenderness)**2)
    calc_limit = 4.71*sqrt(E/Fy)

    if critical_slenderness <= calc_limit:
        F_cr = (0.658**(Fy/F_e))*Fy
    else: 
        F_cr = 0.877*F_e
    lambda_comp = lambda_r*sqrt(Fy/F_cr)
    c1 = 0.2
    c2 = 1.38
    Fel_b = ((c2*(lambda_r/b_t_rel))**2)*Fy
    Fel_h = Fel_b = ((c2*(lambda_r/h_t_rel))**2)*Fy
    if b_t_rel <= lambda_comp:
        b_e = Width
    else:
        b_e = Width*(1-c1*sqrt(Fel_b/F_cr))*sqrt(Fel_b/F_cr)
        
    if h_t_rel <= lambda_comp:
        h_e = Height
    else:
        h_e = Height*(1-c1*sqrt(Fel_h/F_cr))*sqrt(Fel_h/F_cr)

    A_e = A_g - (2*(Width-b_e)*t + 2*(Height-h_e)*t)
    
    Pn_LB = A_e*F_cr
    Pn_FLB = A_g*F_cr
    P_n = min(Pn_LB, Pn_FLB)
    P_allow = P_n/1.67

    return P_allow

def getHSSFlexuralStrength(hss_section:str, Fy:float, Length:float, orientation:str):
    '''Calculation of HSS section bending capacity.
    
    Args:
        hss_section (str): HSS section shape name.
        Fy (float): Yield stress (ksi).
        Length (float): Total member length (ft).
        orientation (str): Member orientation ("strong" or "weak").
    
    Returns:
        forallpeople.Physical: Allowable bending strength of HSS section.
    
    '''
    # DIFFERENT LIMIT STATES TO CHECK ACCORDING TO THE AISC 14TH EDITION MANUAL:
    #       YIELDING
    #       FLANGE LOCAL BUCKLING 
    #       WEB LOCAL BUCKLING
    # Get relevant section data:
    # Note on Lateral Torsional Buckling:
    # Very long rectangular HSS bent about the major axis are subject tolateral-torsional buckling; however, the Specification provides no strength equation for this limit state since beam deflection will control for all reasonable cases
    
    orientation = orientation.lower()
    if orientation != 'strong' and orientation != 'weak':
        return 'ERROR: Incorrect section orientation'
        
    E = 29000*ksi
    Length= Length*ft
    Fy = Fy*ksi
    hss_section_path = os.path.dirname(__file__) + '\HSS_SECTION_DATA_JAMBS.csv'
    rawPath_hss_section = r'' + hss_section_path
    list_sections = []
    try:
        with open(rawPath_hss_section, 'r', encoding="utf-8-sig") as csv_file:
            csv_reader = csv.DictReader(csv_file)
            keys = csv_reader.fieldnames
            for line in csv_reader:
                list_sections.append(line['SECTION'])
            if hss_section not in list_sections:
                return 'ERROR: Section not found in csv file'
            else:
                sectionLocation = list_sections.index(hss_section)
                csv_file.seek(0)
                next(csv_reader)
                for count, line in enumerate(csv_reader):
                    if count == sectionLocation:
                        A_g = float(line[keys[1]])
                        I_x = float(line[keys[2]])
                        r_x = float(line[keys[3]])
                        I_y = float(line[keys[4]])
                        r_y = float(line[keys[5]])
                        Height = float(line[keys[6]])
                        Width = float(line[keys[7]])
                        b_t_rel = float(line[keys[8]])
                        h_t_rel = float(line[keys[9]])
                        t = float(line[keys[10]])
                        Z_x = float(line[keys[11]])
                        S_x = float(line[keys[12]])
                        J = float(line[keys[13]])
                        Z_y = float(line[keys[15]])
                        S_y = float(line[keys[16]])
    except FileNotFoundError as error:
        return 'ERROR: File not found - ' + error
    if hss_section in list_sections:
        A_g = A_g*inch**2
        I_x = I_x*inch**4
        r_x = r_x*inch
        I_y = I_y*inch**4
        r_y = r_y*inch
        Height = Height*inch
        Width = Width*inch
        t = t*inch
        Z_x = Z_x*inch**3
        S_x = S_x*inch**3
        J = J*inch**4
        Z_y = Z_y*inch**3
        S_y = S_y*inch**3
    # Check section type based on limits in AISC Manual
    lambda_p_flanges = 1.12*sqrt(E/Fy)
    lambda_r_flanges = 1.40*sqrt(E/Fy)
    lambda_p_web = 2.42*sqrt(E/Fy)
    lambda_r_web = 5.70*sqrt(E/Fy)
    
    # check flanges: 
    
    if orientation == 'strong':
        if b_t_rel < lambda_p_flanges:
                flange_class = 'compact'
        elif b_t_rel <= lambda_r_flanges and b_t_rel >= lambda_p_flanges:
            flange_class = 'non-compact'
        else:
            flange_class = 'slender'
            
        # check webs:
        if h_t_rel < lambda_p_web:
            web_class = 'compact'
        elif h_t_rel <= lambda_r_web and h_t_rel >= lambda_p_web:
            web_class = 'non-compact'
        else:
            web_class = 'slender'
        # Start flexural calculations:
        M_p = Z_x*Fy # Yielding

        # flange local buckling: 
        b = Width-3*t
        h = Height-3*t
        M_n_flange = 99999999*kip*ft
        M_n_web = 99999999*kip*ft

        if flange_class == 'non-compact':
            # compared_value = M_p*(M_p-Fy*S_x)*(3.57*(b_t_rel)*sqrt(Fy/E)-4.0)
            # print(compared_value)
            M_n_flange = min(M_p-(M_p-Fy*S_x)*(3.57*(b_t_rel)*sqrt(Fy/E)-4.0), M_p)
        elif flange_class == 'slender':
            b_e = min(1.92*t*sqrt(E/Fy)*(1-(0.38/b_t_rel)*sqrt(E/Fy)), b)
            c_y = Height  - ((b_e*t*(Height-0.5*t) + Height*t*Height*0.5*2 + t*(Width-t*2)*(0.5*t))/(b_e*t + 2*(Height*t) + t*(Width-t*2)))
            I_x_effect = I_x - t*(b-b_e)*((0.5*h-0.5*t)**2)-(((t**3)*(b-b_e))/12)
            S_e = I_x_effect/c_y    
            M_n_flange = Fy*S_e
            
        # web local buckling: 
        if web_class == 'non-compact':
            M_n_web = min(M_p - (M_p-Fy*S_x)*(0.305*(h_t_rel)*sqrt(Fy/E)-0.738), M_p)
        elif web_class == 'slender':
            return 'ERROR: HSS section has a slender web - check calculations'
            
        M_n_final = min(M_p, M_n_flange, M_n_web)
        M_allow = M_n_final/1.67
        M_allow = M_allow.to('kipft')
        return M_allow
    else:
        # weak axis orientation
        if h_t_rel < lambda_p_flanges:
                flange_class = 'compact'
        elif h_t_rel <= lambda_r_flanges and b_t_rel >= lambda_p_flanges:
            flange_class = 'non-compact'
        else:
            flange_class = 'slender'
            
        # check webs:
        if b_t_rel < lambda_p_web:
            web_class = 'compact'
        elif b_t_rel <= lambda_r_web and b_t_rel >= lambda_p_web:
            web_class = 'non-compact'
        else:
            web_class = 'slender'
        # Start flexural calculations:
        M_p = Z_y*Fy # Yielding

        # flange local buckling: 
        b = Width-3*t
        h = Height-3*t
        M_n_flange = 99999999*kip*ft
        M_n_web = 99999999*kip*ft

        if flange_class == 'non-compact':
            M_n_flange = min(M_p-(M_p-Fy*S_y)*(3.57*(h_t_rel)*sqrt(Fy/E)-4.0), M_p)
        elif flange_class == 'slender':
            h_e = min(1.92*t*sqrt(E/Fy)*(1-(0.38/h_t_rel)*sqrt(E/Fy)), h)
            c_x = Width  - ((h_e*t*(Width-0.5*t) + Width*t*Width*0.5*2 + t*(Height-t*2)*(0.5*t))/(h_e*t + 2*(Width*t) + t*(Height-t*2)))
            I_y_effect = I_y-t*(h-h_e)*((0.5*Width-0.5*t)**2)-(((t**3)*(h-h_e))/(12))
            S_e = I_y_effect/c_x
            M_n_flange = Fy*S_e
            
        # web local buckling: 
        if web_class == 'non-compact':
            M_n_web = min(M_p - (M_p-Fy*S_y)*(0.305*(b_t_rel)*sqrt(Fy/E)-0.738), M_p)
        elif web_class == 'slender':
            return 'ERROR: HSS section has a slender web - check calculations'
            
        M_n_final = min(M_p, M_n_flange, M_n_web)
        M_allow = M_n_final/1.67
        M_allow = M_allow.to('kipft')
        return M_allow

def getHSSShearStrength(hss_section:str, Fy:float, orientation:str):
    '''Calculation of HSS section shear capacity.
    
    Args:
        hss_section (str): HSS section shape name.
        Fy (float): Yield stress (ksi).
        orientation (str): Member orientation ("strong" or "weak").
    
    Returns:
        forallpeople.Physical: Allowable shear strength of HSS section.
    
    '''
    orientation = orientation.lower()
    if orientation != 'strong' and orientation != 'weak':
        return 'ERROR: Incorrect section orientation'
    E = 29000*ksi
    Fy = Fy*ksi
    hss_section_path = os.path.dirname(__file__) + '\HSS_SECTION_DATA_JAMBS.csv'
    rawPath_hss_section = r'' + hss_section_path
    list_sections = []
    try:
        with open(rawPath_hss_section, 'r', encoding="utf-8-sig") as csv_file:
            csv_reader = csv.DictReader(csv_file)
            keys = csv_reader.fieldnames
            for line in csv_reader:
                list_sections.append(line['SECTION'])
            if hss_section not in list_sections:
                return 'ERROR: Section not found in csv file'
            else:
                sectionLocation = list_sections.index(hss_section)
                csv_file.seek(0)
                next(csv_reader)
                for count, line in enumerate(csv_reader):
                    if count == sectionLocation:
                        A_g = float(line[keys[1]])
                        I_x = float(line[keys[2]])
                        r_x = float(line[keys[3]])
                        I_y = float(line[keys[4]])
                        r_y = float(line[keys[5]])
                        Height = float(line[keys[6]])
                        Width = float(line[keys[7]])
                        b_t_rel = float(line[keys[8]])
                        h_t_rel = float(line[keys[9]])
                        t = float(line[keys[10]])
                        Z_x = float(line[keys[11]])
                        S_x = float(line[keys[12]])
                        J = float(line[keys[13]])
                        Z_y = float(line[keys[15]])
                        S_y = float(line[keys[16]])
    except FileNotFoundError as error:
        return 'ERROR: File not found - ' + error
    if hss_section in list_sections:
        A_g = A_g*inch**2
        I_x = I_x*inch**4
        r_x = r_x*inch
        I_y = I_y*inch**4
        r_y = r_y*inch
        Height = Height*inch
        Width = Width*inch
        t = t*inch
        Z_x = Z_x*inch**3
        S_x = S_x*inch**3
        J = J*inch**4
        Z_y = Z_y*inch**3
        S_y = S_y*inch**3
    
    h = Height - 3*t
    b = Width - 3*t
    k_v = 5
    if orientation == 'strong':
        A_w = 2*h*t    
        if h_t_rel <= 1.10*sqrt((k_v*E)/Fy):
            C_v = 1.0
        elif 1.10*sqrt((k_v*E)/Fy) < h_t_rel <= 1.37*sqrt((k_v*E)/Fy):
            C_v = (1.10*sqrt((k_v*E)/Fy))/(h_t_rel)
        elif h_t_rel > 1.37*sqrt((k_v*E)/Fy):
            C_v = (1.51*k_v*E)/(((h_t_rel)**2)*Fy)
    else:
        A_w = 2*b*t
        if b_t_rel <= 1.10*sqrt((k_v*E)/Fy):
            C_v = 1.0
        elif 1.10*sqrt((k_v*E)/Fy) < b_t_rel <= 1.37*sqrt((k_v*E)/Fy):
            C_v = (1.10*sqrt((k_v*E)/Fy))/(b_t_rel)
        elif b_t_rel > 1.37*sqrt((k_v*E)/Fy):
            C_v = (1.51*k_v*E)/(((b_t_rel)**2)*Fy)

    V_n = 0.6*Fy*A_w*C_v
    V_allow = V_n/1.67
    V_allow = V_allow.to('kip')
    
    return V_allow

