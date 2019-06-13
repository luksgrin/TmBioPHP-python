import sys
import math

def tm_Base_Stacking(
    primer_seq,  #Primer sequence only containing ATGC
    conc_primer, #Primer concentration in nM
    conc_salt,   #Salt concentration in mM
    conc_mg):	 #Mg2+ concentration in mM
    
    H = 0
    S = 0

    # from table at http://www.ncbi.nlm.nih.gov/pmc/articles/PMC19045/table/T2/ (SantaLucia, 1998)
    # enthalpy values
    H_values = {
        'AA' : -7.9,
        'AC' : -8.4,
        'AG' : -7.8,
        'AT' : -7.2,
        'CA' : -8.5,
        'CC' : -8.0,
        'CG' :-10.6,
        'CT' : -7.8,
        'GA' : -8.2,
        'GC' : -9.8,
        'GG' : -8.0,
        'GT' : -8.4,
        'TA' : -7.2,
        'TC' : -8.2,
        'TG' : -8.5,
        'TT' : -7.9}

    # entropy values
    S_values = {
        'AA' :-22.2,
        'AC' :-22.4,
        'AG' :-21.0,
        'AT' :-20.4,
        'CA' :-22.7,
        'CC' :-19.9,
        'CG' :-27.2,
        'CT' :-21.0,
        'GA' :-22.2,
        'GC' :-24.4,
        'GG' :-19.9,
        'GT' :-22.4,
        'TA' :-21.3,
        'TC' :-22.2,
        'TG' :-22.7,
        'TT' :-22.2}
        
    # effect on entropy by salt correction; von Ahsen et al 1999
    # Increase of stability due to presence of Mg;
    salt_effect = (conc_salt/1e3)+((conc_mg/1e3)*140)
    
    # effect on entropy
    S += (0.368*(len(primer_seq)-1)*math.log(salt_effect))
                 
    # terminal corrections. Santalucia 1998
    firstnuc = primer_seq[0]
    if (firstnuc == 'G') or (firstnuc == 'C'):
        H += 0.1
        S += 2.8
    if (firstnuc == 'A') or (firstnuc == 'T'):
        H += 2.3
        S += 4.1
    
    lastnuc = primer_seq[-1]
    if (firstnuc == 'G') or (firstnuc == 'C'):
        H += 0.1
        S += 2.8
    if (firstnuc == 'A') or (firstnuc == 'T'):
        H += 2.3
        S += 4.1
    
    # compute new H and S based on sequence. Santalucia 1998
    for n in range(len(primer_seq)-1):
        sub_seq = primer_seq[n: n + 2]
        H += H_values[sub_seq]
        S += S_values[sub_seq]
    
    Tm = ((1e3*H)/(S + (1.987*math.log(conc_primer/2e9)))) - 273.15
    
    Tm = round(Tm, 2) 
    H = round(H, 2)
    S = round(S, 2)
    
    return Tm, H, S
    
print('No degenerate primers are allowed\n'
    'Calculations are performed assuming primer length spans between 6 and 50 nucleotides')

Primer_seq = input('Primer: ')
Primer_conc = float(input('Primer concentration (in nM): '))
Salt_conc = float(input('Salt concentration (in mM): '))
Mg_conc = float(input('Magensium2+ concentration (in mM): '))

(Tm,
H,
S) = tm_Base_Stacking(
    Primer_seq,
    Primer_conc,
    Salt_conc,
    Mg_conc)
    
print('Melting temperature: ' + str(Tm) + 'ºC')
print('Enthalpy: ' + str(H))
print('Entropy: ' + str(S))