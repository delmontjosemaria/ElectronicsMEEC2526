VREFP = 10.0
VREFN = 0.0

RESOLUTION = 12 

#Existe o C6, que está antes do CB mas não faz parte do conjunto híbrido
M = 7
L = 6

#Condensador unitário
C_u = 1 

#Condensador de ponte
C_b = 2*C_u

weightsMSB = [C_u*2**i for i in range(M, 2*M)]
weightsMSB = weightsMSB[::-1]
weightC6 = 2**(M-1)*C_u
weightsLSB = [C_b*(2**L - 1)/2**(i-L) for i in range(1, L)]
weightsLSB = weightsLSB[::-1]
weights = weightsMSB + [weightC6] + weightsLSB

VLSB = (VREFP-VREFN)/2**RESOLUTION

def comparator(v_xp:float, v_xn:float):
    return 1 if v_xp > v_xn else 0

def sarAdcActivity(vip:float, vin:float):
    preword = []
    bitCount = 0
    
    #fase de sampling, todos os condensadores estão carregados
    vxp = vip
    vxn = vin
    out = 2
    
    while (bitCount < RESOLUTION):
        #primeira comparação, logo à entrada do sinal
        out = comparator(vxp, vxn)
        
        #comparações híbridas, deve haver redundância dos bits a serem verificados
        if bitCount < M:
            if out == 1:
                vxp-=weights[bitCount]*VLSB
            elif out == 0:
                vxn-=weights[bitCount]*VLSB
            else:
                break
        else:
            #faz a comparação dos condensadores para o monotónico, mas vistas a partir do CB (sub-radix) exceto o C6. 
            if out == 1:
                vxp-=weightsLSB[bitCount]*VLSB
            elif out == 0:
                vxn-=weightsLSB[bitCount]*VLSB
            else:
                break
            
        preword.append(out)
        bitCount += 1
    
    
    
    return preword
    
    
def decActivity(preword:list):
    return preword
    
