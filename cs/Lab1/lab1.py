VREF = 10

# ADC
codOut = []
capacitors = [a for a in range(1, 7)]

def adcOutput(vin:float, vref:float, resolution:int):
    vmin = -vref/2
    vmax = -vmin
    mostSigBits = []
    leastSigBits = []
    threshold = vmax if vin > 0 else vmin
    mostSigBits.append(1 if threshold == vmax else 0)
    for bitn in range(resolution):
        if bitn == len(capacitors) - 1:
            bitn = 1
            pass
        if len(mostSigBits) < len(capacitors):
            mostSigBits.append(1 if vin >= threshold else 0)
            threshold = threshold * 2**-capacitors[bitn]
        else:
            leastSigBits.append(1 if vin >= threshold else 0)
            threshold = threshold * 2**-capacitors[bitn]
    return mostSigBits + leastSigBits

codOut = adcOutput(6.5, VREF, 11)