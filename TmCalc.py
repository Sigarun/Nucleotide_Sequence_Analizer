from NucleicAcid import NucleicAcid
import math
import numbers

class TmCalculator(NucleicAcid):
    """
    Nucleic acid melting temperature calculator.
    For now only the nearest neighbour method is available, basic and salt
    adjusted are to be added.
    The class uses sequence, primer concentration(nM) and Na ions concentration(mM)
    to calculate Tm value(s) along with the basic thermodinamic propherties
    of a sequence.
    """

    #Thermodynamic values charateristic for nucleotid pairs
    #(deltaH, deltaS, deltaG(37degreeC))
    TD_VALS = {
        "AA":(8.0, 21.9, 1.2),
        "TT":(8.0, 21.9, 1.2),
        "AT":(5.6, 15.2, 0.9),
        "TA":(6.6, 18.4, 0.9),
        "CA":(8.2, 21.0, 1.7),
        "TG":(8.2, 21.0, 1.7),
        "GT":(9.4, 25.5, 1.5),
        "AC":(9.4, 25.5, 1.5),
        "CT":(6.6, 16.4, 1.5),
        "AG":(6.6, 16.4, 1.5),
        "GA":(8.8, 23.5, 1.5),
        "TC":(8.8, 23.5, 1.5),
        "CG":(11.8, 29.0, 2.8),
        "GC":(10.5, 26.4, 2.3),
        "GG":(10.9, 28.4, 2.1),
        "CC":(10.9, 28.4, 2.1)
    }

    def __init__(self, seq, seq_type, primer_conc = 50, na_conc = 50):
        super().__init__(seq, seq_type)
        if (not isinstance(primer_conc, numbers.Number) and
            not isinstance(na_conc, numbers.Number)):
            raise TypeError("Concentrations have to be of numeric type")
        self.primer_conc = primer_conc
        self.na_conc = na_conc

    def nearest_neighbour(self):
        """
        Nearest neighbour method.
        """
        salt_adj = 16.6 * math.log10(self.na_conc*0.001)
        self.rlnk = 1.987 * math.log(1/(self.primer_conc*0.000000001))

        vals = self.__match_td_vals()
        zipped = list(zip(*vals))
        self.delta_H, self.delta_S, self.delta_G = (sum(zipped[0]),
                                                    sum(zipped[1]),
                                                    sum(zipped[2]))
        #print("deltaH = %f, deltaS = %f, deltaG = %f" % (delta_H, delta_S, delta_G))
        #print(rlnk)
        #print(salt_adj)

        t = 1000 * ((self.delta_H - 3.4) / (self.delta_S + self.rlnk))
        self.tm = (t + salt_adj) - 273.15

    def __get_pairs(self):
        """
        Creates a list of pairs. For seq 'ACTG' - AC, CT, TG
        """
        pairs = [self.seq[i:i+2] for i in range(len(self.seq)-1)]
        return pairs

    def __match_td_vals(self):
        """
        Matches input pairs with their thermodynamic values
        and returns them as a list of tuples - (dH, dS, dG)
        """
        pairs = self.__get_pairs()
        vals = []
        for pair in pairs:
            try:
                vals.append(self.TD_VALS[pair])
            except KeyError:
                pass
        return vals
