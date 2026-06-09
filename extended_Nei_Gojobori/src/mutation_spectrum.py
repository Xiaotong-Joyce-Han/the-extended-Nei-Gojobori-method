class MutationSpectrum:

    def __init__(self, rates: list[float]):

        if len(rates) != 6:
            raise ValueError(f"Expected number of parameters in mutation spectrum is 6, but got {len(rates)}")
        
        self.ATtoTA = rates[0]
        self.ATtoCG = rates[1]
        self.ATtoGC = rates[2]
        self.CGtoTA = rates[3]
        self.CGtoAT = rates[4]
        self.CGtoGC = rates[5]

        self.check_validity()

    
    def check_validity(self):
        AT_valid = self.ATtoTA >= 0 and self.ATtoCG >= 0 and self.ATtoGC >= 0 and (self.ATtoTA + self.ATtoCG + self.ATtoGC > 0)
        CG_valid = self.CGtoTA >= 0 and self.CGtoAT >= 0 and self.CGtoGC >= 0 and (self.CGtoTA + self.CGtoAT + self.CGtoGC > 0)
        valid = AT_valid and CG_valid

        if not valid:
            raise ValueError("Mutation spectrum specified is invalid.")
        
    
    def __str__(self):
        return (
        "Mutation spectrum:\n"
        f"  AT -> TA relative mutation rate: {self.ATtoTA}\n"
        f"  AT -> CG relative mutation rate: {self.ATtoCG}\n"
        f"  AT -> GC relative mutation rate: {self.ATtoGC}\n"
        f"  CG -> TA relative mutation rate: {self.CGtoTA}\n"
        f"  CG -> AT relative mutation rate: {self.CGtoAT}\n"
        f"  CG -> GC relative mutation rate: {self.CGtoGC}"
    ) 


    def get_mutation_rate(self, start_base, target_base):

        if start_base == target_base:
            return 0.0

        if start_base == "A":
            if target_base == "T":
                return self.ATtoTA
            elif target_base == "C":
                return self.ATtoCG
            elif target_base == "G":
                return self.ATtoGC

        elif start_base == "T":
            if target_base == "A":
                return self.ATtoTA
            elif target_base == "G":
                return self.ATtoCG
            elif target_base == "C":
                return self.ATtoGC

        elif start_base == "C":
            if target_base == "A":
                return self.CGtoAT
            elif target_base == "T":
                return self.CGtoTA
            elif target_base == "G":
                return self.CGtoGC

        elif start_base == "G":
            if target_base == "T":
                return self.CGtoAT
            elif target_base == "A":
                return self.CGtoTA
            elif target_base == "C":
                return self.CGtoGC
            
        raise ValueError(
            f"Invalid mutation type: {start_base} -> {target_base}. "
            "Bases must be one of A, T, C, G."
        )