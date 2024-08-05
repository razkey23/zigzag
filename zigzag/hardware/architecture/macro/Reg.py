class Reg:
    def __init__(self, tech_param):
        self.tech_param = tech_param
    
    def get_1b_reg_area(self) -> float:
        """area of 1b DFF"""
        return self.tech_param["dff_area"]
    
    # May need to decouple as well
    def get_1b_reg_energy(self) -> float:
        """energy of 1b DFF"""
        return self.tech_param["dff_cap"] * (self.tech_param["vdd"] ** 2)  # unit: pJ