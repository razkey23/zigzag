import math
class Adder:
    def __init__(self, tech_param):
        self.tech_param = tech_param

    def get_1b_adder_energy(self) -> float:
        """energy of 1b full adder
        Assume a 1b adder has 3 ND2 gate and 2 XOR2 gate"""
        adder_cap = 3 * self.tech_param["nd2_cap"] + 2 * self.tech_param["xor2_cap"]
        return adder_cap * (self.tech_param["vdd"] ** 2)  # unit: pJ
    
    def get_1b_adder_energy_half_activated(self) -> float:
        """energy of 1b full adder when 1 input is 0"""
        adder_cap = 2 * self.tech_param["xor2_cap"]
        return adder_cap * (self.tech_param["vdd"] ** 2)  # unit: pJ

    def get_1b_adder_area(self) -> float:
        adder_area = 3 * self.tech_param["nd2_area"] + 2 * self.tech_param["xor2_area"]
        return adder_area

    def get_1b_adder_dly_in2sum(self) -> float:
        adder_dly = 2 * self.tech_param["xor2_dly"]
        return adder_dly

    def get_1b_adder_dly_in2cout(self) -> float:
        adder_dly = self.tech_param["xor2_dly"] + 2 * self.tech_param["nd2_dly"]
        return adder_dly
    
    def get_1b_adder_dly_cin2cout(self) -> float:
        adder_dly = 2 * self.tech_param["nd2_dly"]
        return adder_dly