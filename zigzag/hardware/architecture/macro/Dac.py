import math

class Dac:
    def __init__(self, tech_param):
        self.tech_param = tech_param

    def get_dac_cost(self, bit_serial_precision) -> tuple[float, float, float]:
        dac_area = 0
        dac_delay = 0
        if bit_serial_precision == 1:
            dac_energy = 0
        else:
            k0 = 50e-3
            dac_energy = k0 * bit_serial_precision * self.tech_param["vdd"] ** 2
        return dac_area, dac_delay, dac_energy