class Adc:
    def __init__(self, tech_param,adc_resolution):
        self.tech_param = tech_param
        self.adc_resolution = adc_resolution
    
    def get_adc_cost(self, bitline_dim_size) -> tuple[float, float, float]:
        if self.adc_resolution == 0:
            adc_area = 0
            adc_delay = 0
            adc_energy = 0
        else:
            k1 = -0.0369
            k2 = 1.206
            adc_area = 10 ** (k1 * self.adc_resolution + k2) * 2**self.adc_resolution * (10**-6)
            k3 = 0.00653
            k4 = 0.640
            adc_delay = self.adc_resolution * (k3 * bitline_dim_size + k4)
            k5 = 100
            k6 = 0.001
            adc_energy = (k5 * self.adc_resolution + k6 * 4**self.adc_resolution) * self.tech_param["vdd"] ** 2 / 1000
        return adc_area, adc_delay, adc_energy