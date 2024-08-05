from zigzag.hardware.architecture.get_cacti_cost import get_cacti_cost

class Cell:
    def __init__(self, tech_param):
        self.tech_param = tech_param
        #self.cells_size = 
        #self.weight_precision = weight_precision

    
    def get_1b_multiplier_dly(self) -> float:
        """delay of 1b multiplier
        1b mult includes 1 NOR gate, which is assumed as the same cost of ND2 gate"""
        return self.tech_param["nd2_dly"]
    
    def get_1b_multiplier_energy(self) -> float:
        """energy of 1b multiplier
        1b mult includes 1 NOR gate, which is assumed as the same cost of ND2 gate
        why 0.5: considering weight stays constant during multiplication"""
        return 0.5 * self.tech_param["nd2_cap"] * (self.tech_param["vdd"] ** 2)  # unit: pJ
    
    def get_1b_multiplier_area(self) -> float:
        """area of 1b multiplier
        1b mult includes 1 NOR gate, which is assumed as the same cost of ND2 gate"""
        return self.tech_param["nd2_area"]
    
    @staticmethod
    def get_single_cell_array_cost_from_cacti(tech_node: float, wordline_dim_size: float, bitline_dim_size: float, cells_size: float, weight_precision: int) -> tuple[float, float, float, float]:
        cell_array_size = wordline_dim_size * bitline_dim_size * cells_size / 8
        array_bw = wordline_dim_size * weight_precision

        cacti_path = "zigzag/cacti/cacti_master"
        access_time, area, r_cost, w_cost = get_cacti_cost(
            cacti_path=cacti_path,
            tech_node=tech_node,
            mem_type="sram",
            mem_size_in_byte=cell_array_size,
            bw=array_bw,
        )
        return access_time, area, r_cost, w_cost
