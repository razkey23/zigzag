import math
import logging as _logging

from zigzag.workload.layer_node import LayerNode
from zigzag.mapping.Mapping import Mapping
from zigzag.datatypes import OADimension
from zigzag.hardware.architecture._ImcUnit import _ImcUnit
from zigzag.utils import json_repr_handler

class _ImcArray(_ImcUnit):
    """definition of an Analog/Digital In-SRAM-Computing (A/DIMC) core
    constraint:
        -- activation precision must be in the power of 2.
        -- bit_serial_precision must be in the power of 2.
    """

    def __init__(
        self,
        is_analog_imc: bool,
        bit_serial_precision: int,
        input_precision: list[int],
        adc_resolution: int,
        cells_size: int,
        cells_area: float | None,
        dimension_sizes: dict[OADimension, int],
        auto_cost_extraction: bool = False,
    ):
        super().__init__(
            is_analog_imc=is_analog_imc,
            bit_serial_precision=bit_serial_precision,
            input_precision=input_precision,
            adc_resolution=adc_resolution,
            cells_size=cells_size,
            cells_area=cells_area,
            dimension_sizes=dimension_sizes,
            auto_cost_extraction=auto_cost_extraction,
        )

        self.get_area() # No need for mapping
        self.get_tclk() # No need for mapping
        self.tops_peak, self.topsw_peak, self.topsmm2_peak = self.get_macro_level_peak_performance()
    
    def get_area(self):

        # Getting Area of Cells
        if self.auto_cost_extraction:
            cost_per_bank = self.cell.get_single_cell_array_cost_from_cacti(
                tech_node=self.tech_param["tech_node"],
                wordline_dim_size=self.wordline_dim_size,
                bitline_dim_size=self.bitline_dim_size,
                cells_size=self.cells_size,
                weight_precision=self.weight_precision,
            )
            area_cells = cost_per_bank[1] * self.nb_of_banks
            self.cells_w_cost = cost_per_bank[3] / self.wordline_dim_size
        else:
            assert (
                self.cells_area is not None
            ), "No cells area given by user yet `auto_cost_extraction` is set to `false`"
            area_cells = self.total_unit_count * self.cells_area
        
        # Getting Area of ADCs, DACs if existent
        if self.is_aimc:
            area_dacs = self.dac.get_dac_cost(self.bit_serial_precision)[0] * self.bitline_dim_size * self.nb_of_banks
            area_adcs = self.adc.get_adc_cost(self.bitline_dim_size)[0] * self.weight_precision * self.wordline_dim_size * self.nb_of_banks
        else:
            area_dacs = 0
            area_adcs = 0
        
        # Multiplier Array (Extra gates for multiplication - based on 6T-Sram Cell)
        if self.is_aimc:
            nb_of_1b_multiplier = (
                self.weight_precision * self.wordline_dim_size * self.bitline_dim_size * self.nb_of_banks
            )
        else:
            nb_of_1b_multiplier = (
                self.bit_serial_precision
                * self.weight_precision
                * self.wordline_dim_size
                * self.bitline_dim_size
                * self.nb_of_banks
            )
        area_mults = self.cell.get_1b_multiplier_area() * nb_of_1b_multiplier

        '''
        Area of adder trees without place values
        If it is analog the output is handled by the ADC 
        If it is digital -> we need to sum #bitline numbers of size weight_precision. 
        The result of output will be weight_precision + depth_of_the_tree 
        To calculate the amount of 1b adders analytically ->  
        number_of_inputs * (adder_input_precision+1) - adder_input_precision-depth-1
        How many of those trees do we have(?)  bit_serial_precision*wordlines*nb_of_banks

        '''
        if self.is_aimc:
            adder_output_precision_regular = 0
            area_adders_regular = 0
        else:
            adder_input_precision_regular = self.weight_precision
            nb_inputs_of_adder_regular = self.bitline_dim_size  # the number of inputs of the adder tree
            adder_depth_regular = math.log2(nb_inputs_of_adder_regular)
            assert (
                adder_depth_regular % 1 == 0
            ), f"The number of inputs [{nb_inputs_of_adder_regular}] for the adder tree is not in the power of 2."
            adder_depth_regular = int(adder_depth_regular)  # float -> int for simplicity
            adder_output_precision_regular = adder_input_precision_regular + adder_depth_regular
            nb_of_1b_adder_per_tree_regular = nb_inputs_of_adder_regular * (adder_input_precision_regular + 1) - (
                adder_input_precision_regular + adder_depth_regular + 1
            )  # nb of 1b adders in a single adder tree
            nb_of_adder_trees = self.bit_serial_precision * self.wordline_dim_size * self.nb_of_banks
            area_adders_regular = self.adder.get_1b_adder_area() * nb_of_1b_adder_per_tree_regular * nb_of_adder_trees
        
        # area of adder trees with place values (type: RCA)
        if self.is_aimc:
            nb_inputs_of_adder_pv = self.weight_precision
            input_precision_pv = self.adc.adc_resolution
        else:
            nb_inputs_of_adder_pv = self.bit_serial_precision
            input_precision_pv = adder_output_precision_regular

        if nb_inputs_of_adder_pv == 1:
            nb_of_1b_adder_per_tree_pv = 0
        else:
            adder_depth_pv = math.log2(nb_inputs_of_adder_pv)
            assert (
                adder_depth_pv % 1 == 0
            ), f"The value [{nb_inputs_of_adder_pv}] of [weight_precision] is not in the power of 2."
            adder_depth_pv = int(adder_depth_pv)  # float -> int for simplicity
            nb_of_1b_adder_per_tree_pv = input_precision_pv * (nb_inputs_of_adder_pv - 1) + nb_inputs_of_adder_pv * (
                adder_depth_pv - 0.5
            )  # nb of 1b adders in a single place-value adder tree
        
        nb_of_adder_trees_pv = self.wordline_dim_size * self.nb_of_banks
        area_adders_pv = self.adder.get_1b_adder_area() * nb_of_1b_adder_per_tree_pv * nb_of_adder_trees_pv

        # area of accumulators (adder type: RCA)
        if self.bit_serial_precision == self.activation_precision:
            area_accumulators = 0
        else:
            if self.is_aimc:
                accumulator_output_precision = (
                    self.activation_precision + self.adc_resolution + self.weight_precision
                )  # output precision from adders_pv + required shifted bits
            else:
                accumulator_output_precision = (
                    self.activation_precision + math.log2(self.bitline_dim_size) + self.weight_precision
                )  # output precision from adders_pv + required shifted bits
            nb_of_1b_adder_accumulator = accumulator_output_precision * self.wordline_dim_size * self.nb_of_banks
            nb_of_1b_reg_accumulator = nb_of_1b_adder_accumulator  # number of regs in an accumulator
            area_accumulators = (
                self.adder.get_1b_adder_area() * nb_of_1b_adder_accumulator
                + self.reg.get_1b_reg_area() * nb_of_1b_reg_accumulator
            )
        # total logic area of imc macros (exclude cells)
        self.area_breakdown = {
            "cells": area_cells,
            "dacs": area_dacs,
            "adcs": area_adcs,
            "mults": area_mults,
            "adders_regular": area_adders_regular,
            "adders_pv": area_adders_pv,
            "accumulators": area_accumulators,
        }
        self.area = sum([v for v in self.area_breakdown.values()])

    def get_tclk(self):
        """! get clock cycle time of imc macros (worst path: dacs -> mults -> adcs -> adders -> accumulators)"""
        # delay of cells
        dly_cells = 0  # cells are not on critical paths

        # delay of dacs
        if self.is_aimc:
            dly_dacs = self.dac.get_dac_cost(self.bit_serial_precision)[1]
        else:
            dly_dacs = 0

        # delay of adcs
        if self.is_aimc:
            dly_adcs = self.adc.get_adc_cost(self.bitline_dim_size)[1]
        else:
            dly_adcs = 0

        # delay of multipliers
        dly_mults = self.cell.get_1b_multiplier_dly()

        # delay of regular adder trees without place value (type: RCA)
        # worst path: in-to-sum -> in-to-sum -> ... -> in-to-cout -> cin-to-cout -> ... -> cin-to-cout
        if self.is_aimc:
            adder_output_precision_regular = 0
            dly_adders_regular = 0
        else:
            adder_input_precision_regular = self.weight_precision
            nb_inputs_of_adder_regular = self.bitline_dim_size  # the number of inputs of the adder tree
            adder_depth_regular = math.log2(nb_inputs_of_adder_regular)
            assert (
                adder_depth_regular % 1 == 0
            ), f"The number of inputs [{nb_inputs_of_adder_regular}] for the adder tree is not in the power of 2."
            adder_depth_regular = int(adder_depth_regular)  # float -> int for simplicity
            adder_output_precision_regular = adder_input_precision_regular + adder_depth_regular
            dly_adders_regular = (
                (adder_depth_regular - 1) * self.adder.get_1b_adder_dly_in2sum()
                + self.adder.get_1b_adder_dly_in2cout()
                + (adder_output_precision_regular - 1 - 1) * self.adder.get_1b_adder_dly_cin2cout()
            )

        # delay of adder trees with place value (type: RCA)
        # worst path: in-to-sum -> in-to-sum -> ... -> in-to-cout -> cin-to-cout -> ... -> cin-to-cout
        if self.is_aimc:
            nb_inputs_of_adder_pv = self.weight_precision
            input_precision_pv = self.adc_resolution
        else:
            nb_inputs_of_adder_pv = self.bit_serial_precision
            input_precision_pv = adder_output_precision_regular

        if nb_inputs_of_adder_pv == 1:
            adder_pv_output_precision = input_precision_pv
            dly_adders_pv = 0
        else:
            adder_depth_pv = math.log2(nb_inputs_of_adder_pv)
            adder_depth_pv = int(adder_depth_pv)  # float -> int for simplicity
            adder_pv_output_precision = nb_inputs_of_adder_pv + input_precision_pv  # output precision from adders_pv
            dly_adders_pv = (
                (adder_depth_pv - 1) * self.adder.get_1b_adder_dly_in2sum()
                + self.adder.get_1b_adder_dly_in2cout()
                + (adder_pv_output_precision - input_precision_pv - 1) * self.adder.get_1b_adder_dly_cin2cout()
            )

        # delay of accumulators (adder type: RCA)
        if self.bit_serial_precision == self.activation_precision:
            dly_accumulators = 0
        else:
            accumulator_input_precision = adder_pv_output_precision
            if self.is_aimc:
                accumulator_output_precision = (
                    self.activation_precision + self.adc_resolution + self.weight_precision
                )  # output precision from adders_pv + required shifted bits
            else:
                accumulator_output_precision = (
                    self.activation_precision + math.log2(self.bitline_dim_size) + self.weight_precision
                )  # output precision from adders_pv + required shifted bits
            assert accumulator_input_precision < accumulator_output_precision, (
                f"accumulator_input_precision {accumulator_input_precision} must be smaller than "
                f"accumulator_output_precision {accumulator_output_precision}"
            )
            dly_accumulators = (
                self.adder.get_1b_adder_dly_in2cout()
                + (accumulator_output_precision - accumulator_input_precision - 1) * self.adder.get_1b_adder_dly_cin2cout()
            )

        # total delay of imc
        self.tclk_breakdown = {
            "cells": dly_cells,
            "dacs": dly_dacs,
            "adcs": dly_adcs,
            "mults": dly_mults,
            "adders_regular": dly_adders_regular,
            "adders_pv": dly_adders_pv,
            "accumulators": dly_accumulators,
        }
        self.tclk = sum([v for v in self.tclk_breakdown.values()])


    def get_peak_energy_single_cycle(self) -> dict[str, float]:
        """! macro-level one-cycle energy of imc arrays (fully utilization, no weight updating)
        (components: cells, mults, adders, adders_pv, accumulators. Not include input/output regs)
        """
        # energy of local bitline precharging during weight stationary in cells
        energy_local_bl_precharging = 0

        # energy of DACs
        if self.is_aimc:
            energy_dacs = self.dac.get_dac_cost(self.bit_serial_precision)[2] * self.bitline_dim_size * self.nb_of_banks
        else:
            energy_dacs = 0

        # energy of ADCs
        if self.is_aimc:
            energy_adcs = self.adc.get_adc_cost(self.bitline_dim_size)[2] * self.weight_precision * self.wordline_dim_size * self.nb_of_banks
        else:
            energy_adcs = 0

        # energy of multiplier array
        if self.is_aimc:
            nb_of_1b_multiplier = (
                self.weight_precision * self.wordline_dim_size * self.bitline_dim_size * self.nb_of_banks
            )
        else:
            nb_of_1b_multiplier = (
                self.bit_serial_precision
                * self.weight_precision
                * self.wordline_dim_size
                * self.bitline_dim_size
                * self.nb_of_banks
            )
        energy_mults = self.cell.get_1b_multiplier_energy() * nb_of_1b_multiplier

        # energy of analog bitline addition, type: voltage-based
        if self.is_aimc:
            energy_analog_bl_addition = (
                (self.tech_param["bl_cap"] * (self.tech_param["vdd"] ** 2) * self.weight_precision)
                * self.wordline_dim_size
                * self.bitline_dim_size
                * self.nb_of_banks
            )
        else:
            energy_analog_bl_addition = 0

        # energy of regular adder trees without place values (type: RCA)
        if self.is_aimc:
            adder_output_precision_regular = 0
            energy_adders_regular = 0
        else:
            adder_input_precision_regular = self.weight_precision
            nb_inputs_of_adder_regular = self.bitline_dim_size  # the number of inputs of the adder tree
            adder_depth_regular = math.log2(nb_inputs_of_adder_regular)
            assert (
                adder_depth_regular % 1 == 0
            ), f"The number of inputs [{nb_inputs_of_adder_regular}] for the adder tree is not in the power of 2."
            adder_depth_regular = int(adder_depth_regular)  # float -> int for simplicity
            adder_output_precision_regular = adder_input_precision_regular + adder_depth_regular
            nb_of_1b_adder_per_tree_regular = nb_inputs_of_adder_regular * (adder_input_precision_regular + 1) - (
                adder_input_precision_regular + adder_depth_regular + 1
            )  # nb of 1b adders in a single adder tree
            nb_of_adder_trees = self.bit_serial_precision * self.wordline_dim_size * self.nb_of_banks
            energy_adders_regular = self.adder.get_1b_adder_energy() * nb_of_1b_adder_per_tree_regular * nb_of_adder_trees

        # energy of adder trees with place values (type: RCA)
        if self.is_aimc:
            nb_inputs_of_adder_pv = self.weight_precision
            input_precision_pv = self.adc_resolution
        else:
            nb_inputs_of_adder_pv = self.bit_serial_precision
            input_precision_pv = adder_output_precision_regular

        if nb_inputs_of_adder_pv == 1:
            nb_of_1b_adder_per_tree_pv = 0
        else:
            nb_of_1b_adder_per_tree_pv = input_precision_pv * (nb_inputs_of_adder_pv - 1) + nb_inputs_of_adder_pv * (
                math.log2(nb_inputs_of_adder_pv) - 0.5
            )
        nb_of_adder_trees_pv = self.wordline_dim_size * self.nb_of_banks
        energy_adders_pv = self.adder.get_1b_adder_energy() * nb_of_1b_adder_per_tree_pv * nb_of_adder_trees_pv

        # energy of accumulators (adder type: RCA)
        if self.bit_serial_precision == self.activation_precision:
            energy_accumulators = 0
        else:
            if self.is_aimc:
                accumulator_output_precision = (
                    self.activation_precision + self.adc_resolution + self.weight_precision
                )  # output precision from adders_pv + required shifted bits
            else:
                accumulator_output_precision = (
                    self.activation_precision + math.log2(self.bitline_dim_size) + self.weight_precision
                )  # output precision from adders_pv + required shifted bits
            nb_of_1b_adder_accumulator = accumulator_output_precision * self.wordline_dim_size * self.nb_of_banks
            nb_of_1b_reg_accumulator = nb_of_1b_adder_accumulator  # number of regs in an accumulator
            energy_accumulators = (
                self.adder.get_1b_adder_energy() * nb_of_1b_adder_accumulator
                + self.reg.get_1b_reg_energy() * nb_of_1b_reg_accumulator
            )

        peak_energy_breakdown = {  # unit: pJ (the unit is borrowed from CACTI)
            "local_bl_precharging": energy_local_bl_precharging,
            "dacs": energy_dacs,
            "adcs": energy_adcs,
            "mults": energy_mults,
            "analog_bl_addition": energy_analog_bl_addition,
            "adders_regular": energy_adders_regular,
            "adders_pv": energy_adders_pv,
            "accumulators": energy_accumulators,
        }
        return peak_energy_breakdown
    
    def get_macro_level_peak_performance(self) -> tuple[float, float, float]:
        """! macro-level peak performance of imc arrays (fully utilization, no weight updating)"""
        nb_of_macs_per_cycle = (
            self.wordline_dim_size
            * self.bitline_dim_size
            / (self.activation_precision / self.bit_serial_precision)
            * self.nb_of_banks
        )

        clock_cycle_period = self.tclk  # unit: ns
        peak_energy_per_cycle = sum([v for v in self.get_peak_energy_single_cycle().values()])  # unit: pJ
        imc_area = self.area  # unit: mm^2

        tops_peak = nb_of_macs_per_cycle * 2 / clock_cycle_period / 1000
        topsw_peak = nb_of_macs_per_cycle * 2 / peak_energy_per_cycle
        topsmm2_peak = tops_peak / imc_area

        logger = _logging.getLogger(__name__)
        imc_type_info = f"Current macro-level peak performance ({'analog' if self.is_aimc else 'digital'} imc):"
        peak_performance_info = f"TOP/s: {tops_peak}, TOP/s/W: {topsw_peak}, TOP/s/mm^2: {topsmm2_peak}"
        logger.info(imc_type_info)
        logger.info(peak_performance_info)

        return tops_peak, topsw_peak, topsmm2_peak
    
    
    def get_energy_for_a_layer(self, layer: LayerNode, mapping: Mapping) -> dict[str, float]:
        # parameter extraction
        (
            mapped_rows_total_per_macro,
            _,
            mapped_cols_per_macro,
            macro_activation_times,  # normalized to only one imc macro (bank)
        ) = self.get_mapped_oa_dim(layer, self.wl_dim, self.bl_dim)
        self.mapped_rows_total_per_macro = mapped_rows_total_per_macro

        # energy of local bitline precharging during weight stationary in cells
        energy_local_bl_precharging, self.mapped_group_depth = self.get_precharge_energy(
            self.tech_param, layer, mapping
        )

        # energy of DACs
        if self.is_aimc:
            energy_dacs = (
                self.dac.get_dac_cost(self.bit_serial_precision)[2]
                * mapped_rows_total_per_macro
                * (self.activation_precision / self.bit_serial_precision)
                * macro_activation_times
            )
        else:
            energy_dacs = 0

        # energy of ADCs
        if self.is_aimc:
            energy_adcs = (
                self.adc.get_adc_cost(self.bitline_dim_size)[2]
                * self.weight_precision
                * mapped_cols_per_macro
                * (self.activation_precision / self.bit_serial_precision)
                * macro_activation_times
            )
        else:
            energy_adcs = 0

        # energy of multiplier array
        if self.is_aimc:
            nb_of_active_1b_multiplier_per_macro = (
                self.weight_precision * self.wordline_dim_size * mapped_rows_total_per_macro
            )
        else:
            nb_of_active_1b_multiplier_per_macro = (
                self.bit_serial_precision * self.weight_precision * self.wordline_dim_size * self.bitline_dim_size
            )

        energy_mults = (
            self.cell.get_1b_multiplier_energy()
            * nb_of_active_1b_multiplier_per_macro
            * (self.activation_precision / self.bit_serial_precision)
            * macro_activation_times
        )

        # energy of analog bitline addition, type: voltage-based
        if self.is_aimc:
            energy_analog_bl_addition = (
                (self.tech_param["bl_cap"] * (self.tech_param["vdd"] ** 2) * self.weight_precision)
                * mapped_cols_per_macro
                * self.bitline_dim_size
                * (self.activation_precision / self.bit_serial_precision)
                * macro_activation_times
            )
        else:
            energy_analog_bl_addition = 0

        # energy of regular adder trees without place values (type: RCA)
        if self.is_aimc:
            adder_output_precision_regular = 0
            energy_adders_regular = 0
        else:
            adder_input_precision_regular = self.weight_precision
            nb_inputs_of_adder_regular = self.bitline_dim_size  # the number of inputs of the adder tree
            adder_depth_regular = math.log2(nb_inputs_of_adder_regular)
            adder_depth_regular = int(adder_depth_regular)  # float -> int for simplicity
            adder_output_precision_regular = adder_input_precision_regular + adder_depth_regular

            nb_of_active_adder_trees_per_macro = self.bit_serial_precision * mapped_cols_per_macro
            energy_adders_per_tree_regular = self.get_regular_adder_trees_energy(
                adder_input_precision=adder_input_precision_regular,
                active_inputs_number=mapped_rows_total_per_macro,
                physical_inputs_number=self.bitline_dim_size,
            )
            energy_adders_regular = (
                energy_adders_per_tree_regular
                * nb_of_active_adder_trees_per_macro
                * (self.activation_precision / self.bit_serial_precision)
                * macro_activation_times
            )

        # energy of adder trees with place values (type: RCA)
        if self.is_aimc:
            nb_inputs_of_adder_pv = self.weight_precision
            input_precision_pv = self.adc_resolution
        else:
            nb_inputs_of_adder_pv = self.bit_serial_precision
            input_precision_pv = adder_output_precision_regular

        if nb_inputs_of_adder_pv == 1:
            nb_of_1b_adder_per_tree_pv = 0
        else:
            nb_of_1b_adder_per_tree_pv = input_precision_pv * (nb_inputs_of_adder_pv - 1) + nb_inputs_of_adder_pv * (
                math.log2(nb_inputs_of_adder_pv) - 0.5
            )
        energy_adders_pv = (
            self.adder.get_1b_adder_energy()
            * nb_of_1b_adder_per_tree_pv
            * mapped_cols_per_macro
            * (self.activation_precision / self.bit_serial_precision)
            * macro_activation_times
        )

        # energy of accumulators (adder type: RCA)
        if self.bit_serial_precision == self.activation_precision:
            energy_accumulators = 0
        else:
            if self.is_aimc:
                accumulator_output_precision = (
                    self.activation_precision + self.adc_resolution + self.weight_precision
                )  # output precision from adders_pv + required shifted bits
            else:
                accumulator_output_precision = (
                    self.activation_precision + math.log2(self.bitline_dim_size) + self.weight_precision
                )  # output precision from adders_pv + required shifted bits

            energy_accumulators = (
                (self.adder.get_1b_adder_energy() + self.reg.get_1b_reg_energy())
                * accumulator_output_precision
                * mapped_cols_per_macro
                * (self.activation_precision / self.bit_serial_precision)
                * macro_activation_times
            )

        self.energy_breakdown = {  # unit: pJ (the unit borrowed from CACTI)
            "local_bl_precharging": energy_local_bl_precharging,
            "dacs": energy_dacs,
            "adcs": energy_adcs,
            "mults": energy_mults,
            "analog_bl_addition": energy_analog_bl_addition,
            "adders_regular": energy_adders_regular,
            "adders_pv": energy_adders_pv,
            "accumulators": energy_accumulators,
        }
        self.energy = sum([v for v in self.energy_breakdown.values()])
        return self.energy_breakdown

    def __jsonrepr__(self):
        return json_repr_handler({"operational_unit: ImcArray, dimensions": self.dimension_sizes})