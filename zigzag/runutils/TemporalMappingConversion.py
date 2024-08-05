import logging
from typing import Any

from zigzag.opt.loma.MemoryAllocator import MemoryAllocator
from zigzag.hardware.architecture.Accelerator import Accelerator
from zigzag.mapping.SpatialMappingInternal import SpatialMappingInternal
from zigzag.mapping.TemporalMapping import TemporalMapping
from zigzag.workload.layer_attributes import LayerTemporalOrdering
from zigzag.workload.layer_node import LayerNode


class TemporalMappingConversion():

    def __init__(
        self,
        *,
        accelerator: Accelerator,
        layer: LayerNode,
        spatial_mapping: SpatialMappingInternal,
        **kwargs: Any,
    ):
        
        self.layer = layer
        self.spatial_mapping = spatial_mapping
        self.accelerator = accelerator

    
    def run(self):
        temporal_mapping = self.convert_user_temporal_mapping(self.layer.temporal_ordering)
        return temporal_mapping
    
    def convert_user_temporal_mapping(self, user_temporal_mapping: LayerTemporalOrdering) -> TemporalMapping:
        """!
        # TODO move to `LayerTemporalOrdering`, fix types.
        """
        spatial_mapping = self.spatial_mapping
        layer = self.layer
        layer_dim_sizes = layer.layer_dim_sizes
        for i, utm in list(enumerate(user_temporal_mapping.data))[::-1]:
            if utm[0] not in layer_dim_sizes.layer_dims:
                logger.warning(
                    f"Supplied temporal ordering {utm} for layer {layer} thrown out because loop not present in the layer"
                )
                del user_temporal_mapping[i]

        converted_mapping = []
        for dim, size in user_temporal_mapping:
            if size == "all":
                size = layer_dim_sizes[dim]
                size_already = 1
                for dim_already, size_already_sub in converted_mapping + spatial_mapping.spatial_loop_dim_size:
                    if dim_already == dim:
                        size_already *= size_already_sub
                size //= size_already
            converted_mapping.append((dim, size))
        allocator = MemoryAllocator(self.accelerator, layer, spatial_mapping, converted_mapping)

        temporal_mapping = allocator.run()  # allocate this ordering to the memories
        return temporal_mapping
