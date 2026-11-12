"""Preserved native-PyG backend isolation; call before importing PyG."""
from contracts import ContractError as ValidationError
def configure_native_pyg():
    """Use native PyTorch CPU operators; exclude incompatible optional binaries.

    Process-local configuration only. Does not alter installed packages or files.
    This must be applied before importing torch_geometric.
    """
    import sys
    if 'torch_geometric' in sys.modules and not globals().get('_NATIVE_PYG_CONFIGURED',False):raise ValidationError('Configure native PyG before importing it')
    disabled=['pyg_lib','torch_scatter','torch_sparse','torch_cluster','torch_spline_conv','torch_frame','intel_extension_for_pytorch']
    for name in disabled:
        if name in sys.modules and sys.modules[name] is not None:raise ValidationError('Optional compiled GNN backend already loaded: '+name)
        sys.modules[name]=None
    globals()['_NATIVE_PYG_CONFIGURED']=True
    return disabled
