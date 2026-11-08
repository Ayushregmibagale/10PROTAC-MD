"""XTC reader that never writes offset/lock sidecars beside an original file.

Uses MDAnalysis's XDR offset hooks; qualified against MDAnalysis 2.10.0.
"""
from MDAnalysis.coordinates.XTC import XTCReader


class ReadOnlyXTCReader(XTCReader):
    def _load_offsets(self):
        return self._read_offsets(store=False)

    def _read_offsets(self,store=False):
        return super()._read_offsets(store=False)
