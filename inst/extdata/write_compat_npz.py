import io as _io
import pickle as _pkl
import zipfile as _zf

import numpy as _np


def _write_compat_npz(path, sample, quality, binsize):
    """Write NPZ with numpy 1.x/2.x cross-compatible pickle payloads."""
    with _zf.ZipFile(path, "w", compression=_zf.ZIP_DEFLATED) as zf:
        for key, val in [
            ("sample", sample),
            ("quality", quality),
            ("binsize", _np.array(binsize)),
        ]:
            arr = _np.asarray(val)
            buf = _io.BytesIO()
            if arr.dtype == object:
                _np.lib.format.write_array_header_2_0(
                    buf, _np.lib.format.header_data_from_array_1_0(arr)
                )
                pkl = _pkl.dumps(arr, protocol=2)
                pkl = pkl.replace(b"numpy._core.", b"numpy.core.")
                buf.write(pkl)
            else:
                _np.lib.format.write_array(buf, arr, allow_pickle=False)
            zf.writestr(key + ".npy", buf.getvalue())
