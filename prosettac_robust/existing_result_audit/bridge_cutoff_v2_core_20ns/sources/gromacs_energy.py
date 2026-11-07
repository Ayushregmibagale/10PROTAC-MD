"""Read named GROMACS energy columns without assuming request order."""
import io
import re
import numpy as np
from protac_qc import ValidationError


def named_energy_series(text,required):
    legends={}
    numeric=[]
    for line in text.splitlines():
        match=re.match(r'\s*@\s+s(\d+)\s+legend\s+"([^"]+)"',line)
        if match:
            index=int(match[1]);name=match[2]
            if index in legends or name in legends.values():
                raise ValidationError('Duplicate GROMACS energy column legend')
            legends[index]=name
        stripped=line.strip()
        if stripped and not stripped.startswith(('#','@')):numeric.append(stripped)
    if not numeric:raise ValidationError('Energy series has no numeric samples')
    try:values=np.loadtxt(io.StringIO('\n'.join(numeric)),ndmin=2)
    except ValueError as exc:raise ValidationError('Malformed GROMACS energy samples') from exc
    if values.shape[1]<2 or set(legends)!=set(range(values.shape[1]-1)):
        raise ValidationError('Every energy column requires an explicit indexed legend')
    if not np.isfinite(values).all() or np.any(np.diff(values[:,0])<=0):
        raise ValidationError('Energy samples are nonfinite or time is not increasing')
    by_name={name:values[:,index+1] for index,name in legends.items()}
    if set(required)-set(by_name):raise ValidationError('Required energy column is missing')
    for name in ['Temperature','Density']:
        if name in by_name and np.any(by_name[name]<=0):
            raise ValidationError('Nonpositive '+name+' in energy series')
    return values[:,0],{name:by_name[name] for name in required}
