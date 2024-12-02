import numpy as np
import os
import requests

_default_cache_dir = os.sep.join((os.path.expanduser("~"), ".adas_data"))


try:
    import xarray as xr
    __HAS_XARRAY__ = True
except ImportError:
    __HAS_XARRAY__ = False


def _auto_decode(s):
    try:
        return int(s)
    except ValueError:
        try: 
            return float(s)
        except ValueError:
            return s.strip()


def assure_directory(dirname):
    if not os.path.exists(dirname):
        parent_dir = dirname[: dirname.rfind(os.sep)]
        assure_directory(parent_dir)
        os.mkdir(dirname)


def search_download(filename, adas_version=None):
    r'''
    Search the url and download the file.
    Returns the url and its contents
    '''
    directory1 = filename[:filename.find('_')].replace('#', '][')
    directory2 = filename.replace('#', '][')

    if adas_version is None:
        versions = range(20)
    else:
        versions = [adas_version]
    
    content = None
    for version in versions:
        url = 'https://open.adas.ac.uk/download/adf{0:s}/{1:s}/{2:s}.dat'.format(
            str(version).zfill(2), 
            directory1, directory2
        )
        response = requests.get(url)
        if b'An error has occured. Please try again or contact us if the problem' not in response.content:
            content = response.content.decode('utf-8')
            break
    
    if content is None:
        raise FileNotFoundError('ADAS data {} was not found'.format(dataname))
    return url, content
    

def load(dataname, force_download=False, return_xr=False):
    cache_dir = os.sep.join((_default_cache_dir, "openadas"))
    assure_directory(cache_dir)

    filename = os.sep.join((cache_dir, dataname + '.dat'))
    if not os.path.exists(filename) or force_download:
        url, data = search_download(dataname)
        # download the file into cache directory
        with open(filename, 'w') as f:
            f.writelines(data)

    return _read_file(filename, dataname, return_xr=return_xr)


def _read_file(filename, dataname, return_xr):
    if dataname[:3] == 'pec':
        return _read_pec(filename, return_xr)
    raise NotImplementedError(
        'Trying to read {}, but filetype {} is not supported.'.format(
            filename, dataname
        ))
        
def _read_pec(filename, return_xr):
    '''read pec file'''
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    line = lines[0]
    n_blocks = int(line[:6])
    title = line[6:].strip()
    
    data = [] if return_xr else {}
    lines = lines[1:]
    for i in range(n_blocks):
        # header
        iline = 0
        header = lines[iline]
        emission = header[:11].strip()
        ndens = int(header[11:15])
        ntemp = int(header[15:19])
        coords = {'line': emission}
        for item in header[19:].split('/')[1:]:
            key, val = item.split('=')
            coords[key.strip()] = _auto_decode(val)
        
        ne = []
        while iline < len(lines):
            iline += 1
            line = lines[iline]
            ne = ne + [float(l) for l in line.split(' ') if len(l.strip()) > 1]
            if len(ne) >= ndens:
                break
        Te = []
        while iline < len(lines):
            iline += 1
            line = lines[iline]
            Te = Te + [float(l) for l in line.split(' ') if len(l.strip()) > 1]
            if len(Te) >= ntemp:
                break
        pec = []
        while iline < len(lines):
            iline += 1
            line = lines[iline]
            pec = pec + [float(l) for l in line.split(' ') if len(l.strip()) > 1]
            if len(pec) >= ntemp * ndens:
                break
        coords['ne'] = 'ne', np.array(ne) * 1e6, {'units': 'm-3'}
        coords['Te'] = 'Te', Te, {'units': 'eV'}
        pec = np.array(pec).reshape(ndens, ntemp) * 1e-6
        if return_xr:
            data.append(xr.DataArray(
                pec, dims=['ne', 'Te'], coords=coords,
                attrs={'units': 'm3/s'}, name='photon emission coefficient'
            ))
        else:
            data[header] = (pec, ne, Te)

        lines = lines[iline + 1:]

    if return_xr:
        data = xr.concat(data, dim='index')
        stack_keys = ['line', 'TYPE']
        if data['TYPE'].ndim == 0:
            data['TYPE'] = 'index', [data['TYPE'].item()] * data.sizes['index']
        return data.set_index(index=stack_keys).unstack()

    return data

