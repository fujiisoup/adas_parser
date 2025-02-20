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


def floatF(s):
    return float(s.replace('D', 'E'))

def assure_directory(dirname):
    if not os.path.exists(dirname):
        parent_dir = dirname[: dirname.rfind(os.sep)]
        assure_directory(parent_dir)
        os.mkdir(dirname)

def SLJ_to_str(S, L, J):
    return '{:s}{:s}{:s}'.format(S.strip(), 'SPDFGHIJKLMN'[int(L)], J.strip())


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
    elif dataname[:3] == 'qcx':
        return _read_qcx(filename, return_xr)
    elif dataname[:3] == 'rrc':
        return _read_rrc(filename, return_xr)
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
        if data['TYPE'].ndim != 0:
            return data.set_index(index=['line', 'TYPE']).unstack()
        else:
            return data.swap_dims(index='line').expand_dims('TYPE')            
    return data


def _read_qcx(filename, return_xr):
    '''read qcx file'''

    if not return_xr:
        raise NotImplementedError('requires xarray to return xarray dataset')

    with open(filename, 'r') as f:
        lines = f.readlines()
    
    line = lines[0]
    receiver = line[:10].replace(' ', '')
    donor = line[10:line.find('/')].replace(' ', '')

    line = lines[1]
    n_energy = int(line[:line.find('/')])
    line = lines[2]
    nmin = int(line[:line.find('/')])
    line = lines[3]
    nmax = int(line[:line.find('/')])
    
    line = lines[4]
    energy = [
        floatF(l) for l in line[:line.find('/')].split(' ') 
        if len(l.strip()) > 0
    ]
    
    cs_n, cs_nl, cs_nlm = [], [], []
    for line in lines[12:]:
        try:
            n = int(line[:4])
        except ValueError:
            break
        try:
            l = int(line[4:7])
            try: 
                m = int(line[7:10])
                cs_nlm.append(xr.DataArray(
                    [floatF(l) for l in line[11:].split(' ') if len(l.strip()) > 0],
                    dims=['energy'], coords={'energy': energy, 'n': n, 'l': l, 'm': m}
                ))
            except ValueError:
                cs_nl.append(xr.DataArray(
                    [floatF(l) for l in line[11:].split(' ') if len(l.strip()) > 0],
                    dims=['energy'], coords={'energy': energy, 'n': n, 'l': l}
                ))
        except ValueError:
            cs_n.append(xr.DataArray(
                [floatF(l) for l in line[11:].split(' ') if len(l.strip()) > 0],
                dims=['energy'], coords={'energy': energy, 'n': n}
            ))
    data = xr.Dataset({}, coords={'energy': energy})
    if len(cs_nlm) > 0:
        data['cross_section_nlm'] = xr.concat(cs_nlm, dim='nlm').set_index(nlm=['n', 'l', 'm']).unstack()
    if len(cs_nl) > 0:
        data['cross_section_nl'] = xr.concat(cs_nl, dim='nl').set_index(nl=['n', 'l']).unstack()
    if len(cs_n) > 0:
        data['cross_section_n'] = xr.concat(cs_n, dim='n')
        
    return data


def _read_rrc(filename, return_xr):
    '''read rrc file'''
    with open(filename, 'r') as f:
        lines = f.readlines()

    seq = lines[0][5:7].strip()
    nucchg = int(lines[0][21:])

    for i in range(len(lines)):
        if "PARENT TERM INDEXING" in lines[i]:
            break              
    line = lines[i]
    n_upper = int(line[line.find('NPRNTI=') + 7:])
    lines = lines[i + 4:]
    upper_term = []
    upper_energy = []
    for i in range(n_upper):
        line = lines[i]
        term = line[7:30].strip() + '(' + SLJ_to_str(line[32], line[34], line[36:40]) + ')'
        upper_term.append(term)
        upper_energy.append(floatF(lines[i][42:53]))

    for i in range(len(lines)):
        if "LS RESOLVED TERM INDEXING" in lines[i]:
            break              
    line = lines[i]
    n_lower = int(line[line.find('NTRM=') + 5:])
    lines = lines[i + 4:]
    lower_term = []
    lower_energy = []
    for i in range(n_lower):
        line = lines[i]
        term = line[7:30].strip() + '(' + SLJ_to_str(line[32], line[34], line[36:40]) + ')'
        lower_term.append(term)
        lower_energy.append(floatF(lines[i][42:53]))
            
    lines = lines[n_lower:]
    rates = []
    for PRTI in range(1, n_upper + 1):
        for i in range(len(lines)):
            if "PRTI= {}".format(PRTI) in lines[i] or "PRTI={}".format(PRTI) in lines[i]:
                break
        line = lines[i + 2]
        Te = [floatF(l) for l in line[line.find('=')+1:].split(' ') if len(l) > 0]
        lines = lines[i + 4:]
        rate1 = []
        for line in lines:
            try:
                dest = int(line[:7].strip())
            except ValueError:
                break
            r = [floatF(l) for l in line[7:].split(' ') if len(l) > 0]
            rate1.append((dest, r))
        rates.append((Te, rate1))

    if not return_xr:
        return rates, lower_term, lower_energy, upper_term, upper_energy
    # return xr
    data = []
    for i, (Te, rate1) in enumerate(rates):
        for dest, r in rate1:
            data.append(xr.DataArray(
                r, dims=['Te'], 
                coords={
                    'Te': Te, 
                    'lower_term': lower_term[dest - 1], 'lower_energy': lower_energy[dest - 1],
                    'upper_term': upper_term[i], 'upper_energy': upper_energy[i],
                }
            ))
    data = xr.concat(data, dim='index')
    data = data.set_index(index=['upper_term', 'lower_term']).unstack()
    return data

