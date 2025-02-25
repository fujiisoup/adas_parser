import os
import numpy as np
import pytest
import adas_parser

test_dir = os.sep.join((os.path.dirname(__file__), "example_data"))
test_files = [os.sep.join((test_dir, f)) for f in os.listdir(test_dir) if f[-4:] == '.dat']


@pytest.mark.parametrize(('dataname', 'path'), [
    ('pec93#c_pju#c0', 'https://open.adas.ac.uk/download/adf15/pec93][c/pec93][c_pju][c0.dat'),
    ('qef97#he_2s-s_kvi#c6', 'https://open.adas.ac.uk/download/adf12/qef97][he/qef97][he_2s-s_kvi][c6.dat'), 
    ('qef93#h_c6', 'https://open.adas.ac.uk/download/adf12/qef93][h/qef93][h_c6.dat'),
    ('rrc96#b_c1ls', 'https://open.adas.ac.uk/download/adf08/rrc96][b/rrc96][b_c1ls.dat'),
    ('rrc93##_c6ls', 'https://open.adas.ac.uk/download/adf08/rrc93][][/rrc93][][_c6ls.dat'),
])
def _test_find_url(dataname, path):
    actual, _ = adas_parser.search_download(dataname)
    assert actual == path


@pytest.mark.parametrize('dataname', [
    'szd93#c_c5',
    'rrc93##_c6ls',
    'rrc96#b_c1ls',
    'rrc96#be_c2ls',
    'qcx#h0_ory#h1',
    'qcx#h0_gyt#c6',
    #'qef97#he_2s-s_kvi#c6', 
    #'qef93#h_c6', 
])
def test_load(dataname):
    adas_parser.load(dataname, return_xr=True)


@pytest.mark.parametrize('filename', [f for f in test_files if 'pec' in f])
def test_read_pec(filename):
    datafile = filename[filename.rfind(os.sep) + 1:]
    data = adas_parser._read_file(filename, datafile, return_xr=True)
    assert 'line' in data.dims
    assert 'TYPE' in data.dims
    
    if datafile == 'pec93#c_llr#c0.dat':
        assert '1212.5 A' in data['line']
        assert 'EXCIT' in data['TYPE']
        
@pytest.mark.parametrize('filename', [f for f in test_files if 'qcx' in f])
def test_read_qcx(filename):
    datafile = filename[filename.rfind(os.sep) + 1:]
    data = adas_parser._read_file(filename, datafile, return_xr=True)

    if datafile == 'qcx#h0_gyt#c6.dat':
        assert len(data['energy']) == 27
        assert data['energy'].attrs['units'] == 'eV/amu'
        assert data['cross_section_n'].attrs['units'] == 'm^2'
        assert 1 in data['n']
        #assert np.allclose(
        #    data['cross_section_total'], data['cross_section_n'].sum('n'), 
        #    atol=1e-20, rtol=0.2
        #)


@pytest.mark.parametrize('filename', [f for f in test_files if 'rrc' in f])
def test_read_rrc(filename):
    datafile = filename[filename.rfind(os.sep) + 1:]
    data = adas_parser._read_file(filename, datafile, return_xr=True)
    print(data)

    if datafile == 'rrc96#b_c1ls.dat':
        assert '2S2 2P2(3P4.0)' in data['lower_term']


@pytest.mark.parametrize('filename', [f for f in test_files if 'szd' in f])
def test_read_rrc(filename):
    datafile = filename[filename.rfind(os.sep) + 1:]
    data = adas_parser._read_file(filename, datafile, return_xr=True)
    print(data)

    if datafile == 'szd93#c_c2.dat':
        pass
        assert len(data['Te']) == 24
