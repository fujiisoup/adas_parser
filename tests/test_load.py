import os
import pytest
import adas_parser

test_dir = os.sep.join((os.path.dirname(__file__), "example_data"))
test_files = [os.sep.join((test_dir, f)) for f in os.listdir(test_dir) if f[-4:] == '.dat']


@pytest.mark.parametrize(('dataname', 'path'), [
    ('pec93#c_pju#c0', 'https://open.adas.ac.uk/download/adf15/pec93][c/pec93][c_pju][c0.dat'),
    ('qef97#he_2s-s_kvi#c6', 'https://open.adas.ac.uk/download/adf12/qef97][he/qef97][he_2s-s_kvi][c6.dat'), 
    ('qef93#h_c6', 'https://open.adas.ac.uk/download/adf12/qef93][h/qef93][h_c6.dat'),
    ('rrc96#b_c1ls', 'https://open.adas.ac.uk/download/adf08/rrc96][b/rrc96][b_c1ls.dat'),
])
def _test_find_url(dataname, path):
    actual, _ = adas_parser.search_download(dataname)
    assert actual == path


@pytest.mark.parametrize('dataname', [
    'qef97#he_2s-s_kvi#c6', 
    'qef93#h_c6', 
    'rrc96#b_c1ls',
])
def _test_load(dataname):
    adas_parser.load(dataname, force_download=True)


@pytest.mark.parametrize('filename', test_files)
def test_read(filename):
    datafile = filename[filename.rfind(os.sep) + 1:]
    if 'pec' in filename:
        data = adas_parser._read_file(filename, datafile, return_xr=True)
        print(data)
        raise ValueError
    else:
        with pytest.raises(NotImplementedError) as e_info:
            adas_parser._read_file(filename, datafile, return_xr=True)
        