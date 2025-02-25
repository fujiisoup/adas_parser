import os
import pytest
import adas_parser


test_dir = os.sep.join((os.path.dirname(__file__), "example_data"))


def test_rate():
    datafile = 'qcx#h0_gyt#c6'
    filename = test_dir + '/qcx#h0_gyt#c6.dat'
    data = adas_parser._read_file(filename, datafile, return_xr=True)
    ratecoef = adas_parser.utils.calculate_rate(
        data['energy'], data['cross_section_n'], mass=1,
        temperature=1e3
    )
    assert 'n' in ratecoef.coords

    ratecoef = adas_parser.utils.calculate_rate(
        data['energy'], data['cross_section_n'], mass=1,
        temperature=[10, 1e3]
    )
    ratecoef = adas_parser.utils.calculate_rate(
        data['energy'], data['cross_section_nlm'], mass=1,
        temperature=[10, 1e3],
    )
