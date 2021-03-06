import os
import tempfile
from itertools import chain
from pathlib import Path

from nose.tools import ok_
import numpy as np
from numpy.testing import assert_allclose, assert_array_almost_equal
from tqdm import tqdm

from morphio import Morphology, Option, RawDataError, set_maximum_warnings

from morph_tool import diff
from morph_tool.converter import (contour2centroid, contourcenter, get_sides,
                                  make_convex, convert)
from morph_tool.neuron_surface import get_NEURON_surface

_path = os.path.dirname(os.path.abspath(__file__))


def _get_surface(filename, extension):
    if extension == 'h5':
        raise NotImplementedError
    if extension == 'asc':
        return get_NEURON_surface(filename)
    if extension == 'swc':
        return get_NEURON_surface(filename)
    raise NotImplementedError


def _walk(folder):
    path = Path(folder)
    return chain(path.rglob('*.asc'), path.rglob('*.ASC'))


def assert_conversion_works(input_file):
    ext_in = input_file.lower().split('.')[-1]
    name = os.path.basename(input_file).split('.')[0]
    for ext_out in ['asc', 'h5', 'swc']:
        output_file = os.path.join('/tmp', name + '.' + ext_out)
        convert(input_file, output_file)
        input = Morphology(input_file)
        output = Morphology(output_file)
        diff_result = diff(input, output)
        ok_(not diff_result,
            'Difference between {} and {}: {}'.format(input_file, output_file, diff_result.info))
        try:
            if ext_in == ext_out or {ext_in, ext_out} == {'asc', 'h5'}:
                assert_array_almost_equal(Morphology(input_file).soma.points,
                                          Morphology(output_file).soma.points)
            else:
                surf1 = _get_surface(input_file, ext_in)
                surf2 = _get_surface(output_file, ext_out)
                try:
                    assert_allclose(surf1, surf2, 5e-1)
                except Exception as e:
                    raise e
        except NotImplementedError:
            pass


def test_run_converter():
    for ext_in in ['asc', 'swc', 'h5']:
        assert_conversion_works(os.path.join(_path, 'circle_contour.' + ext_in))

    # # real_neuron.asc prints a lot of warning that we do not care about
    set_maximum_warnings(0)
    assert_conversion_works(os.path.join(_path, 'real_neuron.asc'))
    assert_conversion_works(os.path.join(_path, 'real_neuron2.asc'))


def test_sides():
    points = np.array(
        [[-1.7688164200706389e+00, -1.4035942635744076e+01, 1.6111242490058810e+01],
         [-1.9805765734067151e+00, -1.3413675329177915e+01, 1.6111242490058810e+01],
         [-2.1923367267427913e+00, -1.2791408022611753e+01, 1.6111242490058810e+01],
         [-2.4040968800788676e+00, -1.2169140716045591e+01, 1.6111242490058810e+01],
         [-2.6158570334149438e+00, -1.1546873409479430e+01, 1.6111242490058810e+01],
         [-2.8276171867510200e+00, -1.0924606102913271e+01, 1.6111242490058810e+01],
         [-3.0393773400870980e+00, -1.0302338796347110e+01, 1.6111242490058810e+01],
         [-3.2511374934231743e+00, -9.6800714897809481e+00, 1.6111242490058810e+01],
         [-3.4628976467592505e+00, -9.0578041832147864e+00, 1.6111242490058810e+01],
         [-3.7017396160667762e+00, -8.4459907311487985e+00, 1.6111242490058810e+01],
         [-3.9630058061317470e+00, -7.8428332547871467e+00, 1.6111242490058810e+01],
         [-4.2242719961967197e+00, -7.2396757784254930e+00, 1.6111242490058810e+01],
         [-4.4855381862616905e+00, -6.6365183020638412e+00, 1.6111242490058810e+01],
         [-4.7468043763266632e+00, -6.0333608257021876e+00, 1.6111242490058810e+01],
         [-5.0080705663916341e+00, -5.4302033493405357e+00, 1.6111242490058810e+01],
         [-5.2693367564566049e+00, -4.8270458729788839e+00, 1.6111242490058810e+01],
         [-5.5306029465215776e+00, -4.2238883966172303e+00, 1.6111242490058810e+01],
         [-5.7918691365865484e+00, -3.6207309202555766e+00, 1.6111242490058810e+01],
         [-6.0531353266515211e+00, -3.0175734438939248e+00, 1.6111242490058810e+01],
         [-6.1046479396507927e+00, -2.3740798099891141e+00, 1.6111242490058810e+01],
         [-6.0857208749253520e+00, -1.7170404431450077e+00, 1.6111242490058810e+01],
         [-6.0667938101999113e+00, -1.0600010763009031e+00, 1.6111242490058810e+01],
         [-6.0478667454744706e+00, -4.0296170945679677e-01, 1.6111242490058810e+01],
         [-6.0289396807490299e+00,  2.5407765738730781e-01, 1.6111242490058810e+01],
         [-6.0100126160235892e+00,  9.1111702423141416e-01, 1.6111242490058810e+01],
         [-5.9910855512981485e+00,  1.5681563910755196e+00, 1.6111242490058810e+01],
         [-5.9721584865727078e+00,  2.2251957579196251e+00, 1.6111242490058810e+01],
         [-5.9532314218472671e+00,  2.8822351247637270e+00, 1.6111242490058810e+01],
         [-5.9343043571218264e+00,  3.5392744916078325e+00, 1.6111242490058810e+01],
         [-5.9153772923963857e+00,  4.1963138584519379e+00, 1.6111242490058810e+01],
         [-5.7261235238967796e+00,  4.8127359022788401e+00, 1.6111242490058810e+01],
         [-5.4465103123668825e+00,  5.4076101875694054e+00, 1.6111242490058810e+01],
         [-5.1668971008369855e+00,  6.0024844728599716e+00, 1.6111242490058810e+01],
         [-4.8872838893070885e+00,  6.5973587581505377e+00, 1.6111242490058810e+01],
         [-4.6076706777771932e+00,  7.1922330434411030e+00, 1.6111242490058810e+01],
         [-4.3280574662472961e+00,  7.7871073287316683e+00, 1.6111242490058810e+01],
         [-4.0484442547173991e+00,  8.3819816140222336e+00, 1.6111242490058810e+01],
         [-3.7688310431875021e+00,  8.9768558993127989e+00, 1.6111242490058810e+01],
         [-3.4892178316576086e+00,  9.5717301846033624e+00, 1.6111242490058810e+01],
         [-3.1613982024135154e+00,  1.0140192740314404e+01, 1.6111242490058810e+01],
         [-2.8105466055139026e+00,  1.0696036351697153e+01, 1.6111242490058810e+01],
         [-2.4596950086142897e+00,  1.1251879963079901e+01, 1.6111242490058810e+01],
         [-2.1088434117146786e+00,  1.1807723574462651e+01, 1.6111242490058810e+01],
         [-1.7579918148150657e+00,  1.2363567185845399e+01, 1.6111242490058810e+01],
         [-1.3167999933561259e+00,  1.2841068364987382e+01, 1.0043666102714512e+01],
         [-8.1621047312024508e-01,  1.3267060261145087e+01, -1.3274671529877935e-02],
         [-3.1562095288436609e-01,  1.3693052157302789e+01, -1.0070215445774267e+01],
         [1.8496856735151468e-01,  1.4119044053460494e+01, -2.0127156220018659e+01],
         [6.8555808758739367e-01,  1.4545035949618196e+01, -3.0184096994263047e+01],
         [1.1279120542519880e+00,  1.4483170152149899e+01, -3.5118757052177514e+01],
         [1.5141635392149047e+00,  1.3951316080316314e+01, -3.5118757052177514e+01],
         [1.9004150241778213e+00,  1.3419462008482725e+01, -3.5118757052177514e+01],
         [2.2866665091407370e+00,  1.2887607936649140e+01, -3.5118757052177514e+01],
         [2.6729179941036536e+00,  1.2355753864815551e+01, -3.5118757052177514e+01],
         [2.9944306329982515e+00,  1.1786229474791574e+01, -3.5118757052177514e+01],
         [3.2659706121101779e+00,  1.1187626932944383e+01, -3.5118757052177514e+01],
         [3.5375105912221052e+00,  1.0589024391097190e+01, -3.5118757052177514e+01],
         [3.8090505703340316e+00,  9.9904218492499997e+00, -3.5118757052177514e+01],
         [4.0805905494459580e+00,  9.3918193074028089e+00, -3.5118757052177514e+01],
         [4.3521305285578853e+00,  8.7932167655556182e+00, -3.5118757052177514e+01],
         [4.6236705076698126e+00,  8.1946142237084274e+00, -3.5118757052177514e+01],
         [4.8952104867817390e+00,  7.5960116818612349e+00, -3.5118757052177514e+01],
         [5.1667504658936654e+00,  6.9974091400140441e+00, -3.5118757052177514e+01],
         [5.4382904450055927e+00,  6.3988065981668525e+00, -3.5118757052177514e+01],
         [5.7098304241175200e+00,  5.8002040563196617e+00, -3.5118757052177514e+01],
         [5.9749870117084587e+00,  5.1996941203991636e+00, -3.5008457815034504e+01],
         [6.0763711363460793e+00,  4.5502480313097866e+00, -3.2068318410733596e+01],
         [6.1777552609836981e+00,  3.9008019422204088e+00, -2.9128179006432678e+01],
         [6.2791393856213187e+00,  3.2513558531310309e+00, -2.6188039602131770e+01],
         [6.3805235102589375e+00,  2.6019097640416540e+00, -2.3247900197830859e+01],
         [6.4819076348965581e+00,  1.9524636749522752e+00, -2.0307760793529948e+01],
         [6.5832917595341769e+00,  1.3030175858628983e+00, -1.7367621389229040e+01],
         [6.6052532110057420e+00,  6.6238986975010050e-01, -1.5108756823295682e+01],
         [6.3638773835365630e+00,  5.1000735666420383e-02, -1.5108756823295682e+01],
         [6.1225015560673830e+00, -5.6038839841726151e-01, -1.5108756823295682e+01],
         [5.8811257285982057e+00, -1.1717775325009363e+00, -1.5108756823295682e+01],
         [5.6397499011290266e+00, -1.7831666665846164e+00, -1.5108756823295682e+01],
         [5.3983740736598467e+00, -2.3945558006682983e+00, -1.5108756823295682e+01],
         [5.1569982461906676e+00, -3.0059449347519784e+00, -1.5108756823295682e+01],
         [4.8916435100596036e+00, -3.6060313556985477e+00, -1.4099929155922805e+01],
         [4.5736686561140907e+00, -4.1813147255436647e+00, -1.0877296370887610e+01],
         [4.2556938021685768e+00, -4.7565980953887816e+00, -7.6546635858524148e+00],
         [3.9377189482230639e+00, -5.3318814652338968e+00, -4.4320308008172180e+00],
         [3.6197440942775509e+00, -5.9071648350790138e+00, -1.2093980157820230e+00],
         [3.3017692403320371e+00, -6.4824482049241308e+00, 2.0132347692531738e+00],
         [2.9837943863865242e+00, -7.0577315747692477e+00, 5.2358675542883688e+00],
         [2.6658195324410112e+00, -7.6330149446143647e+00, 8.4585003393235638e+00],
         [2.3478446784954983e+00, -8.2082983144594817e+00, 1.1681133124358759e+01],
         [2.0298698245499853e+00, -8.7835816843045986e+00, 1.4903765909393954e+01],
         [1.9150644177674048e+00, -9.4192742747462059e+00, 1.5701242642646701e+01],
         [1.8670682816299768e+00, -1.0074831544767596e+01, 1.5701242642646701e+01],
         [1.8190721454925480e+00, -1.0730388814788990e+01, 1.5701242642646701e+01],
         [1.7710760093551201e+00, -1.1385946084810383e+01, 1.5701242642646701e+01],
         [1.7230798732176922e+00, -1.2041503354831777e+01, 1.5701242642646701e+01],
         [1.6750837370802643e+00, -1.2697060624853171e+01, 1.5701242642646701e+01],
         [1.6270876009428363e+00, -1.3352617894874561e+01, 1.5701242642646701e+01],
         [1.5790914648054084e+00, -1.4008175164895954e+01, 1.5701242642646701e+01],
         [1.4344589949807371e+00, -1.4613679164675789e+01, 1.5701242642646701e+01],
         [9.2670062502220141e-01, -1.5031100296267239e+01, 1.5701242642646701e+01],
         [4.1894225506367277e-01, -1.5448521427858685e+01, 1.5701242642646701e+01],
         [-8.8816114894857634e-02, -1.5865942559450131e+01, 1.5701242642646701e+01]])

    # Value taken from NRN software
    expected = {
        'sides': [
            np.array([-19.0193, -18.8413, -18.6633, -18.4853, -18.3126, -18.1414,
                      -17.9701, -17.7989, -17.6276, -17.4564, -17.2851, -17.1139,
                      -16.1746, -12.8996, -9.62467, -6.34972, -3.07477, 0.20018,
                      3.47513, 6.75008, 10.025, 13.3, 14.4544, 14.6425, 14.8307,
                      15.0188, 15.2069, 15.395, 15.5832, 17.9028, 20.8653, 23.8279,
                      26.7904, 29.7529, 32.7154, 35.678, 35.8957, 36.0065, 36.1173,
                      36.2281, 36.3389, 36.4497, 36.5605, 36.6713, 36.7821, 36.8929,
                      37.0037, 37.1, 37.1775, 37.2549, 37.3324, 37.4099
                      ]),
            np.array([-19.1959, -19.0705, -18.9451, -18.8198, -18.6944, -18.569,
                      -18.4437, -18.3183, -18.1929, -18.0741, -17.9607, -17.8473,
                      -17.7338, -17.6204, -17.507, -17.3936, -17.2802, -17.1667,
                      -17.0533, -16.8995, -16.732, -16.5646, -16.3972, -16.2298,
                      -16.0624, -15.8949, -15.7275, -15.5601, -15.3927, -15.2253,
                      -15.0434, -14.8539, -14.6644, -14.4749, -14.2854, -14.0958,
                      -13.9063, -13.7168, -13.5273, -13.3374, -13.1474, -12.9574,
                      -12.7674, -12.5773, -6.58583, 3.2201, 13.026, 22.832, 32.6379
                      ])
        ],
        'rads': [
            np.array([-15.1552, -14.9033, -14.6515, -14.3997, -13.8624, -13.2491,
                      -12.6358, -12.0225, -11.4092, -10.7959, -10.1826, -9.56925,
                      -8.99437, -8.53634, -8.07832, -7.62029, -7.16226, -6.70424,
                      -6.24621, -5.78818, -5.33016, -4.87213, -4.37508, -3.86024,
                      -3.3454, -2.83057, -2.31573, -1.80089, -1.28605, -0.6667,
                      -0.0158241, 0.635052, 1.28593, 1.9368, 2.58768, 3.23856,
                      3.89021, 4.54189, 5.19358, 5.84526, 6.49695, 7.14863, 7.80031,
                      8.452, 9.10368, 9.75537, 10.407, 11.0454, 11.6666, 12.2878,
                      12.9089, 13.5301
                      ]),
            np.array([
                -12.9159, -12.2589, -11.602, -10.945, -10.2881, -9.63111,
                -8.97416, -8.31721, -7.66025, -7.00543, -6.35238, -5.69932,
                -5.04626, -4.39321, -3.74015, -3.0871, -2.43404, -1.78098,
                -1.12793, -0.497242, 0.12593, 0.749103, 1.37228, 1.99545,
                2.61862, 3.24179, 3.86496, 4.48814, 5.11131, 5.73448, 6.26928,
                6.75721, 7.24513, 7.73305, 8.22097, 8.70889, 9.19681, 9.68474,
                10.1727, 10.6213, 11.0512, 11.481, 11.9109, 12.3408, 12.6694,
                12.9315, 13.1936, 13.4557, 13.7178
            ])
        ]
    }

    sides, rads = get_sides(points=points,
                            major=[0.144518, 0.250647, -0.957231],
                            minor=[-0.290655, 0.956828, 0])
    assert_array_almost_equal(sides[0],
                              expected['sides'][0], decimal=4)
    assert_array_almost_equal(sides[1],
                              expected['sides'][1], decimal=4)

    assert_array_almost_equal(rads[0], expected['rads'][0], decimal=4)
    assert_array_almost_equal(rads[1], expected['rads'][1], decimal=4)

    convex_sides, convex_rads = make_convex(expected['sides'],
                                            expected['rads'])

    expected = {
        'rads': [[-15.1552, -14.9033, -14.6515, -14.3997, -13.8624, -13.2491,
                  -12.6358, -12.0225, -11.4092, -10.7959, -10.1826, -9.56925,
                  -8.99437, -8.53634, -8.07832, -7.62029, -7.16226, -6.70424,
                  -6.24621, -5.78818, -5.33016, -4.87213, -4.37508, -3.86024,
                  -3.3454, -2.83057, -2.31573, -1.80089, -1.28605, -0.6667,
                  -0.0158241, 0.635052, 1.28593, 1.9368, 2.58768, 3.23856, 3.89021,
                  4.54189, 5.19358, 5.84526, 6.49695, 7.14863, 7.80031, 8.452,
                  9.10368, 9.75537, 10.407, 11.0454, 11.6666, 12.2878, 12.9089,
                  13.5301
                  ],
                 [-12.9159, -12.2589, -11.602, -10.945, -10.2881, -9.63111,
                  -8.97416, -8.31721, -7.66025, -7.00543, -6.35238,
                  -5.69932, -5.04626, -4.39321, -3.74015, -3.0871, -2.43404,
                  -1.78098, -1.12793, -0.497242, 0.12593, 0.749103, 1.37228,
                  1.99545, 2.61862, 3.24179, 3.86496, 4.48814, 5.11131,
                  5.73448, 6.26928, 6.75721, 7.24513, 7.73305, 8.22097,
                  8.70889, 9.19681, 9.68474, 10.1727, 10.6213, 11.0512,
                  11.481, 11.9109, 12.3408, 12.6694, 12.9315, 13.1936,
                  13.4557, 13.7178
                  ]]
    }

    assert_array_almost_equal(convex_rads[0], expected['rads'][0])
    assert_array_almost_equal(convex_rads[1], expected['rads'][1])


def test_full():
    points = np.array([[-11.25, -22.78,  10.66],
                       [-13.04, -17.52,  10.66],
                       [-15.6, -11.61,  10.66],
                       [-15.39,  -4.32,  10.66],
                       [-12.88,   1.02,  10.66],
                       [-11.1,   3.84,  10.66],
                       [-8.55,   6.01, -40.57],
                       [-6.64,   3.38, -40.57],
                       [-3.51,  -3.52, -40.57],
                       [-2.82,  -7.94, -20.56],
                       [-4.49, -12.17, -20.56],
                       [-7.53, -17.67,  10.25],
                       [-7.94, -23.27,  10.25],
                       [-9.57, -24.61,  10.25]])

    xyz, diameters = contour2centroid(*contourcenter(points))

    expected_xyzd = np.array([[-12.2372, -13.524, 12.8036, 6.92794],
                              [-11.8297, -12.8172, 10.1041, 10.9596],
                              [-11.4221, -12.1103, 7.40457, 19.0127],
                              [-11.0145, -11.4034, 4.70503, 20.6648],
                              [-10.607, -10.6966, 2.0055, 20.4251],
                              [-10.1994, -9.98972, -0.694028, 20.1399],
                              [-9.79184, -9.28285, -3.39356, 19.8209],
                              [-9.38428, -8.57599, -6.09309, 19.5018],
                              [-8.97672, -7.86913, -8.79262, 19.1828],
                              [-8.56916, -7.16227, -11.4922, 18.8638],
                              [-8.16159, -6.45541, -14.1917, 18.5447],
                              [-7.75403, -5.74855, -16.8912, 18.2257],
                              [-7.34647, -5.04169, -19.5907, 16.7483],
                              [-6.9389, -4.33483, -22.2903, 14.0655],
                              [-6.53134, -3.62797, -24.9898, 13.5066],
                              [-6.12378, -2.9211, -27.6893, 12.9623],
                              [-5.71622, -2.21424, -30.3889, 12.4181],
                              [-5.30865, -1.50738, -33.0884, 11.8739],
                              [-4.90109, -0.800521, -35.7879, 11.3297],
                              [-4.49353, -0.0936601, -38.4875, 10.7354],
                              [-4.08597, 0.613201, -41.187, 5.77215]])

    assert_allclose(xyz, expected_xyzd[:, [0, 1, 2]], rtol=1e-2)
    assert_allclose(diameters, expected_xyzd[:, 3], rtol=4e-2)

def test_swc_1pt_soma_to_asc():
    input_file = os.path.join(_path, 'single-point-soma.swc')
    output_file = os.path.join(_path, 'tmp.asc')
    convert(input_file, output_file)
    assert_allclose(_get_surface(input_file, 'swc'),
                    _get_surface(output_file, 'asc'),
                    rtol=0.1)


def test_3pts_cylinder_to_asc():
    input_file = os.path.join(_path, 'soma_three_points_cylinder.swc')
    output_file = os.path.join(_path, 'test_3pts.asc')
    convert(input_file, output_file)


def test_same_conversion_as_asciitoh5():
    '''Check that the morph-tool conversion to produces the same h5 files as asciitoh5 converter.

    repo_base/02-morphology-repository-sanitized contains valid morphologies under ASC format
    repo_base/03-morphology-repository-sanitized-asciitoh5ed
        contains the expected conversion output

    Note: asciitoh5 also reorder the neurites but not 'morph-tool convert'
          so we compare the reordered files
    '''
    repo_base = '/gpfs/bbp.cscs.ch/project/proj30/jenkins/morph-tool/converter'

    with tempfile.TemporaryDirectory() as folder:
        sanitized = Path(repo_base, '02-morphology-repository-sanitized')
        for path in tqdm(list(_walk(sanitized))):
            morphio_morph = Path(folder, path.stem + '.h5')
            convert(path, str(morphio_morph))
            asciitoh5_morph = Path(repo_base, '03-morphology-repository-sanitized-asciitoh5ed', path.stem + '.h5')
            diff_result = diff(Morphology(morphio_morph, Option.nrn_order),
                               asciitoh5_morph)
            ok_(not diff_result, 'mismatch:\n{}\n{}\n{}'.format(path, asciitoh5_morph, diff_result.info))
