import numpy as np
import matplotlib.pyplot as plt

p_array = np.array([
    0, 0.0252086, 0.0503568, 0.0754346, 0.100432, 0.12534, 0.150147, 0.174845,
    0.199424, 0.223874, 0.248186, 0.27235, 0.296357, 0.320197, 0.343862,
    0.367342, 0.390629, 0.374539, 0.358332, 0.342015, 0.325593, 0.309074,
    0.292464, 0.275771, 0.259001, 0.242161, 0.225258, 0.208298, 0.191289,
    0.174238, 0.157151, 0.140036, 0.122898, 0.105746, 0.0885862, 0.0714249,
    0.0542695, 0.0371266, 0.0200031, 0.00290586, -0.0141584, -0.0311829,
    -0.0481609, -0.0650858, -0.0819508, -0.0987494, -0.115475, -0.132121,
    -0.148681, -0.165148, -0.125572, -0.0860254, -0.046523, -0.007081,
    0.0322851, 0.0715596, 0.110727, 0.149772, 0.188679, 0.227433, 0.266019,
    0.304421, 0.342625, 0.380615, 0.418378, 0.455898, 0.493161, 0.530152,
    0.545187, 0.559974, 0.574507, 0.588782, 0.602793, 0.616535, 0.630003,
    0.643193, 0.656099, 0.668718, 0.681044, 0.693073, 0.704801, 0.716224,
    0.727338, 0.738139, 0.748624, 0.758788, 0.768629, 0.778144, 0.787328,
    0.796179, 0.804695, 0.812872, 0.820707, 0.828199, 0.835346, 0.842144,
    0.848592, 0.854688, 0.86043, 0.865817, 0.870847, 0.830301, 0.789504,
    0.748473, 0.707225, 0.665777, 0.624146, 0.582348, 0.540401, 0.498322,
    0.456128, 0.413836, 0.371463, 0.329026, 0.286543, 0.24403, 0.201504,
    0.158983, 0.116483, 0.101495, 0.0864963, 0.071493, 0.056491, 0.0414965,
    0.0265154, 0.0115536, -0.00338295, -0.0182883, -0.0331565, -0.0479818,
    -0.0627583, -0.0774802, -0.0921417, -0.106737, -0.121261, -0.135707,
    -0.15007, -0.164344, -0.178524, -0.192605, -0.20658, -0.220445, -0.234194,
    -0.247822, -0.261324, -0.274695, -0.287929, -0.301021, -0.313967,
    -0.326762, -0.353208, -0.37946, -0.405508, -0.431342, -0.456952, -0.482328,
    -0.507461, -0.53234, -0.556958, -0.581303, -0.605368, -0.629142, -0.652618,
    -0.675786, -0.698638, -0.721165, -0.743358, -0.765211, -0.786714,
    -0.800165, -0.813269, -0.826022, -0.83842, -0.850458, -0.862131, -0.873438,
    -0.884372, -0.894931, -0.905112, -0.91491, -0.924324, -0.933349, -0.941983,
    -0.950223, -0.958068, -0.965514, -0.972559, -0.979201, -0.985439, -0.99127,
    -0.996694, -1.00171, -1.00631, -1.01051, -1.01429, -1.01766, -1.02061,
    -1.02315, -1.02528, -1.02699, -1.02829, -1.02626, -1.02382, -1.02098,
    -1.01774, -1.01409, -1.01005, -1.00561, -1.00078, -0.995559, -0.989951,
    -0.983958, -0.977584, -0.970832, -0.963706, -0.956209, -0.948344,
    -0.940117, -0.93153, -0.922588, -0.898186, -0.873474, -0.848462, -0.823161,
    -0.797581, -0.771734, -0.745631, -0.719281, -0.692697, -0.665889,
    -0.638868, -0.611646, -0.584235, -0.556644, -0.528886, -0.500973,
    -0.472914, -0.444723, -0.416411, -0.387989, -0.359468, -0.330861,
    -0.302179, -0.273433, -0.244635, -0.215798, -0.186931, -0.158048,
    -0.129159, -0.100277, -0.0714114, -0.0597965, -0.0481809, -0.0365692,
    -0.0249661, -0.0133762, -0.00180417, 0.0097455, 0.0134138, 0.0170694,
    0.0207109, 0.0243369, 0.0279459, 0.0315365, 0.0351073, 0.038657, 0.0421841,
    0.0456873, 0.0491653, 0.0526166, 0.0560401, 0.0594343, 0.062798, 0.0661298,
    0.0694286, 0.0726931, 0.075922, 0.0791141, 0.0822681, 0.0853831, 0.0884576,
    0.0914907, 0.0944812, 0.0974279, 0.10033, 0.0676305, 0.0349695, 0.00235971,
    -0.0301858, -0.0626543, -0.0950328, -0.127309, -0.159469, -0.191502,
    -0.223394, -0.255133, -0.286706, -0.318103, -0.349309, -0.380313,
    -0.411103, -0.383224, -0.355248, -0.327185, -0.299048, -0.270847,
    -0.242595, -0.214302, -0.18598, -0.15764, -0.129294, -0.100952, -0.0726275,
    -0.0443301, -0.0160716, 0.0121369, 0.0402841, 0.0683591, 0.0963506,
    0.124248, 0.147139, 0.169926, 0.192599, 0.21515, 0.23757, 0.25985,
    0.281982, 0.303957, 0.325767, 0.347402, 0.368856, 0.39012, 0.411185,
    0.432044, 0.452688, 0.473111, 0.493303, 0.513258, 0.532968, 0.552426,
    0.571624, 0.590555, 0.609213, 0.627589, 0.645679, 0.663474, 0.680968,
    0.698155, 0.715029, 0.731584, 0.747813, 0.768808, 0.789455, 0.809744,
    0.787594, 0.765174, 0.742493, 0.719561, 0.696387, 0.67298, 0.649352,
    0.625512, 0.601469, 0.577233, 0.552816, 0.528226, 0.503475, 0.478572,
    0.453527, 0.428351, 0.403054, 0.377647, 0.352139, 0.326542, 0.300866,
    0.27512, 0.249316, 0.223464, 0.197574, 0.171658, 0.145724, 0.119784,
    0.0938475, 0.0679257, 0.0420285, 0.0161662, -0.00965082, -0.0354124,
    0.00160297, 0.0385438, 0.0753954, 0.112143, 0.148773, 0.18527, 0.22162,
    0.257809, 0.293822, 0.329646, 0.365267, 0.400671, 0.435844, 0.444748,
    0.453458, 0.461968, 0.470278, 0.478382, 0.486279, 0.493966, 0.501441,
    0.5087, 0.515741, 0.522562, 0.529161, 0.535535, 0.541682, 0.547601,
    0.553289, 0.558744, 0.563966, 0.568951, 0.573699, 0.578209, 0.616833,
    0.655134, 0.693096, 0.730705, 0.767948, 0.804809, 0.841275, 0.877331,
    0.912966, 0.948164, 0.982913, 1.0172, 1.05101, 1.08434, 1.11716, 1.14947,
    1.18126, 1.21251, 1.24322, 1.27337, 1.30295, 1.33194, 1.36035, 1.38816,
    1.41536, 1.44194, 1.46789, 1.4932, 1.51786, 1.53877, 1.55903, 1.57862,
    1.59754, 1.61579, 1.63335, 1.65023, 1.66641, 1.6819, 1.69668, 1.71075,
    1.72411, 1.73676, 1.74869, 1.75989, 1.77037, 1.78012, 1.78914, 1.79743,
    1.80498, 1.8118, 1.81788, 1.78304, 1.74756, 1.71146, 1.67474, 1.63743,
    1.59953, 1.56108, 1.52207, 1.48254, 1.44249, 1.40195, 1.36092, 1.31944,
    1.27751, 1.23515, 1.19239, 1.14923, 1.1057, 1.06182, 1.0176, 0.973062,
    0.928223, 0.883104, 0.837721, 0.792094, 0.746242, 0.700184, 0.653937,
    0.584316, 0.5146, 0.444818, 0.374997, 0.305166, 0.235353, 0.165585,
    0.0958903, 0.0262966, -0.0431686, -0.112478, -0.181603, -0.250519,
    -0.319196, -0.387608, -0.455729, -0.501729, -0.547437, -0.592835,
    -0.637905, -0.68263, -0.726993, -0.770977, -0.814564, -0.85774, -0.900486,
    -0.942787, -0.984627, -1.02599, -1.06686, -1.10722, -1.14706, -1.18636,
    -1.22511, -1.26329, -1.30089, -1.3379, -1.37429, -1.41007, -1.44521,
    -1.4797, -1.51353, -1.54669, -1.59564, -1.64385, -1.69131, -1.73799,
    -1.78389, -1.82899, -1.87326, -1.9167, -1.90356, -1.88968, -1.87508,
    -1.85976, -1.84373, -1.82699, -1.80955, -1.79143, -1.77263, -1.75315,
    -1.73302, -1.71223, -1.6908, -1.70252, -1.71353, -1.72383, -1.73343,
    -1.74231, -1.75048, -1.75793, -1.76467, -1.77069, -1.77599, -1.78057,
    -1.78442, -1.78756, -1.78998, -1.79569, -1.80067, -1.80493, -1.80845,
    -1.81125, -1.81331, -1.81465, -1.81526, -1.81514, -1.8143, -1.81273,
    -1.81044, -1.80744, -1.80372, -1.79928, -1.79413, -1.78828, -1.80716,
    -1.82527, -1.84262, -1.8592, -1.875, -1.87433, -1.87292, -1.87076,
    -1.86786, -1.86422, -1.85984, -1.85473, -1.84888, -1.84231, -1.83501,
    -1.827, -1.81827, -1.80883, -1.79869, -1.78785, -1.77632, -1.7641, -1.7512,
    -1.73763, -1.72339, -1.70849, -1.66144, -1.61382, -1.56565, -1.51695,
    -1.46775, -1.41805, -1.36789, -1.31728, -1.26625, -1.21481, -1.16299,
    -1.11081, -1.05829, -1.00545, -0.952316, -0.898907, -0.845245, -0.791352,
    -0.737251, -0.682963, -0.628511, -0.573917, -0.519202, -0.464389, -0.4095,
    -0.354557, -0.299583, -0.244598, -0.189625, -0.133071, -0.0765754,
    -0.0201626, 0.0361455, 0.0923267, 0.148359, 0.20422, 0.259887, 0.31534,
    0.370556, 0.425513, 0.48019, 0.534567, 0.588621, 0.642332, 0.695679,
    0.748641, 0.801198, 0.841917, 0.882219, 0.922087, 0.961507, 1.00046,
    1.03894, 1.07693, 1.11441, 1.15137, 1.1878, 1.22368, 1.259, 1.29374,
    1.3279, 1.36146, 1.39441, 1.42674, 1.45843, 1.48947, 1.51986, 1.54958,
    1.57862, 1.60698, 1.63463, 1.66157, 1.6878, 1.7133, 1.73807, 1.76209,
    1.78536, 1.80787, 1.82961, 1.81735, 1.80439, 1.79074, 1.77639, 1.76137,
    1.74567, 1.72931, 1.71229, 1.69462, 1.6763, 1.65736, 1.63778, 1.6176,
    1.5968, 1.57541, 1.55344, 1.53088, 1.50776, 1.4813, 1.45429, 1.42676,
    1.39872, 1.37017, 1.34113, 1.31161, 1.28163, 1.25119, 1.22032, 1.18902,
    1.1573, 1.12519, 1.09269, 1.05982, 1.0266, 0.993027, 0.959126, 0.92491,
    0.890392, 0.855588, 0.820512, 0.785177, 0.7496, 0.713794, 0.677774,
    0.641555, 0.605153, 0.568581, 0.531855, 0.494991, 0.458001, 0.458074,
    0.457963, 0.45767, 0.457194, 0.456537, 0.455698, 0.454679, 0.45348,
    0.452103, 0.450547, 0.448815, 0.446907, 0.444824, 0.442567, 0.440138,
    0.437538, 0.434769, 0.431831, 0.428726, 0.414576, 0.400289, 0.38587,
    0.371325, 0.356661, 0.341884, 0.327, 0.312015, 0.296936, 0.281768,
    0.266517, 0.25119, 0.235794, 0.220334, 0.204817, 0.18925, 0.173637,
    0.157987, 0.142305, 0.126597, 0.11087, 0.0951297, 0.0793832, 0.0636364,
    0.0478956, 0.0321672, 0.0164573, 0.000772263, -0.0148818, -0.0304986,
    -0.046072, -0.0826654, -0.119153, -0.15552, -0.191752, -0.227835,
    -0.263754, -0.299497, -0.335049, -0.370395, -0.405523, -0.440419,
    -0.475069, -0.50946, -0.543579, -0.577412, -0.610947, -0.644171, -0.677071,
    -0.709635, -0.74185, -0.773705, -0.805186, -0.836283, -0.866984, -0.897276,
    -0.92715, -0.931834, -0.936135, -0.940054, -0.94359, -0.946742, -0.949508,
    -0.95189, -0.953887, -0.955499, -0.956725, -0.957567, -0.958025, -0.958098,
    -0.957789, -0.957098, -0.956025, -0.954573, -0.952742, -0.950534,
    -0.947951, -0.944994, -0.941665, -0.937967, -0.933901, -0.933889,
    -0.933504, -0.932747, -0.931618, -0.93012, -0.928253, -0.926018, -0.923418,
    -0.920454, -0.917129, -0.913443, -0.9094, -0.905001, -0.90025, -0.895148,
    -0.889699, -0.883905, -0.87777, -0.871296, -0.864487, -0.857346, -0.849877,
    -0.842083, -0.833968, -0.825536, -0.816791, -0.785411, -0.753781,
    -0.721913, -0.68982, -0.657515, -0.625012, -0.592324, -0.559465, -0.526448,
    -0.493287, -0.459994, -0.426585, -0.393071, -0.359468, -0.325788,
    -0.292045, -0.258253, -0.224425, -0.190576, -0.156717, -0.122864,
    -0.0890291, -0.0552263, -0.021469, 0.0174767, 0.0563375, 0.0950982,
    0.133744, 0.172258, 0.210627, 0.248835, 0.286867, 0.324709, 0.362345,
    0.399761, 0.436943, 0.473876, 0.510546, 0.546938, 0.583039, 0.618835,
    0.654312, 0.689457, 0.724256, 0.758697, 0.792765, 0.826448, 0.859734,
    0.89261, 0.880337, 0.867737, 0.854815, 0.841577, 0.82803, 0.814179,
    0.80003, 0.78559, 0.770864, 0.75586, 0.740584, 0.725043, 0.709243,
    0.693191, 0.676894, 0.660359, 0.643593, 0.626604, 0.609398, 0.591983,
    0.574366, 0.556555, 0.538557, 0.52038, 0.502032, 0.500629, 0.499029,
    0.497233, 0.495242, 0.493057, 0.49068, 0.488111, 0.485352, 0.482404,
    0.47927, 0.475951, 0.472448, 0.468763, 0.464898, 0.460855, 0.456636,
    0.452243, 0.447679, 0.471328, 0.494743, 0.517913, 0.540829, 0.563484,
    0.585868, 0.607973, 0.629792, 0.651315, 0.672534, 0.693443, 0.714033,
    0.734296, 0.754225, 0.773814, 0.793053, 0.811938, 0.83046, 0.848613,
    0.866391, 0.883788, 0.900796, 0.917411, 0.933625, 0.949434, 0.964833,
    0.979814, 0.994375, 1.00851, 1.02221, 1.03548
])

x_array = np.array([
    0, -0.000213073, -0.000851779, -0.00191501, -0.00340151, -0.00530981,
    -0.00763833, -0.0103853, -0.0135488, -0.0171266, -0.0211167, -0.0255164,
    -0.0303233, -0.0355347, -0.0411476, -0.0471589, -0.0535656, -0.0603642,
    -0.0675512, -0.075123, -0.0830759, -0.0914057, -0.100109, -0.10918,
    -0.118616, -0.128412, -0.138563, -0.149065, -0.159912, -0.1711, -0.182623,
    -0.194477, -0.206655, -0.219152, -0.231963, -0.244503, -0.256188,
    -0.267015, -0.276981, -0.286085, -0.294324, -0.301697, -0.308203, -0.31384,
    -0.318609, -0.322508, -0.325539, -0.327702, -0.328997, -0.329426, -0.32899,
    -0.328232, -0.32769, -0.327365, -0.327257, -0.327365, -0.327689, -0.328228,
    -0.328982, -0.32995, -0.33113, -0.332523, -0.334127, -0.335941, -0.337317,
    -0.33761, -0.336822, -0.334954, -0.33201, -0.327994, -0.322908, -0.316758,
    -0.309547, -0.301281, -0.291965, -0.281605, -0.270207, -0.257778,
    -0.244324, -0.229854, -0.214375, -0.197895, -0.180423, -0.161968,
    -0.142539, -0.122146, -0.100799, -0.078508, -0.0552843, -0.0311388,
    -0.00608319, 0.0198708, 0.0467111, 0.0744251, 0.103, 0.132423, 0.162679,
    0.193267, 0.223684, 0.253917, 0.283955, 0.313787, 0.3434, 0.372784,
    0.401927, 0.43108, 0.460492, 0.490152, 0.520047, 0.550165, 0.580494,
    0.61102, 0.641731, 0.672614, 0.703658, 0.734848, 0.766173, 0.79762,
    0.829175, 0.860826, 0.89256, 0.924364, 0.956226, 0.988132, 1.02007,
    1.05203, 1.08399, 1.11595, 1.14788, 1.17979, 1.21165, 1.24345, 1.27512,
    1.30659, 1.33785, 1.36888, 1.39968, 1.43023, 1.46052, 1.49053, 1.52026,
    1.5497, 1.57883, 1.60764, 1.63613, 1.66427, 1.69207, 1.71951, 1.74657,
    1.77326, 1.79956, 1.82527, 1.85019, 1.87456, 1.89862, 1.92237, 1.94579,
    1.96887, 1.99161, 2.01399, 2.03602, 2.05767, 2.07895, 2.09985, 2.12035,
    2.14045, 2.16014, 2.17933, 2.19791, 2.21588, 2.23322, 2.24995, 2.26604,
    2.28149, 2.29631, 2.31048, 2.324, 2.33687, 2.34908, 2.36063, 2.37151,
    2.38173, 2.39127, 2.40015, 2.40835, 2.41587, 2.42272, 2.42883, 2.43417,
    2.43872, 2.4425, 2.4455, 2.44771, 2.44915, 2.44981, 2.4497, 2.4488,
    2.44714, 2.4447, 2.4415, 2.43753, 2.43286, 2.42754, 2.4216, 2.41502,
    2.40781, 2.39998, 2.39153, 2.38246, 2.37279, 2.36251, 2.35162, 2.34015,
    2.32809, 2.31544, 2.30238, 2.28906, 2.27551, 2.26171, 2.24768, 2.23342,
    2.21894, 2.20425, 2.18934, 2.17424, 2.15894, 2.14345, 2.12777, 2.11193,
    2.09591, 2.07973, 2.06339, 2.04691, 2.03028, 2.01352, 1.99663, 1.97961,
    1.96249, 1.94526, 1.92792, 1.9105, 1.89299, 1.8754, 1.85773, 1.84001,
    1.82223, 1.80439, 1.78652, 1.76861, 1.75067, 1.7323, 1.7131, 1.69308,
    1.67225, 1.65061, 1.62818, 1.60496, 1.58105, 1.5565, 1.53134, 1.50558,
    1.47923, 1.4523, 1.42479, 1.39674, 1.36813, 1.339, 1.30935, 1.27919,
    1.24854, 1.21741, 1.18581, 1.15377, 1.12128, 1.08837, 1.05504, 1.02132,
    0.987221, 0.95275, 0.917926, 0.882763, 0.847275, 0.811808, 0.776706,
    0.741982, 0.70765, 0.673721, 0.64021, 0.607128, 0.574488, 0.542302,
    0.510581, 0.479339, 0.448586, 0.418333, 0.388592, 0.359374, 0.330689,
    0.302547, 0.274959, 0.247935, 0.221209, 0.194518, 0.167872, 0.141281,
    0.114756, 0.088307, 0.0619453, 0.0356808, 0.00952379, -0.0165155,
    -0.0424268, -0.0682001, -0.0938254, -0.119027, -0.14353, -0.167326,
    -0.190407, -0.212765, -0.234392, -0.255282, -0.275428, -0.294822,
    -0.313459, -0.331333, -0.348439, -0.36477, -0.380322, -0.39509, -0.40907,
    -0.422258, -0.434651, -0.446244, -0.457035, -0.467021, -0.476199,
    -0.484569, -0.492127, -0.498873, -0.504805, -0.509924, -0.514227,
    -0.517716, -0.52039, -0.52225, -0.523297, -0.523715, -0.523689, -0.523221,
    -0.52231, -0.520958, -0.519167, -0.516938, -0.514273, -0.511174, -0.507643,
    -0.503683, -0.499295, -0.494482, -0.489247, -0.483594, -0.477524,
    -0.471427, -0.465691, -0.460317, -0.455307, -0.450662, -0.446382, -0.44247,
    -0.438925, -0.435749, -0.432942, -0.430505, -0.428437, -0.426738, -0.42541,
    -0.424451, -0.423861, -0.42364, -0.423786, -0.4243, -0.42518, -0.426425,
    -0.428035, -0.430006, -0.432339, -0.435031, -0.438081, -0.441486,
    -0.445245, -0.449356, -0.453815, -0.458621, -0.463771, -0.469262,
    -0.475092, -0.481257, -0.487487, -0.493514, -0.499336, -0.504951,
    -0.510357, -0.515551, -0.520534, -0.525302, -0.529855, -0.534191,
    -0.538308, -0.542224, -0.545955, -0.549499, -0.552855, -0.556023,
    -0.559001, -0.561789, -0.564386, -0.566792, -0.569005, -0.571026,
    -0.572853, -0.574487, -0.575928, -0.577174, -0.578226, -0.579085,
    -0.579749, -0.580219, -0.580496, -0.580579, -0.580469, -0.580167,
    -0.579672, -0.578986, -0.578109, -0.577041, -0.575785, -0.574559,
    -0.573583, -0.572858, -0.572382, -0.572157, -0.57218, -0.572452, -0.572973,
    -0.573741, -0.574755, -0.576015, -0.577521, -0.579269, -0.581261,
    -0.583493, -0.585729, -0.58773, -0.589497, -0.59103, -0.592328, -0.59339,
    -0.594218, -0.594811, -0.59517, -0.595295, -0.595185, -0.594843, -0.594268,
    -0.593461, -0.592423, -0.591155, -0.589657, -0.587931, -0.585979, -0.5838,
    -0.581396, -0.57877, -0.575922, -0.572727, -0.56906, -0.564924, -0.56032,
    -0.555253, -0.549724, -0.543737, -0.537296, -0.530403, -0.523063,
    -0.515279, -0.507055, -0.498395, -0.489304, -0.479786, -0.469847,
    -0.459489, -0.44872, -0.437543, -0.425964, -0.413989, -0.401623, -0.388872,
    -0.375741, -0.362237, -0.348366, -0.334134, -0.319581, -0.30475, -0.289645,
    -0.274275, -0.258645, -0.242762, -0.226633, -0.210266, -0.193666,
    -0.176842, -0.1598, -0.142547, -0.125092, -0.10744, -0.0896007, -0.0720047,
    -0.0550832, -0.0388415, -0.023285, -0.00841826, 0.00575397, 0.0192275,
    0.0319983, 0.0440626, 0.0554172, 0.0660588, 0.0759846, 0.0851921,
    0.0936791, 0.101444, 0.108484, 0.114799, 0.120387, 0.125248, 0.12938,
    0.132785, 0.135462, 0.137411, 0.138634, 0.13913, 0.138902, 0.138187,
    0.137224, 0.136014, 0.134557, 0.132855, 0.130908, 0.128719, 0.126288,
    0.123616, 0.120582, 0.117062, 0.113059, 0.108781, 0.104438, 0.100032,
    0.0955653, 0.0910387, 0.0864546, 0.0818149, 0.0771216, 0.0723768,
    0.0675823, 0.0627403, 0.0578526, 0.0529215, 0.0478715, 0.0426274,
    0.0371917, 0.031567, 0.0257559, 0.0197611, 0.0135853, 0.00723141,
    0.000702308, -0.00599906, -0.0128697, -0.0199064, -0.0271062, -0.0344658,
    -0.0419819, -0.0496513, -0.057343, -0.0649266, -0.072399, -0.0797577,
    -0.0869998, -0.0941228, -0.101124, -0.108001, -0.114751, -0.121371,
    -0.12786, -0.134215, -0.140434, -0.146514, -0.152454, -0.158251, -0.163903,
    -0.169408, -0.174839, -0.180268, -0.185693, -0.191111, -0.196521,
    -0.201919, -0.207305, -0.212676, -0.21803, -0.223365, -0.228678, -0.233968,
    -0.239233, -0.24447, -0.249677, -0.254853, -0.260219, -0.265996, -0.27218,
    -0.27877, -0.28576, -0.293148, -0.30093, -0.309102, -0.317659, -0.326598,
    -0.335915, -0.345605, -0.355662, -0.366083, -0.376863, -0.387996,
    -0.399478, -0.411303, -0.423466, -0.435961, -0.448782, -0.461924,
    -0.475381, -0.489147, -0.503215, -0.517581, -0.532236, -0.547175,
    -0.562391, -0.577878, -0.593628, -0.609636, -0.625894, -0.641938,
    -0.657268, -0.671841, -0.685652, -0.698698, -0.710974, -0.722478,
    -0.733206, -0.743155, -0.752324, -0.76071, -0.768311, -0.775126, -0.781154,
    -0.786393, -0.791137, -0.795677, -0.800011, -0.804139, -0.808059,
    -0.811771, -0.815272, -0.818562, -0.821641, -0.824507, -0.827159,
    -0.829598, -0.831822, -0.83383, -0.83533, -0.836025, -0.835919, -0.835011,
    -0.833306, -0.830803, -0.827508, -0.823421, -0.818547, -0.812888,
    -0.806449, -0.799235, -0.791248, -0.782495, -0.772979, -0.762707,
    -0.751684, -0.739916, -0.727409, -0.714169, -0.700204, -0.68552, -0.670392,
    -0.655093, -0.63963, -0.62401, -0.60824, -0.592324, -0.576272, -0.560088,
    -0.54378, -0.527354, -0.510818, -0.494178, -0.477441, -0.460614, -0.443704,
    -0.426718, -0.409662, -0.392544, -0.37537, -0.358148, -0.340885, -0.323587,
    -0.306261, -0.288915, -0.271556, -0.25419, -0.236824, -0.219322, -0.201548,
    -0.183509, -0.165213, -0.146667, -0.12788, -0.10886, -0.0896148,
    -0.0701521, -0.0504805, -0.0306083, -0.0105437, 0.00970477, 0.0301287,
    0.0507196, 0.0714689, 0.092368, 0.113408, 0.134581, 0.155877, 0.177288,
    0.198806, 0.220421, 0.242049, 0.263608, 0.285089, 0.306484, 0.327783,
    0.348979, 0.370064, 0.391029, 0.411867, 0.432568, 0.453125, 0.473531,
    0.493777, 0.513856, 0.533759, 0.55348, 0.57301, 0.592343, 0.61147,
    0.630386, 0.649082, 0.667552, 0.685788, 0.703784, 0.721533, 0.739012,
    0.756197, 0.773081, 0.789659, 0.805924, 0.82187, 0.837493, 0.852786,
    0.867743, 0.88236, 0.896631, 0.910552, 0.924117, 0.937322, 0.950162,
    0.962633, 0.97473, 0.98645, 0.997788, 1.00874, 1.0193, 1.02948, 1.03925,
    1.04863, 1.05761, 1.066, 1.07363, 1.08049, 1.08659, 1.09192, 1.09649,
    1.10028, 1.10331, 1.10557, 1.10706, 1.10779, 1.10775, 1.10695, 1.10539,
    1.10307, 1.09999, 1.09616, 1.09157, 1.08624, 1.08017, 1.07335, 1.0658,
    1.05752, 1.04851, 1.03877, 1.02833, 1.01719, 1.0054, 0.992963, 0.979879,
    0.966159, 0.951807, 0.936833, 0.921242, 0.905043, 0.888243, 0.87085,
    0.852871, 0.834316, 0.815194, 0.795511, 0.775279, 0.754506, 0.7332,
    0.711373, 0.689033, 0.666191, 0.642857, 0.619041, 0.594753, 0.569864,
    0.544243, 0.517903, 0.490856, 0.463113, 0.434687, 0.405592, 0.375839,
    0.345443, 0.314416, 0.282772, 0.250525, 0.217689, 0.184279, 0.150309,
    0.115794, 0.0807476, 0.0451862, 0.00912463, -0.0274217, -0.0644372,
    -0.101906, -0.139813, -0.178141, -0.216874, -0.255997, -0.295492,
    -0.335344, -0.375067, -0.414179, -0.452666, -0.490513, -0.527707,
    -0.564234, -0.60008, -0.635233, -0.66968, -0.703409, -0.736408, -0.768664,
    -0.800167, -0.830905, -0.860868, -0.890046, -0.918427, -0.946003,
    -0.972764, -0.998701, -1.02381, -1.04807, -1.07193, -1.09584, -1.11979,
    -1.14376, -1.16775, -1.19175, -1.21574, -1.23972, -1.26367, -1.2876,
    -1.31148, -1.33532, -1.35909, -1.38279, -1.40642, -1.42996, -1.4534,
    -1.47673, -1.49995, -1.52304, -1.546, -1.56882, -1.59149, -1.614, -1.63633,
    -1.65849, -1.68047, -1.70189, -1.7224, -1.74199, -1.76065, -1.77838,
    -1.79518, -1.81103, -1.82593, -1.83988, -1.85288, -1.86492, -1.876,
    -1.88612, -1.89527, -1.90345, -1.91067, -1.91692, -1.92219, -1.9265,
    -1.92984, -1.9322, -1.9336, -1.93404, -1.93348, -1.93193, -1.92937,
    -1.92581, -1.92126, -1.91572, -1.90919, -1.90168, -1.89319, -1.88373,
    -1.8733, -1.86191, -1.84957, -1.83628, -1.82206, -1.8069, -1.79081,
    -1.77382, -1.75591, -1.7371, -1.71741, -1.69709, -1.67643, -1.65541,
    -1.63406, -1.61238, -1.59038, -1.56808, -1.54548, -1.52258, -1.49941,
    -1.47597, -1.45227, -1.42831, -1.40412, -1.3797, -1.35506, -1.33021,
    -1.30516, -1.27966, -1.25343, -1.22651, -1.19889, -1.17059, -1.14162,
    -1.112, -1.08173, -1.05084, -1.01934
])

print(sum(p_array*p_array)/1000)
