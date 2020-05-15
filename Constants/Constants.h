#ifndef Constants_H
#define Constants_H

#include <vector>
#include "../Configuration.h"

namespace Constants
{

const fp ePhi = 27.2113962;

const std::vector<std::vector<fp>> phi = {
{13.5984},
{24.5874, 54.4178},
{5.39171, 75.6401, 122.454},
{9.3227, 18.2111, 153.896, 217.719},
{8.29802, 25.1548, 37.9306, 259.372, 340.226},
{11.2603, 24.3832, 47.8878, 64.4935, 392.091, 489.993},
{14.5341, 29.6013, 47.4453, 77.4735, 97.8901, 552.067, 667.046},
{13.6181, 35.1211, 54.9355, 77.4135, 113.899, 138.119, 739.327, 871.41},
{17.4228, 34.9708, 62.708, 87.175, 114.249, 157.163, 185.187, 953.898, 1103.12},
{21.5645, 40.963, 63.4233, 97.19, 126.247, 157.934, 207.271, 239.097, 1195.81, 1362.2},
{5.13908, 47.2864, 71.62, 98.936, 138.404, 172.23, 208.504, 264.192, 299.856, 1465.13, 1648.7},
{7.64624, 15.0353, 80.1436, 109.265, 141.33, 186.76, 225.02, 265.924, 327.99, 367.489, 1761.8, 1962.66},
{5.98577, 18.8285, 28.4476, 119.992, 153.825, 190.49, 241.76, 284.64, 330.21, 398.65, 442.005, 2085.98, 2304.14},
{8.15168, 16.3458, 33.493, 45.1418, 166.767, 205.279, 246.57, 303.59, 351.28, 401.38, 476.273, 523.415, 2437.66, 2673.18},
{10.4867, 19.7695, 30.2026, 51.4439, 65.0251, 220.43, 263.57, 309.6, 372.31, 424.4, 479.44, 560.62, 611.741, 2816.91, 3069.84},
{10.36, 23.3379, 34.86, 47.222, 72.5945, 88.0529, 280.954, 328.794, 379.84, 447.7, 504.55, 564.41, 651.96, 706.994, 3223.78, 3494.19},
{12.9676, 23.8136, 39.8, 53.24, 67.68, 96.94, 114.201, 348.306, 400.851, 456.7, 530, 591.58, 656.3, 750.23, 809.198, 3658.34, 3946.29},
{15.7596, 27.6297, 40.735, 59.58, 74.84, 91.29, 124.41, 143.457, 422.6, 479.76, 540.4, 619, 685.5, 755.13, 855.5, 918.375, 4120.67, 4426.22},
{4.34066, 31.625, 45.8031, 60.917, 82.66, 99.44, 117.56, 154.87, 175.817, 503.67, 565.6, 631.1, 714.7, 786.3, 860.92, 967.7, 1034.54, 4610.87, 4934.05},
{6.11316, 11.8717, 50.9132, 67.2732, 84.34, 108.78, 127.21, 147.24, 188.54, 211.275, 591.6, 658.2, 728.6, 817.2, 894, 973.7, 1086.8, 1157.73, 5128.86, 5469.86},
{6.56149, 12.7998, 24.7568, 73.4894, 91.95, 110.68, 137.99, 158.08, 180.03, 225.18, 249.798, 687.36, 757.7, 833.2, 926.5, 1008.6, 1093.5, 1213.1, 1287.96, 5674.9, 6033.75},
{6.82812, 13.5755, 27.4917, 43.2672, 99.299, 119.533, 140.68, 170.5, 192.1, 215.92, 265.07, 291.5, 787.67, 864, 944.5, 1042.5, 1130.2, 1220.3, 1346.3, 1425.26, 6249.02, 6625.81},
{6.74619, 14.634, 29.3111, 46.709, 65.2816, 128.125, 150.72, 173.55, 206, 230.5, 254.8, 308.5, 336.274, 896, 977.2, 1062.9, 1165.2, 1258.9, 1354.2, 1486.7, 1569.66, 6851.31, 7246.12},
{6.76651, 16.4863, 30.959, 49.16, 69.46, 90.6349, 160.29, 184.76, 209.5, 244.5, 270.8, 296.7, 354.7, 384.163, 1011.6, 1097.2, 1188, 1294.8, 1394.5, 1495.1, 1634.1, 1721.18, 7481.86, 7894.8},
{7.43404, 15.64, 33.668, 51.21, 72.41, 95.604, 119.203, 195.5, 221.89, 248.6, 286.1, 314.4, 343.6, 402.95, 435.172, 1133.7, 1224.1, 1320.3, 1430.9, 1537.2, 1643.2, 1788.7, 1879.87, 8140.79, 8571.95},
{7.90247, 16.1992, 30.651, 54.91, 75, 98.985, 124.976, 151.06, 233.6, 262.1, 290.9, 330.8, 361, 392.2, 456.2, 489.312, 1262.7, 1357.8, 1460, 1575.6, 1687, 1798.4, 1950.4, 2045.76, 8828.19, 9277.68},
{7.88101, 17.0844, 33.5, 51.27, 79.5, 102, 128.9, 157.8, 186.14, 275.4, 305.32, 336.1, 378.5, 410, 441.1, 511.96, 546.588, 1397.2, 1504.5, 1606, 1724, 1844, 1960.8, 2119.4, 2218.88, 9544.18, 10012.1},
{7.63988, 18.1688, 35.187, 54.92, 76.06, 108, 132, 162, 193.2, 224.7, 319.5, 351.6, 384.5, 429.3, 462.8, 495.4, 571.07, 607.02, 1541, 1646, 1758, 1880, 2008.1, 2130.5, 2295.6, 2399.26, 10288.9, 10775.4},
{7.72638, 20.2924, 36.841, 57.38, 79.8, 103, 139, 166, 198, 232.2, 265.33, 367, 401, 436, 483.1, 518.7, 552.8, 632.5, 670.608, 1690.5, 1800, 1918, 2044, 2179.4, 2307.3, 2479.1, 2586.95, 11062.4, 11567.6},
{9.3942, 17.9644, 39.7233, 59.573, 82.6, 108, 133.9, 173.9, 203, 238, 274.4, 310.8, 417.6, 453.4, 490.6, 540, 577.8, 613.3, 697.5, 737.366, 1846.8, 1961, 2085, 2214, 2358, 2491.5, 2669.9, 2782, 11864.9, 12388.9},
{5.9993, 20.5151, 30.7258, 63.241, 86.01, 112.7, 140.8, 169.9, 211, 244, 280, 319, 356, 471.2, 508.8, 548.3, 599.8, 640, 677, 765.7, 807.308, 2010, 2129, 2258, 2391, 2543.9, 2683, 2868, 2984.43, 12696.6, 13239.5},
{7.89944, 15.9346, 34.0576, 45.7155, 90.5, 115.9, 144.9, 176.4, 212.5, 252.1, 286, 326, 367, 407, 527.9, 567.3, 609.1, 662.8, 706.7, 744, 837.1, 880.44, 2180.1, 2304, 2439, 2575, 2737.1, 2881.9, 3074, 3194.29, 13557.4, 14119.4},
{9.78855, 18.5892, 28.349, 50.15, 62.77, 121.19, 147, 180, 213, 247, 296, 333, 375, 418, 460, 587.6, 628.8, 672.9, 728.9, 774, 814, 911.7, 956.79, 2356.9, 2486, 2626, 2766, 2938, 3088.1, 3287, 3411.64, 14447.7, 15028.9},
{9.75239, 21.196, 31.697, 42.947, 68.3, 81.83, 155.327, 184, 219, 255, 291, 342.9, 383, 426, 473, 517, 650.5, 693.4, 739.8, 798, 845.8, 887, 989.6, 1036.36, 2540.7, 2674, 2820, 2964, 3146, 3301.8, 3507, 3636.53, 15367.5, 15968.1},
{11.8138, 21.591, 34.871, 47.782, 59.595, 87.39, 103.03, 192.61, 224, 261, 301, 338, 393, 436, 481, 530, 577, 716.3, 761, 809.8, 870, 920.8, 963, 1070.6, 1119.17, 2731.4, 2869, 3021, 3169, 3361, 3523.1, 3735, 3868.99, 16317, 16937.1},
{13.9996, 24.3598, 35.838, 50.85, 64.69, 78.49, 109.13, 125.802, 233, 268, 308, 350, 391, 446, 492, 540, 591, 640, 785, 831.6, 882.8, 945, 999, 1042, 1155, 1205.23, 2928.9, 3072, 3228, 3380, 3584, 3752, 3971, 4109.08, 17296.4, 17936.2},
{4.17713, 27.2895, 39.247, 52.2, 68.44, 82.9, 98.67, 132.79, 150.628, 277.12, 313.1, 356, 400, 443, 502, 550, 601, 654, 706, 857, 905.3, 958.9, 1024, 1080, 1125, 1242.5, 1294.57, 3133.3, 3281, 3443, 3600, 3815, 3988, 4214, 4356.86, 18305.9, 18965.5},
{5.69487, 11.0303, 42.8835, 56.28, 70.7, 88, 104, 121.21, 158.33, 177.3, 324.07, 362, 408, 454, 499, 562, 612, 665, 722, 774, 932, 982.1, 1038, 1105, 1165, 1211, 1333.4, 1387.19, 3344.7, 3497, 3664, 3830, 4053, 4232, 4465, 4612.4, 19345.6, 20025.2},
{6.21726, 12.2236, 20.5244, 60.6072, 75.35, 91.39, 110.02, 127, 145.64, 185.7, 205.814, 374.04, 414, 463, 512, 559, 624, 677, 733, 790, 847, 1010, 1061.9, 1120.2, 1190, 1253, 1300, 1427.6, 1483.12, 3562.9, 3720, 3892, 4060, 4299, 4484, 4724, 4875.73, 20415.7, 21115.5},
{6.63412, 13.13, 23.17, 34.4184, 80.348, 96.38, 112, 133.7, 153, 172.02, 214.9, 236.252, 426, 470, 520, 573, 622, 690, 745, 803, 863, 922, 1092, 1144.7, 1205.4, 1277, 1344, 1392, 1525.1, 1582.37, 3788, 3950, 4127, 4300, 4553, 4744, 4991, 5146.94, 21516.5, 22236.7},
{6.75885, 14.32, 25.04, 37.611, 50.5728, 102.069, 119.1, 136, 159.2, 180, 200.28, 246.1, 268.59, 482.5, 530, 581, 636, 688, 758, 816, 877, 940, 1000, 1176, 1230.6, 1293.7, 1368, 1439, 1488, 1625.9, 1684.97, 4020.1, 4187, 4369, 4540, 4815, 5011, 5265, 5426.06, 22648, 23388.8},
{7.09243, 16.16, 27.13, 40.33, 54.417, 68.827, 125.638, 143.6, 164.12, 186.3, 209.3, 230.28, 279.1, 302.6, 544, 591, 646, 702, 758, 829, 890, 953, 1019, 1082, 1263, 1319.6, 1385.1, 1462, 1537, 1587, 1730.1, 1790.93, 4259, 4430, 4618, 4800, 5084, 5287, 5548, 5713.19, 23810.7, 24572.2},
{7.11938, 15.26, 29.55, 41, 57, 72, 88, 150, 169, 189.9, 214, 239, 262.08, 311, 338.55, 604, 655, 713, 773, 829, 904, 968, 1032, 1102, 1166, 1354, 1411.6, 1479.5, 1559, 1638, 1689, 1838, 1900.28, 4505, 4681, 4874, 5060, 5361, 5570, 5838, 6008.39, 25004.5, 25787},
{7.3605, 16.76, 28.47, 45, 59, 76, 93, 110, 178.41, 198, 219.9, 245, 271, 295.9, 348, 376.25, 670, 723, 784, 845, 905, 981, 1048, 1115, 1187, 1253, 1447, 1506.7, 1577, 1659, 1743, 1794, 1949, 2013.04, 4758, 4939, 5136, 5330, 5647, 5861, 6137, 6311.72, 26229.9, 27033.5},
{7.4589, 18.08, 31.06, 42, 63, 80, 97, 115.1, 135, 207.51, 228, 252.1, 277, 306, 331.58, 389.3, 415.97, 739, 794, 857, 921, 984, 1061, 1131, 1202, 1274, 1344, 1544, 1604.9, 1677.6, 1763, 1851, 1903, 2063, 2129.22, 5018, 5203, 5406, 5600, 5940, 6161, 6444, 6623.26, 27487, 28312},
{8.33684, 19.43, 32.93, 46, 61, 84.1, 101, 120, 141, 159.9, 238.57, 260, 286, 311, 342, 369.1, 427, 457.5, 810, 869, 933, 1000, 1065, 1145, 1218, 1290, 1366, 1438, 1644, 1706.2, 1781.3, 1869, 1962, 2016, 2181, 2248.87, 5284, 5475, 5683, 5880, 6242, 6469, 6759, 6943.1, 28776, 29622.6},
{7.57623, 21.4844, 34.8, 49, 65, 82, 106, 125, 145.1, 167, 188, 271.46, 294, 321, 347, 381, 408.43, 469, 500.87, 885, 946, 1013, 1082, 1149, 1231, 1308, 1382, 1460, 1535, 1747, 1810.5, 1888, 1979, 2077, 2131, 2302, 2371.99, 5558, 5753, 5966, 6170, 6551, 6785, 7082, 7271.3, 30097.3, 30965.7},
{8.99382, 16.9083, 37.468, 51, 67.9, 87, 105, 130.1, 150, 173, 195, 218, 305, 329, 358, 385, 421, 452.6, 513, 546.19, 963, 1026, 1095, 1167, 1237, 1320, 1401, 1477, 1558, 1635, 1852, 1917.9, 1998, 2091, 2195, 2250, 2427, 2498.62, 5839, 6039, 6257, 6460, 6869, 7109, 7414, 7607.95, 31451.1, 32341.5},
{5.78636, 18.8704, 28.0441, 55.45, 69.3, 90, 109, 130.1, 156, 178, 201, 226, 249, 341, 368, 396, 425, 462, 497.1, 560, 593.38, 1043, 1109, 1181, 1255, 1328, 1413, 1496, 1575, 1659, 1738, 1961, 2028.5, 2111, 2207, 2317, 2373, 2555, 2628.77, 6126, 6331, 6554, 6770, 7196, 7442, 7754, 7953.14, 32837.6, 33750.3},
{7.34392, 14.6331, 30.506, 40.74, 77.03, 94, 112.9, 135, 156, 184, 208, 232, 258, 282, 379, 407, 437, 466, 506, 537, 608, 642.35, 1127, 1195, 1269, 1347, 1421, 1508, 1596, 1676, 1763, 1844, 2074, 2142.1, 2227, 2326, 2443, 2499, 2687, 2762.49, 6421, 6631, 6859, 7080, 7531, 7790, 8103, 8306.95, 34257.1, 35192.4},
{8.60839, 16.626, 25.3235, 43.804, 55, 99.51, 117, 139, 162, 185, 214, 238, 265, 292, 317, 420, 447, 479, 510, 552, 584, 657, 693.26, 1214, 1285, 1360, 1441, 1518, 1606, 1698, 1781, 1869, 1954, 2190, 2266, 2349, 2428, 2567, 2654, 2815, 2900, 6714, 6929, 7167, 7390, 7887, 8140, 8455, 8669.48, 35710, 36668.1},
{9.00966, 18.6, 27.84, 37.4155, 59.3, 69.1, 124.2, 143, 167, 191.1, 215, 245, 272, 299, 328, 354, 461, 491, 522, 555, 599, 633, 709, 746.12, 1304, 1377, 1455, 1538, 1618, 1707, 1803, 1889, 1979, 2066, 2309, 2386, 2472, 2552, 2700, 2788, 2954, 3041, 7022, 7243, 7485, 7714, 8240, 8499, 8821, 9040.83, 37196.5, 38177.6},
{10.4513, 19.1313, 29.57, 40.357, 51.52, 74.4, 87.61, 150.81, 171, 197, 220.9, 247, 279, 307, 335, 365, 393, 505, 535, 569, 601, 649, 683, 762, 800.8, 1397, 1472, 1553, 1639, 1720, 1812, 1911, 1999, 2093, 2181, 2431, 2510, 2598, 2680, 2836, 2926, 3096, 3185.5, 7337, 7563, 7811, 8044, 8601, 8867, 9196, 9421.1, 38717, 39721.4},
{12.1298, 20.975, 31.05, 42.2, 54.1, 66.703, 91.6, 105.978, 179.84, 202, 229.02, 255, 281, 314, 343, 374, 404, 434, 549, 582, 616, 650, 700, 736, 818, 857, 1493, 1571, 1653, 1742, 1826, 1919, 2023, 2113, 2209, 2300, 2556, 2637, 2726, 2811, 2975, 3068, 3243, 3333.8, 7660, 7889, 8144, 8382, 8971, 9243, 9581, 9810.37, 40271.7, 41299.7},
{3.89391, 23.1575, 33.195, 43, 56, 69.1, 82.9, 110.1, 125.61, 213.3, 233, 261, 289, 316, 352, 382, 413, 445, 476, 597, 629, 666, 700, 753, 791, 875, 916.1, 1592, 1672, 1757, 1848, 1936, 2029, 2137, 2230, 2329, 2422, 2683, 2767, 2859, 2945, 3118, 3214, 3392, 3485, 7989, 8224, 8484, 8726, 9350, 9629, 9974, 10208.8, 41861.1, 42913},
{5.21166, 10.0038, 35.8438, 47, 58, 71, 86, 101, 130.5, 146.52, 241, 267.1, 296, 325, 354, 390, 422, 455, 488, 520, 646, 679, 717, 752, 809, 846, 935, 976.62, 1695, 1776, 1864, 1958, 2047, 2142, 2256, 2349, 2452, 2547, 2814, 2901, 2994, 3081, 3266, 3363, 3546, 3640, 8326, 8565, 8831, 9077, 9739, 10023, 10376, 10616.4, 43485.4, 44561.5},
{5.5769, 11.185, 19.1773, 49.95, 61.6, 74, 88, 105, 119, 151.4, 168.77, 275, 303, 332, 364, 393, 431, 464, 498, 533, 566, 696, 731, 770, 806, 865, 906, 995, 1039.09, 1800, 1884, 1974, 2069, 2162, 2259, 2377, 2473, 2577, 2674, 2950, 3036, 3133, 3222, 3416, 3515, 3704, 3800, 8669, 8914, 9184, 9437, 10136, 10426, 10789, 11033.4, 45145, 46245.6},
{5.5386, 10.956, 20.1974, 36.906, 65.55, 77.6, 91, 106, 125, 140, 172, 192.24, 312, 340, 371, 403, 435, 472, 509, 543, 579, 613, 749, 785, 824, 862, 924, 965, 1060, 1103.5, 1908, 1994, 2087, 2185, 2280, 2378, 2500, 2600, 2706, 2806, 3087, 3176, 3274, 3366, 3570, 3672, 3865, 3963, 9020, 9269, 9545, 9803, 10542, 10840, 11210, 11459.9, 46840.3, 47965.7},
{5.4702, 10.631, 21.6237, 38.981, 57.53, 82, 97, 112, 131, 148, 162, 196, 217.02, 350, 378, 412, 445, 478, 516, 554, 590, 627, 663, 803, 840, 880, 920, 985, 1028, 1124, 1169.9, 2019, 2108, 2202, 2304, 2400, 2501, 2628, 2729, 2838, 2941, 3227, 3319, 3419, 3512, 3729, 3832, 4030, 4130, 9378, 9632, 9913, 10175, 10959, 11262, 11641, 11895.9, 48571.7, 49722.2},
{5.525, 10.783, 22.09, 40.6, 60, 84, 99, 114, 136, 152, 168, 195, 221, 243, 389, 420, 453, 489, 522, 562, 602, 638, 678, 714, 859, 896, 939, 978, 1049, 1092, 1191, 1238.42, 2134, 2224, 2321, 2425, 2525, 2627, 2758, 2861, 2974, 3078, 3371, 3465, 3567, 3662, 3891, 3997, 4198, 4302, 9742, 10002, 10288, 10555, 11384, 11694, 12082, 12341.7, 50339.6, 51515.6},
{5.577, 10.938, 22.44, 41.17, 61.7, 85, 101, 116, 138, 155, 174, 202, 229, 248, 269, 430, 462, 497, 534, 569, 609, 651, 689, 730, 767, 916, 956, 998, 1040, 1113, 1158, 1261, 1308.7, 2251, 2344, 2443, 2549, 2652, 2755, 2892, 2997, 3112, 3219, 3519, 3613, 3718, 3816, 4056, 4166, 4371, 4476, 10115, 10378, 10671, 10942, 11819, 12136, 12532, 12797.3, 52144.3, 53346.1},
{5.64371, 11.078, 23.55, 41.64, 62.7, 87, 103, 118, 141, 158, 179, 208, 237, 257, 276, 306.5, 474, 506, 543, 581, 617, 658, 702, 742, 782, 822, 976, 1016, 1060, 1103, 1180, 1226, 1332, 1381.56, 2371, 2466, 2569, 2676, 2782, 2887, 3028, 3137, 3253, 3363, 3669, 3766, 3873, 3971, 4227, 4337, 4548, 4655, 10494, 10762, 11060, 11337, 12264, 12588, 12992, 13262.9, 53986.1, 55214.2},
{5.67038, 11.24, 24.84, 42.94, 63.2, 89, 105, 120, 144, 161, 183, 213, 243, 263, 281, 311, 344.4, 518, 553, 590, 630, 667, 709, 755, 795, 838, 879, 1037, 1078, 1124, 1167, 1249, 1296, 1406, 1456.06, 2495, 2591, 2697, 2807, 2914, 3022, 3168, 3279, 3398, 3510, 3823, 3921, 4031, 4131, 4400, 4513, 4729, 4838, 10880, 11153, 11457, 11739, 12718, 13050, 13462, 13738.6, 55865.9, 57120.6},
{6.1498, 12.076, 20.54, 44.44, 64.8, 89, 106, 123, 144, 165, 183, 213, 246, 268, 288, 319, 352, 384.4, 565, 601, 639, 680, 719, 761, 810, 851, 895, 937, 1100, 1142, 1189, 1233, 1321, 1368, 1481, 1532.3, 2621, 2720, 2827, 2941, 3050, 3160, 3312, 3424, 3546, 3660, 3980, 4080, 4191, 4294, 4578, 4693, 4914, 5025, 11273, 11552, 11861, 12147, 13183, 13521, 13943, 14224.6, 57783.9, 59065.5},
{5.8638, 11.513, 21.82, 39.33, 66.5, 90, 108, 125, 143, 168, 186, 216, 250, 273, 294, 325, 358, 393, 426.6, 613, 651, 690, 732, 772, 816, 866, 909, 954, 997, 1165, 1208, 1256, 1301, 1393, 1443, 1559, 1610.4, 2750, 2852, 2961, 3078, 3189, 3300, 3458, 3573, 3698, 3814, 4139, 4242, 4355, 4460, 4760, 4877, 5103, 5217, 11673, 11957, 12272, 12563, 13658, 14003, 14434, 14721, 59739.3, 61049.7},
{5.93905, 11.647, 22.89, 41.23, 62.1, 93, 110, 127, 152, 170, 192, 224, 259, 279, 300, 332, 366, 399, 431, 464.9, 664, 702, 743, 786, 827, 872, 924, 969, 1014, 1059, 1232, 1275, 1325, 1371, 1468, 1520, 1638, 1691.7, 2882, 2987, 3098, 3217, 3331, 3445, 3607, 3725, 3852, 3970, 4303, 4407, 4523, 4629, 4945, 5066, 5296, 5412, 12081, 12370, 12690, 12986, 14144, 14495, 14936, 15228.1, 61736.6, 63073.5},
{6.0215, 11.781, 22.79, 42.52, 63.9, 95, 112, 129, 155, 173, 197, 229, 263, 284, 305, 340, 373, 408, 441, 475, 510, 715, 755, 797, 842, 885, 929, 985, 1029, 1077, 1122, 1300, 1346, 1395, 1443, 1545, 1598, 1719, 1773.6, 3018, 3125, 3238, 3359, 3476, 3592, 3760, 3880, 4009, 4131, 4469, 4576, 4693, 4802, 5135, 5258, 5494, 5611, 12495, 12790, 13116, 13417, 14639, 14998, 15448, 15745.8, 63772.4, 65136.8},
{6.1077, 11.916, 22.7, 42.42, 65.1, 96, 114, 131, 158, 177, 201, 235, 268, 290, 311, 345, 381, 415, 450, 486, 520, 555, 770, 810, 853, 899, 943, 989, 1046, 1092, 1142, 1188, 1370, 1416, 1468, 1516, 1625, 1678, 1803, 1858.5, 3157, 3265, 3381, 3505, 3624, 3742, 3916, 4038, 4170, 4294, 4639, 4748, 4866, 4978, 5329, 5455, 5695, 5815, 12918, 13217, 13548, 13855, 15146, 15511, 15971, 16274.6, 65848.2, 67241.8},
{6.18431, 12.065, 23.66, 42.41, 65.4, 98, 116, 133, 160, 180, 205, 239, 274, 295, 317, 352, 387, 424, 460, 496, 530, 570, 603, 825, 866, 911, 958, 1004, 1050, 1110, 1157, 1207, 1255, 1442, 1490, 1542, 1591, 1706, 1761, 1889, 1945.2, 3298, 3409, 3528, 3653, 3775, 3895, 4075, 4199, 4335, 4461, 4812, 4922, 5044, 5157, 5527, 5656, 5901, 6023, 13347, 13651, 13988, 14300, 15663, 16036, 16510, 16814.3, 67965.3, 69387.3},
{6.25416, 12.1792, 25.053, 43.61, 65.6, 99, 117, 135, 163, 182, 209, 244, 279, 301, 324, 360, 396, 431, 469, 505, 540, 580, 610, 651, 882, 924, 971, 1019, 1065, 1114, 1175, 1224, 1275, 1324, 1516, 1564, 1618, 1668, 1789, 1845, 1978, 2036.4, 3443, 3555, 3677, 3805, 3929, 4051, 4238, 4364, 4502, 4630, 4988, 5101, 5224, 5339, 5731, 5860, 6111, 6236, 13784, 14093, 14435, 14752, 16191, 16570, 17050, 17365.4, 70123, 71574.8},
{5.42587, 14.13, 20.9594, 45.249, 66.8, 98, 117, 136, 159, 185, 205, 238, 276, 305, 328, 361, 399, 438, 476, 520, 560, 600, 630, 670, 713, 941, 985, 1032, 1081, 1130, 1178, 1242, 1292, 1345, 1395, 1591, 1641, 1696, 1747, 1875, 1933, 2067, 2125.5, 3590, 3706, 3828, 3960, 4086, 4211, 4403, 4532, 4673, 4803, 5168, 5282, 5408, 5525, 5937, 6070, 6326, 6452, 14228, 14542, 14890, 15211, 16730, 17120, 17610, 17928, 72322.9, 73804.8},
{6.82507, 14.61, 22.55, 33.37, 68.37, 98, 118, 137, 157, 187, 209, 230, 270, 310, 334, 359, 399, 440, 481, 520, 570, 610, 650, 690, 730, 772, 1002, 1047, 1094, 1146, 1195, 1245, 1311, 1362, 1417, 1467, 1669, 1719, 1776, 1827, 1963, 2022, 2159, 2218.9, 3741, 3858, 3984, 4118, 4246, 4372, 4573, 4703, 4846, 4980, 5350, 5468, 5595, 5713, 6149, 6284, 6545, 6674, 14678, 14999, 15351, 15680, 17280, 17680, 18180, 18502.3, 74565.9, 76077.8},
{7.54957, 16.2, 23.1, 35, 48.272, 94.01, 119, 139, 159, 180, 213, 235, 262, 304, 338, 363, 396, 439, 482, 530, 570, 610, 660, 700, 750, 790, 832, 1064, 1110, 1160, 1211, 1262, 1313, 1382, 1434, 1490, 1542, 1748, 1799, 1857, 1910, 2053, 2113, 2254, 2314.7, 3898.7, 4014, 4143, 4278, 4410, 4537, 4745, 4877, 5024, 5159, 5537, 5655, 5785, 5907, 6364, 6502, 6769, 6900, 15137, 15461, 15820, 16150, 17840, 18250, 18760, 19088.5, 76852, 78394.7},
{7.86403, 16.37, 26, 38.2, 51.6, 64.77, 122.01, 141.2, 160.2, 179, 208.9, 231.6, 258.3, 290.7, 325.3, 361.9, 387.9, 420.7, 462.1, 502.6, 543.4, 594.5, 640.6, 685.6, 734.1, 784.4, 833.4, 881.4, 1132.2, 1180, 1230.4, 1283.4, 1335.1, 1386.8, 1459.9, 1512.4, 1569.1, 1621.7, 1829.8, 1882.9, 1940.6, 1994.8, 2149.1, 2210, 2354.5, 2414.1, 4057, 4180, 4309, 4446, 4578, 4709, 4927, 5063, 5209, 5348, 5719, 5840, 5970, 6093, 6596, 6735, 7000, 7130, 15566, 15896, 16252, 16588, 18476, 18872, 19362, 19686.7, 79181.9, 80755.6},
{7.83352, 16.6, 27, 39.1, 51.9, 67, 82.71, 144.4, 165, 187, 208, 236, 268, 291, 330, 377, 403, 429, 476, 520, 570, 620, 670, 720, 760, 810, 860, 910, 953, 1194, 1242, 1294, 1349, 1402, 1454, 1530, 1583, 1641, 1696, 1912, 1966, 2025, 2080, 2240, 2302, 2450, 2514.5, 4214, 4335, 4468, 4609, 4745, 4877, 5099, 5236, 5388, 5528, 5919, 6042, 6176, 6300, 6810, 6952, 7230, 7366, 16080, 16410, 16780, 17120, 19000, 19420, 19950, 20297.4, 81556.9, 83162.3},
{8.43823, 17, 25, 41, 55, 70.1, 85.1, 102.02, 168.7, 190, 213, 235, 269, 298, 322, 367, 410, 436, 470, 520, 570, 620, 670, 720, 770, 820, 870, 920, 970, 1015, 1262, 1311, 1364, 1420, 1474, 1528, 1606, 1660, 1720, 1776, 1996, 2052, 2112, 2168, 2336, 2400, 2552, 2615.5, 4374, 4501, 4635, 4779, 4917, 5052, 5280, 5421, 5575, 5717, 6115, 6240, 6376, 6503, 7039, 7185, 7468, 7610, 16560, 16900, 17270, 17620, 19600, 20030, 20570, 20920.6, 83976.2, 85614.4},
{8.96702, 17, 28, 40, 57, 72, 89, 105, 122.7, 194.8, 217, 240, 264, 303, 329, 356, 407, 445, 472, 510, 560, 610, 670, 720, 770, 820, 870, 920, 980, 1030, 1080, 1331, 1381, 1436, 1493, 1548, 1603, 1684, 1739, 1801, 1857, 2083, 2139, 2201, 2258, 2435, 2500, 2656, 2720.4, 4540, 4668, 4806, 4952, 5092, 5229, 5466, 5609, 5765, 5910, 6315, 6441, 6580, 6708, 7274, 7421, 7710, 7850, 17040, 17390, 17770, 18120, 20210, 20650, 21200, 21556.6, 86438.9, 88113.3},
{8.95883, 18.56, 29, 43, 56, 75, 91, 109, 126, 144.9, 220.4, 245, 269, 293, 332, 358, 392, 445, 479, 507, 550, 610, 660, 710, 760, 820, 870, 930, 980, 1040, 1090, 1140, 1402, 1454, 1509, 1567, 1624, 1680, 1763, 1821, 1883, 1941, 2171, 2228, 2291, 2350, 2536, 2603, 2762, 2827.8, 4715, 4839, 4980, 5128, 5270, 5410, 5654, 5800, 5959, 6106, 6517, 6646, 6787, 6918, 7512, 7660, 7960, 8100, 17540, 17890, 18280, 18630, 20840, 21280, 21840, 22205.7, 88955.2, 90659.7},
{9.22555, 20.203, 30, 45, 60, 74, 94, 112, 130.1, 149, 168.2, 248, 275, 299, 324, 365, 392, 433, 487, 520, 550, 600, 650, 710, 760, 820, 870, 930, 990, 1040, 1100, 1150, 1210, 1475, 1527, 1584, 1644, 1702, 1758, 1845, 1904, 1967, 2026, 2261, 2320, 2383, 2443, 2640, 2708, 2870, 2941, 4888, 5013, 5156, 5307, 5452, 5594, 5846, 5994, 6156, 6305, 6724, 6854, 6997, 7130, 7760, 7910, 8210, 8360, 18040, 18400, 18790, 19150, 21470, 21920, 22500, 22868.1, 91515.8, 93254.3},
{10.4375, 18.7569, 34.46, 48.55, 61.2, 76.6, 93, 113.9, 134, 153, 173, 192.7, 276.9, 307, 332, 357, 402, 429, 477, 530, 560, 590, 650, 710, 760, 820, 880, 930, 990, 1050, 1110, 1160, 1220, 1280, 1549, 1603, 1661, 1723, 1780, 1839, 1928, 1989, 2052, 2113, 2354, 2412, 2478, 2539, 2745, 2815, 2981, 3049.9, 5055, 5191, 5335, 5490, 5636, 5780, 6041, 6192, 6356, 6508, 6933, 7066, 7211, 7350, 8010, 8160, 8470, 8620, 18550, 18910, 19310, 19680, 22120, 22580, 23170, 23544.1, 94124.7, 95897.7},
{6.10829, 20.4283, 29.852, 51.14, 62.6, 80, 97.9, 116, 135, 158, 177, 198, 218.3, 306.9, 340, 366, 392, 439, 467, 520, 570, 600, 640, 700, 760, 820, 880, 930, 990, 1060, 1110, 1170, 1230, 1290, 1350, 1625, 1681, 1740, 1802, 1862, 1920, 2014, 2075, 2140, 2202, 2447, 2508, 2574, 2635, 2854, 2925, 3094, 3164.7, 5234, 5371, 5518, 5674, 5824, 5969, 6241, 6392, 6560, 6714, 7146, 7281, 7430, 7570, 8260, 8420, 8730, 8880, 19070, 19440, 19840, 20210, 22780, 23250, 23850, 24234.1, 96783.2, 98591.6},
{7.41668, 15.0325, 31.9373, 42.3326, 68.8, 82.9, 100.1, 120, 138, 158, 182, 203, 224, 245.1, 338.1, 374, 401, 427, 478, 507, 570, 610, 650, 690, 750, 810, 870, 930, 990, 1050, 1120, 1180, 1240, 1300, 1360, 1430, 1704, 1760, 1819, 1884, 1945, 2004, 2101, 2163, 2230, 2292, 2543, 2605, 2671, 2735, 2965, 3036, 3211, 3282.1, 5414, 5555, 5703, 5862, 6015, 6162, 6442, 6597, 6767, 6924, 7362, 7500, 7650, 7790, 8520, 8680, 9000, 9150, 19590, 19970, 20380, 20750, 23460, 23940, 24550, 24938.2, 99491.9, 101336},
{7.28552, 16.703, 25.563, 45.37, 54.856, 88.4, 103, 122, 143, 161.1, 183, 208, 229, 252, 272.6, 370.2, 409, 436, 464, 520, 550, 620, 660, 690, 750, 810, 870, 930, 990, 1060, 1120, 1180, 1250, 1310, 1380, 1440, 1500, 1784, 1840, 1902, 1967, 2029, 2090, 2190, 2253, 2321, 2385, 2641, 2703, 2771, 2835, 3078, 3151, 3329, 3401.8, 5599, 5740, 5892, 6054, 6208, 6358, 6648, 6804, 6977, 7137, 7580, 7720, 7870, 8010, 8780, 8950, 9270, 9430, 20130, 20500, 20920, 21300, 24150, 24640, 25260, 25656.9, 102252, 104133},
{8.414, 19.3, 27.3, 36, 57, 69.1, 108, 125, 146.1, 166, 186, 209, 235, 257, 281, 304, 416, 444, 473, 502, 560, 590, 670, 700, 740, 800, 870, 930, 990, 1060, 1120, 1180, 1250, 1320, 1380, 1440, 1510, 1570, 1865, 1923, 1986, 2052, 2115, 2177, 2281, 2345, 2414, 2480, 2740, 2803, 2873, 2938, 3194, 3268, 3450, 3524.2, 5785, 5930, 6084, 6248, 6405, 6557, 6856, 7015, 7191, 7350, 7810, 7950, 8100, 8240, 9050, 9220, 9550, 9710, 20670, 21050, 21470, 21860, 24860, 25360, 25990, 26390.4, 105064, 106983},
{9.31751, 17.88, 26.58, 39.65, 50.39, 72, 85.1, 130.1, 149, 169, 192.1, 212, 236, 263, 287, 311, 335, 452, 481, 510, 540, 600, 630, 720, 750, 790, 860, 920, 990, 1050, 1120, 1180, 1250, 1320, 1380, 1450, 1510, 1590, 1650, 1948, 2007, 2071, 2139, 2203, 2266, 2373, 2439, 2510, 2576, 2841, 2905, 2977, 3042, 3312, 3388, 3573, 3649, 5976, 6122, 6279, 6445, 6604, 6759, 7068, 7230, 7410, 7570, 8030, 8180, 8330, 8480, 9330, 9500, 9830, 9990, 21210, 21600, 22030, 22420, 25580, 26090, 26730, 27139, 107923, 109886},
{10.7485, 21.4, 29.4, 36.9, 52.9, 64, 88, 102, 154, 173.9, 195, 218, 240, 264, 293, 317, 342, 367, 488, 520, 550, 580, 640, 680, 760, 800, 850, 920, 980, 1050, 1110, 1180, 1250, 1310, 1390, 1460, 1520, 1590, 1660, 1720, 2033, 2094, 2158, 2227, 2293, 2357, 2467, 2535, 2606, 2674, 2944, 3010, 3082, 3149, 3433, 3510, 3699, 3777, 6169, 6318, 6476, 6646, 6807, 6964, 7283, 7450, 7630, 7800, 8260, 8410, 8570, 8710, 9610, 9780, 10120, 10290, 21770, 22160, 22600, 22990, 26310, 26830, 27490, 27903.1, 110842, 112844},
{4.07274, 22.4, 33.5, 39.1, 50, 67, 80, 106, 120, 179, 200, 222.1, 245, 269, 293, 324, 349, 375, 400, 530, 560, 590, 620, 690, 720, 810, 850, 910, 980, 1040, 1110, 1180, 1250, 1320, 1380, 1460, 1530, 1600, 1670, 1740, 1810, 2119, 2182, 2247, 2317, 2384, 2450, 2564, 2631, 2706, 2774, 3049, 3115, 3190, 3257, 3556, 3635, 3828, 3907, 6365, 6516, 6678, 6849, 7013, 7172, 7500, 7670, 7850, 8020, 8500, 8640, 8800, 8950, 9890, 10070, 10420, 10590, 22330, 22730, 23170, 23570, 27060, 27590, 28260, 28683.4, 113817, 115859},
{5.27842, 10.1472, 31, 41, 52.9, 64, 82, 97, 124, 140, 204.9, 227, 250, 274, 299, 324, 356, 382, 409, 435, 570, 600, 630, 660, 740, 770, 860, 900, 970, 1040, 1110, 1180, 1250, 1320, 1390, 1460, 1530, 1610, 1680, 1750, 1820, 1880, 2208, 2271, 2338, 2409, 2477, 2544, 2662, 2731, 2806, 2876, 3155, 3224, 3298, 3368, 3682, 3762, 3959, 4040, 6565, 6718, 6881, 7056, 7222, 7380, 7720, 7890, 8080, 8250, 8730, 8880, 9040, 9200, 10190, 10360, 10720, 10890, 22900, 23300, 23750, 24160, 27830, 28370, 29050, 29479.8, 116849, 118931},
{5.38023, 11.75, 17.431, 44.8, 55, 67, 79, 98.9, 113.9, 143.9, 161.1, 233, 255, 279, 305, 330, 355, 390, 416, 444, 470, 610, 640, 670, 710, 780, 820, 920, 950, 1030, 1100, 1170, 1240, 1310, 1380, 1460, 1530, 1610, 1680, 1750, 1820, 1900, 1970, 2298, 2362, 2430, 2503, 2572, 2639, 2762, 2833, 2908, 2980, 3264, 3334, 3409, 3479, 3811, 3893, 4093, 4175, 6767, 6923, 7088, 7265, 7430, 7600, 7950, 8120, 8310, 8480, 8970, 9120, 9290, 9440, 10480, 10660, 11030, 11200, 23480, 23890, 24340, 24760, 28610, 29160, 29850, 30293.1, 119939, 122063},
{6.3067, 12.1, 18.32, 28.648, 58, 69.1, 82, 95, 118, 133, 165, 181, 262, 285, 310, 336, 362, 389, 424, 451, 480, 508, 650, 680, 720, 750, 830, 870, 970, 1010, 1090, 1160, 1240, 1310, 1380, 1460, 1530, 1600, 1680, 1760, 1830, 1910, 1980, 2060, 2390, 2455, 2524, 2598, 2669, 2737, 2864, 2935, 3013, 3086, 3375, 3445, 3522, 3593, 3943, 4025, 4230, 4313, 6972, 7130, 7299, 7480, 7650, 7810, 8180, 8350, 8550, 8720, 9220, 9370, 9540, 9690, 10790, 10970, 11340, 11510, 24060, 24480, 24940, 25360, 29410, 29970, 30680, 31122.8, 123086, 125253},
{5.89, 11.9, 18.6, 30.9, 44.3, 72, 85.1, 98.9, 111, 137, 153, 187, 203, 292, 316, 342, 369, 395, 423, 460, 488, 518, 546, 690, 720, 760, 790, 880, 920, 1020, 1060, 1150, 1220, 1300, 1370, 1450, 1520, 1600, 1670, 1760, 1830, 1910, 1980, 2060, 2130, 2483, 2550, 2620, 2696, 2766, 2837, 2968, 3040, 3119, 3193, 3488, 3558, 3637, 3709, 4077, 4161, 4370, 4454, 7181, 7341, 7510, 7690, 7870, 8040, 8410, 8590, 8780, 8960, 9460, 9620, 9790, 9950, 11100, 11290, 11660, 11840, 24660, 25080, 25540, 25970, 30230, 30800, 31520, 31971.6, 126297, 128507},
{6.19405, 11.6, 19.8, 36.7, 46, 62, 89, 101, 116, 128.9, 158, 173, 210, 227, 323, 348, 375, 402, 431, 458, 497, 525, 557, 585, 730, 770, 800, 840, 930, 970, 1070, 1110, 1210, 1290, 1370, 1440, 1520, 1590, 1670, 1750, 1830, 1910, 1990, 2070, 2140, 2220, 2578, 2646, 2718, 2794, 2867, 2938, 3073, 3147, 3228, 3301, 3602, 3675, 3753, 3827, 4214, 4299, 4513, 4598, 7393, 7550, 7730, 7910, 8090, 8260, 8650, 8830, 9030, 9210, 9720, 9870, 10040, 10200, 11410, 11600, 11990, 12160, 25260, 25680, 26150, 26590, 31060, 31640, 32400, 32836.5, 129570, 131821},
{6.26554, 11.5, 19.7, 33.8, 48, 65, 92, 107, 121, 136, 151, 179, 196, 233, 252, 355, 382, 408, 438, 466, 495, 535, 565, 596, 626, 770, 810, 850, 880, 980, 1020, 1130, 1170, 1280, 1360, 1430, 1510, 1590, 1670, 1740, 1820, 1910, 1990, 2070, 2140, 2230, 2310, 2675, 2745, 2817, 2894, 2969, 3041, 3181, 3255, 3338, 3413, 3718, 3792, 3872, 3947, 4353, 4441, 4658, 4744, 7610, 7770, 7950, 8130, 8310, 8480, 8890, 9070, 9270, 9450, 9970, 10130, 10300, 10470, 11730, 11930, 12320, 12500, 25870, 26300, 26770, 27210, 31910, 32500, 33300, 33722.2, 132902, 135202},
{6.02576, 11.5, 21.1, 35, 49, 80, 95, 109, 124, 139, 159, 179, 200, 219, 258, 278, 389, 416, 444, 474, 503, 532, 575, 605, 637, 668, 820, 850, 890, 930, 1030, 1070, 1180, 1220, 1340, 1420, 1500, 1580, 1660, 1740, 1820, 1890, 1990, 2070, 2150, 2230, 2310, 2390, 2774, 2844, 2918, 2997, 3072, 3146, 3290, 3366, 3449, 3527, 3836, 3911, 3993, 4068, 4496, 4585, 4807, 4890, 7830, 7990, 8170, 8360, 8540, 8710, 9130, 9310, 9520, 9700, 10230, 10390, 10570, 10730, 12060, 12260, 12660, 12840, 26480, 26920, 27400, 27840, 32800, 33400, 34100, 34625.8, 136299, 138646},
{5.97381, 11.7, 21.7, 36.8, 50, 67.9, 95, 110, 125, 141, 163, 184, 206, 225, 242, 284, 305, 424, 451, 481, 511, 541, 571, 616, 646, 680, 711, 870, 900, 940, 980, 1090, 1130, 1240, 1280, 1410, 1490, 1570, 1650, 1730, 1820, 1900, 1980, 2070, 2160, 2240, 2320, 2410, 2480, 2874, 2946, 3021, 3101, 3178, 3251, 3402, 3479, 3563, 3641, 3956, 4033, 4115, 4191, 4642, 4733, 4960, 5050, 8040, 8210, 8390, 8590, 8770, 8950, 9380, 9560, 9770, 9960, 10490, 10650, 10830, 11000, 12400, 12600, 13000, 13190, 27110, 27550, 28040, 28500, 33700, 34300, 35100, 35549.4, 139770, 142161},
{5.99141, 12.4, 20.1, 37.7, 51, 69.1, 97, 112, 128, 144, 167, 190, 213, 235, 253, 272, 311, 332, 460, 489, 518, 550, 580, 611, 657, 689, 723, 755, 910, 950, 990, 1030, 1140, 1180, 1300, 1340, 1480, 1560, 1650, 1730, 1810, 1890, 1980, 2060, 2160, 2240, 2320, 2410, 2490, 2580, 2976, 3050, 3125, 3207, 3284, 3360, 3515, 3593, 3679, 3758, 4078, 4156, 4239, 4317, 4791, 4880, 5110, 5200, 8270, 8440, 8620, 8820, 9000, 9180, 9630, 9820, 10020, 10220, 10760, 10920, 11100, 11270, 12740, 12950, 13350, 13550, 27740, 28180, 28700, 29100, 34600, 35200, 36000, 36493, 143300, 145743},
{6.19785, 11.9, 21.6, 36, 56, 70.1, 90, 114, 130, 147, 171, 195, 218, 240, 259, 279, 303, 339, 361, 497, 526, 557, 590, 621, 652, 700, 733, 768, 800, 960, 1000, 1040, 1080, 1200, 1240, 1360, 1410, 1550, 1630, 1720, 1800, 1890, 1970, 2050, 2140, 2240, 2320, 2410, 2490, 2580, 2670, 3080, 3154, 3232, 3315, 3393, 3469, 3630, 3709, 3797, 3877, 4202, 4281, 4365, 4445, 4940, 5040, 5270, 5360, 8500, 8670, 8850, 9050, 9240, 9420, 9880, 10070, 10280, 10480, 11020, 11190, 11380, 11550, 13090, 13300, 13720, 13910, 28380, 28800, 29300, 29800, 35500, 36200, 37000, 37457.6, 146905, 149398},
{6.28166, 12, 22.4, 37.7, 51.9, 75, 91, 112.9, 133, 152, 178, 201, 225, 247, 265, 286, 310, 334, 368, 390, 536, 566, 597, 630, 662, 695, 744, 778, 814, 847, 1010, 1050, 1090, 1120, 1250, 1300, 1420, 1470, 1620, 1700, 1790, 1880, 1960, 2050, 2130, 2220, 2320, 2410, 2490, 2580, 2670, 2750, 3186, 3261, 3340, 3424, 3503, 3581, 3747, 3828, 3915, 3998, 4329, 4407, 4494, 4570, 5100, 5190, 5430, 5520, 8730, 8900, 9090, 9290, 9480, 9660, 10140, 10330, 10550, 10740, 11300, 11470, 11650, 11820, 13450, 13660, 14080, 14280, 29000, 29500, 30000, 30500, 36500, 37100, 37900, 38443.5, 150579, 153124},
{6.36758, 12.2, 22.7, 38.8, 54.1, 71, 97, 112.9, 137, 157, 180, 206, 231, 252, 270, 294, 317, 342, 367, 398, 421, 576, 606, 638, 672, 705, 738, 790, 824, 861, 895, 1060, 1100, 1140, 1180, 1310, 1360, 1480, 1530, 1690, 1780, 1870, 1950, 2040, 2130, 2220, 2300, 2410, 2490, 2580, 2680, 2760, 2850, 3294, 3370, 3449, 3535, 3616, 3694, 3866, 3947, 4038, 4120, 4456, 4537, 4620, 4700, 5260, 5350, 5600, 5690, 8960, 9140, 9330, 9530, 9720, 9910, 10400, 10590, 10810, 11010, 11570, 11740, 11930, 12110, 13810, 14030, 14460, 14700, 29700, 30100, 30700, 31100, 37400, 38100, 38900, 39451.4, 154328, 156926},
{6.5, 12.4, 23.2, 39.3, 55, 74, 93, 120, 136, 162, 185, 209, 237, 257, 276, 300, 326, 351, 377, 402, 430, 453, 616, 647, 680, 716, 749, 782, 837, 871, 909, 944, 1110, 1150, 1190, 1230, 1370, 1420, 1550, 1600, 1770, 1850, 1940, 2030, 2120, 2210, 2300, 2390, 2490, 2590, 2680, 2760, 2850, 2950, 3403, 3480, 3561, 3647, 3730, 3810, 3986, 4070, 4160, 4245, 4586, 4670, 4760, 4840, 5420, 5510, 5760, 5860, 9200, 9370, 9570, 9770, 9970, 10160, 10660, 10860, 11080, 11280, 11850, 12020, 12220, 12390, 14180, 14400, 14800, 15000, 30300, 30800, 31300, 31800, 38400, 39100, 40000, 40482.2, 158152, 160804},
{6.58, 12.4, 24.3, 40, 54.1, 76, 96, 115.1, 143.9, 162, 187, 215, 240, 260, 282, 307, 334, 360, 386, 412, 438, 462, 486, 659, 690, 723, 760, 794, 828, 885, 920, 958, 994, 1160, 1210, 1250, 1290, 1430, 1480, 1620, 1660, 1840, 1930, 2020, 2110, 2200, 2290, 2390, 2480, 2580, 2680, 2760, 2860, 2950, 3050, 3513, 3592, 3675, 3762, 3845, 3926, 4109, 4194, 4286, 4371, 4720, 4800, 4890, 4970, 5580, 5680, 5930, 6030, 9430, 9620, 9810, 10020, 10220, 10410, 10930, 11130, 11350, 11560, 12130, 12310, 12500, 12680, 14560, 14800, 15200, 15400, 31000, 31500, 32000, 32500, 39500, 40100, 41000, 41548, 0, 164764},
{6.66, 12.93, 25.8, 41.5, 60, 74, 97, 119, 140, 170, 187, 216, 246, 267, 285, 312, 341, 367, 394, 422, 448, 475, 496, 520, 701, 734, 768, 805, 840, 875, 934, 969, 1010, 1045, 1220, 1260, 1300, 1350, 1500, 1550, 1680, 1730, 1920, 2010, 2110, 2200, 2290, 2380, 2470, 2570, 2680, 2760, 2860, 2950, 3050, 3140, 3627, 3705, 3790, 3878, 3962, 4045, 4234, 4320, 4413, 4500, 4850, 4930, 5030, 5110, 5750, 5850, 6110, 6210, 9680, 9860, 10060, 10270, 10470, 10660, 11200, 11410, 11630, 11840, 12420, 12600, 12800, 12980, 15000, 15200, 15600, 15800, 31700, 32200, 32700, 33200, 40500, 41200, 42100, 42632, 0, 168806},
{4.96, 14.54, 21.8, 43.6, 56, 80, 96, 121, 143, 165, 197, 216, 244, 269, 290, 322, 344, 374, 403, 431, 459, 487, 510, 540, 560, 745, 779, 814, 852, 888, 922, 985, 1020, 1061, 1098, 1280, 1320, 1360, 1410, 1570, 1620, 1760, 1810, 2010, 2100, 2190, 2290, 2380, 2470, 2570, 2670, 2780, 2860, 2960, 3060, 3150, 3250, 3741, 3821, 3906, 3996, 4082, 4165, 4360, 4448, 4540, 4630, 4990, 5070, 5160, 5250, 5920, 6030, 6290, 6390, 9920, 10110, 10310, 10520, 10720, 10920, 11470, 11680, 11910, 12120, 12710, 12890, 13090, 13300, 15300, 15600, 16000, 16200, 32400, 32900, 33400, 33900, 41600, 42300, 43200, 43759, 0, 172930},
{6.02, 14.35, 23.84, 31.87, 64, 77, 102, 119, 146.1, 169, 193, 225, 244, 275, 0, 791, 825, 860, 899, 936, 972, 1036, 1073, 1114, 1151, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3857, 3938, 4025, 4116, 4203, 4287, 4489, 4580, 4670, 4760, 5130, 5210, 5300, 5390, 6100, 6200, 6470, 6570, 10170, 10360, 10560, 10780, 10980, 11180, 11750, 11960, 12200, 12410, 13010, 13190, 13400, 13600, 15800, 16000, 16400, 16700, 33100, 33600, 34100, 34600, 42700, 43400, 44300, 0, 0, 177148},
{6.8, 14, 23.1, 33, 43, 86, 98.9, 126, 145.1, 172, 196, 220.9, 254, 274, 307, 0, 838, 872, 908, 948, 985, 1022, 1089, 1126, 1168, 1207, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3975, 4057, 4145, 4237, 4326, 4411, 4620, 4710, 4810, 4900, 5260, 5350, 5450, 5530, 6280, 6380, 6650, 6760, 10420, 10610, 10820, 11040, 11240, 11440, 12040, 12250, 12480, 12700, 13300, 13500, 13700, 13900, 16200, 16400, 16900, 17100, 33800, 34300, 34800, 35300, 43800, 44500, 45400, 0, 0, 181444},
{7.8, 17.1, 25.8, 35.5, 47.2, 59.3, 109, 122, 152, 170, 200, 224, 251, 285, 306, 339, 0, 885, 921, 958, 998, 1036, 1073, 1143, 1181, 1223, 1263, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4095, 4178, 4267, 4360, 4450, 4540, 4750, 4840, 4940, 5030, 5410, 5490, 5590, 5680, 6460, 6570, 6840, 6950, 10680, 10870, 11080, 11300, 11510, 11710, 12320, 12540, 12780, 12990, 13600, 13800, 14000, 14200, 16600, 16800, 17300, 17500, 34500, 35000, 35600, 36100, 44900, 45700, 46600, 0, 0, 185839},
{7.7, 17.5, 26.7, 37.3, 49, 62.1, 74.9, 134, 148, 178, 198, 228, 255, 281, 318, 337, 374, 0, 934, 969, 1008, 1049, 1088, 1126, 1197, 1237, 1280, 1320, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4216, 4301, 4390, 4486, 4580, 4660, 4890, 4980, 5080, 5170, 5550, 5640, 5740, 5830, 6650, 6760, 7040, 7140, 10930, 11130, 11340, 11560, 11780, 11980, 12610, 12830, 13070, 13300, 13900, 14100, 14300, 14500, 17000, 17300, 17700, 18000, 35200, 35700, 36300, 36800, 46100, 46900, 47800, 0, 0, 190331},
{7.6, 18.2, 29.3, 37.7, 51.2, 64, 78.1, 91.7, 159.9, 173.9, 206.1, 227, 258, 285, 314, 351, 371, 409, 0, 984, 1020, 1060, 1101, 1140, 1180, 1253, 1294, 1338, 1379, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4339, 4425, 4516, 4610, 4700, 4790, 5020, 5110, 5220, 5310, 5700, 5780, 5880, 5980, 6840, 6950, 7230, 7340, 11200, 11390, 11610, 11830, 12040, 12250, 12910, 13130, 13400, 13600, 14200, 14400, 14600, 14800, 17500, 17700, 18200, 18400, 35900, 36400, 37000, 37500, 47300, 48100, 49000, 0, 0, 194917},
{50, 0, 0, 94, 109, 187, 202, 235.9, 257, 289, 318, 346, 386, 406, 445, 0, 1035, 1072, 1112, 1154, 1195, 1234, 1311, 1352, 1397, 1439, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4464, 4551, 4640, 4740, 4830, 4920, 5160, 5250, 5360, 5450, 5840, 5930, 6030, 6130, 7030, 7150, 7430, 7550, 11460, 11660, 11870, 12100, 12320, 12530, 13200, 13400, 13700, 13900, 14500, 14700, 14900, 15100, 17900, 18200, 18700, 18900, 36700, 37200, 37800, 38300, 48500, 49400, 50300, 0, 0, 199606},
{65, 0, 0, 112.9, 128, 216, 231, 266, 288, 322, 352, 380, 422, 442, 483, 0, 1087, 1125, 1165, 1208, 1250, 1290, 1369, 1412, 1457, 1500, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4590, 4680, 4770, 4870, 4960, 5060, 5300, 5400, 5500, 5600, 5990, 6080, 6190, 6280, 7230, 7350, 7640, 7750, 11730, 11930, 12140, 12380, 12600, 12810, 13500, 13700, 14000, 14200, 14800, 15000, 15300, 15500, 18400, 18600, 19100, 19400, 37400, 37900, 38500, 39100, 49800, 50700, 51600, 0, 0, 204400}
};

fp invI05_log(fp log_z); // Обратная функция к функции Ферми-Дирака I_{1/2}(log(z))
fp I15_invI05(fp z); // Функция Ферми-Дирака I_{3/2}(I_{1/2}^{-1}(z))

fp mu_log(fp T, fp log_V, fp log_xe); // mu(log(xe))

fp I05_mu_T(fp T, fp V, fp xe); // Величина I_{1/2}(mu/T)
fp I15_mu_T(fp T, fp V, fp xe); // Величина I_{3/2}(mu/T)

} // namespace Constants

#endif // Constants_H
