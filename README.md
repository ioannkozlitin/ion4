# ion4
<p>
Код для расчета состава и термодинамики равновесной плазмы по модели Саха в двухтемпературном случае.
</p>
Формат вызова: <br>
./ion4 файл_сценария выходной_файл_в_формате_матлаб
<br><br>
Пример сценария:<br>
<pre>
#Компоненты смеси (воздух)
Z: 7 8 18
#Концентрации компонент смеси
x: 0.7811 0.2095 0.0094
#Диапазон по плоности lg(г/см^3): min max step
lgRho: -6 0 1
#Диапазон по ионной температуре lg(эВ): min max step
lgT: 0 2 1
#Диапазон отличия электронной температуры от ионной lg(эВ): minDiff maxDiff
lgtDiffT: -1 1
#Коэффициент объема ионных остовов, отвечающий за ионизацию сжатием
#Для горячей неплотной плазмы можно взять 0 (быстрее счет), но штатное значение 0.6
rCoeff: 0
</pre>
Строки, начинающиеся с '#' - комментарии.<br><br>

Выходной файл имеет матлабовский формат.<br>
Пример выходного файла (результат расчета по сценарию, данному выше):<br>
<pre>
Z=[7 8 18 ];
x=[7.811000e-01 2.095000e-01 9.400000e-03 ];
lgV = [2.215817e+00 3.215817e+00 4.215817e+00 5.215817e+00 6.215817e+00 7.215817e+00 8.215817e+00 ];
lgRho = [0.000000e+00 -1.000000e+00 -2.000000e+00 -3.000000e+00 -4.000000e+00 -5.000000e+00 -6.000000e+00 ];
xe = [
3.000000e+00 2.000000e+00 4.590227e+00 4.965391e+00 5.487965e+00 5.888209e+00 6.185414e+00 6.411676e+00 6.588385e+00 
2.000000e+00 2.000000e+00 5.234056e+00 5.765987e+00 6.621536e+00 7.145200e+00 7.275491e+00 7.292184e+00 7.293908e+00 
1.000000e+00 2.000000e+00 5.237720e+00 5.715334e+00 7.047153e+00 7.294091e+00 7.294100e+00 7.294100e+00 7.294100e+00 
2.000000e+00 1.000000e+00 2.123432e+00 2.060824e+00 2.511006e+00 3.006991e+00 3.475794e+00 3.880790e+00 4.208517e+00 
1.000000e+00 1.000000e+00 9.612152e-01 1.645863e+00 2.559408e+00 3.339312e+00 4.175685e+00 4.796221e+00 5.090338e+00 
0.000000e+00 1.000000e+00 8.787052e-01 1.773649e+00 2.777722e+00 3.027423e+00 4.027850e+00 4.936688e+00 5.081493e+00 
1.000000e+00 0.000000e+00 4.642311e-01 2.496761e-01 2.926984e-01 3.478783e-01 4.116328e-01 4.838896e-01 5.646976e-01 
0.000000e+00 0.000000e+00 9.639073e-04 9.418266e-04 3.047413e-03 9.620909e-03 3.004842e-02 9.136946e-02 2.570572e-01 
-1.000000e+00 0.000000e+00 6.253150e-08 4.764632e-06 4.164164e-05 3.402499e-04 2.758309e-03 2.217559e-02 1.595887e-01 
];
P = [
3.000000e+00 2.000000e+00 3.524576e-01 3.366686e-02 3.463925e-03 3.552338e-04 3.618729e-05 3.669313e-06 3.708821e-07 
2.000000e+00 2.000000e+00 1.413205e-01 1.514180e-02 1.704149e-03 1.821113e-04 1.850231e-05 1.853962e-06 1.854347e-07 
1.000000e+00 2.000000e+00 1.209379e-01 1.301503e-02 1.598092e-03 1.653181e-04 1.653170e-05 1.653169e-06 1.653169e-07 
2.000000e+00 1.000000e+00 4.080459e-02 2.799929e-03 2.805505e-04 2.908753e-05 3.012962e-06 3.103463e-07 3.176732e-08 
1.000000e+00 1.000000e+00 1.004302e-02 6.161182e-04 7.978085e-05 9.703642e-06 1.157196e-06 1.295918e-07 1.361675e-08 
0.000000e+00 1.000000e+00 5.163653e-03 4.347185e-04 6.449493e-05 6.993807e-06 9.229197e-07 1.126102e-07 1.158475e-08 
1.000000e+00 0.000000e+00 9.227319e-03 2.497004e-04 2.319993e-05 2.315404e-06 2.328007e-07 2.344000e-08 2.362051e-09 
0.000000e+00 0.000000e+00 2.633074e-03 2.463305e-05 2.263300e-06 2.259365e-07 2.303185e-08 2.440099e-09 2.810524e-10 
-1.000000e+00 0.000000e+00 2.643184e-04 2.461202e-06 2.257389e-07 2.245457e-08 2.297674e-09 2.731620e-10 5.803875e-11 
];
E = [
3.000000e+00 2.000000e+00 9.633418e+01 1.024389e+02 1.108632e+02 1.183763e+02 1.247879e+02 1.302313e+02 1.348409e+02 
2.000000e+00 2.000000e+00 4.711427e+01 5.973169e+01 8.430702e+01 1.010004e+02 1.054820e+02 1.060678e+02 1.061284e+02 
1.000000e+00 2.000000e+00 4.080355e+01 5.289199e+01 9.198638e+01 1.011740e+02 1.011741e+02 1.011740e+02 1.011740e+02 
2.000000e+00 1.000000e+00 9.356834e+00 9.247907e+00 1.036593e+01 1.170594e+01 1.309296e+01 1.440344e+01 1.555989e+01 
1.000000e+00 1.000000e+00 1.739530e+00 2.817033e+00 4.694733e+00 6.886787e+00 9.880101e+00 1.243386e+01 1.383666e+01 
0.000000e+00 1.000000e+00 1.039894e+00 2.411217e+00 4.600119e+00 5.251932e+00 8.662150e+00 1.251281e+01 1.325818e+01 
1.000000e+00 0.000000e+00 8.848345e-01 7.048072e-01 7.326387e-01 7.695129e-01 8.133562e-01 8.646829e-01 9.242456e-01 
0.000000e+00 0.000000e+00 5.568882e-02 5.566645e-02 5.687896e-02 6.066511e-02 7.243595e-02 1.078162e-01 2.037287e-01 
-1.000000e+00 0.000000e+00 5.512431e-03 5.515044e-03 5.535532e-03 5.701437e-03 7.044894e-03 1.783303e-02 9.418542e-02 
];
S = [
3.000000e+00 2.000000e+00 5.559473e+01 7.148519e+01 8.998439e+01 1.093783e+02 1.291683e+02 1.491800e+02 1.693244e+02 
2.000000e+00 2.000000e+00 5.441597e+01 7.280609e+01 9.605865e+01 1.188550e+02 1.390236e+02 1.582638e+02 1.773765e+02 
1.000000e+00 2.000000e+00 5.041310e+01 6.872746e+01 9.548753e+01 1.166307e+02 1.357285e+02 1.548264e+02 1.739243e+02 
2.000000e+00 1.000000e+00 3.089423e+01 3.884321e+01 4.910021e+01 6.160098e+01 7.593194e+01 9.147133e+01 1.076223e+02 
1.000000e+00 1.000000e+00 2.323923e+01 3.281020e+01 4.517465e+01 6.024844e+01 7.934895e+01 9.900527e+01 1.165480e+02 
0.000000e+00 1.000000e+00 1.893104e+01 2.906848e+01 4.230191e+01 5.320090e+01 7.343805e+01 9.695570e+01 1.128251e+02 
1.000000e+00 0.000000e+00 1.973168e+01 2.374417e+01 2.700775e+01 3.047872e+01 3.426987e+01 3.844209e+01 4.305156e+01 
0.000000e+00 0.000000e+00 1.397278e+01 1.864866e+01 2.107526e+01 2.350235e+01 2.616742e+01 2.956026e+01 3.484537e+01 
-1.000000e+00 0.000000e+00 1.050158e+01 1.517826e+01 1.756856e+01 1.988618e+01 2.223976e+01 2.490802e+01 2.949544e+01 
];
Si = [
3.000000e+00 2.000000e+00 2.853876e+01 3.084268e+01 3.294752e+01 3.503649e+01 3.713620e+01 3.924836e+01 4.137177e+01 
2.000000e+00 2.000000e+00 2.390031e+01 2.643987e+01 2.848179e+01 3.002544e+01 3.195328e+01 3.417374e+01 3.646394e+01 
1.000000e+00 2.000000e+00 1.987235e+01 2.271798e+01 2.400594e+01 2.610052e+01 2.840302e+01 3.070560e+01 3.300819e+01 
2.000000e+00 1.000000e+00 2.462323e+01 2.733523e+01 2.971761e+01 3.200155e+01 3.421764e+01 3.638836e+01 3.853835e+01 
1.000000e+00 1.000000e+00 2.004252e+01 2.325140e+01 2.546481e+01 2.772746e+01 3.000100e+01 3.194470e+01 3.395715e+01 
0.000000e+00 1.000000e+00 1.593658e+01 1.888881e+01 2.113694e+01 2.342057e+01 2.569201e+01 2.807362e+01 3.036891e+01 
1.000000e+00 0.000000e+00 1.937235e+01 2.269444e+01 2.513140e+01 2.750572e+01 2.987326e+01 3.223775e+01 3.459806e+01 
0.000000e+00 0.000000e+00 1.396789e+01 1.863947e+01 2.104182e+01 2.338561e+01 2.576783e+01 2.823643e+01 3.079493e+01 
-1.000000e+00 0.000000e+00 1.050158e+01 1.517819e+01 1.756792e+01 1.988092e+01 2.219650e+01 2.455533e+01 2.690474e+01 
];
Se = [
3.000000e+00 2.000000e+00 2.705597e+01 4.064250e+01 5.703687e+01 7.434177e+01 9.203208e+01 1.099317e+02 1.279526e+02 
2.000000e+00 2.000000e+00 3.051566e+01 4.636621e+01 6.757686e+01 8.882955e+01 1.070703e+02 1.240901e+02 1.409126e+02 
1.000000e+00 2.000000e+00 3.054074e+01 4.600948e+01 7.148158e+01 9.053014e+01 1.073255e+02 1.241208e+02 1.409161e+02 
2.000000e+00 1.000000e+00 6.271002e+00 1.150798e+01 1.938260e+01 2.959942e+01 4.171430e+01 5.508296e+01 6.908392e+01 
1.000000e+00 1.000000e+00 3.196713e+00 9.558803e+00 1.970984e+01 3.252098e+01 4.934795e+01 6.706058e+01 8.259089e+01 
0.000000e+00 1.000000e+00 2.994465e+00 1.017967e+01 2.116497e+01 2.978033e+01 4.774604e+01 6.888208e+01 8.245621e+01 
1.000000e+00 0.000000e+00 3.593376e-01 1.049725e+00 1.876350e+00 2.972994e+00 4.396609e+00 6.204338e+00 8.453501e+00 
0.000000e+00 0.000000e+00 4.886175e-03 9.189738e-03 3.343759e-02 1.167368e-01 3.995893e-01 1.323829e+00 4.050438e+00 
-1.000000e+00 0.000000e+00 9.192512e-07 7.167844e-05 6.356751e-04 5.265590e-03 4.326787e-02 3.526953e-01 2.590707e+00 
];
</pre>
<p>
Здесь xe - степень ионизации, P - давление в а.е., E - энергия в а.е., Si - ионная компонента энтропии, Se - электронная компонента, S - их сумма.
Первые два столбца любой строки в этих таблицах - десятичные логарифмы электронной и ионной температур в эВ соответственно. Все остальные элементы строки - величины соответствующей 
термодинамической функции при разных плотностях, заданных в строке lgRho. Строка lgV - десятичные логарифмы объемов атомных ячеек для этих плотностей в а.е.
</p>
