<!-- research-repo-tools/run-report/v2 -->
# la-stack complete benchmark measurements

Run: `634bc9527c2caf0a431a9bdfac082c190e1094f0b29eb6e3f28c4b96961cdec5`

Evidence: [run payload](performance-runs/runs/634bc9527c2caf0a431a9bdfac082c190e1094f0b29eb6e3f28c4b96961cdec5/run.json) · [evidence envelope](performance-runs/runs/634bc9527c2caf0a431a9bdfac082c190e1094f0b29eb6e3f28c4b96961cdec5/evidence.json).

Intervals are marginal timing intervals. No significance or ratio interval is inferred.

- la-stack baseline: measured in **baseline** phase (v0.4.6), revision `0215e65763f8eecb4bcad8bf29de4aca67531aa1`.
- la-stack current: measured in **current** phase (v0.5.0), revision `51868ea4725c1a4790bba4c59341d2113b18b2bb`.
- nalgebra: measured in **baseline** phase (v0.4.6), revision `0215e65763f8eecb4bcad8bf29de4aca67531aa1`.
- faer: measured in **baseline** phase (v0.4.6), revision `0215e65763f8eecb4bcad8bf29de4aca67531aa1`.

## Mean (ns)

| Case | la-stack baseline | la-stack current | nalgebra | faer |
| --- | --- | --- | --- | --- |
| d16/la&#95;stack&#95;det | 381.877 ns [380.637, 383.185] (0.95 confidence) | 402.901 ns [400.348, 405.538] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;det&#95;from&#95;ldlt | 2.88903 ns [2.88509, 2.89272] (0.95 confidence) | 2.96685 ns [2.96414, 2.96959] (0.95 confidence) | 1.69417 ns [1.68903, 1.69915] (0.95 confidence) | 4.38111 ns [4.35599, 4.40639] (0.95 confidence) |
| d16/la&#95;stack&#95;det&#95;from&#95;lu | 3.27537 ns [3.27145, 3.27908] (0.95 confidence) | 3.36911 ns [3.35909, 3.38387] (0.95 confidence) | 1.77288 ns [1.76874, 1.77689] (0.95 confidence) | 4.84201 ns [4.80326, 4.88737] (0.95 confidence) |
| d16/la&#95;stack&#95;det&#95;via&#95;lu | 361.076 ns [360.149, 362.072] (0.95 confidence) | 378.232 ns [376.384, 380.127] (0.95 confidence) | 449.519 ns [447.78, 451.056] (0.95 confidence) | 670.718 ns [668.126, 673.372] (0.95 confidence) |
| d16/la&#95;stack&#95;dot | 2.33211 ns [2.3303, 2.33385] (0.95 confidence) | 2.38358 ns [2.37016, 2.40181] (0.95 confidence) | 1.84951 ns [1.84743, 1.85161] (0.95 confidence) | 3.33075 ns [3.32201, 3.34079] (0.95 confidence) |
| d16/la&#95;stack&#95;inf&#95;norm | 32.3675 ns [32.3089, 32.4231] (0.95 confidence) | 33.5932 ns [33.5291, 33.6597] (0.95 confidence) | 32.4876 ns [32.4397, 32.5321] (0.95 confidence) | 32.4244 ns [32.3467, 32.5071] (0.95 confidence) |
| d16/la&#95;stack&#95;ldlt | 339.413 ns [339.103, 339.73] (0.95 confidence) | 346.101 ns [345.399, 346.847] (0.95 confidence) | 421.755 ns [421.085, 422.399] (0.95 confidence) | 418.304 ns [417.574, 419.051] (0.95 confidence) |
| d16/la&#95;stack&#95;ldlt&#95;solve | 382.168 ns [381.957, 382.379] (0.95 confidence) | 391.674 ns [390.916, 392.472] (0.95 confidence) | 627.153 ns [626.113, 628.155] (0.95 confidence) | 612.374 ns [611.444, 613.316] (0.95 confidence) |
| d16/la&#95;stack&#95;lu | 351.056 ns [349.975, 352.238] (0.95 confidence) | 368.085 ns [366.506, 369.652] (0.95 confidence) | 466.726 ns [465.278, 468.463] (0.95 confidence) | 642.229 ns [638.899, 645.695] (0.95 confidence) |
| d16/la&#95;stack&#95;lu&#95;solve | 657.445 ns [645.98, 669.183] (0.95 confidence) | 705.624 ns [692.434, 718.691] (0.95 confidence) | 591.452 ns [589.116, 594.421] (0.95 confidence) | 937.4 ns [933.866, 940.955] (0.95 confidence) |
| d16/la&#95;stack&#95;norm2&#95;sq | 2.11274 ns [2.10163, 2.12181] (0.95 confidence) | 2.14769 ns [2.14117, 2.15295] (0.95 confidence) | 1.73252 ns [1.73063, 1.73425] (0.95 confidence) | 4.09505 ns [4.09295, 4.09779] (0.95 confidence) |
| d16/la&#95;stack&#95;solve&#95;from&#95;ldlt | 27.7874 ns [27.7543, 27.8255] (0.95 confidence) | 29.0296 ns [28.9327, 29.1361] (0.95 confidence) | 119.383 ns [119.111, 119.643] (0.95 confidence) | 182.77 ns [182.572, 182.977] (0.95 confidence) |
| d16/la&#95;stack&#95;solve&#95;from&#95;lu | 191.094 ns [190.45, 191.748] (0.95 confidence) | 209.039 ns [208.68, 209.397] (0.95 confidence) | 93.5301 ns [93.2916, 93.7761] (0.95 confidence) | 257.724 ns [257.076, 258.442] (0.95 confidence) |
| d2/la&#95;stack&#95;det | 0.572913 ns [0.572679, 0.573159] (0.95 confidence) | 0.574026 ns [0.573741, 0.574325] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;det&#95;from&#95;ldlt | 0.500962 ns [0.498793, 0.503394] (0.95 confidence) | 0.499362 ns [0.497568, 0.501539] (0.95 confidence) | 0.438693 ns [0.437689, 0.440203] (0.95 confidence) | 0.561681 ns [0.558912, 0.564336] (0.95 confidence) |
| d2/la&#95;stack&#95;det&#95;from&#95;lu | 0.520734 ns [0.520625, 0.520855] (0.95 confidence) | 0.522398 ns [0.521702, 0.523272] (0.95 confidence) | 0.454669 ns [0.452677, 0.457018] (0.95 confidence) | 0.744466 ns [0.744105, 0.744866] (0.95 confidence) |
| d2/la&#95;stack&#95;det&#95;via&#95;lu | 1.77012 ns [1.76894, 1.7713] (0.95 confidence) | 1.74177 ns [1.73791, 1.74545] (0.95 confidence) | 0.823722 ns [0.822379, 0.825152] (0.95 confidence) | 111.409 ns [110.955, 111.884] (0.95 confidence) |
| d2/la&#95;stack&#95;dot | 0.613022 ns [0.608747, 0.618254] (0.95 confidence) | 0.609238 ns [0.603744, 0.615701] (0.95 confidence) | 0.598959 ns [0.596814, 0.60111] (0.95 confidence) | 2.49854 ns [2.48847, 2.51019] (0.95 confidence) |
| d2/la&#95;stack&#95;inf&#95;norm | 0.703036 ns [0.702728, 0.70336] (0.95 confidence) | 0.703414 ns [0.70302, 0.703828] (0.95 confidence) | 0.559914 ns [0.559665, 0.560193] (0.95 confidence) | 0.702043 ns [0.69971, 0.704562] (0.95 confidence) |
| d2/la&#95;stack&#95;ldlt | 6.89019 ns [6.86966, 6.91184] (0.95 confidence) | 6.57264 ns [6.55792, 6.58782] (0.95 confidence) | 1.72147 ns [1.71961, 1.72315] (0.95 confidence) | 87.8513 ns [87.5879, 88.1188] (0.95 confidence) |
| d2/la&#95;stack&#95;ldlt&#95;solve | 9.63446 ns [9.59263, 9.68741] (0.95 confidence) | 9.49552 ns [9.47364, 9.52069] (0.95 confidence) | 2.61452 ns [2.60732, 2.62125] (0.95 confidence) | 126.156 ns [125.859, 126.469] (0.95 confidence) |
| d2/la&#95;stack&#95;lu | 1.55973 ns [1.55555, 1.56384] (0.95 confidence) | 1.56503 ns [1.56276, 1.56702] (0.95 confidence) | 1.73762 ns [1.72118, 1.75451] (0.95 confidence) | 110.571 ns [109.243, 111.965] (0.95 confidence) |
| d2/la&#95;stack&#95;lu&#95;solve | 2.03852 ns [2.02892, 2.04945] (0.95 confidence) | 2.1184 ns [2.10452, 2.13304] (0.95 confidence) | 4.56865 ns [4.555, 4.58322] (0.95 confidence) | 155.728 ns [155.154, 156.321] (0.95 confidence) |
| d2/la&#95;stack&#95;norm2&#95;sq | 0.395341 ns [0.395056, 0.395603] (0.95 confidence) | 0.394679 ns [0.393798, 0.395596] (0.95 confidence) | 0.388882 ns [0.388576, 0.389166] (0.95 confidence) | 4.12739 ns [4.12475, 4.13006] (0.95 confidence) |
| d2/la&#95;stack&#95;solve&#95;from&#95;ldlt | 1.2274 ns [1.22166, 1.2324] (0.95 confidence) | 1.22675 ns [1.22129, 1.23174] (0.95 confidence) | 1.24962 ns [1.24721, 1.25192] (0.95 confidence) | 39.4451 ns [39.3679, 39.523] (0.95 confidence) |
| d2/la&#95;stack&#95;solve&#95;from&#95;lu | 1.24875 ns [1.24645, 1.25126] (0.95 confidence) | 1.26313 ns [1.26011, 1.26619] (0.95 confidence) | 2.90591 ns [2.8988, 2.91369] (0.95 confidence) | 49.2788 ns [49.1287, 49.4294] (0.95 confidence) |
| d3/la&#95;stack&#95;det | 1.35771 ns [1.35126, 1.36442] (0.95 confidence) | 1.35469 ns [1.34796, 1.36169] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;det&#95;from&#95;ldlt | 0.614083 ns [0.612319, 0.616346] (0.95 confidence) | 0.595951 ns [0.595668, 0.596257] (0.95 confidence) | 0.498063 ns [0.497033, 0.499209] (0.95 confidence) | 0.764778 ns [0.76241, 0.767909] (0.95 confidence) |
| d3/la&#95;stack&#95;det&#95;from&#95;lu | 0.711127 ns [0.710008, 0.712378] (0.95 confidence) | 0.694031 ns [0.693802, 0.694265] (0.95 confidence) | 0.484274 ns [0.483243, 0.485335] (0.95 confidence) | 1.01339 ns [1.01104, 1.01612] (0.95 confidence) |
| d3/la&#95;stack&#95;det&#95;via&#95;lu | 8.82548 ns [8.7993, 8.8528] (0.95 confidence) | 8.71392 ns [8.69608, 8.73389] (0.95 confidence) | 12.7489 ns [12.7001, 12.7992] (0.95 confidence) | 138.965 ns [138.652, 139.286] (0.95 confidence) |
| d3/la&#95;stack&#95;dot | 0.685432 ns [0.683514, 0.687459] (0.95 confidence) | 0.673633 ns [0.672134, 0.675244] (0.95 confidence) | 0.708431 ns [0.706012, 0.711061] (0.95 confidence) | 3.0376 ns [3.01401, 3.06103] (0.95 confidence) |
| d3/la&#95;stack&#95;inf&#95;norm | 1.388 ns [1.38622, 1.38986] (0.95 confidence) | 1.37242 ns [1.36981, 1.37579] (0.95 confidence) | 1.19893 ns [1.19719, 1.20096] (0.95 confidence) | 1.34545 ns [1.34132, 1.34945] (0.95 confidence) |
| d3/la&#95;stack&#95;ldlt | 13.4432 ns [13.3942, 13.4942] (0.95 confidence) | 13.0217 ns [12.9829, 13.0597] (0.95 confidence) | 3.95291 ns [3.94565, 3.96002] (0.95 confidence) | 99.436 ns [99.2032, 99.6837] (0.95 confidence) |
| d3/la&#95;stack&#95;ldlt&#95;solve | 9.86949 ns [9.83823, 9.90216] (0.95 confidence) | 9.81578 ns [9.77303, 9.86134] (0.95 confidence) | 5.73265 ns [5.71507, 5.74869] (0.95 confidence) | 147.806 ns [146.995, 149.019] (0.95 confidence) |
| d3/la&#95;stack&#95;lu | 8.69776 ns [8.64958, 8.75174] (0.95 confidence) | 8.37318 ns [8.35947, 8.38655] (0.95 confidence) | 12.8778 ns [12.7868, 12.9756] (0.95 confidence) | 139.028 ns [138.531, 139.552] (0.95 confidence) |
| d3/la&#95;stack&#95;lu&#95;solve | 9.9778 ns [9.96433, 9.9906] (0.95 confidence) | 9.91611 ns [9.89981, 9.93134] (0.95 confidence) | 33.086 ns [32.9543, 33.2231] (0.95 confidence) | 191.944 ns [191.471, 192.41] (0.95 confidence) |
| d3/la&#95;stack&#95;norm2&#95;sq | 0.509357 ns [0.507526, 0.511601] (0.95 confidence) | 0.4988 ns [0.497818, 0.500049] (0.95 confidence) | 0.45501 ns [0.452661, 0.458448] (0.95 confidence) | 4.26834 ns [4.25797, 4.27851] (0.95 confidence) |
| d3/la&#95;stack&#95;solve&#95;from&#95;ldlt | 1.84441 ns [1.82821, 1.86838] (0.95 confidence) | 1.78465 ns [1.78277, 1.78675] (0.95 confidence) | 3.03299 ns [3.01328, 3.06015] (0.95 confidence) | 41.4997 ns [41.3809, 41.6231] (0.95 confidence) |
| d3/la&#95;stack&#95;solve&#95;from&#95;lu | 2.12983 ns [2.12464, 2.13515] (0.95 confidence) | 2.06385 ns [2.06024, 2.06793] (0.95 confidence) | 4.70164 ns [4.68838, 4.72005] (0.95 confidence) | 54.0274 ns [53.7538, 54.346] (0.95 confidence) |
| d32/la&#95;stack&#95;det | 2178.59 ns [2173.51, 2183.86] (0.95 confidence) | 2252.39 ns [2247.91, 2256.88] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;det&#95;from&#95;ldlt | 7.78595 ns [7.77144, 7.79968] (0.95 confidence) | 8.16646 ns [8.15442, 8.17878] (0.95 confidence) | 3.07986 ns [3.07157, 3.08763] (0.95 confidence) | 8.61234 ns [8.58365, 8.64235] (0.95 confidence) |
| d32/la&#95;stack&#95;det&#95;from&#95;lu | 9.12572 ns [9.10623, 9.14645] (0.95 confidence) | 9.57734 ns [9.54482, 9.61145] (0.95 confidence) | 3.06712 ns [3.06069, 3.07261] (0.95 confidence) | 8.77433 ns [8.74152, 8.80883] (0.95 confidence) |
| d32/la&#95;stack&#95;det&#95;via&#95;lu | 2043.73 ns [2040.32, 2047.24] (0.95 confidence) | 2106.39 ns [2100.31, 2113.01] (0.95 confidence) | 2268.7 ns [2263.43, 2274.04] (0.95 confidence) | 2389.34 ns [2384.5, 2394.1] (0.95 confidence) |
| d32/la&#95;stack&#95;dot | 4.0179 ns [4.01365, 4.02204] (0.95 confidence) | 4.09956 ns [4.09462, 4.1045] (0.95 confidence) | 4.00272 ns [3.99739, 4.00814] (0.95 confidence) | 4.89614 ns [4.88595, 4.90669] (0.95 confidence) |
| d32/la&#95;stack&#95;inf&#95;norm | 127.564 ns [127.21, 127.884] (0.95 confidence) | 133.649 ns [133.334, 134.06] (0.95 confidence) | 160.034 ns [158.875, 161.173] (0.95 confidence) | 162.59 ns [162.29, 162.884] (0.95 confidence) |
| d32/la&#95;stack&#95;ldlt | 2392.54 ns [2386.27, 2402.48] (0.95 confidence) | 2430.23 ns [2426.93, 2433.6] (0.95 confidence) | 2117.9 ns [2097.96, 2145.6] (0.95 confidence) | 1400.06 ns [1395.44, 1404.88] (0.95 confidence) |
| d32/la&#95;stack&#95;ldlt&#95;solve | 2710.85 ns [2707.16, 2714.57] (0.95 confidence) | 2778.27 ns [2771.53, 2785.45] (0.95 confidence) | 2765.06 ns [2748.24, 2781.59] (0.95 confidence) | 1864.35 ns [1861.11, 1867.91] (0.95 confidence) |
| d32/la&#95;stack&#95;lu | 2136.08 ns [2131.23, 2140.81] (0.95 confidence) | 2218.43 ns [2211.08, 2225.47] (0.95 confidence) | 2123.07 ns [2119.46, 2126.64] (0.95 confidence) | 2367.71 ns [2362.3, 2373.17] (0.95 confidence) |
| d32/la&#95;stack&#95;lu&#95;solve | 2764.68 ns [2760.11, 2769.3] (0.95 confidence) | 2835.01 ns [2827.86, 2842.8] (0.95 confidence) | 2693.92 ns [2678.98, 2707.72] (0.95 confidence) | 3030.93 ns [3023.73, 3038.55] (0.95 confidence) |
| d32/la&#95;stack&#95;norm2&#95;sq | 3.96169 ns [3.95641, 3.96658] (0.95 confidence) | 4.11365 ns [4.10097, 4.12831] (0.95 confidence) | 3.61741 ns [3.60959, 3.62457] (0.95 confidence) | 4.2178 ns [4.21225, 4.22402] (0.95 confidence) |
| d32/la&#95;stack&#95;solve&#95;from&#95;ldlt | 311.146 ns [310.312, 312.038] (0.95 confidence) | 320.497 ns [319.642, 321.43] (0.95 confidence) | 547.331 ns [546.18, 548.448] (0.95 confidence) | 469.58 ns [469.101, 470.078] (0.95 confidence) |
| d32/la&#95;stack&#95;solve&#95;from&#95;lu | 670.442 ns [667.718, 673.547] (0.95 confidence) | 736.257 ns [734.399, 738.2] (0.95 confidence) | 328.251 ns [327.657, 328.798] (0.95 confidence) | 639.747 ns [638.454, 641.443] (0.95 confidence) |
| d4/la&#95;stack&#95;det | 2.60575 ns [2.59285, 2.62314] (0.95 confidence) | 2.54743 ns [2.53963, 2.55628] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;det&#95;from&#95;ldlt | 0.792381 ns [0.790149, 0.794662] (0.95 confidence) | 0.771954 ns [0.771298, 0.772628] (0.95 confidence) | 0.504207 ns [0.501823, 0.506781] (0.95 confidence) | 1.00565 ns [1.00376, 1.00786] (0.95 confidence) |
| d4/la&#95;stack&#95;det&#95;from&#95;lu | 0.888149 ns [0.886686, 0.889733] (0.95 confidence) | 0.872037 ns [0.87126, 0.872845] (0.95 confidence) | 0.594308 ns [0.592273, 0.596502] (0.95 confidence) | 1.25897 ns [1.25659, 1.26189] (0.95 confidence) |
| d4/la&#95;stack&#95;det&#95;via&#95;lu | 14.5843 ns [14.5431, 14.6307] (0.95 confidence) | 14.3903 ns [14.3579, 14.4264] (0.95 confidence) | 31.3107 ns [31.2583, 31.366] (0.95 confidence) | 165.205 ns [164.673, 165.745] (0.95 confidence) |
| d4/la&#95;stack&#95;dot | 0.69749 ns [0.694191, 0.701078] (0.95 confidence) | 0.687464 ns [0.685363, 0.689699] (0.95 confidence) | 0.670387 ns [0.668655, 0.672222] (0.95 confidence) | 2.69579 ns [2.67877, 2.71198] (0.95 confidence) |
| d4/la&#95;stack&#95;inf&#95;norm | 2.39947 ns [2.39632, 2.40289] (0.95 confidence) | 2.37035 ns [2.36783, 2.37291] (0.95 confidence) | 2.08687 ns [2.0825, 2.09242] (0.95 confidence) | 2.10548 ns [2.09838, 2.11288] (0.95 confidence) |
| d4/la&#95;stack&#95;ldlt | 18.2527 ns [18.2261, 18.2822] (0.95 confidence) | 17.8838 ns [17.8541, 17.9128] (0.95 confidence) | 8.89904 ns [8.41541, 9.46346] (0.95 confidence) | 124.877 ns [124.327, 125.456] (0.95 confidence) |
| d4/la&#95;stack&#95;ldlt&#95;solve | 18.1806 ns [18.1449, 18.2196] (0.95 confidence) | 18.0042 ns [17.9668, 18.0427] (0.95 confidence) | 11.5421 ns [11.5154, 11.5682] (0.95 confidence) | 166.424 ns [166.006, 166.846] (0.95 confidence) |
| d4/la&#95;stack&#95;lu | 13.9015 ns [13.8778, 13.9256] (0.95 confidence) | 13.781 ns [13.7598, 13.8022] (0.95 confidence) | 30.1319 ns [30.0809, 30.1859] (0.95 confidence) | 158.856 ns [158.369, 159.35] (0.95 confidence) |
| d4/la&#95;stack&#95;lu&#95;solve | 22.1298 ns [22.0973, 22.1616] (0.95 confidence) | 21.7429 ns [21.7173, 21.7676] (0.95 confidence) | 54.9585 ns [54.8726, 55.0434] (0.95 confidence) | 222.052 ns [221.391, 222.733] (0.95 confidence) |
| d4/la&#95;stack&#95;norm2&#95;sq | 0.462637 ns [0.461733, 0.463576] (0.95 confidence) | 0.486285 ns [0.481086, 0.49246] (0.95 confidence) | 0.545274 ns [0.543567, 0.547246] (0.95 confidence) | 4.22054 ns [4.21341, 4.22791] (0.95 confidence) |
| d4/la&#95;stack&#95;solve&#95;from&#95;ldlt | 2.55737 ns [2.54503, 2.57192] (0.95 confidence) | 2.52901 ns [2.51247, 2.5479] (0.95 confidence) | 5.81531 ns [5.80181, 5.82826] (0.95 confidence) | 41.7363 ns [41.6474, 41.8243] (0.95 confidence) |
| d4/la&#95;stack&#95;solve&#95;from&#95;lu | 4.04715 ns [4.03604, 4.0589] (0.95 confidence) | 3.9379 ns [3.9308, 3.94505] (0.95 confidence) | 5.83199 ns [5.82423, 5.84028] (0.95 confidence) | 55.4299 ns [55.2979, 55.5682] (0.95 confidence) |
| d5/la&#95;stack&#95;det | 40.9326 ns [40.8844, 40.9807] (0.95 confidence) | 41.1163 ns [41.0652, 41.1706] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;det&#95;from&#95;ldlt | 0.877508 ns [0.877122, 0.877986] (0.95 confidence) | 0.881074 ns [0.879439, 0.883348] (0.95 confidence) | 0.58162 ns [0.580574, 0.582497] (0.95 confidence) | 1.23603 ns [1.23446, 1.2376] (0.95 confidence) |
| d5/la&#95;stack&#95;det&#95;from&#95;lu | 1.03418 ns [1.03344, 1.03503] (0.95 confidence) | 1.03383 ns [1.0329, 1.03483] (0.95 confidence) | 0.670221 ns [0.666679, 0.673956] (0.95 confidence) | 1.4881 ns [1.48622, 1.49048] (0.95 confidence) |
| d5/la&#95;stack&#95;det&#95;via&#95;lu | 33.4257 ns [33.3587, 33.4966] (0.95 confidence) | 33.6302 ns [33.5526, 33.7071] (0.95 confidence) | 51.4905 ns [51.3641, 51.6174] (0.95 confidence) | 191.54 ns [191.07, 192.02] (0.95 confidence) |
| d5/la&#95;stack&#95;dot | 0.836779 ns [0.833711, 0.83957] (0.95 confidence) | 0.831733 ns [0.827993, 0.835329] (0.95 confidence) | 0.764813 ns [0.763295, 0.766305] (0.95 confidence) | 2.68543 ns [2.66798, 2.70487] (0.95 confidence) |
| d5/la&#95;stack&#95;inf&#95;norm | 3.46401 ns [3.46135, 3.4665] (0.95 confidence) | 3.46346 ns [3.4591, 3.46757] (0.95 confidence) | 3.2708 ns [3.26629, 3.2751] (0.95 confidence) | 3.15991 ns [3.15843, 3.16134] (0.95 confidence) |
| d5/la&#95;stack&#95;ldlt | 24.9692 ns [24.9211, 25.0139] (0.95 confidence) | 25.0969 ns [25.0624, 25.1293] (0.95 confidence) | 25.349 ns [24.6129, 26.0987] (0.95 confidence) | 137.428 ns [137.072, 137.782] (0.95 confidence) |
| d5/la&#95;stack&#95;ldlt&#95;solve | 30.391 ns [30.341, 30.4413] (0.95 confidence) | 30.5802 ns [30.5263, 30.6359] (0.95 confidence) | 40.6196 ns [40.3922, 40.9033] (0.95 confidence) | 200.976 ns [200.255, 201.734] (0.95 confidence) |
| d5/la&#95;stack&#95;lu | 33.4374 ns [33.3708, 33.506] (0.95 confidence) | 33.4032 ns [33.3242, 33.4834] (0.95 confidence) | 50.6301 ns [50.2232, 51.1681] (0.95 confidence) | 180.438 ns [180.051, 180.816] (0.95 confidence) |
| d5/la&#95;stack&#95;lu&#95;solve | 48.4078 ns [48.2116, 48.5877] (0.95 confidence) | 49.1281 ns [49.0529, 49.2046] (0.95 confidence) | 74.7195 ns [74.4329, 74.997] (0.95 confidence) | 272.634 ns [271.856, 273.405] (0.95 confidence) |
| d5/la&#95;stack&#95;norm2&#95;sq | 0.526476 ns [0.526056, 0.526905] (0.95 confidence) | 0.526933 ns [0.526405, 0.52746] (0.95 confidence) | 0.503805 ns [0.503496, 0.504142] (0.95 confidence) | 4.13517 ns [4.13249, 4.13792] (0.95 confidence) |
| d5/la&#95;stack&#95;solve&#95;from&#95;ldlt | 3.94011 ns [3.93095, 3.94896] (0.95 confidence) | 3.90871 ns [3.90269, 3.91439] (0.95 confidence) | 9.15833 ns [9.04235, 9.33769] (0.95 confidence) | 60.2768 ns [60.118, 60.4392] (0.95 confidence) |
| d5/la&#95;stack&#95;solve&#95;from&#95;lu | 6.0822 ns [6.05978, 6.10429] (0.95 confidence) | 6.13798 ns [6.12943, 6.1468] (0.95 confidence) | 8.98041 ns [8.96122, 8.99972] (0.95 confidence) | 84.2034 ns [84.0547, 84.3511] (0.95 confidence) |
| d64/la&#95;stack&#95;det | 15387.5 ns [15288.9, 15484.2] (0.95 confidence) | 15537.1 ns [15450.8, 15624] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;det&#95;from&#95;ldlt | 23.0039 ns [22.9162, 23.0988] (0.95 confidence) | 24.3583 ns [24.3062, 24.4099] (0.95 confidence) | 8.42777 ns [8.41132, 8.44438] (0.95 confidence) | 20.9039 ns [20.8362, 20.976] (0.95 confidence) |
| d64/la&#95;stack&#95;det&#95;from&#95;lu | 23.1735 ns [23.1279, 23.2159] (0.95 confidence) | 24.3531 ns [24.3013, 24.4068] (0.95 confidence) | 8.62853 ns [8.61019, 8.64531] (0.95 confidence) | 21.2206 ns [21.1872, 21.2522] (0.95 confidence) |
| d64/la&#95;stack&#95;det&#95;via&#95;lu | 14581.7 ns [14536.4, 14627.3] (0.95 confidence) | 14845.3 ns [14796.5, 14894] (0.95 confidence) | 13016.8 ns [13002.1, 13031.3] (0.95 confidence) | 11082.9 ns [11056.6, 11112.3] (0.95 confidence) |
| d64/la&#95;stack&#95;dot | 10.8549 ns [10.8255, 10.8815] (0.95 confidence) | 11.5028 ns [11.4821, 11.523] (0.95 confidence) | 7.98078 ns [7.96332, 7.99756] (0.95 confidence) | 8.20819 ns [8.1802, 8.23787] (0.95 confidence) |
| d64/la&#95;stack&#95;inf&#95;norm | 626.923 ns [625.576, 628.207] (0.95 confidence) | 676.668 ns [674.536, 679.046] (0.95 confidence) | 1248.06 ns [1244.39, 1251.4] (0.95 confidence) | 1563.64 ns [1560.04, 1567.18] (0.95 confidence) |
| d64/la&#95;stack&#95;ldlt | 16846.1 ns [16801.6, 16898.6] (0.95 confidence) | 17169.9 ns [17148.4, 17192.6] (0.95 confidence) | 10743.9 ns [10734.7, 10752.7] (0.95 confidence) | 7015.69 ns [7001.32, 7034.73] (0.95 confidence) |
| d64/la&#95;stack&#95;ldlt&#95;solve | 19083.2 ns [19057.2, 19108.6] (0.95 confidence) | 19515.9 ns [19474.6, 19564.5] (0.95 confidence) | 12295.2 ns [12278, 12311.8] (0.95 confidence) | 8289.7 ns [8276.38, 8307.26] (0.95 confidence) |
| d64/la&#95;stack&#95;lu | 14704.2 ns [14683.9, 14725.4] (0.95 confidence) | 14949.5 ns [14919.7, 14979.8] (0.95 confidence) | 13518 ns [13474.3, 13555.6] (0.95 confidence) | 11156.9 ns [11118.2, 11197.3] (0.95 confidence) |
| d64/la&#95;stack&#95;lu&#95;solve | 17976.7 ns [17911.4, 18039.2] (0.95 confidence) | 18307.6 ns [18202.4, 18412.5] (0.95 confidence) | 13841.5 ns [13828.4, 13854.6] (0.95 confidence) | 12773.7 ns [12732.2, 12819.9] (0.95 confidence) |
| d64/la&#95;stack&#95;norm2&#95;sq | 10.8698 ns [10.8185, 10.9261] (0.95 confidence) | 11.4462 ns [11.4135, 11.4806] (0.95 confidence) | 8.49885 ns [8.48628, 8.51061] (0.95 confidence) | 6.38979 ns [6.33342, 6.44895] (0.95 confidence) |
| d64/la&#95;stack&#95;solve&#95;from&#95;ldlt | 1074.91 ns [1072.26, 1077.89] (0.95 confidence) | 1105.08 ns [1102.15, 1108.29] (0.95 confidence) | 1285.03 ns [1282.74, 1287.24] (0.95 confidence) | 1255.36 ns [1253.06, 1257.76] (0.95 confidence) |
| d64/la&#95;stack&#95;solve&#95;from&#95;lu | 2531.46 ns [2526.7, 2536.23] (0.95 confidence) | 2719.87 ns [2713.15, 2726.53] (0.95 confidence) | 790.936 ns [790.432, 791.429] (0.95 confidence) | 1684.65 ns [1682.86, 1686.39] (0.95 confidence) |
| d8/la&#95;stack&#95;det | 95.1736 ns [94.8209, 95.5046] (0.95 confidence) | 95.4649 ns [95.0827, 95.8222] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;det&#95;from&#95;ldlt | 1.31426 ns [1.31128, 1.31763] (0.95 confidence) | 1.30756 ns [1.30637, 1.30868] (0.95 confidence) | 0.906331 ns [0.901297, 0.911379] (0.95 confidence) | 2.00406 ns [1.98489, 2.02829] (0.95 confidence) |
| d8/la&#95;stack&#95;det&#95;from&#95;ldlt&#95;balanced&#95;range | 8.77584 ns [8.75176, 8.79992] (0.95 confidence) | 9.02872 ns [9.00819, 9.04959] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;det&#95;from&#95;lu | 1.46542 ns [1.46401, 1.46677] (0.95 confidence) | 1.47263 ns [1.46917, 1.47648] (0.95 confidence) | 0.988416 ns [0.984047, 0.992586] (0.95 confidence) | 2.31004 ns [2.27776, 2.34696] (0.95 confidence) |
| d8/la&#95;stack&#95;det&#95;from&#95;lu&#95;balanced&#95;range | 9.37731 ns [9.34858, 9.40608] (0.95 confidence) | 9.53277 ns [9.38008, 9.72639] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;det&#95;via&#95;lu | 89.5977 ns [89.0677, 90.1101] (0.95 confidence) | 91.1588 ns [90.6694, 91.6617] (0.95 confidence) | 140.441 ns [140.02, 140.873] (0.95 confidence) | 268.759 ns [267.679, 269.843] (0.95 confidence) |
| d8/la&#95;stack&#95;dot | 0.921649 ns [0.919907, 0.923464] (0.95 confidence) | 0.940921 ns [0.938419, 0.943604] (0.95 confidence) | 1.33039 ns [1.32865, 1.33213] (0.95 confidence) | 2.57592 ns [2.55781, 2.59931] (0.95 confidence) |
| d8/la&#95;stack&#95;inf&#95;norm | 8.40687 ns [8.38488, 8.42637] (0.95 confidence) | 8.43383 ns [8.41901, 8.44747] (0.95 confidence) | 8.0359 ns [8.02179, 8.04917] (0.95 confidence) | 8.09335 ns [8.08499, 8.10121] (0.95 confidence) |
| d8/la&#95;stack&#95;ldlt | 90.7807 ns [90.5801, 90.9944] (0.95 confidence) | 90.5629 ns [90.4554, 90.6769] (0.95 confidence) | 101.992 ns [101.338, 103.027] (0.95 confidence) | 195.409 ns [195.078, 195.749] (0.95 confidence) |
| d8/la&#95;stack&#95;ldlt&#95;ill&#95;conditioned | 90.3923 ns [90.2556, 90.5439] (0.95 confidence) | 90.2935 ns [90.2018, 90.3852] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;ldlt&#95;solve | 103.719 ns [103.505, 103.96] (0.95 confidence) | 105.515 ns [105.221, 105.816] (0.95 confidence) | 146.11 ns [144.922, 147.301] (0.95 confidence) | 266.909 ns [266.634, 267.191] (0.95 confidence) |
| d8/la&#95;stack&#95;lu | 89.8157 ns [89.4433, 90.1605] (0.95 confidence) | 90.2377 ns [90.0297, 90.4463] (0.95 confidence) | 133.03 ns [132.41, 133.614] (0.95 confidence) | 258.431 ns [257.648, 259.234] (0.95 confidence) |
| d8/la&#95;stack&#95;lu&#95;ill&#95;conditioned | 89.7279 ns [89.2354, 90.2294] (0.95 confidence) | 89.0322 ns [88.591, 89.4492] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;lu&#95;pivoting | 94.4715 ns [94.2755, 94.6681] (0.95 confidence) | 95.5521 ns [95.3155, 95.7964] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;lu&#95;solve | 130.536 ns [130.281, 130.8] (0.95 confidence) | 133.72 ns [133.437, 134.014] (0.95 confidence) | 166.473 ns [166.13, 166.819] (0.95 confidence) | 374.434 ns [373.338, 375.529] (0.95 confidence) |
| d8/la&#95;stack&#95;norm2&#95;sq | 0.671224 ns [0.668953, 0.673608] (0.95 confidence) | 0.66535 ns [0.664149, 0.66646] (0.95 confidence) | 0.938623 ns [0.935849, 0.941726] (0.95 confidence) | 4.15172 ns [4.14727, 4.15633] (0.95 confidence) |
| d8/la&#95;stack&#95;solve&#95;from&#95;ldlt | 8.20154 ns [8.19397, 8.20861] (0.95 confidence) | 8.19407 ns [8.18184, 8.2057] (0.95 confidence) | 21.7461 ns [21.6501, 21.8547] (0.95 confidence) | 69.7129 ns [69.6386, 69.7896] (0.95 confidence) |
| d8/la&#95;stack&#95;solve&#95;from&#95;lu | 13.3991 ns [13.3776, 13.4204] (0.95 confidence) | 13.4357 ns [13.4092, 13.4667] (0.95 confidence) | 16.1861 ns [16.1545, 16.2149] (0.95 confidence) | 102.928 ns [102.7, 103.162] (0.95 confidence) |
| exact&#95;d2/det | 0.438193 ns [0.437655, 0.438771] (0.95 confidence) | 0.460065 ns [0.452834, 0.470007] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;direct | 0.441632 ns [0.44026, 0.44306] (0.95 confidence) | 0.449694 ns [0.447109, 0.452912] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;direct&#95;with&#95;errbound | 1.74211 ns [1.73988, 1.7443] (0.95 confidence) | 1.77617 ns [1.77383, 1.77893] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;errbound | 1.59072 ns [1.58958, 1.59192] (0.95 confidence) | 1.62371 ns [1.6209, 1.62705] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;exact | 105.56 ns [105.366, 105.711] (0.95 confidence) | 107.55 ns [106.786, 108.249] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;exact&#95;f64&#95;result | 66.7252 ns [66.5659, 66.8851] (0.95 confidence) | 69.1351 ns [68.9925, 69.2744] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;exact&#95;rounded&#95;f64 | 68.598 ns [68.4678, 68.7306] (0.95 confidence) | 71.2031 ns [70.9945, 71.4249] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;sign&#95;exact | 2.98341 ns [2.978, 2.98886] (0.95 confidence) | 3.06474 ns [3.05713, 3.07371] (0.95 confidence) | — | — |
| exact&#95;d2/solve&#95;exact | 7101.53 ns [7090.38, 7112.74] (0.95 confidence) | 7308.61 ns [7282.98, 7340.25] (0.95 confidence) | — | — |
| exact&#95;d2/solve&#95;exact&#95;f64&#95;result | 7239.05 ns [7223.45, 7256.51] (0.95 confidence) | 7492.75 ns [7479.79, 7506.44] (0.95 confidence) | — | — |
| exact&#95;d2/solve&#95;exact&#95;rounded&#95;f64 | 7533.28 ns [7513.86, 7554.41] (0.95 confidence) | 7669.97 ns [7653.92, 7686.35] (0.95 confidence) | — | — |
| exact&#95;d3/det | 0.830687 ns [0.828107, 0.833417] (0.95 confidence) | 0.84932 ns [0.846959, 0.85175] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;direct | 0.822643 ns [0.821397, 0.824092] (0.95 confidence) | 0.846636 ns [0.844202, 0.849169] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;direct&#95;with&#95;errbound | 3.32045 ns [3.31512, 3.32655] (0.95 confidence) | 3.4032 ns [3.39678, 3.40996] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;errbound | 3.2997 ns [3.29731, 3.30252] (0.95 confidence) | 3.37902 ns [3.37142, 3.38805] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;exact | 334.11 ns [333.523, 334.726] (0.95 confidence) | 347.209 ns [346.755, 347.654] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;exact&#95;f64&#95;result | 302.755 ns [301.839, 303.884] (0.95 confidence) | 308.088 ns [307.635, 308.556] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;exact&#95;rounded&#95;f64 | 303.78 ns [303.492, 304.073] (0.95 confidence) | 311.269 ns [310.686, 311.827] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;sign&#95;exact | 4.71189 ns [4.70753, 4.71633] (0.95 confidence) | 4.80312 ns [4.79668, 4.80958] (0.95 confidence) | — | — |
| exact&#95;d3/solve&#95;exact | 31219.5 ns [31095.6, 31360.1] (0.95 confidence) | 31410.7 ns [31353.4, 31469.3] (0.95 confidence) | — | — |
| exact&#95;d3/solve&#95;exact&#95;f64&#95;result | 30895 ns [30833.8, 30955.8] (0.95 confidence) | 33106.6 ns [32764.1, 33600.4] (0.95 confidence) | — | — |
| exact&#95;d3/solve&#95;exact&#95;rounded&#95;f64 | 31328.3 ns [31283.9, 31372.9] (0.95 confidence) | 31914.8 ns [31842.1, 31989.9] (0.95 confidence) | — | — |
| exact&#95;d4/det | 3.22311 ns [3.21642, 3.23298] (0.95 confidence) | 3.259 ns [3.24465, 3.276] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;direct | 2.35399 ns [2.35274, 2.3554] (0.95 confidence) | 2.37516 ns [2.36789, 2.38652] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;direct&#95;with&#95;errbound | 6.77204 ns [6.76653, 6.77706] (0.95 confidence) | 6.81027 ns [6.80147, 6.8206] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;errbound | 6.75024 ns [6.74157, 6.75886] (0.95 confidence) | 6.79247 ns [6.77999, 6.80542] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;exact | 997.174 ns [994.009, 1000.41] (0.95 confidence) | 994.8 ns [991.392, 998.362] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;exact&#95;f64&#95;result | 952.917 ns [948.439, 957.382] (0.95 confidence) | 941.657 ns [938.485, 944.856] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;exact&#95;rounded&#95;f64 | 994.725 ns [990.92, 998.48] (0.95 confidence) | 954.803 ns [950.651, 959.021] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;sign&#95;exact | 7.93276 ns [7.92856, 7.93724] (0.95 confidence) | 7.98889 ns [7.97367, 8.00617] (0.95 confidence) | — | — |
| exact&#95;d4/solve&#95;exact | 79448 ns [79250.3, 79659.6] (0.95 confidence) | 79419.2 ns [79293.5, 79547.5] (0.95 confidence) | — | — |
| exact&#95;d4/solve&#95;exact&#95;f64&#95;result | 79484.3 ns [79394.2, 79575.5] (0.95 confidence) | 79460.4 ns [79329.4, 79591.9] (0.95 confidence) | — | — |
| exact&#95;d4/solve&#95;exact&#95;rounded&#95;f64 | 80357.2 ns [80175.3, 80565.9] (0.95 confidence) | 80490.9 ns [80334.1, 80650.4] (0.95 confidence) | — | — |
| exact&#95;d5/det | 27.7094 ns [27.6543, 27.7667] (0.95 confidence) | 27.9543 ns [27.8822, 28.0263] (0.95 confidence) | — | — |
| exact&#95;d5/det&#95;exact | 3444.94 ns [3439.22, 3450.7] (0.95 confidence) | 3493.37 ns [3487.82, 3498.91] (0.95 confidence) | — | — |
| exact&#95;d5/det&#95;exact&#95;f64&#95;result | 3398.76 ns [3391.83, 3405.8] (0.95 confidence) | 3435.43 ns [3428.97, 3442.55] (0.95 confidence) | — | — |
| exact&#95;d5/det&#95;exact&#95;rounded&#95;f64 | 3376.17 ns [3371.34, 3381] (0.95 confidence) | 3463.3 ns [3457.91, 3468.72] (0.95 confidence) | — | — |
| exact&#95;d5/det&#95;sign&#95;exact | 3364.89 ns [3359.61, 3370.23] (0.95 confidence) | 3443.45 ns [3436.91, 3450.32] (0.95 confidence) | — | — |
| exact&#95;d5/solve&#95;exact | 156623 ns [156400, 156852] (0.95 confidence) | 157629 ns [157389, 157872] (0.95 confidence) | — | — |
| exact&#95;d5/solve&#95;exact&#95;f64&#95;result | 156857 ns [156635, 157081] (0.95 confidence) | 158041 ns [157742, 158345] (0.95 confidence) | — | — |
| exact&#95;d5/solve&#95;exact&#95;rounded&#95;f64 | 157493 ns [157276, 157717] (0.95 confidence) | 159404 ns [158863, 160068] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/det&#95;exact | 1064.71 ns [1057.77, 1071.42] (0.95 confidence) | 1024.85 ns [1019.28, 1030.83] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/det&#95;sign&#95;exact | 8.07959 ns [8.06018, 8.09947] (0.95 confidence) | 7.91835 ns [7.91319, 7.92302] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/solve&#95;exact | 59655.4 ns [59580.3, 59728.7] (0.95 confidence) | 59527.9 ns [59410.5, 59646] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/solve&#95;exact&#95;f64&#95;result | 59918.6 ns [59834.9, 60005.5] (0.95 confidence) | 59465.3 ns [59360.6, 59568.6] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/solve&#95;exact&#95;rounded&#95;f64 | 60493.4 ns [60401.8, 60583.5] (0.95 confidence) | 60077.2 ns [59962, 60190.8] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/det&#95;exact | 3596 ns [3590.09, 3601.87] (0.95 confidence) | 3705.67 ns [3699.39, 3712.02] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/det&#95;sign&#95;exact | 3583.33 ns [3571.22, 3595.49] (0.95 confidence) | 3676.59 ns [3665.63, 3688.39] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/solve&#95;exact | 123691 ns [123407, 124000] (0.95 confidence) | 122501 ns [122253, 122748] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/solve&#95;exact&#95;f64&#95;result | 122495 ns [122283, 122711] (0.95 confidence) | 122956 ns [122437, 123748] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/solve&#95;exact&#95;rounded&#95;f64 | 123567 ns [123372, 123762] (0.95 confidence) | 131935 ns [130386, 133861] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/det&#95;exact | 320.927 ns [320.384, 321.487] (0.95 confidence) | 322.039 ns [321.459, 322.66] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/det&#95;sign&#95;exact | 313.249 ns [312.481, 314.043] (0.95 confidence) | 312.551 ns [312.197, 312.919] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/solve&#95;exact | 92432.9 ns [92308.6, 92560.9] (0.95 confidence) | 93469.4 ns [93382, 93557.2] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/solve&#95;exact&#95;f64&#95;result | 92313.6 ns [92182.9, 92444.5] (0.95 confidence) | 94078.7 ns [93922.8, 94245.2] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/solve&#95;exact&#95;rounded&#95;f64 | 93037.8 ns [92919, 93158.3] (0.95 confidence) | 95379.1 ns [95112.2, 95743.5] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/det&#95;exact | 330.658 ns [330.155, 331.218] (0.95 confidence) | 331.575 ns [331.123, 332.053] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/det&#95;sign&#95;exact | 382.485 ns [381.592, 383.421] (0.95 confidence) | 383.104 ns [381.834, 384.285] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/solve&#95;exact | 2639.39 ns [2635.26, 2643.52] (0.95 confidence) | 2658.68 ns [2653.54, 2663.92] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/solve&#95;exact&#95;f64&#95;result | 2674.76 ns [2667.73, 2682.9] (0.95 confidence) | 2655.98 ns [2649.46, 2662.29] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/solve&#95;exact&#95;rounded&#95;f64 | 2662.97 ns [2658.73, 2667.28] (0.95 confidence) | 2672.55 ns [2668.49, 2677.22] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/det&#95;exact | 3075.83 ns [3073.67, 3078.22] (0.95 confidence) | 3117.36 ns [3113.97, 3121.09] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/det&#95;sign&#95;exact | 150.183 ns [149.231, 151.209] (0.95 confidence) | 145.468 ns [145.085, 145.856] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/solve&#95;exact | 73943.1 ns [73801.6, 74084.5] (0.95 confidence) | 74677.1 ns [74538.4, 74813.9] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/solve&#95;exact&#95;f64&#95;result | 74789.7 ns [74640.9, 74950] (0.95 confidence) | 75858.3 ns [75719.8, 76004.6] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/solve&#95;exact&#95;rounded&#95;f64 | 75116.3 ns [74965.5, 75269.8] (0.95 confidence) | 76457.9 ns [76242.4, 76681.3] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/det&#95;exact | 8953.7 ns [8829.82, 9114.2] (0.95 confidence) | 8965.63 ns [8933.22, 8999.98] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/det&#95;sign&#95;exact | 225.109 ns [224.834, 225.397] (0.95 confidence) | 225.832 ns [225.61, 226.053] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/solve&#95;exact | 236413 ns [235835, 237009] (0.95 confidence) | 232370 ns [232162, 232582] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/solve&#95;exact&#95;f64&#95;result | 233709 ns [233360, 234053] (0.95 confidence) | 233079 ns [232683, 233453] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/solve&#95;exact&#95;rounded&#95;f64 | 232855 ns [232546, 233185] (0.95 confidence) | 234107 ns [233724, 234489] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/det&#95;exact | 21561.5 ns [21442.9, 21690.8] (0.95 confidence) | 21894.3 ns [21805.9, 21989.7] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/det&#95;sign&#95;exact | 443.25 ns [442.946, 443.556] (0.95 confidence) | 445.317 ns [444.541, 446.14] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/solve&#95;exact | 514994 ns [514483, 515517] (0.95 confidence) | 524619 ns [523548, 525744] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/solve&#95;exact&#95;f64&#95;result | 515398 ns [514737, 516102] (0.95 confidence) | 522851 ns [522158, 523528] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/solve&#95;exact&#95;rounded&#95;f64 | 516704 ns [516028, 517373] (0.95 confidence) | 525305 ns [524419, 526135] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/det&#95;exact | 55370.6 ns [55120.8, 55649.7] (0.95 confidence) | 56219.7 ns [55937, 56522.6] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/det&#95;sign&#95;exact | 54984.3 ns [54776, 55197.6] (0.95 confidence) | 55193.7 ns [54994.6, 55404.9] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/solve&#95;exact | 1.01308e+06 ns [1.01078e+06, 1.01545e+06] (0.95 confidence) | 1.02002e+06 ns [1.01876e+06, 1.02122e+06] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/solve&#95;exact&#95;f64&#95;result | 1.01018e+06 ns [1.00909e+06, 1.01127e+06] (0.95 confidence) | 1.02136e+06 ns [1.01858e+06, 1.02434e+06] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/solve&#95;exact&#95;rounded&#95;f64 | 1.01063e+06 ns [1.00927e+06, 1.01192e+06] (0.95 confidence) | 1.01802e+06 ns [1.01571e+06, 1.0201e+06] (0.95 confidence) | — | — |
| rational&#95;input&#95;d2/det&#95;big&#95;rational&#95;gaussian | 1101.74 ns [1100.24, 1103.35] (0.95 confidence) | 1106.92 ns [1104.44, 1110.22] (0.95 confidence) | — | — |
| rational&#95;input&#95;d2/det&#95;row&#95;cleared&#95;bareiss | 178.104 ns [177.862, 178.353] (0.95 confidence) | 180.318 ns [179.918, 180.746] (0.95 confidence) | — | — |
| rational&#95;input&#95;d2/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 77.0368 ns [76.9432, 77.1334] (0.95 confidence) | 77.9081 ns [77.743, 78.0791] (0.95 confidence) | — | — |
| rational&#95;input&#95;d2/solve&#95;big&#95;rational&#95;gaussian | 2105.43 ns [2101.63, 2109.5] (0.95 confidence) | 2116.51 ns [2113.51, 2119.55] (0.95 confidence) | — | — |
| rational&#95;input&#95;d2/solve&#95;row&#95;cleared&#95;bareiss | 926.427 ns [925.508, 927.351] (0.95 confidence) | 931.53 ns [930.581, 932.578] (0.95 confidence) | — | — |
| rational&#95;input&#95;d3/det&#95;big&#95;rational&#95;gaussian | 4614.22 ns [4602.18, 4626.62] (0.95 confidence) | 4628.07 ns [4615.68, 4641.06] (0.95 confidence) | — | — |
| rational&#95;input&#95;d3/det&#95;row&#95;cleared&#95;bareiss | 562.68 ns [561.488, 564.029] (0.95 confidence) | 568.154 ns [567.323, 569.022] (0.95 confidence) | — | — |
| rational&#95;input&#95;d3/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 362.361 ns [361.303, 363.502] (0.95 confidence) | 356.702 ns [355.781, 357.705] (0.95 confidence) | — | — |
| rational&#95;input&#95;d3/solve&#95;big&#95;rational&#95;gaussian | 7974.74 ns [7953.05, 7996.88] (0.95 confidence) | 8054.63 ns [8025.89, 8084.71] (0.95 confidence) | — | — |
| rational&#95;input&#95;d3/solve&#95;row&#95;cleared&#95;bareiss | 2402.19 ns [2397.85, 2406.79] (0.95 confidence) | 2425.42 ns [2414.75, 2439.47] (0.95 confidence) | — | — |
| rational&#95;input&#95;d4/det&#95;big&#95;rational&#95;gaussian | 13576.4 ns [13527.7, 13624] (0.95 confidence) | 13416.5 ns [13370.5, 13462.4] (0.95 confidence) | — | — |
| rational&#95;input&#95;d4/det&#95;row&#95;cleared&#95;bareiss | 1140.83 ns [1138.74, 1142.92] (0.95 confidence) | 1146.18 ns [1144.82, 1147.56] (0.95 confidence) | — | — |
| rational&#95;input&#95;d4/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 761.432 ns [759.768, 763.137] (0.95 confidence) | 757.253 ns [756.353, 758.199] (0.95 confidence) | — | — |
| rational&#95;input&#95;d4/solve&#95;big&#95;rational&#95;gaussian | 23725.2 ns [23664.2, 23789.3] (0.95 confidence) | 23719.6 ns [23657.8, 23780.9] (0.95 confidence) | — | — |
| rational&#95;input&#95;d4/solve&#95;row&#95;cleared&#95;bareiss | 5500.13 ns [5405.68, 5602.14] (0.95 confidence) | 5911.87 ns [5768.75, 6065.51] (0.95 confidence) | — | — |
| rational&#95;input&#95;d5/det&#95;big&#95;rational&#95;gaussian | 23205.6 ns [23111.8, 23314.3] (0.95 confidence) | 23137.7 ns [23099.9, 23177] (0.95 confidence) | — | — |
| rational&#95;input&#95;d5/det&#95;row&#95;cleared&#95;bareiss | 2055.53 ns [2051.42, 2059.61] (0.95 confidence) | 2077.81 ns [2075.3, 2080.37] (0.95 confidence) | — | — |
| rational&#95;input&#95;d5/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 1540.19 ns [1537.47, 1542.93] (0.95 confidence) | 1572.43 ns [1570.28, 1574.59] (0.95 confidence) | — | — |
| rational&#95;input&#95;d5/solve&#95;big&#95;rational&#95;gaussian | 39427.5 ns [39267.9, 39593.4] (0.95 confidence) | 39299.3 ns [39219.9, 39378.8] (0.95 confidence) | — | — |
| rational&#95;input&#95;d5/solve&#95;row&#95;cleared&#95;bareiss | 9350.24 ns [9331.51, 9369.2] (0.95 confidence) | 9464.07 ns [9448.98, 9479.42] (0.95 confidence) | — | — |
| rational&#95;input&#95;d6/det&#95;big&#95;rational&#95;gaussian | 55982.3 ns [55847.8, 56113.5] (0.95 confidence) | 56994.4 ns [56875, 57112.8] (0.95 confidence) | — | — |
| rational&#95;input&#95;d6/det&#95;row&#95;cleared&#95;bareiss | 3601.81 ns [3595.16, 3608.46] (0.95 confidence) | 3640.98 ns [3631.9, 3652.73] (0.95 confidence) | — | — |
| rational&#95;input&#95;d6/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 2947.81 ns [2942.41, 2953.28] (0.95 confidence) | 2983.67 ns [2978.91, 2988.4] (0.95 confidence) | — | — |
| rational&#95;input&#95;d6/solve&#95;big&#95;rational&#95;gaussian | 89277.2 ns [89057.2, 89501.4] (0.95 confidence) | 90181.9 ns [89971.9, 90397.4] (0.95 confidence) | — | — |
| rational&#95;input&#95;d6/solve&#95;row&#95;cleared&#95;bareiss | 17424.9 ns [17367.6, 17488.8] (0.95 confidence) | 17803.2 ns [17772.1, 17834.1] (0.95 confidence) | — | — |
| rational&#95;input&#95;d7/det&#95;big&#95;rational&#95;gaussian | 108611 ns [108354, 108873] (0.95 confidence) | 110246 ns [109941, 110600] (0.95 confidence) | — | — |
| rational&#95;input&#95;d7/det&#95;row&#95;cleared&#95;bareiss | 5850.32 ns [5839.08, 5863.38] (0.95 confidence) | 5889.54 ns [5872.79, 5907.38] (0.95 confidence) | — | — |
| rational&#95;input&#95;d7/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 5003.98 ns [4992.98, 5015.6] (0.95 confidence) | 5000.2 ns [4991.23, 5009.15] (0.95 confidence) | — | — |
| rational&#95;input&#95;d7/solve&#95;big&#95;rational&#95;gaussian | 159309 ns [158679, 160020] (0.95 confidence) | 162486 ns [162175, 162790] (0.95 confidence) | — | — |
| rational&#95;input&#95;d7/solve&#95;row&#95;cleared&#95;bareiss | 28193.8 ns [28128.1, 28260.2] (0.95 confidence) | 28918.5 ns [28872.2, 28964] (0.95 confidence) | — | — |
| rational&#95;input&#95;d8/det&#95;big&#95;rational&#95;gaussian | 187569 ns [187131, 188019] (0.95 confidence) | 190172 ns [189753, 190586] (0.95 confidence) | — | — |
| rational&#95;input&#95;d8/det&#95;row&#95;cleared&#95;bareiss | 9020.66 ns [8897.13, 9189.34] (0.95 confidence) | 8883.24 ns [8871.22, 8894.84] (0.95 confidence) | — | — |
| rational&#95;input&#95;d8/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 7985.48 ns [7971.39, 7999.63] (0.95 confidence) | 7942.5 ns [7923.98, 7962.36] (0.95 confidence) | — | — |
| rational&#95;input&#95;d8/solve&#95;big&#95;rational&#95;gaussian | 262705 ns [262251, 263168] (0.95 confidence) | 268790 ns [268225, 269354] (0.95 confidence) | — | — |
| rational&#95;input&#95;d8/solve&#95;row&#95;cleared&#95;bareiss | 40881 ns [40819.1, 40944.1] (0.95 confidence) | 41790.7 ns [41637, 42001.2] (0.95 confidence) | — | — |

## Median (ns)

| Case | la-stack baseline | la-stack current | nalgebra | faer |
| --- | --- | --- | --- | --- |
| d16/la&#95;stack&#95;det | 379.272 ns [378.451, 380.366] (0.95 confidence) | 401.093 ns [398.153, 404.658] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;det&#95;from&#95;ldlt | 2.89094 ns [2.88425, 2.90038] (0.95 confidence) | 2.96543 ns [2.96339, 2.96871] (0.95 confidence) | 1.70051 ns [1.69525, 1.70445] (0.95 confidence) | 4.3533 ns [4.33624, 4.36853] (0.95 confidence) |
| d16/la&#95;stack&#95;det&#95;from&#95;lu | 3.28159 ns [3.27565, 3.28649] (0.95 confidence) | 3.35471 ns [3.35195, 3.35876] (0.95 confidence) | 1.77532 ns [1.76963, 1.77906] (0.95 confidence) | 4.82318 ns [4.76454, 4.8634] (0.95 confidence) |
| d16/la&#95;stack&#95;det&#95;via&#95;lu | 360.619 ns [359.716, 361.689] (0.95 confidence) | 377.718 ns [374.568, 381.126] (0.95 confidence) | 451.034 ns [449.913, 451.762] (0.95 confidence) | 666.409 ns [663.273, 673.084] (0.95 confidence) |
| d16/la&#95;stack&#95;dot | 2.33315 ns [2.33218, 2.33422] (0.95 confidence) | 2.36585 ns [2.36324, 2.36963] (0.95 confidence) | 1.84921 ns [1.84831, 1.84993] (0.95 confidence) | 3.31767 ns [3.31433, 3.32792] (0.95 confidence) |
| d16/la&#95;stack&#95;inf&#95;norm | 32.4514 ns [32.3366, 32.4761] (0.95 confidence) | 33.5202 ns [33.4712, 33.5765] (0.95 confidence) | 32.5791 ns [32.5201, 32.605] (0.95 confidence) | 32.2307 ns [32.2168, 32.2844] (0.95 confidence) |
| d16/la&#95;stack&#95;ldlt | 338.907 ns [338.69, 339.485] (0.95 confidence) | 345.474 ns [344.655, 346.116] (0.95 confidence) | 422.509 ns [421.62, 422.946] (0.95 confidence) | 418.131 ns [417.556, 418.894] (0.95 confidence) |
| d16/la&#95;stack&#95;ldlt&#95;solve | 382.109 ns [381.947, 382.505] (0.95 confidence) | 390.634 ns [390.107, 391.056] (0.95 confidence) | 628.019 ns [626.887, 629.04] (0.95 confidence) | 612.089 ns [610.79, 613.119] (0.95 confidence) |
| d16/la&#95;stack&#95;lu | 350.086 ns [348.492, 351.41] (0.95 confidence) | 368.462 ns [364.863, 371.298] (0.95 confidence) | 464.597 ns [463.776, 465.843] (0.95 confidence) | 633.175 ns [631.494, 640.576] (0.95 confidence) |
| d16/la&#95;stack&#95;lu&#95;solve | 627.457 ns [613.21, 635.397] (0.95 confidence) | 750.027 ns [649.969, 761.453] (0.95 confidence) | 589.011 ns [586.474, 592.341] (0.95 confidence) | 938.943 ns [932.02, 944.127] (0.95 confidence) |
| d16/la&#95;stack&#95;norm2&#95;sq | 2.12017 ns [2.11718, 2.1233] (0.95 confidence) | 2.14796 ns [2.14589, 2.15195] (0.95 confidence) | 1.73558 ns [1.73353, 1.7364] (0.95 confidence) | 4.09384 ns [4.09297, 4.09492] (0.95 confidence) |
| d16/la&#95;stack&#95;solve&#95;from&#95;ldlt | 27.7905 ns [27.7687, 27.8075] (0.95 confidence) | 29.016 ns [28.9085, 29.1021] (0.95 confidence) | 119.414 ns [118.977, 119.984] (0.95 confidence) | 182.626 ns [182.299, 182.82] (0.95 confidence) |
| d16/la&#95;stack&#95;solve&#95;from&#95;lu | 191.074 ns [190.623, 191.465] (0.95 confidence) | 209.188 ns [208.682, 209.517] (0.95 confidence) | 93.6387 ns [93.394, 93.7612] (0.95 confidence) | 257.154 ns [256.905, 257.7] (0.95 confidence) |
| d2/la&#95;stack&#95;det | 0.572876 ns [0.572777, 0.57299] (0.95 confidence) | 0.573987 ns [0.573647, 0.574175] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;det&#95;from&#95;ldlt | 0.496238 ns [0.496192, 0.496347] (0.95 confidence) | 0.496604 ns [0.496449, 0.496839] (0.95 confidence) | 0.438546 ns [0.438, 0.4387] (0.95 confidence) | 0.570374 ns [0.569774, 0.57083] (0.95 confidence) |
| d2/la&#95;stack&#95;det&#95;from&#95;lu | 0.520624 ns [0.520547, 0.520724] (0.95 confidence) | 0.520987 ns [0.520651, 0.52145] (0.95 confidence) | 0.45059 ns [0.450332, 0.451806] (0.95 confidence) | 0.743834 ns [0.743654, 0.744151] (0.95 confidence) |
| d2/la&#95;stack&#95;det&#95;via&#95;lu | 1.77086 ns [1.76992, 1.77177] (0.95 confidence) | 1.7435 ns [1.7374, 1.75313] (0.95 confidence) | 0.822672 ns [0.821112, 0.823663] (0.95 confidence) | 110.718 ns [110.261, 111.408] (0.95 confidence) |
| d2/la&#95;stack&#95;dot | 0.612921 ns [0.612366, 0.613432] (0.95 confidence) | 0.610449 ns [0.606927, 0.613366] (0.95 confidence) | 0.598043 ns [0.593009, 0.599721] (0.95 confidence) | 2.47931 ns [2.47867, 2.47984] (0.95 confidence) |
| d2/la&#95;stack&#95;inf&#95;norm | 0.703132 ns [0.702861, 0.703304] (0.95 confidence) | 0.703434 ns [0.703269, 0.703603] (0.95 confidence) | 0.559833 ns [0.559717, 0.559959] (0.95 confidence) | 0.695855 ns [0.695335, 0.696852] (0.95 confidence) |
| d2/la&#95;stack&#95;ldlt | 6.89114 ns [6.87985, 6.91103] (0.95 confidence) | 6.55722 ns [6.54552, 6.57972] (0.95 confidence) | 1.72375 ns [1.72274, 1.7245] (0.95 confidence) | 87.7743 ns [87.4516, 88.1834] (0.95 confidence) |
| d2/la&#95;stack&#95;ldlt&#95;solve | 9.60901 ns [9.59092, 9.63616] (0.95 confidence) | 9.49351 ns [9.47354, 9.50791] (0.95 confidence) | 2.62395 ns [2.61761, 2.62944] (0.95 confidence) | 125.726 ns [125.433, 126.029] (0.95 confidence) |
| d2/la&#95;stack&#95;lu | 1.56569 ns [1.56284, 1.5679] (0.95 confidence) | 1.56693 ns [1.56549, 1.56814] (0.95 confidence) | 1.68189 ns [1.67891, 1.68622] (0.95 confidence) | 108.948 ns [107.575, 110.523] (0.95 confidence) |
| d2/la&#95;stack&#95;lu&#95;solve | 2.01731 ns [2.01633, 2.01805] (0.95 confidence) | 2.10359 ns [2.08908, 2.11912] (0.95 confidence) | 4.53894 ns [4.53141, 4.54661] (0.95 confidence) | 154.45 ns [153.658, 155.846] (0.95 confidence) |
| d2/la&#95;stack&#95;norm2&#95;sq | 0.395429 ns [0.395177, 0.395569] (0.95 confidence) | 0.395034 ns [0.394428, 0.395421] (0.95 confidence) | 0.389089 ns [0.388885, 0.389318] (0.95 confidence) | 4.1265 ns [4.12534, 4.12826] (0.95 confidence) |
| d2/la&#95;stack&#95;solve&#95;from&#95;ldlt | 1.23867 ns [1.2383, 1.23904] (0.95 confidence) | 1.23932 ns [1.23805, 1.24032] (0.95 confidence) | 1.24972 ns [1.24834, 1.25148] (0.95 confidence) | 39.463 ns [39.3332, 39.5216] (0.95 confidence) |
| d2/la&#95;stack&#95;solve&#95;from&#95;lu | 1.24468 ns [1.24094, 1.24904] (0.95 confidence) | 1.26331 ns [1.25802, 1.26633] (0.95 confidence) | 2.894 ns [2.89321, 2.89545] (0.95 confidence) | 49.3097 ns [49.1857, 49.4274] (0.95 confidence) |
| d3/la&#95;stack&#95;det | 1.33913 ns [1.33765, 1.34325] (0.95 confidence) | 1.33991 ns [1.33668, 1.34289] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;det&#95;from&#95;ldlt | 0.612948 ns [0.611813, 0.614127] (0.95 confidence) | 0.595467 ns [0.595158, 0.595737] (0.95 confidence) | 0.496208 ns [0.495686, 0.49675] (0.95 confidence) | 0.760555 ns [0.759415, 0.763128] (0.95 confidence) |
| d3/la&#95;stack&#95;det&#95;from&#95;lu | 0.709547 ns [0.708629, 0.710965] (0.95 confidence) | 0.693946 ns [0.693803, 0.694034] (0.95 confidence) | 0.483905 ns [0.48194, 0.485024] (0.95 confidence) | 1.01064 ns [1.00857, 1.01164] (0.95 confidence) |
| d3/la&#95;stack&#95;det&#95;via&#95;lu | 8.79971 ns [8.76753, 8.84205] (0.95 confidence) | 8.68592 ns [8.67467, 8.69586] (0.95 confidence) | 12.61 ns [12.603, 12.6328] (0.95 confidence) | 138.806 ns [138.328, 139.246] (0.95 confidence) |
| d3/la&#95;stack&#95;dot | 0.681535 ns [0.680378, 0.683296] (0.95 confidence) | 0.671023 ns [0.669925, 0.672805] (0.95 confidence) | 0.706024 ns [0.703341, 0.70755] (0.95 confidence) | 3.0367 ns [3.02574, 3.04921] (0.95 confidence) |
| d3/la&#95;stack&#95;inf&#95;norm | 1.38741 ns [1.3859, 1.38884] (0.95 confidence) | 1.36927 ns [1.36835, 1.37046] (0.95 confidence) | 1.19754 ns [1.19605, 1.19858] (0.95 confidence) | 1.34749 ns [1.34389, 1.3515] (0.95 confidence) |
| d3/la&#95;stack&#95;ldlt | 13.4093 ns [13.3568, 13.4652] (0.95 confidence) | 13.0128 ns [12.9757, 13.0532] (0.95 confidence) | 3.95925 ns [3.95144, 3.97094] (0.95 confidence) | 99.0856 ns [98.9316, 99.2756] (0.95 confidence) |
| d3/la&#95;stack&#95;ldlt&#95;solve | 9.80848 ns [9.78745, 9.84487] (0.95 confidence) | 9.74358 ns [9.71499, 9.76146] (0.95 confidence) | 5.75807 ns [5.73394, 5.76789] (0.95 confidence) | 147.084 ns [146.664, 147.499] (0.95 confidence) |
| d3/la&#95;stack&#95;lu | 8.68803 ns [8.62311, 8.7214] (0.95 confidence) | 8.38217 ns [8.36225, 8.39415] (0.95 confidence) | 12.7922 ns [12.7034, 12.8895] (0.95 confidence) | 139.091 ns [138.443, 139.316] (0.95 confidence) |
| d3/la&#95;stack&#95;lu&#95;solve | 9.98135 ns [9.97281, 9.9906] (0.95 confidence) | 9.93879 ns [9.92881, 9.95117] (0.95 confidence) | 32.9342 ns [32.8554, 33.057] (0.95 confidence) | 191.85 ns [191.206, 192.814] (0.95 confidence) |
| d3/la&#95;stack&#95;norm2&#95;sq | 0.507066 ns [0.505861, 0.508281] (0.95 confidence) | 0.498095 ns [0.496577, 0.499413] (0.95 confidence) | 0.452066 ns [0.451062, 0.453544] (0.95 confidence) | 4.26639 ns [4.25754, 4.28818] (0.95 confidence) |
| d3/la&#95;stack&#95;solve&#95;from&#95;ldlt | 1.82795 ns [1.82373, 1.83132] (0.95 confidence) | 1.78095 ns [1.78017, 1.78173] (0.95 confidence) | 3.00925 ns [2.99857, 3.02052] (0.95 confidence) | 41.4065 ns [41.2795, 41.5182] (0.95 confidence) |
| d3/la&#95;stack&#95;solve&#95;from&#95;lu | 2.13084 ns [2.12446, 2.13602] (0.95 confidence) | 2.05702 ns [2.05559, 2.05822] (0.95 confidence) | 4.69852 ns [4.68887, 4.70656] (0.95 confidence) | 53.7473 ns [53.5314, 53.9728] (0.95 confidence) |
| d32/la&#95;stack&#95;det | 2175.6 ns [2167.52, 2183.12] (0.95 confidence) | 2253.83 ns [2248.47, 2257.61] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;det&#95;from&#95;ldlt | 7.80178 ns [7.78797, 7.81116] (0.95 confidence) | 8.16296 ns [8.14651, 8.1759] (0.95 confidence) | 3.0888 ns [3.08514, 3.09175] (0.95 confidence) | 8.60059 ns [8.5721, 8.63071] (0.95 confidence) |
| d32/la&#95;stack&#95;det&#95;from&#95;lu | 9.09263 ns [9.08045, 9.11029] (0.95 confidence) | 9.55091 ns [9.52895, 9.59579] (0.95 confidence) | 3.07434 ns [3.06716, 3.07811] (0.95 confidence) | 8.78609 ns [8.73977, 8.82193] (0.95 confidence) |
| d32/la&#95;stack&#95;det&#95;via&#95;lu | 2042.01 ns [2039.57, 2044.43] (0.95 confidence) | 2097.15 ns [2093.25, 2104.1] (0.95 confidence) | 2269.4 ns [2260.57, 2275.1] (0.95 confidence) | 2391.84 ns [2384.89, 2394.43] (0.95 confidence) |
| d32/la&#95;stack&#95;dot | 4.02217 ns [4.0174, 4.02843] (0.95 confidence) | 4.09795 ns [4.09615, 4.10263] (0.95 confidence) | 4.00127 ns [3.9967, 4.0067] (0.95 confidence) | 4.87016 ns [4.86241, 4.87819] (0.95 confidence) |
| d32/la&#95;stack&#95;inf&#95;norm | 127.933 ns [127.73, 128.284] (0.95 confidence) | 133.153 ns [133.06, 133.385] (0.95 confidence) | 160.906 ns [158.963, 162.054] (0.95 confidence) | 162.784 ns [162.45, 162.981] (0.95 confidence) |
| d32/la&#95;stack&#95;ldlt | 2384.2 ns [2382.57, 2386.09] (0.95 confidence) | 2431.22 ns [2428.58, 2433.78] (0.95 confidence) | 2092.05 ns [2082.88, 2102.77] (0.95 confidence) | 1399.49 ns [1390.69, 1408.7] (0.95 confidence) |
| d32/la&#95;stack&#95;ldlt&#95;solve | 2707.43 ns [2704.76, 2712.19] (0.95 confidence) | 2768.09 ns [2762.52, 2774.48] (0.95 confidence) | 2810.35 ns [2784.92, 2816.19] (0.95 confidence) | 1861.98 ns [1859.82, 1865.27] (0.95 confidence) |
| d32/la&#95;stack&#95;lu | 2138.35 ns [2134.27, 2143.88] (0.95 confidence) | 2226.92 ns [2218.65, 2230.7] (0.95 confidence) | 2120.89 ns [2118.28, 2124.36] (0.95 confidence) | 2366.95 ns [2363.16, 2370.37] (0.95 confidence) |
| d32/la&#95;stack&#95;lu&#95;solve | 2762.99 ns [2758.9, 2767.72] (0.95 confidence) | 2829.62 ns [2819.66, 2837.22] (0.95 confidence) | 2722.74 ns [2720.98, 2729.21] (0.95 confidence) | 3023.87 ns [3017.73, 3030.19] (0.95 confidence) |
| d32/la&#95;stack&#95;norm2&#95;sq | 3.9658 ns [3.96092, 3.97436] (0.95 confidence) | 4.1013 ns [4.08692, 4.10916] (0.95 confidence) | 3.62609 ns [3.61571, 3.63556] (0.95 confidence) | 4.21289 ns [4.21225, 4.21392] (0.95 confidence) |
| d32/la&#95;stack&#95;solve&#95;from&#95;ldlt | 311.264 ns [310.487, 312.295] (0.95 confidence) | 320.496 ns [319.713, 321.007] (0.95 confidence) | 548.658 ns [546.83, 549.753] (0.95 confidence) | 469.109 ns [468.904, 469.465] (0.95 confidence) |
| d32/la&#95;stack&#95;solve&#95;from&#95;lu | 666.961 ns [665.864, 668.213] (0.95 confidence) | 735.74 ns [734.025, 737.377] (0.95 confidence) | 328.48 ns [328.138, 329.019] (0.95 confidence) | 639.316 ns [638.049, 640.345] (0.95 confidence) |
| d4/la&#95;stack&#95;det | 2.59079 ns [2.58535, 2.60559] (0.95 confidence) | 2.5437 ns [2.52985, 2.55742] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;det&#95;from&#95;ldlt | 0.787562 ns [0.786422, 0.791756] (0.95 confidence) | 0.77159 ns [0.770293, 0.772868] (0.95 confidence) | 0.499227 ns [0.498689, 0.4997] (0.95 confidence) | 1.00377 ns [1.00289, 1.00498] (0.95 confidence) |
| d4/la&#95;stack&#95;det&#95;from&#95;lu | 0.886058 ns [0.884854, 0.887513] (0.95 confidence) | 0.870765 ns [0.870079, 0.872809] (0.95 confidence) | 0.590961 ns [0.588626, 0.592034] (0.95 confidence) | 1.25584 ns [1.25455, 1.2569] (0.95 confidence) |
| d4/la&#95;stack&#95;det&#95;via&#95;lu | 14.5493 ns [14.5287, 14.5816] (0.95 confidence) | 14.3681 ns [14.3312, 14.3953] (0.95 confidence) | 31.2562 ns [31.1699, 31.3552] (0.95 confidence) | 164.743 ns [164.348, 165.425] (0.95 confidence) |
| d4/la&#95;stack&#95;dot | 0.691713 ns [0.688447, 0.693459] (0.95 confidence) | 0.68448 ns [0.681956, 0.688119] (0.95 confidence) | 0.666976 ns [0.665682, 0.668425] (0.95 confidence) | 2.73962 ns [2.7005, 2.74687] (0.95 confidence) |
| d4/la&#95;stack&#95;inf&#95;norm | 2.39495 ns [2.39335, 2.39861] (0.95 confidence) | 2.37107 ns [2.36854, 2.374] (0.95 confidence) | 2.08273 ns [2.08038, 2.08506] (0.95 confidence) | 2.08911 ns [2.08339, 2.10665] (0.95 confidence) |
| d4/la&#95;stack&#95;ldlt | 18.2285 ns [18.2088, 18.2493] (0.95 confidence) | 17.9138 ns [17.8725, 17.9542] (0.95 confidence) | 7.92229 ns [7.9045, 7.95379] (0.95 confidence) | 124.238 ns [123.806, 124.89] (0.95 confidence) |
| d4/la&#95;stack&#95;ldlt&#95;solve | 18.1373 ns [18.0834, 18.1765] (0.95 confidence) | 17.9777 ns [17.9146, 18.0244] (0.95 confidence) | 11.5452 ns [11.4995, 11.5868] (0.95 confidence) | 166.26 ns [165.641, 166.955] (0.95 confidence) |
| d4/la&#95;stack&#95;lu | 13.8975 ns [13.8774, 13.9233] (0.95 confidence) | 13.7862 ns [13.7637, 13.8035] (0.95 confidence) | 30.0526 ns [30.0148, 30.1072] (0.95 confidence) | 158.857 ns [157.98, 159.484] (0.95 confidence) |
| d4/la&#95;stack&#95;lu&#95;solve | 22.1488 ns [22.1185, 22.1688] (0.95 confidence) | 21.7341 ns [21.7035, 21.776] (0.95 confidence) | 55.0106 ns [54.9217, 55.0615] (0.95 confidence) | 221.887 ns [221.678, 222.673] (0.95 confidence) |
| d4/la&#95;stack&#95;norm2&#95;sq | 0.461115 ns [0.460671, 0.462203] (0.95 confidence) | 0.477166 ns [0.474712, 0.481989] (0.95 confidence) | 0.54413 ns [0.543035, 0.545759] (0.95 confidence) | 4.21735 ns [4.20825, 4.22393] (0.95 confidence) |
| d4/la&#95;stack&#95;solve&#95;from&#95;ldlt | 2.54384 ns [2.53093, 2.55676] (0.95 confidence) | 2.49308 ns [2.48506, 2.50136] (0.95 confidence) | 5.82202 ns [5.81192, 5.83516] (0.95 confidence) | 41.7311 ns [41.6519, 41.8151] (0.95 confidence) |
| d4/la&#95;stack&#95;solve&#95;from&#95;lu | 4.0318 ns [4.02105, 4.04677] (0.95 confidence) | 3.9445 ns [3.93473, 3.94977] (0.95 confidence) | 5.83637 ns [5.82832, 5.84141] (0.95 confidence) | 55.3476 ns [55.2688, 55.4685] (0.95 confidence) |
| d5/la&#95;stack&#95;det | 40.9715 ns [40.8612, 41.0168] (0.95 confidence) | 41.1138 ns [41.0703, 41.1471] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;det&#95;from&#95;ldlt | 0.877386 ns [0.877116, 0.877632] (0.95 confidence) | 0.878577 ns [0.877941, 0.879622] (0.95 confidence) | 0.583159 ns [0.582517, 0.58354] (0.95 confidence) | 1.23782 ns [1.23634, 1.23875] (0.95 confidence) |
| d5/la&#95;stack&#95;det&#95;from&#95;lu | 1.03418 ns [1.03346, 1.03456] (0.95 confidence) | 1.033 ns [1.03162, 1.03452] (0.95 confidence) | 0.66136 ns [0.657915, 0.668257] (0.95 confidence) | 1.48662 ns [1.48628, 1.48713] (0.95 confidence) |
| d5/la&#95;stack&#95;det&#95;via&#95;lu | 33.4972 ns [33.4588, 33.5131] (0.95 confidence) | 33.6162 ns [33.5415, 33.6868] (0.95 confidence) | 51.4294 ns [51.3175, 51.524] (0.95 confidence) | 191.177 ns [190.659, 191.624] (0.95 confidence) |
| d5/la&#95;stack&#95;dot | 0.843249 ns [0.8426, 0.843671] (0.95 confidence) | 0.842959 ns [0.839884, 0.843939] (0.95 confidence) | 0.764448 ns [0.763663, 0.76533] (0.95 confidence) | 2.63811 ns [2.63731, 2.64365] (0.95 confidence) |
| d5/la&#95;stack&#95;inf&#95;norm | 3.46955 ns [3.46784, 3.47043] (0.95 confidence) | 3.46893 ns [3.46542, 3.47116] (0.95 confidence) | 3.2796 ns [3.27431, 3.28241] (0.95 confidence) | 3.16084 ns [3.15981, 3.16145] (0.95 confidence) |
| d5/la&#95;stack&#95;ldlt | 25.0049 ns [24.944, 25.0526] (0.95 confidence) | 25.1376 ns [25.11, 25.154] (0.95 confidence) | 25.2354 ns [24.3772, 26.3912] (0.95 confidence) | 137.376 ns [136.651, 137.838] (0.95 confidence) |
| d5/la&#95;stack&#95;ldlt&#95;solve | 30.4296 ns [30.3614, 30.4628] (0.95 confidence) | 30.5698 ns [30.4626, 30.6452] (0.95 confidence) | 40.2172 ns [40.2036, 40.2558] (0.95 confidence) | 199.907 ns [199.061, 200.506] (0.95 confidence) |
| d5/la&#95;stack&#95;lu | 33.5064 ns [33.317, 33.5916] (0.95 confidence) | 33.2397 ns [33.1367, 33.6301] (0.95 confidence) | 50.0776 ns [49.8946, 50.5527] (0.95 confidence) | 180.626 ns [180.122, 181.097] (0.95 confidence) |
| d5/la&#95;stack&#95;lu&#95;solve | 48.7153 ns [48.6718, 48.7504] (0.95 confidence) | 49.1031 ns [49.0318, 49.1659] (0.95 confidence) | 74.629 ns [74.3476, 74.9906] (0.95 confidence) | 272.864 ns [272.009, 273.512] (0.95 confidence) |
| d5/la&#95;stack&#95;norm2&#95;sq | 0.525661 ns [0.525259, 0.526131] (0.95 confidence) | 0.527524 ns [0.525701, 0.528065] (0.95 confidence) | 0.503771 ns [0.503572, 0.503944] (0.95 confidence) | 4.13532 ns [4.13289, 4.13641] (0.95 confidence) |
| d5/la&#95;stack&#95;solve&#95;from&#95;ldlt | 3.94677 ns [3.93389, 3.96169] (0.95 confidence) | 3.91732 ns [3.91421, 3.92183] (0.95 confidence) | 9.03899 ns [9.00572, 9.05456] (0.95 confidence) | 60.1795 ns [60.1129, 60.3863] (0.95 confidence) |
| d5/la&#95;stack&#95;solve&#95;from&#95;lu | 6.08472 ns [6.06153, 6.11231] (0.95 confidence) | 6.12656 ns [6.11897, 6.14238] (0.95 confidence) | 8.9966 ns [8.97481, 9.01403] (0.95 confidence) | 84.2009 ns [84.0597, 84.3299] (0.95 confidence) |
| d64/la&#95;stack&#95;det | 15419.7 ns [15317.1, 15517.3] (0.95 confidence) | 15607.4 ns [15525.1, 15681.8] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;det&#95;from&#95;ldlt | 22.9949 ns [22.8735, 23.0413] (0.95 confidence) | 24.3353 ns [24.2589, 24.4077] (0.95 confidence) | 8.43186 ns [8.40615, 8.45378] (0.95 confidence) | 20.8183 ns [20.7795, 20.8946] (0.95 confidence) |
| d64/la&#95;stack&#95;det&#95;from&#95;lu | 23.231 ns [23.1874, 23.2604] (0.95 confidence) | 24.3422 ns [24.295, 24.3796] (0.95 confidence) | 8.64464 ns [8.63183, 8.65888] (0.95 confidence) | 21.2372 ns [21.2039, 21.2809] (0.95 confidence) |
| d64/la&#95;stack&#95;det&#95;via&#95;lu | 14596.6 ns [14537.5, 14694.7] (0.95 confidence) | 14869.9 ns [14823.1, 14944.5] (0.95 confidence) | 13020.4 ns [13001.5, 13032.7] (0.95 confidence) | 11036.6 ns [11027.5, 11059] (0.95 confidence) |
| d64/la&#95;stack&#95;dot | 10.8944 ns [10.8602, 10.9125] (0.95 confidence) | 11.5095 ns [11.4922, 11.5295] (0.95 confidence) | 7.99232 ns [7.96211, 8.00956] (0.95 confidence) | 8.17384 ns [8.10427, 8.23606] (0.95 confidence) |
| d64/la&#95;stack&#95;inf&#95;norm | 627.793 ns [626.723, 629.574] (0.95 confidence) | 674.478 ns [672.691, 676.417] (0.95 confidence) | 1251.53 ns [1248.47, 1255.67] (0.95 confidence) | 1568.14 ns [1562.6, 1570.88] (0.95 confidence) |
| d64/la&#95;stack&#95;ldlt | 16774.7 ns [16761.3, 16786.9] (0.95 confidence) | 17153.1 ns [17130, 17175.4] (0.95 confidence) | 10747.5 ns [10740.8, 10757.6] (0.95 confidence) | 7001.9 ns [6997.07, 7009.98] (0.95 confidence) |
| d64/la&#95;stack&#95;ldlt&#95;solve | 19075.4 ns [19037.7, 19105.9] (0.95 confidence) | 19478.7 ns [19456.8, 19523] (0.95 confidence) | 12299.4 ns [12291.5, 12313.2] (0.95 confidence) | 8273.44 ns [8262.32, 8283.77] (0.95 confidence) |
| d64/la&#95;stack&#95;lu | 14682.2 ns [14663.6, 14702.9] (0.95 confidence) | 14924.1 ns [14891.1, 14997.4] (0.95 confidence) | 13604.1 ns [13575, 13616.4] (0.95 confidence) | 11098.4 ns [11037.4, 11167.1] (0.95 confidence) |
| d64/la&#95;stack&#95;lu&#95;solve | 18014.7 ns [17849.1, 18166.4] (0.95 confidence) | 18280.4 ns [18138, 18489.9] (0.95 confidence) | 13835.1 ns [13822.9, 13852.2] (0.95 confidence) | 12678.8 ns [12660.2, 12699.3] (0.95 confidence) |
| d64/la&#95;stack&#95;norm2&#95;sq | 10.8165 ns [10.8016, 10.8318] (0.95 confidence) | 11.3915 ns [11.3645, 11.4378] (0.95 confidence) | 8.52183 ns [8.50728, 8.5323] (0.95 confidence) | 6.21011 ns [6.20027, 6.26274] (0.95 confidence) |
| d64/la&#95;stack&#95;solve&#95;from&#95;ldlt | 1074.14 ns [1072.88, 1076.54] (0.95 confidence) | 1103.32 ns [1099.25, 1105.99] (0.95 confidence) | 1289.2 ns [1286.32, 1290.23] (0.95 confidence) | 1252.76 ns [1249.61, 1255.61] (0.95 confidence) |
| d64/la&#95;stack&#95;solve&#95;from&#95;lu | 2530.99 ns [2524.36, 2541.55] (0.95 confidence) | 2719.71 ns [2714.93, 2728.93] (0.95 confidence) | 791.158 ns [790.924, 791.33] (0.95 confidence) | 1685.24 ns [1683.32, 1686.69] (0.95 confidence) |
| d8/la&#95;stack&#95;det | 95.5337 ns [95.422, 95.6859] (0.95 confidence) | 96.0299 ns [95.7172, 96.1538] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;det&#95;from&#95;ldlt | 1.30998 ns [1.30856, 1.31208] (0.95 confidence) | 1.30758 ns [1.30593, 1.30941] (0.95 confidence) | 0.913542 ns [0.896885, 0.918783] (0.95 confidence) | 1.97881 ns [1.97447, 1.98117] (0.95 confidence) |
| d8/la&#95;stack&#95;det&#95;from&#95;ldlt&#95;balanced&#95;range | 8.73914 ns [8.68084, 8.82264] (0.95 confidence) | 9.02925 ns [8.99874, 9.0541] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;det&#95;from&#95;lu | 1.46489 ns [1.46432, 1.46598] (0.95 confidence) | 1.46681 ns [1.46544, 1.46859] (0.95 confidence) | 0.997928 ns [0.990364, 0.998579] (0.95 confidence) | 2.23371 ns [2.23129, 2.23987] (0.95 confidence) |
| d8/la&#95;stack&#95;det&#95;from&#95;lu&#95;balanced&#95;range | 9.40991 ns [9.36982, 9.41723] (0.95 confidence) | 9.36441 ns [9.19918, 9.46503] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;det&#95;via&#95;lu | 90.0003 ns [89.4881, 90.7238] (0.95 confidence) | 91.1377 ns [90.9825, 91.3288] (0.95 confidence) | 140.38 ns [140.256, 140.682] (0.95 confidence) | 268.382 ns [267.16, 270.555] (0.95 confidence) |
| d8/la&#95;stack&#95;dot | 0.920693 ns [0.919859, 0.921825] (0.95 confidence) | 0.938248 ns [0.935731, 0.940223] (0.95 confidence) | 1.33079 ns [1.32954, 1.3322] (0.95 confidence) | 2.54972 ns [2.54228, 2.57066] (0.95 confidence) |
| d8/la&#95;stack&#95;inf&#95;norm | 8.43703 ns [8.42209, 8.45037] (0.95 confidence) | 8.46336 ns [8.44575, 8.47159] (0.95 confidence) | 8.05139 ns [8.03854, 8.06795] (0.95 confidence) | 8.10969 ns [8.10081, 8.11662] (0.95 confidence) |
| d8/la&#95;stack&#95;ldlt | 90.4337 ns [90.2549, 90.6032] (0.95 confidence) | 90.4841 ns [90.4109, 90.6392] (0.95 confidence) | 101.46 ns [100.991, 102.101] (0.95 confidence) | 195.344 ns [195.102, 195.453] (0.95 confidence) |
| d8/la&#95;stack&#95;ldlt&#95;ill&#95;conditioned | 90.2188 ns [90.1175, 90.3155] (0.95 confidence) | 90.2188 ns [90.1527, 90.4098] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;ldlt&#95;solve | 103.441 ns [103.298, 103.708] (0.95 confidence) | 105.555 ns [105.215, 105.799] (0.95 confidence) | 147.45 ns [141.922, 149.357] (0.95 confidence) | 266.801 ns [266.3, 267.128] (0.95 confidence) |
| d8/la&#95;stack&#95;lu | 90.0189 ns [89.8607, 90.3413] (0.95 confidence) | 90.1784 ns [89.9849, 90.362] (0.95 confidence) | 133.545 ns [132.756, 133.916] (0.95 confidence) | 257.893 ns [257.102, 259.004] (0.95 confidence) |
| d8/la&#95;stack&#95;lu&#95;ill&#95;conditioned | 89.6288 ns [88.9996, 90.2728] (0.95 confidence) | 89.6705 ns [89.4697, 89.882] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;lu&#95;pivoting | 94.4911 ns [94.3366, 94.8131] (0.95 confidence) | 95.2265 ns [95.0854, 95.5085] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;lu&#95;solve | 130.471 ns [130.304, 130.733] (0.95 confidence) | 133.665 ns [133.426, 133.861] (0.95 confidence) | 166.317 ns [166.02, 166.604] (0.95 confidence) | 374.912 ns [372.746, 376.096] (0.95 confidence) |
| d8/la&#95;stack&#95;norm2&#95;sq | 0.668475 ns [0.667505, 0.669308] (0.95 confidence) | 0.66656 ns [0.665109, 0.668505] (0.95 confidence) | 0.932707 ns [0.931802, 0.933498] (0.95 confidence) | 4.14741 ns [4.14266, 4.15215] (0.95 confidence) |
| d8/la&#95;stack&#95;solve&#95;from&#95;ldlt | 8.21381 ns [8.20853, 8.21669] (0.95 confidence) | 8.20724 ns [8.19176, 8.21578] (0.95 confidence) | 21.6993 ns [21.6244, 21.7363] (0.95 confidence) | 69.6593 ns [69.6017, 69.6999] (0.95 confidence) |
| d8/la&#95;stack&#95;solve&#95;from&#95;lu | 13.3948 ns [13.3776, 13.4128] (0.95 confidence) | 13.4109 ns [13.3987, 13.4503] (0.95 confidence) | 16.2178 ns [16.186, 16.2658] (0.95 confidence) | 102.443 ns [102.246, 102.805] (0.95 confidence) |
| exact&#95;d2/det | 0.437607 ns [0.437261, 0.43802] (0.95 confidence) | 0.450079 ns [0.449262, 0.451482] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;direct | 0.438443 ns [0.437972, 0.439661] (0.95 confidence) | 0.44589 ns [0.444961, 0.44701] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;direct&#95;with&#95;errbound | 1.74342 ns [1.74005, 1.74601] (0.95 confidence) | 1.77356 ns [1.77148, 1.77481] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;errbound | 1.58811 ns [1.5874, 1.58956] (0.95 confidence) | 1.61984 ns [1.61872, 1.62252] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;exact | 105.681 ns [105.555, 105.714] (0.95 confidence) | 108.898 ns [108.622, 109.16] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;exact&#95;f64&#95;result | 66.6816 ns [66.3984, 67.107] (0.95 confidence) | 69.2162 ns [69.1264, 69.3796] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;exact&#95;rounded&#95;f64 | 68.5379 ns [68.3688, 68.7476] (0.95 confidence) | 70.9349 ns [70.8286, 71.0535] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;sign&#95;exact | 2.97719 ns [2.97444, 2.98201] (0.95 confidence) | 3.05619 ns [3.0517, 3.05976] (0.95 confidence) | — | — |
| exact&#95;d2/solve&#95;exact | 7101.58 ns [7076.55, 7114.39] (0.95 confidence) | 7287.66 ns [7269.82, 7303.67] (0.95 confidence) | — | — |
| exact&#95;d2/solve&#95;exact&#95;f64&#95;result | 7211.4 ns [7198.57, 7221.83] (0.95 confidence) | 7486.5 ns [7474.23, 7503.96] (0.95 confidence) | — | — |
| exact&#95;d2/solve&#95;exact&#95;rounded&#95;f64 | 7521.46 ns [7503.04, 7541.24] (0.95 confidence) | 7667.74 ns [7653.73, 7681.71] (0.95 confidence) | — | — |
| exact&#95;d3/det | 0.824679 ns [0.82185, 0.826576] (0.95 confidence) | 0.844887 ns [0.84331, 0.847052] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;direct | 0.820217 ns [0.819999, 0.820613] (0.95 confidence) | 0.841538 ns [0.839892, 0.843676] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;direct&#95;with&#95;errbound | 3.31432 ns [3.31121, 3.31644] (0.95 confidence) | 3.39753 ns [3.38846, 3.40442] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;errbound | 3.2961 ns [3.29439, 3.29739] (0.95 confidence) | 3.37007 ns [3.36717, 3.37234] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;exact | 333.2 ns [332.537, 334.062] (0.95 confidence) | 347.228 ns [346.674, 347.694] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;exact&#95;f64&#95;result | 301.192 ns [300.546, 301.931] (0.95 confidence) | 307.89 ns [307.366, 308.632] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;exact&#95;rounded&#95;f64 | 303.616 ns [303.313, 304.146] (0.95 confidence) | 311.716 ns [310.916, 312.184] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;sign&#95;exact | 4.70867 ns [4.70694, 4.71219] (0.95 confidence) | 4.80631 ns [4.79645, 4.81409] (0.95 confidence) | — | — |
| exact&#95;d3/solve&#95;exact | 31151.7 ns [30971.3, 31260.6] (0.95 confidence) | 31421.8 ns [31340, 31484.2] (0.95 confidence) | — | — |
| exact&#95;d3/solve&#95;exact&#95;f64&#95;result | 30868.5 ns [30801.3, 30901] (0.95 confidence) | 32665 ns [32367.3, 32932.5] (0.95 confidence) | — | — |
| exact&#95;d3/solve&#95;exact&#95;rounded&#95;f64 | 31301.8 ns [31224.5, 31382.9] (0.95 confidence) | 31889.1 ns [31809.7, 31976.6] (0.95 confidence) | — | — |
| exact&#95;d4/det | 3.21436 ns [3.21379, 3.21571] (0.95 confidence) | 3.22868 ns [3.22382, 3.23502] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;direct | 2.35275 ns [2.35217, 2.35392] (0.95 confidence) | 2.36516 ns [2.36235, 2.36771] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;direct&#95;with&#95;errbound | 6.77784 ns [6.77585, 6.77884] (0.95 confidence) | 6.79239 ns [6.78985, 6.79734] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;errbound | 6.75278 ns [6.74407, 6.75666] (0.95 confidence) | 6.78473 ns [6.76819, 6.79724] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;exact | 995.982 ns [992.447, 998.628] (0.95 confidence) | 994.463 ns [989.647, 997.38] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;exact&#95;f64&#95;result | 953.334 ns [946.67, 956.673] (0.95 confidence) | 942.457 ns [939.248, 946.725] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;exact&#95;rounded&#95;f64 | 997.858 ns [992.58, 1001.8] (0.95 confidence) | 953.063 ns [947.833, 958.796] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;sign&#95;exact | 7.92726 ns [7.92438, 7.92858] (0.95 confidence) | 7.96825 ns [7.95719, 7.9766] (0.95 confidence) | — | — |
| exact&#95;d4/solve&#95;exact | 79180.8 ns [79084.4, 79417.8] (0.95 confidence) | 79433 ns [79240.9, 79558.3] (0.95 confidence) | — | — |
| exact&#95;d4/solve&#95;exact&#95;f64&#95;result | 79524.9 ns [79375.2, 79623.1] (0.95 confidence) | 79595.9 ns [79338.6, 79742.8] (0.95 confidence) | — | — |
| exact&#95;d4/solve&#95;exact&#95;rounded&#95;f64 | 80130.8 ns [80036.9, 80288.4] (0.95 confidence) | 80440.7 ns [80332.6, 80657.4] (0.95 confidence) | — | — |
| exact&#95;d5/det | 27.7235 ns [27.6692, 27.766] (0.95 confidence) | 27.9531 ns [27.866, 28.1329] (0.95 confidence) | — | — |
| exact&#95;d5/det&#95;exact | 3443.38 ns [3437.65, 3453.31] (0.95 confidence) | 3496.21 ns [3485.51, 3503.35] (0.95 confidence) | — | — |
| exact&#95;d5/det&#95;exact&#95;f64&#95;result | 3399.92 ns [3390.22, 3405.48] (0.95 confidence) | 3425.94 ns [3422.05, 3429.93] (0.95 confidence) | — | — |
| exact&#95;d5/det&#95;exact&#95;rounded&#95;f64 | 3371.59 ns [3365.66, 3384.72] (0.95 confidence) | 3464.17 ns [3457.26, 3472] (0.95 confidence) | — | — |
| exact&#95;d5/det&#95;sign&#95;exact | 3366.38 ns [3356.93, 3373.04] (0.95 confidence) | 3443.89 ns [3431.94, 3449.15] (0.95 confidence) | — | — |
| exact&#95;d5/solve&#95;exact | 156461 ns [156212, 156628] (0.95 confidence) | 157713 ns [157301, 158021] (0.95 confidence) | — | — |
| exact&#95;d5/solve&#95;exact&#95;f64&#95;result | 156745 ns [156557, 157013] (0.95 confidence) | 158057 ns [157661, 158457] (0.95 confidence) | — | — |
| exact&#95;d5/solve&#95;exact&#95;rounded&#95;f64 | 157433 ns [157193, 157646] (0.95 confidence) | 158819 ns [158422, 159037] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/det&#95;exact | 1076.5 ns [1068.35, 1087.49] (0.95 confidence) | 1013.79 ns [1010.21, 1019.36] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/det&#95;sign&#95;exact | 8.05776 ns [8.03223, 8.08877] (0.95 confidence) | 7.92417 ns [7.92002, 7.92615] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/solve&#95;exact | 59638.5 ns [59580.5, 59746.1] (0.95 confidence) | 59591.1 ns [59495.9, 59686.2] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/solve&#95;exact&#95;f64&#95;result | 59811.3 ns [59764.6, 59928.3] (0.95 confidence) | 59446.4 ns [59318, 59581.2] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/solve&#95;exact&#95;rounded&#95;f64 | 60476.9 ns [60337.1, 60612.8] (0.95 confidence) | 60139.9 ns [60076.1, 60207.2] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/det&#95;exact | 3596.17 ns [3588.46, 3605.22] (0.95 confidence) | 3703.79 ns [3696.35, 3713.21] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/det&#95;sign&#95;exact | 3590.06 ns [3575.21, 3605.64] (0.95 confidence) | 3670.86 ns [3651.49, 3686.03] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/solve&#95;exact | 123340 ns [123059, 123696] (0.95 confidence) | 122365 ns [122128, 122554] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/solve&#95;exact&#95;f64&#95;result | 122512 ns [122299, 122707] (0.95 confidence) | 122484 ns [122326, 122754] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/solve&#95;exact&#95;rounded&#95;f64 | 123587 ns [123413, 123844] (0.95 confidence) | 129890 ns [129489, 130259] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/det&#95;exact | 320.36 ns [319.794, 321.303] (0.95 confidence) | 321.245 ns [320.736, 322.073] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/det&#95;sign&#95;exact | 313.534 ns [312.574, 314.258] (0.95 confidence) | 312.382 ns [312.027, 312.8] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/solve&#95;exact | 92374.7 ns [92255.2, 92526.3] (0.95 confidence) | 93423.1 ns [93319.1, 93553.1] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/solve&#95;exact&#95;f64&#95;result | 92292.6 ns [92189.6, 92467.7] (0.95 confidence) | 93917.6 ns [93809.9, 94023.6] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/solve&#95;exact&#95;rounded&#95;f64 | 92993.5 ns [92788.4, 93187.4] (0.95 confidence) | 95011.1 ns [94877.7, 95275] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/det&#95;exact | 330.396 ns [329.774, 330.732] (0.95 confidence) | 331.719 ns [331.215, 332.041] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/det&#95;sign&#95;exact | 381.617 ns [380.891, 382.738] (0.95 confidence) | 384.296 ns [383.549, 385.741] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/solve&#95;exact | 2640.57 ns [2634.42, 2648.9] (0.95 confidence) | 2661.62 ns [2657.04, 2664.72] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/solve&#95;exact&#95;f64&#95;result | 2666.04 ns [2662.54, 2670.71] (0.95 confidence) | 2664.43 ns [2658.27, 2668.62] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/solve&#95;exact&#95;rounded&#95;f64 | 2661.42 ns [2656.71, 2666.63] (0.95 confidence) | 2669.25 ns [2665.89, 2672.32] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/det&#95;exact | 3074.61 ns [3073.48, 3076.5] (0.95 confidence) | 3113.29 ns [3110.98, 3116.17] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/det&#95;sign&#95;exact | 149.951 ns [148.821, 150.311] (0.95 confidence) | 145.447 ns [145.019, 145.8] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/solve&#95;exact | 73898.7 ns [73785.4, 74041.7] (0.95 confidence) | 74691.2 ns [74552.8, 74879.2] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/solve&#95;exact&#95;f64&#95;result | 74636 ns [74528.3, 74739.7] (0.95 confidence) | 75679.1 ns [75588.4, 75805] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/solve&#95;exact&#95;rounded&#95;f64 | 74998.9 ns [74909.1, 75146.1] (0.95 confidence) | 76119.2 ns [75828.6, 76475.4] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/det&#95;exact | 8735.9 ns [8706.82, 8747.5] (0.95 confidence) | 8916.02 ns [8890.29, 8952.93] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/det&#95;sign&#95;exact | 225.041 ns [224.843, 225.166] (0.95 confidence) | 225.856 ns [225.386, 226.057] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/solve&#95;exact | 236603 ns [234555, 237500] (0.95 confidence) | 232425 ns [232059, 232620] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/solve&#95;exact&#95;f64&#95;result | 233964 ns [233377, 234178] (0.95 confidence) | 233273 ns [232820, 233741] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/solve&#95;exact&#95;rounded&#95;f64 | 232670 ns [232372, 232989] (0.95 confidence) | 234259 ns [233846, 234719] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/det&#95;exact | 21465.1 ns [21329.4, 21524.8] (0.95 confidence) | 21767.8 ns [21688.3, 21846.1] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/det&#95;sign&#95;exact | 443.083 ns [442.756, 443.478] (0.95 confidence) | 444.827 ns [443.47, 445.771] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/solve&#95;exact | 514860 ns [514153, 515277] (0.95 confidence) | 524102 ns [523441, 525063] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/solve&#95;exact&#95;f64&#95;result | 515182 ns [514699, 515836] (0.95 confidence) | 522636 ns [522137, 522990] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/solve&#95;exact&#95;rounded&#95;f64 | 516974 ns [515876, 517511] (0.95 confidence) | 525773 ns [524838, 527087] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/det&#95;exact | 55095.3 ns [54756.3, 55363.3] (0.95 confidence) | 55825.8 ns [55684.2, 56051.1] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/det&#95;sign&#95;exact | 54917.3 ns [54608.2, 55126.3] (0.95 confidence) | 55009.3 ns [54885, 55329.4] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/solve&#95;exact | 1.01098e+06 ns [1.00922e+06, 1.01757e+06] (0.95 confidence) | 1.02105e+06 ns [1.01928e+06, 1.0227e+06] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/solve&#95;exact&#95;f64&#95;result | 1.0099e+06 ns [1.00872e+06, 1.01108e+06] (0.95 confidence) | 1.02094e+06 ns [1.01881e+06, 1.02267e+06] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/solve&#95;exact&#95;rounded&#95;f64 | 1.0115e+06 ns [1.01054e+06, 1.0123e+06] (0.95 confidence) | 1.01892e+06 ns [1.01752e+06, 1.02096e+06] (0.95 confidence) | — | — |
| rational&#95;input&#95;d2/det&#95;big&#95;rational&#95;gaussian | 1101.62 ns [1099.34, 1102.96] (0.95 confidence) | 1104.07 ns [1102.04, 1107.44] (0.95 confidence) | — | — |
| rational&#95;input&#95;d2/det&#95;row&#95;cleared&#95;bareiss | 178.295 ns [177.825, 178.437] (0.95 confidence) | 179.604 ns [179.355, 180.059] (0.95 confidence) | — | — |
| rational&#95;input&#95;d2/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 76.9844 ns [76.8428, 77.153] (0.95 confidence) | 77.8888 ns [77.737, 78.0925] (0.95 confidence) | — | — |
| rational&#95;input&#95;d2/solve&#95;big&#95;rational&#95;gaussian | 2101.82 ns [2099, 2106.11] (0.95 confidence) | 2115.58 ns [2112.65, 2118.49] (0.95 confidence) | — | — |
| rational&#95;input&#95;d2/solve&#95;row&#95;cleared&#95;bareiss | 926.348 ns [924.941, 927.762] (0.95 confidence) | 931.625 ns [930.645, 932.404] (0.95 confidence) | — | — |
| rational&#95;input&#95;d3/det&#95;big&#95;rational&#95;gaussian | 4607.84 ns [4591.84, 4618.04] (0.95 confidence) | 4618.71 ns [4606.42, 4629.37] (0.95 confidence) | — | — |
| rational&#95;input&#95;d3/det&#95;row&#95;cleared&#95;bareiss | 561.084 ns [560.655, 561.946] (0.95 confidence) | 567.794 ns [566.479, 568.868] (0.95 confidence) | — | — |
| rational&#95;input&#95;d3/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 361.636 ns [360.405, 362.795] (0.95 confidence) | 356.108 ns [355.037, 356.669] (0.95 confidence) | — | — |
| rational&#95;input&#95;d3/solve&#95;big&#95;rational&#95;gaussian | 7954.99 ns [7935.36, 7982.04] (0.95 confidence) | 8041.07 ns [8011.87, 8073.69] (0.95 confidence) | — | — |
| rational&#95;input&#95;d3/solve&#95;row&#95;cleared&#95;bareiss | 2399.62 ns [2397.75, 2401.91] (0.95 confidence) | 2415.2 ns [2410.25, 2420.77] (0.95 confidence) | — | — |
| rational&#95;input&#95;d4/det&#95;big&#95;rational&#95;gaussian | 13591.5 ns [13511.7, 13635.6] (0.95 confidence) | 13455.6 ns [13419.5, 13507.7] (0.95 confidence) | — | — |
| rational&#95;input&#95;d4/det&#95;row&#95;cleared&#95;bareiss | 1142.12 ns [1138.56, 1143.7] (0.95 confidence) | 1145.98 ns [1143.13, 1147.91] (0.95 confidence) | — | — |
| rational&#95;input&#95;d4/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 761.999 ns [759.33, 763.949] (0.95 confidence) | 756.302 ns [755.292, 757.201] (0.95 confidence) | — | — |
| rational&#95;input&#95;d4/solve&#95;big&#95;rational&#95;gaussian | 23663.8 ns [23615.8, 23716.4] (0.95 confidence) | 23719.3 ns [23645.2, 23790.2] (0.95 confidence) | — | — |
| rational&#95;input&#95;d4/solve&#95;row&#95;cleared&#95;bareiss | 5289.48 ns [5280.14, 5298.04] (0.95 confidence) | 5772.1 ns [5438.53, 5833.49] (0.95 confidence) | — | — |
| rational&#95;input&#95;d5/det&#95;big&#95;rational&#95;gaussian | 23230.8 ns [23150.3, 23279.4] (0.95 confidence) | 23137 ns [23102.7, 23164.8] (0.95 confidence) | — | — |
| rational&#95;input&#95;d5/det&#95;row&#95;cleared&#95;bareiss | 2055.03 ns [2051.28, 2060.1] (0.95 confidence) | 2077 ns [2074.49, 2079.9] (0.95 confidence) | — | — |
| rational&#95;input&#95;d5/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 1538.17 ns [1535.12, 1542.43] (0.95 confidence) | 1571.96 ns [1570.15, 1575.37] (0.95 confidence) | — | — |
| rational&#95;input&#95;d5/solve&#95;big&#95;rational&#95;gaussian | 39295.2 ns [39086.4, 39466.5] (0.95 confidence) | 39284.3 ns [39220, 39359.5] (0.95 confidence) | — | — |
| rational&#95;input&#95;d5/solve&#95;row&#95;cleared&#95;bareiss | 9337.06 ns [9322.82, 9367.19] (0.95 confidence) | 9466.7 ns [9442.65, 9484.9] (0.95 confidence) | — | — |
| rational&#95;input&#95;d6/det&#95;big&#95;rational&#95;gaussian | 56010.1 ns [55774.6, 56174.1] (0.95 confidence) | 56994.6 ns [56766, 57163.9] (0.95 confidence) | — | — |
| rational&#95;input&#95;d6/det&#95;row&#95;cleared&#95;bareiss | 3600.36 ns [3590.3, 3610.64] (0.95 confidence) | 3634.44 ns [3625.5, 3643.79] (0.95 confidence) | — | — |
| rational&#95;input&#95;d6/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 2941.38 ns [2937.33, 2948.15] (0.95 confidence) | 2981.62 ns [2978.38, 2987.35] (0.95 confidence) | — | — |
| rational&#95;input&#95;d6/solve&#95;big&#95;rational&#95;gaussian | 89102.5 ns [88929.7, 89313.1] (0.95 confidence) | 90044.2 ns [89810.5, 90232.2] (0.95 confidence) | — | — |
| rational&#95;input&#95;d6/solve&#95;row&#95;cleared&#95;bareiss | 17413.4 ns [17310.6, 17491.9] (0.95 confidence) | 17813.1 ns [17775.4, 17853] (0.95 confidence) | — | — |
| rational&#95;input&#95;d7/det&#95;big&#95;rational&#95;gaussian | 108117 ns [107947, 108683] (0.95 confidence) | 110131 ns [109794, 110370] (0.95 confidence) | — | — |
| rational&#95;input&#95;d7/det&#95;row&#95;cleared&#95;bareiss | 5842.88 ns [5824.17, 5856.91] (0.95 confidence) | 5874.12 ns [5853.19, 5886.59] (0.95 confidence) | — | — |
| rational&#95;input&#95;d7/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 5000.46 ns [4988.05, 5011.27] (0.95 confidence) | 4997.34 ns [4990.86, 5005.29] (0.95 confidence) | — | — |
| rational&#95;input&#95;d7/solve&#95;big&#95;rational&#95;gaussian | 158783 ns [158480, 159372] (0.95 confidence) | 162794 ns [162278, 163066] (0.95 confidence) | — | — |
| rational&#95;input&#95;d7/solve&#95;row&#95;cleared&#95;bareiss | 28160.2 ns [28040.7, 28265.3] (0.95 confidence) | 28929.7 ns [28904, 28967.5] (0.95 confidence) | — | — |
| rational&#95;input&#95;d8/det&#95;big&#95;rational&#95;gaussian | 187037 ns [186707, 187403] (0.95 confidence) | 190295 ns [189907, 190514] (0.95 confidence) | — | — |
| rational&#95;input&#95;d8/det&#95;row&#95;cleared&#95;bareiss | 8821.33 ns [8793.62, 8847.83] (0.95 confidence) | 8890.39 ns [8876.42, 8901.57] (0.95 confidence) | — | — |
| rational&#95;input&#95;d8/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 7984.02 ns [7972.89, 8000.1] (0.95 confidence) | 7938.35 ns [7921.87, 7948.98] (0.95 confidence) | — | — |
| rational&#95;input&#95;d8/solve&#95;big&#95;rational&#95;gaussian | 261996 ns [261648, 262655] (0.95 confidence) | 268285 ns [267897, 268994] (0.95 confidence) | — | — |
| rational&#95;input&#95;d8/solve&#95;row&#95;cleared&#95;bareiss | 40878.2 ns [40802.8, 40944] (0.95 confidence) | 41625.7 ns [41560.5, 41757.9] (0.95 confidence) | — | — |

## Recorded measurement provenance

These observations belong to the measured phases. Rerendering does not recapture the publisher's environment.

### Baseline phase

| Field | Recorded value |
| --- | --- |
| Release | v0.4.6 |
| Revision | 0215e65763f8eecb4bcad8bf29de4aca67531aa1 |
| Measurement mode | tag |
| Original source SHA-256 | d5054a85bd6b163fa142cb94ffc01a217d033893dd9ea8b5e422a442dcbc5074 |
| Measured source SHA-256 | c4756ecf0cde859b365b4c082e12fccf66f21b7e4752c323eb10e9a74bc4f6fa |
| Common harness SHA-256 | bdb17ee72a8f6482ab66fbd92781bba3137a02c528302068f1c3119dc547b17d |
| Fingerprint schema | research-repo-tools/files/v1 |
| OS | macOS-27.0.1-arm64-arm-64bit-Mach-O |
| Architecture | arm64 |
| CPU | Apple M4 Max (arm64) |
| Physical cores | 16 |
| Logical threads | 16 |
| Memory (bytes) | 68719476736 |
| Tool: cargo | cargo 1.99.0 (5f94df478 2026-08-27) |
| Tool: just | just 1.58.0 |
| Tool: rustc | rustc 1.99.0 (b940084d7 2026-09-28)<br>binary: rustc<br>commit-hash: b940084d7eb6a299eb4bfeb8e34901bc051e7ac4<br>commit-date: 2026-09-28<br>host: aarch64-apple-darwin<br>release: 1.99.0<br>LLVM version: 23.1.1 |
| Dependency: criterion | 0.8.2 |
| Dependency: faer | 0.24.4 |
| Dependency: nalgebra | 0.35.0 |
| Environment overrides (names only) | CARGO&#95;TARGET&#95;DIR, CRITERION&#95;HOME, RUSTFLAGS, RUSTUP&#95;TOOLCHAIN |
| Measurement command (argv) | &#91;&quot;just&quot;, &quot;--justfile&quot;, &quot;tooling/performance.just&quot;, &quot;--working-directory&quot;, &quot;.&quot;, &quot;baseline-all&quot;&#93; |
| Independent correctness gate | passed |
| Gate command (argv) | &#91;&quot;just&quot;, &quot;--justfile&quot;, &quot;tooling/performance.just&quot;, &quot;--working-directory&quot;, &quot;.&quot;, &quot;gate&quot;&#93; |

### Current phase

| Field | Recorded value |
| --- | --- |
| Release | v0.5.0 |
| Revision | 51868ea4725c1a4790bba4c59341d2113b18b2bb |
| Measurement mode | working-tree |
| Original source SHA-256 | 1d876cd1476f52283a2f9dbf1bdbf46baec54b86834115f1eddefee601ad6b68 |
| Measured source SHA-256 | 1d876cd1476f52283a2f9dbf1bdbf46baec54b86834115f1eddefee601ad6b68 |
| Common harness SHA-256 | bdb17ee72a8f6482ab66fbd92781bba3137a02c528302068f1c3119dc547b17d |
| Fingerprint schema | research-repo-tools/files/v1 |
| OS | macOS-27.0.1-arm64-arm-64bit-Mach-O |
| Architecture | arm64 |
| CPU | Apple M4 Max (arm64) |
| Physical cores | 16 |
| Logical threads | 16 |
| Memory (bytes) | 68719476736 |
| Tool: cargo | cargo 1.99.0 (5f94df478 2026-08-27) |
| Tool: just | just 1.58.0 |
| Tool: rustc | rustc 1.99.0 (b940084d7 2026-09-28)<br>binary: rustc<br>commit-hash: b940084d7eb6a299eb4bfeb8e34901bc051e7ac4<br>commit-date: 2026-09-28<br>host: aarch64-apple-darwin<br>release: 1.99.0<br>LLVM version: 23.1.1 |
| Dependency: criterion | 0.8.2 |
| Dependency: faer | 0.24.4 |
| Dependency: nalgebra | 0.35.0 |
| Environment overrides (names only) | CARGO&#95;TARGET&#95;DIR, CRITERION&#95;HOME, RUSTFLAGS, RUSTUP&#95;TOOLCHAIN |
| Measurement command (argv) | &#91;&quot;just&quot;, &quot;--justfile&quot;, &quot;tooling/performance.just&quot;, &quot;--working-directory&quot;, &quot;.&quot;, &quot;current-all&quot;&#93; |
| Independent correctness gate | passed |
| Gate command (argv) | &#91;&quot;just&quot;, &quot;--justfile&quot;, &quot;tooling/performance.just&quot;, &quot;--working-directory&quot;, &quot;.&quot;, &quot;gate&quot;&#93; |
