<!-- research-repo-tools/run-report/v2 -->
# Benchmark timings

Run: `26cb584886988dc4a95c717e49f27147bfc68a3edf9b495e944ccf751d0cd413`

Evidence: [run payload](run.json) · [evidence envelope](evidence.json).

Intervals are marginal timing intervals. No significance or ratio interval is inferred.

- la-stack baseline: measured in **baseline** phase (v0.5.0), revision `5d8f5b4ea1c0bd708cd748b86f426960b7597ab7`.
- la-stack current: measured in **current** phase (v0.5.1), revision `eb6f071607c058527f6cb6a9bf86bb307c93ab02`.
- nalgebra: measured in **baseline** phase (v0.5.0), revision `5d8f5b4ea1c0bd708cd748b86f426960b7597ab7`.
- faer: measured in **baseline** phase (v0.5.0), revision `5d8f5b4ea1c0bd708cd748b86f426960b7597ab7`.

## Mean (ns)

| Case | la-stack baseline | la-stack current | nalgebra | faer |
| --- | --- | --- | --- | --- |
| d16/la&#95;stack&#95;det | 432.469 ns [430.85, 434.114] (0.95 confidence) | 437.002 ns [434.599, 439.443] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;det&#95;from&#95;ldlt | 2.90166 ns [2.8958, 2.90769] (0.95 confidence) | 2.91988 ns [2.91681, 2.92322] (0.95 confidence) | 1.68498 ns [1.67984, 1.68998] (0.95 confidence) | 4.36875 ns [4.35419, 4.38386] (0.95 confidence) |
| d16/la&#95;stack&#95;det&#95;from&#95;lu | 3.29887 ns [3.29433, 3.30389] (0.95 confidence) | 3.30704 ns [3.30356, 3.31067] (0.95 confidence) | 1.76792 ns [1.76184, 1.7738] (0.95 confidence) | 4.73194 ns [4.71044, 4.75406] (0.95 confidence) |
| d16/la&#95;stack&#95;det&#95;via&#95;lu | 414.274 ns [411.546, 417.013] (0.95 confidence) | 421.585 ns [419.219, 423.987] (0.95 confidence) | 456.433 ns [454.978, 457.827] (0.95 confidence) | 638.639 ns [635.078, 642.406] (0.95 confidence) |
| d16/la&#95;stack&#95;dot | 2.3743 ns [2.3492, 2.40641] (0.95 confidence) | 2.35592 ns [2.34842, 2.3642] (0.95 confidence) | 1.85114 ns [1.84813, 1.85426] (0.95 confidence) | 3.26553 ns [3.25954, 3.27282] (0.95 confidence) |
| d16/la&#95;stack&#95;inf&#95;norm | 32.0114 ns [31.9385, 32.0688] (0.95 confidence) | 32.959 ns [32.9022, 33.0249] (0.95 confidence) | 32.2007 ns [32.1547, 32.2431] (0.95 confidence) | 32.1679 ns [32.1187, 32.2437] (0.95 confidence) |
| d16/la&#95;stack&#95;ldlt | 341.417 ns [340.871, 342.047] (0.95 confidence) | 341.817 ns [341.382, 342.26] (0.95 confidence) | 416.84 ns [415.6, 418.099] (0.95 confidence) | 406.349 ns [405.247, 407.521] (0.95 confidence) |
| d16/la&#95;stack&#95;ldlt&#95;solve | 384.197 ns [383.703, 384.758] (0.95 confidence) | 384.37 ns [384, 384.763] (0.95 confidence) | 640.119 ns [638.248, 642.154] (0.95 confidence) | 603.036 ns [601.461, 604.667] (0.95 confidence) |
| d16/la&#95;stack&#95;lu | 409.202 ns [407.245, 411.202] (0.95 confidence) | 399.34 ns [397.163, 401.571] (0.95 confidence) | 462.711 ns [461.337, 464.122] (0.95 confidence) | 618.121 ns [614.263, 622.061] (0.95 confidence) |
| d16/la&#95;stack&#95;lu&#95;solve | 657.55 ns [655.399, 659.732] (0.95 confidence) | 665.676 ns [663.951, 667.409] (0.95 confidence) | 591.906 ns [590.665, 593.121] (0.95 confidence) | 876.838 ns [873.041, 880.783] (0.95 confidence) |
| d16/la&#95;stack&#95;norm2&#95;sq | 2.11166 ns [2.10732, 2.11616] (0.95 confidence) | 2.11555 ns [2.11069, 2.12022] (0.95 confidence) | 1.73143 ns [1.72795, 1.73516] (0.95 confidence) | 4.10491 ns [4.09957, 4.11065] (0.95 confidence) |
| d16/la&#95;stack&#95;solve&#95;from&#95;ldlt | 27.8592 ns [27.8277, 27.893] (0.95 confidence) | 27.9577 ns [27.9166, 28.003] (0.95 confidence) | 119.104 ns [118.781, 119.449] (0.95 confidence) | 180.01 ns [179.756, 180.273] (0.95 confidence) |
| d16/la&#95;stack&#95;solve&#95;from&#95;lu | 190.929 ns [190.353, 191.617] (0.95 confidence) | 201.063 ns [200.672, 201.451] (0.95 confidence) | 93.7524 ns [93.5974, 93.9088] (0.95 confidence) | 246.609 ns [246.109, 247.113] (0.95 confidence) |
| d2/la&#95;stack&#95;det | 0.590869 ns [0.589738, 0.592083] (0.95 confidence) | 0.581027 ns [0.580077, 0.582046] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;det&#95;from&#95;ldlt | 0.510214 ns [0.507717, 0.512974] (0.95 confidence) | 0.510046 ns [0.507561, 0.512768] (0.95 confidence) | 0.446264 ns [0.445735, 0.446802] (0.95 confidence) | 0.575151 ns [0.572903, 0.577217] (0.95 confidence) |
| d2/la&#95;stack&#95;det&#95;from&#95;lu | 0.560553 ns [0.559914, 0.561257] (0.95 confidence) | 0.553751 ns [0.552483, 0.555069] (0.95 confidence) | 0.470425 ns [0.468357, 0.472816] (0.95 confidence) | 0.772088 ns [0.770581, 0.773706] (0.95 confidence) |
| d2/la&#95;stack&#95;det&#95;via&#95;lu | 1.79031 ns [1.78515, 1.79569] (0.95 confidence) | 1.77958 ns [1.7735, 1.78576] (0.95 confidence) | 0.848914 ns [0.847146, 0.850911] (0.95 confidence) | 88.4573 ns [88.1649, 88.7537] (0.95 confidence) |
| d2/la&#95;stack&#95;dot | 0.602951 ns [0.598412, 0.608423] (0.95 confidence) | 0.621421 ns [0.616888, 0.626594] (0.95 confidence) | 0.605373 ns [0.603622, 0.60732] (0.95 confidence) | 2.45289 ns [2.44599, 2.46009] (0.95 confidence) |
| d2/la&#95;stack&#95;inf&#95;norm | 0.714153 ns [0.713474, 0.71494] (0.95 confidence) | 0.714906 ns [0.713251, 0.716698] (0.95 confidence) | 0.571747 ns [0.571097, 0.572463] (0.95 confidence) | 0.711529 ns [0.710137, 0.713079] (0.95 confidence) |
| d2/la&#95;stack&#95;ldlt | 7.41642 ns [7.39426, 7.43885] (0.95 confidence) | 7.01843 ns [6.98417, 7.05275] (0.95 confidence) | 1.77697 ns [1.77451, 1.77957] (0.95 confidence) | 75.3138 ns [74.9406, 75.6948] (0.95 confidence) |
| d2/la&#95;stack&#95;ldlt&#95;solve | 10.8502 ns [10.8156, 10.8863] (0.95 confidence) | 10.3573 ns [10.3067, 10.409] (0.95 confidence) | 2.81236 ns [2.80708, 2.81743] (0.95 confidence) | 115.097 ns [114.805, 115.393] (0.95 confidence) |
| d2/la&#95;stack&#95;lu | 1.60051 ns [1.5963, 1.60468] (0.95 confidence) | 1.58253 ns [1.57907, 1.58587] (0.95 confidence) | 1.66035 ns [1.65767, 1.66333] (0.95 confidence) | 84.5829 ns [84.2704, 84.9011] (0.95 confidence) |
| d2/la&#95;stack&#95;lu&#95;solve | 2.11536 ns [2.10737, 2.12377] (0.95 confidence) | 2.0934 ns [2.08389, 2.10314] (0.95 confidence) | 4.67941 ns [4.67473, 4.68429] (0.95 confidence) | 134.298 ns [133.957, 134.646] (0.95 confidence) |
| d2/la&#95;stack&#95;norm2&#95;sq | 0.397284 ns [0.396598, 0.398006] (0.95 confidence) | 0.401982 ns [0.400667, 0.403351] (0.95 confidence) | 0.398959 ns [0.397257, 0.40113] (0.95 confidence) | 4.20422 ns [4.19792, 4.21097] (0.95 confidence) |
| d2/la&#95;stack&#95;solve&#95;from&#95;ldlt | 1.26348 ns [1.2574, 1.26989] (0.95 confidence) | 1.25824 ns [1.25426, 1.262] (0.95 confidence) | 1.29493 ns [1.29232, 1.29782] (0.95 confidence) | 35.0087 ns [34.8983, 35.1253] (0.95 confidence) |
| d2/la&#95;stack&#95;solve&#95;from&#95;lu | 1.29843 ns [1.29535, 1.30161] (0.95 confidence) | 1.29294 ns [1.287, 1.29985] (0.95 confidence) | 3.00986 ns [2.99598, 3.02624] (0.95 confidence) | 42.9869 ns [42.8791, 43.0971] (0.95 confidence) |
| d3/la&#95;stack&#95;det | 1.35654 ns [1.34947, 1.36389] (0.95 confidence) | 1.36898 ns [1.36205, 1.37623] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;det&#95;from&#95;ldlt | 0.595074 ns [0.594486, 0.595847] (0.95 confidence) | 0.605429 ns [0.604267, 0.606636] (0.95 confidence) | 0.486194 ns [0.48555, 0.486894] (0.95 confidence) | 0.743291 ns [0.742724, 0.743957] (0.95 confidence) |
| d3/la&#95;stack&#95;det&#95;from&#95;lu | 0.694443 ns [0.69384, 0.695162] (0.95 confidence) | 0.704082 ns [0.702824, 0.705367] (0.95 confidence) | 0.475397 ns [0.47495, 0.475924] (0.95 confidence) | 0.993279 ns [0.991798, 0.99504] (0.95 confidence) |
| d3/la&#95;stack&#95;det&#95;via&#95;lu | 8.80459 ns [8.78623, 8.82453] (0.95 confidence) | 8.94119 ns [8.89718, 8.99881] (0.95 confidence) | 13.3016 ns [13.2505, 13.3549] (0.95 confidence) | 112.754 ns [112.344, 113.174] (0.95 confidence) |
| d3/la&#95;stack&#95;dot | 0.668047 ns [0.666654, 0.669535] (0.95 confidence) | 0.677372 ns [0.675675, 0.679158] (0.95 confidence) | 0.685014 ns [0.68365, 0.686638] (0.95 confidence) | 2.70561 ns [2.67013, 2.74218] (0.95 confidence) |
| d3/la&#95;stack&#95;inf&#95;norm | 1.36475 ns [1.36296, 1.36712] (0.95 confidence) | 1.37796 ns [1.37597, 1.37997] (0.95 confidence) | 1.17786 ns [1.17687, 1.17896] (0.95 confidence) | 1.34143 ns [1.33314, 1.35018] (0.95 confidence) |
| d3/la&#95;stack&#95;ldlt | 13.3476 ns [13.2991, 13.3968] (0.95 confidence) | 13.6916 ns [13.6368, 13.748] (0.95 confidence) | 3.94829 ns [3.94235, 3.95429] (0.95 confidence) | 84.305 ns [84.1388, 84.478] (0.95 confidence) |
| d3/la&#95;stack&#95;ldlt&#95;solve | 9.87846 ns [9.8374, 9.92295] (0.95 confidence) | 9.87467 ns [9.82948, 9.92171] (0.95 confidence) | 5.72641 ns [5.71402, 5.73905] (0.95 confidence) | 125.076 ns [124.892, 125.286] (0.95 confidence) |
| d3/la&#95;stack&#95;lu | 8.41008 ns [8.39866, 8.42146] (0.95 confidence) | 8.45659 ns [8.4397, 8.47299] (0.95 confidence) | 12.3022 ns [12.26, 12.3453] (0.95 confidence) | 107.261 ns [107.058, 107.47] (0.95 confidence) |
| d3/la&#95;stack&#95;lu&#95;solve | 9.90226 ns [9.87697, 9.93391] (0.95 confidence) | 10.1174 ns [10.0838, 10.1539] (0.95 confidence) | 32.8967 ns [32.7525, 33.0512] (0.95 confidence) | 155.301 ns [155.005, 155.611] (0.95 confidence) |
| d3/la&#95;stack&#95;norm2&#95;sq | 0.445269 ns [0.444449, 0.446226] (0.95 confidence) | 0.462042 ns [0.456977, 0.468409] (0.95 confidence) | 0.444969 ns [0.444348, 0.445567] (0.95 confidence) | 4.1323 ns [4.12855, 4.13642] (0.95 confidence) |
| d3/la&#95;stack&#95;solve&#95;from&#95;ldlt | 1.78285 ns [1.78078, 1.78517] (0.95 confidence) | 1.81893 ns [1.81132, 1.8286] (0.95 confidence) | 2.83441 ns [2.82922, 2.83921] (0.95 confidence) | 35.2794 ns [35.0231, 35.5519] (0.95 confidence) |
| d3/la&#95;stack&#95;solve&#95;from&#95;lu | 2.04759 ns [2.04383, 2.052] (0.95 confidence) | 2.09995 ns [2.09392, 2.10609] (0.95 confidence) | 4.55087 ns [4.54295, 4.55885] (0.95 confidence) | 43.652 ns [43.464, 43.8506] (0.95 confidence) |
| d32/la&#95;stack&#95;det | 2178.76 ns [2174.46, 2182.97] (0.95 confidence) | 2226.47 ns [2220.08, 2233.52] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;det&#95;from&#95;ldlt | 7.68009 ns [7.66713, 7.69229] (0.95 confidence) | 8.04439 ns [8.01937, 8.07767] (0.95 confidence) | 2.94124 ns [2.93636, 2.94609] (0.95 confidence) | 8.26594 ns [8.25378, 8.27749] (0.95 confidence) |
| d32/la&#95;stack&#95;det&#95;from&#95;lu | 9.01043 ns [8.99001, 9.03247] (0.95 confidence) | 9.34347 ns [9.30653, 9.38186] (0.95 confidence) | 3.05002 ns [3.04597, 3.05398] (0.95 confidence) | 8.65583 ns [8.63648, 8.67459] (0.95 confidence) |
| d32/la&#95;stack&#95;det&#95;via&#95;lu | 2045.67 ns [2042.74, 2048.63] (0.95 confidence) | 2068.16 ns [2064.79, 2071.52] (0.95 confidence) | 2265.99 ns [2260.38, 2271.6] (0.95 confidence) | 2320.76 ns [2317.61, 2324.09] (0.95 confidence) |
| d32/la&#95;stack&#95;dot | 3.98436 ns [3.97782, 3.99137] (0.95 confidence) | 4.08448 ns [4.07581, 4.09513] (0.95 confidence) | 3.98898 ns [3.9845, 3.9933] (0.95 confidence) | 4.84085 ns [4.83736, 4.84497] (0.95 confidence) |
| d32/la&#95;stack&#95;inf&#95;norm | 129.881 ns [129.603, 130.186] (0.95 confidence) | 131.582 ns [131.324, 131.871] (0.95 confidence) | 161.517 ns [160.399, 162.632] (0.95 confidence) | 164.558 ns [164.182, 164.934] (0.95 confidence) |
| d32/la&#95;stack&#95;ldlt | 2376.09 ns [2374.46, 2377.76] (0.95 confidence) | 2405.75 ns [2401.49, 2410.48] (0.95 confidence) | 2069.61 ns [2063.28, 2075.37] (0.95 confidence) | 1362.44 ns [1358.04, 1367.02] (0.95 confidence) |
| d32/la&#95;stack&#95;ldlt&#95;solve | 2693.12 ns [2689.48, 2696.96] (0.95 confidence) | 2754.43 ns [2745.44, 2763.88] (0.95 confidence) | 2800 ns [2786.21, 2813.19] (0.95 confidence) | 1866.72 ns [1861.86, 1871.76] (0.95 confidence) |
| d32/la&#95;stack&#95;lu | 2234.77 ns [2222.62, 2247.29] (0.95 confidence) | 2277.76 ns [2261.76, 2294.07] (0.95 confidence) | 2136.83 ns [2133.34, 2140.41] (0.95 confidence) | 2273.44 ns [2269.98, 2276.94] (0.95 confidence) |
| d32/la&#95;stack&#95;lu&#95;solve | 2740.83 ns [2735.54, 2745.6] (0.95 confidence) | 2814.24 ns [2807.44, 2822.15] (0.95 confidence) | 2625.48 ns [2601.68, 2648.63] (0.95 confidence) | 2904.46 ns [2898.17, 2911.03] (0.95 confidence) |
| d32/la&#95;stack&#95;norm2&#95;sq | 3.91906 ns [3.91386, 3.92407] (0.95 confidence) | 4.04294 ns [4.03471, 4.05352] (0.95 confidence) | 3.5769 ns [3.57217, 3.58149] (0.95 confidence) | 4.21873 ns [4.2125, 4.22578] (0.95 confidence) |
| d32/la&#95;stack&#95;solve&#95;from&#95;ldlt | 310.464 ns [309.302, 311.795] (0.95 confidence) | 316.933 ns [315.869, 318.171] (0.95 confidence) | 559.272 ns [557.325, 561.362] (0.95 confidence) | 466.93 ns [466.212, 467.641] (0.95 confidence) |
| d32/la&#95;stack&#95;solve&#95;from&#95;lu | 658.466 ns [656.148, 661.196] (0.95 confidence) | 697.653 ns [696.063, 699.269] (0.95 confidence) | 327.239 ns [326.551, 328.038] (0.95 confidence) | 627.965 ns [626.998, 628.93] (0.95 confidence) |
| d4/la&#95;stack&#95;det | 2.55124 ns [2.54593, 2.55717] (0.95 confidence) | 2.58731 ns [2.57665, 2.59814] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;det&#95;from&#95;ldlt | 0.768516 ns [0.767923, 0.769216] (0.95 confidence) | 0.77732 ns [0.776232, 0.778421] (0.95 confidence) | 0.48641 ns [0.48378, 0.489175] (0.95 confidence) | 0.993584 ns [0.992183, 0.9953] (0.95 confidence) |
| d4/la&#95;stack&#95;det&#95;from&#95;lu | 0.869856 ns [0.86857, 0.871561] (0.95 confidence) | 0.876722 ns [0.875299, 0.878225] (0.95 confidence) | 0.592761 ns [0.590054, 0.59564] (0.95 confidence) | 1.24049 ns [1.23895, 1.24221] (0.95 confidence) |
| d4/la&#95;stack&#95;det&#95;via&#95;lu | 14.2103 ns [14.1878, 14.2334] (0.95 confidence) | 14.433 ns [14.4124, 14.454] (0.95 confidence) | 32.1832 ns [32.1137, 32.2566] (0.95 confidence) | 132.871 ns [132.594, 133.172] (0.95 confidence) |
| d4/la&#95;stack&#95;dot | 0.683454 ns [0.681931, 0.685063] (0.95 confidence) | 0.688374 ns [0.68628, 0.690556] (0.95 confidence) | 0.656873 ns [0.655282, 0.658617] (0.95 confidence) | 2.49008 ns [2.48666, 2.49318] (0.95 confidence) |
| d4/la&#95;stack&#95;inf&#95;norm | 2.34543 ns [2.34265, 2.34804] (0.95 confidence) | 2.38532 ns [2.38064, 2.39056] (0.95 confidence) | 2.02266 ns [2.01992, 2.02527] (0.95 confidence) | 2.05452 ns [2.04944, 2.06015] (0.95 confidence) |
| d4/la&#95;stack&#95;ldlt | 18.3556 ns [18.2888, 18.4202] (0.95 confidence) | 18.1601 ns [18.1223, 18.2039] (0.95 confidence) | 8.34482 ns [7.94179, 8.8199] (0.95 confidence) | 105.553 ns [105.263, 105.858] (0.95 confidence) |
| d4/la&#95;stack&#95;ldlt&#95;solve | 17.8342 ns [17.8042, 17.8661] (0.95 confidence) | 18.1142 ns [17.9995, 18.2758] (0.95 confidence) | 10.6812 ns [10.6631, 10.7003] (0.95 confidence) | 147.74 ns [147.374, 148.141] (0.95 confidence) |
| d4/la&#95;stack&#95;lu | 13.7308 ns [13.7147, 13.7468] (0.95 confidence) | 13.9444 ns [13.9083, 13.9811] (0.95 confidence) | 29.6785 ns [29.6526, 29.7078] (0.95 confidence) | 127.729 ns [127.446, 128.022] (0.95 confidence) |
| d4/la&#95;stack&#95;lu&#95;solve | 21.7783 ns [21.747, 21.8085] (0.95 confidence) | 22.4324 ns [22.3586, 22.5161] (0.95 confidence) | 52.9834 ns [52.863, 53.1145] (0.95 confidence) | 180.555 ns [180.21, 180.91] (0.95 confidence) |
| d4/la&#95;stack&#95;norm2&#95;sq | 0.467722 ns [0.465841, 0.469641] (0.95 confidence) | 0.461323 ns [0.459937, 0.462879] (0.95 confidence) | 0.541316 ns [0.537781, 0.546089] (0.95 confidence) | 4.13045 ns [4.12802, 4.13298] (0.95 confidence) |
| d4/la&#95;stack&#95;solve&#95;from&#95;ldlt | 2.48106 ns [2.47432, 2.48881] (0.95 confidence) | 2.52077 ns [2.50806, 2.53411] (0.95 confidence) | 5.36961 ns [5.35828, 5.38113] (0.95 confidence) | 38.1218 ns [38.0318, 38.2184] (0.95 confidence) |
| d4/la&#95;stack&#95;solve&#95;from&#95;lu | 3.92242 ns [3.91404, 3.93106] (0.95 confidence) | 3.96815 ns [3.95874, 3.97748] (0.95 confidence) | 5.70857 ns [5.69901, 5.71793] (0.95 confidence) | 46.7236 ns [46.6336, 46.8143] (0.95 confidence) |
| d5/la&#95;stack&#95;det | 40.5626 ns [40.4776, 40.641] (0.95 confidence) | 41.0548 ns [40.976, 41.1364] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;det&#95;from&#95;ldlt | 0.875943 ns [0.875524, 0.876385] (0.95 confidence) | 0.881993 ns [0.88059, 0.883704] (0.95 confidence) | 0.577311 ns [0.575677, 0.578849] (0.95 confidence) | 1.22278 ns [1.22115, 1.22429] (0.95 confidence) |
| d5/la&#95;stack&#95;det&#95;from&#95;lu | 1.02983 ns [1.0292, 1.03043] (0.95 confidence) | 1.03599 ns [1.03485, 1.03721] (0.95 confidence) | 0.661631 ns [0.65847, 0.665089] (0.95 confidence) | 1.48064 ns [1.47929, 1.48183] (0.95 confidence) |
| d5/la&#95;stack&#95;det&#95;via&#95;lu | 33.3326 ns [33.2636, 33.4063] (0.95 confidence) | 33.5825 ns [33.4853, 33.6888] (0.95 confidence) | 50.8016 ns [50.6699, 50.932] (0.95 confidence) | 171.129 ns [170.768, 171.51] (0.95 confidence) |
| d5/la&#95;stack&#95;dot | 0.829773 ns [0.825979, 0.833375] (0.95 confidence) | 0.834933 ns [0.830421, 0.839274] (0.95 confidence) | 0.763618 ns [0.762249, 0.76501] (0.95 confidence) | 2.64238 ns [2.63714, 2.64823] (0.95 confidence) |
| d5/la&#95;stack&#95;inf&#95;norm | 3.44363 ns [3.43631, 3.45162] (0.95 confidence) | 3.65047 ns [3.58593, 3.75601] (0.95 confidence) | 3.22854 ns [3.22463, 3.23224] (0.95 confidence) | 3.15102 ns [3.14822, 3.15352] (0.95 confidence) |
| d5/la&#95;stack&#95;ldlt | 24.9443 ns [24.909, 24.9819] (0.95 confidence) | 25.4779 ns [25.4368, 25.5205] (0.95 confidence) | 23.3595 ns [22.776, 23.9525] (0.95 confidence) | 125.769 ns [125.433, 126.108] (0.95 confidence) |
| d5/la&#95;stack&#95;ldlt&#95;solve | 30.3116 ns [30.2458, 30.3809] (0.95 confidence) | 30.4785 ns [30.4198, 30.5411] (0.95 confidence) | 40.3623 ns [40.1731, 40.608] (0.95 confidence) | 184.387 ns [184.024, 184.746] (0.95 confidence) |
| d5/la&#95;stack&#95;lu | 33.516 ns [33.4495, 33.5834] (0.95 confidence) | 33.6877 ns [33.5797, 33.8048] (0.95 confidence) | 49.6134 ns [49.4432, 49.783] (0.95 confidence) | 162.691 ns [162.411, 162.977] (0.95 confidence) |
| d5/la&#95;stack&#95;lu&#95;solve | 48.6994 ns [48.5885, 48.8095] (0.95 confidence) | 48.7733 ns [48.7081, 48.8338] (0.95 confidence) | 73.3974 ns [73.1831, 73.6133] (0.95 confidence) | 267.463 ns [258.808, 282.377] (0.95 confidence) |
| d5/la&#95;stack&#95;norm2&#95;sq | 0.525 ns [0.524521, 0.525507] (0.95 confidence) | 0.525996 ns [0.525322, 0.526869] (0.95 confidence) | 0.502362 ns [0.502014, 0.502753] (0.95 confidence) | 4.11399 ns [4.11151, 4.11646] (0.95 confidence) |
| d5/la&#95;stack&#95;solve&#95;from&#95;ldlt | 3.8729 ns [3.86687, 3.87899] (0.95 confidence) | 3.93807 ns [3.92238, 3.95582] (0.95 confidence) | 8.83912 ns [8.83055, 8.84903] (0.95 confidence) | 57.2651 ns [57.1884, 57.3422] (0.95 confidence) |
| d5/la&#95;stack&#95;solve&#95;from&#95;lu | 6.08649 ns [6.0788, 6.09462] (0.95 confidence) | 6.13959 ns [6.12707, 6.1538] (0.95 confidence) | 8.83838 ns [8.82537, 8.85105] (0.95 confidence) | 78.0131 ns [77.8362, 78.1879] (0.95 confidence) |
| d64/la&#95;stack&#95;det | 15049.4 ns [14981.1, 15117.4] (0.95 confidence) | 14827.6 ns [14785.2, 14875.4] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;det&#95;from&#95;ldlt | 23.5902 ns [23.5189, 23.6648] (0.95 confidence) | 23.691 ns [23.641, 23.7465] (0.95 confidence) | 8.6143 ns [8.58258, 8.64927] (0.95 confidence) | 21.1847 ns [21.1282, 21.2413] (0.95 confidence) |
| d64/la&#95;stack&#95;det&#95;from&#95;lu | 23.8456 ns [23.7874, 23.9062] (0.95 confidence) | 23.8288 ns [23.7853, 23.8742] (0.95 confidence) | 8.86397 ns [8.83947, 8.88813] (0.95 confidence) | 22.5292 ns [22.386, 22.6794] (0.95 confidence) |
| d64/la&#95;stack&#95;det&#95;via&#95;lu | 14601.3 ns [14477.8, 14742.2] (0.95 confidence) | 14403.3 ns [14339.3, 14471.8] (0.95 confidence) | 13621.6 ns [13594.7, 13651] (0.95 confidence) | 11249.4 ns [11213.9, 11285.9] (0.95 confidence) |
| d64/la&#95;stack&#95;dot | 10.9747 ns [10.952, 10.9974] (0.95 confidence) | 11.2887 ns [11.2595, 11.3201] (0.95 confidence) | 8.0296 ns [8.00277, 8.05639] (0.95 confidence) | 8.09611 ns [8.08582, 8.10758] (0.95 confidence) |
| d64/la&#95;stack&#95;inf&#95;norm | 663.278 ns [662.471, 664.147] (0.95 confidence) | 644.586 ns [643.187, 645.984] (0.95 confidence) | 1342.8 ns [1339.89, 1346.03] (0.95 confidence) | 1606.65 ns [1600.76, 1613.09] (0.95 confidence) |
| d64/la&#95;stack&#95;ldlt | 17230.5 ns [17193.5, 17269.3] (0.95 confidence) | 16883.2 ns [16862.5, 16906.3] (0.95 confidence) | 10914.4 ns [10840, 10994.7] (0.95 confidence) | 7202.42 ns [7172.92, 7234.92] (0.95 confidence) |
| d64/la&#95;stack&#95;ldlt&#95;solve | 19119.2 ns [19085.8, 19153.4] (0.95 confidence) | 19323 ns [19282.7, 19364.8] (0.95 confidence) | 12248.7 ns [12232.9, 12265.2] (0.95 confidence) | 8376.88 ns [8362.22, 8393.55] (0.95 confidence) |
| d64/la&#95;stack&#95;lu | 15530.4 ns [15479.9, 15582.8] (0.95 confidence) | 15224.6 ns [15190.8, 15261.1] (0.95 confidence) | 12839.5 ns [12787.8, 12894.7] (0.95 confidence) | 11163 ns [11141.4, 11185.6] (0.95 confidence) |
| d64/la&#95;stack&#95;lu&#95;solve | 17440.1 ns [17346.3, 17556.5] (0.95 confidence) | 17342.3 ns [17295.7, 17393.7] (0.95 confidence) | 14468.6 ns [14438.2, 14501.7] (0.95 confidence) | 12986.9 ns [12949.9, 13023.7] (0.95 confidence) |
| d64/la&#95;stack&#95;norm2&#95;sq | 11.4621 ns [11.3945, 11.5309] (0.95 confidence) | 11.0976 ns [11.0725, 11.1214] (0.95 confidence) | 8.94431 ns [8.90637, 8.98568] (0.95 confidence) | 6.40725 ns [6.36537, 6.45128] (0.95 confidence) |
| d64/la&#95;stack&#95;solve&#95;from&#95;ldlt | 1121.49 ns [1117.84, 1125.87] (0.95 confidence) | 1100.18 ns [1097.91, 1102.58] (0.95 confidence) | 1325 ns [1321.67, 1328.81] (0.95 confidence) | 1276.12 ns [1269.81, 1283.06] (0.95 confidence) |
| d64/la&#95;stack&#95;solve&#95;from&#95;lu | 2664.88 ns [2649.22, 2680.47] (0.95 confidence) | 2618.03 ns [2610.79, 2625.88] (0.95 confidence) | 815.944 ns [813.85, 818.129] (0.95 confidence) | 1707.09 ns [1703.71, 1710.59] (0.95 confidence) |
| d8/la&#95;stack&#95;det | 93.4319 ns [93.0029, 93.8484] (0.95 confidence) | 96.7167 ns [96.0907, 97.3599] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;det&#95;from&#95;ldlt | 1.32003 ns [1.31814, 1.32212] (0.95 confidence) | 1.32934 ns [1.32635, 1.3324] (0.95 confidence) | 0.94603 ns [0.928354, 0.969663] (0.95 confidence) | 2.04618 ns [2.02106, 2.0763] (0.95 confidence) |
| d8/la&#95;stack&#95;det&#95;from&#95;ldlt&#95;balanced&#95;range | 8.91368 ns [8.89804, 8.92844] (0.95 confidence) | 8.93074 ns [8.91084, 8.9503] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;det&#95;from&#95;lu | 1.48156 ns [1.47986, 1.48345] (0.95 confidence) | 1.50098 ns [1.49806, 1.50396] (0.95 confidence) | 0.998838 ns [0.994562, 1.00287] (0.95 confidence) | 2.33274 ns [2.2976, 2.37236] (0.95 confidence) |
| d8/la&#95;stack&#95;det&#95;from&#95;lu&#95;balanced&#95;range | 9.23268 ns [9.21257, 9.25348] (0.95 confidence) | 9.29322 ns [9.2702, 9.3174] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;det&#95;via&#95;lu | 90.3591 ns [89.7671, 90.9393] (0.95 confidence) | 88.5765 ns [87.9204, 89.2792] (0.95 confidence) | 137.682 ns [137.457, 137.905] (0.95 confidence) | 243.108 ns [242.667, 243.552] (0.95 confidence) |
| d8/la&#95;stack&#95;dot | 0.941413 ns [0.939238, 0.94378] (0.95 confidence) | 0.928519 ns [0.92603, 0.931247] (0.95 confidence) | 1.38149 ns [1.37268, 1.39128] (0.95 confidence) | 2.67403 ns [2.65094, 2.69868] (0.95 confidence) |
| d8/la&#95;stack&#95;inf&#95;norm | 8.48411 ns [8.46971, 8.49943] (0.95 confidence) | 8.56628 ns [8.55053, 8.58604] (0.95 confidence) | 8.086 ns [8.07027, 8.10178] (0.95 confidence) | 8.14067 ns [8.12845, 8.15271] (0.95 confidence) |
| d8/la&#95;stack&#95;ldlt | 91.717 ns [91.1439, 92.5411] (0.95 confidence) | 91.5053 ns [91.3342, 91.6824] (0.95 confidence) | 104.495 ns [103.83, 105.268] (0.95 confidence) | 186.414 ns [185.972, 186.858] (0.95 confidence) |
| d8/la&#95;stack&#95;ldlt&#95;ill&#95;conditioned | 90.7398 ns [90.6416, 90.8375] (0.95 confidence) | 91.2669 ns [91.1005, 91.4395] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;ldlt&#95;solve | 104.31 ns [103.961, 104.69] (0.95 confidence) | 104.888 ns [104.544, 105.274] (0.95 confidence) | 147.287 ns [145.988, 148.603] (0.95 confidence) | 263.89 ns [262.768, 265.108] (0.95 confidence) |
| d8/la&#95;stack&#95;lu | 88.9922 ns [88.5584, 89.4008] (0.95 confidence) | 90.6949 ns [90.2205, 91.1759] (0.95 confidence) | 133.305 ns [132.852, 133.702] (0.95 confidence) | 231.746 ns [230.781, 232.785] (0.95 confidence) |
| d8/la&#95;stack&#95;lu&#95;ill&#95;conditioned | 89.3846 ns [88.9696, 89.7872] (0.95 confidence) | 90.0925 ns [89.6298, 90.5498] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;lu&#95;pivoting | 96.7127 ns [96.3621, 97.0772] (0.95 confidence) | 95.6654 ns [95.4061, 95.9271] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;lu&#95;solve | 130.625 ns [130.283, 131.008] (0.95 confidence) | 134.478 ns [133.489, 135.995] (0.95 confidence) | 165.573 ns [165.151, 166.045] (0.95 confidence) | 344.186 ns [342.441, 346.1] (0.95 confidence) |
| d8/la&#95;stack&#95;norm2&#95;sq | 0.657056 ns [0.655199, 0.658876] (0.95 confidence) | 0.65978 ns [0.657268, 0.662373] (0.95 confidence) | 0.941667 ns [0.940616, 0.942819] (0.95 confidence) | 4.17722 ns [4.17268, 4.18212] (0.95 confidence) |
| d8/la&#95;stack&#95;solve&#95;from&#95;ldlt | 8.23761 ns [8.22317, 8.25198] (0.95 confidence) | 8.22696 ns [8.20965, 8.24439] (0.95 confidence) | 21.9222 ns [21.8563, 21.997] (0.95 confidence) | 68.6883 ns [68.5996, 68.7819] (0.95 confidence) |
| d8/la&#95;stack&#95;solve&#95;from&#95;lu | 13.6434 ns [13.5957, 13.7003] (0.95 confidence) | 13.5593 ns [13.5302, 13.5893] (0.95 confidence) | 16.3146 ns [16.2753, 16.3525] (0.95 confidence) | 96.5555 ns [96.3241, 96.7952] (0.95 confidence) |
| exact&#95;d2/det | 0.464907 ns [0.45892, 0.471229] (0.95 confidence) | 0.438622 ns [0.437109, 0.44045] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;direct | 0.461442 ns [0.456221, 0.466904] (0.95 confidence) | 0.436723 ns [0.435053, 0.439191] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;direct&#95;with&#95;errbound | 1.76793 ns [1.7644, 1.77151] (0.95 confidence) | 1.74735 ns [1.7449, 1.75002] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;errbound | 1.61701 ns [1.61379, 1.6208] (0.95 confidence) | 1.60287 ns [1.59941, 1.60668] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;exact | 86.0365 ns [85.7999, 86.2716] (0.95 confidence) | 88.6251 ns [88.4409, 88.8313] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;exact&#95;f64&#95;result | 61.0576 ns [60.8153, 61.3381] (0.95 confidence) | 65.7729 ns [65.5712, 66.0239] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;exact&#95;rounded&#95;f64 | 62.2699 ns [62.1371, 62.4026] (0.95 confidence) | 68.1264 ns [67.9205, 68.3313] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;sign&#95;exact | 3.03592 ns [3.03191, 3.04033] (0.95 confidence) | 2.99344 ns [2.99019, 2.99695] (0.95 confidence) | — | — |
| exact&#95;d2/solve&#95;exact | 7537.1 ns [7475.86, 7598.94] (0.95 confidence) | 3018.79 ns [3013.79, 3024.01] (0.95 confidence) | — | — |
| exact&#95;d2/solve&#95;exact&#95;f64&#95;result | 7292.51 ns [7283.1, 7302.31] (0.95 confidence) | 3162.04 ns [3156.26, 3167.84] (0.95 confidence) | — | — |
| exact&#95;d2/solve&#95;exact&#95;rounded&#95;f64 | 7525.39 ns [7507.21, 7545.57] (0.95 confidence) | 3340.97 ns [3335.36, 3346.73] (0.95 confidence) | — | — |
| exact&#95;d3/det | 0.845755 ns [0.843419, 0.848107] (0.95 confidence) | 0.829088 ns [0.82724, 0.831074] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;direct | 0.833231 ns [0.831605, 0.835012] (0.95 confidence) | 0.835991 ns [0.833528, 0.838534] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;direct&#95;with&#95;errbound | 3.36017 ns [3.35539, 3.36518] (0.95 confidence) | 3.34331 ns [3.33932, 3.34744] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;errbound | 3.34429 ns [3.34085, 3.3481] (0.95 confidence) | 3.34289 ns [3.33442, 3.35338] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;exact | 292.325 ns [291.789, 292.926] (0.95 confidence) | 302.6 ns [302.175, 303.049] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;exact&#95;f64&#95;result | 266.492 ns [265.909, 267.122] (0.95 confidence) | 278.371 ns [277.979, 278.784] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;exact&#95;rounded&#95;f64 | 271.485 ns [270.855, 272.19] (0.95 confidence) | 282.931 ns [282.19, 283.718] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;sign&#95;exact | 4.79567 ns [4.7893, 4.80252] (0.95 confidence) | 4.75107 ns [4.74541, 4.75719] (0.95 confidence) | — | — |
| exact&#95;d3/solve&#95;exact | 30383.4 ns [30339.2, 30430.7] (0.95 confidence) | 9174.5 ns [9142.84, 9205.56] (0.95 confidence) | — | — |
| exact&#95;d3/solve&#95;exact&#95;f64&#95;result | 30435.4 ns [30358.5, 30511.4] (0.95 confidence) | 9291.59 ns [9266.41, 9317.15] (0.95 confidence) | — | — |
| exact&#95;d3/solve&#95;exact&#95;rounded&#95;f64 | 30825.5 ns [30771.2, 30882.2] (0.95 confidence) | 9621.4 ns [9581.51, 9666.6] (0.95 confidence) | — | — |
| exact&#95;d4/det | 3.28635 ns [3.28048, 3.29261] (0.95 confidence) | 3.23493 ns [3.23169, 3.23836] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;direct | 2.44448 ns [2.4335, 2.45702] (0.95 confidence) | 2.3744 ns [2.37144, 2.37765] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;direct&#95;with&#95;errbound | 6.86216 ns [6.84796, 6.87715] (0.95 confidence) | 6.8604 ns [6.84425, 6.87811] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;errbound | 6.80357 ns [6.79007, 6.81857] (0.95 confidence) | 6.82216 ns [6.81282, 6.83201] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;exact | 820.307 ns [815.294, 825.391] (0.95 confidence) | 852.843 ns [850.565, 855.164] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;exact&#95;f64&#95;result | 804.866 ns [798.701, 811.207] (0.95 confidence) | 811.722 ns [809.784, 813.748] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;exact&#95;rounded&#95;f64 | 780.739 ns [777.45, 784.306] (0.95 confidence) | 849.567 ns [840.736, 863.589] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;sign&#95;exact | 8.01255 ns [8.00282, 8.02293] (0.95 confidence) | 8.09025 ns [8.06866, 8.11287] (0.95 confidence) | — | — |
| exact&#95;d4/solve&#95;exact | 79076.2 ns [78837.3, 79318.4] (0.95 confidence) | 18444.1 ns [18406.8, 18481.3] (0.95 confidence) | — | — |
| exact&#95;d4/solve&#95;exact&#95;f64&#95;result | 79069.8 ns [78920.9, 79220.3] (0.95 confidence) | 18609.3 ns [18574.3, 18644.9] (0.95 confidence) | — | — |
| exact&#95;d4/solve&#95;exact&#95;rounded&#95;f64 | 79528.3 ns [79388.4, 79680.7] (0.95 confidence) | 19199.3 ns [19175.1, 19224.1] (0.95 confidence) | — | — |
| exact&#95;d5/det | 28.0197 ns [27.9559, 28.0857] (0.95 confidence) | 27.541 ns [27.5043, 27.5801] (0.95 confidence) | — | — |
| exact&#95;d5/det&#95;exact | 2945.78 ns [2940.09, 2951.63] (0.95 confidence) | 2877.53 ns [2871.03, 2884.36] (0.95 confidence) | — | — |
| exact&#95;d5/det&#95;exact&#95;f64&#95;result | 2902.11 ns [2895.72, 2909.22] (0.95 confidence) | 2838.26 ns [2832.07, 2844.78] (0.95 confidence) | — | — |
| exact&#95;d5/det&#95;exact&#95;rounded&#95;f64 | 2930.32 ns [2924.26, 2936.49] (0.95 confidence) | 2894.38 ns [2885.57, 2903.09] (0.95 confidence) | — | — |
| exact&#95;d5/det&#95;sign&#95;exact | 2907.28 ns [2899.71, 2915.84] (0.95 confidence) | 2852.69 ns [2847.17, 2858.4] (0.95 confidence) | — | — |
| exact&#95;d5/solve&#95;exact | 155228 ns [154979, 155486] (0.95 confidence) | 29946.2 ns [29871.4, 30024.3] (0.95 confidence) | — | — |
| exact&#95;d5/solve&#95;exact&#95;f64&#95;result | 155448 ns [155242, 155666] (0.95 confidence) | 30219.4 ns [30119.8, 30321.7] (0.95 confidence) | — | — |
| exact&#95;d5/solve&#95;exact&#95;rounded&#95;f64 | 156528 ns [156244, 156839] (0.95 confidence) | 31044.4 ns [30949, 31144.4] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/det&#95;exact | 891.058 ns [887.612, 894.765] (0.95 confidence) | 920.824 ns [918.883, 922.932] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/det&#95;sign&#95;exact | 8.03062 ns [8.01891, 8.0425] (0.95 confidence) | 7.94933 ns [7.93841, 7.96289] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/solve&#95;exact | 58847.5 ns [58697.6, 58999.6] (0.95 confidence) | 15244.7 ns [15183.9, 15306.6] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/solve&#95;exact&#95;f64&#95;result | 58994.6 ns [58835.1, 59212.5] (0.95 confidence) | 15746.7 ns [15669.6, 15829] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/solve&#95;exact&#95;rounded&#95;f64 | 59434.9 ns [59277.4, 59605] (0.95 confidence) | 15822.8 ns [15790.9, 15855.2] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/det&#95;exact | 3112.34 ns [3103.53, 3121.59] (0.95 confidence) | 3005.41 ns [2999.04, 3011.71] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/det&#95;sign&#95;exact | 3084.84 ns [3078.17, 3091.58] (0.95 confidence) | 2920.05 ns [2914.93, 2925.26] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/solve&#95;exact | 123307 ns [122931, 123735] (0.95 confidence) | 24944.5 ns [24871.9, 25020.4] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/solve&#95;exact&#95;f64&#95;result | 122084 ns [121814, 122360] (0.95 confidence) | 25183 ns [25127.4, 25239.9] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/solve&#95;exact&#95;rounded&#95;f64 | 123499 ns [123173, 123840] (0.95 confidence) | 25915.2 ns [25870.7, 25960.7] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/det&#95;exact | 311.839 ns [310.664, 313.034] (0.95 confidence) | 323.307 ns [322.517, 324.173] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/det&#95;sign&#95;exact | 287.356 ns [286.531, 288.211] (0.95 confidence) | 299.469 ns [298.25, 300.804] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/solve&#95;exact | 95674.4 ns [95431.2, 95932.8] (0.95 confidence) | 9124.89 ns [9098.44, 9155.92] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/solve&#95;exact&#95;f64&#95;result | 95695.7 ns [95423.5, 96056.2] (0.95 confidence) | 9180.66 ns [9167.34, 9197.66] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/solve&#95;exact&#95;rounded&#95;f64 | 95947.4 ns [95764.7, 96139.2] (0.95 confidence) | 9550.61 ns [9536.34, 9565.4] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/det&#95;exact | 297.329 ns [296.509, 298.201] (0.95 confidence) | 185.434 ns [185.047, 185.863] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/det&#95;sign&#95;exact | 323.245 ns [322.205, 324.332] (0.95 confidence) | 168.077 ns [167.735, 168.434] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/solve&#95;exact | 2476.26 ns [2470.38, 2482.22] (0.95 confidence) | 988.409 ns [985.691, 991.339] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/solve&#95;exact&#95;f64&#95;result | 2512.44 ns [2502.16, 2523.84] (0.95 confidence) | 1013.99 ns [1012.38, 1015.63] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/solve&#95;exact&#95;rounded&#95;f64 | 2516.79 ns [2511.2, 2522.2] (0.95 confidence) | 1015.92 ns [1012.72, 1019.25] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/det&#95;exact | 3105 ns [3100.9, 3109.17] (0.95 confidence) | 3464.77 ns [3457.91, 3471.56] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/det&#95;sign&#95;exact | 148.313 ns [147.831, 148.806] (0.95 confidence) | 144.304 ns [143.768, 144.909] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/solve&#95;exact | 75398.3 ns [75234.5, 75567] (0.95 confidence) | 35843.4 ns [35739.9, 35954.7] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/solve&#95;exact&#95;f64&#95;result | 76281.8 ns [76124.8, 76440.9] (0.95 confidence) | 36496.1 ns [36430.2, 36566.2] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/solve&#95;exact&#95;rounded&#95;f64 | 76464.5 ns [76303.4, 76628.3] (0.95 confidence) | 36770.7 ns [36681.3, 36872.1] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/det&#95;exact | 8942.22 ns [8931.27, 8953.46] (0.95 confidence) | 9414.13 ns [9403.91, 9424.97] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/det&#95;sign&#95;exact | 228.254 ns [227.729, 228.888] (0.95 confidence) | 225.823 ns [225.61, 226.056] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/solve&#95;exact | 235284 ns [235052, 235542] (0.95 confidence) | 81151.3 ns [81042.1, 81260.5] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/solve&#95;exact&#95;f64&#95;result | 236745 ns [236408, 237100] (0.95 confidence) | 81956.8 ns [81834.2, 82094.9] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/solve&#95;exact&#95;rounded&#95;f64 | 236941 ns [236667, 237238] (0.95 confidence) | 81960.2 ns [81866.9, 82059.1] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/det&#95;exact | 22390.6 ns [22303.4, 22487.1] (0.95 confidence) | 23197.9 ns [23134.6, 23262.5] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/det&#95;sign&#95;exact | 459.599 ns [458.585, 460.815] (0.95 confidence) | 444.786 ns [444.035, 445.643] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/solve&#95;exact | 559687 ns [558956, 560522] (0.95 confidence) | 144803 ns [144613, 145004] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/solve&#95;exact&#95;f64&#95;result | 560387 ns [559863, 560933] (0.95 confidence) | 145684 ns [145515, 145862] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/solve&#95;exact&#95;rounded&#95;f64 | 557999 ns [557228, 558830] (0.95 confidence) | 146546 ns [146402, 146693] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/det&#95;exact | 56764.7 ns [56506.1, 57052.2] (0.95 confidence) | 55924.2 ns [55843.4, 56009.9] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/det&#95;sign&#95;exact | 56003 ns [55800, 56225.7] (0.95 confidence) | 54997.7 ns [54895.1, 55106] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/solve&#95;exact | 1.11323e+06 ns [1.11161e+06, 1.11501e+06] (0.95 confidence) | 239488 ns [238616, 240500] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/solve&#95;exact&#95;f64&#95;result | 1.1119e+06 ns [1.11028e+06, 1.11349e+06] (0.95 confidence) | 238779 ns [238398, 239195] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/solve&#95;exact&#95;rounded&#95;f64 | 1.11509e+06 ns [1.11124e+06, 1.11931e+06] (0.95 confidence) | 240106 ns [239644, 240636] (0.95 confidence) | — | — |
| rational&#95;input&#95;d2/det&#95;big&#95;rational&#95;gaussian | 1133 ns [1128.22, 1138.26] (0.95 confidence) | 1113.9 ns [1112.29, 1115.58] (0.95 confidence) | — | — |
| rational&#95;input&#95;d2/det&#95;row&#95;cleared&#95;bareiss | 180.27 ns [179.87, 180.69] (0.95 confidence) | 184.677 ns [182.891, 186.533] (0.95 confidence) | — | — |
| rational&#95;input&#95;d2/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 77.4542 ns [77.3284, 77.5811] (0.95 confidence) | 77.5408 ns [77.3608, 77.7628] (0.95 confidence) | — | — |
| rational&#95;input&#95;d2/solve&#95;big&#95;rational&#95;gaussian | 2156.08 ns [2149.25, 2164.12] (0.95 confidence) | 2134.63 ns [2129.79, 2140.59] (0.95 confidence) | — | — |
| rational&#95;input&#95;d2/solve&#95;row&#95;cleared&#95;bareiss | 936.195 ns [934.232, 938.594] (0.95 confidence) | 612.613 ns [611.602, 613.711] (0.95 confidence) | — | — |
| rational&#95;input&#95;d3/det&#95;big&#95;rational&#95;gaussian | 4773.7 ns [4762.67, 4785.17] (0.95 confidence) | 4757.51 ns [4743.32, 4772.24] (0.95 confidence) | — | — |
| rational&#95;input&#95;d3/det&#95;row&#95;cleared&#95;bareiss | 570.351 ns [568.925, 571.835] (0.95 confidence) | 571.875 ns [570.466, 573.34] (0.95 confidence) | — | — |
| rational&#95;input&#95;d3/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 361.694 ns [360.876, 362.579] (0.95 confidence) | 371.405 ns [368.669, 374.253] (0.95 confidence) | — | — |
| rational&#95;input&#95;d3/solve&#95;big&#95;rational&#95;gaussian | 8252.29 ns [8230.16, 8274.62] (0.95 confidence) | 8156.96 ns [8139.08, 8175.63] (0.95 confidence) | — | — |
| rational&#95;input&#95;d3/solve&#95;row&#95;cleared&#95;bareiss | 2504.33 ns [2479.3, 2538.74] (0.95 confidence) | 1245.75 ns [1243.53, 1248.3] (0.95 confidence) | — | — |
| rational&#95;input&#95;d4/det&#95;big&#95;rational&#95;gaussian | 14030.8 ns [13993.7, 14067.5] (0.95 confidence) | 13763.4 ns [13735.3, 13791.8] (0.95 confidence) | — | — |
| rational&#95;input&#95;d4/det&#95;row&#95;cleared&#95;bareiss | 1174.21 ns [1171.29, 1177.42] (0.95 confidence) | 1169.83 ns [1166.78, 1173.78] (0.95 confidence) | — | — |
| rational&#95;input&#95;d4/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 776.357 ns [774.662, 778.148] (0.95 confidence) | 781.737 ns [779.595, 783.941] (0.95 confidence) | — | — |
| rational&#95;input&#95;d4/solve&#95;big&#95;rational&#95;gaussian | 24157.7 ns [24100.8, 24216.8] (0.95 confidence) | 24446.8 ns [24139.4, 24971.8] (0.95 confidence) | — | — |
| rational&#95;input&#95;d4/solve&#95;row&#95;cleared&#95;bareiss | 6282.84 ns [6181.89, 6385.94] (0.95 confidence) | 2191.78 ns [2189.08, 2194.59] (0.95 confidence) | — | — |
| rational&#95;input&#95;d5/det&#95;big&#95;rational&#95;gaussian | 23789.7 ns [23703.9, 23876.5] (0.95 confidence) | 23922.4 ns [23867.9, 23978.9] (0.95 confidence) | — | — |
| rational&#95;input&#95;d5/det&#95;row&#95;cleared&#95;bareiss | 2058.93 ns [2055.11, 2063.08] (0.95 confidence) | 2079.25 ns [2075.91, 2082.75] (0.95 confidence) | — | — |
| rational&#95;input&#95;d5/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 1554.91 ns [1551.3, 1558.57] (0.95 confidence) | 1566.96 ns [1563.79, 1570.23] (0.95 confidence) | — | — |
| rational&#95;input&#95;d5/solve&#95;big&#95;rational&#95;gaussian | 40182.1 ns [40060.5, 40304.3] (0.95 confidence) | 39652.6 ns [39591.6, 39716.5] (0.95 confidence) | — | — |
| rational&#95;input&#95;d5/solve&#95;row&#95;cleared&#95;bareiss | 9578.22 ns [9557.99, 9599.99] (0.95 confidence) | 3391.9 ns [3384.43, 3400.11] (0.95 confidence) | — | — |
| rational&#95;input&#95;d6/det&#95;big&#95;rational&#95;gaussian | 56861.2 ns [56750.8, 56973.9] (0.95 confidence) | 58123.1 ns [57997.4, 58265.5] (0.95 confidence) | — | — |
| rational&#95;input&#95;d6/det&#95;row&#95;cleared&#95;bareiss | 3543.09 ns [3537.04, 3549.27] (0.95 confidence) | 3529.54 ns [3521.29, 3538.61] (0.95 confidence) | — | — |
| rational&#95;input&#95;d6/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 2859.34 ns [2853.28, 2865.5] (0.95 confidence) | 2915.53 ns [2900.21, 2935.73] (0.95 confidence) | — | — |
| rational&#95;input&#95;d6/solve&#95;big&#95;rational&#95;gaussian | 91090.5 ns [90854.2, 91346.7] (0.95 confidence) | 96016.6 ns [95344.3, 96732.7] (0.95 confidence) | — | — |
| rational&#95;input&#95;d6/solve&#95;row&#95;cleared&#95;bareiss | 17872.1 ns [17839.8, 17906.8] (0.95 confidence) | 5774.91 ns [5762.19, 5788.64] (0.95 confidence) | — | — |
| rational&#95;input&#95;d7/det&#95;big&#95;rational&#95;gaussian | 110964 ns [110577, 111406] (0.95 confidence) | 112441 ns [112208, 112697] (0.95 confidence) | — | — |
| rational&#95;input&#95;d7/det&#95;row&#95;cleared&#95;bareiss | 5614.79 ns [5605.19, 5624.57] (0.95 confidence) | 5595.66 ns [5586.14, 5605.54] (0.95 confidence) | — | — |
| rational&#95;input&#95;d7/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 4666.84 ns [4660.3, 4673.55] (0.95 confidence) | 4760.44 ns [4752.12, 4768.72] (0.95 confidence) | — | — |
| rational&#95;input&#95;d7/solve&#95;big&#95;rational&#95;gaussian | 162196 ns [161885, 162525] (0.95 confidence) | 169179 ns [168176, 170257] (0.95 confidence) | — | — |
| rational&#95;input&#95;d7/solve&#95;row&#95;cleared&#95;bareiss | 27619.4 ns [27571.4, 27666.7] (0.95 confidence) | 10132.3 ns [10112.8, 10153.3] (0.95 confidence) | — | — |
| rational&#95;input&#95;d8/det&#95;big&#95;rational&#95;gaussian | 188248 ns [187806, 188689] (0.95 confidence) | 194331 ns [193766, 194947] (0.95 confidence) | — | — |
| rational&#95;input&#95;d8/det&#95;row&#95;cleared&#95;bareiss | 8097.63 ns [8078.87, 8117.44] (0.95 confidence) | 8237.34 ns [8219.85, 8256.1] (0.95 confidence) | — | — |
| rational&#95;input&#95;d8/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 7241.89 ns [7229.39, 7254.24] (0.95 confidence) | 7386.17 ns [7369.49, 7403.23] (0.95 confidence) | — | — |
| rational&#95;input&#95;d8/solve&#95;big&#95;rational&#95;gaussian | 268381 ns [267729, 269058] (0.95 confidence) | 272901 ns [272515, 273317] (0.95 confidence) | — | — |
| rational&#95;input&#95;d8/solve&#95;row&#95;cleared&#95;bareiss | 40862.4 ns [40797.8, 40931] (0.95 confidence) | 15394.7 ns [15364.1, 15426.3] (0.95 confidence) | — | — |

## Median (ns)

| Case | la-stack baseline | la-stack current | nalgebra | faer |
| --- | --- | --- | --- | --- |
| d16/la&#95;stack&#95;det | 431.355 ns [429.554, 433.409] (0.95 confidence) | 435.958 ns [434.147, 439.205] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;det&#95;from&#95;ldlt | 2.89883 ns [2.88973, 2.90569] (0.95 confidence) | 2.91858 ns [2.9156, 2.92263] (0.95 confidence) | 1.68913 ns [1.68305, 1.69439] (0.95 confidence) | 4.36104 ns [4.35315, 4.3734] (0.95 confidence) |
| d16/la&#95;stack&#95;det&#95;from&#95;lu | 3.29569 ns [3.29208, 3.29911] (0.95 confidence) | 3.30504 ns [3.30178, 3.31021] (0.95 confidence) | 1.77695 ns [1.76877, 1.78267] (0.95 confidence) | 4.6818 ns [4.65863, 4.76215] (0.95 confidence) |
| d16/la&#95;stack&#95;det&#95;via&#95;lu | 414.441 ns [409.712, 417.837] (0.95 confidence) | 419.372 ns [416.471, 424.374] (0.95 confidence) | 456.7 ns [455.176, 459.057] (0.95 confidence) | 637.56 ns [633.11, 642.581] (0.95 confidence) |
| d16/la&#95;stack&#95;dot | 2.34307 ns [2.33946, 2.34768] (0.95 confidence) | 2.34721 ns [2.3409, 2.35083] (0.95 confidence) | 1.84897 ns [1.84649, 1.85204] (0.95 confidence) | 3.25757 ns [3.25467, 3.2636] (0.95 confidence) |
| d16/la&#95;stack&#95;inf&#95;norm | 32.0776 ns [32.0146, 32.1088] (0.95 confidence) | 32.8858 ns [32.8579, 32.9355] (0.95 confidence) | 32.2176 ns [32.1884, 32.2615] (0.95 confidence) | 32.1229 ns [32.1084, 32.1492] (0.95 confidence) |
| d16/la&#95;stack&#95;ldlt | 340.841 ns [340.049, 341.613] (0.95 confidence) | 341.567 ns [340.675, 342.421] (0.95 confidence) | 417.765 ns [414.879, 418.891] (0.95 confidence) | 404.975 ns [403.504, 406.537] (0.95 confidence) |
| d16/la&#95;stack&#95;ldlt&#95;solve | 383.51 ns [383.128, 383.934] (0.95 confidence) | 384.304 ns [383.895, 384.793] (0.95 confidence) | 637.749 ns [636.778, 639.595] (0.95 confidence) | 602.783 ns [601.444, 605.333] (0.95 confidence) |
| d16/la&#95;stack&#95;lu | 407.869 ns [405.9, 411.207] (0.95 confidence) | 397.951 ns [394.518, 399.731] (0.95 confidence) | 462.218 ns [460.297, 463.499] (0.95 confidence) | 614.536 ns [610.04, 618.949] (0.95 confidence) |
| d16/la&#95;stack&#95;lu&#95;solve | 657.049 ns [654.101, 659.92] (0.95 confidence) | 664.211 ns [662.337, 667.556] (0.95 confidence) | 592.628 ns [590.988, 594.368] (0.95 confidence) | 871.515 ns [869.043, 873.723] (0.95 confidence) |
| d16/la&#95;stack&#95;norm2&#95;sq | 2.10975 ns [2.10679, 2.11445] (0.95 confidence) | 2.11572 ns [2.11167, 2.12005] (0.95 confidence) | 1.72908 ns [1.7266, 1.73194] (0.95 confidence) | 4.10049 ns [4.09225, 4.10847] (0.95 confidence) |
| d16/la&#95;stack&#95;solve&#95;from&#95;ldlt | 27.8216 ns [27.7944, 27.8683] (0.95 confidence) | 27.934 ns [27.9019, 27.9761] (0.95 confidence) | 118.985 ns [118.77, 119.395] (0.95 confidence) | 179.969 ns [179.588, 180.171] (0.95 confidence) |
| d16/la&#95;stack&#95;solve&#95;from&#95;lu | 190.734 ns [190.126, 191.305] (0.95 confidence) | 200.99 ns [200.711, 201.414] (0.95 confidence) | 93.6653 ns [93.4961, 93.806] (0.95 confidence) | 247.034 ns [246.643, 247.187] (0.95 confidence) |
| d2/la&#95;stack&#95;det | 0.589071 ns [0.58771, 0.590455] (0.95 confidence) | 0.580101 ns [0.579276, 0.581777] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;det&#95;from&#95;ldlt | 0.504966 ns [0.504045, 0.505835] (0.95 confidence) | 0.506699 ns [0.504713, 0.508559] (0.95 confidence) | 0.44602 ns [0.445693, 0.446374] (0.95 confidence) | 0.578853 ns [0.578537, 0.579115] (0.95 confidence) |
| d2/la&#95;stack&#95;det&#95;from&#95;lu | 0.560109 ns [0.559738, 0.560774] (0.95 confidence) | 0.553998 ns [0.553136, 0.554527] (0.95 confidence) | 0.467319 ns [0.465846, 0.46907] (0.95 confidence) | 0.771338 ns [0.769893, 0.77228] (0.95 confidence) |
| d2/la&#95;stack&#95;det&#95;via&#95;lu | 1.78382 ns [1.77935, 1.79074] (0.95 confidence) | 1.78332 ns [1.77429, 1.78663] (0.95 confidence) | 0.846381 ns [0.845059, 0.848121] (0.95 confidence) | 88.3003 ns [88.0676, 88.6565] (0.95 confidence) |
| d2/la&#95;stack&#95;dot | 0.594018 ns [0.592489, 0.59575] (0.95 confidence) | 0.621878 ns [0.620068, 0.623559] (0.95 confidence) | 0.601261 ns [0.600803, 0.60237] (0.95 confidence) | 2.43825 ns [2.43058, 2.45104] (0.95 confidence) |
| d2/la&#95;stack&#95;inf&#95;norm | 0.712774 ns [0.712264, 0.713288] (0.95 confidence) | 0.712689 ns [0.711605, 0.715323] (0.95 confidence) | 0.570926 ns [0.570514, 0.571607] (0.95 confidence) | 0.708844 ns [0.70798, 0.709461] (0.95 confidence) |
| d2/la&#95;stack&#95;ldlt | 7.40784 ns [7.39823, 7.41923] (0.95 confidence) | 7.01983 ns [6.97531, 7.04159] (0.95 confidence) | 1.77699 ns [1.77357, 1.77972] (0.95 confidence) | 75.2039 ns [74.6867, 75.6538] (0.95 confidence) |
| d2/la&#95;stack&#95;ldlt&#95;solve | 10.8182 ns [10.792, 10.8503] (0.95 confidence) | 10.2846 ns [10.2338, 10.4365] (0.95 confidence) | 2.81043 ns [2.80574, 2.81848] (0.95 confidence) | 114.949 ns [114.708, 115.549] (0.95 confidence) |
| d2/la&#95;stack&#95;lu | 1.60072 ns [1.59544, 1.60676] (0.95 confidence) | 1.58602 ns [1.58198, 1.58856] (0.95 confidence) | 1.65908 ns [1.65715, 1.66273] (0.95 confidence) | 84.4946 ns [84.1955, 84.9341] (0.95 confidence) |
| d2/la&#95;stack&#95;lu&#95;solve | 2.10013 ns [2.09514, 2.10964] (0.95 confidence) | 2.08604 ns [2.07011, 2.10324] (0.95 confidence) | 4.67756 ns [4.66983, 4.68294] (0.95 confidence) | 134.191 ns [133.728, 134.733] (0.95 confidence) |
| d2/la&#95;stack&#95;norm2&#95;sq | 0.396282 ns [0.395772, 0.397039] (0.95 confidence) | 0.401382 ns [0.400177, 0.402878] (0.95 confidence) | 0.396352 ns [0.39566, 0.397058] (0.95 confidence) | 4.18786 ns [4.18646, 4.19618] (0.95 confidence) |
| d2/la&#95;stack&#95;solve&#95;from&#95;ldlt | 1.27006 ns [1.26598, 1.27132] (0.95 confidence) | 1.26198 ns [1.25696, 1.26744] (0.95 confidence) | 1.29209 ns [1.28845, 1.29526] (0.95 confidence) | 34.9752 ns [34.816, 35.0867] (0.95 confidence) |
| d2/la&#95;stack&#95;solve&#95;from&#95;lu | 1.29707 ns [1.29403, 1.30025] (0.95 confidence) | 1.29248 ns [1.28377, 1.29784] (0.95 confidence) | 2.9914 ns [2.98343, 3.00375] (0.95 confidence) | 42.9799 ns [42.824, 43.151] (0.95 confidence) |
| d3/la&#95;stack&#95;det | 1.3452 ns [1.33846, 1.34912] (0.95 confidence) | 1.35317 ns [1.34688, 1.36196] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;det&#95;from&#95;ldlt | 0.594023 ns [0.59398, 0.594086] (0.95 confidence) | 0.604845 ns [0.602909, 0.606165] (0.95 confidence) | 0.485039 ns [0.484764, 0.485341] (0.95 confidence) | 0.742339 ns [0.742234, 0.742411] (0.95 confidence) |
| d3/la&#95;stack&#95;det&#95;from&#95;lu | 0.693126 ns [0.692947, 0.693409] (0.95 confidence) | 0.703788 ns [0.702007, 0.70615] (0.95 confidence) | 0.474913 ns [0.474849, 0.474982] (0.95 confidence) | 0.990024 ns [0.989822, 0.990483] (0.95 confidence) |
| d3/la&#95;stack&#95;det&#95;via&#95;lu | 8.77423 ns [8.75813, 8.8006] (0.95 confidence) | 8.92086 ns [8.8869, 8.95247] (0.95 confidence) | 13.1366 ns [13.108, 13.2644] (0.95 confidence) | 112.513 ns [111.955, 113.031] (0.95 confidence) |
| d3/la&#95;stack&#95;dot | 0.664367 ns [0.66369, 0.665379] (0.95 confidence) | 0.676309 ns [0.674337, 0.678252] (0.95 confidence) | 0.683378 ns [0.682914, 0.683667] (0.95 confidence) | 2.72049 ns [2.55765, 2.72198] (0.95 confidence) |
| d3/la&#95;stack&#95;inf&#95;norm | 1.36208 ns [1.36199, 1.3622] (0.95 confidence) | 1.3786 ns [1.37657, 1.38027] (0.95 confidence) | 1.17687 ns [1.17677, 1.17696] (0.95 confidence) | 1.33943 ns [1.32765, 1.3461] (0.95 confidence) |
| d3/la&#95;stack&#95;ldlt | 13.3075 ns [13.2508, 13.3931] (0.95 confidence) | 13.6252 ns [13.5597, 13.6886] (0.95 confidence) | 3.95167 ns [3.94481, 3.95668] (0.95 confidence) | 84.2451 ns [84.0766, 84.3537] (0.95 confidence) |
| d3/la&#95;stack&#95;ldlt&#95;solve | 9.79675 ns [9.77665, 9.83857] (0.95 confidence) | 9.84898 ns [9.77585, 9.93792] (0.95 confidence) | 5.72724 ns [5.71149, 5.73425] (0.95 confidence) | 124.904 ns [124.791, 125.083] (0.95 confidence) |
| d3/la&#95;stack&#95;lu | 8.40997 ns [8.39921, 8.42334] (0.95 confidence) | 8.46337 ns [8.44305, 8.4874] (0.95 confidence) | 12.324 ns [12.1421, 12.4366] (0.95 confidence) | 106.969 ns [106.816, 107.43] (0.95 confidence) |
| d3/la&#95;stack&#95;lu&#95;solve | 9.92161 ns [9.91361, 9.9287] (0.95 confidence) | 10.1267 ns [10.0964, 10.1596] (0.95 confidence) | 32.6591 ns [32.5867, 32.7761] (0.95 confidence) | 154.854 ns [154.671, 155.247] (0.95 confidence) |
| d3/la&#95;stack&#95;norm2&#95;sq | 0.44411 ns [0.443798, 0.444522] (0.95 confidence) | 0.454869 ns [0.453578, 0.457309] (0.95 confidence) | 0.445103 ns [0.444813, 0.445429] (0.95 confidence) | 4.12934 ns [4.1273, 4.13111] (0.95 confidence) |
| d3/la&#95;stack&#95;solve&#95;from&#95;ldlt | 1.77795 ns [1.77767, 1.77892] (0.95 confidence) | 1.80842 ns [1.80526, 1.81416] (0.95 confidence) | 2.83668 ns [2.83173, 2.84155] (0.95 confidence) | 35.3107 ns [35.0784, 35.4992] (0.95 confidence) |
| d3/la&#95;stack&#95;solve&#95;from&#95;lu | 2.0399 ns [2.0391, 2.04149] (0.95 confidence) | 2.09998 ns [2.08953, 2.10686] (0.95 confidence) | 4.54803 ns [4.54059, 4.55782] (0.95 confidence) | 43.3184 ns [43.0145, 43.5935] (0.95 confidence) |
| d32/la&#95;stack&#95;det | 2182 ns [2175.58, 2187.04] (0.95 confidence) | 2223.08 ns [2218.7, 2227.55] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;det&#95;from&#95;ldlt | 7.69107 ns [7.67645, 7.7027] (0.95 confidence) | 8.03731 ns [8.00739, 8.05051] (0.95 confidence) | 2.94299 ns [2.9393, 2.94622] (0.95 confidence) | 8.27568 ns [8.25958, 8.29245] (0.95 confidence) |
| d32/la&#95;stack&#95;det&#95;from&#95;lu | 8.977 ns [8.97317, 8.99046] (0.95 confidence) | 9.3028 ns [9.25894, 9.3643] (0.95 confidence) | 3.05277 ns [3.04703, 3.05856] (0.95 confidence) | 8.67523 ns [8.64858, 8.69204] (0.95 confidence) |
| d32/la&#95;stack&#95;det&#95;via&#95;lu | 2045.35 ns [2041.17, 2046.43] (0.95 confidence) | 2065.25 ns [2061.39, 2071.94] (0.95 confidence) | 2269.24 ns [2261.85, 2276.18] (0.95 confidence) | 2319.85 ns [2315.48, 2326.17] (0.95 confidence) |
| d32/la&#95;stack&#95;dot | 3.98429 ns [3.97772, 3.9934] (0.95 confidence) | 4.08529 ns [4.07284, 4.09232] (0.95 confidence) | 3.98987 ns [3.98694, 3.99393] (0.95 confidence) | 4.83929 ns [4.83707, 4.84133] (0.95 confidence) |
| d32/la&#95;stack&#95;inf&#95;norm | 129.791 ns [129.619, 130.037] (0.95 confidence) | 131.099 ns [130.947, 131.372] (0.95 confidence) | 161.098 ns [159.761, 163.155] (0.95 confidence) | 164.596 ns [164.063, 164.812] (0.95 confidence) |
| d32/la&#95;stack&#95;ldlt | 2377.19 ns [2375.64, 2378.38] (0.95 confidence) | 2397.51 ns [2395.36, 2400.46] (0.95 confidence) | 2071.76 ns [2064.14, 2078.73] (0.95 confidence) | 1359.7 ns [1356.07, 1364.35] (0.95 confidence) |
| d32/la&#95;stack&#95;ldlt&#95;solve | 2688.18 ns [2685.55, 2692.42] (0.95 confidence) | 2740.05 ns [2729.13, 2752.71] (0.95 confidence) | 2822.03 ns [2812.14, 2831.33] (0.95 confidence) | 1863.74 ns [1859.28, 1870.81] (0.95 confidence) |
| d32/la&#95;stack&#95;lu | 2198.74 ns [2193.65, 2206.81] (0.95 confidence) | 2248.36 ns [2232.33, 2267.13] (0.95 confidence) | 2135.75 ns [2131.1, 2138.6] (0.95 confidence) | 2269.67 ns [2266.36, 2274.27] (0.95 confidence) |
| d32/la&#95;stack&#95;lu&#95;solve | 2741.43 ns [2735.77, 2745.53] (0.95 confidence) | 2808.97 ns [2803.87, 2818.34] (0.95 confidence) | 2690.31 ns [2670.38, 2705.41] (0.95 confidence) | 2896.79 ns [2891.24, 2903.54] (0.95 confidence) |
| d32/la&#95;stack&#95;norm2&#95;sq | 3.92321 ns [3.91686, 3.92944] (0.95 confidence) | 4.03852 ns [4.03006, 4.0452] (0.95 confidence) | 3.58156 ns [3.57374, 3.58809] (0.95 confidence) | 4.20525 ns [4.20511, 4.20554] (0.95 confidence) |
| d32/la&#95;stack&#95;solve&#95;from&#95;ldlt | 308.973 ns [308.386, 309.722] (0.95 confidence) | 315.168 ns [314.772, 316.089] (0.95 confidence) | 555.61 ns [554.315, 557.585] (0.95 confidence) | 466.763 ns [466.289, 467.357] (0.95 confidence) |
| d32/la&#95;stack&#95;solve&#95;from&#95;lu | 656.305 ns [653.956, 659.365] (0.95 confidence) | 697.839 ns [696.782, 699.084] (0.95 confidence) | 327.023 ns [326.009, 327.433] (0.95 confidence) | 627.953 ns [627.26, 628.939] (0.95 confidence) |
| d4/la&#95;stack&#95;det | 2.55314 ns [2.55229, 2.55395] (0.95 confidence) | 2.58978 ns [2.57461, 2.60745] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;det&#95;from&#95;ldlt | 0.767416 ns [0.767354, 0.767469] (0.95 confidence) | 0.777446 ns [0.775474, 0.778829] (0.95 confidence) | 0.483954 ns [0.479339, 0.487772] (0.95 confidence) | 0.991797 ns [0.99148, 0.992335] (0.95 confidence) |
| d4/la&#95;stack&#95;det&#95;from&#95;lu | 0.867575 ns [0.867409, 0.867884] (0.95 confidence) | 0.875126 ns [0.873052, 0.877122] (0.95 confidence) | 0.587088 ns [0.583528, 0.592817] (0.95 confidence) | 1.23724 ns [1.23702, 1.23742] (0.95 confidence) |
| d4/la&#95;stack&#95;det&#95;via&#95;lu | 14.201 ns [14.1908, 14.2162] (0.95 confidence) | 14.4292 ns [14.3915, 14.451] (0.95 confidence) | 32.0753 ns [31.9999, 32.3066] (0.95 confidence) | 132.607 ns [132.3, 132.863] (0.95 confidence) |
| d4/la&#95;stack&#95;dot | 0.679943 ns [0.679453, 0.680522] (0.95 confidence) | 0.684804 ns [0.682462, 0.686139] (0.95 confidence) | 0.653494 ns [0.653246, 0.653695] (0.95 confidence) | 2.48786 ns [2.48554, 2.49104] (0.95 confidence) |
| d4/la&#95;stack&#95;inf&#95;norm | 2.34899 ns [2.34712, 2.35047] (0.95 confidence) | 2.38909 ns [2.38083, 2.39082] (0.95 confidence) | 2.02545 ns [2.02312, 2.02753] (0.95 confidence) | 2.04242 ns [2.04234, 2.04263] (0.95 confidence) |
| d4/la&#95;stack&#95;ldlt | 18.5227 ns [18.4585, 18.5709] (0.95 confidence) | 18.1184 ns [18.0769, 18.1652] (0.95 confidence) | 7.65423 ns [7.64598, 7.6731] (0.95 confidence) | 105.356 ns [105.087, 105.612] (0.95 confidence) |
| d4/la&#95;stack&#95;ldlt&#95;solve | 17.7808 ns [17.7637, 17.7899] (0.95 confidence) | 17.9666 ns [17.9392, 17.9894] (0.95 confidence) | 10.6718 ns [10.6558, 10.6911] (0.95 confidence) | 147.393 ns [147.018, 147.715] (0.95 confidence) |
| d4/la&#95;stack&#95;lu | 13.7463 ns [13.7391, 13.7544] (0.95 confidence) | 13.9351 ns [13.8969, 13.9667] (0.95 confidence) | 29.644 ns [29.6321, 29.6551] (0.95 confidence) | 127.412 ns [127.175, 127.693] (0.95 confidence) |
| d4/la&#95;stack&#95;lu&#95;solve | 21.8137 ns [21.772, 21.8502] (0.95 confidence) | 22.3772 ns [22.2802, 22.459] (0.95 confidence) | 53.0406 ns [52.9385, 53.1043] (0.95 confidence) | 180.446 ns [179.949, 180.87] (0.95 confidence) |
| d4/la&#95;stack&#95;norm2&#95;sq | 0.467956 ns [0.464323, 0.470412] (0.95 confidence) | 0.459498 ns [0.458269, 0.460469] (0.95 confidence) | 0.536845 ns [0.535356, 0.538889] (0.95 confidence) | 4.13013 ns [4.12802, 4.13245] (0.95 confidence) |
| d4/la&#95;stack&#95;solve&#95;from&#95;ldlt | 2.46601 ns [2.46472, 2.4675] (0.95 confidence) | 2.49014 ns [2.4848, 2.49479] (0.95 confidence) | 5.37537 ns [5.35887, 5.38824] (0.95 confidence) | 38.0564 ns [37.9689, 38.1541] (0.95 confidence) |
| d4/la&#95;stack&#95;solve&#95;from&#95;lu | 3.91906 ns [3.90746, 3.93628] (0.95 confidence) | 3.97481 ns [3.96095, 3.989] (0.95 confidence) | 5.71357 ns [5.70559, 5.71961] (0.95 confidence) | 46.7493 ns [46.555, 46.8708] (0.95 confidence) |
| d5/la&#95;stack&#95;det | 40.5954 ns [40.4881, 40.6886] (0.95 confidence) | 41.0598 ns [40.8595, 41.1037] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;det&#95;from&#95;ldlt | 0.875865 ns [0.875681, 0.876044] (0.95 confidence) | 0.879522 ns [0.878943, 0.880689] (0.95 confidence) | 0.579098 ns [0.577543, 0.580635] (0.95 confidence) | 1.22362 ns [1.22259, 1.2253] (0.95 confidence) |
| d5/la&#95;stack&#95;det&#95;from&#95;lu | 1.02922 ns [1.02904, 1.02934] (0.95 confidence) | 1.03525 ns [1.03377, 1.03593] (0.95 confidence) | 0.653228 ns [0.652906, 0.656333] (0.95 confidence) | 1.48377 ns [1.48143, 1.48413] (0.95 confidence) |
| d5/la&#95;stack&#95;det&#95;via&#95;lu | 33.4002 ns [33.2745, 33.4551] (0.95 confidence) | 33.5907 ns [33.4759, 33.6353] (0.95 confidence) | 50.8343 ns [50.6228, 51.0608] (0.95 confidence) | 170.689 ns [170.243, 171.072] (0.95 confidence) |
| d5/la&#95;stack&#95;dot | 0.841179 ns [0.838923, 0.842128] (0.95 confidence) | 0.844329 ns [0.84204, 0.852426] (0.95 confidence) | 0.762545 ns [0.761668, 0.764225] (0.95 confidence) | 2.63178 ns [2.63138, 2.63254] (0.95 confidence) |
| d5/la&#95;stack&#95;inf&#95;norm | 3.44479 ns [3.43898, 3.4515] (0.95 confidence) | 3.5738 ns [3.56102, 3.58815] (0.95 confidence) | 3.23395 ns [3.2273, 3.23911] (0.95 confidence) | 3.15523 ns [3.15319, 3.15562] (0.95 confidence) |
| d5/la&#95;stack&#95;ldlt | 24.93 ns [24.906, 24.9826] (0.95 confidence) | 25.4578 ns [25.4307, 25.486] (0.95 confidence) | 22.7947 ns [21.6891, 23.9297] (0.95 confidence) | 125.613 ns [125.245, 126.079] (0.95 confidence) |
| d5/la&#95;stack&#95;ldlt&#95;solve | 30.2591 ns [30.2128, 30.3451] (0.95 confidence) | 30.4307 ns [30.3607, 30.4917] (0.95 confidence) | 40.0729 ns [40.0625, 40.0837] (0.95 confidence) | 184.375 ns [184.027, 184.871] (0.95 confidence) |
| d5/la&#95;stack&#95;lu | 33.5307 ns [33.484, 33.5691] (0.95 confidence) | 33.6523 ns [33.5058, 33.7452] (0.95 confidence) | 49.6228 ns [49.3422, 49.8434] (0.95 confidence) | 162.362 ns [162.013, 162.886] (0.95 confidence) |
| d5/la&#95;stack&#95;lu&#95;solve | 48.7029 ns [48.5833, 48.8509] (0.95 confidence) | 48.8084 ns [48.777, 48.8493] (0.95 confidence) | 73.0617 ns [72.815, 73.8988] (0.95 confidence) | 257.861 ns [256.801, 259.666] (0.95 confidence) |
| d5/la&#95;stack&#95;norm2&#95;sq | 0.523922 ns [0.523551, 0.524413] (0.95 confidence) | 0.524923 ns [0.524289, 0.526392] (0.95 confidence) | 0.50218 ns [0.502061, 0.502439] (0.95 confidence) | 4.11539 ns [4.11368, 4.11693] (0.95 confidence) |
| d5/la&#95;stack&#95;solve&#95;from&#95;ldlt | 3.86968 ns [3.86118, 3.87836] (0.95 confidence) | 3.93299 ns [3.91599, 3.93909] (0.95 confidence) | 8.83247 ns [8.82543, 8.84179] (0.95 confidence) | 57.1795 ns [57.1152, 57.2948] (0.95 confidence) |
| d5/la&#95;stack&#95;solve&#95;from&#95;lu | 6.08215 ns [6.07524, 6.08969] (0.95 confidence) | 6.12587 ns [6.11194, 6.1346] (0.95 confidence) | 8.8573 ns [8.83692, 8.86672] (0.95 confidence) | 78.1927 ns [77.8529, 78.45] (0.95 confidence) |
| d64/la&#95;stack&#95;det | 15099 ns [15011.2, 15199] (0.95 confidence) | 14814.9 ns [14781.9, 14851.7] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;det&#95;from&#95;ldlt | 23.4905 ns [23.3812, 23.5796] (0.95 confidence) | 23.6613 ns [23.6266, 23.6998] (0.95 confidence) | 8.57614 ns [8.54886, 8.59401] (0.95 confidence) | 21.2024 ns [21.1229, 21.2797] (0.95 confidence) |
| d64/la&#95;stack&#95;det&#95;from&#95;lu | 23.8202 ns [23.7592, 23.8701] (0.95 confidence) | 23.8007 ns [23.7531, 23.8604] (0.95 confidence) | 8.85422 ns [8.83878, 8.87139] (0.95 confidence) | 22.3252 ns [22.2381, 22.4876] (0.95 confidence) |
| d64/la&#95;stack&#95;det&#95;via&#95;lu | 14373 ns [14321.7, 14423.3] (0.95 confidence) | 14279.6 ns [14264, 14301.2] (0.95 confidence) | 13582.4 ns [13563.7, 13616.6] (0.95 confidence) | 11190.7 ns [11146.1, 11258.7] (0.95 confidence) |
| d64/la&#95;stack&#95;dot | 10.9769 ns [10.955, 10.9995] (0.95 confidence) | 11.2595 ns [11.2297, 11.2853] (0.95 confidence) | 8.04063 ns [8.01696, 8.06084] (0.95 confidence) | 8.08276 ns [8.07748, 8.09128] (0.95 confidence) |
| d64/la&#95;stack&#95;inf&#95;norm | 662.324 ns [661.753, 662.849] (0.95 confidence) | 644.322 ns [642.754, 645.656] (0.95 confidence) | 1338 ns [1336.57, 1340.01] (0.95 confidence) | 1594.46 ns [1591.27, 1598.52] (0.95 confidence) |
| d64/la&#95;stack&#95;ldlt | 17220.6 ns [17166.6, 17244.9] (0.95 confidence) | 16856.6 ns [16848.1, 16876.3] (0.95 confidence) | 10734.5 ns [10714.4, 10767] (0.95 confidence) | 7155.37 ns [7134.26, 7170.86] (0.95 confidence) |
| d64/la&#95;stack&#95;ldlt&#95;solve | 19129.7 ns [19078.3, 19167.5] (0.95 confidence) | 19315.9 ns [19280.2, 19362.2] (0.95 confidence) | 12244.3 ns [12226.7, 12260.1] (0.95 confidence) | 8358.63 ns [8350.02, 8376.27] (0.95 confidence) |
| d64/la&#95;stack&#95;lu | 15522.6 ns [15484.8, 15552.6] (0.95 confidence) | 15190.5 ns [15164.8, 15236.7] (0.95 confidence) | 12743.8 ns [12703, 12757.6] (0.95 confidence) | 11152.5 ns [11142.7, 11181.3] (0.95 confidence) |
| d64/la&#95;stack&#95;lu&#95;solve | 17300.2 ns [17256.8, 17337.8] (0.95 confidence) | 17282 ns [17249.2, 17323.3] (0.95 confidence) | 14432 ns [14407.6, 14464.9] (0.95 confidence) | 12989.7 ns [12951.4, 13040.3] (0.95 confidence) |
| d64/la&#95;stack&#95;norm2&#95;sq | 11.4593 ns [11.3548, 11.5663] (0.95 confidence) | 11.1089 ns [11.0989, 11.1307] (0.95 confidence) | 8.93083 ns [8.88963, 8.96048] (0.95 confidence) | 6.41099 ns [6.26565, 6.43075] (0.95 confidence) |
| d64/la&#95;stack&#95;solve&#95;from&#95;ldlt | 1115.47 ns [1112.76, 1120.19] (0.95 confidence) | 1098.07 ns [1094.66, 1102.1] (0.95 confidence) | 1321.98 ns [1318.65, 1325.22] (0.95 confidence) | 1270.44 ns [1262.76, 1276.09] (0.95 confidence) |
| d64/la&#95;stack&#95;solve&#95;from&#95;lu | 2669.11 ns [2624.77, 2712.47] (0.95 confidence) | 2611.06 ns [2605.12, 2621.88] (0.95 confidence) | 813.067 ns [810.889, 817.177] (0.95 confidence) | 1706.62 ns [1703.31, 1709.81] (0.95 confidence) |
| d8/la&#95;stack&#95;det | 94.4297 ns [93.8917, 94.7198] (0.95 confidence) | 96.5096 ns [95.8312, 96.8879] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;det&#95;from&#95;ldlt | 1.31833 ns [1.3161, 1.32035] (0.95 confidence) | 1.32512 ns [1.32139, 1.33136] (0.95 confidence) | 0.925975 ns [0.920341, 0.935258] (0.95 confidence) | 1.99876 ns [1.99712, 2.00292] (0.95 confidence) |
| d8/la&#95;stack&#95;det&#95;from&#95;ldlt&#95;balanced&#95;range | 8.92545 ns [8.91765, 8.93489] (0.95 confidence) | 8.93058 ns [8.92071, 8.9554] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;det&#95;from&#95;lu | 1.47907 ns [1.47769, 1.48124] (0.95 confidence) | 1.50032 ns [1.49594, 1.5065] (0.95 confidence) | 1.00522 ns [1.00379, 1.00679] (0.95 confidence) | 2.25317 ns [2.25041, 2.2595] (0.95 confidence) |
| d8/la&#95;stack&#95;det&#95;from&#95;lu&#95;balanced&#95;range | 9.21361 ns [9.18956, 9.23661] (0.95 confidence) | 9.23756 ns [9.22415, 9.29228] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;det&#95;via&#95;lu | 90.6719 ns [90.0715, 91.2476] (0.95 confidence) | 88.3343 ns [87.2426, 89.3322] (0.95 confidence) | 137.8 ns [137.575, 137.945] (0.95 confidence) | 243.153 ns [242.647, 243.861] (0.95 confidence) |
| d8/la&#95;stack&#95;dot | 0.940209 ns [0.939374, 0.941332] (0.95 confidence) | 0.927324 ns [0.924937, 0.929722] (0.95 confidence) | 1.36474 ns [1.36288, 1.36772] (0.95 confidence) | 2.61538 ns [2.60534, 2.62536] (0.95 confidence) |
| d8/la&#95;stack&#95;inf&#95;norm | 8.48378 ns [8.47377, 8.49436] (0.95 confidence) | 8.55112 ns [8.54346, 8.56367] (0.95 confidence) | 8.08498 ns [8.06452, 8.09826] (0.95 confidence) | 8.13716 ns [8.12354, 8.14716] (0.95 confidence) |
| d8/la&#95;stack&#95;ldlt | 90.9843 ns [90.5604, 91.4939] (0.95 confidence) | 91.3422 ns [91.0941, 91.6294] (0.95 confidence) | 104.214 ns [103.545, 104.746] (0.95 confidence) | 186.501 ns [186.112, 187.032] (0.95 confidence) |
| d8/la&#95;stack&#95;ldlt&#95;ill&#95;conditioned | 90.7093 ns [90.5635, 90.8494] (0.95 confidence) | 91.1547 ns [91.0196, 91.3356] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;ldlt&#95;solve | 103.694 ns [103.47, 104.097] (0.95 confidence) | 104.466 ns [104.038, 104.878] (0.95 confidence) | 149.394 ns [143.506, 150.493] (0.95 confidence) | 261.824 ns [261.002, 263.752] (0.95 confidence) |
| d8/la&#95;stack&#95;lu | 89.6213 ns [89.4622, 89.806] (0.95 confidence) | 90.4649 ns [89.9379, 91.1991] (0.95 confidence) | 133.806 ns [133.458, 134.063] (0.95 confidence) | 230.049 ns [229.557, 231.408] (0.95 confidence) |
| d8/la&#95;stack&#95;lu&#95;ill&#95;conditioned | 89.9146 ns [89.1801, 90.2248] (0.95 confidence) | 90.4408 ns [89.551, 90.8016] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;lu&#95;pivoting | 96.3611 ns [96.1327, 96.7414] (0.95 confidence) | 95.4818 ns [95.1756, 95.8917] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;lu&#95;solve | 130.313 ns [129.926, 130.665] (0.95 confidence) | 133.645 ns [133.022, 134.192] (0.95 confidence) | 165.323 ns [164.976, 165.544] (0.95 confidence) | 343.132 ns [341.19, 344.925] (0.95 confidence) |
| d8/la&#95;stack&#95;norm2&#95;sq | 0.656376 ns [0.654619, 0.659467] (0.95 confidence) | 0.659491 ns [0.656381, 0.661156] (0.95 confidence) | 0.940483 ns [0.939902, 0.941394] (0.95 confidence) | 4.17276 ns [4.16914, 4.17632] (0.95 confidence) |
| d8/la&#95;stack&#95;solve&#95;from&#95;ldlt | 8.23446 ns [8.21941, 8.26455] (0.95 confidence) | 8.22664 ns [8.2065, 8.25128] (0.95 confidence) | 21.8696 ns [21.7953, 21.929] (0.95 confidence) | 68.6217 ns [68.5514, 68.6858] (0.95 confidence) |
| d8/la&#95;stack&#95;solve&#95;from&#95;lu | 13.6594 ns [13.632, 13.6841] (0.95 confidence) | 13.5354 ns [13.5057, 13.5614] (0.95 confidence) | 16.3345 ns [16.2846, 16.3732] (0.95 confidence) | 96.1579 ns [95.9001, 96.758] (0.95 confidence) |
| exact&#95;d2/det | 0.447061 ns [0.445505, 0.448842] (0.95 confidence) | 0.436911 ns [0.435776, 0.438239] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;direct | 0.446168 ns [0.44433, 0.448788] (0.95 confidence) | 0.434805 ns [0.434504, 0.435449] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;direct&#95;with&#95;errbound | 1.76545 ns [1.76255, 1.76855] (0.95 confidence) | 1.74506 ns [1.7421, 1.74715] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;errbound | 1.61193 ns [1.60877, 1.61366] (0.95 confidence) | 1.59569 ns [1.59216, 1.59849] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;exact | 86.2495 ns [85.8077, 86.4668] (0.95 confidence) | 88.4997 ns [88.3078, 88.7312] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;exact&#95;f64&#95;result | 60.781 ns [60.6175, 60.9605] (0.95 confidence) | 65.494 ns [65.4462, 65.5977] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;exact&#95;rounded&#95;f64 | 62.2712 ns [62.079, 62.445] (0.95 confidence) | 68.1834 ns [67.9814, 68.3201] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;sign&#95;exact | 3.03183 ns [3.02998, 3.03508] (0.95 confidence) | 2.98818 ns [2.98512, 2.99248] (0.95 confidence) | — | — |
| exact&#95;d2/solve&#95;exact | 7348.2 ns [7316.64, 7627.86] (0.95 confidence) | 3016.03 ns [3010.21, 3024.57] (0.95 confidence) | — | — |
| exact&#95;d2/solve&#95;exact&#95;f64&#95;result | 7291.2 ns [7279.36, 7304.6] (0.95 confidence) | 3156.23 ns [3152.6, 3170.15] (0.95 confidence) | — | — |
| exact&#95;d2/solve&#95;exact&#95;rounded&#95;f64 | 7526.79 ns [7499.52, 7543.58] (0.95 confidence) | 3341.95 ns [3335.58, 3346.4] (0.95 confidence) | — | — |
| exact&#95;d3/det | 0.8458 ns [0.841489, 0.848934] (0.95 confidence) | 0.825558 ns [0.824, 0.827777] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;direct | 0.828476 ns [0.828314, 0.828806] (0.95 confidence) | 0.831781 ns [0.829159, 0.834748] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;direct&#95;with&#95;errbound | 3.35589 ns [3.3517, 3.3589] (0.95 confidence) | 3.34142 ns [3.33546, 3.34569] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;errbound | 3.33839 ns [3.33525, 3.34047] (0.95 confidence) | 3.33136 ns [3.32587, 3.33689] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;exact | 291.945 ns [291.461, 292.362] (0.95 confidence) | 302.343 ns [302.011, 302.631] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;exact&#95;f64&#95;result | 265.915 ns [264.984, 266.697] (0.95 confidence) | 277.97 ns [277.365, 278.315] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;exact&#95;rounded&#95;f64 | 271.053 ns [270.821, 271.521] (0.95 confidence) | 282.347 ns [281.351, 282.815] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;sign&#95;exact | 4.78565 ns [4.77962, 4.79698] (0.95 confidence) | 4.74248 ns [4.73657, 4.75222] (0.95 confidence) | — | — |
| exact&#95;d3/solve&#95;exact | 30337.2 ns [30306.6, 30360.4] (0.95 confidence) | 9194.78 ns [9159.27, 9246.59] (0.95 confidence) | — | — |
| exact&#95;d3/solve&#95;exact&#95;f64&#95;result | 30433.3 ns [30315.1, 30501.7] (0.95 confidence) | 9280.51 ns [9250.27, 9324.82] (0.95 confidence) | — | — |
| exact&#95;d3/solve&#95;exact&#95;rounded&#95;f64 | 30809.1 ns [30758.5, 30861.8] (0.95 confidence) | 9597.44 ns [9559.69, 9639.71] (0.95 confidence) | — | — |
| exact&#95;d4/det | 3.28729 ns [3.27918, 3.2968] (0.95 confidence) | 3.23082 ns [3.22817, 3.23488] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;direct | 2.42967 ns [2.42463, 2.43468] (0.95 confidence) | 2.36794 ns [2.36618, 2.37203] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;direct&#95;with&#95;errbound | 6.85923 ns [6.8384, 6.86431] (0.95 confidence) | 6.82974 ns [6.81826, 6.8489] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;errbound | 6.79328 ns [6.78353, 6.80274] (0.95 confidence) | 6.81918 ns [6.80793, 6.82838] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;exact | 811.978 ns [804.182, 830.689] (0.95 confidence) | 853.407 ns [850.642, 854.307] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;exact&#95;f64&#95;result | 793.076 ns [783.916, 817.067] (0.95 confidence) | 811.592 ns [807.415, 814.233] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;exact&#95;rounded&#95;f64 | 777.651 ns [775.258, 781.116] (0.95 confidence) | 839.889 ns [835.828, 844.656] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;sign&#95;exact | 7.99175 ns [7.98779, 8.00063] (0.95 confidence) | 8.06773 ns [8.04558, 8.0854] (0.95 confidence) | — | — |
| exact&#95;d4/solve&#95;exact | 79011.1 ns [78417, 79339.6] (0.95 confidence) | 18473.2 ns [18429.9, 18500.8] (0.95 confidence) | — | — |
| exact&#95;d4/solve&#95;exact&#95;f64&#95;result | 78971.9 ns [78834.5, 79203.3] (0.95 confidence) | 18602.9 ns [18579, 18634.7] (0.95 confidence) | — | — |
| exact&#95;d4/solve&#95;exact&#95;rounded&#95;f64 | 79421.8 ns [79255.7, 79534] (0.95 confidence) | 19189.3 ns [19165.4, 19210.5] (0.95 confidence) | — | — |
| exact&#95;d5/det | 28.0208 ns [27.9559, 28.1076] (0.95 confidence) | 27.492 ns [27.4697, 27.5499] (0.95 confidence) | — | — |
| exact&#95;d5/det&#95;exact | 2940.38 ns [2936.46, 2947.21] (0.95 confidence) | 2871.87 ns [2860.42, 2881.17] (0.95 confidence) | — | — |
| exact&#95;d5/det&#95;exact&#95;f64&#95;result | 2895.01 ns [2891.48, 2899.19] (0.95 confidence) | 2825.57 ns [2821.31, 2832.48] (0.95 confidence) | — | — |
| exact&#95;d5/det&#95;exact&#95;rounded&#95;f64 | 2927.16 ns [2920.96, 2934.79] (0.95 confidence) | 2898.08 ns [2887.16, 2916.11] (0.95 confidence) | — | — |
| exact&#95;d5/det&#95;sign&#95;exact | 2901.26 ns [2894.54, 2910.68] (0.95 confidence) | 2848.02 ns [2843.84, 2853.08] (0.95 confidence) | — | — |
| exact&#95;d5/solve&#95;exact | 155170 ns [154910, 155363] (0.95 confidence) | 29911.1 ns [29841.4, 30037.4] (0.95 confidence) | — | — |
| exact&#95;d5/solve&#95;exact&#95;f64&#95;result | 155426 ns [155150, 155794] (0.95 confidence) | 30179.8 ns [30082.2, 30333.5] (0.95 confidence) | — | — |
| exact&#95;d5/solve&#95;exact&#95;rounded&#95;f64 | 156296 ns [155976, 156565] (0.95 confidence) | 31008.2 ns [30903.7, 31094.6] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/det&#95;exact | 886.169 ns [882.918, 890.396] (0.95 confidence) | 918.798 ns [916.54, 920.971] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/det&#95;sign&#95;exact | 8.02681 ns [8.01052, 8.03637] (0.95 confidence) | 7.93198 ns [7.92195, 7.93987] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/solve&#95;exact | 58849.7 ns [58608.2, 59031.7] (0.95 confidence) | 15232.3 ns [15183.6, 15313.7] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/solve&#95;exact&#95;f64&#95;result | 58908.4 ns [58763.4, 59018] (0.95 confidence) | 15606.1 ns [15571.2, 15639.9] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/solve&#95;exact&#95;rounded&#95;f64 | 59340.4 ns [59147.2, 59494.3] (0.95 confidence) | 15815.1 ns [15778.8, 15842.7] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/det&#95;exact | 3109.48 ns [3098.42, 3117.94] (0.95 confidence) | 3006.3 ns [2997.79, 3016.29] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/det&#95;sign&#95;exact | 3081.58 ns [3072.34, 3091.41] (0.95 confidence) | 2919 ns [2913.88, 2926.19] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/solve&#95;exact | 123197 ns [122852, 123553] (0.95 confidence) | 24903.4 ns [24849, 24948.5] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/solve&#95;exact&#95;f64&#95;result | 121815 ns [121626, 122186] (0.95 confidence) | 25144.3 ns [25080.6, 25204.6] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/solve&#95;exact&#95;rounded&#95;f64 | 123372 ns [122999, 123693] (0.95 confidence) | 25911.8 ns [25849.8, 25988.7] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/det&#95;exact | 311.991 ns [311.48, 312.654] (0.95 confidence) | 322.854 ns [322.226, 323.314] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/det&#95;sign&#95;exact | 287.26 ns [286.009, 288.378] (0.95 confidence) | 298.219 ns [296.831, 299.471] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/solve&#95;exact | 95480.2 ns [95244.8, 95775.3] (0.95 confidence) | 9092.68 ns [9069.45, 9115.61] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/solve&#95;exact&#95;f64&#95;result | 95462.4 ns [95313.1, 95677.7] (0.95 confidence) | 9165.22 ns [9161.29, 9171.24] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/solve&#95;exact&#95;rounded&#95;f64 | 95845.9 ns [95604.6, 96090] (0.95 confidence) | 9544.94 ns [9517.59, 9569.46] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/det&#95;exact | 296.858 ns [295.444, 297.515] (0.95 confidence) | 185.153 ns [184.825, 185.518] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/det&#95;sign&#95;exact | 321.698 ns [321.176, 322.72] (0.95 confidence) | 167.592 ns [167.296, 168.096] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/solve&#95;exact | 2472.61 ns [2466.18, 2480.96] (0.95 confidence) | 987.024 ns [983.077, 988.592] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/solve&#95;exact&#95;f64&#95;result | 2494.61 ns [2490.02, 2502.66] (0.95 confidence) | 1014.36 ns [1011.89, 1016.14] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/solve&#95;exact&#95;rounded&#95;f64 | 2517.85 ns [2512.63, 2521.55] (0.95 confidence) | 1014.35 ns [1008.42, 1018.01] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/det&#95;exact | 3108 ns [3101.57, 3110.89] (0.95 confidence) | 3470.69 ns [3461.42, 3479.63] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/det&#95;sign&#95;exact | 147.929 ns [147.371, 148.55] (0.95 confidence) | 143.384 ns [143.117, 143.599] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/solve&#95;exact | 75185.6 ns [75070, 75278.3] (0.95 confidence) | 35625 ns [35557.1, 35661.9] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/solve&#95;exact&#95;f64&#95;result | 76152.1 ns [75926, 76304.1] (0.95 confidence) | 36439.6 ns [36378.1, 36490.1] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/solve&#95;exact&#95;rounded&#95;f64 | 76193.4 ns [76029.7, 76711.2] (0.95 confidence) | 36649.3 ns [36544.2, 36793.5] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/det&#95;exact | 8937.02 ns [8928.87, 8951.56] (0.95 confidence) | 9408.85 ns [9399.03, 9416.69] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/det&#95;sign&#95;exact | 227.702 ns [227.337, 227.941] (0.95 confidence) | 225.451 ns [225.335, 225.569] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/solve&#95;exact | 235100 ns [234917, 235393] (0.95 confidence) | 81283.3 ns [80965.5, 81348.2] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/solve&#95;exact&#95;f64&#95;result | 236378 ns [236066, 236634] (0.95 confidence) | 81756 ns [81694.5, 81871.8] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/solve&#95;exact&#95;rounded&#95;f64 | 236568 ns [236431, 236871] (0.95 confidence) | 81895 ns [81837, 82009.9] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/det&#95;exact | 22304.4 ns [22257.5, 22382.6] (0.95 confidence) | 23096.5 ns [23052.7, 23217.6] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/det&#95;sign&#95;exact | 458.212 ns [457.679, 458.856] (0.95 confidence) | 443.819 ns [443.388, 444.319] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/solve&#95;exact | 558971 ns [558553, 559607] (0.95 confidence) | 144832 ns [144561, 144984] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/solve&#95;exact&#95;f64&#95;result | 560222 ns [559444, 561033] (0.95 confidence) | 145657 ns [145515, 145777] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/solve&#95;exact&#95;rounded&#95;f64 | 556988 ns [556105, 558198] (0.95 confidence) | 146541 ns [146435, 146671] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/det&#95;exact | 56497.4 ns [56262.1, 56699.9] (0.95 confidence) | 55785.3 ns [55721.7, 55853.3] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/det&#95;sign&#95;exact | 55844.3 ns [55626.3, 56004.5] (0.95 confidence) | 54821.9 ns [54720.9, 54962.5] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/solve&#95;exact | 1.11253e+06 ns [1.11113e+06, 1.11508e+06] (0.95 confidence) | 238283 ns [237402, 238931] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/solve&#95;exact&#95;f64&#95;result | 1.11324e+06 ns [1.10959e+06, 1.11532e+06] (0.95 confidence) | 238497 ns [238038, 238820] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/solve&#95;exact&#95;rounded&#95;f64 | 1.10892e+06 ns [1.10482e+06, 1.11365e+06] (0.95 confidence) | 239783 ns [239562, 240425] (0.95 confidence) | — | — |
| rational&#95;input&#95;d2/det&#95;big&#95;rational&#95;gaussian | 1124.88 ns [1121.92, 1126.54] (0.95 confidence) | 1112.35 ns [1110.98, 1114.29] (0.95 confidence) | — | — |
| rational&#95;input&#95;d2/det&#95;row&#95;cleared&#95;bareiss | 179.956 ns [179.493, 180.513] (0.95 confidence) | 180.64 ns [179.437, 182.352] (0.95 confidence) | — | — |
| rational&#95;input&#95;d2/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 77.4413 ns [77.228, 77.609] (0.95 confidence) | 77.414 ns [77.2865, 77.5376] (0.95 confidence) | — | — |
| rational&#95;input&#95;d2/solve&#95;big&#95;rational&#95;gaussian | 2149.64 ns [2142.37, 2154.45] (0.95 confidence) | 2128.03 ns [2123.92, 2132.87] (0.95 confidence) | — | — |
| rational&#95;input&#95;d2/solve&#95;row&#95;cleared&#95;bareiss | 934.448 ns [932.297, 936.493] (0.95 confidence) | 611.085 ns [610.542, 612.449] (0.95 confidence) | — | — |
| rational&#95;input&#95;d3/det&#95;big&#95;rational&#95;gaussian | 4760.34 ns [4749.4, 4771.38] (0.95 confidence) | 4749.5 ns [4728.58, 4781.11] (0.95 confidence) | — | — |
| rational&#95;input&#95;d3/det&#95;row&#95;cleared&#95;bareiss | 569.43 ns [568.248, 571.027] (0.95 confidence) | 570.294 ns [568.207, 572.797] (0.95 confidence) | — | — |
| rational&#95;input&#95;d3/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 360.808 ns [359.832, 361.702] (0.95 confidence) | 368.327 ns [364.907, 372.358] (0.95 confidence) | — | — |
| rational&#95;input&#95;d3/solve&#95;big&#95;rational&#95;gaussian | 8238.12 ns [8214.12, 8274.56] (0.95 confidence) | 8154.58 ns [8137.35, 8161.65] (0.95 confidence) | — | — |
| rational&#95;input&#95;d3/solve&#95;row&#95;cleared&#95;bareiss | 2470.93 ns [2455.06, 2487.43] (0.95 confidence) | 1242.7 ns [1241.43, 1244.97] (0.95 confidence) | — | — |
| rational&#95;input&#95;d4/det&#95;big&#95;rational&#95;gaussian | 14063.7 ns [14006.3, 14108.8] (0.95 confidence) | 13746.5 ns [13709.8, 13770.1] (0.95 confidence) | — | — |
| rational&#95;input&#95;d4/det&#95;row&#95;cleared&#95;bareiss | 1170.02 ns [1168.15, 1174.38] (0.95 confidence) | 1167.28 ns [1164.99, 1169.25] (0.95 confidence) | — | — |
| rational&#95;input&#95;d4/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 775.91 ns [773.321, 777.669] (0.95 confidence) | 780.759 ns [778.043, 784.302] (0.95 confidence) | — | — |
| rational&#95;input&#95;d4/solve&#95;big&#95;rational&#95;gaussian | 24087.4 ns [24045, 24167] (0.95 confidence) | 24061.4 ns [23974, 24165.9] (0.95 confidence) | — | — |
| rational&#95;input&#95;d4/solve&#95;row&#95;cleared&#95;bareiss | 6053.42 ns [5977.7, 6424.57] (0.95 confidence) | 2192.29 ns [2188.82, 2193.59] (0.95 confidence) | — | — |
| rational&#95;input&#95;d5/det&#95;big&#95;rational&#95;gaussian | 23772.3 ns [23706.2, 23852.9] (0.95 confidence) | 23913.2 ns [23827.5, 24005] (0.95 confidence) | — | — |
| rational&#95;input&#95;d5/det&#95;row&#95;cleared&#95;bareiss | 2055.73 ns [2053.82, 2058.3] (0.95 confidence) | 2075.87 ns [2072.09, 2078.93] (0.95 confidence) | — | — |
| rational&#95;input&#95;d5/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 1555.74 ns [1551.2, 1558.2] (0.95 confidence) | 1565.28 ns [1560.95, 1568.14] (0.95 confidence) | — | — |
| rational&#95;input&#95;d5/solve&#95;big&#95;rational&#95;gaussian | 40244.3 ns [40116.4, 40361.6] (0.95 confidence) | 39605.7 ns [39546.7, 39664.3] (0.95 confidence) | — | — |
| rational&#95;input&#95;d5/solve&#95;row&#95;cleared&#95;bareiss | 9572.44 ns [9545.8, 9586.93] (0.95 confidence) | 3389.42 ns [3374.06, 3399.43] (0.95 confidence) | — | — |
| rational&#95;input&#95;d6/det&#95;big&#95;rational&#95;gaussian | 56898 ns [56735.2, 56988.2] (0.95 confidence) | 58101.2 ns [57920.6, 58199.8] (0.95 confidence) | — | — |
| rational&#95;input&#95;d6/det&#95;row&#95;cleared&#95;bareiss | 3543.68 ns [3536.61, 3551.44] (0.95 confidence) | 3520.67 ns [3516.19, 3527.03] (0.95 confidence) | — | — |
| rational&#95;input&#95;d6/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 2860.43 ns [2844.73, 2866.39] (0.95 confidence) | 2907.82 ns [2886.86, 2924.07] (0.95 confidence) | — | — |
| rational&#95;input&#95;d6/solve&#95;big&#95;rational&#95;gaussian | 90819.1 ns [90512.9, 91105.1] (0.95 confidence) | 94520.2 ns [94151.1, 94771.5] (0.95 confidence) | — | — |
| rational&#95;input&#95;d6/solve&#95;row&#95;cleared&#95;bareiss | 17846.2 ns [17800.5, 17876.3] (0.95 confidence) | 5763.11 ns [5752.48, 5776.2] (0.95 confidence) | — | — |
| rational&#95;input&#95;d7/det&#95;big&#95;rational&#95;gaussian | 110709 ns [110499, 110885] (0.95 confidence) | 112144 ns [111938, 112464] (0.95 confidence) | — | — |
| rational&#95;input&#95;d7/det&#95;row&#95;cleared&#95;bareiss | 5611.26 ns [5594.35, 5620.66] (0.95 confidence) | 5594.87 ns [5586.85, 5602.95] (0.95 confidence) | — | — |
| rational&#95;input&#95;d7/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 4658.33 ns [4652.07, 4669.75] (0.95 confidence) | 4764.17 ns [4755.87, 4770.82] (0.95 confidence) | — | — |
| rational&#95;input&#95;d7/solve&#95;big&#95;rational&#95;gaussian | 162045 ns [161702, 162312] (0.95 confidence) | 167832 ns [167087, 168309] (0.95 confidence) | — | — |
| rational&#95;input&#95;d7/solve&#95;row&#95;cleared&#95;bareiss | 27621.5 ns [27601.4, 27675.2] (0.95 confidence) | 10111.2 ns [10097.6, 10128.9] (0.95 confidence) | — | — |
| rational&#95;input&#95;d8/det&#95;big&#95;rational&#95;gaussian | 188217 ns [187912, 188633] (0.95 confidence) | 193407 ns [193086, 194082] (0.95 confidence) | — | — |
| rational&#95;input&#95;d8/det&#95;row&#95;cleared&#95;bareiss | 8066.77 ns [8051.61, 8090.2] (0.95 confidence) | 8221.63 ns [8211.07, 8253.95] (0.95 confidence) | — | — |
| rational&#95;input&#95;d8/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 7240.46 ns [7227.54, 7253.83] (0.95 confidence) | 7377.06 ns [7363.19, 7401.44] (0.95 confidence) | — | — |
| rational&#95;input&#95;d8/solve&#95;big&#95;rational&#95;gaussian | 268034 ns [267183, 268878] (0.95 confidence) | 272292 ns [272046, 272658] (0.95 confidence) | — | — |
| rational&#95;input&#95;d8/solve&#95;row&#95;cleared&#95;bareiss | 40814.3 ns [40749.5, 40889.6] (0.95 confidence) | 15364.7 ns [15329.1, 15406.2] (0.95 confidence) | — | — |

## Recorded measurement provenance

These observations belong to the measured phases. Rerendering does not recapture the publisher's environment.

### Baseline phase

| Field | Recorded value |
| --- | --- |
| Release | v0.5.0 |
| Revision | 5d8f5b4ea1c0bd708cd748b86f426960b7597ab7 |
| Measurement mode | tag |
| Original source SHA-256 | 1d876cd1476f52283a2f9dbf1bdbf46baec54b86834115f1eddefee601ad6b68 |
| Measured source SHA-256 | 36efc0af7c97d4a2508f9ad26b692dff69296e1cb9a94878d323b767f98d78c7 |
| Common harness SHA-256 | 2a5139b9eae822eb41751dfb2003eb74f76b192b79c1acbcba9923927f3e8ee3 |
| Fingerprint schema | research-repo-tools/files/v1 |
| OS | macOS-27.0.1-arm64-arm-64bit-Mach-O |
| Architecture | arm64 |
| CPU | Apple M4 Max (arm64) |
| Physical cores | 16 |
| Logical threads | 16 |
| Memory (bytes) | 68719476736 |
| Tool: cargo | cargo 1.99.1 (5f94df478 2026-08-27) |
| Tool: just | just 1.59.0 |
| Tool: rustc | rustc 1.99.1 (58cf8f9c3 2026-10-13)<br>binary: rustc<br>commit-hash: 58cf8f9c35b355199727627d07ae34a4e2846d67<br>commit-date: 2026-10-13<br>host: aarch64-apple-darwin<br>release: 1.99.1<br>LLVM version: 23.1.1 |
| Dependency: criterion | 0.8.2 |
| Dependency: faer | 0.24.4 |
| Dependency: nalgebra | 0.35.0 |
| Environment overrides (names only) | CARGO&#95;ENCODED&#95;RUSTFLAGS, CARGO&#95;TARGET&#95;DIR, CRITERION&#95;HOME, RUSTC, RUSTC&#95;WORKSPACE&#95;WRAPPER, RUSTC&#95;WRAPPER, RUSTFLAGS, RUSTUP&#95;TOOLCHAIN |
| Measurement command (argv) | &#91;&quot;just&quot;, &quot;--justfile&quot;, &quot;tooling/performance.just&quot;, &quot;--working-directory&quot;, &quot;.&quot;, &quot;baseline-all&quot;&#93; |
| Independent correctness gate | passed |
| Gate command (argv) | &#91;&quot;just&quot;, &quot;--justfile&quot;, &quot;tooling/performance.just&quot;, &quot;--working-directory&quot;, &quot;.&quot;, &quot;gate&quot;&#93; |

### Current phase

| Field | Recorded value |
| --- | --- |
| Release | v0.5.1 |
| Revision | eb6f071607c058527f6cb6a9bf86bb307c93ab02 |
| Measurement mode | working-tree |
| Original source SHA-256 | 8edc38caa565097a74ebe12ad80493067168403c4e7bc8c4cadc1d653152892f |
| Measured source SHA-256 | 8edc38caa565097a74ebe12ad80493067168403c4e7bc8c4cadc1d653152892f |
| Common harness SHA-256 | 2a5139b9eae822eb41751dfb2003eb74f76b192b79c1acbcba9923927f3e8ee3 |
| Fingerprint schema | research-repo-tools/files/v1 |
| OS | macOS-27.0.1-arm64-arm-64bit-Mach-O |
| Architecture | arm64 |
| CPU | Apple M4 Max (arm64) |
| Physical cores | 16 |
| Logical threads | 16 |
| Memory (bytes) | 68719476736 |
| Tool: cargo | cargo 1.99.1 (5f94df478 2026-08-27) |
| Tool: just | just 1.59.0 |
| Tool: rustc | rustc 1.99.1 (58cf8f9c3 2026-10-13)<br>binary: rustc<br>commit-hash: 58cf8f9c35b355199727627d07ae34a4e2846d67<br>commit-date: 2026-10-13<br>host: aarch64-apple-darwin<br>release: 1.99.1<br>LLVM version: 23.1.1 |
| Dependency: criterion | 0.8.2 |
| Dependency: faer | 0.24.4 |
| Dependency: nalgebra | 0.35.0 |
| Environment overrides (names only) | CARGO&#95;ENCODED&#95;RUSTFLAGS, CARGO&#95;TARGET&#95;DIR, CRITERION&#95;HOME, RUSTC, RUSTC&#95;WORKSPACE&#95;WRAPPER, RUSTC&#95;WRAPPER, RUSTFLAGS, RUSTUP&#95;TOOLCHAIN |
| Measurement command (argv) | &#91;&quot;just&quot;, &quot;--justfile&quot;, &quot;tooling/performance.just&quot;, &quot;--working-directory&quot;, &quot;.&quot;, &quot;current-all&quot;&#93; |
| Independent correctness gate | passed |
| Gate command (argv) | &#91;&quot;just&quot;, &quot;--justfile&quot;, &quot;tooling/performance.just&quot;, &quot;--working-directory&quot;, &quot;.&quot;, &quot;gate&quot;&#93; |
