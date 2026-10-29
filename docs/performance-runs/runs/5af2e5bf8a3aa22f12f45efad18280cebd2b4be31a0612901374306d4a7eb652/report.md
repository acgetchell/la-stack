<!-- research-repo-tools/run-report/v2 -->
# Benchmark timings

Run: `5af2e5bf8a3aa22f12f45efad18280cebd2b4be31a0612901374306d4a7eb652`

Evidence: [run payload](run.json) · [evidence envelope](evidence.json).

Intervals are marginal timing intervals. No significance or ratio interval is inferred.

- la-stack baseline: measured in **baseline** phase (v0.5.1), revision `80a67d598b8f5788c2c8ffc12a0e8677c306909e`.
- la-stack current: measured in **current** phase (v0.5.2), revision `b5156e222add90fc62191db55dcccb40875a5e41`.
- nalgebra: measured in **baseline** phase (v0.5.1), revision `80a67d598b8f5788c2c8ffc12a0e8677c306909e`.
- faer: measured in **baseline** phase (v0.5.1), revision `80a67d598b8f5788c2c8ffc12a0e8677c306909e`.

## Mean (ns)

| Case | la-stack baseline | la-stack current | nalgebra | faer |
| --- | --- | --- | --- | --- |
| angle&#95;d3/dense/borrowed&#95;slice | 24.1045 ns [24.0501, 24.1614] (0.95 confidence) | 23.9402 ns [23.8958, 23.9893] (0.95 confidence) | — | — |
| angle&#95;d3/dense/construct&#95;then&#95;angle | 23.0417 ns [22.9929, 23.0959] (0.95 confidence) | 23.0252 ns [22.9932, 23.0579] (0.95 confidence) | — | — |
| angle&#95;d3/dense/prepared&#95;vector | 23.0082 ns [22.9623, 23.0559] (0.95 confidence) | 22.7546 ns [22.7206, 22.787] (0.95 confidence) | — | — |
| angle&#95;d3/dense&#95;antipodal/borrowed&#95;slice | 17.0516 ns [17.0204, 17.0846] (0.95 confidence) | 16.8696 ns [16.8404, 16.9011] (0.95 confidence) | — | — |
| angle&#95;d3/dense&#95;antipodal/construct&#95;then&#95;angle | 15.8358 ns [15.7998, 15.8775] (0.95 confidence) | 15.9013 ns [15.8198, 15.9932] (0.95 confidence) | — | — |
| angle&#95;d3/dense&#95;antipodal/prepared&#95;vector | 16.181 ns [16.1433, 16.2192] (0.95 confidence) | 15.8392 ns [15.8111, 15.8695] (0.95 confidence) | — | — |
| angle&#95;d3/dense&#95;identical/borrowed&#95;slice | 14.3098 ns [14.2929, 14.3277] (0.95 confidence) | 14.3596 ns [14.3462, 14.3736] (0.95 confidence) | — | — |
| angle&#95;d3/dense&#95;identical/construct&#95;then&#95;angle | 13.6055 ns [13.578, 13.634] (0.95 confidence) | 13.5407 ns [13.5254, 13.5568] (0.95 confidence) | — | — |
| angle&#95;d3/dense&#95;identical/prepared&#95;vector | 13.6832 ns [13.647, 13.7256] (0.95 confidence) | 13.5043 ns [13.4844, 13.5255] (0.95 confidence) | — | — |
| angle&#95;d3/dense&#95;varied/borrowed&#95;slice | 26.432 ns [26.3864, 26.4784] (0.95 confidence) | 26.1615 ns [26.1246, 26.1997] (0.95 confidence) | — | — |
| angle&#95;d3/dense&#95;varied/construct&#95;then&#95;angle | 25.5848 ns [25.5264, 25.6478] (0.95 confidence) | 25.173 ns [25.1194, 25.2296] (0.95 confidence) | — | — |
| angle&#95;d3/dense&#95;varied/prepared&#95;vector | 25.1523 ns [25.0805, 25.2391] (0.95 confidence) | 24.6357 ns [24.5963, 24.6749] (0.95 confidence) | — | — |
| angle&#95;d3/maximum&#95;small&#95;angle/borrowed&#95;slice | 10.8533 ns [10.8324, 10.8747] (0.95 confidence) | 10.9019 ns [10.8695, 10.9368] (0.95 confidence) | — | — |
| angle&#95;d3/maximum&#95;small&#95;angle/construct&#95;then&#95;angle | 10.4951 ns [10.4613, 10.5296] (0.95 confidence) | 10.0919 ns [10.0719, 10.1134] (0.95 confidence) | — | — |
| angle&#95;d3/maximum&#95;small&#95;angle/prepared&#95;vector | 9.60406 ns [9.5873, 9.62261] (0.95 confidence) | 9.50976 ns [9.49754, 9.52272] (0.95 confidence) | — | — |
| angle&#95;d3/mixed&#95;scale/borrowed&#95;slice | 12.2969 ns [12.278, 12.3176] (0.95 confidence) | 12.3259 ns [12.286, 12.3703] (0.95 confidence) | — | — |
| angle&#95;d3/mixed&#95;scale/construct&#95;then&#95;angle | 11.8417 ns [11.8227, 11.8609] (0.95 confidence) | 11.3521 ns [11.3252, 11.3801] (0.95 confidence) | — | — |
| angle&#95;d3/mixed&#95;scale/prepared&#95;vector | 10.7291 ns [10.7169, 10.7421] (0.95 confidence) | 10.6583 ns [10.6452, 10.6721] (0.95 confidence) | — | — |
| angle&#95;d3/near&#95;antipodal/borrowed&#95;slice | 16.5935 ns [16.5753, 16.6126] (0.95 confidence) | 16.4667 ns [16.4448, 16.4892] (0.95 confidence) | — | — |
| angle&#95;d3/near&#95;antipodal/construct&#95;then&#95;angle | 15.4403 ns [15.4218, 15.46] (0.95 confidence) | 15.2032 ns [15.1821, 15.2257] (0.95 confidence) | — | — |
| angle&#95;d3/near&#95;antipodal/prepared&#95;vector | 15.1386 ns [15.1229, 15.1548] (0.95 confidence) | 14.9188 ns [14.8966, 14.9435] (0.95 confidence) | — | — |
| angle&#95;d3/near&#95;parallel/borrowed&#95;slice | 10.8373 ns [10.8199, 10.8555] (0.95 confidence) | 10.8549 ns [10.8254, 10.8868] (0.95 confidence) | — | — |
| angle&#95;d3/near&#95;parallel/construct&#95;then&#95;angle | 10.4133 ns [10.3893, 10.4372] (0.95 confidence) | 10.135 ns [10.109, 10.1624] (0.95 confidence) | — | — |
| angle&#95;d3/near&#95;parallel/prepared&#95;vector | 9.6709 ns [9.65129, 9.691] (0.95 confidence) | 9.46647 ns [9.44972, 9.48356] (0.95 confidence) | — | — |
| angle&#95;d3/sparse&#95;adjacent/borrowed&#95;slice | 13.8735 ns [13.8449, 13.905] (0.95 confidence) | 13.6401 ns [13.5976, 13.6866] (0.95 confidence) | — | — |
| angle&#95;d3/sparse&#95;adjacent/construct&#95;then&#95;angle | 13.4274 ns [13.3917, 13.4683] (0.95 confidence) | 13.0165 ns [12.9977, 13.0361] (0.95 confidence) | — | — |
| angle&#95;d3/sparse&#95;adjacent/prepared&#95;vector | 12.0474 ns [12.0302, 12.0663] (0.95 confidence) | 11.7578 ns [11.7398, 11.7775] (0.95 confidence) | — | — |
| angle&#95;d3/sparse&#95;alternating/borrowed&#95;slice | 14.1659 ns [14.1227, 14.2145] (0.95 confidence) | 13.9385 ns [13.8958, 13.9837] (0.95 confidence) | — | — |
| angle&#95;d3/sparse&#95;alternating/construct&#95;then&#95;angle | 13.8047 ns [13.7309, 13.9274] (0.95 confidence) | 13.0442 ns [13.0221, 13.0683] (0.95 confidence) | — | — |
| angle&#95;d3/sparse&#95;alternating/prepared&#95;vector | 12.4264 ns [12.375, 12.5004] (0.95 confidence) | 12.0153 ns [11.9988, 12.0339] (0.95 confidence) | — | — |
| angle&#95;d3/sparse&#95;antipodal/borrowed&#95;slice | 11.2173 ns [11.1977, 11.2401] (0.95 confidence) | 10.8043 ns [10.7724, 10.8421] (0.95 confidence) | — | — |
| angle&#95;d3/sparse&#95;antipodal/construct&#95;then&#95;angle | 10.4847 ns [10.4582, 10.5126] (0.95 confidence) | 9.99354 ns [9.9732, 10.015] (0.95 confidence) | — | — |
| angle&#95;d3/sparse&#95;antipodal/prepared&#95;vector | 9.0665 ns [9.05273, 9.08216] (0.95 confidence) | 8.76192 ns [8.74995, 8.77497] (0.95 confidence) | — | — |
| angle&#95;d3/sparse&#95;identical/borrowed&#95;slice | 8.88425 ns [8.86004, 8.91047] (0.95 confidence) | 9.24627 ns [9.17461, 9.32039] (0.95 confidence) | — | — |
| angle&#95;d3/sparse&#95;identical/construct&#95;then&#95;angle | 8.35069 ns [8.33183, 8.37074] (0.95 confidence) | 8.41849 ns [8.3717, 8.47281] (0.95 confidence) | — | — |
| angle&#95;d3/sparse&#95;identical/prepared&#95;vector | 7.08777 ns [7.07204, 7.10549] (0.95 confidence) | 7.11441 ns [7.08871, 7.14033] (0.95 confidence) | — | — |
| angle&#95;d3/subnormal/borrowed&#95;slice | 10.8422 ns [10.8163, 10.8747] (0.95 confidence) | 10.9951 ns [10.9398, 11.0538] (0.95 confidence) | — | — |
| angle&#95;d3/subnormal/construct&#95;then&#95;angle | 10.3714 ns [10.3406, 10.4055] (0.95 confidence) | 10.1557 ns [10.1319, 10.18] (0.95 confidence) | — | — |
| angle&#95;d3/subnormal/prepared&#95;vector | 9.65915 ns [9.62169, 9.70342] (0.95 confidence) | 9.47147 ns [9.45516, 9.488] (0.95 confidence) | — | — |
| angle&#95;d3/subnormal&#95;axes/borrowed&#95;slice | 14.0469 ns [14.0331, 14.0617] (0.95 confidence) | 13.9781 ns [13.9239, 14.0392] (0.95 confidence) | — | — |
| angle&#95;d3/subnormal&#95;axes/construct&#95;then&#95;angle | 13.572 ns [13.5519, 13.5937] (0.95 confidence) | 13.0656 ns [13.0432, 13.0888] (0.95 confidence) | — | — |
| angle&#95;d3/subnormal&#95;axes/prepared&#95;vector | 12.2804 ns [12.2515, 12.3144] (0.95 confidence) | 12.0054 ns [11.9873, 12.027] (0.95 confidence) | — | — |
| angle&#95;d4/dense/borrowed&#95;slice | 34.6587 ns [34.6089, 34.7072] (0.95 confidence) | 34.3516 ns [34.2921, 34.4133] (0.95 confidence) | — | — |
| angle&#95;d4/dense/construct&#95;then&#95;angle | 33.2409 ns [33.1938, 33.2935] (0.95 confidence) | 32.9266 ns [32.874, 32.9828] (0.95 confidence) | — | — |
| angle&#95;d4/dense/prepared&#95;vector | 33.0183 ns [32.9826, 33.0568] (0.95 confidence) | 32.6268 ns [32.5669, 32.6851] (0.95 confidence) | — | — |
| angle&#95;d4/dense&#95;antipodal/borrowed&#95;slice | 26.8127 ns [26.7897, 26.8374] (0.95 confidence) | 26.7719 ns [26.7357, 26.8105] (0.95 confidence) | — | — |
| angle&#95;d4/dense&#95;antipodal/construct&#95;then&#95;angle | 25.1949 ns [25.1544, 25.2376] (0.95 confidence) | 25.1685 ns [25.1381, 25.2001] (0.95 confidence) | — | — |
| angle&#95;d4/dense&#95;antipodal/prepared&#95;vector | 25.0349 ns [25.0108, 25.0613] (0.95 confidence) | 24.9528 ns [24.9131, 24.9951] (0.95 confidence) | — | — |
| angle&#95;d4/dense&#95;identical/borrowed&#95;slice | 24.4717 ns [24.4196, 24.5298] (0.95 confidence) | 24.4284 ns [24.3722, 24.498] (0.95 confidence) | — | — |
| angle&#95;d4/dense&#95;identical/construct&#95;then&#95;angle | 22.7958 ns [22.7668, 22.8288] (0.95 confidence) | 22.7814 ns [22.7559, 22.8072] (0.95 confidence) | — | — |
| angle&#95;d4/dense&#95;identical/prepared&#95;vector | 22.5829 ns [22.5613, 22.6052] (0.95 confidence) | 22.2958 ns [22.2611, 22.3307] (0.95 confidence) | — | — |
| angle&#95;d4/dense&#95;varied/borrowed&#95;slice | 36.6672 ns [36.5288, 36.8334] (0.95 confidence) | 35.6879 ns [35.6237, 35.7534] (0.95 confidence) | — | — |
| angle&#95;d4/dense&#95;varied/construct&#95;then&#95;angle | 34.6009 ns [34.5603, 34.643] (0.95 confidence) | 34.0731 ns [34.0193, 34.1272] (0.95 confidence) | — | — |
| angle&#95;d4/dense&#95;varied/prepared&#95;vector | 34.2943 ns [34.1923, 34.397] (0.95 confidence) | 33.4258 ns [33.3666, 33.4842] (0.95 confidence) | — | — |
| angle&#95;d4/maximum&#95;small&#95;angle/borrowed&#95;slice | 12.2799 ns [12.2193, 12.3469] (0.95 confidence) | 12.5577 ns [12.488, 12.628] (0.95 confidence) | — | — |
| angle&#95;d4/maximum&#95;small&#95;angle/construct&#95;then&#95;angle | 11.3221 ns [11.2711, 11.3809] (0.95 confidence) | 11.4326 ns [11.3875, 11.4785] (0.95 confidence) | — | — |
| angle&#95;d4/maximum&#95;small&#95;angle/prepared&#95;vector | 10.0218 ns [9.99792, 10.0484] (0.95 confidence) | 9.99594 ns [9.96445, 10.0327] (0.95 confidence) | — | — |
| angle&#95;d4/mixed&#95;scale/borrowed&#95;slice | 13.7528 ns [13.7171, 13.791] (0.95 confidence) | 13.9773 ns [13.9067, 14.053] (0.95 confidence) | — | — |
| angle&#95;d4/mixed&#95;scale/construct&#95;then&#95;angle | 12.7731 ns [12.7242, 12.8264] (0.95 confidence) | 12.6081 ns [12.537, 12.6839] (0.95 confidence) | — | — |
| angle&#95;d4/mixed&#95;scale/prepared&#95;vector | 11.4763 ns [11.423, 11.5388] (0.95 confidence) | 11.1324 ns [11.1102, 11.1554] (0.95 confidence) | — | — |
| angle&#95;d4/near&#95;antipodal/borrowed&#95;slice | 17.8531 ns [17.8351, 17.8746] (0.95 confidence) | 17.7696 ns [17.7461, 17.7948] (0.95 confidence) | — | — |
| angle&#95;d4/near&#95;antipodal/construct&#95;then&#95;angle | 16.2484 ns [16.2287, 16.2699] (0.95 confidence) | 16.495 ns [16.2587, 16.916] (0.95 confidence) | — | — |
| angle&#95;d4/near&#95;antipodal/prepared&#95;vector | 15.8352 ns [15.8197, 15.8516] (0.95 confidence) | 15.5896 ns [15.5701, 15.6107] (0.95 confidence) | — | — |
| angle&#95;d4/near&#95;parallel/borrowed&#95;slice | 12.1235 ns [12.1021, 12.1461] (0.95 confidence) | 12.7498 ns [12.658, 12.8462] (0.95 confidence) | — | — |
| angle&#95;d4/near&#95;parallel/construct&#95;then&#95;angle | 11.1614 ns [11.1283, 11.1966] (0.95 confidence) | 11.2694 ns [11.2192, 11.3222] (0.95 confidence) | — | — |
| angle&#95;d4/near&#95;parallel/prepared&#95;vector | 10.0354 ns [10.0039, 10.0703] (0.95 confidence) | 9.97099 ns [9.94791, 9.99683] (0.95 confidence) | — | — |
| angle&#95;d4/sparse&#95;adjacent/borrowed&#95;slice | 15.7631 ns [15.7388, 15.7915] (0.95 confidence) | 15.7319 ns [15.6923, 15.7724] (0.95 confidence) | — | — |
| angle&#95;d4/sparse&#95;adjacent/construct&#95;then&#95;angle | 14.9177 ns [14.885, 14.9507] (0.95 confidence) | 14.723 ns [14.6854, 14.7612] (0.95 confidence) | — | — |
| angle&#95;d4/sparse&#95;adjacent/prepared&#95;vector | 13.506 ns [13.4942, 13.5186] (0.95 confidence) | 13.2752 ns [13.252, 13.3008] (0.95 confidence) | — | — |
| angle&#95;d4/sparse&#95;alternating/borrowed&#95;slice | 16.5045 ns [16.388, 16.6313] (0.95 confidence) | 15.9947 ns [15.9442, 16.0471] (0.95 confidence) | — | — |
| angle&#95;d4/sparse&#95;alternating/construct&#95;then&#95;angle | 15.1339 ns [15.097, 15.1718] (0.95 confidence) | 14.9851 ns [14.9408, 15.0311] (0.95 confidence) | — | — |
| angle&#95;d4/sparse&#95;alternating/prepared&#95;vector | 13.8278 ns [13.8022, 13.8546] (0.95 confidence) | 13.5873 ns [13.5601, 13.6172] (0.95 confidence) | — | — |
| angle&#95;d4/sparse&#95;antipodal/borrowed&#95;slice | 12.0208 ns [12.0066, 12.0374] (0.95 confidence) | 11.9661 ns [11.9205, 12.0122] (0.95 confidence) | — | — |
| angle&#95;d4/sparse&#95;antipodal/construct&#95;then&#95;angle | 11.138 ns [11.1053, 11.173] (0.95 confidence) | 10.8367 ns [10.8049, 10.8702] (0.95 confidence) | — | — |
| angle&#95;d4/sparse&#95;antipodal/prepared&#95;vector | 9.74935 ns [9.73519, 9.769] (0.95 confidence) | 9.49584 ns [9.4804, 9.51372] (0.95 confidence) | — | — |
| angle&#95;d4/sparse&#95;identical/borrowed&#95;slice | 10.279 ns [10.2573, 10.3041] (0.95 confidence) | 10.7151 ns [10.6455, 10.787] (0.95 confidence) | — | — |
| angle&#95;d4/sparse&#95;identical/construct&#95;then&#95;angle | 9.3929 ns [9.36473, 9.42185] (0.95 confidence) | 9.45617 ns [9.41779, 9.49576] (0.95 confidence) | — | — |
| angle&#95;d4/sparse&#95;identical/prepared&#95;vector | 8.07517 ns [8.04286, 8.11024] (0.95 confidence) | 7.99897 ns [7.98802, 8.01184] (0.95 confidence) | — | — |
| angle&#95;d4/subnormal/borrowed&#95;slice | 12.1725 ns [12.1436, 12.2034] (0.95 confidence) | 12.7121 ns [12.6435, 12.7832] (0.95 confidence) | — | — |
| angle&#95;d4/subnormal/construct&#95;then&#95;angle | 11.2063 ns [11.1622, 11.2528] (0.95 confidence) | 11.4972 ns [11.4609, 11.5332] (0.95 confidence) | — | — |
| angle&#95;d4/subnormal/prepared&#95;vector | 10.0101 ns [9.97626, 10.0517] (0.95 confidence) | 10.2786 ns [10.2223, 10.3367] (0.95 confidence) | — | — |
| angle&#95;d4/subnormal&#95;axes/borrowed&#95;slice | 15.5001 ns [15.4471, 15.5565] (0.95 confidence) | 15.4432 ns [15.3854, 15.5024] (0.95 confidence) | — | — |
| angle&#95;d4/subnormal&#95;axes/construct&#95;then&#95;angle | 14.5309 ns [14.4837, 14.5804] (0.95 confidence) | 14.2424 ns [14.1916, 14.2955] (0.95 confidence) | — | — |
| angle&#95;d4/subnormal&#95;axes/prepared&#95;vector | 13.3663 ns [13.3111, 13.4269] (0.95 confidence) | 12.8402 ns [12.8131, 12.8709] (0.95 confidence) | — | — |
| angle&#95;d5/dense/borrowed&#95;slice | 48.4748 ns [48.4078, 48.56] (0.95 confidence) | 48.9578 ns [48.863, 49.0654] (0.95 confidence) | — | — |
| angle&#95;d5/dense/construct&#95;then&#95;angle | 46.175 ns [46.1218, 46.2322] (0.95 confidence) | 46.3174 ns [46.2381, 46.4006] (0.95 confidence) | — | — |
| angle&#95;d5/dense/prepared&#95;vector | 45.7642 ns [45.7194, 45.8118] (0.95 confidence) | 45.9059 ns [45.8185, 46.008] (0.95 confidence) | — | — |
| angle&#95;d5/dense&#95;antipodal/borrowed&#95;slice | 37.277 ns [37.2257, 37.3312] (0.95 confidence) | 37.3871 ns [37.3407, 37.438] (0.95 confidence) | — | — |
| angle&#95;d5/dense&#95;antipodal/construct&#95;then&#95;angle | 35.2251 ns [35.1809, 35.2754] (0.95 confidence) | 35.2511 ns [35.2001, 35.3081] (0.95 confidence) | — | — |
| angle&#95;d5/dense&#95;antipodal/prepared&#95;vector | 35.5177 ns [35.4297, 35.6129] (0.95 confidence) | 35.181 ns [35.0776, 35.3032] (0.95 confidence) | — | — |
| angle&#95;d5/dense&#95;identical/borrowed&#95;slice | 35.867 ns [35.8353, 35.9016] (0.95 confidence) | 36.5627 ns [36.3416, 36.8001] (0.95 confidence) | — | — |
| angle&#95;d5/dense&#95;identical/construct&#95;then&#95;angle | 34.2926 ns [34.2277, 34.3673] (0.95 confidence) | 34.282 ns [34.2316, 34.3351] (0.95 confidence) | — | — |
| angle&#95;d5/dense&#95;identical/prepared&#95;vector | 33.602 ns [33.5424, 33.6698] (0.95 confidence) | 33.6337 ns [33.5797, 33.6842] (0.95 confidence) | — | — |
| angle&#95;d5/dense&#95;varied/borrowed&#95;slice | 70.6778 ns [70.5695, 70.7874] (0.95 confidence) | 70.6334 ns [70.4849, 70.7783] (0.95 confidence) | — | — |
| angle&#95;d5/dense&#95;varied/construct&#95;then&#95;angle | 68.1036 ns [67.9793, 68.2296] (0.95 confidence) | 67.9214 ns [67.813, 68.0311] (0.95 confidence) | — | — |
| angle&#95;d5/dense&#95;varied/prepared&#95;vector | 66.1363 ns [66.0195, 66.2627] (0.95 confidence) | 65.9772 ns [65.8368, 66.1262] (0.95 confidence) | — | — |
| angle&#95;d5/maximum&#95;small&#95;angle/borrowed&#95;slice | 13.9096 ns [13.8313, 13.9924] (0.95 confidence) | 14.6108 ns [14.5402, 14.6836] (0.95 confidence) | — | — |
| angle&#95;d5/maximum&#95;small&#95;angle/construct&#95;then&#95;angle | 13.3112 ns [13.2542, 13.3694] (0.95 confidence) | 13.2199 ns [13.1481, 13.2923] (0.95 confidence) | — | — |
| angle&#95;d5/maximum&#95;small&#95;angle/prepared&#95;vector | 11.1173 ns [11.0832, 11.1515] (0.95 confidence) | 11.1645 ns [11.1277, 11.2019] (0.95 confidence) | — | — |
| angle&#95;d5/mixed&#95;scale/borrowed&#95;slice | 15.2647 ns [15.2128, 15.3199] (0.95 confidence) | 15.79 ns [15.7458, 15.833] (0.95 confidence) | — | — |
| angle&#95;d5/mixed&#95;scale/construct&#95;then&#95;angle | 14.4863 ns [14.4386, 14.5318] (0.95 confidence) | 14.145 ns [14.0826, 14.2111] (0.95 confidence) | — | — |
| angle&#95;d5/mixed&#95;scale/prepared&#95;vector | 12.5965 ns [12.5666, 12.6266] (0.95 confidence) | 12.4114 ns [12.3832, 12.439] (0.95 confidence) | — | — |
| angle&#95;d5/near&#95;antipodal/borrowed&#95;slice | 19.6325 ns [19.609, 19.6583] (0.95 confidence) | 19.5123 ns [19.4326, 19.6191] (0.95 confidence) | — | — |
| angle&#95;d5/near&#95;antipodal/construct&#95;then&#95;angle | 17.5209 ns [17.4995, 17.5447] (0.95 confidence) | 17.0629 ns [17.0443, 17.0836] (0.95 confidence) | — | — |
| angle&#95;d5/near&#95;antipodal/prepared&#95;vector | 16.6939 ns [16.6635, 16.7269] (0.95 confidence) | 16.6996 ns [16.6508, 16.7579] (0.95 confidence) | — | — |
| angle&#95;d5/near&#95;parallel/borrowed&#95;slice | 13.7301 ns [13.6735, 13.7887] (0.95 confidence) | 14.6257 ns [14.5699, 14.6794] (0.95 confidence) | — | — |
| angle&#95;d5/near&#95;parallel/construct&#95;then&#95;angle | 13.1497 ns [13.0995, 13.1995] (0.95 confidence) | 12.8261 ns [12.7804, 12.8718] (0.95 confidence) | — | — |
| angle&#95;d5/near&#95;parallel/prepared&#95;vector | 11.2645 ns [11.2031, 11.3286] (0.95 confidence) | 11.2142 ns [11.1728, 11.2574] (0.95 confidence) | — | — |
| angle&#95;d5/sparse&#95;adjacent/borrowed&#95;slice | 17.5096 ns [17.418, 17.6102] (0.95 confidence) | 17.6729 ns [17.622, 17.729] (0.95 confidence) | — | — |
| angle&#95;d5/sparse&#95;adjacent/construct&#95;then&#95;angle | 16.8418 ns [16.7843, 16.9006] (0.95 confidence) | 16.3834 ns [16.3399, 16.4283] (0.95 confidence) | — | — |
| angle&#95;d5/sparse&#95;adjacent/prepared&#95;vector | 14.8607 ns [14.8307, 14.8948] (0.95 confidence) | 14.6244 ns [14.5913, 14.6593] (0.95 confidence) | — | — |
| angle&#95;d5/sparse&#95;alternating/borrowed&#95;slice | 28.2011 ns [28.1661, 28.2384] (0.95 confidence) | 27.7565 ns [27.7295, 27.7853] (0.95 confidence) | — | — |
| angle&#95;d5/sparse&#95;alternating/construct&#95;then&#95;angle | 26.7053 ns [26.6834, 26.7308] (0.95 confidence) | 26.7088 ns [26.6873, 26.7333] (0.95 confidence) | — | — |
| angle&#95;d5/sparse&#95;alternating/prepared&#95;vector | 26.4798 ns [26.4248, 26.5446] (0.95 confidence) | 26.2169 ns [26.1923, 26.2435] (0.95 confidence) | — | — |
| angle&#95;d5/sparse&#95;antipodal/borrowed&#95;slice | 13.0305 ns [13.0067, 13.0567] (0.95 confidence) | 13.2083 ns [13.1717, 13.2447] (0.95 confidence) | — | — |
| angle&#95;d5/sparse&#95;antipodal/construct&#95;then&#95;angle | 12.5805 ns [12.5424, 12.6175] (0.95 confidence) | 11.9564 ns [11.9182, 11.9955] (0.95 confidence) | — | — |
| angle&#95;d5/sparse&#95;antipodal/prepared&#95;vector | 10.5258 ns [10.5131, 10.5405] (0.95 confidence) | 10.2772 ns [10.2683, 10.2871] (0.95 confidence) | — | — |
| angle&#95;d5/sparse&#95;identical/borrowed&#95;slice | 11.5494 ns [11.5248, 11.5769] (0.95 confidence) | 12.0665 ns [12.0299, 12.1062] (0.95 confidence) | — | — |
| angle&#95;d5/sparse&#95;identical/construct&#95;then&#95;angle | 10.8967 ns [10.8538, 10.9415] (0.95 confidence) | 10.7036 ns [10.6632, 10.746] (0.95 confidence) | — | — |
| angle&#95;d5/sparse&#95;identical/prepared&#95;vector | 9.05904 ns [9.03597, 9.08555] (0.95 confidence) | 9.10516 ns [9.08377, 9.12868] (0.95 confidence) | — | — |
| angle&#95;d5/subnormal/borrowed&#95;slice | 13.6459 ns [13.5905, 13.7048] (0.95 confidence) | 14.5245 ns [14.4769, 14.5707] (0.95 confidence) | — | — |
| angle&#95;d5/subnormal/construct&#95;then&#95;angle | 13.05 ns [12.9898, 13.1102] (0.95 confidence) | 12.7104 ns [12.6692, 12.7541] (0.95 confidence) | — | — |
| angle&#95;d5/subnormal/prepared&#95;vector | 11.1296 ns [11.1011, 11.1575] (0.95 confidence) | 11.1162 ns [11.0649, 11.1701] (0.95 confidence) | — | — |
| angle&#95;d5/subnormal&#95;axes/borrowed&#95;slice | 17.7337 ns [17.495, 18.0624] (0.95 confidence) | 17.5241 ns [17.4492, 17.5964] (0.95 confidence) | — | — |
| angle&#95;d5/subnormal&#95;axes/construct&#95;then&#95;angle | 16.4393 ns [16.373, 16.5072] (0.95 confidence) | 15.9618 ns [15.8881, 16.0432] (0.95 confidence) | — | — |
| angle&#95;d5/subnormal&#95;axes/prepared&#95;vector | 14.3732 ns [14.3056, 14.447] (0.95 confidence) | 14.2549 ns [14.1749, 14.3469] (0.95 confidence) | — | — |
| angle&#95;d6/dense/borrowed&#95;slice | 62.7279 ns [62.6351, 62.8269] (0.95 confidence) | 64.0222 ns [62.9407, 65.8334] (0.95 confidence) | — | — |
| angle&#95;d6/dense/construct&#95;then&#95;angle | 60.1158 ns [60.0397, 60.1952] (0.95 confidence) | 60.4122 ns [60.1122, 60.8198] (0.95 confidence) | — | — |
| angle&#95;d6/dense/prepared&#95;vector | 59.7743 ns [59.6956, 59.8513] (0.95 confidence) | 59.7492 ns [59.655, 59.853] (0.95 confidence) | — | — |
| angle&#95;d6/dense&#95;antipodal/borrowed&#95;slice | 49.1961 ns [49.1446, 49.2457] (0.95 confidence) | 49.1832 ns [49.1253, 49.2495] (0.95 confidence) | — | — |
| angle&#95;d6/dense&#95;antipodal/construct&#95;then&#95;angle | 46.7109 ns [46.6518, 46.7733] (0.95 confidence) | 46.9354 ns [46.8033, 47.0769] (0.95 confidence) | — | — |
| angle&#95;d6/dense&#95;antipodal/prepared&#95;vector | 46.215 ns [46.1674, 46.2655] (0.95 confidence) | 46.2264 ns [46.1757, 46.2815] (0.95 confidence) | — | — |
| angle&#95;d6/dense&#95;identical/borrowed&#95;slice | 48.0533 ns [48.0063, 48.0996] (0.95 confidence) | 48.6885 ns [48.5877, 48.802] (0.95 confidence) | — | — |
| angle&#95;d6/dense&#95;identical/construct&#95;then&#95;angle | 45.5259 ns [45.4848, 45.5673] (0.95 confidence) | 46.1879 ns [46.0095, 46.4399] (0.95 confidence) | — | — |
| angle&#95;d6/dense&#95;identical/prepared&#95;vector | 45.0639 ns [44.9997, 45.135] (0.95 confidence) | 45.1379 ns [45.0988, 45.18] (0.95 confidence) | — | — |
| angle&#95;d6/dense&#95;varied/borrowed&#95;slice | 95.1667 ns [94.9292, 95.4244] (0.95 confidence) | 95.6313 ns [95.3863, 95.8835] (0.95 confidence) | — | — |
| angle&#95;d6/dense&#95;varied/construct&#95;then&#95;angle | 92.8007 ns [92.6025, 92.9998] (0.95 confidence) | 93.2442 ns [93.0675, 93.4199] (0.95 confidence) | — | — |
| angle&#95;d6/dense&#95;varied/prepared&#95;vector | 91.4998 ns [91.3383, 91.6539] (0.95 confidence) | 92.487 ns [92.1318, 92.8928] (0.95 confidence) | — | — |
| angle&#95;d6/maximum&#95;small&#95;angle/borrowed&#95;slice | 15.1556 ns [15.1195, 15.1903] (0.95 confidence) | 15.8449 ns [15.7859, 15.9002] (0.95 confidence) | — | — |
| angle&#95;d6/maximum&#95;small&#95;angle/construct&#95;then&#95;angle | 14.2136 ns [14.1633, 14.2639] (0.95 confidence) | 14.397 ns [14.353, 14.4423] (0.95 confidence) | — | — |
| angle&#95;d6/maximum&#95;small&#95;angle/prepared&#95;vector | 12.3136 ns [12.2727, 12.3576] (0.95 confidence) | 12.3154 ns [12.2843, 12.3501] (0.95 confidence) | — | — |
| angle&#95;d6/mixed&#95;scale/borrowed&#95;slice | 16.7797 ns [16.7354, 16.8261] (0.95 confidence) | 17.1853 ns [17.1264, 17.247] (0.95 confidence) | — | — |
| angle&#95;d6/mixed&#95;scale/construct&#95;then&#95;angle | 15.8454 ns [15.7823, 15.9097] (0.95 confidence) | 15.6388 ns [15.5978, 15.6795] (0.95 confidence) | — | — |
| angle&#95;d6/mixed&#95;scale/prepared&#95;vector | 13.8212 ns [13.7828, 13.8609] (0.95 confidence) | 13.5798 ns [13.5403, 13.6238] (0.95 confidence) | — | — |
| angle&#95;d6/near&#95;antipodal/borrowed&#95;slice | 21.1817 ns [21.1565, 21.2094] (0.95 confidence) | 21.0872 ns [21.053, 21.1249] (0.95 confidence) | — | — |
| angle&#95;d6/near&#95;antipodal/construct&#95;then&#95;angle | 18.7242 ns [18.6718, 18.7803] (0.95 confidence) | 18.6054 ns [18.5799, 18.6338] (0.95 confidence) | — | — |
| angle&#95;d6/near&#95;antipodal/prepared&#95;vector | 17.4224 ns [17.3877, 17.46] (0.95 confidence) | 17.2592 ns [17.2348, 17.2846] (0.95 confidence) | — | — |
| angle&#95;d6/near&#95;parallel/borrowed&#95;slice | 15.1873 ns [15.1282, 15.2505] (0.95 confidence) | 16.0966 ns [16.0066, 16.2169] (0.95 confidence) | — | — |
| angle&#95;d6/near&#95;parallel/construct&#95;then&#95;angle | 14.4749 ns [14.4093, 14.5413] (0.95 confidence) | 14.5092 ns [14.4636, 14.5559] (0.95 confidence) | — | — |
| angle&#95;d6/near&#95;parallel/prepared&#95;vector | 12.256 ns [12.21, 12.3053] (0.95 confidence) | 12.2875 ns [12.2537, 12.3225] (0.95 confidence) | — | — |
| angle&#95;d6/sparse&#95;adjacent/borrowed&#95;slice | 18.9542 ns [18.924, 18.9858] (0.95 confidence) | 19.2159 ns [19.1493, 19.289] (0.95 confidence) | — | — |
| angle&#95;d6/sparse&#95;adjacent/construct&#95;then&#95;angle | 17.9807 ns [17.9354, 18.0274] (0.95 confidence) | 17.9341 ns [17.8642, 18.0102] (0.95 confidence) | — | — |
| angle&#95;d6/sparse&#95;adjacent/prepared&#95;vector | 16.1153 ns [16.0867, 16.1475] (0.95 confidence) | 15.9217 ns [15.8777, 15.9679] (0.95 confidence) | — | — |
| angle&#95;d6/sparse&#95;alternating/borrowed&#95;slice | 29.4678 ns [29.4339, 29.5046] (0.95 confidence) | 29.5055 ns [29.4609, 29.5531] (0.95 confidence) | — | — |
| angle&#95;d6/sparse&#95;alternating/construct&#95;then&#95;angle | 27.6402 ns [27.609, 27.6766] (0.95 confidence) | 27.8166 ns [27.7736, 27.8673] (0.95 confidence) | — | — |
| angle&#95;d6/sparse&#95;alternating/prepared&#95;vector | 27.0643 ns [27.0253, 27.1052] (0.95 confidence) | 27.2811 ns [27.233, 27.3386] (0.95 confidence) | — | — |
| angle&#95;d6/sparse&#95;antipodal/borrowed&#95;slice | 13.9534 ns [13.9391, 13.9692] (0.95 confidence) | 14.2799 ns [14.2399, 14.3185] (0.95 confidence) | — | — |
| angle&#95;d6/sparse&#95;antipodal/construct&#95;then&#95;angle | 13.1905 ns [13.1617, 13.2215] (0.95 confidence) | 13.226 ns [13.186, 13.2681] (0.95 confidence) | — | — |
| angle&#95;d6/sparse&#95;antipodal/prepared&#95;vector | 11.2208 ns [11.2068, 11.2366] (0.95 confidence) | 11.0536 ns [11.0387, 11.0711] (0.95 confidence) | — | — |
| angle&#95;d6/sparse&#95;identical/borrowed&#95;slice | 12.7827 ns [12.7522, 12.8168] (0.95 confidence) | 13.3388 ns [13.2727, 13.4116] (0.95 confidence) | — | — |
| angle&#95;d6/sparse&#95;identical/construct&#95;then&#95;angle | 11.8593 ns [11.8239, 11.8956] (0.95 confidence) | 12.0685 ns [12.0122, 12.1295] (0.95 confidence) | — | — |
| angle&#95;d6/sparse&#95;identical/prepared&#95;vector | 10.0041 ns [9.98771, 10.0217] (0.95 confidence) | 10.0926 ns [10.0746, 10.112] (0.95 confidence) | — | — |
| angle&#95;d6/subnormal/borrowed&#95;slice | 15.3455 ns [15.2746, 15.4241] (0.95 confidence) | 16.147 ns [16.0508, 16.2465] (0.95 confidence) | — | — |
| angle&#95;d6/subnormal/construct&#95;then&#95;angle | 14.3842 ns [14.3075, 14.468] (0.95 confidence) | 14.555 ns [14.5134, 14.5972] (0.95 confidence) | — | — |
| angle&#95;d6/subnormal/prepared&#95;vector | 12.287 ns [12.2464, 12.3311] (0.95 confidence) | 12.3456 ns [12.3137, 12.381] (0.95 confidence) | — | — |
| angle&#95;d6/subnormal&#95;axes/borrowed&#95;slice | 18.4952 ns [18.4211, 18.5777] (0.95 confidence) | 19.2836 ns [19.2132, 19.355] (0.95 confidence) | — | — |
| angle&#95;d6/subnormal&#95;axes/construct&#95;then&#95;angle | 17.522 ns [17.4657, 17.5773] (0.95 confidence) | 17.6781 ns [17.6176, 17.7404] (0.95 confidence) | — | — |
| angle&#95;d6/subnormal&#95;axes/prepared&#95;vector | 15.4473 ns [15.4214, 15.4721] (0.95 confidence) | 15.4914 ns [15.4459, 15.5387] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d2/rounded&#95;array | 170.392 ns [169.564, 171.242] (0.95 confidence) | 166.05 ns [165.771, 166.336] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d2/rounded&#95;result | 166.98 ns [166.469, 167.497] (0.95 confidence) | 166.756 ns [166.347, 167.183] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d2/strict&#95;array | 169.299 ns [168.869, 169.749] (0.95 confidence) | 170.045 ns [169.592, 170.443] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d2/strict&#95;result | 168.623 ns [168.224, 169.017] (0.95 confidence) | 172.416 ns [172.008, 172.842] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d3/rounded&#95;array | 172.651 ns [172.215, 173.148] (0.95 confidence) | 169.422 ns [169.114, 169.767] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d3/rounded&#95;result | 171.013 ns [170.753, 171.289] (0.95 confidence) | 169.806 ns [169.485, 170.152] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d3/strict&#95;array | 176.471 ns [176.036, 176.936] (0.95 confidence) | 171.347 ns [170.894, 171.835] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d3/strict&#95;result | 173.182 ns [172.643, 173.743] (0.95 confidence) | 169.218 ns [168.987, 169.448] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d4/rounded&#95;array | 204.379 ns [203.607, 205.222] (0.95 confidence) | 177.205 ns [176.863, 177.556] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d4/rounded&#95;result | 203.991 ns [203.383, 204.596] (0.95 confidence) | 176.654 ns [176.347, 176.967] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d4/strict&#95;array | 199.789 ns [199.338, 200.244] (0.95 confidence) | 176.829 ns [176.22, 177.467] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d4/strict&#95;result | 200.616 ns [199.524, 201.665] (0.95 confidence) | 175.042 ns [174.634, 175.461] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d5/rounded&#95;array | 185.06 ns [184.389, 185.762] (0.95 confidence) | 181.257 ns [180.899, 181.648] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d5/rounded&#95;result | 182.5 ns [182.101, 182.916] (0.95 confidence) | 181.639 ns [181.355, 181.932] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d5/strict&#95;array | 180.438 ns [179.912, 180.977] (0.95 confidence) | 180.288 ns [179.829, 180.763] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d5/strict&#95;result | 178.689 ns [178.367, 179.037] (0.95 confidence) | 176.702 ns [176.287, 177.154] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d2/rounded&#95;array | 12.3012 ns [12.2843, 12.3203] (0.95 confidence) | 11.5436 ns [11.536, 11.5509] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d2/rounded&#95;result | 12.2427 ns [12.2308, 12.2562] (0.95 confidence) | 11.5466 ns [11.5381, 11.5546] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d2/strict&#95;array | 9.74414 ns [9.7168, 9.76839] (0.95 confidence) | 8.76123 ns [8.6843, 8.84099] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d2/strict&#95;result | 8.62866 ns [8.48955, 8.80032] (0.95 confidence) | 7.85678 ns [7.83213, 7.88373] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d3/rounded&#95;array | 18.0124 ns [17.9976, 18.0295] (0.95 confidence) | 17.8624 ns [17.8477, 17.8763] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d3/rounded&#95;result | 18.015 ns [18.0001, 18.0314] (0.95 confidence) | 17.8668 ns [17.85, 17.8833] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d3/strict&#95;array | 13.8707 ns [13.801, 13.937] (0.95 confidence) | 12.4616 ns [12.3684, 12.5629] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d3/strict&#95;result | 12.0871 ns [12.0178, 12.157] (0.95 confidence) | 11.4507 ns [11.3886, 11.5189] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d4/rounded&#95;array | 23.889 ns [23.8614, 23.9232] (0.95 confidence) | 23.7089 ns [23.6901, 23.7273] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d4/rounded&#95;result | 23.9005 ns [23.8715, 23.9332] (0.95 confidence) | 23.7171 ns [23.6989, 23.7376] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d4/strict&#95;array | 18.3098 ns [18.2712, 18.3499] (0.95 confidence) | 17.0195 ns [16.8423, 17.1994] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d4/strict&#95;result | 15.7371 ns [15.6386, 15.8418] (0.95 confidence) | 15.5784 ns [15.4634, 15.6965] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d5/rounded&#95;array | 29.158 ns [29.1437, 29.175] (0.95 confidence) | 28.4198 ns [28.4081, 28.4318] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d5/rounded&#95;result | 28.8497 ns [28.8262, 28.8769] (0.95 confidence) | 28.454 ns [28.4313, 28.4816] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d5/strict&#95;array | 22.3636 ns [22.1233, 22.5987] (0.95 confidence) | 21.8926 ns [21.635, 22.1452] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d5/strict&#95;result | 20.4076 ns [20.2294, 20.5847] (0.95 confidence) | 19.524 ns [19.3696, 19.6842] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d2/rounded&#95;array | 172.68 ns [172.374, 173.007] (0.95 confidence) | 170.079 ns [169.859, 170.277] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d2/rounded&#95;result | 172.37 ns [172.179, 172.574] (0.95 confidence) | 170.191 ns [169.957, 170.416] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d2/strict&#95;array | 14.1546 ns [14.1107, 14.2041] (0.95 confidence) | 13.8756 ns [13.7737, 13.9785] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d2/strict&#95;result | 13.1086 ns [13.0045, 13.2219] (0.95 confidence) | 12.4382 ns [12.3231, 12.5519] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d3/rounded&#95;array | 178.74 ns [178.559, 178.925] (0.95 confidence) | 176.663 ns [176.45, 176.875] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d3/rounded&#95;result | 178.726 ns [178.526, 178.939] (0.95 confidence) | 177.719 ns [177.501, 177.937] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d3/strict&#95;array | 18.7087 ns [18.6292, 18.7937] (0.95 confidence) | 17.5993 ns [17.5359, 17.6571] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d3/strict&#95;result | 17.6381 ns [17.5218, 17.7582] (0.95 confidence) | 16.5062 ns [16.4157, 16.5998] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d4/rounded&#95;array | 185.582 ns [185.33, 185.854] (0.95 confidence) | 183.496 ns [183.333, 183.657] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d4/rounded&#95;result | 185.061 ns [184.742, 185.394] (0.95 confidence) | 184.97 ns [184.732, 185.219] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d4/strict&#95;array | 23.4685 ns [23.3873, 23.5508] (0.95 confidence) | 22.0265 ns [21.9716, 22.0875] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d4/strict&#95;result | 22.0218 ns [21.9052, 22.1374] (0.95 confidence) | 20.4761 ns [20.3178, 20.644] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d5/rounded&#95;array | 190.62 ns [190.354, 190.896] (0.95 confidence) | 187.773 ns [187.517, 188.009] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d5/rounded&#95;result | 190.637 ns [190.408, 190.876] (0.95 confidence) | 188.528 ns [188.285, 188.777] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d5/strict&#95;array | 28.2002 ns [28.1104, 28.3034] (0.95 confidence) | 26.7858 ns [26.6273, 26.9596] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d5/strict&#95;result | 26.2401 ns [26.1054, 26.3766] (0.95 confidence) | 24.5318 ns [24.325, 24.7463] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d2/rounded&#95;array | 171.972 ns [171.701, 172.286] (0.95 confidence) | 169.767 ns [169.612, 169.923] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d2/rounded&#95;result | 172.608 ns [172.286, 172.963] (0.95 confidence) | 169.996 ns [169.748, 170.225] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d2/strict&#95;array | 13.4184 ns [13.3675, 13.474] (0.95 confidence) | 12.7836 ns [12.7441, 12.8251] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d2/strict&#95;result | 12.7797 ns [12.7195, 12.8458] (0.95 confidence) | 12.5816 ns [12.4918, 12.6713] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d3/rounded&#95;array | 178.693 ns [178.481, 178.907] (0.95 confidence) | 176.754 ns [176.543, 176.961] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d3/rounded&#95;result | 178.237 ns [177.98, 178.518] (0.95 confidence) | 176.717 ns [176.522, 176.916] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d3/strict&#95;array | 17.8458 ns [17.7732, 17.9217] (0.95 confidence) | 16.8605 ns [16.7777, 16.9466] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d3/strict&#95;result | 16.8692 ns [16.7866, 16.9522] (0.95 confidence) | 16.4516 ns [16.3834, 16.5201] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d4/rounded&#95;array | 195.116 ns [193.022, 197.269] (0.95 confidence) | 183.67 ns [183.472, 183.872] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d4/rounded&#95;result | 184.343 ns [184.149, 184.544] (0.95 confidence) | 183.912 ns [183.666, 184.167] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d4/strict&#95;array | 22.4088 ns [22.3081, 22.5127] (0.95 confidence) | 21.5053 ns [21.4166, 21.5939] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d4/strict&#95;result | 20.8778 ns [20.6988, 21.0763] (0.95 confidence) | 20.1175 ns [20.0179, 20.2193] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d5/rounded&#95;array | 189.539 ns [189.321, 189.769] (0.95 confidence) | 191.686 ns [190.706, 192.745] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d5/rounded&#95;result | 189.276 ns [189.055, 189.532] (0.95 confidence) | 188.439 ns [188.22, 188.656] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d5/strict&#95;array | 26.8756 ns [26.8046, 26.9529] (0.95 confidence) | 24.8754 ns [24.709, 25.0467] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d5/strict&#95;result | 25.0684 ns [24.9289, 25.2123] (0.95 confidence) | 23.4251 ns [23.3074, 23.5525] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d2/rounded&#95;array | 12.0458 ns [12.0334, 12.0599] (0.95 confidence) | 11.5707 ns [11.5626, 11.5784] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d2/rounded&#95;result | 12.3429 ns [12.2909, 12.408] (0.95 confidence) | 11.5789 ns [11.5647, 11.5935] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d2/strict&#95;array | 57.0819 ns [56.9498, 57.2171] (0.95 confidence) | 55.7445 ns [55.638, 55.8524] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d2/strict&#95;result | 12.1698 ns [12.1513, 12.1932] (0.95 confidence) | 11.7804 ns [11.7594, 11.8022] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d3/rounded&#95;array | 18.0527 ns [18.0384, 18.0686] (0.95 confidence) | 18.2314 ns [18.1344, 18.3361] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d3/rounded&#95;result | 18.0498 ns [18.0298, 18.0708] (0.95 confidence) | 17.9396 ns [17.9248, 17.9545] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d3/strict&#95;array | 60.9706 ns [60.8656, 61.0802] (0.95 confidence) | 61.166 ns [60.9928, 61.3382] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d3/strict&#95;result | 15.9575 ns [15.9192, 16.0012] (0.95 confidence) | 15.4527 ns [15.4062, 15.499] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d4/rounded&#95;array | 24.001 ns [23.9732, 24.033] (0.95 confidence) | 23.7818 ns [23.763, 23.7987] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d4/rounded&#95;result | 24.0047 ns [23.9817, 24.0309] (0.95 confidence) | 23.8104 ns [23.793, 23.8272] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d4/strict&#95;array | 64.0852 ns [63.9269, 64.2514] (0.95 confidence) | 63.626 ns [63.4298, 63.8309] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d4/strict&#95;result | 19.2955 ns [19.2423, 19.3484] (0.95 confidence) | 19.4171 ns [19.393, 19.4409] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d5/rounded&#95;array | 28.5346 ns [28.5149, 28.5559] (0.95 confidence) | 28.3576 ns [28.3392, 28.3744] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d5/rounded&#95;result | 29.098 ns [29.0858, 29.1121] (0.95 confidence) | 28.3483 ns [28.3275, 28.3695] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d5/strict&#95;array | 69.1011 ns [68.9445, 69.2618] (0.95 confidence) | 68.0842 ns [67.8709, 68.3042] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d5/strict&#95;result | 23.1101 ns [23.0622, 23.1622] (0.95 confidence) | 22.3496 ns [22.2707, 22.4333] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d2/rounded&#95;array | 167.591 ns [167.123, 168.085] (0.95 confidence) | 165.741 ns [165.567, 165.92] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d2/rounded&#95;result | 166.714 ns [166.426, 167.031] (0.95 confidence) | 166.616 ns [166.261, 166.981] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d2/strict&#95;array | 177.078 ns [176.654, 177.537] (0.95 confidence) | 175.091 ns [174.775, 175.415] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d2/strict&#95;result | 175.91 ns [175.495, 176.352] (0.95 confidence) | 174.732 ns [174.252, 175.214] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d3/rounded&#95;array | 174.442 ns [174.195, 174.693] (0.95 confidence) | 172.022 ns [171.706, 172.344] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d3/rounded&#95;result | 172.713 ns [172.432, 173.017] (0.95 confidence) | 172.397 ns [172.043, 172.76] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d3/strict&#95;array | 182.468 ns [182.05, 182.894] (0.95 confidence) | 181.108 ns [180.392, 181.853] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d3/strict&#95;result | 178.923 ns [178.523, 179.344] (0.95 confidence) | 180.705 ns [180.305, 181.109] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d4/rounded&#95;array | 184.503 ns [184.033, 184.983] (0.95 confidence) | 179.648 ns [179.356, 179.947] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d4/rounded&#95;result | 181.853 ns [181.514, 182.214] (0.95 confidence) | 178.557 ns [178.232, 178.893] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d4/strict&#95;array | 188.566 ns [188.047, 189.09] (0.95 confidence) | 186.052 ns [185.205, 186.951] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d4/strict&#95;result | 186.622 ns [186.086, 187.157] (0.95 confidence) | 182.64 ns [181.999, 183.298] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d5/rounded&#95;array | 185.054 ns [184.594, 185.55] (0.95 confidence) | 183.919 ns [183.628, 184.225] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d5/rounded&#95;result | 183.442 ns [183.253, 183.659] (0.95 confidence) | 184.806 ns [184.327, 185.289] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d5/strict&#95;array | 189.034 ns [188.533, 189.551] (0.95 confidence) | 188.771 ns [188.05, 189.541] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d5/strict&#95;result | 187.894 ns [187.206, 188.599] (0.95 confidence) | 183.409 ns [182.962, 183.899] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d2/rounded&#95;array | 199.016 ns [198.54, 199.524] (0.95 confidence) | 201.942 ns [201.335, 202.559] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d2/rounded&#95;result | 200.679 ns [200.114, 201.275] (0.95 confidence) | 202.54 ns [202.209, 202.878] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d2/strict&#95;array | 13484.6 ns [13462.3, 13506.5] (0.95 confidence) | 14062.6 ns [14034.7, 14090.3] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d2/strict&#95;result | 201.182 ns [200.908, 201.462] (0.95 confidence) | 205.838 ns [205.481, 206.201] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d3/rounded&#95;array | 205.242 ns [204.82, 205.686] (0.95 confidence) | 206.25 ns [205.895, 206.615] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d3/rounded&#95;result | 205.094 ns [204.789, 205.412] (0.95 confidence) | 208.063 ns [207.703, 208.436] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d3/strict&#95;array | 13417.5 ns [13394.3, 13440.5] (0.95 confidence) | 14066.8 ns [14046, 14088] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d3/strict&#95;result | 205.491 ns [205.152, 205.831] (0.95 confidence) | 208.616 ns [208.208, 209.032] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d4/rounded&#95;array | 211.17 ns [210.781, 211.59] (0.95 confidence) | 216.003 ns [215.567, 216.483] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d4/rounded&#95;result | 209.551 ns [209.219, 209.922] (0.95 confidence) | 216.945 ns [216.506, 217.395] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d4/strict&#95;array | 13426.2 ns [13402, 13451.3] (0.95 confidence) | 14103.3 ns [14082.7, 14125.1] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d4/strict&#95;result | 207.48 ns [207.035, 207.955] (0.95 confidence) | 209.534 ns [209.082, 209.994] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d5/rounded&#95;array | 216.062 ns [215.692, 216.469] (0.95 confidence) | 221.593 ns [221.201, 222.009] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d5/rounded&#95;result | 217.777 ns [217.484, 218.088] (0.95 confidence) | 221.686 ns [221.376, 222.001] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d5/strict&#95;array | 13464 ns [13443, 13485.2] (0.95 confidence) | 14137.3 ns [14117, 14157.8] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d5/strict&#95;result | 213.618 ns [213.165, 214.079] (0.95 confidence) | 215.465 ns [215.134, 215.815] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d2/rounded&#95;array | 160.601 ns [160.251, 160.96] (0.95 confidence) | 159.058 ns [158.376, 159.78] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d2/rounded&#95;result | 163.189 ns [162.735, 163.645] (0.95 confidence) | 161.182 ns [160.712, 161.657] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d2/strict&#95;array | 3843.41 ns [3826.97, 3859.65] (0.95 confidence) | 3870.75 ns [3848.27, 3893.33] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d2/strict&#95;result | 166.476 ns [166.06, 166.914] (0.95 confidence) | 171.895 ns [170.566, 173.382] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d3/rounded&#95;array | 169.235 ns [168.726, 169.764] (0.95 confidence) | 165.144 ns [164.654, 165.683] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d3/rounded&#95;result | 168.082 ns [167.543, 168.634] (0.95 confidence) | 163.41 ns [162.844, 164.025] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d3/strict&#95;array | 3807.96 ns [3789.27, 3826.58] (0.95 confidence) | 3766.27 ns [3747.21, 3785.32] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d3/strict&#95;result | 167.737 ns [167.198, 168.289] (0.95 confidence) | 163.482 ns [162.817, 164.172] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d4/rounded&#95;array | 174.731 ns [174.296, 175.163] (0.95 confidence) | 169.263 ns [168.806, 169.777] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d4/rounded&#95;result | 174.488 ns [174.064, 174.93] (0.95 confidence) | 177.837 ns [176.702, 179.003] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d4/strict&#95;array | 3774.96 ns [3760.76, 3789.51] (0.95 confidence) | 3882.32 ns [3862.33, 3901.05] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d4/strict&#95;result | 172.381 ns [171.846, 172.933] (0.95 confidence) | 175.108 ns [174.05, 176.135] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d5/rounded&#95;array | 179.669 ns [179.061, 180.314] (0.95 confidence) | 176.521 ns [175.917, 177.149] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d5/rounded&#95;result | 181.502 ns [180.926, 182.097] (0.95 confidence) | 179.443 ns [178.776, 180.082] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d5/strict&#95;array | 3841.81 ns [3826.02, 3857.45] (0.95 confidence) | 3805.97 ns [3782.68, 3828.55] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d5/strict&#95;result | 177.39 ns [176.838, 177.952] (0.95 confidence) | 170.897 ns [170.2, 171.654] (0.95 confidence) | — | — |
| column&#95;basis/2x2/full&#95;rank | — | 4.88257 ns [4.87257, 4.89253] (0.95 confidence) | — | — |
| column&#95;basis/2x2/large | — | 4.8849 ns [4.87736, 4.89223] (0.95 confidence) | — | — |
| column&#95;basis/2x2/near&#95;dependent | — | 4.92458 ns [4.91475, 4.93428] (0.95 confidence) | — | — |
| column&#95;basis/2x2/rank&#95;one | — | 4.6824 ns [4.674, 4.69084] (0.95 confidence) | — | — |
| column&#95;basis/2x2/small | — | 4.90264 ns [4.89366, 4.91177] (0.95 confidence) | — | — |
| column&#95;basis/2x3/shared&#95;face | — | 9.95374 ns [9.92415, 9.9861] (0.95 confidence) | — | — |
| column&#95;basis/3x3/full&#95;rank | — | 15.9834 ns [15.2334, 16.7484] (0.95 confidence) | — | — |
| column&#95;basis/3x3/large | — | 16.0139 ns [15.2628, 16.7724] (0.95 confidence) | — | — |
| column&#95;basis/3x3/near&#95;dependent | — | 16.0682 ns [15.3238, 16.8259] (0.95 confidence) | — | — |
| column&#95;basis/3x3/rank&#95;one | — | 9.08962 ns [9.02781, 9.14775] (0.95 confidence) | — | — |
| column&#95;basis/3x3/small | — | 16.0336 ns [15.2952, 16.7893] (0.95 confidence) | — | — |
| column&#95;basis/3x4/shared&#95;face | — | 17.5518 ns [17.4345, 17.7641] (0.95 confidence) | — | — |
| column&#95;basis/4x4/full&#95;rank | — | 28.9656 ns [28.8499, 29.0768] (0.95 confidence) | — | — |
| column&#95;basis/4x4/large | — | 29.1335 ns [29.0159, 29.2485] (0.95 confidence) | — | — |
| column&#95;basis/4x4/near&#95;dependent | — | 29.316 ns [29.2027, 29.425] (0.95 confidence) | — | — |
| column&#95;basis/4x4/rank&#95;one | — | 14.5866 ns [14.5642, 14.6089] (0.95 confidence) | — | — |
| column&#95;basis/4x4/small | — | 29.008 ns [28.8901, 29.122] (0.95 confidence) | — | — |
| column&#95;basis/4x5/shared&#95;face | — | 39.0366 ns [38.4126, 39.7215] (0.95 confidence) | — | — |
| column&#95;basis/4x6/full&#95;rank | — | 36.1218 ns [36.0257, 36.2419] (0.95 confidence) | — | — |
| column&#95;basis/4x6/large | — | 36.0076 ns [35.9194, 36.1008] (0.95 confidence) | — | — |
| column&#95;basis/4x6/near&#95;dependent | — | 35.8012 ns [35.7464, 35.8593] (0.95 confidence) | — | — |
| column&#95;basis/4x6/rank&#95;one | — | 20.841 ns [20.7714, 20.912] (0.95 confidence) | — | — |
| column&#95;basis/4x6/small | — | 36.0106 ns [35.9206, 36.1033] (0.95 confidence) | — | — |
| column&#95;basis/5x2/full&#95;rank | — | 13.2202 ns [13.1959, 13.2428] (0.95 confidence) | — | — |
| column&#95;basis/5x2/large | — | 13.1889 ns [13.1534, 13.2173] (0.95 confidence) | — | — |
| column&#95;basis/5x2/near&#95;dependent | — | 13.1491 ns [13.1179, 13.1808] (0.95 confidence) | — | — |
| column&#95;basis/5x2/rank&#95;one | — | 11.0112 ns [10.9824, 11.0391] (0.95 confidence) | — | — |
| column&#95;basis/5x2/small | — | 13.1971 ns [13.1679, 13.2229] (0.95 confidence) | — | — |
| column&#95;basis/5x5/full&#95;rank | — | 52.8605 ns [52.1535, 53.6254] (0.95 confidence) | — | — |
| column&#95;basis/5x5/large | — | 52.7295 ns [52.0859, 53.4317] (0.95 confidence) | — | — |
| column&#95;basis/5x5/near&#95;dependent | — | 52.2591 ns [51.6788, 52.8921] (0.95 confidence) | — | — |
| column&#95;basis/5x5/rank&#95;one | — | 23.7991 ns [23.7555, 23.8434] (0.95 confidence) | — | — |
| column&#95;basis/5x5/small | — | 52.6512 ns [52.0101, 53.3537] (0.95 confidence) | — | — |
| column&#95;basis/5x6/shared&#95;face | — | 62.548 ns [62.0449, 63.1651] (0.95 confidence) | — | — |
| column&#95;basis/5x8/full&#95;rank | — | 61.8334 ns [61.7605, 61.9023] (0.95 confidence) | — | — |
| column&#95;basis/5x8/large | — | 61.8893 ns [61.8052, 61.9762] (0.95 confidence) | — | — |
| column&#95;basis/5x8/near&#95;dependent | — | 61.7574 ns [61.6148, 61.8784] (0.95 confidence) | — | — |
| column&#95;basis/5x8/rank&#95;one | — | 37.3205 ns [37.0004, 37.6547] (0.95 confidence) | — | — |
| column&#95;basis/5x8/small | — | 61.6267 ns [61.5526, 61.6967] (0.95 confidence) | — | — |
| column&#95;basis/6x10/full&#95;rank | — | 97.781 ns [97.6578, 97.8871] (0.95 confidence) | — | — |
| column&#95;basis/6x10/large | — | 97.6072 ns [97.4777, 97.724] (0.95 confidence) | — | — |
| column&#95;basis/6x10/near&#95;dependent | — | 97.0417 ns [96.9517, 97.1314] (0.95 confidence) | — | — |
| column&#95;basis/6x10/rank&#95;one | — | 56.2444 ns [55.9439, 56.5895] (0.95 confidence) | — | — |
| column&#95;basis/6x10/small | — | 97.7249 ns [97.6385, 97.8097] (0.95 confidence) | — | — |
| column&#95;basis/6x7/shared&#95;face | — | 97.9685 ns [97.5884, 98.4143] (0.95 confidence) | — | — |
| column&#95;basis/7x12/full&#95;rank | — | 142.523 ns [142.381, 142.66] (0.95 confidence) | — | — |
| column&#95;basis/7x12/large | — | 142.564 ns [142.424, 142.694] (0.95 confidence) | — | — |
| column&#95;basis/7x12/near&#95;dependent | — | 142.947 ns [142.811, 143.08] (0.95 confidence) | — | — |
| column&#95;basis/7x12/rank&#95;one | — | 78.8884 ns [78.771, 79.0115] (0.95 confidence) | — | — |
| column&#95;basis/7x12/small | — | 142.774 ns [142.653, 142.893] (0.95 confidence) | — | — |
| column&#95;basis/8x14/full&#95;rank | — | 197.252 ns [196.992, 197.522] (0.95 confidence) | — | — |
| column&#95;basis/8x14/large | — | 197.452 ns [197.215, 197.7] (0.95 confidence) | — | — |
| column&#95;basis/8x14/near&#95;dependent | — | 196.834 ns [196.592, 197.075] (0.95 confidence) | — | — |
| column&#95;basis/8x14/rank&#95;one | — | 103.784 ns [103.667, 103.898] (0.95 confidence) | — | — |
| column&#95;basis/8x14/small | — | 196.966 ns [196.724, 197.207] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;det | 401.505 ns [399.452, 403.579] (0.95 confidence) | 394.579 ns [392.711, 396.419] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;det&#95;from&#95;ldlt | 2.92718 ns [2.92269, 2.93175] (0.95 confidence) | 2.92441 ns [2.91848, 2.93117] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;det&#95;from&#95;lu | 3.32316 ns [3.3178, 3.32899] (0.95 confidence) | 3.29999 ns [3.29483, 3.30502] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;det&#95;via&#95;lu | 383.617 ns [381.482, 385.79] (0.95 confidence) | 429.269 ns [421.759, 437.126] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;dot | 2.36948 ns [2.36346, 2.37626] (0.95 confidence) | 2.34403 ns [2.34046, 2.34781] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;inf&#95;norm | 33.2143 ns [33.1727, 33.2579] (0.95 confidence) | 32.675 ns [32.6292, 32.7217] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;ldlt | 341.957 ns [341.595, 342.333] (0.95 confidence) | 344.624 ns [344.2, 345.061] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;ldlt&#95;solve | 384.48 ns [384.042, 384.952] (0.95 confidence) | 386.434 ns [385.903, 387.009] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;lu | 365.047 ns [363.793, 366.323] (0.95 confidence) | 427.067 ns [418.824, 435.46] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;lu&#95;solve | 632.693 ns [630.115, 634.955] (0.95 confidence) | 651.943 ns [644.193, 660.122] (0.95 confidence) | 590.01 ns [588.792, 591.221] (0.95 confidence) | 897.202 ns [894.001, 900.335] (0.95 confidence) |
| d16/la&#95;stack&#95;norm2 | 10.5313 ns [10.4793, 10.5842] (0.95 confidence) | 10.266 ns [10.2241, 10.3084] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;norm2&#95;scenario&#95;descending | 9.66632 ns [9.64182, 9.69241] (0.95 confidence) | 9.47427 ns [9.45376, 9.4948] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale | 9.62115 ns [9.60033, 9.64204] (0.95 confidence) | 9.51531 ns [9.48855, 9.54352] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;norm2&#95;scenario&#95;sparse | 5.01572 ns [4.9672, 5.06417] (0.95 confidence) | 4.58404 ns [4.52182, 4.65206] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range | 8.59527 ns [8.577, 8.61565] (0.95 confidence) | 8.33755 ns [8.32246, 8.35245] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;norm2&#95;sq | 2.13962 ns [2.13449, 2.1453] (0.95 confidence) | 2.10753 ns [2.10166, 2.11299] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;solve&#95;from&#95;ldlt | 28.1731 ns [28.1315, 28.2201] (0.95 confidence) | 28.0484 ns [28.0031, 28.0987] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;solve&#95;from&#95;lu | 195.209 ns [194.39, 196.099] (0.95 confidence) | 192.433 ns [192.017, 192.862] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;det | 0.577927 ns [0.577342, 0.578536] (0.95 confidence) | 0.578185 ns [0.577421, 0.579024] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;det&#95;from&#95;ldlt | 0.50415 ns [0.502261, 0.506367] (0.95 confidence) | 0.505999 ns [0.503212, 0.50904] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;det&#95;from&#95;lu | 0.553248 ns [0.552166, 0.554379] (0.95 confidence) | 0.524862 ns [0.524214, 0.525551] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;det&#95;via&#95;lu | 1.77119 ns [1.767, 1.77524] (0.95 confidence) | 1.60827 ns [1.60187, 1.61553] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;dot | 0.619244 ns [0.614712, 0.624991] (0.95 confidence) | 0.620084 ns [0.616312, 0.624781] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;inf&#95;norm | 0.710091 ns [0.709148, 0.711038] (0.95 confidence) | 0.709184 ns [0.708359, 0.710046] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;ldlt | 6.98083 ns [6.943, 7.02667] (0.95 confidence) | 6.60056 ns [6.58097, 6.62109] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;ldlt&#95;solve | 9.89485 ns [9.87043, 9.91901] (0.95 confidence) | 9.5078 ns [9.47976, 9.53841] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;lu | 1.56934 ns [1.56575, 1.57296] (0.95 confidence) | 1.3956 ns [1.39147, 1.3996] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;lu&#95;solve | 2.06392 ns [2.05763, 2.07044] (0.95 confidence) | 2.10425 ns [2.08425, 2.13147] (0.95 confidence) | 4.52904 ns [4.51928, 4.53913] (0.95 confidence) | 128.199 ns [127.796, 128.61] (0.95 confidence) |
| d2/la&#95;stack&#95;norm2 | 1.24505 ns [1.23841, 1.25169] (0.95 confidence) | 1.16686 ns [1.16283, 1.17139] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;norm2&#95;scenario&#95;descending | 1.36443 ns [1.36166, 1.36717] (0.95 confidence) | 1.35423 ns [1.35219, 1.35617] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale | 1.36518 ns [1.36178, 1.36923] (0.95 confidence) | 1.35592 ns [1.35324, 1.35866] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;norm2&#95;scenario&#95;sparse | 1.18063 ns [1.17447, 1.18716] (0.95 confidence) | 1.19333 ns [1.18548, 1.20115] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range | 1.36611 ns [1.36338, 1.36876] (0.95 confidence) | 1.35356 ns [1.3507, 1.35681] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;norm2&#95;sq | 0.404569 ns [0.403328, 0.405844] (0.95 confidence) | 0.396026 ns [0.395117, 0.397012] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;solve&#95;from&#95;ldlt | 1.24497 ns [1.24094, 1.24884] (0.95 confidence) | 1.23694 ns [1.23242, 1.24103] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;solve&#95;from&#95;lu | 1.28825 ns [1.28404, 1.29266] (0.95 confidence) | 1.28538 ns [1.28178, 1.28905] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;det | 1.36856 ns [1.36193, 1.37549] (0.95 confidence) | 1.36514 ns [1.35856, 1.372] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;det&#95;from&#95;ldlt | 0.597733 ns [0.597247, 0.598286] (0.95 confidence) | 0.600638 ns [0.600084, 0.601227] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;det&#95;from&#95;lu | 0.739396 ns [0.713635, 0.786969] (0.95 confidence) | 0.700524 ns [0.699824, 0.701249] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;det&#95;via&#95;lu | 8.86495 ns [8.85156, 8.87879] (0.95 confidence) | 8.78666 ns [8.77033, 8.80327] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;dot | 0.674627 ns [0.672942, 0.676436] (0.95 confidence) | 0.673466 ns [0.672473, 0.674482] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;inf&#95;norm | 1.37541 ns [1.37353, 1.37733] (0.95 confidence) | 1.38582 ns [1.38154, 1.39039] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;ldlt | 13.1977 ns [13.1527, 13.242] (0.95 confidence) | 8.58149 ns [8.57069, 8.59251] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;ldlt&#95;solve | 9.70153 ns [9.67131, 9.73329] (0.95 confidence) | 9.51881 ns [9.45754, 9.58302] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;lu | 8.45095 ns [8.43597, 8.46475] (0.95 confidence) | 8.41833 ns [8.3997, 8.43714] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;lu&#95;solve | 9.91313 ns [9.89443, 9.9341] (0.95 confidence) | 9.97148 ns [9.95479, 9.98806] (0.95 confidence) | 32.7352 ns [32.5918, 32.8851] (0.95 confidence) | 156.859 ns [156.399, 157.352] (0.95 confidence) |
| d3/la&#95;stack&#95;norm2 | 1.99895 ns [1.95259, 2.04674] (0.95 confidence) | 1.86497 ns [1.83002, 1.9035] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;norm2&#95;scenario&#95;descending | 1.8569 ns [1.8534, 1.86057] (0.95 confidence) | 1.84839 ns [1.84455, 1.85186] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale | 1.85007 ns [1.84511, 1.85474] (0.95 confidence) | 1.85177 ns [1.84797, 1.85581] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;norm2&#95;scenario&#95;sparse | 1.43959 ns [1.42994, 1.44919] (0.95 confidence) | 1.44272 ns [1.43273, 1.45285] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range | 1.85159 ns [1.84773, 1.85546] (0.95 confidence) | 1.85149 ns [1.84643, 1.85701] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;norm2&#95;sq | 0.446082 ns [0.444676, 0.447574] (0.95 confidence) | 0.441906 ns [0.441178, 0.442667] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;solve&#95;from&#95;ldlt | 1.80032 ns [1.79793, 1.80278] (0.95 confidence) | 1.8038 ns [1.80048, 1.80736] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;solve&#95;from&#95;lu | 2.09024 ns [2.08608, 2.09449] (0.95 confidence) | 2.0826 ns [2.07708, 2.08899] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;det | 2223.47 ns [2218.19, 2229.2] (0.95 confidence) | 2039.47 ns [2033.86, 2045.21] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;det&#95;from&#95;ldlt | 8.08993 ns [8.07215, 8.10934] (0.95 confidence) | 7.84333 ns [7.829, 7.85763] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;det&#95;from&#95;lu | 9.40397 ns [9.37389, 9.43531] (0.95 confidence) | 9.20567 ns [9.18108, 9.23167] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;det&#95;via&#95;lu | 2077.69 ns [2074.16, 2081.27] (0.95 confidence) | 2062.43 ns [2059.18, 2065.84] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;dot | 4.10092 ns [4.09644, 4.10586] (0.95 confidence) | 4.07495 ns [4.06487, 4.08563] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;inf&#95;norm | 131.587 ns [131.406, 131.777] (0.95 confidence) | 129.303 ns [129.086, 129.518] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;ldlt | 2436.46 ns [2431.91, 2442.22] (0.95 confidence) | 2417.13 ns [2413.88, 2420.6] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;ldlt&#95;solve | 2773.58 ns [2768.21, 2779.39] (0.95 confidence) | 2742.1 ns [2738, 2746.43] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;lu | 2306.51 ns [2292.16, 2320.56] (0.95 confidence) | 2280.71 ns [2267.53, 2293.66] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;lu&#95;solve | 2849.33 ns [2841.53, 2857.4] (0.95 confidence) | 2773.62 ns [2769.08, 2778.25] (0.95 confidence) | 2505.09 ns [2498.04, 2513.43] (0.95 confidence) | 2998.54 ns [2994.31, 3002.69] (0.95 confidence) |
| d32/la&#95;stack&#95;norm2 | 24.5641 ns [24.5054, 24.6219] (0.95 confidence) | 23.4958 ns [23.4403, 23.5489] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;norm2&#95;scenario&#95;descending | 19.8339 ns [19.7921, 19.876] (0.95 confidence) | 18.7457 ns [18.7116, 18.7786] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale | 19.5843 ns [19.5298, 19.6378] (0.95 confidence) | 18.7721 ns [18.7143, 18.8406] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;norm2&#95;scenario&#95;sparse | 9.83985 ns [9.78822, 9.88459] (0.95 confidence) | 9.31095 ns [9.28395, 9.33881] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range | 17.1601 ns [17.1275, 17.1928] (0.95 confidence) | 16.3343 ns [16.31, 16.3583] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;norm2&#95;sq | 4.08044 ns [4.07242, 4.08869] (0.95 confidence) | 4.01813 ns [4.01203, 4.02462] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;solve&#95;from&#95;ldlt | 320.007 ns [319.183, 320.923] (0.95 confidence) | 314.988 ns [314.404, 315.605] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;solve&#95;from&#95;lu | 718.721 ns [717.053, 720.513] (0.95 confidence) | 673.44 ns [671.973, 675.055] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;det | 2.57282 ns [2.56652, 2.57904] (0.95 confidence) | 2.54882 ns [2.54299, 2.55464] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;det&#95;from&#95;ldlt | 0.776064 ns [0.77502, 0.777213] (0.95 confidence) | 0.772655 ns [0.771977, 0.77344] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;det&#95;from&#95;lu | 0.879167 ns [0.877982, 0.880412] (0.95 confidence) | 0.871302 ns [0.870688, 0.871978] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;det&#95;via&#95;lu | 14.2874 ns [14.255, 14.3306] (0.95 confidence) | 14.7155 ns [14.6838, 14.7469] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;dot | 0.693512 ns [0.690822, 0.696353] (0.95 confidence) | 0.682345 ns [0.68139, 0.683357] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;inf&#95;norm | 2.36304 ns [2.35969, 2.36641] (0.95 confidence) | 2.36342 ns [2.361, 2.36595] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;ldlt | 17.9412 ns [17.8992, 17.9845] (0.95 confidence) | 18.1312 ns [18.098, 18.1641] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;ldlt&#95;solve | 18.0398 ns [18.0027, 18.0798] (0.95 confidence) | 16.6332 ns [16.6107, 16.6562] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;lu | 13.8586 ns [13.8372, 13.8801] (0.95 confidence) | 14.3499 ns [14.3266, 14.3736] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;lu&#95;solve | 21.9562 ns [21.8992, 22.0168] (0.95 confidence) | 23.734 ns [23.6968, 23.7747] (0.95 confidence) | 52.8804 ns [52.7601, 53.0026] (0.95 confidence) | 183.201 ns [182.744, 183.667] (0.95 confidence) |
| d4/la&#95;stack&#95;norm2 | 2.29123 ns [2.26814, 2.31925] (0.95 confidence) | 2.30625 ns [2.27492, 2.34171] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;norm2&#95;scenario&#95;descending | 2.40647 ns [2.40017, 2.41256] (0.95 confidence) | 2.40681 ns [2.40227, 2.41146] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale | 2.40055 ns [2.3948, 2.40576] (0.95 confidence) | 2.39859 ns [2.3949, 2.40236] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;norm2&#95;scenario&#95;sparse | 1.59906 ns [1.55604, 1.64632] (0.95 confidence) | 1.61345 ns [1.56977, 1.66089] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range | 2.16765 ns [2.16359, 2.17177] (0.95 confidence) | 2.15988 ns [2.1557, 2.16421] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;norm2&#95;sq | 0.458678 ns [0.457928, 0.459462] (0.95 confidence) | 0.465777 ns [0.464161, 0.467538] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;solve&#95;from&#95;ldlt | 2.53042 ns [2.522, 2.53914] (0.95 confidence) | 2.49087 ns [2.48456, 2.49793] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;solve&#95;from&#95;lu | 3.97217 ns [3.96361, 3.98091] (0.95 confidence) | 4.03493 ns [4.03125, 4.03899] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;det | 41.1572 ns [41.0685, 41.2468] (0.95 confidence) | 28.882 ns [28.8445, 28.9194] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;det&#95;from&#95;ldlt | 0.883946 ns [0.882646, 0.885304] (0.95 confidence) | 0.880258 ns [0.879597, 0.880984] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;det&#95;from&#95;lu | 1.0409 ns [1.03936, 1.04263] (0.95 confidence) | 1.03794 ns [1.03715, 1.03879] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;det&#95;via&#95;lu | 34.0376 ns [33.9658, 34.109] (0.95 confidence) | 34.4327 ns [34.3697, 34.4919] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;dot | 0.819936 ns [0.815995, 0.823916] (0.95 confidence) | 0.842805 ns [0.838781, 0.846159] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;inf&#95;norm | 3.48409 ns [3.47912, 3.48888] (0.95 confidence) | 3.48483 ns [3.48051, 3.48952] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;ldlt | 25.1194 ns [25.0745, 25.1668] (0.95 confidence) | 24.5301 ns [24.4968, 24.565] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;ldlt&#95;solve | 30.7216 ns [30.6584, 30.7866] (0.95 confidence) | 30.3116 ns [30.2638, 30.3615] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;lu | 33.844 ns [33.7557, 33.9284] (0.95 confidence) | 34.2245 ns [34.1584, 34.288] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;lu&#95;solve | 49.0016 ns [48.9399, 49.0688] (0.95 confidence) | 47.1819 ns [47.041, 47.3268] (0.95 confidence) | 74.7814 ns [74.4543, 75.129] (0.95 confidence) | 251.145 ns [250.673, 251.685] (0.95 confidence) |
| d5/la&#95;stack&#95;norm2 | 2.92892 ns [2.90265, 2.95579] (0.95 confidence) | 2.74456 ns [2.7231, 2.77096] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;norm2&#95;scenario&#95;descending | 3.08162 ns [3.07537, 3.08766] (0.95 confidence) | 3.18555 ns [3.18058, 3.19074] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale | 3.09913 ns [3.09166, 3.10779] (0.95 confidence) | 3.17466 ns [3.1674, 3.18131] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;norm2&#95;scenario&#95;sparse | 1.75682 ns [1.74233, 1.77677] (0.95 confidence) | 1.7313 ns [1.72825, 1.73441] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range | 2.69057 ns [2.68019, 2.70087] (0.95 confidence) | 3.15595 ns [3.14949, 3.16244] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;norm2&#95;sq | 0.530002 ns [0.529164, 0.530923] (0.95 confidence) | 0.528587 ns [0.527824, 0.529399] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;solve&#95;from&#95;ldlt | 3.90445 ns [3.89516, 3.91374] (0.95 confidence) | 3.9733 ns [3.96428, 3.98234] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;solve&#95;from&#95;lu | 6.18384 ns [6.17035, 6.19785] (0.95 confidence) | 6.09741 ns [6.09002, 6.10571] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;det | 14850.7 ns [14824.9, 14877.1] (0.95 confidence) | 14387.2 ns [14367.5, 14407.4] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;det&#95;from&#95;ldlt | 23.5372 ns [23.4912, 23.5848] (0.95 confidence) | 23.1635 ns [23.1068, 23.2165] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;det&#95;from&#95;lu | 23.6203 ns [23.5814, 23.6584] (0.95 confidence) | 23.2679 ns [23.2237, 23.3115] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;det&#95;via&#95;lu | 14298 ns [14280.6, 14315.1] (0.95 confidence) | 14256.3 ns [14236.3, 14277.6] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;dot | 11.1093 ns [11.0905, 11.1282] (0.95 confidence) | 11.0116 ns [10.9888, 11.0351] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;inf&#95;norm | 637.15 ns [635.707, 638.616] (0.95 confidence) | 631.864 ns [629.878, 634.021] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;ldlt | 16885.7 ns [16871.8, 16900.2] (0.95 confidence) | 16952.6 ns [16930.3, 16978.3] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;ldlt&#95;solve | 19187 ns [19130.7, 19250.7] (0.95 confidence) | 19163.5 ns [19134, 19192.9] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;lu | 15108.9 ns [15084.2, 15135.2] (0.95 confidence) | 15190.9 ns [15173.1, 15209.3] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;lu&#95;solve | 17192.8 ns [17165.9, 17220] (0.95 confidence) | 17527.9 ns [17451.7, 17605.5] (0.95 confidence) | 14548.5 ns [14519.8, 14579.6] (0.95 confidence) | 12861.4 ns [12845.4, 12877.9] (0.95 confidence) |
| d64/la&#95;stack&#95;norm2 | 59.5714 ns [59.3671, 59.7957] (0.95 confidence) | 57.8026 ns [57.6277, 57.9871] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;norm2&#95;scenario&#95;descending | 37.6369 ns [37.5508, 37.7264] (0.95 confidence) | 37.0365 ns [36.9239, 37.1599] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale | 37.3847 ns [37.3039, 37.4586] (0.95 confidence) | 36.917 ns [36.8359, 37] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;norm2&#95;scenario&#95;sparse | 17.3884 ns [17.3503, 17.4256] (0.95 confidence) | 17.1363 ns [17.0936, 17.1786] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range | 32.9493 ns [32.8579, 33.0376] (0.95 confidence) | 32.5757 ns [32.4979, 32.6506] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;norm2&#95;sq | 11.004 ns [10.9772, 11.032] (0.95 confidence) | 10.9397 ns [10.9096, 10.9717] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;solve&#95;from&#95;ldlt | 1096.41 ns [1094.52, 1098.36] (0.95 confidence) | 1027.43 ns [1025.49, 1029.44] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;solve&#95;from&#95;lu | 2615.26 ns [2605.68, 2625.84] (0.95 confidence) | 2559.14 ns [2550.98, 2568.73] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;det | 98.0726 ns [97.5683, 98.5805] (0.95 confidence) | 83.6659 ns [83.3796, 83.9485] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;det&#95;from&#95;ldlt | 1.31139 ns [1.30998, 1.31288] (0.95 confidence) | 1.31504 ns [1.31288, 1.31743] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;det&#95;from&#95;ldlt&#95;balanced&#95;range | 8.8479 ns [8.82393, 8.87186] (0.95 confidence) | 9.05012 ns [9.02776, 9.0733] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;det&#95;from&#95;lu | 1.47457 ns [1.47267, 1.47655] (0.95 confidence) | 1.47212 ns [1.47053, 1.47384] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;det&#95;from&#95;lu&#95;balanced&#95;range | 9.08828 ns [9.05797, 9.11985] (0.95 confidence) | 9.02866 ns [9.00587, 9.05233] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;det&#95;via&#95;lu | 90.8195 ns [90.2678, 91.3531] (0.95 confidence) | 92.9329 ns [92.5633, 93.3043] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;dot | 0.938013 ns [0.936592, 0.939457] (0.95 confidence) | 0.936178 ns [0.933581, 0.938855] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;inf&#95;norm | 8.47079 ns [8.458, 8.48266] (0.95 confidence) | 8.53422 ns [8.51667, 8.55455] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;ldlt | 91.0232 ns [90.8866, 91.1755] (0.95 confidence) | 91.7358 ns [91.3468, 92.1373] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;ldlt&#95;ill&#95;conditioned | 90.7641 ns [90.635, 90.8997] (0.95 confidence) | 93.1496 ns [92.1977, 94.499] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;ldlt&#95;solve | 104.304 ns [104.146, 104.467] (0.95 confidence) | 103.208 ns [102.902, 103.516] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;lu | 91.4551 ns [91.0942, 91.8053] (0.95 confidence) | 92.0077 ns [91.78, 92.2432] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;lu&#95;ill&#95;conditioned | 91.2176 ns [90.785, 91.6457] (0.95 confidence) | 90.8655 ns [90.516, 91.2148] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;lu&#95;pivoting | 95.7447 ns [95.4358, 96.0542] (0.95 confidence) | 94.4171 ns [94.2135, 94.6338] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;lu&#95;solve | 133.234 ns [132.867, 133.616] (0.95 confidence) | 132.273 ns [132.041, 132.506] (0.95 confidence) | 167.351 ns [166.824, 167.952] (0.95 confidence) | 342.616 ns [341.418, 343.841] (0.95 confidence) |
| d8/la&#95;stack&#95;norm2 | 4.68254 ns [4.66509, 4.70177] (0.95 confidence) | 4.65191 ns [4.6359, 4.671] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;norm2&#95;scenario&#95;descending | 4.84877 ns [4.84062, 4.85688] (0.95 confidence) | 4.84058 ns [4.83072, 4.85087] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale | 4.84462 ns [4.83486, 4.854] (0.95 confidence) | 4.81887 ns [4.81008, 4.82753] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;norm2&#95;scenario&#95;sparse | 2.62987 ns [2.5784, 2.68532] (0.95 confidence) | 2.4915 ns [2.47451, 2.51296] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range | 4.19056 ns [4.17478, 4.2054] (0.95 confidence) | 4.2651 ns [4.25808, 4.27184] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;norm2&#95;sq | 0.65899 ns [0.657492, 0.660437] (0.95 confidence) | 0.662481 ns [0.660877, 0.664079] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;solve&#95;from&#95;ldlt | 8.25672 ns [8.24375, 8.26921] (0.95 confidence) | 8.15193 ns [8.14108, 8.16325] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;solve&#95;from&#95;lu | 13.5422 ns [13.4902, 13.605] (0.95 confidence) | 13.4835 ns [13.4642, 13.5025] (0.95 confidence) | — | — |
| exact&#95;d2/det | 0.454783 ns [0.449938, 0.460081] (0.95 confidence) | 0.44612 ns [0.442103, 0.450547] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;direct | 0.501216 ns [0.500267, 0.502231] (0.95 confidence) | 0.4498 ns [0.445638, 0.454428] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;direct&#95;with&#95;errbound | 1.75259 ns [1.74939, 1.75597] (0.95 confidence) | 1.73695 ns [1.73469, 1.7393] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;errbound | 1.60091 ns [1.59771, 1.60465] (0.95 confidence) | 1.58529 ns [1.58458, 1.58606] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;exact | 87.7348 ns [87.6428, 87.8303] (0.95 confidence) | 88.1087 ns [87.8851, 88.3368] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;exact&#95;f64&#95;result | 64.8769 ns [64.785, 64.9712] (0.95 confidence) | 68.0866 ns [67.9467, 68.206] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;exact&#95;rounded&#95;f64 | 66.4164 ns [66.2951, 66.5423] (0.95 confidence) | 71.214 ns [71.1125, 71.3124] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;sign&#95;exact | 2.98731 ns [2.98419, 2.99088] (0.95 confidence) | 2.9648 ns [2.96183, 2.96748] (0.95 confidence) | — | — |
| exact&#95;d2/solve&#95;exact | 2957.41 ns [2953.08, 2961.76] (0.95 confidence) | 2943.3 ns [2937.55, 2949.38] (0.95 confidence) | — | — |
| exact&#95;d2/solve&#95;exact&#95;f64&#95;result | 3138.2 ns [3133.23, 3143.07] (0.95 confidence) | 3093.55 ns [3088, 3099.39] (0.95 confidence) | — | — |
| exact&#95;d2/solve&#95;exact&#95;rounded&#95;f64 | 3292.8 ns [3275.83, 3309.53] (0.95 confidence) | 3287.76 ns [3274.48, 3302.1] (0.95 confidence) | — | — |
| exact&#95;d3/det | 0.832426 ns [0.829634, 0.835379] (0.95 confidence) | 0.823044 ns [0.821101, 0.825177] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;direct | 0.8344 ns [0.831679, 0.837257] (0.95 confidence) | 0.822441 ns [0.821003, 0.824116] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;direct&#95;with&#95;errbound | 3.33451 ns [3.32918, 3.3402] (0.95 confidence) | 3.29644 ns [3.29189, 3.30097] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;errbound | 3.31456 ns [3.31157, 3.31773] (0.95 confidence) | 3.29182 ns [3.28861, 3.2954] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;exact | 300.943 ns [300.474, 301.456] (0.95 confidence) | 302.017 ns [301.693, 302.346] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;exact&#95;f64&#95;result | 277.874 ns [277.458, 278.314] (0.95 confidence) | 278.518 ns [278.035, 279.01] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;exact&#95;rounded&#95;f64 | 280.904 ns [280.47, 281.343] (0.95 confidence) | 280.418 ns [280.015, 280.823] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;sign&#95;exact | 4.73424 ns [4.72968, 4.73894] (0.95 confidence) | 4.69626 ns [4.69205, 4.70011] (0.95 confidence) | — | — |
| exact&#95;d3/solve&#95;exact | 8948.44 ns [8927.94, 8968.77] (0.95 confidence) | 9160.16 ns [9136.9, 9181.69] (0.95 confidence) | — | — |
| exact&#95;d3/solve&#95;exact&#95;f64&#95;result | 9073.82 ns [9051.95, 9095.51] (0.95 confidence) | 9278.55 ns [9258.51, 9296.71] (0.95 confidence) | — | — |
| exact&#95;d3/solve&#95;exact&#95;rounded&#95;f64 | 9401.04 ns [9374.44, 9427.99] (0.95 confidence) | 9772.61 ns [9725.35, 9820.59] (0.95 confidence) | — | — |
| exact&#95;d4/det | 3.22894 ns [3.22448, 3.23375] (0.95 confidence) | 3.21689 ns [3.20988, 3.22468] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;direct | 2.3663 ns [2.36356, 2.36926] (0.95 confidence) | 2.34837 ns [2.34681, 2.35007] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;direct&#95;with&#95;errbound | 6.80833 ns [6.80052, 6.81657] (0.95 confidence) | 6.75414 ns [6.74618, 6.7617] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;errbound | 6.81916 ns [6.79787, 6.8425] (0.95 confidence) | 6.72067 ns [6.70989, 6.73074] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;exact | 850.627 ns [849.237, 852.048] (0.95 confidence) | 831.217 ns [829.112, 833.448] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;exact&#95;f64&#95;result | 816.963 ns [815.766, 818.18] (0.95 confidence) | 791.784 ns [790.125, 793.514] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;exact&#95;rounded&#95;f64 | 834.214 ns [832.558, 835.919] (0.95 confidence) | 807.041 ns [804.962, 809.218] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;sign&#95;exact | 8.1593 ns [8.1508, 8.1687] (0.95 confidence) | 7.90718 ns [7.90232, 7.91222] (0.95 confidence) | — | — |
| exact&#95;d4/solve&#95;exact | 18648.6 ns [18609.7, 18686.3] (0.95 confidence) | 18573.1 ns [18525, 18620.3] (0.95 confidence) | — | — |
| exact&#95;d4/solve&#95;exact&#95;f64&#95;result | 18689.9 ns [18662, 18718.3] (0.95 confidence) | 18689.6 ns [18645.3, 18732.5] (0.95 confidence) | — | — |
| exact&#95;d4/solve&#95;exact&#95;rounded&#95;f64 | 19305.3 ns [19246.7, 19366.5] (0.95 confidence) | 19394.7 ns [19363.7, 19426.6] (0.95 confidence) | — | — |
| exact&#95;d5/det | 27.5967 ns [27.4998, 27.6937] (0.95 confidence) | 28.3426 ns [28.3091, 28.3749] (0.95 confidence) | — | — |
| exact&#95;d5/det&#95;exact | 2921.01 ns [2915.43, 2926.64] (0.95 confidence) | 2881.55 ns [2871.83, 2891.96] (0.95 confidence) | — | — |
| exact&#95;d5/det&#95;exact&#95;f64&#95;result | 2883.28 ns [2878.02, 2888.52] (0.95 confidence) | 2825.06 ns [2820.16, 2830.1] (0.95 confidence) | — | — |
| exact&#95;d5/det&#95;exact&#95;rounded&#95;f64 | 2930 ns [2922.49, 2937.73] (0.95 confidence) | 2844.06 ns [2836.01, 2852.45] (0.95 confidence) | — | — |
| exact&#95;d5/det&#95;sign&#95;exact | 2906.18 ns [2900.91, 2911.53] (0.95 confidence) | 2815.92 ns [2809.45, 2822.69] (0.95 confidence) | — | — |
| exact&#95;d5/solve&#95;exact | 30259.3 ns [30167.9, 30352.1] (0.95 confidence) | 30411.7 ns [30338.1, 30480.9] (0.95 confidence) | — | — |
| exact&#95;d5/solve&#95;exact&#95;f64&#95;result | 30465.7 ns [30355, 30586] (0.95 confidence) | 30628.7 ns [30571.4, 30684.7] (0.95 confidence) | — | — |
| exact&#95;d5/solve&#95;exact&#95;rounded&#95;f64 | 31428.6 ns [31326.7, 31534.1] (0.95 confidence) | 31392.3 ns [31300.1, 31488] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/det&#95;exact | 947.804 ns [946.243, 949.458] (0.95 confidence) | 917.843 ns [916.083, 919.679] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/det&#95;sign&#95;exact | 8.12303 ns [8.11785, 8.12898] (0.95 confidence) | 7.90733 ns [7.90357, 7.91126] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/solve&#95;exact | 15286.7 ns [15241.8, 15332.2] (0.95 confidence) | 15286.1 ns [15247.4, 15323.2] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/solve&#95;exact&#95;f64&#95;result | 15297.9 ns [15256.2, 15339.3] (0.95 confidence) | 15522.9 ns [15474.4, 15570.1] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/solve&#95;exact&#95;rounded&#95;f64 | 15724.5 ns [15683.7, 15764.6] (0.95 confidence) | 15778.1 ns [15734.3, 15820.2] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/det&#95;exact | 2951.11 ns [2945.78, 2956.66] (0.95 confidence) | 2961.73 ns [2954.88, 2969.28] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/det&#95;sign&#95;exact | 2898.84 ns [2894.24, 2903.55] (0.95 confidence) | 2852.01 ns [2844.11, 2860.27] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/solve&#95;exact | 24708.6 ns [24634.8, 24784.5] (0.95 confidence) | 24779.3 ns [24727, 24830.2] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/solve&#95;exact&#95;f64&#95;result | 24773.5 ns [24728.5, 24819.3] (0.95 confidence) | 25127.1 ns [25093.3, 25159.4] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/solve&#95;exact&#95;rounded&#95;f64 | 25723.2 ns [25660.8, 25788.5] (0.95 confidence) | 25743.8 ns [25683.7, 25802] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/det&#95;exact | 317.532 ns [316.795, 318.282] (0.95 confidence) | 318.32 ns [314.976, 323.691] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/det&#95;sign&#95;exact | 290.579 ns [290.186, 290.981] (0.95 confidence) | 283.747 ns [283.081, 284.441] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/solve&#95;exact | 9040.44 ns [9025.74, 9055.23] (0.95 confidence) | 8884.54 ns [8869.57, 8899.37] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/solve&#95;exact&#95;f64&#95;result | 9177.84 ns [9164.94, 9190.58] (0.95 confidence) | 8994.26 ns [8976.83, 9009.5] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/solve&#95;exact&#95;rounded&#95;f64 | 9487.28 ns [9476.57, 9498.2] (0.95 confidence) | 9386.34 ns [9369.54, 9403.07] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/det&#95;exact | 183.473 ns [183.263, 183.698] (0.95 confidence) | 187.471 ns [187.307, 187.63] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/det&#95;sign&#95;exact | 163.575 ns [163.341, 163.815] (0.95 confidence) | 164.073 ns [163.927, 164.231] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/solve&#95;exact | 975.223 ns [973.541, 976.973] (0.95 confidence) | 966.438 ns [964.262, 968.704] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/solve&#95;exact&#95;f64&#95;result | 998.164 ns [996.601, 999.721] (0.95 confidence) | 1001.38 ns [999.34, 1003.36] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/solve&#95;exact&#95;rounded&#95;f64 | 987.851 ns [985.967, 989.822] (0.95 confidence) | 981.634 ns [980.439, 982.847] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/det&#95;exact | 3483.41 ns [3473, 3494.14] (0.95 confidence) | 3440.59 ns [3429.32, 3452.35] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/det&#95;sign&#95;exact | 152.073 ns [151.586, 152.582] (0.95 confidence) | 145.485 ns [145.151, 145.829] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/solve&#95;exact | 37311.9 ns [37204.4, 37421.9] (0.95 confidence) | 35224.2 ns [35132.3, 35320.5] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/solve&#95;exact&#95;f64&#95;result | 38113.9 ns [38020.9, 38208.2] (0.95 confidence) | 36087 ns [36022.4, 36150.8] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/solve&#95;exact&#95;rounded&#95;f64 | 38399.1 ns [38304.3, 38494] (0.95 confidence) | 36170.8 ns [36069.8, 36276] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/det&#95;exact | 9657.87 ns [9645.76, 9669.8] (0.95 confidence) | 9390.2 ns [9369.84, 9412.13] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/det&#95;sign&#95;exact | 231.778 ns [231.604, 231.969] (0.95 confidence) | 224.052 ns [223.517, 224.447] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/solve&#95;exact | 84886.5 ns [84763.3, 85009.7] (0.95 confidence) | 79297.2 ns [79157.5, 79438.9] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/solve&#95;exact&#95;f64&#95;result | 85654.2 ns [85536, 85778.5] (0.95 confidence) | 79844.5 ns [79703.9, 79984.3] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/solve&#95;exact&#95;rounded&#95;f64 | 86166.3 ns [86066.3, 86271.2] (0.95 confidence) | 80616.5 ns [80504.1, 80729.8] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/det&#95;exact | 23781.3 ns [23736.9, 23827.4] (0.95 confidence) | 22680.4 ns [22638.4, 22722.1] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/det&#95;sign&#95;exact | 464.973 ns [464.057, 465.938] (0.95 confidence) | 439.352 ns [438.814, 439.888] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/solve&#95;exact | 152394 ns [152176, 152617] (0.95 confidence) | 142382 ns [142106, 142682] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/solve&#95;exact&#95;f64&#95;result | 154145 ns [153914, 154382] (0.95 confidence) | 141900 ns [141689, 142109] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/solve&#95;exact&#95;rounded&#95;f64 | 154880 ns [154690, 155085] (0.95 confidence) | 143209 ns [143001, 143416] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/det&#95;exact | 58687.7 ns [58514.6, 58879.2] (0.95 confidence) | 59073.3 ns [58988.4, 59159] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/det&#95;sign&#95;exact | 57436.8 ns [57232.5, 57683.9] (0.95 confidence) | 57750.9 ns [57623.6, 57869] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/solve&#95;exact | 256085 ns [255351, 256856] (0.95 confidence) | 238998 ns [238480, 239537] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/solve&#95;exact&#95;f64&#95;result | 258505 ns [257485, 259558] (0.95 confidence) | 241036 ns [240606, 241484] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/solve&#95;exact&#95;rounded&#95;f64 | 264500 ns [264106, 264908] (0.95 confidence) | 240958 ns [240572, 241340] (0.95 confidence) | — | — |
| gram/1x2/dependent/la&#95;stack | 0.691288 ns [0.689963, 0.692682] (0.95 confidence) | 0.697442 ns [0.69574, 0.69928] (0.95 confidence) | — | — |
| gram/1x2/mixed&#95;scale/la&#95;stack | 0.695483 ns [0.69209, 0.700302] (0.95 confidence) | 0.696512 ns [0.693883, 0.699234] (0.95 confidence) | — | — |
| gram/1x2/near&#95;dependent/la&#95;stack | 0.691187 ns [0.689977, 0.692456] (0.95 confidence) | 0.700072 ns [0.698122, 0.702311] (0.95 confidence) | — | — |
| gram/1x2/orthogonal/la&#95;stack | 0.691756 ns [0.689819, 0.694365] (0.95 confidence) | 0.696204 ns [0.694678, 0.697804] (0.95 confidence) | — | — |
| gram/2x2/dependent/la&#95;stack | 1.3672 ns [1.3629, 1.37176] (0.95 confidence) | 1.32653 ns [1.32019, 1.33283] (0.95 confidence) | — | — |
| gram/2x2/mixed&#95;scale/la&#95;stack | 1.34559 ns [1.34126, 1.35029] (0.95 confidence) | 1.32978 ns [1.32369, 1.33596] (0.95 confidence) | — | — |
| gram/2x2/near&#95;dependent/la&#95;stack | 1.36747 ns [1.36214, 1.37388] (0.95 confidence) | 1.33113 ns [1.32463, 1.33763] (0.95 confidence) | — | — |
| gram/2x2/orthogonal/la&#95;stack | 1.34542 ns [1.34283, 1.34795] (0.95 confidence) | 1.32963 ns [1.32303, 1.33645] (0.95 confidence) | — | — |
| gram/2x3/dependent/la&#95;stack | 1.32209 ns [1.31967, 1.32492] (0.95 confidence) | 1.26651 ns [1.265, 1.26804] (0.95 confidence) | — | — |
| gram/2x3/mixed&#95;scale/la&#95;stack | 1.32285 ns [1.32091, 1.32487] (0.95 confidence) | 1.30008 ns [1.29428, 1.30695] (0.95 confidence) | — | — |
| gram/2x3/near&#95;dependent/la&#95;stack | 1.31996 ns [1.31811, 1.32188] (0.95 confidence) | 1.26573 ns [1.26408, 1.26744] (0.95 confidence) | — | — |
| gram/2x3/orthogonal/la&#95;stack | 1.32644 ns [1.32033, 1.3344] (0.95 confidence) | 1.27025 ns [1.26709, 1.27489] (0.95 confidence) | — | — |
| gram/3x3/dependent/la&#95;stack | 3.04469 ns [3.04057, 3.04936] (0.95 confidence) | 3.00552 ns [3.0025, 3.00901] (0.95 confidence) | — | — |
| gram/3x3/mixed&#95;scale/la&#95;stack | 3.04444 ns [3.03927, 3.05057] (0.95 confidence) | 3.00983 ns [3.00621, 3.01384] (0.95 confidence) | — | — |
| gram/3x3/near&#95;dependent/la&#95;stack | 3.04325 ns [3.03994, 3.04673] (0.95 confidence) | 3.00828 ns [3.00492, 3.01205] (0.95 confidence) | — | — |
| gram/3x3/orthogonal/la&#95;stack | 3.05093 ns [3.04187, 3.06297] (0.95 confidence) | 3.00498 ns [3.00254, 3.00757] (0.95 confidence) | — | — |
| gram/3x4/dependent/la&#95;stack | 3.07897 ns [3.07494, 3.08291] (0.95 confidence) | 3.04852 ns [3.04638, 3.05088] (0.95 confidence) | — | — |
| gram/3x4/mixed&#95;scale/la&#95;stack | 3.08218 ns [3.07806, 3.08662] (0.95 confidence) | 3.05321 ns [3.04993, 3.05682] (0.95 confidence) | — | — |
| gram/3x4/near&#95;dependent/la&#95;stack | 3.0777 ns [3.07452, 3.08141] (0.95 confidence) | 3.04582 ns [3.04303, 3.04839] (0.95 confidence) | — | — |
| gram/3x4/orthogonal/la&#95;stack | 3.08003 ns [3.07665, 3.08367] (0.95 confidence) | 3.04934 ns [3.04611, 3.05275] (0.95 confidence) | — | — |
| gram/4x4/dependent/la&#95;stack | 7.67246 ns [7.66091, 7.68432] (0.95 confidence) | 7.61777 ns [7.56963, 7.6654] (0.95 confidence) | — | — |
| gram/4x4/mixed&#95;scale/la&#95;stack | 7.66536 ns [7.65486, 7.67587] (0.95 confidence) | 7.6063 ns [7.5596, 7.65289] (0.95 confidence) | — | — |
| gram/4x4/near&#95;dependent/la&#95;stack | 7.66213 ns [7.64935, 7.67509] (0.95 confidence) | 7.62139 ns [7.57592, 7.6668] (0.95 confidence) | — | — |
| gram/4x4/orthogonal/la&#95;stack | 7.69969 ns [7.68636, 7.71378] (0.95 confidence) | 7.65077 ns [7.6017, 7.70113] (0.95 confidence) | — | — |
| gram/4x5/dependent/la&#95;stack | 9.78864 ns [9.70945, 9.87217] (0.95 confidence) | 9.39848 ns [9.36328, 9.43504] (0.95 confidence) | — | — |
| gram/4x5/mixed&#95;scale/la&#95;stack | 9.55299 ns [9.5191, 9.59315] (0.95 confidence) | 9.4738 ns [9.46017, 9.48627] (0.95 confidence) | — | — |
| gram/4x5/near&#95;dependent/la&#95;stack | 9.56836 ns [9.51644, 9.62765] (0.95 confidence) | 9.45292 ns [9.43377, 9.47019] (0.95 confidence) | — | — |
| gram/4x5/orthogonal/la&#95;stack | 9.58809 ns [9.54383, 9.63772] (0.95 confidence) | 9.46572 ns [9.45144, 9.47877] (0.95 confidence) | — | — |
| gram/5x5/dependent/la&#95;stack | 18.2504 ns [18.2024, 18.2981] (0.95 confidence) | 18.3151 ns [18.2419, 18.3825] (0.95 confidence) | — | — |
| gram/5x5/mixed&#95;scale/la&#95;stack | 18.3387 ns [18.2775, 18.4019] (0.95 confidence) | 18.2585 ns [18.1807, 18.3317] (0.95 confidence) | — | — |
| gram/5x5/near&#95;dependent/la&#95;stack | 18.0695 ns [18.0137, 18.1255] (0.95 confidence) | 18.3131 ns [18.2592, 18.3645] (0.95 confidence) | — | — |
| gram/5x5/orthogonal/la&#95;stack | 18.1654 ns [18.1055, 18.2262] (0.95 confidence) | 18.1255 ns [18.0278, 18.221] (0.95 confidence) | — | — |
| gram/5x6/dependent/la&#95;stack | 19.1514 ns [19.0881, 19.2139] (0.95 confidence) | 19.1628 ns [19.1205, 19.2022] (0.95 confidence) | — | — |
| gram/5x6/mixed&#95;scale/la&#95;stack | 19.2473 ns [19.1988, 19.2931] (0.95 confidence) | 19.0764 ns [19.0195, 19.1323] (0.95 confidence) | — | — |
| gram/5x6/near&#95;dependent/la&#95;stack | 19.2474 ns [19.1849, 19.3094] (0.95 confidence) | 18.9812 ns [18.9059, 19.0543] (0.95 confidence) | — | — |
| gram/5x6/orthogonal/la&#95;stack | 19.1732 ns [19.1349, 19.2106] (0.95 confidence) | 18.9669 ns [18.9106, 19.0214] (0.95 confidence) | — | — |
| gram/6x6/dependent/la&#95;stack | 29.621 ns [29.509, 29.7337] (0.95 confidence) | 29.4488 ns [29.3042, 29.5939] (0.95 confidence) | — | — |
| gram/6x6/mixed&#95;scale/la&#95;stack | 29.7246 ns [29.5978, 29.8994] (0.95 confidence) | 29.395 ns [29.2341, 29.5617] (0.95 confidence) | — | — |
| gram/6x6/near&#95;dependent/la&#95;stack | 29.635 ns [29.5146, 29.7718] (0.95 confidence) | 29.6386 ns [29.4797, 29.7982] (0.95 confidence) | — | — |
| gram/6x6/orthogonal/la&#95;stack | 29.8517 ns [29.7781, 29.9219] (0.95 confidence) | 29.5701 ns [29.3904, 29.751] (0.95 confidence) | — | — |
| gram/6x7/dependent/la&#95;stack | 32.5496 ns [32.4445, 32.6512] (0.95 confidence) | 31.8265 ns [31.7207, 31.9295] (0.95 confidence) | — | — |
| gram/6x7/mixed&#95;scale/la&#95;stack | 32.5706 ns [32.5058, 32.6294] (0.95 confidence) | 32.0235 ns [31.9107, 32.1314] (0.95 confidence) | — | — |
| gram/6x7/near&#95;dependent/la&#95;stack | 32.2997 ns [32.1856, 32.4058] (0.95 confidence) | 32.0398 ns [31.9286, 32.1514] (0.95 confidence) | — | — |
| gram/6x7/orthogonal/la&#95;stack | 32.4785 ns [32.4239, 32.5379] (0.95 confidence) | 32.0921 ns [31.9655, 32.2201] (0.95 confidence) | — | — |
| gram/7x7/dependent/la&#95;stack | 42.3001 ns [42.0257, 42.5821] (0.95 confidence) | 39.6052 ns [39.2246, 40.0266] (0.95 confidence) | — | — |
| gram/7x7/mixed&#95;scale/la&#95;stack | 41.7736 ns [41.5319, 42.0296] (0.95 confidence) | 39.5013 ns [39.2418, 39.7893] (0.95 confidence) | — | — |
| gram/7x7/near&#95;dependent/la&#95;stack | 42.2881 ns [42.0704, 42.5201] (0.95 confidence) | 40.248 ns [39.9407, 40.5814] (0.95 confidence) | — | — |
| gram/7x7/orthogonal/la&#95;stack | 43.8727 ns [43.6079, 44.1347] (0.95 confidence) | 39.0846 ns [38.9457, 39.2351] (0.95 confidence) | — | — |
| gram/7x8/dependent/la&#95;stack | 42.3153 ns [42.1614, 42.4703] (0.95 confidence) | 40.1367 ns [40.0724, 40.2091] (0.95 confidence) | — | — |
| gram/7x8/mixed&#95;scale/la&#95;stack | 43.0446 ns [42.7668, 43.3353] (0.95 confidence) | 40.2136 ns [40.1065, 40.3405] (0.95 confidence) | — | — |
| gram/7x8/near&#95;dependent/la&#95;stack | 42.6745 ns [42.3372, 43.0259] (0.95 confidence) | 40.4871 ns [40.2769, 40.7304] (0.95 confidence) | — | — |
| gram/7x8/orthogonal/la&#95;stack | 43.1579 ns [42.7922, 43.5365] (0.95 confidence) | 40.5194 ns [40.2236, 40.8562] (0.95 confidence) | — | — |
| gram/8x8/dependent/la&#95;stack | 52.6214 ns [52.4312, 52.8224] (0.95 confidence) | 51.8777 ns [51.7744, 51.9961] (0.95 confidence) | — | — |
| gram/8x8/mixed&#95;scale/la&#95;stack | 53.0731 ns [52.8391, 53.317] (0.95 confidence) | 52.2938 ns [52.1048, 52.5007] (0.95 confidence) | — | — |
| gram/8x8/near&#95;dependent/la&#95;stack | 52.2506 ns [52.1169, 52.4019] (0.95 confidence) | 52.1104 ns [51.944, 52.2938] (0.95 confidence) | — | — |
| gram/8x8/orthogonal/la&#95;stack | 53.7736 ns [53.588, 53.9632] (0.95 confidence) | 51.9483 ns [51.8043, 52.1091] (0.95 confidence) | — | — |
| interval&#95;d3/dense&#95;inconclusive | 43.9686 ns [43.705, 44.2459] (0.95 confidence) | 41.32 ns [41.2477, 41.3951] (0.95 confidence) | — | — |
| interval&#95;d3/dense&#95;sign | 40.6015 ns [40.2983, 40.9148] (0.95 confidence) | 38.6591 ns [38.2112, 39.1133] (0.95 confidence) | — | — |
| interval&#95;d3/lifted&#95;assembly | 8.37993 ns [8.37128, 8.39027] (0.95 confidence) | 8.34388 ns [8.33595, 8.35185] (0.95 confidence) | — | — |
| interval&#95;d3/lifted&#95;assembly&#95;and&#95;sign | 37.5666 ns [37.3081, 37.8451] (0.95 confidence) | 36.9108 ns [36.7559, 37.07] (0.95 confidence) | — | — |
| interval&#95;d3/sparse&#95;sign | 32.4448 ns [32.1129, 32.7905] (0.95 confidence) | 31.6599 ns [31.3849, 31.9405] (0.95 confidence) | — | — |
| interval&#95;d4/dense&#95;inconclusive | 108.767 ns [107.82, 109.758] (0.95 confidence) | 102.565 ns [101.879, 103.429] (0.95 confidence) | — | — |
| interval&#95;d4/dense&#95;sign | 103.922 ns [103.138, 104.726] (0.95 confidence) | 99.9014 ns [99.229, 100.584] (0.95 confidence) | — | — |
| interval&#95;d4/lifted&#95;assembly | 13.128 ns [13.1032, 13.1517] (0.95 confidence) | 13.2624 ns [13.248, 13.2766] (0.95 confidence) | — | — |
| interval&#95;d4/lifted&#95;assembly&#95;and&#95;sign | 85.3568 ns [84.893, 85.8597] (0.95 confidence) | 81.9549 ns [81.6115, 82.3181] (0.95 confidence) | — | — |
| interval&#95;d4/sparse&#95;sign | 71.2349 ns [70.7125, 71.7901] (0.95 confidence) | 68.4665 ns [68.2196, 68.7179] (0.95 confidence) | — | — |
| interval&#95;d5/dense&#95;inconclusive | 245.46 ns [244.651, 246.29] (0.95 confidence) | 237.666 ns [236.576, 238.811] (0.95 confidence) | — | — |
| interval&#95;d5/dense&#95;sign | 235.465 ns [233.821, 237.22] (0.95 confidence) | 224.915 ns [223.741, 226.136] (0.95 confidence) | — | — |
| interval&#95;d5/lifted&#95;assembly | 21.6352 ns [21.6125, 21.6579] (0.95 confidence) | 21.7896 ns [21.704, 21.8795] (0.95 confidence) | — | — |
| interval&#95;d5/lifted&#95;assembly&#95;and&#95;sign | 197.243 ns [195.988, 198.545] (0.95 confidence) | 194.888 ns [193.453, 196.314] (0.95 confidence) | — | — |
| interval&#95;d5/sparse&#95;sign | 163.278 ns [162.142, 164.405] (0.95 confidence) | 152.512 ns [152.106, 152.922] (0.95 confidence) | — | — |
| interval&#95;d6/dense&#95;inconclusive | 565.855 ns [563.652, 568.133] (0.95 confidence) | 541.679 ns [539.842, 543.484] (0.95 confidence) | — | — |
| interval&#95;d6/dense&#95;sign | 561.571 ns [559.127, 564.004] (0.95 confidence) | 539.508 ns [537.135, 541.866] (0.95 confidence) | — | — |
| interval&#95;d6/lifted&#95;assembly | 29.9073 ns [29.6892, 30.1265] (0.95 confidence) | 30.8582 ns [30.7814, 30.941] (0.95 confidence) | — | — |
| interval&#95;d6/lifted&#95;assembly&#95;and&#95;sign | 406.906 ns [405.556, 408.241] (0.95 confidence) | 385.747 ns [384.712, 386.712] (0.95 confidence) | — | — |
| interval&#95;d6/sparse&#95;sign | 379.632 ns [378.373, 380.93] (0.95 confidence) | 361.311 ns [359.685, 362.958] (0.95 confidence) | — | — |
| interval&#95;d7/dense&#95;inconclusive | 1121.25 ns [1117.06, 1125.55] (0.95 confidence) | 1104.9 ns [1099.67, 1110.08] (0.95 confidence) | — | — |
| interval&#95;d7/dense&#95;sign | 1103 ns [1097.8, 1108.38] (0.95 confidence) | 1066.95 ns [1063.32, 1070.59] (0.95 confidence) | — | — |
| interval&#95;d7/lifted&#95;assembly | 38.1994 ns [38.1049, 38.3367] (0.95 confidence) | 37.2406 ns [37.1846, 37.2942] (0.95 confidence) | — | — |
| interval&#95;d7/lifted&#95;assembly&#95;and&#95;sign | 835.048 ns [831.681, 838.411] (0.95 confidence) | 789.738 ns [788.229, 791.177] (0.95 confidence) | — | — |
| interval&#95;d7/sparse&#95;sign | 786.57 ns [784.251, 788.882] (0.95 confidence) | 773.622 ns [770.628, 776.476] (0.95 confidence) | — | — |
| interval&#95;det&#95;sign/d4&#95;conclusive&#95;lifted | 71.8087 ns [71.1605, 72.471] (0.95 confidence) | 69.4616 ns [68.918, 70.0258] (0.95 confidence) | — | — |
| interval&#95;det&#95;sign/d4&#95;inconclusive&#95;lifted | 68.9024 ns [68.4715, 69.3327] (0.95 confidence) | 65.217 ns [64.8922, 65.525] (0.95 confidence) | — | — |
| interval&#95;det&#95;sign/d7&#95;conclusive&#95;lifted | 788.945 ns [785.864, 791.984] (0.95 confidence) | 757.958 ns [756.398, 759.439] (0.95 confidence) | — | — |
| interval&#95;division/negative | — | 5.20326 ns [5.19776, 5.20906] (0.95 confidence) | — | — |
| interval&#95;division/positive | — | 5.21728 ns [5.21171, 5.22257] (0.95 confidence) | — | — |
| interval&#95;division/reciprocal | — | 10.0324 ns [10.0101, 10.0549] (0.95 confidence) | — | — |
| interval&#95;division/subnormal | — | 3.87826 ns [3.86276, 3.89458] (0.95 confidence) | — | — |
| interval&#95;division/tiny&#95;ratio | — | 6.10193 ns [6.09268, 6.111] (0.95 confidence) | — | — |
| interval&#95;scalar/cancellation/add | 2.15785 ns [2.15519, 2.16106] (0.95 confidence) | 2.11879 ns [2.11624, 2.12124] (0.95 confidence) | — | — |
| interval&#95;scalar/cancellation/multiply | 1.88011 ns [1.87677, 1.88382] (0.95 confidence) | 1.87582 ns [1.87225, 1.87945] (0.95 confidence) | — | — |
| interval&#95;scalar/cancellation/square | 1.90357 ns [1.90162, 1.90549] (0.95 confidence) | 1.86397 ns [1.86106, 1.86675] (0.95 confidence) | — | — |
| interval&#95;scalar/exact/add | 2.15951 ns [2.15676, 2.16269] (0.95 confidence) | 2.18616 ns [2.1679, 2.20556] (0.95 confidence) | — | — |
| interval&#95;scalar/exact/multiply | 2.04157 ns [2.02921, 2.05479] (0.95 confidence) | 2.08714 ns [2.07369, 2.09994] (0.95 confidence) | — | — |
| interval&#95;scalar/exact/square | 1.91258 ns [1.90686, 1.92011] (0.95 confidence) | 1.86132 ns [1.85941, 1.86319] (0.95 confidence) | — | — |
| interval&#95;scalar/point/add | 2.19269 ns [2.18977, 2.19574] (0.95 confidence) | 2.18621 ns [2.18071, 2.19329] (0.95 confidence) | — | — |
| interval&#95;scalar/point/multiply | 2.67442 ns [2.6604, 2.68912] (0.95 confidence) | 2.61265 ns [2.5982, 2.62836] (0.95 confidence) | — | — |
| interval&#95;scalar/point/square | 2.14572 ns [2.1429, 2.14881] (0.95 confidence) | 2.12102 ns [2.1182, 2.12364] (0.95 confidence) | — | — |
| interval&#95;scalar/subnormal/add | 2.42888 ns [2.42665, 2.43139] (0.95 confidence) | 2.40313 ns [2.39996, 2.40633] (0.95 confidence) | — | — |
| interval&#95;scalar/subnormal/multiply | 3.03723 ns [3.02865, 3.04661] (0.95 confidence) | 2.85319 ns [2.84238, 2.86434] (0.95 confidence) | — | — |
| interval&#95;scalar/subnormal/square | 2.3995 ns [2.39695, 2.40249] (0.95 confidence) | 2.37825 ns [2.37499, 2.38206] (0.95 confidence) | — | — |
| interval&#95;scalar/zero/add | 0.950002 ns [0.948839, 0.951264] (0.95 confidence) | 0.938252 ns [0.937252, 0.939229] (0.95 confidence) | — | — |
| interval&#95;scalar/zero/multiply | 1.79525 ns [1.79229, 1.79833] (0.95 confidence) | 1.76518 ns [1.76203, 1.76847] (0.95 confidence) | — | — |
| interval&#95;scalar/zero/square | 0.894211 ns [0.89267, 0.896079] (0.95 confidence) | 0.881382 ns [0.880431, 0.882315] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d1/dense | — | 227.941 ns [227.485, 228.373] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d1/near&#95;singular | — | 291.937 ns [290.538, 294.476] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d1/range&#95;exhausted | — | 155.447 ns [154.984, 155.976] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d1/singular | — | 51.4406 ns [51.21, 51.6885] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d1/uncertain&#95;dense | — | 372.578 ns [372.043, 373.104] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d2/dense | — | 2892.28 ns [2886.26, 2897.91] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d2/near&#95;singular | — | 1593.16 ns [1590.19, 1596.13] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d2/range&#95;exhausted | — | 253.499 ns [253.358, 253.635] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d2/singular | — | 167.885 ns [167.66, 168.121] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d2/uncertain&#95;dense | — | 3983.73 ns [3975.89, 3991.21] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d3/dense | — | 7421.99 ns [7410.59, 7433.04] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d3/near&#95;singular | — | 2023.02 ns [2012.21, 2035.39] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d3/range&#95;exhausted | — | 409.94 ns [408.819, 410.771] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d3/singular | — | 322.298 ns [321.764, 322.852] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d3/uncertain&#95;dense | — | 9058.11 ns [9044.09, 9071.42] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d4/dense | — | 14522.1 ns [14503.6, 14540] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d4/near&#95;singular | — | 2588.86 ns [2585.49, 2592.07] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d4/range&#95;exhausted | — | 756.616 ns [755.693, 757.64] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d4/singular | — | 576.519 ns [575.925, 577.065] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d4/uncertain&#95;dense | — | 17140.9 ns [17109, 17169.9] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d5/dense | — | 24996.6 ns [24950.3, 25044.3] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d5/near&#95;singular | — | 3438.63 ns [3434.33, 3442.85] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d5/range&#95;exhausted | — | 1248.94 ns [1247.35, 1250.44] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d5/singular | — | 1015.89 ns [1014.55, 1017.2] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d5/uncertain&#95;dense | — | 28450.9 ns [28392.8, 28507.3] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d6/dense | — | 39232.3 ns [39178.3, 39284.3] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d6/near&#95;singular | — | 4318.33 ns [4312.88, 4323.8] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d6/range&#95;exhausted | — | 1578.92 ns [1576.88, 1581.01] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d6/singular | — | 1314.21 ns [1312.33, 1316.09] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d6/uncertain&#95;dense | — | 44829.6 ns [44678.1, 44989.1] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;batch/origin/construction | 1.94486 ns [1.93713, 1.95218] (0.95 confidence) | 2.05399 ns [2.0508, 2.05703] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;batch/origin/prepared | 12.7947 ns [12.7809, 12.8101] (0.95 confidence) | 13.8026 ns [13.791, 13.8135] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;batch/origin/with&#95;construction | 14.3408 ns [14.3274, 14.3559] (0.95 confidence) | 15.3831 ns [15.356, 15.4147] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;batch/translated/construction | 1.9617 ns [1.95322, 1.96995] (0.95 confidence) | 2.05699 ns [2.05422, 2.05961] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;batch/translated/prepared | 42.0754 ns [41.9704, 42.19] (0.95 confidence) | 40.4104 ns [40.2969, 40.5309] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;batch/translated/with&#95;construction | 43.8221 ns [43.7105, 43.9446] (0.95 confidence) | 42.2964 ns [42.1355, 42.4739] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/cancellation/difference | 16.5123 ns [16.4882, 16.5357] (0.95 confidence) | 10.4083 ns [10.3252, 10.5093] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/cancellation/difference&#95;endpoints | 17.5592 ns [17.5317, 17.5879] (0.95 confidence) | 11.7016 ns [11.5773, 11.8435] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/cancellation/dot | 5.52474 ns [5.4929, 5.56104] (0.95 confidence) | 6.12315 ns [6.07935, 6.16898] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/cancellation/dot&#95;endpoints | 6.50064 ns [6.49095, 6.51083] (0.95 confidence) | 7.067 ns [7.04716, 7.0881] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/dense/difference | 15.9986 ns [15.9718, 16.0286] (0.95 confidence) | 10.4372 ns [10.3493, 10.5385] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/dense/difference&#95;endpoints | 19.6168 ns [19.5721, 19.6614] (0.95 confidence) | 13.1 ns [12.9823, 13.2353] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/dense/dot | 5.22143 ns [5.19469, 5.25177] (0.95 confidence) | 5.54044 ns [5.52529, 5.55845] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/dense/dot&#95;endpoints | 5.99465 ns [5.97616, 6.01545] (0.95 confidence) | 6.45964 ns [6.43459, 6.48713] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/dense&#95;inexact/difference | 16.0261 ns [15.9987, 16.0552] (0.95 confidence) | 10.4648 ns [10.3878, 10.555] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/dense&#95;inexact/difference&#95;endpoints | 17.1933 ns [17.1662, 17.22] (0.95 confidence) | 11.7622 ns [11.6651, 11.8671] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/dense&#95;inexact/dot | 5.22662 ns [5.19757, 5.25958] (0.95 confidence) | 5.52435 ns [5.51347, 5.53583] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/dense&#95;inexact/dot&#95;endpoints | 5.97822 ns [5.96647, 5.99094] (0.95 confidence) | 6.54018 ns [6.51324, 6.57] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/overflow/dot | 1.2435 ns [1.2424, 1.2447] (0.95 confidence) | 1.07248 ns [1.07056, 1.07452] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/sparse/difference | 6.49866 ns [6.49072, 6.50751] (0.95 confidence) | 4.76367 ns [4.74552, 4.78535] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/sparse/difference&#95;endpoints | 7.2662 ns [7.25236, 7.28002] (0.95 confidence) | 5.53412 ns [5.52682, 5.54128] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/sparse/dot | 2.23743 ns [2.2354, 2.2395] (0.95 confidence) | 2.22388 ns [2.22224, 2.22544] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/sparse/dot&#95;endpoints | 1.99118 ns [1.98945, 1.99313] (0.95 confidence) | 2.22379 ns [2.22255, 2.22501] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/underflow/difference | 5.88818 ns [5.88146, 5.895] (0.95 confidence) | 4.45819 ns [4.45517, 4.4619] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/underflow/difference&#95;endpoints | 5.95928 ns [5.95247, 5.96679] (0.95 confidence) | 4.60642 ns [4.60228, 4.61113] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/underflow/dot | 1.49373 ns [1.4916, 1.49606] (0.95 confidence) | 1.48147 ns [1.47991, 1.48322] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/underflow/dot&#95;endpoints | 1.54992 ns [1.54551, 1.55409] (0.95 confidence) | 1.48311 ns [1.48232, 1.48391] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;batch/origin/construction | 2.05169 ns [2.04501, 2.05849] (0.95 confidence) | 2.24248 ns [2.23402, 2.25071] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;batch/origin/prepared | 24.7524 ns [24.6945, 24.8148] (0.95 confidence) | 24.095 ns [24.0145, 24.178] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;batch/origin/with&#95;construction | 28.8884 ns [28.7876, 28.9933] (0.95 confidence) | 28.3683 ns [28.2808, 28.4556] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;batch/translated/construction | 2.0603 ns [2.05482, 2.06598] (0.95 confidence) | 2.25088 ns [2.24358, 2.25789] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;batch/translated/prepared | 82.5715 ns [82.3103, 82.862] (0.95 confidence) | 68.3414 ns [68.1145, 68.6044] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;batch/translated/with&#95;construction | 86.276 ns [85.9826, 86.5867] (0.95 confidence) | 71.5975 ns [71.2435, 71.9914] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/cancellation/difference | 17.553 ns [17.3493, 17.8927] (0.95 confidence) | 11.1764 ns [11.0498, 11.3151] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/cancellation/difference&#95;endpoints | 18.6115 ns [18.5738, 18.6468] (0.95 confidence) | 11.8862 ns [11.8594, 11.9134] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/cancellation/dot | 7.01298 ns [6.96746, 7.05707] (0.95 confidence) | 6.05106 ns [6.0411, 6.05951] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/cancellation/dot&#95;endpoints | 7.91425 ns [7.89235, 7.93777] (0.95 confidence) | 6.69282 ns [6.64551, 6.74723] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/dense/difference | 23.6575 ns [23.6132, 23.7025] (0.95 confidence) | 14.2959 ns [14.2659, 14.3396] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/dense/difference&#95;endpoints | 26.9968 ns [26.9339, 27.0582] (0.95 confidence) | 16.6533 ns [16.5494, 16.769] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/dense/dot | 9.62576 ns [9.60793, 9.64382] (0.95 confidence) | 8.045 ns [8.00619, 8.0872] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/dense/dot&#95;endpoints | 10.4202 ns [10.3931, 10.4479] (0.95 confidence) | 8.8761 ns [8.78278, 8.97174] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/dense&#95;inexact/difference | 23.7226 ns [23.6823, 23.7629] (0.95 confidence) | 14.4558 ns [14.3704, 14.553] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/dense&#95;inexact/difference&#95;endpoints | 25.3342 ns [25.2789, 25.392] (0.95 confidence) | 16.2561 ns [16.1449, 16.3753] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/dense&#95;inexact/dot | 9.54145 ns [9.50886, 9.57517] (0.95 confidence) | 7.96981 ns [7.94303, 8.00144] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/dense&#95;inexact/dot&#95;endpoints | 11.2277 ns [11.1996, 11.257] (0.95 confidence) | 9.25224 ns [9.21928, 9.28712] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/overflow/dot | 3.09871 ns [3.09107, 3.10662] (0.95 confidence) | 1.00235 ns [1.00016, 1.00523] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/sparse/difference | 8.45972 ns [8.45233, 8.46792] (0.95 confidence) | 5.4623 ns [5.45595, 5.46993] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/sparse/difference&#95;endpoints | 9.20638 ns [9.19563, 9.21664] (0.95 confidence) | 6.37373 ns [6.36515, 6.38194] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/sparse/dot | 2.56394 ns [2.54876, 2.58007] (0.95 confidence) | 2.86103 ns [2.83441, 2.88855] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/sparse/dot&#95;endpoints | 2.79407 ns [2.78223, 2.80756] (0.95 confidence) | 2.80918 ns [2.78919, 2.83018] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/underflow/difference | 7.67836 ns [7.67164, 7.68548] (0.95 confidence) | 4.71338 ns [4.70504, 4.72388] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/underflow/difference&#95;endpoints | 7.75033 ns [7.73856, 7.76453] (0.95 confidence) | 4.96771 ns [4.96646, 4.96903] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/underflow/dot | 3.45187 ns [3.44527, 3.45927] (0.95 confidence) | 3.14983 ns [3.1442, 3.15585] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/underflow/dot&#95;endpoints | 3.47106 ns [3.46641, 3.47585] (0.95 confidence) | 3.15457 ns [3.14951, 3.15959] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4/dot&#95;bounded&#95;inconclusive | 7.25319 ns [7.22225, 7.289] (0.95 confidence) | 6.32547 ns [6.26036, 6.39664] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4/dot&#95;bounded&#95;well&#95;separated | 11.4806 ns [11.4215, 11.5451] (0.95 confidence) | 9.39528 ns [9.34615, 9.45374] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4/dot&#95;difference&#95;bounded&#95;inconclusive | 31.3311 ns [31.2725, 31.3874] (0.95 confidence) | 21.0049 ns [20.8484, 21.1689] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4/dot&#95;difference&#95;bounded&#95;well&#95;separated | 27.8096 ns [27.4429, 28.3147] (0.95 confidence) | 16.4615 ns [16.3751, 16.5646] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4/dot&#95;plain&#95;well&#95;separated | 0.681686 ns [0.680637, 0.682756] (0.95 confidence) | 0.681934 ns [0.680306, 0.683653] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;batch/origin/construction | 2.62783 ns [2.62183, 2.6337] (0.95 confidence) | 2.34742 ns [2.34093, 2.35397] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;batch/origin/prepared | 40.2892 ns [40.222, 40.361] (0.95 confidence) | 34.9296 ns [34.8214, 35.0331] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;batch/origin/with&#95;construction | 46.8609 ns [46.809, 46.9119] (0.95 confidence) | 42.0714 ns [41.8533, 42.2992] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;batch/translated/construction | 2.94963 ns [2.94141, 2.95795] (0.95 confidence) | 2.34395 ns [2.3378, 2.35012] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;batch/translated/prepared | 132.947 ns [132.7, 133.188] (0.95 confidence) | 107.182 ns [106.958, 107.371] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;batch/translated/with&#95;construction | 139.312 ns [138.998, 139.64] (0.95 confidence) | 113.721 ns [113.254, 114.259] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/cancellation/difference | 18.184 ns [18.1441, 18.224] (0.95 confidence) | 11.7395 ns [11.6355, 11.8574] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/cancellation/difference&#95;endpoints | 19.8011 ns [19.7729, 19.8301] (0.95 confidence) | 12.9636 ns [12.8339, 13.105] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/cancellation/dot | 7.46986 ns [7.4468, 7.49647] (0.95 confidence) | 6.35969 ns [6.30525, 6.42148] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/cancellation/dot&#95;endpoints | 8.31456 ns [8.26713, 8.36457] (0.95 confidence) | 6.92672 ns [6.91861, 6.93443] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/dense/difference | 31.304 ns [31.2205, 31.4014] (0.95 confidence) | 18.5546 ns [18.4454, 18.6798] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/dense/difference&#95;endpoints | 35.0265 ns [34.9113, 35.1552] (0.95 confidence) | 20.7828 ns [20.7599, 20.8041] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/dense/dot | 11.4752 ns [11.4315, 11.5241] (0.95 confidence) | 9.61001 ns [9.54748, 9.68222] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/dense/dot&#95;endpoints | 13.7059 ns [13.5964, 13.8521] (0.95 confidence) | 11.1758 ns [11.1133, 11.2453] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/dense&#95;inexact/difference | 31.3867 ns [31.3106, 31.463] (0.95 confidence) | 18.8116 ns [18.5667, 19.2074] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/dense&#95;inexact/difference&#95;endpoints | 32.999 ns [32.905, 33.0971] (0.95 confidence) | 20.1724 ns [20.1092, 20.2472] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/dense&#95;inexact/dot | 11.4643 ns [11.4275, 11.5064] (0.95 confidence) | 9.87758 ns [9.78521, 9.97716] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/dense&#95;inexact/dot&#95;endpoints | 12.9519 ns [12.9285, 12.9753] (0.95 confidence) | 10.3229 ns [10.2486, 10.4069] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/overflow/dot | 3.02529 ns [3.01946, 3.03108] (0.95 confidence) | 0.997613 ns [0.996478, 0.998728] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/sparse/difference | 10.3817 ns [10.3563, 10.4088] (0.95 confidence) | 6.23896 ns [6.21215, 6.2709] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/sparse/difference&#95;endpoints | 11.0616 ns [11.0422, 11.0816] (0.95 confidence) | 7.40526 ns [7.39598, 7.41434] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/sparse/dot | 3.20807 ns [3.2027, 3.21319] (0.95 confidence) | 3.22407 ns [3.18888, 3.26249] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/sparse/dot&#95;endpoints | 3.40666 ns [3.40127, 3.41175] (0.95 confidence) | 3.23424 ns [3.22235, 3.24685] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/underflow/difference | 9.40378 ns [9.393, 9.41345] (0.95 confidence) | 5.68948 ns [5.68632, 5.69242] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/underflow/difference&#95;endpoints | 9.5079 ns [9.4925, 9.52197] (0.95 confidence) | 6.03071 ns [6.00265, 6.06138] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/underflow/dot | 4.59693 ns [4.58846, 4.60675] (0.95 confidence) | 3.86676 ns [3.85372, 3.88578] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/underflow/dot&#95;endpoints | 4.15747 ns [4.15295, 4.16173] (0.95 confidence) | 3.88086 ns [3.87257, 3.88954] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;batch/origin/construction | 3.08404 ns [3.07949, 3.08852] (0.95 confidence) | 2.87331 ns [2.86737, 2.87897] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;batch/origin/prepared | 55.1543 ns [55.0541, 55.2407] (0.95 confidence) | 46.2404 ns [46.1476, 46.3405] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;batch/origin/with&#95;construction | 66.5465 ns [66.4599, 66.63] (0.95 confidence) | 58.2986 ns [58.09, 58.5157] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;batch/translated/construction | 3.49129 ns [3.46948, 3.51421] (0.95 confidence) | 3.21946 ns [3.21347, 3.22532] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;batch/translated/prepared | 180.428 ns [180.097, 180.752] (0.95 confidence) | 147.144 ns [146.129, 148.235] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;batch/translated/with&#95;construction | 190.574 ns [189.911, 191.254] (0.95 confidence) | 155.894 ns [155.048, 156.863] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/cancellation/difference | 19.8848 ns [19.8568, 19.9123] (0.95 confidence) | 12.4082 ns [12.3571, 12.4652] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/cancellation/difference&#95;endpoints | 22.2038 ns [22.1674, 22.2395] (0.95 confidence) | 13.5516 ns [13.453, 13.66] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/cancellation/dot | 7.47088 ns [7.4484, 7.49469] (0.95 confidence) | 6.7673 ns [6.71645, 6.82334] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/cancellation/dot&#95;endpoints | 8.6089 ns [8.56588, 8.65811] (0.95 confidence) | 7.31332 ns [7.26061, 7.37154] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/dense/difference | 38.3278 ns [38.2533, 38.4001] (0.95 confidence) | 22.8999 ns [22.7578, 23.0611] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/dense/difference&#95;endpoints | 39.947 ns [39.8562, 40.0364] (0.95 confidence) | 24.1543 ns [24.0884, 24.2297] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/dense/dot | 14.0795 ns [14.0328, 14.1267] (0.95 confidence) | 11.8912 ns [11.7693, 12.0212] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/dense/dot&#95;endpoints | 15.2008 ns [15.1482, 15.256] (0.95 confidence) | 12.2635 ns [12.2118, 12.3215] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/dense&#95;inexact/difference | 38.5468 ns [38.4536, 38.6495] (0.95 confidence) | 22.8686 ns [22.7683, 22.985] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/dense&#95;inexact/difference&#95;endpoints | 43.7183 ns [43.6348, 43.7999] (0.95 confidence) | 25.4875 ns [25.345, 25.6411] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/dense&#95;inexact/dot | 14.152 ns [14.0829, 14.2227] (0.95 confidence) | 11.7909 ns [11.6989, 11.891] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/dense&#95;inexact/dot&#95;endpoints | 16.4904 ns [16.4557, 16.5245] (0.95 confidence) | 13.2503 ns [13.1726, 13.3357] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/overflow/dot | 2.99418 ns [2.98901, 2.99917] (0.95 confidence) | 0.998251 ns [0.99683, 0.999709] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/sparse/difference | 11.8847 ns [11.8674, 11.9015] (0.95 confidence) | 7.24318 ns [7.22331, 7.26597] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/sparse/difference&#95;endpoints | 12.8203 ns [12.8026, 12.837] (0.95 confidence) | 8.35253 ns [8.34043, 8.36453] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/sparse/dot | 4.10202 ns [4.09781, 4.10606] (0.95 confidence) | 3.68348 ns [3.66071, 3.7093] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/sparse/dot&#95;endpoints | 4.27448 ns [4.26883, 4.28012] (0.95 confidence) | 3.75049 ns [3.74274, 3.75912] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/underflow/difference | 11.131 ns [11.1185, 11.1439] (0.95 confidence) | 7.58527 ns [7.56604, 7.61009] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/underflow/difference&#95;endpoints | 11.2087 ns [11.1959, 11.2213] (0.95 confidence) | 7.85098 ns [7.83177, 7.8716] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/underflow/dot | 4.87682 ns [4.86607, 4.88872] (0.95 confidence) | 4.57564 ns [4.56414, 4.58644] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/underflow/dot&#95;endpoints | 4.87035 ns [4.86516, 4.87548] (0.95 confidence) | 4.58856 ns [4.57642, 4.60173] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;batch/origin/construction | 3.4638 ns [3.46059, 3.46704] (0.95 confidence) | 3.0836 ns [3.07961, 3.08752] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;batch/origin/prepared | 78.6625 ns [78.519, 78.7987] (0.95 confidence) | 60.7993 ns [60.6755, 60.9269] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;batch/origin/with&#95;construction | 101.492 ns [101.374, 101.606] (0.95 confidence) | 80.1084 ns [79.9421, 80.2678] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;batch/translated/construction | 3.61253 ns [3.60859, 3.61636] (0.95 confidence) | 3.54713 ns [3.54275, 3.55144] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;batch/translated/prepared | 253.802 ns [253.348, 254.235] (0.95 confidence) | 202.166 ns [201.317, 203.128] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;batch/translated/with&#95;construction | 276.935 ns [276.443, 277.437] (0.95 confidence) | 228.017 ns [226.688, 229.451] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/cancellation/difference | 22.6738 ns [22.6376, 22.708] (0.95 confidence) | 12.9897 ns [12.8771, 13.1071] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/cancellation/difference&#95;endpoints | 24.521 ns [24.4609, 24.5853] (0.95 confidence) | 14.1945 ns [14.1188, 14.2769] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/cancellation/dot | 8.0872 ns [8.05858, 8.11923] (0.95 confidence) | 6.92801 ns [6.89658, 6.96661] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/cancellation/dot&#95;endpoints | 9.49459 ns [9.48097, 9.50779] (0.95 confidence) | 7.71687 ns [7.6655, 7.77284] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/dense/difference | 48.0731 ns [47.9743, 48.1826] (0.95 confidence) | 28.6961 ns [28.6412, 28.7503] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/dense/difference&#95;endpoints | 53.5182 ns [53.4372, 53.5953] (0.95 confidence) | 32.1277 ns [32.0639, 32.1932] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/dense/dot | 16.5929 ns [16.5427, 16.6444] (0.95 confidence) | 13.9468 ns [13.836, 14.0646] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/dense/dot&#95;endpoints | 19.3972 ns [19.3446, 19.4458] (0.95 confidence) | 15.1317 ns [15.0375, 15.2413] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/dense&#95;inexact/difference | 48.4965 ns [48.0595, 49.1208] (0.95 confidence) | 28.6659 ns [28.6319, 28.6987] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/dense&#95;inexact/difference&#95;endpoints | 49.7462 ns [49.6587, 49.8286] (0.95 confidence) | 31.3564 ns [31.2913, 31.4186] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/dense&#95;inexact/dot | 16.801 ns [16.7451, 16.8586] (0.95 confidence) | 13.9913 ns [13.8882, 14.1021] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/dense&#95;inexact/dot&#95;endpoints | 19.2068 ns [19.1345, 19.2776] (0.95 confidence) | 15.1294 ns [15.0552, 15.2148] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/overflow/dot | 3.05078 ns [3.04394, 3.05748] (0.95 confidence) | 0.99979 ns [0.996842, 1.00392] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/sparse/difference | 13.6097 ns [13.5886, 13.6297] (0.95 confidence) | 8.07792 ns [8.06679, 8.08811] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/sparse/difference&#95;endpoints | 14.6416 ns [14.6179, 14.6653] (0.95 confidence) | 9.40877 ns [9.38724, 9.42891] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/sparse/dot | 5.01579 ns [5.00902, 5.02247] (0.95 confidence) | 4.09967 ns [4.09649, 4.10271] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/sparse/dot&#95;endpoints | 5.21585 ns [5.20648, 5.22573] (0.95 confidence) | 4.17883 ns [4.17511, 4.18235] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/underflow/difference | 12.6704 ns [12.5957, 12.7688] (0.95 confidence) | 8.57361 ns [8.55731, 8.59074] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/underflow/difference&#95;endpoints | 13.2603 ns [13.2178, 13.3048] (0.95 confidence) | 8.85114 ns [8.82632, 8.87718] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/underflow/dot | 5.62272 ns [5.61626, 5.62892] (0.95 confidence) | 5.29493 ns [5.28415, 5.30605] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/underflow/dot&#95;endpoints | 6.03973 ns [6.02924, 6.05008] (0.95 confidence) | 5.34642 ns [5.33946, 5.3532] (0.95 confidence) | — | — |
| polynomial/coincident&#95;lifted/order4/determinant | — | 13077.2 ns [13044.9, 13112.2] (0.95 confidence) | — | — |
| polynomial/coincident&#95;lifted/order5/determinant | — | 103607 ns [103373, 103843] (0.95 confidence) | — | — |
| polynomial/coincident&#95;lifted/order6/determinant | — | 895718 ns [894322, 897167] (0.95 confidence) | — | — |
| polynomial/coincident&#95;lifted/order7/determinant | — | 8.48317e+06 ns [8.46553e+06, 8.50083e+06] (0.95 confidence) | — | — |
| polynomial/coincident&#95;lifted/order8/determinant | — | 8.48955e+07 ns [8.47658e+07, 8.50265e+07] (0.95 confidence) | — | — |
| public&#95;d2/det&#95;bound | 1.95377 ns [1.95017, 1.95755] (0.95 confidence) | 1.95877 ns [1.9564, 1.96106] (0.95 confidence) | — | — |
| public&#95;d2/det&#95;direct | 1.12822 ns [1.12654, 1.12988] (0.95 confidence) | 1.12342 ns [1.12159, 1.12514] (0.95 confidence) | — | — |
| public&#95;d2/det&#95;with&#95;bound | 2.11925 ns [2.1147, 2.12386] (0.95 confidence) | 2.11587 ns [2.11194, 2.11965] (0.95 confidence) | — | — |
| public&#95;d2/first&#95;asymmetry | 1.37976 ns [1.37772, 1.38174] (0.95 confidence) | 1.37375 ns [1.372, 1.3755] (0.95 confidence) | — | — |
| public&#95;d2/interval&#95;det | 8.24584 ns [8.17254, 8.32106] (0.95 confidence) | 7.90585 ns [7.8554, 7.9622] (0.95 confidence) | — | — |
| public&#95;d2/interval&#95;matrix&#95;access | 6.28581 ns [6.2609, 6.31052] (0.95 confidence) | 6.23261 ns [6.21977, 6.24492] (0.95 confidence) | — | — |
| public&#95;d2/interval&#95;matrix&#95;from&#95;matrix | 2.5871 ns [2.58131, 2.59281] (0.95 confidence) | 2.56199 ns [2.55587, 2.56807] (0.95 confidence) | — | — |
| public&#95;d2/interval&#95;matrix&#95;from&#95;points | 2.34132 ns [2.33414, 2.34872] (0.95 confidence) | 2.24331 ns [2.23415, 2.25274] (0.95 confidence) | — | — |
| public&#95;d2/interval&#95;matrix&#95;from&#95;rows | 3.5445 ns [3.53872, 3.55028] (0.95 confidence) | 3.49968 ns [3.49319, 3.50609] (0.95 confidence) | — | — |
| public&#95;d2/interval&#95;matrix&#95;identity | 2.07515 ns [2.07153, 2.07863] (0.95 confidence) | 2.06627 ns [2.06328, 2.06927] (0.95 confidence) | — | — |
| public&#95;d2/interval&#95;matrix&#95;set | 3.77827 ns [3.76936, 3.78677] (0.95 confidence) | 3.721 ns [3.71436, 3.72734] (0.95 confidence) | — | — |
| public&#95;d2/interval&#95;matrix&#95;zero | 2.10014 ns [2.09579, 2.10426] (0.95 confidence) | 2.07286 ns [2.0696, 2.07607] (0.95 confidence) | — | — |
| public&#95;d2/is&#95;symmetric | 0.757044 ns [0.756134, 0.758055] (0.95 confidence) | 0.757249 ns [0.755479, 0.758665] (0.95 confidence) | — | — |
| public&#95;d2/matrix&#95;access | 4.67996 ns [4.66853, 4.69057] (0.95 confidence) | 4.68666 ns [4.67646, 4.69644] (0.95 confidence) | — | — |
| public&#95;d2/matrix&#95;from&#95;rows | 2.89635 ns [2.89037, 2.90502] (0.95 confidence) | 2.88372 ns [2.87858, 2.8884] (0.95 confidence) | — | — |
| public&#95;d2/matrix&#95;identity | 1.2318 ns [1.2296, 1.23393] (0.95 confidence) | 1.23409 ns [1.23184, 1.23625] (0.95 confidence) | — | — |
| public&#95;d2/matrix&#95;set | 2.08514 ns [2.07562, 2.0946] (0.95 confidence) | 2.08285 ns [2.07458, 2.09124] (0.95 confidence) | — | — |
| public&#95;d2/matrix&#95;zero | 1.15113 ns [1.1493, 1.15288] (0.95 confidence) | 1.14521 ns [1.14371, 1.14667] (0.95 confidence) | — | — |
| public&#95;d2/vector&#95;access | 1.60799 ns [1.60434, 1.61154] (0.95 confidence) | 1.5987 ns [1.5951, 1.60215] (0.95 confidence) | — | — |
| public&#95;d2/vector&#95;new | 1.87326 ns [1.8526, 1.8933] (0.95 confidence) | 1.88124 ns [1.86174, 1.89973] (0.95 confidence) | — | — |
| public&#95;d2/vector&#95;zero | 0.585701 ns [0.584702, 0.586653] (0.95 confidence) | 0.580976 ns [0.579954, 0.581982] (0.95 confidence) | — | — |
| public&#95;d3/det&#95;bound | 3.64868 ns [3.64506, 3.65256] (0.95 confidence) | 3.64765 ns [3.64337, 3.65162] (0.95 confidence) | — | — |
| public&#95;d3/det&#95;direct | 1.28478 ns [1.28401, 1.28557] (0.95 confidence) | 1.2638 ns [1.26285, 1.26486] (0.95 confidence) | — | — |
| public&#95;d3/det&#95;with&#95;bound | 3.67426 ns [3.67055, 3.67715] (0.95 confidence) | 3.66577 ns [3.66246, 3.66899] (0.95 confidence) | — | — |
| public&#95;d3/first&#95;asymmetry | 3.62519 ns [3.58524, 3.66582] (0.95 confidence) | 3.56803 ns [3.52965, 3.6072] (0.95 confidence) | — | — |
| public&#95;d3/interval&#95;det | 24.3271 ns [24.089, 24.5586] (0.95 confidence) | 23.5558 ns [23.3121, 23.787] (0.95 confidence) | — | — |
| public&#95;d3/interval&#95;matrix&#95;access | 9.75697 ns [9.73849, 9.77496] (0.95 confidence) | 9.75189 ns [9.73149, 9.77286] (0.95 confidence) | — | — |
| public&#95;d3/interval&#95;matrix&#95;from&#95;matrix | 5.14947 ns [5.13734, 5.16042] (0.95 confidence) | 5.0934 ns [5.08034, 5.10568] (0.95 confidence) | — | — |
| public&#95;d3/interval&#95;matrix&#95;from&#95;points | 6.37092 ns [6.35337, 6.38978] (0.95 confidence) | 6.27796 ns [6.26691, 6.28884] (0.95 confidence) | — | — |
| public&#95;d3/interval&#95;matrix&#95;from&#95;rows | 7.34223 ns [7.32765, 7.35656] (0.95 confidence) | 7.22097 ns [7.20736, 7.23456] (0.95 confidence) | — | — |
| public&#95;d3/interval&#95;matrix&#95;identity | 4.29289 ns [4.28582, 4.29978] (0.95 confidence) | 4.20375 ns [4.19791, 4.20952] (0.95 confidence) | — | — |
| public&#95;d3/interval&#95;matrix&#95;set | 7.79423 ns [7.77338, 7.81604] (0.95 confidence) | 7.70667 ns [7.68056, 7.73512] (0.95 confidence) | — | — |
| public&#95;d3/interval&#95;matrix&#95;zero | 4.30012 ns [4.27414, 4.33344] (0.95 confidence) | 4.21647 ns [4.21144, 4.2216] (0.95 confidence) | — | — |
| public&#95;d3/is&#95;symmetric | 1.33726 ns [1.33577, 1.33859] (0.95 confidence) | 1.49974 ns [1.49354, 1.50702] (0.95 confidence) | — | — |
| public&#95;d3/matrix&#95;access | 6.56973 ns [6.53061, 6.60652] (0.95 confidence) | 6.54463 ns [6.52762, 6.5614] (0.95 confidence) | — | — |
| public&#95;d3/matrix&#95;from&#95;rows | 5.73643 ns [5.72917, 5.74424] (0.95 confidence) | 5.65926 ns [5.64849, 5.67023] (0.95 confidence) | — | — |
| public&#95;d3/matrix&#95;identity | 2.67303 ns [2.66882, 2.67711] (0.95 confidence) | 2.65408 ns [2.64998, 2.65804] (0.95 confidence) | — | — |
| public&#95;d3/matrix&#95;set | 4.20847 ns [4.19335, 4.2237] (0.95 confidence) | 4.16062 ns [4.14534, 4.17607] (0.95 confidence) | — | — |
| public&#95;d3/matrix&#95;zero | 2.38886 ns [2.383, 2.39513] (0.95 confidence) | 2.37966 ns [2.37319, 2.38693] (0.95 confidence) | — | — |
| public&#95;d3/vector&#95;access | 2.07414 ns [2.06896, 2.07889] (0.95 confidence) | 2.04925 ns [2.04449, 2.05368] (0.95 confidence) | — | — |
| public&#95;d3/vector&#95;new | 3.96454 ns [3.93816, 3.98999] (0.95 confidence) | 3.64105 ns [3.63144, 3.65102] (0.95 confidence) | — | — |
| public&#95;d3/vector&#95;zero | 0.946373 ns [0.944706, 0.948002] (0.95 confidence) | 0.933342 ns [0.931479, 0.935158] (0.95 confidence) | — | — |
| public&#95;d4/det&#95;bound | 6.29051 ns [6.27166, 6.31055] (0.95 confidence) | 6.19631 ns [6.18009, 6.21318] (0.95 confidence) | — | — |
| public&#95;d4/det&#95;direct | 2.49293 ns [2.48989, 2.49659] (0.95 confidence) | 2.4714 ns [2.46855, 2.47434] (0.95 confidence) | — | — |
| public&#95;d4/det&#95;with&#95;bound | 6.27796 ns [6.26034, 6.29694] (0.95 confidence) | 6.20083 ns [6.18139, 6.22194] (0.95 confidence) | — | — |
| public&#95;d4/first&#95;asymmetry | 3.68036 ns [3.64569, 3.71685] (0.95 confidence) | 3.67914 ns [3.64128, 3.71732] (0.95 confidence) | — | — |
| public&#95;d4/interval&#95;det | 51.953 ns [51.5603, 52.3831] (0.95 confidence) | 51.6147 ns [51.1325, 52.1066] (0.95 confidence) | — | — |
| public&#95;d4/interval&#95;matrix&#95;access | 14.3958 ns [14.3428, 14.4489] (0.95 confidence) | 13.4976 ns [13.4583, 13.5371] (0.95 confidence) | — | — |
| public&#95;d4/interval&#95;matrix&#95;from&#95;matrix | 8.02359 ns [8.00723, 8.0397] (0.95 confidence) | 7.7925 ns [7.7788, 7.80565] (0.95 confidence) | — | — |
| public&#95;d4/interval&#95;matrix&#95;from&#95;points | 11.0625 ns [11.0441, 11.0825] (0.95 confidence) | 10.7524 ns [10.7315, 10.7714] (0.95 confidence) | — | — |
| public&#95;d4/interval&#95;matrix&#95;from&#95;rows | 11.4045 ns [11.3769, 11.4314] (0.95 confidence) | 10.9647 ns [10.9461, 10.9827] (0.95 confidence) | — | — |
| public&#95;d4/interval&#95;matrix&#95;identity | 6.63256 ns [6.61661, 6.64873] (0.95 confidence) | 6.34731 ns [6.33685, 6.35774] (0.95 confidence) | — | — |
| public&#95;d4/interval&#95;matrix&#95;set | 11.5971 ns [11.5746, 11.6193] (0.95 confidence) | 11.25 ns [11.2257, 11.2729] (0.95 confidence) | — | — |
| public&#95;d4/interval&#95;matrix&#95;zero | 6.68908 ns [6.67556, 6.70317] (0.95 confidence) | 6.46948 ns [6.46053, 6.47879] (0.95 confidence) | — | — |
| public&#95;d4/is&#95;symmetric | 2.124 ns [2.12078, 2.12823] (0.95 confidence) | 2.10252 ns [2.10147, 2.1035] (0.95 confidence) | — | — |
| public&#95;d4/matrix&#95;access | 9.28962 ns [9.26743, 9.31066] (0.95 confidence) | 9.22861 ns [9.20267, 9.25536] (0.95 confidence) | — | — |
| public&#95;d4/matrix&#95;from&#95;rows | 8.88423 ns [8.79975, 8.98852] (0.95 confidence) | 10.2323 ns [10.2163, 10.2459] (0.95 confidence) | — | — |
| public&#95;d4/matrix&#95;identity | 4.46969 ns [4.37108, 4.63253] (0.95 confidence) | 4.02802 ns [4.02058, 4.03575] (0.95 confidence) | — | — |
| public&#95;d4/matrix&#95;set | 7.57202 ns [7.54096, 7.60606] (0.95 confidence) | 7.20277 ns [7.1877, 7.21742] (0.95 confidence) | — | — |
| public&#95;d4/matrix&#95;zero | 3.97911 ns [3.96942, 3.98902] (0.95 confidence) | 3.78052 ns [3.77577, 3.78532] (0.95 confidence) | — | — |
| public&#95;d4/vector&#95;access | 2.65357 ns [2.64766, 2.6593] (0.95 confidence) | 2.56403 ns [2.55727, 2.57045] (0.95 confidence) | — | — |
| public&#95;d4/vector&#95;new | 2.9744 ns [2.96747, 2.98196] (0.95 confidence) | 2.82001 ns [2.8119, 2.82766] (0.95 confidence) | — | — |
| public&#95;d4/vector&#95;zero | 1.19069 ns [1.18828, 1.1931] (0.95 confidence) | 1.13964 ns [1.13765, 1.14152] (0.95 confidence) | — | — |
| public&#95;d5/det&#95;bound | 1.35472 ns [1.35206, 1.35749] (0.95 confidence) | 1.31596 ns [1.3123, 1.32064] (0.95 confidence) | — | — |
| public&#95;d5/det&#95;direct | 1.35495 ns [1.35279, 1.35711] (0.95 confidence) | 1.31277 ns [1.31181, 1.31393] (0.95 confidence) | — | — |
| public&#95;d5/det&#95;with&#95;bound | 0.862362 ns [0.860245, 0.864447] (0.95 confidence) | 0.835709 ns [0.833938, 0.837407] (0.95 confidence) | — | — |
| public&#95;d5/first&#95;asymmetry | 3.81202 ns [3.80838, 3.81611] (0.95 confidence) | 3.7817 ns [3.77808, 3.78602] (0.95 confidence) | — | — |
| public&#95;d5/interval&#95;det | 119.07 ns [118.404, 119.785] (0.95 confidence) | 122.559 ns [121.897, 123.212] (0.95 confidence) | — | — |
| public&#95;d5/interval&#95;matrix&#95;access | 20.1275 ns [19.9942, 20.2618] (0.95 confidence) | 18.2687 ns [18.2484, 18.2883] (0.95 confidence) | — | — |
| public&#95;d5/interval&#95;matrix&#95;from&#95;matrix | 10.7274 ns [10.69, 10.7655] (0.95 confidence) | 10.5108 ns [10.4753, 10.5461] (0.95 confidence) | — | — |
| public&#95;d5/interval&#95;matrix&#95;from&#95;points | 18.754 ns [18.7224, 18.7879] (0.95 confidence) | 18.2354 ns [18.1795, 18.2844] (0.95 confidence) | — | — |
| public&#95;d5/interval&#95;matrix&#95;from&#95;rows | 14.5922 ns [14.32, 14.866] (0.95 confidence) | 14.8194 ns [14.5344, 15.1006] (0.95 confidence) | — | — |
| public&#95;d5/interval&#95;matrix&#95;identity | 8.41685 ns [8.34958, 8.48387] (0.95 confidence) | 8.28581 ns [8.21927, 8.35171] (0.95 confidence) | — | — |
| public&#95;d5/interval&#95;matrix&#95;set | 14.0795 ns [13.8517, 14.3125] (0.95 confidence) | 14.9438 ns [14.6837, 15.1945] (0.95 confidence) | — | — |
| public&#95;d5/interval&#95;matrix&#95;zero | 9.45761 ns [9.44076, 9.47739] (0.95 confidence) | 9.23654 ns [9.22236, 9.25051] (0.95 confidence) | — | — |
| public&#95;d5/is&#95;symmetric | 3.71043 ns [3.70738, 3.7139] (0.95 confidence) | 3.68883 ns [3.68651, 3.69145] (0.95 confidence) | — | — |
| public&#95;d5/matrix&#95;access | 11.8187 ns [11.7635, 11.8747] (0.95 confidence) | 12.0413 ns [11.9964, 12.0859] (0.95 confidence) | — | — |
| public&#95;d5/matrix&#95;from&#95;rows | 11.7698 ns [11.7523, 11.7858] (0.95 confidence) | 11.6853 ns [11.6751, 11.6953] (0.95 confidence) | — | — |
| public&#95;d5/matrix&#95;identity | 6.94175 ns [6.90709, 6.98105] (0.95 confidence) | 6.65531 ns [6.64172, 6.6691] (0.95 confidence) | — | — |
| public&#95;d5/matrix&#95;set | 10.6188 ns [10.5963, 10.6418] (0.95 confidence) | 10.3227 ns [10.3019, 10.343] (0.95 confidence) | — | — |
| public&#95;d5/matrix&#95;zero | 5.81831 ns [5.80716, 5.82946] (0.95 confidence) | 5.60699 ns [5.59912, 5.61567] (0.95 confidence) | — | — |
| public&#95;d5/vector&#95;access | 3.11282 ns [3.10413, 3.12162] (0.95 confidence) | 3.03666 ns [3.03154, 3.04162] (0.95 confidence) | — | — |
| public&#95;d5/vector&#95;new | 4.85133 ns [4.83621, 4.86717] (0.95 confidence) | 4.82995 ns [4.81858, 4.84103] (0.95 confidence) | — | — |
| public&#95;d5/vector&#95;zero | 1.54095 ns [1.53076, 1.55098] (0.95 confidence) | 1.48448 ns [1.47417, 1.4947] (0.95 confidence) | — | — |
| public&#95;exact&#95;d2/matrix&#95;access | 0.788751 ns [0.782776, 0.795715] (0.95 confidence) | 0.752188 ns [0.749156, 0.755004] (0.95 confidence) | — | — |
| public&#95;exact&#95;d2/matrix&#95;from&#95;fn | 135.192 ns [134.902, 135.491] (0.95 confidence) | 130.341 ns [130.047, 130.605] (0.95 confidence) | — | — |
| public&#95;exact&#95;d2/matrix&#95;from&#95;rows | 138.559 ns [138.26, 138.869] (0.95 confidence) | 133.771 ns [133.529, 133.996] (0.95 confidence) | — | — |
| public&#95;exact&#95;d2/matrix&#95;into&#95;rows | 67.017 ns [66.8741, 67.1682] (0.95 confidence) | 64.3394 ns [64.218, 64.4604] (0.95 confidence) | — | — |
| public&#95;exact&#95;d2/matrix&#95;zero | 2.7553 ns [2.74835, 2.76366] (0.95 confidence) | 2.68555 ns [2.68287, 2.68814] (0.95 confidence) | — | — |
| public&#95;exact&#95;d2/vector&#95;access | 0.784764 ns [0.781496, 0.78821] (0.95 confidence) | 0.754558 ns [0.750821, 0.758999] (0.95 confidence) | — | — |
| public&#95;exact&#95;d2/vector&#95;from&#95;fn | 38.4462 ns [38.3354, 38.5491] (0.95 confidence) | 36.5897 ns [36.5376, 36.6365] (0.95 confidence) | — | — |
| public&#95;exact&#95;d2/vector&#95;into&#95;array | 16.1368 ns [16.0544, 16.2299] (0.95 confidence) | 15.3632 ns [15.3322, 15.3929] (0.95 confidence) | — | — |
| public&#95;exact&#95;d2/vector&#95;new | 36.745 ns [36.6787, 36.8161] (0.95 confidence) | 35.8566 ns [35.7999, 35.9088] (0.95 confidence) | — | — |
| public&#95;exact&#95;d2/vector&#95;zero | 1.37967 ns [1.37795, 1.38156] (0.95 confidence) | 1.35263 ns [1.35103, 1.35425] (0.95 confidence) | — | — |
| public&#95;exact&#95;d3/matrix&#95;access | 0.781708 ns [0.779328, 0.783924] (0.95 confidence) | 0.755792 ns [0.753803, 0.757724] (0.95 confidence) | — | — |
| public&#95;exact&#95;d3/matrix&#95;from&#95;fn | 232.995 ns [232.658, 233.326] (0.95 confidence) | 229.125 ns [228.464, 229.858] (0.95 confidence) | — | — |
| public&#95;exact&#95;d3/matrix&#95;from&#95;rows | 272.728 ns [272.302, 273.168] (0.95 confidence) | 265.195 ns [264.839, 265.534] (0.95 confidence) | — | — |
| public&#95;exact&#95;d3/matrix&#95;into&#95;rows | 125.867 ns [125.632, 126.085] (0.95 confidence) | 122.694 ns [122.489, 122.892] (0.95 confidence) | — | — |
| public&#95;exact&#95;d3/matrix&#95;zero | 6.3648 ns [6.3545, 6.37528] (0.95 confidence) | 6.15475 ns [6.14516, 6.16381] (0.95 confidence) | — | — |
| public&#95;exact&#95;d3/vector&#95;access | 0.779482 ns [0.777292, 0.781589] (0.95 confidence) | 0.763888 ns [0.75945, 0.769186] (0.95 confidence) | — | — |
| public&#95;exact&#95;d3/vector&#95;from&#95;fn | 52.8937 ns [52.8213, 52.9668] (0.95 confidence) | 51.7105 ns [51.6264, 51.7924] (0.95 confidence) | — | — |
| public&#95;exact&#95;d3/vector&#95;into&#95;array | 23.1618 ns [23.0647, 23.2657] (0.95 confidence) | 22.3022 ns [22.2614, 22.3408] (0.95 confidence) | — | — |
| public&#95;exact&#95;d3/vector&#95;new | 52.6085 ns [52.5095, 52.715] (0.95 confidence) | 51.4494 ns [51.3502, 51.5426] (0.95 confidence) | — | — |
| public&#95;exact&#95;d3/vector&#95;zero | 2.05879 ns [2.05606, 2.06169] (0.95 confidence) | 2.01641 ns [2.01364, 2.01901] (0.95 confidence) | — | — |
| public&#95;exact&#95;d4/matrix&#95;access | 0.782031 ns [0.779701, 0.784428] (0.95 confidence) | 0.75746 ns [0.755263, 0.75949] (0.95 confidence) | — | — |
| public&#95;exact&#95;d4/matrix&#95;from&#95;fn | 291.467 ns [290.974, 291.948] (0.95 confidence) | 283.462 ns [283.068, 283.841] (0.95 confidence) | — | — |
| public&#95;exact&#95;d4/matrix&#95;from&#95;rows | 352.411 ns [351.849, 353.001] (0.95 confidence) | 343.819 ns [343.365, 344.253] (0.95 confidence) | — | — |
| public&#95;exact&#95;d4/matrix&#95;into&#95;rows | 133.577 ns [133.381, 133.775] (0.95 confidence) | 129.618 ns [129.199, 130.021] (0.95 confidence) | — | — |
| public&#95;exact&#95;d4/matrix&#95;zero | 15.6443 ns [15.5515, 15.7381] (0.95 confidence) | 15.3827 ns [15.2906, 15.4735] (0.95 confidence) | — | — |
| public&#95;exact&#95;d4/vector&#95;access | 0.7813 ns [0.778935, 0.783587] (0.95 confidence) | 0.757103 ns [0.754651, 0.759269] (0.95 confidence) | — | — |
| public&#95;exact&#95;d4/vector&#95;from&#95;fn | 67.5729 ns [67.4477, 67.697] (0.95 confidence) | 65.708 ns [65.6143, 65.7986] (0.95 confidence) | — | — |
| public&#95;exact&#95;d4/vector&#95;into&#95;array | 30.479 ns [30.4233, 30.5366] (0.95 confidence) | 29.5455 ns [29.5022, 29.5872] (0.95 confidence) | — | — |
| public&#95;exact&#95;d4/vector&#95;new | 66.6948 ns [66.5769, 66.8157] (0.95 confidence) | 64.9294 ns [64.8419, 65.0073] (0.95 confidence) | — | — |
| public&#95;exact&#95;d4/vector&#95;zero | 3.482 ns [3.4604, 3.50373] (0.95 confidence) | 3.45487 ns [3.43426, 3.47449] (0.95 confidence) | — | — |
| public&#95;exact&#95;d5/matrix&#95;access | 0.781558 ns [0.778321, 0.785216] (0.95 confidence) | 0.757158 ns [0.754775, 0.759215] (0.95 confidence) | — | — |
| public&#95;exact&#95;d5/matrix&#95;from&#95;fn | 439.838 ns [438.307, 441.384] (0.95 confidence) | 437.999 ns [437.021, 438.964] (0.95 confidence) | — | — |
| public&#95;exact&#95;d5/matrix&#95;from&#95;rows | 635.388 ns [634.356, 636.476] (0.95 confidence) | 621.381 ns [620.503, 622.163] (0.95 confidence) | — | — |
| public&#95;exact&#95;d5/matrix&#95;into&#95;rows | 309.243 ns [308.516, 309.991] (0.95 confidence) | 300.884 ns [300.227, 301.533] (0.95 confidence) | — | — |
| public&#95;exact&#95;d5/matrix&#95;zero | 22.2433 ns [22.2003, 22.2878] (0.95 confidence) | 21.7327 ns [21.6822, 21.7795] (0.95 confidence) | — | — |
| public&#95;exact&#95;d5/vector&#95;access | 0.780676 ns [0.778649, 0.782664] (0.95 confidence) | 0.75976 ns [0.757263, 0.761969] (0.95 confidence) | — | — |
| public&#95;exact&#95;d5/vector&#95;from&#95;fn | 83.6678 ns [83.5454, 83.7898] (0.95 confidence) | 81.7707 ns [81.6511, 81.8899] (0.95 confidence) | — | — |
| public&#95;exact&#95;d5/vector&#95;into&#95;array | 38.8157 ns [38.7492, 38.8841] (0.95 confidence) | 38.0498 ns [37.9746, 38.126] (0.95 confidence) | — | — |
| public&#95;exact&#95;d5/vector&#95;new | 83.62 ns [83.5143, 83.7261] (0.95 confidence) | 81.9606 ns [81.7969, 82.1143] (0.95 confidence) | — | — |
| public&#95;exact&#95;d5/vector&#95;zero | 4.29848 ns [4.27095, 4.32711] (0.95 confidence) | 4.19785 ns [4.16831, 4.22765] (0.95 confidence) | — | — |
| public&#95;exact&#95;scalar/rounded | 5.61009 ns [5.60404, 5.61654] (0.95 confidence) | 5.55237 ns [5.54787, 5.55668] (0.95 confidence) | — | — |
| public&#95;exact&#95;scalar/strict | 4.30833 ns [4.28897, 4.32931] (0.95 confidence) | 4.2207 ns [4.20439, 4.2384] (0.95 confidence) | — | — |
| public&#95;scalar/interval&#95;access | 1.57792 ns [1.57419, 1.58125] (0.95 confidence) | 1.57165 ns [1.56801, 1.57489] (0.95 confidence) | — | — |
| public&#95;scalar/interval&#95;negate | 0.729915 ns [0.729162, 0.730646] (0.95 confidence) | 0.72489 ns [0.723997, 0.725708] (0.95 confidence) | — | — |
| public&#95;scalar/interval&#95;new | 0.892185 ns [0.891272, 0.893047] (0.95 confidence) | 0.890087 ns [0.888923, 0.891185] (0.95 confidence) | — | — |
| public&#95;scalar/interval&#95;point | 0.582883 ns [0.580283, 0.586448] (0.95 confidence) | 0.578593 ns [0.576767, 0.580794] (0.95 confidence) | — | — |
| public&#95;scalar/interval&#95;subtraction | 1.85768 ns [1.85506, 1.86032] (0.95 confidence) | 1.84809 ns [1.8455, 1.85035] (0.95 confidence) | — | — |
| public&#95;scalar/tolerance | 0.96891 ns [0.965305, 0.971955] (0.95 confidence) | 0.961973 ns [0.958009, 0.965532] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d14/prepared | — | 1782.07 ns [1774.67, 1789.73] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d14/slice | — | 627.702 ns [626.348, 628.935] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d2/prepared | — | 256.954 ns [256.616, 257.37] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d2/slice | — | 96.4336 ns [96.2242, 96.6708] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d22/prepared | — | 2799.49 ns [2787.93, 2811.45] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d22/slice | — | 976.336 ns [974.764, 977.899] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d3/prepared | — | 380.256 ns [379.518, 380.994] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d3/slice | — | 125.716 ns [125.503, 125.932] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d4/prepared | — | 509.87 ns [508.072, 511.821] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d4/slice | — | 185.024 ns [184.702, 185.354] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d5/prepared | — | 627.588 ns [625.367, 629.837] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d5/slice | — | 210.951 ns [210.531, 211.373] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d6/prepared | — | 775.32 ns [771.482, 779.114] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d6/slice | — | 272.502 ns [271.549, 273.497] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d7/prepared | — | 874.139 ns [871.28, 877.078] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d7/slice | — | 299.145 ns [298.679, 299.623] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d8/prepared | — | 1019.75 ns [1014.98, 1024.74] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d8/slice | — | 356.394 ns [355.79, 356.999] (0.95 confidence) | — | — |
| rational&#95;dot/large/d14/prepared | — | 2836.65 ns [2832.67, 2840.5] (0.95 confidence) | — | — |
| rational&#95;dot/large/d14/slice | — | 1247.48 ns [1244.1, 1251] (0.95 confidence) | — | — |
| rational&#95;dot/large/d2/prepared | — | 333.795 ns [333.185, 334.411] (0.95 confidence) | — | — |
| rational&#95;dot/large/d2/slice | — | 122.064 ns [121.836, 122.296] (0.95 confidence) | — | — |
| rational&#95;dot/large/d22/prepared | — | 4536.8 ns [4528.67, 4545.16] (0.95 confidence) | — | — |
| rational&#95;dot/large/d22/slice | — | 2040.56 ns [2029.8, 2052.64] (0.95 confidence) | — | — |
| rational&#95;dot/large/d3/prepared | — | 535.532 ns [534.118, 536.976] (0.95 confidence) | — | — |
| rational&#95;dot/large/d3/slice | — | 208.547 ns [208.237, 208.866] (0.95 confidence) | — | — |
| rational&#95;dot/large/d4/prepared | — | 716.56 ns [715.061, 718.093] (0.95 confidence) | — | — |
| rational&#95;dot/large/d4/slice | — | 291.827 ns [290.475, 293.227] (0.95 confidence) | — | — |
| rational&#95;dot/large/d5/prepared | — | 915.709 ns [913.931, 917.54] (0.95 confidence) | — | — |
| rational&#95;dot/large/d5/slice | — | 371.873 ns [370.945, 372.843] (0.95 confidence) | — | — |
| rational&#95;dot/large/d6/prepared | — | 1117.91 ns [1116.23, 1119.6] (0.95 confidence) | — | — |
| rational&#95;dot/large/d6/slice | — | 467.657 ns [466.175, 469.24] (0.95 confidence) | — | — |
| rational&#95;dot/large/d7/prepared | — | 1329.37 ns [1326.95, 1331.85] (0.95 confidence) | — | — |
| rational&#95;dot/large/d7/slice | — | 568.946 ns [567.925, 569.997] (0.95 confidence) | — | — |
| rational&#95;dot/large/d8/prepared | — | 1523.99 ns [1521.09, 1526.96] (0.95 confidence) | — | — |
| rational&#95;dot/large/d8/slice | — | 643.316 ns [642.262, 644.382] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d14/prepared | — | 2375.69 ns [2372.27, 2378.9] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d14/slice | — | 1028.29 ns [1025.82, 1030.83] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d2/prepared | — | 271.552 ns [270.965, 272.258] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d2/slice | — | 92.2384 ns [92.112, 92.3666] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d22/prepared | — | 3809.82 ns [3803.93, 3815.7] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d22/slice | — | 1688.88 ns [1684.11, 1693.79] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d3/prepared | — | 438.351 ns [437.841, 438.843] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d3/slice | — | 168.653 ns [168.248, 169.061] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d4/prepared | — | 588.689 ns [587.938, 589.343] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d4/slice | — | 226.249 ns [225.84, 226.664] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d5/prepared | — | 758.046 ns [757.131, 758.969] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d5/slice | — | 300.054 ns [299.484, 300.62] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d6/prepared | — | 924.641 ns [923.509, 925.763] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d6/slice | — | 376.718 ns [375.856, 377.597] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d7/prepared | — | 1109.81 ns [1108.05, 1111.61] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d7/slice | — | 475.087 ns [473.891, 476.209] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d8/prepared | — | 1261.19 ns [1259.42, 1263.05] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d8/slice | — | 530.164 ns [529.108, 531.199] (0.95 confidence) | — | — |
| rational&#95;input&#95;d2/det&#95;row&#95;cleared&#95;bareiss | 179.021 ns [178.66, 179.408] (0.95 confidence) | 178.581 ns [178.371, 178.791] (0.95 confidence) | — | — |
| rational&#95;input&#95;d2/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 77.3422 ns [77.1896, 77.5057] (0.95 confidence) | 76.6864 ns [76.5773, 76.7956] (0.95 confidence) | — | — |
| rational&#95;input&#95;d2/solve&#95;row&#95;cleared&#95;bareiss | 615.768 ns [615.16, 616.442] (0.95 confidence) | 623.588 ns [622.768, 624.438] (0.95 confidence) | — | — |
| rational&#95;input&#95;d3/det&#95;row&#95;cleared&#95;bareiss | 569.289 ns [568.31, 570.319] (0.95 confidence) | 590.315 ns [588.762, 591.851] (0.95 confidence) | — | — |
| rational&#95;input&#95;d3/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 356.333 ns [354.904, 357.892] (0.95 confidence) | 389.615 ns [388.381, 390.833] (0.95 confidence) | — | — |
| rational&#95;input&#95;d3/solve&#95;row&#95;cleared&#95;bareiss | 1230.88 ns [1228.96, 1232.9] (0.95 confidence) | 1237.92 ns [1236.24, 1239.57] (0.95 confidence) | — | — |
| rational&#95;input&#95;d4/det&#95;row&#95;cleared&#95;bareiss | 1161.49 ns [1159.51, 1163.88] (0.95 confidence) | 1219.64 ns [1218.01, 1221.15] (0.95 confidence) | — | — |
| rational&#95;input&#95;d4/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 761.912 ns [760.419, 763.448] (0.95 confidence) | 824.138 ns [820.924, 827.121] (0.95 confidence) | — | — |
| rational&#95;input&#95;d4/solve&#95;row&#95;cleared&#95;bareiss | 2284.33 ns [2249.53, 2322.09] (0.95 confidence) | 2346.93 ns [2323.21, 2370.72] (0.95 confidence) | — | — |
| rational&#95;input&#95;d5/det&#95;row&#95;cleared&#95;bareiss | 2051.97 ns [2049.43, 2054.59] (0.95 confidence) | 2069.46 ns [2065.43, 2073.61] (0.95 confidence) | — | — |
| rational&#95;input&#95;d5/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 1534.56 ns [1531.6, 1537.64] (0.95 confidence) | 1569.19 ns [1565.76, 1572.85] (0.95 confidence) | — | — |
| rational&#95;input&#95;d5/solve&#95;row&#95;cleared&#95;bareiss | 3351.04 ns [3343.37, 3360.12] (0.95 confidence) | 3342.48 ns [3338.6, 3346.4] (0.95 confidence) | — | — |
| rational&#95;input&#95;d6/det&#95;row&#95;cleared&#95;bareiss | 3510.87 ns [3505.63, 3516.44] (0.95 confidence) | 3469.89 ns [3464.93, 3475.14] (0.95 confidence) | — | — |
| rational&#95;input&#95;d6/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 2842.43 ns [2837.84, 2847.38] (0.95 confidence) | 2855.33 ns [2849.98, 2860.86] (0.95 confidence) | — | — |
| rational&#95;input&#95;d6/solve&#95;row&#95;cleared&#95;bareiss | 5763.7 ns [5752.46, 5775.5] (0.95 confidence) | 5762.91 ns [5753.18, 5772.59] (0.95 confidence) | — | — |
| rational&#95;input&#95;d7/det&#95;row&#95;cleared&#95;bareiss | 5529.53 ns [5522.44, 5536.89] (0.95 confidence) | 5565.9 ns [5555.18, 5577.63] (0.95 confidence) | — | — |
| rational&#95;input&#95;d7/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 4657.58 ns [4651.71, 4663.59] (0.95 confidence) | 4671.96 ns [4666.68, 4677.2] (0.95 confidence) | — | — |
| rational&#95;input&#95;d7/solve&#95;row&#95;cleared&#95;bareiss | 9983.3 ns [9967.02, 10000] (0.95 confidence) | 9918.07 ns [9902.09, 9933.77] (0.95 confidence) | — | — |
| rational&#95;input&#95;d8/det&#95;row&#95;cleared&#95;bareiss | 8028.54 ns [8015.57, 8042.33] (0.95 confidence) | 8094.75 ns [8082.37, 8107.21] (0.95 confidence) | — | — |
| rational&#95;input&#95;d8/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 7181.85 ns [7169.31, 7195.37] (0.95 confidence) | 7290.74 ns [7279.22, 7301.85] (0.95 confidence) | — | — |
| rational&#95;input&#95;d8/solve&#95;row&#95;cleared&#95;bareiss | 15199.2 ns [15171.8, 15228.7] (0.95 confidence) | 14946.5 ns [14924.4, 14968.9] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d2/det&#95;row&#95;cleared&#95;bareiss | 42286.6 ns [42223.2, 42351.7] (0.95 confidence) | 42055.2 ns [41970.2, 42149.4] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d2/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 952.135 ns [951.064, 953.336] (0.95 confidence) | 949.878 ns [947.547, 952.438] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d2/solve&#95;row&#95;cleared&#95;bareiss | 3055.23 ns [3052.17, 3058.71] (0.95 confidence) | 3057.7 ns [3054.16, 3061.09] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d3/det&#95;row&#95;cleared&#95;bareiss | 70075.1 ns [69949.3, 70228.4] (0.95 confidence) | 68976.6 ns [68906.2, 69047.7] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d3/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 4210.55 ns [4205.83, 4215.53] (0.95 confidence) | 4185.57 ns [4179.56, 4191.7] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d3/solve&#95;row&#95;cleared&#95;bareiss | 11376 ns [11364.5, 11388.8] (0.95 confidence) | 11355.8 ns [11339.2, 11371.7] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d4/det&#95;row&#95;cleared&#95;bareiss | 111003 ns [110886, 111133] (0.95 confidence) | 110119 ns [110006, 110232] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d4/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 9466.42 ns [9443.74, 9493.38] (0.95 confidence) | 9362.49 ns [9346.62, 9378.72] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d4/solve&#95;row&#95;cleared&#95;bareiss | 33727.1 ns [33664.5, 33792.7] (0.95 confidence) | 33538.8 ns [33485, 33587] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d5/det&#95;row&#95;cleared&#95;bareiss | 187316 ns [187016, 187660] (0.95 confidence) | 186093 ns [185869, 186318] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d5/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 43043.2 ns [42989.5, 43100.9] (0.95 confidence) | 42398.6 ns [42350.1, 42442.6] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d5/solve&#95;row&#95;cleared&#95;bareiss | 69474.6 ns [69372.9, 69593.8] (0.95 confidence) | 67936.2 ns [67811.7, 68041.2] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d6/det&#95;row&#95;cleared&#95;bareiss | 302605 ns [302187, 303129] (0.95 confidence) | 299734 ns [299394, 300043] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d6/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 104294 ns [104157, 104435] (0.95 confidence) | 102496 ns [102346, 102633] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d6/solve&#95;row&#95;cleared&#95;bareiss | 161176 ns [160935, 161413] (0.95 confidence) | 158830 ns [158651, 158999] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d7/det&#95;row&#95;cleared&#95;bareiss | 503808 ns [503227, 504450] (0.95 confidence) | 496800 ns [496134, 497430] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d7/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 233499 ns [233152, 233869] (0.95 confidence) | 228252 ns [227833, 228608] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d7/solve&#95;row&#95;cleared&#95;bareiss | 339463 ns [339005, 339930] (0.95 confidence) | 331513 ns [330905, 332033] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d8/det&#95;row&#95;cleared&#95;bareiss | 822453 ns [821779, 823156] (0.95 confidence) | 811834 ns [810897, 812666] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d8/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 468132 ns [467392, 468867] (0.95 confidence) | 459388 ns [458703, 460068] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d8/solve&#95;row&#95;cleared&#95;bareiss | 650462 ns [649605, 651383] (0.95 confidence) | 636508 ns [635713, 637278] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d2/det&#95;row&#95;cleared&#95;bareiss | 9023.99 ns [8986.97, 9062.18] (0.95 confidence) | 8796.76 ns [8771.33, 8821.79] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d2/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 652.636 ns [651.155, 654.121] (0.95 confidence) | 631.58 ns [629.9, 633.434] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d2/solve&#95;row&#95;cleared&#95;bareiss | 2078.47 ns [2074.95, 2081.99] (0.95 confidence) | 2072.93 ns [2069.07, 2077.02] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d3/det&#95;row&#95;cleared&#95;bareiss | 13614.3 ns [13587.6, 13641.7] (0.95 confidence) | 13210.9 ns [13175.8, 13246.9] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d3/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 2794.49 ns [2787.32, 2801.92] (0.95 confidence) | 2719.98 ns [2711.9, 2728.38] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d3/solve&#95;row&#95;cleared&#95;bareiss | 5748.61 ns [5732.61, 5764.76] (0.95 confidence) | 5719.01 ns [5705.73, 5732.78] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d4/det&#95;row&#95;cleared&#95;bareiss | 22321.4 ns [22287.3, 22356] (0.95 confidence) | 21713.4 ns [21674.6, 21753.1] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d4/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 5978.67 ns [5957.97, 5999.87] (0.95 confidence) | 5876 ns [5856, 5896.28] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d4/solve&#95;row&#95;cleared&#95;bareiss | 13233 ns [13181.8, 13285.6] (0.95 confidence) | 12586 ns [12553.4, 12618.8] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d5/det&#95;row&#95;cleared&#95;bareiss | 33488.9 ns [33449, 33534.5] (0.95 confidence) | 33406.5 ns [33337.1, 33479.2] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d5/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 13332.9 ns [13301.4, 13365.8] (0.95 confidence) | 13074.3 ns [13042.2, 13107.9] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d5/solve&#95;row&#95;cleared&#95;bareiss | 20687.4 ns [20639.3, 20738.1] (0.95 confidence) | 20372 ns [20338.1, 20406.3] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d6/det&#95;row&#95;cleared&#95;bareiss | 52074.2 ns [51996.4, 52155.5] (0.95 confidence) | 51667.7 ns [51582.7, 51753.6] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d6/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 26377.7 ns [26312.3, 26444.1] (0.95 confidence) | 26002.1 ns [25945.7, 26061.7] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d6/solve&#95;row&#95;cleared&#95;bareiss | 39617.9 ns [39541.8, 39695.6] (0.95 confidence) | 38556.9 ns [38478.8, 38636.6] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d7/det&#95;row&#95;cleared&#95;bareiss | 83405.8 ns [83284.1, 83539.7] (0.95 confidence) | 82020.9 ns [81926.5, 82108.3] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d7/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 50605.6 ns [50541.2, 50674.3] (0.95 confidence) | 49487.1 ns [49330.2, 49640.5] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d7/solve&#95;row&#95;cleared&#95;bareiss | 72278.6 ns [72172.7, 72390.7] (0.95 confidence) | 70552.1 ns [70432.2, 70673.1] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d8/det&#95;row&#95;cleared&#95;bareiss | 132094 ns [131952, 132241] (0.95 confidence) | 129741 ns [129585, 129891] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d8/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 91300.8 ns [91157.8, 91455.3] (0.95 confidence) | 88819 ns [88636.2, 89007.7] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d8/solve&#95;row&#95;cleared&#95;bareiss | 126848 ns [126705, 126998] (0.95 confidence) | 125410 ns [125147, 125691] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d2/adapter | 1053.52 ns [1043.79, 1063.51] (0.95 confidence) | 1090.71 ns [1072.56, 1109.8] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d2/construct | 290.527 ns [289.937, 291.148] (0.95 confidence) | 286.3 ns [285.774, 286.828] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d2/prepared | 637.998 ns [636.894, 639.14] (0.95 confidence) | 647.728 ns [646.166, 649.36] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d3/adapter | 1982.36 ns [1973.72, 1993.68] (0.95 confidence) | 1886.8 ns [1881.92, 1890.89] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d3/construct | 642.408 ns [641.141, 643.737] (0.95 confidence) | 632.996 ns [631.922, 634.066] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d3/prepared | 1177.42 ns [1174.81, 1180.23] (0.95 confidence) | 1146.27 ns [1144.34, 1148.23] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d4/adapter | 3389.75 ns [3382.88, 3396.67] (0.95 confidence) | 3318.92 ns [3314.34, 3323.5] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d4/construct | 1060.79 ns [1057.91, 1063.7] (0.95 confidence) | 1037.72 ns [1035.83, 1039.57] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d4/prepared | 2056.24 ns [2050.99, 2061.75] (0.95 confidence) | 2078.89 ns [2073.07, 2084.75] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d5/adapter | 5566.53 ns [5555.86, 5578.45] (0.95 confidence) | 5397.14 ns [5381.71, 5410.37] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d5/construct | 1969.13 ns [1962.5, 1976.07] (0.95 confidence) | 1920.27 ns [1915.45, 1925.32] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d5/prepared | 3287.11 ns [3276.63, 3297.23] (0.95 confidence) | 3136.47 ns [3130.51, 3142.51] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d6/adapter | 7787.09 ns [7765.59, 7812.06] (0.95 confidence) | 7537.04 ns [7526.68, 7546.86] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d6/construct | 2547.16 ns [2541.4, 2552.95] (0.95 confidence) | 2494.65 ns [2489.39, 2500.22] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d6/prepared | 4723.52 ns [4712.29, 4734.82] (0.95 confidence) | 4645.83 ns [4639.66, 4651.96] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d7/adapter | 10671.6 ns [10657.1, 10686.5] (0.95 confidence) | 10237.4 ns [10222.9, 10251.7] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d7/construct | 3497.2 ns [3485.69, 3509.38] (0.95 confidence) | 3426.81 ns [3421.91, 3431.77] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d7/prepared | 6559.42 ns [6542.64, 6577.18] (0.95 confidence) | 6437.09 ns [6424.01, 6450.23] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d8/adapter | 13868.4 ns [13848.4, 13888.8] (0.95 confidence) | 13278.9 ns [13261.2, 13296] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d8/construct | 4106.44 ns [4099.78, 4113.05] (0.95 confidence) | 4006.46 ns [3999.35, 4013.36] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d8/prepared | 9620.86 ns [9605.13, 9637.05] (0.95 confidence) | 8973.52 ns [8956.33, 8990.93] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d2/adapter | 1052.64 ns [1047.28, 1058.34] (0.95 confidence) | 1114.62 ns [1097.81, 1131.72] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d2/construct | 331.181 ns [330.13, 332.268] (0.95 confidence) | 323.581 ns [322.614, 324.471] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d2/prepared | 615.312 ns [614.344, 616.373] (0.95 confidence) | 623.349 ns [622.294, 624.327] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d3/adapter | 2096.02 ns [2093.7, 2098.39] (0.95 confidence) | 2075.44 ns [2072.21, 2078.25] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d3/construct | 683.25 ns [680.501, 685.867] (0.95 confidence) | 658.492 ns [655.828, 660.999] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d3/prepared | 1279.57 ns [1277.55, 1281.7] (0.95 confidence) | 1276.1 ns [1274.69, 1277.47] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d4/adapter | 3504.54 ns [3497.83, 3511.79] (0.95 confidence) | 3472.18 ns [3467.15, 3477.11] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d4/construct | 1151.17 ns [1147.96, 1154.3] (0.95 confidence) | 1121.98 ns [1117.45, 1126.42] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d4/prepared | 2161.39 ns [2158.98, 2163.81] (0.95 confidence) | 2151.52 ns [2145.73, 2156.45] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d5/adapter | 5629.61 ns [5621.04, 5638.74] (0.95 confidence) | 5521.78 ns [5512.95, 5530.71] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d5/construct | 1932.17 ns [1927.81, 1936.43] (0.95 confidence) | 1891.74 ns [1883.43, 1900.54] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d5/prepared | 3437.7 ns [3432.49, 3442.94] (0.95 confidence) | 3319.55 ns [3315.57, 3323.42] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d6/adapter | 7993.21 ns [7980.63, 8006.28] (0.95 confidence) | 7861 ns [7844.13, 7878.14] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d6/construct | 2796.04 ns [2786.72, 2805.33] (0.95 confidence) | 2689.39 ns [2673.71, 2704.75] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d6/prepared | 4927.14 ns [4918.18, 4936.65] (0.95 confidence) | 4808.87 ns [4802.73, 4814.84] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d7/adapter | 11798.4 ns [11776.6, 11821.8] (0.95 confidence) | 11589.4 ns [11565.6, 11613.1] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d7/construct | 3951.25 ns [3937.93, 3964.86] (0.95 confidence) | 3849.28 ns [3832.08, 3865.97] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d7/prepared | 7401.28 ns [7390.39, 7412.67] (0.95 confidence) | 7376.07 ns [7359.4, 7393.38] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d8/adapter | 16390.4 ns [16369.5, 16412.8] (0.95 confidence) | 15830.9 ns [15805.5, 15854.9] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d8/construct | 5234.83 ns [5213.97, 5255.17] (0.95 confidence) | 5019.94 ns [4983.76, 5055.46] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d8/prepared | 10541.9 ns [10527.8, 10556.6] (0.95 confidence) | 10433.3 ns [10412.8, 10453.9] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d2/adapter | 458.956 ns [453.25, 464.77] (0.95 confidence) | 429.652 ns [427.615, 431.851] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d2/construct | 227.557 ns [227.259, 227.875] (0.95 confidence) | 223.365 ns [222.893, 223.883] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d2/prepared | 152.894 ns [152.663, 153.134] (0.95 confidence) | 151.109 ns [150.813, 151.401] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d3/adapter | 737.293 ns [736.411, 738.25] (0.95 confidence) | 727.58 ns [726.412, 728.784] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d3/construct | 413.412 ns [412.8, 414.043] (0.95 confidence) | 404.426 ns [403.599, 405.245] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d3/prepared | 250.742 ns [250.294, 251.245] (0.95 confidence) | 247.148 ns [246.725, 247.572] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d4/adapter | 1033.21 ns [1031.93, 1034.51] (0.95 confidence) | 1011.56 ns [1010.1, 1013.02] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d4/construct | 611.902 ns [610.784, 613.054] (0.95 confidence) | 595.13 ns [593.492, 596.799] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d4/prepared | 323.025 ns [322.586, 323.458] (0.95 confidence) | 321.741 ns [321.285, 322.191] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d5/adapter | 1387.24 ns [1385.16, 1389.35] (0.95 confidence) | 1354.39 ns [1351.03, 1357.62] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d5/construct | 849.127 ns [846.822, 851.485] (0.95 confidence) | 839.765 ns [836.19, 843.591] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d5/prepared | 425.163 ns [424.352, 426.119] (0.95 confidence) | 415.03 ns [414.224, 415.814] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d6/adapter | 1763.1 ns [1759.51, 1766.84] (0.95 confidence) | 1728.79 ns [1726.52, 1730.96] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d6/construct | 1116.31 ns [1114.01, 1118.89] (0.95 confidence) | 1096.6 ns [1095.19, 1098.02] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d6/prepared | 514.891 ns [514.123, 515.65] (0.95 confidence) | 499.039 ns [498.197, 499.851] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d7/adapter | 2204.8 ns [2199.68, 2210.42] (0.95 confidence) | 2150.69 ns [2147.19, 2154.13] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d7/construct | 1457.9 ns [1454.78, 1461.07] (0.95 confidence) | 1402.75 ns [1400.24, 1405.12] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d7/prepared | 603.939 ns [602.624, 605.58] (0.95 confidence) | 596.728 ns [595.815, 597.643] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d8/adapter | 2521.62 ns [2515.67, 2527.56] (0.95 confidence) | 2462.39 ns [2458.56, 2466.04] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d8/construct | 1743.96 ns [1737.4, 1750.6] (0.95 confidence) | 1680.54 ns [1677.93, 1683.19] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d8/prepared | 656.806 ns [650.007, 666.335] (0.95 confidence) | 635.407 ns [634.624, 636.182] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d2/adapter | 603.613 ns [597.877, 609.451] (0.95 confidence) | 562.604 ns [561.163, 564.135] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d2/construct | 308.988 ns [307.742, 310.187] (0.95 confidence) | 303.18 ns [301.846, 304.443] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d2/prepared | 189.656 ns [189.251, 190.106] (0.95 confidence) | 188.448 ns [188.14, 188.749] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d3/adapter | 923.012 ns [921.624, 924.507] (0.95 confidence) | 900.276 ns [898.877, 901.697] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d3/construct | 537.414 ns [534.717, 539.868] (0.95 confidence) | 521.415 ns [518.846, 523.801] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d3/prepared | 282.271 ns [281.773, 282.87] (0.95 confidence) | 278.5 ns [278.011, 278.984] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d4/adapter | 1292.73 ns [1290.52, 1295.08] (0.95 confidence) | 1248.05 ns [1246.03, 1250.09] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d4/construct | 789.175 ns [786.058, 792.185] (0.95 confidence) | 750.974 ns [746.853, 755.034] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d4/prepared | 367.644 ns [367.048, 368.303] (0.95 confidence) | 367.243 ns [366.622, 367.878] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d5/adapter | 1715.7 ns [1712.66, 1718.9] (0.95 confidence) | 1662.03 ns [1659.67, 1664.49] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d5/construct | 1072.55 ns [1067.5, 1077.4] (0.95 confidence) | 1058.18 ns [1052.72, 1063.23] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d5/prepared | 485.341 ns [484.082, 486.736] (0.95 confidence) | 479.478 ns [478.083, 480.945] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d6/adapter | 2163.97 ns [2159.46, 2168.58] (0.95 confidence) | 2068.75 ns [2065.12, 2072.5] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d6/construct | 1408.7 ns [1401.82, 1415.01] (0.95 confidence) | 1364.29 ns [1354.94, 1373.18] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d6/prepared | 579.295 ns [577.478, 581.244] (0.95 confidence) | 562.872 ns [561.4, 564.376] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d7/adapter | 2658.76 ns [2651.98, 2665.64] (0.95 confidence) | 2550.48 ns [2546.37, 2554.78] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d7/construct | 1802.94 ns [1793.78, 1811.7] (0.95 confidence) | 1751.04 ns [1742.2, 1759.63] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d7/prepared | 683.299 ns [680.909, 685.81] (0.95 confidence) | 667.841 ns [666.171, 669.607] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d8/adapter | 3077.74 ns [3071.26, 3084.45] (0.95 confidence) | 2977.44 ns [2972.63, 2982.12] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d8/construct | 2144.63 ns [2134.97, 2153.84] (0.95 confidence) | 2065.21 ns [2052.89, 2076.52] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d8/prepared | 777.379 ns [774.065, 781.02] (0.95 confidence) | 756.408 ns [754.138, 758.477] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d2/adapter | 548.364 ns [542.914, 553.915] (0.95 confidence) | 517.674 ns [513.685, 522.065] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d2/construct | 244.598 ns [243.968, 245.277] (0.95 confidence) | 240.979 ns [240.028, 241.983] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d2/prepared | 217.686 ns [217.35, 218.057] (0.95 confidence) | 217.389 ns [217.13, 217.648] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d3/adapter | 1028.54 ns [1026.25, 1030.82] (0.95 confidence) | 1025.96 ns [1024.11, 1027.84] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d3/construct | 487.825 ns [486.761, 488.94] (0.95 confidence) | 477.231 ns [476.325, 478.175] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d3/prepared | 457.707 ns [456.604, 459.082] (0.95 confidence) | 450.343 ns [449.723, 450.939] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d4/adapter | 1614.05 ns [1611.87, 1616.39] (0.95 confidence) | 1586.76 ns [1584.38, 1589.05] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d4/construct | 825.919 ns [824.496, 827.347] (0.95 confidence) | 799.759 ns [798.433, 801.12] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d4/prepared | 645.81 ns [644.783, 647.024] (0.95 confidence) | 650.079 ns [649.322, 650.811] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d5/adapter | 2444.22 ns [2440.04, 2448.5] (0.95 confidence) | 2406.23 ns [2400.82, 2411.95] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d5/construct | 1249.68 ns [1245.86, 1253.57] (0.95 confidence) | 1225.18 ns [1221.54, 1229.06] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d5/prepared | 992.484 ns [990.157, 994.759] (0.95 confidence) | 967.638 ns [966.639, 968.606] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d6/adapter | 3316.25 ns [3308.23, 3324.43] (0.95 confidence) | 3237.41 ns [3232.81, 3242.05] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d6/construct | 1736.96 ns [1732.27, 1741.86] (0.95 confidence) | 1707.47 ns [1702.64, 1712.66] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d6/prepared | 1323.63 ns [1321.38, 1325.94] (0.95 confidence) | 1305.78 ns [1303.97, 1307.59] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d7/adapter | 4305.71 ns [4298.16, 4313.74] (0.95 confidence) | 4247.89 ns [4242.65, 4252.96] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d7/construct | 2308.15 ns [2303.39, 2313.22] (0.95 confidence) | 2274.75 ns [2270.01, 2279.68] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d7/prepared | 1743.54 ns [1740.78, 1746.43] (0.95 confidence) | 1726.69 ns [1724.39, 1729.06] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d8/adapter | 5314.42 ns [5304.29, 5324.84] (0.95 confidence) | 5174.94 ns [5163.02, 5186.41] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d8/construct | 2950.83 ns [2943.83, 2958.28] (0.95 confidence) | 2924.35 ns [2914.59, 2934.97] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d8/prepared | 2157.01 ns [2143.77, 2175.19] (0.95 confidence) | 2102.79 ns [2099.2, 2106.29] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d2/adapter | 678.793 ns [673.524, 684.196] (0.95 confidence) | 633.362 ns [632.523, 634.187] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d2/construct | 348.319 ns [346.709, 349.908] (0.95 confidence) | 341.88 ns [340.378, 343.367] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d2/prepared | 216.269 ns [215.959, 216.594] (0.95 confidence) | 212.652 ns [212.356, 212.953] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d3/adapter | 1164.45 ns [1162.79, 1166.14] (0.95 confidence) | 1136.76 ns [1135.26, 1138.26] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d3/construct | 661.341 ns [657.921, 664.583] (0.95 confidence) | 648.006 ns [644.112, 651.632] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d3/prepared | 372.34 ns [371.829, 372.882] (0.95 confidence) | 363.856 ns [363.286, 364.413] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d4/adapter | 1837.74 ns [1835.69, 1839.84] (0.95 confidence) | 1781.94 ns [1779.79, 1784.08] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d4/construct | 1056.12 ns [1050.22, 1061.79] (0.95 confidence) | 1039.84 ns [1035.66, 1043.8] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d4/prepared | 590.131 ns [589.469, 590.847] (0.95 confidence) | 578.821 ns [577.729, 579.844] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d5/adapter | 2657.03 ns [2650.98, 2663.31] (0.95 confidence) | 2557.89 ns [2553.7, 2562.12] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d5/construct | 1595.94 ns [1586.87, 1604.22] (0.95 confidence) | 1559.91 ns [1553.07, 1566.41] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d5/prepared | 823.61 ns [822.017, 825.305] (0.95 confidence) | 803.917 ns [802.195, 805.627] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d6/adapter | 3608.52 ns [3601.91, 3615.55] (0.95 confidence) | 3493.39 ns [3487.91, 3498.49] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d6/construct | 2320.82 ns [2312.33, 2328.88] (0.95 confidence) | 2232.98 ns [2218.76, 2246.11] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d6/prepared | 1043.03 ns [1040.62, 1045.5] (0.95 confidence) | 1018.99 ns [1016.44, 1021.61] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d7/adapter | 4767.19 ns [4758.72, 4775.81] (0.95 confidence) | 4563.24 ns [4557.82, 4568.51] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d7/construct | 3061.54 ns [3048.41, 3074.18] (0.95 confidence) | 2968.97 ns [2957.86, 2979.61] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d7/prepared | 1396.11 ns [1393.1, 1399.25] (0.95 confidence) | 1355.18 ns [1352.77, 1357.61] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d8/adapter | 5922.14 ns [5909.38, 5935.41] (0.95 confidence) | 5702.17 ns [5694.79, 5709.49] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d8/construct | 3938.44 ns [3922.28, 3954.49] (0.95 confidence) | 3812.37 ns [3794.8, 3828.93] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d8/prepared | 1650.86 ns [1645.96, 1656.12] (0.95 confidence) | 1614.26 ns [1610.76, 1617.81] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d2/adapter | 548.998 ns [543.764, 554.336] (0.95 confidence) | 517.035 ns [512.973, 521.53] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d2/construct | 245.714 ns [245.14, 246.326] (0.95 confidence) | 239.25 ns [238.897, 239.634] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d2/prepared | 217.537 ns [217.257, 217.815] (0.95 confidence) | 217.059 ns [216.743, 217.365] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d3/adapter | 1040.72 ns [1038.23, 1043.24] (0.95 confidence) | 1028.09 ns [1026.23, 1029.9] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d3/construct | 493.856 ns [492.317, 495.417] (0.95 confidence) | 482.98 ns [481.68, 484.311] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d3/prepared | 457.754 ns [457.12, 458.379] (0.95 confidence) | 447.262 ns [446.549, 447.985] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d4/adapter | 1618.49 ns [1616.2, 1620.9] (0.95 confidence) | 1593.88 ns [1591.63, 1596.15] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d4/construct | 833.201 ns [831.128, 835.497] (0.95 confidence) | 804.023 ns [802.036, 806.02] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d4/prepared | 647.168 ns [646.233, 648.157] (0.95 confidence) | 651.143 ns [650.46, 651.808] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d5/adapter | 2450.28 ns [2444.3, 2456.38] (0.95 confidence) | 2396.03 ns [2391.58, 2400.54] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d5/construct | 1256.85 ns [1251.97, 1262.14] (0.95 confidence) | 1220.48 ns [1216.35, 1224.98] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d5/prepared | 995.197 ns [992.982, 997.369] (0.95 confidence) | 963.025 ns [961.053, 964.878] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d6/adapter | 3300.12 ns [3293.79, 3306.74] (0.95 confidence) | 3260.58 ns [3254.9, 3266.37] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d6/construct | 1739.56 ns [1736.3, 1742.86] (0.95 confidence) | 1716.03 ns [1710.37, 1721.87] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d6/prepared | 1326.69 ns [1324.18, 1329.29] (0.95 confidence) | 1309.28 ns [1306.81, 1311.91] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d7/adapter | 4325.35 ns [4319.08, 4331.82] (0.95 confidence) | 4259.3 ns [4253.72, 4265.15] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d7/construct | 2305.02 ns [2299.65, 2310.6] (0.95 confidence) | 2248.37 ns [2241.43, 2255.59] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d7/prepared | 1747.75 ns [1744.84, 1750.75] (0.95 confidence) | 1725.2 ns [1722.38, 1728.06] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d8/adapter | 5322.37 ns [5313, 5332.04] (0.95 confidence) | 5202.13 ns [5194.2, 5209.95] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d8/construct | 2942.23 ns [2936.04, 2948.46] (0.95 confidence) | 2881.28 ns [2875.54, 2887.02] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d8/prepared | 2136.2 ns [2132.99, 2139.47] (0.95 confidence) | 2124.49 ns [2121.74, 2127.16] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d2/adapter | 678.582 ns [673.047, 684.157] (0.95 confidence) | 632.531 ns [631.763, 633.3] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d2/construct | 349.929 ns [348.386, 351.42] (0.95 confidence) | 342.849 ns [341.499, 344.147] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d2/prepared | 216.254 ns [215.884, 216.636] (0.95 confidence) | 212.725 ns [212.43, 213.019] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d3/adapter | 1161.83 ns [1160.32, 1163.4] (0.95 confidence) | 1140.29 ns [1138.14, 1142.48] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d3/construct | 667.208 ns [664.286, 670.159] (0.95 confidence) | 652.026 ns [649.464, 654.517] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d3/prepared | 374.069 ns [373.398, 374.814] (0.95 confidence) | 361.451 ns [360.791, 362.1] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d4/adapter | 1829 ns [1827.09, 1830.95] (0.95 confidence) | 1779.47 ns [1776.47, 1782.62] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d4/construct | 1058.48 ns [1052.7, 1063.96] (0.95 confidence) | 1024.53 ns [1018.5, 1030.24] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d4/prepared | 588.327 ns [587.563, 589.145] (0.95 confidence) | 581.562 ns [580.397, 582.699] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d5/adapter | 2642.61 ns [2637.64, 2647.64] (0.95 confidence) | 2571.08 ns [2565.76, 2576.25] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d5/construct | 1602.94 ns [1592.57, 1612.38] (0.95 confidence) | 1562.31 ns [1554.41, 1569.88] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d5/prepared | 824.228 ns [822.55, 825.976] (0.95 confidence) | 802.342 ns [800.935, 803.74] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d6/adapter | 3617.04 ns [3610.71, 3622.93] (0.95 confidence) | 3489.45 ns [3484.13, 3494.6] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d6/construct | 2318.22 ns [2308.9, 2327.89] (0.95 confidence) | 2226.18 ns [2213.84, 2237.39] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d6/prepared | 1041.39 ns [1038.89, 1044.04] (0.95 confidence) | 1014.8 ns [1012.25, 1017.44] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d7/adapter | 4768.94 ns [4759.61, 4778.49] (0.95 confidence) | 4591.83 ns [4583.39, 4601] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d7/construct | 3037.26 ns [3017.16, 3055.87] (0.95 confidence) | 2991.25 ns [2981.97, 3000.42] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d7/prepared | 1392.6 ns [1389.58, 1395.77] (0.95 confidence) | 1357.58 ns [1355.27, 1359.97] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d8/adapter | 5949.06 ns [5937.23, 5962.49] (0.95 confidence) | 5745.66 ns [5738.89, 5752.62] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d8/construct | 3931.88 ns [3909.63, 3952.46] (0.95 confidence) | 3809.59 ns [3789.19, 3828.73] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d8/prepared | 1666.48 ns [1662.49, 1670.53] (0.95 confidence) | 1629.14 ns [1625.25, 1633.1] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d2/adapter | 554.825 ns [549.466, 560.229] (0.95 confidence) | 522.125 ns [518.234, 526.408] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d2/construct | 239.658 ns [239.278, 240.041] (0.95 confidence) | 233.528 ns [233.165, 233.899] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d2/prepared | 217.648 ns [217.362, 217.967] (0.95 confidence) | 217.564 ns [217.243, 217.872] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d3/adapter | 1018.43 ns [1017.04, 1019.85] (0.95 confidence) | 996.283 ns [993.518, 998.57] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d3/construct | 480.202 ns [478.791, 481.761] (0.95 confidence) | 469.291 ns [468.508, 470.087] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d3/prepared | 445.166 ns [444.349, 446.035] (0.95 confidence) | 433.937 ns [433.417, 434.46] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d4/adapter | 1650.77 ns [1647.43, 1654.48] (0.95 confidence) | 1610.77 ns [1607.68, 1613.8] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d4/construct | 823.859 ns [821.642, 826.123] (0.95 confidence) | 797.739 ns [796.317, 799.155] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d4/prepared | 674.603 ns [673.604, 675.65] (0.95 confidence) | 680.569 ns [679.283, 681.755] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d5/adapter | 2382.02 ns [2376.97, 2387.19] (0.95 confidence) | 2327.81 ns [2322.22, 2333.44] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d5/construct | 1212.27 ns [1207.91, 1216.94] (0.95 confidence) | 1184.95 ns [1181.09, 1189.26] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d5/prepared | 972.604 ns [970.603, 974.78] (0.95 confidence) | 942.946 ns [941.876, 944.015] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d6/adapter | 3374.83 ns [3368.91, 3380.9] (0.95 confidence) | 3316.72 ns [3311.2, 3322.27] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d6/construct | 1798.43 ns [1794.13, 1802.76] (0.95 confidence) | 1766.62 ns [1761.38, 1771.84] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d6/prepared | 1355.75 ns [1353.2, 1358.92] (0.95 confidence) | 1340.27 ns [1338.92, 1341.57] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d7/adapter | 4407.49 ns [4400.37, 4414.76] (0.95 confidence) | 4295.26 ns [4287.82, 4302.74] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d7/construct | 2371.72 ns [2364.73, 2378.89] (0.95 confidence) | 2318.76 ns [2310.93, 2326.83] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d7/prepared | 1768.33 ns [1765.09, 1771.7] (0.95 confidence) | 1749.19 ns [1746.59, 1751.71] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d8/adapter | 5494.29 ns [5481.28, 5507.28] (0.95 confidence) | 5382.23 ns [5370.2, 5394.51] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d8/construct | 3136.44 ns [3127.05, 3145.86] (0.95 confidence) | 3038.21 ns [3029.73, 3046.55] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d8/prepared | 2107.38 ns [2104.98, 2109.87] (0.95 confidence) | 2087.75 ns [2084.49, 2091.05] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d2/adapter | 713.345 ns [708.189, 718.621] (0.95 confidence) | 666.809 ns [665.833, 667.752] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d2/construct | 374.894 ns [372.857, 376.779] (0.95 confidence) | 366.013 ns [364.076, 367.802] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d2/prepared | 231.1 ns [230.064, 232.234] (0.95 confidence) | 225.296 ns [224.982, 225.623] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d3/adapter | 1190.27 ns [1188.59, 1191.98] (0.95 confidence) | 1162.93 ns [1161.58, 1164.24] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d3/construct | 689.822 ns [686.926, 692.467] (0.95 confidence) | 674.799 ns [671.302, 678.047] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d3/prepared | 389.486 ns [388.543, 390.536] (0.95 confidence) | 381.942 ns [381.425, 382.443] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d4/adapter | 1788.5 ns [1783.91, 1793.11] (0.95 confidence) | 1738.61 ns [1735.71, 1741.65] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d4/construct | 1007.77 ns [1001.02, 1014.07] (0.95 confidence) | 970.597 ns [964.894, 975.801] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d4/prepared | 620.96 ns [620.135, 621.841] (0.95 confidence) | 613.586 ns [612.81, 614.373] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d5/adapter | 2544.79 ns [2537.72, 2552.07] (0.95 confidence) | 2459.31 ns [2455.11, 2463.72] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d5/construct | 1449.37 ns [1443.92, 1454.63] (0.95 confidence) | 1413.68 ns [1405.32, 1421.28] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d5/prepared | 887.931 ns [886.455, 889.478] (0.95 confidence) | 866.572 ns [865.476, 867.666] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d6/adapter | 3446.05 ns [3437.55, 3454.93] (0.95 confidence) | 3318.37 ns [3312.17, 3324.67] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d6/construct | 2031.92 ns [2023.26, 2040.27] (0.95 confidence) | 1965.3 ns [1954.85, 1974.93] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d6/prepared | 1149.45 ns [1147.67, 1151.28] (0.95 confidence) | 1118.86 ns [1116.84, 1121.02] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d7/adapter | 4568.77 ns [4560.81, 4576.94] (0.95 confidence) | 4366.59 ns [4358.73, 4374.79] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d7/construct | 2782.5 ns [2772, 2792.47] (0.95 confidence) | 2685.35 ns [2671.94, 2697.92] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d7/prepared | 1482.09 ns [1479.25, 1485.08] (0.95 confidence) | 1452.17 ns [1450.62, 1453.71] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d8/adapter | 5826.47 ns [5811.2, 5841.78] (0.95 confidence) | 5759.81 ns [5736.13, 5784.26] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d8/construct | 3605.24 ns [3588.26, 3622.35] (0.95 confidence) | 3522.74 ns [3509.06, 3536.56] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d8/prepared | 1892.2 ns [1889.56, 1894.87] (0.95 confidence) | 1853.73 ns [1849.77, 1858.02] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d2/adapter | 1041.23 ns [1032.37, 1050.46] (0.95 confidence) | 1075.11 ns [1057.69, 1093.29] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d2/construct | 293.648 ns [293.112, 294.183] (0.95 confidence) | 287.203 ns [286.595, 287.815] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d2/prepared | 618.178 ns [617.18, 619.248] (0.95 confidence) | 629.269 ns [628.312, 630.29] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d3/adapter | 1806.22 ns [1804.08, 1808.48] (0.95 confidence) | 1754.07 ns [1750.91, 1757.27] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d3/construct | 579.765 ns [578.517, 581.04] (0.95 confidence) | 572.241 ns [570.72, 573.788] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d3/prepared | 1102.53 ns [1096.81, 1109.83] (0.95 confidence) | 1069.87 ns [1068.23, 1071.58] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d4/adapter | 2717.83 ns [2713.85, 2722.03] (0.95 confidence) | 2674.49 ns [2670.95, 2677.96] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d4/construct | 810.246 ns [807.65, 812.889] (0.95 confidence) | 782.146 ns [779.945, 784.354] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d4/prepared | 1742.66 ns [1739.68, 1745.63] (0.95 confidence) | 1737.94 ns [1733.81, 1742.16] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d5/adapter | 4213.51 ns [4207.98, 4219.18] (0.95 confidence) | 4151.84 ns [4146.79, 4158.22] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d5/construct | 1326.62 ns [1322.86, 1330.5] (0.95 confidence) | 1313.83 ns [1309.9, 1317.84] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d5/prepared | 2688.12 ns [2676.45, 2700.92] (0.95 confidence) | 2606.4 ns [2602.63, 2610.18] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d6/adapter | 5662.16 ns [5645.53, 5681.96] (0.95 confidence) | 5577.91 ns [5560.18, 5597.99] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d6/construct | 1704.57 ns [1700.57, 1708.51] (0.95 confidence) | 1655.54 ns [1651.11, 1660.07] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d6/prepared | 3675.27 ns [3669.77, 3680.87] (0.95 confidence) | 3607.11 ns [3603.01, 3611.01] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d7/adapter | 7939.63 ns [7925.33, 7954.73] (0.95 confidence) | 7705.28 ns [7694.61, 7716.07] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d7/construct | 2231.79 ns [2223.56, 2240.42] (0.95 confidence) | 2147.55 ns [2141.27, 2154.05] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d7/prepared | 5334.93 ns [5324.62, 5345.79] (0.95 confidence) | 5187.83 ns [5182.09, 5193.51] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d8/adapter | 9603.44 ns [9583.21, 9626.49] (0.95 confidence) | 9403.26 ns [9391.73, 9414.56] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d8/construct | 2428.43 ns [2422, 2434.86] (0.95 confidence) | 2353.22 ns [2347.75, 2358.7] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d8/prepared | 7049.64 ns [7036.28, 7065.16] (0.95 confidence) | 6856.97 ns [6847.29, 6866.12] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d2/adapter | 1176.62 ns [1170.58, 1183.15] (0.95 confidence) | 1218.14 ns [1202.45, 1234.31] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d2/construct | 479.144 ns [477.467, 480.832] (0.95 confidence) | 468.185 ns [465.329, 470.683] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d2/prepared | 571.673 ns [570.855, 572.553] (0.95 confidence) | 579.856 ns [578.959, 580.764] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d3/adapter | 1933.28 ns [1930.34, 1936.61] (0.95 confidence) | 1884.78 ns [1881.8, 1887.83] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d3/construct | 750.776 ns [748.655, 752.826] (0.95 confidence) | 727.395 ns [725.414, 729.316] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d3/prepared | 1027.22 ns [1025.8, 1028.68] (0.95 confidence) | 1019.2 ns [1017.94, 1020.47] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d4/adapter | 3164 ns [3159.67, 3168.87] (0.95 confidence) | 3061.35 ns [3054.01, 3067.43] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d4/construct | 1271.86 ns [1265.82, 1277.23] (0.95 confidence) | 1241.53 ns [1237.96, 1244.83] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d4/prepared | 1661.63 ns [1658.94, 1664.41] (0.95 confidence) | 1639.83 ns [1637.57, 1642.04] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d5/adapter | 4425.53 ns [4419.69, 4431.25] (0.95 confidence) | 4282.79 ns [4276.37, 4288.88] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d5/construct | 1671.28 ns [1665.31, 1676.68] (0.95 confidence) | 1624.53 ns [1617.69, 1630.54] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d5/prepared | 2486.95 ns [2482.47, 2491.74] (0.95 confidence) | 2431.78 ns [2429.04, 2434.52] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d6/adapter | 6380 ns [6182.24, 6754.85] (0.95 confidence) | 5836.16 ns [5824.97, 5848.02] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d6/construct | 2272.2 ns [2263.06, 2280.69] (0.95 confidence) | 2219.58 ns [2211.92, 2226.76] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d6/prepared | 3485.14 ns [3479.72, 3490.91] (0.95 confidence) | 3395.43 ns [3391.44, 3399.42] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d7/adapter | 8585.96 ns [8565.56, 8608.96] (0.95 confidence) | 8434.43 ns [8422.88, 8445.52] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d7/construct | 2925.44 ns [2915.18, 2935.93] (0.95 confidence) | 2859.35 ns [2844.75, 2874.73] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d7/prepared | 5403.54 ns [5393.55, 5413.7] (0.95 confidence) | 5366.49 ns [5358.64, 5374.28] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d8/adapter | 10597.7 ns [10583, 10612.6] (0.95 confidence) | 10322.7 ns [10308.7, 10336.3] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d8/construct | 3436.79 ns [3428.25, 3445.67] (0.95 confidence) | 3326.76 ns [3318.57, 3334.26] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d8/prepared | 6826.19 ns [6814.63, 6838.22] (0.95 confidence) | 6717.26 ns [6707.82, 6726.6] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d2/adapter | 555.172 ns [549.806, 560.623] (0.95 confidence) | 521.417 ns [517.369, 525.844] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d2/construct | 241.43 ns [240.947, 241.947] (0.95 confidence) | 236.185 ns [235.764, 236.605] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d2/prepared | 218.468 ns [218.202, 218.752] (0.95 confidence) | 217.812 ns [217.507, 218.109] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d3/adapter | 1012.78 ns [1011.49, 1014.11] (0.95 confidence) | 999.341 ns [997.213, 1001.5] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d3/construct | 477.312 ns [476.19, 478.491] (0.95 confidence) | 470.523 ns [469.822, 471.209] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d3/prepared | 443.964 ns [443.505, 444.439] (0.95 confidence) | 436.436 ns [435.917, 436.977] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d4/adapter | 1659.44 ns [1656.64, 1662.3] (0.95 confidence) | 1621.13 ns [1618.47, 1623.77] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d4/construct | 823.452 ns [821.601, 825.311] (0.95 confidence) | 802.89 ns [801.578, 804.171] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d4/prepared | 674.55 ns [673.589, 675.433] (0.95 confidence) | 679.577 ns [678.472, 680.607] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d5/adapter | 2358.64 ns [2354.61, 2363.01] (0.95 confidence) | 2313.74 ns [2310.56, 2317.05] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d5/construct | 1217.71 ns [1212.8, 1222.82] (0.95 confidence) | 1193.61 ns [1188.82, 1198.86] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d5/prepared | 963.863 ns [962.343, 965.508] (0.95 confidence) | 941.643 ns [940.15, 942.949] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d6/adapter | 3382.77 ns [3374.02, 3391.62] (0.95 confidence) | 3339.19 ns [3333.77, 3344.44] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d6/construct | 1810.69 ns [1806.09, 1815.47] (0.95 confidence) | 1810.81 ns [1792.61, 1835.2] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d6/prepared | 1360.27 ns [1357.76, 1362.97] (0.95 confidence) | 1350.55 ns [1346.93, 1354.4] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d7/adapter | 4383.51 ns [4375.52, 4392.19] (0.95 confidence) | 4300.35 ns [4293.45, 4307.39] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d7/construct | 2385.57 ns [2378.66, 2392.58] (0.95 confidence) | 2317.73 ns [2311.07, 2324.71] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d7/prepared | 1776.49 ns [1774.27, 1778.85] (0.95 confidence) | 1754.98 ns [1752.34, 1757.49] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d8/adapter | 5464.54 ns [5451.82, 5477.77] (0.95 confidence) | 5415.02 ns [5402.91, 5427.36] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d8/construct | 3146.69 ns [3136.62, 3156.5] (0.95 confidence) | 3037.91 ns [3029.75, 3046.07] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d8/prepared | 2110.21 ns [2106.72, 2113.87] (0.95 confidence) | 2092.23 ns [2089.28, 2095.16] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d2/adapter | 720.139 ns [715.525, 724.918] (0.95 confidence) | 669.92 ns [668.637, 671.112] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d2/construct | 373.132 ns [371.112, 375.049] (0.95 confidence) | 370.255 ns [368.034, 372.299] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d2/prepared | 230.022 ns [229.658, 230.409] (0.95 confidence) | 225.9 ns [225.59, 226.214] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d3/adapter | 1190.65 ns [1188.8, 1192.67] (0.95 confidence) | 1167.27 ns [1165.53, 1169.13] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d3/construct | 689.602 ns [685.796, 693.202] (0.95 confidence) | 674.71 ns [670.865, 678.287] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d3/prepared | 390.04 ns [389.452, 390.637] (0.95 confidence) | 382.178 ns [381.301, 383.153] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d4/adapter | 1798.76 ns [1795.32, 1802.27] (0.95 confidence) | 1738.73 ns [1736.44, 1740.98] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d4/construct | 1006.68 ns [1000.88, 1012.1] (0.95 confidence) | 975.326 ns [969.739, 980.487] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d4/prepared | 623.921 ns [622.649, 625.423] (0.95 confidence) | 624.822 ns [622.686, 626.973] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d5/adapter | 2522.67 ns [2518.06, 2527.62] (0.95 confidence) | 2468.05 ns [2463.44, 2472.74] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d5/construct | 1442.94 ns [1436.28, 1449.36] (0.95 confidence) | 1404.2 ns [1393.58, 1414.33] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d5/prepared | 882.186 ns [881.038, 883.361] (0.95 confidence) | 867.546 ns [866.491, 868.587] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d6/adapter | 3430.49 ns [3424.91, 3436.37] (0.95 confidence) | 3301.95 ns [3296.03, 3307.88] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d6/construct | 2046.37 ns [2039.12, 2053.87] (0.95 confidence) | 1978.78 ns [1969.51, 1987.18] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d6/prepared | 1146.57 ns [1144.07, 1149.51] (0.95 confidence) | 1112.32 ns [1111.01, 1113.59] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d7/adapter | 4559.43 ns [4551.3, 4567.61] (0.95 confidence) | 4382.92 ns [4374.68, 4391.29] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d7/construct | 2788.7 ns [2780.93, 2796.15] (0.95 confidence) | 2657.91 ns [2646.29, 2668.65] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d7/prepared | 1482.39 ns [1479.05, 1485.88] (0.95 confidence) | 1454.68 ns [1452.96, 1456.36] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d8/adapter | 5831.76 ns [5819.76, 5843.92] (0.95 confidence) | 5620.58 ns [5611.16, 5630.38] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d8/construct | 3671.77 ns [3644.82, 3701.74] (0.95 confidence) | 3513.48 ns [3499.15, 3527.56] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d8/prepared | 1892.4 ns [1887.75, 1897.45] (0.95 confidence) | 1840.37 ns [1837.21, 1843.65] (0.95 confidence) | — | — |

## Median (ns)

| Case | la-stack baseline | la-stack current | nalgebra | faer |
| --- | --- | --- | --- | --- |
| angle&#95;d3/dense/borrowed&#95;slice | 24.0879 ns [24.053, 24.1391] (0.95 confidence) | 23.9313 ns [23.8771, 23.9811] (0.95 confidence) | — | — |
| angle&#95;d3/dense/construct&#95;then&#95;angle | 23.0042 ns [22.9629, 23.0563] (0.95 confidence) | 22.9978 ns [22.9758, 23.0552] (0.95 confidence) | — | — |
| angle&#95;d3/dense/prepared&#95;vector | 22.9881 ns [22.9486, 23.0432] (0.95 confidence) | 22.7778 ns [22.7348, 22.8036] (0.95 confidence) | — | — |
| angle&#95;d3/dense&#95;antipodal/borrowed&#95;slice | 17.0061 ns [16.9812, 17.0607] (0.95 confidence) | 16.8446 ns [16.8037, 16.8776] (0.95 confidence) | — | — |
| angle&#95;d3/dense&#95;antipodal/construct&#95;then&#95;angle | 15.7858 ns [15.7291, 15.8203] (0.95 confidence) | 15.7541 ns [15.7306, 15.7791] (0.95 confidence) | — | — |
| angle&#95;d3/dense&#95;antipodal/prepared&#95;vector | 16.1846 ns [16.113, 16.246] (0.95 confidence) | 15.8541 ns [15.8071, 15.8841] (0.95 confidence) | — | — |
| angle&#95;d3/dense&#95;identical/borrowed&#95;slice | 14.2942 ns [14.2751, 14.3067] (0.95 confidence) | 14.3557 ns [14.3439, 14.3659] (0.95 confidence) | — | — |
| angle&#95;d3/dense&#95;identical/construct&#95;then&#95;angle | 13.5795 ns [13.5569, 13.6171] (0.95 confidence) | 13.5413 ns [13.5179, 13.5588] (0.95 confidence) | — | — |
| angle&#95;d3/dense&#95;identical/prepared&#95;vector | 13.6323 ns [13.611, 13.6611] (0.95 confidence) | 13.5009 ns [13.4802, 13.5241] (0.95 confidence) | — | — |
| angle&#95;d3/dense&#95;varied/borrowed&#95;slice | 26.4423 ns [26.3942, 26.4704] (0.95 confidence) | 26.1458 ns [26.109, 26.1882] (0.95 confidence) | — | — |
| angle&#95;d3/dense&#95;varied/construct&#95;then&#95;angle | 25.5293 ns [25.4808, 25.5818] (0.95 confidence) | 25.1511 ns [25.0971, 25.2075] (0.95 confidence) | — | — |
| angle&#95;d3/dense&#95;varied/prepared&#95;vector | 25.0846 ns [25.057, 25.1147] (0.95 confidence) | 24.6637 ns [24.6111, 24.7016] (0.95 confidence) | — | — |
| angle&#95;d3/maximum&#95;small&#95;angle/borrowed&#95;slice | 10.815 ns [10.7932, 10.8542] (0.95 confidence) | 10.8295 ns [10.8137, 10.8635] (0.95 confidence) | — | — |
| angle&#95;d3/maximum&#95;small&#95;angle/construct&#95;then&#95;angle | 10.4992 ns [10.45, 10.5457] (0.95 confidence) | 10.0515 ns [10.0461, 10.0866] (0.95 confidence) | — | — |
| angle&#95;d3/maximum&#95;small&#95;angle/prepared&#95;vector | 9.58074 ns [9.56357, 9.59163] (0.95 confidence) | 9.50805 ns [9.49195, 9.52284] (0.95 confidence) | — | — |
| angle&#95;d3/mixed&#95;scale/borrowed&#95;slice | 12.2727 ns [12.2609, 12.285] (0.95 confidence) | 12.2582 ns [12.2209, 12.3027] (0.95 confidence) | — | — |
| angle&#95;d3/mixed&#95;scale/construct&#95;then&#95;angle | 11.8298 ns [11.8077, 11.8629] (0.95 confidence) | 11.3027 ns [11.2718, 11.3389] (0.95 confidence) | — | — |
| angle&#95;d3/mixed&#95;scale/prepared&#95;vector | 10.7225 ns [10.7044, 10.737] (0.95 confidence) | 10.6562 ns [10.633, 10.6728] (0.95 confidence) | — | — |
| angle&#95;d3/near&#95;antipodal/borrowed&#95;slice | 16.5843 ns [16.5725, 16.5962] (0.95 confidence) | 16.4534 ns [16.4265, 16.4784] (0.95 confidence) | — | — |
| angle&#95;d3/near&#95;antipodal/construct&#95;then&#95;angle | 15.438 ns [15.4043, 15.4659] (0.95 confidence) | 15.1937 ns [15.1708, 15.2104] (0.95 confidence) | — | — |
| angle&#95;d3/near&#95;antipodal/prepared&#95;vector | 15.1397 ns [15.129, 15.15] (0.95 confidence) | 14.9241 ns [14.905, 14.9382] (0.95 confidence) | — | — |
| angle&#95;d3/near&#95;parallel/borrowed&#95;slice | 10.811 ns [10.7989, 10.8366] (0.95 confidence) | 10.7896 ns [10.767, 10.809] (0.95 confidence) | — | — |
| angle&#95;d3/near&#95;parallel/construct&#95;then&#95;angle | 10.4483 ns [10.3909, 10.4697] (0.95 confidence) | 10.1072 ns [10.0655, 10.1795] (0.95 confidence) | — | — |
| angle&#95;d3/near&#95;parallel/prepared&#95;vector | 9.67713 ns [9.64223, 9.69931] (0.95 confidence) | 9.46857 ns [9.42515, 9.48964] (0.95 confidence) | — | — |
| angle&#95;d3/sparse&#95;adjacent/borrowed&#95;slice | 13.8229 ns [13.8061, 13.8465] (0.95 confidence) | 13.5705 ns [13.5417, 13.6494] (0.95 confidence) | — | — |
| angle&#95;d3/sparse&#95;adjacent/construct&#95;then&#95;angle | 13.3547 ns [13.3277, 13.437] (0.95 confidence) | 13.0015 ns [12.9827, 13.0138] (0.95 confidence) | — | — |
| angle&#95;d3/sparse&#95;adjacent/prepared&#95;vector | 12.0284 ns [12.0044, 12.0524] (0.95 confidence) | 11.7487 ns [11.7256, 11.7635] (0.95 confidence) | — | — |
| angle&#95;d3/sparse&#95;alternating/borrowed&#95;slice | 14.0866 ns [14.0689, 14.107] (0.95 confidence) | 13.9024 ns [13.8207, 13.9395] (0.95 confidence) | — | — |
| angle&#95;d3/sparse&#95;alternating/construct&#95;then&#95;angle | 13.7328 ns [13.6974, 13.7698] (0.95 confidence) | 13.0134 ns [13.0042, 13.028] (0.95 confidence) | — | — |
| angle&#95;d3/sparse&#95;alternating/prepared&#95;vector | 12.3594 ns [12.3287, 12.407] (0.95 confidence) | 12.0133 ns [11.9959, 12.0278] (0.95 confidence) | — | — |
| angle&#95;d3/sparse&#95;antipodal/borrowed&#95;slice | 11.2083 ns [11.1918, 11.2203] (0.95 confidence) | 10.7586 ns [10.7443, 10.7697] (0.95 confidence) | — | — |
| angle&#95;d3/sparse&#95;antipodal/construct&#95;then&#95;angle | 10.4673 ns [10.429, 10.5173] (0.95 confidence) | 9.96194 ns [9.94187, 9.99628] (0.95 confidence) | — | — |
| angle&#95;d3/sparse&#95;antipodal/prepared&#95;vector | 9.0543 ns [9.04153, 9.06567] (0.95 confidence) | 8.7555 ns [8.74616, 8.76475] (0.95 confidence) | — | — |
| angle&#95;d3/sparse&#95;identical/borrowed&#95;slice | 8.87129 ns [8.82802, 8.91111] (0.95 confidence) | 9.17295 ns [9.03425, 9.36623] (0.95 confidence) | — | — |
| angle&#95;d3/sparse&#95;identical/construct&#95;then&#95;angle | 8.31464 ns [8.30173, 8.33224] (0.95 confidence) | 8.40837 ns [8.35572, 8.45124] (0.95 confidence) | — | — |
| angle&#95;d3/sparse&#95;identical/prepared&#95;vector | 7.05948 ns [7.05296, 7.07245] (0.95 confidence) | 7.0937 ns [7.03918, 7.17723] (0.95 confidence) | — | — |
| angle&#95;d3/subnormal/borrowed&#95;slice | 10.8001 ns [10.7908, 10.8074] (0.95 confidence) | 10.8662 ns [10.8227, 10.9329] (0.95 confidence) | — | — |
| angle&#95;d3/subnormal/construct&#95;then&#95;angle | 10.2929 ns [10.2813, 10.3237] (0.95 confidence) | 10.1498 ns [10.0773, 10.2098] (0.95 confidence) | — | — |
| angle&#95;d3/subnormal/prepared&#95;vector | 9.60335 ns [9.57702, 9.62634] (0.95 confidence) | 9.4553 ns [9.43672, 9.50595] (0.95 confidence) | — | — |
| angle&#95;d3/subnormal&#95;axes/borrowed&#95;slice | 14.0384 ns [14.0293, 14.052] (0.95 confidence) | 13.8873 ns [13.8603, 13.9196] (0.95 confidence) | — | — |
| angle&#95;d3/subnormal&#95;axes/construct&#95;then&#95;angle | 13.5365 ns [13.5236, 13.5611] (0.95 confidence) | 13.0367 ns [13.02, 13.0491] (0.95 confidence) | — | — |
| angle&#95;d3/subnormal&#95;axes/prepared&#95;vector | 12.2277 ns [12.22, 12.2422] (0.95 confidence) | 11.9959 ns [11.9803, 12.0068] (0.95 confidence) | — | — |
| angle&#95;d4/dense/borrowed&#95;slice | 34.656 ns [34.6223, 34.6799] (0.95 confidence) | 34.3115 ns [34.2593, 34.4044] (0.95 confidence) | — | — |
| angle&#95;d4/dense/construct&#95;then&#95;angle | 33.1815 ns [33.1539, 33.2213] (0.95 confidence) | 32.9013 ns [32.8496, 32.9507] (0.95 confidence) | — | — |
| angle&#95;d4/dense/prepared&#95;vector | 32.9852 ns [32.9541, 33.0164] (0.95 confidence) | 32.6172 ns [32.5738, 32.6951] (0.95 confidence) | — | — |
| angle&#95;d4/dense&#95;antipodal/borrowed&#95;slice | 26.794 ns [26.7815, 26.8156] (0.95 confidence) | 26.7761 ns [26.7436, 26.8143] (0.95 confidence) | — | — |
| angle&#95;d4/dense&#95;antipodal/construct&#95;then&#95;angle | 25.1597 ns [25.1232, 25.194] (0.95 confidence) | 25.1651 ns [25.1231, 25.2078] (0.95 confidence) | — | — |
| angle&#95;d4/dense&#95;antipodal/prepared&#95;vector | 25.0312 ns [25.0097, 25.0446] (0.95 confidence) | 24.9378 ns [24.9031, 24.9739] (0.95 confidence) | — | — |
| angle&#95;d4/dense&#95;identical/borrowed&#95;slice | 24.3806 ns [24.3381, 24.4148] (0.95 confidence) | 24.3791 ns [24.3446, 24.4053] (0.95 confidence) | — | — |
| angle&#95;d4/dense&#95;identical/construct&#95;then&#95;angle | 22.7818 ns [22.7643, 22.7986] (0.95 confidence) | 22.7864 ns [22.7594, 22.8094] (0.95 confidence) | — | — |
| angle&#95;d4/dense&#95;identical/prepared&#95;vector | 22.5904 ns [22.5548, 22.6098] (0.95 confidence) | 22.2985 ns [22.2541, 22.3387] (0.95 confidence) | — | — |
| angle&#95;d4/dense&#95;varied/borrowed&#95;slice | 36.5218 ns [36.3208, 36.6146] (0.95 confidence) | 35.635 ns [35.5924, 35.7532] (0.95 confidence) | — | — |
| angle&#95;d4/dense&#95;varied/construct&#95;then&#95;angle | 34.6056 ns [34.5447, 34.6324] (0.95 confidence) | 34.0692 ns [34.0193, 34.1603] (0.95 confidence) | — | — |
| angle&#95;d4/dense&#95;varied/prepared&#95;vector | 34.1474 ns [34.0888, 34.2669] (0.95 confidence) | 33.4474 ns [33.3832, 33.5101] (0.95 confidence) | — | — |
| angle&#95;d4/maximum&#95;small&#95;angle/borrowed&#95;slice | 12.1422 ns [12.0988, 12.181] (0.95 confidence) | 12.5233 ns [12.4705, 12.6103] (0.95 confidence) | — | — |
| angle&#95;d4/maximum&#95;small&#95;angle/construct&#95;then&#95;angle | 11.2778 ns [11.2623, 11.2978] (0.95 confidence) | 11.4427 ns [11.3874, 11.493] (0.95 confidence) | — | — |
| angle&#95;d4/maximum&#95;small&#95;angle/prepared&#95;vector | 9.97579 ns [9.95848, 10.0015] (0.95 confidence) | 9.95326 ns [9.93527, 9.97679] (0.95 confidence) | — | — |
| angle&#95;d4/mixed&#95;scale/borrowed&#95;slice | 13.6584 ns [13.6438, 13.6855] (0.95 confidence) | 13.8642 ns [13.7865, 13.9311] (0.95 confidence) | — | — |
| angle&#95;d4/mixed&#95;scale/construct&#95;then&#95;angle | 12.7393 ns [12.6563, 12.7856] (0.95 confidence) | 12.4882 ns [12.4146, 12.5427] (0.95 confidence) | — | — |
| angle&#95;d4/mixed&#95;scale/prepared&#95;vector | 11.378 ns [11.362, 11.4145] (0.95 confidence) | 11.1002 ns [11.0783, 11.1369] (0.95 confidence) | — | — |
| angle&#95;d4/near&#95;antipodal/borrowed&#95;slice | 17.8299 ns [17.8157, 17.8413] (0.95 confidence) | 17.7537 ns [17.7356, 17.7865] (0.95 confidence) | — | — |
| angle&#95;d4/near&#95;antipodal/construct&#95;then&#95;angle | 16.2281 ns [16.2057, 16.2535] (0.95 confidence) | 16.1921 ns [16.1363, 16.2872] (0.95 confidence) | — | — |
| angle&#95;d4/near&#95;antipodal/prepared&#95;vector | 15.8191 ns [15.8104, 15.8258] (0.95 confidence) | 15.568 ns [15.55, 15.6042] (0.95 confidence) | — | — |
| angle&#95;d4/near&#95;parallel/borrowed&#95;slice | 12.079 ns [12.0575, 12.1251] (0.95 confidence) | 12.6153 ns [12.545, 12.8416] (0.95 confidence) | — | — |
| angle&#95;d4/near&#95;parallel/construct&#95;then&#95;angle | 11.1172 ns [11.038, 11.2122] (0.95 confidence) | 11.204 ns [11.1599, 11.2465] (0.95 confidence) | — | — |
| angle&#95;d4/near&#95;parallel/prepared&#95;vector | 9.97161 ns [9.95737, 9.99387] (0.95 confidence) | 9.95015 ns [9.92939, 9.9607] (0.95 confidence) | — | — |
| angle&#95;d4/sparse&#95;adjacent/borrowed&#95;slice | 15.7304 ns [15.7207, 15.7389] (0.95 confidence) | 15.6538 ns [15.6061, 15.7685] (0.95 confidence) | — | — |
| angle&#95;d4/sparse&#95;adjacent/construct&#95;then&#95;angle | 14.9491 ns [14.8677, 14.9851] (0.95 confidence) | 14.7386 ns [14.6699, 14.7769] (0.95 confidence) | — | — |
| angle&#95;d4/sparse&#95;adjacent/prepared&#95;vector | 13.4909 ns [13.482, 13.5035] (0.95 confidence) | 13.2594 ns [13.2398, 13.2728] (0.95 confidence) | — | — |
| angle&#95;d4/sparse&#95;alternating/borrowed&#95;slice | 16.3661 ns [16.2149, 16.4881] (0.95 confidence) | 15.9453 ns [15.8918, 16.0471] (0.95 confidence) | — | — |
| angle&#95;d4/sparse&#95;alternating/construct&#95;then&#95;angle | 15.1135 ns [15.0606, 15.1968] (0.95 confidence) | 14.9343 ns [14.8893, 15.0002] (0.95 confidence) | — | — |
| angle&#95;d4/sparse&#95;alternating/prepared&#95;vector | 13.7929 ns [13.768, 13.8421] (0.95 confidence) | 13.5574 ns [13.5391, 13.573] (0.95 confidence) | — | — |
| angle&#95;d4/sparse&#95;antipodal/borrowed&#95;slice | 12.0043 ns [11.9972, 12.0096] (0.95 confidence) | 11.8817 ns [11.8016, 12.0162] (0.95 confidence) | — | — |
| angle&#95;d4/sparse&#95;antipodal/construct&#95;then&#95;angle | 11.0543 ns [11.0294, 11.1297] (0.95 confidence) | 10.7754 ns [10.7542, 10.7939] (0.95 confidence) | — | — |
| angle&#95;d4/sparse&#95;antipodal/prepared&#95;vector | 9.73619 ns [9.73044, 9.74006] (0.95 confidence) | 9.48821 ns [9.47339, 9.50502] (0.95 confidence) | — | — |
| angle&#95;d4/sparse&#95;identical/borrowed&#95;slice | 10.2378 ns [10.2324, 10.254] (0.95 confidence) | 10.6816 ns [10.6188, 10.7341] (0.95 confidence) | — | — |
| angle&#95;d4/sparse&#95;identical/construct&#95;then&#95;angle | 9.42229 ns [9.28811, 9.46929] (0.95 confidence) | 9.44571 ns [9.40612, 9.48925] (0.95 confidence) | — | — |
| angle&#95;d4/sparse&#95;identical/prepared&#95;vector | 7.98672 ns [7.97871, 8.02187] (0.95 confidence) | 7.99104 ns [7.98051, 8.00166] (0.95 confidence) | — | — |
| angle&#95;d4/subnormal/borrowed&#95;slice | 12.1235 ns [12.0937, 12.1786] (0.95 confidence) | 12.7157 ns [12.5575, 12.8383] (0.95 confidence) | — | — |
| angle&#95;d4/subnormal/construct&#95;then&#95;angle | 11.1323 ns [11.0625, 11.1975] (0.95 confidence) | 11.5185 ns [11.476, 11.5538] (0.95 confidence) | — | — |
| angle&#95;d4/subnormal/prepared&#95;vector | 9.94724 ns [9.9371, 9.9593] (0.95 confidence) | 10.2349 ns [10.1352, 10.29] (0.95 confidence) | — | — |
| angle&#95;d4/subnormal&#95;axes/borrowed&#95;slice | 15.404 ns [15.359, 15.4599] (0.95 confidence) | 15.4124 ns [15.2684, 15.492] (0.95 confidence) | — | — |
| angle&#95;d4/subnormal&#95;axes/construct&#95;then&#95;angle | 14.4777 ns [14.4551, 14.5172] (0.95 confidence) | 14.1978 ns [14.1651, 14.2689] (0.95 confidence) | — | — |
| angle&#95;d4/subnormal&#95;axes/prepared&#95;vector | 13.2735 ns [13.2411, 13.3481] (0.95 confidence) | 12.8013 ns [12.7826, 12.8136] (0.95 confidence) | — | — |
| angle&#95;d5/dense/borrowed&#95;slice | 48.3883 ns [48.3617, 48.4319] (0.95 confidence) | 48.8853 ns [48.7979, 48.9584] (0.95 confidence) | — | — |
| angle&#95;d5/dense/construct&#95;then&#95;angle | 46.1089 ns [46.0813, 46.1409] (0.95 confidence) | 46.3286 ns [46.2686, 46.3619] (0.95 confidence) | — | — |
| angle&#95;d5/dense/prepared&#95;vector | 45.7418 ns [45.7195, 45.7651] (0.95 confidence) | 45.8803 ns [45.7972, 45.9187] (0.95 confidence) | — | — |
| angle&#95;d5/dense&#95;antipodal/borrowed&#95;slice | 37.2245 ns [37.1895, 37.2494] (0.95 confidence) | 37.3592 ns [37.321, 37.3974] (0.95 confidence) | — | — |
| angle&#95;d5/dense&#95;antipodal/construct&#95;then&#95;angle | 35.1753 ns [35.1512, 35.21] (0.95 confidence) | 35.1687 ns [35.1363, 35.2024] (0.95 confidence) | — | — |
| angle&#95;d5/dense&#95;antipodal/prepared&#95;vector | 35.4377 ns [35.3642, 35.5358] (0.95 confidence) | 35.0346 ns [34.9806, 35.0823] (0.95 confidence) | — | — |
| angle&#95;d5/dense&#95;identical/borrowed&#95;slice | 35.8387 ns [35.8206, 35.8492] (0.95 confidence) | 36.0542 ns [35.9839, 36.0988] (0.95 confidence) | — | — |
| angle&#95;d5/dense&#95;identical/construct&#95;then&#95;angle | 34.2026 ns [34.1722, 34.2743] (0.95 confidence) | 34.2474 ns [34.2053, 34.2999] (0.95 confidence) | — | — |
| angle&#95;d5/dense&#95;identical/prepared&#95;vector | 33.5264 ns [33.4957, 33.5564] (0.95 confidence) | 33.6363 ns [33.5972, 33.6678] (0.95 confidence) | — | — |
| angle&#95;d5/dense&#95;varied/borrowed&#95;slice | 70.6968 ns [70.6172, 70.7631] (0.95 confidence) | 70.7355 ns [70.5995, 70.9089] (0.95 confidence) | — | — |
| angle&#95;d5/dense&#95;varied/construct&#95;then&#95;angle | 68.0814 ns [67.9812, 68.1902] (0.95 confidence) | 67.9341 ns [67.8082, 68.0779] (0.95 confidence) | — | — |
| angle&#95;d5/dense&#95;varied/prepared&#95;vector | 66.0852 ns [65.9581, 66.1965] (0.95 confidence) | 65.9582 ns [65.8551, 66.0359] (0.95 confidence) | — | — |
| angle&#95;d5/maximum&#95;small&#95;angle/borrowed&#95;slice | 13.9205 ns [13.7333, 13.9613] (0.95 confidence) | 14.7046 ns [14.5815, 14.7206] (0.95 confidence) | — | — |
| angle&#95;d5/maximum&#95;small&#95;angle/construct&#95;then&#95;angle | 13.2574 ns [13.2326, 13.3268] (0.95 confidence) | 13.215 ns [13.0873, 13.2677] (0.95 confidence) | — | — |
| angle&#95;d5/maximum&#95;small&#95;angle/prepared&#95;vector | 11.1534 ns [11.0502, 11.2226] (0.95 confidence) | 11.2292 ns [11.1262, 11.2489] (0.95 confidence) | — | — |
| angle&#95;d5/mixed&#95;scale/borrowed&#95;slice | 15.1895 ns [15.1316, 15.3384] (0.95 confidence) | 15.8805 ns [15.7893, 15.9074] (0.95 confidence) | — | — |
| angle&#95;d5/mixed&#95;scale/construct&#95;then&#95;angle | 14.5183 ns [14.4732, 14.6028] (0.95 confidence) | 14.0542 ns [13.9954, 14.1912] (0.95 confidence) | — | — |
| angle&#95;d5/mixed&#95;scale/prepared&#95;vector | 12.5885 ns [12.5372, 12.666] (0.95 confidence) | 12.4583 ns [12.4451, 12.4725] (0.95 confidence) | — | — |
| angle&#95;d5/near&#95;antipodal/borrowed&#95;slice | 19.6198 ns [19.6086, 19.633] (0.95 confidence) | 19.3927 ns [19.3589, 19.453] (0.95 confidence) | — | — |
| angle&#95;d5/near&#95;antipodal/construct&#95;then&#95;angle | 17.4969 ns [17.4794, 17.5141] (0.95 confidence) | 17.0435 ns [17.0349, 17.0555] (0.95 confidence) | — | — |
| angle&#95;d5/near&#95;antipodal/prepared&#95;vector | 16.6754 ns [16.653, 16.6955] (0.95 confidence) | 16.6437 ns [16.6199, 16.6775] (0.95 confidence) | — | — |
| angle&#95;d5/near&#95;parallel/borrowed&#95;slice | 13.6309 ns [13.5722, 13.7096] (0.95 confidence) | 14.7107 ns [14.6967, 14.7202] (0.95 confidence) | — | — |
| angle&#95;d5/near&#95;parallel/construct&#95;then&#95;angle | 13.2195 ns [13.1238, 13.2407] (0.95 confidence) | 12.8078 ns [12.7698, 12.904] (0.95 confidence) | — | — |
| angle&#95;d5/near&#95;parallel/prepared&#95;vector | 11.2494 ns [11.1718, 11.2594] (0.95 confidence) | 11.2527 ns [11.2258, 11.2642] (0.95 confidence) | — | — |
| angle&#95;d5/sparse&#95;adjacent/borrowed&#95;slice | 17.288 ns [17.2748, 17.3208] (0.95 confidence) | 17.5763 ns [17.5363, 17.6537] (0.95 confidence) | — | — |
| angle&#95;d5/sparse&#95;adjacent/construct&#95;then&#95;angle | 16.8192 ns [16.7787, 16.9581] (0.95 confidence) | 16.3453 ns [16.3088, 16.3765] (0.95 confidence) | — | — |
| angle&#95;d5/sparse&#95;adjacent/prepared&#95;vector | 14.8019 ns [14.7872, 14.8266] (0.95 confidence) | 14.5401 ns [14.5217, 14.5713] (0.95 confidence) | — | — |
| angle&#95;d5/sparse&#95;alternating/borrowed&#95;slice | 28.1789 ns [28.1578, 28.208] (0.95 confidence) | 27.7222 ns [27.7068, 27.7357] (0.95 confidence) | — | — |
| angle&#95;d5/sparse&#95;alternating/construct&#95;then&#95;angle | 26.6723 ns [26.6569, 26.692] (0.95 confidence) | 26.6881 ns [26.6776, 26.6984] (0.95 confidence) | — | — |
| angle&#95;d5/sparse&#95;alternating/prepared&#95;vector | 26.3885 ns [26.3639, 26.4072] (0.95 confidence) | 26.191 ns [26.1728, 26.2026] (0.95 confidence) | — | — |
| angle&#95;d5/sparse&#95;antipodal/borrowed&#95;slice | 12.983 ns [12.9743, 12.9985] (0.95 confidence) | 13.2125 ns [13.2077, 13.2199] (0.95 confidence) | — | — |
| angle&#95;d5/sparse&#95;antipodal/construct&#95;then&#95;angle | 12.6743 ns [12.5526, 12.7042] (0.95 confidence) | 11.9117 ns [11.8303, 11.9879] (0.95 confidence) | — | — |
| angle&#95;d5/sparse&#95;antipodal/prepared&#95;vector | 10.5054 ns [10.5001, 10.5135] (0.95 confidence) | 10.2662 ns [10.2572, 10.2725] (0.95 confidence) | — | — |
| angle&#95;d5/sparse&#95;identical/borrowed&#95;slice | 11.4943 ns [11.4809, 11.507] (0.95 confidence) | 11.9907 ns [11.9805, 12.0157] (0.95 confidence) | — | — |
| angle&#95;d5/sparse&#95;identical/construct&#95;then&#95;angle | 10.8441 ns [10.7735, 10.9293] (0.95 confidence) | 10.7178 ns [10.6298, 10.7331] (0.95 confidence) | — | — |
| angle&#95;d5/sparse&#95;identical/prepared&#95;vector | 9.01029 ns [9.00559, 9.01775] (0.95 confidence) | 9.07047 ns [9.04956, 9.09495] (0.95 confidence) | — | — |
| angle&#95;d5/subnormal/borrowed&#95;slice | 13.5669 ns [13.4839, 13.7227] (0.95 confidence) | 14.6281 ns [14.5251, 14.6594] (0.95 confidence) | — | — |
| angle&#95;d5/subnormal/construct&#95;then&#95;angle | 13.0531 ns [12.9617, 13.1621] (0.95 confidence) | 12.7002 ns [12.689, 12.7148] (0.95 confidence) | — | — |
| angle&#95;d5/subnormal/prepared&#95;vector | 11.1689 ns [11.0861, 11.2075] (0.95 confidence) | 11.1348 ns [11, 11.2141] (0.95 confidence) | — | — |
| angle&#95;d5/subnormal&#95;axes/borrowed&#95;slice | 17.5173 ns [17.3397, 17.6332] (0.95 confidence) | 17.6632 ns [17.5584, 17.6969] (0.95 confidence) | — | — |
| angle&#95;d5/subnormal&#95;axes/construct&#95;then&#95;angle | 16.4684 ns [16.4519, 16.4988] (0.95 confidence) | 15.9284 ns [15.8242, 15.9728] (0.95 confidence) | — | — |
| angle&#95;d5/subnormal&#95;axes/prepared&#95;vector | 14.2686 ns [14.2402, 14.3161] (0.95 confidence) | 14.1928 ns [14.0739, 14.2453] (0.95 confidence) | — | — |
| angle&#95;d6/dense/borrowed&#95;slice | 62.6647 ns [62.6234, 62.7485] (0.95 confidence) | 62.8104 ns [62.7539, 62.9215] (0.95 confidence) | — | — |
| angle&#95;d6/dense/construct&#95;then&#95;angle | 60.1026 ns [60.044, 60.1605] (0.95 confidence) | 60.128 ns [60.036, 60.2076] (0.95 confidence) | — | — |
| angle&#95;d6/dense/prepared&#95;vector | 59.8348 ns [59.7577, 59.9062] (0.95 confidence) | 59.7571 ns [59.6729, 59.8381] (0.95 confidence) | — | — |
| angle&#95;d6/dense&#95;antipodal/borrowed&#95;slice | 49.2077 ns [49.1826, 49.2723] (0.95 confidence) | 49.1375 ns [49.1035, 49.1759] (0.95 confidence) | — | — |
| angle&#95;d6/dense&#95;antipodal/construct&#95;then&#95;angle | 46.6901 ns [46.6438, 46.7744] (0.95 confidence) | 46.5981 ns [46.5417, 46.7419] (0.95 confidence) | — | — |
| angle&#95;d6/dense&#95;antipodal/prepared&#95;vector | 46.2149 ns [46.1521, 46.2581] (0.95 confidence) | 46.1783 ns [46.1355, 46.2108] (0.95 confidence) | — | — |
| angle&#95;d6/dense&#95;identical/borrowed&#95;slice | 48.0553 ns [48.0105, 48.0876] (0.95 confidence) | 48.6202 ns [48.5535, 48.6731] (0.95 confidence) | — | — |
| angle&#95;d6/dense&#95;identical/construct&#95;then&#95;angle | 45.5389 ns [45.5031, 45.5804] (0.95 confidence) | 45.9969 ns [45.9567, 46.0883] (0.95 confidence) | — | — |
| angle&#95;d6/dense&#95;identical/prepared&#95;vector | 45.012 ns [44.9789, 45.1047] (0.95 confidence) | 45.1169 ns [45.0874, 45.1567] (0.95 confidence) | — | — |
| angle&#95;d6/dense&#95;varied/borrowed&#95;slice | 95.131 ns [94.9293, 95.2889] (0.95 confidence) | 95.7104 ns [95.4165, 95.8065] (0.95 confidence) | — | — |
| angle&#95;d6/dense&#95;varied/construct&#95;then&#95;angle | 92.7727 ns [92.6179, 92.9979] (0.95 confidence) | 93.1571 ns [92.9881, 93.4324] (0.95 confidence) | — | — |
| angle&#95;d6/dense&#95;varied/prepared&#95;vector | 91.5717 ns [91.3878, 91.7587] (0.95 confidence) | 92.1362 ns [91.928, 92.2566] (0.95 confidence) | — | — |
| angle&#95;d6/maximum&#95;small&#95;angle/borrowed&#95;slice | 15.1904 ns [15.1718, 15.1983] (0.95 confidence) | 15.9706 ns [15.8844, 15.9963] (0.95 confidence) | — | — |
| angle&#95;d6/maximum&#95;small&#95;angle/construct&#95;then&#95;angle | 14.2118 ns [14.1995, 14.2341] (0.95 confidence) | 14.3792 ns [14.289, 14.4212] (0.95 confidence) | — | — |
| angle&#95;d6/maximum&#95;small&#95;angle/prepared&#95;vector | 12.2232 ns [12.2134, 12.2719] (0.95 confidence) | 12.256 ns [12.2435, 12.2878] (0.95 confidence) | — | — |
| angle&#95;d6/mixed&#95;scale/borrowed&#95;slice | 16.7446 ns [16.7282, 16.7636] (0.95 confidence) | 17.2025 ns [17.1427, 17.2849] (0.95 confidence) | — | — |
| angle&#95;d6/mixed&#95;scale/construct&#95;then&#95;angle | 15.7944 ns [15.7216, 15.9482] (0.95 confidence) | 15.6041 ns [15.5794, 15.68] (0.95 confidence) | — | — |
| angle&#95;d6/mixed&#95;scale/prepared&#95;vector | 13.7596 ns [13.737, 13.7941] (0.95 confidence) | 13.5698 ns [13.5516, 13.5838] (0.95 confidence) | — | — |
| angle&#95;d6/near&#95;antipodal/borrowed&#95;slice | 21.1486 ns [21.1278, 21.1747] (0.95 confidence) | 21.0516 ns [21.0337, 21.0923] (0.95 confidence) | — | — |
| angle&#95;d6/near&#95;antipodal/construct&#95;then&#95;angle | 18.6321 ns [18.5774, 18.7122] (0.95 confidence) | 18.5846 ns [18.5584, 18.5989] (0.95 confidence) | — | — |
| angle&#95;d6/near&#95;antipodal/prepared&#95;vector | 17.3832 ns [17.3452, 17.4104] (0.95 confidence) | 17.2198 ns [17.194, 17.262] (0.95 confidence) | — | — |
| angle&#95;d6/near&#95;parallel/borrowed&#95;slice | 15.1722 ns [15.1587, 15.1787] (0.95 confidence) | 16.0034 ns [15.9873, 16.0291] (0.95 confidence) | — | — |
| angle&#95;d6/near&#95;parallel/construct&#95;then&#95;angle | 14.4783 ns [14.4323, 14.5347] (0.95 confidence) | 14.4786 ns [14.3847, 14.5464] (0.95 confidence) | — | — |
| angle&#95;d6/near&#95;parallel/prepared&#95;vector | 12.2019 ns [12.1938, 12.2109] (0.95 confidence) | 12.2288 ns [12.1971, 12.2718] (0.95 confidence) | — | — |
| angle&#95;d6/sparse&#95;adjacent/borrowed&#95;slice | 18.9213 ns [18.8824, 18.9687] (0.95 confidence) | 19.1039 ns [19.0677, 19.19] (0.95 confidence) | — | — |
| angle&#95;d6/sparse&#95;adjacent/construct&#95;then&#95;angle | 17.9229 ns [17.8243, 18.0052] (0.95 confidence) | 17.8411 ns [17.7795, 17.9659] (0.95 confidence) | — | — |
| angle&#95;d6/sparse&#95;adjacent/prepared&#95;vector | 16.0628 ns [16.0523, 16.0827] (0.95 confidence) | 15.8207 ns [15.7605, 15.934] (0.95 confidence) | — | — |
| angle&#95;d6/sparse&#95;alternating/borrowed&#95;slice | 29.443 ns [29.4079, 29.4589] (0.95 confidence) | 29.4867 ns [29.4454, 29.5116] (0.95 confidence) | — | — |
| angle&#95;d6/sparse&#95;alternating/construct&#95;then&#95;angle | 27.6185 ns [27.5987, 27.6372] (0.95 confidence) | 27.7882 ns [27.7555, 27.8241] (0.95 confidence) | — | — |
| angle&#95;d6/sparse&#95;alternating/prepared&#95;vector | 27.0362 ns [27.0165, 27.0557] (0.95 confidence) | 27.2606 ns [27.2142, 27.2741] (0.95 confidence) | — | — |
| angle&#95;d6/sparse&#95;antipodal/borrowed&#95;slice | 13.9363 ns [13.926, 13.9442] (0.95 confidence) | 14.3105 ns [14.2883, 14.3359] (0.95 confidence) | — | — |
| angle&#95;d6/sparse&#95;antipodal/construct&#95;then&#95;angle | 13.1698 ns [13.1658, 13.1804] (0.95 confidence) | 13.17 ns [13.1111, 13.2353] (0.95 confidence) | — | — |
| angle&#95;d6/sparse&#95;antipodal/prepared&#95;vector | 11.2125 ns [11.2046, 11.22] (0.95 confidence) | 11.0459 ns [11.0342, 11.058] (0.95 confidence) | — | — |
| angle&#95;d6/sparse&#95;identical/borrowed&#95;slice | 12.7177 ns [12.7025, 12.7446] (0.95 confidence) | 13.3133 ns [13.2877, 13.3306] (0.95 confidence) | — | — |
| angle&#95;d6/sparse&#95;identical/construct&#95;then&#95;angle | 11.8454 ns [11.7656, 11.8842] (0.95 confidence) | 12.0484 ns [11.9705, 12.119] (0.95 confidence) | — | — |
| angle&#95;d6/sparse&#95;identical/prepared&#95;vector | 9.97315 ns [9.9664, 9.9873] (0.95 confidence) | 10.0604 ns [10.0505, 10.0765] (0.95 confidence) | — | — |
| angle&#95;d6/subnormal/borrowed&#95;slice | 15.2634 ns [15.2388, 15.3174] (0.95 confidence) | 16.0641 ns [16.0024, 16.2084] (0.95 confidence) | — | — |
| angle&#95;d6/subnormal/construct&#95;then&#95;angle | 14.28 ns [14.2335, 14.339] (0.95 confidence) | 14.5698 ns [14.4742, 14.6161] (0.95 confidence) | — | — |
| angle&#95;d6/subnormal/prepared&#95;vector | 12.2481 ns [12.2416, 12.2646] (0.95 confidence) | 12.3158 ns [12.3018, 12.3274] (0.95 confidence) | — | — |
| angle&#95;d6/subnormal&#95;axes/borrowed&#95;slice | 18.4297 ns [18.4108, 18.458] (0.95 confidence) | 19.2845 ns [19.2361, 19.3496] (0.95 confidence) | — | — |
| angle&#95;d6/subnormal&#95;axes/construct&#95;then&#95;angle | 17.5007 ns [17.4691, 17.5641] (0.95 confidence) | 17.6434 ns [17.5799, 17.7156] (0.95 confidence) | — | — |
| angle&#95;d6/subnormal&#95;axes/prepared&#95;vector | 15.4567 ns [15.4498, 15.4706] (0.95 confidence) | 15.4722 ns [15.4422, 15.5063] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d2/rounded&#95;array | 169.378 ns [168.824, 170.092] (0.95 confidence) | 165.891 ns [165.508, 166.351] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d2/rounded&#95;result | 167.023 ns [166.129, 167.811] (0.95 confidence) | 166.161 ns [165.779, 166.677] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d2/strict&#95;array | 169.012 ns [168.591, 169.516] (0.95 confidence) | 169.97 ns [169.802, 170.141] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d2/strict&#95;result | 168.86 ns [168.162, 169.377] (0.95 confidence) | 171.941 ns [171.574, 172.415] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d3/rounded&#95;array | 172.14 ns [171.758, 172.466] (0.95 confidence) | 168.979 ns [168.795, 169.175] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d3/rounded&#95;result | 170.788 ns [170.516, 170.972] (0.95 confidence) | 169.317 ns [169.15, 169.57] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d3/strict&#95;array | 176.169 ns [175.656, 176.688] (0.95 confidence) | 170.422 ns [170.219, 170.836] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d3/strict&#95;result | 173.056 ns [172.268, 173.679] (0.95 confidence) | 169.149 ns [168.947, 169.412] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d4/rounded&#95;array | 203.478 ns [203.057, 204.274] (0.95 confidence) | 176.864 ns [176.589, 177.26] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d4/rounded&#95;result | 203.84 ns [203.297, 204.813] (0.95 confidence) | 176.545 ns [176.193, 176.759] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d4/strict&#95;array | 199.587 ns [199.083, 200.349] (0.95 confidence) | 175.283 ns [174.999, 176.718] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d4/strict&#95;result | 202.405 ns [200.435, 203.337] (0.95 confidence) | 174.597 ns [174.074, 175.425] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d5/rounded&#95;array | 183.846 ns [183.203, 184.509] (0.95 confidence) | 180.881 ns [180.603, 181.173] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d5/rounded&#95;result | 182.25 ns [181.772, 182.947] (0.95 confidence) | 181.382 ns [181.159, 181.806] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d5/strict&#95;array | 180.066 ns [179.629, 180.797] (0.95 confidence) | 179.773 ns [179.13, 180.759] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;below&#95;overflow&#95;d5/strict&#95;result | 178.218 ns [177.865, 178.629] (0.95 confidence) | 176.314 ns [175.965, 176.604] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d2/rounded&#95;array | 12.2877 ns [12.2651, 12.301] (0.95 confidence) | 11.5532 ns [11.5511, 11.5554] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d2/rounded&#95;result | 12.2297 ns [12.2177, 12.2387] (0.95 confidence) | 11.5555 ns [11.5518, 11.5608] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d2/strict&#95;array | 9.75463 ns [9.73719, 9.76636] (0.95 confidence) | 8.59831 ns [8.42609, 8.88843] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d2/strict&#95;result | 8.48447 ns [8.42224, 8.53416] (0.95 confidence) | 7.79616 ns [7.79276, 7.84795] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d3/rounded&#95;array | 17.9951 ns [17.9867, 18.0086] (0.95 confidence) | 17.8624 ns [17.8565, 17.8665] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d3/rounded&#95;result | 17.998 ns [17.9866, 18.0097] (0.95 confidence) | 17.8559 ns [17.8529, 17.861] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d3/strict&#95;array | 13.9358 ns [13.8799, 13.9876] (0.95 confidence) | 12.3034 ns [12.1734, 12.343] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d3/strict&#95;result | 12.1254 ns [11.8604, 12.23] (0.95 confidence) | 11.279 ns [11.2614, 11.3355] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d4/rounded&#95;array | 23.862 ns [23.843, 23.8789] (0.95 confidence) | 23.7126 ns [23.7016, 23.7238] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d4/rounded&#95;result | 23.8515 ns [23.8417, 23.8735] (0.95 confidence) | 23.7075 ns [23.6874, 23.7281] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d4/strict&#95;array | 18.262 ns [18.1977, 18.3442] (0.95 confidence) | 16.8885 ns [16.2384, 17.2784] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d4/strict&#95;result | 15.4468 ns [15.4356, 15.5168] (0.95 confidence) | 15.3426 ns [15.3374, 15.3659] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d5/rounded&#95;array | 29.1374 ns [29.1294, 29.1527] (0.95 confidence) | 28.4286 ns [28.4204, 28.4346] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d5/rounded&#95;result | 28.8038 ns [28.7938, 28.8244] (0.95 confidence) | 28.4442 ns [28.4386, 28.4562] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d5/strict&#95;array | 23.0406 ns [22.1819, 23.1708] (0.95 confidence) | 22.1504 ns [21.8057, 22.9275] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;dyadic&#95;d5/strict&#95;result | 20.2851 ns [19.804, 21.1798] (0.95 confidence) | 19.0669 ns [19.0424, 19.2579] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d2/rounded&#95;array | 172.403 ns [172.13, 172.69] (0.95 confidence) | 170.125 ns [169.937, 170.299] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d2/rounded&#95;result | 172.288 ns [171.972, 172.439] (0.95 confidence) | 170.194 ns [169.862, 170.453] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d2/strict&#95;array | 14.0716 ns [14.0511, 14.1004] (0.95 confidence) | 13.9124 ns [13.7865, 14.0183] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d2/strict&#95;result | 12.9003 ns [12.7774, 12.9901] (0.95 confidence) | 12.4771 ns [12.3361, 12.6766] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d3/rounded&#95;array | 178.738 ns [178.463, 178.915] (0.95 confidence) | 176.595 ns [176.438, 176.803] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d3/rounded&#95;result | 178.474 ns [178.318, 178.706] (0.95 confidence) | 177.705 ns [177.482, 177.992] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d3/strict&#95;array | 18.6058 ns [18.512, 18.7083] (0.95 confidence) | 17.7051 ns [17.5992, 17.7709] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d3/strict&#95;result | 17.5039 ns [17.4299, 17.6555] (0.95 confidence) | 16.3976 ns [16.2739, 16.631] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d4/rounded&#95;array | 185.35 ns [185.193, 185.505] (0.95 confidence) | 183.491 ns [183.281, 183.798] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d4/rounded&#95;result | 184.627 ns [184.343, 185.022] (0.95 confidence) | 184.648 ns [184.413, 184.963] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d4/strict&#95;array | 23.4507 ns [23.3272, 23.535] (0.95 confidence) | 21.9484 ns [21.9277, 21.9632] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d4/strict&#95;result | 22.0545 ns [21.9821, 22.2184] (0.95 confidence) | 20.0631 ns [19.9017, 20.6028] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d5/rounded&#95;array | 190.642 ns [190.302, 190.826] (0.95 confidence) | 187.828 ns [187.622, 187.926] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d5/rounded&#95;result | 190.569 ns [190.24, 190.723] (0.95 confidence) | 188.331 ns [188.103, 188.722] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d5/strict&#95;array | 28.0909 ns [28.0181, 28.1515] (0.95 confidence) | 26.5173 ns [26.3884, 26.6547] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;min&#95;subnormal&#95;d5/strict&#95;result | 26.1995 ns [26.0678, 26.2965] (0.95 confidence) | 24.2077 ns [23.8395, 24.6699] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d2/rounded&#95;array | 171.749 ns [171.398, 171.962] (0.95 confidence) | 169.682 ns [169.564, 169.795] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d2/rounded&#95;result | 172.256 ns [171.889, 172.522] (0.95 confidence) | 170.012 ns [169.723, 170.264] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d2/strict&#95;array | 13.3091 ns [13.2762, 13.4215] (0.95 confidence) | 12.7135 ns [12.6793, 12.8035] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d2/strict&#95;result | 12.663 ns [12.6449, 12.6838] (0.95 confidence) | 12.6491 ns [12.4276, 12.7339] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d3/rounded&#95;array | 178.648 ns [178.452, 178.777] (0.95 confidence) | 176.841 ns [176.474, 176.963] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d3/rounded&#95;result | 178.049 ns [177.759, 178.169] (0.95 confidence) | 176.69 ns [176.435, 176.931] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d3/strict&#95;array | 17.7578 ns [17.6344, 17.8569] (0.95 confidence) | 16.8292 ns [16.732, 16.9119] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d3/strict&#95;result | 16.7631 ns [16.6985, 16.9878] (0.95 confidence) | 16.5154 ns [16.4249, 16.58] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d4/rounded&#95;array | 189.212 ns [187.154, 193.553] (0.95 confidence) | 183.592 ns [183.417, 183.896] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d4/rounded&#95;result | 184.216 ns [184.006, 184.414] (0.95 confidence) | 183.685 ns [183.369, 183.935] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d4/strict&#95;array | 22.4192 ns [22.319, 22.488] (0.95 confidence) | 21.5198 ns [21.4178, 21.6459] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d4/strict&#95;result | 20.5687 ns [20.4332, 20.6915] (0.95 confidence) | 20.0344 ns [19.8551, 20.1657] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d5/rounded&#95;array | 189.492 ns [189.236, 189.693] (0.95 confidence) | 189.489 ns [189.33, 190.122] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d5/rounded&#95;result | 189.224 ns [188.88, 189.517] (0.95 confidence) | 188.353 ns [188.029, 188.672] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d5/strict&#95;array | 26.7939 ns [26.7243, 26.8464] (0.95 confidence) | 24.2758 ns [24.2539, 25.0146] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;negative&#95;underflow&#95;d5/strict&#95;result | 25.0873 ns [24.9481, 25.1604] (0.95 confidence) | 23.2088 ns [23.1368, 23.2788] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d2/rounded&#95;array | 12.0307 ns [12.0251, 12.0393] (0.95 confidence) | 11.5762 ns [11.5687, 11.5807] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d2/rounded&#95;result | 12.3277 ns [12.302, 12.3506] (0.95 confidence) | 11.5755 ns [11.567, 11.5835] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d2/strict&#95;array | 56.9936 ns [56.85, 57.1777] (0.95 confidence) | 55.6388 ns [55.5246, 55.8321] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d2/strict&#95;result | 12.1476 ns [12.1408, 12.1589] (0.95 confidence) | 11.7443 ns [11.7336, 11.7526] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d3/rounded&#95;array | 18.0387 ns [18.0315, 18.0472] (0.95 confidence) | 17.9734 ns [17.9493, 17.9901] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d3/rounded&#95;result | 18.0348 ns [18.0103, 18.0518] (0.95 confidence) | 17.938 ns [17.9253, 17.9495] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d3/strict&#95;array | 60.8547 ns [60.7606, 60.9546] (0.95 confidence) | 61.2297 ns [61.0908, 61.4146] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d3/strict&#95;result | 15.9397 ns [15.9186, 15.965] (0.95 confidence) | 15.487 ns [15.4532, 15.5009] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d4/rounded&#95;array | 23.9691 ns [23.9588, 23.9814] (0.95 confidence) | 23.8014 ns [23.7898, 23.8084] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d4/rounded&#95;result | 23.9768 ns [23.9652, 24.0001] (0.95 confidence) | 23.8213 ns [23.8065, 23.8302] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d4/strict&#95;array | 63.9741 ns [63.8285, 64.1235] (0.95 confidence) | 63.3716 ns [63.2268, 63.6956] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d4/strict&#95;result | 19.3051 ns [19.2161, 19.3963] (0.95 confidence) | 19.4249 ns [19.3849, 19.4597] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d5/rounded&#95;array | 28.5349 ns [28.5144, 28.5457] (0.95 confidence) | 28.3595 ns [28.3543, 28.3655] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d5/rounded&#95;result | 29.0914 ns [29.0748, 29.1029] (0.95 confidence) | 28.3465 ns [28.3403, 28.3542] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d5/strict&#95;array | 68.9066 ns [68.7961, 69.1164] (0.95 confidence) | 67.9443 ns [67.7418, 68.2523] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;non&#95;dyadic&#95;d5/strict&#95;result | 23.021 ns [22.9906, 23.0554] (0.95 confidence) | 22.1602 ns [22.1442, 22.2904] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d2/rounded&#95;array | 167.16 ns [166.495, 167.82] (0.95 confidence) | 165.63 ns [165.503, 165.817] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d2/rounded&#95;result | 166.313 ns [166.003, 166.575] (0.95 confidence) | 166.212 ns [165.656, 166.717] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d2/strict&#95;array | 176.613 ns [176.177, 177.298] (0.95 confidence) | 175.11 ns [174.736, 175.344] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d2/strict&#95;result | 175.644 ns [175.157, 176.204] (0.95 confidence) | 174.546 ns [173.934, 175.284] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d3/rounded&#95;array | 174.254 ns [174.019, 174.782] (0.95 confidence) | 171.629 ns [171.329, 172.227] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d3/rounded&#95;result | 172.336 ns [172.173, 172.721] (0.95 confidence) | 171.925 ns [171.562, 172.165] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d3/strict&#95;array | 182.255 ns [181.63, 182.658] (0.95 confidence) | 180.373 ns [179.841, 180.935] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d3/strict&#95;result | 178.338 ns [178.095, 178.611] (0.95 confidence) | 180.463 ns [180.042, 181.026] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d4/rounded&#95;array | 184.546 ns [183.576, 185.054] (0.95 confidence) | 179.536 ns [179.077, 179.82] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d4/rounded&#95;result | 181.728 ns [181.15, 182.086] (0.95 confidence) | 178.393 ns [178.271, 178.774] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d4/strict&#95;array | 188.424 ns [187.754, 189.424] (0.95 confidence) | 184.35 ns [183.688, 185.185] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d4/strict&#95;result | 186.59 ns [186.196, 187.531] (0.95 confidence) | 182.734 ns [182.059, 183.431] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d5/rounded&#95;array | 184.726 ns [184.153, 185.159] (0.95 confidence) | 183.621 ns [183.421, 183.9] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d5/rounded&#95;result | 183.272 ns [183.089, 183.523] (0.95 confidence) | 184.407 ns [183.952, 185.004] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d5/strict&#95;array | 188.928 ns [188.217, 189.537] (0.95 confidence) | 187.49 ns [186.976, 187.935] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d5/strict&#95;result | 187.615 ns [185.771, 188.796] (0.95 confidence) | 182.794 ns [182.647, 183.086] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d2/rounded&#95;array | 198.021 ns [197.827, 198.582] (0.95 confidence) | 201.167 ns [200.818, 202.081] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d2/rounded&#95;result | 200.614 ns [199.806, 201.225] (0.95 confidence) | 202.466 ns [202.077, 202.858] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d2/strict&#95;array | 13494.7 ns [13471.7, 13510.8] (0.95 confidence) | 14055.4 ns [14034.7, 14096.2] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d2/strict&#95;result | 201.084 ns [200.544, 201.388] (0.95 confidence) | 205.64 ns [205.06, 206.348] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d3/rounded&#95;array | 204.541 ns [204.071, 205.306] (0.95 confidence) | 206.076 ns [205.64, 206.468] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d3/rounded&#95;result | 204.661 ns [204.367, 205.003] (0.95 confidence) | 207.965 ns [207.5, 208.207] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d3/strict&#95;array | 13435.4 ns [13413.8, 13456.6] (0.95 confidence) | 14050.9 ns [14035.2, 14083.6] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d3/strict&#95;result | 205.719 ns [205.211, 206.288] (0.95 confidence) | 207.971 ns [207.611, 208.3] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d4/rounded&#95;array | 210.716 ns [210.481, 211.171] (0.95 confidence) | 215.464 ns [215.097, 215.768] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d4/rounded&#95;result | 209.186 ns [208.994, 209.537] (0.95 confidence) | 216.426 ns [215.986, 217.337] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d4/strict&#95;array | 13436.1 ns [13389.7, 13452.1] (0.95 confidence) | 14070.9 ns [14050.8, 14102.1] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d4/strict&#95;result | 206.736 ns [206.455, 207.022] (0.95 confidence) | 209.7 ns [209.033, 210.144] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d5/rounded&#95;array | 215.696 ns [215.48, 215.966] (0.95 confidence) | 221.015 ns [220.663, 221.66] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d5/rounded&#95;result | 217.709 ns [217.451, 217.959] (0.95 confidence) | 221.522 ns [221.035, 221.951] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d5/strict&#95;array | 13471.1 ns [13454.9, 13487.7] (0.95 confidence) | 14118.9 ns [14102.4, 14158.2] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide1024&#95;d5/strict&#95;result | 213.251 ns [212.714, 213.721] (0.95 confidence) | 215.15 ns [214.948, 215.484] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d2/rounded&#95;array | 160.53 ns [159.872, 161.019] (0.95 confidence) | 157.591 ns [157.072, 158.76] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d2/rounded&#95;result | 162.982 ns [161.965, 163.637] (0.95 confidence) | 161.389 ns [160.559, 162.04] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d2/strict&#95;array | 3860.63 ns [3828.75, 3878.25] (0.95 confidence) | 3876.37 ns [3847.66, 3899.19] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d2/strict&#95;result | 166.29 ns [165.833, 166.784] (0.95 confidence) | 171.431 ns [169.683, 172.491] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d3/rounded&#95;array | 168.741 ns [168.232, 169.66] (0.95 confidence) | 164.563 ns [164.324, 164.741] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d3/rounded&#95;result | 167.753 ns [167.113, 168.442] (0.95 confidence) | 162.917 ns [162.537, 163.152] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d3/strict&#95;array | 3834.76 ns [3778.53, 3848.65] (0.95 confidence) | 3771 ns [3755.86, 3780.67] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d3/strict&#95;result | 167.201 ns [166.481, 168.518] (0.95 confidence) | 162.825 ns [162.468, 163.137] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d4/rounded&#95;array | 174.906 ns [174.395, 175.179] (0.95 confidence) | 169.066 ns [168.761, 169.308] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d4/rounded&#95;result | 174.283 ns [173.787, 174.791] (0.95 confidence) | 176.042 ns [174.316, 178.073] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d4/strict&#95;array | 3769.04 ns [3744.89, 3784.82] (0.95 confidence) | 3909.96 ns [3886.34, 3923.73] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d4/strict&#95;result | 171.467 ns [171.084, 172.221] (0.95 confidence) | 176.382 ns [174.201, 177.346] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d5/rounded&#95;array | 179.384 ns [178.553, 179.87] (0.95 confidence) | 175.73 ns [175.196, 176.771] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d5/rounded&#95;result | 181.01 ns [180.462, 181.735] (0.95 confidence) | 180.065 ns [179.584, 180.441] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d5/strict&#95;array | 3844.34 ns [3825.6, 3876.81] (0.95 confidence) | 3829.63 ns [3786.12, 3858.66] (0.95 confidence) | — | — |
| canonical&#95;conversion&#95;wide256&#95;d5/strict&#95;result | 176.803 ns [176.017, 178.471] (0.95 confidence) | 169.455 ns [169.233, 170.272] (0.95 confidence) | — | — |
| column&#95;basis/2x2/full&#95;rank | — | 4.87551 ns [4.86459, 4.89447] (0.95 confidence) | — | — |
| column&#95;basis/2x2/large | — | 4.88566 ns [4.87448, 4.89693] (0.95 confidence) | — | — |
| column&#95;basis/2x2/near&#95;dependent | — | 4.9299 ns [4.91518, 4.94344] (0.95 confidence) | — | — |
| column&#95;basis/2x2/rank&#95;one | — | 4.68395 ns [4.67386, 4.68901] (0.95 confidence) | — | — |
| column&#95;basis/2x2/small | — | 4.89218 ns [4.88609, 4.90233] (0.95 confidence) | — | — |
| column&#95;basis/2x3/shared&#95;face | — | 9.91512 ns [9.8953, 9.93527] (0.95 confidence) | — | — |
| column&#95;basis/3x3/full&#95;rank | — | 13.3705 ns [13.3477, 13.4102] (0.95 confidence) | — | — |
| column&#95;basis/3x3/large | — | 13.3977 ns [13.3753, 13.4356] (0.95 confidence) | — | — |
| column&#95;basis/3x3/near&#95;dependent | — | 13.4564 ns [13.4291, 13.4799] (0.95 confidence) | — | — |
| column&#95;basis/3x3/rank&#95;one | — | 9.20716 ns [9.16816, 9.23003] (0.95 confidence) | — | — |
| column&#95;basis/3x3/small | — | 13.4676 ns [13.4213, 13.5087] (0.95 confidence) | — | — |
| column&#95;basis/3x4/shared&#95;face | — | 17.4341 ns [17.4185, 17.4646] (0.95 confidence) | — | — |
| column&#95;basis/4x4/full&#95;rank | — | 29.1789 ns [28.9557, 29.328] (0.95 confidence) | — | — |
| column&#95;basis/4x4/large | — | 29.2383 ns [28.9913, 29.4059] (0.95 confidence) | — | — |
| column&#95;basis/4x4/near&#95;dependent | — | 29.3965 ns [29.3006, 29.5259] (0.95 confidence) | — | — |
| column&#95;basis/4x4/rank&#95;one | — | 14.5871 ns [14.5632, 14.6004] (0.95 confidence) | — | — |
| column&#95;basis/4x4/small | — | 29.1842 ns [28.9803, 29.2867] (0.95 confidence) | — | — |
| column&#95;basis/4x5/shared&#95;face | — | 37.497 ns [37.2911, 38.6217] (0.95 confidence) | — | — |
| column&#95;basis/4x6/full&#95;rank | — | 36.0818 ns [36.0167, 36.1567] (0.95 confidence) | — | — |
| column&#95;basis/4x6/large | — | 35.9234 ns [35.8813, 36.0778] (0.95 confidence) | — | — |
| column&#95;basis/4x6/near&#95;dependent | — | 35.7584 ns [35.6844, 35.8516] (0.95 confidence) | — | — |
| column&#95;basis/4x6/rank&#95;one | — | 20.7802 ns [20.7217, 20.8337] (0.95 confidence) | — | — |
| column&#95;basis/4x6/small | — | 36.0075 ns [35.9145, 36.103] (0.95 confidence) | — | — |
| column&#95;basis/5x2/full&#95;rank | — | 13.2423 ns [13.2054, 13.2648] (0.95 confidence) | — | — |
| column&#95;basis/5x2/large | — | 13.203 ns [13.1897, 13.2249] (0.95 confidence) | — | — |
| column&#95;basis/5x2/near&#95;dependent | — | 13.1183 ns [13.0821, 13.1711] (0.95 confidence) | — | — |
| column&#95;basis/5x2/rank&#95;one | — | 11.0461 ns [10.9974, 11.0867] (0.95 confidence) | — | — |
| column&#95;basis/5x2/small | — | 13.1912 ns [13.1568, 13.2164] (0.95 confidence) | — | — |
| column&#95;basis/5x5/full&#95;rank | — | 51.3191 ns [51.1665, 51.4943] (0.95 confidence) | — | — |
| column&#95;basis/5x5/large | — | 51.353 ns [51.2417, 51.4358] (0.95 confidence) | — | — |
| column&#95;basis/5x5/near&#95;dependent | — | 51.1935 ns [51.019, 51.3286] (0.95 confidence) | — | — |
| column&#95;basis/5x5/rank&#95;one | — | 23.7538 ns [23.721, 23.8057] (0.95 confidence) | — | — |
| column&#95;basis/5x5/small | — | 51.291 ns [51.2142, 51.4959] (0.95 confidence) | — | — |
| column&#95;basis/5x6/shared&#95;face | — | 61.9312 ns [61.709, 62.2505] (0.95 confidence) | — | — |
| column&#95;basis/5x8/full&#95;rank | — | 61.8642 ns [61.7462, 61.9333] (0.95 confidence) | — | — |
| column&#95;basis/5x8/large | — | 61.8465 ns [61.7682, 61.9365] (0.95 confidence) | — | — |
| column&#95;basis/5x8/near&#95;dependent | — | 61.7506 ns [61.6986, 61.8052] (0.95 confidence) | — | — |
| column&#95;basis/5x8/rank&#95;one | — | 37.3936 ns [36.5934, 37.6122] (0.95 confidence) | — | — |
| column&#95;basis/5x8/small | — | 61.6877 ns [61.5963, 61.769] (0.95 confidence) | — | — |
| column&#95;basis/6x10/full&#95;rank | — | 97.8628 ns [97.7894, 97.9471] (0.95 confidence) | — | — |
| column&#95;basis/6x10/large | — | 97.6644 ns [97.5779, 97.7826] (0.95 confidence) | — | — |
| column&#95;basis/6x10/near&#95;dependent | — | 97.0547 ns [96.9728, 97.153] (0.95 confidence) | — | — |
| column&#95;basis/6x10/rank&#95;one | — | 55.7397 ns [55.6868, 55.9252] (0.95 confidence) | — | — |
| column&#95;basis/6x10/small | — | 97.7712 ns [97.6431, 97.8339] (0.95 confidence) | — | — |
| column&#95;basis/6x7/shared&#95;face | — | 97.7082 ns [97.3705, 98.0467] (0.95 confidence) | — | — |
| column&#95;basis/7x12/full&#95;rank | — | 142.518 ns [142.36, 142.673] (0.95 confidence) | — | — |
| column&#95;basis/7x12/large | — | 142.67 ns [142.512, 142.779] (0.95 confidence) | — | — |
| column&#95;basis/7x12/near&#95;dependent | — | 143.146 ns [142.91, 143.333] (0.95 confidence) | — | — |
| column&#95;basis/7x12/rank&#95;one | — | 78.8162 ns [78.6967, 78.9118] (0.95 confidence) | — | — |
| column&#95;basis/7x12/small | — | 142.866 ns [142.632, 143.04] (0.95 confidence) | — | — |
| column&#95;basis/8x14/full&#95;rank | — | 197.067 ns [196.679, 197.228] (0.95 confidence) | — | — |
| column&#95;basis/8x14/large | — | 197.406 ns [197.149, 197.678] (0.95 confidence) | — | — |
| column&#95;basis/8x14/near&#95;dependent | — | 196.915 ns [196.619, 197.087] (0.95 confidence) | — | — |
| column&#95;basis/8x14/rank&#95;one | — | 103.852 ns [103.782, 103.928] (0.95 confidence) | — | — |
| column&#95;basis/8x14/small | — | 196.882 ns [196.598, 197.267] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;det | 400.475 ns [397.654, 403.477] (0.95 confidence) | 396.167 ns [394.043, 397.739] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;det&#95;from&#95;ldlt | 2.92958 ns [2.92534, 2.93206] (0.95 confidence) | 2.92139 ns [2.91715, 2.92599] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;det&#95;from&#95;lu | 3.32053 ns [3.31634, 3.3237] (0.95 confidence) | 3.29684 ns [3.29292, 3.30236] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;det&#95;via&#95;lu | 383.544 ns [380.335, 385.779] (0.95 confidence) | 414.369 ns [407.784, 419.823] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;dot | 2.36457 ns [2.35992, 2.36894] (0.95 confidence) | 2.34045 ns [2.33499, 2.34976] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;inf&#95;norm | 33.2031 ns [33.1403, 33.2379] (0.95 confidence) | 32.6599 ns [32.6293, 32.6991] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;ldlt | 341.89 ns [341.43, 342.129] (0.95 confidence) | 344.564 ns [344.081, 344.976] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;ldlt&#95;solve | 383.805 ns [383.422, 384.206] (0.95 confidence) | 386.388 ns [385.841, 386.863] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;lu | 365.836 ns [364.153, 366.937] (0.95 confidence) | 406.273 ns [403.113, 409.5] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;lu&#95;solve | 633.163 ns [631.635, 634.907] (0.95 confidence) | 631.9 ns [628.777, 638.18] (0.95 confidence) | 589.967 ns [588.536, 591.405] (0.95 confidence) | 900.746 ns [897.01, 902.971] (0.95 confidence) |
| d16/la&#95;stack&#95;norm2 | 10.5581 ns [10.4871, 10.6205] (0.95 confidence) | 10.2881 ns [10.1609, 10.3694] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;norm2&#95;scenario&#95;descending | 9.67065 ns [9.64515, 9.68597] (0.95 confidence) | 9.48122 ns [9.44822, 9.5047] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale | 9.61287 ns [9.59129, 9.64811] (0.95 confidence) | 9.49346 ns [9.46517, 9.53399] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;norm2&#95;scenario&#95;sparse | 5.05144 ns [4.96285, 5.1092] (0.95 confidence) | 4.42447 ns [4.40139, 4.4529] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range | 8.58392 ns [8.57245, 8.5956] (0.95 confidence) | 8.33589 ns [8.31522, 8.3675] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;norm2&#95;sq | 2.13582 ns [2.13158, 2.14127] (0.95 confidence) | 2.10971 ns [2.10409, 2.11526] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;solve&#95;from&#95;ldlt | 28.1214 ns [28.1, 28.1602] (0.95 confidence) | 28.0231 ns [27.9648, 28.056] (0.95 confidence) | — | — |
| d16/la&#95;stack&#95;solve&#95;from&#95;lu | 194.418 ns [193.64, 195.02] (0.95 confidence) | 192.263 ns [191.781, 192.872] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;det | 0.577625 ns [0.576939, 0.578455] (0.95 confidence) | 0.576858 ns [0.576023, 0.578265] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;det&#95;from&#95;ldlt | 0.50138 ns [0.50102, 0.501804] (0.95 confidence) | 0.499765 ns [0.498951, 0.501182] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;det&#95;from&#95;lu | 0.552009 ns [0.550757, 0.553512] (0.95 confidence) | 0.525001 ns [0.524099, 0.525605] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;det&#95;via&#95;lu | 1.77187 ns [1.76941, 1.77368] (0.95 confidence) | 1.5981 ns [1.59579, 1.60064] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;dot | 0.616439 ns [0.61485, 0.618042] (0.95 confidence) | 0.6185 ns [0.617968, 0.619101] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;inf&#95;norm | 0.709753 ns [0.708699, 0.710384] (0.95 confidence) | 0.708593 ns [0.707998, 0.709845] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;ldlt | 6.93545 ns [6.9252, 6.94391] (0.95 confidence) | 6.60203 ns [6.57931, 6.61506] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;ldlt&#95;solve | 9.87096 ns [9.85544, 9.91457] (0.95 confidence) | 9.50609 ns [9.48448, 9.52406] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;lu | 1.57001 ns [1.56776, 1.57254] (0.95 confidence) | 1.39947 ns [1.39661, 1.40236] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;lu&#95;solve | 2.06081 ns [2.05195, 2.06951] (0.95 confidence) | 2.06256 ns [2.04953, 2.09431] (0.95 confidence) | 4.52663 ns [4.51758, 4.53267] (0.95 confidence) | 128.125 ns [127.445, 128.545] (0.95 confidence) |
| d2/la&#95;stack&#95;norm2 | 1.24263 ns [1.23885, 1.24911] (0.95 confidence) | 1.16126 ns [1.15994, 1.16262] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;norm2&#95;scenario&#95;descending | 1.36543 ns [1.36341, 1.36826] (0.95 confidence) | 1.35695 ns [1.35287, 1.35765] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale | 1.3638 ns [1.36068, 1.36911] (0.95 confidence) | 1.3553 ns [1.35379, 1.35818] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;norm2&#95;scenario&#95;sparse | 1.16722 ns [1.16379, 1.16902] (0.95 confidence) | 1.19988 ns [1.16606, 1.22116] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range | 1.36678 ns [1.36374, 1.37196] (0.95 confidence) | 1.35398 ns [1.35006, 1.35599] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;norm2&#95;sq | 0.404173 ns [0.402338, 0.40648] (0.95 confidence) | 0.394665 ns [0.393848, 0.395844] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;solve&#95;from&#95;ldlt | 1.2457 ns [1.24222, 1.24895] (0.95 confidence) | 1.24439 ns [1.24261, 1.24524] (0.95 confidence) | — | — |
| d2/la&#95;stack&#95;solve&#95;from&#95;lu | 1.28669 ns [1.28277, 1.29053] (0.95 confidence) | 1.28359 ns [1.28103, 1.28718] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;det | 1.35695 ns [1.35124, 1.3647] (0.95 confidence) | 1.35047 ns [1.34662, 1.36171] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;det&#95;from&#95;ldlt | 0.59751 ns [0.597163, 0.597679] (0.95 confidence) | 0.600426 ns [0.600053, 0.600826] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;det&#95;from&#95;lu | 0.71243 ns [0.710531, 0.715423] (0.95 confidence) | 0.700692 ns [0.700007, 0.701188] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;det&#95;via&#95;lu | 8.86687 ns [8.8499, 8.87945] (0.95 confidence) | 8.78512 ns [8.76191, 8.80381] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;dot | 0.671561 ns [0.670471, 0.672575] (0.95 confidence) | 0.673201 ns [0.672208, 0.67457] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;inf&#95;norm | 1.37394 ns [1.3719, 1.37625] (0.95 confidence) | 1.37344 ns [1.37039, 1.38758] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;ldlt | 13.1726 ns [13.1387, 13.2355] (0.95 confidence) | 8.57796 ns [8.56247, 8.59357] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;ldlt&#95;solve | 9.64226 ns [9.62827, 9.66186] (0.95 confidence) | 9.42105 ns [9.39126, 9.45339] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;lu | 8.44357 ns [8.43179, 8.45519] (0.95 confidence) | 8.4172 ns [8.39453, 8.4354] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;lu&#95;solve | 9.89299 ns [9.87485, 9.90579] (0.95 confidence) | 9.98278 ns [9.96279, 10.0083] (0.95 confidence) | 32.5058 ns [32.3949, 32.6237] (0.95 confidence) | 156.551 ns [156.195, 157.093] (0.95 confidence) |
| d3/la&#95;stack&#95;norm2 | 1.95482 ns [1.86494, 2.00939] (0.95 confidence) | 1.78608 ns [1.77909, 1.79287] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;norm2&#95;scenario&#95;descending | 1.85721 ns [1.85238, 1.85984] (0.95 confidence) | 1.85068 ns [1.84616, 1.85519] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale | 1.85213 ns [1.84822, 1.85795] (0.95 confidence) | 1.85131 ns [1.84539, 1.8542] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;norm2&#95;scenario&#95;sparse | 1.43864 ns [1.40636, 1.47307] (0.95 confidence) | 1.46052 ns [1.43373, 1.47018] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range | 1.85221 ns [1.84873, 1.8567] (0.95 confidence) | 1.84853 ns [1.84243, 1.85253] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;norm2&#95;sq | 0.44311 ns [0.441832, 0.444752] (0.95 confidence) | 0.441195 ns [0.44033, 0.441916] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;solve&#95;from&#95;ldlt | 1.79878 ns [1.79623, 1.80177] (0.95 confidence) | 1.79838 ns [1.79638, 1.80019] (0.95 confidence) | — | — |
| d3/la&#95;stack&#95;solve&#95;from&#95;lu | 2.09039 ns [2.08569, 2.09414] (0.95 confidence) | 2.07388 ns [2.0717, 2.0763] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;det | 2219.73 ns [2214.84, 2222.86] (0.95 confidence) | 2040.12 ns [2032.03, 2044.49] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;det&#95;from&#95;ldlt | 8.08075 ns [8.06844, 8.09245] (0.95 confidence) | 7.8484 ns [7.82629, 7.85882] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;det&#95;from&#95;lu | 9.33744 ns [9.3095, 9.38465] (0.95 confidence) | 9.16064 ns [9.14285, 9.17761] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;det&#95;via&#95;lu | 2077.12 ns [2070.61, 2080.62] (0.95 confidence) | 2059.15 ns [2054.3, 2064.37] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;dot | 4.0975 ns [4.09436, 4.10111] (0.95 confidence) | 4.05523 ns [4.04856, 4.06384] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;inf&#95;norm | 131.515 ns [131.353, 131.676] (0.95 confidence) | 129.469 ns [129.178, 129.69] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;ldlt | 2435.31 ns [2429.34, 2438.19] (0.95 confidence) | 2415.42 ns [2412.27, 2418.57] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;ldlt&#95;solve | 2766.55 ns [2762.23, 2772.97] (0.95 confidence) | 2738.41 ns [2733.31, 2742.8] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;lu | 2322.96 ns [2292.69, 2354.57] (0.95 confidence) | 2284.71 ns [2273.05, 2325.63] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;lu&#95;solve | 2848.06 ns [2835.22, 2858.39] (0.95 confidence) | 2772.06 ns [2766.99, 2779.6] (0.95 confidence) | 2499.64 ns [2493.29, 2506.42] (0.95 confidence) | 2998.91 ns [2995.26, 3003.01] (0.95 confidence) |
| d32/la&#95;stack&#95;norm2 | 24.5467 ns [24.4978, 24.6199] (0.95 confidence) | 23.5342 ns [23.4784, 23.5988] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;norm2&#95;scenario&#95;descending | 19.8387 ns [19.7855, 19.8717] (0.95 confidence) | 18.761 ns [18.7244, 18.8219] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale | 19.5915 ns [19.5335, 19.6313] (0.95 confidence) | 18.7325 ns [18.7024, 18.7784] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;norm2&#95;scenario&#95;sparse | 9.84375 ns [9.81066, 9.86317] (0.95 confidence) | 9.2909 ns [9.27892, 9.31593] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range | 17.1293 ns [17.1013, 17.1897] (0.95 confidence) | 16.3483 ns [16.3123, 16.3747] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;norm2&#95;sq | 4.07092 ns [4.06396, 4.08385] (0.95 confidence) | 4.01816 ns [4.00858, 4.02468] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;solve&#95;from&#95;ldlt | 319.109 ns [318.64, 319.543] (0.95 confidence) | 314.7 ns [313.751, 315.227] (0.95 confidence) | — | — |
| d32/la&#95;stack&#95;solve&#95;from&#95;lu | 717.625 ns [715.93, 718.883] (0.95 confidence) | 672.538 ns [671.207, 673.725] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;det | 2.57718 ns [2.56736, 2.58353] (0.95 confidence) | 2.55777 ns [2.53254, 2.56702] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;det&#95;from&#95;ldlt | 0.774725 ns [0.77374, 0.775741] (0.95 confidence) | 0.772085 ns [0.771759, 0.772327] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;det&#95;from&#95;lu | 0.877497 ns [0.875988, 0.879378] (0.95 confidence) | 0.871482 ns [0.871054, 0.871644] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;det&#95;via&#95;lu | 14.2642 ns [14.2434, 14.3046] (0.95 confidence) | 14.7236 ns [14.6956, 14.7622] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;dot | 0.687603 ns [0.68587, 0.691345] (0.95 confidence) | 0.681486 ns [0.680242, 0.68257] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;inf&#95;norm | 2.36308 ns [2.36075, 2.36595] (0.95 confidence) | 2.36308 ns [2.36124, 2.36465] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;ldlt | 17.9143 ns [17.859, 17.9821] (0.95 confidence) | 18.1543 ns [18.0874, 18.1918] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;ldlt&#95;solve | 18.0206 ns [17.9586, 18.0755] (0.95 confidence) | 16.6302 ns [16.6122, 16.6538] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;lu | 13.8673 ns [13.8352, 13.8846] (0.95 confidence) | 14.361 ns [14.3391, 14.3813] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;lu&#95;solve | 21.9827 ns [21.9489, 22.0132] (0.95 confidence) | 23.7355 ns [23.6905, 23.7801] (0.95 confidence) | 52.8361 ns [52.7137, 52.9963] (0.95 confidence) | 183.09 ns [182.58, 183.503] (0.95 confidence) |
| d4/la&#95;stack&#95;norm2 | 2.2491 ns [2.24373, 2.25946] (0.95 confidence) | 2.23538 ns [2.22784, 2.25104] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;norm2&#95;scenario&#95;descending | 2.4125 ns [2.40519, 2.41875] (0.95 confidence) | 2.40758 ns [2.40373, 2.41159] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale | 2.40785 ns [2.40256, 2.41215] (0.95 confidence) | 2.39852 ns [2.39489, 2.40408] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;norm2&#95;scenario&#95;sparse | 1.49545 ns [1.4911, 1.49814] (0.95 confidence) | 1.49438 ns [1.49216, 1.49697] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range | 2.16798 ns [2.16297, 2.17033] (0.95 confidence) | 2.1605 ns [2.15458, 2.16447] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;norm2&#95;sq | 0.457462 ns [0.45688, 0.459312] (0.95 confidence) | 0.46516 ns [0.463941, 0.466585] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;solve&#95;from&#95;ldlt | 2.52402 ns [2.51318, 2.53296] (0.95 confidence) | 2.47768 ns [2.47568, 2.47945] (0.95 confidence) | — | — |
| d4/la&#95;stack&#95;solve&#95;from&#95;lu | 3.9732 ns [3.95865, 3.98715] (0.95 confidence) | 4.03234 ns [4.02887, 4.03571] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;det | 41.1898 ns [41.1131, 41.2546] (0.95 confidence) | 28.873 ns [28.8513, 28.9087] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;det&#95;from&#95;ldlt | 0.88214 ns [0.881264, 0.884422] (0.95 confidence) | 0.880188 ns [0.87963, 0.880902] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;det&#95;from&#95;lu | 1.03831 ns [1.03709, 1.04088] (0.95 confidence) | 1.03743 ns [1.03668, 1.03795] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;det&#95;via&#95;lu | 34.0393 ns [33.9885, 34.1047] (0.95 confidence) | 34.5134 ns [34.4911, 34.5356] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;dot | 0.813817 ns [0.80907, 0.821305] (0.95 confidence) | 0.846226 ns [0.845676, 0.846786] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;inf&#95;norm | 3.48393 ns [3.47859, 3.48912] (0.95 confidence) | 3.48404 ns [3.47969, 3.48697] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;ldlt | 25.088 ns [25.0291, 25.1187] (0.95 confidence) | 24.5233 ns [24.5014, 24.5592] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;ldlt&#95;solve | 30.6662 ns [30.595, 30.7694] (0.95 confidence) | 30.2983 ns [30.2532, 30.3627] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;lu | 33.9328 ns [33.7964, 34.0331] (0.95 confidence) | 34.2958 ns [34.2494, 34.3176] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;lu&#95;solve | 48.9417 ns [48.8679, 49.0075] (0.95 confidence) | 47.0332 ns [46.9411, 47.2694] (0.95 confidence) | 74.3705 ns [74.1427, 74.7764] (0.95 confidence) | 250.816 ns [250.442, 251.241] (0.95 confidence) |
| d5/la&#95;stack&#95;norm2 | 2.98015 ns [2.95753, 2.99206] (0.95 confidence) | 2.71271 ns [2.70665, 2.71794] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;norm2&#95;scenario&#95;descending | 3.08399 ns [3.07901, 3.09136] (0.95 confidence) | 3.1851 ns [3.18239, 3.18823] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale | 3.09205 ns [3.08546, 3.09867] (0.95 confidence) | 3.17755 ns [3.17328, 3.18161] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;norm2&#95;scenario&#95;sparse | 1.7401 ns [1.73806, 1.74273] (0.95 confidence) | 1.73229 ns [1.72914, 1.73528] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range | 2.69958 ns [2.68882, 2.70962] (0.95 confidence) | 3.15624 ns [3.14966, 3.16738] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;norm2&#95;sq | 0.52868 ns [0.527629, 0.529971] (0.95 confidence) | 0.527569 ns [0.526732, 0.528412] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;solve&#95;from&#95;ldlt | 3.90215 ns [3.88585, 3.92128] (0.95 confidence) | 3.97755 ns [3.96298, 3.98931] (0.95 confidence) | — | — |
| d5/la&#95;stack&#95;solve&#95;from&#95;lu | 6.16628 ns [6.15322, 6.18557] (0.95 confidence) | 6.09185 ns [6.0862, 6.09642] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;det | 14850.8 ns [14761, 14912] (0.95 confidence) | 14411.1 ns [14376.7, 14425.9] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;det&#95;from&#95;ldlt | 23.5272 ns [23.4886, 23.5627] (0.95 confidence) | 23.1743 ns [23.1304, 23.2651] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;det&#95;from&#95;lu | 23.6266 ns [23.5766, 23.6576] (0.95 confidence) | 23.2611 ns [23.2259, 23.2994] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;det&#95;via&#95;lu | 14301.7 ns [14285.4, 14318.9] (0.95 confidence) | 14232.5 ns [14217.7, 14261.2] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;dot | 11.1068 ns [11.0866, 11.1216] (0.95 confidence) | 11.0042 ns [10.9885, 11.0236] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;inf&#95;norm | 637.909 ns [634.625, 639.696] (0.95 confidence) | 630.117 ns [628.54, 631.393] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;ldlt | 16884.4 ns [16876.5, 16896.4] (0.95 confidence) | 16941.1 ns [16919.4, 16959.7] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;ldlt&#95;solve | 19159.4 ns [19030, 19192.4] (0.95 confidence) | 19186.6 ns [19125.6, 19228.9] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;lu | 15088 ns [15056.4, 15128.1] (0.95 confidence) | 15183.1 ns [15170.6, 15198.5] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;lu&#95;solve | 17206.3 ns [17157.5, 17230.8] (0.95 confidence) | 17505.6 ns [17263.4, 17749.8] (0.95 confidence) | 14495.6 ns [14481, 14519.4] (0.95 confidence) | 12840.2 ns [12828.6, 12881.9] (0.95 confidence) |
| d64/la&#95;stack&#95;norm2 | 59.5681 ns [59.4065, 59.67] (0.95 confidence) | 57.684 ns [57.4383, 57.8988] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;norm2&#95;scenario&#95;descending | 37.6669 ns [37.5408, 37.7362] (0.95 confidence) | 37.0344 ns [36.9153, 37.0848] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale | 37.4383 ns [37.3678, 37.5004] (0.95 confidence) | 36.9123 ns [36.8246, 37.0126] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;norm2&#95;scenario&#95;sparse | 17.4115 ns [17.3765, 17.4555] (0.95 confidence) | 17.1403 ns [17.107, 17.2072] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range | 33.0102 ns [32.9249, 33.0882] (0.95 confidence) | 32.6075 ns [32.5132, 32.6998] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;norm2&#95;sq | 10.9827 ns [10.9652, 11.0071] (0.95 confidence) | 10.9275 ns [10.8959, 10.9604] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;solve&#95;from&#95;ldlt | 1094.97 ns [1092.84, 1096.95] (0.95 confidence) | 1026.15 ns [1025.34, 1027.97] (0.95 confidence) | — | — |
| d64/la&#95;stack&#95;solve&#95;from&#95;lu | 2609.75 ns [2600.43, 2617.82] (0.95 confidence) | 2551.85 ns [2545.6, 2556.24] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;det | 97.8025 ns [97.4755, 98.6699] (0.95 confidence) | 83.6088 ns [83.4681, 83.9718] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;det&#95;from&#95;ldlt | 1.31119 ns [1.30999, 1.31251] (0.95 confidence) | 1.31412 ns [1.31222, 1.31604] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;det&#95;from&#95;ldlt&#95;balanced&#95;range | 8.84089 ns [8.80905, 8.88341] (0.95 confidence) | 9.00993 ns [8.99633, 9.04097] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;det&#95;from&#95;lu | 1.47266 ns [1.47101, 1.47466] (0.95 confidence) | 1.47102 ns [1.46991, 1.47223] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;det&#95;from&#95;lu&#95;balanced&#95;range | 9.02862 ns [9.01063, 9.11903] (0.95 confidence) | 8.97458 ns [8.96075, 9.01066] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;det&#95;via&#95;lu | 90.9294 ns [90.379, 91.9684] (0.95 confidence) | 92.993 ns [92.6753, 93.4541] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;dot | 0.937296 ns [0.936154, 0.939237] (0.95 confidence) | 0.936585 ns [0.932864, 0.940633] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;inf&#95;norm | 8.48274 ns [8.4731, 8.492] (0.95 confidence) | 8.53003 ns [8.51663, 8.53651] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;ldlt | 90.9087 ns [90.7524, 91.0579] (0.95 confidence) | 91.1194 ns [90.3935, 91.9193] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;ldlt&#95;ill&#95;conditioned | 90.6948 ns [90.5823, 90.8167] (0.95 confidence) | 91.9402 ns [90.7683, 92.4929] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;ldlt&#95;solve | 104.126 ns [103.885, 104.376] (0.95 confidence) | 103.621 ns [102.237, 104.2] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;lu | 91.2908 ns [90.9827, 91.8646] (0.95 confidence) | 91.8047 ns [91.3936, 92.2277] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;lu&#95;ill&#95;conditioned | 91.0945 ns [90.7316, 91.7782] (0.95 confidence) | 90.9006 ns [90.5614, 91.0754] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;lu&#95;pivoting | 95.6736 ns [95.2054, 96.1071] (0.95 confidence) | 94.1377 ns [93.9387, 94.4925] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;lu&#95;solve | 133.158 ns [132.483, 133.773] (0.95 confidence) | 132.096 ns [131.87, 132.513] (0.95 confidence) | 166.8 ns [166.498, 167.23] (0.95 confidence) | 342.479 ns [340.695, 343.771] (0.95 confidence) |
| d8/la&#95;stack&#95;norm2 | 4.6549 ns [4.64624, 4.66174] (0.95 confidence) | 4.63565 ns [4.62702, 4.64578] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;norm2&#95;scenario&#95;descending | 4.85115 ns [4.84303, 4.85774] (0.95 confidence) | 4.83837 ns [4.82806, 4.84839] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale | 4.84693 ns [4.84032, 4.85873] (0.95 confidence) | 4.81456 ns [4.80496, 4.82881] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;norm2&#95;scenario&#95;sparse | 2.4801 ns [2.46772, 2.50059] (0.95 confidence) | 2.47345 ns [2.46979, 2.47773] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range | 4.19425 ns [4.17986, 4.2204] (0.95 confidence) | 4.26646 ns [4.26213, 4.27681] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;norm2&#95;sq | 0.660081 ns [0.658018, 0.662154] (0.95 confidence) | 0.663766 ns [0.660676, 0.665466] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;solve&#95;from&#95;ldlt | 8.25904 ns [8.24338, 8.26755] (0.95 confidence) | 8.15204 ns [8.14491, 8.16144] (0.95 confidence) | — | — |
| d8/la&#95;stack&#95;solve&#95;from&#95;lu | 13.4672 ns [13.4421, 13.4873] (0.95 confidence) | 13.4879 ns [13.4684, 13.5178] (0.95 confidence) | — | — |
| exact&#95;d2/det | 0.443679 ns [0.442438, 0.44545] (0.95 confidence) | 0.437289 ns [0.436058, 0.438567] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;direct | 0.499977 ns [0.498816, 0.501041] (0.95 confidence) | 0.440568 ns [0.439745, 0.441307] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;direct&#95;with&#95;errbound | 1.75004 ns [1.74739, 1.75306] (0.95 confidence) | 1.7363 ns [1.734, 1.73904] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;errbound | 1.59409 ns [1.59314, 1.59612] (0.95 confidence) | 1.58445 ns [1.58399, 1.58491] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;exact | 87.6626 ns [87.5841, 87.7443] (0.95 confidence) | 87.8599 ns [87.6289, 88.0253] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;exact&#95;f64&#95;result | 64.8572 ns [64.7395, 64.9237] (0.95 confidence) | 68.1011 ns [67.9669, 68.2989] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;exact&#95;rounded&#95;f64 | 66.3488 ns [66.2122, 66.4581] (0.95 confidence) | 71.2966 ns [71.1855, 71.3897] (0.95 confidence) | — | — |
| exact&#95;d2/det&#95;sign&#95;exact | 2.98316 ns [2.98158, 2.98591] (0.95 confidence) | 2.96916 ns [2.9685, 2.96993] (0.95 confidence) | — | — |
| exact&#95;d2/solve&#95;exact | 2957.2 ns [2952.91, 2961.23] (0.95 confidence) | 2935.15 ns [2931.6, 2938.03] (0.95 confidence) | — | — |
| exact&#95;d2/solve&#95;exact&#95;f64&#95;result | 3139.59 ns [3133.98, 3146.61] (0.95 confidence) | 3086.88 ns [3082.5, 3096.87] (0.95 confidence) | — | — |
| exact&#95;d2/solve&#95;exact&#95;rounded&#95;f64 | 3314.02 ns [3302.92, 3328.75] (0.95 confidence) | 3260.26 ns [3254.14, 3263.61] (0.95 confidence) | — | — |
| exact&#95;d3/det | 0.82411 ns [0.82314, 0.826621] (0.95 confidence) | 0.819392 ns [0.818897, 0.819927] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;direct | 0.830312 ns [0.825973, 0.833761] (0.95 confidence) | 0.819617 ns [0.819367, 0.820099] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;direct&#95;with&#95;errbound | 3.32751 ns [3.32202, 3.33611] (0.95 confidence) | 3.29231 ns [3.28626, 3.29942] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;errbound | 3.31065 ns [3.30839, 3.31656] (0.95 confidence) | 3.29024 ns [3.2876, 3.29195] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;exact | 300.239 ns [300.005, 300.893] (0.95 confidence) | 301.639 ns [301.399, 302.038] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;exact&#95;f64&#95;result | 277.682 ns [277.124, 278.275] (0.95 confidence) | 278.305 ns [277.225, 278.743] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;exact&#95;rounded&#95;f64 | 280.825 ns [280.262, 281.421] (0.95 confidence) | 280.392 ns [279.771, 280.858] (0.95 confidence) | — | — |
| exact&#95;d3/det&#95;sign&#95;exact | 4.73185 ns [4.7286, 4.73725] (0.95 confidence) | 4.69999 ns [4.69853, 4.70147] (0.95 confidence) | — | — |
| exact&#95;d3/solve&#95;exact | 8959.5 ns [8925.77, 8983.71] (0.95 confidence) | 9177.06 ns [9155.27, 9191] (0.95 confidence) | — | — |
| exact&#95;d3/solve&#95;exact&#95;f64&#95;result | 9073.11 ns [9052.45, 9107.74] (0.95 confidence) | 9295.48 ns [9263.12, 9311.62] (0.95 confidence) | — | — |
| exact&#95;d3/solve&#95;exact&#95;rounded&#95;f64 | 9385 ns [9357.29, 9423.09] (0.95 confidence) | 9712.92 ns [9666.38, 9786.92] (0.95 confidence) | — | — |
| exact&#95;d4/det | 3.22455 ns [3.22095, 3.22732] (0.95 confidence) | 3.20562 ns [3.2044, 3.2085] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;direct | 2.36258 ns [2.3609, 2.36516] (0.95 confidence) | 2.3476 ns [2.34706, 2.34791] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;direct&#95;with&#95;errbound | 6.80264 ns [6.79118, 6.81258] (0.95 confidence) | 6.76565 ns [6.76018, 6.76784] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;errbound | 6.79631 ns [6.77655, 6.8132] (0.95 confidence) | 6.7272 ns [6.71827, 6.73528] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;exact | 849.78 ns [848.838, 851.669] (0.95 confidence) | 827.717 ns [826.227, 828.905] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;exact&#95;f64&#95;result | 816.878 ns [814.982, 818.393] (0.95 confidence) | 790.724 ns [789.037, 792.132] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;exact&#95;rounded&#95;f64 | 833.945 ns [831.964, 835.694] (0.95 confidence) | 806.199 ns [804.186, 807.928] (0.95 confidence) | — | — |
| exact&#95;d4/det&#95;sign&#95;exact | 8.14947 ns [8.14417, 8.15266] (0.95 confidence) | 7.91046 ns [7.9081, 7.91345] (0.95 confidence) | — | — |
| exact&#95;d4/solve&#95;exact | 18661.7 ns [18620, 18753.5] (0.95 confidence) | 18596.9 ns [18566.3, 18654.6] (0.95 confidence) | — | — |
| exact&#95;d4/solve&#95;exact&#95;f64&#95;result | 18699 ns [18677.6, 18720.1] (0.95 confidence) | 18696.9 ns [18664.6, 18777.3] (0.95 confidence) | — | — |
| exact&#95;d4/solve&#95;exact&#95;rounded&#95;f64 | 19279.5 ns [19198.3, 19356.3] (0.95 confidence) | 19374.3 ns [19355.1, 19411.8] (0.95 confidence) | — | — |
| exact&#95;d5/det | 27.7042 ns [27.5857, 27.781] (0.95 confidence) | 28.3684 ns [28.3093, 28.4114] (0.95 confidence) | — | — |
| exact&#95;d5/det&#95;exact | 2928.1 ns [2920.03, 2934] (0.95 confidence) | 2863.85 ns [2858.06, 2879.64] (0.95 confidence) | — | — |
| exact&#95;d5/det&#95;exact&#95;f64&#95;result | 2884.1 ns [2877.28, 2890.81] (0.95 confidence) | 2822.07 ns [2817.91, 2827.32] (0.95 confidence) | — | — |
| exact&#95;d5/det&#95;exact&#95;rounded&#95;f64 | 2919.86 ns [2914.71, 2929.87] (0.95 confidence) | 2830.34 ns [2823.56, 2835.71] (0.95 confidence) | — | — |
| exact&#95;d5/det&#95;sign&#95;exact | 2904.44 ns [2900.51, 2911.67] (0.95 confidence) | 2806.7 ns [2804.1, 2815.87] (0.95 confidence) | — | — |
| exact&#95;d5/solve&#95;exact | 30279.5 ns [30214.6, 30367.5] (0.95 confidence) | 30478 ns [30421.1, 30559.8] (0.95 confidence) | — | — |
| exact&#95;d5/solve&#95;exact&#95;f64&#95;result | 30388.5 ns [30294.2, 30505.7] (0.95 confidence) | 30612.7 ns [30558.9, 30679] (0.95 confidence) | — | — |
| exact&#95;d5/solve&#95;exact&#95;rounded&#95;f64 | 31362.9 ns [31261.6, 31574.8] (0.95 confidence) | 31279.5 ns [31203.4, 31354.2] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/det&#95;exact | 946.209 ns [944.792, 948.64] (0.95 confidence) | 916.54 ns [913.355, 919.521] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/det&#95;sign&#95;exact | 8.12176 ns [8.11783, 8.12548] (0.95 confidence) | 7.91046 ns [7.90934, 7.91133] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/solve&#95;exact | 15283.7 ns [15237.5, 15337.2] (0.95 confidence) | 15327.5 ns [15264.3, 15362.7] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/solve&#95;exact&#95;f64&#95;result | 15319.6 ns [15248.4, 15348.9] (0.95 confidence) | 15546.8 ns [15504.4, 15617.1] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;4x4/solve&#95;exact&#95;rounded&#95;f64 | 15762.2 ns [15736.1, 15785.8] (0.95 confidence) | 15848.6 ns [15810, 15882.4] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/det&#95;exact | 2948.96 ns [2941.39, 2955.37] (0.95 confidence) | 2955.89 ns [2949.82, 2962.57] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/det&#95;sign&#95;exact | 2895.85 ns [2890, 2903.04] (0.95 confidence) | 2842.84 ns [2832.41, 2851.67] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/solve&#95;exact | 24667.6 ns [24618.5, 24701.6] (0.95 confidence) | 24781.4 ns [24711.2, 24837.2] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/solve&#95;exact&#95;f64&#95;result | 24757.8 ns [24693.1, 24802.9] (0.95 confidence) | 25136.8 ns [25109.1, 25171.2] (0.95 confidence) | — | — |
| exact&#95;hilbert&#95;5x5/solve&#95;exact&#95;rounded&#95;f64 | 25678.4 ns [25604.5, 25714.3] (0.95 confidence) | 25770.4 ns [25672.6, 25858.1] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/det&#95;exact | 317.618 ns [316.346, 318.505] (0.95 confidence) | 313.867 ns [312.205, 316.276] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/det&#95;sign&#95;exact | 290.197 ns [289.848, 290.943] (0.95 confidence) | 282.423 ns [281.867, 283.086] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/solve&#95;exact | 9033.56 ns [9018.39, 9065.29] (0.95 confidence) | 8879.99 ns [8870.17, 8892.12] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/solve&#95;exact&#95;f64&#95;result | 9174.28 ns [9164.49, 9185.1] (0.95 confidence) | 9003.24 ns [8987.3, 9021.08] (0.95 confidence) | — | — |
| exact&#95;large&#95;entries&#95;3x3/solve&#95;exact&#95;rounded&#95;f64 | 9484.44 ns [9475.99, 9492.04] (0.95 confidence) | 9384.71 ns [9367.01, 9403.03] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/det&#95;exact | 183.201 ns [183.079, 183.49] (0.95 confidence) | 187.488 ns [187.362, 187.666] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/det&#95;sign&#95;exact | 163.513 ns [163.18, 163.718] (0.95 confidence) | 163.947 ns [163.792, 164.104] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/solve&#95;exact | 974.857 ns [973.16, 976.296] (0.95 confidence) | 964.384 ns [962.376, 968.01] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/solve&#95;exact&#95;f64&#95;result | 998.111 ns [995.929, 999.608] (0.95 confidence) | 1003.88 ns [1001.89, 1005.62] (0.95 confidence) | — | — |
| exact&#95;near&#95;singular&#95;3x3/solve&#95;exact&#95;rounded&#95;f64 | 986.581 ns [984.651, 989.052] (0.95 confidence) | 980.728 ns [979.225, 982.815] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/det&#95;exact | 3485.12 ns [3477.31, 3490.73] (0.95 confidence) | 3418.69 ns [3412.47, 3430.15] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/det&#95;sign&#95;exact | 151.878 ns [151.205, 152.341] (0.95 confidence) | 145.328 ns [144.915, 145.743] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/solve&#95;exact | 37248.2 ns [37170.3, 37369.5] (0.95 confidence) | 35075.2 ns [34992, 35245.7] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/solve&#95;exact&#95;f64&#95;result | 38038.7 ns [37990.5, 38125.9] (0.95 confidence) | 36092 ns [36063.4, 36135.9] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d2/solve&#95;exact&#95;rounded&#95;f64 | 38410.2 ns [38275.2, 38516.3] (0.95 confidence) | 35985.3 ns [35861.6, 36110.7] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/det&#95;exact | 9664.84 ns [9655.49, 9678.2] (0.95 confidence) | 9367.29 ns [9339.14, 9383.44] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/det&#95;sign&#95;exact | 231.519 ns [231.41, 231.695] (0.95 confidence) | 224.71 ns [224.547, 224.737] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/solve&#95;exact | 84925.6 ns [84663.8, 85097.6] (0.95 confidence) | 79115.9 ns [78919.3, 79509.5] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/solve&#95;exact&#95;f64&#95;result | 85586.4 ns [85480.1, 85774.2] (0.95 confidence) | 79919.8 ns [79654.1, 80134] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d3/solve&#95;exact&#95;rounded&#95;f64 | 86112.3 ns [86012.8, 86237.1] (0.95 confidence) | 80599.1 ns [80461.2, 80725.6] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/det&#95;exact | 23734 ns [23681, 23761.9] (0.95 confidence) | 22670.9 ns [22616.9, 22768.4] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/det&#95;sign&#95;exact | 464.103 ns [462.676, 465.864] (0.95 confidence) | 439.395 ns [439.081, 439.748] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/solve&#95;exact | 152306 ns [152014, 152463] (0.95 confidence) | 142242 ns [141970, 142520] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/solve&#95;exact&#95;f64&#95;result | 154023 ns [153682, 154344] (0.95 confidence) | 141811 ns [141586, 142314] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d4/solve&#95;exact&#95;rounded&#95;f64 | 154793 ns [154632, 154971] (0.95 confidence) | 143282 ns [143036, 143384] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/det&#95;exact | 58486.6 ns [58294.5, 58639.9] (0.95 confidence) | 59084.4 ns [59009.3, 59191.1] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/det&#95;sign&#95;exact | 57291.4 ns [57100.4, 57521.3] (0.95 confidence) | 57761 ns [57619.6, 57869.5] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/solve&#95;exact | 255271 ns [254756, 256664] (0.95 confidence) | 238297 ns [237811, 238681] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/solve&#95;exact&#95;f64&#95;result | 257775 ns [256632, 258987] (0.95 confidence) | 240906 ns [240594, 241328] (0.95 confidence) | — | — |
| exact&#95;random&#95;corpus&#95;d5/solve&#95;exact&#95;rounded&#95;f64 | 264550 ns [263746, 264704] (0.95 confidence) | 240700 ns [240202, 241376] (0.95 confidence) | — | — |
| gram/1x2/dependent/la&#95;stack | 0.688217 ns [0.68765, 0.688682] (0.95 confidence) | 0.694062 ns [0.693371, 0.695676] (0.95 confidence) | — | — |
| gram/1x2/mixed&#95;scale/la&#95;stack | 0.690096 ns [0.689019, 0.690537] (0.95 confidence) | 0.694036 ns [0.692216, 0.697763] (0.95 confidence) | — | — |
| gram/1x2/near&#95;dependent/la&#95;stack | 0.688782 ns [0.688077, 0.689137] (0.95 confidence) | 0.697215 ns [0.694952, 0.698906] (0.95 confidence) | — | — |
| gram/1x2/orthogonal/la&#95;stack | 0.688122 ns [0.687837, 0.688696] (0.95 confidence) | 0.693527 ns [0.692235, 0.69475] (0.95 confidence) | — | — |
| gram/2x2/dependent/la&#95;stack | 1.3657 ns [1.36206, 1.36886] (0.95 confidence) | 1.32851 ns [1.32154, 1.33759] (0.95 confidence) | — | — |
| gram/2x2/mixed&#95;scale/la&#95;stack | 1.34548 ns [1.3416, 1.34823] (0.95 confidence) | 1.33173 ns [1.32113, 1.33881] (0.95 confidence) | — | — |
| gram/2x2/near&#95;dependent/la&#95;stack | 1.36164 ns [1.35767, 1.3686] (0.95 confidence) | 1.32964 ns [1.32014, 1.33824] (0.95 confidence) | — | — |
| gram/2x2/orthogonal/la&#95;stack | 1.34751 ns [1.34427, 1.35038] (0.95 confidence) | 1.33146 ns [1.31992, 1.3425] (0.95 confidence) | — | — |
| gram/2x3/dependent/la&#95;stack | 1.32146 ns [1.31889, 1.32368] (0.95 confidence) | 1.26532 ns [1.26345, 1.26701] (0.95 confidence) | — | — |
| gram/2x3/mixed&#95;scale/la&#95;stack | 1.32217 ns [1.32043, 1.32399] (0.95 confidence) | 1.30569 ns [1.28722, 1.31013] (0.95 confidence) | — | — |
| gram/2x3/near&#95;dependent/la&#95;stack | 1.31985 ns [1.31813, 1.32111] (0.95 confidence) | 1.26324 ns [1.26135, 1.26425] (0.95 confidence) | — | — |
| gram/2x3/orthogonal/la&#95;stack | 1.31788 ns [1.31681, 1.32096] (0.95 confidence) | 1.26695 ns [1.26474, 1.27078] (0.95 confidence) | — | — |
| gram/3x3/dependent/la&#95;stack | 3.04131 ns [3.03792, 3.04537] (0.95 confidence) | 3.00348 ns [3.00095, 3.00685] (0.95 confidence) | — | — |
| gram/3x3/mixed&#95;scale/la&#95;stack | 3.03873 ns [3.03611, 3.04196] (0.95 confidence) | 3.0055 ns [3.00316, 3.00788] (0.95 confidence) | — | — |
| gram/3x3/near&#95;dependent/la&#95;stack | 3.03986 ns [3.03721, 3.04313] (0.95 confidence) | 3.00495 ns [3.00312, 3.00938] (0.95 confidence) | — | — |
| gram/3x3/orthogonal/la&#95;stack | 3.03926 ns [3.03228, 3.04449] (0.95 confidence) | 3.00411 ns [3.00057, 3.00623] (0.95 confidence) | — | — |
| gram/3x4/dependent/la&#95;stack | 3.07849 ns [3.07602, 3.07997] (0.95 confidence) | 3.0487 ns [3.0467, 3.05091] (0.95 confidence) | — | — |
| gram/3x4/mixed&#95;scale/la&#95;stack | 3.07783 ns [3.07502, 3.08121] (0.95 confidence) | 3.05111 ns [3.04656, 3.05354] (0.95 confidence) | — | — |
| gram/3x4/near&#95;dependent/la&#95;stack | 3.07785 ns [3.07631, 3.0792] (0.95 confidence) | 3.04756 ns [3.04539, 3.0492] (0.95 confidence) | — | — |
| gram/3x4/orthogonal/la&#95;stack | 3.08027 ns [3.07837, 3.08252] (0.95 confidence) | 3.04854 ns [3.04619, 3.05132] (0.95 confidence) | — | — |
| gram/4x4/dependent/la&#95;stack | 7.67728 ns [7.66077, 7.68784] (0.95 confidence) | 7.51104 ns [7.44533, 7.81867] (0.95 confidence) | — | — |
| gram/4x4/mixed&#95;scale/la&#95;stack | 7.6665 ns [7.64964, 7.68983] (0.95 confidence) | 7.49418 ns [7.45441, 7.79816] (0.95 confidence) | — | — |
| gram/4x4/near&#95;dependent/la&#95;stack | 7.66139 ns [7.64998, 7.67863] (0.95 confidence) | 7.49123 ns [7.45499, 7.81016] (0.95 confidence) | — | — |
| gram/4x4/orthogonal/la&#95;stack | 7.69578 ns [7.67735, 7.70371] (0.95 confidence) | 7.50641 ns [7.47916, 7.83406] (0.95 confidence) | — | — |
| gram/4x5/dependent/la&#95;stack | 9.65489 ns [9.59153, 9.71598] (0.95 confidence) | 9.44447 ns [9.41627, 9.46087] (0.95 confidence) | — | — |
| gram/4x5/mixed&#95;scale/la&#95;stack | 9.51485 ns [9.50711, 9.52545] (0.95 confidence) | 9.49321 ns [9.48481, 9.49742] (0.95 confidence) | — | — |
| gram/4x5/near&#95;dependent/la&#95;stack | 9.50886 ns [9.49807, 9.51946] (0.95 confidence) | 9.48421 ns [9.46483, 9.49632] (0.95 confidence) | — | — |
| gram/4x5/orthogonal/la&#95;stack | 9.53377 ns [9.52453, 9.54297] (0.95 confidence) | 9.4848 ns [9.47905, 9.49089] (0.95 confidence) | — | — |
| gram/5x5/dependent/la&#95;stack | 18.2206 ns [18.1646, 18.2857] (0.95 confidence) | 18.4022 ns [18.3876, 18.4107] (0.95 confidence) | — | — |
| gram/5x5/mixed&#95;scale/la&#95;stack | 18.3668 ns [18.3012, 18.4162] (0.95 confidence) | 18.3363 ns [18.3172, 18.3758] (0.95 confidence) | — | — |
| gram/5x5/near&#95;dependent/la&#95;stack | 18.0665 ns [18.0163, 18.1384] (0.95 confidence) | 18.2766 ns [18.2582, 18.342] (0.95 confidence) | — | — |
| gram/5x5/orthogonal/la&#95;stack | 18.1869 ns [18.139, 18.248] (0.95 confidence) | 18.3129 ns [18.2477, 18.3869] (0.95 confidence) | — | — |
| gram/5x6/dependent/la&#95;stack | 19.2391 ns [19.1272, 19.288] (0.95 confidence) | 19.2217 ns [19.1547, 19.262] (0.95 confidence) | — | — |
| gram/5x6/mixed&#95;scale/la&#95;stack | 19.2749 ns [19.2525, 19.3233] (0.95 confidence) | 19.0876 ns [19.0441, 19.1213] (0.95 confidence) | — | — |
| gram/5x6/near&#95;dependent/la&#95;stack | 19.2643 ns [19.206, 19.3404] (0.95 confidence) | 19.122 ns [19.0545, 19.1658] (0.95 confidence) | — | — |
| gram/5x6/orthogonal/la&#95;stack | 19.2181 ns [19.1769, 19.2385] (0.95 confidence) | 19.0179 ns [18.9349, 19.0722] (0.95 confidence) | — | — |
| gram/6x6/dependent/la&#95;stack | 29.669 ns [29.457, 29.7913] (0.95 confidence) | 29.0209 ns [28.8655, 29.7152] (0.95 confidence) | — | — |
| gram/6x6/mixed&#95;scale/la&#95;stack | 29.7003 ns [29.6398, 29.8245] (0.95 confidence) | 28.9178 ns [28.8143, 29.0101] (0.95 confidence) | — | — |
| gram/6x6/near&#95;dependent/la&#95;stack | 29.6557 ns [29.5101, 29.7169] (0.95 confidence) | 29.6721 ns [28.897, 30.2522] (0.95 confidence) | — | — |
| gram/6x6/orthogonal/la&#95;stack | 29.9119 ns [29.8787, 29.9658] (0.95 confidence) | 28.941 ns [28.8429, 29.5409] (0.95 confidence) | — | — |
| gram/6x7/dependent/la&#95;stack | 32.6161 ns [32.5437, 32.7134] (0.95 confidence) | 31.81 ns [31.7472, 31.8796] (0.95 confidence) | — | — |
| gram/6x7/mixed&#95;scale/la&#95;stack | 32.5786 ns [32.5297, 32.6685] (0.95 confidence) | 31.894 ns [31.8708, 32.0718] (0.95 confidence) | — | — |
| gram/6x7/near&#95;dependent/la&#95;stack | 32.4376 ns [32.3558, 32.4729] (0.95 confidence) | 31.9548 ns [31.8528, 32.3402] (0.95 confidence) | — | — |
| gram/6x7/orthogonal/la&#95;stack | 32.4649 ns [32.3892, 32.5104] (0.95 confidence) | 32.0095 ns [31.6173, 32.3394] (0.95 confidence) | — | — |
| gram/7x7/dependent/la&#95;stack | 42.0611 ns [41.8819, 42.2704] (0.95 confidence) | 38.6924 ns [38.6605, 38.7602] (0.95 confidence) | — | — |
| gram/7x7/mixed&#95;scale/la&#95;stack | 41.3303 ns [41.2189, 41.6089] (0.95 confidence) | 38.8711 ns [38.8203, 39.0061] (0.95 confidence) | — | — |
| gram/7x7/near&#95;dependent/la&#95;stack | 41.9281 ns [41.8665, 42.0856] (0.95 confidence) | 39.6471 ns [39.2629, 40.1935] (0.95 confidence) | — | — |
| gram/7x7/orthogonal/la&#95;stack | 44.3947 ns [43.7069, 44.6195] (0.95 confidence) | 38.7539 ns [38.689, 38.8093] (0.95 confidence) | — | — |
| gram/7x8/dependent/la&#95;stack | 42.2354 ns [41.9852, 42.529] (0.95 confidence) | 40.0351 ns [40.006, 40.112] (0.95 confidence) | — | — |
| gram/7x8/mixed&#95;scale/la&#95;stack | 42.6812 ns [42.3636, 42.8276] (0.95 confidence) | 40.0253 ns [39.9925, 40.0585] (0.95 confidence) | — | — |
| gram/7x8/near&#95;dependent/la&#95;stack | 42.2239 ns [41.8721, 42.6539] (0.95 confidence) | 40.1305 ns [40.0172, 40.1773] (0.95 confidence) | — | — |
| gram/7x8/orthogonal/la&#95;stack | 42.7217 ns [42.0732, 43.0329] (0.95 confidence) | 39.9929 ns [39.9307, 40.0751] (0.95 confidence) | — | — |
| gram/8x8/dependent/la&#95;stack | 52.1925 ns [52.1166, 52.2838] (0.95 confidence) | 51.7686 ns [51.671, 51.8677] (0.95 confidence) | — | — |
| gram/8x8/mixed&#95;scale/la&#95;stack | 52.55 ns [52.3887, 52.9233] (0.95 confidence) | 51.937 ns [51.8631, 52.0236] (0.95 confidence) | — | — |
| gram/8x8/near&#95;dependent/la&#95;stack | 52.045 ns [52.011, 52.0954] (0.95 confidence) | 51.9214 ns [51.8504, 52.0065] (0.95 confidence) | — | — |
| gram/8x8/orthogonal/la&#95;stack | 53.771 ns [53.4799, 54.0122] (0.95 confidence) | 51.7652 ns [51.7001, 51.8108] (0.95 confidence) | — | — |
| interval&#95;d3/dense&#95;inconclusive | 43.6777 ns [43.4999, 43.893] (0.95 confidence) | 41.3044 ns [41.2669, 41.3169] (0.95 confidence) | — | — |
| interval&#95;d3/dense&#95;sign | 40.3511 ns [40.0195, 40.6957] (0.95 confidence) | 38.6382 ns [38.2364, 39.0316] (0.95 confidence) | — | — |
| interval&#95;d3/lifted&#95;assembly | 8.36435 ns [8.3608, 8.37158] (0.95 confidence) | 8.34372 ns [8.33228, 8.3518] (0.95 confidence) | — | — |
| interval&#95;d3/lifted&#95;assembly&#95;and&#95;sign | 37.2083 ns [37.0373, 37.4007] (0.95 confidence) | 36.7324 ns [36.5291, 37.0782] (0.95 confidence) | — | — |
| interval&#95;d3/sparse&#95;sign | 32.3931 ns [31.9353, 32.7281] (0.95 confidence) | 31.4919 ns [31.2831, 31.8013] (0.95 confidence) | — | — |
| interval&#95;d4/dense&#95;inconclusive | 108.142 ns [106.647, 109.601] (0.95 confidence) | 102.425 ns [101.754, 102.813] (0.95 confidence) | — | — |
| interval&#95;d4/dense&#95;sign | 103.766 ns [102.826, 104.437] (0.95 confidence) | 99.9475 ns [99.0867, 100.445] (0.95 confidence) | — | — |
| interval&#95;d4/lifted&#95;assembly | 13.1186 ns [13.107, 13.1326] (0.95 confidence) | 13.2665 ns [13.254, 13.2767] (0.95 confidence) | — | — |
| interval&#95;d4/lifted&#95;assembly&#95;and&#95;sign | 85.012 ns [84.4391, 85.69] (0.95 confidence) | 81.851 ns [81.674, 82.0605] (0.95 confidence) | — | — |
| interval&#95;d4/sparse&#95;sign | 70.4079 ns [69.9608, 71.0626] (0.95 confidence) | 68.4328 ns [68.112, 68.5465] (0.95 confidence) | — | — |
| interval&#95;d5/dense&#95;inconclusive | 245.717 ns [244.554, 246.779] (0.95 confidence) | 236.734 ns [235.534, 237.402] (0.95 confidence) | — | — |
| interval&#95;d5/dense&#95;sign | 233.052 ns [231.256, 236.012] (0.95 confidence) | 224.167 ns [222.97, 225.236] (0.95 confidence) | — | — |
| interval&#95;d5/lifted&#95;assembly | 21.6333 ns [21.6145, 21.6442] (0.95 confidence) | 21.6985 ns [21.5522, 21.8389] (0.95 confidence) | — | — |
| interval&#95;d5/lifted&#95;assembly&#95;and&#95;sign | 196.702 ns [194.099, 198.18] (0.95 confidence) | 194.712 ns [192.118, 197.618] (0.95 confidence) | — | — |
| interval&#95;d5/sparse&#95;sign | 162.889 ns [160.753, 165.113] (0.95 confidence) | 152.433 ns [152.079, 152.681] (0.95 confidence) | — | — |
| interval&#95;d6/dense&#95;inconclusive | 564.455 ns [562.876, 567.92] (0.95 confidence) | 543.583 ns [541.755, 545.37] (0.95 confidence) | — | — |
| interval&#95;d6/dense&#95;sign | 559.805 ns [558.35, 564.537] (0.95 confidence) | 540.338 ns [537.608, 541.253] (0.95 confidence) | — | — |
| interval&#95;d6/lifted&#95;assembly | 29.825 ns [29.4025, 30.6451] (0.95 confidence) | 30.7535 ns [30.6853, 30.8279] (0.95 confidence) | — | — |
| interval&#95;d6/lifted&#95;assembly&#95;and&#95;sign | 407.166 ns [405.455, 408.686] (0.95 confidence) | 387.086 ns [386.071, 387.815] (0.95 confidence) | — | — |
| interval&#95;d6/sparse&#95;sign | 378.308 ns [377.261, 380.151] (0.95 confidence) | 361.278 ns [359.893, 362.105] (0.95 confidence) | — | — |
| interval&#95;d7/dense&#95;inconclusive | 1117.33 ns [1112.95, 1125.97] (0.95 confidence) | 1107.81 ns [1099.47, 1111.9] (0.95 confidence) | — | — |
| interval&#95;d7/dense&#95;sign | 1101.14 ns [1094.18, 1108.77] (0.95 confidence) | 1067.55 ns [1060.51, 1072.87] (0.95 confidence) | — | — |
| interval&#95;d7/lifted&#95;assembly | 38.1579 ns [38.1159, 38.1925] (0.95 confidence) | 37.3233 ns [37.2833, 37.3812] (0.95 confidence) | — | — |
| interval&#95;d7/lifted&#95;assembly&#95;and&#95;sign | 832.766 ns [829.383, 839.704] (0.95 confidence) | 792.012 ns [790.496, 792.97] (0.95 confidence) | — | — |
| interval&#95;d7/sparse&#95;sign | 788.877 ns [785.238, 790.196] (0.95 confidence) | 777.791 ns [774.317, 779.602] (0.95 confidence) | — | — |
| interval&#95;det&#95;sign/d4&#95;conclusive&#95;lifted | 71.2195 ns [70.4772, 72.2118] (0.95 confidence) | 68.8496 ns [68.4276, 69.4791] (0.95 confidence) | — | — |
| interval&#95;det&#95;sign/d4&#95;inconclusive&#95;lifted | 68.7962 ns [68.4381, 69.2044] (0.95 confidence) | 65.5913 ns [65.3951, 65.9149] (0.95 confidence) | — | — |
| interval&#95;det&#95;sign/d7&#95;conclusive&#95;lifted | 790.309 ns [786.742, 793.3] (0.95 confidence) | 759.483 ns [757.97, 760.993] (0.95 confidence) | — | — |
| interval&#95;division/negative | — | 5.20122 ns [5.19661, 5.20428] (0.95 confidence) | — | — |
| interval&#95;division/positive | — | 5.21446 ns [5.21009, 5.223] (0.95 confidence) | — | — |
| interval&#95;division/reciprocal | — | 10.0225 ns [9.99006, 10.0526] (0.95 confidence) | — | — |
| interval&#95;division/subnormal | — | 3.83682 ns [3.82606, 3.88788] (0.95 confidence) | — | — |
| interval&#95;division/tiny&#95;ratio | — | 6.11068 ns [6.10368, 6.11761] (0.95 confidence) | — | — |
| interval&#95;scalar/cancellation/add | 2.15436 ns [2.15334, 2.15514] (0.95 confidence) | 2.11986 ns [2.11704, 2.12293] (0.95 confidence) | — | — |
| interval&#95;scalar/cancellation/multiply | 1.87897 ns [1.87213, 1.88336] (0.95 confidence) | 1.87515 ns [1.86626, 1.88136] (0.95 confidence) | — | — |
| interval&#95;scalar/cancellation/square | 1.90296 ns [1.90161, 1.90421] (0.95 confidence) | 1.86607 ns [1.86356, 1.86788] (0.95 confidence) | — | — |
| interval&#95;scalar/exact/add | 2.15796 ns [2.1553, 2.16026] (0.95 confidence) | 2.13202 ns [2.12893, 2.13555] (0.95 confidence) | — | — |
| interval&#95;scalar/exact/multiply | 2.01307 ns [2.00582, 2.02236] (0.95 confidence) | 2.09247 ns [2.08459, 2.11716] (0.95 confidence) | — | — |
| interval&#95;scalar/exact/square | 1.90559 ns [1.9038, 1.90737] (0.95 confidence) | 1.86367 ns [1.86075, 1.86548] (0.95 confidence) | — | — |
| interval&#95;scalar/point/add | 2.18918 ns [2.18518, 2.19274] (0.95 confidence) | 2.17817 ns [2.1762, 2.18302] (0.95 confidence) | — | — |
| interval&#95;scalar/point/multiply | 2.66183 ns [2.64974, 2.68212] (0.95 confidence) | 2.5785 ns [2.57152, 2.58862] (0.95 confidence) | — | — |
| interval&#95;scalar/point/square | 2.14533 ns [2.14202, 2.14614] (0.95 confidence) | 2.12382 ns [2.11939, 2.12789] (0.95 confidence) | — | — |
| interval&#95;scalar/subnormal/add | 2.42644 ns [2.42473, 2.42769] (0.95 confidence) | 2.40163 ns [2.4001, 2.40436] (0.95 confidence) | — | — |
| interval&#95;scalar/subnormal/multiply | 3.02735 ns [3.01317, 3.03925] (0.95 confidence) | 2.84418 ns [2.83234, 2.85112] (0.95 confidence) | — | — |
| interval&#95;scalar/subnormal/square | 2.39543 ns [2.39453, 2.39653] (0.95 confidence) | 2.37575 ns [2.37389, 2.37671] (0.95 confidence) | — | — |
| interval&#95;scalar/zero/add | 0.948981 ns [0.947427, 0.950704] (0.95 confidence) | 0.939264 ns [0.937517, 0.940362] (0.95 confidence) | — | — |
| interval&#95;scalar/zero/multiply | 1.79168 ns [1.78894, 1.7953] (0.95 confidence) | 1.76617 ns [1.76121, 1.7691] (0.95 confidence) | — | — |
| interval&#95;scalar/zero/square | 0.892078 ns [0.890487, 0.893158] (0.95 confidence) | 0.88179 ns [0.880943, 0.882104] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d1/dense | — | 228.171 ns [227.935, 228.449] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d1/near&#95;singular | — | 290.58 ns [290.41, 290.79] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d1/range&#95;exhausted | — | 154.653 ns [154.561, 154.917] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d1/singular | — | 51.0022 ns [50.9241, 51.0911] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d1/uncertain&#95;dense | — | 373.058 ns [372.495, 373.413] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d2/dense | — | 2898.98 ns [2892.42, 2904.75] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d2/near&#95;singular | — | 1593.61 ns [1588.65, 1596.19] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d2/range&#95;exhausted | — | 253.57 ns [253.481, 253.618] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d2/singular | — | 167.945 ns [167.718, 168.072] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d2/uncertain&#95;dense | — | 3993.07 ns [3984.14, 4003.14] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d3/dense | — | 7436.99 ns [7426.99, 7442.29] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d3/near&#95;singular | — | 2021.37 ns [2019.02, 2025.14] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d3/range&#95;exhausted | — | 411.269 ns [410.388, 411.702] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d3/singular | — | 322.511 ns [321.863, 322.696] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d3/uncertain&#95;dense | — | 9084.15 ns [9063.08, 9095.07] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d4/dense | — | 14537.9 ns [14522.6, 14558.8] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d4/near&#95;singular | — | 2592.62 ns [2586.02, 2595.86] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d4/range&#95;exhausted | — | 757.256 ns [756.306, 757.753] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d4/singular | — | 577.434 ns [576.907, 577.611] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d4/uncertain&#95;dense | — | 17160.7 ns [17125.4, 17210.8] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d5/dense | — | 24971 ns [24921.6, 25011.2] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d5/near&#95;singular | — | 3441.72 ns [3434.4, 3445.23] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d5/range&#95;exhausted | — | 1248.86 ns [1247.82, 1252.04] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d5/singular | — | 1015.14 ns [1013.73, 1017.47] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d5/uncertain&#95;dense | — | 28461.1 ns [28425.6, 28511.7] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d6/dense | — | 39282.1 ns [39234.9, 39326.4] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d6/near&#95;singular | — | 4323.4 ns [4312.04, 4328.06] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d6/range&#95;exhausted | — | 1577.63 ns [1576.28, 1579.71] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d6/singular | — | 1312.98 ns [1311.06, 1315.71] (0.95 confidence) | — | — |
| interval&#95;solve&#95;d6/uncertain&#95;dense | — | 44586.3 ns [44340.3, 44782.6] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;batch/origin/construction | 1.95764 ns [1.95212, 1.96343] (0.95 confidence) | 2.05561 ns [2.05276, 2.05955] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;batch/origin/prepared | 12.7746 ns [12.7689, 12.7877] (0.95 confidence) | 13.8189 ns [13.812, 13.8238] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;batch/origin/with&#95;construction | 14.3201 ns [14.3135, 14.3287] (0.95 confidence) | 15.3389 ns [15.3355, 15.3437] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;batch/translated/construction | 1.96983 ns [1.96512, 1.97431] (0.95 confidence) | 2.05868 ns [2.05592, 2.06048] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;batch/translated/prepared | 42.0752 ns [41.9876, 42.1221] (0.95 confidence) | 40.3291 ns [40.2074, 40.4052] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;batch/translated/with&#95;construction | 43.7079 ns [43.6371, 43.742] (0.95 confidence) | 42.1336 ns [42.0777, 42.1853] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/cancellation/difference | 16.5131 ns [16.4885, 16.5354] (0.95 confidence) | 10.2895 ns [10.2876, 10.2917] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/cancellation/difference&#95;endpoints | 17.5547 ns [17.5259, 17.5748] (0.95 confidence) | 11.4428 ns [11.413, 11.497] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/cancellation/dot | 5.46448 ns [5.46055, 5.46947] (0.95 confidence) | 6.00467 ns [6.00062, 6.01007] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/cancellation/dot&#95;endpoints | 6.49568 ns [6.48663, 6.50723] (0.95 confidence) | 7.03352 ns [7.02085, 7.04773] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/dense/difference | 15.997 ns [15.9744, 16.0233] (0.95 confidence) | 10.2916 ns [10.2795, 10.311] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/dense/difference&#95;endpoints | 19.603 ns [19.5704, 19.6694] (0.95 confidence) | 12.9642 ns [12.877, 13.1434] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/dense/dot | 5.17149 ns [5.16783, 5.17996] (0.95 confidence) | 5.52486 ns [5.51593, 5.52924] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/dense/dot&#95;endpoints | 5.9708 ns [5.96235, 5.97677] (0.95 confidence) | 6.4219 ns [6.41058, 6.43278] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/dense&#95;inexact/difference | 16.0284 ns [15.9993, 16.0517] (0.95 confidence) | 10.3453 ns [10.3313, 10.3542] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/dense&#95;inexact/difference&#95;endpoints | 17.203 ns [17.1726, 17.2344] (0.95 confidence) | 11.5247 ns [11.4952, 11.5494] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/dense&#95;inexact/dot | 5.16956 ns [5.16381, 5.18197] (0.95 confidence) | 5.52608 ns [5.50131, 5.53086] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/dense&#95;inexact/dot&#95;endpoints | 5.96391 ns [5.9589, 5.97112] (0.95 confidence) | 6.49137 ns [6.48339, 6.50486] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/overflow/dot | 1.24268 ns [1.24189, 1.24351] (0.95 confidence) | 1.07081 ns [1.06955, 1.07171] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/sparse/difference | 6.48799 ns [6.48463, 6.49261] (0.95 confidence) | 4.73494 ns [4.732, 4.73814] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/sparse/difference&#95;endpoints | 7.28506 ns [7.2322, 7.30063] (0.95 confidence) | 5.53558 ns [5.5228, 5.55072] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/sparse/dot | 2.23767 ns [2.23539, 2.23918] (0.95 confidence) | 2.22623 ns [2.22483, 2.22631] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/sparse/dot&#95;endpoints | 1.98997 ns [1.98871, 1.9906] (0.95 confidence) | 2.22627 ns [2.22499, 2.22638] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/underflow/difference | 5.8908 ns [5.88564, 5.89546] (0.95 confidence) | 4.45626 ns [4.45491, 4.45787] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/underflow/difference&#95;endpoints | 5.95316 ns [5.94987, 5.95744] (0.95 confidence) | 4.60089 ns [4.59955, 4.60306] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/underflow/dot | 1.49232 ns [1.49107, 1.49318] (0.95 confidence) | 1.48308 ns [1.48053, 1.48422] (0.95 confidence) | — | — |
| linear&#95;form&#95;d2&#95;prepared/underflow/dot&#95;endpoints | 1.55381 ns [1.55122, 1.55799] (0.95 confidence) | 1.4842 ns [1.48381, 1.48429] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;batch/origin/construction | 2.04723 ns [2.0417, 2.05634] (0.95 confidence) | 2.26383 ns [2.23912, 2.27037] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;batch/origin/prepared | 24.6288 ns [24.5766, 24.6846] (0.95 confidence) | 23.9219 ns [23.843, 24.0564] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;batch/origin/with&#95;construction | 28.6987 ns [28.5649, 28.8529] (0.95 confidence) | 28.2668 ns [28.1989, 28.4733] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;batch/translated/construction | 2.05769 ns [2.05111, 2.0635] (0.95 confidence) | 2.2681 ns [2.24845, 2.27505] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;batch/translated/prepared | 82.2446 ns [82.086, 82.3434] (0.95 confidence) | 68.1149 ns [68.0075, 68.1918] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;batch/translated/with&#95;construction | 85.7979 ns [85.5925, 86.2242] (0.95 confidence) | 71.0537 ns [71.0012, 71.1221] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/cancellation/difference | 17.3751 ns [17.3499, 17.4097] (0.95 confidence) | 10.8743 ns [10.8681, 10.8792] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/cancellation/difference&#95;endpoints | 18.6192 ns [18.5896, 18.6561] (0.95 confidence) | 11.8557 ns [11.826, 11.8924] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/cancellation/dot | 7.07716 ns [7.06399, 7.0942] (0.95 confidence) | 6.05745 ns [6.05371, 6.06194] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/cancellation/dot&#95;endpoints | 7.87255 ns [7.86055, 7.9256] (0.95 confidence) | 6.60638 ns [6.59024, 6.63641] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/dense/difference | 23.6731 ns [23.621, 23.7117] (0.95 confidence) | 14.2878 ns [14.2691, 14.315] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/dense/difference&#95;endpoints | 26.982 ns [26.9309, 27.0189] (0.95 confidence) | 16.4802 ns [16.4533, 16.4988] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/dense/dot | 9.63411 ns [9.61363, 9.64974] (0.95 confidence) | 7.98182 ns [7.95075, 7.98477] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/dense/dot&#95;endpoints | 10.4091 ns [10.3639, 10.4421] (0.95 confidence) | 8.55832 ns [8.49535, 9.2153] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/dense&#95;inexact/difference | 23.7184 ns [23.6739, 23.7692] (0.95 confidence) | 14.3481 ns [14.3064, 14.3605] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/dense&#95;inexact/difference&#95;endpoints | 25.3482 ns [25.291, 25.3832] (0.95 confidence) | 15.9915 ns [15.9538, 16.0145] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/dense&#95;inexact/dot | 9.52741 ns [9.48898, 9.56149] (0.95 confidence) | 7.95183 ns [7.92919, 7.97192] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/dense&#95;inexact/dot&#95;endpoints | 11.1805 ns [11.1664, 11.2236] (0.95 confidence) | 9.19935 ns [9.18163, 9.21633] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/overflow/dot | 3.09807 ns [3.09287, 3.10595] (0.95 confidence) | 1.00061 ns [0.999054, 1.00264] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/sparse/difference | 8.4509 ns [8.44749, 8.45716] (0.95 confidence) | 5.45126 ns [5.44975, 5.45482] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/sparse/difference&#95;endpoints | 9.21662 ns [9.2087, 9.22196] (0.95 confidence) | 6.38118 ns [6.37328, 6.39877] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/sparse/dot | 2.52497 ns [2.52097, 2.53306] (0.95 confidence) | 2.82018 ns [2.76048, 2.94905] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/sparse/dot&#95;endpoints | 2.78413 ns [2.77286, 2.79924] (0.95 confidence) | 2.74303 ns [2.74244, 2.74435] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/underflow/difference | 7.67494 ns [7.66781, 7.68247] (0.95 confidence) | 4.70284 ns [4.70221, 4.70332] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/underflow/difference&#95;endpoints | 7.74478 ns [7.73944, 7.75115] (0.95 confidence) | 4.96773 ns [4.9669, 4.96826] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/underflow/dot | 3.45048 ns [3.44488, 3.45505] (0.95 confidence) | 3.1525 ns [3.14547, 3.1566] (0.95 confidence) | — | — |
| linear&#95;form&#95;d3&#95;prepared/underflow/dot&#95;endpoints | 3.47158 ns [3.46518, 3.47752] (0.95 confidence) | 3.15478 ns [3.15053, 3.15931] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4/dot&#95;bounded&#95;inconclusive | 7.21617 ns [7.20781, 7.22163] (0.95 confidence) | 6.22429 ns [6.17994, 6.2657] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4/dot&#95;bounded&#95;well&#95;separated | 11.3947 ns [11.3527, 11.4154] (0.95 confidence) | 9.34147 ns [9.33208, 9.35244] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4/dot&#95;difference&#95;bounded&#95;inconclusive | 31.4313 ns [31.3502, 31.4881] (0.95 confidence) | 20.581 ns [20.5465, 20.6254] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4/dot&#95;difference&#95;bounded&#95;well&#95;separated | 27.4542 ns [27.4075, 27.5341] (0.95 confidence) | 16.3667 ns [16.3197, 16.3859] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4/dot&#95;plain&#95;well&#95;separated | 0.679843 ns [0.678655, 0.6816] (0.95 confidence) | 0.678392 ns [0.677916, 0.681955] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;batch/origin/construction | 2.62588 ns [2.62046, 2.64329] (0.95 confidence) | 2.34627 ns [2.33189, 2.36199] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;batch/origin/prepared | 40.2474 ns [40.1936, 40.3072] (0.95 confidence) | 34.9296 ns [34.8982, 35.0586] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;batch/origin/with&#95;construction | 46.8533 ns [46.8077, 46.9116] (0.95 confidence) | 42.0101 ns [41.6295, 42.1819] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;batch/translated/construction | 2.93381 ns [2.92906, 2.9486] (0.95 confidence) | 2.34344 ns [2.33339, 2.35924] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;batch/translated/prepared | 133.086 ns [132.81, 133.27] (0.95 confidence) | 107.258 ns [107.143, 107.524] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;batch/translated/with&#95;construction | 139.013 ns [138.688, 139.368] (0.95 confidence) | 113.133 ns [113.021, 113.258] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/cancellation/difference | 18.1805 ns [18.1566, 18.2203] (0.95 confidence) | 11.5377 ns [11.4934, 11.5963] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/cancellation/difference&#95;endpoints | 19.7901 ns [19.7684, 19.8252] (0.95 confidence) | 12.7465 ns [12.7081, 12.7708] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/cancellation/dot | 7.44414 ns [7.42769, 7.45715] (0.95 confidence) | 6.26037 ns [6.24677, 6.2651] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/cancellation/dot&#95;endpoints | 8.19466 ns [8.18437, 8.21799] (0.95 confidence) | 6.93066 ns [6.92249, 6.93794] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/dense/difference | 31.2988 ns [31.1979, 31.3472] (0.95 confidence) | 18.42 ns [18.3777, 18.4343] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/dense/difference&#95;endpoints | 35.0061 ns [34.8977, 35.0634] (0.95 confidence) | 20.8224 ns [20.792, 20.8338] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/dense/dot | 11.4357 ns [11.4136, 11.445] (0.95 confidence) | 9.52813 ns [9.52001, 9.54033] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/dense/dot&#95;endpoints | 13.5532 ns [13.5044, 13.5825] (0.95 confidence) | 11.1165 ns [11.0807, 11.1461] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/dense&#95;inexact/difference | 31.4597 ns [31.3589, 31.5227] (0.95 confidence) | 18.4214 ns [18.3823, 18.4575] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/dense&#95;inexact/difference&#95;endpoints | 32.9333 ns [32.8339, 33.0308] (0.95 confidence) | 20.0763 ns [20.04, 20.1273] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/dense&#95;inexact/dot | 11.4249 ns [11.4045, 11.4441] (0.95 confidence) | 9.61746 ns [9.57981, 9.72216] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/dense&#95;inexact/dot&#95;endpoints | 12.9509 ns [12.9277, 12.974] (0.95 confidence) | 10.2156 ns [10.1842, 10.2381] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/overflow/dot | 3.02682 ns [3.01612, 3.03786] (0.95 confidence) | 0.99782 ns [0.996907, 0.999865] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/sparse/difference | 10.3236 ns [10.3092, 10.3601] (0.95 confidence) | 6.19055 ns [6.1849, 6.1953] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/sparse/difference&#95;endpoints | 11.087 ns [11.0613, 11.0943] (0.95 confidence) | 7.39863 ns [7.3948, 7.40844] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/sparse/dot | 3.20814 ns [3.20589, 3.21147] (0.95 confidence) | 3.11788 ns [3.11732, 3.12586] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/sparse/dot&#95;endpoints | 3.40923 ns [3.40714, 3.41164] (0.95 confidence) | 3.20377 ns [3.19947, 3.2107] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/underflow/difference | 9.41151 ns [9.40364, 9.41666] (0.95 confidence) | 5.6906 ns [5.68995, 5.69237] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/underflow/difference&#95;endpoints | 9.52246 ns [9.50998, 9.5302] (0.95 confidence) | 5.95253 ns [5.9453, 5.97283] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/underflow/dot | 4.60016 ns [4.59187, 4.60296] (0.95 confidence) | 3.85826 ns [3.85371, 3.86538] (0.95 confidence) | — | — |
| linear&#95;form&#95;d4&#95;prepared/underflow/dot&#95;endpoints | 4.16018 ns [4.15384, 4.16305] (0.95 confidence) | 3.8802 ns [3.87343, 3.89116] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;batch/origin/construction | 3.08483 ns [3.07832, 3.08827] (0.95 confidence) | 2.88478 ns [2.87615, 2.88794] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;batch/origin/prepared | 55.2354 ns [55.2054, 55.2704] (0.95 confidence) | 46.1407 ns [46.1339, 46.1483] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;batch/origin/with&#95;construction | 66.5723 ns [66.5051, 66.6437] (0.95 confidence) | 57.6837 ns [57.648, 58.3603] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;batch/translated/construction | 3.42943 ns [3.42208, 3.46142] (0.95 confidence) | 3.22488 ns [3.21661, 3.23656] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;batch/translated/prepared | 180.664 ns [180.405, 180.836] (0.95 confidence) | 144.946 ns [144.561, 145.214] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;batch/translated/with&#95;construction | 189.404 ns [189.143, 189.946] (0.95 confidence) | 154.838 ns [154.61, 154.975] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/cancellation/difference | 19.8938 ns [19.8779, 19.9389] (0.95 confidence) | 12.3023 ns [12.2497, 12.4026] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/cancellation/difference&#95;endpoints | 22.2129 ns [22.1497, 22.2507] (0.95 confidence) | 13.3216 ns [13.2779, 13.4001] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/cancellation/dot | 7.4775 ns [7.45963, 7.49103] (0.95 confidence) | 6.66802 ns [6.63204, 6.6924] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/cancellation/dot&#95;endpoints | 8.56899 ns [8.55217, 8.57757] (0.95 confidence) | 7.2151 ns [7.2059, 7.22597] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/dense/difference | 38.3619 ns [38.2415, 38.4125] (0.95 confidence) | 22.6649 ns [22.6257, 22.7377] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/dense/difference&#95;endpoints | 40.0011 ns [39.9, 40.131] (0.95 confidence) | 24.1401 ns [24.0954, 24.1638] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/dense/dot | 14.1278 ns [14.086, 14.1729] (0.95 confidence) | 11.5709 ns [11.5602, 11.5919] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/dense/dot&#95;endpoints | 15.1225 ns [15.0697, 15.1576] (0.95 confidence) | 12.2147 ns [12.1998, 12.2245] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/dense&#95;inexact/difference | 38.5807 ns [38.4877, 38.6255] (0.95 confidence) | 22.7475 ns [22.6996, 22.774] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/dense&#95;inexact/difference&#95;endpoints | 43.7804 ns [43.6847, 43.8342] (0.95 confidence) | 25.2227 ns [25.1813, 25.2623] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/dense&#95;inexact/dot | 14.106 ns [14.0793, 14.1332] (0.95 confidence) | 11.5937 ns [11.576, 11.6077] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/dense&#95;inexact/dot&#95;endpoints | 16.5418 ns [16.4594, 16.6009] (0.95 confidence) | 13.0849 ns [13.0729, 13.1072] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/overflow/dot | 2.99546 ns [2.98948, 3.00368] (0.95 confidence) | 0.997011 ns [0.994281, 0.999465] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/sparse/difference | 11.8954 ns [11.8753, 11.9162] (0.95 confidence) | 7.21053 ns [7.20753, 7.21432] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/sparse/difference&#95;endpoints | 12.8397 ns [12.8257, 12.856] (0.95 confidence) | 8.35547 ns [8.33652, 8.38226] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/sparse/dot | 4.1031 ns [4.09894, 4.10918] (0.95 confidence) | 3.62528 ns [3.61972, 3.63653] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/sparse/dot&#95;endpoints | 4.27209 ns [4.26798, 4.28243] (0.95 confidence) | 3.73959 ns [3.73713, 3.74325] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/underflow/difference | 11.1308 ns [11.1175, 11.1475] (0.95 confidence) | 7.58073 ns [7.57634, 7.58838] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/underflow/difference&#95;endpoints | 11.2173 ns [11.2022, 11.2269] (0.95 confidence) | 7.82203 ns [7.81129, 7.82462] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/underflow/dot | 4.86774 ns [4.86219, 4.87362] (0.95 confidence) | 4.57918 ns [4.57669, 4.58063] (0.95 confidence) | — | — |
| linear&#95;form&#95;d5&#95;prepared/underflow/dot&#95;endpoints | 4.87561 ns [4.86651, 4.88279] (0.95 confidence) | 4.5784 ns [4.5701, 4.58584] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;batch/origin/construction | 3.46353 ns [3.45982, 3.46755] (0.95 confidence) | 3.08543 ns [3.07724, 3.08933] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;batch/origin/prepared | 78.7062 ns [78.561, 78.7732] (0.95 confidence) | 60.6998 ns [60.5898, 60.7498] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;batch/origin/with&#95;construction | 101.433 ns [101.317, 101.652] (0.95 confidence) | 79.9721 ns [79.8775, 80.0872] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;batch/translated/construction | 3.61175 ns [3.6082, 3.61769] (0.95 confidence) | 3.54985 ns [3.54674, 3.55224] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;batch/translated/prepared | 253.906 ns [253.445, 254.621] (0.95 confidence) | 201.205 ns [200.988, 201.34] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;batch/translated/with&#95;construction | 276.775 ns [276.294, 277.224] (0.95 confidence) | 225.303 ns [224.983, 225.804] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/cancellation/difference | 22.7307 ns [22.6556, 22.7879] (0.95 confidence) | 12.6548 ns [12.6329, 12.6946] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/cancellation/difference&#95;endpoints | 24.4986 ns [24.4421, 24.6098] (0.95 confidence) | 13.9891 ns [13.963, 14.0394] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/cancellation/dot | 8.07428 ns [8.03552, 8.10877] (0.95 confidence) | 6.89947 ns [6.8828, 6.92326] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/cancellation/dot&#95;endpoints | 9.49806 ns [9.47757, 9.51853] (0.95 confidence) | 7.60571 ns [7.5935, 7.63157] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/dense/difference | 48.0631 ns [47.8902, 48.17] (0.95 confidence) | 28.6891 ns [28.6005, 28.7541] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/dense/difference&#95;endpoints | 53.6343 ns [53.5327, 53.708] (0.95 confidence) | 32.1441 ns [32.0631, 32.2255] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/dense/dot | 16.5236 ns [16.5168, 16.5793] (0.95 confidence) | 13.6775 ns [13.67, 13.7073] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/dense/dot&#95;endpoints | 19.4328 ns [19.4178, 19.4897] (0.95 confidence) | 14.9982 ns [14.9739, 15.0317] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/dense&#95;inexact/difference | 48.0673 ns [47.9679, 48.2276] (0.95 confidence) | 28.717 ns [28.6267, 28.7563] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/dense&#95;inexact/difference&#95;endpoints | 49.84 ns [49.7186, 49.94] (0.95 confidence) | 31.4739 ns [31.3303, 31.5295] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/dense&#95;inexact/dot | 16.7695 ns [16.7345, 16.8259] (0.95 confidence) | 13.7623 ns [13.7316, 13.7775] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/dense&#95;inexact/dot&#95;endpoints | 19.359 ns [19.0892, 19.4532] (0.95 confidence) | 15.0403 ns [15.0177, 15.0485] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/overflow/dot | 3.05125 ns [3.0445, 3.0655] (0.95 confidence) | 0.996395 ns [0.994067, 0.99948] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/sparse/difference | 13.6286 ns [13.6173, 13.6412] (0.95 confidence) | 8.08297 ns [8.07587, 8.0901] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/sparse/difference&#95;endpoints | 14.6575 ns [14.6297, 14.6782] (0.95 confidence) | 9.43665 ns [9.42472, 9.4536] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/sparse/dot | 5.02101 ns [5.01252, 5.02646] (0.95 confidence) | 4.10064 ns [4.09728, 4.10498] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/sparse/dot&#95;endpoints | 5.21993 ns [5.21469, 5.22835] (0.95 confidence) | 4.18337 ns [4.17688, 4.18497] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/underflow/difference | 12.5878 ns [12.5672, 12.608] (0.95 confidence) | 8.55042 ns [8.54282, 8.55557] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/underflow/difference&#95;endpoints | 13.2039 ns [13.1877, 13.2217] (0.95 confidence) | 8.82416 ns [8.80345, 8.83588] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/underflow/dot | 5.63023 ns [5.62236, 5.63986] (0.95 confidence) | 5.30179 ns [5.29528, 5.30354] (0.95 confidence) | — | — |
| linear&#95;form&#95;d6&#95;prepared/underflow/dot&#95;endpoints | 6.04373 ns [6.03616, 6.05342] (0.95 confidence) | 5.35698 ns [5.34605, 5.36259] (0.95 confidence) | — | — |
| polynomial/coincident&#95;lifted/order4/determinant | — | 13030.6 ns [13009.9, 13058.2] (0.95 confidence) | — | — |
| polynomial/coincident&#95;lifted/order5/determinant | — | 103419 ns [103216, 103771] (0.95 confidence) | — | — |
| polynomial/coincident&#95;lifted/order6/determinant | — | 894694 ns [893193, 895993] (0.95 confidence) | — | — |
| polynomial/coincident&#95;lifted/order7/determinant | — | 8.48315e+06 ns [8.45394e+06, 8.49939e+06] (0.95 confidence) | — | — |
| polynomial/coincident&#95;lifted/order8/determinant | — | 8.48231e+07 ns [8.47308e+07, 8.50018e+07] (0.95 confidence) | — | — |
| public&#95;d2/det&#95;bound | 1.95425 ns [1.95244, 1.95537] (0.95 confidence) | 1.95895 ns [1.9551, 1.9616] (0.95 confidence) | — | — |
| public&#95;d2/det&#95;direct | 1.12917 ns [1.12633, 1.13029] (0.95 confidence) | 1.12452 ns [1.1231, 1.12541] (0.95 confidence) | — | — |
| public&#95;d2/det&#95;with&#95;bound | 2.11797 ns [2.11285, 2.12247] (0.95 confidence) | 2.11745 ns [2.11292, 2.12205] (0.95 confidence) | — | — |
| public&#95;d2/first&#95;asymmetry | 1.38189 ns [1.37921, 1.38291] (0.95 confidence) | 1.37444 ns [1.37279, 1.37567] (0.95 confidence) | — | — |
| public&#95;d2/interval&#95;det | 8.09126 ns [7.9378, 8.31869] (0.95 confidence) | 7.91504 ns [7.89827, 7.92484] (0.95 confidence) | — | — |
| public&#95;d2/interval&#95;matrix&#95;access | 6.26419 ns [6.24464, 6.35671] (0.95 confidence) | 6.24263 ns [6.22866, 6.25651] (0.95 confidence) | — | — |
| public&#95;d2/interval&#95;matrix&#95;from&#95;matrix | 2.58664 ns [2.57931, 2.60184] (0.95 confidence) | 2.55824 ns [2.55089, 2.57772] (0.95 confidence) | — | — |
| public&#95;d2/interval&#95;matrix&#95;from&#95;points | 2.33149 ns [2.32566, 2.33818] (0.95 confidence) | 2.23866 ns [2.21935, 2.25256] (0.95 confidence) | — | — |
| public&#95;d2/interval&#95;matrix&#95;from&#95;rows | 3.54607 ns [3.53484, 3.5527] (0.95 confidence) | 3.50303 ns [3.49156, 3.51644] (0.95 confidence) | — | — |
| public&#95;d2/interval&#95;matrix&#95;identity | 2.07738 ns [2.07349, 2.08032] (0.95 confidence) | 2.06674 ns [2.06412, 2.06892] (0.95 confidence) | — | — |
| public&#95;d2/interval&#95;matrix&#95;set | 3.7843 ns [3.77236, 3.79075] (0.95 confidence) | 3.72489 ns [3.7185, 3.73218] (0.95 confidence) | — | — |
| public&#95;d2/interval&#95;matrix&#95;zero | 2.10363 ns [2.09917, 2.10815] (0.95 confidence) | 2.07562 ns [2.06765, 2.08104] (0.95 confidence) | — | — |
| public&#95;d2/is&#95;symmetric | 0.755677 ns [0.755472, 0.755938] (0.95 confidence) | 0.756872 ns [0.756307, 0.757859] (0.95 confidence) | — | — |
| public&#95;d2/matrix&#95;access | 4.6886 ns [4.68114, 4.6978] (0.95 confidence) | 4.69217 ns [4.68179, 4.70555] (0.95 confidence) | — | — |
| public&#95;d2/matrix&#95;from&#95;rows | 2.89757 ns [2.89459, 2.90035] (0.95 confidence) | 2.89067 ns [2.88556, 2.89601] (0.95 confidence) | — | — |
| public&#95;d2/matrix&#95;identity | 1.23312 ns [1.2305, 1.2354] (0.95 confidence) | 1.23545 ns [1.23269, 1.23842] (0.95 confidence) | — | — |
| public&#95;d2/matrix&#95;set | 2.08447 ns [2.07161, 2.09592] (0.95 confidence) | 2.08003 ns [2.07281, 2.08918] (0.95 confidence) | — | — |
| public&#95;d2/matrix&#95;zero | 1.1521 ns [1.14921, 1.15535] (0.95 confidence) | 1.14562 ns [1.14385, 1.1472] (0.95 confidence) | — | — |
| public&#95;d2/vector&#95;access | 1.60971 ns [1.60638, 1.61372] (0.95 confidence) | 1.59905 ns [1.59566, 1.60338] (0.95 confidence) | — | — |
| public&#95;d2/vector&#95;new | 1.93678 ns [1.89933, 1.9412] (0.95 confidence) | 1.92229 ns [1.91485, 1.927] (0.95 confidence) | — | — |
| public&#95;d2/vector&#95;zero | 0.586456 ns [0.585882, 0.58723] (0.95 confidence) | 0.581542 ns [0.579749, 0.582719] (0.95 confidence) | — | — |
| public&#95;d3/det&#95;bound | 3.64688 ns [3.6447, 3.64996] (0.95 confidence) | 3.64757 ns [3.64249, 3.65373] (0.95 confidence) | — | — |
| public&#95;d3/det&#95;direct | 1.28439 ns [1.28379, 1.28475] (0.95 confidence) | 1.26236 ns [1.26183, 1.26315] (0.95 confidence) | — | — |
| public&#95;d3/det&#95;with&#95;bound | 3.67755 ns [3.67678, 3.67822] (0.95 confidence) | 3.66781 ns [3.66309, 3.6749] (0.95 confidence) | — | — |
| public&#95;d3/first&#95;asymmetry | 3.56833 ns [3.5443, 3.60893] (0.95 confidence) | 3.50134 ns [3.45955, 3.55895] (0.95 confidence) | — | — |
| public&#95;d3/interval&#95;det | 24.4149 ns [24.1767, 24.6827] (0.95 confidence) | 23.9405 ns [23.6131, 24.1218] (0.95 confidence) | — | — |
| public&#95;d3/interval&#95;matrix&#95;access | 9.76427 ns [9.73542, 9.79009] (0.95 confidence) | 9.76897 ns [9.73797, 9.78462] (0.95 confidence) | — | — |
| public&#95;d3/interval&#95;matrix&#95;from&#95;matrix | 5.16557 ns [5.14964, 5.17248] (0.95 confidence) | 5.11108 ns [5.09366, 5.12121] (0.95 confidence) | — | — |
| public&#95;d3/interval&#95;matrix&#95;from&#95;points | 6.3544 ns [6.34282, 6.36653] (0.95 confidence) | 6.27696 ns [6.26707, 6.28564] (0.95 confidence) | — | — |
| public&#95;d3/interval&#95;matrix&#95;from&#95;rows | 7.34437 ns [7.33049, 7.35934] (0.95 confidence) | 7.21547 ns [7.19611, 7.23171] (0.95 confidence) | — | — |
| public&#95;d3/interval&#95;matrix&#95;identity | 4.29613 ns [4.28949, 4.30178] (0.95 confidence) | 4.20755 ns [4.20112, 4.21326] (0.95 confidence) | — | — |
| public&#95;d3/interval&#95;matrix&#95;set | 7.79235 ns [7.77189, 7.8281] (0.95 confidence) | 7.6999 ns [7.67101, 7.74614] (0.95 confidence) | — | — |
| public&#95;d3/interval&#95;matrix&#95;zero | 4.27136 ns [4.26442, 4.28598] (0.95 confidence) | 4.21592 ns [4.2072, 4.22161] (0.95 confidence) | — | — |
| public&#95;d3/is&#95;symmetric | 1.339 ns [1.33749, 1.3402] (0.95 confidence) | 1.48588 ns [1.48466, 1.48759] (0.95 confidence) | — | — |
| public&#95;d3/matrix&#95;access | 6.63971 ns [6.61471, 6.6523] (0.95 confidence) | 6.54981 ns [6.532, 6.5693] (0.95 confidence) | — | — |
| public&#95;d3/matrix&#95;from&#95;rows | 5.7303 ns [5.72054, 5.73856] (0.95 confidence) | 5.65897 ns [5.64087, 5.67208] (0.95 confidence) | — | — |
| public&#95;d3/matrix&#95;identity | 2.67776 ns [2.67267, 2.68084] (0.95 confidence) | 2.65578 ns [2.65071, 2.66104] (0.95 confidence) | — | — |
| public&#95;d3/matrix&#95;set | 4.2145 ns [4.19309, 4.22338] (0.95 confidence) | 4.16835 ns [4.14388, 4.17953] (0.95 confidence) | — | — |
| public&#95;d3/matrix&#95;zero | 2.38912 ns [2.38387, 2.39689] (0.95 confidence) | 2.37447 ns [2.36949, 2.38071] (0.95 confidence) | — | — |
| public&#95;d3/vector&#95;access | 2.07756 ns [2.07132, 2.08318] (0.95 confidence) | 2.05052 ns [2.04407, 2.05662] (0.95 confidence) | — | — |
| public&#95;d3/vector&#95;new | 3.98657 ns [3.97069, 4.01596] (0.95 confidence) | 3.64318 ns [3.62949, 3.64718] (0.95 confidence) | — | — |
| public&#95;d3/vector&#95;zero | 0.946915 ns [0.945268, 0.948583] (0.95 confidence) | 0.934631 ns [0.932032, 0.936264] (0.95 confidence) | — | — |
| public&#95;d4/det&#95;bound | 6.24464 ns [6.23236, 6.26586] (0.95 confidence) | 6.1848 ns [6.18273, 6.18949] (0.95 confidence) | — | — |
| public&#95;d4/det&#95;direct | 2.48869 ns [2.48757, 2.49036] (0.95 confidence) | 2.47297 ns [2.46852, 2.47438] (0.95 confidence) | — | — |
| public&#95;d4/det&#95;with&#95;bound | 6.23827 ns [6.23048, 6.24361] (0.95 confidence) | 6.18261 ns [6.17356, 6.18442] (0.95 confidence) | — | — |
| public&#95;d4/first&#95;asymmetry | 3.62038 ns [3.59441, 3.65551] (0.95 confidence) | 3.61452 ns [3.58072, 3.6613] (0.95 confidence) | — | — |
| public&#95;d4/interval&#95;det | 51.539 ns [51.1941, 51.8543] (0.95 confidence) | 50.655 ns [50.3573, 51.1149] (0.95 confidence) | — | — |
| public&#95;d4/interval&#95;matrix&#95;access | 14.4137 ns [14.1949, 14.5647] (0.95 confidence) | 13.5203 ns [13.4075, 13.5797] (0.95 confidence) | — | — |
| public&#95;d4/interval&#95;matrix&#95;from&#95;matrix | 8.02271 ns [8.01004, 8.03448] (0.95 confidence) | 7.80369 ns [7.78684, 7.81999] (0.95 confidence) | — | — |
| public&#95;d4/interval&#95;matrix&#95;from&#95;points | 11.0549 ns [11.0468, 11.071] (0.95 confidence) | 10.7687 ns [10.7475, 10.7936] (0.95 confidence) | — | — |
| public&#95;d4/interval&#95;matrix&#95;from&#95;rows | 11.4026 ns [11.3778, 11.4389] (0.95 confidence) | 10.9662 ns [10.9521, 11.0052] (0.95 confidence) | — | — |
| public&#95;d4/interval&#95;matrix&#95;identity | 6.63401 ns [6.61322, 6.65297] (0.95 confidence) | 6.33923 ns [6.32332, 6.36298] (0.95 confidence) | — | — |
| public&#95;d4/interval&#95;matrix&#95;set | 11.5963 ns [11.5722, 11.6193] (0.95 confidence) | 11.2784 ns [11.252, 11.298] (0.95 confidence) | — | — |
| public&#95;d4/interval&#95;matrix&#95;zero | 6.68162 ns [6.67137, 6.69728] (0.95 confidence) | 6.46484 ns [6.45279, 6.47488] (0.95 confidence) | — | — |
| public&#95;d4/is&#95;symmetric | 2.11935 ns [2.1183, 2.12068] (0.95 confidence) | 2.10272 ns [2.10246, 2.10309] (0.95 confidence) | — | — |
| public&#95;d4/matrix&#95;access | 9.30408 ns [9.28874, 9.32736] (0.95 confidence) | 9.2135 ns [9.18459, 9.24482] (0.95 confidence) | — | — |
| public&#95;d4/matrix&#95;from&#95;rows | 8.7149 ns [8.69736, 8.72936] (0.95 confidence) | 10.2377 ns [10.2289, 10.2539] (0.95 confidence) | — | — |
| public&#95;d4/matrix&#95;identity | 4.30407 ns [4.29126, 4.36429] (0.95 confidence) | 4.02212 ns [4.01759, 4.03066] (0.95 confidence) | — | — |
| public&#95;d4/matrix&#95;set | 7.52521 ns [7.50998, 7.5448] (0.95 confidence) | 7.20915 ns [7.19567, 7.22534] (0.95 confidence) | — | — |
| public&#95;d4/matrix&#95;zero | 3.97624 ns [3.9679, 3.9903] (0.95 confidence) | 3.7809 ns [3.77106, 3.78694] (0.95 confidence) | — | — |
| public&#95;d4/vector&#95;access | 2.65694 ns [2.64768, 2.66204] (0.95 confidence) | 2.56494 ns [2.56037, 2.57503] (0.95 confidence) | — | — |
| public&#95;d4/vector&#95;new | 2.96576 ns [2.95929, 2.97282] (0.95 confidence) | 2.81823 ns [2.80958, 2.82741] (0.95 confidence) | — | — |
| public&#95;d4/vector&#95;zero | 1.19209 ns [1.1898, 1.19378] (0.95 confidence) | 1.14148 ns [1.13882, 1.14246] (0.95 confidence) | — | — |
| public&#95;d5/det&#95;bound | 1.35434 ns [1.35129, 1.35698] (0.95 confidence) | 1.31094 ns [1.31041, 1.31154] (0.95 confidence) | — | — |
| public&#95;d5/det&#95;direct | 1.35665 ns [1.35454, 1.35817] (0.95 confidence) | 1.31114 ns [1.31078, 1.31189] (0.95 confidence) | — | — |
| public&#95;d5/det&#95;with&#95;bound | 0.860717 ns [0.857805, 0.864817] (0.95 confidence) | 0.835909 ns [0.834723, 0.838044] (0.95 confidence) | — | — |
| public&#95;d5/first&#95;asymmetry | 3.80618 ns [3.80402, 3.80952] (0.95 confidence) | 3.77854 ns [3.77541, 3.78171] (0.95 confidence) | — | — |
| public&#95;d5/interval&#95;det | 118.138 ns [117.423, 118.662] (0.95 confidence) | 122.724 ns [121.814, 123.48] (0.95 confidence) | — | — |
| public&#95;d5/interval&#95;matrix&#95;access | 20.2676 ns [19.5172, 20.4898] (0.95 confidence) | 18.2863 ns [18.2489, 18.3012] (0.95 confidence) | — | — |
| public&#95;d5/interval&#95;matrix&#95;from&#95;matrix | 10.7524 ns [10.6281, 10.8313] (0.95 confidence) | 10.4821 ns [10.3999, 10.6367] (0.95 confidence) | — | — |
| public&#95;d5/interval&#95;matrix&#95;from&#95;points | 18.7428 ns [18.7211, 18.7826] (0.95 confidence) | 18.3015 ns [18.2478, 18.3221] (0.95 confidence) | — | — |
| public&#95;d5/interval&#95;matrix&#95;from&#95;rows | 14.4726 ns [14.0187, 15.1463] (0.95 confidence) | 15.6187 ns [14.7541, 16.0549] (0.95 confidence) | — | — |
| public&#95;d5/interval&#95;matrix&#95;identity | 8.59183 ns [8.10371, 8.69373] (0.95 confidence) | 8.39727 ns [7.97745, 8.57111] (0.95 confidence) | — | — |
| public&#95;d5/interval&#95;matrix&#95;set | 13.7304 ns [13.3479, 14.134] (0.95 confidence) | 15.7054 ns [14.9857, 15.9636] (0.95 confidence) | — | — |
| public&#95;d5/interval&#95;matrix&#95;zero | 9.45084 ns [9.44119, 9.46629] (0.95 confidence) | 9.25117 ns [9.2158, 9.27348] (0.95 confidence) | — | — |
| public&#95;d5/is&#95;symmetric | 3.70791 ns [3.7062, 3.70919] (0.95 confidence) | 3.68636 ns [3.68601, 3.68696] (0.95 confidence) | — | — |
| public&#95;d5/matrix&#95;access | 11.8913 ns [11.6059, 12.008] (0.95 confidence) | 12.0077 ns [11.9604, 12.1498] (0.95 confidence) | — | — |
| public&#95;d5/matrix&#95;from&#95;rows | 11.7826 ns [11.7604, 11.791] (0.95 confidence) | 11.7007 ns [11.6832, 11.713] (0.95 confidence) | — | — |
| public&#95;d5/matrix&#95;identity | 6.90549 ns [6.88292, 6.93271] (0.95 confidence) | 6.67087 ns [6.63139, 6.68503] (0.95 confidence) | — | — |
| public&#95;d5/matrix&#95;set | 10.6169 ns [10.5873, 10.6393] (0.95 confidence) | 10.3242 ns [10.2998, 10.3555] (0.95 confidence) | — | — |
| public&#95;d5/matrix&#95;zero | 5.82588 ns [5.8149, 5.83072] (0.95 confidence) | 5.6095 ns [5.60443, 5.61542] (0.95 confidence) | — | — |
| public&#95;d5/vector&#95;access | 3.11107 ns [3.10221, 3.12675] (0.95 confidence) | 3.04164 ns [3.03167, 3.04916] (0.95 confidence) | — | — |
| public&#95;d5/vector&#95;new | 4.83949 ns [4.82131, 4.85392] (0.95 confidence) | 4.84082 ns [4.8305, 4.85153] (0.95 confidence) | — | — |
| public&#95;d5/vector&#95;zero | 1.55075 ns [1.53741, 1.56388] (0.95 confidence) | 1.49081 ns [1.47346, 1.4996] (0.95 confidence) | — | — |
| public&#95;exact&#95;d2/matrix&#95;access | 0.783143 ns [0.781318, 0.785793] (0.95 confidence) | 0.754298 ns [0.75108, 0.757038] (0.95 confidence) | — | — |
| public&#95;exact&#95;d2/matrix&#95;from&#95;fn | 135.36 ns [135.184, 135.535] (0.95 confidence) | 130.679 ns [130.37, 130.893] (0.95 confidence) | — | — |
| public&#95;exact&#95;d2/matrix&#95;from&#95;rows | 138.642 ns [138.295, 138.855] (0.95 confidence) | 134.021 ns [133.834, 134.241] (0.95 confidence) | — | — |
| public&#95;exact&#95;d2/matrix&#95;into&#95;rows | 66.9789 ns [66.8326, 67.161] (0.95 confidence) | 64.2995 ns [64.129, 64.5517] (0.95 confidence) | — | — |
| public&#95;exact&#95;d2/matrix&#95;zero | 2.74706 ns [2.74245, 2.75283] (0.95 confidence) | 2.68829 ns [2.6845, 2.69262] (0.95 confidence) | — | — |
| public&#95;exact&#95;d2/vector&#95;access | 0.785684 ns [0.782087, 0.78743] (0.95 confidence) | 0.753859 ns [0.751912, 0.757063] (0.95 confidence) | — | — |
| public&#95;exact&#95;d2/vector&#95;from&#95;fn | 38.6122 ns [38.5562, 38.6718] (0.95 confidence) | 36.6532 ns [36.5512, 36.7093] (0.95 confidence) | — | — |
| public&#95;exact&#95;d2/vector&#95;into&#95;array | 16.0256 ns [16.0008, 16.0846] (0.95 confidence) | 15.3789 ns [15.3337, 15.4315] (0.95 confidence) | — | — |
| public&#95;exact&#95;d2/vector&#95;new | 36.7692 ns [36.6675, 36.8277] (0.95 confidence) | 35.922 ns [35.8463, 35.9597] (0.95 confidence) | — | — |
| public&#95;exact&#95;d2/vector&#95;zero | 1.37804 ns [1.37719, 1.37907] (0.95 confidence) | 1.3533 ns [1.35183, 1.35417] (0.95 confidence) | — | — |
| public&#95;exact&#95;d3/matrix&#95;access | 0.78168 ns [0.779555, 0.784657] (0.95 confidence) | 0.755696 ns [0.753848, 0.758525] (0.95 confidence) | — | — |
| public&#95;exact&#95;d3/matrix&#95;from&#95;fn | 233.054 ns [232.763, 233.487] (0.95 confidence) | 228.414 ns [227.856, 228.971] (0.95 confidence) | — | — |
| public&#95;exact&#95;d3/matrix&#95;from&#95;rows | 272.828 ns [272.458, 273.097] (0.95 confidence) | 265.746 ns [265.011, 265.98] (0.95 confidence) | — | — |
| public&#95;exact&#95;d3/matrix&#95;into&#95;rows | 125.906 ns [125.658, 126.155] (0.95 confidence) | 122.775 ns [122.589, 122.971] (0.95 confidence) | — | — |
| public&#95;exact&#95;d3/matrix&#95;zero | 6.3692 ns [6.3577, 6.37923] (0.95 confidence) | 6.16586 ns [6.15687, 6.17716] (0.95 confidence) | — | — |
| public&#95;exact&#95;d3/vector&#95;access | 0.779548 ns [0.778962, 0.781404] (0.95 confidence) | 0.761848 ns [0.759412, 0.764693] (0.95 confidence) | — | — |
| public&#95;exact&#95;d3/vector&#95;from&#95;fn | 52.8798 ns [52.7761, 52.9877] (0.95 confidence) | 51.7165 ns [51.6551, 51.7931] (0.95 confidence) | — | — |
| public&#95;exact&#95;d3/vector&#95;into&#95;array | 23.0064 ns [22.9529, 23.054] (0.95 confidence) | 22.3356 ns [22.2877, 22.3746] (0.95 confidence) | — | — |
| public&#95;exact&#95;d3/vector&#95;new | 52.547 ns [52.4434, 52.642] (0.95 confidence) | 51.398 ns [51.325, 51.5158] (0.95 confidence) | — | — |
| public&#95;exact&#95;d3/vector&#95;zero | 2.05819 ns [2.05596, 2.06139] (0.95 confidence) | 2.01699 ns [2.01497, 2.01994] (0.95 confidence) | — | — |
| public&#95;exact&#95;d4/matrix&#95;access | 0.782333 ns [0.779796, 0.784365] (0.95 confidence) | 0.758367 ns [0.755138, 0.760587] (0.95 confidence) | — | — |
| public&#95;exact&#95;d4/matrix&#95;from&#95;fn | 291.499 ns [291.015, 292.22] (0.95 confidence) | 283.814 ns [283.428, 284.213] (0.95 confidence) | — | — |
| public&#95;exact&#95;d4/matrix&#95;from&#95;rows | 352.211 ns [351.591, 353.002] (0.95 confidence) | 344.167 ns [343.652, 344.493] (0.95 confidence) | — | — |
| public&#95;exact&#95;d4/matrix&#95;into&#95;rows | 133.509 ns [133.384, 133.763] (0.95 confidence) | 129.673 ns [129.369, 129.918] (0.95 confidence) | — | — |
| public&#95;exact&#95;d4/matrix&#95;zero | 15.5436 ns [15.4905, 15.674] (0.95 confidence) | 15.4061 ns [15.2825, 15.6589] (0.95 confidence) | — | — |
| public&#95;exact&#95;d4/vector&#95;access | 0.782022 ns [0.779766, 0.784321] (0.95 confidence) | 0.758791 ns [0.756932, 0.760432] (0.95 confidence) | — | — |
| public&#95;exact&#95;d4/vector&#95;from&#95;fn | 67.5261 ns [67.4714, 67.7017] (0.95 confidence) | 65.8041 ns [65.7334, 65.8556] (0.95 confidence) | — | — |
| public&#95;exact&#95;d4/vector&#95;into&#95;array | 30.481 ns [30.4196, 30.5202] (0.95 confidence) | 29.5823 ns [29.5279, 29.6265] (0.95 confidence) | — | — |
| public&#95;exact&#95;d4/vector&#95;new | 66.7524 ns [66.6796, 66.8537] (0.95 confidence) | 65.0457 ns [65.0075, 65.0837] (0.95 confidence) | — | — |
| public&#95;exact&#95;d4/vector&#95;zero | 3.47711 ns [3.44192, 3.52718] (0.95 confidence) | 3.48719 ns [3.45224, 3.50349] (0.95 confidence) | — | — |
| public&#95;exact&#95;d5/matrix&#95;access | 0.779626 ns [0.777094, 0.782142] (0.95 confidence) | 0.759472 ns [0.757495, 0.761201] (0.95 confidence) | — | — |
| public&#95;exact&#95;d5/matrix&#95;from&#95;fn | 440.428 ns [435.45, 443.028] (0.95 confidence) | 437.464 ns [436.456, 438.985] (0.95 confidence) | — | — |
| public&#95;exact&#95;d5/matrix&#95;from&#95;rows | 634.235 ns [633.452, 635.543] (0.95 confidence) | 621.57 ns [620.777, 622.853] (0.95 confidence) | — | — |
| public&#95;exact&#95;d5/matrix&#95;into&#95;rows | 308.987 ns [307.874, 310.76] (0.95 confidence) | 300.476 ns [299.454, 302.465] (0.95 confidence) | — | — |
| public&#95;exact&#95;d5/matrix&#95;zero | 22.2349 ns [22.1563, 22.2926] (0.95 confidence) | 21.789 ns [21.7099, 21.8579] (0.95 confidence) | — | — |
| public&#95;exact&#95;d5/vector&#95;access | 0.782101 ns [0.779456, 0.783233] (0.95 confidence) | 0.761072 ns [0.759043, 0.763089] (0.95 confidence) | — | — |
| public&#95;exact&#95;d5/vector&#95;from&#95;fn | 83.686 ns [83.5321, 83.8326] (0.95 confidence) | 81.7078 ns [81.5155, 81.9331] (0.95 confidence) | — | — |
| public&#95;exact&#95;d5/vector&#95;into&#95;array | 38.8168 ns [38.716, 38.8851] (0.95 confidence) | 37.9701 ns [37.8721, 38.1891] (0.95 confidence) | — | — |
| public&#95;exact&#95;d5/vector&#95;new | 83.6116 ns [83.5086, 83.7174] (0.95 confidence) | 81.9602 ns [81.7457, 82.1717] (0.95 confidence) | — | — |
| public&#95;exact&#95;d5/vector&#95;zero | 4.26042 ns [4.23619, 4.29255] (0.95 confidence) | 4.18914 ns [4.13479, 4.23614] (0.95 confidence) | — | — |
| public&#95;exact&#95;scalar/rounded | 5.60579 ns [5.59882, 5.61421] (0.95 confidence) | 5.55788 ns [5.5554, 5.55858] (0.95 confidence) | — | — |
| public&#95;exact&#95;scalar/strict | 4.26417 ns [4.25523, 4.26877] (0.95 confidence) | 4.18783 ns [4.18003, 4.19987] (0.95 confidence) | — | — |
| public&#95;scalar/interval&#95;access | 1.58163 ns [1.5803, 1.58525] (0.95 confidence) | 1.57685 ns [1.57546, 1.57924] (0.95 confidence) | — | — |
| public&#95;scalar/interval&#95;negate | 0.730185 ns [0.729651, 0.730984] (0.95 confidence) | 0.725859 ns [0.724932, 0.726396] (0.95 confidence) | — | — |
| public&#95;scalar/interval&#95;new | 0.893286 ns [0.892219, 0.893886] (0.95 confidence) | 0.89168 ns [0.890555, 0.892681] (0.95 confidence) | — | — |
| public&#95;scalar/interval&#95;point | 0.580283 ns [0.579782, 0.581068] (0.95 confidence) | 0.577213 ns [0.576399, 0.578628] (0.95 confidence) | — | — |
| public&#95;scalar/interval&#95;subtraction | 1.85835 ns [1.85672, 1.85946] (0.95 confidence) | 1.85109 ns [1.84874, 1.85307] (0.95 confidence) | — | — |
| public&#95;scalar/tolerance | 0.974368 ns [0.972578, 0.975068] (0.95 confidence) | 0.967521 ns [0.965368, 0.968941] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d14/prepared | — | 1762.79 ns [1760.75, 1769.04] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d14/slice | — | 629.731 ns [628.927, 630.242] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d2/prepared | — | 256.674 ns [256.34, 256.98] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d2/slice | — | 95.98 ns [95.8994, 96.1598] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d22/prepared | — | 2775.77 ns [2763.8, 2790.25] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d22/slice | — | 976.636 ns [975.615, 977.444] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d3/prepared | — | 380.418 ns [379.83, 380.834] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d3/slice | — | 125.337 ns [125.217, 125.703] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d4/prepared | — | 505.812 ns [504.819, 507.39] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d4/slice | — | 184.332 ns [184.169, 184.797] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d5/prepared | — | 625.923 ns [621.104, 630.612] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d5/slice | — | 211.11 ns [210.334, 211.803] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d6/prepared | — | 784.346 ns [762.099, 790.461] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d6/slice | — | 271.629 ns [270.499, 272.554] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d7/prepared | — | 866.572 ns [865.211, 870.445] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d7/slice | — | 298.356 ns [297.889, 298.977] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d8/prepared | — | 1006.29 ns [1004.04, 1009.8] (0.95 confidence) | — | — |
| rational&#95;dot/cancellation/d8/slice | — | 355.545 ns [355.049, 356.416] (0.95 confidence) | — | — |
| rational&#95;dot/large/d14/prepared | — | 2840.04 ns [2834.05, 2843.34] (0.95 confidence) | — | — |
| rational&#95;dot/large/d14/slice | — | 1245.07 ns [1238.78, 1250.86] (0.95 confidence) | — | — |
| rational&#95;dot/large/d2/prepared | — | 333.628 ns [332.773, 334.394] (0.95 confidence) | — | — |
| rational&#95;dot/large/d2/slice | — | 121.779 ns [121.39, 122.282] (0.95 confidence) | — | — |
| rational&#95;dot/large/d22/prepared | — | 4533.48 ns [4527.04, 4542.62] (0.95 confidence) | — | — |
| rational&#95;dot/large/d22/slice | — | 2012.84 ns [2008.85, 2022.74] (0.95 confidence) | — | — |
| rational&#95;dot/large/d3/prepared | — | 533.984 ns [532.825, 535.517] (0.95 confidence) | — | — |
| rational&#95;dot/large/d3/slice | — | 208.4 ns [207.961, 208.828] (0.95 confidence) | — | — |
| rational&#95;dot/large/d4/prepared | — | 716.192 ns [714.957, 716.658] (0.95 confidence) | — | — |
| rational&#95;dot/large/d4/slice | — | 292.6 ns [289.947, 293.889] (0.95 confidence) | — | — |
| rational&#95;dot/large/d5/prepared | — | 913.582 ns [911.793, 915.883] (0.95 confidence) | — | — |
| rational&#95;dot/large/d5/slice | — | 371.385 ns [368.996, 373.128] (0.95 confidence) | — | — |
| rational&#95;dot/large/d6/prepared | — | 1117.32 ns [1115.82, 1119.62] (0.95 confidence) | — | — |
| rational&#95;dot/large/d6/slice | — | 465.641 ns [464.977, 467.113] (0.95 confidence) | — | — |
| rational&#95;dot/large/d7/prepared | — | 1328.36 ns [1325.33, 1330.71] (0.95 confidence) | — | — |
| rational&#95;dot/large/d7/slice | — | 568.497 ns [567.429, 569.48] (0.95 confidence) | — | — |
| rational&#95;dot/large/d8/prepared | — | 1522.79 ns [1518.67, 1525.93] (0.95 confidence) | — | — |
| rational&#95;dot/large/d8/slice | — | 643.256 ns [641.807, 644.44] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d14/prepared | — | 2378.18 ns [2375.09, 2381.22] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d14/slice | — | 1026.12 ns [1023.23, 1028.64] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d2/prepared | — | 270.685 ns [270.261, 270.997] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d2/slice | — | 92.2164 ns [92.0467, 92.4231] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d22/prepared | — | 3809.77 ns [3805.59, 3812.84] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d22/slice | — | 1687.01 ns [1675.17, 1692.46] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d3/prepared | — | 438.269 ns [437.825, 439.15] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d3/slice | — | 168.463 ns [168.001, 169.062] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d4/prepared | — | 588.849 ns [588.494, 589.714] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d4/slice | — | 226.3 ns [225.526, 226.83] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d5/prepared | — | 757.6 ns [756.805, 759.146] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d5/slice | — | 300.221 ns [299.508, 300.776] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d6/prepared | — | 924.26 ns [922.982, 925.751] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d6/slice | — | 376.53 ns [375.862, 377.967] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d7/prepared | — | 1109.61 ns [1107.53, 1111.91] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d7/slice | — | 475.063 ns [474.112, 475.999] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d8/prepared | — | 1260.19 ns [1258.33, 1262.08] (0.95 confidence) | — | — |
| rational&#95;dot/non&#95;dyadic/d8/slice | — | 529.913 ns [528.557, 531.384] (0.95 confidence) | — | — |
| rational&#95;input&#95;d2/det&#95;row&#95;cleared&#95;bareiss | 178.923 ns [178.422, 179.482] (0.95 confidence) | 178.758 ns [178.419, 178.944] (0.95 confidence) | — | — |
| rational&#95;input&#95;d2/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 77.1264 ns [76.9536, 77.3438] (0.95 confidence) | 76.7382 ns [76.5073, 76.8367] (0.95 confidence) | — | — |
| rational&#95;input&#95;d2/solve&#95;row&#95;cleared&#95;bareiss | 615.301 ns [614.705, 615.826] (0.95 confidence) | 623.22 ns [622.759, 624.253] (0.95 confidence) | — | — |
| rational&#95;input&#95;d3/det&#95;row&#95;cleared&#95;bareiss | 569.182 ns [567.574, 569.919] (0.95 confidence) | 590.543 ns [589.085, 592.147] (0.95 confidence) | — | — |
| rational&#95;input&#95;d3/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 355.7 ns [353.625, 357.161] (0.95 confidence) | 389.335 ns [388.75, 390.638] (0.95 confidence) | — | — |
| rational&#95;input&#95;d3/solve&#95;row&#95;cleared&#95;bareiss | 1229.61 ns [1228.42, 1231.17] (0.95 confidence) | 1238 ns [1236.54, 1239.38] (0.95 confidence) | — | — |
| rational&#95;input&#95;d4/det&#95;row&#95;cleared&#95;bareiss | 1160.28 ns [1158.7, 1161.85] (0.95 confidence) | 1220.5 ns [1219.1, 1221.75] (0.95 confidence) | — | — |
| rational&#95;input&#95;d4/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 761.552 ns [760.623, 763.382] (0.95 confidence) | 829.312 ns [826.847, 831.528] (0.95 confidence) | — | — |
| rational&#95;input&#95;d4/solve&#95;row&#95;cleared&#95;bareiss | 2201.27 ns [2197.95, 2207.22] (0.95 confidence) | 2256.25 ns [2240.86, 2446.2] (0.95 confidence) | — | — |
| rational&#95;input&#95;d5/det&#95;row&#95;cleared&#95;bareiss | 2051.13 ns [2045.93, 2055.6] (0.95 confidence) | 2068.44 ns [2061.32, 2071.29] (0.95 confidence) | — | — |
| rational&#95;input&#95;d5/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 1535.72 ns [1528.67, 1539.47] (0.95 confidence) | 1564.41 ns [1562.18, 1567.45] (0.95 confidence) | — | — |
| rational&#95;input&#95;d5/solve&#95;row&#95;cleared&#95;bareiss | 3338.71 ns [3335.24, 3343] (0.95 confidence) | 3340.11 ns [3335.52, 3347] (0.95 confidence) | — | — |
| rational&#95;input&#95;d6/det&#95;row&#95;cleared&#95;bareiss | 3509.11 ns [3501.83, 3516.94] (0.95 confidence) | 3463.1 ns [3457.23, 3468.44] (0.95 confidence) | — | — |
| rational&#95;input&#95;d6/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 2838.61 ns [2834, 2842.66] (0.95 confidence) | 2847.15 ns [2841.38, 2853.69] (0.95 confidence) | — | — |
| rational&#95;input&#95;d6/solve&#95;row&#95;cleared&#95;bareiss | 5753.87 ns [5744.02, 5765.17] (0.95 confidence) | 5772.1 ns [5752.47, 5782.24] (0.95 confidence) | — | — |
| rational&#95;input&#95;d7/det&#95;row&#95;cleared&#95;bareiss | 5527 ns [5519.96, 5533.26] (0.95 confidence) | 5558.85 ns [5550.64, 5576.89] (0.95 confidence) | — | — |
| rational&#95;input&#95;d7/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 4655.7 ns [4650.08, 4663.42] (0.95 confidence) | 4675.89 ns [4667.92, 4680.16] (0.95 confidence) | — | — |
| rational&#95;input&#95;d7/solve&#95;row&#95;cleared&#95;bareiss | 9978.45 ns [9954.61, 9997.16] (0.95 confidence) | 9920.65 ns [9892.13, 9948.21] (0.95 confidence) | — | — |
| rational&#95;input&#95;d8/det&#95;row&#95;cleared&#95;bareiss | 8012.69 ns [8004.51, 8025.9] (0.95 confidence) | 8086.5 ns [8078.96, 8102.83] (0.95 confidence) | — | — |
| rational&#95;input&#95;d8/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 7166.59 ns [7155.52, 7180.28] (0.95 confidence) | 7287.49 ns [7282.65, 7305.28] (0.95 confidence) | — | — |
| rational&#95;input&#95;d8/solve&#95;row&#95;cleared&#95;bareiss | 15185.7 ns [15158.7, 15204.2] (0.95 confidence) | 14928 ns [14906.2, 14965.2] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d2/det&#95;row&#95;cleared&#95;bareiss | 42267 ns [42202.9, 42346.1] (0.95 confidence) | 41922.2 ns [41888.2, 41990.5] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d2/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 951.63 ns [950.688, 952.166] (0.95 confidence) | 946.652 ns [945.109, 948.833] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d2/solve&#95;row&#95;cleared&#95;bareiss | 3051.72 ns [3049, 3053.88] (0.95 confidence) | 3057.01 ns [3054.61, 3060.27] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d3/det&#95;row&#95;cleared&#95;bareiss | 69939.3 ns [69887.7, 70012.3] (0.95 confidence) | 68925.9 ns [68866.8, 69030] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d3/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 4208.2 ns [4201, 4213.22] (0.95 confidence) | 4183.77 ns [4174.63, 4188.78] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d3/solve&#95;row&#95;cleared&#95;bareiss | 11370.3 ns [11356.8, 11385.8] (0.95 confidence) | 11363.9 ns [11334, 11385.6] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d4/det&#95;row&#95;cleared&#95;bareiss | 110967 ns [110820, 111054] (0.95 confidence) | 110103 ns [109988, 110255] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d4/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 9436.04 ns [9426.5, 9448.55] (0.95 confidence) | 9344.38 ns [9334.52, 9359.37] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d4/solve&#95;row&#95;cleared&#95;bareiss | 33642.9 ns [33611, 33685.8] (0.95 confidence) | 33556.6 ns [33485.3, 33608.8] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d5/det&#95;row&#95;cleared&#95;bareiss | 186816 ns [186720, 187068] (0.95 confidence) | 186139 ns [185987, 186399] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d5/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 43051.9 ns [42971.6, 43127] (0.95 confidence) | 42449.3 ns [42409.5, 42493.1] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d5/solve&#95;row&#95;cleared&#95;bareiss | 69412.6 ns [69351.7, 69483.1] (0.95 confidence) | 68030.2 ns [67892, 68155.1] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d6/det&#95;row&#95;cleared&#95;bareiss | 302005 ns [301872, 302366] (0.95 confidence) | 300064 ns [299722, 300312] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d6/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 104280 ns [104155, 104375] (0.95 confidence) | 102755 ns [102562, 102872] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d6/solve&#95;row&#95;cleared&#95;bareiss | 161358 ns [161086, 161507] (0.95 confidence) | 158996 ns [158755, 159213] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d7/det&#95;row&#95;cleared&#95;bareiss | 503208 ns [502913, 503422] (0.95 confidence) | 497239 ns [496569, 498132] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d7/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 233505 ns [233132, 233725] (0.95 confidence) | 228622 ns [228378, 228842] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d7/solve&#95;row&#95;cleared&#95;bareiss | 339805 ns [339559, 339945] (0.95 confidence) | 332411 ns [331486, 332784] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d8/det&#95;row&#95;cleared&#95;bareiss | 822180 ns [821795, 822385] (0.95 confidence) | 812397 ns [811791, 813040] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d8/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 467967 ns [467642, 468742] (0.95 confidence) | 459684 ns [458853, 460378] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide1024&#95;d8/solve&#95;row&#95;cleared&#95;bareiss | 649629 ns [648958, 650037] (0.95 confidence) | 636835 ns [635428, 637639] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d2/det&#95;row&#95;cleared&#95;bareiss | 9015.5 ns [8967.05, 9062.49] (0.95 confidence) | 8814.2 ns [8784.21, 8846.61] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d2/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 651.941 ns [649.711, 654.246] (0.95 confidence) | 628.944 ns [628.011, 629.863] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d2/solve&#95;row&#95;cleared&#95;bareiss | 2079.31 ns [2075.39, 2083.36] (0.95 confidence) | 2069.67 ns [2066.18, 2074.58] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d3/det&#95;row&#95;cleared&#95;bareiss | 13596.2 ns [13561.5, 13640] (0.95 confidence) | 13192.2 ns [13164.4, 13226.4] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d3/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 2790.11 ns [2778.39, 2797.41] (0.95 confidence) | 2706.28 ns [2700.52, 2719.65] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d3/solve&#95;row&#95;cleared&#95;bareiss | 5753.47 ns [5731.98, 5759.75] (0.95 confidence) | 5712.4 ns [5706.59, 5716.7] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d4/det&#95;row&#95;cleared&#95;bareiss | 22305.1 ns [22261.7, 22372] (0.95 confidence) | 21660.4 ns [21618.9, 21736.1] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d4/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 5947.46 ns [5921.12, 5985.78] (0.95 confidence) | 5850.43 ns [5807.3, 5919.89] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d4/solve&#95;row&#95;cleared&#95;bareiss | 13175.8 ns [13109.9, 13269.7] (0.95 confidence) | 12571.6 ns [12525.1, 12602.9] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d5/det&#95;row&#95;cleared&#95;bareiss | 33430.9 ns [33397, 33486.9] (0.95 confidence) | 33413.1 ns [33269.6, 33482.8] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d5/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 13288.9 ns [13248.2, 13345.1] (0.95 confidence) | 13055.8 ns [13033, 13074] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d5/solve&#95;row&#95;cleared&#95;bareiss | 20590.5 ns [20545.9, 20670.7] (0.95 confidence) | 20369.9 ns [20341, 20397.7] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d6/det&#95;row&#95;cleared&#95;bareiss | 51904.6 ns [51832.1, 52192.9] (0.95 confidence) | 51631.2 ns [51501.5, 51695.3] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d6/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 26370.3 ns [26314.5, 26457.7] (0.95 confidence) | 25971.4 ns [25915.1, 26011.2] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d6/solve&#95;row&#95;cleared&#95;bareiss | 39678.3 ns [39567.7, 39748.7] (0.95 confidence) | 38452.3 ns [38386.3, 38611.5] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d7/det&#95;row&#95;cleared&#95;bareiss | 83257 ns [83175.2, 83372.7] (0.95 confidence) | 82111.7 ns [82049.3, 82166.6] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d7/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 50573.1 ns [50518.2, 50614.6] (0.95 confidence) | 49590.7 ns [49252.9, 49767.8] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d7/solve&#95;row&#95;cleared&#95;bareiss | 72181.4 ns [72143.7, 72241.1] (0.95 confidence) | 70519.6 ns [70388.1, 70595.7] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d8/det&#95;row&#95;cleared&#95;bareiss | 132001 ns [131825, 132152] (0.95 confidence) | 129874 ns [129642, 130041] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d8/det&#95;sign&#95;row&#95;cleared&#95;bareiss | 91290.7 ns [91078.4, 91375.3] (0.95 confidence) | 88555 ns [88421.9, 88686.5] (0.95 confidence) | — | — |
| rational&#95;input&#95;wide256&#95;d8/solve&#95;row&#95;cleared&#95;bareiss | 126767 ns [126731, 126848] (0.95 confidence) | 125277 ns [125092, 125514] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d2/adapter | 1021.89 ns [1014.81, 1069.25] (0.95 confidence) | 1045.7 ns [1031.57, 1078.14] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d2/construct | 289.79 ns [288.969, 290.52] (0.95 confidence) | 286.421 ns [285.797, 286.899] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d2/prepared | 637.701 ns [636.63, 638.858] (0.95 confidence) | 647.691 ns [645.86, 649.291] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d3/adapter | 1971.72 ns [1965.64, 1976.72] (0.95 confidence) | 1889.31 ns [1887.4, 1892.15] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d3/construct | 641.362 ns [639.894, 642.359] (0.95 confidence) | 632.271 ns [631.29, 634.266] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d3/prepared | 1174.82 ns [1171.98, 1176.45] (0.95 confidence) | 1146.84 ns [1142.64, 1148.78] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d4/adapter | 3386.83 ns [3377.59, 3400.51] (0.95 confidence) | 3319.53 ns [3312.69, 3325.11] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d4/construct | 1055.62 ns [1052.11, 1065.67] (0.95 confidence) | 1040.05 ns [1036.93, 1041.73] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d4/prepared | 2053.28 ns [2048.29, 2058.99] (0.95 confidence) | 2082.62 ns [2071.67, 2087.45] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d5/adapter | 5564.94 ns [5549.58, 5572.66] (0.95 confidence) | 5399.36 ns [5392, 5411.89] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d5/construct | 1956.13 ns [1951.25, 1963.88] (0.95 confidence) | 1919.21 ns [1912.03, 1922.92] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d5/prepared | 3295.74 ns [3286.49, 3304.29] (0.95 confidence) | 3130.63 ns [3123.72, 3137.39] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d6/adapter | 7762.99 ns [7756.36, 7777.01] (0.95 confidence) | 7541.46 ns [7536.42, 7550.09] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d6/construct | 2544.52 ns [2541.11, 2551.58] (0.95 confidence) | 2489.67 ns [2484.34, 2492.95] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d6/prepared | 4720.18 ns [4710.14, 4727.41] (0.95 confidence) | 4645.34 ns [4637.66, 4650.89] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d7/adapter | 10668.9 ns [10651.1, 10679.2] (0.95 confidence) | 10245.8 ns [10216.5, 10267.4] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d7/construct | 3490.83 ns [3475.87, 3500.22] (0.95 confidence) | 3424.18 ns [3420.38, 3429.89] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d7/prepared | 6533.56 ns [6519.63, 6561.91] (0.95 confidence) | 6435.59 ns [6413.28, 6452.83] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d8/adapter | 13867.5 ns [13850.6, 13882.4] (0.95 confidence) | 13280.1 ns [13264.7, 13301.7] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d8/construct | 4103.89 ns [4096.87, 4110.28] (0.95 confidence) | 4003.1 ns [3996.73, 4021.65] (0.95 confidence) | — | — |
| rational&#95;solve/dense/dyadic&#95;d8/prepared | 9619.44 ns [9603.5, 9632.91] (0.95 confidence) | 8969.29 ns [8943.95, 8998.03] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d2/adapter | 1043.36 ns [1032.73, 1059.15] (0.95 confidence) | 1110.16 ns [1056.83, 1132.96] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d2/construct | 330.933 ns [329.948, 331.94] (0.95 confidence) | 324.449 ns [323.517, 325.211] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d2/prepared | 615.289 ns [614.258, 616.1] (0.95 confidence) | 623.275 ns [622.058, 624.632] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d3/adapter | 2095.29 ns [2094.14, 2097.97] (0.95 confidence) | 2076.69 ns [2073.52, 2079.22] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d3/construct | 684.629 ns [681.081, 687.018] (0.95 confidence) | 661.009 ns [657.749, 665.676] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d3/prepared | 1278.33 ns [1276.12, 1279.93] (0.95 confidence) | 1276.82 ns [1274.89, 1278.77] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d4/adapter | 3497.42 ns [3493.34, 3503.5] (0.95 confidence) | 3474.5 ns [3468.72, 3477.38] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d4/construct | 1150.67 ns [1148.1, 1154.77] (0.95 confidence) | 1122.71 ns [1117.25, 1126.74] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d4/prepared | 2161.14 ns [2158.71, 2163.98] (0.95 confidence) | 2153.75 ns [2150.3, 2157.41] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d5/adapter | 5623.72 ns [5616.89, 5633.08] (0.95 confidence) | 5517.11 ns [5507.69, 5525.79] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d5/construct | 1934.11 ns [1929.5, 1936.52] (0.95 confidence) | 1893.54 ns [1885.66, 1899.17] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d5/prepared | 3438.83 ns [3433.67, 3443.01] (0.95 confidence) | 3322.38 ns [3318.68, 3326.04] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d6/adapter | 7987.12 ns [7975.16, 8000.27] (0.95 confidence) | 7851.1 ns [7839.75, 7864.14] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d6/construct | 2790.69 ns [2781.81, 2806.75] (0.95 confidence) | 2700.56 ns [2690.91, 2712.94] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d6/prepared | 4925.71 ns [4913.76, 4934.49] (0.95 confidence) | 4813.25 ns [4804.2, 4818.61] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d7/adapter | 11771.4 ns [11753.5, 11792.3] (0.95 confidence) | 11592.7 ns [11535.1, 11650.3] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d7/construct | 3944.11 ns [3935.07, 3961.32] (0.95 confidence) | 3863.59 ns [3837.33, 3893.22] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d7/prepared | 7397.66 ns [7387.98, 7411.11] (0.95 confidence) | 7370.48 ns [7338.9, 7380.61] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d8/adapter | 16370.4 ns [16347.1, 16405.6] (0.95 confidence) | 15848.7 ns [15808.8, 15884.3] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d8/construct | 5253.24 ns [5219.9, 5271.32] (0.95 confidence) | 5070.7 ns [5014.78, 5097.65] (0.95 confidence) | — | — |
| rational&#95;solve/dense/rational&#95;d8/prepared | 10539.8 ns [10519.2, 10559.8] (0.95 confidence) | 10435.9 ns [10397.9, 10475.1] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d2/adapter | 437.22 ns [433.734, 478.532] (0.95 confidence) | 425.406 ns [424.838, 425.989] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d2/construct | 227.308 ns [227.114, 227.6] (0.95 confidence) | 223.366 ns [222.645, 223.791] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d2/prepared | 152.776 ns [152.348, 153.166] (0.95 confidence) | 151.425 ns [150.582, 151.734] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d3/adapter | 737.061 ns [736.005, 738.42] (0.95 confidence) | 727.744 ns [726.164, 728.497] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d3/construct | 412.856 ns [412.206, 414.018] (0.95 confidence) | 404.214 ns [403.532, 405.228] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d3/prepared | 250.306 ns [249.918, 250.699] (0.95 confidence) | 247.255 ns [246.488, 247.81] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d4/adapter | 1032.14 ns [1031.07, 1034.33] (0.95 confidence) | 1011.38 ns [1009.96, 1013.42] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d4/construct | 611.593 ns [610.674, 613.104] (0.95 confidence) | 594.422 ns [592.92, 596.307] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d4/prepared | 323.072 ns [322.639, 323.57] (0.95 confidence) | 321.782 ns [321.308, 322.177] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d5/adapter | 1386.89 ns [1383.47, 1389.89] (0.95 confidence) | 1353.62 ns [1348.58, 1359.52] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d5/construct | 848.353 ns [844.888, 851.225] (0.95 confidence) | 833.876 ns [831.515, 836.534] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d5/prepared | 425.064 ns [424.162, 425.572] (0.95 confidence) | 415.933 ns [414.652, 416.876] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d6/adapter | 1760.58 ns [1757.9, 1763.6] (0.95 confidence) | 1729.44 ns [1725.85, 1732.39] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d6/construct | 1114.04 ns [1112.72, 1115.76] (0.95 confidence) | 1097.13 ns [1095.1, 1098.84] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d6/prepared | 515.006 ns [513.863, 515.756] (0.95 confidence) | 499.341 ns [498.514, 499.98] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d7/adapter | 2200.41 ns [2194.34, 2207.68] (0.95 confidence) | 2152.46 ns [2149.16, 2158.71] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d7/construct | 1457.19 ns [1452.65, 1460.96] (0.95 confidence) | 1402.32 ns [1401.27, 1405.17] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d7/prepared | 602.406 ns [601.991, 603.85] (0.95 confidence) | 597.304 ns [596.048, 598.205] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d8/adapter | 2525.31 ns [2510.13, 2533.7] (0.95 confidence) | 2463.81 ns [2459.77, 2467.8] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d8/construct | 1747.95 ns [1729.35, 1755.66] (0.95 confidence) | 1680.48 ns [1679.01, 1683.5] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/dyadic&#95;d8/prepared | 647.223 ns [646.137, 649.111] (0.95 confidence) | 635.292 ns [634.586, 636.428] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d2/adapter | 583.359 ns [579.584, 623.946] (0.95 confidence) | 561.361 ns [560.579, 562.815] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d2/construct | 309.97 ns [308.268, 311.051] (0.95 confidence) | 305.407 ns [303.923, 306.417] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d2/prepared | 189.491 ns [189.206, 189.869] (0.95 confidence) | 188.665 ns [188.377, 188.926] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d3/adapter | 922.525 ns [920.481, 924.321] (0.95 confidence) | 900.667 ns [897.58, 902.26] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d3/construct | 539.931 ns [538.173, 541.408] (0.95 confidence) | 524.223 ns [522.174, 527.062] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d3/prepared | 282.072 ns [281.75, 282.454] (0.95 confidence) | 278.716 ns [277.836, 279.382] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d4/adapter | 1291.29 ns [1289.29, 1292.79] (0.95 confidence) | 1246.5 ns [1244.01, 1249.95] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d4/construct | 791.143 ns [788.197, 795.289] (0.95 confidence) | 754.668 ns [747.711, 758.473] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d4/prepared | 367.341 ns [366.886, 367.884] (0.95 confidence) | 367.05 ns [366.284, 367.727] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d5/adapter | 1713.07 ns [1710.98, 1718.46] (0.95 confidence) | 1660.49 ns [1657.95, 1662.74] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d5/construct | 1078.28 ns [1071.82, 1083.11] (0.95 confidence) | 1067.1 ns [1062.54, 1072.18] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d5/prepared | 483.438 ns [482.684, 484.234] (0.95 confidence) | 477.524 ns [476.067, 478.571] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d6/adapter | 2162.07 ns [2156.37, 2168.16] (0.95 confidence) | 2067.48 ns [2062.08, 2070.81] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d6/construct | 1413.45 ns [1408.91, 1418.94] (0.95 confidence) | 1379.24 ns [1372.49, 1385.94] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d6/prepared | 576.035 ns [574.744, 577.343] (0.95 confidence) | 561.874 ns [560.696, 563.241] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d7/adapter | 2657.16 ns [2648.34, 2667.16] (0.95 confidence) | 2544.75 ns [2541.49, 2551.23] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d7/construct | 1806.43 ns [1802.86, 1812.64] (0.95 confidence) | 1753.57 ns [1746.28, 1766.37] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d7/prepared | 678.616 ns [676.612, 679.599] (0.95 confidence) | 666.053 ns [664.78, 666.929] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d8/adapter | 3071.34 ns [3063.86, 3078.33] (0.95 confidence) | 2976.31 ns [2969.36, 2981.67] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d8/construct | 2148.64 ns [2135.07, 2160.02] (0.95 confidence) | 2084.5 ns [2072.93, 2092.27] (0.95 confidence) | — | — |
| rational&#95;solve/diagonal/rational&#95;d8/prepared | 769.792 ns [768.983, 772.048] (0.95 confidence) | 755.467 ns [754.062, 757.14] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d2/adapter | 531.12 ns [523.818, 566.979] (0.95 confidence) | 508.777 ns [508.203, 509.497] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d2/construct | 243.393 ns [242.956, 243.891] (0.95 confidence) | 238.464 ns [238.253, 239.089] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d2/prepared | 217.311 ns [217.107, 217.668] (0.95 confidence) | 217.301 ns [216.991, 217.717] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d3/adapter | 1027.78 ns [1025.09, 1032.17] (0.95 confidence) | 1024.86 ns [1022.67, 1027.65] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d3/construct | 487.176 ns [485.463, 488.215] (0.95 confidence) | 476.473 ns [475.616, 477.712] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d3/prepared | 456.656 ns [456.062, 457.101] (0.95 confidence) | 450.801 ns [450.135, 451.691] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d4/adapter | 1612.35 ns [1611.11, 1613.42] (0.95 confidence) | 1588.37 ns [1586.38, 1590.36] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d4/construct | 824.84 ns [823.197, 828.043] (0.95 confidence) | 798.981 ns [798.051, 800.577] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d4/prepared | 644.892 ns [644.157, 645.715] (0.95 confidence) | 650.233 ns [649.47, 650.983] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d5/adapter | 2442.95 ns [2437.74, 2448.13] (0.95 confidence) | 2404.79 ns [2400.78, 2409.76] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d5/construct | 1246.31 ns [1242.41, 1252.32] (0.95 confidence) | 1220.87 ns [1217.88, 1222.87] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d5/prepared | 991.778 ns [988.987, 993.879] (0.95 confidence) | 967.726 ns [966.363, 968.55] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d6/adapter | 3309.69 ns [3292.57, 3321.34] (0.95 confidence) | 3237.53 ns [3233.39, 3241.46] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d6/construct | 1731.27 ns [1728.23, 1735.96] (0.95 confidence) | 1701.11 ns [1699.21, 1704.19] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d6/prepared | 1322.45 ns [1319.71, 1325.51] (0.95 confidence) | 1305.97 ns [1302.97, 1308.56] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d7/adapter | 4300.27 ns [4290.65, 4304.66] (0.95 confidence) | 4248.04 ns [4243.78, 4252.61] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d7/construct | 2303.88 ns [2300.75, 2308.16] (0.95 confidence) | 2269.26 ns [2265.93, 2275.69] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d7/prepared | 1742.59 ns [1739.59, 1744.4] (0.95 confidence) | 1725.44 ns [1722.61, 1728.25] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d8/adapter | 5317.89 ns [5294.63, 5327.24] (0.95 confidence) | 5167.92 ns [5157.87, 5184.16] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d8/construct | 2947.77 ns [2936.68, 2955.08] (0.95 confidence) | 2913.75 ns [2907.63, 2924.21] (0.95 confidence) | — | — |
| rational&#95;solve/lower/dyadic&#95;d8/prepared | 2140.73 ns [2137.76, 2144.78] (0.95 confidence) | 2102.13 ns [2100.35, 2106.21] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d2/adapter | 659.011 ns [657.272, 700.159] (0.95 confidence) | 633.722 ns [632.698, 634.275] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d2/construct | 349.448 ns [346.756, 350.241] (0.95 confidence) | 342.254 ns [340.137, 343.792] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d2/prepared | 216.338 ns [215.969, 216.591] (0.95 confidence) | 212.82 ns [212.298, 213.152] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d3/adapter | 1163.13 ns [1161.91, 1165.18] (0.95 confidence) | 1137.24 ns [1135.29, 1139.5] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d3/construct | 665.421 ns [661.089, 668.649] (0.95 confidence) | 653.393 ns [650.887, 656.089] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d3/prepared | 372.033 ns [371.421, 372.402] (0.95 confidence) | 363.694 ns [362.995, 364.781] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d4/adapter | 1836.83 ns [1835.32, 1838.73] (0.95 confidence) | 1781.32 ns [1779.65, 1783.71] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d4/construct | 1059.49 ns [1052.81, 1063.73] (0.95 confidence) | 1044.2 ns [1039.78, 1050.92] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d4/prepared | 589.58 ns [588.895, 590.153] (0.95 confidence) | 578.188 ns [577.383, 579.697] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d5/adapter | 2653.81 ns [2647.43, 2660.42] (0.95 confidence) | 2557.93 ns [2551.49, 2562.47] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d5/construct | 1601.38 ns [1590.59, 1614.11] (0.95 confidence) | 1569.85 ns [1558.5, 1576.75] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d5/prepared | 821.123 ns [819.77, 822.214] (0.95 confidence) | 804.259 ns [802.17, 805.751] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d6/adapter | 3606.23 ns [3599.37, 3612.76] (0.95 confidence) | 3493.95 ns [3489.52, 3500.09] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d6/construct | 2328.81 ns [2315.13, 2335.68] (0.95 confidence) | 2257.95 ns [2251.8, 2263.96] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d6/prepared | 1039.32 ns [1037.2, 1041.62] (0.95 confidence) | 1015.56 ns [1012.92, 1018.58] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d7/adapter | 4761.58 ns [4747.42, 4773.9] (0.95 confidence) | 4560.35 ns [4556.77, 4566.33] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d7/construct | 3067.48 ns [3049.06, 3081.04] (0.95 confidence) | 2973.08 ns [2961.89, 2983.04] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d7/prepared | 1393.14 ns [1389.6, 1396.65] (0.95 confidence) | 1355.46 ns [1352.57, 1359.78] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d8/adapter | 5914.4 ns [5902.08, 5940.02] (0.95 confidence) | 5703.18 ns [5690.64, 5713.52] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d8/construct | 3943.46 ns [3934.52, 3961.18] (0.95 confidence) | 3835.37 ns [3821.21, 3844.3] (0.95 confidence) | — | — |
| rational&#95;solve/lower/rational&#95;d8/prepared | 1642.6 ns [1639.12, 1645.06] (0.95 confidence) | 1608.26 ns [1605.16, 1614.97] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d2/adapter | 539.402 ns [526.663, 569.803] (0.95 confidence) | 507.984 ns [507.494, 508.998] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d2/construct | 244.662 ns [244.219, 245.206] (0.95 confidence) | 239.281 ns [239.001, 239.586] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d2/prepared | 217.707 ns [217.359, 218.001] (0.95 confidence) | 217.062 ns [216.831, 217.246] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d3/adapter | 1040.62 ns [1038.15, 1041.95] (0.95 confidence) | 1028.99 ns [1026.52, 1031.13] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d3/construct | 491.883 ns [489.768, 494.388] (0.95 confidence) | 482.383 ns [481.306, 483.661] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d3/prepared | 458.047 ns [457.476, 458.545] (0.95 confidence) | 447.108 ns [445.752, 448.08] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d4/adapter | 1618.17 ns [1613.58, 1621.1] (0.95 confidence) | 1594.09 ns [1590.82, 1595.7] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d4/construct | 832.769 ns [830.644, 833.969] (0.95 confidence) | 802.841 ns [798.827, 807.723] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d4/prepared | 646.393 ns [645.333, 648.247] (0.95 confidence) | 651.045 ns [650.455, 652.081] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d5/adapter | 2450 ns [2443.45, 2457.2] (0.95 confidence) | 2395.42 ns [2385.79, 2400.82] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d5/construct | 1248.55 ns [1242.73, 1254.94] (0.95 confidence) | 1214.38 ns [1212.5, 1215.75] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d5/prepared | 995.378 ns [993.843, 997.691] (0.95 confidence) | 965.219 ns [963.955, 966.284] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d6/adapter | 3292.79 ns [3286.62, 3301.89] (0.95 confidence) | 3256.18 ns [3249.77, 3266.89] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d6/construct | 1739.59 ns [1734.55, 1744.73] (0.95 confidence) | 1705.44 ns [1697.65, 1712.83] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d6/prepared | 1322.76 ns [1320.36, 1326.46] (0.95 confidence) | 1305.73 ns [1304.59, 1309.06] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d7/adapter | 4320.43 ns [4314.01, 4330.63] (0.95 confidence) | 4254.3 ns [4250.32, 4261.16] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d7/construct | 2302.49 ns [2297.17, 2309.92] (0.95 confidence) | 2243.46 ns [2240.25, 2248.52] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d7/prepared | 1745.67 ns [1742.4, 1748.68] (0.95 confidence) | 1724.86 ns [1720.82, 1727.27] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d8/adapter | 5311.16 ns [5300.65, 5328.57] (0.95 confidence) | 5201 ns [5195.74, 5207.21] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d8/construct | 2942.17 ns [2929.57, 2948.35] (0.95 confidence) | 2881.11 ns [2875.51, 2886.08] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/dyadic&#95;d8/prepared | 2136.61 ns [2130.72, 2138.95] (0.95 confidence) | 2126.84 ns [2123.06, 2131.02] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d2/adapter | 659.282 ns [655.085, 702.346] (0.95 confidence) | 632.716 ns [631.418, 633.692] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d2/construct | 351.958 ns [348.815, 353.686] (0.95 confidence) | 343.363 ns [342.163, 344.882] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d2/prepared | 216.157 ns [215.749, 216.364] (0.95 confidence) | 212.62 ns [212.268, 212.934] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d3/adapter | 1161.06 ns [1159.01, 1162.64] (0.95 confidence) | 1137.97 ns [1135.2, 1142.48] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d3/construct | 666.34 ns [664.283, 669.782] (0.95 confidence) | 652.546 ns [649.152, 656.133] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d3/prepared | 373.782 ns [372.965, 374.237] (0.95 confidence) | 361.557 ns [360.939, 362.028] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d4/adapter | 1828.49 ns [1826.49, 1830.88] (0.95 confidence) | 1778.05 ns [1774.93, 1779.69] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d4/construct | 1062.51 ns [1057.45, 1067.48] (0.95 confidence) | 1034.45 ns [1025.32, 1039.27] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d4/prepared | 587.466 ns [586.973, 588.225] (0.95 confidence) | 582.939 ns [580.571, 584.592] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d5/adapter | 2637.73 ns [2630.14, 2650.74] (0.95 confidence) | 2568.07 ns [2561.51, 2576.8] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d5/construct | 1620.63 ns [1605.24, 1626.27] (0.95 confidence) | 1573.46 ns [1562.58, 1581.36] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d5/prepared | 822.473 ns [820.481, 824.181] (0.95 confidence) | 802.73 ns [800.518, 804.11] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d6/adapter | 3621.67 ns [3616.1, 3627.93] (0.95 confidence) | 3490.81 ns [3488.1, 3496.02] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d6/construct | 2314.14 ns [2305.67, 2320.3] (0.95 confidence) | 2242.74 ns [2232.99, 2248.51] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d6/prepared | 1037.1 ns [1035.54, 1039.04] (0.95 confidence) | 1012.19 ns [1009.14, 1014.6] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d7/adapter | 4759.7 ns [4753.79, 4773.56] (0.95 confidence) | 4590.85 ns [4580.22, 4599.72] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d7/construct | 3040.38 ns [3023.75, 3068.78] (0.95 confidence) | 2990.16 ns [2983.38, 2998.16] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d7/prepared | 1388.8 ns [1385.52, 1393.37] (0.95 confidence) | 1356.02 ns [1353.81, 1358.21] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d8/adapter | 5942.42 ns [5927.63, 5953.68] (0.95 confidence) | 5743.57 ns [5735.18, 5750.77] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d8/construct | 3956.18 ns [3922.26, 3975.8] (0.95 confidence) | 3837.38 ns [3815.32, 3859.08] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;lower/rational&#95;d8/prepared | 1662.97 ns [1654.67, 1674.23] (0.95 confidence) | 1629.9 ns [1616.56, 1636.1] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d2/adapter | 536.257 ns [530.014, 575.187] (0.95 confidence) | 513.804 ns [512.674, 515.021] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d2/construct | 239.843 ns [239.267, 240.298] (0.95 confidence) | 233.283 ns [232.948, 233.545] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d2/prepared | 217.527 ns [217.036, 217.825] (0.95 confidence) | 217.834 ns [217.384, 218.281] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d3/adapter | 1018.09 ns [1016.88, 1019.47] (0.95 confidence) | 997.639 ns [995.659, 999.489] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d3/construct | 479.595 ns [477.32, 480.953] (0.95 confidence) | 468.947 ns [468.099, 469.897] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d3/prepared | 444.223 ns [443.339, 445.528] (0.95 confidence) | 434.164 ns [433.333, 434.862] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d4/adapter | 1647.84 ns [1644.87, 1652.76] (0.95 confidence) | 1613.18 ns [1607.54, 1616.83] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d4/construct | 824.285 ns [822.294, 825.418] (0.95 confidence) | 798.121 ns [796.411, 799.904] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d4/prepared | 674.72 ns [674.075, 675.325] (0.95 confidence) | 680.518 ns [679.699, 681.852] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d5/adapter | 2377.16 ns [2370.47, 2385.31] (0.95 confidence) | 2324.31 ns [2318.63, 2329.79] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d5/construct | 1204.7 ns [1202.59, 1207.31] (0.95 confidence) | 1179.49 ns [1177.42, 1181.89] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d5/prepared | 969.752 ns [968.445, 972.744] (0.95 confidence) | 942.489 ns [941.98, 943.455] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d6/adapter | 3369.96 ns [3362, 3375.61] (0.95 confidence) | 3314.52 ns [3307.93, 3322.51] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d6/construct | 1798.68 ns [1793.67, 1803.82] (0.95 confidence) | 1761.88 ns [1757.01, 1766.71] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d6/prepared | 1353.76 ns [1350.93, 1356.46] (0.95 confidence) | 1341.66 ns [1340.04, 1342.84] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d7/adapter | 4402.07 ns [4394.38, 4415.03] (0.95 confidence) | 4292.64 ns [4288.14, 4303.02] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d7/construct | 2367.67 ns [2360.87, 2378.01] (0.95 confidence) | 2311.68 ns [2303.36, 2320.38] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d7/prepared | 1765.82 ns [1761.93, 1770.6] (0.95 confidence) | 1750.15 ns [1747.58, 1753.72] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d8/adapter | 5496.7 ns [5473.47, 5515.07] (0.95 confidence) | 5368.7 ns [5355.86, 5387.32] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d8/construct | 3138.65 ns [3122.63, 3151.18] (0.95 confidence) | 3036.19 ns [3024.27, 3047.36] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/dyadic&#95;d8/prepared | 2105.29 ns [2103.1, 2108.23] (0.95 confidence) | 2085.29 ns [2082.25, 2090.79] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d2/adapter | 694.915 ns [691.558, 734.212] (0.95 confidence) | 667.577 ns [666.172, 668.336] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d2/construct | 377.923 ns [377.064, 379.483] (0.95 confidence) | 367.913 ns [366.243, 369.91] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d2/prepared | 229.27 ns [228.763, 229.863] (0.95 confidence) | 225.143 ns [224.726, 225.842] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d3/adapter | 1190.9 ns [1187.75, 1193.24] (0.95 confidence) | 1163.44 ns [1161.93, 1164.77] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d3/construct | 691.989 ns [690.371, 693.729] (0.95 confidence) | 680.248 ns [677.336, 682.901] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d3/prepared | 388.802 ns [387.721, 389.759] (0.95 confidence) | 382.166 ns [381.409, 382.809] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d4/adapter | 1785.94 ns [1776.9, 1795.47] (0.95 confidence) | 1737.64 ns [1735.4, 1740.17] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d4/construct | 1014.69 ns [1007.56, 1020.57] (0.95 confidence) | 980.256 ns [976.961, 983.437] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d4/prepared | 620.72 ns [620.042, 621.569] (0.95 confidence) | 613.184 ns [612.799, 613.917] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d5/adapter | 2536.31 ns [2529.13, 2544.71] (0.95 confidence) | 2455.91 ns [2450.88, 2459.51] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d5/construct | 1449.24 ns [1444.53, 1457.52] (0.95 confidence) | 1420.68 ns [1416.36, 1427.23] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d5/prepared | 887.266 ns [885.956, 888.639] (0.95 confidence) | 866.652 ns [864.506, 868.007] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d6/adapter | 3440.77 ns [3426.55, 3447.72] (0.95 confidence) | 3314.09 ns [3309.11, 3320.15] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d6/construct | 2036.38 ns [2026.33, 2045.66] (0.95 confidence) | 1977.63 ns [1966.07, 1983.63] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d6/prepared | 1147.21 ns [1145.04, 1149.56] (0.95 confidence) | 1117.12 ns [1115.04, 1119.09] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d7/adapter | 4564.39 ns [4556.85, 4578.58] (0.95 confidence) | 4363.67 ns [4356.25, 4368.24] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d7/construct | 2789.14 ns [2780.13, 2801.28] (0.95 confidence) | 2698.03 ns [2687.72, 2702.51] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d7/prepared | 1479.64 ns [1476.55, 1482.72] (0.95 confidence) | 1452.78 ns [1450.3, 1453.99] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d8/adapter | 5809.84 ns [5793.92, 5843.5] (0.95 confidence) | 5717.21 ns [5696.57, 5741.25] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d8/construct | 3584.64 ns [3568.93, 3621.26] (0.95 confidence) | 3502.99 ns [3491.75, 3511.54] (0.95 confidence) | — | — |
| rational&#95;solve/permuted&#95;upper/rational&#95;d8/prepared | 1890.18 ns [1888.33, 1894.19] (0.95 confidence) | 1851.25 ns [1847.87, 1855.57] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d2/adapter | 1030.04 ns [1002.02, 1064.96] (0.95 confidence) | 1030.39 ns [1024.54, 1060.94] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d2/construct | 293.791 ns [292.94, 294.651] (0.95 confidence) | 286.901 ns [286.037, 287.913] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d2/prepared | 617.015 ns [616.172, 617.813] (0.95 confidence) | 629.567 ns [628.322, 630.168] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d3/adapter | 1805.03 ns [1802.68, 1806.61] (0.95 confidence) | 1753.26 ns [1750.32, 1757.7] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d3/construct | 580.313 ns [577.304, 581.034] (0.95 confidence) | 571.457 ns [569.617, 573.145] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d3/prepared | 1092.46 ns [1088.69, 1100.38] (0.95 confidence) | 1069.38 ns [1068.47, 1071.13] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d4/adapter | 2715.19 ns [2710.35, 2718.46] (0.95 confidence) | 2676.13 ns [2669.26, 2680.23] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d4/construct | 808.343 ns [804.799, 812.819] (0.95 confidence) | 781.676 ns [777.474, 785.112] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d4/prepared | 1740.91 ns [1737.62, 1747.23] (0.95 confidence) | 1735.09 ns [1728.47, 1739.72] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d5/adapter | 4215 ns [4202.3, 4219.34] (0.95 confidence) | 4153.12 ns [4149.15, 4156.47] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d5/construct | 1323.05 ns [1317.21, 1331.71] (0.95 confidence) | 1309.06 ns [1306.35, 1314.07] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d5/prepared | 2674.47 ns [2662.98, 2683.64] (0.95 confidence) | 2607.19 ns [2602.8, 2610.33] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d6/adapter | 5638.41 ns [5631.77, 5646.73] (0.95 confidence) | 5555.36 ns [5548.44, 5563.5] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d6/construct | 1703.29 ns [1701.29, 1712] (0.95 confidence) | 1651.65 ns [1646.16, 1659.21] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d6/prepared | 3672.65 ns [3667.8, 3680.61] (0.95 confidence) | 3609.62 ns [3604.28, 3613.62] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d7/adapter | 7930.95 ns [7915.68, 7952.9] (0.95 confidence) | 7701.43 ns [7683.98, 7716.92] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d7/construct | 2223.37 ns [2214.59, 2238.31] (0.95 confidence) | 2140.71 ns [2136.72, 2144.82] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d7/prepared | 5324.83 ns [5316.92, 5335.38] (0.95 confidence) | 5190.13 ns [5181.1, 5199.55] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d8/adapter | 9576.44 ns [9562.86, 9582.89] (0.95 confidence) | 9403.17 ns [9387.65, 9419.44] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d8/construct | 2426.27 ns [2416.42, 2440.97] (0.95 confidence) | 2348.45 ns [2344.76, 2355.76] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/dyadic&#95;d8/prepared | 7041.27 ns [7030.03, 7047.89] (0.95 confidence) | 6858.79 ns [6850.26, 6870.09] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d2/adapter | 1169.33 ns [1160.32, 1179.14] (0.95 confidence) | 1170.66 ns [1164.02, 1226.27] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d2/construct | 479.576 ns [477.593, 481.592] (0.95 confidence) | 472.396 ns [471.469, 473.221] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d2/prepared | 570.8 ns [569.986, 571.395] (0.95 confidence) | 579.281 ns [578.433, 580.184] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d3/adapter | 1929.94 ns [1928.14, 1932.45] (0.95 confidence) | 1885.12 ns [1882.48, 1887.34] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d3/construct | 752.066 ns [749.329, 754.726] (0.95 confidence) | 728.489 ns [726.129, 730.001] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d3/prepared | 1026.18 ns [1025, 1028.07] (0.95 confidence) | 1019.44 ns [1017.93, 1020.54] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d4/adapter | 3163.65 ns [3159.31, 3166.27] (0.95 confidence) | 3067.91 ns [3063.59, 3071.61] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d4/construct | 1277.36 ns [1272.28, 1281.42] (0.95 confidence) | 1242.28 ns [1239.85, 1245.22] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d4/prepared | 1660.51 ns [1656.46, 1665.93] (0.95 confidence) | 1639.22 ns [1636.49, 1643.8] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d5/adapter | 4425.24 ns [4420.93, 4430.64] (0.95 confidence) | 4281.26 ns [4275.55, 4291.75] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d5/construct | 1671.72 ns [1667.36, 1675.91] (0.95 confidence) | 1632.85 ns [1628.65, 1636.03] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d5/prepared | 2480.29 ns [2476.76, 2485.12] (0.95 confidence) | 2431.97 ns [2429.02, 2434.35] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d6/adapter | 6194.18 ns [6170.5, 6217.82] (0.95 confidence) | 5827.35 ns [5813.81, 5842] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d6/construct | 2273.93 ns [2266.76, 2279.37] (0.95 confidence) | 2221.45 ns [2215.19, 2224.45] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d6/prepared | 3481.84 ns [3476.81, 3486.35] (0.95 confidence) | 3395.34 ns [3392.25, 3401.52] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d7/adapter | 8559.26 ns [8551.96, 8576.94] (0.95 confidence) | 8431.94 ns [8426.05, 8443.19] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d7/construct | 2915.96 ns [2906.03, 2929.8] (0.95 confidence) | 2835.9 ns [2826.43, 2852.06] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d7/prepared | 5395.3 ns [5383.08, 5415.65] (0.95 confidence) | 5364.65 ns [5356.5, 5374.09] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d8/adapter | 10600.6 ns [10573.9, 10612.7] (0.95 confidence) | 10325.8 ns [10315.3, 10336.5] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d8/construct | 3429.27 ns [3424.34, 3441.2] (0.95 confidence) | 3336.21 ns [3323.91, 3344.01] (0.95 confidence) | — | — |
| rational&#95;solve/sparse/rational&#95;d8/prepared | 6824.57 ns [6813.09, 6837.94] (0.95 confidence) | 6716.79 ns [6710.05, 6728.19] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d2/adapter | 534.623 ns [531.357, 573.504] (0.95 confidence) | 512.631 ns [511.88, 513.741] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d2/construct | 241.317 ns [241.053, 241.704] (0.95 confidence) | 235.943 ns [235.164, 236.688] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d2/prepared | 218.145 ns [218.034, 218.457] (0.95 confidence) | 217.966 ns [217.663, 218.284] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d3/adapter | 1011.58 ns [1010.42, 1012.5] (0.95 confidence) | 999.847 ns [997.11, 1002.65] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d3/construct | 475.695 ns [474.69, 477.696] (0.95 confidence) | 470.892 ns [470.093, 471.684] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d3/prepared | 443.964 ns [443.391, 444.58] (0.95 confidence) | 436.544 ns [435.869, 436.937] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d4/adapter | 1657.49 ns [1653.92, 1662.95] (0.95 confidence) | 1619.87 ns [1616.71, 1624.97] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d4/construct | 822.512 ns [819.75, 825.887] (0.95 confidence) | 803.333 ns [801.401, 805.662] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d4/prepared | 674.634 ns [674.019, 675.296] (0.95 confidence) | 680.528 ns [679.915, 681.061] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d5/adapter | 2355.28 ns [2350.54, 2359.32] (0.95 confidence) | 2313.34 ns [2310.06, 2315.11] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d5/construct | 1210.04 ns [1205.49, 1215.49] (0.95 confidence) | 1185.43 ns [1184.11, 1188.24] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d5/prepared | 963.068 ns [961.548, 964.172] (0.95 confidence) | 942.26 ns [940.633, 942.961] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d6/adapter | 3381.24 ns [3368.48, 3394.89] (0.95 confidence) | 3341.48 ns [3334.03, 3345.25] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d6/construct | 1804.5 ns [1799.87, 1811.24] (0.95 confidence) | 1800.87 ns [1782.23, 1810.32] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d6/prepared | 1358.13 ns [1355.24, 1362.27] (0.95 confidence) | 1347.27 ns [1344.4, 1351.29] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d7/adapter | 4371.86 ns [4364.79, 4385.02] (0.95 confidence) | 4299.26 ns [4295.05, 4303.62] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d7/construct | 2385.55 ns [2378.61, 2391.18] (0.95 confidence) | 2307.93 ns [2301.82, 2316.16] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d7/prepared | 1774.19 ns [1771.7, 1777.03] (0.95 confidence) | 1754.09 ns [1751.95, 1757.97] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d8/adapter | 5448.86 ns [5436.93, 5461.96] (0.95 confidence) | 5409.56 ns [5389.23, 5435.19] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d8/construct | 3157.75 ns [3148.54, 3164.39] (0.95 confidence) | 3033.04 ns [3023.03, 3048.14] (0.95 confidence) | — | — |
| rational&#95;solve/upper/dyadic&#95;d8/prepared | 2107.05 ns [2103.52, 2109.39] (0.95 confidence) | 2091.75 ns [2087.71, 2093.93] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d2/adapter | 721.561 ns [706.677, 735.123] (0.95 confidence) | 670.364 ns [669.25, 671.586] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d2/construct | 375.691 ns [374.179, 377.299] (0.95 confidence) | 373.254 ns [372.434, 374.022] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d2/prepared | 229.586 ns [229.182, 229.967] (0.95 confidence) | 225.816 ns [225.332, 226.35] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d3/adapter | 1189.19 ns [1187.12, 1191.95] (0.95 confidence) | 1166.85 ns [1165.12, 1168.21] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d3/construct | 692.435 ns [688.83, 694.744] (0.95 confidence) | 681.714 ns [678.229, 683.138] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d3/prepared | 390.267 ns [389.251, 390.918] (0.95 confidence) | 381.873 ns [381.302, 382.534] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d4/adapter | 1794.83 ns [1790.67, 1801.73] (0.95 confidence) | 1738.62 ns [1736.51, 1741.23] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d4/construct | 1012.43 ns [1009.08, 1016.34] (0.95 confidence) | 984.97 ns [981.226, 986.594] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d4/prepared | 623.139 ns [622.124, 624.028] (0.95 confidence) | 626.06 ns [618.077, 632.238] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d5/adapter | 2517.48 ns [2514.53, 2522.75] (0.95 confidence) | 2464.53 ns [2457.58, 2471.73] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d5/construct | 1447.89 ns [1440.34, 1452.44] (0.95 confidence) | 1418.7 ns [1412.91, 1422.5] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d5/prepared | 880.677 ns [879.756, 882.268] (0.95 confidence) | 867.639 ns [866.416, 868.254] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d6/adapter | 3427.31 ns [3422.38, 3436.8] (0.95 confidence) | 3302.36 ns [3295.52, 3308.17] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d6/construct | 2045.49 ns [2036.64, 2053.84] (0.95 confidence) | 1985.2 ns [1982.15, 1989.75] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d6/prepared | 1144.09 ns [1142.17, 1147.02] (0.95 confidence) | 1113.77 ns [1111.6, 1115.15] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d7/adapter | 4558.23 ns [4546.96, 4567.53] (0.95 confidence) | 4378.57 ns [4367.84, 4384.82] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d7/construct | 2791.18 ns [2785.01, 2797.91] (0.95 confidence) | 2673.11 ns [2660.87, 2678.36] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d7/prepared | 1479.49 ns [1475.21, 1482.46] (0.95 confidence) | 1455.68 ns [1453.71, 1456.77] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d8/adapter | 5823.13 ns [5806.47, 5853.15] (0.95 confidence) | 5615.62 ns [5605.14, 5620.13] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d8/construct | 3655.46 ns [3635.49, 3695.73] (0.95 confidence) | 3506.09 ns [3493.26, 3514.5] (0.95 confidence) | — | — |
| rational&#95;solve/upper/rational&#95;d8/prepared | 1885.47 ns [1882.52, 1889.51] (0.95 confidence) | 1837.85 ns [1835.49, 1840.01] (0.95 confidence) | — | — |

## Recorded measurement provenance

These observations belong to the measured phases. Rerendering does not recapture the publisher's environment.

### Baseline phase

| Field | Recorded value |
| --- | --- |
| Release | v0.5.1 |
| Revision | 80a67d598b8f5788c2c8ffc12a0e8677c306909e |
| Measurement mode | tag |
| Original source SHA-256 | 8edc38caa565097a74ebe12ad80493067168403c4e7bc8c4cadc1d653152892f |
| Measured source SHA-256 | cadc6c3c45e5b5a6ba194bbd9d58882ed8bba1e5f9e714dfbdfb9ea3053f416b |
| Common harness SHA-256 | e07a93b7f4dfd1ad1ad6bd25ab5c9f0134d654b8aebe075826709c7782c5be96 |
| Fingerprint schema | research-repo-tools/files/v1 |
| OS | macOS-27.2-arm64-arm-64bit-Mach-O |
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
| Release | v0.5.2 |
| Revision | b5156e222add90fc62191db55dcccb40875a5e41 |
| Measurement mode | working-tree |
| Original source SHA-256 | 80e684f9e40f0b89417aa4a8317b22e37e0f594aef8043b250656f2c6cb3990f |
| Measured source SHA-256 | 80e684f9e40f0b89417aa4a8317b22e37e0f594aef8043b250656f2c6cb3990f |
| Common harness SHA-256 | e07a93b7f4dfd1ad1ad6bd25ab5c9f0134d654b8aebe075826709c7782c5be96 |
| Fingerprint schema | research-repo-tools/files/v1 |
| OS | macOS-27.2-arm64-arm-64bit-Mach-O |
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
| Gate command (argv) | &#91;&quot;just&quot;, &quot;--justfile&quot;, &quot;tooling/performance.just&quot;, &quot;--working-directory&quot;, &quot;.&quot;, &quot;gate-public&quot;&#93; |
