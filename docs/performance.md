<!-- research-repo-tools/run-report/v3 -->
# la-stack complete benchmark measurements

Run: `2e55909621d0aacbfa7670f678574ae750311ba29353e6e9ebb11461414876a2`

Evidence: [run payload](performance-runs/runs/2e55909621d0aacbfa7670f678574ae750311ba29353e6e9ebb11461414876a2/run.json) ·
[evidence envelope](performance-runs/runs/2e55909621d0aacbfa7670f678574ae750311ba29353e6e9ebb11461414876a2/evidence.json).

Acquisition schedule: sequential; trials: 1; seed: None.
Each trial is retained separately. Within-phase samples are not independent trials; no pooled significance claim is made.

## Trial 1: baseline → current

Collection: 1 segment(s); 0 abandoned attempt(s). Validated abandoned phases are retained and excluded from accepted pairs; failed-output capture status is
recorded.
These are completed campaign slots; recovery does not establish uninterrupted collection, independent trials, or absence of optional-stopping bias.

- Collection segment `153c7233feac457685a9928eefd5a4ab` started: 1794363023829085000 (Unix ns).
- core-default/baseline: 1794363363666136000–1794369348978131000 (Unix ns).
- core-default/current: 1794369373996376000–1794375390051283000 (Unix ns).
- core-exact/baseline: 1794375599189457000–1794381303655942000 (Unix ns).
- core-exact/current: 1794381323045760000–1794387030693435000 (Unix ns).
- lu-solve/baseline: 1794387245833543000–1794387475775313000 (Unix ns).
- lu-solve/current: 1794387480574368000–1794387712280044000 (Unix ns).
- column-basis/baseline: 1794387790639980000–1794388305248041000 (Unix ns).
- column-basis/current: 1794388308379823000–1794388825718022000 (Unix ns).
- interval-solve/baseline: 1794388910264788000–1794389242861147000 (Unix ns).
- interval-solve/current: 1794389245760375000–1794389576029583000 (Unix ns).
- polynomial/baseline: 1794389657781193000–1794389712546254000 (Unix ns).
- polynomial/current: 1794389715333844000–1794389769737601000 (Unix ns).
- rational-dot/baseline: 1794389847700706000–1794390352057462000 (Unix ns).
- rational-dot/current: 1794390355296443000–1794390860904864000 (Unix ns).

Ratios are point estimates within one accepted paired suite. Marginal timing intervals do not establish significance.

## core-default

Semantic contract: la-stack-core-default-v1. Collection session: `a78ec5e5cce54fd682def359579174dd`.
Collection interval (Unix ns): 1794363023834533000–1794375423278495000.
Acquisition: newly collected.

- baseline: available — inside declared release interval; release v0.5.2; revision `eb5e44e52399401311e676faad590e6cce1610e5`.
  Source SHA-256: `bdd47c0e9bc3126dc3a32a98e0f1a688aaa52e7c6a3868009f1adfb0735f367e`; harness SHA-256:
  `0fedb76c77086b1cf8a698034ac4a1f3539fe06feda597672e2d11f01704f9cb`; scientific gate: passed.
- current: available — inside declared release interval; release v0.5.3; revision `be9277bf04400dccbe1115da69c7fda2141c4a39`.
  Source SHA-256: `793169b4c5085624d7c0e773f1d28504c017e9b78332f250150d19df87391a32`; harness SHA-256:
  `0fedb76c77086b1cf8a698034ac4a1f3539fe06feda597672e2d11f01704f9cb`; scientific gate: passed.

### Mean (ns)

<table>
<thead>
<tr>
<th>Workload</th>
<th>Baseline</th>
<th>Current</th>
<th>Baseline / current</th>
</tr>
</thead>
<tbody>
<tr>
<td>angle&#95;d3/dense/borrowed&#95;slice</td>
<td>23.7297 ns [23.6991, 23.7609] (0.95 confidence)</td>
<td>24.2735 ns [24.2418, 24.3068] (0.95 confidence)</td>
<td>0.977597</td>
</tr>
<tr>
<td>angle&#95;d3/dense/construct&#95;then&#95;angle</td>
<td>22.8851 ns [22.8533, 22.9175] (0.95 confidence)</td>
<td>23.3835 ns [23.353, 23.4139] (0.95 confidence)</td>
<td>0.978685</td>
</tr>
<tr>
<td>angle&#95;d3/dense/prepared&#95;vector</td>
<td>22.5945 ns [22.5504, 22.6366] (0.95 confidence)</td>
<td>23.2686 ns [23.2339, 23.3157] (0.95 confidence)</td>
<td>0.971029</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;antipodal/borrowed&#95;slice</td>
<td>16.7085 ns [16.6956, 16.7222] (0.95 confidence)</td>
<td>17.0881 ns [17.0422, 17.1368] (0.95 confidence)</td>
<td>0.977782</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;antipodal/construct&#95;then&#95;angle</td>
<td>15.5815 ns [15.5658, 15.5983] (0.95 confidence)</td>
<td>16.1254 ns [16.0828, 16.1731] (0.95 confidence)</td>
<td>0.966271</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;antipodal/prepared&#95;vector</td>
<td>15.7254 ns [15.7072, 15.7441] (0.95 confidence)</td>
<td>16.1762 ns [16.1517, 16.2015] (0.95 confidence)</td>
<td>0.972132</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;identical/borrowed&#95;slice</td>
<td>14.2498 ns [14.2399, 14.2616] (0.95 confidence)</td>
<td>14.4507 ns [14.4313, 14.4717] (0.95 confidence)</td>
<td>0.986101</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;identical/construct&#95;then&#95;angle</td>
<td>13.4964 ns [13.4826, 13.5136] (0.95 confidence)</td>
<td>13.6636 ns [13.6288, 13.7023] (0.95 confidence)</td>
<td>0.987765</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;identical/prepared&#95;vector</td>
<td>13.3954 ns [13.3775, 13.4126] (0.95 confidence)</td>
<td>13.5926 ns [13.5705, 13.6171] (0.95 confidence)</td>
<td>0.98549</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;varied/borrowed&#95;slice</td>
<td>26.2188 ns [26.1298, 26.3124] (0.95 confidence)</td>
<td>28.4253 ns [28.1603, 28.6933] (0.95 confidence)</td>
<td>0.922376</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;varied/construct&#95;then&#95;angle</td>
<td>25.0626 ns [25.0321, 25.0938] (0.95 confidence)</td>
<td>25.6528 ns [25.6256, 25.6807] (0.95 confidence)</td>
<td>0.976993</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;varied/prepared&#95;vector</td>
<td>24.5602 ns [24.5284, 24.5914] (0.95 confidence)</td>
<td>25.4805 ns [25.4276, 25.5358] (0.95 confidence)</td>
<td>0.963882</td>
</tr>
<tr>
<td>angle&#95;d3/maximum&#95;small&#95;angle/borrowed&#95;slice</td>
<td>10.6709 ns [10.6496, 10.6933] (0.95 confidence)</td>
<td>11.3958 ns [11.3465, 11.4492] (0.95 confidence)</td>
<td>0.936388</td>
</tr>
<tr>
<td>angle&#95;d3/maximum&#95;small&#95;angle/construct&#95;then&#95;angle</td>
<td>9.84091 ns [9.82443, 9.85787] (0.95 confidence)</td>
<td>10.266 ns [10.235, 10.2968] (0.95 confidence)</td>
<td>0.958596</td>
</tr>
<tr>
<td>angle&#95;d3/maximum&#95;small&#95;angle/prepared&#95;vector</td>
<td>9.2913 ns [9.2788, 9.30531] (0.95 confidence)</td>
<td>9.58301 ns [9.55836, 9.60767] (0.95 confidence)</td>
<td>0.96956</td>
</tr>
<tr>
<td>angle&#95;d3/mixed&#95;scale/borrowed&#95;slice</td>
<td>12.1984 ns [12.1624, 12.2379] (0.95 confidence)</td>
<td>12.7041 ns [12.6482, 12.7634] (0.95 confidence)</td>
<td>0.960198</td>
</tr>
<tr>
<td>angle&#95;d3/mixed&#95;scale/construct&#95;then&#95;angle</td>
<td>11.2355 ns [11.2175, 11.2546] (0.95 confidence)</td>
<td>11.6765 ns [11.6294, 11.7309] (0.95 confidence)</td>
<td>0.962234</td>
</tr>
<tr>
<td>angle&#95;d3/mixed&#95;scale/prepared&#95;vector</td>
<td>10.5198 ns [10.5116, 10.5284] (0.95 confidence)</td>
<td>10.8455 ns [10.8282, 10.864] (0.95 confidence)</td>
<td>0.969974</td>
</tr>
<tr>
<td>angle&#95;d3/near&#95;antipodal/borrowed&#95;slice</td>
<td>16.3337 ns [16.317, 16.3512] (0.95 confidence)</td>
<td>16.8247 ns [16.7974, 16.8545] (0.95 confidence)</td>
<td>0.970822</td>
</tr>
<tr>
<td>angle&#95;d3/near&#95;antipodal/construct&#95;then&#95;angle</td>
<td>15.1054 ns [15.0892, 15.1237] (0.95 confidence)</td>
<td>15.47 ns [15.4487, 15.4915] (0.95 confidence)</td>
<td>0.976435</td>
</tr>
<tr>
<td>angle&#95;d3/near&#95;antipodal/prepared&#95;vector</td>
<td>14.906 ns [14.8744, 14.9396] (0.95 confidence)</td>
<td>15.6536 ns [15.618, 15.6891] (0.95 confidence)</td>
<td>0.952236</td>
</tr>
<tr>
<td>angle&#95;d3/near&#95;parallel/borrowed&#95;slice</td>
<td>10.7256 ns [10.6894, 10.7639] (0.95 confidence)</td>
<td>11.2186 ns [11.175, 11.2658] (0.95 confidence)</td>
<td>0.956048</td>
</tr>
<tr>
<td>angle&#95;d3/near&#95;parallel/construct&#95;then&#95;angle</td>
<td>9.83974 ns [9.82072, 9.85944] (0.95 confidence)</td>
<td>10.2731 ns [10.2564, 10.2898] (0.95 confidence)</td>
<td>0.957813</td>
</tr>
<tr>
<td>angle&#95;d3/near&#95;parallel/prepared&#95;vector</td>
<td>9.26364 ns [9.25229, 9.27541] (0.95 confidence)</td>
<td>9.70263 ns [9.68946, 9.71595] (0.95 confidence)</td>
<td>0.954756</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;adjacent/borrowed&#95;slice</td>
<td>13.5898 ns [13.564, 13.6167] (0.95 confidence)</td>
<td>13.821 ns [13.7772, 13.8668] (0.95 confidence)</td>
<td>0.983274</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;adjacent/construct&#95;then&#95;angle</td>
<td>12.9824 ns [12.9565, 13.0089] (0.95 confidence)</td>
<td>13.0939 ns [13.0669, 13.1234] (0.95 confidence)</td>
<td>0.991484</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;adjacent/prepared&#95;vector</td>
<td>11.6502 ns [11.6417, 11.66] (0.95 confidence)</td>
<td>11.8678 ns [11.8444, 11.8944] (0.95 confidence)</td>
<td>0.981661</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;alternating/borrowed&#95;slice</td>
<td>13.9208 ns [13.8714, 13.9733] (0.95 confidence)</td>
<td>14.1064 ns [14.0458, 14.1735] (0.95 confidence)</td>
<td>0.986844</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;alternating/construct&#95;then&#95;angle</td>
<td>12.9478 ns [12.9218, 12.9752] (0.95 confidence)</td>
<td>13.1005 ns [13.0761, 13.1261] (0.95 confidence)</td>
<td>0.988345</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;alternating/prepared&#95;vector</td>
<td>11.891 ns [11.8823, 11.9011] (0.95 confidence)</td>
<td>12.0702 ns [12.0509, 12.0907] (0.95 confidence)</td>
<td>0.985156</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;antipodal/borrowed&#95;slice</td>
<td>10.7306 ns [10.6942, 10.7708] (0.95 confidence)</td>
<td>11.0649 ns [11.0119, 11.1208] (0.95 confidence)</td>
<td>0.969788</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;antipodal/construct&#95;then&#95;angle</td>
<td>10.004 ns [9.97862, 10.0295] (0.95 confidence)</td>
<td>10.3665 ns [10.3273, 10.4089] (0.95 confidence)</td>
<td>0.965037</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;antipodal/prepared&#95;vector</td>
<td>8.93321 ns [8.92412, 8.94242] (0.95 confidence)</td>
<td>9.11866 ns [9.10048, 9.13708] (0.95 confidence)</td>
<td>0.979663</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;identical/borrowed&#95;slice</td>
<td>8.66292 ns [8.62461, 8.70294] (0.95 confidence)</td>
<td>8.74886 ns [8.70636, 8.79419] (0.95 confidence)</td>
<td>0.990177</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;identical/construct&#95;then&#95;angle</td>
<td>7.79546 ns [7.77048, 7.82095] (0.95 confidence)</td>
<td>7.91748 ns [7.89119, 7.94474] (0.95 confidence)</td>
<td>0.984588</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;identical/prepared&#95;vector</td>
<td>6.70629 ns [6.70121, 6.71235] (0.95 confidence)</td>
<td>6.82254 ns [6.80921, 6.83661] (0.95 confidence)</td>
<td>0.982961</td>
</tr>
<tr>
<td>angle&#95;d3/subnormal/borrowed&#95;slice</td>
<td>10.7451 ns [10.6992, 10.7952] (0.95 confidence)</td>
<td>11.078 ns [11.0383, 11.1191] (0.95 confidence)</td>
<td>0.969941</td>
</tr>
<tr>
<td>angle&#95;d3/subnormal/construct&#95;then&#95;angle</td>
<td>9.83504 ns [9.81929, 9.85108] (0.95 confidence)</td>
<td>10.1583 ns [10.1361, 10.1806] (0.95 confidence)</td>
<td>0.968181</td>
</tr>
<tr>
<td>angle&#95;d3/subnormal/prepared&#95;vector</td>
<td>9.27391 ns [9.26419, 9.28399] (0.95 confidence)</td>
<td>9.55225 ns [9.52581, 9.58428] (0.95 confidence)</td>
<td>0.970862</td>
</tr>
<tr>
<td>angle&#95;d3/subnormal&#95;axes/borrowed&#95;slice</td>
<td>13.8185 ns [13.7786, 13.8614] (0.95 confidence)</td>
<td>14.9466 ns [14.8638, 15.0292] (0.95 confidence)</td>
<td>0.924525</td>
</tr>
<tr>
<td>angle&#95;d3/subnormal&#95;axes/construct&#95;then&#95;angle</td>
<td>13.0528 ns [13.0297, 13.0758] (0.95 confidence)</td>
<td>13.4958 ns [13.443, 13.5524] (0.95 confidence)</td>
<td>0.967178</td>
</tr>
<tr>
<td>angle&#95;d3/subnormal&#95;axes/prepared&#95;vector</td>
<td>11.9166 ns [11.9003, 11.9351] (0.95 confidence)</td>
<td>12.5921 ns [12.5129, 12.6763] (0.95 confidence)</td>
<td>0.946354</td>
</tr>
<tr>
<td>angle&#95;d4/dense/borrowed&#95;slice</td>
<td>34.0769 ns [34.0328, 34.1217] (0.95 confidence)</td>
<td>35.3288 ns [35.2846, 35.379] (0.95 confidence)</td>
<td>0.964562</td>
</tr>
<tr>
<td>angle&#95;d4/dense/construct&#95;then&#95;angle</td>
<td>32.7881 ns [32.746, 32.8327] (0.95 confidence)</td>
<td>33.411 ns [33.3717, 33.4525] (0.95 confidence)</td>
<td>0.981355</td>
</tr>
<tr>
<td>angle&#95;d4/dense/prepared&#95;vector</td>
<td>32.3368 ns [32.2859, 32.3847] (0.95 confidence)</td>
<td>33.6078 ns [33.5613, 33.6593] (0.95 confidence)</td>
<td>0.962182</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;antipodal/borrowed&#95;slice</td>
<td>26.5858 ns [26.5517, 26.6206] (0.95 confidence)</td>
<td>27.4413 ns [27.3737, 27.5184] (0.95 confidence)</td>
<td>0.968823</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;antipodal/construct&#95;then&#95;angle</td>
<td>24.9429 ns [24.9162, 24.9688] (0.95 confidence)</td>
<td>25.9598 ns [25.9153, 26.0031] (0.95 confidence)</td>
<td>0.960827</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;antipodal/prepared&#95;vector</td>
<td>24.667 ns [24.6318, 24.7025] (0.95 confidence)</td>
<td>26.5642 ns [26.4521, 26.6776] (0.95 confidence)</td>
<td>0.928582</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;identical/borrowed&#95;slice</td>
<td>24.1047 ns [24.0814, 24.1291] (0.95 confidence)</td>
<td>24.8927 ns [24.8316, 24.961] (0.95 confidence)</td>
<td>0.968343</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;identical/construct&#95;then&#95;angle</td>
<td>22.5749 ns [22.5458, 22.6039] (0.95 confidence)</td>
<td>24.1353 ns [24.0036, 24.2707] (0.95 confidence)</td>
<td>0.935346</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;identical/prepared&#95;vector</td>
<td>22.0887 ns [22.0633, 22.114] (0.95 confidence)</td>
<td>22.7044 ns [22.643, 22.7706] (0.95 confidence)</td>
<td>0.972884</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;varied/borrowed&#95;slice</td>
<td>35.4198 ns [35.3712, 35.4677] (0.95 confidence)</td>
<td>36.5325 ns [36.4718, 36.5998] (0.95 confidence)</td>
<td>0.969543</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;varied/construct&#95;then&#95;angle</td>
<td>33.933 ns [33.877, 33.9889] (0.95 confidence)</td>
<td>35.293 ns [35.2123, 35.3733] (0.95 confidence)</td>
<td>0.961466</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;varied/prepared&#95;vector</td>
<td>33.2825 ns [33.2251, 33.3402] (0.95 confidence)</td>
<td>34.2088 ns [34.1706, 34.2486] (0.95 confidence)</td>
<td>0.972922</td>
</tr>
<tr>
<td>angle&#95;d4/maximum&#95;small&#95;angle/borrowed&#95;slice</td>
<td>12.3988 ns [12.3302, 12.4667] (0.95 confidence)</td>
<td>12.6896 ns [12.6151, 12.7667] (0.95 confidence)</td>
<td>0.977082</td>
</tr>
<tr>
<td>angle&#95;d4/maximum&#95;small&#95;angle/construct&#95;then&#95;angle</td>
<td>11.1288 ns [11.0774, 11.1815] (0.95 confidence)</td>
<td>11.5431 ns [11.4982, 11.5883] (0.95 confidence)</td>
<td>0.964106</td>
</tr>
<tr>
<td>angle&#95;d4/maximum&#95;small&#95;angle/prepared&#95;vector</td>
<td>9.84237 ns [9.82371, 9.86454] (0.95 confidence)</td>
<td>10.1406 ns [10.1149, 10.169] (0.95 confidence)</td>
<td>0.970592</td>
</tr>
<tr>
<td>angle&#95;d4/mixed&#95;scale/borrowed&#95;slice</td>
<td>13.9381 ns [13.88, 13.9987] (0.95 confidence)</td>
<td>14.3119 ns [14.2448, 14.3812] (0.95 confidence)</td>
<td>0.973885</td>
</tr>
<tr>
<td>angle&#95;d4/mixed&#95;scale/construct&#95;then&#95;angle</td>
<td>12.543 ns [12.5035, 12.5828] (0.95 confidence)</td>
<td>12.9845 ns [12.94, 13.0309] (0.95 confidence)</td>
<td>0.965996</td>
</tr>
<tr>
<td>angle&#95;d4/mixed&#95;scale/prepared&#95;vector</td>
<td>11.1669 ns [11.125, 11.2115] (0.95 confidence)</td>
<td>11.4129 ns [11.3819, 11.4462] (0.95 confidence)</td>
<td>0.978445</td>
</tr>
<tr>
<td>angle&#95;d4/near&#95;antipodal/borrowed&#95;slice</td>
<td>17.6236 ns [17.6086, 17.6401] (0.95 confidence)</td>
<td>18.0458 ns [18.0256, 18.07] (0.95 confidence)</td>
<td>0.976607</td>
</tr>
<tr>
<td>angle&#95;d4/near&#95;antipodal/construct&#95;then&#95;angle</td>
<td>15.9866 ns [15.9698, 16.0047] (0.95 confidence)</td>
<td>16.3667 ns [16.35, 16.3834] (0.95 confidence)</td>
<td>0.976779</td>
</tr>
<tr>
<td>angle&#95;d4/near&#95;antipodal/prepared&#95;vector</td>
<td>15.4428 ns [15.424, 15.4615] (0.95 confidence)</td>
<td>15.8313 ns [15.8135, 15.8509] (0.95 confidence)</td>
<td>0.975461</td>
</tr>
<tr>
<td>angle&#95;d4/near&#95;parallel/borrowed&#95;slice</td>
<td>12.4252 ns [12.3557, 12.4949] (0.95 confidence)</td>
<td>12.5811 ns [12.5028, 12.6611] (0.95 confidence)</td>
<td>0.987611</td>
</tr>
<tr>
<td>angle&#95;d4/near&#95;parallel/construct&#95;then&#95;angle</td>
<td>11.1674 ns [11.1181, 11.218] (0.95 confidence)</td>
<td>11.6249 ns [11.5771, 11.6738] (0.95 confidence)</td>
<td>0.960642</td>
</tr>
<tr>
<td>angle&#95;d4/near&#95;parallel/prepared&#95;vector</td>
<td>9.87755 ns [9.85356, 9.90457] (0.95 confidence)</td>
<td>10.1608 ns [10.1251, 10.2017] (0.95 confidence)</td>
<td>0.972121</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;adjacent/borrowed&#95;slice</td>
<td>15.6733 ns [15.6276, 15.7196] (0.95 confidence)</td>
<td>16.4586 ns [16.403, 16.5175] (0.95 confidence)</td>
<td>0.952286</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;adjacent/construct&#95;then&#95;angle</td>
<td>14.5683 ns [14.5389, 14.5975] (0.95 confidence)</td>
<td>15.4753 ns [15.4408, 15.5124] (0.95 confidence)</td>
<td>0.941388</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;adjacent/prepared&#95;vector</td>
<td>13.1267 ns [13.1147, 13.1392] (0.95 confidence)</td>
<td>13.7842 ns [13.7644, 13.8053] (0.95 confidence)</td>
<td>0.952297</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;alternating/borrowed&#95;slice</td>
<td>15.9996 ns [15.9485, 16.051] (0.95 confidence)</td>
<td>16.9462 ns [16.8873, 17.0048] (0.95 confidence)</td>
<td>0.944143</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;alternating/construct&#95;then&#95;angle</td>
<td>14.9589 ns [14.9296, 14.9898] (0.95 confidence)</td>
<td>15.7241 ns [15.6923, 15.7584] (0.95 confidence)</td>
<td>0.951331</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;alternating/prepared&#95;vector</td>
<td>13.4326 ns [13.4101, 13.4573] (0.95 confidence)</td>
<td>14.02 ns [13.988, 14.0519] (0.95 confidence)</td>
<td>0.958103</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;antipodal/borrowed&#95;slice</td>
<td>11.9597 ns [11.9161, 12.0029] (0.95 confidence)</td>
<td>12.2876 ns [12.2316, 12.346] (0.95 confidence)</td>
<td>0.973308</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;antipodal/construct&#95;then&#95;angle</td>
<td>10.8864 ns [10.8551, 10.9181] (0.95 confidence)</td>
<td>11.3333 ns [11.3087, 11.3576] (0.95 confidence)</td>
<td>0.960571</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;antipodal/prepared&#95;vector</td>
<td>9.38741 ns [9.37751, 9.39821] (0.95 confidence)</td>
<td>9.65882 ns [9.6446, 9.67376] (0.95 confidence)</td>
<td>0.9719</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;identical/borrowed&#95;slice</td>
<td>10.2336 ns [10.179, 10.2872] (0.95 confidence)</td>
<td>10.5302 ns [10.4559, 10.6066] (0.95 confidence)</td>
<td>0.971832</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;identical/construct&#95;then&#95;angle</td>
<td>9.16167 ns [9.13022, 9.19293] (0.95 confidence)</td>
<td>9.54983 ns [9.52337, 9.57646] (0.95 confidence)</td>
<td>0.959354</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;identical/prepared&#95;vector</td>
<td>7.7168 ns [7.69759, 7.74246] (0.95 confidence)</td>
<td>7.96239 ns [7.93534, 7.99006] (0.95 confidence)</td>
<td>0.969155</td>
</tr>
<tr>
<td>angle&#95;d4/subnormal/borrowed&#95;slice</td>
<td>12.3767 ns [12.3108, 12.4451] (0.95 confidence)</td>
<td>12.7102 ns [12.631, 12.793] (0.95 confidence)</td>
<td>0.973758</td>
</tr>
<tr>
<td>angle&#95;d4/subnormal/construct&#95;then&#95;angle</td>
<td>11.0808 ns [11.0379, 11.1259] (0.95 confidence)</td>
<td>11.5058 ns [11.457, 11.5549] (0.95 confidence)</td>
<td>0.963059</td>
</tr>
<tr>
<td>angle&#95;d4/subnormal/prepared&#95;vector</td>
<td>9.88515 ns [9.85737, 9.91656] (0.95 confidence)</td>
<td>10.2042 ns [10.1742, 10.2422] (0.95 confidence)</td>
<td>0.968737</td>
</tr>
<tr>
<td>angle&#95;d4/subnormal&#95;axes/borrowed&#95;slice</td>
<td>15.4906 ns [15.4369, 15.5463] (0.95 confidence)</td>
<td>16.5324 ns [16.4398, 16.6267] (0.95 confidence)</td>
<td>0.936981</td>
</tr>
<tr>
<td>angle&#95;d4/subnormal&#95;axes/construct&#95;then&#95;angle</td>
<td>14.3042 ns [14.2637, 14.3476] (0.95 confidence)</td>
<td>15.0088 ns [14.9421, 15.0807] (0.95 confidence)</td>
<td>0.953058</td>
</tr>
<tr>
<td>angle&#95;d4/subnormal&#95;axes/prepared&#95;vector</td>
<td>12.8079 ns [12.7656, 12.8532] (0.95 confidence)</td>
<td>13.1641 ns [13.1201, 13.2116] (0.95 confidence)</td>
<td>0.97294</td>
</tr>
<tr>
<td>angle&#95;d5/dense/borrowed&#95;slice</td>
<td>47.7096 ns [47.6483, 47.773] (0.95 confidence)</td>
<td>49.7189 ns [49.4227, 50.0452] (0.95 confidence)</td>
<td>0.959587</td>
</tr>
<tr>
<td>angle&#95;d5/dense/construct&#95;then&#95;angle</td>
<td>45.241 ns [45.1934, 45.2862] (0.95 confidence)</td>
<td>46.7221 ns [46.6638, 46.7863] (0.95 confidence)</td>
<td>0.968301</td>
</tr>
<tr>
<td>angle&#95;d5/dense/prepared&#95;vector</td>
<td>44.8147 ns [44.7473, 44.8833] (0.95 confidence)</td>
<td>46.3215 ns [46.2738, 46.3842] (0.95 confidence)</td>
<td>0.967471</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;antipodal/borrowed&#95;slice</td>
<td>36.7206 ns [36.6791, 36.7632] (0.95 confidence)</td>
<td>37.6524 ns [37.6114, 37.6984] (0.95 confidence)</td>
<td>0.975253</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;antipodal/construct&#95;then&#95;angle</td>
<td>34.8026 ns [34.7553, 34.8474] (0.95 confidence)</td>
<td>36.1498 ns [35.9498, 36.3676] (0.95 confidence)</td>
<td>0.962731</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;antipodal/prepared&#95;vector</td>
<td>34.1694 ns [34.115, 34.2234] (0.95 confidence)</td>
<td>35.1194 ns [35.0939, 35.1486] (0.95 confidence)</td>
<td>0.972948</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;identical/borrowed&#95;slice</td>
<td>35.3536 ns [35.3022, 35.4061] (0.95 confidence)</td>
<td>37.0623 ns [36.9776, 37.1515] (0.95 confidence)</td>
<td>0.953897</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;identical/construct&#95;then&#95;angle</td>
<td>33.5875 ns [33.4825, 33.7065] (0.95 confidence)</td>
<td>36.3521 ns [35.9904, 36.7167] (0.95 confidence)</td>
<td>0.92395</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;identical/prepared&#95;vector</td>
<td>32.9016 ns [32.8655, 32.9372] (0.95 confidence)</td>
<td>34.796 ns [34.6673, 34.9346] (0.95 confidence)</td>
<td>0.945558</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;varied/borrowed&#95;slice</td>
<td>68.1352 ns [68, 68.2651] (0.95 confidence)</td>
<td>74.3455 ns [73.783, 75.0806] (0.95 confidence)</td>
<td>0.916467</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;varied/construct&#95;then&#95;angle</td>
<td>65.7938 ns [65.6639, 65.9163] (0.95 confidence)</td>
<td>70.5233 ns [70.4137, 70.6356] (0.95 confidence)</td>
<td>0.932937</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;varied/prepared&#95;vector</td>
<td>63.6491 ns [63.5343, 63.7598] (0.95 confidence)</td>
<td>69.4213 ns [68.6652, 70.6353] (0.95 confidence)</td>
<td>0.916852</td>
</tr>
<tr>
<td>angle&#95;d5/maximum&#95;small&#95;angle/borrowed&#95;slice</td>
<td>14.228 ns [14.1747, 14.2796] (0.95 confidence)</td>
<td>14.556 ns [14.4806, 14.6319] (0.95 confidence)</td>
<td>0.977464</td>
</tr>
<tr>
<td>angle&#95;d5/maximum&#95;small&#95;angle/construct&#95;then&#95;angle</td>
<td>12.6719 ns [12.6116, 12.733] (0.95 confidence)</td>
<td>13.0489 ns [12.9713, 13.1306] (0.95 confidence)</td>
<td>0.971105</td>
</tr>
<tr>
<td>angle&#95;d5/maximum&#95;small&#95;angle/prepared&#95;vector</td>
<td>10.8565 ns [10.8221, 10.8927] (0.95 confidence)</td>
<td>11.199 ns [11.159, 11.24] (0.95 confidence)</td>
<td>0.969415</td>
</tr>
<tr>
<td>angle&#95;d5/mixed&#95;scale/borrowed&#95;slice</td>
<td>15.7886 ns [15.7516, 15.8231] (0.95 confidence)</td>
<td>16.5744 ns [16.4275, 16.7304] (0.95 confidence)</td>
<td>0.95259</td>
</tr>
<tr>
<td>angle&#95;d5/mixed&#95;scale/construct&#95;then&#95;angle</td>
<td>13.9986 ns [13.9555, 14.0409] (0.95 confidence)</td>
<td>15.3632 ns [15.1698, 15.5934] (0.95 confidence)</td>
<td>0.911176</td>
</tr>
<tr>
<td>angle&#95;d5/mixed&#95;scale/prepared&#95;vector</td>
<td>12.3109 ns [12.2806, 12.3407] (0.95 confidence)</td>
<td>12.8588 ns [12.8007, 12.9196] (0.95 confidence)</td>
<td>0.957392</td>
</tr>
<tr>
<td>angle&#95;d5/near&#95;antipodal/borrowed&#95;slice</td>
<td>19.1549 ns [19.1175, 19.1876] (0.95 confidence)</td>
<td>19.6936 ns [19.64, 19.7529] (0.95 confidence)</td>
<td>0.972645</td>
</tr>
<tr>
<td>angle&#95;d5/near&#95;antipodal/construct&#95;then&#95;angle</td>
<td>16.9824 ns [16.9656, 17.001] (0.95 confidence)</td>
<td>17.4111 ns [17.3717, 17.4536] (0.95 confidence)</td>
<td>0.975382</td>
</tr>
<tr>
<td>angle&#95;d5/near&#95;antipodal/prepared&#95;vector</td>
<td>16.2335 ns [16.2126, 16.2552] (0.95 confidence)</td>
<td>16.7854 ns [16.7466, 16.8256] (0.95 confidence)</td>
<td>0.967125</td>
</tr>
<tr>
<td>angle&#95;d5/near&#95;parallel/borrowed&#95;slice</td>
<td>14.2514 ns [14.2014, 14.2988] (0.95 confidence)</td>
<td>14.5612 ns [14.4872, 14.6356] (0.95 confidence)</td>
<td>0.97872</td>
</tr>
<tr>
<td>angle&#95;d5/near&#95;parallel/construct&#95;then&#95;angle</td>
<td>12.6683 ns [12.6194, 12.7173] (0.95 confidence)</td>
<td>13.0602 ns [12.9971, 13.1239] (0.95 confidence)</td>
<td>0.969991</td>
</tr>
<tr>
<td>angle&#95;d5/near&#95;parallel/prepared&#95;vector</td>
<td>10.8753 ns [10.8502, 10.9007] (0.95 confidence)</td>
<td>11.1744 ns [11.1317, 11.2191] (0.95 confidence)</td>
<td>0.973237</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;adjacent/borrowed&#95;slice</td>
<td>17.4285 ns [17.3845, 17.4794] (0.95 confidence)</td>
<td>18.3518 ns [18.1999, 18.5111] (0.95 confidence)</td>
<td>0.949686</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;adjacent/construct&#95;then&#95;angle</td>
<td>16.2344 ns [16.18, 16.2925] (0.95 confidence)</td>
<td>16.7535 ns [16.7124, 16.7948] (0.95 confidence)</td>
<td>0.969013</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;adjacent/prepared&#95;vector</td>
<td>14.44 ns [14.416, 14.4666] (0.95 confidence)</td>
<td>14.9318 ns [14.8348, 15.0391] (0.95 confidence)</td>
<td>0.967064</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;alternating/borrowed&#95;slice</td>
<td>27.4414 ns [27.4124, 27.4696] (0.95 confidence)</td>
<td>28.109 ns [28.0821, 28.1386] (0.95 confidence)</td>
<td>0.976248</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;alternating/construct&#95;then&#95;angle</td>
<td>26.437 ns [26.4106, 26.4651] (0.95 confidence)</td>
<td>27.1732 ns [27.1228, 27.2296] (0.95 confidence)</td>
<td>0.97291</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;alternating/prepared&#95;vector</td>
<td>25.8953 ns [25.8634, 25.9283] (0.95 confidence)</td>
<td>26.4826 ns [26.464, 26.5057] (0.95 confidence)</td>
<td>0.977823</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;antipodal/borrowed&#95;slice</td>
<td>13.0971 ns [13.0806, 13.1143] (0.95 confidence)</td>
<td>13.6837 ns [13.6318, 13.7399] (0.95 confidence)</td>
<td>0.957128</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;antipodal/construct&#95;then&#95;angle</td>
<td>11.9274 ns [11.8981, 11.9569] (0.95 confidence)</td>
<td>12.3061 ns [12.2742, 12.3377] (0.95 confidence)</td>
<td>0.969227</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;antipodal/prepared&#95;vector</td>
<td>10.1538 ns [10.1439, 10.1635] (0.95 confidence)</td>
<td>10.566 ns [10.5327, 10.6009] (0.95 confidence)</td>
<td>0.960981</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;identical/borrowed&#95;slice</td>
<td>11.723 ns [11.683, 11.767] (0.95 confidence)</td>
<td>12.0953 ns [12.036, 12.1601] (0.95 confidence)</td>
<td>0.96922</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;identical/construct&#95;then&#95;angle</td>
<td>10.4386 ns [10.3956, 10.4816] (0.95 confidence)</td>
<td>10.9181 ns [10.8659, 10.9721] (0.95 confidence)</td>
<td>0.956078</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;identical/prepared&#95;vector</td>
<td>8.74639 ns [8.72788, 8.7667] (0.95 confidence)</td>
<td>8.92012 ns [8.90356, 8.93825] (0.95 confidence)</td>
<td>0.980524</td>
</tr>
<tr>
<td>angle&#95;d5/subnormal/borrowed&#95;slice</td>
<td>14.3479 ns [14.2901, 14.4048] (0.95 confidence)</td>
<td>14.6097 ns [14.5467, 14.6739] (0.95 confidence)</td>
<td>0.982084</td>
</tr>
<tr>
<td>angle&#95;d5/subnormal/construct&#95;then&#95;angle</td>
<td>12.6892 ns [12.6312, 12.7481] (0.95 confidence)</td>
<td>12.7046 ns [12.6357, 12.7722] (0.95 confidence)</td>
<td>0.998786</td>
</tr>
<tr>
<td>angle&#95;d5/subnormal/prepared&#95;vector</td>
<td>10.9472 ns [10.9103, 10.9858] (0.95 confidence)</td>
<td>11.2011 ns [11.1602, 11.2437] (0.95 confidence)</td>
<td>0.977335</td>
</tr>
<tr>
<td>angle&#95;d5/subnormal&#95;axes/borrowed&#95;slice</td>
<td>17.428 ns [17.384, 17.47] (0.95 confidence)</td>
<td>17.8258 ns [17.7587, 17.8919] (0.95 confidence)</td>
<td>0.97768</td>
</tr>
<tr>
<td>angle&#95;d5/subnormal&#95;axes/construct&#95;then&#95;angle</td>
<td>15.7627 ns [15.7075, 15.8174] (0.95 confidence)</td>
<td>16.3494 ns [16.2838, 16.4167] (0.95 confidence)</td>
<td>0.964116</td>
</tr>
<tr>
<td>angle&#95;d5/subnormal&#95;axes/prepared&#95;vector</td>
<td>14.0215 ns [13.9783, 14.0645] (0.95 confidence)</td>
<td>14.4619 ns [14.3962, 14.5297] (0.95 confidence)</td>
<td>0.96955</td>
</tr>
<tr>
<td>angle&#95;d6/dense/borrowed&#95;slice</td>
<td>62.759 ns [62.6581, 62.8611] (0.95 confidence)</td>
<td>63.6723 ns [63.6333, 63.7186] (0.95 confidence)</td>
<td>0.985658</td>
</tr>
<tr>
<td>angle&#95;d6/dense/construct&#95;then&#95;angle</td>
<td>59.6193 ns [59.5297, 59.7121] (0.95 confidence)</td>
<td>60.9514 ns [60.9193, 60.9891] (0.95 confidence)</td>
<td>0.978145</td>
</tr>
<tr>
<td>angle&#95;d6/dense/prepared&#95;vector</td>
<td>59.1116 ns [59.022, 59.1965] (0.95 confidence)</td>
<td>60.5478 ns [60.4993, 60.6086] (0.95 confidence)</td>
<td>0.976281</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;antipodal/borrowed&#95;slice</td>
<td>48.7574 ns [48.6944, 48.82] (0.95 confidence)</td>
<td>49.6893 ns [49.6529, 49.7326] (0.95 confidence)</td>
<td>0.981246</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;antipodal/construct&#95;then&#95;angle</td>
<td>46.1565 ns [46.0971, 46.2165] (0.95 confidence)</td>
<td>46.9196 ns [46.8928, 46.9516] (0.95 confidence)</td>
<td>0.983736</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;antipodal/prepared&#95;vector</td>
<td>45.6289 ns [45.5637, 45.6906] (0.95 confidence)</td>
<td>46.5988 ns [46.565, 46.6397] (0.95 confidence)</td>
<td>0.979184</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;identical/borrowed&#95;slice</td>
<td>47.6227 ns [47.5548, 47.6897] (0.95 confidence)</td>
<td>48.6527 ns [48.6205, 48.693] (0.95 confidence)</td>
<td>0.978829</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;identical/construct&#95;then&#95;angle</td>
<td>45.145 ns [45.0883, 45.2024] (0.95 confidence)</td>
<td>45.9302 ns [45.8905, 45.9834] (0.95 confidence)</td>
<td>0.982906</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;identical/prepared&#95;vector</td>
<td>44.4939 ns [44.4416, 44.5472] (0.95 confidence)</td>
<td>45.4333 ns [45.4082, 45.4637] (0.95 confidence)</td>
<td>0.979322</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;varied/borrowed&#95;slice</td>
<td>93.3557 ns [93.1138, 93.6085] (0.95 confidence)</td>
<td>99.0541 ns [98.8307, 99.2987] (0.95 confidence)</td>
<td>0.942471</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;varied/construct&#95;then&#95;angle</td>
<td>91.153 ns [90.9095, 91.4005] (0.95 confidence)</td>
<td>96.5849 ns [96.3829, 96.8082] (0.95 confidence)</td>
<td>0.943761</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;varied/prepared&#95;vector</td>
<td>89.5233 ns [89.3526, 89.6906] (0.95 confidence)</td>
<td>95.0904 ns [94.9924, 95.2015] (0.95 confidence)</td>
<td>0.941455</td>
</tr>
<tr>
<td>angle&#95;d6/maximum&#95;small&#95;angle/borrowed&#95;slice</td>
<td>15.5901 ns [15.5602, 15.6172] (0.95 confidence)</td>
<td>22.4098 ns [21.8392, 23.0641] (0.95 confidence)</td>
<td>0.695682</td>
</tr>
<tr>
<td>angle&#95;d6/maximum&#95;small&#95;angle/construct&#95;then&#95;angle</td>
<td>13.6978 ns [13.6616, 13.7373] (0.95 confidence)</td>
<td>22.1677 ns [21.6359, 22.7143] (0.95 confidence)</td>
<td>0.617917</td>
</tr>
<tr>
<td>angle&#95;d6/maximum&#95;small&#95;angle/prepared&#95;vector</td>
<td>11.9528 ns [11.9228, 11.9836] (0.95 confidence)</td>
<td>17.2217 ns [16.7542, 17.7389] (0.95 confidence)</td>
<td>0.694052</td>
</tr>
<tr>
<td>angle&#95;d6/mixed&#95;scale/borrowed&#95;slice</td>
<td>17.0672 ns [17.0356, 17.0997] (0.95 confidence)</td>
<td>17.4423 ns [17.2858, 17.6214] (0.95 confidence)</td>
<td>0.978499</td>
</tr>
<tr>
<td>angle&#95;d6/mixed&#95;scale/construct&#95;then&#95;angle</td>
<td>15.3544 ns [15.2874, 15.4261] (0.95 confidence)</td>
<td>25.733 ns [24.6855, 26.8126] (0.95 confidence)</td>
<td>0.596684</td>
</tr>
<tr>
<td>angle&#95;d6/mixed&#95;scale/prepared&#95;vector</td>
<td>13.406 ns [13.383, 13.4292] (0.95 confidence)</td>
<td>17.4928 ns [16.8514, 18.1478] (0.95 confidence)</td>
<td>0.766374</td>
</tr>
<tr>
<td>angle&#95;d6/near&#95;antipodal/borrowed&#95;slice</td>
<td>20.7838 ns [20.7476, 20.8279] (0.95 confidence)</td>
<td>24.1453 ns [23.5075, 24.8018] (0.95 confidence)</td>
<td>0.86078</td>
</tr>
<tr>
<td>angle&#95;d6/near&#95;antipodal/construct&#95;then&#95;angle</td>
<td>18.3834 ns [18.3704, 18.3983] (0.95 confidence)</td>
<td>30.988 ns [29.9577, 32.0551] (0.95 confidence)</td>
<td>0.593245</td>
</tr>
<tr>
<td>angle&#95;d6/near&#95;antipodal/prepared&#95;vector</td>
<td>17.0593 ns [17.034, 17.0861] (0.95 confidence)</td>
<td>28.9751 ns [28.104, 29.9762] (0.95 confidence)</td>
<td>0.588757</td>
</tr>
<tr>
<td>angle&#95;d6/near&#95;parallel/borrowed&#95;slice</td>
<td>15.5631 ns [15.5215, 15.6023] (0.95 confidence)</td>
<td>25.2883 ns [24.4199, 26.2945] (0.95 confidence)</td>
<td>0.615426</td>
</tr>
<tr>
<td>angle&#95;d6/near&#95;parallel/construct&#95;then&#95;angle</td>
<td>13.7163 ns [13.6686, 13.7668] (0.95 confidence)</td>
<td>23.6085 ns [23.0169, 24.2317] (0.95 confidence)</td>
<td>0.580992</td>
</tr>
<tr>
<td>angle&#95;d6/near&#95;parallel/prepared&#95;vector</td>
<td>11.9182 ns [11.8908, 11.9439] (0.95 confidence)</td>
<td>19.9832 ns [19.147, 21.0048] (0.95 confidence)</td>
<td>0.59641</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;adjacent/borrowed&#95;slice</td>
<td>19.1539 ns [19.082, 19.2301] (0.95 confidence)</td>
<td>19.4284 ns [19.3442, 19.518] (0.95 confidence)</td>
<td>0.98587</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;adjacent/construct&#95;then&#95;angle</td>
<td>18.052 ns [17.8993, 18.2685] (0.95 confidence)</td>
<td>20.6769 ns [20.1249, 21.3191] (0.95 confidence)</td>
<td>0.873051</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;adjacent/prepared&#95;vector</td>
<td>15.8411 ns [15.7906, 15.8948] (0.95 confidence)</td>
<td>15.9815 ns [15.9378, 16.0312] (0.95 confidence)</td>
<td>0.991216</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;alternating/borrowed&#95;slice</td>
<td>29.0033 ns [28.9811, 29.0249] (0.95 confidence)</td>
<td>33.7994 ns [32.8229, 34.817] (0.95 confidence)</td>
<td>0.858101</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;alternating/construct&#95;then&#95;angle</td>
<td>27.5898 ns [27.5546, 27.6263] (0.95 confidence)</td>
<td>35.3133 ns [33.9082, 37.0562] (0.95 confidence)</td>
<td>0.781286</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;alternating/prepared&#95;vector</td>
<td>26.8592 ns [26.8184, 26.9024] (0.95 confidence)</td>
<td>27.7299 ns [27.6723, 27.7877] (0.95 confidence)</td>
<td>0.9686</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;antipodal/borrowed&#95;slice</td>
<td>14.0947 ns [14.0735, 14.1181] (0.95 confidence)</td>
<td>18.3023 ns [17.8545, 18.8122] (0.95 confidence)</td>
<td>0.770105</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;antipodal/construct&#95;then&#95;angle</td>
<td>12.7896 ns [12.7336, 12.8491] (0.95 confidence)</td>
<td>16.9355 ns [16.5521, 17.3303] (0.95 confidence)</td>
<td>0.755199</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;antipodal/prepared&#95;vector</td>
<td>10.8757 ns [10.8635, 10.8883] (0.95 confidence)</td>
<td>12.61 ns [12.3623, 12.8683] (0.95 confidence)</td>
<td>0.862466</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;identical/borrowed&#95;slice</td>
<td>13.1925 ns [13.12, 13.2656] (0.95 confidence)</td>
<td>13.6538 ns [13.5715, 13.7381] (0.95 confidence)</td>
<td>0.966215</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;identical/construct&#95;then&#95;angle</td>
<td>11.5046 ns [11.4167, 11.6002] (0.95 confidence)</td>
<td>12.0447 ns [11.9772, 12.1182] (0.95 confidence)</td>
<td>0.955164</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;identical/prepared&#95;vector</td>
<td>9.50836 ns [9.47745, 9.54177] (0.95 confidence)</td>
<td>9.74155 ns [9.70308, 9.78318] (0.95 confidence)</td>
<td>0.976062</td>
</tr>
<tr>
<td>angle&#95;d6/subnormal/borrowed&#95;slice</td>
<td>15.5265 ns [15.4827, 15.5688] (0.95 confidence)</td>
<td>15.8124 ns [15.7574, 15.8696] (0.95 confidence)</td>
<td>0.981922</td>
</tr>
<tr>
<td>angle&#95;d6/subnormal/construct&#95;then&#95;angle</td>
<td>13.7678 ns [13.696, 13.867] (0.95 confidence)</td>
<td>14.0969 ns [14.0423, 14.1529] (0.95 confidence)</td>
<td>0.976659</td>
</tr>
<tr>
<td>angle&#95;d6/subnormal/prepared&#95;vector</td>
<td>11.9682 ns [11.9477, 11.9902] (0.95 confidence)</td>
<td>12.2918 ns [12.2116, 12.3795] (0.95 confidence)</td>
<td>0.973672</td>
</tr>
<tr>
<td>angle&#95;d6/subnormal&#95;axes/borrowed&#95;slice</td>
<td>18.7551 ns [18.7118, 18.7977] (0.95 confidence)</td>
<td>25.3652 ns [24.5874, 26.2553] (0.95 confidence)</td>
<td>0.739401</td>
</tr>
<tr>
<td>angle&#95;d6/subnormal&#95;axes/construct&#95;then&#95;angle</td>
<td>16.9599 ns [16.9205, 17.0028] (0.95 confidence)</td>
<td>23.4595 ns [22.9025, 24.0363] (0.95 confidence)</td>
<td>0.722943</td>
</tr>
<tr>
<td>angle&#95;d6/subnormal&#95;axes/prepared&#95;vector</td>
<td>15.1372 ns [15.0992, 15.1777] (0.95 confidence)</td>
<td>20.6367 ns [19.9911, 21.3838] (0.95 confidence)</td>
<td>0.733508</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;det</td>
<td>374.855 ns [371.374, 378.615] (0.95 confidence)</td>
<td>395.224 ns [387.437, 403.539] (0.95 confidence)</td>
<td>0.948461</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>2.89288 ns [2.88848, 2.89713] (0.95 confidence)</td>
<td>2.95701 ns [2.95412, 2.96014] (0.95 confidence)</td>
<td>0.978312</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>3.28242 ns [3.27947, 3.28559] (0.95 confidence)</td>
<td>3.35285 ns [3.34819, 3.35804] (0.95 confidence)</td>
<td>0.978993</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>444.311 ns [434.888, 453.512] (0.95 confidence)</td>
<td>483.809 ns [475.198, 492.262] (0.95 confidence)</td>
<td>0.91836</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;dot</td>
<td>2.32946 ns [2.32698, 2.33208] (0.95 confidence)</td>
<td>2.37151 ns [2.36747, 2.3757] (0.95 confidence)</td>
<td>0.982269</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;inf&#95;norm</td>
<td>32.3368 ns [32.2941, 32.3795] (0.95 confidence)</td>
<td>33.0075 ns [32.9516, 33.0676] (0.95 confidence)</td>
<td>0.97968</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;ldlt</td>
<td>340.861 ns [340.45, 341.268] (0.95 confidence)</td>
<td>347.722 ns [347.223, 348.228] (0.95 confidence)</td>
<td>0.980269</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;ldlt&#95;solve</td>
<td>382.64 ns [382.317, 383.007] (0.95 confidence)</td>
<td>399.681 ns [397.002, 402.535] (0.95 confidence)</td>
<td>0.957361</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;lu</td>
<td>385.034 ns [381.073, 389.388] (0.95 confidence)</td>
<td>392.92 ns [389.484, 396.245] (0.95 confidence)</td>
<td>0.979929</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;norm2</td>
<td>10.0129 ns [9.98617, 10.0403] (0.95 confidence)</td>
<td>10.6307 ns [10.5979, 10.6643] (0.95 confidence)</td>
<td>0.941885</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>9.39058 ns [9.36687, 9.41416] (0.95 confidence)</td>
<td>9.85871 ns [9.83984, 9.8782] (0.95 confidence)</td>
<td>0.952516</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>9.38819 ns [9.36825, 9.40784] (0.95 confidence)</td>
<td>9.75952 ns [9.73839, 9.78055] (0.95 confidence)</td>
<td>0.961952</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>4.85249 ns [4.78713, 4.91603] (0.95 confidence)</td>
<td>4.75097 ns [4.68104, 4.8249] (0.95 confidence)</td>
<td>1.02137</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>8.2441 ns [8.22874, 8.2587] (0.95 confidence)</td>
<td>8.5696 ns [8.55109, 8.58958] (0.95 confidence)</td>
<td>0.962017</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;norm2&#95;sq</td>
<td>2.08521 ns [2.07268, 2.09511] (0.95 confidence)</td>
<td>2.11561 ns [2.1073, 2.12322] (0.95 confidence)</td>
<td>0.985631</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>27.758 ns [27.7391, 27.7791] (0.95 confidence)</td>
<td>29.1701 ns [28.9472, 29.4054] (0.95 confidence)</td>
<td>0.951591</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>190.198 ns [189.8, 190.59] (0.95 confidence)</td>
<td>211.514 ns [210.252, 212.916] (0.95 confidence)</td>
<td>0.899224</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;det</td>
<td>0.572636 ns [0.572213, 0.573133] (0.95 confidence)</td>
<td>0.57237 ns [0.572008, 0.572793] (0.95 confidence)</td>
<td>1.00047</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>0.502827 ns [0.499905, 0.506032] (0.95 confidence)</td>
<td>0.511473 ns [0.507812, 0.515463] (0.95 confidence)</td>
<td>0.983096</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>0.520444 ns [0.520098, 0.520874] (0.95 confidence)</td>
<td>0.521134 ns [0.520447, 0.521967] (0.95 confidence)</td>
<td>0.998676</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>1.61143 ns [1.60215, 1.62202] (0.95 confidence)</td>
<td>1.59466 ns [1.59099, 1.59875] (0.95 confidence)</td>
<td>1.01052</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;dot</td>
<td>0.612739 ns [0.607959, 0.618659] (0.95 confidence)</td>
<td>0.638917 ns [0.625654, 0.661144] (0.95 confidence)</td>
<td>0.959028</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;inf&#95;norm</td>
<td>0.702943 ns [0.702449, 0.703538] (0.95 confidence)</td>
<td>0.705748 ns [0.704624, 0.707024] (0.95 confidence)</td>
<td>0.996025</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;ldlt</td>
<td>6.90915 ns [6.8729, 6.94778] (0.95 confidence)</td>
<td>6.50137 ns [6.48792, 6.51593] (0.95 confidence)</td>
<td>1.06272</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;ldlt&#95;solve</td>
<td>10.1798 ns [10.051, 10.314] (0.95 confidence)</td>
<td>9.33263 ns [9.31921, 9.34749] (0.95 confidence)</td>
<td>1.09077</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;lu</td>
<td>1.38314 ns [1.38024, 1.38592] (0.95 confidence)</td>
<td>1.39254 ns [1.38926, 1.39626] (0.95 confidence)</td>
<td>0.993243</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;norm2</td>
<td>1.15876 ns [1.15441, 1.16354] (0.95 confidence)</td>
<td>1.30788 ns [1.27994, 1.35705] (0.95 confidence)</td>
<td>0.885985</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>1.34324 ns [1.34123, 1.34522] (0.95 confidence)</td>
<td>1.44751 ns [1.4373, 1.45825] (0.95 confidence)</td>
<td>0.927963</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>1.33949 ns [1.33786, 1.341] (0.95 confidence)</td>
<td>1.40758 ns [1.40534, 1.40983] (0.95 confidence)</td>
<td>0.951628</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>1.18126 ns [1.17392, 1.18878] (0.95 confidence)</td>
<td>1.20492 ns [1.19389, 1.21731] (0.95 confidence)</td>
<td>0.980367</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>1.34006 ns [1.33775, 1.34244] (0.95 confidence)</td>
<td>1.36272 ns [1.36151, 1.36405] (0.95 confidence)</td>
<td>0.983373</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;norm2&#95;sq</td>
<td>0.394321 ns [0.393786, 0.394842] (0.95 confidence)</td>
<td>0.406735 ns [0.403596, 0.410762] (0.95 confidence)</td>
<td>0.96948</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>1.2228 ns [1.21505, 1.22946] (0.95 confidence)</td>
<td>1.21673 ns [1.2091, 1.22367] (0.95 confidence)</td>
<td>1.00499</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>1.2525 ns [1.24938, 1.25581] (0.95 confidence)</td>
<td>1.25864 ns [1.25529, 1.26209] (0.95 confidence)</td>
<td>0.995122</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;det</td>
<td>1.35565 ns [1.34959, 1.36209] (0.95 confidence)</td>
<td>1.35378 ns [1.34759, 1.36042] (0.95 confidence)</td>
<td>1.00139</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>0.594526 ns [0.594178, 0.594959] (0.95 confidence)</td>
<td>0.594964 ns [0.594458, 0.595565] (0.95 confidence)</td>
<td>0.999265</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>0.693738 ns [0.693315, 0.694255] (0.95 confidence)</td>
<td>0.694516 ns [0.693859, 0.695298] (0.95 confidence)</td>
<td>0.99888</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>8.70273 ns [8.68872, 8.7167] (0.95 confidence)</td>
<td>8.72527 ns [8.71319, 8.73765] (0.95 confidence)</td>
<td>0.997416</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;dot</td>
<td>0.66759 ns [0.666833, 0.66838] (0.95 confidence)</td>
<td>0.667112 ns [0.666423, 0.667835] (0.95 confidence)</td>
<td>1.00072</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;inf&#95;norm</td>
<td>1.36308 ns [1.3622, 1.36408] (0.95 confidence)</td>
<td>1.37069 ns [1.36924, 1.37235] (0.95 confidence)</td>
<td>0.994451</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;ldlt</td>
<td>8.47433 ns [8.45691, 8.4883] (0.95 confidence)</td>
<td>8.50994 ns [8.50321, 8.51719] (0.95 confidence)</td>
<td>0.995815</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;ldlt&#95;solve</td>
<td>9.33281 ns [9.32062, 9.34446] (0.95 confidence)</td>
<td>9.35159 ns [9.33677, 9.36608] (0.95 confidence)</td>
<td>0.997992</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;lu</td>
<td>8.33935 ns [8.32353, 8.3549] (0.95 confidence)</td>
<td>8.39784 ns [8.38012, 8.415] (0.95 confidence)</td>
<td>0.993035</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;norm2</td>
<td>1.84484 ns [1.81055, 1.8822] (0.95 confidence)</td>
<td>1.91848 ns [1.8731, 1.96847] (0.95 confidence)</td>
<td>0.961614</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>1.83388 ns [1.83106, 1.83664] (0.95 confidence)</td>
<td>1.86963 ns [1.86693, 1.87254] (0.95 confidence)</td>
<td>0.98088</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>1.83132 ns [1.82859, 1.83403] (0.95 confidence)</td>
<td>1.87003 ns [1.8653, 1.87516] (0.95 confidence)</td>
<td>0.979301</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>1.41185 ns [1.40319, 1.42076] (0.95 confidence)</td>
<td>1.42171 ns [1.41163, 1.43217] (0.95 confidence)</td>
<td>0.993064</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>1.83234 ns [1.82914, 1.83535] (0.95 confidence)</td>
<td>1.85886 ns [1.85577, 1.86187] (0.95 confidence)</td>
<td>0.985734</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;norm2&#95;sq</td>
<td>0.444011 ns [0.443572, 0.444487] (0.95 confidence)</td>
<td>0.449659 ns [0.448114, 0.451345] (0.95 confidence)</td>
<td>0.987439</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>1.78035 ns [1.77857, 1.7825] (0.95 confidence)</td>
<td>1.78032 ns [1.77885, 1.7821] (0.95 confidence)</td>
<td>1.00002</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>2.04743 ns [2.04439, 2.05075] (0.95 confidence)</td>
<td>2.06053 ns [2.0565, 2.06509] (0.95 confidence)</td>
<td>0.993642</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;det</td>
<td>2000.19 ns [1996.51, 2003.91] (0.95 confidence)</td>
<td>2033.3 ns [2029.28, 2037.38] (0.95 confidence)</td>
<td>0.983713</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>7.7604 ns [7.74486, 7.77505] (0.95 confidence)</td>
<td>7.96951 ns [7.95397, 7.98519] (0.95 confidence)</td>
<td>0.973761</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>9.10293 ns [9.08373, 9.12382] (0.95 confidence)</td>
<td>9.2461 ns [9.22249, 9.27204] (0.95 confidence)</td>
<td>0.984515</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>2038.57 ns [2035.19, 2042.15] (0.95 confidence)</td>
<td>2077.69 ns [2072.41, 2083.68] (0.95 confidence)</td>
<td>0.981171</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;dot</td>
<td>4.01633 ns [4.01131, 4.02148] (0.95 confidence)</td>
<td>4.05481 ns [4.05009, 4.05996] (0.95 confidence)</td>
<td>0.990511</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;inf&#95;norm</td>
<td>127.977 ns [127.761, 128.191] (0.95 confidence)</td>
<td>130.984 ns [130.766, 131.218] (0.95 confidence)</td>
<td>0.97704</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;ldlt</td>
<td>2393.54 ns [2391.48, 2395.84] (0.95 confidence)</td>
<td>2424.93 ns [2421.51, 2428.77] (0.95 confidence)</td>
<td>0.987054</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;ldlt&#95;solve</td>
<td>2726.63 ns [2722.03, 2731.47] (0.95 confidence)</td>
<td>2761.06 ns [2755.59, 2766.83] (0.95 confidence)</td>
<td>0.987532</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;lu</td>
<td>2227.12 ns [2213.3, 2241.39] (0.95 confidence)</td>
<td>2232.3 ns [2218.92, 2246.7] (0.95 confidence)</td>
<td>0.997677</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;norm2</td>
<td>23.1163 ns [23.0511, 23.1907] (0.95 confidence)</td>
<td>23.7551 ns [23.7015, 23.8091] (0.95 confidence)</td>
<td>0.973108</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>18.5428 ns [18.5087, 18.5746] (0.95 confidence)</td>
<td>18.9495 ns [18.9023, 18.994] (0.95 confidence)</td>
<td>0.978536</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>18.4703 ns [18.4383, 18.5] (0.95 confidence)</td>
<td>18.9568 ns [18.9237, 18.9907] (0.95 confidence)</td>
<td>0.974336</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>9.19407 ns [9.1708, 9.2166] (0.95 confidence)</td>
<td>9.37401 ns [9.33696, 9.41187] (0.95 confidence)</td>
<td>0.980804</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>16.1136 ns [16.0772, 16.1456] (0.95 confidence)</td>
<td>16.6299 ns [16.582, 16.6804] (0.95 confidence)</td>
<td>0.968953</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;norm2&#95;sq</td>
<td>3.98908 ns [3.97402, 4.00755] (0.95 confidence)</td>
<td>4.0087 ns [4.00524, 4.01236] (0.95 confidence)</td>
<td>0.995108</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>312.223 ns [311.419, 313.139] (0.95 confidence)</td>
<td>317.199 ns [316.394, 318.145] (0.95 confidence)</td>
<td>0.984315</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>665.579 ns [663.99, 667.261] (0.95 confidence)</td>
<td>693.557 ns [690.504, 697.229] (0.95 confidence)</td>
<td>0.959659</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;det</td>
<td>2.54562 ns [2.5404, 2.55068] (0.95 confidence)</td>
<td>2.53778 ns [2.53114, 2.54436] (0.95 confidence)</td>
<td>1.00309</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>0.767976 ns [0.767443, 0.76866] (0.95 confidence)</td>
<td>0.768048 ns [0.767605, 0.768622] (0.95 confidence)</td>
<td>0.999906</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>0.867199 ns [0.866688, 0.867816] (0.95 confidence)</td>
<td>0.876146 ns [0.874217, 0.878292] (0.95 confidence)</td>
<td>0.989788</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>14.5563 ns [14.5288, 14.5857] (0.95 confidence)</td>
<td>14.5679 ns [14.542, 14.5961] (0.95 confidence)</td>
<td>0.9992</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;dot</td>
<td>0.682833 ns [0.681801, 0.683961] (0.95 confidence)</td>
<td>0.680809 ns [0.679894, 0.68176] (0.95 confidence)</td>
<td>1.00297</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;inf&#95;norm</td>
<td>2.35001 ns [2.34782, 2.35232] (0.95 confidence)</td>
<td>2.35202 ns [2.35037, 2.35399] (0.95 confidence)</td>
<td>0.999145</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;ldlt</td>
<td>18.0612 ns [18.0319, 18.0899] (0.95 confidence)</td>
<td>18.1801 ns [18.0999, 18.2728] (0.95 confidence)</td>
<td>0.993457</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;ldlt&#95;solve</td>
<td>16.5245 ns [16.5045, 16.5452] (0.95 confidence)</td>
<td>16.5566 ns [16.5295, 16.584] (0.95 confidence)</td>
<td>0.998065</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;lu</td>
<td>14.2179 ns [14.2021, 14.2352] (0.95 confidence)</td>
<td>14.209 ns [14.1899, 14.2302] (0.95 confidence)</td>
<td>1.00063</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;norm2</td>
<td>2.39025 ns [2.36503, 2.41517] (0.95 confidence)</td>
<td>2.32069 ns [2.28601, 2.36995] (0.95 confidence)</td>
<td>1.02997</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>2.39957 ns [2.39649, 2.40266] (0.95 confidence)</td>
<td>2.417 ns [2.41325, 2.42059] (0.95 confidence)</td>
<td>0.992788</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>2.39088 ns [2.38669, 2.39506] (0.95 confidence)</td>
<td>2.41059 ns [2.40629, 2.41492] (0.95 confidence)</td>
<td>0.991823</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>1.82327 ns [1.75898, 1.88673] (0.95 confidence)</td>
<td>1.60703 ns [1.56114, 1.6566] (0.95 confidence)</td>
<td>1.13456</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>2.13926 ns [2.13544, 2.14308] (0.95 confidence)</td>
<td>2.15514 ns [2.15159, 2.1587] (0.95 confidence)</td>
<td>0.992631</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;norm2&#95;sq</td>
<td>0.460274 ns [0.459021, 0.461597] (0.95 confidence)</td>
<td>0.459757 ns [0.458163, 0.461513] (0.95 confidence)</td>
<td>1.00112</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>2.48364 ns [2.4737, 2.49471] (0.95 confidence)</td>
<td>2.49432 ns [2.48625, 2.50341] (0.95 confidence)</td>
<td>0.995719</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>4.01104 ns [4.00513, 4.01779] (0.95 confidence)</td>
<td>4.03184 ns [4.02207, 4.04314] (0.95 confidence)</td>
<td>0.994839</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;det</td>
<td>28.9008 ns [28.8612, 28.9415] (0.95 confidence)</td>
<td>28.8199 ns [28.7977, 28.8445] (0.95 confidence)</td>
<td>1.0028</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>0.877292 ns [0.876613, 0.878072] (0.95 confidence)</td>
<td>0.876909 ns [0.876397, 0.877525] (0.95 confidence)</td>
<td>1.00044</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>1.03347 ns [1.03263, 1.03454] (0.95 confidence)</td>
<td>1.03369 ns [1.03296, 1.03458] (0.95 confidence)</td>
<td>0.999793</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>34.1699 ns [34.1208, 34.2204] (0.95 confidence)</td>
<td>34.1926 ns [34.1243, 34.259] (0.95 confidence)</td>
<td>0.999337</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;dot</td>
<td>0.822208 ns [0.818422, 0.825951] (0.95 confidence)</td>
<td>0.811684 ns [0.806458, 0.816819] (0.95 confidence)</td>
<td>1.01297</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;inf&#95;norm</td>
<td>3.46287 ns [3.46005, 3.46581] (0.95 confidence)</td>
<td>3.46864 ns [3.46408, 3.47335] (0.95 confidence)</td>
<td>0.998336</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;ldlt</td>
<td>24.5145 ns [24.4682, 24.567] (0.95 confidence)</td>
<td>24.5344 ns [24.5069, 24.5642] (0.95 confidence)</td>
<td>0.99919</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;ldlt&#95;solve</td>
<td>30.1646 ns [30.1265, 30.2008] (0.95 confidence)</td>
<td>30.1447 ns [30.0983, 30.1971] (0.95 confidence)</td>
<td>1.00066</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;lu</td>
<td>33.8969 ns [33.8535, 33.9412] (0.95 confidence)</td>
<td>33.8196 ns [33.7419, 33.8966] (0.95 confidence)</td>
<td>1.00229</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;norm2</td>
<td>2.83544 ns [2.80589, 2.86674] (0.95 confidence)</td>
<td>2.86038 ns [2.83635, 2.88433] (0.95 confidence)</td>
<td>0.991281</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>3.15941 ns [3.15465, 3.16432] (0.95 confidence)</td>
<td>3.16529 ns [3.16083, 3.16987] (0.95 confidence)</td>
<td>0.998143</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>3.15278 ns [3.14774, 3.15786] (0.95 confidence)</td>
<td>3.16511 ns [3.15941, 3.17092] (0.95 confidence)</td>
<td>0.996105</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>1.72702 ns [1.72497, 1.72908] (0.95 confidence)</td>
<td>1.73028 ns [1.72806, 1.73253] (0.95 confidence)</td>
<td>0.998116</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>3.13618 ns [3.13038, 3.14202] (0.95 confidence)</td>
<td>3.15069 ns [3.14423, 3.15718] (0.95 confidence)</td>
<td>0.995394</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;norm2&#95;sq</td>
<td>0.525516 ns [0.525067, 0.525989] (0.95 confidence)</td>
<td>0.525579 ns [0.52512, 0.526074] (0.95 confidence)</td>
<td>0.999881</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>3.96456 ns [3.95597, 3.97304] (0.95 confidence)</td>
<td>3.96582 ns [3.95827, 3.97306] (0.95 confidence)</td>
<td>0.999681</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>6.07768 ns [6.06931, 6.08689] (0.95 confidence)</td>
<td>6.06864 ns [6.06188, 6.07608] (0.95 confidence)</td>
<td>1.00149</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;det</td>
<td>14821.4 ns [14801.9, 14841.5] (0.95 confidence)</td>
<td>14332.6 ns [14304.4, 14370.7] (0.95 confidence)</td>
<td>1.0341</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>23.0114 ns [22.9765, 23.0444] (0.95 confidence)</td>
<td>23.5652 ns [23.5197, 23.6122] (0.95 confidence)</td>
<td>0.976499</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>23.0469 ns [23.0085, 23.0846] (0.95 confidence)</td>
<td>23.6266 ns [23.5801, 23.6708] (0.95 confidence)</td>
<td>0.975463</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>14870.8 ns [14804.5, 14943.4] (0.95 confidence)</td>
<td>14657.4 ns [14552.4, 14767.3] (0.95 confidence)</td>
<td>1.01456</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;dot</td>
<td>10.9006 ns [10.8856, 10.9145] (0.95 confidence)</td>
<td>11.2112 ns [11.1832, 11.2407] (0.95 confidence)</td>
<td>0.972296</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;inf&#95;norm</td>
<td>625.33 ns [624.026, 626.568] (0.95 confidence)</td>
<td>654.654 ns [653.294, 656.042] (0.95 confidence)</td>
<td>0.955207</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;ldlt</td>
<td>16833 ns [16817, 16851] (0.95 confidence)</td>
<td>17309.6 ns [17291.2, 17328.6] (0.95 confidence)</td>
<td>0.972467</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;ldlt&#95;solve</td>
<td>19139.7 ns [19108.8, 19169.9] (0.95 confidence)</td>
<td>19174.3 ns [19145.5, 19205.9] (0.95 confidence)</td>
<td>0.998198</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;lu</td>
<td>15292.4 ns [15272.5, 15315.9] (0.95 confidence)</td>
<td>15249.6 ns [15196.6, 15304.5] (0.95 confidence)</td>
<td>1.0028</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;norm2</td>
<td>57.501 ns [57.3692, 57.6357] (0.95 confidence)</td>
<td>60.4853 ns [60.025, 60.9673] (0.95 confidence)</td>
<td>0.950661</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>36.4755 ns [36.3788, 36.5696] (0.95 confidence)</td>
<td>39.5594 ns [39.4934, 39.6247] (0.95 confidence)</td>
<td>0.922042</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>36.3637 ns [36.2863, 36.4399] (0.95 confidence)</td>
<td>39.7902 ns [39.7365, 39.8411] (0.95 confidence)</td>
<td>0.913887</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>17.0586 ns [17.0227, 17.0937] (0.95 confidence)</td>
<td>18.1269 ns [18.0902, 18.1667] (0.95 confidence)</td>
<td>0.941062</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>32.4171 ns [32.3529, 32.4746] (0.95 confidence)</td>
<td>34.1155 ns [34.0563, 34.1785] (0.95 confidence)</td>
<td>0.950216</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;norm2&#95;sq</td>
<td>10.7461 ns [10.7231, 10.7681] (0.95 confidence)</td>
<td>11.07 ns [11.0499, 11.0909] (0.95 confidence)</td>
<td>0.97074</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>1019.11 ns [1017.26, 1021.09] (0.95 confidence)</td>
<td>1038.2 ns [1036.27, 1040.2] (0.95 confidence)</td>
<td>0.981612</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>2517.45 ns [2512.2, 2522.9] (0.95 confidence)</td>
<td>2633.77 ns [2625.96, 2642.35] (0.95 confidence)</td>
<td>0.955837</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;det</td>
<td>82.7954 ns [82.5989, 82.9921] (0.95 confidence)</td>
<td>82.8004 ns [82.5406, 83.0459] (0.95 confidence)</td>
<td>0.99994</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>1.31653 ns [1.31447, 1.31869] (0.95 confidence)</td>
<td>1.30663 ns [1.30537, 1.30793] (0.95 confidence)</td>
<td>1.00758</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;det&#95;from&#95;ldlt&#95;balanced&#95;range</td>
<td>8.67015 ns [8.65279, 8.68807] (0.95 confidence)</td>
<td>8.89678 ns [8.87977, 8.91454] (0.95 confidence)</td>
<td>0.974526</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>1.4712 ns [1.46835, 1.47426] (0.95 confidence)</td>
<td>1.47004 ns [1.46789, 1.47242] (0.95 confidence)</td>
<td>1.00079</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;det&#95;from&#95;lu&#95;balanced&#95;range</td>
<td>8.96025 ns [8.93841, 8.98285] (0.95 confidence)</td>
<td>11.2109 ns [11.1831, 11.2389] (0.95 confidence)</td>
<td>0.799245</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>93.5384 ns [92.9492, 94.1484] (0.95 confidence)</td>
<td>92.2086 ns [91.9807, 92.4447] (0.95 confidence)</td>
<td>1.01442</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;dot</td>
<td>0.929499 ns [0.926431, 0.932605] (0.95 confidence)</td>
<td>0.933639 ns [0.93087, 0.936421] (0.95 confidence)</td>
<td>0.995565</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;inf&#95;norm</td>
<td>8.43254 ns [8.42241, 8.44226] (0.95 confidence)</td>
<td>8.51129 ns [8.49464, 8.52764] (0.95 confidence)</td>
<td>0.990748</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;ldlt</td>
<td>90.8734 ns [90.5785, 91.1783] (0.95 confidence)</td>
<td>91.0659 ns [90.726, 91.4203] (0.95 confidence)</td>
<td>0.997886</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;ldlt&#95;ill&#95;conditioned</td>
<td>91.8546 ns [91.3084, 92.4254] (0.95 confidence)</td>
<td>92.7181 ns [92.2068, 93.2608] (0.95 confidence)</td>
<td>0.990687</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;ldlt&#95;solve</td>
<td>103.143 ns [102.875, 103.41] (0.95 confidence)</td>
<td>103.173 ns [102.903, 103.443] (0.95 confidence)</td>
<td>0.999711</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;lu</td>
<td>90.2174 ns [89.8798, 90.551] (0.95 confidence)</td>
<td>90.8863 ns [90.5822, 91.1941] (0.95 confidence)</td>
<td>0.99264</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;lu&#95;ill&#95;conditioned</td>
<td>89.4184 ns [89.1453, 89.6919] (0.95 confidence)</td>
<td>92.3313 ns [91.9995, 92.6728] (0.95 confidence)</td>
<td>0.968452</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;lu&#95;pivoting</td>
<td>94.4288 ns [94.0769, 94.812] (0.95 confidence)</td>
<td>96.1626 ns [95.8595, 96.4847] (0.95 confidence)</td>
<td>0.981971</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;norm2</td>
<td>4.59407 ns [4.57462, 4.61525] (0.95 confidence)</td>
<td>4.70306 ns [4.68112, 4.72642] (0.95 confidence)</td>
<td>0.976825</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>4.82523 ns [4.81579, 4.83474] (0.95 confidence)</td>
<td>4.8736 ns [4.86336, 4.88431] (0.95 confidence)</td>
<td>0.990075</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>4.8184 ns [4.80539, 4.83004] (0.95 confidence)</td>
<td>4.88617 ns [4.86974, 4.90377] (0.95 confidence)</td>
<td>0.986131</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>2.80707 ns [2.74298, 2.87062] (0.95 confidence)</td>
<td>2.76914 ns [2.70565, 2.83476] (0.95 confidence)</td>
<td>1.0137</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>4.20888 ns [4.20034, 4.21691] (0.95 confidence)</td>
<td>4.27181 ns [4.2627, 4.28124] (0.95 confidence)</td>
<td>0.98527</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;norm2&#95;sq</td>
<td>0.662231 ns [0.660744, 0.663713] (0.95 confidence)</td>
<td>0.658078 ns [0.656342, 0.659836] (0.95 confidence)</td>
<td>1.00631</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>8.09859 ns [8.08944, 8.10837] (0.95 confidence)</td>
<td>8.11805 ns [8.1047, 8.13282] (0.95 confidence)</td>
<td>0.997603</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>13.8087 ns [13.7351, 13.8856] (0.95 confidence)</td>
<td>13.8011 ns [13.6834, 13.9313] (0.95 confidence)</td>
<td>1.00056</td>
</tr>
<tr>
<td>gram/1x2/dependent/la&#95;stack</td>
<td>0.686554 ns [0.685225, 0.687957] (0.95 confidence)</td>
<td>0.694393 ns [0.693023, 0.695861] (0.95 confidence)</td>
<td>0.988711</td>
</tr>
<tr>
<td>gram/1x2/mixed&#95;scale/la&#95;stack</td>
<td>0.686637 ns [0.685284, 0.688061] (0.95 confidence)</td>
<td>0.696152 ns [0.694301, 0.698227] (0.95 confidence)</td>
<td>0.986333</td>
</tr>
<tr>
<td>gram/1x2/near&#95;dependent/la&#95;stack</td>
<td>0.686876 ns [0.685533, 0.688295] (0.95 confidence)</td>
<td>0.697503 ns [0.69576, 0.699339] (0.95 confidence)</td>
<td>0.984765</td>
</tr>
<tr>
<td>gram/1x2/orthogonal/la&#95;stack</td>
<td>0.686984 ns [0.685506, 0.688595] (0.95 confidence)</td>
<td>0.694483 ns [0.693102, 0.695942] (0.95 confidence)</td>
<td>0.989202</td>
</tr>
<tr>
<td>gram/2x2/dependent/la&#95;stack</td>
<td>1.30757 ns [1.30165, 1.31338] (0.95 confidence)</td>
<td>1.32323 ns [1.31753, 1.32889] (0.95 confidence)</td>
<td>0.98816</td>
</tr>
<tr>
<td>gram/2x2/mixed&#95;scale/la&#95;stack</td>
<td>1.32194 ns [1.31388, 1.33062] (0.95 confidence)</td>
<td>1.32392 ns [1.31837, 1.32941] (0.95 confidence)</td>
<td>0.9985</td>
</tr>
<tr>
<td>gram/2x2/near&#95;dependent/la&#95;stack</td>
<td>1.30769 ns [1.30213, 1.31318] (0.95 confidence)</td>
<td>1.3249 ns [1.31905, 1.33073] (0.95 confidence)</td>
<td>0.987012</td>
</tr>
<tr>
<td>gram/2x2/orthogonal/la&#95;stack</td>
<td>1.31071 ns [1.30455, 1.31706] (0.95 confidence)</td>
<td>1.32868 ns [1.32148, 1.3363] (0.95 confidence)</td>
<td>0.986477</td>
</tr>
<tr>
<td>gram/2x3/dependent/la&#95;stack</td>
<td>1.27534 ns [1.26836, 1.28647] (0.95 confidence)</td>
<td>1.30911 ns [1.30015, 1.32229] (0.95 confidence)</td>
<td>0.974204</td>
</tr>
<tr>
<td>gram/2x3/mixed&#95;scale/la&#95;stack</td>
<td>1.26807 ns [1.26643, 1.26979] (0.95 confidence)</td>
<td>1.30854 ns [1.30549, 1.31198] (0.95 confidence)</td>
<td>0.969066</td>
</tr>
<tr>
<td>gram/2x3/near&#95;dependent/la&#95;stack</td>
<td>1.26778 ns [1.26604, 1.26966] (0.95 confidence)</td>
<td>1.2922 ns [1.2901, 1.29438] (0.95 confidence)</td>
<td>0.981102</td>
</tr>
<tr>
<td>gram/2x3/orthogonal/la&#95;stack</td>
<td>1.2682 ns [1.26642, 1.27006] (0.95 confidence)</td>
<td>1.29317 ns [1.28575, 1.30517] (0.95 confidence)</td>
<td>0.980687</td>
</tr>
<tr>
<td>gram/3x3/dependent/la&#95;stack</td>
<td>3.00633 ns [3.00251, 3.01056] (0.95 confidence)</td>
<td>3.02683 ns [3.02342, 3.03043] (0.95 confidence)</td>
<td>0.993227</td>
</tr>
<tr>
<td>gram/3x3/mixed&#95;scale/la&#95;stack</td>
<td>3.01606 ns [3.00986, 3.02279] (0.95 confidence)</td>
<td>3.0317 ns [3.02806, 3.03558] (0.95 confidence)</td>
<td>0.994841</td>
</tr>
<tr>
<td>gram/3x3/near&#95;dependent/la&#95;stack</td>
<td>3.00439 ns [2.99987, 3.00967] (0.95 confidence)</td>
<td>3.02677 ns [3.02305, 3.03085] (0.95 confidence)</td>
<td>0.992604</td>
</tr>
<tr>
<td>gram/3x3/orthogonal/la&#95;stack</td>
<td>3.00323 ns [2.99983, 3.00656] (0.95 confidence)</td>
<td>3.03405 ns [3.02871, 3.04029] (0.95 confidence)</td>
<td>0.989839</td>
</tr>
<tr>
<td>gram/3x4/dependent/la&#95;stack</td>
<td>3.05165 ns [3.04839, 3.05519] (0.95 confidence)</td>
<td>4.32553 ns [4.14038, 4.53897] (0.95 confidence)</td>
<td>0.705498</td>
</tr>
<tr>
<td>gram/3x4/mixed&#95;scale/la&#95;stack</td>
<td>3.05325 ns [3.04831, 3.05959] (0.95 confidence)</td>
<td>3.1599 ns [3.14361, 3.17778] (0.95 confidence)</td>
<td>0.966249</td>
</tr>
<tr>
<td>gram/3x4/near&#95;dependent/la&#95;stack</td>
<td>3.04916 ns [3.04481, 3.05263] (0.95 confidence)</td>
<td>3.78715 ns [3.71379, 3.87394] (0.95 confidence)</td>
<td>0.805134</td>
</tr>
<tr>
<td>gram/3x4/orthogonal/la&#95;stack</td>
<td>3.05082 ns [3.04768, 3.0542] (0.95 confidence)</td>
<td>3.08891 ns [3.08478, 3.09363] (0.95 confidence)</td>
<td>0.987666</td>
</tr>
<tr>
<td>gram/4x4/dependent/la&#95;stack</td>
<td>7.58099 ns [7.53587, 7.62571] (0.95 confidence)</td>
<td>7.72053 ns [7.67555, 7.76514] (0.95 confidence)</td>
<td>0.981926</td>
</tr>
<tr>
<td>gram/4x4/mixed&#95;scale/la&#95;stack</td>
<td>7.59536 ns [7.54908, 7.64121] (0.95 confidence)</td>
<td>8.94192 ns [8.71522, 9.18554] (0.95 confidence)</td>
<td>0.849411</td>
</tr>
<tr>
<td>gram/4x4/near&#95;dependent/la&#95;stack</td>
<td>7.57485 ns [7.53054, 7.61929] (0.95 confidence)</td>
<td>7.84145 ns [7.78823, 7.89686] (0.95 confidence)</td>
<td>0.966001</td>
</tr>
<tr>
<td>gram/4x4/orthogonal/la&#95;stack</td>
<td>7.60858 ns [7.562, 7.655] (0.95 confidence)</td>
<td>7.75992 ns [7.7095, 7.81125] (0.95 confidence)</td>
<td>0.980498</td>
</tr>
<tr>
<td>gram/4x5/dependent/la&#95;stack</td>
<td>9.48698 ns [9.46718, 9.50804] (0.95 confidence)</td>
<td>11.7891 ns [11.3922, 12.2443] (0.95 confidence)</td>
<td>0.804722</td>
</tr>
<tr>
<td>gram/4x5/mixed&#95;scale/la&#95;stack</td>
<td>9.43186 ns [9.40679, 9.45623] (0.95 confidence)</td>
<td>13.0738 ns [12.4995, 13.7294] (0.95 confidence)</td>
<td>0.721433</td>
</tr>
<tr>
<td>gram/4x5/near&#95;dependent/la&#95;stack</td>
<td>9.47989 ns [9.46575, 9.49384] (0.95 confidence)</td>
<td>13.7036 ns [12.8616, 14.6524] (0.95 confidence)</td>
<td>0.691779</td>
</tr>
<tr>
<td>gram/4x5/orthogonal/la&#95;stack</td>
<td>9.40964 ns [9.38368, 9.43475] (0.95 confidence)</td>
<td>10.2945 ns [10.1004, 10.5045] (0.95 confidence)</td>
<td>0.914048</td>
</tr>
<tr>
<td>gram/5x5/dependent/la&#95;stack</td>
<td>18.423 ns [18.3892, 18.458] (0.95 confidence)</td>
<td>18.5407 ns [18.4842, 18.5974] (0.95 confidence)</td>
<td>0.993649</td>
</tr>
<tr>
<td>gram/5x5/mixed&#95;scale/la&#95;stack</td>
<td>18.4884 ns [18.4101, 18.5601] (0.95 confidence)</td>
<td>18.5925 ns [18.5013, 18.6888] (0.95 confidence)</td>
<td>0.994402</td>
</tr>
<tr>
<td>gram/5x5/near&#95;dependent/la&#95;stack</td>
<td>18.4414 ns [18.3577, 18.5195] (0.95 confidence)</td>
<td>18.6126 ns [18.5459, 18.6794] (0.95 confidence)</td>
<td>0.990804</td>
</tr>
<tr>
<td>gram/5x5/orthogonal/la&#95;stack</td>
<td>18.497 ns [18.4193, 18.5668] (0.95 confidence)</td>
<td>18.9381 ns [18.7742, 19.132] (0.95 confidence)</td>
<td>0.976709</td>
</tr>
<tr>
<td>gram/5x6/dependent/la&#95;stack</td>
<td>19.283 ns [19.2409, 19.3237] (0.95 confidence)</td>
<td>29.6311 ns [28.5998, 30.9413] (0.95 confidence)</td>
<td>0.650768</td>
</tr>
<tr>
<td>gram/5x6/mixed&#95;scale/la&#95;stack</td>
<td>19.2153 ns [19.1626, 19.2646] (0.95 confidence)</td>
<td>19.4247 ns [19.3678, 19.4812] (0.95 confidence)</td>
<td>0.989224</td>
</tr>
<tr>
<td>gram/5x6/near&#95;dependent/la&#95;stack</td>
<td>19.0856 ns [19.0248, 19.1437] (0.95 confidence)</td>
<td>26.4399 ns [25.9339, 27.026] (0.95 confidence)</td>
<td>0.721851</td>
</tr>
<tr>
<td>gram/5x6/orthogonal/la&#95;stack</td>
<td>19.1951 ns [19.1491, 19.239] (0.95 confidence)</td>
<td>21.4264 ns [21.3004, 21.5507] (0.95 confidence)</td>
<td>0.895863</td>
</tr>
<tr>
<td>gram/6x6/dependent/la&#95;stack</td>
<td>29.3586 ns [29.1769, 29.5422] (0.95 confidence)</td>
<td>39.6842 ns [38.5064, 40.9225] (0.95 confidence)</td>
<td>0.739806</td>
</tr>
<tr>
<td>gram/6x6/mixed&#95;scale/la&#95;stack</td>
<td>29.3076 ns [29.1188, 29.502] (0.95 confidence)</td>
<td>42.0471 ns [40.6239, 43.6641] (0.95 confidence)</td>
<td>0.697018</td>
</tr>
<tr>
<td>gram/6x6/near&#95;dependent/la&#95;stack</td>
<td>29.4143 ns [29.2518, 29.5826] (0.95 confidence)</td>
<td>40.836 ns [39.5378, 42.3194] (0.95 confidence)</td>
<td>0.720304</td>
</tr>
<tr>
<td>gram/6x6/orthogonal/la&#95;stack</td>
<td>29.6191 ns [29.4263, 29.8146] (0.95 confidence)</td>
<td>39.0287 ns [37.5893, 40.6234] (0.95 confidence)</td>
<td>0.758904</td>
</tr>
<tr>
<td>gram/6x7/dependent/la&#95;stack</td>
<td>32.1442 ns [32.0192, 32.268] (0.95 confidence)</td>
<td>32.0059 ns [31.9068, 32.1055] (0.95 confidence)</td>
<td>1.00432</td>
</tr>
<tr>
<td>gram/6x7/mixed&#95;scale/la&#95;stack</td>
<td>32.2011 ns [32.0698, 32.3306] (0.95 confidence)</td>
<td>32.2415 ns [32.1377, 32.3453] (0.95 confidence)</td>
<td>0.998749</td>
</tr>
<tr>
<td>gram/6x7/near&#95;dependent/la&#95;stack</td>
<td>32.2633 ns [32.1733, 32.3536] (0.95 confidence)</td>
<td>32.0789 ns [31.9549, 32.1998] (0.95 confidence)</td>
<td>1.00575</td>
</tr>
<tr>
<td>gram/6x7/orthogonal/la&#95;stack</td>
<td>32.1922 ns [32.0496, 32.3366] (0.95 confidence)</td>
<td>32.2025 ns [32.068, 32.3365] (0.95 confidence)</td>
<td>0.999681</td>
</tr>
<tr>
<td>gram/7x7/dependent/la&#95;stack</td>
<td>39.2631 ns [39.1028, 39.4409] (0.95 confidence)</td>
<td>72.6031 ns [70.5955, 74.8363] (0.95 confidence)</td>
<td>0.540791</td>
</tr>
<tr>
<td>gram/7x7/mixed&#95;scale/la&#95;stack</td>
<td>39.8351 ns [39.5311, 40.1773] (0.95 confidence)</td>
<td>46.3954 ns [45.555, 47.3694] (0.95 confidence)</td>
<td>0.858599</td>
</tr>
<tr>
<td>gram/7x7/near&#95;dependent/la&#95;stack</td>
<td>39.8875 ns [39.5095, 40.2992] (0.95 confidence)</td>
<td>42.5067 ns [42.2216, 42.7983] (0.95 confidence)</td>
<td>0.938382</td>
</tr>
<tr>
<td>gram/7x7/orthogonal/la&#95;stack</td>
<td>39.6303 ns [39.3504, 39.9531] (0.95 confidence)</td>
<td>57.2046 ns [55.2902, 59.1841] (0.95 confidence)</td>
<td>0.692782</td>
</tr>
<tr>
<td>gram/7x8/dependent/la&#95;stack</td>
<td>40.8438 ns [40.5009, 41.2292] (0.95 confidence)</td>
<td>40.708 ns [40.5916, 40.8396] (0.95 confidence)</td>
<td>1.00334</td>
</tr>
<tr>
<td>gram/7x8/mixed&#95;scale/la&#95;stack</td>
<td>41.5158 ns [41.0854, 41.9744] (0.95 confidence)</td>
<td>42.3546 ns [41.9865, 42.7281] (0.95 confidence)</td>
<td>0.980194</td>
</tr>
<tr>
<td>gram/7x8/near&#95;dependent/la&#95;stack</td>
<td>40.6504 ns [40.3654, 40.9739] (0.95 confidence)</td>
<td>41.7561 ns [41.4546, 42.0732] (0.95 confidence)</td>
<td>0.97352</td>
</tr>
<tr>
<td>gram/7x8/orthogonal/la&#95;stack</td>
<td>41.9241 ns [41.388, 42.4794] (0.95 confidence)</td>
<td>41.7167 ns [41.4869, 41.9623] (0.95 confidence)</td>
<td>1.00497</td>
</tr>
<tr>
<td>gram/8x8/dependent/la&#95;stack</td>
<td>51.9391 ns [51.8348, 52.0555] (0.95 confidence)</td>
<td>52.7757 ns [52.5576, 52.9994] (0.95 confidence)</td>
<td>0.984148</td>
</tr>
<tr>
<td>gram/8x8/mixed&#95;scale/la&#95;stack</td>
<td>51.8843 ns [51.8238, 51.956] (0.95 confidence)</td>
<td>51.9988 ns [51.8512, 52.1657] (0.95 confidence)</td>
<td>0.997799</td>
</tr>
<tr>
<td>gram/8x8/near&#95;dependent/la&#95;stack</td>
<td>52.0055 ns [51.8718, 52.1541] (0.95 confidence)</td>
<td>51.8606 ns [51.7306, 52.0078] (0.95 confidence)</td>
<td>1.00279</td>
</tr>
<tr>
<td>gram/8x8/orthogonal/la&#95;stack</td>
<td>52.2229 ns [52.0611, 52.3935] (0.95 confidence)</td>
<td>51.9693 ns [51.8208, 52.1325] (0.95 confidence)</td>
<td>1.00488</td>
</tr>
<tr>
<td>interval&#95;d3/dense&#95;inconclusive</td>
<td>41.6783 ns [41.5836, 41.7803] (0.95 confidence)</td>
<td>43.6174 ns [43.3436, 43.9154] (0.95 confidence)</td>
<td>0.955543</td>
</tr>
<tr>
<td>interval&#95;d3/dense&#95;sign</td>
<td>37.7162 ns [37.4009, 38.0312] (0.95 confidence)</td>
<td>39.3183 ns [38.915, 39.7136] (0.95 confidence)</td>
<td>0.959253</td>
</tr>
<tr>
<td>interval&#95;d3/lifted&#95;assembly</td>
<td>8.36173 ns [8.35232, 8.37184] (0.95 confidence)</td>
<td>8.41746 ns [8.40484, 8.43187] (0.95 confidence)</td>
<td>0.993379</td>
</tr>
<tr>
<td>interval&#95;d3/lifted&#95;assembly&#95;and&#95;sign</td>
<td>37.05 ns [36.9169, 37.1855] (0.95 confidence)</td>
<td>37.545 ns [37.3393, 37.7584] (0.95 confidence)</td>
<td>0.986817</td>
</tr>
<tr>
<td>interval&#95;d3/sparse&#95;sign</td>
<td>31.4658 ns [31.3054, 31.616] (0.95 confidence)</td>
<td>31.0511 ns [30.8608, 31.2459] (0.95 confidence)</td>
<td>1.01335</td>
</tr>
<tr>
<td>interval&#95;d4/dense&#95;inconclusive</td>
<td>100.222 ns [99.6724, 100.792] (0.95 confidence)</td>
<td>106.362 ns [105.409, 107.355] (0.95 confidence)</td>
<td>0.94227</td>
</tr>
<tr>
<td>interval&#95;d4/dense&#95;sign</td>
<td>100.924 ns [100.238, 101.611] (0.95 confidence)</td>
<td>109.425 ns [108.572, 110.264] (0.95 confidence)</td>
<td>0.922314</td>
</tr>
<tr>
<td>interval&#95;d4/lifted&#95;assembly</td>
<td>13.3071 ns [13.2941, 13.322] (0.95 confidence)</td>
<td>13.611 ns [13.5266, 13.6969] (0.95 confidence)</td>
<td>0.977672</td>
</tr>
<tr>
<td>interval&#95;d4/lifted&#95;assembly&#95;and&#95;sign</td>
<td>82.3934 ns [81.9953, 82.7892] (0.95 confidence)</td>
<td>85.271 ns [84.8317, 85.722] (0.95 confidence)</td>
<td>0.966254</td>
</tr>
<tr>
<td>interval&#95;d4/sparse&#95;sign</td>
<td>68.2494 ns [67.8808, 68.6229] (0.95 confidence)</td>
<td>70.9735 ns [70.4201, 71.5601] (0.95 confidence)</td>
<td>0.961617</td>
</tr>
<tr>
<td>interval&#95;d5/dense&#95;inconclusive</td>
<td>237.366 ns [236.688, 238.05] (0.95 confidence)</td>
<td>238.519 ns [237.723, 239.337] (0.95 confidence)</td>
<td>0.995169</td>
</tr>
<tr>
<td>interval&#95;d5/dense&#95;sign</td>
<td>223.692 ns [222.739, 224.688] (0.95 confidence)</td>
<td>234.131 ns [231.985, 236.383] (0.95 confidence)</td>
<td>0.955414</td>
</tr>
<tr>
<td>interval&#95;d5/lifted&#95;assembly</td>
<td>21.6842 ns [21.6275, 21.7408] (0.95 confidence)</td>
<td>21.9298 ns [21.8634, 21.9968] (0.95 confidence)</td>
<td>0.988803</td>
</tr>
<tr>
<td>interval&#95;d5/lifted&#95;assembly&#95;and&#95;sign</td>
<td>186.842 ns [186.28, 187.416] (0.95 confidence)</td>
<td>198.179 ns [197.188, 199.212] (0.95 confidence)</td>
<td>0.942791</td>
</tr>
<tr>
<td>interval&#95;d5/sparse&#95;sign</td>
<td>153.751 ns [153.243, 154.251] (0.95 confidence)</td>
<td>166.712 ns [165.379, 168.029] (0.95 confidence)</td>
<td>0.922257</td>
</tr>
<tr>
<td>interval&#95;d6/dense&#95;inconclusive</td>
<td>535.598 ns [533.961, 537.269] (0.95 confidence)</td>
<td>545.482 ns [541.712, 549.333] (0.95 confidence)</td>
<td>0.98188</td>
</tr>
<tr>
<td>interval&#95;d6/dense&#95;sign</td>
<td>524.469 ns [522.532, 526.387] (0.95 confidence)</td>
<td>527.15 ns [525.032, 529.463] (0.95 confidence)</td>
<td>0.994914</td>
</tr>
<tr>
<td>interval&#95;d6/lifted&#95;assembly</td>
<td>30.9939 ns [30.9459, 31.0405] (0.95 confidence)</td>
<td>31.2 ns [31.1463, 31.2555] (0.95 confidence)</td>
<td>0.993394</td>
</tr>
<tr>
<td>interval&#95;d6/lifted&#95;assembly&#95;and&#95;sign</td>
<td>387.07 ns [385.205, 388.984] (0.95 confidence)</td>
<td>397.046 ns [395.412, 398.72] (0.95 confidence)</td>
<td>0.974876</td>
</tr>
<tr>
<td>interval&#95;d6/sparse&#95;sign</td>
<td>364.96 ns [363.744, 366.186] (0.95 confidence)</td>
<td>376.408 ns [373.995, 378.916] (0.95 confidence)</td>
<td>0.969586</td>
</tr>
<tr>
<td>interval&#95;d7/dense&#95;inconclusive</td>
<td>1163.92 ns [1156.87, 1170.72] (0.95 confidence)</td>
<td>1156.59 ns [1151.03, 1162.1] (0.95 confidence)</td>
<td>1.00634</td>
</tr>
<tr>
<td>interval&#95;d7/dense&#95;sign</td>
<td>1136.73 ns [1131.77, 1141.69] (0.95 confidence)</td>
<td>1126.29 ns [1119.88, 1133.14] (0.95 confidence)</td>
<td>1.00927</td>
</tr>
<tr>
<td>interval&#95;d7/lifted&#95;assembly</td>
<td>37.4663 ns [37.3994, 37.5328] (0.95 confidence)</td>
<td>38.0795 ns [37.9823, 38.1792] (0.95 confidence)</td>
<td>0.983896</td>
</tr>
<tr>
<td>interval&#95;d7/lifted&#95;assembly&#95;and&#95;sign</td>
<td>779.866 ns [777.773, 782.078] (0.95 confidence)</td>
<td>806.783 ns [804.027, 809.674] (0.95 confidence)</td>
<td>0.966638</td>
</tr>
<tr>
<td>interval&#95;d7/sparse&#95;sign</td>
<td>752.181 ns [750.082, 754.375] (0.95 confidence)</td>
<td>770.621 ns [767.991, 773.416] (0.95 confidence)</td>
<td>0.976072</td>
</tr>
<tr>
<td>interval&#95;det&#95;sign/d4&#95;conclusive&#95;lifted</td>
<td>69.034 ns [68.6, 69.4653] (0.95 confidence)</td>
<td>67.6038 ns [67.2468, 67.9929] (0.95 confidence)</td>
<td>1.02115</td>
</tr>
<tr>
<td>interval&#95;det&#95;sign/d4&#95;inconclusive&#95;lifted</td>
<td>68.2205 ns [67.8361, 68.605] (0.95 confidence)</td>
<td>66.1255 ns [65.7481, 66.5309] (0.95 confidence)</td>
<td>1.03168</td>
</tr>
<tr>
<td>interval&#95;det&#95;sign/d7&#95;conclusive&#95;lifted</td>
<td>751.326 ns [748.829, 754.008] (0.95 confidence)</td>
<td>760.285 ns [757.01, 763.824] (0.95 confidence)</td>
<td>0.988216</td>
</tr>
<tr>
<td>interval&#95;scalar/cancellation/add</td>
<td>2.12984 ns [2.12651, 2.13302] (0.95 confidence)</td>
<td>2.12215 ns [2.11914, 2.12505] (0.95 confidence)</td>
<td>1.00362</td>
</tr>
<tr>
<td>interval&#95;scalar/cancellation/multiply</td>
<td>1.88453 ns [1.88076, 1.88832] (0.95 confidence)</td>
<td>1.88495 ns [1.88142, 1.88848] (0.95 confidence)</td>
<td>0.999774</td>
</tr>
<tr>
<td>interval&#95;scalar/cancellation/square</td>
<td>1.8727 ns [1.86965, 1.87572] (0.95 confidence)</td>
<td>1.87396 ns [1.86896, 1.87928] (0.95 confidence)</td>
<td>0.999327</td>
</tr>
<tr>
<td>interval&#95;scalar/exact/add</td>
<td>2.12653 ns [2.12335, 2.12945] (0.95 confidence)</td>
<td>2.09524 ns [2.09207, 2.09843] (0.95 confidence)</td>
<td>1.01493</td>
</tr>
<tr>
<td>interval&#95;scalar/exact/multiply</td>
<td>2.09355 ns [2.08437, 2.10219] (0.95 confidence)</td>
<td>2.07472 ns [2.06259, 2.08647] (0.95 confidence)</td>
<td>1.00907</td>
</tr>
<tr>
<td>interval&#95;scalar/exact/square</td>
<td>1.87366 ns [1.87099, 1.87632] (0.95 confidence)</td>
<td>1.8473 ns [1.84388, 1.85049] (0.95 confidence)</td>
<td>1.01427</td>
</tr>
<tr>
<td>interval&#95;scalar/point/add</td>
<td>2.17759 ns [2.17579, 2.17953] (0.95 confidence)</td>
<td>2.163 ns [2.15885, 2.16689] (0.95 confidence)</td>
<td>1.00675</td>
</tr>
<tr>
<td>interval&#95;scalar/point/multiply</td>
<td>2.60903 ns [2.59706, 2.62204] (0.95 confidence)</td>
<td>2.59323 ns [2.58385, 2.6033] (0.95 confidence)</td>
<td>1.00609</td>
</tr>
<tr>
<td>interval&#95;scalar/point/square</td>
<td>2.12777 ns [2.12496, 2.13049] (0.95 confidence)</td>
<td>2.08336 ns [2.07964, 2.08703] (0.95 confidence)</td>
<td>1.02131</td>
</tr>
<tr>
<td>interval&#95;scalar/subnormal/add</td>
<td>2.41479 ns [2.40801, 2.42275] (0.95 confidence)</td>
<td>2.40145 ns [2.39778, 2.40511] (0.95 confidence)</td>
<td>1.00555</td>
</tr>
<tr>
<td>interval&#95;scalar/subnormal/multiply</td>
<td>2.84981 ns [2.83715, 2.86344] (0.95 confidence)</td>
<td>2.92398 ns [2.90637, 2.94176] (0.95 confidence)</td>
<td>0.974632</td>
</tr>
<tr>
<td>interval&#95;scalar/subnormal/square</td>
<td>2.37699 ns [2.37459, 2.37974] (0.95 confidence)</td>
<td>2.37132 ns [2.36882, 2.37396] (0.95 confidence)</td>
<td>1.00239</td>
</tr>
<tr>
<td>interval&#95;scalar/zero/add</td>
<td>0.941051 ns [0.939642, 0.942748] (0.95 confidence)</td>
<td>0.933392 ns [0.932044, 0.93471] (0.95 confidence)</td>
<td>1.00821</td>
</tr>
<tr>
<td>interval&#95;scalar/zero/multiply</td>
<td>1.76802 ns [1.76541, 1.77067] (0.95 confidence)</td>
<td>1.75359 ns [1.75073, 1.75639] (0.95 confidence)</td>
<td>1.00823</td>
</tr>
<tr>
<td>interval&#95;scalar/zero/square</td>
<td>0.881135 ns [0.879756, 0.882269] (0.95 confidence)</td>
<td>0.879828 ns [0.879159, 0.880493] (0.95 confidence)</td>
<td>1.00149</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;batch/origin/construction</td>
<td>2.07257 ns [2.06937, 2.07575] (0.95 confidence)</td>
<td>2.08368 ns [2.07967, 2.08758] (0.95 confidence)</td>
<td>0.99467</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;batch/origin/prepared</td>
<td>13.8301 ns [13.8161, 13.8461] (0.95 confidence)</td>
<td>13.8611 ns [13.8468, 13.8768] (0.95 confidence)</td>
<td>0.997757</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;batch/origin/with&#95;construction</td>
<td>15.3626 ns [15.3421, 15.3866] (0.95 confidence)</td>
<td>15.3538 ns [15.3394, 15.3701] (0.95 confidence)</td>
<td>1.00057</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;batch/translated/construction</td>
<td>2.07343 ns [2.07041, 2.07627] (0.95 confidence)</td>
<td>2.07408 ns [2.06947, 2.07828] (0.95 confidence)</td>
<td>0.999683</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;batch/translated/prepared</td>
<td>41.072 ns [40.9562, 41.1973] (0.95 confidence)</td>
<td>41.1985 ns [41.0954, 41.306] (0.95 confidence)</td>
<td>0.996931</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;batch/translated/with&#95;construction</td>
<td>42.8472 ns [42.7219, 42.9813] (0.95 confidence)</td>
<td>43.2018 ns [43.0685, 43.34] (0.95 confidence)</td>
<td>0.991791</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/cancellation/difference</td>
<td>10.8479 ns [10.765, 10.9452] (0.95 confidence)</td>
<td>10.9383 ns [10.8492, 11.0374] (0.95 confidence)</td>
<td>0.991728</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/cancellation/difference&#95;endpoints</td>
<td>11.892 ns [11.7547, 12.0478] (0.95 confidence)</td>
<td>11.6706 ns [11.5966, 11.7621] (0.95 confidence)</td>
<td>1.01897</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/cancellation/dot</td>
<td>5.92874 ns [5.92112, 5.93669] (0.95 confidence)</td>
<td>5.98938 ns [5.96624, 6.01554] (0.95 confidence)</td>
<td>0.989876</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/cancellation/dot&#95;endpoints</td>
<td>7.04162 ns [7.02504, 7.05866] (0.95 confidence)</td>
<td>7.17826 ns [7.15518, 7.20157] (0.95 confidence)</td>
<td>0.980964</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense/difference</td>
<td>10.4051 ns [10.3525, 10.4703] (0.95 confidence)</td>
<td>10.5554 ns [10.4655, 10.6549] (0.95 confidence)</td>
<td>0.985761</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense/difference&#95;endpoints</td>
<td>12.7901 ns [12.6929, 12.8989] (0.95 confidence)</td>
<td>12.9066 ns [12.7995, 13.0206] (0.95 confidence)</td>
<td>0.990968</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense/dot</td>
<td>5.52712 ns [5.51519, 5.53997] (0.95 confidence)</td>
<td>5.56744 ns [5.55265, 5.58302] (0.95 confidence)</td>
<td>0.992757</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense/dot&#95;endpoints</td>
<td>6.5175 ns [6.49582, 6.54077] (0.95 confidence)</td>
<td>6.67251 ns [6.63847, 6.70904] (0.95 confidence)</td>
<td>0.976769</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense&#95;inexact/difference</td>
<td>10.4593 ns [10.3928, 10.5392] (0.95 confidence)</td>
<td>10.6023 ns [10.4903, 10.7264] (0.95 confidence)</td>
<td>0.986512</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense&#95;inexact/difference&#95;endpoints</td>
<td>11.7852 ns [11.6915, 11.8897] (0.95 confidence)</td>
<td>12.0621 ns [11.9377, 12.1927] (0.95 confidence)</td>
<td>0.977043</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense&#95;inexact/dot</td>
<td>5.55083 ns [5.53588, 5.5664] (0.95 confidence)</td>
<td>5.59001 ns [5.57896, 5.60129] (0.95 confidence)</td>
<td>0.992991</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense&#95;inexact/dot&#95;endpoints</td>
<td>6.72995 ns [6.69989, 6.76075] (0.95 confidence)</td>
<td>6.60576 ns [6.58192, 6.63138] (0.95 confidence)</td>
<td>1.0188</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/overflow/dot</td>
<td>1.07179 ns [1.07004, 1.0736] (0.95 confidence)</td>
<td>1.07177 ns [1.06959, 1.07407] (0.95 confidence)</td>
<td>1.00002</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/sparse/difference</td>
<td>4.75658 ns [4.74197, 4.77571] (0.95 confidence)</td>
<td>4.7523 ns [4.73801, 4.77113] (0.95 confidence)</td>
<td>1.0009</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/sparse/difference&#95;endpoints</td>
<td>5.59961 ns [5.58371, 5.61565] (0.95 confidence)</td>
<td>5.61839 ns [5.60677, 5.63028] (0.95 confidence)</td>
<td>0.996658</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/sparse/dot</td>
<td>2.22702 ns [2.22512, 2.22923] (0.95 confidence)</td>
<td>2.23463 ns [2.23183, 2.23773] (0.95 confidence)</td>
<td>0.996595</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/sparse/dot&#95;endpoints</td>
<td>2.22835 ns [2.22645, 2.2305] (0.95 confidence)</td>
<td>2.23039 ns [2.22853, 2.23259] (0.95 confidence)</td>
<td>0.999082</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/underflow/difference</td>
<td>4.464 ns [4.46012, 4.46877] (0.95 confidence)</td>
<td>4.46454 ns [4.46077, 4.46901] (0.95 confidence)</td>
<td>0.99988</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/underflow/difference&#95;endpoints</td>
<td>4.60958 ns [4.60587, 4.61396] (0.95 confidence)</td>
<td>4.6162 ns [4.60655, 4.62833] (0.95 confidence)</td>
<td>0.998566</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/underflow/dot</td>
<td>1.48357 ns [1.48211, 1.48516] (0.95 confidence)</td>
<td>1.48487 ns [1.48384, 1.48608] (0.95 confidence)</td>
<td>0.999122</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/underflow/dot&#95;endpoints</td>
<td>1.48505 ns [1.48384, 1.48637] (0.95 confidence)</td>
<td>1.48626 ns [1.48491, 1.48775] (0.95 confidence)</td>
<td>0.999182</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;batch/origin/construction</td>
<td>2.26741 ns [2.25885, 2.27566] (0.95 confidence)</td>
<td>2.26282 ns [2.25443, 2.27109] (0.95 confidence)</td>
<td>1.00202</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;batch/origin/prepared</td>
<td>23.5348 ns [23.4362, 23.6458] (0.95 confidence)</td>
<td>24.0021 ns [23.9288, 24.0812] (0.95 confidence)</td>
<td>0.980531</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;batch/origin/with&#95;construction</td>
<td>27.9835 ns [27.88, 28.0907] (0.95 confidence)</td>
<td>27.6646 ns [27.6119, 27.7207] (0.95 confidence)</td>
<td>1.01153</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;batch/translated/construction</td>
<td>2.27584 ns [2.2679, 2.28361] (0.95 confidence)</td>
<td>2.26378 ns [2.25523, 2.2723] (0.95 confidence)</td>
<td>1.00533</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;batch/translated/prepared</td>
<td>68.5552 ns [68.4729, 68.6347] (0.95 confidence)</td>
<td>70.3717 ns [70.0274, 70.7248] (0.95 confidence)</td>
<td>0.974188</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;batch/translated/with&#95;construction</td>
<td>72.0859 ns [71.7275, 72.4814] (0.95 confidence)</td>
<td>72.1918 ns [71.9316, 72.4702] (0.95 confidence)</td>
<td>0.998533</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/cancellation/difference</td>
<td>11.2361 ns [11.1347, 11.3445] (0.95 confidence)</td>
<td>11.3363 ns [11.2221, 11.4608] (0.95 confidence)</td>
<td>0.991166</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/cancellation/difference&#95;endpoints</td>
<td>12.2106 ns [12.0995, 12.334] (0.95 confidence)</td>
<td>12.0787 ns [11.9825, 12.1941] (0.95 confidence)</td>
<td>1.01092</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/cancellation/dot</td>
<td>6.04934 ns [6.03471, 6.07088] (0.95 confidence)</td>
<td>6.25588 ns [6.20524, 6.31167] (0.95 confidence)</td>
<td>0.966986</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/cancellation/dot&#95;endpoints</td>
<td>6.68969 ns [6.64259, 6.74362] (0.95 confidence)</td>
<td>6.7322 ns [6.69936, 6.76837] (0.95 confidence)</td>
<td>0.993686</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense/difference</td>
<td>14.3698 ns [14.3382, 14.4041] (0.95 confidence)</td>
<td>14.5754 ns [14.48, 14.6828] (0.95 confidence)</td>
<td>0.985889</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense/difference&#95;endpoints</td>
<td>16.7801 ns [16.6628, 16.9076] (0.95 confidence)</td>
<td>16.666 ns [16.6437, 16.6899] (0.95 confidence)</td>
<td>1.00684</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense/dot</td>
<td>8.01941 ns [8.00493, 8.03349] (0.95 confidence)</td>
<td>8.10194 ns [8.06826, 8.1379] (0.95 confidence)</td>
<td>0.989813</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense/dot&#95;endpoints</td>
<td>8.62254 ns [8.56023, 8.69177] (0.95 confidence)</td>
<td>8.52613 ns [8.50419, 8.54995] (0.95 confidence)</td>
<td>1.01131</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense&#95;inexact/difference</td>
<td>14.5298 ns [14.442, 14.6313] (0.95 confidence)</td>
<td>14.6419 ns [14.5747, 14.7154] (0.95 confidence)</td>
<td>0.992342</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense&#95;inexact/difference&#95;endpoints</td>
<td>16.1563 ns [16.115, 16.2058] (0.95 confidence)</td>
<td>16.2958 ns [16.2126, 16.3888] (0.95 confidence)</td>
<td>0.991436</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense&#95;inexact/dot</td>
<td>8.08881 ns [8.05784, 8.12315] (0.95 confidence)</td>
<td>8.11891 ns [8.08975, 8.15135] (0.95 confidence)</td>
<td>0.996292</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense&#95;inexact/dot&#95;endpoints</td>
<td>9.62985 ns [9.56887, 9.69643] (0.95 confidence)</td>
<td>9.59879 ns [9.55957, 9.64159] (0.95 confidence)</td>
<td>1.00324</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/overflow/dot</td>
<td>0.998045 ns [0.996675, 0.999457] (0.95 confidence)</td>
<td>0.998778 ns [0.997121, 1.00056] (0.95 confidence)</td>
<td>0.999267</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/sparse/difference</td>
<td>5.46259 ns [5.45669, 5.46912] (0.95 confidence)</td>
<td>5.4777 ns [5.47038, 5.48691] (0.95 confidence)</td>
<td>0.997241</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/sparse/difference&#95;endpoints</td>
<td>6.36893 ns [6.35778, 6.38027] (0.95 confidence)</td>
<td>6.45672 ns [6.44136, 6.47175] (0.95 confidence)</td>
<td>0.986404</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/sparse/dot</td>
<td>2.74773 ns [2.73557, 2.76187] (0.95 confidence)</td>
<td>2.82947 ns [2.80645, 2.85302] (0.95 confidence)</td>
<td>0.971111</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/sparse/dot&#95;endpoints</td>
<td>2.80792 ns [2.78754, 2.82971] (0.95 confidence)</td>
<td>2.82123 ns [2.80115, 2.84231] (0.95 confidence)</td>
<td>0.995281</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/underflow/difference</td>
<td>4.7077 ns [4.70515, 4.71091] (0.95 confidence)</td>
<td>4.71013 ns [4.70627, 4.71495] (0.95 confidence)</td>
<td>0.999483</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/underflow/difference&#95;endpoints</td>
<td>4.9931 ns [4.98122, 5.00672] (0.95 confidence)</td>
<td>5.1341 ns [5.11337, 5.15415] (0.95 confidence)</td>
<td>0.972538</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/underflow/dot</td>
<td>3.16979 ns [3.1632, 3.17645] (0.95 confidence)</td>
<td>3.20029 ns [3.18431, 3.22686] (0.95 confidence)</td>
<td>0.990471</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/underflow/dot&#95;endpoints</td>
<td>3.24022 ns [3.23101, 3.24958] (0.95 confidence)</td>
<td>3.19701 ns [3.19165, 3.20275] (0.95 confidence)</td>
<td>1.01352</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4/dot&#95;bounded&#95;inconclusive</td>
<td>6.31165 ns [6.23781, 6.39288] (0.95 confidence)</td>
<td>6.23999 ns [6.2011, 6.28481] (0.95 confidence)</td>
<td>1.01148</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4/dot&#95;bounded&#95;well&#95;separated</td>
<td>9.56412 ns [9.51191, 9.62502] (0.95 confidence)</td>
<td>9.89911 ns [9.78455, 10.0194] (0.95 confidence)</td>
<td>0.966159</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4/dot&#95;difference&#95;bounded&#95;inconclusive</td>
<td>20.6885 ns [20.6562, 20.7194] (0.95 confidence)</td>
<td>21.0529 ns [20.9513, 21.167] (0.95 confidence)</td>
<td>0.982692</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4/dot&#95;difference&#95;bounded&#95;well&#95;separated</td>
<td>16.7313 ns [16.6366, 16.8404] (0.95 confidence)</td>
<td>17.6837 ns [17.2195, 18.2384] (0.95 confidence)</td>
<td>0.946142</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4/dot&#95;plain&#95;well&#95;separated</td>
<td>0.683334 ns [0.681647, 0.685138] (0.95 confidence)</td>
<td>0.683854 ns [0.68214, 0.685732] (0.95 confidence)</td>
<td>0.999239</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;batch/origin/construction</td>
<td>2.35147 ns [2.34546, 2.35748] (0.95 confidence)</td>
<td>2.36933 ns [2.36262, 2.37594] (0.95 confidence)</td>
<td>0.992466</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;batch/origin/prepared</td>
<td>34.8106 ns [34.6845, 34.9446] (0.95 confidence)</td>
<td>35.4746 ns [35.3044, 35.6504] (0.95 confidence)</td>
<td>0.981282</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;batch/origin/with&#95;construction</td>
<td>41.933 ns [41.7975, 42.0668] (0.95 confidence)</td>
<td>42.3385 ns [42.1681, 42.5148] (0.95 confidence)</td>
<td>0.990423</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;batch/translated/construction</td>
<td>2.35254 ns [2.34647, 2.35865] (0.95 confidence)</td>
<td>2.37367 ns [2.36671, 2.38063] (0.95 confidence)</td>
<td>0.9911</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;batch/translated/prepared</td>
<td>108.937 ns [108.779, 109.096] (0.95 confidence)</td>
<td>111.426 ns [110.694, 112.205] (0.95 confidence)</td>
<td>0.977668</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;batch/translated/with&#95;construction</td>
<td>115.702 ns [115.179, 116.282] (0.95 confidence)</td>
<td>116.682 ns [116.046, 117.363] (0.95 confidence)</td>
<td>0.991598</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/cancellation/difference</td>
<td>11.9422 ns [11.8363, 12.0567] (0.95 confidence)</td>
<td>11.8907 ns [11.7806, 12.0097] (0.95 confidence)</td>
<td>1.00433</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/cancellation/difference&#95;endpoints</td>
<td>12.9838 ns [12.8706, 13.107] (0.95 confidence)</td>
<td>12.8108 ns [12.7527, 12.8704] (0.95 confidence)</td>
<td>1.01351</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/cancellation/dot</td>
<td>6.39527 ns [6.3569, 6.44032] (0.95 confidence)</td>
<td>6.333 ns [6.29554, 6.37572] (0.95 confidence)</td>
<td>1.00983</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/cancellation/dot&#95;endpoints</td>
<td>6.91863 ns [6.91009, 6.92751] (0.95 confidence)</td>
<td>7.02145 ns [6.97679, 7.07157] (0.95 confidence)</td>
<td>0.985356</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense/difference</td>
<td>18.7945 ns [18.6902, 18.9098] (0.95 confidence)</td>
<td>18.8509 ns [18.7189, 18.9949] (0.95 confidence)</td>
<td>0.997011</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense/difference&#95;endpoints</td>
<td>21.1613 ns [21.0901, 21.2429] (0.95 confidence)</td>
<td>21.0593 ns [20.9621, 21.1693] (0.95 confidence)</td>
<td>1.00485</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense/dot</td>
<td>9.73801 ns [9.6637, 9.82242] (0.95 confidence)</td>
<td>9.73918 ns [9.66777, 9.82005] (0.95 confidence)</td>
<td>0.999879</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense/dot&#95;endpoints</td>
<td>11.1794 ns [11.0954, 11.2695] (0.95 confidence)</td>
<td>11.231 ns [11.2166, 11.2446] (0.95 confidence)</td>
<td>0.995403</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense&#95;inexact/difference</td>
<td>18.916 ns [18.8063, 19.0355] (0.95 confidence)</td>
<td>19.0916 ns [18.9533, 19.2384] (0.95 confidence)</td>
<td>0.990805</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense&#95;inexact/difference&#95;endpoints</td>
<td>20.3946 ns [20.2815, 20.5192] (0.95 confidence)</td>
<td>20.3819 ns [20.2802, 20.4948] (0.95 confidence)</td>
<td>1.00062</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense&#95;inexact/dot</td>
<td>9.80356 ns [9.73469, 9.88256] (0.95 confidence)</td>
<td>9.98102 ns [9.88749, 10.0784] (0.95 confidence)</td>
<td>0.982221</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense&#95;inexact/dot&#95;endpoints</td>
<td>10.3516 ns [10.3358, 10.3687] (0.95 confidence)</td>
<td>10.3596 ns [10.3376, 10.3821] (0.95 confidence)</td>
<td>0.999226</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/overflow/dot</td>
<td>1.00229 ns [1.00038, 1.00423] (0.95 confidence)</td>
<td>1.00009 ns [0.998305, 1.00224] (0.95 confidence)</td>
<td>1.0022</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/sparse/difference</td>
<td>6.2159 ns [6.20891, 6.22411] (0.95 confidence)</td>
<td>6.20861 ns [6.19999, 6.21941] (0.95 confidence)</td>
<td>1.00117</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/sparse/difference&#95;endpoints</td>
<td>7.49102 ns [7.4791, 7.5024] (0.95 confidence)</td>
<td>7.44052 ns [7.42699, 7.45449] (0.95 confidence)</td>
<td>1.00679</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/sparse/dot</td>
<td>3.18973 ns [3.17727, 3.20223] (0.95 confidence)</td>
<td>3.1876 ns [3.17008, 3.20699] (0.95 confidence)</td>
<td>1.00067</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/sparse/dot&#95;endpoints</td>
<td>3.24935 ns [3.23787, 3.26231] (0.95 confidence)</td>
<td>3.2239 ns [3.21534, 3.23348] (0.95 confidence)</td>
<td>1.00789</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/underflow/difference</td>
<td>5.70562 ns [5.69771, 5.71627] (0.95 confidence)</td>
<td>5.7018 ns [5.69259, 5.71707] (0.95 confidence)</td>
<td>1.00067</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/underflow/difference&#95;endpoints</td>
<td>6.04947 ns [6.01027, 6.09313] (0.95 confidence)</td>
<td>6.04133 ns [6.01649, 6.06825] (0.95 confidence)</td>
<td>1.00135</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/underflow/dot</td>
<td>3.89893 ns [3.89493, 3.90269] (0.95 confidence)</td>
<td>3.90204 ns [3.89635, 3.90768] (0.95 confidence)</td>
<td>0.999203</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/underflow/dot&#95;endpoints</td>
<td>3.91469 ns [3.90706, 3.92278] (0.95 confidence)</td>
<td>3.91815 ns [3.90667, 3.93317] (0.95 confidence)</td>
<td>0.999117</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;batch/origin/construction</td>
<td>2.90408 ns [2.89819, 2.90973] (0.95 confidence)</td>
<td>2.92388 ns [2.91953, 2.9284] (0.95 confidence)</td>
<td>0.993225</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;batch/origin/prepared</td>
<td>46.9513 ns [46.6971, 47.2286] (0.95 confidence)</td>
<td>47.0958 ns [46.8841, 47.3399] (0.95 confidence)</td>
<td>0.996932</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;batch/origin/with&#95;construction</td>
<td>56.6745 ns [56.5915, 56.7643] (0.95 confidence)</td>
<td>58.3173 ns [58.1535, 58.4775] (0.95 confidence)</td>
<td>0.971831</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;batch/translated/construction</td>
<td>3.2622 ns [3.2567, 3.26754] (0.95 confidence)</td>
<td>3.26554 ns [3.25979, 3.27114] (0.95 confidence)</td>
<td>0.998978</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;batch/translated/prepared</td>
<td>145.53 ns [145.324, 145.738] (0.95 confidence)</td>
<td>147.159 ns [146.485, 147.933] (0.95 confidence)</td>
<td>0.988929</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;batch/translated/with&#95;construction</td>
<td>158.477 ns [157.451, 159.612] (0.95 confidence)</td>
<td>157.86 ns [157.359, 158.417] (0.95 confidence)</td>
<td>1.0039</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/cancellation/difference</td>
<td>12.366 ns [12.262, 12.48] (0.95 confidence)</td>
<td>12.2501 ns [12.1634, 12.3495] (0.95 confidence)</td>
<td>1.00946</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/cancellation/difference&#95;endpoints</td>
<td>13.2627 ns [13.2194, 13.3072] (0.95 confidence)</td>
<td>13.4686 ns [13.3564, 13.591] (0.95 confidence)</td>
<td>0.984711</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/cancellation/dot</td>
<td>6.74787 ns [6.69223, 6.80804] (0.95 confidence)</td>
<td>6.71234 ns [6.66675, 6.7653] (0.95 confidence)</td>
<td>1.00529</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/cancellation/dot&#95;endpoints</td>
<td>7.29322 ns [7.24745, 7.34526] (0.95 confidence)</td>
<td>7.4328 ns [7.37818, 7.49074] (0.95 confidence)</td>
<td>0.981222</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense/difference</td>
<td>23.0432 ns [22.9148, 23.1879] (0.95 confidence)</td>
<td>23.0735 ns [22.9512, 23.2102] (0.95 confidence)</td>
<td>0.998687</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense/difference&#95;endpoints</td>
<td>24.222 ns [24.1857, 24.2571] (0.95 confidence)</td>
<td>24.2532 ns [24.2141, 24.2928] (0.95 confidence)</td>
<td>0.998713</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense/dot</td>
<td>12.0636 ns [11.9281, 12.206] (0.95 confidence)</td>
<td>11.9981 ns [11.8729, 12.1294] (0.95 confidence)</td>
<td>1.00546</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense/dot&#95;endpoints</td>
<td>12.4035 ns [12.3474, 12.4667] (0.95 confidence)</td>
<td>12.4372 ns [12.386, 12.4955] (0.95 confidence)</td>
<td>0.997289</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense&#95;inexact/difference</td>
<td>22.9642 ns [22.8655, 23.0818] (0.95 confidence)</td>
<td>23.3279 ns [23.1905, 23.4755] (0.95 confidence)</td>
<td>0.984407</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense&#95;inexact/difference&#95;endpoints</td>
<td>25.4994 ns [25.3785, 25.6342] (0.95 confidence)</td>
<td>25.7273 ns [25.5838, 25.8832] (0.95 confidence)</td>
<td>0.991139</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense&#95;inexact/dot</td>
<td>11.7826 ns [11.7136, 11.8609] (0.95 confidence)</td>
<td>11.9201 ns [11.8265, 12.0224] (0.95 confidence)</td>
<td>0.988468</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense&#95;inexact/dot&#95;endpoints</td>
<td>13.478 ns [13.3777, 13.5868] (0.95 confidence)</td>
<td>13.6352 ns [13.5267, 13.7494] (0.95 confidence)</td>
<td>0.988472</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/overflow/dot</td>
<td>0.996116 ns [0.994764, 0.997553] (0.95 confidence)</td>
<td>1.00049 ns [0.998712, 1.00235] (0.95 confidence)</td>
<td>0.995624</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/sparse/difference</td>
<td>7.21477 ns [7.2068, 7.22352] (0.95 confidence)</td>
<td>7.21624 ns [7.20475, 7.22818] (0.95 confidence)</td>
<td>0.999796</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/sparse/difference&#95;endpoints</td>
<td>8.3687 ns [8.35875, 8.37892] (0.95 confidence)</td>
<td>8.35154 ns [8.34086, 8.36247] (0.95 confidence)</td>
<td>1.00205</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/sparse/dot</td>
<td>3.68962 ns [3.67043, 3.71038] (0.95 confidence)</td>
<td>3.66276 ns [3.65066, 3.67589] (0.95 confidence)</td>
<td>1.00733</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/sparse/dot&#95;endpoints</td>
<td>3.77165 ns [3.76074, 3.78362] (0.95 confidence)</td>
<td>3.77063 ns [3.76171, 3.78017] (0.95 confidence)</td>
<td>1.00027</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/underflow/difference</td>
<td>7.69628 ns [7.68241, 7.71061] (0.95 confidence)</td>
<td>7.74628 ns [7.7199, 7.77698] (0.95 confidence)</td>
<td>0.993545</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/underflow/difference&#95;endpoints</td>
<td>8.17004 ns [8.15511, 8.18553] (0.95 confidence)</td>
<td>8.36952 ns [8.33986, 8.39916] (0.95 confidence)</td>
<td>0.976166</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/underflow/dot</td>
<td>4.62296 ns [4.61672, 4.6292] (0.95 confidence)</td>
<td>4.63848 ns [4.63041, 4.64609] (0.95 confidence)</td>
<td>0.996654</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/underflow/dot&#95;endpoints</td>
<td>4.61712 ns [4.60891, 4.62494] (0.95 confidence)</td>
<td>4.65893 ns [4.64984, 4.6691] (0.95 confidence)</td>
<td>0.991025</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;batch/origin/construction</td>
<td>3.10441 ns [3.09908, 3.10959] (0.95 confidence)</td>
<td>3.1123 ns [3.10548, 3.12215] (0.95 confidence)</td>
<td>0.997466</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;batch/origin/prepared</td>
<td>61.2179 ns [61.0688, 61.3787] (0.95 confidence)</td>
<td>61.4993 ns [61.2847, 61.726] (0.95 confidence)</td>
<td>0.995424</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;batch/origin/with&#95;construction</td>
<td>82.2115 ns [81.8953, 82.5436] (0.95 confidence)</td>
<td>81.3002 ns [81.1264, 81.4867] (0.95 confidence)</td>
<td>1.01121</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;batch/translated/construction</td>
<td>3.59019 ns [3.58707, 3.59338] (0.95 confidence)</td>
<td>3.58847 ns [3.58244, 3.59445] (0.95 confidence)</td>
<td>1.00048</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;batch/translated/prepared</td>
<td>205.236 ns [204.465, 206.098] (0.95 confidence)</td>
<td>207.339 ns [206.339, 208.372] (0.95 confidence)</td>
<td>0.989856</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;batch/translated/with&#95;construction</td>
<td>231.992 ns [230.811, 233.246] (0.95 confidence)</td>
<td>229.236 ns [228.212, 230.375] (0.95 confidence)</td>
<td>1.01202</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/cancellation/difference</td>
<td>12.7864 ns [12.7041, 12.8786] (0.95 confidence)</td>
<td>12.6224 ns [12.5675, 12.6846] (0.95 confidence)</td>
<td>1.01299</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/cancellation/difference&#95;endpoints</td>
<td>14.0663 ns [14.019, 14.1258] (0.95 confidence)</td>
<td>14.1341 ns [14.0767, 14.1983] (0.95 confidence)</td>
<td>0.995207</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/cancellation/dot</td>
<td>7.02234 ns [6.96024, 7.09031] (0.95 confidence)</td>
<td>7.19805 ns [7.12814, 7.2722] (0.95 confidence)</td>
<td>0.975589</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/cancellation/dot&#95;endpoints</td>
<td>7.6217 ns [7.59732, 7.65381] (0.95 confidence)</td>
<td>7.75012 ns [7.71199, 7.79083] (0.95 confidence)</td>
<td>0.98343</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense/difference</td>
<td>29.0351 ns [28.9877, 29.0826] (0.95 confidence)</td>
<td>29.5819 ns [29.2572, 29.9791] (0.95 confidence)</td>
<td>0.981516</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense/difference&#95;endpoints</td>
<td>32.3852 ns [32.327, 32.4421] (0.95 confidence)</td>
<td>32.616 ns [32.5405, 32.6986] (0.95 confidence)</td>
<td>0.992922</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense/dot</td>
<td>13.8853 ns [13.8281, 13.952] (0.95 confidence)</td>
<td>13.98 ns [13.8924, 14.0808] (0.95 confidence)</td>
<td>0.993227</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense/dot&#95;endpoints</td>
<td>15.3473 ns [15.2737, 15.4315] (0.95 confidence)</td>
<td>15.2087 ns [15.186, 15.2313] (0.95 confidence)</td>
<td>1.00911</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense&#95;inexact/difference</td>
<td>29.0994 ns [29.0563, 29.1423] (0.95 confidence)</td>
<td>29.8979 ns [29.4113, 30.5138] (0.95 confidence)</td>
<td>0.973291</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense&#95;inexact/difference&#95;endpoints</td>
<td>31.6566 ns [31.5833, 31.7274] (0.95 confidence)</td>
<td>31.879 ns [31.8221, 31.9348] (0.95 confidence)</td>
<td>0.993023</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense&#95;inexact/dot</td>
<td>13.898 ns [13.8391, 13.966] (0.95 confidence)</td>
<td>14.0992 ns [13.9958, 14.2104] (0.95 confidence)</td>
<td>0.985728</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense&#95;inexact/dot&#95;endpoints</td>
<td>15.4588 ns [15.3543, 15.5762] (0.95 confidence)</td>
<td>15.3926 ns [15.3199, 15.4746] (0.95 confidence)</td>
<td>1.0043</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/overflow/dot</td>
<td>1.00038 ns [0.998023, 1.0032] (0.95 confidence)</td>
<td>1.00098 ns [0.999833, 1.00217] (0.95 confidence)</td>
<td>0.999396</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/sparse/difference</td>
<td>8.2384 ns [8.20555, 8.27328] (0.95 confidence)</td>
<td>8.18599 ns [8.17165, 8.20098] (0.95 confidence)</td>
<td>1.0064</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/sparse/difference&#95;endpoints</td>
<td>9.38851 ns [9.36654, 9.41169] (0.95 confidence)</td>
<td>9.44922 ns [9.42666, 9.47115] (0.95 confidence)</td>
<td>0.993575</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/sparse/dot</td>
<td>4.11362 ns [4.10792, 4.12063] (0.95 confidence)</td>
<td>4.11907 ns [4.11298, 4.12596] (0.95 confidence)</td>
<td>0.998675</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/sparse/dot&#95;endpoints</td>
<td>4.17944 ns [4.17518, 4.1841] (0.95 confidence)</td>
<td>4.18663 ns [4.18102, 4.19125] (0.95 confidence)</td>
<td>0.998281</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/underflow/difference</td>
<td>8.63701 ns [8.62019, 8.65584] (0.95 confidence)</td>
<td>8.69364 ns [8.66655, 8.72893] (0.95 confidence)</td>
<td>0.993487</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/underflow/difference&#95;endpoints</td>
<td>8.99286 ns [8.96101, 9.02453] (0.95 confidence)</td>
<td>8.94974 ns [8.91306, 8.98708] (0.95 confidence)</td>
<td>1.00482</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/underflow/dot</td>
<td>5.32761 ns [5.31866, 5.33636] (0.95 confidence)</td>
<td>5.37242 ns [5.36089, 5.38439] (0.95 confidence)</td>
<td>0.99166</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/underflow/dot&#95;endpoints</td>
<td>5.39601 ns [5.38952, 5.40251] (0.95 confidence)</td>
<td>5.42105 ns [5.41318, 5.42875] (0.95 confidence)</td>
<td>0.99538</td>
</tr>
<tr>
<td>public&#95;d2/det&#95;bound</td>
<td>1.95663 ns [1.95359, 1.95971] (0.95 confidence)</td>
<td>1.96698 ns [1.9648, 1.96944] (0.95 confidence)</td>
<td>0.994737</td>
</tr>
<tr>
<td>public&#95;d2/det&#95;direct</td>
<td>1.13571 ns [1.13345, 1.13815] (0.95 confidence)</td>
<td>1.13417 ns [1.13241, 1.13595] (0.95 confidence)</td>
<td>1.00136</td>
</tr>
<tr>
<td>public&#95;d2/det&#95;with&#95;bound</td>
<td>2.12177 ns [2.11698, 2.12684] (0.95 confidence)</td>
<td>2.11954 ns [2.11623, 2.12279] (0.95 confidence)</td>
<td>1.00105</td>
</tr>
<tr>
<td>public&#95;d2/first&#95;asymmetry</td>
<td>1.38137 ns [1.37969, 1.38296] (0.95 confidence)</td>
<td>1.38871 ns [1.38647, 1.39116] (0.95 confidence)</td>
<td>0.994714</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;det</td>
<td>7.83106 ns [7.78978, 7.87308] (0.95 confidence)</td>
<td>7.9092 ns [7.85764, 7.96317] (0.95 confidence)</td>
<td>0.990121</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;access</td>
<td>6.36266 ns [6.35093, 6.37437] (0.95 confidence)</td>
<td>6.36852 ns [6.35477, 6.38176] (0.95 confidence)</td>
<td>0.99908</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;from&#95;matrix</td>
<td>2.60218 ns [2.59741, 2.60699] (0.95 confidence)</td>
<td>2.60341 ns [2.59814, 2.60867] (0.95 confidence)</td>
<td>0.999529</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;from&#95;points</td>
<td>2.25862 ns [2.24985, 2.26771] (0.95 confidence)</td>
<td>2.25701 ns [2.24838, 2.26586] (0.95 confidence)</td>
<td>1.00071</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;from&#95;rows</td>
<td>3.57968 ns [3.57317, 3.58615] (0.95 confidence)</td>
<td>3.57148 ns [3.56461, 3.57798] (0.95 confidence)</td>
<td>1.0023</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;identity</td>
<td>2.10118 ns [2.09626, 2.10626] (0.95 confidence)</td>
<td>2.10849 ns [2.10341, 2.11393] (0.95 confidence)</td>
<td>0.99653</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;set</td>
<td>3.82329 ns [3.81572, 3.83054] (0.95 confidence)</td>
<td>3.81506 ns [3.80636, 3.82361] (0.95 confidence)</td>
<td>1.00216</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;zero</td>
<td>2.1129 ns [2.1085, 2.11725] (0.95 confidence)</td>
<td>2.12446 ns [2.12044, 2.12865] (0.95 confidence)</td>
<td>0.994558</td>
</tr>
<tr>
<td>public&#95;d2/is&#95;symmetric</td>
<td>0.759106 ns [0.756601, 0.762752] (0.95 confidence)</td>
<td>0.760366 ns [0.759421, 0.761445] (0.95 confidence)</td>
<td>0.998343</td>
</tr>
<tr>
<td>public&#95;d2/matrix&#95;access</td>
<td>4.76113 ns [4.7487, 4.77351] (0.95 confidence)</td>
<td>4.76655 ns [4.75658, 4.77612] (0.95 confidence)</td>
<td>0.998862</td>
</tr>
<tr>
<td>public&#95;d2/matrix&#95;from&#95;rows</td>
<td>2.91122 ns [2.9069, 2.91542] (0.95 confidence)</td>
<td>2.92121 ns [2.9177, 2.92461] (0.95 confidence)</td>
<td>0.996582</td>
</tr>
<tr>
<td>public&#95;d2/matrix&#95;identity</td>
<td>1.25794 ns [1.25491, 1.26089] (0.95 confidence)</td>
<td>1.2591 ns [1.25589, 1.26222] (0.95 confidence)</td>
<td>0.99908</td>
</tr>
<tr>
<td>public&#95;d2/matrix&#95;set</td>
<td>2.11758 ns [2.10853, 2.12657] (0.95 confidence)</td>
<td>2.10445 ns [2.09628, 2.11259] (0.95 confidence)</td>
<td>1.00624</td>
</tr>
<tr>
<td>public&#95;d2/matrix&#95;zero</td>
<td>1.16319 ns [1.16081, 1.16558] (0.95 confidence)</td>
<td>1.16598 ns [1.16343, 1.16849] (0.95 confidence)</td>
<td>0.997602</td>
</tr>
<tr>
<td>public&#95;d2/vector&#95;access</td>
<td>1.6226 ns [1.61801, 1.62776] (0.95 confidence)</td>
<td>1.631 ns [1.62743, 1.63472] (0.95 confidence)</td>
<td>0.994847</td>
</tr>
<tr>
<td>public&#95;d2/vector&#95;new</td>
<td>1.91323 ns [1.89864, 1.92658] (0.95 confidence)</td>
<td>1.89729 ns [1.87861, 1.91494] (0.95 confidence)</td>
<td>1.0084</td>
</tr>
<tr>
<td>public&#95;d2/vector&#95;zero</td>
<td>0.591217 ns [0.589914, 0.592467] (0.95 confidence)</td>
<td>0.59487 ns [0.593369, 0.596387] (0.95 confidence)</td>
<td>0.993859</td>
</tr>
<tr>
<td>public&#95;d3/det&#95;bound</td>
<td>3.66705 ns [3.66341, 3.67141] (0.95 confidence)</td>
<td>3.65708 ns [3.65291, 3.66138] (0.95 confidence)</td>
<td>1.00272</td>
</tr>
<tr>
<td>public&#95;d3/det&#95;direct</td>
<td>1.29753 ns [1.28939, 1.30941] (0.95 confidence)</td>
<td>1.29289 ns [1.29087, 1.29527] (0.95 confidence)</td>
<td>1.00359</td>
</tr>
<tr>
<td>public&#95;d3/det&#95;with&#95;bound</td>
<td>3.67438 ns [3.67053, 3.67825] (0.95 confidence)</td>
<td>3.66812 ns [3.66359, 3.673] (0.95 confidence)</td>
<td>1.00171</td>
</tr>
<tr>
<td>public&#95;d3/first&#95;asymmetry</td>
<td>3.59173 ns [3.55874, 3.62663] (0.95 confidence)</td>
<td>3.60537 ns [3.57417, 3.63806] (0.95 confidence)</td>
<td>0.996216</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;det</td>
<td>23.2804 ns [23.0141, 23.5449] (0.95 confidence)</td>
<td>23.7368 ns [23.4962, 23.9699] (0.95 confidence)</td>
<td>0.980772</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;access</td>
<td>10.0095 ns [9.98622, 10.0391] (0.95 confidence)</td>
<td>10.1642 ns [10.1385, 10.1923] (0.95 confidence)</td>
<td>0.984783</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;from&#95;matrix</td>
<td>5.17688 ns [5.16629, 5.18761] (0.95 confidence)</td>
<td>5.20121 ns [5.18725, 5.21569] (0.95 confidence)</td>
<td>0.995323</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;from&#95;points</td>
<td>6.3279 ns [6.3176, 6.33762] (0.95 confidence)</td>
<td>6.41434 ns [6.39991, 6.42976] (0.95 confidence)</td>
<td>0.986524</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;from&#95;rows</td>
<td>7.37206 ns [7.35535, 7.38813] (0.95 confidence)</td>
<td>7.42261 ns [7.40035, 7.44448] (0.95 confidence)</td>
<td>0.993189</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;identity</td>
<td>4.33255 ns [4.32363, 4.34169] (0.95 confidence)</td>
<td>4.55988 ns [4.52194, 4.59862] (0.95 confidence)</td>
<td>0.950145</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;set</td>
<td>7.8115 ns [7.79265, 7.82983] (0.95 confidence)</td>
<td>7.98022 ns [7.92787, 8.03982] (0.95 confidence)</td>
<td>0.978857</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;zero</td>
<td>4.32797 ns [4.32073, 4.33563] (0.95 confidence)</td>
<td>4.36584 ns [4.35608, 4.37566] (0.95 confidence)</td>
<td>0.991326</td>
</tr>
<tr>
<td>public&#95;d3/is&#95;symmetric</td>
<td>1.49686 ns [1.49388, 1.50015] (0.95 confidence)</td>
<td>1.34464 ns [1.34237, 1.34684] (0.95 confidence)</td>
<td>1.1132</td>
</tr>
<tr>
<td>public&#95;d3/matrix&#95;access</td>
<td>6.68839 ns [6.66674, 6.70996] (0.95 confidence)</td>
<td>6.96002 ns [6.92456, 6.99553] (0.95 confidence)</td>
<td>0.960973</td>
</tr>
<tr>
<td>public&#95;d3/matrix&#95;from&#95;rows</td>
<td>5.71821 ns [5.70504, 5.73195] (0.95 confidence)</td>
<td>5.74111 ns [5.7295, 5.75305] (0.95 confidence)</td>
<td>0.996011</td>
</tr>
<tr>
<td>public&#95;d3/matrix&#95;identity</td>
<td>2.71203 ns [2.70598, 2.71838] (0.95 confidence)</td>
<td>2.72135 ns [2.71617, 2.72646] (0.95 confidence)</td>
<td>0.996575</td>
</tr>
<tr>
<td>public&#95;d3/matrix&#95;set</td>
<td>4.24606 ns [4.23142, 4.26049] (0.95 confidence)</td>
<td>4.34244 ns [4.32169, 4.36321] (0.95 confidence)</td>
<td>0.977805</td>
</tr>
<tr>
<td>public&#95;d3/matrix&#95;zero</td>
<td>2.41589 ns [2.41069, 2.42111] (0.95 confidence)</td>
<td>2.43452 ns [2.42969, 2.43941] (0.95 confidence)</td>
<td>0.992347</td>
</tr>
<tr>
<td>public&#95;d3/vector&#95;access</td>
<td>2.0971 ns [2.09404, 2.1001] (0.95 confidence)</td>
<td>2.10487 ns [2.10132, 2.10826] (0.95 confidence)</td>
<td>0.996306</td>
</tr>
<tr>
<td>public&#95;d3/vector&#95;new</td>
<td>3.65984 ns [3.64735, 3.67339] (0.95 confidence)</td>
<td>3.69119 ns [3.66372, 3.72186] (0.95 confidence)</td>
<td>0.991508</td>
</tr>
<tr>
<td>public&#95;d3/vector&#95;zero</td>
<td>0.958735 ns [0.957142, 0.960298] (0.95 confidence)</td>
<td>0.967774 ns [0.965054, 0.970472] (0.95 confidence)</td>
<td>0.990659</td>
</tr>
<tr>
<td>public&#95;d4/det&#95;bound</td>
<td>6.21432 ns [6.19992, 6.23037] (0.95 confidence)</td>
<td>6.85308 ns [6.7866, 6.91548] (0.95 confidence)</td>
<td>0.906793</td>
</tr>
<tr>
<td>public&#95;d4/det&#95;direct</td>
<td>2.47792 ns [2.47577, 2.48042] (0.95 confidence)</td>
<td>2.76436 ns [2.74583, 2.78155] (0.95 confidence)</td>
<td>0.896383</td>
</tr>
<tr>
<td>public&#95;d4/det&#95;with&#95;bound</td>
<td>6.2079 ns [6.1947, 6.22246] (0.95 confidence)</td>
<td>6.27901 ns [6.25526, 6.30438] (0.95 confidence)</td>
<td>0.988674</td>
</tr>
<tr>
<td>public&#95;d4/first&#95;asymmetry</td>
<td>3.67001 ns [3.63795, 3.70357] (0.95 confidence)</td>
<td>3.68008 ns [3.64809, 3.71334] (0.95 confidence)</td>
<td>0.997265</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;det</td>
<td>49.0357 ns [48.714, 49.4048] (0.95 confidence)</td>
<td>50.6666 ns [50.121, 51.2536] (0.95 confidence)</td>
<td>0.967809</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;access</td>
<td>13.642 ns [13.6017, 13.6817] (0.95 confidence)</td>
<td>13.7785 ns [13.74, 13.8169] (0.95 confidence)</td>
<td>0.990089</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;from&#95;matrix</td>
<td>7.92668 ns [7.91391, 7.939] (0.95 confidence)</td>
<td>8.08358 ns [8.05685, 8.11349] (0.95 confidence)</td>
<td>0.980591</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;from&#95;points</td>
<td>10.8812 ns [10.8612, 10.8994] (0.95 confidence)</td>
<td>11.0209 ns [11.0062, 11.0353] (0.95 confidence)</td>
<td>0.987324</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;from&#95;rows</td>
<td>11.2044 ns [11.1847, 11.2239] (0.95 confidence)</td>
<td>11.4335 ns [11.4089, 11.4574] (0.95 confidence)</td>
<td>0.979968</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;identity</td>
<td>6.5333 ns [6.5189, 6.54741] (0.95 confidence)</td>
<td>6.72461 ns [6.68166, 6.77739] (0.95 confidence)</td>
<td>0.971551</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;set</td>
<td>11.4603 ns [11.4366, 11.4832] (0.95 confidence)</td>
<td>11.6313 ns [11.6097, 11.6535] (0.95 confidence)</td>
<td>0.985301</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;zero</td>
<td>6.62551 ns [6.61557, 6.63596] (0.95 confidence)</td>
<td>6.87306 ns [6.82653, 6.9261] (0.95 confidence)</td>
<td>0.963984</td>
</tr>
<tr>
<td>public&#95;d4/is&#95;symmetric</td>
<td>2.10936 ns [2.10755, 2.11142] (0.95 confidence)</td>
<td>2.11679 ns [2.11351, 2.12034] (0.95 confidence)</td>
<td>0.996491</td>
</tr>
<tr>
<td>public&#95;d4/matrix&#95;access</td>
<td>9.42128 ns [9.4025, 9.43848] (0.95 confidence)</td>
<td>9.75422 ns [9.72526, 9.78597] (0.95 confidence)</td>
<td>0.965868</td>
</tr>
<tr>
<td>public&#95;d4/matrix&#95;from&#95;rows</td>
<td>10.3244 ns [10.3137, 10.3358] (0.95 confidence)</td>
<td>10.5594 ns [10.5267, 10.5943] (0.95 confidence)</td>
<td>0.977747</td>
</tr>
<tr>
<td>public&#95;d4/matrix&#95;identity</td>
<td>4.13039 ns [4.12001, 4.14119] (0.95 confidence)</td>
<td>4.40613 ns [4.37055, 4.44223] (0.95 confidence)</td>
<td>0.93742</td>
</tr>
<tr>
<td>public&#95;d4/matrix&#95;set</td>
<td>7.36402 ns [7.34744, 7.38032] (0.95 confidence)</td>
<td>7.63308 ns [7.61104, 7.65503] (0.95 confidence)</td>
<td>0.964751</td>
</tr>
<tr>
<td>public&#95;d4/matrix&#95;zero</td>
<td>3.8696 ns [3.86343, 3.876] (0.95 confidence)</td>
<td>3.97335 ns [3.95883, 3.98822] (0.95 confidence)</td>
<td>0.973889</td>
</tr>
<tr>
<td>public&#95;d4/vector&#95;access</td>
<td>2.62395 ns [2.6184, 2.6291] (0.95 confidence)</td>
<td>2.68185 ns [2.67424, 2.68991] (0.95 confidence)</td>
<td>0.978412</td>
</tr>
<tr>
<td>public&#95;d4/vector&#95;new</td>
<td>2.82418 ns [2.7995, 2.84699] (0.95 confidence)</td>
<td>2.85455 ns [2.82425, 2.88264] (0.95 confidence)</td>
<td>0.989361</td>
</tr>
<tr>
<td>public&#95;d4/vector&#95;zero</td>
<td>1.17665 ns [1.17455, 1.17882] (0.95 confidence)</td>
<td>1.20928 ns [1.20635, 1.2123] (0.95 confidence)</td>
<td>0.973019</td>
</tr>
<tr>
<td>public&#95;d5/det&#95;bound</td>
<td>1.34042 ns [1.33789, 1.34313] (0.95 confidence)</td>
<td>1.3597 ns [1.35709, 1.36237] (0.95 confidence)</td>
<td>0.98582</td>
</tr>
<tr>
<td>public&#95;d5/det&#95;direct</td>
<td>1.33682 ns [1.33445, 1.33929] (0.95 confidence)</td>
<td>1.37835 ns [1.37163, 1.38576] (0.95 confidence)</td>
<td>0.969874</td>
</tr>
<tr>
<td>public&#95;d5/det&#95;with&#95;bound</td>
<td>0.856137 ns [0.854077, 0.858175] (0.95 confidence)</td>
<td>0.867654 ns [0.865948, 0.869314] (0.95 confidence)</td>
<td>0.986727</td>
</tr>
<tr>
<td>public&#95;d5/first&#95;asymmetry</td>
<td>3.79095 ns [3.78764, 3.79474] (0.95 confidence)</td>
<td>3.7929 ns [3.78855, 3.79797] (0.95 confidence)</td>
<td>0.999486</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;det</td>
<td>125.697 ns [125.093, 126.291] (0.95 confidence)</td>
<td>135.893 ns [134.145, 137.675] (0.95 confidence)</td>
<td>0.924968</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;access</td>
<td>16.3934 ns [16.176, 16.6128] (0.95 confidence)</td>
<td>16.1987 ns [16.0289, 16.3804] (0.95 confidence)</td>
<td>1.01202</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;from&#95;matrix</td>
<td>10.6754 ns [10.6442, 10.7062] (0.95 confidence)</td>
<td>10.7401 ns [10.7092, 10.7712] (0.95 confidence)</td>
<td>0.993974</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;from&#95;points</td>
<td>18.4856 ns [18.459, 18.5121] (0.95 confidence)</td>
<td>18.6114 ns [18.5789, 18.6461] (0.95 confidence)</td>
<td>0.993238</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;from&#95;rows</td>
<td>12.7277 ns [12.6886, 12.771] (0.95 confidence)</td>
<td>12.7395 ns [12.7094, 12.7734] (0.95 confidence)</td>
<td>0.999077</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;identity</td>
<td>8.36595 ns [8.30095, 8.4314] (0.95 confidence)</td>
<td>8.52994 ns [8.46352, 8.59637] (0.95 confidence)</td>
<td>0.980776</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;set</td>
<td>13.0304 ns [12.9688, 13.0936] (0.95 confidence)</td>
<td>13.2614 ns [13.1693, 13.3579] (0.95 confidence)</td>
<td>0.982588</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;zero</td>
<td>9.35655 ns [9.34208, 9.37174] (0.95 confidence)</td>
<td>9.66481 ns [9.62602, 9.70629] (0.95 confidence)</td>
<td>0.968104</td>
</tr>
<tr>
<td>public&#95;d5/is&#95;symmetric</td>
<td>3.69614 ns [3.69355, 3.69917] (0.95 confidence)</td>
<td>3.69895 ns [3.69545, 3.70294] (0.95 confidence)</td>
<td>0.999241</td>
</tr>
<tr>
<td>public&#95;d5/matrix&#95;access</td>
<td>12.1725 ns [12.1255, 12.2193] (0.95 confidence)</td>
<td>12.2775 ns [12.2341, 12.3209] (0.95 confidence)</td>
<td>0.991452</td>
</tr>
<tr>
<td>public&#95;d5/matrix&#95;from&#95;rows</td>
<td>11.79 ns [11.7742, 11.8061] (0.95 confidence)</td>
<td>11.9163 ns [11.8912, 11.9465] (0.95 confidence)</td>
<td>0.989401</td>
</tr>
<tr>
<td>public&#95;d5/matrix&#95;identity</td>
<td>6.80095 ns [6.78319, 6.81918] (0.95 confidence)</td>
<td>6.93749 ns [6.91795, 6.95708] (0.95 confidence)</td>
<td>0.980319</td>
</tr>
<tr>
<td>public&#95;d5/matrix&#95;set</td>
<td>10.4823 ns [10.4597, 10.5051] (0.95 confidence)</td>
<td>10.7711 ns [10.606, 11.0485] (0.95 confidence)</td>
<td>0.973182</td>
</tr>
<tr>
<td>public&#95;d5/matrix&#95;zero</td>
<td>5.70823 ns [5.69737, 5.71922] (0.95 confidence)</td>
<td>5.85968 ns [5.8467, 5.87287] (0.95 confidence)</td>
<td>0.974153</td>
</tr>
<tr>
<td>public&#95;d5/vector&#95;access</td>
<td>3.09438 ns [3.08873, 3.09981] (0.95 confidence)</td>
<td>3.09694 ns [3.09089, 3.10261] (0.95 confidence)</td>
<td>0.999173</td>
</tr>
<tr>
<td>public&#95;d5/vector&#95;new</td>
<td>4.85431 ns [4.84603, 4.86247] (0.95 confidence)</td>
<td>4.87153 ns [4.85761, 4.88744] (0.95 confidence)</td>
<td>0.996465</td>
</tr>
<tr>
<td>public&#95;d5/vector&#95;zero</td>
<td>1.5021 ns [1.49286, 1.51118] (0.95 confidence)</td>
<td>1.53946 ns [1.52998, 1.54912] (0.95 confidence)</td>
<td>0.975732</td>
</tr>
<tr>
<td>public&#95;scalar/interval&#95;access</td>
<td>1.59123 ns [1.58641, 1.59615] (0.95 confidence)</td>
<td>1.5922 ns [1.58748, 1.59706] (0.95 confidence)</td>
<td>0.99939</td>
</tr>
<tr>
<td>public&#95;scalar/interval&#95;negate</td>
<td>0.734152 ns [0.733212, 0.735135] (0.95 confidence)</td>
<td>0.734679 ns [0.733536, 0.735897] (0.95 confidence)</td>
<td>0.999283</td>
</tr>
<tr>
<td>public&#95;scalar/interval&#95;new</td>
<td>0.893035 ns [0.892049, 0.894031] (0.95 confidence)</td>
<td>0.893437 ns [0.892804, 0.894077] (0.95 confidence)</td>
<td>0.999549</td>
</tr>
<tr>
<td>public&#95;scalar/interval&#95;point</td>
<td>0.58228 ns [0.581627, 0.582952] (0.95 confidence)</td>
<td>0.590174 ns [0.586155, 0.59477] (0.95 confidence)</td>
<td>0.986624</td>
</tr>
<tr>
<td>public&#95;scalar/interval&#95;subtraction</td>
<td>1.86871 ns [1.86367, 1.8761] (0.95 confidence)</td>
<td>1.87726 ns [1.87262, 1.88218] (0.95 confidence)</td>
<td>0.995444</td>
</tr>
<tr>
<td>public&#95;scalar/tolerance</td>
<td>0.975598 ns [0.971658, 0.979046] (0.95 confidence)</td>
<td>0.97846 ns [0.97536, 0.981202] (0.95 confidence)</td>
<td>0.997075</td>
</tr>
</tbody>
</table>

### Median (ns)

<table>
<thead>
<tr>
<th>Workload</th>
<th>Baseline</th>
<th>Current</th>
<th>Baseline / current</th>
</tr>
</thead>
<tbody>
<tr>
<td>angle&#95;d3/dense/borrowed&#95;slice</td>
<td>23.7216 ns [23.6954, 23.777] (0.95 confidence)</td>
<td>24.2331 ns [24.2043, 24.2929] (0.95 confidence)</td>
<td>0.978892</td>
</tr>
<tr>
<td>angle&#95;d3/dense/construct&#95;then&#95;angle</td>
<td>22.9059 ns [22.8408, 22.9301] (0.95 confidence)</td>
<td>23.3629 ns [23.3445, 23.388] (0.95 confidence)</td>
<td>0.98044</td>
</tr>
<tr>
<td>angle&#95;d3/dense/prepared&#95;vector</td>
<td>22.5975 ns [22.5706, 22.6399] (0.95 confidence)</td>
<td>23.2547 ns [23.2199, 23.2781] (0.95 confidence)</td>
<td>0.97174</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;antipodal/borrowed&#95;slice</td>
<td>16.7075 ns [16.699, 16.7151] (0.95 confidence)</td>
<td>17.0525 ns [16.9883, 17.101] (0.95 confidence)</td>
<td>0.979769</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;antipodal/construct&#95;then&#95;angle</td>
<td>15.5421 ns [15.5342, 15.6028] (0.95 confidence)</td>
<td>16.1209 ns [16.0608, 16.1613] (0.95 confidence)</td>
<td>0.964098</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;antipodal/prepared&#95;vector</td>
<td>15.7179 ns [15.6955, 15.7351] (0.95 confidence)</td>
<td>16.1643 ns [16.1398, 16.1929] (0.95 confidence)</td>
<td>0.972386</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;identical/borrowed&#95;slice</td>
<td>14.2448 ns [14.2417, 14.2485] (0.95 confidence)</td>
<td>14.4289 ns [14.415, 14.4365] (0.95 confidence)</td>
<td>0.98724</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;identical/construct&#95;then&#95;angle</td>
<td>13.4831 ns [13.479, 13.4893] (0.95 confidence)</td>
<td>13.636 ns [13.6087, 13.6546] (0.95 confidence)</td>
<td>0.988786</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;identical/prepared&#95;vector</td>
<td>13.3974 ns [13.3869, 13.4081] (0.95 confidence)</td>
<td>13.5631 ns [13.5412, 13.5793] (0.95 confidence)</td>
<td>0.987788</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;varied/borrowed&#95;slice</td>
<td>26.0424 ns [26.0252, 26.0886] (0.95 confidence)</td>
<td>28.2035 ns [27.7448, 28.6696] (0.95 confidence)</td>
<td>0.923376</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;varied/construct&#95;then&#95;angle</td>
<td>25.0637 ns [25.0285, 25.1031] (0.95 confidence)</td>
<td>25.638 ns [25.6146, 25.6663] (0.95 confidence)</td>
<td>0.977601</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;varied/prepared&#95;vector</td>
<td>24.5745 ns [24.5599, 24.6028] (0.95 confidence)</td>
<td>25.4037 ns [25.3797, 25.4957] (0.95 confidence)</td>
<td>0.96736</td>
</tr>
<tr>
<td>angle&#95;d3/maximum&#95;small&#95;angle/borrowed&#95;slice</td>
<td>10.6224 ns [10.6028, 10.6444] (0.95 confidence)</td>
<td>11.3288 ns [11.2871, 11.3609] (0.95 confidence)</td>
<td>0.937646</td>
</tr>
<tr>
<td>angle&#95;d3/maximum&#95;small&#95;angle/construct&#95;then&#95;angle</td>
<td>9.81408 ns [9.78924, 9.85403] (0.95 confidence)</td>
<td>10.2713 ns [10.2469, 10.3088] (0.95 confidence)</td>
<td>0.955486</td>
</tr>
<tr>
<td>angle&#95;d3/maximum&#95;small&#95;angle/prepared&#95;vector</td>
<td>9.27639 ns [9.25612, 9.29112] (0.95 confidence)</td>
<td>9.62336 ns [9.52821, 9.64421] (0.95 confidence)</td>
<td>0.963945</td>
</tr>
<tr>
<td>angle&#95;d3/mixed&#95;scale/borrowed&#95;slice</td>
<td>12.0937 ns [12.085, 12.1417] (0.95 confidence)</td>
<td>12.6226 ns [12.597, 12.6841] (0.95 confidence)</td>
<td>0.958097</td>
</tr>
<tr>
<td>angle&#95;d3/mixed&#95;scale/construct&#95;then&#95;angle</td>
<td>11.177 ns [11.1699, 11.2092] (0.95 confidence)</td>
<td>11.6246 ns [11.5769, 11.6781] (0.95 confidence)</td>
<td>0.961495</td>
</tr>
<tr>
<td>angle&#95;d3/mixed&#95;scale/prepared&#95;vector</td>
<td>10.5114 ns [10.5065, 10.5271] (0.95 confidence)</td>
<td>10.8358 ns [10.8201, 10.8505] (0.95 confidence)</td>
<td>0.970063</td>
</tr>
<tr>
<td>angle&#95;d3/near&#95;antipodal/borrowed&#95;slice</td>
<td>16.335 ns [16.3271, 16.341] (0.95 confidence)</td>
<td>16.7933 ns [16.7732, 16.8245] (0.95 confidence)</td>
<td>0.97271</td>
</tr>
<tr>
<td>angle&#95;d3/near&#95;antipodal/construct&#95;then&#95;angle</td>
<td>15.0904 ns [15.0829, 15.107] (0.95 confidence)</td>
<td>15.4648 ns [15.4364, 15.4971] (0.95 confidence)</td>
<td>0.975793</td>
</tr>
<tr>
<td>angle&#95;d3/near&#95;antipodal/prepared&#95;vector</td>
<td>14.887 ns [14.8436, 14.9054] (0.95 confidence)</td>
<td>15.6669 ns [15.6225, 15.699] (0.95 confidence)</td>
<td>0.950221</td>
</tr>
<tr>
<td>angle&#95;d3/near&#95;parallel/borrowed&#95;slice</td>
<td>10.6155 ns [10.5922, 10.6994] (0.95 confidence)</td>
<td>11.1814 ns [11.1368, 11.2418] (0.95 confidence)</td>
<td>0.949388</td>
</tr>
<tr>
<td>angle&#95;d3/near&#95;parallel/construct&#95;then&#95;angle</td>
<td>9.81133 ns [9.78214, 9.85502] (0.95 confidence)</td>
<td>10.2575 ns [10.2475, 10.2986] (0.95 confidence)</td>
<td>0.956503</td>
</tr>
<tr>
<td>angle&#95;d3/near&#95;parallel/prepared&#95;vector</td>
<td>9.24811 ns [9.24016, 9.26855] (0.95 confidence)</td>
<td>9.70411 ns [9.68052, 9.71917] (0.95 confidence)</td>
<td>0.953009</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;adjacent/borrowed&#95;slice</td>
<td>13.6075 ns [13.6042, 13.6118] (0.95 confidence)</td>
<td>13.7874 ns [13.7473, 13.8264] (0.95 confidence)</td>
<td>0.986949</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;adjacent/construct&#95;then&#95;angle</td>
<td>12.9177 ns [12.8763, 13.0472] (0.95 confidence)</td>
<td>13.0549 ns [13.0356, 13.0619] (0.95 confidence)</td>
<td>0.989493</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;adjacent/prepared&#95;vector</td>
<td>11.6374 ns [11.6355, 11.6394] (0.95 confidence)</td>
<td>11.8379 ns [11.809, 11.8648] (0.95 confidence)</td>
<td>0.983064</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;alternating/borrowed&#95;slice</td>
<td>13.8569 ns [13.8521, 13.8696] (0.95 confidence)</td>
<td>14.0366 ns [13.9898, 14.0979] (0.95 confidence)</td>
<td>0.987199</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;alternating/construct&#95;then&#95;angle</td>
<td>12.8704 ns [12.8683, 12.8793] (0.95 confidence)</td>
<td>13.066 ns [13.0475, 13.1114] (0.95 confidence)</td>
<td>0.98503</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;alternating/prepared&#95;vector</td>
<td>11.8829 ns [11.8799, 11.8845] (0.95 confidence)</td>
<td>12.0505 ns [12.0323, 12.0762] (0.95 confidence)</td>
<td>0.986092</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;antipodal/borrowed&#95;slice</td>
<td>10.6428 ns [10.6404, 10.6479] (0.95 confidence)</td>
<td>11.0102 ns [10.97, 11.075] (0.95 confidence)</td>
<td>0.966626</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;antipodal/construct&#95;then&#95;angle</td>
<td>10.0213 ns [9.9085, 10.0807] (0.95 confidence)</td>
<td>10.2968 ns [10.2521, 10.3806] (0.95 confidence)</td>
<td>0.973248</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;antipodal/prepared&#95;vector</td>
<td>8.93716 ns [8.93539, 8.93965] (0.95 confidence)</td>
<td>9.11559 ns [9.07936, 9.15954] (0.95 confidence)</td>
<td>0.980426</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;identical/borrowed&#95;slice</td>
<td>8.66008 ns [8.65702, 8.66317] (0.95 confidence)</td>
<td>8.72456 ns [8.70782, 8.7582] (0.95 confidence)</td>
<td>0.992609</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;identical/construct&#95;then&#95;angle</td>
<td>7.73278 ns [7.69899, 7.86064] (0.95 confidence)</td>
<td>7.91639 ns [7.84952, 7.96057] (0.95 confidence)</td>
<td>0.976807</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;identical/prepared&#95;vector</td>
<td>6.70197 ns [6.70102, 6.70293] (0.95 confidence)</td>
<td>6.80658 ns [6.78311, 6.82855] (0.95 confidence)</td>
<td>0.984631</td>
</tr>
<tr>
<td>angle&#95;d3/subnormal/borrowed&#95;slice</td>
<td>10.6257 ns [10.5963, 10.6635] (0.95 confidence)</td>
<td>11.0446 ns [10.9908, 11.1105] (0.95 confidence)</td>
<td>0.962076</td>
</tr>
<tr>
<td>angle&#95;d3/subnormal/construct&#95;then&#95;angle</td>
<td>9.80408 ns [9.77744, 9.84474] (0.95 confidence)</td>
<td>10.1628 ns [10.1331, 10.1798] (0.95 confidence)</td>
<td>0.964701</td>
</tr>
<tr>
<td>angle&#95;d3/subnormal/prepared&#95;vector</td>
<td>9.25867 ns [9.24881, 9.27094] (0.95 confidence)</td>
<td>9.52505 ns [9.50991, 9.54521] (0.95 confidence)</td>
<td>0.972034</td>
</tr>
<tr>
<td>angle&#95;d3/subnormal&#95;axes/borrowed&#95;slice</td>
<td>13.7133 ns [13.6978, 13.7542] (0.95 confidence)</td>
<td>14.8963 ns [14.841, 15.0575] (0.95 confidence)</td>
<td>0.92058</td>
</tr>
<tr>
<td>angle&#95;d3/subnormal&#95;axes/construct&#95;then&#95;angle</td>
<td>13.1087 ns [13.0624, 13.1371] (0.95 confidence)</td>
<td>13.4295 ns [13.4052, 13.5017] (0.95 confidence)</td>
<td>0.976112</td>
</tr>
<tr>
<td>angle&#95;d3/subnormal&#95;axes/prepared&#95;vector</td>
<td>11.8819 ns [11.8802, 11.8854] (0.95 confidence)</td>
<td>12.3865 ns [12.3522, 12.4469] (0.95 confidence)</td>
<td>0.959262</td>
</tr>
<tr>
<td>angle&#95;d4/dense/borrowed&#95;slice</td>
<td>34.1112 ns [34.0272, 34.1597] (0.95 confidence)</td>
<td>35.2923 ns [35.2392, 35.3213] (0.95 confidence)</td>
<td>0.966532</td>
</tr>
<tr>
<td>angle&#95;d4/dense/construct&#95;then&#95;angle</td>
<td>32.765 ns [32.7305, 32.8156] (0.95 confidence)</td>
<td>33.4145 ns [33.3572, 33.4524] (0.95 confidence)</td>
<td>0.980564</td>
</tr>
<tr>
<td>angle&#95;d4/dense/prepared&#95;vector</td>
<td>32.3477 ns [32.2929, 32.3816] (0.95 confidence)</td>
<td>33.5556 ns [33.5063, 33.6174] (0.95 confidence)</td>
<td>0.964004</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;antipodal/borrowed&#95;slice</td>
<td>26.6095 ns [26.5816, 26.6224] (0.95 confidence)</td>
<td>27.3767 ns [27.3002, 27.4395] (0.95 confidence)</td>
<td>0.971976</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;antipodal/construct&#95;then&#95;angle</td>
<td>24.9545 ns [24.9255, 24.9788] (0.95 confidence)</td>
<td>26.0569 ns [26.0187, 26.0753] (0.95 confidence)</td>
<td>0.957691</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;antipodal/prepared&#95;vector</td>
<td>24.6866 ns [24.6212, 24.7475] (0.95 confidence)</td>
<td>26.5459 ns [26.4599, 26.6702] (0.95 confidence)</td>
<td>0.929961</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;identical/borrowed&#95;slice</td>
<td>24.1149 ns [24.0824, 24.1385] (0.95 confidence)</td>
<td>24.8486 ns [24.7793, 24.8828] (0.95 confidence)</td>
<td>0.970473</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;identical/construct&#95;then&#95;angle</td>
<td>22.5833 ns [22.5711, 22.598] (0.95 confidence)</td>
<td>24.0491 ns [23.8814, 24.221] (0.95 confidence)</td>
<td>0.93905</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;identical/prepared&#95;vector</td>
<td>22.1023 ns [22.0645, 22.1269] (0.95 confidence)</td>
<td>22.6215 ns [22.5362, 22.6778] (0.95 confidence)</td>
<td>0.977048</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;varied/borrowed&#95;slice</td>
<td>35.4566 ns [35.4083, 35.499] (0.95 confidence)</td>
<td>36.4458 ns [36.3875, 36.4837] (0.95 confidence)</td>
<td>0.972859</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;varied/construct&#95;then&#95;angle</td>
<td>33.95 ns [33.8659, 34.0147] (0.95 confidence)</td>
<td>35.3611 ns [35.2133, 35.4613] (0.95 confidence)</td>
<td>0.960096</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;varied/prepared&#95;vector</td>
<td>33.2881 ns [33.2609, 33.3277] (0.95 confidence)</td>
<td>34.1928 ns [34.1568, 34.2279] (0.95 confidence)</td>
<td>0.97354</td>
</tr>
<tr>
<td>angle&#95;d4/maximum&#95;small&#95;angle/borrowed&#95;slice</td>
<td>12.3833 ns [12.2159, 12.6031] (0.95 confidence)</td>
<td>12.6213 ns [12.5371, 12.7341] (0.95 confidence)</td>
<td>0.981138</td>
</tr>
<tr>
<td>angle&#95;d4/maximum&#95;small&#95;angle/construct&#95;then&#95;angle</td>
<td>11.1117 ns [10.9687, 11.1478] (0.95 confidence)</td>
<td>11.52 ns [11.4783, 11.5972] (0.95 confidence)</td>
<td>0.964549</td>
</tr>
<tr>
<td>angle&#95;d4/maximum&#95;small&#95;angle/prepared&#95;vector</td>
<td>9.80742 ns [9.80009, 9.81995] (0.95 confidence)</td>
<td>10.1074 ns [10.0788, 10.128] (0.95 confidence)</td>
<td>0.970322</td>
</tr>
<tr>
<td>angle&#95;d4/mixed&#95;scale/borrowed&#95;slice</td>
<td>13.9531 ns [13.8293, 14.0852] (0.95 confidence)</td>
<td>14.2986 ns [14.2565, 14.3359] (0.95 confidence)</td>
<td>0.975841</td>
</tr>
<tr>
<td>angle&#95;d4/mixed&#95;scale/construct&#95;then&#95;angle</td>
<td>12.5028 ns [12.4605, 12.5919] (0.95 confidence)</td>
<td>12.9522 ns [12.8992, 12.9874] (0.95 confidence)</td>
<td>0.965301</td>
</tr>
<tr>
<td>angle&#95;d4/mixed&#95;scale/prepared&#95;vector</td>
<td>11.0833 ns [11.0445, 11.139] (0.95 confidence)</td>
<td>11.3916 ns [11.3652, 11.4154] (0.95 confidence)</td>
<td>0.972936</td>
</tr>
<tr>
<td>angle&#95;d4/near&#95;antipodal/borrowed&#95;slice</td>
<td>17.6154 ns [17.6071, 17.624] (0.95 confidence)</td>
<td>18.0245 ns [18.0153, 18.0425] (0.95 confidence)</td>
<td>0.977306</td>
</tr>
<tr>
<td>angle&#95;d4/near&#95;antipodal/construct&#95;then&#95;angle</td>
<td>15.9822 ns [15.9488, 15.9984] (0.95 confidence)</td>
<td>16.3714 ns [16.3503, 16.3857] (0.95 confidence)</td>
<td>0.976232</td>
</tr>
<tr>
<td>angle&#95;d4/near&#95;antipodal/prepared&#95;vector</td>
<td>15.4338 ns [15.4214, 15.4551] (0.95 confidence)</td>
<td>15.8186 ns [15.8081, 15.8332] (0.95 confidence)</td>
<td>0.975674</td>
</tr>
<tr>
<td>angle&#95;d4/near&#95;parallel/borrowed&#95;slice</td>
<td>12.5011 ns [12.2832, 12.6139] (0.95 confidence)</td>
<td>12.5108 ns [12.4169, 12.6364] (0.95 confidence)</td>
<td>0.999226</td>
</tr>
<tr>
<td>angle&#95;d4/near&#95;parallel/construct&#95;then&#95;angle</td>
<td>11.1347 ns [11.0311, 11.1422] (0.95 confidence)</td>
<td>11.6198 ns [11.5692, 11.6868] (0.95 confidence)</td>
<td>0.958252</td>
</tr>
<tr>
<td>angle&#95;d4/near&#95;parallel/prepared&#95;vector</td>
<td>9.82509 ns [9.81174, 9.85677] (0.95 confidence)</td>
<td>10.1082 ns [10.0859, 10.1279] (0.95 confidence)</td>
<td>0.971989</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;adjacent/borrowed&#95;slice</td>
<td>15.6951 ns [15.5541, 15.8268] (0.95 confidence)</td>
<td>16.4505 ns [16.419, 16.4965] (0.95 confidence)</td>
<td>0.954079</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;adjacent/construct&#95;then&#95;angle</td>
<td>14.6124 ns [14.5663, 14.622] (0.95 confidence)</td>
<td>15.4333 ns [15.3991, 15.4898] (0.95 confidence)</td>
<td>0.946813</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;adjacent/prepared&#95;vector</td>
<td>13.1277 ns [13.1251, 13.1293] (0.95 confidence)</td>
<td>13.7623 ns [13.7377, 13.7935] (0.95 confidence)</td>
<td>0.953888</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;alternating/borrowed&#95;slice</td>
<td>16.0498 ns [15.9442, 16.0819] (0.95 confidence)</td>
<td>16.9195 ns [16.8609, 17.0052] (0.95 confidence)</td>
<td>0.948596</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;alternating/construct&#95;then&#95;angle</td>
<td>14.896 ns [14.8784, 14.9285] (0.95 confidence)</td>
<td>15.699 ns [15.6749, 15.728] (0.95 confidence)</td>
<td>0.948854</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;alternating/prepared&#95;vector</td>
<td>13.3783 ns [13.3755, 13.391] (0.95 confidence)</td>
<td>14.0386 ns [14.0256, 14.0754] (0.95 confidence)</td>
<td>0.952963</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;antipodal/borrowed&#95;slice</td>
<td>12.0605 ns [11.987, 12.0906] (0.95 confidence)</td>
<td>12.2886 ns [12.1627, 12.3586] (0.95 confidence)</td>
<td>0.981437</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;antipodal/construct&#95;then&#95;angle</td>
<td>10.8845 ns [10.8532, 10.8879] (0.95 confidence)</td>
<td>11.3649 ns [11.321, 11.3912] (0.95 confidence)</td>
<td>0.957729</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;antipodal/prepared&#95;vector</td>
<td>9.39887 ns [9.37988, 9.40086] (0.95 confidence)</td>
<td>9.64409 ns [9.62061, 9.6734] (0.95 confidence)</td>
<td>0.974573</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;identical/borrowed&#95;slice</td>
<td>10.3622 ns [10.2744, 10.3904] (0.95 confidence)</td>
<td>10.5313 ns [10.4544, 10.6143] (0.95 confidence)</td>
<td>0.983943</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;identical/construct&#95;then&#95;angle</td>
<td>9.15826 ns [9.15486, 9.16493] (0.95 confidence)</td>
<td>9.56566 ns [9.50978, 9.59923] (0.95 confidence)</td>
<td>0.95741</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;identical/prepared&#95;vector</td>
<td>7.67363 ns [7.67156, 7.68369] (0.95 confidence)</td>
<td>7.97824 ns [7.90135, 8.00155] (0.95 confidence)</td>
<td>0.96182</td>
</tr>
<tr>
<td>angle&#95;d4/subnormal/borrowed&#95;slice</td>
<td>12.2733 ns [12.1732, 12.3798] (0.95 confidence)</td>
<td>12.6982 ns [12.5685, 12.82] (0.95 confidence)</td>
<td>0.966544</td>
</tr>
<tr>
<td>angle&#95;d4/subnormal/construct&#95;then&#95;angle</td>
<td>10.978 ns [10.9308, 11.0859] (0.95 confidence)</td>
<td>11.509 ns [11.4468, 11.5934] (0.95 confidence)</td>
<td>0.953867</td>
</tr>
<tr>
<td>angle&#95;d4/subnormal/prepared&#95;vector</td>
<td>9.84067 ns [9.81284, 9.86688] (0.95 confidence)</td>
<td>10.1817 ns [10.1596, 10.2009] (0.95 confidence)</td>
<td>0.966508</td>
</tr>
<tr>
<td>angle&#95;d4/subnormal&#95;axes/borrowed&#95;slice</td>
<td>15.4307 ns [15.3837, 15.4998] (0.95 confidence)</td>
<td>16.4265 ns [16.3504, 16.5476] (0.95 confidence)</td>
<td>0.939379</td>
</tr>
<tr>
<td>angle&#95;d4/subnormal&#95;axes/construct&#95;then&#95;angle</td>
<td>14.2423 ns [14.1658, 14.3438] (0.95 confidence)</td>
<td>14.9556 ns [14.8987, 15.0022] (0.95 confidence)</td>
<td>0.952307</td>
</tr>
<tr>
<td>angle&#95;d4/subnormal&#95;axes/prepared&#95;vector</td>
<td>12.7102 ns [12.6564, 12.7964] (0.95 confidence)</td>
<td>13.1087 ns [13.0496, 13.1562] (0.95 confidence)</td>
<td>0.969597</td>
</tr>
<tr>
<td>angle&#95;d5/dense/borrowed&#95;slice</td>
<td>47.7304 ns [47.6782, 47.7495] (0.95 confidence)</td>
<td>49.0117 ns [48.7682, 49.816] (0.95 confidence)</td>
<td>0.973858</td>
</tr>
<tr>
<td>angle&#95;d5/dense/construct&#95;then&#95;angle</td>
<td>45.259 ns [45.2186, 45.2948] (0.95 confidence)</td>
<td>46.6861 ns [46.6082, 46.7346] (0.95 confidence)</td>
<td>0.969431</td>
</tr>
<tr>
<td>angle&#95;d5/dense/prepared&#95;vector</td>
<td>44.8206 ns [44.7447, 44.8891] (0.95 confidence)</td>
<td>46.2646 ns [46.2419, 46.2903] (0.95 confidence)</td>
<td>0.968788</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;antipodal/borrowed&#95;slice</td>
<td>36.7556 ns [36.6704, 36.7981] (0.95 confidence)</td>
<td>37.5898 ns [37.5674, 37.6239] (0.95 confidence)</td>
<td>0.977807</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;antipodal/construct&#95;then&#95;angle</td>
<td>34.8673 ns [34.8081, 34.8864] (0.95 confidence)</td>
<td>35.6135 ns [35.5775, 35.6603] (0.95 confidence)</td>
<td>0.979049</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;antipodal/prepared&#95;vector</td>
<td>34.1653 ns [34.116, 34.2351] (0.95 confidence)</td>
<td>35.0864 ns [35.0626, 35.1068] (0.95 confidence)</td>
<td>0.973749</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;identical/borrowed&#95;slice</td>
<td>35.3528 ns [35.2928, 35.4121] (0.95 confidence)</td>
<td>37.0795 ns [37.0255, 37.1319] (0.95 confidence)</td>
<td>0.953432</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;identical/construct&#95;then&#95;angle</td>
<td>33.4438 ns [33.3743, 33.5277] (0.95 confidence)</td>
<td>36.436 ns [34.6168, 37.9745] (0.95 confidence)</td>
<td>0.917877</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;identical/prepared&#95;vector</td>
<td>32.9229 ns [32.875, 32.9473] (0.95 confidence)</td>
<td>34.6725 ns [34.543, 34.7711] (0.95 confidence)</td>
<td>0.94954</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;varied/borrowed&#95;slice</td>
<td>68.3369 ns [68.1146, 68.5017] (0.95 confidence)</td>
<td>73.7419 ns [73.5926, 73.9405] (0.95 confidence)</td>
<td>0.926704</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;varied/construct&#95;then&#95;angle</td>
<td>65.9484 ns [65.7368, 66.1148] (0.95 confidence)</td>
<td>70.5048 ns [70.3532, 70.6824] (0.95 confidence)</td>
<td>0.935375</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;varied/prepared&#95;vector</td>
<td>63.7527 ns [63.5641, 63.8941] (0.95 confidence)</td>
<td>68.7437 ns [68.5983, 68.9136] (0.95 confidence)</td>
<td>0.927396</td>
</tr>
<tr>
<td>angle&#95;d5/maximum&#95;small&#95;angle/borrowed&#95;slice</td>
<td>14.3485 ns [14.333, 14.3548] (0.95 confidence)</td>
<td>14.6416 ns [14.5751, 14.6737] (0.95 confidence)</td>
<td>0.979981</td>
</tr>
<tr>
<td>angle&#95;d5/maximum&#95;small&#95;angle/construct&#95;then&#95;angle</td>
<td>12.6437 ns [12.6157, 12.6519] (0.95 confidence)</td>
<td>13.0445 ns [12.9738, 13.1265] (0.95 confidence)</td>
<td>0.969272</td>
</tr>
<tr>
<td>angle&#95;d5/maximum&#95;small&#95;angle/prepared&#95;vector</td>
<td>10.909 ns [10.8528, 10.9125] (0.95 confidence)</td>
<td>11.1675 ns [11.1311, 11.2097] (0.95 confidence)</td>
<td>0.97685</td>
</tr>
<tr>
<td>angle&#95;d5/mixed&#95;scale/borrowed&#95;slice</td>
<td>15.8392 ns [15.8353, 15.8448] (0.95 confidence)</td>
<td>16.4862 ns [16.2204, 16.6951] (0.95 confidence)</td>
<td>0.960752</td>
</tr>
<tr>
<td>angle&#95;d5/mixed&#95;scale/construct&#95;then&#95;angle</td>
<td>13.9796 ns [13.9231, 14.0892] (0.95 confidence)</td>
<td>15.2285 ns [14.9571, 15.3378] (0.95 confidence)</td>
<td>0.917988</td>
</tr>
<tr>
<td>angle&#95;d5/mixed&#95;scale/prepared&#95;vector</td>
<td>12.3928 ns [12.3084, 12.409] (0.95 confidence)</td>
<td>12.8407 ns [12.7754, 12.9028] (0.95 confidence)</td>
<td>0.965115</td>
</tr>
<tr>
<td>angle&#95;d5/near&#95;antipodal/borrowed&#95;slice</td>
<td>19.1715 ns [19.1501, 19.1907] (0.95 confidence)</td>
<td>19.6337 ns [19.5726, 19.709] (0.95 confidence)</td>
<td>0.976456</td>
</tr>
<tr>
<td>angle&#95;d5/near&#95;antipodal/construct&#95;then&#95;angle</td>
<td>16.9704 ns [16.9571, 16.9882] (0.95 confidence)</td>
<td>17.369 ns [17.3305, 17.4] (0.95 confidence)</td>
<td>0.97705</td>
</tr>
<tr>
<td>angle&#95;d5/near&#95;antipodal/prepared&#95;vector</td>
<td>16.2226 ns [16.2117, 16.2466] (0.95 confidence)</td>
<td>16.7474 ns [16.7114, 16.7952] (0.95 confidence)</td>
<td>0.968662</td>
</tr>
<tr>
<td>angle&#95;d5/near&#95;parallel/borrowed&#95;slice</td>
<td>14.3485 ns [14.3379, 14.351] (0.95 confidence)</td>
<td>14.6321 ns [14.5564, 14.67] (0.95 confidence)</td>
<td>0.980617</td>
</tr>
<tr>
<td>angle&#95;d5/near&#95;parallel/construct&#95;then&#95;angle</td>
<td>12.6465 ns [12.6142, 12.7352] (0.95 confidence)</td>
<td>13.077 ns [13.0039, 13.1558] (0.95 confidence)</td>
<td>0.967077</td>
</tr>
<tr>
<td>angle&#95;d5/near&#95;parallel/prepared&#95;vector</td>
<td>10.9096 ns [10.9014, 10.9135] (0.95 confidence)</td>
<td>11.1508 ns [11.1308, 11.2121] (0.95 confidence)</td>
<td>0.978367</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;adjacent/borrowed&#95;slice</td>
<td>17.3375 ns [17.3342, 17.3501] (0.95 confidence)</td>
<td>18.0094 ns [17.9355, 18.0623] (0.95 confidence)</td>
<td>0.962692</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;adjacent/construct&#95;then&#95;angle</td>
<td>16.1735 ns [16.1504, 16.3043] (0.95 confidence)</td>
<td>16.7345 ns [16.7027, 16.8015] (0.95 confidence)</td>
<td>0.966476</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;adjacent/prepared&#95;vector</td>
<td>14.392 ns [14.389, 14.3958] (0.95 confidence)</td>
<td>14.7096 ns [14.6788, 14.7407] (0.95 confidence)</td>
<td>0.978407</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;alternating/borrowed&#95;slice</td>
<td>27.4846 ns [27.4434, 27.4995] (0.95 confidence)</td>
<td>28.0515 ns [28.0408, 28.073] (0.95 confidence)</td>
<td>0.979792</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;alternating/construct&#95;then&#95;angle</td>
<td>26.4642 ns [26.4289, 26.4762] (0.95 confidence)</td>
<td>27.0643 ns [27.0217, 27.0926] (0.95 confidence)</td>
<td>0.977826</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;alternating/prepared&#95;vector</td>
<td>25.9338 ns [25.8944, 25.9527] (0.95 confidence)</td>
<td>26.4721 ns [26.452, 26.4797] (0.95 confidence)</td>
<td>0.979667</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;antipodal/borrowed&#95;slice</td>
<td>13.11 ns [13.0904, 13.1129] (0.95 confidence)</td>
<td>13.6507 ns [13.5909, 13.706] (0.95 confidence)</td>
<td>0.960394</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;antipodal/construct&#95;then&#95;angle</td>
<td>11.9183 ns [11.9097, 11.9241] (0.95 confidence)</td>
<td>12.3146 ns [12.2705, 12.3609] (0.95 confidence)</td>
<td>0.967816</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;antipodal/prepared&#95;vector</td>
<td>10.1695 ns [10.157, 10.1749] (0.95 confidence)</td>
<td>10.5476 ns [10.5195, 10.5693] (0.95 confidence)</td>
<td>0.964149</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;identical/borrowed&#95;slice</td>
<td>11.6382 ns [11.6318, 11.6543] (0.95 confidence)</td>
<td>12.0281 ns [12.0022, 12.0856] (0.95 confidence)</td>
<td>0.967585</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;identical/construct&#95;then&#95;angle</td>
<td>10.4183 ns [10.407, 10.4639] (0.95 confidence)</td>
<td>10.905 ns [10.8153, 10.9394] (0.95 confidence)</td>
<td>0.955369</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;identical/prepared&#95;vector</td>
<td>8.69878 ns [8.69016, 8.71866] (0.95 confidence)</td>
<td>8.89129 ns [8.87205, 8.9131] (0.95 confidence)</td>
<td>0.978348</td>
</tr>
<tr>
<td>angle&#95;d5/subnormal/borrowed&#95;slice</td>
<td>14.3596 ns [14.3513, 14.3947] (0.95 confidence)</td>
<td>14.6979 ns [14.6963, 14.7016] (0.95 confidence)</td>
<td>0.976984</td>
</tr>
<tr>
<td>angle&#95;d5/subnormal/construct&#95;then&#95;angle</td>
<td>12.6491 ns [12.5561, 12.7818] (0.95 confidence)</td>
<td>12.7935 ns [12.6762, 12.7977] (0.95 confidence)</td>
<td>0.98871</td>
</tr>
<tr>
<td>angle&#95;d5/subnormal/prepared&#95;vector</td>
<td>10.9153 ns [10.9112, 10.9295] (0.95 confidence)</td>
<td>11.1781 ns [11.1763, 11.1805] (0.95 confidence)</td>
<td>0.976496</td>
</tr>
<tr>
<td>angle&#95;d5/subnormal&#95;axes/borrowed&#95;slice</td>
<td>17.5305 ns [17.4984, 17.5667] (0.95 confidence)</td>
<td>17.8659 ns [17.8088, 17.9238] (0.95 confidence)</td>
<td>0.981227</td>
</tr>
<tr>
<td>angle&#95;d5/subnormal&#95;axes/construct&#95;then&#95;angle</td>
<td>15.7942 ns [15.6554, 15.8358] (0.95 confidence)</td>
<td>16.3094 ns [16.2388, 16.4061] (0.95 confidence)</td>
<td>0.968415</td>
</tr>
<tr>
<td>angle&#95;d5/subnormal&#95;axes/prepared&#95;vector</td>
<td>14.0864 ns [14.0335, 14.1293] (0.95 confidence)</td>
<td>14.4492 ns [14.4031, 14.5307] (0.95 confidence)</td>
<td>0.974894</td>
</tr>
<tr>
<td>angle&#95;d6/dense/borrowed&#95;slice</td>
<td>62.6908 ns [62.5952, 62.8393] (0.95 confidence)</td>
<td>63.6192 ns [63.6103, 63.6326] (0.95 confidence)</td>
<td>0.985407</td>
</tr>
<tr>
<td>angle&#95;d6/dense/construct&#95;then&#95;angle</td>
<td>59.5797 ns [59.497, 59.6655] (0.95 confidence)</td>
<td>60.9019 ns [60.883, 60.9241] (0.95 confidence)</td>
<td>0.978289</td>
</tr>
<tr>
<td>angle&#95;d6/dense/prepared&#95;vector</td>
<td>59.1339 ns [59.0133, 59.2115] (0.95 confidence)</td>
<td>60.485 ns [60.4699, 60.5] (0.95 confidence)</td>
<td>0.977663</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;antipodal/borrowed&#95;slice</td>
<td>48.8058 ns [48.68, 48.8977] (0.95 confidence)</td>
<td>49.6328 ns [49.6248, 49.6381] (0.95 confidence)</td>
<td>0.983337</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;antipodal/construct&#95;then&#95;angle</td>
<td>46.2381 ns [46.1285, 46.2824] (0.95 confidence)</td>
<td>46.878 ns [46.8698, 46.895] (0.95 confidence)</td>
<td>0.986351</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;antipodal/prepared&#95;vector</td>
<td>45.6532 ns [45.6045, 45.7289] (0.95 confidence)</td>
<td>46.5688 ns [46.5555, 46.5813] (0.95 confidence)</td>
<td>0.980338</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;identical/borrowed&#95;slice</td>
<td>47.6671 ns [47.57, 47.7333] (0.95 confidence)</td>
<td>48.62 ns [48.5994, 48.6354] (0.95 confidence)</td>
<td>0.980401</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;identical/construct&#95;then&#95;angle</td>
<td>45.1697 ns [45.0931, 45.2548] (0.95 confidence)</td>
<td>45.89 ns [45.8772, 45.9045] (0.95 confidence)</td>
<td>0.984303</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;identical/prepared&#95;vector</td>
<td>44.5166 ns [44.458, 44.5632] (0.95 confidence)</td>
<td>45.4069 ns [45.3947, 45.4173] (0.95 confidence)</td>
<td>0.980392</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;varied/borrowed&#95;slice</td>
<td>93.4086 ns [93.1297, 93.5652] (0.95 confidence)</td>
<td>98.524 ns [98.4684, 98.5973] (0.95 confidence)</td>
<td>0.94808</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;varied/construct&#95;then&#95;angle</td>
<td>91.0468 ns [90.9523, 91.1804] (0.95 confidence)</td>
<td>96.0863 ns [96.0167, 96.1964] (0.95 confidence)</td>
<td>0.947553</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;varied/prepared&#95;vector</td>
<td>89.5966 ns [89.3026, 89.9077] (0.95 confidence)</td>
<td>94.9694 ns [94.9295, 95.0057] (0.95 confidence)</td>
<td>0.943425</td>
</tr>
<tr>
<td>angle&#95;d6/maximum&#95;small&#95;angle/borrowed&#95;slice</td>
<td>15.615 ns [15.6052, 15.6261] (0.95 confidence)</td>
<td>22.4492 ns [21.7679, 23.0512] (0.95 confidence)</td>
<td>0.69557</td>
</tr>
<tr>
<td>angle&#95;d6/maximum&#95;small&#95;angle/construct&#95;then&#95;angle</td>
<td>13.6325 ns [13.629, 13.6373] (0.95 confidence)</td>
<td>22.0261 ns [21.4398, 22.5548] (0.95 confidence)</td>
<td>0.618925</td>
</tr>
<tr>
<td>angle&#95;d6/maximum&#95;small&#95;angle/prepared&#95;vector</td>
<td>11.9265 ns [11.9184, 11.9433] (0.95 confidence)</td>
<td>17.1102 ns [16.6218, 17.3862] (0.95 confidence)</td>
<td>0.69704</td>
</tr>
<tr>
<td>angle&#95;d6/mixed&#95;scale/borrowed&#95;slice</td>
<td>17.0742 ns [17.0695, 17.0758] (0.95 confidence)</td>
<td>17.2777 ns [17.2631, 17.2952] (0.95 confidence)</td>
<td>0.988222</td>
</tr>
<tr>
<td>angle&#95;d6/mixed&#95;scale/construct&#95;then&#95;angle</td>
<td>15.2204 ns [15.1461, 15.3123] (0.95 confidence)</td>
<td>26.3898 ns [25.7243, 26.9199] (0.95 confidence)</td>
<td>0.576753</td>
</tr>
<tr>
<td>angle&#95;d6/mixed&#95;scale/prepared&#95;vector</td>
<td>13.3862 ns [13.3802, 13.3912] (0.95 confidence)</td>
<td>15.9598 ns [15.1965, 18.7087] (0.95 confidence)</td>
<td>0.83874</td>
</tr>
<tr>
<td>angle&#95;d6/near&#95;antipodal/borrowed&#95;slice</td>
<td>20.745 ns [20.729, 20.7591] (0.95 confidence)</td>
<td>21.9449 ns [21.3167, 24.2397] (0.95 confidence)</td>
<td>0.945322</td>
</tr>
<tr>
<td>angle&#95;d6/near&#95;antipodal/construct&#95;then&#95;angle</td>
<td>18.3742 ns [18.3694, 18.3831] (0.95 confidence)</td>
<td>31.4305 ns [28.7281, 32.4116] (0.95 confidence)</td>
<td>0.584599</td>
</tr>
<tr>
<td>angle&#95;d6/near&#95;antipodal/prepared&#95;vector</td>
<td>17.033 ns [16.9999, 17.0637] (0.95 confidence)</td>
<td>28.3067 ns [27.6686, 28.9515] (0.95 confidence)</td>
<td>0.601729</td>
</tr>
<tr>
<td>angle&#95;d6/near&#95;parallel/borrowed&#95;slice</td>
<td>15.6082 ns [15.602, 15.6153] (0.95 confidence)</td>
<td>24.3517 ns [23.9004, 24.9239] (0.95 confidence)</td>
<td>0.640947</td>
</tr>
<tr>
<td>angle&#95;d6/near&#95;parallel/construct&#95;then&#95;angle</td>
<td>13.6343 ns [13.6277, 13.6395] (0.95 confidence)</td>
<td>23.477 ns [22.736, 23.9365] (0.95 confidence)</td>
<td>0.580752</td>
</tr>
<tr>
<td>angle&#95;d6/near&#95;parallel/prepared&#95;vector</td>
<td>11.9254 ns [11.9182, 11.935] (0.95 confidence)</td>
<td>19.0504 ns [18.6899, 19.439] (0.95 confidence)</td>
<td>0.625992</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;adjacent/borrowed&#95;slice</td>
<td>19.0483 ns [18.8961, 19.1976] (0.95 confidence)</td>
<td>19.2255 ns [19.1558, 19.3717] (0.95 confidence)</td>
<td>0.990782</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;adjacent/construct&#95;then&#95;angle</td>
<td>17.9182 ns [17.8398, 17.9922] (0.95 confidence)</td>
<td>19.5395 ns [18.8553, 20.7] (0.95 confidence)</td>
<td>0.917024</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;adjacent/prepared&#95;vector</td>
<td>15.7083 ns [15.6687, 15.806] (0.95 confidence)</td>
<td>15.8973 ns [15.8619, 15.9319] (0.95 confidence)</td>
<td>0.988112</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;alternating/borrowed&#95;slice</td>
<td>29.0321 ns [29.003, 29.048] (0.95 confidence)</td>
<td>30.6755 ns [29.4883, 35.1567] (0.95 confidence)</td>
<td>0.946426</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;alternating/construct&#95;then&#95;angle</td>
<td>27.5696 ns [27.5279, 27.5943] (0.95 confidence)</td>
<td>34.0764 ns [32.4428, 34.8887] (0.95 confidence)</td>
<td>0.809052</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;alternating/prepared&#95;vector</td>
<td>26.839 ns [26.8054, 26.8686] (0.95 confidence)</td>
<td>27.774 ns [27.5763, 27.8729] (0.95 confidence)</td>
<td>0.966334</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;antipodal/borrowed&#95;slice</td>
<td>14.101 ns [14.0848, 14.1045] (0.95 confidence)</td>
<td>17.9408 ns [17.8557, 18.4677] (0.95 confidence)</td>
<td>0.785971</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;antipodal/construct&#95;then&#95;angle</td>
<td>12.6237 ns [12.6198, 12.6286] (0.95 confidence)</td>
<td>16.7434 ns [16.2085, 17.219] (0.95 confidence)</td>
<td>0.753949</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;antipodal/prepared&#95;vector</td>
<td>10.8908 ns [10.88, 10.8937] (0.95 confidence)</td>
<td>12.5617 ns [12.0246, 12.9362] (0.95 confidence)</td>
<td>0.866984</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;identical/borrowed&#95;slice</td>
<td>13.1415 ns [12.9498, 13.3283] (0.95 confidence)</td>
<td>13.6352 ns [13.5192, 13.7575] (0.95 confidence)</td>
<td>0.963795</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;identical/construct&#95;then&#95;angle</td>
<td>11.3861 ns [11.2427, 11.4551] (0.95 confidence)</td>
<td>12.0242 ns [11.9911, 12.0343] (0.95 confidence)</td>
<td>0.946938</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;identical/prepared&#95;vector</td>
<td>9.42047 ns [9.41676, 9.44686] (0.95 confidence)</td>
<td>9.65863 ns [9.64748, 9.68789] (0.95 confidence)</td>
<td>0.975343</td>
</tr>
<tr>
<td>angle&#95;d6/subnormal/borrowed&#95;slice</td>
<td>15.588 ns [15.5833, 15.5903] (0.95 confidence)</td>
<td>15.8003 ns [15.7877, 15.8149] (0.95 confidence)</td>
<td>0.986562</td>
</tr>
<tr>
<td>angle&#95;d6/subnormal/construct&#95;then&#95;angle</td>
<td>13.6297 ns [13.6235, 13.6358] (0.95 confidence)</td>
<td>14.031 ns [14.025, 14.0566] (0.95 confidence)</td>
<td>0.971403</td>
</tr>
<tr>
<td>angle&#95;d6/subnormal/prepared&#95;vector</td>
<td>11.9221 ns [11.9163, 11.9282] (0.95 confidence)</td>
<td>12.1649 ns [12.1059, 12.2039] (0.95 confidence)</td>
<td>0.980044</td>
</tr>
<tr>
<td>angle&#95;d6/subnormal&#95;axes/borrowed&#95;slice</td>
<td>18.7987 ns [18.7847, 18.805] (0.95 confidence)</td>
<td>24.8784 ns [24.1527, 25.5673] (0.95 confidence)</td>
<td>0.755625</td>
</tr>
<tr>
<td>angle&#95;d6/subnormal&#95;axes/construct&#95;then&#95;angle</td>
<td>16.862 ns [16.8579, 16.88] (0.95 confidence)</td>
<td>23.6125 ns [22.9664, 24.0363] (0.95 confidence)</td>
<td>0.714115</td>
</tr>
<tr>
<td>angle&#95;d6/subnormal&#95;axes/prepared&#95;vector</td>
<td>15.1129 ns [15.1092, 15.1158] (0.95 confidence)</td>
<td>20.315 ns [19.1786, 20.6758] (0.95 confidence)</td>
<td>0.74393</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;det</td>
<td>377.552 ns [369.134, 381.492] (0.95 confidence)</td>
<td>384.27 ns [378.618, 389.466] (0.95 confidence)</td>
<td>0.982516</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>2.89789 ns [2.89298, 2.90051] (0.95 confidence)</td>
<td>2.95155 ns [2.9489, 2.95342] (0.95 confidence)</td>
<td>0.981821</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>3.28393 ns [3.28235, 3.28482] (0.95 confidence)</td>
<td>3.34751 ns [3.34146, 3.35409] (0.95 confidence)</td>
<td>0.981007</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>457.836 ns [455.295, 482.884] (0.95 confidence)</td>
<td>503.216 ns [501.575, 504.248] (0.95 confidence)</td>
<td>0.90982</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;dot</td>
<td>2.33086 ns [2.32773, 2.33171] (0.95 confidence)</td>
<td>2.36913 ns [2.36364, 2.37451] (0.95 confidence)</td>
<td>0.983847</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;inf&#95;norm</td>
<td>32.3401 ns [32.2999, 32.386] (0.95 confidence)</td>
<td>32.9401 ns [32.9158, 32.9977] (0.95 confidence)</td>
<td>0.981786</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;ldlt</td>
<td>341.059 ns [340.752, 341.495] (0.95 confidence)</td>
<td>347.824 ns [346.943, 348.629] (0.95 confidence)</td>
<td>0.980549</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;ldlt&#95;solve</td>
<td>382.436 ns [382.036, 382.662] (0.95 confidence)</td>
<td>393.872 ns [391.402, 397.036] (0.95 confidence)</td>
<td>0.970965</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;lu</td>
<td>385.224 ns [375.686, 390.708] (0.95 confidence)</td>
<td>400.282 ns [390.349, 403.558] (0.95 confidence)</td>
<td>0.96238</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;norm2</td>
<td>10.0053 ns [9.97855, 10.0311] (0.95 confidence)</td>
<td>10.6663 ns [10.563, 10.7015] (0.95 confidence)</td>
<td>0.938031</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>9.39719 ns [9.37267, 9.42305] (0.95 confidence)</td>
<td>9.84612 ns [9.8314, 9.86163] (0.95 confidence)</td>
<td>0.954405</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>9.39775 ns [9.37432, 9.41287] (0.95 confidence)</td>
<td>9.76235 ns [9.7452, 9.78115] (0.95 confidence)</td>
<td>0.962653</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>5.00698 ns [4.99314, 5.02267] (0.95 confidence)</td>
<td>4.53038 ns [4.51313, 4.66966] (0.95 confidence)</td>
<td>1.1052</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>8.25053 ns [8.23858, 8.26485] (0.95 confidence)</td>
<td>8.57718 ns [8.55043, 8.59331] (0.95 confidence)</td>
<td>0.961916</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;norm2&#95;sq</td>
<td>2.09408 ns [2.09218, 2.09659] (0.95 confidence)</td>
<td>2.11733 ns [2.11321, 2.13171] (0.95 confidence)</td>
<td>0.98902</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>27.7519 ns [27.7408, 27.761] (0.95 confidence)</td>
<td>28.7559 ns [28.4739, 28.8898] (0.95 confidence)</td>
<td>0.965085</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>190.232 ns [189.938, 190.521] (0.95 confidence)</td>
<td>209.303 ns [208.619, 209.87] (0.95 confidence)</td>
<td>0.908886</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;det</td>
<td>0.572129 ns [0.572047, 0.572242] (0.95 confidence)</td>
<td>0.572004 ns [0.571849, 0.572147] (0.95 confidence)</td>
<td>1.00022</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>0.495541 ns [0.495467, 0.495768] (0.95 confidence)</td>
<td>0.504191 ns [0.502319, 0.506689] (0.95 confidence)</td>
<td>0.982845</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>0.519938 ns [0.51986, 0.519986] (0.95 confidence)</td>
<td>0.519812 ns [0.519733, 0.519918] (0.95 confidence)</td>
<td>1.00024</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>1.58622 ns [1.5846, 1.60619] (0.95 confidence)</td>
<td>1.58676 ns [1.58516, 1.58846] (0.95 confidence)</td>
<td>0.99966</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;dot</td>
<td>0.612704 ns [0.612144, 0.61295] (0.95 confidence)</td>
<td>0.626202 ns [0.623019, 0.627986] (0.95 confidence)</td>
<td>0.978444</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;inf&#95;norm</td>
<td>0.702335 ns [0.702282, 0.702549] (0.95 confidence)</td>
<td>0.702802 ns [0.702537, 0.703594] (0.95 confidence)</td>
<td>0.999336</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;ldlt</td>
<td>6.90258 ns [6.85292, 6.94268] (0.95 confidence)</td>
<td>6.49491 ns [6.48933, 6.50544] (0.95 confidence)</td>
<td>1.06277</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;ldlt&#95;solve</td>
<td>9.83535 ns [9.81688, 9.86058] (0.95 confidence)</td>
<td>9.31411 ns [9.29477, 9.33318] (0.95 confidence)</td>
<td>1.05596</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;lu</td>
<td>1.38617 ns [1.38358, 1.38857] (0.95 confidence)</td>
<td>1.39436 ns [1.3897, 1.39574] (0.95 confidence)</td>
<td>0.994124</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;norm2</td>
<td>1.15197 ns [1.15148, 1.15227] (0.95 confidence)</td>
<td>1.27712 ns [1.27189, 1.28598] (0.95 confidence)</td>
<td>0.902005</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>1.34428 ns [1.34261, 1.34501] (0.95 confidence)</td>
<td>1.42106 ns [1.4141, 1.43729] (0.95 confidence)</td>
<td>0.945969</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>1.34152 ns [1.33971, 1.34337] (0.95 confidence)</td>
<td>1.4088 ns [1.40576, 1.41046] (0.95 confidence)</td>
<td>0.952242</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>1.1562 ns [1.15207, 1.20297] (0.95 confidence)</td>
<td>1.19821 ns [1.16722, 1.21505] (0.95 confidence)</td>
<td>0.964946</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>1.3408 ns [1.33743, 1.34337] (0.95 confidence)</td>
<td>1.36057 ns [1.36021, 1.36097] (0.95 confidence)</td>
<td>0.985464</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;norm2&#95;sq</td>
<td>0.39448 ns [0.394295, 0.394779] (0.95 confidence)</td>
<td>0.40369 ns [0.40223, 0.404784] (0.95 confidence)</td>
<td>0.977186</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>1.23714 ns [1.23696, 1.23734] (0.95 confidence)</td>
<td>1.23537 ns [1.23445, 1.23653] (0.95 confidence)</td>
<td>1.00143</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>1.24733 ns [1.24271, 1.25286] (0.95 confidence)</td>
<td>1.2568 ns [1.25198, 1.25986] (0.95 confidence)</td>
<td>0.992471</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;det</td>
<td>1.33627 ns [1.33579, 1.33875] (0.95 confidence)</td>
<td>1.33861 ns [1.33572, 1.3447] (0.95 confidence)</td>
<td>0.998251</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>0.59397 ns [0.593943, 0.594031] (0.95 confidence)</td>
<td>0.59407 ns [0.593991, 0.59419] (0.95 confidence)</td>
<td>0.999833</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>0.692998 ns [0.692906, 0.693114] (0.95 confidence)</td>
<td>0.693101 ns [0.692993, 0.693362] (0.95 confidence)</td>
<td>0.999851</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>8.69851 ns [8.67891, 8.72332] (0.95 confidence)</td>
<td>8.7203 ns [8.70739, 8.73362] (0.95 confidence)</td>
<td>0.997502</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;dot</td>
<td>0.668063 ns [0.66694, 0.668358] (0.95 confidence)</td>
<td>0.665927 ns [0.66495, 0.667496] (0.95 confidence)</td>
<td>1.00321</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;inf&#95;norm</td>
<td>1.36214 ns [1.36195, 1.36239] (0.95 confidence)</td>
<td>1.36953 ns [1.3691, 1.36996] (0.95 confidence)</td>
<td>0.994606</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;ldlt</td>
<td>8.48654 ns [8.48098, 8.49401] (0.95 confidence)</td>
<td>8.50315 ns [8.49465, 8.51539] (0.95 confidence)</td>
<td>0.998047</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;ldlt&#95;solve</td>
<td>9.34543 ns [9.3374, 9.35461] (0.95 confidence)</td>
<td>9.37449 ns [9.34477, 9.39446] (0.95 confidence)</td>
<td>0.996901</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;lu</td>
<td>8.34886 ns [8.33511, 8.36538] (0.95 confidence)</td>
<td>8.42532 ns [8.41164, 8.4395] (0.95 confidence)</td>
<td>0.990924</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;norm2</td>
<td>1.77727 ns [1.75554, 1.77813] (0.95 confidence)</td>
<td>1.81397 ns [1.81064, 1.81906] (0.95 confidence)</td>
<td>0.979772</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>1.83517 ns [1.83363, 1.83686] (0.95 confidence)</td>
<td>1.86877 ns [1.86512, 1.87052] (0.95 confidence)</td>
<td>0.98202</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>1.83156 ns [1.82869, 1.83459] (0.95 confidence)</td>
<td>1.8654 ns [1.86386, 1.86801] (0.95 confidence)</td>
<td>0.981863</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>1.38058 ns [1.37903, 1.41518] (0.95 confidence)</td>
<td>1.38698 ns [1.38515, 1.39801] (0.95 confidence)</td>
<td>0.995386</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>1.83511 ns [1.83199, 1.8376] (0.95 confidence)</td>
<td>1.86052 ns [1.85567, 1.86395] (0.95 confidence)</td>
<td>0.986343</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;norm2&#95;sq</td>
<td>0.443882 ns [0.443573, 0.444035] (0.95 confidence)</td>
<td>0.446768 ns [0.446336, 0.447862] (0.95 confidence)</td>
<td>0.993541</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>1.77718 ns [1.77671, 1.77769] (0.95 confidence)</td>
<td>1.77785 ns [1.77736, 1.77826] (0.95 confidence)</td>
<td>0.99962</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>2.04115 ns [2.04037, 2.04194] (0.95 confidence)</td>
<td>2.05257 ns [2.05105, 2.0544] (0.95 confidence)</td>
<td>0.994436</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;det</td>
<td>1999.6 ns [1996.89, 2003.18] (0.95 confidence)</td>
<td>2032.64 ns [2028.65, 2036.19] (0.95 confidence)</td>
<td>0.983749</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>7.76642 ns [7.73965, 7.79611] (0.95 confidence)</td>
<td>7.96098 ns [7.95055, 7.98124] (0.95 confidence)</td>
<td>0.975561</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>9.08697 ns [9.06991, 9.10132] (0.95 confidence)</td>
<td>9.20112 ns [9.18527, 9.22044] (0.95 confidence)</td>
<td>0.987594</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>2033.79 ns [2030.79, 2037.03] (0.95 confidence)</td>
<td>2072.4 ns [2065.2, 2078.18] (0.95 confidence)</td>
<td>0.981367</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;dot</td>
<td>4.0178 ns [4.01143, 4.02364] (0.95 confidence)</td>
<td>4.05049 ns [4.04806, 4.05348] (0.95 confidence)</td>
<td>0.991929</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;inf&#95;norm</td>
<td>128.016 ns [127.797, 128.306] (0.95 confidence)</td>
<td>131.036 ns [130.868, 131.185] (0.95 confidence)</td>
<td>0.976952</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;ldlt</td>
<td>2392.43 ns [2391.22, 2393.31] (0.95 confidence)</td>
<td>2421.9 ns [2420.01, 2425.3] (0.95 confidence)</td>
<td>0.987835</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;ldlt&#95;solve</td>
<td>2721.35 ns [2715.58, 2725.06] (0.95 confidence)</td>
<td>2759.58 ns [2750.37, 2764.06] (0.95 confidence)</td>
<td>0.986148</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;lu</td>
<td>2192.82 ns [2186.28, 2203.51] (0.95 confidence)</td>
<td>2205.48 ns [2199.23, 2214.68] (0.95 confidence)</td>
<td>0.994263</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;norm2</td>
<td>23.1297 ns [23.0296, 23.2135] (0.95 confidence)</td>
<td>23.7672 ns [23.7158, 23.8188] (0.95 confidence)</td>
<td>0.973178</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>18.6041 ns [18.5503, 18.6143] (0.95 confidence)</td>
<td>18.9695 ns [18.9385, 19.0084] (0.95 confidence)</td>
<td>0.980737</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>18.4847 ns [18.4586, 18.527] (0.95 confidence)</td>
<td>18.9715 ns [18.9246, 19.0113] (0.95 confidence)</td>
<td>0.974342</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>9.22776 ns [9.20786, 9.23135] (0.95 confidence)</td>
<td>9.29407 ns [9.26704, 9.42591] (0.95 confidence)</td>
<td>0.992866</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>16.1531 ns [16.1229, 16.199] (0.95 confidence)</td>
<td>16.5877 ns [16.5422, 16.6172] (0.95 confidence)</td>
<td>0.973795</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;norm2&#95;sq</td>
<td>3.96809 ns [3.96021, 3.97316] (0.95 confidence)</td>
<td>4.00447 ns [4.00306, 4.00849] (0.95 confidence)</td>
<td>0.990916</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>311.571 ns [310.672, 311.949] (0.95 confidence)</td>
<td>316.057 ns [315.457, 316.553] (0.95 confidence)</td>
<td>0.985804</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>663.431 ns [662.498, 665.227] (0.95 confidence)</td>
<td>690.165 ns [687.785, 692.693] (0.95 confidence)</td>
<td>0.961264</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;det</td>
<td>2.54938 ns [2.54858, 2.55168] (0.95 confidence)</td>
<td>2.54816 ns [2.54713, 2.54944] (0.95 confidence)</td>
<td>1.00048</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>0.767163 ns [0.767108, 0.767247] (0.95 confidence)</td>
<td>0.767298 ns [0.76723, 0.767362] (0.95 confidence)</td>
<td>0.999825</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>0.866288 ns [0.866224, 0.866348] (0.95 confidence)</td>
<td>0.871575 ns [0.871106, 0.872237] (0.95 confidence)</td>
<td>0.993934</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>14.5497 ns [14.5264, 14.5905] (0.95 confidence)</td>
<td>14.5341 ns [14.5172, 14.5727] (0.95 confidence)</td>
<td>1.00107</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;dot</td>
<td>0.681411 ns [0.679407, 0.68373] (0.95 confidence)</td>
<td>0.678744 ns [0.678199, 0.681011] (0.95 confidence)</td>
<td>1.00393</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;inf&#95;norm</td>
<td>2.35196 ns [2.3504, 2.3522] (0.95 confidence)</td>
<td>2.35155 ns [2.35106, 2.35182] (0.95 confidence)</td>
<td>1.00017</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;ldlt</td>
<td>18.0875 ns [18.0585, 18.1374] (0.95 confidence)</td>
<td>18.1026 ns [18.0691, 18.1416] (0.95 confidence)</td>
<td>0.999168</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;ldlt&#95;solve</td>
<td>16.5367 ns [16.5277, 16.5412] (0.95 confidence)</td>
<td>16.5718 ns [16.514, 16.6066] (0.95 confidence)</td>
<td>0.99788</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;lu</td>
<td>14.2122 ns [14.1965, 14.2342] (0.95 confidence)</td>
<td>14.1977 ns [14.1941, 14.2113] (0.95 confidence)</td>
<td>1.00102</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;norm2</td>
<td>2.46706 ns [2.40526, 2.48349] (0.95 confidence)</td>
<td>2.2685 ns [2.25333, 2.27027] (0.95 confidence)</td>
<td>1.08753</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>2.40125 ns [2.39756, 2.40632] (0.95 confidence)</td>
<td>2.41745 ns [2.41226, 2.42313] (0.95 confidence)</td>
<td>0.993297</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>2.39094 ns [2.3881, 2.39436] (0.95 confidence)</td>
<td>2.40852 ns [2.40066, 2.41875] (0.95 confidence)</td>
<td>0.992701</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>1.85564 ns [1.48878, 2.13753] (0.95 confidence)</td>
<td>1.48637 ns [1.48532, 1.48722] (0.95 confidence)</td>
<td>1.24844</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>2.14034 ns [2.13313, 2.14397] (0.95 confidence)</td>
<td>2.15475 ns [2.15266, 2.15709] (0.95 confidence)</td>
<td>0.99331</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;norm2&#95;sq</td>
<td>0.458391 ns [0.45785, 0.459244] (0.95 confidence)</td>
<td>0.458063 ns [0.457538, 0.460053] (0.95 confidence)</td>
<td>1.00072</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>2.46202 ns [2.4608, 2.46274] (0.95 confidence)</td>
<td>2.47726 ns [2.47502, 2.48201] (0.95 confidence)</td>
<td>0.993846</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>4.00902 ns [4.00575, 4.01151] (0.95 confidence)</td>
<td>4.02943 ns [4.02519, 4.03134] (0.95 confidence)</td>
<td>0.994933</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;det</td>
<td>28.8842 ns [28.8496, 28.928] (0.95 confidence)</td>
<td>28.8007 ns [28.7867, 28.8278] (0.95 confidence)</td>
<td>1.0029</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>0.876327 ns [0.876099, 0.876441] (0.95 confidence)</td>
<td>0.87639 ns [0.876189, 0.876563] (0.95 confidence)</td>
<td>0.999928</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>1.0321 ns [1.0319, 1.03236] (0.95 confidence)</td>
<td>1.03237 ns [1.03195, 1.03272] (0.95 confidence)</td>
<td>0.999736</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>34.1687 ns [34.1378, 34.2274] (0.95 confidence)</td>
<td>34.2748 ns [34.2429, 34.3059] (0.95 confidence)</td>
<td>0.996906</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;dot</td>
<td>0.830836 ns [0.805988, 0.836416] (0.95 confidence)</td>
<td>0.805363 ns [0.802246, 0.821759] (0.95 confidence)</td>
<td>1.03163</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;inf&#95;norm</td>
<td>3.46521 ns [3.46403, 3.46567] (0.95 confidence)</td>
<td>3.46904 ns [3.46675, 3.47035] (0.95 confidence)</td>
<td>0.998893</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;ldlt</td>
<td>24.4768 ns [24.4496, 24.4927] (0.95 confidence)</td>
<td>24.5352 ns [24.5176, 24.5524] (0.95 confidence)</td>
<td>0.997617</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;ldlt&#95;solve</td>
<td>30.2181 ns [30.1879, 30.2566] (0.95 confidence)</td>
<td>30.1372 ns [30.1034, 30.1676] (0.95 confidence)</td>
<td>1.00269</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;lu</td>
<td>33.8342 ns [33.7764, 33.9033] (0.95 confidence)</td>
<td>33.9319 ns [33.7407, 34.0166] (0.95 confidence)</td>
<td>0.997119</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;norm2</td>
<td>2.78325 ns [2.71182, 2.91223] (0.95 confidence)</td>
<td>2.93128 ns [2.7931, 2.9627] (0.95 confidence)</td>
<td>0.9495</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>3.15807 ns [3.15461, 3.16276] (0.95 confidence)</td>
<td>3.16579 ns [3.16162, 3.16873] (0.95 confidence)</td>
<td>0.997562</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>3.15283 ns [3.14698, 3.1572] (0.95 confidence)</td>
<td>3.16157 ns [3.15722, 3.17081] (0.95 confidence)</td>
<td>0.997237</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>1.72907 ns [1.72574, 1.73082] (0.95 confidence)</td>
<td>1.73142 ns [1.72837, 1.73241] (0.95 confidence)</td>
<td>0.998644</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>3.13977 ns [3.13349, 3.14464] (0.95 confidence)</td>
<td>3.15054 ns [3.14351, 3.15396] (0.95 confidence)</td>
<td>0.996583</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;norm2&#95;sq</td>
<td>0.524456 ns [0.524038, 0.525036] (0.95 confidence)</td>
<td>0.52466 ns [0.524258, 0.525374] (0.95 confidence)</td>
<td>0.999613</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>3.97389 ns [3.95891, 3.9863] (0.95 confidence)</td>
<td>3.97394 ns [3.96309, 3.97855] (0.95 confidence)</td>
<td>0.999987</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>6.06678 ns [6.06182, 6.07257] (0.95 confidence)</td>
<td>6.06033 ns [6.05514, 6.06615] (0.95 confidence)</td>
<td>1.00106</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;det</td>
<td>14818.5 ns [14800, 14838] (0.95 confidence)</td>
<td>14333.1 ns [14276.5, 14354.8] (0.95 confidence)</td>
<td>1.03387</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>23.0369 ns [23.0066, 23.0667] (0.95 confidence)</td>
<td>23.5658 ns [23.5384, 23.595] (0.95 confidence)</td>
<td>0.977557</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>23.056 ns [23.0268, 23.0899] (0.95 confidence)</td>
<td>23.6488 ns [23.6101, 23.6886] (0.95 confidence)</td>
<td>0.974934</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>14737.1 ns [14710.6, 14760.1] (0.95 confidence)</td>
<td>14333.9 ns [14300.1, 14389.1] (0.95 confidence)</td>
<td>1.02813</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;dot</td>
<td>10.925 ns [10.9084, 10.9376] (0.95 confidence)</td>
<td>11.2052 ns [11.1759, 11.2264] (0.95 confidence)</td>
<td>0.974997</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;inf&#95;norm</td>
<td>626.498 ns [625.343, 627.433] (0.95 confidence)</td>
<td>654.575 ns [653.523, 655.916] (0.95 confidence)</td>
<td>0.957106</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;ldlt</td>
<td>16828.3 ns [16811.1, 16845.8] (0.95 confidence)</td>
<td>17305.6 ns [17274.5, 17318.6] (0.95 confidence)</td>
<td>0.972421</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;ldlt&#95;solve</td>
<td>19166.9 ns [19105.5, 19197.6] (0.95 confidence)</td>
<td>19168.9 ns [19143, 19203.2] (0.95 confidence)</td>
<td>0.999895</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;lu</td>
<td>15280.8 ns [15265.9, 15291.7] (0.95 confidence)</td>
<td>15121.2 ns [15074.6, 15221.4] (0.95 confidence)</td>
<td>1.01056</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;norm2</td>
<td>57.579 ns [57.3485, 57.7494] (0.95 confidence)</td>
<td>59.3197 ns [59.1191, 60.2646] (0.95 confidence)</td>
<td>0.970656</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>36.5496 ns [36.4448, 36.6602] (0.95 confidence)</td>
<td>39.6441 ns [39.4689, 39.7012] (0.95 confidence)</td>
<td>0.921942</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>36.4419 ns [36.2314, 36.5335] (0.95 confidence)</td>
<td>39.7963 ns [39.7571, 39.8576] (0.95 confidence)</td>
<td>0.915711</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>17.064 ns [17.0316, 17.1186] (0.95 confidence)</td>
<td>18.0769 ns [18.0526, 18.1532] (0.95 confidence)</td>
<td>0.943964</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>32.5246 ns [32.4201, 32.5852] (0.95 confidence)</td>
<td>34.1111 ns [34.0674, 34.1544] (0.95 confidence)</td>
<td>0.953489</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;norm2&#95;sq</td>
<td>10.7521 ns [10.7282, 10.7797] (0.95 confidence)</td>
<td>11.0605 ns [11.036, 11.093] (0.95 confidence)</td>
<td>0.972114</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>1017.7 ns [1015.63, 1020.12] (0.95 confidence)</td>
<td>1036.38 ns [1034.63, 1038.72] (0.95 confidence)</td>
<td>0.981981</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>2517.5 ns [2513.92, 2520.17] (0.95 confidence)</td>
<td>2631.47 ns [2624.95, 2639.85] (0.95 confidence)</td>
<td>0.956689</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;det</td>
<td>82.868 ns [82.3697, 82.9824] (0.95 confidence)</td>
<td>82.9638 ns [82.7774, 83.3394] (0.95 confidence)</td>
<td>0.998846</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>1.31758 ns [1.31648, 1.31802] (0.95 confidence)</td>
<td>1.30528 ns [1.30434, 1.30718] (0.95 confidence)</td>
<td>1.00943</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;det&#95;from&#95;ldlt&#95;balanced&#95;range</td>
<td>8.66089 ns [8.65948, 8.66331] (0.95 confidence)</td>
<td>8.89311 ns [8.85191, 8.90858] (0.95 confidence)</td>
<td>0.973887</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>1.46562 ns [1.46328, 1.46709] (0.95 confidence)</td>
<td>1.46767 ns [1.46608, 1.4688] (0.95 confidence)</td>
<td>0.998603</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;det&#95;from&#95;lu&#95;balanced&#95;range</td>
<td>8.9098 ns [8.90616, 8.9245] (0.95 confidence)</td>
<td>11.2185 ns [11.2034, 11.2393] (0.95 confidence)</td>
<td>0.794205</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>92.7072 ns [92.3066, 93.4104] (0.95 confidence)</td>
<td>92.1068 ns [91.9193, 92.4603] (0.95 confidence)</td>
<td>1.00652</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;dot</td>
<td>0.931383 ns [0.925274, 0.934744] (0.95 confidence)</td>
<td>0.935043 ns [0.932987, 0.938274] (0.95 confidence)</td>
<td>0.996086</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;inf&#95;norm</td>
<td>8.44404 ns [8.42618, 8.45928] (0.95 confidence)</td>
<td>8.51908 ns [8.50652, 8.52786] (0.95 confidence)</td>
<td>0.991192</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;ldlt</td>
<td>90.9539 ns [89.859, 91.3991] (0.95 confidence)</td>
<td>91.1351 ns [89.7375, 91.5908] (0.95 confidence)</td>
<td>0.998012</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;ldlt&#95;ill&#95;conditioned</td>
<td>90.7713 ns [89.7363, 92.016] (0.95 confidence)</td>
<td>92.1541 ns [91.0903, 92.7969] (0.95 confidence)</td>
<td>0.984995</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;ldlt&#95;solve</td>
<td>103.56 ns [102.259, 104.077] (0.95 confidence)</td>
<td>103.242 ns [102.352, 104.137] (0.95 confidence)</td>
<td>1.00308</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;lu</td>
<td>90.2517 ns [90.0624, 90.5806] (0.95 confidence)</td>
<td>90.5462 ns [90.3236, 91.0277] (0.95 confidence)</td>
<td>0.996747</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;lu&#95;ill&#95;conditioned</td>
<td>89.7649 ns [89.675, 89.8994] (0.95 confidence)</td>
<td>92.0866 ns [91.7234, 92.4472] (0.95 confidence)</td>
<td>0.974787</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;lu&#95;pivoting</td>
<td>93.8967 ns [93.5466, 94.6072] (0.95 confidence)</td>
<td>95.8871 ns [95.4773, 96.2629] (0.95 confidence)</td>
<td>0.979242</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;norm2</td>
<td>4.57895 ns [4.55036, 4.60343] (0.95 confidence)</td>
<td>4.66988 ns [4.66131, 4.68635] (0.95 confidence)</td>
<td>0.980527</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>4.82048 ns [4.81419, 4.82777] (0.95 confidence)</td>
<td>4.86745 ns [4.8589, 4.88095] (0.95 confidence)</td>
<td>0.990351</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>4.82121 ns [4.81107, 4.83147] (0.95 confidence)</td>
<td>4.86929 ns [4.85811, 4.88408] (0.95 confidence)</td>
<td>0.990126</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>2.76971 ns [2.49563, 3.12693] (0.95 confidence)</td>
<td>2.55371 ns [2.5084, 2.739] (0.95 confidence)</td>
<td>1.08458</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>4.21264 ns [4.2035, 4.2216] (0.95 confidence)</td>
<td>4.26873 ns [4.26332, 4.27917] (0.95 confidence)</td>
<td>0.98686</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;norm2&#95;sq</td>
<td>0.663461 ns [0.66132, 0.664177] (0.95 confidence)</td>
<td>0.656503 ns [0.654255, 0.660009] (0.95 confidence)</td>
<td>1.0106</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>8.11136 ns [8.0935, 8.11507] (0.95 confidence)</td>
<td>8.11594 ns [8.09985, 8.12444] (0.95 confidence)</td>
<td>0.999436</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>13.6736 ns [13.6479, 13.727] (0.95 confidence)</td>
<td>13.5068 ns [13.4918, 13.525] (0.95 confidence)</td>
<td>1.01235</td>
</tr>
<tr>
<td>gram/1x2/dependent/la&#95;stack</td>
<td>0.684229 ns [0.683789, 0.684859] (0.95 confidence)</td>
<td>0.69171 ns [0.691382, 0.692041] (0.95 confidence)</td>
<td>0.989186</td>
</tr>
<tr>
<td>gram/1x2/mixed&#95;scale/la&#95;stack</td>
<td>0.684608 ns [0.684184, 0.685] (0.95 confidence)</td>
<td>0.692254 ns [0.691899, 0.692956] (0.95 confidence)</td>
<td>0.988956</td>
</tr>
<tr>
<td>gram/1x2/near&#95;dependent/la&#95;stack</td>
<td>0.684604 ns [0.684214, 0.684993] (0.95 confidence)</td>
<td>0.693132 ns [0.692208, 0.695322] (0.95 confidence)</td>
<td>0.987697</td>
</tr>
<tr>
<td>gram/1x2/orthogonal/la&#95;stack</td>
<td>0.684063 ns [0.683726, 0.684896] (0.95 confidence)</td>
<td>0.69207 ns [0.691695, 0.692527] (0.95 confidence)</td>
<td>0.988431</td>
</tr>
<tr>
<td>gram/2x2/dependent/la&#95;stack</td>
<td>1.31232 ns [1.3003, 1.32008] (0.95 confidence)</td>
<td>1.32701 ns [1.31472, 1.33478] (0.95 confidence)</td>
<td>0.988929</td>
</tr>
<tr>
<td>gram/2x2/mixed&#95;scale/la&#95;stack</td>
<td>1.32591 ns [1.32103, 1.33303] (0.95 confidence)</td>
<td>1.32759 ns [1.31542, 1.33596] (0.95 confidence)</td>
<td>0.998737</td>
</tr>
<tr>
<td>gram/2x2/near&#95;dependent/la&#95;stack</td>
<td>1.31031 ns [1.30009, 1.31763] (0.95 confidence)</td>
<td>1.32621 ns [1.31523, 1.33688] (0.95 confidence)</td>
<td>0.988015</td>
</tr>
<tr>
<td>gram/2x2/orthogonal/la&#95;stack</td>
<td>1.31369 ns [1.2996, 1.32253] (0.95 confidence)</td>
<td>1.32662 ns [1.31556, 1.33756] (0.95 confidence)</td>
<td>0.99025</td>
</tr>
<tr>
<td>gram/2x3/dependent/la&#95;stack</td>
<td>1.26487 ns [1.26392, 1.26855] (0.95 confidence)</td>
<td>1.29944 ns [1.29533, 1.30481] (0.95 confidence)</td>
<td>0.973397</td>
</tr>
<tr>
<td>gram/2x3/mixed&#95;scale/la&#95;stack</td>
<td>1.26599 ns [1.26435, 1.26981] (0.95 confidence)</td>
<td>1.31015 ns [1.30397, 1.31308] (0.95 confidence)</td>
<td>0.966293</td>
</tr>
<tr>
<td>gram/2x3/near&#95;dependent/la&#95;stack</td>
<td>1.26591 ns [1.26355, 1.26894] (0.95 confidence)</td>
<td>1.29008 ns [1.28688, 1.29431] (0.95 confidence)</td>
<td>0.981268</td>
</tr>
<tr>
<td>gram/2x3/orthogonal/la&#95;stack</td>
<td>1.266 ns [1.26391, 1.26984] (0.95 confidence)</td>
<td>1.28448 ns [1.2829, 1.28872] (0.95 confidence)</td>
<td>0.985613</td>
</tr>
<tr>
<td>gram/3x3/dependent/la&#95;stack</td>
<td>3.00337 ns [3.00057, 3.00792] (0.95 confidence)</td>
<td>3.02629 ns [3.02338, 3.02975] (0.95 confidence)</td>
<td>0.992425</td>
</tr>
<tr>
<td>gram/3x3/mixed&#95;scale/la&#95;stack</td>
<td>3.00576 ns [3.00184, 3.0126] (0.95 confidence)</td>
<td>3.03267 ns [3.02743, 3.03642] (0.95 confidence)</td>
<td>0.991126</td>
</tr>
<tr>
<td>gram/3x3/near&#95;dependent/la&#95;stack</td>
<td>3.00198 ns [2.99852, 3.0064] (0.95 confidence)</td>
<td>3.02695 ns [3.02281, 3.03025] (0.95 confidence)</td>
<td>0.991749</td>
</tr>
<tr>
<td>gram/3x3/orthogonal/la&#95;stack</td>
<td>3.00342 ns [2.99982, 3.00712] (0.95 confidence)</td>
<td>3.02956 ns [3.02647, 3.0349] (0.95 confidence)</td>
<td>0.991373</td>
</tr>
<tr>
<td>gram/3x4/dependent/la&#95;stack</td>
<td>3.05166 ns [3.04861, 3.05323] (0.95 confidence)</td>
<td>4.11139 ns [4.03087, 4.22397] (0.95 confidence)</td>
<td>0.742244</td>
</tr>
<tr>
<td>gram/3x4/mixed&#95;scale/la&#95;stack</td>
<td>3.05183 ns [3.04777, 3.05404] (0.95 confidence)</td>
<td>3.12804 ns [3.11861, 3.13323] (0.95 confidence)</td>
<td>0.975635</td>
</tr>
<tr>
<td>gram/3x4/near&#95;dependent/la&#95;stack</td>
<td>3.0501 ns [3.04691, 3.05319] (0.95 confidence)</td>
<td>3.65582 ns [3.6472, 3.66783] (0.95 confidence)</td>
<td>0.834314</td>
</tr>
<tr>
<td>gram/3x4/orthogonal/la&#95;stack</td>
<td>3.04911 ns [3.04597, 3.05242] (0.95 confidence)</td>
<td>3.08841 ns [3.08148, 3.09308] (0.95 confidence)</td>
<td>0.987276</td>
</tr>
<tr>
<td>gram/4x4/dependent/la&#95;stack</td>
<td>7.50056 ns [7.40922, 7.77175] (0.95 confidence)</td>
<td>7.8628 ns [7.52518, 7.90086] (0.95 confidence)</td>
<td>0.953929</td>
</tr>
<tr>
<td>gram/4x4/mixed&#95;scale/la&#95;stack</td>
<td>7.641 ns [7.43389, 7.7669] (0.95 confidence)</td>
<td>8.56386 ns [8.37038, 8.89095] (0.95 confidence)</td>
<td>0.892237</td>
</tr>
<tr>
<td>gram/4x4/near&#95;dependent/la&#95;stack</td>
<td>7.53183 ns [7.40221, 7.73521] (0.95 confidence)</td>
<td>7.86048 ns [7.71498, 7.88698] (0.95 confidence)</td>
<td>0.958189</td>
</tr>
<tr>
<td>gram/4x4/orthogonal/la&#95;stack</td>
<td>7.50282 ns [7.42251, 7.76078] (0.95 confidence)</td>
<td>7.86457 ns [7.63127, 7.8992] (0.95 confidence)</td>
<td>0.954002</td>
</tr>
<tr>
<td>gram/4x5/dependent/la&#95;stack</td>
<td>9.49516 ns [9.48113, 9.50528] (0.95 confidence)</td>
<td>11.3335 ns [10.9379, 11.8234] (0.95 confidence)</td>
<td>0.837795</td>
</tr>
<tr>
<td>gram/4x5/mixed&#95;scale/la&#95;stack</td>
<td>9.47728 ns [9.45932, 9.50018] (0.95 confidence)</td>
<td>12.2973 ns [11.7177, 12.9628] (0.95 confidence)</td>
<td>0.770683</td>
</tr>
<tr>
<td>gram/4x5/near&#95;dependent/la&#95;stack</td>
<td>9.50009 ns [9.48814, 9.50614] (0.95 confidence)</td>
<td>12.2642 ns [11.8829, 12.7166] (0.95 confidence)</td>
<td>0.774619</td>
</tr>
<tr>
<td>gram/4x5/orthogonal/la&#95;stack</td>
<td>9.47384 ns [9.4517, 9.48941] (0.95 confidence)</td>
<td>9.77945 ns [9.75866, 9.82576] (0.95 confidence)</td>
<td>0.96875</td>
</tr>
<tr>
<td>gram/5x5/dependent/la&#95;stack</td>
<td>18.4158 ns [18.3699, 18.4389] (0.95 confidence)</td>
<td>18.5274 ns [18.4753, 18.6162] (0.95 confidence)</td>
<td>0.993977</td>
</tr>
<tr>
<td>gram/5x5/mixed&#95;scale/la&#95;stack</td>
<td>18.5882 ns [18.5251, 18.6493] (0.95 confidence)</td>
<td>18.5758 ns [18.4259, 18.7005] (0.95 confidence)</td>
<td>1.00067</td>
</tr>
<tr>
<td>gram/5x5/near&#95;dependent/la&#95;stack</td>
<td>18.5812 ns [18.5012, 18.6411] (0.95 confidence)</td>
<td>18.6063 ns [18.5495, 18.6781] (0.95 confidence)</td>
<td>0.998652</td>
</tr>
<tr>
<td>gram/5x5/orthogonal/la&#95;stack</td>
<td>18.6111 ns [18.5693, 18.6442] (0.95 confidence)</td>
<td>18.8591 ns [18.8089, 18.8984] (0.95 confidence)</td>
<td>0.986847</td>
</tr>
<tr>
<td>gram/5x6/dependent/la&#95;stack</td>
<td>19.2937 ns [19.2542, 19.3723] (0.95 confidence)</td>
<td>28.8584 ns [27.883, 29.6033] (0.95 confidence)</td>
<td>0.668564</td>
</tr>
<tr>
<td>gram/5x6/mixed&#95;scale/la&#95;stack</td>
<td>19.273 ns [19.2338, 19.3226] (0.95 confidence)</td>
<td>19.4757 ns [19.3522, 19.5855] (0.95 confidence)</td>
<td>0.98959</td>
</tr>
<tr>
<td>gram/5x6/near&#95;dependent/la&#95;stack</td>
<td>19.198 ns [19.1354, 19.2415] (0.95 confidence)</td>
<td>25.4334 ns [25.2068, 25.9288] (0.95 confidence)</td>
<td>0.754835</td>
</tr>
<tr>
<td>gram/5x6/orthogonal/la&#95;stack</td>
<td>19.2193 ns [19.176, 19.2667] (0.95 confidence)</td>
<td>21.4348 ns [21.2767, 21.544] (0.95 confidence)</td>
<td>0.896641</td>
</tr>
<tr>
<td>gram/6x6/dependent/la&#95;stack</td>
<td>28.9829 ns [28.7761, 29.16] (0.95 confidence)</td>
<td>39.0206 ns [37.9749, 40.0332] (0.95 confidence)</td>
<td>0.742759</td>
</tr>
<tr>
<td>gram/6x6/mixed&#95;scale/la&#95;stack</td>
<td>28.7987 ns [28.7066, 28.9625] (0.95 confidence)</td>
<td>40.0538 ns [39.0992, 41.0343] (0.95 confidence)</td>
<td>0.719</td>
</tr>
<tr>
<td>gram/6x6/near&#95;dependent/la&#95;stack</td>
<td>28.8337 ns [28.7914, 29.3202] (0.95 confidence)</td>
<td>39.5323 ns [38.7857, 40.1376] (0.95 confidence)</td>
<td>0.729371</td>
</tr>
<tr>
<td>gram/6x6/orthogonal/la&#95;stack</td>
<td>29.1075 ns [28.8232, 29.7981] (0.95 confidence)</td>
<td>37.6856 ns [35.8543, 38.6281] (0.95 confidence)</td>
<td>0.772376</td>
</tr>
<tr>
<td>gram/6x7/dependent/la&#95;stack</td>
<td>31.9708 ns [31.8342, 32.5803] (0.95 confidence)</td>
<td>31.8465 ns [31.7844, 31.9219] (0.95 confidence)</td>
<td>1.0039</td>
</tr>
<tr>
<td>gram/6x7/mixed&#95;scale/la&#95;stack</td>
<td>32.0129 ns [31.8568, 32.6535] (0.95 confidence)</td>
<td>32.3463 ns [31.9654, 32.4993] (0.95 confidence)</td>
<td>0.989695</td>
</tr>
<tr>
<td>gram/6x7/near&#95;dependent/la&#95;stack</td>
<td>31.9829 ns [31.9146, 32.6476] (0.95 confidence)</td>
<td>32.0185 ns [31.8023, 32.5294] (0.95 confidence)</td>
<td>0.998887</td>
</tr>
<tr>
<td>gram/6x7/orthogonal/la&#95;stack</td>
<td>31.8877 ns [31.5693, 32.583] (0.95 confidence)</td>
<td>32.4421 ns [31.5858, 32.7629] (0.95 confidence)</td>
<td>0.982912</td>
</tr>
<tr>
<td>gram/7x7/dependent/la&#95;stack</td>
<td>38.9444 ns [38.821, 39.194] (0.95 confidence)</td>
<td>71.1962 ns [69.2428, 72.1129] (0.95 confidence)</td>
<td>0.547002</td>
</tr>
<tr>
<td>gram/7x7/mixed&#95;scale/la&#95;stack</td>
<td>39.3985 ns [39.2382, 39.5797] (0.95 confidence)</td>
<td>45.1211 ns [44.8914, 45.399] (0.95 confidence)</td>
<td>0.873174</td>
</tr>
<tr>
<td>gram/7x7/near&#95;dependent/la&#95;stack</td>
<td>39.1492 ns [39.0734, 39.2164] (0.95 confidence)</td>
<td>42.3281 ns [41.8826, 42.6294] (0.95 confidence)</td>
<td>0.924898</td>
</tr>
<tr>
<td>gram/7x7/orthogonal/la&#95;stack</td>
<td>39.1878 ns [39.0817, 39.3313] (0.95 confidence)</td>
<td>55.5328 ns [53.6728, 59.3069] (0.95 confidence)</td>
<td>0.705669</td>
</tr>
<tr>
<td>gram/7x8/dependent/la&#95;stack</td>
<td>40.1163 ns [40.0207, 40.2913] (0.95 confidence)</td>
<td>40.5952 ns [40.4791, 40.7718] (0.95 confidence)</td>
<td>0.988201</td>
</tr>
<tr>
<td>gram/7x8/mixed&#95;scale/la&#95;stack</td>
<td>40.322 ns [40.1409, 40.6188] (0.95 confidence)</td>
<td>41.7811 ns [41.0962, 43.056] (0.95 confidence)</td>
<td>0.965078</td>
</tr>
<tr>
<td>gram/7x8/near&#95;dependent/la&#95;stack</td>
<td>40.0791 ns [39.9618, 40.3296] (0.95 confidence)</td>
<td>40.9713 ns [40.892, 41.2455] (0.95 confidence)</td>
<td>0.978224</td>
</tr>
<tr>
<td>gram/7x8/orthogonal/la&#95;stack</td>
<td>40.0611 ns [39.8425, 40.7587] (0.95 confidence)</td>
<td>41.2898 ns [41.1586, 41.4859] (0.95 confidence)</td>
<td>0.970243</td>
</tr>
<tr>
<td>gram/8x8/dependent/la&#95;stack</td>
<td>51.7562 ns [51.713, 51.8508] (0.95 confidence)</td>
<td>52.5845 ns [51.9911, 53.1356] (0.95 confidence)</td>
<td>0.984248</td>
</tr>
<tr>
<td>gram/8x8/mixed&#95;scale/la&#95;stack</td>
<td>51.8232 ns [51.7812, 51.851] (0.95 confidence)</td>
<td>51.7701 ns [51.7297, 51.8398] (0.95 confidence)</td>
<td>1.00103</td>
</tr>
<tr>
<td>gram/8x8/near&#95;dependent/la&#95;stack</td>
<td>51.7969 ns [51.7593, 51.8397] (0.95 confidence)</td>
<td>51.7119 ns [51.6549, 51.8337] (0.95 confidence)</td>
<td>1.00165</td>
</tr>
<tr>
<td>gram/8x8/orthogonal/la&#95;stack</td>
<td>51.8968 ns [51.8231, 51.9921] (0.95 confidence)</td>
<td>51.722 ns [51.6143, 51.7972] (0.95 confidence)</td>
<td>1.00338</td>
</tr>
<tr>
<td>interval&#95;d3/dense&#95;inconclusive</td>
<td>41.635 ns [41.4991, 41.7007] (0.95 confidence)</td>
<td>43.1819 ns [42.943, 43.3697] (0.95 confidence)</td>
<td>0.964177</td>
</tr>
<tr>
<td>interval&#95;d3/dense&#95;sign</td>
<td>37.6097 ns [37.1925, 38.2345] (0.95 confidence)</td>
<td>39.4633 ns [38.9232, 39.8069] (0.95 confidence)</td>
<td>0.953031</td>
</tr>
<tr>
<td>interval&#95;d3/lifted&#95;assembly</td>
<td>8.35286 ns [8.34584, 8.3618] (0.95 confidence)</td>
<td>8.40164 ns [8.39496, 8.41053] (0.95 confidence)</td>
<td>0.994194</td>
</tr>
<tr>
<td>interval&#95;d3/lifted&#95;assembly&#95;and&#95;sign</td>
<td>36.9994 ns [36.8524, 37.256] (0.95 confidence)</td>
<td>37.2917 ns [37.0857, 37.5728] (0.95 confidence)</td>
<td>0.992162</td>
</tr>
<tr>
<td>interval&#95;d3/sparse&#95;sign</td>
<td>31.7756 ns [31.633, 31.845] (0.95 confidence)</td>
<td>30.9618 ns [30.6802, 31.2864] (0.95 confidence)</td>
<td>1.02628</td>
</tr>
<tr>
<td>interval&#95;d4/dense&#95;inconclusive</td>
<td>99.6798 ns [99.3223, 100.163] (0.95 confidence)</td>
<td>105.26 ns [104.365, 106.521] (0.95 confidence)</td>
<td>0.946985</td>
</tr>
<tr>
<td>interval&#95;d4/dense&#95;sign</td>
<td>100.402 ns [99.8368, 101.762] (0.95 confidence)</td>
<td>109.682 ns [108.499, 110.743] (0.95 confidence)</td>
<td>0.915387</td>
</tr>
<tr>
<td>interval&#95;d4/lifted&#95;assembly</td>
<td>13.3029 ns [13.289, 13.3145] (0.95 confidence)</td>
<td>13.4764 ns [13.4487, 13.5219] (0.95 confidence)</td>
<td>0.987122</td>
</tr>
<tr>
<td>interval&#95;d4/lifted&#95;assembly&#95;and&#95;sign</td>
<td>82.3838 ns [81.7697, 83.2433] (0.95 confidence)</td>
<td>84.907 ns [84.3925, 85.583] (0.95 confidence)</td>
<td>0.970283</td>
</tr>
<tr>
<td>interval&#95;d4/sparse&#95;sign</td>
<td>68.0468 ns [67.4648, 68.5741] (0.95 confidence)</td>
<td>70.185 ns [69.4682, 70.9479] (0.95 confidence)</td>
<td>0.969534</td>
</tr>
<tr>
<td>interval&#95;d5/dense&#95;inconclusive</td>
<td>237.338 ns [236.518, 237.94] (0.95 confidence)</td>
<td>237.555 ns [236.452, 239.15] (0.95 confidence)</td>
<td>0.999086</td>
</tr>
<tr>
<td>interval&#95;d5/dense&#95;sign</td>
<td>221.918 ns [221.023, 223.676] (0.95 confidence)</td>
<td>230.419 ns [228.457, 232.065] (0.95 confidence)</td>
<td>0.963107</td>
</tr>
<tr>
<td>interval&#95;d5/lifted&#95;assembly</td>
<td>21.7081 ns [21.4916, 21.8416] (0.95 confidence)</td>
<td>22.0266 ns [21.8003, 22.0856] (0.95 confidence)</td>
<td>0.985539</td>
</tr>
<tr>
<td>interval&#95;d5/lifted&#95;assembly&#95;and&#95;sign</td>
<td>186.655 ns [185.669, 187.722] (0.95 confidence)</td>
<td>197.433 ns [195.795, 198.409] (0.95 confidence)</td>
<td>0.945409</td>
</tr>
<tr>
<td>interval&#95;d5/sparse&#95;sign</td>
<td>154.057 ns [153.46, 154.443] (0.95 confidence)</td>
<td>166.74 ns [165.169, 168.607] (0.95 confidence)</td>
<td>0.923938</td>
</tr>
<tr>
<td>interval&#95;d6/dense&#95;inconclusive</td>
<td>535.577 ns [533.532, 537.152] (0.95 confidence)</td>
<td>537.631 ns [535.231, 542.519] (0.95 confidence)</td>
<td>0.996179</td>
</tr>
<tr>
<td>interval&#95;d6/dense&#95;sign</td>
<td>524.891 ns [522.984, 527.096] (0.95 confidence)</td>
<td>527.354 ns [525.11, 528.906] (0.95 confidence)</td>
<td>0.995328</td>
</tr>
<tr>
<td>interval&#95;d6/lifted&#95;assembly</td>
<td>31.0141 ns [30.968, 31.051] (0.95 confidence)</td>
<td>31.2335 ns [31.1774, 31.2579] (0.95 confidence)</td>
<td>0.992974</td>
</tr>
<tr>
<td>interval&#95;d6/lifted&#95;assembly&#95;and&#95;sign</td>
<td>385.455 ns [384.042, 387.327] (0.95 confidence)</td>
<td>396.889 ns [395.064, 397.882] (0.95 confidence)</td>
<td>0.971191</td>
</tr>
<tr>
<td>interval&#95;d6/sparse&#95;sign</td>
<td>365.124 ns [363.519, 366.762] (0.95 confidence)</td>
<td>374.913 ns [372.932, 376.281] (0.95 confidence)</td>
<td>0.973892</td>
</tr>
<tr>
<td>interval&#95;d7/dense&#95;inconclusive</td>
<td>1176.37 ns [1167.94, 1180.32] (0.95 confidence)</td>
<td>1155.71 ns [1149.25, 1164.46] (0.95 confidence)</td>
<td>1.01788</td>
</tr>
<tr>
<td>interval&#95;d7/dense&#95;sign</td>
<td>1141.22 ns [1132.74, 1145.3] (0.95 confidence)</td>
<td>1115.23 ns [1110.25, 1123.59] (0.95 confidence)</td>
<td>1.0233</td>
</tr>
<tr>
<td>interval&#95;d7/lifted&#95;assembly</td>
<td>37.5019 ns [37.4174, 37.5314] (0.95 confidence)</td>
<td>38.0318 ns [37.9599, 38.09] (0.95 confidence)</td>
<td>0.986067</td>
</tr>
<tr>
<td>interval&#95;d7/lifted&#95;assembly&#95;and&#95;sign</td>
<td>777.718 ns [776.184, 779.122] (0.95 confidence)</td>
<td>806.018 ns [799.349, 809.778] (0.95 confidence)</td>
<td>0.964889</td>
</tr>
<tr>
<td>interval&#95;d7/sparse&#95;sign</td>
<td>749.928 ns [748.272, 751.237] (0.95 confidence)</td>
<td>768.668 ns [766.364, 770.058] (0.95 confidence)</td>
<td>0.97562</td>
</tr>
<tr>
<td>interval&#95;det&#95;sign/d4&#95;conclusive&#95;lifted</td>
<td>68.9407 ns [68.4857, 69.6622] (0.95 confidence)</td>
<td>67.0797 ns [66.8191, 67.3668] (0.95 confidence)</td>
<td>1.02774</td>
</tr>
<tr>
<td>interval&#95;det&#95;sign/d4&#95;inconclusive&#95;lifted</td>
<td>68.3062 ns [67.5004, 68.9004] (0.95 confidence)</td>
<td>65.7997 ns [65.3996, 66.1733] (0.95 confidence)</td>
<td>1.03809</td>
</tr>
<tr>
<td>interval&#95;det&#95;sign/d7&#95;conclusive&#95;lifted</td>
<td>747.271 ns [745.77, 749.17] (0.95 confidence)</td>
<td>754.729 ns [753.013, 756.398] (0.95 confidence)</td>
<td>0.990119</td>
</tr>
<tr>
<td>interval&#95;scalar/cancellation/add</td>
<td>2.13084 ns [2.12784, 2.13519] (0.95 confidence)</td>
<td>2.12391 ns [2.1218, 2.12577] (0.95 confidence)</td>
<td>1.00326</td>
</tr>
<tr>
<td>interval&#95;scalar/cancellation/multiply</td>
<td>1.8833 ns [1.87486, 1.89242] (0.95 confidence)</td>
<td>1.88661 ns [1.88013, 1.89052] (0.95 confidence)</td>
<td>0.998247</td>
</tr>
<tr>
<td>interval&#95;scalar/cancellation/square</td>
<td>1.87321 ns [1.87049, 1.87501] (0.95 confidence)</td>
<td>1.86908 ns [1.86538, 1.87393] (0.95 confidence)</td>
<td>1.00221</td>
</tr>
<tr>
<td>interval&#95;scalar/exact/add</td>
<td>2.12996 ns [2.12478, 2.13244] (0.95 confidence)</td>
<td>2.09857 ns [2.09221, 2.101] (0.95 confidence)</td>
<td>1.01495</td>
</tr>
<tr>
<td>interval&#95;scalar/exact/multiply</td>
<td>2.10838 ns [2.08958, 2.11662] (0.95 confidence)</td>
<td>2.09027 ns [2.07247, 2.10658] (0.95 confidence)</td>
<td>1.00867</td>
</tr>
<tr>
<td>interval&#95;scalar/exact/square</td>
<td>1.87429 ns [1.87119, 1.87653] (0.95 confidence)</td>
<td>1.84931 ns [1.84713, 1.85147] (0.95 confidence)</td>
<td>1.01351</td>
</tr>
<tr>
<td>interval&#95;scalar/point/add</td>
<td>2.1766 ns [2.17517, 2.17802] (0.95 confidence)</td>
<td>2.16712 ns [2.16248, 2.17218] (0.95 confidence)</td>
<td>1.00438</td>
</tr>
<tr>
<td>interval&#95;scalar/point/multiply</td>
<td>2.58476 ns [2.57841, 2.60239] (0.95 confidence)</td>
<td>2.57734 ns [2.56809, 2.58866] (0.95 confidence)</td>
<td>1.00288</td>
</tr>
<tr>
<td>interval&#95;scalar/point/square</td>
<td>2.12967 ns [2.12655, 2.1323] (0.95 confidence)</td>
<td>2.08281 ns [2.08046, 2.08873] (0.95 confidence)</td>
<td>1.0225</td>
</tr>
<tr>
<td>interval&#95;scalar/subnormal/add</td>
<td>2.40486 ns [2.40336, 2.40569] (0.95 confidence)</td>
<td>2.40141 ns [2.39922, 2.40252] (0.95 confidence)</td>
<td>1.00144</td>
</tr>
<tr>
<td>interval&#95;scalar/subnormal/multiply</td>
<td>2.82492 ns [2.81988, 2.83888] (0.95 confidence)</td>
<td>2.95377 ns [2.88166, 2.97064] (0.95 confidence)</td>
<td>0.956376</td>
</tr>
<tr>
<td>interval&#95;scalar/subnormal/square</td>
<td>2.37632 ns [2.37518, 2.37748] (0.95 confidence)</td>
<td>2.37243 ns [2.37069, 2.37324] (0.95 confidence)</td>
<td>1.00164</td>
</tr>
<tr>
<td>interval&#95;scalar/zero/add</td>
<td>0.939894 ns [0.939378, 0.941275] (0.95 confidence)</td>
<td>0.933311 ns [0.932303, 0.935459] (0.95 confidence)</td>
<td>1.00705</td>
</tr>
<tr>
<td>interval&#95;scalar/zero/multiply</td>
<td>1.76726 ns [1.76514, 1.77002] (0.95 confidence)</td>
<td>1.75643 ns [1.74882, 1.75869] (0.95 confidence)</td>
<td>1.00617</td>
</tr>
<tr>
<td>interval&#95;scalar/zero/square</td>
<td>0.882122 ns [0.881187, 0.882494] (0.95 confidence)</td>
<td>0.880718 ns [0.879529, 0.88131] (0.95 confidence)</td>
<td>1.00159</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;batch/origin/construction</td>
<td>2.07152 ns [2.0689, 2.07494] (0.95 confidence)</td>
<td>2.08439 ns [2.08005, 2.08849] (0.95 confidence)</td>
<td>0.993826</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;batch/origin/prepared</td>
<td>13.8219 ns [13.8198, 13.8242] (0.95 confidence)</td>
<td>13.8345 ns [13.8255, 13.8523] (0.95 confidence)</td>
<td>0.999083</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;batch/origin/with&#95;construction</td>
<td>15.3355 ns [15.3308, 15.3422] (0.95 confidence)</td>
<td>15.3371 ns [15.3293, 15.3434] (0.95 confidence)</td>
<td>0.999898</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;batch/translated/construction</td>
<td>2.07389 ns [2.07188, 2.07632] (0.95 confidence)</td>
<td>2.07594 ns [2.07265, 2.07839] (0.95 confidence)</td>
<td>0.999011</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;batch/translated/prepared</td>
<td>40.9044 ns [40.8444, 40.9544] (0.95 confidence)</td>
<td>41.1115 ns [41.0046, 41.1934] (0.95 confidence)</td>
<td>0.994962</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;batch/translated/with&#95;construction</td>
<td>42.6453 ns [42.5735, 42.7954] (0.95 confidence)</td>
<td>43.0572 ns [42.901, 43.1645] (0.95 confidence)</td>
<td>0.990434</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/cancellation/difference</td>
<td>10.7001 ns [10.6985, 10.7051] (0.95 confidence)</td>
<td>10.7021 ns [10.6985, 10.719] (0.95 confidence)</td>
<td>0.999812</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/cancellation/difference&#95;endpoints</td>
<td>11.6173 ns [11.5742, 11.6652] (0.95 confidence)</td>
<td>11.5346 ns [11.5169, 11.5871] (0.95 confidence)</td>
<td>1.00717</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/cancellation/dot</td>
<td>5.91941 ns [5.91527, 5.92305] (0.95 confidence)</td>
<td>5.96751 ns [5.93936, 5.99333] (0.95 confidence)</td>
<td>0.991941</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/cancellation/dot&#95;endpoints</td>
<td>7.04489 ns [7.02518, 7.05544] (0.95 confidence)</td>
<td>7.15772 ns [7.11534, 7.21953] (0.95 confidence)</td>
<td>0.984236</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense/difference</td>
<td>10.3354 ns [10.3298, 10.3407] (0.95 confidence)</td>
<td>10.3528 ns [10.3433, 10.3774] (0.95 confidence)</td>
<td>0.998326</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense/difference&#95;endpoints</td>
<td>12.6227 ns [12.5986, 12.6419] (0.95 confidence)</td>
<td>12.5751 ns [12.5615, 12.5987] (0.95 confidence)</td>
<td>1.00378</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense/dot</td>
<td>5.51592 ns [5.50429, 5.52738] (0.95 confidence)</td>
<td>5.56122 ns [5.55192, 5.57774] (0.95 confidence)</td>
<td>0.991854</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense/dot&#95;endpoints</td>
<td>6.5011 ns [6.48707, 6.50977] (0.95 confidence)</td>
<td>6.60232 ns [6.57301, 6.63038] (0.95 confidence)</td>
<td>0.98467</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense&#95;inexact/difference</td>
<td>10.3515 ns [10.3422, 10.3565] (0.95 confidence)</td>
<td>10.3393 ns [10.3244, 10.3464] (0.95 confidence)</td>
<td>1.00117</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense&#95;inexact/difference&#95;endpoints</td>
<td>11.6002 ns [11.589, 11.6093] (0.95 confidence)</td>
<td>11.6542 ns [11.6105, 11.8614] (0.95 confidence)</td>
<td>0.995361</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense&#95;inexact/dot</td>
<td>5.53229 ns [5.51129, 5.55794] (0.95 confidence)</td>
<td>5.59089 ns [5.58485, 5.59688] (0.95 confidence)</td>
<td>0.989518</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense&#95;inexact/dot&#95;endpoints</td>
<td>6.70165 ns [6.65882, 6.74659] (0.95 confidence)</td>
<td>6.57862 ns [6.55817, 6.59453] (0.95 confidence)</td>
<td>1.0187</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/overflow/dot</td>
<td>1.06914 ns [1.06842, 1.07004] (0.95 confidence)</td>
<td>1.06914 ns [1.06753, 1.07072] (0.95 confidence)</td>
<td>1</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/sparse/difference</td>
<td>4.73766 ns [4.73385, 4.7417] (0.95 confidence)</td>
<td>4.73061 ns [4.72871, 4.73251] (0.95 confidence)</td>
<td>1.00149</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/sparse/difference&#95;endpoints</td>
<td>5.60267 ns [5.58722, 5.61992] (0.95 confidence)</td>
<td>5.61838 ns [5.60347, 5.63137] (0.95 confidence)</td>
<td>0.997202</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/sparse/dot</td>
<td>2.22647 ns [2.22572, 2.22683] (0.95 confidence)</td>
<td>2.2295 ns [2.22848, 2.23182] (0.95 confidence)</td>
<td>0.998643</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/sparse/dot&#95;endpoints</td>
<td>2.22681 ns [2.22645, 2.2271] (0.95 confidence)</td>
<td>2.22813 ns [2.2278, 2.22903] (0.95 confidence)</td>
<td>0.999408</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/underflow/difference</td>
<td>4.45845 ns [4.45734, 4.45931] (0.95 confidence)</td>
<td>4.45961 ns [4.45864, 4.46063] (0.95 confidence)</td>
<td>0.999739</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/underflow/difference&#95;endpoints</td>
<td>4.60513 ns [4.60452, 4.60726] (0.95 confidence)</td>
<td>4.60251 ns [4.60071, 4.60374] (0.95 confidence)</td>
<td>1.00057</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/underflow/dot</td>
<td>1.48431 ns [1.48401, 1.48453] (0.95 confidence)</td>
<td>1.48455 ns [1.48419, 1.48477] (0.95 confidence)</td>
<td>0.999832</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/underflow/dot&#95;endpoints</td>
<td>1.48468 ns [1.48429, 1.48489] (0.95 confidence)</td>
<td>1.48483 ns [1.48464, 1.485] (0.95 confidence)</td>
<td>0.999894</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;batch/origin/construction</td>
<td>2.29138 ns [2.27493, 2.29623] (0.95 confidence)</td>
<td>2.27551 ns [2.26042, 2.28564] (0.95 confidence)</td>
<td>1.00698</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;batch/origin/prepared</td>
<td>23.2831 ns [23.2732, 23.3209] (0.95 confidence)</td>
<td>23.8129 ns [23.781, 23.8738] (0.95 confidence)</td>
<td>0.977753</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;batch/origin/with&#95;construction</td>
<td>27.795 ns [27.6156, 28.0572] (0.95 confidence)</td>
<td>27.5526 ns [27.4939, 27.6073] (0.95 confidence)</td>
<td>1.0088</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;batch/translated/construction</td>
<td>2.28814 ns [2.27161, 2.29574] (0.95 confidence)</td>
<td>2.27997 ns [2.25537, 2.28992] (0.95 confidence)</td>
<td>1.00358</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;batch/translated/prepared</td>
<td>68.5892 ns [68.464, 68.7469] (0.95 confidence)</td>
<td>69.6704 ns [69.2702, 70.5888] (0.95 confidence)</td>
<td>0.984482</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;batch/translated/with&#95;construction</td>
<td>71.4484 ns [71.3786, 71.5062] (0.95 confidence)</td>
<td>71.6419 ns [71.5292, 71.8792] (0.95 confidence)</td>
<td>0.9973</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/cancellation/difference</td>
<td>10.954 ns [10.9269, 11.0555] (0.95 confidence)</td>
<td>11.0407 ns [10.9513, 11.2173] (0.95 confidence)</td>
<td>0.99215</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/cancellation/difference&#95;endpoints</td>
<td>12.0356 ns [11.9733, 12.0911] (0.95 confidence)</td>
<td>11.937 ns [11.8929, 11.9876] (0.95 confidence)</td>
<td>1.00825</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/cancellation/dot</td>
<td>6.03556 ns [6.02678, 6.04784] (0.95 confidence)</td>
<td>6.14793 ns [6.13156, 6.16787] (0.95 confidence)</td>
<td>0.981721</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/cancellation/dot&#95;endpoints</td>
<td>6.60185 ns [6.58872, 6.63914] (0.95 confidence)</td>
<td>6.68865 ns [6.66838, 6.70632] (0.95 confidence)</td>
<td>0.987023</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense/difference</td>
<td>14.3755 ns [14.3548, 14.3928] (0.95 confidence)</td>
<td>14.421 ns [14.4074, 14.4306] (0.95 confidence)</td>
<td>0.996849</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense/difference&#95;endpoints</td>
<td>16.5617 ns [16.5184, 16.5758] (0.95 confidence)</td>
<td>16.6687 ns [16.6481, 16.6834] (0.95 confidence)</td>
<td>0.99358</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense/dot</td>
<td>8.02696 ns [8.00441, 8.04124] (0.95 confidence)</td>
<td>8.08536 ns [8.05384, 8.11515] (0.95 confidence)</td>
<td>0.992778</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense/dot&#95;endpoints</td>
<td>8.50859 ns [8.4897, 8.55285] (0.95 confidence)</td>
<td>8.50071 ns [8.49015, 8.50737] (0.95 confidence)</td>
<td>1.00093</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense&#95;inexact/difference</td>
<td>14.4063 ns [14.3971, 14.4143] (0.95 confidence)</td>
<td>14.467 ns [14.4493, 14.5296] (0.95 confidence)</td>
<td>0.995802</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense&#95;inexact/difference&#95;endpoints</td>
<td>16.0909 ns [16.0613, 16.126] (0.95 confidence)</td>
<td>16.1282 ns [16.1089, 16.1548] (0.95 confidence)</td>
<td>0.997683</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense&#95;inexact/dot</td>
<td>8.06307 ns [8.04099, 8.08116] (0.95 confidence)</td>
<td>8.08654 ns [8.06816, 8.11039] (0.95 confidence)</td>
<td>0.997097</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense&#95;inexact/dot&#95;endpoints</td>
<td>9.52249 ns [9.49348, 9.55228] (0.95 confidence)</td>
<td>9.55609 ns [9.53034, 9.59582] (0.95 confidence)</td>
<td>0.996484</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/overflow/dot</td>
<td>0.996468 ns [0.994012, 0.998825] (0.95 confidence)</td>
<td>0.997016 ns [0.994736, 0.998419] (0.95 confidence)</td>
<td>0.999451</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/sparse/difference</td>
<td>5.45456 ns [5.45019, 5.45825] (0.95 confidence)</td>
<td>5.46403 ns [5.4614, 5.4703] (0.95 confidence)</td>
<td>0.998267</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/sparse/difference&#95;endpoints</td>
<td>6.37902 ns [6.35013, 6.38763] (0.95 confidence)</td>
<td>6.47811 ns [6.44222, 6.49248] (0.95 confidence)</td>
<td>0.984704</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/sparse/dot</td>
<td>2.72531 ns [2.72461, 2.72701] (0.95 confidence)</td>
<td>2.74687 ns [2.72948, 2.83701] (0.95 confidence)</td>
<td>0.992153</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/sparse/dot&#95;endpoints</td>
<td>2.74395 ns [2.74291, 2.74587] (0.95 confidence)</td>
<td>2.74757 ns [2.74444, 2.77809] (0.95 confidence)</td>
<td>0.998682</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/underflow/difference</td>
<td>4.70393 ns [4.70338, 4.70501] (0.95 confidence)</td>
<td>4.70375 ns [4.70315, 4.7046] (0.95 confidence)</td>
<td>1.00004</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/underflow/difference&#95;endpoints</td>
<td>4.96948 ns [4.96893, 4.97076] (0.95 confidence)</td>
<td>5.1956 ns [5.17362, 5.1969] (0.95 confidence)</td>
<td>0.956478</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/underflow/dot</td>
<td>3.16995 ns [3.16486, 3.17578] (0.95 confidence)</td>
<td>3.18622 ns [3.17932, 3.1926] (0.95 confidence)</td>
<td>0.994894</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/underflow/dot&#95;endpoints</td>
<td>3.25276 ns [3.24051, 3.25528] (0.95 confidence)</td>
<td>3.1958 ns [3.19235, 3.19994] (0.95 confidence)</td>
<td>1.01782</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4/dot&#95;bounded&#95;inconclusive</td>
<td>6.13272 ns [6.12013, 6.142] (0.95 confidence)</td>
<td>6.18361 ns [6.1724, 6.19292] (0.95 confidence)</td>
<td>0.99177</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4/dot&#95;bounded&#95;well&#95;separated</td>
<td>9.49785 ns [9.47855, 9.50645] (0.95 confidence)</td>
<td>9.5294 ns [9.49334, 9.75219] (0.95 confidence)</td>
<td>0.996689</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4/dot&#95;difference&#95;bounded&#95;inconclusive</td>
<td>20.6859 ns [20.6459, 20.7287] (0.95 confidence)</td>
<td>20.8701 ns [20.8014, 20.9866] (0.95 confidence)</td>
<td>0.991172</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4/dot&#95;difference&#95;bounded&#95;well&#95;separated</td>
<td>16.5904 ns [16.5781, 16.6077] (0.95 confidence)</td>
<td>16.6183 ns [16.5687, 16.6877] (0.95 confidence)</td>
<td>0.998317</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4/dot&#95;plain&#95;well&#95;separated</td>
<td>0.681474 ns [0.678887, 0.682175] (0.95 confidence)</td>
<td>0.679748 ns [0.679395, 0.680327] (0.95 confidence)</td>
<td>1.00254</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;batch/origin/construction</td>
<td>2.35595 ns [2.34622, 2.36429] (0.95 confidence)</td>
<td>2.37139 ns [2.36528, 2.37897] (0.95 confidence)</td>
<td>0.993489</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;batch/origin/prepared</td>
<td>34.717 ns [34.5891, 34.878] (0.95 confidence)</td>
<td>35.4356 ns [35.2775, 35.6399] (0.95 confidence)</td>
<td>0.979722</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;batch/origin/with&#95;construction</td>
<td>41.9689 ns [41.8337, 42.1889] (0.95 confidence)</td>
<td>42.3231 ns [42.0715, 42.4351] (0.95 confidence)</td>
<td>0.99163</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;batch/translated/construction</td>
<td>2.35129 ns [2.34207, 2.358] (0.95 confidence)</td>
<td>2.37583 ns [2.36557, 2.38601] (0.95 confidence)</td>
<td>0.989673</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;batch/translated/prepared</td>
<td>109.019 ns [108.879, 109.092] (0.95 confidence)</td>
<td>109.467 ns [109.312, 109.773] (0.95 confidence)</td>
<td>0.995912</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;batch/translated/with&#95;construction</td>
<td>114.927 ns [114.646, 115.239] (0.95 confidence)</td>
<td>115.073 ns [114.818, 115.504] (0.95 confidence)</td>
<td>0.998735</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/cancellation/difference</td>
<td>11.6913 ns [11.6506, 11.8921] (0.95 confidence)</td>
<td>11.685 ns [11.6116, 11.7096] (0.95 confidence)</td>
<td>1.00054</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/cancellation/difference&#95;endpoints</td>
<td>12.7568 ns [12.7248, 12.7956] (0.95 confidence)</td>
<td>12.7277 ns [12.691, 12.8158] (0.95 confidence)</td>
<td>1.00229</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/cancellation/dot</td>
<td>6.33429 ns [6.32613, 6.34807] (0.95 confidence)</td>
<td>6.27657 ns [6.26535, 6.29579] (0.95 confidence)</td>
<td>1.0092</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/cancellation/dot&#95;endpoints</td>
<td>6.91652 ns [6.90428, 6.922] (0.95 confidence)</td>
<td>6.939 ns [6.92545, 6.96078] (0.95 confidence)</td>
<td>0.99676</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense/difference</td>
<td>18.6116 ns [18.5836, 18.6485] (0.95 confidence)</td>
<td>18.5474 ns [18.5141, 18.5837] (0.95 confidence)</td>
<td>1.00346</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense/difference&#95;endpoints</td>
<td>21.0995 ns [21.0785, 21.1231] (0.95 confidence)</td>
<td>20.941 ns [20.9142, 20.9572] (0.95 confidence)</td>
<td>1.00757</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense/dot</td>
<td>9.63862 ns [9.61725, 9.6715] (0.95 confidence)</td>
<td>9.6114 ns [9.57944, 9.6504] (0.95 confidence)</td>
<td>1.00283</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense/dot&#95;endpoints</td>
<td>11.0327 ns [10.9453, 11.1388] (0.95 confidence)</td>
<td>11.2501 ns [11.2302, 11.2596] (0.95 confidence)</td>
<td>0.980679</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense&#95;inexact/difference</td>
<td>18.677 ns [18.6315, 18.7024] (0.95 confidence)</td>
<td>18.6417 ns [18.5985, 19.0384] (0.95 confidence)</td>
<td>1.00189</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense&#95;inexact/difference&#95;endpoints</td>
<td>20.1926 ns [20.1732, 20.2382] (0.95 confidence)</td>
<td>20.2061 ns [20.1859, 20.2462] (0.95 confidence)</td>
<td>0.999331</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense&#95;inexact/dot</td>
<td>9.69772 ns [9.6887, 9.70797] (0.95 confidence)</td>
<td>9.7092 ns [9.69113, 9.77681] (0.95 confidence)</td>
<td>0.998817</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense&#95;inexact/dot&#95;endpoints</td>
<td>10.3547 ns [10.3368, 10.3691] (0.95 confidence)</td>
<td>10.3608 ns [10.3424, 10.3722] (0.95 confidence)</td>
<td>0.99941</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/overflow/dot</td>
<td>1.00295 ns [1.00041, 1.00549] (0.95 confidence)</td>
<td>0.999184 ns [0.996709, 1.00052] (0.95 confidence)</td>
<td>1.00377</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/sparse/difference</td>
<td>6.20857 ns [6.20326, 6.21283] (0.95 confidence)</td>
<td>6.19708 ns [6.19202, 6.20056] (0.95 confidence)</td>
<td>1.00186</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/sparse/difference&#95;endpoints</td>
<td>7.48958 ns [7.48452, 7.49319] (0.95 confidence)</td>
<td>7.42523 ns [7.40587, 7.45718] (0.95 confidence)</td>
<td>1.00867</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/sparse/dot</td>
<td>3.20484 ns [3.15104, 3.21995] (0.95 confidence)</td>
<td>3.12474 ns [3.11944, 3.19853] (0.95 confidence)</td>
<td>1.02563</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/sparse/dot&#95;endpoints</td>
<td>3.23246 ns [3.23039, 3.23518] (0.95 confidence)</td>
<td>3.20557 ns [3.19981, 3.22093] (0.95 confidence)</td>
<td>1.00839</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/underflow/difference</td>
<td>5.69722 ns [5.69483, 5.69968] (0.95 confidence)</td>
<td>5.69191 ns [5.69103, 5.69308] (0.95 confidence)</td>
<td>1.00093</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/underflow/difference&#95;endpoints</td>
<td>5.95124 ns [5.94871, 5.95594] (0.95 confidence)</td>
<td>5.95806 ns [5.94276, 6.03897] (0.95 confidence)</td>
<td>0.998855</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/underflow/dot</td>
<td>3.90652 ns [3.90104, 3.90975] (0.95 confidence)</td>
<td>3.90601 ns [3.89915, 3.91024] (0.95 confidence)</td>
<td>1.00013</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/underflow/dot&#95;endpoints</td>
<td>3.91365 ns [3.90502, 3.91654] (0.95 confidence)</td>
<td>3.90673 ns [3.89944, 3.91276] (0.95 confidence)</td>
<td>1.00177</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;batch/origin/construction</td>
<td>2.90948 ns [2.90236, 2.91577] (0.95 confidence)</td>
<td>2.926 ns [2.92336, 2.92846] (0.95 confidence)</td>
<td>0.994356</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;batch/origin/prepared</td>
<td>46.3272 ns [46.1313, 46.7369] (0.95 confidence)</td>
<td>46.8455 ns [46.7581, 46.9868] (0.95 confidence)</td>
<td>0.988937</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;batch/origin/with&#95;construction</td>
<td>56.4829 ns [56.4383, 56.5745] (0.95 confidence)</td>
<td>58.3742 ns [58.2361, 58.4731] (0.95 confidence)</td>
<td>0.967601</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;batch/translated/construction</td>
<td>3.26714 ns [3.26054, 3.27518] (0.95 confidence)</td>
<td>3.27107 ns [3.26246, 3.27642] (0.95 confidence)</td>
<td>0.998799</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;batch/translated/prepared</td>
<td>145.403 ns [145.306, 145.84] (0.95 confidence)</td>
<td>146.233 ns [146.094, 146.418] (0.95 confidence)</td>
<td>0.994325</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;batch/translated/with&#95;construction</td>
<td>156.542 ns [156.25, 156.726] (0.95 confidence)</td>
<td>157.084 ns [156.878, 157.331] (0.95 confidence)</td>
<td>0.996553</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/cancellation/difference</td>
<td>12.1572 ns [12.1324, 12.1938] (0.95 confidence)</td>
<td>12.0917 ns [12.0837, 12.1182] (0.95 confidence)</td>
<td>1.00542</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/cancellation/difference&#95;endpoints</td>
<td>13.2458 ns [13.2248, 13.2993] (0.95 confidence)</td>
<td>13.2826 ns [13.234, 13.3363] (0.95 confidence)</td>
<td>0.997227</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/cancellation/dot</td>
<td>6.60165 ns [6.59242, 6.62481] (0.95 confidence)</td>
<td>6.62936 ns [6.61498, 6.65726] (0.95 confidence)</td>
<td>0.99582</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/cancellation/dot&#95;endpoints</td>
<td>7.21546 ns [7.20589, 7.23004] (0.95 confidence)</td>
<td>7.32661 ns [7.26474, 7.35456] (0.95 confidence)</td>
<td>0.98483</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense/difference</td>
<td>22.8315 ns [22.7883, 22.8653] (0.95 confidence)</td>
<td>22.8556 ns [22.8336, 22.914] (0.95 confidence)</td>
<td>0.998948</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense/difference&#95;endpoints</td>
<td>24.234 ns [24.1918, 24.2541] (0.95 confidence)</td>
<td>24.2706 ns [24.2356, 24.2886] (0.95 confidence)</td>
<td>0.99849</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense/dot</td>
<td>11.6377 ns [11.6129, 11.6821] (0.95 confidence)</td>
<td>11.6607 ns [11.65, 11.7099] (0.95 confidence)</td>
<td>0.998027</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense/dot&#95;endpoints</td>
<td>12.3236 ns [12.3057, 12.3383] (0.95 confidence)</td>
<td>12.378 ns [12.3644, 12.4029] (0.95 confidence)</td>
<td>0.995605</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense&#95;inexact/difference</td>
<td>22.8572 ns [22.8311, 22.88] (0.95 confidence)</td>
<td>23.0161 ns [22.9778, 23.044] (0.95 confidence)</td>
<td>0.9931</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense&#95;inexact/difference&#95;endpoints</td>
<td>25.3199 ns [25.2659, 25.3734] (0.95 confidence)</td>
<td>25.418 ns [25.389, 25.4773] (0.95 confidence)</td>
<td>0.996143</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense&#95;inexact/dot</td>
<td>11.6891 ns [11.6617, 11.6979] (0.95 confidence)</td>
<td>11.7022 ns [11.6645, 11.7322] (0.95 confidence)</td>
<td>0.998883</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense&#95;inexact/dot&#95;endpoints</td>
<td>13.2363 ns [13.2219, 13.2626] (0.95 confidence)</td>
<td>13.3147 ns [13.2973, 13.3482] (0.95 confidence)</td>
<td>0.994108</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/overflow/dot</td>
<td>0.994167 ns [0.992614, 0.995279] (0.95 confidence)</td>
<td>0.999834 ns [0.997773, 1.00188] (0.95 confidence)</td>
<td>0.994332</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/sparse/difference</td>
<td>7.2102 ns [7.20756, 7.21498] (0.95 confidence)</td>
<td>7.21538 ns [7.20869, 7.22104] (0.95 confidence)</td>
<td>0.999282</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/sparse/difference&#95;endpoints</td>
<td>8.37067 ns [8.34733, 8.38948] (0.95 confidence)</td>
<td>8.34408 ns [8.33145, 8.36214] (0.95 confidence)</td>
<td>1.00319</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/sparse/dot</td>
<td>3.63274 ns [3.62453, 3.66723] (0.95 confidence)</td>
<td>3.63848 ns [3.63125, 3.65012] (0.95 confidence)</td>
<td>0.998421</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/sparse/dot&#95;endpoints</td>
<td>3.74578 ns [3.74095, 3.75712] (0.95 confidence)</td>
<td>3.75081 ns [3.74668, 3.76219] (0.95 confidence)</td>
<td>0.998658</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/underflow/difference</td>
<td>7.68173 ns [7.66676, 7.70563] (0.95 confidence)</td>
<td>7.72672 ns [7.71315, 7.73816] (0.95 confidence)</td>
<td>0.994177</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/underflow/difference&#95;endpoints</td>
<td>8.15954 ns [8.15105, 8.16446] (0.95 confidence)</td>
<td>8.39012 ns [8.37311, 8.4203] (0.95 confidence)</td>
<td>0.972517</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/underflow/dot</td>
<td>4.62012 ns [4.61281, 4.62611] (0.95 confidence)</td>
<td>4.64726 ns [4.63784, 4.6528] (0.95 confidence)</td>
<td>0.99416</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/underflow/dot&#95;endpoints</td>
<td>4.62062 ns [4.61463, 4.62659] (0.95 confidence)</td>
<td>4.65546 ns [4.64721, 4.65954] (0.95 confidence)</td>
<td>0.992516</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;batch/origin/construction</td>
<td>3.10651 ns [3.10201, 3.11101] (0.95 confidence)</td>
<td>3.1121 ns [3.10656, 3.11485] (0.95 confidence)</td>
<td>0.998204</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;batch/origin/prepared</td>
<td>60.9242 ns [60.8085, 61.1343] (0.95 confidence)</td>
<td>60.9984 ns [60.8603, 61.2142] (0.95 confidence)</td>
<td>0.998783</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;batch/origin/with&#95;construction</td>
<td>81.9307 ns [81.4382, 82.7448] (0.95 confidence)</td>
<td>80.977 ns [80.8887, 81.2723] (0.95 confidence)</td>
<td>1.01178</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;batch/translated/construction</td>
<td>3.5917 ns [3.58771, 3.59549] (0.95 confidence)</td>
<td>3.59155 ns [3.5881, 3.59563] (0.95 confidence)</td>
<td>1.00004</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;batch/translated/prepared</td>
<td>204.227 ns [203.703, 204.63] (0.95 confidence)</td>
<td>206.204 ns [204.273, 207.895] (0.95 confidence)</td>
<td>0.990411</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;batch/translated/with&#95;construction</td>
<td>230.032 ns [228.806, 231.195] (0.95 confidence)</td>
<td>227.561 ns [227.012, 228.006] (0.95 confidence)</td>
<td>1.01086</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/cancellation/difference</td>
<td>12.6174 ns [12.5962, 12.6413] (0.95 confidence)</td>
<td>12.5249 ns [12.5113, 12.5367] (0.95 confidence)</td>
<td>1.00739</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/cancellation/difference&#95;endpoints</td>
<td>14.0021 ns [13.9931, 14.0164] (0.95 confidence)</td>
<td>14.0287 ns [14.0084, 14.0461] (0.95 confidence)</td>
<td>0.998101</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/cancellation/dot</td>
<td>6.90846 ns [6.89136, 6.9243] (0.95 confidence)</td>
<td>7.05878 ns [6.98221, 7.09927] (0.95 confidence)</td>
<td>0.978705</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/cancellation/dot&#95;endpoints</td>
<td>7.59503 ns [7.59208, 7.60006] (0.95 confidence)</td>
<td>7.64698 ns [7.62876, 7.69179] (0.95 confidence)</td>
<td>0.993206</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense/difference</td>
<td>29.0409 ns [28.999, 29.0775] (0.95 confidence)</td>
<td>29.1737 ns [29.049, 29.2236] (0.95 confidence)</td>
<td>0.995448</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense/difference&#95;endpoints</td>
<td>32.4653 ns [32.3431, 32.5039] (0.95 confidence)</td>
<td>32.5806 ns [32.5462, 32.6396] (0.95 confidence)</td>
<td>0.996459</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense/dot</td>
<td>13.83 ns [13.8104, 13.8389] (0.95 confidence)</td>
<td>13.8503 ns [13.821, 13.8675] (0.95 confidence)</td>
<td>0.998536</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense/dot&#95;endpoints</td>
<td>15.2514 ns [15.2227, 15.2659] (0.95 confidence)</td>
<td>15.2345 ns [15.1893, 15.2537] (0.95 confidence)</td>
<td>1.00111</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense&#95;inexact/difference</td>
<td>29.0829 ns [29.0543, 29.1176] (0.95 confidence)</td>
<td>29.1489 ns [29.0821, 29.243] (0.95 confidence)</td>
<td>0.997736</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense&#95;inexact/difference&#95;endpoints</td>
<td>31.6463 ns [31.5923, 31.7009] (0.95 confidence)</td>
<td>31.9465 ns [31.8482, 32.0016] (0.95 confidence)</td>
<td>0.990602</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense&#95;inexact/dot</td>
<td>13.8366 ns [13.82, 13.8499] (0.95 confidence)</td>
<td>13.8648 ns [13.8359, 13.9007] (0.95 confidence)</td>
<td>0.997962</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense&#95;inexact/dot&#95;endpoints</td>
<td>15.2667 ns [15.245, 15.2884] (0.95 confidence)</td>
<td>15.2557 ns [15.2344, 15.3408] (0.95 confidence)</td>
<td>1.00072</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/overflow/dot</td>
<td>0.996479 ns [0.995208, 0.999786] (0.95 confidence)</td>
<td>1.00055 ns [0.999182, 1.0021] (0.95 confidence)</td>
<td>0.995931</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/sparse/difference</td>
<td>8.17102 ns [8.14882, 8.18237] (0.95 confidence)</td>
<td>8.1829 ns [8.17088, 8.19252] (0.95 confidence)</td>
<td>0.998549</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/sparse/difference&#95;endpoints</td>
<td>9.32462 ns [9.31832, 9.36817] (0.95 confidence)</td>
<td>9.49495 ns [9.47546, 9.51595] (0.95 confidence)</td>
<td>0.98206</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/sparse/dot</td>
<td>4.108 ns [4.10456, 4.11065] (0.95 confidence)</td>
<td>4.11447 ns [4.11149, 4.11776] (0.95 confidence)</td>
<td>0.998428</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/sparse/dot&#95;endpoints</td>
<td>4.17631 ns [4.17479, 4.17692] (0.95 confidence)</td>
<td>4.18986 ns [4.1875, 4.19191] (0.95 confidence)</td>
<td>0.996768</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/underflow/difference</td>
<td>8.63523 ns [8.62149, 8.64722] (0.95 confidence)</td>
<td>8.67195 ns [8.66222, 8.67941] (0.95 confidence)</td>
<td>0.995765</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/underflow/difference&#95;endpoints</td>
<td>8.94743 ns [8.92993, 9.03814] (0.95 confidence)</td>
<td>8.94217 ns [8.91777, 8.95607] (0.95 confidence)</td>
<td>1.00059</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/underflow/dot</td>
<td>5.33089 ns [5.32493, 5.33677] (0.95 confidence)</td>
<td>5.37087 ns [5.35432, 5.3847] (0.95 confidence)</td>
<td>0.992556</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/underflow/dot&#95;endpoints</td>
<td>5.39812 ns [5.38924, 5.40795] (0.95 confidence)</td>
<td>5.42274 ns [5.41357, 5.43196] (0.95 confidence)</td>
<td>0.99546</td>
</tr>
<tr>
<td>public&#95;d2/det&#95;bound</td>
<td>1.95664 ns [1.95217, 1.95847] (0.95 confidence)</td>
<td>1.96514 ns [1.96354, 1.96708] (0.95 confidence)</td>
<td>0.995675</td>
</tr>
<tr>
<td>public&#95;d2/det&#95;direct</td>
<td>1.13468 ns [1.13325, 1.13537] (0.95 confidence)</td>
<td>1.13465 ns [1.13358, 1.13625] (0.95 confidence)</td>
<td>1.00003</td>
</tr>
<tr>
<td>public&#95;d2/det&#95;with&#95;bound</td>
<td>2.12183 ns [2.1143, 2.12778] (0.95 confidence)</td>
<td>2.11814 ns [2.11333, 2.12301] (0.95 confidence)</td>
<td>1.00174</td>
</tr>
<tr>
<td>public&#95;d2/first&#95;asymmetry</td>
<td>1.3834 ns [1.38226, 1.3841] (0.95 confidence)</td>
<td>1.38564 ns [1.38442, 1.3869] (0.95 confidence)</td>
<td>0.998383</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;det</td>
<td>7.7934 ns [7.71027, 7.91615] (0.95 confidence)</td>
<td>7.9311 ns [7.8421, 7.94839] (0.95 confidence)</td>
<td>0.982638</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;access</td>
<td>6.36566 ns [6.34568, 6.37962] (0.95 confidence)</td>
<td>6.37664 ns [6.36432, 6.38653] (0.95 confidence)</td>
<td>0.998278</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;from&#95;matrix</td>
<td>2.5997 ns [2.59268, 2.6124] (0.95 confidence)</td>
<td>2.6046 ns [2.59733, 2.61173] (0.95 confidence)</td>
<td>0.998119</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;from&#95;points</td>
<td>2.24632 ns [2.23206, 2.26519] (0.95 confidence)</td>
<td>2.25494 ns [2.2339, 2.26778] (0.95 confidence)</td>
<td>0.996181</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;from&#95;rows</td>
<td>3.5809 ns [3.57694, 3.58891] (0.95 confidence)</td>
<td>3.58125 ns [3.5736, 3.58615] (0.95 confidence)</td>
<td>0.999903</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;identity</td>
<td>2.10121 ns [2.09323, 2.10916] (0.95 confidence)</td>
<td>2.10873 ns [2.10491, 2.1147] (0.95 confidence)</td>
<td>0.996433</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;set</td>
<td>3.8324 ns [3.82261, 3.83859] (0.95 confidence)</td>
<td>3.81699 ns [3.80735, 3.8273] (0.95 confidence)</td>
<td>1.00404</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;zero</td>
<td>2.11183 ns [2.1071, 2.11879] (0.95 confidence)</td>
<td>2.12562 ns [2.11983, 2.12962] (0.95 confidence)</td>
<td>0.993514</td>
</tr>
<tr>
<td>public&#95;d2/is&#95;symmetric</td>
<td>0.757479 ns [0.756805, 0.757957] (0.95 confidence)</td>
<td>0.759338 ns [0.758783, 0.759927] (0.95 confidence)</td>
<td>0.997552</td>
</tr>
<tr>
<td>public&#95;d2/matrix&#95;access</td>
<td>4.76058 ns [4.74843, 4.77499] (0.95 confidence)</td>
<td>4.7767 ns [4.76554, 4.78224] (0.95 confidence)</td>
<td>0.996625</td>
</tr>
<tr>
<td>public&#95;d2/matrix&#95;from&#95;rows</td>
<td>2.913 ns [2.90965, 2.91584] (0.95 confidence)</td>
<td>2.92204 ns [2.91696, 2.92716] (0.95 confidence)</td>
<td>0.996907</td>
</tr>
<tr>
<td>public&#95;d2/matrix&#95;identity</td>
<td>1.25949 ns [1.25732, 1.2625] (0.95 confidence)</td>
<td>1.2613 ns [1.25843, 1.26469] (0.95 confidence)</td>
<td>0.998568</td>
</tr>
<tr>
<td>public&#95;d2/matrix&#95;set</td>
<td>2.12336 ns [2.10641, 2.13044] (0.95 confidence)</td>
<td>2.10148 ns [2.09364, 2.11055] (0.95 confidence)</td>
<td>1.01041</td>
</tr>
<tr>
<td>public&#95;d2/matrix&#95;zero</td>
<td>1.16224 ns [1.15933, 1.16537] (0.95 confidence)</td>
<td>1.1658 ns [1.1636, 1.16918] (0.95 confidence)</td>
<td>0.99695</td>
</tr>
<tr>
<td>public&#95;d2/vector&#95;access</td>
<td>1.61911 ns [1.61569, 1.62371] (0.95 confidence)</td>
<td>1.63055 ns [1.62709, 1.63432] (0.95 confidence)</td>
<td>0.992987</td>
</tr>
<tr>
<td>public&#95;d2/vector&#95;new</td>
<td>1.94095 ns [1.93805, 1.94432] (0.95 confidence)</td>
<td>1.94043 ns [1.93127, 1.94987] (0.95 confidence)</td>
<td>1.00027</td>
</tr>
<tr>
<td>public&#95;d2/vector&#95;zero</td>
<td>0.591446 ns [0.590549, 0.59294] (0.95 confidence)</td>
<td>0.594682 ns [0.592713, 0.596628] (0.95 confidence)</td>
<td>0.994557</td>
</tr>
<tr>
<td>public&#95;d3/det&#95;bound</td>
<td>3.66565 ns [3.6644, 3.66757] (0.95 confidence)</td>
<td>3.65778 ns [3.64737, 3.66423] (0.95 confidence)</td>
<td>1.00215</td>
</tr>
<tr>
<td>public&#95;d3/det&#95;direct</td>
<td>1.28711 ns [1.28673, 1.28802] (0.95 confidence)</td>
<td>1.28927 ns [1.28757, 1.29206] (0.95 confidence)</td>
<td>0.998328</td>
</tr>
<tr>
<td>public&#95;d3/det&#95;with&#95;bound</td>
<td>3.67744 ns [3.67235, 3.68009] (0.95 confidence)</td>
<td>3.66553 ns [3.6613, 3.67118] (0.95 confidence)</td>
<td>1.00325</td>
</tr>
<tr>
<td>public&#95;d3/first&#95;asymmetry</td>
<td>3.52765 ns [3.50953, 3.56052] (0.95 confidence)</td>
<td>3.55373 ns [3.5328, 3.58731] (0.95 confidence)</td>
<td>0.992662</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;det</td>
<td>23.6946 ns [23.1601, 24.0108] (0.95 confidence)</td>
<td>24.1336 ns [23.8254, 24.311] (0.95 confidence)</td>
<td>0.981811</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;access</td>
<td>10.0263 ns [9.9967, 10.0342] (0.95 confidence)</td>
<td>10.162 ns [10.1461, 10.182] (0.95 confidence)</td>
<td>0.986641</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;from&#95;matrix</td>
<td>5.18056 ns [5.16745, 5.18971] (0.95 confidence)</td>
<td>5.19655 ns [5.18722, 5.20697] (0.95 confidence)</td>
<td>0.996924</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;from&#95;points</td>
<td>6.3321 ns [6.31716, 6.33924] (0.95 confidence)</td>
<td>6.41037 ns [6.39382, 6.41986] (0.95 confidence)</td>
<td>0.98779</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;from&#95;rows</td>
<td>7.37786 ns [7.36395, 7.39881] (0.95 confidence)</td>
<td>7.42294 ns [7.40042, 7.45617] (0.95 confidence)</td>
<td>0.993927</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;identity</td>
<td>4.33362 ns [4.32213, 4.34066] (0.95 confidence)</td>
<td>4.47483 ns [4.439, 4.56631] (0.95 confidence)</td>
<td>0.968443</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;set</td>
<td>7.82044 ns [7.80034, 7.84705] (0.95 confidence)</td>
<td>7.93497 ns [7.89658, 7.95434] (0.95 confidence)</td>
<td>0.985567</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;zero</td>
<td>4.32661 ns [4.32013, 4.33625] (0.95 confidence)</td>
<td>4.37214 ns [4.34871, 4.38037] (0.95 confidence)</td>
<td>0.989585</td>
</tr>
<tr>
<td>public&#95;d3/is&#95;symmetric</td>
<td>1.49044 ns [1.48743, 1.49315] (0.95 confidence)</td>
<td>1.34355 ns [1.34255, 1.34523] (0.95 confidence)</td>
<td>1.10933</td>
</tr>
<tr>
<td>public&#95;d3/matrix&#95;access</td>
<td>6.69919 ns [6.66677, 6.71384] (0.95 confidence)</td>
<td>6.97709 ns [6.9349, 7.01896] (0.95 confidence)</td>
<td>0.96017</td>
</tr>
<tr>
<td>public&#95;d3/matrix&#95;from&#95;rows</td>
<td>5.70964 ns [5.68825, 5.73631] (0.95 confidence)</td>
<td>5.73515 ns [5.71534, 5.74878] (0.95 confidence)</td>
<td>0.995552</td>
</tr>
<tr>
<td>public&#95;d3/matrix&#95;identity</td>
<td>2.7134 ns [2.70193, 2.72176] (0.95 confidence)</td>
<td>2.72064 ns [2.7163, 2.73016] (0.95 confidence)</td>
<td>0.997339</td>
</tr>
<tr>
<td>public&#95;d3/matrix&#95;set</td>
<td>4.25673 ns [4.2347, 4.2715] (0.95 confidence)</td>
<td>4.34403 ns [4.318, 4.36906] (0.95 confidence)</td>
<td>0.979904</td>
</tr>
<tr>
<td>public&#95;d3/matrix&#95;zero</td>
<td>2.41291 ns [2.40647, 2.42039] (0.95 confidence)</td>
<td>2.43336 ns [2.42714, 2.44021] (0.95 confidence)</td>
<td>0.991599</td>
</tr>
<tr>
<td>public&#95;d3/vector&#95;access</td>
<td>2.09708 ns [2.09259, 2.10164] (0.95 confidence)</td>
<td>2.10683 ns [2.10284, 2.11069] (0.95 confidence)</td>
<td>0.995369</td>
</tr>
<tr>
<td>public&#95;d3/vector&#95;new</td>
<td>3.65061 ns [3.63946, 3.66285] (0.95 confidence)</td>
<td>3.64449 ns [3.6351, 3.65915] (0.95 confidence)</td>
<td>1.00168</td>
</tr>
<tr>
<td>public&#95;d3/vector&#95;zero</td>
<td>0.959549 ns [0.957472, 0.961661] (0.95 confidence)</td>
<td>0.968004 ns [0.966127, 0.970576] (0.95 confidence)</td>
<td>0.991266</td>
</tr>
<tr>
<td>public&#95;d4/det&#95;bound</td>
<td>6.1915 ns [6.1852, 6.19664] (0.95 confidence)</td>
<td>7.00967 ns [6.98388, 7.02617] (0.95 confidence)</td>
<td>0.88328</td>
</tr>
<tr>
<td>public&#95;d4/det&#95;direct</td>
<td>2.47648 ns [2.47504, 2.47766] (0.95 confidence)</td>
<td>2.80647 ns [2.79169, 2.81609] (0.95 confidence)</td>
<td>0.88242</td>
</tr>
<tr>
<td>public&#95;d4/det&#95;with&#95;bound</td>
<td>6.1938 ns [6.19142, 6.19521] (0.95 confidence)</td>
<td>6.20103 ns [6.19678, 6.20556] (0.95 confidence)</td>
<td>0.998835</td>
</tr>
<tr>
<td>public&#95;d4/first&#95;asymmetry</td>
<td>3.62292 ns [3.6024, 3.64156] (0.95 confidence)</td>
<td>3.65411 ns [3.60472, 3.68347] (0.95 confidence)</td>
<td>0.991463</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;det</td>
<td>48.3828 ns [48.2087, 48.6253] (0.95 confidence)</td>
<td>49.3042 ns [48.6513, 51.0911] (0.95 confidence)</td>
<td>0.981313</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;access</td>
<td>13.5995 ns [13.5306, 13.7773] (0.95 confidence)</td>
<td>13.8055 ns [13.661, 13.896] (0.95 confidence)</td>
<td>0.985073</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;from&#95;matrix</td>
<td>7.94081 ns [7.92476, 7.95452] (0.95 confidence)</td>
<td>8.05866 ns [8.04007, 8.07749] (0.95 confidence)</td>
<td>0.985376</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;from&#95;points</td>
<td>10.8896 ns [10.8737, 10.9027] (0.95 confidence)</td>
<td>11.0298 ns [11.0115, 11.0441] (0.95 confidence)</td>
<td>0.987292</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;from&#95;rows</td>
<td>11.2131 ns [11.1817, 11.2282] (0.95 confidence)</td>
<td>11.4555 ns [11.4303, 11.4704] (0.95 confidence)</td>
<td>0.978843</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;identity</td>
<td>6.53603 ns [6.51539, 6.54871] (0.95 confidence)</td>
<td>6.66892 ns [6.64473, 6.69073] (0.95 confidence)</td>
<td>0.980073</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;set</td>
<td>11.4694 ns [11.4416, 11.5012] (0.95 confidence)</td>
<td>11.6451 ns [11.6143, 11.6584] (0.95 confidence)</td>
<td>0.984916</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;zero</td>
<td>6.62053 ns [6.61243, 6.64061] (0.95 confidence)</td>
<td>6.78746 ns [6.7567, 6.80812] (0.95 confidence)</td>
<td>0.975405</td>
</tr>
<tr>
<td>public&#95;d4/is&#95;symmetric</td>
<td>2.10627 ns [2.10574, 2.10664] (0.95 confidence)</td>
<td>2.10835 ns [2.10702, 2.11084] (0.95 confidence)</td>
<td>0.999016</td>
</tr>
<tr>
<td>public&#95;d4/matrix&#95;access</td>
<td>9.44024 ns [9.41604, 9.45878] (0.95 confidence)</td>
<td>9.72293 ns [9.69976, 9.73757] (0.95 confidence)</td>
<td>0.970925</td>
</tr>
<tr>
<td>public&#95;d4/matrix&#95;from&#95;rows</td>
<td>10.3309 ns [10.3165, 10.3365] (0.95 confidence)</td>
<td>10.4825 ns [10.4636, 10.5216] (0.95 confidence)</td>
<td>0.985536</td>
</tr>
<tr>
<td>public&#95;d4/matrix&#95;identity</td>
<td>4.12042 ns [4.11093, 4.13012] (0.95 confidence)</td>
<td>4.33765 ns [4.284, 4.47742] (0.95 confidence)</td>
<td>0.94992</td>
</tr>
<tr>
<td>public&#95;d4/matrix&#95;set</td>
<td>7.36531 ns [7.35086, 7.38621] (0.95 confidence)</td>
<td>7.63363 ns [7.62135, 7.65548] (0.95 confidence)</td>
<td>0.964851</td>
</tr>
<tr>
<td>public&#95;d4/matrix&#95;zero</td>
<td>3.87272 ns [3.86636, 3.87808] (0.95 confidence)</td>
<td>3.96141 ns [3.95578, 3.97897] (0.95 confidence)</td>
<td>0.977611</td>
</tr>
<tr>
<td>public&#95;d4/vector&#95;access</td>
<td>2.62714 ns [2.62292, 2.63186] (0.95 confidence)</td>
<td>2.67881 ns [2.67184, 2.68401] (0.95 confidence)</td>
<td>0.980712</td>
</tr>
<tr>
<td>public&#95;d4/vector&#95;new</td>
<td>2.8288 ns [2.82171, 2.83659] (0.95 confidence)</td>
<td>2.85693 ns [2.83743, 2.86533] (0.95 confidence)</td>
<td>0.990153</td>
</tr>
<tr>
<td>public&#95;d4/vector&#95;zero</td>
<td>1.17901 ns [1.17583, 1.18094] (0.95 confidence)</td>
<td>1.20705 ns [1.20424, 1.20951] (0.95 confidence)</td>
<td>0.97677</td>
</tr>
<tr>
<td>public&#95;d5/det&#95;bound</td>
<td>1.3395 ns [1.33841, 1.34107] (0.95 confidence)</td>
<td>1.35904 ns [1.3557, 1.36239] (0.95 confidence)</td>
<td>0.985619</td>
</tr>
<tr>
<td>public&#95;d5/det&#95;direct</td>
<td>1.33753 ns [1.33415, 1.34082] (0.95 confidence)</td>
<td>1.36889 ns [1.36439, 1.37485] (0.95 confidence)</td>
<td>0.977093</td>
</tr>
<tr>
<td>public&#95;d5/det&#95;with&#95;bound</td>
<td>0.856813 ns [0.85442, 0.858827] (0.95 confidence)</td>
<td>0.868672 ns [0.866184, 0.870306] (0.95 confidence)</td>
<td>0.986349</td>
</tr>
<tr>
<td>public&#95;d5/first&#95;asymmetry</td>
<td>3.78565 ns [3.78332, 3.7871] (0.95 confidence)</td>
<td>3.78622 ns [3.7833, 3.78859] (0.95 confidence)</td>
<td>0.99985</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;det</td>
<td>125.547 ns [124.859, 126.736] (0.95 confidence)</td>
<td>133.391 ns [131.845, 137.855] (0.95 confidence)</td>
<td>0.941196</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;access</td>
<td>16.2943 ns [15.8506, 16.4095] (0.95 confidence)</td>
<td>15.9283 ns [15.792, 16.3053] (0.95 confidence)</td>
<td>1.02298</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;from&#95;matrix</td>
<td>10.6943 ns [10.6274, 10.7548] (0.95 confidence)</td>
<td>10.7694 ns [10.7102, 10.8119] (0.95 confidence)</td>
<td>0.993022</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;from&#95;points</td>
<td>18.4907 ns [18.4749, 18.5202] (0.95 confidence)</td>
<td>18.5823 ns [18.5577, 18.6246] (0.95 confidence)</td>
<td>0.995067</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;from&#95;rows</td>
<td>12.6471 ns [12.641, 12.6962] (0.95 confidence)</td>
<td>12.709 ns [12.6708, 12.736] (0.95 confidence)</td>
<td>0.995126</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;identity</td>
<td>8.39408 ns [8.0575, 8.66217] (0.95 confidence)</td>
<td>8.55546 ns [8.30833, 8.76291] (0.95 confidence)</td>
<td>0.981137</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;set</td>
<td>12.8861 ns [12.8658, 13.0514] (0.95 confidence)</td>
<td>13.0077 ns [12.95, 13.2063] (0.95 confidence)</td>
<td>0.990645</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;zero</td>
<td>9.34898 ns [9.33168, 9.36607] (0.95 confidence)</td>
<td>9.59338 ns [9.56535, 9.63931] (0.95 confidence)</td>
<td>0.974524</td>
</tr>
<tr>
<td>public&#95;d5/is&#95;symmetric</td>
<td>3.69191 ns [3.69117, 3.69297] (0.95 confidence)</td>
<td>3.69208 ns [3.691, 3.69335] (0.95 confidence)</td>
<td>0.999953</td>
</tr>
<tr>
<td>public&#95;d5/matrix&#95;access</td>
<td>12.1451 ns [12.0202, 12.3277] (0.95 confidence)</td>
<td>12.3022 ns [12.1327, 12.4112] (0.95 confidence)</td>
<td>0.987232</td>
</tr>
<tr>
<td>public&#95;d5/matrix&#95;from&#95;rows</td>
<td>11.7925 ns [11.7794, 11.8023] (0.95 confidence)</td>
<td>11.8705 ns [11.8627, 11.8806] (0.95 confidence)</td>
<td>0.993432</td>
</tr>
<tr>
<td>public&#95;d5/matrix&#95;identity</td>
<td>6.80234 ns [6.7697, 6.828] (0.95 confidence)</td>
<td>6.92828 ns [6.90834, 6.96514] (0.95 confidence)</td>
<td>0.981823</td>
</tr>
<tr>
<td>public&#95;d5/matrix&#95;set</td>
<td>10.4954 ns [10.4623, 10.5237] (0.95 confidence)</td>
<td>10.6193 ns [10.5849, 10.6446] (0.95 confidence)</td>
<td>0.988337</td>
</tr>
<tr>
<td>public&#95;d5/matrix&#95;zero</td>
<td>5.70736 ns [5.693, 5.72254] (0.95 confidence)</td>
<td>5.86593 ns [5.85234, 5.87711] (0.95 confidence)</td>
<td>0.972969</td>
</tr>
<tr>
<td>public&#95;d5/vector&#95;access</td>
<td>3.1004 ns [3.09257, 3.10462] (0.95 confidence)</td>
<td>3.10303 ns [3.09797, 3.10622] (0.95 confidence)</td>
<td>0.999155</td>
</tr>
<tr>
<td>public&#95;d5/vector&#95;new</td>
<td>4.86258 ns [4.85467, 4.86847] (0.95 confidence)</td>
<td>4.86773 ns [4.85818, 4.87585] (0.95 confidence)</td>
<td>0.998942</td>
</tr>
<tr>
<td>public&#95;d5/vector&#95;zero</td>
<td>1.51493 ns [1.50675, 1.52515] (0.95 confidence)</td>
<td>1.54872 ns [1.53332, 1.55728] (0.95 confidence)</td>
<td>0.978185</td>
</tr>
<tr>
<td>public&#95;scalar/interval&#95;access</td>
<td>1.59363 ns [1.58853, 1.59584] (0.95 confidence)</td>
<td>1.59533 ns [1.59273, 1.59761] (0.95 confidence)</td>
<td>0.998936</td>
</tr>
<tr>
<td>public&#95;scalar/interval&#95;negate</td>
<td>0.733931 ns [0.732952, 0.734547] (0.95 confidence)</td>
<td>0.734398 ns [0.733084, 0.734993] (0.95 confidence)</td>
<td>0.999365</td>
</tr>
<tr>
<td>public&#95;scalar/interval&#95;new</td>
<td>0.893183 ns [0.892522, 0.894376] (0.95 confidence)</td>
<td>0.893373 ns [0.89285, 0.894503] (0.95 confidence)</td>
<td>0.999787</td>
</tr>
<tr>
<td>public&#95;scalar/interval&#95;point</td>
<td>0.582607 ns [0.581797, 0.583414] (0.95 confidence)</td>
<td>0.583195 ns [0.582302, 0.584175] (0.95 confidence)</td>
<td>0.998992</td>
</tr>
<tr>
<td>public&#95;scalar/interval&#95;subtraction</td>
<td>1.86419 ns [1.8621, 1.86615] (0.95 confidence)</td>
<td>1.87035 ns [1.86576, 1.87564] (0.95 confidence)</td>
<td>0.996711</td>
</tr>
<tr>
<td>public&#95;scalar/tolerance</td>
<td>0.979754 ns [0.977365, 0.982309] (0.95 confidence)</td>
<td>0.981141 ns [0.978186, 0.984157] (0.95 confidence)</td>
<td>0.998585</td>
</tr>
</tbody>
</table>

## core-exact

Semantic contract: la-stack-core-exact-v1. Collection session: `f846a217f1914228a0116b2849d1b360`.
Collection interval (Unix ns): 1794375424898021000–1794387068765900000.
Acquisition: newly collected.

- baseline: available — inside declared release interval; release v0.5.2; revision `eb5e44e52399401311e676faad590e6cce1610e5`.
  Source SHA-256: `bdd47c0e9bc3126dc3a32a98e0f1a688aaa52e7c6a3868009f1adfb0735f367e`; harness SHA-256:
  `0fedb76c77086b1cf8a698034ac4a1f3539fe06feda597672e2d11f01704f9cb`; scientific gate: passed.
- current: available — inside declared release interval; release v0.5.3; revision `be9277bf04400dccbe1115da69c7fda2141c4a39`.
  Source SHA-256: `793169b4c5085624d7c0e773f1d28504c017e9b78332f250150d19df87391a32`; harness SHA-256:
  `0fedb76c77086b1cf8a698034ac4a1f3539fe06feda597672e2d11f01704f9cb`; scientific gate: passed.

### Mean (ns)

<table>
<thead>
<tr>
<th>Workload</th>
<th>Baseline</th>
<th>Current</th>
<th>Baseline / current</th>
</tr>
</thead>
<tbody>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d2/rounded&#95;array</td>
<td>170.837 ns [170.283, 171.462] (0.95 confidence)</td>
<td>168.083 ns [167.549, 168.638] (0.95 confidence)</td>
<td>1.01639</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d2/rounded&#95;result</td>
<td>170.742 ns [170.212, 171.296] (0.95 confidence)</td>
<td>169.237 ns [168.564, 169.939] (0.95 confidence)</td>
<td>1.00889</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d2/strict&#95;array</td>
<td>175.229 ns [174.65, 175.837] (0.95 confidence)</td>
<td>171.103 ns [170.723, 171.495] (0.95 confidence)</td>
<td>1.02412</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d2/strict&#95;result</td>
<td>178.763 ns [178.139, 179.429] (0.95 confidence)</td>
<td>174.233 ns [173.371, 175.148] (0.95 confidence)</td>
<td>1.026</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d3/rounded&#95;array</td>
<td>180.657 ns [179.088, 182.316] (0.95 confidence)</td>
<td>177.818 ns [177.009, 178.653] (0.95 confidence)</td>
<td>1.01597</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d3/rounded&#95;result</td>
<td>176.993 ns [176.24, 177.803] (0.95 confidence)</td>
<td>173.937 ns [173.558, 174.325] (0.95 confidence)</td>
<td>1.01757</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d3/strict&#95;array</td>
<td>180.657 ns [179.907, 181.474] (0.95 confidence)</td>
<td>176.343 ns [175.645, 177.077] (0.95 confidence)</td>
<td>1.02446</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d3/strict&#95;result</td>
<td>175.064 ns [174.542, 175.603] (0.95 confidence)</td>
<td>176.081 ns [175.515, 176.65] (0.95 confidence)</td>
<td>0.994223</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d4/rounded&#95;array</td>
<td>185.885 ns [185.47, 186.311] (0.95 confidence)</td>
<td>190.444 ns [188.506, 192.399] (0.95 confidence)</td>
<td>0.976061</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d4/rounded&#95;result</td>
<td>185.475 ns [184.95, 185.993] (0.95 confidence)</td>
<td>188.844 ns [187.164, 190.553] (0.95 confidence)</td>
<td>0.982156</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d4/strict&#95;array</td>
<td>187.45 ns [186.894, 187.993] (0.95 confidence)</td>
<td>192.31 ns [190.319, 194.287] (0.95 confidence)</td>
<td>0.97473</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d4/strict&#95;result</td>
<td>181.995 ns [181.402, 182.573] (0.95 confidence)</td>
<td>190.868 ns [189.375, 192.37] (0.95 confidence)</td>
<td>0.95351</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d5/rounded&#95;array</td>
<td>187.538 ns [185.857, 189.438] (0.95 confidence)</td>
<td>184.701 ns [184.327, 185.093] (0.95 confidence)</td>
<td>1.01536</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d5/rounded&#95;result</td>
<td>186.546 ns [185.966, 187.144] (0.95 confidence)</td>
<td>185.773 ns [185.262, 186.305] (0.95 confidence)</td>
<td>1.00416</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d5/strict&#95;array</td>
<td>183.692 ns [183.222, 184.185] (0.95 confidence)</td>
<td>184.786 ns [184.367, 185.217] (0.95 confidence)</td>
<td>0.994079</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d5/strict&#95;result</td>
<td>181.015 ns [180.578, 181.462] (0.95 confidence)</td>
<td>181.317 ns [180.85, 181.79] (0.95 confidence)</td>
<td>0.998335</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d2/rounded&#95;array</td>
<td>11.546 ns [11.5183, 11.5739] (0.95 confidence)</td>
<td>11.5869 ns [11.5752, 11.6005] (0.95 confidence)</td>
<td>0.996466</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d2/rounded&#95;result</td>
<td>11.5012 ns [11.4856, 11.519] (0.95 confidence)</td>
<td>11.5828 ns [11.5731, 11.5936] (0.95 confidence)</td>
<td>0.992953</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d2/strict&#95;array</td>
<td>9.11253 ns [9.04553, 9.18023] (0.95 confidence)</td>
<td>8.52379 ns [8.50165, 8.54733] (0.95 confidence)</td>
<td>1.06907</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d2/strict&#95;result</td>
<td>8.22977 ns [8.18416, 8.27532] (0.95 confidence)</td>
<td>7.85116 ns [7.8348, 7.86841] (0.95 confidence)</td>
<td>1.04822</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d3/rounded&#95;array</td>
<td>17.6867 ns [17.6702, 17.7062] (0.95 confidence)</td>
<td>18.1119 ns [18.0716, 18.1577] (0.95 confidence)</td>
<td>0.976525</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d3/rounded&#95;result</td>
<td>17.7021 ns [17.6818, 17.7268] (0.95 confidence)</td>
<td>18.2218 ns [18.1865, 18.2577] (0.95 confidence)</td>
<td>0.971483</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d3/strict&#95;array</td>
<td>13.4689 ns [13.3551, 13.5845] (0.95 confidence)</td>
<td>13.5051 ns [13.417, 13.5945] (0.95 confidence)</td>
<td>0.997319</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d3/strict&#95;result</td>
<td>12.0232 ns [11.9534, 12.0977] (0.95 confidence)</td>
<td>11.9232 ns [11.8633, 11.9849] (0.95 confidence)</td>
<td>1.00839</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d4/rounded&#95;array</td>
<td>23.2367 ns [23.204, 23.2725] (0.95 confidence)</td>
<td>24.1501 ns [24.0908, 24.2157] (0.95 confidence)</td>
<td>0.962178</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d4/rounded&#95;result</td>
<td>23.2964 ns [23.267, 23.327] (0.95 confidence)</td>
<td>24.0477 ns [24.0057, 24.0918] (0.95 confidence)</td>
<td>0.968758</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d4/strict&#95;array</td>
<td>18.0589 ns [17.9558, 18.1608] (0.95 confidence)</td>
<td>18.2513 ns [18.1375, 18.3606] (0.95 confidence)</td>
<td>0.989455</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d4/strict&#95;result</td>
<td>15.7821 ns [15.6667, 15.8994] (0.95 confidence)</td>
<td>16.1853 ns [16.0937, 16.2767] (0.95 confidence)</td>
<td>0.975084</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d5/rounded&#95;array</td>
<td>30.0767 ns [29.9795, 30.1764] (0.95 confidence)</td>
<td>30.7689 ns [30.4598, 31.0776] (0.95 confidence)</td>
<td>0.977503</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d5/rounded&#95;result</td>
<td>29.4838 ns [29.4462, 29.5242] (0.95 confidence)</td>
<td>30.4636 ns [30.1921, 30.7391] (0.95 confidence)</td>
<td>0.967836</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d5/strict&#95;array</td>
<td>23.0423 ns [22.8865, 23.1938] (0.95 confidence)</td>
<td>23.9589 ns [23.6787, 24.2455] (0.95 confidence)</td>
<td>0.961741</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d5/strict&#95;result</td>
<td>20.1035 ns [19.9564, 20.2519] (0.95 confidence)</td>
<td>20.8087 ns [20.5743, 21.0435] (0.95 confidence)</td>
<td>0.966112</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d2/rounded&#95;array</td>
<td>177.292 ns [176.65, 178.034] (0.95 confidence)</td>
<td>174.087 ns [173.465, 174.723] (0.95 confidence)</td>
<td>1.01841</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d2/rounded&#95;result</td>
<td>178.231 ns [177.586, 178.919] (0.95 confidence)</td>
<td>172.691 ns [172.345, 173.07] (0.95 confidence)</td>
<td>1.03208</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d2/strict&#95;array</td>
<td>13.7905 ns [13.6805, 13.9011] (0.95 confidence)</td>
<td>14.2317 ns [14.1388, 14.3259] (0.95 confidence)</td>
<td>0.969001</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d2/strict&#95;result</td>
<td>13.5488 ns [13.4438, 13.6526] (0.95 confidence)</td>
<td>12.9119 ns [12.8229, 13.0033] (0.95 confidence)</td>
<td>1.04933</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d3/rounded&#95;array</td>
<td>180.549 ns [180.226, 180.896] (0.95 confidence)</td>
<td>181.142 ns [180.795, 181.534] (0.95 confidence)</td>
<td>0.996725</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d3/rounded&#95;result</td>
<td>179.3 ns [178.973, 179.694] (0.95 confidence)</td>
<td>180.737 ns [180.495, 180.986] (0.95 confidence)</td>
<td>0.992048</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d3/strict&#95;array</td>
<td>18.2514 ns [18.1423, 18.36] (0.95 confidence)</td>
<td>18.5127 ns [18.4177, 18.6094] (0.95 confidence)</td>
<td>0.985887</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d3/strict&#95;result</td>
<td>17.8323 ns [17.6893, 17.981] (0.95 confidence)</td>
<td>17.164 ns [16.9943, 17.3377] (0.95 confidence)</td>
<td>1.03894</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d4/rounded&#95;array</td>
<td>186.173 ns [185.863, 186.509] (0.95 confidence)</td>
<td>197.999 ns [196.013, 199.99] (0.95 confidence)</td>
<td>0.940276</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d4/rounded&#95;result</td>
<td>186.203 ns [185.882, 186.552] (0.95 confidence)</td>
<td>195.789 ns [193.737, 197.878] (0.95 confidence)</td>
<td>0.951044</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d4/strict&#95;array</td>
<td>23.2713 ns [23.1256, 23.4177] (0.95 confidence)</td>
<td>25.2074 ns [24.8747, 25.5458] (0.95 confidence)</td>
<td>0.923191</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d4/strict&#95;result</td>
<td>21.9624 ns [21.8207, 22.1117] (0.95 confidence)</td>
<td>22.9827 ns [22.674, 23.2963] (0.95 confidence)</td>
<td>0.955603</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d5/rounded&#95;array</td>
<td>191.876 ns [191.497, 192.3] (0.95 confidence)</td>
<td>192.67 ns [192.3, 193.056] (0.95 confidence)</td>
<td>0.995875</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d5/rounded&#95;result</td>
<td>194.789 ns [193.557, 196.278] (0.95 confidence)</td>
<td>201.61 ns [199.489, 203.808] (0.95 confidence)</td>
<td>0.966168</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d5/strict&#95;array</td>
<td>27.6162 ns [27.3789, 27.8611] (0.95 confidence)</td>
<td>27.4688 ns [27.3577, 27.5893] (0.95 confidence)</td>
<td>1.00537</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d5/strict&#95;result</td>
<td>26.4244 ns [26.2699, 26.5764] (0.95 confidence)</td>
<td>26.7777 ns [26.4996, 27.0599] (0.95 confidence)</td>
<td>0.986808</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d2/rounded&#95;array</td>
<td>177.591 ns [177.322, 177.872] (0.95 confidence)</td>
<td>174.876 ns [174.338, 175.422] (0.95 confidence)</td>
<td>1.01553</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d2/rounded&#95;result</td>
<td>176.092 ns [175.659, 176.582] (0.95 confidence)</td>
<td>171.637 ns [171.421, 171.867] (0.95 confidence)</td>
<td>1.02596</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d2/strict&#95;array</td>
<td>13.7282 ns [13.6615, 13.7978] (0.95 confidence)</td>
<td>13.4495 ns [13.3884, 13.5141] (0.95 confidence)</td>
<td>1.02072</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d2/strict&#95;result</td>
<td>13.301 ns [13.1992, 13.4082] (0.95 confidence)</td>
<td>12.6859 ns [12.6248, 12.7513] (0.95 confidence)</td>
<td>1.04849</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d3/rounded&#95;array</td>
<td>178.407 ns [178.073, 178.768] (0.95 confidence)</td>
<td>179.743 ns [179.43, 180.063] (0.95 confidence)</td>
<td>0.992564</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d3/rounded&#95;result</td>
<td>179.224 ns [178.885, 179.589] (0.95 confidence)</td>
<td>182.423 ns [181.71, 183.188] (0.95 confidence)</td>
<td>0.982464</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d3/strict&#95;array</td>
<td>17.5205 ns [17.429, 17.6137] (0.95 confidence)</td>
<td>17.7132 ns [17.6289, 17.8003] (0.95 confidence)</td>
<td>0.989125</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d3/strict&#95;result</td>
<td>16.7573 ns [16.6455, 16.8714] (0.95 confidence)</td>
<td>17.0238 ns [16.8852, 17.1642] (0.95 confidence)</td>
<td>0.984349</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d4/rounded&#95;array</td>
<td>190.47 ns [189.932, 191.054] (0.95 confidence)</td>
<td>195.983 ns [194.227, 197.774] (0.95 confidence)</td>
<td>0.97187</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d4/rounded&#95;result</td>
<td>187.867 ns [187.199, 188.625] (0.95 confidence)</td>
<td>198.572 ns [196.572, 200.559] (0.95 confidence)</td>
<td>0.946089</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d4/strict&#95;array</td>
<td>22.528 ns [22.4105, 22.6439] (0.95 confidence)</td>
<td>23.4926 ns [23.2293, 23.7503] (0.95 confidence)</td>
<td>0.958942</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d4/strict&#95;result</td>
<td>21.6347 ns [21.5038, 21.767] (0.95 confidence)</td>
<td>21.3357 ns [21.0633, 21.6114] (0.95 confidence)</td>
<td>1.01401</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d5/rounded&#95;array</td>
<td>193.581 ns [192.871, 194.462] (0.95 confidence)</td>
<td>190.51 ns [190.289, 190.756] (0.95 confidence)</td>
<td>1.01612</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d5/rounded&#95;result</td>
<td>191.062 ns [190.831, 191.31] (0.95 confidence)</td>
<td>190.912 ns [190.656, 191.19] (0.95 confidence)</td>
<td>1.00079</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d5/strict&#95;array</td>
<td>26.6534 ns [26.5202, 26.7905] (0.95 confidence)</td>
<td>26.7651 ns [26.5774, 26.9569] (0.95 confidence)</td>
<td>0.995826</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d5/strict&#95;result</td>
<td>25.8933 ns [25.6099, 26.1918] (0.95 confidence)</td>
<td>24.3202 ns [24.1817, 24.465] (0.95 confidence)</td>
<td>1.06468</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d2/rounded&#95;array</td>
<td>11.4075 ns [11.3921, 11.4244] (0.95 confidence)</td>
<td>11.6269 ns [11.6167, 11.6395] (0.95 confidence)</td>
<td>0.981131</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d2/rounded&#95;result</td>
<td>11.3782 ns [11.3675, 11.39] (0.95 confidence)</td>
<td>11.6166 ns [11.6076, 11.627] (0.95 confidence)</td>
<td>0.979482</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d2/strict&#95;array</td>
<td>57.2098 ns [56.9894, 57.441] (0.95 confidence)</td>
<td>56.1455 ns [55.9581, 56.3386] (0.95 confidence)</td>
<td>1.01896</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d2/strict&#95;result</td>
<td>12.2625 ns [12.2412, 12.2839] (0.95 confidence)</td>
<td>11.8609 ns [11.8347, 11.8882] (0.95 confidence)</td>
<td>1.03386</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d3/rounded&#95;array</td>
<td>17.7137 ns [17.6971, 17.7327] (0.95 confidence)</td>
<td>18.1337 ns [18.1054, 18.1643] (0.95 confidence)</td>
<td>0.976839</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d3/rounded&#95;result</td>
<td>17.7207 ns [17.7016, 17.742] (0.95 confidence)</td>
<td>18.1263 ns [18.0747, 18.2062] (0.95 confidence)</td>
<td>0.977621</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d3/strict&#95;array</td>
<td>60.8557 ns [60.7256, 60.9847] (0.95 confidence)</td>
<td>60.601 ns [60.3716, 60.8634] (0.95 confidence)</td>
<td>1.0042</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d3/strict&#95;result</td>
<td>15.9155 ns [15.8618, 15.9708] (0.95 confidence)</td>
<td>15.7749 ns [15.7219, 15.8292] (0.95 confidence)</td>
<td>1.00892</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d4/rounded&#95;array</td>
<td>23.2383 ns [23.1909, 23.2964] (0.95 confidence)</td>
<td>24.0379 ns [24.0158, 24.0642] (0.95 confidence)</td>
<td>0.966735</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d4/rounded&#95;result</td>
<td>23.3069 ns [23.2711, 23.3444] (0.95 confidence)</td>
<td>24.0115 ns [23.9719, 24.0564] (0.95 confidence)</td>
<td>0.970655</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d4/strict&#95;array</td>
<td>64.9929 ns [64.7985, 65.186] (0.95 confidence)</td>
<td>64.6415 ns [64.3153, 64.973] (0.95 confidence)</td>
<td>1.00543</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d4/strict&#95;result</td>
<td>19.5146 ns [19.4396, 19.5873] (0.95 confidence)</td>
<td>19.7315 ns [19.6794, 19.7837] (0.95 confidence)</td>
<td>0.989005</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d5/rounded&#95;array</td>
<td>29.5692 ns [29.5303, 29.6135] (0.95 confidence)</td>
<td>30.2521 ns [29.9505, 30.5581] (0.95 confidence)</td>
<td>0.977427</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d5/rounded&#95;result</td>
<td>29.725 ns [29.6683, 29.7844] (0.95 confidence)</td>
<td>30.3239 ns [30.0257, 30.6265] (0.95 confidence)</td>
<td>0.98025</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d5/strict&#95;array</td>
<td>70.3297 ns [70.0242, 70.6401] (0.95 confidence)</td>
<td>73.0567 ns [72.2568, 73.8403] (0.95 confidence)</td>
<td>0.962672</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d5/strict&#95;result</td>
<td>23.8905 ns [23.7882, 23.9927] (0.95 confidence)</td>
<td>23.9901 ns [23.7407, 24.2438] (0.95 confidence)</td>
<td>0.995849</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d2/rounded&#95;array</td>
<td>173.488 ns [172.878, 174.113] (0.95 confidence)</td>
<td>171.953 ns [171.402, 172.601] (0.95 confidence)</td>
<td>1.00892</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d2/rounded&#95;result</td>
<td>172.272 ns [171.724, 172.823] (0.95 confidence)</td>
<td>169.13 ns [168.742, 169.527] (0.95 confidence)</td>
<td>1.01858</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d2/strict&#95;array</td>
<td>180.827 ns [180.356, 181.308] (0.95 confidence)</td>
<td>181.304 ns [180.708, 181.916] (0.95 confidence)</td>
<td>0.997368</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d2/strict&#95;result</td>
<td>182.311 ns [181.733, 182.929] (0.95 confidence)</td>
<td>176.962 ns [176.586, 177.366] (0.95 confidence)</td>
<td>1.03023</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d3/rounded&#95;array</td>
<td>176.699 ns [176.168, 177.248] (0.95 confidence)</td>
<td>177.977 ns [177.491, 178.482] (0.95 confidence)</td>
<td>0.992814</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d3/rounded&#95;result</td>
<td>181.897 ns [181.349, 182.467] (0.95 confidence)</td>
<td>178.84 ns [178.305, 179.392] (0.95 confidence)</td>
<td>1.01709</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d3/strict&#95;array</td>
<td>181.655 ns [181.051, 182.273] (0.95 confidence)</td>
<td>186.203 ns [185.646, 186.776] (0.95 confidence)</td>
<td>0.975574</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d3/strict&#95;result</td>
<td>183.645 ns [183.085, 184.218] (0.95 confidence)</td>
<td>184.71 ns [184.052, 185.387] (0.95 confidence)</td>
<td>0.994238</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d4/rounded&#95;array</td>
<td>187.556 ns [186.962, 188.174] (0.95 confidence)</td>
<td>190.152 ns [188.276, 192.081] (0.95 confidence)</td>
<td>0.986347</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d4/rounded&#95;result</td>
<td>188.047 ns [187.348, 188.775] (0.95 confidence)</td>
<td>190.823 ns [188.954, 192.751] (0.95 confidence)</td>
<td>0.985451</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d4/strict&#95;array</td>
<td>196.624 ns [195.395, 198.072] (0.95 confidence)</td>
<td>200.318 ns [198.324, 202.318] (0.95 confidence)</td>
<td>0.981562</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d4/strict&#95;result</td>
<td>191.217 ns [190.6, 191.855] (0.95 confidence)</td>
<td>199.745 ns [197.842, 201.642] (0.95 confidence)</td>
<td>0.957303</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d5/rounded&#95;array</td>
<td>188.192 ns [187.368, 189.131] (0.95 confidence)</td>
<td>188.087 ns [187.673, 188.517] (0.95 confidence)</td>
<td>1.00055</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d5/rounded&#95;result</td>
<td>189.193 ns [188.071, 190.477] (0.95 confidence)</td>
<td>189.201 ns [188.807, 189.607] (0.95 confidence)</td>
<td>0.999957</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d5/strict&#95;array</td>
<td>197.731 ns [195.372, 200.216] (0.95 confidence)</td>
<td>191.36 ns [190.816, 191.931] (0.95 confidence)</td>
<td>1.03329</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d5/strict&#95;result</td>
<td>188.256 ns [187.116, 189.665] (0.95 confidence)</td>
<td>188.701 ns [188.21, 189.231] (0.95 confidence)</td>
<td>0.997641</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d2/rounded&#95;array</td>
<td>208.018 ns [207.469, 208.612] (0.95 confidence)</td>
<td>198.44 ns [197.848, 199.05] (0.95 confidence)</td>
<td>1.04827</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d2/rounded&#95;result</td>
<td>208.215 ns [207.593, 208.868] (0.95 confidence)</td>
<td>196.648 ns [196.328, 196.995] (0.95 confidence)</td>
<td>1.05882</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d2/strict&#95;array</td>
<td>13773.9 ns [13678.2, 13874.3] (0.95 confidence)</td>
<td>13259.7 ns [13234.9, 13284.3] (0.95 confidence)</td>
<td>1.03878</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d2/strict&#95;result</td>
<td>213.01 ns [212.319, 213.709] (0.95 confidence)</td>
<td>201.061 ns [200.474, 201.659] (0.95 confidence)</td>
<td>1.05943</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d3/rounded&#95;array</td>
<td>215.619 ns [215.028, 216.218] (0.95 confidence)</td>
<td>209.341 ns [208.643, 210.057] (0.95 confidence)</td>
<td>1.02999</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d3/rounded&#95;result</td>
<td>211.728 ns [211.308, 212.156] (0.95 confidence)</td>
<td>209.262 ns [208.822, 209.722] (0.95 confidence)</td>
<td>1.01179</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d3/strict&#95;array</td>
<td>13505.7 ns [13472.4, 13538] (0.95 confidence)</td>
<td>13370.9 ns [13335.4, 13408.7] (0.95 confidence)</td>
<td>1.01008</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d3/strict&#95;result</td>
<td>210.084 ns [209.728, 210.448] (0.95 confidence)</td>
<td>207.775 ns [207.315, 208.234] (0.95 confidence)</td>
<td>1.01111</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d4/rounded&#95;array</td>
<td>215.868 ns [215.012, 216.802] (0.95 confidence)</td>
<td>227.014 ns [224.622, 229.414] (0.95 confidence)</td>
<td>0.950902</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d4/rounded&#95;result</td>
<td>213.582 ns [213.047, 214.148] (0.95 confidence)</td>
<td>228.174 ns [226.412, 229.972] (0.95 confidence)</td>
<td>0.936047</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d4/strict&#95;array</td>
<td>13416.6 ns [13387.8, 13448.2] (0.95 confidence)</td>
<td>14045.3 ns [13930.6, 14162.8] (0.95 confidence)</td>
<td>0.955239</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d4/strict&#95;result</td>
<td>213.083 ns [212.281, 213.921] (0.95 confidence)</td>
<td>223.239 ns [221.143, 225.384] (0.95 confidence)</td>
<td>0.954506</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d5/rounded&#95;array</td>
<td>224.506 ns [223.896, 225.121] (0.95 confidence)</td>
<td>232.674 ns [230.309, 235.033] (0.95 confidence)</td>
<td>0.964892</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d5/rounded&#95;result</td>
<td>225.179 ns [224.645, 225.727] (0.95 confidence)</td>
<td>231.849 ns [229.538, 234.185] (0.95 confidence)</td>
<td>0.97123</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d5/strict&#95;array</td>
<td>13433.6 ns [13405, 13463.1] (0.95 confidence)</td>
<td>14234 ns [14105.9, 14363.9] (0.95 confidence)</td>
<td>0.943768</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d5/strict&#95;result</td>
<td>220.959 ns [220.355, 221.579] (0.95 confidence)</td>
<td>224.447 ns [222.147, 226.766] (0.95 confidence)</td>
<td>0.984459</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d2/rounded&#95;array</td>
<td>169.492 ns [168.651, 170.376] (0.95 confidence)</td>
<td>159.964 ns [159.468, 160.48] (0.95 confidence)</td>
<td>1.05956</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d2/rounded&#95;result</td>
<td>169.115 ns [168.255, 170.002] (0.95 confidence)</td>
<td>158.96 ns [158.501, 159.448] (0.95 confidence)</td>
<td>1.06388</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d2/strict&#95;array</td>
<td>3876.59 ns [3853.06, 3900.46] (0.95 confidence)</td>
<td>3915.89 ns [3900.04, 3932.64] (0.95 confidence)</td>
<td>0.989965</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d2/strict&#95;result</td>
<td>165.194 ns [164.571, 165.829] (0.95 confidence)</td>
<td>162.628 ns [162.094, 163.176] (0.95 confidence)</td>
<td>1.01578</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d3/rounded&#95;array</td>
<td>172.593 ns [172.069, 173.131] (0.95 confidence)</td>
<td>172.82 ns [172.249, 173.401] (0.95 confidence)</td>
<td>0.998684</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d3/rounded&#95;result</td>
<td>172.304 ns [171.796, 172.843] (0.95 confidence)</td>
<td>171.466 ns [170.709, 172.296] (0.95 confidence)</td>
<td>1.00489</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d3/strict&#95;array</td>
<td>3924.15 ns [3899.69, 3949.08] (0.95 confidence)</td>
<td>3859.86 ns [3839.49, 3881.43] (0.95 confidence)</td>
<td>1.01666</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d3/strict&#95;result</td>
<td>172.476 ns [172.036, 172.933] (0.95 confidence)</td>
<td>171.392 ns [170.802, 171.983] (0.95 confidence)</td>
<td>1.00632</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d4/rounded&#95;array</td>
<td>181.01 ns [179.565, 182.602] (0.95 confidence)</td>
<td>188.611 ns [186.964, 190.247] (0.95 confidence)</td>
<td>0.959703</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d4/rounded&#95;result</td>
<td>180.671 ns [179.583, 181.876] (0.95 confidence)</td>
<td>188.437 ns [186.801, 190.058] (0.95 confidence)</td>
<td>0.958789</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d4/strict&#95;array</td>
<td>3916.7 ns [3881.65, 3952.81] (0.95 confidence)</td>
<td>4075.19 ns [4042.24, 4108.92] (0.95 confidence)</td>
<td>0.961108</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d4/strict&#95;result</td>
<td>176.154 ns [175.601, 176.713] (0.95 confidence)</td>
<td>177.539 ns [176.456, 178.7] (0.95 confidence)</td>
<td>0.992199</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d5/rounded&#95;array</td>
<td>186.689 ns [186.093, 187.312] (0.95 confidence)</td>
<td>193.225 ns [191.249, 195.206] (0.95 confidence)</td>
<td>0.966176</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d5/rounded&#95;result</td>
<td>186.691 ns [186.066, 187.342] (0.95 confidence)</td>
<td>193.044 ns [191.047, 195.034] (0.95 confidence)</td>
<td>0.967092</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d5/strict&#95;array</td>
<td>3890.42 ns [3871.02, 3910.08] (0.95 confidence)</td>
<td>4184.68 ns [4142.65, 4226.57] (0.95 confidence)</td>
<td>0.929681</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d5/strict&#95;result</td>
<td>180.771 ns [180.049, 181.573] (0.95 confidence)</td>
<td>185.835 ns [184.08, 187.643] (0.95 confidence)</td>
<td>0.972751</td>
</tr>
<tr>
<td>exact&#95;d2/det</td>
<td>0.44951 ns [0.444411, 0.45511] (0.95 confidence)</td>
<td>0.449913 ns [0.444738, 0.455476] (0.95 confidence)</td>
<td>0.999106</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;direct</td>
<td>0.446875 ns [0.442821, 0.451308] (0.95 confidence)</td>
<td>0.45032 ns [0.445258, 0.455824] (0.95 confidence)</td>
<td>0.992351</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;direct&#95;with&#95;errbound</td>
<td>1.73886 ns [1.73633, 1.7415] (0.95 confidence)</td>
<td>1.737 ns [1.73471, 1.73935] (0.95 confidence)</td>
<td>1.00107</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;errbound</td>
<td>1.59048 ns [1.58882, 1.59235] (0.95 confidence)</td>
<td>1.59023 ns [1.58881, 1.59184] (0.95 confidence)</td>
<td>1.00016</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;exact</td>
<td>89.1715 ns [88.8973, 89.4526] (0.95 confidence)</td>
<td>87.9554 ns [87.7943, 88.1212] (0.95 confidence)</td>
<td>1.01383</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;exact&#95;f64&#95;result</td>
<td>67.9972 ns [67.875, 68.1313] (0.95 confidence)</td>
<td>67.6063 ns [67.4786, 67.7477] (0.95 confidence)</td>
<td>1.00578</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;exact&#95;rounded&#95;f64</td>
<td>70.782 ns [70.5011, 71.1414] (0.95 confidence)</td>
<td>70.1203 ns [70.0281, 70.2159] (0.95 confidence)</td>
<td>1.00944</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;sign&#95;exact</td>
<td>3.00747 ns [3.00286, 3.01312] (0.95 confidence)</td>
<td>2.97177 ns [2.96943, 2.97453] (0.95 confidence)</td>
<td>1.01201</td>
</tr>
<tr>
<td>exact&#95;d2/solve&#95;exact</td>
<td>2973.02 ns [2966.98, 2979.21] (0.95 confidence)</td>
<td>2978.28 ns [2974.67, 2982.06] (0.95 confidence)</td>
<td>0.998232</td>
</tr>
<tr>
<td>exact&#95;d2/solve&#95;exact&#95;f64&#95;result</td>
<td>3133.98 ns [3125.22, 3143.14] (0.95 confidence)</td>
<td>3120.3 ns [3115.08, 3125.48] (0.95 confidence)</td>
<td>1.00438</td>
</tr>
<tr>
<td>exact&#95;d2/solve&#95;exact&#95;rounded&#95;f64</td>
<td>3359.49 ns [3353.7, 3365.35] (0.95 confidence)</td>
<td>3285.06 ns [3271.99, 3300.62] (0.95 confidence)</td>
<td>1.02266</td>
</tr>
<tr>
<td>exact&#95;d3/det</td>
<td>0.853092 ns [0.850418, 0.855816] (0.95 confidence)</td>
<td>0.822565 ns [0.820839, 0.824521] (0.95 confidence)</td>
<td>1.03711</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;direct</td>
<td>0.859921 ns [0.855884, 0.864117] (0.95 confidence)</td>
<td>0.824769 ns [0.822555, 0.827149] (0.95 confidence)</td>
<td>1.04262</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;direct&#95;with&#95;errbound</td>
<td>3.35103 ns [3.34407, 3.35826] (0.95 confidence)</td>
<td>3.31536 ns [3.30865, 3.32282] (0.95 confidence)</td>
<td>1.01076</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;errbound</td>
<td>3.29662 ns [3.29389, 3.29977] (0.95 confidence)</td>
<td>3.29204 ns [3.2896, 3.29502] (0.95 confidence)</td>
<td>1.00139</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;exact</td>
<td>304.962 ns [304.512, 305.429] (0.95 confidence)</td>
<td>301.654 ns [301.249, 302.073] (0.95 confidence)</td>
<td>1.01097</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;exact&#95;f64&#95;result</td>
<td>280.084 ns [279.471, 280.716] (0.95 confidence)</td>
<td>278.143 ns [277.684, 278.601] (0.95 confidence)</td>
<td>1.00698</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;exact&#95;rounded&#95;f64</td>
<td>287.21 ns [286.572, 287.858] (0.95 confidence)</td>
<td>280.705 ns [280.335, 281.071] (0.95 confidence)</td>
<td>1.02317</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;sign&#95;exact</td>
<td>4.75368 ns [4.74891, 4.75894] (0.95 confidence)</td>
<td>4.71105 ns [4.70782, 4.71478] (0.95 confidence)</td>
<td>1.00905</td>
</tr>
<tr>
<td>exact&#95;d3/solve&#95;exact</td>
<td>9135.83 ns [9115.92, 9155.79] (0.95 confidence)</td>
<td>9112.45 ns [9092.74, 9132.57] (0.95 confidence)</td>
<td>1.00257</td>
</tr>
<tr>
<td>exact&#95;d3/solve&#95;exact&#95;f64&#95;result</td>
<td>9329.17 ns [9302.44, 9355.99] (0.95 confidence)</td>
<td>9289.14 ns [9262.45, 9314.16] (0.95 confidence)</td>
<td>1.00431</td>
</tr>
<tr>
<td>exact&#95;d3/solve&#95;exact&#95;rounded&#95;f64</td>
<td>9602.46 ns [9574.01, 9631.06] (0.95 confidence)</td>
<td>9595.92 ns [9580.49, 9611.46] (0.95 confidence)</td>
<td>1.00068</td>
</tr>
<tr>
<td>exact&#95;d4/det</td>
<td>3.24386 ns [3.23868, 3.24984] (0.95 confidence)</td>
<td>4.7951 ns [4.79155, 4.79918] (0.95 confidence)</td>
<td>0.676495</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;direct</td>
<td>2.37628 ns [2.3722, 2.3808] (0.95 confidence)</td>
<td>2.35076 ns [2.34914, 2.35279] (0.95 confidence)</td>
<td>1.01086</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;direct&#95;with&#95;errbound</td>
<td>6.88109 ns [6.81274, 6.97839] (0.95 confidence)</td>
<td>6.78338 ns [6.77778, 6.79018] (0.95 confidence)</td>
<td>1.0144</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;errbound</td>
<td>6.78969 ns [6.77137, 6.81025] (0.95 confidence)</td>
<td>6.77236 ns [6.75913, 6.78733] (0.95 confidence)</td>
<td>1.00256</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;exact</td>
<td>855.179 ns [852.716, 857.716] (0.95 confidence)</td>
<td>826.113 ns [825.326, 826.929] (0.95 confidence)</td>
<td>1.03518</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;exact&#95;f64&#95;result</td>
<td>814.218 ns [811.575, 816.964] (0.95 confidence)</td>
<td>794.431 ns [793.485, 795.371] (0.95 confidence)</td>
<td>1.02491</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;exact&#95;rounded&#95;f64</td>
<td>835.595 ns [831.981, 839.513] (0.95 confidence)</td>
<td>808.133 ns [806.26, 809.984] (0.95 confidence)</td>
<td>1.03398</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;sign&#95;exact</td>
<td>8.14874 ns [8.0888, 8.22362] (0.95 confidence)</td>
<td>7.92653 ns [7.91945, 7.93451] (0.95 confidence)</td>
<td>1.02803</td>
</tr>
<tr>
<td>exact&#95;d4/solve&#95;exact</td>
<td>18960 ns [18864.1, 19057.1] (0.95 confidence)</td>
<td>18981.7 ns [18940.2, 19022.7] (0.95 confidence)</td>
<td>0.998854</td>
</tr>
<tr>
<td>exact&#95;d4/solve&#95;exact&#95;f64&#95;result</td>
<td>18806.7 ns [18748.9, 18866.4] (0.95 confidence)</td>
<td>19136 ns [19103.8, 19169] (0.95 confidence)</td>
<td>0.982791</td>
</tr>
<tr>
<td>exact&#95;d4/solve&#95;exact&#95;rounded&#95;f64</td>
<td>19582.8 ns [19525, 19643.1] (0.95 confidence)</td>
<td>19726.1 ns [19694.3, 19758] (0.95 confidence)</td>
<td>0.992736</td>
</tr>
<tr>
<td>exact&#95;d5/det</td>
<td>28.5136 ns [28.4806, 28.5482] (0.95 confidence)</td>
<td>31.9676 ns [31.9271, 32.012] (0.95 confidence)</td>
<td>0.891954</td>
</tr>
<tr>
<td>exact&#95;d5/det&#95;exact</td>
<td>2868.6 ns [2865.2, 2872.13] (0.95 confidence)</td>
<td>2870.55 ns [2866.22, 2875.16] (0.95 confidence)</td>
<td>0.99932</td>
</tr>
<tr>
<td>exact&#95;d5/det&#95;exact&#95;f64&#95;result</td>
<td>2867.24 ns [2860.46, 2874.63] (0.95 confidence)</td>
<td>2844.85 ns [2841.11, 2848.86] (0.95 confidence)</td>
<td>1.00787</td>
</tr>
<tr>
<td>exact&#95;d5/det&#95;exact&#95;rounded&#95;f64</td>
<td>2892.82 ns [2884.91, 2902.04] (0.95 confidence)</td>
<td>2855.57 ns [2851.04, 2860.29] (0.95 confidence)</td>
<td>1.01304</td>
</tr>
<tr>
<td>exact&#95;d5/det&#95;sign&#95;exact</td>
<td>2835.49 ns [2829.25, 2842.51] (0.95 confidence)</td>
<td>2825.2 ns [2821.28, 2829.4] (0.95 confidence)</td>
<td>1.00364</td>
</tr>
<tr>
<td>exact&#95;d5/solve&#95;exact</td>
<td>30079.2 ns [29973, 30189.6] (0.95 confidence)</td>
<td>30536.6 ns [30459.5, 30615] (0.95 confidence)</td>
<td>0.98502</td>
</tr>
<tr>
<td>exact&#95;d5/solve&#95;exact&#95;f64&#95;result</td>
<td>30045.6 ns [29945.7, 30152.6] (0.95 confidence)</td>
<td>30864.1 ns [30774.7, 30965.6] (0.95 confidence)</td>
<td>0.973482</td>
</tr>
<tr>
<td>exact&#95;d5/solve&#95;exact&#95;rounded&#95;f64</td>
<td>30611.8 ns [30514, 30709.4] (0.95 confidence)</td>
<td>31217.3 ns [31111.4, 31318.8] (0.95 confidence)</td>
<td>0.980604</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;4x4/det&#95;exact</td>
<td>936.341 ns [935.145, 937.565] (0.95 confidence)</td>
<td>930.382 ns [928.489, 932.395] (0.95 confidence)</td>
<td>1.00641</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;4x4/det&#95;sign&#95;exact</td>
<td>7.97097 ns [7.96124, 7.98191] (0.95 confidence)</td>
<td>8.03548 ns [8.0229, 8.04888] (0.95 confidence)</td>
<td>0.991972</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;4x4/solve&#95;exact</td>
<td>15317.5 ns [15266.9, 15368.7] (0.95 confidence)</td>
<td>15586.2 ns [15550.8, 15621.4] (0.95 confidence)</td>
<td>0.982762</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;4x4/solve&#95;exact&#95;f64&#95;result</td>
<td>15487.7 ns [15457.8, 15518.6] (0.95 confidence)</td>
<td>15714.5 ns [15686.4, 15743] (0.95 confidence)</td>
<td>0.985572</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;4x4/solve&#95;exact&#95;rounded&#95;f64</td>
<td>15777 ns [15736.6, 15818.7] (0.95 confidence)</td>
<td>16117.1 ns [16080.2, 16153.4] (0.95 confidence)</td>
<td>0.978896</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;5x5/det&#95;exact</td>
<td>2951.43 ns [2947.07, 2955.73] (0.95 confidence)</td>
<td>2993.4 ns [2985.36, 3002.26] (0.95 confidence)</td>
<td>0.985978</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;5x5/det&#95;sign&#95;exact</td>
<td>2898.25 ns [2893.14, 2903.58] (0.95 confidence)</td>
<td>2967.32 ns [2961.63, 2973.11] (0.95 confidence)</td>
<td>0.976723</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;5x5/solve&#95;exact</td>
<td>24922 ns [24828.4, 25017.9] (0.95 confidence)</td>
<td>25324.8 ns [25272.8, 25376.5] (0.95 confidence)</td>
<td>0.984095</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;5x5/solve&#95;exact&#95;f64&#95;result</td>
<td>25044.2 ns [24956.3, 25137.5] (0.95 confidence)</td>
<td>25851 ns [25804.9, 25898.9] (0.95 confidence)</td>
<td>0.968791</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;5x5/solve&#95;exact&#95;rounded&#95;f64</td>
<td>25744 ns [25657, 25831.4] (0.95 confidence)</td>
<td>26783.6 ns [26701, 26870.9] (0.95 confidence)</td>
<td>0.961185</td>
</tr>
<tr>
<td>exact&#95;large&#95;entries&#95;3x3/det&#95;exact</td>
<td>318.753 ns [317.872, 319.644] (0.95 confidence)</td>
<td>323.313 ns [322.166, 324.465] (0.95 confidence)</td>
<td>0.985895</td>
</tr>
<tr>
<td>exact&#95;large&#95;entries&#95;3x3/det&#95;sign&#95;exact</td>
<td>287.529 ns [287.017, 288.062] (0.95 confidence)</td>
<td>295.479 ns [294.616, 296.374] (0.95 confidence)</td>
<td>0.973095</td>
</tr>
<tr>
<td>exact&#95;large&#95;entries&#95;3x3/solve&#95;exact</td>
<td>9108.53 ns [9089.25, 9129.95] (0.95 confidence)</td>
<td>9218.89 ns [9197.47, 9241.54] (0.95 confidence)</td>
<td>0.988029</td>
</tr>
<tr>
<td>exact&#95;large&#95;entries&#95;3x3/solve&#95;exact&#95;f64&#95;result</td>
<td>9287.5 ns [9267.91, 9308.22] (0.95 confidence)</td>
<td>9388.52 ns [9364.07, 9415.07] (0.95 confidence)</td>
<td>0.98924</td>
</tr>
<tr>
<td>exact&#95;large&#95;entries&#95;3x3/solve&#95;exact&#95;rounded&#95;f64</td>
<td>9526.02 ns [9511.54, 9541.14] (0.95 confidence)</td>
<td>9659.13 ns [9644.03, 9674.84] (0.95 confidence)</td>
<td>0.98622</td>
</tr>
<tr>
<td>exact&#95;near&#95;singular&#95;3x3/det&#95;exact</td>
<td>185.539 ns [185.313, 185.79] (0.95 confidence)</td>
<td>188.596 ns [188.291, 188.994] (0.95 confidence)</td>
<td>0.983794</td>
</tr>
<tr>
<td>exact&#95;near&#95;singular&#95;3x3/det&#95;sign&#95;exact</td>
<td>162.484 ns [162.24, 162.735] (0.95 confidence)</td>
<td>164.31 ns [164.06, 164.575] (0.95 confidence)</td>
<td>0.988885</td>
</tr>
<tr>
<td>exact&#95;near&#95;singular&#95;3x3/solve&#95;exact</td>
<td>968.752 ns [967.713, 969.828] (0.95 confidence)</td>
<td>1001.17 ns [998.684, 1003.69] (0.95 confidence)</td>
<td>0.967621</td>
</tr>
<tr>
<td>exact&#95;near&#95;singular&#95;3x3/solve&#95;exact&#95;f64&#95;result</td>
<td>996.469 ns [994.819, 998.095] (0.95 confidence)</td>
<td>1029.62 ns [1026.81, 1033.07] (0.95 confidence)</td>
<td>0.967804</td>
</tr>
<tr>
<td>exact&#95;near&#95;singular&#95;3x3/solve&#95;exact&#95;rounded&#95;f64</td>
<td>993.499 ns [989.777, 997.518] (0.95 confidence)</td>
<td>1017.3 ns [1014.98, 1019.65] (0.95 confidence)</td>
<td>0.976607</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d2/det&#95;exact</td>
<td>3561.94 ns [3553.79, 3570.17] (0.95 confidence)</td>
<td>3486.21 ns [3480.15, 3492.59] (0.95 confidence)</td>
<td>1.02172</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d2/det&#95;sign&#95;exact</td>
<td>152.364 ns [151.811, 152.931] (0.95 confidence)</td>
<td>143.164 ns [142.994, 143.345] (0.95 confidence)</td>
<td>1.06426</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d2/solve&#95;exact</td>
<td>37209.7 ns [37124.2, 37292.4] (0.95 confidence)</td>
<td>35952.6 ns [35881.7, 36028.4] (0.95 confidence)</td>
<td>1.03496</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d2/solve&#95;exact&#95;f64&#95;result</td>
<td>38507.8 ns [38319.4, 38704] (0.95 confidence)</td>
<td>36459.7 ns [36410.2, 36510.8] (0.95 confidence)</td>
<td>1.05617</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d2/solve&#95;exact&#95;rounded&#95;f64</td>
<td>37680 ns [37583.8, 37778.1] (0.95 confidence)</td>
<td>36989.7 ns [36927.5, 37054.5] (0.95 confidence)</td>
<td>1.01866</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d3/det&#95;exact</td>
<td>9614.56 ns [9582.24, 9650.68] (0.95 confidence)</td>
<td>9318.53 ns [9308.65, 9328.57] (0.95 confidence)</td>
<td>1.03177</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d3/det&#95;sign&#95;exact</td>
<td>230.337 ns [229.918, 230.784] (0.95 confidence)</td>
<td>225.604 ns [225.364, 225.882] (0.95 confidence)</td>
<td>1.02098</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d3/solve&#95;exact</td>
<td>84537.5 ns [84292.7, 84798.3] (0.95 confidence)</td>
<td>81192.8 ns [81050.3, 81339.7] (0.95 confidence)</td>
<td>1.04119</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d3/solve&#95;exact&#95;f64&#95;result</td>
<td>85110.2 ns [84883.3, 85351.6] (0.95 confidence)</td>
<td>81446.4 ns [81322.1, 81574.4] (0.95 confidence)</td>
<td>1.04498</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d3/solve&#95;exact&#95;rounded&#95;f64</td>
<td>86529.8 ns [86113.4, 86974] (0.95 confidence)</td>
<td>81793.2 ns [81682.2, 81906] (0.95 confidence)</td>
<td>1.05791</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d4/det&#95;exact</td>
<td>23657.4 ns [23596.6, 23719.5] (0.95 confidence)</td>
<td>22354 ns [22300.8, 22409.4] (0.95 confidence)</td>
<td>1.05831</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d4/det&#95;sign&#95;exact</td>
<td>463.576 ns [462.656, 464.549] (0.95 confidence)</td>
<td>443.515 ns [442.98, 444.102] (0.95 confidence)</td>
<td>1.04523</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d4/solve&#95;exact</td>
<td>150777 ns [150383, 151195] (0.95 confidence)</td>
<td>144402 ns [144136, 144688] (0.95 confidence)</td>
<td>1.04415</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d4/solve&#95;exact&#95;f64&#95;result</td>
<td>148586 ns [148195, 148997] (0.95 confidence)</td>
<td>145791 ns [145434, 146195] (0.95 confidence)</td>
<td>1.01917</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d4/solve&#95;exact&#95;rounded&#95;f64</td>
<td>150341 ns [149557, 151351] (0.95 confidence)</td>
<td>146087 ns [145839, 146344] (0.95 confidence)</td>
<td>1.02913</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d5/det&#95;exact</td>
<td>58262.6 ns [58015.7, 58522.7] (0.95 confidence)</td>
<td>56156.4 ns [55959.8, 56379.9] (0.95 confidence)</td>
<td>1.03751</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d5/det&#95;sign&#95;exact</td>
<td>57299.6 ns [56961.2, 57665.1] (0.95 confidence)</td>
<td>55085 ns [54971.1, 55203.3] (0.95 confidence)</td>
<td>1.0402</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d5/solve&#95;exact</td>
<td>244302 ns [243541, 245121] (0.95 confidence)</td>
<td>239192 ns [238497, 240000] (0.95 confidence)</td>
<td>1.02136</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d5/solve&#95;exact&#95;f64&#95;result</td>
<td>242787 ns [242242, 243337] (0.95 confidence)</td>
<td>239148 ns [238776, 239536] (0.95 confidence)</td>
<td>1.01522</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d5/solve&#95;exact&#95;rounded&#95;f64</td>
<td>245668 ns [245069, 246319] (0.95 confidence)</td>
<td>240608 ns [240179, 241067] (0.95 confidence)</td>
<td>1.02103</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/matrix&#95;access</td>
<td>0.779782 ns [0.777862, 0.781628] (0.95 confidence)</td>
<td>0.799685 ns [0.791835, 0.808644] (0.95 confidence)</td>
<td>0.975111</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/matrix&#95;from&#95;fn</td>
<td>135.869 ns [135.629, 136.127] (0.95 confidence)</td>
<td>136.269 ns [136.004, 136.545] (0.95 confidence)</td>
<td>0.997071</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/matrix&#95;from&#95;rows</td>
<td>138.603 ns [138.337, 138.848] (0.95 confidence)</td>
<td>140.025 ns [139.739, 140.315] (0.95 confidence)</td>
<td>0.989842</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/matrix&#95;into&#95;rows</td>
<td>66.9787 ns [66.8265, 67.1438] (0.95 confidence)</td>
<td>67.7196 ns [67.5688, 67.8654] (0.95 confidence)</td>
<td>0.98906</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/matrix&#95;zero</td>
<td>2.73949 ns [2.73674, 2.74236] (0.95 confidence)</td>
<td>2.7511 ns [2.7466, 2.75586] (0.95 confidence)</td>
<td>0.995779</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/vector&#95;access</td>
<td>0.785975 ns [0.781748, 0.7907] (0.95 confidence)</td>
<td>0.79091 ns [0.788075, 0.793691] (0.95 confidence)</td>
<td>0.99376</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/vector&#95;from&#95;fn</td>
<td>37.4018 ns [37.3518, 37.4514] (0.95 confidence)</td>
<td>37.7872 ns [37.7163, 37.8636] (0.95 confidence)</td>
<td>0.989799</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/vector&#95;into&#95;array</td>
<td>15.9061 ns [15.877, 15.9364] (0.95 confidence)</td>
<td>16.1116 ns [16.0837, 16.1407] (0.95 confidence)</td>
<td>0.987245</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/vector&#95;new</td>
<td>36.8934 ns [36.8417, 36.947] (0.95 confidence)</td>
<td>37.1092 ns [37.0502, 37.1702] (0.95 confidence)</td>
<td>0.994186</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/vector&#95;zero</td>
<td>1.37306 ns [1.37092, 1.37556] (0.95 confidence)</td>
<td>1.38003 ns [1.3778, 1.38239] (0.95 confidence)</td>
<td>0.994943</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/matrix&#95;access</td>
<td>0.778246 ns [0.776074, 0.780348] (0.95 confidence)</td>
<td>0.788873 ns [0.786631, 0.7913] (0.95 confidence)</td>
<td>0.986529</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/matrix&#95;from&#95;fn</td>
<td>232.893 ns [232.591, 233.2] (0.95 confidence)</td>
<td>235.663 ns [235.102, 236.25] (0.95 confidence)</td>
<td>0.988244</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/matrix&#95;from&#95;rows</td>
<td>270.903 ns [270.593, 271.189] (0.95 confidence)</td>
<td>273.887 ns [273.426, 274.335] (0.95 confidence)</td>
<td>0.989102</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/matrix&#95;into&#95;rows</td>
<td>125.475 ns [125.309, 125.645] (0.95 confidence)</td>
<td>126.445 ns [126.208, 126.713] (0.95 confidence)</td>
<td>0.992329</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/matrix&#95;zero</td>
<td>6.3446 ns [6.33487, 6.35415] (0.95 confidence)</td>
<td>6.42486 ns [6.41266, 6.43772] (0.95 confidence)</td>
<td>0.987509</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/vector&#95;access</td>
<td>0.781559 ns [0.779307, 0.783766] (0.95 confidence)</td>
<td>0.793679 ns [0.789306, 0.79969] (0.95 confidence)</td>
<td>0.984729</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/vector&#95;from&#95;fn</td>
<td>52.8438 ns [52.7665, 52.924] (0.95 confidence)</td>
<td>53.2278 ns [53.1362, 53.3277] (0.95 confidence)</td>
<td>0.992785</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/vector&#95;into&#95;array</td>
<td>22.9249 ns [22.8942, 22.9555] (0.95 confidence)</td>
<td>23.1117 ns [23.0594, 23.1649] (0.95 confidence)</td>
<td>0.991919</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/vector&#95;new</td>
<td>52.3599 ns [52.2944, 52.4274] (0.95 confidence)</td>
<td>52.7564 ns [52.6621, 52.8535] (0.95 confidence)</td>
<td>0.992484</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/vector&#95;zero</td>
<td>2.04814 ns [2.04555, 2.05068] (0.95 confidence)</td>
<td>2.05957 ns [2.05658, 2.06273] (0.95 confidence)</td>
<td>0.994449</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/matrix&#95;access</td>
<td>0.780009 ns [0.777591, 0.782334] (0.95 confidence)</td>
<td>0.786287 ns [0.783064, 0.790014] (0.95 confidence)</td>
<td>0.992016</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/matrix&#95;from&#95;fn</td>
<td>289.668 ns [289.07, 290.254] (0.95 confidence)</td>
<td>293.417 ns [292.824, 294.023] (0.95 confidence)</td>
<td>0.987222</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/matrix&#95;from&#95;rows</td>
<td>351.945 ns [351.342, 352.576] (0.95 confidence)</td>
<td>353.819 ns [353.155, 354.506] (0.95 confidence)</td>
<td>0.994704</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/matrix&#95;into&#95;rows</td>
<td>133.391 ns [133.155, 133.624] (0.95 confidence)</td>
<td>134.221 ns [134.024, 134.422] (0.95 confidence)</td>
<td>0.993819</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/matrix&#95;zero</td>
<td>15.2281 ns [15.1419, 15.3181] (0.95 confidence)</td>
<td>15.6312 ns [15.5773, 15.6859] (0.95 confidence)</td>
<td>0.974209</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/vector&#95;access</td>
<td>0.780215 ns [0.777666, 0.782868] (0.95 confidence)</td>
<td>0.790809 ns [0.787536, 0.79429] (0.95 confidence)</td>
<td>0.986603</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/vector&#95;from&#95;fn</td>
<td>67.521 ns [67.4387, 67.604] (0.95 confidence)</td>
<td>68.1967 ns [68.0536, 68.3558] (0.95 confidence)</td>
<td>0.990093</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/vector&#95;into&#95;array</td>
<td>30.4226 ns [30.3892, 30.4552] (0.95 confidence)</td>
<td>30.6536 ns [30.6057, 30.7006] (0.95 confidence)</td>
<td>0.992462</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/vector&#95;new</td>
<td>66.5071 ns [66.4149, 66.5967] (0.95 confidence)</td>
<td>67.0376 ns [66.9181, 67.154] (0.95 confidence)</td>
<td>0.992087</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/vector&#95;zero</td>
<td>3.33796 ns [3.32621, 3.35008] (0.95 confidence)</td>
<td>3.35923 ns [3.34893, 3.37036] (0.95 confidence)</td>
<td>0.993668</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/matrix&#95;access</td>
<td>0.782256 ns [0.779143, 0.785837] (0.95 confidence)</td>
<td>0.804639 ns [0.801482, 0.807972] (0.95 confidence)</td>
<td>0.972183</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/matrix&#95;from&#95;fn</td>
<td>440.967 ns [439.891, 442.049] (0.95 confidence)</td>
<td>440.972 ns [439.577, 442.383] (0.95 confidence)</td>
<td>0.999988</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/matrix&#95;from&#95;rows</td>
<td>631.351 ns [630.179, 632.632] (0.95 confidence)</td>
<td>633.659 ns [632.28, 635.104] (0.95 confidence)</td>
<td>0.996358</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/matrix&#95;into&#95;rows</td>
<td>305.936 ns [305.158, 306.72] (0.95 confidence)</td>
<td>307.319 ns [306.604, 308.052] (0.95 confidence)</td>
<td>0.995499</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/matrix&#95;zero</td>
<td>22.3066 ns [22.2743, 22.339] (0.95 confidence)</td>
<td>22.8939 ns [22.843, 22.9471] (0.95 confidence)</td>
<td>0.974345</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/vector&#95;access</td>
<td>0.777262 ns [0.774688, 0.779722] (0.95 confidence)</td>
<td>0.817594 ns [0.807919, 0.828795] (0.95 confidence)</td>
<td>0.95067</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/vector&#95;from&#95;fn</td>
<td>82.8874 ns [82.7949, 82.9841] (0.95 confidence)</td>
<td>83.6639 ns [83.5359, 83.7968] (0.95 confidence)</td>
<td>0.990718</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/vector&#95;into&#95;array</td>
<td>38.8144 ns [38.7452, 38.8844] (0.95 confidence)</td>
<td>39.0394 ns [38.9622, 39.1183] (0.95 confidence)</td>
<td>0.994235</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/vector&#95;new</td>
<td>83.0896 ns [82.9861, 83.1986] (0.95 confidence)</td>
<td>83.4279 ns [83.3018, 83.557] (0.95 confidence)</td>
<td>0.995944</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/vector&#95;zero</td>
<td>4.32568 ns [4.30166, 4.34923] (0.95 confidence)</td>
<td>4.31444 ns [4.29482, 4.33433] (0.95 confidence)</td>
<td>1.0026</td>
</tr>
<tr>
<td>public&#95;exact&#95;scalar/rounded</td>
<td>5.58609 ns [5.58042, 5.59234] (0.95 confidence)</td>
<td>5.61703 ns [5.60843, 5.62617] (0.95 confidence)</td>
<td>0.99449</td>
</tr>
<tr>
<td>public&#95;exact&#95;scalar/strict</td>
<td>4.32772 ns [4.30696, 4.34915] (0.95 confidence)</td>
<td>4.53163 ns [4.50748, 4.56147] (0.95 confidence)</td>
<td>0.955002</td>
</tr>
<tr>
<td>rational&#95;input&#95;d2/det&#95;row&#95;cleared&#95;bareiss</td>
<td>179.664 ns [178.691, 180.825] (0.95 confidence)</td>
<td>179.598 ns [179.295, 179.906] (0.95 confidence)</td>
<td>1.00037</td>
</tr>
<tr>
<td>rational&#95;input&#95;d2/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>77.3566 ns [77.1969, 77.526] (0.95 confidence)</td>
<td>77.8297 ns [77.6866, 77.9783] (0.95 confidence)</td>
<td>0.993922</td>
</tr>
<tr>
<td>rational&#95;input&#95;d2/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>624.761 ns [623.846, 625.74] (0.95 confidence)</td>
<td>625.603 ns [624.593, 626.581] (0.95 confidence)</td>
<td>0.998653</td>
</tr>
<tr>
<td>rational&#95;input&#95;d3/det&#95;row&#95;cleared&#95;bareiss</td>
<td>596.656 ns [595.673, 597.675] (0.95 confidence)</td>
<td>598.862 ns [597.875, 599.834] (0.95 confidence)</td>
<td>0.996317</td>
</tr>
<tr>
<td>rational&#95;input&#95;d3/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>390.656 ns [389.5, 391.899] (0.95 confidence)</td>
<td>396.812 ns [395.712, 397.947] (0.95 confidence)</td>
<td>0.984486</td>
</tr>
<tr>
<td>rational&#95;input&#95;d3/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>1262.53 ns [1260.86, 1264.35] (0.95 confidence)</td>
<td>1268.81 ns [1266.89, 1270.76] (0.95 confidence)</td>
<td>0.995045</td>
</tr>
<tr>
<td>rational&#95;input&#95;d4/det&#95;row&#95;cleared&#95;bareiss</td>
<td>1234.22 ns [1229.56, 1239.17] (0.95 confidence)</td>
<td>1242.23 ns [1240.36, 1244.2] (0.95 confidence)</td>
<td>0.993549</td>
</tr>
<tr>
<td>rational&#95;input&#95;d4/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>821.58 ns [820.064, 823.125] (0.95 confidence)</td>
<td>830.17 ns [828.72, 831.635] (0.95 confidence)</td>
<td>0.989653</td>
</tr>
<tr>
<td>rational&#95;input&#95;d4/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>2242.61 ns [2228.09, 2261.47] (0.95 confidence)</td>
<td>2277.63 ns [2263.61, 2295.79] (0.95 confidence)</td>
<td>0.984623</td>
</tr>
<tr>
<td>rational&#95;input&#95;d5/det&#95;row&#95;cleared&#95;bareiss</td>
<td>2105.24 ns [2101.55, 2109.05] (0.95 confidence)</td>
<td>2099.79 ns [2094.53, 2105.41] (0.95 confidence)</td>
<td>1.0026</td>
</tr>
<tr>
<td>rational&#95;input&#95;d5/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>1576.24 ns [1573.49, 1579.12] (0.95 confidence)</td>
<td>1585.99 ns [1581.24, 1590.99] (0.95 confidence)</td>
<td>0.993857</td>
</tr>
<tr>
<td>rational&#95;input&#95;d5/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>3368.46 ns [3364.88, 3372.08] (0.95 confidence)</td>
<td>3385.96 ns [3379.83, 3392.58] (0.95 confidence)</td>
<td>0.994834</td>
</tr>
<tr>
<td>rational&#95;input&#95;d6/det&#95;row&#95;cleared&#95;bareiss</td>
<td>3520.35 ns [3512.66, 3528.93] (0.95 confidence)</td>
<td>3537.4 ns [3529.73, 3545.45] (0.95 confidence)</td>
<td>0.99518</td>
</tr>
<tr>
<td>rational&#95;input&#95;d6/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>2858.05 ns [2851.63, 2864.66] (0.95 confidence)</td>
<td>2882.86 ns [2876.87, 2889.19] (0.95 confidence)</td>
<td>0.991392</td>
</tr>
<tr>
<td>rational&#95;input&#95;d6/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>5743 ns [5735.95, 5750.45] (0.95 confidence)</td>
<td>5835.65 ns [5823.74, 5847.74] (0.95 confidence)</td>
<td>0.984124</td>
</tr>
<tr>
<td>rational&#95;input&#95;d7/det&#95;row&#95;cleared&#95;bareiss</td>
<td>5475.94 ns [5469.73, 5482.46] (0.95 confidence)</td>
<td>5554.33 ns [5544.05, 5565.94] (0.95 confidence)</td>
<td>0.985887</td>
</tr>
<tr>
<td>rational&#95;input&#95;d7/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>4689.11 ns [4681.61, 4697.13] (0.95 confidence)</td>
<td>4732.26 ns [4724.03, 4740.8] (0.95 confidence)</td>
<td>0.990883</td>
</tr>
<tr>
<td>rational&#95;input&#95;d7/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>9901.04 ns [9888.44, 9913.82] (0.95 confidence)</td>
<td>10084.2 ns [10067.6, 10101.1] (0.95 confidence)</td>
<td>0.981837</td>
</tr>
<tr>
<td>rational&#95;input&#95;d8/det&#95;row&#95;cleared&#95;bareiss</td>
<td>8157.1 ns [8142.05, 8172.94] (0.95 confidence)</td>
<td>8181.86 ns [8173.69, 8190.31] (0.95 confidence)</td>
<td>0.996974</td>
</tr>
<tr>
<td>rational&#95;input&#95;d8/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>7333.18 ns [7324.68, 7341.94] (0.95 confidence)</td>
<td>7397.08 ns [7384.12, 7411.29] (0.95 confidence)</td>
<td>0.991362</td>
</tr>
<tr>
<td>rational&#95;input&#95;d8/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>15136 ns [15112.4, 15161.6] (0.95 confidence)</td>
<td>15391.7 ns [15353.6, 15431.1] (0.95 confidence)</td>
<td>0.983388</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d2/det&#95;row&#95;cleared&#95;bareiss</td>
<td>42297.1 ns [42213.5, 42385] (0.95 confidence)</td>
<td>43469.2 ns [43370.4, 43572.3] (0.95 confidence)</td>
<td>0.973038</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d2/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>954.341 ns [953.003, 955.772] (0.95 confidence)</td>
<td>984.917 ns [982.456, 987.476] (0.95 confidence)</td>
<td>0.968956</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d2/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>3037.31 ns [3034.13, 3040.67] (0.95 confidence)</td>
<td>3263.59 ns [3231.99, 3299.4] (0.95 confidence)</td>
<td>0.930666</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d3/det&#95;row&#95;cleared&#95;bareiss</td>
<td>69199.7 ns [69102, 69304.4] (0.95 confidence)</td>
<td>70837.3 ns [70699.1, 70977.6] (0.95 confidence)</td>
<td>0.976882</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d3/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>4198.22 ns [4193.33, 4203.46] (0.95 confidence)</td>
<td>4275.19 ns [4268.26, 4282.42] (0.95 confidence)</td>
<td>0.981996</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d3/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>11432.8 ns [11420.3, 11446.2] (0.95 confidence)</td>
<td>11608.7 ns [11593.5, 11625.2] (0.95 confidence)</td>
<td>0.984842</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d4/det&#95;row&#95;cleared&#95;bareiss</td>
<td>110157 ns [110045, 110277] (0.95 confidence)</td>
<td>112028 ns [111915, 112147] (0.95 confidence)</td>
<td>0.9833</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d4/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>9398.59 ns [9387.9, 9409.6] (0.95 confidence)</td>
<td>9596.2 ns [9575.98, 9617.47] (0.95 confidence)</td>
<td>0.979408</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d4/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>33669.8 ns [33637.2, 33703.2] (0.95 confidence)</td>
<td>34553.2 ns [34499.2, 34607.9] (0.95 confidence)</td>
<td>0.974434</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d5/det&#95;row&#95;cleared&#95;bareiss</td>
<td>187705 ns [187389, 188063] (0.95 confidence)</td>
<td>188664 ns [188516, 188825] (0.95 confidence)</td>
<td>0.994919</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d5/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>43129.3 ns [43087.9, 43173.4] (0.95 confidence)</td>
<td>43792.3 ns [43744.5, 43846.1] (0.95 confidence)</td>
<td>0.984861</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d5/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>70629.7 ns [70485.7, 70807.9] (0.95 confidence)</td>
<td>70778.3 ns [70726, 70837.8] (0.95 confidence)</td>
<td>0.997901</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d6/det&#95;row&#95;cleared&#95;bareiss</td>
<td>300980 ns [300707, 301286] (0.95 confidence)</td>
<td>308645 ns [308313, 309057] (0.95 confidence)</td>
<td>0.975163</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d6/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>104660 ns [104560, 104774] (0.95 confidence)</td>
<td>106753 ns [106577, 106939] (0.95 confidence)</td>
<td>0.980398</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d6/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>162147 ns [161939, 162358] (0.95 confidence)</td>
<td>166510 ns [166347, 166679] (0.95 confidence)</td>
<td>0.973798</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d7/det&#95;row&#95;cleared&#95;bareiss</td>
<td>503230 ns [502527, 504039] (0.95 confidence)</td>
<td>515882 ns [515293, 516538] (0.95 confidence)</td>
<td>0.975474</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d7/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>234003 ns [233642, 234434] (0.95 confidence)</td>
<td>239917 ns [239700, 240146] (0.95 confidence)</td>
<td>0.97535</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d7/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>340276 ns [339855, 340743] (0.95 confidence)</td>
<td>351506 ns [350946, 352104] (0.95 confidence)</td>
<td>0.968052</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d8/det&#95;row&#95;cleared&#95;bareiss</td>
<td>824027 ns [823317, 824846] (0.95 confidence)</td>
<td>836306 ns [835394, 837301] (0.95 confidence)</td>
<td>0.985318</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d8/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>469097 ns [468673, 469576] (0.95 confidence)</td>
<td>476956 ns [476429, 477516] (0.95 confidence)</td>
<td>0.983523</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d8/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>652935 ns [651969, 654038] (0.95 confidence)</td>
<td>668570 ns [667163, 670032] (0.95 confidence)</td>
<td>0.976614</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d2/det&#95;row&#95;cleared&#95;bareiss</td>
<td>9062.72 ns [9032.42, 9096.22] (0.95 confidence)</td>
<td>9251.57 ns [9224.26, 9280.18] (0.95 confidence)</td>
<td>0.979588</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d2/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>698.7 ns [690.908, 706.655] (0.95 confidence)</td>
<td>660.946 ns [659.481, 662.459] (0.95 confidence)</td>
<td>1.05712</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d2/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>2066.36 ns [2062.22, 2070.74] (0.95 confidence)</td>
<td>2124.91 ns [2116.63, 2133.49] (0.95 confidence)</td>
<td>0.972444</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d3/det&#95;row&#95;cleared&#95;bareiss</td>
<td>13443.9 ns [13413.5, 13476.6] (0.95 confidence)</td>
<td>13701.1 ns [13673.2, 13728.9] (0.95 confidence)</td>
<td>0.981226</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d3/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>2765.41 ns [2758.24, 2773.18] (0.95 confidence)</td>
<td>2798.22 ns [2791.24, 2805.4] (0.95 confidence)</td>
<td>0.988275</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d3/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>5782.83 ns [5770.65, 5795.67] (0.95 confidence)</td>
<td>5908.3 ns [5891.97, 5925.5] (0.95 confidence)</td>
<td>0.978764</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d4/det&#95;row&#95;cleared&#95;bareiss</td>
<td>21694.8 ns [21658.6, 21733.6] (0.95 confidence)</td>
<td>22230.3 ns [22196, 22265.4] (0.95 confidence)</td>
<td>0.975913</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d4/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>5914.91 ns [5900.48, 5930.41] (0.95 confidence)</td>
<td>6091.92 ns [6070.13, 6114.42] (0.95 confidence)</td>
<td>0.970943</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d4/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>12753.8 ns [12713.2, 12796.8] (0.95 confidence)</td>
<td>13020.9 ns [12989.2, 13053.8] (0.95 confidence)</td>
<td>0.979483</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d5/det&#95;row&#95;cleared&#95;bareiss</td>
<td>33634.5 ns [33585.8, 33690.5] (0.95 confidence)</td>
<td>34404 ns [34357.8, 34457.1] (0.95 confidence)</td>
<td>0.977633</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d5/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>13467.9 ns [13434.9, 13503.2] (0.95 confidence)</td>
<td>13510.4 ns [13480.5, 13541.4] (0.95 confidence)</td>
<td>0.996854</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d5/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>20950.8 ns [20901.4, 21003] (0.95 confidence)</td>
<td>21177.6 ns [21128.3, 21229.3] (0.95 confidence)</td>
<td>0.98929</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d6/det&#95;row&#95;cleared&#95;bareiss</td>
<td>52562 ns [52444.5, 52696.5] (0.95 confidence)</td>
<td>53074.3 ns [52972.8, 53204.6] (0.95 confidence)</td>
<td>0.990349</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d6/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>27461.8 ns [27244.6, 27686.9] (0.95 confidence)</td>
<td>26518.5 ns [26485.1, 26553.8] (0.95 confidence)</td>
<td>1.03557</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d6/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>39628 ns [39559.4, 39703.1] (0.95 confidence)</td>
<td>40258.8 ns [40166.6, 40352.2] (0.95 confidence)</td>
<td>0.984331</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d7/det&#95;row&#95;cleared&#95;bareiss</td>
<td>83577.9 ns [83424.2, 83755.9] (0.95 confidence)</td>
<td>84535.6 ns [84446.3, 84629.9] (0.95 confidence)</td>
<td>0.988671</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d7/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>50215.3 ns [50153.4, 50281.1] (0.95 confidence)</td>
<td>51264.5 ns [51131.8, 51405.6] (0.95 confidence)</td>
<td>0.979534</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d7/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>72237.9 ns [72098.3, 72395.5] (0.95 confidence)</td>
<td>73235.4 ns [73168, 73308.6] (0.95 confidence)</td>
<td>0.986379</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d8/det&#95;row&#95;cleared&#95;bareiss</td>
<td>131582 ns [131453, 131727] (0.95 confidence)</td>
<td>134356 ns [134093, 134633] (0.95 confidence)</td>
<td>0.979354</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d8/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>90804 ns [90560.9, 91075.5] (0.95 confidence)</td>
<td>92129.6 ns [91993, 92269] (0.95 confidence)</td>
<td>0.985612</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d8/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>127272 ns [127152, 127401] (0.95 confidence)</td>
<td>130728 ns [130380, 131085] (0.95 confidence)</td>
<td>0.97356</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d2/adapter</td>
<td>1085.21 ns [1070.51, 1100.9] (0.95 confidence)</td>
<td>1078.16 ns [1063.45, 1093.84] (0.95 confidence)</td>
<td>1.00653</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d2/construct</td>
<td>292.748 ns [292.038, 293.464] (0.95 confidence)</td>
<td>292.645 ns [291.959, 293.347] (0.95 confidence)</td>
<td>1.00035</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d2/prepared</td>
<td>651.721 ns [650.215, 653.288] (0.95 confidence)</td>
<td>650.751 ns [649.089, 652.559] (0.95 confidence)</td>
<td>1.00149</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d3/adapter</td>
<td>1930.84 ns [1928.08, 1933.69] (0.95 confidence)</td>
<td>1944.81 ns [1941.6, 1948.06] (0.95 confidence)</td>
<td>0.992818</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d3/construct</td>
<td>650.285 ns [647.915, 652.7] (0.95 confidence)</td>
<td>642.874 ns [641.407, 644.436] (0.95 confidence)</td>
<td>1.01153</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d3/prepared</td>
<td>1153.93 ns [1151.84, 1156.32] (0.95 confidence)</td>
<td>1156.37 ns [1154.91, 1157.9] (0.95 confidence)</td>
<td>0.997895</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d4/adapter</td>
<td>3491.09 ns [3479.22, 3504.43] (0.95 confidence)</td>
<td>3366.98 ns [3361.18, 3372.88] (0.95 confidence)</td>
<td>1.03686</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d4/construct</td>
<td>1186.96 ns [1151.51, 1225.64] (0.95 confidence)</td>
<td>1159.07 ns [1125.2, 1195.8] (0.95 confidence)</td>
<td>1.02406</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d4/prepared</td>
<td>2161.47 ns [2151.56, 2171.45] (0.95 confidence)</td>
<td>2111.51 ns [2103.23, 2120.79] (0.95 confidence)</td>
<td>1.02366</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d5/adapter</td>
<td>5667.23 ns [5644.97, 5693.8] (0.95 confidence)</td>
<td>5625.74 ns [5604.63, 5654.33] (0.95 confidence)</td>
<td>1.00738</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d5/construct</td>
<td>2052.87 ns [2043.49, 2061.86] (0.95 confidence)</td>
<td>1991.3 ns [1982.16, 2000.44] (0.95 confidence)</td>
<td>1.03092</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d5/prepared</td>
<td>3291.06 ns [3284.56, 3297.58] (0.95 confidence)</td>
<td>3169.69 ns [3161.25, 3178.44] (0.95 confidence)</td>
<td>1.03829</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d6/adapter</td>
<td>7759.12 ns [7748.53, 7770.44] (0.95 confidence)</td>
<td>7794.85 ns [7780.1, 7810.4] (0.95 confidence)</td>
<td>0.995417</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d6/construct</td>
<td>2608.79 ns [2598.82, 2618.85] (0.95 confidence)</td>
<td>2666.76 ns [2651.92, 2682.24] (0.95 confidence)</td>
<td>0.978259</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d6/prepared</td>
<td>4755.69 ns [4750.46, 4761.11] (0.95 confidence)</td>
<td>4760.52 ns [4749.71, 4771.65] (0.95 confidence)</td>
<td>0.998985</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d7/adapter</td>
<td>10483.5 ns [10464.8, 10502.4] (0.95 confidence)</td>
<td>10495.6 ns [10480.4, 10510.6] (0.95 confidence)</td>
<td>0.998841</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d7/construct</td>
<td>3548.84 ns [3516.98, 3595.33] (0.95 confidence)</td>
<td>3571.94 ns [3544.88, 3613.92] (0.95 confidence)</td>
<td>0.993531</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d7/prepared</td>
<td>6720.31 ns [6689.47, 6751.93] (0.95 confidence)</td>
<td>6541 ns [6524.18, 6558.71] (0.95 confidence)</td>
<td>1.02741</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d8/adapter</td>
<td>13641.3 ns [13618.9, 13664.3] (0.95 confidence)</td>
<td>13583.3 ns [13558, 13611.3] (0.95 confidence)</td>
<td>1.00427</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d8/construct</td>
<td>4126.39 ns [4112.24, 4140.34] (0.95 confidence)</td>
<td>4166.66 ns [4155.69, 4177.77] (0.95 confidence)</td>
<td>0.990336</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d8/prepared</td>
<td>9329.67 ns [9300.68, 9360.42] (0.95 confidence)</td>
<td>9141.35 ns [9115, 9168.62] (0.95 confidence)</td>
<td>1.0206</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d2/adapter</td>
<td>1071.63 ns [1063.34, 1080.4] (0.95 confidence)</td>
<td>1075.18 ns [1066.7, 1084.09] (0.95 confidence)</td>
<td>0.996702</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d2/construct</td>
<td>328.869 ns [328.119, 329.591] (0.95 confidence)</td>
<td>328.274 ns [327.239, 329.233] (0.95 confidence)</td>
<td>1.00181</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d2/prepared</td>
<td>633.586 ns [632.546, 634.639] (0.95 confidence)</td>
<td>630.965 ns [629.866, 632.133] (0.95 confidence)</td>
<td>1.00415</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d3/adapter</td>
<td>2112.69 ns [2109.56, 2115.93] (0.95 confidence)</td>
<td>2113.13 ns [2110.67, 2115.66] (0.95 confidence)</td>
<td>0.999795</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d3/construct</td>
<td>683.186 ns [680.863, 685.544] (0.95 confidence)</td>
<td>684.219 ns [682.012, 686.378] (0.95 confidence)</td>
<td>0.998489</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d3/prepared</td>
<td>1299.79 ns [1297.96, 1301.73] (0.95 confidence)</td>
<td>1305.99 ns [1304.72, 1307.35] (0.95 confidence)</td>
<td>0.995249</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d4/adapter</td>
<td>3742.74 ns [3719.1, 3766.62] (0.95 confidence)</td>
<td>3535.96 ns [3532.58, 3539.59] (0.95 confidence)</td>
<td>1.05848</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d4/construct</td>
<td>1285.63 ns [1247.6, 1327.12] (0.95 confidence)</td>
<td>1246.53 ns [1209.36, 1286.69] (0.95 confidence)</td>
<td>1.03137</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d4/prepared</td>
<td>2281.27 ns [2275.41, 2287.35] (0.95 confidence)</td>
<td>2204.14 ns [2201.05, 2207.2] (0.95 confidence)</td>
<td>1.03499</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d5/adapter</td>
<td>5711.13 ns [5705.02, 5717.44] (0.95 confidence)</td>
<td>5645.36 ns [5635.92, 5655.51] (0.95 confidence)</td>
<td>1.01165</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d5/construct</td>
<td>1995.22 ns [1988.29, 2002.08] (0.95 confidence)</td>
<td>1956.59 ns [1948.93, 1964.52] (0.95 confidence)</td>
<td>1.01974</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d5/prepared</td>
<td>3482.24 ns [3476.76, 3488.04] (0.95 confidence)</td>
<td>3474.49 ns [3463.36, 3487.25] (0.95 confidence)</td>
<td>1.00223</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d6/adapter</td>
<td>8106.05 ns [8097.79, 8114.69] (0.95 confidence)</td>
<td>8065.61 ns [8050.46, 8085.26] (0.95 confidence)</td>
<td>1.00501</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d6/construct</td>
<td>2890.37 ns [2881.4, 2899.18] (0.95 confidence)</td>
<td>2807.55 ns [2793.44, 2821.05] (0.95 confidence)</td>
<td>1.0295</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d6/prepared</td>
<td>5072.81 ns [5057.09, 5089.65] (0.95 confidence)</td>
<td>4949.4 ns [4942.77, 4956.75] (0.95 confidence)</td>
<td>1.02493</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d7/adapter</td>
<td>11675.7 ns [11664.1, 11687.9] (0.95 confidence)</td>
<td>11910.7 ns [11883.5, 11940.1] (0.95 confidence)</td>
<td>0.98027</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d7/construct</td>
<td>3938.78 ns [3904.01, 3987.4] (0.95 confidence)</td>
<td>4028.39 ns [4000.22, 4067.44] (0.95 confidence)</td>
<td>0.977755</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d7/prepared</td>
<td>7386.74 ns [7377.55, 7396.37] (0.95 confidence)</td>
<td>7510.85 ns [7492.54, 7529.83] (0.95 confidence)</td>
<td>0.983475</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d8/adapter</td>
<td>16179.8 ns [16153.5, 16206.5] (0.95 confidence)</td>
<td>16187.8 ns [16161.8, 16214] (0.95 confidence)</td>
<td>0.999505</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d8/construct</td>
<td>5195.81 ns [5171.1, 5219.96] (0.95 confidence)</td>
<td>5229.12 ns [5205.28, 5252.35] (0.95 confidence)</td>
<td>0.99363</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d8/prepared</td>
<td>10556.9 ns [10543.6, 10570] (0.95 confidence)</td>
<td>10605.4 ns [10590.1, 10621.8] (0.95 confidence)</td>
<td>0.995425</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d2/adapter</td>
<td>454.548 ns [447.184, 462.343] (0.95 confidence)</td>
<td>455.263 ns [448.01, 463.112] (0.95 confidence)</td>
<td>0.998429</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d2/construct</td>
<td>225.677 ns [225.329, 226.032] (0.95 confidence)</td>
<td>233.237 ns [231.496, 235.116] (0.95 confidence)</td>
<td>0.967586</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d2/prepared</td>
<td>153.995 ns [153.768, 154.232] (0.95 confidence)</td>
<td>157.588 ns [157.114, 158.069] (0.95 confidence)</td>
<td>0.9772</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d3/adapter</td>
<td>750.863 ns [749.267, 752.572] (0.95 confidence)</td>
<td>740.314 ns [739.102, 741.566] (0.95 confidence)</td>
<td>1.01425</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d3/construct</td>
<td>411.312 ns [410.383, 412.252] (0.95 confidence)</td>
<td>411.864 ns [411.063, 412.691] (0.95 confidence)</td>
<td>0.99866</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d3/prepared</td>
<td>250.512 ns [249.872, 251.166] (0.95 confidence)</td>
<td>250.183 ns [249.772, 250.64] (0.95 confidence)</td>
<td>1.00131</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d4/adapter</td>
<td>1039.27 ns [1037.6, 1041] (0.95 confidence)</td>
<td>1031.78 ns [1029.8, 1034] (0.95 confidence)</td>
<td>1.00726</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d4/construct</td>
<td>636.878 ns [626.952, 647.621] (0.95 confidence)</td>
<td>634.259 ns [625.151, 644.117] (0.95 confidence)</td>
<td>1.00413</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d4/prepared</td>
<td>329.798 ns [329.291, 330.319] (0.95 confidence)</td>
<td>325.819 ns [325.248, 326.43] (0.95 confidence)</td>
<td>1.01221</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d5/adapter</td>
<td>1432.34 ns [1429.65, 1434.99] (0.95 confidence)</td>
<td>1393.91 ns [1391.44, 1396.42] (0.95 confidence)</td>
<td>1.02757</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d5/construct</td>
<td>881.104 ns [878.325, 883.916] (0.95 confidence)</td>
<td>857.379 ns [853.943, 860.798] (0.95 confidence)</td>
<td>1.02767</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d5/prepared</td>
<td>442.523 ns [439.906, 446.814] (0.95 confidence)</td>
<td>422.429 ns [422.02, 422.849] (0.95 confidence)</td>
<td>1.04757</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d6/adapter</td>
<td>1847.66 ns [1835.5, 1861.01] (0.95 confidence)</td>
<td>1760.24 ns [1755.35, 1765.53] (0.95 confidence)</td>
<td>1.04966</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d6/construct</td>
<td>1130.12 ns [1127.56, 1132.76] (0.95 confidence)</td>
<td>1117.49 ns [1114, 1121.09] (0.95 confidence)</td>
<td>1.01131</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d6/prepared</td>
<td>532.164 ns [528.086, 536.576] (0.95 confidence)</td>
<td>514.208 ns [512.808, 515.78] (0.95 confidence)</td>
<td>1.03492</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d7/adapter</td>
<td>2202.93 ns [2198.16, 2207.95] (0.95 confidence)</td>
<td>2218.78 ns [2213.76, 2223.97] (0.95 confidence)</td>
<td>0.992859</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d7/construct</td>
<td>1477.96 ns [1468.28, 1489.32] (0.95 confidence)</td>
<td>1446.68 ns [1440.03, 1454.77] (0.95 confidence)</td>
<td>1.02162</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d7/prepared</td>
<td>612.205 ns [611.154, 613.319] (0.95 confidence)</td>
<td>608.342 ns [607.02, 609.743] (0.95 confidence)</td>
<td>1.00635</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d8/adapter</td>
<td>2516.32 ns [2512.71, 2519.91] (0.95 confidence)</td>
<td>2519.83 ns [2515.53, 2524.32] (0.95 confidence)</td>
<td>0.998609</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d8/construct</td>
<td>1736.21 ns [1731.79, 1740.62] (0.95 confidence)</td>
<td>1741.1 ns [1734.45, 1747.73] (0.95 confidence)</td>
<td>0.997188</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d8/prepared</td>
<td>646.41 ns [645.717, 647.128] (0.95 confidence)</td>
<td>650.503 ns [649.553, 651.474] (0.95 confidence)</td>
<td>0.993708</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d2/adapter</td>
<td>596.583 ns [588.964, 604.572] (0.95 confidence)</td>
<td>595.626 ns [587.801, 604.058] (0.95 confidence)</td>
<td>1.00161</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d2/construct</td>
<td>307.637 ns [306.468, 308.792] (0.95 confidence)</td>
<td>312.427 ns [310.565, 314.215] (0.95 confidence)</td>
<td>0.98467</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d2/prepared</td>
<td>190.281 ns [189.938, 190.622] (0.95 confidence)</td>
<td>199.659 ns [197.934, 201.414] (0.95 confidence)</td>
<td>0.95303</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d3/adapter</td>
<td>917.222 ns [915.57, 918.971] (0.95 confidence)</td>
<td>915.647 ns [913.675, 917.859] (0.95 confidence)</td>
<td>1.00172</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d3/construct</td>
<td>536.285 ns [534.414, 538.155] (0.95 confidence)</td>
<td>537.437 ns [535.822, 538.966] (0.95 confidence)</td>
<td>0.997857</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d3/prepared</td>
<td>277.257 ns [276.886, 277.638] (0.95 confidence)</td>
<td>277.969 ns [277.508, 278.453] (0.95 confidence)</td>
<td>0.997438</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d4/adapter</td>
<td>1277.76 ns [1275.91, 1279.68] (0.95 confidence)</td>
<td>1276.61 ns [1273.11, 1280.31] (0.95 confidence)</td>
<td>1.0009</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d4/construct</td>
<td>807.23 ns [797.338, 817.771] (0.95 confidence)</td>
<td>812.513 ns [802.328, 823.485] (0.95 confidence)</td>
<td>0.993498</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d4/prepared</td>
<td>373.469 ns [372.825, 374.141] (0.95 confidence)</td>
<td>373.455 ns [372.403, 374.6] (0.95 confidence)</td>
<td>1.00004</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d5/adapter</td>
<td>1742.57 ns [1740.3, 1744.93] (0.95 confidence)</td>
<td>1735.18 ns [1730.38, 1740.05] (0.95 confidence)</td>
<td>1.00426</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d5/construct</td>
<td>1128.86 ns [1123.21, 1134.82] (0.95 confidence)</td>
<td>1110.54 ns [1104.79, 1116.6] (0.95 confidence)</td>
<td>1.01649</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d5/prepared</td>
<td>509.166 ns [506.625, 511.789] (0.95 confidence)</td>
<td>479.759 ns [478.884, 480.718] (0.95 confidence)</td>
<td>1.0613</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d6/adapter</td>
<td>2173.94 ns [2169.24, 2178.71] (0.95 confidence)</td>
<td>2121.58 ns [2113.94, 2130.33] (0.95 confidence)</td>
<td>1.02468</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d6/construct</td>
<td>1438.78 ns [1433.69, 1443.77] (0.95 confidence)</td>
<td>1432.91 ns [1424.94, 1441.14] (0.95 confidence)</td>
<td>1.00409</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d6/prepared</td>
<td>580.952 ns [578.42, 583.523] (0.95 confidence)</td>
<td>572.536 ns [571.121, 574.05] (0.95 confidence)</td>
<td>1.0147</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d7/adapter</td>
<td>2646.99 ns [2640.68, 2653.49] (0.95 confidence)</td>
<td>2636.38 ns [2630.03, 2643.12] (0.95 confidence)</td>
<td>1.00403</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d7/construct</td>
<td>1846.46 ns [1834.18, 1859.43] (0.95 confidence)</td>
<td>1815.24 ns [1806.52, 1825.5] (0.95 confidence)</td>
<td>1.0172</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d7/prepared</td>
<td>681.978 ns [678.98, 685.472] (0.95 confidence)</td>
<td>684.379 ns [681.93, 686.927] (0.95 confidence)</td>
<td>0.996493</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d8/adapter</td>
<td>3055 ns [3049.57, 3060.45] (0.95 confidence)</td>
<td>3044.04 ns [3037.84, 3050.36] (0.95 confidence)</td>
<td>1.0036</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d8/construct</td>
<td>2136.72 ns [2126.22, 2146.51] (0.95 confidence)</td>
<td>2152.79 ns [2141.79, 2162.83] (0.95 confidence)</td>
<td>0.992537</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d8/prepared</td>
<td>763.131 ns [760.985, 765.408] (0.95 confidence)</td>
<td>775.928 ns [773.805, 778.123] (0.95 confidence)</td>
<td>0.983507</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d2/adapter</td>
<td>545.505 ns [538.646, 552.876] (0.95 confidence)</td>
<td>543.759 ns [536.943, 550.931] (0.95 confidence)</td>
<td>1.00321</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d2/construct</td>
<td>241.451 ns [240.955, 241.99] (0.95 confidence)</td>
<td>242.139 ns [241.594, 242.747] (0.95 confidence)</td>
<td>0.997157</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d2/prepared</td>
<td>223.021 ns [222.648, 223.419] (0.95 confidence)</td>
<td>221.286 ns [220.941, 221.65] (0.95 confidence)</td>
<td>1.00784</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d3/adapter</td>
<td>1049.95 ns [1047.1, 1052.91] (0.95 confidence)</td>
<td>1041.96 ns [1039.87, 1044.07] (0.95 confidence)</td>
<td>1.00767</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d3/construct</td>
<td>507.765 ns [505.623, 509.915] (0.95 confidence)</td>
<td>487.785 ns [486.497, 489.124] (0.95 confidence)</td>
<td>1.04096</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d3/prepared</td>
<td>463.041 ns [461.706, 464.442] (0.95 confidence)</td>
<td>452.405 ns [451.845, 452.978] (0.95 confidence)</td>
<td>1.02351</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d4/adapter</td>
<td>1642.58 ns [1640.01, 1645.13] (0.95 confidence)</td>
<td>1617.62 ns [1615.41, 1619.86] (0.95 confidence)</td>
<td>1.01543</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d4/construct</td>
<td>921.169 ns [896.633, 947.844] (0.95 confidence)</td>
<td>885.723 ns [862.456, 911.017] (0.95 confidence)</td>
<td>1.04002</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d4/prepared</td>
<td>686.223 ns [684.582, 687.844] (0.95 confidence)</td>
<td>661.043 ns [660.29, 661.823] (0.95 confidence)</td>
<td>1.03809</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d5/adapter</td>
<td>2449.97 ns [2445.52, 2454.4] (0.95 confidence)</td>
<td>2450.52 ns [2445.84, 2455.28] (0.95 confidence)</td>
<td>0.999774</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d5/construct</td>
<td>1298.59 ns [1293.53, 1303.51] (0.95 confidence)</td>
<td>1254.54 ns [1250.29, 1258.93] (0.95 confidence)</td>
<td>1.03511</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d5/prepared</td>
<td>1000.63 ns [998.429, 1002.91] (0.95 confidence)</td>
<td>984.233 ns [982.196, 986.38] (0.95 confidence)</td>
<td>1.01666</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d6/adapter</td>
<td>3402.89 ns [3398.13, 3407.64] (0.95 confidence)</td>
<td>3346.05 ns [3333.52, 3360.98] (0.95 confidence)</td>
<td>1.01699</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d6/construct</td>
<td>1792.2 ns [1785.63, 1798.77] (0.95 confidence)</td>
<td>1770.12 ns [1763.72, 1776.55] (0.95 confidence)</td>
<td>1.01247</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d6/prepared</td>
<td>1421.62 ns [1412.06, 1431.24] (0.95 confidence)</td>
<td>1344.07 ns [1341.25, 1347.01] (0.95 confidence)</td>
<td>1.0577</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d7/adapter</td>
<td>4279.9 ns [4274.45, 4285.57] (0.95 confidence)</td>
<td>4341.33 ns [4335.44, 4347.68] (0.95 confidence)</td>
<td>0.985851</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d7/construct</td>
<td>2313.19 ns [2295.35, 2340.69] (0.95 confidence)</td>
<td>2351.63 ns [2334.97, 2377.96] (0.95 confidence)</td>
<td>0.983654</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d7/prepared</td>
<td>1743.67 ns [1741.7, 1745.71] (0.95 confidence)</td>
<td>1754.46 ns [1751.89, 1757.1] (0.95 confidence)</td>
<td>0.993846</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d8/adapter</td>
<td>5287.89 ns [5281.28, 5294.6] (0.95 confidence)</td>
<td>5326.62 ns [5314.78, 5339.13] (0.95 confidence)</td>
<td>0.992729</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d8/construct</td>
<td>2934.39 ns [2926.69, 2942.47] (0.95 confidence)</td>
<td>3109.86 ns [3086.52, 3135.01] (0.95 confidence)</td>
<td>0.943575</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d8/prepared</td>
<td>2124.66 ns [2121.86, 2127.58] (0.95 confidence)</td>
<td>2136.09 ns [2132.99, 2139.29] (0.95 confidence)</td>
<td>0.994651</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d2/adapter</td>
<td>669.891 ns [662.794, 677.372] (0.95 confidence)</td>
<td>667.702 ns [660.463, 675.472] (0.95 confidence)</td>
<td>1.00328</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d2/construct</td>
<td>348.487 ns [347.228, 349.782] (0.95 confidence)</td>
<td>349.858 ns [348.643, 351.055] (0.95 confidence)</td>
<td>0.996081</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d2/prepared</td>
<td>215.696 ns [215.342, 216.054] (0.95 confidence)</td>
<td>217.36 ns [216.881, 217.881] (0.95 confidence)</td>
<td>0.992345</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d3/adapter</td>
<td>1170.08 ns [1168.39, 1171.9] (0.95 confidence)</td>
<td>1156.66 ns [1155.04, 1158.38] (0.95 confidence)</td>
<td>1.0116</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d3/construct</td>
<td>664.642 ns [661.939, 667.256] (0.95 confidence)</td>
<td>659.896 ns [656.58, 663.115] (0.95 confidence)</td>
<td>1.00719</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d3/prepared</td>
<td>369.414 ns [368.476, 370.407] (0.95 confidence)</td>
<td>365.35 ns [364.771, 365.951] (0.95 confidence)</td>
<td>1.01112</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d4/adapter</td>
<td>1883.72 ns [1880.71, 1886.7] (0.95 confidence)</td>
<td>1813.14 ns [1810.85, 1815.54] (0.95 confidence)</td>
<td>1.03893</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d4/construct</td>
<td>1162.06 ns [1135.48, 1190.7] (0.95 confidence)</td>
<td>1127.24 ns [1101.51, 1155.05] (0.95 confidence)</td>
<td>1.03089</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d4/prepared</td>
<td>613.059 ns [612.214, 613.922] (0.95 confidence)</td>
<td>589.471 ns [588.851, 590.104] (0.95 confidence)</td>
<td>1.04001</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d5/adapter</td>
<td>2622.72 ns [2619.11, 2626.39] (0.95 confidence)</td>
<td>2603.63 ns [2598.19, 2609.35] (0.95 confidence)</td>
<td>1.00733</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d5/construct</td>
<td>1602.64 ns [1595.93, 1609.01] (0.95 confidence)</td>
<td>1630.76 ns [1618.85, 1645.15] (0.95 confidence)</td>
<td>0.982756</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d5/prepared</td>
<td>823.207 ns [821.823, 824.695] (0.95 confidence)</td>
<td>813.669 ns [812.346, 815.067] (0.95 confidence)</td>
<td>1.01172</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d6/adapter</td>
<td>3691.14 ns [3674.67, 3708.77] (0.95 confidence)</td>
<td>3574.68 ns [3566.12, 3584.65] (0.95 confidence)</td>
<td>1.03258</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d6/construct</td>
<td>2343.7 ns [2335.69, 2350.94] (0.95 confidence)</td>
<td>2322.94 ns [2313.93, 2331.68] (0.95 confidence)</td>
<td>1.00894</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d6/prepared</td>
<td>1060.72 ns [1057.64, 1063.92] (0.95 confidence)</td>
<td>1041.55 ns [1038.81, 1044.29] (0.95 confidence)</td>
<td>1.01841</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d7/adapter</td>
<td>4740.49 ns [4719.67, 4768.24] (0.95 confidence)</td>
<td>4713.3 ns [4701.2, 4726.54] (0.95 confidence)</td>
<td>1.00577</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d7/construct</td>
<td>3053.67 ns [3036.06, 3079.43] (0.95 confidence)</td>
<td>3039.74 ns [3014.78, 3070.47] (0.95 confidence)</td>
<td>1.00458</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d7/prepared</td>
<td>1394.41 ns [1387.94, 1401.44] (0.95 confidence)</td>
<td>1385.59 ns [1382.07, 1389.31] (0.95 confidence)</td>
<td>1.00636</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d8/adapter</td>
<td>5839.38 ns [5829.72, 5849.27] (0.95 confidence)</td>
<td>5917.97 ns [5844.21, 6050.18] (0.95 confidence)</td>
<td>0.986719</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d8/construct</td>
<td>3937.71 ns [3922.28, 3953.01] (0.95 confidence)</td>
<td>3941.41 ns [3929.3, 3953.17] (0.95 confidence)</td>
<td>0.999063</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d8/prepared</td>
<td>1630.51 ns [1626.71, 1634.36] (0.95 confidence)</td>
<td>1634.99 ns [1630.86, 1639.51] (0.95 confidence)</td>
<td>0.997261</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d2/adapter</td>
<td>552.073 ns [544.384, 560.144] (0.95 confidence)</td>
<td>546.469 ns [539.468, 553.91] (0.95 confidence)</td>
<td>1.01025</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d2/construct</td>
<td>244.137 ns [243.604, 244.704] (0.95 confidence)</td>
<td>243.811 ns [243.213, 244.473] (0.95 confidence)</td>
<td>1.00134</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d2/prepared</td>
<td>223.737 ns [223.196, 224.309] (0.95 confidence)</td>
<td>220.611 ns [220.229, 221] (0.95 confidence)</td>
<td>1.01417</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d3/adapter</td>
<td>1045.34 ns [1042.39, 1048.44] (0.95 confidence)</td>
<td>1044.36 ns [1042.69, 1045.99] (0.95 confidence)</td>
<td>1.00094</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d3/construct</td>
<td>492.259 ns [490.552, 494.025] (0.95 confidence)</td>
<td>487.097 ns [485.856, 488.454] (0.95 confidence)</td>
<td>1.0106</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d3/prepared</td>
<td>456.871 ns [455.948, 457.801] (0.95 confidence)</td>
<td>451.529 ns [450.875, 452.206] (0.95 confidence)</td>
<td>1.01183</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d4/adapter</td>
<td>1658.43 ns [1653.22, 1664.13] (0.95 confidence)</td>
<td>1623.27 ns [1621.15, 1625.4] (0.95 confidence)</td>
<td>1.02166</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d4/construct</td>
<td>906.371 ns [881.744, 933.265] (0.95 confidence)</td>
<td>888.589 ns [864.517, 914.627] (0.95 confidence)</td>
<td>1.02001</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d4/prepared</td>
<td>677.678 ns [676.695, 678.685] (0.95 confidence)</td>
<td>660.411 ns [659.662, 661.202] (0.95 confidence)</td>
<td>1.02615</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d5/adapter</td>
<td>2529.01 ns [2523.08, 2534.97] (0.95 confidence)</td>
<td>2451.19 ns [2449.07, 2453.3] (0.95 confidence)</td>
<td>1.03175</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d5/construct</td>
<td>1306.29 ns [1302.02, 1310.57] (0.95 confidence)</td>
<td>1255.9 ns [1250.39, 1261.58] (0.95 confidence)</td>
<td>1.04012</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d5/prepared</td>
<td>1029.88 ns [1027.75, 1032.14] (0.95 confidence)</td>
<td>985.794 ns [984.547, 987.101] (0.95 confidence)</td>
<td>1.04473</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d6/adapter</td>
<td>3427.55 ns [3423.04, 3432.1] (0.95 confidence)</td>
<td>3320.5 ns [3315.48, 3325.68] (0.95 confidence)</td>
<td>1.03224</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d6/construct</td>
<td>1800.85 ns [1794.05, 1807.61] (0.95 confidence)</td>
<td>1748.65 ns [1741.95, 1755.73] (0.95 confidence)</td>
<td>1.02986</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d6/prepared</td>
<td>1373.3 ns [1371.48, 1375.22] (0.95 confidence)</td>
<td>1348.92 ns [1343.66, 1354.83] (0.95 confidence)</td>
<td>1.01807</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d7/adapter</td>
<td>4295 ns [4286.83, 4303.36] (0.95 confidence)</td>
<td>4375.65 ns [4362.21, 4390.37] (0.95 confidence)</td>
<td>0.981569</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d7/construct</td>
<td>2295.64 ns [2278.96, 2322.19] (0.95 confidence)</td>
<td>2312.72 ns [2295.38, 2338.9] (0.95 confidence)</td>
<td>0.992615</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d7/prepared</td>
<td>1746.41 ns [1744.51, 1748.42] (0.95 confidence)</td>
<td>1765.29 ns [1760.08, 1770.99] (0.95 confidence)</td>
<td>0.989309</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d8/adapter</td>
<td>5287.1 ns [5279.72, 5294.85] (0.95 confidence)</td>
<td>5308.37 ns [5298.04, 5319.24] (0.95 confidence)</td>
<td>0.995994</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d8/construct</td>
<td>2954.57 ns [2944.82, 2964.25] (0.95 confidence)</td>
<td>2979.78 ns [2970.81, 2989.01] (0.95 confidence)</td>
<td>0.991539</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d8/prepared</td>
<td>2139.56 ns [2136.16, 2142.95] (0.95 confidence)</td>
<td>2168.07 ns [2164.14, 2172.15] (0.95 confidence)</td>
<td>0.986854</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d2/adapter</td>
<td>671.851 ns [664.611, 679.481] (0.95 confidence)</td>
<td>684.208 ns [676.096, 692.823] (0.95 confidence)</td>
<td>0.981941</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d2/construct</td>
<td>350.456 ns [349.002, 351.9] (0.95 confidence)</td>
<td>351.466 ns [349.915, 353.012] (0.95 confidence)</td>
<td>0.997127</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d2/prepared</td>
<td>217.478 ns [216.978, 217.988] (0.95 confidence)</td>
<td>216.837 ns [216.459, 217.216] (0.95 confidence)</td>
<td>1.00296</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d3/adapter</td>
<td>1155.51 ns [1153.21, 1157.95] (0.95 confidence)</td>
<td>1151.82 ns [1150.29, 1153.39] (0.95 confidence)</td>
<td>1.0032</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d3/construct</td>
<td>664.907 ns [662.364, 667.404] (0.95 confidence)</td>
<td>664.954 ns [662.551, 667.228] (0.95 confidence)</td>
<td>0.999929</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d3/prepared</td>
<td>367.777 ns [367.148, 368.421] (0.95 confidence)</td>
<td>364.751 ns [364.067, 365.443] (0.95 confidence)</td>
<td>1.0083</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d4/adapter</td>
<td>1846.49 ns [1843.5, 1849.54] (0.95 confidence)</td>
<td>1820.31 ns [1817.12, 1823.44] (0.95 confidence)</td>
<td>1.01438</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d4/construct</td>
<td>1138.96 ns [1112.18, 1168.19] (0.95 confidence)</td>
<td>1123.96 ns [1097.55, 1152.67] (0.95 confidence)</td>
<td>1.01335</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d4/prepared</td>
<td>603.255 ns [602.367, 604.171] (0.95 confidence)</td>
<td>593.941 ns [592.936, 595.008] (0.95 confidence)</td>
<td>1.01568</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d5/adapter</td>
<td>2713.66 ns [2708.73, 2718.71] (0.95 confidence)</td>
<td>2604.96 ns [2601.07, 2609.08] (0.95 confidence)</td>
<td>1.04173</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d5/construct</td>
<td>1647.68 ns [1641.38, 1653.71] (0.95 confidence)</td>
<td>1605.01 ns [1598.35, 1611.54] (0.95 confidence)</td>
<td>1.02659</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d5/prepared</td>
<td>838.654 ns [836.773, 840.595] (0.95 confidence)</td>
<td>815.34 ns [814.007, 816.781] (0.95 confidence)</td>
<td>1.02859</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d6/adapter</td>
<td>3663.16 ns [3657.47, 3669.03] (0.95 confidence)</td>
<td>3558.86 ns [3550.38, 3568.06] (0.95 confidence)</td>
<td>1.02931</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d6/construct</td>
<td>2384.99 ns [2374.86, 2394.91] (0.95 confidence)</td>
<td>2334.72 ns [2325.77, 2343.28] (0.95 confidence)</td>
<td>1.02153</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d6/prepared</td>
<td>1068.12 ns [1064.53, 1071.84] (0.95 confidence)</td>
<td>1049.75 ns [1045.99, 1053.58] (0.95 confidence)</td>
<td>1.0175</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d7/adapter</td>
<td>4688.49 ns [4681.59, 4695.46] (0.95 confidence)</td>
<td>4720.63 ns [4707.11, 4734.77] (0.95 confidence)</td>
<td>0.993193</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d7/construct</td>
<td>3017.19 ns [2994.4, 3047.33] (0.95 confidence)</td>
<td>3086.28 ns [3065.69, 3115.09] (0.95 confidence)</td>
<td>0.977611</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d7/prepared</td>
<td>1373.36 ns [1370.28, 1376.57] (0.95 confidence)</td>
<td>1379.57 ns [1376.44, 1382.8] (0.95 confidence)</td>
<td>0.995497</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d8/adapter</td>
<td>5858.17 ns [5848.64, 5867.9] (0.95 confidence)</td>
<td>5861.97 ns [5852.5, 5871.76] (0.95 confidence)</td>
<td>0.999352</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d8/construct</td>
<td>3891.2 ns [3878.62, 3903.69] (0.95 confidence)</td>
<td>3957.59 ns [3941, 3973.37] (0.95 confidence)</td>
<td>0.983226</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d8/prepared</td>
<td>1636.55 ns [1632.77, 1640.57] (0.95 confidence)</td>
<td>1652.59 ns [1647.96, 1657.38] (0.95 confidence)</td>
<td>0.990291</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d2/adapter</td>
<td>548.498 ns [541.98, 555.458] (0.95 confidence)</td>
<td>550.698 ns [543.848, 557.982] (0.95 confidence)</td>
<td>0.996006</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d2/construct</td>
<td>237.175 ns [236.715, 237.638] (0.95 confidence)</td>
<td>237.853 ns [237.319, 238.427] (0.95 confidence)</td>
<td>0.997151</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d2/prepared</td>
<td>222.605 ns [222.22, 223.016] (0.95 confidence)</td>
<td>222.108 ns [221.784, 222.438] (0.95 confidence)</td>
<td>1.00223</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d3/adapter</td>
<td>1014.33 ns [1012.62, 1016.13] (0.95 confidence)</td>
<td>1019.21 ns [1017.61, 1020.84] (0.95 confidence)</td>
<td>0.995217</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d3/construct</td>
<td>479.223 ns [478.116, 480.34] (0.95 confidence)</td>
<td>473.24 ns [471.643, 474.923] (0.95 confidence)</td>
<td>1.01264</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d3/prepared</td>
<td>441.899 ns [440.691, 443.267] (0.95 confidence)</td>
<td>438.852 ns [438.294, 439.424] (0.95 confidence)</td>
<td>1.00694</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d4/adapter</td>
<td>1677.13 ns [1673.88, 1680.39] (0.95 confidence)</td>
<td>1656.28 ns [1654.03, 1658.67] (0.95 confidence)</td>
<td>1.01259</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d4/construct</td>
<td>892.023 ns [872.168, 913.611] (0.95 confidence)</td>
<td>872.591 ns [852.92, 893.919] (0.95 confidence)</td>
<td>1.02227</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d4/prepared</td>
<td>718.815 ns [715.561, 722.171] (0.95 confidence)</td>
<td>692.355 ns [691.513, 693.232] (0.95 confidence)</td>
<td>1.03822</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d5/adapter</td>
<td>2419.32 ns [2411.5, 2427.37] (0.95 confidence)</td>
<td>2384.34 ns [2380.55, 2388.38] (0.95 confidence)</td>
<td>1.01467</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d5/construct</td>
<td>1239.91 ns [1234.12, 1245.67] (0.95 confidence)</td>
<td>1214.17 ns [1208.81, 1219.97] (0.95 confidence)</td>
<td>1.02119</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d5/prepared</td>
<td>964.865 ns [963.111, 966.747] (0.95 confidence)</td>
<td>954.508 ns [953.411, 955.669] (0.95 confidence)</td>
<td>1.01085</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d6/adapter</td>
<td>3658.46 ns [3615.63, 3708.7] (0.95 confidence)</td>
<td>3455.49 ns [3444.97, 3466.62] (0.95 confidence)</td>
<td>1.05874</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d6/construct</td>
<td>1892.61 ns [1884.21, 1901.21] (0.95 confidence)</td>
<td>1799.13 ns [1793.62, 1804.89] (0.95 confidence)</td>
<td>1.05196</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d6/prepared</td>
<td>1413.24 ns [1409.47, 1417.34] (0.95 confidence)</td>
<td>1399.14 ns [1395.82, 1402.4] (0.95 confidence)</td>
<td>1.01008</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d7/adapter</td>
<td>4370.38 ns [4363.75, 4377.06] (0.95 confidence)</td>
<td>4413.74 ns [4406.94, 4420.55] (0.95 confidence)</td>
<td>0.990176</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d7/construct</td>
<td>2361.01 ns [2345.61, 2384.27] (0.95 confidence)</td>
<td>2364.72 ns [2349.71, 2388.02] (0.95 confidence)</td>
<td>0.99843</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d7/prepared</td>
<td>1775.79 ns [1773.5, 1778.18] (0.95 confidence)</td>
<td>1794.32 ns [1791.21, 1797.49] (0.95 confidence)</td>
<td>0.989672</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d8/adapter</td>
<td>5457.73 ns [5450.28, 5465.68] (0.95 confidence)</td>
<td>5764.32 ns [5723.51, 5804.63] (0.95 confidence)</td>
<td>0.946812</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d8/construct</td>
<td>3088.83 ns [3077.34, 3100.78] (0.95 confidence)</td>
<td>3166.91 ns [3153.98, 3179.62] (0.95 confidence)</td>
<td>0.975345</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d8/prepared</td>
<td>2108.46 ns [2106.13, 2110.95] (0.95 confidence)</td>
<td>2135.11 ns [2130.46, 2140.37] (0.95 confidence)</td>
<td>0.987518</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d2/adapter</td>
<td>706.439 ns [699.405, 713.914] (0.95 confidence)</td>
<td>705.931 ns [698.793, 713.568] (0.95 confidence)</td>
<td>1.00072</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d2/construct</td>
<td>373.226 ns [371.753, 374.647] (0.95 confidence)</td>
<td>374.471 ns [372.248, 376.5] (0.95 confidence)</td>
<td>0.996674</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d2/prepared</td>
<td>230.484 ns [230.066, 230.922] (0.95 confidence)</td>
<td>232.948 ns [232.46, 233.451] (0.95 confidence)</td>
<td>0.98942</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d3/adapter</td>
<td>1187.63 ns [1185.73, 1189.61] (0.95 confidence)</td>
<td>1182.09 ns [1180.33, 1183.93] (0.95 confidence)</td>
<td>1.00469</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d3/construct</td>
<td>695.239 ns [692.695, 697.546] (0.95 confidence)</td>
<td>687.486 ns [684.263, 690.275] (0.95 confidence)</td>
<td>1.01128</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d3/prepared</td>
<td>386.14 ns [385.328, 387.011] (0.95 confidence)</td>
<td>381.591 ns [380.972, 382.204] (0.95 confidence)</td>
<td>1.01192</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d4/adapter</td>
<td>1817.68 ns [1813.78, 1821.69] (0.95 confidence)</td>
<td>1759.45 ns [1757.49, 1761.49] (0.95 confidence)</td>
<td>1.03309</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d4/construct</td>
<td>1090.25 ns [1069.42, 1112.64] (0.95 confidence)</td>
<td>1051.3 ns [1029.79, 1074.36] (0.95 confidence)</td>
<td>1.03704</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d4/prepared</td>
<td>637.002 ns [635.751, 638.342] (0.95 confidence)</td>
<td>624.001 ns [623.21, 624.812] (0.95 confidence)</td>
<td>1.02083</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d5/adapter</td>
<td>2521.52 ns [2514.42, 2529.44] (0.95 confidence)</td>
<td>2502.14 ns [2498.06, 2506.65] (0.95 confidence)</td>
<td>1.00775</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d5/construct</td>
<td>1455.82 ns [1447.18, 1463.84] (0.95 confidence)</td>
<td>1445.94 ns [1438.39, 1452.92] (0.95 confidence)</td>
<td>1.00683</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d5/prepared</td>
<td>883.279 ns [881.872, 884.779] (0.95 confidence)</td>
<td>876.38 ns [875.37, 877.481] (0.95 confidence)</td>
<td>1.00787</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d6/adapter</td>
<td>3531.14 ns [3521.11, 3541.77] (0.95 confidence)</td>
<td>3397.57 ns [3390.61, 3404.7] (0.95 confidence)</td>
<td>1.03931</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d6/construct</td>
<td>2097.3 ns [2089.87, 2104.67] (0.95 confidence)</td>
<td>2050.17 ns [2040.08, 2059.7] (0.95 confidence)</td>
<td>1.02299</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d6/prepared</td>
<td>1172.73 ns [1170.13, 1175.47] (0.95 confidence)</td>
<td>1139.64 ns [1137.4, 1142] (0.95 confidence)</td>
<td>1.02903</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d7/adapter</td>
<td>4445.71 ns [4437.72, 4453.92] (0.95 confidence)</td>
<td>4522.44 ns [4506.11, 4539.72] (0.95 confidence)</td>
<td>0.983033</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d7/construct</td>
<td>2775.73 ns [2757.77, 2799.97] (0.95 confidence)</td>
<td>2833.87 ns [2817.99, 2857.81] (0.95 confidence)</td>
<td>0.979483</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d7/prepared</td>
<td>1470.06 ns [1468.05, 1472.19] (0.95 confidence)</td>
<td>1478.61 ns [1475.59, 1481.74] (0.95 confidence)</td>
<td>0.994219</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d8/adapter</td>
<td>5773.84 ns [5758.41, 5789.5] (0.95 confidence)</td>
<td>5751.11 ns [5736.83, 5766.55] (0.95 confidence)</td>
<td>1.00395</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d8/construct</td>
<td>3567.85 ns [3557.47, 3578.03] (0.95 confidence)</td>
<td>3615.74 ns [3601.36, 3630.31] (0.95 confidence)</td>
<td>0.986754</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d8/prepared</td>
<td>1870.18 ns [1867.92, 1872.49] (0.95 confidence)</td>
<td>1887.13 ns [1883.78, 1890.64] (0.95 confidence)</td>
<td>0.991016</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d2/adapter</td>
<td>1058 ns [1043.92, 1072.95] (0.95 confidence)</td>
<td>1063.83 ns [1049.44, 1079.25] (0.95 confidence)</td>
<td>0.994518</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d2/construct</td>
<td>295.764 ns [294.696, 296.979] (0.95 confidence)</td>
<td>295.107 ns [294.477, 295.759] (0.95 confidence)</td>
<td>1.00222</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d2/prepared</td>
<td>636.166 ns [635.349, 637.004] (0.95 confidence)</td>
<td>631.189 ns [630.129, 632.291] (0.95 confidence)</td>
<td>1.00789</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d3/adapter</td>
<td>1806.37 ns [1802.59, 1810.3] (0.95 confidence)</td>
<td>1782.51 ns [1778.98, 1786.13] (0.95 confidence)</td>
<td>1.01338</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d3/construct</td>
<td>585.52 ns [584.022, 587.055] (0.95 confidence)</td>
<td>582.377 ns [581.369, 583.408] (0.95 confidence)</td>
<td>1.0054</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d3/prepared</td>
<td>1091.34 ns [1089.43, 1093.32] (0.95 confidence)</td>
<td>1087.43 ns [1085.22, 1089.73] (0.95 confidence)</td>
<td>1.0036</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d4/adapter</td>
<td>2878.24 ns [2868.17, 2888.67] (0.95 confidence)</td>
<td>2744.81 ns [2740.49, 2749.41] (0.95 confidence)</td>
<td>1.04861</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d4/construct</td>
<td>886.125 ns [863.276, 910.945] (0.95 confidence)</td>
<td>856.386 ns [836.646, 877.81] (0.95 confidence)</td>
<td>1.03473</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d4/prepared</td>
<td>1998.91 ns [1983.55, 2016.12] (0.95 confidence)</td>
<td>1783.56 ns [1781.08, 1786.1] (0.95 confidence)</td>
<td>1.12074</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d5/adapter</td>
<td>4301.74 ns [4293.81, 4310.23] (0.95 confidence)</td>
<td>4212.13 ns [4207.92, 4216.45] (0.95 confidence)</td>
<td>1.02127</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d5/construct</td>
<td>1377.03 ns [1371.39, 1382.74] (0.95 confidence)</td>
<td>1342.31 ns [1338.01, 1346.69] (0.95 confidence)</td>
<td>1.02587</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d5/prepared</td>
<td>2707.91 ns [2703.08, 2712.9] (0.95 confidence)</td>
<td>2631.43 ns [2626.95, 2636.04] (0.95 confidence)</td>
<td>1.02906</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d6/adapter</td>
<td>5800.95 ns [5788.99, 5812.7] (0.95 confidence)</td>
<td>5697.51 ns [5683.17, 5712.67] (0.95 confidence)</td>
<td>1.01815</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d6/construct</td>
<td>1761.68 ns [1755.54, 1767.74] (0.95 confidence)</td>
<td>1700.45 ns [1695.1, 1705.86] (0.95 confidence)</td>
<td>1.03601</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d6/prepared</td>
<td>3769.45 ns [3759.8, 3779.56] (0.95 confidence)</td>
<td>3689.57 ns [3682.19, 3697.84] (0.95 confidence)</td>
<td>1.02165</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d7/adapter</td>
<td>7830.81 ns [7817.06, 7845.35] (0.95 confidence)</td>
<td>7870.46 ns [7858.62, 7883.86] (0.95 confidence)</td>
<td>0.994962</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d7/construct</td>
<td>2204.14 ns [2190.88, 2222.03] (0.95 confidence)</td>
<td>2249.04 ns [2235.26, 2268.67] (0.95 confidence)</td>
<td>0.980033</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d7/prepared</td>
<td>5267.29 ns [5261.39, 5273.58] (0.95 confidence)</td>
<td>5304.65 ns [5296.02, 5313.45] (0.95 confidence)</td>
<td>0.992956</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d8/adapter</td>
<td>9516.21 ns [9505.65, 9527.18] (0.95 confidence)</td>
<td>9572.15 ns [9558.89, 9586.41] (0.95 confidence)</td>
<td>0.994156</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d8/construct</td>
<td>2414.99 ns [2407.56, 2422.54] (0.95 confidence)</td>
<td>2435.83 ns [2428.92, 2442.8] (0.95 confidence)</td>
<td>0.991445</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d8/prepared</td>
<td>6997.65 ns [6989.07, 7006.52] (0.95 confidence)</td>
<td>7052.73 ns [7034.33, 7071.66] (0.95 confidence)</td>
<td>0.99219</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d2/adapter</td>
<td>1190.64 ns [1182.99, 1198.69] (0.95 confidence)</td>
<td>1195.55 ns [1187.92, 1203.58] (0.95 confidence)</td>
<td>0.995887</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d2/construct</td>
<td>477.943 ns [476.814, 479.067] (0.95 confidence)</td>
<td>479.22 ns [476.411, 481.743] (0.95 confidence)</td>
<td>0.997334</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d2/prepared</td>
<td>585.862 ns [584.94, 586.84] (0.95 confidence)</td>
<td>587.333 ns [586.405, 588.315] (0.95 confidence)</td>
<td>0.997496</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d3/adapter</td>
<td>1918.6 ns [1915.81, 1921.48] (0.95 confidence)</td>
<td>1920.13 ns [1917.62, 1922.78] (0.95 confidence)</td>
<td>0.999199</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d3/construct</td>
<td>752.746 ns [750.622, 754.834] (0.95 confidence)</td>
<td>750.303 ns [748.003, 752.504] (0.95 confidence)</td>
<td>1.00326</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d3/prepared</td>
<td>1033.21 ns [1031.44, 1035.15] (0.95 confidence)</td>
<td>1034.64 ns [1033.38, 1035.98] (0.95 confidence)</td>
<td>0.998614</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d4/adapter</td>
<td>3269.21 ns [3250.34, 3290.04] (0.95 confidence)</td>
<td>3159.85 ns [3155.57, 3164.51] (0.95 confidence)</td>
<td>1.03461</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d4/construct</td>
<td>1372.1 ns [1346.63, 1399.3] (0.95 confidence)</td>
<td>1343.85 ns [1317.53, 1372.37] (0.95 confidence)</td>
<td>1.02102</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d4/prepared</td>
<td>1735.47 ns [1722.58, 1749.6] (0.95 confidence)</td>
<td>1683.8 ns [1681.34, 1686.35] (0.95 confidence)</td>
<td>1.03069</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d5/adapter</td>
<td>4455.82 ns [4448.72, 4463.15] (0.95 confidence)</td>
<td>4423.47 ns [4416.84, 4430.79] (0.95 confidence)</td>
<td>1.00731</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d5/construct</td>
<td>1739.88 ns [1733.26, 1746.55] (0.95 confidence)</td>
<td>1678.08 ns [1672.81, 1683.3] (0.95 confidence)</td>
<td>1.03683</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d5/prepared</td>
<td>2583.02 ns [2577.95, 2588.11] (0.95 confidence)</td>
<td>2494.5 ns [2491.12, 2498.05] (0.95 confidence)</td>
<td>1.03549</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d6/adapter</td>
<td>6173.16 ns [6163.85, 6182.95] (0.95 confidence)</td>
<td>6103.79 ns [6087.99, 6120.61] (0.95 confidence)</td>
<td>1.01137</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d6/construct</td>
<td>2365.26 ns [2357.27, 2373.17] (0.95 confidence)</td>
<td>2276.93 ns [2270.31, 2283.62] (0.95 confidence)</td>
<td>1.0388</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d6/prepared</td>
<td>3570.22 ns [3563.67, 3577.08] (0.95 confidence)</td>
<td>3493.12 ns [3488.01, 3498.34] (0.95 confidence)</td>
<td>1.02207</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d7/adapter</td>
<td>8560.35 ns [8548.42, 8573.3] (0.95 confidence)</td>
<td>8668.9 ns [8650.21, 8689.19] (0.95 confidence)</td>
<td>0.987478</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d7/construct</td>
<td>2899.48 ns [2884.65, 2921.22] (0.95 confidence)</td>
<td>2918.55 ns [2902.34, 2940.78] (0.95 confidence)</td>
<td>0.993467</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d7/prepared</td>
<td>5425.08 ns [5420.21, 5430.55] (0.95 confidence)</td>
<td>5455.26 ns [5448.19, 5462.56] (0.95 confidence)</td>
<td>0.994467</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d8/adapter</td>
<td>10488.7 ns [10478.9, 10500] (0.95 confidence)</td>
<td>10575.7 ns [10561.5, 10590.1] (0.95 confidence)</td>
<td>0.991775</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d8/construct</td>
<td>3425.9 ns [3416.87, 3435.07] (0.95 confidence)</td>
<td>3455.24 ns [3445.89, 3464.4] (0.95 confidence)</td>
<td>0.991507</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d8/prepared</td>
<td>6870.52 ns [6860.05, 6881.55] (0.95 confidence)</td>
<td>6939.58 ns [6920.35, 6959.19] (0.95 confidence)</td>
<td>0.990049</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d2/adapter</td>
<td>547.905 ns [541.229, 555.079] (0.95 confidence)</td>
<td>548.702 ns [541.668, 556.209] (0.95 confidence)</td>
<td>0.998546</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d2/construct</td>
<td>238.249 ns [237.715, 238.802] (0.95 confidence)</td>
<td>239.884 ns [239.524, 240.246] (0.95 confidence)</td>
<td>0.993185</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d2/prepared</td>
<td>223.292 ns [222.913, 223.685] (0.95 confidence)</td>
<td>221.708 ns [221.352, 222.093] (0.95 confidence)</td>
<td>1.00714</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d3/adapter</td>
<td>1016.16 ns [1014.09, 1018.53] (0.95 confidence)</td>
<td>1019.27 ns [1017.04, 1021.92] (0.95 confidence)</td>
<td>0.996947</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d3/construct</td>
<td>480.123 ns [478.779, 481.459] (0.95 confidence)</td>
<td>477.511 ns [476.418, 478.654] (0.95 confidence)</td>
<td>1.00547</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d3/prepared</td>
<td>440.31 ns [439.754, 440.888] (0.95 confidence)</td>
<td>441.883 ns [441.398, 442.375] (0.95 confidence)</td>
<td>0.99644</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d4/adapter</td>
<td>1707.82 ns [1704.28, 1711.44] (0.95 confidence)</td>
<td>1647.51 ns [1644.88, 1650.37] (0.95 confidence)</td>
<td>1.03661</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d4/construct</td>
<td>896.248 ns [875.885, 918.633] (0.95 confidence)</td>
<td>874.547 ns [855.468, 895.088] (0.95 confidence)</td>
<td>1.02481</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d4/prepared</td>
<td>721.649 ns [720.24, 723.127] (0.95 confidence)</td>
<td>690.309 ns [689.434, 691.26] (0.95 confidence)</td>
<td>1.0454</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d5/adapter</td>
<td>2409.46 ns [2403.54, 2415.27] (0.95 confidence)</td>
<td>2370.66 ns [2364.93, 2376.52] (0.95 confidence)</td>
<td>1.01636</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d5/construct</td>
<td>1291.93 ns [1287.16, 1296.54] (0.95 confidence)</td>
<td>1233.83 ns [1227.38, 1240.49] (0.95 confidence)</td>
<td>1.04709</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d5/prepared</td>
<td>990.798 ns [989.024, 992.592] (0.95 confidence)</td>
<td>961.056 ns [958.208, 964.087] (0.95 confidence)</td>
<td>1.03095</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d6/adapter</td>
<td>3465.65 ns [3461.11, 3470.3] (0.95 confidence)</td>
<td>3432.79 ns [3423.72, 3442.43] (0.95 confidence)</td>
<td>1.00957</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d6/construct</td>
<td>1879.73 ns [1872.82, 1886.4] (0.95 confidence)</td>
<td>1829.56 ns [1823.3, 1835.74] (0.95 confidence)</td>
<td>1.02742</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d6/prepared</td>
<td>1441.55 ns [1435.19, 1448.33] (0.95 confidence)</td>
<td>1371.51 ns [1369.73, 1373.4] (0.95 confidence)</td>
<td>1.05107</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d7/adapter</td>
<td>4381.95 ns [4373.09, 4390.83] (0.95 confidence)</td>
<td>4458.99 ns [4446.32, 4472.09] (0.95 confidence)</td>
<td>0.982721</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d7/construct</td>
<td>2396.55 ns [2378.11, 2422.67] (0.95 confidence)</td>
<td>2429.53 ns [2413.39, 2454] (0.95 confidence)</td>
<td>0.986427</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d7/prepared</td>
<td>1780.98 ns [1776.41, 1786.35] (0.95 confidence)</td>
<td>1803.54 ns [1797.83, 1809.51] (0.95 confidence)</td>
<td>0.987493</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d8/adapter</td>
<td>5452.03 ns [5444.31, 5460.28] (0.95 confidence)</td>
<td>5476.84 ns [5466.81, 5487.21] (0.95 confidence)</td>
<td>0.995471</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d8/construct</td>
<td>3079.71 ns [3067.81, 3092.15] (0.95 confidence)</td>
<td>3102.78 ns [3092.58, 3113.64] (0.95 confidence)</td>
<td>0.992565</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d8/prepared</td>
<td>2109.29 ns [2105.73, 2113.06] (0.95 confidence)</td>
<td>2130.22 ns [2126.65, 2133.82] (0.95 confidence)</td>
<td>0.990173</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d2/adapter</td>
<td>704.911 ns [698.007, 712.271] (0.95 confidence)</td>
<td>709.334 ns [702.007, 717.208] (0.95 confidence)</td>
<td>0.993763</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d2/construct</td>
<td>377.61 ns [376.201, 378.983] (0.95 confidence)</td>
<td>382.691 ns [381.272, 384.014] (0.95 confidence)</td>
<td>0.986723</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d2/prepared</td>
<td>230.067 ns [229.608, 230.581] (0.95 confidence)</td>
<td>231.181 ns [230.726, 231.665] (0.95 confidence)</td>
<td>0.99518</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d3/adapter</td>
<td>1187.97 ns [1185.63, 1190.38] (0.95 confidence)</td>
<td>1182.84 ns [1180.96, 1184.79] (0.95 confidence)</td>
<td>1.00434</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d3/construct</td>
<td>701.973 ns [698.677, 704.904] (0.95 confidence)</td>
<td>695.411 ns [692.968, 697.797] (0.95 confidence)</td>
<td>1.00944</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d3/prepared</td>
<td>385.522 ns [384.836, 386.283] (0.95 confidence)</td>
<td>384.303 ns [383.611, 385.011] (0.95 confidence)</td>
<td>1.00317</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d4/adapter</td>
<td>1841 ns [1832.21, 1850.58] (0.95 confidence)</td>
<td>1770.14 ns [1765.76, 1774.82] (0.95 confidence)</td>
<td>1.04003</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d4/construct</td>
<td>1092.63 ns [1071.78, 1115.09] (0.95 confidence)</td>
<td>1059.08 ns [1039.41, 1080.35] (0.95 confidence)</td>
<td>1.03168</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d4/prepared</td>
<td>637.775 ns [636.148, 639.454] (0.95 confidence)</td>
<td>621.634 ns [620.651, 622.624] (0.95 confidence)</td>
<td>1.02597</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d5/adapter</td>
<td>2662.54 ns [2646.45, 2679.72] (0.95 confidence)</td>
<td>2526.86 ns [2520.44, 2533.74] (0.95 confidence)</td>
<td>1.05369</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d5/construct</td>
<td>1505.56 ns [1499.52, 1511.59] (0.95 confidence)</td>
<td>1466.39 ns [1460.27, 1472.53] (0.95 confidence)</td>
<td>1.02672</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d5/prepared</td>
<td>926.298 ns [923.618, 929.172] (0.95 confidence)</td>
<td>878.611 ns [876.39, 880.919] (0.95 confidence)</td>
<td>1.05428</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d6/adapter</td>
<td>3546.04 ns [3536.76, 3555.53] (0.95 confidence)</td>
<td>3418.99 ns [3406.94, 3431.63] (0.95 confidence)</td>
<td>1.03716</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d6/construct</td>
<td>2181.51 ns [2167.15, 2197.97] (0.95 confidence)</td>
<td>2023.68 ns [2016.92, 2030.63] (0.95 confidence)</td>
<td>1.07799</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d6/prepared</td>
<td>1185.46 ns [1182.5, 1188.68] (0.95 confidence)</td>
<td>1134.87 ns [1133.05, 1136.8] (0.95 confidence)</td>
<td>1.04458</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d7/adapter</td>
<td>4483.51 ns [4475.35, 4491.8] (0.95 confidence)</td>
<td>4493.01 ns [4483.66, 4502.71] (0.95 confidence)</td>
<td>0.997885</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d7/construct</td>
<td>2778.01 ns [2761.39, 2800.31] (0.95 confidence)</td>
<td>2854.55 ns [2838.48, 2873.97] (0.95 confidence)</td>
<td>0.973186</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d7/prepared</td>
<td>1476.26 ns [1473.4, 1479.19] (0.95 confidence)</td>
<td>1470.77 ns [1468.97, 1472.66] (0.95 confidence)</td>
<td>1.00373</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d8/adapter</td>
<td>5739.66 ns [5729.33, 5750.28] (0.95 confidence)</td>
<td>5753.01 ns [5740.11, 5766.26] (0.95 confidence)</td>
<td>0.99768</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d8/construct</td>
<td>3569.37 ns [3560.53, 3578.28] (0.95 confidence)</td>
<td>3614.05 ns [3602.5, 3626.16] (0.95 confidence)</td>
<td>0.987636</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d8/prepared</td>
<td>1869.19 ns [1865.36, 1873.05] (0.95 confidence)</td>
<td>1878.44 ns [1874.47, 1882.53] (0.95 confidence)</td>
<td>0.995078</td>
</tr>
</tbody>
</table>

### Median (ns)

<table>
<thead>
<tr>
<th>Workload</th>
<th>Baseline</th>
<th>Current</th>
<th>Baseline / current</th>
</tr>
</thead>
<tbody>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d2/rounded&#95;array</td>
<td>170.267 ns [169.769, 170.704] (0.95 confidence)</td>
<td>167.564 ns [167.197, 168.57] (0.95 confidence)</td>
<td>1.01614</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d2/rounded&#95;result</td>
<td>170.54 ns [170.027, 170.826] (0.95 confidence)</td>
<td>168.589 ns [168.064, 169.492] (0.95 confidence)</td>
<td>1.01157</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d2/strict&#95;array</td>
<td>174.763 ns [173.92, 175.849] (0.95 confidence)</td>
<td>170.806 ns [170.308, 171.216] (0.95 confidence)</td>
<td>1.02317</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d2/strict&#95;result</td>
<td>178.279 ns [177.626, 179.11] (0.95 confidence)</td>
<td>173.428 ns [172.221, 174.444] (0.95 confidence)</td>
<td>1.02797</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d3/rounded&#95;array</td>
<td>177.643 ns [176.163, 179.26] (0.95 confidence)</td>
<td>177.887 ns [176.08, 178.707] (0.95 confidence)</td>
<td>0.998629</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d3/rounded&#95;result</td>
<td>176.48 ns [175.979, 177.136] (0.95 confidence)</td>
<td>173.704 ns [173.316, 174.257] (0.95 confidence)</td>
<td>1.01598</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d3/strict&#95;array</td>
<td>179.981 ns [179.133, 180.577] (0.95 confidence)</td>
<td>175.173 ns [174.57, 176.187] (0.95 confidence)</td>
<td>1.02745</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d3/strict&#95;result</td>
<td>174.966 ns [174.46, 175.295] (0.95 confidence)</td>
<td>175.907 ns [175.486, 176.875] (0.95 confidence)</td>
<td>0.994649</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d4/rounded&#95;array</td>
<td>185.51 ns [185.131, 186.273] (0.95 confidence)</td>
<td>189.909 ns [184.399, 194.093] (0.95 confidence)</td>
<td>0.976836</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d4/rounded&#95;result</td>
<td>185.484 ns [184.617, 186.412] (0.95 confidence)</td>
<td>185.625 ns [184.313, 188.263] (0.95 confidence)</td>
<td>0.99924</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d4/strict&#95;array</td>
<td>187.663 ns [187.032, 188.238] (0.95 confidence)</td>
<td>196.112 ns [190.752, 198.317] (0.95 confidence)</td>
<td>0.956915</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d4/strict&#95;result</td>
<td>182.587 ns [181.266, 183.136] (0.95 confidence)</td>
<td>190.786 ns [188.298, 192.075] (0.95 confidence)</td>
<td>0.957023</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d5/rounded&#95;array</td>
<td>184.417 ns [183.74, 185.031] (0.95 confidence)</td>
<td>184.212 ns [183.784, 184.687] (0.95 confidence)</td>
<td>1.00111</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d5/rounded&#95;result</td>
<td>186.443 ns [185.066, 187.027] (0.95 confidence)</td>
<td>185.309 ns [184.795, 186.022] (0.95 confidence)</td>
<td>1.00612</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d5/strict&#95;array</td>
<td>183.468 ns [182.662, 184.106] (0.95 confidence)</td>
<td>184.653 ns [183.982, 185.118] (0.95 confidence)</td>
<td>0.99358</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d5/strict&#95;result</td>
<td>180.862 ns [180.36, 181.479] (0.95 confidence)</td>
<td>181.137 ns [180.779, 181.791] (0.95 confidence)</td>
<td>0.998479</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d2/rounded&#95;array</td>
<td>11.545 ns [11.4873, 11.6048] (0.95 confidence)</td>
<td>11.5719 ns [11.5681, 11.5769] (0.95 confidence)</td>
<td>0.997678</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d2/rounded&#95;result</td>
<td>11.4779 ns [11.4712, 11.484] (0.95 confidence)</td>
<td>11.5674 ns [11.5653, 11.5732] (0.95 confidence)</td>
<td>0.992266</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d2/strict&#95;array</td>
<td>9.0759 ns [8.96383, 9.16632] (0.95 confidence)</td>
<td>8.47547 ns [8.45159, 8.5183] (0.95 confidence)</td>
<td>1.07084</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d2/strict&#95;result</td>
<td>8.28168 ns [8.27113, 8.28831] (0.95 confidence)</td>
<td>7.82516 ns [7.80684, 7.88143] (0.95 confidence)</td>
<td>1.05834</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d3/rounded&#95;array</td>
<td>17.6553 ns [17.6504, 17.6655] (0.95 confidence)</td>
<td>18.042 ns [18.019, 18.0558] (0.95 confidence)</td>
<td>0.978566</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d3/rounded&#95;result</td>
<td>17.6679 ns [17.6581, 17.681] (0.95 confidence)</td>
<td>18.2281 ns [18.1121, 18.2878] (0.95 confidence)</td>
<td>0.969266</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d3/strict&#95;array</td>
<td>13.4109 ns [13.3516, 13.427] (0.95 confidence)</td>
<td>13.406 ns [13.24, 13.6235] (0.95 confidence)</td>
<td>1.00036</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d3/strict&#95;result</td>
<td>11.8118 ns [11.8046, 11.8359] (0.95 confidence)</td>
<td>11.8563 ns [11.7692, 11.9381] (0.95 confidence)</td>
<td>0.996247</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d4/rounded&#95;array</td>
<td>23.2124 ns [23.1678, 23.2372] (0.95 confidence)</td>
<td>24.1443 ns [24.0325, 24.1782] (0.95 confidence)</td>
<td>0.961404</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d4/rounded&#95;result</td>
<td>23.2769 ns [23.2481, 23.3023] (0.95 confidence)</td>
<td>23.9608 ns [23.9131, 24.1255] (0.95 confidence)</td>
<td>0.97146</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d4/strict&#95;array</td>
<td>18.0248 ns [17.8089, 18.1516] (0.95 confidence)</td>
<td>18.2378 ns [18.0731, 18.6103] (0.95 confidence)</td>
<td>0.988325</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d4/strict&#95;result</td>
<td>15.6387 ns [15.5278, 15.7692] (0.95 confidence)</td>
<td>16.2175 ns [16.0503, 16.2943] (0.95 confidence)</td>
<td>0.964307</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d5/rounded&#95;array</td>
<td>29.9287 ns [29.8704, 30.1126] (0.95 confidence)</td>
<td>30.7722 ns [30.1771, 31.4234] (0.95 confidence)</td>
<td>0.972589</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d5/rounded&#95;result</td>
<td>29.4416 ns [29.4013, 29.4882] (0.95 confidence)</td>
<td>30.3059 ns [29.7782, 30.8672] (0.95 confidence)</td>
<td>0.971481</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d5/strict&#95;array</td>
<td>23.3039 ns [22.878, 23.4237] (0.95 confidence)</td>
<td>23.6394 ns [23.3749, 24.1314] (0.95 confidence)</td>
<td>0.985808</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d5/strict&#95;result</td>
<td>19.9016 ns [19.8172, 20.1184] (0.95 confidence)</td>
<td>21.0316 ns [20.5584, 21.2689] (0.95 confidence)</td>
<td>0.946273</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d2/rounded&#95;array</td>
<td>176.154 ns [175.659, 176.918] (0.95 confidence)</td>
<td>173.427 ns [172.52, 174.474] (0.95 confidence)</td>
<td>1.01572</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d2/rounded&#95;result</td>
<td>177.021 ns [176.259, 177.908] (0.95 confidence)</td>
<td>172.256 ns [171.803, 172.572] (0.95 confidence)</td>
<td>1.02766</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d2/strict&#95;array</td>
<td>13.852 ns [13.6203, 13.9647] (0.95 confidence)</td>
<td>14.2051 ns [14.078, 14.3352] (0.95 confidence)</td>
<td>0.975143</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d2/strict&#95;result</td>
<td>13.5657 ns [13.4667, 13.6702] (0.95 confidence)</td>
<td>12.8597 ns [12.7515, 13.0306] (0.95 confidence)</td>
<td>1.0549</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d3/rounded&#95;array</td>
<td>180.084 ns [179.657, 180.494] (0.95 confidence)</td>
<td>180.723 ns [180.539, 181.042] (0.95 confidence)</td>
<td>0.996467</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d3/rounded&#95;result</td>
<td>178.775 ns [178.703, 178.949] (0.95 confidence)</td>
<td>180.566 ns [180.29, 181.019] (0.95 confidence)</td>
<td>0.990084</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d3/strict&#95;array</td>
<td>18.2517 ns [18.1245, 18.4488] (0.95 confidence)</td>
<td>18.4747 ns [18.3846, 18.5458] (0.95 confidence)</td>
<td>0.987929</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d3/strict&#95;result</td>
<td>17.8128 ns [17.6939, 17.8911] (0.95 confidence)</td>
<td>17.0662 ns [16.9666, 17.2187] (0.95 confidence)</td>
<td>1.04375</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d4/rounded&#95;array</td>
<td>185.954 ns [185.694, 186.108] (0.95 confidence)</td>
<td>197.543 ns [192.097, 200.645] (0.95 confidence)</td>
<td>0.941333</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d4/rounded&#95;result</td>
<td>185.802 ns [185.509, 185.993] (0.95 confidence)</td>
<td>191.204 ns [188.709, 198.579] (0.95 confidence)</td>
<td>0.97175</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d4/strict&#95;array</td>
<td>23.1134 ns [23.0268, 23.4194] (0.95 confidence)</td>
<td>25.1784 ns [24.5094, 25.7475] (0.95 confidence)</td>
<td>0.917986</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d4/strict&#95;result</td>
<td>21.8755 ns [21.6878, 21.9909] (0.95 confidence)</td>
<td>22.8644 ns [22.3874, 23.2757] (0.95 confidence)</td>
<td>0.95675</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d5/rounded&#95;array</td>
<td>191.437 ns [191.252, 191.557] (0.95 confidence)</td>
<td>192.569 ns [192.09, 192.967] (0.95 confidence)</td>
<td>0.994122</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d5/rounded&#95;result</td>
<td>192.132 ns [191.624, 193.309] (0.95 confidence)</td>
<td>196.634 ns [194.54, 199.155] (0.95 confidence)</td>
<td>0.977106</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d5/strict&#95;array</td>
<td>27.3025 ns [27.2026, 27.6397] (0.95 confidence)</td>
<td>27.2687 ns [27.2102, 27.3955] (0.95 confidence)</td>
<td>1.00124</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d5/strict&#95;result</td>
<td>26.4007 ns [26.2474, 26.5858] (0.95 confidence)</td>
<td>26.6276 ns [26.0695, 27.1] (0.95 confidence)</td>
<td>0.991481</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d2/rounded&#95;array</td>
<td>177.229 ns [176.995, 177.674] (0.95 confidence)</td>
<td>174.763 ns [174.134, 175.454] (0.95 confidence)</td>
<td>1.01411</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d2/rounded&#95;result</td>
<td>175.574 ns [174.854, 175.931] (0.95 confidence)</td>
<td>171.417 ns [171.268, 171.683] (0.95 confidence)</td>
<td>1.02425</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d2/strict&#95;array</td>
<td>13.6433 ns [13.595, 13.7211] (0.95 confidence)</td>
<td>13.3682 ns [13.3154, 13.4787] (0.95 confidence)</td>
<td>1.02057</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d2/strict&#95;result</td>
<td>13.2136 ns [13.1219, 13.3278] (0.95 confidence)</td>
<td>12.668 ns [12.599, 12.7136] (0.95 confidence)</td>
<td>1.04307</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d3/rounded&#95;array</td>
<td>177.878 ns [177.727, 178.108] (0.95 confidence)</td>
<td>179.448 ns [178.957, 180.165] (0.95 confidence)</td>
<td>0.99125</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d3/rounded&#95;result</td>
<td>178.618 ns [178.409, 178.887] (0.95 confidence)</td>
<td>181.141 ns [180.51, 182.283] (0.95 confidence)</td>
<td>0.986071</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d3/strict&#95;array</td>
<td>17.5328 ns [17.4196, 17.6198] (0.95 confidence)</td>
<td>17.6718 ns [17.5876, 17.7341] (0.95 confidence)</td>
<td>0.992131</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d3/strict&#95;result</td>
<td>16.6947 ns [16.566, 16.7958] (0.95 confidence)</td>
<td>16.9554 ns [16.7853, 17.1205] (0.95 confidence)</td>
<td>0.984622</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d4/rounded&#95;array</td>
<td>190.106 ns [189.544, 190.625] (0.95 confidence)</td>
<td>193.714 ns [192.177, 197.298] (0.95 confidence)</td>
<td>0.981376</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d4/rounded&#95;result</td>
<td>186.84 ns [186.476, 187.303] (0.95 confidence)</td>
<td>199.65 ns [194.779, 203.42] (0.95 confidence)</td>
<td>0.935836</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d4/strict&#95;array</td>
<td>22.5403 ns [22.4072, 22.715] (0.95 confidence)</td>
<td>23.344 ns [23.0782, 24.0956] (0.95 confidence)</td>
<td>0.965571</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d4/strict&#95;result</td>
<td>21.63 ns [21.5317, 21.7022] (0.95 confidence)</td>
<td>21.3579 ns [20.879, 21.656] (0.95 confidence)</td>
<td>1.01274</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d5/rounded&#95;array</td>
<td>192.637 ns [192.2, 193.154] (0.95 confidence)</td>
<td>190.37 ns [190.206, 190.453] (0.95 confidence)</td>
<td>1.01191</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d5/rounded&#95;result</td>
<td>190.649 ns [190.449, 190.99] (0.95 confidence)</td>
<td>190.638 ns [190.299, 190.881] (0.95 confidence)</td>
<td>1.00006</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d5/strict&#95;array</td>
<td>26.5072 ns [26.3414, 26.8552] (0.95 confidence)</td>
<td>26.6981 ns [26.4818, 26.8196] (0.95 confidence)</td>
<td>0.992851</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d5/strict&#95;result</td>
<td>25.3567 ns [25.2643, 25.5294] (0.95 confidence)</td>
<td>24.2218 ns [24.0813, 24.3841] (0.95 confidence)</td>
<td>1.04685</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d2/rounded&#95;array</td>
<td>11.3864 ns [11.3708, 11.4115] (0.95 confidence)</td>
<td>11.6127 ns [11.6108, 11.6264] (0.95 confidence)</td>
<td>0.98051</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d2/rounded&#95;result</td>
<td>11.361 ns [11.3468, 11.3833] (0.95 confidence)</td>
<td>11.6101 ns [11.6039, 11.6122] (0.95 confidence)</td>
<td>0.978546</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d2/strict&#95;array</td>
<td>57.1123 ns [56.9555, 57.2923] (0.95 confidence)</td>
<td>55.8697 ns [55.6726, 56.4314] (0.95 confidence)</td>
<td>1.02224</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d2/strict&#95;result</td>
<td>12.2616 ns [12.2548, 12.2787] (0.95 confidence)</td>
<td>11.8021 ns [11.7827, 11.8295] (0.95 confidence)</td>
<td>1.03893</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d3/rounded&#95;array</td>
<td>17.6963 ns [17.6687, 17.7131] (0.95 confidence)</td>
<td>18.0737 ns [18.0541, 18.1454] (0.95 confidence)</td>
<td>0.979122</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d3/rounded&#95;result</td>
<td>17.6833 ns [17.6704, 17.6948] (0.95 confidence)</td>
<td>18.0528 ns [18.0284, 18.0783] (0.95 confidence)</td>
<td>0.979529</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d3/strict&#95;array</td>
<td>60.8535 ns [60.672, 61.0463] (0.95 confidence)</td>
<td>60.3612 ns [60.1613, 60.6588] (0.95 confidence)</td>
<td>1.00816</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d3/strict&#95;result</td>
<td>15.8612 ns [15.7409, 16.0143] (0.95 confidence)</td>
<td>15.7053 ns [15.6669, 15.7885] (0.95 confidence)</td>
<td>1.00993</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d4/rounded&#95;array</td>
<td>23.1613 ns [23.1275, 23.2249] (0.95 confidence)</td>
<td>24.0046 ns [23.9951, 24.0177] (0.95 confidence)</td>
<td>0.964869</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d4/rounded&#95;result</td>
<td>23.2459 ns [23.2053, 23.307] (0.95 confidence)</td>
<td>23.9351 ns [23.9056, 23.9676] (0.95 confidence)</td>
<td>0.971205</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d4/strict&#95;array</td>
<td>65.1003 ns [64.746, 65.3777] (0.95 confidence)</td>
<td>64.2679 ns [63.6743, 65.0856] (0.95 confidence)</td>
<td>1.01295</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d4/strict&#95;result</td>
<td>19.5113 ns [19.4377, 19.6708] (0.95 confidence)</td>
<td>19.7169 ns [19.63, 19.8108] (0.95 confidence)</td>
<td>0.989571</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d5/rounded&#95;array</td>
<td>29.5107 ns [29.4701, 29.5553] (0.95 confidence)</td>
<td>29.6969 ns [29.2643, 30.8893] (0.95 confidence)</td>
<td>0.993729</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d5/rounded&#95;result</td>
<td>29.7035 ns [29.6263, 29.7579] (0.95 confidence)</td>
<td>30.2167 ns [29.3554, 30.9773] (0.95 confidence)</td>
<td>0.983018</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d5/strict&#95;array</td>
<td>70.3099 ns [69.5562, 70.7972] (0.95 confidence)</td>
<td>73.4043 ns [72.0186, 74.9856] (0.95 confidence)</td>
<td>0.957844</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d5/strict&#95;result</td>
<td>23.9291 ns [23.8228, 23.9976] (0.95 confidence)</td>
<td>23.5875 ns [23.2422, 24.1889] (0.95 confidence)</td>
<td>1.01449</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d2/rounded&#95;array</td>
<td>173.298 ns [172.517, 173.882] (0.95 confidence)</td>
<td>171.492 ns [171.25, 171.834] (0.95 confidence)</td>
<td>1.01053</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d2/rounded&#95;result</td>
<td>171.803 ns [171.032, 173.043] (0.95 confidence)</td>
<td>168.934 ns [168.572, 169.598] (0.95 confidence)</td>
<td>1.01698</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d2/strict&#95;array</td>
<td>180.748 ns [179.813, 181.586] (0.95 confidence)</td>
<td>180.632 ns [180.116, 181.811] (0.95 confidence)</td>
<td>1.00064</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d2/strict&#95;result</td>
<td>182.077 ns [181.676, 182.397] (0.95 confidence)</td>
<td>176.662 ns [176.408, 176.957] (0.95 confidence)</td>
<td>1.03065</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d3/rounded&#95;array</td>
<td>176.129 ns [175.705, 176.688] (0.95 confidence)</td>
<td>177.78 ns [177.292, 178.303] (0.95 confidence)</td>
<td>0.990714</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d3/rounded&#95;result</td>
<td>181.542 ns [181.047, 182.302] (0.95 confidence)</td>
<td>178.527 ns [177.959, 179.06] (0.95 confidence)</td>
<td>1.01689</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d3/strict&#95;array</td>
<td>180.973 ns [180.294, 181.906] (0.95 confidence)</td>
<td>185.6 ns [184.732, 186.815] (0.95 confidence)</td>
<td>0.975074</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d3/strict&#95;result</td>
<td>183.517 ns [182.848, 183.899] (0.95 confidence)</td>
<td>184.531 ns [183.483, 185.018] (0.95 confidence)</td>
<td>0.994508</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d4/rounded&#95;array</td>
<td>187.034 ns [186.355, 187.576] (0.95 confidence)</td>
<td>186.852 ns [183.709, 191.033] (0.95 confidence)</td>
<td>1.00097</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d4/rounded&#95;result</td>
<td>187.098 ns [186.444, 188.226] (0.95 confidence)</td>
<td>187.547 ns [185.09, 189.663] (0.95 confidence)</td>
<td>0.997608</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d4/strict&#95;array</td>
<td>194.349 ns [193.559, 195.629] (0.95 confidence)</td>
<td>201.294 ns [195.759, 204.853] (0.95 confidence)</td>
<td>0.965501</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d4/strict&#95;result</td>
<td>190.847 ns [190.34, 191.371] (0.95 confidence)</td>
<td>199.876 ns [197.713, 202.704] (0.95 confidence)</td>
<td>0.954828</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d5/rounded&#95;array</td>
<td>187.251 ns [186.866, 187.507] (0.95 confidence)</td>
<td>187.636 ns [187.194, 188.124] (0.95 confidence)</td>
<td>0.997951</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d5/rounded&#95;result</td>
<td>187.85 ns [187.123, 188.751] (0.95 confidence)</td>
<td>188.974 ns [188.376, 189.524] (0.95 confidence)</td>
<td>0.994049</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d5/strict&#95;array</td>
<td>191.312 ns [190.333, 191.981] (0.95 confidence)</td>
<td>191.277 ns [190.682, 191.618] (0.95 confidence)</td>
<td>1.00018</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d5/strict&#95;result</td>
<td>186.908 ns [186.378, 187.238] (0.95 confidence)</td>
<td>188.326 ns [187.817, 188.923] (0.95 confidence)</td>
<td>0.992474</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d2/rounded&#95;array</td>
<td>207.561 ns [207.109, 208.089] (0.95 confidence)</td>
<td>198.338 ns [197.371, 199.126] (0.95 confidence)</td>
<td>1.0465</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d2/rounded&#95;result</td>
<td>207.498 ns [206.978, 208.245] (0.95 confidence)</td>
<td>196.196 ns [195.954, 196.477] (0.95 confidence)</td>
<td>1.0576</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d2/strict&#95;array</td>
<td>13613.1 ns [13493.6, 13709.2] (0.95 confidence)</td>
<td>13269.1 ns [13245.4, 13304.2] (0.95 confidence)</td>
<td>1.02593</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d2/strict&#95;result</td>
<td>213.253 ns [212.282, 214.119] (0.95 confidence)</td>
<td>201.102 ns [199.642, 202.183] (0.95 confidence)</td>
<td>1.06042</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d3/rounded&#95;array</td>
<td>215.156 ns [214.304, 216.354] (0.95 confidence)</td>
<td>209.645 ns [208.194, 210.409] (0.95 confidence)</td>
<td>1.02629</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d3/rounded&#95;result</td>
<td>211.56 ns [210.92, 212.451] (0.95 confidence)</td>
<td>208.722 ns [208.319, 209.439] (0.95 confidence)</td>
<td>1.0136</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d3/strict&#95;array</td>
<td>13528.3 ns [13498.6, 13558.6] (0.95 confidence)</td>
<td>13341.8 ns [13295.6, 13377.4] (0.95 confidence)</td>
<td>1.01398</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d3/strict&#95;result</td>
<td>210.202 ns [210.036, 210.373] (0.95 confidence)</td>
<td>207.834 ns [207.317, 208.367] (0.95 confidence)</td>
<td>1.01139</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d4/rounded&#95;array</td>
<td>214.893 ns [214.269, 215.801] (0.95 confidence)</td>
<td>227.316 ns [219.461, 232.174] (0.95 confidence)</td>
<td>0.945349</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d4/rounded&#95;result</td>
<td>213.173 ns [212.802, 213.442] (0.95 confidence)</td>
<td>226.031 ns [224.683, 229.217] (0.95 confidence)</td>
<td>0.943116</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d4/strict&#95;array</td>
<td>13408.9 ns [13376, 13424] (0.95 confidence)</td>
<td>13950.7 ns [13831.6, 14034.3] (0.95 confidence)</td>
<td>0.961161</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d4/strict&#95;result</td>
<td>212.41 ns [211.522, 213.217] (0.95 confidence)</td>
<td>219.594 ns [218.221, 222.748] (0.95 confidence)</td>
<td>0.967285</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d5/rounded&#95;array</td>
<td>224.018 ns [223.767, 225.218] (0.95 confidence)</td>
<td>230.314 ns [225.867, 236.797] (0.95 confidence)</td>
<td>0.972664</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d5/rounded&#95;result</td>
<td>224.785 ns [224.219, 225.498] (0.95 confidence)</td>
<td>231.974 ns [226.299, 235.147] (0.95 confidence)</td>
<td>0.969006</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d5/strict&#95;array</td>
<td>13420.6 ns [13378.4, 13455.9] (0.95 confidence)</td>
<td>14163.2 ns [13899.1, 14419.7] (0.95 confidence)</td>
<td>0.947572</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d5/strict&#95;result</td>
<td>221.042 ns [219.764, 221.847] (0.95 confidence)</td>
<td>220.83 ns [218.659, 226.423] (0.95 confidence)</td>
<td>1.00096</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d2/rounded&#95;array</td>
<td>169.098 ns [168.261, 170.153] (0.95 confidence)</td>
<td>159.822 ns [159.467, 160.357] (0.95 confidence)</td>
<td>1.05804</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d2/rounded&#95;result</td>
<td>168.381 ns [167.484, 169.441] (0.95 confidence)</td>
<td>158.079 ns [157.908, 158.743] (0.95 confidence)</td>
<td>1.06517</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d2/strict&#95;array</td>
<td>3873.27 ns [3813.78, 3903.7] (0.95 confidence)</td>
<td>3904.71 ns [3893.21, 3911.7] (0.95 confidence)</td>
<td>0.991947</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d2/strict&#95;result</td>
<td>165.211 ns [164.634, 165.654] (0.95 confidence)</td>
<td>161.88 ns [161.266, 163.005] (0.95 confidence)</td>
<td>1.02058</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d3/rounded&#95;array</td>
<td>172.274 ns [171.702, 173.179] (0.95 confidence)</td>
<td>172.356 ns [171.758, 173.19] (0.95 confidence)</td>
<td>0.999527</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d3/rounded&#95;result</td>
<td>171.776 ns [171.067, 172.517] (0.95 confidence)</td>
<td>171.242 ns [170.508, 171.665] (0.95 confidence)</td>
<td>1.00312</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d3/strict&#95;array</td>
<td>3913.38 ns [3852.08, 3956] (0.95 confidence)</td>
<td>3849.93 ns [3831.83, 3872.56] (0.95 confidence)</td>
<td>1.01648</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d3/strict&#95;result</td>
<td>172.427 ns [171.988, 172.805] (0.95 confidence)</td>
<td>171.194 ns [170.447, 172.287] (0.95 confidence)</td>
<td>1.0072</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d4/rounded&#95;array</td>
<td>179.348 ns [178.417, 180.287] (0.95 confidence)</td>
<td>188.622 ns [185.545, 191.259] (0.95 confidence)</td>
<td>0.950835</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d4/rounded&#95;result</td>
<td>179.644 ns [178.431, 180.57] (0.95 confidence)</td>
<td>189.389 ns [186.485, 192.444] (0.95 confidence)</td>
<td>0.948543</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d4/strict&#95;array</td>
<td>3876.04 ns [3812.61, 3968.41] (0.95 confidence)</td>
<td>4041.74 ns [3984.55, 4085.37] (0.95 confidence)</td>
<td>0.959002</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d4/strict&#95;result</td>
<td>176.297 ns [175.424, 177.153] (0.95 confidence)</td>
<td>176.422 ns [175.478, 176.86] (0.95 confidence)</td>
<td>0.999294</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d5/rounded&#95;array</td>
<td>185.93 ns [185.434, 186.731] (0.95 confidence)</td>
<td>193.206 ns [189.456, 198.605] (0.95 confidence)</td>
<td>0.962343</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d5/rounded&#95;result</td>
<td>186.398 ns [185.72, 187.177] (0.95 confidence)</td>
<td>194.451 ns [189.808, 198.011] (0.95 confidence)</td>
<td>0.958584</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d5/strict&#95;array</td>
<td>3881.88 ns [3862.9, 3905.13] (0.95 confidence)</td>
<td>4114.3 ns [4079.06, 4227.93] (0.95 confidence)</td>
<td>0.943508</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d5/strict&#95;result</td>
<td>180.445 ns [179.524, 181.355] (0.95 confidence)</td>
<td>184.372 ns [180.818, 187.682] (0.95 confidence)</td>
<td>0.978701</td>
</tr>
<tr>
<td>exact&#95;d2/det</td>
<td>0.437439 ns [0.437091, 0.437794] (0.95 confidence)</td>
<td>0.437915 ns [0.437368, 0.439197] (0.95 confidence)</td>
<td>0.998913</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;direct</td>
<td>0.437328 ns [0.437017, 0.438059] (0.95 confidence)</td>
<td>0.438738 ns [0.437975, 0.439416] (0.95 confidence)</td>
<td>0.996785</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;direct&#95;with&#95;errbound</td>
<td>1.73673 ns [1.7335, 1.74015] (0.95 confidence)</td>
<td>1.73466 ns [1.73162, 1.73822] (0.95 confidence)</td>
<td>1.0012</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;errbound</td>
<td>1.5873 ns [1.58544, 1.58975] (0.95 confidence)</td>
<td>1.58795 ns [1.58689, 1.59026] (0.95 confidence)</td>
<td>0.999588</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;exact</td>
<td>88.9472 ns [88.5107, 89.2936] (0.95 confidence)</td>
<td>87.9631 ns [87.8524, 88.1314] (0.95 confidence)</td>
<td>1.01119</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;exact&#95;f64&#95;result</td>
<td>67.7832 ns [67.671, 67.9359] (0.95 confidence)</td>
<td>67.3893 ns [67.3036, 67.4757] (0.95 confidence)</td>
<td>1.00584</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;exact&#95;rounded&#95;f64</td>
<td>70.4246 ns [70.2498, 70.5916] (0.95 confidence)</td>
<td>70.0029 ns [69.8266, 70.2306] (0.95 confidence)</td>
<td>1.00602</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;sign&#95;exact</td>
<td>3.00141 ns [2.99846, 3.00529] (0.95 confidence)</td>
<td>2.96913 ns [2.96891, 2.96963] (0.95 confidence)</td>
<td>1.01087</td>
</tr>
<tr>
<td>exact&#95;d2/solve&#95;exact</td>
<td>2968.41 ns [2958.38, 2979.71] (0.95 confidence)</td>
<td>2977.13 ns [2972.54, 2979.98] (0.95 confidence)</td>
<td>0.997072</td>
</tr>
<tr>
<td>exact&#95;d2/solve&#95;exact&#95;f64&#95;result</td>
<td>3131.09 ns [3122.66, 3137.25] (0.95 confidence)</td>
<td>3121.32 ns [3108.84, 3129.78] (0.95 confidence)</td>
<td>1.00313</td>
</tr>
<tr>
<td>exact&#95;d2/solve&#95;exact&#95;rounded&#95;f64</td>
<td>3357.99 ns [3351.85, 3363.05] (0.95 confidence)</td>
<td>3268.85 ns [3263.54, 3272.84] (0.95 confidence)</td>
<td>1.02727</td>
</tr>
<tr>
<td>exact&#95;d3/det</td>
<td>0.852501 ns [0.846696, 0.85727] (0.95 confidence)</td>
<td>0.818606 ns [0.81849, 0.818743] (0.95 confidence)</td>
<td>1.04141</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;direct</td>
<td>0.855756 ns [0.851763, 0.860336] (0.95 confidence)</td>
<td>0.818661 ns [0.818536, 0.818786] (0.95 confidence)</td>
<td>1.04531</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;direct&#95;with&#95;errbound</td>
<td>3.34905 ns [3.33871, 3.35649] (0.95 confidence)</td>
<td>3.30414 ns [3.30086, 3.30863] (0.95 confidence)</td>
<td>1.01359</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;errbound</td>
<td>3.29189 ns [3.29113, 3.29344] (0.95 confidence)</td>
<td>3.28904 ns [3.28776, 3.29022] (0.95 confidence)</td>
<td>1.00087</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;exact</td>
<td>304.715 ns [304.325, 305.598] (0.95 confidence)</td>
<td>301.352 ns [300.995, 302.138] (0.95 confidence)</td>
<td>1.01116</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;exact&#95;f64&#95;result</td>
<td>280.098 ns [279.243, 280.771] (0.95 confidence)</td>
<td>278.007 ns [277.473, 278.69] (0.95 confidence)</td>
<td>1.00752</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;exact&#95;rounded&#95;f64</td>
<td>286.97 ns [285.881, 287.995] (0.95 confidence)</td>
<td>281.086 ns [280.664, 281.538] (0.95 confidence)</td>
<td>1.02093</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;sign&#95;exact</td>
<td>4.74964 ns [4.74299, 4.75436] (0.95 confidence)</td>
<td>4.70809 ns [4.70604, 4.70932] (0.95 confidence)</td>
<td>1.00883</td>
</tr>
<tr>
<td>exact&#95;d3/solve&#95;exact</td>
<td>9134.44 ns [9113.87, 9158.76] (0.95 confidence)</td>
<td>9099.49 ns [9079.59, 9121.03] (0.95 confidence)</td>
<td>1.00384</td>
</tr>
<tr>
<td>exact&#95;d3/solve&#95;exact&#95;f64&#95;result</td>
<td>9342.79 ns [9300.5, 9371.68] (0.95 confidence)</td>
<td>9301.47 ns [9276.85, 9314.83] (0.95 confidence)</td>
<td>1.00444</td>
</tr>
<tr>
<td>exact&#95;d3/solve&#95;exact&#95;rounded&#95;f64</td>
<td>9592.31 ns [9565.83, 9620.18] (0.95 confidence)</td>
<td>9594.1 ns [9571.92, 9620.69] (0.95 confidence)</td>
<td>0.999814</td>
</tr>
<tr>
<td>exact&#95;d4/det</td>
<td>3.23577 ns [3.23154, 3.24051] (0.95 confidence)</td>
<td>4.79016 ns [4.78906, 4.79147] (0.95 confidence)</td>
<td>0.675505</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;direct</td>
<td>2.36939 ns [2.36522, 2.37139] (0.95 confidence)</td>
<td>2.34854 ns [2.34813, 2.34908] (0.95 confidence)</td>
<td>1.00888</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;direct&#95;with&#95;errbound</td>
<td>6.79501 ns [6.78543, 6.80768] (0.95 confidence)</td>
<td>6.77661 ns [6.77401, 6.77916] (0.95 confidence)</td>
<td>1.00272</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;errbound</td>
<td>6.75481 ns [6.74641, 6.77596] (0.95 confidence)</td>
<td>6.7611 ns [6.75076, 6.76679] (0.95 confidence)</td>
<td>0.999071</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;exact</td>
<td>853.169 ns [851.603, 856.24] (0.95 confidence)</td>
<td>826.173 ns [824.571, 827.294] (0.95 confidence)</td>
<td>1.03267</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;exact&#95;f64&#95;result</td>
<td>812.145 ns [809.85, 813.153] (0.95 confidence)</td>
<td>794.508 ns [793.001, 796.117] (0.95 confidence)</td>
<td>1.0222</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;exact&#95;rounded&#95;f64</td>
<td>830.59 ns [827.38, 835.316] (0.95 confidence)</td>
<td>807.691 ns [805.869, 809.605] (0.95 confidence)</td>
<td>1.02835</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;sign&#95;exact</td>
<td>8.0822 ns [8.05789, 8.10897] (0.95 confidence)</td>
<td>7.9137 ns [7.9128, 7.9166] (0.95 confidence)</td>
<td>1.02129</td>
</tr>
<tr>
<td>exact&#95;d4/solve&#95;exact</td>
<td>18954.1 ns [18748, 19046.5] (0.95 confidence)</td>
<td>19028.7 ns [18947.9, 19070.1] (0.95 confidence)</td>
<td>0.996081</td>
</tr>
<tr>
<td>exact&#95;d4/solve&#95;exact&#95;f64&#95;result</td>
<td>18781 ns [18715.4, 18815.9] (0.95 confidence)</td>
<td>19122.2 ns [19074.8, 19165.9] (0.95 confidence)</td>
<td>0.982158</td>
</tr>
<tr>
<td>exact&#95;d4/solve&#95;exact&#95;rounded&#95;f64</td>
<td>19536.7 ns [19490.1, 19632.7] (0.95 confidence)</td>
<td>19738.8 ns [19702.9, 19766.9] (0.95 confidence)</td>
<td>0.98976</td>
</tr>
<tr>
<td>exact&#95;d5/det</td>
<td>28.4751 ns [28.4482, 28.5041] (0.95 confidence)</td>
<td>31.9244 ns [31.8596, 31.9551] (0.95 confidence)</td>
<td>0.891953</td>
</tr>
<tr>
<td>exact&#95;d5/det&#95;exact</td>
<td>2864.3 ns [2862.79, 2870.03] (0.95 confidence)</td>
<td>2864.16 ns [2862.36, 2866.66] (0.95 confidence)</td>
<td>1.00005</td>
</tr>
<tr>
<td>exact&#95;d5/det&#95;exact&#95;f64&#95;result</td>
<td>2861.92 ns [2853.86, 2865.81] (0.95 confidence)</td>
<td>2842.32 ns [2838.58, 2845.56] (0.95 confidence)</td>
<td>1.0069</td>
</tr>
<tr>
<td>exact&#95;d5/det&#95;exact&#95;rounded&#95;f64</td>
<td>2890.65 ns [2884.27, 2895.8] (0.95 confidence)</td>
<td>2852.09 ns [2848.04, 2857.98] (0.95 confidence)</td>
<td>1.01352</td>
</tr>
<tr>
<td>exact&#95;d5/det&#95;sign&#95;exact</td>
<td>2829.68 ns [2825.4, 2833.79] (0.95 confidence)</td>
<td>2819.63 ns [2817.4, 2822.95] (0.95 confidence)</td>
<td>1.00356</td>
</tr>
<tr>
<td>exact&#95;d5/solve&#95;exact</td>
<td>30008.1 ns [29857.6, 30115.2] (0.95 confidence)</td>
<td>30589.5 ns [30452.7, 30665.8] (0.95 confidence)</td>
<td>0.980994</td>
</tr>
<tr>
<td>exact&#95;d5/solve&#95;exact&#95;f64&#95;result</td>
<td>30059.1 ns [29949, 30160.6] (0.95 confidence)</td>
<td>30773.9 ns [30676.5, 30866.2] (0.95 confidence)</td>
<td>0.976771</td>
</tr>
<tr>
<td>exact&#95;d5/solve&#95;exact&#95;rounded&#95;f64</td>
<td>30610.4 ns [30479.2, 30714.7] (0.95 confidence)</td>
<td>31419 ns [31294.4, 31498.3] (0.95 confidence)</td>
<td>0.974264</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;4x4/det&#95;exact</td>
<td>935.428 ns [934.389, 937.753] (0.95 confidence)</td>
<td>927.04 ns [926.151, 928.477] (0.95 confidence)</td>
<td>1.00905</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;4x4/det&#95;sign&#95;exact</td>
<td>7.95459 ns [7.95072, 7.96164] (0.95 confidence)</td>
<td>8.01669 ns [8.01509, 8.0216] (0.95 confidence)</td>
<td>0.992253</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;4x4/solve&#95;exact</td>
<td>15318.4 ns [15274.2, 15366.2] (0.95 confidence)</td>
<td>15607.3 ns [15545.1, 15640.8] (0.95 confidence)</td>
<td>0.981493</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;4x4/solve&#95;exact&#95;f64&#95;result</td>
<td>15471.9 ns [15433.8, 15514.1] (0.95 confidence)</td>
<td>15714.2 ns [15663, 15743.6] (0.95 confidence)</td>
<td>0.984582</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;4x4/solve&#95;exact&#95;rounded&#95;f64</td>
<td>15768.6 ns [15736.2, 15807.6] (0.95 confidence)</td>
<td>16120.1 ns [16077.7, 16160.4] (0.95 confidence)</td>
<td>0.978194</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;5x5/det&#95;exact</td>
<td>2952.27 ns [2946.58, 2958.05] (0.95 confidence)</td>
<td>2987.78 ns [2978.51, 3000.29] (0.95 confidence)</td>
<td>0.988115</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;5x5/det&#95;sign&#95;exact</td>
<td>2891.2 ns [2884.74, 2897.25] (0.95 confidence)</td>
<td>2965.85 ns [2959.1, 2970.87] (0.95 confidence)</td>
<td>0.974831</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;5x5/solve&#95;exact</td>
<td>24798 ns [24742.1, 24916.1] (0.95 confidence)</td>
<td>25321 ns [25245.4, 25400] (0.95 confidence)</td>
<td>0.979345</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;5x5/solve&#95;exact&#95;f64&#95;result</td>
<td>24931.9 ns [24827.4, 25014.6] (0.95 confidence)</td>
<td>25807.2 ns [25767.5, 25869.3] (0.95 confidence)</td>
<td>0.966085</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;5x5/solve&#95;exact&#95;rounded&#95;f64</td>
<td>25711.6 ns [25631.1, 25812.8] (0.95 confidence)</td>
<td>26684.4 ns [26601.1, 26779.4] (0.95 confidence)</td>
<td>0.963546</td>
</tr>
<tr>
<td>exact&#95;large&#95;entries&#95;3x3/det&#95;exact</td>
<td>318.858 ns [318.371, 319.78] (0.95 confidence)</td>
<td>322.96 ns [322.341, 324.04] (0.95 confidence)</td>
<td>0.987299</td>
</tr>
<tr>
<td>exact&#95;large&#95;entries&#95;3x3/det&#95;sign&#95;exact</td>
<td>287.332 ns [286.806, 287.909] (0.95 confidence)</td>
<td>294.472 ns [293.538, 295.622] (0.95 confidence)</td>
<td>0.975752</td>
</tr>
<tr>
<td>exact&#95;large&#95;entries&#95;3x3/solve&#95;exact</td>
<td>9098.33 ns [9077.79, 9118.89] (0.95 confidence)</td>
<td>9183.27 ns [9161.85, 9208.87] (0.95 confidence)</td>
<td>0.990751</td>
</tr>
<tr>
<td>exact&#95;large&#95;entries&#95;3x3/solve&#95;exact&#95;f64&#95;result</td>
<td>9269.46 ns [9262.39, 9287.43] (0.95 confidence)</td>
<td>9354.59 ns [9336.78, 9369.88] (0.95 confidence)</td>
<td>0.9909</td>
</tr>
<tr>
<td>exact&#95;large&#95;entries&#95;3x3/solve&#95;exact&#95;rounded&#95;f64</td>
<td>9513.19 ns [9495.38, 9534.89] (0.95 confidence)</td>
<td>9645.53 ns [9630.2, 9661.1] (0.95 confidence)</td>
<td>0.98628</td>
</tr>
<tr>
<td>exact&#95;near&#95;singular&#95;3x3/det&#95;exact</td>
<td>185.258 ns [185.079, 185.436] (0.95 confidence)</td>
<td>188.289 ns [188.045, 188.595] (0.95 confidence)</td>
<td>0.983903</td>
</tr>
<tr>
<td>exact&#95;near&#95;singular&#95;3x3/det&#95;sign&#95;exact</td>
<td>162.262 ns [162.004, 162.535] (0.95 confidence)</td>
<td>164.154 ns [163.851, 164.267] (0.95 confidence)</td>
<td>0.988479</td>
</tr>
<tr>
<td>exact&#95;near&#95;singular&#95;3x3/solve&#95;exact</td>
<td>968.507 ns [966.566, 969.988] (0.95 confidence)</td>
<td>998.773 ns [996.049, 1001.08] (0.95 confidence)</td>
<td>0.969697</td>
</tr>
<tr>
<td>exact&#95;near&#95;singular&#95;3x3/solve&#95;exact&#95;f64&#95;result</td>
<td>996.61 ns [995.05, 998.942] (0.95 confidence)</td>
<td>1026.67 ns [1024.22, 1029.49] (0.95 confidence)</td>
<td>0.970716</td>
</tr>
<tr>
<td>exact&#95;near&#95;singular&#95;3x3/solve&#95;exact&#95;rounded&#95;f64</td>
<td>988.536 ns [985.927, 990.288] (0.95 confidence)</td>
<td>1015.87 ns [1014.19, 1018.19] (0.95 confidence)</td>
<td>0.973091</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d2/det&#95;exact</td>
<td>3558.9 ns [3545.02, 3572.6] (0.95 confidence)</td>
<td>3481.11 ns [3472.55, 3488.34] (0.95 confidence)</td>
<td>1.02234</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d2/det&#95;sign&#95;exact</td>
<td>152.008 ns [151.305, 152.431] (0.95 confidence)</td>
<td>142.89 ns [142.752, 142.991] (0.95 confidence)</td>
<td>1.06381</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d2/solve&#95;exact</td>
<td>37283.2 ns [37219, 37337] (0.95 confidence)</td>
<td>35910.5 ns [35807.6, 35962.8] (0.95 confidence)</td>
<td>1.03822</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d2/solve&#95;exact&#95;f64&#95;result</td>
<td>38202.5 ns [38116.9, 38387] (0.95 confidence)</td>
<td>36446.5 ns [36392.8, 36489] (0.95 confidence)</td>
<td>1.04818</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d2/solve&#95;exact&#95;rounded&#95;f64</td>
<td>37679.8 ns [37472.7, 37833] (0.95 confidence)</td>
<td>36911.7 ns [36872.4, 36997.4] (0.95 confidence)</td>
<td>1.02081</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d3/det&#95;exact</td>
<td>9591.8 ns [9554, 9610.18] (0.95 confidence)</td>
<td>9315.15 ns [9305.72, 9320.93] (0.95 confidence)</td>
<td>1.0297</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d3/det&#95;sign&#95;exact</td>
<td>229.657 ns [229.33, 230.237] (0.95 confidence)</td>
<td>225.254 ns [225.151, 225.436] (0.95 confidence)</td>
<td>1.01955</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d3/solve&#95;exact</td>
<td>84516.5 ns [84104.8, 84685.9] (0.95 confidence)</td>
<td>81042.7 ns [80908.3, 81318.3] (0.95 confidence)</td>
<td>1.04286</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d3/solve&#95;exact&#95;f64&#95;result</td>
<td>84906.1 ns [84713.9, 85118.8] (0.95 confidence)</td>
<td>81395.5 ns [81248.6, 81513.5] (0.95 confidence)</td>
<td>1.04313</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d3/solve&#95;exact&#95;rounded&#95;f64</td>
<td>86010.2 ns [85504.5, 86239.2] (0.95 confidence)</td>
<td>81775.9 ns [81620.8, 81890.4] (0.95 confidence)</td>
<td>1.05178</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d4/det&#95;exact</td>
<td>23645.2 ns [23535.5, 23704.2] (0.95 confidence)</td>
<td>22276.9 ns [22214.1, 22357.8] (0.95 confidence)</td>
<td>1.06142</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d4/det&#95;sign&#95;exact</td>
<td>463.034 ns [462.132, 463.996] (0.95 confidence)</td>
<td>443.064 ns [442.546, 443.797] (0.95 confidence)</td>
<td>1.04507</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d4/solve&#95;exact</td>
<td>151057 ns [150622, 151291] (0.95 confidence)</td>
<td>143885 ns [143707, 144209] (0.95 confidence)</td>
<td>1.04985</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d4/solve&#95;exact&#95;f64&#95;result</td>
<td>148249 ns [147877, 148705] (0.95 confidence)</td>
<td>145343 ns [144918, 145970] (0.95 confidence)</td>
<td>1.01999</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d4/solve&#95;exact&#95;rounded&#95;f64</td>
<td>149319 ns [148849, 149942] (0.95 confidence)</td>
<td>145753 ns [145432, 146099] (0.95 confidence)</td>
<td>1.02447</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d5/det&#95;exact</td>
<td>57981.6 ns [57787.6, 58242.7] (0.95 confidence)</td>
<td>55739.8 ns [55697.1, 55864.1] (0.95 confidence)</td>
<td>1.04022</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d5/det&#95;sign&#95;exact</td>
<td>56725.3 ns [56455.7, 56904.9] (0.95 confidence)</td>
<td>54902.8 ns [54833, 55119.4] (0.95 confidence)</td>
<td>1.0332</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d5/solve&#95;exact</td>
<td>243238 ns [242808, 243970] (0.95 confidence)</td>
<td>238101 ns [237626, 238775] (0.95 confidence)</td>
<td>1.02157</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d5/solve&#95;exact&#95;f64&#95;result</td>
<td>242779 ns [241915, 243714] (0.95 confidence)</td>
<td>238721 ns [238253, 239261] (0.95 confidence)</td>
<td>1.017</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d5/solve&#95;exact&#95;rounded&#95;f64</td>
<td>245396 ns [244945, 245832] (0.95 confidence)</td>
<td>240317 ns [239839, 240867] (0.95 confidence)</td>
<td>1.02114</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/matrix&#95;access</td>
<td>0.780832 ns [0.778257, 0.782005] (0.95 confidence)</td>
<td>0.788367 ns [0.785733, 0.791456] (0.95 confidence)</td>
<td>0.990442</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/matrix&#95;from&#95;fn</td>
<td>135.779 ns [135.634, 135.959] (0.95 confidence)</td>
<td>136.404 ns [136.074, 136.518] (0.95 confidence)</td>
<td>0.995417</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/matrix&#95;from&#95;rows</td>
<td>138.738 ns [138.456, 138.961] (0.95 confidence)</td>
<td>139.927 ns [139.726, 140.164] (0.95 confidence)</td>
<td>0.9915</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/matrix&#95;into&#95;rows</td>
<td>66.9196 ns [66.8246, 67.0234] (0.95 confidence)</td>
<td>67.69 ns [67.5786, 67.8245] (0.95 confidence)</td>
<td>0.988619</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/matrix&#95;zero</td>
<td>2.73779 ns [2.73581, 2.7394] (0.95 confidence)</td>
<td>2.74523 ns [2.74285, 2.74863] (0.95 confidence)</td>
<td>0.997291</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/vector&#95;access</td>
<td>0.78334 ns [0.780676, 0.785152] (0.95 confidence)</td>
<td>0.792045 ns [0.789494, 0.793734] (0.95 confidence)</td>
<td>0.98901</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/vector&#95;from&#95;fn</td>
<td>37.4186 ns [37.351, 37.4836] (0.95 confidence)</td>
<td>37.7053 ns [37.6674, 37.7677] (0.95 confidence)</td>
<td>0.992396</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/vector&#95;into&#95;array</td>
<td>15.8952 ns [15.8814, 15.9183] (0.95 confidence)</td>
<td>16.0985 ns [16.0785, 16.1229] (0.95 confidence)</td>
<td>0.987377</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/vector&#95;new</td>
<td>36.8861 ns [36.823, 36.9343] (0.95 confidence)</td>
<td>37.0685 ns [37.0235, 37.1028] (0.95 confidence)</td>
<td>0.995078</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/vector&#95;zero</td>
<td>1.37087 ns [1.36856, 1.3725] (0.95 confidence)</td>
<td>1.3768 ns [1.37565, 1.37861] (0.95 confidence)</td>
<td>0.995691</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/matrix&#95;access</td>
<td>0.778346 ns [0.776336, 0.780757] (0.95 confidence)</td>
<td>0.788392 ns [0.785713, 0.790379] (0.95 confidence)</td>
<td>0.987258</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/matrix&#95;from&#95;fn</td>
<td>232.99 ns [232.506, 233.302] (0.95 confidence)</td>
<td>235.124 ns [234.679, 235.623] (0.95 confidence)</td>
<td>0.990925</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/matrix&#95;from&#95;rows</td>
<td>271.242 ns [270.865, 271.52] (0.95 confidence)</td>
<td>273.865 ns [273.272, 274.288] (0.95 confidence)</td>
<td>0.990421</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/matrix&#95;into&#95;rows</td>
<td>125.549 ns [125.392, 125.661] (0.95 confidence)</td>
<td>126.235 ns [126.077, 126.455] (0.95 confidence)</td>
<td>0.994563</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/matrix&#95;zero</td>
<td>6.34513 ns [6.33771, 6.35376] (0.95 confidence)</td>
<td>6.41617 ns [6.40514, 6.43077] (0.95 confidence)</td>
<td>0.988929</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/vector&#95;access</td>
<td>0.782352 ns [0.780409, 0.785227] (0.95 confidence)</td>
<td>0.790876 ns [0.788407, 0.793618] (0.95 confidence)</td>
<td>0.989222</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/vector&#95;from&#95;fn</td>
<td>52.8455 ns [52.782, 52.9274] (0.95 confidence)</td>
<td>53.1031 ns [53.03, 53.1856] (0.95 confidence)</td>
<td>0.995149</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/vector&#95;into&#95;array</td>
<td>22.9374 ns [22.9032, 22.9695] (0.95 confidence)</td>
<td>23.0728 ns [23.0522, 23.0937] (0.95 confidence)</td>
<td>0.994132</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/vector&#95;new</td>
<td>52.3116 ns [52.2706, 52.3405] (0.95 confidence)</td>
<td>52.7036 ns [52.6178, 52.8055] (0.95 confidence)</td>
<td>0.992562</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/vector&#95;zero</td>
<td>2.04744 ns [2.04568, 2.05043] (0.95 confidence)</td>
<td>2.05773 ns [2.05483, 2.05971] (0.95 confidence)</td>
<td>0.994997</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/matrix&#95;access</td>
<td>0.780269 ns [0.778648, 0.783117] (0.95 confidence)</td>
<td>0.786655 ns [0.78466, 0.788565] (0.95 confidence)</td>
<td>0.991881</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/matrix&#95;from&#95;fn</td>
<td>289.753 ns [288.907, 290.376] (0.95 confidence)</td>
<td>293.205 ns [292.691, 293.958] (0.95 confidence)</td>
<td>0.988225</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/matrix&#95;from&#95;rows</td>
<td>351.576 ns [351.308, 351.921] (0.95 confidence)</td>
<td>353.549 ns [352.854, 354.68] (0.95 confidence)</td>
<td>0.994418</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/matrix&#95;into&#95;rows</td>
<td>133.474 ns [133.168, 133.741] (0.95 confidence)</td>
<td>134.211 ns [133.941, 134.353] (0.95 confidence)</td>
<td>0.994503</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/matrix&#95;zero</td>
<td>15.1976 ns [15.0702, 15.3092] (0.95 confidence)</td>
<td>15.6226 ns [15.576, 15.6779] (0.95 confidence)</td>
<td>0.972793</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/vector&#95;access</td>
<td>0.780073 ns [0.777877, 0.781819] (0.95 confidence)</td>
<td>0.788752 ns [0.786615, 0.792164] (0.95 confidence)</td>
<td>0.988997</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/vector&#95;from&#95;fn</td>
<td>67.4894 ns [67.4223, 67.683] (0.95 confidence)</td>
<td>68.138 ns [68.0536, 68.3035] (0.95 confidence)</td>
<td>0.990482</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/vector&#95;into&#95;array</td>
<td>30.4201 ns [30.3895, 30.4767] (0.95 confidence)</td>
<td>30.6726 ns [30.6218, 30.7044] (0.95 confidence)</td>
<td>0.991768</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/vector&#95;new</td>
<td>66.4669 ns [66.4042, 66.5798] (0.95 confidence)</td>
<td>67.0391 ns [66.932, 67.1278] (0.95 confidence)</td>
<td>0.991464</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/vector&#95;zero</td>
<td>3.32833 ns [3.31854, 3.34185] (0.95 confidence)</td>
<td>3.34941 ns [3.33965, 3.36212] (0.95 confidence)</td>
<td>0.993706</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/matrix&#95;access</td>
<td>0.78194 ns [0.779246, 0.783857] (0.95 confidence)</td>
<td>0.802427 ns [0.798151, 0.805735] (0.95 confidence)</td>
<td>0.974468</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/matrix&#95;from&#95;fn</td>
<td>441.077 ns [438.322, 442.788] (0.95 confidence)</td>
<td>440.222 ns [438.209, 442.183] (0.95 confidence)</td>
<td>1.00194</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/matrix&#95;from&#95;rows</td>
<td>630.548 ns [629.286, 631.824] (0.95 confidence)</td>
<td>631.724 ns [630.162, 633.295] (0.95 confidence)</td>
<td>0.998138</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/matrix&#95;into&#95;rows</td>
<td>306.895 ns [303.62, 307.741] (0.95 confidence)</td>
<td>307.967 ns [306.336, 308.567] (0.95 confidence)</td>
<td>0.996519</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/matrix&#95;zero</td>
<td>22.3146 ns [22.2734, 22.335] (0.95 confidence)</td>
<td>22.8993 ns [22.8095, 22.9418] (0.95 confidence)</td>
<td>0.974465</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/vector&#95;access</td>
<td>0.777722 ns [0.775793, 0.780029] (0.95 confidence)</td>
<td>0.802511 ns [0.798397, 0.806196] (0.95 confidence)</td>
<td>0.969111</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/vector&#95;from&#95;fn</td>
<td>82.8907 ns [82.7175, 83.0472] (0.95 confidence)</td>
<td>83.5793 ns [83.3502, 83.7711] (0.95 confidence)</td>
<td>0.991761</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/vector&#95;into&#95;array</td>
<td>38.8386 ns [38.6853, 38.9403] (0.95 confidence)</td>
<td>39.0188 ns [38.9477, 39.0935] (0.95 confidence)</td>
<td>0.995383</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/vector&#95;new</td>
<td>83.1168 ns [82.954, 83.1925] (0.95 confidence)</td>
<td>83.3766 ns [83.284, 83.531] (0.95 confidence)</td>
<td>0.996884</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/vector&#95;zero</td>
<td>4.338 ns [4.29246, 4.36616] (0.95 confidence)</td>
<td>4.30881 ns [4.28333, 4.34674] (0.95 confidence)</td>
<td>1.00678</td>
</tr>
<tr>
<td>public&#95;exact&#95;scalar/rounded</td>
<td>5.58263 ns [5.57564, 5.58756] (0.95 confidence)</td>
<td>5.60008 ns [5.59298, 5.60954] (0.95 confidence)</td>
<td>0.996885</td>
</tr>
<tr>
<td>public&#95;exact&#95;scalar/strict</td>
<td>4.2923 ns [4.24732, 4.33184] (0.95 confidence)</td>
<td>4.49576 ns [4.48863, 4.50274] (0.95 confidence)</td>
<td>0.954744</td>
</tr>
<tr>
<td>rational&#95;input&#95;d2/det&#95;row&#95;cleared&#95;bareiss</td>
<td>178.509 ns [178.111, 178.962] (0.95 confidence)</td>
<td>179.707 ns [179.038, 180.139] (0.95 confidence)</td>
<td>0.993334</td>
</tr>
<tr>
<td>rational&#95;input&#95;d2/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>77.1861 ns [77.0208, 77.3878] (0.95 confidence)</td>
<td>77.6917 ns [77.5901, 77.9512] (0.95 confidence)</td>
<td>0.993492</td>
</tr>
<tr>
<td>rational&#95;input&#95;d2/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>624.118 ns [623.509, 625.37] (0.95 confidence)</td>
<td>626.125 ns [625.217, 626.859] (0.95 confidence)</td>
<td>0.996795</td>
</tr>
<tr>
<td>rational&#95;input&#95;d3/det&#95;row&#95;cleared&#95;bareiss</td>
<td>596.63 ns [595.961, 597.726] (0.95 confidence)</td>
<td>599.424 ns [597.655, 600.388] (0.95 confidence)</td>
<td>0.995339</td>
</tr>
<tr>
<td>rational&#95;input&#95;d3/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>389.337 ns [388.505, 390.429] (0.95 confidence)</td>
<td>395.54 ns [394.798, 396.278] (0.95 confidence)</td>
<td>0.984317</td>
</tr>
<tr>
<td>rational&#95;input&#95;d3/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>1261.79 ns [1260.21, 1263.41] (0.95 confidence)</td>
<td>1268.86 ns [1266.89, 1271.24] (0.95 confidence)</td>
<td>0.994432</td>
</tr>
<tr>
<td>rational&#95;input&#95;d4/det&#95;row&#95;cleared&#95;bareiss</td>
<td>1229.18 ns [1225.63, 1232.71] (0.95 confidence)</td>
<td>1240.48 ns [1239.15, 1242.05] (0.95 confidence)</td>
<td>0.990889</td>
</tr>
<tr>
<td>rational&#95;input&#95;d4/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>821.615 ns [820.343, 823.159] (0.95 confidence)</td>
<td>829.51 ns [828.561, 831.182] (0.95 confidence)</td>
<td>0.990481</td>
</tr>
<tr>
<td>rational&#95;input&#95;d4/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>2226.78 ns [2223.08, 2230.3] (0.95 confidence)</td>
<td>2265.15 ns [2261.57, 2269.85] (0.95 confidence)</td>
<td>0.98306</td>
</tr>
<tr>
<td>rational&#95;input&#95;d5/det&#95;row&#95;cleared&#95;bareiss</td>
<td>2103.79 ns [2099.35, 2108.25] (0.95 confidence)</td>
<td>2092.96 ns [2090.32, 2099.85] (0.95 confidence)</td>
<td>1.00518</td>
</tr>
<tr>
<td>rational&#95;input&#95;d5/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>1576.35 ns [1572.2, 1578.36] (0.95 confidence)</td>
<td>1581.1 ns [1574.36, 1586.78] (0.95 confidence)</td>
<td>0.996992</td>
</tr>
<tr>
<td>rational&#95;input&#95;d5/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>3369.27 ns [3364.07, 3373.38] (0.95 confidence)</td>
<td>3384.32 ns [3376.65, 3388.34] (0.95 confidence)</td>
<td>0.995555</td>
</tr>
<tr>
<td>rational&#95;input&#95;d6/det&#95;row&#95;cleared&#95;bareiss</td>
<td>3513.61 ns [3508.18, 3522.52] (0.95 confidence)</td>
<td>3524.89 ns [3521.09, 3536.05] (0.95 confidence)</td>
<td>0.996798</td>
</tr>
<tr>
<td>rational&#95;input&#95;d6/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>2848.03 ns [2841.25, 2862.54] (0.95 confidence)</td>
<td>2877.88 ns [2868.98, 2885.91] (0.95 confidence)</td>
<td>0.989625</td>
</tr>
<tr>
<td>rational&#95;input&#95;d6/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>5737.9 ns [5732.65, 5743.77] (0.95 confidence)</td>
<td>5840.38 ns [5806.49, 5859.38] (0.95 confidence)</td>
<td>0.982454</td>
</tr>
<tr>
<td>rational&#95;input&#95;d7/det&#95;row&#95;cleared&#95;bareiss</td>
<td>5469.42 ns [5464.03, 5477.95] (0.95 confidence)</td>
<td>5551.86 ns [5537.11, 5560.48] (0.95 confidence)</td>
<td>0.985151</td>
</tr>
<tr>
<td>rational&#95;input&#95;d7/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>4677.58 ns [4672.49, 4689.9] (0.95 confidence)</td>
<td>4726.28 ns [4713.78, 4735.36] (0.95 confidence)</td>
<td>0.989696</td>
</tr>
<tr>
<td>rational&#95;input&#95;d7/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>9901.57 ns [9877.45, 9921.59] (0.95 confidence)</td>
<td>10083 ns [10054.1, 10111.3] (0.95 confidence)</td>
<td>0.982006</td>
</tr>
<tr>
<td>rational&#95;input&#95;d8/det&#95;row&#95;cleared&#95;bareiss</td>
<td>8137.11 ns [8126.3, 8169.05] (0.95 confidence)</td>
<td>8175.65 ns [8167.9, 8183.22] (0.95 confidence)</td>
<td>0.995286</td>
</tr>
<tr>
<td>rational&#95;input&#95;d8/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>7332.59 ns [7320.56, 7339.13] (0.95 confidence)</td>
<td>7381.15 ns [7370.04, 7390.56] (0.95 confidence)</td>
<td>0.993421</td>
</tr>
<tr>
<td>rational&#95;input&#95;d8/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>15123.6 ns [15103.4, 15147.4] (0.95 confidence)</td>
<td>15338.1 ns [15297.5, 15407.9] (0.95 confidence)</td>
<td>0.986015</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d2/det&#95;row&#95;cleared&#95;bareiss</td>
<td>42193.3 ns [42148.4, 42283.3] (0.95 confidence)</td>
<td>43373.1 ns [43319.9, 43488.7] (0.95 confidence)</td>
<td>0.972799</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d2/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>953.357 ns [952.775, 954.097] (0.95 confidence)</td>
<td>982.735 ns [979.064, 986.42] (0.95 confidence)</td>
<td>0.970106</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d2/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>3036.93 ns [3030.89, 3039.19] (0.95 confidence)</td>
<td>3201.15 ns [3190.81, 3214.89] (0.95 confidence)</td>
<td>0.948697</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d3/det&#95;row&#95;cleared&#95;bareiss</td>
<td>69081.3 ns [69004.4, 69143.5] (0.95 confidence)</td>
<td>70822.5 ns [70736.9, 71072.7] (0.95 confidence)</td>
<td>0.975414</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d3/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>4194.1 ns [4187.09, 4196.47] (0.95 confidence)</td>
<td>4270.85 ns [4266.26, 4277.33] (0.95 confidence)</td>
<td>0.98203</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d3/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>11415.5 ns [11412.3, 11434.2] (0.95 confidence)</td>
<td>11583.6 ns [11573.9, 11599.2] (0.95 confidence)</td>
<td>0.985485</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d4/det&#95;row&#95;cleared&#95;bareiss</td>
<td>110073 ns [109943, 110159] (0.95 confidence)</td>
<td>111981 ns [111689, 112140] (0.95 confidence)</td>
<td>0.982964</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d4/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>9390.15 ns [9372.25, 9401.79] (0.95 confidence)</td>
<td>9556.83 ns [9548.5, 9583.45] (0.95 confidence)</td>
<td>0.982559</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d4/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>33628.4 ns [33589.7, 33706.2] (0.95 confidence)</td>
<td>34445.6 ns [34388, 34587.5] (0.95 confidence)</td>
<td>0.976277</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d5/det&#95;row&#95;cleared&#95;bareiss</td>
<td>187418 ns [186999, 187628] (0.95 confidence)</td>
<td>188437 ns [188291, 188560] (0.95 confidence)</td>
<td>0.994592</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d5/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>43111.7 ns [43057.9, 43159.7] (0.95 confidence)</td>
<td>43735 ns [43681.6, 43778.6] (0.95 confidence)</td>
<td>0.985748</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d5/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>70425 ns [70347.3, 70574.2] (0.95 confidence)</td>
<td>70726.5 ns [70691.9, 70775.4] (0.95 confidence)</td>
<td>0.995737</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d6/det&#95;row&#95;cleared&#95;bareiss</td>
<td>300467 ns [300344, 300773] (0.95 confidence)</td>
<td>308258 ns [308125, 308563] (0.95 confidence)</td>
<td>0.974728</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d6/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>104553 ns [104515, 104617] (0.95 confidence)</td>
<td>106281 ns [106204, 106460] (0.95 confidence)</td>
<td>0.983744</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d6/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>162004 ns [161863, 162137] (0.95 confidence)</td>
<td>166249 ns [166190, 166329] (0.95 confidence)</td>
<td>0.974466</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d7/det&#95;row&#95;cleared&#95;bareiss</td>
<td>501825 ns [501673, 502316] (0.95 confidence)</td>
<td>514813 ns [514395, 515500] (0.95 confidence)</td>
<td>0.974771</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d7/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>233429 ns [233218, 233760] (0.95 confidence)</td>
<td>239566 ns [239474, 239702] (0.95 confidence)</td>
<td>0.97438</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d7/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>340392 ns [340181, 340551] (0.95 confidence)</td>
<td>350433 ns [349862, 351007] (0.95 confidence)</td>
<td>0.971348</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d8/det&#95;row&#95;cleared&#95;bareiss</td>
<td>822776 ns [822603, 823553] (0.95 confidence)</td>
<td>834807 ns [834543, 835134] (0.95 confidence)</td>
<td>0.985587</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d8/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>468867 ns [468613, 469098] (0.95 confidence)</td>
<td>476072 ns [475893, 476551] (0.95 confidence)</td>
<td>0.984867</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d8/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>651457 ns [651232, 652209] (0.95 confidence)</td>
<td>670533 ns [663211, 672054] (0.95 confidence)</td>
<td>0.971551</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d2/det&#95;row&#95;cleared&#95;bareiss</td>
<td>9061.41 ns [9026.32, 9105.22] (0.95 confidence)</td>
<td>9222.18 ns [9205.34, 9257.58] (0.95 confidence)</td>
<td>0.982567</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d2/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>687.722 ns [669.611, 720.967] (0.95 confidence)</td>
<td>659.164 ns [657.674, 662.486] (0.95 confidence)</td>
<td>1.04333</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d2/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>2058.7 ns [2055.08, 2067.03] (0.95 confidence)</td>
<td>2111.84 ns [2100.29, 2122.39] (0.95 confidence)</td>
<td>0.97484</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d3/det&#95;row&#95;cleared&#95;bareiss</td>
<td>13431.7 ns [13395.4, 13454.2] (0.95 confidence)</td>
<td>13698.6 ns [13651.1, 13739.3] (0.95 confidence)</td>
<td>0.980513</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d3/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>2755.2 ns [2746.52, 2766.39] (0.95 confidence)</td>
<td>2792.73 ns [2785.81, 2801.5] (0.95 confidence)</td>
<td>0.986559</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d3/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>5777.77 ns [5765.61, 5793.28] (0.95 confidence)</td>
<td>5895.42 ns [5875.02, 5909.92] (0.95 confidence)</td>
<td>0.980043</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d4/det&#95;row&#95;cleared&#95;bareiss</td>
<td>21672.4 ns [21634.3, 21712.4] (0.95 confidence)</td>
<td>22203.4 ns [22177, 22233.8] (0.95 confidence)</td>
<td>0.976086</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d4/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>5891.07 ns [5870.96, 5918.87] (0.95 confidence)</td>
<td>6077.23 ns [6051.01, 6092.83] (0.95 confidence)</td>
<td>0.969368</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d4/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>12744.5 ns [12727.9, 12760.7] (0.95 confidence)</td>
<td>13010.2 ns [12957.5, 13040.3] (0.95 confidence)</td>
<td>0.979577</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d5/det&#95;row&#95;cleared&#95;bareiss</td>
<td>33557.2 ns [33536.8, 33585.9] (0.95 confidence)</td>
<td>34357.4 ns [34331.2, 34386.1] (0.95 confidence)</td>
<td>0.97671</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d5/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>13460.7 ns [13413.1, 13485.5] (0.95 confidence)</td>
<td>13500.9 ns [13454.5, 13535.5] (0.95 confidence)</td>
<td>0.997023</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d5/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>20917.4 ns [20890.1, 20987.1] (0.95 confidence)</td>
<td>21092.1 ns [21050.8, 21185.3] (0.95 confidence)</td>
<td>0.99172</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d6/det&#95;row&#95;cleared&#95;bareiss</td>
<td>52414.1 ns [52323, 52508.3] (0.95 confidence)</td>
<td>53006.2 ns [52919.6, 53071.3] (0.95 confidence)</td>
<td>0.98883</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d6/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>27061.9 ns [26817.9, 27604.9] (0.95 confidence)</td>
<td>26473.5 ns [26442.1, 26519.4] (0.95 confidence)</td>
<td>1.02222</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d6/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>39541.3 ns [39467.9, 39584.7] (0.95 confidence)</td>
<td>40094.4 ns [40013.8, 40371.6] (0.95 confidence)</td>
<td>0.986204</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d7/det&#95;row&#95;cleared&#95;bareiss</td>
<td>83309.7 ns [83201.7, 83415.2] (0.95 confidence)</td>
<td>84414.4 ns [84324.8, 84555.2] (0.95 confidence)</td>
<td>0.986914</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d7/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>50119.1 ns [50072.4, 50178.9] (0.95 confidence)</td>
<td>50993.8 ns [50904.8, 51386.5] (0.95 confidence)</td>
<td>0.982847</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d7/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>72024.1 ns [71952.1, 72104] (0.95 confidence)</td>
<td>73180.2 ns [73126.7, 73233.1] (0.95 confidence)</td>
<td>0.984202</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d8/det&#95;row&#95;cleared&#95;bareiss</td>
<td>131370 ns [131281, 131498] (0.95 confidence)</td>
<td>134131 ns [133645, 134372] (0.95 confidence)</td>
<td>0.979419</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d8/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>90378.5 ns [90266.7, 90486.5] (0.95 confidence)</td>
<td>91988 ns [91809.4, 92168.8] (0.95 confidence)</td>
<td>0.982503</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d8/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>127141 ns [127042, 127281] (0.95 confidence)</td>
<td>130524 ns [130073, 131022] (0.95 confidence)</td>
<td>0.974075</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d2/adapter</td>
<td>1064.02 ns [1038.76, 1073.6] (0.95 confidence)</td>
<td>1058.05 ns [1028.63, 1068.61] (0.95 confidence)</td>
<td>1.00564</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d2/construct</td>
<td>292.824 ns [291.458, 294.008] (0.95 confidence)</td>
<td>292.412 ns [291.197, 293.427] (0.95 confidence)</td>
<td>1.00141</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d2/prepared</td>
<td>650.816 ns [650.077, 652.436] (0.95 confidence)</td>
<td>649.654 ns [649.192, 651.851] (0.95 confidence)</td>
<td>1.00179</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d3/adapter</td>
<td>1928.96 ns [1926.06, 1931.54] (0.95 confidence)</td>
<td>1942.28 ns [1938.11, 1947.25] (0.95 confidence)</td>
<td>0.993141</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d3/construct</td>
<td>650.528 ns [642.886, 655.029] (0.95 confidence)</td>
<td>640.456 ns [639.736, 642.517] (0.95 confidence)</td>
<td>1.01573</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d3/prepared</td>
<td>1152.3 ns [1151.15, 1154.16] (0.95 confidence)</td>
<td>1154.24 ns [1152.99, 1156.86] (0.95 confidence)</td>
<td>0.998313</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d4/adapter</td>
<td>3477.02 ns [3472.42, 3481.63] (0.95 confidence)</td>
<td>3367.17 ns [3357.86, 3372.32] (0.95 confidence)</td>
<td>1.03263</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d4/construct</td>
<td>1101.28 ns [1098.82, 1110.3] (0.95 confidence)</td>
<td>1081.24 ns [1077.38, 1090.35] (0.95 confidence)</td>
<td>1.01853</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d4/prepared</td>
<td>2158.45 ns [2144.04, 2172.3] (0.95 confidence)</td>
<td>2103.06 ns [2096.71, 2111.55] (0.95 confidence)</td>
<td>1.02634</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d5/adapter</td>
<td>5635.19 ns [5629.11, 5645.92] (0.95 confidence)</td>
<td>5598.15 ns [5583.65, 5630.29] (0.95 confidence)</td>
<td>1.00662</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d5/construct</td>
<td>2067.26 ns [2055.76, 2075.43] (0.95 confidence)</td>
<td>1987.72 ns [1965.61, 2004.94] (0.95 confidence)</td>
<td>1.04002</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d5/prepared</td>
<td>3289.73 ns [3282.17, 3298.58] (0.95 confidence)</td>
<td>3166.96 ns [3152.21, 3178.83] (0.95 confidence)</td>
<td>1.03876</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d6/adapter</td>
<td>7753.41 ns [7748.36, 7766.7] (0.95 confidence)</td>
<td>7791.1 ns [7767.85, 7803.11] (0.95 confidence)</td>
<td>0.995163</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d6/construct</td>
<td>2599.51 ns [2585.06, 2620.48] (0.95 confidence)</td>
<td>2683.82 ns [2654.49, 2694.41] (0.95 confidence)</td>
<td>0.968588</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d6/prepared</td>
<td>4752.31 ns [4745.79, 4757.17] (0.95 confidence)</td>
<td>4751.62 ns [4733.45, 4767.65] (0.95 confidence)</td>
<td>1.00015</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d7/adapter</td>
<td>10471.6 ns [10453.6, 10501.5] (0.95 confidence)</td>
<td>10504 ns [10480.5, 10519.9] (0.95 confidence)</td>
<td>0.996912</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d7/construct</td>
<td>3483.91 ns [3473.23, 3529.64] (0.95 confidence)</td>
<td>3559.33 ns [3521.7, 3580.27] (0.95 confidence)</td>
<td>0.978811</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d7/prepared</td>
<td>6708.72 ns [6677.73, 6746.45] (0.95 confidence)</td>
<td>6525.52 ns [6510.47, 6550.38] (0.95 confidence)</td>
<td>1.02807</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d8/adapter</td>
<td>13633.3 ns [13620.3, 13645.5] (0.95 confidence)</td>
<td>13577 ns [13544.4, 13596.6] (0.95 confidence)</td>
<td>1.00415</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d8/construct</td>
<td>4128.93 ns [4096.32, 4154.12] (0.95 confidence)</td>
<td>4164.63 ns [4142.48, 4183.26] (0.95 confidence)</td>
<td>0.991427</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d8/prepared</td>
<td>9320.48 ns [9257.02, 9361.85] (0.95 confidence)</td>
<td>9120.64 ns [9089.53, 9155.8] (0.95 confidence)</td>
<td>1.02191</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d2/adapter</td>
<td>1046.6 ns [1043.37, 1055.89] (0.95 confidence)</td>
<td>1051.03 ns [1047.01, 1054.99] (0.95 confidence)</td>
<td>0.995788</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d2/construct</td>
<td>329.525 ns [328.613, 330.417] (0.95 confidence)</td>
<td>329.244 ns [328.557, 329.931] (0.95 confidence)</td>
<td>1.00085</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d2/prepared</td>
<td>632.808 ns [631.732, 633.989] (0.95 confidence)</td>
<td>629.503 ns [628.843, 630.768] (0.95 confidence)</td>
<td>1.00525</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d3/adapter</td>
<td>2107.53 ns [2105.27, 2113.93] (0.95 confidence)</td>
<td>2111.07 ns [2108.2, 2115.28] (0.95 confidence)</td>
<td>0.998325</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d3/construct</td>
<td>683.972 ns [679.34, 685.665] (0.95 confidence)</td>
<td>685.014 ns [681.673, 687.513] (0.95 confidence)</td>
<td>0.998478</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d3/prepared</td>
<td>1298 ns [1295.92, 1299.66] (0.95 confidence)</td>
<td>1305.39 ns [1304.02, 1306.77] (0.95 confidence)</td>
<td>0.994338</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d4/adapter</td>
<td>3765.49 ns [3660.04, 3798.48] (0.95 confidence)</td>
<td>3532.35 ns [3530.92, 3536.63] (0.95 confidence)</td>
<td>1.066</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d4/construct</td>
<td>1194.54 ns [1189.1, 1202.12] (0.95 confidence)</td>
<td>1158.61 ns [1151.23, 1165.73] (0.95 confidence)</td>
<td>1.03101</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d4/prepared</td>
<td>2278.93 ns [2274.97, 2283.62] (0.95 confidence)</td>
<td>2204.52 ns [2202.6, 2206.8] (0.95 confidence)</td>
<td>1.03375</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d5/adapter</td>
<td>5705.12 ns [5697.24, 5714.99] (0.95 confidence)</td>
<td>5631.77 ns [5627.15, 5641.4] (0.95 confidence)</td>
<td>1.01302</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d5/construct</td>
<td>1995.95 ns [1987.99, 2005.22] (0.95 confidence)</td>
<td>1953.02 ns [1944.57, 1959.36] (0.95 confidence)</td>
<td>1.02198</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d5/prepared</td>
<td>3477.34 ns [3472.05, 3481.42] (0.95 confidence)</td>
<td>3459.86 ns [3452.09, 3467.15] (0.95 confidence)</td>
<td>1.00505</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d6/adapter</td>
<td>8101.41 ns [8089.87, 8112.83] (0.95 confidence)</td>
<td>8050.17 ns [8038.26, 8066.93] (0.95 confidence)</td>
<td>1.00637</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d6/construct</td>
<td>2895.71 ns [2882.71, 2905.85] (0.95 confidence)</td>
<td>2814.75 ns [2804.71, 2836.6] (0.95 confidence)</td>
<td>1.02876</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d6/prepared</td>
<td>5051.48 ns [5036.95, 5058.32] (0.95 confidence)</td>
<td>4945.59 ns [4938.69, 4951.92] (0.95 confidence)</td>
<td>1.02141</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d7/adapter</td>
<td>11664.8 ns [11654.2, 11682.2] (0.95 confidence)</td>
<td>11855.1 ns [11833.8, 11888.2] (0.95 confidence)</td>
<td>0.983946</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d7/construct</td>
<td>3941.64 ns [3927.9, 3954.59] (0.95 confidence)</td>
<td>4004.17 ns [3990.22, 4029.85] (0.95 confidence)</td>
<td>0.984383</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d7/prepared</td>
<td>7377.7 ns [7364.2, 7390.91] (0.95 confidence)</td>
<td>7485.7 ns [7457.5, 7530.3] (0.95 confidence)</td>
<td>0.985572</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d8/adapter</td>
<td>16197.5 ns [16173.6, 16216.3] (0.95 confidence)</td>
<td>16175.5 ns [16150, 16220.9] (0.95 confidence)</td>
<td>1.00136</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d8/construct</td>
<td>5199.58 ns [5173.11, 5237.39] (0.95 confidence)</td>
<td>5239.43 ns [5223.7, 5265.17] (0.95 confidence)</td>
<td>0.992394</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d8/prepared</td>
<td>10558.6 ns [10541.7, 10571.4] (0.95 confidence)</td>
<td>10587.2 ns [10577, 10603.4] (0.95 confidence)</td>
<td>0.997295</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d2/adapter</td>
<td>431.379 ns [430.733, 432.562] (0.95 confidence)</td>
<td>433.641 ns [432.746, 435.007] (0.95 confidence)</td>
<td>0.994783</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d2/construct</td>
<td>225.927 ns [225.535, 226.236] (0.95 confidence)</td>
<td>228.958 ns [227.923, 231.34] (0.95 confidence)</td>
<td>0.986759</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d2/prepared</td>
<td>153.88 ns [153.53, 154.152] (0.95 confidence)</td>
<td>157.691 ns [157.168, 158.155] (0.95 confidence)</td>
<td>0.975835</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d3/adapter</td>
<td>748.283 ns [747.464, 750.031] (0.95 confidence)</td>
<td>739.706 ns [738.707, 741.026] (0.95 confidence)</td>
<td>1.01159</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d3/construct</td>
<td>412.119 ns [410.111, 413.255] (0.95 confidence)</td>
<td>411.852 ns [410.707, 412.652] (0.95 confidence)</td>
<td>1.00065</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d3/prepared</td>
<td>250.897 ns [249.17, 251.641] (0.95 confidence)</td>
<td>250.011 ns [249.596, 250.346] (0.95 confidence)</td>
<td>1.00355</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d4/adapter</td>
<td>1037.84 ns [1036.57, 1039.7] (0.95 confidence)</td>
<td>1030.63 ns [1028.47, 1032.17] (0.95 confidence)</td>
<td>1.007</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d4/construct</td>
<td>617.709 ns [613.227, 621.469] (0.95 confidence)</td>
<td>616.555 ns [612.954, 619.926] (0.95 confidence)</td>
<td>1.00187</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d4/prepared</td>
<td>329.513 ns [329.108, 329.736] (0.95 confidence)</td>
<td>325.516 ns [324.936, 325.907] (0.95 confidence)</td>
<td>1.01228</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d5/adapter</td>
<td>1432.38 ns [1428.21, 1434.03] (0.95 confidence)</td>
<td>1394.73 ns [1389.96, 1398.1] (0.95 confidence)</td>
<td>1.02699</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d5/construct</td>
<td>881.087 ns [876.748, 883.731] (0.95 confidence)</td>
<td>859.495 ns [853.807, 862.759] (0.95 confidence)</td>
<td>1.02512</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d5/prepared</td>
<td>439.884 ns [439.045, 440.708] (0.95 confidence)</td>
<td>422.242 ns [421.721, 422.5] (0.95 confidence)</td>
<td>1.04178</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d6/adapter</td>
<td>1829.97 ns [1823.51, 1834.06] (0.95 confidence)</td>
<td>1756.1 ns [1751.2, 1760.05] (0.95 confidence)</td>
<td>1.04206</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d6/construct</td>
<td>1127.94 ns [1123.53, 1132.57] (0.95 confidence)</td>
<td>1119.07 ns [1110.17, 1124.76] (0.95 confidence)</td>
<td>1.00792</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d6/prepared</td>
<td>524.944 ns [522.953, 527.245] (0.95 confidence)</td>
<td>512.302 ns [511.238, 513.105] (0.95 confidence)</td>
<td>1.02468</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d7/adapter</td>
<td>2203.03 ns [2195.39, 2207.41] (0.95 confidence)</td>
<td>2215.47 ns [2207.17, 2222.25] (0.95 confidence)</td>
<td>0.994383</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d7/construct</td>
<td>1463.98 ns [1462.31, 1471.25] (0.95 confidence)</td>
<td>1440.89 ns [1432.61, 1448.97] (0.95 confidence)</td>
<td>1.01603</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d7/prepared</td>
<td>612.081 ns [610.141, 612.918] (0.95 confidence)</td>
<td>607.441 ns [605.846, 609.242] (0.95 confidence)</td>
<td>1.00764</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d8/adapter</td>
<td>2514.61 ns [2510.51, 2519.18] (0.95 confidence)</td>
<td>2519.16 ns [2511.45, 2521.81] (0.95 confidence)</td>
<td>0.998196</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d8/construct</td>
<td>1736.03 ns [1731.98, 1741.16] (0.95 confidence)</td>
<td>1742.36 ns [1734.8, 1753.73] (0.95 confidence)</td>
<td>0.996369</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d8/prepared</td>
<td>646.07 ns [645.333, 647.266] (0.95 confidence)</td>
<td>650.073 ns [649.074, 651.696] (0.95 confidence)</td>
<td>0.993843</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d2/adapter</td>
<td>573.815 ns [572.82, 575.304] (0.95 confidence)</td>
<td>573.281 ns [572.02, 575.531] (0.95 confidence)</td>
<td>1.00093</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d2/construct</td>
<td>308.609 ns [307.032, 309.824] (0.95 confidence)</td>
<td>313.503 ns [312.678, 315.754] (0.95 confidence)</td>
<td>0.984391</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d2/prepared</td>
<td>190.415 ns [189.945, 190.708] (0.95 confidence)</td>
<td>194.408 ns [193.752, 195.372] (0.95 confidence)</td>
<td>0.979458</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d3/adapter</td>
<td>916.19 ns [914.101, 917.605] (0.95 confidence)</td>
<td>912.585 ns [911.738, 914.592] (0.95 confidence)</td>
<td>1.00395</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d3/construct</td>
<td>535.153 ns [533.654, 537.534] (0.95 confidence)</td>
<td>538.218 ns [537.38, 539.936] (0.95 confidence)</td>
<td>0.994304</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d3/prepared</td>
<td>277.299 ns [276.877, 277.652] (0.95 confidence)</td>
<td>277.685 ns [277.188, 278.208] (0.95 confidence)</td>
<td>0.998611</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d4/adapter</td>
<td>1275.79 ns [1274.12, 1278.9] (0.95 confidence)</td>
<td>1274.09 ns [1270.31, 1276.6] (0.95 confidence)</td>
<td>1.00133</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d4/construct</td>
<td>788.578 ns [781.236, 793.558] (0.95 confidence)</td>
<td>797.401 ns [790.052, 805.062] (0.95 confidence)</td>
<td>0.988936</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d4/prepared</td>
<td>372.814 ns [372.299, 373.63] (0.95 confidence)</td>
<td>372.057 ns [371.029, 373.366] (0.95 confidence)</td>
<td>1.00203</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d5/adapter</td>
<td>1740.6 ns [1736.69, 1744.49] (0.95 confidence)</td>
<td>1731.09 ns [1724.79, 1740.8] (0.95 confidence)</td>
<td>1.00549</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d5/construct</td>
<td>1126.65 ns [1120.76, 1133.14] (0.95 confidence)</td>
<td>1107.06 ns [1103.38, 1113.13] (0.95 confidence)</td>
<td>1.01769</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d5/prepared</td>
<td>503.425 ns [501.81, 509.803] (0.95 confidence)</td>
<td>478.565 ns [478.244, 478.955] (0.95 confidence)</td>
<td>1.05195</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d6/adapter</td>
<td>2172.99 ns [2166.42, 2177.85] (0.95 confidence)</td>
<td>2112.19 ns [2105.51, 2118.53] (0.95 confidence)</td>
<td>1.02879</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d6/construct</td>
<td>1440.22 ns [1433.8, 1445.58] (0.95 confidence)</td>
<td>1427.85 ns [1419.61, 1436.8] (0.95 confidence)</td>
<td>1.00866</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d6/prepared</td>
<td>583.103 ns [575.853, 587.153] (0.95 confidence)</td>
<td>571.207 ns [569.878, 571.987] (0.95 confidence)</td>
<td>1.02083</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d7/adapter</td>
<td>2643.32 ns [2632.34, 2652.15] (0.95 confidence)</td>
<td>2630.81 ns [2624.27, 2640.98] (0.95 confidence)</td>
<td>1.00475</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d7/construct</td>
<td>1842.56 ns [1833.97, 1854.08] (0.95 confidence)</td>
<td>1814.96 ns [1809.01, 1822.14] (0.95 confidence)</td>
<td>1.01521</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d7/prepared</td>
<td>676.106 ns [674.365, 677.72] (0.95 confidence)</td>
<td>683.242 ns [680.148, 687.12] (0.95 confidence)</td>
<td>0.989556</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d8/adapter</td>
<td>3054.18 ns [3047.85, 3060.74] (0.95 confidence)</td>
<td>3040.7 ns [3033.29, 3051.43] (0.95 confidence)</td>
<td>1.00443</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d8/construct</td>
<td>2149.48 ns [2138.25, 2157.44] (0.95 confidence)</td>
<td>2164.67 ns [2158.55, 2172.14] (0.95 confidence)</td>
<td>0.992982</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d8/prepared</td>
<td>759.485 ns [758.537, 760.992] (0.95 confidence)</td>
<td>774.578 ns [772.115, 778.124] (0.95 confidence)</td>
<td>0.980515</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d2/adapter</td>
<td>527.923 ns [525.307, 534.429] (0.95 confidence)</td>
<td>526.471 ns [522.113, 529.131] (0.95 confidence)</td>
<td>1.00276</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d2/construct</td>
<td>240.908 ns [240.417, 241.395] (0.95 confidence)</td>
<td>241.55 ns [240.925, 241.924] (0.95 confidence)</td>
<td>0.997343</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d2/prepared</td>
<td>222.634 ns [222.349, 223.028] (0.95 confidence)</td>
<td>221.053 ns [220.853, 221.452] (0.95 confidence)</td>
<td>1.00715</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d3/adapter</td>
<td>1050.27 ns [1045.69, 1053.61] (0.95 confidence)</td>
<td>1040.02 ns [1038.58, 1043.78] (0.95 confidence)</td>
<td>1.00986</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d3/construct</td>
<td>505.89 ns [501.808, 507.985] (0.95 confidence)</td>
<td>486.442 ns [485.045, 488.116] (0.95 confidence)</td>
<td>1.03998</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d3/prepared</td>
<td>462.832 ns [461.088, 464.28] (0.95 confidence)</td>
<td>452.493 ns [451.642, 453.053] (0.95 confidence)</td>
<td>1.02285</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d4/adapter</td>
<td>1644.4 ns [1640.12, 1647.48] (0.95 confidence)</td>
<td>1616.59 ns [1613.9, 1618.5] (0.95 confidence)</td>
<td>1.0172</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d4/construct</td>
<td>862.535 ns [860.82, 866.544] (0.95 confidence)</td>
<td>832.064 ns [829.731, 835.724] (0.95 confidence)</td>
<td>1.03662</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d4/prepared</td>
<td>687.334 ns [685.55, 689.556] (0.95 confidence)</td>
<td>660.207 ns [659.172, 661.308] (0.95 confidence)</td>
<td>1.04109</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d5/adapter</td>
<td>2449.25 ns [2444.26, 2455.27] (0.95 confidence)</td>
<td>2449.32 ns [2440.1, 2455.23] (0.95 confidence)</td>
<td>0.99997</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d5/construct</td>
<td>1304.67 ns [1299.13, 1308.79] (0.95 confidence)</td>
<td>1251.15 ns [1247.52, 1254.35] (0.95 confidence)</td>
<td>1.04278</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d5/prepared</td>
<td>999.415 ns [998.064, 1001.87] (0.95 confidence)</td>
<td>984.001 ns [981.735, 985.253] (0.95 confidence)</td>
<td>1.01566</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d6/adapter</td>
<td>3404.98 ns [3398.36, 3409.83] (0.95 confidence)</td>
<td>3333.03 ns [3324.59, 3342.65] (0.95 confidence)</td>
<td>1.02159</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d6/construct</td>
<td>1800.52 ns [1776.21, 1808.31] (0.95 confidence)</td>
<td>1780.27 ns [1755.34, 1786.08] (0.95 confidence)</td>
<td>1.01137</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d6/prepared</td>
<td>1430.28 ns [1383.54, 1446.3] (0.95 confidence)</td>
<td>1342.43 ns [1339.55, 1345.08] (0.95 confidence)</td>
<td>1.06544</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d7/adapter</td>
<td>4278.01 ns [4269.59, 4282.1] (0.95 confidence)</td>
<td>4337.41 ns [4328.87, 4343.3] (0.95 confidence)</td>
<td>0.986305</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d7/construct</td>
<td>2286.58 ns [2277.75, 2302.12] (0.95 confidence)</td>
<td>2347.3 ns [2336.55, 2357.73] (0.95 confidence)</td>
<td>0.974131</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d7/prepared</td>
<td>1743.75 ns [1741.5, 1745.55] (0.95 confidence)</td>
<td>1752.49 ns [1748.49, 1757.47] (0.95 confidence)</td>
<td>0.995011</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d8/adapter</td>
<td>5282.71 ns [5271.74, 5296.76] (0.95 confidence)</td>
<td>5322.53 ns [5313.4, 5329.08] (0.95 confidence)</td>
<td>0.992519</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d8/construct</td>
<td>2922.74 ns [2915.03, 2930.28] (0.95 confidence)</td>
<td>3087 ns [3063.28, 3106.09] (0.95 confidence)</td>
<td>0.946791</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d8/prepared</td>
<td>2121.49 ns [2118.39, 2125.6] (0.95 confidence)</td>
<td>2133.98 ns [2130.82, 2138.23] (0.95 confidence)</td>
<td>0.99415</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d2/adapter</td>
<td>649.34 ns [647.178, 652.321] (0.95 confidence)</td>
<td>646.512 ns [645.402, 649.45] (0.95 confidence)</td>
<td>1.00437</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d2/construct</td>
<td>348.42 ns [346.802, 350.037] (0.95 confidence)</td>
<td>350.225 ns [349.393, 351.79] (0.95 confidence)</td>
<td>0.994846</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d2/prepared</td>
<td>215.854 ns [215.608, 216.095] (0.95 confidence)</td>
<td>216.918 ns [216.556, 217.267] (0.95 confidence)</td>
<td>0.995096</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d3/adapter</td>
<td>1168.37 ns [1167.09, 1170.26] (0.95 confidence)</td>
<td>1155.79 ns [1154.35, 1157.59] (0.95 confidence)</td>
<td>1.01088</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d3/construct</td>
<td>666.764 ns [663.847, 669.519] (0.95 confidence)</td>
<td>662.75 ns [660.162, 666.669] (0.95 confidence)</td>
<td>1.00606</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d3/prepared</td>
<td>367.927 ns [367.343, 368.488] (0.95 confidence)</td>
<td>365.313 ns [364.361, 365.96] (0.95 confidence)</td>
<td>1.00716</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d4/adapter</td>
<td>1884.13 ns [1879.45, 1885.98] (0.95 confidence)</td>
<td>1813.04 ns [1811.05, 1814.7] (0.95 confidence)</td>
<td>1.03921</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d4/construct</td>
<td>1103.68 ns [1097.4, 1109.73] (0.95 confidence)</td>
<td>1069.2 ns [1064.17, 1079.35] (0.95 confidence)</td>
<td>1.03225</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d4/prepared</td>
<td>612.639 ns [612.129, 613.431] (0.95 confidence)</td>
<td>589.352 ns [588.932, 589.803] (0.95 confidence)</td>
<td>1.03951</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d5/adapter</td>
<td>2621.41 ns [2617.83, 2623.86] (0.95 confidence)</td>
<td>2594.28 ns [2590.38, 2602.38] (0.95 confidence)</td>
<td>1.01046</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d5/construct</td>
<td>1608.52 ns [1602.77, 1614.02] (0.95 confidence)</td>
<td>1625.57 ns [1615.22, 1633.84] (0.95 confidence)</td>
<td>0.98951</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d5/prepared</td>
<td>820.846 ns [819.581, 822.826] (0.95 confidence)</td>
<td>812.177 ns [811.252, 813.374] (0.95 confidence)</td>
<td>1.01067</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d6/adapter</td>
<td>3659.06 ns [3645.84, 3675.61] (0.95 confidence)</td>
<td>3570.8 ns [3562.43, 3578.75] (0.95 confidence)</td>
<td>1.02472</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d6/construct</td>
<td>2353.98 ns [2349.22, 2357.71] (0.95 confidence)</td>
<td>2325.91 ns [2320.22, 2333.33] (0.95 confidence)</td>
<td>1.01207</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d6/prepared</td>
<td>1054.06 ns [1051.57, 1056.13] (0.95 confidence)</td>
<td>1040.68 ns [1037.72, 1044.06] (0.95 confidence)</td>
<td>1.01286</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d7/adapter</td>
<td>4712.4 ns [4702.65, 4736.69] (0.95 confidence)</td>
<td>4689.86 ns [4687.2, 4701.93] (0.95 confidence)</td>
<td>1.0048</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d7/construct</td>
<td>3050.13 ns [3039.19, 3059.56] (0.95 confidence)</td>
<td>3055.61 ns [3039.24, 3063.16] (0.95 confidence)</td>
<td>0.998205</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d7/prepared</td>
<td>1384.27 ns [1378.29, 1389.17] (0.95 confidence)</td>
<td>1383.02 ns [1379.71, 1387.75] (0.95 confidence)</td>
<td>1.00091</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d8/adapter</td>
<td>5837.7 ns [5830.72, 5849.41] (0.95 confidence)</td>
<td>5843.74 ns [5831.25, 5861.71] (0.95 confidence)</td>
<td>0.998967</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d8/construct</td>
<td>3940.6 ns [3931.65, 3953.21] (0.95 confidence)</td>
<td>3943.38 ns [3926.88, 3952.38] (0.95 confidence)</td>
<td>0.999294</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d8/prepared</td>
<td>1627.66 ns [1622.75, 1633.48] (0.95 confidence)</td>
<td>1627.73 ns [1623.97, 1632.54] (0.95 confidence)</td>
<td>0.999956</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d2/adapter</td>
<td>532.213 ns [526.527, 537.558] (0.95 confidence)</td>
<td>529.658 ns [523.065, 535.313] (0.95 confidence)</td>
<td>1.00482</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d2/construct</td>
<td>243.646 ns [243.1, 244.197] (0.95 confidence)</td>
<td>243.152 ns [242.551, 243.527] (0.95 confidence)</td>
<td>1.00204</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d2/prepared</td>
<td>223.309 ns [222.627, 223.724] (0.95 confidence)</td>
<td>220.666 ns [220.048, 221.208] (0.95 confidence)</td>
<td>1.01198</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d3/adapter</td>
<td>1042.54 ns [1039.76, 1045.97] (0.95 confidence)</td>
<td>1044.7 ns [1042.41, 1047.63] (0.95 confidence)</td>
<td>0.997931</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d3/construct</td>
<td>490.074 ns [488.133, 493.611] (0.95 confidence)</td>
<td>485.413 ns [484.696, 486.764] (0.95 confidence)</td>
<td>1.0096</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d3/prepared</td>
<td>456.729 ns [456.067, 457.746] (0.95 confidence)</td>
<td>450.992 ns [450.633, 451.746] (0.95 confidence)</td>
<td>1.01272</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d4/adapter</td>
<td>1655.09 ns [1651.41, 1661.72] (0.95 confidence)</td>
<td>1623.48 ns [1620.13, 1625.47] (0.95 confidence)</td>
<td>1.01947</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d4/construct</td>
<td>846.04 ns [843.755, 849.129] (0.95 confidence)</td>
<td>836.919 ns [825.901, 851.705] (0.95 confidence)</td>
<td>1.0109</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d4/prepared</td>
<td>677.401 ns [676.624, 678.627] (0.95 confidence)</td>
<td>659.866 ns [659.319, 660.468] (0.95 confidence)</td>
<td>1.02657</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d5/adapter</td>
<td>2528.97 ns [2522.68, 2538.27] (0.95 confidence)</td>
<td>2451.45 ns [2449.03, 2453.31] (0.95 confidence)</td>
<td>1.03162</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d5/construct</td>
<td>1303.91 ns [1296.02, 1313.71] (0.95 confidence)</td>
<td>1246.25 ns [1241.51, 1252.78] (0.95 confidence)</td>
<td>1.04626</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d5/prepared</td>
<td>1029.15 ns [1027.15, 1030.66] (0.95 confidence)</td>
<td>984.715 ns [983.633, 986.023] (0.95 confidence)</td>
<td>1.04513</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d6/adapter</td>
<td>3427.48 ns [3422.6, 3432.45] (0.95 confidence)</td>
<td>3320.55 ns [3312.49, 3323.41] (0.95 confidence)</td>
<td>1.0322</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d6/construct</td>
<td>1805.21 ns [1793.65, 1812.45] (0.95 confidence)</td>
<td>1739.2 ns [1731.48, 1752.87] (0.95 confidence)</td>
<td>1.03795</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d6/prepared</td>
<td>1373.71 ns [1371.32, 1375.67] (0.95 confidence)</td>
<td>1338.26 ns [1335.87, 1343.1] (0.95 confidence)</td>
<td>1.02649</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d7/adapter</td>
<td>4291.13 ns [4279.93, 4306.77] (0.95 confidence)</td>
<td>4356.03 ns [4346.27, 4365.8] (0.95 confidence)</td>
<td>0.985102</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d7/construct</td>
<td>2273.24 ns [2264.71, 2287.73] (0.95 confidence)</td>
<td>2288.65 ns [2279.41, 2306.48] (0.95 confidence)</td>
<td>0.993268</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d7/prepared</td>
<td>1745.82 ns [1742.56, 1748.32] (0.95 confidence)</td>
<td>1757.58 ns [1754.67, 1760.55] (0.95 confidence)</td>
<td>0.99331</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d8/adapter</td>
<td>5275.52 ns [5264.7, 5290.5] (0.95 confidence)</td>
<td>5291.13 ns [5283, 5305.25] (0.95 confidence)</td>
<td>0.997048</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d8/construct</td>
<td>2964.93 ns [2936.3, 2982.58] (0.95 confidence)</td>
<td>2977.44 ns [2961.35, 2993.7] (0.95 confidence)</td>
<td>0.995797</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d8/prepared</td>
<td>2139.95 ns [2132.36, 2147.68] (0.95 confidence)</td>
<td>2163.85 ns [2160.64, 2168.03] (0.95 confidence)</td>
<td>0.988957</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d2/adapter</td>
<td>651.752 ns [648.924, 660.928] (0.95 confidence)</td>
<td>664.138 ns [661.725, 665.444] (0.95 confidence)</td>
<td>0.98135</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d2/construct</td>
<td>349.965 ns [348.184, 351.927] (0.95 confidence)</td>
<td>351.686 ns [350.497, 353.533] (0.95 confidence)</td>
<td>0.995106</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d2/prepared</td>
<td>217.033 ns [216.517, 218.011] (0.95 confidence)</td>
<td>216.962 ns [216.454, 217.402] (0.95 confidence)</td>
<td>1.00033</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d3/adapter</td>
<td>1153.82 ns [1150.53, 1156.47] (0.95 confidence)</td>
<td>1152.3 ns [1150.13, 1153.29] (0.95 confidence)</td>
<td>1.00132</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d3/construct</td>
<td>666.143 ns [664.71, 667.854] (0.95 confidence)</td>
<td>666.414 ns [664.213, 670.253] (0.95 confidence)</td>
<td>0.999593</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d3/prepared</td>
<td>367.213 ns [366.538, 368.418] (0.95 confidence)</td>
<td>364.822 ns [363.951, 365.756] (0.95 confidence)</td>
<td>1.00655</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d4/adapter</td>
<td>1844.65 ns [1842.28, 1848.29] (0.95 confidence)</td>
<td>1820.78 ns [1818.34, 1823.64] (0.95 confidence)</td>
<td>1.01311</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d4/construct</td>
<td>1081.18 ns [1076.51, 1085.74] (0.95 confidence)</td>
<td>1065.77 ns [1061.58, 1071.2] (0.95 confidence)</td>
<td>1.01445</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d4/prepared</td>
<td>603.114 ns [601.942, 604.037] (0.95 confidence)</td>
<td>593.675 ns [592.457, 594.358] (0.95 confidence)</td>
<td>1.0159</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d5/adapter</td>
<td>2714.62 ns [2709.82, 2717.02] (0.95 confidence)</td>
<td>2605.4 ns [2599.25, 2608.94] (0.95 confidence)</td>
<td>1.04192</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d5/construct</td>
<td>1649.31 ns [1640.86, 1653.56] (0.95 confidence)</td>
<td>1606.99 ns [1596.58, 1614.3] (0.95 confidence)</td>
<td>1.02634</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d5/prepared</td>
<td>837.039 ns [834.056, 839.541] (0.95 confidence)</td>
<td>814.028 ns [812.411, 815.618] (0.95 confidence)</td>
<td>1.02827</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d6/adapter</td>
<td>3658.09 ns [3653.45, 3666.2] (0.95 confidence)</td>
<td>3550.2 ns [3543.74, 3561] (0.95 confidence)</td>
<td>1.03039</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d6/construct</td>
<td>2384.54 ns [2374.37, 2394.37] (0.95 confidence)</td>
<td>2339.99 ns [2331.04, 2345.58] (0.95 confidence)</td>
<td>1.01903</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d6/prepared</td>
<td>1063.05 ns [1059.17, 1069.82] (0.95 confidence)</td>
<td>1046.53 ns [1040.4, 1052.85] (0.95 confidence)</td>
<td>1.01578</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d7/adapter</td>
<td>4687.04 ns [4682.13, 4692.95] (0.95 confidence)</td>
<td>4705.43 ns [4692.44, 4727.42] (0.95 confidence)</td>
<td>0.996092</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d7/construct</td>
<td>3016.65 ns [3009.21, 3030.01] (0.95 confidence)</td>
<td>3079.93 ns [3063.8, 3091.07] (0.95 confidence)</td>
<td>0.979455</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d7/prepared</td>
<td>1368.47 ns [1366.44, 1374.9] (0.95 confidence)</td>
<td>1378.02 ns [1375.02, 1381.87] (0.95 confidence)</td>
<td>0.993073</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d8/adapter</td>
<td>5856.42 ns [5848.91, 5867.11] (0.95 confidence)</td>
<td>5849.61 ns [5843.93, 5859.72] (0.95 confidence)</td>
<td>1.00116</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d8/construct</td>
<td>3894.06 ns [3878.9, 3912.15] (0.95 confidence)</td>
<td>3963.35 ns [3951.27, 3978.52] (0.95 confidence)</td>
<td>0.982516</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d8/prepared</td>
<td>1630.1 ns [1627.42, 1632.55] (0.95 confidence)</td>
<td>1646.67 ns [1642.06, 1654.56] (0.95 confidence)</td>
<td>0.98994</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d2/adapter</td>
<td>532.654 ns [526.088, 535.252] (0.95 confidence)</td>
<td>532.774 ns [527.596, 539.072] (0.95 confidence)</td>
<td>0.999774</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d2/construct</td>
<td>237.146 ns [236.711, 237.816] (0.95 confidence)</td>
<td>237.726 ns [237.078, 238.222] (0.95 confidence)</td>
<td>0.997561</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d2/prepared</td>
<td>222.355 ns [221.908, 222.723] (0.95 confidence)</td>
<td>222.186 ns [221.652, 222.42] (0.95 confidence)</td>
<td>1.00076</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d3/adapter</td>
<td>1013.42 ns [1011.38, 1014.71] (0.95 confidence)</td>
<td>1019.8 ns [1017.34, 1021.6] (0.95 confidence)</td>
<td>0.993739</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d3/construct</td>
<td>479.036 ns [477.653, 480.868] (0.95 confidence)</td>
<td>471.089 ns [470.045, 472.094] (0.95 confidence)</td>
<td>1.01687</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d3/prepared</td>
<td>440.647 ns [439.389, 442.318] (0.95 confidence)</td>
<td>438.658 ns [438.036, 439.405] (0.95 confidence)</td>
<td>1.00453</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d4/adapter</td>
<td>1676.59 ns [1671.81, 1683.55] (0.95 confidence)</td>
<td>1653.93 ns [1651.6, 1657.47] (0.95 confidence)</td>
<td>1.0137</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d4/construct</td>
<td>846.4 ns [842.654, 850.39] (0.95 confidence)</td>
<td>829.325 ns [821.496, 836.087] (0.95 confidence)</td>
<td>1.02059</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d4/prepared</td>
<td>712.395 ns [709.283, 718.34] (0.95 confidence)</td>
<td>691.943 ns [690.963, 693.134] (0.95 confidence)</td>
<td>1.02956</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d5/adapter</td>
<td>2416.21 ns [2407.05, 2426.94] (0.95 confidence)</td>
<td>2378.96 ns [2375.3, 2383] (0.95 confidence)</td>
<td>1.01566</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d5/construct</td>
<td>1237.52 ns [1226.75, 1246.57] (0.95 confidence)</td>
<td>1203.49 ns [1200.42, 1208.94] (0.95 confidence)</td>
<td>1.02828</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d5/prepared</td>
<td>962.027 ns [960.715, 964.503] (0.95 confidence)</td>
<td>953.397 ns [952.463, 954.05] (0.95 confidence)</td>
<td>1.00905</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d6/adapter</td>
<td>3595.42 ns [3565.37, 3620.5] (0.95 confidence)</td>
<td>3445.79 ns [3435.66, 3454.72] (0.95 confidence)</td>
<td>1.04342</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d6/construct</td>
<td>1895.73 ns [1891.17, 1905.09] (0.95 confidence)</td>
<td>1791.13 ns [1784.38, 1797.26] (0.95 confidence)</td>
<td>1.0584</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d6/prepared</td>
<td>1411.1 ns [1407.31, 1413.3] (0.95 confidence)</td>
<td>1399.9 ns [1394.89, 1406.82] (0.95 confidence)</td>
<td>1.008</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d7/adapter</td>
<td>4373.35 ns [4361.74, 4382.94] (0.95 confidence)</td>
<td>4417.69 ns [4407.28, 4421.78] (0.95 confidence)</td>
<td>0.989963</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d7/construct</td>
<td>2350.78 ns [2340.4, 2357.83] (0.95 confidence)</td>
<td>2342.21 ns [2338.08, 2348.21] (0.95 confidence)</td>
<td>1.00366</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d7/prepared</td>
<td>1773.29 ns [1770.34, 1777.53] (0.95 confidence)</td>
<td>1795.02 ns [1788.05, 1798.4] (0.95 confidence)</td>
<td>0.987897</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d8/adapter</td>
<td>5446.53 ns [5440.28, 5464.37] (0.95 confidence)</td>
<td>5894.36 ns [5636.77, 5930.6] (0.95 confidence)</td>
<td>0.924023</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d8/construct</td>
<td>3070.84 ns [3061.33, 3080.46] (0.95 confidence)</td>
<td>3176.6 ns [3161.23, 3197.45] (0.95 confidence)</td>
<td>0.966706</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d8/prepared</td>
<td>2107.62 ns [2104.64, 2108.98] (0.95 confidence)</td>
<td>2131.34 ns [2128.2, 2135.58] (0.95 confidence)</td>
<td>0.988868</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d2/adapter</td>
<td>686.109 ns [683.126, 687.936] (0.95 confidence)</td>
<td>686.879 ns [684.754, 689.136] (0.95 confidence)</td>
<td>0.998878</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d2/construct</td>
<td>374.61 ns [372.686, 375.356] (0.95 confidence)</td>
<td>376.941 ns [375.926, 378.491] (0.95 confidence)</td>
<td>0.993817</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d2/prepared</td>
<td>230.246 ns [229.965, 230.54] (0.95 confidence)</td>
<td>232.558 ns [232.278, 233.489] (0.95 confidence)</td>
<td>0.990061</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d3/adapter</td>
<td>1184.82 ns [1183.29, 1187.64] (0.95 confidence)</td>
<td>1181.82 ns [1178.92, 1183.98] (0.95 confidence)</td>
<td>1.00253</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d3/construct</td>
<td>695.846 ns [694.593, 698.324] (0.95 confidence)</td>
<td>691.581 ns [689.946, 692.675] (0.95 confidence)</td>
<td>1.00617</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d3/prepared</td>
<td>385.008 ns [384.356, 386.333] (0.95 confidence)</td>
<td>381.618 ns [380.868, 382.212] (0.95 confidence)</td>
<td>1.00888</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d4/adapter</td>
<td>1814.26 ns [1810.09, 1818.1] (0.95 confidence)</td>
<td>1759.3 ns [1756.39, 1761.24] (0.95 confidence)</td>
<td>1.03124</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d4/construct</td>
<td>1044.67 ns [1035.78, 1057.08] (0.95 confidence)</td>
<td>1002.63 ns [996.52, 1013.71] (0.95 confidence)</td>
<td>1.04193</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d4/prepared</td>
<td>636.228 ns [635.278, 637.697] (0.95 confidence)</td>
<td>624.01 ns [622.702, 624.641] (0.95 confidence)</td>
<td>1.01958</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d5/adapter</td>
<td>2509.85 ns [2503.24, 2513.07] (0.95 confidence)</td>
<td>2495.93 ns [2493.75, 2499.72] (0.95 confidence)</td>
<td>1.00558</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d5/construct</td>
<td>1459.87 ns [1449.46, 1469.77] (0.95 confidence)</td>
<td>1449.6 ns [1446.33, 1454.04] (0.95 confidence)</td>
<td>1.00708</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d5/prepared</td>
<td>881.76 ns [880.239, 882.938] (0.95 confidence)</td>
<td>875.341 ns [874.351, 876.431] (0.95 confidence)</td>
<td>1.00733</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d6/adapter</td>
<td>3528.55 ns [3513.05, 3537.18] (0.95 confidence)</td>
<td>3392.97 ns [3383.49, 3404.88] (0.95 confidence)</td>
<td>1.03996</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d6/construct</td>
<td>2100.5 ns [2082.99, 2108.51] (0.95 confidence)</td>
<td>2046.86 ns [2038.52, 2054.77] (0.95 confidence)</td>
<td>1.02621</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d6/prepared</td>
<td>1168.58 ns [1166.5, 1173.33] (0.95 confidence)</td>
<td>1135.94 ns [1134.37, 1138.34] (0.95 confidence)</td>
<td>1.02874</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d7/adapter</td>
<td>4437.14 ns [4426.93, 4444.63] (0.95 confidence)</td>
<td>4503.76 ns [4478.04, 4533.07] (0.95 confidence)</td>
<td>0.985207</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d7/construct</td>
<td>2748.8 ns [2742.81, 2770.7] (0.95 confidence)</td>
<td>2823.28 ns [2817.92, 2833.74] (0.95 confidence)</td>
<td>0.973619</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d7/prepared</td>
<td>1468.13 ns [1466.33, 1471.17] (0.95 confidence)</td>
<td>1477.75 ns [1472.68, 1480.06] (0.95 confidence)</td>
<td>0.993488</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d8/adapter</td>
<td>5763.35 ns [5746.25, 5794.59] (0.95 confidence)</td>
<td>5746.41 ns [5730.07, 5758.18] (0.95 confidence)</td>
<td>1.00295</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d8/construct</td>
<td>3571.95 ns [3561.71, 3579.09] (0.95 confidence)</td>
<td>3607.19 ns [3587.79, 3623.85] (0.95 confidence)</td>
<td>0.990229</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d8/prepared</td>
<td>1869.29 ns [1867.47, 1871.54] (0.95 confidence)</td>
<td>1883.44 ns [1879.71, 1889.34] (0.95 confidence)</td>
<td>0.992487</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d2/adapter</td>
<td>1036.4 ns [1009.88, 1042.9] (0.95 confidence)</td>
<td>1035.32 ns [1015.22, 1050.67] (0.95 confidence)</td>
<td>1.00105</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d2/construct</td>
<td>296.002 ns [295.018, 296.526] (0.95 confidence)</td>
<td>294.806 ns [294.098, 295.587] (0.95 confidence)</td>
<td>1.00406</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d2/prepared</td>
<td>635.547 ns [634.314, 636.639] (0.95 confidence)</td>
<td>630.187 ns [629.261, 631.98] (0.95 confidence)</td>
<td>1.00851</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d3/adapter</td>
<td>1805.12 ns [1801.23, 1807.66] (0.95 confidence)</td>
<td>1782.88 ns [1777.36, 1786.37] (0.95 confidence)</td>
<td>1.01247</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d3/construct</td>
<td>585.553 ns [584.462, 586.851] (0.95 confidence)</td>
<td>582.329 ns [581.222, 583.304] (0.95 confidence)</td>
<td>1.00554</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d3/prepared</td>
<td>1090.75 ns [1087.17, 1092.97] (0.95 confidence)</td>
<td>1084.11 ns [1082.73, 1087.49] (0.95 confidence)</td>
<td>1.00612</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d4/adapter</td>
<td>2866.1 ns [2854.86, 2877.85] (0.95 confidence)</td>
<td>2742.43 ns [2737.23, 2750.21] (0.95 confidence)</td>
<td>1.04509</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d4/construct</td>
<td>847.484 ns [829.778, 857.942] (0.95 confidence)</td>
<td>810.156 ns [803.063, 820.224] (0.95 confidence)</td>
<td>1.04607</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d4/prepared</td>
<td>1987.41 ns [1977.87, 2003.06] (0.95 confidence)</td>
<td>1780.87 ns [1778.53, 1785.49] (0.95 confidence)</td>
<td>1.11598</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d5/adapter</td>
<td>4292.56 ns [4286.15, 4299.67] (0.95 confidence)</td>
<td>4211.31 ns [4203.25, 4214.79] (0.95 confidence)</td>
<td>1.01929</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d5/construct</td>
<td>1378.46 ns [1372.31, 1386.18] (0.95 confidence)</td>
<td>1339.9 ns [1336.11, 1346.44] (0.95 confidence)</td>
<td>1.02878</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d5/prepared</td>
<td>2706.84 ns [2701.26, 2710.55] (0.95 confidence)</td>
<td>2629.96 ns [2624.23, 2633.74] (0.95 confidence)</td>
<td>1.02924</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d6/adapter</td>
<td>5803.77 ns [5788.3, 5818.03] (0.95 confidence)</td>
<td>5682.75 ns [5662.9, 5706.29] (0.95 confidence)</td>
<td>1.0213</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d6/construct</td>
<td>1764.12 ns [1756.06, 1773.9] (0.95 confidence)</td>
<td>1696.14 ns [1689.11, 1705.42] (0.95 confidence)</td>
<td>1.04008</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d6/prepared</td>
<td>3757.04 ns [3749.48, 3764.01] (0.95 confidence)</td>
<td>3679.18 ns [3673.12, 3686.45] (0.95 confidence)</td>
<td>1.02116</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d7/adapter</td>
<td>7831.53 ns [7802.7, 7844.05] (0.95 confidence)</td>
<td>7859.89 ns [7846.15, 7870.64] (0.95 confidence)</td>
<td>0.996391</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d7/construct</td>
<td>2187.75 ns [2176.02, 2196.38] (0.95 confidence)</td>
<td>2232.42 ns [2223.69, 2249.58] (0.95 confidence)</td>
<td>0.979989</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d7/prepared</td>
<td>5265.25 ns [5259.07, 5270.25] (0.95 confidence)</td>
<td>5300.71 ns [5293.38, 5307.75] (0.95 confidence)</td>
<td>0.993312</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d8/adapter</td>
<td>9505.9 ns [9490.99, 9522.4] (0.95 confidence)</td>
<td>9558.04 ns [9539.01, 9575.93] (0.95 confidence)</td>
<td>0.994544</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d8/construct</td>
<td>2411.51 ns [2401.57, 2426.35] (0.95 confidence)</td>
<td>2441.17 ns [2424.57, 2449.56] (0.95 confidence)</td>
<td>0.987851</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d8/prepared</td>
<td>6988.31 ns [6972.83, 7003.92] (0.95 confidence)</td>
<td>7013.41 ns [7000.06, 7071.39] (0.95 confidence)</td>
<td>0.996421</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d2/adapter</td>
<td>1175.35 ns [1165.97, 1184.73] (0.95 confidence)</td>
<td>1174.91 ns [1172.14, 1178.86] (0.95 confidence)</td>
<td>1.00037</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d2/construct</td>
<td>478.605 ns [477.396, 479.378] (0.95 confidence)</td>
<td>481.065 ns [478.624, 482.092] (0.95 confidence)</td>
<td>0.994885</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d2/prepared</td>
<td>585.259 ns [584.853, 586.33] (0.95 confidence)</td>
<td>587.241 ns [585.847, 588.101] (0.95 confidence)</td>
<td>0.996624</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d3/adapter</td>
<td>1916.39 ns [1912.54, 1920.68] (0.95 confidence)</td>
<td>1918.11 ns [1916.01, 1921.11] (0.95 confidence)</td>
<td>0.999101</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d3/construct</td>
<td>752.145 ns [749.482, 754.802] (0.95 confidence)</td>
<td>750.529 ns [747.517, 753.197] (0.95 confidence)</td>
<td>1.00215</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d3/prepared</td>
<td>1031.3 ns [1029.61, 1032.5] (0.95 confidence)</td>
<td>1033.74 ns [1032.75, 1034.99] (0.95 confidence)</td>
<td>0.997644</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d4/adapter</td>
<td>3248.67 ns [3238.53, 3261.07] (0.95 confidence)</td>
<td>3156.87 ns [3154.99, 3160.74] (0.95 confidence)</td>
<td>1.02908</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d4/construct</td>
<td>1311.56 ns [1304.08, 1323.84] (0.95 confidence)</td>
<td>1279.37 ns [1275.37, 1285.49] (0.95 confidence)</td>
<td>1.02516</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d4/prepared</td>
<td>1718.78 ns [1707.99, 1724.41] (0.95 confidence)</td>
<td>1684.46 ns [1679.99, 1686.57] (0.95 confidence)</td>
<td>1.02037</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d5/adapter</td>
<td>4448.34 ns [4441.61, 4458.83] (0.95 confidence)</td>
<td>4416.32 ns [4411.69, 4423.89] (0.95 confidence)</td>
<td>1.00725</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d5/construct</td>
<td>1737.1 ns [1730.29, 1744.62] (0.95 confidence)</td>
<td>1677.39 ns [1673.26, 1682.66] (0.95 confidence)</td>
<td>1.0356</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d5/prepared</td>
<td>2582.18 ns [2578.95, 2586.92] (0.95 confidence)</td>
<td>2491.82 ns [2488.03, 2494.97] (0.95 confidence)</td>
<td>1.03626</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d6/adapter</td>
<td>6164.82 ns [6152.51, 6181.22] (0.95 confidence)</td>
<td>6099.39 ns [6066.52, 6117.93] (0.95 confidence)</td>
<td>1.01073</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d6/construct</td>
<td>2376.09 ns [2362.78, 2382.36] (0.95 confidence)</td>
<td>2266.71 ns [2262.51, 2281.69] (0.95 confidence)</td>
<td>1.04825</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d6/prepared</td>
<td>3565.2 ns [3551.79, 3575.52] (0.95 confidence)</td>
<td>3490.93 ns [3482.44, 3499.35] (0.95 confidence)</td>
<td>1.02127</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d7/adapter</td>
<td>8540.39 ns [8530.94, 8558.54] (0.95 confidence)</td>
<td>8655.62 ns [8638.85, 8669.09] (0.95 confidence)</td>
<td>0.986688</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d7/construct</td>
<td>2890.47 ns [2877.76, 2900.4] (0.95 confidence)</td>
<td>2910.22 ns [2897.92, 2930.53] (0.95 confidence)</td>
<td>0.993213</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d7/prepared</td>
<td>5422.82 ns [5416.8, 5428.49] (0.95 confidence)</td>
<td>5452.35 ns [5445.15, 5455.53] (0.95 confidence)</td>
<td>0.994584</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d8/adapter</td>
<td>10486.6 ns [10476.6, 10492.2] (0.95 confidence)</td>
<td>10578.9 ns [10551.9, 10602.9] (0.95 confidence)</td>
<td>0.991281</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d8/construct</td>
<td>3422.99 ns [3406.92, 3434.86] (0.95 confidence)</td>
<td>3455.24 ns [3444.02, 3469.94] (0.95 confidence)</td>
<td>0.990667</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d8/prepared</td>
<td>6856.76 ns [6850.5, 6879.61] (0.95 confidence)</td>
<td>6926.23 ns [6910.15, 6949.65] (0.95 confidence)</td>
<td>0.98997</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d2/adapter</td>
<td>530.963 ns [526.643, 535.512] (0.95 confidence)</td>
<td>531.846 ns [525.934, 536.529] (0.95 confidence)</td>
<td>0.998339</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d2/construct</td>
<td>238.408 ns [237.296, 239.188] (0.95 confidence)</td>
<td>239.798 ns [239.457, 240.475] (0.95 confidence)</td>
<td>0.994205</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d2/prepared</td>
<td>222.908 ns [222.561, 223.494] (0.95 confidence)</td>
<td>221.539 ns [221.126, 221.852] (0.95 confidence)</td>
<td>1.00618</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d3/adapter</td>
<td>1013.13 ns [1011.83, 1014.84] (0.95 confidence)</td>
<td>1017.97 ns [1016.76, 1019.48] (0.95 confidence)</td>
<td>0.995245</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d3/construct</td>
<td>480.387 ns [478.338, 483.059] (0.95 confidence)</td>
<td>476.627 ns [475.935, 477.751] (0.95 confidence)</td>
<td>1.00789</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d3/prepared</td>
<td>439.924 ns [439.362, 440.497] (0.95 confidence)</td>
<td>441.522 ns [440.982, 442.18] (0.95 confidence)</td>
<td>0.996381</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d4/adapter</td>
<td>1708.37 ns [1706.06, 1712.56] (0.95 confidence)</td>
<td>1644.77 ns [1642.42, 1647.44] (0.95 confidence)</td>
<td>1.03867</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d4/construct</td>
<td>852.114 ns [845.458, 860.234] (0.95 confidence)</td>
<td>829.392 ns [824.935, 834.065] (0.95 confidence)</td>
<td>1.0274</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d4/prepared</td>
<td>721.728 ns [720.599, 722.76] (0.95 confidence)</td>
<td>689.895 ns [689.157, 690.7] (0.95 confidence)</td>
<td>1.04614</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d5/adapter</td>
<td>2410.8 ns [2403.87, 2418.49] (0.95 confidence)</td>
<td>2361.62 ns [2358.49, 2373.56] (0.95 confidence)</td>
<td>1.02082</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d5/construct</td>
<td>1293.19 ns [1287.85, 1301.4] (0.95 confidence)</td>
<td>1227.71 ns [1212.72, 1239.74] (0.95 confidence)</td>
<td>1.05334</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d5/prepared</td>
<td>991.437 ns [988.529, 992.944] (0.95 confidence)</td>
<td>955.233 ns [952.585, 959.697] (0.95 confidence)</td>
<td>1.0379</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d6/adapter</td>
<td>3461.06 ns [3454.34, 3471.78] (0.95 confidence)</td>
<td>3435.56 ns [3426.63, 3441.16] (0.95 confidence)</td>
<td>1.00742</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d6/construct</td>
<td>1888.65 ns [1879, 1896.24] (0.95 confidence)</td>
<td>1835.96 ns [1820.46, 1844.29] (0.95 confidence)</td>
<td>1.0287</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d6/prepared</td>
<td>1429.02 ns [1426.47, 1434.21] (0.95 confidence)</td>
<td>1370.29 ns [1368.57, 1372.07] (0.95 confidence)</td>
<td>1.04286</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d7/adapter</td>
<td>4381.07 ns [4370.02, 4390.43] (0.95 confidence)</td>
<td>4449.05 ns [4433.22, 4462.42] (0.95 confidence)</td>
<td>0.984721</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d7/construct</td>
<td>2397.64 ns [2364.74, 2417.37] (0.95 confidence)</td>
<td>2423.36 ns [2403.43, 2437.96] (0.95 confidence)</td>
<td>0.989385</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d7/prepared</td>
<td>1775.19 ns [1773.27, 1777.6] (0.95 confidence)</td>
<td>1797.44 ns [1790.95, 1804.58] (0.95 confidence)</td>
<td>0.987623</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d8/adapter</td>
<td>5443.64 ns [5435.62, 5452.76] (0.95 confidence)</td>
<td>5463.3 ns [5457.38, 5473.12] (0.95 confidence)</td>
<td>0.996401</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d8/construct</td>
<td>3063.4 ns [3055.6, 3075.28] (0.95 confidence)</td>
<td>3090.55 ns [3080.5, 3094.55] (0.95 confidence)</td>
<td>0.991217</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d8/prepared</td>
<td>2109.11 ns [2104.31, 2111.48] (0.95 confidence)</td>
<td>2131.29 ns [2127.41, 2133.43] (0.95 confidence)</td>
<td>0.989591</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d2/adapter</td>
<td>684.938 ns [682.777, 688.075] (0.95 confidence)</td>
<td>688.701 ns [686.85, 691.926] (0.95 confidence)</td>
<td>0.994537</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d2/construct</td>
<td>377.735 ns [376.279, 379.42] (0.95 confidence)</td>
<td>383.696 ns [382.18, 384.87] (0.95 confidence)</td>
<td>0.984463</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d2/prepared</td>
<td>229.656 ns [229.138, 230.162] (0.95 confidence)</td>
<td>230.971 ns [230.638, 231.836] (0.95 confidence)</td>
<td>0.994307</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d3/adapter</td>
<td>1186.6 ns [1184.56, 1189.77] (0.95 confidence)</td>
<td>1181.66 ns [1178.4, 1184.52] (0.95 confidence)</td>
<td>1.00418</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d3/construct</td>
<td>704.378 ns [702.582, 706.413] (0.95 confidence)</td>
<td>697.535 ns [694.119, 700.216] (0.95 confidence)</td>
<td>1.00981</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d3/prepared</td>
<td>385.094 ns [384.23, 385.866] (0.95 confidence)</td>
<td>384.303 ns [383.691, 385.023] (0.95 confidence)</td>
<td>1.00206</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d4/adapter</td>
<td>1828.83 ns [1823.69, 1841.18] (0.95 confidence)</td>
<td>1762.55 ns [1758.61, 1767.06] (0.95 confidence)</td>
<td>1.0376</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d4/construct</td>
<td>1051.19 ns [1045.14, 1060.24] (0.95 confidence)</td>
<td>1018.22 ns [1007.65, 1029.66] (0.95 confidence)</td>
<td>1.03239</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d4/prepared</td>
<td>635.291 ns [634.366, 637.909] (0.95 confidence)</td>
<td>621.535 ns [620.358, 622.705] (0.95 confidence)</td>
<td>1.02213</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d5/adapter</td>
<td>2651.05 ns [2629.32, 2667.26] (0.95 confidence)</td>
<td>2525.56 ns [2516.44, 2530.28] (0.95 confidence)</td>
<td>1.04969</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d5/construct</td>
<td>1507.68 ns [1496.24, 1515.78] (0.95 confidence)</td>
<td>1461.32 ns [1458.38, 1466.91] (0.95 confidence)</td>
<td>1.03173</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d5/prepared</td>
<td>925.091 ns [923.14, 926.85] (0.95 confidence)</td>
<td>875.854 ns [874.514, 877.755] (0.95 confidence)</td>
<td>1.05622</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d6/adapter</td>
<td>3541.31 ns [3535.45, 3559.51] (0.95 confidence)</td>
<td>3411.18 ns [3395.78, 3428.79] (0.95 confidence)</td>
<td>1.03815</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d6/construct</td>
<td>2165.71 ns [2149.07, 2173.83] (0.95 confidence)</td>
<td>2015.16 ns [2008.56, 2023.96] (0.95 confidence)</td>
<td>1.07471</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d6/prepared</td>
<td>1181.7 ns [1179.58, 1184.16] (0.95 confidence)</td>
<td>1131.78 ns [1130.57, 1133.91] (0.95 confidence)</td>
<td>1.04411</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d7/adapter</td>
<td>4480.86 ns [4464.88, 4494.54] (0.95 confidence)</td>
<td>4488.14 ns [4467.88, 4505.16] (0.95 confidence)</td>
<td>0.998378</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d7/construct</td>
<td>2778.72 ns [2768.44, 2789.46] (0.95 confidence)</td>
<td>2842.85 ns [2829.3, 2862.13] (0.95 confidence)</td>
<td>0.977443</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d7/prepared</td>
<td>1472.92 ns [1469.95, 1478.46] (0.95 confidence)</td>
<td>1468.98 ns [1467.14, 1472.3] (0.95 confidence)</td>
<td>1.00268</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d8/adapter</td>
<td>5736.36 ns [5718.14, 5746.02] (0.95 confidence)</td>
<td>5736.43 ns [5716.1, 5768.26] (0.95 confidence)</td>
<td>0.999986</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d8/construct</td>
<td>3564.16 ns [3556.09, 3573.04] (0.95 confidence)</td>
<td>3602.86 ns [3590.57, 3611.37] (0.95 confidence)</td>
<td>0.989259</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d8/prepared</td>
<td>1865.2 ns [1861.57, 1871.87] (0.95 confidence)</td>
<td>1878.32 ns [1873.89, 1884.77] (0.95 confidence)</td>
<td>0.993015</td>
</tr>
</tbody>
</table>

## lu-solve

Semantic contract: la-stack-lu-solve-v1. Collection session: `b90b2cc4e14749b9914c75a3a05586f6`.
Collection interval (Unix ns): 1794387070544411000–1794387718718775000.
Acquisition: newly collected.

- baseline: available — inside declared release interval; release v0.5.2; revision `eb5e44e52399401311e676faad590e6cce1610e5`.
  Source SHA-256: `bdd47c0e9bc3126dc3a32a98e0f1a688aaa52e7c6a3868009f1adfb0735f367e`; harness SHA-256:
  `0fedb76c77086b1cf8a698034ac4a1f3539fe06feda597672e2d11f01704f9cb`; scientific gate: passed.
- current: available — inside declared release interval; release v0.5.3; revision `be9277bf04400dccbe1115da69c7fda2141c4a39`.
  Source SHA-256: `793169b4c5085624d7c0e773f1d28504c017e9b78332f250150d19df87391a32`; harness SHA-256:
  `0fedb76c77086b1cf8a698034ac4a1f3539fe06feda597672e2d11f01704f9cb`; scientific gate: passed.

### Mean (ns)

| Workload | Baseline | Current | Baseline / current |
| --- | --- | --- | --- |
| d16/faer&#95;lu&#95;solve | 891.634 ns [888.601, 894.728] (0.95 confidence) | 889.57 ns [886.341, 892.813] (0.95 confidence) | 1.00232 |
| d16/la&#95;stack&#95;lu&#95;solve | 643.276 ns [637.738, 649.088] (0.95 confidence) | 637.181 ns [632.308, 642.165] (0.95 confidence) | 1.00957 |
| d16/nalgebra&#95;lu&#95;solve | 581.77 ns [580.915, 582.657] (0.95 confidence) | 582.674 ns [581.544, 584.02] (0.95 confidence) | 0.998448 |
| d2/faer&#95;lu&#95;solve | 128.602 ns [128.203, 129.019] (0.95 confidence) | 131.13 ns [130.789, 131.478] (0.95 confidence) | 0.980726 |
| d2/la&#95;stack&#95;lu&#95;solve | 2.06411 ns [2.05621, 2.07343] (0.95 confidence) | 2.10215 ns [2.09053, 2.11465] (0.95 confidence) | 0.981904 |
| d2/nalgebra&#95;lu&#95;solve | 4.55386 ns [4.54805, 4.56037] (0.95 confidence) | 4.59959 ns [4.59367, 4.60624] (0.95 confidence) | 0.990056 |
| d3/faer&#95;lu&#95;solve | 158.251 ns [157.824, 158.755] (0.95 confidence) | 164.427 ns [163.774, 165.139] (0.95 confidence) | 0.962436 |
| d3/la&#95;stack&#95;lu&#95;solve | 9.9268 ns [9.90829, 9.94509] (0.95 confidence) | 10.0648 ns [10.0397, 10.0895] (0.95 confidence) | 0.986289 |
| d3/nalgebra&#95;lu&#95;solve | 34.1425 ns [33.9956, 34.2966] (0.95 confidence) | 35.4423 ns [35.2888, 35.6045] (0.95 confidence) | 0.963324 |
| d32/faer&#95;lu&#95;solve | 2953.64 ns [2948.27, 2959.08] (0.95 confidence) | 2978.41 ns [2973.22, 2983.68] (0.95 confidence) | 0.991682 |
| d32/la&#95;stack&#95;lu&#95;solve | 2768.72 ns [2765.38, 2772.21] (0.95 confidence) | 2788.72 ns [2785.16, 2792.38] (0.95 confidence) | 0.992828 |
| d32/nalgebra&#95;lu&#95;solve | 2738.79 ns [2733.74, 2743.27] (0.95 confidence) | 2700.6 ns [2684.77, 2715.38] (0.95 confidence) | 1.01414 |
| d4/faer&#95;lu&#95;solve | 183.588 ns [183.193, 184.001] (0.95 confidence) | 189.543 ns [189.21, 189.893] (0.95 confidence) | 0.968579 |
| d4/la&#95;stack&#95;lu&#95;solve | 23.6434 ns [23.6152, 23.6734] (0.95 confidence) | 24.0238 ns [24.0089, 24.039] (0.95 confidence) | 0.984168 |
| d4/nalgebra&#95;lu&#95;solve | 53.6516 ns [53.5576, 53.7465] (0.95 confidence) | 55.4183 ns [55.3366, 55.5018] (0.95 confidence) | 0.96812 |
| d5/faer&#95;lu&#95;solve | 250.473 ns [250.003, 250.96] (0.95 confidence) | 252.711 ns [252.19, 253.222] (0.95 confidence) | 0.991146 |
| d5/la&#95;stack&#95;lu&#95;solve | 47.5538 ns [47.3448, 47.7868] (0.95 confidence) | 48.1208 ns [47.9547, 48.2945] (0.95 confidence) | 0.988218 |
| d5/nalgebra&#95;lu&#95;solve | 75.7539 ns [75.5586, 75.9539] (0.95 confidence) | 78.2899 ns [78.1445, 78.4381] (0.95 confidence) | 0.967607 |
| d64/faer&#95;lu&#95;solve | 12912.1 ns [12873.3, 12953.3] (0.95 confidence) | 12975.5 ns [12918.9, 13041] (0.95 confidence) | 0.995118 |
| d64/la&#95;stack&#95;lu&#95;solve | 18627.3 ns [18586.5, 18671.3] (0.95 confidence) | 18092.8 ns [18061.7, 18126.1] (0.95 confidence) | 1.02954 |
| d64/nalgebra&#95;lu&#95;solve | 14456.6 ns [14431.7, 14482.3] (0.95 confidence) | 14374 ns [14356.3, 14393.6] (0.95 confidence) | 1.00575 |
| d8/faer&#95;lu&#95;solve | 335.339 ns [334.466, 336.272] (0.95 confidence) | 338.001 ns [337.317, 338.689] (0.95 confidence) | 0.992125 |
| d8/la&#95;stack&#95;lu&#95;solve | 131.397 ns [131.167, 131.632] (0.95 confidence) | 133.622 ns [133.354, 133.911] (0.95 confidence) | 0.98335 |
| d8/nalgebra&#95;lu&#95;solve | 167.971 ns [167.614, 168.328] (0.95 confidence) | 169.128 ns [168.673, 169.633] (0.95 confidence) | 0.993159 |

### Median (ns)

| Workload | Baseline | Current | Baseline / current |
| --- | --- | --- | --- |
| d16/faer&#95;lu&#95;solve | 891.125 ns [888.358, 895.089] (0.95 confidence) | 889.228 ns [885.972, 894.046] (0.95 confidence) | 1.00213 |
| d16/la&#95;stack&#95;lu&#95;solve | 638.756 ns [635.257, 642.895] (0.95 confidence) | 631.082 ns [627.847, 636.268] (0.95 confidence) | 1.01216 |
| d16/nalgebra&#95;lu&#95;solve | 580.644 ns [580.113, 581.847] (0.95 confidence) | 581.252 ns [580.272, 582.309] (0.95 confidence) | 0.998953 |
| d2/faer&#95;lu&#95;solve | 128.098 ns [127.923, 128.505] (0.95 confidence) | 130.897 ns [130.436, 131.469] (0.95 confidence) | 0.978617 |
| d2/la&#95;stack&#95;lu&#95;solve | 2.04872 ns [2.0471, 2.05208] (0.95 confidence) | 2.06884 ns [2.0664, 2.07466] (0.95 confidence) | 0.990277 |
| d2/nalgebra&#95;lu&#95;solve | 4.54511 ns [4.54048, 4.54985] (0.95 confidence) | 4.5918 ns [4.58866, 4.59692] (0.95 confidence) | 0.989833 |
| d3/faer&#95;lu&#95;solve | 157.652 ns [157.502, 158.003] (0.95 confidence) | 163.46 ns [162.918, 164.545] (0.95 confidence) | 0.964467 |
| d3/la&#95;stack&#95;lu&#95;solve | 9.93158 ns [9.91324, 9.95553] (0.95 confidence) | 10.1094 ns [10.0641, 10.135] (0.95 confidence) | 0.982413 |
| d3/nalgebra&#95;lu&#95;solve | 33.8755 ns [33.8276, 33.9858] (0.95 confidence) | 35.0464 ns [35.0098, 35.1322] (0.95 confidence) | 0.96659 |
| d32/faer&#95;lu&#95;solve | 2949.45 ns [2945.04, 2958.89] (0.95 confidence) | 2974.82 ns [2969.33, 2979.13] (0.95 confidence) | 0.991473 |
| d32/la&#95;stack&#95;lu&#95;solve | 2766.49 ns [2763.42, 2771.15] (0.95 confidence) | 2787.86 ns [2781.93, 2792.44] (0.95 confidence) | 0.992337 |
| d32/nalgebra&#95;lu&#95;solve | 2737.51 ns [2733.92, 2740.01] (0.95 confidence) | 2729.77 ns [2726.39, 2731.83] (0.95 confidence) | 1.00284 |
| d4/faer&#95;lu&#95;solve | 183.115 ns [182.666, 183.737] (0.95 confidence) | 189.329 ns [189.14, 189.77] (0.95 confidence) | 0.967177 |
| d4/la&#95;stack&#95;lu&#95;solve | 23.6092 ns [23.5848, 23.6336] (0.95 confidence) | 24.0133 ns [24.0048, 24.0359] (0.95 confidence) | 0.983172 |
| d4/nalgebra&#95;lu&#95;solve | 53.6451 ns [53.5895, 53.7487] (0.95 confidence) | 55.4432 ns [55.2079, 55.6073] (0.95 confidence) | 0.967569 |
| d5/faer&#95;lu&#95;solve | 250.156 ns [249.643, 250.61] (0.95 confidence) | 252.753 ns [252.337, 253.319] (0.95 confidence) | 0.989726 |
| d5/la&#95;stack&#95;lu&#95;solve | 47.491 ns [47.1864, 47.7292] (0.95 confidence) | 47.8723 ns [47.6686, 48.1696] (0.95 confidence) | 0.992036 |
| d5/nalgebra&#95;lu&#95;solve | 75.5898 ns [75.4762, 75.7332] (0.95 confidence) | 78.1378 ns [77.9205, 78.4518] (0.95 confidence) | 0.967392 |
| d64/faer&#95;lu&#95;solve | 12829.2 ns [12815.6, 12858.5] (0.95 confidence) | 13017.3 ns [12958.9, 13039] (0.95 confidence) | 0.985549 |
| d64/la&#95;stack&#95;lu&#95;solve | 18629.1 ns [18599.8, 18663.8] (0.95 confidence) | 18063.6 ns [18036.7, 18106.7] (0.95 confidence) | 1.03131 |
| d64/nalgebra&#95;lu&#95;solve | 14441.1 ns [14397.9, 14483.2] (0.95 confidence) | 14343.6 ns [14332.8, 14361.8] (0.95 confidence) | 1.0068 |
| d8/faer&#95;lu&#95;solve | 334.227 ns [333.627, 335.329] (0.95 confidence) | 338.191 ns [337.146, 338.919] (0.95 confidence) | 0.988277 |
| d8/la&#95;stack&#95;lu&#95;solve | 131.221 ns [130.934, 131.421] (0.95 confidence) | 133.256 ns [133.148, 133.486] (0.95 confidence) | 0.984728 |
| d8/nalgebra&#95;lu&#95;solve | 167.883 ns [167.682, 168.585] (0.95 confidence) | 168.299 ns [168.089, 168.876] (0.95 confidence) | 0.997531 |

## column-basis

Semantic contract: la-stack-column-basis-v1. Collection session: `04f8194525b24da4a67a76c579c6fb93`.
Collection interval (Unix ns): 1794387719673835000–1794388833760030000.
Acquisition: newly collected.

- baseline: available — inside declared release interval; release v0.5.2; revision `eb5e44e52399401311e676faad590e6cce1610e5`.
  Source SHA-256: `bdd47c0e9bc3126dc3a32a98e0f1a688aaa52e7c6a3868009f1adfb0735f367e`; harness SHA-256:
  `0fedb76c77086b1cf8a698034ac4a1f3539fe06feda597672e2d11f01704f9cb`; scientific gate: passed.
- current: available — inside declared release interval; release v0.5.3; revision `be9277bf04400dccbe1115da69c7fda2141c4a39`.
  Source SHA-256: `793169b4c5085624d7c0e773f1d28504c017e9b78332f250150d19df87391a32`; harness SHA-256:
  `0fedb76c77086b1cf8a698034ac4a1f3539fe06feda597672e2d11f01704f9cb`; scientific gate: passed.

### Mean (ns)

| Workload | Baseline | Current | Baseline / current |
| --- | --- | --- | --- |
| column&#95;basis/2x2/full&#95;rank | 4.95857 ns [4.94946, 4.96799] (0.95 confidence) | 4.99445 ns [4.98573, 5.00323] (0.95 confidence) | 0.992817 |
| column&#95;basis/2x2/large | 4.98151 ns [4.97293, 4.99025] (0.95 confidence) | 5.00686 ns [4.9981, 5.01589] (0.95 confidence) | 0.994936 |
| column&#95;basis/2x2/near&#95;dependent | 5.04413 ns [5.02498, 5.07084] (0.95 confidence) | 5.0563 ns [5.04757, 5.06471] (0.95 confidence) | 0.997593 |
| column&#95;basis/2x2/rank&#95;one | 4.78021 ns [4.77374, 4.78664] (0.95 confidence) | 4.78727 ns [4.77683, 4.79761] (0.95 confidence) | 0.998525 |
| column&#95;basis/2x2/small | 4.96337 ns [4.95738, 4.96946] (0.95 confidence) | 4.9956 ns [4.98903, 5.00215] (0.95 confidence) | 0.993547 |
| column&#95;basis/2x3/shared&#95;face | 10.3585 ns [10.3239, 10.3935] (0.95 confidence) | 10.4999 ns [10.4678, 10.5341] (0.95 confidence) | 0.98653 |
| column&#95;basis/3x3/full&#95;rank | 15.2904 ns [14.6792, 15.9768] (0.95 confidence) | 15.4427 ns [14.8236, 16.1161] (0.95 confidence) | 0.990137 |
| column&#95;basis/3x3/large | 15.3024 ns [14.6966, 15.9841] (0.95 confidence) | 15.423 ns [14.7972, 16.1092] (0.95 confidence) | 0.992182 |
| column&#95;basis/3x3/near&#95;dependent | 15.365 ns [14.7614, 16.0365] (0.95 confidence) | 15.3527 ns [14.7428, 16.0341] (0.95 confidence) | 1.0008 |
| column&#95;basis/3x3/rank&#95;one | 9.36361 ns [9.34725, 9.38186] (0.95 confidence) | 9.31781 ns [9.25137, 9.38276] (0.95 confidence) | 1.00492 |
| column&#95;basis/3x3/small | 15.3093 ns [14.7011, 15.9909] (0.95 confidence) | 15.3122 ns [14.701, 15.9965] (0.95 confidence) | 0.999808 |
| column&#95;basis/3x4/shared&#95;face | 18.2592 ns [17.9073, 18.7222] (0.95 confidence) | 18.058 ns [17.9556, 18.2225] (0.95 confidence) | 1.01114 |
| column&#95;basis/4x4/full&#95;rank | 30.4948 ns [30.4303, 30.5586] (0.95 confidence) | 30.6305 ns [30.5463, 30.7148] (0.95 confidence) | 0.995569 |
| column&#95;basis/4x4/large | 30.6106 ns [30.55, 30.671] (0.95 confidence) | 30.505 ns [30.4198, 30.5891] (0.95 confidence) | 1.00346 |
| column&#95;basis/4x4/near&#95;dependent | 30.64 ns [30.5667, 30.7141] (0.95 confidence) | 30.7826 ns [30.7019, 30.8635] (0.95 confidence) | 0.995368 |
| column&#95;basis/4x4/rank&#95;one | 14.6506 ns [14.6195, 14.6844] (0.95 confidence) | 14.8209 ns [14.7667, 14.8776] (0.95 confidence) | 0.988506 |
| column&#95;basis/4x4/small | 30.4533 ns [30.3814, 30.5251] (0.95 confidence) | 30.5759 ns [30.4808, 30.679] (0.95 confidence) | 0.995989 |
| column&#95;basis/4x5/shared&#95;face | 40.8523 ns [40.374, 41.3276] (0.95 confidence) | 41.1964 ns [40.6945, 41.6869] (0.95 confidence) | 0.991649 |
| column&#95;basis/4x6/full&#95;rank | 36.622 ns [36.5532, 36.6899] (0.95 confidence) | 37.1737 ns [37.1003, 37.2452] (0.95 confidence) | 0.985158 |
| column&#95;basis/4x6/large | 36.4041 ns [36.3257, 36.4845] (0.95 confidence) | 37.3302 ns [37.2617, 37.4006] (0.95 confidence) | 0.975192 |
| column&#95;basis/4x6/near&#95;dependent | 36.7806 ns [36.7067, 36.8578] (0.95 confidence) | 37.3993 ns [37.2968, 37.5025] (0.95 confidence) | 0.983457 |
| column&#95;basis/4x6/rank&#95;one | 20.7955 ns [20.7711, 20.8201] (0.95 confidence) | 21.1325 ns [21.0928, 21.1718] (0.95 confidence) | 0.984053 |
| column&#95;basis/4x6/small | 36.7549 ns [36.6923, 36.8175] (0.95 confidence) | 37.391 ns [37.3296, 37.455] (0.95 confidence) | 0.98299 |
| column&#95;basis/5x2/full&#95;rank | 13.4934 ns [13.471, 13.5155] (0.95 confidence) | 13.4423 ns [13.4125, 13.4735] (0.95 confidence) | 1.0038 |
| column&#95;basis/5x2/large | 13.5124 ns [13.4858, 13.5399] (0.95 confidence) | 13.6001 ns [13.5634, 13.6372] (0.95 confidence) | 0.99355 |
| column&#95;basis/5x2/near&#95;dependent | 13.493 ns [13.4623, 13.5327] (0.95 confidence) | 13.4692 ns [13.4327, 13.5054] (0.95 confidence) | 1.00177 |
| column&#95;basis/5x2/rank&#95;one | 11.0391 ns [11.0145, 11.0648] (0.95 confidence) | 11.3118 ns [11.2215, 11.4151] (0.95 confidence) | 0.975894 |
| column&#95;basis/5x2/small | 13.4317 ns [13.401, 13.4625] (0.95 confidence) | 13.7135 ns [13.6815, 13.746] (0.95 confidence) | 0.979451 |
| column&#95;basis/5x5/full&#95;rank | 55.0425 ns [54.2222, 55.8859] (0.95 confidence) | 55.1392 ns [54.3182, 55.9923] (0.95 confidence) | 0.998246 |
| column&#95;basis/5x5/large | 54.9789 ns [54.162, 55.8202] (0.95 confidence) | 57.8369 ns [56.8538, 58.8504] (0.95 confidence) | 0.950585 |
| column&#95;basis/5x5/near&#95;dependent | 54.5977 ns [53.8586, 55.3624] (0.95 confidence) | 54.682 ns [53.9417, 55.4407] (0.95 confidence) | 0.998459 |
| column&#95;basis/5x5/rank&#95;one | 23.9293 ns [23.8999, 23.9616] (0.95 confidence) | 24.2538 ns [24.1941, 24.3154] (0.95 confidence) | 0.98662 |
| column&#95;basis/5x5/small | 54.9416 ns [54.1288, 55.7853] (0.95 confidence) | 56.6418 ns [55.7997, 57.5207] (0.95 confidence) | 0.969984 |
| column&#95;basis/5x6/shared&#95;face | 62.6146 ns [62.4525, 62.7859] (0.95 confidence) | 63.5088 ns [63.2853, 63.7463] (0.95 confidence) | 0.98592 |
| column&#95;basis/5x8/full&#95;rank | 62.9981 ns [62.9215, 63.0804] (0.95 confidence) | 63.9928 ns [63.8701, 64.1228] (0.95 confidence) | 0.984456 |
| column&#95;basis/5x8/large | 63.1311 ns [63.0355, 63.2322] (0.95 confidence) | 63.3481 ns [63.288, 63.4113] (0.95 confidence) | 0.996575 |
| column&#95;basis/5x8/near&#95;dependent | 62.8182 ns [62.742, 62.902] (0.95 confidence) | 63.1838 ns [63.119, 63.2535] (0.95 confidence) | 0.994213 |
| column&#95;basis/5x8/rank&#95;one | 38.684 ns [38.0735, 39.313] (0.95 confidence) | 39.6855 ns [39.0429, 40.332] (0.95 confidence) | 0.974765 |
| column&#95;basis/5x8/small | 63.0828 ns [63.0099, 63.1563] (0.95 confidence) | 63.3481 ns [63.2795, 63.4248] (0.95 confidence) | 0.995812 |
| column&#95;basis/6x10/full&#95;rank | 97.7059 ns [97.5639, 97.8583] (0.95 confidence) | 98.6351 ns [98.4825, 98.7998] (0.95 confidence) | 0.99058 |
| column&#95;basis/6x10/large | 97.3925 ns [97.268, 97.5251] (0.95 confidence) | 98.5616 ns [98.4301, 98.7064] (0.95 confidence) | 0.988139 |
| column&#95;basis/6x10/near&#95;dependent | 97.5668 ns [97.4457, 97.6957] (0.95 confidence) | 99.0418 ns [98.8748, 99.214] (0.95 confidence) | 0.985107 |
| column&#95;basis/6x10/rank&#95;one | 56.3396 ns [55.8984, 56.8133] (0.95 confidence) | 57.8989 ns [57.3823, 58.42] (0.95 confidence) | 0.973069 |
| column&#95;basis/6x10/small | 97.3596 ns [97.2333, 97.49] (0.95 confidence) | 99.0984 ns [98.8937, 99.3065] (0.95 confidence) | 0.982454 |
| column&#95;basis/6x7/shared&#95;face | 100.151 ns [99.8376, 100.457] (0.95 confidence) | 98.2549 ns [97.5946, 99.0073] (0.95 confidence) | 1.01929 |
| column&#95;basis/7x12/full&#95;rank | 143.048 ns [142.855, 143.263] (0.95 confidence) | 144.49 ns [144.357, 144.652] (0.95 confidence) | 0.990019 |
| column&#95;basis/7x12/large | 143.082 ns [142.901, 143.29] (0.95 confidence) | 144.486 ns [144.303, 144.697] (0.95 confidence) | 0.99028 |
| column&#95;basis/7x12/near&#95;dependent | 142.895 ns [142.758, 143.045] (0.95 confidence) | 144.264 ns [144.155, 144.387] (0.95 confidence) | 0.990514 |
| column&#95;basis/7x12/rank&#95;one | 79.2708 ns [79.2027, 79.3459] (0.95 confidence) | 80.0876 ns [80.0346, 80.1437] (0.95 confidence) | 0.9898 |
| column&#95;basis/7x12/small | 143.036 ns [142.893, 143.189] (0.95 confidence) | 144.445 ns [144.305, 144.603] (0.95 confidence) | 0.990247 |
| column&#95;basis/8x14/full&#95;rank | 198.64 ns [198.406, 198.9] (0.95 confidence) | 201.456 ns [201.106, 201.837] (0.95 confidence) | 0.986022 |
| column&#95;basis/8x14/large | 198.65 ns [198.44, 198.874] (0.95 confidence) | 201.076 ns [200.807, 201.4] (0.95 confidence) | 0.987935 |
| column&#95;basis/8x14/near&#95;dependent | 198.912 ns [198.651, 199.191] (0.95 confidence) | 200.419 ns [200.182, 200.68] (0.95 confidence) | 0.99248 |
| column&#95;basis/8x14/rank&#95;one | 105.297 ns [105.161, 105.452] (0.95 confidence) | 106.495 ns [106.351, 106.653] (0.95 confidence) | 0.988754 |
| column&#95;basis/8x14/small | 198.999 ns [198.75, 199.262] (0.95 confidence) | 200.877 ns [200.65, 201.125] (0.95 confidence) | 0.990653 |

### Median (ns)

| Workload | Baseline | Current | Baseline / current |
| --- | --- | --- | --- |
| column&#95;basis/2x2/full&#95;rank | 4.94999 ns [4.94137, 4.96854] (0.95 confidence) | 4.98621 ns [4.97854, 4.99824] (0.95 confidence) | 0.992737 |
| column&#95;basis/2x2/large | 4.98642 ns [4.97372, 4.99642] (0.95 confidence) | 5.00171 ns [4.99105, 5.00785] (0.95 confidence) | 0.996943 |
| column&#95;basis/2x2/near&#95;dependent | 5.03073 ns [5.01487, 5.04577] (0.95 confidence) | 5.05517 ns [5.05127, 5.06435] (0.95 confidence) | 0.995166 |
| column&#95;basis/2x2/rank&#95;one | 4.78449 ns [4.7796, 4.78924] (0.95 confidence) | 4.7967 ns [4.77312, 4.81859] (0.95 confidence) | 0.997453 |
| column&#95;basis/2x2/small | 4.96134 ns [4.95566, 4.9709] (0.95 confidence) | 4.9943 ns [4.98693, 5.00153] (0.95 confidence) | 0.993401 |
| column&#95;basis/2x3/shared&#95;face | 10.3229 ns [10.2815, 10.383] (0.95 confidence) | 10.438 ns [10.4008, 10.4693] (0.95 confidence) | 0.988978 |
| column&#95;basis/3x3/full&#95;rank | 13.6927 ns [13.6684, 13.7119] (0.95 confidence) | 13.7592 ns [13.7231, 13.8992] (0.95 confidence) | 0.995166 |
| column&#95;basis/3x3/large | 13.6862 ns [13.669, 13.7075] (0.95 confidence) | 13.7725 ns [13.7288, 13.843] (0.95 confidence) | 0.993739 |
| column&#95;basis/3x3/near&#95;dependent | 13.7446 ns [13.7067, 13.797] (0.95 confidence) | 13.7423 ns [13.7142, 13.7698] (0.95 confidence) | 1.00017 |
| column&#95;basis/3x3/rank&#95;one | 9.33828 ns [9.32787, 9.35347] (0.95 confidence) | 9.36702 ns [9.34791, 9.38098] (0.95 confidence) | 0.996932 |
| column&#95;basis/3x3/small | 13.6921 ns [13.6759, 13.7079] (0.95 confidence) | 13.6985 ns [13.6748, 13.7263] (0.95 confidence) | 0.999534 |
| column&#95;basis/3x4/shared&#95;face | 17.8349 ns [17.8187, 17.8606] (0.95 confidence) | 17.908 ns [17.8878, 17.9406] (0.95 confidence) | 0.995922 |
| column&#95;basis/4x4/full&#95;rank | 30.5398 ns [30.4437, 30.5982] (0.95 confidence) | 30.6733 ns [30.5705, 30.7514] (0.95 confidence) | 0.99565 |
| column&#95;basis/4x4/large | 30.657 ns [30.6095, 30.7013] (0.95 confidence) | 30.5197 ns [30.3916, 30.6299] (0.95 confidence) | 1.0045 |
| column&#95;basis/4x4/near&#95;dependent | 30.6852 ns [30.6117, 30.7677] (0.95 confidence) | 30.8269 ns [30.7238, 30.9033] (0.95 confidence) | 0.995403 |
| column&#95;basis/4x4/rank&#95;one | 14.6553 ns [14.6301, 14.6947] (0.95 confidence) | 14.8005 ns [14.7482, 14.841] (0.95 confidence) | 0.990195 |
| column&#95;basis/4x4/small | 30.5041 ns [30.4184, 30.5523] (0.95 confidence) | 30.5008 ns [30.4009, 30.6735] (0.95 confidence) | 1.00011 |
| column&#95;basis/4x5/shared&#95;face | 41.4533 ns [39.9637, 42.5294] (0.95 confidence) | 42.6174 ns [41.236, 43.0788] (0.95 confidence) | 0.972686 |
| column&#95;basis/4x6/full&#95;rank | 36.6456 ns [36.5931, 36.7164] (0.95 confidence) | 37.1871 ns [37.123, 37.2802] (0.95 confidence) | 0.985439 |
| column&#95;basis/4x6/large | 36.3472 ns [36.2304, 36.5052] (0.95 confidence) | 37.2435 ns [37.1764, 37.3334] (0.95 confidence) | 0.975936 |
| column&#95;basis/4x6/near&#95;dependent | 36.7602 ns [36.6719, 36.7961] (0.95 confidence) | 37.3672 ns [37.2159, 37.4571] (0.95 confidence) | 0.983755 |
| column&#95;basis/4x6/rank&#95;one | 20.7941 ns [20.7528, 20.8233] (0.95 confidence) | 21.1553 ns [21.0636, 21.1897] (0.95 confidence) | 0.982925 |
| column&#95;basis/4x6/small | 36.7483 ns [36.658, 36.8251] (0.95 confidence) | 37.3533 ns [37.3125, 37.4285] (0.95 confidence) | 0.983803 |
| column&#95;basis/5x2/full&#95;rank | 13.4857 ns [13.465, 13.534] (0.95 confidence) | 13.4284 ns [13.3931, 13.4756] (0.95 confidence) | 1.00427 |
| column&#95;basis/5x2/large | 13.4882 ns [13.4657, 13.5042] (0.95 confidence) | 13.5951 ns [13.5672, 13.634] (0.95 confidence) | 0.992135 |
| column&#95;basis/5x2/near&#95;dependent | 13.4721 ns [13.45, 13.4885] (0.95 confidence) | 13.4846 ns [13.405, 13.5277] (0.95 confidence) | 0.999071 |
| column&#95;basis/5x2/rank&#95;one | 11.033 ns [11.0181, 11.0549] (0.95 confidence) | 11.1766 ns [11.0957, 11.324] (0.95 confidence) | 0.987149 |
| column&#95;basis/5x2/small | 13.4355 ns [13.4041, 13.4666] (0.95 confidence) | 13.7079 ns [13.6627, 13.7499] (0.95 confidence) | 0.980132 |
| column&#95;basis/5x5/full&#95;rank | 52.4387 ns [52.3653, 52.5899] (0.95 confidence) | 52.5759 ns [52.4207, 52.8182] (0.95 confidence) | 0.99739 |
| column&#95;basis/5x5/large | 52.3718 ns [52.1968, 52.6746] (0.95 confidence) | 57.1276 ns [54.4495, 58.8424] (0.95 confidence) | 0.916752 |
| column&#95;basis/5x5/near&#95;dependent | 52.6094 ns [52.44, 53.2882] (0.95 confidence) | 52.5045 ns [52.4037, 52.8806] (0.95 confidence) | 1.002 |
| column&#95;basis/5x5/rank&#95;one | 23.9258 ns [23.9011, 23.9488] (0.95 confidence) | 24.3013 ns [24.1774, 24.3635] (0.95 confidence) | 0.984548 |
| column&#95;basis/5x5/small | 52.2136 ns [52.0536, 52.5622] (0.95 confidence) | 54.4598 ns [54.1695, 54.8149] (0.95 confidence) | 0.958755 |
| column&#95;basis/5x6/shared&#95;face | 62.5242 ns [62.3277, 62.7171] (0.95 confidence) | 63.2916 ns [63.0003, 63.6337] (0.95 confidence) | 0.987875 |
| column&#95;basis/5x8/full&#95;rank | 62.9484 ns [62.8439, 63.0524] (0.95 confidence) | 63.9084 ns [63.8054, 63.9959] (0.95 confidence) | 0.984978 |
| column&#95;basis/5x8/large | 63.057 ns [62.9259, 63.1571] (0.95 confidence) | 63.3914 ns [63.2935, 63.4258] (0.95 confidence) | 0.994726 |
| column&#95;basis/5x8/near&#95;dependent | 62.6851 ns [62.6462, 62.7628] (0.95 confidence) | 63.1368 ns [63.0954, 63.208] (0.95 confidence) | 0.992846 |
| column&#95;basis/5x8/rank&#95;one | 38.294 ns [37.541, 38.7991] (0.95 confidence) | 39.6516 ns [38.3108, 40.5396] (0.95 confidence) | 0.965762 |
| column&#95;basis/5x8/small | 63.0795 ns [62.9598, 63.1363] (0.95 confidence) | 63.3458 ns [63.2923, 63.3744] (0.95 confidence) | 0.995797 |
| column&#95;basis/6x10/full&#95;rank | 97.5746 ns [97.4285, 97.6794] (0.95 confidence) | 98.3174 ns [98.2363, 98.4038] (0.95 confidence) | 0.992445 |
| column&#95;basis/6x10/large | 97.2404 ns [97.1781, 97.2845] (0.95 confidence) | 98.3635 ns [98.2184, 98.4636] (0.95 confidence) | 0.988582 |
| column&#95;basis/6x10/near&#95;dependent | 97.4051 ns [97.3029, 97.5137] (0.95 confidence) | 98.8136 ns [98.6663, 99.0853] (0.95 confidence) | 0.985746 |
| column&#95;basis/6x10/rank&#95;one | 55.5749 ns [55.3204, 55.7977] (0.95 confidence) | 57.619 ns [56.6751, 58.5022] (0.95 confidence) | 0.964524 |
| column&#95;basis/6x10/small | 97.2463 ns [97.0602, 97.4376] (0.95 confidence) | 99.1046 ns [98.4609, 99.5201] (0.95 confidence) | 0.981249 |
| column&#95;basis/6x7/shared&#95;face | 100.292 ns [100.065, 100.582] (0.95 confidence) | 97.4249 ns [96.85, 97.8915] (0.95 confidence) | 1.02943 |
| column&#95;basis/7x12/full&#95;rank | 142.912 ns [142.689, 143.05] (0.95 confidence) | 144.326 ns [144.251, 144.45] (0.95 confidence) | 0.990207 |
| column&#95;basis/7x12/large | 142.901 ns [142.723, 143.043] (0.95 confidence) | 144.208 ns [144.055, 144.3] (0.95 confidence) | 0.990935 |
| column&#95;basis/7x12/near&#95;dependent | 142.845 ns [142.677, 143.072] (0.95 confidence) | 144.132 ns [144.066, 144.2] (0.95 confidence) | 0.991066 |
| column&#95;basis/7x12/rank&#95;one | 79.2412 ns [79.1535, 79.3086] (0.95 confidence) | 80.0584 ns [80.0289, 80.0977] (0.95 confidence) | 0.989793 |
| column&#95;basis/7x12/small | 142.923 ns [142.812, 143.072] (0.95 confidence) | 144.246 ns [144.167, 144.329] (0.95 confidence) | 0.99083 |
| column&#95;basis/8x14/full&#95;rank | 198.367 ns [198.212, 198.645] (0.95 confidence) | 200.722 ns [200.511, 201.312] (0.95 confidence) | 0.988267 |
| column&#95;basis/8x14/large | 198.46 ns [198.263, 198.7] (0.95 confidence) | 200.855 ns [200.58, 201.07] (0.95 confidence) | 0.988077 |
| column&#95;basis/8x14/near&#95;dependent | 198.612 ns [198.326, 198.795] (0.95 confidence) | 200.221 ns [200.017, 200.501] (0.95 confidence) | 0.991963 |
| column&#95;basis/8x14/rank&#95;one | 105.103 ns [104.98, 105.172] (0.95 confidence) | 106.267 ns [106.222, 106.349] (0.95 confidence) | 0.98904 |
| column&#95;basis/8x14/small | 198.684 ns [198.455, 199.044] (0.95 confidence) | 200.618 ns [200.468, 200.927] (0.95 confidence) | 0.990364 |

## interval-solve

Semantic contract: la-stack-interval-solve-v1. Collection session: `fbdfe7817d0845f993432cd0c6ef1f5b`.
Collection interval (Unix ns): 1794388835057450000–1794389581938382000.
Acquisition: newly collected.

- baseline: available — inside declared release interval; release v0.5.2; revision `eb5e44e52399401311e676faad590e6cce1610e5`.
  Source SHA-256: `bdd47c0e9bc3126dc3a32a98e0f1a688aaa52e7c6a3868009f1adfb0735f367e`; harness SHA-256:
  `0fedb76c77086b1cf8a698034ac4a1f3539fe06feda597672e2d11f01704f9cb`; scientific gate: passed.
- current: available — inside declared release interval; release v0.5.3; revision `be9277bf04400dccbe1115da69c7fda2141c4a39`.
  Source SHA-256: `793169b4c5085624d7c0e773f1d28504c017e9b78332f250150d19df87391a32`; harness SHA-256:
  `0fedb76c77086b1cf8a698034ac4a1f3539fe06feda597672e2d11f01704f9cb`; scientific gate: passed.

### Mean (ns)

| Workload | Baseline | Current | Baseline / current |
| --- | --- | --- | --- |
| interval&#95;division/negative | 5.30969 ns [5.30091, 5.31876] (0.95 confidence) | 5.26081 ns [5.25094, 5.27156] (0.95 confidence) | 1.00929 |
| interval&#95;division/positive | 5.32575 ns [5.31602, 5.336] (0.95 confidence) | 5.29187 ns [5.28, 5.30426] (0.95 confidence) | 1.0064 |
| interval&#95;division/reciprocal | 10.074 ns [10.0522, 10.0956] (0.95 confidence) | 10.0423 ns [10.0289, 10.0561] (0.95 confidence) | 1.00316 |
| interval&#95;division/subnormal | 3.88365 ns [3.8777, 3.8899] (0.95 confidence) | 3.84365 ns [3.8371, 3.85033] (0.95 confidence) | 1.01041 |
| interval&#95;division/tiny&#95;ratio | 6.1485 ns [6.138, 6.15913] (0.95 confidence) | 6.11212 ns [6.10288, 6.12177] (0.95 confidence) | 1.00595 |
| interval&#95;solve&#95;d1/dense | 232.956 ns [231.334, 235.272] (0.95 confidence) | 230.35 ns [229.931, 230.779] (0.95 confidence) | 1.01131 |
| interval&#95;solve&#95;d1/near&#95;singular | 298.006 ns [296.545, 300.552] (0.95 confidence) | 296.05 ns [294.345, 298.631] (0.95 confidence) | 1.0066 |
| interval&#95;solve&#95;d1/range&#95;exhausted | 157.471 ns [156.981, 158.031] (0.95 confidence) | 155.548 ns [155.06, 156.105] (0.95 confidence) | 1.01236 |
| interval&#95;solve&#95;d1/singular | 53.3616 ns [53.1276, 53.6197] (0.95 confidence) | 52.7329 ns [52.5021, 52.9818] (0.95 confidence) | 1.01192 |
| interval&#95;solve&#95;d1/uncertain&#95;dense | 380.465 ns [379.882, 381.122] (0.95 confidence) | 376.305 ns [375.956, 376.677] (0.95 confidence) | 1.01106 |
| interval&#95;solve&#95;d2/dense | 3101.97 ns [3099.02, 3105.33] (0.95 confidence) | 3006.66 ns [3000.72, 3013.53] (0.95 confidence) | 1.0317 |
| interval&#95;solve&#95;d2/near&#95;singular | 1678.01 ns [1675.28, 1680.87] (0.95 confidence) | 1653.09 ns [1650.63, 1655.53] (0.95 confidence) | 1.01507 |
| interval&#95;solve&#95;d2/range&#95;exhausted | 258.838 ns [258.395, 259.354] (0.95 confidence) | 255.872 ns [255.67, 256.094] (0.95 confidence) | 1.01159 |
| interval&#95;solve&#95;d2/singular | 172.816 ns [172.533, 173.115] (0.95 confidence) | 170.139 ns [169.937, 170.346] (0.95 confidence) | 1.01573 |
| interval&#95;solve&#95;d2/uncertain&#95;dense | 4244.66 ns [4241.82, 4248.2] (0.95 confidence) | 4123.7 ns [4116.55, 4130.93] (0.95 confidence) | 1.02933 |
| interval&#95;solve&#95;d3/dense | 7896.77 ns [7883.53, 7911.41] (0.95 confidence) | 7655.91 ns [7643.94, 7667.79] (0.95 confidence) | 1.03146 |
| interval&#95;solve&#95;d3/near&#95;singular | 2010.78 ns [2006.84, 2015.02] (0.95 confidence) | 1962.06 ns [1959.37, 1964.76] (0.95 confidence) | 1.02483 |
| interval&#95;solve&#95;d3/range&#95;exhausted | 418.507 ns [418.129, 418.924] (0.95 confidence) | 413.265 ns [412.875, 413.736] (0.95 confidence) | 1.01268 |
| interval&#95;solve&#95;d3/singular | 346.011 ns [343.138, 348.992] (0.95 confidence) | 326.579 ns [326.104, 327.059] (0.95 confidence) | 1.0595 |
| interval&#95;solve&#95;d3/uncertain&#95;dense | 9637.69 ns [9620.95, 9654.52] (0.95 confidence) | 9334.83 ns [9320.06, 9349.64] (0.95 confidence) | 1.03244 |
| interval&#95;solve&#95;d4/dense | 15241.8 ns [15215.6, 15270.7] (0.95 confidence) | 14958.2 ns [14929.4, 14987] (0.95 confidence) | 1.01896 |
| interval&#95;solve&#95;d4/near&#95;singular | 2673.6 ns [2670.44, 2677] (0.95 confidence) | 2634.32 ns [2631.14, 2637.75] (0.95 confidence) | 1.01491 |
| interval&#95;solve&#95;d4/range&#95;exhausted | 773.36 ns [772.823, 774.002] (0.95 confidence) | 764.069 ns [763.243, 764.993] (0.95 confidence) | 1.01216 |
| interval&#95;solve&#95;d4/singular | 592.9 ns [592.432, 593.451] (0.95 confidence) | 584.809 ns [584.11, 585.632] (0.95 confidence) | 1.01384 |
| interval&#95;solve&#95;d4/uncertain&#95;dense | 18049.9 ns [18024.5, 18077.4] (0.95 confidence) | 17677.5 ns [17644.9, 17712.2] (0.95 confidence) | 1.02107 |
| interval&#95;solve&#95;d5/dense | 25918.6 ns [25885.4, 25954.6] (0.95 confidence) | 25604 ns [25566.3, 25641.2] (0.95 confidence) | 1.01229 |
| interval&#95;solve&#95;d5/near&#95;singular | 3514.71 ns [3510.36, 3519.17] (0.95 confidence) | 3477.79 ns [3473.02, 3482.94] (0.95 confidence) | 1.01061 |
| interval&#95;solve&#95;d5/range&#95;exhausted | 1276.89 ns [1275.4, 1278.42] (0.95 confidence) | 1263.21 ns [1261.7, 1264.72] (0.95 confidence) | 1.01083 |
| interval&#95;solve&#95;d5/singular | 1039.38 ns [1037.99, 1040.79] (0.95 confidence) | 1026.87 ns [1025.3, 1028.45] (0.95 confidence) | 1.01218 |
| interval&#95;solve&#95;d5/uncertain&#95;dense | 30007.1 ns [29961.4, 30054.4] (0.95 confidence) | 29443.8 ns [29377.5, 29511.2] (0.95 confidence) | 1.01913 |
| interval&#95;solve&#95;d6/dense | 41458.5 ns [41315.3, 41605.7] (0.95 confidence) | 40604.3 ns [40548.9, 40658.5] (0.95 confidence) | 1.02104 |
| interval&#95;solve&#95;d6/near&#95;singular | 4446.66 ns [4441.92, 4451.74] (0.95 confidence) | 4403.03 ns [4397.07, 4409.66] (0.95 confidence) | 1.00991 |
| interval&#95;solve&#95;d6/range&#95;exhausted | 1657 ns [1646.89, 1668.12] (0.95 confidence) | 1600.65 ns [1598.68, 1602.78] (0.95 confidence) | 1.0352 |
| interval&#95;solve&#95;d6/singular | 1382.95 ns [1374.61, 1391.93] (0.95 confidence) | 1332.87 ns [1331.05, 1334.74] (0.95 confidence) | 1.03757 |
| interval&#95;solve&#95;d6/uncertain&#95;dense | 46902.5 ns [46774.2, 47030] (0.95 confidence) | 46071.5 ns [45903.7, 46241.3] (0.95 confidence) | 1.01804 |

### Median (ns)

| Workload | Baseline | Current | Baseline / current |
| --- | --- | --- | --- |
| interval&#95;division/negative | 5.28774 ns [5.28512, 5.31845] (0.95 confidence) | 5.24517 ns [5.2355, 5.26328] (0.95 confidence) | 1.00812 |
| interval&#95;division/positive | 5.30577 ns [5.29903, 5.31641] (0.95 confidence) | 5.28445 ns [5.2625, 5.29934] (0.95 confidence) | 1.00404 |
| interval&#95;division/reciprocal | 10.0965 ns [10.0624, 10.1381] (0.95 confidence) | 10.0329 ns [10.0121, 10.0545] (0.95 confidence) | 1.00635 |
| interval&#95;division/subnormal | 3.88278 ns [3.86635, 3.89045] (0.95 confidence) | 3.84179 ns [3.83227, 3.85065] (0.95 confidence) | 1.01067 |
| interval&#95;division/tiny&#95;ratio | 6.15399 ns [6.12828, 6.16302] (0.95 confidence) | 6.10978 ns [6.09622, 6.12098] (0.95 confidence) | 1.00724 |
| interval&#95;solve&#95;d1/dense | 232.229 ns [231.684, 233.456] (0.95 confidence) | 229.65 ns [229.367, 230.173] (0.95 confidence) | 1.01123 |
| interval&#95;solve&#95;d1/near&#95;singular | 296.514 ns [296.335, 297.206] (0.95 confidence) | 294.116 ns [293.682, 294.573] (0.95 confidence) | 1.00815 |
| interval&#95;solve&#95;d1/range&#95;exhausted | 156.663 ns [156.447, 157.43] (0.95 confidence) | 154.578 ns [154.335, 155.216] (0.95 confidence) | 1.01349 |
| interval&#95;solve&#95;d1/singular | 52.7512 ns [52.7414, 52.779] (0.95 confidence) | 52.3138 ns [52.233, 52.3992] (0.95 confidence) | 1.00836 |
| interval&#95;solve&#95;d1/uncertain&#95;dense | 379.617 ns [379.376, 380.006] (0.95 confidence) | 376.161 ns [375.904, 376.394] (0.95 confidence) | 1.00919 |
| interval&#95;solve&#95;d2/dense | 3098.14 ns [3096.54, 3100.1] (0.95 confidence) | 3005.07 ns [3000.19, 3009.36] (0.95 confidence) | 1.03097 |
| interval&#95;solve&#95;d2/near&#95;singular | 1676.05 ns [1672.76, 1678.35] (0.95 confidence) | 1653.73 ns [1650.42, 1656.16] (0.95 confidence) | 1.0135 |
| interval&#95;solve&#95;d2/range&#95;exhausted | 258.491 ns [258.188, 258.729] (0.95 confidence) | 255.656 ns [255.485, 255.853] (0.95 confidence) | 1.01109 |
| interval&#95;solve&#95;d2/singular | 172.249 ns [171.951, 172.373] (0.95 confidence) | 170.094 ns [169.893, 170.201] (0.95 confidence) | 1.01267 |
| interval&#95;solve&#95;d2/uncertain&#95;dense | 4241.5 ns [4240.08, 4243.21] (0.95 confidence) | 4125.37 ns [4120.12, 4131.38] (0.95 confidence) | 1.02815 |
| interval&#95;solve&#95;d3/dense | 7879.17 ns [7869.71, 7886.37] (0.95 confidence) | 7657.63 ns [7647.28, 7669.34] (0.95 confidence) | 1.02893 |
| interval&#95;solve&#95;d3/near&#95;singular | 2006.08 ns [2003.22, 2008.9] (0.95 confidence) | 1960.9 ns [1958.19, 1963.3] (0.95 confidence) | 1.02304 |
| interval&#95;solve&#95;d3/range&#95;exhausted | 417.655 ns [417.577, 418.042] (0.95 confidence) | 412.852 ns [412.674, 413.206] (0.95 confidence) | 1.01164 |
| interval&#95;solve&#95;d3/singular | 337.741 ns [336.606, 339.549] (0.95 confidence) | 326.805 ns [326.278, 327.076] (0.95 confidence) | 1.03346 |
| interval&#95;solve&#95;d3/uncertain&#95;dense | 9628.58 ns [9617.02, 9647.42] (0.95 confidence) | 9331.74 ns [9321.61, 9353.74] (0.95 confidence) | 1.03181 |
| interval&#95;solve&#95;d4/dense | 15217.7 ns [15197.2, 15258.4] (0.95 confidence) | 14945.3 ns [14914.6, 14985.7] (0.95 confidence) | 1.01823 |
| interval&#95;solve&#95;d4/near&#95;singular | 2668.05 ns [2666.75, 2672.44] (0.95 confidence) | 2630.68 ns [2629.02, 2633.48] (0.95 confidence) | 1.01421 |
| interval&#95;solve&#95;d4/range&#95;exhausted | 772.768 ns [772.388, 773.435] (0.95 confidence) | 763.352 ns [762.941, 764.071] (0.95 confidence) | 1.01234 |
| interval&#95;solve&#95;d4/singular | 591.923 ns [591.721, 592.464] (0.95 confidence) | 584.327 ns [584.066, 584.877] (0.95 confidence) | 1.013 |
| interval&#95;solve&#95;d4/uncertain&#95;dense | 18035.4 ns [17999.2, 18065.4] (0.95 confidence) | 17666.2 ns [17638.7, 17692.6] (0.95 confidence) | 1.0209 |
| interval&#95;solve&#95;d5/dense | 25894.8 ns [25871.7, 25906.3] (0.95 confidence) | 25597.3 ns [25530.2, 25663] (0.95 confidence) | 1.01162 |
| interval&#95;solve&#95;d5/near&#95;singular | 3512.64 ns [3505.29, 3519.79] (0.95 confidence) | 3474.35 ns [3468.98, 3481.17] (0.95 confidence) | 1.01102 |
| interval&#95;solve&#95;d5/range&#95;exhausted | 1278.72 ns [1274.09, 1279.47] (0.95 confidence) | 1264.42 ns [1261.15, 1265.95] (0.95 confidence) | 1.01132 |
| interval&#95;solve&#95;d5/singular | 1041.97 ns [1034.89, 1042.7] (0.95 confidence) | 1029.11 ns [1022.62, 1030.3] (0.95 confidence) | 1.0125 |
| interval&#95;solve&#95;d5/uncertain&#95;dense | 29911.1 ns [29883, 30010.4] (0.95 confidence) | 29378 ns [29316.5, 29469.2] (0.95 confidence) | 1.01815 |
| interval&#95;solve&#95;d6/dense | 41275.6 ns [41130.5, 41433.1] (0.95 confidence) | 40630.3 ns [40563.9, 40703.3] (0.95 confidence) | 1.01588 |
| interval&#95;solve&#95;d6/near&#95;singular | 4443.28 ns [4438.1, 4448.03] (0.95 confidence) | 4399.3 ns [4392.75, 4404.93] (0.95 confidence) | 1.01 |
| interval&#95;solve&#95;d6/range&#95;exhausted | 1639.82 ns [1633.02, 1645.28] (0.95 confidence) | 1598.76 ns [1597.42, 1601.33] (0.95 confidence) | 1.02568 |
| interval&#95;solve&#95;d6/singular | 1367.45 ns [1363.07, 1373.83] (0.95 confidence) | 1333.85 ns [1327.83, 1336.1] (0.95 confidence) | 1.02519 |
| interval&#95;solve&#95;d6/uncertain&#95;dense | 47079.6 ns [46759.5, 47196.2] (0.95 confidence) | 46213.1 ns [45539.9, 46433.9] (0.95 confidence) | 1.01875 |

## polynomial

Semantic contract: la-stack-polynomial-v1. Collection session: `983ca2ec2cab417d9c17695135d56ac4`.
Collection interval (Unix ns): 1794389582882250000–1794389775283775000.
Acquisition: newly collected.

- baseline: available — inside declared release interval; release v0.5.2; revision `eb5e44e52399401311e676faad590e6cce1610e5`.
  Source SHA-256: `bdd47c0e9bc3126dc3a32a98e0f1a688aaa52e7c6a3868009f1adfb0735f367e`; harness SHA-256:
  `0fedb76c77086b1cf8a698034ac4a1f3539fe06feda597672e2d11f01704f9cb`; scientific gate: passed.
- current: available — inside declared release interval; release v0.5.3; revision `be9277bf04400dccbe1115da69c7fda2141c4a39`.
  Source SHA-256: `793169b4c5085624d7c0e773f1d28504c017e9b78332f250150d19df87391a32`; harness SHA-256:
  `0fedb76c77086b1cf8a698034ac4a1f3539fe06feda597672e2d11f01704f9cb`; scientific gate: passed.

### Mean (ns)

<table>
<thead>
<tr>
<th>Workload</th>
<th>Baseline</th>
<th>Current</th>
<th>Baseline / current</th>
</tr>
</thead>
<tbody>
<tr>
<td>polynomial/coincident&#95;lifted/order4/determinant</td>
<td>13177.5 ns [13152.1, 13203.8] (0.95 confidence)</td>
<td>13276.8 ns [13244.2, 13309.1] (0.95 confidence)</td>
<td>0.992522</td>
</tr>
<tr>
<td>polynomial/coincident&#95;lifted/order5/determinant</td>
<td>105371 ns [105158, 105596] (0.95 confidence)</td>
<td>105648 ns [105398, 105897] (0.95 confidence)</td>
<td>0.997383</td>
</tr>
<tr>
<td>polynomial/coincident&#95;lifted/order6/determinant</td>
<td>911964 ns [910343, 913657] (0.95 confidence)</td>
<td>907804 ns [906325, 909331] (0.95 confidence)</td>
<td>1.00458</td>
</tr>
<tr>
<td>polynomial/coincident&#95;lifted/order7/determinant</td>
<td>8.47281e+06 ns [8.4609e+06, 8.48571e+06] (0.95 confidence)</td>
<td>8.50226e+06 ns [8.48922e+06, 8.51574e+06] (0.95 confidence)</td>
<td>0.996536</td>
</tr>
<tr>
<td>polynomial/coincident&#95;lifted/order8/determinant</td>
<td>8.57603e+07 ns [8.56128e+07, 8.59145e+07] (0.95 confidence)</td>
<td>8.54967e+07 ns [8.53594e+07, 8.56404e+07] (0.95 confidence)</td>
<td>1.00308</td>
</tr>
</tbody>
</table>

### Median (ns)

<table>
<thead>
<tr>
<th>Workload</th>
<th>Baseline</th>
<th>Current</th>
<th>Baseline / current</th>
</tr>
</thead>
<tbody>
<tr>
<td>polynomial/coincident&#95;lifted/order4/determinant</td>
<td>13137.6 ns [13114.2, 13164.8] (0.95 confidence)</td>
<td>13323.6 ns [13279.9, 13360.2] (0.95 confidence)</td>
<td>0.986038</td>
</tr>
<tr>
<td>polynomial/coincident&#95;lifted/order5/determinant</td>
<td>105061 ns [104749, 105261] (0.95 confidence)</td>
<td>105686 ns [105456, 105946] (0.95 confidence)</td>
<td>0.994086</td>
</tr>
<tr>
<td>polynomial/coincident&#95;lifted/order6/determinant</td>
<td>910025 ns [908453, 912123] (0.95 confidence)</td>
<td>908740 ns [904475, 909957] (0.95 confidence)</td>
<td>1.00141</td>
</tr>
<tr>
<td>polynomial/coincident&#95;lifted/order7/determinant</td>
<td>8.45672e+06 ns [8.44622e+06, 8.46157e+06] (0.95 confidence)</td>
<td>8.48525e+06 ns [8.47358e+06, 8.49841e+06] (0.95 confidence)</td>
<td>0.996638</td>
</tr>
<tr>
<td>polynomial/coincident&#95;lifted/order8/determinant</td>
<td>8.5454e+07 ns [8.5325e+07, 8.56305e+07] (0.95 confidence)</td>
<td>8.53886e+07 ns [8.51802e+07, 8.56451e+07] (0.95 confidence)</td>
<td>1.00077</td>
</tr>
</tbody>
</table>

## rational-dot

Semantic contract: la-stack-rational-dot-v1. Collection session: `a7af2e83c9854881a2c4e81500d7877a`.
Collection interval (Unix ns): 1794389776178361000–1794390867180518000.
Acquisition: newly collected.

- baseline: available — inside declared release interval; release v0.5.2; revision `eb5e44e52399401311e676faad590e6cce1610e5`.
  Source SHA-256: `bdd47c0e9bc3126dc3a32a98e0f1a688aaa52e7c6a3868009f1adfb0735f367e`; harness SHA-256:
  `0fedb76c77086b1cf8a698034ac4a1f3539fe06feda597672e2d11f01704f9cb`; scientific gate: passed.
- current: available — inside declared release interval; release v0.5.3; revision `be9277bf04400dccbe1115da69c7fda2141c4a39`.
  Source SHA-256: `793169b4c5085624d7c0e773f1d28504c017e9b78332f250150d19df87391a32`; harness SHA-256:
  `0fedb76c77086b1cf8a698034ac4a1f3539fe06feda597672e2d11f01704f9cb`; scientific gate: passed.

### Mean (ns)

| Workload | Baseline | Current | Baseline / current |
| --- | --- | --- | --- |
| rational&#95;dot/cancellation/d14/prepared | 1868.51 ns [1863.8, 1874.54] (0.95 confidence) | 1829.66 ns [1822.12, 1837.04] (0.95 confidence) | 1.02124 |
| rational&#95;dot/cancellation/d14/slice | 647.635 ns [646.068, 649.321] (0.95 confidence) | 646.245 ns [644.45, 648.042] (0.95 confidence) | 1.00215 |
| rational&#95;dot/cancellation/d2/prepared | 268.271 ns [267.937, 268.613] (0.95 confidence) | 272.118 ns [271.203, 272.983] (0.95 confidence) | 0.985864 |
| rational&#95;dot/cancellation/d2/slice | 96.8204 ns [96.7245, 96.9215] (0.95 confidence) | 98.0236 ns [97.9104, 98.1424] (0.95 confidence) | 0.987725 |
| rational&#95;dot/cancellation/d22/prepared | 2930.13 ns [2923.61, 2936.95] (0.95 confidence) | 2857.39 ns [2844.86, 2869.83] (0.95 confidence) | 1.02546 |
| rational&#95;dot/cancellation/d22/slice | 1000.68 ns [999.269, 1002.15] (0.95 confidence) | 994.448 ns [993.107, 995.82] (0.95 confidence) | 1.00627 |
| rational&#95;dot/cancellation/d3/prepared | 382.328 ns [381.714, 382.953] (0.95 confidence) | 385.75 ns [384.782, 386.69] (0.95 confidence) | 0.991129 |
| rational&#95;dot/cancellation/d3/slice | 126.724 ns [126.525, 126.927] (0.95 confidence) | 128.147 ns [127.987, 128.321] (0.95 confidence) | 0.988894 |
| rational&#95;dot/cancellation/d4/prepared | 532.839 ns [531.334, 534.403] (0.95 confidence) | 534.499 ns [532.481, 536.455] (0.95 confidence) | 0.996896 |
| rational&#95;dot/cancellation/d4/slice | 188.093 ns [187.662, 188.551] (0.95 confidence) | 189.219 ns [188.817, 189.66] (0.95 confidence) | 0.99405 |
| rational&#95;dot/cancellation/d5/prepared | 644.424 ns [643.254, 645.668] (0.95 confidence) | 654.454 ns [652.416, 656.546] (0.95 confidence) | 0.984674 |
| rational&#95;dot/cancellation/d5/slice | 213.324 ns [212.847, 213.85] (0.95 confidence) | 216.896 ns [216.473, 217.326] (0.95 confidence) | 0.983533 |
| rational&#95;dot/cancellation/d6/prepared | 808.091 ns [806.332, 809.876] (0.95 confidence) | 817.222 ns [811.345, 823.268] (0.95 confidence) | 0.988827 |
| rational&#95;dot/cancellation/d6/slice | 276.914 ns [275.952, 277.956] (0.95 confidence) | 278.598 ns [277.436, 279.885] (0.95 confidence) | 0.993954 |
| rational&#95;dot/cancellation/d7/prepared | 919.328 ns [916.963, 921.74] (0.95 confidence) | 920.054 ns [917.52, 922.658] (0.95 confidence) | 0.999211 |
| rational&#95;dot/cancellation/d7/slice | 307.59 ns [306.951, 308.286] (0.95 confidence) | 309.748 ns [308.664, 310.929] (0.95 confidence) | 0.993033 |
| rational&#95;dot/cancellation/d8/prepared | 1061.63 ns [1057.17, 1065.92] (0.95 confidence) | 1063.22 ns [1059.56, 1066.74] (0.95 confidence) | 0.998504 |
| rational&#95;dot/cancellation/d8/slice | 364.715 ns [364.033, 365.453] (0.95 confidence) | 372.366 ns [371.33, 373.419] (0.95 confidence) | 0.979451 |
| rational&#95;dot/large/d14/prepared | 2977.51 ns [2968.12, 2987.29] (0.95 confidence) | 3001.14 ns [2990.75, 3011.85] (0.95 confidence) | 0.992127 |
| rational&#95;dot/large/d14/slice | 1298.07 ns [1293.66, 1302.68] (0.95 confidence) | 1298.66 ns [1294.22, 1303.3] (0.95 confidence) | 0.999551 |
| rational&#95;dot/large/d2/prepared | 340.137 ns [339.569, 340.713] (0.95 confidence) | 342.462 ns [341.59, 343.352] (0.95 confidence) | 0.993212 |
| rational&#95;dot/large/d2/slice | 122.96 ns [122.611, 123.332] (0.95 confidence) | 125.324 ns [124.997, 125.663] (0.95 confidence) | 0.981133 |
| rational&#95;dot/large/d22/prepared | 4731.83 ns [4721.34, 4742.62] (0.95 confidence) | 4789.87 ns [4777.95, 4802.39] (0.95 confidence) | 0.987884 |
| rational&#95;dot/large/d22/slice | 2108.37 ns [2100.53, 2116.22] (0.95 confidence) | 2121.39 ns [2115.46, 2127.39] (0.95 confidence) | 0.993859 |
| rational&#95;dot/large/d3/prepared | 555.13 ns [552.984, 557.379] (0.95 confidence) | 553.917 ns [552.41, 555.485] (0.95 confidence) | 1.00219 |
| rational&#95;dot/large/d3/slice | 215.402 ns [214.948, 215.87] (0.95 confidence) | 218.444 ns [217.555, 219.456] (0.95 confidence) | 0.986074 |
| rational&#95;dot/large/d4/prepared | 741.624 ns [739.698, 743.662] (0.95 confidence) | 738.906 ns [737.51, 740.347] (0.95 confidence) | 1.00368 |
| rational&#95;dot/large/d4/slice | 294.881 ns [293.58, 296.23] (0.95 confidence) | 307.91 ns [306.206, 309.627] (0.95 confidence) | 0.957686 |
| rational&#95;dot/large/d5/prepared | 955.17 ns [952.033, 958.542] (0.95 confidence) | 952.624 ns [950.552, 954.761] (0.95 confidence) | 1.00267 |
| rational&#95;dot/large/d5/slice | 400.465 ns [398.763, 402.208] (0.95 confidence) | 393.547 ns [391.898, 395.3] (0.95 confidence) | 1.01758 |
| rational&#95;dot/large/d6/prepared | 1179.37 ns [1176.66, 1182.1] (0.95 confidence) | 1162.48 ns [1160.15, 1164.89] (0.95 confidence) | 1.01453 |
| rational&#95;dot/large/d6/slice | 492.188 ns [490.568, 493.835] (0.95 confidence) | 495.539 ns [493.785, 497.348] (0.95 confidence) | 0.993238 |
| rational&#95;dot/large/d7/prepared | 1402.55 ns [1398.04, 1407.25] (0.95 confidence) | 1399.43 ns [1396.7, 1402.25] (0.95 confidence) | 1.00223 |
| rational&#95;dot/large/d7/slice | 599.765 ns [597.578, 602.011] (0.95 confidence) | 596.971 ns [594.672, 599.422] (0.95 confidence) | 1.00468 |
| rational&#95;dot/large/d8/prepared | 1583.39 ns [1580.9, 1585.97] (0.95 confidence) | 1604.22 ns [1600.18, 1608.46] (0.95 confidence) | 0.987021 |
| rational&#95;dot/large/d8/slice | 673.788 ns [671.405, 676.282] (0.95 confidence) | 704.055 ns [699.972, 708.372] (0.95 confidence) | 0.95701 |
| rational&#95;dot/non&#95;dyadic/d14/prepared | 2409.34 ns [2405.64, 2412.9] (0.95 confidence) | 2445.42 ns [2429.35, 2462.65] (0.95 confidence) | 0.985246 |
| rational&#95;dot/non&#95;dyadic/d14/slice | 1052.52 ns [1048.81, 1056.26] (0.95 confidence) | 1018.62 ns [1016.31, 1020.85] (0.95 confidence) | 1.03328 |
| rational&#95;dot/non&#95;dyadic/d2/prepared | 275.149 ns [274.523, 275.835] (0.95 confidence) | 277.218 ns [276.592, 277.882] (0.95 confidence) | 0.992537 |
| rational&#95;dot/non&#95;dyadic/d2/slice | 92.4407 ns [92.2982, 92.5853] (0.95 confidence) | 95.5251 ns [95.2083, 95.8544] (0.95 confidence) | 0.967711 |
| rational&#95;dot/non&#95;dyadic/d22/prepared | 3876.23 ns [3868.88, 3883.81] (0.95 confidence) | 3875.37 ns [3866.95, 3884] (0.95 confidence) | 1.00022 |
| rational&#95;dot/non&#95;dyadic/d22/slice | 1710.65 ns [1705.95, 1715.49] (0.95 confidence) | 1712.38 ns [1705.93, 1718.84] (0.95 confidence) | 0.998991 |
| rational&#95;dot/non&#95;dyadic/d3/prepared | 447.392 ns [446.876, 447.903] (0.95 confidence) | 448.44 ns [447.752, 449.145] (0.95 confidence) | 0.997663 |
| rational&#95;dot/non&#95;dyadic/d3/slice | 166.987 ns [166.72, 167.256] (0.95 confidence) | 170.725 ns [170.353, 171.115] (0.95 confidence) | 0.978101 |
| rational&#95;dot/non&#95;dyadic/d4/prepared | 599.132 ns [598.202, 600.141] (0.95 confidence) | 599.973 ns [598.349, 601.92] (0.95 confidence) | 0.998599 |
| rational&#95;dot/non&#95;dyadic/d4/slice | 228.693 ns [228.142, 229.287] (0.95 confidence) | 230.166 ns [229.624, 230.747] (0.95 confidence) | 0.993604 |
| rational&#95;dot/non&#95;dyadic/d5/prepared | 770.855 ns [769.358, 772.551] (0.95 confidence) | 775.66 ns [773.317, 778.219] (0.95 confidence) | 0.993804 |
| rational&#95;dot/non&#95;dyadic/d5/slice | 303.844 ns [302.894, 304.83] (0.95 confidence) | 310.141 ns [309.305, 311.067] (0.95 confidence) | 0.979696 |
| rational&#95;dot/non&#95;dyadic/d6/prepared | 953.401 ns [951.441, 955.352] (0.95 confidence) | 954.849 ns [952.678, 957.121] (0.95 confidence) | 0.998483 |
| rational&#95;dot/non&#95;dyadic/d6/slice | 386.631 ns [385.703, 387.569] (0.95 confidence) | 391.032 ns [389.85, 392.217] (0.95 confidence) | 0.988744 |
| rational&#95;dot/non&#95;dyadic/d7/prepared | 1142.9 ns [1140.69, 1145.2] (0.95 confidence) | 1134.91 ns [1133.66, 1136.24] (0.95 confidence) | 1.00704 |
| rational&#95;dot/non&#95;dyadic/d7/slice | 482.737 ns [481.697, 483.816] (0.95 confidence) | 487.846 ns [486.53, 489.172] (0.95 confidence) | 0.989527 |
| rational&#95;dot/non&#95;dyadic/d8/prepared | 1287.98 ns [1285.35, 1290.77] (0.95 confidence) | 1284.49 ns [1282.76, 1286.35] (0.95 confidence) | 1.00271 |
| rational&#95;dot/non&#95;dyadic/d8/slice | 541.923 ns [540.476, 543.424] (0.95 confidence) | 539.175 ns [537.684, 540.724] (0.95 confidence) | 1.0051 |

### Median (ns)

| Workload | Baseline | Current | Baseline / current |
| --- | --- | --- | --- |
| rational&#95;dot/cancellation/d14/prepared | 1863.31 ns [1861.34, 1867.04] (0.95 confidence) | 1831.53 ns [1829.21, 1838.65] (0.95 confidence) | 1.01735 |
| rational&#95;dot/cancellation/d14/slice | 645.235 ns [643.925, 648.081] (0.95 confidence) | 646.492 ns [643.689, 649.726] (0.95 confidence) | 0.998055 |
| rational&#95;dot/cancellation/d2/prepared | 268.06 ns [267.766, 268.537] (0.95 confidence) | 272.3 ns [271.656, 273.136] (0.95 confidence) | 0.984429 |
| rational&#95;dot/cancellation/d2/slice | 96.7527 ns [96.6464, 96.8219] (0.95 confidence) | 97.9773 ns [97.838, 98.1229] (0.95 confidence) | 0.987501 |
| rational&#95;dot/cancellation/d22/prepared | 2920.9 ns [2913.96, 2930.66] (0.95 confidence) | 2879.48 ns [2855.17, 2894.49] (0.95 confidence) | 1.01438 |
| rational&#95;dot/cancellation/d22/slice | 999.613 ns [998.225, 1001.16] (0.95 confidence) | 994.363 ns [992.339, 995.739] (0.95 confidence) | 1.00528 |
| rational&#95;dot/cancellation/d3/prepared | 381.793 ns [381.043, 383.123] (0.95 confidence) | 386.851 ns [386.352, 387.435] (0.95 confidence) | 0.986925 |
| rational&#95;dot/cancellation/d3/slice | 126.668 ns [126.419, 126.869] (0.95 confidence) | 128.04 ns [127.846, 128.189] (0.95 confidence) | 0.989286 |
| rational&#95;dot/cancellation/d4/prepared | 531.89 ns [530.686, 533.436] (0.95 confidence) | 537.001 ns [535.981, 537.801] (0.95 confidence) | 0.990483 |
| rational&#95;dot/cancellation/d4/slice | 187.73 ns [187.26, 188.231] (0.95 confidence) | 188.675 ns [188.256, 189.188] (0.95 confidence) | 0.994991 |
| rational&#95;dot/cancellation/d5/prepared | 643.333 ns [642.443, 644.004] (0.95 confidence) | 654.146 ns [653.224, 655.748] (0.95 confidence) | 0.983469 |
| rational&#95;dot/cancellation/d5/slice | 213.387 ns [213.039, 213.748] (0.95 confidence) | 216.98 ns [216.218, 217.299] (0.95 confidence) | 0.983438 |
| rational&#95;dot/cancellation/d6/prepared | 806.724 ns [805.128, 809.388] (0.95 confidence) | 809.215 ns [806.603, 813.061] (0.95 confidence) | 0.996921 |
| rational&#95;dot/cancellation/d6/slice | 275.329 ns [274.667, 276.229] (0.95 confidence) | 276.594 ns [275.93, 277.517] (0.95 confidence) | 0.99543 |
| rational&#95;dot/cancellation/d7/prepared | 917.142 ns [915.154, 919.729] (0.95 confidence) | 920.915 ns [919.242, 922.838] (0.95 confidence) | 0.995903 |
| rational&#95;dot/cancellation/d7/slice | 306.673 ns [306.142, 307.286] (0.95 confidence) | 307.543 ns [307.066, 309.02] (0.95 confidence) | 0.997172 |
| rational&#95;dot/cancellation/d8/prepared | 1065.38 ns [1062.95, 1068.32] (0.95 confidence) | 1066.85 ns [1062.28, 1070.31] (0.95 confidence) | 0.998622 |
| rational&#95;dot/cancellation/d8/slice | 364.03 ns [363.236, 365.036] (0.95 confidence) | 372.931 ns [371.408, 374.094] (0.95 confidence) | 0.976133 |
| rational&#95;dot/large/d14/prepared | 2968 ns [2958.99, 2981.19] (0.95 confidence) | 3007.52 ns [2985.66, 3015.82] (0.95 confidence) | 0.98686 |
| rational&#95;dot/large/d14/slice | 1293.59 ns [1288.24, 1299.38] (0.95 confidence) | 1298.21 ns [1294.35, 1302.18] (0.95 confidence) | 0.996442 |
| rational&#95;dot/large/d2/prepared | 340.074 ns [339.503, 340.739] (0.95 confidence) | 341.661 ns [340.785, 343.07] (0.95 confidence) | 0.995355 |
| rational&#95;dot/large/d2/slice | 122.52 ns [122.247, 122.841] (0.95 confidence) | 124.825 ns [124.628, 125.431] (0.95 confidence) | 0.981528 |
| rational&#95;dot/large/d22/prepared | 4735.56 ns [4723.37, 4745.63] (0.95 confidence) | 4790.83 ns [4774.52, 4802.77] (0.95 confidence) | 0.988464 |
| rational&#95;dot/large/d22/slice | 2110.7 ns [2094.64, 2120.18] (0.95 confidence) | 2121.79 ns [2110.61, 2129.26] (0.95 confidence) | 0.994773 |
| rational&#95;dot/large/d3/prepared | 553.271 ns [550.845, 555.477] (0.95 confidence) | 552.083 ns [550.496, 553.294] (0.95 confidence) | 1.00215 |
| rational&#95;dot/large/d3/slice | 215.4 ns [214.833, 215.808] (0.95 confidence) | 217.334 ns [216.464, 218.167] (0.95 confidence) | 0.991104 |
| rational&#95;dot/large/d4/prepared | 739.528 ns [738.069, 741.521] (0.95 confidence) | 737.913 ns [736.317, 739.384] (0.95 confidence) | 1.00219 |
| rational&#95;dot/large/d4/slice | 293.31 ns [291.376, 296.198] (0.95 confidence) | 307.603 ns [305.797, 309.515] (0.95 confidence) | 0.953537 |
| rational&#95;dot/large/d5/prepared | 951.234 ns [948.856, 955.112] (0.95 confidence) | 951.186 ns [947.923, 953.605] (0.95 confidence) | 1.00005 |
| rational&#95;dot/large/d5/slice | 400.261 ns [398.017, 401.907] (0.95 confidence) | 392.632 ns [391.283, 393.812] (0.95 confidence) | 1.01943 |
| rational&#95;dot/large/d6/prepared | 1179.79 ns [1175.53, 1183.18] (0.95 confidence) | 1161.57 ns [1158.82, 1165.65] (0.95 confidence) | 1.01569 |
| rational&#95;dot/large/d6/slice | 492.903 ns [488.072, 495.344] (0.95 confidence) | 494.569 ns [493.086, 495.824] (0.95 confidence) | 0.996631 |
| rational&#95;dot/large/d7/prepared | 1396.51 ns [1393.73, 1402.68] (0.95 confidence) | 1398.77 ns [1393.75, 1402.39] (0.95 confidence) | 0.99839 |
| rational&#95;dot/large/d7/slice | 599.702 ns [596.652, 601.995] (0.95 confidence) | 594.83 ns [592.982, 596.2] (0.95 confidence) | 1.00819 |
| rational&#95;dot/large/d8/prepared | 1583.75 ns [1580.26, 1585.17] (0.95 confidence) | 1599.78 ns [1593.73, 1603.97] (0.95 confidence) | 0.989979 |
| rational&#95;dot/large/d8/slice | 671.461 ns [669.33, 673.538] (0.95 confidence) | 701.41 ns [697.449, 704.864] (0.95 confidence) | 0.957302 |
| rational&#95;dot/non&#95;dyadic/d14/prepared | 2408.21 ns [2405.8, 2412.1] (0.95 confidence) | 2409.29 ns [2402.8, 2416.44] (0.95 confidence) | 0.999553 |
| rational&#95;dot/non&#95;dyadic/d14/slice | 1052.27 ns [1046.99, 1056.99] (0.95 confidence) | 1018.53 ns [1017.01, 1020.67] (0.95 confidence) | 1.03312 |
| rational&#95;dot/non&#95;dyadic/d2/prepared | 274.766 ns [274.234, 275.391] (0.95 confidence) | 276.507 ns [275.845, 277.024] (0.95 confidence) | 0.993705 |
| rational&#95;dot/non&#95;dyadic/d2/slice | 92.3515 ns [92.2079, 92.5161] (0.95 confidence) | 95.1917 ns [94.9351, 95.5183] (0.95 confidence) | 0.970164 |
| rational&#95;dot/non&#95;dyadic/d22/prepared | 3874.05 ns [3866.66, 3880.36] (0.95 confidence) | 3876.06 ns [3866.25, 3882.45] (0.95 confidence) | 0.999482 |
| rational&#95;dot/non&#95;dyadic/d22/slice | 1710.64 ns [1699.05, 1716.34] (0.95 confidence) | 1714.49 ns [1706.07, 1723.6] (0.95 confidence) | 0.997755 |
| rational&#95;dot/non&#95;dyadic/d3/prepared | 447.567 ns [447, 448.077] (0.95 confidence) | 447.906 ns [446.856, 448.716] (0.95 confidence) | 0.999244 |
| rational&#95;dot/non&#95;dyadic/d3/slice | 167.1 ns [166.758, 167.413] (0.95 confidence) | 170.574 ns [170.031, 170.994] (0.95 confidence) | 0.979632 |
| rational&#95;dot/non&#95;dyadic/d4/prepared | 598.763 ns [597.664, 599.504] (0.95 confidence) | 598.143 ns [597.102, 599.341] (0.95 confidence) | 1.00104 |
| rational&#95;dot/non&#95;dyadic/d4/slice | 228.402 ns [227.636, 229.014] (0.95 confidence) | 229.525 ns [228.931, 230.518] (0.95 confidence) | 0.995103 |
| rational&#95;dot/non&#95;dyadic/d5/prepared | 769.704 ns [768.358, 771.423] (0.95 confidence) | 772.173 ns [770.698, 774.001] (0.95 confidence) | 0.996803 |
| rational&#95;dot/non&#95;dyadic/d5/slice | 303.446 ns [302.104, 304.453] (0.95 confidence) | 309.342 ns [308.521, 310.599] (0.95 confidence) | 0.980942 |
| rational&#95;dot/non&#95;dyadic/d6/prepared | 953.662 ns [951.977, 955.058] (0.95 confidence) | 953.694 ns [950.94, 955.292] (0.95 confidence) | 0.999967 |
| rational&#95;dot/non&#95;dyadic/d6/slice | 386.219 ns [385.321, 387.793] (0.95 confidence) | 391.056 ns [389.072, 392.764] (0.95 confidence) | 0.987629 |
| rational&#95;dot/non&#95;dyadic/d7/prepared | 1140.5 ns [1138.39, 1143.8] (0.95 confidence) | 1133.75 ns [1132.87, 1135.38] (0.95 confidence) | 1.00595 |
| rational&#95;dot/non&#95;dyadic/d7/slice | 482.115 ns [480.459, 483.209] (0.95 confidence) | 487.358 ns [485.942, 490.479] (0.95 confidence) | 0.989243 |
| rational&#95;dot/non&#95;dyadic/d8/prepared | 1285.02 ns [1283.14, 1287.18] (0.95 confidence) | 1283.58 ns [1281.6, 1285.32] (0.95 confidence) | 1.00113 |
| rational&#95;dot/non&#95;dyadic/d8/slice | 540.345 ns [539.26, 542.131] (0.95 confidence) | 537.062 ns [535.613, 538.614] (0.95 confidence) | 1.00611 |
