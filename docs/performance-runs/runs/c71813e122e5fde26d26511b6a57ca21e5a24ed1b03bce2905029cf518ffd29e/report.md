<!-- research-repo-tools/run-report/v3 -->
# Benchmark timings

Run: `c71813e122e5fde26d26511b6a57ca21e5a24ed1b03bce2905029cf518ffd29e`

Evidence: [run payload](run.json) · [evidence envelope](evidence.json).

Acquisition schedule: sequential; trials: 1; seed: None.
Each trial is retained separately. Within-phase samples are not independent trials; no pooled significance claim is made.

## Trial 1: baseline → current

Collection: 1 segment(s); 0 abandoned attempt(s). Validated abandoned phases are retained and excluded from accepted pairs; failed-output capture status is
recorded.
These are completed campaign slots; recovery does not establish uninterrupted collection, independent trials, or absence of optional-stopping bias.

- Collection segment `081b4a0390324a1794091ed46711f042` started: 1794894043471548000 (Unix ns).
- core-default/baseline: 1794894450174032000–1794900434796846000 (Unix ns).
- core-default/current: 1794900453469856000–1794906468236730000 (Unix ns).
- core-exact/baseline: 1794906689492048000–1794912395331095000 (Unix ns).
- core-exact/current: 1794912416049951000–1794918193584729000 (Unix ns).
- lu-solve/baseline: 1794918408505207000–1794918637085957000 (Unix ns).
- lu-solve/current: 1794918640134901000–1794918867105624000 (Unix ns).
- column-basis/baseline: 1794918953409635000–1794919464918980000 (Unix ns).
- column-basis/current: 1794919468206735000–1794919983671173000 (Unix ns).
- interval-solve/baseline: 1794920077740109000–1794920414287723000 (Unix ns).
- interval-solve/current: 1794920417938338000–1794920752532128000 (Unix ns).
- polynomial/baseline: 1794920839246935000–1794920894724822000 (Unix ns).
- polynomial/current: 1794920897926951000–1794920954437757000 (Unix ns).
- rational-dot/baseline: 1794921046815527000–1794921551633299000 (Unix ns).
- rational-dot/current: 1794921555080963000–1794922063151143000 (Unix ns).

Ratios are point estimates within one accepted paired suite. Marginal timing intervals do not establish significance.

## core-default

Semantic contract: la-stack-core-default-v1. Collection session: `f74705328ddc415a93bd01622df222a0`.
Collection interval (Unix ns): 1794894043476571000–1794906508361461000.
Acquisition: newly collected.

- baseline: available — inside declared release interval; release v0.5.3; revision `7165ed562f2e4611d201e574fc1f600a8ad8b33e`.
  Source SHA-256: `073502f95ae58344aa776b54303ae24d1e6cf2ce60d1437d950bcbe64d0b4dbf`; harness SHA-256:
  `ba63ffff14e6124d5aca8df516e6836a606ad10809bbe3b1ab13dd5fd7ceedbf`; scientific gate: passed.
- current: available — inside declared release interval; release v0.5.4; revision `531f75eb58eb8ac810ad5d89bf1d3e66a301e84e`.
  Source SHA-256: `50d72eb42299ad56e00e6838b34b10824c560a3eb5ae9c9a05e48859fd3cfe46`; harness SHA-256:
  `ba63ffff14e6124d5aca8df516e6836a606ad10809bbe3b1ab13dd5fd7ceedbf`; scientific gate: passed.

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
<td>23.6191 ns [23.5916, 23.6456] (0.95 confidence)</td>
<td>24.6319 ns [24.5833, 24.6825] (0.95 confidence)</td>
<td>0.958884</td>
</tr>
<tr>
<td>angle&#95;d3/dense/construct&#95;then&#95;angle</td>
<td>22.6903 ns [22.6539, 22.7218] (0.95 confidence)</td>
<td>23.6939 ns [23.6292, 23.7621] (0.95 confidence)</td>
<td>0.957644</td>
</tr>
<tr>
<td>angle&#95;d3/dense/prepared&#95;vector</td>
<td>22.3849 ns [22.3485, 22.4173] (0.95 confidence)</td>
<td>23.6495 ns [23.5743, 23.7294] (0.95 confidence)</td>
<td>0.946529</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;antipodal/borrowed&#95;slice</td>
<td>16.6957 ns [16.6724, 16.725] (0.95 confidence)</td>
<td>22.7954 ns [21.9715, 23.7682] (0.95 confidence)</td>
<td>0.732416</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;antipodal/construct&#95;then&#95;angle</td>
<td>15.5943 ns [15.5779, 15.6111] (0.95 confidence)</td>
<td>20.2723 ns [19.7169, 20.8655] (0.95 confidence)</td>
<td>0.769243</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;antipodal/prepared&#95;vector</td>
<td>15.6148 ns [15.5942, 15.6342] (0.95 confidence)</td>
<td>21.0544 ns [20.4696, 21.6822] (0.95 confidence)</td>
<td>0.741641</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;identical/borrowed&#95;slice</td>
<td>14.1703 ns [14.1474, 14.1901] (0.95 confidence)</td>
<td>17.9503 ns [17.5073, 18.4045] (0.95 confidence)</td>
<td>0.789419</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;identical/construct&#95;then&#95;angle</td>
<td>13.4634 ns [13.4553, 13.4714] (0.95 confidence)</td>
<td>18.5847 ns [18.0411, 19.182] (0.95 confidence)</td>
<td>0.724438</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;identical/prepared&#95;vector</td>
<td>13.333 ns [13.3146, 13.3519] (0.95 confidence)</td>
<td>14.6218 ns [14.3445, 14.9256] (0.95 confidence)</td>
<td>0.911856</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;varied/borrowed&#95;slice</td>
<td>25.7516 ns [25.7165, 25.785] (0.95 confidence)</td>
<td>26.6952 ns [26.689, 26.7027] (0.95 confidence)</td>
<td>0.964652</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;varied/construct&#95;then&#95;angle</td>
<td>24.8007 ns [24.773, 24.8267] (0.95 confidence)</td>
<td>25.8022 ns [25.7371, 25.8823] (0.95 confidence)</td>
<td>0.961185</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;varied/prepared&#95;vector</td>
<td>24.4338 ns [24.3829, 24.4964] (0.95 confidence)</td>
<td>25.3279 ns [25.316, 25.3414] (0.95 confidence)</td>
<td>0.964699</td>
</tr>
<tr>
<td>angle&#95;d3/maximum&#95;small&#95;angle/borrowed&#95;slice</td>
<td>11.626 ns [11.6179, 11.6338] (0.95 confidence)</td>
<td>11.6651 ns [11.6473, 11.6843] (0.95 confidence)</td>
<td>0.99665</td>
</tr>
<tr>
<td>angle&#95;d3/maximum&#95;small&#95;angle/construct&#95;then&#95;angle</td>
<td>10.4313 ns [10.4264, 10.4364] (0.95 confidence)</td>
<td>10.3577 ns [10.3494, 10.3667] (0.95 confidence)</td>
<td>1.0071</td>
</tr>
<tr>
<td>angle&#95;d3/maximum&#95;small&#95;angle/prepared&#95;vector</td>
<td>9.74664 ns [9.73776, 9.75566] (0.95 confidence)</td>
<td>9.74635 ns [9.73676, 9.7559] (0.95 confidence)</td>
<td>1.00003</td>
</tr>
<tr>
<td>angle&#95;d3/mixed&#95;scale/borrowed&#95;slice</td>
<td>12.9724 ns [12.9612, 12.9846] (0.95 confidence)</td>
<td>13.0496 ns [13.0392, 13.0596] (0.95 confidence)</td>
<td>0.994084</td>
</tr>
<tr>
<td>angle&#95;d3/mixed&#95;scale/construct&#95;then&#95;angle</td>
<td>11.6401 ns [11.6265, 11.654] (0.95 confidence)</td>
<td>11.5405 ns [11.5221, 11.5602] (0.95 confidence)</td>
<td>1.00863</td>
</tr>
<tr>
<td>angle&#95;d3/mixed&#95;scale/prepared&#95;vector</td>
<td>11.042 ns [11.0313, 11.0528] (0.95 confidence)</td>
<td>11.0217 ns [11.0096, 11.0335] (0.95 confidence)</td>
<td>1.00185</td>
</tr>
<tr>
<td>angle&#95;d3/near&#95;antipodal/borrowed&#95;slice</td>
<td>17.2048 ns [17.1864, 17.2227] (0.95 confidence)</td>
<td>17.0954 ns [17.0729, 17.1167] (0.95 confidence)</td>
<td>1.0064</td>
</tr>
<tr>
<td>angle&#95;d3/near&#95;antipodal/construct&#95;then&#95;angle</td>
<td>15.8046 ns [15.7932, 15.8142] (0.95 confidence)</td>
<td>15.5863 ns [15.5408, 15.6231] (0.95 confidence)</td>
<td>1.01401</td>
</tr>
<tr>
<td>angle&#95;d3/near&#95;antipodal/prepared&#95;vector</td>
<td>15.3133 ns [15.2817, 15.3474] (0.95 confidence)</td>
<td>15.281 ns [15.2338, 15.3288] (0.95 confidence)</td>
<td>1.00211</td>
</tr>
<tr>
<td>angle&#95;d3/near&#95;parallel/borrowed&#95;slice</td>
<td>11.6619 ns [11.6531, 11.6714] (0.95 confidence)</td>
<td>11.6352 ns [11.6249, 11.6456] (0.95 confidence)</td>
<td>1.00229</td>
</tr>
<tr>
<td>angle&#95;d3/near&#95;parallel/construct&#95;then&#95;angle</td>
<td>10.433 ns [10.4273, 10.4388] (0.95 confidence)</td>
<td>10.3542 ns [10.3466, 10.3609] (0.95 confidence)</td>
<td>1.0076</td>
</tr>
<tr>
<td>angle&#95;d3/near&#95;parallel/prepared&#95;vector</td>
<td>9.74128 ns [9.73193, 9.75071] (0.95 confidence)</td>
<td>9.73968 ns [9.7295, 9.74993] (0.95 confidence)</td>
<td>1.00016</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;adjacent/borrowed&#95;slice</td>
<td>14.007 ns [13.9815, 14.033] (0.95 confidence)</td>
<td>15.1981 ns [15.0992, 15.3015] (0.95 confidence)</td>
<td>0.921632</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;adjacent/construct&#95;then&#95;angle</td>
<td>13.0201 ns [12.9951, 13.0457] (0.95 confidence)</td>
<td>14.6887 ns [14.3208, 15.0987] (0.95 confidence)</td>
<td>0.886402</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;adjacent/prepared&#95;vector</td>
<td>12.0552 ns [12.0396, 12.0707] (0.95 confidence)</td>
<td>15.8621 ns [15.3854, 16.3948] (0.95 confidence)</td>
<td>0.76</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;alternating/borrowed&#95;slice</td>
<td>14.2034 ns [14.174, 14.2347] (0.95 confidence)</td>
<td>15.3187 ns [15.2008, 15.436] (0.95 confidence)</td>
<td>0.927191</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;alternating/construct&#95;then&#95;angle</td>
<td>13.3709 ns [13.3134, 13.4373] (0.95 confidence)</td>
<td>14.5467 ns [14.2299, 14.8839] (0.95 confidence)</td>
<td>0.919174</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;alternating/prepared&#95;vector</td>
<td>12.3497 ns [12.3329, 12.366] (0.95 confidence)</td>
<td>16.5265 ns [16.0362, 17.0892] (0.95 confidence)</td>
<td>0.747267</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;antipodal/borrowed&#95;slice</td>
<td>10.9234 ns [10.8962, 10.9532] (0.95 confidence)</td>
<td>10.936 ns [10.9124, 10.9631] (0.95 confidence)</td>
<td>0.998855</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;antipodal/construct&#95;then&#95;angle</td>
<td>10.1776 ns [10.1511, 10.2056] (0.95 confidence)</td>
<td>10.2307 ns [10.203, 10.2603] (0.95 confidence)</td>
<td>0.994807</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;antipodal/prepared&#95;vector</td>
<td>8.88674 ns [8.87749, 8.89531] (0.95 confidence)</td>
<td>8.91164 ns [8.90842, 8.91555] (0.95 confidence)</td>
<td>0.997206</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;identical/borrowed&#95;slice</td>
<td>8.69509 ns [8.67845, 8.71412] (0.95 confidence)</td>
<td>8.9077 ns [8.86236, 8.95744] (0.95 confidence)</td>
<td>0.976132</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;identical/construct&#95;then&#95;angle</td>
<td>8.08129 ns [8.05925, 8.10347] (0.95 confidence)</td>
<td>8.70374 ns [8.55169, 8.871] (0.95 confidence)</td>
<td>0.928485</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;identical/prepared&#95;vector</td>
<td>6.94585 ns [6.9416, 6.94987] (0.95 confidence)</td>
<td>8.9891 ns [8.71188, 9.3126] (0.95 confidence)</td>
<td>0.772698</td>
</tr>
<tr>
<td>angle&#95;d3/subnormal/borrowed&#95;slice</td>
<td>11.661 ns [11.6487, 11.6745] (0.95 confidence)</td>
<td>11.6663 ns [11.6539, 11.6804] (0.95 confidence)</td>
<td>0.999546</td>
</tr>
<tr>
<td>angle&#95;d3/subnormal/construct&#95;then&#95;angle</td>
<td>10.4388 ns [10.4339, 10.4439] (0.95 confidence)</td>
<td>10.3604 ns [10.3525, 10.3685] (0.95 confidence)</td>
<td>1.00757</td>
</tr>
<tr>
<td>angle&#95;d3/subnormal/prepared&#95;vector</td>
<td>9.74744 ns [9.73913, 9.75576] (0.95 confidence)</td>
<td>9.75183 ns [9.73831, 9.76768] (0.95 confidence)</td>
<td>0.999549</td>
</tr>
<tr>
<td>angle&#95;d3/subnormal&#95;axes/borrowed&#95;slice</td>
<td>14.1955 ns [14.1789, 14.2129] (0.95 confidence)</td>
<td>14.2706 ns [14.2404, 14.3032] (0.95 confidence)</td>
<td>0.994736</td>
</tr>
<tr>
<td>angle&#95;d3/subnormal&#95;axes/construct&#95;then&#95;angle</td>
<td>13.1957 ns [13.1813, 13.2113] (0.95 confidence)</td>
<td>13.2395 ns [13.2124, 13.2677] (0.95 confidence)</td>
<td>0.996689</td>
</tr>
<tr>
<td>angle&#95;d3/subnormal&#95;axes/prepared&#95;vector</td>
<td>12.3577 ns [12.3353, 12.3803] (0.95 confidence)</td>
<td>12.4119 ns [12.3856, 12.4415] (0.95 confidence)</td>
<td>0.995637</td>
</tr>
<tr>
<td>angle&#95;d4/dense/borrowed&#95;slice</td>
<td>33.5269 ns [33.4843, 33.5659] (0.95 confidence)</td>
<td>33.4374 ns [33.3878, 33.4805] (0.95 confidence)</td>
<td>1.00268</td>
</tr>
<tr>
<td>angle&#95;d4/dense/construct&#95;then&#95;angle</td>
<td>32.366 ns [32.3283, 32.4017] (0.95 confidence)</td>
<td>32.1909 ns [32.1495, 32.2294] (0.95 confidence)</td>
<td>1.00544</td>
</tr>
<tr>
<td>angle&#95;d4/dense/prepared&#95;vector</td>
<td>31.9905 ns [31.9252, 32.0725] (0.95 confidence)</td>
<td>32.0635 ns [32.0167, 32.1117] (0.95 confidence)</td>
<td>0.997724</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;antipodal/borrowed&#95;slice</td>
<td>26.2405 ns [26.2137, 26.2657] (0.95 confidence)</td>
<td>26.1523 ns [26.1086, 26.2002] (0.95 confidence)</td>
<td>1.00337</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;antipodal/construct&#95;then&#95;angle</td>
<td>24.8397 ns [24.8117, 24.8666] (0.95 confidence)</td>
<td>24.7438 ns [24.7158, 24.7704] (0.95 confidence)</td>
<td>1.00388</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;antipodal/prepared&#95;vector</td>
<td>24.4057 ns [24.3684, 24.439] (0.95 confidence)</td>
<td>24.454 ns [24.4274, 24.4783] (0.95 confidence)</td>
<td>0.998023</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;identical/borrowed&#95;slice</td>
<td>23.9057 ns [23.8816, 23.9282] (0.95 confidence)</td>
<td>23.7824 ns [23.7577, 23.8065] (0.95 confidence)</td>
<td>1.00519</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;identical/construct&#95;then&#95;angle</td>
<td>22.3986 ns [22.375, 22.4216] (0.95 confidence)</td>
<td>22.468 ns [22.4347, 22.5042] (0.95 confidence)</td>
<td>0.996909</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;identical/prepared&#95;vector</td>
<td>21.9048 ns [21.8581, 21.9565] (0.95 confidence)</td>
<td>21.8752 ns [21.8389, 21.9125] (0.95 confidence)</td>
<td>1.00135</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;varied/borrowed&#95;slice</td>
<td>34.9173 ns [34.8697, 34.9628] (0.95 confidence)</td>
<td>34.9903 ns [34.9433, 35.0362] (0.95 confidence)</td>
<td>0.997913</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;varied/construct&#95;then&#95;angle</td>
<td>33.7019 ns [33.6437, 33.7588] (0.95 confidence)</td>
<td>33.7705 ns [33.7277, 33.8121] (0.95 confidence)</td>
<td>0.997966</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;varied/prepared&#95;vector</td>
<td>32.7854 ns [32.7065, 32.8538] (0.95 confidence)</td>
<td>32.7618 ns [32.7086, 32.8119] (0.95 confidence)</td>
<td>1.00072</td>
</tr>
<tr>
<td>angle&#95;d4/maximum&#95;small&#95;angle/borrowed&#95;slice</td>
<td>13.029 ns [12.9945, 13.0671] (0.95 confidence)</td>
<td>12.7511 ns [12.7084, 12.7982] (0.95 confidence)</td>
<td>1.02179</td>
</tr>
<tr>
<td>angle&#95;d4/maximum&#95;small&#95;angle/construct&#95;then&#95;angle</td>
<td>11.6179 ns [11.5787, 11.6583] (0.95 confidence)</td>
<td>11.1865 ns [11.1662, 11.2079] (0.95 confidence)</td>
<td>1.03856</td>
</tr>
<tr>
<td>angle&#95;d4/maximum&#95;small&#95;angle/prepared&#95;vector</td>
<td>10.5574 ns [10.5473, 10.5677] (0.95 confidence)</td>
<td>10.3105 ns [10.2958, 10.3268] (0.95 confidence)</td>
<td>1.02395</td>
</tr>
<tr>
<td>angle&#95;d4/mixed&#95;scale/borrowed&#95;slice</td>
<td>14.1705 ns [14.1428, 14.2027] (0.95 confidence)</td>
<td>14.2029 ns [14.1652, 14.2465] (0.95 confidence)</td>
<td>0.997724</td>
</tr>
<tr>
<td>angle&#95;d4/mixed&#95;scale/construct&#95;then&#95;angle</td>
<td>12.73 ns [12.6763, 12.786] (0.95 confidence)</td>
<td>12.5561 ns [12.5324, 12.5811] (0.95 confidence)</td>
<td>1.01386</td>
</tr>
<tr>
<td>angle&#95;d4/mixed&#95;scale/prepared&#95;vector</td>
<td>11.7692 ns [11.7506, 11.7879] (0.95 confidence)</td>
<td>11.7517 ns [11.7262, 11.7798] (0.95 confidence)</td>
<td>1.00148</td>
</tr>
<tr>
<td>angle&#95;d4/near&#95;antipodal/borrowed&#95;slice</td>
<td>18.1458 ns [18.1388, 18.1534] (0.95 confidence)</td>
<td>18.1012 ns [18.0796, 18.1205] (0.95 confidence)</td>
<td>1.00246</td>
</tr>
<tr>
<td>angle&#95;d4/near&#95;antipodal/construct&#95;then&#95;angle</td>
<td>16.5569 ns [16.5474, 16.5664] (0.95 confidence)</td>
<td>16.5156 ns [16.5018, 16.5293] (0.95 confidence)</td>
<td>1.0025</td>
</tr>
<tr>
<td>angle&#95;d4/near&#95;antipodal/prepared&#95;vector</td>
<td>16.0653 ns [16.0554, 16.0756] (0.95 confidence)</td>
<td>15.9976 ns [15.9839, 16.0106] (0.95 confidence)</td>
<td>1.00423</td>
</tr>
<tr>
<td>angle&#95;d4/near&#95;parallel/borrowed&#95;slice</td>
<td>12.9383 ns [12.8969, 12.9807] (0.95 confidence)</td>
<td>12.6573 ns [12.6288, 12.6877] (0.95 confidence)</td>
<td>1.0222</td>
</tr>
<tr>
<td>angle&#95;d4/near&#95;parallel/construct&#95;then&#95;angle</td>
<td>11.2214 ns [11.1936, 11.2506] (0.95 confidence)</td>
<td>11.1915 ns [11.1607, 11.2252] (0.95 confidence)</td>
<td>1.00267</td>
</tr>
<tr>
<td>angle&#95;d4/near&#95;parallel/prepared&#95;vector</td>
<td>10.5603 ns [10.5498, 10.5712] (0.95 confidence)</td>
<td>10.3076 ns [10.2984, 10.3171] (0.95 confidence)</td>
<td>1.02451</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;adjacent/borrowed&#95;slice</td>
<td>15.7246 ns [15.6783, 15.7712] (0.95 confidence)</td>
<td>15.8556 ns [15.7926, 15.9201] (0.95 confidence)</td>
<td>0.991738</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;adjacent/construct&#95;then&#95;angle</td>
<td>14.8081 ns [14.7655, 14.8517] (0.95 confidence)</td>
<td>14.6802 ns [14.6598, 14.7014] (0.95 confidence)</td>
<td>1.00871</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;adjacent/prepared&#95;vector</td>
<td>13.3689 ns [13.3592, 13.3791] (0.95 confidence)</td>
<td>13.3696 ns [13.3604, 13.378] (0.95 confidence)</td>
<td>0.999948</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;alternating/borrowed&#95;slice</td>
<td>15.9397 ns [15.89, 15.9912] (0.95 confidence)</td>
<td>16.0018 ns [15.9531, 16.0523] (0.95 confidence)</td>
<td>0.996119</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;alternating/construct&#95;then&#95;angle</td>
<td>14.8977 ns [14.8754, 14.9215] (0.95 confidence)</td>
<td>14.9478 ns [14.9164, 14.981] (0.95 confidence)</td>
<td>0.996648</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;alternating/prepared&#95;vector</td>
<td>13.6248 ns [13.614, 13.6362] (0.95 confidence)</td>
<td>13.6431 ns [13.6264, 13.6619] (0.95 confidence)</td>
<td>0.998662</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;antipodal/borrowed&#95;slice</td>
<td>12.0311 ns [11.9967, 12.0659] (0.95 confidence)</td>
<td>12.0545 ns [12.0212, 12.0881] (0.95 confidence)</td>
<td>0.998054</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;antipodal/construct&#95;then&#95;angle</td>
<td>10.9538 ns [10.9225, 10.9865] (0.95 confidence)</td>
<td>11.035 ns [10.9911, 11.0873] (0.95 confidence)</td>
<td>0.99264</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;antipodal/prepared&#95;vector</td>
<td>9.58397 ns [9.57329, 9.59429] (0.95 confidence)</td>
<td>9.60452 ns [9.59446, 9.61407] (0.95 confidence)</td>
<td>0.99786</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;identical/borrowed&#95;slice</td>
<td>10.2761 ns [10.228, 10.3238] (0.95 confidence)</td>
<td>10.2726 ns [10.2206, 10.3241] (0.95 confidence)</td>
<td>1.00033</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;identical/construct&#95;then&#95;angle</td>
<td>9.1445 ns [9.10557, 9.18411] (0.95 confidence)</td>
<td>9.17413 ns [9.13751, 9.21007] (0.95 confidence)</td>
<td>0.996771</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;identical/prepared&#95;vector</td>
<td>7.68934 ns [7.68013, 7.69921] (0.95 confidence)</td>
<td>7.69515 ns [7.68402, 7.70676] (0.95 confidence)</td>
<td>0.999245</td>
</tr>
<tr>
<td>angle&#95;d4/subnormal/borrowed&#95;slice</td>
<td>12.7736 ns [12.7271, 12.8228] (0.95 confidence)</td>
<td>12.6895 ns [12.6553, 12.7279] (0.95 confidence)</td>
<td>1.00662</td>
</tr>
<tr>
<td>angle&#95;d4/subnormal/construct&#95;then&#95;angle</td>
<td>11.1613 ns [11.1307, 11.196] (0.95 confidence)</td>
<td>11.3018 ns [11.2456, 11.3626] (0.95 confidence)</td>
<td>0.987562</td>
</tr>
<tr>
<td>angle&#95;d4/subnormal/prepared&#95;vector</td>
<td>10.3313 ns [10.313, 10.3525] (0.95 confidence)</td>
<td>10.317 ns [10.3013, 10.3344] (0.95 confidence)</td>
<td>1.00139</td>
</tr>
<tr>
<td>angle&#95;d4/subnormal&#95;axes/borrowed&#95;slice</td>
<td>16.2116 ns [16.1315, 16.2922] (0.95 confidence)</td>
<td>15.715 ns [15.6576, 15.7821] (0.95 confidence)</td>
<td>1.0316</td>
</tr>
<tr>
<td>angle&#95;d4/subnormal&#95;axes/construct&#95;then&#95;angle</td>
<td>15.5269 ns [15.4045, 15.6544] (0.95 confidence)</td>
<td>14.4597 ns [14.3991, 14.5235] (0.95 confidence)</td>
<td>1.07381</td>
</tr>
<tr>
<td>angle&#95;d4/subnormal&#95;axes/prepared&#95;vector</td>
<td>13.4469 ns [13.3121, 13.6145] (0.95 confidence)</td>
<td>13.0393 ns [13.0071, 13.074] (0.95 confidence)</td>
<td>1.03125</td>
</tr>
<tr>
<td>angle&#95;d5/dense/borrowed&#95;slice</td>
<td>47.7267 ns [47.6657, 47.7816] (0.95 confidence)</td>
<td>47.0784 ns [47.0051, 47.1552] (0.95 confidence)</td>
<td>1.01377</td>
</tr>
<tr>
<td>angle&#95;d5/dense/construct&#95;then&#95;angle</td>
<td>45.3421 ns [45.289, 45.3923] (0.95 confidence)</td>
<td>44.9372 ns [44.8853, 44.9865] (0.95 confidence)</td>
<td>1.00901</td>
</tr>
<tr>
<td>angle&#95;d5/dense/prepared&#95;vector</td>
<td>44.9307 ns [44.8705, 44.9878] (0.95 confidence)</td>
<td>44.1311 ns [44.0617, 44.1959] (0.95 confidence)</td>
<td>1.01812</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;antipodal/borrowed&#95;slice</td>
<td>36.1699 ns [36.1348, 36.2041] (0.95 confidence)</td>
<td>36.0093 ns [35.96, 36.0606] (0.95 confidence)</td>
<td>1.00446</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;antipodal/construct&#95;then&#95;angle</td>
<td>34.2868 ns [34.2416, 34.3331] (0.95 confidence)</td>
<td>34.2123 ns [34.1769, 34.2467] (0.95 confidence)</td>
<td>1.00218</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;antipodal/prepared&#95;vector</td>
<td>34.1877 ns [34.1483, 34.2257] (0.95 confidence)</td>
<td>33.8926 ns [33.8568, 33.9269] (0.95 confidence)</td>
<td>1.00871</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;identical/borrowed&#95;slice</td>
<td>34.9007 ns [34.8602, 34.9413] (0.95 confidence)</td>
<td>34.7863 ns [34.7241, 34.8587] (0.95 confidence)</td>
<td>1.00329</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;identical/construct&#95;then&#95;angle</td>
<td>33.1249 ns [33.0852, 33.1636] (0.95 confidence)</td>
<td>32.8842 ns [32.848, 32.9185] (0.95 confidence)</td>
<td>1.00732</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;identical/prepared&#95;vector</td>
<td>32.7468 ns [32.7051, 32.787] (0.95 confidence)</td>
<td>32.7182 ns [32.6066, 32.8475] (0.95 confidence)</td>
<td>1.00088</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;varied/borrowed&#95;slice</td>
<td>68.2853 ns [68.1243, 68.4688] (0.95 confidence)</td>
<td>67.1009 ns [66.9026, 67.3045] (0.95 confidence)</td>
<td>1.01765</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;varied/construct&#95;then&#95;angle</td>
<td>65.7956 ns [65.6498, 65.9234] (0.95 confidence)</td>
<td>64.5051 ns [64.3749, 64.6276] (0.95 confidence)</td>
<td>1.02001</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;varied/prepared&#95;vector</td>
<td>63.9238 ns [63.7853, 64.049] (0.95 confidence)</td>
<td>62.8011 ns [62.682, 62.913] (0.95 confidence)</td>
<td>1.01788</td>
</tr>
<tr>
<td>angle&#95;d5/maximum&#95;small&#95;angle/borrowed&#95;slice</td>
<td>14.2103 ns [14.1678, 14.2509] (0.95 confidence)</td>
<td>14.2075 ns [14.1674, 14.2455] (0.95 confidence)</td>
<td>1.0002</td>
</tr>
<tr>
<td>angle&#95;d5/maximum&#95;small&#95;angle/construct&#95;then&#95;angle</td>
<td>12.6245 ns [12.5817, 12.6682] (0.95 confidence)</td>
<td>12.4734 ns [12.443, 12.5068] (0.95 confidence)</td>
<td>1.01211</td>
</tr>
<tr>
<td>angle&#95;d5/maximum&#95;small&#95;angle/prepared&#95;vector</td>
<td>11.1744 ns [11.1665, 11.1818] (0.95 confidence)</td>
<td>11.1995 ns [11.1915, 11.2088] (0.95 confidence)</td>
<td>0.997759</td>
</tr>
<tr>
<td>angle&#95;d5/mixed&#95;scale/borrowed&#95;slice</td>
<td>15.8141 ns [15.7956, 15.83] (0.95 confidence)</td>
<td>15.7265 ns [15.6906, 15.7604] (0.95 confidence)</td>
<td>1.00557</td>
</tr>
<tr>
<td>angle&#95;d5/mixed&#95;scale/construct&#95;then&#95;angle</td>
<td>14.2249 ns [14.1758, 14.2744] (0.95 confidence)</td>
<td>14.074 ns [14.0367, 14.1132] (0.95 confidence)</td>
<td>1.01072</td>
</tr>
<tr>
<td>angle&#95;d5/mixed&#95;scale/prepared&#95;vector</td>
<td>12.6869 ns [12.6795, 12.6945] (0.95 confidence)</td>
<td>12.7103 ns [12.6954, 12.7276] (0.95 confidence)</td>
<td>0.998156</td>
</tr>
<tr>
<td>angle&#95;d5/near&#95;antipodal/borrowed&#95;slice</td>
<td>19.8478 ns [19.8179, 19.8753] (0.95 confidence)</td>
<td>19.8097 ns [19.7856, 19.8333] (0.95 confidence)</td>
<td>1.00192</td>
</tr>
<tr>
<td>angle&#95;d5/near&#95;antipodal/construct&#95;then&#95;angle</td>
<td>17.6725 ns [17.6613, 17.6831] (0.95 confidence)</td>
<td>17.676 ns [17.6648, 17.6866] (0.95 confidence)</td>
<td>0.999801</td>
</tr>
<tr>
<td>angle&#95;d5/near&#95;antipodal/prepared&#95;vector</td>
<td>16.8015 ns [16.7744, 16.8284] (0.95 confidence)</td>
<td>16.7885 ns [16.7556, 16.8256] (0.95 confidence)</td>
<td>1.00077</td>
</tr>
<tr>
<td>angle&#95;d5/near&#95;parallel/borrowed&#95;slice</td>
<td>14.1822 ns [14.1334, 14.2286] (0.95 confidence)</td>
<td>14.2336 ns [14.1899, 14.2754] (0.95 confidence)</td>
<td>0.99639</td>
</tr>
<tr>
<td>angle&#95;d5/near&#95;parallel/construct&#95;then&#95;angle</td>
<td>12.6613 ns [12.6131, 12.7107] (0.95 confidence)</td>
<td>12.4922 ns [12.455, 12.5307] (0.95 confidence)</td>
<td>1.01354</td>
</tr>
<tr>
<td>angle&#95;d5/near&#95;parallel/prepared&#95;vector</td>
<td>11.2022 ns [11.1896, 11.2191] (0.95 confidence)</td>
<td>11.1837 ns [11.1734, 11.1941] (0.95 confidence)</td>
<td>1.00166</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;adjacent/borrowed&#95;slice</td>
<td>17.4226 ns [17.3891, 17.4603] (0.95 confidence)</td>
<td>17.3744 ns [17.3331, 17.4127] (0.95 confidence)</td>
<td>1.00277</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;adjacent/construct&#95;then&#95;angle</td>
<td>16.3601 ns [16.3175, 16.405] (0.95 confidence)</td>
<td>16.2813 ns [16.2376, 16.3278] (0.95 confidence)</td>
<td>1.00484</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;adjacent/prepared&#95;vector</td>
<td>14.6343 ns [14.6227, 14.6466] (0.95 confidence)</td>
<td>14.6521 ns [14.6332, 14.6727] (0.95 confidence)</td>
<td>0.998784</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;alternating/borrowed&#95;slice</td>
<td>28.1884 ns [28.1501, 28.2234] (0.95 confidence)</td>
<td>28.1227 ns [28.0889, 28.1537] (0.95 confidence)</td>
<td>1.00234</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;alternating/construct&#95;then&#95;angle</td>
<td>27.1155 ns [27.0735, 27.1643] (0.95 confidence)</td>
<td>27.0145 ns [26.9562, 27.0709] (0.95 confidence)</td>
<td>1.00374</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;alternating/prepared&#95;vector</td>
<td>26.4136 ns [26.38, 26.4448] (0.95 confidence)</td>
<td>26.4571 ns [26.4255, 26.4883] (0.95 confidence)</td>
<td>0.998358</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;antipodal/borrowed&#95;slice</td>
<td>13.3268 ns [13.3027, 13.3516] (0.95 confidence)</td>
<td>13.2539 ns [13.2288, 13.2761] (0.95 confidence)</td>
<td>1.0055</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;antipodal/construct&#95;then&#95;angle</td>
<td>11.9503 ns [11.9192, 11.9831] (0.95 confidence)</td>
<td>12.0676 ns [12.0181, 12.113] (0.95 confidence)</td>
<td>0.990279</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;antipodal/prepared&#95;vector</td>
<td>10.3896 ns [10.3703, 10.4122] (0.95 confidence)</td>
<td>10.3206 ns [10.3042, 10.3369] (0.95 confidence)</td>
<td>1.00668</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;identical/borrowed&#95;slice</td>
<td>11.7651 ns [11.7196, 11.8138] (0.95 confidence)</td>
<td>11.7849 ns [11.7328, 11.8375] (0.95 confidence)</td>
<td>0.998325</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;identical/construct&#95;then&#95;angle</td>
<td>10.6241 ns [10.5865, 10.6622] (0.95 confidence)</td>
<td>10.5082 ns [10.4837, 10.5334] (0.95 confidence)</td>
<td>1.01103</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;identical/prepared&#95;vector</td>
<td>8.7473 ns [8.72899, 8.76668] (0.95 confidence)</td>
<td>8.73996 ns [8.72131, 8.75961] (0.95 confidence)</td>
<td>1.00084</td>
</tr>
<tr>
<td>angle&#95;d5/subnormal/borrowed&#95;slice</td>
<td>14.1608 ns [14.1176, 14.2028] (0.95 confidence)</td>
<td>14.2438 ns [14.2022, 14.2834] (0.95 confidence)</td>
<td>0.994174</td>
</tr>
<tr>
<td>angle&#95;d5/subnormal/construct&#95;then&#95;angle</td>
<td>12.6195 ns [12.569, 12.6715] (0.95 confidence)</td>
<td>12.4758 ns [12.4427, 12.5116] (0.95 confidence)</td>
<td>1.01152</td>
</tr>
<tr>
<td>angle&#95;d5/subnormal/prepared&#95;vector</td>
<td>11.1908 ns [11.1838, 11.1976] (0.95 confidence)</td>
<td>11.1878 ns [11.1789, 11.1969] (0.95 confidence)</td>
<td>1.00027</td>
</tr>
<tr>
<td>angle&#95;d5/subnormal&#95;axes/borrowed&#95;slice</td>
<td>17.5451 ns [17.5044, 17.5859] (0.95 confidence)</td>
<td>17.4325 ns [17.376, 17.4898] (0.95 confidence)</td>
<td>1.00646</td>
</tr>
<tr>
<td>angle&#95;d5/subnormal&#95;axes/construct&#95;then&#95;angle</td>
<td>15.8395 ns [15.8038, 15.8762] (0.95 confidence)</td>
<td>15.7992 ns [15.747, 15.8543] (0.95 confidence)</td>
<td>1.00256</td>
</tr>
<tr>
<td>angle&#95;d5/subnormal&#95;axes/prepared&#95;vector</td>
<td>14.0943 ns [14.0472, 14.1454] (0.95 confidence)</td>
<td>14.1227 ns [14.1007, 14.1446] (0.95 confidence)</td>
<td>0.997983</td>
</tr>
<tr>
<td>angle&#95;d6/dense/borrowed&#95;slice</td>
<td>61.1727 ns [61.0467, 61.3252] (0.95 confidence)</td>
<td>61.0514 ns [60.9578, 61.1445] (0.95 confidence)</td>
<td>1.00199</td>
</tr>
<tr>
<td>angle&#95;d6/dense/construct&#95;then&#95;angle</td>
<td>58.7899 ns [58.7232, 58.8525] (0.95 confidence)</td>
<td>58.7123 ns [58.5918, 58.8476] (0.95 confidence)</td>
<td>1.00132</td>
</tr>
<tr>
<td>angle&#95;d6/dense/prepared&#95;vector</td>
<td>57.8693 ns [57.7629, 57.9922] (0.95 confidence)</td>
<td>57.8061 ns [57.7137, 57.8932] (0.95 confidence)</td>
<td>1.00109</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;antipodal/borrowed&#95;slice</td>
<td>48.0267 ns [47.941, 48.1201] (0.95 confidence)</td>
<td>47.9023 ns [47.8186, 47.971] (0.95 confidence)</td>
<td>1.0026</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;antipodal/construct&#95;then&#95;angle</td>
<td>45.5231 ns [45.4751, 45.5696] (0.95 confidence)</td>
<td>45.371 ns [45.2402, 45.5585] (0.95 confidence)</td>
<td>1.00335</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;antipodal/prepared&#95;vector</td>
<td>44.9155 ns [44.8574, 44.9673] (0.95 confidence)</td>
<td>44.8055 ns [44.7338, 44.8668] (0.95 confidence)</td>
<td>1.00245</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;identical/borrowed&#95;slice</td>
<td>46.8549 ns [46.8077, 46.8995] (0.95 confidence)</td>
<td>46.7352 ns [46.66, 46.8136] (0.95 confidence)</td>
<td>1.00256</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;identical/construct&#95;then&#95;angle</td>
<td>44.343 ns [44.3012, 44.3818] (0.95 confidence)</td>
<td>44.1749 ns [44.1038, 44.236] (0.95 confidence)</td>
<td>1.00381</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;identical/prepared&#95;vector</td>
<td>43.6673 ns [43.6065, 43.7237] (0.95 confidence)</td>
<td>43.7582 ns [43.7047, 43.8146] (0.95 confidence)</td>
<td>0.997923</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;varied/borrowed&#95;slice</td>
<td>90.9937 ns [90.7422, 91.2585] (0.95 confidence)</td>
<td>90.7588 ns [90.4571, 91.1008] (0.95 confidence)</td>
<td>1.00259</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;varied/construct&#95;then&#95;angle</td>
<td>88.7591 ns [88.5631, 88.9591] (0.95 confidence)</td>
<td>88.9377 ns [88.7323, 89.136] (0.95 confidence)</td>
<td>0.997991</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;varied/prepared&#95;vector</td>
<td>87.9048 ns [87.7421, 88.0568] (0.95 confidence)</td>
<td>87.9603 ns [87.7501, 88.18] (0.95 confidence)</td>
<td>0.999369</td>
</tr>
<tr>
<td>angle&#95;d6/maximum&#95;small&#95;angle/borrowed&#95;slice</td>
<td>15.6658 ns [15.6167, 15.7181] (0.95 confidence)</td>
<td>15.5423 ns [15.5214, 15.5604] (0.95 confidence)</td>
<td>1.00794</td>
</tr>
<tr>
<td>angle&#95;d6/maximum&#95;small&#95;angle/construct&#95;then&#95;angle</td>
<td>13.9099 ns [13.8746, 13.9462] (0.95 confidence)</td>
<td>13.8294 ns [13.7871, 13.8736] (0.95 confidence)</td>
<td>1.00582</td>
</tr>
<tr>
<td>angle&#95;d6/maximum&#95;small&#95;angle/prepared&#95;vector</td>
<td>12.0592 ns [11.996, 12.1316] (0.95 confidence)</td>
<td>11.9924 ns [11.9538, 12.0447] (0.95 confidence)</td>
<td>1.00557</td>
</tr>
<tr>
<td>angle&#95;d6/mixed&#95;scale/borrowed&#95;slice</td>
<td>17.0732 ns [17.0461, 17.0994] (0.95 confidence)</td>
<td>17.0375 ns [17.0012, 17.0732] (0.95 confidence)</td>
<td>1.0021</td>
</tr>
<tr>
<td>angle&#95;d6/mixed&#95;scale/construct&#95;then&#95;angle</td>
<td>15.27 ns [15.2348, 15.3079] (0.95 confidence)</td>
<td>15.2552 ns [15.216, 15.2961] (0.95 confidence)</td>
<td>1.00096</td>
</tr>
<tr>
<td>angle&#95;d6/mixed&#95;scale/prepared&#95;vector</td>
<td>13.4091 ns [13.3787, 13.4422] (0.95 confidence)</td>
<td>13.4115 ns [13.3881, 13.4377] (0.95 confidence)</td>
<td>0.999815</td>
</tr>
<tr>
<td>angle&#95;d6/near&#95;antipodal/borrowed&#95;slice</td>
<td>21.3178 ns [21.2959, 21.3386] (0.95 confidence)</td>
<td>21.2087 ns [21.1874, 21.2295] (0.95 confidence)</td>
<td>1.00514</td>
</tr>
<tr>
<td>angle&#95;d6/near&#95;antipodal/construct&#95;then&#95;angle</td>
<td>18.8336 ns [18.8071, 18.8575] (0.95 confidence)</td>
<td>18.8281 ns [18.811, 18.844] (0.95 confidence)</td>
<td>1.00029</td>
</tr>
<tr>
<td>angle&#95;d6/near&#95;antipodal/prepared&#95;vector</td>
<td>17.7254 ns [17.7097, 17.7405] (0.95 confidence)</td>
<td>17.6601 ns [17.637, 17.6829] (0.95 confidence)</td>
<td>1.0037</td>
</tr>
<tr>
<td>angle&#95;d6/near&#95;parallel/borrowed&#95;slice</td>
<td>15.5688 ns [15.5509, 15.5844] (0.95 confidence)</td>
<td>15.5598 ns [15.5287, 15.589] (0.95 confidence)</td>
<td>1.00058</td>
</tr>
<tr>
<td>angle&#95;d6/near&#95;parallel/construct&#95;then&#95;angle</td>
<td>13.8473 ns [13.8093, 13.8872] (0.95 confidence)</td>
<td>13.8306 ns [13.7901, 13.8722] (0.95 confidence)</td>
<td>1.00121</td>
</tr>
<tr>
<td>angle&#95;d6/near&#95;parallel/prepared&#95;vector</td>
<td>11.9469 ns [11.9346, 11.9604] (0.95 confidence)</td>
<td>11.9461 ns [11.9369, 11.9556] (0.95 confidence)</td>
<td>1.00007</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;adjacent/borrowed&#95;slice</td>
<td>19.3309 ns [19.2476, 19.421] (0.95 confidence)</td>
<td>19.1524 ns [19.1011, 19.2099] (0.95 confidence)</td>
<td>1.00932</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;adjacent/construct&#95;then&#95;angle</td>
<td>17.7706 ns [17.7109, 17.833] (0.95 confidence)</td>
<td>17.7536 ns [17.7045, 17.8037] (0.95 confidence)</td>
<td>1.00096</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;adjacent/prepared&#95;vector</td>
<td>15.8072 ns [15.7523, 15.8615] (0.95 confidence)</td>
<td>15.7969 ns [15.7644, 15.8301] (0.95 confidence)</td>
<td>1.00065</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;alternating/borrowed&#95;slice</td>
<td>29.8506 ns [29.8212, 29.8792] (0.95 confidence)</td>
<td>29.653 ns [29.6151, 29.69] (0.95 confidence)</td>
<td>1.00666</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;alternating/construct&#95;then&#95;angle</td>
<td>28.1209 ns [28.0954, 28.1443] (0.95 confidence)</td>
<td>28.0516 ns [27.9908, 28.1055] (0.95 confidence)</td>
<td>1.00247</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;alternating/prepared&#95;vector</td>
<td>27.211 ns [27.1775, 27.2442] (0.95 confidence)</td>
<td>27.0883 ns [27.049, 27.1241] (0.95 confidence)</td>
<td>1.00453</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;antipodal/borrowed&#95;slice</td>
<td>14.2784 ns [14.2609, 14.2946] (0.95 confidence)</td>
<td>14.2793 ns [14.2571, 14.3014] (0.95 confidence)</td>
<td>0.999937</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;antipodal/construct&#95;then&#95;angle</td>
<td>13.0152 ns [12.9831, 13.0487] (0.95 confidence)</td>
<td>13.0718 ns [13.034, 13.1113] (0.95 confidence)</td>
<td>0.99567</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;antipodal/prepared&#95;vector</td>
<td>11.0537 ns [11.0403, 11.0661] (0.95 confidence)</td>
<td>10.9989 ns [10.9855, 11.0115] (0.95 confidence)</td>
<td>1.00498</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;identical/borrowed&#95;slice</td>
<td>13.041 ns [12.9869, 13.0981] (0.95 confidence)</td>
<td>13.0598 ns [13.0017, 13.121] (0.95 confidence)</td>
<td>0.99856</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;identical/construct&#95;then&#95;angle</td>
<td>11.7303 ns [11.6813, 11.7788] (0.95 confidence)</td>
<td>11.6461 ns [11.6063, 11.6876] (0.95 confidence)</td>
<td>1.00724</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;identical/prepared&#95;vector</td>
<td>9.6753 ns [9.66353, 9.68795] (0.95 confidence)</td>
<td>9.67578 ns [9.65505, 9.69489] (0.95 confidence)</td>
<td>0.99995</td>
</tr>
<tr>
<td>angle&#95;d6/subnormal/borrowed&#95;slice</td>
<td>15.5699 ns [15.5271, 15.6121] (0.95 confidence)</td>
<td>15.5197 ns [15.4502, 15.5967] (0.95 confidence)</td>
<td>1.00324</td>
</tr>
<tr>
<td>angle&#95;d6/subnormal/construct&#95;then&#95;angle</td>
<td>13.7711 ns [13.7432, 13.8004] (0.95 confidence)</td>
<td>13.7983 ns [13.7582, 13.8404] (0.95 confidence)</td>
<td>0.998034</td>
</tr>
<tr>
<td>angle&#95;d6/subnormal/prepared&#95;vector</td>
<td>11.9608 ns [11.944, 11.9795] (0.95 confidence)</td>
<td>11.9538 ns [11.9334, 11.9763] (0.95 confidence)</td>
<td>1.00059</td>
</tr>
<tr>
<td>angle&#95;d6/subnormal&#95;axes/borrowed&#95;slice</td>
<td>19.0063 ns [18.9556, 19.0561] (0.95 confidence)</td>
<td>18.9521 ns [18.9027, 19.0011] (0.95 confidence)</td>
<td>1.00286</td>
</tr>
<tr>
<td>angle&#95;d6/subnormal&#95;axes/construct&#95;then&#95;angle</td>
<td>17.056 ns [17.0248, 17.0891] (0.95 confidence)</td>
<td>16.9883 ns [16.9497, 17.0285] (0.95 confidence)</td>
<td>1.00399</td>
</tr>
<tr>
<td>angle&#95;d6/subnormal&#95;axes/prepared&#95;vector</td>
<td>15.0859 ns [15.0568, 15.1097] (0.95 confidence)</td>
<td>15.1518 ns [15.1165, 15.1958] (0.95 confidence)</td>
<td>0.995647</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;det</td>
<td>391.684 ns [389.504, 393.724] (0.95 confidence)</td>
<td>398.35 ns [396.017, 400.664] (0.95 confidence)</td>
<td>0.983267</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>2.86973 ns [2.8656, 2.87374] (0.95 confidence)</td>
<td>2.94379 ns [2.93923, 2.94911] (0.95 confidence)</td>
<td>0.974843</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>3.26726 ns [3.26354, 3.27077] (0.95 confidence)</td>
<td>3.44469 ns [3.42698, 3.46321] (0.95 confidence)</td>
<td>0.948491</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>412.416 ns [407.184, 417.876] (0.95 confidence)</td>
<td>429.352 ns [421.822, 437.121] (0.95 confidence)</td>
<td>0.960554</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;dot</td>
<td>2.32641 ns [2.3243, 2.32845] (0.95 confidence)</td>
<td>2.7103 ns [2.64177, 2.78454] (0.95 confidence)</td>
<td>0.858361</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;inf&#95;norm</td>
<td>32.4255 ns [32.3644, 32.4908] (0.95 confidence)</td>
<td>45.75 ns [44.7763, 46.7913] (0.95 confidence)</td>
<td>0.708753</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;ldlt</td>
<td>341.879 ns [341.549, 342.192] (0.95 confidence)</td>
<td>359.998 ns [358.333, 361.71] (0.95 confidence)</td>
<td>0.949668</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;ldlt&#95;solve</td>
<td>382.579 ns [382.339, 382.796] (0.95 confidence)</td>
<td>501.869 ns [489.271, 514.967] (0.95 confidence)</td>
<td>0.762308</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;lu</td>
<td>436.44 ns [428.45, 444.35] (0.95 confidence)</td>
<td>523.726 ns [504.483, 546.232] (0.95 confidence)</td>
<td>0.833336</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;norm2</td>
<td>9.98409 ns [9.94708, 10.0214] (0.95 confidence)</td>
<td>13.8797 ns [13.4084, 14.3935] (0.95 confidence)</td>
<td>0.71933</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>9.25205 ns [9.23116, 9.27233] (0.95 confidence)</td>
<td>12.2532 ns [11.9509, 12.5949] (0.95 confidence)</td>
<td>0.755073</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>9.20438 ns [9.185, 9.22325] (0.95 confidence)</td>
<td>12.2501 ns [11.9743, 12.5679] (0.95 confidence)</td>
<td>0.751374</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>4.72069 ns [4.65912, 4.78042] (0.95 confidence)</td>
<td>6.24027 ns [6.11081, 6.37593] (0.95 confidence)</td>
<td>0.756488</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>8.23076 ns [8.17742, 8.29916] (0.95 confidence)</td>
<td>11.1961 ns [10.9881, 11.4206] (0.95 confidence)</td>
<td>0.735148</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;norm2&#95;sq</td>
<td>2.07295 ns [2.05975, 2.08381] (0.95 confidence)</td>
<td>2.88837 ns [2.7964, 2.99063] (0.95 confidence)</td>
<td>0.717689</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>27.7257 ns [27.7027, 27.7473] (0.95 confidence)</td>
<td>28.9772 ns [28.8564, 29.1059] (0.95 confidence)</td>
<td>0.956811</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>188.275 ns [187.955, 188.58] (0.95 confidence)</td>
<td>253.651 ns [247.736, 260.203] (0.95 confidence)</td>
<td>0.742258</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;det</td>
<td>0.573078 ns [0.572499, 0.574016] (0.95 confidence)</td>
<td>0.590489 ns [0.58856, 0.592388] (0.95 confidence)</td>
<td>0.970514</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>0.503753 ns [0.501194, 0.506576] (0.95 confidence)</td>
<td>0.5208 ns [0.518097, 0.523744] (0.95 confidence)</td>
<td>0.967267</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>0.519794 ns [0.519609, 0.520014] (0.95 confidence)</td>
<td>0.538361 ns [0.537946, 0.538776] (0.95 confidence)</td>
<td>0.96551</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>1.6074 ns [1.59662, 1.61947] (0.95 confidence)</td>
<td>1.58718 ns [1.58615, 1.58824] (0.95 confidence)</td>
<td>1.01274</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;dot</td>
<td>0.602431 ns [0.59763, 0.60807] (0.95 confidence)</td>
<td>0.635015 ns [0.630049, 0.641155] (0.95 confidence)</td>
<td>0.948688</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;inf&#95;norm</td>
<td>0.70995 ns [0.708756, 0.711209] (0.95 confidence)</td>
<td>0.713905 ns [0.712993, 0.714947] (0.95 confidence)</td>
<td>0.99446</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;ldlt</td>
<td>6.53809 ns [6.52044, 6.55632] (0.95 confidence)</td>
<td>5.22534 ns [5.186, 5.27056] (0.95 confidence)</td>
<td>1.25123</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;ldlt&#95;solve</td>
<td>9.35495 ns [9.3327, 9.3818] (0.95 confidence)</td>
<td>10.929 ns [10.8668, 10.9931] (0.95 confidence)</td>
<td>0.855973</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;lu</td>
<td>1.39335 ns [1.38958, 1.39703] (0.95 confidence)</td>
<td>1.41995 ns [1.40785, 1.43482] (0.95 confidence)</td>
<td>0.981264</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;norm2</td>
<td>1.16128 ns [1.15483, 1.16767] (0.95 confidence)</td>
<td>1.21219 ns [1.20567, 1.21914] (0.95 confidence)</td>
<td>0.957996</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>1.29993 ns [1.29736, 1.30233] (0.95 confidence)</td>
<td>1.41489 ns [1.41125, 1.41899] (0.95 confidence)</td>
<td>0.91875</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>1.29858 ns [1.29631, 1.30071] (0.95 confidence)</td>
<td>1.40731 ns [1.40575, 1.40868] (0.95 confidence)</td>
<td>0.922737</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>1.12322 ns [1.1216, 1.12462] (0.95 confidence)</td>
<td>1.205 ns [1.19643, 1.21368] (0.95 confidence)</td>
<td>0.93213</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>1.30077 ns [1.29825, 1.30315] (0.95 confidence)</td>
<td>1.38207 ns [1.38009, 1.38411] (0.95 confidence)</td>
<td>0.941175</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;norm2&#95;sq</td>
<td>0.393521 ns [0.392783, 0.394261] (0.95 confidence)</td>
<td>0.406928 ns [0.405998, 0.40782] (0.95 confidence)</td>
<td>0.967055</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>1.22152 ns [1.21358, 1.22811] (0.95 confidence)</td>
<td>1.26211 ns [1.25513, 1.26804] (0.95 confidence)</td>
<td>0.967836</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>1.25523 ns [1.25145, 1.25921] (0.95 confidence)</td>
<td>1.33524 ns [1.32558, 1.34562] (0.95 confidence)</td>
<td>0.940076</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;det</td>
<td>1.3543 ns [1.34844, 1.36054] (0.95 confidence)</td>
<td>1.37243 ns [1.36579, 1.37935] (0.95 confidence)</td>
<td>0.986789</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>0.59388 ns [0.593753, 0.594063] (0.95 confidence)</td>
<td>0.602854 ns [0.602352, 0.603406] (0.95 confidence)</td>
<td>0.985115</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>0.692647 ns [0.69237, 0.692891] (0.95 confidence)</td>
<td>0.704409 ns [0.703818, 0.705045] (0.95 confidence)</td>
<td>0.983302</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>8.82219 ns [8.77281, 8.88999] (0.95 confidence)</td>
<td>8.85175 ns [8.83623, 8.86761] (0.95 confidence)</td>
<td>0.99666</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;dot</td>
<td>0.666034 ns [0.665465, 0.666605] (0.95 confidence)</td>
<td>0.678469 ns [0.677082, 0.680104] (0.95 confidence)</td>
<td>0.981672</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;inf&#95;norm</td>
<td>1.36541 ns [1.36373, 1.36718] (0.95 confidence)</td>
<td>1.38113 ns [1.37983, 1.38276] (0.95 confidence)</td>
<td>0.988618</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;ldlt</td>
<td>8.40549 ns [8.39498, 8.41567] (0.95 confidence)</td>
<td>8.63834 ns [8.62853, 8.64859] (0.95 confidence)</td>
<td>0.973045</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;ldlt&#95;solve</td>
<td>9.23759 ns [9.22165, 9.25325] (0.95 confidence)</td>
<td>9.45039 ns [9.43011, 9.47116] (0.95 confidence)</td>
<td>0.977482</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;lu</td>
<td>8.18656 ns [8.1741, 8.19889] (0.95 confidence)</td>
<td>8.50011 ns [8.48412, 8.51579] (0.95 confidence)</td>
<td>0.963112</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;norm2</td>
<td>1.78876 ns [1.7602, 1.82116] (0.95 confidence)</td>
<td>1.87509 ns [1.83128, 1.92351] (0.95 confidence)</td>
<td>0.95396</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>1.80804 ns [1.80423, 1.8115] (0.95 confidence)</td>
<td>1.88712 ns [1.88424, 1.89075] (0.95 confidence)</td>
<td>0.958095</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>1.79832 ns [1.79302, 1.80275] (0.95 confidence)</td>
<td>1.88671 ns [1.88487, 1.88875] (0.95 confidence)</td>
<td>0.953148</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>1.38292 ns [1.37064, 1.39823] (0.95 confidence)</td>
<td>1.41875 ns [1.41041, 1.42774] (0.95 confidence)</td>
<td>0.974743</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>1.8195 ns [1.81696, 1.82195] (0.95 confidence)</td>
<td>1.91538 ns [1.9114, 1.91963] (0.95 confidence)</td>
<td>0.949944</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;norm2&#95;sq</td>
<td>0.444106 ns [0.443526, 0.44478] (0.95 confidence)</td>
<td>0.447199 ns [0.446664, 0.447707] (0.95 confidence)</td>
<td>0.993082</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>1.78497 ns [1.78239, 1.78787] (0.95 confidence)</td>
<td>1.81461 ns [1.81202, 1.81734] (0.95 confidence)</td>
<td>0.983668</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>2.05849 ns [2.05403, 2.06351] (0.95 confidence)</td>
<td>2.11214 ns [2.10543, 2.11892] (0.95 confidence)</td>
<td>0.974601</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;det</td>
<td>1983.66 ns [1980.04, 1987.49] (0.95 confidence)</td>
<td>3150.39 ns [3050.78, 3266.7] (0.95 confidence)</td>
<td>0.629655</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>7.69958 ns [7.68266, 7.7152] (0.95 confidence)</td>
<td>12.1318 ns [11.6465, 12.6343] (0.95 confidence)</td>
<td>0.634661</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>9.01955 ns [8.99851, 9.04137] (0.95 confidence)</td>
<td>9.49764 ns [9.41252, 9.59844] (0.95 confidence)</td>
<td>0.949661</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>2025.13 ns [2023.37, 2026.98] (0.95 confidence)</td>
<td>3008.27 ns [2947.2, 3071.75] (0.95 confidence)</td>
<td>0.673188</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;dot</td>
<td>3.9855 ns [3.98156, 3.98938] (0.95 confidence)</td>
<td>4.08336 ns [4.08191, 4.08496] (0.95 confidence)</td>
<td>0.976033</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;inf&#95;norm</td>
<td>126.974 ns [126.784, 127.157] (0.95 confidence)</td>
<td>131.04 ns [130.885, 131.198] (0.95 confidence)</td>
<td>0.968973</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;ldlt</td>
<td>2386.17 ns [2383.94, 2388.14] (0.95 confidence)</td>
<td>3529.04 ns [3481.13, 3582.37] (0.95 confidence)</td>
<td>0.676154</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;ldlt&#95;solve</td>
<td>2726.74 ns [2722.68, 2730.86] (0.95 confidence)</td>
<td>4047.02 ns [4004.84, 4089.32] (0.95 confidence)</td>
<td>0.673766</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;lu</td>
<td>2243.7 ns [2225.68, 2260.9] (0.95 confidence)</td>
<td>3438.6 ns [3368.19, 3512.27] (0.95 confidence)</td>
<td>0.652504</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;norm2</td>
<td>22.7597 ns [22.7175, 22.8009] (0.95 confidence)</td>
<td>24.8957 ns [24.8718, 24.9217] (0.95 confidence)</td>
<td>0.914201</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>18.5802 ns [18.5231, 18.6456] (0.95 confidence)</td>
<td>19.5911 ns [19.5692, 19.6138] (0.95 confidence)</td>
<td>0.948397</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>18.3796 ns [18.3482, 18.4096] (0.95 confidence)</td>
<td>19.5443 ns [19.5208, 19.5678] (0.95 confidence)</td>
<td>0.94041</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>9.08769 ns [9.06673, 9.10756] (0.95 confidence)</td>
<td>9.54477 ns [9.52097, 9.56902] (0.95 confidence)</td>
<td>0.952112</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>15.9759 ns [15.9467, 16.0036] (0.95 confidence)</td>
<td>16.8617 ns [16.8397, 16.8834] (0.95 confidence)</td>
<td>0.947469</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;norm2&#95;sq</td>
<td>3.92623 ns [3.92244, 3.92995] (0.95 confidence)</td>
<td>4.05008 ns [4.04184, 4.05941] (0.95 confidence)</td>
<td>0.96942</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>309.516 ns [309.042, 310.003] (0.95 confidence)</td>
<td>378.11 ns [364.962, 391.651] (0.95 confidence)</td>
<td>0.818588</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>652.945 ns [651.309, 654.681] (0.95 confidence)</td>
<td>1082.66 ns [1067.88, 1097.82] (0.95 confidence)</td>
<td>0.603091</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;det</td>
<td>2.54914 ns [2.54568, 2.55229] (0.95 confidence)</td>
<td>2.57477 ns [2.56823, 2.58126] (0.95 confidence)</td>
<td>0.990046</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>0.767517 ns [0.767024, 0.768345] (0.95 confidence)</td>
<td>0.795214 ns [0.793524, 0.796908] (0.95 confidence)</td>
<td>0.965171</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>0.866243 ns [0.865972, 0.866629] (0.95 confidence)</td>
<td>0.885531 ns [0.88344, 0.887893] (0.95 confidence)</td>
<td>0.978218</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>14.551 ns [14.5312, 14.5707] (0.95 confidence)</td>
<td>14.7305 ns [14.7041, 14.7561] (0.95 confidence)</td>
<td>0.987815</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;dot</td>
<td>0.680169 ns [0.679465, 0.680901] (0.95 confidence)</td>
<td>0.732166 ns [0.725292, 0.738966] (0.95 confidence)</td>
<td>0.928981</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;inf&#95;norm</td>
<td>2.33872 ns [2.33454, 2.34231] (0.95 confidence)</td>
<td>2.37985 ns [2.37867, 2.38094] (0.95 confidence)</td>
<td>0.982717</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;ldlt</td>
<td>17.8649 ns [17.8384, 17.8911] (0.95 confidence)</td>
<td>19.175 ns [19.0409, 19.3187] (0.95 confidence)</td>
<td>0.931678</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;ldlt&#95;solve</td>
<td>16.4911 ns [16.4649, 16.516] (0.95 confidence)</td>
<td>16.311 ns [16.2182, 16.4116] (0.95 confidence)</td>
<td>1.01104</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;lu</td>
<td>14.1844 ns [14.1592, 14.211] (0.95 confidence)</td>
<td>14.3557 ns [14.3352, 14.3769] (0.95 confidence)</td>
<td>0.988067</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;norm2</td>
<td>2.37299 ns [2.34809, 2.39789] (0.95 confidence)</td>
<td>2.44152 ns [2.41444, 2.46801] (0.95 confidence)</td>
<td>0.971933</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>2.37061 ns [2.36605, 2.37454] (0.95 confidence)</td>
<td>2.4743 ns [2.47082, 2.47951] (0.95 confidence)</td>
<td>0.958093</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>2.38624 ns [2.37929, 2.39387] (0.95 confidence)</td>
<td>2.52009 ns [2.50392, 2.5378] (0.95 confidence)</td>
<td>0.946889</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>1.53618 ns [1.50322, 1.57393] (0.95 confidence)</td>
<td>1.65713 ns [1.60418, 1.71393] (0.95 confidence)</td>
<td>0.927013</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>2.11669 ns [2.11316, 2.12011] (0.95 confidence)</td>
<td>2.30135 ns [2.2759, 2.32766] (0.95 confidence)</td>
<td>0.91976</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;norm2&#95;sq</td>
<td>0.458723 ns [0.457298, 0.460326] (0.95 confidence)</td>
<td>0.457908 ns [0.4573, 0.458561] (0.95 confidence)</td>
<td>1.00178</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>2.5255 ns [2.50881, 2.54266] (0.95 confidence)</td>
<td>2.66734 ns [2.64997, 2.68469] (0.95 confidence)</td>
<td>0.946822</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>4.00625 ns [4.0041, 4.00836] (0.95 confidence)</td>
<td>4.10176 ns [4.09108, 4.1122] (0.95 confidence)</td>
<td>0.976716</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;det</td>
<td>28.5189 ns [28.4948, 28.5437] (0.95 confidence)</td>
<td>29.6685 ns [29.6184, 29.7244] (0.95 confidence)</td>
<td>0.961249</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>0.875684 ns [0.875382, 0.876001] (0.95 confidence)</td>
<td>0.906381 ns [0.904346, 0.908477] (0.95 confidence)</td>
<td>0.966133</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>1.03161 ns [1.03118, 1.03212] (0.95 confidence)</td>
<td>1.05021 ns [1.04755, 1.05301] (0.95 confidence)</td>
<td>0.982287</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>33.9362 ns [33.8887, 33.9822] (0.95 confidence)</td>
<td>34.1695 ns [34.1509, 34.1884] (0.95 confidence)</td>
<td>0.993171</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;dot</td>
<td>0.83753 ns [0.835123, 0.839645] (0.95 confidence)</td>
<td>0.828302 ns [0.823507, 0.833046] (0.95 confidence)</td>
<td>1.01114</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;inf&#95;norm</td>
<td>3.47144 ns [3.45803, 3.48686] (0.95 confidence)</td>
<td>3.51744 ns [3.51426, 3.52089] (0.95 confidence)</td>
<td>0.986924</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;ldlt</td>
<td>24.4036 ns [24.3617, 24.4486] (0.95 confidence)</td>
<td>24.2925 ns [24.2291, 24.3563] (0.95 confidence)</td>
<td>1.00457</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;ldlt&#95;solve</td>
<td>30.166 ns [30.0727, 30.2746] (0.95 confidence)</td>
<td>28.8382 ns [28.8017, 28.8763] (0.95 confidence)</td>
<td>1.04604</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;lu</td>
<td>33.5855 ns [33.5316, 33.6422] (0.95 confidence)</td>
<td>33.4626 ns [33.43, 33.4968] (0.95 confidence)</td>
<td>1.00367</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;norm2</td>
<td>2.84419 ns [2.81802, 2.86982] (0.95 confidence)</td>
<td>2.94423 ns [2.91466, 2.97477] (0.95 confidence)</td>
<td>0.966022</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>3.13276 ns [3.12469, 3.14454] (0.95 confidence)</td>
<td>3.18139 ns [3.17825, 3.18489] (0.95 confidence)</td>
<td>0.984712</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>3.12279 ns [3.1174, 3.12796] (0.95 confidence)</td>
<td>3.1891 ns [3.18391, 3.19497] (0.95 confidence)</td>
<td>0.979206</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>1.7148 ns [1.71283, 1.7167] (0.95 confidence)</td>
<td>1.80231 ns [1.77863, 1.83146] (0.95 confidence)</td>
<td>0.951444</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>3.09549 ns [3.08887, 3.10136] (0.95 confidence)</td>
<td>2.65993 ns [2.65318, 2.66649] (0.95 confidence)</td>
<td>1.16375</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;norm2&#95;sq</td>
<td>0.525335 ns [0.524972, 0.525709] (0.95 confidence)</td>
<td>0.532895 ns [0.532092, 0.533737] (0.95 confidence)</td>
<td>0.985814</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>3.95182 ns [3.94389, 3.95942] (0.95 confidence)</td>
<td>3.95244 ns [3.93961, 3.9661] (0.95 confidence)</td>
<td>0.999845</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>6.0888 ns [6.08072, 6.09747] (0.95 confidence)</td>
<td>6.30068 ns [6.26394, 6.34028] (0.95 confidence)</td>
<td>0.966372</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;det</td>
<td>14378.5 ns [14355.4, 14403.2] (0.95 confidence)</td>
<td>14532.8 ns [14517.9, 14547.6] (0.95 confidence)</td>
<td>0.989386</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>23.1553 ns [22.9938, 23.3441] (0.95 confidence)</td>
<td>23.7147 ns [23.6771, 23.7521] (0.95 confidence)</td>
<td>0.97641</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>22.9523 ns [22.7837, 23.1362] (0.95 confidence)</td>
<td>23.7985 ns [23.7586, 23.84] (0.95 confidence)</td>
<td>0.964442</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>14343.8 ns [14318.4, 14370.4] (0.95 confidence)</td>
<td>14391.8 ns [14372.2, 14410.8] (0.95 confidence)</td>
<td>0.996671</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;dot</td>
<td>10.8454 ns [10.8223, 10.8675] (0.95 confidence)</td>
<td>11.2581 ns [11.2357, 11.2807] (0.95 confidence)</td>
<td>0.963338</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;inf&#95;norm</td>
<td>619.511 ns [618.288, 620.662] (0.95 confidence)</td>
<td>858.835 ns [833.178, 885.212] (0.95 confidence)</td>
<td>0.721339</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;ldlt</td>
<td>16872.4 ns [16856, 16888.1] (0.95 confidence)</td>
<td>17053.5 ns [17034.2, 17071.6] (0.95 confidence)</td>
<td>0.989379</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;ldlt&#95;solve</td>
<td>18943.1 ns [18916.2, 18969.3] (0.95 confidence)</td>
<td>19219.9 ns [19176.1, 19266.7] (0.95 confidence)</td>
<td>0.9856</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;lu</td>
<td>14946 ns [14925.3, 14967] (0.95 confidence)</td>
<td>15257 ns [15243.2, 15270.7] (0.95 confidence)</td>
<td>0.979612</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;norm2</td>
<td>57.0012 ns [56.9076, 57.0969] (0.95 confidence)</td>
<td>61.6936 ns [61.2342, 62.1684] (0.95 confidence)</td>
<td>0.92394</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>36.0511 ns [35.9855, 36.1172] (0.95 confidence)</td>
<td>38.4514 ns [38.3931, 38.51] (0.95 confidence)</td>
<td>0.937575</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>35.9015 ns [35.8456, 35.9563] (0.95 confidence)</td>
<td>41.582 ns [41.0365, 42.1781] (0.95 confidence)</td>
<td>0.86339</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>16.9744 ns [16.8939, 17.0918] (0.95 confidence)</td>
<td>18.6607 ns [18.6164, 18.7034] (0.95 confidence)</td>
<td>0.909635</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>32.1362 ns [32.0729, 32.196] (0.95 confidence)</td>
<td>39.2354 ns [38.3505, 40.1498] (0.95 confidence)</td>
<td>0.819063</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;norm2&#95;sq</td>
<td>10.7 ns [10.6766, 10.7248] (0.95 confidence)</td>
<td>11.2856 ns [11.228, 11.3481] (0.95 confidence)</td>
<td>0.948115</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>1014.4 ns [1011.62, 1018.04] (0.95 confidence)</td>
<td>1026.12 ns [1020.76, 1033.69] (0.95 confidence)</td>
<td>0.988578</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>2496.17 ns [2487.64, 2508.59] (0.95 confidence)</td>
<td>2707.49 ns [2702.14, 2713.11] (0.95 confidence)</td>
<td>0.921949</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;det</td>
<td>81.5576 ns [81.2869, 81.8042] (0.95 confidence)</td>
<td>83.0821 ns [82.5178, 83.6687] (0.95 confidence)</td>
<td>0.981651</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>1.30359 ns [1.30249, 1.30464] (0.95 confidence)</td>
<td>1.33515 ns [1.33475, 1.33559] (0.95 confidence)</td>
<td>0.976366</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;det&#95;from&#95;ldlt&#95;balanced&#95;range</td>
<td>8.57256 ns [8.55379, 8.58902] (0.95 confidence)</td>
<td>8.90404 ns [8.86592, 8.94866] (0.95 confidence)</td>
<td>0.962772</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>1.4627 ns [1.46114, 1.4642] (0.95 confidence)</td>
<td>1.58573 ns [1.57066, 1.60455] (0.95 confidence)</td>
<td>0.922413</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;det&#95;from&#95;lu&#95;balanced&#95;range</td>
<td>10.9809 ns [10.9562, 11.0055] (0.95 confidence)</td>
<td>11.1034 ns [11.0806, 11.1275] (0.95 confidence)</td>
<td>0.988966</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>89.415 ns [89.0393, 89.7797] (0.95 confidence)</td>
<td>91.5974 ns [91.0112, 92.2037] (0.95 confidence)</td>
<td>0.976174</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;dot</td>
<td>0.932825 ns [0.931183, 0.934473] (0.95 confidence)</td>
<td>0.945822 ns [0.942299, 0.949183] (0.95 confidence)</td>
<td>0.986259</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;inf&#95;norm</td>
<td>8.34857 ns [8.33811, 8.35842] (0.95 confidence)</td>
<td>8.58166 ns [8.57983, 8.58391] (0.95 confidence)</td>
<td>0.972838</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;ldlt</td>
<td>90.9195 ns [90.5288, 91.3144] (0.95 confidence)</td>
<td>92.275 ns [92.0114, 92.5517] (0.95 confidence)</td>
<td>0.985311</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;ldlt&#95;ill&#95;conditioned</td>
<td>90.4759 ns [90.2289, 90.7241] (0.95 confidence)</td>
<td>92.3699 ns [91.9908, 92.7727] (0.95 confidence)</td>
<td>0.979495</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;ldlt&#95;solve</td>
<td>102.518 ns [102.282, 102.749] (0.95 confidence)</td>
<td>104.744 ns [104.405, 105.08] (0.95 confidence)</td>
<td>0.978746</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;lu</td>
<td>89.6454 ns [89.071, 90.2158] (0.95 confidence)</td>
<td>90.8255 ns [90.1411, 91.5464] (0.95 confidence)</td>
<td>0.987007</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;lu&#95;ill&#95;conditioned</td>
<td>88.7934 ns [88.3849, 89.1739] (0.95 confidence)</td>
<td>89.9498 ns [89.5464, 90.3523] (0.95 confidence)</td>
<td>0.987143</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;lu&#95;pivoting</td>
<td>94.1729 ns [93.9437, 94.4047] (0.95 confidence)</td>
<td>95.1672 ns [94.9652, 95.3813] (0.95 confidence)</td>
<td>0.989552</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;norm2</td>
<td>4.46732 ns [4.45186, 4.48134] (0.95 confidence)</td>
<td>4.7318 ns [4.71974, 4.74515] (0.95 confidence)</td>
<td>0.944106</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>4.7304 ns [4.72019, 4.74053] (0.95 confidence)</td>
<td>5.01063 ns [5.00421, 5.01728] (0.95 confidence)</td>
<td>0.944072</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>4.7429 ns [4.73298, 4.75289] (0.95 confidence)</td>
<td>5.01744 ns [5.01014, 5.02509] (0.95 confidence)</td>
<td>0.945282</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>3.01616 ns [2.97141, 3.05648] (0.95 confidence)</td>
<td>3.07294 ns [3.01334, 3.12908] (0.95 confidence)</td>
<td>0.981523</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>4.17319 ns [4.16651, 4.17973] (0.95 confidence)</td>
<td>4.37293 ns [4.37045, 4.37544] (0.95 confidence)</td>
<td>0.954323</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;norm2&#95;sq</td>
<td>0.664128 ns [0.663156, 0.665041] (0.95 confidence)</td>
<td>0.669806 ns [0.668592, 0.671001] (0.95 confidence)</td>
<td>0.991523</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>8.09418 ns [8.08847, 8.10012] (0.95 confidence)</td>
<td>8.28953 ns [8.26933, 8.31069] (0.95 confidence)</td>
<td>0.976434</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>13.3082 ns [13.2882, 13.3254] (0.95 confidence)</td>
<td>14.0527 ns [13.9826, 14.1276] (0.95 confidence)</td>
<td>0.947022</td>
</tr>
<tr>
<td>gram/1x2/dependent/la&#95;stack</td>
<td>0.682578 ns [0.681125, 0.684076] (0.95 confidence)</td>
<td>0.68167 ns [0.680348, 0.683062] (0.95 confidence)</td>
<td>1.00133</td>
</tr>
<tr>
<td>gram/1x2/mixed&#95;scale/la&#95;stack</td>
<td>0.682629 ns [0.680801, 0.684383] (0.95 confidence)</td>
<td>0.680767 ns [0.679501, 0.682086] (0.95 confidence)</td>
<td>1.00274</td>
</tr>
<tr>
<td>gram/1x2/near&#95;dependent/la&#95;stack</td>
<td>0.683183 ns [0.681804, 0.684628] (0.95 confidence)</td>
<td>0.680946 ns [0.679532, 0.682439] (0.95 confidence)</td>
<td>1.00329</td>
</tr>
<tr>
<td>gram/1x2/orthogonal/la&#95;stack</td>
<td>0.683873 ns [0.682559, 0.685276] (0.95 confidence)</td>
<td>0.682397 ns [0.68083, 0.684185] (0.95 confidence)</td>
<td>1.00216</td>
</tr>
<tr>
<td>gram/2x2/dependent/la&#95;stack</td>
<td>1.31021 ns [1.30305, 1.31876] (0.95 confidence)</td>
<td>1.30045 ns [1.29562, 1.3052] (0.95 confidence)</td>
<td>1.00751</td>
</tr>
<tr>
<td>gram/2x2/mixed&#95;scale/la&#95;stack</td>
<td>1.30285 ns [1.29698, 1.30858] (0.95 confidence)</td>
<td>1.30349 ns [1.29752, 1.30936] (0.95 confidence)</td>
<td>0.999515</td>
</tr>
<tr>
<td>gram/2x2/near&#95;dependent/la&#95;stack</td>
<td>1.30529 ns [1.29996, 1.31059] (0.95 confidence)</td>
<td>1.30293 ns [1.2976, 1.30817] (0.95 confidence)</td>
<td>1.00181</td>
</tr>
<tr>
<td>gram/2x2/orthogonal/la&#95;stack</td>
<td>1.30438 ns [1.29899, 1.30963] (0.95 confidence)</td>
<td>1.30257 ns [1.29722, 1.30776] (0.95 confidence)</td>
<td>1.00139</td>
</tr>
<tr>
<td>gram/2x3/dependent/la&#95;stack</td>
<td>1.26812 ns [1.26637, 1.26985] (0.95 confidence)</td>
<td>1.26566 ns [1.26355, 1.26808] (0.95 confidence)</td>
<td>1.00194</td>
</tr>
<tr>
<td>gram/2x3/mixed&#95;scale/la&#95;stack</td>
<td>1.26387 ns [1.26221, 1.26554] (0.95 confidence)</td>
<td>1.26461 ns [1.26256, 1.26678] (0.95 confidence)</td>
<td>0.999415</td>
</tr>
<tr>
<td>gram/2x3/near&#95;dependent/la&#95;stack</td>
<td>1.2667 ns [1.26517, 1.26825] (0.95 confidence)</td>
<td>1.2642 ns [1.26242, 1.26598] (0.95 confidence)</td>
<td>1.00198</td>
</tr>
<tr>
<td>gram/2x3/orthogonal/la&#95;stack</td>
<td>1.2723 ns [1.26778, 1.27735] (0.95 confidence)</td>
<td>1.26319 ns [1.26125, 1.26501] (0.95 confidence)</td>
<td>1.00721</td>
</tr>
<tr>
<td>gram/3x3/dependent/la&#95;stack</td>
<td>2.99297 ns [2.98826, 2.9978] (0.95 confidence)</td>
<td>2.97394 ns [2.96951, 2.97859] (0.95 confidence)</td>
<td>1.0064</td>
</tr>
<tr>
<td>gram/3x3/mixed&#95;scale/la&#95;stack</td>
<td>2.98612 ns [2.98306, 2.98911] (0.95 confidence)</td>
<td>2.97393 ns [2.9703, 2.97748] (0.95 confidence)</td>
<td>1.0041</td>
</tr>
<tr>
<td>gram/3x3/near&#95;dependent/la&#95;stack</td>
<td>2.9759 ns [2.97243, 2.9793] (0.95 confidence)</td>
<td>2.99527 ns [2.98104, 3.01135] (0.95 confidence)</td>
<td>0.993532</td>
</tr>
<tr>
<td>gram/3x3/orthogonal/la&#95;stack</td>
<td>2.98605 ns [2.98021, 2.99162] (0.95 confidence)</td>
<td>2.97285 ns [2.96855, 2.97693] (0.95 confidence)</td>
<td>1.00444</td>
</tr>
<tr>
<td>gram/3x4/dependent/la&#95;stack</td>
<td>3.03711 ns [3.0332, 3.04102] (0.95 confidence)</td>
<td>3.03628 ns [3.0326, 3.04021] (0.95 confidence)</td>
<td>1.00027</td>
</tr>
<tr>
<td>gram/3x4/mixed&#95;scale/la&#95;stack</td>
<td>3.03978 ns [3.03671, 3.04261] (0.95 confidence)</td>
<td>3.0685 ns [3.05978, 3.07746] (0.95 confidence)</td>
<td>0.990641</td>
</tr>
<tr>
<td>gram/3x4/near&#95;dependent/la&#95;stack</td>
<td>3.03467 ns [3.03104, 3.03794] (0.95 confidence)</td>
<td>3.03041 ns [3.02562, 3.03552] (0.95 confidence)</td>
<td>1.0014</td>
</tr>
<tr>
<td>gram/3x4/orthogonal/la&#95;stack</td>
<td>3.04024 ns [3.03729, 3.04301] (0.95 confidence)</td>
<td>3.02176 ns [3.0171, 3.02606] (0.95 confidence)</td>
<td>1.00612</td>
</tr>
<tr>
<td>gram/4x4/dependent/la&#95;stack</td>
<td>7.51501 ns [7.47247, 7.55763] (0.95 confidence)</td>
<td>7.4737 ns [7.43146, 7.51575] (0.95 confidence)</td>
<td>1.00553</td>
</tr>
<tr>
<td>gram/4x4/mixed&#95;scale/la&#95;stack</td>
<td>7.49296 ns [7.45136, 7.53455] (0.95 confidence)</td>
<td>7.48904 ns [7.4428, 7.53764] (0.95 confidence)</td>
<td>1.00052</td>
</tr>
<tr>
<td>gram/4x4/near&#95;dependent/la&#95;stack</td>
<td>7.49142 ns [7.44872, 7.53378] (0.95 confidence)</td>
<td>7.48457 ns [7.44227, 7.52719] (0.95 confidence)</td>
<td>1.00091</td>
</tr>
<tr>
<td>gram/4x4/orthogonal/la&#95;stack</td>
<td>7.50964 ns [7.46785, 7.55155] (0.95 confidence)</td>
<td>7.45919 ns [7.41575, 7.50187] (0.95 confidence)</td>
<td>1.00676</td>
</tr>
<tr>
<td>gram/4x5/dependent/la&#95;stack</td>
<td>9.43283 ns [9.41932, 9.44538] (0.95 confidence)</td>
<td>9.37057 ns [9.35455, 9.38551] (0.95 confidence)</td>
<td>1.00664</td>
</tr>
<tr>
<td>gram/4x5/mixed&#95;scale/la&#95;stack</td>
<td>9.4391 ns [9.41601, 9.46405] (0.95 confidence)</td>
<td>9.31743 ns [9.2903, 9.34311] (0.95 confidence)</td>
<td>1.01306</td>
</tr>
<tr>
<td>gram/4x5/near&#95;dependent/la&#95;stack</td>
<td>9.39514 ns [9.37949, 9.40987] (0.95 confidence)</td>
<td>9.3629 ns [9.33859, 9.38412] (0.95 confidence)</td>
<td>1.00344</td>
</tr>
<tr>
<td>gram/4x5/orthogonal/la&#95;stack</td>
<td>9.42203 ns [9.40909, 9.43523] (0.95 confidence)</td>
<td>9.30345 ns [9.27321, 9.33255] (0.95 confidence)</td>
<td>1.01275</td>
</tr>
<tr>
<td>gram/5x5/dependent/la&#95;stack</td>
<td>18.2382 ns [18.1661, 18.3052] (0.95 confidence)</td>
<td>18.1427 ns [18.0653, 18.216] (0.95 confidence)</td>
<td>1.00526</td>
</tr>
<tr>
<td>gram/5x5/mixed&#95;scale/la&#95;stack</td>
<td>18.1815 ns [18.0849, 18.2729] (0.95 confidence)</td>
<td>18.3702 ns [18.2977, 18.4344] (0.95 confidence)</td>
<td>0.989725</td>
</tr>
<tr>
<td>gram/5x5/near&#95;dependent/la&#95;stack</td>
<td>18.2596 ns [18.202, 18.3119] (0.95 confidence)</td>
<td>18.2641 ns [18.1803, 18.3411] (0.95 confidence)</td>
<td>0.999753</td>
</tr>
<tr>
<td>gram/5x5/orthogonal/la&#95;stack</td>
<td>18.3817 ns [18.3392, 18.4213] (0.95 confidence)</td>
<td>18.3004 ns [18.2342, 18.3601] (0.95 confidence)</td>
<td>1.00444</td>
</tr>
<tr>
<td>gram/5x6/dependent/la&#95;stack</td>
<td>19.1052 ns [19.0508, 19.1559] (0.95 confidence)</td>
<td>19.0428 ns [18.9974, 19.0839] (0.95 confidence)</td>
<td>1.00327</td>
</tr>
<tr>
<td>gram/5x6/mixed&#95;scale/la&#95;stack</td>
<td>18.9901 ns [18.9149, 19.0613] (0.95 confidence)</td>
<td>18.9538 ns [18.9059, 18.9979] (0.95 confidence)</td>
<td>1.00192</td>
</tr>
<tr>
<td>gram/5x6/near&#95;dependent/la&#95;stack</td>
<td>19.1485 ns [19.109, 19.1845] (0.95 confidence)</td>
<td>19.0148 ns [18.9546, 19.08] (0.95 confidence)</td>
<td>1.00703</td>
</tr>
<tr>
<td>gram/5x6/orthogonal/la&#95;stack</td>
<td>18.9908 ns [18.9283, 19.0502] (0.95 confidence)</td>
<td>19.0222 ns [18.9617, 19.0822] (0.95 confidence)</td>
<td>0.998345</td>
</tr>
<tr>
<td>gram/6x6/dependent/la&#95;stack</td>
<td>29.3685 ns [29.2049, 29.5354] (0.95 confidence)</td>
<td>29.2156 ns [29.0583, 29.378] (0.95 confidence)</td>
<td>1.00523</td>
</tr>
<tr>
<td>gram/6x6/mixed&#95;scale/la&#95;stack</td>
<td>29.4024 ns [29.2284, 29.5798] (0.95 confidence)</td>
<td>29.1408 ns [28.9635, 29.3209] (0.95 confidence)</td>
<td>1.00898</td>
</tr>
<tr>
<td>gram/6x6/near&#95;dependent/la&#95;stack</td>
<td>29.4381 ns [29.2726, 29.6058] (0.95 confidence)</td>
<td>29.2616 ns [29.1039, 29.4232] (0.95 confidence)</td>
<td>1.00603</td>
</tr>
<tr>
<td>gram/6x6/orthogonal/la&#95;stack</td>
<td>29.5158 ns [29.3137, 29.7202] (0.95 confidence)</td>
<td>29.5091 ns [29.3322, 29.6891] (0.95 confidence)</td>
<td>1.00023</td>
</tr>
<tr>
<td>gram/6x7/dependent/la&#95;stack</td>
<td>32.0229 ns [31.9115, 32.1305] (0.95 confidence)</td>
<td>31.8304 ns [31.7323, 31.9269] (0.95 confidence)</td>
<td>1.00605</td>
</tr>
<tr>
<td>gram/6x7/mixed&#95;scale/la&#95;stack</td>
<td>32.1056 ns [31.9982, 32.2198] (0.95 confidence)</td>
<td>31.598 ns [31.4679, 31.7282] (0.95 confidence)</td>
<td>1.01606</td>
</tr>
<tr>
<td>gram/6x7/near&#95;dependent/la&#95;stack</td>
<td>32.0422 ns [31.9388, 32.1425] (0.95 confidence)</td>
<td>31.845 ns [31.7391, 31.9473] (0.95 confidence)</td>
<td>1.00619</td>
</tr>
<tr>
<td>gram/6x7/orthogonal/la&#95;stack</td>
<td>32.0387 ns [31.8964, 32.1826] (0.95 confidence)</td>
<td>32.1105 ns [31.9345, 32.3034] (0.95 confidence)</td>
<td>0.997764</td>
</tr>
<tr>
<td>gram/7x7/dependent/la&#95;stack</td>
<td>40.2045 ns [40.014, 40.4306] (0.95 confidence)</td>
<td>40.1644 ns [40.02, 40.336] (0.95 confidence)</td>
<td>1.001</td>
</tr>
<tr>
<td>gram/7x7/mixed&#95;scale/la&#95;stack</td>
<td>40.5847 ns [40.3496, 40.8475] (0.95 confidence)</td>
<td>40.4786 ns [40.2307, 40.7543] (0.95 confidence)</td>
<td>1.00262</td>
</tr>
<tr>
<td>gram/7x7/near&#95;dependent/la&#95;stack</td>
<td>40.7341 ns [40.3606, 41.1327] (0.95 confidence)</td>
<td>40.9667 ns [40.5984, 41.3495] (0.95 confidence)</td>
<td>0.994323</td>
</tr>
<tr>
<td>gram/7x7/orthogonal/la&#95;stack</td>
<td>40.4439 ns [40.1594, 40.7616] (0.95 confidence)</td>
<td>40.6523 ns [40.388, 40.9425] (0.95 confidence)</td>
<td>0.994874</td>
</tr>
<tr>
<td>gram/7x8/dependent/la&#95;stack</td>
<td>41.6452 ns [41.38, 41.9431] (0.95 confidence)</td>
<td>42.7477 ns [42.4171, 43.0873] (0.95 confidence)</td>
<td>0.974208</td>
</tr>
<tr>
<td>gram/7x8/mixed&#95;scale/la&#95;stack</td>
<td>42.134 ns [41.7668, 42.5248] (0.95 confidence)</td>
<td>41.8839 ns [41.6506, 42.1477] (0.95 confidence)</td>
<td>1.00597</td>
</tr>
<tr>
<td>gram/7x8/near&#95;dependent/la&#95;stack</td>
<td>41.3276 ns [41.1634, 41.521] (0.95 confidence)</td>
<td>42.0438 ns [41.7624, 42.3415] (0.95 confidence)</td>
<td>0.982967</td>
</tr>
<tr>
<td>gram/7x8/orthogonal/la&#95;stack</td>
<td>41.9343 ns [41.6078, 42.2885] (0.95 confidence)</td>
<td>41.26 ns [41.1768, 41.3404] (0.95 confidence)</td>
<td>1.01634</td>
</tr>
<tr>
<td>gram/8x8/dependent/la&#95;stack</td>
<td>51.9098 ns [51.7659, 52.073] (0.95 confidence)</td>
<td>52.1054 ns [51.9092, 52.3095] (0.95 confidence)</td>
<td>0.996246</td>
</tr>
<tr>
<td>gram/8x8/mixed&#95;scale/la&#95;stack</td>
<td>52.0772 ns [51.9251, 52.2436] (0.95 confidence)</td>
<td>51.7293 ns [51.57, 51.9081] (0.95 confidence)</td>
<td>1.00673</td>
</tr>
<tr>
<td>gram/8x8/near&#95;dependent/la&#95;stack</td>
<td>51.979 ns [51.8433, 52.1336] (0.95 confidence)</td>
<td>51.5699 ns [51.4998, 51.6308] (0.95 confidence)</td>
<td>1.00793</td>
</tr>
<tr>
<td>gram/8x8/orthogonal/la&#95;stack</td>
<td>52.2853 ns [52.0845, 52.4975] (0.95 confidence)</td>
<td>51.7252 ns [51.6281, 51.8408] (0.95 confidence)</td>
<td>1.01083</td>
</tr>
<tr>
<td>interval&#95;d3/dense&#95;inconclusive</td>
<td>41.2767 ns [41.1844, 41.3701] (0.95 confidence)</td>
<td>41.2741 ns [41.1454, 41.4186] (0.95 confidence)</td>
<td>1.00006</td>
</tr>
<tr>
<td>interval&#95;d3/dense&#95;sign</td>
<td>35.7748 ns [35.5818, 35.9768] (0.95 confidence)</td>
<td>35.6414 ns [35.3415, 35.9667] (0.95 confidence)</td>
<td>1.00374</td>
</tr>
<tr>
<td>interval&#95;d3/lifted&#95;assembly</td>
<td>8.3619 ns [8.35134, 8.37352] (0.95 confidence)</td>
<td>8.34687 ns [8.3396, 8.35364] (0.95 confidence)</td>
<td>1.0018</td>
</tr>
<tr>
<td>interval&#95;d3/lifted&#95;assembly&#95;and&#95;sign</td>
<td>37.9832 ns [37.8557, 38.1253] (0.95 confidence)</td>
<td>37.6689 ns [37.5362, 37.7983] (0.95 confidence)</td>
<td>1.00834</td>
</tr>
<tr>
<td>interval&#95;d3/sparse&#95;sign</td>
<td>30.9902 ns [30.7383, 31.2353] (0.95 confidence)</td>
<td>31.2528 ns [31.1147, 31.3835] (0.95 confidence)</td>
<td>0.991595</td>
</tr>
<tr>
<td>interval&#95;d4/dense&#95;inconclusive</td>
<td>100.037 ns [99.3807, 100.712] (0.95 confidence)</td>
<td>101.382 ns [100.743, 102.013] (0.95 confidence)</td>
<td>0.986732</td>
</tr>
<tr>
<td>interval&#95;d4/dense&#95;sign</td>
<td>97.1614 ns [96.6075, 97.7277] (0.95 confidence)</td>
<td>99.0904 ns [98.5172, 99.6513] (0.95 confidence)</td>
<td>0.980533</td>
</tr>
<tr>
<td>interval&#95;d4/lifted&#95;assembly</td>
<td>13.2154 ns [13.1967, 13.2335] (0.95 confidence)</td>
<td>13.1328 ns [13.1005, 13.1605] (0.95 confidence)</td>
<td>1.00629</td>
</tr>
<tr>
<td>interval&#95;d4/lifted&#95;assembly&#95;and&#95;sign</td>
<td>81.5099 ns [81.3615, 81.656] (0.95 confidence)</td>
<td>81.5751 ns [81.1104, 82.0736] (0.95 confidence)</td>
<td>0.9992</td>
</tr>
<tr>
<td>interval&#95;d4/sparse&#95;sign</td>
<td>70.5952 ns [69.8388, 71.3672] (0.95 confidence)</td>
<td>67.6613 ns [67.3249, 68.0119] (0.95 confidence)</td>
<td>1.04336</td>
</tr>
<tr>
<td>interval&#95;d5/dense&#95;inconclusive</td>
<td>236.103 ns [235.18, 237.042] (0.95 confidence)</td>
<td>236.485 ns [235.405, 237.596] (0.95 confidence)</td>
<td>0.998384</td>
</tr>
<tr>
<td>interval&#95;d5/dense&#95;sign</td>
<td>224.122 ns [223.152, 225.138] (0.95 confidence)</td>
<td>227.299 ns [226.097, 228.517] (0.95 confidence)</td>
<td>0.986024</td>
</tr>
<tr>
<td>interval&#95;d5/lifted&#95;assembly</td>
<td>21.4495 ns [21.3983, 21.5009] (0.95 confidence)</td>
<td>21.434 ns [21.3837, 21.4835] (0.95 confidence)</td>
<td>1.00072</td>
</tr>
<tr>
<td>interval&#95;d5/lifted&#95;assembly&#95;and&#95;sign</td>
<td>192.183 ns [191.009, 193.389] (0.95 confidence)</td>
<td>189.711 ns [188.524, 190.916] (0.95 confidence)</td>
<td>1.01303</td>
</tr>
<tr>
<td>interval&#95;d5/sparse&#95;sign</td>
<td>155.12 ns [154.168, 156.141] (0.95 confidence)</td>
<td>151.64 ns [151.064, 152.224] (0.95 confidence)</td>
<td>1.02295</td>
</tr>
<tr>
<td>interval&#95;d6/dense&#95;inconclusive</td>
<td>543.906 ns [540.352, 547.553] (0.95 confidence)</td>
<td>525.28 ns [523.488, 527.114] (0.95 confidence)</td>
<td>1.03546</td>
</tr>
<tr>
<td>interval&#95;d6/dense&#95;sign</td>
<td>513.365 ns [512.101, 514.675] (0.95 confidence)</td>
<td>516.858 ns [514.576, 519.262] (0.95 confidence)</td>
<td>0.99324</td>
</tr>
<tr>
<td>interval&#95;d6/lifted&#95;assembly</td>
<td>30.5456 ns [30.505, 30.5852] (0.95 confidence)</td>
<td>30.4692 ns [30.3796, 30.5452] (0.95 confidence)</td>
<td>1.00251</td>
</tr>
<tr>
<td>interval&#95;d6/lifted&#95;assembly&#95;and&#95;sign</td>
<td>377.6 ns [376.651, 378.582] (0.95 confidence)</td>
<td>376.409 ns [375.143, 377.591] (0.95 confidence)</td>
<td>1.00317</td>
</tr>
<tr>
<td>interval&#95;d6/sparse&#95;sign</td>
<td>364.042 ns [362.007, 366.167] (0.95 confidence)</td>
<td>363.4 ns [361.465, 365.416] (0.95 confidence)</td>
<td>1.00177</td>
</tr>
<tr>
<td>interval&#95;d7/dense&#95;inconclusive</td>
<td>1098.07 ns [1091.33, 1104.87] (0.95 confidence)</td>
<td>1099.37 ns [1090.21, 1110.67] (0.95 confidence)</td>
<td>0.998818</td>
</tr>
<tr>
<td>interval&#95;d7/dense&#95;sign</td>
<td>1091.68 ns [1084.41, 1098.88] (0.95 confidence)</td>
<td>1053.82 ns [1048.5, 1059.21] (0.95 confidence)</td>
<td>1.03592</td>
</tr>
<tr>
<td>interval&#95;d7/lifted&#95;assembly</td>
<td>37.8794 ns [37.784, 37.9743] (0.95 confidence)</td>
<td>37.145 ns [37.0674, 37.2236] (0.95 confidence)</td>
<td>1.01977</td>
</tr>
<tr>
<td>interval&#95;d7/lifted&#95;assembly&#95;and&#95;sign</td>
<td>784.808 ns [783.053, 786.55] (0.95 confidence)</td>
<td>767.662 ns [766.333, 768.993] (0.95 confidence)</td>
<td>1.02233</td>
</tr>
<tr>
<td>interval&#95;d7/sparse&#95;sign</td>
<td>741.873 ns [740.095, 743.647] (0.95 confidence)</td>
<td>741.447 ns [738.535, 744.305] (0.95 confidence)</td>
<td>1.00057</td>
</tr>
<tr>
<td>interval&#95;det&#95;sign/d4&#95;conclusive&#95;lifted</td>
<td>70.6195 ns [70.2112, 71.0413] (0.95 confidence)</td>
<td>66.1521 ns [65.9054, 66.4039] (0.95 confidence)</td>
<td>1.06753</td>
</tr>
<tr>
<td>interval&#95;det&#95;sign/d4&#95;inconclusive&#95;lifted</td>
<td>67.3171 ns [66.8706, 67.7708] (0.95 confidence)</td>
<td>64.9882 ns [64.7042, 65.3087] (0.95 confidence)</td>
<td>1.03584</td>
</tr>
<tr>
<td>interval&#95;det&#95;sign/d7&#95;conclusive&#95;lifted</td>
<td>735.617 ns [733.839, 737.574] (0.95 confidence)</td>
<td>741.729 ns [739.602, 743.858] (0.95 confidence)</td>
<td>0.99176</td>
</tr>
<tr>
<td>interval&#95;scalar/cancellation/add</td>
<td>2.10278 ns [2.10061, 2.10489] (0.95 confidence)</td>
<td>2.09915 ns [2.09635, 2.10191] (0.95 confidence)</td>
<td>1.00173</td>
</tr>
<tr>
<td>interval&#95;scalar/cancellation/multiply</td>
<td>1.86815 ns [1.86426, 1.87205] (0.95 confidence)</td>
<td>1.87655 ns [1.86306, 1.89268] (0.95 confidence)</td>
<td>0.995524</td>
</tr>
<tr>
<td>interval&#95;scalar/cancellation/square</td>
<td>1.8492 ns [1.84648, 1.85171] (0.95 confidence)</td>
<td>1.84204 ns [1.83791, 1.84552] (0.95 confidence)</td>
<td>1.00389</td>
</tr>
<tr>
<td>interval&#95;scalar/exact/add</td>
<td>2.10359 ns [2.10129, 2.10578] (0.95 confidence)</td>
<td>2.09602 ns [2.093, 2.09887] (0.95 confidence)</td>
<td>1.00361</td>
</tr>
<tr>
<td>interval&#95;scalar/exact/multiply</td>
<td>2.06695 ns [2.05532, 2.07866] (0.95 confidence)</td>
<td>2.07932 ns [2.06544, 2.0929] (0.95 confidence)</td>
<td>0.994051</td>
</tr>
<tr>
<td>interval&#95;scalar/exact/square</td>
<td>1.85267 ns [1.85029, 1.85496] (0.95 confidence)</td>
<td>1.83689 ns [1.83186, 1.84106] (0.95 confidence)</td>
<td>1.00859</td>
</tr>
<tr>
<td>interval&#95;scalar/point/add</td>
<td>2.16825 ns [2.16648, 2.16993] (0.95 confidence)</td>
<td>2.16462 ns [2.15992, 2.16982] (0.95 confidence)</td>
<td>1.00168</td>
</tr>
<tr>
<td>interval&#95;scalar/point/multiply</td>
<td>2.57676 ns [2.57083, 2.58316] (0.95 confidence)</td>
<td>2.58536 ns [2.57584, 2.59594] (0.95 confidence)</td>
<td>0.996671</td>
</tr>
<tr>
<td>interval&#95;scalar/point/square</td>
<td>2.10414 ns [2.10152, 2.10674] (0.95 confidence)</td>
<td>2.09625 ns [2.09316, 2.09902] (0.95 confidence)</td>
<td>1.00376</td>
</tr>
<tr>
<td>interval&#95;scalar/subnormal/add</td>
<td>2.40035 ns [2.39642, 2.40488] (0.95 confidence)</td>
<td>2.39284 ns [2.38858, 2.39731] (0.95 confidence)</td>
<td>1.00314</td>
</tr>
<tr>
<td>interval&#95;scalar/subnormal/multiply</td>
<td>2.85655 ns [2.84382, 2.86957] (0.95 confidence)</td>
<td>2.85885 ns [2.84704, 2.87083] (0.95 confidence)</td>
<td>0.999199</td>
</tr>
<tr>
<td>interval&#95;scalar/subnormal/square</td>
<td>2.36529 ns [2.36257, 2.36778] (0.95 confidence)</td>
<td>2.36584 ns [2.3626, 2.36873] (0.95 confidence)</td>
<td>0.999766</td>
</tr>
<tr>
<td>interval&#95;scalar/zero/add</td>
<td>0.933207 ns [0.931905, 0.934456] (0.95 confidence)</td>
<td>0.931825 ns [0.930365, 0.933385] (0.95 confidence)</td>
<td>1.00148</td>
</tr>
<tr>
<td>interval&#95;scalar/zero/multiply</td>
<td>1.74264 ns [1.73979, 1.74527] (0.95 confidence)</td>
<td>1.74142 ns [1.73866, 1.74405] (0.95 confidence)</td>
<td>1.0007</td>
</tr>
<tr>
<td>interval&#95;scalar/zero/square</td>
<td>0.876898 ns [0.87603, 0.877738] (0.95 confidence)</td>
<td>0.874881 ns [0.874172, 0.875563] (0.95 confidence)</td>
<td>1.00231</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;batch/origin/construction</td>
<td>2.05389 ns [2.05063, 2.05702] (0.95 confidence)</td>
<td>2.05269 ns [2.04737, 2.0583] (0.95 confidence)</td>
<td>1.00058</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;batch/origin/prepared</td>
<td>13.8164 ns [13.7973, 13.8368] (0.95 confidence)</td>
<td>13.8307 ns [13.8134, 13.8489] (0.95 confidence)</td>
<td>0.998962</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;batch/origin/with&#95;construction</td>
<td>15.2937 ns [15.2784, 15.3081] (0.95 confidence)</td>
<td>15.3204 ns [15.305, 15.3351] (0.95 confidence)</td>
<td>0.998261</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;batch/translated/construction</td>
<td>2.04998 ns [2.04647, 2.05341] (0.95 confidence)</td>
<td>2.04798 ns [2.04401, 2.05182] (0.95 confidence)</td>
<td>1.00098</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;batch/translated/prepared</td>
<td>40.5787 ns [40.4176, 40.7841] (0.95 confidence)</td>
<td>40.3567 ns [40.2523, 40.4643] (0.95 confidence)</td>
<td>1.0055</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;batch/translated/with&#95;construction</td>
<td>42.3089 ns [42.1953, 42.4317] (0.95 confidence)</td>
<td>42.2322 ns [42.0973, 42.3742] (0.95 confidence)</td>
<td>1.00182</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/cancellation/difference</td>
<td>10.4423 ns [10.3547, 10.5474] (0.95 confidence)</td>
<td>10.4905 ns [10.3873, 10.6097] (0.95 confidence)</td>
<td>0.995408</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/cancellation/difference&#95;endpoints</td>
<td>11.8123 ns [11.6766, 11.9638] (0.95 confidence)</td>
<td>11.7821 ns [11.6424, 11.9367] (0.95 confidence)</td>
<td>1.00256</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/cancellation/dot</td>
<td>5.9389 ns [5.92768, 5.95154] (0.95 confidence)</td>
<td>5.90255 ns [5.8775, 5.93028] (0.95 confidence)</td>
<td>1.00616</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/cancellation/dot&#95;endpoints</td>
<td>6.99773 ns [6.98233, 7.01339] (0.95 confidence)</td>
<td>6.98192 ns [6.96885, 6.9951] (0.95 confidence)</td>
<td>1.00226</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense/difference</td>
<td>10.3867 ns [10.3258, 10.4557] (0.95 confidence)</td>
<td>10.2778 ns [10.2178, 10.3501] (0.95 confidence)</td>
<td>1.0106</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense/difference&#95;endpoints</td>
<td>12.529 ns [12.4681, 12.5992] (0.95 confidence)</td>
<td>12.5417 ns [12.4584, 12.6343] (0.95 confidence)</td>
<td>0.998994</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense/dot</td>
<td>5.50716 ns [5.49696, 5.51766] (0.95 confidence)</td>
<td>5.44194 ns [5.4332, 5.45036] (0.95 confidence)</td>
<td>1.01199</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense/dot&#95;endpoints</td>
<td>6.44325 ns [6.42734, 6.46092] (0.95 confidence)</td>
<td>6.4071 ns [6.3882, 6.42647] (0.95 confidence)</td>
<td>1.00564</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense&#95;inexact/difference</td>
<td>10.4385 ns [10.356, 10.5331] (0.95 confidence)</td>
<td>10.4307 ns [10.3322, 10.5398] (0.95 confidence)</td>
<td>1.00075</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense&#95;inexact/difference&#95;endpoints</td>
<td>11.6341 ns [11.5421, 11.7335] (0.95 confidence)</td>
<td>11.4383 ns [11.4257, 11.4506] (0.95 confidence)</td>
<td>1.01711</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense&#95;inexact/dot</td>
<td>5.51369 ns [5.50088, 5.53041] (0.95 confidence)</td>
<td>5.5009 ns [5.49105, 5.5109] (0.95 confidence)</td>
<td>1.00232</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense&#95;inexact/dot&#95;endpoints</td>
<td>6.47965 ns [6.45795, 6.50321] (0.95 confidence)</td>
<td>6.50576 ns [6.48345, 6.52969] (0.95 confidence)</td>
<td>0.995985</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/overflow/dot</td>
<td>1.07103 ns [1.06931, 1.07287] (0.95 confidence)</td>
<td>1.07241 ns [1.07081, 1.07411] (0.95 confidence)</td>
<td>0.998715</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/sparse/difference</td>
<td>4.82726 ns [4.78954, 4.86885] (0.95 confidence)</td>
<td>4.92439 ns [4.88802, 4.96799] (0.95 confidence)</td>
<td>0.980275</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/sparse/difference&#95;endpoints</td>
<td>5.55793 ns [5.54307, 5.57347] (0.95 confidence)</td>
<td>5.57058 ns [5.5523, 5.58962] (0.95 confidence)</td>
<td>0.997728</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/sparse/dot</td>
<td>2.2212 ns [2.21931, 2.22296] (0.95 confidence)</td>
<td>2.22212 ns [2.22055, 2.22364] (0.95 confidence)</td>
<td>0.999587</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/sparse/dot&#95;endpoints</td>
<td>2.22258 ns [2.22093, 2.22418] (0.95 confidence)</td>
<td>2.22101 ns [2.21849, 2.22307] (0.95 confidence)</td>
<td>1.00071</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/underflow/difference</td>
<td>4.45547 ns [4.45303, 4.4581] (0.95 confidence)</td>
<td>4.47291 ns [4.46474, 4.48188] (0.95 confidence)</td>
<td>0.996102</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/underflow/difference&#95;endpoints</td>
<td>4.60353 ns [4.59956, 4.60785] (0.95 confidence)</td>
<td>4.60903 ns [4.60383, 4.61486] (0.95 confidence)</td>
<td>0.998805</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/underflow/dot</td>
<td>1.47627 ns [1.4746, 1.47789] (0.95 confidence)</td>
<td>1.47532 ns [1.47361, 1.477] (0.95 confidence)</td>
<td>1.00065</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/underflow/dot&#95;endpoints</td>
<td>1.48258 ns [1.48191, 1.48322] (0.95 confidence)</td>
<td>1.48244 ns [1.48141, 1.48341] (0.95 confidence)</td>
<td>1.00009</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;batch/origin/construction</td>
<td>2.31911 ns [2.29162, 2.35836] (0.95 confidence)</td>
<td>2.24773 ns [2.23929, 2.25795] (0.95 confidence)</td>
<td>1.03176</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;batch/origin/prepared</td>
<td>24.7934 ns [24.6427, 24.9536] (0.95 confidence)</td>
<td>23.3008 ns [23.2397, 23.3677] (0.95 confidence)</td>
<td>1.06405</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;batch/origin/with&#95;construction</td>
<td>28.4576 ns [28.2989, 28.6217] (0.95 confidence)</td>
<td>27.6926 ns [27.6092, 27.7801] (0.95 confidence)</td>
<td>1.02763</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;batch/translated/construction</td>
<td>2.31279 ns [2.28857, 2.34616] (0.95 confidence)</td>
<td>2.22531 ns [2.2169, 2.23365] (0.95 confidence)</td>
<td>1.03931</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;batch/translated/prepared</td>
<td>70.4239 ns [70.1361, 70.7263] (0.95 confidence)</td>
<td>68.975 ns [68.6293, 69.3342] (0.95 confidence)</td>
<td>1.02101</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;batch/translated/with&#95;construction</td>
<td>72.7753 ns [72.4005, 73.1852] (0.95 confidence)</td>
<td>71.8662 ns [71.4612, 72.305] (0.95 confidence)</td>
<td>1.01265</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/cancellation/difference</td>
<td>11.077 ns [10.9896, 11.1784] (0.95 confidence)</td>
<td>10.9455 ns [10.9107, 10.9964] (0.95 confidence)</td>
<td>1.01202</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/cancellation/difference&#95;endpoints</td>
<td>12.3465 ns [12.2179, 12.4889] (0.95 confidence)</td>
<td>12.1853 ns [12.0644, 12.3223] (0.95 confidence)</td>
<td>1.01323</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/cancellation/dot</td>
<td>6.14417 ns [6.08718, 6.20827] (0.95 confidence)</td>
<td>6.0626 ns [6.02462, 6.10724] (0.95 confidence)</td>
<td>1.01345</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/cancellation/dot&#95;endpoints</td>
<td>6.5853 ns [6.56638, 6.60468] (0.95 confidence)</td>
<td>6.81176 ns [6.74464, 6.88204] (0.95 confidence)</td>
<td>0.966755</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense/difference</td>
<td>14.6808 ns [14.5564, 14.8098] (0.95 confidence)</td>
<td>14.3256 ns [14.2386, 14.4251] (0.95 confidence)</td>
<td>1.0248</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense/difference&#95;endpoints</td>
<td>16.6326 ns [16.5438, 16.7314] (0.95 confidence)</td>
<td>16.8537 ns [16.7453, 16.9669] (0.95 confidence)</td>
<td>0.986883</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense/dot</td>
<td>8.04597 ns [8.01025, 8.08442] (0.95 confidence)</td>
<td>7.88502 ns [7.85731, 7.91427] (0.95 confidence)</td>
<td>1.02041</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense/dot&#95;endpoints</td>
<td>8.4362 ns [8.38609, 8.49375] (0.95 confidence)</td>
<td>8.54698 ns [8.47975, 8.61845] (0.95 confidence)</td>
<td>0.987038</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense&#95;inexact/difference</td>
<td>14.2754 ns [14.2282, 14.3377] (0.95 confidence)</td>
<td>14.384 ns [14.298, 14.4837] (0.95 confidence)</td>
<td>0.992453</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense&#95;inexact/difference&#95;endpoints</td>
<td>16.1046 ns [16.0124, 16.2057] (0.95 confidence)</td>
<td>15.9139 ns [15.8468, 15.9899] (0.95 confidence)</td>
<td>1.01199</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense&#95;inexact/dot</td>
<td>7.99225 ns [7.97686, 8.00832] (0.95 confidence)</td>
<td>7.99989 ns [7.95096, 8.05163] (0.95 confidence)</td>
<td>0.999045</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense&#95;inexact/dot&#95;endpoints</td>
<td>9.37811 ns [9.33103, 9.42828] (0.95 confidence)</td>
<td>9.37876 ns [9.33516, 9.42657] (0.95 confidence)</td>
<td>0.99993</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/overflow/dot</td>
<td>0.999465 ns [0.998781, 1.00015] (0.95 confidence)</td>
<td>0.996919 ns [0.996008, 0.997849] (0.95 confidence)</td>
<td>1.00255</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/sparse/difference</td>
<td>5.45125 ns [5.44705, 5.45563] (0.95 confidence)</td>
<td>5.47839 ns [5.47143, 5.48661] (0.95 confidence)</td>
<td>0.995046</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/sparse/difference&#95;endpoints</td>
<td>6.41523 ns [6.4055, 6.42442] (0.95 confidence)</td>
<td>6.39055 ns [6.37213, 6.40885] (0.95 confidence)</td>
<td>1.00386</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/sparse/dot</td>
<td>2.71924 ns [2.71707, 2.72113] (0.95 confidence)</td>
<td>2.79124 ns [2.77148, 2.81229] (0.95 confidence)</td>
<td>0.974204</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/sparse/dot&#95;endpoints</td>
<td>2.76691 ns [2.75207, 2.78406] (0.95 confidence)</td>
<td>2.93586 ns [2.91429, 2.95669] (0.95 confidence)</td>
<td>0.942452</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/underflow/difference</td>
<td>4.71248 ns [4.70732, 4.71891] (0.95 confidence)</td>
<td>4.70413 ns [4.70161, 4.70656] (0.95 confidence)</td>
<td>1.00177</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/underflow/difference&#95;endpoints</td>
<td>5.05867 ns [5.03722, 5.08089] (0.95 confidence)</td>
<td>5.02823 ns [5.00933, 5.04821] (0.95 confidence)</td>
<td>1.00605</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/underflow/dot</td>
<td>3.12727 ns [3.1223, 3.13215] (0.95 confidence)</td>
<td>3.12643 ns [3.12051, 3.13322] (0.95 confidence)</td>
<td>1.00027</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/underflow/dot&#95;endpoints</td>
<td>3.16826 ns [3.1609, 3.17659] (0.95 confidence)</td>
<td>3.11707 ns [3.1094, 3.12382] (0.95 confidence)</td>
<td>1.01642</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4/dot&#95;bounded&#95;inconclusive</td>
<td>6.23804 ns [6.20487, 6.27689] (0.95 confidence)</td>
<td>6.31584 ns [6.2872, 6.34662] (0.95 confidence)</td>
<td>0.987681</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4/dot&#95;bounded&#95;well&#95;separated</td>
<td>9.65214 ns [9.59818, 9.71237] (0.95 confidence)</td>
<td>9.74712 ns [9.66637, 9.83418] (0.95 confidence)</td>
<td>0.990255</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4/dot&#95;difference&#95;bounded&#95;inconclusive</td>
<td>21.0788 ns [20.9926, 21.1785] (0.95 confidence)</td>
<td>21.0847 ns [21.0614, 21.1101] (0.95 confidence)</td>
<td>0.999721</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4/dot&#95;difference&#95;bounded&#95;well&#95;separated</td>
<td>17.0841 ns [16.9241, 17.2575] (0.95 confidence)</td>
<td>16.8704 ns [16.8001, 16.9533] (0.95 confidence)</td>
<td>1.01267</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4/dot&#95;plain&#95;well&#95;separated</td>
<td>0.685632 ns [0.683816, 0.687599] (0.95 confidence)</td>
<td>0.689734 ns [0.68753, 0.692137] (0.95 confidence)</td>
<td>0.994054</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;batch/origin/construction</td>
<td>2.3838 ns [2.37633, 2.39144] (0.95 confidence)</td>
<td>2.39347 ns [2.38512, 2.40253] (0.95 confidence)</td>
<td>0.99596</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;batch/origin/prepared</td>
<td>35.4566 ns [35.2655, 35.6555] (0.95 confidence)</td>
<td>36.1346 ns [35.9723, 36.2982] (0.95 confidence)</td>
<td>0.981237</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;batch/origin/with&#95;construction</td>
<td>42.2127 ns [42.015, 42.4175] (0.95 confidence)</td>
<td>42.6756 ns [42.4747, 42.8816] (0.95 confidence)</td>
<td>0.989154</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;batch/translated/construction</td>
<td>2.37647 ns [2.36995, 2.38309] (0.95 confidence)</td>
<td>2.40241 ns [2.39542, 2.40914] (0.95 confidence)</td>
<td>0.989201</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;batch/translated/prepared</td>
<td>112.409 ns [111.888, 112.962] (0.95 confidence)</td>
<td>112.25 ns [111.86, 112.661] (0.95 confidence)</td>
<td>1.00141</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;batch/translated/with&#95;construction</td>
<td>117.078 ns [116.719, 117.483] (0.95 confidence)</td>
<td>117.887 ns [117.337, 118.472] (0.95 confidence)</td>
<td>0.99313</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/cancellation/difference</td>
<td>11.9039 ns [11.8262, 11.988] (0.95 confidence)</td>
<td>12.0415 ns [11.9426, 12.1454] (0.95 confidence)</td>
<td>0.988573</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/cancellation/difference&#95;endpoints</td>
<td>13.3579 ns [13.2061, 13.514] (0.95 confidence)</td>
<td>13.0284 ns [12.9374, 13.1289] (0.95 confidence)</td>
<td>1.02529</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/cancellation/dot</td>
<td>6.6078 ns [6.56811, 6.64942] (0.95 confidence)</td>
<td>6.41641 ns [6.37875, 6.45979] (0.95 confidence)</td>
<td>1.02983</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/cancellation/dot&#95;endpoints</td>
<td>7.1415 ns [7.10787, 7.17818] (0.95 confidence)</td>
<td>7.06671 ns [7.0277, 7.11119] (0.95 confidence)</td>
<td>1.01058</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense/difference</td>
<td>19.0432 ns [18.9131, 19.182] (0.95 confidence)</td>
<td>18.3501 ns [18.2627, 18.4478] (0.95 confidence)</td>
<td>1.03777</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense/difference&#95;endpoints</td>
<td>21.5243 ns [21.3745, 21.6831] (0.95 confidence)</td>
<td>21.0657 ns [20.9272, 21.2358] (0.95 confidence)</td>
<td>1.02177</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense/dot</td>
<td>9.86211 ns [9.80434, 9.92736] (0.95 confidence)</td>
<td>9.60998 ns [9.54271, 9.68515] (0.95 confidence)</td>
<td>1.02624</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense/dot&#95;endpoints</td>
<td>11.363 ns [11.2918, 11.4378] (0.95 confidence)</td>
<td>11.0049 ns [10.9438, 11.0707] (0.95 confidence)</td>
<td>1.03253</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense&#95;inexact/difference</td>
<td>19.1299 ns [19.0077, 19.2651] (0.95 confidence)</td>
<td>19.0693 ns [18.9577, 19.1886] (0.95 confidence)</td>
<td>1.00318</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense&#95;inexact/difference&#95;endpoints</td>
<td>20.854 ns [20.7291, 20.988] (0.95 confidence)</td>
<td>20.7739 ns [20.659, 20.8955] (0.95 confidence)</td>
<td>1.00385</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense&#95;inexact/dot</td>
<td>9.8838 ns [9.81114, 9.96117] (0.95 confidence)</td>
<td>9.77345 ns [9.691, 9.86351] (0.95 confidence)</td>
<td>1.01129</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense&#95;inexact/dot&#95;endpoints</td>
<td>10.6052 ns [10.536, 10.6821] (0.95 confidence)</td>
<td>10.5181 ns [10.4283, 10.6179] (0.95 confidence)</td>
<td>1.00829</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/overflow/dot</td>
<td>1.0013 ns [1.00003, 1.00271] (0.95 confidence)</td>
<td>1.01063 ns [1.00863, 1.01274] (0.95 confidence)</td>
<td>0.990768</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/sparse/difference</td>
<td>6.43651 ns [6.41454, 6.45902] (0.95 confidence)</td>
<td>6.29831 ns [6.28773, 6.31236] (0.95 confidence)</td>
<td>1.02194</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/sparse/difference&#95;endpoints</td>
<td>7.49994 ns [7.47071, 7.52963] (0.95 confidence)</td>
<td>7.45148 ns [7.43231, 7.47156] (0.95 confidence)</td>
<td>1.0065</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/sparse/dot</td>
<td>3.2253 ns [3.2087, 3.24427] (0.95 confidence)</td>
<td>3.23677 ns [3.22184, 3.25283] (0.95 confidence)</td>
<td>0.996456</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/sparse/dot&#95;endpoints</td>
<td>3.26348 ns [3.25544, 3.27218] (0.95 confidence)</td>
<td>3.24784 ns [3.24088, 3.25619] (0.95 confidence)</td>
<td>1.00482</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/underflow/difference</td>
<td>5.79085 ns [5.76844, 5.8178] (0.95 confidence)</td>
<td>5.77357 ns [5.76343, 5.78542] (0.95 confidence)</td>
<td>1.00299</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/underflow/difference&#95;endpoints</td>
<td>6.24327 ns [6.2347, 6.25125] (0.95 confidence)</td>
<td>6.21123 ns [6.18498, 6.23765] (0.95 confidence)</td>
<td>1.00516</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/underflow/dot</td>
<td>3.96806 ns [3.95841, 3.97854] (0.95 confidence)</td>
<td>3.95502 ns [3.94939, 3.96059] (0.95 confidence)</td>
<td>1.0033</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/underflow/dot&#95;endpoints</td>
<td>3.97767 ns [3.96938, 3.98593] (0.95 confidence)</td>
<td>3.97017 ns [3.95992, 3.98045] (0.95 confidence)</td>
<td>1.00189</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;batch/origin/construction</td>
<td>2.87478 ns [2.86638, 2.88323] (0.95 confidence)</td>
<td>2.94497 ns [2.92828, 2.96275] (0.95 confidence)</td>
<td>0.976166</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;batch/origin/prepared</td>
<td>47.5431 ns [47.299, 47.8065] (0.95 confidence)</td>
<td>49.2709 ns [48.7566, 49.8107] (0.95 confidence)</td>
<td>0.964934</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;batch/origin/with&#95;construction</td>
<td>59.854 ns [59.5035, 60.2079] (0.95 confidence)</td>
<td>60.8134 ns [60.3642, 61.2786] (0.95 confidence)</td>
<td>0.984224</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;batch/translated/construction</td>
<td>3.28267 ns [3.27697, 3.28834] (0.95 confidence)</td>
<td>3.31657 ns [3.30782, 3.32913] (0.95 confidence)</td>
<td>0.989779</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;batch/translated/prepared</td>
<td>151.434 ns [150.186, 152.729] (0.95 confidence)</td>
<td>148.094 ns [147.814, 148.382] (0.95 confidence)</td>
<td>1.02255</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;batch/translated/with&#95;construction</td>
<td>158.359 ns [157.538, 159.284] (0.95 confidence)</td>
<td>161.127 ns [160.08, 162.258] (0.95 confidence)</td>
<td>0.982821</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/cancellation/difference</td>
<td>12.8103 ns [12.6912, 12.9324] (0.95 confidence)</td>
<td>12.5887 ns [12.4786, 12.7059] (0.95 confidence)</td>
<td>1.0176</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/cancellation/difference&#95;endpoints</td>
<td>13.7146 ns [13.6037, 13.8304] (0.95 confidence)</td>
<td>13.6938 ns [13.6149, 13.775] (0.95 confidence)</td>
<td>1.00152</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/cancellation/dot</td>
<td>6.98034 ns [6.90855, 7.05486] (0.95 confidence)</td>
<td>6.89508 ns [6.84698, 6.94576] (0.95 confidence)</td>
<td>1.01237</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/cancellation/dot&#95;endpoints</td>
<td>7.40052 ns [7.3575, 7.44737] (0.95 confidence)</td>
<td>7.34485 ns [7.31919, 7.37465] (0.95 confidence)</td>
<td>1.00758</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense/difference</td>
<td>23.8169 ns [23.6334, 24.0028] (0.95 confidence)</td>
<td>23.2424 ns [23.1343, 23.3636] (0.95 confidence)</td>
<td>1.02472</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense/difference&#95;endpoints</td>
<td>25.1066 ns [24.925, 25.295] (0.95 confidence)</td>
<td>24.6995 ns [24.5848, 24.8234] (0.95 confidence)</td>
<td>1.01648</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense/dot</td>
<td>12.0207 ns [11.9406, 12.1051] (0.95 confidence)</td>
<td>11.9755 ns [11.8985, 12.0613] (0.95 confidence)</td>
<td>1.00377</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense/dot&#95;endpoints</td>
<td>12.5544 ns [12.4963, 12.6185] (0.95 confidence)</td>
<td>12.6627 ns [12.5876, 12.7434] (0.95 confidence)</td>
<td>0.991446</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense&#95;inexact/difference</td>
<td>23.7841 ns [23.5664, 24.0104] (0.95 confidence)</td>
<td>23.5171 ns [23.3935, 23.6558] (0.95 confidence)</td>
<td>1.01136</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense&#95;inexact/difference&#95;endpoints</td>
<td>26.2449 ns [26.083, 26.4096] (0.95 confidence)</td>
<td>26.1696 ns [26.0423, 26.3037] (0.95 confidence)</td>
<td>1.00288</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense&#95;inexact/dot</td>
<td>12.0347 ns [11.9477, 12.128] (0.95 confidence)</td>
<td>11.794 ns [11.7615, 11.8338] (0.95 confidence)</td>
<td>1.02041</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense&#95;inexact/dot&#95;endpoints</td>
<td>13.6729 ns [13.5914, 13.7578] (0.95 confidence)</td>
<td>13.6006 ns [13.5218, 13.6857] (0.95 confidence)</td>
<td>1.00532</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/overflow/dot</td>
<td>1.00337 ns [1.00218, 1.00461] (0.95 confidence)</td>
<td>1.01338 ns [1.011, 1.01586] (0.95 confidence)</td>
<td>0.990122</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/sparse/difference</td>
<td>7.31984 ns [7.30238, 7.33818] (0.95 confidence)</td>
<td>7.39094 ns [7.35858, 7.42756] (0.95 confidence)</td>
<td>0.99038</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/sparse/difference&#95;endpoints</td>
<td>8.46101 ns [8.44244, 8.48162] (0.95 confidence)</td>
<td>8.50704 ns [8.48715, 8.52852] (0.95 confidence)</td>
<td>0.994589</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/sparse/dot</td>
<td>3.81286 ns [3.7743, 3.85448] (0.95 confidence)</td>
<td>3.69359 ns [3.6777, 3.71098] (0.95 confidence)</td>
<td>1.03229</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/sparse/dot&#95;endpoints</td>
<td>3.79981 ns [3.78317, 3.81869] (0.95 confidence)</td>
<td>3.83199 ns [3.82058, 3.84458] (0.95 confidence)</td>
<td>0.991604</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/underflow/difference</td>
<td>7.78785 ns [7.76907, 7.80752] (0.95 confidence)</td>
<td>7.85671 ns [7.83124, 7.8858] (0.95 confidence)</td>
<td>0.991236</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/underflow/difference&#95;endpoints</td>
<td>8.35431 ns [8.32959, 8.37915] (0.95 confidence)</td>
<td>8.44528 ns [8.41179, 8.48061] (0.95 confidence)</td>
<td>0.989228</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/underflow/dot</td>
<td>4.66525 ns [4.65326, 4.68027] (0.95 confidence)</td>
<td>4.74982 ns [4.73661, 4.76391] (0.95 confidence)</td>
<td>0.982196</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/underflow/dot&#95;endpoints</td>
<td>4.65589 ns [4.64737, 4.66478] (0.95 confidence)</td>
<td>4.71441 ns [4.70729, 4.72137] (0.95 confidence)</td>
<td>0.987588</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;batch/origin/construction</td>
<td>3.12633 ns [3.12138, 3.13242] (0.95 confidence)</td>
<td>3.16258 ns [3.15654, 3.16983] (0.95 confidence)</td>
<td>0.988537</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;batch/origin/prepared</td>
<td>62.7622 ns [62.4142, 63.1287] (0.95 confidence)</td>
<td>63.708 ns [63.3632, 64.06] (0.95 confidence)</td>
<td>0.985153</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;batch/origin/with&#95;construction</td>
<td>82.9628 ns [82.5661, 83.3764] (0.95 confidence)</td>
<td>84.2044 ns [83.8705, 84.5412] (0.95 confidence)</td>
<td>0.985255</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;batch/translated/construction</td>
<td>3.607 ns [3.60321, 3.61081] (0.95 confidence)</td>
<td>3.62526 ns [3.621, 3.6308] (0.95 confidence)</td>
<td>0.994965</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;batch/translated/prepared</td>
<td>209.245 ns [208.446, 210.106] (0.95 confidence)</td>
<td>211.107 ns [210.249, 212.033] (0.95 confidence)</td>
<td>0.991182</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;batch/translated/with&#95;construction</td>
<td>232.275 ns [231.539, 233.079] (0.95 confidence)</td>
<td>235.722 ns [234.6, 237.16] (0.95 confidence)</td>
<td>0.985377</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/cancellation/difference</td>
<td>12.894 ns [12.8178, 12.9774] (0.95 confidence)</td>
<td>12.9319 ns [12.8695, 13.0002] (0.95 confidence)</td>
<td>0.997074</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/cancellation/difference&#95;endpoints</td>
<td>14.4712 ns [14.3574, 14.5897] (0.95 confidence)</td>
<td>14.2142 ns [14.1514, 14.286] (0.95 confidence)</td>
<td>1.01808</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/cancellation/dot</td>
<td>6.98307 ns [6.9449, 7.0271] (0.95 confidence)</td>
<td>6.99078 ns [6.95995, 7.0266] (0.95 confidence)</td>
<td>0.998898</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/cancellation/dot&#95;endpoints</td>
<td>7.82314 ns [7.78698, 7.86178] (0.95 confidence)</td>
<td>7.80408 ns [7.76888, 7.84295] (0.95 confidence)</td>
<td>1.00244</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense/difference</td>
<td>29.3697 ns [29.2956, 29.4633] (0.95 confidence)</td>
<td>29.4886 ns [29.4449, 29.5323] (0.95 confidence)</td>
<td>0.995968</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense/difference&#95;endpoints</td>
<td>32.7869 ns [32.7243, 32.8497] (0.95 confidence)</td>
<td>32.8956 ns [32.8317, 32.9598] (0.95 confidence)</td>
<td>0.996693</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense/dot</td>
<td>14.0751 ns [13.9911, 14.1663] (0.95 confidence)</td>
<td>14.2209 ns [14.1505, 14.2967] (0.95 confidence)</td>
<td>0.989748</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense/dot&#95;endpoints</td>
<td>15.4349 ns [15.3735, 15.5074] (0.95 confidence)</td>
<td>15.7078 ns [15.6204, 15.8037] (0.95 confidence)</td>
<td>0.982629</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense&#95;inexact/difference</td>
<td>29.3661 ns [29.3141, 29.4173] (0.95 confidence)</td>
<td>29.6652 ns [29.6099, 29.7212] (0.95 confidence)</td>
<td>0.989916</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense&#95;inexact/difference&#95;endpoints</td>
<td>32.0627 ns [31.996, 32.1287] (0.95 confidence)</td>
<td>32.2081 ns [32.1518, 32.2618] (0.95 confidence)</td>
<td>0.995485</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense&#95;inexact/dot</td>
<td>14.0411 ns [13.9763, 14.1131] (0.95 confidence)</td>
<td>14.0377 ns [13.9903, 14.0902] (0.95 confidence)</td>
<td>1.00024</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense&#95;inexact/dot&#95;endpoints</td>
<td>15.7155 ns [15.5827, 15.86] (0.95 confidence)</td>
<td>15.9835 ns [15.8534, 16.1179] (0.95 confidence)</td>
<td>0.983234</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/overflow/dot</td>
<td>1.00266 ns [1.00096, 1.00447] (0.95 confidence)</td>
<td>1.01543 ns [1.00924, 1.02547] (0.95 confidence)</td>
<td>0.987427</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/sparse/difference</td>
<td>8.16483 ns [8.15056, 8.18011] (0.95 confidence)</td>
<td>8.21671 ns [8.19927, 8.23573] (0.95 confidence)</td>
<td>0.993685</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/sparse/difference&#95;endpoints</td>
<td>9.44762 ns [9.42316, 9.47249] (0.95 confidence)</td>
<td>9.57896 ns [9.55324, 9.60664] (0.95 confidence)</td>
<td>0.98629</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/sparse/dot</td>
<td>4.12989 ns [4.12502, 4.13487] (0.95 confidence)</td>
<td>4.20746 ns [4.17648, 4.24361] (0.95 confidence)</td>
<td>0.981563</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/sparse/dot&#95;endpoints</td>
<td>4.20093 ns [4.19521, 4.20726] (0.95 confidence)</td>
<td>4.26227 ns [4.25403, 4.27085] (0.95 confidence)</td>
<td>0.985609</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/underflow/difference</td>
<td>8.73037 ns [8.70033, 8.76266] (0.95 confidence)</td>
<td>8.75112 ns [8.7042, 8.79765] (0.95 confidence)</td>
<td>0.997629</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/underflow/difference&#95;endpoints</td>
<td>9.17232 ns [9.14485, 9.19795] (0.95 confidence)</td>
<td>9.25919 ns [9.21136, 9.30663] (0.95 confidence)</td>
<td>0.990618</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/underflow/dot</td>
<td>5.37825 ns [5.36931, 5.38739] (0.95 confidence)</td>
<td>5.51529 ns [5.4889, 5.54715] (0.95 confidence)</td>
<td>0.975153</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/underflow/dot&#95;endpoints</td>
<td>5.44061 ns [5.42536, 5.45373] (0.95 confidence)</td>
<td>5.54183 ns [5.53574, 5.548] (0.95 confidence)</td>
<td>0.981734</td>
</tr>
<tr>
<td>public&#95;d2/det&#95;bound</td>
<td>1.95114 ns [1.94737, 1.9551] (0.95 confidence)</td>
<td>1.94636 ns [1.94426, 1.94855] (0.95 confidence)</td>
<td>1.00246</td>
</tr>
<tr>
<td>public&#95;d2/det&#95;direct</td>
<td>1.14089 ns [1.13958, 1.14223] (0.95 confidence)</td>
<td>1.14428 ns [1.14219, 1.14679] (0.95 confidence)</td>
<td>0.997035</td>
</tr>
<tr>
<td>public&#95;d2/det&#95;with&#95;bound</td>
<td>2.15879 ns [2.15077, 2.16546] (0.95 confidence)</td>
<td>2.13283 ns [2.12484, 2.13955] (0.95 confidence)</td>
<td>1.01217</td>
</tr>
<tr>
<td>public&#95;d2/first&#95;asymmetry</td>
<td>1.39105 ns [1.3892, 1.3931] (0.95 confidence)</td>
<td>1.40094 ns [1.39897, 1.40315] (0.95 confidence)</td>
<td>0.992938</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;det</td>
<td>7.89607 ns [7.84956, 7.94459] (0.95 confidence)</td>
<td>7.94507 ns [7.89366, 7.99926] (0.95 confidence)</td>
<td>0.993833</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;access</td>
<td>6.42205 ns [6.40769, 6.43568] (0.95 confidence)</td>
<td>6.53599 ns [6.51808, 6.55902] (0.95 confidence)</td>
<td>0.982567</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;from&#95;matrix</td>
<td>2.61934 ns [2.61267, 2.62587] (0.95 confidence)</td>
<td>2.65352 ns [2.64783, 2.65906] (0.95 confidence)</td>
<td>0.987119</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;from&#95;points</td>
<td>2.27525 ns [2.26604, 2.28471] (0.95 confidence)</td>
<td>2.29198 ns [2.28282, 2.30144] (0.95 confidence)</td>
<td>0.992699</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;from&#95;rows</td>
<td>3.6149 ns [3.60747, 3.62196] (0.95 confidence)</td>
<td>3.68068 ns [3.67218, 3.69002] (0.95 confidence)</td>
<td>0.98213</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;identity</td>
<td>2.13433 ns [2.12886, 2.13989] (0.95 confidence)</td>
<td>2.17186 ns [2.16694, 2.17672] (0.95 confidence)</td>
<td>0.982718</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;set</td>
<td>3.85296 ns [3.84394, 3.86197] (0.95 confidence)</td>
<td>3.93916 ns [3.93208, 3.94631] (0.95 confidence)</td>
<td>0.978119</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;zero</td>
<td>2.14323 ns [2.1393, 2.14714] (0.95 confidence)</td>
<td>2.19583 ns [2.19051, 2.20133] (0.95 confidence)</td>
<td>0.976047</td>
</tr>
<tr>
<td>public&#95;d2/is&#95;symmetric</td>
<td>0.7639 ns [0.76183, 0.765899] (0.95 confidence)</td>
<td>0.76932 ns [0.767445, 0.771977] (0.95 confidence)</td>
<td>0.992956</td>
</tr>
<tr>
<td>public&#95;d2/matrix&#95;access</td>
<td>4.79356 ns [4.78257, 4.80421] (0.95 confidence)</td>
<td>4.89528 ns [4.88512, 4.90515] (0.95 confidence)</td>
<td>0.97922</td>
</tr>
<tr>
<td>public&#95;d2/matrix&#95;from&#95;rows</td>
<td>2.93643 ns [2.93255, 2.94028] (0.95 confidence)</td>
<td>2.99958 ns [2.99369, 3.00582] (0.95 confidence)</td>
<td>0.978949</td>
</tr>
<tr>
<td>public&#95;d2/matrix&#95;identity</td>
<td>1.2786 ns [1.27511, 1.28271] (0.95 confidence)</td>
<td>1.3047 ns [1.30201, 1.30732] (0.95 confidence)</td>
<td>0.979998</td>
</tr>
<tr>
<td>public&#95;d2/matrix&#95;set</td>
<td>2.13106 ns [2.1218, 2.1403] (0.95 confidence)</td>
<td>2.15921 ns [2.15285, 2.16569] (0.95 confidence)</td>
<td>0.986966</td>
</tr>
<tr>
<td>public&#95;d2/matrix&#95;zero</td>
<td>1.18119 ns [1.17804, 1.18455] (0.95 confidence)</td>
<td>1.20852 ns [1.20573, 1.21135] (0.95 confidence)</td>
<td>0.977379</td>
</tr>
<tr>
<td>public&#95;d2/vector&#95;access</td>
<td>1.64163 ns [1.63784, 1.64507] (0.95 confidence)</td>
<td>1.66233 ns [1.65847, 1.66609] (0.95 confidence)</td>
<td>0.987545</td>
</tr>
<tr>
<td>public&#95;d2/vector&#95;new</td>
<td>1.96022 ns [1.94946, 1.97009] (0.95 confidence)</td>
<td>1.82872 ns [1.80772, 1.85067] (0.95 confidence)</td>
<td>1.07191</td>
</tr>
<tr>
<td>public&#95;d2/vector&#95;zero</td>
<td>0.601984 ns [0.600975, 0.602987] (0.95 confidence)</td>
<td>0.610871 ns [0.609895, 0.611825] (0.95 confidence)</td>
<td>0.985453</td>
</tr>
<tr>
<td>public&#95;d3/det&#95;bound</td>
<td>3.679 ns [3.67335, 3.68546] (0.95 confidence)</td>
<td>3.70726 ns [3.70309, 3.7116] (0.95 confidence)</td>
<td>0.992377</td>
</tr>
<tr>
<td>public&#95;d3/det&#95;direct</td>
<td>1.29698 ns [1.29228, 1.30401] (0.95 confidence)</td>
<td>1.34224 ns [1.32049, 1.36717] (0.95 confidence)</td>
<td>0.966286</td>
</tr>
<tr>
<td>public&#95;d3/det&#95;with&#95;bound</td>
<td>3.6785 ns [3.67302, 3.68451] (0.95 confidence)</td>
<td>3.73279 ns [3.72306, 3.74777] (0.95 confidence)</td>
<td>0.985455</td>
</tr>
<tr>
<td>public&#95;d3/first&#95;asymmetry</td>
<td>3.62336 ns [3.58348, 3.6642] (0.95 confidence)</td>
<td>3.71284 ns [3.66987, 3.75716] (0.95 confidence)</td>
<td>0.975899</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;det</td>
<td>23.1004 ns [22.9416, 23.2444] (0.95 confidence)</td>
<td>22.7312 ns [22.5195, 22.9422] (0.95 confidence)</td>
<td>1.01624</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;access</td>
<td>9.99397 ns [9.97141, 10.0144] (0.95 confidence)</td>
<td>10.3027 ns [10.2848, 10.3205] (0.95 confidence)</td>
<td>0.970033</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;from&#95;matrix</td>
<td>5.11543 ns [5.10345, 5.12654] (0.95 confidence)</td>
<td>5.28587 ns [5.27684, 5.29496] (0.95 confidence)</td>
<td>0.967756</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;from&#95;points</td>
<td>6.316 ns [6.30733, 6.32409] (0.95 confidence)</td>
<td>6.49024 ns [6.47983, 6.50079] (0.95 confidence)</td>
<td>0.973153</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;from&#95;rows</td>
<td>7.32929 ns [7.31313, 7.34507] (0.95 confidence)</td>
<td>7.58448 ns [7.56888, 7.59976] (0.95 confidence)</td>
<td>0.966354</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;identity</td>
<td>4.70753 ns [4.67544, 4.74307] (0.95 confidence)</td>
<td>4.47421 ns [4.46415, 4.4843] (0.95 confidence)</td>
<td>1.05215</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;set</td>
<td>7.83892 ns [7.79821, 7.89393] (0.95 confidence)</td>
<td>8.10395 ns [8.08135, 8.12739] (0.95 confidence)</td>
<td>0.967296</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;zero</td>
<td>4.29144 ns [4.27796, 4.30714] (0.95 confidence)</td>
<td>4.49838 ns [4.48667, 4.51175] (0.95 confidence)</td>
<td>0.953998</td>
</tr>
<tr>
<td>public&#95;d3/is&#95;symmetric</td>
<td>1.34132 ns [1.33949, 1.34311] (0.95 confidence)</td>
<td>1.35089 ns [1.34844, 1.35354] (0.95 confidence)</td>
<td>0.992914</td>
</tr>
<tr>
<td>public&#95;d3/matrix&#95;access</td>
<td>6.76914 ns [6.74316, 6.7958] (0.95 confidence)</td>
<td>6.88154 ns [6.86312, 6.89952] (0.95 confidence)</td>
<td>0.983666</td>
</tr>
<tr>
<td>public&#95;d3/matrix&#95;from&#95;rows</td>
<td>5.80572 ns [5.79118, 5.82132] (0.95 confidence)</td>
<td>5.8589 ns [5.84008, 5.87891] (0.95 confidence)</td>
<td>0.990923</td>
</tr>
<tr>
<td>public&#95;d3/matrix&#95;identity</td>
<td>2.75901 ns [2.74973, 2.76888] (0.95 confidence)</td>
<td>2.79535 ns [2.78829, 2.80294] (0.95 confidence)</td>
<td>0.987001</td>
</tr>
<tr>
<td>public&#95;d3/matrix&#95;set</td>
<td>4.27838 ns [4.26546, 4.29129] (0.95 confidence)</td>
<td>4.35024 ns [4.33453, 4.36775] (0.95 confidence)</td>
<td>0.98348</td>
</tr>
<tr>
<td>public&#95;d3/matrix&#95;zero</td>
<td>2.46107 ns [2.45541, 2.46712] (0.95 confidence)</td>
<td>2.51535 ns [2.50962, 2.52115] (0.95 confidence)</td>
<td>0.978421</td>
</tr>
<tr>
<td>public&#95;d3/vector&#95;access</td>
<td>2.08496 ns [2.08098, 2.08892] (0.95 confidence)</td>
<td>2.15007 ns [2.14597, 2.15426] (0.95 confidence)</td>
<td>0.969718</td>
</tr>
<tr>
<td>public&#95;d3/vector&#95;new</td>
<td>3.65905 ns [3.6358, 3.68573] (0.95 confidence)</td>
<td>3.71035 ns [3.69094, 3.73266] (0.95 confidence)</td>
<td>0.986174</td>
</tr>
<tr>
<td>public&#95;d3/vector&#95;zero</td>
<td>0.95297 ns [0.951182, 0.954739] (0.95 confidence)</td>
<td>0.988988 ns [0.984371, 0.995653] (0.95 confidence)</td>
<td>0.963581</td>
</tr>
<tr>
<td>public&#95;d4/det&#95;bound</td>
<td>6.39577 ns [6.37833, 6.41673] (0.95 confidence)</td>
<td>6.62511 ns [6.60142, 6.64902] (0.95 confidence)</td>
<td>0.965384</td>
</tr>
<tr>
<td>public&#95;d4/det&#95;direct</td>
<td>2.47212 ns [2.4707, 2.47355] (0.95 confidence)</td>
<td>2.50799 ns [2.50502, 2.51116] (0.95 confidence)</td>
<td>0.985698</td>
</tr>
<tr>
<td>public&#95;d4/det&#95;with&#95;bound</td>
<td>6.45419 ns [6.43298, 6.47646] (0.95 confidence)</td>
<td>6.63225 ns [6.60846, 6.6559] (0.95 confidence)</td>
<td>0.973151</td>
</tr>
<tr>
<td>public&#95;d4/first&#95;asymmetry</td>
<td>3.72093 ns [3.68567, 3.75723] (0.95 confidence)</td>
<td>3.78543 ns [3.75266, 3.81922] (0.95 confidence)</td>
<td>0.982961</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;det</td>
<td>49.6771 ns [49.2148, 50.1726] (0.95 confidence)</td>
<td>54.1887 ns [53.7543, 54.6548] (0.95 confidence)</td>
<td>0.916743</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;access</td>
<td>13.5073 ns [13.465, 13.5486] (0.95 confidence)</td>
<td>14.0645 ns [14.0232, 14.1071] (0.95 confidence)</td>
<td>0.960384</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;from&#95;matrix</td>
<td>7.80035 ns [7.78665, 7.8131] (0.95 confidence)</td>
<td>8.18562 ns [8.16921, 8.20234] (0.95 confidence)</td>
<td>0.952934</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;from&#95;points</td>
<td>10.8374 ns [10.8134, 10.8647] (0.95 confidence)</td>
<td>11.1776 ns [11.1475, 11.2097] (0.95 confidence)</td>
<td>0.969568</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;from&#95;rows</td>
<td>11.033 ns [11.0151, 11.0506] (0.95 confidence)</td>
<td>11.6038 ns [11.5818, 11.6257] (0.95 confidence)</td>
<td>0.950811</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;identity</td>
<td>6.42719 ns [6.40675, 6.45287] (0.95 confidence)</td>
<td>6.89123 ns [6.8717, 6.91134] (0.95 confidence)</td>
<td>0.932663</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;set</td>
<td>11.231 ns [11.2143, 11.2465] (0.95 confidence)</td>
<td>11.941 ns [11.9168, 11.966] (0.95 confidence)</td>
<td>0.940547</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;zero</td>
<td>6.48848 ns [6.47335, 6.50579] (0.95 confidence)</td>
<td>7.08644 ns [7.04487, 7.12907] (0.95 confidence)</td>
<td>0.915619</td>
</tr>
<tr>
<td>public&#95;d4/is&#95;symmetric</td>
<td>2.10906 ns [2.10316, 2.11951] (0.95 confidence)</td>
<td>2.13028 ns [2.12808, 2.13258] (0.95 confidence)</td>
<td>0.990039</td>
</tr>
<tr>
<td>public&#95;d4/matrix&#95;access</td>
<td>9.3021 ns [9.28295, 9.32066] (0.95 confidence)</td>
<td>9.73495 ns [9.71045, 9.76019] (0.95 confidence)</td>
<td>0.955537</td>
</tr>
<tr>
<td>public&#95;d4/matrix&#95;from&#95;rows</td>
<td>10.2693 ns [10.2603, 10.278] (0.95 confidence)</td>
<td>10.5457 ns [10.5345, 10.5571] (0.95 confidence)</td>
<td>0.973784</td>
</tr>
<tr>
<td>public&#95;d4/matrix&#95;identity</td>
<td>4.1111 ns [4.1024, 4.11975] (0.95 confidence)</td>
<td>4.28306 ns [4.27347, 4.29294] (0.95 confidence)</td>
<td>0.959851</td>
</tr>
<tr>
<td>public&#95;d4/matrix&#95;set</td>
<td>7.32682 ns [7.31298, 7.3398] (0.95 confidence)</td>
<td>7.65385 ns [7.63594, 7.67208] (0.95 confidence)</td>
<td>0.957272</td>
</tr>
<tr>
<td>public&#95;d4/matrix&#95;zero</td>
<td>3.8433 ns [3.83698, 3.8496] (0.95 confidence)</td>
<td>3.98957 ns [3.97823, 4.00241] (0.95 confidence)</td>
<td>0.963335</td>
</tr>
<tr>
<td>public&#95;d4/vector&#95;access</td>
<td>2.59029 ns [2.5845, 2.59555] (0.95 confidence)</td>
<td>2.69899 ns [2.69399, 2.70386] (0.95 confidence)</td>
<td>0.959728</td>
</tr>
<tr>
<td>public&#95;d4/vector&#95;new</td>
<td>2.7841 ns [2.75818, 2.80706] (0.95 confidence)</td>
<td>2.90478 ns [2.88384, 2.92612] (0.95 confidence)</td>
<td>0.958456</td>
</tr>
<tr>
<td>public&#95;d4/vector&#95;zero</td>
<td>1.15288 ns [1.15113, 1.15459] (0.95 confidence)</td>
<td>1.21304 ns [1.21024, 1.2159] (0.95 confidence)</td>
<td>0.950404</td>
</tr>
<tr>
<td>public&#95;d5/det&#95;bound</td>
<td>1.31074 ns [1.31029, 1.31127] (0.95 confidence)</td>
<td>1.41783 ns [1.41198, 1.42425] (0.95 confidence)</td>
<td>0.924474</td>
</tr>
<tr>
<td>public&#95;d5/det&#95;direct</td>
<td>1.31633 ns [1.31208, 1.32375] (0.95 confidence)</td>
<td>1.41265 ns [1.40597, 1.42044] (0.95 confidence)</td>
<td>0.931817</td>
</tr>
<tr>
<td>public&#95;d5/det&#95;with&#95;bound</td>
<td>0.843058 ns [0.841348, 0.844763] (0.95 confidence)</td>
<td>0.898503 ns [0.896125, 0.900957] (0.95 confidence)</td>
<td>0.938291</td>
</tr>
<tr>
<td>public&#95;d5/first&#95;asymmetry</td>
<td>3.80354 ns [3.7988, 3.80883] (0.95 confidence)</td>
<td>3.85648 ns [3.85046, 3.86277] (0.95 confidence)</td>
<td>0.986273</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;det</td>
<td>124.161 ns [123.459, 124.886] (0.95 confidence)</td>
<td>131.119 ns [130.443, 131.775] (0.95 confidence)</td>
<td>0.946933</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;access</td>
<td>15.4533 ns [15.3499, 15.5715] (0.95 confidence)</td>
<td>15.7427 ns [15.5543, 15.9466] (0.95 confidence)</td>
<td>0.981621</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;from&#95;matrix</td>
<td>11.1501 ns [11.1095, 11.1915] (0.95 confidence)</td>
<td>10.8727 ns [10.8387, 10.9069] (0.95 confidence)</td>
<td>1.02552</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;from&#95;points</td>
<td>19.3303 ns [19.3118, 19.3493] (0.95 confidence)</td>
<td>19.0386 ns [19.0093, 19.0697] (0.95 confidence)</td>
<td>1.01532</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;from&#95;rows</td>
<td>15.8279 ns [15.6511, 15.9839] (0.95 confidence)</td>
<td>15.2415 ns [15.0168, 15.4586] (0.95 confidence)</td>
<td>1.03847</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;identity</td>
<td>8.73271 ns [8.65785, 8.80761] (0.95 confidence)</td>
<td>8.18986 ns [8.17588, 8.20553] (0.95 confidence)</td>
<td>1.06628</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;set</td>
<td>13.9892 ns [13.7545, 14.2294] (0.95 confidence)</td>
<td>14.0363 ns [13.8192, 14.2598] (0.95 confidence)</td>
<td>0.996642</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;zero</td>
<td>9.8935 ns [9.87967, 9.90797] (0.95 confidence)</td>
<td>9.65042 ns [9.63488, 9.66589] (0.95 confidence)</td>
<td>1.02519</td>
</tr>
<tr>
<td>public&#95;d5/is&#95;symmetric</td>
<td>3.68602 ns [3.68358, 3.68937] (0.95 confidence)</td>
<td>3.73626 ns [3.73186, 3.74116] (0.95 confidence)</td>
<td>0.986555</td>
</tr>
<tr>
<td>public&#95;d5/matrix&#95;access</td>
<td>12.0422 ns [11.9952, 12.0885] (0.95 confidence)</td>
<td>12.5899 ns [12.5463, 12.6346] (0.95 confidence)</td>
<td>0.956499</td>
</tr>
<tr>
<td>public&#95;d5/matrix&#95;from&#95;rows</td>
<td>11.7666 ns [11.7395, 11.7957] (0.95 confidence)</td>
<td>12.0375 ns [12.014, 12.0649] (0.95 confidence)</td>
<td>0.977502</td>
</tr>
<tr>
<td>public&#95;d5/matrix&#95;identity</td>
<td>6.6445 ns [6.63047, 6.65845] (0.95 confidence)</td>
<td>7.11644 ns [7.0927, 7.14183] (0.95 confidence)</td>
<td>0.933684</td>
</tr>
<tr>
<td>public&#95;d5/matrix&#95;set</td>
<td>10.3085 ns [10.2853, 10.3309] (0.95 confidence)</td>
<td>10.9044 ns [10.8862, 10.9223] (0.95 confidence)</td>
<td>0.945354</td>
</tr>
<tr>
<td>public&#95;d5/matrix&#95;zero</td>
<td>5.61284 ns [5.60384, 5.62192] (0.95 confidence)</td>
<td>6.01221 ns [5.99553, 6.02999] (0.95 confidence)</td>
<td>0.933573</td>
</tr>
<tr>
<td>public&#95;d5/vector&#95;access</td>
<td>3.03988 ns [3.03232, 3.04731] (0.95 confidence)</td>
<td>3.19712 ns [3.1899, 3.20445] (0.95 confidence)</td>
<td>0.950819</td>
</tr>
<tr>
<td>public&#95;d5/vector&#95;new</td>
<td>4.82633 ns [4.81152, 4.83904] (0.95 confidence)</td>
<td>4.91636 ns [4.90468, 4.92805] (0.95 confidence)</td>
<td>0.981687</td>
</tr>
<tr>
<td>public&#95;d5/vector&#95;zero</td>
<td>1.47931 ns [1.46986, 1.48837] (0.95 confidence)</td>
<td>1.61404 ns [1.60454, 1.62355] (0.95 confidence)</td>
<td>0.916522</td>
</tr>
<tr>
<td>public&#95;scalar/interval&#95;access</td>
<td>1.60978 ns [1.60446, 1.61524] (0.95 confidence)</td>
<td>1.6147 ns [1.61115, 1.61789] (0.95 confidence)</td>
<td>0.996953</td>
</tr>
<tr>
<td>public&#95;scalar/interval&#95;negate</td>
<td>0.737897 ns [0.737179, 0.738586] (0.95 confidence)</td>
<td>0.744525 ns [0.743047, 0.746207] (0.95 confidence)</td>
<td>0.991098</td>
</tr>
<tr>
<td>public&#95;scalar/interval&#95;new</td>
<td>0.896676 ns [0.895427, 0.898077] (0.95 confidence)</td>
<td>0.900118 ns [0.899492, 0.900745] (0.95 confidence)</td>
<td>0.996176</td>
</tr>
<tr>
<td>public&#95;scalar/interval&#95;point</td>
<td>0.59006 ns [0.588009, 0.592412] (0.95 confidence)</td>
<td>0.588506 ns [0.588026, 0.589023] (0.95 confidence)</td>
<td>1.00264</td>
</tr>
<tr>
<td>public&#95;scalar/interval&#95;subtraction</td>
<td>1.89251 ns [1.88755, 1.89775] (0.95 confidence)</td>
<td>1.90154 ns [1.89624, 1.90987] (0.95 confidence)</td>
<td>0.995256</td>
</tr>
<tr>
<td>public&#95;scalar/tolerance</td>
<td>0.989852 ns [0.988303, 0.991335] (0.95 confidence)</td>
<td>1.00463 ns [1.00301, 1.00626] (0.95 confidence)</td>
<td>0.985291</td>
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
<td>23.6383 ns [23.5924, 23.6712] (0.95 confidence)</td>
<td>24.5487 ns [24.4859, 24.6166] (0.95 confidence)</td>
<td>0.962916</td>
</tr>
<tr>
<td>angle&#95;d3/dense/construct&#95;then&#95;angle</td>
<td>22.7139 ns [22.695, 22.7373] (0.95 confidence)</td>
<td>23.6371 ns [23.6014, 23.6905] (0.95 confidence)</td>
<td>0.960941</td>
</tr>
<tr>
<td>angle&#95;d3/dense/prepared&#95;vector</td>
<td>22.4223 ns [22.4013, 22.4584] (0.95 confidence)</td>
<td>23.5494 ns [23.474, 23.7458] (0.95 confidence)</td>
<td>0.952137</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;antipodal/borrowed&#95;slice</td>
<td>16.7019 ns [16.6794, 16.7096] (0.95 confidence)</td>
<td>22.2659 ns [21.2658, 22.8578] (0.95 confidence)</td>
<td>0.750109</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;antipodal/construct&#95;then&#95;angle</td>
<td>15.562 ns [15.5536, 15.596] (0.95 confidence)</td>
<td>20.0299 ns [19.0835, 20.7088] (0.95 confidence)</td>
<td>0.776941</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;antipodal/prepared&#95;vector</td>
<td>15.6408 ns [15.6163, 15.6509] (0.95 confidence)</td>
<td>20.8114 ns [20.2269, 21.6587] (0.95 confidence)</td>
<td>0.751548</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;identical/borrowed&#95;slice</td>
<td>14.2009 ns [14.1856, 14.2061] (0.95 confidence)</td>
<td>17.8948 ns [17.3898, 18.4653] (0.95 confidence)</td>
<td>0.79358</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;identical/construct&#95;then&#95;angle</td>
<td>13.4729 ns [13.4643, 13.4765] (0.95 confidence)</td>
<td>18.2128 ns [17.7298, 18.519] (0.95 confidence)</td>
<td>0.739748</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;identical/prepared&#95;vector</td>
<td>13.3278 ns [13.3169, 13.3432] (0.95 confidence)</td>
<td>13.8909 ns [13.8079, 13.938] (0.95 confidence)</td>
<td>0.959462</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;varied/borrowed&#95;slice</td>
<td>25.782 ns [25.7172, 25.8216] (0.95 confidence)</td>
<td>26.6881 ns [26.6852, 26.6919] (0.95 confidence)</td>
<td>0.966047</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;varied/construct&#95;then&#95;angle</td>
<td>24.8311 ns [24.8043, 24.862] (0.95 confidence)</td>
<td>25.6874 ns [25.6541, 25.7577] (0.95 confidence)</td>
<td>0.966664</td>
</tr>
<tr>
<td>angle&#95;d3/dense&#95;varied/prepared&#95;vector</td>
<td>24.4303 ns [24.3782, 24.4709] (0.95 confidence)</td>
<td>25.3098 ns [25.2989, 25.318] (0.95 confidence)</td>
<td>0.965253</td>
</tr>
<tr>
<td>angle&#95;d3/maximum&#95;small&#95;angle/borrowed&#95;slice</td>
<td>11.6307 ns [11.6214, 11.6414] (0.95 confidence)</td>
<td>11.6541 ns [11.6438, 11.6603] (0.95 confidence)</td>
<td>0.99799</td>
</tr>
<tr>
<td>angle&#95;d3/maximum&#95;small&#95;angle/construct&#95;then&#95;angle</td>
<td>10.4296 ns [10.4247, 10.4367] (0.95 confidence)</td>
<td>10.3577 ns [10.351, 10.3661] (0.95 confidence)</td>
<td>1.00694</td>
</tr>
<tr>
<td>angle&#95;d3/maximum&#95;small&#95;angle/prepared&#95;vector</td>
<td>9.73881 ns [9.72348, 9.74609] (0.95 confidence)</td>
<td>9.7401 ns [9.7269, 9.76801] (0.95 confidence)</td>
<td>0.999868</td>
</tr>
<tr>
<td>angle&#95;d3/mixed&#95;scale/borrowed&#95;slice</td>
<td>12.9742 ns [12.9688, 12.9786] (0.95 confidence)</td>
<td>13.0603 ns [13.0535, 13.0683] (0.95 confidence)</td>
<td>0.993411</td>
</tr>
<tr>
<td>angle&#95;d3/mixed&#95;scale/construct&#95;then&#95;angle</td>
<td>11.6123 ns [11.6025, 11.6685] (0.95 confidence)</td>
<td>11.5119 ns [11.4988, 11.552] (0.95 confidence)</td>
<td>1.00872</td>
</tr>
<tr>
<td>angle&#95;d3/mixed&#95;scale/prepared&#95;vector</td>
<td>11.0269 ns [11.013, 11.0577] (0.95 confidence)</td>
<td>11.0057 ns [10.9911, 11.0332] (0.95 confidence)</td>
<td>1.00192</td>
</tr>
<tr>
<td>angle&#95;d3/near&#95;antipodal/borrowed&#95;slice</td>
<td>17.2169 ns [17.1951, 17.2411] (0.95 confidence)</td>
<td>17.1292 ns [17.0975, 17.1357] (0.95 confidence)</td>
<td>1.00512</td>
</tr>
<tr>
<td>angle&#95;d3/near&#95;antipodal/construct&#95;then&#95;angle</td>
<td>15.8133 ns [15.8043, 15.8185] (0.95 confidence)</td>
<td>15.6385 ns [15.6179, 15.6516] (0.95 confidence)</td>
<td>1.01118</td>
</tr>
<tr>
<td>angle&#95;d3/near&#95;antipodal/prepared&#95;vector</td>
<td>15.284 ns [15.2766, 15.3043] (0.95 confidence)</td>
<td>15.2024 ns [15.1745, 15.3165] (0.95 confidence)</td>
<td>1.00537</td>
</tr>
<tr>
<td>angle&#95;d3/near&#95;parallel/borrowed&#95;slice</td>
<td>11.6526 ns [11.6455, 11.662] (0.95 confidence)</td>
<td>11.6465 ns [11.6224, 11.6561] (0.95 confidence)</td>
<td>1.00052</td>
</tr>
<tr>
<td>angle&#95;d3/near&#95;parallel/construct&#95;then&#95;angle</td>
<td>10.4304 ns [10.4201, 10.4411] (0.95 confidence)</td>
<td>10.3548 ns [10.3459, 10.3611] (0.95 confidence)</td>
<td>1.0073</td>
</tr>
<tr>
<td>angle&#95;d3/near&#95;parallel/prepared&#95;vector</td>
<td>9.72887 ns [9.71883, 9.75006] (0.95 confidence)</td>
<td>9.73257 ns [9.72298, 9.74892] (0.95 confidence)</td>
<td>0.99962</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;adjacent/borrowed&#95;slice</td>
<td>13.97 ns [13.9554, 14.0086] (0.95 confidence)</td>
<td>15.0813 ns [15.0256, 15.2076] (0.95 confidence)</td>
<td>0.926315</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;adjacent/construct&#95;then&#95;angle</td>
<td>12.9663 ns [12.9329, 13.0042] (0.95 confidence)</td>
<td>14.0259 ns [13.6343, 14.1474] (0.95 confidence)</td>
<td>0.924457</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;adjacent/prepared&#95;vector</td>
<td>12.064 ns [12.0455, 12.0805] (0.95 confidence)</td>
<td>15.4574 ns [15.0732, 15.8618] (0.95 confidence)</td>
<td>0.780468</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;alternating/borrowed&#95;slice</td>
<td>14.146 ns [14.1249, 14.1993] (0.95 confidence)</td>
<td>15.482 ns [15.2192, 15.6073] (0.95 confidence)</td>
<td>0.913704</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;alternating/construct&#95;then&#95;angle</td>
<td>13.2396 ns [13.2164, 13.343] (0.95 confidence)</td>
<td>13.6231 ns [13.4964, 13.9059] (0.95 confidence)</td>
<td>0.971843</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;alternating/prepared&#95;vector</td>
<td>12.384 ns [12.3637, 12.3877] (0.95 confidence)</td>
<td>16.219 ns [15.7528, 16.923] (0.95 confidence)</td>
<td>0.763549</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;antipodal/borrowed&#95;slice</td>
<td>10.8864 ns [10.8833, 10.9002] (0.95 confidence)</td>
<td>10.8864 ns [10.8851, 10.8877] (0.95 confidence)</td>
<td>1</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;antipodal/construct&#95;then&#95;angle</td>
<td>10.1338 ns [10.1165, 10.1577] (0.95 confidence)</td>
<td>10.1482 ns [10.1431, 10.204] (0.95 confidence)</td>
<td>0.998584</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;antipodal/prepared&#95;vector</td>
<td>8.90137 ns [8.88643, 8.90553] (0.95 confidence)</td>
<td>8.90708 ns [8.90652, 8.90779] (0.95 confidence)</td>
<td>0.999359</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;identical/borrowed&#95;slice</td>
<td>8.65831 ns [8.65799, 8.6598] (0.95 confidence)</td>
<td>8.77244 ns [8.76809, 8.80723] (0.95 confidence)</td>
<td>0.986991</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;identical/construct&#95;then&#95;angle</td>
<td>8.08594 ns [7.99101, 8.16512] (0.95 confidence)</td>
<td>8.39363 ns [8.35997, 8.48883] (0.95 confidence)</td>
<td>0.963342</td>
</tr>
<tr>
<td>angle&#95;d3/sparse&#95;identical/prepared&#95;vector</td>
<td>6.94734 ns [6.94673, 6.94811] (0.95 confidence)</td>
<td>8.69779 ns [8.55146, 8.99391] (0.95 confidence)</td>
<td>0.798748</td>
</tr>
<tr>
<td>angle&#95;d3/subnormal/borrowed&#95;slice</td>
<td>11.6516 ns [11.6439, 11.6653] (0.95 confidence)</td>
<td>11.6573 ns [11.6504, 11.6669] (0.95 confidence)</td>
<td>0.999506</td>
</tr>
<tr>
<td>angle&#95;d3/subnormal/construct&#95;then&#95;angle</td>
<td>10.4385 ns [10.4341, 10.4412] (0.95 confidence)</td>
<td>10.3639 ns [10.3555, 10.368] (0.95 confidence)</td>
<td>1.00719</td>
</tr>
<tr>
<td>angle&#95;d3/subnormal/prepared&#95;vector</td>
<td>9.74016 ns [9.72962, 9.75431] (0.95 confidence)</td>
<td>9.73624 ns [9.72353, 9.75722] (0.95 confidence)</td>
<td>1.0004</td>
</tr>
<tr>
<td>angle&#95;d3/subnormal&#95;axes/borrowed&#95;slice</td>
<td>14.1889 ns [14.1815, 14.199] (0.95 confidence)</td>
<td>14.2211 ns [14.2068, 14.2313] (0.95 confidence)</td>
<td>0.997737</td>
</tr>
<tr>
<td>angle&#95;d3/subnormal&#95;axes/construct&#95;then&#95;angle</td>
<td>13.1711 ns [13.1701, 13.1761] (0.95 confidence)</td>
<td>13.2356 ns [13.1301, 13.2977] (0.95 confidence)</td>
<td>0.995129</td>
</tr>
<tr>
<td>angle&#95;d3/subnormal&#95;axes/prepared&#95;vector</td>
<td>12.367 ns [12.3011, 12.3936] (0.95 confidence)</td>
<td>12.3908 ns [12.3706, 12.4128] (0.95 confidence)</td>
<td>0.998083</td>
</tr>
<tr>
<td>angle&#95;d4/dense/borrowed&#95;slice</td>
<td>33.5586 ns [33.5001, 33.6265] (0.95 confidence)</td>
<td>33.463 ns [33.4107, 33.5007] (0.95 confidence)</td>
<td>1.00286</td>
</tr>
<tr>
<td>angle&#95;d4/dense/construct&#95;then&#95;angle</td>
<td>32.4337 ns [32.3808, 32.4522] (0.95 confidence)</td>
<td>32.2133 ns [32.1787, 32.2635] (0.95 confidence)</td>
<td>1.00684</td>
</tr>
<tr>
<td>angle&#95;d4/dense/prepared&#95;vector</td>
<td>32.002 ns [31.9678, 32.0382] (0.95 confidence)</td>
<td>32.0777 ns [32.0459, 32.0952] (0.95 confidence)</td>
<td>0.997641</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;antipodal/borrowed&#95;slice</td>
<td>26.2669 ns [26.239, 26.2973] (0.95 confidence)</td>
<td>26.1585 ns [26.1326, 26.1793] (0.95 confidence)</td>
<td>1.00414</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;antipodal/construct&#95;then&#95;angle</td>
<td>24.8586 ns [24.8246, 24.8904] (0.95 confidence)</td>
<td>24.7479 ns [24.7165, 24.8131] (0.95 confidence)</td>
<td>1.00447</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;antipodal/prepared&#95;vector</td>
<td>24.4317 ns [24.3943, 24.451] (0.95 confidence)</td>
<td>24.4684 ns [24.443, 24.493] (0.95 confidence)</td>
<td>0.998502</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;identical/borrowed&#95;slice</td>
<td>23.9218 ns [23.9013, 23.9511] (0.95 confidence)</td>
<td>23.8044 ns [23.7779, 23.8247] (0.95 confidence)</td>
<td>1.00493</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;identical/construct&#95;then&#95;angle</td>
<td>22.409 ns [22.3732, 22.44] (0.95 confidence)</td>
<td>22.4403 ns [22.4217, 22.4887] (0.95 confidence)</td>
<td>0.998608</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;identical/prepared&#95;vector</td>
<td>21.877 ns [21.8398, 21.9039] (0.95 confidence)</td>
<td>21.8546 ns [21.8212, 21.8865] (0.95 confidence)</td>
<td>1.00102</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;varied/borrowed&#95;slice</td>
<td>34.9559 ns [34.9262, 34.984] (0.95 confidence)</td>
<td>35.0145 ns [34.9519, 35.0796] (0.95 confidence)</td>
<td>0.998327</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;varied/construct&#95;then&#95;angle</td>
<td>33.7585 ns [33.6581, 33.7935] (0.95 confidence)</td>
<td>33.8205 ns [33.7395, 33.8571] (0.95 confidence)</td>
<td>0.998165</td>
</tr>
<tr>
<td>angle&#95;d4/dense&#95;varied/prepared&#95;vector</td>
<td>32.8359 ns [32.8016, 32.8788] (0.95 confidence)</td>
<td>32.7997 ns [32.746, 32.8811] (0.95 confidence)</td>
<td>1.0011</td>
</tr>
<tr>
<td>angle&#95;d4/maximum&#95;small&#95;angle/borrowed&#95;slice</td>
<td>12.9654 ns [12.9491, 12.9966] (0.95 confidence)</td>
<td>12.6829 ns [12.6525, 12.7195] (0.95 confidence)</td>
<td>1.02227</td>
</tr>
<tr>
<td>angle&#95;d4/maximum&#95;small&#95;angle/construct&#95;then&#95;angle</td>
<td>11.6601 ns [11.5545, 11.6809] (0.95 confidence)</td>
<td>11.1465 ns [11.1387, 11.1543] (0.95 confidence)</td>
<td>1.04608</td>
</tr>
<tr>
<td>angle&#95;d4/maximum&#95;small&#95;angle/prepared&#95;vector</td>
<td>10.5619 ns [10.5292, 10.5809] (0.95 confidence)</td>
<td>10.2926 ns [10.2731, 10.3186] (0.95 confidence)</td>
<td>1.02616</td>
</tr>
<tr>
<td>angle&#95;d4/mixed&#95;scale/borrowed&#95;slice</td>
<td>14.137 ns [14.13, 14.1471] (0.95 confidence)</td>
<td>14.1382 ns [14.1281, 14.1548] (0.95 confidence)</td>
<td>0.999916</td>
</tr>
<tr>
<td>angle&#95;d4/mixed&#95;scale/construct&#95;then&#95;angle</td>
<td>12.6242 ns [12.6188, 12.681] (0.95 confidence)</td>
<td>12.6011 ns [12.5012, 12.618] (0.95 confidence)</td>
<td>1.00183</td>
</tr>
<tr>
<td>angle&#95;d4/mixed&#95;scale/prepared&#95;vector</td>
<td>11.7643 ns [11.727, 11.7916] (0.95 confidence)</td>
<td>11.7183 ns [11.6874, 11.7659] (0.95 confidence)</td>
<td>1.00392</td>
</tr>
<tr>
<td>angle&#95;d4/near&#95;antipodal/borrowed&#95;slice</td>
<td>18.1395 ns [18.136, 18.1451] (0.95 confidence)</td>
<td>18.1286 ns [18.0962, 18.1464] (0.95 confidence)</td>
<td>1.0006</td>
</tr>
<tr>
<td>angle&#95;d4/near&#95;antipodal/construct&#95;then&#95;angle</td>
<td>16.5505 ns [16.5397, 16.5599] (0.95 confidence)</td>
<td>16.5082 ns [16.4993, 16.5272] (0.95 confidence)</td>
<td>1.00256</td>
</tr>
<tr>
<td>angle&#95;d4/near&#95;antipodal/prepared&#95;vector</td>
<td>16.0499 ns [16.037, 16.0666] (0.95 confidence)</td>
<td>16.0024 ns [15.987, 16.0189] (0.95 confidence)</td>
<td>1.00297</td>
</tr>
<tr>
<td>angle&#95;d4/near&#95;parallel/borrowed&#95;slice</td>
<td>12.9203 ns [12.878, 12.9671] (0.95 confidence)</td>
<td>12.6422 ns [12.6177, 12.6584] (0.95 confidence)</td>
<td>1.022</td>
</tr>
<tr>
<td>angle&#95;d4/near&#95;parallel/construct&#95;then&#95;angle</td>
<td>11.1772 ns [11.1519, 11.1905] (0.95 confidence)</td>
<td>11.1436 ns [11.1217, 11.1497] (0.95 confidence)</td>
<td>1.00301</td>
</tr>
<tr>
<td>angle&#95;d4/near&#95;parallel/prepared&#95;vector</td>
<td>10.5491 ns [10.5368, 10.5758] (0.95 confidence)</td>
<td>10.3023 ns [10.2866, 10.3243] (0.95 confidence)</td>
<td>1.02395</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;adjacent/borrowed&#95;slice</td>
<td>15.8114 ns [15.6223, 15.8672] (0.95 confidence)</td>
<td>15.8584 ns [15.7675, 15.8772] (0.95 confidence)</td>
<td>0.997033</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;adjacent/construct&#95;then&#95;angle</td>
<td>14.6861 ns [14.6584, 14.8016] (0.95 confidence)</td>
<td>14.6416 ns [14.6342, 14.6494] (0.95 confidence)</td>
<td>1.00304</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;adjacent/prepared&#95;vector</td>
<td>13.3778 ns [13.3696, 13.3818] (0.95 confidence)</td>
<td>13.3816 ns [13.3724, 13.3839] (0.95 confidence)</td>
<td>0.999711</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;alternating/borrowed&#95;slice</td>
<td>15.9099 ns [15.7964, 16.0527] (0.95 confidence)</td>
<td>16.0448 ns [15.9843, 16.0813] (0.95 confidence)</td>
<td>0.99159</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;alternating/construct&#95;then&#95;angle</td>
<td>14.8663 ns [14.8647, 14.8708] (0.95 confidence)</td>
<td>14.8803 ns [14.8656, 14.9153] (0.95 confidence)</td>
<td>0.999058</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;alternating/prepared&#95;vector</td>
<td>13.6292 ns [13.6284, 13.63] (0.95 confidence)</td>
<td>13.6297 ns [13.6291, 13.6305] (0.95 confidence)</td>
<td>0.999963</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;antipodal/borrowed&#95;slice</td>
<td>11.9234 ns [11.8892, 12.1085] (0.95 confidence)</td>
<td>12.0254 ns [11.9381, 12.1379] (0.95 confidence)</td>
<td>0.991524</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;antipodal/construct&#95;then&#95;angle</td>
<td>10.8855 ns [10.8848, 10.8946] (0.95 confidence)</td>
<td>10.9762 ns [10.9178, 11.0371] (0.95 confidence)</td>
<td>0.991742</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;antipodal/prepared&#95;vector</td>
<td>9.58785 ns [9.57724, 9.61198] (0.95 confidence)</td>
<td>9.61445 ns [9.60182, 9.6243] (0.95 confidence)</td>
<td>0.997234</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;identical/borrowed&#95;slice</td>
<td>10.3666 ns [10.298, 10.3894] (0.95 confidence)</td>
<td>10.3764 ns [10.2752, 10.3902] (0.95 confidence)</td>
<td>0.999054</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;identical/construct&#95;then&#95;angle</td>
<td>9.14562 ns [9.09658, 9.15555] (0.95 confidence)</td>
<td>9.16331 ns [9.15385, 9.18895] (0.95 confidence)</td>
<td>0.99807</td>
</tr>
<tr>
<td>angle&#95;d4/sparse&#95;identical/prepared&#95;vector</td>
<td>7.67004 ns [7.66922, 7.67255] (0.95 confidence)</td>
<td>7.67361 ns [7.67021, 7.6811] (0.95 confidence)</td>
<td>0.999534</td>
</tr>
<tr>
<td>angle&#95;d4/subnormal/borrowed&#95;slice</td>
<td>12.6693 ns [12.6404, 12.7131] (0.95 confidence)</td>
<td>12.6453 ns [12.6233, 12.6669] (0.95 confidence)</td>
<td>1.0019</td>
</tr>
<tr>
<td>angle&#95;d4/subnormal/construct&#95;then&#95;angle</td>
<td>11.1436 ns [11.0939, 11.1478] (0.95 confidence)</td>
<td>11.1639 ns [11.1353, 11.2013] (0.95 confidence)</td>
<td>0.998177</td>
</tr>
<tr>
<td>angle&#95;d4/subnormal/prepared&#95;vector</td>
<td>10.323 ns [10.2922, 10.3396] (0.95 confidence)</td>
<td>10.3018 ns [10.2895, 10.3258] (0.95 confidence)</td>
<td>1.00206</td>
</tr>
<tr>
<td>angle&#95;d4/subnormal&#95;axes/borrowed&#95;slice</td>
<td>16.166 ns [16.0806, 16.2749] (0.95 confidence)</td>
<td>15.6347 ns [15.6327, 15.6387] (0.95 confidence)</td>
<td>1.03398</td>
</tr>
<tr>
<td>angle&#95;d4/subnormal&#95;axes/construct&#95;then&#95;angle</td>
<td>15.2507 ns [15.1415, 15.5037] (0.95 confidence)</td>
<td>14.3699 ns [14.3487, 14.4606] (0.95 confidence)</td>
<td>1.0613</td>
</tr>
<tr>
<td>angle&#95;d4/subnormal&#95;axes/prepared&#95;vector</td>
<td>13.1372 ns [13.1147, 13.3741] (0.95 confidence)</td>
<td>12.997 ns [12.9527, 13.0488] (0.95 confidence)</td>
<td>1.01078</td>
</tr>
<tr>
<td>angle&#95;d5/dense/borrowed&#95;slice</td>
<td>47.7813 ns [47.6963, 47.8622] (0.95 confidence)</td>
<td>47.1131 ns [46.9868, 47.2291] (0.95 confidence)</td>
<td>1.01418</td>
</tr>
<tr>
<td>angle&#95;d5/dense/construct&#95;then&#95;angle</td>
<td>45.3672 ns [45.3137, 45.4138] (0.95 confidence)</td>
<td>45.0296 ns [44.9652, 45.0524] (0.95 confidence)</td>
<td>1.0075</td>
</tr>
<tr>
<td>angle&#95;d5/dense/prepared&#95;vector</td>
<td>44.9896 ns [44.894, 45.058] (0.95 confidence)</td>
<td>44.1956 ns [44.0867, 44.2851] (0.95 confidence)</td>
<td>1.01796</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;antipodal/borrowed&#95;slice</td>
<td>36.1755 ns [36.1323, 36.2252] (0.95 confidence)</td>
<td>36.0191 ns [35.9719, 36.0393] (0.95 confidence)</td>
<td>1.00434</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;antipodal/construct&#95;then&#95;angle</td>
<td>34.2856 ns [34.2358, 34.321] (0.95 confidence)</td>
<td>34.2359 ns [34.191, 34.278] (0.95 confidence)</td>
<td>1.00145</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;antipodal/prepared&#95;vector</td>
<td>34.2057 ns [34.1838, 34.2629] (0.95 confidence)</td>
<td>33.9016 ns [33.8726, 33.9773] (0.95 confidence)</td>
<td>1.00897</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;identical/borrowed&#95;slice</td>
<td>34.9172 ns [34.8518, 34.9805] (0.95 confidence)</td>
<td>34.8124 ns [34.7823, 34.8285] (0.95 confidence)</td>
<td>1.00301</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;identical/construct&#95;then&#95;angle</td>
<td>33.1376 ns [33.0914, 33.1853] (0.95 confidence)</td>
<td>32.9014 ns [32.8819, 32.9481] (0.95 confidence)</td>
<td>1.00718</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;identical/prepared&#95;vector</td>
<td>32.7586 ns [32.6976, 32.81] (0.95 confidence)</td>
<td>32.6232 ns [32.5587, 32.6621] (0.95 confidence)</td>
<td>1.00415</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;varied/borrowed&#95;slice</td>
<td>68.4538 ns [68.2148, 68.5417] (0.95 confidence)</td>
<td>67.0815 ns [66.8824, 67.2573] (0.95 confidence)</td>
<td>1.02046</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;varied/construct&#95;then&#95;angle</td>
<td>65.8836 ns [65.7725, 66.0297] (0.95 confidence)</td>
<td>64.7305 ns [64.5203, 64.8077] (0.95 confidence)</td>
<td>1.01781</td>
</tr>
<tr>
<td>angle&#95;d5/dense&#95;varied/prepared&#95;vector</td>
<td>64.0366 ns [63.8333, 64.2415] (0.95 confidence)</td>
<td>62.977 ns [62.8211, 63.0272] (0.95 confidence)</td>
<td>1.01683</td>
</tr>
<tr>
<td>angle&#95;d5/maximum&#95;small&#95;angle/borrowed&#95;slice</td>
<td>14.3239 ns [14.3129, 14.3473] (0.95 confidence)</td>
<td>14.3233 ns [14.2826, 14.3474] (0.95 confidence)</td>
<td>1.00004</td>
</tr>
<tr>
<td>angle&#95;d5/maximum&#95;small&#95;angle/construct&#95;then&#95;angle</td>
<td>12.6399 ns [12.5223, 12.6449] (0.95 confidence)</td>
<td>12.3979 ns [12.3972, 12.4038] (0.95 confidence)</td>
<td>1.01952</td>
</tr>
<tr>
<td>angle&#95;d5/maximum&#95;small&#95;angle/prepared&#95;vector</td>
<td>11.1826 ns [11.1806, 11.1885] (0.95 confidence)</td>
<td>11.1936 ns [11.1885, 11.198] (0.95 confidence)</td>
<td>0.999022</td>
</tr>
<tr>
<td>angle&#95;d5/mixed&#95;scale/borrowed&#95;slice</td>
<td>15.8357 ns [15.8277, 15.8378] (0.95 confidence)</td>
<td>15.808 ns [15.7904, 15.837] (0.95 confidence)</td>
<td>1.00175</td>
</tr>
<tr>
<td>angle&#95;d5/mixed&#95;scale/construct&#95;then&#95;angle</td>
<td>14.1781 ns [14.1566, 14.3454] (0.95 confidence)</td>
<td>14.1091 ns [14.029, 14.1293] (0.95 confidence)</td>
<td>1.00489</td>
</tr>
<tr>
<td>angle&#95;d5/mixed&#95;scale/prepared&#95;vector</td>
<td>12.6806 ns [12.673, 12.6905] (0.95 confidence)</td>
<td>12.6908 ns [12.6881, 12.6941] (0.95 confidence)</td>
<td>0.999197</td>
</tr>
<tr>
<td>angle&#95;d5/near&#95;antipodal/borrowed&#95;slice</td>
<td>19.8835 ns [19.8521, 19.9109] (0.95 confidence)</td>
<td>19.8175 ns [19.7875, 19.8534] (0.95 confidence)</td>
<td>1.00333</td>
</tr>
<tr>
<td>angle&#95;d5/near&#95;antipodal/construct&#95;then&#95;angle</td>
<td>17.6874 ns [17.679, 17.6918] (0.95 confidence)</td>
<td>17.6841 ns [17.6764, 17.6959] (0.95 confidence)</td>
<td>1.00019</td>
</tr>
<tr>
<td>angle&#95;d5/near&#95;antipodal/prepared&#95;vector</td>
<td>16.7835 ns [16.7241, 16.8688] (0.95 confidence)</td>
<td>16.7835 ns [16.7274, 16.8127] (0.95 confidence)</td>
<td>1.00001</td>
</tr>
<tr>
<td>angle&#95;d5/near&#95;parallel/borrowed&#95;slice</td>
<td>14.3138 ns [14.2997, 14.3342] (0.95 confidence)</td>
<td>14.3387 ns [14.3099, 14.3478] (0.95 confidence)</td>
<td>0.998262</td>
</tr>
<tr>
<td>angle&#95;d5/near&#95;parallel/construct&#95;then&#95;angle</td>
<td>12.6432 ns [12.5324, 12.661] (0.95 confidence)</td>
<td>12.4101 ns [12.3984, 12.4392] (0.95 confidence)</td>
<td>1.01878</td>
</tr>
<tr>
<td>angle&#95;d5/near&#95;parallel/prepared&#95;vector</td>
<td>11.1876 ns [11.1856, 11.1947] (0.95 confidence)</td>
<td>11.1869 ns [11.1831, 11.1903] (0.95 confidence)</td>
<td>1.00006</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;adjacent/borrowed&#95;slice</td>
<td>17.3632 ns [17.3491, 17.3716] (0.95 confidence)</td>
<td>17.3532 ns [17.3476, 17.3714] (0.95 confidence)</td>
<td>1.00058</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;adjacent/construct&#95;then&#95;angle</td>
<td>16.3165 ns [16.2346, 16.4091] (0.95 confidence)</td>
<td>16.1572 ns [16.1563, 16.177] (0.95 confidence)</td>
<td>1.00986</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;adjacent/prepared&#95;vector</td>
<td>14.6401 ns [14.6358, 14.6408] (0.95 confidence)</td>
<td>14.6391 ns [14.6358, 14.6402] (0.95 confidence)</td>
<td>1.00007</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;alternating/borrowed&#95;slice</td>
<td>28.1857 ns [28.1599, 28.2409] (0.95 confidence)</td>
<td>28.1415 ns [28.1235, 28.1776] (0.95 confidence)</td>
<td>1.00157</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;alternating/construct&#95;then&#95;angle</td>
<td>27.1314 ns [27.1097, 27.1552] (0.95 confidence)</td>
<td>27.0116 ns [26.9871, 27.0826] (0.95 confidence)</td>
<td>1.00443</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;alternating/prepared&#95;vector</td>
<td>26.4511 ns [26.4054, 26.4772] (0.95 confidence)</td>
<td>26.4504 ns [26.4191, 26.5129] (0.95 confidence)</td>
<td>1.00003</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;antipodal/borrowed&#95;slice</td>
<td>13.3365 ns [13.3111, 13.3505] (0.95 confidence)</td>
<td>13.2833 ns [13.2445, 13.3048] (0.95 confidence)</td>
<td>1.004</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;antipodal/construct&#95;then&#95;angle</td>
<td>11.911 ns [11.9044, 11.9141] (0.95 confidence)</td>
<td>12.0966 ns [12.0538, 12.1109] (0.95 confidence)</td>
<td>0.984657</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;antipodal/prepared&#95;vector</td>
<td>10.3822 ns [10.374, 10.3942] (0.95 confidence)</td>
<td>10.3227 ns [10.3055, 10.3353] (0.95 confidence)</td>
<td>1.00576</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;identical/borrowed&#95;slice</td>
<td>11.6305 ns [11.6272, 11.6621] (0.95 confidence)</td>
<td>11.6496 ns [11.6286, 11.7165] (0.95 confidence)</td>
<td>0.998363</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;identical/construct&#95;then&#95;angle</td>
<td>10.6498 ns [10.5068, 10.6648] (0.95 confidence)</td>
<td>10.4491 ns [10.4145, 10.5393] (0.95 confidence)</td>
<td>1.0192</td>
</tr>
<tr>
<td>angle&#95;d5/sparse&#95;identical/prepared&#95;vector</td>
<td>8.69457 ns [8.683, 8.7145] (0.95 confidence)</td>
<td>8.6824 ns [8.68133, 8.70043] (0.95 confidence)</td>
<td>1.0014</td>
</tr>
<tr>
<td>angle&#95;d5/subnormal/borrowed&#95;slice</td>
<td>14.3051 ns [14.2258, 14.3241] (0.95 confidence)</td>
<td>14.3466 ns [14.3267, 14.3476] (0.95 confidence)</td>
<td>0.997106</td>
</tr>
<tr>
<td>angle&#95;d5/subnormal/construct&#95;then&#95;angle</td>
<td>12.6143 ns [12.4917, 12.6448] (0.95 confidence)</td>
<td>12.4002 ns [12.3983, 12.4093] (0.95 confidence)</td>
<td>1.01727</td>
</tr>
<tr>
<td>angle&#95;d5/subnormal/prepared&#95;vector</td>
<td>11.1864 ns [11.1843, 11.1942] (0.95 confidence)</td>
<td>11.1847 ns [11.1822, 11.1892] (0.95 confidence)</td>
<td>1.00015</td>
</tr>
<tr>
<td>angle&#95;d5/subnormal&#95;axes/borrowed&#95;slice</td>
<td>17.5637 ns [17.5555, 17.5645] (0.95 confidence)</td>
<td>17.5205 ns [17.4443, 17.5611] (0.95 confidence)</td>
<td>1.00246</td>
</tr>
<tr>
<td>angle&#95;d5/subnormal&#95;axes/construct&#95;then&#95;angle</td>
<td>15.8338 ns [15.7798, 15.8651] (0.95 confidence)</td>
<td>15.6342 ns [15.6192, 15.7686] (0.95 confidence)</td>
<td>1.01277</td>
</tr>
<tr>
<td>angle&#95;d5/subnormal&#95;axes/prepared&#95;vector</td>
<td>14.0816 ns [13.9492, 14.1281] (0.95 confidence)</td>
<td>14.1306 ns [14.1297, 14.1323] (0.95 confidence)</td>
<td>0.996532</td>
</tr>
<tr>
<td>angle&#95;d6/dense/borrowed&#95;slice</td>
<td>61.2063 ns [61.13, 61.3721] (0.95 confidence)</td>
<td>61.0878 ns [61.0133, 61.1909] (0.95 confidence)</td>
<td>1.00194</td>
</tr>
<tr>
<td>angle&#95;d6/dense/construct&#95;then&#95;angle</td>
<td>58.866 ns [58.8001, 58.9019] (0.95 confidence)</td>
<td>58.5914 ns [58.5025, 58.6679] (0.95 confidence)</td>
<td>1.00469</td>
</tr>
<tr>
<td>angle&#95;d6/dense/prepared&#95;vector</td>
<td>57.9158 ns [57.7826, 58.0011] (0.95 confidence)</td>
<td>57.9723 ns [57.7659, 58.0287] (0.95 confidence)</td>
<td>0.999024</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;antipodal/borrowed&#95;slice</td>
<td>47.9949 ns [47.9415, 48.0445] (0.95 confidence)</td>
<td>48.0253 ns [47.9451, 48.0554] (0.95 confidence)</td>
<td>0.999368</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;antipodal/construct&#95;then&#95;angle</td>
<td>45.5516 ns [45.5094, 45.586] (0.95 confidence)</td>
<td>45.3473 ns [45.3203, 45.3638] (0.95 confidence)</td>
<td>1.0045</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;antipodal/prepared&#95;vector</td>
<td>44.9783 ns [44.911, 45.0283] (0.95 confidence)</td>
<td>44.8838 ns [44.8072, 44.9404] (0.95 confidence)</td>
<td>1.00211</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;identical/borrowed&#95;slice</td>
<td>46.9019 ns [46.8663, 46.9348] (0.95 confidence)</td>
<td>46.7591 ns [46.7108, 46.8171] (0.95 confidence)</td>
<td>1.00305</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;identical/construct&#95;then&#95;angle</td>
<td>44.4029 ns [44.3379, 44.4234] (0.95 confidence)</td>
<td>44.2638 ns [44.2283, 44.3037] (0.95 confidence)</td>
<td>1.00314</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;identical/prepared&#95;vector</td>
<td>43.7301 ns [43.651, 43.7894] (0.95 confidence)</td>
<td>43.7554 ns [43.6969, 43.8591] (0.95 confidence)</td>
<td>0.999424</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;varied/borrowed&#95;slice</td>
<td>90.8163 ns [90.6295, 91.1596] (0.95 confidence)</td>
<td>90.7091 ns [90.4408, 90.8741] (0.95 confidence)</td>
<td>1.00118</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;varied/construct&#95;then&#95;angle</td>
<td>88.6482 ns [88.4253, 88.9732] (0.95 confidence)</td>
<td>89.0695 ns [88.7252, 89.2039] (0.95 confidence)</td>
<td>0.995271</td>
</tr>
<tr>
<td>angle&#95;d6/dense&#95;varied/prepared&#95;vector</td>
<td>88.0694 ns [87.8535, 88.2218] (0.95 confidence)</td>
<td>87.8889 ns [87.7148, 88.1232] (0.95 confidence)</td>
<td>1.00205</td>
</tr>
<tr>
<td>angle&#95;d6/maximum&#95;small&#95;angle/borrowed&#95;slice</td>
<td>15.5901 ns [15.5873, 15.5993] (0.95 confidence)</td>
<td>15.5859 ns [15.5586, 15.587] (0.95 confidence)</td>
<td>1.00027</td>
</tr>
<tr>
<td>angle&#95;d6/maximum&#95;small&#95;angle/construct&#95;then&#95;angle</td>
<td>13.859 ns [13.8563, 13.8803] (0.95 confidence)</td>
<td>13.8519 ns [13.7444, 13.8563] (0.95 confidence)</td>
<td>1.00051</td>
</tr>
<tr>
<td>angle&#95;d6/maximum&#95;small&#95;angle/prepared&#95;vector</td>
<td>11.9449 ns [11.9294, 11.9572] (0.95 confidence)</td>
<td>11.94 ns [11.9291, 11.9484] (0.95 confidence)</td>
<td>1.00041</td>
</tr>
<tr>
<td>angle&#95;d6/mixed&#95;scale/borrowed&#95;slice</td>
<td>17.0735 ns [17.0723, 17.0757] (0.95 confidence)</td>
<td>17.0715 ns [17.052, 17.0739] (0.95 confidence)</td>
<td>1.00011</td>
</tr>
<tr>
<td>angle&#95;d6/mixed&#95;scale/construct&#95;then&#95;angle</td>
<td>15.251 ns [15.177, 15.3144] (0.95 confidence)</td>
<td>15.1687 ns [15.1222, 15.3107] (0.95 confidence)</td>
<td>1.00543</td>
</tr>
<tr>
<td>angle&#95;d6/mixed&#95;scale/prepared&#95;vector</td>
<td>13.3704 ns [13.3681, 13.3787] (0.95 confidence)</td>
<td>13.3851 ns [13.3765, 13.3944] (0.95 confidence)</td>
<td>0.998907</td>
</tr>
<tr>
<td>angle&#95;d6/near&#95;antipodal/borrowed&#95;slice</td>
<td>21.3383 ns [21.3163, 21.3599] (0.95 confidence)</td>
<td>21.2256 ns [21.1905, 21.2529] (0.95 confidence)</td>
<td>1.00531</td>
</tr>
<tr>
<td>angle&#95;d6/near&#95;antipodal/construct&#95;then&#95;angle</td>
<td>18.8522 ns [18.8413, 18.8589] (0.95 confidence)</td>
<td>18.8416 ns [18.8357, 18.8471] (0.95 confidence)</td>
<td>1.00056</td>
</tr>
<tr>
<td>angle&#95;d6/near&#95;antipodal/prepared&#95;vector</td>
<td>17.7453 ns [17.7106, 17.7672] (0.95 confidence)</td>
<td>17.6729 ns [17.6383, 17.7136] (0.95 confidence)</td>
<td>1.0041</td>
</tr>
<tr>
<td>angle&#95;d6/near&#95;parallel/borrowed&#95;slice</td>
<td>15.5868 ns [15.5838, 15.5887] (0.95 confidence)</td>
<td>15.5875 ns [15.5854, 15.5894] (0.95 confidence)</td>
<td>0.99995</td>
</tr>
<tr>
<td>angle&#95;d6/near&#95;parallel/construct&#95;then&#95;angle</td>
<td>13.8556 ns [13.8535, 13.857] (0.95 confidence)</td>
<td>13.8358 ns [13.7618, 13.8584] (0.95 confidence)</td>
<td>1.00143</td>
</tr>
<tr>
<td>angle&#95;d6/near&#95;parallel/prepared&#95;vector</td>
<td>11.9398 ns [11.9309, 11.9519] (0.95 confidence)</td>
<td>11.9431 ns [11.9322, 11.9536] (0.95 confidence)</td>
<td>0.999724</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;adjacent/borrowed&#95;slice</td>
<td>19.1564 ns [19.1032, 19.2119] (0.95 confidence)</td>
<td>19.0459 ns [19.0209, 19.0865] (0.95 confidence)</td>
<td>1.00581</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;adjacent/construct&#95;then&#95;angle</td>
<td>17.6918 ns [17.6482, 17.7653] (0.95 confidence)</td>
<td>17.6694 ns [17.6564, 17.8005] (0.95 confidence)</td>
<td>1.00127</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;adjacent/prepared&#95;vector</td>
<td>15.6997 ns [15.6552, 15.8101] (0.95 confidence)</td>
<td>15.7559 ns [15.7124, 15.866] (0.95 confidence)</td>
<td>0.996432</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;alternating/borrowed&#95;slice</td>
<td>29.889 ns [29.8399, 29.9292] (0.95 confidence)</td>
<td>29.658 ns [29.6102, 29.6931] (0.95 confidence)</td>
<td>1.00779</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;alternating/construct&#95;then&#95;angle</td>
<td>28.1514 ns [28.1155, 28.178] (0.95 confidence)</td>
<td>28.0855 ns [28.0478, 28.1129] (0.95 confidence)</td>
<td>1.00235</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;alternating/prepared&#95;vector</td>
<td>27.2257 ns [27.1783, 27.2816] (0.95 confidence)</td>
<td>27.1 ns [27.0697, 27.1236] (0.95 confidence)</td>
<td>1.00464</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;antipodal/borrowed&#95;slice</td>
<td>14.2992 ns [14.2753, 14.3251] (0.95 confidence)</td>
<td>14.3014 ns [14.2561, 14.3215] (0.95 confidence)</td>
<td>0.999851</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;antipodal/construct&#95;then&#95;angle</td>
<td>13.0316 ns [13.012, 13.0389] (0.95 confidence)</td>
<td>13.0207 ns [13.0048, 13.036] (0.95 confidence)</td>
<td>1.00084</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;antipodal/prepared&#95;vector</td>
<td>11.0563 ns [11.0462, 11.0743] (0.95 confidence)</td>
<td>11.0066 ns [10.9869, 11.0273] (0.95 confidence)</td>
<td>1.00451</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;identical/borrowed&#95;slice</td>
<td>12.876 ns [12.8643, 12.9691] (0.95 confidence)</td>
<td>12.8709 ns [12.8644, 12.9239] (0.95 confidence)</td>
<td>1.00039</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;identical/construct&#95;then&#95;angle</td>
<td>11.7123 ns [11.6306, 11.8754] (0.95 confidence)</td>
<td>11.6287 ns [11.6181, 11.642] (0.95 confidence)</td>
<td>1.00719</td>
</tr>
<tr>
<td>angle&#95;d6/sparse&#95;identical/prepared&#95;vector</td>
<td>9.65962 ns [9.65692, 9.66306] (0.95 confidence)</td>
<td>9.65999 ns [9.65758, 9.66438] (0.95 confidence)</td>
<td>0.999962</td>
</tr>
<tr>
<td>angle&#95;d6/subnormal/borrowed&#95;slice</td>
<td>15.5868 ns [15.5796, 15.5915] (0.95 confidence)</td>
<td>15.547 ns [15.5217, 15.5869] (0.95 confidence)</td>
<td>1.00256</td>
</tr>
<tr>
<td>angle&#95;d6/subnormal/construct&#95;then&#95;angle</td>
<td>13.8258 ns [13.7109, 13.8536] (0.95 confidence)</td>
<td>13.8146 ns [13.6983, 13.8513] (0.95 confidence)</td>
<td>1.00081</td>
</tr>
<tr>
<td>angle&#95;d6/subnormal/prepared&#95;vector</td>
<td>11.9479 ns [11.9365, 11.9576] (0.95 confidence)</td>
<td>11.94 ns [11.9297, 11.9489] (0.95 confidence)</td>
<td>1.00067</td>
</tr>
<tr>
<td>angle&#95;d6/subnormal&#95;axes/borrowed&#95;slice</td>
<td>18.9538 ns [18.8169, 19.1527] (0.95 confidence)</td>
<td>18.859 ns [18.8129, 18.9551] (0.95 confidence)</td>
<td>1.00503</td>
</tr>
<tr>
<td>angle&#95;d6/subnormal&#95;axes/construct&#95;then&#95;angle</td>
<td>17.0769 ns [17.0496, 17.0782] (0.95 confidence)</td>
<td>16.9777 ns [16.8973, 17.0716] (0.95 confidence)</td>
<td>1.00585</td>
</tr>
<tr>
<td>angle&#95;d6/subnormal&#95;axes/prepared&#95;vector</td>
<td>15.1088 ns [15.1048, 15.1096] (0.95 confidence)</td>
<td>15.1114 ns [15.1095, 15.116] (0.95 confidence)</td>
<td>0.999829</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;det</td>
<td>395.788 ns [394.772, 396.203] (0.95 confidence)</td>
<td>399.853 ns [396.948, 403.292] (0.95 confidence)</td>
<td>0.989833</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>2.87018 ns [2.86629, 2.87658] (0.95 confidence)</td>
<td>2.93702 ns [2.93594, 2.94027] (0.95 confidence)</td>
<td>0.977242</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>3.27243 ns [3.26878, 3.2743] (0.95 confidence)</td>
<td>3.43345 ns [3.42742, 3.4468] (0.95 confidence)</td>
<td>0.953103</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>406.318 ns [405.124, 408.132] (0.95 confidence)</td>
<td>415.374 ns [408.509, 421.23] (0.95 confidence)</td>
<td>0.978199</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;dot</td>
<td>2.32688 ns [2.32637, 2.3296] (0.95 confidence)</td>
<td>2.63688 ns [2.43107, 2.77079] (0.95 confidence)</td>
<td>0.882435</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;inf&#95;norm</td>
<td>32.3338 ns [32.2839, 32.4265] (0.95 confidence)</td>
<td>45.0675 ns [44.2991, 46.6219] (0.95 confidence)</td>
<td>0.717451</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;ldlt</td>
<td>342.022 ns [341.409, 342.585] (0.95 confidence)</td>
<td>359.478 ns [356.612, 361.71] (0.95 confidence)</td>
<td>0.951441</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;ldlt&#95;solve</td>
<td>382.649 ns [382.337, 383.084] (0.95 confidence)</td>
<td>499.852 ns [491.709, 510.715] (0.95 confidence)</td>
<td>0.765525</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;lu</td>
<td>457.195 ns [404.227, 474.007] (0.95 confidence)</td>
<td>505.975 ns [498.092, 523.058] (0.95 confidence)</td>
<td>0.903592</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;norm2</td>
<td>9.9445 ns [9.87134, 10.1083] (0.95 confidence)</td>
<td>13.3845 ns [13.1587, 13.9337] (0.95 confidence)</td>
<td>0.742984</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>9.26137 ns [9.24269, 9.27948] (0.95 confidence)</td>
<td>12.229 ns [12.0359, 12.3577] (0.95 confidence)</td>
<td>0.757328</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>9.21526 ns [9.18561, 9.2344] (0.95 confidence)</td>
<td>12.2049 ns [11.9168, 12.3464] (0.95 confidence)</td>
<td>0.755046</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>4.89545 ns [4.84812, 4.91617] (0.95 confidence)</td>
<td>6.18841 ns [6.05891, 6.34033] (0.95 confidence)</td>
<td>0.791067</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>8.18397 ns [8.15239, 8.20281] (0.95 confidence)</td>
<td>10.9882 ns [10.8312, 11.3374] (0.95 confidence)</td>
<td>0.744793</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;norm2&#95;sq</td>
<td>2.09015 ns [2.08623, 2.09332] (0.95 confidence)</td>
<td>2.83759 ns [2.73855, 2.8932] (0.95 confidence)</td>
<td>0.736594</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>27.7326 ns [27.7237, 27.7489] (0.95 confidence)</td>
<td>28.7841 ns [28.6853, 28.8314] (0.95 confidence)</td>
<td>0.963469</td>
</tr>
<tr>
<td>d16/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>188.604 ns [188.27, 189.024] (0.95 confidence)</td>
<td>242.472 ns [238.263, 252.653] (0.95 confidence)</td>
<td>0.777837</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;det</td>
<td>0.572737 ns [0.57239, 0.573051] (0.95 confidence)</td>
<td>0.592601 ns [0.590712, 0.595267] (0.95 confidence)</td>
<td>0.966479</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>0.498842 ns [0.495496, 0.501311] (0.95 confidence)</td>
<td>0.514233 ns [0.514063, 0.51478] (0.95 confidence)</td>
<td>0.970071</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>0.519517 ns [0.519495, 0.519561] (0.95 confidence)</td>
<td>0.538578 ns [0.538493, 0.538683] (0.95 confidence)</td>
<td>0.964608</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>1.58267 ns [1.58128, 1.5841] (0.95 confidence)</td>
<td>1.5858 ns [1.5847, 1.58655] (0.95 confidence)</td>
<td>0.998031</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;dot</td>
<td>0.603376 ns [0.591942, 0.609279] (0.95 confidence)</td>
<td>0.634941 ns [0.63443, 0.635456] (0.95 confidence)</td>
<td>0.950286</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;inf&#95;norm</td>
<td>0.710366 ns [0.710219, 0.710806] (0.95 confidence)</td>
<td>0.712506 ns [0.712009, 0.713423] (0.95 confidence)</td>
<td>0.996996</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;ldlt</td>
<td>6.52787 ns [6.51763, 6.54709] (0.95 confidence)</td>
<td>5.21965 ns [5.19837, 5.23946] (0.95 confidence)</td>
<td>1.25063</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;ldlt&#95;solve</td>
<td>9.3538 ns [9.34649, 9.35723] (0.95 confidence)</td>
<td>10.8933 ns [10.8379, 10.9713] (0.95 confidence)</td>
<td>0.858674</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;lu</td>
<td>1.39355 ns [1.39297, 1.39957] (0.95 confidence)</td>
<td>1.40129 ns [1.39444, 1.41941] (0.95 confidence)</td>
<td>0.994478</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;norm2</td>
<td>1.17175 ns [1.15045, 1.18127] (0.95 confidence)</td>
<td>1.19412 ns [1.19377, 1.19497] (0.95 confidence)</td>
<td>0.98127</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>1.30167 ns [1.29845, 1.30651] (0.95 confidence)</td>
<td>1.41009 ns [1.40946, 1.41062] (0.95 confidence)</td>
<td>0.923109</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>1.30153 ns [1.29893, 1.30503] (0.95 confidence)</td>
<td>1.40935 ns [1.409, 1.40992] (0.95 confidence)</td>
<td>0.923502</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>1.12599 ns [1.12437, 1.12703] (0.95 confidence)</td>
<td>1.17911 ns [1.17172, 1.22044] (0.95 confidence)</td>
<td>0.954948</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>1.30314 ns [1.30142, 1.30419] (0.95 confidence)</td>
<td>1.38006 ns [1.37848, 1.38302] (0.95 confidence)</td>
<td>0.944263</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;norm2&#95;sq</td>
<td>0.394162 ns [0.393495, 0.394407] (0.95 confidence)</td>
<td>0.408793 ns [0.408382, 0.408966] (0.95 confidence)</td>
<td>0.964209</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>1.23487 ns [1.23348, 1.23668] (0.95 confidence)</td>
<td>1.26927 ns [1.26337, 1.27623] (0.95 confidence)</td>
<td>0.972898</td>
</tr>
<tr>
<td>d2/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>1.24955 ns [1.24585, 1.25485] (0.95 confidence)</td>
<td>1.32329 ns [1.31491, 1.33093] (0.95 confidence)</td>
<td>0.944278</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;det</td>
<td>1.33461 ns [1.33441, 1.33614] (0.95 confidence)</td>
<td>1.36002 ns [1.35403, 1.37225] (0.95 confidence)</td>
<td>0.98132</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>0.593729 ns [0.593713, 0.593745] (0.95 confidence)</td>
<td>0.601717 ns [0.601481, 0.602209] (0.95 confidence)</td>
<td>0.986724</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>0.69264 ns [0.692621, 0.692663] (0.95 confidence)</td>
<td>0.703566 ns [0.702937, 0.704513] (0.95 confidence)</td>
<td>0.98447</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>8.77476 ns [8.74856, 8.80286] (0.95 confidence)</td>
<td>8.8405 ns [8.83053, 8.85785] (0.95 confidence)</td>
<td>0.992564</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;dot</td>
<td>0.665643 ns [0.664461, 0.667266] (0.95 confidence)</td>
<td>0.677299 ns [0.675947, 0.678477] (0.95 confidence)</td>
<td>0.98279</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;inf&#95;norm</td>
<td>1.36186 ns [1.36179, 1.36206] (0.95 confidence)</td>
<td>1.37865 ns [1.37852, 1.37887] (0.95 confidence)</td>
<td>0.987826</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;ldlt</td>
<td>8.4107 ns [8.39399, 8.42756] (0.95 confidence)</td>
<td>8.6318 ns [8.61819, 8.64338] (0.95 confidence)</td>
<td>0.974384</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;ldlt&#95;solve</td>
<td>9.24038 ns [9.23171, 9.24434] (0.95 confidence)</td>
<td>9.42895 ns [9.39402, 9.47049] (0.95 confidence)</td>
<td>0.98</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;lu</td>
<td>8.18299 ns [8.16846, 8.19682] (0.95 confidence)</td>
<td>8.52234 ns [8.48226, 8.53144] (0.95 confidence)</td>
<td>0.960181</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;norm2</td>
<td>1.72792 ns [1.71189, 1.74871] (0.95 confidence)</td>
<td>1.77422 ns [1.77377, 1.77489] (0.95 confidence)</td>
<td>0.973903</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>1.81101 ns [1.80775, 1.81365] (0.95 confidence)</td>
<td>1.8855 ns [1.88417, 1.88679] (0.95 confidence)</td>
<td>0.96049</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>1.80336 ns [1.79941, 1.80773] (0.95 confidence)</td>
<td>1.88635 ns [1.885, 1.88779] (0.95 confidence)</td>
<td>0.956006</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>1.3632 ns [1.3624, 1.36478] (0.95 confidence)</td>
<td>1.39523 ns [1.39511, 1.3954] (0.95 confidence)</td>
<td>0.977038</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>1.82209 ns [1.81833, 1.82465] (0.95 confidence)</td>
<td>1.91214 ns [1.90884, 1.91557] (0.95 confidence)</td>
<td>0.952905</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;norm2&#95;sq</td>
<td>0.443475 ns [0.443383, 0.443652] (0.95 confidence)</td>
<td>0.446996 ns [0.446281, 0.448515] (0.95 confidence)</td>
<td>0.992122</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>1.77893 ns [1.77771, 1.77987] (0.95 confidence)</td>
<td>1.81273 ns [1.80827, 1.81692] (0.95 confidence)</td>
<td>0.981354</td>
</tr>
<tr>
<td>d3/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>2.04876 ns [2.04847, 2.05074] (0.95 confidence)</td>
<td>2.11026 ns [2.0888, 2.12607] (0.95 confidence)</td>
<td>0.970855</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;det</td>
<td>1983.87 ns [1979.76, 1987.68] (0.95 confidence)</td>
<td>3070.39 ns [2987.8, 3128.66] (0.95 confidence)</td>
<td>0.646128</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>7.70947 ns [7.69553, 7.72802] (0.95 confidence)</td>
<td>12.1552 ns [11.2226, 12.6761] (0.95 confidence)</td>
<td>0.634255</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>8.99916 ns [8.9862, 9.01331] (0.95 confidence)</td>
<td>9.30227 ns [9.2808, 9.45948] (0.95 confidence)</td>
<td>0.967415</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>2022.48 ns [2021.24, 2024.93] (0.95 confidence)</td>
<td>3000.48 ns [2935.33, 3056.92] (0.95 confidence)</td>
<td>0.674052</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;dot</td>
<td>3.98813 ns [3.98262, 3.98994] (0.95 confidence)</td>
<td>4.08122 ns [4.07905, 4.08259] (0.95 confidence)</td>
<td>0.977192</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;inf&#95;norm</td>
<td>127.149 ns [126.778, 127.484] (0.95 confidence)</td>
<td>130.909 ns [130.745, 131.129] (0.95 confidence)</td>
<td>0.971272</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;ldlt</td>
<td>2387.43 ns [2386.81, 2388.46] (0.95 confidence)</td>
<td>3490.46 ns [3466.29, 3516.61] (0.95 confidence)</td>
<td>0.683987</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;ldlt&#95;solve</td>
<td>2725.99 ns [2718.98, 2733.63] (0.95 confidence)</td>
<td>4054.77 ns [4008.6, 4106.23] (0.95 confidence)</td>
<td>0.672291</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;lu</td>
<td>2274.14 ns [2239.1, 2303.85] (0.95 confidence)</td>
<td>3433.49 ns [3369.52, 3498.16] (0.95 confidence)</td>
<td>0.662341</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;norm2</td>
<td>22.7844 ns [22.7114, 22.8603] (0.95 confidence)</td>
<td>24.8611 ns [24.8521, 24.8679] (0.95 confidence)</td>
<td>0.916468</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>18.5788 ns [18.5393, 18.6085] (0.95 confidence)</td>
<td>19.5895 ns [19.5749, 19.5925] (0.95 confidence)</td>
<td>0.948406</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>18.3995 ns [18.3406, 18.4293] (0.95 confidence)</td>
<td>19.5688 ns [19.5355, 19.5914] (0.95 confidence)</td>
<td>0.940247</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>9.11271 ns [9.08577, 9.13313] (0.95 confidence)</td>
<td>9.51544 ns [9.49635, 9.57847] (0.95 confidence)</td>
<td>0.957676</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>16.0115 ns [15.9755, 16.0276] (0.95 confidence)</td>
<td>16.8424 ns [16.8185, 16.8705] (0.95 confidence)</td>
<td>0.950664</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;norm2&#95;sq</td>
<td>3.92801 ns [3.92525, 3.93369] (0.95 confidence)</td>
<td>4.03888 ns [4.0333, 4.04115] (0.95 confidence)</td>
<td>0.972549</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>309.753 ns [308.903, 310.222] (0.95 confidence)</td>
<td>332.077 ns [325.204, 354.638] (0.95 confidence)</td>
<td>0.932775</td>
</tr>
<tr>
<td>d32/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>651.535 ns [649.804, 653.661] (0.95 confidence)</td>
<td>1079.43 ns [1064.79, 1091.55] (0.95 confidence)</td>
<td>0.603591</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;det</td>
<td>2.55171 ns [2.55015, 2.55497] (0.95 confidence)</td>
<td>2.58214 ns [2.57955, 2.58869] (0.95 confidence)</td>
<td>0.988218</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>0.76691 ns [0.766854, 0.766978] (0.95 confidence)</td>
<td>0.796264 ns [0.795006, 0.797129] (0.95 confidence)</td>
<td>0.963134</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>0.865898 ns [0.865841, 0.866004] (0.95 confidence)</td>
<td>0.880956 ns [0.880614, 0.881476] (0.95 confidence)</td>
<td>0.982908</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>14.5498 ns [14.5211, 14.574] (0.95 confidence)</td>
<td>14.7562 ns [14.7097, 14.7769] (0.95 confidence)</td>
<td>0.986012</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;dot</td>
<td>0.678324 ns [0.678071, 0.679073] (0.95 confidence)</td>
<td>0.742086 ns [0.710388, 0.755315] (0.95 confidence)</td>
<td>0.914078</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;inf&#95;norm</td>
<td>2.34352 ns [2.3408, 2.347] (0.95 confidence)</td>
<td>2.37977 ns [2.37962, 2.38008] (0.95 confidence)</td>
<td>0.984766</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;ldlt</td>
<td>17.8769 ns [17.8372, 17.9302] (0.95 confidence)</td>
<td>19.028 ns [18.9004, 19.0813] (0.95 confidence)</td>
<td>0.939509</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;ldlt&#95;solve</td>
<td>16.5324 ns [16.4971, 16.551] (0.95 confidence)</td>
<td>16.162 ns [16.097, 16.2168] (0.95 confidence)</td>
<td>1.02291</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;lu</td>
<td>14.1988 ns [14.1803, 14.2114] (0.95 confidence)</td>
<td>14.3807 ns [14.3618, 14.4073] (0.95 confidence)</td>
<td>0.987355</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;norm2</td>
<td>2.44103 ns [2.37792, 2.45179] (0.95 confidence)</td>
<td>2.53989 ns [2.45177, 2.54006] (0.95 confidence)</td>
<td>0.96108</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>2.37589 ns [2.37011, 2.38025] (0.95 confidence)</td>
<td>2.46805 ns [2.46788, 2.46822] (0.95 confidence)</td>
<td>0.962658</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>2.38541 ns [2.37841, 2.38683] (0.95 confidence)</td>
<td>2.47101 ns [2.4692, 2.47589] (0.95 confidence)</td>
<td>0.965356</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>1.47276 ns [1.46999, 1.4741] (0.95 confidence)</td>
<td>1.50287 ns [1.50248, 1.50358] (0.95 confidence)</td>
<td>0.979966</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>2.11932 ns [2.11358, 2.12154] (0.95 confidence)</td>
<td>2.21241 ns [2.20815, 2.2203] (0.95 confidence)</td>
<td>0.957922</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;norm2&#95;sq</td>
<td>0.457011 ns [0.456875, 0.457244] (0.95 confidence)</td>
<td>0.456252 ns [0.456079, 0.456553] (0.95 confidence)</td>
<td>1.00166</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>2.4615 ns [2.46009, 2.46734] (0.95 confidence)</td>
<td>2.67994 ns [2.66152, 2.70158] (0.95 confidence)</td>
<td>0.918488</td>
</tr>
<tr>
<td>d4/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>4.00706 ns [4.00438, 4.00973] (0.95 confidence)</td>
<td>4.10903 ns [4.09832, 4.12098] (0.95 confidence)</td>
<td>0.975186</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;det</td>
<td>28.5008 ns [28.4657, 28.5291] (0.95 confidence)</td>
<td>29.6596 ns [29.5923, 29.6783] (0.95 confidence)</td>
<td>0.960929</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>0.875565 ns [0.875313, 0.875815] (0.95 confidence)</td>
<td>0.905966 ns [0.903468, 0.908587] (0.95 confidence)</td>
<td>0.966444</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>1.03103 ns [1.03096, 1.03125] (0.95 confidence)</td>
<td>1.04649 ns [1.04524, 1.04792] (0.95 confidence)</td>
<td>0.98523</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>34.0321 ns [33.9418, 34.072] (0.95 confidence)</td>
<td>34.1885 ns [34.1748, 34.195] (0.95 confidence)</td>
<td>0.995424</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;dot</td>
<td>0.84194 ns [0.840582, 0.842292] (0.95 confidence)</td>
<td>0.821313 ns [0.816159, 0.844488] (0.95 confidence)</td>
<td>1.02511</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;inf&#95;norm</td>
<td>3.45412 ns [3.44998, 3.46314] (0.95 confidence)</td>
<td>3.51254 ns [3.51105, 3.51466] (0.95 confidence)</td>
<td>0.983369</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;ldlt</td>
<td>24.3917 ns [24.336, 24.4671] (0.95 confidence)</td>
<td>24.2796 ns [24.1231, 24.4375] (0.95 confidence)</td>
<td>1.00462</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;ldlt&#95;solve</td>
<td>30.0866 ns [30.0182, 30.0994] (0.95 confidence)</td>
<td>28.8552 ns [28.796, 28.9185] (0.95 confidence)</td>
<td>1.04267</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;lu</td>
<td>33.581 ns [33.5705, 33.5914] (0.95 confidence)</td>
<td>33.4432 ns [33.411, 33.4718] (0.95 confidence)</td>
<td>1.00412</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;norm2</td>
<td>2.93685 ns [2.87234, 2.94272] (0.95 confidence)</td>
<td>2.975 ns [2.87766, 3.04803] (0.95 confidence)</td>
<td>0.987177</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>3.13262 ns [3.12837, 3.13498] (0.95 confidence)</td>
<td>3.17453 ns [3.17324, 3.17738] (0.95 confidence)</td>
<td>0.986797</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>3.12946 ns [3.11904, 3.13326] (0.95 confidence)</td>
<td>3.17682 ns [3.17386, 3.18247] (0.95 confidence)</td>
<td>0.985093</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>1.71493 ns [1.71299, 1.71805] (0.95 confidence)</td>
<td>1.75764 ns [1.75541, 1.76361] (0.95 confidence)</td>
<td>0.975703</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>3.10057 ns [3.09551, 3.1064] (0.95 confidence)</td>
<td>2.66186 ns [2.65525, 2.66959] (0.95 confidence)</td>
<td>1.16481</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;norm2&#95;sq</td>
<td>0.524443 ns [0.52436, 0.524655] (0.95 confidence)</td>
<td>0.531349 ns [0.530745, 0.53281] (0.95 confidence)</td>
<td>0.987004</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>3.9621 ns [3.94854, 3.9716] (0.95 confidence)</td>
<td>3.94199 ns [3.93066, 3.95312] (0.95 confidence)</td>
<td>1.0051</td>
</tr>
<tr>
<td>d5/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>6.08064 ns [6.07096, 6.08713] (0.95 confidence)</td>
<td>6.25693 ns [6.23242, 6.29425] (0.95 confidence)</td>
<td>0.971824</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;det</td>
<td>14347.6 ns [14339.7, 14368.5] (0.95 confidence)</td>
<td>14529.4 ns [14508.9, 14545.6] (0.95 confidence)</td>
<td>0.98749</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>22.9062 ns [22.7412, 22.9998] (0.95 confidence)</td>
<td>23.7528 ns [23.672, 23.7632] (0.95 confidence)</td>
<td>0.964358</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>22.7588 ns [22.6881, 22.813] (0.95 confidence)</td>
<td>23.804 ns [23.7466, 23.8689] (0.95 confidence)</td>
<td>0.956095</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>14326.2 ns [14306.3, 14339.8] (0.95 confidence)</td>
<td>14396 ns [14362.9, 14439.9] (0.95 confidence)</td>
<td>0.995149</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;dot</td>
<td>10.8635 ns [10.8287, 10.8864] (0.95 confidence)</td>
<td>11.2526 ns [11.2295, 11.2688] (0.95 confidence)</td>
<td>0.965416</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;inf&#95;norm</td>
<td>620.228 ns [618.85, 621.288] (0.95 confidence)</td>
<td>848.677 ns [795.056, 871.222] (0.95 confidence)</td>
<td>0.730818</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;ldlt</td>
<td>16885.4 ns [16866, 16895.3] (0.95 confidence)</td>
<td>17074.5 ns [17050.3, 17091.4] (0.95 confidence)</td>
<td>0.98892</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;ldlt&#95;solve</td>
<td>18933.6 ns [18904.4, 18993.7] (0.95 confidence)</td>
<td>19237.7 ns [19138.9, 19271.7] (0.95 confidence)</td>
<td>0.984195</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;lu</td>
<td>14930.1 ns [14918.7, 14953.6] (0.95 confidence)</td>
<td>15253.5 ns [15242, 15270.6] (0.95 confidence)</td>
<td>0.9788</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;norm2</td>
<td>56.9506 ns [56.8312, 57.0868] (0.95 confidence)</td>
<td>61.8844 ns [61.2847, 62.3146] (0.95 confidence)</td>
<td>0.920274</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>36.037 ns [35.947, 36.1733] (0.95 confidence)</td>
<td>38.4379 ns [38.3903, 38.5195] (0.95 confidence)</td>
<td>0.937536</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>35.9241 ns [35.8419, 35.9792] (0.95 confidence)</td>
<td>40.0864 ns [40.0813, 40.0917] (0.95 confidence)</td>
<td>0.896165</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>16.9182 ns [16.8714, 16.9642] (0.95 confidence)</td>
<td>18.7786 ns [18.7008, 18.7798] (0.95 confidence)</td>
<td>0.900932</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>32.1889 ns [32.0647, 32.2869] (0.95 confidence)</td>
<td>39.9761 ns [38.1551, 40.923] (0.95 confidence)</td>
<td>0.805203</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;norm2&#95;sq</td>
<td>10.6835 ns [10.6558, 10.7034] (0.95 confidence)</td>
<td>11.2377 ns [11.203, 11.2811] (0.95 confidence)</td>
<td>0.950679</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>1013.34 ns [1011.5, 1013.95] (0.95 confidence)</td>
<td>1018.47 ns [1017.39, 1020.31] (0.95 confidence)</td>
<td>0.994963</td>
</tr>
<tr>
<td>d64/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>2492.65 ns [2486.66, 2496.49] (0.95 confidence)</td>
<td>2708.09 ns [2701.82, 2711.52] (0.95 confidence)</td>
<td>0.920448</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;det</td>
<td>81.7378 ns [81.5784, 81.9427] (0.95 confidence)</td>
<td>82.8783 ns [82.4098, 83.327] (0.95 confidence)</td>
<td>0.986239</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;det&#95;from&#95;ldlt</td>
<td>1.3035 ns [1.30206, 1.30493] (0.95 confidence)</td>
<td>1.33418 ns [1.33389, 1.33545] (0.95 confidence)</td>
<td>0.977005</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;det&#95;from&#95;ldlt&#95;balanced&#95;range</td>
<td>8.57129 ns [8.56349, 8.57667] (0.95 confidence)</td>
<td>8.83354 ns [8.81694, 8.85208] (0.95 confidence)</td>
<td>0.970312</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;det&#95;from&#95;lu</td>
<td>1.46214 ns [1.46185, 1.46289] (0.95 confidence)</td>
<td>1.56803 ns [1.55725, 1.57656] (0.95 confidence)</td>
<td>0.932468</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;det&#95;from&#95;lu&#95;balanced&#95;range</td>
<td>10.925 ns [10.8865, 11.0294] (0.95 confidence)</td>
<td>11.0271 ns [11.0202, 11.0498] (0.95 confidence)</td>
<td>0.99074</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;det&#95;via&#95;lu</td>
<td>89.4929 ns [88.9916, 90.0221] (0.95 confidence)</td>
<td>91.2184 ns [90.4118, 91.8218] (0.95 confidence)</td>
<td>0.981084</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;dot</td>
<td>0.932747 ns [0.931426, 0.934931] (0.95 confidence)</td>
<td>0.95199 ns [0.946786, 0.956131] (0.95 confidence)</td>
<td>0.979786</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;inf&#95;norm</td>
<td>8.36855 ns [8.35605, 8.37265] (0.95 confidence)</td>
<td>8.57885 ns [8.57844, 8.57933] (0.95 confidence)</td>
<td>0.975487</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;ldlt</td>
<td>90.3272 ns [89.4885, 91.3502] (0.95 confidence)</td>
<td>91.7333 ns [91.4548, 92.28] (0.95 confidence)</td>
<td>0.984672</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;ldlt&#95;ill&#95;conditioned</td>
<td>90.6874 ns [89.6459, 91.1185] (0.95 confidence)</td>
<td>91.3964 ns [90.9643, 92.714] (0.95 confidence)</td>
<td>0.992242</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;ldlt&#95;solve</td>
<td>102.095 ns [101.846, 103.344] (0.95 confidence)</td>
<td>104.94 ns [104.179, 105.529] (0.95 confidence)</td>
<td>0.97289</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;lu</td>
<td>89.9048 ns [89.6683, 90.0517] (0.95 confidence)</td>
<td>90.4285 ns [89.45, 91.1565] (0.95 confidence)</td>
<td>0.99421</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;lu&#95;ill&#95;conditioned</td>
<td>89.7633 ns [89.3494, 89.9307] (0.95 confidence)</td>
<td>90.671 ns [89.5761, 90.8653] (0.95 confidence)</td>
<td>0.989989</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;lu&#95;pivoting</td>
<td>94.1266 ns [93.8136, 94.2213] (0.95 confidence)</td>
<td>95.0717 ns [94.9176, 95.2174] (0.95 confidence)</td>
<td>0.990058</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;norm2</td>
<td>4.47619 ns [4.46182, 4.48067] (0.95 confidence)</td>
<td>4.75701 ns [4.74981, 4.75749] (0.95 confidence)</td>
<td>0.940966</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;norm2&#95;scenario&#95;descending</td>
<td>4.73652 ns [4.72531, 4.74823] (0.95 confidence)</td>
<td>5.00888 ns [5.00032, 5.01524] (0.95 confidence)</td>
<td>0.945625</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;norm2&#95;scenario&#95;repeated&#95;scale</td>
<td>4.74129 ns [4.73559, 4.75144] (0.95 confidence)</td>
<td>5.0124 ns [5.00419, 5.01951] (0.95 confidence)</td>
<td>0.945911</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;norm2&#95;scenario&#95;sparse</td>
<td>3.09619 ns [3.08446, 3.10505] (0.95 confidence)</td>
<td>3.25455 ns [3.25063, 3.25508] (0.95 confidence)</td>
<td>0.951343</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;norm2&#95;scenario&#95;wide&#95;dynamic&#95;range</td>
<td>4.17979 ns [4.17116, 4.18459] (0.95 confidence)</td>
<td>4.37159 ns [4.36782, 4.375] (0.95 confidence)</td>
<td>0.956125</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;norm2&#95;sq</td>
<td>0.664573 ns [0.664013, 0.666974] (0.95 confidence)</td>
<td>0.671153 ns [0.668506, 0.671688] (0.95 confidence)</td>
<td>0.990196</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;solve&#95;from&#95;ldlt</td>
<td>8.1 ns [8.09787, 8.10189] (0.95 confidence)</td>
<td>8.26301 ns [8.24185, 8.29232] (0.95 confidence)</td>
<td>0.980273</td>
</tr>
<tr>
<td>d8/la&#95;stack&#95;solve&#95;from&#95;lu</td>
<td>13.3244 ns [13.3094, 13.3379] (0.95 confidence)</td>
<td>13.982 ns [13.9614, 14.0029] (0.95 confidence)</td>
<td>0.952971</td>
</tr>
<tr>
<td>gram/1x2/dependent/la&#95;stack</td>
<td>0.681487 ns [0.680645, 0.682537] (0.95 confidence)</td>
<td>0.680557 ns [0.679396, 0.681628] (0.95 confidence)</td>
<td>1.00137</td>
</tr>
<tr>
<td>gram/1x2/mixed&#95;scale/la&#95;stack</td>
<td>0.682354 ns [0.680684, 0.683103] (0.95 confidence)</td>
<td>0.679452 ns [0.678535, 0.680405] (0.95 confidence)</td>
<td>1.00427</td>
</tr>
<tr>
<td>gram/1x2/near&#95;dependent/la&#95;stack</td>
<td>0.682061 ns [0.68137, 0.682864] (0.95 confidence)</td>
<td>0.679842 ns [0.678202, 0.680794] (0.95 confidence)</td>
<td>1.00326</td>
</tr>
<tr>
<td>gram/1x2/orthogonal/la&#95;stack</td>
<td>0.682565 ns [0.681498, 0.683294] (0.95 confidence)</td>
<td>0.680389 ns [0.679467, 0.681176] (0.95 confidence)</td>
<td>1.0032</td>
</tr>
<tr>
<td>gram/2x2/dependent/la&#95;stack</td>
<td>1.31097 ns [1.29889, 1.31514] (0.95 confidence)</td>
<td>1.3026 ns [1.29461, 1.31078] (0.95 confidence)</td>
<td>1.00643</td>
</tr>
<tr>
<td>gram/2x2/mixed&#95;scale/la&#95;stack</td>
<td>1.30996 ns [1.29679, 1.31389] (0.95 confidence)</td>
<td>1.31005 ns [1.29193, 1.31979] (0.95 confidence)</td>
<td>0.999928</td>
</tr>
<tr>
<td>gram/2x2/near&#95;dependent/la&#95;stack</td>
<td>1.30998 ns [1.29386, 1.31446] (0.95 confidence)</td>
<td>1.30999 ns [1.29715, 1.31698] (0.95 confidence)</td>
<td>0.99999</td>
</tr>
<tr>
<td>gram/2x2/orthogonal/la&#95;stack</td>
<td>1.30921 ns [1.29754, 1.31456] (0.95 confidence)</td>
<td>1.30806 ns [1.2961, 1.31167] (0.95 confidence)</td>
<td>1.00088</td>
</tr>
<tr>
<td>gram/2x3/dependent/la&#95;stack</td>
<td>1.26903 ns [1.26742, 1.26989] (0.95 confidence)</td>
<td>1.26251 ns [1.26157, 1.26445] (0.95 confidence)</td>
<td>1.00517</td>
</tr>
<tr>
<td>gram/2x3/mixed&#95;scale/la&#95;stack</td>
<td>1.26111 ns [1.26026, 1.2639] (0.95 confidence)</td>
<td>1.26165 ns [1.26028, 1.2636] (0.95 confidence)</td>
<td>0.999575</td>
</tr>
<tr>
<td>gram/2x3/near&#95;dependent/la&#95;stack</td>
<td>1.26572 ns [1.26388, 1.26801] (0.95 confidence)</td>
<td>1.26145 ns [1.26004, 1.26438] (0.95 confidence)</td>
<td>1.00338</td>
</tr>
<tr>
<td>gram/2x3/orthogonal/la&#95;stack</td>
<td>1.26354 ns [1.26185, 1.26852] (0.95 confidence)</td>
<td>1.26211 ns [1.26068, 1.26379] (0.95 confidence)</td>
<td>1.00113</td>
</tr>
<tr>
<td>gram/3x3/dependent/la&#95;stack</td>
<td>2.98674 ns [2.98228, 2.99563] (0.95 confidence)</td>
<td>2.97276 ns [2.96979, 2.97658] (0.95 confidence)</td>
<td>1.0047</td>
</tr>
<tr>
<td>gram/3x3/mixed&#95;scale/la&#95;stack</td>
<td>2.9873 ns [2.98295, 2.99127] (0.95 confidence)</td>
<td>2.97679 ns [2.96907, 2.98044] (0.95 confidence)</td>
<td>1.00353</td>
</tr>
<tr>
<td>gram/3x3/near&#95;dependent/la&#95;stack</td>
<td>2.97814 ns [2.97351, 2.98117] (0.95 confidence)</td>
<td>2.97504 ns [2.9697, 2.98031] (0.95 confidence)</td>
<td>1.00104</td>
</tr>
<tr>
<td>gram/3x3/orthogonal/la&#95;stack</td>
<td>2.98622 ns [2.98078, 2.99084] (0.95 confidence)</td>
<td>2.97498 ns [2.96881, 2.97896] (0.95 confidence)</td>
<td>1.00378</td>
</tr>
<tr>
<td>gram/3x4/dependent/la&#95;stack</td>
<td>3.04101 ns [3.03674, 3.04316] (0.95 confidence)</td>
<td>3.03699 ns [3.03526, 3.03923] (0.95 confidence)</td>
<td>1.00132</td>
</tr>
<tr>
<td>gram/3x4/mixed&#95;scale/la&#95;stack</td>
<td>3.03944 ns [3.03824, 3.04603] (0.95 confidence)</td>
<td>3.0601 ns [3.04896, 3.07818] (0.95 confidence)</td>
<td>0.993246</td>
</tr>
<tr>
<td>gram/3x4/near&#95;dependent/la&#95;stack</td>
<td>3.03749 ns [3.03445, 3.0415] (0.95 confidence)</td>
<td>3.02941 ns [3.02524, 3.03567] (0.95 confidence)</td>
<td>1.00267</td>
</tr>
<tr>
<td>gram/3x4/orthogonal/la&#95;stack</td>
<td>3.04158 ns [3.03949, 3.04425] (0.95 confidence)</td>
<td>3.02549 ns [3.02161, 3.02925] (0.95 confidence)</td>
<td>1.00532</td>
</tr>
<tr>
<td>gram/4x4/dependent/la&#95;stack</td>
<td>7.39853 ns [7.34977, 7.69908] (0.95 confidence)</td>
<td>7.37018 ns [7.3193, 7.63578] (0.95 confidence)</td>
<td>1.00385</td>
</tr>
<tr>
<td>gram/4x4/mixed&#95;scale/la&#95;stack</td>
<td>7.39049 ns [7.34744, 7.68011] (0.95 confidence)</td>
<td>7.38289 ns [7.33099, 7.65174] (0.95 confidence)</td>
<td>1.00103</td>
</tr>
<tr>
<td>gram/4x4/near&#95;dependent/la&#95;stack</td>
<td>7.38722 ns [7.32935, 7.66531] (0.95 confidence)</td>
<td>7.39199 ns [7.34196, 7.66985] (0.95 confidence)</td>
<td>0.999354</td>
</tr>
<tr>
<td>gram/4x4/orthogonal/la&#95;stack</td>
<td>7.4042 ns [7.34415, 7.68637] (0.95 confidence)</td>
<td>7.3527 ns [7.31069, 7.62563] (0.95 confidence)</td>
<td>1.00701</td>
</tr>
<tr>
<td>gram/4x5/dependent/la&#95;stack</td>
<td>9.44469 ns [9.4326, 9.45167] (0.95 confidence)</td>
<td>9.39146 ns [9.38047, 9.41248] (0.95 confidence)</td>
<td>1.00567</td>
</tr>
<tr>
<td>gram/4x5/mixed&#95;scale/la&#95;stack</td>
<td>9.43324 ns [9.41763, 9.45163] (0.95 confidence)</td>
<td>9.37149 ns [9.35313, 9.38336] (0.95 confidence)</td>
<td>1.00659</td>
</tr>
<tr>
<td>gram/4x5/near&#95;dependent/la&#95;stack</td>
<td>9.4104 ns [9.38877, 9.42777] (0.95 confidence)</td>
<td>9.37575 ns [9.36026, 9.39196] (0.95 confidence)</td>
<td>1.0037</td>
</tr>
<tr>
<td>gram/4x5/orthogonal/la&#95;stack</td>
<td>9.43848 ns [9.42344, 9.44675] (0.95 confidence)</td>
<td>9.36887 ns [9.33821, 9.39229] (0.95 confidence)</td>
<td>1.00743</td>
</tr>
<tr>
<td>gram/5x5/dependent/la&#95;stack</td>
<td>18.33 ns [18.2801, 18.3921] (0.95 confidence)</td>
<td>18.3027 ns [18.2637, 18.335] (0.95 confidence)</td>
<td>1.00149</td>
</tr>
<tr>
<td>gram/5x5/mixed&#95;scale/la&#95;stack</td>
<td>18.4031 ns [18.3412, 18.4469] (0.95 confidence)</td>
<td>18.4993 ns [18.4478, 18.5227] (0.95 confidence)</td>
<td>0.994801</td>
</tr>
<tr>
<td>gram/5x5/near&#95;dependent/la&#95;stack</td>
<td>18.3054 ns [18.2759, 18.3404] (0.95 confidence)</td>
<td>18.398 ns [18.3658, 18.4368] (0.95 confidence)</td>
<td>0.994971</td>
</tr>
<tr>
<td>gram/5x5/orthogonal/la&#95;stack</td>
<td>18.4424 ns [18.3938, 18.4973] (0.95 confidence)</td>
<td>18.3761 ns [18.3367, 18.434] (0.95 confidence)</td>
<td>1.00361</td>
</tr>
<tr>
<td>gram/5x6/dependent/la&#95;stack</td>
<td>19.161 ns [19.0917, 19.19] (0.95 confidence)</td>
<td>19.1067 ns [19.0825, 19.1319] (0.95 confidence)</td>
<td>1.00284</td>
</tr>
<tr>
<td>gram/5x6/mixed&#95;scale/la&#95;stack</td>
<td>19.1701 ns [19.0983, 19.2127] (0.95 confidence)</td>
<td>19.0354 ns [19.0003, 19.0711] (0.95 confidence)</td>
<td>1.00708</td>
</tr>
<tr>
<td>gram/5x6/near&#95;dependent/la&#95;stack</td>
<td>19.1567 ns [19.137, 19.1969] (0.95 confidence)</td>
<td>19.0572 ns [19.0196, 19.0988] (0.95 confidence)</td>
<td>1.00522</td>
</tr>
<tr>
<td>gram/5x6/orthogonal/la&#95;stack</td>
<td>19.0886 ns [19.0323, 19.1302] (0.95 confidence)</td>
<td>19.0651 ns [19.01, 19.1108] (0.95 confidence)</td>
<td>1.00123</td>
</tr>
<tr>
<td>gram/6x6/dependent/la&#95;stack</td>
<td>28.8129 ns [28.7215, 29.3277] (0.95 confidence)</td>
<td>28.7601 ns [28.6836, 28.9128] (0.95 confidence)</td>
<td>1.00184</td>
</tr>
<tr>
<td>gram/6x6/mixed&#95;scale/la&#95;stack</td>
<td>28.9258 ns [28.7092, 29.7781] (0.95 confidence)</td>
<td>28.6502 ns [28.5635, 28.7601] (0.95 confidence)</td>
<td>1.00962</td>
</tr>
<tr>
<td>gram/6x6/near&#95;dependent/la&#95;stack</td>
<td>29.4478 ns [28.7568, 29.7121] (0.95 confidence)</td>
<td>28.7083 ns [28.6662, 29.3316] (0.95 confidence)</td>
<td>1.02576</td>
</tr>
<tr>
<td>gram/6x6/orthogonal/la&#95;stack</td>
<td>28.8684 ns [28.785, 30.4444] (0.95 confidence)</td>
<td>28.9892 ns [28.8827, 29.9298] (0.95 confidence)</td>
<td>0.995832</td>
</tr>
<tr>
<td>gram/6x7/dependent/la&#95;stack</td>
<td>31.8939 ns [31.7256, 32.4372] (0.95 confidence)</td>
<td>31.7063 ns [31.659, 31.9684] (0.95 confidence)</td>
<td>1.00592</td>
</tr>
<tr>
<td>gram/6x7/mixed&#95;scale/la&#95;stack</td>
<td>31.9793 ns [31.8758, 32.316] (0.95 confidence)</td>
<td>31.6694 ns [31.5444, 31.7173] (0.95 confidence)</td>
<td>1.00979</td>
</tr>
<tr>
<td>gram/6x7/near&#95;dependent/la&#95;stack</td>
<td>31.944 ns [31.7126, 32.3743] (0.95 confidence)</td>
<td>31.7526 ns [31.7068, 31.9434] (0.95 confidence)</td>
<td>1.00603</td>
</tr>
<tr>
<td>gram/6x7/orthogonal/la&#95;stack</td>
<td>31.5769 ns [31.4478, 32.5141] (0.95 confidence)</td>
<td>31.8139 ns [31.4186, 32.6061] (0.95 confidence)</td>
<td>0.992552</td>
</tr>
<tr>
<td>gram/7x7/dependent/la&#95;stack</td>
<td>40.0685 ns [39.8111, 40.1393] (0.95 confidence)</td>
<td>40.0461 ns [39.8727, 40.1802] (0.95 confidence)</td>
<td>1.00056</td>
</tr>
<tr>
<td>gram/7x7/mixed&#95;scale/la&#95;stack</td>
<td>40.365 ns [40.1288, 40.4459] (0.95 confidence)</td>
<td>40.228 ns [40.0785, 40.31] (0.95 confidence)</td>
<td>1.00341</td>
</tr>
<tr>
<td>gram/7x7/near&#95;dependent/la&#95;stack</td>
<td>39.814 ns [39.6787, 40.1292] (0.95 confidence)</td>
<td>39.9258 ns [39.8615, 40.0887] (0.95 confidence)</td>
<td>0.9972</td>
</tr>
<tr>
<td>gram/7x7/orthogonal/la&#95;stack</td>
<td>40.0013 ns [39.8227, 40.1501] (0.95 confidence)</td>
<td>40.6232 ns [39.9245, 40.7828] (0.95 confidence)</td>
<td>0.984691</td>
</tr>
<tr>
<td>gram/7x8/dependent/la&#95;stack</td>
<td>41.2976 ns [41.2652, 41.3809] (0.95 confidence)</td>
<td>42.4064 ns [42.1543, 42.8265] (0.95 confidence)</td>
<td>0.973853</td>
</tr>
<tr>
<td>gram/7x8/mixed&#95;scale/la&#95;stack</td>
<td>41.4011 ns [41.2665, 41.4374] (0.95 confidence)</td>
<td>41.7018 ns [41.3915, 41.8946] (0.95 confidence)</td>
<td>0.99279</td>
</tr>
<tr>
<td>gram/7x8/near&#95;dependent/la&#95;stack</td>
<td>41.1279 ns [41.0161, 41.3322] (0.95 confidence)</td>
<td>41.5185 ns [41.4251, 41.8595] (0.95 confidence)</td>
<td>0.990591</td>
</tr>
<tr>
<td>gram/7x8/orthogonal/la&#95;stack</td>
<td>41.1327 ns [41.0518, 41.2656] (0.95 confidence)</td>
<td>41.2821 ns [41.1988, 41.3475] (0.95 confidence)</td>
<td>0.996381</td>
</tr>
<tr>
<td>gram/8x8/dependent/la&#95;stack</td>
<td>51.6781 ns [51.6383, 51.7303] (0.95 confidence)</td>
<td>51.7105 ns [51.6273, 51.7996] (0.95 confidence)</td>
<td>0.999374</td>
</tr>
<tr>
<td>gram/8x8/mixed&#95;scale/la&#95;stack</td>
<td>51.7696 ns [51.7272, 51.8555] (0.95 confidence)</td>
<td>51.5458 ns [51.4416, 51.6644] (0.95 confidence)</td>
<td>1.00434</td>
</tr>
<tr>
<td>gram/8x8/near&#95;dependent/la&#95;stack</td>
<td>51.7764 ns [51.7071, 51.8455] (0.95 confidence)</td>
<td>51.6106 ns [51.5597, 51.6676] (0.95 confidence)</td>
<td>1.00321</td>
</tr>
<tr>
<td>gram/8x8/orthogonal/la&#95;stack</td>
<td>51.7818 ns [51.7017, 51.8508] (0.95 confidence)</td>
<td>51.6202 ns [51.606, 51.6408] (0.95 confidence)</td>
<td>1.00313</td>
</tr>
<tr>
<td>interval&#95;d3/dense&#95;inconclusive</td>
<td>41.1097 ns [41.0615, 41.2089] (0.95 confidence)</td>
<td>41.1077 ns [41.0197, 41.2516] (0.95 confidence)</td>
<td>1.00005</td>
</tr>
<tr>
<td>interval&#95;d3/dense&#95;sign</td>
<td>35.6857 ns [35.316, 35.7734] (0.95 confidence)</td>
<td>35.035 ns [34.8003, 35.3013] (0.95 confidence)</td>
<td>1.01857</td>
</tr>
<tr>
<td>interval&#95;d3/lifted&#95;assembly</td>
<td>8.36119 ns [8.35385, 8.36642] (0.95 confidence)</td>
<td>8.35787 ns [8.34706, 8.36352] (0.95 confidence)</td>
<td>1.0004</td>
</tr>
<tr>
<td>interval&#95;d3/lifted&#95;assembly&#95;and&#95;sign</td>
<td>37.8972 ns [37.7335, 38.0135] (0.95 confidence)</td>
<td>37.7926 ns [37.6363, 37.9136] (0.95 confidence)</td>
<td>1.00277</td>
</tr>
<tr>
<td>interval&#95;d3/sparse&#95;sign</td>
<td>31.3017 ns [31.0767, 31.4583] (0.95 confidence)</td>
<td>31.3833 ns [31.2895, 31.4699] (0.95 confidence)</td>
<td>0.9974</td>
</tr>
<tr>
<td>interval&#95;d4/dense&#95;inconclusive</td>
<td>99.836 ns [99.0972, 100.379] (0.95 confidence)</td>
<td>101.643 ns [101.091, 102.597] (0.95 confidence)</td>
<td>0.98222</td>
</tr>
<tr>
<td>interval&#95;d4/dense&#95;sign</td>
<td>96.4784 ns [95.9026, 97.1418] (0.95 confidence)</td>
<td>98.9517 ns [98.5302, 99.8874] (0.95 confidence)</td>
<td>0.975005</td>
</tr>
<tr>
<td>interval&#95;d4/lifted&#95;assembly</td>
<td>13.2158 ns [13.1946, 13.2406] (0.95 confidence)</td>
<td>13.1725 ns [13.1629, 13.1854] (0.95 confidence)</td>
<td>1.00329</td>
</tr>
<tr>
<td>interval&#95;d4/lifted&#95;assembly&#95;and&#95;sign</td>
<td>81.4894 ns [81.2546, 81.627] (0.95 confidence)</td>
<td>81.0787 ns [80.6327, 81.3619] (0.95 confidence)</td>
<td>1.00507</td>
</tr>
<tr>
<td>interval&#95;d4/sparse&#95;sign</td>
<td>69.5583 ns [68.8466, 71.243] (0.95 confidence)</td>
<td>67.508 ns [67.0811, 67.8957] (0.95 confidence)</td>
<td>1.03037</td>
</tr>
<tr>
<td>interval&#95;d5/dense&#95;inconclusive</td>
<td>235.48 ns [234.23, 236.996] (0.95 confidence)</td>
<td>235.687 ns [234.327, 236.845] (0.95 confidence)</td>
<td>0.999122</td>
</tr>
<tr>
<td>interval&#95;d5/dense&#95;sign</td>
<td>223.555 ns [222.136, 224.295] (0.95 confidence)</td>
<td>227.991 ns [224.866, 229.891] (0.95 confidence)</td>
<td>0.980545</td>
</tr>
<tr>
<td>interval&#95;d5/lifted&#95;assembly</td>
<td>21.3936 ns [21.3163, 21.5424] (0.95 confidence)</td>
<td>21.337 ns [21.2769, 21.5756] (0.95 confidence)</td>
<td>1.00265</td>
</tr>
<tr>
<td>interval&#95;d5/lifted&#95;assembly&#95;and&#95;sign</td>
<td>192.291 ns [190.972, 192.98] (0.95 confidence)</td>
<td>188.14 ns [186.211, 190.19] (0.95 confidence)</td>
<td>1.02206</td>
</tr>
<tr>
<td>interval&#95;d5/sparse&#95;sign</td>
<td>153.309 ns [152.261, 155.015] (0.95 confidence)</td>
<td>151.27 ns [150.714, 151.801] (0.95 confidence)</td>
<td>1.01348</td>
</tr>
<tr>
<td>interval&#95;d6/dense&#95;inconclusive</td>
<td>539.22 ns [536.758, 543.951] (0.95 confidence)</td>
<td>524.283 ns [521.836, 526.768] (0.95 confidence)</td>
<td>1.02849</td>
</tr>
<tr>
<td>interval&#95;d6/dense&#95;sign</td>
<td>512.514 ns [511.062, 514.392] (0.95 confidence)</td>
<td>513.524 ns [511.893, 517.406] (0.95 confidence)</td>
<td>0.998033</td>
</tr>
<tr>
<td>interval&#95;d6/lifted&#95;assembly</td>
<td>30.5373 ns [30.5056, 30.5802] (0.95 confidence)</td>
<td>30.5396 ns [30.4867, 30.603] (0.95 confidence)</td>
<td>0.999925</td>
</tr>
<tr>
<td>interval&#95;d6/lifted&#95;assembly&#95;and&#95;sign</td>
<td>376.591 ns [375.839, 377.774] (0.95 confidence)</td>
<td>377.561 ns [376.121, 378.458] (0.95 confidence)</td>
<td>0.99743</td>
</tr>
<tr>
<td>interval&#95;d6/sparse&#95;sign</td>
<td>361.485 ns [359.374, 364.095] (0.95 confidence)</td>
<td>361.952 ns [359.227, 364.429] (0.95 confidence)</td>
<td>0.998708</td>
</tr>
<tr>
<td>interval&#95;d7/dense&#95;inconclusive</td>
<td>1094.04 ns [1082.11, 1112.29] (0.95 confidence)</td>
<td>1089.6 ns [1083.05, 1098.42] (0.95 confidence)</td>
<td>1.00408</td>
</tr>
<tr>
<td>interval&#95;d7/dense&#95;sign</td>
<td>1093.28 ns [1080.18, 1111.59] (0.95 confidence)</td>
<td>1049.05 ns [1043.49, 1056.36] (0.95 confidence)</td>
<td>1.04216</td>
</tr>
<tr>
<td>interval&#95;d7/lifted&#95;assembly</td>
<td>37.9483 ns [37.7098, 38.1172] (0.95 confidence)</td>
<td>37.1771 ns [37.108, 37.2334] (0.95 confidence)</td>
<td>1.02075</td>
</tr>
<tr>
<td>interval&#95;d7/lifted&#95;assembly&#95;and&#95;sign</td>
<td>784.466 ns [782.074, 786.878] (0.95 confidence)</td>
<td>767.51 ns [766.786, 769.88] (0.95 confidence)</td>
<td>1.02209</td>
</tr>
<tr>
<td>interval&#95;d7/sparse&#95;sign</td>
<td>742.266 ns [739.385, 744.961] (0.95 confidence)</td>
<td>740.445 ns [737.453, 742.755] (0.95 confidence)</td>
<td>1.00246</td>
</tr>
<tr>
<td>interval&#95;det&#95;sign/d4&#95;conclusive&#95;lifted</td>
<td>70.4132 ns [70.0404, 70.754] (0.95 confidence)</td>
<td>65.9754 ns [65.7773, 66.2901] (0.95 confidence)</td>
<td>1.06726</td>
</tr>
<tr>
<td>interval&#95;det&#95;sign/d4&#95;inconclusive&#95;lifted</td>
<td>67.1653 ns [66.6651, 67.4921] (0.95 confidence)</td>
<td>64.6679 ns [64.5836, 64.9085] (0.95 confidence)</td>
<td>1.03862</td>
</tr>
<tr>
<td>interval&#95;det&#95;sign/d7&#95;conclusive&#95;lifted</td>
<td>736.556 ns [733.493, 737.552] (0.95 confidence)</td>
<td>742.221 ns [740.208, 743.723] (0.95 confidence)</td>
<td>0.992367</td>
</tr>
<tr>
<td>interval&#95;scalar/cancellation/add</td>
<td>2.10355 ns [2.10092, 2.10691] (0.95 confidence)</td>
<td>2.10021 ns [2.09786, 2.10196] (0.95 confidence)</td>
<td>1.00159</td>
</tr>
<tr>
<td>interval&#95;scalar/cancellation/multiply</td>
<td>1.86626 ns [1.86153, 1.87623] (0.95 confidence)</td>
<td>1.85883 ns [1.84922, 1.86698] (0.95 confidence)</td>
<td>1.004</td>
</tr>
<tr>
<td>interval&#95;scalar/cancellation/square</td>
<td>1.85242 ns [1.84918, 1.85511] (0.95 confidence)</td>
<td>1.84166 ns [1.83917, 1.84986] (0.95 confidence)</td>
<td>1.00584</td>
</tr>
<tr>
<td>interval&#95;scalar/exact/add</td>
<td>2.10466 ns [2.10282, 2.10874] (0.95 confidence)</td>
<td>2.09887 ns [2.09542, 2.101] (0.95 confidence)</td>
<td>1.00276</td>
</tr>
<tr>
<td>interval&#95;scalar/exact/multiply</td>
<td>2.07506 ns [2.05276, 2.0891] (0.95 confidence)</td>
<td>2.09624 ns [2.07386, 2.11394] (0.95 confidence)</td>
<td>0.989894</td>
</tr>
<tr>
<td>interval&#95;scalar/exact/square</td>
<td>1.85492 ns [1.84961, 1.85856] (0.95 confidence)</td>
<td>1.841 ns [1.83813, 1.84362] (0.95 confidence)</td>
<td>1.00756</td>
</tr>
<tr>
<td>interval&#95;scalar/point/add</td>
<td>2.17129 ns [2.17022, 2.17235] (0.95 confidence)</td>
<td>2.16497 ns [2.16088, 2.17045] (0.95 confidence)</td>
<td>1.00292</td>
</tr>
<tr>
<td>interval&#95;scalar/point/multiply</td>
<td>2.56719 ns [2.55807, 2.57205] (0.95 confidence)</td>
<td>2.56992 ns [2.56013, 2.57713] (0.95 confidence)</td>
<td>0.998937</td>
</tr>
<tr>
<td>interval&#95;scalar/point/square</td>
<td>2.10548 ns [2.104, 2.10711] (0.95 confidence)</td>
<td>2.10043 ns [2.094, 2.10246] (0.95 confidence)</td>
<td>1.0024</td>
</tr>
<tr>
<td>interval&#95;scalar/subnormal/add</td>
<td>2.40126 ns [2.39938, 2.40221] (0.95 confidence)</td>
<td>2.39344 ns [2.39011, 2.3976] (0.95 confidence)</td>
<td>1.00327</td>
</tr>
<tr>
<td>interval&#95;scalar/subnormal/multiply</td>
<td>2.83864 ns [2.82765, 2.86563] (0.95 confidence)</td>
<td>2.84784 ns [2.83129, 2.8684] (0.95 confidence)</td>
<td>0.996768</td>
</tr>
<tr>
<td>interval&#95;scalar/subnormal/square</td>
<td>2.37053 ns [2.36557, 2.37219] (0.95 confidence)</td>
<td>2.36943 ns [2.36553, 2.37231] (0.95 confidence)</td>
<td>1.00046</td>
</tr>
<tr>
<td>interval&#95;scalar/zero/add</td>
<td>0.933287 ns [0.93207, 0.934387] (0.95 confidence)</td>
<td>0.931723 ns [0.930066, 0.932874] (0.95 confidence)</td>
<td>1.00168</td>
</tr>
<tr>
<td>interval&#95;scalar/zero/multiply</td>
<td>1.74332 ns [1.73981, 1.74685] (0.95 confidence)</td>
<td>1.74139 ns [1.73745, 1.74555] (0.95 confidence)</td>
<td>1.00111</td>
</tr>
<tr>
<td>interval&#95;scalar/zero/square</td>
<td>0.87802 ns [0.875854, 0.879296] (0.95 confidence)</td>
<td>0.875606 ns [0.87496, 0.876032] (0.95 confidence)</td>
<td>1.00276</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;batch/origin/construction</td>
<td>2.0552 ns [2.05193, 2.05913] (0.95 confidence)</td>
<td>2.04919 ns [2.04368, 2.05696] (0.95 confidence)</td>
<td>1.00293</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;batch/origin/prepared</td>
<td>13.8219 ns [13.8125, 13.8302] (0.95 confidence)</td>
<td>13.8208 ns [13.8167, 13.839] (0.95 confidence)</td>
<td>1.00008</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;batch/origin/with&#95;construction</td>
<td>15.3234 ns [15.3035, 15.3264] (0.95 confidence)</td>
<td>15.3236 ns [15.3197, 15.3292] (0.95 confidence)</td>
<td>0.999982</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;batch/translated/construction</td>
<td>2.05076 ns [2.04579, 2.05479] (0.95 confidence)</td>
<td>2.04986 ns [2.04603, 2.05246] (0.95 confidence)</td>
<td>1.00044</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;batch/translated/prepared</td>
<td>40.4141 ns [40.3331, 40.4853] (0.95 confidence)</td>
<td>40.2164 ns [40.1114, 40.2812] (0.95 confidence)</td>
<td>1.00492</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;batch/translated/with&#95;construction</td>
<td>42.2661 ns [42.1976, 42.3471] (0.95 confidence)</td>
<td>42.1614 ns [41.956, 42.3164] (0.95 confidence)</td>
<td>1.00248</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/cancellation/difference</td>
<td>10.3145 ns [10.3108, 10.3205] (0.95 confidence)</td>
<td>10.2979 ns [10.2898, 10.3358] (0.95 confidence)</td>
<td>1.00162</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/cancellation/difference&#95;endpoints</td>
<td>11.4836 ns [11.4398, 11.6226] (0.95 confidence)</td>
<td>11.4822 ns [11.429, 11.5232] (0.95 confidence)</td>
<td>1.00012</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/cancellation/dot</td>
<td>5.94459 ns [5.93094, 5.95126] (0.95 confidence)</td>
<td>5.89861 ns [5.87038, 5.92303] (0.95 confidence)</td>
<td>1.00779</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/cancellation/dot&#95;endpoints</td>
<td>6.99481 ns [6.98416, 7.00487] (0.95 confidence)</td>
<td>6.97745 ns [6.96751, 6.98809] (0.95 confidence)</td>
<td>1.00249</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense/difference</td>
<td>10.288 ns [10.279, 10.3073] (0.95 confidence)</td>
<td>10.2064 ns [10.1942, 10.2261] (0.95 confidence)</td>
<td>1.008</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense/difference&#95;endpoints</td>
<td>12.4604 ns [12.4316, 12.4749] (0.95 confidence)</td>
<td>12.398 ns [12.3699, 12.4204] (0.95 confidence)</td>
<td>1.00504</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense/dot</td>
<td>5.50318 ns [5.49701, 5.51317] (0.95 confidence)</td>
<td>5.44595 ns [5.43578, 5.4525] (0.95 confidence)</td>
<td>1.01051</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense/dot&#95;endpoints</td>
<td>6.43784 ns [6.42888, 6.44093] (0.95 confidence)</td>
<td>6.40395 ns [6.39394, 6.41421] (0.95 confidence)</td>
<td>1.00529</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense&#95;inexact/difference</td>
<td>10.3074 ns [10.2975, 10.3172] (0.95 confidence)</td>
<td>10.2479 ns [10.237, 10.2661] (0.95 confidence)</td>
<td>1.00581</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense&#95;inexact/difference&#95;endpoints</td>
<td>11.4948 ns [11.479, 11.5099] (0.95 confidence)</td>
<td>11.442 ns [11.4246, 11.45] (0.95 confidence)</td>
<td>1.00461</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense&#95;inexact/dot</td>
<td>5.5113 ns [5.49616, 5.51566] (0.95 confidence)</td>
<td>5.49857 ns [5.48763, 5.51048] (0.95 confidence)</td>
<td>1.00231</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/dense&#95;inexact/dot&#95;endpoints</td>
<td>6.45673 ns [6.4332, 6.47361] (0.95 confidence)</td>
<td>6.46966 ns [6.45504, 6.49104] (0.95 confidence)</td>
<td>0.998001</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/overflow/dot</td>
<td>1.0687 ns [1.06839, 1.06916] (0.95 confidence)</td>
<td>1.06942 ns [1.06886, 1.06977] (0.95 confidence)</td>
<td>0.999329</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/sparse/difference</td>
<td>4.73105 ns [4.72849, 4.73478] (0.95 confidence)</td>
<td>4.86982 ns [4.86625, 4.87513] (0.95 confidence)</td>
<td>0.971503</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/sparse/difference&#95;endpoints</td>
<td>5.55379 ns [5.54931, 5.56027] (0.95 confidence)</td>
<td>5.56721 ns [5.5533, 5.58803] (0.95 confidence)</td>
<td>0.997591</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/sparse/dot</td>
<td>2.22588 ns [2.2214, 2.22628] (0.95 confidence)</td>
<td>2.22607 ns [2.22267, 2.2263] (0.95 confidence)</td>
<td>0.999913</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/sparse/dot&#95;endpoints</td>
<td>2.22621 ns [2.22425, 2.22635] (0.95 confidence)</td>
<td>2.22617 ns [2.22217, 2.2263] (0.95 confidence)</td>
<td>1.00002</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/underflow/difference</td>
<td>4.45617 ns [4.45377, 4.45736] (0.95 confidence)</td>
<td>4.45849 ns [4.45601, 4.45952] (0.95 confidence)</td>
<td>0.999479</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/underflow/difference&#95;endpoints</td>
<td>4.60181 ns [4.59938, 4.60311] (0.95 confidence)</td>
<td>4.60216 ns [4.59878, 4.60531] (0.95 confidence)</td>
<td>0.999925</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/underflow/dot</td>
<td>1.47822 ns [1.4745, 1.48066] (0.95 confidence)</td>
<td>1.4757 ns [1.47234, 1.47874] (0.95 confidence)</td>
<td>1.0017</td>
</tr>
<tr>
<td>linear&#95;form&#95;d2&#95;prepared/underflow/dot&#95;endpoints</td>
<td>1.48407 ns [1.48277, 1.48421] (0.95 confidence)</td>
<td>1.48413 ns [1.4831, 1.48421] (0.95 confidence)</td>
<td>0.999959</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;batch/origin/construction</td>
<td>2.29563 ns [2.28365, 2.30027] (0.95 confidence)</td>
<td>2.25766 ns [2.25322, 2.26012] (0.95 confidence)</td>
<td>1.01681</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;batch/origin/prepared</td>
<td>24.6501 ns [24.4957, 24.797] (0.95 confidence)</td>
<td>23.2519 ns [23.2423, 23.2538] (0.95 confidence)</td>
<td>1.06013</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;batch/origin/with&#95;construction</td>
<td>28.4208 ns [28.007, 28.6767] (0.95 confidence)</td>
<td>27.4876 ns [27.4621, 27.5451] (0.95 confidence)</td>
<td>1.03395</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;batch/translated/construction</td>
<td>2.3029 ns [2.29406, 2.30948] (0.95 confidence)</td>
<td>2.22787 ns [2.21025, 2.2512] (0.95 confidence)</td>
<td>1.03368</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;batch/translated/prepared</td>
<td>70.035 ns [69.932, 70.2371] (0.95 confidence)</td>
<td>68.0131 ns [67.7889, 68.7687] (0.95 confidence)</td>
<td>1.02973</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;batch/translated/with&#95;construction</td>
<td>72.203 ns [71.9869, 72.4448] (0.95 confidence)</td>
<td>71.0447 ns [70.9223, 71.1131] (0.95 confidence)</td>
<td>1.0163</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/cancellation/difference</td>
<td>10.9256 ns [10.9224, 10.9296] (0.95 confidence)</td>
<td>10.9106 ns [10.9063, 10.914] (0.95 confidence)</td>
<td>1.00137</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/cancellation/difference&#95;endpoints</td>
<td>12.1146 ns [12.0663, 12.1523] (0.95 confidence)</td>
<td>11.9929 ns [11.9206, 12.0707] (0.95 confidence)</td>
<td>1.01015</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/cancellation/dot</td>
<td>6.02758 ns [6.01671, 6.03622] (0.95 confidence)</td>
<td>6.01352 ns [5.98705, 6.02586] (0.95 confidence)</td>
<td>1.00234</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/cancellation/dot&#95;endpoints</td>
<td>6.57293 ns [6.56418, 6.5777] (0.95 confidence)</td>
<td>6.6297 ns [6.59787, 6.67365] (0.95 confidence)</td>
<td>0.991436</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense/difference</td>
<td>14.2613 ns [14.2418, 14.3109] (0.95 confidence)</td>
<td>14.2061 ns [14.1822, 14.2237] (0.95 confidence)</td>
<td>1.00389</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense/difference&#95;endpoints</td>
<td>16.4962 ns [16.4741, 16.5211] (0.95 confidence)</td>
<td>16.5483 ns [16.4816, 16.753] (0.95 confidence)</td>
<td>0.996853</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense/dot</td>
<td>7.97493 ns [7.96338, 7.98758] (0.95 confidence)</td>
<td>7.84146 ns [7.82368, 7.90643] (0.95 confidence)</td>
<td>1.01702</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense/dot&#95;endpoints</td>
<td>8.37167 ns [8.35643, 8.37749] (0.95 confidence)</td>
<td>8.38871 ns [8.36888, 8.40131] (0.95 confidence)</td>
<td>0.997968</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense&#95;inexact/difference</td>
<td>14.2494 ns [14.2305, 14.2757] (0.95 confidence)</td>
<td>14.2483 ns [14.2234, 14.266] (0.95 confidence)</td>
<td>1.00008</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense&#95;inexact/difference&#95;endpoints</td>
<td>15.9192 ns [15.8987, 15.9549] (0.95 confidence)</td>
<td>15.8334 ns [15.8079, 15.8563] (0.95 confidence)</td>
<td>1.00542</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense&#95;inexact/dot</td>
<td>8.00504 ns [7.99672, 8.01127] (0.95 confidence)</td>
<td>7.9455 ns [7.91287, 7.99322] (0.95 confidence)</td>
<td>1.00749</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/dense&#95;inexact/dot&#95;endpoints</td>
<td>9.28952 ns [9.24835, 9.42108] (0.95 confidence)</td>
<td>9.30917 ns [9.25198, 9.40796] (0.95 confidence)</td>
<td>0.997889</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/overflow/dot</td>
<td>0.999574 ns [0.998639, 1.00005] (0.95 confidence)</td>
<td>0.99652 ns [0.995572, 0.997727] (0.95 confidence)</td>
<td>1.00306</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/sparse/difference</td>
<td>5.44911 ns [5.44694, 5.45097] (0.95 confidence)</td>
<td>5.47226 ns [5.46339, 5.48129] (0.95 confidence)</td>
<td>0.99577</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/sparse/difference&#95;endpoints</td>
<td>6.42163 ns [6.41359, 6.42382] (0.95 confidence)</td>
<td>6.4004 ns [6.35648, 6.41552] (0.95 confidence)</td>
<td>1.00332</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/sparse/dot</td>
<td>2.72196 ns [2.7212, 2.72216] (0.95 confidence)</td>
<td>2.72429 ns [2.72226, 2.74979] (0.95 confidence)</td>
<td>0.999147</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/sparse/dot&#95;endpoints</td>
<td>2.74206 ns [2.74192, 2.74239] (0.95 confidence)</td>
<td>2.96821 ns [2.96241, 2.96864] (0.95 confidence)</td>
<td>0.92381</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/underflow/difference</td>
<td>4.70308 ns [4.70205, 4.70424] (0.95 confidence)</td>
<td>4.70294 ns [4.70158, 4.70397] (0.95 confidence)</td>
<td>1.00003</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/underflow/difference&#95;endpoints</td>
<td>4.97165 ns [4.97023, 5.03997] (0.95 confidence)</td>
<td>4.97285 ns [4.97044, 4.97857] (0.95 confidence)</td>
<td>0.999758</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/underflow/dot</td>
<td>3.13088 ns [3.12573, 3.13431] (0.95 confidence)</td>
<td>3.13201 ns [3.12429, 3.13407] (0.95 confidence)</td>
<td>0.999638</td>
</tr>
<tr>
<td>linear&#95;form&#95;d3&#95;prepared/underflow/dot&#95;endpoints</td>
<td>3.16335 ns [3.15921, 3.16689] (0.95 confidence)</td>
<td>3.12481 ns [3.11852, 3.13104] (0.95 confidence)</td>
<td>1.01233</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4/dot&#95;bounded&#95;inconclusive</td>
<td>6.19639 ns [6.19001, 6.20605] (0.95 confidence)</td>
<td>6.28295 ns [6.27385, 6.29402] (0.95 confidence)</td>
<td>0.986222</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4/dot&#95;bounded&#95;well&#95;separated</td>
<td>9.53848 ns [9.50843, 9.63863] (0.95 confidence)</td>
<td>9.5681 ns [9.55671, 9.58591] (0.95 confidence)</td>
<td>0.996905</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4/dot&#95;difference&#95;bounded&#95;inconclusive</td>
<td>20.9611 ns [20.9271, 20.9887] (0.95 confidence)</td>
<td>21.0524 ns [21.0408, 21.0821] (0.95 confidence)</td>
<td>0.995663</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4/dot&#95;difference&#95;bounded&#95;well&#95;separated</td>
<td>16.6654 ns [16.6306, 16.7192] (0.95 confidence)</td>
<td>16.742 ns [16.7273, 16.7669] (0.95 confidence)</td>
<td>0.99543</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4/dot&#95;plain&#95;well&#95;separated</td>
<td>0.682772 ns [0.68174, 0.684286] (0.95 confidence)</td>
<td>0.684917 ns [0.683723, 0.687081] (0.95 confidence)</td>
<td>0.996868</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;batch/origin/construction</td>
<td>2.3771 ns [2.36615, 2.39811] (0.95 confidence)</td>
<td>2.39596 ns [2.3845, 2.40832] (0.95 confidence)</td>
<td>0.99213</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;batch/origin/prepared</td>
<td>35.4302 ns [35.1538, 35.5842] (0.95 confidence)</td>
<td>36.1771 ns [36.0146, 36.3497] (0.95 confidence)</td>
<td>0.979355</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;batch/origin/with&#95;construction</td>
<td>42.0214 ns [41.8178, 42.3286] (0.95 confidence)</td>
<td>42.6551 ns [42.2894, 42.8396] (0.95 confidence)</td>
<td>0.985143</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;batch/translated/construction</td>
<td>2.36555 ns [2.35974, 2.38211] (0.95 confidence)</td>
<td>2.41093 ns [2.40677, 2.4197] (0.95 confidence)</td>
<td>0.981181</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;batch/translated/prepared</td>
<td>111.669 ns [111.202, 112.239] (0.95 confidence)</td>
<td>111.845 ns [111.286, 112.18] (0.95 confidence)</td>
<td>0.998428</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;batch/translated/with&#95;construction</td>
<td>116.588 ns [116.481, 116.755] (0.95 confidence)</td>
<td>116.961 ns [116.49, 117.594] (0.95 confidence)</td>
<td>0.996809</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/cancellation/difference</td>
<td>11.7837 ns [11.7126, 11.8479] (0.95 confidence)</td>
<td>11.9197 ns [11.7578, 12.0686] (0.95 confidence)</td>
<td>0.988587</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/cancellation/difference&#95;endpoints</td>
<td>12.9863 ns [12.9084, 13.1764] (0.95 confidence)</td>
<td>12.906 ns [12.8824, 12.9672] (0.95 confidence)</td>
<td>1.00622</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/cancellation/dot</td>
<td>6.57086 ns [6.53801, 6.60812] (0.95 confidence)</td>
<td>6.35822 ns [6.34018, 6.36982] (0.95 confidence)</td>
<td>1.03344</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/cancellation/dot&#95;endpoints</td>
<td>7.09107 ns [7.06706, 7.12381] (0.95 confidence)</td>
<td>6.99345 ns [6.98621, 7.01995] (0.95 confidence)</td>
<td>1.01396</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense/difference</td>
<td>18.6993 ns [18.6512, 18.7583] (0.95 confidence)</td>
<td>18.2389 ns [18.2032, 18.2739] (0.95 confidence)</td>
<td>1.02525</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense/difference&#95;endpoints</td>
<td>21.2181 ns [21.0597, 21.354] (0.95 confidence)</td>
<td>20.8375 ns [20.7779, 20.9391] (0.95 confidence)</td>
<td>1.01827</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense/dot</td>
<td>9.76413 ns [9.74049, 9.80526] (0.95 confidence)</td>
<td>9.50022 ns [9.48459, 9.52221] (0.95 confidence)</td>
<td>1.02778</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense/dot&#95;endpoints</td>
<td>11.3194 ns [11.2763, 11.3614] (0.95 confidence)</td>
<td>11.0156 ns [10.9733, 11.0379] (0.95 confidence)</td>
<td>1.02758</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense&#95;inexact/difference</td>
<td>18.8402 ns [18.8161, 18.8988] (0.95 confidence)</td>
<td>18.886 ns [18.8083, 18.9757] (0.95 confidence)</td>
<td>0.997573</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense&#95;inexact/difference&#95;endpoints</td>
<td>20.5396 ns [20.4849, 20.653] (0.95 confidence)</td>
<td>20.5769 ns [20.4702, 20.7475] (0.95 confidence)</td>
<td>0.998188</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense&#95;inexact/dot</td>
<td>9.71435 ns [9.69418, 9.75667] (0.95 confidence)</td>
<td>9.58426 ns [9.56461, 9.614] (0.95 confidence)</td>
<td>1.01357</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/dense&#95;inexact/dot&#95;endpoints</td>
<td>10.4541 ns [10.4408, 10.4915] (0.95 confidence)</td>
<td>10.3127 ns [10.2913, 10.3522] (0.95 confidence)</td>
<td>1.01372</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/overflow/dot</td>
<td>1.00049 ns [0.998626, 1.00154] (0.95 confidence)</td>
<td>1.00781 ns [1.00597, 1.01006] (0.95 confidence)</td>
<td>0.992733</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/sparse/difference</td>
<td>6.42232 ns [6.39466, 6.45769] (0.95 confidence)</td>
<td>6.29216 ns [6.28517, 6.29828] (0.95 confidence)</td>
<td>1.02069</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/sparse/difference&#95;endpoints</td>
<td>7.48959 ns [7.44284, 7.52531] (0.95 confidence)</td>
<td>7.41295 ns [7.39292, 7.44321] (0.95 confidence)</td>
<td>1.01034</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/sparse/dot</td>
<td>3.21912 ns [3.17999, 3.23837] (0.95 confidence)</td>
<td>3.22608 ns [3.20777, 3.24572] (0.95 confidence)</td>
<td>0.997843</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/sparse/dot&#95;endpoints</td>
<td>3.2487 ns [3.24553, 3.25621] (0.95 confidence)</td>
<td>3.24137 ns [3.23428, 3.24629] (0.95 confidence)</td>
<td>1.00226</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/underflow/difference</td>
<td>5.74614 ns [5.73329, 5.76812] (0.95 confidence)</td>
<td>5.76465 ns [5.75711, 5.7707] (0.95 confidence)</td>
<td>0.996789</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/underflow/difference&#95;endpoints</td>
<td>6.24114 ns [6.23643, 6.24675] (0.95 confidence)</td>
<td>6.24871 ns [6.19653, 6.2731] (0.95 confidence)</td>
<td>0.998789</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/underflow/dot</td>
<td>3.96697 ns [3.9573, 3.97186] (0.95 confidence)</td>
<td>3.9533 ns [3.94619, 3.95728] (0.95 confidence)</td>
<td>1.00346</td>
</tr>
<tr>
<td>linear&#95;form&#95;d4&#95;prepared/underflow/dot&#95;endpoints</td>
<td>3.97866 ns [3.97115, 3.98269] (0.95 confidence)</td>
<td>3.96648 ns [3.95618, 3.97664] (0.95 confidence)</td>
<td>1.00307</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;batch/origin/construction</td>
<td>2.86613 ns [2.85772, 2.88505] (0.95 confidence)</td>
<td>2.9304 ns [2.91253, 2.94281] (0.95 confidence)</td>
<td>0.978067</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;batch/origin/prepared</td>
<td>47.2579 ns [47.022, 47.4763] (0.95 confidence)</td>
<td>48.4073 ns [47.7824, 49.3192] (0.95 confidence)</td>
<td>0.976256</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;batch/origin/with&#95;construction</td>
<td>59.278 ns [58.972, 59.8347] (0.95 confidence)</td>
<td>60.0483 ns [59.5653, 60.773] (0.95 confidence)</td>
<td>0.987172</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;batch/translated/construction</td>
<td>3.28613 ns [3.27393, 3.29277] (0.95 confidence)</td>
<td>3.31292 ns [3.30698, 3.31731] (0.95 confidence)</td>
<td>0.991913</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;batch/translated/prepared</td>
<td>147.466 ns [146.917, 150.364] (0.95 confidence)</td>
<td>147.884 ns [147.712, 148.054] (0.95 confidence)</td>
<td>0.997174</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;batch/translated/with&#95;construction</td>
<td>157.131 ns [156.822, 157.42] (0.95 confidence)</td>
<td>158.167 ns [157.887, 158.724] (0.95 confidence)</td>
<td>0.993451</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/cancellation/difference</td>
<td>12.6137 ns [12.4051, 12.9539] (0.95 confidence)</td>
<td>12.374 ns [12.2446, 12.5593] (0.95 confidence)</td>
<td>1.01937</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/cancellation/difference&#95;endpoints</td>
<td>13.5033 ns [13.3965, 13.6233] (0.95 confidence)</td>
<td>13.6217 ns [13.5374, 13.7483] (0.95 confidence)</td>
<td>0.99131</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/cancellation/dot</td>
<td>6.77983 ns [6.69537, 6.9628] (0.95 confidence)</td>
<td>6.80961 ns [6.77612, 6.85507] (0.95 confidence)</td>
<td>0.995627</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/cancellation/dot&#95;endpoints</td>
<td>7.31495 ns [7.29516, 7.34682] (0.95 confidence)</td>
<td>7.30902 ns [7.28496, 7.33775] (0.95 confidence)</td>
<td>1.00081</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense/difference</td>
<td>23.6518 ns [23.2894, 23.9791] (0.95 confidence)</td>
<td>23.075 ns [23.0322, 23.1167] (0.95 confidence)</td>
<td>1.02499</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense/difference&#95;endpoints</td>
<td>24.5835 ns [24.5217, 24.6678] (0.95 confidence)</td>
<td>24.5235 ns [24.471, 24.5644] (0.95 confidence)</td>
<td>1.00245</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense/dot</td>
<td>11.8055 ns [11.7659, 11.8481] (0.95 confidence)</td>
<td>11.808 ns [11.7835, 11.8579] (0.95 confidence)</td>
<td>0.999791</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense/dot&#95;endpoints</td>
<td>12.4632 ns [12.4414, 12.4912] (0.95 confidence)</td>
<td>12.4957 ns [12.4808, 12.5241] (0.95 confidence)</td>
<td>0.997397</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense&#95;inexact/difference</td>
<td>23.0815 ns [23.014, 23.1743] (0.95 confidence)</td>
<td>23.2526 ns [23.2104, 23.4028] (0.95 confidence)</td>
<td>0.992644</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense&#95;inexact/difference&#95;endpoints</td>
<td>25.8421 ns [25.7153, 26.2065] (0.95 confidence)</td>
<td>25.9684 ns [25.8577, 26.0769] (0.95 confidence)</td>
<td>0.995134</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense&#95;inexact/dot</td>
<td>11.8276 ns [11.7716, 11.8578] (0.95 confidence)</td>
<td>11.7771 ns [11.7396, 11.8069] (0.95 confidence)</td>
<td>1.00429</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/dense&#95;inexact/dot&#95;endpoints</td>
<td>13.4906 ns [13.4086, 13.6636] (0.95 confidence)</td>
<td>13.4641 ns [13.4354, 13.5089] (0.95 confidence)</td>
<td>1.00197</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/overflow/dot</td>
<td>1.00244 ns [1.00112, 1.0039] (0.95 confidence)</td>
<td>1.01274 ns [1.00925, 1.01469] (0.95 confidence)</td>
<td>0.989835</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/sparse/difference</td>
<td>7.27731 ns [7.26928, 7.29506] (0.95 confidence)</td>
<td>7.34811 ns [7.31796, 7.39655] (0.95 confidence)</td>
<td>0.990365</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/sparse/difference&#95;endpoints</td>
<td>8.44813 ns [8.42376, 8.47343] (0.95 confidence)</td>
<td>8.48781 ns [8.47135, 8.49864] (0.95 confidence)</td>
<td>0.995325</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/sparse/dot</td>
<td>3.74502 ns [3.69841, 3.77736] (0.95 confidence)</td>
<td>3.65576 ns [3.65392, 3.66434] (0.95 confidence)</td>
<td>1.02441</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/sparse/dot&#95;endpoints</td>
<td>3.77928 ns [3.76957, 3.78334] (0.95 confidence)</td>
<td>3.82065 ns [3.80912, 3.83745] (0.95 confidence)</td>
<td>0.989174</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/underflow/difference</td>
<td>7.76213 ns [7.75034, 7.78078] (0.95 confidence)</td>
<td>7.82391 ns [7.80858, 7.84414] (0.95 confidence)</td>
<td>0.992103</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/underflow/difference&#95;endpoints</td>
<td>8.35735 ns [8.30683, 8.40819] (0.95 confidence)</td>
<td>8.4297 ns [8.36898, 8.48793] (0.95 confidence)</td>
<td>0.991417</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/underflow/dot</td>
<td>4.66395 ns [4.65624, 4.66811] (0.95 confidence)</td>
<td>4.73154 ns [4.72605, 4.74202] (0.95 confidence)</td>
<td>0.985716</td>
</tr>
<tr>
<td>linear&#95;form&#95;d5&#95;prepared/underflow/dot&#95;endpoints</td>
<td>4.65272 ns [4.64565, 4.66214] (0.95 confidence)</td>
<td>4.7225 ns [4.71222, 4.72933] (0.95 confidence)</td>
<td>0.985223</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;batch/origin/construction</td>
<td>3.12364 ns [3.12047, 3.12962] (0.95 confidence)</td>
<td>3.15126 ns [3.14681, 3.15968] (0.95 confidence)</td>
<td>0.991236</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;batch/origin/prepared</td>
<td>62.2839 ns [61.845, 62.7274] (0.95 confidence)</td>
<td>63.6574 ns [63.1406, 64.147] (0.95 confidence)</td>
<td>0.978423</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;batch/origin/with&#95;construction</td>
<td>82.2614 ns [81.9278, 82.9114] (0.95 confidence)</td>
<td>84.4803 ns [83.4949, 84.7162] (0.95 confidence)</td>
<td>0.973734</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;batch/translated/construction</td>
<td>3.60834 ns [3.60513, 3.61181] (0.95 confidence)</td>
<td>3.62063 ns [3.61898, 3.62232] (0.95 confidence)</td>
<td>0.996606</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;batch/translated/prepared</td>
<td>207.647 ns [207.049, 208.03] (0.95 confidence)</td>
<td>209.77 ns [209.255, 210.514] (0.95 confidence)</td>
<td>0.989876</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;batch/translated/with&#95;construction</td>
<td>231.5 ns [230.952, 231.833] (0.95 confidence)</td>
<td>234.351 ns [233.596, 235.168] (0.95 confidence)</td>
<td>0.987834</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/cancellation/difference</td>
<td>12.7219 ns [12.6929, 12.7661] (0.95 confidence)</td>
<td>12.8108 ns [12.7622, 12.8706] (0.95 confidence)</td>
<td>0.993057</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/cancellation/difference&#95;endpoints</td>
<td>14.1417 ns [14.0696, 14.3208] (0.95 confidence)</td>
<td>14.0908 ns [14.0593, 14.1453] (0.95 confidence)</td>
<td>1.00361</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/cancellation/dot</td>
<td>6.9179 ns [6.90945, 6.93367] (0.95 confidence)</td>
<td>6.93243 ns [6.92507, 6.94312] (0.95 confidence)</td>
<td>0.997904</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/cancellation/dot&#95;endpoints</td>
<td>7.74173 ns [7.7244, 7.7912] (0.95 confidence)</td>
<td>7.74324 ns [7.71485, 7.76405] (0.95 confidence)</td>
<td>0.999805</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense/difference</td>
<td>29.3384 ns [29.2883, 29.3959] (0.95 confidence)</td>
<td>29.4837 ns [29.4196, 29.5359] (0.95 confidence)</td>
<td>0.995073</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense/difference&#95;endpoints</td>
<td>32.8187 ns [32.7173, 32.8778] (0.95 confidence)</td>
<td>32.8676 ns [32.8359, 32.9413] (0.95 confidence)</td>
<td>0.998512</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense/dot</td>
<td>13.8932 ns [13.8652, 13.9429] (0.95 confidence)</td>
<td>14.0848 ns [14.0641, 14.1252] (0.95 confidence)</td>
<td>0.986395</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense/dot&#95;endpoints</td>
<td>15.3696 ns [15.3384, 15.3935] (0.95 confidence)</td>
<td>15.5605 ns [15.5192, 15.5985] (0.95 confidence)</td>
<td>0.987734</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense&#95;inexact/difference</td>
<td>29.3818 ns [29.3346, 29.4311] (0.95 confidence)</td>
<td>29.6306 ns [29.5881, 29.7285] (0.95 confidence)</td>
<td>0.991605</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense&#95;inexact/difference&#95;endpoints</td>
<td>32.1313 ns [32.0521, 32.1888] (0.95 confidence)</td>
<td>32.2784 ns [32.2036, 32.3155] (0.95 confidence)</td>
<td>0.995443</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense&#95;inexact/dot</td>
<td>13.971 ns [13.9101, 13.9974] (0.95 confidence)</td>
<td>14.0126 ns [13.9847, 14.0406] (0.95 confidence)</td>
<td>0.997029</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/dense&#95;inexact/dot&#95;endpoints</td>
<td>15.4053 ns [15.3671, 15.4701] (0.95 confidence)</td>
<td>15.6122 ns [15.5542, 15.7493] (0.95 confidence)</td>
<td>0.986751</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/overflow/dot</td>
<td>1.00219 ns [0.999386, 1.0037] (0.95 confidence)</td>
<td>1.01005 ns [1.00791, 1.01158] (0.95 confidence)</td>
<td>0.992223</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/sparse/difference</td>
<td>8.15165 ns [8.13687, 8.16904] (0.95 confidence)</td>
<td>8.18368 ns [8.1632, 8.2142] (0.95 confidence)</td>
<td>0.996086</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/sparse/difference&#95;endpoints</td>
<td>9.42126 ns [9.36846, 9.50884] (0.95 confidence)</td>
<td>9.59333 ns [9.56905, 9.61888] (0.95 confidence)</td>
<td>0.982063</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/sparse/dot</td>
<td>4.12612 ns [4.12356, 4.13095] (0.95 confidence)</td>
<td>4.14305 ns [4.13856, 4.15275] (0.95 confidence)</td>
<td>0.995915</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/sparse/dot&#95;endpoints</td>
<td>4.19577 ns [4.1921, 4.19913] (0.95 confidence)</td>
<td>4.25016 ns [4.24199, 4.25895] (0.95 confidence)</td>
<td>0.987203</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/underflow/difference</td>
<td>8.70112 ns [8.69302, 8.71169] (0.95 confidence)</td>
<td>8.74621 ns [8.73371, 8.76377] (0.95 confidence)</td>
<td>0.994844</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/underflow/difference&#95;endpoints</td>
<td>9.18663 ns [9.17934, 9.20024] (0.95 confidence)</td>
<td>9.24031 ns [9.22392, 9.27777] (0.95 confidence)</td>
<td>0.994191</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/underflow/dot</td>
<td>5.37868 ns [5.36805, 5.38945] (0.95 confidence)</td>
<td>5.481 ns [5.47325, 5.48605] (0.95 confidence)</td>
<td>0.981331</td>
</tr>
<tr>
<td>linear&#95;form&#95;d6&#95;prepared/underflow/dot&#95;endpoints</td>
<td>5.45555 ns [5.43747, 5.4618] (0.95 confidence)</td>
<td>5.54068 ns [5.53635, 5.54572] (0.95 confidence)</td>
<td>0.984635</td>
</tr>
<tr>
<td>public&#95;d2/det&#95;bound</td>
<td>1.94387 ns [1.94087, 1.94892] (0.95 confidence)</td>
<td>1.94307 ns [1.9412, 1.94799] (0.95 confidence)</td>
<td>1.00041</td>
</tr>
<tr>
<td>public&#95;d2/det&#95;direct</td>
<td>1.14085 ns [1.13983, 1.14173] (0.95 confidence)</td>
<td>1.14286 ns [1.14237, 1.14335] (0.95 confidence)</td>
<td>0.998247</td>
</tr>
<tr>
<td>public&#95;d2/det&#95;with&#95;bound</td>
<td>2.16788 ns [2.16135, 2.17128] (0.95 confidence)</td>
<td>2.1386 ns [2.13047, 2.14537] (0.95 confidence)</td>
<td>1.01369</td>
</tr>
<tr>
<td>public&#95;d2/first&#95;asymmetry</td>
<td>1.38944 ns [1.38811, 1.39045] (0.95 confidence)</td>
<td>1.40027 ns [1.39885, 1.40159] (0.95 confidence)</td>
<td>0.992268</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;det</td>
<td>7.84608 ns [7.77509, 7.97999] (0.95 confidence)</td>
<td>7.96887 ns [7.84865, 7.99613] (0.95 confidence)</td>
<td>0.984591</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;access</td>
<td>6.44154 ns [6.42179, 6.45069] (0.95 confidence)</td>
<td>6.53482 ns [6.5265, 6.5443] (0.95 confidence)</td>
<td>0.985726</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;from&#95;matrix</td>
<td>2.61535 ns [2.6084, 2.63174] (0.95 confidence)</td>
<td>2.65526 ns [2.64711, 2.66434] (0.95 confidence)</td>
<td>0.984967</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;from&#95;points</td>
<td>2.27099 ns [2.25725, 2.28007] (0.95 confidence)</td>
<td>2.28497 ns [2.27103, 2.29821] (0.95 confidence)</td>
<td>0.993882</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;from&#95;rows</td>
<td>3.61651 ns [3.60889, 3.62578] (0.95 confidence)</td>
<td>3.68006 ns [3.67036, 3.68529] (0.95 confidence)</td>
<td>0.98273</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;identity</td>
<td>2.13798 ns [2.12896, 2.14302] (0.95 confidence)</td>
<td>2.17299 ns [2.17049, 2.17865] (0.95 confidence)</td>
<td>0.983886</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;set</td>
<td>3.8557 ns [3.84567, 3.86322] (0.95 confidence)</td>
<td>3.9397 ns [3.93306, 3.94521] (0.95 confidence)</td>
<td>0.978678</td>
</tr>
<tr>
<td>public&#95;d2/interval&#95;matrix&#95;zero</td>
<td>2.1447 ns [2.13994, 2.14951] (0.95 confidence)</td>
<td>2.1954 ns [2.19011, 2.2003] (0.95 confidence)</td>
<td>0.976908</td>
</tr>
<tr>
<td>public&#95;d2/is&#95;symmetric</td>
<td>0.761818 ns [0.761074, 0.762429] (0.95 confidence)</td>
<td>0.767699 ns [0.766779, 0.76845] (0.95 confidence)</td>
<td>0.992339</td>
</tr>
<tr>
<td>public&#95;d2/matrix&#95;access</td>
<td>4.79435 ns [4.78627, 4.80345] (0.95 confidence)</td>
<td>4.89564 ns [4.88798, 4.90277] (0.95 confidence)</td>
<td>0.97931</td>
</tr>
<tr>
<td>public&#95;d2/matrix&#95;from&#95;rows</td>
<td>2.93818 ns [2.93326, 2.94086] (0.95 confidence)</td>
<td>2.99304 ns [2.98637, 2.99734] (0.95 confidence)</td>
<td>0.981671</td>
</tr>
<tr>
<td>public&#95;d2/matrix&#95;identity</td>
<td>1.27748 ns [1.27327, 1.28322] (0.95 confidence)</td>
<td>1.30842 ns [1.30484, 1.31087] (0.95 confidence)</td>
<td>0.976351</td>
</tr>
<tr>
<td>public&#95;d2/matrix&#95;set</td>
<td>2.13376 ns [2.12321, 2.14786] (0.95 confidence)</td>
<td>2.1578 ns [2.15088, 2.16229] (0.95 confidence)</td>
<td>0.988861</td>
</tr>
<tr>
<td>public&#95;d2/matrix&#95;zero</td>
<td>1.17912 ns [1.17724, 1.18466] (0.95 confidence)</td>
<td>1.20899 ns [1.20506, 1.21104] (0.95 confidence)</td>
<td>0.975292</td>
</tr>
<tr>
<td>public&#95;d2/vector&#95;access</td>
<td>1.64303 ns [1.64091, 1.64514] (0.95 confidence)</td>
<td>1.66563 ns [1.66275, 1.66876] (0.95 confidence)</td>
<td>0.986433</td>
</tr>
<tr>
<td>public&#95;d2/vector&#95;new</td>
<td>1.96761 ns [1.96279, 1.97347] (0.95 confidence)</td>
<td>1.76165 ns [1.75567, 1.77235] (0.95 confidence)</td>
<td>1.11691</td>
</tr>
<tr>
<td>public&#95;d2/vector&#95;zero</td>
<td>0.602799 ns [0.601253, 0.603701] (0.95 confidence)</td>
<td>0.61155 ns [0.60965, 0.612843] (0.95 confidence)</td>
<td>0.98569</td>
</tr>
<tr>
<td>public&#95;d3/det&#95;bound</td>
<td>3.67422 ns [3.66832, 3.67938] (0.95 confidence)</td>
<td>3.70568 ns [3.70067, 3.71089] (0.95 confidence)</td>
<td>0.99151</td>
</tr>
<tr>
<td>public&#95;d3/det&#95;direct</td>
<td>1.29244 ns [1.29167, 1.29318] (0.95 confidence)</td>
<td>1.30279 ns [1.30093, 1.30469] (0.95 confidence)</td>
<td>0.992053</td>
</tr>
<tr>
<td>public&#95;d3/det&#95;with&#95;bound</td>
<td>3.66996 ns [3.66776, 3.67326] (0.95 confidence)</td>
<td>3.72518 ns [3.72202, 3.72946] (0.95 confidence)</td>
<td>0.985178</td>
</tr>
<tr>
<td>public&#95;d3/first&#95;asymmetry</td>
<td>3.53398 ns [3.50563, 3.57523] (0.95 confidence)</td>
<td>3.61288 ns [3.57889, 3.65383] (0.95 confidence)</td>
<td>0.978161</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;det</td>
<td>23.3439 ns [23.2256, 23.4651] (0.95 confidence)</td>
<td>22.7102 ns [22.5299, 22.9688] (0.95 confidence)</td>
<td>1.0279</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;access</td>
<td>10.0103 ns [9.99165, 10.0337] (0.95 confidence)</td>
<td>10.314 ns [10.294, 10.3287] (0.95 confidence)</td>
<td>0.970556</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;from&#95;matrix</td>
<td>5.12023 ns [5.11033, 5.13293] (0.95 confidence)</td>
<td>5.28562 ns [5.27803, 5.29839] (0.95 confidence)</td>
<td>0.968709</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;from&#95;points</td>
<td>6.31868 ns [6.3144, 6.32586] (0.95 confidence)</td>
<td>6.48973 ns [6.48225, 6.49827] (0.95 confidence)</td>
<td>0.973643</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;from&#95;rows</td>
<td>7.34075 ns [7.3236, 7.36019] (0.95 confidence)</td>
<td>7.59014 ns [7.57116, 7.60582] (0.95 confidence)</td>
<td>0.967143</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;identity</td>
<td>4.67122 ns [4.62911, 4.75314] (0.95 confidence)</td>
<td>4.47024 ns [4.46211, 4.48431] (0.95 confidence)</td>
<td>1.04496</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;set</td>
<td>7.81482 ns [7.79877, 7.83672] (0.95 confidence)</td>
<td>8.0944 ns [8.07371, 8.12466] (0.95 confidence)</td>
<td>0.96546</td>
</tr>
<tr>
<td>public&#95;d3/interval&#95;matrix&#95;zero</td>
<td>4.28669 ns [4.27133, 4.29701] (0.95 confidence)</td>
<td>4.49554 ns [4.48664, 4.50325] (0.95 confidence)</td>
<td>0.953543</td>
</tr>
<tr>
<td>public&#95;d3/is&#95;symmetric</td>
<td>1.3421 ns [1.34113, 1.34262] (0.95 confidence)</td>
<td>1.35068 ns [1.34828, 1.35276] (0.95 confidence)</td>
<td>0.993642</td>
</tr>
<tr>
<td>public&#95;d3/matrix&#95;access</td>
<td>6.76559 ns [6.73349, 6.79272] (0.95 confidence)</td>
<td>6.89216 ns [6.86763, 6.91185] (0.95 confidence)</td>
<td>0.981636</td>
</tr>
<tr>
<td>public&#95;d3/matrix&#95;from&#95;rows</td>
<td>5.79534 ns [5.77799, 5.80943] (0.95 confidence)</td>
<td>5.84901 ns [5.82594, 5.87455] (0.95 confidence)</td>
<td>0.990823</td>
</tr>
<tr>
<td>public&#95;d3/matrix&#95;identity</td>
<td>2.75548 ns [2.74371, 2.76573] (0.95 confidence)</td>
<td>2.79092 ns [2.78592, 2.8027] (0.95 confidence)</td>
<td>0.987303</td>
</tr>
<tr>
<td>public&#95;d3/matrix&#95;set</td>
<td>4.27928 ns [4.26063, 4.29164] (0.95 confidence)</td>
<td>4.34506 ns [4.32808, 4.35437] (0.95 confidence)</td>
<td>0.984861</td>
</tr>
<tr>
<td>public&#95;d3/matrix&#95;zero</td>
<td>2.45734 ns [2.45415, 2.46974] (0.95 confidence)</td>
<td>2.51635 ns [2.50813, 2.52054] (0.95 confidence)</td>
<td>0.976551</td>
</tr>
<tr>
<td>public&#95;d3/vector&#95;access</td>
<td>2.08637 ns [2.0798, 2.0907] (0.95 confidence)</td>
<td>2.14876 ns [2.1449, 2.15169] (0.95 confidence)</td>
<td>0.970965</td>
</tr>
<tr>
<td>public&#95;d3/vector&#95;new</td>
<td>3.63155 ns [3.61943, 3.64055] (0.95 confidence)</td>
<td>3.68724 ns [3.67299, 3.69597] (0.95 confidence)</td>
<td>0.984898</td>
</tr>
<tr>
<td>public&#95;d3/vector&#95;zero</td>
<td>0.954174 ns [0.95187, 0.956787] (0.95 confidence)</td>
<td>0.986501 ns [0.984875, 0.987578] (0.95 confidence)</td>
<td>0.96723</td>
</tr>
<tr>
<td>public&#95;d4/det&#95;bound</td>
<td>6.37762 ns [6.37032, 6.38708] (0.95 confidence)</td>
<td>6.62964 ns [6.52631, 6.68941] (0.95 confidence)</td>
<td>0.961987</td>
</tr>
<tr>
<td>public&#95;d4/det&#95;direct</td>
<td>2.47364 ns [2.47316, 2.47378] (0.95 confidence)</td>
<td>2.50505 ns [2.50237, 2.50893] (0.95 confidence)</td>
<td>0.987463</td>
</tr>
<tr>
<td>public&#95;d4/det&#95;with&#95;bound</td>
<td>6.41158 ns [6.39391, 6.42253] (0.95 confidence)</td>
<td>6.66738 ns [6.55005, 6.71849] (0.95 confidence)</td>
<td>0.961633</td>
</tr>
<tr>
<td>public&#95;d4/first&#95;asymmetry</td>
<td>3.64802 ns [3.62173, 3.68528] (0.95 confidence)</td>
<td>3.70142 ns [3.69141, 3.72898] (0.95 confidence)</td>
<td>0.985572</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;det</td>
<td>48.6644 ns [48.2495, 49.2578] (0.95 confidence)</td>
<td>53.7948 ns [53.6482, 54.0948] (0.95 confidence)</td>
<td>0.90463</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;access</td>
<td>13.4584 ns [13.4196, 13.6411] (0.95 confidence)</td>
<td>14.0886 ns [13.976, 14.145] (0.95 confidence)</td>
<td>0.955269</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;from&#95;matrix</td>
<td>7.80561 ns [7.79517, 7.82008] (0.95 confidence)</td>
<td>8.17428 ns [8.16722, 8.18694] (0.95 confidence)</td>
<td>0.954898</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;from&#95;points</td>
<td>10.8415 ns [10.8158, 10.8572] (0.95 confidence)</td>
<td>11.1681 ns [11.1376, 11.1881] (0.95 confidence)</td>
<td>0.970761</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;from&#95;rows</td>
<td>11.0311 ns [11.0117, 11.0479] (0.95 confidence)</td>
<td>11.6107 ns [11.5881, 11.6295] (0.95 confidence)</td>
<td>0.950082</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;identity</td>
<td>6.41912 ns [6.40275, 6.42916] (0.95 confidence)</td>
<td>6.88525 ns [6.8654, 6.90667] (0.95 confidence)</td>
<td>0.932301</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;set</td>
<td>11.2494 ns [11.2374, 11.2616] (0.95 confidence)</td>
<td>11.9395 ns [11.9231, 11.9548] (0.95 confidence)</td>
<td>0.942206</td>
</tr>
<tr>
<td>public&#95;d4/interval&#95;matrix&#95;zero</td>
<td>6.4747 ns [6.46616, 6.49034] (0.95 confidence)</td>
<td>7.01475 ns [6.96243, 7.20333] (0.95 confidence)</td>
<td>0.923013</td>
</tr>
<tr>
<td>public&#95;d4/is&#95;symmetric</td>
<td>2.10281 ns [2.10263, 2.10302] (0.95 confidence)</td>
<td>2.12996 ns [2.12782, 2.13131] (0.95 confidence)</td>
<td>0.987253</td>
</tr>
<tr>
<td>public&#95;d4/matrix&#95;access</td>
<td>9.3151 ns [9.29236, 9.33634] (0.95 confidence)</td>
<td>9.7251 ns [9.71024, 9.75224] (0.95 confidence)</td>
<td>0.957842</td>
</tr>
<tr>
<td>public&#95;d4/matrix&#95;from&#95;rows</td>
<td>10.2734 ns [10.2621, 10.2841] (0.95 confidence)</td>
<td>10.543 ns [10.5326, 10.5631] (0.95 confidence)</td>
<td>0.974431</td>
</tr>
<tr>
<td>public&#95;d4/matrix&#95;identity</td>
<td>4.10883 ns [4.09324, 4.12258] (0.95 confidence)</td>
<td>4.27582 ns [4.26675, 4.28489] (0.95 confidence)</td>
<td>0.960946</td>
</tr>
<tr>
<td>public&#95;d4/matrix&#95;set</td>
<td>7.33761 ns [7.32816, 7.3541] (0.95 confidence)</td>
<td>7.64868 ns [7.63992, 7.67442] (0.95 confidence)</td>
<td>0.95933</td>
</tr>
<tr>
<td>public&#95;d4/matrix&#95;zero</td>
<td>3.84687 ns [3.83638, 3.85274] (0.95 confidence)</td>
<td>3.99165 ns [3.9812, 4.00092] (0.95 confidence)</td>
<td>0.963728</td>
</tr>
<tr>
<td>public&#95;d4/vector&#95;access</td>
<td>2.59668 ns [2.59187, 2.60151] (0.95 confidence)</td>
<td>2.70282 ns [2.69674, 2.70756] (0.95 confidence)</td>
<td>0.960729</td>
</tr>
<tr>
<td>public&#95;d4/vector&#95;new</td>
<td>2.80767 ns [2.80121, 2.82598] (0.95 confidence)</td>
<td>2.87781 ns [2.86117, 2.89678] (0.95 confidence)</td>
<td>0.975625</td>
</tr>
<tr>
<td>public&#95;d4/vector&#95;zero</td>
<td>1.15445 ns [1.15152, 1.15645] (0.95 confidence)</td>
<td>1.21223 ns [1.20848, 1.21526] (0.95 confidence)</td>
<td>0.952338</td>
</tr>
<tr>
<td>public&#95;d5/det&#95;bound</td>
<td>1.31001 ns [1.30991, 1.31016] (0.95 confidence)</td>
<td>1.41115 ns [1.40774, 1.41368] (0.95 confidence)</td>
<td>0.928327</td>
</tr>
<tr>
<td>public&#95;d5/det&#95;direct</td>
<td>1.31041 ns [1.31019, 1.31076] (0.95 confidence)</td>
<td>1.40497 ns [1.39824, 1.4084] (0.95 confidence)</td>
<td>0.93269</td>
</tr>
<tr>
<td>public&#95;d5/det&#95;with&#95;bound</td>
<td>0.843094 ns [0.840616, 0.845509] (0.95 confidence)</td>
<td>0.898635 ns [0.895712, 0.900988] (0.95 confidence)</td>
<td>0.938194</td>
</tr>
<tr>
<td>public&#95;d5/first&#95;asymmetry</td>
<td>3.79534 ns [3.79111, 3.80453] (0.95 confidence)</td>
<td>3.86037 ns [3.85161, 3.86491] (0.95 confidence)</td>
<td>0.983156</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;det</td>
<td>122.694 ns [122.165, 123.268] (0.95 confidence)</td>
<td>131.622 ns [130.696, 132.313] (0.95 confidence)</td>
<td>0.932167</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;access</td>
<td>15.3045 ns [15.1818, 15.6116] (0.95 confidence)</td>
<td>15.3425 ns [15.2142, 15.6131] (0.95 confidence)</td>
<td>0.997526</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;from&#95;matrix</td>
<td>11.2064 ns [11.1129, 11.2282] (0.95 confidence)</td>
<td>10.9004 ns [10.8299, 10.9345] (0.95 confidence)</td>
<td>1.02807</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;from&#95;points</td>
<td>19.3 ns [19.2874, 19.3172] (0.95 confidence)</td>
<td>19.027 ns [19.0018, 19.0542] (0.95 confidence)</td>
<td>1.01435</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;from&#95;rows</td>
<td>16.1425 ns [16.0231, 16.1651] (0.95 confidence)</td>
<td>15.5384 ns [15.1745, 15.8695] (0.95 confidence)</td>
<td>1.03888</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;identity</td>
<td>8.99665 ns [8.39897, 9.04777] (0.95 confidence)</td>
<td>8.16754 ns [8.16014, 8.17734] (0.95 confidence)</td>
<td>1.10151</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;set</td>
<td>13.3761 ns [13.0966, 14.0127] (0.95 confidence)</td>
<td>13.6667 ns [13.2054, 13.9465] (0.95 confidence)</td>
<td>0.978739</td>
</tr>
<tr>
<td>public&#95;d5/interval&#95;matrix&#95;zero</td>
<td>9.90295 ns [9.89906, 9.90675] (0.95 confidence)</td>
<td>9.64802 ns [9.63295, 9.66546] (0.95 confidence)</td>
<td>1.02642</td>
</tr>
<tr>
<td>public&#95;d5/is&#95;symmetric</td>
<td>3.68584 ns [3.68557, 3.68609] (0.95 confidence)</td>
<td>3.73198 ns [3.72887, 3.73475] (0.95 confidence)</td>
<td>0.987636</td>
</tr>
<tr>
<td>public&#95;d5/matrix&#95;access</td>
<td>11.9483 ns [11.9187, 12.1921] (0.95 confidence)</td>
<td>12.6112 ns [12.4791, 12.7039] (0.95 confidence)</td>
<td>0.947439</td>
</tr>
<tr>
<td>public&#95;d5/matrix&#95;from&#95;rows</td>
<td>11.7525 ns [11.7322, 11.7628] (0.95 confidence)</td>
<td>12.0078 ns [11.9929, 12.0211] (0.95 confidence)</td>
<td>0.97874</td>
</tr>
<tr>
<td>public&#95;d5/matrix&#95;identity</td>
<td>6.66273 ns [6.62334, 6.67512] (0.95 confidence)</td>
<td>7.11069 ns [7.08195, 7.14026] (0.95 confidence)</td>
<td>0.937001</td>
</tr>
<tr>
<td>public&#95;d5/matrix&#95;set</td>
<td>10.3218 ns [10.3004, 10.3517] (0.95 confidence)</td>
<td>10.9117 ns [10.8847, 10.9341] (0.95 confidence)</td>
<td>0.945937</td>
</tr>
<tr>
<td>public&#95;d5/matrix&#95;zero</td>
<td>5.60791 ns [5.59867, 5.62688] (0.95 confidence)</td>
<td>5.99943 ns [5.98744, 6.00687] (0.95 confidence)</td>
<td>0.934742</td>
</tr>
<tr>
<td>public&#95;d5/vector&#95;access</td>
<td>3.04147 ns [3.03412, 3.04814] (0.95 confidence)</td>
<td>3.19567 ns [3.18844, 3.20136] (0.95 confidence)</td>
<td>0.951749</td>
</tr>
<tr>
<td>public&#95;d5/vector&#95;new</td>
<td>4.84713 ns [4.83327, 4.85243] (0.95 confidence)</td>
<td>4.92635 ns [4.91247, 4.93489] (0.95 confidence)</td>
<td>0.983917</td>
</tr>
<tr>
<td>public&#95;d5/vector&#95;zero</td>
<td>1.49217 ns [1.47789, 1.50115] (0.95 confidence)</td>
<td>1.61857 ns [1.60729, 1.62493] (0.95 confidence)</td>
<td>0.921905</td>
</tr>
<tr>
<td>public&#95;scalar/interval&#95;access</td>
<td>1.61245 ns [1.61014, 1.61437] (0.95 confidence)</td>
<td>1.61824 ns [1.6171, 1.62069] (0.95 confidence)</td>
<td>0.996421</td>
</tr>
<tr>
<td>public&#95;scalar/interval&#95;negate</td>
<td>0.738279 ns [0.737671, 0.738996] (0.95 confidence)</td>
<td>0.742832 ns [0.741921, 0.74442] (0.95 confidence)</td>
<td>0.993871</td>
</tr>
<tr>
<td>public&#95;scalar/interval&#95;new</td>
<td>0.896569 ns [0.896044, 0.896978] (0.95 confidence)</td>
<td>0.900083 ns [0.89945, 0.900638] (0.95 confidence)</td>
<td>0.996097</td>
</tr>
<tr>
<td>public&#95;scalar/interval&#95;point</td>
<td>0.58625 ns [0.585719, 0.587218] (0.95 confidence)</td>
<td>0.587801 ns [0.587558, 0.588097] (0.95 confidence)</td>
<td>0.99736</td>
</tr>
<tr>
<td>public&#95;scalar/interval&#95;subtraction</td>
<td>1.88641 ns [1.8832, 1.89039] (0.95 confidence)</td>
<td>1.8947 ns [1.893, 1.89768] (0.95 confidence)</td>
<td>0.995623</td>
</tr>
<tr>
<td>public&#95;scalar/tolerance</td>
<td>0.990502 ns [0.989109, 0.992242] (0.95 confidence)</td>
<td>1.00527 ns [1.00353, 1.00713] (0.95 confidence)</td>
<td>0.985311</td>
</tr>
</tbody>
</table>

## core-exact

Semantic contract: la-stack-core-exact-v1. Collection session: `6886cf37ab85432a8fe07f82593c9754`.
Collection interval (Unix ns): 1794906510083407000–1794918235914777000.
Acquisition: newly collected.

- baseline: available — inside declared release interval; release v0.5.3; revision `7165ed562f2e4611d201e574fc1f600a8ad8b33e`.
  Source SHA-256: `073502f95ae58344aa776b54303ae24d1e6cf2ce60d1437d950bcbe64d0b4dbf`; harness SHA-256:
  `ba63ffff14e6124d5aca8df516e6836a606ad10809bbe3b1ab13dd5fd7ceedbf`; scientific gate: passed.
- current: available — inside declared release interval; release v0.5.4; revision `531f75eb58eb8ac810ad5d89bf1d3e66a301e84e`.
  Source SHA-256: `50d72eb42299ad56e00e6838b34b10824c560a3eb5ae9c9a05e48859fd3cfe46`; harness SHA-256:
  `ba63ffff14e6124d5aca8df516e6836a606ad10809bbe3b1ab13dd5fd7ceedbf`; scientific gate: passed.

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
<td>168.803 ns [168.417, 169.213] (0.95 confidence)</td>
<td>163.719 ns [163.228, 164.237] (0.95 confidence)</td>
<td>1.03105</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d2/rounded&#95;result</td>
<td>170.81 ns [170.259, 171.369] (0.95 confidence)</td>
<td>164.702 ns [164.287, 165.124] (0.95 confidence)</td>
<td>1.03709</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d2/strict&#95;array</td>
<td>173.924 ns [173.185, 174.686] (0.95 confidence)</td>
<td>167.14 ns [166.719, 167.648] (0.95 confidence)</td>
<td>1.04059</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d2/strict&#95;result</td>
<td>173.076 ns [172.504, 173.659] (0.95 confidence)</td>
<td>169.801 ns [169.289, 170.311] (0.95 confidence)</td>
<td>1.01929</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d3/rounded&#95;array</td>
<td>177.123 ns [176.519, 177.748] (0.95 confidence)</td>
<td>169.636 ns [169.152, 170.124] (0.95 confidence)</td>
<td>1.04414</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d3/rounded&#95;result</td>
<td>176.915 ns [176.386, 177.454] (0.95 confidence)</td>
<td>173.264 ns [172.712, 173.826] (0.95 confidence)</td>
<td>1.02107</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d3/strict&#95;array</td>
<td>177.238 ns [176.584, 177.916] (0.95 confidence)</td>
<td>170.246 ns [169.908, 170.587] (0.95 confidence)</td>
<td>1.04107</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d3/strict&#95;result</td>
<td>176.378 ns [175.847, 176.949] (0.95 confidence)</td>
<td>169.754 ns [169.247, 170.275] (0.95 confidence)</td>
<td>1.03902</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d4/rounded&#95;array</td>
<td>178.442 ns [177.869, 179.017] (0.95 confidence)</td>
<td>175.088 ns [174.696, 175.487] (0.95 confidence)</td>
<td>1.01916</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d4/rounded&#95;result</td>
<td>177.569 ns [177.136, 178.014] (0.95 confidence)</td>
<td>177.104 ns [176.635, 177.587] (0.95 confidence)</td>
<td>1.00263</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d4/strict&#95;array</td>
<td>177.503 ns [176.994, 178.045] (0.95 confidence)</td>
<td>174.521 ns [174.118, 174.932] (0.95 confidence)</td>
<td>1.01709</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d4/strict&#95;result</td>
<td>176.776 ns [176.345, 177.237] (0.95 confidence)</td>
<td>173.974 ns [173.462, 174.54] (0.95 confidence)</td>
<td>1.01611</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d5/rounded&#95;array</td>
<td>185.628 ns [185.124, 186.13] (0.95 confidence)</td>
<td>181.719 ns [181.183, 182.246] (0.95 confidence)</td>
<td>1.02151</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d5/rounded&#95;result</td>
<td>182.932 ns [182.552, 183.315] (0.95 confidence)</td>
<td>181.768 ns [181.309, 182.227] (0.95 confidence)</td>
<td>1.0064</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d5/strict&#95;array</td>
<td>182.563 ns [182.096, 183.042] (0.95 confidence)</td>
<td>181.807 ns [181.203, 182.443] (0.95 confidence)</td>
<td>1.00416</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d5/strict&#95;result</td>
<td>180.087 ns [179.43, 180.763] (0.95 confidence)</td>
<td>178.643 ns [178.268, 179.029] (0.95 confidence)</td>
<td>1.00809</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d2/rounded&#95;array</td>
<td>11.7657 ns [11.7425, 11.7924] (0.95 confidence)</td>
<td>11.5404 ns [11.5325, 11.5479] (0.95 confidence)</td>
<td>1.01952</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d2/rounded&#95;result</td>
<td>11.7746 ns [11.7513, 11.8022] (0.95 confidence)</td>
<td>11.5334 ns [11.5229, 11.5429] (0.95 confidence)</td>
<td>1.02092</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d2/strict&#95;array</td>
<td>8.76726 ns [8.71211, 8.82606] (0.95 confidence)</td>
<td>8.69433 ns [8.60533, 8.78973] (0.95 confidence)</td>
<td>1.00839</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d2/strict&#95;result</td>
<td>8.20456 ns [8.16399, 8.24655] (0.95 confidence)</td>
<td>8.01097 ns [7.97088, 8.05231] (0.95 confidence)</td>
<td>1.02417</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d3/rounded&#95;array</td>
<td>18.1443 ns [18.0967, 18.2045] (0.95 confidence)</td>
<td>17.8673 ns [17.847, 17.8867] (0.95 confidence)</td>
<td>1.01551</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d3/rounded&#95;result</td>
<td>18.058 ns [18.0426, 18.0752] (0.95 confidence)</td>
<td>17.8628 ns [17.8411, 17.8854] (0.95 confidence)</td>
<td>1.01093</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d3/strict&#95;array</td>
<td>12.9637 ns [12.8537, 13.0777] (0.95 confidence)</td>
<td>12.2597 ns [12.1827, 12.3477] (0.95 confidence)</td>
<td>1.05742</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d3/strict&#95;result</td>
<td>11.4914 ns [11.4143, 11.5723] (0.95 confidence)</td>
<td>11.8049 ns [11.7354, 11.875] (0.95 confidence)</td>
<td>0.973444</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d4/rounded&#95;array</td>
<td>23.7375 ns [23.719, 23.7626] (0.95 confidence)</td>
<td>23.6618 ns [23.6353, 23.6866] (0.95 confidence)</td>
<td>1.0032</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d4/rounded&#95;result</td>
<td>23.7475 ns [23.7339, 23.7632] (0.95 confidence)</td>
<td>23.6553 ns [23.6159, 23.6991] (0.95 confidence)</td>
<td>1.0039</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d4/strict&#95;array</td>
<td>16.4923 ns [16.4289, 16.559] (0.95 confidence)</td>
<td>17.0097 ns [16.8696, 17.1542] (0.95 confidence)</td>
<td>0.969584</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d4/strict&#95;result</td>
<td>15.3745 ns [15.2943, 15.4593] (0.95 confidence)</td>
<td>15.6357 ns [15.5416, 15.7323] (0.95 confidence)</td>
<td>0.983296</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d5/rounded&#95;array</td>
<td>28.4071 ns [28.394, 28.4231] (0.95 confidence)</td>
<td>28.3566 ns [28.3306, 28.3812] (0.95 confidence)</td>
<td>1.00178</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d5/rounded&#95;result</td>
<td>28.4319 ns [28.4022, 28.4575] (0.95 confidence)</td>
<td>28.3337 ns [28.2917, 28.3696] (0.95 confidence)</td>
<td>1.00347</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d5/strict&#95;array</td>
<td>21.2019 ns [21.0126, 21.3966] (0.95 confidence)</td>
<td>20.6092 ns [20.46, 20.7771] (0.95 confidence)</td>
<td>1.02876</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d5/strict&#95;result</td>
<td>18.9488 ns [18.8285, 19.0817] (0.95 confidence)</td>
<td>19.4563 ns [19.3223, 19.5886] (0.95 confidence)</td>
<td>0.973915</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d2/rounded&#95;array</td>
<td>181.158 ns [180.776, 181.576] (0.95 confidence)</td>
<td>169.008 ns [168.796, 169.213] (0.95 confidence)</td>
<td>1.07189</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d2/rounded&#95;result</td>
<td>180.536 ns [180.174, 180.93] (0.95 confidence)</td>
<td>169.378 ns [168.979, 169.73] (0.95 confidence)</td>
<td>1.06588</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d2/strict&#95;array</td>
<td>14.1148 ns [13.9809, 14.2509] (0.95 confidence)</td>
<td>13.1051 ns [13.0169, 13.1945] (0.95 confidence)</td>
<td>1.07705</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d2/strict&#95;result</td>
<td>12.7275 ns [12.6609, 12.7965] (0.95 confidence)</td>
<td>12.2806 ns [12.2075, 12.3549] (0.95 confidence)</td>
<td>1.03639</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d3/rounded&#95;array</td>
<td>179.516 ns [179.224, 179.831] (0.95 confidence)</td>
<td>176.043 ns [175.766, 176.319] (0.95 confidence)</td>
<td>1.01972</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d3/rounded&#95;result</td>
<td>180.308 ns [179.942, 180.687] (0.95 confidence)</td>
<td>175.371 ns [175.12, 175.619] (0.95 confidence)</td>
<td>1.02815</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d3/strict&#95;array</td>
<td>18.999 ns [18.8439, 19.1487] (0.95 confidence)</td>
<td>18.2251 ns [18.1004, 18.3493] (0.95 confidence)</td>
<td>1.04246</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d3/strict&#95;result</td>
<td>17.3245 ns [17.2064, 17.4438] (0.95 confidence)</td>
<td>16.5734 ns [16.4827, 16.6669] (0.95 confidence)</td>
<td>1.04532</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d4/rounded&#95;array</td>
<td>191.098 ns [190.869, 191.343] (0.95 confidence)</td>
<td>183.429 ns [183.03, 183.845] (0.95 confidence)</td>
<td>1.04181</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d4/rounded&#95;result</td>
<td>191.515 ns [191.264, 191.777] (0.95 confidence)</td>
<td>181.939 ns [181.619, 182.294] (0.95 confidence)</td>
<td>1.05263</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d4/strict&#95;array</td>
<td>23.6651 ns [23.5398, 23.79] (0.95 confidence)</td>
<td>22.8766 ns [22.7381, 23.0157] (0.95 confidence)</td>
<td>1.03446</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d4/strict&#95;result</td>
<td>20.9173 ns [20.7755, 21.0647] (0.95 confidence)</td>
<td>21.1753 ns [21.0493, 21.3001] (0.95 confidence)</td>
<td>0.987815</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d5/rounded&#95;array</td>
<td>188.868 ns [188.403, 189.428] (0.95 confidence)</td>
<td>187.567 ns [187.197, 187.941] (0.95 confidence)</td>
<td>1.00694</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d5/rounded&#95;result</td>
<td>188.695 ns [188.422, 188.977] (0.95 confidence)</td>
<td>186.38 ns [186.117, 186.637] (0.95 confidence)</td>
<td>1.01242</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d5/strict&#95;array</td>
<td>27.8753 ns [27.753, 27.9988] (0.95 confidence)</td>
<td>27.0024 ns [26.88, 27.1337] (0.95 confidence)</td>
<td>1.03233</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d5/strict&#95;result</td>
<td>24.854 ns [24.6903, 25.0235] (0.95 confidence)</td>
<td>24.2753 ns [24.0905, 24.4646] (0.95 confidence)</td>
<td>1.02384</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d2/rounded&#95;array</td>
<td>180.767 ns [180.4, 181.16] (0.95 confidence)</td>
<td>169.959 ns [169.643, 170.273] (0.95 confidence)</td>
<td>1.06359</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d2/rounded&#95;result</td>
<td>181.105 ns [180.609, 181.725] (0.95 confidence)</td>
<td>168.417 ns [168.174, 168.662] (0.95 confidence)</td>
<td>1.07534</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d2/strict&#95;array</td>
<td>13.5995 ns [13.5157, 13.6858] (0.95 confidence)</td>
<td>13.2183 ns [13.0816, 13.3623] (0.95 confidence)</td>
<td>1.02884</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d2/strict&#95;result</td>
<td>12.6307 ns [12.5619, 12.7028] (0.95 confidence)</td>
<td>12.3413 ns [12.2502, 12.4335] (0.95 confidence)</td>
<td>1.02345</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d3/rounded&#95;array</td>
<td>180.416 ns [180.002, 180.864] (0.95 confidence)</td>
<td>174.609 ns [174.309, 174.887] (0.95 confidence)</td>
<td>1.03326</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d3/rounded&#95;result</td>
<td>179.046 ns [178.751, 179.346] (0.95 confidence)</td>
<td>176.103 ns [175.853, 176.356] (0.95 confidence)</td>
<td>1.01671</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d3/strict&#95;array</td>
<td>17.9692 ns [17.857, 18.0795] (0.95 confidence)</td>
<td>17.215 ns [17.101, 17.3272] (0.95 confidence)</td>
<td>1.04381</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d3/strict&#95;result</td>
<td>16.7368 ns [16.6435, 16.829] (0.95 confidence)</td>
<td>16.3229 ns [16.213, 16.4362] (0.95 confidence)</td>
<td>1.02536</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d4/rounded&#95;array</td>
<td>191.358 ns [191.089, 191.635] (0.95 confidence)</td>
<td>180.874 ns [180.637, 181.108] (0.95 confidence)</td>
<td>1.05796</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d4/rounded&#95;result</td>
<td>192.177 ns [191.871, 192.531] (0.95 confidence)</td>
<td>180.94 ns [180.702, 181.173] (0.95 confidence)</td>
<td>1.0621</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d4/strict&#95;array</td>
<td>22.0584 ns [21.9789, 22.1383] (0.95 confidence)</td>
<td>21.4644 ns [21.3843, 21.5445] (0.95 confidence)</td>
<td>1.02767</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d4/strict&#95;result</td>
<td>21.2461 ns [21.1006, 21.4035] (0.95 confidence)</td>
<td>20.4725 ns [20.3311, 20.6103] (0.95 confidence)</td>
<td>1.03779</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d5/rounded&#95;array</td>
<td>188.229 ns [187.981, 188.48] (0.95 confidence)</td>
<td>185.354 ns [185.041, 185.663] (0.95 confidence)</td>
<td>1.01551</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d5/rounded&#95;result</td>
<td>188.476 ns [188.125, 188.853] (0.95 confidence)</td>
<td>185.719 ns [185.422, 186.009] (0.95 confidence)</td>
<td>1.01485</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d5/strict&#95;array</td>
<td>26.2719 ns [26.1381, 26.4113] (0.95 confidence)</td>
<td>25.1257 ns [24.9683, 25.2943] (0.95 confidence)</td>
<td>1.04562</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d5/strict&#95;result</td>
<td>24.0639 ns [23.9613, 24.1726] (0.95 confidence)</td>
<td>24.1301 ns [24.0222, 24.2343] (0.95 confidence)</td>
<td>0.997257</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d2/rounded&#95;array</td>
<td>11.7976 ns [11.7794, 11.8178] (0.95 confidence)</td>
<td>11.5748 ns [11.5615, 11.5924] (0.95 confidence)</td>
<td>1.01925</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d2/rounded&#95;result</td>
<td>11.7242 ns [11.7079, 11.7432] (0.95 confidence)</td>
<td>11.5667 ns [11.5548, 11.578] (0.95 confidence)</td>
<td>1.01362</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d2/strict&#95;array</td>
<td>56.472 ns [56.355, 56.5957] (0.95 confidence)</td>
<td>54.8654 ns [54.7306, 54.9999] (0.95 confidence)</td>
<td>1.02928</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d2/strict&#95;result</td>
<td>11.9506 ns [11.9232, 11.9787] (0.95 confidence)</td>
<td>11.8611 ns [11.8356, 11.8866] (0.95 confidence)</td>
<td>1.00755</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d3/rounded&#95;array</td>
<td>18.1596 ns [18.1304, 18.192] (0.95 confidence)</td>
<td>17.9625 ns [17.9454, 17.9797] (0.95 confidence)</td>
<td>1.01097</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d3/rounded&#95;result</td>
<td>18.1689 ns [18.1428, 18.2027] (0.95 confidence)</td>
<td>17.9306 ns [17.9109, 17.9494] (0.95 confidence)</td>
<td>1.01329</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d3/strict&#95;array</td>
<td>60.6232 ns [60.4784, 60.772] (0.95 confidence)</td>
<td>59.2829 ns [59.1386, 59.4448] (0.95 confidence)</td>
<td>1.02261</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d3/strict&#95;result</td>
<td>15.5596 ns [15.5117, 15.6062] (0.95 confidence)</td>
<td>15.496 ns [15.4672, 15.5244] (0.95 confidence)</td>
<td>1.0041</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d4/rounded&#95;array</td>
<td>23.8265 ns [23.8101, 23.8465] (0.95 confidence)</td>
<td>23.7253 ns [23.6781, 23.7646] (0.95 confidence)</td>
<td>1.00427</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d4/rounded&#95;result</td>
<td>23.8463 ns [23.8257, 23.8786] (0.95 confidence)</td>
<td>23.7634 ns [23.736, 23.79] (0.95 confidence)</td>
<td>1.00349</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d4/strict&#95;array</td>
<td>63.2999 ns [63.0987, 63.5059] (0.95 confidence)</td>
<td>62.5591 ns [62.3174, 62.8086] (0.95 confidence)</td>
<td>1.01184</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d4/strict&#95;result</td>
<td>19.5322 ns [19.4925, 19.5736] (0.95 confidence)</td>
<td>19.4448 ns [19.404, 19.4852] (0.95 confidence)</td>
<td>1.0045</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d5/rounded&#95;array</td>
<td>28.3417 ns [28.3305, 28.3535] (0.95 confidence)</td>
<td>28.3222 ns [28.2825, 28.3591] (0.95 confidence)</td>
<td>1.00069</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d5/rounded&#95;result</td>
<td>28.366 ns [28.3543, 28.3794] (0.95 confidence)</td>
<td>28.2855 ns [28.2585, 28.3104] (0.95 confidence)</td>
<td>1.00285</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d5/strict&#95;array</td>
<td>66.1842 ns [65.9856, 66.397] (0.95 confidence)</td>
<td>65.0798 ns [64.8529, 65.2996] (0.95 confidence)</td>
<td>1.01697</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d5/strict&#95;result</td>
<td>22.4628 ns [22.3609, 22.5675] (0.95 confidence)</td>
<td>22.7212 ns [22.6435, 22.8001] (0.95 confidence)</td>
<td>0.98863</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d2/rounded&#95;array</td>
<td>171.024 ns [170.508, 171.536] (0.95 confidence)</td>
<td>165.682 ns [165.365, 166.01] (0.95 confidence)</td>
<td>1.03224</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d2/rounded&#95;result</td>
<td>171.785 ns [171.177, 172.418] (0.95 confidence)</td>
<td>167.184 ns [166.811, 167.541] (0.95 confidence)</td>
<td>1.02752</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d2/strict&#95;array</td>
<td>180.037 ns [179.497, 180.595] (0.95 confidence)</td>
<td>173.633 ns [173.285, 173.985] (0.95 confidence)</td>
<td>1.03688</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d2/strict&#95;result</td>
<td>180.514 ns [179.986, 181.06] (0.95 confidence)</td>
<td>172.244 ns [171.808, 172.678] (0.95 confidence)</td>
<td>1.04801</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d3/rounded&#95;array</td>
<td>178.036 ns [177.455, 178.652] (0.95 confidence)</td>
<td>171.95 ns [171.499, 172.407] (0.95 confidence)</td>
<td>1.03539</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d3/rounded&#95;result</td>
<td>178.675 ns [178.052, 179.356] (0.95 confidence)</td>
<td>171.068 ns [170.772, 171.369] (0.95 confidence)</td>
<td>1.04446</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d3/strict&#95;array</td>
<td>185.272 ns [184.627, 185.946] (0.95 confidence)</td>
<td>178.68 ns [178.175, 179.241] (0.95 confidence)</td>
<td>1.03689</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d3/strict&#95;result</td>
<td>182.382 ns [181.816, 182.975] (0.95 confidence)</td>
<td>176.486 ns [176.016, 176.924] (0.95 confidence)</td>
<td>1.03341</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d4/rounded&#95;array</td>
<td>181.936 ns [181.369, 182.538] (0.95 confidence)</td>
<td>178.452 ns [178.05, 178.866] (0.95 confidence)</td>
<td>1.01952</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d4/rounded&#95;result</td>
<td>180.369 ns [179.833, 180.909] (0.95 confidence)</td>
<td>178.031 ns [177.631, 178.425] (0.95 confidence)</td>
<td>1.01313</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d4/strict&#95;array</td>
<td>184.581 ns [184.091, 185.084] (0.95 confidence)</td>
<td>182.124 ns [181.575, 182.81] (0.95 confidence)</td>
<td>1.01349</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d4/strict&#95;result</td>
<td>182.345 ns [181.817, 182.892] (0.95 confidence)</td>
<td>179.188 ns [178.839, 179.542] (0.95 confidence)</td>
<td>1.01762</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d5/rounded&#95;array</td>
<td>187.652 ns [187.099, 188.192] (0.95 confidence)</td>
<td>183.574 ns [183.277, 183.865] (0.95 confidence)</td>
<td>1.02221</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d5/rounded&#95;result</td>
<td>187.704 ns [187.205, 188.215] (0.95 confidence)</td>
<td>183.269 ns [182.926, 183.614] (0.95 confidence)</td>
<td>1.0242</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d5/strict&#95;array</td>
<td>190.901 ns [190.273, 191.537] (0.95 confidence)</td>
<td>184.247 ns [183.811, 184.684] (0.95 confidence)</td>
<td>1.03612</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d5/strict&#95;result</td>
<td>184.757 ns [184.4, 185.139] (0.95 confidence)</td>
<td>182.672 ns [182.391, 182.963] (0.95 confidence)</td>
<td>1.01141</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d2/rounded&#95;array</td>
<td>204.741 ns [204.114, 205.42] (0.95 confidence)</td>
<td>198.36 ns [197.814, 198.918] (0.95 confidence)</td>
<td>1.03217</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d2/rounded&#95;result</td>
<td>204.407 ns [204.017, 204.802] (0.95 confidence)</td>
<td>195.848 ns [195.538, 196.156] (0.95 confidence)</td>
<td>1.0437</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d2/strict&#95;array</td>
<td>13504.9 ns [13469.7, 13540.2] (0.95 confidence)</td>
<td>13171.3 ns [13146.9, 13195.9] (0.95 confidence)</td>
<td>1.02533</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d2/strict&#95;result</td>
<td>205.317 ns [204.809, 205.839] (0.95 confidence)</td>
<td>198.314 ns [197.781, 198.872] (0.95 confidence)</td>
<td>1.03532</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d3/rounded&#95;array</td>
<td>209.741 ns [209.263, 210.241] (0.95 confidence)</td>
<td>201.82 ns [201.544, 202.081] (0.95 confidence)</td>
<td>1.03925</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d3/rounded&#95;result</td>
<td>210.83 ns [210.139, 211.72] (0.95 confidence)</td>
<td>205.529 ns [205.108, 205.944] (0.95 confidence)</td>
<td>1.02579</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d3/strict&#95;array</td>
<td>13496.6 ns [13465.4, 13528.6] (0.95 confidence)</td>
<td>13309.4 ns [13270.5, 13348.3] (0.95 confidence)</td>
<td>1.01407</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d3/strict&#95;result</td>
<td>209.495 ns [209.113, 209.888] (0.95 confidence)</td>
<td>203.738 ns [202.962, 204.545] (0.95 confidence)</td>
<td>1.02825</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d4/rounded&#95;array</td>
<td>210.981 ns [210.418, 211.606] (0.95 confidence)</td>
<td>210.362 ns [210.02, 210.698] (0.95 confidence)</td>
<td>1.00294</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d4/rounded&#95;result</td>
<td>211.704 ns [211.107, 212.313] (0.95 confidence)</td>
<td>211.168 ns [210.784, 211.549] (0.95 confidence)</td>
<td>1.00254</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d4/strict&#95;array</td>
<td>13360.1 ns [13331.3, 13389.5] (0.95 confidence)</td>
<td>13268.7 ns [13247.3, 13290.3] (0.95 confidence)</td>
<td>1.00689</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d4/strict&#95;result</td>
<td>206.516 ns [205.965, 207.108] (0.95 confidence)</td>
<td>212.689 ns [211.768, 213.656] (0.95 confidence)</td>
<td>0.970977</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d5/rounded&#95;array</td>
<td>216.566 ns [216.24, 216.885] (0.95 confidence)</td>
<td>214.722 ns [214.203, 215.25] (0.95 confidence)</td>
<td>1.00858</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d5/rounded&#95;result</td>
<td>216.172 ns [215.93, 216.427] (0.95 confidence)</td>
<td>213.609 ns [213.255, 213.961] (0.95 confidence)</td>
<td>1.012</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d5/strict&#95;array</td>
<td>13295.5 ns [13267.3, 13324.1] (0.95 confidence)</td>
<td>13357.5 ns [13331.4, 13383.9] (0.95 confidence)</td>
<td>0.995355</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d5/strict&#95;result</td>
<td>214.251 ns [213.472, 215.033] (0.95 confidence)</td>
<td>211.412 ns [210.757, 212.078] (0.95 confidence)</td>
<td>1.01343</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d2/rounded&#95;array</td>
<td>165.543 ns [165.033, 166.055] (0.95 confidence)</td>
<td>159.256 ns [158.749, 159.777] (0.95 confidence)</td>
<td>1.03948</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d2/rounded&#95;result</td>
<td>165.643 ns [165.18, 166.14] (0.95 confidence)</td>
<td>156.796 ns [156.168, 157.449] (0.95 confidence)</td>
<td>1.05642</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d2/strict&#95;array</td>
<td>3991.06 ns [3972.22, 4009.62] (0.95 confidence)</td>
<td>3873.98 ns [3862.21, 3885.15] (0.95 confidence)</td>
<td>1.03022</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d2/strict&#95;result</td>
<td>165.746 ns [165.136, 166.35] (0.95 confidence)</td>
<td>160.83 ns [160.249, 161.413] (0.95 confidence)</td>
<td>1.03057</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d3/rounded&#95;array</td>
<td>172.645 ns [171.757, 173.595] (0.95 confidence)</td>
<td>167.883 ns [167.231, 168.523] (0.95 confidence)</td>
<td>1.02836</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d3/rounded&#95;result</td>
<td>172.345 ns [171.769, 172.94] (0.95 confidence)</td>
<td>166.583 ns [166.001, 167.181] (0.95 confidence)</td>
<td>1.03458</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d3/strict&#95;array</td>
<td>3970.98 ns [3953.79, 3987.83] (0.95 confidence)</td>
<td>3871.76 ns [3852.69, 3890.09] (0.95 confidence)</td>
<td>1.02563</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d3/strict&#95;result</td>
<td>172.383 ns [171.781, 172.989] (0.95 confidence)</td>
<td>169.401 ns [168.864, 169.929] (0.95 confidence)</td>
<td>1.0176</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d4/rounded&#95;array</td>
<td>171.215 ns [170.426, 172.012] (0.95 confidence)</td>
<td>171.38 ns [170.782, 171.965] (0.95 confidence)</td>
<td>0.999036</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d4/rounded&#95;result</td>
<td>174.317 ns [173.207, 175.42] (0.95 confidence)</td>
<td>170.132 ns [169.651, 170.634] (0.95 confidence)</td>
<td>1.0246</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d4/strict&#95;array</td>
<td>3964.98 ns [3953.82, 3976.23] (0.95 confidence)</td>
<td>3896.53 ns [3885.1, 3908.12] (0.95 confidence)</td>
<td>1.01757</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d4/strict&#95;result</td>
<td>176.202 ns [174.772, 177.68] (0.95 confidence)</td>
<td>170.176 ns [169.587, 170.768] (0.95 confidence)</td>
<td>1.03541</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d5/rounded&#95;array</td>
<td>180.362 ns [179.712, 180.983] (0.95 confidence)</td>
<td>178.9 ns [178.37, 179.438] (0.95 confidence)</td>
<td>1.00817</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d5/rounded&#95;result</td>
<td>178.475 ns [177.952, 179.005] (0.95 confidence)</td>
<td>176.203 ns [175.548, 176.878] (0.95 confidence)</td>
<td>1.01289</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d5/strict&#95;array</td>
<td>3908.24 ns [3892.08, 3922.95] (0.95 confidence)</td>
<td>3880.23 ns [3860.59, 3899.33] (0.95 confidence)</td>
<td>1.00722</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d5/strict&#95;result</td>
<td>176.547 ns [175.525, 177.605] (0.95 confidence)</td>
<td>174.592 ns [174.057, 175.121] (0.95 confidence)</td>
<td>1.0112</td>
</tr>
<tr>
<td>exact&#95;d2/det</td>
<td>0.449255 ns [0.445354, 0.453483] (0.95 confidence)</td>
<td>0.447541 ns [0.442846, 0.452659] (0.95 confidence)</td>
<td>1.00383</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;direct</td>
<td>0.448425 ns [0.443941, 0.45335] (0.95 confidence)</td>
<td>0.448068 ns [0.443287, 0.453303] (0.95 confidence)</td>
<td>1.0008</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;direct&#95;with&#95;errbound</td>
<td>1.75135 ns [1.7485, 1.75437] (0.95 confidence)</td>
<td>1.72957 ns [1.72702, 1.73214] (0.95 confidence)</td>
<td>1.0126</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;errbound</td>
<td>1.59969 ns [1.59826, 1.60125] (0.95 confidence)</td>
<td>1.58271 ns [1.58211, 1.58333] (0.95 confidence)</td>
<td>1.01073</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;exact</td>
<td>90.2703 ns [90.1181, 90.4338] (0.95 confidence)</td>
<td>89.3689 ns [89.2006, 89.5443] (0.95 confidence)</td>
<td>1.01009</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;exact&#95;f64&#95;result</td>
<td>70.8119 ns [70.6873, 70.9401] (0.95 confidence)</td>
<td>71.3181 ns [71.2077, 71.4244] (0.95 confidence)</td>
<td>0.992902</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;exact&#95;rounded&#95;f64</td>
<td>72.1141 ns [72.0258, 72.2083] (0.95 confidence)</td>
<td>72.9844 ns [72.815, 73.1569] (0.95 confidence)</td>
<td>0.988076</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;sign&#95;exact</td>
<td>3.00237 ns [2.99814, 3.00734] (0.95 confidence)</td>
<td>2.94885 ns [2.94512, 2.9524] (0.95 confidence)</td>
<td>1.01815</td>
</tr>
<tr>
<td>exact&#95;d2/solve&#95;exact</td>
<td>3010.79 ns [3006.86, 3014.94] (0.95 confidence)</td>
<td>2963.29 ns [2959.14, 2967.41] (0.95 confidence)</td>
<td>1.01603</td>
</tr>
<tr>
<td>exact&#95;d2/solve&#95;exact&#95;f64&#95;result</td>
<td>3172.71 ns [3168.79, 3176.79] (0.95 confidence)</td>
<td>3111.76 ns [3107.43, 3115.94] (0.95 confidence)</td>
<td>1.01959</td>
</tr>
<tr>
<td>exact&#95;d2/solve&#95;exact&#95;rounded&#95;f64</td>
<td>3327.59 ns [3320.7, 3334.92] (0.95 confidence)</td>
<td>3258.28 ns [3253.16, 3263.11] (0.95 confidence)</td>
<td>1.02127</td>
</tr>
<tr>
<td>exact&#95;d3/det</td>
<td>0.829398 ns [0.82789, 0.831087] (0.95 confidence)</td>
<td>0.822886 ns [0.820719, 0.825234] (0.95 confidence)</td>
<td>1.00791</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;direct</td>
<td>0.836207 ns [0.834113, 0.838435] (0.95 confidence)</td>
<td>0.818228 ns [0.817938, 0.818533] (0.95 confidence)</td>
<td>1.02197</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;direct&#95;with&#95;errbound</td>
<td>3.31281 ns [3.30911, 3.31687] (0.95 confidence)</td>
<td>3.27819 ns [3.27305, 3.28396] (0.95 confidence)</td>
<td>1.01056</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;errbound</td>
<td>3.26121 ns [3.25719, 3.26573] (0.95 confidence)</td>
<td>3.21987 ns [3.21614, 3.22339] (0.95 confidence)</td>
<td>1.01284</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;exact</td>
<td>309.602 ns [309.022, 310.269] (0.95 confidence)</td>
<td>302.833 ns [302.279, 303.436] (0.95 confidence)</td>
<td>1.02235</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;exact&#95;f64&#95;result</td>
<td>283.344 ns [282.893, 283.814] (0.95 confidence)</td>
<td>276.564 ns [276.272, 276.853] (0.95 confidence)</td>
<td>1.02452</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;exact&#95;rounded&#95;f64</td>
<td>285.45 ns [285.077, 285.847] (0.95 confidence)</td>
<td>280.668 ns [280.182, 281.142] (0.95 confidence)</td>
<td>1.01704</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;sign&#95;exact</td>
<td>4.76021 ns [4.75515, 4.76578] (0.95 confidence)</td>
<td>4.69179 ns [4.68328, 4.70281] (0.95 confidence)</td>
<td>1.01458</td>
</tr>
<tr>
<td>exact&#95;d3/solve&#95;exact</td>
<td>9210.3 ns [9185.47, 9236.6] (0.95 confidence)</td>
<td>9228.34 ns [9208.82, 9247.25] (0.95 confidence)</td>
<td>0.998045</td>
</tr>
<tr>
<td>exact&#95;d3/solve&#95;exact&#95;f64&#95;result</td>
<td>9424.32 ns [9404.06, 9446.25] (0.95 confidence)</td>
<td>9313.39 ns [9277.85, 9345.48] (0.95 confidence)</td>
<td>1.01191</td>
</tr>
<tr>
<td>exact&#95;d3/solve&#95;exact&#95;rounded&#95;f64</td>
<td>9705.51 ns [9682.02, 9728.4] (0.95 confidence)</td>
<td>9908.23 ns [9860.56, 9959.44] (0.95 confidence)</td>
<td>0.97954</td>
</tr>
<tr>
<td>exact&#95;d4/det</td>
<td>4.79131 ns [4.78554, 4.79808] (0.95 confidence)</td>
<td>4.7468 ns [4.74078, 4.75312] (0.95 confidence)</td>
<td>1.00938</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;direct</td>
<td>2.38295 ns [2.37775, 2.3889] (0.95 confidence)</td>
<td>2.34516 ns [2.34386, 2.34632] (0.95 confidence)</td>
<td>1.01611</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;direct&#95;with&#95;errbound</td>
<td>6.83703 ns [6.82991, 6.84515] (0.95 confidence)</td>
<td>6.76662 ns [6.75828, 6.7774] (0.95 confidence)</td>
<td>1.01041</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;errbound</td>
<td>6.80473 ns [6.79507, 6.81478] (0.95 confidence)</td>
<td>6.70808 ns [6.6985, 6.71742] (0.95 confidence)</td>
<td>1.01441</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;exact</td>
<td>859.048 ns [857.475, 860.638] (0.95 confidence)</td>
<td>852.806 ns [851.444, 854.142] (0.95 confidence)</td>
<td>1.00732</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;exact&#95;f64&#95;result</td>
<td>819.954 ns [818.008, 821.921] (0.95 confidence)</td>
<td>814.099 ns [812.759, 815.449] (0.95 confidence)</td>
<td>1.00719</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;exact&#95;rounded&#95;f64</td>
<td>839.83 ns [837.834, 841.842] (0.95 confidence)</td>
<td>829.481 ns [827.95, 830.99] (0.95 confidence)</td>
<td>1.01248</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;sign&#95;exact</td>
<td>7.99013 ns [7.98188, 8.00058] (0.95 confidence)</td>
<td>7.89361 ns [7.88705, 7.89996] (0.95 confidence)</td>
<td>1.01223</td>
</tr>
<tr>
<td>exact&#95;d4/solve&#95;exact</td>
<td>18991 ns [18951.1, 19033.7] (0.95 confidence)</td>
<td>19431.1 ns [19355.2, 19508.4] (0.95 confidence)</td>
<td>0.977353</td>
</tr>
<tr>
<td>exact&#95;d4/solve&#95;exact&#95;f64&#95;result</td>
<td>19240.7 ns [19202.3, 19279.5] (0.95 confidence)</td>
<td>19464.1 ns [19385.6, 19541.3] (0.95 confidence)</td>
<td>0.988518</td>
</tr>
<tr>
<td>exact&#95;d4/solve&#95;exact&#95;rounded&#95;f64</td>
<td>19656.1 ns [19630, 19681.8] (0.95 confidence)</td>
<td>19612.1 ns [19565.6, 19659.2] (0.95 confidence)</td>
<td>1.00224</td>
</tr>
<tr>
<td>exact&#95;d5/det</td>
<td>32.3073 ns [32.2585, 32.3595] (0.95 confidence)</td>
<td>31.8163 ns [31.7767, 31.8538] (0.95 confidence)</td>
<td>1.01543</td>
</tr>
<tr>
<td>exact&#95;d5/det&#95;exact</td>
<td>2927.72 ns [2923.15, 2933.28] (0.95 confidence)</td>
<td>2863.65 ns [2857.44, 2869.91] (0.95 confidence)</td>
<td>1.02237</td>
</tr>
<tr>
<td>exact&#95;d5/det&#95;exact&#95;f64&#95;result</td>
<td>2874.5 ns [2870.74, 2878.43] (0.95 confidence)</td>
<td>2829.98 ns [2823.63, 2836.38] (0.95 confidence)</td>
<td>1.01573</td>
</tr>
<tr>
<td>exact&#95;d5/det&#95;exact&#95;rounded&#95;f64</td>
<td>2908.79 ns [2902.68, 2915.32] (0.95 confidence)</td>
<td>2807.86 ns [2801.38, 2814.61] (0.95 confidence)</td>
<td>1.03595</td>
</tr>
<tr>
<td>exact&#95;d5/det&#95;sign&#95;exact</td>
<td>2877.13 ns [2870.67, 2884.42] (0.95 confidence)</td>
<td>2815.46 ns [2810.01, 2820.78] (0.95 confidence)</td>
<td>1.02191</td>
</tr>
<tr>
<td>exact&#95;d5/solve&#95;exact</td>
<td>31349 ns [31245.7, 31465.3] (0.95 confidence)</td>
<td>31301.3 ns [31062.7, 31552.9] (0.95 confidence)</td>
<td>1.00153</td>
</tr>
<tr>
<td>exact&#95;d5/solve&#95;exact&#95;f64&#95;result</td>
<td>31443.3 ns [31352.9, 31564.3] (0.95 confidence)</td>
<td>31195.5 ns [30994, 31416.7] (0.95 confidence)</td>
<td>1.00794</td>
</tr>
<tr>
<td>exact&#95;d5/solve&#95;exact&#95;rounded&#95;f64</td>
<td>31476.8 ns [31398.5, 31570.7] (0.95 confidence)</td>
<td>31765.7 ns [31562.6, 31985.1] (0.95 confidence)</td>
<td>0.990904</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;4x4/det&#95;exact</td>
<td>931.038 ns [929.607, 932.535] (0.95 confidence)</td>
<td>932.866 ns [930.482, 935.75] (0.95 confidence)</td>
<td>0.99804</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;4x4/det&#95;sign&#95;exact</td>
<td>7.93424 ns [7.92582, 7.94268] (0.95 confidence)</td>
<td>7.89388 ns [7.88643, 7.90083] (0.95 confidence)</td>
<td>1.00511</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;4x4/solve&#95;exact</td>
<td>15742.4 ns [15711.1, 15774.4] (0.95 confidence)</td>
<td>15968.7 ns [15926.5, 16009.8] (0.95 confidence)</td>
<td>0.985827</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;4x4/solve&#95;exact&#95;f64&#95;result</td>
<td>15958.9 ns [15926.1, 15996] (0.95 confidence)</td>
<td>16136.5 ns [16108, 16165.6] (0.95 confidence)</td>
<td>0.988994</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;4x4/solve&#95;exact&#95;rounded&#95;f64</td>
<td>16190.7 ns [16157.8, 16224.6] (0.95 confidence)</td>
<td>16469.6 ns [16426.8, 16512.8] (0.95 confidence)</td>
<td>0.983063</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;5x5/det&#95;exact</td>
<td>2945.96 ns [2942.26, 2949.66] (0.95 confidence)</td>
<td>2899.66 ns [2893.96, 2905.3] (0.95 confidence)</td>
<td>1.01597</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;5x5/det&#95;sign&#95;exact</td>
<td>2898.83 ns [2893.67, 2904.22] (0.95 confidence)</td>
<td>2836.7 ns [2831.67, 2841.47] (0.95 confidence)</td>
<td>1.0219</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;5x5/solve&#95;exact</td>
<td>25238 ns [25185.3, 25295.1] (0.95 confidence)</td>
<td>25948.9 ns [25838.5, 26062.8] (0.95 confidence)</td>
<td>0.972604</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;5x5/solve&#95;exact&#95;f64&#95;result</td>
<td>25534.4 ns [25484, 25590] (0.95 confidence)</td>
<td>26148.6 ns [26036.7, 26264.7] (0.95 confidence)</td>
<td>0.97651</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;5x5/solve&#95;exact&#95;rounded&#95;f64</td>
<td>26010.8 ns [25960.2, 26062.6] (0.95 confidence)</td>
<td>26703.6 ns [26575.4, 26838.6] (0.95 confidence)</td>
<td>0.974056</td>
</tr>
<tr>
<td>exact&#95;large&#95;entries&#95;3x3/det&#95;exact</td>
<td>321.385 ns [320.803, 321.968] (0.95 confidence)</td>
<td>319.457 ns [318.809, 320.105] (0.95 confidence)</td>
<td>1.00603</td>
</tr>
<tr>
<td>exact&#95;large&#95;entries&#95;3x3/det&#95;sign&#95;exact</td>
<td>289.542 ns [288.635, 290.511] (0.95 confidence)</td>
<td>290.418 ns [290.034, 290.796] (0.95 confidence)</td>
<td>0.996984</td>
</tr>
<tr>
<td>exact&#95;large&#95;entries&#95;3x3/solve&#95;exact</td>
<td>9000.92 ns [8986.86, 9015.63] (0.95 confidence)</td>
<td>8829.71 ns [8814.24, 8845.57] (0.95 confidence)</td>
<td>1.01939</td>
</tr>
<tr>
<td>exact&#95;large&#95;entries&#95;3x3/solve&#95;exact&#95;f64&#95;result</td>
<td>9164.87 ns [9148.58, 9181.18] (0.95 confidence)</td>
<td>8959.11 ns [8948.36, 8969.59] (0.95 confidence)</td>
<td>1.02297</td>
</tr>
<tr>
<td>exact&#95;large&#95;entries&#95;3x3/solve&#95;exact&#95;rounded&#95;f64</td>
<td>9455.19 ns [9440.69, 9469.81] (0.95 confidence)</td>
<td>9306.22 ns [9293.93, 9318.11] (0.95 confidence)</td>
<td>1.01601</td>
</tr>
<tr>
<td>exact&#95;near&#95;singular&#95;3x3/det&#95;exact</td>
<td>189.641 ns [189.379, 189.911] (0.95 confidence)</td>
<td>187.746 ns [187.564, 187.918] (0.95 confidence)</td>
<td>1.0101</td>
</tr>
<tr>
<td>exact&#95;near&#95;singular&#95;3x3/det&#95;sign&#95;exact</td>
<td>164.541 ns [164.323, 164.765] (0.95 confidence)</td>
<td>164.571 ns [164.445, 164.695] (0.95 confidence)</td>
<td>0.999818</td>
</tr>
<tr>
<td>exact&#95;near&#95;singular&#95;3x3/solve&#95;exact</td>
<td>996.135 ns [989.345, 1005.7] (0.95 confidence)</td>
<td>959.36 ns [957.876, 960.774] (0.95 confidence)</td>
<td>1.03833</td>
</tr>
<tr>
<td>exact&#95;near&#95;singular&#95;3x3/solve&#95;exact&#95;f64&#95;result</td>
<td>1010.08 ns [1008.43, 1011.8] (0.95 confidence)</td>
<td>986.967 ns [984.749, 989.191] (0.95 confidence)</td>
<td>1.02342</td>
</tr>
<tr>
<td>exact&#95;near&#95;singular&#95;3x3/solve&#95;exact&#95;rounded&#95;f64</td>
<td>993.502 ns [991.931, 995.077] (0.95 confidence)</td>
<td>979.17 ns [977.587, 980.659] (0.95 confidence)</td>
<td>1.01464</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d2/det&#95;exact</td>
<td>3526.85 ns [3519.67, 3534.32] (0.95 confidence)</td>
<td>3451.79 ns [3440.27, 3463.33] (0.95 confidence)</td>
<td>1.02175</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d2/det&#95;sign&#95;exact</td>
<td>147.907 ns [147.485, 148.333] (0.95 confidence)</td>
<td>142.976 ns [142.71, 143.244] (0.95 confidence)</td>
<td>1.03448</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d2/solve&#95;exact</td>
<td>36851.5 ns [36762, 36943.4] (0.95 confidence)</td>
<td>35844 ns [35704.7, 36027.7] (0.95 confidence)</td>
<td>1.02811</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d2/solve&#95;exact&#95;f64&#95;result</td>
<td>37452.9 ns [37336, 37588.6] (0.95 confidence)</td>
<td>36313.4 ns [36240.7, 36388.1] (0.95 confidence)</td>
<td>1.03138</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d2/solve&#95;exact&#95;rounded&#95;f64</td>
<td>37633.3 ns [37529.9, 37745.8] (0.95 confidence)</td>
<td>36710.5 ns [36632.7, 36790.2] (0.95 confidence)</td>
<td>1.02514</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d3/det&#95;exact</td>
<td>9582.23 ns [9559.41, 9607.31] (0.95 confidence)</td>
<td>9281.03 ns [9271.15, 9290.7] (0.95 confidence)</td>
<td>1.03245</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d3/det&#95;sign&#95;exact</td>
<td>233.609 ns [232.96, 234.307] (0.95 confidence)</td>
<td>227.931 ns [227.541, 228.334] (0.95 confidence)</td>
<td>1.02491</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d3/solve&#95;exact</td>
<td>83084.9 ns [82905.8, 83283.1] (0.95 confidence)</td>
<td>79851.3 ns [79751.4, 79947.2] (0.95 confidence)</td>
<td>1.04049</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d3/solve&#95;exact&#95;f64&#95;result</td>
<td>83507.6 ns [83337.4, 83694.1] (0.95 confidence)</td>
<td>80508.4 ns [80309.8, 80763.8] (0.95 confidence)</td>
<td>1.03725</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d3/solve&#95;exact&#95;rounded&#95;f64</td>
<td>84399.2 ns [84188.3, 84634.5] (0.95 confidence)</td>
<td>80824.7 ns [80705.5, 80941.8] (0.95 confidence)</td>
<td>1.04423</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d4/det&#95;exact</td>
<td>23241.1 ns [23160.9, 23324.2] (0.95 confidence)</td>
<td>22304.6 ns [22258.4, 22355.8] (0.95 confidence)</td>
<td>1.04199</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d4/det&#95;sign&#95;exact</td>
<td>447.521 ns [446.72, 448.41] (0.95 confidence)</td>
<td>438.942 ns [438.073, 439.764] (0.95 confidence)</td>
<td>1.01954</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d4/solve&#95;exact</td>
<td>150062 ns [149672, 150479] (0.95 confidence)</td>
<td>143347 ns [143175, 143513] (0.95 confidence)</td>
<td>1.04685</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d4/solve&#95;exact&#95;f64&#95;result</td>
<td>150338 ns [149966, 150748] (0.95 confidence)</td>
<td>144318 ns [144128, 144510] (0.95 confidence)</td>
<td>1.04171</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d4/solve&#95;exact&#95;rounded&#95;f64</td>
<td>154891 ns [154304, 155475] (0.95 confidence)</td>
<td>145634 ns [145060, 146395] (0.95 confidence)</td>
<td>1.06356</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d5/det&#95;exact</td>
<td>57476.8 ns [57202.7, 57775.1] (0.95 confidence)</td>
<td>55673.3 ns [55493.1, 55864] (0.95 confidence)</td>
<td>1.03239</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d5/det&#95;sign&#95;exact</td>
<td>56380.2 ns [56184.1, 56590.5] (0.95 confidence)</td>
<td>54406.2 ns [54276.3, 54537.2] (0.95 confidence)</td>
<td>1.03628</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d5/solve&#95;exact</td>
<td>249318 ns [248574, 250121] (0.95 confidence)</td>
<td>236798 ns [236477, 237124] (0.95 confidence)</td>
<td>1.05287</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d5/solve&#95;exact&#95;f64&#95;result</td>
<td>249192 ns [248488, 249942] (0.95 confidence)</td>
<td>237145 ns [236846, 237436] (0.95 confidence)</td>
<td>1.0508</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d5/solve&#95;exact&#95;rounded&#95;f64</td>
<td>249842 ns [249217, 250500] (0.95 confidence)</td>
<td>238420 ns [238164, 238681] (0.95 confidence)</td>
<td>1.04791</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/matrix&#95;access</td>
<td>0.751529 ns [0.749058, 0.75377] (0.95 confidence)</td>
<td>0.747059 ns [0.744651, 0.749268] (0.95 confidence)</td>
<td>1.00598</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/matrix&#95;from&#95;fn</td>
<td>131.062 ns [130.848, 131.26] (0.95 confidence)</td>
<td>130.202 ns [129.977, 130.408] (0.95 confidence)</td>
<td>1.00661</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/matrix&#95;from&#95;rows</td>
<td>133.558 ns [133.323, 133.779] (0.95 confidence)</td>
<td>133.09 ns [132.877, 133.297] (0.95 confidence)</td>
<td>1.00352</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/matrix&#95;into&#95;rows</td>
<td>63.9908 ns [63.9037, 64.0756] (0.95 confidence)</td>
<td>63.5106 ns [63.3737, 63.6474] (0.95 confidence)</td>
<td>1.00756</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/matrix&#95;zero</td>
<td>2.68748 ns [2.68251, 2.69315] (0.95 confidence)</td>
<td>2.67832 ns [2.67536, 2.68103] (0.95 confidence)</td>
<td>1.00342</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/vector&#95;access</td>
<td>0.756667 ns [0.752529, 0.761265] (0.95 confidence)</td>
<td>0.748351 ns [0.746071, 0.750484] (0.95 confidence)</td>
<td>1.01111</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/vector&#95;from&#95;fn</td>
<td>36.5738 ns [36.5215, 36.623] (0.95 confidence)</td>
<td>36.3829 ns [36.3162, 36.4462] (0.95 confidence)</td>
<td>1.00525</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/vector&#95;into&#95;array</td>
<td>15.2952 ns [15.2527, 15.3342] (0.95 confidence)</td>
<td>15.2534 ns [15.2221, 15.2832] (0.95 confidence)</td>
<td>1.00274</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/vector&#95;new</td>
<td>35.87 ns [35.8099, 35.9277] (0.95 confidence)</td>
<td>35.774 ns [35.7171, 35.8291] (0.95 confidence)</td>
<td>1.00268</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/vector&#95;zero</td>
<td>1.35184 ns [1.35073, 1.35291] (0.95 confidence)</td>
<td>1.34655 ns [1.34504, 1.34795] (0.95 confidence)</td>
<td>1.00393</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/matrix&#95;access</td>
<td>0.754236 ns [0.751877, 0.756363] (0.95 confidence)</td>
<td>0.748654 ns [0.745531, 0.752012] (0.95 confidence)</td>
<td>1.00746</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/matrix&#95;from&#95;fn</td>
<td>227.19 ns [226.829, 227.528] (0.95 confidence)</td>
<td>226.518 ns [226.131, 226.88] (0.95 confidence)</td>
<td>1.00297</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/matrix&#95;from&#95;rows</td>
<td>264.675 ns [264.267, 265.093] (0.95 confidence)</td>
<td>262.656 ns [262.092, 263.16] (0.95 confidence)</td>
<td>1.00769</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/matrix&#95;into&#95;rows</td>
<td>122.417 ns [122.202, 122.623] (0.95 confidence)</td>
<td>121.691 ns [121.449, 121.919] (0.95 confidence)</td>
<td>1.00596</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/matrix&#95;zero</td>
<td>6.20649 ns [6.1894, 6.22252] (0.95 confidence)</td>
<td>6.15085 ns [6.13273, 6.17613] (0.95 confidence)</td>
<td>1.00905</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/vector&#95;access</td>
<td>0.754942 ns [0.752503, 0.757226] (0.95 confidence)</td>
<td>0.749935 ns [0.746409, 0.753198] (0.95 confidence)</td>
<td>1.00668</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/vector&#95;from&#95;fn</td>
<td>51.6957 ns [51.6308, 51.7551] (0.95 confidence)</td>
<td>51.4255 ns [51.3512, 51.4943] (0.95 confidence)</td>
<td>1.00525</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/vector&#95;into&#95;array</td>
<td>22.2144 ns [22.1744, 22.2525] (0.95 confidence)</td>
<td>22.1676 ns [22.1263, 22.2072] (0.95 confidence)</td>
<td>1.00211</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/vector&#95;new</td>
<td>51.3898 ns [51.2938, 51.5167] (0.95 confidence)</td>
<td>51.2305 ns [51.1185, 51.3574] (0.95 confidence)</td>
<td>1.00311</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/vector&#95;zero</td>
<td>2.01399 ns [2.01082, 2.01701] (0.95 confidence)</td>
<td>2.00192 ns [1.99905, 2.00465] (0.95 confidence)</td>
<td>1.00603</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/matrix&#95;access</td>
<td>0.751811 ns [0.749268, 0.754166] (0.95 confidence)</td>
<td>0.756091 ns [0.753499, 0.758559] (0.95 confidence)</td>
<td>0.994339</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/matrix&#95;from&#95;fn</td>
<td>283.609 ns [283.161, 284.027] (0.95 confidence)</td>
<td>285.79 ns [285.124, 286.473] (0.95 confidence)</td>
<td>0.99237</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/matrix&#95;from&#95;rows</td>
<td>341.977 ns [341.401, 342.528] (0.95 confidence)</td>
<td>358.904 ns [356.891, 361.239] (0.95 confidence)</td>
<td>0.952837</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/matrix&#95;into&#95;rows</td>
<td>129.546 ns [129.304, 129.778] (0.95 confidence)</td>
<td>129.943 ns [129.699, 130.163] (0.95 confidence)</td>
<td>0.996947</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/matrix&#95;zero</td>
<td>15.2142 ns [15.1344, 15.2964] (0.95 confidence)</td>
<td>15.614 ns [15.5266, 15.702] (0.95 confidence)</td>
<td>0.974392</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/vector&#95;access</td>
<td>0.754794 ns [0.752279, 0.757307] (0.95 confidence)</td>
<td>0.753732 ns [0.751386, 0.755874] (0.95 confidence)</td>
<td>1.00141</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/vector&#95;from&#95;fn</td>
<td>65.7253 ns [65.6388, 65.8008] (0.95 confidence)</td>
<td>65.6807 ns [65.531, 65.7965] (0.95 confidence)</td>
<td>1.00068</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/vector&#95;into&#95;array</td>
<td>29.5384 ns [29.4935, 29.5809] (0.95 confidence)</td>
<td>29.5248 ns [29.4765, 29.569] (0.95 confidence)</td>
<td>1.00046</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/vector&#95;new</td>
<td>64.8987 ns [64.7799, 65.0304] (0.95 confidence)</td>
<td>65.0267 ns [64.8634, 65.2353] (0.95 confidence)</td>
<td>0.998032</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/vector&#95;zero</td>
<td>3.28527 ns [3.26721, 3.30459] (0.95 confidence)</td>
<td>3.38094 ns [3.36049, 3.40135] (0.95 confidence)</td>
<td>0.971703</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/matrix&#95;access</td>
<td>0.756624 ns [0.752221, 0.761855] (0.95 confidence)</td>
<td>0.75288 ns [0.750334, 0.755095] (0.95 confidence)</td>
<td>1.00497</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/matrix&#95;from&#95;fn</td>
<td>429.038 ns [427.452, 430.479] (0.95 confidence)</td>
<td>430.446 ns [429.141, 431.751] (0.95 confidence)</td>
<td>0.996728</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/matrix&#95;from&#95;rows</td>
<td>620.102 ns [619.185, 621.007] (0.95 confidence)</td>
<td>617.669 ns [616.708, 618.615] (0.95 confidence)</td>
<td>1.00394</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/matrix&#95;into&#95;rows</td>
<td>301.924 ns [301.055, 302.875] (0.95 confidence)</td>
<td>299.72 ns [299.009, 300.41] (0.95 confidence)</td>
<td>1.00735</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/matrix&#95;zero</td>
<td>21.4612 ns [21.4228, 21.4971] (0.95 confidence)</td>
<td>21.5768 ns [21.5224, 21.633] (0.95 confidence)</td>
<td>0.994639</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/vector&#95;access</td>
<td>0.749079 ns [0.744471, 0.753981] (0.95 confidence)</td>
<td>0.753508 ns [0.750182, 0.757303] (0.95 confidence)</td>
<td>0.994122</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/vector&#95;from&#95;fn</td>
<td>81.8699 ns [81.7542, 81.9831] (0.95 confidence)</td>
<td>81.4132 ns [81.2832, 81.5377] (0.95 confidence)</td>
<td>1.00561</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/vector&#95;into&#95;array</td>
<td>37.9196 ns [37.8491, 37.9877] (0.95 confidence)</td>
<td>37.8679 ns [37.7922, 37.9425] (0.95 confidence)</td>
<td>1.00137</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/vector&#95;new</td>
<td>81.9742 ns [81.8559, 82.0904] (0.95 confidence)</td>
<td>81.5544 ns [81.4449, 81.66] (0.95 confidence)</td>
<td>1.00515</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/vector&#95;zero</td>
<td>4.23138 ns [4.21499, 4.24758] (0.95 confidence)</td>
<td>4.25998 ns [4.23258, 4.28889] (0.95 confidence)</td>
<td>0.993288</td>
</tr>
<tr>
<td>public&#95;exact&#95;scalar/rounded</td>
<td>5.54531 ns [5.53979, 5.55034] (0.95 confidence)</td>
<td>5.54759 ns [5.53952, 5.55639] (0.95 confidence)</td>
<td>0.99959</td>
</tr>
<tr>
<td>public&#95;exact&#95;scalar/strict</td>
<td>4.21264 ns [4.19771, 4.22887] (0.95 confidence)</td>
<td>4.26761 ns [4.24519, 4.29077] (0.95 confidence)</td>
<td>0.98712</td>
</tr>
<tr>
<td>rational&#95;input&#95;d2/det&#95;row&#95;cleared&#95;bareiss</td>
<td>179.184 ns [178.732, 179.58] (0.95 confidence)</td>
<td>177.82 ns [177.567, 178.078] (0.95 confidence)</td>
<td>1.00767</td>
</tr>
<tr>
<td>rational&#95;input&#95;d2/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>76.6023 ns [76.4682, 76.754] (0.95 confidence)</td>
<td>76.437 ns [76.3327, 76.5281] (0.95 confidence)</td>
<td>1.00216</td>
</tr>
<tr>
<td>rational&#95;input&#95;d2/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>624.687 ns [623.292, 626.275] (0.95 confidence)</td>
<td>617.227 ns [616.363, 618.073] (0.95 confidence)</td>
<td>1.01209</td>
</tr>
<tr>
<td>rational&#95;input&#95;d3/det&#95;row&#95;cleared&#95;bareiss</td>
<td>593.489 ns [592.472, 594.458] (0.95 confidence)</td>
<td>586.124 ns [584.807, 587.434] (0.95 confidence)</td>
<td>1.01257</td>
</tr>
<tr>
<td>rational&#95;input&#95;d3/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>369.366 ns [368.312, 370.456] (0.95 confidence)</td>
<td>369.984 ns [369.01, 371.078] (0.95 confidence)</td>
<td>0.998329</td>
</tr>
<tr>
<td>rational&#95;input&#95;d3/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>1263.39 ns [1261.08, 1265.74] (0.95 confidence)</td>
<td>1239.88 ns [1237.95, 1241.72] (0.95 confidence)</td>
<td>1.01896</td>
</tr>
<tr>
<td>rational&#95;input&#95;d4/det&#95;row&#95;cleared&#95;bareiss</td>
<td>1223.66 ns [1222.41, 1224.85] (0.95 confidence)</td>
<td>1207.56 ns [1206.03, 1209.03] (0.95 confidence)</td>
<td>1.01333</td>
</tr>
<tr>
<td>rational&#95;input&#95;d4/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>819.806 ns [816.904, 822.447] (0.95 confidence)</td>
<td>804.985 ns [801.786, 808.08] (0.95 confidence)</td>
<td>1.01841</td>
</tr>
<tr>
<td>rational&#95;input&#95;d4/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>2239.48 ns [2237.54, 2241.42] (0.95 confidence)</td>
<td>2216.19 ns [2213.63, 2218.63] (0.95 confidence)</td>
<td>1.01051</td>
</tr>
<tr>
<td>rational&#95;input&#95;d5/det&#95;row&#95;cleared&#95;bareiss</td>
<td>2111.89 ns [2108.55, 2115.31] (0.95 confidence)</td>
<td>2051.33 ns [2047.88, 2054.85] (0.95 confidence)</td>
<td>1.02952</td>
</tr>
<tr>
<td>rational&#95;input&#95;d5/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>1562.05 ns [1559.79, 1564.38] (0.95 confidence)</td>
<td>1551.24 ns [1547.67, 1555.07] (0.95 confidence)</td>
<td>1.00697</td>
</tr>
<tr>
<td>rational&#95;input&#95;d5/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>3373.95 ns [3370.11, 3377.83] (0.95 confidence)</td>
<td>3340.08 ns [3335.91, 3344.21] (0.95 confidence)</td>
<td>1.01014</td>
</tr>
<tr>
<td>rational&#95;input&#95;d6/det&#95;row&#95;cleared&#95;bareiss</td>
<td>3534.49 ns [3529.21, 3539.86] (0.95 confidence)</td>
<td>3441.8 ns [3434.81, 3448.63] (0.95 confidence)</td>
<td>1.02693</td>
</tr>
<tr>
<td>rational&#95;input&#95;d6/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>2830.57 ns [2827.05, 2834.17] (0.95 confidence)</td>
<td>2810.14 ns [2805.15, 2815.36] (0.95 confidence)</td>
<td>1.00727</td>
</tr>
<tr>
<td>rational&#95;input&#95;d6/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>5701.4 ns [5690.33, 5713.15] (0.95 confidence)</td>
<td>5634.52 ns [5626.52, 5642.76] (0.95 confidence)</td>
<td>1.01187</td>
</tr>
<tr>
<td>rational&#95;input&#95;d7/det&#95;row&#95;cleared&#95;bareiss</td>
<td>5561.44 ns [5551.53, 5572.24] (0.95 confidence)</td>
<td>5518.8 ns [5508.72, 5527.79] (0.95 confidence)</td>
<td>1.00773</td>
</tr>
<tr>
<td>rational&#95;input&#95;d7/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>4696.04 ns [4670.69, 4724.56] (0.95 confidence)</td>
<td>4595.3 ns [4589.71, 4600.94] (0.95 confidence)</td>
<td>1.02192</td>
</tr>
<tr>
<td>rational&#95;input&#95;d7/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>9856.3 ns [9844.68, 9867.41] (0.95 confidence)</td>
<td>9823.52 ns [9808.37, 9838.74] (0.95 confidence)</td>
<td>1.00334</td>
</tr>
<tr>
<td>rational&#95;input&#95;d8/det&#95;row&#95;cleared&#95;bareiss</td>
<td>8090.91 ns [8083.13, 8098.62] (0.95 confidence)</td>
<td>8110.89 ns [8095.03, 8126.33] (0.95 confidence)</td>
<td>0.997536</td>
</tr>
<tr>
<td>rational&#95;input&#95;d8/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>7260.69 ns [7250.92, 7270.38] (0.95 confidence)</td>
<td>7231.61 ns [7222.95, 7240.22] (0.95 confidence)</td>
<td>1.00402</td>
</tr>
<tr>
<td>rational&#95;input&#95;d8/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>15113.4 ns [15076.1, 15150.6] (0.95 confidence)</td>
<td>14971.4 ns [14948.6, 14993.9] (0.95 confidence)</td>
<td>1.00948</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d2/det&#95;row&#95;cleared&#95;bareiss</td>
<td>42145 ns [42094.7, 42197.2] (0.95 confidence)</td>
<td>41928.8 ns [41876.2, 41982.3] (0.95 confidence)</td>
<td>1.00516</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d2/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>942.862 ns [942.022, 943.674] (0.95 confidence)</td>
<td>919.683 ns [918.338, 921.003] (0.95 confidence)</td>
<td>1.0252</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d2/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>3038.11 ns [3033.98, 3042.28] (0.95 confidence)</td>
<td>2999.94 ns [2994.8, 3004.44] (0.95 confidence)</td>
<td>1.01272</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d3/det&#95;row&#95;cleared&#95;bareiss</td>
<td>69061.8 ns [69002.4, 69124.9] (0.95 confidence)</td>
<td>68848.2 ns [68770.9, 68923.5] (0.95 confidence)</td>
<td>1.0031</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d3/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>4195 ns [4185.86, 4205.31] (0.95 confidence)</td>
<td>4150.08 ns [4142.76, 4157.46] (0.95 confidence)</td>
<td>1.01082</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d3/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>11360.2 ns [11346.4, 11373.8] (0.95 confidence)</td>
<td>11124.7 ns [11112.2, 11136.6] (0.95 confidence)</td>
<td>1.02117</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d4/det&#95;row&#95;cleared&#95;bareiss</td>
<td>110136 ns [109997, 110277] (0.95 confidence)</td>
<td>109835 ns [109726, 109937] (0.95 confidence)</td>
<td>1.00274</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d4/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>9428.13 ns [9414.96, 9441.48] (0.95 confidence)</td>
<td>9297.44 ns [9283.84, 9311.49] (0.95 confidence)</td>
<td>1.01406</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d4/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>33422.7 ns [33370, 33474.3] (0.95 confidence)</td>
<td>32986.8 ns [32806.4, 33229.7] (0.95 confidence)</td>
<td>1.01322</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d5/det&#95;row&#95;cleared&#95;bareiss</td>
<td>185962 ns [185819, 186114] (0.95 confidence)</td>
<td>184465 ns [184255, 184674] (0.95 confidence)</td>
<td>1.00811</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d5/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>42825 ns [42692.4, 42937.7] (0.95 confidence)</td>
<td>41772.7 ns [41700.2, 41843.6] (0.95 confidence)</td>
<td>1.02519</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d5/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>69038.7 ns [68940.8, 69131.8] (0.95 confidence)</td>
<td>67469.6 ns [67371.9, 67563.4] (0.95 confidence)</td>
<td>1.02326</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d6/det&#95;row&#95;cleared&#95;bareiss</td>
<td>301791 ns [301343, 302230] (0.95 confidence)</td>
<td>296791 ns [296366, 297197] (0.95 confidence)</td>
<td>1.01685</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d6/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>104728 ns [104523, 104929] (0.95 confidence)</td>
<td>101040 ns [100864, 101210] (0.95 confidence)</td>
<td>1.0365</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d6/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>160514 ns [160249, 160782] (0.95 confidence)</td>
<td>155980 ns [155584, 156448] (0.95 confidence)</td>
<td>1.02907</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d7/det&#95;row&#95;cleared&#95;bareiss</td>
<td>500896 ns [500165, 501653] (0.95 confidence)</td>
<td>493423 ns [492862, 493978] (0.95 confidence)</td>
<td>1.01515</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d7/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>232631 ns [232230, 233032] (0.95 confidence)</td>
<td>223979 ns [223557, 224377] (0.95 confidence)</td>
<td>1.03863</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d7/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>337301 ns [336707, 337887] (0.95 confidence)</td>
<td>327617 ns [327174, 328051] (0.95 confidence)</td>
<td>1.02956</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d8/det&#95;row&#95;cleared&#95;bareiss</td>
<td>821638 ns [820795, 822484] (0.95 confidence)</td>
<td>806098 ns [805335, 806841] (0.95 confidence)</td>
<td>1.01928</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d8/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>466313 ns [464701, 468035] (0.95 confidence)</td>
<td>452486 ns [451811, 453168] (0.95 confidence)</td>
<td>1.03056</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d8/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>646144 ns [645047, 647208] (0.95 confidence)</td>
<td>630464 ns [629514, 631398] (0.95 confidence)</td>
<td>1.02487</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d2/det&#95;row&#95;cleared&#95;bareiss</td>
<td>8782.11 ns [8763.14, 8801.14] (0.95 confidence)</td>
<td>8809.58 ns [8788.84, 8830.95] (0.95 confidence)</td>
<td>0.996882</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d2/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>645.504 ns [643.531, 647.57] (0.95 confidence)</td>
<td>648.316 ns [646.576, 650.143] (0.95 confidence)</td>
<td>0.995664</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d2/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>2071.44 ns [2063.84, 2079.02] (0.95 confidence)</td>
<td>2025.63 ns [2020.94, 2030.34] (0.95 confidence)</td>
<td>1.02261</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d3/det&#95;row&#95;cleared&#95;bareiss</td>
<td>13206.3 ns [13189.8, 13223.6] (0.95 confidence)</td>
<td>13305.2 ns [13281.3, 13329.8] (0.95 confidence)</td>
<td>0.992561</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d3/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>2759.13 ns [2752.87, 2765.54] (0.95 confidence)</td>
<td>2727.17 ns [2718.51, 2735.84] (0.95 confidence)</td>
<td>1.01172</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d3/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>5791.01 ns [5779.75, 5802.25] (0.95 confidence)</td>
<td>5741.98 ns [5730.41, 5753.87] (0.95 confidence)</td>
<td>1.00854</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d4/det&#95;row&#95;cleared&#95;bareiss</td>
<td>21987.8 ns [21942.6, 22033] (0.95 confidence)</td>
<td>21756.5 ns [21722.1, 21791.4] (0.95 confidence)</td>
<td>1.01063</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d4/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>5965.15 ns [5949.15, 5981.56] (0.95 confidence)</td>
<td>5818.1 ns [5805.53, 5831.13] (0.95 confidence)</td>
<td>1.02528</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d4/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>12862 ns [12834.3, 12891.3] (0.95 confidence)</td>
<td>12523.6 ns [12494.9, 12551.9] (0.95 confidence)</td>
<td>1.02702</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d5/det&#95;row&#95;cleared&#95;bareiss</td>
<td>33538.4 ns [33499.8, 33578] (0.95 confidence)</td>
<td>33542 ns [33466.1, 33628.2] (0.95 confidence)</td>
<td>0.999891</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d5/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>13405.2 ns [13381, 13429.9] (0.95 confidence)</td>
<td>13089.2 ns [13056.8, 13122] (0.95 confidence)</td>
<td>1.02414</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d5/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>20806.8 ns [20773.4, 20839.2] (0.95 confidence)</td>
<td>20315.3 ns [20267.1, 20364] (0.95 confidence)</td>
<td>1.02419</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d6/det&#95;row&#95;cleared&#95;bareiss</td>
<td>52080.1 ns [52009.1, 52152.8] (0.95 confidence)</td>
<td>51786.8 ns [51686, 51887.4] (0.95 confidence)</td>
<td>1.00566</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d6/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>26183.2 ns [26131.8, 26245.3] (0.95 confidence)</td>
<td>25619.7 ns [25573.2, 25665.6] (0.95 confidence)</td>
<td>1.02199</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d6/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>39503.9 ns [39438.1, 39568.3] (0.95 confidence)</td>
<td>38387.8 ns [38309.9, 38461.9] (0.95 confidence)</td>
<td>1.02907</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d7/det&#95;row&#95;cleared&#95;bareiss</td>
<td>83172.8 ns [83053.2, 83292.2] (0.95 confidence)</td>
<td>81984.1 ns [81869.8, 82101.1] (0.95 confidence)</td>
<td>1.0145</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d7/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>50174.5 ns [50084.3, 50264.9] (0.95 confidence)</td>
<td>48713.2 ns [48622.4, 48804.3] (0.95 confidence)</td>
<td>1.03</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d7/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>72284.2 ns [72181.2, 72384.9] (0.95 confidence)</td>
<td>70085.5 ns [69941.5, 70213.6] (0.95 confidence)</td>
<td>1.03137</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d8/det&#95;row&#95;cleared&#95;bareiss</td>
<td>131248 ns [131082, 131415] (0.95 confidence)</td>
<td>129567 ns [129362, 129771] (0.95 confidence)</td>
<td>1.01297</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d8/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>90260.6 ns [90108.6, 90412.5] (0.95 confidence)</td>
<td>89103.8 ns [88981.1, 89220.4] (0.95 confidence)</td>
<td>1.01298</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d8/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>127188 ns [126868, 127513] (0.95 confidence)</td>
<td>123161 ns [122970, 123344] (0.95 confidence)</td>
<td>1.0327</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d2/adapter</td>
<td>1069.15 ns [1053.43, 1085.78] (0.95 confidence)</td>
<td>1050.86 ns [1035.93, 1066.73] (0.95 confidence)</td>
<td>1.0174</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d2/construct</td>
<td>292.153 ns [291.684, 292.621] (0.95 confidence)</td>
<td>287.242 ns [286.781, 287.708] (0.95 confidence)</td>
<td>1.0171</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d2/prepared</td>
<td>651.986 ns [650.248, 653.742] (0.95 confidence)</td>
<td>637.315 ns [635.927, 638.711] (0.95 confidence)</td>
<td>1.02302</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d3/adapter</td>
<td>1964.28 ns [1960.13, 1968.82] (0.95 confidence)</td>
<td>1930.37 ns [1926.66, 1934.03] (0.95 confidence)</td>
<td>1.01756</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d3/construct</td>
<td>651.102 ns [649.393, 652.878] (0.95 confidence)</td>
<td>635.788 ns [634.499, 637.176] (0.95 confidence)</td>
<td>1.02409</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d3/prepared</td>
<td>1172.22 ns [1170.77, 1173.66] (0.95 confidence)</td>
<td>1143.14 ns [1141.44, 1144.72] (0.95 confidence)</td>
<td>1.02543</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d4/adapter</td>
<td>3295.62 ns [3290.87, 3300.13] (0.95 confidence)</td>
<td>3322.59 ns [3313.87, 3332.75] (0.95 confidence)</td>
<td>0.991884</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d4/construct</td>
<td>1144.45 ns [1108.54, 1183.57] (0.95 confidence)</td>
<td>1138.17 ns [1103.72, 1175.58] (0.95 confidence)</td>
<td>1.00552</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d4/prepared</td>
<td>2066.83 ns [2061.12, 2072.73] (0.95 confidence)</td>
<td>2114.1 ns [2109.89, 2118.34] (0.95 confidence)</td>
<td>0.977638</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d5/adapter</td>
<td>5398.01 ns [5387.47, 5408.81] (0.95 confidence)</td>
<td>5369 ns [5358.23, 5379.99] (0.95 confidence)</td>
<td>1.0054</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d5/construct</td>
<td>1931.37 ns [1926.04, 1937.22] (0.95 confidence)</td>
<td>1927.86 ns [1922.75, 1933.29] (0.95 confidence)</td>
<td>1.00182</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d5/prepared</td>
<td>3157.72 ns [3152.06, 3163.31] (0.95 confidence)</td>
<td>3152.58 ns [3144.64, 3160.55] (0.95 confidence)</td>
<td>1.00163</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d6/adapter</td>
<td>7473.47 ns [7463.27, 7483.72] (0.95 confidence)</td>
<td>7672.14 ns [7664.83, 7679.67] (0.95 confidence)</td>
<td>0.974105</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d6/construct</td>
<td>2514.69 ns [2506.96, 2523] (0.95 confidence)</td>
<td>2529.39 ns [2524.06, 2535.32] (0.95 confidence)</td>
<td>0.994188</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d6/prepared</td>
<td>4626.1 ns [4618.75, 4633.16] (0.95 confidence)</td>
<td>4725.34 ns [4717.18, 4734.09] (0.95 confidence)</td>
<td>0.978998</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d7/adapter</td>
<td>10252.6 ns [10238.7, 10266.1] (0.95 confidence)</td>
<td>10309.5 ns [10287.7, 10334] (0.95 confidence)</td>
<td>0.994478</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d7/construct</td>
<td>3700.43 ns [3595.65, 3814.55] (0.95 confidence)</td>
<td>3693.38 ns [3582.98, 3813.07] (0.95 confidence)</td>
<td>1.00191</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d7/prepared</td>
<td>6401.06 ns [6387.77, 6414.28] (0.95 confidence)</td>
<td>6567.72 ns [6550.18, 6584.74] (0.95 confidence)</td>
<td>0.974623</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d8/adapter</td>
<td>13258.4 ns [13237.6, 13279.4] (0.95 confidence)</td>
<td>13419.7 ns [13397.2, 13440.3] (0.95 confidence)</td>
<td>0.98798</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d8/construct</td>
<td>4043.41 ns [4034.94, 4051.97] (0.95 confidence)</td>
<td>4053.69 ns [4046.49, 4060.62] (0.95 confidence)</td>
<td>0.997463</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d8/prepared</td>
<td>8962.36 ns [8941.48, 8984.86] (0.95 confidence)</td>
<td>9220.66 ns [9208, 9233.22] (0.95 confidence)</td>
<td>0.971987</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d2/adapter</td>
<td>1069.22 ns [1060.8, 1078.2] (0.95 confidence)</td>
<td>1045.53 ns [1037.39, 1054.13] (0.95 confidence)</td>
<td>1.02266</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d2/construct</td>
<td>328.471 ns [327.343, 329.573] (0.95 confidence)</td>
<td>322.339 ns [321.295, 323.298] (0.95 confidence)</td>
<td>1.01902</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d2/prepared</td>
<td>634.424 ns [633.601, 635.274] (0.95 confidence)</td>
<td>621.335 ns [620.04, 622.625] (0.95 confidence)</td>
<td>1.02107</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d3/adapter</td>
<td>2116.06 ns [2113.3, 2118.67] (0.95 confidence)</td>
<td>2072.62 ns [2070.26, 2074.81] (0.95 confidence)</td>
<td>1.02096</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d3/construct</td>
<td>683.075 ns [680.706, 685.394] (0.95 confidence)</td>
<td>661.778 ns [659.297, 664.168] (0.95 confidence)</td>
<td>1.03218</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d3/prepared</td>
<td>1299.86 ns [1298.48, 1301.27] (0.95 confidence)</td>
<td>1275.04 ns [1273.54, 1276.48] (0.95 confidence)</td>
<td>1.01947</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d4/adapter</td>
<td>3506.54 ns [3502.78, 3510.36] (0.95 confidence)</td>
<td>3473.56 ns [3468.69, 3478.15] (0.95 confidence)</td>
<td>1.0095</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d4/construct</td>
<td>1225.23 ns [1187.15, 1267.05] (0.95 confidence)</td>
<td>1200.45 ns [1162.13, 1242.01] (0.95 confidence)</td>
<td>1.02065</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d4/prepared</td>
<td>2192.55 ns [2187.53, 2198.77] (0.95 confidence)</td>
<td>2168.69 ns [2165.39, 2171.83] (0.95 confidence)</td>
<td>1.011</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d5/adapter</td>
<td>5567.37 ns [5561.52, 5572.99] (0.95 confidence)</td>
<td>5503.67 ns [5495.06, 5512.13] (0.95 confidence)</td>
<td>1.01157</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d5/construct</td>
<td>1886.81 ns [1879.48, 1893.7] (0.95 confidence)</td>
<td>1857.77 ns [1848.09, 1867.03] (0.95 confidence)</td>
<td>1.01564</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d5/prepared</td>
<td>3410.79 ns [3406.76, 3414.61] (0.95 confidence)</td>
<td>3370.27 ns [3365.83, 3374.6] (0.95 confidence)</td>
<td>1.01202</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d6/adapter</td>
<td>7896.6 ns [7887.56, 7905.6] (0.95 confidence)</td>
<td>7876.89 ns [7862.1, 7891.41] (0.95 confidence)</td>
<td>1.0025</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d6/construct</td>
<td>2718.71 ns [2706.12, 2731.25] (0.95 confidence)</td>
<td>2711.58 ns [2696.7, 2725.78] (0.95 confidence)</td>
<td>1.00263</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d6/prepared</td>
<td>4863.3 ns [4857.86, 4868.72] (0.95 confidence)</td>
<td>4861.53 ns [4849.51, 4878.33] (0.95 confidence)</td>
<td>1.00037</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d7/adapter</td>
<td>11566.4 ns [11552, 11580.1] (0.95 confidence)</td>
<td>11423.3 ns [11403.6, 11443.7] (0.95 confidence)</td>
<td>1.01252</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d7/construct</td>
<td>4095.49 ns [3976.7, 4224.04] (0.95 confidence)</td>
<td>4023.02 ns [3908.7, 4148.51] (0.95 confidence)</td>
<td>1.01801</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d7/prepared</td>
<td>7355.88 ns [7333.71, 7378.83] (0.95 confidence)</td>
<td>7319.14 ns [7309.83, 7327.94] (0.95 confidence)</td>
<td>1.00502</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d8/adapter</td>
<td>15786.4 ns [15756, 15816.5] (0.95 confidence)</td>
<td>15794.2 ns [15758.7, 15828] (0.95 confidence)</td>
<td>0.999509</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d8/construct</td>
<td>5104.65 ns [5079.91, 5128.54] (0.95 confidence)</td>
<td>5072.82 ns [5041.77, 5103.39] (0.95 confidence)</td>
<td>1.00627</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d8/prepared</td>
<td>10297.7 ns [10276.9, 10316.9] (0.95 confidence)</td>
<td>10389.2 ns [10370.2, 10408.6] (0.95 confidence)</td>
<td>0.991187</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d2/adapter</td>
<td>448.2 ns [440.824, 456.112] (0.95 confidence)</td>
<td>444.685 ns [437.575, 452.366] (0.95 confidence)</td>
<td>1.0079</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d2/construct</td>
<td>225.715 ns [225.386, 226.05] (0.95 confidence)</td>
<td>222.083 ns [221.796, 222.366] (0.95 confidence)</td>
<td>1.01635</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d2/prepared</td>
<td>151.36 ns [151.15, 151.567] (0.95 confidence)</td>
<td>149.324 ns [149.071, 149.571] (0.95 confidence)</td>
<td>1.01364</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d3/adapter</td>
<td>739.501 ns [738.659, 740.35] (0.95 confidence)</td>
<td>728.521 ns [727.515, 729.511] (0.95 confidence)</td>
<td>1.01507</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d3/construct</td>
<td>414.247 ns [413.741, 414.751] (0.95 confidence)</td>
<td>404.429 ns [403.816, 405.015] (0.95 confidence)</td>
<td>1.02428</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d3/prepared</td>
<td>246.816 ns [246.469, 247.171] (0.95 confidence)</td>
<td>243.325 ns [242.939, 243.73] (0.95 confidence)</td>
<td>1.01435</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d4/adapter</td>
<td>1033.35 ns [1031.83, 1034.92] (0.95 confidence)</td>
<td>1029.76 ns [1026.02, 1033.64] (0.95 confidence)</td>
<td>1.00348</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d4/construct</td>
<td>630.57 ns [621.096, 640.846] (0.95 confidence)</td>
<td>620.731 ns [611.345, 630.797] (0.95 confidence)</td>
<td>1.01585</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d4/prepared</td>
<td>322.606 ns [322.169, 323.035] (0.95 confidence)</td>
<td>321.293 ns [320.874, 321.697] (0.95 confidence)</td>
<td>1.00409</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d5/adapter</td>
<td>1370.21 ns [1368.18, 1372.31] (0.95 confidence)</td>
<td>1359.85 ns [1357.91, 1361.79] (0.95 confidence)</td>
<td>1.00762</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d5/construct</td>
<td>848.034 ns [844.465, 851.754] (0.95 confidence)</td>
<td>830.895 ns [829.527, 832.291] (0.95 confidence)</td>
<td>1.02063</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d5/prepared</td>
<td>417.19 ns [416.754, 417.636] (0.95 confidence)</td>
<td>416.765 ns [416.3, 417.209] (0.95 confidence)</td>
<td>1.00102</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d6/adapter</td>
<td>1733.6 ns [1730.81, 1736.5] (0.95 confidence)</td>
<td>1730.13 ns [1727.3, 1733.15] (0.95 confidence)</td>
<td>1.00201</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d6/construct</td>
<td>1105.58 ns [1102.21, 1109.1] (0.95 confidence)</td>
<td>1091.71 ns [1089.45, 1094.04] (0.95 confidence)</td>
<td>1.01271</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d6/prepared</td>
<td>503.093 ns [502.297, 503.892] (0.95 confidence)</td>
<td>501.409 ns [500.223, 502.539] (0.95 confidence)</td>
<td>1.00336</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d7/adapter</td>
<td>2149.56 ns [2146.99, 2152.08] (0.95 confidence)</td>
<td>2181.4 ns [2179.8, 2182.93] (0.95 confidence)</td>
<td>0.985405</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d7/construct</td>
<td>1441.87 ns [1425.46, 1459.83] (0.95 confidence)</td>
<td>1463.29 ns [1446.58, 1481.48] (0.95 confidence)</td>
<td>0.985366</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d7/prepared</td>
<td>599.52 ns [598.765, 600.239] (0.95 confidence)</td>
<td>603.41 ns [602.623, 604.288] (0.95 confidence)</td>
<td>0.993553</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d8/adapter</td>
<td>2451.75 ns [2446.83, 2456.66] (0.95 confidence)</td>
<td>2455.65 ns [2451.69, 2459.61] (0.95 confidence)</td>
<td>0.998413</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d8/construct</td>
<td>1688.21 ns [1683.65, 1692.85] (0.95 confidence)</td>
<td>1693.44 ns [1689.93, 1697.09] (0.95 confidence)</td>
<td>0.99691</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d8/prepared</td>
<td>658.99 ns [658.01, 659.977] (0.95 confidence)</td>
<td>633.807 ns [633.112, 634.487] (0.95 confidence)</td>
<td>1.03973</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d2/adapter</td>
<td>595.251 ns [587.928, 603.011] (0.95 confidence)</td>
<td>585.393 ns [577.698, 593.519] (0.95 confidence)</td>
<td>1.01684</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d2/construct</td>
<td>304.595 ns [303.11, 305.996] (0.95 confidence)</td>
<td>299.513 ns [297.832, 301.114] (0.95 confidence)</td>
<td>1.01697</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d2/prepared</td>
<td>189.015 ns [188.635, 189.399] (0.95 confidence)</td>
<td>187.516 ns [187.021, 187.996] (0.95 confidence)</td>
<td>1.00799</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d3/adapter</td>
<td>916.749 ns [915.576, 917.975] (0.95 confidence)</td>
<td>897.274 ns [896.018, 898.45] (0.95 confidence)</td>
<td>1.0217</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d3/construct</td>
<td>534.82 ns [532.215, 537.258] (0.95 confidence)</td>
<td>521.054 ns [518.431, 523.509] (0.95 confidence)</td>
<td>1.02642</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d3/prepared</td>
<td>279.056 ns [278.708, 279.413] (0.95 confidence)</td>
<td>273.427 ns [272.951, 273.852] (0.95 confidence)</td>
<td>1.02059</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d4/adapter</td>
<td>1278.66 ns [1277.16, 1280.17] (0.95 confidence)</td>
<td>1245.85 ns [1244.26, 1247.42] (0.95 confidence)</td>
<td>1.02633</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d4/construct</td>
<td>793.667 ns [783.167, 804.869] (0.95 confidence)</td>
<td>774.827 ns [764.484, 785.831] (0.95 confidence)</td>
<td>1.02432</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d4/prepared</td>
<td>367.221 ns [366.597, 367.871] (0.95 confidence)</td>
<td>366.654 ns [366.156, 367.154] (0.95 confidence)</td>
<td>1.00155</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d5/adapter</td>
<td>1694.51 ns [1692.3, 1696.73] (0.95 confidence)</td>
<td>1675.81 ns [1672.24, 1678.63] (0.95 confidence)</td>
<td>1.01116</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d5/construct</td>
<td>1070.41 ns [1064.95, 1075.26] (0.95 confidence)</td>
<td>1068.77 ns [1062.9, 1074.34] (0.95 confidence)</td>
<td>1.00153</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d5/prepared</td>
<td>472.213 ns [471.648, 472.803] (0.95 confidence)</td>
<td>479.423 ns [478.419, 480.436] (0.95 confidence)</td>
<td>0.98496</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d6/adapter</td>
<td>2116.1 ns [2111.24, 2121.2] (0.95 confidence)</td>
<td>2068.13 ns [2063.02, 2072.86] (0.95 confidence)</td>
<td>1.02319</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d6/construct</td>
<td>1431.93 ns [1425.48, 1438.02] (0.95 confidence)</td>
<td>1355.06 ns [1348.53, 1361.14] (0.95 confidence)</td>
<td>1.05673</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d6/prepared</td>
<td>565.975 ns [564.553, 567.446] (0.95 confidence)</td>
<td>564.171 ns [563.287, 565.035] (0.95 confidence)</td>
<td>1.0032</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d7/adapter</td>
<td>2590.07 ns [2584.06, 2596.16] (0.95 confidence)</td>
<td>2571.26 ns [2568.11, 2574.27] (0.95 confidence)</td>
<td>1.00732</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d7/construct</td>
<td>1789.94 ns [1770.14, 1810.87] (0.95 confidence)</td>
<td>1779.78 ns [1762.66, 1798.06] (0.95 confidence)</td>
<td>1.0057</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d7/prepared</td>
<td>662.587 ns [661.595, 663.571] (0.95 confidence)</td>
<td>670.335 ns [669.136, 671.674] (0.95 confidence)</td>
<td>0.988442</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d8/adapter</td>
<td>2979.38 ns [2972.75, 2986.88] (0.95 confidence)</td>
<td>2976.95 ns [2971.75, 2982.29] (0.95 confidence)</td>
<td>1.00082</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d8/construct</td>
<td>2078.2 ns [2064.53, 2091.17] (0.95 confidence)</td>
<td>2083.17 ns [2073.15, 2092.43] (0.95 confidence)</td>
<td>0.997615</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d8/prepared</td>
<td>750.084 ns [747.768, 752.567] (0.95 confidence)</td>
<td>759.872 ns [758.397, 761.422] (0.95 confidence)</td>
<td>0.987118</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d2/adapter</td>
<td>540.418 ns [533.385, 547.874] (0.95 confidence)</td>
<td>530.639 ns [523.864, 537.836] (0.95 confidence)</td>
<td>1.01843</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d2/construct</td>
<td>243.874 ns [243.004, 244.797] (0.95 confidence)</td>
<td>237.463 ns [237.029, 237.95] (0.95 confidence)</td>
<td>1.027</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d2/prepared</td>
<td>219.958 ns [219.648, 220.272] (0.95 confidence)</td>
<td>215.434 ns [215.116, 215.742] (0.95 confidence)</td>
<td>1.021</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d3/adapter</td>
<td>1038.33 ns [1036.19, 1040.57] (0.95 confidence)</td>
<td>1022.62 ns [1020.77, 1024.35] (0.95 confidence)</td>
<td>1.01536</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d3/construct</td>
<td>491.93 ns [490.683, 493.174] (0.95 confidence)</td>
<td>475.252 ns [474.314, 476.255] (0.95 confidence)</td>
<td>1.03509</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d3/prepared</td>
<td>454.11 ns [453.182, 454.955] (0.95 confidence)</td>
<td>449.432 ns [448.776, 450.077] (0.95 confidence)</td>
<td>1.01041</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d4/adapter</td>
<td>1622.77 ns [1619.6, 1626.65] (0.95 confidence)</td>
<td>1601.07 ns [1598.61, 1603.33] (0.95 confidence)</td>
<td>1.01356</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d4/construct</td>
<td>896.578 ns [871.129, 923.923] (0.95 confidence)</td>
<td>869.293 ns [844.996, 895.581] (0.95 confidence)</td>
<td>1.03139</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d4/prepared</td>
<td>662.096 ns [661.306, 662.894] (0.95 confidence)</td>
<td>660.344 ns [658.889, 661.723] (0.95 confidence)</td>
<td>1.00265</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d5/adapter</td>
<td>2419.62 ns [2416.03, 2423.21] (0.95 confidence)</td>
<td>2390.76 ns [2387.85, 2393.7] (0.95 confidence)</td>
<td>1.01207</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d5/construct</td>
<td>1255.75 ns [1253.5, 1258.02] (0.95 confidence)</td>
<td>1225.61 ns [1223.25, 1227.95] (0.95 confidence)</td>
<td>1.02459</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d5/prepared</td>
<td>984.063 ns [982.684, 985.445] (0.95 confidence)</td>
<td>980.312 ns [978.38, 982.049] (0.95 confidence)</td>
<td>1.00383</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d6/adapter</td>
<td>3253.58 ns [3249.67, 3257.6] (0.95 confidence)</td>
<td>3306.36 ns [3301.77, 3310.73] (0.95 confidence)</td>
<td>0.984036</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d6/construct</td>
<td>1697.23 ns [1693.28, 1701.57] (0.95 confidence)</td>
<td>1716.51 ns [1711.02, 1722.33] (0.95 confidence)</td>
<td>0.988771</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d6/prepared</td>
<td>1323.83 ns [1321.38, 1326.3] (0.95 confidence)</td>
<td>1329.79 ns [1328.18, 1331.22] (0.95 confidence)</td>
<td>0.995516</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d7/adapter</td>
<td>4291.49 ns [4284.38, 4298.27] (0.95 confidence)</td>
<td>4297.8 ns [4293.8, 4301.85] (0.95 confidence)</td>
<td>0.998534</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d7/construct</td>
<td>2414.14 ns [2347.27, 2486.72] (0.95 confidence)</td>
<td>2433.23 ns [2365.09, 2506.64] (0.95 confidence)</td>
<td>0.992154</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d7/prepared</td>
<td>1744.02 ns [1742, 1746.06] (0.95 confidence)</td>
<td>1746.66 ns [1744.42, 1748.95] (0.95 confidence)</td>
<td>0.998493</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d8/adapter</td>
<td>5153.02 ns [5145.79, 5159.99] (0.95 confidence)</td>
<td>5249.16 ns [5242.64, 5255.49] (0.95 confidence)</td>
<td>0.981686</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d8/construct</td>
<td>2878.54 ns [2871.64, 2887.52] (0.95 confidence)</td>
<td>2890.09 ns [2884.64, 2895.81] (0.95 confidence)</td>
<td>0.996003</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d8/prepared</td>
<td>2103.11 ns [2099.96, 2106.29] (0.95 confidence)</td>
<td>2120.31 ns [2116.71, 2124.01] (0.95 confidence)</td>
<td>0.991888</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d2/adapter</td>
<td>670.075 ns [662.734, 677.885] (0.95 confidence)</td>
<td>657.962 ns [650.797, 665.547] (0.95 confidence)</td>
<td>1.01841</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d2/construct</td>
<td>342.253 ns [340.443, 344.018] (0.95 confidence)</td>
<td>345.136 ns [343.331, 346.93] (0.95 confidence)</td>
<td>0.991647</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d2/prepared</td>
<td>213.479 ns [213.101, 213.885] (0.95 confidence)</td>
<td>212.55 ns [212.252, 212.842] (0.95 confidence)</td>
<td>1.00437</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d3/adapter</td>
<td>1155.69 ns [1153.82, 1157.51] (0.95 confidence)</td>
<td>1129.34 ns [1127.61, 1131] (0.95 confidence)</td>
<td>1.02333</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d3/construct</td>
<td>664.352 ns [661.652, 666.997] (0.95 confidence)</td>
<td>652.7 ns [649.867, 655.366] (0.95 confidence)</td>
<td>1.01785</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d3/prepared</td>
<td>368.206 ns [367.694, 368.718] (0.95 confidence)</td>
<td>361.423 ns [360.919, 361.899] (0.95 confidence)</td>
<td>1.01877</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d4/adapter</td>
<td>1838.53 ns [1836.08, 1841.05] (0.95 confidence)</td>
<td>1783.04 ns [1780.78, 1785.32] (0.95 confidence)</td>
<td>1.03112</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d4/construct</td>
<td>1125.23 ns [1099.19, 1153.52] (0.95 confidence)</td>
<td>1099.46 ns [1074.01, 1126.96] (0.95 confidence)</td>
<td>1.02344</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d4/prepared</td>
<td>592.822 ns [592.039, 593.635] (0.95 confidence)</td>
<td>587.023 ns [586.31, 587.717] (0.95 confidence)</td>
<td>1.00988</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d5/adapter</td>
<td>2631.67 ns [2626.13, 2637.73] (0.95 confidence)</td>
<td>2571.56 ns [2567.49, 2575.35] (0.95 confidence)</td>
<td>1.02338</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d5/construct</td>
<td>1568.88 ns [1557.98, 1579.32] (0.95 confidence)</td>
<td>1554.58 ns [1544.81, 1563.76] (0.95 confidence)</td>
<td>1.0092</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d5/prepared</td>
<td>815.397 ns [813.99, 816.861] (0.95 confidence)</td>
<td>810.074 ns [808.298, 811.922] (0.95 confidence)</td>
<td>1.00657</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d6/adapter</td>
<td>3545.32 ns [3540.65, 3550.04] (0.95 confidence)</td>
<td>3505.89 ns [3501.17, 3510.57] (0.95 confidence)</td>
<td>1.01125</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d6/construct</td>
<td>2250.06 ns [2238.17, 2261.09] (0.95 confidence)</td>
<td>2261.6 ns [2253.44, 2268.98] (0.95 confidence)</td>
<td>0.994895</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d6/prepared</td>
<td>1030.36 ns [1027.12, 1033.79] (0.95 confidence)</td>
<td>1036.39 ns [1033.41, 1039.4] (0.95 confidence)</td>
<td>0.994189</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d7/adapter</td>
<td>4686.28 ns [4678.54, 4694.61] (0.95 confidence)</td>
<td>4644.8 ns [4640.14, 4649.43] (0.95 confidence)</td>
<td>1.00893</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d7/construct</td>
<td>3149.47 ns [3069.36, 3236.88] (0.95 confidence)</td>
<td>3156.59 ns [3083.35, 3236.7] (0.95 confidence)</td>
<td>0.997745</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d7/prepared</td>
<td>1354.85 ns [1351.81, 1357.89] (0.95 confidence)</td>
<td>1368.48 ns [1366.17, 1370.89] (0.95 confidence)</td>
<td>0.99004</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d8/adapter</td>
<td>5760.84 ns [5750.19, 5771.38] (0.95 confidence)</td>
<td>5711.94 ns [5702.5, 5721.93] (0.95 confidence)</td>
<td>1.00856</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d8/construct</td>
<td>3801.21 ns [3781.1, 3820.28] (0.95 confidence)</td>
<td>3825.27 ns [3807.3, 3843.33] (0.95 confidence)</td>
<td>0.993711</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d8/prepared</td>
<td>1614.51 ns [1609.85, 1619.37] (0.95 confidence)</td>
<td>1629.31 ns [1624.42, 1634.43] (0.95 confidence)</td>
<td>0.990915</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d2/adapter</td>
<td>541.04 ns [534.082, 548.441] (0.95 confidence)</td>
<td>533.187 ns [526.389, 540.381] (0.95 confidence)</td>
<td>1.01473</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d2/construct</td>
<td>240.538 ns [240.041, 241.109] (0.95 confidence)</td>
<td>239.103 ns [238.759, 239.467] (0.95 confidence)</td>
<td>1.006</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d2/prepared</td>
<td>220.575 ns [220.157, 221.002] (0.95 confidence)</td>
<td>214.646 ns [214, 215.162] (0.95 confidence)</td>
<td>1.02762</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d3/adapter</td>
<td>1043.73 ns [1040.38, 1047.46] (0.95 confidence)</td>
<td>1030.44 ns [1028.49, 1032.37] (0.95 confidence)</td>
<td>1.01289</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d3/construct</td>
<td>489.59 ns [488.904, 490.312] (0.95 confidence)</td>
<td>478.728 ns [477.285, 480.21] (0.95 confidence)</td>
<td>1.02269</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d3/prepared</td>
<td>451.348 ns [450.801, 451.902] (0.95 confidence)</td>
<td>447.399 ns [446.872, 447.918] (0.95 confidence)</td>
<td>1.00883</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d4/adapter</td>
<td>1602.16 ns [1599.99, 1604.4] (0.95 confidence)</td>
<td>1606.52 ns [1604.15, 1608.92] (0.95 confidence)</td>
<td>0.997287</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d4/construct</td>
<td>881.403 ns [856.825, 908.161] (0.95 confidence)</td>
<td>866.935 ns [842.212, 894.071] (0.95 confidence)</td>
<td>1.01669</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d4/prepared</td>
<td>655.514 ns [654.821, 656.249] (0.95 confidence)</td>
<td>656.026 ns [654.371, 657.819] (0.95 confidence)</td>
<td>0.999219</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d5/adapter</td>
<td>2418.64 ns [2414.2, 2423.3] (0.95 confidence)</td>
<td>2391.55 ns [2384.4, 2397.24] (0.95 confidence)</td>
<td>1.01133</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d5/construct</td>
<td>1229.74 ns [1225.64, 1234.08] (0.95 confidence)</td>
<td>1222.14 ns [1219.37, 1225.25] (0.95 confidence)</td>
<td>1.00622</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d5/prepared</td>
<td>975.437 ns [974.528, 976.305] (0.95 confidence)</td>
<td>980.874 ns [979.298, 982.43] (0.95 confidence)</td>
<td>0.994457</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d6/adapter</td>
<td>3256.73 ns [3252.08, 3262.04] (0.95 confidence)</td>
<td>3326.12 ns [3320.91, 3331.01] (0.95 confidence)</td>
<td>0.979141</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d6/construct</td>
<td>1696.68 ns [1694.02, 1699.71] (0.95 confidence)</td>
<td>1723.86 ns [1718.75, 1729.36] (0.95 confidence)</td>
<td>0.984232</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d6/prepared</td>
<td>1322.64 ns [1320.39, 1324.89] (0.95 confidence)</td>
<td>1336.15 ns [1334.53, 1337.8] (0.95 confidence)</td>
<td>0.989891</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d7/adapter</td>
<td>4306.51 ns [4300.87, 4312.06] (0.95 confidence)</td>
<td>4287.26 ns [4281.09, 4293.48] (0.95 confidence)</td>
<td>1.00449</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d7/construct</td>
<td>2400.35 ns [2333.56, 2473.09] (0.95 confidence)</td>
<td>2435.05 ns [2363.81, 2512.23] (0.95 confidence)</td>
<td>0.985748</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d7/prepared</td>
<td>1744.18 ns [1741.94, 1746.36] (0.95 confidence)</td>
<td>1736.94 ns [1733.84, 1739.87] (0.95 confidence)</td>
<td>1.00417</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d8/adapter</td>
<td>5199.54 ns [5184.57, 5216.85] (0.95 confidence)</td>
<td>5269.74 ns [5262.11, 5276.66] (0.95 confidence)</td>
<td>0.986679</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d8/construct</td>
<td>2883.53 ns [2875.94, 2891.15] (0.95 confidence)</td>
<td>2871.83 ns [2865.35, 2878.29] (0.95 confidence)</td>
<td>1.00407</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d8/prepared</td>
<td>2125.55 ns [2121.91, 2129.69] (0.95 confidence)</td>
<td>2140.29 ns [2136.92, 2143.8] (0.95 confidence)</td>
<td>0.993115</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d2/adapter</td>
<td>669.353 ns [661.498, 677.691] (0.95 confidence)</td>
<td>661.233 ns [654.25, 668.681] (0.95 confidence)</td>
<td>1.01228</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d2/construct</td>
<td>350.475 ns [349.071, 351.884] (0.95 confidence)</td>
<td>345.689 ns [344.278, 347.064] (0.95 confidence)</td>
<td>1.01385</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d2/prepared</td>
<td>215.072 ns [214.781, 215.358] (0.95 confidence)</td>
<td>221.206 ns [219.889, 222.594] (0.95 confidence)</td>
<td>0.972269</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d3/adapter</td>
<td>1153.9 ns [1151.19, 1156.92] (0.95 confidence)</td>
<td>1125.89 ns [1124.52, 1127.2] (0.95 confidence)</td>
<td>1.02487</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d3/construct</td>
<td>663.243 ns [660.25, 666.118] (0.95 confidence)</td>
<td>649.327 ns [646.081, 652.433] (0.95 confidence)</td>
<td>1.02143</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d3/prepared</td>
<td>365.562 ns [364.995, 366.136] (0.95 confidence)</td>
<td>361.428 ns [360.923, 361.933] (0.95 confidence)</td>
<td>1.01144</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d4/adapter</td>
<td>1823.35 ns [1819.97, 1826.98] (0.95 confidence)</td>
<td>1785.69 ns [1782.73, 1788.73] (0.95 confidence)</td>
<td>1.02109</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d4/construct</td>
<td>1128.93 ns [1100.64, 1159.25] (0.95 confidence)</td>
<td>1106.6 ns [1081.76, 1133.27] (0.95 confidence)</td>
<td>1.02018</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d4/prepared</td>
<td>597.033 ns [595.948, 598.116] (0.95 confidence)</td>
<td>589.271 ns [588.572, 589.947] (0.95 confidence)</td>
<td>1.01317</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d5/adapter</td>
<td>2597.55 ns [2593.64, 2601.48] (0.95 confidence)</td>
<td>2562.49 ns [2558.82, 2565.95] (0.95 confidence)</td>
<td>1.01368</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d5/construct</td>
<td>1564.82 ns [1555.68, 1573.27] (0.95 confidence)</td>
<td>1560.09 ns [1552.81, 1566.56] (0.95 confidence)</td>
<td>1.00303</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d5/prepared</td>
<td>811.369 ns [810.031, 812.748] (0.95 confidence)</td>
<td>808.391 ns [806.872, 809.943] (0.95 confidence)</td>
<td>1.00368</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d6/adapter</td>
<td>3540.22 ns [3535.31, 3545.13] (0.95 confidence)</td>
<td>3532.67 ns [3528.16, 3537.03] (0.95 confidence)</td>
<td>1.00214</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d6/construct</td>
<td>2246.46 ns [2235.91, 2255.96] (0.95 confidence)</td>
<td>2259.34 ns [2247.02, 2270.11] (0.95 confidence)</td>
<td>0.994303</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d6/prepared</td>
<td>1023.36 ns [1020.35, 1026.42] (0.95 confidence)</td>
<td>1034.14 ns [1031.19, 1037.17] (0.95 confidence)</td>
<td>0.989575</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d7/adapter</td>
<td>4672 ns [4664.06, 4679.85] (0.95 confidence)</td>
<td>4585.93 ns [4580.22, 4591.44] (0.95 confidence)</td>
<td>1.01877</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d7/construct</td>
<td>3111.7 ns [3044.76, 3185.41] (0.95 confidence)</td>
<td>3127.9 ns [3055.5, 3208.07] (0.95 confidence)</td>
<td>0.994819</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d7/prepared</td>
<td>1361.27 ns [1356.66, 1366.44] (0.95 confidence)</td>
<td>1369.64 ns [1367.18, 1372.17] (0.95 confidence)</td>
<td>0.993886</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d8/adapter</td>
<td>5797.8 ns [5787.54, 5808.22] (0.95 confidence)</td>
<td>5728.2 ns [5719.97, 5736.35] (0.95 confidence)</td>
<td>1.01215</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d8/construct</td>
<td>3794.16 ns [3770.88, 3815.62] (0.95 confidence)</td>
<td>3820.64 ns [3809.17, 3831.49] (0.95 confidence)</td>
<td>0.993068</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d8/prepared</td>
<td>1626.98 ns [1622.51, 1631.64] (0.95 confidence)</td>
<td>1638.41 ns [1634.17, 1642.72] (0.95 confidence)</td>
<td>0.993022</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d2/adapter</td>
<td>546.323 ns [539.468, 553.594] (0.95 confidence)</td>
<td>540.067 ns [533.436, 547.121] (0.95 confidence)</td>
<td>1.01158</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d2/construct</td>
<td>239.793 ns [239.403, 240.235] (0.95 confidence)</td>
<td>233.634 ns [233.281, 233.973] (0.95 confidence)</td>
<td>1.02636</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d2/prepared</td>
<td>220.459 ns [220.156, 220.781] (0.95 confidence)</td>
<td>216.185 ns [215.785, 216.515] (0.95 confidence)</td>
<td>1.01977</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d3/adapter</td>
<td>1015.82 ns [1014.41, 1017.28] (0.95 confidence)</td>
<td>1007.99 ns [1006.53, 1009.38] (0.95 confidence)</td>
<td>1.00776</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d3/construct</td>
<td>475.032 ns [474.077, 476.002] (0.95 confidence)</td>
<td>467.862 ns [466.645, 469.236] (0.95 confidence)</td>
<td>1.01533</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d3/prepared</td>
<td>442.239 ns [441.525, 443.012] (0.95 confidence)</td>
<td>432.552 ns [431.829, 433.308] (0.95 confidence)</td>
<td>1.02239</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d4/adapter</td>
<td>1648.88 ns [1645.69, 1652.1] (0.95 confidence)</td>
<td>1633.28 ns [1630.95, 1635.59] (0.95 confidence)</td>
<td>1.00955</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d4/construct</td>
<td>874.597 ns [854.809, 896.072] (0.95 confidence)</td>
<td>859.115 ns [839.317, 880.533] (0.95 confidence)</td>
<td>1.01802</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d4/prepared</td>
<td>691.441 ns [690.643, 692.233] (0.95 confidence)</td>
<td>690.373 ns [688.737, 692.011] (0.95 confidence)</td>
<td>1.00155</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d5/adapter</td>
<td>2329.92 ns [2324.51, 2335.59] (0.95 confidence)</td>
<td>2333.83 ns [2328.72, 2338.95] (0.95 confidence)</td>
<td>0.998327</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d5/construct</td>
<td>1198.19 ns [1193.37, 1203.49] (0.95 confidence)</td>
<td>1200.79 ns [1195.58, 1206.35] (0.95 confidence)</td>
<td>0.997833</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d5/prepared</td>
<td>946.636 ns [945.78, 947.507] (0.95 confidence)</td>
<td>943.053 ns [941.848, 944.215] (0.95 confidence)</td>
<td>1.0038</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d6/adapter</td>
<td>3322.35 ns [3314.82, 3329.77] (0.95 confidence)</td>
<td>3396.11 ns [3391.92, 3400.38] (0.95 confidence)</td>
<td>0.978282</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d6/construct</td>
<td>1764.23 ns [1759.13, 1769.68] (0.95 confidence)</td>
<td>1780.9 ns [1776.82, 1785.36] (0.95 confidence)</td>
<td>0.990641</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d6/prepared</td>
<td>1358.68 ns [1357.56, 1359.78] (0.95 confidence)</td>
<td>1363.71 ns [1361.92, 1365.39] (0.95 confidence)</td>
<td>0.996311</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d7/adapter</td>
<td>4296.6 ns [4291.21, 4301.89] (0.95 confidence)</td>
<td>4386.52 ns [4376.46, 4400.85] (0.95 confidence)</td>
<td>0.9795</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d7/construct</td>
<td>2449.27 ns [2386.88, 2515.99] (0.95 confidence)</td>
<td>2476.21 ns [2412.36, 2545.58] (0.95 confidence)</td>
<td>0.989121</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d7/prepared</td>
<td>1773.67 ns [1771.71, 1775.56] (0.95 confidence)</td>
<td>1773.38 ns [1771.25, 1775.73] (0.95 confidence)</td>
<td>1.00017</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d8/adapter</td>
<td>5397.8 ns [5386.56, 5409.14] (0.95 confidence)</td>
<td>5381.35 ns [5371.02, 5391.86] (0.95 confidence)</td>
<td>1.00306</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d8/construct</td>
<td>3029.28 ns [3023.39, 3035.58] (0.95 confidence)</td>
<td>3044.26 ns [3035.94, 3052.96] (0.95 confidence)</td>
<td>0.995079</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d8/prepared</td>
<td>2081.08 ns [2078.06, 2084.28] (0.95 confidence)</td>
<td>2081.62 ns [2078.85, 2084.38] (0.95 confidence)</td>
<td>0.999736</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d2/adapter</td>
<td>705.075 ns [697.716, 712.969] (0.95 confidence)</td>
<td>691.797 ns [685.031, 698.883] (0.95 confidence)</td>
<td>1.01919</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d2/construct</td>
<td>373.648 ns [371.992, 375.155] (0.95 confidence)</td>
<td>367.084 ns [365.499, 368.566] (0.95 confidence)</td>
<td>1.01788</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d2/prepared</td>
<td>228.892 ns [228.423, 229.354] (0.95 confidence)</td>
<td>226.09 ns [225.735, 226.443] (0.95 confidence)</td>
<td>1.01239</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d3/adapter</td>
<td>1183.05 ns [1181.3, 1184.84] (0.95 confidence)</td>
<td>1156.18 ns [1154.08, 1158.21] (0.95 confidence)</td>
<td>1.02325</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d3/construct</td>
<td>685.323 ns [683.035, 687.449] (0.95 confidence)</td>
<td>671.479 ns [668.681, 674.092] (0.95 confidence)</td>
<td>1.02062</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d3/prepared</td>
<td>385.025 ns [384.471, 385.585] (0.95 confidence)</td>
<td>378.52 ns [377.849, 379.239] (0.95 confidence)</td>
<td>1.01719</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d4/adapter</td>
<td>1771.08 ns [1769.12, 1773.16] (0.95 confidence)</td>
<td>1730.18 ns [1727.72, 1732.64] (0.95 confidence)</td>
<td>1.02364</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d4/construct</td>
<td>1056.27 ns [1034.7, 1079.55] (0.95 confidence)</td>
<td>1015.54 ns [994.78, 1037.73] (0.95 confidence)</td>
<td>1.04011</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d4/prepared</td>
<td>624.89 ns [623.968, 625.855] (0.95 confidence)</td>
<td>619.238 ns [618.348, 620.098] (0.95 confidence)</td>
<td>1.00913</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d5/adapter</td>
<td>2495.46 ns [2493.01, 2497.77] (0.95 confidence)</td>
<td>2469.49 ns [2465.96, 2472.91] (0.95 confidence)</td>
<td>1.01052</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d5/construct</td>
<td>1430.17 ns [1422.37, 1437.04] (0.95 confidence)</td>
<td>1408.32 ns [1400.05, 1415.87] (0.95 confidence)</td>
<td>1.01552</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d5/prepared</td>
<td>869.025 ns [868.021, 870.088] (0.95 confidence)</td>
<td>874.367 ns [872.907, 875.944] (0.95 confidence)</td>
<td>0.993891</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d6/adapter</td>
<td>3366.07 ns [3360.21, 3371.95] (0.95 confidence)</td>
<td>3357.95 ns [3351.67, 3364.43] (0.95 confidence)</td>
<td>1.00242</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d6/construct</td>
<td>1999.47 ns [1989.49, 2008.85] (0.95 confidence)</td>
<td>2011.57 ns [2004.08, 2018.95] (0.95 confidence)</td>
<td>0.993983</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d6/prepared</td>
<td>1125.7 ns [1124.16, 1127.27] (0.95 confidence)</td>
<td>1137.32 ns [1135.64, 1139.01] (0.95 confidence)</td>
<td>0.989789</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d7/adapter</td>
<td>4429.87 ns [4422.28, 4437.41] (0.95 confidence)</td>
<td>4388.1 ns [4381.42, 4394.01] (0.95 confidence)</td>
<td>1.00952</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d7/construct</td>
<td>2821.5 ns [2766.89, 2881.34] (0.95 confidence)</td>
<td>2812.04 ns [2757.65, 2871.44] (0.95 confidence)</td>
<td>1.00337</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d7/prepared</td>
<td>1454.64 ns [1452.66, 1456.45] (0.95 confidence)</td>
<td>1466.78 ns [1464.86, 1468.68] (0.95 confidence)</td>
<td>0.991726</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d8/adapter</td>
<td>5614.21 ns [5606.87, 5621.36] (0.95 confidence)</td>
<td>5635.81 ns [5627.67, 5643.91] (0.95 confidence)</td>
<td>0.996167</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d8/construct</td>
<td>3469.69 ns [3454.43, 3483.8] (0.95 confidence)</td>
<td>3463.72 ns [3450.2, 3476.8] (0.95 confidence)</td>
<td>1.00172</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d8/prepared</td>
<td>1849.63 ns [1846.73, 1852.61] (0.95 confidence)</td>
<td>1862.7 ns [1859.65, 1865.73] (0.95 confidence)</td>
<td>0.992982</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d2/adapter</td>
<td>1056.22 ns [1041.45, 1071.99] (0.95 confidence)</td>
<td>1037.55 ns [1023.41, 1052.77] (0.95 confidence)</td>
<td>1.01799</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d2/construct</td>
<td>294.148 ns [293.678, 294.625] (0.95 confidence)</td>
<td>286.736 ns [286.29, 287.201] (0.95 confidence)</td>
<td>1.02585</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d2/prepared</td>
<td>635.846 ns [634.869, 636.987] (0.95 confidence)</td>
<td>628.08 ns [627.295, 628.846] (0.95 confidence)</td>
<td>1.01236</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d3/adapter</td>
<td>1815.1 ns [1811.87, 1818.46] (0.95 confidence)</td>
<td>1782.59 ns [1779.88, 1785.3] (0.95 confidence)</td>
<td>1.01824</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d3/construct</td>
<td>583.448 ns [582.424, 584.466] (0.95 confidence)</td>
<td>570.236 ns [569.118, 571.402] (0.95 confidence)</td>
<td>1.02317</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d3/prepared</td>
<td>1101.74 ns [1100.55, 1102.94] (0.95 confidence)</td>
<td>1076.08 ns [1074.66, 1077.53] (0.95 confidence)</td>
<td>1.02385</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d4/adapter</td>
<td>2689.56 ns [2686.05, 2693.03] (0.95 confidence)</td>
<td>2677.48 ns [2673.27, 2681.67] (0.95 confidence)</td>
<td>1.00451</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d4/construct</td>
<td>853.536 ns [833.446, 875.069] (0.95 confidence)</td>
<td>839.309 ns [818.917, 861.33] (0.95 confidence)</td>
<td>1.01695</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d4/prepared</td>
<td>1740.82 ns [1738.57, 1743.26] (0.95 confidence)</td>
<td>1726.05 ns [1723.61, 1728.38] (0.95 confidence)</td>
<td>1.00855</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d5/adapter</td>
<td>4126.93 ns [4119.38, 4134.48] (0.95 confidence)</td>
<td>4136.66 ns [4132.09, 4140.98] (0.95 confidence)</td>
<td>0.997648</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d5/construct</td>
<td>1315.45 ns [1311.43, 1319.77] (0.95 confidence)</td>
<td>1320.05 ns [1314.61, 1325.75] (0.95 confidence)</td>
<td>0.996518</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d5/prepared</td>
<td>2613.24 ns [2609.06, 2617.46] (0.95 confidence)</td>
<td>2615.22 ns [2611.11, 2619.33] (0.95 confidence)</td>
<td>0.999245</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d6/adapter</td>
<td>5501.22 ns [5493.8, 5508.68] (0.95 confidence)</td>
<td>5577.09 ns [5567.11, 5585.82] (0.95 confidence)</td>
<td>0.986395</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d6/construct</td>
<td>1681.15 ns [1674.23, 1688.4] (0.95 confidence)</td>
<td>1677.12 ns [1673.4, 1680.98] (0.95 confidence)</td>
<td>1.00241</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d6/prepared</td>
<td>3598.21 ns [3594.28, 3602.03] (0.95 confidence)</td>
<td>3627.96 ns [3622.83, 3633.2] (0.95 confidence)</td>
<td>0.9918</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d7/adapter</td>
<td>7714.84 ns [7704.91, 7724.95] (0.95 confidence)</td>
<td>7721.59 ns [7710.02, 7733.2] (0.95 confidence)</td>
<td>0.999127</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d7/construct</td>
<td>2273.15 ns [2223.88, 2327.09] (0.95 confidence)</td>
<td>2268.65 ns [2226.51, 2314.93] (0.95 confidence)</td>
<td>1.00198</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d7/prepared</td>
<td>5211.71 ns [5202.97, 5220.41] (0.95 confidence)</td>
<td>5189.78 ns [5180.22, 5201.02] (0.95 confidence)</td>
<td>1.00422</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d8/adapter</td>
<td>9400.38 ns [9388.37, 9411.81] (0.95 confidence)</td>
<td>9342.26 ns [9332.54, 9351.77] (0.95 confidence)</td>
<td>1.00622</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d8/construct</td>
<td>2363.38 ns [2358.38, 2368.52] (0.95 confidence)</td>
<td>2374.81 ns [2368.67, 2381.05] (0.95 confidence)</td>
<td>0.995185</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d8/prepared</td>
<td>6892.71 ns [6884.36, 6900.86] (0.95 confidence)</td>
<td>6849.86 ns [6840.92, 6858.41] (0.95 confidence)</td>
<td>1.00626</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d2/adapter</td>
<td>1191.85 ns [1183.92, 1200.22] (0.95 confidence)</td>
<td>1168.63 ns [1161.42, 1176.17] (0.95 confidence)</td>
<td>1.01987</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d2/construct</td>
<td>475.732 ns [473.154, 478.045] (0.95 confidence)</td>
<td>466.661 ns [464.266, 468.748] (0.95 confidence)</td>
<td>1.01944</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d2/prepared</td>
<td>586.524 ns [585.737, 587.354] (0.95 confidence)</td>
<td>577.577 ns [576.6, 578.546] (0.95 confidence)</td>
<td>1.01549</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d3/adapter</td>
<td>1922.73 ns [1920.2, 1925.06] (0.95 confidence)</td>
<td>1888.39 ns [1884.99, 1892.27] (0.95 confidence)</td>
<td>1.01818</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d3/construct</td>
<td>740.519 ns [737.918, 743.08] (0.95 confidence)</td>
<td>729.674 ns [727.602, 731.576] (0.95 confidence)</td>
<td>1.01486</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d3/prepared</td>
<td>1029.39 ns [1028.44, 1030.34] (0.95 confidence)</td>
<td>1014.7 ns [1013.23, 1016.11] (0.95 confidence)</td>
<td>1.01448</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d4/adapter</td>
<td>3107.58 ns [3103.05, 3112.36] (0.95 confidence)</td>
<td>3046.24 ns [3040.45, 3051.35] (0.95 confidence)</td>
<td>1.02014</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d4/construct</td>
<td>1326.71 ns [1299.34, 1356.35] (0.95 confidence)</td>
<td>1302.65 ns [1277.69, 1329.51] (0.95 confidence)</td>
<td>1.01847</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d4/prepared</td>
<td>1643.9 ns [1641.77, 1646.19] (0.95 confidence)</td>
<td>1623.81 ns [1621.6, 1625.79] (0.95 confidence)</td>
<td>1.01237</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d5/adapter</td>
<td>4323.87 ns [4318.02, 4329.42] (0.95 confidence)</td>
<td>4282.28 ns [4276.64, 4287.98] (0.95 confidence)</td>
<td>1.00971</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d5/construct</td>
<td>1639.06 ns [1634.05, 1643.89] (0.95 confidence)</td>
<td>1642.05 ns [1636.32, 1647.88] (0.95 confidence)</td>
<td>0.998182</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d5/prepared</td>
<td>2450.19 ns [2447.41, 2453.03] (0.95 confidence)</td>
<td>2438.48 ns [2434.77, 2442.28] (0.95 confidence)</td>
<td>1.0048</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d6/adapter</td>
<td>5877.09 ns [5867.64, 5886.41] (0.95 confidence)</td>
<td>5915.33 ns [5909.91, 5920.57] (0.95 confidence)</td>
<td>0.993536</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d6/construct</td>
<td>2235.06 ns [2226.68, 2243.14] (0.95 confidence)</td>
<td>2247.55 ns [2242.69, 2252.36] (0.95 confidence)</td>
<td>0.994443</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d6/prepared</td>
<td>3406.39 ns [3401.66, 3411.24] (0.95 confidence)</td>
<td>3434.54 ns [3430.82, 3438.29] (0.95 confidence)</td>
<td>0.991807</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d7/adapter</td>
<td>8447.16 ns [8433.73, 8459.41] (0.95 confidence)</td>
<td>8422.74 ns [8409.58, 8435.45] (0.95 confidence)</td>
<td>1.0029</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d7/construct</td>
<td>2928.88 ns [2878.95, 2983.05] (0.95 confidence)</td>
<td>2909.07 ns [2859.34, 2963.12] (0.95 confidence)</td>
<td>1.00681</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d7/prepared</td>
<td>5381.95 ns [5369.09, 5398.32] (0.95 confidence)</td>
<td>5359.03 ns [5346.8, 5368.32] (0.95 confidence)</td>
<td>1.00428</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d8/adapter</td>
<td>10333 ns [10320.5, 10344.9] (0.95 confidence)</td>
<td>10315.8 ns [10295.2, 10336.7] (0.95 confidence)</td>
<td>1.00166</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d8/construct</td>
<td>3353.66 ns [3347.12, 3360.25] (0.95 confidence)</td>
<td>3337.73 ns [3329.01, 3345.96] (0.95 confidence)</td>
<td>1.00477</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d8/prepared</td>
<td>6745.6 ns [6731.55, 6757.72] (0.95 confidence)</td>
<td>6703.37 ns [6691.5, 6715.55] (0.95 confidence)</td>
<td>1.0063</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d2/adapter</td>
<td>538.907 ns [531.687, 546.483] (0.95 confidence)</td>
<td>538.805 ns [531.767, 546.259] (0.95 confidence)</td>
<td>1.00019</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d2/construct</td>
<td>239.191 ns [238.848, 239.525] (0.95 confidence)</td>
<td>234.636 ns [234.321, 234.943] (0.95 confidence)</td>
<td>1.01941</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d2/prepared</td>
<td>217.499 ns [217.241, 217.77] (0.95 confidence)</td>
<td>216.021 ns [215.727, 216.363] (0.95 confidence)</td>
<td>1.00684</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d3/adapter</td>
<td>1015.12 ns [1013.13, 1017.18] (0.95 confidence)</td>
<td>1008.82 ns [1006.86, 1010.69] (0.95 confidence)</td>
<td>1.00624</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d3/construct</td>
<td>476.56 ns [475.623, 477.551] (0.95 confidence)</td>
<td>471.531 ns [470.771, 472.301] (0.95 confidence)</td>
<td>1.01067</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d3/prepared</td>
<td>440.667 ns [440.129, 441.264] (0.95 confidence)</td>
<td>434.625 ns [434.122, 435.095] (0.95 confidence)</td>
<td>1.0139</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d4/adapter</td>
<td>1640.1 ns [1637.85, 1642.39] (0.95 confidence)</td>
<td>1644.38 ns [1641.63, 1647.08] (0.95 confidence)</td>
<td>0.997398</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d4/construct</td>
<td>876.285 ns [856.115, 898.208] (0.95 confidence)</td>
<td>861.674 ns [842.888, 882.229] (0.95 confidence)</td>
<td>1.01696</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d4/prepared</td>
<td>688.202 ns [687.427, 688.958] (0.95 confidence)</td>
<td>691.652 ns [689.722, 693.614] (0.95 confidence)</td>
<td>0.995012</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d5/adapter</td>
<td>2346.73 ns [2342.61, 2350.89] (0.95 confidence)</td>
<td>2333.45 ns [2330.75, 2336.19] (0.95 confidence)</td>
<td>1.00569</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d5/construct</td>
<td>1228.01 ns [1224.41, 1231.83] (0.95 confidence)</td>
<td>1212.14 ns [1206.67, 1217.82] (0.95 confidence)</td>
<td>1.01309</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d5/prepared</td>
<td>952.984 ns [951.787, 954.199] (0.95 confidence)</td>
<td>946.259 ns [945.254, 947.269] (0.95 confidence)</td>
<td>1.00711</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d6/adapter</td>
<td>3327.64 ns [3323.06, 3332.28] (0.95 confidence)</td>
<td>3396.47 ns [3391.39, 3401.42] (0.95 confidence)</td>
<td>0.979735</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d6/construct</td>
<td>1763.33 ns [1757.55, 1769.33] (0.95 confidence)</td>
<td>1764.14 ns [1759.7, 1768.95] (0.95 confidence)</td>
<td>0.999542</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d6/prepared</td>
<td>1359.86 ns [1358.5, 1361.16] (0.95 confidence)</td>
<td>1362.85 ns [1360.75, 1364.75] (0.95 confidence)</td>
<td>0.997802</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d7/adapter</td>
<td>4318.04 ns [4310.79, 4325.26] (0.95 confidence)</td>
<td>4390.66 ns [4387.07, 4394.25] (0.95 confidence)</td>
<td>0.983459</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d7/construct</td>
<td>2495.97 ns [2431.16, 2566.22] (0.95 confidence)</td>
<td>2480.9 ns [2415.36, 2550.01] (0.95 confidence)</td>
<td>1.00607</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d7/prepared</td>
<td>1776.52 ns [1774.46, 1778.52] (0.95 confidence)</td>
<td>1787.05 ns [1784.34, 1790.04] (0.95 confidence)</td>
<td>0.994104</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d8/adapter</td>
<td>5361.52 ns [5351.11, 5372.27] (0.95 confidence)</td>
<td>5380.94 ns [5374.11, 5388.06] (0.95 confidence)</td>
<td>0.996392</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d8/construct</td>
<td>3049.05 ns [3038.93, 3059.63] (0.95 confidence)</td>
<td>3053.52 ns [3041.91, 3065.24] (0.95 confidence)</td>
<td>0.998538</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d8/prepared</td>
<td>2080.53 ns [2077.89, 2083.2] (0.95 confidence)</td>
<td>2084.44 ns [2081.51, 2087.24] (0.95 confidence)</td>
<td>0.998125</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d2/adapter</td>
<td>701.144 ns [693.825, 708.957] (0.95 confidence)</td>
<td>692.717 ns [685.778, 700.091] (0.95 confidence)</td>
<td>1.01217</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d2/construct</td>
<td>368.938 ns [366.999, 370.782] (0.95 confidence)</td>
<td>369.116 ns [366.959, 371.142] (0.95 confidence)</td>
<td>0.999517</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d2/prepared</td>
<td>226.454 ns [226.012, 226.925] (0.95 confidence)</td>
<td>225.701 ns [225.319, 226.075] (0.95 confidence)</td>
<td>1.00334</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d3/adapter</td>
<td>1189.45 ns [1188.07, 1190.79] (0.95 confidence)</td>
<td>1156.42 ns [1154.67, 1158.06] (0.95 confidence)</td>
<td>1.02856</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d3/construct</td>
<td>690.107 ns [686.901, 693.134] (0.95 confidence)</td>
<td>674.776 ns [671.112, 678.105] (0.95 confidence)</td>
<td>1.02272</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d3/prepared</td>
<td>384.526 ns [383.898, 385.148] (0.95 confidence)</td>
<td>380.465 ns [379.286, 381.787] (0.95 confidence)</td>
<td>1.01067</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d4/adapter</td>
<td>1782.92 ns [1780.53, 1785.35] (0.95 confidence)</td>
<td>1740.61 ns [1737.95, 1743.13] (0.95 confidence)</td>
<td>1.02431</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d4/construct</td>
<td>1052.97 ns [1033.03, 1074.55] (0.95 confidence)</td>
<td>1041.57 ns [1019.9, 1065.22] (0.95 confidence)</td>
<td>1.01094</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d4/prepared</td>
<td>626.164 ns [625.018, 627.334] (0.95 confidence)</td>
<td>624.953 ns [624.111, 625.776] (0.95 confidence)</td>
<td>1.00194</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d5/adapter</td>
<td>2523.92 ns [2518.06, 2529.88] (0.95 confidence)</td>
<td>2465.47 ns [2460.01, 2470.83] (0.95 confidence)</td>
<td>1.02371</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d5/construct</td>
<td>1433.95 ns [1427.54, 1439.98] (0.95 confidence)</td>
<td>1406.69 ns [1397.56, 1414.96] (0.95 confidence)</td>
<td>1.01938</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d5/prepared</td>
<td>869.399 ns [868.393, 870.444] (0.95 confidence)</td>
<td>870.271 ns [869.047, 871.55] (0.95 confidence)</td>
<td>0.998998</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d6/adapter</td>
<td>3368.76 ns [3363.37, 3374.25] (0.95 confidence)</td>
<td>3350.99 ns [3345.18, 3356.83] (0.95 confidence)</td>
<td>1.0053</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d6/construct</td>
<td>1988.85 ns [1980.56, 1996.91] (0.95 confidence)</td>
<td>1978.76 ns [1973.13, 1984.48] (0.95 confidence)</td>
<td>1.0051</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d6/prepared</td>
<td>1123.36 ns [1121.7, 1125.05] (0.95 confidence)</td>
<td>1136.55 ns [1134.45, 1138.79] (0.95 confidence)</td>
<td>0.988394</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d7/adapter</td>
<td>4461 ns [4453.13, 4468.67] (0.95 confidence)</td>
<td>4439.28 ns [4433.28, 4445.42] (0.95 confidence)</td>
<td>1.00489</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d7/construct</td>
<td>2839.2 ns [2779.75, 2904.1] (0.95 confidence)</td>
<td>2848.64 ns [2786.67, 2916.16] (0.95 confidence)</td>
<td>0.996685</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d7/prepared</td>
<td>1469.45 ns [1465.12, 1474.2] (0.95 confidence)</td>
<td>1467.58 ns [1465.96, 1469.3] (0.95 confidence)</td>
<td>1.00128</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d8/adapter</td>
<td>5635.41 ns [5624.07, 5648.55] (0.95 confidence)</td>
<td>5642.4 ns [5632.35, 5652.2] (0.95 confidence)</td>
<td>0.998761</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d8/construct</td>
<td>3501.33 ns [3487.77, 3513.88] (0.95 confidence)</td>
<td>3494.98 ns [3484.58, 3505.06] (0.95 confidence)</td>
<td>1.00182</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d8/prepared</td>
<td>1852.31 ns [1849.95, 1854.64] (0.95 confidence)</td>
<td>1864.93 ns [1862.09, 1867.73] (0.95 confidence)</td>
<td>0.993233</td>
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
<td>168.584 ns [167.986, 168.935] (0.95 confidence)</td>
<td>163.201 ns [162.783, 163.697] (0.95 confidence)</td>
<td>1.03299</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d2/rounded&#95;result</td>
<td>170.641 ns [170.036, 171.365] (0.95 confidence)</td>
<td>164.306 ns [164.003, 164.809] (0.95 confidence)</td>
<td>1.03856</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d2/strict&#95;array</td>
<td>173.556 ns [172.251, 174.555] (0.95 confidence)</td>
<td>167.341 ns [166.916, 167.707] (0.95 confidence)</td>
<td>1.03714</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d2/strict&#95;result</td>
<td>172.543 ns [171.801, 173.542] (0.95 confidence)</td>
<td>169.622 ns [169.016, 170.51] (0.95 confidence)</td>
<td>1.01722</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d3/rounded&#95;array</td>
<td>176.447 ns [175.906, 177.302] (0.95 confidence)</td>
<td>169.096 ns [168.564, 170.524] (0.95 confidence)</td>
<td>1.04347</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d3/rounded&#95;result</td>
<td>176.745 ns [175.857, 177.525] (0.95 confidence)</td>
<td>172.978 ns [172.265, 173.44] (0.95 confidence)</td>
<td>1.02178</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d3/strict&#95;array</td>
<td>176.631 ns [176.176, 177.42] (0.95 confidence)</td>
<td>170.356 ns [169.862, 170.62] (0.95 confidence)</td>
<td>1.03683</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d3/strict&#95;result</td>
<td>175.966 ns [175.643, 176.486] (0.95 confidence)</td>
<td>169.411 ns [168.761, 170.167] (0.95 confidence)</td>
<td>1.03869</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d4/rounded&#95;array</td>
<td>178.127 ns [177.46, 179.745] (0.95 confidence)</td>
<td>174.702 ns [174.096, 175.272] (0.95 confidence)</td>
<td>1.0196</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d4/rounded&#95;result</td>
<td>177.28 ns [176.503, 177.813] (0.95 confidence)</td>
<td>176.74 ns [176.184, 177.318] (0.95 confidence)</td>
<td>1.00305</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d4/strict&#95;array</td>
<td>176.444 ns [175.863, 177.231] (0.95 confidence)</td>
<td>174.133 ns [173.826, 174.564] (0.95 confidence)</td>
<td>1.01327</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d4/strict&#95;result</td>
<td>176.198 ns [175.936, 176.531] (0.95 confidence)</td>
<td>173.425 ns [173.076, 174.222] (0.95 confidence)</td>
<td>1.01599</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d5/rounded&#95;array</td>
<td>185.926 ns [185.225, 186.405] (0.95 confidence)</td>
<td>181.295 ns [180.857, 182.266] (0.95 confidence)</td>
<td>1.02554</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d5/rounded&#95;result</td>
<td>182.523 ns [182.006, 183.658] (0.95 confidence)</td>
<td>181.925 ns [181.394, 182.335] (0.95 confidence)</td>
<td>1.00329</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d5/strict&#95;array</td>
<td>182.53 ns [182.08, 182.765] (0.95 confidence)</td>
<td>181.455 ns [180.738, 182.416] (0.95 confidence)</td>
<td>1.00593</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;below&#95;overflow&#95;d5/strict&#95;result</td>
<td>178.724 ns [178.478, 179.565] (0.95 confidence)</td>
<td>178.416 ns [178.184, 178.788] (0.95 confidence)</td>
<td>1.00172</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d2/rounded&#95;array</td>
<td>11.734 ns [11.7201, 11.7477] (0.95 confidence)</td>
<td>11.5526 ns [11.5469, 11.5552] (0.95 confidence)</td>
<td>1.0157</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d2/rounded&#95;result</td>
<td>11.7381 ns [11.7283, 11.7535] (0.95 confidence)</td>
<td>11.5539 ns [11.5459, 11.5569] (0.95 confidence)</td>
<td>1.01594</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d2/strict&#95;array</td>
<td>8.72057 ns [8.6396, 8.75761] (0.95 confidence)</td>
<td>8.42613 ns [8.41133, 8.44143] (0.95 confidence)</td>
<td>1.03494</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d2/strict&#95;result</td>
<td>8.12464 ns [8.08833, 8.21402] (0.95 confidence)</td>
<td>7.91059 ns [7.89658, 7.99753] (0.95 confidence)</td>
<td>1.02706</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d3/rounded&#95;array</td>
<td>18.0556 ns [18.0342, 18.0762] (0.95 confidence)</td>
<td>17.8625 ns [17.854, 17.8724] (0.95 confidence)</td>
<td>1.01081</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d3/rounded&#95;result</td>
<td>18.0365 ns [18.0285, 18.0493] (0.95 confidence)</td>
<td>17.8649 ns [17.8583, 17.8699] (0.95 confidence)</td>
<td>1.00961</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d3/strict&#95;array</td>
<td>12.7161 ns [12.5754, 12.9025] (0.95 confidence)</td>
<td>12.1249 ns [12.1236, 12.1295] (0.95 confidence)</td>
<td>1.04876</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d3/strict&#95;result</td>
<td>11.365 ns [11.3233, 11.5823] (0.95 confidence)</td>
<td>11.7383 ns [11.6273, 11.9617] (0.95 confidence)</td>
<td>0.968198</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d4/rounded&#95;array</td>
<td>23.7132 ns [23.7091, 23.7185] (0.95 confidence)</td>
<td>23.703 ns [23.6502, 23.7098] (0.95 confidence)</td>
<td>1.00043</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d4/rounded&#95;result</td>
<td>23.7199 ns [23.7166, 23.7276] (0.95 confidence)</td>
<td>23.6553 ns [23.615, 23.6957] (0.95 confidence)</td>
<td>1.00273</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d4/strict&#95;array</td>
<td>16.3397 ns [16.3334, 16.3469] (0.95 confidence)</td>
<td>17.0825 ns [16.6892, 17.2125] (0.95 confidence)</td>
<td>0.956521</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d4/strict&#95;result</td>
<td>15.3238 ns [15.2662, 15.3603] (0.95 confidence)</td>
<td>15.4495 ns [15.3442, 15.6642] (0.95 confidence)</td>
<td>0.991864</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d5/rounded&#95;array</td>
<td>28.4036 ns [28.3985, 28.406] (0.95 confidence)</td>
<td>28.4015 ns [28.3639, 28.4097] (0.95 confidence)</td>
<td>1.00008</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d5/rounded&#95;result</td>
<td>28.4202 ns [28.4124, 28.4313] (0.95 confidence)</td>
<td>28.399 ns [28.3313, 28.4131] (0.95 confidence)</td>
<td>1.00074</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d5/strict&#95;array</td>
<td>20.8894 ns [20.5304, 21.3771] (0.95 confidence)</td>
<td>20.366 ns [20.3082, 20.4094] (0.95 confidence)</td>
<td>1.0257</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;dyadic&#95;d5/strict&#95;result</td>
<td>18.7951 ns [18.6251, 18.8001] (0.95 confidence)</td>
<td>19.7967 ns [19.075, 19.8908] (0.95 confidence)</td>
<td>0.949403</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d2/rounded&#95;array</td>
<td>180.707 ns [180.386, 180.997] (0.95 confidence)</td>
<td>169.005 ns [168.725, 169.199] (0.95 confidence)</td>
<td>1.06924</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d2/rounded&#95;result</td>
<td>180.28 ns [179.845, 180.698] (0.95 confidence)</td>
<td>169.34 ns [168.894, 169.782] (0.95 confidence)</td>
<td>1.0646</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d2/strict&#95;array</td>
<td>14.1245 ns [13.9868, 14.2728] (0.95 confidence)</td>
<td>13.0196 ns [12.7975, 13.3501] (0.95 confidence)</td>
<td>1.08486</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d2/strict&#95;result</td>
<td>12.675 ns [12.6177, 12.7327] (0.95 confidence)</td>
<td>12.218 ns [12.1953, 12.2905] (0.95 confidence)</td>
<td>1.03741</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d3/rounded&#95;array</td>
<td>179.061 ns [178.852, 179.486] (0.95 confidence)</td>
<td>176.188 ns [175.78, 176.571] (0.95 confidence)</td>
<td>1.01631</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d3/rounded&#95;result</td>
<td>180.022 ns [179.668, 180.509] (0.95 confidence)</td>
<td>175.42 ns [175.147, 175.553] (0.95 confidence)</td>
<td>1.02624</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d3/strict&#95;array</td>
<td>19.0595 ns [18.813, 19.2734] (0.95 confidence)</td>
<td>18.1843 ns [18.0743, 18.3131] (0.95 confidence)</td>
<td>1.04813</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d3/strict&#95;result</td>
<td>17.2229 ns [17.1016, 17.4989] (0.95 confidence)</td>
<td>16.5004 ns [16.3633, 16.6843] (0.95 confidence)</td>
<td>1.04379</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d4/rounded&#95;array</td>
<td>190.832 ns [190.631, 191.08] (0.95 confidence)</td>
<td>183.395 ns [182.979, 183.83] (0.95 confidence)</td>
<td>1.04055</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d4/rounded&#95;result</td>
<td>191.272 ns [190.907, 191.542] (0.95 confidence)</td>
<td>181.917 ns [181.685, 182.22] (0.95 confidence)</td>
<td>1.05142</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d4/strict&#95;array</td>
<td>23.7158 ns [23.5034, 23.7879] (0.95 confidence)</td>
<td>22.7465 ns [22.6046, 23.1862] (0.95 confidence)</td>
<td>1.04261</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d4/strict&#95;result</td>
<td>20.7434 ns [20.2895, 21.1002] (0.95 confidence)</td>
<td>21.3119 ns [21.193, 21.4473] (0.95 confidence)</td>
<td>0.973326</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d5/rounded&#95;array</td>
<td>188.107 ns [187.783, 188.493] (0.95 confidence)</td>
<td>187.454 ns [187.038, 187.74] (0.95 confidence)</td>
<td>1.00349</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d5/rounded&#95;result</td>
<td>188.219 ns [188.002, 188.759] (0.95 confidence)</td>
<td>186.564 ns [186.28, 186.799] (0.95 confidence)</td>
<td>1.00887</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d5/strict&#95;array</td>
<td>27.9083 ns [27.777, 28.1609] (0.95 confidence)</td>
<td>26.7108 ns [26.6798, 26.8032] (0.95 confidence)</td>
<td>1.04483</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;min&#95;subnormal&#95;d5/strict&#95;result</td>
<td>24.7 ns [24.5916, 24.9031] (0.95 confidence)</td>
<td>24.2005 ns [23.478, 24.6662] (0.95 confidence)</td>
<td>1.02064</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d2/rounded&#95;array</td>
<td>180.596 ns [180.053, 181.032] (0.95 confidence)</td>
<td>170.044 ns [169.61, 170.721] (0.95 confidence)</td>
<td>1.06205</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d2/rounded&#95;result</td>
<td>180.604 ns [180.122, 180.69] (0.95 confidence)</td>
<td>168.592 ns [168.162, 168.836] (0.95 confidence)</td>
<td>1.07125</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d2/strict&#95;array</td>
<td>13.5199 ns [13.459, 13.7216] (0.95 confidence)</td>
<td>13.0706 ns [12.9416, 13.2096] (0.95 confidence)</td>
<td>1.03437</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d2/strict&#95;result</td>
<td>12.5184 ns [12.4476, 12.6099] (0.95 confidence)</td>
<td>12.1827 ns [12.1235, 12.3912] (0.95 confidence)</td>
<td>1.02756</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d3/rounded&#95;array</td>
<td>180.118 ns [179.483, 180.298] (0.95 confidence)</td>
<td>174.776 ns [174.276, 175.009] (0.95 confidence)</td>
<td>1.03057</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d3/rounded&#95;result</td>
<td>178.911 ns [178.545, 179.349] (0.95 confidence)</td>
<td>176.048 ns [175.982, 176.275] (0.95 confidence)</td>
<td>1.01626</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d3/strict&#95;array</td>
<td>18.0437 ns [17.9839, 18.1504] (0.95 confidence)</td>
<td>17.1959 ns [17.0671, 17.4524] (0.95 confidence)</td>
<td>1.0493</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d3/strict&#95;result</td>
<td>16.7632 ns [16.6552, 16.9432] (0.95 confidence)</td>
<td>16.2167 ns [16.0097, 16.4334] (0.95 confidence)</td>
<td>1.0337</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d4/rounded&#95;array</td>
<td>191.23 ns [190.853, 191.701] (0.95 confidence)</td>
<td>180.842 ns [180.528, 181.223] (0.95 confidence)</td>
<td>1.05744</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d4/rounded&#95;result</td>
<td>191.864 ns [191.536, 192.066] (0.95 confidence)</td>
<td>180.851 ns [180.666, 181.25] (0.95 confidence)</td>
<td>1.06089</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d4/strict&#95;array</td>
<td>21.9628 ns [21.8974, 22.1144] (0.95 confidence)</td>
<td>21.3794 ns [21.3124, 21.5155] (0.95 confidence)</td>
<td>1.02729</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d4/strict&#95;result</td>
<td>21.253 ns [21.0633, 21.4093] (0.95 confidence)</td>
<td>20.7213 ns [20.4145, 20.8287] (0.95 confidence)</td>
<td>1.02566</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d5/rounded&#95;array</td>
<td>188.146 ns [187.919, 188.684] (0.95 confidence)</td>
<td>185.259 ns [184.828, 185.741] (0.95 confidence)</td>
<td>1.01558</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d5/rounded&#95;result</td>
<td>188.261 ns [187.948, 188.534] (0.95 confidence)</td>
<td>185.846 ns [185.496, 186.109] (0.95 confidence)</td>
<td>1.01299</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d5/strict&#95;array</td>
<td>26.2131 ns [25.8897, 26.5021] (0.95 confidence)</td>
<td>24.8071 ns [24.6507, 25.1521] (0.95 confidence)</td>
<td>1.05668</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;negative&#95;underflow&#95;d5/strict&#95;result</td>
<td>24.0074 ns [23.912, 24.0991] (0.95 confidence)</td>
<td>24.319 ns [24.1192, 24.4239] (0.95 confidence)</td>
<td>0.987186</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d2/rounded&#95;array</td>
<td>11.7702 ns [11.7602, 11.7803] (0.95 confidence)</td>
<td>11.5749 ns [11.5675, 11.5847] (0.95 confidence)</td>
<td>1.01687</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d2/rounded&#95;result</td>
<td>11.6969 ns [11.6881, 11.7024] (0.95 confidence)</td>
<td>11.5784 ns [11.5704, 11.5859] (0.95 confidence)</td>
<td>1.01023</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d2/strict&#95;array</td>
<td>56.4589 ns [56.3511, 56.5754] (0.95 confidence)</td>
<td>54.7786 ns [54.6539, 54.8806] (0.95 confidence)</td>
<td>1.03067</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d2/strict&#95;result</td>
<td>11.91 ns [11.8628, 12.0215] (0.95 confidence)</td>
<td>11.8605 ns [11.7817, 11.9074] (0.95 confidence)</td>
<td>1.00417</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d3/rounded&#95;array</td>
<td>18.1032 ns [18.0895, 18.1177] (0.95 confidence)</td>
<td>17.9642 ns [17.9564, 17.9705] (0.95 confidence)</td>
<td>1.00774</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d3/rounded&#95;result</td>
<td>18.1233 ns [18.1169, 18.1387] (0.95 confidence)</td>
<td>17.9463 ns [17.9253, 17.9549] (0.95 confidence)</td>
<td>1.00986</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d3/strict&#95;array</td>
<td>60.4605 ns [60.3071, 60.6253] (0.95 confidence)</td>
<td>59.1318 ns [58.985, 59.2249] (0.95 confidence)</td>
<td>1.02247</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d3/strict&#95;result</td>
<td>15.6172 ns [15.5972, 15.6567] (0.95 confidence)</td>
<td>15.5007 ns [15.4762, 15.5305] (0.95 confidence)</td>
<td>1.00752</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d4/rounded&#95;array</td>
<td>23.8113 ns [23.807, 23.818] (0.95 confidence)</td>
<td>23.7783 ns [23.7395, 23.8083] (0.95 confidence)</td>
<td>1.00139</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d4/rounded&#95;result</td>
<td>23.8192 ns [23.8117, 23.8321] (0.95 confidence)</td>
<td>23.7968 ns [23.742, 23.8312] (0.95 confidence)</td>
<td>1.00094</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d4/strict&#95;array</td>
<td>63.0875 ns [62.8624, 63.3715] (0.95 confidence)</td>
<td>62.4699 ns [62.1941, 62.8129] (0.95 confidence)</td>
<td>1.00989</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d4/strict&#95;result</td>
<td>19.4672 ns [19.4331, 19.5154] (0.95 confidence)</td>
<td>19.442 ns [19.4031, 19.5037] (0.95 confidence)</td>
<td>1.0013</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d5/rounded&#95;array</td>
<td>28.3348 ns [28.3279, 28.3395] (0.95 confidence)</td>
<td>28.3403 ns [28.3316, 28.3481] (0.95 confidence)</td>
<td>0.999803</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d5/rounded&#95;result</td>
<td>28.3542 ns [28.3443, 28.3629] (0.95 confidence)</td>
<td>28.3388 ns [28.3077, 28.346] (0.95 confidence)</td>
<td>1.00054</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d5/strict&#95;array</td>
<td>65.9816 ns [65.8492, 66.2297] (0.95 confidence)</td>
<td>64.8969 ns [64.7726, 65.1455] (0.95 confidence)</td>
<td>1.01672</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;non&#95;dyadic&#95;d5/strict&#95;result</td>
<td>22.4043 ns [22.3402, 22.4377] (0.95 confidence)</td>
<td>22.6312 ns [22.5075, 22.9263] (0.95 confidence)</td>
<td>0.989975</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d2/rounded&#95;array</td>
<td>170.854 ns [170.398, 171.71] (0.95 confidence)</td>
<td>165.544 ns [165.172, 165.923] (0.95 confidence)</td>
<td>1.03208</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d2/rounded&#95;result</td>
<td>171.42 ns [170.315, 172.177] (0.95 confidence)</td>
<td>167.625 ns [167.246, 167.842] (0.95 confidence)</td>
<td>1.02264</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d2/strict&#95;array</td>
<td>179.601 ns [179.038, 180.371] (0.95 confidence)</td>
<td>173.121 ns [172.841, 173.677] (0.95 confidence)</td>
<td>1.03743</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d2/strict&#95;result</td>
<td>180.315 ns [179.828, 180.882] (0.95 confidence)</td>
<td>172.297 ns [171.22, 172.936] (0.95 confidence)</td>
<td>1.04654</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d3/rounded&#95;array</td>
<td>177.212 ns [176.807, 177.797] (0.95 confidence)</td>
<td>172.07 ns [171.208, 172.559] (0.95 confidence)</td>
<td>1.02988</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d3/rounded&#95;result</td>
<td>177.677 ns [177.206, 178.233] (0.95 confidence)</td>
<td>171.16 ns [170.765, 171.486] (0.95 confidence)</td>
<td>1.03807</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d3/strict&#95;array</td>
<td>184.678 ns [184.199, 185.4] (0.95 confidence)</td>
<td>178.889 ns [177.705, 179.234] (0.95 confidence)</td>
<td>1.03236</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d3/strict&#95;result</td>
<td>181.53 ns [180.892, 182.355] (0.95 confidence)</td>
<td>176.659 ns [176.152, 176.992] (0.95 confidence)</td>
<td>1.02757</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d4/rounded&#95;array</td>
<td>181.75 ns [181.192, 182.166] (0.95 confidence)</td>
<td>178.247 ns [177.724, 178.81] (0.95 confidence)</td>
<td>1.01965</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d4/rounded&#95;result</td>
<td>180.585 ns [179.526, 181.387] (0.95 confidence)</td>
<td>178.403 ns [177.776, 179.207] (0.95 confidence)</td>
<td>1.01223</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d4/strict&#95;array</td>
<td>183.987 ns [183.36, 184.867] (0.95 confidence)</td>
<td>181.645 ns [181.234, 182.175] (0.95 confidence)</td>
<td>1.01289</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d4/strict&#95;result</td>
<td>181.538 ns [180.756, 182.633] (0.95 confidence)</td>
<td>179.213 ns [178.775, 179.812] (0.95 confidence)</td>
<td>1.01297</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d5/rounded&#95;array</td>
<td>187.755 ns [186.899, 188.668] (0.95 confidence)</td>
<td>183.677 ns [183.154, 183.975] (0.95 confidence)</td>
<td>1.0222</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d5/rounded&#95;result</td>
<td>188.055 ns [187.497, 188.53] (0.95 confidence)</td>
<td>183.097 ns [182.518, 183.888] (0.95 confidence)</td>
<td>1.02708</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d5/strict&#95;array</td>
<td>190.611 ns [189.735, 191.87] (0.95 confidence)</td>
<td>183.803 ns [183.194, 184.402] (0.95 confidence)</td>
<td>1.03704</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;overflow&#95;midpoint&#95;d5/strict&#95;result</td>
<td>184.19 ns [184.054, 184.873] (0.95 confidence)</td>
<td>182.524 ns [182.342, 182.807] (0.95 confidence)</td>
<td>1.00913</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d2/rounded&#95;array</td>
<td>204.02 ns [203.413, 204.807] (0.95 confidence)</td>
<td>197.936 ns [196.92, 198.672] (0.95 confidence)</td>
<td>1.03074</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d2/rounded&#95;result</td>
<td>204.48 ns [204.026, 204.855] (0.95 confidence)</td>
<td>195.819 ns [195.67, 196.052] (0.95 confidence)</td>
<td>1.04422</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d2/strict&#95;array</td>
<td>13506.1 ns [13469.1, 13545.2] (0.95 confidence)</td>
<td>13172.2 ns [13154.5, 13194.4] (0.95 confidence)</td>
<td>1.02535</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d2/strict&#95;result</td>
<td>204.977 ns [204.162, 205.911] (0.95 confidence)</td>
<td>197.123 ns [196.88, 198.01] (0.95 confidence)</td>
<td>1.03984</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d3/rounded&#95;array</td>
<td>209.512 ns [209.072, 209.963] (0.95 confidence)</td>
<td>201.916 ns [201.509, 202.291] (0.95 confidence)</td>
<td>1.03762</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d3/rounded&#95;result</td>
<td>210.376 ns [209.844, 210.896] (0.95 confidence)</td>
<td>205.225 ns [205.022, 205.63] (0.95 confidence)</td>
<td>1.0251</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d3/strict&#95;array</td>
<td>13506.1 ns [13467.3, 13525.3] (0.95 confidence)</td>
<td>13317.9 ns [13278.6, 13347.4] (0.95 confidence)</td>
<td>1.01413</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d3/strict&#95;result</td>
<td>209.244 ns [208.684, 209.755] (0.95 confidence)</td>
<td>201.956 ns [201.437, 202.675] (0.95 confidence)</td>
<td>1.03608</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d4/rounded&#95;array</td>
<td>209.873 ns [209.622, 210.193] (0.95 confidence)</td>
<td>210.271 ns [210.031, 210.445] (0.95 confidence)</td>
<td>0.99811</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d4/rounded&#95;result</td>
<td>211.321 ns [210.661, 211.856] (0.95 confidence)</td>
<td>211.049 ns [210.516, 211.558] (0.95 confidence)</td>
<td>1.00128</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d4/strict&#95;array</td>
<td>13360.3 ns [13309.3, 13403.3] (0.95 confidence)</td>
<td>13275.5 ns [13252.9, 13300.2] (0.95 confidence)</td>
<td>1.00638</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d4/strict&#95;result</td>
<td>205.904 ns [205.089, 206.358] (0.95 confidence)</td>
<td>211.488 ns [210.719, 212.108] (0.95 confidence)</td>
<td>0.973597</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d5/rounded&#95;array</td>
<td>216.723 ns [216.162, 217.294] (0.95 confidence)</td>
<td>214.586 ns [213.424, 215.382] (0.95 confidence)</td>
<td>1.00996</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d5/rounded&#95;result</td>
<td>216.002 ns [215.802, 216.328] (0.95 confidence)</td>
<td>213.659 ns [213.053, 214.3] (0.95 confidence)</td>
<td>1.01097</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d5/strict&#95;array</td>
<td>13283.5 ns [13254, 13305.4] (0.95 confidence)</td>
<td>13349.2 ns [13314.7, 13386] (0.95 confidence)</td>
<td>0.995081</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide1024&#95;d5/strict&#95;result</td>
<td>213.518 ns [211.521, 215.951] (0.95 confidence)</td>
<td>210.849 ns [210.12, 211.94] (0.95 confidence)</td>
<td>1.01266</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d2/rounded&#95;array</td>
<td>165.53 ns [164.886, 166.01] (0.95 confidence)</td>
<td>158.774 ns [158.445, 159.061] (0.95 confidence)</td>
<td>1.04256</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d2/rounded&#95;result</td>
<td>165.126 ns [164.584, 166.144] (0.95 confidence)</td>
<td>156.464 ns [156.051, 156.862] (0.95 confidence)</td>
<td>1.05536</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d2/strict&#95;array</td>
<td>4012.81 ns [3981.48, 4032.39] (0.95 confidence)</td>
<td>3877.29 ns [3864.49, 3898.69] (0.95 confidence)</td>
<td>1.03495</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d2/strict&#95;result</td>
<td>165.736 ns [165.391, 166.225] (0.95 confidence)</td>
<td>160.634 ns [159.583, 161.653] (0.95 confidence)</td>
<td>1.03176</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d3/rounded&#95;array</td>
<td>172.023 ns [171.022, 173.031] (0.95 confidence)</td>
<td>168.504 ns [167.434, 169.655] (0.95 confidence)</td>
<td>1.02088</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d3/rounded&#95;result</td>
<td>172.361 ns [171.462, 172.653] (0.95 confidence)</td>
<td>165.377 ns [164.918, 167.727] (0.95 confidence)</td>
<td>1.04223</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d3/strict&#95;array</td>
<td>3986.69 ns [3947.77, 4004.07] (0.95 confidence)</td>
<td>3880.79 ns [3870.62, 3908.87] (0.95 confidence)</td>
<td>1.02729</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d3/strict&#95;result</td>
<td>172.307 ns [171.715, 172.884] (0.95 confidence)</td>
<td>169.462 ns [169.033, 170.094] (0.95 confidence)</td>
<td>1.01679</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d4/rounded&#95;array</td>
<td>170.744 ns [169.313, 171.677] (0.95 confidence)</td>
<td>171.455 ns [170.974, 172.181] (0.95 confidence)</td>
<td>0.995851</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d4/rounded&#95;result</td>
<td>174.783 ns [174.28, 176.737] (0.95 confidence)</td>
<td>169.992 ns [169.41, 170.268] (0.95 confidence)</td>
<td>1.02819</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d4/strict&#95;array</td>
<td>3963 ns [3940.22, 3979.98] (0.95 confidence)</td>
<td>3898.74 ns [3882.99, 3905.57] (0.95 confidence)</td>
<td>1.01648</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d4/strict&#95;result</td>
<td>174.585 ns [172.707, 175.851] (0.95 confidence)</td>
<td>170.622 ns [170.159, 171.072] (0.95 confidence)</td>
<td>1.02322</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d5/rounded&#95;array</td>
<td>181.119 ns [180.569, 181.624] (0.95 confidence)</td>
<td>178.282 ns [177.614, 179.033] (0.95 confidence)</td>
<td>1.01591</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d5/rounded&#95;result</td>
<td>178.021 ns [177.453, 178.791] (0.95 confidence)</td>
<td>175.699 ns [174.597, 176.694] (0.95 confidence)</td>
<td>1.01321</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d5/strict&#95;array</td>
<td>3933.12 ns [3912.44, 3943.15] (0.95 confidence)</td>
<td>3888.86 ns [3876.25, 3905.75] (0.95 confidence)</td>
<td>1.01138</td>
</tr>
<tr>
<td>canonical&#95;conversion&#95;wide256&#95;d5/strict&#95;result</td>
<td>177.025 ns [175.712, 177.737] (0.95 confidence)</td>
<td>174.642 ns [174.097, 175.224] (0.95 confidence)</td>
<td>1.01365</td>
</tr>
<tr>
<td>exact&#95;d2/det</td>
<td>0.440662 ns [0.438562, 0.444248] (0.95 confidence)</td>
<td>0.436954 ns [0.436421, 0.437256] (0.95 confidence)</td>
<td>1.00848</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;direct</td>
<td>0.437997 ns [0.436532, 0.441574] (0.95 confidence)</td>
<td>0.436919 ns [0.436685, 0.437284] (0.95 confidence)</td>
<td>1.00247</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;direct&#95;with&#95;errbound</td>
<td>1.74747 ns [1.74551, 1.75157] (0.95 confidence)</td>
<td>1.72991 ns [1.72908, 1.73156] (0.95 confidence)</td>
<td>1.01015</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;errbound</td>
<td>1.59761 ns [1.59633, 1.59878] (0.95 confidence)</td>
<td>1.58313 ns [1.58303, 1.58316] (0.95 confidence)</td>
<td>1.00915</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;exact</td>
<td>90.169 ns [89.9706, 90.3425] (0.95 confidence)</td>
<td>89.2496 ns [89.1382, 89.5308] (0.95 confidence)</td>
<td>1.0103</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;exact&#95;f64&#95;result</td>
<td>70.699 ns [70.5412, 70.9039] (0.95 confidence)</td>
<td>71.4265 ns [71.3262, 71.5602] (0.95 confidence)</td>
<td>0.989815</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;exact&#95;rounded&#95;f64</td>
<td>72.0505 ns [71.9832, 72.106] (0.95 confidence)</td>
<td>73.0551 ns [72.8223, 73.2335] (0.95 confidence)</td>
<td>0.986249</td>
</tr>
<tr>
<td>exact&#95;d2/det&#95;sign&#95;exact</td>
<td>2.99501 ns [2.99267, 2.99722] (0.95 confidence)</td>
<td>2.95251 ns [2.94729, 2.95539] (0.95 confidence)</td>
<td>1.0144</td>
</tr>
<tr>
<td>exact&#95;d2/solve&#95;exact</td>
<td>3010.02 ns [3006.93, 3015.64] (0.95 confidence)</td>
<td>2960.92 ns [2957.12, 2968.02] (0.95 confidence)</td>
<td>1.01658</td>
</tr>
<tr>
<td>exact&#95;d2/solve&#95;exact&#95;f64&#95;result</td>
<td>3171.21 ns [3165.85, 3177.8] (0.95 confidence)</td>
<td>3110.91 ns [3106.68, 3115.49] (0.95 confidence)</td>
<td>1.01938</td>
</tr>
<tr>
<td>exact&#95;d2/solve&#95;exact&#95;rounded&#95;f64</td>
<td>3322.2 ns [3313.33, 3327.79] (0.95 confidence)</td>
<td>3262.18 ns [3254.21, 3265.3] (0.95 confidence)</td>
<td>1.0184</td>
</tr>
<tr>
<td>exact&#95;d3/det</td>
<td>0.82615 ns [0.82558, 0.826793] (0.95 confidence)</td>
<td>0.81827 ns [0.818213, 0.818371] (0.95 confidence)</td>
<td>1.00963</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;direct</td>
<td>0.832009 ns [0.830986, 0.834565] (0.95 confidence)</td>
<td>0.818215 ns [0.818096, 0.818266] (0.95 confidence)</td>
<td>1.01686</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;direct&#95;with&#95;errbound</td>
<td>3.3062 ns [3.30342, 3.30984] (0.95 confidence)</td>
<td>3.27244 ns [3.27001, 3.27385] (0.95 confidence)</td>
<td>1.01032</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;errbound</td>
<td>3.25603 ns [3.25202, 3.26141] (0.95 confidence)</td>
<td>3.21745 ns [3.21655, 3.21902] (0.95 confidence)</td>
<td>1.01199</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;exact</td>
<td>309.316 ns [308.891, 309.877] (0.95 confidence)</td>
<td>302.669 ns [302.118, 303.082] (0.95 confidence)</td>
<td>1.02196</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;exact&#95;f64&#95;result</td>
<td>283.216 ns [282.73, 283.797] (0.95 confidence)</td>
<td>276.433 ns [276.141, 276.922] (0.95 confidence)</td>
<td>1.02454</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;exact&#95;rounded&#95;f64</td>
<td>285.151 ns [284.765, 285.586] (0.95 confidence)</td>
<td>280.907 ns [280.051, 281.377] (0.95 confidence)</td>
<td>1.01511</td>
</tr>
<tr>
<td>exact&#95;d3/det&#95;sign&#95;exact</td>
<td>4.75301 ns [4.74681, 4.75956] (0.95 confidence)</td>
<td>4.69498 ns [4.68518, 4.7001] (0.95 confidence)</td>
<td>1.01236</td>
</tr>
<tr>
<td>exact&#95;d3/solve&#95;exact</td>
<td>9211.56 ns [9188.81, 9230.72] (0.95 confidence)</td>
<td>9236.7 ns [9204.64, 9268.71] (0.95 confidence)</td>
<td>0.997278</td>
</tr>
<tr>
<td>exact&#95;d3/solve&#95;exact&#95;f64&#95;result</td>
<td>9419.52 ns [9392.66, 9444.5] (0.95 confidence)</td>
<td>9333.06 ns [9291.78, 9361.81] (0.95 confidence)</td>
<td>1.00926</td>
</tr>
<tr>
<td>exact&#95;d3/solve&#95;exact&#95;rounded&#95;f64</td>
<td>9707.99 ns [9695.37, 9731.39] (0.95 confidence)</td>
<td>9868.17 ns [9824.2, 9919.58] (0.95 confidence)</td>
<td>0.983768</td>
</tr>
<tr>
<td>exact&#95;d4/det</td>
<td>4.78303 ns [4.78044, 4.78677] (0.95 confidence)</td>
<td>4.73798 ns [4.73625, 4.73961] (0.95 confidence)</td>
<td>1.00951</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;direct</td>
<td>2.37025 ns [2.36901, 2.37276] (0.95 confidence)</td>
<td>2.34691 ns [2.34578, 2.34752] (0.95 confidence)</td>
<td>1.00994</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;direct&#95;with&#95;errbound</td>
<td>6.82586 ns [6.82085, 6.83212] (0.95 confidence)</td>
<td>6.76795 ns [6.76441, 6.77134] (0.95 confidence)</td>
<td>1.00856</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;errbound</td>
<td>6.80414 ns [6.79529, 6.81298] (0.95 confidence)</td>
<td>6.71246 ns [6.69935, 6.72632] (0.95 confidence)</td>
<td>1.01366</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;exact</td>
<td>858.4 ns [857.272, 860.716] (0.95 confidence)</td>
<td>853.771 ns [852.288, 855.044] (0.95 confidence)</td>
<td>1.00542</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;exact&#95;f64&#95;result</td>
<td>819.452 ns [818.3, 821.163] (0.95 confidence)</td>
<td>813.359 ns [810.666, 816.134] (0.95 confidence)</td>
<td>1.00749</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;exact&#95;rounded&#95;f64</td>
<td>839.691 ns [836.798, 841.939] (0.95 confidence)</td>
<td>830.766 ns [828.065, 832.591] (0.95 confidence)</td>
<td>1.01074</td>
</tr>
<tr>
<td>exact&#95;d4/det&#95;sign&#95;exact</td>
<td>7.97623 ns [7.97239, 7.98179] (0.95 confidence)</td>
<td>7.90821 ns [7.89855, 7.90983] (0.95 confidence)</td>
<td>1.0086</td>
</tr>
<tr>
<td>exact&#95;d4/solve&#95;exact</td>
<td>18994.2 ns [18955, 19021.9] (0.95 confidence)</td>
<td>19367.8 ns [19293.3, 19490.9] (0.95 confidence)</td>
<td>0.980707</td>
</tr>
<tr>
<td>exact&#95;d4/solve&#95;exact&#95;f64&#95;result</td>
<td>19232.4 ns [19180.3, 19252.7] (0.95 confidence)</td>
<td>19419.8 ns [19300.2, 19556.9] (0.95 confidence)</td>
<td>0.990349</td>
</tr>
<tr>
<td>exact&#95;d4/solve&#95;exact&#95;rounded&#95;f64</td>
<td>19645.6 ns [19622.7, 19672.2] (0.95 confidence)</td>
<td>19554.5 ns [19523.7, 19624.5] (0.95 confidence)</td>
<td>1.00466</td>
</tr>
<tr>
<td>exact&#95;d5/det</td>
<td>32.2692 ns [32.2081, 32.3201] (0.95 confidence)</td>
<td>31.8374 ns [31.8061, 31.8662] (0.95 confidence)</td>
<td>1.01356</td>
</tr>
<tr>
<td>exact&#95;d5/det&#95;exact</td>
<td>2922.24 ns [2918.04, 2926.89] (0.95 confidence)</td>
<td>2858.25 ns [2848.4, 2871.8] (0.95 confidence)</td>
<td>1.02239</td>
</tr>
<tr>
<td>exact&#95;d5/det&#95;exact&#95;f64&#95;result</td>
<td>2870.84 ns [2866.79, 2876.14] (0.95 confidence)</td>
<td>2823.65 ns [2817.46, 2839.22] (0.95 confidence)</td>
<td>1.01671</td>
</tr>
<tr>
<td>exact&#95;d5/det&#95;exact&#95;rounded&#95;f64</td>
<td>2899.38 ns [2894.89, 2904] (0.95 confidence)</td>
<td>2802.17 ns [2798.46, 2807.32] (0.95 confidence)</td>
<td>1.03469</td>
</tr>
<tr>
<td>exact&#95;d5/det&#95;sign&#95;exact</td>
<td>2872.43 ns [2864.45, 2878.13] (0.95 confidence)</td>
<td>2813.82 ns [2807.68, 2818.47] (0.95 confidence)</td>
<td>1.02083</td>
</tr>
<tr>
<td>exact&#95;d5/solve&#95;exact</td>
<td>31285.1 ns [31191.4, 31358.6] (0.95 confidence)</td>
<td>31041.9 ns [30769.6, 31314.3] (0.95 confidence)</td>
<td>1.00784</td>
</tr>
<tr>
<td>exact&#95;d5/solve&#95;exact&#95;f64&#95;result</td>
<td>31428.7 ns [31307.1, 31489.8] (0.95 confidence)</td>
<td>30804.4 ns [30738.3, 30892.4] (0.95 confidence)</td>
<td>1.02027</td>
</tr>
<tr>
<td>exact&#95;d5/solve&#95;exact&#95;rounded&#95;f64</td>
<td>31428.1 ns [31354.3, 31492] (0.95 confidence)</td>
<td>31385.8 ns [31303, 31469.7] (0.95 confidence)</td>
<td>1.00135</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;4x4/det&#95;exact</td>
<td>929.218 ns [928.616, 930.648] (0.95 confidence)</td>
<td>930.253 ns [928.761, 931.265] (0.95 confidence)</td>
<td>0.998887</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;4x4/det&#95;sign&#95;exact</td>
<td>7.91756 ns [7.91272, 7.93077] (0.95 confidence)</td>
<td>7.90808 ns [7.89341, 7.91101] (0.95 confidence)</td>
<td>1.0012</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;4x4/solve&#95;exact</td>
<td>15730.4 ns [15703.3, 15765.2] (0.95 confidence)</td>
<td>16030.6 ns [15954.1, 16070.3] (0.95 confidence)</td>
<td>0.98127</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;4x4/solve&#95;exact&#95;f64&#95;result</td>
<td>15962.7 ns [15921.2, 15984.4] (0.95 confidence)</td>
<td>16136.6 ns [16106.9, 16172.8] (0.95 confidence)</td>
<td>0.989223</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;4x4/solve&#95;exact&#95;rounded&#95;f64</td>
<td>16180.8 ns [16129.6, 16225] (0.95 confidence)</td>
<td>16449 ns [16398.4, 16517.9] (0.95 confidence)</td>
<td>0.983694</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;5x5/det&#95;exact</td>
<td>2946.77 ns [2941.44, 2951.81] (0.95 confidence)</td>
<td>2903.27 ns [2892.52, 2911.69] (0.95 confidence)</td>
<td>1.01498</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;5x5/det&#95;sign&#95;exact</td>
<td>2895.42 ns [2889.32, 2901.07] (0.95 confidence)</td>
<td>2839.11 ns [2833.39, 2844.35] (0.95 confidence)</td>
<td>1.01983</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;5x5/solve&#95;exact</td>
<td>25240.4 ns [25200.5, 25291.8] (0.95 confidence)</td>
<td>25922.1 ns [25866.9, 25956.9] (0.95 confidence)</td>
<td>0.973702</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;5x5/solve&#95;exact&#95;f64&#95;result</td>
<td>25518.6 ns [25478.4, 25573.1] (0.95 confidence)</td>
<td>26092.6 ns [26003, 26136.5] (0.95 confidence)</td>
<td>0.978001</td>
</tr>
<tr>
<td>exact&#95;hilbert&#95;5x5/solve&#95;exact&#95;rounded&#95;f64</td>
<td>26011.1 ns [25937.3, 26042.2] (0.95 confidence)</td>
<td>26561.2 ns [26473.6, 26636.6] (0.95 confidence)</td>
<td>0.979287</td>
</tr>
<tr>
<td>exact&#95;large&#95;entries&#95;3x3/det&#95;exact</td>
<td>321.463 ns [320.618, 322.388] (0.95 confidence)</td>
<td>319.372 ns [318.457, 320.563] (0.95 confidence)</td>
<td>1.00655</td>
</tr>
<tr>
<td>exact&#95;large&#95;entries&#95;3x3/det&#95;sign&#95;exact</td>
<td>288.815 ns [287.749, 289.799] (0.95 confidence)</td>
<td>290.35 ns [289.838, 290.835] (0.95 confidence)</td>
<td>0.994714</td>
</tr>
<tr>
<td>exact&#95;large&#95;entries&#95;3x3/solve&#95;exact</td>
<td>8993.64 ns [8980.4, 9007.14] (0.95 confidence)</td>
<td>8835.63 ns [8824.05, 8842.12] (0.95 confidence)</td>
<td>1.01788</td>
</tr>
<tr>
<td>exact&#95;large&#95;entries&#95;3x3/solve&#95;exact&#95;f64&#95;result</td>
<td>9168.37 ns [9147.68, 9177.71] (0.95 confidence)</td>
<td>8968.73 ns [8951.71, 8976.67] (0.95 confidence)</td>
<td>1.02226</td>
</tr>
<tr>
<td>exact&#95;large&#95;entries&#95;3x3/solve&#95;exact&#95;rounded&#95;f64</td>
<td>9447.64 ns [9436.42, 9467] (0.95 confidence)</td>
<td>9319.76 ns [9312.19, 9330.77] (0.95 confidence)</td>
<td>1.01372</td>
</tr>
<tr>
<td>exact&#95;near&#95;singular&#95;3x3/det&#95;exact</td>
<td>189.537 ns [189.203, 190.001] (0.95 confidence)</td>
<td>187.904 ns [187.697, 188.014] (0.95 confidence)</td>
<td>1.00869</td>
</tr>
<tr>
<td>exact&#95;near&#95;singular&#95;3x3/det&#95;sign&#95;exact</td>
<td>164.431 ns [164.183, 164.75] (0.95 confidence)</td>
<td>164.677 ns [164.477, 164.795] (0.95 confidence)</td>
<td>0.99851</td>
</tr>
<tr>
<td>exact&#95;near&#95;singular&#95;3x3/solve&#95;exact</td>
<td>987.682 ns [984.222, 994.627] (0.95 confidence)</td>
<td>959.574 ns [958.609, 961.643] (0.95 confidence)</td>
<td>1.02929</td>
</tr>
<tr>
<td>exact&#95;near&#95;singular&#95;3x3/solve&#95;exact&#95;f64&#95;result</td>
<td>1008.5 ns [1006.67, 1010.55] (0.95 confidence)</td>
<td>986.406 ns [983.508, 988.198] (0.95 confidence)</td>
<td>1.02239</td>
</tr>
<tr>
<td>exact&#95;near&#95;singular&#95;3x3/solve&#95;exact&#95;rounded&#95;f64</td>
<td>992.355 ns [990.863, 994.286] (0.95 confidence)</td>
<td>981.003 ns [978.681, 982.183] (0.95 confidence)</td>
<td>1.01157</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d2/det&#95;exact</td>
<td>3523.29 ns [3511.9, 3535.46] (0.95 confidence)</td>
<td>3440.46 ns [3432.47, 3457.31] (0.95 confidence)</td>
<td>1.02407</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d2/det&#95;sign&#95;exact</td>
<td>147.976 ns [147.627, 148.459] (0.95 confidence)</td>
<td>142.615 ns [142.319, 143.124] (0.95 confidence)</td>
<td>1.03759</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d2/solve&#95;exact</td>
<td>36818 ns [36705, 36922.4] (0.95 confidence)</td>
<td>35714.3 ns [35664.8, 35745.5] (0.95 confidence)</td>
<td>1.0309</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d2/solve&#95;exact&#95;f64&#95;result</td>
<td>37300.9 ns [37204.2, 37479.6] (0.95 confidence)</td>
<td>36304.6 ns [36180.7, 36372.6] (0.95 confidence)</td>
<td>1.02744</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d2/solve&#95;exact&#95;rounded&#95;f64</td>
<td>37579.6 ns [37394.7, 37715.9] (0.95 confidence)</td>
<td>36682 ns [36561.2, 36750] (0.95 confidence)</td>
<td>1.02447</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d3/det&#95;exact</td>
<td>9553.25 ns [9535.43, 9569.91] (0.95 confidence)</td>
<td>9280.66 ns [9274.22, 9300.16] (0.95 confidence)</td>
<td>1.02937</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d3/det&#95;sign&#95;exact</td>
<td>232.79 ns [232.074, 233.617] (0.95 confidence)</td>
<td>227.712 ns [227.132, 228.108] (0.95 confidence)</td>
<td>1.0223</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d3/solve&#95;exact</td>
<td>83000 ns [82807, 83082.5] (0.95 confidence)</td>
<td>79946.6 ns [79744.1, 80044.6] (0.95 confidence)</td>
<td>1.03819</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d3/solve&#95;exact&#95;f64&#95;result</td>
<td>83275.7 ns [83142.7, 83404.6] (0.95 confidence)</td>
<td>80353.6 ns [80250.5, 80493.9] (0.95 confidence)</td>
<td>1.03637</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d3/solve&#95;exact&#95;rounded&#95;f64</td>
<td>84210.8 ns [83976.4, 84377.9] (0.95 confidence)</td>
<td>80870.7 ns [80734, 80988] (0.95 confidence)</td>
<td>1.0413</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d4/det&#95;exact</td>
<td>23184.4 ns [23111.1, 23249.6] (0.95 confidence)</td>
<td>22305.8 ns [22250.4, 22350.9] (0.95 confidence)</td>
<td>1.03939</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d4/det&#95;sign&#95;exact</td>
<td>446.48 ns [446.109, 447.212] (0.95 confidence)</td>
<td>438.767 ns [438.272, 439.485] (0.95 confidence)</td>
<td>1.01758</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d4/solve&#95;exact</td>
<td>149657 ns [149335, 150007] (0.95 confidence)</td>
<td>143370 ns [143123, 143609] (0.95 confidence)</td>
<td>1.04385</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d4/solve&#95;exact&#95;f64&#95;result</td>
<td>149882 ns [149607, 150179] (0.95 confidence)</td>
<td>144201 ns [143991, 144432] (0.95 confidence)</td>
<td>1.0394</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d4/solve&#95;exact&#95;rounded&#95;f64</td>
<td>154828 ns [154272, 155505] (0.95 confidence)</td>
<td>145579 ns [145320, 145733] (0.95 confidence)</td>
<td>1.06353</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d5/det&#95;exact</td>
<td>57169.2 ns [56819.7, 57412.4] (0.95 confidence)</td>
<td>55349.1 ns [55190.5, 55745.7] (0.95 confidence)</td>
<td>1.03288</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d5/det&#95;sign&#95;exact</td>
<td>56141.7 ns [55987.9, 56313.8] (0.95 confidence)</td>
<td>54243.8 ns [54013.2, 54480.7] (0.95 confidence)</td>
<td>1.03499</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d5/solve&#95;exact</td>
<td>248282 ns [247451, 249019] (0.95 confidence)</td>
<td>236687 ns [236067, 237384] (0.95 confidence)</td>
<td>1.04899</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d5/solve&#95;exact&#95;f64&#95;result</td>
<td>248915 ns [248302, 249648] (0.95 confidence)</td>
<td>236917 ns [236672, 237443] (0.95 confidence)</td>
<td>1.05064</td>
</tr>
<tr>
<td>exact&#95;random&#95;corpus&#95;d5/solve&#95;exact&#95;rounded&#95;f64</td>
<td>249040 ns [248602, 249827] (0.95 confidence)</td>
<td>238511 ns [237924, 238744] (0.95 confidence)</td>
<td>1.04414</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/matrix&#95;access</td>
<td>0.752881 ns [0.751182, 0.755628] (0.95 confidence)</td>
<td>0.749177 ns [0.747176, 0.750297] (0.95 confidence)</td>
<td>1.00494</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/matrix&#95;from&#95;fn</td>
<td>131.292 ns [131.053, 131.399] (0.95 confidence)</td>
<td>130.506 ns [130.203, 130.713] (0.95 confidence)</td>
<td>1.00602</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/matrix&#95;from&#95;rows</td>
<td>133.929 ns [133.379, 134.196] (0.95 confidence)</td>
<td>133.219 ns [132.736, 133.838] (0.95 confidence)</td>
<td>1.00533</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/matrix&#95;into&#95;rows</td>
<td>64.0189 ns [63.9085, 64.1671] (0.95 confidence)</td>
<td>63.6174 ns [63.4266, 63.7348] (0.95 confidence)</td>
<td>1.00631</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/matrix&#95;zero</td>
<td>2.6892 ns [2.68522, 2.69305] (0.95 confidence)</td>
<td>2.68124 ns [2.67933, 2.68315] (0.95 confidence)</td>
<td>1.00297</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/vector&#95;access</td>
<td>0.757049 ns [0.755139, 0.758964] (0.95 confidence)</td>
<td>0.750505 ns [0.748711, 0.751659] (0.95 confidence)</td>
<td>1.00872</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/vector&#95;from&#95;fn</td>
<td>36.6485 ns [36.5672, 36.6934] (0.95 confidence)</td>
<td>36.4306 ns [36.358, 36.4874] (0.95 confidence)</td>
<td>1.00598</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/vector&#95;into&#95;array</td>
<td>15.3314 ns [15.2971, 15.3683] (0.95 confidence)</td>
<td>15.2669 ns [15.2165, 15.3185] (0.95 confidence)</td>
<td>1.00423</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/vector&#95;new</td>
<td>35.9345 ns [35.8423, 36.0083] (0.95 confidence)</td>
<td>35.8058 ns [35.7094, 35.937] (0.95 confidence)</td>
<td>1.0036</td>
</tr>
<tr>
<td>public&#95;exact&#95;d2/vector&#95;zero</td>
<td>1.35364 ns [1.35201, 1.35496] (0.95 confidence)</td>
<td>1.34906 ns [1.3474, 1.34951] (0.95 confidence)</td>
<td>1.0034</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/matrix&#95;access</td>
<td>0.756793 ns [0.753797, 0.758451] (0.95 confidence)</td>
<td>0.747524 ns [0.746006, 0.750705] (0.95 confidence)</td>
<td>1.0124</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/matrix&#95;from&#95;fn</td>
<td>227.811 ns [227.153, 228.05] (0.95 confidence)</td>
<td>227.03 ns [226.443, 227.402] (0.95 confidence)</td>
<td>1.00344</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/matrix&#95;from&#95;rows</td>
<td>264.754 ns [264.181, 265.354] (0.95 confidence)</td>
<td>262.871 ns [262.6, 263.566] (0.95 confidence)</td>
<td>1.00716</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/matrix&#95;into&#95;rows</td>
<td>122.635 ns [122.073, 122.83] (0.95 confidence)</td>
<td>121.881 ns [121.486, 122.233] (0.95 confidence)</td>
<td>1.00618</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/matrix&#95;zero</td>
<td>6.222 ns [6.19667, 6.2306] (0.95 confidence)</td>
<td>6.14776 ns [6.13124, 6.15977] (0.95 confidence)</td>
<td>1.01208</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/vector&#95;access</td>
<td>0.755927 ns [0.754405, 0.757654] (0.95 confidence)</td>
<td>0.751778 ns [0.748226, 0.753623] (0.95 confidence)</td>
<td>1.00552</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/vector&#95;from&#95;fn</td>
<td>51.7425 ns [51.696, 51.8209] (0.95 confidence)</td>
<td>51.5327 ns [51.4278, 51.5847] (0.95 confidence)</td>
<td>1.00407</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/vector&#95;into&#95;array</td>
<td>22.2378 ns [22.1772, 22.3219] (0.95 confidence)</td>
<td>22.2115 ns [22.174, 22.2507] (0.95 confidence)</td>
<td>1.00119</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/vector&#95;new</td>
<td>51.4016 ns [51.343, 51.4324] (0.95 confidence)</td>
<td>51.1656 ns [51.1166, 51.2268] (0.95 confidence)</td>
<td>1.00461</td>
</tr>
<tr>
<td>public&#95;exact&#95;d3/vector&#95;zero</td>
<td>2.01457 ns [2.0117, 2.01797] (0.95 confidence)</td>
<td>2.00524 ns [2.00223, 2.00726] (0.95 confidence)</td>
<td>1.00465</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/matrix&#95;access</td>
<td>0.753548 ns [0.750771, 0.757605] (0.95 confidence)</td>
<td>0.75543 ns [0.75281, 0.758106] (0.95 confidence)</td>
<td>0.997509</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/matrix&#95;from&#95;fn</td>
<td>283.884 ns [283.613, 284.278] (0.95 confidence)</td>
<td>285.294 ns [284.834, 286.46] (0.95 confidence)</td>
<td>0.995056</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/matrix&#95;from&#95;rows</td>
<td>342.709 ns [341.969, 343.256] (0.95 confidence)</td>
<td>357.947 ns [354.309, 361.015] (0.95 confidence)</td>
<td>0.957429</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/matrix&#95;into&#95;rows</td>
<td>129.753 ns [129.497, 130.057] (0.95 confidence)</td>
<td>130.246 ns [130.009, 130.565] (0.95 confidence)</td>
<td>0.99622</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/matrix&#95;zero</td>
<td>15.1534 ns [15.0392, 15.2416] (0.95 confidence)</td>
<td>15.5442 ns [15.4621, 15.7621] (0.95 confidence)</td>
<td>0.974855</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/vector&#95;access</td>
<td>0.754366 ns [0.752125, 0.756903] (0.95 confidence)</td>
<td>0.754779 ns [0.753171, 0.755929] (0.95 confidence)</td>
<td>0.999453</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/vector&#95;from&#95;fn</td>
<td>65.7834 ns [65.7238, 65.8573] (0.95 confidence)</td>
<td>65.8668 ns [65.8025, 65.8919] (0.95 confidence)</td>
<td>0.998734</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/vector&#95;into&#95;array</td>
<td>29.61 ns [29.5423, 29.6579] (0.95 confidence)</td>
<td>29.6046 ns [29.5276, 29.6484] (0.95 confidence)</td>
<td>1.00018</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/vector&#95;new</td>
<td>64.9513 ns [64.8673, 65.0245] (0.95 confidence)</td>
<td>64.9722 ns [64.9111, 65.039] (0.95 confidence)</td>
<td>0.999679</td>
</tr>
<tr>
<td>public&#95;exact&#95;d4/vector&#95;zero</td>
<td>3.25081 ns [3.2431, 3.2669] (0.95 confidence)</td>
<td>3.38742 ns [3.3579, 3.4109] (0.95 confidence)</td>
<td>0.95967</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/matrix&#95;access</td>
<td>0.753978 ns [0.752362, 0.756477] (0.95 confidence)</td>
<td>0.755062 ns [0.751649, 0.75656] (0.95 confidence)</td>
<td>0.998564</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/matrix&#95;from&#95;fn</td>
<td>430.18 ns [428.079, 431.818] (0.95 confidence)</td>
<td>431.46 ns [427.586, 433.432] (0.95 confidence)</td>
<td>0.997034</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/matrix&#95;from&#95;rows</td>
<td>619.769 ns [618.813, 620.529] (0.95 confidence)</td>
<td>617.828 ns [616.809, 618.722] (0.95 confidence)</td>
<td>1.00314</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/matrix&#95;into&#95;rows</td>
<td>301.544 ns [300.339, 302.579] (0.95 confidence)</td>
<td>299.758 ns [298.451, 300.853] (0.95 confidence)</td>
<td>1.00596</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/matrix&#95;zero</td>
<td>21.5126 ns [21.4664, 21.5525] (0.95 confidence)</td>
<td>21.5942 ns [21.5047, 21.6502] (0.95 confidence)</td>
<td>0.996223</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/vector&#95;access</td>
<td>0.749983 ns [0.748414, 0.752059] (0.95 confidence)</td>
<td>0.751513 ns [0.749943, 0.754607] (0.95 confidence)</td>
<td>0.997965</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/vector&#95;from&#95;fn</td>
<td>81.8524 ns [81.7989, 82.0266] (0.95 confidence)</td>
<td>81.3726 ns [81.2948, 81.536] (0.95 confidence)</td>
<td>1.0059</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/vector&#95;into&#95;array</td>
<td>37.9326 ns [37.8915, 37.9682] (0.95 confidence)</td>
<td>37.84 ns [37.7326, 37.9744] (0.95 confidence)</td>
<td>1.00245</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/vector&#95;new</td>
<td>81.9319 ns [81.807, 82.0616] (0.95 confidence)</td>
<td>81.5746 ns [81.4565, 81.8208] (0.95 confidence)</td>
<td>1.00438</td>
</tr>
<tr>
<td>public&#95;exact&#95;d5/vector&#95;zero</td>
<td>4.23462 ns [4.21453, 4.24695] (0.95 confidence)</td>
<td>4.25737 ns [4.20793, 4.27825] (0.95 confidence)</td>
<td>0.994655</td>
</tr>
<tr>
<td>public&#95;exact&#95;scalar/rounded</td>
<td>5.5532 ns [5.54579, 5.55513] (0.95 confidence)</td>
<td>5.55743 ns [5.54671, 5.56066] (0.95 confidence)</td>
<td>0.999237</td>
</tr>
<tr>
<td>public&#95;exact&#95;scalar/strict</td>
<td>4.184 ns [4.17334, 4.19725] (0.95 confidence)</td>
<td>4.2527 ns [4.1755, 4.35369] (0.95 confidence)</td>
<td>0.983846</td>
</tr>
<tr>
<td>rational&#95;input&#95;d2/det&#95;row&#95;cleared&#95;bareiss</td>
<td>179.348 ns [179.117, 179.514] (0.95 confidence)</td>
<td>177.913 ns [177.563, 178.061] (0.95 confidence)</td>
<td>1.00807</td>
</tr>
<tr>
<td>rational&#95;input&#95;d2/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>76.4459 ns [76.307, 76.5091] (0.95 confidence)</td>
<td>76.4504 ns [76.3907, 76.5657] (0.95 confidence)</td>
<td>0.999941</td>
</tr>
<tr>
<td>rational&#95;input&#95;d2/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>622.931 ns [622.417, 623.64] (0.95 confidence)</td>
<td>617.845 ns [616.34, 619.403] (0.95 confidence)</td>
<td>1.00823</td>
</tr>
<tr>
<td>rational&#95;input&#95;d3/det&#95;row&#95;cleared&#95;bareiss</td>
<td>593.946 ns [593.274, 594.98] (0.95 confidence)</td>
<td>587.285 ns [585.535, 588.278] (0.95 confidence)</td>
<td>1.01134</td>
</tr>
<tr>
<td>rational&#95;input&#95;d3/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>368.446 ns [367.087, 369.359] (0.95 confidence)</td>
<td>368.486 ns [367.986, 369.732] (0.95 confidence)</td>
<td>0.999891</td>
</tr>
<tr>
<td>rational&#95;input&#95;d3/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>1262.73 ns [1260.03, 1266.5] (0.95 confidence)</td>
<td>1240.02 ns [1238.15, 1242.09] (0.95 confidence)</td>
<td>1.01831</td>
</tr>
<tr>
<td>rational&#95;input&#95;d4/det&#95;row&#95;cleared&#95;bareiss</td>
<td>1224.6 ns [1222.77, 1226.02] (0.95 confidence)</td>
<td>1208.81 ns [1206.67, 1210.71] (0.95 confidence)</td>
<td>1.01306</td>
</tr>
<tr>
<td>rational&#95;input&#95;d4/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>823.657 ns [822.459, 824.714] (0.95 confidence)</td>
<td>809.386 ns [805.941, 812.034] (0.95 confidence)</td>
<td>1.01763</td>
</tr>
<tr>
<td>rational&#95;input&#95;d4/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>2239.92 ns [2237.4, 2241.69] (0.95 confidence)</td>
<td>2217.36 ns [2214.59, 2220.3] (0.95 confidence)</td>
<td>1.01017</td>
</tr>
<tr>
<td>rational&#95;input&#95;d5/det&#95;row&#95;cleared&#95;bareiss</td>
<td>2109.81 ns [2106.27, 2112.62] (0.95 confidence)</td>
<td>2049.24 ns [2045.17, 2055.92] (0.95 confidence)</td>
<td>1.02956</td>
</tr>
<tr>
<td>rational&#95;input&#95;d5/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>1561.52 ns [1556.47, 1564.03] (0.95 confidence)</td>
<td>1547.16 ns [1545.56, 1550.05] (0.95 confidence)</td>
<td>1.00928</td>
</tr>
<tr>
<td>rational&#95;input&#95;d5/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>3374.19 ns [3369.75, 3376.45] (0.95 confidence)</td>
<td>3341.09 ns [3336.46, 3348.91] (0.95 confidence)</td>
<td>1.00991</td>
</tr>
<tr>
<td>rational&#95;input&#95;d6/det&#95;row&#95;cleared&#95;bareiss</td>
<td>3528.31 ns [3519.7, 3538.76] (0.95 confidence)</td>
<td>3436.7 ns [3429.31, 3442.18] (0.95 confidence)</td>
<td>1.02666</td>
</tr>
<tr>
<td>rational&#95;input&#95;d6/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>2827.81 ns [2825.39, 2831.86] (0.95 confidence)</td>
<td>2805.2 ns [2795.92, 2812.89] (0.95 confidence)</td>
<td>1.00806</td>
</tr>
<tr>
<td>rational&#95;input&#95;d6/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>5699.4 ns [5677.38, 5716.49] (0.95 confidence)</td>
<td>5634.1 ns [5626.11, 5639.42] (0.95 confidence)</td>
<td>1.01159</td>
</tr>
<tr>
<td>rational&#95;input&#95;d7/det&#95;row&#95;cleared&#95;bareiss</td>
<td>5547.98 ns [5541.64, 5558.96] (0.95 confidence)</td>
<td>5519.95 ns [5513.68, 5531.65] (0.95 confidence)</td>
<td>1.00508</td>
</tr>
<tr>
<td>rational&#95;input&#95;d7/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>4648.5 ns [4643.76, 4656.09] (0.95 confidence)</td>
<td>4597.77 ns [4588.8, 4603.12] (0.95 confidence)</td>
<td>1.01103</td>
</tr>
<tr>
<td>rational&#95;input&#95;d7/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>9855.63 ns [9847.2, 9872.56] (0.95 confidence)</td>
<td>9810.61 ns [9797.02, 9839.12] (0.95 confidence)</td>
<td>1.00459</td>
</tr>
<tr>
<td>rational&#95;input&#95;d8/det&#95;row&#95;cleared&#95;bareiss</td>
<td>8097.42 ns [8083.13, 8103.01] (0.95 confidence)</td>
<td>8118.07 ns [8092.15, 8146.57] (0.95 confidence)</td>
<td>0.997456</td>
</tr>
<tr>
<td>rational&#95;input&#95;d8/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>7263.88 ns [7247.13, 7273.49] (0.95 confidence)</td>
<td>7228.18 ns [7217.35, 7240.58] (0.95 confidence)</td>
<td>1.00494</td>
</tr>
<tr>
<td>rational&#95;input&#95;d8/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>15095.9 ns [15064.3, 15143.8] (0.95 confidence)</td>
<td>14961.9 ns [14952.1, 15000.2] (0.95 confidence)</td>
<td>1.00896</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d2/det&#95;row&#95;cleared&#95;bareiss</td>
<td>42103.3 ns [42034.7, 42187.6] (0.95 confidence)</td>
<td>41941.2 ns [41881.7, 41953.5] (0.95 confidence)</td>
<td>1.00387</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d2/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>943.264 ns [942.788, 943.904] (0.95 confidence)</td>
<td>919.332 ns [918.318, 921.772] (0.95 confidence)</td>
<td>1.02603</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d2/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>3036.77 ns [3032.2, 3041.08] (0.95 confidence)</td>
<td>3003.33 ns [2998.27, 3007.54] (0.95 confidence)</td>
<td>1.01113</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d3/det&#95;row&#95;cleared&#95;bareiss</td>
<td>68981.4 ns [68945, 69053] (0.95 confidence)</td>
<td>68880.9 ns [68780.1, 68981.7] (0.95 confidence)</td>
<td>1.00146</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d3/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>4182.64 ns [4177.78, 4190.66] (0.95 confidence)</td>
<td>4146.84 ns [4137.76, 4160.91] (0.95 confidence)</td>
<td>1.00863</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d3/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>11345.9 ns [11333, 11365.9] (0.95 confidence)</td>
<td>11136.7 ns [11118.5, 11148.2] (0.95 confidence)</td>
<td>1.01879</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d4/det&#95;row&#95;cleared&#95;bareiss</td>
<td>109960 ns [109787, 110243] (0.95 confidence)</td>
<td>109991 ns [109868, 110080] (0.95 confidence)</td>
<td>0.999716</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d4/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>9421.04 ns [9396.41, 9450.85] (0.95 confidence)</td>
<td>9293.41 ns [9286.4, 9306.65] (0.95 confidence)</td>
<td>1.01373</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d4/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>33485 ns [33423.4, 33527.8] (0.95 confidence)</td>
<td>32895.6 ns [32838.7, 32926.8] (0.95 confidence)</td>
<td>1.01792</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d5/det&#95;row&#95;cleared&#95;bareiss</td>
<td>185812 ns [185739, 186012] (0.95 confidence)</td>
<td>184590 ns [184107, 184821] (0.95 confidence)</td>
<td>1.00662</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d5/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>42908.3 ns [42777.8, 43013.8] (0.95 confidence)</td>
<td>41802.9 ns [41680, 41918.8] (0.95 confidence)</td>
<td>1.02644</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d5/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>69044.3 ns [68968.5, 69183.7] (0.95 confidence)</td>
<td>67570.7 ns [67405.7, 67722.9] (0.95 confidence)</td>
<td>1.02181</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d6/det&#95;row&#95;cleared&#95;bareiss</td>
<td>301738 ns [301183, 302390] (0.95 confidence)</td>
<td>296821 ns [296416, 297246] (0.95 confidence)</td>
<td>1.01656</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d6/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>104738 ns [104624, 104954] (0.95 confidence)</td>
<td>101069 ns [100929, 101239] (0.95 confidence)</td>
<td>1.03631</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d6/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>160249 ns [160160, 160662] (0.95 confidence)</td>
<td>155899 ns [155646, 156091] (0.95 confidence)</td>
<td>1.0279</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d7/det&#95;row&#95;cleared&#95;bareiss</td>
<td>501123 ns [500872, 501642] (0.95 confidence)</td>
<td>493028 ns [492327, 493966] (0.95 confidence)</td>
<td>1.01642</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d7/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>232579 ns [232067, 233066] (0.95 confidence)</td>
<td>224397 ns [223886, 224663] (0.95 confidence)</td>
<td>1.03646</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d7/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>337471 ns [337041, 338044] (0.95 confidence)</td>
<td>327709 ns [327112, 328670] (0.95 confidence)</td>
<td>1.02979</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d8/det&#95;row&#95;cleared&#95;bareiss</td>
<td>820753 ns [819943, 821683] (0.95 confidence)</td>
<td>806102 ns [805054, 807447] (0.95 confidence)</td>
<td>1.01817</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d8/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>463569 ns [462454, 464768] (0.95 confidence)</td>
<td>452446 ns [451326, 453565] (0.95 confidence)</td>
<td>1.02458</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide1024&#95;d8/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>646100 ns [644980, 647115] (0.95 confidence)</td>
<td>630696 ns [629008, 632391] (0.95 confidence)</td>
<td>1.02442</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d2/det&#95;row&#95;cleared&#95;bareiss</td>
<td>8774.79 ns [8751.15, 8805.73] (0.95 confidence)</td>
<td>8796.71 ns [8779.78, 8820.94] (0.95 confidence)</td>
<td>0.997508</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d2/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>643.344 ns [641.735, 646.037] (0.95 confidence)</td>
<td>646.735 ns [644.528, 648.869] (0.95 confidence)</td>
<td>0.994756</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d2/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>2075.07 ns [2061.33, 2087.85] (0.95 confidence)</td>
<td>2022.26 ns [2017.01, 2027.32] (0.95 confidence)</td>
<td>1.02612</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d3/det&#95;row&#95;cleared&#95;bareiss</td>
<td>13194.2 ns [13173.9, 13205.2] (0.95 confidence)</td>
<td>13289.5 ns [13270.8, 13316.3] (0.95 confidence)</td>
<td>0.992831</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d3/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>2751.4 ns [2744.38, 2762.44] (0.95 confidence)</td>
<td>2727.27 ns [2700.04, 2752.91] (0.95 confidence)</td>
<td>1.00884</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d3/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>5793.43 ns [5776.91, 5802.86] (0.95 confidence)</td>
<td>5742.14 ns [5725.49, 5756.62] (0.95 confidence)</td>
<td>1.00893</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d4/det&#95;row&#95;cleared&#95;bareiss</td>
<td>22008.5 ns [21971, 22041.7] (0.95 confidence)</td>
<td>21742 ns [21713, 21777] (0.95 confidence)</td>
<td>1.01226</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d4/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>5958.24 ns [5931.84, 5982.08] (0.95 confidence)</td>
<td>5810.67 ns [5804.13, 5824.45] (0.95 confidence)</td>
<td>1.0254</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d4/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>12838.3 ns [12814.1, 12859.7] (0.95 confidence)</td>
<td>12547.3 ns [12514.1, 12566.1] (0.95 confidence)</td>
<td>1.02319</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d5/det&#95;row&#95;cleared&#95;bareiss</td>
<td>33514.9 ns [33483.9, 33539.9] (0.95 confidence)</td>
<td>33510 ns [33458, 33560.5] (0.95 confidence)</td>
<td>1.00015</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d5/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>13392.1 ns [13380.5, 13426.2] (0.95 confidence)</td>
<td>13046.9 ns [13009.1, 13104.6] (0.95 confidence)</td>
<td>1.02646</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d5/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>20824.2 ns [20793.3, 20859.2] (0.95 confidence)</td>
<td>20260.7 ns [20221.1, 20318.9] (0.95 confidence)</td>
<td>1.02781</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d6/det&#95;row&#95;cleared&#95;bareiss</td>
<td>52030.7 ns [51938.3, 52124.3] (0.95 confidence)</td>
<td>51717.8 ns [51589.4, 51906.9] (0.95 confidence)</td>
<td>1.00605</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d6/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>26161.2 ns [26121.1, 26212] (0.95 confidence)</td>
<td>25617 ns [25574, 25684.2] (0.95 confidence)</td>
<td>1.02125</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d6/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>39501.5 ns [39416.4, 39573] (0.95 confidence)</td>
<td>38385.1 ns [38292.3, 38465.3] (0.95 confidence)</td>
<td>1.02908</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d7/det&#95;row&#95;cleared&#95;bareiss</td>
<td>83142.3 ns [83015.2, 83255.1] (0.95 confidence)</td>
<td>81995.2 ns [81846.6, 82088.3] (0.95 confidence)</td>
<td>1.01399</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d7/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>50156.9 ns [50104.8, 50254] (0.95 confidence)</td>
<td>48683.1 ns [48510.2, 48803.2] (0.95 confidence)</td>
<td>1.03027</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d7/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>72335.2 ns [72156.6, 72476.8] (0.95 confidence)</td>
<td>70105.2 ns [70036.5, 70328.8] (0.95 confidence)</td>
<td>1.03181</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d8/det&#95;row&#95;cleared&#95;bareiss</td>
<td>131374 ns [131089, 131471] (0.95 confidence)</td>
<td>129560 ns [129268, 129689] (0.95 confidence)</td>
<td>1.014</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d8/det&#95;sign&#95;row&#95;cleared&#95;bareiss</td>
<td>90340.3 ns [90200.1, 90486.8] (0.95 confidence)</td>
<td>89248.2 ns [89080.1, 89311] (0.95 confidence)</td>
<td>1.01224</td>
</tr>
<tr>
<td>rational&#95;input&#95;wide256&#95;d8/solve&#95;row&#95;cleared&#95;bareiss</td>
<td>126843 ns [126614, 127241] (0.95 confidence)</td>
<td>123295 ns [123026, 123509] (0.95 confidence)</td>
<td>1.02878</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d2/adapter</td>
<td>1023.07 ns [1018.64, 1026.73] (0.95 confidence)</td>
<td>1004.86 ns [1002.81, 1010] (0.95 confidence)</td>
<td>1.01812</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d2/construct</td>
<td>292.076 ns [291.301, 293.078] (0.95 confidence)</td>
<td>287.243 ns [286.359, 287.732] (0.95 confidence)</td>
<td>1.01682</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d2/prepared</td>
<td>652.568 ns [649.519, 654.727] (0.95 confidence)</td>
<td>635.393 ns [633.792, 639.39] (0.95 confidence)</td>
<td>1.02703</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d3/adapter</td>
<td>1959.33 ns [1955.5, 1966.31] (0.95 confidence)</td>
<td>1931.06 ns [1926.28, 1934.84] (0.95 confidence)</td>
<td>1.01464</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d3/construct</td>
<td>647.607 ns [646.559, 651.332] (0.95 confidence)</td>
<td>634.964 ns [633.542, 636.62] (0.95 confidence)</td>
<td>1.01991</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d3/prepared</td>
<td>1172 ns [1170.52, 1173.49] (0.95 confidence)</td>
<td>1143.5 ns [1142.25, 1145.02] (0.95 confidence)</td>
<td>1.02492</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d4/adapter</td>
<td>3299.59 ns [3290.71, 3303.18] (0.95 confidence)</td>
<td>3316.4 ns [3311.99, 3326.3] (0.95 confidence)</td>
<td>0.994931</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d4/construct</td>
<td>1052.44 ns [1043.69, 1066] (0.95 confidence)</td>
<td>1045.86 ns [1041.76, 1052.84] (0.95 confidence)</td>
<td>1.00629</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d4/prepared</td>
<td>2057.26 ns [2049.96, 2068.68] (0.95 confidence)</td>
<td>2114.81 ns [2111.1, 2118.88] (0.95 confidence)</td>
<td>0.972791</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d5/adapter</td>
<td>5396.08 ns [5385.27, 5406.18] (0.95 confidence)</td>
<td>5367.27 ns [5352.92, 5379.38] (0.95 confidence)</td>
<td>1.00537</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d5/construct</td>
<td>1927.28 ns [1921.14, 1931.49] (0.95 confidence)</td>
<td>1926.37 ns [1922.85, 1928.07] (0.95 confidence)</td>
<td>1.00047</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d5/prepared</td>
<td>3159.17 ns [3151.04, 3164.81] (0.95 confidence)</td>
<td>3153.48 ns [3143.05, 3165.09] (0.95 confidence)</td>
<td>1.0018</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d6/adapter</td>
<td>7471.63 ns [7459.16, 7485.66] (0.95 confidence)</td>
<td>7665.61 ns [7657.65, 7678.65] (0.95 confidence)</td>
<td>0.974694</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d6/construct</td>
<td>2505.46 ns [2501.44, 2507.96] (0.95 confidence)</td>
<td>2528.91 ns [2524.91, 2533.98] (0.95 confidence)</td>
<td>0.990729</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d6/prepared</td>
<td>4627.95 ns [4618.6, 4635.2] (0.95 confidence)</td>
<td>4717.03 ns [4711.04, 4728.2] (0.95 confidence)</td>
<td>0.981115</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d7/adapter</td>
<td>10265.4 ns [10245.3, 10280.1] (0.95 confidence)</td>
<td>10305.7 ns [10281.2, 10320] (0.95 confidence)</td>
<td>0.996091</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d7/construct</td>
<td>3460.77 ns [3441.85, 3496.44] (0.95 confidence)</td>
<td>3440.92 ns [3435.01, 3449.07] (0.95 confidence)</td>
<td>1.00577</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d7/prepared</td>
<td>6402.08 ns [6383.19, 6417.66] (0.95 confidence)</td>
<td>6578.61 ns [6557.09, 6594.06] (0.95 confidence)</td>
<td>0.973167</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d8/adapter</td>
<td>13247.9 ns [13233.7, 13272.4] (0.95 confidence)</td>
<td>13417.9 ns [13404.3, 13438.2] (0.95 confidence)</td>
<td>0.987333</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d8/construct</td>
<td>4042.24 ns [4033.58, 4050.53] (0.95 confidence)</td>
<td>4055.78 ns [4049.86, 4062.64] (0.95 confidence)</td>
<td>0.996661</td>
</tr>
<tr>
<td>rational&#95;solve/dense/dyadic&#95;d8/prepared</td>
<td>8951.38 ns [8935.27, 8973.06] (0.95 confidence)</td>
<td>9224.58 ns [9210.84, 9236.98] (0.95 confidence)</td>
<td>0.970383</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d2/adapter</td>
<td>1044.76 ns [1041.47, 1053.12] (0.95 confidence)</td>
<td>1021.8 ns [1019.75, 1025.85] (0.95 confidence)</td>
<td>1.02246</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d2/construct</td>
<td>328.182 ns [327.392, 329.094] (0.95 confidence)</td>
<td>323.505 ns [322.659, 324.458] (0.95 confidence)</td>
<td>1.01446</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d2/prepared</td>
<td>633.765 ns [633.121, 634.781] (0.95 confidence)</td>
<td>620.943 ns [619.542, 623.122] (0.95 confidence)</td>
<td>1.02065</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d3/adapter</td>
<td>2117.51 ns [2115.09, 2119.64] (0.95 confidence)</td>
<td>2073.99 ns [2071.39, 2075.9] (0.95 confidence)</td>
<td>1.02098</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d3/construct</td>
<td>683.064 ns [680.558, 686.197] (0.95 confidence)</td>
<td>665.361 ns [662.421, 667.888] (0.95 confidence)</td>
<td>1.02661</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d3/prepared</td>
<td>1299.73 ns [1297.35, 1301.4] (0.95 confidence)</td>
<td>1276.41 ns [1273.99, 1278.25] (0.95 confidence)</td>
<td>1.01827</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d4/adapter</td>
<td>3504.37 ns [3500.35, 3512.02] (0.95 confidence)</td>
<td>3476.28 ns [3471.03, 3482.07] (0.95 confidence)</td>
<td>1.00808</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d4/construct</td>
<td>1129.01 ns [1125.05, 1137.43] (0.95 confidence)</td>
<td>1104.67 ns [1102.17, 1107.9] (0.95 confidence)</td>
<td>1.02204</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d4/prepared</td>
<td>2188.19 ns [2185.52, 2189.47] (0.95 confidence)</td>
<td>2167.81 ns [2165.24, 2174.23] (0.95 confidence)</td>
<td>1.0094</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d5/adapter</td>
<td>5572.64 ns [5564.12, 5577.6] (0.95 confidence)</td>
<td>5508.96 ns [5495.5, 5517.6] (0.95 confidence)</td>
<td>1.01156</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d5/construct</td>
<td>1890.11 ns [1885.19, 1895.42] (0.95 confidence)</td>
<td>1869.52 ns [1858.84, 1874.59] (0.95 confidence)</td>
<td>1.01102</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d5/prepared</td>
<td>3412.97 ns [3409.44, 3417.01] (0.95 confidence)</td>
<td>3372.22 ns [3366.81, 3374.44] (0.95 confidence)</td>
<td>1.01208</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d6/adapter</td>
<td>7897.68 ns [7885.42, 7911.25] (0.95 confidence)</td>
<td>7885.78 ns [7849.24, 7907.89] (0.95 confidence)</td>
<td>1.00151</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d6/construct</td>
<td>2717.76 ns [2707.24, 2732.84] (0.95 confidence)</td>
<td>2725.14 ns [2713.05, 2745.79] (0.95 confidence)</td>
<td>0.997289</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d6/prepared</td>
<td>4864.74 ns [4856.61, 4870.35] (0.95 confidence)</td>
<td>4850.37 ns [4846.6, 4852.8] (0.95 confidence)</td>
<td>1.00296</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d7/adapter</td>
<td>11571.5 ns [11559.3, 11591.3] (0.95 confidence)</td>
<td>11418.2 ns [11390.4, 11442.2] (0.95 confidence)</td>
<td>1.01343</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d7/construct</td>
<td>3856.03 ns [3836.44, 3876.06] (0.95 confidence)</td>
<td>3789.99 ns [3754.44, 3816.52] (0.95 confidence)</td>
<td>1.01742</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d7/prepared</td>
<td>7310.46 ns [7288.17, 7338.63] (0.95 confidence)</td>
<td>7326.59 ns [7316.34, 7331.94] (0.95 confidence)</td>
<td>0.997799</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d8/adapter</td>
<td>15803.7 ns [15728.9, 15840.8] (0.95 confidence)</td>
<td>15832 ns [15780.5, 15862.7] (0.95 confidence)</td>
<td>0.99821</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d8/construct</td>
<td>5138.56 ns [5087.08, 5172.65] (0.95 confidence)</td>
<td>5087.43 ns [5045.89, 5110.96] (0.95 confidence)</td>
<td>1.01005</td>
</tr>
<tr>
<td>rational&#95;solve/dense/rational&#95;d8/prepared</td>
<td>10299.5 ns [10276.2, 10334.1] (0.95 confidence)</td>
<td>10391.8 ns [10348.8, 10415.5] (0.95 confidence)</td>
<td>0.991124</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d2/adapter</td>
<td>426.399 ns [424.787, 427.245] (0.95 confidence)</td>
<td>423.973 ns [422.789, 424.753] (0.95 confidence)</td>
<td>1.00572</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d2/construct</td>
<td>225.363 ns [225.121, 225.713] (0.95 confidence)</td>
<td>222.093 ns [221.733, 222.42] (0.95 confidence)</td>
<td>1.01472</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d2/prepared</td>
<td>151.436 ns [151.274, 151.624] (0.95 confidence)</td>
<td>149.522 ns [149.258, 149.688] (0.95 confidence)</td>
<td>1.0128</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d3/adapter</td>
<td>739.918 ns [739.1, 740.779] (0.95 confidence)</td>
<td>728.601 ns [727.028, 729.828] (0.95 confidence)</td>
<td>1.01553</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d3/construct</td>
<td>414.461 ns [413.649, 414.999] (0.95 confidence)</td>
<td>404.682 ns [403.881, 405.542] (0.95 confidence)</td>
<td>1.02416</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d3/prepared</td>
<td>246.651 ns [246.122, 247.291] (0.95 confidence)</td>
<td>243.458 ns [242.886, 243.965] (0.95 confidence)</td>
<td>1.01312</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d4/adapter</td>
<td>1033.18 ns [1032.27, 1033.94] (0.95 confidence)</td>
<td>1026.66 ns [1024.88, 1030.29] (0.95 confidence)</td>
<td>1.00635</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d4/construct</td>
<td>606.699 ns [604.787, 609.319] (0.95 confidence)</td>
<td>596.958 ns [592.423, 599.859] (0.95 confidence)</td>
<td>1.01632</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d4/prepared</td>
<td>322.772 ns [322.041, 323.129] (0.95 confidence)</td>
<td>321.563 ns [320.993, 321.997] (0.95 confidence)</td>
<td>1.00376</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d5/adapter</td>
<td>1368 ns [1366.99, 1369.59] (0.95 confidence)</td>
<td>1359.54 ns [1357.99, 1361.77] (0.95 confidence)</td>
<td>1.00623</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d5/construct</td>
<td>842.364 ns [839.566, 845.38] (0.95 confidence)</td>
<td>829.778 ns [828.51, 831.895] (0.95 confidence)</td>
<td>1.01517</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d5/prepared</td>
<td>416.867 ns [416.543, 417.205] (0.95 confidence)</td>
<td>416.859 ns [416.519, 417.386] (0.95 confidence)</td>
<td>1.00002</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d6/adapter</td>
<td>1731.77 ns [1727.57, 1735.68] (0.95 confidence)</td>
<td>1730.61 ns [1725.77, 1734.16] (0.95 confidence)</td>
<td>1.00067</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d6/construct</td>
<td>1102.09 ns [1099.35, 1105.05] (0.95 confidence)</td>
<td>1089.94 ns [1087.67, 1091.46] (0.95 confidence)</td>
<td>1.01115</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d6/prepared</td>
<td>503.63 ns [502.321, 504.2] (0.95 confidence)</td>
<td>503.515 ns [502.528, 504.161] (0.95 confidence)</td>
<td>1.00023</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d7/adapter</td>
<td>2150.95 ns [2148.85, 2152.88] (0.95 confidence)</td>
<td>2182.67 ns [2181.43, 2183.94] (0.95 confidence)</td>
<td>0.985469</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d7/construct</td>
<td>1402.65 ns [1399.57, 1410.47] (0.95 confidence)</td>
<td>1424.28 ns [1419.35, 1427.33] (0.95 confidence)</td>
<td>0.984813</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d7/prepared</td>
<td>599.796 ns [599.015, 600.551] (0.95 confidence)</td>
<td>602.411 ns [601.808, 603.854] (0.95 confidence)</td>
<td>0.995658</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d8/adapter</td>
<td>2452.39 ns [2442.04, 2461.66] (0.95 confidence)</td>
<td>2455.71 ns [2449.56, 2463.76] (0.95 confidence)</td>
<td>0.998648</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d8/construct</td>
<td>1682.87 ns [1678.77, 1690.16] (0.95 confidence)</td>
<td>1692.07 ns [1688.18, 1698.08] (0.95 confidence)</td>
<td>0.994561</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/dyadic&#95;d8/prepared</td>
<td>658.853 ns [658.229, 659.921] (0.95 confidence)</td>
<td>633.907 ns [632.997, 634.463] (0.95 confidence)</td>
<td>1.03935</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d2/adapter</td>
<td>576.414 ns [572.828, 583.909] (0.95 confidence)</td>
<td>561.544 ns [559.867, 565.681] (0.95 confidence)</td>
<td>1.02648</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d2/construct</td>
<td>305.67 ns [304.165, 307.931] (0.95 confidence)</td>
<td>301.452 ns [299.229, 302.972] (0.95 confidence)</td>
<td>1.01399</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d2/prepared</td>
<td>188.679 ns [188.271, 189.549] (0.95 confidence)</td>
<td>187.569 ns [186.771, 187.91] (0.95 confidence)</td>
<td>1.00592</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d3/adapter</td>
<td>916.744 ns [915.725, 917.636] (0.95 confidence)</td>
<td>897.258 ns [896.697, 898.657] (0.95 confidence)</td>
<td>1.02172</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d3/construct</td>
<td>537.099 ns [534.121, 539.4] (0.95 confidence)</td>
<td>525.214 ns [521.448, 526.935] (0.95 confidence)</td>
<td>1.02263</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d3/prepared</td>
<td>278.994 ns [278.714, 279.323] (0.95 confidence)</td>
<td>273.554 ns [273.217, 273.842] (0.95 confidence)</td>
<td>1.01989</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d4/adapter</td>
<td>1279.63 ns [1277.2, 1281.27] (0.95 confidence)</td>
<td>1246.41 ns [1244.14, 1247.63] (0.95 confidence)</td>
<td>1.02666</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d4/construct</td>
<td>777.72 ns [772.699, 781.676] (0.95 confidence)</td>
<td>761.476 ns [756.491, 765.02] (0.95 confidence)</td>
<td>1.02133</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d4/prepared</td>
<td>367.068 ns [366.298, 367.607] (0.95 confidence)</td>
<td>366.891 ns [366.121, 367.532] (0.95 confidence)</td>
<td>1.00048</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d5/adapter</td>
<td>1693.5 ns [1691.42, 1695.22] (0.95 confidence)</td>
<td>1678.7 ns [1676.1, 1680.26] (0.95 confidence)</td>
<td>1.00881</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d5/construct</td>
<td>1076.85 ns [1074.04, 1078.18] (0.95 confidence)</td>
<td>1072.51 ns [1069.5, 1076.4] (0.95 confidence)</td>
<td>1.00405</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d5/prepared</td>
<td>471.743 ns [471.227, 472.36] (0.95 confidence)</td>
<td>479.626 ns [478.14, 481.104] (0.95 confidence)</td>
<td>0.983566</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d6/adapter</td>
<td>2112.77 ns [2105.2, 2117.03] (0.95 confidence)</td>
<td>2068.15 ns [2063.62, 2071.82] (0.95 confidence)</td>
<td>1.02157</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d6/construct</td>
<td>1434.64 ns [1427.35, 1439.62] (0.95 confidence)</td>
<td>1363.32 ns [1357.27, 1370.27] (0.95 confidence)</td>
<td>1.05232</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d6/prepared</td>
<td>564.171 ns [563.239, 565.751] (0.95 confidence)</td>
<td>564.916 ns [562.782, 565.677] (0.95 confidence)</td>
<td>0.998682</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d7/adapter</td>
<td>2583.45 ns [2576.41, 2590.73] (0.95 confidence)</td>
<td>2572.93 ns [2568.3, 2575.25] (0.95 confidence)</td>
<td>1.00409</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d7/construct</td>
<td>1763.55 ns [1754.12, 1771.51] (0.95 confidence)</td>
<td>1769.15 ns [1756.23, 1774.39] (0.95 confidence)</td>
<td>0.996835</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d7/prepared</td>
<td>662.365 ns [660.613, 663.495] (0.95 confidence)</td>
<td>668.82 ns [667.837, 669.662] (0.95 confidence)</td>
<td>0.990348</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d8/adapter</td>
<td>2979.2 ns [2971.04, 2983.19] (0.95 confidence)</td>
<td>2968.05 ns [2963.48, 2975.78] (0.95 confidence)</td>
<td>1.00376</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d8/construct</td>
<td>2099.1 ns [2087.06, 2107.89] (0.95 confidence)</td>
<td>2096.45 ns [2090.06, 2100.52] (0.95 confidence)</td>
<td>1.00127</td>
</tr>
<tr>
<td>rational&#95;solve/diagonal/rational&#95;d8/prepared</td>
<td>745.979 ns [745.323, 747.917] (0.95 confidence)</td>
<td>757.942 ns [756.513, 759.675] (0.95 confidence)</td>
<td>0.984216</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d2/adapter</td>
<td>520.886 ns [516.361, 527.371] (0.95 confidence)</td>
<td>511.832 ns [509.509, 514.875] (0.95 confidence)</td>
<td>1.01769</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d2/construct</td>
<td>242.106 ns [241.469, 242.772] (0.95 confidence)</td>
<td>236.997 ns [236.784, 237.226] (0.95 confidence)</td>
<td>1.02156</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d2/prepared</td>
<td>220.041 ns [219.713, 220.303] (0.95 confidence)</td>
<td>215.659 ns [215.404, 215.982] (0.95 confidence)</td>
<td>1.02032</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d3/adapter</td>
<td>1036.55 ns [1034.96, 1040.1] (0.95 confidence)</td>
<td>1024.39 ns [1021.46, 1026.41] (0.95 confidence)</td>
<td>1.01188</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d3/construct</td>
<td>491.787 ns [490.005, 494.039] (0.95 confidence)</td>
<td>474.491 ns [473.848, 475.532] (0.95 confidence)</td>
<td>1.03645</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d3/prepared</td>
<td>454.734 ns [453.966, 455.297] (0.95 confidence)</td>
<td>450.057 ns [449.185, 450.693] (0.95 confidence)</td>
<td>1.01039</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d4/adapter</td>
<td>1620.03 ns [1618.58, 1621.77] (0.95 confidence)</td>
<td>1603.92 ns [1601.53, 1605.17] (0.95 confidence)</td>
<td>1.01004</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d4/construct</td>
<td>829.246 ns [826.973, 835.692] (0.95 confidence)</td>
<td>806.553 ns [804.69, 810.394] (0.95 confidence)</td>
<td>1.02814</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d4/prepared</td>
<td>661.915 ns [660.616, 663.07] (0.95 confidence)</td>
<td>659.754 ns [658.106, 662.859] (0.95 confidence)</td>
<td>1.00328</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d5/adapter</td>
<td>2420 ns [2417.68, 2424.88] (0.95 confidence)</td>
<td>2389.73 ns [2385.91, 2394.94] (0.95 confidence)</td>
<td>1.01267</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d5/construct</td>
<td>1256.38 ns [1252.03, 1259.03] (0.95 confidence)</td>
<td>1225.19 ns [1223.67, 1227.42] (0.95 confidence)</td>
<td>1.02546</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d5/prepared</td>
<td>983.346 ns [982.346, 985.026] (0.95 confidence)</td>
<td>981.252 ns [979.71, 983.06] (0.95 confidence)</td>
<td>1.00213</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d6/adapter</td>
<td>3250.57 ns [3247.56, 3256.3] (0.95 confidence)</td>
<td>3307.91 ns [3300.73, 3312.49] (0.95 confidence)</td>
<td>0.982668</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d6/construct</td>
<td>1690.62 ns [1688.92, 1694.9] (0.95 confidence)</td>
<td>1710.61 ns [1708, 1713.24] (0.95 confidence)</td>
<td>0.988312</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d6/prepared</td>
<td>1322.36 ns [1320.76, 1325.39] (0.95 confidence)</td>
<td>1331.01 ns [1329.54, 1332.23] (0.95 confidence)</td>
<td>0.993506</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d7/adapter</td>
<td>4294.48 ns [4289.02, 4299.63] (0.95 confidence)</td>
<td>4298.23 ns [4292.74, 4302.06] (0.95 confidence)</td>
<td>0.999127</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d7/construct</td>
<td>2255.21 ns [2250.09, 2266.67] (0.95 confidence)</td>
<td>2274.57 ns [2269.72, 2278.75] (0.95 confidence)</td>
<td>0.991488</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d7/prepared</td>
<td>1744.05 ns [1742.11, 1746.66] (0.95 confidence)</td>
<td>1743.56 ns [1741.62, 1746.86] (0.95 confidence)</td>
<td>1.00028</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d8/adapter</td>
<td>5155.78 ns [5144.9, 5161.57] (0.95 confidence)</td>
<td>5250.54 ns [5244.95, 5255.84] (0.95 confidence)</td>
<td>0.981953</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d8/construct</td>
<td>2875.18 ns [2870.25, 2881.58] (0.95 confidence)</td>
<td>2885.56 ns [2882.52, 2891.24] (0.95 confidence)</td>
<td>0.996403</td>
</tr>
<tr>
<td>rational&#95;solve/lower/dyadic&#95;d8/prepared</td>
<td>2102.32 ns [2098.62, 2104.6] (0.95 confidence)</td>
<td>2115.27 ns [2111.29, 2120.96] (0.95 confidence)</td>
<td>0.993875</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d2/adapter</td>
<td>649.071 ns [647.931, 650.222] (0.95 confidence)</td>
<td>637.21 ns [636.347, 638.884] (0.95 confidence)</td>
<td>1.01861</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d2/construct</td>
<td>344.836 ns [342.932, 346.106] (0.95 confidence)</td>
<td>345.244 ns [343.809, 348.621] (0.95 confidence)</td>
<td>0.998818</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d2/prepared</td>
<td>213.178 ns [212.71, 213.657] (0.95 confidence)</td>
<td>212.651 ns [212.263, 212.947] (0.95 confidence)</td>
<td>1.00248</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d3/adapter</td>
<td>1156.56 ns [1154.96, 1158.91] (0.95 confidence)</td>
<td>1130.73 ns [1128.06, 1132.8] (0.95 confidence)</td>
<td>1.02284</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d3/construct</td>
<td>666.528 ns [664.322, 668.634] (0.95 confidence)</td>
<td>656.39 ns [653.428, 659.243] (0.95 confidence)</td>
<td>1.01544</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d3/prepared</td>
<td>368.236 ns [367.666, 368.963] (0.95 confidence)</td>
<td>361.298 ns [360.789, 361.865] (0.95 confidence)</td>
<td>1.0192</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d4/adapter</td>
<td>1838.14 ns [1835.76, 1839.88] (0.95 confidence)</td>
<td>1781.52 ns [1779.04, 1785.43] (0.95 confidence)</td>
<td>1.03178</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d4/construct</td>
<td>1061.07 ns [1057.26, 1074.62] (0.95 confidence)</td>
<td>1051.1 ns [1046.82, 1058.24] (0.95 confidence)</td>
<td>1.00948</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d4/prepared</td>
<td>592.314 ns [591.358, 593.15] (0.95 confidence)</td>
<td>587.296 ns [586.422, 588.126] (0.95 confidence)</td>
<td>1.00854</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d5/adapter</td>
<td>2629.05 ns [2622.01, 2634.94] (0.95 confidence)</td>
<td>2571.18 ns [2568.63, 2574.03] (0.95 confidence)</td>
<td>1.02251</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d5/construct</td>
<td>1578.76 ns [1569.14, 1588.88] (0.95 confidence)</td>
<td>1570.81 ns [1556.33, 1580.91] (0.95 confidence)</td>
<td>1.00506</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d5/prepared</td>
<td>814.075 ns [811.948, 815.679] (0.95 confidence)</td>
<td>807.661 ns [805.643, 810.879] (0.95 confidence)</td>
<td>1.00794</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d6/adapter</td>
<td>3543.5 ns [3539.59, 3548.56] (0.95 confidence)</td>
<td>3503.92 ns [3498.13, 3510.93] (0.95 confidence)</td>
<td>1.0113</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d6/construct</td>
<td>2265.42 ns [2252.37, 2272.16] (0.95 confidence)</td>
<td>2272.39 ns [2267.13, 2275.89] (0.95 confidence)</td>
<td>0.996933</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d6/prepared</td>
<td>1023.17 ns [1021.16, 1025.71] (0.95 confidence)</td>
<td>1031.57 ns [1028.6, 1037.32] (0.95 confidence)</td>
<td>0.991858</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d7/adapter</td>
<td>4684.7 ns [4673.44, 4695.93] (0.95 confidence)</td>
<td>4643.63 ns [4638.54, 4651.44] (0.95 confidence)</td>
<td>1.00885</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d7/construct</td>
<td>2997.62 ns [2976.62, 3012.05] (0.95 confidence)</td>
<td>3025.35 ns [3014.9, 3031.98] (0.95 confidence)</td>
<td>0.990833</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d7/prepared</td>
<td>1352.23 ns [1349.5, 1354.02] (0.95 confidence)</td>
<td>1364.76 ns [1363.75, 1367.64] (0.95 confidence)</td>
<td>0.99082</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d8/adapter</td>
<td>5754.25 ns [5740.16, 5770.77] (0.95 confidence)</td>
<td>5706.37 ns [5699.68, 5718.52] (0.95 confidence)</td>
<td>1.00839</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d8/construct</td>
<td>3842.75 ns [3820.73, 3858.17] (0.95 confidence)</td>
<td>3841.45 ns [3822.1, 3854.1] (0.95 confidence)</td>
<td>1.00034</td>
</tr>
<tr>
<td>rational&#95;solve/lower/rational&#95;d8/prepared</td>
<td>1603.88 ns [1600.38, 1608.76] (0.95 confidence)</td>
<td>1621.7 ns [1617.83, 1627.88] (0.95 confidence)</td>
<td>0.989012</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d2/adapter</td>
<td>521.93 ns [518.105, 524.933] (0.95 confidence)</td>
<td>512.741 ns [511.401, 516.991] (0.95 confidence)</td>
<td>1.01792</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d2/construct</td>
<td>239.92 ns [239.58, 240.234] (0.95 confidence)</td>
<td>238.96 ns [238.69, 239.291] (0.95 confidence)</td>
<td>1.00402</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d2/prepared</td>
<td>220.461 ns [219.94, 221.202] (0.95 confidence)</td>
<td>214.945 ns [214.681, 215.469] (0.95 confidence)</td>
<td>1.02566</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d3/adapter</td>
<td>1039.74 ns [1035.26, 1044.56] (0.95 confidence)</td>
<td>1031.19 ns [1029.08, 1033.93] (0.95 confidence)</td>
<td>1.00829</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d3/construct</td>
<td>489.27 ns [488.403, 489.788] (0.95 confidence)</td>
<td>477.037 ns [475.208, 479.244] (0.95 confidence)</td>
<td>1.02564</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d3/prepared</td>
<td>451.331 ns [450.701, 451.899] (0.95 confidence)</td>
<td>447.781 ns [447.006, 448.353] (0.95 confidence)</td>
<td>1.00793</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d4/adapter</td>
<td>1600.4 ns [1596.12, 1602.57] (0.95 confidence)</td>
<td>1606.86 ns [1603.56, 1608.11] (0.95 confidence)</td>
<td>0.995977</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d4/construct</td>
<td>819.007 ns [814.48, 822.851] (0.95 confidence)</td>
<td>801.771 ns [799.763, 804.933] (0.95 confidence)</td>
<td>1.0215</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d4/prepared</td>
<td>654.378 ns [653.843, 655.902] (0.95 confidence)</td>
<td>655.427 ns [653.995, 656.7] (0.95 confidence)</td>
<td>0.998399</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d5/adapter</td>
<td>2411.03 ns [2408.73, 2414.55] (0.95 confidence)</td>
<td>2395.27 ns [2390.27, 2398.64] (0.95 confidence)</td>
<td>1.00658</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d5/construct</td>
<td>1226.96 ns [1222.22, 1229.74] (0.95 confidence)</td>
<td>1220.66 ns [1218.25, 1222.32] (0.95 confidence)</td>
<td>1.00516</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d5/prepared</td>
<td>976.065 ns [975.175, 977.207] (0.95 confidence)</td>
<td>981.648 ns [978.678, 983.603] (0.95 confidence)</td>
<td>0.994313</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d6/adapter</td>
<td>3254.8 ns [3249.4, 3259.45] (0.95 confidence)</td>
<td>3322.56 ns [3318.63, 3327.16] (0.95 confidence)</td>
<td>0.979605</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d6/construct</td>
<td>1695.04 ns [1693.78, 1697.59] (0.95 confidence)</td>
<td>1713.94 ns [1711.51, 1717.82] (0.95 confidence)</td>
<td>0.988969</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d6/prepared</td>
<td>1320.88 ns [1319.21, 1323.21] (0.95 confidence)</td>
<td>1335.22 ns [1333.39, 1336.58] (0.95 confidence)</td>
<td>0.989261</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d7/adapter</td>
<td>4306.74 ns [4300.86, 4311.76] (0.95 confidence)</td>
<td>4286.79 ns [4282.51, 4290.09] (0.95 confidence)</td>
<td>1.00465</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d7/construct</td>
<td>2241.81 ns [2237.42, 2249.81] (0.95 confidence)</td>
<td>2266.64 ns [2262.6, 2272.22] (0.95 confidence)</td>
<td>0.989047</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d7/prepared</td>
<td>1744.67 ns [1742.2, 1746.55] (0.95 confidence)</td>
<td>1737.22 ns [1733.44, 1739.93] (0.95 confidence)</td>
<td>1.00429</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d8/adapter</td>
<td>5180.37 ns [5171.3, 5191.45] (0.95 confidence)</td>
<td>5270.5 ns [5265, 5276.94] (0.95 confidence)</td>
<td>0.982898</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d8/construct</td>
<td>2876.54 ns [2870.17, 2882.58] (0.95 confidence)</td>
<td>2868.17 ns [2864.81, 2873.15] (0.95 confidence)</td>
<td>1.00292</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/dyadic&#95;d8/prepared</td>
<td>2125.96 ns [2120.21, 2128.02] (0.95 confidence)</td>
<td>2137.99 ns [2134.72, 2140.78] (0.95 confidence)</td>
<td>0.994373</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d2/adapter</td>
<td>645.21 ns [643.236, 651.93] (0.95 confidence)</td>
<td>640.051 ns [639.035, 642.668] (0.95 confidence)</td>
<td>1.00806</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d2/construct</td>
<td>350.296 ns [348.513, 352.933] (0.95 confidence)</td>
<td>346.894 ns [344.7, 348.312] (0.95 confidence)</td>
<td>1.00981</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d2/prepared</td>
<td>215.068 ns [214.737, 215.507] (0.95 confidence)</td>
<td>221.166 ns [219.961, 221.884] (0.95 confidence)</td>
<td>0.972426</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d3/adapter</td>
<td>1152.98 ns [1150.24, 1154.17] (0.95 confidence)</td>
<td>1126.72 ns [1124.66, 1128.87] (0.95 confidence)</td>
<td>1.02331</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d3/construct</td>
<td>665.952 ns [661.587, 668.511] (0.95 confidence)</td>
<td>654.135 ns [648.541, 657.28] (0.95 confidence)</td>
<td>1.01807</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d3/prepared</td>
<td>365.555 ns [364.729, 366.275] (0.95 confidence)</td>
<td>361.362 ns [360.894, 361.828] (0.95 confidence)</td>
<td>1.0116</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d4/adapter</td>
<td>1821.47 ns [1817.08, 1824.06] (0.95 confidence)</td>
<td>1786.27 ns [1782.33, 1788.69] (0.95 confidence)</td>
<td>1.01971</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d4/construct</td>
<td>1063.99 ns [1057.68, 1080.05] (0.95 confidence)</td>
<td>1051.92 ns [1046.67, 1056.35] (0.95 confidence)</td>
<td>1.01148</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d4/prepared</td>
<td>597.852 ns [596.939, 598.568] (0.95 confidence)</td>
<td>589.494 ns [588.553, 590.509] (0.95 confidence)</td>
<td>1.01418</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d5/adapter</td>
<td>2596.86 ns [2592.68, 2601.15] (0.95 confidence)</td>
<td>2562.31 ns [2559.74, 2565.43] (0.95 confidence)</td>
<td>1.01348</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d5/construct</td>
<td>1576.6 ns [1569.54, 1583.8] (0.95 confidence)</td>
<td>1568.19 ns [1561.52, 1572.61] (0.95 confidence)</td>
<td>1.00536</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d5/prepared</td>
<td>808.891 ns [807.807, 810.169] (0.95 confidence)</td>
<td>807.532 ns [805.826, 809.082] (0.95 confidence)</td>
<td>1.00168</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d6/adapter</td>
<td>3538.41 ns [3534.67, 3544.45] (0.95 confidence)</td>
<td>3532.72 ns [3526.87, 3540.96] (0.95 confidence)</td>
<td>1.00161</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d6/construct</td>
<td>2261.42 ns [2256.84, 2266.89] (0.95 confidence)</td>
<td>2273.52 ns [2267.53, 2279.63] (0.95 confidence)</td>
<td>0.99468</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d6/prepared</td>
<td>1017.36 ns [1014.75, 1020.45] (0.95 confidence)</td>
<td>1029.43 ns [1026.82, 1031.82] (0.95 confidence)</td>
<td>0.988273</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d7/adapter</td>
<td>4670.67 ns [4665.18, 4678.15] (0.95 confidence)</td>
<td>4586.5 ns [4583.6, 4591.99] (0.95 confidence)</td>
<td>1.01835</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d7/construct</td>
<td>2986.97 ns [2975.88, 3006.61] (0.95 confidence)</td>
<td>2984.86 ns [2973.01, 2992.31] (0.95 confidence)</td>
<td>1.00071</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d7/prepared</td>
<td>1355.32 ns [1353.16, 1357.41] (0.95 confidence)</td>
<td>1366.59 ns [1363.1, 1370.64] (0.95 confidence)</td>
<td>0.991754</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d8/adapter</td>
<td>5796.24 ns [5780.17, 5808.47] (0.95 confidence)</td>
<td>5726.86 ns [5717.12, 5733.95] (0.95 confidence)</td>
<td>1.01212</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d8/construct</td>
<td>3825.28 ns [3819.18, 3844.74] (0.95 confidence)</td>
<td>3832.59 ns [3824.17, 3842.26] (0.95 confidence)</td>
<td>0.998093</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;lower/rational&#95;d8/prepared</td>
<td>1619.66 ns [1614.09, 1623.9] (0.95 confidence)</td>
<td>1635.39 ns [1629.43, 1645.81] (0.95 confidence)</td>
<td>0.990377</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d2/adapter</td>
<td>526.953 ns [524.337, 533.378] (0.95 confidence)</td>
<td>520.913 ns [517.508, 527.227] (0.95 confidence)</td>
<td>1.0116</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d2/construct</td>
<td>239.743 ns [239.373, 240.064] (0.95 confidence)</td>
<td>233.836 ns [233.512, 234.178] (0.95 confidence)</td>
<td>1.02526</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d2/prepared</td>
<td>220.129 ns [219.802, 220.648] (0.95 confidence)</td>
<td>216.428 ns [216.173, 216.726] (0.95 confidence)</td>
<td>1.0171</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d3/adapter</td>
<td>1014.87 ns [1014.26, 1016.4] (0.95 confidence)</td>
<td>1008.6 ns [1007.03, 1010.31] (0.95 confidence)</td>
<td>1.00622</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d3/construct</td>
<td>474.873 ns [473.703, 475.799] (0.95 confidence)</td>
<td>467.634 ns [466.577, 468.286] (0.95 confidence)</td>
<td>1.01548</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d3/prepared</td>
<td>441.73 ns [440.935, 442.263] (0.95 confidence)</td>
<td>432.173 ns [431.553, 432.725] (0.95 confidence)</td>
<td>1.02211</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d4/adapter</td>
<td>1647.76 ns [1644.43, 1651.34] (0.95 confidence)</td>
<td>1634.05 ns [1631.68, 1636.31] (0.95 confidence)</td>
<td>1.00839</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d4/construct</td>
<td>825.368 ns [822.396, 834.55] (0.95 confidence)</td>
<td>812.518 ns [809.988, 816.674] (0.95 confidence)</td>
<td>1.01581</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d4/prepared</td>
<td>691.522 ns [690.9, 692.543] (0.95 confidence)</td>
<td>688.923 ns [686.619, 691.016] (0.95 confidence)</td>
<td>1.00377</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d5/adapter</td>
<td>2325.32 ns [2320.19, 2330.15] (0.95 confidence)</td>
<td>2327.31 ns [2324.73, 2330.89] (0.95 confidence)</td>
<td>0.999143</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d5/construct</td>
<td>1189.6 ns [1186.85, 1194.39] (0.95 confidence)</td>
<td>1190.57 ns [1187.98, 1191.83] (0.95 confidence)</td>
<td>0.999181</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d5/prepared</td>
<td>946.481 ns [945.686, 947.475] (0.95 confidence)</td>
<td>943.516 ns [942.043, 944.706] (0.95 confidence)</td>
<td>1.00314</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d6/adapter</td>
<td>3315.31 ns [3308.75, 3321.91] (0.95 confidence)</td>
<td>3393.71 ns [3389.25, 3400.05] (0.95 confidence)</td>
<td>0.976898</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d6/construct</td>
<td>1758.76 ns [1754.54, 1762.94] (0.95 confidence)</td>
<td>1779.07 ns [1775.7, 1780.87] (0.95 confidence)</td>
<td>0.988588</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d6/prepared</td>
<td>1360.06 ns [1358.29, 1360.84] (0.95 confidence)</td>
<td>1365.72 ns [1363.71, 1367.16] (0.95 confidence)</td>
<td>0.995859</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d7/adapter</td>
<td>4298.44 ns [4291.24, 4303.71] (0.95 confidence)</td>
<td>4381.5 ns [4374.98, 4390.68] (0.95 confidence)</td>
<td>0.981043</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d7/construct</td>
<td>2306.65 ns [2302.34, 2316.8] (0.95 confidence)</td>
<td>2330.18 ns [2323.04, 2333.07] (0.95 confidence)</td>
<td>0.989903</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d7/prepared</td>
<td>1775.18 ns [1773.48, 1776.5] (0.95 confidence)</td>
<td>1771.76 ns [1769.56, 1774.52] (0.95 confidence)</td>
<td>1.00193</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d8/adapter</td>
<td>5394.56 ns [5374.8, 5414.04] (0.95 confidence)</td>
<td>5367.83 ns [5356.45, 5382.23] (0.95 confidence)</td>
<td>1.00498</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d8/construct</td>
<td>3026.63 ns [3022.73, 3031.7] (0.95 confidence)</td>
<td>3036.43 ns [3029.34, 3040.6] (0.95 confidence)</td>
<td>0.996773</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/dyadic&#95;d8/prepared</td>
<td>2079.44 ns [2075.54, 2082.63] (0.95 confidence)</td>
<td>2081.63 ns [2078.69, 2085.84] (0.95 confidence)</td>
<td>0.998951</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d2/adapter</td>
<td>683.579 ns [681.696, 686.227] (0.95 confidence)</td>
<td>672.743 ns [671.242, 675.012] (0.95 confidence)</td>
<td>1.01611</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d2/construct</td>
<td>375.625 ns [374.454, 376.103] (0.95 confidence)</td>
<td>369.534 ns [367.816, 370.528] (0.95 confidence)</td>
<td>1.01648</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d2/prepared</td>
<td>229.297 ns [228.654, 229.968] (0.95 confidence)</td>
<td>226.061 ns [225.565, 226.536] (0.95 confidence)</td>
<td>1.01431</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d3/adapter</td>
<td>1182.39 ns [1180.18, 1185.26] (0.95 confidence)</td>
<td>1157.17 ns [1154.59, 1158.29] (0.95 confidence)</td>
<td>1.0218</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d3/construct</td>
<td>687.584 ns [685.671, 689.649] (0.95 confidence)</td>
<td>677.065 ns [674.606, 678.49] (0.95 confidence)</td>
<td>1.01554</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d3/prepared</td>
<td>385.232 ns [384.397, 385.783] (0.95 confidence)</td>
<td>377.704 ns [377.358, 378.448] (0.95 confidence)</td>
<td>1.01993</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d4/adapter</td>
<td>1768.86 ns [1767.2, 1771.04] (0.95 confidence)</td>
<td>1730.6 ns [1728.39, 1734.71] (0.95 confidence)</td>
<td>1.02211</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d4/construct</td>
<td>1004.94 ns [1000.42, 1025.22] (0.95 confidence)</td>
<td>983.017 ns [976.155, 986.32] (0.95 confidence)</td>
<td>1.02231</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d4/prepared</td>
<td>624.376 ns [623.034, 625.446] (0.95 confidence)</td>
<td>619.581 ns [618.138, 620.792] (0.95 confidence)</td>
<td>1.00774</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d5/adapter</td>
<td>2495.69 ns [2492.89, 2497.44] (0.95 confidence)</td>
<td>2469.04 ns [2465.64, 2471.69] (0.95 confidence)</td>
<td>1.0108</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d5/construct</td>
<td>1441.19 ns [1436.94, 1444.17] (0.95 confidence)</td>
<td>1417.91 ns [1412.8, 1424.78] (0.95 confidence)</td>
<td>1.01641</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d5/prepared</td>
<td>868.547 ns [867.187, 869.555] (0.95 confidence)</td>
<td>873.484 ns [870.827, 874.665] (0.95 confidence)</td>
<td>0.994348</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d6/adapter</td>
<td>3363.75 ns [3355.98, 3373.11] (0.95 confidence)</td>
<td>3355.25 ns [3349.24, 3363.04] (0.95 confidence)</td>
<td>1.00254</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d6/construct</td>
<td>2005.48 ns [1995.35, 2013.91] (0.95 confidence)</td>
<td>2011.01 ns [2002.19, 2020.5] (0.95 confidence)</td>
<td>0.997249</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d6/prepared</td>
<td>1124.88 ns [1123.49, 1126.4] (0.95 confidence)</td>
<td>1136.92 ns [1134.29, 1139.59] (0.95 confidence)</td>
<td>0.989418</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d7/adapter</td>
<td>4430.95 ns [4422.81, 4437.31] (0.95 confidence)</td>
<td>4387.15 ns [4383.79, 4394.73] (0.95 confidence)</td>
<td>1.00998</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d7/construct</td>
<td>2701.59 ns [2695.62, 2711.64] (0.95 confidence)</td>
<td>2712.33 ns [2704.51, 2722.56] (0.95 confidence)</td>
<td>0.99604</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d7/prepared</td>
<td>1455.86 ns [1454.57, 1458.49] (0.95 confidence)</td>
<td>1465.86 ns [1464.57, 1468.3] (0.95 confidence)</td>
<td>0.993177</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d8/adapter</td>
<td>5617.51 ns [5612.62, 5623.6] (0.95 confidence)</td>
<td>5637.12 ns [5622.79, 5643.57] (0.95 confidence)</td>
<td>0.996521</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d8/construct</td>
<td>3474.56 ns [3466.69, 3489.79] (0.95 confidence)</td>
<td>3471.53 ns [3458.46, 3488.43] (0.95 confidence)</td>
<td>1.00087</td>
</tr>
<tr>
<td>rational&#95;solve/permuted&#95;upper/rational&#95;d8/prepared</td>
<td>1849.72 ns [1846.68, 1852.23] (0.95 confidence)</td>
<td>1862.19 ns [1858.77, 1866.47] (0.95 confidence)</td>
<td>0.993305</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d2/adapter</td>
<td>1010.1 ns [1007.97, 1018.09] (0.95 confidence)</td>
<td>995.89 ns [994.707, 998.736] (0.95 confidence)</td>
<td>1.01427</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d2/construct</td>
<td>294.141 ns [293.438, 294.602] (0.95 confidence)</td>
<td>286.561 ns [285.824, 287.085] (0.95 confidence)</td>
<td>1.02645</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d2/prepared</td>
<td>635.22 ns [634.015, 636.07] (0.95 confidence)</td>
<td>628.889 ns [628.151, 629.916] (0.95 confidence)</td>
<td>1.01007</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d3/adapter</td>
<td>1814.85 ns [1810.89, 1819.75] (0.95 confidence)</td>
<td>1780.57 ns [1778.05, 1782.9] (0.95 confidence)</td>
<td>1.01925</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d3/construct</td>
<td>584.026 ns [582.92, 585.011] (0.95 confidence)</td>
<td>570.444 ns [569.219, 571.68] (0.95 confidence)</td>
<td>1.02381</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d3/prepared</td>
<td>1101.95 ns [1100.57, 1103.21] (0.95 confidence)</td>
<td>1075.24 ns [1073.49, 1077.89] (0.95 confidence)</td>
<td>1.02484</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d4/adapter</td>
<td>2690.46 ns [2682.36, 2696.26] (0.95 confidence)</td>
<td>2677.59 ns [2672.67, 2683.54] (0.95 confidence)</td>
<td>1.00481</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d4/construct</td>
<td>803.655 ns [796.666, 820.045] (0.95 confidence)</td>
<td>786.828 ns [784.939, 789.984] (0.95 confidence)</td>
<td>1.02139</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d4/prepared</td>
<td>1739.08 ns [1736.5, 1740.75] (0.95 confidence)</td>
<td>1726.72 ns [1724.73, 1729.53] (0.95 confidence)</td>
<td>1.00716</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d5/adapter</td>
<td>4126.49 ns [4114.31, 4137.46] (0.95 confidence)</td>
<td>4141.67 ns [4135.9, 4144.51] (0.95 confidence)</td>
<td>0.996334</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d5/construct</td>
<td>1312.69 ns [1308.59, 1315.87] (0.95 confidence)</td>
<td>1307.76 ns [1306.15, 1314.33] (0.95 confidence)</td>
<td>1.00377</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d5/prepared</td>
<td>2608.33 ns [2604.12, 2618.67] (0.95 confidence)</td>
<td>2617.68 ns [2609.78, 2622.02] (0.95 confidence)</td>
<td>0.996427</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d6/adapter</td>
<td>5493.61 ns [5490.2, 5502.67] (0.95 confidence)</td>
<td>5585.76 ns [5575.21, 5590.92] (0.95 confidence)</td>
<td>0.983503</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d6/construct</td>
<td>1664.35 ns [1658.62, 1676.12] (0.95 confidence)</td>
<td>1676.64 ns [1674.77, 1679.14] (0.95 confidence)</td>
<td>0.992667</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d6/prepared</td>
<td>3601.6 ns [3596.84, 3605.37] (0.95 confidence)</td>
<td>3624.39 ns [3621.07, 3629.56] (0.95 confidence)</td>
<td>0.993712</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d7/adapter</td>
<td>7710.3 ns [7703.42, 7721.68] (0.95 confidence)</td>
<td>7730.46 ns [7709.08, 7736.6] (0.95 confidence)</td>
<td>0.997392</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d7/construct</td>
<td>2146.06 ns [2141.95, 2173.09] (0.95 confidence)</td>
<td>2170.41 ns [2163.18, 2181.21] (0.95 confidence)</td>
<td>0.988781</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d7/prepared</td>
<td>5210.61 ns [5199.57, 5219.82] (0.95 confidence)</td>
<td>5192.74 ns [5184.17, 5199.32] (0.95 confidence)</td>
<td>1.00344</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d8/adapter</td>
<td>9400.88 ns [9387.15, 9417.26] (0.95 confidence)</td>
<td>9348.29 ns [9334.32, 9357.1] (0.95 confidence)</td>
<td>1.00563</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d8/construct</td>
<td>2358.36 ns [2354.09, 2362.22] (0.95 confidence)</td>
<td>2374.3 ns [2366.74, 2378.72] (0.95 confidence)</td>
<td>0.99329</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/dyadic&#95;d8/prepared</td>
<td>6901.21 ns [6892.78, 6914.49] (0.95 confidence)</td>
<td>6852.61 ns [6842.82, 6860.45] (0.95 confidence)</td>
<td>1.00709</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d2/adapter</td>
<td>1170.17 ns [1166.8, 1176.16] (0.95 confidence)</td>
<td>1147.3 ns [1145.51, 1151.26] (0.95 confidence)</td>
<td>1.01993</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d2/construct</td>
<td>478.794 ns [477.544, 480.14] (0.95 confidence)</td>
<td>470.238 ns [469.246, 471.521] (0.95 confidence)</td>
<td>1.0182</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d2/prepared</td>
<td>586.391 ns [585.355, 587.31] (0.95 confidence)</td>
<td>577.975 ns [576.26, 579.276] (0.95 confidence)</td>
<td>1.01456</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d3/adapter</td>
<td>1923.21 ns [1921.03, 1925.07] (0.95 confidence)</td>
<td>1888.87 ns [1886.73, 1892.01] (0.95 confidence)</td>
<td>1.01818</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d3/construct</td>
<td>740.809 ns [738.152, 743.689] (0.95 confidence)</td>
<td>731.945 ns [730.037, 733.825] (0.95 confidence)</td>
<td>1.01211</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d3/prepared</td>
<td>1029.3 ns [1028.5, 1030.17] (0.95 confidence)</td>
<td>1014.62 ns [1012.81, 1016.47] (0.95 confidence)</td>
<td>1.01447</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d4/adapter</td>
<td>3102.88 ns [3098.2, 3106.1] (0.95 confidence)</td>
<td>3050.3 ns [3045.76, 3054.42] (0.95 confidence)</td>
<td>1.01724</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d4/construct</td>
<td>1262.73 ns [1258.63, 1270.99] (0.95 confidence)</td>
<td>1252.63 ns [1247.85, 1257.12] (0.95 confidence)</td>
<td>1.00806</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d4/prepared</td>
<td>1640.87 ns [1639.25, 1642.21] (0.95 confidence)</td>
<td>1624.55 ns [1622.43, 1626.28] (0.95 confidence)</td>
<td>1.01005</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d5/adapter</td>
<td>4327.75 ns [4323.9, 4331.91] (0.95 confidence)</td>
<td>4279.09 ns [4271.38, 4287.84] (0.95 confidence)</td>
<td>1.01137</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d5/construct</td>
<td>1643.51 ns [1636.98, 1646.98] (0.95 confidence)</td>
<td>1639.43 ns [1636.23, 1642.3] (0.95 confidence)</td>
<td>1.00249</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d5/prepared</td>
<td>2449.01 ns [2445.1, 2452.77] (0.95 confidence)</td>
<td>2434.27 ns [2430.65, 2439.28] (0.95 confidence)</td>
<td>1.00606</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d6/adapter</td>
<td>5874.21 ns [5864.09, 5885.61] (0.95 confidence)</td>
<td>5919.61 ns [5913.35, 5924.17] (0.95 confidence)</td>
<td>0.992331</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d6/construct</td>
<td>2231.61 ns [2225.21, 2241.64] (0.95 confidence)</td>
<td>2246.52 ns [2243.09, 2249.78] (0.95 confidence)</td>
<td>0.993362</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d6/prepared</td>
<td>3404.73 ns [3402.23, 3409.91] (0.95 confidence)</td>
<td>3432.65 ns [3426.58, 3439.23] (0.95 confidence)</td>
<td>0.991866</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d7/adapter</td>
<td>8446.31 ns [8439.91, 8463.82] (0.95 confidence)</td>
<td>8423.63 ns [8412.02, 8433.62] (0.95 confidence)</td>
<td>1.00269</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d7/construct</td>
<td>2819.12 ns [2809.51, 2830.27] (0.95 confidence)</td>
<td>2804.6 ns [2793.8, 2814.75] (0.95 confidence)</td>
<td>1.00518</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d7/prepared</td>
<td>5377.26 ns [5370.05, 5386.3] (0.95 confidence)</td>
<td>5364.61 ns [5360.77, 5372.4] (0.95 confidence)</td>
<td>1.00236</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d8/adapter</td>
<td>10344.4 ns [10324.5, 10359.1] (0.95 confidence)</td>
<td>10317.7 ns [10302, 10332.3] (0.95 confidence)</td>
<td>1.00258</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d8/construct</td>
<td>3355.43 ns [3346.04, 3367.16] (0.95 confidence)</td>
<td>3342.77 ns [3336.83, 3354.49] (0.95 confidence)</td>
<td>1.00379</td>
</tr>
<tr>
<td>rational&#95;solve/sparse/rational&#95;d8/prepared</td>
<td>6750.83 ns [6739.2, 6764.79] (0.95 confidence)</td>
<td>6691.93 ns [6668.55, 6717.54] (0.95 confidence)</td>
<td>1.0088</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d2/adapter</td>
<td>518.629 ns [515.232, 528.072] (0.95 confidence)</td>
<td>521.121 ns [516.032, 524.123] (0.95 confidence)</td>
<td>0.995218</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d2/construct</td>
<td>239.359 ns [238.858, 239.642] (0.95 confidence)</td>
<td>234.787 ns [234.301, 235.186] (0.95 confidence)</td>
<td>1.01947</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d2/prepared</td>
<td>217.245 ns [216.945, 217.596] (0.95 confidence)</td>
<td>215.965 ns [215.693, 216.245] (0.95 confidence)</td>
<td>1.00593</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d3/adapter</td>
<td>1013.41 ns [1011, 1016.91] (0.95 confidence)</td>
<td>1010.58 ns [1009.11, 1012.11] (0.95 confidence)</td>
<td>1.0028</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d3/construct</td>
<td>476.297 ns [475.459, 477.516] (0.95 confidence)</td>
<td>471.216 ns [470.312, 472.244] (0.95 confidence)</td>
<td>1.01078</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d3/prepared</td>
<td>440.268 ns [439.754, 440.832] (0.95 confidence)</td>
<td>435.176 ns [434.48, 435.589] (0.95 confidence)</td>
<td>1.0117</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d4/adapter</td>
<td>1639.99 ns [1637.34, 1642.73] (0.95 confidence)</td>
<td>1645.94 ns [1642.07, 1649.67] (0.95 confidence)</td>
<td>0.996381</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d4/construct</td>
<td>827.983 ns [824.952, 831.817] (0.95 confidence)</td>
<td>816.493 ns [810.796, 829.545] (0.95 confidence)</td>
<td>1.01407</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d4/prepared</td>
<td>688.781 ns [687.837, 689.491] (0.95 confidence)</td>
<td>690.375 ns [687.665, 694.344] (0.95 confidence)</td>
<td>0.99769</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d5/adapter</td>
<td>2347.22 ns [2341.3, 2351.98] (0.95 confidence)</td>
<td>2333.2 ns [2330.12, 2336.34] (0.95 confidence)</td>
<td>1.00601</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d5/construct</td>
<td>1223.66 ns [1220.48, 1228.1] (0.95 confidence)</td>
<td>1197.3 ns [1194.68, 1205.58] (0.95 confidence)</td>
<td>1.02202</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d5/prepared</td>
<td>953.326 ns [951.857, 954.024] (0.95 confidence)</td>
<td>946.054 ns [945.046, 947.001] (0.95 confidence)</td>
<td>1.00769</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d6/adapter</td>
<td>3324.46 ns [3318.68, 3330.18] (0.95 confidence)</td>
<td>3398.22 ns [3390.92, 3402.13] (0.95 confidence)</td>
<td>0.978295</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d6/construct</td>
<td>1752.96 ns [1748.28, 1760.45] (0.95 confidence)</td>
<td>1760.04 ns [1757.21, 1763.82] (0.95 confidence)</td>
<td>0.995976</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d6/prepared</td>
<td>1360.48 ns [1358.74, 1361.8] (0.95 confidence)</td>
<td>1364.12 ns [1363.1, 1366.14] (0.95 confidence)</td>
<td>0.997333</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d7/adapter</td>
<td>4313.82 ns [4306.41, 4325.01] (0.95 confidence)</td>
<td>4389.54 ns [4384.77, 4393.35] (0.95 confidence)</td>
<td>0.982749</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d7/construct</td>
<td>2353.43 ns [2345.2, 2362.94] (0.95 confidence)</td>
<td>2332.08 ns [2329.65, 2336.89] (0.95 confidence)</td>
<td>1.00916</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d7/prepared</td>
<td>1778.04 ns [1775.56, 1779.67] (0.95 confidence)</td>
<td>1784.97 ns [1782.99, 1786.74] (0.95 confidence)</td>
<td>0.99612</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d8/adapter</td>
<td>5347.66 ns [5342.08, 5359.65] (0.95 confidence)</td>
<td>5386.33 ns [5372.02, 5391.12] (0.95 confidence)</td>
<td>0.99282</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d8/construct</td>
<td>3029.07 ns [3023.6, 3036.71] (0.95 confidence)</td>
<td>3030.07 ns [3021.52, 3042.19] (0.95 confidence)</td>
<td>0.999669</td>
</tr>
<tr>
<td>rational&#95;solve/upper/dyadic&#95;d8/prepared</td>
<td>2079.14 ns [2076.85, 2081.73] (0.95 confidence)</td>
<td>2086.45 ns [2082.25, 2089.31] (0.95 confidence)</td>
<td>0.996495</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d2/adapter</td>
<td>679.787 ns [678.283, 682.623] (0.95 confidence)</td>
<td>673.806 ns [671.78, 676.002] (0.95 confidence)</td>
<td>1.00888</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d2/construct</td>
<td>371.58 ns [369.663, 372.783] (0.95 confidence)</td>
<td>372.594 ns [370.856, 374.714] (0.95 confidence)</td>
<td>0.997278</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d2/prepared</td>
<td>226.113 ns [225.751, 226.883] (0.95 confidence)</td>
<td>225.596 ns [225.223, 226.148] (0.95 confidence)</td>
<td>1.00229</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d3/adapter</td>
<td>1190.67 ns [1188.87, 1191.58] (0.95 confidence)</td>
<td>1158.09 ns [1156.93, 1160.24] (0.95 confidence)</td>
<td>1.02814</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d3/construct</td>
<td>692.707 ns [691.752, 694.242] (0.95 confidence)</td>
<td>679.89 ns [676.57, 682.953] (0.95 confidence)</td>
<td>1.01885</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d3/prepared</td>
<td>384.626 ns [383.81, 385.391] (0.95 confidence)</td>
<td>379.094 ns [378.481, 379.766] (0.95 confidence)</td>
<td>1.01459</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d4/adapter</td>
<td>1782.71 ns [1781.17, 1784.43] (0.95 confidence)</td>
<td>1740.24 ns [1738.34, 1743.1] (0.95 confidence)</td>
<td>1.02441</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d4/construct</td>
<td>1008.39 ns [1003.81, 1012.58] (0.95 confidence)</td>
<td>996.891 ns [990.586, 1001.17] (0.95 confidence)</td>
<td>1.01153</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d4/prepared</td>
<td>625.528 ns [624.039, 627.208] (0.95 confidence)</td>
<td>625.342 ns [624.054, 626.489] (0.95 confidence)</td>
<td>1.0003</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d5/adapter</td>
<td>2523.05 ns [2516.9, 2529.09] (0.95 confidence)</td>
<td>2464.16 ns [2458.55, 2470.75] (0.95 confidence)</td>
<td>1.0239</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d5/construct</td>
<td>1436.27 ns [1434.19, 1439.21] (0.95 confidence)</td>
<td>1419.32 ns [1416.54, 1422.44] (0.95 confidence)</td>
<td>1.01195</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d5/prepared</td>
<td>868.945 ns [867.332, 870.052] (0.95 confidence)</td>
<td>870.065 ns [868.672, 870.852] (0.95 confidence)</td>
<td>0.998712</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d6/adapter</td>
<td>3365.14 ns [3357.02, 3373.82] (0.95 confidence)</td>
<td>3350.67 ns [3344.97, 3359.61] (0.95 confidence)</td>
<td>1.00432</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d6/construct</td>
<td>1988.31 ns [1984.44, 1993.68] (0.95 confidence)</td>
<td>1975.54 ns [1972.41, 1979.97] (0.95 confidence)</td>
<td>1.00647</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d6/prepared</td>
<td>1123.24 ns [1119.62, 1125.21] (0.95 confidence)</td>
<td>1136.62 ns [1135.39, 1137.73] (0.95 confidence)</td>
<td>0.988231</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d7/adapter</td>
<td>4461.1 ns [4452.76, 4466.36] (0.95 confidence)</td>
<td>4436.66 ns [4431.25, 4444.17] (0.95 confidence)</td>
<td>1.00551</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d7/construct</td>
<td>2708.54 ns [2703.26, 2717.87] (0.95 confidence)</td>
<td>2715.38 ns [2709.11, 2722.64] (0.95 confidence)</td>
<td>0.99748</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d7/prepared</td>
<td>1461.56 ns [1458.46, 1465.31] (0.95 confidence)</td>
<td>1465.72 ns [1464.69, 1467.67] (0.95 confidence)</td>
<td>0.997166</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d8/adapter</td>
<td>5625.36 ns [5620.53, 5642.26] (0.95 confidence)</td>
<td>5644.13 ns [5636.4, 5654.42] (0.95 confidence)</td>
<td>0.996674</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d8/construct</td>
<td>3503.46 ns [3489.08, 3514.83] (0.95 confidence)</td>
<td>3497.53 ns [3489.1, 3505.66] (0.95 confidence)</td>
<td>1.0017</td>
</tr>
<tr>
<td>rational&#95;solve/upper/rational&#95;d8/prepared</td>
<td>1853.3 ns [1850.44, 1856.48] (0.95 confidence)</td>
<td>1866 ns [1862.5, 1867.95] (0.95 confidence)</td>
<td>0.993194</td>
</tr>
</tbody>
</table>

## lu-solve

Semantic contract: la-stack-lu-solve-v1. Collection session: `537edaa10b434ff88fae9345fd8d0260`.
Collection interval (Unix ns): 1794918237502732000–1794918874570519000.
Acquisition: newly collected.

- baseline: available — inside declared release interval; release v0.5.3; revision `7165ed562f2e4611d201e574fc1f600a8ad8b33e`.
  Source SHA-256: `073502f95ae58344aa776b54303ae24d1e6cf2ce60d1437d950bcbe64d0b4dbf`; harness SHA-256:
  `ba63ffff14e6124d5aca8df516e6836a606ad10809bbe3b1ab13dd5fd7ceedbf`; scientific gate: passed.
- current: available — inside declared release interval; release v0.5.4; revision `531f75eb58eb8ac810ad5d89bf1d3e66a301e84e`.
  Source SHA-256: `50d72eb42299ad56e00e6838b34b10824c560a3eb5ae9c9a05e48859fd3cfe46`; harness SHA-256:
  `ba63ffff14e6124d5aca8df516e6836a606ad10809bbe3b1ab13dd5fd7ceedbf`; scientific gate: passed.

### Mean (ns)

| Workload | Baseline | Current | Baseline / current |
| --- | --- | --- | --- |
| d16/faer&#95;lu&#95;solve | 869.828 ns [868.213, 871.432] (0.95 confidence) | 868.144 ns [866.47, 869.843] (0.95 confidence) | 1.00194 |
| d16/la&#95;stack&#95;lu&#95;solve | 649.62 ns [643.409, 656.046] (0.95 confidence) | 646.963 ns [640.761, 653.414] (0.95 confidence) | 1.00411 |
| d16/nalgebra&#95;lu&#95;solve | 570.337 ns [569.637, 570.991] (0.95 confidence) | 564.66 ns [563.864, 565.392] (0.95 confidence) | 1.01005 |
| d2/faer&#95;lu&#95;solve | 127.11 ns [126.907, 127.316] (0.95 confidence) | 128.077 ns [127.61, 128.549] (0.95 confidence) | 0.992449 |
| d2/la&#95;stack&#95;lu&#95;solve | 2.07694 ns [2.06394, 2.09078] (0.95 confidence) | 2.05833 ns [2.04736, 2.07029] (0.95 confidence) | 1.00904 |
| d2/nalgebra&#95;lu&#95;solve | 4.54168 ns [4.53387, 4.54994] (0.95 confidence) | 4.47262 ns [4.45864, 4.48834] (0.95 confidence) | 1.01544 |
| d3/faer&#95;lu&#95;solve | 157.552 ns [157.293, 157.817] (0.95 confidence) | 155.882 ns [155.543, 156.233] (0.95 confidence) | 1.01071 |
| d3/la&#95;stack&#95;lu&#95;solve | 9.94496 ns [9.92704, 9.96212] (0.95 confidence) | 9.80129 ns [9.7868, 9.81598] (0.95 confidence) | 1.01466 |
| d3/nalgebra&#95;lu&#95;solve | 32.8507 ns [32.6918, 33.0151] (0.95 confidence) | 31.6962 ns [31.5534, 31.8406] (0.95 confidence) | 1.03642 |
| d32/faer&#95;lu&#95;solve | 2939.33 ns [2933.65, 2945.65] (0.95 confidence) | 2920.89 ns [2914.33, 2928.6] (0.95 confidence) | 1.00631 |
| d32/la&#95;stack&#95;lu&#95;solve | 2713.07 ns [2709.36, 2716.79] (0.95 confidence) | 2705.88 ns [2701.65, 2710] (0.95 confidence) | 1.00266 |
| d32/nalgebra&#95;lu&#95;solve | 2684.13 ns [2674.13, 2694.54] (0.95 confidence) | 2665.88 ns [2656.4, 2673.88] (0.95 confidence) | 1.00685 |
| d4/faer&#95;lu&#95;solve | 181.202 ns [180.97, 181.442] (0.95 confidence) | 180.759 ns [180.374, 181.145] (0.95 confidence) | 1.00245 |
| d4/la&#95;stack&#95;lu&#95;solve | 23.2935 ns [23.2609, 23.3258] (0.95 confidence) | 23.0641 ns [23.028, 23.098] (0.95 confidence) | 1.00995 |
| d4/nalgebra&#95;lu&#95;solve | 52.0021 ns [51.839, 52.1543] (0.95 confidence) | 49.968 ns [49.8627, 50.0705] (0.95 confidence) | 1.04071 |
| d5/faer&#95;lu&#95;solve | 248.552 ns [247.955, 249.079] (0.95 confidence) | 248.29 ns [247.699, 248.824] (0.95 confidence) | 1.00106 |
| d5/la&#95;stack&#95;lu&#95;solve | 47.561 ns [47.3879, 47.7307] (0.95 confidence) | 46.9038 ns [46.7417, 47.0677] (0.95 confidence) | 1.01401 |
| d5/nalgebra&#95;lu&#95;solve | 73.0724 ns [72.8535, 73.2944] (0.95 confidence) | 71.2598 ns [71.0808, 71.4414] (0.95 confidence) | 1.02544 |
| d64/faer&#95;lu&#95;solve | 12677.4 ns [12654.8, 12701.2] (0.95 confidence) | 12523.6 ns [12495.6, 12552.8] (0.95 confidence) | 1.01228 |
| d64/la&#95;stack&#95;lu&#95;solve | 17804.8 ns [17734.6, 17871.8] (0.95 confidence) | 17826 ns [17767.9, 17883.8] (0.95 confidence) | 0.998812 |
| d64/nalgebra&#95;lu&#95;solve | 14333.2 ns [14315.8, 14349.6] (0.95 confidence) | 14384.7 ns [14323, 14450.5] (0.95 confidence) | 0.996418 |
| d8/faer&#95;lu&#95;solve | 338.208 ns [337.334, 339.101] (0.95 confidence) | 330.951 ns [330.18, 331.836] (0.95 confidence) | 1.02193 |
| d8/la&#95;stack&#95;lu&#95;solve | 129.467 ns [129.27, 129.659] (0.95 confidence) | 128.351 ns [128.109, 128.582] (0.95 confidence) | 1.00869 |
| d8/nalgebra&#95;lu&#95;solve | 163.725 ns [163.264, 164.327] (0.95 confidence) | 163.131 ns [162.874, 163.379] (0.95 confidence) | 1.00364 |

### Median (ns)

| Workload | Baseline | Current | Baseline / current |
| --- | --- | --- | --- |
| d16/faer&#95;lu&#95;solve | 870.621 ns [869.275, 872.121] (0.95 confidence) | 868.305 ns [865.979, 870.438] (0.95 confidence) | 1.00267 |
| d16/la&#95;stack&#95;lu&#95;solve | 650.728 ns [644.065, 653.174] (0.95 confidence) | 648.762 ns [640.587, 650.581] (0.95 confidence) | 1.00303 |
| d16/nalgebra&#95;lu&#95;solve | 571.205 ns [570.733, 571.738] (0.95 confidence) | 565.132 ns [564.443, 565.645] (0.95 confidence) | 1.01075 |
| d2/faer&#95;lu&#95;solve | 127.178 ns [126.866, 127.274] (0.95 confidence) | 127.938 ns [127.501, 128.618] (0.95 confidence) | 0.994058 |
| d2/la&#95;stack&#95;lu&#95;solve | 2.0402 ns [2.03662, 2.04702] (0.95 confidence) | 2.03568 ns [2.03521, 2.03768] (0.95 confidence) | 1.00222 |
| d2/nalgebra&#95;lu&#95;solve | 4.53023 ns [4.52293, 4.53838] (0.95 confidence) | 4.45772 ns [4.45145, 4.46162] (0.95 confidence) | 1.01627 |
| d3/faer&#95;lu&#95;solve | 157.36 ns [157.122, 157.735] (0.95 confidence) | 155.795 ns [155.273, 156.394] (0.95 confidence) | 1.01004 |
| d3/la&#95;stack&#95;lu&#95;solve | 9.961 ns [9.92805, 9.99721] (0.95 confidence) | 9.80694 ns [9.78888, 9.82592] (0.95 confidence) | 1.01571 |
| d3/nalgebra&#95;lu&#95;solve | 32.5441 ns [32.466, 32.6572] (0.95 confidence) | 31.3746 ns [31.2597, 32.2559] (0.95 confidence) | 1.03728 |
| d32/faer&#95;lu&#95;solve | 2935.09 ns [2928.56, 2940.8] (0.95 confidence) | 2910.33 ns [2906.48, 2915.35] (0.95 confidence) | 1.00851 |
| d32/la&#95;stack&#95;lu&#95;solve | 2711.57 ns [2706.63, 2715.56] (0.95 confidence) | 2707.39 ns [2698.58, 2714.77] (0.95 confidence) | 1.00154 |
| d32/nalgebra&#95;lu&#95;solve | 2689.27 ns [2682.69, 2694.48] (0.95 confidence) | 2675.44 ns [2668.25, 2681.65] (0.95 confidence) | 1.00517 |
| d4/faer&#95;lu&#95;solve | 181.258 ns [180.834, 181.465] (0.95 confidence) | 180.604 ns [180.112, 181.261] (0.95 confidence) | 1.00362 |
| d4/la&#95;stack&#95;lu&#95;solve | 23.3115 ns [23.2599, 23.3393] (0.95 confidence) | 23.0854 ns [23.0361, 23.1269] (0.95 confidence) | 1.00979 |
| d4/nalgebra&#95;lu&#95;solve | 52.0524 ns [51.9979, 52.196] (0.95 confidence) | 50.005 ns [49.8659, 50.167] (0.95 confidence) | 1.04094 |
| d5/faer&#95;lu&#95;solve | 248.945 ns [248.274, 249.433] (0.95 confidence) | 248.634 ns [248.008, 249.091] (0.95 confidence) | 1.00125 |
| d5/la&#95;stack&#95;lu&#95;solve | 47.9257 ns [47.3619, 48.1396] (0.95 confidence) | 46.8005 ns [46.6284, 46.996] (0.95 confidence) | 1.02404 |
| d5/nalgebra&#95;lu&#95;solve | 72.7499 ns [72.4063, 73.5066] (0.95 confidence) | 71.0579 ns [70.8848, 71.2155] (0.95 confidence) | 1.02381 |
| d64/faer&#95;lu&#95;solve | 12676.9 ns [12660.1, 12685.1] (0.95 confidence) | 12540.4 ns [12526.9, 12563.7] (0.95 confidence) | 1.01089 |
| d64/la&#95;stack&#95;lu&#95;solve | 17864.5 ns [17707, 17945] (0.95 confidence) | 17915.3 ns [17788.1, 17947.7] (0.95 confidence) | 0.997164 |
| d64/nalgebra&#95;lu&#95;solve | 14344.1 ns [14321.2, 14366.1] (0.95 confidence) | 14234.3 ns [14216, 14273.7] (0.95 confidence) | 1.00772 |
| d8/faer&#95;lu&#95;solve | 338.034 ns [336.773, 339.285] (0.95 confidence) | 330.342 ns [329.317, 331.311] (0.95 confidence) | 1.02329 |
| d8/la&#95;stack&#95;lu&#95;solve | 129.515 ns [129.216, 129.772] (0.95 confidence) | 128.594 ns [128.299, 128.732] (0.95 confidence) | 1.00716 |
| d8/nalgebra&#95;lu&#95;solve | 163.308 ns [163.189, 163.801] (0.95 confidence) | 163.316 ns [163.059, 163.603] (0.95 confidence) | 0.999951 |

## column-basis

Semantic contract: la-stack-column-basis-v1. Collection session: `9bf9c769dfe84e3ea4a3b8861bc924e8`.
Collection interval (Unix ns): 1794918875581667000–1794919996636100000.
Acquisition: newly collected.

- baseline: available — inside declared release interval; release v0.5.3; revision `7165ed562f2e4611d201e574fc1f600a8ad8b33e`.
  Source SHA-256: `073502f95ae58344aa776b54303ae24d1e6cf2ce60d1437d950bcbe64d0b4dbf`; harness SHA-256:
  `ba63ffff14e6124d5aca8df516e6836a606ad10809bbe3b1ab13dd5fd7ceedbf`; scientific gate: passed.
- current: available — inside declared release interval; release v0.5.4; revision `531f75eb58eb8ac810ad5d89bf1d3e66a301e84e`.
  Source SHA-256: `50d72eb42299ad56e00e6838b34b10824c560a3eb5ae9c9a05e48859fd3cfe46`; harness SHA-256:
  `ba63ffff14e6124d5aca8df516e6836a606ad10809bbe3b1ab13dd5fd7ceedbf`; scientific gate: passed.

### Mean (ns)

| Workload | Baseline | Current | Baseline / current |
| --- | --- | --- | --- |
| column&#95;basis/2x2/full&#95;rank | 4.95966 ns [4.94994, 4.9696] (0.95 confidence) | 4.83783 ns [4.8286, 4.84687] (0.95 confidence) | 1.02518 |
| column&#95;basis/2x2/large | 4.9199 ns [4.91243, 4.92747] (0.95 confidence) | 4.87045 ns [4.86067, 4.87991] (0.95 confidence) | 1.01015 |
| column&#95;basis/2x2/near&#95;dependent | 4.96558 ns [4.95636, 4.97498] (0.95 confidence) | 4.90075 ns [4.88876, 4.91695] (0.95 confidence) | 1.01323 |
| column&#95;basis/2x2/rank&#95;one | 4.73316 ns [4.72474, 4.74156] (0.95 confidence) | 4.64712 ns [4.63737, 4.65656] (0.95 confidence) | 1.01851 |
| column&#95;basis/2x2/small | 4.92644 ns [4.91662, 4.93611] (0.95 confidence) | 4.84088 ns [4.83417, 4.84737] (0.95 confidence) | 1.01767 |
| column&#95;basis/2x3/shared&#95;face | 9.9587 ns [9.91878, 9.99947] (0.95 confidence) | 10.2796 ns [10.2445, 10.3177] (0.95 confidence) | 0.968783 |
| column&#95;basis/3x3/full&#95;rank | 14.8223 ns [14.2363, 15.479] (0.95 confidence) | 15.8602 ns [15.1224, 16.6175] (0.95 confidence) | 0.934564 |
| column&#95;basis/3x3/large | 14.7988 ns [14.2192, 15.4476] (0.95 confidence) | 15.8881 ns [15.155, 16.6315] (0.95 confidence) | 0.931442 |
| column&#95;basis/3x3/near&#95;dependent | 14.802 ns [14.2231, 15.4579] (0.95 confidence) | 15.9199 ns [15.1901, 16.6643] (0.95 confidence) | 0.929781 |
| column&#95;basis/3x3/rank&#95;one | 9.04429 ns [8.96469, 9.12088] (0.95 confidence) | 9.13933 ns [9.08734, 9.18681] (0.95 confidence) | 0.9896 |
| column&#95;basis/3x3/small | 14.8081 ns [14.2152, 15.4697] (0.95 confidence) | 16.035 ns [15.3053, 16.7884] (0.95 confidence) | 0.923485 |
| column&#95;basis/3x4/shared&#95;face | 17.616 ns [17.3365, 17.974] (0.95 confidence) | 18.2426 ns [18.0663, 18.4473] (0.95 confidence) | 0.965649 |
| column&#95;basis/4x4/full&#95;rank | 29.74 ns [29.6323, 29.8366] (0.95 confidence) | 29.3611 ns [29.2671, 29.4505] (0.95 confidence) | 1.01291 |
| column&#95;basis/4x4/large | 29.7475 ns [29.6478, 29.8425] (0.95 confidence) | 29.5345 ns [29.4386, 29.6295] (0.95 confidence) | 1.00721 |
| column&#95;basis/4x4/near&#95;dependent | 30.0041 ns [29.9133, 30.1023] (0.95 confidence) | 29.5666 ns [29.4749, 29.6567] (0.95 confidence) | 1.0148 |
| column&#95;basis/4x4/rank&#95;one | 14.5818 ns [14.5552, 14.6086] (0.95 confidence) | 14.5347 ns [14.5115, 14.5606] (0.95 confidence) | 1.00324 |
| column&#95;basis/4x4/small | 30.5895 ns [30.3986, 30.7912] (0.95 confidence) | 29.4893 ns [29.3911, 29.5807] (0.95 confidence) | 1.03731 |
| column&#95;basis/4x5/shared&#95;face | 39.2914 ns [38.7475, 39.8458] (0.95 confidence) | 40.1177 ns [39.4597, 40.8371] (0.95 confidence) | 0.979403 |
| column&#95;basis/4x6/full&#95;rank | 36.085 ns [36.0105, 36.157] (0.95 confidence) | 36.2625 ns [36.1819, 36.3421] (0.95 confidence) | 0.995106 |
| column&#95;basis/4x6/large | 35.9689 ns [35.894, 36.0416] (0.95 confidence) | 36.081 ns [36.0216, 36.1417] (0.95 confidence) | 0.996891 |
| column&#95;basis/4x6/near&#95;dependent | 36.0344 ns [35.9809, 36.0865] (0.95 confidence) | 36.4621 ns [36.3943, 36.534] (0.95 confidence) | 0.988271 |
| column&#95;basis/4x6/rank&#95;one | 20.7442 ns [20.7185, 20.7696] (0.95 confidence) | 20.8474 ns [20.8082, 20.8882] (0.95 confidence) | 0.995047 |
| column&#95;basis/4x6/small | 35.9703 ns [35.9263, 36.0148] (0.95 confidence) | 36.3011 ns [36.2041, 36.4244] (0.95 confidence) | 0.990889 |
| column&#95;basis/5x2/full&#95;rank | 13.1694 ns [13.1447, 13.1936] (0.95 confidence) | 13.0633 ns [13.0382, 13.0874] (0.95 confidence) | 1.00812 |
| column&#95;basis/5x2/large | 13.2034 ns [13.1674, 13.2426] (0.95 confidence) | 13.4528 ns [13.3967, 13.5102] (0.95 confidence) | 0.981456 |
| column&#95;basis/5x2/near&#95;dependent | 13.1424 ns [13.1215, 13.1635] (0.95 confidence) | 13.4016 ns [13.3416, 13.4663] (0.95 confidence) | 0.980663 |
| column&#95;basis/5x2/rank&#95;one | 10.7887 ns [10.7695, 10.8089] (0.95 confidence) | 10.8183 ns [10.7928, 10.8433] (0.95 confidence) | 0.99726 |
| column&#95;basis/5x2/small | 13.1216 ns [13.0958, 13.1476] (0.95 confidence) | 13.4264 ns [13.3759, 13.481] (0.95 confidence) | 0.977301 |
| column&#95;basis/5x5/full&#95;rank | 52.5872 ns [51.9226, 53.2966] (0.95 confidence) | 52.1985 ns [51.5571, 52.8858] (0.95 confidence) | 1.00745 |
| column&#95;basis/5x5/large | 52.769 ns [52.0922, 53.5005] (0.95 confidence) | 52.0027 ns [51.3586, 52.7007] (0.95 confidence) | 1.01474 |
| column&#95;basis/5x5/near&#95;dependent | 52.3087 ns [51.6981, 52.968] (0.95 confidence) | 51.6005 ns [51.0191, 52.2305] (0.95 confidence) | 1.01372 |
| column&#95;basis/5x5/rank&#95;one | 23.8217 ns [23.7701, 23.8739] (0.95 confidence) | 23.5659 ns [23.5412, 23.5889] (0.95 confidence) | 1.01085 |
| column&#95;basis/5x5/small | 52.4758 ns [51.8014, 53.204] (0.95 confidence) | 52.0799 ns [51.4504, 52.7672] (0.95 confidence) | 1.0076 |
| column&#95;basis/5x6/shared&#95;face | 60.8743 ns [60.4889, 61.372] (0.95 confidence) | 62.9388 ns [62.7606, 63.1213] (0.95 confidence) | 0.967197 |
| column&#95;basis/5x8/full&#95;rank | 61.7964 ns [61.63, 61.9267] (0.95 confidence) | 62.4461 ns [62.2623, 62.667] (0.95 confidence) | 0.989595 |
| column&#95;basis/5x8/large | 61.739 ns [61.6477, 61.8339] (0.95 confidence) | 62.6662 ns [62.5125, 62.8294] (0.95 confidence) | 0.985204 |
| column&#95;basis/5x8/near&#95;dependent | 61.4838 ns [61.3889, 61.5816] (0.95 confidence) | 62.3579 ns [62.2407, 62.4808] (0.95 confidence) | 0.985982 |
| column&#95;basis/5x8/rank&#95;one | 37.2181 ns [36.7865, 37.6707] (0.95 confidence) | 37.7746 ns [37.402, 38.164] (0.95 confidence) | 0.985266 |
| column&#95;basis/5x8/small | 61.552 ns [61.4825, 61.6199] (0.95 confidence) | 62.2116 ns [62.1026, 62.3217] (0.95 confidence) | 0.989397 |
| column&#95;basis/6x10/full&#95;rank | 96.0565 ns [95.9288, 96.1835] (0.95 confidence) | 96.9912 ns [96.8325, 97.1518] (0.95 confidence) | 0.990364 |
| column&#95;basis/6x10/large | 96.1839 ns [96.0284, 96.3458] (0.95 confidence) | 98.0617 ns [97.8366, 98.2928] (0.95 confidence) | 0.980851 |
| column&#95;basis/6x10/near&#95;dependent | 96.1643 ns [96.0765, 96.2501] (0.95 confidence) | 97.7876 ns [97.5723, 98.0249] (0.95 confidence) | 0.983399 |
| column&#95;basis/6x10/rank&#95;one | 53.925 ns [53.7632, 54.1045] (0.95 confidence) | 58.7271 ns [58.266, 59.2029] (0.95 confidence) | 0.91823 |
| column&#95;basis/6x10/small | 96.1567 ns [96.0305, 96.2841] (0.95 confidence) | 98.0299 ns [97.7665, 98.3266] (0.95 confidence) | 0.980891 |
| column&#95;basis/6x7/shared&#95;face | 96.6532 ns [96.4011, 96.9061] (0.95 confidence) | 98.6426 ns [98.2462, 99.0516] (0.95 confidence) | 0.979831 |
| column&#95;basis/7x12/full&#95;rank | 140.707 ns [140.484, 140.922] (0.95 confidence) | 143.395 ns [143.118, 143.694] (0.95 confidence) | 0.981252 |
| column&#95;basis/7x12/large | 140.582 ns [140.443, 140.714] (0.95 confidence) | 145.715 ns [145.147, 146.328] (0.95 confidence) | 0.964771 |
| column&#95;basis/7x12/near&#95;dependent | 140.793 ns [140.642, 140.941] (0.95 confidence) | 144.116 ns [143.758, 144.496] (0.95 confidence) | 0.97694 |
| column&#95;basis/7x12/rank&#95;one | 78.7261 ns [78.6093, 78.8375] (0.95 confidence) | 79.4069 ns [79.2218, 79.6319] (0.95 confidence) | 0.991426 |
| column&#95;basis/7x12/small | 140.757 ns [140.585, 140.921] (0.95 confidence) | 146.867 ns [146.503, 147.223] (0.95 confidence) | 0.958398 |
| column&#95;basis/8x14/full&#95;rank | 196.772 ns [196.219, 197.461] (0.95 confidence) | 203.039 ns [202.594, 203.522] (0.95 confidence) | 0.969133 |
| column&#95;basis/8x14/large | 195.928 ns [195.69, 196.159] (0.95 confidence) | 199.719 ns [199.441, 200.033] (0.95 confidence) | 0.981021 |
| column&#95;basis/8x14/near&#95;dependent | 195.57 ns [195.344, 195.789] (0.95 confidence) | 201.405 ns [200.988, 201.852] (0.95 confidence) | 0.971028 |
| column&#95;basis/8x14/rank&#95;one | 104.168 ns [103.996, 104.331] (0.95 confidence) | 106.453 ns [106.203, 106.716] (0.95 confidence) | 0.978534 |
| column&#95;basis/8x14/small | 195.743 ns [195.434, 196.042] (0.95 confidence) | 199.71 ns [199.425, 200.024] (0.95 confidence) | 0.980132 |

### Median (ns)

| Workload | Baseline | Current | Baseline / current |
| --- | --- | --- | --- |
| column&#95;basis/2x2/full&#95;rank | 4.95972 ns [4.93853, 4.97863] (0.95 confidence) | 4.828 ns [4.81965, 4.84477] (0.95 confidence) | 1.02728 |
| column&#95;basis/2x2/large | 4.91984 ns [4.90302, 4.93503] (0.95 confidence) | 4.87367 ns [4.86572, 4.88439] (0.95 confidence) | 1.00947 |
| column&#95;basis/2x2/near&#95;dependent | 4.95989 ns [4.94601, 4.97383] (0.95 confidence) | 4.90007 ns [4.88242, 4.90703] (0.95 confidence) | 1.01221 |
| column&#95;basis/2x2/rank&#95;one | 4.72453 ns [4.71839, 4.74314] (0.95 confidence) | 4.66211 ns [4.64763, 4.67584] (0.95 confidence) | 1.01339 |
| column&#95;basis/2x2/small | 4.93176 ns [4.91481, 4.94781] (0.95 confidence) | 4.84097 ns [4.8376, 4.84868] (0.95 confidence) | 1.01875 |
| column&#95;basis/2x3/shared&#95;face | 9.92427 ns [9.8661, 10.011] (0.95 confidence) | 10.2223 ns [10.1985, 10.2592] (0.95 confidence) | 0.970842 |
| column&#95;basis/3x3/full&#95;rank | 13.4359 ns [13.4009, 13.4463] (0.95 confidence) | 13.2823 ns [13.2655, 13.317] (0.95 confidence) | 1.01156 |
| column&#95;basis/3x3/large | 13.4073 ns [13.3922, 13.4274] (0.95 confidence) | 13.3262 ns [13.2935, 13.3672] (0.95 confidence) | 1.00609 |
| column&#95;basis/3x3/near&#95;dependent | 13.4109 ns [13.3999, 13.4384] (0.95 confidence) | 13.3493 ns [13.33, 13.3901] (0.95 confidence) | 1.00462 |
| column&#95;basis/3x3/rank&#95;one | 9.21872 ns [9.19132, 9.25387] (0.95 confidence) | 9.22944 ns [9.19742, 9.25089] (0.95 confidence) | 0.998838 |
| column&#95;basis/3x3/small | 13.3908 ns [13.3779, 13.4071] (0.95 confidence) | 13.3981 ns [13.3348, 13.8797] (0.95 confidence) | 0.999457 |
| column&#95;basis/3x4/shared&#95;face | 17.2899 ns [17.259, 17.3092] (0.95 confidence) | 17.9827 ns [17.9255, 18.025] (0.95 confidence) | 0.961476 |
| column&#95;basis/4x4/full&#95;rank | 29.8678 ns [29.7804, 29.9401] (0.95 confidence) | 29.456 ns [29.3135, 29.5504] (0.95 confidence) | 1.01398 |
| column&#95;basis/4x4/large | 29.8693 ns [29.6764, 29.9901] (0.95 confidence) | 29.6029 ns [29.4237, 29.7062] (0.95 confidence) | 1.009 |
| column&#95;basis/4x4/near&#95;dependent | 30.0937 ns [29.9171, 30.1605] (0.95 confidence) | 29.674 ns [29.4125, 29.7567] (0.95 confidence) | 1.01414 |
| column&#95;basis/4x4/rank&#95;one | 14.5872 ns [14.5686, 14.6074] (0.95 confidence) | 14.5335 ns [14.5246, 14.5513] (0.95 confidence) | 1.0037 |
| column&#95;basis/4x4/small | 30.1426 ns [30.0732, 30.3138] (0.95 confidence) | 29.5351 ns [29.4466, 29.6498] (0.95 confidence) | 1.02057 |
| column&#95;basis/4x5/shared&#95;face | 39.203 ns [37.1602, 40.0315] (0.95 confidence) | 38.6027 ns [38.3401, 39.2796] (0.95 confidence) | 1.01555 |
| column&#95;basis/4x6/full&#95;rank | 36.0887 ns [36.0366, 36.1908] (0.95 confidence) | 36.3179 ns [36.2309, 36.402] (0.95 confidence) | 0.993689 |
| column&#95;basis/4x6/large | 35.962 ns [35.9102, 36.0061] (0.95 confidence) | 36.0589 ns [36.0067, 36.1255] (0.95 confidence) | 0.997313 |
| column&#95;basis/4x6/near&#95;dependent | 36.0354 ns [36.0033, 36.067] (0.95 confidence) | 36.4271 ns [36.3696, 36.4943] (0.95 confidence) | 0.989247 |
| column&#95;basis/4x6/rank&#95;one | 20.7476 ns [20.728, 20.7717] (0.95 confidence) | 20.803 ns [20.7673, 20.8387] (0.95 confidence) | 0.997336 |
| column&#95;basis/4x6/small | 35.9573 ns [35.9182, 36.0021] (0.95 confidence) | 36.1986 ns [36.1536, 36.3103] (0.95 confidence) | 0.993332 |
| column&#95;basis/5x2/full&#95;rank | 13.1847 ns [13.1434, 13.2086] (0.95 confidence) | 13.0909 ns [13.0683, 13.1091] (0.95 confidence) | 1.00717 |
| column&#95;basis/5x2/large | 13.1727 ns [13.1311, 13.199] (0.95 confidence) | 13.4099 ns [13.3385, 13.4664] (0.95 confidence) | 0.982305 |
| column&#95;basis/5x2/near&#95;dependent | 13.132 ns [13.1034, 13.1452] (0.95 confidence) | 13.3438 ns [13.3036, 13.387] (0.95 confidence) | 0.984128 |
| column&#95;basis/5x2/rank&#95;one | 10.7709 ns [10.7507, 10.8037] (0.95 confidence) | 10.8218 ns [10.7992, 10.8543] (0.95 confidence) | 0.995297 |
| column&#95;basis/5x2/small | 13.1291 ns [13.0954, 13.1458] (0.95 confidence) | 13.3518 ns [13.3159, 13.4071] (0.95 confidence) | 0.983313 |
| column&#95;basis/5x5/full&#95;rank | 50.9585 ns [50.8514, 51.1216] (0.95 confidence) | 50.6185 ns [50.5683, 50.8305] (0.95 confidence) | 1.00672 |
| column&#95;basis/5x5/large | 51.165 ns [51.1234, 51.2459] (0.95 confidence) | 50.6252 ns [50.547, 50.7127] (0.95 confidence) | 1.01066 |
| column&#95;basis/5x5/near&#95;dependent | 51.1714 ns [50.9183, 51.2728] (0.95 confidence) | 50.5924 ns [50.2168, 50.8168] (0.95 confidence) | 1.01144 |
| column&#95;basis/5x5/rank&#95;one | 23.7696 ns [23.7246, 23.8589] (0.95 confidence) | 23.5956 ns [23.5625, 23.6211] (0.95 confidence) | 1.00738 |
| column&#95;basis/5x5/small | 50.8553 ns [50.7915, 50.9804] (0.95 confidence) | 50.6678 ns [50.5617, 50.8088] (0.95 confidence) | 1.0037 |
| column&#95;basis/5x6/shared&#95;face | 60.4123 ns [60.3215, 60.5546] (0.95 confidence) | 62.9381 ns [62.5966, 63.1168] (0.95 confidence) | 0.959869 |
| column&#95;basis/5x8/full&#95;rank | 61.8813 ns [61.813, 61.9866] (0.95 confidence) | 62.3144 ns [62.1512, 62.4456] (0.95 confidence) | 0.99305 |
| column&#95;basis/5x8/large | 61.6304 ns [61.5548, 61.7825] (0.95 confidence) | 62.4919 ns [62.2962, 62.6273] (0.95 confidence) | 0.986213 |
| column&#95;basis/5x8/near&#95;dependent | 61.47 ns [61.4077, 61.515] (0.95 confidence) | 62.2266 ns [62.153, 62.3935] (0.95 confidence) | 0.987841 |
| column&#95;basis/5x8/rank&#95;one | 36.5204 ns [35.7479, 37.613] (0.95 confidence) | 37.3592 ns [36.7562, 37.7487] (0.95 confidence) | 0.977546 |
| column&#95;basis/5x8/small | 61.58 ns [61.5306, 61.6685] (0.95 confidence) | 62.1971 ns [62.0765, 62.3133] (0.95 confidence) | 0.990079 |
| column&#95;basis/6x10/full&#95;rank | 96.0995 ns [95.9504, 96.2873] (0.95 confidence) | 96.8954 ns [96.6494, 97.0958] (0.95 confidence) | 0.991787 |
| column&#95;basis/6x10/large | 96.1287 ns [95.9586, 96.3412] (0.95 confidence) | 98.0269 ns [97.8452, 98.1698] (0.95 confidence) | 0.980636 |
| column&#95;basis/6x10/near&#95;dependent | 96.2388 ns [96.1023, 96.315] (0.95 confidence) | 97.3912 ns [97.2862, 97.5652] (0.95 confidence) | 0.988167 |
| column&#95;basis/6x10/rank&#95;one | 53.7379 ns [53.623, 53.8915] (0.95 confidence) | 58.6324 ns [58.0538, 58.9741] (0.95 confidence) | 0.916522 |
| column&#95;basis/6x10/small | 96.1424 ns [96.025, 96.211] (0.95 confidence) | 97.4988 ns [97.3491, 97.9195] (0.95 confidence) | 0.986088 |
| column&#95;basis/6x7/shared&#95;face | 96.6206 ns [96.3175, 96.8951] (0.95 confidence) | 98.7602 ns [97.9653, 99.2936] (0.95 confidence) | 0.978336 |
| column&#95;basis/7x12/full&#95;rank | 140.81 ns [140.609, 140.989] (0.95 confidence) | 142.885 ns [142.713, 143.108] (0.95 confidence) | 0.985475 |
| column&#95;basis/7x12/large | 140.668 ns [140.576, 140.799] (0.95 confidence) | 145.586 ns [145.019, 146.149] (0.95 confidence) | 0.966221 |
| column&#95;basis/7x12/near&#95;dependent | 140.836 ns [140.665, 140.95] (0.95 confidence) | 143.923 ns [143.535, 144.302] (0.95 confidence) | 0.978555 |
| column&#95;basis/7x12/rank&#95;one | 78.7799 ns [78.6954, 78.9366] (0.95 confidence) | 79.1241 ns [79.027, 79.2376] (0.95 confidence) | 0.99565 |
| column&#95;basis/7x12/small | 140.83 ns [140.627, 140.998] (0.95 confidence) | 147.649 ns [146.696, 147.954] (0.95 confidence) | 0.953811 |
| column&#95;basis/8x14/full&#95;rank | 196.058 ns [195.837, 196.404] (0.95 confidence) | 202.992 ns [202.571, 203.152] (0.95 confidence) | 0.96584 |
| column&#95;basis/8x14/large | 195.887 ns [195.721, 196.169] (0.95 confidence) | 199.377 ns [199.218, 199.569] (0.95 confidence) | 0.982497 |
| column&#95;basis/8x14/near&#95;dependent | 195.656 ns [195.536, 195.787] (0.95 confidence) | 200.998 ns [200.489, 201.656] (0.95 confidence) | 0.973422 |
| column&#95;basis/8x14/rank&#95;one | 104.38 ns [104.089, 104.526] (0.95 confidence) | 105.893 ns [105.62, 106.305] (0.95 confidence) | 0.98571 |
| column&#95;basis/8x14/small | 195.736 ns [195.491, 195.929] (0.95 confidence) | 199.342 ns [199.073, 199.643] (0.95 confidence) | 0.981914 |

## interval-solve

Semantic contract: la-stack-interval-solve-v1. Collection session: `08822c214acf491d9fca8e6bbed69f09`.
Collection interval (Unix ns): 1794919997752899000–1794920759549830000.
Acquisition: newly collected.

- baseline: available — inside declared release interval; release v0.5.3; revision `7165ed562f2e4611d201e574fc1f600a8ad8b33e`.
  Source SHA-256: `073502f95ae58344aa776b54303ae24d1e6cf2ce60d1437d950bcbe64d0b4dbf`; harness SHA-256:
  `ba63ffff14e6124d5aca8df516e6836a606ad10809bbe3b1ab13dd5fd7ceedbf`; scientific gate: passed.
- current: available — inside declared release interval; release v0.5.4; revision `531f75eb58eb8ac810ad5d89bf1d3e66a301e84e`.
  Source SHA-256: `50d72eb42299ad56e00e6838b34b10824c560a3eb5ae9c9a05e48859fd3cfe46`; harness SHA-256:
  `ba63ffff14e6124d5aca8df516e6836a606ad10809bbe3b1ab13dd5fd7ceedbf`; scientific gate: passed.

### Mean (ns)

| Workload | Baseline | Current | Baseline / current |
| --- | --- | --- | --- |
| interval&#95;division/negative | 5.31922 ns [5.29621, 5.34287] (0.95 confidence) | 5.23637 ns [5.22784, 5.24491] (0.95 confidence) | 1.01582 |
| interval&#95;division/positive | 5.35988 ns [5.34828, 5.372] (0.95 confidence) | 5.34889 ns [5.33967, 5.35826] (0.95 confidence) | 1.00205 |
| interval&#95;division/reciprocal | 10.517 ns [10.4941, 10.5408] (0.95 confidence) | 10.1222 ns [10.0998, 10.145] (0.95 confidence) | 1.03901 |
| interval&#95;division/subnormal | 4.15534 ns [4.13616, 4.17697] (0.95 confidence) | 3.99778 ns [3.98699, 4.00844] (0.95 confidence) | 1.03941 |
| interval&#95;division/tiny&#95;ratio | 6.53195 ns [6.51935, 6.54575] (0.95 confidence) | 6.26725 ns [6.25691, 6.27763] (0.95 confidence) | 1.04224 |
| interval&#95;solve&#95;d1/dense | 241.149 ns [239.877, 242.995] (0.95 confidence) | 232.699 ns [231.598, 234.579] (0.95 confidence) | 1.03631 |
| interval&#95;solve&#95;d1/near&#95;singular | 312.193 ns [310.583, 313.852] (0.95 confidence) | 296.229 ns [295.819, 296.673] (0.95 confidence) | 1.05389 |
| interval&#95;solve&#95;d1/range&#95;exhausted | 164.925 ns [164.315, 165.599] (0.95 confidence) | 157.089 ns [156.602, 157.646] (0.95 confidence) | 1.04988 |
| interval&#95;solve&#95;d1/singular | 54.9183 ns [54.6623, 55.186] (0.95 confidence) | 53.3387 ns [53.0759, 53.6227] (0.95 confidence) | 1.02962 |
| interval&#95;solve&#95;d1/uncertain&#95;dense | 396.769 ns [395.324, 398.338] (0.95 confidence) | 381.14 ns [380.396, 381.988] (0.95 confidence) | 1.04101 |
| interval&#95;solve&#95;d2/dense | 3231.12 ns [3217.64, 3244.96] (0.95 confidence) | 3058.59 ns [3051.35, 3065.68] (0.95 confidence) | 1.05641 |
| interval&#95;solve&#95;d2/near&#95;singular | 1768.19 ns [1757.24, 1779.1] (0.95 confidence) | 1663.68 ns [1660.65, 1666.86] (0.95 confidence) | 1.06282 |
| interval&#95;solve&#95;d2/range&#95;exhausted | 268.535 ns [267.8, 269.261] (0.95 confidence) | 255.059 ns [254.914, 255.234] (0.95 confidence) | 1.05283 |
| interval&#95;solve&#95;d2/singular | 178.774 ns [178.283, 179.231] (0.95 confidence) | 173.421 ns [173.061, 173.876] (0.95 confidence) | 1.03087 |
| interval&#95;solve&#95;d2/uncertain&#95;dense | 4421.08 ns [4409.7, 4431.96] (0.95 confidence) | 4224.34 ns [4217, 4231.69] (0.95 confidence) | 1.04657 |
| interval&#95;solve&#95;d3/dense | 8044.21 ns [8022.65, 8067.5] (0.95 confidence) | 7717.5 ns [7706.25, 7728.32] (0.95 confidence) | 1.04233 |
| interval&#95;solve&#95;d3/near&#95;singular | 2050.98 ns [2048.06, 2054.06] (0.95 confidence) | 1986.8 ns [1984.32, 1989.24] (0.95 confidence) | 1.0323 |
| interval&#95;solve&#95;d3/range&#95;exhausted | 429.042 ns [428.13, 429.977] (0.95 confidence) | 413.157 ns [412.843, 413.501] (0.95 confidence) | 1.03845 |
| interval&#95;solve&#95;d3/singular | 342.504 ns [341.551, 343.53] (0.95 confidence) | 328.024 ns [327.63, 328.502] (0.95 confidence) | 1.04414 |
| interval&#95;solve&#95;d3/uncertain&#95;dense | 9773.27 ns [9757.75, 9790.55] (0.95 confidence) | 9421.36 ns [9407.31, 9435.35] (0.95 confidence) | 1.03735 |
| interval&#95;solve&#95;d4/dense | 15616.1 ns [15573.7, 15659.4] (0.95 confidence) | 14859.8 ns [14836.4, 14882.9] (0.95 confidence) | 1.0509 |
| interval&#95;solve&#95;d4/near&#95;singular | 2714.65 ns [2710.74, 2718.74] (0.95 confidence) | 2617.81 ns [2614.67, 2621.13] (0.95 confidence) | 1.03699 |
| interval&#95;solve&#95;d4/range&#95;exhausted | 788.033 ns [787.234, 788.892] (0.95 confidence) | 774.991 ns [772.593, 777.605] (0.95 confidence) | 1.01683 |
| interval&#95;solve&#95;d4/singular | 605.891 ns [604.84, 607.051] (0.95 confidence) | 584.579 ns [583.986, 585.19] (0.95 confidence) | 1.03646 |
| interval&#95;solve&#95;d4/uncertain&#95;dense | 18803.3 ns [18738, 18872.5] (0.95 confidence) | 17579.5 ns [17548.4, 17610.6] (0.95 confidence) | 1.06962 |
| interval&#95;solve&#95;d5/dense | 26736.7 ns [26660.3, 26833.8] (0.95 confidence) | 25844.7 ns [25770.5, 25936] (0.95 confidence) | 1.03451 |
| interval&#95;solve&#95;d5/near&#95;singular | 3578.39 ns [3571.74, 3585.53] (0.95 confidence) | 3487.56 ns [3482.55, 3492.82] (0.95 confidence) | 1.02605 |
| interval&#95;solve&#95;d5/range&#95;exhausted | 1292.28 ns [1290.09, 1294.64] (0.95 confidence) | 1271.09 ns [1269, 1273.5] (0.95 confidence) | 1.01667 |
| interval&#95;solve&#95;d5/singular | 1055.76 ns [1053.32, 1058.25] (0.95 confidence) | 1034.06 ns [1032.52, 1035.66] (0.95 confidence) | 1.02098 |
| interval&#95;solve&#95;d5/uncertain&#95;dense | 30821 ns [30741, 30907.8] (0.95 confidence) | 29514.3 ns [29447.9, 29580.1] (0.95 confidence) | 1.04427 |
| interval&#95;solve&#95;d6/dense | 41534.3 ns [41453.9, 41617.9] (0.95 confidence) | 40611.4 ns [40502.5, 40731.6] (0.95 confidence) | 1.02272 |
| interval&#95;solve&#95;d6/near&#95;singular | 4449.16 ns [4442.73, 4455.59] (0.95 confidence) | 4407.96 ns [4402.29, 4413.94] (0.95 confidence) | 1.00935 |
| interval&#95;solve&#95;d6/range&#95;exhausted | 1620.6 ns [1616.16, 1626.18] (0.95 confidence) | 1612.68 ns [1609.76, 1615.86] (0.95 confidence) | 1.00491 |
| interval&#95;solve&#95;d6/singular | 1360.94 ns [1358.11, 1364.02] (0.95 confidence) | 1339.89 ns [1337.81, 1342.04] (0.95 confidence) | 1.01571 |
| interval&#95;solve&#95;d6/uncertain&#95;dense | 47305.4 ns [47204.5, 47410.1] (0.95 confidence) | 45992.8 ns [45854.5, 46138.8] (0.95 confidence) | 1.02854 |

### Median (ns)

| Workload | Baseline | Current | Baseline / current |
| --- | --- | --- | --- |
| interval&#95;division/negative | 5.26442 ns [5.25076, 5.28368] (0.95 confidence) | 5.23531 ns [5.22491, 5.24522] (0.95 confidence) | 1.00556 |
| interval&#95;division/positive | 5.35254 ns [5.34476, 5.36995] (0.95 confidence) | 5.3505 ns [5.33438, 5.36361] (0.95 confidence) | 1.00038 |
| interval&#95;division/reciprocal | 10.4833 ns [10.459, 10.5169] (0.95 confidence) | 10.1282 ns [10.0906, 10.1652] (0.95 confidence) | 1.03507 |
| interval&#95;division/subnormal | 4.13905 ns [4.12841, 4.15677] (0.95 confidence) | 4.00696 ns [4.00255, 4.01711] (0.95 confidence) | 1.03296 |
| interval&#95;division/tiny&#95;ratio | 6.51612 ns [6.50778, 6.53195] (0.95 confidence) | 6.26887 ns [6.25185, 6.27462] (0.95 confidence) | 1.03944 |
| interval&#95;solve&#95;d1/dense | 239.148 ns [238.814, 239.947] (0.95 confidence) | 231.654 ns [231.541, 231.849] (0.95 confidence) | 1.03235 |
| interval&#95;solve&#95;d1/near&#95;singular | 310.689 ns [307.914, 313.062] (0.95 confidence) | 295.701 ns [295.215, 296.335] (0.95 confidence) | 1.05069 |
| interval&#95;solve&#95;d1/range&#95;exhausted | 164.66 ns [164.451, 164.919] (0.95 confidence) | 156.162 ns [155.877, 156.782] (0.95 confidence) | 1.05441 |
| interval&#95;solve&#95;d1/singular | 54.8088 ns [53.8721, 55.4504] (0.95 confidence) | 52.6921 ns [52.6502, 52.8265] (0.95 confidence) | 1.04017 |
| interval&#95;solve&#95;d1/uncertain&#95;dense | 393.827 ns [392.183, 395.914] (0.95 confidence) | 380.241 ns [379.561, 380.782] (0.95 confidence) | 1.03573 |
| interval&#95;solve&#95;d2/dense | 3234.43 ns [3224.62, 3241.28] (0.95 confidence) | 3060.57 ns [3055.61, 3065.01] (0.95 confidence) | 1.05681 |
| interval&#95;solve&#95;d2/near&#95;singular | 1772.93 ns [1736.22, 1795.62] (0.95 confidence) | 1662.12 ns [1659.27, 1664.68] (0.95 confidence) | 1.06667 |
| interval&#95;solve&#95;d2/range&#95;exhausted | 269.818 ns [268.702, 270.51] (0.95 confidence) | 254.99 ns [254.826, 255.079] (0.95 confidence) | 1.05815 |
| interval&#95;solve&#95;d2/singular | 180.1 ns [179.379, 180.309] (0.95 confidence) | 173.257 ns [173.019, 173.45] (0.95 confidence) | 1.0395 |
| interval&#95;solve&#95;d2/uncertain&#95;dense | 4452.05 ns [4441.48, 4456.86] (0.95 confidence) | 4224.08 ns [4217.66, 4231.41] (0.95 confidence) | 1.05397 |
| interval&#95;solve&#95;d3/dense | 7996.7 ns [7985.34, 8010.45] (0.95 confidence) | 7728.75 ns [7715.63, 7734.9] (0.95 confidence) | 1.03467 |
| interval&#95;solve&#95;d3/near&#95;singular | 2049.21 ns [2045.65, 2052.29] (0.95 confidence) | 1985.88 ns [1983.9, 1987.67] (0.95 confidence) | 1.03189 |
| interval&#95;solve&#95;d3/range&#95;exhausted | 427.95 ns [427.459, 428.984] (0.95 confidence) | 412.899 ns [412.771, 413.036] (0.95 confidence) | 1.03645 |
| interval&#95;solve&#95;d3/singular | 341.131 ns [340.69, 341.645] (0.95 confidence) | 327.764 ns [327.448, 328.034] (0.95 confidence) | 1.04078 |
| interval&#95;solve&#95;d3/uncertain&#95;dense | 9754.45 ns [9744.48, 9763.11] (0.95 confidence) | 9434 ns [9412.83, 9438.25] (0.95 confidence) | 1.03397 |
| interval&#95;solve&#95;d4/dense | 15566.6 ns [15537.1, 15632.5] (0.95 confidence) | 14863.3 ns [14839, 14884.7] (0.95 confidence) | 1.04732 |
| interval&#95;solve&#95;d4/near&#95;singular | 2710.75 ns [2705.14, 2713.69] (0.95 confidence) | 2616.35 ns [2614.14, 2618.14] (0.95 confidence) | 1.03608 |
| interval&#95;solve&#95;d4/range&#95;exhausted | 787.147 ns [786.432, 787.976] (0.95 confidence) | 771.897 ns [770.127, 774.072] (0.95 confidence) | 1.01976 |
| interval&#95;solve&#95;d4/singular | 604.235 ns [603.387, 605.026] (0.95 confidence) | 584.353 ns [584.009, 584.851] (0.95 confidence) | 1.03402 |
| interval&#95;solve&#95;d4/uncertain&#95;dense | 18726.2 ns [18640.6, 18793.3] (0.95 confidence) | 17578.6 ns [17510.1, 17642.3] (0.95 confidence) | 1.06528 |
| interval&#95;solve&#95;d5/dense | 26705.1 ns [26615.3, 26780] (0.95 confidence) | 25757.8 ns [25695.1, 25823] (0.95 confidence) | 1.03678 |
| interval&#95;solve&#95;d5/near&#95;singular | 3567.58 ns [3564.19, 3571.67] (0.95 confidence) | 3479.24 ns [3474.09, 3489.23] (0.95 confidence) | 1.02539 |
| interval&#95;solve&#95;d5/range&#95;exhausted | 1290.56 ns [1288.71, 1292.21] (0.95 confidence) | 1270.09 ns [1269.05, 1271.48] (0.95 confidence) | 1.01612 |
| interval&#95;solve&#95;d5/singular | 1053.15 ns [1051.02, 1055.18] (0.95 confidence) | 1034.89 ns [1032.58, 1036.7] (0.95 confidence) | 1.01764 |
| interval&#95;solve&#95;d5/uncertain&#95;dense | 30722.6 ns [30656.7, 30815.3] (0.95 confidence) | 29443.6 ns [29401.8, 29503.4] (0.95 confidence) | 1.04344 |
| interval&#95;solve&#95;d6/dense | 41492.1 ns [41400, 41615.2] (0.95 confidence) | 40519.5 ns [40445.8, 40637.4] (0.95 confidence) | 1.024 |
| interval&#95;solve&#95;d6/near&#95;singular | 4451.61 ns [4441.41, 4460.45] (0.95 confidence) | 4400.73 ns [4394.58, 4409.01] (0.95 confidence) | 1.01156 |
| interval&#95;solve&#95;d6/range&#95;exhausted | 1614.26 ns [1611.96, 1619.69] (0.95 confidence) | 1610.35 ns [1606.88, 1612.9] (0.95 confidence) | 1.00243 |
| interval&#95;solve&#95;d6/singular | 1358.72 ns [1356.71, 1362.33] (0.95 confidence) | 1340.11 ns [1334.95, 1341.89] (0.95 confidence) | 1.01389 |
| interval&#95;solve&#95;d6/uncertain&#95;dense | 47211.9 ns [47105, 47308.3] (0.95 confidence) | 45875.7 ns [45748.9, 46059.9] (0.95 confidence) | 1.02913 |

## polynomial

Semantic contract: la-stack-polynomial-v1. Collection session: `7c82055c7c4f4bba928b13a6692f41a5`.
Collection interval (Unix ns): 1794920760634717000–1794920961309851000.
Acquisition: newly collected.

- baseline: available — inside declared release interval; release v0.5.3; revision `7165ed562f2e4611d201e574fc1f600a8ad8b33e`.
  Source SHA-256: `073502f95ae58344aa776b54303ae24d1e6cf2ce60d1437d950bcbe64d0b4dbf`; harness SHA-256:
  `ba63ffff14e6124d5aca8df516e6836a606ad10809bbe3b1ab13dd5fd7ceedbf`; scientific gate: passed.
- current: available — inside declared release interval; release v0.5.4; revision `531f75eb58eb8ac810ad5d89bf1d3e66a301e84e`.
  Source SHA-256: `50d72eb42299ad56e00e6838b34b10824c560a3eb5ae9c9a05e48859fd3cfe46`; harness SHA-256:
  `ba63ffff14e6124d5aca8df516e6836a606ad10809bbe3b1ab13dd5fd7ceedbf`; scientific gate: passed.

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
<td>13510.1 ns [13477.9, 13542.7] (0.95 confidence)</td>
<td>13976.4 ns [13932.2, 14025.3] (0.95 confidence)</td>
<td>0.966634</td>
</tr>
<tr>
<td>polynomial/coincident&#95;lifted/order5/determinant</td>
<td>107192 ns [106955, 107430] (0.95 confidence)</td>
<td>109075 ns [108720, 109432] (0.95 confidence)</td>
<td>0.982744</td>
</tr>
<tr>
<td>polynomial/coincident&#95;lifted/order6/determinant</td>
<td>922336 ns [920547, 924134] (0.95 confidence)</td>
<td>942104 ns [939866, 944338] (0.95 confidence)</td>
<td>0.979017</td>
</tr>
<tr>
<td>polynomial/coincident&#95;lifted/order7/determinant</td>
<td>8.63594e+06 ns [8.61915e+06, 8.65341e+06] (0.95 confidence)</td>
<td>8.9215e+06 ns [8.90037e+06, 8.94272e+06] (0.95 confidence)</td>
<td>0.967991</td>
</tr>
<tr>
<td>polynomial/coincident&#95;lifted/order8/determinant</td>
<td>8.79634e+07 ns [8.76769e+07, 8.82866e+07] (0.95 confidence)</td>
<td>8.85879e+07 ns [8.82445e+07, 8.89316e+07] (0.95 confidence)</td>
<td>0.99295</td>
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
<td>13500.9 ns [13476.9, 13545] (0.95 confidence)</td>
<td>13978.6 ns [13934.6, 14021.7] (0.95 confidence)</td>
<td>0.96583</td>
</tr>
<tr>
<td>polynomial/coincident&#95;lifted/order5/determinant</td>
<td>107140 ns [106809, 107478] (0.95 confidence)</td>
<td>108899 ns [108504, 109547] (0.95 confidence)</td>
<td>0.983848</td>
</tr>
<tr>
<td>polynomial/coincident&#95;lifted/order6/determinant</td>
<td>921898 ns [919443, 924374] (0.95 confidence)</td>
<td>941792 ns [939335, 944702] (0.95 confidence)</td>
<td>0.978876</td>
</tr>
<tr>
<td>polynomial/coincident&#95;lifted/order7/determinant</td>
<td>8.63632e+06 ns [8.6196e+06, 8.6474e+06] (0.95 confidence)</td>
<td>8.92533e+06 ns [8.90181e+06, 8.94206e+06] (0.95 confidence)</td>
<td>0.96762</td>
</tr>
<tr>
<td>polynomial/coincident&#95;lifted/order8/determinant</td>
<td>8.78232e+07 ns [8.72825e+07, 8.80884e+07] (0.95 confidence)</td>
<td>8.84672e+07 ns [8.79208e+07, 8.90537e+07] (0.95 confidence)</td>
<td>0.99272</td>
</tr>
</tbody>
</table>

## rational-dot

Semantic contract: la-stack-rational-dot-v1. Collection session: `82322b5f71a3444481968f1405ba74b3`.
Collection interval (Unix ns): 1794920962362176000–1794922073158949000.
Acquisition: newly collected.

- baseline: available — inside declared release interval; release v0.5.3; revision `7165ed562f2e4611d201e574fc1f600a8ad8b33e`.
  Source SHA-256: `073502f95ae58344aa776b54303ae24d1e6cf2ce60d1437d950bcbe64d0b4dbf`; harness SHA-256:
  `ba63ffff14e6124d5aca8df516e6836a606ad10809bbe3b1ab13dd5fd7ceedbf`; scientific gate: passed.
- current: available — inside declared release interval; release v0.5.4; revision `531f75eb58eb8ac810ad5d89bf1d3e66a301e84e`.
  Source SHA-256: `50d72eb42299ad56e00e6838b34b10824c560a3eb5ae9c9a05e48859fd3cfe46`; harness SHA-256:
  `ba63ffff14e6124d5aca8df516e6836a606ad10809bbe3b1ab13dd5fd7ceedbf`; scientific gate: passed.

### Mean (ns)

| Workload | Baseline | Current | Baseline / current |
| --- | --- | --- | --- |
| rational&#95;dot/cancellation/d14/prepared | 1822.79 ns [1814.29, 1831.21] (0.95 confidence) | 1889.21 ns [1883.69, 1895.07] (0.95 confidence) | 0.964841 |
| rational&#95;dot/cancellation/d14/slice | 628.65 ns [626.941, 630.312] (0.95 confidence) | 684.511 ns [682.253, 686.797] (0.95 confidence) | 0.918393 |
| rational&#95;dot/cancellation/d2/prepared | 316.34 ns [312.816, 320.65] (0.95 confidence) | 267.308 ns [266.064, 268.545] (0.95 confidence) | 1.18343 |
| rational&#95;dot/cancellation/d2/slice | 103.728 ns [103.265, 104.216] (0.95 confidence) | 98.1341 ns [97.8665, 98.4177] (0.95 confidence) | 1.05701 |
| rational&#95;dot/cancellation/d22/prepared | 2928.17 ns [2920.2, 2936.34] (0.95 confidence) | 2983.6 ns [2972.34, 2994.07] (0.95 confidence) | 0.981422 |
| rational&#95;dot/cancellation/d22/slice | 1001.07 ns [995.691, 1006.81] (0.95 confidence) | 1023.07 ns [1021.36, 1024.82] (0.95 confidence) | 0.978496 |
| rational&#95;dot/cancellation/d3/prepared | 399.94 ns [398.707, 401.184] (0.95 confidence) | 400.227 ns [398.274, 402.673] (0.95 confidence) | 0.999282 |
| rational&#95;dot/cancellation/d3/slice | 135.849 ns [135.083, 136.644] (0.95 confidence) | 129.894 ns [129.571, 130.228] (0.95 confidence) | 1.04584 |
| rational&#95;dot/cancellation/d4/prepared | 555.869 ns [552.437, 559.341] (0.95 confidence) | 536.585 ns [534.185, 539.01] (0.95 confidence) | 1.03594 |
| rational&#95;dot/cancellation/d4/slice | 197.398 ns [196.714, 198.119] (0.95 confidence) | 194.432 ns [193.915, 195.009] (0.95 confidence) | 1.01525 |
| rational&#95;dot/cancellation/d5/prepared | 699.659 ns [690.945, 708.443] (0.95 confidence) | 660.821 ns [658.615, 663.034] (0.95 confidence) | 1.05877 |
| rational&#95;dot/cancellation/d5/slice | 224.16 ns [223.405, 224.91] (0.95 confidence) | 221.969 ns [221.27, 222.72] (0.95 confidence) | 1.00987 |
| rational&#95;dot/cancellation/d6/prepared | 832.578 ns [830.914, 834.324] (0.95 confidence) | 811.406 ns [808.705, 813.968] (0.95 confidence) | 1.02609 |
| rational&#95;dot/cancellation/d6/slice | 284.835 ns [284.034, 285.652] (0.95 confidence) | 278.371 ns [277.46, 279.301] (0.95 confidence) | 1.02322 |
| rational&#95;dot/cancellation/d7/prepared | 928.076 ns [924.522, 931.555] (0.95 confidence) | 935.469 ns [930.411, 940.984] (0.95 confidence) | 0.992097 |
| rational&#95;dot/cancellation/d7/slice | 303.411 ns [302.871, 303.999] (0.95 confidence) | 312.482 ns [311.43, 313.581] (0.95 confidence) | 0.970969 |
| rational&#95;dot/cancellation/d8/prepared | 1058.89 ns [1054.44, 1063.32] (0.95 confidence) | 1071.51 ns [1067.42, 1075.44] (0.95 confidence) | 0.988223 |
| rational&#95;dot/cancellation/d8/slice | 362.628 ns [361.347, 364.061] (0.95 confidence) | 379.587 ns [377.747, 381.433] (0.95 confidence) | 0.955323 |
| rational&#95;dot/large/d14/prepared | 2909.03 ns [2901.21, 2917.56] (0.95 confidence) | 3102.33 ns [3084.62, 3120.26] (0.95 confidence) | 0.937691 |
| rational&#95;dot/large/d14/slice | 1275 ns [1271.31, 1278.85] (0.95 confidence) | 1333.62 ns [1325.24, 1342.57] (0.95 confidence) | 0.95605 |
| rational&#95;dot/large/d2/prepared | 367.006 ns [365.121, 368.93] (0.95 confidence) | 379.906 ns [377.906, 381.889] (0.95 confidence) | 0.966043 |
| rational&#95;dot/large/d2/slice | 138.741 ns [137.944, 139.635] (0.95 confidence) | 135.138 ns [133.992, 136.569] (0.95 confidence) | 1.02666 |
| rational&#95;dot/large/d22/prepared | 4697.02 ns [4686.7, 4707.44] (0.95 confidence) | 4858.72 ns [4845.62, 4872.35] (0.95 confidence) | 0.966719 |
| rational&#95;dot/large/d22/slice | 2083.36 ns [2077.86, 2088.96] (0.95 confidence) | 2157.39 ns [2150.86, 2164.1] (0.95 confidence) | 0.965683 |
| rational&#95;dot/large/d3/prepared | 578.264 ns [576.431, 580.182] (0.95 confidence) | 586.365 ns [582.952, 589.834] (0.95 confidence) | 0.986185 |
| rational&#95;dot/large/d3/slice | 227.507 ns [226.417, 228.644] (0.95 confidence) | 222 ns [220.91, 223.187] (0.95 confidence) | 1.02481 |
| rational&#95;dot/large/d4/prepared | 773.819 ns [771.347, 776.357] (0.95 confidence) | 773.919 ns [771.975, 775.915] (0.95 confidence) | 0.999871 |
| rational&#95;dot/large/d4/slice | 368.509 ns [366.17, 370.767] (0.95 confidence) | 310.123 ns [308.934, 311.325] (0.95 confidence) | 1.18827 |
| rational&#95;dot/large/d5/prepared | 982.395 ns [980.601, 984.17] (0.95 confidence) | 964.451 ns [962.506, 966.421] (0.95 confidence) | 1.01861 |
| rational&#95;dot/large/d5/slice | 449.336 ns [446.213, 453.168] (0.95 confidence) | 396.697 ns [395.23, 398.179] (0.95 confidence) | 1.13269 |
| rational&#95;dot/large/d6/prepared | 1194.99 ns [1191.2, 1198.86] (0.95 confidence) | 1184.82 ns [1181.57, 1188.2] (0.95 confidence) | 1.00859 |
| rational&#95;dot/large/d6/slice | 500.878 ns [498.885, 502.939] (0.95 confidence) | 502.076 ns [500.295, 503.891] (0.95 confidence) | 0.997614 |
| rational&#95;dot/large/d7/prepared | 1364.97 ns [1362.6, 1367.45] (0.95 confidence) | 1407.45 ns [1404.6, 1410.34] (0.95 confidence) | 0.969819 |
| rational&#95;dot/large/d7/slice | 607.714 ns [605.319, 610.186] (0.95 confidence) | 616.202 ns [613.653, 618.768] (0.95 confidence) | 0.986225 |
| rational&#95;dot/large/d8/prepared | 1557.82 ns [1554.21, 1561.63] (0.95 confidence) | 1609.17 ns [1604.46, 1614.17] (0.95 confidence) | 0.968089 |
| rational&#95;dot/large/d8/slice | 675.076 ns [672.289, 677.908] (0.95 confidence) | 674.889 ns [671.871, 677.905] (0.95 confidence) | 1.00028 |
| rational&#95;dot/non&#95;dyadic/d14/prepared | 2399.73 ns [2393.38, 2406.53] (0.95 confidence) | 2540.48 ns [2530.15, 2551.47] (0.95 confidence) | 0.944596 |
| rational&#95;dot/non&#95;dyadic/d14/slice | 1038.35 ns [1035.42, 1041.32] (0.95 confidence) | 1049.61 ns [1044.51, 1054.78] (0.95 confidence) | 0.989276 |
| rational&#95;dot/non&#95;dyadic/d2/prepared | 323.119 ns [321.431, 324.815] (0.95 confidence) | 276.84 ns [276.042, 277.682] (0.95 confidence) | 1.16717 |
| rational&#95;dot/non&#95;dyadic/d2/slice | 105.075 ns [104.038, 106.194] (0.95 confidence) | 92.9724 ns [92.8067, 93.1386] (0.95 confidence) | 1.13017 |
| rational&#95;dot/non&#95;dyadic/d22/prepared | 3852.64 ns [3845.61, 3859.78] (0.95 confidence) | 3993.37 ns [3984.38, 4002.55] (0.95 confidence) | 0.964758 |
| rational&#95;dot/non&#95;dyadic/d22/slice | 1700.03 ns [1695.88, 1704.47] (0.95 confidence) | 1765.6 ns [1760.13, 1771.26] (0.95 confidence) | 0.962862 |
| rational&#95;dot/non&#95;dyadic/d3/prepared | 508.368 ns [506.123, 510.6] (0.95 confidence) | 458.805 ns [457.345, 460.382] (0.95 confidence) | 1.10803 |
| rational&#95;dot/non&#95;dyadic/d3/slice | 188.408 ns [185.904, 190.931] (0.95 confidence) | 177.496 ns [176.726, 178.336] (0.95 confidence) | 1.06148 |
| rational&#95;dot/non&#95;dyadic/d4/prepared | 618.411 ns [615.922, 621.828] (0.95 confidence) | 634.993 ns [631.224, 638.719] (0.95 confidence) | 0.973886 |
| rational&#95;dot/non&#95;dyadic/d4/slice | 242.766 ns [241.698, 243.88] (0.95 confidence) | 232.931 ns [232.238, 233.626] (0.95 confidence) | 1.04222 |
| rational&#95;dot/non&#95;dyadic/d5/prepared | 893.137 ns [886.618, 899.021] (0.95 confidence) | 802.503 ns [800.897, 804.276] (0.95 confidence) | 1.11294 |
| rational&#95;dot/non&#95;dyadic/d5/slice | 320.577 ns [319.213, 322.051] (0.95 confidence) | 310.995 ns [310.102, 311.854] (0.95 confidence) | 1.03081 |
| rational&#95;dot/non&#95;dyadic/d6/prepared | 1057.39 ns [1047.13, 1067.08] (0.95 confidence) | 962.052 ns [959.321, 965.102] (0.95 confidence) | 1.0991 |
| rational&#95;dot/non&#95;dyadic/d6/slice | 399.022 ns [397.856, 400.229] (0.95 confidence) | 392.826 ns [391.809, 393.857] (0.95 confidence) | 1.01577 |
| rational&#95;dot/non&#95;dyadic/d7/prepared | 1145.87 ns [1143.84, 1148.03] (0.95 confidence) | 1154.31 ns [1152.15, 1156.7] (0.95 confidence) | 0.992689 |
| rational&#95;dot/non&#95;dyadic/d7/slice | 493.507 ns [489.036, 498.396] (0.95 confidence) | 490.225 ns [488.908, 491.57] (0.95 confidence) | 1.0067 |
| rational&#95;dot/non&#95;dyadic/d8/prepared | 1289.04 ns [1285.19, 1292.65] (0.95 confidence) | 1296.42 ns [1294.62, 1298.3] (0.95 confidence) | 0.994312 |
| rational&#95;dot/non&#95;dyadic/d8/slice | 544.348 ns [542.163, 546.822] (0.95 confidence) | 549.159 ns [547.549, 550.804] (0.95 confidence) | 0.991239 |

### Median (ns)

| Workload | Baseline | Current | Baseline / current |
| --- | --- | --- | --- |
| rational&#95;dot/cancellation/d14/prepared | 1841.63 ns [1795.71, 1850.68] (0.95 confidence) | 1888.69 ns [1879.67, 1894.73] (0.95 confidence) | 0.975085 |
| rational&#95;dot/cancellation/d14/slice | 628.993 ns [628.032, 630.093] (0.95 confidence) | 683.224 ns [680.899, 685.599] (0.95 confidence) | 0.920625 |
| rational&#95;dot/cancellation/d2/prepared | 314.655 ns [312.185, 316.85] (0.95 confidence) | 269.92 ns [264.892, 271.543] (0.95 confidence) | 1.16574 |
| rational&#95;dot/cancellation/d2/slice | 103.074 ns [102.653, 103.646] (0.95 confidence) | 98.0167 ns [97.6244, 98.3045] (0.95 confidence) | 1.0516 |
| rational&#95;dot/cancellation/d22/prepared | 2927.02 ns [2919.16, 2931.99] (0.95 confidence) | 3002.56 ns [2989.64, 3010.1] (0.95 confidence) | 0.974843 |
| rational&#95;dot/cancellation/d22/slice | 987.854 ns [986.34, 989.373] (0.95 confidence) | 1022.64 ns [1021.09, 1023.95] (0.95 confidence) | 0.965988 |
| rational&#95;dot/cancellation/d3/prepared | 399.474 ns [398.878, 401.163] (0.95 confidence) | 399.314 ns [398.202, 400.085] (0.95 confidence) | 1.0004 |
| rational&#95;dot/cancellation/d3/slice | 134.722 ns [134.347, 135.51] (0.95 confidence) | 129.634 ns [129.235, 130.069] (0.95 confidence) | 1.03925 |
| rational&#95;dot/cancellation/d4/prepared | 556.388 ns [554.236, 560.396] (0.95 confidence) | 537.215 ns [535.093, 538.995] (0.95 confidence) | 1.03569 |
| rational&#95;dot/cancellation/d4/slice | 197.048 ns [196.229, 197.413] (0.95 confidence) | 194.062 ns [193.628, 194.526] (0.95 confidence) | 1.01539 |
| rational&#95;dot/cancellation/d5/prepared | 674.381 ns [667.982, 725.766] (0.95 confidence) | 659.984 ns [655.41, 664.697] (0.95 confidence) | 1.02181 |
| rational&#95;dot/cancellation/d5/slice | 224.481 ns [223.195, 225.375] (0.95 confidence) | 221.944 ns [221.132, 222.42] (0.95 confidence) | 1.01143 |
| rational&#95;dot/cancellation/d6/prepared | 831.997 ns [830.606, 833.696] (0.95 confidence) | 814.703 ns [812.56, 816.593] (0.95 confidence) | 1.02123 |
| rational&#95;dot/cancellation/d6/slice | 284.677 ns [283.532, 285.472] (0.95 confidence) | 277.429 ns [276.684, 278.965] (0.95 confidence) | 1.02612 |
| rational&#95;dot/cancellation/d7/prepared | 930.919 ns [927.678, 934.203] (0.95 confidence) | 929.099 ns [927.641, 931.723] (0.95 confidence) | 1.00196 |
| rational&#95;dot/cancellation/d7/slice | 302.918 ns [302.162, 303.596] (0.95 confidence) | 311.122 ns [310.033, 312.461] (0.95 confidence) | 0.97363 |
| rational&#95;dot/cancellation/d8/prepared | 1062.53 ns [1058.76, 1064.71] (0.95 confidence) | 1077.15 ns [1075.64, 1081.35] (0.95 confidence) | 0.98643 |
| rational&#95;dot/cancellation/d8/slice | 360.692 ns [359.675, 361.164] (0.95 confidence) | 381.254 ns [378.081, 383.758] (0.95 confidence) | 0.946068 |
| rational&#95;dot/large/d14/prepared | 2899.79 ns [2894.09, 2906.65] (0.95 confidence) | 3085.85 ns [3069.09, 3132.94] (0.95 confidence) | 0.939706 |
| rational&#95;dot/large/d14/slice | 1272.4 ns [1267.62, 1277.68] (0.95 confidence) | 1322.38 ns [1314.39, 1333.01] (0.95 confidence) | 0.962202 |
| rational&#95;dot/large/d2/prepared | 368.119 ns [366.192, 369.195] (0.95 confidence) | 380.362 ns [378.436, 382.964] (0.95 confidence) | 0.967815 |
| rational&#95;dot/large/d2/slice | 138.122 ns [137.208, 139.12] (0.95 confidence) | 133.645 ns [132.881, 134.088] (0.95 confidence) | 1.0335 |
| rational&#95;dot/large/d22/prepared | 4696.64 ns [4686.05, 4708.61] (0.95 confidence) | 4843.84 ns [4834.14, 4865.27] (0.95 confidence) | 0.96961 |
| rational&#95;dot/large/d22/slice | 2078.21 ns [2071.16, 2090.5] (0.95 confidence) | 2151.32 ns [2145.3, 2161.89] (0.95 confidence) | 0.966019 |
| rational&#95;dot/large/d3/prepared | 576.9 ns [574.892, 578.581] (0.95 confidence) | 589.824 ns [584.05, 592.97] (0.95 confidence) | 0.978089 |
| rational&#95;dot/large/d3/slice | 225.635 ns [224.899, 226.651] (0.95 confidence) | 221.356 ns [219.956, 222.448] (0.95 confidence) | 1.01933 |
| rational&#95;dot/large/d4/prepared | 772.127 ns [769.006, 775.021] (0.95 confidence) | 772.855 ns [770.543, 775.144] (0.95 confidence) | 0.999058 |
| rational&#95;dot/large/d4/slice | 370.154 ns [367.159, 372.919] (0.95 confidence) | 309.881 ns [307.853, 310.93] (0.95 confidence) | 1.19451 |
| rational&#95;dot/large/d5/prepared | 982.513 ns [980.46, 984.532] (0.95 confidence) | 965.401 ns [961.956, 967.532] (0.95 confidence) | 1.01773 |
| rational&#95;dot/large/d5/slice | 446.974 ns [444.334, 450.267] (0.95 confidence) | 397.214 ns [395.573, 398.455] (0.95 confidence) | 1.12527 |
| rational&#95;dot/large/d6/prepared | 1193.91 ns [1188.71, 1198.7] (0.95 confidence) | 1181.97 ns [1178.59, 1184.79] (0.95 confidence) | 1.0101 |
| rational&#95;dot/large/d6/slice | 499.651 ns [496.989, 502.701] (0.95 confidence) | 501.118 ns [499.684, 503.196] (0.95 confidence) | 0.997073 |
| rational&#95;dot/large/d7/prepared | 1363.11 ns [1361.64, 1366.21] (0.95 confidence) | 1405.97 ns [1403.76, 1409.39] (0.95 confidence) | 0.96952 |
| rational&#95;dot/large/d7/slice | 606.131 ns [604.072, 608.859] (0.95 confidence) | 615.089 ns [611.615, 619.289] (0.95 confidence) | 0.985436 |
| rational&#95;dot/large/d8/prepared | 1553.97 ns [1549.82, 1559.83] (0.95 confidence) | 1605.4 ns [1600.18, 1611.96] (0.95 confidence) | 0.967966 |
| rational&#95;dot/large/d8/slice | 673.458 ns [671.43, 675.708] (0.95 confidence) | 677.326 ns [673.735, 679.847] (0.95 confidence) | 0.994289 |
| rational&#95;dot/non&#95;dyadic/d14/prepared | 2391.83 ns [2386.51, 2396.29] (0.95 confidence) | 2531.83 ns [2525.71, 2537.18] (0.95 confidence) | 0.944702 |
| rational&#95;dot/non&#95;dyadic/d14/slice | 1036.81 ns [1032.27, 1040.22] (0.95 confidence) | 1050.04 ns [1035.31, 1058.26] (0.95 confidence) | 0.987397 |
| rational&#95;dot/non&#95;dyadic/d2/prepared | 323.546 ns [322.059, 325.08] (0.95 confidence) | 275.374 ns [274.733, 276.928] (0.95 confidence) | 1.17493 |
| rational&#95;dot/non&#95;dyadic/d2/slice | 103.985 ns [102.369, 105.352] (0.95 confidence) | 92.8809 ns [92.7316, 93.1052] (0.95 confidence) | 1.11955 |
| rational&#95;dot/non&#95;dyadic/d22/prepared | 3852.32 ns [3843.33, 3860.42] (0.95 confidence) | 3989.44 ns [3978.9, 4001.91] (0.95 confidence) | 0.96563 |
| rational&#95;dot/non&#95;dyadic/d22/slice | 1695.42 ns [1689.54, 1701.15] (0.95 confidence) | 1759.93 ns [1749.17, 1771.1] (0.95 confidence) | 0.963347 |
| rational&#95;dot/non&#95;dyadic/d3/prepared | 505.897 ns [503.492, 511.705] (0.95 confidence) | 458.319 ns [457.271, 459.115] (0.95 confidence) | 1.10381 |
| rational&#95;dot/non&#95;dyadic/d3/slice | 180.031 ns [178.045, 197.903] (0.95 confidence) | 176.288 ns [175.786, 176.638] (0.95 confidence) | 1.02123 |
| rational&#95;dot/non&#95;dyadic/d4/prepared | 617.876 ns [615.623, 619.547] (0.95 confidence) | 636.176 ns [631.665, 640.081] (0.95 confidence) | 0.971233 |
| rational&#95;dot/non&#95;dyadic/d4/slice | 242.104 ns [241.389, 243.187] (0.95 confidence) | 232.919 ns [231.802, 233.898] (0.95 confidence) | 1.03944 |
| rational&#95;dot/non&#95;dyadic/d5/prepared | 893.773 ns [891.941, 896.953] (0.95 confidence) | 800.837 ns [799.742, 802.682] (0.95 confidence) | 1.11605 |
| rational&#95;dot/non&#95;dyadic/d5/slice | 318.866 ns [318.102, 319.768] (0.95 confidence) | 311.989 ns [310.29, 312.996] (0.95 confidence) | 1.02204 |
| rational&#95;dot/non&#95;dyadic/d6/prepared | 1083.59 ns [1080.3, 1085.41] (0.95 confidence) | 958.951 ns [956.412, 960.757] (0.95 confidence) | 1.12997 |
| rational&#95;dot/non&#95;dyadic/d6/slice | 398.655 ns [396.842, 399.828] (0.95 confidence) | 392.726 ns [391.087, 393.876] (0.95 confidence) | 1.0151 |
| rational&#95;dot/non&#95;dyadic/d7/prepared | 1144.76 ns [1143.33, 1147.29] (0.95 confidence) | 1153.6 ns [1151.39, 1155.49] (0.95 confidence) | 0.992339 |
| rational&#95;dot/non&#95;dyadic/d7/slice | 484.099 ns [483.48, 485.861] (0.95 confidence) | 489.462 ns [488.62, 491.125] (0.95 confidence) | 0.989044 |
| rational&#95;dot/non&#95;dyadic/d8/prepared | 1288.84 ns [1286.82, 1293.13] (0.95 confidence) | 1295.96 ns [1294.17, 1297.83] (0.95 confidence) | 0.994503 |
| rational&#95;dot/non&#95;dyadic/d8/slice | 543.675 ns [539.856, 546.19] (0.95 confidence) | 548.85 ns [546.532, 550.193] (0.95 confidence) | 0.990573 |
