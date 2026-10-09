# Double estimator for the generator: pre-registered report

## Outcome table

| Criterion | Verdict |
|---|---|
| D-1 | PASS |
| D-2 | FAIL |
| D-3 | PASS |
| D-4 | PASS |
| D-5 | FAIL |
| D-6 | EXPLORATORY |

## Criteria (verbatim)

### D-1

> For P1 and MR, in cells (3,6) and (4,6), |mean generator bias| of `double` <= 0.05 x that of `raw`, at every d.

**Verdict: PASS**

```json
{
  "details": [
    {
      "cell": [
        3,
        6
      ],
      "dimension": 20,
      "double_abs_bias": 0.734795826938473,
      "passed": true,
      "pde_id": "P1",
      "ratio": 0.024601644990218798,
      "raw_abs_bias": 29.867751820279317
    },
    {
      "cell": [
        4,
        6
      ],
      "dimension": 20,
      "double_abs_bias": 0.5620813779597679,
      "passed": true,
      "pde_id": "P1",
      "ratio": 0.0003340420343652255,
      "raw_abs_bias": 1682.6666111882648
    },
    {
      "cell": [
        3,
        6
      ],
      "dimension": 100,
      "double_abs_bias": 0.6941429132591488,
      "passed": true,
      "pde_id": "P1",
      "ratio": 0.00048041481045315237,
      "raw_abs_bias": 1444.8824186007023
    },
    {
      "cell": [
        4,
        6
      ],
      "dimension": 100,
      "double_abs_bias": 0.5332557068598712,
      "passed": true,
      "pde_id": "P1",
      "ratio": 1.9434662686393364e-08,
      "raw_abs_bias": 27438382.4131517
    },
    {
      "cell": [
        3,
        6
      ],
      "dimension": 400,
      "double_abs_bias": 0.7315134928322118,
      "passed": true,
      "pde_id": "P1",
      "ratio": 9.420654711656088e-06,
      "raw_abs_bias": 77649.96332230678
    },
    {
      "cell": [
        4,
        6
      ],
      "dimension": 400,
      "double_abs_bias": 0.9219871236811116,
      "passed": true,
      "pde_id": "P1",
      "ratio": 2.913705326988244e-12,
      "raw_abs_bias": 316431148730.51526
    },
    {
      "cell": [
        3,
        6
      ],
      "dimension": 100,
      "double_abs_bias": 2.6711293412363393,
      "passed": true,
      "pde_id": "MR",
      "ratio": 0.0002482228824657319,
      "raw_abs_bias": 10761.011695225556
    },
    {
      "cell": [
        4,
        6
      ],
      "dimension": 100,
      "double_abs_bias": 2.1559003260806597,
      "passed": true,
      "pde_id": "MR",
      "ratio": 4.505621874527096e-09,
      "raw_abs_bias": 478491179.71244323
    },
    {
      "cell": [
        3,
        6
      ],
      "dimension": 400,
      "double_abs_bias": 1.563632167093017,
      "passed": true,
      "pde_id": "MR",
      "ratio": 4.740244717446048e-06,
      "raw_abs_bias": 329863.17380159895
    },
    {
      "cell": [
        4,
        6
      ],
      "dimension": 400,
      "double_abs_bias": 8.440904448243902,
      "passed": true,
      "pde_id": "MR",
      "ratio": 3.904637049274607e-12,
      "raw_abs_bias": 2161764164434.1902
    }
  ],
  "missing": [],
  "verdict": "PASS"
}
```

### D-2

> For P1 and MR, `double` skill <= 0.1 x `raw` skill in cells (3,6) and (4,6) at every d, and `double` skill at the largest d <= 1.5 x its value at the smallest d in each of those cells.

**Verdict: FAIL**

```json
{
  "dimension_growth": [
    {
      "cell": [
        3,
        6
      ],
      "growth": 18.027908218893785,
      "largest_d_skill": 8.779866986465283,
      "passed": false,
      "pde_id": "P1",
      "smallest_d_skill": 0.4870152920605465
    },
    {
      "cell": [
        4,
        6
      ],
      "growth": 1699.6820046198939,
      "largest_d_skill": 2526.0780219801068,
      "passed": false,
      "pde_id": "P1",
      "smallest_d_skill": 1.4862062521777555
    },
    {
      "cell": [
        3,
        6
      ],
      "growth": 4.1523109232699245,
      "largest_d_skill": 12.135265005364337,
      "passed": false,
      "pde_id": "MR",
      "smallest_d_skill": 2.922532832827431
    },
    {
      "cell": [
        4,
        6
      ],
      "growth": 38.493725314344054,
      "largest_d_skill": 5997.723024339303,
      "passed": false,
      "pde_id": "MR",
      "smallest_d_skill": 155.81040741994255
    }
  ],
  "missing": [],
  "pointwise": [
    {
      "cell": [
        3,
        6
      ],
      "dimension": 20,
      "double_skill": 0.4870152920605465,
      "passed": true,
      "pde_id": "P1",
      "ratio": 0.008933603139712998,
      "raw_skill": 54.51499069794054
    },
    {
      "cell": [
        3,
        6
      ],
      "dimension": 100,
      "double_skill": 1.6676206608378337,
      "passed": true,
      "pde_id": "P1",
      "ratio": 0.00048199051485358057,
      "raw_skill": 3459.861987832737
    },
    {
      "cell": [
        3,
        6
      ],
      "dimension": 400,
      "double_skill": 8.779866986465283,
      "passed": true,
      "pde_id": "P1",
      "ratio": 4.651833344890065e-05,
      "raw_skill": 188739.9297335054
    },
    {
      "cell": [
        4,
        6
      ],
      "dimension": 20,
      "double_skill": 1.4862062521777555,
      "passed": true,
      "pde_id": "P1",
      "ratio": 1.2100041983114059e-05,
      "raw_skill": 122826.5368212604
    },
    {
      "cell": [
        4,
        6
      ],
      "dimension": 100,
      "double_skill": 43.99962770385012,
      "passed": true,
      "pde_id": "P1",
      "ratio": 2.8627570507704718e-08,
      "raw_skill": 1536966879.2539775
    },
    {
      "cell": [
        4,
        6
      ],
      "dimension": 400,
      "double_skill": 2526.0780219801068,
      "passed": true,
      "pde_id": "P1",
      "ratio": 1.8543260449918965e-10,
      "raw_skill": 13622620621667.135
    },
    {
      "cell": [
        3,
        6
      ],
      "dimension": 100,
      "double_skill": 2.922532832827431,
      "passed": true,
      "pde_id": "MR",
      "ratio": 0.00039738660870432355,
      "raw_skill": 7354.3817753606
    },
    {
      "cell": [
        3,
        6
      ],
      "dimension": 400,
      "double_skill": 12.135265005364337,
      "passed": true,
      "pde_id": "MR",
      "ratio": 4.4652973217055364e-05,
      "raw_skill": 271768.3533944214
    },
    {
      "cell": [
        4,
        6
      ],
      "dimension": 100,
      "double_skill": 155.81040741994255,
      "passed": true,
      "pde_id": "MR",
      "ratio": 2.6476686599077898e-08,
      "raw_skill": 5884815187.764807
    },
    {
      "cell": [
        4,
        6
      ],
      "dimension": 400,
      "double_skill": 5997.723024339303,
      "passed": true,
      "pde_id": "MR",
      "ratio": 1.6594625396404287e-10,
      "raw_skill": 36142563517214.95
    }
  ],
  "verdict": "FAIL"
}
```

### D-3

> `double_path` skill <= `path` skill in cells (3,6) and (4,6) for at least 75% of the (PDE, d, cell) combinations over P1 and MR.

**Verdict: PASS**

```json
{
  "combinations": 10,
  "details": [
    {
      "cell": [
        3,
        6
      ],
      "dimension": 20,
      "double_path_skill": 0.20372474706016663,
      "passed": true,
      "path_skill": 1.4469982713368428,
      "pde_id": "P1"
    },
    {
      "cell": [
        4,
        6
      ],
      "dimension": 20,
      "double_path_skill": 0.38644459025424377,
      "passed": true,
      "path_skill": 44.573736364787074,
      "pde_id": "P1"
    },
    {
      "cell": [
        3,
        6
      ],
      "dimension": 100,
      "double_path_skill": 0.3684825673978145,
      "passed": true,
      "path_skill": 7.292648491104707,
      "pde_id": "P1"
    },
    {
      "cell": [
        4,
        6
      ],
      "dimension": 100,
      "double_path_skill": 2.678214487568524,
      "passed": true,
      "path_skill": 4984.437754020386,
      "pde_id": "P1"
    },
    {
      "cell": [
        3,
        6
      ],
      "dimension": 400,
      "double_path_skill": 0.6853662415589847,
      "passed": true,
      "path_skill": 28.458059467038872,
      "pde_id": "P1"
    },
    {
      "cell": [
        4,
        6
      ],
      "dimension": 400,
      "double_path_skill": 16.47648455207607,
      "passed": true,
      "path_skill": 258411.6351287691,
      "pde_id": "P1"
    },
    {
      "cell": [
        3,
        6
      ],
      "dimension": 100,
      "double_path_skill": 0.5888786564239403,
      "passed": true,
      "path_skill": 12.294436187285765,
      "pde_id": "MR"
    },
    {
      "cell": [
        4,
        6
      ],
      "dimension": 100,
      "double_path_skill": 7.758222029983775,
      "passed": true,
      "path_skill": 15434.367533334389,
      "pde_id": "MR"
    },
    {
      "cell": [
        3,
        6
      ],
      "dimension": 400,
      "double_path_skill": 0.8031552870668769,
      "passed": true,
      "path_skill": 34.159467342086984,
      "pde_id": "MR"
    },
    {
      "cell": [
        4,
        6
      ],
      "dimension": 400,
      "double_path_skill": 30.609082037764146,
      "passed": true,
      "path_skill": 532705.410914412,
      "pde_id": "MR"
    }
  ],
  "missing": [],
  "success_fraction": 1.0,
  "successes": 10,
  "verdict": "PASS"
}
```

### D-4

> For P1 and MR at d=100 and d=400, the frontier of the better of `double`/`double_path` is <= 0.8 x the `raw` frontier at >= 8 of 10 cost levels.

**Verdict: PASS**

```json
{
  "problems": {
    "MR_d100": {
      "levels": [
        {
          "cost_level": 9600.0,
          "double_method": "double_path",
          "double_skill": 0.21056949019643256,
          "level_index": 0,
          "passed": true,
          "ratio": 0.02759233993885674,
          "raw_skill": 7.631447374997704
        },
        {
          "cost_level": 16515.455598736746,
          "double_method": "double_path",
          "double_skill": 0.21056949019643256,
          "level_index": 1,
          "passed": true,
          "ratio": 0.02759233993885674,
          "raw_skill": 7.631447374997704
        },
        {
          "cost_level": 28412.528503525547,
          "double_method": "double_path",
          "double_skill": 0.17633691206600838,
          "level_index": 2,
          "passed": true,
          "ratio": 0.04845754603347272,
          "raw_skill": 3.6389979786471485
        },
        {
          "cost_level": 48879.77634873113,
          "double_method": "double_path",
          "double_skill": 0.160043253058667,
          "level_index": 3,
          "passed": true,
          "ratio": 0.09294686101988604,
          "raw_skill": 1.721879053284281
        },
        {
          "cost_level": 84090.80999621362,
          "double_method": "double_path",
          "double_skill": 0.15236811102251915,
          "level_index": 4,
          "passed": true,
          "ratio": 0.19257335211354365,
          "raw_skill": 0.7912211598865507
        },
        {
          "cost_level": 144666.46237023704,
          "double_method": "double_path",
          "double_skill": 0.15236811102251915,
          "level_index": 5,
          "passed": true,
          "ratio": 0.19257335211354365,
          "raw_skill": 0.7912211598865507
        },
        {
          "cost_level": 248878.38915645552,
          "double_method": "double_path",
          "double_skill": 0.1489682815882403,
          "level_index": 6,
          "passed": true,
          "ratio": 0.45075293017467405,
          "raw_skill": 0.3304876610131246
        },
        {
          "cost_level": 428160.4151665178,
          "double_method": "double_path",
          "double_skill": 0.1489682815882403,
          "level_index": 7,
          "passed": true,
          "ratio": 0.45075293017467405,
          "raw_skill": 0.3304876610131246
        },
        {
          "cost_level": 736590.0339395113,
          "double_method": "double_path",
          "double_skill": 0.1489682815882403,
          "level_index": 8,
          "passed": true,
          "ratio": 0.45075293017467405,
          "raw_skill": 0.3304876610131246
        },
        {
          "cost_level": 1267200.0,
          "double_method": "double_path",
          "double_skill": 0.1489682815882403,
          "level_index": 9,
          "passed": true,
          "ratio": 0.45075293017467405,
          "raw_skill": 0.3304876610131246
        }
      ],
      "verdict": "PASS",
      "wins": 10
    },
    "MR_d400": {
      "levels": [
        {
          "cost_level": 9600.0,
          "double_method": "double_path",
          "double_skill": 0.1719990410223505,
          "level_index": 0,
          "passed": true,
          "ratio": 0.00734470161076538,
          "raw_skill": 23.41811146830603
        },
        {
          "cost_level": 16515.455598736746,
          "double_method": "double_path",
          "double_skill": 0.1719990410223505,
          "level_index": 1,
          "passed": true,
          "ratio": 0.00734470161076538,
          "raw_skill": 23.41811146830603
        },
        {
          "cost_level": 28412.528503525547,
          "double_method": "double_path",
          "double_skill": 0.1392817796988881,
          "level_index": 2,
          "passed": true,
          "ratio": 0.01231459157666459,
          "raw_skill": 11.31030443289883
        },
        {
          "cost_level": 48879.77634873113,
          "double_method": "double_path",
          "double_skill": 0.12458599676297406,
          "level_index": 3,
          "passed": true,
          "ratio": 0.022537245267852155,
          "raw_skill": 5.528004655506301
        },
        {
          "cost_level": 84090.80999621362,
          "double_method": "double_path",
          "double_skill": 0.11822314028696498,
          "level_index": 4,
          "passed": true,
          "ratio": 0.04367940738195474,
          "raw_skill": 2.706610445814026
        },
        {
          "cost_level": 144666.46237023704,
          "double_method": "double_path",
          "double_skill": 0.11822314028696498,
          "level_index": 5,
          "passed": true,
          "ratio": 0.04367940738195474,
          "raw_skill": 2.706610445814026
        },
        {
          "cost_level": 248878.38915645552,
          "double_method": "double_path",
          "double_skill": 0.11558920165536424,
          "level_index": 6,
          "passed": true,
          "ratio": 0.08923500895923175,
          "raw_skill": 1.2953346786592779
        },
        {
          "cost_level": 428160.4151665178,
          "double_method": "double_path",
          "double_skill": 0.11558920165536424,
          "level_index": 7,
          "passed": true,
          "ratio": 0.08923500895923175,
          "raw_skill": 1.2953346786592779
        },
        {
          "cost_level": 736590.0339395113,
          "double_method": "double_path",
          "double_skill": 0.11558920165536424,
          "level_index": 8,
          "passed": true,
          "ratio": 0.08923500895923175,
          "raw_skill": 1.2953346786592779
        },
        {
          "cost_level": 1267200.0,
          "double_method": "double_path",
          "double_skill": 0.11558920165536424,
          "level_index": 9,
          "passed": true,
          "ratio": 0.08923500895923175,
          "raw_skill": 1.2953346786592779
        }
      ],
      "verdict": "PASS",
      "wins": 10
    },
    "P1_d100": {
      "levels": [
        {
          "cost_level": 9600.0,
          "double_method": "double_path",
          "double_skill": 0.1807553900257024,
          "level_index": 0,
          "passed": true,
          "ratio": 0.03423163566429853,
          "raw_skill": 5.280360886004026
        },
        {
          "cost_level": 16515.455598736746,
          "double_method": "double_path",
          "double_skill": 0.1807553900257024,
          "level_index": 1,
          "passed": true,
          "ratio": 0.03423163566429853,
          "raw_skill": 5.280360886004026
        },
        {
          "cost_level": 28412.528503525547,
          "double_method": "double_path",
          "double_skill": 0.146807652596202,
          "level_index": 2,
          "passed": true,
          "ratio": 0.05879091021538902,
          "raw_skill": 2.4971148100675915
        },
        {
          "cost_level": 48879.77634873113,
          "double_method": "double_path",
          "double_skill": 0.13338670556397264,
          "level_index": 3,
          "passed": true,
          "ratio": 0.11373927091355689,
          "raw_skill": 1.1727409934370692
        },
        {
          "cost_level": 84090.80999621362,
          "double_method": "double_path",
          "double_skill": 0.12741883783979713,
          "level_index": 4,
          "passed": true,
          "ratio": 0.23952091332504843,
          "raw_skill": 0.5319737473899822
        },
        {
          "cost_level": 144666.46237023704,
          "double_method": "double_path",
          "double_skill": 0.12741883783979713,
          "level_index": 5,
          "passed": true,
          "ratio": 0.23952091332504843,
          "raw_skill": 0.5319737473899822
        },
        {
          "cost_level": 248878.38915645552,
          "double_method": "double_path",
          "double_skill": 0.12424554394186468,
          "level_index": 6,
          "passed": true,
          "ratio": 0.5787682491195558,
          "raw_skill": 0.21467235656218478
        },
        {
          "cost_level": 428160.4151665178,
          "double_method": "double_path",
          "double_skill": 0.12424554394186468,
          "level_index": 7,
          "passed": true,
          "ratio": 0.5787682491195558,
          "raw_skill": 0.21467235656218478
        },
        {
          "cost_level": 736590.0339395113,
          "double_method": "double_path",
          "double_skill": 0.12424554394186468,
          "level_index": 8,
          "passed": true,
          "ratio": 0.5787682491195558,
          "raw_skill": 0.21467235656218478
        },
        {
          "cost_level": 1267200.0,
          "double_method": "double_path",
          "double_skill": 0.09385827602969858,
          "level_index": 9,
          "passed": true,
          "ratio": 0.43721640518960053,
          "raw_skill": 0.21467235656218478
        }
      ],
      "verdict": "PASS",
      "wins": 10
    },
    "P1_d400": {
      "levels": [
        {
          "cost_level": 9600.0,
          "double_method": "double_path",
          "double_skill": 0.18107807235689,
          "level_index": 0,
          "passed": true,
          "ratio": 0.008629168580808524,
          "raw_skill": 20.98441705723677
        },
        {
          "cost_level": 16515.455598736746,
          "double_method": "double_path",
          "double_skill": 0.18107807235689,
          "level_index": 1,
          "passed": true,
          "ratio": 0.008629168580808524,
          "raw_skill": 20.98441705723677
        },
        {
          "cost_level": 28412.528503525547,
          "double_method": "double_path",
          "double_skill": 0.14932226407720436,
          "level_index": 2,
          "passed": true,
          "ratio": 0.014715044370663562,
          "raw_skill": 10.147591832947413
        },
        {
          "cost_level": 48879.77634873113,
          "double_method": "double_path",
          "double_skill": 0.13509961554205857,
          "level_index": 3,
          "passed": true,
          "ratio": 0.027253900383477112,
          "raw_skill": 4.9570745339615225
        },
        {
          "cost_level": 84090.80999621362,
          "double_method": "double_path",
          "double_skill": 0.1283135691706779,
          "level_index": 4,
          "passed": true,
          "ratio": 0.053120198804536785,
          "raw_skill": 2.4155325480393186
        },
        {
          "cost_level": 144666.46237023704,
          "double_method": "double_path",
          "double_skill": 0.1283135691706779,
          "level_index": 5,
          "passed": true,
          "ratio": 0.053120198804536785,
          "raw_skill": 2.4155325480393186
        },
        {
          "cost_level": 248878.38915645552,
          "double_method": "double_path",
          "double_skill": 0.12710402376837965,
          "level_index": 6,
          "passed": true,
          "ratio": 0.11111809780507702,
          "raw_skill": 1.1438642874479825
        },
        {
          "cost_level": 428160.4151665178,
          "double_method": "double_path",
          "double_skill": 0.12710402376837965,
          "level_index": 7,
          "passed": true,
          "ratio": 0.11111809780507702,
          "raw_skill": 1.1438642874479825
        },
        {
          "cost_level": 736590.0339395113,
          "double_method": "double_path",
          "double_skill": 0.12710402376837965,
          "level_index": 8,
          "passed": true,
          "ratio": 0.11111809780507702,
          "raw_skill": 1.1438642874479825
        },
        {
          "cost_level": 1267200.0,
          "double_method": "double_path",
          "double_skill": 0.12710402376837965,
          "level_index": 9,
          "passed": true,
          "ratio": 0.11111809780507702,
          "raw_skill": 1.1438642874479825
        }
      ],
      "verdict": "PASS",
      "wins": 10
    }
  },
  "verdict": "PASS"
}
```

### D-5

> On C2, |bias(`double`)| <= 0.1 x |bias(`raw`)| for convex and flip, at cells (2,32), (3,6), every d.

**Verdict: FAIL**

```json
{
  "details": [
    {
      "cell": [
        2,
        32
      ],
      "dimension": 20,
      "double_abs_bias": 0.9838558984127184,
      "passed": false,
      "pde_id": "C2-convex",
      "ratio": 0.8620256891343859,
      "raw_abs_bias": 1.1413301376211535
    },
    {
      "cell": [
        3,
        6
      ],
      "dimension": 20,
      "double_abs_bias": 0.7255714128301431,
      "passed": true,
      "pde_id": "C2-convex",
      "ratio": 0.02267345340397788,
      "raw_abs_bias": 32.0009219549611
    },
    {
      "cell": [
        2,
        32
      ],
      "dimension": 100,
      "double_abs_bias": 0.9942403563905462,
      "passed": false,
      "pde_id": "C2-convex",
      "ratio": 0.10483893307391703,
      "raw_abs_bias": 9.48350319141033
    },
    {
      "cell": [
        3,
        6
      ],
      "dimension": 100,
      "double_abs_bias": 0.7452452860626414,
      "passed": true,
      "pde_id": "C2-convex",
      "ratio": 0.0005018377985380861,
      "raw_abs_bias": 1485.0321921418247
    },
    {
      "cell": [
        2,
        32
      ],
      "dimension": 400,
      "double_abs_bias": 0.9790400521329218,
      "passed": true,
      "pde_id": "C2-convex",
      "ratio": 0.025396176998368748,
      "raw_abs_bias": 38.55068627832479
    },
    {
      "cell": [
        3,
        6
      ],
      "dimension": 400,
      "double_abs_bias": 0.8353490281615326,
      "passed": true,
      "pde_id": "C2-convex",
      "ratio": 1.148128032153557e-05,
      "raw_abs_bias": 72757.48041746343
    },
    {
      "cell": [
        2,
        32
      ],
      "dimension": 20,
      "double_abs_bias": 0.9708989857347724,
      "passed": false,
      "pde_id": "C2-flip",
      "ratio": 0.4708437071723035,
      "raw_abs_bias": 2.062040908575795
    },
    {
      "cell": [
        3,
        6
      ],
      "dimension": 20,
      "double_abs_bias": 0.7170443509735736,
      "passed": true,
      "pde_id": "C2-flip",
      "ratio": 0.07558651225753876,
      "raw_abs_bias": 9.486406100210793
    },
    {
      "cell": [
        2,
        32
      ],
      "dimension": 100,
      "double_abs_bias": 0.9895142529420452,
      "passed": false,
      "pde_id": "C2-flip",
      "ratio": 0.16176200117958053,
      "raw_abs_bias": 6.117099477791037
    },
    {
      "cell": [
        3,
        6
      ],
      "dimension": 100,
      "double_abs_bias": 0.7250485152576411,
      "passed": true,
      "pde_id": "C2-flip",
      "ratio": 0.004408345474105214,
      "raw_abs_bias": 164.47180002488537
    },
    {
      "cell": [
        2,
        32
      ],
      "dimension": 400,
      "double_abs_bias": 0.9774395295072591,
      "passed": true,
      "pde_id": "C2-flip",
      "ratio": 0.04715675671601598,
      "raw_abs_bias": 20.727454506541342
    },
    {
      "cell": [
        3,
        6
      ],
      "dimension": 400,
      "double_abs_bias": 0.7353198623740705,
      "passed": true,
      "pde_id": "C2-flip",
      "ratio": 9.060579441702606e-05,
      "raw_abs_bias": 8115.594229985516
    }
  ],
  "missing": [],
  "verdict": "FAIL"
}
```

### D-6

> Exploratory P4: report bias and skill of the Double-Q norm form; no verdict.

**Verdict: EXPLORATORY**

```json
{
  "details": [
    {
      "cell": [
        2,
        32
      ],
      "dimension": 20,
      "mean_generator_bias": 0.08581670567509078,
      "mean_generator_rmse": 0.19675007975890038,
      "method": "double",
      "skill": 0.22072611261009945
    },
    {
      "cell": [
        2,
        32
      ],
      "dimension": 20,
      "mean_generator_bias": -0.05385404176753765,
      "mean_generator_rmse": 0.12139108695653937,
      "method": "double_path",
      "skill": 0.20495611481892623
    },
    {
      "cell": [
        3,
        6
      ],
      "dimension": 20,
      "mean_generator_bias": 0.25637191076463517,
      "mean_generator_rmse": 0.45040375461799814,
      "method": "double",
      "skill": 0.4917462594025236
    },
    {
      "cell": [
        3,
        6
      ],
      "dimension": 20,
      "mean_generator_bias": -0.01762235911421178,
      "mean_generator_rmse": 0.21999506049881665,
      "method": "double_path",
      "skill": 0.4551829298691546
    },
    {
      "cell": [
        3,
        10
      ],
      "dimension": 20,
      "mean_generator_bias": 0.20795434073788138,
      "mean_generator_rmse": 0.3640594780530711,
      "method": "double",
      "skill": 0.2952666788671946
    },
    {
      "cell": [
        3,
        10
      ],
      "dimension": 20,
      "mean_generator_bias": -0.03030328054513568,
      "mean_generator_rmse": 0.17198257434579747,
      "method": "double_path",
      "skill": 0.29538634114375995
    },
    {
      "cell": [
        4,
        6
      ],
      "dimension": 20,
      "mean_generator_bias": 0.27777828229909785,
      "mean_generator_rmse": 0.46610831249742934,
      "method": "double",
      "skill": 0.32161370542884965
    },
    {
      "cell": [
        4,
        6
      ],
      "dimension": 20,
      "mean_generator_bias": -0.01428354455491409,
      "mean_generator_rmse": 0.19621749930263324,
      "method": "double_path",
      "skill": 0.27155273452580786
    },
    {
      "cell": [
        2,
        32
      ],
      "dimension": 100,
      "mean_generator_bias": 0.23105182773997118,
      "mean_generator_rmse": 0.30580588201331116,
      "method": "double",
      "skill": 0.41875818846917784
    },
    {
      "cell": [
        2,
        32
      ],
      "dimension": 100,
      "mean_generator_bias": -0.05257804923096051,
      "mean_generator_rmse": 0.11910536148057076,
      "method": "double_path",
      "skill": 0.19257684515746754
    },
    {
      "cell": [
        3,
        6
      ],
      "dimension": 100,
      "mean_generator_bias": 0.3780636319665788,
      "mean_generator_rmse": 0.5379923483529396,
      "method": "double",
      "skill": 0.5640975498652687
    },
    {
      "cell": [
        3,
        6
      ],
      "dimension": 100,
      "mean_generator_bias": -0.012046945938171863,
      "mean_generator_rmse": 0.22182781267648735,
      "method": "double_path",
      "skill": 0.5615759534313547
    },
    {
      "cell": [
        3,
        10
      ],
      "dimension": 100,
      "mean_generator_bias": 0.34817143396601624,
      "mean_generator_rmse": 0.4732083157890877,
      "method": "double",
      "skill": 0.37737121147510005
    },
    {
      "cell": [
        3,
        10
      ],
      "dimension": 100,
      "mean_generator_bias": -0.026582145362309385,
      "mean_generator_rmse": 0.1724752746946127,
      "method": "double_path",
      "skill": 0.4413420905662056
    },
    {
      "cell": [
        4,
        6
      ],
      "dimension": 100,
      "mean_generator_bias": 0.408625283761433,
      "mean_generator_rmse": 0.5627680465047954,
      "method": "double",
      "skill": 0.4052556848007434
    },
    {
      "cell": [
        4,
        6
      ],
      "dimension": 100,
      "mean_generator_bias": -0.008238029277677347,
      "mean_generator_rmse": 0.1996334010407419,
      "method": "double_path",
      "skill": 0.39874331959636106
    }
  ],
  "missing": [],
  "verdict": "EXPLORATORY"
}
```

## Bias--variance separation

The table reports the requested generator RMSE and across-repetition skill spread. Missing rows mean that the corresponding priority block was not completed.

| PDE | d | cell | method | mean bias | generator RMSE | mean skill | skill SD |
|---|---:|---|---|---:|---:|---:|---:|
| MR | 100 | (3,6) | double | 2.67113 | 53.8997 | 2.92253 | 0.200756 |
| MR | 100 | (3,6) | double_path | 2.63486 | 9.20301 | 0.588879 | 0.0582625 |
| MR | 100 | (3,6) | path | -14.0062 | 102.115 | 12.2944 | 0.349825 |
| MR | 100 | (3,6) | raw | -10761 | 70452.3 | 7354.38 | 235.723 |
| MR | 100 | (4,6) | double | 2.1559 | 655.706 | 155.81 | 39.7815 |
| MR | 100 | (4,6) | double_path | 2.08444 | 32.6274 | 7.75822 | 2.2484 |
| MR | 100 | (4,6) | path | -1064.5 | 52397.8 | 15434.4 | 984.334 |
| MR | 100 | (4,6) | raw | -4.78491e+08 | 2.59065e+10 | 5.88482e+09 | 3.64734e+08 |
| MR | 400 | (3,6) | double | 1.56363 | 173.622 | 12.1353 | 1.39243 |
| MR | 400 | (3,6) | double_path | 1.80092 | 11.0665 | 0.803155 | 0.0726128 |
| MR | 400 | (3,6) | path | -36.3762 | 269.275 | 34.1595 | 0.523662 |
| MR | 400 | (3,6) | raw | -329863 | 2.43685e+06 | 271768 | 12746 |
| MR | 400 | (4,6) | double | 8.4409 | 27185.6 | 5997.72 | 1984.77 |
| MR | 400 | (4,6) | double_path | 1.4485 | 119.832 | 30.6091 | 8.53823 |
| MR | 400 | (4,6) | path | -28160.9 | 1.65235e+06 | 532705 | 43108.9 |
| MR | 400 | (4,6) | raw | -2.16176e+12 | 1.23727e+14 | 3.61426e+13 | 6.55349e+12 |
| P1 | 20 | (3,6) | double | 0.734796 | 4.69245 | 0.487015 | 0.0368496 |
| P1 | 20 | (3,6) | double_path | 0.736603 | 1.40763 | 0.203725 | 0.0123077 |
| P1 | 20 | (3,6) | path | 0.207598 | 4.14317 | 1.447 | 0.0571138 |
| P1 | 20 | (3,6) | raw | -29.8678 | 178.933 | 54.515 | 2.59966 |
| P1 | 20 | (4,6) | double | 0.562081 | 4.86473 | 1.48621 | 0.223174 |
| P1 | 20 | (4,6) | double_path | 0.56251 | 1.2556 | 0.386445 | 0.0276538 |
| P1 | 20 | (4,6) | path | -0.801183 | 46.3075 | 44.5737 | 5.04947 |
| P1 | 20 | (4,6) | raw | -1682.67 | 166800 | 122827 | 97931.9 |
| P1 | 100 | (3,6) | double | 0.694143 | 10.276 | 1.66762 | 0.137668 |
| P1 | 100 | (3,6) | double_path | 0.718266 | 1.96858 | 0.368483 | 0.0320701 |
| P1 | 100 | (3,6) | path | -2.19624 | 19.2903 | 7.29265 | 0.0974455 |
| P1 | 100 | (3,6) | raw | -1444.88 | 10584.4 | 3459.86 | 88.9328 |
| P1 | 100 | (4,6) | double | 0.533256 | 60.3982 | 43.9996 | 11.8358 |
| P1 | 100 | (4,6) | double_path | 0.549127 | 3.79241 | 2.67821 | 0.709692 |
| P1 | 100 | (4,6) | path | -94.6636 | 4891.1 | 4984.44 | 366.209 |
| P1 | 100 | (4,6) | raw | -2.74384e+07 | 1.85776e+09 | 1.53697e+09 | 1.66291e+08 |
| P1 | 400 | (3,6) | double | 0.731513 | 43.5709 | 8.77987 | 1.23511 |
| P1 | 400 | (3,6) | double_path | 0.73322 | 3.3755 | 0.685366 | 0.0365857 |
| P1 | 400 | (3,6) | path | -10.9741 | 76.9948 | 28.4581 | 0.807163 |
| P1 | 400 | (3,6) | raw | -77650 | 576030 | 188740 | 5182.24 |
| P1 | 400 | (4,6) | double | -0.921987 | 3547.07 | 2526.08 | 496.3 |
| P1 | 400 | (4,6) | double_path | 0.567843 | 21.55 | 16.4765 | 2.34916 |
| P1 | 400 | (4,6) | path | -5814.73 | 272073 | 258412 | 17801.2 |
| P1 | 400 | (4,6) | raw | -3.16431e+11 | 1.7464e+13 | 1.36226e+13 | 9.76101e+11 |

## Failures and deviations

- No failed criterion is suppressed. Criteria without every required cell and all 10 repetitions are marked `NOT EVALUATED (COMPUTE)`.
- Criteria currently not evaluated: none.
- The two recursive estimates are used separately only inside the generator. Their arithmetic average is the state used for diagnostics and for every non-generator use, as pre-registered.
- The primary RNG uses the exact registered six-part SeedSequence. Auxiliary trees use a separately seeded recursive spawn tree, so their extra draws cannot shift the primary stream.
- For P4, `lambda` in the registered Double-Q formula is implemented as the equation's effective z-coordinate coefficient `lambda_f/sigma`, because this codebase stores `z=sigma*grad(u)`.
- The MR `centre` control retains the established ambient box-hull centre from `effdim_mlp.py`; the certified reference is `sub_box`.

## Provenance and totals

- Frozen preregistration commit: `ab7337e7b7aec9357f0b1746933c846fb0ed5836`.
- Analysis code commit: `1b977ef411776da391e57cfe62d6f63830baaaf6`.
- Result-producing code commits: `1b977ef411776da391e57cfe62d6f63830baaaf6`.
- Completed PDEs: C2-cancel, C2-convex, C2-flip, MR, P1, P4.
- Completed rows: 6080 / 6080 unique registered tasks.
- Aggregate worker wall time: 280956 seconds.
- Non-finite test-state values: 0.
- Non-finite generator values: 0.
- Test points: 1,200 fixed points per PDE/dimension with a fixed 20% validation split; verdicts use the 960-point test subset.
- Arithmetic: float64 throughout.

## Figures

- `results/double_estimator/figures/skill_vs_d_per_cell.png`
- `results/double_estimator/figures/generator_bias_vs_d.png`
- `results/double_estimator/figures/equal_cost_frontiers.png`

## Verification

- The focused new-and-dependent test suite passed: 22 tests passed (`double_estimator`, `mechanism_suite`, and `expert_iteration`).
- Artifact audit: 6,080 unique task rows, no duplicates, exactly 10 repetitions per registered cell, no missing tasks, and no non-finite state or generator values.
- RNG audit: all 960 registered paired groups share the required primary-draw fingerprints for `raw`/`double` and `path`/`double_path`; no primary/auxiliary fingerprint collision was found.
- Test-point audit: every PDE/dimension group uses one fixed test-point fingerprint, as registered.
- The three generated figures were visually inspected for clipping, unreadable labels, and malformed log axes; no rendering issue was found.
- A broader legacy test run also encounters one existing Round-2 mismatch in `test_projection_uses_interval_without_inward_margin`. Both the implementation and that test are unchanged from the base commit, so no out-of-scope repair was attempted.
