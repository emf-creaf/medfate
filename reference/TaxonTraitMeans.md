# Parameter average values

Internal data set with parameter averages for taxonomic families. This
is used by input initialization functions to provide suitable parameter
values when missing from species parameter tables.

## Format

Data frame `trait_family_means` has taxonomic families in rows and
parameter names as columns.

## Source

Same sources as
[`SpParamsMED`](https://emf-creaf.github.io/medfate/reference/SpParams.md)

## See also

[`SpParamsMED`](https://emf-creaf.github.io/medfate/reference/SpParams.md),
[`spwbInput`](https://emf-creaf.github.io/medfate/reference/modelInput.md)

## Examples

``` r
medfate::trait_family_means
#>                    WoodDensity LeafDensity     WoodC       Ptlp    LeafPI0
#> Acanthaceae         0.62347293  0.29727615 0.4108250 -2.7400000 -3.3950000
#> Achariaceae         0.61001599  0.36224725        NA -2.4200000         NA
#> Achatocarpaceae     0.87000000          NA        NA         NA         NA
#> Acoraceae           0.21506871  0.10000000        NA         NA         NA
#> Actinidiaceae       0.37145115  0.45142686        NA -0.7897536         NA
#> Aextoxicaceae       0.56666667          NA        NA         NA         NA
#> Aizoaceae           0.13097800  0.08241249        NA         NA         NA
#> Akaniaceae          0.57133685          NA        NA         NA         NA
#> Alismataceae        0.10037916  0.17703114        NA         NA         NA
#> Alstroemeriaceae            NA          NA        NA         NA         NA
#> Altingiaceae        0.53151501  0.68339710 0.4435000 -1.2487962         NA
#> Alzateaceae         0.66180271          NA        NA         NA         NA
#> Amaranthaceae       0.68154337  0.17128218 0.4199317 -1.5450000 -1.2679576
#> Amaryllidaceae      0.12438751  0.12749970 0.4508400         NA         NA
#> Amborellaceae               NA          NA        NA         NA         NA
#> Anacardiaceae       0.53134629  0.51158932 0.4580588 -2.8013368 -2.4428297
#> Anemiaceae                  NA  0.44850000        NA -1.6450000 -1.3700000
#> Anisophylleaceae    0.57468312          NA        NA         NA         NA
#> Annonaceae          0.55383625  0.41099949 0.4683333 -1.9433333 -2.5250000
#> Aphloiaceae         0.62052000  0.46270602        NA         NA         NA
#> Apiaceae            0.56064111  0.23171368 0.4382558 -1.5040000 -1.1040000
#> Apocynaceae         0.54159387  0.32835059 0.4928286 -0.7800000 -0.6000000
#> Aquifoliaceae       0.58119398  0.43610051 0.4400000 -1.9144202 -2.3000000
#> Araceae             0.06372500  0.17849835        NA         NA         NA
#> Araliaceae          0.41671464  0.31681165 0.4517143 -1.5407671 -0.7673038
#> Araucariaceae       0.42218480  0.38589002        NA         NA         NA
#> Arecaceae           0.55204229  0.42568643 0.4556000 -3.8100000 -3.4000000
#> Aristolochiaceae    0.28046343  0.26814343        NA         NA         NA
#> Asparagaceae        0.65496705  0.18761443 0.4631000         NA -1.7922006
#> Asphodelaceae       0.41942475  0.55721844        NA         NA -1.3794753
#> Aspleniaceae        0.30000000  0.34159513        NA -1.3591176 -1.2400000
#> Asteliaceae                 NA          NA        NA         NA         NA
#> Asteraceae          0.52236825  0.22487138 0.4291509 -1.5967910 -1.4119700
#> Asteropeiaceae      0.81717121          NA        NA         NA         NA
#> Atherospermataceae  0.53917542  0.24008723        NA -1.7600000 -1.3400000
#> Athyriaceae                 NA  0.13005000        NA         NA         NA
#> Aulacomniaceae              NA          NA        NA         NA         NA
#> Austrobaileyaceae   0.36334000  0.26209917        NA         NA         NA
#> Balanopaceae        0.68270079          NA        NA         NA         NA
#> Balsaminaceae       0.16107319  0.38966452        NA         NA         NA
#> Basellaceae         0.20000000          NA        NA         NA         NA
#> Begoniaceae         0.05056471  0.16267894        NA         NA         NA
#> Berberidaceae       0.67673862  0.39692512        NA -1.8233333 -1.7100000
#> Berberidopsidaceae  0.33000000          NA        NA         NA         NA
#> Betulaceae          0.55730771  0.49528026 0.4750059 -1.9169583 -1.5015389
#> Bignoniaceae        0.52135729  0.35958354 0.4692962 -1.8134743 -1.9900000
#> Bixaceae            0.40709211  0.48683360        NA -1.4500000         NA
#> Blechnaceae         0.37060000  0.31479929        NA -1.3400000         NA
#> Bonnetiaceae        0.78381722  0.27978289        NA         NA         NA
#> Boraginaceae        0.47932589  0.37598065 0.4536333 -1.6333333 -1.2177273
#> Brassicaceae        0.46620780  0.24471597 0.4244838 -1.8700000 -1.4800000
#> Bromeliaceae        0.22321054  0.17755239        NA         NA         NA
#> Brunelliaceae       0.43945805          NA        NA         NA         NA
#> Bruniaceae          0.56436667          NA        NA         NA         NA
#> Bryaceae                    NA          NA        NA         NA         NA
#> Burseraceae         0.54281394  0.44904656 0.4703312 -1.8108379 -1.4350000
#> Butomaceae          0.10083260          NA        NA         NA         NA
#> Buxaceae            0.73163747  0.27062945        NA -2.4059098  2.3219014
#> Cabombaceae                 NA          NA        NA         NA         NA
#> Cactaceae           0.60300000  0.15317264        NA         NA         NA
#> Calomniaceae                NA          NA        NA         NA         NA
#> Calophyllaceae      0.60755955  0.39397650 0.4729500         NA         NA
#>                      LeafEPS    LeafAF       Gswmin     Gswmax     Nleaf
#> Acanthaceae        23.230000 0.1270000           NA 0.06000000 25.836769
#> Achariaceae               NA        NA           NA         NA 21.945455
#> Achatocarpaceae           NA        NA           NA         NA        NA
#> Acoraceae                 NA        NA           NA         NA 18.000000
#> Actinidiaceae             NA        NA           NA         NA 21.635986
#> Aextoxicaceae             NA        NA           NA         NA 10.121739
#> Aizoaceae                 NA        NA           NA         NA 14.800000
#> Akaniaceae                NA        NA           NA         NA        NA
#> Alismataceae              NA        NA           NA         NA 31.388709
#> Alstroemeriaceae          NA        NA           NA         NA 18.126667
#> Altingiaceae              NA        NA 0.0031009001 0.53939633 16.490558
#> Alzateaceae               NA        NA           NA         NA  9.802997
#> Amaranthaceae      15.476547 0.6930000 0.0113780532 0.92710526 25.940024
#> Amaryllidaceae            NA        NA 0.0022841333         NA 27.680470
#> Amborellaceae             NA        NA           NA         NA        NA
#> Anacardiaceae      12.760000        NA 0.0036981857 0.06956795 18.853621
#> Anemiaceae         20.040000 0.2550000           NA         NA 21.440000
#> Anisophylleaceae          NA        NA           NA         NA        NA
#> Annonaceae         34.215000        NA 0.0053616476 0.19400000 23.501445
#> Aphloiaceae               NA        NA           NA         NA        NA
#> Apiaceae           23.410000 0.8380000 0.0080862140         NA 26.373905
#> Apocynaceae        14.970000 0.7930000 0.0030300000 1.49093750 21.079302
#> Aquifoliaceae      20.700000 0.4000000 0.0027715622 0.04925700 14.636919
#> Araceae                   NA        NA 0.0042445099         NA 19.245495
#> Araliaceae         12.394746 0.8730000 0.0023960466 1.33950000 18.090127
#> Araucariaceae             NA        NA 0.0018785800         NA 17.428378
#> Arecaceae          73.400000 0.2000000 0.0008583501         NA 19.231886
#> Aristolochiaceae          NA        NA 0.0035895556         NA 31.770806
#> Asparagaceae       22.810356        NA 0.0010977161         NA 24.786812
#> Asphodelaceae       7.425622        NA 0.0055481213         NA 13.623111
#> Aspleniaceae       35.300000        NA 0.0015831553 0.12950000 25.551582
#> Asteliaceae               NA        NA           NA         NA  7.205741
#> Asteraceae          9.586856 0.5153333 0.0113591271 0.32161918 22.667037
#> Asteropeiaceae            NA        NA           NA         NA        NA
#> Atherospermataceae  8.380000        NA           NA 0.55458333 15.392333
#> Athyriaceae               NA        NA           NA         NA 28.101444
#> Aulacomniaceae            NA        NA           NA         NA  8.000000
#> Austrobaileyaceae         NA        NA           NA         NA        NA
#> Balanopaceae              NA        NA           NA         NA        NA
#> Balsaminaceae             NA        NA 0.0126573016         NA 38.195871
#> Basellaceae               NA        NA           NA         NA        NA
#> Begoniaceae               NA        NA           NA         NA 28.277827
#> Berberidaceae             NA        NA 0.0016500000         NA 21.510712
#> Berberidopsidaceae        NA        NA           NA         NA        NA
#> Betulaceae          3.306000        NA 0.0027260202 0.13051410 24.889296
#> Bignoniaceae       17.610000 0.1770000 0.0036919354 0.17275000 25.851301
#> Bixaceae                  NA        NA 0.0027600000 0.44000000 19.032367
#> Blechnaceae               NA        NA           NA         NA 11.745388
#> Bonnetiaceae              NA        NA           NA         NA 10.808790
#> Boraginaceae       12.128750        NA 0.0060134744 0.18000000 25.567842
#> Brassicaceae              NA        NA 0.0119382443         NA 38.586343
#> Bromeliaceae              NA        NA           NA         NA  9.750073
#> Brunelliaceae             NA        NA           NA         NA 16.166559
#> Bruniaceae                NA        NA           NA         NA  8.848485
#> Bryaceae                  NA        NA           NA         NA 16.050000
#> Burseraceae        14.980000        NA 0.0041680952         NA 18.934166
#> Butomaceae                NA        NA           NA         NA 42.600000
#> Buxaceae           13.741165        NA 0.0015158428         NA 19.119108
#> Cabombaceae               NA        NA           NA         NA 19.500000
#> Cactaceae           8.700000        NA 0.0026724592 0.14900000 15.314267
#> Calomniaceae              NA        NA           NA         NA  6.280000
#> Calophyllaceae            NA        NA           NA 0.33000000 13.429887
#>                    Nsapwood Nfineroot Kmax_stemxylem VCstem_P50      Al2As
#> Acanthaceae              NA  7.805074   7.539139e-01 -5.2100000  1620.9420
#> Achariaceae              NA        NA             NA -3.1500000  8300.0000
#> Achatocarpaceae          NA        NA             NA         NA         NA
#> Acoraceae                NA        NA             NA         NA         NA
#> Actinidiaceae            NA 24.905621   1.920380e+00 -1.1500000  7057.3160
#> Aextoxicaceae            NA        NA   3.800000e-01         NA 66666.6667
#> Aizoaceae                NA        NA             NA         NA         NA
#> Akaniaceae               NA        NA             NA         NA         NA
#> Alismataceae             NA 35.985750             NA         NA         NA
#> Alstroemeriaceae         NA        NA             NA         NA         NA
#> Altingiaceae             NA 10.723118   1.715415e+00 -3.0834551  6700.0869
#> Alzateaceae              NA        NA   3.216863e+00         NA         NA
#> Amaranthaceae            NA 10.971150   1.735153e-01 -2.5064151   973.0000
#> Amaryllidaceae           NA 11.188333             NA         NA         NA
#> Amborellaceae            NA        NA   4.370000e-01 -2.9985149  4260.0000
#> Anacardiaceae            NA 16.859303   4.257418e+00 -4.5493004 22515.7911
#> Anemiaceae               NA        NA             NA         NA         NA
#> Anisophylleaceae         NA        NA             NA         NA         NA
#> Annonaceae               NA 20.282433   5.142811e+00 -1.8571897  5488.1003
#> Aphloiaceae              NA        NA             NA         NA         NA
#> Apiaceae                 NA 10.599666   5.200000e-01 -5.7000000  2174.4000
#> Apocynaceae              NA 13.757041   1.047348e+01 -1.7282466  4389.3290
#> Aquifoliaceae            NA 10.155075   7.937288e-01 -3.6679851  9508.1433
#> Araceae                  NA        NA             NA         NA         NA
#> Araliaceae               NA 16.451735   1.363396e+00 -2.8447366  4828.2712
#> Araucariaceae            NA 13.000000   3.106368e-04 -2.6875469  3865.5204
#> Arecaceae                NA 12.267625   7.780000e+00 -1.8100000  5821.8050
#> Aristolochiaceae         NA        NA             NA         NA         NA
#> Asparagaceae             NA 12.888227   3.450000e-03 -2.1080000         NA
#> Asphodelaceae            NA        NA             NA         NA         NA
#> Aspleniaceae             NA 11.177640             NA         NA         NA
#> Asteliaceae              NA  8.305000             NA         NA         NA
#> Asteraceae               NA 12.814810   5.923996e-01 -3.0368525  2248.5141
#> Asteropeiaceae           NA        NA             NA         NA         NA
#> Atherospermataceae       NA 23.200000   5.508007e-01 -3.3100000  3320.1431
#> Athyriaceae              NA        NA             NA         NA         NA
#> Aulacomniaceae           NA        NA             NA         NA         NA
#> Austrobaileyaceae        NA        NA             NA -0.4990000         NA
#> Balanopaceae             NA        NA             NA         NA         NA
#> Balsaminaceae            NA 17.328412             NA         NA         NA
#> Basellaceae              NA        NA             NA         NA         NA
#> Begoniaceae              NA        NA             NA         NA         NA
#> Berberidaceae            NA 20.252051   1.400000e-01 -4.5130097  2200.5358
#> Berberidopsidaceae       NA        NA             NA         NA         NA
#> Betulaceae         12.53400 12.843458   2.256390e+00 -2.6037822  5727.7371
#> Bignoniaceae       12.96411 19.420246   2.249788e+00 -1.8361890  3932.1200
#> Bixaceae                 NA        NA   7.750000e+00 -1.4416163  1757.4692
#> Blechnaceae              NA  9.353304   1.215000e+01         NA         NA
#> Bonnetiaceae             NA        NA             NA         NA         NA
#> Boraginaceae             NA 17.421666   5.798530e+00 -3.2601286  9427.6448
#> Brassicaceae             NA 18.495733   1.617493e-02 -3.1543750         NA
#> Bromeliaceae             NA        NA             NA         NA         NA
#> Brunelliaceae            NA        NA   1.097763e+00         NA         NA
#> Bruniaceae               NA        NA   3.702603e-01 -3.0449976         NA
#> Bryaceae                 NA        NA             NA         NA         NA
#> Burseraceae              NA 12.769375   2.593880e+00 -2.3795462 12457.0704
#> Butomaceae               NA        NA             NA         NA         NA
#> Buxaceae                 NA        NA   1.393277e-01 -8.0000000         NA
#> Cabombaceae              NA        NA             NA         NA         NA
#> Cactaceae                NA        NA   9.700000e-01 -2.0000000  1479.2899
#> Calomniaceae             NA        NA             NA         NA         NA
#> Calophyllaceae           NA        NA   7.182893e+00 -0.9723964   186.1386
#>                    conduit2sapwood
#> Acanthaceae              0.6300000
#> Achariaceae                     NA
#> Achatocarpaceae                 NA
#> Acoraceae                       NA
#> Actinidiaceae                   NA
#> Aextoxicaceae                   NA
#> Aizoaceae                       NA
#> Akaniaceae                      NA
#> Alismataceae                    NA
#> Alstroemeriaceae                NA
#> Altingiaceae             0.8015000
#> Alzateaceae                     NA
#> Amaranthaceae                   NA
#> Amaryllidaceae                  NA
#> Amborellaceae                   NA
#> Anacardiaceae            0.7105300
#> Anemiaceae                      NA
#> Anisophylleaceae                NA
#> Annonaceae               0.5603750
#> Aphloiaceae                     NA
#> Apiaceae                        NA
#> Apocynaceae              0.7112500
#> Aquifoliaceae            0.6528500
#> Araceae                         NA
#> Araliaceae               0.7785000
#> Araucariaceae            0.9375000
#> Arecaceae                       NA
#> Aristolochiaceae                NA
#> Asparagaceae                    NA
#> Asphodelaceae                   NA
#> Aspleniaceae                    NA
#> Asteliaceae                     NA
#> Asteraceae               0.7098857
#> Asteropeiaceae                  NA
#> Atherospermataceae       0.7560000
#> Athyriaceae                     NA
#> Aulacomniaceae                  NA
#> Austrobaileyaceae               NA
#> Balanopaceae                    NA
#> Balsaminaceae                   NA
#> Basellaceae                     NA
#> Begoniaceae                     NA
#> Berberidaceae                   NA
#> Berberidopsidaceae              NA
#> Betulaceae               0.8305000
#> Bignoniaceae             0.6358032
#> Bixaceae                 0.7070000
#> Blechnaceae                     NA
#> Bonnetiaceae                    NA
#> Boraginaceae             0.6426000
#> Brassicaceae                    NA
#> Bromeliaceae                    NA
#> Brunelliaceae                   NA
#> Bruniaceae                      NA
#> Bryaceae                        NA
#> Burseraceae              0.8252848
#> Butomaceae                      NA
#> Buxaceae                 0.8330000
#> Cabombaceae                     NA
#> Cactaceae                0.3574375
#> Calomniaceae                    NA
#> Calophyllaceae           0.7335000
#>  [ reached 'max' / getOption("max.print") -- omitted 316 rows ]
```
