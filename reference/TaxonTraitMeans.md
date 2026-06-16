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
#>                    WoodDensity LeafDensity     WoodC    LeafPI0   LeafEPS
#> Acanthaceae         0.59819866  0.23109017 0.4108250 -3.3950000 23.230000
#> Achariaceae         0.59291465  0.36224725        NA         NA        NA
#> Achatocarpaceae     0.87000000          NA        NA         NA        NA
#> Acoraceae           0.21506871  0.10000000        NA         NA        NA
#> Actinidiaceae       0.36619865  0.38537880        NA         NA        NA
#> Aextoxicaceae       0.56654642          NA        NA         NA        NA
#> Aizoaceae           0.13097800  0.08241249        NA         NA        NA
#> Akaniaceae          0.56430643          NA        NA         NA        NA
#> Alismataceae        0.10037916  0.17703114        NA         NA        NA
#> Alstroemeriaceae            NA          NA        NA         NA        NA
#> Altingiaceae        0.57931303  0.68339710 0.4435000         NA        NA
#> Alzateaceae         0.66180271          NA        NA         NA        NA
#> Amaranthaceae       0.44919143  0.19925998 0.4199317 -1.4810943  5.641365
#> Amaryllidaceae      0.12438751  0.12749970 0.4508400         NA        NA
#> Amborellaceae               NA          NA        NA         NA        NA
#> Anacardiaceae       0.55251214  0.42726402 0.4580588 -2.0183556 12.760000
#> Anemiaceae                  NA  0.44850000        NA -1.3700000 20.040000
#> Anisophylleaceae    0.59253100          NA        NA         NA        NA
#> Annonaceae          0.53275262  0.36479085 0.4683333 -2.5250000 34.215000
#> Aphloiaceae         0.61962314  0.46270602        NA         NA        NA
#> Apiaceae            0.26837295  0.32662503 0.4382558 -1.1040000 23.410000
#> Apocynaceae         0.54132595  0.29921186 0.4928286 -1.4950000 17.955000
#> Aquifoliaceae       0.57615551  0.36690261 0.4400000 -2.3000000 20.700000
#> Araceae             0.08849923  0.17849835        NA         NA        NA
#> Araliaceae          0.41094541  0.28151142 0.4517143 -1.4535028 11.387712
#> Araucariaceae       0.46044971  0.36739509        NA         NA        NA
#> Arecaceae           0.57661947  0.42339049 0.4556000 -3.4000000 73.400000
#> Aristolochiaceae    0.24502011  0.26814343        NA         NA        NA
#> Asparagaceae        0.28530769  0.18761443 0.4631000 -1.7922006 22.810356
#> Asphodelaceae       0.41615533  0.55721844        NA -1.3794753  7.425622
#> Aspleniaceae        0.30000000  0.34159513        NA -1.2400000 35.300000
#> Asteliaceae                 NA          NA        NA         NA        NA
#> Asteraceae          0.40106594  0.23275728 0.4291509 -1.3018058  6.022339
#> Asteropeiaceae      0.81479790          NA        NA         NA        NA
#> Atherospermataceae  0.49540018  0.24008723        NA -1.3400000  8.380000
#> Athyriaceae                 NA  0.13005000        NA         NA        NA
#> Aulacomniaceae              NA          NA        NA         NA        NA
#> Austrobaileyaceae   0.36334000  0.26209917        NA         NA        NA
#> Balanopaceae        0.72927837          NA        NA         NA        NA
#> Balsaminaceae       0.16107319  0.38966452        NA         NA        NA
#> Basellaceae         0.20000000          NA        NA         NA        NA
#> Begoniaceae         0.05056471  0.16267894        NA         NA        NA
#> Berberidaceae       0.66892669  0.39692512        NA -1.7100000        NA
#> Berberidopsidaceae  0.33000000          NA        NA         NA        NA
#> Betulaceae          0.57153544  0.49433425 0.4750059 -1.3583989  4.458000
#> Bignoniaceae        0.56821892  0.35517753 0.4692962 -1.9900000 17.610000
#> Bixaceae            0.32813746  0.36302045        NA         NA        NA
#> Blechnaceae         0.37060000  0.31479929        NA         NA        NA
#> Bonnetiaceae        0.79986944  0.27978289        NA         NA        NA
#> Boraginaceae        0.46228583  0.30176908 0.4536333 -1.5210000 11.797333
#> Brassicaceae        0.20428225  0.24471597 0.4244838 -1.4200000        NA
#> Bromeliaceae        0.22321054  0.17755239        NA         NA        NA
#> Brunelliaceae       0.37862442          NA        NA         NA        NA
#> Bruniaceae          0.56395714          NA        NA         NA        NA
#> Bryaceae                    NA          NA        NA         NA        NA
#> Burseraceae         0.51952671  0.42079782 0.4703312 -1.4350000 14.980000
#> Butomaceae          0.10083260          NA        NA         NA        NA
#> Buxaceae            0.60345652  0.27062945        NA         NA        NA
#> Cabombaceae                 NA          NA        NA         NA        NA
#> Cactaceae           0.57064278  0.15317264        NA         NA  8.700000
#> Calomniaceae                NA          NA        NA         NA        NA
#> Calophyllaceae      0.61022831  0.39397650 0.4729500         NA        NA
#> Calycanthaceae      0.62946278  0.36542389        NA         NA        NA
#> Calyceraceae        0.16399000  0.16893665        NA         NA        NA
#> Campanulaceae       0.22357033  0.22938075 0.4489120 -0.8345534  3.195661
#> Canellaceae         0.66653305          NA        NA         NA        NA
#>                    LeafAF       Gswmin    Gswmax     Nleaf Nsapwood Nfineroot
#> Acanthaceae        0.1270           NA 0.2500000 28.166211       NA  7.805074
#> Achariaceae            NA           NA        NA 22.638234       NA        NA
#> Achatocarpaceae        NA           NA        NA        NA       NA        NA
#> Acoraceae              NA           NA        NA 18.000000       NA        NA
#> Actinidiaceae          NA           NA        NA 22.057101       NA 24.905621
#> Aextoxicaceae          NA           NA        NA  9.652280       NA        NA
#> Aizoaceae              NA           NA        NA 14.800000       NA        NA
#> Akaniaceae             NA           NA        NA        NA       NA        NA
#> Alismataceae           NA           NA        NA 24.279413       NA 35.985750
#> Alstroemeriaceae       NA           NA        NA 18.484000       NA        NA
#> Altingiaceae           NA 0.0031009001 0.5000000 15.704509       NA 10.723118
#> Alzateaceae            NA           NA        NA  9.802997       NA        NA
#> Amaranthaceae      0.6930 0.0113780532 0.0750000 25.143267       NA 10.971150
#> Amaryllidaceae         NA 0.0093742733        NA 27.807973       NA 11.188333
#> Amborellaceae          NA           NA        NA        NA       NA        NA
#> Anacardiaceae          NA 0.0043356163 0.2017136 17.851435       NA 16.859303
#> Anemiaceae         0.2550           NA        NA 21.440000       NA        NA
#> Anisophylleaceae       NA           NA        NA        NA       NA        NA
#> Annonaceae             NA 0.0045999511        NA 23.182512       NA 20.282433
#> Aphloiaceae            NA           NA        NA        NA       NA        NA
#> Apiaceae           0.8380 0.0051617034        NA 25.885090       NA 10.599666
#> Apocynaceae        0.7930 0.0030300000        NA 21.091263       NA 13.757041
#> Aquifoliaceae      0.4000 0.0017826592 0.2600000 15.239217       NA 10.155075
#> Araceae                NA 0.0042445099        NA 19.587660       NA        NA
#> Araliaceae         0.5620 0.0005668349        NA 17.527488       NA 16.451735
#> Araucariaceae          NA 0.0018785800        NA 18.412875       NA 13.000000
#> Arecaceae          0.2000 0.0009988400        NA 17.618260       NA 12.267625
#> Aristolochiaceae       NA 0.0035895556        NA 31.050092       NA        NA
#> Asparagaceae           NA 0.0010977161        NA 24.878019       NA 12.888227
#> Asphodelaceae          NA 0.0055481213        NA 10.880839       NA        NA
#> Aspleniaceae           NA 0.0045607354        NA 25.425769       NA 11.177640
#> Asteliaceae            NA           NA        NA  7.205741       NA  8.305000
#> Asteraceae         0.4066 0.0107116006 0.1976901 21.823946       NA 12.814810
#> Asteropeiaceae         NA           NA        NA        NA       NA        NA
#> Atherospermataceae     NA           NA        NA 14.776944       NA 23.200000
#> Athyriaceae            NA           NA        NA 27.001204       NA        NA
#> Aulacomniaceae         NA           NA        NA  8.000000       NA        NA
#> Austrobaileyaceae      NA           NA        NA        NA       NA        NA
#> Balanopaceae           NA           NA        NA        NA       NA        NA
#> Balsaminaceae          NA 0.0126573016        NA 35.893386       NA 17.328412
#> Basellaceae            NA           NA        NA        NA       NA        NA
#> Begoniaceae            NA           NA        NA 27.734506       NA        NA
#> Berberidaceae          NA 0.0016500000        NA 20.459189       NA 20.252051
#> Berberidopsidaceae     NA           NA        NA        NA       NA        NA
#> Betulaceae             NA 0.0031703165 0.3536968 24.816045 12.53400 12.843458
#> Bignoniaceae       0.1770 0.0026880454        NA 25.194014 12.96411 19.420246
#> Bixaceae               NA 0.0027600000        NA 20.220367       NA        NA
#> Blechnaceae            NA           NA        NA 11.907284       NA  9.353304
#> Bonnetiaceae           NA           NA        NA 10.472527       NA        NA
#> Boraginaceae           NA 0.0060134744 0.2900000 25.911014       NA 17.421666
#> Brassicaceae           NA 0.0113630344        NA 39.093972       NA 18.495733
#> Bromeliaceae           NA           NA        NA  8.659122       NA        NA
#> Brunelliaceae          NA           NA        NA 16.166559       NA        NA
#> Bruniaceae             NA           NA        NA  7.785422       NA        NA
#> Bryaceae               NA           NA        NA 16.050000       NA        NA
#> Burseraceae            NA 0.0027100000        NA 18.524136       NA 12.769375
#> Butomaceae             NA           NA        NA 42.600000       NA        NA
#> Buxaceae               NA           NA        NA 19.346119       NA        NA
#> Cabombaceae            NA           NA        NA 19.500000       NA        NA
#> Cactaceae              NA 0.0014808049        NA 16.597289       NA        NA
#> Calomniaceae           NA           NA        NA  6.280000       NA        NA
#> Calophyllaceae         NA           NA 0.1350000 13.474516       NA        NA
#> Calycanthaceae         NA           NA        NA 21.297115       NA 19.029600
#> Calyceraceae           NA           NA        NA 38.773333       NA        NA
#> Campanulaceae          NA 0.0059264620        NA 27.418234       NA 11.512892
#> Canellaceae            NA           NA        NA        NA       NA        NA
#>                    Kmax_stemxylem VCstem_P50      Al2As conduit2sapwood
#> Acanthaceae             3.6524786 -5.7450000  1803.6645       0.6300000
#> Achariaceae                    NA -3.1500000  8300.0000              NA
#> Achatocarpaceae                NA         NA         NA              NA
#> Acoraceae                      NA         NA         NA              NA
#> Actinidiaceae           1.4453860 -0.9833333 11543.4188              NA
#> Aextoxicaceae           0.3800000         NA 66666.6667              NA
#> Aizoaceae                      NA         NA         NA              NA
#> Akaniaceae                     NA         NA         NA              NA
#> Alismataceae                   NA         NA         NA              NA
#> Alstroemeriaceae               NA         NA         NA              NA
#> Altingiaceae            1.3055515 -2.7349020  6472.3731       0.8015000
#> Alzateaceae             3.2168629         NA         NA              NA
#> Amaranthaceae           0.1564636 -3.0284442   973.0677              NA
#> Amaryllidaceae                 NA         NA         NA              NA
#> Amborellaceae           0.5400000 -3.0000000  4257.6596              NA
#> Anacardiaceae           3.5289937 -2.6108245 21764.4502       0.7105300
#> Anemiaceae                     NA         NA         NA              NA
#> Anisophylleaceae               NA         NA         NA              NA
#> Annonaceae              4.3383977 -3.1516340  9203.3824       0.5603750
#> Aphloiaceae                    NA         NA         NA              NA
#> Apiaceae                0.5150000 -5.7000000  1389.4314              NA
#> Apocynaceae             1.9545389 -2.2883044 20941.4964       0.7112500
#> Aquifoliaceae           0.7222116 -4.2944620  7414.6991       0.6528500
#> Araceae                        NA         NA         NA              NA
#> Araliaceae              1.6541078 -2.0065789  5044.9277       0.7785000
#> Araucariaceae           0.7378000 -2.7808696  3860.5770       0.9375000
#> Arecaceae               4.0884057 -1.8100000  9560.5163              NA
#> Aristolochiaceae       45.0000000         NA         NA              NA
#> Asparagaceae            0.0034500 -2.1080000         NA              NA
#> Asphodelaceae                  NA         NA         NA              NA
#> Aspleniaceae                   NA         NA         NA              NA
#> Asteliaceae                    NA         NA         NA              NA
#> Asteraceae              2.4028983 -3.3353531  2465.6724       0.7098857
#> Asteropeiaceae                 NA         NA         NA              NA
#> Atherospermataceae      0.5500000 -3.3118889  5080.5040       0.7560000
#> Athyriaceae                    NA         NA         NA              NA
#> Aulacomniaceae                 NA         NA         NA              NA
#> Austrobaileyaceae       2.3000000 -0.4990000 15384.6154              NA
#> Balanopaceae                   NA         NA         NA              NA
#> Balsaminaceae                  NA         NA         NA              NA
#> Basellaceae                    NA         NA         NA              NA
#> Begoniaceae                    NA         NA         NA              NA
#> Berberidaceae           0.2374195 -4.5000000  3461.9656              NA
#> Berberidopsidaceae             NA         NA         NA              NA
#> Betulaceae              2.6700672 -2.0812191  6262.0017       0.8305000
#> Bignoniaceae            1.9256324 -1.5351814 15502.7818       0.6325556
#> Bixaceae                7.7500000 -1.4403556  1758.3128       0.7070000
#> Blechnaceae                    NA         NA         NA              NA
#> Bonnetiaceae                   NA         NA         NA              NA
#> Boraginaceae            4.9482106 -2.3958000  9807.0622       0.6300000
#> Brassicaceae                   NA         NA         NA              NA
#> Bromeliaceae                   NA         NA         NA              NA
#> Brunelliaceae           1.0977628         NA         NA              NA
#> Bruniaceae              0.3807548 -3.0448784         NA              NA
#> Bryaceae                       NA         NA         NA              NA
#> Burseraceae             4.1829227 -1.5985357 16869.9644       0.8204286
#> Butomaceae                     NA         NA         NA              NA
#> Buxaceae                       NA -8.0000000         NA       0.8330000
#> Cabombaceae                    NA         NA         NA              NA
#> Cactaceae               1.6672222 -2.0000000  2339.3221       0.3574375
#> Calomniaceae                   NA         NA         NA              NA
#> Calophyllaceae          2.3932233 -1.7505034   186.1379       0.7335000
#> Calycanthaceae          0.4100000 -1.9780000  3400.6803       0.6535500
#> Calyceraceae                   NA         NA         NA              NA
#> Campanulaceae                  NA         NA         NA              NA
#> Canellaceae                    NA -1.7880000         NA              NA
#>  [ reached 'max' / getOption("max.print") -- omitted 312 rows ]
```
