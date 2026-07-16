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
#>                    WoodDensity LeafDensity FineRootDensity     WoodC       Ptlp
#> Acanthaceae         0.62347293  0.29727615       0.3884868 0.4108250 -2.7400000
#> Achariaceae         0.61001599  0.36224725              NA        NA -2.4200000
#> Achatocarpaceae     0.87000000          NA              NA        NA         NA
#> Acoraceae           0.21506871  0.10000000       0.0800000        NA         NA
#> Actinidiaceae       0.37145115  0.45142686       0.1078125        NA -0.7897536
#> Aextoxicaceae       0.56666667          NA              NA        NA         NA
#> Aizoaceae           0.13097800  0.08241249              NA        NA         NA
#> Akaniaceae          0.57133685          NA              NA        NA         NA
#> Alismataceae        0.10037916  0.17703114       0.1921740        NA         NA
#> Alstroemeriaceae            NA          NA              NA        NA         NA
#> Altingiaceae        0.53151501  0.68339710       0.3722962 0.4435000 -1.2487962
#> Alzateaceae         0.66180271          NA              NA        NA         NA
#> Amaranthaceae       0.68154337  0.17128218       0.4149995 0.4199317 -1.5450000
#> Amaryllidaceae      0.12438751  0.12749970       0.1049986 0.4508400         NA
#> Amborellaceae               NA          NA              NA        NA         NA
#> Anacardiaceae       0.53134629  0.51158932       0.3659151 0.4580588 -2.8013368
#> Anemiaceae                  NA  0.44850000              NA        NA -1.6450000
#> Anisophylleaceae    0.57468312          NA              NA        NA         NA
#> Annonaceae          0.55383625  0.41099949       0.2443377 0.4683333 -1.9433333
#> Aphloiaceae         0.62052000  0.46270602              NA        NA         NA
#> Apiaceae            0.56064111  0.23171368       0.1438487 0.4382558 -1.5040000
#> Apocynaceae         0.54159387  0.32835059       0.2621171 0.4928286 -0.7800000
#> Aquifoliaceae       0.58119398  0.43610051       0.2106023 0.4400000 -1.9144202
#> Araceae             0.06372500  0.17849835       0.1974646        NA         NA
#> Araliaceae          0.41671464  0.31681165       0.1395793 0.4517143 -1.5407671
#> Araucariaceae       0.42218480  0.38589002       0.2321983        NA         NA
#> Arecaceae           0.55204229  0.42568643       0.2593264 0.4556000 -3.8100000
#> Aristolochiaceae    0.28046343  0.26814343              NA        NA         NA
#> Asparagaceae        0.65496705  0.18761443       0.1345246 0.4631000         NA
#> Asphodelaceae       0.41942475  0.55721844              NA        NA         NA
#> Aspleniaceae        0.30000000  0.34159513       0.2339205        NA -1.3700000
#> Asteliaceae                 NA          NA       0.0900000        NA         NA
#> Asteraceae          0.52236825  0.22487138       0.1770911 0.4291509 -1.5967910
#> Asteropeiaceae      0.81717121          NA              NA        NA         NA
#> Atherospermataceae  0.53917542  0.24008723       0.1125824        NA -1.7600000
#> Athyriaceae                 NA  0.13005000              NA        NA         NA
#> Aulacomniaceae              NA          NA              NA        NA         NA
#> Austrobaileyaceae   0.36334000  0.26209917              NA        NA         NA
#> Balanopaceae        0.68270079          NA              NA        NA         NA
#> Balsaminaceae       0.16107319  0.38966452              NA        NA         NA
#> Basellaceae         0.20000000          NA              NA        NA         NA
#> Begoniaceae         0.05056471  0.16267894              NA        NA         NA
#> Berberidaceae       0.67673862  0.39692512       0.2867463        NA -1.8233333
#> Berberidopsidaceae  0.33000000          NA              NA        NA         NA
#> Betulaceae          0.55730771  0.49528026       0.3900151 0.4750059 -1.9169583
#> Bignoniaceae        0.52135729  0.35958354       0.2034424 0.4692962 -1.8134743
#> Bixaceae            0.40709211  0.48683360              NA        NA -1.4500000
#> Blechnaceae         0.37060000  0.31479929       0.3100000        NA         NA
#> Bonnetiaceae        0.78381722  0.27978289              NA        NA         NA
#> Boraginaceae        0.47932589  0.37598065       0.2088358 0.4536333 -1.6333333
#> Brassicaceae        0.46620780  0.24471597       0.3473000 0.4244838 -1.8700000
#> Bromeliaceae        0.22321054  0.17755239              NA        NA         NA
#> Brunelliaceae       0.43945805          NA              NA        NA         NA
#> Bruniaceae          0.56436667          NA              NA        NA         NA
#> Bryaceae                    NA          NA              NA        NA         NA
#>                       LeafPI0   LeafEPS    LeafAF       Gswmin     Gswmax
#> Acanthaceae        -3.3950000 23.230000 0.1270000           NA 0.25000000
#> Achariaceae                NA        NA        NA           NA         NA
#> Achatocarpaceae            NA        NA        NA           NA         NA
#> Acoraceae                  NA        NA        NA           NA         NA
#> Actinidiaceae              NA        NA        NA           NA         NA
#> Aextoxicaceae              NA        NA        NA           NA         NA
#> Aizoaceae                  NA        NA        NA           NA         NA
#> Akaniaceae                 NA        NA        NA           NA         NA
#> Alismataceae               NA        NA        NA           NA         NA
#> Alstroemeriaceae           NA        NA        NA           NA         NA
#> Altingiaceae               NA        NA        NA 0.0031009001 0.50000000
#> Alzateaceae                NA        NA        NA           NA         NA
#> Amaranthaceae      -1.2679576 15.476547 0.6930000 0.0113780532 0.92710526
#> Amaryllidaceae             NA        NA        NA 0.0022841333         NA
#> Amborellaceae              NA        NA        NA           NA         NA
#> Anacardiaceae      -2.4428297 12.760000        NA 0.0036981857 0.06956795
#> Anemiaceae         -1.3700000 20.040000 0.2550000           NA         NA
#> Anisophylleaceae           NA        NA        NA           NA         NA
#> Annonaceae         -2.5250000 34.215000        NA 0.0053616476         NA
#> Aphloiaceae                NA        NA        NA           NA         NA
#> Apiaceae           -1.1040000 23.410000 0.8380000 0.0080862140         NA
#> Apocynaceae        -0.6000000 14.970000 0.7930000 0.0030300000 1.49093750
#> Aquifoliaceae      -2.3000000 20.700000 0.4000000 0.0027715622 0.26000000
#> Araceae                    NA        NA        NA 0.0042445099         NA
#> Araliaceae         -0.7673038 12.394746 0.8730000 0.0023960466 1.33950000
#> Araucariaceae              NA        NA        NA 0.0018785800         NA
#> Arecaceae          -3.4000000 73.400000 0.2000000 0.0008583501         NA
#> Aristolochiaceae           NA        NA        NA 0.0035895556         NA
#> Asparagaceae       -1.7922006 22.810356        NA 0.0010977161         NA
#> Asphodelaceae      -1.3794753  7.425622        NA 0.0055481213         NA
#> Aspleniaceae       -1.2400000 35.300000        NA 0.0015831553         NA
#> Asteliaceae                NA        NA        NA           NA         NA
#> Asteraceae         -1.4119700  9.586856 0.5153333 0.0113591271 0.32161918
#> Asteropeiaceae             NA        NA        NA           NA         NA
#> Atherospermataceae -1.3400000  8.380000        NA           NA 0.55458333
#> Athyriaceae                NA        NA        NA           NA         NA
#> Aulacomniaceae             NA        NA        NA           NA         NA
#> Austrobaileyaceae          NA        NA        NA           NA         NA
#> Balanopaceae               NA        NA        NA           NA         NA
#> Balsaminaceae              NA        NA        NA 0.0126573016         NA
#> Basellaceae                NA        NA        NA           NA         NA
#> Begoniaceae                NA        NA        NA           NA         NA
#> Berberidaceae      -1.7100000        NA        NA 0.0016500000         NA
#> Berberidopsidaceae         NA        NA        NA           NA         NA
#> Betulaceae         -1.5015389  3.306000        NA 0.0027260202 0.13051410
#> Bignoniaceae       -1.9900000 17.610000 0.1770000 0.0036919354         NA
#> Bixaceae                   NA        NA        NA 0.0027600000         NA
#> Blechnaceae                NA        NA        NA           NA         NA
#> Bonnetiaceae               NA        NA        NA           NA         NA
#> Boraginaceae       -1.2177273 12.128750        NA 0.0060134744 0.29000000
#> Brassicaceae       -1.4800000        NA        NA 0.0119382443         NA
#> Bromeliaceae               NA        NA        NA           NA         NA
#> Brunelliaceae              NA        NA        NA           NA         NA
#> Bruniaceae                 NA        NA        NA           NA         NA
#> Bryaceae                   NA        NA        NA           NA         NA
#>                        Nleaf Nsapwood Nfineroot Kmax_stemxylem VCstem_P50
#> Acanthaceae        25.836769       NA  7.805074   7.539139e-01  -5.210000
#> Achariaceae        21.945455       NA        NA             NA  -3.150000
#> Achatocarpaceae           NA       NA        NA             NA         NA
#> Acoraceae          18.000000       NA        NA             NA         NA
#> Actinidiaceae      21.635986       NA 24.905621   1.920380e+00  -1.150000
#> Aextoxicaceae      10.121739       NA        NA   3.800000e-01         NA
#> Aizoaceae          14.800000       NA        NA             NA         NA
#> Akaniaceae                NA       NA        NA             NA         NA
#> Alismataceae       31.388709       NA 35.985750             NA         NA
#> Alstroemeriaceae   18.126667       NA        NA             NA         NA
#> Altingiaceae       16.490558       NA 10.723118   1.715415e+00  -3.083455
#> Alzateaceae         9.802997       NA        NA   3.216863e+00         NA
#> Amaranthaceae      25.940024       NA 10.971150   1.734250e-01  -2.506415
#> Amaryllidaceae     27.680470       NA 11.188333             NA         NA
#> Amborellaceae             NA       NA        NA   4.370000e-01  -3.000000
#> Anacardiaceae      18.853621       NA 16.859303   4.257418e+00  -4.549300
#> Anemiaceae         21.440000       NA        NA             NA         NA
#> Anisophylleaceae          NA       NA        NA             NA         NA
#> Annonaceae         23.501445       NA 20.282433   5.142811e+00  -1.857190
#> Aphloiaceae               NA       NA        NA             NA         NA
#> Apiaceae           26.373905       NA 10.599666   5.200000e-01  -5.700000
#> Apocynaceae        21.079302       NA 13.757041   1.047348e+01  -1.728247
#> Aquifoliaceae      14.636919       NA 10.155075   7.937288e-01  -3.667985
#> Araceae            19.245495       NA        NA             NA         NA
#> Araliaceae         18.090127       NA 16.451735   1.363396e+00  -2.844737
#> Araucariaceae      17.428378       NA 13.000000   3.106368e-04  -2.687547
#> Arecaceae          19.231886       NA 12.267625   7.780000e+00  -1.810000
#> Aristolochiaceae   31.770806       NA        NA             NA         NA
#> Asparagaceae       24.786812       NA 12.888227             NA         NA
#> Asphodelaceae      13.623111       NA        NA             NA         NA
#> Aspleniaceae       25.551582       NA 11.177640             NA         NA
#> Asteliaceae         7.205741       NA  8.305000             NA         NA
#> Asteraceae         22.667037       NA 12.814810   5.923996e-01  -3.036852
#> Asteropeiaceae            NA       NA        NA             NA         NA
#> Atherospermataceae 15.392333       NA 23.200000   5.500000e-01  -3.310000
#> Athyriaceae        28.101444       NA        NA             NA         NA
#> Aulacomniaceae      8.000000       NA        NA             NA         NA
#> Austrobaileyaceae         NA       NA        NA             NA  -0.499000
#> Balanopaceae              NA       NA        NA             NA         NA
#> Balsaminaceae      38.195871       NA 17.328412             NA         NA
#> Basellaceae               NA       NA        NA             NA         NA
#> Begoniaceae        28.277827       NA        NA             NA         NA
#> Berberidaceae      21.510712       NA 20.252051   1.400000e-01  -4.500000
#> Berberidopsidaceae        NA       NA        NA             NA         NA
#> Betulaceae         24.889296 12.53400 12.843458   2.256390e+00  -2.603782
#> Bignoniaceae       25.851301 12.96411 19.420246   2.249788e+00  -1.836189
#> Bixaceae           19.032367       NA        NA   7.750000e+00  -1.441616
#> Blechnaceae        11.745388       NA  9.353304             NA         NA
#> Bonnetiaceae       10.808790       NA        NA             NA         NA
#> Boraginaceae       25.567842       NA 17.421666   5.798530e+00  -3.260129
#> Brassicaceae       38.586343       NA 18.495733             NA         NA
#> Bromeliaceae        9.750073       NA        NA             NA         NA
#> Brunelliaceae      16.166559       NA        NA   1.097763e+00         NA
#> Bruniaceae          8.848485       NA        NA   3.687956e-01  -3.045000
#> Bryaceae           16.050000       NA        NA             NA         NA
#>                        Al2As conduit2sapwood       SRL
#> Acanthaceae         1620.942       0.6300000  2627.503
#> Achariaceae         8300.000              NA        NA
#> Achatocarpaceae           NA              NA        NA
#> Acoraceae                 NA              NA        NA
#> Actinidiaceae       7057.316              NA 22452.377
#> Aextoxicaceae      66666.667              NA        NA
#> Aizoaceae                 NA              NA        NA
#> Akaniaceae                NA              NA        NA
#> Alismataceae              NA              NA  2999.158
#> Alstroemeriaceae          NA              NA        NA
#> Altingiaceae        6696.120       0.8015000  3083.189
#> Alzateaceae               NA              NA        NA
#> Amaranthaceae        973.000              NA  8633.782
#> Amaryllidaceae            NA              NA  3499.583
#> Amborellaceae       4260.000              NA        NA
#> Anacardiaceae      22671.398       0.7105300  4570.389
#> Anemiaceae                NA              NA        NA
#> Anisophylleaceae          NA              NA        NA
#> Annonaceae          5488.100       0.5603750  1299.778
#> Aphloiaceae               NA              NA        NA
#> Apiaceae            2174.400              NA 12803.648
#> Apocynaceae         4389.329       0.7112500  3436.547
#> Aquifoliaceae       9495.144       0.6528500  4087.608
#> Araceae                   NA              NA  3630.397
#> Araliaceae          4828.271       0.7785000  5159.123
#> Araucariaceae       3865.385       0.9375000  1342.881
#> Arecaceae           5821.805              NA  2819.556
#> Aristolochiaceae          NA              NA        NA
#> Asparagaceae              NA              NA  4119.521
#> Asphodelaceae             NA              NA  2535.250
#> Aspleniaceae              NA              NA  7377.539
#> Asteliaceae               NA              NA  1370.000
#> Asteraceae          2248.514       0.7098857 10301.279
#> Asteropeiaceae            NA              NA        NA
#> Atherospermataceae  3320.143       0.7560000  1050.719
#> Athyriaceae               NA              NA        NA
#> Aulacomniaceae            NA              NA        NA
#> Austrobaileyaceae         NA              NA        NA
#> Balanopaceae              NA              NA        NA
#> Balsaminaceae             NA              NA  3984.138
#> Basellaceae               NA              NA        NA
#> Begoniaceae               NA              NA        NA
#> Berberidaceae       2208.689              NA  2475.989
#> Berberidopsidaceae        NA              NA        NA
#> Betulaceae          5677.893       0.8305000  6367.772
#> Bignoniaceae        3932.120       0.6358032  1091.807
#> Bixaceae            1757.469       0.7070000        NA
#> Blechnaceae               NA              NA  5767.422
#> Bonnetiaceae              NA              NA        NA
#> Boraginaceae        9420.097       0.6426000 17122.563
#> Brassicaceae              NA              NA 35563.773
#> Bromeliaceae              NA              NA        NA
#> Brunelliaceae             NA              NA        NA
#> Bruniaceae                NA              NA        NA
#> Bryaceae                  NA              NA        NA
#>  [ reached 'max' / getOption("max.print") -- omitted 324 rows ]
```
