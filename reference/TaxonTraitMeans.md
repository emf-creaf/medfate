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
#>                    LeafDensity WoodDensity   LeafPI0   LeafEPS LeafAF    Gswmax
#> Acanthaceae         0.29590033   0.5684693 -3.395000 23.230000 0.1270 0.2500000
#> Achariaceae         0.25524948   0.6052036        NA        NA     NA        NA
#> Acoraceae           0.10000000          NA        NA        NA     NA        NA
#> Actinidiaceae       0.44397760   0.4092320        NA        NA     NA        NA
#> Adoxaceae           0.39161942   0.5157416 -1.560000 12.790000     NA        NA
#> Aextoxicaceae               NA   0.5666667        NA        NA     NA        NA
#> Aizoaceae           0.08241249          NA        NA        NA     NA        NA
#> Akaniaceae                  NA   0.5547825        NA        NA     NA        NA
#> Alismataceae        0.19142078          NA        NA        NA     NA        NA
#> Alstroemeriaceae            NA          NA        NA        NA     NA        NA
#> Altingiaceae        0.68339710   0.6010948        NA        NA     NA 0.5000000
#> Amaranthaceae       0.20458309   0.6315739 -2.250000        NA     NA 0.0750000
#> Amaryllidaceae      0.12207038          NA        NA        NA     NA        NA
#> Amborellaceae               NA          NA        NA        NA     NA        NA
#> Amphorogynaceae             NA   0.6097400        NA        NA     NA        NA
#> Anacardiaceae       0.45680882   0.5685583 -1.700000 12.760000     NA 0.3148333
#> Anisophylleaceae            NA   0.6734780        NA        NA     NA        NA
#> Annonaceae          0.37578113   0.5642062 -2.160000 23.710000     NA        NA
#> Aphloiaceae         0.46270602   0.6205200        NA        NA     NA        NA
#> Apiaceae            0.28747559   0.2561785        NA        NA     NA        NA
#> Apocynaceae         0.29811485   0.5683635 -2.390000 20.940000     NA        NA
#> Aptandraceae        0.33184918   0.7076756        NA        NA     NA        NA
#> Aquifoliaceae       0.48126262   0.5579305 -2.300000 20.700000 0.4000 0.2600000
#> Araceae             0.16597459          NA        NA        NA     NA        NA
#> Araliaceae          0.31188214   0.4142687 -1.528503 11.231462 0.4065        NA
#> Araucariaceae       0.34479018   0.4641456        NA        NA     NA        NA
#> Arecaceae           0.43935890   0.5913967 -3.400000 73.400000 0.2000        NA
#> Aristolochiaceae    0.26814343   0.2900000        NA        NA     NA        NA
#> Asparagaceae        0.14250995   0.4254258        NA        NA     NA        NA
#> Asphodelaceae       0.55124717   0.3951990        NA        NA     NA        NA
#> Aspleniaceae        0.26445566          NA -1.240000 35.300000     NA        NA
#> Asteraceae          0.25110926   0.4822337 -1.471389 14.491429 0.2435 0.1275000
#> Asteropeiaceae              NA   0.7554862        NA        NA     NA        NA
#> Atherospermataceae  0.24513566   0.4767621 -1.340000  8.380000     NA        NA
#> Athyriaceae         0.22706154          NA        NA        NA     NA        NA
#> Aulacomniaceae              NA          NA        NA        NA     NA        NA
#> Austrobaileyaceae   0.26209917          NA        NA        NA     NA        NA
#> Balanopaceae                NA   0.7348976        NA        NA     NA        NA
#> Balsaminaceae       0.39480943          NA        NA        NA     NA        NA
#> Begoniaceae         0.16267894          NA        NA        NA     NA        NA
#> Berberidaceae       0.35994527   0.7028850        NA        NA     NA        NA
#> Betulaceae          0.44465276   0.5381493 -1.246984  5.498667     NA 0.3965000
#> Bignoniaceae        0.35725900   0.6256030 -1.990000 17.610000 0.1770        NA
#> Bixaceae            0.44000000   0.3546357        NA        NA     NA        NA
#> Blechnaceae         0.31479929          NA        NA        NA     NA        NA
#> Bonnetiaceae        0.27978289   0.8400000        NA        NA     NA        NA
#> Boraginaceae        0.28504624   0.4987559 -1.140000        NA     NA        NA
#> Brassicaceae        0.22676049   0.4516377 -1.485000  7.710000 0.2320        NA
#> Bromeliaceae        0.16771395          NA        NA        NA     NA        NA
#> Brunelliaceae               NA   0.3112500        NA        NA     NA        NA
#> Bruniaceae                  NA   0.5636500        NA        NA     NA        NA
#> Burseraceae         0.42267872   0.5205008 -1.435000 14.980000     NA        NA
#> Butomaceae                  NA          NA        NA        NA     NA        NA
#> Buxaceae            0.27062945   0.7314511        NA        NA     NA        NA
#> Cabombaceae                 NA          NA        NA        NA     NA        NA
#> Cactaceae           0.18193720   0.6187500        NA  8.700000     NA        NA
#> Calophyllaceae      0.40650924   0.6067855        NA        NA     NA 0.1350000
#> Calycanthaceae      0.32040734   0.6500550        NA        NA     NA        NA
#> Calyceraceae        0.16893665          NA        NA        NA     NA        NA
#> Campanulaceae       0.23963177          NA        NA        NA     NA        NA
#> Canellaceae                 NA   0.6808255        NA        NA     NA        NA
#> Cannabaceae         0.35297198   0.5502170 -1.580000  5.180000     NA 0.3300000
#>                          Gswmin     Nleaf  Nsapwood Nfineroot     WoodC
#> Acanthaceae                  NA 27.122403        NA  5.580000 0.4242167
#> Achariaceae                  NA 22.162815        NA        NA        NA
#> Acoraceae                    NA 18.000000        NA        NA        NA
#> Actinidiaceae                NA 20.183321        NA        NA        NA
#> Adoxaceae                    NA 20.566364        NA 11.013750        NA
#> Aextoxicaceae                NA  9.621429        NA        NA        NA
#> Aizoaceae                    NA 14.800000        NA        NA        NA
#> Akaniaceae                   NA        NA        NA        NA        NA
#> Alismataceae                 NA 27.287834        NA        NA        NA
#> Alstroemeriaceae             NA 19.402857        NA        NA        NA
#> Altingiaceae                 NA 15.462547        NA  7.550000 0.4435000
#> Amaranthaceae                NA 23.740694        NA 12.411441 0.4197662
#> Amaryllidaceae     1.804000e-02 28.734076        NA 11.440000 0.4508400
#> Amborellaceae                NA        NA        NA        NA        NA
#> Amphorogynaceae              NA 25.236766        NA        NA        NA
#> Anacardiaceae      1.222994e-02 17.914439        NA 10.539736 0.4566867
#> Anisophylleaceae             NA        NA        NA        NA        NA
#> Annonaceae                   NA 23.593695        NA 23.830391 0.4726000
#> Aphloiaceae                  NA        NA        NA        NA        NA
#> Apiaceae           1.587302e-03 24.746891        NA 10.787434 0.4413076
#> Apocynaceae                  NA 21.754119        NA 16.845574 0.4924750
#> Aptandraceae                 NA 28.435524        NA        NA        NA
#> Aquifoliaceae      5.740000e-04 14.645818        NA 14.324450 0.4400000
#> Araceae                      NA 22.255347        NA        NA        NA
#> Araliaceae         3.902458e-04 18.093569        NA 20.500000 0.4540667
#> Araucariaceae      1.878580e-03 12.622035        NA 13.000000        NA
#> Arecaceae          4.500000e-04 18.349681        NA 13.216667 0.4556000
#> Aristolochiaceae   3.589556e-03 31.726525        NA        NA        NA
#> Asparagaceae                 NA 22.875870        NA        NA 0.4631000
#> Asphodelaceae                NA 12.352222        NA        NA        NA
#> Aspleniaceae       8.200000e-03 28.260000        NA        NA        NA
#> Asteraceae         1.052882e-02 22.036543        NA  9.346060 0.4375956
#> Asteropeiaceae               NA        NA        NA        NA        NA
#> Atherospermataceae           NA 17.869722        NA 23.200000        NA
#> Athyriaceae                  NA 26.901037        NA        NA        NA
#> Aulacomniaceae               NA  8.000000        NA        NA        NA
#> Austrobaileyaceae            NA        NA        NA        NA        NA
#> Balanopaceae                 NA        NA        NA        NA        NA
#> Balsaminaceae      1.265730e-02 36.003950        NA        NA        NA
#> Begoniaceae                  NA 34.200000        NA        NA        NA
#> Berberidaceae      1.650000e-03 17.997372        NA 21.087500        NA
#> Betulaceae         2.973057e-03 24.194029 14.505150 13.430951 0.4721962
#> Bignoniaceae                 NA 23.830902  6.575973 19.915047 0.4618077
#> Bixaceae                     NA 25.531043        NA        NA        NA
#> Blechnaceae                  NA 11.749051        NA        NA        NA
#> Bonnetiaceae                 NA  9.800000        NA        NA        NA
#> Boraginaceae       5.740212e-03 23.032555        NA        NA 0.4153500
#> Brassicaceae       1.066000e-02 34.316821        NA 18.501306 0.4274775
#> Bromeliaceae                 NA  9.559591        NA        NA        NA
#> Brunelliaceae                NA 21.006960        NA        NA        NA
#> Bruniaceae                   NA  7.781667        NA        NA        NA
#> Burseraceae                  NA 19.028712        NA 10.964427 0.4548544
#> Butomaceae                   NA 42.600000        NA        NA        NA
#> Buxaceae                     NA 22.691830        NA        NA        NA
#> Cabombaceae                  NA 19.500000        NA        NA        NA
#> Cactaceae          2.433862e-05 16.988641        NA        NA        NA
#> Calophyllaceae               NA 12.424790        NA        NA 0.4663000
#> Calycanthaceae               NA 17.400000        NA        NA        NA
#> Calyceraceae                 NA 43.000000        NA        NA        NA
#> Campanulaceae                NA 27.682767        NA  5.721013 0.4498475
#> Canellaceae                  NA        NA        NA        NA        NA
#> Cannabaceae                  NA 28.941055        NA        NA 0.4683583
#>                    Kmax_stemxylem       WUE        P50       Al2As
#> Acanthaceae                    NA        NA         NA  3070.43202
#> Achariaceae                    NA  3.662630         NA          NA
#> Acoraceae                      NA        NA         NA          NA
#> Actinidiaceae                  NA        NA         NA          NA
#> Adoxaceae              4.05355750  3.117284 -3.0384833  4889.57295
#> Aextoxicaceae                  NA        NA         NA          NA
#> Aizoaceae                      NA        NA         NA          NA
#> Akaniaceae                     NA        NA         NA          NA
#> Alismataceae                   NA        NA         NA          NA
#> Alstroemeriaceae               NA        NA         NA          NA
#> Altingiaceae           0.50500000        NA -2.0370147  5129.90090
#> Amaranthaceae          0.08190000 13.478385 -2.4252402   973.13500
#> Amaryllidaceae                 NA 15.949005         NA          NA
#> Amborellaceae          0.54000000        NA -3.0000000  4255.31915
#> Amphorogynaceae                NA        NA         NA          NA
#> Anacardiaceae          4.07720162  2.339255 -2.6535235 19581.70525
#> Anisophylleaceae               NA        NA         NA          NA
#> Annonaceae             5.27066667        NA -2.5276068 10266.65654
#> Aphloiaceae                    NA        NA         NA          NA
#> Apiaceae               0.51500000  3.099757 -5.7000000    81.15013
#> Apocynaceae            2.56516667  9.082218 -2.4864334 19766.80134
#> Aptandraceae                   NA        NA         NA          NA
#> Aquifoliaceae          0.22557572        NA -3.6437782  4886.36183
#> Araceae                        NA        NA         NA          NA
#> Araliaceae             1.68090114        NA -1.6530859  3928.39877
#> Araucariaceae          0.73225000        NA -2.6183226  3846.15385
#> Arecaceae                      NA  3.962648 -1.8100000  5492.53731
#> Aristolochiaceae               NA        NA         NA          NA
#> Asparagaceae                   NA 13.248101 -1.6960000          NA
#> Asphodelaceae                  NA        NA         NA          NA
#> Aspleniaceae                   NA        NA         NA          NA
#> Asteraceae             0.49865714  3.724090 -3.2565860  2421.76507
#> Asteropeiaceae                 NA        NA         NA          NA
#> Atherospermataceae             NA        NA -3.0063333  2435.95630
#> Athyriaceae                    NA        NA         NA          NA
#> Aulacomniaceae                 NA        NA         NA          NA
#> Austrobaileyaceae      2.30000000        NA -0.4990000 15384.61538
#> Balanopaceae                   NA        NA         NA          NA
#> Balsaminaceae                  NA  2.842391         NA          NA
#> Begoniaceae                    NA        NA         NA          NA
#> Berberidaceae          0.08733333  4.039735 -4.5000000   570.06271
#> Betulaceae             2.87333740        NA -2.1017591  6158.73417
#> Bignoniaceae           2.10149000        NA -0.8616667 12439.80827
#> Bixaceae                       NA        NA         NA 12274.01424
#> Blechnaceae                    NA        NA         NA          NA
#> Bonnetiaceae                   NA        NA         NA          NA
#> Boraginaceae                   NA  2.819216 -3.5677066          NA
#> Brassicaceae                   NA  3.632075         NA          NA
#> Bromeliaceae                   NA        NA         NA          NA
#> Brunelliaceae                  NA        NA         NA          NA
#> Bruniaceae             0.25150000        NA -3.3883558          NA
#> Burseraceae            3.40500000  4.390857 -1.3054970 12218.54705
#> Butomaceae                     NA        NA         NA          NA
#> Buxaceae                       NA        NA -8.0000000          NA
#> Cabombaceae                    NA        NA         NA          NA
#> Cactaceae              1.86880952        NA -1.2875000  2554.61304
#> Calophyllaceae         0.89828571        NA -1.5400000  2662.52127
#> Calycanthaceae                 NA        NA -1.2808475          NA
#> Calyceraceae                   NA        NA         NA          NA
#> Campanulaceae                  NA  4.790244         NA          NA
#> Canellaceae                    NA        NA -0.2320000          NA
#> Cannabaceae            4.33961052        NA -1.5325304 29406.04715
#>                    conduit2sapwood
#> Acanthaceae              0.6300000
#> Achariaceae                     NA
#> Acoraceae                       NA
#> Actinidiaceae                   NA
#> Adoxaceae                       NA
#> Aextoxicaceae                   NA
#> Aizoaceae                       NA
#> Akaniaceae                      NA
#> Alismataceae                    NA
#> Alstroemeriaceae                NA
#> Altingiaceae             0.8226667
#> Amaranthaceae                   NA
#> Amaryllidaceae                  NA
#> Amborellaceae                   NA
#> Amphorogynaceae                 NA
#> Anacardiaceae            0.7155889
#> Anisophylleaceae                NA
#> Annonaceae               0.5685000
#> Aphloiaceae                     NA
#> Apiaceae                        NA
#> Apocynaceae              0.7112500
#> Aptandraceae                    NA
#> Aquifoliaceae            0.6528500
#> Araceae                         NA
#> Araliaceae               0.7785000
#> Araucariaceae            0.9375000
#> Arecaceae                       NA
#> Aristolochiaceae                NA
#> Asparagaceae                    NA
#> Asphodelaceae                   NA
#> Aspleniaceae                    NA
#> Asteraceae               0.7219423
#> Asteropeiaceae                  NA
#> Atherospermataceae       0.7560000
#> Athyriaceae                     NA
#> Aulacomniaceae                  NA
#> Austrobaileyaceae               NA
#> Balanopaceae                    NA
#> Balsaminaceae                   NA
#> Begoniaceae                     NA
#> Berberidaceae                   NA
#> Betulaceae               0.8444000
#> Bignoniaceae             0.6360476
#> Bixaceae                        NA
#> Blechnaceae                     NA
#> Bonnetiaceae                    NA
#> Boraginaceae                    NA
#> Brassicaceae                    NA
#> Bromeliaceae                    NA
#> Brunelliaceae                   NA
#> Bruniaceae                      NA
#> Burseraceae              0.8204286
#> Butomaceae                      NA
#> Buxaceae                 0.8330000
#> Cabombaceae                     NA
#> Cactaceae                0.3636905
#> Calophyllaceae           0.7335000
#> Calycanthaceae           0.6535500
#> Calyceraceae                    NA
#> Campanulaceae                   NA
#> Canellaceae                     NA
#> Cannabaceae              0.7598125
#>  [ reached 'max' / getOption("max.print") -- omitted 317 rows ]
```
