# phamclust

PhamClust is a tool for performing gene phamily based clustering of bacteriophage genomes. It makes use of a novel genome 
similarity index, the proteomic equivalence quotient (PEQ) to cluster genomes according to their global similarity. It was 
[published in mSphere](https://doi.org/10.1128/msystems.00443-23) in 2023.

# Installation

The most straightforward way to install PhamClust is from the Python Package Index (PyPI): `pip install phamclust`

# Usage

Once installed, invoking `phamclust` at the commandline with no arguments will display the help menu:

```
(phamclust) chg60 % phamclust
usage: phamclust [-h] [-g] [-k] [-s] [-sl] [-c] [-cl] [-nr] [-nl] [-m] [-d] [-n] [-r] [-t] infile outdir

Cluster phage genomes using gene content similarity-based metrics.

positional arguments:
  infile                path to a TSV file mapping genomes to phams and translations
  outdir                path to which output files should be written

options:
  -h, --help            show this help message and exit
  -g, --genome-dir      interpret `infile` as a directory of genome FASTA files instead of TSV
  -d, --debug           increase verbosity of logging for debug purposes
  -n, --no-sub          do not perform sub-clustering
  -r, --remove-tmp      remove temporary files (not recommended if repeated runs are planned on the same dataset)
  -t , --threads        number of CPU cores to use [default: 16]

clustering arguments:
  -k , --k-min          minimum cluster size to perform subclustering [default: 6]
  -s , --sub-thresh     similarity threshold to use for sub-clustering [default: 0.6]
  -sl , --sub-linkage   linkage type to use for sub-clustering [default: single]
  -c , --clu-thresh     similarity threshold to use for clustering [default: 0.25]
  -cl , --clu-linkage   linkage type to use for clustering [default: average]
  -nr , --nr-thresh     similarity threshold above which to pre-group very similar genomes that must be clustered together [default: 0.75]
  -nl , --nr-linkage    linkage type to use for pre-grouping very similar genomes [default: complete]
  -m , --metric         relatedness index to use for pairwise genome comparisons [default: peq]
  
heatmap arguments:
  -hc , --heatmap-colors
                        comma-separated list of 2 or 3 colors to use in heatmaps [default: red,yellow,green]
  -hm , --heatmap-midpoint
                        midpoint to use for color gradient in heatmaps [default: same as clustering threshold]

Available metrics:

    Acronym     Name                                Reference
(1) gcs         gene content similarity             https://doi.org/10.1038/nmicrobiol.2017.112
(2) jc          jaccard coefficient                 https://doi.org/10.1111/j.1469-8137.1912.tb05611.x
(3) pocp        percentage of conserved proteins    https://doi.org/10.1128/JB.01688-14
(4) af          alignment fraction                  https://doi.org/10.1093/nar/gkv657
(5) aai         average aminoacid identity          https://doi.org/10.1073/pnas.0409727102
(6) peq         proteomic equivalence quotient      https://doi.org/10.1128/msystems.00443-23
```

PhamClust takes an input filepath and an output filepath as its primary arguments. By default, the input path is 
assumed to be a TSV file mapping genome names/identifiers to pham identifiers and translations.

If the `-g/--genome-dir` argument is provided, PhamClust interprets the input filepath as a directory containing 
one FASTA file per genome, with headers structured as in the example below:

```
>name=Bipper|pham=pham_1|n=1
MTAPLLQSVTADDGNMITVPTLQFTRWLDETRDKVIGADGAPDPVRDPMSAYRYLKGRRSVIEGAARQRPMLRLFDKNMDPIAQIAGERLASVEEMMSDSGQANVVLRYDNWLTDFILHQTKIHEDLHLVVDPNPTNRTWRTRWGGKITGINAKRDSSGIHTLELEAISNRQHAKHMLFASNPVFPPEIQLPKMWVLPGNTRTILSISMFVNLARRFFPLLSIPTNIFNPMAWVNGWGAGLDPLMWPLQVAFVNPLLDQSRLSVLGSSWTDWHTAMDSMLKDAGVLFRAYTWLTEDADTPHTELVDMVRGLGPLQDTVDNLTRPHRNCVVFALEDKSGVQGPTGTAADGVINLIGATADDMITETLFNLDRDGDGETDPIFRKLLGVAPEKPKTIWYDGQFSGIIESEIRRHKGPVKKIHTGGRSPSILNQAQTYAIRYALSQLAQVISYGIGAYQQYGTEGLDNLYQGQLDNTLFAWQAFDDPIRALQTGDMAWQEHFERGSGTAYTLSGIVTLRVGHYKTRAWQGFTVKVVNGRPHAVDVDITLGDRAGFEQGGIIFVDQITAIKRSWSRTEPVTVQLSIGDDQDKEDPAARGLRAIQAVWTTLGMLLGEGTIF
>name=Bipper|pham=pham_37|n=1
MTSPSGVAVAALKGHTKPRLYTPPLAVNCNIWIAPELSCPCGCGLHAGTSWGFDCIDFLTNVLKWQLIPYQRWLYIHALEKGPGGEGFRFKTLVILIARQNGKTQWLRGLGLWRLYLDSRGRSSPDCPAAKTVVIAAQGLEYAEGTLGEVVNDVKECRALKREFLRHRQTNGKHAMLLSGRRSWRAVAANRKGGRSMSVDLAELDELREHHDWLAWNAITPTTQARQYSQNVAASNAGDKRSVVLRSLRDGAMAKILARDTEDTKTGLFEYSAPQDANPLERKYWPMANPALGYLPGHDEDALAAKAEAMADNMAGFVTEHLCQWVDTLLPGVMPMEDWNATTDPESRRAEGAPVYAAVDVSHSRSKAYIAVASRRSDGLLHVEVVAAHRGTDWVVPWFKARPGKFVAVAVQARGCPASDLIEPLTEAGVPVMELGGAELVRGAGGVLFDGIRKHAIWHRPSPALDTAAKGTVSRSLGGDTWVLDRKNSPVDAAPLVACAAAAWAEGQGPMVPDKVPEVHEWPDEEEIAEWEKELDELQ
...
```

If a genome encodes more than one copy of a pham (namely, paralogs), each copy after the first should increment 
the value of `n`.

Six metrics are available for calculating intergenomic similarities (processing speeds estimated in parentheses 
as the number of genome pairs calculated per second on an M1 Macbook Pro):

1. Gene Content Similarity (gcs)                    (>100,000 pairs/second)
2. Jaccard Coefficient (jc)                         (>55,000 pairs/second)
3. Percentage of Conserved Proteins (pocp)          (>25,000 pairs/second)
4. Alignment Fraction (AF)                          (>15,000 pairs/second)
5. Average Aminoacid Identity (aai)                 (~500 pairs/second)
6. Proteomic Equivalence Quotient (peq)             (~475 pairs/second)

Because of how cheap the first four metrics are, the maximum parallelization benefit is seen with just 4 CPU 
cores for datasets consisting of fewer than 2,500 genomes. The last two metrics are considerably more expensive 
to calculate, and will see benefit from as many physical cores are available on your machine (i.e., core count, 
not thread count), as long as there is enough system memory. 

# Heatmap settings

Among the outputs from PhamClust are per-cluster matrix heatmaps that nicely illustrate the pairwise similarities 
within clusters/subclusters.

For reasonably small datasets (those with fewer than 1000 genomes), a heatmap will also be drawn for the complete 
dataset matrix.

Two commandline arguments can be used to alter the appearance of these heatmaps. The `-hc/--heatmap-colors` 
argument allows you to specify either two or three colors that define the heatmap colorscheme. By default, a 
3-point gradient is used, where the lowest similarity genome pairs are in red, intermediate similarity genomes 
pairs are in yellow, and the highest-similarity genome pairs are green. This default behavior is the same as 
invoking phamclust with `-hc red,yellow,green`. Another visually appealing option might be a 2-color gradient
from white (0% similar) to green (100% identical), which could be applied with `-hc white,green`.

Any valid CSS named colors can be used:

                aliceblue, antiquewhite, aqua, aquamarine, azure,
                beige, bisque, black, blanchedalmond, blue,
                blueviolet, brown, burlywood, cadetblue,
                chartreuse, chocolate, coral, cornflowerblue,
                cornsilk, crimson, cyan, darkblue, darkcyan,
                darkgoldenrod, darkgray, darkgrey, darkgreen,
                darkkhaki, darkmagenta, darkolivegreen, darkorange,
                darkorchid, darkred, darksalmon, darkseagreen,
                darkslateblue, darkslategray, darkslategrey,
                darkturquoise, darkviolet, deeppink, deepskyblue,
                dimgray, dimgrey, dodgerblue, firebrick,
                floralwhite, forestgreen, fuchsia, gainsboro,
                ghostwhite, gold, goldenrod, gray, grey, green,
                greenyellow, honeydew, hotpink, indianred, indigo,
                ivory, khaki, lavender, lavenderblush, lawngreen,
                lemonchiffon, lightblue, lightcoral, lightcyan,
                lightgoldenrodyellow, lightgray, lightgrey,
                lightgreen, lightpink, lightsalmon, lightseagreen,
                lightskyblue, lightslategray, lightslategrey,
                lightsteelblue, lightyellow, lime, limegreen,
                linen, magenta, maroon, mediumaquamarine,
                mediumblue, mediumorchid, mediumpurple,
                mediumseagreen, mediumslateblue, mediumspringgreen,
                mediumturquoise, mediumvioletred, midnightblue,
                mintcream, mistyrose, moccasin, navajowhite, navy,
                oldlace, olive, olivedrab, orange, orangered,
                orchid, palegoldenrod, palegreen, paleturquoise,
                palevioletred, papayawhip, peachpuff, peru, pink,
                plum, powderblue, purple, red, rosybrown,
                royalblue, saddlebrown, salmon, sandybrown,
                seagreen, seashell, sienna, silver, skyblue,
                slateblue, slategray, slategrey, snow, springgreen,
                steelblue, tan, teal, thistle, tomato, turquoise,
                violet, wheat, white, whitesmoke, yellow,
                yellowgreen

For 3-point color scales, the `-hm/--heatmap-midpoint` argument can be used to adjust the similarity threshold 
where the middle color goes. For example, to highlight diversity within clusters (i.e., dissimilarity between 
subclusters), a midpoint at the subcluster threshold would maximize contrast between intra-subcluster and 
inter-subcluster similarity values.