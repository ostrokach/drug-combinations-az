



http://folk.uio.no/sigven/az_dream/



Round 1
-------

1A
^^^^
5129062
  - 0.07	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
  
5157022
  - ML with basic parameters (no gridsearch)
  - 0.09	0.21	0.71	0.37	0.65	0.26	0.65	0.5	0.59
  - 71 / 149 global correlation
  - 40 / 149 mean correlation
  - 11 / 149 mean correlation (top 30%)
  - 36 / 149 mean correlation (top 20%)
  - 10 / 149 mean correlation (top 10%)


5177062
  - ML with gridsearch 
  - mean synergy score for pair +/- SE, mean synergy score for cell line +/- SE, chemo
  - 0.08	0.19	0.7	0.23	0.79	0.15	0.84	0.02	1.0
  - 54 / 149 mean correlation
  - 44 / 149 mean correlation (top 30%)
  - 75 / 149 mean correlation (top 20%)
  - 
  
  

1B
^^^^
5129282
  - 0.07	0.0	0.0	0.0	0.0	0.0	0.0	0.0	0.0
  
5157653
  - ML with basic parameters (no gridsearch)
  - 0.08	0.16	0.72	0.19	0.65	0.01	0.62	0.09	0.64
  
5177064
  - ML with gridsearch
  - mean synergy score for pair +/- SE, mean synergy score for cell line +/- SE
  - 0.05	0.17	0.67	0.18	0.74	0.23	0.79	0.22	0.92


2
^^^^



To Do
-----

Evaluation script
  - Figure out the metric that they use for scoring participants.
  - `Concern about synergy scores and evaluation strategy <http://support.sagebase.org/sagebase/topics/concern-about-synergy-scores-and-evaluation-strategy>`_.
  - 



Challange 1
-----------

Leave one out XV, but calculating the means **without** the pair that was left out. 

Take the weighted mean of predictions for all drug pairs.

Add many cell line similarity features.



Challange 2
------------



Info
=====

BRAF_mut -> BRAF
BRAF_V600E -> BRAF
CD19 antibody -> CD19
VEGFR2 -> KDR
NIAP -> NAIP
TNFA -> TNF
NAE2 -> UBA3
TIE2 -> TEK
cMET -> MET
Gamma secretase -> APH*
Proteasome -> PSM*



Methylation -> DNMT*,MGMT


DNA ->  


# use stitch
microtubule -> TUBA*,TUBB*,[MAP10,MAP1A,MAP1B,MAP1LC3A,MAP1LC3B,MAP1LC3B2,MAP1LC3C,MAP1S,MAP2,MAP4,MAP6,MAP7,MAP9,MAPRE1,MAPRE2,MAPRE3,MAPT]


Topotecan -> TOP1,BCG2,HIF1A,EGFR,TP53,ABCB1,CASP3,CYP3A4,ATM,TOP1MT
Gemcitabine -> DNA polymerase
Oxaliplatin -> DNA repair
Cisplatin -> DNA repair

FOLFIRI -> TOP1*, DNA polymerase

CCNH,CDK7,CETN2,DDB1,DDB2,ERCC1,ERCC2,ERCC3,ERCC4,ERCC5,ERCC6,ERCC8,LIG1,MNAT1,MMS19,RAD23A,RAD23B,RPA1,RPA2,TFIIH,XAB2,XPA,XPC




-------

Try to get a list of mutations that each drug targets...




Manual Intervention
====================

BRAF_M	BRAF_mut
BRAF_M2	BRAF_V600E




sorafenib (Braf/VEGFR2 kinase inhibitor)






Teams
-----

`DrugAddicts <https://www.synapse.org/#!Team:3331391>`_
  - `Anna Goldenberg <http://www.cs.toronto.edu/~goldenberg/Anna_Goldenberg/Home.html>`_
  - `Benjamin Haibe-Kains <https://www.pmgenomics.ca/bhklab/people>`_ (and 2 of his post-docs)


`OurTeam <https://www.synapse.org/#!Team:3332165>`_
  - Headed by `Lodewyk Wessels <http://bioinformatics.tudelft.nl/users/lodewyk-wessels>`_.
  - `Magali Michaut <http://people.rez-gif.supelec.fr/mmichaut/welcome.php>`_

Mikhail Zaslavskiy

North Atlantic Dream

GMM

PrecisionHunter

Mushthofa's http://support.sagebase.org/people/mushthofa_7539069



