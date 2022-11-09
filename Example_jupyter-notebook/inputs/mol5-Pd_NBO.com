%chk=mol5-Pd_NBO_NBO.chk
%mem=24GB
%nprocshared=12
# wb97xd/genecp pop=(nbo6read,savenbos) scrf=(solvent=toluene,smd)

mol5-Pd_NBO_NBO

0 1
P  -1.627933  -0.714127  -0.2491
Pd  0.089108  -2.235772  -0.013457
P  1.676885  -0.585214  0.249237
Cl  -1.523562  -3.976556  0.313677
Cl  1.847927  -3.821958  -0.366352
C  -1.077088  0.740414  -1.246103
C  -2.991991  -1.432062  -1.229029
C  -2.340323  -0.081258  1.294799
C  -0.24784  1.713177  -0.717272
C  -1.382071  0.761831  -2.631588
C  -4.330257  -1.181813  -0.938626
C  -2.666027  -2.255544  -2.310401
C  -3.200463  1.018349  1.300325
C  -2.014915  -0.712994  2.493843
C  0.113379  1.723145  0.732141
C  0.287305  2.736686  -1.56064
C  -0.850555  1.709854  -3.455608
C  -5.330464  -1.733003  -1.732758
C  -3.664977  -2.795542  -3.104274
C  -3.730303  1.479022  2.496083
C  -2.541986  -0.24323  3.689845
C  1.013817  0.81346  1.25688
C  -0.497679  2.699031  1.580235
C  1.127073  3.768114  -1.061526
C  0.001751  2.716708  -2.949275
C  -5.000907  -2.532845  -2.817603
C  -3.397224  0.85149  3.691557
C  1.314805  0.849886  2.642838
C  -1.41058  3.668526  1.085376
C  -0.213456  2.692377  2.969204
C  1.64995  4.715307  -1.898034
C  0.570421  3.702304  -3.794704
C  0.711536  1.749733  3.471395
C  -2.003137  4.569687  1.926324
C  -0.854638  3.628144  3.819335
C  1.375471  4.681269  -3.282537
C  -1.728798  4.548119  3.31116
C  3.085977  -1.215035  1.226087
C  2.343347  0.118257  -1.284538
C  2.813909  -2.054464  2.310239
C  4.406563  -0.896423  0.922846
C  2.049642  -0.507308  -2.494992
C  3.132053  1.270241  -1.271585
C  3.846564  -2.541455  3.094809
C  5.441434  -1.394588  1.707755
C  2.536895  0.020243  -3.683691
C  3.623448  1.788177  -2.460245
C  5.16436  -2.20971  2.795737
C  3.321383  1.166604  -3.666967
H  -2.044674  0.017004  -3.050992
H  -4.608495  -0.580322  -0.081986
H  -1.628695  -2.495004  -2.520401
H  -3.436453  1.531535  0.373701
H  -1.352264  -1.57192  2.489153
H  -1.085812  1.699151  -4.515232
H  -6.370533  -1.541309  -1.491505
H  -3.399745  -3.440345  -3.934827
H  -4.387704  2.341223  2.496527
H  -2.281462  -0.735497  4.620349
H  1.362432  3.801458  -0.004951
H  -5.783434  -2.9666  -3.431269
H  -3.803982  1.220604  4.627095
H  2.03177  0.155406  3.058996
H  -1.646235  3.690494  0.028551
H  2.288095  5.4953  -1.496279
H  0.350751  3.667279  -4.857342
H  0.945293  1.750464  4.531399
H  -2.696406  5.302938  1.527828
H  -0.633863  3.603247  4.882004
H  1.803869  5.434038  -3.935804
H  -2.211998  5.263948  3.967678
H  1.792947  -2.350262  2.528036
H  4.645275  -0.283107  0.06276
H  1.441114  -1.405267  -2.505007
H  3.342519  1.777796  -0.335711
H  3.622836  -3.199626  3.927098
H  6.467818  -1.150026  1.456058
H  2.300356  -0.467272  -4.623069
H  4.225488  2.68977  -2.446227
H  5.974122  -2.602348  3.401661
H  3.69719  1.580335  -4.596752

C Cl H P 0
6-311G(d)
****
Pd     0
S   3   1.00
      2.7870000             -1.6102393
      1.9650000              1.8489842
      0.6243000              0.6037492
S   4   1.00
      2.7870000              1.3540775
      1.9650000             -1.6780848
      0.6243000             -0.8559381
      0.1496000              1.0200299
S   1   1.00
      0.0436000              1.0000000
P   3   1.00
      5.9990000             -0.1034910
      1.4430000              0.7456952
      0.5264000              0.3656494
P   2   1.00
      0.7368000              0.0763285
      0.0899000              0.9740065
P   1   1.00
      0.0262000              1.0000000
D   3   1.00
      6.0910000              0.0376146
      1.7190000              0.5200479
      0.6056000              0.5706071
D   1   1.00
      0.1883000              1.0000000
F   1   1.00
      1.4720000              1.0000000
****

Pd 0
LANL2DZ

$nbo bndidx $end


 