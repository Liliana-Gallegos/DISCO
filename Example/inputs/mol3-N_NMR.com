%chk=mol3-N_NBO_NMR.chk
%mem=24GB
%nprocshared=12
# nmr=GIAO mpw1pw91/6-311+G(2d,p)

mol3-N_NBO_NMR

0 1
N  -1.757791  -0.609892  -0.065956
C  -0.437869  -0.282213  -0.026696
C  0.013089  1.047766  -0.02429
C  1.376923  1.292706  -0.003811
C  1.883764  -0.91042  0.019635
C  0.554611  -1.279934  -0.00205
H  -0.679644  1.880496  -0.038742
H  1.730151  2.321623  -0.003034
H  2.652468  -1.67957  0.039437
H  0.280221  -2.330544  -0.002465
H  -1.992532  -1.57657  0.078529
C  -2.822014  0.357764  0.045528
H  -3.776003  -0.168213  0.015967
H  -2.772669  0.927953  0.981688
H  -2.805298  1.064141  -0.790402
N  2.322402  0.353702  0.018689

