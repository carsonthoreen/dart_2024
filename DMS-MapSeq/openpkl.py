#!/usr/bin/env python3

import pickle
import numpy

with open("denatured1__barcode_Bar2_1_mut_cov_PE.pkl","rb") as f:
    mut, cov = pickle.load(f)
print(f"mut len: {len(mut)}")
print(mut[0][20:290])
print(f"cov len: {len(cov)}")
print(cov[0][20:290])
print(len(cov[0]))


