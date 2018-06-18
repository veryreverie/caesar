# ======================================================================
# Shared functions between python plotting scripts.
# ======================================================================
def dblecomplex(displacement):
  if displacement[0]=='-':
    re = float(displacement[:24])
    im = float(displacement[24:-1])
  else:
    re = float(displacement[:23])
    im = float(displacement[23:-1])
  return complex(re,im)

