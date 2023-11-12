import numpy as np

rs = [0, 0, 2, 3, 4, 3, 9, 8]
cs = [0, 4, 1, 3, 5, 1, 7, 8]
vs = [2.0, -1.0, 3.0, 5.0, -1.0, 4.0, 2.0, -8.0]

class COOMatrix:
    def __init__(self, m, n, rs, cs, vs):
      self._m = m
      self._n = n
      self._rs = np.array(rs)
      self._cs = np.array(cs)
      self._vs = np.array(vs)

   # --- m -----------------------------------------------
    @property
    def m(self):
        """Get the 'm' property."""
        return self._m

    @m.setter
    def m(self, value):
        self._m = value

    @m.deleter
    def m(self):
        del self._m
   # --- n -----------------------------------------------
    @property
    def n(self):
        """Get the 'n' property."""
        return self._n

    @n.setter
    def n(self, value):
        self._n = value

    @n.deleter
    def n(self):
        del self._n
   # --- rs -----------------------------------------------
    @property
    def rs(self):
        """Get the 'rs' property."""
        return self._rs

    @rs.setter
    def rs(self, value):
        self._rs = np.array(value)

    @rs.deleter
    def rs(self):
        del self._rs
# --- cs -----------------------------------------------
    @property
    def cs(self):
        """Get the 'cs' property."""
        return self._cs

    @cs.setter
    def cs(self, value):
        self._cs = np.array(value)

    @cs.deleter
    def cs(self):
        del self._cs
# --- vs -----------------------------------------------
    @property
    def vs(self):
        """Get the 'vs' property."""
        return self._vs

    @vs.setter
    def vs(self, value):
        self._vs = np.array(value)

    @vs.deleter
    def vs(self):
        del self._vs

    def apply(self, vec):
        my_ve
        vec_copy = np.array(vec)
        res = np.zeros(np.shape(vec))
        print(f"res: {res}")
        for ii, r in enumerate(rs):
            res[r] += vs[r] * vec_copy[r]
            print(f"{r}: {res[r]}")
            # print(f"ii:{ii} r:{r} res:{res[ii]} += vs: {vs[ii]} * {vec_copy[ii]}")
            # print(f"{ii}: {r} {vs[ii]} x {vec_copy[ii]}")
        # print(f"res[{ii}]={res[ii]}")
        return res
m = 10
n = 10
coo = COOMatrix(m, n, rs, cs, vs)
vec = [1,2,3,4,5,6,7,8,9,10]
ans = coo.apply(vec)