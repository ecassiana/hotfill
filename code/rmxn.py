#! /usr/bin/python -t
# _*_ coding: iso-8859-1 _*_

MODULE_NAME = "rmxn"
MODULE_DESC = "Linear algebra operations on rectangular numeric matrices"
MODULE_VERS = "1.0"

MODULE_COPYRIGHT = "Copyright © 2009 State University of Campinas"

MODULE_INFO = \
  "A library module to perform linear algebra operations on rectangular numeric matrices.\n" \
  "\n" \
  "  Bla bla.\n"

import sys
import rn

def zero_matrix(m,n) :
  # Zero matrix with {m} rows and {n} cols.
  R = [None]*m;
  for i in range(m) :
    R[i] = [0.0]*n;
  return R;
  # ----------------------------------------------------------------------
  
def ident_matrix(m,n) :
  # Identity matrix with {m} rows and {n} cols.
  R = [None]*m;
  for i in range(m) :
    R[i] = [0.0]*n;
    if (i < n) : R[i][i] = 1.0;
  return R;
  # ----------------------------------------------------------------------

def diag_matrix(x) :
  # Square diagonal matrix with {x} as the diagonal.
  n = len(x);
  R = zero_matrix(n,n);
  for i in range(n) :
    R[i][i] = x[i];
  return R;
  # ----------------------------------------------------------------------
  
def add(M,N):
  # Elementwise sum {M+N}.
  m =len(M); assert len(N) == m
  n = len(M[0]); assert len(N[0]) == n
  R = [None]*m
  for i in range(m):
    R[i] = rn.add(M[i],N[i])
  return R
  # ----------------------------------------------------------------------

def sub(M,N):
  # Elementwise difference {M-N}.
  m =len(M); assert len(N) == m
  n = len(M[0]); assert len(N[0]) == n
  R = [None]*m
  for i in range(m):
    R[i] = rn.sub(M[i],N[i])
  return R
  # ----------------------------------------------------------------------

def mul(M,N) :
  # Multiplies the matrices {M} and {N}.
  m = len(M);
  n = len(N[0]);
  R = [None]*m;
  for i in range(m) :
    R[i]= map_row(M[i],N)
  return R;
  # ----------------------------------------------------------------------

def map_row(x,M) :
  # Multiplies the vector {x} by the matrix {M}.
  m = len(x);
  assert len(M) == m, "incompatible {x} lenght and {M} rows";
  n = len(M[0]);
  r = [None] * n;
  for j in range(n) :
    s = 0;
    for i in range(m) :
      s = s + x[i] * M[i][j];
    r[j] = s
  return r;
  # ----------------------------------------------------------------------

def map_col(M,x) :
  # Multiplies the matrix {M} by the vector {x}
  n = len(x);
  m = len(M);
  r = [None] * m;
  for i in range(m) :
    r[i] = rn.dot(M[i], x);
  return r;
  # ----------------------------------------------------------------------

def norm_sqr(M):
  # The sum of the squares of the elements.
  m = len(M)
  sum = 0
  for i in range(m):
    sum += rn.norm_sqr(M[i])
  return sum
  # ----------------------------------------------------------------------

def inv(M):
  n = len(M)
  assert len(M[0]) == n, "non-square matrix"
  if n == 1:
    N = [ [ 1/M[0][0] ] ]
  elif n == 2:
    N = inv_2(M)
  elif n == 3:
    N = inv_3(M)
  else:
    N = inv_gen(M)
  # check_inv(M, N, n*2.0e-4)
  return N
  # ----------------------------------------------------------------------
  
def check_inv(M,N,tol):
  # Checks whether {N} is the inverse of {M}, within a tolerance {tol}
  sys.stderr.write("--- checking result of {inv}---\n")
  n = len(M)
  d = norm_sqr(sub(ident_matrix(n,n), mul(M, N)))
  if d > tol*tol:
    sys.stderr.write("d = %14.9f\n" % d)
    assert False
  return
  # ----------------------------------------------------------------------
  
def inv_2(M):
  assert len(M) == 2
  assert len(M[0]) == 2
  # Invert {M} by hand:
  D = M[0][0]*M[1][1] - M[0][1]*M[1][0]
  N = \
    [ 
      [ +M[1][1]/D, -M[0][1]/D, ],
      [ -M[1][0]/D, +M[0][0]/D, ]
    ]
  return N
  # ----------------------------------------------------------------------
  
def inv_3(M):
  assert len(M) == 3
  assert len(M[0]) == 3
  # Invert {M} by hand:
  N = \
    [ 
      [ + (M[1][1]*M[2][2] - M[1][2]*M[2][1]), 
        - (M[0][1]*M[2][2] - M[0][2]*M[2][1]),
        + (M[0][1]*M[1][2] - M[0][2]*M[1][1]),
      ],
      [ - (M[1][0]*M[2][2] - M[1][2]*M[2][0]), 
        + (M[0][0]*M[2][2] - M[0][2]*M[2][0]),
        - (M[0][0]*M[1][2] - M[0][2]*M[1][0]),
      ],
      [ + (M[1][0]*M[2][1] - M[1][1]*M[2][0]), 
        - (M[0][0]*M[2][1] - M[0][1]*M[2][0]),
        + (M[0][0]*M[1][1] - M[0][1]*M[1][0]),
      ],
    ]
  J = mul(N,M)
  # sys.stderr.write("J = %s\n" % str(J))
  D = J[0][0]
  for i in range(3):
   for j in range(3):
     N[i][j] = N[i][j]/D
  return N
  # ----------------------------------------------------------------------

def inv_gen(M):
  n = len(M)
  assert len(M[0]) == n
  
  # Computes a Householder reflection {R}
  # that takes the last row to a vector {[0,0,..,0,A]}:
  R = householder(M[n-1])
  
  # Applies {R} to {M}:
  Q = mul(M,R)
  
  # Split {Q} into {N} = first {n-1} rows and cols, {x} = last col, {A} = last elem.
  N = [ Q[i][:n-1] for i in range(n-1) ]
  x = [ Q[i][n-1] for i in range(n-1) ]
  A = Q[n-1][n-1]
  assert rn.norm_sqr(Q[n-1][:n-1]) < 1.0e-14
  
  # Reacurse:
  Ninv = inv(N)
  y = rn.scale(-1/A,map_col(Ninv,x))
  
  # Assemble the Q inverse: 
  Qinv = [ Ninv[i] + [y[i]] for i in range(n-1) ]
  Qinv.append([0]*(n-1) + [1/A])
  
  # Apply the householder map to rows:
  Minv = mul(R,Qinv)
  
  return Minv
  # ----------------------------------------------------------------------
  
def householder(u):
  # Given an {n}-vector {u}, returns an {n} by {n} 
  # matrix that is a reflection that takes the row vector
  # {u} to the vector {[0,0,.. 0,A]} where {A} is the squared norm of {u}.
  
  n = len(u)
  A = rn.norm(u)
  v, d = rn.dir(rn.sub(u, [0]*(n-1)+[A]))
  if d < 1.0e-8:
    R = ident_matrix(n,n)
    R[0][0] = -1
  else:
    VV = [ [ 2*v[i]*v[j] for j in range(n) ] for i in range(n) ]
    R = sub(ident_matrix(n,n), VV)
  # check_householder(u, R, 1.0e-7)
  return R
  # ----------------------------------------------------------------------
  
def check_householder(u,R, tol):
  # Checks whether {R} is the householder matrix for {u}, with tolerance {tol}, 
  sys.stderr.write("--- checking result of {householder}---\n")
  n = len(u)
  
  u1 = map_row(u,R)
  err0_sqr = rn.norm_sqr(u1[:n-1])
  errn = u1[n-1] - rn.norm(u)
  if err0_sqr > tol*tol:
    sys.stderr.write("u  = %s\n" % (str(u)))
    sys.stderr.write("u1 = %s\n" % (str(u1)))
    sys.stderr.write("err0_sqr = %.16f\n" % (err0_sqr,errn))
    assert err0_sqr <= tol*tol
  assert abs(errn) < tol
  
  T = mul(R,R)
  assert norm_sqr(sub(ident_matrix(n,n), T)) < tol*tol
  return R
  # ----------------------------------------------------------------------
  check_householder(u, R, 1.0e-7)

def print(wr, pref, M, fmt, sep, suff):
  # Prints lines of {M} with prefix {pref}, elem separator {sep}, suffix {suff}
  m = len(M)
  n = len(M[0])
  for i in range(m):
    wr.write(pref)
    for j in range(n):
      if j > 0: wr.write(sep)
      wr.write(fmt % M[i][j])
    wr.write(suff)
  

def test():
  sys.stderr.write("--- testng {rmxn} module---\n")

  n = 5
  sys.stderr.write("n = %d\n" % n)
  
  # Test {householder}:
  u = [ 1/(i+1) for i in range(n) ]
  R = householder(u)
  check_householder(u, R, 1.0e-7)
  
  # Test inv:
  M = [ [ 1/(i+j+1) for j in range(n) ] for i in range(n) ]
  N = inv(M)
  check_inv(M, N, n*2.0e-4)

  return
  # ----------------------------------------------------------------------

# BOXES

def box_affine(B, A,b):
  # Returns the smallest box that encloses the box {B} mapped by the
  # affine map {A,b}, that is {p --> p*A + b}.
  # If {B} is {None}, returns {None}.
  #
  if B == None:
    Bmap = B
  else:
    Bmap = None
    for xc in B[0][0], B[1][0]:
      for yc in B[0][1], B[1][1]:
        c = (xc, yc)
        cmap = rn.add(map_row(c, A), b)
        Bmap = rn.box_include_point(Bmap, cmap)
    return Bmap
  # ----------------------------------------------------------------------

# ----------------------------------------------------------------------
