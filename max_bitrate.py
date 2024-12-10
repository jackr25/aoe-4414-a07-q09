# max_bitrate.py
#
# Usage: python max_bitrate.py tx_w tx_gain_db freq_hz dist_km rx_gain_db n0_j bw_hz
#  Text explaining script usage
# Parameters:
# tx_w:
# tx_gain_db 
# freq_hz 
# dist_km 
# rx_gain_db 
# n0_j 
# bw_hz
# Output:
#  A description of the script output
#
# Written by Jack Rathert
# Other contributors: None
#


## imports---------------------------------------------------------------------------
import sys
import math as m

## constants variable input----------------------------------------------------------
Ll=0.7943
La = 1
c = 2.99792458e8
## Classes---------------------------------------------------------------------------
class numpy_lite:
    def __init__(self):
        from math import sqrt
        self.sqrt = sqrt 
    def matrix_mult(self, matrix1, matrix2):
        # Wrap flat lists into 1xN matrices
        if isinstance(matrix1[0], (int, float)):
            matrix1 = [[elem] for elem in matrix1]  # Convert to Nx1 matrix
        if isinstance(matrix2[0], (int, float)):
            matrix2 = [[elem] for elem in matrix2]  # Convert to Nx1 matrix
        # Get matrix dimensions
        rowsA, colsA = len(matrix1), len(matrix1[0])
        rowsB, colsB = len(matrix2), len(matrix2[0])
        # Nx1 * 1xN -> NxN (outer product)
        if colsA == 1 and rowsB == 1:
            return [[matrix1[i][0] * matrix2[0][j] for j in range(colsB)] for i in range(rowsA)]
        # 1xN * Nx1 -> scalar (dot product)
        if rowsA == 1 and colsB == 1 and colsA == rowsB:
            return sum(matrix1[0][i] * matrix2[i][0] for i in range(colsA))
        # General matrix multiplication
        if colsA != rowsB:
            raise ValueError(f"Incompatible dimensions: {colsA} != {rowsB}")
        # Perform matrix multiplication
        return [[sum(matrix1[i][k] * matrix2[k][j] for k in range(colsA)) for j in range(colsB)] for i in range(rowsA)]
    def matrix_add(self,list1,list2):
      return [x+y for x,y in zip(list1, list2)]
    def matrix_sub(self,list1,list2):
      return [x-y for x,y in zip(list1,list2)]
    def mag(self,list1):
      return sqrt(sum(x**2 for x in list1)) # type: ignore
    def smul(self,scalar, vector):
      return [scalar * x for x in vector]
    def vecadd(self,vector1,vector2):
      if len(vector1) != len(vector2):
         raise ValueError("Vectors must be of the same length.")
      return [x + y for x, y in zip(vector1, vector2)]
    def vecsub(self,vector1,vector2):
      if len(vector1) != len(vector2):
         raise ValueError("Vectors must be of the same length.")
      return [x - y for x, y in zip(vector1, vector2)]  
    def dotprod(self,vector1,vector2):
      if len(vector1) != len(vector2):
         raise ValueError("Vectors must be of the same length.")
      return sum(x * y for x, y in zip(vector1, vector2))
nplite = numpy_lite()


## Functions--------------------------------------------------------------------------
def calc_C(P,Gt,lam,S,Gr):
  C = P*Ll*Gt*(lam/(4*m.pi*S*1000))**2 *La*Gr
  return C

# Arguments---------------------------------------------------------------------------

# parsing
if len(sys.argv)==8:
  tx_w = float(sys.argv[1])
  tx_gain_db = float(sys.argv[2])
  freq_hz = float(sys.argv[3])
  dist_km = float(sys.argv[4])
  rx_gain_db = float(sys.argv[5])
  n0_j = float(sys.argv[6])
  bw_hz = float(sys.argv[7])
  ...
else:
  print(\
   'Usage: '\
   ' python3 max_bitrate.py tx_w tx_gain_db freq_hz dist_km rx_gain_db n0_j bw_hz'\
  )
  exit()

## Function Calls ----------------------------------------------------------------------

N = n0_j*bw_hz
lam = c/freq_hz
Gt = 10**(tx_gain_db/10)
Gr = 10**(rx_gain_db/10)
C = calc_C(tx_w,Gt,lam, dist_km, Gr)

r_max = bw_hz*m.log2(1+C/N)
print(m.floor(r_max))